clc; clear all; close all
%% example code direct search for high-dimensional mod emultiplers
% original paper title: "Unlocking Mode Programming with Multi-Plane Light
% Conversion using Computer-Generated Hologram Optimisation"
% currently under review

% Authors of this repository: Stefan Rothe and Filipe M. Ferreira
% version: 1, 05/22/2024

% to start with our code, load a initial guess produced with wavefront 
% matching algorithm (WMA). We took the code published in Fontaine, N. K. 
% et al., "Laguerre-Gaussian mode sorter", Nat Comm., 10.1 (2019): 1865.

% their code can be found via
% Fontaine, N. et al. Laguerre–Gaussian mode sorter supplemental 
% information: https://doi.org/10.14264/uql.2019.81 (2019).

%% load data and set flags
% here's our WMA example for 5 target modes and 5 passages. All other 
% parameters like input MFD or plane spacing can be found in our manuscript
load WMA_5modes_5passages.mat;

% here's a direct search (DS) example we ran over 100k iterations 
load workspace_DS_5passages_5targetModes_100k_iterations

% if you only want to see our results, switch flag_produce_new_results=0
% and switch flag_plot=1
flag_produce_new_results=0;     % produce new data? Y:1, N:0
flag_plot=1;                    % plot progression and result? Y:1, N:0

if flag_produce_new_results==1
    % to produce new results, reset the phasemasks
    load("WMA_5modes_5passages.mat","MASKS")
    %% Initialisation
    N_fiberModes=size(lp_modes,1);
    N_targetModes=length(target_modes);
    
    N_passages=size(MASKS,1);
    
    dx=size(MASKS,2);
    dy=size(MASKS,3);
    
    N_changes=0;
    
    MER_avg=-inf;           % initialise average MER
    MER_std=inf;            % initialise std MER
    IL_lim=-5;              % initialise IL lower bound
    
    N_it=1000;              % number of iterations. To not overwhelm from 
                            % the beginning, we set it to 1000. To see 
                            % meainingful results, we recommend at least 
                            % 50k itrerations. In our example code, we 
                            % used 100k.
                            
    store_vals=5;           % store progression into arrays every store_vals 
                            % iterations. We recommend a ratio
                            % N_plot=N_it/store_vals=200

    N_plot=N_it/store_vals;     % ratio between total iterations and length 
                                % of array that stores progression

    index_plot=1;

    MER_vals=zeros(1,N_plot);       % store progression
    MER_std_vals=MER_vals;          % store progression
    IL_vals=MER_std_vals;           % store progression
    
    n_save=N_it;           % save every n_save iterations
    if flag_plot==1
        figure;
    end
    %% start algorithm
    for index_direct_search=1:N_it
        % calculate random phase change
        mask_nr=ceil(N_passages*rand);
        pixel_nr_1=ceil(dx*rand);
        pixel_nr_2=ceil(dy*rand);
    
        phase_val_new=rand*2*pi;
        
        % store default phase value
        phase_val_default=angle(MASKS(mask_nr,pixel_nr_1,pixel_nr_2));
        
        if index_direct_search>1
            % apply new phase value after first iteration --> algorithm needs
            % one initial execution
            MASKS(mask_nr,pixel_nr_1,pixel_nr_2)=single(exp(1i*phase_val_new));
        end
    
        % propagate all modes to last plane
        result=single(zeros(N_targetModes,dy,dx));
        for i = 1:N_targetModes
            field_temp=squeeze(input_fields(i,:,:));
            for plane_index=1:4
                SLM_mask_temp=exp(-1i.*angle(squeeze(MASKS(plane_index,:,:))));
                
                field_temp=field_temp.*SLM_mask_temp;
                field_temp = propagate(field_temp,kernel);
            end
            result(i,:,:)=single(field_temp);
        end
    
        % calculate coupling matrix CM
        CM=zeros(N_targetModes,N_fiberModes);
        for modeIdx=1:N_targetModes
                %get the conjugate field in the forward direction at the last plane
                fieldIn = conj(squeeze(result(modeIdx,:,:)));
                %Apply the phase mask to the field
                fieldIn = fieldIn.*squeeze(exp(1i.*angle(MASKS(end,:,:))));
                %For every mode out
                for modeIdy=1:N_fiberModes
                    %get the field in the backward direction at the last plane
                    fieldOut = squeeze(lp_modes(modeIdy,:,:));
                    %Calculate the overlap integral between the two (fieldIn already
                    %conjugated)
                    CM(modeIdx,modeIdy) = abs(sum(sum(fieldIn.*fieldOut)));
                end
        end
    
        CM_log=10*log10(abs(CM).^2);
    
        if index_direct_search==1
            CM_log_default=CM_log;
        end
    
        S_it=svd(CM);
        S2_it = S_it.^2;
    
        %Insertion loss is the mean singular value squared
        IL_it = 10.*log10(mean(S2_it));
    
        % calculate MER
        MER=zeros(1,size(CM,1));
        for i=1:size(CM,1)
            vecTemp=CM(i,:);
            mode_val=vecTemp(target_modes(i));
            vec_clear=vecTemp;
            vec_clear(target_modes(i))=[];
            MER_temp=20*log10(mode_val/sum(vec_clear));
            MER(i)=MER_temp;
        end
    
        MER_avg_it=mean(MER);
        MER_std_it=std(MER);
    
        % save optimization parameters
        if mod(index_direct_search,store_vals)==0
            MER_vals(index_plot)=MER_avg_it;
            MER_std_vals(index_plot)=MER_std_it;
            IL_vals(index_plot)=IL_it;
            index_plot=index_plot+1;
            if flag_plot==1
                subplot(1,3,1);plot(1:index_plot,MER_vals(1:index_plot),'o','LineWidth',2);xlabel(['iterations \times',num2str(store_vals)]); ylabel('MER [dB]'); set(gca, 'FontSize',18); grid on; drawnow
                subplot(1,3,2);plot(1:index_plot,IL_vals(1:index_plot),'o','LineWidth',2);xlabel(['iterations \times',num2str(store_vals)]); ylabel('IL [dB]'); set(gca, 'FontSize',18); grid on; drawnow
                subplot(1,3,3);plot(1:index_plot,MER_std_vals(1:index_plot),'o','LineWidth',2);xlabel(['iterations \times',num2str(store_vals)]); ylabel('std MER [dB]'); set(gca, 'FontSize',18); grid on; drawnow
            end
            disp(['Complete: ',num2str(index_direct_search/N_it*100,'%.2f'),' %'])
            index_plot=index_plot+1;
        end
    
        % save workspace 
        if mod (index_direct_search, n_save)==0
            savename=sprintf('progression_directSearch_5_targetModes_%d_iterations',index_direct_search);
            save(savename,'-v7.3');
        end
        
        % CHECK COST FUNCTION
        if MER_avg_it>MER_avg &&  MER_std_it<=max(MER_std,1) && IL_it>IL_lim
            MER_avg=MER_avg_it;
            MER_std=MER_std_it;
            IL=IL_it;
            N_changes=N_changes+1;
        else
            MASKS(mask_nr,pixel_nr_1,pixel_nr_2)=single(exp(1i*phase_val_default));
        end
    end
end

%% Plot results
if flag_plot == 1
    figure;
        subplot(1,2,1); imagesc(CM_log); colorbar; caxis([-40 0]); title('Multiplexer coupling matrix final [dB]');
        subplot(1,2,2); imagesc(CM_log_default); colorbar; caxis([-40 0]); title('Multiplexer cupling matrix default [dB]')
    
    if flag_produce_new_results==0
        figure;
            subplot(1,3,1);plot(1:N_plot,MER_vals,'o','LineWidth',2);xlabel(['iterations \times',num2str(store_vals)]); ylabel('MER [dB]'); set(gca, 'FontSize',18); grid on; drawnow
            subplot(1,3,2);plot(1:N_plot,IL_vals,'o','LineWidth',2);xlabel(['iterations \times',num2str(store_vals)]); ylabel('IL [dB]'); set(gca, 'FontSize',18); grid on; drawnow
            subplot(1,3,3);plot(1:N_plot,MER_std_vals,'o','LineWidth',2);xlabel(['iterations \times',num2str(store_vals)]); ylabel('std MER [dB]'); set(gca, 'FontSize',18); grid on; drawnow
    end

    disp(['average MER=',mat2str(mean(MER)),' dB']);

    MER_improvement=MER_vals(end)-MER_vals(1);
    disp(['Improvement on average MER=',mat2str(MER_improvement),' dB']);
    disp(['average IL=',mat2str(IL_it),' dB']);
    
    figure;
    N_plot_result=N_targetModes;
    for i=1:N_plot_result
        subplot(2,N_plot_result,i); imagesc(abs(squeeze(result(i,:,:))));axis image; title(strcat(['Multiplexer output ', num2str(i)]));  
        subplot(2,N_plot_result,N_plot_result+i); imagesc(abs(squeeze(lp_modes(target_modes(i),:,:))));axis image;title(strcat(['Desired output ', num2str(i)]));  
    end
end

% done. Cheers.
