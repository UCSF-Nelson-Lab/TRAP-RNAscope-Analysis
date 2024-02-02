%RNAscopeAnalysis
%Code to load excel data exported from QuPath and analyze in Matlab for
%automated quantification and summary.

%See Ryan et al., 2024 - Figure 4L,M; Figure S5D

%Written by Michael B. Ryan - October 16th, 2023
%% Load QuPath data into struct
%Set path for folder with all animal data - each animal should be its own
%folder. Inside that folder, each slide should have its own excel sheet (one
%for each cell population of interest, ie Slide1_1_D1 refers to Slide 1,
%the first slice, quantification for D1+ cells).

cd('D:/TRAP_RNAscope');
folders = dir; %list folders
RNAdata = [];
%Compile all animal data in struct called RNAdata
for ii = 1:length(folders) %loop through folders
    if strfind(folders(ii).name,'an2') %only analyze folders with animal data in them
        animal_name = folders(ii).name; %define animal name
        cd(sprintf('D:/TRAP_RNAscope/%s', animal_name)) %open animal folder
        files = dir; %files from this animal
        for i = 1:length(files) %loop through files
            if strfind(files(i).name,'Slide') %only analyze files for individual slides
                slide_name = files(i).name;
                slide_data = readcell(sprintf('D:/TRAP_RNAscope/%s/%s', animal_name, slide_name)); %read the data in matlab
                RNAdata.(animal_name).(slide_name(1:end-4)).rawdata = slide_data; %import slide data into cell array
            end
        end
    end
end

%% Normalize and quantify D1-channel and D2-channel for all dMSNs (for all striatum)
animals = fieldnames(RNAdata); %list animals
Summary = [];

for i = 1:length(animals) %loop through mice
    animal_TRAP.D1 = []; %initialize variables
    animal_TRAP.D2 = [];
    animal_unTRAP.D1 = [];
    animal_unTRAP.D2 = [];
    animal_name = animals{i};
    slides = fieldnames(RNAdata.(animal_name)); %create a list of the slices for this animal
    for ii = 1:length(slides) %loop through slides
        slide = slides{ii};
        if strfind(slide, 'D1') %analyze only D1+ cells *Important: you will have to add another loop to do this if you were interested in more than one classification probe. Here we are only analyzing dMSNs (aka D1R+ cells)

            data = RNAdata.(animal_name).(slide).rawdata; %create a data variable with this slice's data
            TRAP_col = find(strcmp(data(1,:),'Name')); %find the column that identifies TRAP vs unTRAP
            Nuc_area_col = find(strcmp(data(1,:),'Nucleus: Area')); %find the column that identifies Nucleus measurements
            Nuc_peri_col = find(strcmp(data(1,:),'Nucleus: Perimeter'));
            AF488_mean_col = find(strcmp(data(1,:),'Nucleus: AF488 mean')); %find column with AF488 (D1) measurements
            AF488_sum_col = find(strcmp(data(1,:),'Nucleus: AF488 sum'));
            AF647_mean_col = find(strcmp(data(1,:),'Nucleus: AF647 mean'));%find column with AF647 (D2) measurements
            AF647_sum_col = find(strcmp(data(1,:),'Nucleus: AF647 sum'));
            TRAP_index = find(strcmp(data(:,TRAP_col),'TRAP')); %index rows with TRAP cells
            unTRAP_index = find(strcmp(data(:,TRAP_col),'unTRAP')); %index rows with unTRAP cells
            
            for i = 1:2 %loop through to quantify D1R and D2R expression for all D1R+ cells (dMSNs)
                TRAP_avg = []; %initialize varibles
                unTRAP_avg = [];
                if i == 1
                    channel_col = AF488_sum_col; %D1 probe channel
                    channel = 'D1';
                else
                    channel_col = AF647_sum_col; %D2 probe channel
                    channel = 'D2';
                    
                end
                %normalize by taking the intensity sum and dividing by nucleus area
                TRAP_area_norm = [data{TRAP_index,channel_col}]./[data{TRAP_index,Nuc_area_col}];
                unTRAP_area_norm = [data{unTRAP_index,channel_col}]./[data{unTRAP_index,Nuc_area_col}];
                %normalize to the mean fluorescence of the slice
                all_cell_mean = mean([TRAP_area_norm,unTRAP_area_norm]);
                mean_normalized_TRAP = TRAP_area_norm/all_cell_mean;
                mean_normalized_unTRAP = unTRAP_area_norm/all_cell_mean;
                %calculate average for each slice
                TRAP_avg = [mean(mean_normalized_TRAP)];
                unTRAP_avg = [mean(mean_normalized_unTRAP)];
                %add slice average to array for animal
                animal_TRAP.(channel) = [animal_TRAP.(channel),TRAP_avg];
                animal_unTRAP.(channel) = [animal_unTRAP.(channel),unTRAP_avg];
                
                %save normalized data for each slice
                RNAdata.(animal_name).(slide).normalized_data.(channel).TRAP = mean_normalized_TRAP;
                RNAdata.(animal_name).(slide).normalized_data.(channel).unTRAP = mean_normalized_unTRAP;
                %save average for TRAP and unTRAP for each slice
                RNAdata.(animal_name).(slide).TRAP_average.(channel) = [TRAP_avg];
                RNAdata.(animal_name).(slide).unTRAP_average.(channel) = [unTRAP_avg];
                %save data for each animal in new Summary struct (will use for later analysis).
                Summary.(animal_name).(channel).TRAP = animal_TRAP.(channel);
                Summary.(animal_name).(channel).TRAP_avg = mean(animal_TRAP.(channel));
                Summary.(animal_name).(channel).unTRAP = animal_unTRAP.(channel);
                Summary.(animal_name).(channel).unTRAP_avg = mean(animal_unTRAP.(channel));
            end
        end
    end
end

%% Plot data for for all striatum
animals = fieldnames(Summary); %create list of animal names
channels = {'D1';'D2'}; %create labels for the probes
all_animal_summary_TRAP_D1 = []; %initialize variables
all_animal_summary_unTRAP_D1 = [];
all_animal_summary_TRAP_D2 = [];
all_animal_summary_unTRAP_D2 = [];

%Plot all mice on one graph (separate graphs for D1 and D2 probes)
for i = 1:2 %loop through D1 and D2 probes
    figure; %new figure for each probe
    for a = 1:length(animals) %loop through mice
        for ii = 1:length(Summary.(animals{a}).(channels{i}).TRAP) %loop through slices
            %plot each slice
            plot([1,2],[Summary.(animals{a}).(channels{i}).TRAP(ii),Summary.(animals{a}).(channels{i}).unTRAP(ii)],'k-o') %Relative expression of probe for TRAP vs unTRAP D1R+ cells
            hold on;
            title(sprintf('%s Summary',(channels{i})))
            %set axis for D1 or D2 probes and save this slice in the appropriate summary array
            if i == 1 %D1R expression
                ylim([0.95, 1.2])
                all_animal_summary_TRAP_D1 = [all_animal_summary_TRAP_D1, Summary.(animals{a}).(channels{i}).TRAP(ii)];
                all_animal_summary_unTRAP_D1 = [all_animal_summary_unTRAP_D1, Summary.(animals{a}).(channels{i}).unTRAP(ii)];
            else %D2R expression
                ylim([0.8,1.3])
                all_animal_summary_TRAP_D2 = [all_animal_summary_TRAP_D2, Summary.(animals{a}).(channels{i}).TRAP(ii)];
                all_animal_summary_unTRAP_D2 = [all_animal_summary_unTRAP_D2, Summary.(animals{a}).(channels{i}).unTRAP(ii)];
            end
        end
    end
end

%% Plot bar graphs for average +/- SEM for each channel (D1R and D2R)
%for D1R expression
x = [1,2];
y = [mean(all_animal_summary_TRAP_D1),mean(all_animal_summary_unTRAP_D1)];
yerr = [[std(all_animal_summary_TRAP_D1)/sqrt(length(all_animal_summary_TRAP_D1))],[std(all_animal_summary_unTRAP_D1)/sqrt(length(all_animal_summary_unTRAP_D1))]];
bar(x,y)
hold on
errorbar(x,y,yerr,yerr)
title('Average Relative D1R expression')
xticklabels({'TRAP dMSN','unTRAP dMSN'})
ylim([0.95,1.2])
hold off

%for D2R expression
x = [1,2];
y = [mean(all_animal_summary_TRAP_D2),mean(all_animal_summary_unTRAP_D2)];
yerr = [[std(all_animal_summary_TRAP_D2)/sqrt(length(all_animal_summary_TRAP_D2))],[std(all_animal_summary_unTRAP_D2)/sqrt(length(all_animal_summary_unTRAP_D2))]];
figure;
bar(x,y)
hold on
errorbar(x,y,yerr,yerr)
title('Average Relative D2R expression')
xticklabels({'TRAP dMSN','unTRAP dMSN'})
ylim([0.8,1.3])
hold off

%% Normalize and quantify D1-channel and D2-channel for all dMSNs (for DLS, DMS, VLS separately)
animals = fieldnames(RNAdata); %list animals
for i = 1:length(animals) %loop through mice
    animal_TRAP.DLS.D1 = []; %initialize variables
    animal_TRAP.DMS.D1 = [];
    animal_TRAP.VLS.D1 = [];
    animal_TRAP.DLS.D2 = [];
    animal_TRAP.DMS.D2 = [];
    animal_TRAP.VLS.D2 = [];
    animal_unTRAP.DLS.D1 = [];
    animal_unTRAP.DMS.D1 = [];
    animal_unTRAP.VLS.D1 = [];
    animal_unTRAP.DLS.D2 = [];
    animal_unTRAP.DMS.D2 = [];
    animal_unTRAP.VLS.D2 = [];
    animal_name = animals{i};
    slides = fieldnames(RNAdata.(animal_name)); %list slides
    for ii = 1:length(slides) %loop through slides
        slide = slides{ii};
        for iii = 1:3 %loop through DLS, DMS, VLS subregions
            if iii == 1
                region = 'DLS'; %create a label for each region
            elseif iii == 2
                region = 'DMS';
            else
                region = 'VLS';
            end
            if strfind(slide, 'D1') %analyze only D1+ cells
                
                data = RNAdata.(animal_name).(slide).rawdata;
                TRAP_col = find(strcmp(data(1,:),'Name')); %find the column that identifies TRAP vs unTRAP
                Nuc_area_col = find(strcmp(data(1,:),'Nucleus: Area')); %find the column that identifies Nucleus parameters
                Nuc_peri_col = find(strcmp(data(1,:),'Nucleus: Perimeter'));
                AF488_mean_col = find(strcmp(data(1,:),'Nucleus: AF488 mean')); %find the column that identifies AF488 (D1) parameters
                AF488_sum_col = find(strcmp(data(1,:),'Nucleus: AF488 sum'));
                AF647_mean_col = find(strcmp(data(1,:),'Nucleus: AF647 mean'));%find the column that identifies AF647 (D2) parameters
                AF647_sum_col = find(strcmp(data(1,:),'Nucleus: AF647 sum'));
                subregion_col = find(strcmp(data(1,:),'Parent'));%find the column that identifies subregion
                
                TRAP_index = find(strcmp(data(:,TRAP_col),'TRAP')); %index all TRAP cells
                unTRAP_index = find(strcmp(data(:,TRAP_col),'unTRAP')); %index all unTRAP cells
                region_index = find(strcmp(data(:,subregion_col),region)); %index subregion cells
                
                for i = 1:2 %loop through to quantify D1 probe and D2 probe for all cells
                    TRAP_avg = [];
                    unTRAP_avg = [];
                    if i == 1
                        channel_col = AF488_sum_col;
                        channel = 'D1';
                    else
                        channel_col = AF647_sum_col;
                        channel = 'D2';
                        
                    end
                    %find region specific cells
                    TRAP_region_index = intersect(TRAP_index, region_index);
                    unTRAP_region_index = intersect(unTRAP_index, region_index);
                    %normalize by taking the intensity sum and dividing by nucleus area
                    TRAP_area_norm = [data{TRAP_region_index,channel_col}]./[data{TRAP_region_index,Nuc_area_col}];
                    unTRAP_area_norm = [data{unTRAP_region_index,channel_col}]./[data{unTRAP_region_index,Nuc_area_col}];
                    %normalize to the mean fluorescence for the slice
                    all_cell_mean = mean([TRAP_area_norm,unTRAP_area_norm]);
                    mean_normalized_TRAP = TRAP_area_norm/all_cell_mean;
                    mean_normalized_unTRAP = unTRAP_area_norm/all_cell_mean;
                    %calculate average for each slice
                    TRAP_avg = [mean(mean_normalized_TRAP)];
                    unTRAP_avg = [mean(mean_normalized_unTRAP)];
                    %add slice average to array for animal
                    animal_TRAP.(region).(channel) = [animal_TRAP.(region).(channel),TRAP_avg];
                    animal_unTRAP.(region).(channel) = [animal_unTRAP.(region).(channel),unTRAP_avg];
                    
                    %save normalized data for each slide
                    RNAdata.(animal_name).(slide).normalized_data.(channel).(region).TRAP = mean_normalized_TRAP;
                    RNAdata.(animal_name).(slide).normalized_data.(channel).(region).unTRAP = mean_normalized_unTRAP;
                    %save average for TRAP and unTRAP for each slice
                    RNAdata.(animal_name).(slide).(region).TRAP_average.(channel) = [TRAP_avg];
                    RNAdata.(animal_name).(slide).(region).unTRAP_average.(channel) = [unTRAP_avg];
                    %save data each animal
                    Summary.(animal_name).(channel).(region).TRAP = animal_TRAP.(region).(channel);
                    Summary.(animal_name).(channel).(region).TRAP_avg = mean(animal_TRAP.(region).(channel));
                    Summary.(animal_name).(channel).(region).unTRAP = animal_unTRAP.(region).(channel);
                    Summary.(animal_name).(channel).(region).unTRAP_avg = mean(animal_unTRAP.(region).(channel));
                end
            end
        end
    end
end
%% Plot data for for striatal subregions
animals = fieldnames(Summary); %list animals
channels = {'D1';'D2'};
regions = {'DLS','DMS','VLS'};
all_animal_summary_TRAP_D1.DLS = [];
all_animal_summary_TRAP_D1.DMS = [];
all_animal_summary_TRAP_D1.VLS = [];
all_animal_summary_unTRAP_D1.DLS = [];
all_animal_summary_unTRAP_D1.DMS = [];
all_animal_summary_unTRAP_D1.VLS = [];
all_animal_summary_TRAP_D2.DLS = [];
all_animal_summary_TRAP_D2.DMS = [];
all_animal_summary_TRAP_D2.VLS = [];
all_animal_summary_unTRAP_D2.DLS = [];
all_animal_summary_unTRAP_D2.DMS = [];
all_animal_summary_unTRAP_D2.VLS = [];

for s = 1:3 %loop through subregions
    for i = 1:2 %loop through D1 and D2 probes
        figure;
        for a = 1:length(animals) %loop through mice
            for ii = 1:length(Summary.(animals{a}).(channels{i}).(regions{s}).TRAP)
                plot([1,2],[Summary.(animals{a}).(channels{i}).(regions{s}).TRAP(ii),Summary.(animals{a}).(channels{i}).(regions{s}).unTRAP(ii)],'k-o');
                hold on;
                
                title(sprintf('%s %s Summary', (regions{s}), (channels{i})))
                if i == 1
                    ylim([0.85, 1.3])
                    all_animal_summary_TRAP_D1.(regions{s}) = [all_animal_summary_TRAP_D1.(regions{s}), Summary.(animals{a}).(channels{i}).(regions{s}).TRAP(ii)];
                    all_animal_summary_unTRAP_D1.(regions{s}) = [all_animal_summary_unTRAP_D1.(regions{s}), Summary.(animals{a}).(channels{i}).(regions{s}).unTRAP(ii)];
                else
                    ylim([0.6,2.2])
                    all_animal_summary_TRAP_D2.(regions{s}) = [all_animal_summary_TRAP_D2.(regions{s}), Summary.(animals{a}).(channels{i}).(regions{s}).TRAP(ii)];
                    all_animal_summary_unTRAP_D2.(regions{s}) = [all_animal_summary_unTRAP_D2.(regions{s}), Summary.(animals{a}).(channels{i}).(regions{s}).unTRAP(ii)];
                end
            end
        end
        hold off
    end
end

%% Plot average +/- SEM for each striatal subregion for D1 probe separate for TRAP and unTRAP
figure;
x = [1,2,3];
DLS = all_animal_summary_TRAP_D1.DLS;
DMS = all_animal_summary_TRAP_D1.DMS;
VLS = all_animal_summary_TRAP_D1.VLS;
y = [mean(DLS),mean(DMS), mean(VLS)];
yerr = [[std(DLS)/sqrt(length(DLS))],[std(DMS)/sqrt(length(DMS))], [std(VLS)/sqrt(length(VLS))]];
bar(x,y)
hold on
errorbar(x,y,yerr,yerr)
title(sprintf('Average Subregion TRAP D1 expression'))
xticklabels({'DLS','DMS','VLS'})
ylim([0.9,1.1])
hold off


figure;
x = [1,2,3];
DLS = all_animal_summary_unTRAP_D1.DLS;
DMS = all_animal_summary_unTRAP_D1.DMS;
VLS = all_animal_summary_unTRAP_D1.VLS;
y = [mean(DLS),mean(DMS), mean(VLS)];
yerr = [[std(DLS)/sqrt(length(DLS))],[std(DMS)/sqrt(length(DMS))], [std(VLS)/sqrt(length(VLS))]];
bar(x,y)
hold on
errorbar(x,y,yerr,yerr)
title(sprintf('Average Subregion unTRAP D1 expression'))
xticklabels({'DLS','DMS','VLS'})
ylim([0.9,1.1])
hold off





