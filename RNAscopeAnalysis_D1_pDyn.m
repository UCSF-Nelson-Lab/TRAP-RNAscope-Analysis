%RNAscopeAnalysis
%Code to load excel data exported from QuPath and analyze in Matlab for
%automated quantification and summary.

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
    if strfind(folders(ii).name,'an') %only analyze folders with animal data in them
        animal_name = folders(ii).name; %define animal name
        cd(sprintf('D:/TRAP_RNAscope/%s', animal_name)) %open animal folder
        files = dir; %files from this animal
        for i = 1:length(files) %loop through files
            if strfind(files(i).name,'pDyn') %only analyze files for individual slides with pDyn probes
                slide_name = files(i).name;
                slide_data = readcell(sprintf('D:/TRAP_RNAscope/%s/%s', animal_name, slide_name)); %read the data in matlab
                pDynRNAdata.(animal_name).(slide_name(1:end-4)).rawdata = slide_data; %import slide data into cell array
            end
        end
    end
end

%% Normalize and quantify pDyn expression for all dMSNs
animals = fieldnames(pDynRNAdata); %list animals
Summary = [];
for i = 1:length(animals) %loop through mice
    animal_TRAP.D1.pDyn = []; %initialize variables
    animal_unTRAP.D1.pDyn = [];
    
    animal_name = animals{i};
    slides = fieldnames(pDynRNAdata.(animal_name)); %create a list of the slices for this animal
    for ii = 1:length(slides) %loop through slides
        slide = slides{ii};
        if strfind(slide, 'D1_pDyn') %analyze D1+ cells
            
            data = pDynRNAdata.(animal_name).(slide).rawdata; %create a data variable with this slice's data
            TRAP_col = find(strcmp(data(1,:),'Name')); %find the column that identifies TRAP vs unTRAP
            Nuc_area_col = find(strcmp(data(1,:),'Nucleus: Area')); %find the column that identifies Nucleus measurements
            Nuc_peri_col = find(strcmp(data(1,:),'Nucleus: Perimeter'));
            AF488_mean_col = find(strcmp(data(1,:),'Nucleus: AF488 mean')); %find column with AF488 (D1) measurements
            AF488_sum_col = find(strcmp(data(1,:),'Nucleus: AF488 sum'));
            AF647_mean_col = find(strcmp(data(1,:),'Nucleus: AF647 mean'));%find column with AF647 (pDyn) measurements
            AF647_sum_col = find(strcmp(data(1,:),'Nucleus: AF647 sum'));
            TRAP_index = find(strcmp(data(:,TRAP_col),'TRAP')); %index rows with TRAP cells
            unTRAP_index = find(strcmp(data(:,TRAP_col),'unTRAP')); %index rows with unTRAP cells
            
            TRAP_avg = []; %initialize varibles
            unTRAP_avg = [];
            
            channel_col = AF647_sum_col; %D1 probe channel
            channel = 'D1';
            probe = 'pDyn';
            
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
            animal_TRAP.(channel).(probe) = [animal_TRAP.(channel).(probe),TRAP_avg];
            animal_unTRAP.(channel).(probe) = [animal_unTRAP.(channel).(probe),unTRAP_avg];
            
            %save normalized data for each slice
            pDynRNAdata.(animal_name).(slide).normalized_data.(channel).(probe).TRAP = mean_normalized_TRAP;
            pDynRNAdata.(animal_name).(slide).normalized_data.(channel).(probe).unTRAP = mean_normalized_unTRAP;
            %save average for TRAP and unTRAP for each slice
            pDynRNAdata.(animal_name).(slide).TRAP_average.(channel).(probe) = [TRAP_avg];
            RNAdata.(animal_name).(slide).unTRAP_average.(channel).(probe) = [unTRAP_avg];
            %save data for each animal in new Summary struct (will use for later analysis)
            pDynSummary.(animal_name).(channel).(probe).TRAP = animal_TRAP.(channel).(probe);
            pDynSummary.(animal_name).(channel).(probe).TRAP_avg = mean(animal_TRAP.(channel).(probe));
            pDynSummary.(animal_name).(channel).(probe).unTRAP = animal_unTRAP.(channel).(probe);
            pDynSummary.(animal_name).(channel).(probe).unTRAP_avg = mean(animal_unTRAP.(channel).(probe));
        end
    end
end
%% Plot pDyn data for all D1R+ cells for each slice
animals = fieldnames(pDynSummary); %create list of animal names
all_animal_summary_TRAP_D1_pDyn = []; %initialize variables
all_animal_summary_unTRAP_D1_pDyn = [];

%Plot all slices for pDyn expression TRAP vs unTRAP
figure;
for a = 1:length(animals) %loop through mice   
    for ii = 1:length(pDynSummary.(animals{a}).D1.pDyn.TRAP) %loop through slices
        %plot each slice
        plot([1,2],[pDynSummary.(animals{a}).D1.pDyn.TRAP(ii),pDynSummary.(animals{a}).D1.pDyn.unTRAP(ii)],'k-o')
        hold on;
        title(sprintf('dMSN pDyn Summary'))
        all_animal_summary_TRAP_D1_pDyn = [all_animal_summary_TRAP_D1_pDyn, pDynSummary.(animals{a}).D1.pDyn.TRAP(ii)];
        all_animal_summary_unTRAP_D1_pDyn = [all_animal_summary_unTRAP_D1_pDyn, pDynSummary.(animals{a}).D1.pDyn.unTRAP(ii)];
    end
end

%% Plot bar graphs for average +/- SEM for each channel for all of striatum
%for pDyn
figure;
x = [1,2];
y = [mean(all_animal_summary_TRAP_D1_pDyn),mean(all_animal_summary_unTRAP_D1_pDyn)];
yerr = [[std(all_animal_summary_TRAP_D1_pDyn)/sqrt(length(all_animal_summary_TRAP_D1_pDyn))],[std(all_animal_summary_unTRAP_D1_pDyn)/sqrt(length(all_animal_summary_unTRAP_D1_pDyn))]];
bar(x,y)
hold on
errorbar(x,y,yerr,yerr)
title('Average Relative pDyn expression')
xticklabels({'TRAP dMSN','unTRAP dMSN'})
ylim([0.8,2.2])
hold off







