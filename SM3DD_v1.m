clc
clearvars
fclose('all');
close all
delete(gcp('nocreate'))
warning('off')

% Check if SM3DD exists
if ~isfolder('.\SM3DD')
    ErrorMessage = "SM3DD folder not detected at root";
    writematrix(ErrorMessage,'.\ErrorMessage_SM3DD folder not detected at root.txt')
    return
end

% Check if SM3DD\Data exists
if ~isfolder('.\SM3DD\Data')
    ErrorMessage = "SM3DD\Data folder not detected at root";
    writematrix(ErrorMessage,'.\ErrorMessage_Data folder not detected within SM3DD.txt')
    return
end

if ~isfile('.\SM3DD\Platform.csv')
    ErrorMessage = ".\SM3DD\Platform.csv does not exist";
    writematrix(ErrorMessage,'.\ErrorMessage_Platform csv file does not exist.txt')
    return
else
    PlatformTable = readtable('.\SM3DD\Platform.csv');
    PlatformAssay = logical(PlatformTable.Platform_Assay_Binary);
    if sum(PlatformAssay) > 1
        ErrorMessage = "More than one Platform Assay flagged as TRUE (1)";
        ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
        writematrix(ErrorMessage,ErrorMessage_file)
        return
    elseif ~any(PlatformAssay)
        ErrorMessage = "No Platform Assay flagged as TRUE (1)";
        ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
        writematrix(ErrorMessage,ErrorMessage_file)
        return
    end
    Z_factor = single(PlatformTable.Z_factor(PlatformAssay));
    if ~(Z_factor > 0)
        ErrorMessage = "Z_factor not numeric";
        ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
        writematrix(ErrorMessage,ErrorMessage_file)
        return
    end
    NumSections = single(PlatformTable.Number_of_Sections(PlatformAssay));
    if ~(NumSections > 0) || ~(floor(NumSections) == NumSections)
        ErrorMessage = "Number_of_Sections not positive integar";
        ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
        writematrix(ErrorMessage,ErrorMessage_file)
        return
    end
    [~,PlatformAssayTableColumns] = size(PlatformTable);
    NumLeadColumns = 11;
    NumControlProbeNameColumns = PlatformAssayTableColumns - NumLeadColumns;
    ControlProbeNames = strings(NumControlProbeNameColumns,1);
    for ControlProbe = 1 : NumControlProbeNameColumns
        ControlProbeNames(ControlProbe) = string(table2cell(PlatformTable(PlatformAssay,NumLeadColumns+ControlProbe)));
    end
    ControlProbeNames(ControlProbeNames == "") = [];
    NumControlProbeNames = length(ControlProbeNames);
    NumCharacters_in_ControlProbeNames = unique(strlength(ControlProbeNames));
    if length(NumCharacters_in_ControlProbeNames) > 1
        ErrorMessage = "Inconsistent number of characters in StartOfControlProbeNames";
        ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
        writematrix(ErrorMessage,ErrorMessage_file)
        return
    end
    ProbeCoordinatesFileExtension = string(PlatformTable.Probe_Coordinate_file_extension(PlatformAssay));
    ExpressionMatrixFileExtension = string(PlatformTable.Expression_Matrix_file_extension(PlatformAssay));
    FirstColumn = PlatformTable.First_Column_with_ProbeName_along_top_row(PlatformAssay);
end

if isfile('.\SM3DD\Probes.mat')
    load('.\SM3DD\Probes.mat')
else
    % Generate ProbeName list from 'first' expression matrix file
    SearchFor = strcat('.\SM3DD\Data\**\*',ExpressionMatrixFileExtension);
    exprMat_files = dir(SearchFor);
    if isempty(exprMat_files)
        ErrorMessage = "Expression matrix file not found";
        ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
        writematrix(ErrorMessage,ErrorMessage_file)
        return
    end
    exprMat_file = strcat(exprMat_files(1).folder,'\',exprMat_files(1).name);
    opts = detectImportOptions(exprMat_file);
    opts.DataLines = [1,1];
    ProbeNames = transpose(readmatrix(exprMat_file,opts,'OutputType','char'));
    % Delete fov & cell_ID
    if FirstColumn > 1
        ProbeNames(1:FirstColumn-1,:) = [];
    end

    % Find NONcontrolProbes
    if NumControlProbeNames > 0
        NONControlProbes_binary = ones(length(ProbeNames),1,'logical');
        for P = 1 : length(ProbeNames)
            ThisProbe =  ProbeNames{P};
            ThisProbeStartsWith = string(ThisProbe(1:min(NumCharacters_in_ControlProbeNames,length(ThisProbe))));
            if sum(ControlProbeNames == ThisProbeStartsWith)
                NONControlProbes_binary(P) = false;
            end
        end
        NONControlProbes = sort(ProbeNames(NONControlProbes_binary));
        ControlProbes = sort(ProbeNames(~NONControlProbes_binary));
    else
        NONControlProbes = sort(ProbeNames);
        ControlProbes = "";
    end
    NumNONControlProbes = length(NONControlProbes);
    NumberOfProbes = length(ProbeNames);
    save('.\SM3DD\Probes.mat',"NumberOfProbes","ProbeNames","NumNONControlProbes","NONControlProbes","ControlProbes")
    writecell(ControlProbes,'.\SM3DD\ControlProbes.csv')
    writecell(NONControlProbes,'.\SM3DD\NONControlProbes.csv')
end

SectionSize = ceil(NumNONControlProbes / NumSections);

if ~isfile('.\SM3DD\Samples to FOVs.csv')
    ErrorMessage = "Samples to FOVs.csv not found";
    ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
    writematrix(ErrorMessage,ErrorMessage_file)
    return
else
    opts = detectImportOptions('.\SM3DD\Samples to FOVs.csv');
    opts.SelectedVariableNames = 1;
    TMAs = readmatrix('.\SM3DD\Samples to FOVs.csv',opts,'OutputType','string');
    opts.SelectedVariableNames = 2;
    FOVs_on_TMAs = readmatrix('.\SM3DD\Samples to FOVs.csv',opts,'OutputType','single');
    opts.SelectedVariableNames = 3;
    Matching_Sample_IDs = readmatrix('.\SM3DD\Samples to FOVs.csv',opts,'OutputType','string');
end

if ~isfile('.\SM3DD\Sample Groups.csv')
    ErrorMessage = "Sample Groups.csv not found";
    ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
    writematrix(ErrorMessage,ErrorMessage_file)
    return
else
    opts = detectImportOptions('.\SM3DD\Sample Groups.csv');
    opts.DataLines = [2,Inf];
    opts.SelectedVariableNames = 1;
    Sample_Group_1 = readmatrix('.\SM3DD\Sample Groups.csv',opts,'OutputType','string');
    opts.SelectedVariableNames = 2;
    Sample_Group_2 = readmatrix('.\SM3DD\Sample Groups.csv',opts,'OutputType','string');
    opts.SelectedVariableNames = 5;
    opts.DataLines = [1,1];
    Sample_Group_1_name = readmatrix('.\SM3DD\Sample Groups.csv',opts,'OutputType','string');
    opts.DataLines = [2,2];
    Sample_Group_2_name = readmatrix('.\SM3DD\Sample Groups.csv',opts,'OutputType','string');
    opts.DataLines = [3,3];
    Paired_comparison_binary = logical(readmatrix('.\SM3DD\Sample Groups.csv',opts,'OutputType','single'));
end

% Generate lists FOVs on TMAs to process
Sample_Group_1(ismissing(Sample_Group_1)) = [];
Sample_Group_2(ismissing(Sample_Group_2)) = [];
Samples_to_process = [Sample_Group_1;Sample_Group_2];
Num_Sample_Group_1 = length(Sample_Group_1);
Num_Sample_Group_2 = length(Sample_Group_2);
Num_Samples_to_process = Num_Sample_Group_1 + Num_Sample_Group_2;
Num_FOVs_per_Sample = zeros(Num_Samples_to_process,1);
for Sample = 1 : Num_Samples_to_process
    Num_FOVs_per_Sample(Sample) = sum(Matching_Sample_IDs == Samples_to_process(Sample));
end
if sum(Num_FOVs_per_Sample == 0) > 0
    NotFound = Samples_to_process(Num_FOVs_per_Sample == 0);
    ErrorMessage = ["Sample(s) listed in Sample Groups.csv not found in Samples to FOVs.csv:";NotFound];
    writematrix(ErrorMessage,'.\ErrorMessage_Sample_Group_mismatch.txt')
    return
end
Num_FOVs_to_process = sum(Num_FOVs_per_Sample);
TMAs_for_FOVs_to_process = strings(Num_FOVs_to_process,1);
FOVs_for_FOVs_to_process = zeros(Num_FOVs_to_process,1,'single');
count = 0;
for Sample = 1 : Num_Samples_to_process
    Matches = Matching_Sample_IDs == Samples_to_process(Sample);
    TMAs_for_FOVs_to_process(count+1:count+Num_FOVs_per_Sample(Sample)) = TMAs(Matches);
    FOVs_for_FOVs_to_process(count+1:count+Num_FOVs_per_Sample(Sample)) = FOVs_on_TMAs(Matches);
    count = count + Num_FOVs_per_Sample(Sample);
end
UniqueTMAs = unique(TMAs_for_FOVs_to_process);
Num_UniqueTMAs = length(UniqueTMAs);

FOVRangesFolderName = '.\SM3DD\FOVRanges';
if ~isfolder(FOVRangesFolderName)
    mkdir(FOVRangesFolderName)
end

% Generate lists of FOV Ranges
for TMA = 1 : Num_UniqueTMAs
    Filename_and_Folder = strcat('.\SM3DD\Data\**\',UniqueTMAs(TMA),ProbeCoordinatesFileExtension);
    TX_file_dir = dir(Filename_and_Folder);
    if isempty(TX_file_dir)
        ErrorMessage = strcat(Filename_and_Folder,' not found');
        writematrix(ErrorMessage,'.\ErrorMessage_ProbeCoordinatesFile not found.txt')
        return
    end
    TX_file = strcat(TX_file_dir.folder,'\',TX_file_dir.name);
    opts = detectImportOptions(TX_file);
    FOVRanges_file = strcat(FOVRangesFolderName,'\',UniqueTMAs(TMA),'.mat');
    if ~isfile(FOVRanges_file)
        opts.SelectedVariableNames = 1; % FOV
        TX_FOVs = readmatrix(TX_file,opts);
        UniqueFOVs = unique(TX_FOVs);
        Max_UniqueFOVs = max(UniqueFOVs);
        MaxDigitsInFOV = ceil(log10((Max_UniqueFOVs + 1)));
        FOVRanges = nan(Max_UniqueFOVs,2);
        for FOV = 1 : Max_UniqueFOVs
            First = find(TX_FOVs == FOV,1,"first");
            if ~isempty(First)
                FOVRanges(FOV,1) = First;
                FOVRanges(FOV,2) = find(TX_FOVs == FOV,1,"last");
            end
        end
        clear TX_FOVs
        FOVRanges = FOVRanges + 1; % opt.DataLines goes from [2,Inf]
        save(FOVRanges_file,"FOVRanges","UniqueFOVs","Max_UniqueFOVs","MaxDigitsInFOV","-v7.3","-nocompression")
    end
end

PCA_options_file = '.\SM3DD\PCA options.csv';
if ~isfile(PCA_options_file)
    ErrorMessage = "No PCA options file";
    ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
    writematrix(ErrorMessage,ErrorMessage_file)
    return
else
    PCA_options_table = readtable(PCA_options_file);
end

StandardisedM3DDFolderName = '.\SM3DD\StandardisedM3DD';
if ~isfolder(StandardisedM3DDFolderName)
    mkdir(StandardisedM3DDFolderName)
end

FOV_level_means_FolderName = '.\SM3DD\FOV_level_means';
if ~isfolder(FOV_level_means_FolderName)
    mkdir(FOV_level_means_FolderName)
end

FOV_level_ProbeCounts_FolderName = '.\SM3DD\FOV_level_ProbeCounts';
if ~isfolder(FOV_level_ProbeCounts_FolderName)
    mkdir(FOV_level_ProbeCounts_FolderName)
end

MeanToProbeValues_FolderName = '.\SM3DD\MeanToProbeValues';
if ~isfolder(MeanToProbeValues_FolderName)
    mkdir(MeanToProbeValues_FolderName)
end

FOV_Completion_FolderName = '.\SM3DD\FOV_Completion';
if ~isfolder(FOV_Completion_FolderName)
    mkdir(FOV_Completion_FolderName)
end

PCA_output_FolderName = '.\SM3DD\PCA_output';
if ~isfolder(PCA_output_FolderName)
    mkdir(PCA_output_FolderName)
end

FROM_probe_PCA_diffs_FolderName = '.\SM3DD\PCA_output\FROM_probe_PCA_diffs';
if ~isfolder(FROM_probe_PCA_diffs_FolderName)
    mkdir(FROM_probe_PCA_diffs_FolderName)
end

PCA_options = logical(PCA_options_table.Binary);
if PCA_options(1)
    FROM_probe_PCA_diffs_folder_IncludeFOVs = strcat(FROM_probe_PCA_diffs_FolderName,'\IncludeFOVs_With_Absent_Transcripts');
    if ~isfolder(FROM_probe_PCA_diffs_folder_IncludeFOVs)
        mkdir(FROM_probe_PCA_diffs_folder_IncludeFOVs)
    end
end
if PCA_options(2)
    FROM_probe_PCA_diffs_folder_ExcludeFOVs = strcat(FROM_probe_PCA_diffs_FolderName,'\ExcludeFOVs_With_Absent_Transcripts');
    if ~isfolder(FROM_probe_PCA_diffs_folder_ExcludeFOVs)
        mkdir(FROM_probe_PCA_diffs_folder_ExcludeFOVs)
    end
end


% Completion Times
FOVCompletionTimes_file = '.\SM3DD\FOVCompletionTimes.csv';
FOVCompletionTimes_sections_file = '.\SM3DD\FOVCompletionTimes_sections.mat';
if ~isfile(FOVCompletionTimes_file)
    FOVCompletionTimes = zeros(Num_FOVs_to_process,1);
    FOVCompletionTimes_sections = zeros(Num_FOVs_to_process,2+NumSections); % LoadDataTime,ProbeBinaryTime,TimePerSection
    writematrix([TMAs_for_FOVs_to_process,FOVs_for_FOVs_to_process,FOVCompletionTimes],FOVCompletionTimes_file)
    save(FOVCompletionTimes_sections_file,"FOVCompletionTimes_sections","FOVCompletionTimes");
else
    load(FOVCompletionTimes_sections_file,"FOVCompletionTimes_sections","FOVCompletionTimes");
end

% calculate SM3DD and FOV level mean SM3DDs
for TMA = 1 : Num_UniqueTMAs
    FOVRanges_file = strcat(FOVRangesFolderName,'\',UniqueTMAs(TMA),'.mat');
    load(FOVRanges_file,"FOVRanges","UniqueFOVs","Max_UniqueFOVs","MaxDigitsInFOV")
    Filename_and_Folder = strcat('.\SM3DD\Data\**\',UniqueTMAs(TMA),'_tx_file.csv');
    TX_file_dir = dir(Filename_and_Folder);
    if isempty(TX_file_dir)
        ErrorMessage = strcat(Filename_and_Folder,' not found');
        writematrix(ErrorMessage,'.\ErrorMessage.txt')
        return
    end
    TX_file = strcat(TX_file_dir.folder,'\',TX_file_dir.name);
    opts = detectImportOptions(TX_file);
    FOVs_for_this_TMA = FOVs_for_FOVs_to_process(TMAs_for_FOVs_to_process == UniqueTMAs(TMA));
    Num_FOVs_for_this_TMA = length(FOVs_for_this_TMA);

    for FOV = 1 : Num_FOVs_for_this_TMA
        Digits_in_current_FOV = ceil(log10((FOVs_for_this_TMA(FOV) + 1)));
        FOV_as_string = string(FOVs_for_this_TMA(FOV));
        while Digits_in_current_FOV < MaxDigitsInFOV
            FOV_as_string = append('0',FOV_as_string);
            Digits_in_current_FOV = Digits_in_current_FOV + 1;
        end
        if isnan(FOVRanges(FOVs_for_this_TMA(FOV),1))
            ErrorMessage = strcat("No data for F",FOV_as_string,"_on TMA_",UniqueTMAs(TMA));
            ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
            writematrix(ErrorMessage,ErrorMessage_file)
            continue
        end
        FOV_Completion_file = strcat(FOV_Completion_FolderName,'\FOV_Completion_',UniqueTMAs(TMA),'_F',FOV_as_string,'.mat');
        if isfile(FOV_Completion_file)
            load(FOV_Completion_file,"FOV_Num_Sections_done")
            if FOV_Num_Sections_done == NumSections
                continue
            else
                SectionStart = FOV_Num_Sections_done + 1;
            end
            MeanToProbeValues_file = strcat(MeanToProbeValues_FolderName,'\MeanToProbeValues_',UniqueTMAs(TMA),'_','F',FOV_as_string,'.mat');
            load(MeanToProbeValues_file,"MeanToProbeValues")
            MSM3DD_file = strcat(FOV_level_means_FolderName,'\MSM3DD_',UniqueTMAs(TMA),'_','F',FOV_as_string,'.mat');
            load(MSM3DD_file,"MSM3DD")
        else
            FOV_StandardisedM3DDFolderName = strcat(StandardisedM3DDFolderName,'\',UniqueTMAs(TMA),'_','F',FOV_as_string);
            Start_P = 1;
            SectionStart = 1;
            mkdir(FOV_StandardisedM3DDFolderName)
            for P = 1 : NumNONControlProbes
                FROM_FOV_StandardisedM3DDFolderName = strcat(FOV_StandardisedM3DDFolderName,'\',NONControlProbes(P));
                mkdir(FROM_FOV_StandardisedM3DDFolderName)
            end
            MeanToProbeValues = nan(1,NumNONControlProbes,'single');
            MeanToProbeValues_file = strcat(MeanToProbeValues_FolderName,'\MeanToProbeValues_',UniqueTMAs(TMA),'_','F',FOV_as_string,'.mat');
            save(MeanToProbeValues_file,"MeanToProbeValues","-v7.3","-nocompression")
            MSM3DD_file = strcat(FOV_level_means_FolderName,'\MSM3DD_',UniqueTMAs(TMA),'_','F',FOV_as_string,'.mat');
            MSM3DD = nan(NumNONControlProbes,NumNONControlProbes,'single');
            save(MSM3DD_file,"MSM3DD","-v7.3","-nocompression")
        end

        % localX, localY, Z
        opts.SelectedVariableNames = [PlatformTable.Local_X_column(PlatformAssay),PlatformTable.Local_Y_column(PlatformAssay),PlatformTable.Local_Z_column(PlatformAssay)];

        LoadDataTimeTIC = tic;

        opts.DataLines = [FOVRanges(FOV,1),FOVRanges(FOV,2)];
        FOV_TX_data = readmatrix(TX_file,opts,'OutputType','single');

        LoadDataTimeTOC = toc(LoadDataTimeTIC);

        ProbeBinaryTIC = tic;

        % Target
        opts.SelectedVariableNames = PlatformTable.ProbeID_column(PlatformAssay);

        FOV_TX_Probes = readmatrix(TX_file,opts,'OutputType','string');
        NumFOVProbes = FOVRanges(FOV,2)-FOVRanges(FOV,1)+1;

        ProbeBinary = zeros(NumFOVProbes,NumNONControlProbes,'logical');

        % Check if there is already a parpool, make one if there isn't
        p = gcp('nocreate');
        if isempty(p)
            parpool("threads");
        end

        parfor P = 1 : NumNONControlProbes
            ProbeBinary(:,P) = FOV_TX_Probes == NONControlProbes(P);
        end
        % delete(gcp('nocreate'))
        clear FOV_TX_Probes
        ControlProbes = ~any(ProbeBinary,2);
        ProbeBinary(ControlProbes,:) = [];
        FOV_TX_data(ControlProbes,:) = [];
        [NumNONControlProbesinFOV,~] = size(FOV_TX_data);
        clear ControlProbes

        ProbeBinaryTOC = toc(ProbeBinaryTIC);

        FOV_TX_data(:,3) = FOV_TX_data(:,3) * Z_factor; % Adjust for Z scaling

        ProbeCounts_file = strcat(FOV_level_ProbeCounts_FolderName,'\ProbeCounts_',UniqueTMAs(TMA),'_','F',FOV_as_string,'.mat');
        if isfile(ProbeCounts_file)
            load(ProbeCounts_file,"ProbeCounts","ZeroProbeCounts")
        else
            ProbeCounts = transpose(sum(ProbeBinary));
            ZeroProbeCounts = ProbeCounts == 0;
            save(ProbeCounts_file,"ProbeCounts","ZeroProbeCounts")
        end

        for Section = SectionStart : NumSections
            SectionTIC = tic;
            if Section == NumSections
                ProbesInSection = NumNONControlProbes - (NumSections-1)*SectionSize;
            else
                ProbesInSection = SectionSize;
            end
            AdjustForSection = (Section - 1) * SectionSize;
            M3DD = nan(NumNONControlProbesinFOV,ProbesInSection);

            % Check if there is already a parpool, make one if there isn't
            p = gcp('nocreate');
            if isempty(p)
                parpool("threads");
            end

            for P = AdjustForSection + 1 : AdjustForSection + ProbesInSection
                if ZeroProbeCounts(P)
                    continue
                end
                TargetXYZs = FOV_TX_data(ProbeBinary(:,P),:);

                % Save as FROM probe coords for use by other scripts
                ProbeCoords_file = strcat(FOV_StandardisedM3DDFolderName,'\',NONControlProbes(P),'\FROM_ProbeCoords.mat');
                save(ProbeCoords_file,"TargetXYZs","-v7.3","-nocompression")

                [NumTargets,~] = size(TargetXYZs);
                OnlyOneTarget = NumTargets == 1;
                M3DD_column = P - AdjustForSection;
                parfor FromProbe = 1 : NumNONControlProbesinFOV
                    Distances_Squared = (TargetXYZs(:,1) - FOV_TX_data(FromProbe,1)) .^2 + (TargetXYZs(:,2) - FOV_TX_data(FromProbe,2)) .^2 + (TargetXYZs(:,3) - FOV_TX_data(FromProbe,3)) .^2;
                    [Min_Distance_Squared, Idx_Min_Distances_Squared] = min(Distances_Squared);
                    if Min_Distance_Squared == 0 % then exclude distance to self
                        if OnlyOneTarget
                            continue
                        end
                        Distances_Squared(Idx_Min_Distances_Squared,:) = [];
                        [Min_Distance_Squared, ~] = min(Distances_Squared);
                    end
                    M3DD(FromProbe,M3DD_column) = sqrt(Min_Distance_Squared);
                end
            end % P finished
            % delete(gcp('nocreate'))
            MeanToProbeValues(AdjustForSection + 1 : AdjustForSection + ProbesInSection) = mean(M3DD,"omitnan");
            SM3DD = M3DD ./ MeanToProbeValues(AdjustForSection + 1 : AdjustForSection + ProbesInSection);

            for FromProbe = 1 : NumNONControlProbes
                if ZeroProbeCounts(FromProbe)
                    continue
                end
                FROM_FOV_StandardisedM3DDFolderName = strcat(FOV_StandardisedM3DDFolderName,'\',NONControlProbes(FromProbe));
                SM3DD_file = strcat(FROM_FOV_StandardisedM3DDFolderName,'\SM3DD_to_Section_',num2str(Section),'.mat');
                from_to_specific_SM3DD = SM3DD(ProbeBinary(:,FromProbe),:);
                save(SM3DD_file,"from_to_specific_SM3DD","-v7.3","-nocompression")
                MSM3DD(FromProbe,AdjustForSection + 1 : AdjustForSection + ProbesInSection) = mean(from_to_specific_SM3DD,"omitnan");
            end
            save(MeanToProbeValues_file,"MeanToProbeValues","-v7.3","-nocompression")



            save(MSM3DD_file,"MSM3DD","-v7.3","-nocompression")
            FOV_Num_Sections_done = Section;
            save(FOV_Completion_file,"FOV_Num_Sections_done")

            SectionTOC = toc(SectionTIC);
            MatchSlide = TMAs_for_FOVs_to_process == UniqueTMAs(TMA);
            MatchFOV = FOVs_for_FOVs_to_process == FOVs_for_this_TMA(FOV);
            MatchBoth = MatchSlide & MatchFOV;
            FOVCompletionTimes_sections(MatchBoth,1) = LoadDataTimeTOC;
            FOVCompletionTimes_sections(MatchBoth,2) = ProbeBinaryTOC;
            FOVCompletionTimes_sections(MatchBoth,2 + Section) = SectionTOC;
            FOVCompletionTimes = sum(FOVCompletionTimes_sections,2);
            save(FOVCompletionTimes_sections_file,"FOVCompletionTimes_sections","FOVCompletionTimes");
            if Section == NumSections
                writematrix([TMAs_for_FOVs_to_process,FOVs_for_FOVs_to_process,FOVCompletionTimes],FOVCompletionTimes_file)
            end
            % Section finished
        end
        % FOV finished
    end
    % TMA finished
end

delete(gcp('nocreate'))

% calculate sample level mean SM3DDs for each Sample Group
Sample_level_means_FolderName = '.\SM3DD\Sample_level_means';
if ~isfolder(Sample_level_means_FolderName)
    mkdir(Sample_level_means_FolderName)
end
Sample_level_mean_SM3DD_file = strcat(Sample_level_means_FolderName,'\Sample_level_mean_SM3DD.mat');
if ~isfile(Sample_level_mean_SM3DD_file)
    Sample_level_mean_SM3DD = nan(Num_Samples_to_process,NumNONControlProbes,NumNONControlProbes,'single');
    count = 0;
    SampleDataExists = ones(Num_Samples_to_process,1,'logical');
    for Sample = 1 : Num_Samples_to_process
        if Num_FOVs_per_Sample(Sample) == 1
            count = count + 1;
            Digits_in_current_FOV = ceil(log10((FOVs_for_FOVs_to_process(count) + 1)));
            FOV_as_string = string(FOVs_for_FOVs_to_process(count));
            while Digits_in_current_FOV < MaxDigitsInFOV
                FOV_as_string = append('0',FOV_as_string);
                Digits_in_current_FOV = Digits_in_current_FOV + 1;
            end
            MSM3DD_file = strcat(FOV_level_means_FolderName,'\MSM3DD_',TMAs_for_FOVs_to_process(count),'_','F',FOV_as_string,'.mat');
            if ~isfile(MSM3DD_file)
                ErrorMessage = strcat("No data for F",num2str(FOVs_for_FOVs_to_process(Sample)),"_on TMA_",TMAs_for_FOVs_to_process(Sample));
                ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
                writematrix(ErrorMessage,ErrorMessage_file)
                SampleDataExists(Sample) = false;
                continue
            end
            load(MSM3DD_file,"MSM3DD")
            Sample_level_mean_SM3DD(Sample,:,:) = MSM3DD;
        else
            Sample_MSM3DDs = nan(Num_FOVs_per_Sample(Sample),NumNONControlProbes,NumNONControlProbes,'single');
            Sample_ProbeCounts = nan(Num_FOVs_per_Sample(Sample),NumNONControlProbes,'single');
            for Sample_FOV = 1 : Num_FOVs_per_Sample(Sample)
                count = count + 1;
                Digits_in_current_FOV = ceil(log10((FOVs_for_FOVs_to_process(count) + 1)));
                FOV_as_string = string(FOVs_for_FOVs_to_process(count));
                while Digits_in_current_FOV < MaxDigitsInFOV
                    FOV_as_string = append('0',FOV_as_string);
                    Digits_in_current_FOV = Digits_in_current_FOV + 1;
                end
                MSM3DD_file = strcat(FOV_level_means_FolderName,'\MSM3DD_',TMAs_for_FOVs_to_process(count),'_','F',FOV_as_string,'.mat');
                if ~isfile(MSM3DD_file)
                    ErrorMessage = strcat("No data for F",num2str(FOVs_for_FOVs_to_process(count)),"_on TMA_",TMAs_for_FOVs_to_process(count));
                    ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
                    writematrix(ErrorMessage,ErrorMessage_file)
                    continue
                end
                load(MSM3DD_file,"MSM3DD")
                Sample_MSM3DDs(Sample_FOV,:,:) = MSM3DD;
                ProbeCounts_file = strcat(FOV_level_ProbeCounts_FolderName,'\ProbeCounts_',TMAs_for_FOVs_to_process(count),'_','F',FOV_as_string,'.mat');
                load(ProbeCounts_file,"ProbeCounts")
                Sample_ProbeCounts(Sample_FOV,:) = ProbeCounts;
            end
            Sample_Total_ProbeCounts = sum(Sample_ProbeCounts);
            Sample_Proportions_ProbeCounts = Sample_ProbeCounts ./ Sample_Total_ProbeCounts;
            Sample_Weighted_MSM3DDs = zeros(Num_FOVs_per_Sample(Sample),NumNONControlProbes,NumNONControlProbes,'single');
            for Sample_FOV = 1 : Num_FOVs_per_Sample(Sample)
                Sample_Weighted_MSM3DDs(Sample_FOV,:,:) = Sample_MSM3DDs(Sample_FOV,:,:) .* Sample_Proportions_ProbeCounts(Sample_FOV,:);
            end
            Sample_level_mean_SM3DD(Sample,:,:) = sum(Sample_Weighted_MSM3DDs);
        end
    end
    save(Sample_level_mean_SM3DD_file,"SampleDataExists","Sample_level_mean_SM3DD","Sample_Group_1_name","Sample_Group_2_name","Sample_Group_1","Sample_Group_2","Num_Sample_Group_1","Num_Sample_Group_2","Num_FOVs_per_Sample","-v7.3","-nocompression")
elseif ~isfile('.\SM3DD\DirectionalLog10P.mat')
    load(Sample_level_mean_SM3DD_file,"SampleDataExists","Sample_level_mean_SM3DD","Sample_Group_1_name","Sample_Group_2_name","Sample_Group_1","Sample_Group_2","Num_Sample_Group_1","Num_Sample_Group_2","Num_FOVs_per_Sample")
end

% t-tests comparing sample level mean SM3DD between Sample Groups
% Ratios are Sample Group 1 SM3DDs / Sample Group 2 SM3DDs
% Ratios > 1 == Closer in Sample Group 2
if ~isfile('.\SM3DD\DirectionalLog10P.mat')
    if Paired_comparison_binary
        [~,P_values] = ttest(Sample_level_mean_SM3DD(1:Num_Sample_Group_1,:,:),Sample_level_mean_SM3DD(Num_Sample_Group_1 + 1:end,:,:));
    else
        [~,P_values] = ttest2(Sample_level_mean_SM3DD(1:Num_Sample_Group_1,:,:),Sample_level_mean_SM3DD(Num_Sample_Group_1 + 1:end,:,:));
    end
    P_values = squeeze(P_values);
    Ratios = squeeze(mean(Sample_level_mean_SM3DD(1:Num_Sample_Group_1,:,:),"omitnan") ./ mean(Sample_level_mean_SM3DD(Num_Sample_Group_1 + 1:end,:,:),"omitnan"));

    % DirectionalLog10P values are positive if "Closer in Sample Group 2"
    NegativeLog10P = -log10(P_values); % All positive values
    DirectionalLog10P = NegativeLog10P;
    DirectionalLog10P(Ratios < 1) = - DirectionalLog10P(Ratios < 1); % Make negative if "Closer in Sample Group 1"

    save('.\SM3DD\DirectionalLog10P.mat',"DirectionalLog10P","-v7.3","-nocompression")
    writematrix([["",transpose(NONControlProbes)];NONControlProbes,string(DirectionalLog10P)],'.\SM3DD\DirectionalLog10P.csv')

    [H,critP] = fdr_bky(P_values);

    FROM = strings(NumNONControlProbes * NumNONControlProbes,1);
    TO = strings(NumNONControlProbes * NumNONControlProbes,1);
    FROM_TO = strings(NumNONControlProbes * NumNONControlProbes,1);
    for P = 1 : NumNONControlProbes
        FROM((P-1)*NumNONControlProbes + 1 : (P-1)*NumNONControlProbes + NumNONControlProbes,1) = NONControlProbes;
        TO((P-1)*NumNONControlProbes + 1 : (P-1)*NumNONControlProbes + NumNONControlProbes,1) = NONControlProbes(P);
    end
    for row = 1 : NumNONControlProbes * NumNONControlProbes
        FROM_TO(row) = strcat(FROM(row),{' to '},TO(row));
    end

    Ratios_column = Ratios(:);
    P_column = P_values(:);
    H_column = H(:);
    NegativeLog10P_column = NegativeLog10P(:);
    DirectionalLog10P_column = DirectionalLog10P(:);

    Headers = {'FROM-probe to TO-probe','FROM-probe','TO-probe',strcat('Ratio_',Sample_Group_1_name,'/',Sample_Group_2_name),'P','Passes FDR control','Negative log10(P)','Directional log10(P)'};
    STATStable_file = strcat('.\SM3DD\STATs_table_',Sample_Group_1_name,'_vs_',Sample_Group_2_name,'.csv');
    writematrix([Headers;[FROM_TO,FROM,TO,Ratios_column,P_column,H_column,NegativeLog10P_column,DirectionalLog10P_column]],STATStable_file)
    critP_file = strcat('.\SM3DD\critP_',Sample_Group_1_name,'_vs_',Sample_Group_2_name,'.csv');
    writematrix(critP,critP_file)
end

% PCA
if any(PCA_options)
    TotalProbecounts_file = strcat(PCA_output_FolderName,'\TotalProbecounts.mat');
    if isfile(TotalProbecounts_file)
        load(TotalProbecounts_file,"Any_Absent","FOVs_With_Absent_Transcripts","TotalProbecounts","TotalProbecounts_ExcludeFOVs","CompositeProbecounts","CompositeProbecounts_ExcludeFOVs","No_transcripts")
    else % Determine Total Probecounts
        CompositeProbecounts = zeros(NumNONControlProbes,Num_FOVs_to_process);
        for F = 1 : Num_FOVs_to_process

            Digits_in_current_FOV = ceil(log10((FOVs_for_FOVs_to_process(F) + 1)));
            FOV_as_string = string(FOVs_for_FOVs_to_process(F));
            while Digits_in_current_FOV < MaxDigitsInFOV
                FOV_as_string = append('0',FOV_as_string);
                Digits_in_current_FOV = Digits_in_current_FOV + 1;
            end

            ProbeCounts_file = strcat(FOV_level_ProbeCounts_FolderName,'\ProbeCounts_',TMAs_for_FOVs_to_process(F),'_F',FOV_as_string,'.mat');
            if ~isfile(ProbeCounts_file)
                ErrorMessage = strcat("No ProbeCounts_file for F",FOV_as_string,"_on TMA_",TMAs_for_FOVs_to_process(F));
                ErrorMessage_file = strcat('.\ErrorMessage_',ErrorMessage,'.txt');
                writematrix(ErrorMessage,ErrorMessage_file)
                continue
            end
            load(ProbeCounts_file,"ProbeCounts")
            CompositeProbecounts(:,F) = ProbeCounts;
        end
        TotalProbecounts = sum(CompositeProbecounts,2);
        No_transcripts = TotalProbecounts == 0;

        FOVs_With_Absent_Transcripts = any(CompositeProbecounts(~No_transcripts,:) == 0);
        Any_Absent = any(FOVs_With_Absent_Transcripts);
        CompositeProbecounts_ExcludeFOVs = CompositeProbecounts;
        if Any_Absent
            CompositeProbecounts_ExcludeFOVs(:,FOVs_With_Absent_Transcripts) = 0;
            TotalProbecounts_ExcludeFOVs = sum(CompositeProbecounts_ExcludeFOVs,2);

        else
            TotalProbecounts_ExcludeFOVs = TotalProbecounts;
        end
        save(TotalProbecounts_file,"Any_Absent","FOVs_With_Absent_Transcripts","TotalProbecounts","TotalProbecounts_ExcludeFOVs","CompositeProbecounts","CompositeProbecounts_ExcludeFOVs","No_transcripts","-v7.3","-nocompression")
    end

    Exclude_Absent = PCA_options(2);
    % Include_Absent is TRUE is told to AND any are absent or Exclude_Absent is FALSE)
    Include_Absent = PCA_options(1) & (Any_Absent | ~Exclude_Absent);

    if Include_Absent
        Num_PCA_dimensions_per_FROM_probe_file = strcat(PCA_output_FolderName,'\Num_PCA_dimensions_per_FROM_probe.mat');
        if isfile(Num_PCA_dimensions_per_FROM_probe_file)
            load(Num_PCA_dimensions_per_FROM_probe_file,"Num_PCA_dimensions_per_FROM_probe","Finished_PCA_for_FROM_probe")
            StartProbe_IncludeFOVs = Finished_PCA_for_FROM_probe + 1;
        else
            StartProbe_IncludeFOVs = 1;
            Num_PCA_dimensions_per_FROM_probe = zeros(NumNONControlProbes,1);
        end
    end

    if Exclude_Absent
        Num_PCA_dimensions_per_FROM_probe_file_ExcludeFOVs = strcat(PCA_output_FolderName,'\Num_PCA_dimensions_per_FROM_probe_ExcludeFOVsWithAbsentTranscripts.mat');
        if isfile(Num_PCA_dimensions_per_FROM_probe_file_ExcludeFOVs)
            load(Num_PCA_dimensions_per_FROM_probe_file_ExcludeFOVs,"Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs","Finished_PCA_for_FROM_probe_ExcludeFOVs");
            StartProbe_ExcludeFOVs = Finished_PCA_for_FROM_probe_ExcludeFOVs + 1;
        else
            StartProbe_ExcludeFOVs = 1;
            Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs = zeros(NumNONControlProbes,1);
        end
    end

    StartProbe = min([StartProbe_IncludeFOVs,StartProbe_ExcludeFOVs]);

    if ~(StartProbe > NumNONControlProbes)
        for FromProbe = StartProbe : NumNONControlProbes
            if No_transcripts(FromProbe)
                continue
            end

            if Include_Absent && ~(FromProbe < StartProbe_IncludeFOVs)
                Num_FROM_probes_Group_1 = sum(CompositeProbecounts(FromProbe,1:sum(Num_FOVs_per_Sample(1:Num_Sample_Group_1))));
                Num_FROM_probes_Group_2 = TotalProbecounts(FromProbe) - Num_FROM_probes_Group_1;
                if Num_FROM_probes_Group_1 == 0 || Num_FROM_probes_Group_2 == 0
                    Diff_in_mean_PCA_scores.(NONControlProbes(FromProbe)) = nan;
                    continue
                end
                Composite_SM3DD_for_P = zeros(TotalProbecounts(FromProbe),NumNONControlProbes,'single');
            end


            if Exclude_Absent && ~(FromProbe < StartProbe_ExcludeFOVs)
                Num_FROM_probes_Group_1_ExcludeFOVs = sum(CompositeProbecounts_ExcludeFOVs(FromProbe,1:sum(Num_FOVs_per_Sample(1:Num_Sample_Group_1))));
                Num_FROM_probes_Group_2_ExcludeFOVs = TotalProbecounts_ExcludeFOVs(FromProbe) - Num_FROM_probes_Group_1_ExcludeFOVs;
                if Num_FROM_probes_Group_1_ExcludeFOVs == 0 || Num_FROM_probes_Group_2_ExcludeFOVs == 0
                    Diff_in_mean_PCA_scores_ExcludeFOVs.(NONControlProbes(FromProbe)) = nan;
                    continue
                end
                Composite_SM3DD_for_P_ExcludeFOVs = zeros(TotalProbecounts_ExcludeFOVs(FromProbe),NumNONControlProbes,'single');
            end

            count = 0;
            count_ExcludeFOVs = 0;
            for F = 1 : Num_FOVs_to_process
                if CompositeProbecounts(FromProbe,F) == 0 || (CompositeProbecounts_ExcludeFOVs(FromProbe,F) == 0 && ~Include_Absent)
                    continue
                end
                Digits_in_current_FOV = ceil(log10((FOVs_for_FOVs_to_process(F) + 1)));
                FOV_as_string = string(FOVs_for_FOVs_to_process(F));
                while Digits_in_current_FOV < MaxDigitsInFOV
                    FOV_as_string = append('0',FOV_as_string);
                    Digits_in_current_FOV = Digits_in_current_FOV + 1;
                end

                for Section = 1 : NumSections
                    FOV_StandardisedM3DDFolderName = strcat(StandardisedM3DDFolderName,'\',TMAs_for_FOVs_to_process(F),'_','F',FOV_as_string);
                    FROM_FOV_StandardisedM3DDFolderName = strcat(FOV_StandardisedM3DDFolderName,'\',NONControlProbes(FromProbe));
                    SM3DD_file = strcat(FROM_FOV_StandardisedM3DDFolderName,'\SM3DD_to_Section_',num2str(Section),'.mat');
                    load(SM3DD_file,"from_to_specific_SM3DD")
                    if Section == NumSections
                        ProbesInSection = NumNONControlProbes - (NumSections-1)*SectionSize;
                    else
                        ProbesInSection = SectionSize;
                    end
                    AdjustForSection = (Section - 1) * SectionSize;

                    if Include_Absent && ~(CompositeProbecounts(FromProbe,F) == 0)
                        Composite_SM3DD_for_P(1 + count_ExcludeFOVs:CompositeProbecounts(FromProbe,F) + count_ExcludeFOVs,1 + AdjustForSection:ProbesInSection + AdjustForSection) = from_to_specific_SM3DD;
                    end
                    if Exclude_Absent && ~(CompositeProbecounts_ExcludeFOVs(FromProbe,F) == 0)
                        Composite_SM3DD_for_P_ExcludeFOVs(1 + count_ExcludeFOVs:CompositeProbecounts_ExcludeFOVs(FromProbe,F) + count_ExcludeFOVs,1 + AdjustForSection:ProbesInSection + AdjustForSection) = from_to_specific_SM3DD;
                    end
                end
                count = count + CompositeProbecounts(FromProbe,F);
                count_ExcludeFOVs = count_ExcludeFOVs + CompositeProbecounts_ExcludeFOVs(FromProbe,F);
            end

            if Include_Absent
                Any_NaN = sum(isnan(Composite_SM3DD_for_P)) > 0;
                [~,Score,~,~,~,~] = pca(Composite_SM3DD_for_P(:,~Any_NaN));

                if ~isempty(Score)
                    Diff_in_mean_PCA_scores = mean(Score(Num_FROM_probes_Group_1 + 1:end,:),"omitnan") - mean(Score(1:Num_FROM_probes_Group_1,:),"omitnan");
                    [~,Num_PCA_dimensions_per_FROM_probe(FromProbe)] = size(Diff_in_mean_PCA_scores);
                else
                    Diff_in_mean_PCA_scores = nan;
                end
                FROM_probe_PCA_diffs_file = strcat(FROM_probe_PCA_diffs_FolderName,'\IncludeFOVs_With_Absent_Transcripts\',string(NONControlProbes(FromProbe)),'.mat');
                save(FROM_probe_PCA_diffs_file,"Diff_in_mean_PCA_scores","-v7.3","-nocompression")
                Finished_PCA_for_FROM_probe = FromProbe;
                save(Num_PCA_dimensions_per_FROM_probe_file,"Num_PCA_dimensions_per_FROM_probe","Finished_PCA_for_FROM_probe","-v7.3","-nocompression")
            end

            if Exclude_Absent
                Any_NaN = sum(isnan(Composite_SM3DD_for_P_ExcludeFOVs)) > 0;
                [~,Score,~,~,~,~] = pca(Composite_SM3DD_for_P_ExcludeFOVs(:,~Any_NaN));

                if ~isempty(Score)
                    Diff_in_mean_PCA_scores = mean(Score(Num_FROM_probes_Group_1_ExcludeFOVs + 1:end,:),"omitnan") - mean(Score(1:Num_FROM_probes_Group_1_ExcludeFOVs,:),"omitnan");
                    [~,Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs(FromProbe)] = size(Diff_in_mean_PCA_scores);
                else
                    Diff_in_mean_PCA_scores = nan;
                end
                FROM_probe_PCA_diffs_file_ExcludeFOVs = strcat(FROM_probe_PCA_diffs_FolderName,'\ExcludeFOVs_With_Absent_Transcripts\',string(NONControlProbes(FromProbe)),'.mat');
                save(FROM_probe_PCA_diffs_file_ExcludeFOVs,"Diff_in_mean_PCA_scores","-v7.3","-nocompression")
                Finished_PCA_for_FROM_probe_ExcludeFOVs = FromProbe;
                save(Num_PCA_dimensions_per_FROM_probe_file_ExcludeFOVs,"Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs","Finished_PCA_for_FROM_probe_ExcludeFOVs","-v7.3","-nocompression")
            end
        end
    end

    % Generate Coomposite PCA diffs
    if Include_Absent
        Num_PCA_dimensions_per_FROM_probe_file = strcat(PCA_output_FolderName,'\Num_PCA_dimensions_per_FROM_probe.mat');
        load(Num_PCA_dimensions_per_FROM_probe_file,"Num_PCA_dimensions_per_FROM_probe")
        Composite_PCA_diffs_csv_file = strcat(PCA_output_FolderName,'\Composite_PCA_diffs_',Sample_Group_1_name,'_vs_',Sample_Group_2_name,'_Include_Absent.csv');
        if isfile(Composite_PCA_diffs_csv_file)
            ErrorMessage = strcat(Composite_PCA_diffs_csv_file,' already exists?!');
            writematrix(ErrorMessage,'.\ErrorMessage.txt')
        else
            Max_PCA_dimensions = max(Num_PCA_dimensions_per_FROM_probe);
            Composite_PCA_diffs = nan(NumNONControlProbes,Max_PCA_dimensions);

            for FromProbe = 1 : NumNONControlProbes
                if ~(Num_PCA_dimensions_per_FROM_probe(FromProbe) == 0)
                    FROM_probe_PCA_diffs_file = strcat(FROM_probe_PCA_diffs_FolderName,'\IncludeFOVs_With_Absent_Transcripts\',string(NONControlProbes(FromProbe)),'.mat');
                    load(FROM_probe_PCA_diffs_file,"Diff_in_mean_PCA_scores")
                    Composite_PCA_diffs(FromProbe,1:Num_PCA_dimensions_per_FROM_probe(FromProbe)) = Diff_in_mean_PCA_scores;
                end
            end
            Composite_PCA_diffs_mat_file = strcat(PCA_output_FolderName,'\Composite_PCA_diffs_',Sample_Group_1_name,'_vs_',Sample_Group_2_name,'_Include_Absent.mat');
            save(Composite_PCA_diffs_mat_file,"Composite_PCA_diffs","-v7.3","-nocompression")
            writematrix([NONControlProbes,string(Composite_PCA_diffs)],Composite_PCA_diffs_csv_file)
        end
    end
    if Exclude_Absent
        Num_PCA_dimensions_per_FROM_probe_file_ExcludeFOVs = strcat(PCA_output_FolderName,'\Num_PCA_dimensions_per_FROM_probe_ExcludeFOVsWithAbsentTranscripts.mat');
        load(Num_PCA_dimensions_per_FROM_probe_file_ExcludeFOVs,"Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs")
        Composite_PCA_diffs_csv_file = strcat(PCA_output_FolderName,'\Composite_PCA_diffs_',Sample_Group_1_name,'_vs_',Sample_Group_2_name,'_Exclude_Absent.csv');
        if isfile(Composite_PCA_diffs_csv_file)
            ErrorMessage = strcat(Composite_PCA_diffs_csv_file,' already exists?!');
            writematrix(ErrorMessage,'.\ErrorMessage.txt')
        else
            Max_PCA_dimensions = max(Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs);
            Composite_PCA_diffs = nan(NumNONControlProbes,Max_PCA_dimensions);

            for FromProbe = 1 : NumNONControlProbes
                if ~(Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs(FromProbe) == 0)
                    FROM_probe_PCA_diffs_file = strcat(FROM_probe_PCA_diffs_FolderName,'\ExcludeFOVs_With_Absent_Transcripts\',string(NONControlProbes(FromProbe)),'.mat');
                    load(FROM_probe_PCA_diffs_file,"Diff_in_mean_PCA_scores")
                    Composite_PCA_diffs(FromProbe,1:Num_PCA_dimensions_per_FROM_probe_ExcludeFOVs(FromProbe)) = Diff_in_mean_PCA_scores;
                end
            end
            Composite_PCA_diffs_mat_file = strcat(PCA_output_FolderName,'\Composite_PCA_diffs_',Sample_Group_1_name,'_vs_',Sample_Group_2_name,'_Exclude_Absent.mat');
            save(Composite_PCA_diffs_mat_file,"Composite_PCA_diffs","-v7.3","-nocompression")
            writematrix([NONControlProbes,string(Composite_PCA_diffs)],Composite_PCA_diffs_csv_file)
        end
    end
end