%load file with names and info of models
[~,infoFile,~] = xlsread('infofile.xlsx');

%get number of rows - number of samples
numRows = size(infoFile, 1) - 1;

modPath = [pwd filesep 'models20']; %location where models are stored - xml/sbml
modelFiles = dir(fullfile(modPath, '*.xml')); %list fo file names .xml 
inputModels = cell(numRows,1); %initialize a cell to contain all the models

% Create a cell array with headers
meta = cell(21,2);
meta{1,1} = 'index';
meta{1,2} = 'model_name';

%load models into a cell
for i=1:numRows
    meta{i+1,1} = i;
    meta{i+1,2} = modelFiles(i).name;
    modelPath = fullfile(modPath, modelFiles(i).name); %path to model location
    modelnames{i,1} = modelFiles(i).name;
    model = readCbModel(modelPath); %read model

    inputModels{i,1} = model; %store models in a cell - a column array sortof
    %nameTagsModels{i,1} = sprintf('M%d_', i);
end
disp("Done with loading the models!");


%%
writecell(meta,'metadata.csv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%diet - for monoculture
diet = readtable('WesternDietAGORA2.txt');
dietConstraints = table2cell(diet);
dietConstraints(:,2) = cellstr(num2str(cell2mat(dietConstraints(:,2))));
dietConstraints;

%community_diet - change all extracellular reactions to uptake (e) -> [u]
com_diet = readtable('com_diet_western.txt');
com_dietConstraints = table2cell(com_diet);
com_dietConstraints(:,2) = cellstr(num2str(cell2mat(com_dietConstraints(:,2))));
com_dietConstraints;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%turn off warnings
warning off;

%monoculture simulation
data = cell(21,2);
data{1,1} = "organism";
data{1,2} = "anaer_growth";
%data{1,3} = "aer_growth";

%monoculture growth rates
for i = 1:numRows
    disp(["Current model:", i]);
    data{i+1,1} = i;
    model = inputModels{i,1};

    %find and set biomass as objective
    biomass_mono = model.rxns(find(strncmp(model.rxns, 'biomass', 7))); %finds the reaction name
    model = changeObjective(model, biomass_mono); 
    
    %set diet
    model = useDiet(model, dietConstraints);

    %FBA without oxygen - anaerobic
    model = changeRxnBounds(model,'EX_o2(e)',0,'l');
    FBA = optimizeCbModel(model,'max');
    data{i+1,2} = FBA.f;

    % Enable uptake of oxygen - aerobic
    %model = changeRxnBounds(model,'EX_o2(e)',-10,'l');
    %FBA = optimizeCbModel(model,'max');
    %data{i+1,3} = FBA.f;
end


%%
writecell(data, 'monoculture_growth.csv');


%%
% Extract keys and values from the cell array
keysArray = data(2:end, 1); % First column as keys
valuesArray = data(2:end, 2); % Second column as values
% Create a map from the keys and values
mono_growth = containers.Map(keysArray, valuesArray);
%for i=1:20
%    mono_growth(i)
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 2 membered communities - making pairs from the given list of species using nchoosek

% pairs from 20species - 20C2 = 190 pairwise models
combinations = nchoosek(1:size(inputModels, 1), 2); %getting the different combinations
pairwiseResults = cell(size(combinations, 1), 1); %initialization
pairwiseNames = cell(size(combinations, 1), 1); %initialization
combinations;


%%
%pairwise culture growth

resultCell = cell(191, 5);
resultCell{1,1} = "idx1";
resultCell{1,2} = "idx2";
resultCell{1,3} = "Common_growth_anaer";
resultCell{1,4} = "M1growth_anaer";
resultCell{1,5} = "M2growth_anaer";
%resultCell{1,6} = "Common_growth_aer";
%resultCell{1,7} = "M1growth_aer";
%resultCell{1,8} = "M2growth_aer";

totalIterations = 190; % Total number of iterations

% Create a waitbar
h = waitbar(0, 'Calculating FBA...');

for i = 1:190
    idx1 = combinations(i,1); %get index of first microbe
    idx2 = combinations(i,2); %get index of second microbe
    resultCell{i+1,1} = idx1;
    resultCell{i+1,2} = idx2;

    models = cell(2,1); %initialization
    nameTagsModels = cell(2,1); %initialization

    model1 = inputModels{idx1, 1}; %access the model from the inputModels and store in model1
    model2 = inputModels{idx2, 1}; %access and store second model
    models{1,1} = model1; %creating a smaller cell(models) to store the pairwise models
    models{2,1} = model2; 
    nameTagsModels{1,1} = sprintf('M%d_', idx1);
    nameTagsModels{2,1} = sprintf('M%d_', idx2);

    [modelCom] = createMultipleSpeciesModel(models, 'nameTagsModels', nameTagsModels);
    [infoCom, indCom] = getMultiSpeciesModelId(modelCom, nameTagsModels);

    % Store the model info and indices
    modelCom.infoCom = infoCom(:);
    modelCom.indCom = indCom(:);
   
    %set biomass rxns of the models as objective functions
    biomass_reaction1 = sprintf('M%d_biomass', idx1);
    M1Biomass = find(contains(modelCom.rxns, biomass_reaction1, 'IgnoreCase', true));
    modelCom.c(M1Biomass,1) = 1; %index based assignment

    %modelCom = changeObjective(modelCom, modelCom.rxns(find(strncmp(modelCom.rxns,'M1_biomass',10))));
    biomass_reaction2 = sprintf('M%d_biomass', idx2);
    M2Biomass = find(contains(modelCom.rxns, biomass_reaction2, 'IgnoreCase', true));
    modelCom.c(M2Biomass,1) = 1; %index based assignment
    
    %modelCom = changeObjective(modelCom, [M1Biomass, M2Biomass]);
    %modelCom = changeObjective(modelCom, modelCom.rxns(find(strncmp(modelCom.rxns,'M2_biomass',10))));
    
    %use diet
    modelCom = useDiet(modelCom, com_dietConstraints);

    %FBA without oxygen - anaerobic
    modelCom = changeRxnBounds(modelCom,'EX_o2[u]',0,'l');  
    FBA = optimizeCbModel(modelCom,'max');
    comgrowth = FBA.f;
    
    %disp(["Comgrowth:",comgrowth)]
    M1growth = FBA.x(M1Biomass);
    M2growth = FBA.x(M2Biomass);
    resultCell{i+1,3} = comgrowth;
    resultCell{i+1,4} = M1growth;
    resultCell{i+1,5} = M2growth;

    % Enable uptake of oxygen - aerobic FBA
    %modelCom = changeRxnBounds(modelCom,'EX_o2[u]',-10,'l');
    %FBA = optimizeCbModel(modelCom,'max');
    %comgrowth = FBA.f;
    %M1growth = FBA.x(M1Biomass);
    %M2growth = FBA.x(M2Biomass);
    %resultCell{i+1,6} = comgrowth;
    %resultCell{i+1,7} = M1growth;
    %resultCell{i+1,8} = M2growth;

     % Update the waitbar
    waitbar(i / totalIterations, h, sprintf('Processing iteration %d of %d...', i, totalIterations));
end
close(h);
%%
writecell(resultCell, 'pairwise_growth.csv')
%%
resultCell = readcell('pairwise_growth.csv');


%%
%pairwise data annotations
% annotations - alpha values and interaction types
resultCell{1,6} = 'M1growth_mono';
resultCell{1,7} = 'M2growth_mono';
resultCell{1,8} = 'alpha1';
resultCell{1,9} = 'alpha2';
resultCell{1,10} = 'interaction_type';
resultCell{1,11} = 'interaction_subclass1';
resultCell{1,12} = 'interaction_subclass2';

for i=1:190
    ind1 = resultCell{i+1,1};
    ind2 = resultCell{i+1,2};
    com_grow1 = resultCell{i+1,4};
    com_grow2 = resultCell{i+1,5};
    wild_grow1 = mono_growth(ind1);
    wild_grow2 = mono_growth(ind2);
 
    alpha1 = (com_grow1 - wild_grow1)*100/wild_grow1;
    alpha2 = (com_grow2 - wild_grow2)*100/wild_grow2;

    resultCell{i+1,6} = wild_grow1;
    resultCell{i+1,7} = wild_grow2;
    resultCell{i+1,8} = alpha1;
    resultCell{i+1,9} = alpha2;

    interactionType = '';
    int_subtype1 = '';
    int_subtype2 = '';

    %update the following to a function - refer to the end section
    % Determine interaction type based on the ranges of alpha1 and alpha2
    if (-10 <= alpha1 && alpha1 <= 10) && (alpha2 <= -10)
        interactionType = 'amensalism';
        int_subtype1 = 'amensal_unnaffected';
        int_subtype2 = 'amensal_affected';
    elseif (alpha1 <= -10) && (-10 <= alpha2 && alpha2 <= 10)
        interactionType = 'amensalism';
        int_subtype1 = 'amensal_affected';
        int_subtype2 = 'amensal_unaffected';
    elseif (-10 <= alpha1 && alpha1 <= 10) && (alpha2 >= 10)
        interactionType = 'commensalism';
        int_subtype1 = 'commensal_unaffected';
        int_subtype2 = 'commensal_taker';
    elseif (alpha1 >= 10) && (-10 <= alpha2 && alpha2 <= 10)
        interactionType = 'commensalism';
        int_subtype1 = 'commensal_taker';
        int_subtype2 = 'commensal_unaffected';
    elseif (alpha1 <= -10) && (alpha2 <= -10)
        interactionType = 'competition';
    elseif (alpha1 >= 10) && (alpha2 >= 10)
        interactionType = 'mutualism';
    elseif (-10 <= alpha1 && alpha1 <= 10) && (-10 <= alpha2 && alpha2 <= 10)
        interactionType = 'neutralism';
    elseif (alpha1 <= -10) && (alpha2 >= 10)
        interactionType = 'parasitism';
        int_subtype1 = 'parasitism_giver';
        int_subtype2 = 'parasitism_taker';
    elseif (alpha1 >= 10) && (alpha2 <= -10)
        interactionType = 'parasitism';        
        int_subtype1 = 'parasitism_taker';
        int_subtype2 = 'parasitism_giver';
    else
        interactionType = 'unknown_interaction';
    end

    resultCell{i+1,10} = interactionType;
    resultCell{i+1,11} = int_subtype1;
    resultCell{i+1,12} = int_subtype2;

end
%%
writecell(resultCell,'pairwise.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 3 membered communities - making trios from the given list of species using nchoosek

% trios from 20species - 20C3 = 1140 pairwise models
combinations_trio = nchoosek(1:size(inputModels, 1), 3); %getting the different combinations
trioResults = cell(size(combinations_trio, 1), 1); %initialization
trioNames = cell(size(combinations_trio, 1), 1); %initialization
combinations_trio;
%%
%trio culture growth
resultCell_trio = cell(1141, 7);
resultCell_trio{1,1} = "idx1";
resultCell_trio{1,2} = "idx2";
resultCell_trio{1,3} = "idx3";

resultCell_trio{1,4} = "Common_growth_anaer";
resultCell_trio{1,5} = "M1growth_anaer";
resultCell_trio{1,6} = "M2growth_anaer";
resultCell_trio{1,7} = "M3growth_anaer";
%resultCell_trio{1,8} = "Common_growth_aer";
%resultCell_trio{1,9} = "M1growth_aer";
%resultCell_trio{1,10} = "M2growth_aer";
%resultCell_trio{1,11} = "M3growth_aer";

totalIterations = 1140; % Total number of iterations

% Create a waitbar
h = waitbar(0, 'Calculating FBA...');

for i = 1:1140
    idx1 = combinations_trio(i,1); %get index of first microbe
    idx2 = combinations_trio(i,2); %get index of second microbe
    idx3 = combinations_trio(i,3);
    resultCell_trio{i+1,1} = idx1;
    resultCell_trio{i+1,2} = idx2;
    resultCell_trio{i+1,3} = idx3;
    
    models = cell(3,1);
    nameTagsModels = cell(3,1);
    model1 = inputModels{idx1, 1}; %access the model from the inputModels and store in model1
    model2 = inputModels{idx2, 1}; %access and store second model
    model3 = inputModels{idx3, 1};
    models{1,1} = model1; %creating a smaller cell(models) to store the pairwise models
    models{2,1} = model2; 
    models{3,1} = model3;
    nameTagsModels{1,1} = sprintf('M%d_', idx1);
    nameTagsModels{2,1} = sprintf('M%d_', idx2);
    nameTagsModels{3,1} = sprintf('M%d_', idx3);

    [modelCom] = createMultipleSpeciesModel(models, 'nameTagsModels', nameTagsModels);
    [infoCom, indCom] = getMultiSpeciesModelId(modelCom, nameTagsModels);

    % Store the model info and indices
    modelCom.infoCom = infoCom(:);
    modelCom.indCom = indCom(:);
   
    %set biomass rxns of the models as objective functions
    biomass_reaction1 = sprintf('M%d_biomass', idx1);
    M1Biomass = find(contains(modelCom.rxns, biomass_reaction1, 'IgnoreCase', true));
    modelCom.c(M1Biomass,1) = 1; %index based assignment
    %modelCom = changeObjective(modelCom, modelCom.rxns(find(strncmp(modelCom.rxns,'M1_biomass',10))));


    biomass_reaction2 = sprintf('M%d_biomass', idx2);
    M2Biomass = find(contains(modelCom.rxns, biomass_reaction2, 'IgnoreCase', true));
    modelCom.c(M2Biomass,1) = 1; %index based assignment
    %modelCom = changeObjective(modelCom, modelCom.rxns(find(strncmp(modelCom.rxns,'M2_biomass',10))));


    biomass_reaction3 = sprintf('M%d_biomass', idx3);
    M3Biomass = find(contains(modelCom.rxns, biomass_reaction3, 'IgnoreCase', true));
    modelCom.c(M3Biomass,1) = 1; %index based assignment
    %modelCom = changeObjective(modelCom, modelCom.rxns(find(strncmp(modelCom.rxns,'M3_biomass',10))));

    %modelCom = changeObjective(modelCom, [M1Biomass, M2Biomass, M3Biomass]);
    
    %use diet
    modelCom = useDiet(modelCom, com_dietConstraints);

    %FBA without oxygen - anaerobic
    modelCom = changeRxnBounds(modelCom,'EX_o2[u]',0,'l');  
    FBA = optimizeCbModel(modelCom,'max');
    comgrowth = FBA.f;
    
    %disp(["Comgrowth:",comgrowth)]
    M1growth = FBA.x(M1Biomass);
    M2growth = FBA.x(M2Biomass);
    M3growth = FBA.x(M3Biomass);
    resultCell_trio{i+1,4} = comgrowth;
    resultCell_trio{i+1,5} = M1growth;
    resultCell_trio{i+1,6} = M2growth;
    resultCell_trio{i+1,7} = M3growth;

    % Enable uptake of oxygen - aerobic FBA
    %modelCom = changeRxnBounds(modelCom,'EX_o2[u]',-10,'l');
    %FBA = optimizeCbModel(modelCom,'max');
    %comgrowth = FBA.f;
    %M1growth = FBA.x(M1Biomass);
    %M2growth = FBA.x(M2Biomass);
    %resultCell{i+1,6} = comgrowth;
    %resultCell{i+1,7} = M1growth;
    %resultCell{i+1,8} = M2growth;

    % Update the waitbar
    waitbar(i / totalIterations, h, sprintf('Processing iteration %d of %d...', i, totalIterations));

end
close(h);
%%
writecell(resultCell_trio, 'trio_growth.csv');
%%
resultCell_trio = readcell("trio_growth.csv");
%%
% trio annotations - interaction types
resultCell_trio{1,8} = 'M1growth_mono';
resultCell_trio{1,9} = 'M2growth_mono';
resultCell_trio{1,10} = 'M3growth_mono';
resultCell_trio{1,11} = 'alpha1';
resultCell_trio{1,12} = 'alpha2';
resultCell_trio{1,13} = 'alpha3';
resultCell_trio{1,14} = 'interaction_type12';
resultCell_trio{1,15} = 'interaction_type13';
resultCell_trio{1,16} = 'interaction_type23';

resultCell_trio{1,17} = 'int_subclass12_1';
resultCell_trio{1,18} = 'int_subclass12_2';
resultCell_trio{1,19} = 'int_subclass13_1';
resultCell_trio{1,20} = 'int_subclass13_3';
resultCell_trio{1,21} = 'int_subclass23_2';
resultCell_trio{1,22} = 'int_subclass23_3';

for i=1:1140

    ind1 = resultCell_trio{i+1,1};
    ind2 = resultCell_trio{i+1,2};
    ind3 = resultCell_trio{i+1,3};
    com_grow1 = resultCell_trio{i+1,5};
    com_grow2 = resultCell_trio{i+1,6};
    com_grow3 = resultCell_trio{i+1,7};
    wild_grow1 = mono_growth(ind1);
    wild_grow2 = mono_growth(ind2);
    wild_grow3 = mono_growth(ind3);
 
    alpha1 = (com_grow1 - wild_grow1)*100/wild_grow1;
    alpha2 = (com_grow2 - wild_grow2)*100/wild_grow2;
    alpha3 = (com_grow3 - wild_grow3)*100/wild_grow3;

    resultCell_trio{i+1,8} = wild_grow1;
    resultCell_trio{i+1,9} = wild_grow2;
    resultCell_trio{i+1,10} = wild_grow3;
    resultCell_trio{i+1,11} = alpha1;
    resultCell_trio{i+1,12} = alpha2;
    resultCell_trio{i+1,13} = alpha3;

    [interactionType, int_subtype1, int_subtype2] = determineInteractionType(alpha1, alpha2);
    resultCell_trio{i+1,14} = interactionType; %'interaction_type12';
    resultCell_trio{i+1,17} = int_subtype1; %'int_subclass12_1';
    resultCell_trio{i+1,18} = int_subtype2; %'int_subclass12_2';

    [interactionType, int_subtype1, int_subtype2] = determineInteractionType(alpha1, alpha3);
    resultCell_trio{i+1,15} = interactionType; %'interaction_type13';
    resultCell_trio{i+1,19} = int_subtype1; %'int_subclass13_1';
    resultCell_trio{i+1,20} = int_subtype2; %'int_subclass13_3';

    [interactionType, int_subtype1, int_subtype2] = determineInteractionType(alpha2, alpha3);
    resultCell_trio{i+1,16} = interactionType; %'interaction_type23';
    resultCell_trio{i+1,21} = int_subtype1; %'int_subclass23_2';
    resultCell_trio{i+1,22} = int_subtype2; %'int_subclass23_3';

end

%%
writecell(resultCell_trio, 'trio.csv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
function [interactionType, int_subtype1, int_subtype2] = determineInteractionType(alpha1, alpha2)
    % Determine interaction type based on the ranges of alpha1 and alpha2
    if (-10 <= alpha1 && alpha1 <= 10) && (alpha2 <= -10)
        interactionType = 'amensalism';
        int_subtype1 = 'amensal_unaffected';
        int_subtype2 = 'amensal_affected';
    elseif (alpha1 <= -10) && (-10 <= alpha2 && alpha2 <= 10)
        interactionType = 'amensalism';
        int_subtype1 = 'amensal_affected';
        int_subtype2 = 'amensal_unaffected';
    elseif (-10 <= alpha1 && alpha1 <= 10) && (alpha2 >= 10)
        interactionType = 'commensalism';
        int_subtype1 = 'commensal_unaffected';
        int_subtype2 = 'commensal_taker';
    elseif (alpha1 >= 10) && (-10 <= alpha2 && alpha2 <= 10)
        interactionType = 'commensalism';
        int_subtype1 = 'commensal_taker';
        int_subtype2 = 'commensal_unaffected';
    elseif (alpha1 <= -10) && (alpha2 <= -10)
        interactionType = 'competition';
        int_subtype1 = '';
        int_subtype2 = '';
    elseif (alpha1 >= 10) && (alpha2 >= 10)
        interactionType = 'mutualism';
        int_subtype1 = '';
        int_subtype2 = '';
    elseif (-10 <= alpha1 && alpha1 <= 10) && (-10 <= alpha2 && alpha2 <= 10)
        interactionType = 'neutralism';
        int_subtype1 = '';
        int_subtype2 = '';
    elseif (alpha1 <= -10) && (alpha2 >= 10)
        interactionType = 'parasitism';
        int_subtype1 = 'parasitism_giver';
        int_subtype2 = 'parasitism_taker';
    elseif (alpha1 >= 10) && (alpha2 <= -10)
        interactionType = 'parasitism';
        int_subtype1 = 'parasitism_taker';
        int_subtype2 = 'parasitism_giver';
    else
        interactionType = 'unknown_interaction';
        int_subtype1 = '';
        int_subtype2 = '';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%