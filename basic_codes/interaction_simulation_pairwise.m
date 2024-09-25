%load file with names and info of models
[~,infoFile,~] = xlsread('infofile.xlsx');

%get number of rows - number of samples
numRows = size(infoFile, 1) - 1;

modPath = [pwd filesep 'redo_pairwise_gayathri']; %location where models are stored - xml/sbml
modelFiles = dir(fullfile(modPath, '*.xml')); %list fo file names .xml 
inputModels = cell(numRows,1); %initialize a cell to contain all the models

%load models into a cell
for i=1:numRows
    modelPath = fullfile(modPath, modelFiles(i).name); %path to model location
    modelnames{i,1} = modelFiles(i).name;
    model = readCbModel(modelPath); %read model
    inputModels{i,1} = model; %store models in a cell - a column array sortof
    %nameTagsModels{i,1} = sprintf('M%d_', i);
end
disp("Done with loading the models!");

%%
%diet - for monoculture
diet = readtable('WesternDietAGORA2.txt');
dietConstraints = table2cell(diet);
dietConstraints(:,2) = cellstr(num2str(cell2mat(dietConstraints(:,2))));
dietConstraints;

%communit_diet - change all extracellular reactions to uptake (e) -> [u]
com_diet = readtable('com_diet_western.txt');
com_dietConstraints = table2cell(com_diet);
com_dietConstraints(:,2) = cellstr(num2str(cell2mat(com_dietConstraints(:,2))));
com_dietConstraints;
%%
%monoculture simulation
data = cell(6,3);
data{1,1} = "organism";
data{1,2} = "anaer_growth";
data{1,3} = "aer_growth";

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
    model = changeRxnBounds(model,'EX_o2(e)',-10,'l');
    FBA = optimizeCbModel(model,'max');
    data{i+1,3} = FBA.f;
end

%display the data table - monoculture
data;

%%

%making pairs from the given list of species using nchoosek

% pairs from 5species - 5C2 = 10 pairwise models
combinations = nchoosek(1:size(inputModels, 1), 2); %getting the different combinations - using numbers 1to5
pairwiseResults = cell(size(combinations, 1), 1); %just initialization
pairwiseNames = cell(size(combinations, 1), 1); %just initialization
combinations;
%%

%pairwise culture growth
resultCell = cell(11, 8);
resultCell{1,1} = "idx1";
resultCell{1,2} = "idx2";
resultCell{1,3} = "Common_growth_anaer";
resultCell{1,4} = "M1growth_anaer";
resultCell{1,5} = "M2growth_anaer";
resultCell{1,6} = "Common_growth_aer";
resultCell{1,7} = "M1growth_aer";
resultCell{1,8} = "M2growth_aer";


for i = 1:10
    idx1 = combinations(i,1); %get index of first microbe
    idx2 = combinations(i,2); %get index of second microbe
    resultCell{i+1,1} = idx1;
    resultCell{i+1,2} = idx2;

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
    rxns = modelCom.rxns();
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
    modelCom = changeRxnBounds(modelCom,'EX_o2[u]',-10,'l');
    FBA = optimizeCbModel(modelCom,'max');
    comgrowth = FBA.f;
    M1growth = FBA.x(M1Biomass);
    M2growth = FBA.x(M2Biomass);
    resultCell{i+1,6} = comgrowth;
    resultCell{i+1,7} = M1growth;
    resultCell{i+1,8} = M2growth;

end
%%
resultCell
%%
%another way

%{
resultCell2 = cell(11, 8);
resultCell2{1,1} = "idx1";
resultCell2{1,2} = "idx2";
resultCell2{1,3} = "Common_growth_anaer";
resultCell2{1,4} = "M1growth_anaer";
resultCell2{1,5} = "M2growth_anaer";
resultCell2{1,6} = "Common_growth_aer";
resultCell2{1,7} = "M1growth_aer";
resultCell2{1,8} = "M2growth_aer";


for i = 1:10
    idx1 = combinations(i,1); %get index of first microbe
    idx2 = combinations(i,2); %get index of second microbe
    resultCell2{i+1,1} = idx1;
    resultCell2{i+1,2} = idx2;

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
    M1Biomass_idx = find(contains(modelCom.rxns, biomass_reaction1, 'IgnoreCase', true));
    M1Biomass = modelCom.rxns(find(strncmp(modelCom.rxns,biomass_reaction1, 10)));
    %modelCom.c(M1Biomass,1) = 1; %index based assignment
    %modelCom = changeObjective(modelCom, M1Biomass);
    
    biomass_reaction2 = sprintf('M%d_biomass', idx2);
    M2Biomass_idx = find(contains(modelCom.rxns, biomass_reaction2, 'IgnoreCase', true));
    M2Biomass = modelCom.rxns(find(strncmp(modelCom.rxns,biomass_reaction2, 10)));

    %modelCom.c(M2Biomass,1) = 1; %index based assignment
    
    modelCom = changeObjective(modelCom, [M1Biomass, M2Biomass]);

    %use diet
    modelCom = useDiet(modelCom, dietConstraints);

    %FBA without oxygen - anaerobic
    FBA = optimizeCbModel(modelCom,'max');
    comgrowth = FBA.f;
    %disp(["Comgrowth:",comgrowth)]
    M1growth = FBA.x(M1Biomass_idx);
    M2growth = FBA.x(M2Biomass_idx);
    resultCell2{i+1,3} = comgrowth;
    resultCell2{i+1,4} = M1growth;
    resultCell2{i+1,5} = M2growth;

    % Enable uptake of oxygen - aerobic FBA
    modelCom = changeRxnBounds(modelCom,'EX_o2(e)',-10,'l');
    FBA = optimizeCbModel(modelCom,'max');
    comgrowth = FBA.f;
    M1growth = FBA.x(M1Biomass_idx);
    M2growth = FBA.x(M2Biomass_idx);
    resultCell2{i+1,6} = comgrowth;
    resultCell2{i+1,7} = M1growth;
    resultCell2{i+1,8} = M2growth;
    
    
   

    %{
    %set medium uptake constraints - default example
    modelCom.lb(indCom.EXcom) = -1000;
    %FBA 
    sol = optimizeCbModel(modelCom,'max');
    %observe growth rate for community and individual model
    Commgrowth = sol.f;
    M1growth = sol.x(M1Biomass);
    M2growth = sol.x(M2Biomass);
    Commgrowth;
    M1growth;
    M2growth;
    %}

end
%}
