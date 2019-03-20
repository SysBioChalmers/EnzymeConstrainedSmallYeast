clc
addpath('sourceCode')
%Load matlab model generated from data/model.xlsx using the RAVEN toolbox.
load('data/model.mat');

%Load specific activity data.
model = mapDataToRxns(model, 'data/RxnAndSA.txt');

%Add a mass constraint metabolite to each enzymatic reaction in the S-matrix
model = addSpecificActivityConstraint(model, 0.5, 0.1, 60);

%Make the S matrix strictly positive.
model = addReversedReactions(model);

%Chose substrate
substrate1 = 'glcIN'; %glcIN 

GAM = 40;

gluIn = 1000;
O2In = 1000;

%Set protein constraint
model.b(end,1) = 0;
model.b(end,2) = 0.107;


%Constrain exchange fluxes
model = setParam(model,'ub',{'glcIN', 'o2IN'},[0, 0]);
model = setParam(model,'ub',{substrate1, 'o2IN'},[gluIn, O2In]);
%model = setParam(model,'ub',{'ethOUT'}, [1000]);
model = setParam(model,'lb',{'ATPX'}, [0.7]);  %0.7 mol/h maintainence

%Values taken from experimental data
model = setParam(model,'lb',{'acOUT'}, [1.3]);
model = setParam(model,'lb',{'glyOUT'}, [1.7]);
model = setParam(model,'ub',{'glyOUT'}, [1.7]);
model = setParam(model,'lb',{'ethOUT'}, [29.6]);
%model = setParam(model,'ub',{'ethOUT'}, [29.6]);


%Objective function
model = setParam(model,'obj',{'GROWTH'}, 1);
model = setParam(model,'ub',{'GROWTH'}, 1000);

%model = setParam(model,'lb',{'ShuttleX'}, [0]);  
%model = setParam(model,'ub',{'ShuttleX'}, [0]);   


%Allow or dissallow uncoupling, preventing uncoupling requires the blocking
%of some reactions
decouple = true;

if decouple 
    model = setParam(model,'lb',{'HDECOUP'}, [0]);  
    model = setParam(model,'ub',{'HDECOUP'}, 1000); 
else
    model = setParam(model,'ub',{'HDECOUP'}, 0);
    model = setParam(model, 'ub', {'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
    model = setParam(model, 'lb', {'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
end

%Reactions of interest
messuredFluxes = {'co2OUT', substrate1, 'ethOUT', 'acOUT', 'glyOUT'};

resultIndex = [];
for i = 1:length(messuredFluxes)
    resultIndex = [resultIndex findIndex(model.rxns, messuredFluxes{i})];
end

%Set up biomass equation
valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, GAM, 1);    

%Run optimization
resX1 = runOptimization(model, valueObject);    
results = resX1(resultIndex);
growthRate = resX1(findIndex(model.rxns, 'GROWTH'));


%Print results
batchAprox = [33.8 19.9 29.6 1.3 1.7];

yield = growthRate/(results(2)*180/1000);

fprintf('Parameters\tModel\tExperiment\tRelation\n', growthRate);
for i = 1:length(messuredFluxes)
   fprintf('%s\t%2.1f\t%2.1f\t%2.1f%%\n', messuredFluxes{i}, results(i), batchAprox(i), 100*results(i)/batchAprox(i));
end
   fprintf('GrowthRate\t%2.3f\t0.40\n', growthRate);
   fprintf('Yield\t%2.2f\t0.11\n', yield);
   
   
%%
massVector = full(model.S(end,:));
allRxns = {};
allFluxes = [];
allMass = [];
allSubs = [];
for i = 1:length(model.rxns)
   fluxAndEnzyme = and(resX1(i)>10^-6, massVector(i)>0);
   if fluxAndEnzyme
       currentMass = resX1(i)*massVector(i);
       currentRxn = split(model.rxns(i), '_back');
        allRxns = [allRxns;currentRxn(1)];
        allFluxes = [allFluxes; resX1(i)];
        allMass = [allMass; currentMass];
        allSubs = [allSubs; model.subSystems(i)];
   end
end

%Normalize all mass
allMass = allMass/sum(allMass);

%Rename IDPH to include IDH
allRxns{contains(allRxns,'IDPH')} = 'IDH+IDPH';

SAData = importdata('data/proteomicsData.txt');
proteomics = SAData.data;

mappedProteomics = zeros(length(allMass),1);

for i = 1:length(mappedProteomics)
    mappedIndx = contains(SAData.textdata, allRxns{i});
    mappedProteomics(i) = proteomics(mappedIndx);
end

mappedProteomics = mappedProteomics/sum(mappedProteomics);

logX = log10(mappedProteomics);
logY = log10(allMass);

hold all


fill([-5 0 0 -5], [-5-1, -1, 1, 1-5], 0.8*[1 1 1], 'facealpha',.5, 'edgecolor', 'none')
fill([-5 0 0 -5], [-5-log10(2), -log10(2), log10(2), + log10(2)-5], 0.8*[1 1 1], 'facealpha',.5, 'edgecolor', 'none')
plot([-5 0], [-5 0], '-', 'linewidth', 2, 'color', [0.6 0.6 0.6])

allSubs(ismember(allSubs,'Oxidative Phosphorylation')) = {'OXPHOS'}; 
subSystemCategories = {'Glycolysis', 'TCA', 'OXPHOS', 'Other'};
fillColors = [106 189 69
              237 34 36
              58 83 164
              100 100 100
            ]/256;
allSubs(not(ismember(allSubs,subSystemCategories))) = {'Other'}; 
        
        
for i = 1:length(allRxns)
    curSub = find(contains(subSystemCategories, allSubs{i}));
    scatter(logX(i), logY(i), 'MarkerFaceColor', fillColors(curSub,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.75)
    if abs(logX(i)-logY(i))>log10(10)
        text(logX(i)+0.1, logY(i), sprintf('%s',allRxns{i}))
    end
end

mdl = fitlm(logX, logY);
slope = mdl.Coefficients.Estimate;
standardError = mdl.Coefficients.SE;
xValues = linspace(min(logX),max(logY));
% plotErrorSlope(xValues,slope,standardError*2)
% plotErrorSlope(xValues,slope,standardError)

[r, p] = corr(logX, logY);
plot(xValues, xValues*slope(2) + slope(1), 'k', 'linewidth', 2)
text(-1, -1.5, sprintf('r=%2.1f\np=%2.1e', r, p), 'color', 'k')

%[r, p] = corr(logX(ismember(allSubs,'Glycolysis')), logY(ismember(allSubs,'Glycolysis')))
%[r, p] = corr(logX(ismember(allSubs,'TCA')), logY(ismember(allSubs,'TCA')))
%[r, p] = corr(logX(ismember(allSubs,'OXPHOS')), logY(ismember(allSubs,'OXPHOS')))
%[r, p] = corr(logX(ismember(allSubs,'Other')), logY(ismember(allSubs,'Other')))

axis equal
xlim([-4.5 -0])
ylim([-4.5 0])
xlabel('log10(proteomics)')
ylabel('log10(model)')

%%
figure()
subMassProteomics = zeros(length(subSystemCategories),1);
subMassModel = zeros(length(subSystemCategories),1);
for i = 1:length(subSystemCategories)
    curRxns = ismember(allSubs, subSystemCategories(i));
    subMassProteomics(i) = sum(mappedProteomics(curRxns));
    subMassModel(i) = sum(allMass(curRxns));    
end

subplot(1,2,1)
pie(subMassProteomics, subSystemCategories)
title('Proteomics')

subplot(1,2,2)
pie(subMassModel, subSystemCategories)
title('Model')
