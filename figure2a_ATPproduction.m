addpath('sourceCode')
%Load matlab model generated from data/model.xlsx using the RAVEN toolbox.
load('data/model.mat');

%Load specific activity data.
model = mapDataToRxns(model, 'data/RxnAndSA.txt');

%Add a mass constraint metabolite to each enzymatic reaction in the S-matrix
model = addSpecificActivityConstraint(model, 0.5, 0.1, 60);

%Make the S matrix strictly positive.
model = addReversedReactions(model);
substrate1 = 'glcIN'; %glcIN 

gluIn = 1000;
O2In = 1000;

expData = importdata('data/chemostatData.txt');
dataTable = expData.data;
%1 'D (h−1)' 2    'Yield (g · g−1)'  3  'qO2'   4 'qCO2'   5 'qglucose'    6 'qethanol'  7  'qacetate'   8 'qpyruvate'   9 'qglycerol'



model.b(end,1) = 0;

model = setParam(model,'ub',{'glcIN', 'o2IN'},[0, 0]);
model = setParam(model,'ub',{substrate1, 'o2IN'},[gluIn, O2In]);
model = setParam(model,'ub',{'ethOUT'}, [1000]);
model = setParam(model,'lb',{'ATPX'}, 0); 
model = setParam(model,'obj',{substrate1}, -1);
model = setParam(model,'ub',{'GROWTH'}, 1000);


messuredFluxes = {'o2IN', 'co2OUT', substrate1, 'ethOUT'};

resultIndex = [];
for i = 1:length(messuredFluxes)
    resultIndex = [resultIndex findIndex(model.rxns, messuredFluxes{i})];
end

growthRate = 5:0.1:35;

proteinContent = linearProteinContent(growthRate);

factor = 6.8;
estimatedATPDemand = factor/0.11 .* proteinContent';

results = zeros(length(growthRate), length(resultIndex));
for i=1:length(growthRate)
    proteinAmount = linearProteinContent(growthRate(i));

    model.b(end,2) = 0.1;

    %Original biomass equation
    valueObject = makeValueObjectWeight(0,  0, 0, 0,0, 0, 0, 0, 1, 1);
     
    model = setParam(model,'obj',{substrate1}, -1);
    model = setParam(model,'lb',{'GROWTH'}, growthRate(i));
    model = setParam(model,'ub',{'GROWTH'}, growthRate(i));
    resX1 = runOptimization(model, valueObject);    
    results(i,:) = resX1(resultIndex);
end

clf
hold all

growthRate(results(:,1) == 0) =[];
results(results(:,1) == 0,:) = [];

color = {'r', 'g', 'b', 'k', 'm'};

for i = 1:(length(messuredFluxes))
    plot(growthRate, results(:,i), [color{i} '-'], 'linewidth', 2)
end
    
    %legend(messuredFluxes)
    xlabel('ATP Production', 'FontSize',15,'FontName', 'Arial')
    ylabel('flux mMol/h/g dw', 'FontSize',15,'FontName', 'Arial')
    set(gca,'XTick', 0:5:40)
    xlim([min(growthRate), max(growthRate)*1.05])
    ylim([0 25])
    legend({'O2',  'CO2'  , 'glucose', 'ethanol'}, 'location', 'nw')
    
    set(gca,'FontSize',15,'FontName', 'Arial')
    
 