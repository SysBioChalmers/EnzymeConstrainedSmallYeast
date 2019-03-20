addpath('sourceCode')
load('data/model.mat');

%Load specific activity data.
model = mapDataToRxns(model, 'data/RxnAndSA.txt');

%Add a mass constraint metabolite to each enzymatic reaction in the S-matrix
model = addSpecificActivityConstraint(model, 0.5, 0.1, 60);

%Make the S matrix strictly positive.
model = addReversedReactions(model);

%Chose substrate
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
model = setParam(model,'lb',{'ATPX'}, [0.7]);  %mol/h maintainence
%Original constraint on acetate
%model = setCostParam(model,'ub',{'acOUT'}, [1000]);
%model = setParam(model,'ub',{'acOUT'}, [0]); 
model = setParam(model,'ub',{'acOUT'}, [0.6]); 


model = setParam(model,'obj',{substrate1}, -1);
model = setParam(model,'ub',{'GROWTH'}, 1000);

%model = setParam(model,'lb',{'ShuttleX'}, [0]);  
%model = setParam(model,'ub',{'ShuttleX'}, [0]);   

%Allow uncoupling
%model = setParam(model,'lb',{'HDECOUP'}, [0]);  
%model = setParam(model,'ub',{'HDECOUP'}, [0]);   

%Prevent uncoupling
model = setParam(model, 'ub', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
model = setParam(model, 'lb', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);



messuredFluxes = {'o2IN', 'co2OUT', substrate1, 'ethOUT', 'acOUT'};

resultIndex = [];
for i = 1:length(messuredFluxes)
    resultIndex = [resultIndex findIndex(model.rxns, messuredFluxes{i})];
end

growthRate = 0.02:0.01:0.41;

proteinContent = linearProteinContent(growthRate);

factor = 6.8;
estimatedATPDemand = factor/0.11 .* proteinContent';

results = zeros(length(growthRate), length(resultIndex));
for i=1:length(growthRate)
    proteinAmount = linearProteinContent(growthRate(i));

    %Protein limitation
    model.b(end,2) = 0.1;

    %Original biomass equation
    valueObject = makeValueObject(3.6856,  0.1966, 0.012, 0.0269, 0.518500, 0.023400, 0.807900, 1.134800, 35, 1);
   
    %Final biomass equation
    carbohydrateContent = 1 - proteinAmount -0.12 -0.006 -0.025 -0.005 -0.01 - 0.05; %0.05 for ash

    valueObject = makeValueObjectWeight(proteinAmount, 0.12, 0.006, 0.025, 0.005, 0.01, 0, carbohydrateContent, estimatedATPDemand(i), 1);    
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
    xlabel('Growthrate', 'FontSize',15,'FontName', 'Arial')
    ylabel('flux mMol/h/g dw', 'FontSize',15,'FontName', 'Arial')
    set(gca,'XTick', 0:0.05:0.45)
    xlim([min(growthRate), max(growthRate)*1.05])
    ylim([0 25])
    legend({'O2',  'CO2'  , 'glucose', 'ethanol', 'acetate'}, 'location', 'nw')
    
    set(gca,'FontSize',15,'FontName', 'Arial')
    
    
 for i = 1:(length(messuredFluxes))
    plot(dataTable(:,1), dataTable(:,i+2), [color{i} 'o'], 'linewidth', 2)
end