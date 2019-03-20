clc
substrate1 = 'glcIN';
%substrate1 = 'acIN';

%Select the condition:
%1=growth 2=growth W/O uncoupling 3=critical dilution rate
conditionNr = 3; 


gluIn = 1000;
O2In = 1000;
proteinLimit = 0.1;
fluxTreshold = 10^-5;
sampling = [1 1.001 inf];


addpath('sourceCode')

%load('data/model');

model = setParam(model,'ub',{'glcIN', 'o2IN'},[0, O2In]);
model = setParam(model,'ub',{substrate1, 'o2IN'},[gluIn, O2In]);


resultIndex = [
    findIndex(model.rxns, 'glcIN');
    findIndex(model.rxns, 'ethOUT');    
    findIndex(model.rxns, 'GROWTH');
    ];

% Objective: Biomass Equation
valueObject = makeValueObject(3.6856,  0.1966, 0.012, 0.0269, 0.518500, 0.023400, 0.807900, 1.134800, 35, 1);

% Objective: ATP maximization
%valueObject = makeValueObject(0,  0, 0, 0, 0, 0, 0, 0, 80, 1);






model.b(end,1) = 0;
model.b(end,2) = proteinLimit;

model = setParam(model,'ub',{substrate1, 'o2IN'},[gluIn, O2In]);
model = setParam(model,'obj',{'GROWTH'}, 1);
model = setParam(model,'lb',{'GROWTH'}, 0);
model = setParam(model,'ub',{'GROWTH'}, 1000);

%Critical = no ethanol or acetate production
%Mu max free (1000) ethanol or acetate production
model = setParam(model,'ub',{'ethOUT'}, 1000);
model = setParam(model,'ub',{'acOUT'}, 1000);
if conditionNr == 2
    %Turn of decoupling
    model = setParam(model, 'ub', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
    model = setParam(model, 'lb', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);

    model = setParam(model,'lb',{'HDECOUP'}, [0]);  
    model = setParam(model,'ub',{'HDECOUP'}, [0]);   
elseif conditionNr == 3
    %Turn of fermentation
    model = setParam(model,'ub',{'ethOUT'}, 0);
    model = setParam(model,'ub',{'acOUT'}, 0);
end

fluxes = runOptimization(model, valueObject);

nonZeroFluxes = find(and(fluxes > fluxTreshold, full(model.S(end,:))'>0));

%hold all
results = zeros(length(sampling), length(nonZeroFluxes));
resultDelta = zeros(length(sampling), length(nonZeroFluxes));
resultsEthanol =  zeros(length(sampling), length(nonZeroFluxes));

for i = 1:length(nonZeroFluxes)
    tmpModel = model;
    currentSA = 0.5 * 60 * model.specificActivity(nonZeroFluxes(i));
    currentSA = currentSA .* sampling;
    allWeights =  1./currentSA;

    for j = 1:length(sampling)
        tmpModel.S(end, nonZeroFluxes(i)) = allWeights(j);
        allflux = runOptimization(tmpModel, valueObject);
        results(j, i) = allflux(resultIndex(3));
        resultsEthanol(j,i) = allflux(resultIndex(2));
    end
    resultDelta(:,i) = currentSA;
end

resultsGrowthLog = log(results);
resultsFluxesLog = log(resultDelta);
fprintf('------------\n\n')
fprintf('c\tSpecific Activity\tFlux\tName\n')
fprintf('------------------\n')

dGH = (results(2, :)./results(1, :) - 1);
dXH = abs(resultsFluxesLog(2, :) - resultsFluxesLog(1, :));
dXO = abs(resultsFluxesLog(3, :) - resultsFluxesLog(1, :));
CH = dGH./(sampling(2)-1);

growthBost = results(3,:)./results(1,:);
ethanolBoost = resultsEthanol(3,:)./resultsEthanol(1,:);

for i = 1:length(model.rxns)
    currentNr = find(i == nonZeroFluxes);
    if currentNr ~= 0
        fluxValue = fluxes(nonZeroFluxes(currentNr));
        specificActivity = model.specificActivity(nonZeroFluxes(currentNr));
        fprintf('%2.4f\t%2.0f\t%2.5f\t%s\n', CH(currentNr), specificActivity, fluxValue, model.rxns{i})
        %fprintf('%s\t%2.3f\n',  model.rxns{i}, CH(currentNr))
    else
        fprintf('\t\t\t%s\n', model.rxns{i})        
    end
end