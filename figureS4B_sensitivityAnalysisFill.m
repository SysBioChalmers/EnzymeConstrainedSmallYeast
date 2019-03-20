substrate1 = 'glcIN'; %glcIN 

gluIn = 1000;
O2In = 1000;
growth = 1000;

addpath('sourceCode')

numberOfTests = 100;
aParam = 0.3%[0.3 0.6 1];

tmpModel = model;
originalS = tmpModel.S;

growthReaction = findIndex(tmpModel.rxns, 'biomassOUT');
ethanolReaction = findIndex(tmpModel.rxns, 'ethOUT');
acReaction = findIndex(tmpModel.rxns, 'acOUT');
glyReaction = findIndex(tmpModel.rxns, 'glyOUT');
dataPositions = [growthReaction, ethanolReaction, acReaction];

%Turn of decoupling
tmpModel = setParam(tmpModel, 'ub', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
tmpModel = setParam(tmpModel, 'lb', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);

valueObject = makeValueObject(3.6856,  0.1966, 0.012, 0.0269, 0.518500, 0.023400, 0.807900, 1.134800, 40, 1);

tmpModel.b(end,1) = 0;
tmpModel.b(end,2) = 0.095;

tmpModel = setParam(tmpModel,'ub',{'glcIN'},0);
tmpModel = setParam(tmpModel,'ub',{substrate1, 'o2IN'},[gluIn, O2In]);
tmpModel = setParam(tmpModel,'obj',{'GROWTH'}, 1);
tmpModel = setParam(tmpModel,'ub',{'ethOUT', 'acOUT'}, [1000, 1000]);
tmpModel = setParam(tmpModel,'obj',{substrate1}, -1);
tmpModel = setParam(tmpModel,'lb',{'GROWTH'}, 0);
tmpModel = setParam(tmpModel,'ub',{'GROWTH'}, 1000);
tmpModel = setParam(tmpModel,'lb',{'ATPX'}, 0);  %0.5 mol/h maintainence
tmpModel.b(end,2) = 0.1;
    
messuredFluxes = {'o2IN', 'ethOUT'};

resultIndex = [];
for i = 1:length(messuredFluxes)
    resultIndex = [resultIndex findIndex(model.rxns, messuredFluxes{i})];
end

growthRates = linspace(0.2, 0.45, 25);
newS = originalS;
kcatRow = full(newS(end, :));

randValues = 2.^(rand(numberOfTests,size(newS,2)) * 2 * aParam - aParam);

results = zeros(numberOfTests, length(growthRates), length(resultIndex));  
for i = 1:numberOfTests
  pertubatedKcat = kcatRow .* randValues(i,:);
  newS(end, :) = pertubatedKcat;
  tmpModel.S = newS;
  results(i,:,:) = runChemostat(tmpModel, growthRates, valueObject, resultIndex);
end


clf
hold all
cutOff = max(growthRates(prctile(results(:,:,1), 10)>0));
medianValues = median(results);

color = {'r', 'g', 'k'};
for i = 1:length(messuredFluxes)
    plot(growthRates, medianValues(:,:,i), [color{i} '--'], 'linewidth', 3)    
end

legend(messuredFluxes, 'location', 'nw')

for i = 1:length(messuredFluxes)
    for j = 1:5
        resultUp = prctile(results(:,:,i), 80 - 5*j);
        resultDown = prctile(results(:,:,i), 20 + 5*j);
        X=[growthRates,fliplr(growthRates)];
        Y=[resultDown, fliplr(resultUp)];
        cVal = 0.95 - 0.08*j;
        fill(X,Y, [cVal cVal cVal], 'edgecolor', 'none')
    end
end



for i = 1:length(messuredFluxes)
    plot(growthRates, medianValues(:,:,i), [color{i} '--'], 'linewidth', 3)  
end

xlabel('Growthrate', 'FontSize',15,'FontName', 'Arial')
ylabel('flux mMol/h/g dw', 'FontSize',15,'FontName', 'Arial')
xlim([min(growthRates), cutOff]);

set(gca,'FontSize',15,'FontName', 'Arial')

hold off
