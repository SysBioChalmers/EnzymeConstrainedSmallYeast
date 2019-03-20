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
growth = 1000;

addpath('sourceCode')

numberOfTests = 1000;
aParam = [0.3 0.6 1];

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

valueObject = makeValueObjectWeight(0.4, 0.10, 0.006, 0.025, 0.005, 0.01, 0.1, 0.40, 40, 1);

tmpModel.b(end,1) = 0;
tmpModel.b(end,2) = 0.1;

tmpModel = setParam(tmpModel,'ub',{'glcIN'},0);
tmpModel = setParam(tmpModel,'ub',{substrate1, 'o2IN'},[gluIn, O2In]);
tmpModel = setParam(tmpModel,'obj',{'GROWTH'}, 1);
tmpModel = setParam(tmpModel,'ub',{'ethOUT', 'acOUT'}, [1000, 1000]);
tmpModel = setParam(tmpModel,'obj',{'GROWTH'}, 1);
tmpModel = setParam(tmpModel,'lb',{'GROWTH'}, 0);
tmpModel = setParam(tmpModel,'ub',{'GROWTH'}, 1000);
%can cause infeasibility
tmpModel = setParam(tmpModel,'lb',{'ATPX'}, 0);  %0.5 mol/h maintainence

result = zeros(length(aParam), numberOfTests, length(dataPositions));


newS = originalS;
kcatRow = full(newS(end, :));
for j = 1:length(aParam)
    randValues = 2.^(rand(numberOfTests,size(newS,2)) * 2 * aParam(j) - aParam(j));

    
   for i = 1:numberOfTests
      
      pertubatedKcat = kcatRow .* randValues(i,:);
      newS(end, :) = pertubatedKcat;
      tmpModel.S = newS;
      resX1 = runOptimization(tmpModel, valueObject);

      result(j, i, :) = resX1(dataPositions);
   end
end

close all
hold all
histogramPostions = linspace(0.25, 0.65, 50);
ethCount = zeros(length(aParam), 1);
acCount = zeros(length(aParam), 1);

for i = 1:length(aParam)
   [numbers positions] = hist(result(i, :,1), histogramPostions);
   plot(positions, numbers, 'o-', 'LineWidth',2);
   ethCount(i) = 100 * sum(result(i,:,2)>0.1)/numberOfTests;
   acCount(i) = 100 * sum(result(i,:,3)>0.1)/numberOfTests;
end

xlabel('Growth rate')
ylabel('Count')

val1 = ['0.2, Eth ' num2str(ethCount(1)) '%'];
val2 = ['0.5, Eth ' num2str(ethCount(2)) '%'];
val3 = ['1, Eth ' num2str(ethCount(3)) '%'];

legend({val1, val2, val3}, 'Location', 'NorthEast')
%%

close all
histogramPostions = linspace(0.2, 0.7, 50);
ethCount = zeros(length(aParam), 1);
acCount = zeros(length(aParam), 1);

maxValue = 0;

for i = 1:length(aParam)
    subplot(length(aParam),1,i);
    [numbers positions] = hist(result(i, :,1), histogramPostions);
    
    bar(positions, numbers, 'FaceColor', [0.36 0.61 0.82]);
    ethCount(i) = 100 * sum(result(i,:,2)>0.1)/numberOfTests;
    acCount(i) = 100 * sum(result(i,:,3)>0.1)/numberOfTests;
    legendString = sprintf('I: %1.1f, Eth %1.2f%%', aParam(i), ethCount(i));
    legend(legendString)
    xlabel('Growth rate', 'FontSize',15,'FontName', 'Arial')
    ylabel('Count', 'FontSize',15,'FontName', 'Arial')
    
    maxValue = max(maxValue, max(numbers));
    ylim([0 maxValue]);
    set(gca,'FontSize',15,'FontName', 'Arial')
   % xlim([min(histogramPostions) max(histogramPostions)]);
end




