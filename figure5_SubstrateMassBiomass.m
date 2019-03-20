addpath('sourceCode')
load('data/model.mat');

%Load specific activity data.
model = mapDataToRxns(model, 'data/RxnAndSA.txt');

%Add a mass constraint metabolite to each enzymatic reaction in the S-matrix
model = addSpecificActivityConstraint(model, 0.5, 0.1, 60);

%Make the S matrix strictly positive.
model = addReversedReactions(model);

%
gluIn = 1000;
O2In = 1000;


addpath('sourceCode')

%Transport reactions moving pyrovate, oxaloacetate, ACCOA, NADH and ATP
%from mit are added to allow biomass to be formed for growth on acetate and
%ethanol (a simplification that creates som weird options for proton
%decoupling), should still be a good estimate
%load('data/modelWithBackTransport.mat')

weightVector = full(model.S(end, :));

interval = 200;
results = zeros(interval, 5);

proteinLimit = 0.1;

model.b(:,1) = 0;
model.b(:,2) = proteinLimit;


model = setParam(model, 'ub', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 1000);
model = setParam(model, 'lb', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, -1000);

growthRate = linspace(0.05, 0.5, interval);

experiments = {'acIN', 'ethIN', 'galIN', 'gluActiveIn', 'glcIN'};
valueObject = makeValueObjectWeight(0.4, 0.10, 0.006, 0.025, 0.005, 0.01, 0.1, 0.40, 40, 1);
%valueObject = makeValueObjectWeight(0, 0, 0, 0, 0, 0, 0, 0, 90, 1);

%model = setParam(model,'lb',{'ATPX'}, 0);

for j = 1:length(experiments)
    model = setParam(model,'ub',{'glcIN', 'galIN', 'ethIN', 'gluActiveIn', 'acIN', 'o2IN'},[0, 0, 0, 0, 0, O2In]);
    model = setParam(model,'lb',{'GROWTH'}, 0);
    model = setParam(model,'ub',{'GROWTH'}, 1000);
    model = setParam(model,'ub', experiments(j), [gluIn]);
    model = setParam(model,'obj', experiments(j), -1);
    
    for i = 1:length(growthRate)
        model = setParam(model,'lb',{'GROWTH'}, growthRate(i));
        model = setParam(model,'ub',{'GROWTH'}, growthRate(i));
        resX = runOptimization(model, valueObject);
        totalWeight = sum(weightVector' .* resX);    
        results(i, j) = totalWeight/growthRate(i);   
    end
end

plotXvalues = growthRate;

clf
hold all

colorOfLines = {'g', 'r', 'b', 'm', 'k'};

for i = 1:length(experiments)
    currentValues = plotXvalues(results(:,i)>0)';
    currentValues = [0; currentValues];
    currentResults = 1./results(results(:, i)>0,i);
    currentResults = [min(currentResults); currentResults];
    plot(currentValues,  currentResults, colorOfLines{i}, 'LineWidth',3)
end

for i = 1:length(experiments)
    currentValues = plotXvalues(results(:,i)>0)';
    currentResults = results(results(:, i)>0,i);
    currentResults = max(1./currentResults);
    plot([currentValues(end) currentValues(end)],  [currentResults 0], [colorOfLines{i} '--'], 'LineWidth',2);
end

plot([0 0.5 0 0], [0 5 5 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth',2);


legend('Acetate', 'Ethanol', 'Galactose', 'Active Glucose', 'Facilitated Glucose', 'Location', 'NW')
xlabel('Growth Rate [h-1]', 'FontSize',15,'FontName', 'Arial')
ylabel('Protein Efficiency [gdw h-1 g-1]', 'FontSize',15,'FontName', 'Arial')
xlim([0, 0.5]);
ylim([0, 5.01]);
%set(gca,'XTick', 0:0.05:0.45)
%set(gca,'YTick', 0:0.025:0.15)

set(gca,'FontSize',15,'FontName', 'Arial')

colorOfLines = [0 1 0
     1 0 0
     0 0 1
     0 0 0];

growthRates = [0.17 0.13 0.28 0.41];
for i = 1:4
    plot(growthRates(i), 0.06, 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colorOfLines(i,:));
end
hold off

for i = 1:length(experiments)
    maxGrowth = growthRate(results(:,i)>0);
    maxGrowth = maxGrowth(end);
    fprintf('%2.2f\t%s\n', maxGrowth, experiments{i});
end

legend boxoff 
