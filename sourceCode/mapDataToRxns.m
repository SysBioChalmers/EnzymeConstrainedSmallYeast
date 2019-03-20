function model = mapDataToRxns(model, rxnFile)
%Load specific activity data from a file
%maps the data to the model reactions.
%for reactions with specific activty -1 set the median
%for reactions without data do nothing
SAData = importdata(rxnFile);
rxns = SAData.textdata;
specificActivity = SAData.data;

results = -ones(length(model.eccodes), 1);
medianRate = median(specificActivity);

for i = 1:length(rxns)
    currentRxn = ismember(model.rxns, rxns{i});
    
    if specificActivity(i) == -1
        results(currentRxn) = medianRate;
    else
        results(currentRxn) = specificActivity(i);
    end
end

model.specificActivity = results;
end

