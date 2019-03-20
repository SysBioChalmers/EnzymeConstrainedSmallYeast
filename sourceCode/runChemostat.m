function results = runChemostat(model, growthRate, valueObject, resultIndex)
    results = zeros(length(growthRate), length(resultIndex));
    for i=1:length(growthRate)
        model = setParam(model,'lb',{'GROWTH'}, growthRate(i));
        model = setParam(model,'ub',{'GROWTH'}, growthRate(i));

        resX1 = runOptimization(model, valueObject);    
        results(i,:) = resX1(resultIndex);
    end
end

