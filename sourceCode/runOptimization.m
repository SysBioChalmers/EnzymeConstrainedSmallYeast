 function resX = runOptimization(model, valueObject)
   biomassLocations.protein = findIndex(model.mets, 'growthProtein[c]');
   biomassLocations.rna = findIndex(model.mets, 'growthRNA[c]');
   biomassLocations.dna = findIndex(model.mets, 'growthDNA[c]');
   biomassLocations.lipid = findIndex(model.mets, 'growthLipid[c]');
   biomassLocations.glycogen = findIndex(model.mets, 'growthGlycogen[c]');
   biomassLocations.trehalose = findIndex(model.mets, 'growthTrehalose[c]');
   biomassLocations.mannan = findIndex(model.mets, 'growthMannan[c]');
   biomassLocations.glucan = findIndex(model.mets, 'growthGlucan[c]');
   biomassLocations.maintain = findIndex(model.mets, 'growthMaintainance[c]');
   biomassLocations.biomass = findIndex(model.mets, 'BIOMASS[c]');
   
   %Set biomass equation
   biomassEquation = findIndex(model.rxns, 'GROWTH');
   growthCol = makeGrowthCol(size(model.S, 1), biomassLocations, valueObject);
   model.S(:,biomassEquation) = growthCol;

   %Run linear optimization
   output = solveLinMin(model);
   if output.x == 0
       resX = zeros(size(model.S,2),1);
   else
       resX = output.x;
   end
end

