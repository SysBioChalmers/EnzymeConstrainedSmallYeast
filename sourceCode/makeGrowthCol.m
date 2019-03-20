function newCol = makeGrowthCol(colLength, positionObject, valueObject)
    %add biomass equation
   newCol = zeros(colLength, 1);
   newCol(positionObject.protein) = -valueObject.protein;
   newCol(positionObject.rna) = -valueObject.rna;
   newCol(positionObject.dna) = -valueObject.dna;
   newCol(positionObject.lipid) = -valueObject.lipid;
   newCol(positionObject.glycogen) = -valueObject.glycogen;
   newCol(positionObject.trehalose) = -valueObject.trehalose;
   newCol(positionObject.mannan) = -valueObject.mannan;
   newCol(positionObject.glucan) = -valueObject.glucan;
   newCol(positionObject.maintain) = -valueObject.maintain;
   newCol(positionObject.biomass) = valueObject.biomass;
end

