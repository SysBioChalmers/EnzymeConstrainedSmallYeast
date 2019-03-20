function valueObject = makeValueObjectWeight(prot, rna, dna, lip, gly, tre, man, glu, main, bio)
   %Set up biomass equation give values in g/g
   valueObject.protein = prot/0.108;
   valueObject.rna = rna/0.324;
   valueObject.dna = dna/0.324;
   valueObject.lipid = lip/0.6;
   valueObject.glycogen = gly/0.162;
   valueObject.trehalose = tre/0.3423;
   valueObject.mannan = man/0.162;
   valueObject.glucan = glu/0.162;
   valueObject.maintain = main;
   valueObject.biomass = bio; 
end

