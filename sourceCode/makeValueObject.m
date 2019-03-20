function valueObject = makeValueObject(prot, rna, dna, lip, gly, tre, man, glu, main, bio)
    %set up biomass equation give values in mmol/g
   valueObject.protein = prot;
   valueObject.rna = rna;
   valueObject.dna = dna;
   valueObject.lipid = lip;
   valueObject.glycogen = gly;
   valueObject.trehalose = tre;
   valueObject.mannan = man;
   valueObject.glucan = glu;
   valueObject.maintain = main;
   valueObject.biomass = bio;
end

