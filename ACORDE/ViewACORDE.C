void ViewACORDE()
{
   gMC->Gsatt("ALIC","seen",0);

   gMC->Gsatt("L3MO","seen",0); // L3 Magnet, Mother
   gMC->Gsatt("L3CO","seen",1); // Coils
   gMC->Gsatt("L3C1","seen",1); // Coils
   gMC->Gsatt("L3YO","seen",1); // Yoke
   gMC->Gsatt("L3DO","seen",0); // return Yoke (DOOR)
   gMC->Gsatt("L3FR","seen",1); // DOOR
   gMC->Gsatt("L3IR","seen",0); // Inner layer
   gMC->Gsatt("L3O1","seen",1); // Door opening
   gMC->Gsatt("L3O2","seen",1); // Door opening

   gMC->Gsatt("ACORDE1", "seen", 0); // ACORDE Mother
   gMC->Gsatt("ACORDE2", "seen", 0); // Module air box
   gMC->Gsatt("ACORDE3", "seen", 1); // Scintillators
   gMC->Gsatt("ACORDE3", "colo", 2); // Scintillators
   gMC->Gsatt("ACORDE4", "seen", 1); // Aluminium frame (long bars)
   gMC->Gsatt("ACORDE4", "colo", 3); //
   gMC->Gsatt("ACORDE5", "seen", 1); // Aluminium frame (short bars)
   gMC->Gsatt("ACORDE5", "colo", 3); //
   gMC->Gsatt("ACORDE6", "seen", 1); // Module support
   gMC->Gsatt("ACORDE6", "colo", 3); //
}
