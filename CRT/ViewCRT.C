void ViewCRT()
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

   gMC->Gsatt("CRT1", "seen", 0); // CRT Mother
   gMC->Gsatt("CRT2", "seen", 0); // Module air box
   gMC->Gsatt("CRT3", "seen", 1); // Scintillators
   gMC->Gsatt("CRT3", "colo", 2); // Scintillators
   gMC->Gsatt("CRT4", "seen", 1); // Aluminium frame (long bars)
   gMC->Gsatt("CRT4", "colo", 3); //
   gMC->Gsatt("CRT5", "seen", 1); // Aluminium frame (short bars)
   gMC->Gsatt("CRT5", "colo", 3); //
   gMC->Gsatt("CRT6", "seen", 1); // Module support
   gMC->Gsatt("CRT6", "colo", 3); //
}
