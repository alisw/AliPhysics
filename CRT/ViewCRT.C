void ViewCRT()
{
   gMC->Gsatt("ALIC","seen",0);
   gMC->Gsatt("L3MO","seen",1); // L3 Magnet
   gMC->Gsatt("CRT1","seen",1); // Scintillators

   // Draw the molasse volumes
   gMC->Gsatt("CMO1","seen",0); // Exactly above the HALL
   gMC->Gsatt("CMO2","seen",0); // Molasse, along the PM25
   gMC->Gsatt("CMO3","seen",0); // molasse along the PGC2
   gMC->Gsatt("CMO4","seen",0); // Molasse, behind the PX24 upper part
   gMC->Gsatt("CMO5","seen",0); // molasse behind px24, lower part
   gMC->Gsatt("CMO6","seen",0); // behind the PX24
   gMC->Gsatt("CMO7","seen",0); // behind the PGC2
   gMC->Gsatt("CMO8","seen",0); // on the right side.
   gMC->Gsatt("CMO9","seen",0); // on the left side.
   gMC->Gsatt("CM10","seen",0); // betwen PX24 & PM25.
   gMC->Gsatt("CM11","seen",0); // betwen PGC2 & PM25.
   gMC->Gsatt("CM12","seen",0); // box above the hall.

}
