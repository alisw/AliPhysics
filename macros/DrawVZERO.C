void DrawVZERO()

{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewVZERO.C");
   gInterpreter->ProcessLine("ViewVZERO()");
   gMC->Gdopt("hide", "off");
   gMC->Gdopt("shad", "off");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->DefaultRange();
   
//   Vues de la partie gauche :

//   gMC->Gdraw("alic", 30, 0, 0, -25, 10, 0.2, 0.2);

//   Vues de la partie droite :

//  gMC->Gdraw("alic", 0, 0, 0, 10, 10, 0.2, 0.2);


    gMC->Gdraw("alic", 90, 0, 0, 90, 9, 0.92, 0.2 );
  
    gMC->Gdhead(1111, "VZERO Detector");


}
