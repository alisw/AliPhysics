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
   
//  right part view  :

//  gMC->Gdraw("alic", 0, 0, 0, 10, 10, 0.2, 0.2);

    gMC->Gdraw("alic", 90, 0, 0, -70, 10, 0.9, 0.2 );
    
    gMC->Gdhead(1111, "VZERO Detector");


}
