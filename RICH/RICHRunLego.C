
RICHRunLego()
{ 
   TControlBar *menu = new TControlBar("vertical","RICH Lego");
   menu->AddButton("Center chamber (40 bins)","gAlice->RunLego(\"Config.C\",40,80,100,40,80,100,0,600,500)","Fast run");
   menu->AddButton("All chambers (120 bins)","gAlice->RunLego(\"Config.C\",120,60,120,120,60,120,0,600,500)","Fast run");
   menu->AddButton("Center chamber (80 bins)","gAlice->RunLego(\"Config.C\",80,80,100,80,80,100,0,510,500)","Slow run");
   menu->AddButton("All chambers hi (200 bins)","gAlice->RunLego(\"Config.C\",200,60,120,200,60,120,0,600,500)","Very slow run");
   menu->AddButton("Central chamber, x-z (40 bins)","gAlice->RunLego(\"Config.C\", 40, 0, 80, 40, 0, 80, 0, 600, 500, new  AliLegoGeneratorXYZ(\"y\"))","Fast run");
   menu->AddButton("Central chamber, x-z (100 bins)","gAlice->RunLego(\"Config.C\", 100, 0, 80, 100, 0, 80, 0, 510, 500, new AliLegoGeneratorXYZ(\"y\"))","Slow run");
   menu->Show();
}



