
RICHRunLego()
{ 
   TControlBar *menu = new TControlBar("vertical","RICH Lego");
   menu->AddButton("Center chamber (40 bins)","gAlice->RunLego(\"Config.C\",40,80,100,40,80,100,0,500,600)","Fast run");
   menu->AddButton("All chambers (120 bins)","gAlice->RunLego(\"Config.C\",120,60,120,120,60,120,0,500,600)","Fast run");
   menu->AddButton("Center chamber (160 bins)","gAlice->RunLego(\"Config.C\",1600,80,100,160,80,100,0,500,600)","Slow run");
   menu->AddButton("All chambers hi (200 bins)","gAlice->RunLego(\"Config.C\",200,60,120,200,60,120,0,500,600)","Very slow run");
   menu->Show();
}



