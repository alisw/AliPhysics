
Int_t events;

mrich(Int_t nev=1)
{ 

  events=nev;
  
   TControlBar *menu = new TControlBar("vertical","RICH diagnostics");
   menu->AddButton("Diagnostics",       ".x RICHpadtest(3,0,events-1","Miscellaneous diagnostics");
   menu->AddButton("Quit",             ".q","Close session");
   
   menu->Show();
}



