
Int_t events;

RICHDiagnostics(Int_t nev=1)
{ 

  events=nev;
  
   TControlBar *menu = new TControlBar("vertical","RICH diagnostics");
   menu->AddButton("Single Ring Hits",".x RICHpadtest.C(1,0,events-1","Hits in central region");
   menu->AddButton("Single Ring Spectra",".x RICHpadtest.C(2,0,events-1","Photon spectra");
   menu->AddButton("Single Ring Statistics",".x RICHpadtest.C(3,0,events-1","Production and clusters");
   menu->AddButton("Single Ring Reconstruction",".x RICHpadtest.C(4,0,events-1","Generated and reconstructed values");
   menu->AddButton("Full Event Hits",".x RICHpadtest.C(5,0,events-1","Hits in seven modules");
   menu->AddButton("Full Event Spectra",".x RICHpadtestC.C","Individual particles' spectra and fluxes");
   menu->AddButton("Full Event Occupancy",".x RICHoccupancy.C","Mean and per chamber occupancy");
   menu->Show();
}



