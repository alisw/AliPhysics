Int_t events;

RICHreco(Int_t nev=1)
{ 

  events=nev;
  
   TControlBar *menu = new TControlBar("vertical","RICH reconstruction");
   menu->AddButton("1D Hough Pattern Recognition",      ".x RICHpatrec.C(0,events-1)","Bari");
   menu->AddButton("3D Hough Pat. Rec. v0 from digits",      ".x RICHdetect.C(0,events-1, 0, 0)","Lisbon");
   menu->AddButton("3D Hough Pat. Rec. v0 from clusters",      ".x RICHdetect.C(0,events-1, 1, 0)","Lisbon");
   menu->AddButton("3D Hough Pat. Rec. v1 from digits",      ".x RICHdetect.C(0,events-1, 0, 1)","Lisbon");
   menu->AddButton("3D Hough Pat. Rec. v1 from clusters",      ".x RICHdetect.C(0,events-1, 1, 1)","Lisbon");
   menu->Show();
}
