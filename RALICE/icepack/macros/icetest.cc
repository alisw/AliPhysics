//////////////////////////////////////////////////////////////////////////
// ROOT macro to show the functionality of data handling in IcePack.
// To run this ROOT macro in batch mode and get the output in a file
// test.log, just type the following at the regular command prompt ($)
//
// $ root -b -q icetest.cc >test.log 
//
// Interactive execution of this macro can be achieved by
// the following command after starting a ROOT session
// 
// $ root> .x icetest.cc
//  
// The commented-out lines at the end of the macro show how
// to obtain a 3D event view of the stored data.
// Obviously such an event display can only be invoked by
// executing this test macro in interactive mode.
//
// NvE 01-dec-2004 Utrecht University
//////////////////////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 
 IceEvent* evt=new IceEvent();
 evt->SetOwner();

 // Amanda module
 IceAOM m;
 m.SetUniqueID(123);
 m.SetNameTitle("OM123","Amanda module");

 Float_t pos[3]={1,2,3};
 m.SetPosition(pos,"car");

 // The starting unique signal ID.
 // It will be increased automatically in this macro
 // when a new signal is created.
 Int_t sid=1;

 AliSignal s;

 s.SetSlotName("ADC",1);
 s.SetSlotName("LE",2);
 s.SetSlotName("TOT",3);

 s.Reset();
 s.SetName("OM123 Hit 1");
 s.SetUniqueID(sid++);
 s.SetSignal(100,"ADC");
 s.SetSignal(-100,"LE");
 s.SetSignal(-1000,"TOT");
 m.AddHit(s);

 s.Reset();
 s.SetName("OM123 Hit 2");
 s.SetUniqueID(sid++);
 s.SetSignal(110,"ADC");
 s.SetSignal(-101,"LE");
 s.SetSignal(1001,"TOT");
 m.AddHit(s);

 s.Reset();
 s.SetName("OM123 Hit 3");
 s.SetUniqueID(sid++);
 s.SetSignal(120,"ADC");
 s.SetSignal(-102,"LE");
 s.SetSignal(-1002,"TOT");
 m.AddHit(s);

 evt->AddDevice(m);

 m.Reset();
 m.SetUniqueID(456);
 m.SetName("OM456");

 pos[0]=-4;
 pos[1]=-5;
 pos[2]=-6;
 m.SetPosition(pos,"car");

 s.Reset();
 s.SetName("OM456 Hit 1");
 s.SetUniqueID(sid++);
 s.SetSignal(20,"ADC");
 s.SetSignal(-200,"LE");
 s.SetSignal(-2000,"TOT");
 m.AddHit(s);

 s.Reset();
 s.SetName("OM456 Hit 2");
 s.SetUniqueID(sid++);
 s.SetSignal(21,"ADC");
 s.SetSignal(-201,"LE");
 s.SetSignal(2001,"TOT");
 m.AddHit(s);

 s.Reset();
 s.SetName("OM456 Hit 3");
 s.SetUniqueID(sid++);
 s.SetSignal(22,"ADC");
 s.SetSignal(-202,"LE");
 s.SetSignal(-2002,"TOT");
 m.AddHit(s);

 // Example of explicit hit selection
 AliSignal* sx=m.GetIdHit(5);
 if (sx)
 {
  cout << " === Hit selection on UniqueID=5 for OM 456 ===" << endl;
  sx->Data();
 }

 evt->AddDevice(m);

 m.Reset();
 m.SetUniqueID(558);
 m.SetName("OM558");

 pos[0]=5;
 pos[1]=5;
 pos[2]=8;
 m.SetPosition(pos,"car");

 s.Reset();
 s.SetName("OM558 Hit 1");
 s.SetUniqueID(sid++);
 s.SetSignal(30,"ADC");
 s.SetSignal(-300,"LE");
 s.SetSignal(-3000,"TOT");
 m.AddHit(s);

 s.Reset();
 s.SetName("OM558 Hit 2");
 s.SetUniqueID(sid++);
 s.SetSignal(31,"ADC");
 s.SetSignal(-301,"LE");
 s.SetSignal(3001,"TOT");
 m.AddHit(s);

 s.Reset();
 s.SetName("OM558 Hit 3");
 s.SetUniqueID(sid++);
 s.SetSignal(32,"ADC");
 s.SetSignal(-302,"LE");
 s.SetSignal(-3002,"TOT");
 m.AddHit(s);

 evt->AddDevice(m);

 // IceCube in-ice DOM
 IceIDOM mid;
 mid.SetUniqueID(958);
 mid.SetNameTitle("OM958","IceCube in-ice module");

 pos[0]=9;
 pos[1]=5;
 pos[2]=8;
 mid.SetPosition(pos,"car");

 s.Reset();
 s.SetName("OM958 Hit 1");
 s.SetUniqueID(sid++);
 s.SetSignal(40,"ADC");
 s.SetSignal(-400,"LE");
 s.SetSignal(-4000,"TOT");
 mid.AddHit(s);

 s.Reset();
 s.SetName("OM958 Hit 2");
 s.SetUniqueID(sid++);
 s.SetSignal(41,"ADC");
 s.SetSignal(-401,"LE");
 s.SetSignal(4001,"TOT");
 mid.AddHit(s);

 s.Reset();
 s.SetName("OM958 Hit 3");
 s.SetUniqueID(sid++);
 s.SetSignal(42,"ADC");
 s.SetSignal(-402,"LE");
 s.SetSignal(-4002,"TOT");
 mid.AddHit(s);

 evt->AddDevice(mid);

 // IceTop DOM
 IceTDOM mtd;
 mtd.SetUniqueID(4958);
 mtd.SetNameTitle("OM4958","IceTop module");

 pos[0]=49;
 pos[1]=5;
 pos[2]=8;
 mtd.SetPosition(pos,"car");

 s.Reset();
 s.SetName("OM4958 Hit 1");
 s.SetUniqueID(sid++);
 s.SetSignal(50,"ADC");
 s.SetSignal(-500,"LE");
 s.SetSignal(-5000,"TOT");
 mtd.AddHit(s);

 s.Reset();
 s.SetName("OM4958 Hit 2");
 s.SetUniqueID(sid++);
 s.SetSignal(51,"ADC");
 s.SetSignal(-501,"LE");
 s.SetSignal(5001,"TOT");
 mtd.AddHit(s);

 s.Reset();
 s.SetName("OM4958 Hit 3");
 s.SetUniqueID(sid++);
 s.SetSignal(52,"ADC");
 s.SetSignal(-502,"LE");
 s.SetSignal(-5002,"TOT");
 mtd.AddHit(s);

 evt->AddDevice(mtd);

 // Provide event data overview
 evt->Data();

 // Select a specific device (i.e. OM) from the event
 AliDevice* dx=(AliDevice*)evt->GetIdDevice(958);
 if (dx)
 {
  cout << " === Device selection on UniqueID=958 for the whole event ===" << endl;
  dx->Data();
 }

 // Select a specific hit from the event
 AliSignal* sx=evt->GetIdHit(5,"IceGOM");
 if (sx)
 {
  cout << " === Hit selection on UniqueID=5 for the whole event via IceGOM ===" << endl;
  sx->Data();
 }

 // Dump all the information for the various stored devices
 cout << endl;
 cout << " ======= devices dump ========" << endl;
 cout << endl;

 Int_t ndev=evt->GetNdevices();
 for (Int_t idev=1; idev<=ndev; idev++)
 {
  cout << " Device number : " << idev << endl;
  IceGOM* om=(IceGOM*)evt->GetDevice(idev);
  if (om) om->Data();
 }

 // Dump all the information for the various stored hits
 cout << endl;
 cout << " ======= Event all hits dump ========" << endl;
 cout << endl;

 // Obtain pointers to the hits for all generic OM's (i.e. IceGOM)
 TObjArray* hits=evt->GetHits("IceGOM");
 Int_t nhits=0;
 if (hits) nhits=hits->GetEntries();
 for (Int_t ih=0; ih<nhits; ih++)
 {
  AliSignal* sx=(AliSignal*)hits->At(ih);
  if (sx) sx->Data();
 }

 // Obtain the minimum and maximum observed TOT value 
 Float_t vmin,vmax;
 evt->GetExtremes("IceGOM",vmin,vmax,"TOT");

 cout << endl;
 cout << " === Extreme values : vmin = " << vmin << " vmax = " << vmax << endl;
 cout << endl;

 // Various hit orderings
 cout << endl;
 cout << " ======= ordered hits w.r.t. decreasing TOT ========" << endl;
 cout << endl;

 TObjArray* ordered=evt->SortHits("IceGOM","TOT",-1);
 nhits=0;
 if (ordered) nhits=ordered->GetEntries();
 for (Int_t i=0; i<nhits; i++)
 {
  AliSignal* sx=(AliSignal*)ordered->At(i);
  if (sx) sx->Data();
 }

 cout << endl;
 cout << " ======= ordered devices from the ordered hit array  ========" << endl;
 cout << endl;

 TObjArray* devs=evt->SortDevices(ordered,0,0);
 ndev=0;
 if (devs) ndev=devs->GetEntries();
 for (Int_t id=0; id<ndev; id++)
 {
  AliDevice* dx=(AliDevice*)devs->At(id);
  if (dx) dx->Data();
 }

 cout << endl;
 cout << " ======= newly ordered devices w.r.t. decreasing ADC ========" << endl;
 cout << endl;

 TObjArray* devs=evt->SortDevices("IceGOM","ADC",-1);
 ndev=0;
 if (devs) ndev=devs->GetEntries();
 for (Int_t id2=0; id2<ndev; id2++)
 {
  AliDevice* dx=(AliDevice*)devs->At(id2);
  if (dx) dx->Data();
 }

 // Example for a 3D signal display of the devices
/***********
 TCanvas* c1=new TCanvas("c1","c1");
 c1->x3d();
 TView* view=new TView(1);
 view->SetRange(-50,-50,-50,50,50,50);
 view->ShowAxis();

 evt->DisplayHits("IceGOM","TOT",1e4,1);
************/
}
