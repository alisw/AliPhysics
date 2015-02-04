/// \file testPed.C
///
/// Simple macro to make pedstal calibration
/// and visualize calibration data

void testPed(Char_t *fileName)
{
   AliRawReaderRoot *rawReader = new AliRawReaderRoot(fileName);
   if ( !rawReader ) return;
   AliTPCCalibPedestal *calib = new AliTPCCalibPedestal;
   printf("Processing data\n");
   Int_t event=0;
   while (rawReader->NextEvent()){
     calib->ProcessEvent(rawReader);
     printf("Processing event\t%d\n",event);
     event++;
   }
   calib->Analyse();
   calib->DumpToFile("PedestalData.root");
   delete rawReader;
   delete calib;

   //
   TFile f("PedestalData.root");
   AliTPCCalibPedestal* ped = (AliTPCCalibPedestal*)f.Get("AliTPCCalibPedestal");
   AliTPCCalPad * pad0 = new AliTPCCalPad(ped->GetCalPadRMS());
   AliTPCCalPad * pad1 = new AliTPCCalPad(ped->GetCalPadPedestal());
   pad0->SetName("Noise");
   pad1->SetName("Pedestal");

   AliTPCPreprocessorOnline preprocesor;
   preprocesor.AddComponent(pad0);
   preprocesor.AddComponent(pad1);
   preprocesor.DumpToFile("CalibTree.root");
   AliTPCCalibViewerGUI::ShowGUI("CalibTree.root");
}


void testPed0(Char_t *fileName)
{
  ///

  AliRawReaderRoot *rawReader = new AliRawReaderRoot(fileName);
  if ( !rawReader ) return;
  AliTPCCalibPedestal *calib = new AliTPCCalibPedestal;
  printf("Processing data\n");
  Int_t event=0;
  while (rawReader->NextEvent()){
    calib->ProcessEvent(rawReader);
    printf("Processing event\t%d\n",event);
    event++;
  }
  calib->Analyse();
  //
  AliTPCCalPad * pad0 = new AliTPCCalPad(calib->GetCalPadRMS());
  AliTPCCalPad * pad1 = new AliTPCCalPad(calib->GetCalPadPedestal());
  pad0->SetName("Noise");
  pad1->SetName("Pedestal");
  //
  AliTPCPreprocessorOnline preprocesor;
  preprocesor.AddComponent(pad0);
  preprocesor.AddComponent(pad1);
  preprocesor.DumpToFile("CalibTree.root");
  AliTPCCalibViewerGUI::ShowGUI("CalibTree.root");
}

