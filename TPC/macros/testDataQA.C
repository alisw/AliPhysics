/*
  Simple macro to make pedstal calibration
  and visualize calibration data

  .L $ALICE_ROOT/TPC/macros/testDataQA.C
  
*/




void testDataQA(Char_t *fileName, Int_t maxevent=10)
{
   AliRawReaderRoot *rawReader = new AliRawReaderRoot(fileName);
   if ( !rawReader ) return;
   AliTPCdataQA *calib = new AliTPCdataQA; 
   calib->SetRangeTime(200,800);
   AliTPCCalibPulser *calibPulser = new AliTPCCalibPulser; 
   calibPulser->SetRangeTime(5,300);
   printf("Processing data\n");
   Int_t event=0;
   while (rawReader->NextEvent()){
     calib->ProcessEvent(rawReader);
     rawReader->Reset();
     calibPulser->ProcessEvent(rawReader);
     printf("Processing event\t%d\n",event);
     event++;
     if (event>maxevent) break;
   }
   calibPulser->Analyse();
   TFile f("dataQA.root","recreate");
   calibPulser->Write("Pulser");
   calib->Write("AliTPCdataQA");
   delete rawReader;
   delete calib;
   delete calibPulser;

   f.Close();
   //
   TFile f2("dataQA.root");
   AliTPCdataQA* ped = (AliTPCdataQA*)f2.Get("AliTPCdataQA");
   AliTPCCalibPulser* pedPulser = (AliTPCCalibPulser*)f2.Get("Pulser");
   
   AliTPCPreprocessorOnline preprocesor;
   preprocesor.AddComponent(ped->GetMaxCharge());
   preprocesor.AddComponent(ped->GetOverThreshold0());
   preprocesor.AddComponent(ped->GetOverThreshold5());
   preprocesor.AddComponent(ped->GetOverThreshold10());
   preprocesor.AddComponent(ped->GetOverThreshold20());
   preprocesor.AddComponent(ped->GetOverThreshold30());
   //
   AliTPCCalPad * pad0 = new AliTPCCalPad(pedPulser->GetCalPadT0());
   AliTPCCalPad * pad1 = new AliTPCCalPad(pedPulser->GetCalPadQ());
   AliTPCCalPad * pad2 = new AliTPCCalPad(pedPulser->GetCalPadRMS());
   pad0->SetName("PulserT0");
   pad1->SetName("PulserQ");
   pad2->SetName("PulserRMS");
   
   preprocesor.AddComponent(pad0);
   preprocesor.AddComponent(pad1);
   preprocesor.AddComponent(pad2);
   
   preprocesor.DumpToFile("CalibTree2.root");
   AliTPCCalibViewerGUI::ShowGUI("CalibTree2.root");
}

