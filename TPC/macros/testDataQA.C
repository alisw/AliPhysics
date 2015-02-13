/// \file testDataQA.C
///
/// Simple macro to make pedstal calibration
/// and visualize calibration data
///
/// ~~~
/// .L $ALICE_ROOT/TPC/macros/testDataQA.C
/// ~~~

void testDataQA(Char_t *fileName, Int_t maxevent=10)
{
   AliRawReaderRoot *rawReader = new AliRawReaderRoot(fileName);
   if ( !rawReader ) return;
   AliTPCdataQA *calib     = new AliTPCdataQA;
   AliTPCdataQA *calibFast = new AliTPCdataQA;
   calib->SetRangeTime(1,1000);
   calibFast->SetRangeTime(1,1000);
   printf("Processing data\n");
   Int_t event=0;
   while (rawReader->NextEvent()){
     calib->ProcessEvent(rawReader);
     rawReader->Reset();
     event++;
     // if you wann to check the handling of Analyse updates uncomment this
//      if(gRandom->Rndm()<0.3) {
//        calib->Analyse();
//        calibFast->Analyse();
//      }
     if (event>maxevent) break;
   }
   calib->Analyse();
   calibFast->Analyse();

   TFile file1("dataQATestAnalyse.root","recreate");
   calib->Write("AliTPCdataQA");
   calibFast->Write("AliTPCdataQAFast");
   delete rawReader;
   //   delete calib;
   delete calibFast;

   file1.Close();
   //
   TFile file2("dataQATestAnalyse.root");
   AliTPCdataQA* cal = (AliTPCdataQA*)file2.Get("AliTPCdataQA");
   AliTPCdataQA* calFast = (AliTPCdataQA*)file2.Get("AliTPCdataQAFast");

   AliTPCPreprocessorOnline preprocesor;
   preprocesor.AddComponent(cal->GetNoThreshold());
   preprocesor.AddComponent(cal->GetOverThreshold10());
   preprocesor.AddComponent(cal->GetOverThreshold20());
   preprocesor.AddComponent(cal->GetOverThreshold30());
   preprocesor.AddComponent(cal->GetNLocalMaxima());
   preprocesor.AddComponent(cal->GetMeanCharge());
   preprocesor.AddComponent(cal->GetMaxCharge());
   preprocesor.AddComponent(cal->GetNTimeBins());
   preprocesor.AddComponent(cal->GetNPads());
   preprocesor.AddComponent(cal->GetTimePosition());

   AliTPCCalPad* calPadNoThr= new AliTPCCalPad(*calFast->GetNoThreshold());
   AliTPCCalPad* calPadThr10= new AliTPCCalPad(*calFast->GetOverThreshold10());
   AliTPCCalPad* calPadThr20= new AliTPCCalPad(*calFast->GetOverThreshold20());
   AliTPCCalPad* calPadThr30= new AliTPCCalPad(*calFast->GetOverThreshold30());
   AliTPCCalPad* calNLocal  = new AliTPCCalPad(*calFast->GetNLocalMaxima());
   AliTPCCalPad* calPadMean = new AliTPCCalPad(*calFast->GetMeanCharge());
   AliTPCCalPad* calPadMax  = new AliTPCCalPad(*calFast->GetMaxCharge());
   AliTPCCalPad* calNTime   = new AliTPCCalPad(*calFast->GetNTimeBins());
   AliTPCCalPad* calNPad    = new AliTPCCalPad(*calFast->GetNPads());
   AliTPCCalPad* calTimePos = new AliTPCCalPad(*calFast->GetTimePosition());
   //
   calPadNoThr ->SetName(Form("%sFast", calPadNoThr->GetName()));
   calPadThr10 ->SetName(Form("%sFast", calPadThr10->GetName()));
   calPadThr20 ->SetName(Form("%sFast", calPadThr20->GetName()));
   calPadThr30 ->SetName(Form("%sFast", calPadThr30->GetName()));
   calNLocal   ->SetName(Form("%sFast", calNLocal  ->GetName()));
   calPadMean  ->SetName(Form("%sFast", calPadMean ->GetName()));
   calPadMax   ->SetName(Form("%sFast", calPadMax  ->GetName()));
   calNTime    ->SetName(Form("%sFast", calNTime   ->GetName()));
   calNPad     ->SetName(Form("%sFast", calNPad    ->GetName()));
   calTimePos  ->SetName(Form("%sFast", calTimePos ->GetName()));

   preprocesor.AddComponent(calPadNoThr );
   preprocesor.AddComponent(calPadThr10 );
   preprocesor.AddComponent(calPadThr20 );
   preprocesor.AddComponent(calPadThr30 );
   preprocesor.AddComponent(calNLocal   );
   preprocesor.AddComponent(calPadMean  );
   preprocesor.AddComponent(calPadMax   );
   preprocesor.AddComponent(calNTime    );
   preprocesor.AddComponent(calNPad     );
   preprocesor.AddComponent(calTimePos  );

   preprocesor.DumpToFile("CalibTreeTestAnalyse.root");
   AliTPCCalibViewerGUI::ShowGUI("CalibTreeTestAnalyse.root");
}

