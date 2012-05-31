/*
  Make default plotsfrom AliTPCdataQA components:

  aliroot -b -q  $ALICE_ROOT/TPC/CalibMacros/CalibQA.C\(121694\);
  
  .L $ALICE_ROOT/TPC/CalibMacros/CalibQA.C
  Int_t run=121694;
  CalibQA(run);
*/

TCut cutNoise="PadNoise.fElements<1.5&&abs(PadNoise.fElements/PadNoise_Median-1)<0.5";
TCut cutTime="abs(TimePosition.fElements-TimePosition_Median)<100";
TCut cutOccu="abs(NoThreshold.fElements/NoThreshold_Median-1)<0.9";
TCut cutAmp="abs(MaxCharge.fElements/MaxCharge_Median-1)<0.99";
TCut cutIROC="sector<36";
TCut cutOROC="sector>=36";


void CalibQA(Int_t run){
  InitOCDB(run);
  MakeTree();
  TCanvas *canvas=0;
  //
  TPostScript *ps = new TPostScript("rawQA.ps", 112);  
  ps->NewPage();
  canvas=DrawOccupancy();
  ps->NewPage();
  canvas->Update();
  ps->Close();
  delete ps;

}

void InitOCDB(Int_t run){
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/scripts/OCDBscan/ConfigOCDBLustre.C");
  //gROOT->LoadMacro("$ALICE_ROOT/TPC/scripts/OCDBscan/ConfigOCDB.C");
  gROOT->Macro("$ALICE_ROOT/TPC/scripts/OCDBscan/NimStyle.C");
  ConfigOCDB(run);
}

void MakeTree(){
  //
  // make summary tree
  //
  AliTPCcalibDB::Instance()->UpdateNonRec();
  AliTPCdataQA* dataQA =   AliTPCcalibDB::Instance()->GetDataQA();
  AliTPCCalPad* gain   =   AliTPCcalibDB::Instance()->GetDedxGainFactor();
  AliTPCCalPad * padNoise = AliTPCcalibDB::Instance()->GetPadNoise();
  AliTPCPreprocessorOnline preprocesor;
  gain->SetName("krGain");
  preprocesor.AddComponent(gain);
  preprocesor.AddComponent(dataQA->GetNPads());
  preprocesor.AddComponent(dataQA->GetNTimeBins());
  preprocesor.AddComponent(dataQA->GetMaxCharge());
  preprocesor.AddComponent(dataQA->GetNoThreshold());
  preprocesor.AddComponent(dataQA->GetNLocalMaxima());
  preprocesor.AddComponent(dataQA->GetTimePosition());
  preprocesor.AddComponent(padNoise);
  preprocesor.DumpToFile("QA.root");
}


TCanvas * DrawOccupancy(){
  TH1::AddDirectory(0);
  gStyle->SetOptStat(0);
  TFile f("QA.root");
  TTree  * tree = (TTree*)f.Get("calPads");
  TLegend *legend=0;
  TProfile * phoccGainIROC=0;
  TProfile * phoccGainOROC=0;
  //
  TCanvas * canvas = new TCanvas("occupancy","occupancy",700,700);
  canvas->Divide(2,2);
  canvas->cd(1);
  tree->Draw("NoThreshold.fElements:gy.fElements:gx.fElements>>hisOccuA(250,-250,250,250,-250,250)",cutNoise+cutTime+cutOccu+"sector%36<18","profcolz");
  canvas->cd(2);
  tree->Draw("NoThreshold.fElements:gy.fElements:gx.fElements>>hisOccuC(250,-250,250,250,-250,250)",cutNoise+cutTime+cutOccu+"sector%36>=18","profcolz");

  canvas->cd(3);
  tree->Draw("NoThreshold.fElements:krGain.fElements>>hoccGainIROC(20,0.7,1.2)",cutNoise+cutTime+cutOccu+cutIROC,"prof");
  legend = new TLegend(0.45,0.15,0.85,0.35, "Raw cluster occupancy");
  phoccGainIROC = (TProfile*)(gROOT->FindObject("hoccGainIROC")->Clone());
  phoccGainIROC->Draw();
  phoccGainIROC->SetTitle("IROC");
  phoccGainIROC->SetName("IROC");
  phoccGainIROC->GetXaxis()->SetTitle("Krypton Amp (a.u.)");
  legend->AddEntry(phoccGainIROC);
  legend->Draw();

  canvas->cd(4);
  tree->Draw("NoThreshold.fElements:krGain.fElements>>hoccGainOROC(20,0.8,1.2)",cutNoise+cutTime+cutOccu+cutOROC+"lx.fElements<200","prof");
  phoccGainOROC = (TProfile*)(gROOT->FindObject("hoccGainOROC")->Clone());
  phoccGainOROC->Draw();
  phoccGainOROC->SetTitle("OROC");
  phoccGainOROC->SetName("OROC");
  phoccGainOROC->GetXaxis()->SetTitle("Krypton Amp (a.u.)");
  legend = new TLegend(0.45,0.15,0.85,0.35, "Raw cluster occupancy");
  legend->AddEntry(phoccGainOROC);
  legend->Draw();
  return canvas;
}

TCanvas * DrawGain(){
  //
  // Compare the amplitude with krypton gain amplitude
  // Similar filtering as in occupancy plot
  //						
  return 0;
}
