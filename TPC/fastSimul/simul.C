/*

gROOT->LoadMacro("$ALICE_ROOT/TPC/fastSimul/AliTPCclusterFast.cxx+");
.L $ALICE_ROOT/TPC/fastSimul/simul.C
//Merge()

TFile f("mergetrack.root");
track = (AliTPCtrackFast*)f.Get("track");
//
// Draw debug stream
//
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
AliXRDPROOFtoolkit tool;
TChain * chain = tool.MakeChain("track.txt","simulTrack",0,200000);
chain->Lookup();
 gProof->Exec("gSystem->Load(\"$ALICE_ROOT/TPC/fastSimul/AliTPCclusterFast_cxx.so\")",kFALSE);

*/

class THnSparse;

void simul(Int_t npoints){ 
  //
  // simulation submit script
  //
  printf("Hallo world\n");
  gRandom->SetSeed(0);
  gROOT->LoadMacro("$ALICE_ROOT/TPC/fastSimul/AliTPCclusterFast.cxx+");
  AliTPCclusterFast::fPRF = new TF1("fprf","gausn",-5,5);
  AliTPCclusterFast::fTRF = new TF1("ftrf","gausn",-5,5);
  AliTPCclusterFast::fPRF->SetParameters(1,0,0.5);
  AliTPCclusterFast::fTRF->SetParameters(1,0,0.5);
  //
  AliTPCtrackFast::Simul("trackerSimul.root",npoints); 
}



void Merge(){
  //
  //
  //
  TString objfile;
  AliTPCtrackFast track0;
  track0.MakeHisto();
  AliTPCtrackFast *track1;
  ifstream in;
  Int_t counter=0;
  in.open("track.txt");
  while(in.good()) {
    in >> objfile;
    if (!objfile.Contains("root")) continue; // protection
    TFile currentFile(objfile.Data());
    printf("Open file:Counter\t%d\tMerging file %s\n",counter,objfile.Data());
    track1=(AliTPCtrackFast)currentFile.Get("track");
    if (!track1) continue;
    track0.Add(*track1);
    counter++;
  } 
  TFile f("mergetrack.root","recreate");
  track0.Write("track");
  f.Close("");
}



void DrawDedxMC(THnSparse * hstat){
  //
  //
  //
  TH1 * hisMean[7];
  TH1 * hisSigma[7];
  TObjArray arr;
  for (Int_t ifrac=0; ifrac<6; ifrac++){
    Float_t frac = 0.5+0.1*Float_t(ifrac);
    hstat->GetAxis(3)->SetRange(ifrac+1,ifrac+1);
    hstat->GetAxis(2)->SetRangeUser(120,160);
    TH2F * his = (TH2F*)hstat->Projection(0,1);
    his->FitSlicesY(0,0,-1,0,"QNR",&arr);
    delete his;
    hisMean[ifrac]  = (TH1*) arr.At(1)->Clone();
    hisSigma[ifrac] = (TH1*) arr.At(2)->Clone();
    arr.SetOwner(kTRUE); arr.Delete();
    //
    hisSigma[ifrac]->Divide(hisMean[ifrac]);
    hisMean[ifrac]->SetMaximum(6);
    hisMean[ifrac]->SetMinimum(0);
    hisSigma[ifrac]->SetMaximum(0.07);
    hisSigma[ifrac]->SetMinimum(0.03);
    //
    hisMean[ifrac]->SetDirectory(0);
    hisSigma[ifrac]->SetDirectory(0);
    hisMean[ifrac]->SetXTitle("N_{prim}");
    hisSigma[ifrac]->SetXTitle("N_{prim}");
    hisMean[ifrac]->SetYTitle("Q/N_{prim}");
    hisSigma[ifrac]->SetYTitle("#sigma_{Q/N_{prim}}/(Q/N_{prim})");
    hisMean[ifrac]->SetMarkerColor(kmicolors[ifrac+1]);
    hisMean[ifrac]->SetMarkerStyle(kmimarkers[ifrac+1]);
    hisSigma[ifrac]->SetMarkerColor(kmicolors[ifrac+1]);
    hisSigma[ifrac]->SetMarkerStyle(kmimarkers[ifrac+1]);
  }
  TCanvas * c = new TCanvas(hstat->GetName(),hstat->GetName(),600,800);
  TLegend *legend = new TLegend(0.55,0.70,0.95,0.95, hstat->GetName());
  c->Divide(1,2);
  for (Int_t ifrac=0; ifrac<6; ifrac++){
    c->cd(1);
    if (ifrac==0) hisMean[ifrac]->Draw();
    legend->AddEntry(hisMean[ifrac],Form("%f",0.5+0.1*ifrac));
    hisMean[ifrac]->Draw("same");
    c->cd(2);
    if (ifrac==0) hisSigma[ifrac]->Draw();
    hisSigma[ifrac]->Draw("same");
  }
  c->Draw();
  legend->Draw();
  TString fname=hstat->GetName();
  fname.ReplaceAll("/","_");
  c->SaveAs(Form("pic/%s.eps",fname.Data()));
  c->SaveAs(Form("pic/%s.gif",fname.Data()));
  c->SaveAs(Form("pic/%s.root",fname.Data()));
}

