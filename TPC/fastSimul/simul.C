/// \file simul.C
/// Macro to generate random tracks and clusters.
/// Fast MC - Geant equivalent used

Int_t  kmarkes[5]={20,21,24,25,23};
Int_t  kcolors[5]={1,2,4,3,6};

void simul(Int_t npoints, Double_t diffFactor){ 
  /// simulation submit script

  printf("Hallo world\n");
  gRandom->SetSeed(0);
  gROOT->LoadMacro("$mcPath/AliTPCclusterFast.cxx+");
 
  AliTPCclusterFast::fPRF = new TF1("fprf","gausn",-5,5);
  AliTPCclusterFast::fTRF = new TF1("ftrf","gausn",-5,5);
  AliTPCclusterFast::fPRF->SetParameters(1,0,0.5);
  AliTPCclusterFast::fTRF->SetParameters(1,0,0.5);
  //
  AliTPCtrackFast::Simul("trackerSimul.root",npoints, diffFactor); 
}



void Merge(){
  ///

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




void DrawdEdxResolExample(){
  /// Example analysis to make an space point resolution study

  TChain * chain  = AliXRDPROOFtoolkit::MakeChain("trackerSimul.list", "simulTrack",0,100); 
  chain->SetCacheSize(10000000000);

  //
  // 1.) Qmax/Qtot as function of the input ionization density
  //
  chain->Draw("tr.CookdEdxDmax(0,0.6,1,0,1,0)/tr.CookdEdxDtot(0,0.6,1,0,1,0):tr.fMNprim>>hisQtotMax(10,10,50)","","prof",10000);
  //
  // 2.) Non linearity due to the truncation Qtot_{60%}/Qtot 100% 
  //
  chain->Draw("tr.CookdEdxDtot(0,0.6,1,0,1,0)/tr.CookdEdxDtot(0,0.99,1,0,1,0):tr.fMNprim>>hisQtot60100(10,10,50)","","prof",10000);
  //
  // 3.) 
  //
  chain->Draw("tr.CookdEdxDtot(0,0.6,1,0,1,0)/tr.fMNprim:tr.fMNprim>>profQtot60(10,10,50)","","prof",10000);
  chain->Draw("tr.CookdEdxDtot(0,0.99,1,0,1,0)/tr.fMNprim:tr.fMNprim>>profQtot100(10,10,50)","","profsame",10000);

}
