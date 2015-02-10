/// \file AliTPCclusterAnalysis.C
///
/// 0. Load libraries
/// ~~~
/// gSystem->Load("libSTAT");
/// .x ~/NimStyle.C
/// .L $ALICE_ROOT/TPC/fastSimul/AliTPCclusterFast.cxx+
/// .L $ALICE_ROOT/TPC/fastSimul/AliTPCclusterAnalysis.C
/// ~~~
///
/// 1. load tree
///
/// ~~~
/// LoadTree();
/// LoadTrack();
/// ~~~

TChain * treeCluster=0;
TChain * treeTrack=0;
AliTPCfastTrack * track =0;

void LoadTree(const char* fname="cluterSimul.root"){
  ///

  treeCluster = new TChain("simul","simul");
  treeCluster->AddFile(fname);
}

void LoadTrack(const char* fname="trackerSimul.root"){
  ///

  treeTrack = new TChain("simulTrack","simulTrack");
  treeTrack->AddFile(fname);
  TFile f(fname);
  
}






void MakeQNormalization(Int_t maxPoints){
  /// Normalize Q to the diffusion and angular effect

  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam;
  TMatrixD covMatrix;
  Int_t npoints;
  TString fstringQ="";
  fstringQ+="fDiff++";
  fstringQ+="abs(fAngleY)++";
  fstringQ+="abs(fAngleZ)++";
  TString *strQmax=0;
  TString *strQtot=0;
  //
  strQmax = toolkit.FitPlane(treeCluster,"s.GetQmax(0.33,2.5,1,1,1)/s.fNtot",fstringQ.Data(), "", chi2,npoints,fitParam,covMatrix,-1,0,50000);
  printf("Qmax norm\t%f\n",TMath::Sqrt(chi2/npoints));
  strQmax->Tokenize("++")->Print();
  treeCluster->SetAlias("qMaxCorr",strQmax->Data());
  //
  strQtot = toolkit.FitPlane(treeCluster,"s.GetQtot(0.33,2.5,1,1,1)/s.fNtot",fstringQ.Data(),"", chi2,npoints,fitParam,covMatrix,-1,0,50000);
  printf("Qmax norm\t%f\n",TMath::Sqrt(chi2/npoints));
  strQtot->Tokenize("++")->Print();
  treeCluster->SetAlias("qTotCorr",strQtot->Data());    
  //
  //
}




//
// Correction tests - normalization of response functions
//

void sum(){
  TF2 f2("f2","AliTPCclusterFast::GaussConvolution(x,y,2,0.5,0.5,0.5)",-3,3,-3,3);
  Float_t sumg=0;
  for (Float_t x=-5; x<5; x+=0.5)  for (Float_t y=-5; y<5; y+=0.5) sumg+=f2->Eval(x,y);
  printf("%f\n", sumg);
}



Double_t testSumGaus(Float_t k0,Float_t k1, Float_t s0, Float_t s1){
  TF2 f2("f2",Form("AliTPCclusterFast::GaussConvolution(x,y,%f,%f,%f,%f)",k0,k1,s0,s1),-2,2,-2,2);
  Float_t sumg=0;
  for (Float_t x=-5; x<5; x+=0.2)  for (Float_t y=-5; y<5; y+=0.2) sumg+=f2.Eval(x,y);
  sumg*=0.2*0.2;
  printf("%f\t%f\t%f\t%f\t%f\n",k0,k1,s0,s1, sumg);
  return sumg;
}

Double_t testSumExp(Float_t s0, Float_t k1){
  TF1 f1("f1",Form("AliTPCclusterFast::GaussExpConvolution(x,%f,%f)",s0,k1),-2,2);
  Float_t sumg=0;
  for (Float_t x=-5; x<5; x+=0.2)  sumg+=f1.Eval(x);
  sumg*=0.2;
  printf("%f\t%f\t%f\t%f\t%f\n",s0,k1, sumg);
  return sumg;
}




