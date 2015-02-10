/// \file CalibTime.C
///
/// Macro to visualize and analyze time dependent variables
///
/// ~~~
/// .L $ALICE_ROOT/TPC/CalibMacros/CalibTime.C
///
/// // 1. Init - Load libraries tools 
/// Init()
///
/// // 2. Load debug streamers 
/// InitChains()
/// ~~~

gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");

TFile f("CalibObjects.root");
AliTPCcalibTime *calibTime = (AliTPCcalibTime *)f->Get("TPCCalib")->FindObject("calibTime");

TPair * addPair= calibTime->GetMapDz()->FindObject(" D0SCO ");
THnSparse* addHist=dynamic_cast<THnSparseF*>(addPair->Value());

TGraph * gr = AliTPCcalibBase::FitSlices(addHist,2,0,100,100);
gr->SetMarkerColor(2);
gr->Draw("same*");
//Make Fit

AliSplineFit fit;
fit.SetGraph(gr)
fit->SetMinPoints(gr->GetN()+1);
fit->InitKnots(gr,2,0,0.001)
fit.SplineFit(0)
TGraph * grfit = fit.MakeGraph(gr->GetX()[0],gr->GetX()[gr->GetN()-1],50000,0);
gr->SetMarkerStyle(25);
//gr->Draw("alp");
grfit->SetLineColor(2);
grfit->Draw("lu");

//
// Chain Based analysis
//


TChain * chainLaser=0, *chainDz=0, *chaindEdx=0; 

void Init(){  
  /// Load neccesary libraries

  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  gROOT->Macro("~/NimStyle.C");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
}

void InitChains(){  
  AliXRDPROOFtoolkit tool;   
  //load debug streams
  chainLaser = tool.MakeChain("time.txt","laserInfo",0,10200);
  chainLaser->Lookup();
  // cosmic dZ - drfit velocity part 
  chainDz = tool.MakeChain("time.txt","cosmicDz",0,10200);
  chainDz->Lookup();
  //
  // cosmic dE - drift velocity part 
  chaindEdx = tool.MakeChain("time.txt","cosmicdEdx",0,10200);
  chaindEdx->Lookup();
  //
  // Set Alias
  //
  chaindEdx->SetAlias("side","(-1+(p0.fP[1]>0)*2)");
  chaindEdx->SetAlias("der","side*(dedx0-dedx1)/(dedx0+dedx1)");
  chaindEdx->SetAlias("derIO","side*(dedx0Out-dedx0In)/(dedx0In+dedx0Out)");
  chaindEdx->SetAlias("dr","(1-abs(p0.fP[1]/250))");
  chaindEdx->SetAlias("isOK","dedx0In>0&&dedx0Out>0&&dedx1In>0&&dedx1Out>0");
  chaindEdx->SetAlias("dedxM","(dedx0+dedx1)*0.5");

}


void MakeTglFitCosmic(){
  /// Fit the z correction factor

  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParamC,fitParamIO;
  TMatrixD covMatrix;  
  //
  // dedx delta
  // 
  TString fstring="";
  fstring+="side++";
  fstring+="p0.fP[3]++";
  //fstring+="p0.fP[3]*dr++";
  TString *strTheta = toolkit.FitPlane(chaindEdx,"der",fstring->Data(), "isOK", chi2,npoints,fitParamC,covMatrix);
  chaindEdx->SetAlias("derc",strTheta.Data());
  strTheta->Tokenize("+")->Print(); 
  //
  // dedx delta
  // 
  fstring="";
  fstring+="side++";
  fstring+="p0.fP[3]++";
  //fstring+="p0.fP[3]*dr++";
  TString *strTheta = toolkit.FitPlane(chaindEdx,"derIO",fstring->Data(), "isOK", chi2,npoints,fitParamIO,covMatrix);
  chaindEdx->SetAlias("derIOc",strTheta.Data());
  strTheta->Tokenize("+")->Print(); 
  //
  // Make Plot
  //
  chaindEdx->Draw("der:p0.fP[3]>>his(20,-0.5,0.5)","","prof");
  chaindEdx->Draw("der-derc:p0.fP[3]>>hisC(20,-0.5,0.5)","","prof");
  his->SetXTitle("tan(#theta)");
  his->SetYTitle("(Q_{u}-Q_{d})/(Q_{u}+Q_{d})");
  his->Draw("");
  hisC->Draw("same");
}


void MakeTglFitCosmic(){
  ///

  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParamC,fitParamIO;
  TMatrixD covMatrix;  
  TString fstring="";
  fstring="";
  fstring+="ppit++";
  fstring+="(trigger==1)++";
  fstring+="(trigger==2)++";
  fstring+="(trigger==4)++";
  fstring+="(trigger==8)++";
  TString *strPress = toolkit.FitPlane(chaindEdx,"(dedx0+dedx1)*0.5",fstring->Data(), "isOK", chi2,npoints,fitParamIO,covMatrix);
  //
  //
  chaindEdx->SetAlias("de",strTheta.Data());
  strTheta->Tokenize("+")->Print(); 


}
