//------------------------------------------------------------------------------
// Implementation of AliComparisonRes class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliComparisonRes.
//
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  //AliComparisonRes * compObj = (AliComparisonRes*)f.Get("AliComparisonRes");
  AliComparisonRes * compObj = (AliComparisonRes*)cOutput->FindObject("AliComparisonRes");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderRes" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Res.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include <iostream>

#include "TFile.h"
#include "TCint.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TGraph.h"

#include "AliESDEvent.h"   
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDVertex.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
#include "AliTracker.h" 

#include "AliMathBase.h"
#include "AliTreeDraw.h" 

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonRes.h" 

using namespace std;

ClassImp(AliComparisonRes)

//_____________________________________________________________________________
AliComparisonRes::AliComparisonRes():
  AliComparisonObject("AliComparisonRes"),

  // Vertices
  fMCVertex(0),  //-> MC primary vertex
  fRecVertex(0), //-> Reconstructed primary vertex

  // global observables
  fPhiTanPtTPC(0),
  fPhiTanPtTPCITS(0),

  // Resolution 

  // TPC only
  fPtResolTPC(0),
  fPtPullTPC(0),
  fPhiResolTanTPC(0),
  fTanResolTanTPC(0),
  fPhiPullTanTPC(0),
  fTanPullTanTPC(0),

  // TPC+ITS
  fPtResolTPCITS(0),
  fPtPullTPCITS(0),
  fPhiResolTanTPCITS(0),
  fTanResolTanTPCITS(0),
  fPhiPullTanTPCITS(0),
  fTanPullTanTPCITS(0),

  //
  // Resolution constrained param
  //
  fCPhiResolTan(0),  // angular resolution -  constrained
  fCTanResolTan(0),  // angular resolution -  constrained
  fCPtResolTan(0),   // pt resolution      -  constrained
  fCPhiPullTan(0),   // angular resolution -  constrained
  fCTanPullTan(0),   // angular resolution -  constrained
  fCPtPullTan(0),    // pt resolution      -  constrained

  //
  // Parametrisation histograms
  //

  f1Pt2ResolS1PtTPC(0),
  f1Pt2ResolS1PtTPCITS(0),
  fYResolS1PtTPC(0),
  fYResolS1PtTPCITS(0),
  fZResolS1PtTPC(0),
  fZResolS1PtTPCITS(0),
  fPhiResolS1PtTPC(0),
  fPhiResolS1PtTPCITS(0),
  fThetaResolS1PtTPC(0),
  fThetaResolS1PtTPCITS(0),

  // constrained
  fC1Pt2ResolS1PtTPC(0),
  fC1Pt2ResolS1PtTPCITS(0),
  fCYResolS1PtTPC(0),
  fCYResolS1PtTPCITS(0),
  fCZResolS1PtTPC(0),
  fCZResolS1PtTPCITS(0),
  fCPhiResolS1PtTPC(0),
  fCPhiResolS1PtTPCITS(0),
  fCThetaResolS1PtTPC(0),
  fCThetaResolS1PtTPCITS(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0)
{
  Init();
}

//_____________________________________________________________________________
AliComparisonRes::~AliComparisonRes(){
  
  // Vertices
  if(fMCVertex) delete fMCVertex; fMCVertex=0;     
  if(fRecVertex) delete fRecVertex; fRecVertex=0;     

  // Global observables
  //
  if(fPhiTanPtTPC) delete fPhiTanPtTPC; fPhiTanPtTPC=0;     
  if(fPhiTanPtTPCITS) delete fPhiTanPtTPCITS; fPhiTanPtTPCITS=0;     

  // Resolution histograms
  if(fPtResolTPC) delete fPtResolTPC; fPtResolTPC=0;     
  if(fPtPullTPC) delete fPtPullTPC; fPtPullTPC=0;     
  if(fPhiResolTanTPC) delete fPhiResolTanTPC; fPhiResolTanTPC=0;     
  if(fTanResolTanTPC) delete fTanResolTanTPC; fTanResolTanTPC=0;     
  if(fPhiPullTanTPC) delete fPhiPullTanTPC; fPhiPullTanTPC=0;     
  if(fTanPullTanTPC) delete fTanPullTanTPC; fTanPullTanTPC=0;     

  if(fPtResolTPCITS) delete fPtResolTPCITS; fPtResolTPCITS=0;     
  if(fPtPullTPCITS) delete fPtPullTPCITS; fPtPullTPCITS=0;     
  if(fPhiResolTanTPCITS) delete fPhiResolTanTPCITS; fPhiResolTanTPCITS=0;     
  if(fTanResolTanTPCITS) delete fTanResolTanTPCITS; fTanResolTanTPCITS=0;     
  if(fPhiPullTanTPCITS) delete fPhiPullTanTPCITS; fPhiPullTanTPCITS=0;     
  if(fTanPullTanTPCITS) delete fTanPullTanTPCITS; fTanPullTanTPCITS=0;     
  
  // Resolution histograms (constrained param)
  if(fCPhiResolTan) delete fCPhiResolTan; fCPhiResolTan=0;
  if(fCTanResolTan) delete fCTanResolTan; fCTanResolTan=0;
  if(fCPtResolTan)  delete fCPtResolTan;  fCPtResolTan=0; 
  if(fCPhiPullTan)  delete fCPhiPullTan;  fCPhiPullTan=0;
  if(fCTanPullTan)  delete fCTanPullTan;  fCTanPullTan=0;
  if(fCPtPullTan)   delete fCPtPullTan;   fCPtPullTan=0;

  // Parametrisation histograms
  // 
  if(f1Pt2ResolS1PtTPC) delete f1Pt2ResolS1PtTPC; f1Pt2ResolS1PtTPC=0;
  if(f1Pt2ResolS1PtTPCITS) delete f1Pt2ResolS1PtTPCITS; f1Pt2ResolS1PtTPCITS=0;
  if(fYResolS1PtTPC) delete fYResolS1PtTPC; fYResolS1PtTPC=0;
  if(fYResolS1PtTPCITS) delete fYResolS1PtTPCITS; fYResolS1PtTPCITS=0;
  if(fZResolS1PtTPC) delete fZResolS1PtTPC; fZResolS1PtTPC=0;
  if(fZResolS1PtTPCITS) delete fZResolS1PtTPCITS; fZResolS1PtTPCITS=0;
  if(fPhiResolS1PtTPC) delete fPhiResolS1PtTPC; fPhiResolS1PtTPC=0;
  if(fPhiResolS1PtTPCITS) delete fPhiResolS1PtTPCITS; fPhiResolS1PtTPCITS=0;
  if(fThetaResolS1PtTPC) delete fThetaResolS1PtTPC; fThetaResolS1PtTPC=0;
  if(fThetaResolS1PtTPCITS) delete fThetaResolS1PtTPCITS; fThetaResolS1PtTPCITS=0;

  // constrained
  if(fC1Pt2ResolS1PtTPC) delete fC1Pt2ResolS1PtTPC; fC1Pt2ResolS1PtTPC=0;
  if(fC1Pt2ResolS1PtTPCITS) delete fC1Pt2ResolS1PtTPCITS; fC1Pt2ResolS1PtTPCITS=0;
  if(fCYResolS1PtTPC) delete fCYResolS1PtTPC; fCYResolS1PtTPC=0;
  if(fCYResolS1PtTPCITS) delete fCYResolS1PtTPCITS; fCYResolS1PtTPCITS=0;
  if(fCZResolS1PtTPC) delete fCZResolS1PtTPC; fCZResolS1PtTPC=0;
  if(fCZResolS1PtTPCITS) delete fCZResolS1PtTPCITS; fCZResolS1PtTPCITS=0;
  if(fCPhiResolS1PtTPC) delete fCPhiResolS1PtTPC; fCPhiResolS1PtTPC=0;
  if(fCPhiResolS1PtTPCITS) delete fCPhiResolS1PtTPCITS; fCPhiResolS1PtTPCITS=0;
  if(fCThetaResolS1PtTPC) delete fCThetaResolS1PtTPC; fCThetaResolS1PtTPC=0;
  if(fCThetaResolS1PtTPCITS) delete fCThetaResolS1PtTPCITS; fCThetaResolS1PtTPCITS=0;

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliComparisonRes::Init(){

  // Control histograms
  fMCVertex = new TH3F("MCVertex","MC vertex Xv:Yv:Zv",100,-0.05,0.05,100,-0.05,0.05,100,-10,10);
  fMCVertex->SetXTitle("Xv [cm]");
  fMCVertex->SetYTitle("Yv [cm]");
  fMCVertex->SetZTitle("Zv [cm]");

  fRecVertex = new TH3F("RecVertex","Rec prim. vertex Xv:Yv:Zv",100,-0.01,0.01,100,-0.01,0.01,100,-0.5,0.5);
  fRecVertex->SetXTitle("Xv [cm]");
  fRecVertex->SetYTitle("Yv [cm]");
  fRecVertex->SetZTitle("Zv [cm]");

  //
  fPhiTanPtTPC = new TH3F("PhiTanPtTPC","phi vs tan#theta vs pt - TPC only",200,-3.5,3.5,50,-2,2,100,0.1,3.0);
  fPhiTanPtTPC->SetXTitle("#phi [rad]");
  fPhiTanPtTPC->SetYTitle("tan#theta");
  fPhiTanPtTPC->SetZTitle("p_{t} [GeV]");

  fPhiTanPtTPCITS = new TH3F("PhiTanPtTPCITS","phi vs tan#theta vs pt - TPC+ITS",200,-3.5,3.5,50,-2,2,100,0.1,3.0);
  fPhiTanPtTPCITS->SetXTitle("#phi [rad]");
  fPhiTanPtTPCITS->SetYTitle("tan#theta");
  fPhiTanPtTPCITS->SetZTitle("p_{t} [GeV]");

  //
  fPtResolTPC = new TH2F("PtResolTPC","pt resol",10, 0.1,3,200,-0.2,0.2);
  fPtResolTPC->SetXTitle("p_{t} [GeV]");
  fPtResolTPC->SetYTitle("#Deltap_{t}/p_{t}");

  fPtPullTPC = new TH2F("fPtPullTPC","pt pull",10, 0.1,3,200,-6,6);
  fPtPullTPC->SetXTitle("p_{t} [GeV]");
  fPtPullTPC->SetYTitle("#Deltap_{t}/#Sigma");

  fPhiResolTanTPC = new TH2F("PhiResolTanTPC","PhiResolTanTPC",50, -2,2,200,-0.025,0.025);   
  fPhiResolTanTPC->SetXTitle("tan(#theta)");
  fPhiResolTanTPC->SetYTitle("#Delta#phi [rad]");

  fTanResolTanTPC = new TH2F("TanResolTanTPC","TanResolTanTPC",50, -2,2,200,-0.025,0.025);
  fTanResolTanTPC->SetXTitle("tan(#theta)");
  fTanResolTanTPC->SetYTitle("#Delta#theta [rad]");

  fPhiPullTanTPC = new TH2F("PhiPullTanTPC","PhiPullTanTPC",50, -2,2,200,-5,5);   
  fPhiPullTanTPC->SetXTitle("Tan(#theta)");
  fPhiPullTanTPC->SetYTitle("#Delta#phi/#Sigma");

  fTanPullTanTPC = new TH2F("TanPullTanTPC","TanPullTanTPC",50, -2,2,200,-5,5);
  fTanPullTanTPC->SetXTitle("Tan(#theta)");
  fTanPullTanTPC->SetYTitle("#Delta#theta/#Sigma");

  fPtResolTPCITS = new TH2F("PtResolTPCITS","pt resol",10, 0.1,3,200,-0.2,0.2);
  fPtResolTPCITS->SetXTitle("p_{t} [GeV]");
  fPtResolTPCITS->SetYTitle("#Deltap_{t}/p_{t}");

  fPtPullTPCITS = new TH2F("fPtPullTPCITS","pt pull",10, 0.1,3,200,-6,6);
  fPtPullTPCITS->SetXTitle("p_{t} [GeV]");
  fPtPullTPCITS->SetYTitle("#Deltap_{t}/#Sigma");

  fPhiResolTanTPCITS = new TH2F("PhiResolTanTPCITS","PhiResolTanTPCITS",50, -2,2,200,-0.025,0.025);   
  fPhiResolTanTPCITS->SetXTitle("tan(#theta)");
  fPhiResolTanTPCITS->SetYTitle("#Delta#phi [rad]");

  fTanResolTanTPCITS = new TH2F("TanResolTanTPCITS","TanResolTanTPCITS",50, -2,2,200,-0.025,0.025);
  fTanResolTanTPCITS->SetXTitle("tan(#theta)");
  fTanResolTanTPCITS->SetYTitle("#Delta#theta [rad]");

  fPhiPullTanTPCITS = new TH2F("PhiPullTanTPCITS","PhiPullTanTPCITS",50, -2,2,200,-5,5);   
  fPhiPullTanTPCITS->SetXTitle("Tan(#theta)");
  fPhiPullTanTPCITS->SetYTitle("#Delta#phi/#Sigma");

  fTanPullTanTPCITS = new TH2F("TanPullTanTPCITS","TanPullTanTPCITS",50, -2,2,200,-5,5);
  fTanPullTanTPCITS->SetXTitle("Tan(#theta)");
  fTanPullTanTPCITS->SetYTitle("#Delta#theta/#Sigma");

  // Resolution constrained param
  fCPhiResolTan = new TH2F("CPhiResolTan","CPhiResolTan",50, -2,2,200,-0.025,0.025);   
  fCPhiResolTan->SetXTitle("tan(#theta)");
  fCPhiResolTan->SetYTitle("#Delta#phi [rad]");

  fCTanResolTan = new TH2F("CTanResolTan","CTanResolTan",50, -2,2,200,-0.025,0.025);
  fCTanResolTan->SetXTitle("tan(#theta)");
  fCTanResolTan->SetYTitle("#Delta#theta [rad]");

  fCPtResolTan=new TH2F("CPtResol","CPtResol",50, -2,2,200,-0.2,0.2);    
  fCPtResolTan->SetXTitle("Tan(#theta)");
  fCPtResolTan->SetYTitle("#Deltap_{t}/p_{t}");

  fCPhiPullTan = new TH2F("CPhiPullTan","CPhiPullTan",50, -2,2,200,-5,5);   
  fCPhiPullTan->SetXTitle("Tan(#theta)");
  fCPhiPullTan->SetYTitle("#Delta#phi/#Sigma");

  fCTanPullTan = new TH2F("CTanPullTan","CTanPullTan",50, -2,2,200,-5,5);
  fCTanPullTan->SetXTitle("Tan(#theta)");
  fCTanPullTan->SetYTitle("#Delta#theta/#Sigma");

  fCPtPullTan=new TH2F("CPtPull","CPtPull",50, -2,2,200,-5,5);    
  fCPtPullTan->SetXTitle("Tan(#theta)");
  fCPtPullTan->SetYTitle("(1/mcp_{t}-1/p_{t})/#Sigma");
  
  //
  // Parametrisation histograms
  // 

  f1Pt2ResolS1PtTPC = new TH2F("f1Pt2ResolS1PtTPC","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs sqrt(1/pt))",100,0,3,200,-0.010,0.010);  
  f1Pt2ResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  f1Pt2ResolS1PtTPC->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  f1Pt2ResolS1PtTPCITS = new TH2F("f1Pt2ResolS1PtTPCITS","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs sqrt(1/pt))",100,0,3,200,-0.010,0.010);  
  f1Pt2ResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  f1Pt2ResolS1PtTPCITS->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  fYResolS1PtTPC = new TH2F("fYResolS1PtTPC","fYResolS1PtTPC",100, 0,3,200,-1.0,1.0);   
  fYResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fYResolS1PtTPC->SetYTitle("#DeltaY");

  fYResolS1PtTPCITS = new TH2F("fYResolS1PtTPCITS","fYResolS1PtTPCITS",100, 0,3,200,-0.05,0.05);   
  fYResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fYResolS1PtTPCITS->SetYTitle("#DeltaY");

  fZResolS1PtTPC = new TH2F("fZResolS1PtTPC","fZResolS1PtTPC",100, 0,3,200,-1.0,1.0);   
  fZResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fZResolS1PtTPC->SetYTitle("#DeltaZ");

  fZResolS1PtTPCITS = new TH2F("fZResolS1PtTPCITS","fZResolS1PtTPCITS",100, 0,3,200,-0.05,0.05);   
  fZResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fZResolS1PtTPCITS->SetYTitle("#DeltaZ");

  fPhiResolS1PtTPC = new TH2F("fPhiResolS1PtTPC","fPhiResolS1PtTPC",100, 0,3,200,-0.025,0.025);   
  fPhiResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fPhiResolS1PtTPC->SetYTitle("#Delta#phi");

  fPhiResolS1PtTPCITS = new TH2F("fPhiResolS1PtTPCITS","fPhiResolS1PtTPCITS",100, 0,3,200,-0.01,0.01);   
  fPhiResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fPhiResolS1PtTPCITS->SetYTitle("#Delta#phi");

  fThetaResolS1PtTPC = new TH2F("fThetaResolS1PtTPC","fThetaResolS1PtTPC",100, 0,3,200,-0.025,0.025);   
  fThetaResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fThetaResolS1PtTPC->SetYTitle("#Delta#theta");

  fThetaResolS1PtTPCITS = new TH2F("fThetaResolS1PtTPCITS","fThetaResolS1PtTPCITS",100, 0,3,200,-0.01,0.01);   
  fThetaResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fThetaResolS1PtTPCITS->SetYTitle("#Delta#theta");
  
  // constrained
  fC1Pt2ResolS1PtTPC = new TH2F("fC1Pt2ResolS1PtTPC","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs 1/pt)",100,0,3,200,-0.010,0.010);  
  fC1Pt2ResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fC1Pt2ResolS1PtTPC->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  fC1Pt2ResolS1PtTPCITS = new TH2F("fC1Pt2ResolS1PtTPCITS","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs 1/pt)",100,0,3,200,-0.010,0.010);  
  fC1Pt2ResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fC1Pt2ResolS1PtTPCITS->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  fCYResolS1PtTPC = new TH2F("fCYResolS1PtTPC","fCYResolS1PtTPC",100, 0,3,200,-1.0,1.0);   
  fCYResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fCYResolS1PtTPC->SetYTitle("#DeltaY");

  fCYResolS1PtTPCITS = new TH2F("fCYResolS1PtTPCITS","fCYResolS1PtTPCITS",100, 0,3,200,-0.01,0.01);   
  fCYResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fCYResolS1PtTPCITS->SetYTitle("#DeltaY");

  fCZResolS1PtTPC = new TH2F("fCZResolS1PtTPC","fCZResolS1PtTPC",100, 0,3,200,-1.0,1.0);   
  fCZResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fCZResolS1PtTPC->SetYTitle("#DeltaZ");

  fCZResolS1PtTPCITS = new TH2F("fCZResolS1PtTPCITS","fCZResolS1PtTPCITS",100, 0,3,200,-0.025,0.025);   
  fCZResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fCZResolS1PtTPCITS->SetYTitle("#DeltaZ");

  fCPhiResolS1PtTPC = new TH2F("fCPhiResolS1PtTPC","fCPhiResolS1PtTPC",100, 0,3,200,-0.025,0.025);   
  fCPhiResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fCPhiResolS1PtTPC->SetYTitle("#Delta#phi");

  fCPhiResolS1PtTPCITS = new TH2F("fCPhiResolS1PtTPCITS","fCPhiResolS1PtTPCITS",100, 0,3,200,-0.003,0.003);   
  fCPhiResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fCPhiResolS1PtTPCITS->SetYTitle("#Delta#phi");

  fCThetaResolS1PtTPC = new TH2F("fCThetaResolS1PtTPC","fCThetaResolS1PtTPC",100, 0,3,200,-0.025,0.025);   
  fCThetaResolS1PtTPC->SetXTitle("#sqrt{1/mcp_{t}}");
  fCThetaResolS1PtTPC->SetYTitle("#Delta#theta");

  fCThetaResolS1PtTPCITS = new TH2F("fCThetaResolS1PtTPCITS","fCThetaResolS1PtTPCITS",100, 0,3,200,-0.005,0.005);   
  fCThetaResolS1PtTPCITS->SetXTitle("#sqrt{1/mcp_{t}}");
  fCThetaResolS1PtTPCITS->SetYTitle("#Delta#theta");

  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderRes","Analysis Resolution Folder");
}

//_____________________________________________________________________________
void AliComparisonRes::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill resolution comparison information 
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  AliESDVertex vertexMC;  // MC primary vertex

  Int_t clusterITS[200];
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t deltaPt, pullPt, deltaPhi,pullPhi, deltaTan, pullTan, delta1Pt2, deltaY1Pt, deltaZ1Pt, deltaPhi1Pt, deltaTheta1Pt; 
  Float_t deltaPtTPC, pullPtTPC, deltaPhiTPC, pullPhiTPC, deltaTanTPC, pullTanTPC, delta1Pt2TPC, deltaY1PtTPC, deltaZ1PtTPC, deltaPhi1PtTPC, deltaTheta1PtTPC; 

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t s1mcpt = TMath::Sqrt(1./infoMC->GetParticle().Pt());
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Float_t eta =  infoMC->GetParticle().Eta();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return; // TPC refit
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;

  // calculate and set prim. vertex
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // Fill MC vertex histo
  fMCVertex->Fill(vertexMC.GetXv(),vertexMC.GetYv(),vertexMC.GetZv());

   
  // Fill Rec vertex histo
  //fRecVertex->Fill(infoRC->GetESDvertex()->GetXv(),infoRC->GetESDvertex()->GetYv(),infoRC->GetESDvertex()->GetZv());
  //printf("Rec vertex xv %f, yv %f, zv %f \n",infoRC->GetESDvertex()->GetXv(),infoRC->GetESDvertex()->GetYv(),infoRC->GetESDvertex()->GetZv());

  deltaPt= (mcpt-infoRC->GetESDtrack()->Pt())/mcpt;  
  pullPt= (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/TMath::Sqrt(infoRC->GetESDtrack()->GetSigma1Pt2());  
  deltaPhi = TMath::ATan2(infoRC->GetESDtrack()->Py(),infoRC->GetESDtrack()->Px())-
                     TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
  pullPhi = deltaPhi/TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaSnp2()); 

  deltaTan = TMath::ATan2(infoRC->GetESDtrack()->Pz(),infoRC->GetESDtrack()->Pt())-
                     TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
  pullTan = deltaTan/TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaTgl2());

  delta1Pt2 = (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/TMath::Power(1+1/mcpt,2);       
  deltaY1Pt = (infoMC->GetParticle().Vy()-infoRC->GetESDtrack()->GetY()) / (0.2+1/mcpt);
  deltaZ1Pt = (infoMC->GetParticle().Vz()-infoRC->GetESDtrack()->GetZ()) / (0.2+1/mcpt);
  deltaPhi1Pt = deltaPhi   / (0.1+1/mcpt);
  deltaTheta1Pt = deltaTan / (0.1+1/mcpt);

  // calculate track parameters at vertex
  const AliExternalTrackParam *innerTPC =  0;
  if ((innerTPC = infoRC->GetESDtrack()->GetTPCInnerParam()) != 0)
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);

      // Fill parametrisation histograms (only TPC track)
      if(bDCAStatus) 
	  {
		  
			deltaPtTPC= (mcpt-innerTPC->Pt())/mcpt;  
			pullPtTPC= (1/mcpt-innerTPC->OneOverPt())/TMath::Sqrt(innerTPC->GetSigma1Pt2());  
			deltaPhiTPC = TMath::ATan2(innerTPC->Py(),innerTPC->Px())-
								TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
                        pullPhiTPC = deltaPhiTPC/TMath::Sqrt(innerTPC->GetSigmaSnp2()); 

			deltaTanTPC = TMath::ATan2(innerTPC->Pz(),innerTPC->Pt())-
								TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
                        pullTanTPC = deltaTanTPC/TMath::Sqrt(infoRC->GetESDtrack()->GetSigmaTgl2());

			delta1Pt2TPC = (1/mcpt-innerTPC->OneOverPt())/TMath::Power(1+1/mcpt,2);       
			deltaY1PtTPC= (infoMC->GetParticle().Vy()-innerTPC->GetY()) / (0.2+1/mcpt);
			deltaZ1PtTPC = (infoMC->GetParticle().Vz()-innerTPC->GetZ()) / (0.2+1/mcpt);
			deltaPhi1PtTPC = deltaPhiTPC   / (0.1+1/mcpt);
			deltaTheta1PtTPC = deltaTanTPC / (0.1+1/mcpt);

	        fPhiTanPtTPC->Fill(TMath::ATan2(innerTPC->Py(),innerTPC->Px()), 
			                   TMath::ATan2(innerTPC->Pz(),innerTPC->Pt()),
							   innerTPC->Pt()); 

			f1Pt2ResolS1PtTPC->Fill(s1mcpt,delta1Pt2TPC);
			fYResolS1PtTPC->Fill(s1mcpt,deltaY1PtTPC);
			fZResolS1PtTPC->Fill(s1mcpt,deltaZ1PtTPC);
			fPhiResolS1PtTPC->Fill(s1mcpt,deltaPhi1PtTPC);
			fThetaResolS1PtTPC->Fill(s1mcpt,deltaTheta1PtTPC);

	        fPtResolTPC->Fill(mcpt,deltaPtTPC);
		    fPtPullTPC->Fill(mcpt,pullPtTPC);
		    fPhiResolTanTPC->Fill(tantheta,deltaPhiTPC);
		    fPhiPullTanTPC->Fill(tantheta,pullPhiTPC);
		    fTanResolTanTPC->Fill(tantheta,deltaTanTPC);
		    fTanPullTanTPC->Fill(tantheta,pullTanTPC);
	  }
	  delete track;
    }
  }

  // TPC and ITS clusters in the system
  if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS()) 
  {
	  fPhiTanPtTPCITS->Fill(TMath::ATan2(infoRC->GetESDtrack()->Py(),infoRC->GetESDtrack()->Px()), 
	                        TMath::ATan2(infoRC->GetESDtrack()->Pz(),infoRC->GetESDtrack()->Pt()), 
	                        infoRC->GetESDtrack()->Pt()); 

      f1Pt2ResolS1PtTPCITS->Fill(s1mcpt,delta1Pt2);
      fYResolS1PtTPCITS->Fill(s1mcpt,deltaY1Pt);
      fZResolS1PtTPCITS->Fill(s1mcpt,deltaZ1Pt);
      fPhiResolS1PtTPCITS->Fill(s1mcpt,deltaPhi1Pt);
      fThetaResolS1PtTPCITS->Fill(s1mcpt,deltaTheta1Pt);

	  fPtResolTPCITS->Fill(mcpt,deltaPt);
	  fPtPullTPCITS->Fill(mcpt,pullPt);

	  fPhiResolTanTPCITS->Fill(tantheta,deltaPhi);
	  fPhiPullTanTPCITS->Fill(tantheta,pullPhi);
	  fTanResolTanTPCITS->Fill(tantheta,deltaTan);
	  fTanPullTanTPCITS->Fill(tantheta,pullTan);
  }
}

//_____________________________________________________________________________
void AliComparisonRes::ProcessConstrained(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill resolution comparison information (constarained parameters) 
  //
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  AliESDVertex vertexMC;  // MC primary vertex

  Int_t clusterITS[200];
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
 
  Float_t deltaPt, pullPt, deltaPhi, pullPhi, deltaTan, pullTan, delta1Pt2, deltaY1Pt, deltaZ1Pt, deltaPhi1Pt, deltaTheta1Pt; 
  Float_t deltaPtTPC, pullPtTPC, deltaPhiTPC, deltaTanTPC, delta1Pt2TPC, deltaY1PtTPC, deltaZ1PtTPC, deltaPhi1PtTPC, deltaTheta1PtTPC; 

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t s1mcpt = TMath::Sqrt(1./infoMC->GetParticle().Pt());
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();
  
  // Check selection cuts
 
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;
  if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;

// calculate and set prim. vertex
  vertexMC.SetXv( infoMC->GetParticle().Vx() - dv[0] );
  vertexMC.SetYv( infoMC->GetParticle().Vy() - dv[1] );
  vertexMC.SetZv( infoMC->GetParticle().Vz() - dv[2] );

  // constrained parameters resolution
  const AliExternalTrackParam * cparam = infoRC->GetESDtrack()->GetConstrainedParam();
  deltaPt= (mcpt-cparam->Pt())/mcpt;  
  pullPt= (1/mcpt-cparam->OneOverPt())/TMath::Sqrt(cparam->GetSigma1Pt2());          
  deltaPhi = TMath::ATan2(cparam->Py(),cparam->Px())-
                     TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
  pullPhi = deltaPhi/TMath::Sqrt(cparam->GetSigmaSnp2()); 
  deltaTan = TMath::ATan2(cparam->Pz(),cparam->Pt())-TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
  pullTan = deltaPhi/TMath::Sqrt(cparam->GetSigmaSnp2()); 


  delta1Pt2 = (1/mcpt-cparam->OneOverPt())/TMath::Power(1+1/mcpt,2);       

  deltaY1Pt = (infoMC->GetParticle().Vy()-cparam->GetY()) / (0.2+1/mcpt);
  deltaZ1Pt = (infoMC->GetParticle().Vz()-cparam->GetZ()) / (0.2+1/mcpt);
  deltaPhi1Pt = deltaPhi   / (0.1+1/mcpt);
  deltaTheta1Pt = deltaTan / (0.1+1/mcpt);

  // calculate track parameters at vertex
  const AliExternalTrackParam *innerTPC =  0;
  if ((innerTPC = infoRC->GetESDtrack()->GetTPCInnerParam()) != 0)
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      Bool_t bDCAStatus = track->PropagateToDCA(&vertexMC,field,kMaxD,dca,cov);

      // Fill parametrisation histograms (only TPC track)
      if(bDCAStatus) 
	  {
		  deltaPtTPC= (mcpt-innerTPC->Pt())/mcpt;  
		  pullPtTPC= (1/mcpt-innerTPC->OneOverPt())/TMath::Sqrt(innerTPC->GetSigma1Pt2());  
		  deltaPhiTPC = TMath::ATan2(innerTPC->Py(),innerTPC->Px())-
								TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());

		  deltaTanTPC = TMath::ATan2(innerTPC->Pz(),innerTPC->Pt())-
								TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());

		  delta1Pt2TPC = (1/mcpt-innerTPC->OneOverPt())/TMath::Power(1+1/mcpt,2);       
		  deltaY1PtTPC= (infoMC->GetParticle().Vy()-innerTPC->GetY()) / (0.2+1/mcpt);
		  deltaZ1PtTPC = (infoMC->GetParticle().Vz()-innerTPC->GetZ()) / (0.2+1/mcpt);
		  deltaPhi1PtTPC = deltaPhiTPC   / (0.1+1/mcpt);
		  deltaTheta1PtTPC = deltaTanTPC / (0.1+1/mcpt);

          fC1Pt2ResolS1PtTPC->Fill(s1mcpt,delta1Pt2TPC);
          fCYResolS1PtTPC->Fill(s1mcpt,deltaY1PtTPC);
          fCZResolS1PtTPC->Fill(s1mcpt,deltaZ1PtTPC);
          fCPhiResolS1PtTPC->Fill(s1mcpt,deltaPhi1PtTPC);
          fCThetaResolS1PtTPC->Fill(s1mcpt,deltaTheta1PtTPC);
	  }
	  delete track;
    }
  }

 // TPC and ITS in the system
  if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS()) 
  {
      fC1Pt2ResolS1PtTPCITS->Fill(s1mcpt,delta1Pt2);
      fCYResolS1PtTPCITS->Fill(s1mcpt,deltaY1Pt);
      fCZResolS1PtTPCITS->Fill(s1mcpt,deltaZ1Pt);
      fCPhiResolS1PtTPCITS->Fill(s1mcpt,deltaPhi1Pt);
      fCThetaResolS1PtTPCITS->Fill(s1mcpt,deltaTheta1Pt);

	  // Fill histograms
	  fCPtResolTan->Fill(tantheta,deltaPt);
	  fCPtPullTan->Fill(tantheta,pullPt);
	  fCPhiResolTan->Fill(tantheta,deltaPhi);
	  fCPhiPullTan->Fill(tantheta,pullPhi);
	  fCTanResolTan->Fill(tantheta,deltaTan);
	  fCTanPullTan->Fill(tantheta,pullTan);
  }
}
 
//_____________________________________________________________________________
void AliComparisonRes::Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  
  // Process comparison information 
  Process(infoMC,infoRC);
  ProcessConstrained(infoMC,infoRC);
}

//_____________________________________________________________________________
TH1F* AliComparisonRes::MakeResol(TH2F * his, Int_t integ, Bool_t type){
  // Create resolution histograms
 
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoI(his,&hism,integ);
  if (type) return hism;
  else 
    return hisr;
}

//_____________________________________________________________________________
void AliComparisonRes::Analyse(){
  // Analyse comparison information and store output histograms
  // in the folder "folderRes"
  //
 
  TH1::AddDirectory(kFALSE);

  AliComparisonRes * comp=this;
  TH1F *hiss=0;
  TObjArray *aFolderObj = new TObjArray;

  // write results in the folder 

  TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  c->cd();
  //
  hiss = comp->MakeResol(comp->fPtResolTPC,1,0);
  hiss->SetXTitle("p_{t} (GeV)");
  hiss->SetYTitle("#sigmap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->SetName("PtResolPtTPC");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fPtResolTPC,1,1);
  hiss->SetXTitle("p_{t} (GeV)");
  hiss->SetYTitle("mean #Deltap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->SetName("PtMeanPtTPC");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fPhiResolTanTPC,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#phi (rad)");
  hiss->Draw();
  hiss->SetName("PhiResolTanTPC");

  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fPhiResolTanTPC,1,1);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("mean #Delta#phi (rad)");
  hiss->Draw();
  hiss->SetName("PhiMeanTanTPC");
  

  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fTanResolTanTPC,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#theta (rad)");
  hiss->Draw();
  hiss->SetName("ThetaResolTanTPC");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fTanResolTanTPC,1,1);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("mean #Delta#theta (rad)");
  hiss->Draw();
  hiss->SetName("ThetaMeanTanTPC");
  
  aFolderObj->Add(hiss);

  //
  //
  hiss = comp->MakeResol(comp->fPtResolTPCITS,1,0);
  hiss->SetXTitle("p_{t}");
  hiss->SetYTitle("#sigmap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->SetName("PtResolPtTPCITS");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fPtResolTPCITS,1,1);
  hiss->SetXTitle("p_{t}");
  hiss->SetYTitle("mean #Deltap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->SetName("PtMeanPtTPCITS");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fPhiResolTanTPCITS,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#phi (rad)");
  hiss->Draw();
  hiss->SetName("PhiResolTanTPCITS");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fPhiResolTanTPCITS,1,1);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("mean #Delta#phi (rad)");
  hiss->Draw();
  hiss->SetName("PhiMeanTanTPCITS");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fTanResolTanTPCITS,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#theta (rad)");
  hiss->Draw();
  hiss->SetName("ThetaResolTanTPCITS");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->fTanResolTanTPCITS,1,1);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("mean #Delta#theta (rad)");
  hiss->Draw();
  hiss->SetName("ThetaMeanTanTPCITS");
  
  aFolderObj->Add(hiss);

  //
  //
  hiss = comp->MakeResol(comp->fCPtResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigmap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->SetName("CptResolTan");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fCPhiResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#phi (rad)");
  hiss->Draw();
  hiss->SetName("CPhiResolTan");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fCTanResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#theta (rad)");
  hiss->Draw();
  hiss->SetName("CThetaResolTan");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fCPtPullTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/#Sigma(1/p_{t})");
  hiss->Draw();
  hiss->SetName("CptPullTan");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fC1Pt2ResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->SetName("C1Pt2ResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fC1Pt2ResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->SetName("C1Pt2ResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fCYResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CYResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fCYResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CYResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fCZResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CZResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fCZResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CZResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fCPhiResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CPhiResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fCPhiResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CPhiResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fCThetaResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CThetaResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fCThetaResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("CThetaResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);

  //
  hiss = comp->MakeResol(comp->f1Pt2ResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->SetName("OnePt2ResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->f1Pt2ResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->SetName("OnePt2ResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fYResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("YResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fYResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("YResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fZResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("ZResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fZResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("ZResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fPhiResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("PhiResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fPhiResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("PhiResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);
  //
  hiss = comp->MakeResol(comp->fThetaResolS1PtTPC,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("ThetaResolS1PtTPC");
  
  aFolderObj->Add(hiss);

  hiss = comp->MakeResol(comp->fThetaResolS1PtTPCITS,1,0);
  hiss->SetXTitle("#sqrt(1/mcp_{t})");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->SetName("ThetaResolS1PtTPCITS");
  
  aFolderObj->Add(hiss);

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliComparisonRes::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliComparisonRes * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}

//_____________________________________________________________________________
Long64_t AliComparisonRes::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
  AliComparisonRes* entry = dynamic_cast<AliComparisonRes*>(obj);
  if (entry == 0) continue; 
  
  fMCVertex->Add(entry->fMCVertex);
  fRecVertex->Add(entry->fRecVertex);

  fPhiTanPtTPC->Add(entry->fPhiTanPtTPC);
  fPhiTanPtTPCITS->Add(entry->fPhiTanPtTPCITS);

  fPtResolTPC->Add(entry->fPtResolTPC);
  fPtPullTPC->Add(entry->fPtPullTPC);
  fPhiResolTanTPC->Add(entry->fPhiResolTanTPC);
  fTanResolTanTPC->Add(entry->fTanResolTanTPC);
  fPhiPullTanTPC->Add(entry->fPhiPullTanTPC);
  fTanPullTanTPC->Add(entry->fTanPullTanTPC);

  fPtResolTPCITS->Add(entry->fPtResolTPCITS);
  fPtPullTPCITS->Add(entry->fPtPullTPCITS);
  fPhiResolTanTPCITS->Add(entry->fPhiResolTanTPCITS);
  fTanResolTanTPCITS->Add(entry->fTanResolTanTPCITS);
  fPhiPullTanTPCITS->Add(entry->fPhiPullTanTPCITS);
  fTanPullTanTPCITS->Add(entry->fTanPullTanTPCITS);

  // Histograms for 1/pt parameterisation
  f1Pt2ResolS1PtTPC->Add(entry->f1Pt2ResolS1PtTPC);
  fYResolS1PtTPC->Add(entry->fYResolS1PtTPC);
  fZResolS1PtTPC->Add(entry->fZResolS1PtTPC);
  fPhiResolS1PtTPC->Add(entry->fPhiResolS1PtTPC);
  fThetaResolS1PtTPC->Add(entry->fThetaResolS1PtTPC);

  f1Pt2ResolS1PtTPCITS->Add(entry->f1Pt2ResolS1PtTPCITS);
  fYResolS1PtTPCITS->Add(entry->fYResolS1PtTPCITS);
  fZResolS1PtTPCITS->Add(entry->fZResolS1PtTPCITS);
  fPhiResolS1PtTPCITS->Add(entry->fPhiResolS1PtTPCITS);
  fThetaResolS1PtTPCITS->Add(entry->fThetaResolS1PtTPCITS);

  // Resolution histograms (constrained param)
  fCPhiResolTan->Add(entry->fCPhiResolTan);
  fCTanResolTan->Add(entry->fCTanResolTan);
  fCPtResolTan->Add(entry->fCPtResolTan);
  fCPhiPullTan->Add(entry->fCPhiPullTan);
  fCTanPullTan->Add(entry->fCTanPullTan);
  fCPtPullTan->Add(entry->fCPtPullTan);

  //  Histograms for 1/pt parameterisation (constrained)
  fC1Pt2ResolS1PtTPC->Add(entry->fC1Pt2ResolS1PtTPC);
  fCYResolS1PtTPC->Add(entry->fCYResolS1PtTPC);
  fCZResolS1PtTPC->Add(entry->fCZResolS1PtTPC);
  fCPhiResolS1PtTPC->Add(entry->fCPhiResolS1PtTPC);
  fCThetaResolS1PtTPC->Add(entry->fCThetaResolS1PtTPC);

  fC1Pt2ResolS1PtTPCITS->Add(entry->fC1Pt2ResolS1PtTPCITS);
  fCYResolS1PtTPCITS->Add(entry->fCYResolS1PtTPCITS);
  fCZResolS1PtTPCITS->Add(entry->fCZResolS1PtTPCITS);
  fCPhiResolS1PtTPCITS->Add(entry->fCPhiResolS1PtTPCITS);
  fCThetaResolS1PtTPCITS->Add(entry->fCThetaResolS1PtTPCITS);

  count++;
  }

return count;
}

//_____________________________________________________________________________
TFolder* AliComparisonRes::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
