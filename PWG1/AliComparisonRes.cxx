//------------------------------------------------------------------------------
// Implementation of AliComparisonRes class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms) are stored in the output picture_res.root file.
//  
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
  //after running analysis, read the file, and get component
  gSystem->Load("libPWG1.so");
  TFile f("Output.root");
  AliComparisonRes * comp = (AliComparisonRes*)f.Get("AliComparisonRes");

  // analyse comparison data (output stored in pictures_res.root)
  comp->Analyse();
  
  // paramtetrisation of the TPC track length (for information only) 
  TF1 fl("fl","((min(250./(abs(x+0.000001)),250)-90))",0,2);  // TPC track length function
  TF1 fl2("fl2","[0]/((min(250./(abs(x+0.000001)),250)-90))^[1]",0,2);
  fl2.SetParameter(1,1);
  fl2.SetParameter(0,1);
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
  TNamed("AliComparisonRes","AliComparisonRes"),

  // Resolution 
  fPtResolLPT(0),        // pt resolution - low pt
  fPtResolHPT(0),        // pt resolution - high pt 
  fPtPullLPT(0),         // pt resolution - low pt
  fPtPullHPT(0),         // pt resolution - high pt 
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

  f1Pt2Resol1PtTPC(0),
  f1Pt2Resol1PtTPCITS(0),
  fYResol1PtTPC(0),
  fYResol1PtTPCITS(0),
  fZResol1PtTPC(0),
  fZResol1PtTPCITS(0),
  fPhiResol1PtTPC(0),
  fPhiResol1PtTPCITS(0),
  fThetaResol1PtTPC(0),
  fThetaResol1PtTPCITS(0),

  // constrained
  fC1Pt2Resol1PtTPC(0),
  fC1Pt2Resol1PtTPCITS(0),
  fCYResol1PtTPC(0),
  fCYResol1PtTPCITS(0),
  fCZResol1PtTPC(0),
  fCZResol1PtTPCITS(0),
  fCPhiResol1PtTPC(0),
  fCPhiResol1PtTPCITS(0),
  fCThetaResol1PtTPC(0),
  fCThetaResol1PtTPCITS(0),

  // vertex
  fVertex(0),
 
  // Cuts 
  fCutsRC(0),  
  fCutsMC(0)  
{
  InitHisto();
  InitCuts();
  
  // vertex (0,0,0)
  fVertex = new AliESDVertex();
  fVertex->SetXv(0.0);
  fVertex->SetYv(0.0);
  fVertex->SetZv(0.0);
}

//_____________________________________________________________________________
AliComparisonRes::~AliComparisonRes(){
  
  // Resolution histograms
  if(fPtResolLPT) delete  fPtResolLPT; fPtResolLPT=0;     
  if(fPtResolHPT) delete  fPtResolHPT; fPtResolHPT=0;    
  if(fPtPullLPT)  delete  fPtPullLPT;  fPtPullLPT=0;    
  if(fPtPullHPT)  delete  fPtPullHPT;  fPtPullHPT=0;   

  // Resolution histograms (constrained param)
  if(fCPhiResolTan) delete fCPhiResolTan; fCPhiResolTan=0;
  if(fCTanResolTan) delete fCTanResolTan; fCTanResolTan=0;
  if(fCPtResolTan)  delete fCPtResolTan;  fCPtResolTan=0; 
  if(fCPhiPullTan)  delete fCPhiPullTan;  fCPhiPullTan=0;
  if(fCTanPullTan)  delete fCTanPullTan;  fCTanPullTan=0;
  if(fCPtPullTan)   delete fCPtPullTan;   fCPtPullTan=0;

  // Parametrisation histograms
  // 
  if(f1Pt2Resol1PtTPC) delete f1Pt2Resol1PtTPC; f1Pt2Resol1PtTPC=0;
  if(f1Pt2Resol1PtTPCITS) delete f1Pt2Resol1PtTPCITS; f1Pt2Resol1PtTPCITS=0;
  if(fYResol1PtTPC) delete fYResol1PtTPC; fYResol1PtTPC=0;
  if(fYResol1PtTPCITS) delete fYResol1PtTPCITS; fYResol1PtTPCITS=0;
  if(fZResol1PtTPC) delete fZResol1PtTPC; fZResol1PtTPC=0;
  if(fZResol1PtTPCITS) delete fZResol1PtTPCITS; fZResol1PtTPCITS=0;
  if(fPhiResol1PtTPC) delete fPhiResol1PtTPC; fPhiResol1PtTPC=0;
  if(fPhiResol1PtTPCITS) delete fPhiResol1PtTPCITS; fPhiResol1PtTPCITS=0;
  if(fThetaResol1PtTPC) delete fThetaResol1PtTPC; fThetaResol1PtTPC=0;
  if(fThetaResol1PtTPCITS) delete fThetaResol1PtTPCITS; fThetaResol1PtTPCITS=0;

  // constrained
  if(fC1Pt2Resol1PtTPC) delete fC1Pt2Resol1PtTPC; fC1Pt2Resol1PtTPC=0;
  if(fC1Pt2Resol1PtTPCITS) delete fC1Pt2Resol1PtTPCITS; fC1Pt2Resol1PtTPCITS=0;
  if(fCYResol1PtTPC) delete fCYResol1PtTPC; fCYResol1PtTPC=0;
  if(fCYResol1PtTPCITS) delete fCYResol1PtTPCITS; fCYResol1PtTPCITS=0;
  if(fCZResol1PtTPC) delete fCZResol1PtTPC; fCZResol1PtTPC=0;
  if(fCZResol1PtTPCITS) delete fCZResol1PtTPCITS; fCZResol1PtTPCITS=0;
  if(fCPhiResol1PtTPC) delete fCPhiResol1PtTPC; fCPhiResol1PtTPC=0;
  if(fCPhiResol1PtTPCITS) delete fCPhiResol1PtTPCITS; fCPhiResol1PtTPCITS=0;
  if(fCThetaResol1PtTPC) delete fCThetaResol1PtTPC; fCThetaResol1PtTPC=0;
  if(fCThetaResol1PtTPCITS) delete fCThetaResol1PtTPCITS; fCThetaResol1PtTPCITS=0;

  if(fVertex) delete fVertex; fVertex=0;

}

//_____________________________________________________________________________
void AliComparisonRes::InitHisto(){

  // Init histograms
  fCPhiResolTan = new TH2F("CPhiResolTan","CPhiResolTan",50, -2,2,200,-0.025,0.025);   
  fCPhiResolTan->SetXTitle("tan(#theta)");
  fCPhiResolTan->SetYTitle("#Delta#phi");

  fCTanResolTan = new TH2F("CTanResolTan","CTanResolTan",50, -2,2,200,-0.025,0.025);
  fCTanResolTan->SetXTitle("tan(#theta)");
  fCTanResolTan->SetYTitle("#Delta#theta");

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

  fPtResolLPT = new TH2F("Pt_resol_lpt","pt resol",10, 0.1,3,200,-0.2,0.2);
  fPtResolLPT->SetXTitle("p_{t}");
  fPtResolLPT->SetYTitle("#Deltap_{t}/p_{t}");

  fPtResolHPT = new TH2F("Pt_resol_hpt","pt resol",10, 2,100,200,-0.3,0.3);  
  fPtResolHPT->SetXTitle("p_{t}");
  fPtResolHPT->SetYTitle("#Deltap_{t}/p_{t}");

  fPtPullLPT = new TH2F("Pt_pull_lpt","pt pull",10, 0.1,3,200,-6,6);
  fPtPullLPT->SetXTitle("p_{t}");
  fPtPullLPT->SetYTitle("#Deltap_{t}/#Sigma");

  fPtPullHPT = new TH2F("Pt_pull_hpt","pt pull",10,2,100,200,-6,6);  
  fPtPullHPT->SetXTitle("p_{t}");
  fPtPullHPT->SetYTitle("#Deltap_{t}/#Sigma");

  //
  // Parametrisation histograms
  // 

  f1Pt2Resol1PtTPC = new TH2F("f1Pt2Resol1PtTPC","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs 1/pt)",100,0,10,200,-0.010,0.010);  
  f1Pt2Resol1PtTPC->SetXTitle("1/mcp_{t}");
  f1Pt2Resol1PtTPC->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  f1Pt2Resol1PtTPCITS = new TH2F("f1Pt2Resol1PtTPCITS","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs 1/pt)",100,0,10,200,-0.010,0.010);  
  f1Pt2Resol1PtTPCITS->SetXTitle("1/mcp_{t}");
  f1Pt2Resol1PtTPCITS->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  fYResol1PtTPC = new TH2F("fYResol1PtTPC","fYResol1PtTPC",100, 0,10,200,-1.0,1.0);   
  fYResol1PtTPC->SetXTitle("1/mcpt");
  fYResol1PtTPC->SetYTitle("#DeltaY");

  fYResol1PtTPCITS = new TH2F("fYResol1PtTPCITS","fYResol1PtTPCITS",100, 0,10,200,-0.05,0.05);   
  fYResol1PtTPCITS->SetXTitle("1/mcpt");
  fYResol1PtTPCITS->SetYTitle("#DeltaY");

  fZResol1PtTPC = new TH2F("fZResol1PtTPC","fZResol1PtTPC",100, 0,10,200,-1.0,1.0);   
  fZResol1PtTPC->SetXTitle("1/mcpt");
  fZResol1PtTPC->SetYTitle("#DeltaZ");

  fZResol1PtTPCITS = new TH2F("fZResol1PtTPCITS","fZResol1PtTPCITS",100, 0,10,200,-0.05,0.05);   
  fZResol1PtTPCITS->SetXTitle("1/mcpt");
  fZResol1PtTPCITS->SetYTitle("#DeltaZ");

  fPhiResol1PtTPC = new TH2F("fPhiResol1PtTPC","fPhiResol1PtTPC",100, 0,10,200,-0.025,0.025);   
  fPhiResol1PtTPC->SetXTitle("1/mcpt");
  fPhiResol1PtTPC->SetYTitle("#Delta#phi");

  fPhiResol1PtTPCITS = new TH2F("fPhiResol1PtTPCITS","fPhiResol1PtTPCITS",100, 0,10,200,-0.01,0.01);   
  fPhiResol1PtTPCITS->SetXTitle("1/mcpt");
  fPhiResol1PtTPCITS->SetYTitle("#Delta#phi");

  fThetaResol1PtTPC = new TH2F("fThetaResol1PtTPC","fThetaResol1PtTPC",100, 0,10,200,-0.025,0.025);   
  fThetaResol1PtTPC->SetXTitle("1/mcpt");
  fThetaResol1PtTPC->SetYTitle("#Delta#theta");

  fThetaResol1PtTPCITS = new TH2F("fThetaResol1PtTPCITS","fThetaResol1PtTPCITS",100, 0,10,200,-0.01,0.01);   
  fThetaResol1PtTPCITS->SetXTitle("1/mcpt");
  fThetaResol1PtTPCITS->SetYTitle("#Delta#theta");
  
  // constrained
  fC1Pt2Resol1PtTPC = new TH2F("fC1Pt2Resol1PtTPC","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs 1/pt)",100,0,10,200,-0.010,0.010);  
  fC1Pt2Resol1PtTPC->SetXTitle("1/mcp_{t}");
  fC1Pt2Resol1PtTPC->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  fC1Pt2Resol1PtTPCITS = new TH2F("fC1Pt2Resol1PtTPCITS","(1/mcpt-1/pt)/(1+1/mcpt)^2 vs 1/pt)",100,0,10,200,-0.010,0.010);  
  fC1Pt2Resol1PtTPCITS->SetXTitle("1/mcp_{t}");
  fC1Pt2Resol1PtTPCITS->SetYTitle("(1/mcp_{t}-1/p_{t})/(1+1/mcp_{t})^2)");

  fCYResol1PtTPC = new TH2F("fCYResol1PtTPC","fCYResol1PtTPC",100, 0,10,200,-1.0,1.0);   
  fCYResol1PtTPC->SetXTitle("1/mcpt");
  fCYResol1PtTPC->SetYTitle("#DeltaY");

  fCYResol1PtTPCITS = new TH2F("fCYResol1PtTPCITS","fCYResol1PtTPCITS",100, 0,10,200,-0.05,0.05);   
  fCYResol1PtTPCITS->SetXTitle("1/mcpt");
  fCYResol1PtTPCITS->SetYTitle("#DeltaY");

  fCZResol1PtTPC = new TH2F("fCZResol1PtTPC","fCZResol1PtTPC",100, 0,10,200,-1.0,1.0);   
  fCZResol1PtTPC->SetXTitle("1/mcpt");
  fCZResol1PtTPC->SetYTitle("#DeltaZ");

  fCZResol1PtTPCITS = new TH2F("fCZResol1PtTPCITS","fCZResol1PtTPCITS",100, 0,10,200,-0.05,0.05);   
  fCZResol1PtTPCITS->SetXTitle("1/mcpt");
  fCZResol1PtTPCITS->SetYTitle("#DeltaZ");

  fCPhiResol1PtTPC = new TH2F("fCPhiResol1PtTPC","fCPhiResol1PtTPC",100, 0,10,200,-0.025,0.025);   
  fCPhiResol1PtTPC->SetXTitle("1/mcpt");
  fCPhiResol1PtTPC->SetYTitle("#Delta#phi");

  fCPhiResol1PtTPCITS = new TH2F("fCPhiResol1PtTPCITS","fCPhiResol1PtTPCITS",100, 0,10,200,-0.01,0.01);   
  fCPhiResol1PtTPCITS->SetXTitle("1/mcpt");
  fCPhiResol1PtTPCITS->SetYTitle("#Delta#phi");

  fCThetaResol1PtTPC = new TH2F("fCThetaResol1PtTPC","fCThetaResol1PtTPC",100, 0,10,200,-0.025,0.025);   
  fCThetaResol1PtTPC->SetXTitle("1/mcpt");
  fCThetaResol1PtTPC->SetYTitle("#Delta#theta");

  fCThetaResol1PtTPCITS = new TH2F("fCThetaResol1PtTPCITS","fCThetaResol1PtTPCITS",100, 0,10,200,-0.01,0.01);   
  fCThetaResol1PtTPCITS->SetXTitle("1/mcpt");
  fCThetaResol1PtTPCITS->SetYTitle("#Delta#theta");
}

//_____________________________________________________________________________
void AliComparisonRes::InitCuts()
{
  // Init cuts 
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");
}

//_____________________________________________________________________________
void AliComparisonRes::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill resolution comparison information 
  AliExternalTrackParam *track = 0;
  Double_t kRadius    = 3.0;      // beam pipe radius
  Double_t kMaxStep   = 5.0;      // max step
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance

  Int_t clusterITS[200];
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  Float_t deltaPt, pullPt, deltaPhi, deltaTan, delta1Pt2, deltaY1Pt, deltaZ1Pt, deltaPhi1Pt, deltaTheta1Pt; 
  Float_t deltaPtTPC, pullPtTPC, deltaPhiTPC, deltaTanTPC, delta1Pt2TPC, deltaY1PtTPC, deltaZ1PtTPC, deltaPhi1PtTPC, deltaTheta1PtTPC; 

  Float_t mcpt = infoMC->GetParticle().Pt();

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();

  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 
  if (!isPrim) return;
  //if (infoRC->GetStatus(1)==0) return;
  if (infoRC->GetStatus(1)!=3) return; // TPC refit
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<fCutsRC->GetMinNClustersTPC()) return;

  // calculate and set prim. vertex
  fVertex->SetXv( infoMC->GetParticle().Vx() - dv[0] );
  fVertex->SetYv( infoMC->GetParticle().Vy() - dv[1] );
  fVertex->SetZv( infoMC->GetParticle().Vz() - dv[2] );

  deltaPt= (mcpt-infoRC->GetESDtrack()->Pt())/mcpt;  
  pullPt= (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/TMath::Sqrt(infoRC->GetESDtrack()->GetSigma1Pt2());  
  deltaPhi = TMath::ATan2(infoRC->GetESDtrack()->Py(),infoRC->GetESDtrack()->Px())-
                     TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());

  deltaTan = TMath::ATan2(infoRC->GetESDtrack()->Pz(),infoRC->GetESDtrack()->Pt())-
                     TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());

  delta1Pt2 = (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/TMath::Power(1+1/mcpt,2);       
  deltaY1Pt = (infoMC->GetParticle().Vy()-infoRC->GetESDtrack()->GetY()) / (0.2+1/mcpt);
  deltaZ1Pt = (infoMC->GetParticle().Vz()-infoRC->GetESDtrack()->GetZ()) / (0.2+1/mcpt);
  deltaPhi1Pt = deltaPhi   / (0.1+1/mcpt);
  deltaTheta1Pt = deltaTan / (0.1+1/mcpt);

  //track parameters at the first measured point (TPC) 
  //Double_t param[5],x,alpha; // [0]-Y [cm],[1]-Z [cm],[2]-sin(phi),[3]-tan(theta),[4]-1/pt [1/GeV]   
  //infoRC->GetESDtrack()->GetInnerExternalParameters(alpha,x,param);
  //const AliExternalTrackParam *innerTPC =  infoRC->GetESDtrack()->GetInnerParam(); 
  //const AliExternalTrackParam *innerTPC =  infoRC->GetESDtrack()->GetTPCInnerParam(); 


  // calculate track parameters at vertex
  const AliExternalTrackParam *innerTPC =  0;
  if ((innerTPC = infoRC->GetESDtrack()->GetTPCInnerParam()) != 0)
  {
    if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
    {
      Bool_t bStatus = AliTracker::PropagateTrackTo(track,kRadius,infoMC->GetMass(),kMaxStep,kTRUE);
      Bool_t bDCAStatus = track->PropagateToDCA(fVertex,field,kMaxD,dca,cov);

      // Fill parametrisation histograms (only TPC track)
      if(bStatus && bDCAStatus) 
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

			f1Pt2Resol1PtTPC->Fill(1/mcpt,delta1Pt2TPC);
			fYResol1PtTPC->Fill(1/mcpt,deltaY1PtTPC);
			fZResol1PtTPC->Fill(1/mcpt,deltaZ1PtTPC);
			fPhiResol1PtTPC->Fill(1/mcpt,deltaPhi1PtTPC);
			fThetaResol1PtTPC->Fill(1/mcpt,deltaTheta1PtTPC);
	  }
	  delete track;
    }
  }

  // TPC and ITS (nb. of clusters >2) in the system
  if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>2) 
  {
      f1Pt2Resol1PtTPCITS->Fill(1/mcpt,delta1Pt2);
      fYResol1PtTPCITS->Fill(1/mcpt,deltaY1Pt);
      fZResol1PtTPCITS->Fill(1/mcpt,deltaZ1Pt);
      fPhiResol1PtTPCITS->Fill(1/mcpt,deltaPhi1Pt);
      fThetaResol1PtTPCITS->Fill(1/mcpt,deltaTheta1Pt);
  }

  // Fill histograms
  fPtResolLPT->Fill(mcpt,deltaPt);
  fPtResolHPT->Fill(mcpt,deltaPt);
  fPtPullLPT->Fill(mcpt,pullPt);
  fPtPullHPT->Fill(mcpt,pullPt);  
}

//_____________________________________________________________________________
void AliComparisonRes::ProcessConstrained(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill resolution comparison information (constarained parameters) 
  //
  AliExternalTrackParam *track = 0;
  Double_t kRadius    = 3.0;      // beam pipe radius
  Double_t kMaxStep   = 5.0;      // max step
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance

  Int_t clusterITS[200];
  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
 
  Float_t deltaPt, pullPt, deltaPhi, pullPhi, deltaTan, pullTan, delta1Pt2, deltaY1Pt, deltaZ1Pt, deltaPhi1Pt, deltaTheta1Pt; 
  Float_t deltaPtTPC, pullPtTPC, deltaPhiTPC, deltaTanTPC, delta1Pt2TPC, deltaY1PtTPC, deltaZ1PtTPC, deltaPhi1PtTPC, deltaTheta1PtTPC; 

  Float_t mcpt = infoMC->GetParticle().Pt();
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
  fVertex->SetXv( infoMC->GetParticle().Vx() - dv[0] );
  fVertex->SetYv( infoMC->GetParticle().Vy() - dv[1] );
  fVertex->SetZv( infoMC->GetParticle().Vz() - dv[2] );

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
      Bool_t bStatus = AliTracker::PropagateTrackTo(track,kRadius,infoMC->GetMass(),kMaxStep,kTRUE);
      Bool_t bDCAStatus = track->PropagateToDCA(fVertex,field,kMaxD,dca,cov);

      // Fill parametrisation histograms (only TPC track)
      if(bStatus && bDCAStatus) 
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

          fC1Pt2Resol1PtTPC->Fill(1/mcpt,delta1Pt2TPC);
          fCYResol1PtTPC->Fill(1/mcpt,deltaY1PtTPC);
          fCZResol1PtTPC->Fill(1/mcpt,deltaZ1PtTPC);
          fCPhiResol1PtTPC->Fill(1/mcpt,deltaPhi1PtTPC);
          fCThetaResol1PtTPC->Fill(1/mcpt,deltaTheta1PtTPC);
	  }
	  delete track;
    }
  }

 // TPC and ITS (nb. of clusters >2) in the system
  if(infoRC->GetESDtrack()->GetITSclusters(clusterITS)>2) 
  {
      fC1Pt2Resol1PtTPCITS->Fill(1/mcpt,delta1Pt2);
      fCYResol1PtTPCITS->Fill(1/mcpt,deltaY1Pt);
      fCZResol1PtTPCITS->Fill(1/mcpt,deltaZ1Pt);
      fCPhiResol1PtTPCITS->Fill(1/mcpt,deltaPhi1Pt);
      fCThetaResol1PtTPCITS->Fill(1/mcpt,deltaTheta1Pt);
  }

  // Fill histograms
  fCPtResolTan->Fill(tantheta,deltaPt);
  fCPtPullTan->Fill(tantheta,pullPt);
  fCPhiResolTan->Fill(tantheta,deltaPhi);
  fCPhiPullTan->Fill(tantheta,pullPhi);
  fCTanResolTan->Fill(tantheta,deltaTan);
  fCTanPullTan->Fill(tantheta,pullTan);
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
  // in the "pictures_res.root" file 
  
  AliComparisonRes * comp=this;
  TH1F *hiss=0;

  TFile *fp = new TFile("pictures_res.root","recreate");
  fp->cd();

  TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  c->cd();
  //
  hiss = comp->MakeResol(comp->fCPtResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigmap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->Write("CptResolTan");
  //
  hiss = comp->MakeResol(comp->fCPhiResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#phi (rad)");
  hiss->Draw();
  hiss->Write("PhiResolTan");
  //
  hiss = comp->MakeResol(comp->fCTanResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#theta (rad)");
  hiss->Draw();
  hiss->Write("ThetaResolTan");
  //
  hiss = comp->MakeResol(comp->fCPtPullTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/#Sigma(1/p_{t})");
  hiss->Draw();
  hiss->Write("CptPullTan");
  //
  hiss = comp->MakeResol(comp->fC1Pt2Resol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->Write("C1Pt2Resol1PtTPC");
  fC1Pt2Resol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fC1Pt2Resol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->Write("C1Pt2Resol1PtTPCITS");
  fC1Pt2Resol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fCYResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CYResol1PtTPC");
  fCYResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fCYResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CYResol1PtTPCITS");
  fCYResol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fCZResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CZResol1PtTPC");
  fCZResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fCZResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CZResol1PtTPCITS");
  fCZResol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fCPhiResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CPhiResol1PtTPC");
  fCPhiResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fCPhiResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CPhiResol1PtTPCITS");
  fCPhiResol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fCThetaResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CThetaResol1PtTPC");
  fCThetaResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fCThetaResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("CThetaResol1PtTPCITS");
  fCThetaResol1PtTPCITS->Write();

  //
  hiss = comp->MakeResol(comp->f1Pt2Resol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->Write("OnePt2Resol1PtTPC");
  f1Pt2Resol1PtTPC->Write();

  hiss = comp->MakeResol(comp->f1Pt2Resol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("1/mcp_{t}-1/p_{t}/(1+1/p_{t})^2");
  hiss->Draw();
  hiss->Write("OnePt2Resol1PtTPCITS");
  f1Pt2Resol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fYResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("YResol1PtTPC");
  fYResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fYResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcy-y)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("YResol1PtTPCITS");
  fYResol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fZResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("ZResol1PtTPC");
  fZResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fZResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mcz-z)/(0.2+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("ZResol1PtTPCITS");
  fZResol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fPhiResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("PhiResol1PtTPC");
  fPhiResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fPhiResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#phi-#phi)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("PhiResol1PtTPCITS");
  fPhiResol1PtTPCITS->Write();
  //
  hiss = comp->MakeResol(comp->fThetaResol1PtTPC,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("ThetaResol1PtTPC");
  fThetaResol1PtTPC->Write();

  hiss = comp->MakeResol(comp->fThetaResol1PtTPCITS,1,0);
  hiss->SetXTitle("1/mcp_{t}");
  hiss->SetYTitle("(mc#theta-#theta)/(0.1+1/mcp_{t})");
  hiss->Draw();
  hiss->Write("ThetaResol1PtTPCITS");
  fThetaResol1PtTPCITS->Write();

  fp->Close();
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

  fPtResolLPT->Add(entry->fPtResolLPT);
  fPtResolHPT->Add(entry->fPtResolHPT);
  fPtPullLPT->Add(entry->fPtPullLPT);
  fPtPullHPT->Add(entry->fPtPullHPT);

  // Histograms for 1/pt parameterisation
  f1Pt2Resol1PtTPC->Add(entry->f1Pt2Resol1PtTPC);
  fYResol1PtTPC->Add(entry->fYResol1PtTPC);
  fZResol1PtTPC->Add(entry->fZResol1PtTPC);
  fPhiResol1PtTPC->Add(entry->fPhiResol1PtTPC);
  fThetaResol1PtTPC->Add(entry->fThetaResol1PtTPC);

  f1Pt2Resol1PtTPCITS->Add(entry->f1Pt2Resol1PtTPCITS);
  fYResol1PtTPCITS->Add(entry->fYResol1PtTPCITS);
  fZResol1PtTPCITS->Add(entry->fZResol1PtTPCITS);
  fPhiResol1PtTPCITS->Add(entry->fPhiResol1PtTPCITS);
  fThetaResol1PtTPCITS->Add(entry->fThetaResol1PtTPCITS);

  // Resolution histograms (constrained param)
  fCPhiResolTan->Add(entry->fCPhiResolTan);
  fCTanResolTan->Add(entry->fCTanResolTan);
  fCPtResolTan->Add(entry->fCPtResolTan);
  fCPhiPullTan->Add(entry->fCPhiPullTan);
  fCTanPullTan->Add(entry->fCTanPullTan);
  fCPtPullTan->Add(entry->fCPtPullTan);

  //  Histograms for 1/pt parameterisation (constrained)
  fC1Pt2Resol1PtTPC->Add(entry->fC1Pt2Resol1PtTPC);
  fCYResol1PtTPC->Add(entry->fCYResol1PtTPC);
  fCZResol1PtTPC->Add(entry->fCZResol1PtTPC);
  fCPhiResol1PtTPC->Add(entry->fCPhiResol1PtTPC);
  fCThetaResol1PtTPC->Add(entry->fCThetaResol1PtTPC);

  fC1Pt2Resol1PtTPCITS->Add(entry->fC1Pt2Resol1PtTPCITS);
  fCYResol1PtTPCITS->Add(entry->fCYResol1PtTPCITS);
  fCZResol1PtTPCITS->Add(entry->fCZResol1PtTPCITS);
  fCPhiResol1PtTPCITS->Add(entry->fCPhiResol1PtTPCITS);
  fCThetaResol1PtTPCITS->Add(entry->fCThetaResol1PtTPCITS);

  count++;
  }

return count;
}
