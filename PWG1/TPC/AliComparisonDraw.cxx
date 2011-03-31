

//
// Comparison draw
// Compare the MC information with the reconstructed 
//

/*
  after running analysis, read the file, and get component
  gSystem->Load("libPWG1.so");
  TFile f("Output.root");
  AliComparisonDraw * comp = (AliComparisonDraw*)f.Get("AliComparisonDraw");
  TF1 fl("fl","((min(250./(abs(x+0.000001)),250)-90))",0,2);  // length function
  TF1 fl2("fl2","[0]/((min(250./(abs(x+0.000001)),250)-90))^[1]",0,2);
  fl2.SetParameter(1,1);
  fl2.SetParameter(0,1);

*/




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
//
// 
#include "AliESDEvent.h"   // new container
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
//
#include "AliMathBase.h"
#include "AliTreeDraw.h" 

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonDraw.h" 


ClassImp(AliComparisonDraw)

Bool_t    AliComparisonDraw::fgBDraw=kFALSE;         //option draw temporary results

AliComparisonDraw::AliComparisonDraw():
  TNamed("ComparisonDraw","ComparisonDraw"),
  fEffTPCPt(0),      // TPC efficiency as function of Pt (tan+-1)
  fEffTPCPtMC(0),    // MC -TPC efficiency as function of Pt (tan+-1)
  fEffTPCPtF(0),     // efficiency for findable tracks
  //
  fEffTPCTan(0),   // TPC efficiency as function of Tan (pt>0.15
  fEffTPCTanMC(0), // MC -TPC efficiency as function of Tan (pt>0.15)
  fEffTPCTanF(0),  // efficiency for findable tracks Tan (pt>0.15)
  //
  fEffTPCPtTan(0),    // TPC efficiency as function of Pt and tan
  fEffTPCPtTanMC(0),  // MC -TPC efficiency as function of Pt and tan
  fEffTPCPtTanF(0),  // TPC efficiency as function of Pt and tan
  //
  // dEdx resolution
  //
  fTPCSignalNormTan(0), // tpc signal normalized to the mean signal - MC
  fTPCSignalNormSPhi(0),   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTPhi(0),   // tpc signal normalized to the mean signal - MC
  //
  fTPCSignalNormTanSPhi(0),   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTanTPhi(0),   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTanSPt(0),   // tpc signal normalized to the mean signal - MC
  //
  //
  fPtResolLPT(0),        // pt resolution - low pt
  fPtResolHPT(0),        // pt resolution - high pt 
  fPtPullLPT(0),         // pt resolution - low pt
  fPtPullHPT(0),         // pt resolution - high pt 
  //
  // Resolution constrained param
  //
  fCPhiResolTan(0),   // angular resolution -  constrained
  fCTanResolTan(0),   // angular resolution -  constrained
  fCPtResolTan(0),    // pt resolution      -  constrained
  fCPhiPullTan(0),   // angular resolution -  constrained
  fCTanPullTan(0),   // angular resolution -  constrained
  fCPtPullTan(0),    // pt resolution      -  constrained
  //
  // DCA resolution
  //
  fD0TanSPtB1(0),   // distance to vertex y  
  fD1TanSPtB1(0),   // distance to vertex z  
  fD0TanSPtL1(0),   // distance to vertex y  
  fD1TanSPtL1(0)   // distance to vertex z  
{
  InitHisto();
}

AliComparisonDraw::AliComparisonDraw(const AliComparisonDraw& draw):
  TNamed(draw.GetName(),draw.GetTitle()),
  fEffTPCPt(draw.fEffTPCPt),      // TPC efficiency as function of Pt (tan+-1)
  fEffTPCPtMC(draw.fEffTPCPtMC),    // MC -TPC efficiency as function of Pt (tan+-1)
  fEffTPCPtF(draw.fEffTPCPtF),     // efficiency for findable tracks
  //
  fEffTPCTan(draw.fEffTPCTan),   // TPC efficiency as function of Tan (pt>0.15
  fEffTPCTanMC(draw.fEffTPCTanMC), // MC -TPC efficiency as function of Tan (pt>0.15)
  fEffTPCTanF(draw.fEffTPCTanF),  // efficiency for findable tracks Tan (pt>0.15)
  //
  fEffTPCPtTan(draw.fEffTPCPtTan),    // TPC efficiency as function of Pt and tan
  fEffTPCPtTanMC(draw.fEffTPCPtTanMC),  // MC -TPC efficiency as function of Pt and tan
  fEffTPCPtTanF(draw.fEffTPCPtTanF),  // TPC efficiency as function of Pt and tan
  //
  // dEdx resolution
  //
  fTPCSignalNormTan(draw.fTPCSignalNormTan), // tpc signal normalized to the mean signal - MC
  fTPCSignalNormSPhi(draw.fTPCSignalNormSPhi),   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTPhi(draw.fTPCSignalNormTPhi),   // tpc signal normalized to the mean signal - MC
  //
  fTPCSignalNormTanSPhi(draw.fTPCSignalNormTanSPhi),   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTanTPhi(draw.fTPCSignalNormTanTPhi),   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTanSPt(draw.fTPCSignalNormTanSPt),   // tpc signal normalized to the mean signal - MC
  //
  //
  fPtResolLPT(draw.fPtResolLPT),        // pt resolution - low pt
  fPtResolHPT(draw.fPtResolHPT),        // pt resolution - high pt 
  fPtPullLPT(draw.fPtPullLPT),         // pt resolution - low pt
  fPtPullHPT(draw.fPtPullHPT),         // pt resolution - high pt 
  //
  // Resolution constrained param
  //
  fCPhiResolTan(draw.fCPhiResolTan),   // angular resolution -  constrained
  fCTanResolTan(draw.fCTanResolTan),   // angular resolution -  constrained
  fCPtResolTan(draw.fCPtResolTan),    // pt resolution      -  constrained
  fCPhiPullTan(draw.fCPhiPullTan),   // angular resolution -  constrained
  fCTanPullTan(draw.fCTanPullTan),   // angular resolution -  constrained
  fCPtPullTan(draw.fCPtPullTan),    // pt resolution      -  constrained
  //
  // DCA resolution
  //
  fD0TanSPtB1(draw.fD0TanSPtB1),   // distance to vertex y  
  fD1TanSPtB1(draw.fD1TanSPtB1),   // distance to vertex z  
  fD0TanSPtL1(draw.fD0TanSPtL1),   // distance to vertex y  
  fD1TanSPtL1(draw.fD1TanSPtL1)   // distance to vertex z  
{
  //
  // copy constructor
  //
}

AliComparisonDraw& AliComparisonDraw::operator=(const AliComparisonDraw& info){
  //
  // assignment operator
  //
  delete this;
  new (this) AliComparisonDraw(info);
  return *this;  
}




AliComparisonDraw::~AliComparisonDraw(){
  //
  //
  //
  delete  fEffTPCPt;      // TPC efficiency as function of Pt (tan+-1)
  delete  fEffTPCPtMC;    // MC -TPC efficiency as function of Pt (tan+-1)
  delete  fEffTPCPtF;     // efficiency for findable tracks
  //
  delete  fEffTPCTan;   // TPC efficiency as function of Tan (pt>0.15
  delete  fEffTPCTanMC; // MC -TPC efficiency as function of Tan (pt>0.15)
  delete  fEffTPCTanF;  // efficiency for findable tracks Tan (pt>0.15)
  //
  delete  fEffTPCPtTan;    // TPC efficiency as function of Pt and tan
  delete  fEffTPCPtTanMC;  // MC -TPC efficiency as function of Pt and tan
  delete  fEffTPCPtTanF;  // TPC efficiency as function of Pt and tan
  //
  // dEdx resolution
  //
  delete  fTPCSignalNormTan; // tpc signal normalized to the mean signal - MC
  delete  fTPCSignalNormSPhi;   // tpc signal normalized to the mean signal - MC
  delete  fTPCSignalNormTPhi;   // tpc signal normalized to the mean signal - MC
  //
  delete  fTPCSignalNormTanSPhi;   // tpc signal normalized to the mean signal - MC
  delete  fTPCSignalNormTanTPhi;   // tpc signal normalized to the mean signal - MC
  delete  fTPCSignalNormTanSPt;   // tpc signal normalized to the mean signal - MC
  //
  //
  delete  fPtResolLPT;        // pt resolution - low pt
  delete  fPtResolHPT;        // pt resolution - high pt 
  delete  fPtPullLPT;         // pt resolution - low pt
  delete  fPtPullHPT;         // pt resolution - high pt 
  //
  // Resolution constrained param
  //
  delete fCPhiResolTan;   // angular resolution -  constrained
  delete fCTanResolTan;   // angular resolution -  constrained
  delete fCPtResolTan;    // pt resolution      -  constrained
  delete fCPhiPullTan;   // angular resolution -  constrained
  delete fCTanPullTan;   // angular resolution -  constrained
  delete fCPtPullTan;    // pt resolution      -  constrained
  //
  // DCA resolution
  //
  delete fD0TanSPtB1;   // distance to vertex y  
  delete fD1TanSPtB1;   // distance to vertex z  
  delete fD0TanSPtL1;   // distance to vertex y  
  delete fD1TanSPtL1;   // distance to vertex z  
 
}




void AliComparisonDraw::InitHisto(){
  //
  //
  // EFFICIENCY
  //  
  // Efficiency as function of pt
  fEffTPCPt = new TProfile("Eff_pt","Eff_Pt",50,0.1,3);            // physical
  fEffTPCPtMC = new TProfile("MC_Eff_pt","MC_Eff_Pt",50,0.1,3);    // MC - particles make more than 50 rowdigits
  fEffTPCPtF = new TProfile("F_Eff_pt","F_Eff_Pt",50,0.1,3);     // tracking - under condition more than 50 rdigits

  // Efficiency as function of pt
  fEffTPCTan = new TProfile("Eff_tan","Eff_tan",50,-2.5,2.5);            // physical
  fEffTPCTanMC = new TProfile("MC_Eff_tan","MC_Eff_tan",50,-2.5,2.5);    // MC - particles make more than 50 rowdigits
  fEffTPCTanF = new TProfile("F_Eff_tan","F_Eff_tan",50,-2.5,2.5);     // tracking - under condition more than 50 rdigits

  fEffTPCPtTan = new TProfile2D("Eff_pt","Eff_Pt",10,0.1,3,20,-2.,2.);
  fEffTPCPtTanMC = new TProfile2D("MC_Eff_pt","MC Eff Pt",10,0.1,3,20, -2.,2.);
  fEffTPCPtTanF = new TProfile2D("MC_Eff_pt","MC Eff Pt",10,0.1,3,20, -2.,2.);
  
  //
  // TPC dEdx
  //
  fTPCSignalNormTan = new TH2F("CdEdxTan","CdEdxTan",50, -2,2,  40,30,70); // tpc signal normalized to the MC
  fTPCSignalNormSPhi   = new TH2F("CdEdxSPhi","CdEdxSPhi",10,0.0,1,40,30,70); // tpc signal normalized to the MC
  fTPCSignalNormTPhi   = new TH2F("CdEdxTPhi","CdEdxTPhi",10,0.0,2,40,30,70); // tpc signal normalized to the MC

  fTPCSignalNormTanSPhi= new TH3F("CdEdxTanSPhi","CdEdxTanSPhi",20, -2,2, 10,0.0 ,1,  40,30,70);   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTanTPhi= new TH3F("CdEdxTanTPhi","CdEdxTanTPhi",20, -2,2, 10,0.0 ,1,  40,30,70);   // tpc signal normalized to the mean signal - MC
  fTPCSignalNormTanSPt= new TH3F("CdEdxTanSPt","CdEdxTanSPt",20, -2,2, 10,0.3 ,3,  40,30,70);   // tpc signal normalized to the mean signal - MC



  //
  // RESOLUTION
  //
  fCPtResolTan = new TH2F("C Pt resol","C pt resol",50, -2,2,200,-0.2,0.2);
  fCPtPullTan = new TH2F("C Pt pull","C pt pull",50, -2,2,200,-5,5);
  //
  fCPhiResolTan = new TH2F("CPhiResolTan","CPhiResolTan",50, -2,2,200,-0.025,0.025);   
  // angular resolution -  constrained
  fCTanResolTan = new TH2F("CTanResolTan","CTanResolTan",50, -2,2,200,-0.025,0.025);
  // angular resolution -  constrained
  fCPtResolTan=new TH2F("CPtResol","CPtResol",50, -2,2,200,-0.2,0.2);;    
  // pt resolution      -  constrained
  fCPhiPullTan = new TH2F("CPhiPullTan","CPhiPullTan",50, -2,2,200,-5,5);   
  // angular resolution -  constrained
  fCTanPullTan = new TH2F("CTanPullTan","CTanPullTan",50, -2,2,200,-5,5);
  // angular resolution -  constrained
  fCPtPullTan=new TH2F("CPtPull","CPtPull",50, -2,2,200,-5,5);    
  // pt resolution      -  constrained
  //
  fPtResolLPT = new TH2F("Pt resol","pt resol",10, 0.1,3,200,-0.2,0.2);
  fPtResolHPT = new TH2F("Pt resol","pt resol",10, 2,100,200,-0.3,0.3);  
  fPtPullLPT = new TH2F("Pt pool","pt pool",10, 0.1,3,200,-6,6);
  fPtPullHPT = new TH2F("Pt pool","pt pool",10, 2,100,200,-6,6);  
  //
  fD0TanSPtB1 = new TH3F("DCAyTanSPt","DCAyTanSPt",20,1,2, 10,0.3,2, 100,-4,4);
  fD1TanSPtB1 = new TH3F("DCAzTanSPt","DCAzTanSPt",20,1,2, 10,0.3,2, 100,-4,4);
  fD0TanSPtL1 = new TH3F("DCAyTanSPt","DCAyTanSPt",20,0,1, 10,0.3,2, 100,-0.1,0.1);
  fD1TanSPtL1 = new TH3F("DCAzTanSPt","DCAzTanSPt",20,0,1, 10,0.3,2, 100, -0.1,0.1);



}

void AliComparisonDraw::ProcessEff(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  // make efficiencies histograms
  //
  Float_t kptcut = 0.15;
  Float_t ktancut=1.;
  Int_t   kmincl =50;
  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Bool_t isPrim = infoMC->GetParticle().R()<0.1 && TMath::Abs(infoMC->GetParticle().Vz())<10;
  //z diamond and 
  
  if (!isPrim) return;

  //pt
  if (TMath::Abs(tantheta)<ktancut){
    fEffTPCPt->Fill(mcpt, infoRC->GetStatus(1)==3);
    fEffTPCPtMC->Fill(mcpt, infoMC->GetRowsWithDigits()>kmincl);
    if (infoMC->GetRowsWithDigits()>kmincl){
      fEffTPCPtF->Fill(mcpt, infoRC->GetStatus(1)==3);
    }
  }

  //theta
  if (TMath::Abs(mcpt)>kptcut){
    fEffTPCTan->Fill(tantheta, infoRC->GetStatus(1)==3);
    fEffTPCTanMC->Fill(tantheta, infoMC->GetRowsWithDigits()>kmincl);
    if (infoMC->GetRowsWithDigits()>kmincl){
      fEffTPCTanF->Fill(tantheta, infoRC->GetStatus(1)==3);
    }
  }
  // 
  // pt-theta
  //
  fEffTPCPtTan->Fill(mcpt,tantheta,infoRC->GetStatus(1)==3);
  fEffTPCPtTanMC->Fill(mcpt,tantheta,infoMC->GetRowsWithDigits()>50); 
  if (infoMC->GetRowsWithDigits()>kmincl){
    fEffTPCPtTanF->Fill(mcpt,tantheta,infoRC->GetStatus(1)==3); 
  }
}


void AliComparisonDraw::ProcessResolConstrained(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  //
  //
  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Bool_t isPrim = infoMC->GetParticle().R()<0.1 && TMath::Abs(infoMC->GetParticle().Vz())<10;
  //z diamond and 
  
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<10) return;
  if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;
  
  //
  // constrained parameters resolution
  //
  const AliExternalTrackParam * cparam = infoRC->GetESDtrack()->GetConstrainedParam();
  Float_t deltaCPt= (mcpt-cparam->Pt())/mcpt;  
  Float_t pullCPt= (1/mcpt-cparam->OneOverPt())/
    TMath::Sqrt(cparam->GetSigma1Pt2());          
  Float_t deltaPhi = TMath::ATan2(cparam->Py(),cparam->Px())-
    TMath::ATan2(infoMC->GetParticle().Py(),infoMC->GetParticle().Px());
  Float_t pullPhi = deltaPhi/TMath::Sqrt(cparam->GetSigmaSnp2()); 

  Float_t deltaTan = TMath::ATan2(cparam->Pz(),cparam->Pt())-
    TMath::ATan2(infoMC->GetParticle().Pz(),infoMC->GetParticle().Pt());
  Float_t pullTan = deltaPhi/TMath::Sqrt(cparam->GetSigmaSnp2()); 

  fCPtResolTan->Fill(tantheta,deltaCPt);
  fCPtPullTan->Fill(tantheta,pullCPt);
  fCPhiResolTan->Fill(tantheta,deltaPhi);
  fCPhiPullTan->Fill(tantheta,pullPhi);
  fCTanResolTan->Fill(tantheta,deltaTan);
  fCTanPullTan->Fill(tantheta,pullTan);

}



void  AliComparisonDraw::ProcessTPCdedx(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  //
  //
  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Bool_t isPrim = infoMC->GetParticle().R()<0.1 && TMath::Abs(infoMC->GetParticle().Vz())<10;
  //z diamond and 
  
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<10) return;
  if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;
  Float_t mprim = infoMC->GetPrim();
  if (mprim>1.4) return;
  if (mprim<0.5) return;
  if (infoRC->GetESDtrack()->GetTPCsignalN()<50) return;
  //
  Float_t ratio = infoRC->GetESDtrack()->GetTPCsignal()/infoMC->GetPrim();
  Float_t sphi =  infoRC->GetESDtrack()->GetInnerParam()->GetSnp();
  Float_t tphi =  sphi/TMath::Sqrt((1.-sphi)*(1.+sphi));


  if (TMath::Abs(infoMC->GetParticle().GetPdgCode())!=211) return;
  if (mcpt>0.5){
    fTPCSignalNormTan->Fill(tantheta,ratio);    //only subset
  }
  if (TMath::Abs(tantheta)<0.5){
    fTPCSignalNormSPhi->Fill(sphi,ratio);        // only subset
    fTPCSignalNormTPhi->Fill(tphi,ratio);        // only subset
  }
  fTPCSignalNormTanSPhi->Fill(tantheta,sphi,ratio);    
  fTPCSignalNormTanTPhi->Fill(tantheta,tphi,ratio);    
  fTPCSignalNormTanSPt->Fill(tantheta,TMath::Sqrt(mcpt),ratio);    
}

void      AliComparisonDraw::ProcessDCA(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  //
  //
  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Bool_t isPrim = infoMC->GetParticle().R()<0.1 && TMath::Abs(infoMC->GetParticle().Vz())<10;
  //z diamond and 
  if (!isPrim) return;
  if (infoRC->GetStatus(1)!=3) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<10) return;
  if (!infoRC->GetESDtrack()->GetConstrainedParam()) return;
  Float_t spt = TMath::Sqrt(mcpt);
  Float_t dca[2],cov[3];
  infoRC->GetESDtrack()->GetImpactParameters(dca,cov);
  Int_t clusterITS[100];
  if (infoRC->GetESDtrack()->GetITSclusters(clusterITS)==0){
    fD0TanSPtB1->Fill(tantheta,spt,dca[0]);
    fD1TanSPtB1->Fill(tantheta,spt,dca[1]);
  }
  fD0TanSPtL1->Fill(tantheta,spt,dca[0]);
  fD1TanSPtL1->Fill(tantheta,spt,dca[1]);  
}




void AliComparisonDraw::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  // 
  //
  ProcessEff(infoMC,infoRC);
  ProcessResolConstrained(infoMC,infoRC);
  ProcessTPCdedx(infoMC, infoRC);
  ProcessDCA(infoMC, infoRC);

  Float_t mcpt = infoMC->GetParticle().Pt();
  Bool_t isPrim = infoMC->GetParticle().R()<0.1 && TMath::Abs(infoMC->GetParticle().Vz())<10;
  //z diamond and 
  
  if (!isPrim) return;
  //
  //
  if (infoRC->GetStatus(1)==0) return;
  if (!infoRC->GetESDtrack()) return;  
  if (infoRC->GetESDtrack()->GetTPCNcls()<10) return;
  //  printf("Pt\t%f\t%f\n",mcpt, infoRC->GetESDtrack()->Pt());
  
  Float_t deltaPt= (mcpt-infoRC->GetESDtrack()->Pt())/mcpt;  
  Float_t poolPt= (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/
    TMath::Sqrt(infoRC->GetESDtrack()->GetSigma1Pt2());  

  fPtResolLPT->Fill(mcpt,deltaPt);
  fPtResolHPT->Fill(mcpt,deltaPt);
  fPtPullLPT->Fill(mcpt,poolPt);
  fPtPullHPT->Fill(mcpt,poolPt);  
}



TH1F* AliComparisonDraw::MakeResol(TH2F * his, Int_t integ, Bool_t type){
  TH1F *hisr, *hism;
  if (!gPad) new TCanvas;
  hisr = AliTreeDraw::CreateResHistoI(his,&hism,integ);
  if (type) return hism;
  else 
    return hisr;
}


TGraph2D * AliComparisonDraw::MakeStat2D(TH3 * his, Int_t delta0, Int_t delta1, Int_t type){
  //
  //
  //
  // delta - number of bins to integrate
  // type - 0 - mean value

  TAxis * xaxis  = his->GetXaxis();
  TAxis * yaxis  = his->GetYaxis();
  //  TAxis * zaxis  = his->GetZaxis();
  Int_t   nbinx  = xaxis->GetNbins();
  Int_t   nbiny  = yaxis->GetNbins();
  char name[1000];
  Int_t icount=0;
  TGraph2D  *graph = new TGraph2D(nbinx*nbiny);
  TF1 f1("f1","gaus");
  for (Int_t ix=0; ix<nbinx;ix++)
    for (Int_t iy=0; iy<nbiny;iy++){
      Float_t xcenter = xaxis->GetBinCenter(ix); 
      Float_t ycenter = yaxis->GetBinCenter(iy); 
      snprintf(name,1000,"%s_%d_%d",his->GetName(), ix,iy);
      TH1 *projection = his->ProjectionZ(name,ix-delta0,ix+delta0,iy-delta1,iy+delta1);
      Float_t stat= 0;
      if (type==0) stat = projection->GetMean();
      if (type==1) stat = projection->GetRMS();
      if (type==2 || type==3){
	TVectorD vec(3);
	AliMathBase::LTM((TH1F*)projection,&vec,0.7);
	if (type==2) stat= vec[1];
	if (type==3) stat= vec[0];	
      }
      if (type==4|| type==5){
	projection->Fit(&f1);
	if (type==4) stat= f1.GetParameter(1);
	if (type==5) stat= f1.GetParameter(2);
      }
      //printf("%d\t%f\t%f\t%f\n", icount,xcenter, ycenter, stat);
      graph->SetPoint(icount,xcenter, ycenter, stat);
      icount++;
    }
  return graph;
}

TGraph * AliComparisonDraw::MakeStat1D(TH3 * his, Int_t delta1, Int_t type){
  //
  //
  //
  // delta - number of bins to integrate
  // type - 0 - mean value

  TAxis * xaxis  = his->GetXaxis();
  TAxis * yaxis  = his->GetYaxis();
  //  TAxis * zaxis  = his->GetZaxis();
  Int_t   nbinx  = xaxis->GetNbins();
  Int_t   nbiny  = yaxis->GetNbins();
  char name[1000];
  Int_t icount=0;
  TGraph  *graph = new TGraph(nbinx);
  TF1 f1("f1","gaus");
  for (Int_t ix=0; ix<nbinx;ix++){
    Float_t xcenter = xaxis->GetBinCenter(ix); 
    //    Float_t ycenter = yaxis->GetBinCenter(iy); 
    snprintf(name,1000,"%s_%d",his->GetName(), ix);
    TH1 *projection = his->ProjectionZ(name,ix-delta1,ix+delta1,0,nbiny);
    Float_t stat= 0;
    if (type==0) stat = projection->GetMean();
    if (type==1) stat = projection->GetRMS();
    if (type==2 || type==3){
      TVectorD vec(3);
	AliMathBase::LTM((TH1F*)projection,&vec,0.7);
	if (type==2) stat= vec[1];
	if (type==3) stat= vec[0];	
    }
    if (type==4|| type==5){
      projection->Fit(&f1);
      if (type==4) stat= f1.GetParameter(1);
      if (type==5) stat= f1.GetParameter(2);
    }
      //printf("%d\t%f\t%f\t%f\n", icount,xcenter, ycenter, stat);
    graph->SetPoint(icount,xcenter, stat);
    icount++;
  }
  return graph;
}

//
// Make derived plots
//

void AliComparisonDraw::MakePlots(){
  //
  //
  //
  AliComparisonDraw * comp=this;

  TFile *fp = new TFile("picutures.root","recreate");
  TH1F *hiss=0;
  //TH1F *hism=0;
  TGraph2D * gr=0, gr2=0;
  TGraph * gr0 = 0;
  TCanvas * c = new TCanvas("Phi resol Tan","Phi resol Tan");
  //
  //
  //
  hiss = comp->MakeResol(comp->fCPtResolTan,1,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigmap_{t}/p_{t}");
  hiss->Draw(); 
  hiss->Write("CptResolTan");
  //
  //
  hiss = comp->MakeResol(comp->fCPhiResolTan,1,0);
  c->cd();
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#phi (rad)");
  hiss->Draw();
  fp->cd();
  hiss->Write("PhiResolTan");
  //
  hiss = comp->MakeResol(comp->fCTanResolTan,1,0);
  c->cd();
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma#theta (rad)");
  hiss->Draw();
  fp->cd();
  hiss->Write("ThetaResolTan");
  //
  //
  hiss = comp->MakeResol(comp->fCTanResolTan,1,0);
  c->cd();
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigmap_{t}/p_{t} ");
  hiss->Draw();
  fp->cd();
  //
  //
  //
  hiss = comp->MakeResol(comp->fTPCSignalNormTan,4,0);
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("#sigma_{dEdx}");
  hiss->Draw();
  fp->cd();
  hiss->Write("TPCdEdxResolTan");
  //
  //
  //
  hiss = comp->MakeResol(comp->fTPCSignalNormTan,4,1); 
  hiss->SetXTitle("Tan(#theta)");
  hiss->SetYTitle("<dEdx>");
  hiss->Draw(); 
  hiss->Write("TPCdEdxMeanTan");
  //
  //
  gr = comp->MakeStat2D(comp->fTPCSignalNormTanSPt,3,1,4);
  gr->GetXaxis()->SetTitle("Tan(#theta)");
  gr->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr->GetZaxis()->SetTitle("<dEdx>");
  gr->Draw("colz"); 
  gr->GetHistogram()->Write("TPCdEdxMeanTanPt");
  //
  //
  gr = comp->MakeStat2D(comp->fTPCSignalNormTanSPt,3,1,5);
  gr->GetXaxis()->SetTitle("Tan(#theta)");
  gr->GetYaxis()->SetTitle("#sqrt{p_{t}(GeV)}");
  gr->GetZaxis()->SetTitle("#sigma_{dEdx}");
  gr->Draw("colz"); 
  gr->GetHistogram()->Write("TPCdEdxMeanTanPt");
  //
  //
  //
  comp->fEffTPCTanF->SetXTitle("Tan(#theta)");
  comp->fEffTPCTanF->SetYTitle("eff_{findable}");
  comp->fEffTPCTanF->Draw();
  comp->fEffTPCTanF->Write("EffTanFindable");
  //
  //
  comp->fEffTPCTan->SetXTitle("Tan(#theta)");
  comp->fEffTPCTan->SetYTitle("eff_{all}");
  comp->fEffTPCTan->Draw();
  comp->fEffTPCTan->Write("EffTanAll");
  //
  //DCA resolution
  //
  gr0 = comp->MakeStat1D(comp->fD0TanSPtB1,2,5);
  gr0->GetXaxis()->SetTitle("Tan(#theta)");
  gr0->GetYaxis()->SetTitle("#sigmaDCA (cm)");
  gPad->Clear();
  gr0->Draw("al*");
  gr->GetHistogram()->Write("DCAResolTan");
  //
  //
  //
  gr = comp->MakeStat2D(comp->fD0TanSPtB1,4,2,5); 
  gr0->GetXaxis()->SetTitle("Tan(#theta)");
  gr0->GetYaxis()->SetTitle("#sigmaDCA (cm)");
  gPad->Clear();
  gr0->Draw("al*");
  gr->GetHistogram()->Write("DCAResolSPTTan");

  fp->Close();


}



