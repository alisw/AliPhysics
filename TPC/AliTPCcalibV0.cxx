
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
#include <TH2F.h>
//
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "AliRunLoader.h"
#include "AliStack.h"



#include <TPDGCode.h>
#include <TStyle.h>
#include "TLinearFitter.h"
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TF1.h"



#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliMathBase.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliTrackPointArray.h"
#include "TCint.h"
#include "AliTPCcalibV0.h"
#include "AliV0.h"





ClassImp(AliTPCcalibV0)


AliTPCcalibV0::AliTPCcalibV0() : 
   AliTPCcalibBase(),
   fStack(0),
   fESD(0),
   fPdg(0),
   fParticles(0),
   fV0s(0),
   fGammas(0),
   fV0type(0),
   fTPCdEdx(0),
   fTPCdEdxPi(0),
   fTPCdEdxEl(0),
   fTPCdEdxP(0)
{
  fPdg = new TDatabasePDG;       
  // create output histograms
  fTPCdEdx   = new TH2F("TPCdEdX",  "dE/dX; BetaGamma; TPC signal (a.u.)", 1000, 0.1, 10000, 300, 0, 300);
  BinLogX(fTPCdEdx); 
  fTPCdEdxPi = new TH2F("TPCdEdXPi","dE/dX; BetaGamma; TPC signal (a.u.)", 1000, 0.1, 10000, 300, 0, 300);
  fTPCdEdxEl = new TH2F("TPCdEdXEl","dE/dX; BetaGamma; TPC signal (a.u.)", 1000, 0.1, 10000, 300, 0, 300);  
  fTPCdEdxP  = new TH2F("TPCdEdXP", "dE/dX; BetaGamma; TPC signal (a.u.)", 1000, 0.1, 10000, 300, 0, 300);        
  
  
}   

AliTPCcalibV0::~AliTPCcalibV0(){
  //
  //
  //
}





void  AliTPCcalibV0::ProcessESD(AliESDEvent *esd, AliStack *stack){
  //
  //
  //
  fESD = esd;
  AliKFParticle::SetField(esd->GetMagneticField());
  MakeV0s();
  if (stack) {
    fStack = stack;
    MakeMC();
  }else{
    fStack =0;
  }
}

void AliTPCcalibV0::MakeMC(){
  //
  // MC comparison 
  //    1. Select interesting particles
  //    2. Assign the recosntructed particles 
  //
  //1. Select interesting particles
  const Float_t kMinP   = 0.2;
  const Float_t kMinPt  = 0.1;
  const Float_t kMaxR   = 0.5;
  const Float_t kMaxTan = 1.2;
  const Float_t kMaxRad = 150;
  //
  if (!fParticles) fParticles = new TObjArray;
  TParticle *part=0;
  //  
  Int_t entries = fStack->GetNtrack();
  for (Int_t ipart=0; ipart<entries; ipart++){
    part = fStack->Particle(ipart);
    if (!part) continue;
    if (part->P()<kMinP) continue;
    if (part->R()>kMaxR) continue;
    if (TMath::Abs(TMath::Tan(part->Theta()-TMath::Pi()*0.5))>kMaxTan) continue;
    Bool_t isInteresting = kFALSE;
    if (part->GetPdgCode()==22) isInteresting =kTRUE;
    if (part->GetPdgCode()==310) isInteresting =kTRUE;
    if (part->GetPdgCode()==111) isInteresting =kTRUE;
    if (TMath::Abs(part->GetPdgCode()==3122)) isInteresting =kTRUE;

    //
    if (!isInteresting) continue;    
    fParticles->AddLast(new TParticle(*part));
  }
  if (fParticles->GetEntries()<1) {
    return;
  }
  //
  //
  //
  Int_t sentries=fParticles->GetEntries();;
  for (Int_t ipart=0; ipart<sentries; ipart++){
    part = (TParticle*)fParticles->At(ipart);
    TParticle *p0 = 0;
    TParticle *p1 = 0;

    Int_t nold =0;
    Int_t nnew =0;
    Int_t id0  = part->GetDaughter(0);
    Int_t id1  = part->GetDaughter(1);    
    if (id0>=fStack->GetNtrack() ) id0*=-1;
    if (id1>=fStack->GetNtrack() ) id1*=-1;
    Bool_t findable = kTRUE;
    if (id0<0 || id1<0) findable = kFALSE;
    Int_t charge =0; 
    if (findable){
      p0 = fStack->Particle(id0);
      if (p0->R()>kMaxRad) findable = kFALSE;
      if (p0->Pt()<kMinPt) findable = kFALSE;
      if (p0->Vz()>250) findable= kFALSE;
      if (TMath::Abs(TMath::Tan(p0->Theta()-TMath::Pi()*0.5))>2) findable=kFALSE;
      if (fPdg->GetParticle(p0->GetPdgCode())==0) findable =kFALSE;
      else
	if (fPdg->GetParticle(p0->GetPdgCode())->Charge()==0) charge++;
	  
      p1 = fStack->Particle(id1);
      if (p1->R()>kMaxRad) findable = kFALSE;
      if (p1->Pt()<kMinPt) findable = kFALSE;
      if (TMath::Abs(p1->Vz())>250) findable= kFALSE;
      if (TMath::Abs(TMath::Tan(p1->Theta()-TMath::Pi()*0.5))>2) findable=kFALSE;
      if (fPdg->GetParticle(p1->GetPdgCode())==0) findable = kFALSE;
      else
	if (fPdg->GetParticle(p1->GetPdgCode())->Charge()==0) charge++;
			  
    }
    //   (*fDebugStream)<<"MC0"<<
    //       "P.="<<part<<
    //       "findable="<<findable<<
    //       "id0="<<id0<<
    //       "id1="<<id1<<
    //       "\n";
    if (!findable) continue;
    Float_t minpt = TMath::Min(p0->Pt(), p1->Pt());
    Int_t type=-1;
    
    //
    // 
    AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
    for (Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++){
      AliESDv0 * v0 = fESD->GetV0(ivertex);
      // select coresponding track
      AliESDtrack * trackN = fESD->GetTrack(v0->GetIndex(0));
      if (TMath::Abs(trackN->GetLabel())!=id0 && TMath::Abs(trackN->GetLabel())!=id1) continue;
      AliESDtrack * trackP = fESD->GetTrack(v0->GetIndex(1));
      if (TMath::Abs(trackP->GetLabel())!=id0 && TMath::Abs(trackP->GetLabel())!=id1) continue;
      TParticle *pn = fStack->Particle(TMath::Abs(trackN->GetLabel()));
      TParticle *pp = fStack->Particle(TMath::Abs(trackP->GetLabel()));
      //
      //
      if ( v0->GetOnFlyStatus()) nnew++;
      if (!v0->GetOnFlyStatus()) nold++;
      if (part->GetPdgCode()==22 && TMath::Abs(pn->GetPdgCode())==11 &&  TMath::Abs(pp->GetPdgCode())==11) 
	type =1;
      if (part->GetPdgCode()==310 && TMath::Abs(pn->GetPdgCode())==211 &&  TMath::Abs(pp->GetPdgCode())==211) 
	type =0;
      if (part->GetPdgCode()==3122){
	if (TMath::Abs(pn->GetPdgCode())==210 ) type=2;
	else type=3;
      }
      AliKFParticle *v0kf       = Fit(primVtx,v0,pn->GetPdgCode(),pp->GetPdgCode());
      v0kf->SetProductionVertex( primVtx );
      AliKFParticle *v0kfc       = Fit(primVtx,v0,pn->GetPdgCode(),pp->GetPdgCode());
      v0kfc->SetProductionVertex( primVtx );
      v0kfc->SetMassConstraint(fPdg->GetParticle(part->GetPdgCode())->Mass());
      Float_t chi2 = v0kf->GetChi2();
      Float_t chi2C = v0kf->GetChi2();
      //
      //
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"MCRC"<<
	  "P.="<<part<<
	  "type="<<type<<
	  "chi2="<<chi2<<
	  "chi2C="<<chi2C<<
	  "minpt="<<minpt<<
	  "id0="<<id0<<
	  "id1="<<id1<<
	  "Pn.="<<pn<<
	  "Pp.="<<pp<<
	  "tn.="<<trackN<<
	  "tp.="<<trackP<<
	  "nold.="<<nold<<
	  "nnew.="<<nnew<<
	  "v0.="<<v0<<
	  "v0kf.="<<v0kf<<
	  "v0kfc.="<<v0kfc<<
	  "\n";     
	delete v0kf;
	delete v0kfc; 
	//
      }
      
      if (cstream){
	(*cstream)<<"MC"<<
	  "P.="<<part<<
	  "charge="<<charge<<
	  "type="<<type<<
	  "minpt="<<minpt<<
	  "id0="<<id0<<
	  "id1="<<id1<<
	  "P0.="<<p0<<
	  "P1.="<<p1<<
	  "nold="<<nold<<
	  "nnew="<<nnew<<
	  "\n";
      }
    }
    fParticles->Delete(); 
  }
}



void AliTPCcalibV0::MakeV0s(){
  //
  //
  //
  const Int_t kMinCluster=40;
  const Float_t kMinR    =0;
  if (! fV0s) fV0s = new TObjArray(10);
  fV0s->Clear();
  //
  // Old V0 finder
  //
  for (Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++){
    AliESDv0 * v0 = fESD->GetV0(ivertex);
    if (v0->GetOnFlyStatus()) continue;
    fV0s->AddLast(v0);
  }
  ProcessV0(0);
  fV0s->Clear(0);
  //
  // MI V0 finder
  //
  for (Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++){
    AliESDv0 * v0 = fESD->GetV0(ivertex);
    if (!v0->GetOnFlyStatus()) continue;
    fV0s->AddLast(v0);
  }
  ProcessV0(1);
  fV0s->Clear();
  return;
  //
  // combinatorial
  //
  Int_t ntracks = fESD->GetNumberOfTracks();
  for (Int_t itrack0=0; itrack0<ntracks; itrack0++){
    AliESDtrack * track0 = fESD->GetTrack(itrack0);
    if (track0->GetSign()>0) continue;
    if ( track0->GetTPCNcls()<kMinCluster) continue;
    if (track0->GetKinkIndex(0)>0) continue;    
    //
    for (Int_t itrack1=0; itrack1<ntracks; itrack1++){
      AliESDtrack * track1 = fESD->GetTrack(itrack1);
      if (track1->GetSign()<0) continue;
      if ( track1->GetTPCNcls()<kMinCluster) continue;
      if (track1->GetKinkIndex(0)>0) continue;
      //
      //      AliExternalTrackParam param0(*track0);
      // AliExternalTrackParam param1(*track1);
      AliV0 vertex;
      vertex.SetParamN(*track0);
      vertex.SetParamP(*track1);
      Float_t xyz[3];
      xyz[0] = fESD->GetPrimaryVertex()->GetXv();
      xyz[1] = fESD->GetPrimaryVertex()->GetYv();
      xyz[2] = fESD->GetPrimaryVertex()->GetZv();
      vertex.Update(xyz);
      if (vertex.GetRr()<kMinR) continue;
      if (vertex.GetDcaV0Daughters()>1.) continue;
      if (vertex.GetDcaV0Daughters()>0.3*vertex.GetRr()) continue;
      // if (vertex.GetPointAngle()<0.9) continue;
      vertex.SetIndex(0,itrack0);
      vertex.SetIndex(1,itrack1);      
      fV0s->AddLast(new AliV0(vertex));
    }
  }
  ProcessV0(2);
  for (Int_t i=0;i<fV0s->GetEntries(); i++) delete fV0s->At(i);
  fV0s->Clear();
}






// void AliTPCcalibV0::ProcessV0(Int_t ftype){
//   //
//   //
//   const Double_t ktimeK0     = 2.684;
//   const Double_t ktimeLambda = 7.89; 
  
  
//   if (! fGammas) fGammas = new TObjArray(10);
//   fGammas->Clear();
//   Int_t nV0s  = fV0s->GetEntries();
//   if (nV0s==0) return;
//   AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
//   //
//   for (Int_t ivertex=0; ivertex<nV0s; ivertex++){
//     AliESDv0 * v0 = (AliESDv0*)fV0s->At(ivertex);
//     AliESDtrack * trackN = fESD->GetTrack(v0->GetIndex(0));
//     AliESDtrack * trackP = fESD->GetTrack(v0->GetIndex(1));
//     // 
//     // 
//     //
//     AliKFParticle *v0K0       = Fit(primVtx,v0,211,211);
//     AliKFParticle *v0Gamma    = Fit(primVtx,v0,11,-11);
//     AliKFParticle *v0Lambda42 = Fit(primVtx,v0,2212,211);
//     AliKFParticle *v0Lambda24 = Fit(primVtx,v0,211,2212);
//     //Set production vertex
//     v0K0->SetProductionVertex( primVtx );
//     v0Gamma->SetProductionVertex( primVtx );
//     v0Lambda42->SetProductionVertex( primVtx );
//     v0Lambda24->SetProductionVertex( primVtx );
//     Double_t massK0, massGamma, massLambda42,massLambda24, massSigma;
//     v0K0->GetMass( massK0,massSigma);
//     v0Gamma->GetMass( massGamma,massSigma);
//     v0Lambda42->GetMass( massLambda42,massSigma);
//     v0Lambda24->GetMass( massLambda24,massSigma);
//     Float_t chi2K0       = v0K0->GetChi2()/v0K0->GetNDF();
//     Float_t chi2Gamma    = v0Gamma->GetChi2()/v0Gamma->GetNDF();
//     Float_t chi2Lambda42 = v0Lambda42->GetChi2()/v0Lambda42->GetNDF();
//     Float_t chi2Lambda24 = v0Lambda24->GetChi2()/v0Lambda24->GetNDF();
//     //
//     // Mass Contrained params
//     //
//     AliKFParticle *v0K0C       = Fit(primVtx,v0,211,211);
//     AliKFParticle *v0GammaC    = Fit(primVtx,v0,11,-11);
//     AliKFParticle *v0Lambda42C = Fit(primVtx,v0,2212,211);
//     AliKFParticle *v0Lambda24C = Fit(primVtx,v0,211,2212);
//     //   
//     v0K0C->SetProductionVertex( primVtx );
//     v0GammaC->SetProductionVertex( primVtx );
//     v0Lambda42C->SetProductionVertex( primVtx );
//     v0Lambda24C->SetProductionVertex( primVtx );

//     v0K0C->SetMassConstraint(fPdg->GetParticle(310)->Mass());
//     v0GammaC->SetMassConstraint(0);
//     v0Lambda42C->SetMassConstraint(fPdg->GetParticle(3122)->Mass());
//     v0Lambda24C->SetMassConstraint(fPdg->GetParticle(-3122)->Mass());
//     //    
//     Double_t timeK0, sigmaTimeK0;  
//     Double_t timeLambda42, sigmaTimeLambda42;  
//     Double_t timeLambda24, sigmaTimeLambda24;  
//     v0K0C->GetLifeTime(timeK0, sigmaTimeK0);
//     //v0K0Gamma->GetLifeTime(timeK0, sigmaTimeK0);
//     v0Lambda42C->GetLifeTime(timeLambda42, sigmaTimeLambda42);
//     v0Lambda24C->GetLifeTime(timeLambda24, sigmaTimeLambda24);
    

//     //
//     Float_t chi2K0C       = v0K0C->GetChi2()/v0K0C->GetNDF();
//     if (chi2K0C<0) chi2K0C=100;
//     Float_t chi2GammaC    = v0GammaC->GetChi2()/v0GammaC->GetNDF();
//     if (chi2GammaC<0) chi2GammaC=100;
//     Float_t chi2Lambda42C = v0Lambda42C->GetChi2()/v0Lambda42C->GetNDF();
//     if (chi2Lambda42C<0) chi2Lambda42C=100;
//     Float_t chi2Lambda24C = v0Lambda24C->GetChi2()/v0Lambda24C->GetNDF();
//     if (chi2Lambda24C<0) chi2Lambda24C=100;
//     //
//     Float_t  minChi2C=99;
//     Int_t   type   =-1;
//     if (chi2K0C<minChi2C) { minChi2C= chi2K0C; type=0;}
//     if (chi2GammaC<minChi2C) { minChi2C= chi2GammaC; type=1;}
//     if (chi2Lambda42C<minChi2C) { minChi2C= chi2Lambda42C; type=2;}
//     if (chi2Lambda24C<minChi2C) { minChi2C= chi2Lambda24C; type=3;}
//     Float_t  minChi2=99;
//     Int_t   type0   =-1;
//     if (chi2K0<minChi2) { minChi2= chi2K0; type0=0;}
//     if (chi2Gamma<minChi2) { minChi2= chi2Gamma; type0=1;}
//     if (chi2Lambda42<minChi2) { minChi2= chi2Lambda42; type0=2;}
//     if (chi2Lambda24<minChi2) { minChi2= chi2Lambda24; type0=3;}
//     Float_t betaGammaP  = trackN->GetP()/fPdg->GetParticle(-2212)->Mass(); 
//     Float_t betaGammaPi = trackN->GetP()/fPdg->GetParticle(-211)->Mass();
//     Float_t betaGammaEl = trackN->GetP()/fPdg->GetParticle(11)->Mass();
//     Float_t dedxTeorP = BetheBlochAleph(betaGammaP);
//     Float_t dedxTeorPi = BetheBlochAleph(betaGammaPi);;
//     Float_t dedxTeorEl = BetheBlochAleph(betaGammaEl);;
//     //
//     //
//     if (minChi2>50) continue;
//     (*fDebugStream)<<"V0"<<
//       "ftype="<<ftype<<
//       "v0.="<<v0<<
//       "trackN.="<<trackN<<
//       "trackP.="<<trackP<<
//       //
//       "dedxTeorP="<<dedxTeorP<<
//       "dedxTeorPi="<<dedxTeorPi<<
//       "dedxTeorEl="<<dedxTeorEl<<
//       //
//       "type="<<type<<
//       "chi2C="<<minChi2C<<
//       "v0K0.="<<v0K0<<
//       "v0Gamma.="<<v0Gamma<<
//       "v0Lambda42.="<<v0Lambda42<<
//       "v0Lambda24.="<<v0Lambda24<<
//       //
//       "chi20K0.="<<chi2K0<<
//       "chi2Gamma.="<<chi2Gamma<<
//       "chi2Lambda42.="<<chi2Lambda42<<
//       "chi2Lambda24.="<<chi2Lambda24<<
//       //
//       "chi20K0c.="<<chi2K0C<<
//       "chi2Gammac.="<<chi2GammaC<<
//       "chi2Lambda42c.="<<chi2Lambda42C<<
//       "chi2Lambda24c.="<<chi2Lambda24C<<
//       //
//       "v0K0C.="<<v0K0C<<
//       "v0GammaC.="<<v0GammaC<<
//       "v0Lambda42C.="<<v0Lambda42C<<
//       "v0Lambda24C.="<<v0Lambda24C<<
//       //
//       "massK0="<<massK0<<
//       "massGamma="<<massGamma<<
//       "massLambda42="<<massLambda42<<
//       "massLambda24="<<massLambda24<<
//       //
//       "timeK0="<<timeK0<<
//       "timeLambda42="<<timeLambda42<<
//       "timeLambda24="<<timeLambda24<<
//       "\n";
//     if (type==1) fGammas->AddLast(v0); 
//     //
//     //
//     //
//     delete v0K0;
//     delete v0Gamma;
//     delete v0Lambda42;
//     delete v0Lambda24;    
//     delete v0K0C;
//     delete v0GammaC;
//     delete v0Lambda42C;
//     delete v0Lambda24C;    
//   }
//   ProcessPI0(); 
// }




void AliTPCcalibV0::ProcessV0(Int_t ftype){
  //
  //
  
  if (! fGammas) fGammas = new TObjArray(10);
  fGammas->Clear();
  Int_t nV0s  = fV0s->GetEntries();
  if (nV0s==0) return;
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  //
  for (Int_t ivertex=0; ivertex<nV0s; ivertex++){
    AliESDv0 * v0 = (AliESDv0*)fV0s->At(ivertex);
    AliESDtrack * trackN = fESD->GetTrack(v0->GetIndex(0)); // negative track
    AliESDtrack * trackP = fESD->GetTrack(v0->GetIndex(1)); // positive track
    
    const AliExternalTrackParam * paramInNeg = trackN->GetInnerParam();
    const AliExternalTrackParam * paramInPos = trackP->GetInnerParam();
  
    if (!paramInPos) continue; // in case the inner paramters do not exist
    if (!paramInNeg) continue;
    // 
    // 
    //
    AliKFParticle *v0K0       = Fit(primVtx,v0,-211,211);
    AliKFParticle *v0Gamma    = Fit(primVtx,v0,11,-11);
    AliKFParticle *v0Lambda42 = Fit(primVtx,v0,-2212,211);
    AliKFParticle *v0Lambda24 = Fit(primVtx,v0,-211,2212);
    //Set production vertex
    v0K0->SetProductionVertex( primVtx );
    v0Gamma->SetProductionVertex( primVtx );
    v0Lambda42->SetProductionVertex( primVtx );
    v0Lambda24->SetProductionVertex( primVtx );
    Double_t massK0, massGamma, massLambda42,massLambda24, massSigma;
    v0K0->GetMass( massK0,massSigma);
    v0Gamma->GetMass( massGamma,massSigma);
    v0Lambda42->GetMass( massLambda42,massSigma);
    v0Lambda24->GetMass( massLambda24,massSigma);
    Float_t chi2K0       = v0K0->GetChi2()/v0K0->GetNDF();
    Float_t chi2Gamma    = v0Gamma->GetChi2()/v0Gamma->GetNDF();
    Float_t chi2Lambda42 = v0Lambda42->GetChi2()/v0Lambda42->GetNDF();
    Float_t chi2Lambda24 = v0Lambda24->GetChi2()/v0Lambda24->GetNDF();
    //
    // Mass Contrained params
    //
    AliKFParticle *v0K0C       = Fit(primVtx,v0,-211,211);
    AliKFParticle *v0GammaC    = Fit(primVtx,v0,11,-11);
    AliKFParticle *v0Lambda42C = Fit(primVtx,v0,-2212,211); //lambdaBar
    AliKFParticle *v0Lambda24C = Fit(primVtx,v0,-211,2212); //lambda
    //   
    v0K0C->SetProductionVertex( primVtx );
    v0GammaC->SetProductionVertex( primVtx );
    v0Lambda42C->SetProductionVertex( primVtx );
    v0Lambda24C->SetProductionVertex( primVtx );

    v0K0C->SetMassConstraint(fPdg->GetParticle(310)->Mass());
    v0GammaC->SetMassConstraint(0);
    v0Lambda42C->SetMassConstraint(fPdg->GetParticle(-3122)->Mass());
    v0Lambda24C->SetMassConstraint(fPdg->GetParticle(3122)->Mass());
    //    
    Double_t timeK0, sigmaTimeK0;  
    Double_t timeLambda42, sigmaTimeLambda42;  
    Double_t timeLambda24, sigmaTimeLambda24;  
    v0K0C->GetLifeTime(timeK0, sigmaTimeK0);
    //v0K0Gamma->GetLifeTime(timeK0, sigmaTimeK0);
    v0Lambda42C->GetLifeTime(timeLambda42, sigmaTimeLambda42);
    v0Lambda24C->GetLifeTime(timeLambda24, sigmaTimeLambda24);
    

    //
    Float_t chi2K0C       = v0K0C->GetChi2()/v0K0C->GetNDF();
    if (chi2K0C<0) chi2K0C=100;
    Float_t chi2GammaC    = v0GammaC->GetChi2()/v0GammaC->GetNDF();
    if (chi2GammaC<0) chi2GammaC=100;
    Float_t chi2Lambda42C = v0Lambda42C->GetChi2()/v0Lambda42C->GetNDF();
    if (chi2Lambda42C<0) chi2Lambda42C=100;
    Float_t chi2Lambda24C = v0Lambda24C->GetChi2()/v0Lambda24C->GetNDF();
    if (chi2Lambda24C<0) chi2Lambda24C=100;
    //
    Float_t  minChi2C=99;
    Int_t   type   =-1;
    if (chi2K0C<minChi2C) { minChi2C= chi2K0C; type=0;}
    if (chi2GammaC<minChi2C) { minChi2C= chi2GammaC; type=1;}
    if (chi2Lambda42C<minChi2C) { minChi2C= chi2Lambda42C; type=2;}
    if (chi2Lambda24C<minChi2C) { minChi2C= chi2Lambda24C; type=3;}
    Float_t  minChi2=99;
    Int_t   type0   =-1;
    if (chi2K0<minChi2) { minChi2= chi2K0; type0=0;}
    if (chi2Gamma<minChi2) { minChi2= chi2Gamma; type0=1;}
    if (chi2Lambda42<minChi2) { minChi2= chi2Lambda42; type0=2;}
    if (chi2Lambda24<minChi2) { minChi2= chi2Lambda24; type0=3;}
    
     // 0 is  negative particle; 1 is positive particle
    Float_t betaGamma0 = 0;
    Float_t betaGamma1 = 0;
    
    switch (type) {
     case 0:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(-211)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(211)->Mass();
      break;
     case 1:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(11)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(-11)->Mass();
      break;
     case 2:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(-2212)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(211)->Mass();
      break;
     case 3:
      betaGamma0 = paramInNeg->GetP()/fPdg->GetParticle(-211)->Mass();
      betaGamma1 = paramInPos->GetP()/fPdg->GetParticle(2212)->Mass();
      break;
    }
 
    // cuts and histogram filling
    Int_t numCand = 0; // number of particle types which have a chi2 < 10*minChi2
        
    if (minChi2C < 2 && ftype == 1) {
     //
     if (chi2K0C < 10*minChi2C) numCand++;
     if (chi2GammaC < 10*minChi2C) numCand++;
     if (chi2Lambda42C < 10*minChi2C) numCand++;
     if (chi2Lambda24C < 10*minChi2C) numCand++;
     //
     if (numCand < 2) {
      if (paramInNeg->GetP() > 0.4) fTPCdEdx->Fill(betaGamma0, trackN->GetTPCsignal());
      if (paramInPos->GetP() > 0.4) fTPCdEdx->Fill(betaGamma1, trackP->GetTPCsignal());
     }    
    }
    
    //
    //
    // write output tree
    if (minChi2>50) continue;
    TTreeSRedirector *cstream = GetDebugStreamer();
    if (cstream){
      (*cstream)<<"V0"<<
	"ftype="<<ftype<<
	"v0.="<<v0<<
	"trackN.="<<trackN<<
	"trackP.="<<trackP<<
	//
	"betaGamma0="<<betaGamma0<<
	"betaGamma1="<<betaGamma1<<
	//
	"type="<<type<<
	"chi2C="<<minChi2C<<
	"v0K0.="<<v0K0<<
	"v0Gamma.="<<v0Gamma<<
	"v0Lambda42.="<<v0Lambda42<<
	"v0Lambda24.="<<v0Lambda24<<
	//
	"chi20K0.="<<chi2K0<<
	"chi2Gamma.="<<chi2Gamma<<
	"chi2Lambda42.="<<chi2Lambda42<<
	"chi2Lambda24.="<<chi2Lambda24<<
	//
	"chi20K0c.="<<chi2K0C<<
	"chi2Gammac.="<<chi2GammaC<<
	"chi2Lambda42c.="<<chi2Lambda42C<<
	"chi2Lambda24c.="<<chi2Lambda24C<<
	//
	"v0K0C.="<<v0K0C<<
	"v0GammaC.="<<v0GammaC<<
	"v0Lambda42C.="<<v0Lambda42C<<
	"v0Lambda24C.="<<v0Lambda24C<<
	//
	"massK0="<<massK0<<
	"massGamma="<<massGamma<<
	"massLambda42="<<massLambda42<<
	"massLambda24="<<massLambda24<<
	//
	"timeK0="<<timeK0<<
	"timeLambda42="<<timeLambda42<<
	"timeLambda24="<<timeLambda24<<
	"\n";
    }
    if (type==1) fGammas->AddLast(v0); 
    //
    //
    //
    delete v0K0;
    delete v0Gamma;
    delete v0Lambda42;
    delete v0Lambda24;    
    delete v0K0C;
    delete v0GammaC;
    delete v0Lambda42C;
    delete v0Lambda24C; 
  }
  ProcessPI0(); 
}



void AliTPCcalibV0::ProcessPI0(){
  //
  //
  //
  Int_t nentries = fGammas->GetEntries();
  if (nentries<2) return;
  // 
  Double_t m0[3], m1[3];
  AliKFVertex primVtx(*(fESD->GetPrimaryVertex()));
  for (Int_t i0=0; i0<nentries; i0++){
    AliESDv0 *v00 = (AliESDv0*)fGammas->At(i0); 
    v00->GetPxPyPz (m0[0], m0[1], m0[2]);
    AliKFParticle *p00 = Fit(primVtx, v00, 11,-11);
    p00->SetProductionVertex( primVtx );
    p00->SetMassConstraint(0);
    //
    for (Int_t i1=i0; i1<nentries; i1++){
      AliESDv0 *v01 = (AliESDv0*)fGammas->At(i1);
      v01->GetPxPyPz (m1[0], m1[1], m1[2]);
      AliKFParticle *p01 = Fit(primVtx, v01, 11,-11);
      p01->SetProductionVertex( primVtx );
      p01->SetMassConstraint(0);
      if (v00->GetIndex(0) != v01->GetIndex(0) && 
	  v00->GetIndex(1) != v01->GetIndex(1)){
	AliKFParticle pi0( *p00,*p01); 
	pi0.SetProductionVertex(primVtx);
	Double_t n1 = TMath::Sqrt (m0[0]*m0[0] + m0[1]*m0[1] + m0[2]*m0[2]);
        Double_t n2 = TMath::Sqrt (m1[0]*m1[0] + m1[1]*m1[1] + m1[2]*m1[2]);
        Double_t mass = TMath::Sqrt(2.*(n1*n2 - (m0[0]*m1[0] + m0[1]*m1[1] + m0[2]*m1[2])));
	TTreeSRedirector *cstream = GetDebugStreamer();
	if (cstream){
	  (*cstream)<<"PI0"<<
	    "v00.="<<v00<<
	    "v01.="<<v01<<
	    "mass="<<mass<<
	    "p00.="<<p00<<
	    "p01.="<<p01<<
	    "pi0.="<<&pi0<<
	    "\n";	
	}
      }
      delete p01;
    }
    delete p00;
  }
}





AliKFParticle * AliTPCcalibV0::Fit(AliKFVertex & /*primVtx*/, AliESDv0 *v0, Int_t PDG1, Int_t PDG2){
  //
  // Make KF Particle
  //
  AliKFParticle p1( *(v0->GetParamN()), PDG1 );
  AliKFParticle p2( *(v0->GetParamP()), PDG2 );
  AliKFParticle *V0 = new AliKFParticle;
  Double_t x, y, z;
  v0->GetXYZ(x,y,z );
  V0->SetVtxGuess(x,y,z);
  *(V0)+=p1;
  *(V0)+=p2;
  return V0;  
}

TH2F * AliTPCcalibV0::GetHistograms() {
  //
  //
  //
 return fTPCdEdx;
}




void AliTPCcalibV0::BinLogX(TH2F *h) {
  //
  //
  //
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Double_t from = axis->GetXmin();
   Double_t to = axis->GetXmax();
   Double_t *new_bins = new Double_t[bins + 1];
   
   new_bins[0] = from;
   Double_t factor = pow(to/from, 1./bins);
  
   for (int i = 1; i <= bins; i++) {
     new_bins[i] = factor * new_bins[i-1];
   }
   axis->Set(bins, new_bins);
   delete new_bins;
   
}




