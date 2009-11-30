/**************************************************************************
 * Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)   *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//----------------------------------------------------------------------------------------------------------------------
//                        class AliResonanceKinkLikeSign
//        Example of an analysis task for producing a like-sign background for resonances having at least one 
//        kaon-kink in their decay products. 
//        Background is calculated from a positive kaon kink and a positive track but other possibilities are feasible.
//----------------------------------------------------------------------------------------------------------------------

#include "AliESDEvent.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "AliESDkink.h"

#include "AliResonanceKinkLikeSign.h"

ClassImp(AliResonanceKinkLikeSign)

//________________________________________________________________________
AliResonanceKinkLikeSign::AliResonanceKinkLikeSign(const char *name) 
  : AliAnalysisTaskSE(name), fDebug(0), fListOfHistos(0), f1(0), f2(0), fPosKaonLikeSign(0), fLikeSignInvmassPt(0), fMaxNSigmaToVertex(0), fMinPtTrackCut(0), fMaxDCAxy(0), fMaxDCAzaxis(0), 
fMinTPCclusters(0),fMaxChi2PerTPCcluster(0), fMaxCov0(0), fMaxCov2(0), fMaxCov5(0) , fMaxCov9(0), fMaxCov14(0), fdaughter1pdg(0), fdaughter2pdg(0), fnbins(0), fnlowx(0), fnhighx(0), floweta(0), fuppereta(0), fminKinkRadius(0), fmaxKinkRadius(0), fminQt(0), fmaxQt(0), fptbins(0), flowpt(0), fupperpt(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliResonanceKinkLikeSign::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  
   f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
   f1->SetParameter(0,0.493677);
   f1->SetParameter(1,0.9127037);
   f1->SetParameter(2,TMath::Pi());

   f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
   f2->SetParameter(0,0.13957018);
   f2->SetParameter(1,0.2731374);
   f2->SetParameter(2,TMath::Pi());
  
   fPosKaonLikeSign=new TH1D("fPosKaonLikeSign"," ", fnbins, fnlowx, fnhighx);
   fLikeSignInvmassPt=new TH2D("fLikeSignInvmassPt"," ", fnbins, fnlowx, fnhighx, fptbins, flowpt, fupperpt);
 
   fListOfHistos=new TList();
   fListOfHistos->Add(fPosKaonLikeSign);
   fListOfHistos->Add(fLikeSignInvmassPt);     

}

//________________________________________________________________________
void AliResonanceKinkLikeSign::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
 
  AliVEvent *event = InputEvent();
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
     Printf("ERROR: Could not retrieve esd");
     return;
  }

   Double_t daughter1Mass, daughter2Mass;
  
   if (fdaughter1pdg==kKPlus)  {
     daughter1Mass=TDatabasePDG::Instance()->GetParticle(fdaughter1pdg)->Mass();
     daughter2Mass=TDatabasePDG::Instance()->GetParticle(fdaughter2pdg)->Mass(); 
   }
  
   if (fdaughter1pdg!=kKPlus)  {
     daughter1Mass=TDatabasePDG::Instance()->GetParticle(fdaughter2pdg)->Mass();
     daughter2Mass=TDatabasePDG::Instance()->GetParticle(fdaughter1pdg)->Mass();   
   }  // to ensure that daughter1Mass has always the kaon mass

  if (!esd) {
    Printf("ERROR: esd not available");
    return;
  }  

  const AliESDVertex* vertex = GetEventVertex(esd);
  if(!vertex) return;
  
  Double_t ptrackpos[3], ptrackneg[3];
  
  TLorentzVector p4pos;
  TLorentzVector p4neg;
  TLorentzVector p4comb; 

  for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
    
    AliESDtrack* trackpos = esd->GetTrack(iTracks);
    if (!trackpos) {
      if (fDebug > 0) Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if (trackpos->GetSign() < 0) continue;
    
    AliExternalTrackParam *tpcTrackpos = (AliExternalTrackParam *)trackpos->GetTPCInnerParam();
    if(!tpcTrackpos) continue;
    ptrackpos[0]=tpcTrackpos->Px();
    ptrackpos[1]=tpcTrackpos->Py();   
    ptrackpos[2]=tpcTrackpos->Pz();  
    
    Bool_t firstLevelAcceptPosTrack=IsAcceptedForKink(esd, vertex, trackpos);
    if(firstLevelAcceptPosTrack==kFALSE) continue;
    
    TVector3 posTrackMom(ptrackpos[0],ptrackpos[1],ptrackpos[2]);
    
    Int_t indexKinkPos=trackpos->GetKinkIndex(0);
    Bool_t posKaonKinkFlag=0;
    if(indexKinkPos<0) posKaonKinkFlag=IsKink(esd, indexKinkPos, posTrackMom);
    if(posKaonKinkFlag==0) continue; 
    if(posKaonKinkFlag==1) p4pos.SetVectM(posTrackMom, daughter1Mass);

    for (Int_t j=0; j<esd->GetNumberOfTracks(); j++) {
        if(iTracks==j) continue;
        AliESDtrack* trackneg=esd->GetTrack(j);
        if (!trackneg) {
          if (fDebug > 0) Printf("ERROR: Could not receive track %d", j);
          continue;
        } 
        if (trackneg->GetSign() < 0) continue;

	AliExternalTrackParam *tpcTrackneg = (AliExternalTrackParam *)trackneg->GetTPCInnerParam();
        if(!tpcTrackneg) continue;
        ptrackneg[0]=tpcTrackneg->Px();
        ptrackneg[1]=tpcTrackneg->Py();   
        ptrackneg[2]=tpcTrackneg->Pz();  
    
        Bool_t firstLevelAcceptNegTrack=IsAcceptedForKink(esd, vertex, trackneg);
        if(firstLevelAcceptNegTrack==kFALSE) continue;	
	
        TVector3 negTrackMom(ptrackneg[0],ptrackneg[1],ptrackneg[2]);

	Bool_t secondLevelAcceptNegTrack=IsAcceptedForTrack(esd, vertex, trackneg);
        if(secondLevelAcceptNegTrack==kFALSE) continue;  
	 	
      	p4neg.SetVectM(negTrackMom, daughter2Mass);  
	
	p4comb=p4pos;
	p4comb+=p4neg;
	
	if(p4comb.Vect().Pt()<=fMinPtTrackCut) continue;	
	
        if((TMath::Abs(p4pos.Vect().Eta())<fuppereta)&&(TMath::Abs(p4neg.Vect().Eta())<fuppereta)&&(p4comb.Vect().Eta()<fuppereta)) {
	
	  fPosKaonLikeSign->Fill(p4comb.M());
	  fLikeSignInvmassPt->Fill(p4comb.M(), p4comb.Vect().Pt());
	
	}
	
      } //inner track loop 

  } //outer track loop 
    
    PostData(1, fListOfHistos);
}      

//________________________________________________________________________
void AliResonanceKinkLikeSign::Terminate(Option_t *) 
{
 
}

//____________________________________________________________________//

Float_t AliResonanceKinkLikeSign::GetSigmaToVertex(AliESDtrack* esdTrack) const {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];

    esdTrack->GetImpactParameters(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    //AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);
  
  if (bRes[0] == 0 || bRes[1] ==0) return -1;
  
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  
  if (TMath::Exp(-d * d / 2) < 1e-10) return 1000;
  
  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  
  return d;
}

//____________________________________________________________________//

const AliESDVertex* AliResonanceKinkLikeSign::GetEventVertex(const AliESDEvent* esd) const

{
  // Get the vertex

  const AliESDVertex* vertex = esd->GetPrimaryVertex();

  if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
  else
  {
     vertex = esd->GetPrimaryVertexSPD();
     if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
     else
     return 0;
  }
}

//________________________________________________________________________

 Bool_t AliResonanceKinkLikeSign::IsAcceptedForKink(AliESDEvent *localesd,
            const AliESDVertex *localvertex, AliESDtrack* localtrack) {
   // Apply the selections for each kink

  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  
  AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)localtrack->GetTPCInnerParam();
  if(!tpcTrack) {
    gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
    dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
    cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
  }
  else {
    gPt = tpcTrack->Pt();
    gPx = tpcTrack->Px();
    gPy = tpcTrack->Py();
    gPz = tpcTrack->Pz();
    tpcTrack->PropagateToDCA(localvertex,
    	       localesd->GetMagneticField(),100.,dca,cov);
  }
  
  if(GetSigmaToVertex(localtrack) > fMaxNSigmaToVertex) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a %lf sigmas to vertex TPC (max. requested: %lf)",   GetSigmaToVertex(localtrack),fMaxNSigmaToVertex);
      return kFALSE;
  }
  
  if(TMath::Abs(dca[0]) > fMaxDCAxy) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a value of dca(xy) (TPC) of %lf (max. requested: %lf)", TMath::Abs(dca[0]), fMaxDCAxy);
      return kFALSE;
  }
    
  if(TMath::Abs(dca[1]) > fMaxDCAzaxis) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a value of dca(z) of %lf (max. requested: %lf)", TMath::Abs(dca[1]), fMaxDCAzaxis);
      return kFALSE;
  }
  
  if(gPt <= fMinPtTrackCut) {
      if (fDebug > 1) Printf("IsAcceptedKink: Track rejected because it has a min value of pt of %lf (min. requested: %lf)", gPt, fMinPtTrackCut);
      return kFALSE;
  } 
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliResonanceKinkLikeSign::IsAcceptedForTrack(AliESDEvent *localesd,                                                                                                                                          const AliESDVertex *localvertex, AliESDtrack *localtrack) {
   // Apply the selections for each track

  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  
  AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)localtrack->GetTPCInnerParam();
  if(!tpcTrack) {
    gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
    dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
    cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
  }
  else {
    gPt = tpcTrack->Pt();
    gPx = tpcTrack->Px();
    gPy = tpcTrack->Py();
    gPz = tpcTrack->Pz();
    tpcTrack->PropagateToDCA(localvertex,
    	       localesd->GetMagneticField(),100.,dca,cov);
  }
  
  Int_t fcls[200];
  Int_t nClustersTPC=localtrack->GetTPCclusters(fcls);
  Float_t chi2perTPCcluster=-1.0;
  if(nClustersTPC!=0) chi2perTPCcluster=(localtrack->GetTPCchi2())/Float_t(nClustersTPC);
  
  Double_t extCov[15];
  localtrack->GetExternalCovariance(extCov);
  
  if((localtrack->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because of no refited in TPC");
      return kFALSE;
  } 

  if(nClustersTPC < fMinTPCclusters) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of nclusters (TPC) of %ld (min. requested: %ld)", nClustersTPC, fMinTPCclusters);
      return kFALSE;
  } 
  
  if(chi2perTPCcluster > fMaxChi2PerTPCcluster) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of chi2perTPCcluster of %lf (max. requested: %lf)", chi2perTPCcluster, fMaxChi2PerTPCcluster);
      return kFALSE;
  } 

  if(extCov[0] > fMaxCov0) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[0] of %lf (max. requested: %lf)", cov[0], fMaxCov0);
      return kFALSE;
  }
  
  if(extCov[2] > fMaxCov2) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[2] of %lf (max. requested: %lf)", cov[2], fMaxCov2);
      return kFALSE;
  }
    
  if(extCov[5] > fMaxCov5) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[5] of %lf (max. requested: %lf)", cov[5], fMaxCov5);
      return kFALSE;
  }  
  
  if(extCov[9] > fMaxCov9) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[9] of %lf (max. requested: %lf)", cov[9], fMaxCov9);
      return kFALSE;
  }  
  
  if(extCov[14] > fMaxCov14) {
      if (fDebug > 1) Printf("IsAccepted: Track rejected because it has a value of cov[14] of %lf (max. requested: %lf)", cov[14], fMaxCov14);
      return kFALSE;
  } 
 
  return kTRUE;

}

//________________________________________________________________________
Bool_t AliResonanceKinkLikeSign::IsKink(AliESDEvent *localesd, Int_t kinkIndex, TVector3 trackMom) 
{
   // Test some kinematical criteria for each kink

	 AliESDkink *kink=localesd->GetKink(TMath::Abs(kinkIndex)-1);
	 const TVector3 motherMfromKink(kink->GetMotherP());
	 const TVector3 daughterMKink(kink->GetDaughterP());
	 Float_t qt=kink->GetQt();

	 Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
	 Double_t maxDecAngpimu=f2->Eval(motherMfromKink.Mag(),0.,0.,0.);

         Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);
	 
	 Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
	 Float_t p1XM= motherMfromKink.Px();
         Float_t p1YM= motherMfromKink.Py();
         Float_t p1ZM= motherMfromKink.Pz();
         Float_t p2XM= daughterMKink.Px();
         Float_t p2YM= daughterMKink.Py();
         Float_t p2ZM= daughterMKink.Pz();
         Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
         Double_t invariantMassKmu= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());

         if((kinkAngle>maxDecAngpimu)&&(qt>fminQt)&&(qt<fmaxQt)&&((kink->GetR()>fminKinkRadius)&&(kink->GetR()<fmaxKinkRadius))&&(TMath::Abs(trackMom.Eta())<fuppereta)&&(invariantMassKmu<0.6)) {

           if(trackMom.Mag()<=1.1) {
		return kTRUE;
           }
	   else 
	   if (kinkAngle<maxDecAngKmu) {
		return kTRUE;
	   }
	 }
	 return kFALSE;
}
