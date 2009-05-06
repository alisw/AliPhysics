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

//----------------------------------------------------------------------------------------------------------------
//                        class AliResonanceKinkLikeSign
//        Example of an analysis task for producing a like-sign background for resonances having at least one 
//        kaon-kink in their decay products. 
//        Background is calculated from a negative kaon kink and a negative track.
//-----------------------------------------------------------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDInputHandler.h"

#include "AliResonanceKinkLikeSign.h"
#include "AliESDkink.h"

ClassImp(AliResonanceKinkLikeSign)

//________________________________________________________________________
AliResonanceKinkLikeSign::AliResonanceKinkLikeSign() 
  : AliAnalysisTaskSE(), fESD(0), fListOfHistos(0), f1(0), f2(0), fNegKaonLikeSign(0)
{
  // Constructor
}

//________________________________________________________________________
AliResonanceKinkLikeSign::AliResonanceKinkLikeSign(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListOfHistos(0), f1(0), f2(0), fNegKaonLikeSign(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliResonanceKinkLikeSign::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches, we want to process only MC
    tree->SetBranchStatus("*",kTRUE );
    tree->SetBranchStatus("*Calo*", kFALSE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliResonanceKinkLikeSign::CreateOutputObjects() 
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
  
   //OpenFile(1);  //uncomment for proof analysis 
   
    // for K*0(892)
   fNegKaonLikeSign=new TH1D("fNegKaonLikeSign"," ", 60,0.6,1.2);
   
   // for phi(1020)
 //  fNegKaonLikeSign=new TH1D("fNegKaonLikeSign"," ", 70,0.99,1.088);
    

   fListOfHistos=new TList();
   fListOfHistos->Add(fNegKaonLikeSign);  

}

//________________________________________________________________________
void AliResonanceKinkLikeSign::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
 
   Float_t daughter1Mass=AliPID::ParticleMass(AliPID::kKaon);
   Float_t daughter2Mass=AliPID::ParticleMass(AliPID::kPion);

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }  

  const AliESDVertex* vertex = GetEventVertex(fESD);
  if(!vertex) return;
  
  Double_t ptrackpos[3], ptrackneg[3];
  
  TLorentzVector p4pos;
  TLorentzVector p4neg;
  TLorentzVector p4comb; 

  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* trackpos = fESD->GetTrack(iTracks);
    if (!trackpos) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    
    if (trackpos->GetSign() > 0) continue;
    
    trackpos->GetPxPyPz(ptrackpos);
    
    Float_t nSigmaToVertex = GetSigmaToVertex(trackpos);      
    
    Float_t bpos[2];
    Float_t bCovpos[3];
    trackpos->GetImpactParameters(bpos,bCovpos);
    
    if (bCovpos[0]<=0 || bCovpos[2]<=0) {
     Printf("Estimated b resolution lower or equal zero!");
     bCovpos[0]=0; bCovpos[2]=0;
    }

    Float_t dcaToVertexXYpos = bpos[0];
    Float_t dcaToVertexZpos = bpos[1];

    if(nSigmaToVertex>=4) continue;
    if((dcaToVertexXYpos>3.0)||(dcaToVertexZpos>3.0)) continue;
    
    TVector3 posTrackMom(ptrackpos[0],ptrackpos[1],ptrackpos[2]);
  
    if(posTrackMom.Perp()<=0.25) continue; 
	
    //Uncomment the following block if the Like Sign is made of K- kink + negative track
    
    //Int_t indexKink=trackpos->GetKinkIndex(0);
    //Int_t kaonKinkFlag=0;
//     if(indexKink<0){
// 		
//         AliESDkink *kink=fESD->GetKink(TMath::Abs(IndexKink)-1);
// 	const TVector3 motherMfromKink(kink->GetMotherP());
// 	const TVector3 daughterMKink(kink->GetDaughterP());
// 	Float_t qt=kink->GetQt();
// 
//         Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
//         Double_t maxDecAngpimu=f2->Eval(motherMfromKink.Mag(),0.,0.,0.);
// 
//         Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);
// 	 
// 	Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
//         Float_t p1XM= motherMfromKink.Px();
//         Float_t p1YM= motherMfromKink.Py();
//         Float_t p1ZM= motherMfromKink.Pz();
//         Float_t p2XM= daughterMKink.Px();
//         Float_t p2YM= daughterMKink.Py();
//         Float_t p2ZM= daughterMKink.Pz();
//         Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
//         Double_t invariantMassKmu= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());
// 
//         if((kinkAngle>maxDecAngpimu)&&(qt>0.05)&&(qt<0.25)&&((kink->GetR()>110.)&&(kink->GetR()<230.))&&(TMath::Abs(posTrackMom.Eta())<1.1)&&(invariantMassKmu<0.6)) {
// 
//           if(posTrackMom.Mag()<=1.1) {
//  	   kaonKinkFlag=1;
//  	  }
// 	  else 
// 	  if (kinkAngle<maxDecAngKmu) {
// 	   kaonKinkFlag=1;	
// 	  }
// 	}
// 
//     }  //End Kink Information   
//     
//     if(kaonKinkFlag==0) continue; 
//     if(kaonKinkFlag==1) p4pos.SetVectM(posTrackMom,daughter1Mass);

      // Comment the following statements till the "for" if the Like Sign of K- kink + negative track is needed

      UInt_t status=trackpos->GetStatus();
      if((status&AliESDtrack::kTPCrefit)==0) continue;
      if(trackpos->GetTPCclusters(0)<50) continue;
      if((trackpos->GetTPCchi2()/trackpos->GetTPCclusters(0))>3.5) continue;
      Double_t extCovPos[15];
      trackpos->GetExternalCovariance(extCovPos);    
      if(extCovPos[0]>2) continue;
      if(extCovPos[2]>2) continue;    
      if(extCovPos[5]>0.5) continue;  
      if(extCovPos[9]>0.5) continue;
      if(extCovPos[14]>2) continue; 

       p4pos.SetVectM(posTrackMom,daughter1Mass);

      for (Int_t j=0; j<fESD->GetNumberOfTracks(); j++) {
        if(iTracks==j) continue;
        AliESDtrack* trackneg=fESD->GetTrack(j);
        if (trackneg->GetSign() > 0) continue;
	
	trackneg->GetPxPyPz(ptrackneg);
        Float_t negSigmaToVertex = GetSigmaToVertex(trackneg);
      
        Float_t bneg[2];
        Float_t bCovneg[3];
        trackneg->GetImpactParameters(bneg,bCovneg);
        if (bCovneg[0]<=0 || bCovneg[2]<=0) {
          Printf("Estimated b resolution lower or equal zero!");
          bCovneg[0]=0; bCovneg[2]=0;
        }

        Float_t dcaToVertexXYneg = bneg[0];
        Float_t dcaToVertexZneg = bneg[1];
    
        if(negSigmaToVertex>=4) continue;
        if((dcaToVertexXYneg>3.0)||(dcaToVertexZneg>3.0)) continue;

        TVector3 negTrackMom(ptrackneg[0],ptrackneg[1],ptrackneg[2]);

        if(negTrackMom.Perp()<=0.25) continue;	
	
	UInt_t statusneg=trackneg->GetStatus();

        if((statusneg&AliESDtrack::kTPCrefit)==0) continue;

        if(trackneg->GetTPCclusters(0)<50) continue;
        if((trackneg->GetTPCchi2()/trackneg->GetTPCclusters(0))>3.5) continue;
       	Double_t extCovneg[15];
        trackneg->GetExternalCovariance(extCovneg);
        if(extCovneg[0]>2) continue;
        if(extCovneg[2]>2) continue;    
        if(extCovneg[5]>0.5) continue;  
        if(extCovneg[9]>0.5) continue;
        if(extCovneg[14]>2) continue;

	p4neg.SetVectM(negTrackMom, daughter2Mass);
	
	p4comb=p4pos;
	p4comb+=p4neg;
	fNegKaonLikeSign->Fill(p4comb.M());
	
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

