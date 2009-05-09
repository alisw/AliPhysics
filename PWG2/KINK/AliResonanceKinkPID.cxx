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
//                        class AliResonanceKinkPID
//        Example of an analysis task for reconstructing resonances having at least one kaon-kink in their decay 
//        products. It provides basic plots as well as plots helping to calculate the corrections.
//        Usage: To analyse a resonance having a kaon in its decay products, one has to modify the integer 
//        variables resonancePDG, daughter1 and daughter2 accordingly as well as daughter1Mass  and daughter2Mass.
//        Also, depending on the analysis mode (ESD or MC), fAnalysisType in the constructor must also be changed 
//-----------------------------------------------------------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TH2D.h"
#include "TParticle.h"
#include <TVector3.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "TF1.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliResonanceKinkPID.h"
#include "AliESDkink.h"
#include "AliStack.h"

ClassImp(AliResonanceKinkPID)

//________________________________________________________________________
AliResonanceKinkPID::AliResonanceKinkPID() 
  : AliAnalysisTaskSE(), fESD(0), fListOfHistos(0), fOpeningAngle(0), fInvariantMass(0), fInvMassTrue(0), fPhiBothKinks(0), fRecPt(0), fRecEta(0), fRecEtaPt(0), 
  fSimPt(0), fSimEta(0), fSimEtaPt(0), fSimPtKink(0), fSimEtaKink(0),  fSimEtaPtKink(0), 
  fhdr(0), fhdz(0), f1(0), f2(0), fAnalysisType("ESD"), fvtxz(0)

{
  // Constructor
}

//________________________________________________________________________
AliResonanceKinkPID::AliResonanceKinkPID(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListOfHistos(0), fOpeningAngle(0), fInvariantMass(0), fInvMassTrue(0), fPhiBothKinks(0), fRecPt(0), fRecEta(0), fRecEtaPt(0),
  fSimPt(0), fSimEta(0), fSimEtaPt(0), fSimPtKink(0), fSimEtaKink(0),  fSimEtaPtKink(0), 
  fhdr(0), fhdz(0), f1(0), f2(0), fAnalysisType("ESD"), fvtxz(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliResonanceKinkPID::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliResonanceKinkPID::CreateOutputObjects() 
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
  
   // OpenFile(1); Uncomment for proof analysis
   fOpeningAngle=new TH1D("fOpeningAngle"," ", 100,-1.0,1.0);

   // for K*0(892)
   fInvariantMass=new TH1D("fInvariantMass"," ", 60,0.6,1.2);
   fInvMassTrue=new TH1D("fInvMassTrue"," ",60,0.6,1.2);
   fPhiBothKinks=new TH1D("fPhiBothKinks"," ", 60,0.6,1.2);  // Not applicable for K*0

   //for phi(1020)
   //fInvariantMass=new TH1D("fInvariantMass"," ", 70,0.99,1.088);
   //fInvMassTrue=new TH1D("fInvMassTrue"," ", 70,0.99,1.088);
   //fPhiBothKinks=new TH1D("fPhiBothKinks"," ", 70,0.99,1.088);

   fRecPt=new TH1D("fRecPt"," ", 50,0.0,5.0);
   fRecEta=new TH1D("fRecEta"," ", 44,-1.1,1.1);
   fRecEtaPt=new TH2D("fRecEtaPt"," ", 50,0.0,5.0, 44,-1.1,1.1); 
   fSimPt=new TH1D("fSimPt"," ", 50,0.0,5.0);
   fSimEta=new TH1D("fSimEta"," ", 44,-1.1,1.1); 
   fSimEtaPt=new TH2D("fSimEtaPt"," ", 50,0.0,5.0, 44,-1.1,1.1);
   fSimPtKink=new TH1D("fSimPtKink"," ", 50,0.0,5.0);
   fSimEtaKink=new TH1D("fSimEtaKink"," ", 44,-1.1,1.1);
   fSimEtaPtKink=new TH2D("fSimEtaPtKink"," ", 50,0.0,5.0, 44,-1.1,1.1);                
   fhdr=new TH1D("fhdr"," ", 100,0.0,5.0);  
   fhdz=new TH1D("fhdz"," ", 100,0.0,5.0);
   fvtxz=new TH1D("fvtxz"," ", 100,-20.0,20.0);
 
   fListOfHistos=new TList();
  
   fListOfHistos->Add(fOpeningAngle);
   fListOfHistos->Add(fInvariantMass);
   fListOfHistos->Add(fInvMassTrue);
   fListOfHistos->Add(fPhiBothKinks);
   fListOfHistos->Add(fRecPt);    
   fListOfHistos->Add(fRecEta);   
   fListOfHistos->Add(fRecEtaPt);    
   fListOfHistos->Add(fSimPt);    
   fListOfHistos->Add(fSimEta);   
   fListOfHistos->Add(fSimEtaPt);     
   fListOfHistos->Add(fSimPtKink);    
   fListOfHistos->Add(fSimEtaKink);   
   fListOfHistos->Add(fSimEtaPtKink);                                                           
   fListOfHistos->Add(fhdr);
   fListOfHistos->Add(fhdz);
   fListOfHistos->Add(fvtxz);

}

//________________________________________________________________________
void AliResonanceKinkPID::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
 
   enum ResonanceType {kPhi=333, kKstar0=313, kLambda1520=3124};
   enum DaughterType {kdaughterPion=211, kdaughterKaon=321, kdaughterProton=2212};
   
   Int_t resonancePDG=kKstar0;
  
   Int_t daughter1=kdaughterKaon;
   Int_t daughter2=kdaughterPion;
   
   Double_t daughter1Mass=TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();
   Double_t daughter2Mass=TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
   
   AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
   if (!eventHandler) {
    Printf("ERROR: Could not retrieve MC event handler");
    return;
   }

  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }  

  AliStack* stack=mcEvent->Stack();

  if(fAnalysisType == "MC") {
    
  for (Int_t iMc = 0; iMc < stack->GetNprimary(); ++iMc)
  {
    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      Printf("UNEXPECTED: particle with label %d not found in stack (mc loop)", iMc);
      continue;
    }

     if(TMath::Abs(particle->GetPdgCode())==resonancePDG) {
       Int_t firstD=particle->GetFirstDaughter();
       Int_t lastD=particle->GetLastDaughter();
       TParticle *daughterParticle1=stack->Particle(firstD);
       TParticle *daughterParticle2=stack->Particle(lastD);
       
       TParticle* kaonFirstDaughter;
       Int_t mcProcessKaonFirstDaughter = -999;
       
       for(Int_t ia=0; ia<daughterParticle1->GetNDaughters(); ia++){
        if ((daughterParticle1->GetFirstDaughter()+ia)!=-1) {
	  kaonFirstDaughter=stack->Particle(daughterParticle1->GetFirstDaughter()+ia);
          mcProcessKaonFirstDaughter=kaonFirstDaughter->GetUniqueID();
        }
       }       
 
       if((daughterParticle1->Pt()>0.25)&&(daughterParticle2->Pt()>0.25)&&(TMath::Abs(daughterParticle1->Eta())<1.1)&&            (TMath::Abs(daughterParticle2->Eta())<1.1)&&(TMath::Abs(particle->Eta())<1.1)) {
         fSimEta->Fill(particle->Eta());
	 fSimPt->Fill(particle->Pt());
	 fSimEtaPt->Fill(particle->Pt(), particle->Eta());
	 if(mcProcessKaonFirstDaughter==4) {
	   fSimPtKink->Fill(particle->Pt());
	   fSimEtaKink->Fill(particle->Eta());
	   fSimEtaPtKink->Fill(particle->Pt(), particle->Eta());
	 }
       }
     }
  } 
  
  } // end fAnalysisType==MC
  else 
  
  if(fAnalysisType == "ESD") {
  
  const AliESDVertex* vertex = GetEventVertex(fESD);
  if(!vertex) return;
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  fvtxz->Fill(vtx[2]);
  
  Double_t ptrackpos[3], ptrackneg[3];
  
  TLorentzVector p4pos, anp4pos;
  TLorentzVector p4neg, anp4neg;
  TLorentzVector p4comb, anp4comb;
  
  
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* trackpos = fESD->GetTrack(iTracks);
    if (!trackpos) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if (trackpos->GetSign() < 0) continue;
    
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
    
    fhdr->Fill(dcaToVertexXYpos);
    fhdz->Fill(dcaToVertexZpos);

    if(nSigmaToVertex>=4) continue;
    if((dcaToVertexXYpos>3.0)||(dcaToVertexZpos>3.0)) continue;
    
    TVector3 posTrackMom(ptrackpos[0],ptrackpos[1],ptrackpos[2]);
  
    if(posTrackMom.Perp()<=0.25) continue; 
	
    TParticle * partpos = stack->Particle(TMath::Abs(trackpos->GetLabel()));
    if (!partpos) continue;
    Int_t pdgpos = partpos->GetPdgCode();
    Int_t mumlabelpos=partpos->GetFirstMother();
    mumlabelpos = TMath::Abs(mumlabelpos);
    TParticle * mumpos=stack->Particle(mumlabelpos);
    if (!mumpos) continue;
    Int_t mumpdgpos = mumpos->GetPdgCode();
    
    Int_t indexKinkPos=trackpos->GetKinkIndex(0);
    Int_t kaonKinkFlag=0;
    if(indexKinkPos<0){
		
        AliESDkink *poskink=fESD->GetKink(TMath::Abs(indexKinkPos)-1);
	const TVector3 motherMfromKinkPos(poskink->GetMotherP());
	const TVector3 daughterMKinkPos(poskink->GetDaughterP());
	Float_t posQt=poskink->GetQt();

        Double_t maxDecAngKmuPos=f1->Eval(motherMfromKinkPos.Mag(),0.,0.,0.);
        Double_t maxDecAngpimuPos=f2->Eval(motherMfromKinkPos.Mag(),0.,0.,0.);

        Float_t kinkAnglePos=TMath::RadToDeg()*poskink->GetAngle(2);
	 
	Float_t energyDaughterMu=TMath::Sqrt(daughterMKinkPos.Mag()*daughterMKinkPos.Mag()+0.105658*0.105658);
        Float_t p1XM= motherMfromKinkPos.Px();
        Float_t p1YM= motherMfromKinkPos.Py();
        Float_t p1ZM= motherMfromKinkPos.Pz();
        Float_t p2XM= daughterMKinkPos.Px();
        Float_t p2YM= daughterMKinkPos.Py();
        Float_t p2ZM= daughterMKinkPos.Pz();
        Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
        Double_t invariantMassKmuPos= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKinkPos.Mag()*motherMfromKinkPos.Mag());

        if((kinkAnglePos>maxDecAngpimuPos)&&(posQt>0.05)&&(posQt<0.25)&&((poskink->GetR()>110.)&&(poskink->GetR()<230.))&&(TMath::Abs(posTrackMom.Eta())<1.1)&&(invariantMassKmuPos<0.6)) {

          if(posTrackMom.Mag()<=1.1) {
 	   kaonKinkFlag=1;
 	  }
	  else 
	  if (kinkAnglePos<maxDecAngKmuPos) {
	   kaonKinkFlag=1;	
	  }
	}

    }  //End Kink Information   
    
    if(kaonKinkFlag==1) anp4pos.SetVectM(posTrackMom,daughter1Mass);
    
    if(indexKinkPos==0) {
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
   
      p4pos.SetVectM(posTrackMom, daughter2Mass);
    
    }
	
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
    
        fhdr->Fill(dcaToVertexXYneg);
        fhdz->Fill(dcaToVertexZneg);

        if(negSigmaToVertex>=4) continue;
        if((dcaToVertexXYneg>3.0)||(dcaToVertexZneg>3.0)) continue;

        TVector3 negTrackMom(ptrackneg[0],ptrackneg[1],ptrackneg[2]);

        if(negTrackMom.Perp()<=0.25) continue;	
	
        TParticle * partneg = stack->Particle(TMath::Abs(trackneg->GetLabel()));
        if (!partneg) continue;
        Int_t pdgneg = partneg->GetPdgCode();
	Int_t mumlabelneg=partneg->GetFirstMother();
        mumlabelneg = TMath::Abs(mumlabelneg);
        TParticle * mumneg=stack->Particle(mumlabelneg);
        if (!mumneg) continue;
        Int_t mumpdgneg = mumneg->GetPdgCode();
	
	Int_t indexKinkNeg=trackneg->GetKinkIndex(0);
	Int_t negKaonKinkFlag=0;
	if(indexKinkNeg<0){
		
	 AliESDkink *negkink=fESD->GetKink(TMath::Abs(indexKinkNeg)-1);
	 const TVector3 motherMfromKinkNeg(negkink->GetMotherP());
	 const TVector3 daughterMKinkNeg(negkink->GetDaughterP());
	 Float_t negQt=negkink->GetQt();

	 Double_t maxDecAngKmuNeg=f1->Eval(motherMfromKinkNeg.Mag(),0.,0.,0.);
	 Double_t maxDecAngpimuNeg=f2->Eval(motherMfromKinkNeg.Mag(),0.,0.,0.);

         Float_t kinkAngleNeg=TMath::RadToDeg()*negkink->GetAngle(2);
	 
	 Float_t energyDaughterMuNeg=TMath::Sqrt(daughterMKinkNeg.Mag()*daughterMKinkNeg.Mag()+0.105658*0.105658);
	 Float_t p1XMNeg= motherMfromKinkNeg.Px();
         Float_t p1YMNeg= motherMfromKinkNeg.Py();
         Float_t p1ZMNeg= motherMfromKinkNeg.Pz();
         Float_t p2XMNeg= daughterMKinkNeg.Px();
         Float_t p2YMNeg= daughterMKinkNeg.Py();
         Float_t p2ZMNeg= daughterMKinkNeg.Pz();
         Float_t p3DaughterNeg=TMath::Sqrt(((p1XMNeg-p2XMNeg)*(p1XMNeg-p2XMNeg))+((p1YMNeg-p2YMNeg)*(p1YMNeg-p2YMNeg))+((p1ZMNeg-p2ZMNeg)*(p1ZMNeg-p2ZMNeg)));
         Double_t invariantMassKmuNeg= TMath::Sqrt((energyDaughterMuNeg+p3DaughterNeg)*(energyDaughterMuNeg+p3DaughterNeg)-motherMfromKinkNeg.Mag()*motherMfromKinkNeg.Mag());

         if((kinkAngleNeg>maxDecAngpimuNeg)&&(negQt>0.05)&&(negQt<0.25)&&((negkink->GetR()>110.)&&(negkink->GetR()<230.))&&(TMath::Abs(negTrackMom.Eta())<1.1)&&(invariantMassKmuNeg<0.6)) {

           if(negTrackMom.Mag()<=1.1) {
 	  	negKaonKinkFlag=1;
           }
	   else 
	   if (kinkAngleNeg<maxDecAngKmuNeg) {
	        negKaonKinkFlag=1;	
	   }
	 }

	}  //End Kink Information   
	
	if(negKaonKinkFlag==1) p4neg.SetVectM(negTrackMom,daughter1Mass);
	
	if(indexKinkNeg==0)  {
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

	  anp4neg.SetVectM(negTrackMom, daughter2Mass);
	
        }
	
	Double_t openingAngle=(ptrackpos[0]*ptrackneg[0]+ptrackpos[1]*ptrackneg[1]+ptrackpos[2]*ptrackneg[2])/(posTrackMom.Mag()*negTrackMom.Mag());

	if((kaonKinkFlag==1)&&(negKaonKinkFlag==1)) {
	 p4comb=anp4pos;
	 p4comb+=p4neg;
	 if(openingAngle>0.6) fPhiBothKinks->Fill(p4comb.M());
	}
		
	if(negKaonKinkFlag==1) {
	  p4comb=p4pos;
          p4comb+=p4neg;
	  fInvariantMass->Fill(p4comb.M());
	  if ((mumpdgpos==(-resonancePDG))&&(mumpdgneg==(-resonancePDG))&&(mumlabelpos==mumlabelneg)
          &&(pdgpos==daughter2)&&(pdgneg==(-daughter1))&&(TMath::Abs(trackpos->GetLabel())>=0)&&(TMath::Abs(trackneg->GetLabel())>=0)&&(mumlabelpos>=0)&&(mumlabelneg>=0)) {
            fOpeningAngle->Fill(openingAngle);
            fInvMassTrue->Fill(p4comb.M());
	    if((TMath::Abs(p4pos.Vect().Eta())<1.1)&&(TMath::Abs(p4neg.Vect().Eta())<1.1)&&(p4comb.Vect().Eta()<1.1)) {
	      fRecPt->Fill(p4comb.Vect().Pt());
	      fRecEta->Fill(p4comb.Vect().Eta());
	      fRecEtaPt->Fill(p4comb.Vect().Perp(),p4comb.Vect().Eta());

	    }

           }
	  
	}
	
	if(kaonKinkFlag==1) {
          anp4comb=anp4pos;
          anp4comb+=anp4neg;  
	  fInvariantMass->Fill(anp4comb.M());
	  if ((mumpdgpos==resonancePDG)&&(mumpdgneg==resonancePDG)&&(mumlabelpos==mumlabelneg)
          &&(pdgpos==daughter1)&&(pdgneg==(-daughter2))&&(TMath::Abs(trackpos->GetLabel())>=0)&&(TMath::Abs(trackneg->GetLabel())>=0)&&(mumlabelpos>=0)  &&(mumlabelneg>=0)) {
            fOpeningAngle->Fill(openingAngle);
            fInvMassTrue->Fill(p4comb.M());
            if((TMath::Abs(anp4neg.Vect().Eta())<1.1)&&(TMath::Abs(anp4pos.Vect().Eta())<1.1)&&(anp4comb.Vect().Eta()<1.1)) {	
	     fRecPt->Fill(anp4comb.Vect().Pt());
	     fRecEta->Fill(anp4comb.Vect().Eta());
	     fRecEtaPt->Fill(anp4comb.Vect().Pt(), anp4comb.Vect().Eta());
	   }

         }

	}
	 
      } //inner track loop 

  } //outer track loop 
  
  } // end fAnalysisType == ESD
  
    PostData(1, fListOfHistos);
}      

//________________________________________________________________________
void AliResonanceKinkPID::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

//  fHistPt = dynamic_cast<TH1F*> (GetOutputData(0));
//  if (!fHistPt) {
//    Printf("ERROR: fHistPt not available");
//    return;
//  }
   
//  TCanvas *c1 = new TCanvas("AliResonanceKinkPID","Pt MC",10,10,510,510);
//  c1->cd(1)->SetLogy();
//  fHistPt->DrawCopy("E");
}

//____________________________________________________________________//

Float_t AliResonanceKinkPID::GetSigmaToVertex(AliESDtrack* esdTrack) const {
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

const AliESDVertex* AliResonanceKinkPID::GetEventVertex(const AliESDEvent* esd) const
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

