/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
//
//
//             Base class for Heavy Flavour Correlations Analysis
//             Single Event and Mixed Event Analysis are implemented
//
//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------

/* $Id$ */

#include <TParticle.h>
#include <TVector3.h>
#include <TChain.h>
#include "TROOT.h"
#include "AliHFCorrelator.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliEventPoolManager.h"
#include "AliReducedParticle.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"

using std::cout;
using std::endl;

//_____________________________________________________
AliHFCorrelator::AliHFCorrelator() :
//
// default constructor
//
TNamed(),
fPoolMgr(0x0),         
fPool(0x0),
fhadcuts(0x0),
fAODEvent(0x0),
fAssociatedTracks(0x0),
fmcArray(0x0),
fReducedPart(0x0),
fD0cand(0x0), 
fhypD0(0), 
fDCharge(0),

fmixing(kFALSE),
fmontecarlo(kFALSE),
fsystem(kFALSE),
fUseReco(kTRUE),
fselect(kUndefined),

fUseImpactParameter(0),
fPIDmode(0),

fNofTracks(0),
fPoolContent(0),

fPhiMin(0),
fPhiMax(0),

fPtTrigger(0),
fPhiTrigger(0),
fEtaTrigger(0),


fDeltaPhi(0),
fDeltaEta(0),
fk0InvMass(0)

{
	// default constructor	
}



//_____________________________________________________
AliHFCorrelator::AliHFCorrelator(const Char_t* name, AliHFAssociatedTrackCuts *cuts, Bool_t ppOrPbPb) :
TNamed(name,"title"),
fPoolMgr(0x0),         
fPool(0x0),
fhadcuts(0x0),
fAODEvent(0x0),
fAssociatedTracks(0x0),
fmcArray(0x0),
fReducedPart(0x0),
fD0cand(0x0), 
fhypD0(0),
fDCharge(0),

fmixing(kFALSE),
fmontecarlo(kFALSE),
fsystem(ppOrPbPb),
fUseReco(kTRUE),
fselect(kUndefined),
fUseImpactParameter(0),
fPIDmode(0),

fNofTracks(0),
fPoolContent(0),

fPhiMin(0),
fPhiMax(0),

fPtTrigger(0),
fPhiTrigger(0),
fEtaTrigger(0),


fDeltaPhi(0),
fDeltaEta(0),
fk0InvMass(0)
{
	fhadcuts = cuts; 
}

//_____________________________________________________
AliHFCorrelator::~AliHFCorrelator() 
{
//
// destructor
//	
	
	if(fPoolMgr)  {delete fPoolMgr; fPoolMgr=0;}       
	if(fPool) {delete fPool; fPool=0;}
	if(fhadcuts) {delete fhadcuts; fhadcuts=0;}
	if(fAODEvent) {delete fAODEvent; fAODEvent=0;}
	if(fAssociatedTracks) {delete fAssociatedTracks; fAssociatedTracks=0;}
	if(fmcArray) {delete fmcArray; fmcArray=0;}
	if(fReducedPart) {delete fReducedPart; fReducedPart=0;}
	if(fD0cand) {delete fD0cand; fD0cand=0;}
	
	
	if(fNofTracks) fNofTracks = 0;
	
	if(fPhiMin) fPhiMin = 0;
	if(fPhiMax) fPhiMax = 0;
	
	if(fPtTrigger) fPtTrigger=0;
	if(fPhiTrigger) fPhiTrigger=0;
	if(fEtaTrigger) fEtaTrigger=0;
	
	if(fDeltaPhi) fDeltaPhi=0;
	if(fDeltaEta) fDeltaEta=0;
	
	if(fk0InvMass) fk0InvMass=0;
	
}
//_____________________________________________________
Bool_t AliHFCorrelator::DefineEventPool(){
	// definition of the Pool Manager for Event Mixing
	

	Int_t MaxNofEvents = fhadcuts->GetMaxNEventsInPool();
	Int_t MinNofTracks = fhadcuts->GetMinNTracksInPool();
	Int_t NofCentBins = fhadcuts->GetNCentPoolBins();
	Double_t * CentBins = fhadcuts->GetCentPoolBins();
	Int_t NofZVrtxBins = fhadcuts->GetNZvtxPoolBins();
	Double_t *ZVrtxBins = fhadcuts->GetZvtxPoolBins();
		
			
	fPoolMgr = new AliEventPoolManager(MaxNofEvents, MinNofTracks, NofCentBins, CentBins, NofZVrtxBins, ZVrtxBins);
	if(!fPoolMgr) return kFALSE;
	return kTRUE;
}
//_____________________________________________________
Bool_t AliHFCorrelator::Initialize(){
	
    //  std::cout << "AliHFCorrelator::Initialize"<< std::endl;
  AliInfo("AliHFCorrelator::Initialize") ;
  if(!fAODEvent){
    AliInfo("No AOD event") ;
    return kFALSE;
  }
    //std::cout << "No AOD event" << std::endl;
	
	AliCentrality *centralityObj = 0;
	Int_t multiplicity = -1;
	Double_t MultipOrCent = -1;
	
	// initialize the pool for event mixing
	if(!fsystem){ // pp
	multiplicity = fAODEvent->GetNTracks();
		MultipOrCent = multiplicity; // convert from Int_t to Double_t
	}
	if(fsystem){ // PbPb
		
		centralityObj = fAODEvent->GetHeader()->GetCentralityP();
		MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
		AliInfo(Form("Centrality is %f", MultipOrCent));
	}
	
	AliAODVertex *vtx = fAODEvent->GetPrimaryVertex();
	Double_t zvertex = vtx->GetZ(); // zvertex
	Double_t * CentBins = fhadcuts->GetCentPoolBins();
	Double_t poolmin=CentBins[0];
	Double_t poolmax=CentBins[fhadcuts->GetNCentPoolBins()];

	
		if(TMath::Abs(zvertex)>=10 || MultipOrCent>poolmax || MultipOrCent < poolmin) {
		if(!fsystem)AliInfo(Form("pp Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",zvertex,MultipOrCent));
		if(fsystem) AliInfo(Form("PbPb Event with Zvertex = %.2f cm and centrality = %.1f  out of pool bounds, SKIPPING",zvertex,MultipOrCent));

			return kFALSE;
		}
	
	fPool = fPoolMgr->GetEventPool(MultipOrCent, zvertex);
	
	if (!fPool){
		AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipOrCent, zvertex));
	    return kFALSE;
	}
	//fPool->PrintInfo();
	return kTRUE;
}

//_____________________________________________________
Bool_t AliHFCorrelator::ProcessEventPool(){
	 // analysis on Mixed Events
	//cout << "AliHFCorrelator::ProcessEventPool"<< endl;
		if(!fmixing) return kFALSE;
		if(!fPool->IsReady()) return kFALSE;
		if(fPool->GetCurrentNEvents()<fhadcuts->GetMinEventsToMix()) return kFALSE;
	//	fPool->PrintInfo();
		fPoolContent = fPool->GetCurrentNEvents();
		
		return kTRUE;
	
}

//_____________________________________________________
Bool_t AliHFCorrelator::ProcessAssociatedTracks(Int_t EventLoopIndex, const TObjArray* associatedTracks){
  // TODO: memory leak needs to be fixed, for every call, a new array
  // is allocated, but the pointer immediately lost. The cleanup is
  // not straightforward as in the case of event mixing the pointer
  // will be an external array which must not be deleted.
	fAssociatedTracks = new TObjArray();

	if(!fmixing){ // analysis on Single Event
		
		
		
		if(fselect==kHadron || fselect ==kKaon)	fAssociatedTracks = AcceptAndReduceTracks(fAODEvent);
		if(fselect==kKZero) {fAssociatedTracks = AcceptAndReduceKZero(fAODEvent);}	
		if(fselect==kElectron && associatedTracks) fAssociatedTracks=new TObjArray(*associatedTracks);
		
	}
	
	if(fmixing) { // analysis on Mixed Events
		
			
		fAssociatedTracks = fPool->GetEvent(EventLoopIndex);
				
				
			

	} // end if mixing
	
	if(!fAssociatedTracks) return kFALSE;
	
	fNofTracks = fAssociatedTracks->GetEntriesFast(); 
		
	return kTRUE;
	
}
//_____________________________________________________
Bool_t AliHFCorrelator::Correlate(Int_t loopindex){

	if(loopindex >= fNofTracks) return kFALSE;
	if(!fAssociatedTracks) return kFALSE;
	
	fReducedPart = (AliReducedParticle*)fAssociatedTracks->At(loopindex);
	

	fDeltaPhi = SetCorrectPhiRange(fPhiTrigger - fReducedPart->Phi());
	
	fDeltaEta = fEtaTrigger - fReducedPart->Eta();

	return kTRUE;
	
}
		
//_____________________________________________________
Bool_t AliHFCorrelator::PoolUpdate(const TObjArray* associatedTracks){

	if(!fmixing) return kFALSE;
	if(!fPool) return kFALSE;
	if(fmixing) { // update the pool for Event Mixing
		TObjArray* objArr = NULL;
		if(fselect==kHadron || fselect==kKaon) objArr = (TObjArray*)AcceptAndReduceTracks(fAODEvent);
		else if(fselect==kKZero) objArr = (TObjArray*)AcceptAndReduceKZero(fAODEvent);
		else if(fselect==kElectron && associatedTracks){
		  objArr = new TObjArray(*associatedTracks);
		}
		else return kFALSE;
		if(objArr->GetEntriesFast()>0) fPool->UpdatePool(objArr); // updating the pool only if there are entries in the array
	}
		
	return kTRUE;
	
}
		
//_____________________________________________________
Double_t AliHFCorrelator::SetCorrectPhiRange(Double_t phi){
	Double_t pi = TMath::Pi();
	
	if(phi<fPhiMin) phi = phi + 2*pi;
	if(phi>fPhiMax) phi = phi - 2*pi;
	
	return phi;
}

//_____________________________________________________
TObjArray*  AliHFCorrelator::AcceptAndReduceTracks(AliAODEvent* inputEvent){

  Double_t weight=1.;
  Int_t nTracks = inputEvent->GetNTracks();
  AliAODVertex * vtx = inputEvent->GetPrimaryVertex();
  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  
  Double_t Bz = inputEvent->GetMagneticField();
	
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  //*******************************************************
  // use reconstruction
  if(fUseReco){
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {
      AliAODTrack* track = inputEvent->GetTrack(iTrack);
      if (!track) continue;
      if(!fhadcuts->IsHadronSelected(track,&vESD,Bz)) continue; // apply ESD level selections
      if(!fhadcuts->Charge(fDCharge,track)) continue; // apply selection on charge, if required

      Double_t pT = track->Pt();
      
      //compute impact parameter
      Double_t d0z0[2],covd0z0[3];
      Double_t d0=-999999.;
      if(fUseImpactParameter) track->PropagateToDCA(vtx,Bz,100,d0z0,covd0z0);
      else d0z0[0] = 1. ; // random number - be careful with the cuts you applied
      
      if(fUseImpactParameter==1) d0 = TMath::Abs(d0z0[0]); // use impact parameter
      if(fUseImpactParameter==2) { // use impact parameter over resolution
	if(TMath::Abs(covd0z0[0])>0.00000001) d0 = TMath::Abs(d0z0[0])/TMath::Sqrt(covd0z0[0]); 
	else d0 = -1.; // if the resoultion is Zero, rejects the track - to be on the safe side
	
      }
      
      if(fmontecarlo) {// THIS TO BE CHECKED
        Int_t hadLabel = track->GetLabel();
	if(hadLabel < 0) continue;	
      }
      
      if(!fhadcuts->CheckHadronKinematic(pT,d0)) continue; // apply kinematic cuts
      Bool_t rejectsoftpi = kTRUE;// TO BE CHECKED: DO WE WANT IT TO kTRUE AS A DEFAULT?
      if(fD0cand && !fmixing) rejectsoftpi = fhadcuts->InvMassDstarRejection(fD0cand,track,fhypD0); // TO BE CHECKED: WHY NOT FOR EM?
      
      
      if(fselect ==kKaon){	
	if(!fhadcuts->CheckKaonCompatibility(track,fmontecarlo,fmcArray,fPIDmode)) continue; // check if it is a Kaon - data and MC
      }
      weight=fhadcuts->GetTrackWeight(pT,track->Eta(),pos[2]);
      tracksClone->Add(new AliReducedParticle(track->Eta(), track->Phi(), pT,track->GetLabel(),track->GetID(),d0,rejectsoftpi,track->Charge(),weight));
    } // end loop on tracks
  } // end if use reconstruction kTRUE
  
  //*******************************************************
  
  //use MC truth
  if(fmontecarlo && !fUseReco){
    
    for (Int_t iPart=0; iPart<fmcArray->GetEntriesFast(); iPart++) { 
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iPart));
      if (!mcPart) {
	AliWarning("MC Particle not found in tree, skipping"); 
	continue;
      }
      if(!mcPart->Charge()) continue; // consider only charged tracks
      
      Int_t PDG =TMath::Abs(mcPart->PdgCode()); 
if(fselect ==kHadron) {if(!((PDG==321)||(PDG==211)||(PDG==2212)||(PDG==13)||(PDG==11))) continue;} // select only if kaon, pion, proton, muon or electron
      else if(fselect ==kKaon) {if(!(PDG==321)) continue;} // select only if kaon
      else if(fselect ==kElectron) {if(!(PDG==11)) continue;} // select only if electron

      Double_t pT = mcPart->Pt();
      Double_t d0 =1; // set 1 fot the moment - no displacement calculation implemented yet
      if(!fhadcuts->CheckHadronKinematic(pT,d0)) continue; // apply kinematic cuts
      
      tracksClone->Add(new AliReducedParticle(mcPart->Eta(), mcPart->Phi(), pT,iPart,-1,d0,kFALSE,mcPart->Charge()));
    }
    
  } // end if use  MC truth
  
  
  return tracksClone;
}

//_____________________________________________________
TObjArray*  AliHFCorrelator::AcceptAndReduceKZero(AliAODEvent* inputEvent){
	
	Int_t nOfVZeros = inputEvent->GetNumberOfV0s();
	TObjArray* KZeroClone = new TObjArray;
	AliAODVertex *vertex1 = (AliAODVertex*)inputEvent->GetPrimaryVertex();

 // use reconstruction	 	
     if(fUseReco){
	Int_t v0label = -1;
	Int_t pdgDgK0toPipi[2] = {211,211};
	Double_t mPDGK0=0.497614;//TDatabasePDG::Instance()->GetParticle(310)->Mass();
	const Int_t size = inputEvent->GetNumberOfV0s();
	Int_t prevnegID[size];
	Int_t prevposID[size];
	for(Int_t iv0 =0; iv0< nOfVZeros; iv0++){// loop on all v0 candidates
		AliAODv0 *v0 = (static_cast<AliAODEvent*> (inputEvent))->GetV0(iv0);
		if(!v0) {
		  AliInfo(Form("is not a v0 at step, %d", iv0)) ;
		  //cout << "is not a v0 at step " << iv0 << endl;
		  continue;
		}
		
		if(!fhadcuts->IsKZeroSelected(v0,vertex1)) continue; // check if it is a k0
	    
		// checks to avoid double counting
		Int_t negID = -999;
		Int_t posID = -998;
		//Int_t a = 0;
		prevnegID[iv0] = -997;
		prevposID[iv0]= -996;
		negID = v0->GetNegID();
		posID = v0->GetPosID();
		Bool_t DoubleCounting = kFALSE;
		
		for(Int_t k=0; k<iv0; k++){
			if(((negID==prevnegID[k])&&(posID==prevposID[k]))||((negID==prevposID[k])&&(posID==prevnegID[k]))){
				DoubleCounting = kTRUE;
				//a=k;
				break;
			}//end if
		} // end for
		
		if(DoubleCounting) continue;
		else{ 
			prevposID[iv0] = posID; 
			prevnegID[iv0] = negID;
		}
		
		if(fmontecarlo)	v0label = v0->MatchToMC(310,fmcArray, 0, pdgDgK0toPipi); //match a K0 short
		Double_t k0pt = v0->Pt();
		Double_t k0eta = v0->Eta();
		Double_t k0Phi = v0->Phi();
	    fk0InvMass = v0->MassK0Short();	
		
		//if (loopindex == 0) {
		//	if(!plotassociation) ((TH2F*)fOutput->FindObject("KZeroSpectra"))->Fill(k0InvMass,k0pt); // spectra for all k0
		//	if(plotassociation) ((TH2F*)fOutput->FindObject("KZeroSpectraifHF"))->Fill(k0InvMass,k0pt); // spectra for k0 in association with a D*
		//}
		// if there are more D* candidates per event, loopindex == 0 makes sure you fill the mass spectra only once!
		
		if(TMath::Abs(fk0InvMass-mPDGK0)>3*0.004) continue; // select candidates within 3 sigma
		KZeroClone->Add(new AliReducedParticle(k0eta,k0Phi,k0pt,v0label));
		
	}
     } // end if use reconstruction kTRUE
	


//*********************************************************************//
     //use MC truth
     if(fmontecarlo && !fUseReco){
		
		for (Int_t iPart=0; iPart<fmcArray->GetEntriesFast(); iPart++) { 
			AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iPart));
			if (!mcPart) {
				AliWarning("MC Particle not found in tree, skipping"); 
				continue;
			}
			
			Int_t PDG =TMath::Abs(mcPart->PdgCode()); 
			if(!(PDG==310)) continue; // select only if k0 short
		
			Double_t pT = mcPart->Pt();
            Double_t d0 =1; // set 1 fot the moment - no displacement calculation implemented yet
			if(!fhadcuts->CheckHadronKinematic(pT,d0)) continue; // apply kinematic cuts
			
			KZeroClone->Add(new AliReducedParticle(mcPart->Eta(), mcPart->Phi(), pT,iPart,-1,d0,kTRUE,mcPart->Charge()));
		}
		
	} // end if use  MC truth



	return KZeroClone;
}
