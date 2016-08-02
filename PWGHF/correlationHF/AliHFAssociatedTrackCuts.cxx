/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// Base class for cuts on Associated tracks for HF Correlation analysis
//
// Author: S.Bjelogrlic (Utrecht) sandro.bjelogrlic@cern.ch
////////////////////////////////////////////////////////////////////////
#include <Riostream.h>
#include "AliHFAssociatedTrackCuts.h"
#include "AliAODPidHF.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TString.h"

using std::cout;
using std::endl;

ClassImp(AliHFAssociatedTrackCuts)

//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts::AliHFAssociatedTrackCuts():
AliAnalysisCuts(),
fESDTrackCuts(0),
fPidObj(0),
  fEffWeights(0),

fTrigEffWeightsvspt(0),
fTrigEffWeightsvsptB(0),
fTrigEffWeights(0),
fTrigEffWeightsB(0),
fPoolMaxNEvents(0), 
fPoolMinNTracks(0), 
fMinEventsToMix(0),
fTargetFracTracks(0),
fNzVtxBins(0), 
fNzVtxBinsDim(0), 
fZvtxBins(0), 
fNCentBins(0), 
fNCentBinsDim(0), 
fCentBins(0),

fNofMCEventType(0),
fMCEventType(0),

fNTrackCuts(0),
fAODTrackCuts(0),
fTrackCutsNames(0),
fNvZeroCuts(0),
fAODvZeroCuts(0),
fvZeroCutsNames(0),
fBit(-1),
fCharge(0),
fDescription("")

{
	//
	//default constructor
	//
	//
	//default constructor
	//
	
}

//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts::AliHFAssociatedTrackCuts(const char* name, const char* title):
AliAnalysisCuts(name,title),
fESDTrackCuts(0),
fPidObj(0),
fEffWeights(0),
fTrigEffWeightsvspt(0),
fTrigEffWeightsvsptB(0),
fTrigEffWeights(0),
fTrigEffWeightsB(0),
fPoolMaxNEvents(0), 
fPoolMinNTracks(0), 
fMinEventsToMix(0),
fTargetFracTracks(0),
fNzVtxBins(0), 
fNzVtxBinsDim(0), 
fZvtxBins(0), 
fNCentBins(0), 
fNCentBinsDim(0), 
fCentBins(0),

fNofMCEventType(0),
fMCEventType(0),

fNTrackCuts(0),
fAODTrackCuts(0),
fTrackCutsNames(0),
fNvZeroCuts(0),
fAODvZeroCuts(0),
fvZeroCutsNames(0),
fBit(-1),
fCharge(0),
fDescription("")

{
	//
	//default constructor
	//
	
}
//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts::AliHFAssociatedTrackCuts(const AliHFAssociatedTrackCuts &source) :
AliAnalysisCuts(source),
fESDTrackCuts(source.fESDTrackCuts),
fPidObj(source.fPidObj),
fEffWeights(source.fEffWeights),
fTrigEffWeightsvspt(source.fTrigEffWeightsvspt),
fTrigEffWeightsvsptB(source.fTrigEffWeightsvsptB),
fTrigEffWeights(source.fTrigEffWeights),
fTrigEffWeightsB(source.fTrigEffWeightsB),

fPoolMaxNEvents(source.fPoolMaxNEvents), 
fPoolMinNTracks(source.fPoolMinNTracks), 
fMinEventsToMix(source.fMinEventsToMix),
fTargetFracTracks(source.fTargetFracTracks),
fNzVtxBins(source.fNzVtxBins), 
fNzVtxBinsDim(source.fNzVtxBinsDim), 
fZvtxBins(source.fZvtxBins), 
fNCentBins(source.fNCentBins), 
fNCentBinsDim(source.fNCentBinsDim), 
fCentBins(source.fCentBins),

fNofMCEventType(source.fNofMCEventType),
fMCEventType(source.fMCEventType),

fNTrackCuts(source.fNTrackCuts),
fAODTrackCuts(source.fAODTrackCuts),
fTrackCutsNames(source.fTrackCutsNames),
fNvZeroCuts(source.fNvZeroCuts),
fAODvZeroCuts(source.fAODvZeroCuts),
fvZeroCutsNames(source.fvZeroCutsNames),
fBit(source.fBit),
fCharge(source.fCharge),
fDescription(source.fDescription)
{
	//
	// copy constructor
	//
	

        AliInfo("AliHFAssociatedTrackCuts::Copy constructor ");
	if(source.fESDTrackCuts) AddTrackCuts(source.fESDTrackCuts);
	if(source.fAODTrackCuts) SetAODTrackCuts(source.fAODTrackCuts);
	if(source.fAODvZeroCuts) SetAODvZeroCuts(source.fAODvZeroCuts);
	if(source.fPidObj) SetPidHF(source.fPidObj);
	if(source.fEffWeights) SetEfficiencyWeightMap(source.fEffWeights);
	if(source.fTrigEffWeightsvspt) SetTriggerEffWeightMapvspt(source.fTrigEffWeightsvspt);
	if(source.fTrigEffWeightsvsptB) SetTriggerEffWeightMapvsptB(source.fTrigEffWeightsvsptB);
	if(source.fTrigEffWeights) SetTriggerEffWeightMap(source.fTrigEffWeights);
	if(source.fTrigEffWeightsB)SetTriggerEffWeightMapB(source.fTrigEffWeightsB);
   
    
}
//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts &AliHFAssociatedTrackCuts::operator=(const AliHFAssociatedTrackCuts &source)
{
	//
	// assignment operator
	//	
	if(&source == this) return *this;
	
	AliAnalysisCuts::operator=(source);
	fESDTrackCuts=source.fESDTrackCuts;
	fPidObj=source.fPidObj;
	fEffWeights=source.fEffWeights;
    fTrigEffWeightsvspt=source.fTrigEffWeightsvspt;
	fTrigEffWeightsvsptB=source.fTrigEffWeightsvsptB;
    fTrigEffWeights=source.fTrigEffWeights;
	fTrigEffWeightsB=source.fTrigEffWeightsB;
	fNTrackCuts=source.fNTrackCuts;
	fAODTrackCuts=source.fAODTrackCuts;
	fTrackCutsNames=source.fTrackCutsNames;
	fNvZeroCuts=source.fNvZeroCuts;
	fAODvZeroCuts=source.fAODvZeroCuts;
	fvZeroCutsNames=source.fvZeroCutsNames;
	fBit=source.fBit;
	fCharge=source.fCharge;
	
	return *this;

}


//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts::~AliHFAssociatedTrackCuts()
{
	if(fESDTrackCuts) {delete fESDTrackCuts; fESDTrackCuts = 0;}
	if(fPidObj) {delete fPidObj; fPidObj = 0;}
	if(fEffWeights){delete fEffWeights;fEffWeights=0;}
    if(fTrigEffWeightsvspt){delete fTrigEffWeightsvspt;fTrigEffWeightsvspt=0;}
	if(fTrigEffWeightsvsptB){delete fTrigEffWeightsvsptB;fTrigEffWeightsvsptB=0;}
    if(fTrigEffWeights){delete fTrigEffWeights;fTrigEffWeights=0;}
	if(fTrigEffWeightsB){delete fTrigEffWeightsB;fTrigEffWeightsB=0;}
	if(fZvtxBins) {delete[] fZvtxBins; fZvtxBins=0;} 
	if(fCentBins) {delete[] fCentBins; fCentBins=0;}
	if(fAODTrackCuts) {delete[] fAODTrackCuts; fAODTrackCuts=0;}
	if(fTrackCutsNames) {delete[] fTrackCutsNames; fTrackCutsNames=0;}
	if(fAODvZeroCuts){delete[] fAODvZeroCuts; fAODvZeroCuts=0;}
	if(fvZeroCutsNames) {delete[] fvZeroCutsNames; fvZeroCutsNames=0;}

	
}
//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::IsInAcceptance()
{
	printf("Careful: method AliHFAssociatedTrackCuts::IsInAcceptance is not implemented yet \n");
	return kFALSE;
}
//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::IsHadronSelected(AliAODTrack * track,const AliESDVertex *primary, Double_t magfield)
{
  
  AliESDtrack esdtrack(track);
  if(primary){// needed to calculate impact parameters
    // needed to calculate the impact parameters
    esdtrack.RelateToVertex(primary,magfield,3.); 
  }
  // set the TPC cluster info
  esdtrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdtrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdtrack.SetTPCPointsF(track->GetTPCNclsF());
  
  if(!fESDTrackCuts->IsSelected(&esdtrack)) return kFALSE;
  
  if(fBit>-1 && !track->TestFilterBit(fBit)) return kFALSE; // check the filter bit
  
  return kTRUE;
	
}

//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::CheckHadronKinematic(Double_t pt, Double_t d0) 
{
	
	
	
    if(pt < fAODTrackCuts[0]) {return kFALSE; cout << "reject min pt" << endl;}
	if(pt > fAODTrackCuts[1]) {return kFALSE; cout << "reject max pt" << endl;}
	if(d0 < fAODTrackCuts[2]) {return kFALSE; cout << "reject min d0" << endl;}
	if(d0 > fAODTrackCuts[3]) {return kFALSE; cout << "reject max d0" << endl;}
	
	return kTRUE;

	
}
//--------------------------------------------------------------------------

Bool_t AliHFAssociatedTrackCuts::Charge(Short_t charge, AliAODTrack* track) 
{// charge is the charge to compare to (for example, a daughter of a D meson)
	
	if(!fCharge) return kTRUE; // if fCharge is set to 0 (no selection on the charge), returns always true
	if(track->Charge()!= fCharge*charge) return kFALSE;
	return kTRUE;
}

//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::CheckKaonCompatibility(AliAODTrack * track, Bool_t useMc, TClonesArray* mcArray, Int_t method)
{
	Bool_t isKaon = kFALSE;
	
	if(useMc) { // on MC
		Int_t hadLabel = track->GetLabel();
		if(hadLabel < 0) return kFALSE;
		AliAODMCParticle* hadron = dynamic_cast<AliAODMCParticle*>(mcArray->At(hadLabel)); 
		if(hadron){
		  Int_t pdg = TMath::Abs(hadron->GetPdgCode()); 
		  if (pdg == 321) isKaon = kTRUE;
		}
	}
	
	if(!useMc) { // on DATA
          switch(method) {
	  case(1): {
		Bool_t isKTPC=kFALSE;
		Bool_t isPiTPC=kFALSE;
		Bool_t isPTPC=kFALSE;
		Bool_t isKTOF=kFALSE;
		Bool_t isPiTOF=kFALSE;
		Bool_t isPTOF=kFALSE;
		
		Bool_t KaonHyp = kFALSE;
		Bool_t PionHyp = kFALSE;
		Bool_t ProtonHyp = kFALSE;
		
		if(fPidObj->CheckStatus(track,"TOF")) {
			isKTOF=fPidObj->IsKaonRaw(track,"TOF");
			isPiTOF=fPidObj->IsPionRaw(track,"TOF");
			isPTOF=fPidObj->IsProtonRaw(track,"TOF");
		}
		if(fPidObj->CheckStatus(track,"TPC")){
			isKTPC=fPidObj->IsKaonRaw(track,"TPC");
			isPiTPC=fPidObj->IsPionRaw(track,"TPC");
			isPTPC=fPidObj->IsProtonRaw(track,"TPC");		
		}
		
		if (isKTOF && isKTPC) KaonHyp = kTRUE;
		if (isPiTOF && isPiTPC) PionHyp = kTRUE;
		if (isPTOF && isPTPC) ProtonHyp = kTRUE;
		
		if(KaonHyp && !PionHyp && !ProtonHyp) isKaon = kTRUE; 
	        break;
	      }
      case(2): {
		if(fPidObj->MakeRawPid(track,3)>=1) isKaon = kTRUE;
		break;
	      }
      }
	}
	
	return isKaon;
	
}
//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::IsKZeroSelected(AliAODv0 *vzero, AliAODVertex *vtx1)
{
	
	if(vzero->DcaV0Daughters()>fAODvZeroCuts[0]) return kFALSE;
	if(vzero->Chi2V0()>fAODvZeroCuts[1]) return kFALSE;
	if(vzero->DecayLength(vtx1) < fAODvZeroCuts[2]) return kFALSE;
	if(vzero->DecayLength(vtx1) > fAODvZeroCuts[3]) return kFALSE;
	if(vzero->OpenAngleV0() > fAODvZeroCuts[4]) return kFALSE;
	if(vzero->Pt() < fAODvZeroCuts[5]) return kFALSE;
	if(TMath::Abs(vzero->Eta()) > fAODvZeroCuts[6]) return kFALSE;

	
	return kTRUE;
}
//--------------------------------------------------------------------------
Bool_t *AliHFAssociatedTrackCuts::IsMCpartFromHF(Int_t label, TClonesArray*mcArray){
  // Check origin in MC
    
  Bool_t isCharmy = kFALSE;
  Bool_t isBeauty = kFALSE;
  Bool_t isD = kFALSE;
  Bool_t isB = kFALSE;
    
     Bool_t *originvect = new Bool_t[5];
    
    originvect[0] = kFALSE; // is from charm
	originvect[1] = kFALSE; // is from beauty
	originvect[2] = kFALSE; // is from D
	originvect[3] = kFALSE; // is from B
    originvect[4] = kFALSE; // did something go wrong? (kTRUE yes, kFALSE no))
    
    //__________________________________

    if (label<0) {originvect[4] = kTRUE; return originvect;}
    
    AliAODMCParticle* mcParticle = dynamic_cast<AliAODMCParticle*>(mcArray->At(label));
    if(!mcParticle) {originvect[4] = kTRUE; return originvect;}
    Int_t pdgCode = -1;
    Int_t  mother = mcParticle->GetMother();
    
    if(mother <0) {originvect[4] = kTRUE; return originvect;}
    
    while (mother >= 0){
        mcParticle = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother));
        if(!mcParticle) {AliError("NO MC PARTICLE"); break;}
        pdgCode =  TMath::Abs(mcParticle->GetPdgCode());
        
        mother = mcParticle->GetMother();
        
        if((pdgCode>=400 && pdgCode <500) || (pdgCode>=4000 && pdgCode<5000 )) isD = kTRUE;
        if((pdgCode>=500 && pdgCode <600) || (pdgCode>=5000 && pdgCode<6000 )) {isD = kFALSE; isB = kTRUE;}
        if(pdgCode == 4) isCharmy = kTRUE;
        if(pdgCode == 5) {isBeauty = kTRUE; isCharmy = kFALSE;}
        if(mother<0) break;
    }
    
	originvect[0] = isCharmy;
	originvect[1] = isBeauty;
	originvect[2] = isD;
	originvect[3] = isB;
 
    
  return originvect;
}

//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::InvMassDstarRejection(AliAODRecoDecayHF2Prong* d, AliAODTrack *track, Int_t hypD0) const {
	//
	// Calculates invmass of track+D0 and rejects if compatible with D*
	// (to remove pions from D*)
	// 
	Double_t nsigma = 3.;
	
	Double_t mD0, mD0bar;
	d->InvMassD0(mD0,mD0bar);
	
	Double_t invmassDstar1 = 0, invmassDstar2 = 0; 
	Double_t e1Pi = d->EProng(0,211), e2K = d->EProng(1,321); //hyp 1 (pi,K) - D0
	Double_t e1K = d->EProng(0,321), e2Pi = d->EProng(1,211); //hyp 2 (K,pi) - D0bar
	Double_t psum2 = (d->Px()+track->Px())*(d->Px()+track->Px())
	+(d->Py()+track->Py())*(d->Py()+track->Py())
	+(d->Pz()+track->Pz())*(d->Pz()+track->Pz());
	
	switch(hypD0) {
		case 1:
			invmassDstar1 = TMath::Sqrt(pow(e1Pi+e2K+track->E(0.1396),2.)-psum2);
			if ((TMath::Abs(invmassDstar1-mD0)-0.14543) < nsigma*800.*pow(10.,-6.)) return kFALSE;
			break;
		case 2:
			invmassDstar2 = TMath::Sqrt(pow(e2Pi+e1K+track->E(0.1396),2.)-psum2);
			if ((TMath::Abs(invmassDstar2-mD0bar)-0.14543) < nsigma*800.*pow(10.,-6.)) return kFALSE;
			break;
		case 3:
			invmassDstar1 = TMath::Sqrt(pow(e1Pi+e2K+track->E(0.1396),2.)-psum2);
			invmassDstar2 = TMath::Sqrt(pow(e2Pi+e1K+track->E(0.1396),2.)-psum2);
			if ((TMath::Abs(invmassDstar1-mD0)-0.14543) < nsigma*800.*pow(10.,-6.)) return kFALSE;
			if ((TMath::Abs(invmassDstar2-mD0bar)-0.14543) < nsigma*800.*pow(10.,-6.)) return kFALSE;
			break;
	}
	
	return kTRUE;
}
//________________________________________________________
void AliHFAssociatedTrackCuts::SetMCEventTypes(Int_t *MCEventTypeArray)
// set the array of event types you want to process in MonteCarlo (gluon splitting, pair production etc.)
{
	if(!fMCEventType) fMCEventType = new Int_t[fNofMCEventType];
	
	for(Int_t k=0; k<fNofMCEventType; k++){
		fMCEventType[k] = MCEventTypeArray[k];
	}
	return;	
}

//________________________________________________________
void AliHFAssociatedTrackCuts::SetAODTrackCuts(Float_t *cutsarray)
{
	if(!fAODTrackCuts) fAODTrackCuts = new Float_t[fNTrackCuts];
	for(Int_t i =0; i<fNTrackCuts; i++){
		fAODTrackCuts[i] = cutsarray[i];
	}
	SetTrackCutsNames();
	return;
}
//________________________________________________________
void AliHFAssociatedTrackCuts::SetTrackCutsNames(/*TString *namearray*/){
	
	fTrackCutsNames = new TString[4];
	fTrackCutsNames[0]= "associated track:: pt min [GeV/c]................: ";
	fTrackCutsNames[1]= "associated track:: pt max [GeV/c]................: ";
	fTrackCutsNames[2]= "associated track:: d0 min [cm]...................: ";
	fTrackCutsNames[3]= "associated track:: d0 max [cm]...................: ";
	

	
	return;
}
//--------------------------------------------------------------------------
void AliHFAssociatedTrackCuts::SetAODvZeroCuts(Float_t *cutsarray)
{
	

	if(!fAODvZeroCuts) fAODvZeroCuts = new Float_t[fNvZeroCuts];
	for(Int_t i =0; i<fNvZeroCuts; i++){
		fAODvZeroCuts[i] = cutsarray[i];
	}
	SetvZeroCutsNames();
	return;
}
//--------------------------------------------------------------------------
void AliHFAssociatedTrackCuts::SetvZeroCutsNames(/*TString *namearray*/){
	
	fvZeroCutsNames = new TString[7];
	fvZeroCutsNames[0] = "vZero:: max DCA between two daughters [cm].......: ";
	fvZeroCutsNames[1] = "vZero:: max fit Chi Square between two daughters.: ";
	fvZeroCutsNames[2] = "vZero:: min decay length [cm]....................: ";
	fvZeroCutsNames[3] = "vZero:: max decay length [cm]....................: ";
	fvZeroCutsNames[4] = "vZero:: max opening angle between daughters [rad]: ";
	fvZeroCutsNames[5] = "vZero:: pt min [Gev/c]...........................: ";
	fvZeroCutsNames[6] = "vZero:: |Eta| range <............................: ";
	
	
	return;
}

//--------------------------------------------------------------------------
void AliHFAssociatedTrackCuts::SetPidAssociated()
{
  //setting PidResponse
  if(fPidObj->GetOldPid()==kFALSE && fPidObj->GetPidResponse()==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
    fPidObj->SetPidResponse(pidResp);
  }
}
//--------------------------------------------------------------------------

void AliHFAssociatedTrackCuts::Print(Option_t *option) const
{
  /// overloaded from TObject: print info
  if (strcmp(option, "parameters")==0) {
    PrintPoolParameters();
    return;
  } else if (strcmp(option, "selectedMC")==0) {
    PrintSelectedMCevents();
    return;
  }
  PrintAll();
}

//--------------------------------------------------------------------------
Int_t AliHFAssociatedTrackCuts::GetPoolBin(Double_t multorcent, Double_t zVtx) const
{
 
    Int_t poolbin = -1;
    Int_t centbin = -1;
    Int_t zvtxbin = -1;
    
    
    if(multorcent <fCentBins[0]) return poolbin;
    if(zVtx <fZvtxBins[0]) return poolbin;
    
    
    for (Int_t i=0;i<fNCentBins;i++){
        if(multorcent<fCentBins[i+1]) {
            centbin=i;
            break;
        }
    }
    
    for (Int_t i=0;i<fNzVtxBins;i++){
        if(zVtx<fZvtxBins[i+1]) {
            zvtxbin=i;
            break;
        }
    }

    if(centbin!=-1 && zvtxbin!=-1) poolbin = centbin + zvtxbin*fNCentBins;
    
    return poolbin;
}
//--------------------------------------------------------------------------
void AliHFAssociatedTrackCuts::PrintAll() const
{
	
	if(fDescription){
	  printf("=================================================");
	  printf("\nAdditional description\n");
	  std::cout << fDescription << std::endl;
	  printf("\n");
	}
	printf("\n=================================================");
	if(fESDTrackCuts){
	  printf("\nCuts for the associated track: \n \n");
	  
	  printf("ITS Refit........................................: %s\n",fESDTrackCuts->GetRequireITSRefit() ? "Yes" : "No");
	  printf("TPC Refit........................................: %s\n",fESDTrackCuts->GetRequireTPCRefit() ? "Yes" : "No");
	  printf("ITS SA...........................................: %s\n",fESDTrackCuts->GetRequireITSStandAlone() ? "Yes" : "No");
	  printf("TPC SA...........................................: %s\n",fESDTrackCuts->GetRequireTPCStandAlone() ? "Yes" : "No");
	  printf("Min number of ITS clusters.......................: %d\n",fESDTrackCuts->GetMinNClustersITS());
	  printf("Min number of TPC clusters.......................: %d\n",fESDTrackCuts->GetMinNClusterTPC());
	  Int_t spd = fESDTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD);
	  if(spd==0) std::cout <<  "SPD..............................................: kOff"  << std::endl;
	  if(spd==1) std::cout <<  "SPD..............................................: kNone"  << std::endl;
	  if(spd==2) std::cout <<  "SPD..............................................: kAny"  << std::endl;
	  if(spd==3) std::cout <<  "SPD..............................................: kFirst"  << std::endl;
	  if(spd==4) std::cout <<  "SPD..............................................: kOnlyFirst"  << std::endl;
	  if(spd==5) std::cout <<  "SPD..............................................: kSecond"  << std::endl;
	  if(spd==6) std::cout <<  "SPD..............................................: kOnlySecond"  << std::endl;
	  if(spd==7) std::cout <<  "SPD..............................................: kBoth"  << std::endl;
	}
	else printf("\nNo Cuts for Associated Tracks\n");
	std::cout <<  "Filter Bit.......................................: " << fBit  << std::endl;
	std::cout <<  "Charge...........................................: " << fCharge  << std::endl;
	
	if(fAODTrackCuts){
	  for(Int_t j=0;j<fNTrackCuts;j++){
	    std::cout << fTrackCutsNames[j] << fAODTrackCuts[j] << std::endl;
	  }
	}

	if(fAODvZeroCuts){
	  printf("\n");
	  printf("=================================================");
	  printf("\nCuts for the K0 candidates: \n \n");
	  for(Int_t k=0;k<fNvZeroCuts;k++){
	    std::cout << fvZeroCutsNames[k] <<  fAODvZeroCuts[k] << std::endl;
	  }
	}
	else printf("\nNo Cuts for the K0 candidates\n");
	std::cout << " " << std::endl;
	PrintPoolParameters();
	PrintSelectedMCevents();

}

//--------------------------------------------------------------------------
void AliHFAssociatedTrackCuts::PrintPoolParameters() const
{   
	printf("=================================================");
	printf("\nEvent Pool settings: \n \n");
	
	printf("Number of zVtx Bins: %d\n", fNzVtxBins);
	printf("\nzVtx Bins:\n");
	//Double_t zVtxbinLims[fNzVtxBins+1] = fNzVtxBins;
	for(Int_t k=0; k<fNzVtxBins; k++){
		printf("Bin %d..............................................: %.1f - %.1f cm\n", k, fZvtxBins[k], fZvtxBins[k+1]);	
	}
	printf("\n");
	printf("\nNumber of Centrality(multiplicity) Bins: %d\n", fNCentBins);
	printf("\nCentrality(multiplicity) Bins:\n");
	for(Int_t k=0; k<fNCentBins; k++){
		printf("Bin %d..............................................: %.1f - %.1f\n", k, fCentBins[k], fCentBins[k+1]);
	}

	
	
}
//--------------------------------------------------------------------------

Double_t AliHFAssociatedTrackCuts::GetTrackWeight(Double_t pt, Double_t eta,Double_t zvtx){
  if(!fEffWeights)return 1.;
  
  Int_t bin=fEffWeights->FindBin(pt,eta,zvtx);
  if(fEffWeights->IsBinUnderflow(bin)||fEffWeights->IsBinOverflow(bin))return 1.;
  return fEffWeights->GetBinContent(bin);

}


//--------------------------------------------------------------------------
Double_t AliHFAssociatedTrackCuts::GetTrigWeight(Double_t pt, Double_t mult){
    
    
    
    if(fTrigEffWeightsvspt){
       Int_t bin=fTrigEffWeightsvspt->FindBin(pt);
        if(fTrigEffWeightsvspt->IsBinUnderflow(bin)||fTrigEffWeightsvspt->IsBinOverflow(bin))return 1.;
        return fTrigEffWeightsvspt->GetBinContent(bin);
        
    }
    
    if(fTrigEffWeights){
        Int_t bin=fTrigEffWeights->FindBin(pt,mult);
        if(fTrigEffWeights->IsBinUnderflow(bin)||fTrigEffWeights->IsBinOverflow(bin))return 1.;
        return fTrigEffWeights->GetBinContent(bin);
        
    }
    
    //if(!fTrigEffWeights && !fTrigEffWeightsvspt)return 1.;
    
    return 1.;
    
}

//--------------------------------------------------------------------------
Double_t AliHFAssociatedTrackCuts::GetTrigWeightB(Double_t pt, Double_t mult){
    
    if(fTrigEffWeightsvsptB){
        Int_t bin=fTrigEffWeightsvsptB->FindBin(pt);
        if(fTrigEffWeightsvsptB->IsBinUnderflow(bin)||fTrigEffWeightsvsptB->IsBinOverflow(bin))return 1.;
        return fTrigEffWeightsvsptB->GetBinContent(bin);
        
    }
    
    if(fTrigEffWeightsB){
        Int_t bin=fTrigEffWeightsB->FindBin(pt,mult);
        if(fTrigEffWeightsB->IsBinUnderflow(bin)||fTrigEffWeightsB->IsBinOverflow(bin))return 1.;
        return fTrigEffWeightsB->GetBinContent(bin);
        
    }
    
 //   if(!fTrigEffWeightsB && !fTrigEffWeightsvsptB)return 1.;
    return 1;
}
//--------------------------------------------------------------------------
void AliHFAssociatedTrackCuts::PrintSelectedMCevents() const
{
	printf("\n=================================================");
	
	printf("\nSelected MC events: \n \n");
	printf("Number of selected events: %d\n",fNofMCEventType);
	
	for(Int_t k=0; k<fNofMCEventType; k++){
	if(fMCEventType[k]==28)	printf("=> Flavour excitation \n");	
	if(fMCEventType[k]==53)	printf("=> Pair creation \n");	
	if(fMCEventType[k]==68)	printf("=> Gluon splitting \n");	
	}
	
	printf("\n");
	for(Int_t k=0; k<fNofMCEventType; k++){
		printf("MC process code %d \n",fMCEventType[k]);		
	}
	
	printf("\n");
	
	
		
	
}


