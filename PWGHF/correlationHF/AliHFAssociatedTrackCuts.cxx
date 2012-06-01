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
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "TString.h"




ClassImp(AliHFAssociatedTrackCuts)

//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts::AliHFAssociatedTrackCuts():
AliAnalysisCuts(),
fESDTrackCuts(0),
fPidObj(0),
fNTrackCuts(0),
fAODTrackCuts(0),
fTrackCutsNames(0),
fNvZeroCuts(0),
fAODvZeroCuts(0),
fvZeroCutsNames(0)

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
fNTrackCuts(0),
fAODTrackCuts(0),
fTrackCutsNames(0),
fNvZeroCuts(0),
fAODvZeroCuts(0),
fvZeroCutsNames(0)

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
fNTrackCuts(source.fNTrackCuts),
fAODTrackCuts(source.fAODTrackCuts),
fTrackCutsNames(source.fTrackCutsNames),
fNvZeroCuts(source.fNvZeroCuts),
fAODvZeroCuts(source.fAODvZeroCuts),
fvZeroCutsNames(source.fvZeroCutsNames)
{
	//
	// copy constructor
	//
	

	cout << "AliHFAssociatedTrackCuts::Copy constructor " << endl;
	if(source.fESDTrackCuts) AddTrackCuts(source.fESDTrackCuts);
	if(source.fAODTrackCuts) SetAODTrackCuts(source.fAODTrackCuts);
	if(source.fAODvZeroCuts) SetAODvZeroCuts(source.fAODvZeroCuts);
	if(source.fPidObj) SetPidHF(source.fPidObj);
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
	fNTrackCuts=source.fNTrackCuts;
	fAODTrackCuts=source.fAODTrackCuts;
	fTrackCutsNames=source.fTrackCutsNames;
	fNvZeroCuts=source.fNvZeroCuts;
	fAODvZeroCuts=source.fAODvZeroCuts;
	fvZeroCutsNames=source.fvZeroCutsNames;
	
	return *this;

}


//--------------------------------------------------------------------------
AliHFAssociatedTrackCuts::~AliHFAssociatedTrackCuts()
{
	if(fESDTrackCuts) {delete fESDTrackCuts; fESDTrackCuts = 0;}
	if(fPidObj) {delete fPidObj; fPidObj = 0;}
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
Bool_t AliHFAssociatedTrackCuts::IsHadronSelected(AliAODTrack * track, AliAODVertex *vtx1, Double_t bz) 
{
	
	
	
	Double_t d0 = -999.;
	Double_t d0err = -998.;
	Double_t VeryBig = 100;
    
	AliESDtrack esdtrack(track);
	Double_t d0z0[2],covd0z0[3];
	
	esdtrack.PropagateToDCA(vtx1,bz,VeryBig,d0z0,covd0z0);
	
	d0 = TMath::Abs(d0z0[0]);
	d0err = TMath::Sqrt(covd0z0[0]);
	
	if(!fESDTrackCuts->IsSelected(&esdtrack)) return kFALSE;
	
	if(track->Pt() < fAODTrackCuts[0]) return kFALSE;
	if(track->Pt() > fAODTrackCuts[1]) return kFALSE;
	if(d0 < fAODTrackCuts[2])  return kFALSE;
	if(d0 > fAODTrackCuts[3])  return kFALSE;
	
	return kTRUE;

	
}
//--------------------------------------------------------------------------
Bool_t AliHFAssociatedTrackCuts::CheckKaonCompatibility(AliAODTrack * track, Bool_t useMc, TClonesArray* mcArray)
{
	Bool_t isKaon = kFALSE;
	
	if(useMc) { // on MC
		Int_t hadLabel = track->GetLabel();
		if(hadLabel < 0) return kFALSE;
		AliAODMCParticle* hadron = dynamic_cast<AliAODMCParticle*>(mcArray->At(hadLabel));
		Int_t pdg = TMath::Abs(hadron->GetPdgCode()); 
		if (pdg == 321) isKaon = kTRUE;
	}
	
	if(!useMc) { // on DATA
		
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
Int_t AliHFAssociatedTrackCuts::IsMCpartFromHF(Int_t label, TClonesArray*mcArray)

{
AliAODMCParticle* MCParticle;
Int_t generationpdgcode = -1;
if (label<0) return 0;
Bool_t isCharmy = kFALSE;
Bool_t isBeauty = kFALSE;
Int_t generation = 0;
while(generationpdgcode!=2212){	// loops back to the collision to check the particle source	
	MCParticle = dynamic_cast<AliAODMCParticle*>(mcArray->At(label));
	generationpdgcode =  TMath::Abs(MCParticle->GetPdgCode());
	label = MCParticle->GetMother();
	generation++;
	if(generationpdgcode == 4) isCharmy = kTRUE;
		if(generationpdgcode == 5) {isBeauty = kTRUE; isCharmy = kFALSE;}
	if(label<0) break;
	
}
if (isCharmy) return 4;
else if (isBeauty) return 5;
else return 0;

	
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
void AliHFAssociatedTrackCuts::PrintAll()
{
	printf("\nCuts for the associated track: \n \n");
	       
	printf("ITS Refit........................................: %s\n",fESDTrackCuts->GetRequireITSRefit() ? "Yes" : "No");
	printf("TPC Refit........................................: %s\n",fESDTrackCuts->GetRequireTPCRefit() ? "Yes" : "No");
	printf("ITS SA...........................................: %s\n",fESDTrackCuts->GetRequireITSStandAlone() ? "Yes" : "No");
	printf("TPC SA...........................................: %s\n",fESDTrackCuts->GetRequireTPCStandAlone() ? "Yes" : "No");
	printf("Min number of ITS clusters.......................: %d\n",fESDTrackCuts->GetMinNClustersITS());
	printf("Min number of TPC clusters.......................: %d\n",fESDTrackCuts->GetMinNClusterTPC());
	Int_t spd = fESDTrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD);
	if(spd==0) cout <<  "SPD..............................................: kOff"  << endl;
	if(spd==1) cout <<  "SPD..............................................: kNone"  << endl;
	if(spd==2) cout <<  "SPD..............................................: kAny"  << endl;
	if(spd==3) cout <<  "SPD..............................................: kFirst"  << endl;
	if(spd==4) cout <<  "SPD..............................................: kOnlyFirst"  << endl;
	if(spd==5) cout <<  "SPD..............................................: kSecond"  << endl;
	if(spd==6) cout <<  "SPD..............................................: kOnlySecond"  << endl;
	if(spd==7) cout <<  "SPD..............................................: kBoth"  << endl;
		
	for(Int_t j=0;j<fNTrackCuts;j++){
		cout << fTrackCutsNames[j] << fAODTrackCuts[j] << endl;
	}
	
	printf("\nCuts for the K0 candidates: \n \n");
	for(Int_t k=0;k<fNvZeroCuts;k++){
		cout << fvZeroCutsNames[k] <<  fAODvZeroCuts[k] << endl;
	}
	cout << " " << endl;
}
