/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
// AliCDMesonUtilsStripped
//
//  Author:
//  Felix Reidt <Felix.Reidt@cern.ch>

#include <TH1.h>
#include <TH2.h>
#include <TGeoMatrix.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TParticle.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDFMD.h"
#include "AliGeomManager.h"
#include "AliITSAlignMille2Module.h"
#include "AliITSsegmentationSPD.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliSPDUtils.h"
#include "AliTriggerAnalysis.h"

#include "AliAODTracklets.h"
#include "AliAODEvent.h"

#include "AliCDMesonBaseStripped.h"
#include "AliCDMesonTracks.h"

#include "AliCDMesonUtilsStripped.h"


//==============================================================================
//------------------------------------------------------------------------------
Bool_t AliCDMesonUtilsStripped::CutEvent(const AliESDEvent *ESDEvent)
{
	//
	// CutEvent
	//

	AliTriggerAnalysis triggerAnalysis;

	// collision vertex cut
	// A cut in XY is implicitly done during the reconstruction by constraining
	// the vertex to the beam diamond.

	// Primary vertex
	Bool_t kpr0 = kTRUE;
	const AliESDVertex *vertex = ESDEvent->GetPrimaryVertexTracks();
	if(vertex->GetNContributors()<1) {
		// SPD vertex
		vertex = ESDEvent->GetPrimaryVertexSPD();
		if(vertex->GetNContributors()<1) {
			// NO GOOD VERTEX, SKIP EVENT
			kpr0 = kFALSE;
		}
	}
	const Bool_t kpriv = kpr0 && (fabs(ESDEvent->GetPrimaryVertex()->GetZ())<10.);
	// 10 is the common value, unit: cm
	if(!kpriv)
		return kFALSE;

	return kTRUE;
}

//------------------------------------------------------------------------------
Bool_t AliCDMesonUtilsStripped::CutEvent(const AliAODEvent *AODEvent)
{
	//
	// Cut Event for AOD Events, to be combined with the ESD Track Cut
	//

	// TODO: no idea about fast or yet, to be thought of

	// Primary vertex
	Bool_t kpr0 = kTRUE;
	const AliAODVertex *vertex = AODEvent->GetPrimaryVertex();
	if(vertex->GetNContributors()<1) {
		// SPD vertex
		vertex = AODEvent->GetPrimaryVertexSPD();
		if(vertex->GetNContributors()<1) {
			// NO GOOD VERTEX, SKIP EVENT
			kpr0 = kFALSE;
		}
	}
	const Bool_t kpriv = kpr0 && (fabs(AODEvent->GetPrimaryVertex()->GetZ())<10.);
	// 10 is the common value, unit: cm

	if(!kpriv)
		return kFALSE;

	return kTRUE;
}


//------------------------------------------------------------------------------
Int_t AliCDMesonUtilsStripped::GetGapConfig(const AliESDEvent *ESDEvent)
{
	//
	// GetGapConfigAndTracks
	//
	// retrieves the gap configuration of a track and returns it as
	// an bit vector
	// kBaseLine ensures, that this event is valid
	// + is equivalent to | in this case
	return AliCDMesonBaseStripped::kBitBaseLine
		+ GetV0(ESDEvent) + GetFMD(ESDEvent) + GetSPD(ESDEvent) + GetTPC(ESDEvent);
}


//==============================================================================
//------------------------------------------------------------------------------
Int_t AliCDMesonUtilsStripped::GetV0(const AliESDEvent * ESDEvent)
{
	//
	//GetV0
	//

	AliTriggerAnalysis triggerAnalysis;
	const Bool_t khw = kFALSE;
	const Bool_t v0A =
		(triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kASide, khw) ==
		 AliTriggerAnalysis::kV0BB);
	const Bool_t v0C =
		(triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kCSide, khw) ==
		 AliTriggerAnalysis::kV0BB);

	return v0A * AliCDMesonBaseStripped::kBitV0A
		+ v0C * AliCDMesonBaseStripped::kBitV0C;
}


//------------------------------------------------------------------------------
Int_t AliCDMesonUtilsStripped::GetFMD(const AliESDEvent *ESDEvent)
{
	//
	// GetFMD
	//

	AliTriggerAnalysis triggerAnalysis;
	triggerAnalysis.SetFMDThreshold(0.3, 0.5); // parameters got from FMD
	const Bool_t fmdA =
		triggerAnalysis.FMDTrigger(ESDEvent, AliTriggerAnalysis::kASide);
	const Bool_t fmdC =
		triggerAnalysis.FMDTrigger(ESDEvent, AliTriggerAnalysis::kCSide);

	return fmdA * AliCDMesonBaseStripped::kBitFMDA
		+ fmdC * AliCDMesonBaseStripped::kBitFMDC;
}


//------------------------------------------------------------------------------
Int_t AliCDMesonUtilsStripped::GetSPD(const AliESDEvent *ESDEvent)
{
	//
	// GetSPD
	//

	Int_t nfoctr[10];
	GetNFO(ESDEvent, "]0.9[", nfoctr);
	// get multiplicity from fastOR and fill corresponding hit maps


	const Int_t ipA = nfoctr[kIPA]; // inner layer A side
	const Int_t ipC = nfoctr[kIPC]; // inner layer C side
	const Int_t opA = nfoctr[kOPA]; // outer layer A side
	const Int_t opC = nfoctr[kOPC]; // outer layer C side

	const Bool_t spdA = ipA + opA; // A side hit?
	const Bool_t spdC = ipC + opC; // C side hit?

	return spdA * AliCDMesonBaseStripped::kBitSPDA
		+ spdC * AliCDMesonBaseStripped::kBitSPDC;
}


//------------------------------------------------------------------------------
Int_t AliCDMesonUtilsStripped::GetTPC(const AliESDEvent * ESDEvent)
{
	//
	//GetTPC
	//

	const Double_t etacut = 0.9;
	Int_t nA = 0;
	Int_t nC = 0;
	for(Int_t itrack = 0; itrack < ESDEvent->GetNumberOfTracks(); itrack++){
		const AliESDtrack* esdtrack = ESDEvent->GetTrack(itrack);
		if( esdtrack->Eta() > etacut ){
			nA ++;
		}
		else if( esdtrack->Eta() < -etacut ){
			nC ++;
		}
	}

	const Bool_t tpcA = nA;
	const Bool_t tpcC = nC;

	return tpcA * AliCDMesonBaseStripped::kBitTPCA
		+ tpcC * AliCDMesonBaseStripped::kBitTPCC;
}

//==============================================================================
//------------------------------------------------------------------------------
void AliCDMesonUtilsStripped::SPDLoadGeom(Int_t run)
{
	// method to get the gGeomanager
	// it is called at the CreatedOutputObject stage
	// to comply with the CAF environment

	AliCDBManager *man = AliCDBManager::Instance();
	// WARNING THE OCDB PATH SHOULD BE ADJUSTED TO THE RUNNING CONDITIONS


	TString cdbpath;
	if (man->IsDefaultStorageSet()) {
		const AliCDBStorage *dsto = man->GetDefaultStorage();
		cdbpath = TString(dsto->GetBaseFolder());
	}
	else {
		man->SetDefaultStorage("raw://");
		cdbpath = "raw://";
	}

	man->SetSpecificStorage("ITS/Align/Data",cdbpath);
	man->SetSpecificStorage("GRP/Geometry/Data",cdbpath);
	man->SetRun(run);

	AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
	if (!obj) {
		printf("AliCDMesonUtilsStripped failed loading geometry object\n");
		return;
	}
	AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
	AliGeomManager::ApplyAlignObjsFromCDB("ITS");
}

//------------------------------------------------------------------------------
Bool_t AliCDMesonUtilsStripped::SPDLoc2Glo(Int_t id, const Double_t *loc,
                                           Double_t *glo)
{
	//
	//SPDLoc2Glo, do not touch
	//

	static TGeoHMatrix mat;
	Int_t vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
	if (vid<0) {
		printf("AliCDMesonUtilsStripped Did not find module with such ID %d\n",id);
		return kFALSE;
	}
	AliITSAlignMille2Module::SensVolMatrix(vid,&mat);
	mat.LocalToMaster(loc,glo);
	return kTRUE;
}


//------------------------------------------------------------------------------
Int_t AliCDMesonUtilsStripped::CheckChipEta(Int_t chipKey,
                                            TString scut,
                                            const Double_t vtxPos[])
{
	//
	//CheckChipEta
	//

	// retrieves the position in eta for a given chip and applies the cut
	// results:
	// 0 <= out of range
	// -1 <= negative pseudo-rapidity position, in range (C-Side)
	// 1 <= positive pseudo-rapidity position, in range (A-Side)
	//
	// scut: "[0.9" or "]0.9", only 3 digits for the value!!


	const Bool_t kincl = (scut[0] == '[');
	const TString cutval = scut(1,3);
	const Double_t etacut = fabs(cutval.Atof());

	//no eta cut, save time
	if(kincl && etacut>=2)
		return kTRUE;

	Int_t etaside = 1;
	//------------------------------- NOT TO TOUCH ------------------------>>
	UInt_t module=999, offchip=999;
	AliSPDUtils::GetOfflineFromOfflineChipKey(chipKey,module,offchip);
	UInt_t hs = AliSPDUtils::GetOnlineHSFromOffline(module);
	if(hs<2) offchip = 4 - offchip; // inversion  in the inner layer...

	const Int_t col[]={
		hs<2? 0 : 31,
		hs<2? 31 : 0,
		hs<2? 31 : 0,
		hs<2? 0 : 31};
	const Int_t aa[]={0, 0, 255, 255};
	const AliITSsegmentationSPD seg;

	for(Int_t ic=0; ic<4; ic++){
		Float_t localchip[3]={0.,0.,0.};
		seg.DetToLocal(aa[ic],col[ic]+32*offchip,localchip[0],localchip[2]);
		// local coordinate of the chip center
		//printf("local coordinates %d %d: %f %f \n",chipKey, ic, localchip[0],localchip[2]);
		const Double_t local[3] = {localchip[0],localchip[1],localchip[2]};
		Double_t glochip[3]={0.,0.,0.};
		if(!SPDLoc2Glo(module,local,glochip)){
			return kFALSE;
		}

		//-------------------------------------------------------------------<<

		const TVector3 pos(glochip[0]-vtxPos[0], glochip[1]-vtxPos[1],
		                   glochip[2]-vtxPos[2]);
		//pos.Print();

		if( kincl && fabs(pos.Eta()) > etacut)
			return kFALSE;

		if(!kincl){
			if(fabs(pos.Eta()) < etacut)
				return kFALSE;
			else if(pos.Eta()<0)
				etaside = -1;
			else
				etaside = 1;
		}
	}

	return etaside;
}


//------------------------------------------------------------------------------
void AliCDMesonUtilsStripped::GetNFO(const AliESDEvent *ESDEvent,
                                     TString etacut, Int_t ctr[])
{
	//
	// GetNFO
	//
	// analyzes the SPD fastOR for a given eta range and returns
	// an array with the number of hits in:

	Int_t ninner=0; // inner layer
	Int_t nouter=0; // outer layer
	Int_t ipA = 0; // inner layer A side
	Int_t ipC = 0; // inner layer C side
	Int_t opA = 0; // outer layer A side
	Int_t opC = 0; // outer layer C side

	const AliMultiplicity *mult = ESDEvent->GetMultiplicity();

	// position of the primary vertex
	Double_t tmp[3] = { 0., 0., 0. };
	ESDEvent->GetPrimaryVertex()->GetXYZ(tmp);
	Double_t vtxPos[3] = { tmp[0], tmp[1], tmp[2] };


	for(Int_t iChipKey=0; iChipKey < 1200; iChipKey++){
		if(mult->TestFastOrFiredChips(iChipKey)){
			// here you check if the FastOr bit is 1 or 0
			const Int_t iseta = CheckChipEta(iChipKey, etacut, vtxPos);
			if(iseta==0)
				continue;

			if(iChipKey<400) {
				ninner++;  // here you count the FastOr bits in the inner layer
				if(iseta>0)
					ipA ++;
				else
					ipC ++;
			}
			else {
				nouter++;  // here you count the FastOr bits in the outer layer
				if(iseta>0)
					opA ++;
				else
					opC ++;
			}
		}
	}

	ctr[kInnerPixel]= ninner;
	ctr[kOuterPixel]= nouter;
	ctr[kIPA]= ipA;
	ctr[kIPC]= ipC;
	ctr[kOPA]= opA;
	ctr[kOPC]= opC;

	return;
}
