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
// AliPBUtils
// for
// AliAnalysisTaskPB
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>

#include <TH1.h>
#include <TH2.h>
#include <TGeoMatrix.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TParticle.h>
#include <TObjArray.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliGeomManager.h"
#include "AliITSAlignMille2Module.h"
#include "AliITSsegmentationSPD.h"
#include "AliMultiplicity.h"
#include "AliSPDUtils.h"
#include "AliTriggerAnalysis.h"
#include "AliVZEROTriggerData.h"
#include "AliVZEROCalibData.h"

#include "AliAODTracklets.h"
#include "AliAODEvent.h"

#include "AliPBBase.h"
#include "AliPBTracks.h"

#include "AliPBUtils.h"


//==============================================================================
//------------------------------------------------------------------------------
void AliPBUtils::GetMassPtCtsOA(Int_t pid,
                                     const TVector3* momenta[], Double_t & mass,
                                     Double_t &pt, Double_t &cts, Double_t &oa)
{
	//
	// Get Mass Pt Cts OA
	//

	Double_t tmpp1[3], tmpp2[3];
	momenta[0]->GetXYZ(tmpp1);
	momenta[1]->GetXYZ(tmpp2);

	Double_t masstrk = -999;
	if(pid==AliPBBase::kBinPionE || pid==AliPBBase::kBinPion ||
	   pid==AliPBBase::kBinSinglePion)
		masstrk = AliPID::ParticleMass(AliPID::kPion);
	else if(pid==AliPBBase::kBinKaonE || pid==AliPBBase::kBinKaon ||
	        pid==AliPBBase::kBinSingleKaon)
		masstrk = AliPID::ParticleMass(AliPID::kKaon);
	else if(pid==AliPBBase::kBinProtonE || pid==AliPBBase::kBinProton ||
	        pid==AliPBBase::kBinSingleProton)
		masstrk = AliPID::ParticleMass(AliPID::kProton);
	else if(pid==AliPBBase::kBinElectronE ||
	        pid==AliPBBase::kBinElectron ||
	        pid==AliPBBase::kBinSingleElectron)
		masstrk = AliPID::ParticleMass(AliPID::kElectron);
	else
		masstrk = AliPID::ParticleMass(AliPID::kPion);
    
  //printf("Particle mass: %f\n",masstrk);

	const TLorentzVector sumv = GetKinematics(tmpp1, tmpp2, masstrk, masstrk,
	                                          cts);
	mass = sumv.M();
	pt = sumv.Pt();

	oa = GetOA(tmpp1, tmpp2);
}


//------------------------------------------------------------------------------
void AliPBUtils::GetPWAinfo(Int_t pid, const AliVTrack *trks[],
    Float_t& theta, Float_t& phi, Float_t& mass,Float_t momentum[])
{
	// transforms the coordiante system to the helicity frame and determines
	// invariant mass and 3-vector of the two track system
	//

	Double_t tmpp1[3], tmpp2[3]; // 3-vectors of the daughter tracks
	trks[0]->GetPxPyPz(tmpp1);
	trks[1]->GetPxPyPz(tmpp2);

	// determine mass of the daugther tracks
	Double_t masstrk = -999;
	if (pid==AliPBBase::kBinPion || pid==AliPBBase::kBinSinglePion)
		masstrk = AliPID::ParticleMass(AliPID::kPion);
	else if(pid==AliPBBase::kBinKaon || pid==AliPBBase::kBinSingleKaon)
		masstrk = AliPID::ParticleMass(AliPID::kKaon);
	else if(pid==AliPBBase::kBinProton ||
	        pid==AliPBBase::kBinSingleProton)
		masstrk = AliPID::ParticleMass(AliPID::kProton);
	else if(pid==AliPBBase::kBinElectron ||
	        pid==AliPBBase::kBinSingleElectron)
		masstrk = AliPID::ParticleMass(AliPID::kElectron);
	else
		masstrk = AliPID::ParticleMass(AliPID::kPion);

	TLorentzVector va, vb; // 4-vectors of the two tracks
	va.SetXYZM(tmpp1[0], tmpp1[1], tmpp1[2], masstrk);
	vb.SetXYZM(tmpp2[0], tmpp2[1], tmpp2[2], masstrk);

	TLorentzVector sumv = va+vb; // 4-vector of the pair
	TVector3 sumXZ(sumv.X(), 0, sumv.Z()); // its projection to the xz-plane

	// mass and momentum
	mass = sumv.M();
	momentum[0] = sumv.X();
	momentum[1] = sumv.Y();
	momentum[2] = sumv.Z();

	// coordinate axis vectors
	const TVector3 x(1.,0,0);
	const TVector3 y(0,1.,0);
	const TVector3 z(0,0,1.);

	const Double_t alpha = -z.Angle(sumXZ);

	sumXZ.Rotate(alpha, y);
	sumv.Rotate(alpha, y);
	va.Rotate(alpha, y);
	vb.Rotate(alpha, y);

	const Double_t beta = z.Angle(sumv.Vect());
	sumv.Rotate(beta, x);
	va.Rotate(beta, x);
	vb.Rotate(beta, x);

	const TVector3 bv = -sumv.BoostVector();
	va.Boost(bv);
	vb.Boost(bv);

	// angles
	theta = va.Theta();
	phi = va.Phi();
}

//------------------------------------------------------------------------------
Int_t AliPBUtils::GetEventType(const AliESDEvent *ESDEvent)
{
	// checks of which type a event is:
	// beam-beam interaction (I), beam-gas (A/C), empty (E)
	//

	TString firedTriggerClasses = ESDEvent->GetFiredTriggerClasses();
  
  printf("<AliPBUtils::GetEventType> firedTriggerClasses: %s\n",
    firedTriggerClasses.Data());

	if (firedTriggerClasses.Contains("CINT1A-ABCE-NOPF-ALL")) { // A
		return AliPBBase::kBinEventA;
	}
	if (firedTriggerClasses.Contains("CINT1C-ABCE-NOPF-ALL")) { // C
		return AliPBBase::kBinEventC;
	}
	if (firedTriggerClasses.Contains("CINT1B-ABCE-NOPF-ALL")) { // I
		return AliPBBase::kBinEventI;
	}
	if (firedTriggerClasses.Contains("CINT1-E-NOPF-ALL")) { // E
		return AliPBBase::kBinEventE;
	}
	if (firedTriggerClasses.Contains("CDG5-E")) { // E
		return AliPBBase::kBinEventE;
	}
	if (firedTriggerClasses.Contains("CDG5-I")) { // I
		return AliPBBase::kBinEventI;
	}
	if (firedTriggerClasses.Contains("CDG5-AC")) { // AC
		return AliPBBase::kBinEventAC;
	}
	return AliPBBase::kBinEventUnknown;
}


//------------------------------------------------------------------------------
void AliPBUtils::SwapTrack(const AliVTrack* trks[])
{
	//
  // swap two esdtracks, needed since only the properties of one a stored and
	// AliRoot sorts them
	//

	TRandom3 tmprd(0);
	if(tmprd.Rndm()>0.5)
		return;

	const AliVTrack* tmpt = trks[0];
	trks[0]=trks[1];
	trks[1]=tmpt;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetCombCh(const AliVTrack * trks[])
{
	//
	// get combination of charges
	//

	const Double_t ch1 = trks[0]->Charge();
	const Double_t ch2 = trks[1]->Charge();

	if(ch1*ch2<0){
		return AliPBBase::kBinPM;
	}
	else 
		return AliPBBase::kBinPPMM;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetCombPID(AliPIDResponse* pid, const AliVTrack* trks[],
                                  Int_t mode, TH2* comb2trkPID /*= 0x0 */)
{
	//
	// Get PID for a two track system (data)
	//

	// printf("<I - AliPBUtils::GetCombPID> mode: %i\n",mode);
  
  if (!pid){
		return AliPBBase::kBinPIDUnknown;
	}

	// get PID for the single tracks
	Int_t kpid[2];
	for(Int_t ii=0; ii<2; ii++){
		kpid[ii] = GetPID(pid, trks[ii], mode);
	}
  //printf("<AliPBUtils::GetCombPID> pids: %i %i\n",kpid[0],kpid[1]);
  
	// PID QA histogram
	if (comb2trkPID) {
		comb2trkPID->Fill(kpid[0], kpid[1]);
	}

	return CombinePID(kpid);
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetCombPID(const TParticle* particles[],
                                  TH2 *comb2trkPID /* = 0x0 */)
{
	//
	// Get PID for a two track system (MC)
	//

	// get PID for the single tracks
	Int_t kpid[2];
	for(Int_t ii=0; ii<2; ii++){
		kpid[ii] = GetPID(particles[ii]->GetPdgCode());
	}

	// PID QA histogram
	if (comb2trkPID) {
		comb2trkPID->Fill(kpid[0], kpid[1]);
	}

	return CombinePID(kpid);
}


//------------------------------------------------------------------------------
void AliPBUtils::GetSPDTrackletMult(const AliESDEvent *ESDEvent,
                                         Int_t& sum, Int_t& forwardA,
                                         Int_t& forwardC, Int_t& central)
{
	// obtain the multiplicity for eta < -.9 (forwardC), eta > 0.9 (fowardA), in
	// the central barrel (central) and its sum from SPD tracklets with the
	// AliMultiplicity class
	// code simular to PWGPP/ITS/AliAnalysisTaskSPD.cxx

	sum = forwardA = forwardC = central = 0; // initialize values

	const AliMultiplicity *mult = ESDEvent->GetMultiplicity();
	for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets();
	     iTracklet++) {
		float_t eta = mult->GetEta(iTracklet);
		if (eta < -0.9) {
			forwardC++;
		}
		else if (eta < 0.9) {
			central++;
		}
		else {
			forwardA++;
		}
		sum++;
	}
}


//------------------------------------------------------------------------------
// AliPBUtils::CutEvent using ESD
Bool_t AliPBUtils::CutEvent(const AliESDEvent *ESDEvent, TH1 *hspd,
                                 TH1 *hpriv, TH1 *hpriVtxPos, TH1 *hpriVtxDist,
                                 TH2 *hfo, TH1* hfochans, Int_t &kfo,
                                 Int_t &nip, Int_t &nop, TH1 *hpriVtxX,
                                 TH1 *hpriVtxY, TH1 *hpriVtxZ)
{
	//
	// CutEvent
	//

	AliTriggerAnalysis triggerAnalysis;
	/*
	//printf("AliPBUtils active triggers: %s\n",
	//       ESDEvent->GetHeader()->GetActiveTriggerInputs().Data());
	//printf("AliPBUtils fired triggers: %s\n",
	//       ESDEvent->GetHeader()->GetFiredTriggerInputs().Data());
	//*/
	//http://alisoft.cern.ch/viewvc/trunk/ANALYSIS/AliTriggerAnalysis.cxx?view=markup&root=AliRoot
	//Int_t AliTriggerAnalysis::SPDFiredChips(const AliESDEvent* aEsd, Int_t origin, Bool_t fillHists, Int_t layer)
	// returns the number of fired chips in the SPD
	// origin = 0 --> aEsd->GetMultiplicity()->GetNumberOfFiredChips() (filled from clusters)
	// origin = 1 --> aEsd->GetMultiplicity()->TestFastOrFiredChips() (from hardware bits)


	// SPD FastOR cut
  // return if no chips have fired 
	const Int_t fastORhw = triggerAnalysis.SPDFiredChips(ESDEvent, 1);
	if (hspd) hspd->Fill(fastORhw);
	if(fastORhw<1)
		return kFALSE;

	if (hfochans) {
		const AliMultiplicity *mult = ESDEvent->GetMultiplicity();

		for (Int_t iChipKey=0; iChipKey < 1200; iChipKey++) {
			if (mult->TestFastOrFiredChips(iChipKey)) {
				hfochans->Fill((Double_t)iChipKey);
			}
		}
	}

	if (kfo) {              // spd trigger study
		Int_t nfoctr[10];
		GetNFO(ESDEvent, "[1.0]", nfoctr, 0x0, 0x0);
		nip = nfoctr[kInnerPixel];
		nop = nfoctr[kOuterPixel];

		if(AliPBBase::GetGapBin("V0", GetV0(ESDEvent)) == AliPBBase::kBinDG) {
			hfo->Fill(nip, nop);
		}

		if(nop<2)
			kfo = 1;
		else
			kfo = nop;

		if(kfo>=10)
			kfo=9;
	}

	// collision vertex cut
	// A cut in XY is implicitly done during the reconstruction by constraining
	// the vertex to the beam diamond.

	// Primary vertex
	Bool_t kpr0 = kTRUE;
	const AliESDVertex *vertex = ESDEvent->GetPrimaryVertexTracks();
	if (!vertex) return kFALSE;
  
  if (vertex->GetNContributors()<1) {
		// SPD vertex
		vertex = ESDEvent->GetPrimaryVertexSPD();
		if(vertex->GetNContributors()<1) {
			// NO VERTEX, SKIP EVENT
			kpr0 = kFALSE;
		}
	}
	//printf("Number of vertex contributions: %i %i\n",
  //  vertex->GetNContributors(),ESDEvent->GetPrimaryVertexSPD());
  //printf("z-position of vertex: %f\n",vertex->GetZ());

	// 10 is the common value, unit: cm
  // use here 15 cm to allow systematic cut studies
  const Bool_t kpriv = kpr0 && (fabs(vertex->GetZ()) < 15.);
	if (hpriv) hpriv->Fill(kpriv);
	if (hpriVtxDist) hpriVtxDist->Fill(vertex->GetZ());
	if (hpriVtxPos && kpr0) hpriVtxPos->Fill(!kpriv);
	if(!kpriv) return kFALSE;

	if(hpriVtxX) hpriVtxX->Fill(vertex->GetX());
	if(hpriVtxY) hpriVtxY->Fill(vertex->GetY());
	if(hpriVtxZ) hpriVtxZ->Fill(vertex->GetZ());
	return kTRUE;

}

//------------------------------------------------------------------------------
// AliPBUtils::CutEvent using AOD
Bool_t AliPBUtils::CutEvent(const AliAODEvent *AODEvent, TH1 *hpriv,
                                 TH1 *hpriVtxX, TH1 *hpriVtxY, TH1 *hpriVtxZ,
                                 TH1 *hpriVtxPos, TH1 *hpriVtxDist)
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
	const Bool_t kpriv = kpr0 && (fabs(vertex->GetZ())<10.);
	// 10 is the common value, unit: cm
	hpriv->Fill(kpriv);
	if (hpriVtxDist) hpriVtxDist->Fill(vertex->GetZ());
	if (hpriVtxPos && kpr0) hpriVtxPos->Fill(!kpriv);
	if(!kpriv) return kFALSE;

	if(hpriVtxX) hpriVtxX->Fill(vertex->GetX());
	if(hpriVtxY) hpriVtxY->Fill(vertex->GetY());
	if(hpriVtxZ) hpriVtxZ->Fill(vertex->GetZ());
	return kTRUE;
}


//------------------------------------------------------------------------------
void AliPBUtils::DoVZEROStudy(const AliESDEvent *ESDEvent,
                                   TObjArray* hists, Int_t run)
{
	//
	//
	// IMPORTANT: the order of the histograms here and in
	// AliPBBase::GetHistVZEROStudies(..) has to match

	const AliESDVZERO* esdV0 = ESDEvent->GetVZEROData();
	if (!esdV0) {
		Printf("ERROR: esd V0  not available");
		return;
	}

	// determine trigger decision
	AliTriggerAnalysis triggerAnalysis;
	//const Bool_t khw = kFALSE;

	//Float_t v0a = (triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kASide, khw) ==
	//             AliTriggerAnalysis::kV0BB) ? 1. : 0.;
	//Float_t v0c = (triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kCSide, khw) ==
	//             AliTriggerAnalysis::kV0BB) ? 1. : 0.;

	// obtain OCDB objects
	AliCDBManager *man = AliCDBManager::Instance();
	TString cdbpath;
	if (man->IsDefaultStorageSet()) {
		const AliCDBStorage *dsto = man->GetDefaultStorage();
		cdbpath = TString(dsto->GetBaseFolder());
	}
	else { //should not be used!
		man->SetDefaultStorage(gSystem->Getenv("TRAIN_CDB_PATH"));
		cdbpath = TString(gSystem->Getenv("TRAIN_CDB_PATH"));
	}
	man->SetSpecificStorage("VZERO/Trigger/Data",cdbpath);
	man->SetSpecificStorage("VZERO/Calib/Data",cdbpath);
	man->SetRun(run);

	AliCDBEntry *ent1 = man->Get("VZERO/Trigger/Data");
	if (!ent1) {
		printf("AliPBUtils failed loading VZERO trigger entry from OCDB\n");
		return;
	}
	AliVZEROTriggerData *trigData = (AliVZEROTriggerData*)ent1->GetObject();
	if (!trigData) {
		printf("AliPBUtils failed loading VZERO trigger data from OCDB\n");
		return;
	}

	AliCDBEntry *ent2 = man->Get("VZERO/Calib/Data");
	if (!ent2) {
		printf("AliPBUtils failed loading VZERO calib entry from OCDB\n");
		return;
	}
	AliVZEROCalibData *calData = (AliVZEROCalibData*)ent2->GetObject();
	if (!calData) {
		printf("AliPBUtils failed loading VZERO calibration data from OCDB\n");
		return;
	}

	// fill histograms
	Int_t pmtHist = 0;
	Int_t iPMT = 0;
	for (iPMT = 0; iPMT < 64; ++iPMT) {
		Int_t board   = AliVZEROCalibData::GetBoardNumber(iPMT);
		Int_t channel = AliVZEROCalibData::GetFEEChannelNumber(iPMT);
		((TH2*)hists->At(iPMT+pmtHist))->Fill(esdV0->GetAdc(iPMT),
		                                      trigData->GetPedestalCut(0, board, channel));
	//                                      calData->GetDiscriThr(iPMT));
	}
	pmtHist = iPMT;
	for (iPMT = 0; iPMT < 64; ++iPMT) {
		((TH2*)hists->At(iPMT+pmtHist))->Fill(esdV0->GetAdc(iPMT),
		                                      esdV0->GetMultiplicity(iPMT));
	}
	/*
	// not used in 2010 for pp
	((TH2*)hists->At(pmtHist++))->Fill((Float_t)esdV0->GetTriggerChargeA(), v0a);
	((TH2*)hists->At(pmtHist++))->Fill((Float_t)esdV0->GetTriggerChargeC(), v0c);
	*/
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetGapConfig(const AliESDEvent *ESDEvent,
      TH2 *hitMapSPDinner, TH2 *hitMapSPDouter, TH2 *hitMapSPDtrklt,
      TH2 *hitMapFMDa, TH2 *hitMapFMDc, TH1 **fmdSums,
      TH1 *mulADA, TH1 *mulADC, TH1 *timeADA, TH1 *timeADC,
      TH2 *TPCGapDCAaSide, TH2 *TPCGapDCAcSide)
{
	//
	// GetGapConfigAndTracks
	//
	// retrieves the gap configuration of an event and returns it as
	// a bit vector
	// kBaseLine ensures, that this event is valid
	// + is equivalent to | in this case
  //
  // only trigger information is used
	Int_t gapConfig = AliPBBase::kBitBaseLine
    + GetV0(ESDEvent)
		+ GetFMD(ESDEvent, hitMapFMDa, hitMapFMDc, fmdSums)
    + GetAD(ESDEvent, mulADA, mulADC, timeADA, timeADC)
		+ GetSPD(ESDEvent, hitMapSPDinner, hitMapSPDouter, hitMapSPDtrklt);
    
	if (gapConfig == AliPBBase::kBitBaseLine) {
		gapConfig += GetTPC(ESDEvent, TPCGapDCAaSide, TPCGapDCAcSide);
	}
	else {
		gapConfig += GetTPC(ESDEvent, 0x0, 0x0);
	}
	if (GetFastORmultiplicity(ESDEvent) > 0) {
		gapConfig += AliPBBase::kBitCentAct;
	}
	
  return gapConfig + GetZDC(ESDEvent);

}


//------------------------------------------------------------------------------
void AliPBUtils::FillEtaPhiMap(const AliVEvent *event,
                               const AliPBTracks* tracks, TH2 *map,
                               TH2 *map_c)
{
	//
	// Fills the eta phi information about all tracks and about the tracks surving
	// the track cuts (CutTrack) into two separate histograms
	//

	// all tracks
	for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); ++itrk) {
		if (AliVParticle *trk = event->GetTrack(itrk)) {
			if (map) {
				map->Fill(trk->Eta(), trk->Phi());
			}
		}
	}

	// tracks that survived the cuts
	for (Int_t itrk = 0; itrk < tracks->GetCombinedTracks(); ++itrk) {
		if (map_c) {
			map_c->Fill(tracks->GetTrack(itrk)->Eta(), tracks->GetTrack(itrk)->Phi());
		}
	}
}


//------------------------------------------------------------------------------
void AliPBUtils::GetMultFMD(const AliESDEvent *ESDEvent, Int_t& fmdA,
                                 Int_t& fmdC, Float_t  *fmdSums)
{
	//
	// Multiplicity seen by FMD
	//
	AliTriggerAnalysis triggerAnalysis;
	triggerAnalysis.SetFMDThreshold(0.3, 0.5); // parameters got from FMD
  fmdA = triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kFMDA);
  fmdC = triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kFMDC);

	if (fmdSums) {
		const AliESDFMD* fmdData =
			(const_cast<AliESDEvent*>(ESDEvent))->GetFMDData();
		for (UShort_t det=1; det<=3; det++) {
			Int_t nRings = (det == 1 ? 1 : 2);
			for (UShort_t ir = 0; ir < nRings; ir++) {
				Char_t   ring = (ir == 0 ? 'I' : 'O');
				UShort_t nsec = (ir == 0 ? 20  : 40);
				UShort_t nstr = (ir == 0 ? 512 : 256);
				for (UShort_t sec =0; sec < nsec;  sec++) {
					for (UShort_t strip = 0; strip < nstr; strip++) {
						Float_t mult = fmdData->Multiplicity(det,ring,sec,strip);
						if (mult == AliESDFMD::kInvalidMult) continue;

						if (det == 1 && ring == 'I') fmdSums[0] += mult;
						else if (det == 2 && ring == 'I') fmdSums[1] += mult;
						else if (det == 2 && ring == 'O') fmdSums[2] += mult;
						else if (det == 3 && ring == 'I') fmdSums[3] += mult;
						else if (det == 3 && ring == 'O') fmdSums[4] += mult;
					}
				}
			}
		}
	}
}


//------------------------------------------------------------------------------
void AliPBUtils::GetMultSPD(const AliESDEvent * ESDEvent, Int_t& spdIA,
                                 Int_t& spdIC, Int_t& spdOA, Int_t& spdOC)
{
	//
	// Retrieves the multiplicity seen by the SPD FastOR
	//

	Int_t nfoctr[10];
	GetNFO(ESDEvent, "]0.9[", nfoctr, 0x0, 0x0);

	spdIA = nfoctr[kIPA];
	spdIC = nfoctr[kIPC];
	spdOA = nfoctr[kOPA];
	spdOC = nfoctr[kOPC];
}

//------------------------------------------------------------------------------
Int_t AliPBUtils::GetV0(const AliESDEvent * ESDEvent)
{
	//
	// GetV0
	//
  // AliPBBase::kBitV0A = (1 << 1) = 2
  // AliPBBase::kBitV0C = (1 << 2) = 4
	// -> V0AND = 6
	// -> V0OR  = [2,4]
  AliTriggerAnalysis triggerAnalysis;
	const Bool_t khw = kFALSE;
	const Bool_t v0A =
		(triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kASide, khw) ==
		 AliTriggerAnalysis::kV0BB);
	const Bool_t v0C =
		(triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kCSide, khw) ==
		 AliTriggerAnalysis::kV0BB);

	// the output can be 0,2,4, or 6!
  Int_t res = v0A * AliPBBase::kBitV0A + v0C * AliPBBase::kBitV0C;
  // printf("<I - AliPBUtils::GetV0> output: %i\n",res);
  
  return res;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetFMD(const AliESDEvent *ESDEvent, TH2 *hitMapFMDa,
                              TH2 *hitMapFMDc, TH1 **fmdSums)
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

	//printf("<I - GetFMD> online : %i (A) : %i (C) \n",fmdA,fmdC);
	//printf("<I - GetFMD> offline: %i (A) : %i (C) \n",
  //  triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kFMDA),
  //  triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kFMDC) );

	// preparations for a charge summation algorithm
	Bool_t hitMaps = (Bool_t)(hitMapFMDa && hitMapFMDc);
	Bool_t calcSum = fmdSums ?
		(fmdSums[0] && fmdSums[1] && fmdSums[2] && fmdSums[3] && fmdSums[4]) :
		kFALSE; // are the histograms defined?
	Float_t sum[] = { 0., 0., 0., 0., 0. };
	// summed multiplicity in the FMD detectors

	// hit map generation
	if (hitMaps || calcSum) {
		const AliESDFMD* fmdData =
			(const_cast<AliESDEvent*>(ESDEvent))->GetFMDData();

		for (UShort_t det=1; det<=3;det++) {
			Int_t nRings = (det == 1 ? 1 : 2);
			for (UShort_t ir = 0; ir < nRings; ir++) {
				Char_t   ring = (ir == 0 ? 'I' : 'O');
				UShort_t nsec = (ir == 0 ? 20  : 40);
				UShort_t nstr = (ir == 0 ? 512 : 256);
				for (UShort_t sec =0; sec < nsec;  sec++) {
					for (UShort_t strip = 0; strip < nstr; strip++) {
						Float_t mult = fmdData->Multiplicity(det,ring,sec,strip);
						if (mult == AliESDFMD::kInvalidMult) continue;

						if (calcSum) {
							if (det == 1 && ring == 'I') sum[0] += mult;
							else if (det == 2 && ring == 'I') sum[1] += mult;
							else if (det == 2 && ring == 'O') sum[2] += mult;
							else if (det == 3 && ring == 'I') sum[3] += mult;
							else if (det == 3 && ring == 'O') sum[4] += mult;
						}

						if (hitMaps) { // care about hit map specific information
							const Float_t eta = fmdData->Eta(det,ring,sec,strip);
							const Float_t phi =
								fmdData->Phi(det,ring,sec,strip) / 180. * TMath::Pi();
							//printf("FR - GetFMD: %f %f %f\n", eta, phi, mult);
							if (eta != AliESDFMD::kInvalidEta) {
								if ((-3.5 < eta) && (eta < -1.5)) {
									hitMapFMDc->Fill(eta, phi, mult); // use mult as weight
								}
								else if ((1.5 < eta) && (eta < 5.5)) {
									hitMapFMDa->Fill(eta, phi, mult); // use mult as weight
								}
							}
						}
					}
				}
			}
		}
	}

	if (calcSum) {
    //printf("DEBUG -- SUM(%f,%f,%f,%f,%f)\n", sum[0], sum[1], sum[2], sum[3],
    //  sum[4]);
		for (UInt_t i = 0; i < 5; i++) { // 
			fmdSums[i]->Fill(sum[i]);
		}
	}

	return fmdA * AliPBBase::kBitFMDA + fmdC * AliPBBase::kBitFMDC;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetAD(const AliESDEvent *ESDEvent,
                        TH1 *mulADA, TH1 *mulADC,
                        TH1 *timeADA, TH1 *timeADC)
{
	//
	// GetAD
	//
	AliESDAD *fADData = ESDEvent->GetADData();
	if (!fADData) return 0;

	// multiplicities
	if (mulADA) mulADA->Fill(fADData->GetMTotADA());
  if (mulADC) mulADC->Fill(fADData->GetMTotADA());
	
	// sum of trigger charge
  if (timeADA) timeADA->Fill(fADData->GetADATime());
	if (timeADC) timeADC->Fill(fADData->GetADCTime());

	// trigger decision
	const Bool_t vADA  = fADData->GetADADecision();
	const Bool_t vADC  = fADData->GetADCDecision();
  
  return (vADA * AliPBBase::kBitADA + vADC * AliPBBase::kBitADC);

}

//------------------------------------------------------------------------------
Int_t AliPBUtils::GetSPD(const AliESDEvent *ESDEvent, TH2 *hitMapSPDinner,
                              TH2 *hitMapSPDouter, TH2 *hitMapSPDtrklt)
{
	//
	// GetSPD
	//

	// get multiplicity from fastOR and fill corresponding hit maps
	Int_t nfoctr[10];
	GetNFO(ESDEvent, "]0.9[", nfoctr, hitMapSPDinner, hitMapSPDouter);

	// fill tracklet hit map
	if (hitMapSPDtrklt) FillSPDtrkltMap(ESDEvent, hitMapSPDtrklt);

	const Int_t ipA = nfoctr[kIPA]; // inner layer A side
	const Int_t ipC = nfoctr[kIPC]; // inner layer C side
	const Int_t opA = nfoctr[kOPA]; // outer layer A side
	const Int_t opC = nfoctr[kOPC]; // outer layer C side

	const Bool_t spdA = ipA + opA; // A side hit?
	const Bool_t spdC = ipC + opC; // C side hit?

	return spdA * AliPBBase::kBitSPDA + spdC * AliPBBase::kBitSPDC;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetTPC(const AliESDEvent * ESDEvent, TH2 *TPCGapDCAaSide,
                              TH2 *TPCGapDCAcSide)
{
	//
	//GetTPC
	//

	const Double_t etacut = 0.9;
	Int_t nA = 0;
	Int_t nC = 0;

	AliESDtrackCuts cuts;
	cuts.SetMaxDCAToVertexXY(0.1);
	cuts.SetMaxDCAToVertexZ(2.);

	for(Int_t itrack = 0; itrack < ESDEvent->GetNumberOfTracks(); itrack++){
		const AliESDtrack* esdtrack = ESDEvent->GetTrack(itrack);
		Float_t b[2];
		Float_t bCov[3];
		esdtrack->GetImpactParameters(b,bCov);

		if (bCov[0]<=0 || bCov[2]<=0) {
			printf("AliPBUtils - Estimated b resolution lower or equal zero!\n");
			bCov[0]=0;
			bCov[2]=0;
		}

		Float_t dcaToVertexXY = b[0];
		Float_t dcaToVertexZ = b[1];
		if (esdtrack->Eta() > etacut) {
			if (cuts.AcceptTrack(esdtrack)) nA++;
			if (TPCGapDCAaSide) {
				TPCGapDCAaSide->Fill(dcaToVertexXY, dcaToVertexZ);
			}
		}
		else if (esdtrack->Eta() < -etacut) {
			if (cuts.AcceptTrack(esdtrack)) nC++;
			if (TPCGapDCAcSide) {
				TPCGapDCAcSide->Fill(dcaToVertexXY, dcaToVertexZ);
			}
		}
	}

	const Bool_t tpcA = nA;
	const Bool_t tpcC = nC;

	return tpcA * AliPBBase::kBitTPCA + tpcC * AliPBBase::kBitTPCC;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::GetZDC(const AliESDEvent * ESDEvent)
{
	//
	//GetZDC
	//

	const Int_t qa = ESDEvent->GetESDZDC()->GetESDQuality();
	Bool_t zdcA = kFALSE, zdcC = kFALSE;
	for(Int_t ii=0; ii<6; ii++){
		if(qa & (1<<ii)){
			if(ii<4) zdcA = kTRUE;
			else zdcC = kTRUE;
		}
	}

	return zdcA * AliPBBase::kBitZDCA + zdcC * AliPBBase::kBitZDCC;
}


//==============================================================================
//------------------------------------------------------------------------------
#ifdef STD_ALIROOT
Int_t AliPBUtils::FMDHitCombinations(const AliESDEvent* ESDEvent,
                                          Int_t side)
{
	//
	// copy of the FMDHitCombinations function originating from AliTriggerAnalysis
	//
	// side == 0 -> A side, side > 0 -> C side

	// workaround for AliESDEvent::GetFMDData is not const!
	const AliESDFMD* fmdData =
		(const_cast<AliESDEvent*>(ESDEvent))->GetFMDData();
	if (!fmdData)
		{
			puts("AliESDFMD not available");
			return -1;
		}

	Int_t detFrom = (side == 0) ? 1 : 3;
	Int_t detTo   = (side == 0) ? 2 : 3;

	Int_t triggers = 0;
	Float_t totalMult = 0;
	for (UShort_t det=detFrom;det<=detTo;det++) {
		Int_t nRings = (det == 1 ? 1 : 2);
		for (UShort_t ir = 0; ir < nRings; ir++) {
			Char_t   ring = (ir == 0 ? 'I' : 'O');
			UShort_t nsec = (ir == 0 ? 20  : 40);
			UShort_t nstr = (ir == 0 ? 512 : 256);
			for (UShort_t sec =0; sec < nsec;  sec++) {
				for (UShort_t strip = 0; strip < nstr; strip++) {
					Float_t mult = fmdData->Multiplicity(det,ring,sec,strip);
					if (mult == AliESDFMD::kInvalidMult) continue;

					if (mult > 0.3) // fFMDLowCut
						totalMult = totalMult + mult;
					else {
						if (totalMult > 0.5) // fFMDHitCut
							triggers++;
					}
				}
			}
		}
	}
	return triggers;
}
#endif


//==============================================================================
//------------------------------------------------------------------------------
void AliPBUtils::SPDLoadGeom(Int_t run)
{
	// method to get the gGeomanager
	// it is called at the CreatedOutputObject stage
	// to comply with the CAF environment

	AliCDBManager *man = AliCDBManager::Instance();

	TString cdbpath;
	if (man->IsDefaultStorageSet()) {
		const AliCDBStorage *dsto = man->GetDefaultStorage();
		cdbpath = TString(dsto->GetBaseFolder());
		//printf("default was set to: %s\n", cdbpath.Data());
	}
	else { //should not be used!
		// man->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
		// would be needed on grid
		man->SetDefaultStorage(gSystem->Getenv("TRAIN_CDB_PATH"));
		cdbpath = TString(gSystem->Getenv("TRAIN_CDB_PATH"));
	}

	man->SetSpecificStorage("ITS/Align/Data",cdbpath);
	man->SetSpecificStorage("GRP/Geometry/Data",cdbpath);
	man->SetRun(run);

	AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
	if (!obj) {
		printf("AliPBUtils failed loading geometry object\n");
		return;
	}
	AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
	AliGeomManager::ApplyAlignObjsFromCDB("ITS");
}

//------------------------------------------------------------------------------
Bool_t AliPBUtils::SPDLoc2Glo(Int_t id, const Double_t *loc,
                                   Double_t *glo)
{
	//
	//SPDLoc2Glo, do not touch
	//

	static TGeoHMatrix mat;
	Int_t vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
	if (vid<0) {
		printf("AliPBUtils Did not find module with such ID %d\n",id);
		return kFALSE;
	}
	AliITSAlignMille2Module::SensVolMatrix(vid,&mat);
	mat.LocalToMaster(loc,glo);
	return kTRUE;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::CheckChipEta(Int_t chipKey, TString scut,
                                    const Double_t vtxPos[],
                                    TH2 *hitMapSPDinner, TH2 *hitMapSPDouter)
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

		if (chipKey < 400) { // inner SPD layer
			if (hitMapSPDinner) {
				Double_t phi = pos.Phi(); // output in the range -Pi +Pi
				if (phi < 0.) phi += TMath::TwoPi(); // remap to the interval [0, TwoPi)
				const Double_t eta = pos.Eta();
				hitMapSPDinner->Fill(eta, phi);
			}
		}
		else {
			if (hitMapSPDouter) { // outer SPD layer
				Double_t phi = pos.Phi(); // output in the range -Pi +Pi
				if (phi < 0.) phi += TMath::TwoPi(); // remap to the interval [0, TwoPi)
				const Double_t eta = pos.Eta();
				hitMapSPDouter->Fill(eta, phi);
			}
		}

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
void AliPBUtils::GetNFO(const AliESDEvent *ESDEvent, TString etacut,
                             Int_t ctr[], TH2 *hitMapSPDinner,
                             TH2 *hitMapSPDouter)
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
			const Int_t iseta = CheckChipEta(iChipKey, etacut, vtxPos, hitMapSPDinner,
			                                 hitMapSPDouter);
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


//--------------------------------------------------------------------------
Int_t AliPBUtils::GetFastORmultiplicity(const AliESDEvent* ESDEvent)
{
	// determine the number of fired fastOR chips in both layers within
	// -0.9 < eta < 0.9
	//

	const AliMultiplicity *mult = ESDEvent->GetMultiplicity();

	// position of the primary vertex
	Double_t tmp[3] = { 0., 0., 0. };
	ESDEvent->GetPrimaryVertex()->GetXYZ(tmp);
	Double_t vtxPos[3] = { tmp[0], tmp[1], tmp[2] };

	Int_t multiplicity = 0;

	for (Int_t iChipKey=0; iChipKey < 1200; iChipKey++) {
		if(mult->TestFastOrFiredChips(iChipKey)){
			// here you check if the FastOr bit is 1 or 0
			const Int_t iseta = CheckChipEta(iChipKey, "[0.9]", vtxPos, 0x0, 0x0);
			if(iseta==0)
				continue;
			else
				++multiplicity;
		}
	}

	return multiplicity;
}


//--------------------------------------------------------------------------
void AliPBUtils::FillSPDtrkltMap(const AliVEvent* event,
                                      TH2 *hitMapSPDtrklt)
{
	//
	// fill eta phi map of SPD tracklets
	//

	if (hitMapSPDtrklt) {
		if (TString(event->ClassName()) == "AliESDEvent") {
			const AliMultiplicity *mult = ((AliESDEvent*)event)->GetMultiplicity();
			for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets();
			     iTracklet++) {
				Double_t eta = mult->GetEta(iTracklet);
				Double_t phi = mult->GetPhi(iTracklet);
				hitMapSPDtrklt->Fill(eta, phi);
			}
		}
		else if (TString(event->ClassName()) == "AliAODEvent") {
			const AliAODTracklets* mult = ((AliAODEvent*)event)->GetTracklets();
			for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets();
			     iTracklet++) {
				Double_t eta = -TMath::Log(TMath::Tan(mult->GetTheta(iTracklet)/2.));
				Double_t phi = mult->GetPhi(iTracklet);
				hitMapSPDtrklt->Fill(eta, phi);
			}
		}
	}
}


//==========================================================================
Int_t AliPBUtils::GetPID(AliPIDResponse *pid, const AliVTrack *trk,
                              Int_t mode /* = 0 */)
{
	// determines PID for ESDs and AODs
	// 
	//

	if (!pid) return AliPBBase::kBinPIDUnknown; // no pid available

	Double_t tpcpion = -999., tpckaon = -999., tpcproton = -999.,
		tpcelectron = -999.;
	Double_t tofpion = -999., tofkaon = -999., tofproton = -999.,
		tofelectron = -999.;
	Double_t itspion = -999., itskaon = -999., itsproton = -999.,
		itselectron = -999.;


	// check whether the track was pure ITS standalone track (has not left TPC)
	const Bool_t kits = !(trk->GetStatus() & AliESDtrack::kTPCout);

	if (kits) { // do ITS pid
		const Double_t nin = 3.;

		itspion     = pid->NumberOfSigmasITS( trk, AliPID::kPion );
		itskaon     = pid->NumberOfSigmasITS( trk, AliPID::kKaon );
		itsproton   = pid->NumberOfSigmasITS( trk, AliPID::kProton );
		itselectron = pid->NumberOfSigmasITS( trk, AliPID::kElectron );

		// printf("PIDa: pi,K,p,e - %e %e %e %e\n",
    //  itspion,itskaon,itsproton,itselectron);
      
    if(fabs(itspion) < nin) // inclusion only
			return AliPBBase::kBinPion;
		else if(fabs(itskaon) < nin)
			return AliPBBase::kBinKaon;
		else if(fabs(itsproton) < nin)
			return AliPBBase::kBinProton;
		else if(fabs(itselectron) < nin)
			return AliPBBase::kBinElectron;
		else{
			return AliPBBase::kBinPIDUnknown;
		}
	}

	tpcpion     = pid->NumberOfSigmasTPC( trk, AliPID::kPion );
	tpckaon     = pid->NumberOfSigmasTPC( trk, AliPID::kKaon );
	tpcproton   = pid->NumberOfSigmasTPC( trk, AliPID::kProton );
	tpcelectron = pid->NumberOfSigmasTPC( trk, AliPID::kElectron );
  
      
	// derive information, whether tof pid is available
	const Bool_t ka = !(trk->GetStatus() & AliESDtrack::kTOFmismatch); 
	const Bool_t kb =  (trk->GetStatus() & AliESDtrack::kTOFpid);
	const Bool_t ktof = ka && kb;

	// retrieve TOF information if available
	if(ktof){
		tofpion     = pid->NumberOfSigmasTOF(trk, AliPID::kPion);
		tofkaon     = pid->NumberOfSigmasTOF(trk, AliPID::kKaon);
		tofproton   = pid->NumberOfSigmasTOF(trk, AliPID::kProton);
		tofelectron = pid->NumberOfSigmasTOF(trk, AliPID::kElectron);
    
      
	}

	if (mode == 0) {
		// 2010 cuts for PID
		//printf("Getting PID with 2010 method ...\n");
    //printf("PID TPC: pi,K,p,e - %e %e %e %e\n",
    //  tpcpion,tpckaon,tpcproton,tpcelectron);
    //printf("PID TOF: pi,K,p,e - %e %e %e %e\n",
    //  tofpion,tofkaon,tofproton,tofelectron);
    
    const Double_t nin = 3.;

		// inclusion + exclusion (either TPC or TOF)
		// inclusive: given particle within nin sigma
		// exclusive: given particle within nin sigma and all other particles outside of nin sigma
    
    if((fabs(tpcpion) < nin && fabs(tpckaon) > nin && fabs(tpcproton) > nin
		    && fabs(tpcelectron) > nin) ||
		   (fabs(tofpion) < nin && fabs(tofkaon) > nin && fabs(tofproton) > nin
		    && fabs(tofelectron) > nin))
			return AliPBBase::kBinPionE;
      
		else if((fabs(tpcpion) > nin && fabs(tpckaon) < nin && fabs(tpcproton) > nin
		         && fabs(tpcelectron) > nin) ||
		        (fabs(tofpion) > nin && fabs(tofkaon) < nin && fabs(tofproton) > nin
		         && fabs(tofelectron) > nin))
			return AliPBBase::kBinKaonE;
      
		else if((fabs(tpcpion) > nin && fabs(tpckaon) > nin && fabs(tpcproton) < nin
		         && fabs(tpcelectron) > nin) ||
		        (fabs(tofpion) > nin && fabs(tofkaon) > nin && fabs(tofproton) < nin
		         && fabs(tofelectron) > nin))
			return AliPBBase::kBinProtonE;
      
		else if((fabs(tpcpion) > nin && fabs(tpckaon) > nin && fabs(tpcproton) > nin
		         && fabs(tpcelectron) < nin) ||
		        (fabs(tofpion) > nin && fabs(tofkaon) > nin && fabs(tofproton) > nin
		         && fabs(tofelectron) < nin))
			return AliPBBase::kBinElectronE;
      
		else if(fabs(tpcpion) < nin || fabs(tofpion) < nin) // inclusion (TPC + TOF)
			return AliPBBase::kBinPion;
      
		else if(fabs(tpckaon) < nin || fabs(tofkaon) < nin)
			return AliPBBase::kBinKaon;
      
		else if(fabs(tpcproton) < nin || fabs(tofproton) < nin)
			return AliPBBase::kBinProton;
      
		else if(fabs(tpcelectron) < nin || fabs(tofelectron) < nin)
			return AliPBBase::kBinElectron;
      
		else{
			return AliPBBase::kBinPIDUnknown;
		}
	}
	else if (mode == 1) {
		// printf("Getting PID with 2010 method ...\n");
		// 2011 cuts for PID in LHC11f
		// TPC: [-3,5] sigma (pion)
		// TOF: 3 sigma for all,
		// only Pion is tuned!
		// ONLY INCLUSION CUTS NO EXCLUSION
		const Double_t nin = 3.;

		// use bayesian method
    
    
    // use old method
    if(tpcpion > -2. && tpcpion < 4. && fabs(tofpion) < nin )
			return AliPBBase::kBinPion;
		else if(fabs(tpckaon) < nin && fabs(tofkaon) < nin)
			return AliPBBase::kBinKaon;
		else if(fabs(tpcproton) < nin &&  fabs(tofproton) < nin)
			return AliPBBase::kBinProton;
		else if(fabs(tpcelectron) < nin &&  fabs(tofelectron) < nin)
			return AliPBBase::kBinElectron;
		else{
			return AliPBBase::kBinPIDUnknown;
		}
    
	}
	return AliPBBase::kBinPIDUnknown;
}


//==========================================================================
Int_t AliPBUtils::GetPID(Int_t pdgCode)
{
	//
	// determine particle type based on PDG code
	//

	if (TMath::Abs(pdgCode) == 211) return AliPBBase::kBinPionE;
	else if (TMath::Abs(pdgCode) == 321) return AliPBBase::kBinKaonE;
	else if (TMath::Abs(pdgCode) == 2212) return AliPBBase::kBinProtonE;
	else if (TMath::Abs(pdgCode) == 11) return AliPBBase::kBinElectronE;
	else return AliPBBase::kBinPIDUnknown;
}


//------------------------------------------------------------------------------
Int_t AliPBUtils::CombinePID(const Int_t pid[])
{
	//
	// combine the PID result
	//

	// determine return value
	if (pid[0] == pid[1]) { // same result for both tracks
		return pid[0];
	}
	// one track identified with exclusion the other only without
	else if ((pid[0] == AliPBBase::kBinPionE &&
	          pid[1] == AliPBBase::kBinPion) ||
	         (pid[1] == AliPBBase::kBinPionE &&
	          pid[0] == AliPBBase::kBinPion)) {
		return AliPBBase::kBinPion;
	}
	else if ((pid[0] == AliPBBase::kBinKaonE &&
	          pid[1] == AliPBBase::kBinKaon) ||
	         (pid[1] == AliPBBase::kBinKaonE &&
	          pid[0] == AliPBBase::kBinKaon)) {
		return AliPBBase::kBinKaon;
	}
	else if ((pid[0] == AliPBBase::kBinProtonE &&
	          pid[1] == AliPBBase::kBinProton) ||
	         (pid[1] == AliPBBase::kBinProtonE &&
	          pid[0] == AliPBBase::kBinProton)) {
		return AliPBBase::kBinProton;
	}
	else if ((pid[0] == AliPBBase::kBinElectronE &&
	          pid[1] == AliPBBase::kBinElectron) ||
	         (pid[1] == AliPBBase::kBinElectronE &&
	          pid[0] == AliPBBase::kBinElectron)) {
		return AliPBBase::kBinElectron;
	}
	// one track identified and one not
	else if (((pid[0] == AliPBBase::kBinPionE ||
	           pid[0] == AliPBBase::kBinPion) &&
	          pid[1] == AliPBBase::kBinPIDUnknown) ||
	         (pid[1] == AliPBBase::kBinPion &&
	          pid[0] == AliPBBase::kBinPIDUnknown)) {
		return AliPBBase::kBinSinglePion;
	}
	else if (((pid[0] == AliPBBase::kBinKaonE ||
	           pid[0] == AliPBBase::kBinKaon)&&
	          pid[1] == AliPBBase::kBinPIDUnknown) ||
	         (pid[1] == AliPBBase::kBinKaon &&
	          pid[0] == AliPBBase::kBinPIDUnknown)) {
		return AliPBBase::kBinSingleKaon;
	}
	else if (((pid[0] == AliPBBase::kBinProtonE ||
	           pid[0] == AliPBBase::kBinProton) &&
	          pid[1] == AliPBBase::kBinPIDUnknown) ||
	         (pid[1] == AliPBBase::kBinProton &&
	          pid[0] == AliPBBase::kBinPIDUnknown)) {
		return AliPBBase::kBinSingleProton;
	}
	else if (((pid[0] == AliPBBase::kBinElectronE ||
	           pid[0] == AliPBBase::kBinElectron) &&
	          pid[1] == AliPBBase::kBinPIDUnknown) ||
	         (pid[1] == AliPBBase::kBinElectron &&
	          pid[0] == AliPBBase::kBinPIDUnknown)) {
		return AliPBBase::kBinSingleElectron;
	}
	else
		return AliPBBase::kBinPIDUnknown;
}


//------------------------------------------------------------------------------
TLorentzVector AliPBUtils::GetKinematics(const Double_t *pa,
                                              const Double_t *pb,
                                              Double_t ma,
                                              Double_t mb, Double_t& cts)
{
	//
	//get kinematics, cts = cos(theta_{#star})
	//

	TLorentzVector va, vb;
	va.SetXYZM(pa[0], pa[1], pa[2], ma);
	vb.SetXYZM(pb[0], pb[1], pb[2], mb);
	const TLorentzVector sumv = va+vb;
  //printf("IVM: %f\n",sumv.Mag());

	const TVector3 bv = -sumv.BoostVector();

	va.Boost(bv);
	vb.Boost(bv);

	// 3-vectors in the restframe of the mother particle
	const TVector3 pra = va.Vect();
	const TVector3 prb = vb.Vect();
	const TVector3 diff = pra - prb; // their difference
  //printf("Transformed into rest frame of mother ...\n");
  //printf("Transformed variables: %f %f %f %f %f\n",
  //  diff.Mag2(),pra.Mag2(),prb.Mag2(),pra.Mag(),prb.Mag());
	
  cts = 0.;
  if (pra.Mag()*prb.Mag() > 0) {
    cts = (diff.Mag2()-pra.Mag2()-prb.Mag2()) / (-2.*pra.Mag()*prb.Mag());
  }
  //printf("cos theta star: %f\n",cts);
	// cosine theta star, calculated according to the law-of-cosines

	return sumv;
}


//------------------------------------------------------------------------------
Double_t AliPBUtils::GetOA(const Double_t *pa, const Double_t *pb)
{
	//
	//cosOpeningAngle
	//

	TVector3 va, vb;
	va.SetXYZ(pa[0], pa[1], pa[2]);
	vb.SetXYZ(pb[0], pb[1], pb[2]);

	const TVector3 ua = va.Unit();
	const TVector3 ub = vb.Unit();

	return ua.Dot(ub);
}
