/*************************************************************************
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
// AliCDMesonUtils
// for
// AliAnalysisTaskCDMeson
//
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>

#ifndef ALICDMESONUTILS_H
#define ALICDMESONUTILS_H

#include <THnSparse.h> // forward declaration not possible


#define STD_ALIROOT
// if this is defined the code should run with a standard aliroot


class TH1;
class TH2;
class TLorentzVector;
class TParticle;
class TVector3;
class TObjArray;

class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVTrack;
class AliPIDResponse;

class AliCDMesonTracks;

class AliCDMesonUtils
{
public:
	enum{
		kInnerPixel = 0,
		kOuterPixel,
		kIPA,
		kIPC,
		kOPA,
		kOPC
	};

	// ESD only
	//---------

	// cuts for ESD analysis
	static Bool_t CutEvent(const AliESDEvent *ESDEvent, TH1 *hspd, TH1 *hpriv,
	                       TH1* hpriVtxPos, TH1* hpriVtxDist, TH2 *hfo,
	                       TH1* hfochans, Int_t &kfo, Int_t &nip, Int_t &nop,
	                       TH1* hpriVtxX, TH1* hpriVtxY, TH1* hpriVtxZ);

	// V0-AND/OR
	static Bool_t V0AND(const AliESDEvent *ESDEvent) {
		// were both V0 detectors A/C firing?
		// (4 = AliCDMesonBase::kBinNG)
		return (GetV0(ESDEvent) == 4);
	}

	static Bool_t V0OR(const AliESDEvent *ESDEvent) {
		// did one of the V0 detectors fire?
		// (1-4 => single or both detectors active)
		return (GetV0(ESDEvent) > 0);
	}

	// VZERO study
	static void DoVZEROStudy(const AliESDEvent *ESDEvent, TObjArray* hists,
	                         const Int_t run);

	// gap determination
	static Int_t GetGapConfig(const AliESDEvent *ESDEvent, TH2 *hitMapSPDinner,
	                          TH2 *hitMapSPDouter, TH2 *hitMapSPDtrklt,
	                          TH2 *hitMapFMDa, TH2 *hitMapFMDc,
	                          TH1 **fmdSums, TH2 *TPCGapDCAaSide,
	                          TH2 *TPCGapDCAcSide);

	static Int_t GetFastORmultiplicity(const AliESDEvent *ESDEvent);

	static void GetSPDTrackletMult(const AliESDEvent *ESDEvent, Int_t& sum,
	                               Int_t& forwardA, Int_t& forwardC,
	                               Int_t& central);

	static void SPDLoadGeom(const Int_t run); // only needed for ESDs, not in AODs

	// functions needed for the empty event study (ESD only)
	static Int_t GetEventType(const AliESDEvent *ESDEvent);
	static void GetMultFMD(const AliESDEvent *ESDEvent, Int_t& fmdA,
	                       Int_t& fmdC, Float_t *fmdSums);
	static void GetMultSPD(const AliESDEvent *ESDEvent, Int_t& spdIA,
	                       Int_t& spdIC, Int_t& spdOA, Int_t& spdOC);


	// AOD only
	//---------
	static Bool_t CutEvent(const AliAODEvent *AODEvent, TH1 *hpriv, TH1* hpriVtxX,
	                       TH1* hpriVtxY, TH1* hpriVtxZ, TH1* hpriVtxPos,
	                       TH1* hpriVtxDist);


	// independent
	//------------
	static void FillEtaPhiMap(const AliVEvent *event,
	                          const AliCDMesonTracks* tracks, TH2 *map,
	                          TH2 *map_c);
	static void SwapTrack(const AliVTrack *trks[]);
	static Int_t GetCombCh(const AliVTrack *trks[]);
	static Int_t GetCombPID(AliPIDResponse *pid, const AliVTrack *trks[],
	                        const Int_t mode, TH2 *comb2trkPID = 0x0);
	static Int_t GetCombPID(const TParticle* particles[], TH2 *comb2trkPID = 0x0);
	static void GetMassPtCtsOA(const Int_t pid, const TVector3* momenta[],
	                           Double_t& mass, Double_t& pt, Double_t& cts,
	                           Double_t& oa);
	static void GetPWAinfo(const Int_t pid, const AliVTrack *trks[],
	                       Float_t& theta, Float_t& phi, Float_t& mass,
	                       Float_t momentum[]);
	static void FillSPDtrkltMap(const AliVEvent *event, TH2 *hitMapSPDtrklt);

private:
	// ESD only
	//---------

	// Gap determination functions
	static Int_t GetV0(const AliESDEvent *ESDEvent);
	static Int_t GetFMD(const AliESDEvent *ESDEvent, TH2 *hitMapFMDa,
	                    TH2 *hitMapFMDc, TH1** fmdSums);
	static Int_t GetSPD(const AliESDEvent *ESDEvent, TH2 *hitMapSPDinner,
	                    TH2 *hitMapSPDouter, TH2 *hitMapSPDtrklt);
	static Int_t GetTPC(const AliESDEvent *ESDEvent, TH2 *TPCGapDCAaSide,
	                    TH2 *TPCGapDCAcSide);
	static Int_t GetZDC(const AliESDEvent *ESDEvent); // not used so far

#ifdef STD_ALIROOT
	// helpers for the FMD gap determination
	static Int_t FMDHitCombinations(const AliESDEvent* ESDEvent, Int_t side);
#endif

	// helpers for the SPD gap determination
	static Bool_t SPDLoc2Glo(const Int_t id, const Double_t *loc, Double_t *glo);
	static Int_t CheckChipEta(const Int_t chipKey, const TString scut,
	                          const Double_t vtxPos[], TH2 *hitMapSPDinner,
	                          TH2 *hitMapSPDouter);
	static void GetNFO(const AliESDEvent *ESDEvent, const TString etacut,
	                   Int_t ctr[], TH2 *hitMapSPDinner, TH2 *hitMapSPDouter);

	// AOD only
	//---------


	// independent
	//----------
	static Int_t GetPID(AliPIDResponse *pid, const AliVTrack *trk,
	                    const Int_t mode = 0);
	static Int_t GetPID(const Int_t pdgCode);
	static Int_t CombinePID(const Int_t pid[]);

	static TLorentzVector GetKinematics(const Double_t *pa, const Double_t *pb,
	                                    const Double_t ma, const Double_t mb,
	                                    Double_t & cts);
	static Double_t GetOA(const Double_t *pa, const Double_t *pb);
};

#endif
