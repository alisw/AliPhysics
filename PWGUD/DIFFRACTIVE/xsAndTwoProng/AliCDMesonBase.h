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
// AliCDMesonBase
// for
// AliAnalysisTaskCDMeson
//
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>

#ifndef ALICDMESONBASE_H
#define ALICDMESONBASE_H

#include "THnSparse.h"

class TH2;
class TH2F;
class TH1F;
class TList;
class TObjArray;

class AliCDMesonBase
{
public:
	enum{
		// combined charged
		kBinPM = 1, // unlike-sign
		kBinPPMM, // like-sign

		// pid results (used for single particles and pairs)
		kBinPionE = 1, // pions with exclusion cut
		kBinPion, // pions without exclusion cut but TPC and TOF
		kBinSinglePion, // one pion identified, the other particle unidentified
		kBinKaonE, // same for kaons ...
		kBinKaon,
		kBinSingleKaon,
		kBinProtonE, // protons
		kBinProton,
		kBinSingleProton,
		kBinElectronE, // and electrons
		kBinElectron,
		kBinSingleElectron,
		kBinPIDUnknown, // could not be identified

		// gap Conditions
		kBinDG = 1, // double gap
		kBinGC, // single gap c side
		kBinGA, // single gap a side
		kBinNG, // no gap


		// event types used for the empty event studies
		kBinEventI = 1, // Beam-Beam Interaction
		kBinEventA, // Beam from A-side
		kBinEventC, // Beam from C-side
		kBinEventAC, // Beam from one side (there is no separation in 2011)
		kBinEventE, // no beam (empty event)
		kBinEventUnknown,

		// StatsFlow histogram entries
		// names for the bins are specified in the .cxx-file
		kBinTotalInput = 0,
		kBinGoodInput,
		kBinV0OR,
		kBinV0AND,
		kBinEventsAfterCuts,
		kBinEventsWithOutPileUp,
		kBinv0Gap,
		kBinv0fmdGap,
		kBinv0fmdspdGap,
		kBinv0fmdspdtpcGap,
		kBinv0fmdspdtpczdcGap,
		kBinfmdGap,
		kBinspdGap,
		kBintpcGap,
		kBintpcspdGap,
		kBintpcspdfmdGap,
		kBintpcspdfmdv0Gap,
		kBinspdfmdGap,
		kBinspdfmdv0Gap,
		kBinTwoTrackEvents,
		kBinThreeTrackEvents,
		kBinPionEvents,
		kBinKaonEvents,
		kBinProtonEvents,
		kBinElectronEvents,
		kBinUnknownPIDEvents,
		kBinResidualTracks,
		kBinResidualTracklets,
		kBinCDonlyEvents,
		kBinLastValue, // used to specify the correct histogram width


		// MC event/process types
		kBinND = 0, // used for data as well, "Non-Diffractive"
		kBinCD, // central diffractive
		kBinSD, // single diffractive
		kBinDD, // double diffractive


		// gap-condition bits used in AliAnalysisTaskCDMeson::fGapRun
		kBitBaseLine = (1<<0),

		kBitV0A  = (1<<1),
		kBitV0C  = (1<<2),
		kBitFMDA = (1<<3),
		kBitFMDC = (1<<4),

		kBitSPDA = (1<<5),
		kBitSPDC = (1<<6),
		kBitTPCA = (1<<7),
		kBitTPCC = (1<<8),

		kBitZDCA   = (1<<9),
		kBitZDCC   = (1<<10),
		kBitCentAct = (1<<11),
		kBitGapMax = (1<<12),


		// analysis task status bits
		// do not change the order in order to be backward compatible!
		kBitConfigurationSet = (1<<0), // if not set everything is active
		kBitImposeZDCGap = (1<<1), // excluding the ZDC gap
		kBitEtaPhiMaps = (1<<2),
		kBitEtaPhiMapsWithCuts = (1<<3),
		kBitStatsFlow = (1<<4),
		kBitMultPerGapHists = (1<<5),
		kBitRmMultPerGapHists = (1<<6),
		kBitTHnMother = (1<<7),
		kBitFastORStudy = (1<<8), // also not active by default
		kBitHitMapSPD = (1<<9),
		kBitHitMapFMD = (1<<10),
		kBitVtxStudies = (1<<11),
		kBitEEStudy = (1<<12),
		kBitPIDStudy = (1<<13),
		kBitMCProcess = (1<<14),
		kBitFMDsum = (1<<15),
		kBitSoftTracks = (1<<16),
		kBitPWAtree = (1<<17),
		kBitTHnMC = (1<<18),
		kBitMultResponseMC = (1<<19),
		kBitMultStudy = (1<<20),
		kBitReadPreprocessedGap = (1<<21), // not active by default as well
		kBitReduceGapEvents = (1<<22),
		kBitVZEROStudy = (1<<23),
		kBitTPCGapStudy = (1<<24),
		kBitAllTrackMass = (1<<25),
		kBitCDonly = (1<<26),
		kBitFastORmultStudy = (1<<27),
		kBitConfigurationVersion = (1<<28) // always set, last bit
	};

	static Int_t GetGapBin(const TString tag, const Int_t gapcg,
	                       Bool_t checkCentralActivity = kTRUE);

	static THnSparseD* GetThnMother(TString name = "CDMeson_Mother");
	static void FillThnMother(THnSparseD *thn, const Double_t vNch,
	                          const Double_t vCombCh, const Double_t vCombPID,
	                          const Double_t vV0, const Double_t vFMD,
	                          const Double_t vSPD, const Double_t vTPC,
	                          const Double_t vMass, const Double_t vPt,
	                          const Double_t vOA, const Double_t vCTS,
	                          const Double_t vDaughterPt,
	                          const Double_t vTrackResiduals,
	                          const Double_t vVertexZ,
	                          const Double_t vProcessType,
	                          const Double_t vVertexCoincidence,
	                          const Double_t vTrkltResiduals);
	static Int_t GetAxisMother(const TString name);

	static THnSparseD* GetThnEmptyEvents();
	static void FillThnEmptyEvents(THnSparseD *thn, const Int_t eventType,
	                               const Int_t multFMDA, const Int_t multFMDC,
	                               const Int_t multSPDIA, const Int_t multSPDIC,
	                               const Int_t multSPDOA, const Int_t multSPDOC,
	                               const Int_t multSPDtrklA,
	                               const Int_t multSPDtrklC, const Int_t fmdSum1I,
	                               const Int_t fmdSum2I, const Int_t fmdSum20,
	                               const Int_t fmdSum3I, const Int_t fmdSum3O/*,
	                               const Int_t multTPC,
	                               const Int_t multTPCdiffVertex */);
	static Int_t GetAxisEmptyEvents(const TString name);

	static THnSparseD* GetThnMultiplicity();
	static void FillThnMultiplicity(THnSparseD *thn, const Double_t vNch,
	                                const Double_t vNsoft,
	                                const Double_t vNcombined,
	                                const Double_t vV0, const Double_t vFMD,
	                                const Double_t vSPD, const Double_t vTPC,
	                                const Double_t vNresidualTracks,
	                                const Double_t vNresidualTracklets,
	                                const Double_t vVertexZ,
	                                const Double_t vVerticesDistance,
	                                const Double_t vProcessType);
	static Int_t GetAxisMultiplicity(const TString name);

	static TH1F* GetHistStatsFlow();
	static TH2F* GetHistPIDStudies(TString name);
	static TObjArray* GetHistVZEROStudies(TList* l);

	static void GetGapTriggers(THnSparseI* gaprun, Int_t gapCondition,
	                           const Int_t run, Double_t& triggers,
	                           Double_t& total);

	static void GetNoGapTriggers(THnSparseI* gaprun, Int_t gapCondition,
	                             const Int_t run, Double_t& triggers,
	                             Double_t& total);
private:
	static void CheckRange(Double_t &var, const Double_t min, const Double_t max);
	static Int_t GetAxis(TString thntit, TString name);
	static TString GetTitleMother();
	static TString GetTitleEmptyEvents();
	static TString GetTitleMultiplicity();

	//== INVARIANT MASS DISTRIBUTIONS (ThnMother) ================================
	//-- event characteristics
	// number of charged primary particles - combined => can be w/ soft tracks
	static const Int_t fgkNNcombined      = 2; // number of bins
	static const Double_t fgkMinNcombined = 2; // lower border
	static const Double_t fgkMaxNcombined = 4; // upper border (#bins + lower border)

	// unlike sign or like sign charged?
	static const Int_t fgkNCombCh      = 2; // kBinPPMM = number of bins
	static const Double_t fgkMinCombCh = 1.; // kBinPM = lower border
	static const Double_t fgkMaxCombCh = 3.; // kBinPPMM + kBinPM = upper border

	// two track PID, take care the enum changed from time to time
	static const Int_t fgkNCombPID      = 13; // kBinPIDUnknown = number of bins
	static const Double_t fgkMinCombPID = 1.; // kBinPionE = lower border
	static const Double_t fgkMaxCombPID = 14.;// ...  = upper border

	// gap configuration, used for different detectors
	static const Int_t fgkNGapConfig      = 4; // kBinNG = number of bins
	static const Double_t fgkMinGapConfig = 1.; // kBinDG = lower border
	static const Double_t fgkMaxGapConfig = 5.; // kBinNG + kBinDG = upper border

	//-- mother kinematics
	// invariant mass of the two-track system
	static const Int_t fgkNMass      = 1024; // number of bins
	static const Double_t fgkMinMass = 0.; // lower border
	static const Double_t fgkMaxMass = 5.12; // upper border

	// transverse momentum of the two-track system
	static const Int_t fgkNMotherPt      = 128;//10; // number of bins
	static const Double_t fgkMinMotherPt = 0.; // lower border
	static const Double_t fgkMaxMotherPt = .64; //6.4 //2; // upper border

	// cosine theta* (opening angle of the two daugther tracks in the
	// centre-of-mass system of the two-track/mother system)
	// **no meaning full variable**
	static const Int_t fgkNCTS      = 2; // number of bins
	static const Double_t fgkMinCTS = -1.; // lower border
	static const Double_t fgkMaxCTS = -0.9; // upper border

	// opening angle in the lab frame
	static const Int_t fgkNOA      = 20; // number of bins
	static const Double_t fgkMinOA = -1.; // lower border
	static const Double_t fgkMaxOA = 1.; // upper border

	//-- daughter kinematics
	// transverse momentum of one of the two daughter particles
	// (randomly selected)
	static const Int_t fgkNDaughterPt      = 128; // number of bins
	static const Double_t fgkMinDaughterPt = 0.; // lower border
	static const Double_t fgkMaxDaughterPt = 6.4; // upper border

	// pseudo rapidity of one of the two daughter particles
	// (randomly selected)
	//static const Int_t fgkNDaughterEta      = 64; // number of bins
	//static const Double_t fgkMinDaughterEta = -1.28; // lower border
	//static const Double_t fgkMaxDaughterEta =  1.28; // upper border

	//-- Event quality information
	// boolean values to reduce output size

	// are there tracks in addition to the ones selected using AliCDMesonTracks
	// (ITSTPC, ITSsa, ITSpureSA) (0 = no, 1 = yes)
	static const Int_t fgkNTrackResiduals      = 2; // number of bins
	static const Double_t fgkMinTrackResiduals = 0.; // lower border
	static const Double_t fgkMaxTrackResiduals = 2.; // upper border

	// vertex with in +/-4cm (0 = no, 1 = yes)
	static const Int_t fgkNVertexZinRng      = 2; // number of bins
	static const Double_t fgkMinVertexZinRng = 0.; // lower border
	static const Double_t fgkMaxVertexZinRng = 2.; // upper border

	// are the vertices from SPD and tracks within 0.5cm? (0 = no, 1 = yes)
	static const Int_t fgkNVertexCoincidence      = 2; // number of bins
	static const Double_t fgkMinVertexCoincidence = 0.; // lower border
	static const Double_t fgkMaxVertexCoincidence = 2.; // upper border

	// are there SPD tracklets which are not assigned to tracks? (0 = no, 1 = yes)
	static const Int_t fgkNTrackletResiduals      = 2; // number of bins
	static const Double_t fgkMinTrackletResiduals = 0.; // lower border
	static const Double_t fgkMaxTrackletResiduals = 2.; // upper border

	//-- MC event information
	static const Int_t fgkNProcessType      = 4; // kBinDD = number of bins
	static const Double_t fgkMinProcessType = 0.; // kBinND = lower border
	static const Double_t fgkMaxProcessType = 4.; // kBinDD = upper border


	//== EMPTY EVENT STUDY (ThnEmptyEvents) ======================================
	// event type
	static const Int_t fgkNEventType      = 5; // kBinEventE = number of bins (5)
	static const Double_t fgkMinEventType = 1.; // kBinEventI = lower border (1)
	static const Double_t fgkMaxEventType = 6.; // kBinEventE+kBinEventI = u. b.

	// multiplicities (reused for different detectors and ways of counting)
	static const Int_t fgkNMult      = 32; // number of bins
	static const Double_t fgkMinMult = 0.; // lower border
	static const Double_t fgkMaxMult = 31.; // upper border

	// multplicities - extended range
	// (reused for different detectors and ways of counting)
	static const Int_t fgkNMultW      = 64; // number of bins
	static const Double_t fgkMinMultW = 0; // lower border
	static const Double_t fgkMaxMultW = 63; // upper border

	//== MULTIPLICITY STUDY (TnnMultiplicity) ====================================
	// number of ITSTPC tracks in event
	static const Int_t fgkNNch      = 51; // number of bins
	static const Double_t fgkMinNch = 0.; // lower border
	static const Double_t fgkMaxNch = 51.; // upper border

	// number of ITS standalone tracks in event
	static const Int_t fgkNNsoft      = 11; // number of bins
	static const Double_t fgkMinNsoft = 0.; // lower border
	static const Double_t fgkMaxNsoft = 11.; // upper border

	// combined multiplicity
	static const Int_t fgkNNcomb      = 61; // number of bins
	static const Double_t fgkMinNcomb = 0.; // lower border
	static const Double_t fgkMaxNcomb = 61.; // upper border

	// gap configuration is reused from THnMother

	// number of residual tracks
	static const Int_t fgkNNresidualTracks      = 11; // number of bins
	static const Double_t fgkMinNresidualTracks = 0.; // lower border
	static const Double_t fgkMaxNresidualTracks = 11.; // upper border

	// number of residual tracklets
	static const Int_t fgkNNresidualTracklets      = 21; // number of bins
	static const Double_t fgkMinNresidualTracklets = 0.; // lower border
	static const Double_t fgkMaxNresidualTracklets = 21.; // upper border

	// vertex z-position
	static const Int_t fgkNVertexZ      = 20; // number of bins
	static const Double_t fgkMinVertexZ = -10.; // lower border
	static const Double_t fgkMaxVertexZ = 10.; // upper border

	// SPD and track vertex distance
	static const Int_t fgkNVerticesDistance      = 10; // number of bins
	static const Double_t fgkMinVerticesDistance = 0.; // lower border
	static const Double_t fgkMaxVerticesDistance = 5.; // upper border

	// MC process type is resued from THnMother
};

#endif
