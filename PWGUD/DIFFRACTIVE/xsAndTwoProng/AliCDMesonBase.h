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
	static const Int_t fgkNNcombined;      // number of bins
	static const Double_t fgkMinNcombined; // lower border
	static const Double_t fgkMaxNcombined; // upper border (#bins + lower border)

	// unlike sign or like sign charged?
	static const Int_t fgkNCombCh;      // kBinPPMM = number of bins
	static const Double_t fgkMinCombCh; // kBinPM = lower border
	static const Double_t fgkMaxCombCh; // kBinPPMM + kBinPM = upper border

	// two track PID, take care the enum changed from time to time
	static const Int_t fgkNCombPID;      // kBinPIDUnknown = number of bins
	static const Double_t fgkMinCombPID; // kBinPionE = lower border
	static const Double_t fgkMaxCombPID; // ...  = upper border

	// gap configuration, used for different detectors
	static const Int_t fgkNGapConfig;      // kBinNG = number of bins
	static const Double_t fgkMinGapConfig; // kBinDG = lower border
	static const Double_t fgkMaxGapConfig; // kBinNG + kBinDG = upper border

	//-- mother kinematics
	// invariant mass of the two-track system
	static const Int_t fgkNMass;      // number of bins
	static const Double_t fgkMinMass; // lower border
	static const Double_t fgkMaxMass; // upper border

	// transverse momentum of the two-track system
	static const Int_t fgkNMotherPt;      // number of bins
	static const Double_t fgkMinMotherPt; // lower border
	static const Double_t fgkMaxMotherPt; // upper border

	// cosine theta* (opening angle of the two daugther tracks in the
	// centre-of-mass system of the two-track/mother system)
	// **no meaning full variable**
	static const Int_t fgkNCTS;      // number of bins
	static const Double_t fgkMinCTS; // lower border
	static const Double_t fgkMaxCTS; // upper border

	// opening angle in the lab frame
	static const Int_t fgkNOA;      // number of bins
	static const Double_t fgkMinOA; // lower border
	static const Double_t fgkMaxOA; // upper border

	//-- daughter kinematics
	// transverse momentum of one of the two daughter particles
	// (randomly selected)
	static const Int_t fgkNDaughterPt;      // number of bins
	static const Double_t fgkMinDaughterPt; // lower border
	static const Double_t fgkMaxDaughterPt; // upper border

	// pseudo rapidity of one of the two daughter particles
	// (randomly selected)
	//static const Int_t fgkNDaughterEta;      // number of bins
	//static const Double_t fgkMinDaughterEta; // lower border
	//static const Double_t fgkMaxDaughterEta; // upper border

	//-- Event quality information
	// boolean values to reduce output size

	// are there tracks in addition to the ones selected using AliCDMesonTracks
	static const Int_t fgkNTrackResiduals;     // number of bins
	static const Double_t fgkMinTrackResiduals; // lower border
	static const Double_t fgkMaxTrackResiduals; // upper border

	// vertex with in +/-4cm (0 = no, 1 = yes)
	static const Int_t fgkNVertexZinRng;      // number of bins
	static const Double_t fgkMinVertexZinRng; // lower border
	static const Double_t fgkMaxVertexZinRng; // upper border

	// are the vertices from SPD and tracks within 0.5cm? (0 = no, 1 = yes)
	static const Int_t fgkNVertexCoincidence;      // number of bins
	static const Double_t fgkMinVertexCoincidence; // lower border
	static const Double_t fgkMaxVertexCoincidence; // upper border

	// are there SPD tracklets which are not assigned to tracks? (0 = no, 1 = yes)
	static const Int_t fgkNTrackletResiduals;      // number of bins
	static const Double_t fgkMinTrackletResiduals; // lower border
	static const Double_t fgkMaxTrackletResiduals; // upper border

	//-- MC event information
	static const Int_t fgkNProcessType;      // kBinDD = number of bins
	static const Double_t fgkMinProcessType; // kBinND = lower border
	static const Double_t fgkMaxProcessType; // kBinDD = upper border


	//== EMPTY EVENT STUDY (ThnEmptyEvents) ======================================
	// event type
	static const Int_t fgkNEventType;      // kBinEventE = number of bins
	static const Double_t fgkMinEventType; // kBinEventI = lower border
	static const Double_t fgkMaxEventType; // kBinEventE+kBinEventI = u. b.

	// multiplicities (reused for different detectors and ways of counting)
	static const Int_t fgkNMult;     // number of bins
	static const Double_t fgkMinMult; // lower border
	static const Double_t fgkMaxMult; // upper border

	// multplicities - extended range
	// (reused for different detectors and ways of counting)
	static const Int_t fgkNMultW;      // number of bins
	static const Double_t fgkMinMultW; // lower border
	static const Double_t fgkMaxMultW; // upper border

	//== MULTIPLICITY STUDY (TnnMultiplicity) ====================================
	// number of ITSTPC tracks in event
	static const Int_t fgkNNch;      // number of bins
	static const Double_t fgkMinNch; // lower border
	static const Double_t fgkMaxNch; // upper border

	// number of ITS standalone tracks in event
	static const Int_t fgkNNsoft;      // number of bins
	static const Double_t fgkMinNsoft; // lower border
	static const Double_t fgkMaxNsoft; // upper border

	// combined multiplicity
	static const Int_t fgkNNcomb;      // number of bins
	static const Double_t fgkMinNcomb; // lower border
	static const Double_t fgkMaxNcomb; // upper border

	// gap configuration is reused from THnMother

	// number of residual tracks
	static const Int_t fgkNNresidualTracks;      // number of bins
	static const Double_t fgkMinNresidualTracks; // lower border
	static const Double_t fgkMaxNresidualTracks; // upper border

	// number of residual tracklets
	static const Int_t fgkNNresidualTracklets;      // number of bins
	static const Double_t fgkMinNresidualTracklets; // lower border
	static const Double_t fgkMaxNresidualTracklets; // upper border

	// vertex z-position
	static const Int_t fgkNVertexZ;      // number of bins
	static const Double_t fgkMinVertexZ; // lower border
	static const Double_t fgkMaxVertexZ; // upper border

	// SPD and track vertex distance
	static const Int_t fgkNVerticesDistance;      // number of bins
	static const Double_t fgkMinVerticesDistance; // lower border
	static const Double_t fgkMaxVerticesDistance; // upper border

	// MC process type is resued from THnMother
};

#endif
