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
// AliCEPBase
// for
// AliAnalysisTaskCEP
//
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  rewritten by
//  Paul Buehler <paul.buehler@oeaw.ac.at>

#ifndef AliCEPBase_H
#define AliCEPBase_H

#include "TH1I.h"

class AliCEPBase : public TObject {

  public:

  // define some constants
  static const Int_t kdumval = -999;

	enum {

		// collission types
		kBinEventUnknown = 0,
		kBinEventI,             // Beam-Beam Interaction
		kBinEventA,             // Beam from A-side
		kBinEventC,             // Beam from C-side
		kBinEventAC,            // Beam from one side (there is no separation in 2011)
		kBinEventE,             // no beam (empty event)

		// track status bits
		kTTBaseLine        = 0,
    kTTTOFBunchCrossing= BIT( 0), // TOFBunchCrossing==0
    kTTTPCScluster     = BIT( 1), // number of TPC shared clusters <= fTPCnclsS(3)
    kTTDCA             = BIT( 2), // DCA to vertex is < 500
    kTTV0              = BIT( 3), // is a daughter of a V0
    kTTITSpure         = BIT( 4), // is an ITS pure track
    kTTZv              = BIT( 5), // |Zv-VtxZ| <= 6
    kTTeta             = BIT( 6), // -0.9<eta<0.9
    kTTAccITSTPC       = BIT( 7), // accepted by ITSTPC criteria
    kTTAccITSSA        = BIT( 8), // accepted by ITSSA  criteria
    kTTFiredChips      = BIT( 9), // passed FiredChips test
    kTTAccTPCOnly      = BIT(10), // passed standard TPCOnly criteria
    kTTSPDHit          = BIT(11), // has at least one SPD hit
    kTTCaloMatch       = BIT(12), // match with calorimeter cluster
    kTTAccV0daughter   = BIT(13), // passed standard V0daughter criteria

    // type of vertex
    kVtxUnknown         = 0,
    kVtxSPD             = BIT(0),  // from SPD tracklets
    kVtxTracks          = BIT(1),  // from tracks
    kVtxErrRes          = BIT(2),  // z-resolution of SPD vertex is out-of-bounds
    kVtxErrDif          = BIT(3),  // difference in z between SPD and track
                                    // vertex is out-of-bounds
    kVtxErrZ            = BIT(4),  // z-position of vertex is out-of-bounds
    kVtxAOD             = BIT(5),  // On AOD only primary vertex is stored

		// StatsFlow histogram entries
		// names for the bins are specified in AliCEPUtils.cxx
		kBinTotalInput = 0,
		kBinGoodInput,
		kBinMCEvent,
		kBinPhysEvent,
		kBinEventCut,
		kBinPhysel,
		kBinPileup,
		kBinClusterCut,
		kBinDGTrigger,
		kBinSharedCluster,
    kBinVtx,
		kBinMBOR,
		kBinMBAND,
		kBinnoV0,
		kBinnoFMD,
		kBinnoAD,
		kBinDG,
		kBinNDG,
		kBinSaved,
		kBinLastValue,

		// definition of bits in AliAnalysisTaskCEP::fCurrentEventCondition
		kETBaseLine       = BIT( 0),
    kETEventCut       = BIT( 1),
    kETPhyssel        = BIT( 2),
    kETPileup         = BIT( 3),
    kETClusterCut     = BIT( 4),
		kETDGTrigger      = BIT( 5),
    kETVtx            = BIT( 6),
		kETMBOR           = BIT( 7),
    kETMBAND          = BIT( 8),
		kETSPDA           = BIT( 9),
		kETSPDC           = BIT(10),
		kETTPCA           = BIT(11),
		kETTPCC           = BIT(12),
    kETV0A            = BIT(13),
		kETV0C            = BIT(14),
		kETFMDA           = BIT(15),
		kETFMDC           = BIT(16),
		kETADA            = BIT(17),
		kETADC            = BIT(18), 
		kETZDCA           = BIT(19), 
		kETZDCC           = BIT(20), 
		kETZDNA           = BIT(21),
		kETZDNC           = BIT(22),
    kETSClusterCut    = BIT(23),

    // MC process types
    kProctypeUnknown = 0,
    kProctypeMB,
    kProctypeND,
    kProctypeEL,
    kProctypeSDA,
    kProctypeSDB,
    kProctypeDD,
    kProctypeCD,

		// analysis task status bits
		// do not change the order in order to be backward compatible!
    
		kBitConfigurationSet      = BIT( 0), // if not set everything is active
    kBitisRun1                = BIT( 1), // is it run1
    kBitSaveAllEvents         = BIT( 2), // save all events
    kBitisMC                  = BIT( 3), // is Monte Carlo
    kBitQArnumStudy           = BIT( 4), // QA as function of rnum
    kBitSPDPileupStudy        = BIT( 5), // SPD pileup study
    kBitnClunTraStudy         = BIT( 6), // cluster vs tracklet study
    kBitVtxStudy              = BIT( 7), // Vtx study
    kBitTrackCutStudy         = BIT( 8), // track cut study
    kBitBBFlagStudy           = BIT( 9), // BBFlag study
    kBitV0Study               = BIT(10), // V0 study
    kBitFMDStudy              = BIT(11), // FMD study
    kBitEMCStudy              = BIT(12), // EMC study
		kBitRawBuffer             = BIT(13), // save a CEPRawEventBuffer
		kBitConfigurationVersion  = BIT(14)  // always set, last bit
	
  };

	ClassDef(AliCEPBase, 1);

};

#endif
