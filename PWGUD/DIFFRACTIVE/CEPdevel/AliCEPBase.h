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
    kTTTOFBunchCrossing= (1<< 0), // TOFBunchCrossing==0
    kTTTPCScluster     = (1<< 1), // number of TPC shared clusters <= fTPCnclsS(3)
    kTTDCA             = (1<< 2), // DCA to vertex is < 500
    kTTV0              = (1<< 3), // is a daughter of a V0
    kTTITSpure         = (1<< 4), // is an ITS pure track
    kTTZv              = (1<< 5), // |Zv-VtxZ| <= 6
    kTTeta             = (1<< 6), // -0.9<eta<0.9
    kTTAccITSTPC       = (1<< 7), // accepted by ITSTPC criteria
    kTTAccITSSA        = (1<< 8), // accepted by ITSSA  criteria
    kTTFiredChips      = (1<< 9), // passed FiredChips test

    // type of vertex
    kVtxUnknown         = 0,
    kVtxSPD             = (1<< 0),  // from ITS
    kVtxTracks          = (1<< 1),  // from tracks
    kVtxErrRes          = (1<< 2),  // z-resolution of SPD vertex is out-of-bounds
    kVtxErrDif          = (1<< 3),  // difference in z between SPD and track
                                    // vertex is out-of-bounds
    kVtxErrZ            = (1<< 4),  // z-position of vertex is out-of-bounds
    kVtxAOD             = (1<< 5),  // On AOD only primary vertex is stored

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
		kETBaseLine       = (1<< 0),
    kETEventCut       = (1<< 1),
    kETPhyssel        = (1<< 2),
    kETPileup         = (1<< 3),
    kETClusterCut     = (1<< 4),
		kETDGTrigger      = (1<< 5),
    kETVtx            = (1<< 6),
		kETMBOR           = (1<< 7),
    kETMBAND          = (1<< 8),
		kETSPDA           = (1<< 9),
		kETSPDC           = (1<<10),
		kETTPCA           = (1<<11),
		kETTPCC           = (1<<12),
    kETV0A            = (1<<13),
		kETV0C            = (1<<14),
		kETFMDA           = (1<<15),
		kETFMDC           = (1<<16),
		kETADA            = (1<<17),
		kETADC            = (1<<18), 
		kETZDCA           = (1<<19), 
		kETZDCC           = (1<<20), 
		kETZDNA           = (1<<21),
		kETZDNC           = (1<<22),
    kETSClusterCut    = (1<<23),

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
    
		kBitConfigurationSet      = (1<< 0), // if not set everything is active
    kBitisRun1                = (1<< 1), // is it run1
    kBitSaveAllEvents         = (1<< 2), // save all events
    kBitisMC                  = (1<< 3), // is Monte Carlo
    kBitQArnumStudy           = (1<< 4), // QA as function of rnum
    kBitSPDPileupStudy        = (1<< 5), // SPD pileup study
    kBitnClunTraStudy         = (1<< 6), // cluster vs tracklet study
    kBitVtxStudy              = (1<< 7), // Vtx study
    kBitTrackCutStudy         = (1<< 8), // track cut study
    kBitBBFlagStudy           = (1<< 9), // BBFlag study
    kBitV0Study               = (1<<10), // V0 study
    kBitFMDStudy              = (1<<11), // FMD study
		kBitConfigurationVersion  = (1<<12)  // always set, last bit
	
  };

	ClassDef(AliCEPBase, 1);

};

#endif
