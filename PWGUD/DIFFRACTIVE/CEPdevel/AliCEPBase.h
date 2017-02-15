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
		kTTUnknown         = 0,
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
    kVtxUnknown = 0,
    kVtxSPD,                // from ITS
    kVtxTracks,             // from tracks

		// StatsFlow histogram entries
		// names for the bins are specified in AliCEPUtils.cxx
		kBinTotalInput = 0,
		kBinGoodInput,
		kBinMCEvent,
		kBinDGTrigger,
		kBinPileup,
		kBinPhysEvent,
		kBinEventCut,
		kBinPhySel,
		kBinMBOR,
		kBinMBAND,
		kBinClusterCut,
		kBinSharedCluster,
		kBinSaved,
		kBinLastValue,        // used to specify the correct histogram width

		// definition of bits in AliAnalysisTaskCEP::fCurrentGapCondition
		kBitBaseLine       = (1<< 0),
		kBitMBOR           = (1<< 1),
    kBitMBAND          = (1<< 2),
		kBitSPDA           = (1<< 3),
		kBitSPDC           = (1<< 4),
		kBitTPCA           = (1<< 5),
		kBitTPCC           = (1<< 6),
    kBitV0A            = (1<< 7),
		kBitV0C            = (1<< 8),
		kBitFMDA           = (1<< 9),
		kBitFMDC           = (1<<10),
		kBitADA            = (1<<11),
		kBitADC            = (1<<12),
		kBitZDCA           = (1<<13),
		kBitZDCC           = (1<<14),
		kBitZDNA           = (1<<15),
		kBitZDNC           = (1<<16),
		kBitDGTrigger      = (1<<17),
    
    // MC process types
    kProctypeUnknown = 0,
    kProctypeND,
    kProctypeSD,
    kProctypeDD,
    kProctypeCD,

		// analysis task status bits
		// do not change the order in order to be backward compatible!
    
		kBitConfigurationSet      = (1<< 0), // if not set everything is active
    kBitisRun1                = (1<< 1), // save all events
    kBitSaveAllEvents         = (1<< 2), // save all events
		kBitConfigurationVersion  = (1<< 3)  // always set, last bit
	
  };

	ClassDef(AliCEPBase, 1);

};

#endif
