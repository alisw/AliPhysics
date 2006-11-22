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


/* $Id$ */

//=========================================================================
// Modified class for JETAN
// Some flags have been changed
// Some are not used for the moment
// Author: magali.estienne@ires.in2p3.fr
//=========================================================================
//  Enumerated types for use in JetFinder classes
//
//*-- Author: Mark Horner (LBL/UCT)
//

#ifndef ALIJETFINDERBGCALCTYPE_T
#define ALIJETFINDERBGCALCTYPE_T

        typedef enum {	kRatio, kCone, kConstant  
	} AliJetFinderBGCalcType_t;
#endif

#ifndef ALIJETFINDERRESETTYPE_T
#define ALIJETFINDERRESETTYPE_T

        typedef enum {  kResetData, kResetTracks, kResetDigits, kResetParameters,
                        kResetAll, kResetPartons, kResetParticles, kResetJets
        } AliJetFinderResetType_t;
#endif

#ifndef  ALIJETFINDERTRACKTYPE_T
#define  ALIJETFINDERTRACKTYPE_T
	typedef enum {	kAllP, kEM, kCharged, kNeutral, kHadron, kChargedHadron, kNoTracks, kEMChargedPi0, kNoNeutronNeutrinoKlong
	} AliJetFinderTrackType_t;
#endif

#ifndef  ALIJETFINDERSMEARINGTYPE_T
#define  ALIJETFINDERSMEARINGTYPE_T
	typedef enum {	kSmear, kEfficiency , kSmearEffic, kPerfectTracks
	} AliJetFinderSmearingType_t;
#endif
	
#ifndef  ALIJETFINDEREMCALTYPE_T
#define  ALIJETFINDEREMCALTYPE_T
typedef enum {	kHits, kTimeCut,kNoHits, kDigits, kClusters 
	} AliJetFinderEmcalType_t;
#endif

#ifndef  ALIJETFINDERFILETYPE_T
#define  ALIJETFINDERFILETYPE_T
	typedef enum {	kHijing,kPythia,kData
	} AliJetFinderFileType_t;
#endif

#ifndef  ALIJETFINDERUA1UNITFLAGTYPE_T
#define  ALIJETFINDERUA1UNITFLAGTYPE_T
 	typedef enum {	kInCurrentJet, kInJet, kOutJet, kBelowMinEt
	} AliJetFinderUnitFlagType_t;
#endif

#ifndef  ALIJETFINDERUA1UNITCUTFLAGTYPE_T
#define  ALIJETFINDERUA1UNITCUTFLAGTYPE_T
 	typedef enum {	kPtSmaller, kPtHigher
	} AliJetFinderUnitCutFlagType_t;
#endif

#ifndef  ALIJETFINDERUA1UNITSIGNALFLAGTYPE_T
#define  ALIJETFINDERUA1UNITSIGNALFLAGTYPE_T
 	typedef enum {	kGood, kBad
	} AliJetFinderUnitSignalFlagType_t;
#endif

#ifndef  ALIJETFINDERUNITDETECTORFLAGTYPE_T
#define  ALIJETFINDERUNITDETECTORFLAGTYPE_T
 	typedef enum {	kTpc, kEmcal, kAll
	} AliJetFinderUnitDetectorFlagType_t;
#endif
