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

//_________________________________________________________________________
//  Enumerated types for use in JetFinder classes
//
//*-- Author: Mark Horner (LBL/UCT)
//

#ifndef ALIEMCALJETFINDERALGOBGCALCTYPE_T
#define ALIEMCALJETFINDERALGOBGCALCTYPE_T

        typedef enum {	kRatio, kCone, kConstant  
	} AliEMCALJetFinderAlgoBGCalcType_t;
#endif

#ifndef ALIEMCALJETFINDERRESETTYPE_T
#define ALIEMCALJETFINDERRESETTYPE_T

        typedef enum {  kResetData, kResetTracks, kResetDigits, kResetParameters,
                        kResetAll, kResetPartons, kResetParticles, kResetJets
        } AliEMCALJetFinderResetType_t;
#endif

#ifndef  ALIEMCALJETFINDERTRACKTYPE_T
#define  ALIEMCALJETFINDERTRACKTYPE_T
	typedef enum {	kAllP, kEM, kCharged, kNeutral, kHadron, kChargedHadron, kNoTracks, kEMChargedPi0, kNoNeutronNeutrinoKlong
	} AliEMCALJetFinderTrackType_t;
#endif

#ifndef  ALIEMCALJETFINDERSMEARINGTYPE_T
#define  ALIEMCALJETFINDERSMEARINGTYPE_T
	typedef enum {	kSmear, kEfficiency , kSmearEffic, kPerfectTracks
	} AliEMCALJetFinderSmearingType_t;
#endif
	
#ifndef  ALIEMCALJETFINDEREMCALTYPE_T
#define  ALIEMCALJETFINDEREMCALTYPE_T
	typedef enum {	kHits, kTimeCut,kNoHits 
	} AliEMCALJetFinderEMCALType_t;
#endif

#ifndef  ALIEMCALJETFINDERFILETYPE_T
#define  ALIEMCALJETFINDERFILETYPE_T
	typedef enum {	kHijing,kPythia,kData
	} AliEMCALJetFinderFileType_t;
#endif

#ifndef  ALIEMCALJETFINDERUA1UNITFLAGTYPE_T
#define  ALIEMCALJETFINDERUA1UNITFLAGTYPE_T
 	typedef enum {	kInCurrentJet, kInJet, kOutJet, kBelowMinEt
	} AliEMCALJetFinderAlgoUA1UnitFlagType_t;
#endif

#ifndef  ALIEMCALJETFINDERUA1FILLUNITFLAGTYPE_T
#define  ALIEMCALJETFINDERUA1FILLUNITFLAGTYPE_T
 	typedef enum {	kFillTracksOnly, kFillDigitsOnly, kFillAll
	} AliEMCALJetFinderAlgoUA1FillUnitFlagType_t;
#endif
