// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawEventHeader                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef ALI_DATE
#include "event.h"
#endif

#include "AliRawEventHeader.h"


ClassImp(AliRawEventHeader)


//______________________________________________________________________________
Bool_t AliRawEventHeader::DataIsSwapped() const
{
   // Returns true if event data is swapped.

#ifdef ALI_DATE
   if (TEST_SYSTEM_ATTRIBUTE(fTypeAttribute, ATTR_EVENT_SWAPPED))
      return kTRUE;
#endif
   return kFALSE;
}

//______________________________________________________________________________
void AliRawEventHeader::Swap()
{
   // Swap header data.

   if (IsSwapped()) {
      fSize         = net2host(fSize);
      fMagic        = net2host(fMagic);
      fHeadLen      = net2host(fHeadLen);
      fVersion      = net2host(fVersion);
      fType         = net2host(fType);
      fRunNb        = net2host(fRunNb);
      for (int i = 0; i < kIdWords; i++)
         fId[i] = net2host(fId[i]);
      for (int i = 0; i < kTriggerWords; i++)
         fTriggerPattern[i] = net2host(fTriggerPattern[i]);
      for (int i = 0; i < kDetectorWords; i++)
         fDetectorPattern[i] = net2host(fDetectorPattern[i]);
      for (int i = 0; i < kAttributeWords; i++)
         fTypeAttribute[i] = net2host(fTypeAttribute[i]);
      fLDCId        = net2host(fLDCId);
      fGDCId        = net2host(fGDCId);
   }
}

//______________________________________________________________________________
UInt_t AliRawEventHeader::GetEventInRun() const
{
   // Get event number in run. Correct for fixed target mode which is used
   // in the Data Challenge Setup.

#ifdef ALI_DATE
   if (!TEST_SYSTEM_ATTRIBUTE(fTypeAttribute, ATTR_ORBIT_BC)) {
      return EVENT_ID_GET_NB_IN_RUN(fId);
   }
#endif
   return 0;
}

//______________________________________________________________________________
const char *AliRawEventHeader::GetTypeName() const
{
   // Get event type as a string.

   switch (GetType()) {
      case kStartOfRun:
         return "START_OF_RUN";
         break;
      case kEndOfRun:
         return "END_OF_RUN";
         break;
      case kStartOfRunFiles:
         return "START_OF_RUN_FILES";
         break;
      case kEndOfRunFiles:
         return "END_OF_RUN_FILES";
         break;
      case kStartOfBurst:
         return "START_OF_BURST";
         break;
      case kEndOfBurst:
         return "END_OF_BURST";
         break;
      case kPhysicsEvent:
         return "PHYSICS_EVENT";
         break;
      case kCalibrationEvent:
         return "CALIBRATION_EVENT";
         break;
      case kFormatError:
         return "EVENT_FORMAT_ERROR";
         break;
      default:
         return "*** UNKNOWN EVENT TYPE ***";
         break;
   }
}
