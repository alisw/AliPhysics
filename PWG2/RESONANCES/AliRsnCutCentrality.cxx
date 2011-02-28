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

////////////////////////////////////////////////////////////////////////////////
//
//  This cut is used to select PbPb events w.r. to their centrality class.
//  It uses the AliCentrality object, and specifically its method to get
//  the centrality percentile.
//  The centrality estimator is decided by the user, in a second string
//  in the constructor after the name.
//
////////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliRsnCutCentrality.h"

ClassImp(AliRsnCutCentrality)

//__________________________________________________________________________________________________
AliRsnCutCentrality::AliRsnCutCentrality(const char *name, const char *est, Double_t min, Double_t max) :
   AliRsnCut(name, AliRsnTarget::kEvent, min, max)
{
//
// Constructor
//

   SetTitle(est);
}

//__________________________________________________________________________________________________
AliRsnCutCentrality::AliRsnCutCentrality(const AliRsnCutCentrality& copy) :
   AliRsnCut(copy)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnCutCentrality& AliRsnCutCentrality::operator=(const AliRsnCutCentrality& copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   
   return (*this);
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutCentrality::IsSelected(TObject *object)
{
//
// Cut checking method.
// Checks current event and compares the percentile centrality
// with the allowed range.
//

   if (!TargetOK(object)) return kFALSE;
   
   AliESDEvent *esd = fEvent->GetRefESD();
   AliAODEvent *aod = fEvent->GetRefAOD();
   
   // esd
   if (esd) {
      AliDebug(AliLog::kDebug + 2, "Centrality for ESDs");
      AliCentrality *centrality = esd->GetCentrality();
      if (centrality) {
         fCutValueD = centrality->GetCentralityPercentile(fTitle.Data());
         return OkRangeD();
      } else {
         AliError("Centrality object is not present");
         return kFALSE;
      }
   }
   else {
      AliError("Currently the implementation works only with ESDs");
      return kFALSE;
   }
}
