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
//  This class contains all code which is used to compute any of the values
//  which can be of interest within a resonance analysis. Besides the obvious
//  invariant mass, it allows to compute other utility values on all possible
//  targets, in order to allow a wide spectrum of binning and checks.
//  When needed, this object can also define a binning in the variable which
//  it is required to compute, which is used for initializing axes of output
//  histograms (see AliRsnFunction).
//  The value computation requires this object to be passed the object whose
//  informations will be used. This object can be of any allowed input type
//  (track, Daughter, event), then this class must inherit from AliRsnTarget.
//  Then, when value computation is attempted, a check on target type is done
//  and computation is successful only if expected target matches that of the
//  passed object.
//  In some cases, the value computation can require a support external object,
//  which must then be passed to this class. It can be of any type inheriting
//  from TObject.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include "AliVVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliPIDResponse.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnDaughterDef.h"
#include "AliRsnDaughterDef.h"

#include "AliRsnValueDaughter.h"

ClassImp(AliRsnValueDaughter)

//_____________________________________________________________________________
AliRsnValueDaughter::AliRsnValueDaughter(const char *name, EType type) :
   AliRsnValue(name, AliRsnTarget::kDaughter),
   fType(type)
{
//
// Constructor
//
}

//_____________________________________________________________________________
AliRsnValueDaughter::AliRsnValueDaughter(const AliRsnValueDaughter& copy) :
   AliRsnValue(copy),
   fType(copy.fType)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnValueDaughter& AliRsnValueDaughter::operator=(const AliRsnValueDaughter& copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

   AliRsnValue::operator=(copy);
   fType = copy.fType;

   return (*this);
}

//_____________________________________________________________________________
const char* AliRsnValueDaughter::GetTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (fType) {
      case kP:           return "SingleTrackPtot";
      case kPt:          return "SingleTrackPt";
      case kPtpc:        return "SingleTrackPtpc";
      case kEta:         return "SingleTrackEta";
      case kITSsignal:   return "SingleTrackITSsignal";
      case kTPCsignal:   return "SingleTrackTPCsignal";
      case kTOFsignal:   return "SingleTrackTOFsignal";
      case kTPCnsigmaPi: return "SingleTrackTPCnsigmaPion";
      case kTPCnsigmaK:  return "SingleTrackTPCnsigmaKaon";
      case kTPCnsigmaP:  return "SingleTrackTPCnsigmaProton";
      case kTOFnsigmaPi: return "SingleTrackTOFnsigmaPion";
      case kTOFnsigmaK:  return "SingleTrackTOFnsigmaKaon";
      case kTOFnsigmaP:  return "SingleTrackTOFnsigmaProton";
      default:           return "Undefined";
   }
}

//_____________________________________________________________________________
Bool_t AliRsnValueDaughter::Eval(TObject *object)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, computes the required value.
// The output of the function tells if computing was successful,
// and the values must be taken with GetValue().
//
   
   // coherence check, which also casts object 
   // to AliRsnTarget data members and returns kFALSE
   // in case the object is NULL
   if (!TargetOK(object)) return kFALSE;

   // set a reference to the mother momentum
   AliVParticle   *ref   = fDaughter->GetRef();
   AliVParticle   *refMC = fDaughter->GetRefMC();
   AliVTrack      *track = fDaughter->Ref2Vtrack();
   if (fUseMCInfo && !refMC) {
      AliError("No MC");
      return kFALSE;
   }
   if (!fUseMCInfo && !ref) {
      AliError("No DATA");
      return kFALSE;
   }
   
   // compute value depending on types in the enumeration
   // if the type does not match any available choice, or if
   // the computation is not doable due to any problem
   // (not initialized support object, wrong values, risk of floating point errors)
   // the method returng kFALSE and sets the computed value to a meaningless number
   switch (fType) {
      case kP:
         fComputedValue = (fUseMCInfo ? refMC->P() : ref->P());
         return kTRUE;
      case kPt:
         fComputedValue = (fUseMCInfo ? refMC->Pt() : ref->Pt());
         return kTRUE;
      case kEta:
         fComputedValue = (fUseMCInfo ? refMC->Eta() : ref->Eta());
         return kTRUE;
      case kPtpc:
         if (track) {
            fComputedValue = track->GetTPCmomentum();
            return kTRUE;
         } else {
            AliWarning("Cannot get TPC momentum for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kITSsignal:
         if (track) {
            fComputedValue = track->GetITSsignal();
            return kTRUE;
         } else {
            AliWarning("Cannot get ITS signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTPCsignal:
         if (track) {
            fComputedValue = track->GetTPCsignal();
            return kTRUE;
         } else {
            AliWarning("Cannot get TPC signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTOFsignal:
         if (track) {
            fComputedValue = track->GetTOFsignal();
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTPCnsigmaPi:
         if (track) {
            AliPIDResponse *pid = fEvent->GetPIDResponse();
            fComputedValue = pid->NumberOfSigmasTPC(track, AliPID::kPion);
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTPCnsigmaK:
         if (track) {
            AliPIDResponse *pid = fEvent->GetPIDResponse();
            fComputedValue = pid->NumberOfSigmasTPC(track, AliPID::kKaon);
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTPCnsigmaP:
         if (track) {
            AliPIDResponse *pid = fEvent->GetPIDResponse();
            fComputedValue = pid->NumberOfSigmasTPC(track, AliPID::kProton);
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTOFnsigmaPi:
         if (track) {
            AliPIDResponse *pid = fEvent->GetPIDResponse();
            fComputedValue = pid->NumberOfSigmasTOF(track, AliPID::kPion);
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTOFnsigmaK:
         if (track) {
            AliPIDResponse *pid = fEvent->GetPIDResponse();
            fComputedValue = pid->NumberOfSigmasTOF(track, AliPID::kKaon);
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTOFnsigmaP:
         if (track) {
            AliPIDResponse *pid = fEvent->GetPIDResponse();
            fComputedValue = pid->NumberOfSigmasTOF(track, AliPID::kProton);
            return kTRUE;
         } else {
            AliWarning("Cannot get TOF signal for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      default:
         AliError(Form("[%s] Invalid value type for this computation", GetName()));
         return kFALSE;
   }
}
