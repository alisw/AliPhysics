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
//  (track, pair, event), then this class must inherit from AliRsnTarget.
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

#include "Riostream.h"
#include "AliVVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"
#include "AliRsnDaughterDef.h"

#include "AliRsnValuePair.h"

ClassImp(AliRsnValuePair)

//_____________________________________________________________________________
AliRsnValuePair::AliRsnValuePair(const char *name, EType type) :
   AliRsnValue(name, AliRsnTarget::kMother),
   fType(type)
{
//
// Constructor
//
}

//_____________________________________________________________________________
AliRsnValuePair::AliRsnValuePair(const AliRsnValuePair& copy) :
   AliRsnValue(copy),
   fType(copy.fType)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnValuePair& AliRsnValuePair::operator=(const AliRsnValuePair& copy)
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
const char* AliRsnValuePair::GetTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (fType) {
      case kPt:           return "PairPt";
      case kPz:           return "PairPz";
      case kInvMass:      return "PairInvMass";
      case kInvMassRes:   return "PairInvMassResolution";
      case kEta:          return "PairEta";
      case kMt:           return "PairMt";
      case kY:            return "PairY";
      case kPtRatio:      return "PairPtRatio";
      case kDipAngle:     return "PairDipAngle";
      case kCosThetaStar: return "PairCosThetaStar";
      default:            return "Undefined";
   }
}

//_____________________________________________________________________________
Bool_t AliRsnValuePair::Eval(TObject *object)
{
//
// Evaluation of the required value.
// In this implementation, fills the member 4-vectors with data
// coming from the object passed as argument, and then returns the value
//

   // coherence check, which also casts object 
   // to AliRsnTarget data members and returns kFALSE
   // in case the object is NULL
   if (!TargetOK(object)) return kFALSE;
   
   // set a reference to the mother momentum
   TLorentzVector &sum   = fMother->Sum(fUseMCInfo);
   TLorentzVector &ref   = fMother->Ref(fUseMCInfo);
   TLorentzVector &p1    = fMother->GetDaughter(0)->P(fUseMCInfo);
   TLorentzVector &p2    = fMother->GetDaughter(1)->P(fUseMCInfo);
   
   // compute value depending on types in the enumeration
   // if the type does not match any available choice, or if
   // the computation is not doable due to any problem
   // (not initialized support object, wrong values, risk of floating point errors)
   // the method returng kFALSE and sets the computed value to a meaningless number
   switch (fType) {
      case kPt:
         fComputedValue = sum.Perp();
         return kTRUE;
      case kInvMass:
         fComputedValue = sum.M();
         return kTRUE;
      case kEta:
         fComputedValue = sum.Eta();
         return kTRUE;
      case kInvMassRes:
         fComputedValue  = fMother->Sum(kFALSE).M() - fMother->Sum(kTRUE).M();
         fComputedValue /= fMother->Sum(kTRUE).M();
         return kTRUE;
      case kMt:
         fComputedValue = ref.Mt();
         return kTRUE;
      case kY:
         fComputedValue = ref.Rapidity();
         return kTRUE;
      case kPtRatio:
         fComputedValue  = TMath::Abs(p1.Perp() - p2.Perp());
         fComputedValue /= TMath::Abs(p1.Perp() + p2.Perp());
         return kTRUE;
      case kDipAngle:
         fComputedValue  = p1.Perp() * p2.Perp() + p1.Z() * p2.Z();
         fComputedValue /= p1.Mag() * p2.Mag();
         return kTRUE;
      case kCosThetaStar:
         fComputedValue = fMother->CosThetaStar();
         return kTRUE;
      default:
         AliError(Form("[%s] Invalid value type for this computation", GetName()));
         return kFALSE;
   }
}
