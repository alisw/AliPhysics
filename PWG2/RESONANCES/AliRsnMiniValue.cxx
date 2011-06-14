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

#include "AliLog.h"

#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"

#include "AliRsnMiniValue.h"

ClassImp(AliRsnMiniValue)

//_____________________________________________________________________________
AliRsnMiniValue::AliRsnMiniValue(EType type, Bool_t useMC) :
   TNamed(ValueName(type, useMC), ""),
   fType(type),
   fUseMCInfo(useMC)
{
//
// Constructor
//
}

//_____________________________________________________________________________
AliRsnMiniValue::AliRsnMiniValue(const AliRsnMiniValue& copy) :
   TNamed(copy),
   fType(copy.fType),
   fUseMCInfo(copy.fUseMCInfo)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnMiniValue& AliRsnMiniValue::operator=(const AliRsnMiniValue& copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

   TNamed::operator=(copy);
   fType = copy.fType;
   fUseMCInfo = copy.fUseMCInfo;

   return (*this);
}

//_____________________________________________________________________________
const char* AliRsnMiniValue::TypeName(EType type)
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (type) {
      case kVz:           return "EventVz";
      case kMult:         return "EventMult";
      case kPlaneAngle:   return "EventPlane";
      case kLeadingPt:    return "EventLeadingPt";
      case kPt:           return "Pt";
      case kPz:           return "Pz";
      case kInvMass:      return "InvMass";
      case kInvMassRes:   return "InvMassResolution";
      case kEta:          return "Eta";
      case kMt:           return "Mt";
      case kY:            return "Y";
      case kPtRatio:      return "PtRatio";
      case kDipAngle:     return "DipAngle";
      case kCosThetaStar: return "CosThetaStar";
      case kAngleLeading: return "AngleToLeading";
      default:            return "Undefined";
   }
}

//_____________________________________________________________________________
Float_t AliRsnMiniValue::Eval(AliRsnMiniPair *pair, AliRsnMiniEvent *event)
{
//
// Evaluation of the required value.
// In this implementation, fills the member 4-vectors with data
// coming from the object passed as argument, and then returns the value
//

   if (!pair && fType > kEventCuts) {
      AliError("Null pair passed!");
      return 1E20;
   }
   
   // compute value depending on types in the enumeration
   // if the type does not match any available choice, or if
   // the computation is not doable due to any problem
   // (not initialized support object, wrong values, risk of floating point errors)
   // the method returng kFALSE and sets the computed value to a meaningless number
   switch (fType) {
      // ---- event values -------------------------------------------------------------------------
      case kVz:
         return event->Vz();
      case kMult:
         return event->Mult();
      case kPlaneAngle:
         return event->Angle();
      case kLeadingPt:
         return 0.0;
      // ---- pair values --------------------------------------------------------------------------
      case kPt: 
         return pair->Pt(fUseMCInfo);
      case kInvMass:
         return pair->InvMass(fUseMCInfo);
      case kEta:
         return pair->Eta(fUseMCInfo);
      case kInvMassRes:
         return pair->InvMassRes();
      case kMt:
         return pair->Mt(fUseMCInfo);
      case kY:
         return pair->Y(fUseMCInfo);
      case kPtRatio:
         return pair->PtRatio(fUseMCInfo);
      case kDipAngle:
         return pair->DipAngle(fUseMCInfo);
      case kCosThetaStar:
         return pair->CosThetaStar(fUseMCInfo);
      case kAngleLeading:
         AliWarning("This method is not yet implemented");
         return 1E20;
      default:
         AliError("Invalid value type");
         return 1E20;
   }
}
