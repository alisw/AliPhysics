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
//  developers: F. Bellini (fbellini@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "AliLog.h"

#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"

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
AliRsnMiniValue::AliRsnMiniValue(const AliRsnMiniValue &copy) :
   TNamed(copy),
   fType(copy.fType),
   fUseMCInfo(copy.fUseMCInfo)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnMiniValue &AliRsnMiniValue::operator=(const AliRsnMiniValue &copy)
{
//
// Assignment operator.
// Works like copy constructor.
//
   TNamed::operator=(copy);
   if (this == &copy)
      return *this;
   fType = copy.fType;
   fUseMCInfo = copy.fUseMCInfo;

   return (*this);
}

//_____________________________________________________________________________
const char *AliRsnMiniValue::TypeName(EType type)
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (type) {
      case kVz:           return "EventVz";
      case kSpherocity:   return "EventSpherocity";
      case kMult:         return "EventMult";
      case kRefMult:      return "EventReferenceMult";
      case kTracklets:    return "EventTracklets";
      case kPlaneAngle:   return "EventPlane";
      case kLeadingPt:    return "EventLeadingPt";
      case kPt:           return "Pt";
      case kPz:           return "Pz";
      case kInvMass:      return "InvMass";
      case kInvMassMother: return "InvMassMother";
      case kInvMassRes:   return "InvMassResolution";
      case kInvMassDiff:  return "InvMassDifference";
      case kEta:          return "Eta";
      case kMt:           return "Mt";
      case kY:            return "Y";
      case kPtRatio:      return "PtRatio";
      case kDipAngle:     return "DipAngle";
      case kCosThetaStar: return "CosThetaStar";
      case kCosThetaStarAbs:    return "CosThetaStarAbs";
      case kCosThetaJackson:    return "CosThetaJackson";
      case kCosThetaTransversity:    return "CosThetaTransversity";
      case kCosThetaToEventPlane:    return "CosThetaToEventPlane";
      case kAngleLeading: return "AngleToLeading";
      case kFirstDaughterPt: return "FirstDaughterPt";
      case kSecondDaughterPt: return "SecondDaughterPt";
      case kFirstDaughterP: return "FirstDaughterP";
      case kSecondDaughterP: return "SecondDaughterP";
      case kDCAproduct:   return "DaughterDCAproduct";
      case kFirstDaughterDCA: return "FirstDaughterDCA";
      case kSecondDaughterDCA: return "SecondDaughterDCA";
      case kNSisters:     return "NumberOfSisters";
      case kPairPtRes:        return "PairPtResolution";
      case kPairYRes:         return "PairYResolution";
      case kPhiV:         return "PhiV";
      case kAsym:         return "PairAsymmetry";
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
   Double_t p3[3]= {0.,0.,0.};
   AliRsnMiniParticle *l;
   TLorentzVector v;
   switch (fType) {
         // ---- event values -------------------------------------------------------------------------
      case kVz:
         return event->Vz();
      case kSpherocity:
         return event->Spherocity();
      case kMult:
         return event->Mult();
      case kRefMult:
         return event->RefMult();
      case kTracklets:
         return event->Tracklets();	 
      case kPlaneAngle:
         return event->Angle();
      case kLeadingPt:
         l = event->LeadingParticle(fUseMCInfo);
         if (l) {
            l->Set4Vector(v,-1.0,fUseMCInfo);
            return v.Pt();
         }
         return 0.0;
      case kPt:
         return pair->Pt(fUseMCInfo);
      case kInvMass:
         return pair->InvMass(fUseMCInfo);
      case kInvMassMother:
         return pair->InvMass(kTRUE);
      case kEta:
         return pair->Eta(fUseMCInfo);
      case kInvMassRes:
         return pair->InvMassRes();
      case kInvMassDiff:
         return pair->InvMassDiff();
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
      case kCosThetaStarAbs:
          return pair->CosThetaStarAbs(fUseMCInfo);           
      case kCosThetaJackson:
         return pair->CosThetaJackson(fUseMCInfo);
      case kCosThetaTransversity:
         return pair->CosThetaTransversity(fUseMCInfo);
      case kCosThetaToEventPlane:
         return pair->CosThetaToEventPlane(event, fUseMCInfo);
      case kAngleLeading:
         l = event->LeadingParticle(fUseMCInfo);
         if (l) {
            l->Set4Vector(v,-1.0,fUseMCInfo);
            Double_t angle = v.Phi() - pair->Sum(fUseMCInfo).Phi();

            //return angle w.r.t. leading particle in the range -pi/2, 3/2pi
            while (angle >= 1.5 * TMath::Pi()) angle -= 2 * TMath::Pi();
            while (angle < -0.5 * TMath::Pi()) angle += 2 * TMath::Pi();
            return angle;
         }
//         AliWarning("This method is not yet implemented");
         return 1E20;
      case kFirstDaughterPt:
         return pair->DaughterPt(0,fUseMCInfo);
      case kSecondDaughterPt:
         return pair->DaughterPt(1,fUseMCInfo);
      case kFirstDaughterP:
         pair->DaughterPxPyPz(0,fUseMCInfo, p3);
         return TMath::Sqrt(p3[0]*p3[0]+p3[1]*p3[1]+p3[2]*p3[2]);
      case kSecondDaughterP:
         pair->DaughterPxPyPz(1,fUseMCInfo, p3);
         return TMath::Sqrt(p3[0]*p3[0]+p3[1]*p3[1]+p3[2]*p3[2]);
      case kDCAproduct:
         return pair->DCAProduct();
      case kFirstDaughterDCA:
         return pair->DaughterDCA(0);
      case kSecondDaughterDCA:
         return pair->DaughterDCA(1);
      case kNSisters:
         return pair->NSisters();
      case kPairPtRes:
         return pair->PairPtRes();
      case kPairYRes:
         return pair->PairYRes();     
      case kPhiV:
         return pair->PhiV(fUseMCInfo);
      case kAsym:
         return pair->PairAsymmetry(fUseMCInfo);
      default:
         AliError("Invalid value type");
         return 1E20;
   }
}
