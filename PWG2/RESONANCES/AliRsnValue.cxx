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

#include <Riostream.h>
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliCentrality.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"
#include "AliRsnDaughterDef.h"

#include "AliRsnValue.h"

ClassImp(AliRsnValue)

//_____________________________________________________________________________
AliRsnValue::AliRsnValue() :
   AliRsnTarget(),
   fComputedValue(0),
   fValueType(kValueTypes),
   fBinArray(0),
   fSupportObject(0x0)
{
//
// Default constructor without arguments.
// Initialize data members to meaningless values.
// This method is provided for ROOT streaming,
// but should never be used directly by a user.
//
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Int_t nbins, Double_t min, Double_t max) :
   AliRsnTarget(name, TargetType(type)),
   fComputedValue(0.0),
   fValueType(type),
   fBinArray(0),
   fSupportObject(0x0)
{
//
// Main constructor (version 1).
// This constructor defines in meaningful way all data members,
// and defined a fixed binnings, subdividing the specified interval
// into that many bins as specified in the integer argument.
// ---
// This method is also the entry point for all instances
// of this class which don't need to do binning (e.g.: TNtuple inputs),
// since arguments 3 to 5 have default values which don't create any
// binning array, in order not to allocate memory when this is useless.
//

   SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Double_t min, Double_t max, Double_t step) :
   AliRsnTarget(name, TargetType(type)),
   fComputedValue(0.0),
   fValueType(type),
   fBinArray(0),
   fSupportObject(0x0)
{
//
// Main constructor (version 2).
// This constructor defines in meaningful way all data members
// and creates enough equal bins of the specified size to cover
// the required interval.
//

   SetBins(min, max, step);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, EValueType type, Int_t nbins, Double_t *array) :
   AliRsnTarget(name, TargetType(type)),
   fComputedValue(0.0),
   fValueType(type),
   fBinArray(0),
   fSupportObject(0x0)
{
//
// Main constructor (version 3).
// This constructor defines in meaningful way all data members
// and creates a set of variable bins delimited by the passed array.
//

   SetBins(nbins, array);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue(const AliRsnValue& copy) :
   AliRsnTarget(copy),
   fComputedValue(copy.fComputedValue),
   fValueType(copy.fValueType),
   fBinArray(copy.fBinArray),
   fSupportObject(copy.fSupportObject)
{
//
// Copy constructor.
// Duplicates the binning array and copies all settings.
//
}

//_____________________________________________________________________________
AliRsnValue& AliRsnValue::operator=(const AliRsnValue& copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

   AliRsnTarget::operator=(copy);

   fComputedValue = copy.fComputedValue;
   fBinArray = copy.fBinArray;
   fSupportObject = copy.fSupportObject;
   fValueType = copy.fValueType;

   return (*this);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t min, Double_t max)
{
//
// Set binning for the axis in equally spaced bins
// where the number of bins, minimum and maximum are given.
//

   if (!nbins) {
      fBinArray.Set(0);
      return;
   }

   fBinArray.Set(nbins + 1);

   Double_t mymax = TMath::Max(min, max);
   Double_t mymin = TMath::Min(min, max);

   Int_t    k = 0;
   Double_t binSize = (mymax - mymin) / ((Double_t)nbins);

   fBinArray[0] = mymin;
   for (k = 1; k <= nbins; k++) fBinArray[k] = fBinArray[k - 1] + binSize;
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Double_t min, Double_t max, Double_t step)
{
//
// Set binning for the axis in equally spaced bins
// where the bin size, minimum and maximum are given.
//

   Double_t dblNbins = TMath::Abs(max - min) / step;
   Int_t    intNbins = ((Int_t)dblNbins) + 1;

   SetBins(intNbins, min, max);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t *array)
{
//
// Set binning for the axis in unequally spaced bins
// using the same way it is done in TAxis
//

   if (!nbins) {
      fBinArray.Set(0);
      return;
   }

   fBinArray.Adopt(nbins, array);
}

//_____________________________________________________________________________
const char* AliRsnValue::GetValueTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (fValueType) {
      case kTrackP:               return "SingleTrackPtot";
      case kTrackPt:              return "SingleTrackPt";
      case kTrackEta:             return "SingleTrackEta";
      case kTrackY:               return "SingleTrackRapidity";
      case kTrackITSsignal:       return "SingleTrackITSsignal";
      case kTrackTPCsignal:       return "SingleTrackTPCsignal";
      case kTrackTOFsignal:       return "SingleTrackTOFsignal";
      case kTrackLength:          return "SingleTrackLength";
      case kPairP1:               return "PairPtotDaughter1";
      case kPairP2:               return "PairPtotDaughter2";
      case kPairP1t:              return "PairPtDaughter1";
      case kPairP2t:              return "PairPtDaughter2";
      case kPairP1z:              return "PairPzDaughter1";
      case kPairP2z:              return "PairPzDaughter2";
      case kPairInvMass:          return "PairInvMass";
      case kPairInvMassMC:        return "PairInvMassMC";
      case kPairInvMassRes:       return "PairInvMassResolution";
      case kPairPt:               return "PairPt";
      case kPairPz:               return "PairPz";
      case kPairEta:              return "PairEta";
      case kPairMt:               return "PairMt";
      case kPairY:                return "PairY";
      case kPairPhi:              return "PairPhi";
      case kPairPhiMC:            return "PairPhiMC";
      case kPairPtRatio:          return "PairPtRatio";
      case kPairDipAngle:         return "PairDipAngle";
      case kPairCosThetaStar:     return "PairCosThetaStar";
      case kPairQInv:             return "PairQInv";
      case kPairAngleToLeading:   return "PairAngleToLeading";
      case kEventLeadingPt:       return "EventLeadingPt";
      case kEventMult:            return "EventMult";
      case kEventMultMC:          return "EventMultMC";
      case kEventMultESDCuts:     return "EventMultESDCuts";
      case kEventVz:              return "EventVz";
      case kEventCentralityV0:    return "EventCentralityV0";
      case kEventCentralityTrack: return "EventCentralityTrack";
      case kEventCentralityCL1:   return "EventCentralityCL1";
      default:                    return "Undefined";
   }
}

//_____________________________________________________________________________
Bool_t AliRsnValue::Eval(TObject *object, Bool_t useMC)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, computes the required value.
// The output of the function tells if computing was successful,
// and the values must be taken with GetValue().
//

   // coherence check 
   // (which also casts object to AliRsnTarget data members)
   if (!TargetOK(object)) return kFALSE;
   if (IsAllNull()) {
      AliError("TargetOK passed but all pointers are NULL");
      return kFALSE;
   }

   // cast the input to the allowed types
   AliESDEvent *esdev  = fgCurrentEvent->GetRefESD();
   AliESDtrack *esdt   = 0x0;
   AliAODTrack *aodt   = 0x0;
   AliAODPid   *pidObj = 0x0;
   
   // conditional initializations
   if (fDaughter) {
      esdt = fDaughter->GetRefESDtrack();
      aodt = fDaughter->GetRefAODtrack();
      if (aodt) pidObj = aodt->GetDetPid();
   }

   // common variables
   TLorentzVector pRec;   // 4-momentum for single track or pair sum (reco)
   TLorentzVector pSim;   // 4-momentum for single track or pair sum (MC)
   TLorentzVector pRec0;  // 4-momentum of first daughter (reco)
   TLorentzVector pSim0;  // 4-momentum of first daughter (MC)
   TLorentzVector pRec1;  // 4-momentum of second daughter (reco)
   TLorentzVector pSim1;  // 4-momentum of second daughter (MC)

   // initialize the above 4-vectors according to the 
   // expected target type (which has been checked above)
   switch (fTargetType) {
      case AliRsnTarget::kDaughter:
         pRec = fDaughter->Psim();
         pSim = fDaughter->Prec();
         break;
      case AliRsnTarget::kMother:
         pRec  = fMother->Sum();
         pSim  = fMother->SumMC();
         pRec0 = fMother->GetDaughter(0)->Prec();
         pRec1 = fMother->GetDaughter(1)->Prec();
         pSim0 = fMother->GetDaughter(0)->Psim();
         pSim1 = fMother->GetDaughter(1)->Psim();
         break;
      case AliRsnTarget::kEvent:
         if (!fgCurrentEvent) {
            AliError(Form("[%s] current event not initialized", GetName()));
            return kFALSE;
         }
         break;
      default:
         AliError(Form("[%s] Wrong type", GetName()));
         return kFALSE;
   }

   // cast the support object to the types which could be needed
   AliRsnPairDef     *pairDef     = 0x0;
   AliRsnDaughterDef *daughterDef = 0x0;
   AliESDpid         *esdPID      = 0x0;
   if (fSupportObject) {
      if (fSupportObject->InheritsFrom(AliRsnPairDef    ::Class())) pairDef     = static_cast<AliRsnPairDef*>(fSupportObject);
      if (fSupportObject->InheritsFrom(AliRsnDaughterDef::Class())) daughterDef = static_cast<AliRsnDaughterDef*>(fSupportObject);
      if (fSupportObject->InheritsFrom(AliESDpid        ::Class())) esdPID      = static_cast<AliESDpid*>(fSupportObject);
   }
   
   // compute value depending on types in the enumeration
   // if the type does not match any available choice, or if
   // the computation is not doable due to any problem
   // (not initialized support object, wrong values, risk of floating point errors)
   // the method returng kFALSE and sets the computed value to a very large number
   switch (fValueType) {
      case kTrackP:
         // single track:
         // total momentum 
         fComputedValue = useMC ? pSim.Mag() : pRec.Mag();
         return kTRUE;
      case kTrackPt:
         // single track:
         // transverse momentum
         fComputedValue = useMC ? pSim.Perp() : pRec.Perp();
         return kTRUE;
      case kTrackEta:
         // single track:
         // pseudo-rapidity
         fComputedValue = useMC ? pSim.Eta() : pRec.Eta();
         return kTRUE;
      case kTrackY:
         // single track:
         // rapidity (requires an AliRsnDaughterDef support)
         if (daughterDef) {
            pRec.SetXYZM(pRec.X(), pRec.Y(), pRec.Z(), daughterDef->GetMass());
            pSim.SetXYZM(pSim.X(), pSim.Y(), pSim.Z(), daughterDef->GetMass());
            fComputedValue = useMC ? pSim.Rapidity() : pRec.Rapidity();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Required a correctly initialized AliRsnDaughterDef support object to compute this value", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         } 
      case kTrackITSsignal:
         // single track:
         // ITS signal (successful only for tracks)
         if (esdt) {
            fComputedValue = esdt->GetITSsignal();
            return kTRUE;
         }
         else if (aodt && pidObj) {
            fComputedValue = pidObj->GetITSsignal();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Detector signals can be computed only on reconstructed tracks", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      case kTrackTPCsignal:
         // single track:
         // TPC signal (successful only for tracks)
         if (esdt) {
            fComputedValue = esdt->GetTPCsignal();
            return kTRUE;
         }
         else if (aodt && pidObj) {
            fComputedValue = pidObj->GetTPCsignal();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Detector signals can be computed only on reconstructed tracks", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      case kTrackTOFsignal:
         // single track:
         // TOF signal (successful only for tracks, for ESD requires an AliESDpid support)
         if (esdt) {
            if (!esdPID) {
               AliError(Form("[%s] Required a correctly initialized AliESDpid support object to compute this value with ESDs", GetName()));
               fComputedValue = fgkVeryBig;
               return kFALSE;
            }
            else {
               fComputedValue = (esdt->GetTOFsignal() - esdPID->GetTOFResponse().GetStartTime(esdt->P()));
               return kTRUE;
            }
         }
         else if (aodt && pidObj) {
            fComputedValue = pidObj->GetTOFsignal();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Detector signals can be computed only on reconstructed tracks", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      case kTrackLength:
         // single tracks:
         // integrated length (computed only on ESDs)
         if (esdt) {
            fComputedValue = esdt->GetIntegratedLength();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Length information not available in AODs", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      //---------------------------------------------------------------------------------------------------------------------
      case kPairP1:
         // pair:
         // momentum of first daughter (which matches definition #1 in pairDef)
         fComputedValue = useMC ? pSim0.Mag() : pRec0.Mag();
         return kTRUE;
      case kPairP2:
         // pair:
         // momentum of second daughter (which matches definition #2 in pairDef)
         fComputedValue = useMC ? pSim1.Mag() : pRec1.Mag();
         return kTRUE;
      case kPairP1t:
         // pair:
         // transverse momentum of first daughter 
         fComputedValue = useMC ? pSim0.Perp() : pRec0.Perp();
         return kTRUE;
      case kPairP2t:
         // pair:
         // transverse momentum of second daughter 
         fComputedValue = useMC ? pSim1.Perp() : pRec1.Perp();
         return kTRUE;
      case kPairP1z:
         // pair:
         // longitudinal momentum of first daughter 
         fComputedValue = useMC ? pSim0.Z() : pRec0.Z();
         return kTRUE;
      case kPairP2z:
         // pair:
         // longitudinal momentum of second daughter 
         fComputedValue = useMC ? pSim1.Z() : pRec1.Z();
         return kTRUE;
      case kPairInvMass:
         // pair:
         // invariant mass
         fComputedValue = useMC ? pSim.M() : pRec.M();
         return kTRUE;
      case kPairInvMassRes:
         // pair:
         // invariant mass resolution (requires MC)
         if (TMath::Abs(pSim.M()) > 0.0) {
            fComputedValue = (pSim.M() - pRec.M()) / pSim.M();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Caught a null MC mass", GetName()));
            return kFALSE;
         }
      case kPairPt:
         // pair:
         // total transverse momentum
         fComputedValue = useMC ? pSim.Perp() : pRec.Perp();
         return kTRUE;
      case kPairEta:
         // pair:
         // pseudo-rapidiry
         fComputedValue = useMC ? pSim.Eta() : pRec.Eta();
         return kTRUE;
      case kPairMt:
         // pair:
         // transverse mass (requires an AliRsnPairDef to get mass hypothesis)
         if (pairDef) {
            pRec.SetXYZM(pRec.X(), pRec.Y(), pRec.Z(), pairDef->GetMotherMass());
            pSim.SetXYZM(pSim.X(), pSim.Y(), pSim.Z(), pairDef->GetMotherMass());
            fComputedValue = useMC ? pSim.Mt() : pRec.Mt();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Required a correctly initialized AliRsnPairDef support object to compute this value", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      case kPairY:
         // pair:
         // rapidity (requires an AliRsnPairDef to get mass hypothesis)
         if (pairDef) {
            pRec.SetXYZM(pRec.X(), pRec.Y(), pRec.Z(), pairDef->GetMotherMass());
            pSim.SetXYZM(pSim.X(), pSim.Y(), pSim.Z(), pairDef->GetMotherMass());
            fComputedValue = useMC ? pSim.Rapidity() : pRec.Rapidity();
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Required a correctly initialized AliRsnPairDef support object to compute this value", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      case kPairPhi:
         // pair:
         // phi angle of total momentum
         fComputedValue = useMC ? pSim.Phi() : pRec.Phi();
         return kTRUE;
      case kPairPtRatio:
         // pair:
         // ratio of relative sum and difference of daughter transverse momenta
         if (useMC) {
            fComputedValue  = TMath::Abs(pSim0.Perp() - pSim1.Perp());
            fComputedValue /= TMath::Abs(pSim0.Perp() + pSim1.Perp());
         } else {
            fComputedValue  = TMath::Abs(pRec0.Perp() - pRec1.Perp());
            fComputedValue /= TMath::Abs(pRec0.Perp() + pRec1.Perp());
         }
         return kTRUE;
      case kPairDipAngle:
         // pair:
         // dip-angle in the transverse-Z plane
         // (used to check conversion electrons)
         if (useMC) {
            fComputedValue  = pSim0.Perp() * pSim1.Perp() + pSim0.Z() * pSim1.Z();
            fComputedValue /= pSim0.Mag() * pSim1.Mag();
         } else {
            fComputedValue  = pRec0.Perp() * pRec1.Perp() + pRec0.Z() * pRec1.Z();
            fComputedValue /= pRec0.Mag() * pRec1.Mag();
         }
         fComputedValue = TMath::Abs(TMath::ACos(fComputedValue));
         return kTRUE;
      case kPairCosThetaStar:
         // pair:
         // cosine of theta star angle
         // (= angle of first daughter to the total momentum, in resonance rest frame)
         fComputedValue = fMother->CosThetaStar(useMC);
         return kTRUE;
      case kPairQInv:
         // pair:
         // Q-invariant
         pSim0 -= pSim1;
         pRec0 -= pRec1;
         fComputedValue = useMC ? pSim0.M() : pRec0.M();
         return kTRUE;
      case kPairAngleToLeading:
         // pair:
         // angle w.r. to leading particle (if any)
         {
            int ID1 = (fMother->GetDaughter(0))->GetID();
            int ID2 = (fMother->GetDaughter(1))->GetID();
            Int_t leadingID = fgCurrentEvent->GetLeadingParticleID();
            if (leadingID == ID1 || leadingID == ID2) {
               fComputedValue = -99.;
               return kFALSE;
            }
            AliRsnDaughter leadingPart = fgCurrentEvent->GetDaughter(leadingID);
            AliVParticle  *ref = leadingPart.GetRef();
            fComputedValue = ref->Phi() - fMother->Sum().Phi();
            //return angle w.r.t. leading particle in the range -pi/2, 3/2pi
            while (fComputedValue >= 1.5 * TMath::Pi()) fComputedValue -= 2 * TMath::Pi();
            while (fComputedValue < -0.5 * TMath::Pi()) fComputedValue += 2 * TMath::Pi();
         }
         return kTRUE;
      //---------------------------------------------------------------------------------------------------------------------
      case kEventMult:
         // event:
         // multiplicity of tracks
         fComputedValue = (Double_t)fgCurrentEvent->GetMultiplicity(0x0);
         return kTRUE;
      case kEventMultMC:
         // event:
         // multiplicity of MC tracks
         fComputedValue = (Double_t)fgCurrentEvent->GetMultiplicityMC();
      case kEventMultESDCuts:
         // event:
         // multiplicity of good quality tracks 
         if (esdev) {
            fComputedValue = AliESDtrackCuts::GetReferenceMultiplicity(esdev, kTRUE);
            return kTRUE;
         }
         else {
            AliError(Form("[%s] Can be computed only on ESDs", GetName()));
            fComputedValue = fgkVeryBig;
            return kFALSE;
         }
      case kEventLeadingPt: 
         // event:
         // transverse momentum of leading particle
         {
            int leadingID = fgCurrentEvent->GetLeadingParticleID(); //fEvent->SelectLeadingParticle(0);
            if (leadingID >= 0) {
               AliRsnDaughter leadingPart = fgCurrentEvent->GetDaughter(leadingID);
               AliVParticle *ref = leadingPart.GetRef();
               fComputedValue = ref->Pt();
            } else fComputedValue = 0;
         }
         return kTRUE;
      case kEventVz:
         // event:
         // Z position of primary vertex
         fComputedValue = fgCurrentEvent->GetRef()->GetPrimaryVertex()->GetZ();
         return kTRUE;
      case kEventCentralityV0:
         // event:
         // centrality using V0 method
         if (esdev) {
            AliCentrality *centrality = esdev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("V0M");
            return kTRUE;
         } else {
            AliError(Form("[%s] Centrality computation is implemented for ESDs only up to now", GetName()));
            return kFALSE;
         }
      case kEventCentralityTrack:
         // event:
         // centrality using tracks method
         if (esdev) {
            AliCentrality *centrality = esdev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("TRK");
            return kTRUE;
         } else {
            AliError(Form("[%s] Centrality computation is implemented for ESDs only up to now", GetName()));
            return kFALSE;
         }
      case kEventCentralityCL1:
         // event:
         // centrality using CL1 method
         if (esdev) {
            AliCentrality *centrality = esdev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("CL1");
            return kTRUE;
         } else {
            AliError(Form("[%s] Centrality computation is implemented for ESDs only up to now", GetName()));
            return kFALSE;
         }
      default:
         AliError(Form("[%s] Invalid value type for this computation", GetName()));
         return kFALSE;
   }
}

//_____________________________________________________________________________
void AliRsnValue::Print(Option_t * /*option */) const
{
//
// Print informations about this object
//

   AliInfo("=== VALUE INFO =================================================");
   AliInfo(Form(" Name                  : %s", GetName()));
   AliInfo(Form(" Type                  : %s", GetValueTypeName()));
   AliInfo(Form(" Current computed value: %f", fComputedValue));
   Int_t i;
   for (i = 0; i < fBinArray.GetSize(); i++) {
      AliInfo(Form(" Bin limit #%d         = %f", i, fBinArray[i]));
   }
   AliInfo(Form(" Support object        : %s", (fSupportObject ? fSupportObject->ClassName() : " NO SUPPORT")));
   AliInfo("=== END VALUE INFO =============================================");
}

//_____________________________________________________________________________
RSNTARGET AliRsnValue::TargetType(EValueType type)
{
//
// This method assigns the target to be expected by this object
// in the computation, depending on its type chosen in the enum.
//

   if (type < kTrackValues)
      return AliRsnTarget::kDaughter;
   else if (type < kPairValues)
      return AliRsnTarget::kMother;
   else
      return AliRsnTarget::kEvent;
}
