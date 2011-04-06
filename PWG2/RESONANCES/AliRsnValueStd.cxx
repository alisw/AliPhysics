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

#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliCentrality.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"
#include "AliRsnDaughterDef.h"

#include "AliRsnValueStd.h"

ClassImp(AliRsnValueStd)

//_____________________________________________________________________________
AliRsnValueStd::AliRsnValueStd() :
   AliRsnValue(),
   fValueType(kValueTypes),
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
AliRsnValueStd::AliRsnValueStd
(const char *name, EValueType type, Int_t nbins, Double_t min, Double_t max) :
   AliRsnValue(name),
   fValueType(type),
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

   SetTargetType(TargetType(type));
   SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnValueStd::AliRsnValueStd
(const char *name, EValueType type, Double_t min, Double_t max, Double_t step) :
   AliRsnValue(name),
   fValueType(type),
   fSupportObject(0x0)
{
//
// Main constructor (version 2).
// This constructor defines in meaningful way all data members
// and creates enough equal bins of the specified size to cover
// the required interval.
//

   SetTargetType(TargetType(type));
   SetBins(min, max, step);
}

//_____________________________________________________________________________
AliRsnValueStd::AliRsnValueStd
(const char *name, EValueType type, Int_t nbins, Double_t *array) :
   AliRsnValue(name),
   fValueType(type),
   fSupportObject(0x0)
{
//
// Main constructor (version 3).
// This constructor defines in meaningful way all data members
// and creates a set of variable bins delimited by the passed array.
//

   SetTargetType(TargetType(type));
   SetBins(nbins, array);
}

//_____________________________________________________________________________
AliRsnValueStd::AliRsnValueStd(const AliRsnValueStd& copy) :
   AliRsnValue(copy),
   fValueType(copy.fValueType),
   fSupportObject(copy.fSupportObject)
{
//
// Copy constructor.
// Duplicates the binning array and copies all settings.
//
}

//_____________________________________________________________________________
AliRsnValueStd& AliRsnValueStd::operator=(const AliRsnValueStd& copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

   AliRsnValue::operator=(copy);

   fSupportObject = copy.fSupportObject;
   fValueType = copy.fValueType;

   return (*this);
}

//_____________________________________________________________________________
const char* AliRsnValueStd::GetValueTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (fValueType) {
      case kTrackP:               return "SingleTrackPtot";
      case kTrackPt:              return "SingleTrackPt";
      case kTrackPtpc:            return "SingleTrackPtpc";
      case kTrackEta:             return "SingleTrackEta";
      case kTrackY:               return "SingleTrackRapidity";
      case kTrackITSsignal:       return "SingleTrackITSsignal";
      case kTrackTPCsignal:       return "SingleTrackTPCsignal";
      case kTrackTOFsignal:       return "SingleTrackTOFsignal";
      case kTrackTOFbeta:         return "SingleTrackTOFbeta";
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
      case kEventMultSPD:         return "EventMultSPD";
      case kEventVz:              return "EventVz";
      case kEventCentralityV0:    return "EventCentralityV0";
      case kEventCentralityTrack: return "EventCentralityTrack";
      case kEventCentralityCL1:   return "EventCentralityCL1";
      default:                    return "Undefined";
   }
}

//_____________________________________________________________________________
Bool_t AliRsnValueStd::Eval(TObject *object, Bool_t useMC)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, computes the required value.
// The output of the function tells if computing was successful,
// and the values must be taken with GetValue().
//

   // utility variables
   Bool_t         success;
   const Double_t fgkVeryBig = 1E20;
   Double_t       time;
   Int_t          leadingID = -1;
   ULong_t        status = 0x0;

   // coherence check, which also casts object 
   // to AliRsnTarget data members and returns kFALSE
   // in case the object is NULL
   if (!TargetOK(object)) return kFALSE;

   // these variables are initialized
   // from the target object, once it
   // is casted to one of the expected
   // types (daughter/mother/event)
   // -- not all are initialized always
   TLorentzVector pRec;          // 4-momentum for single track or pair sum (reco)
   TLorentzVector pSim;          // 4-momentum for single track or pair sum (MC)
   TLorentzVector pRec0;         // 4-momentum of first daughter (reco)
   TLorentzVector pSim0;         // 4-momentum of first daughter (MC)
   TLorentzVector pRec1;         // 4-momentum of second daughter (reco)
   TLorentzVector pSim1;         // 4-momentum of second daughter (MC)
   AliESDEvent   *esdev  = 0x0;  // reference ESD event
   AliAODEvent   *aodev  = 0x0;  // reference AOD event
   AliESDtrack   *esdt   = 0x0;  // reference ESD track
   AliAODTrack   *aodt   = 0x0;  // reference AOD track
   AliAODPid     *pidObj = 0x0;  // reference AOD PID object

   // initialize the above 4-vectors according to the 
   // expected target type (which has been checked above)
   // in particular, the 'fEvent' data member of base AliRsnTarget
   // will be *always* well initialized if the TargetOK() returns kTRUE
   switch (fTargetType) {
      case AliRsnTarget::kDaughter:
         pRec = fDaughter->Prec();
         pSim = fDaughter->Psim();
         esdt = fDaughter->GetRefESDtrack();
         aodt = fDaughter->GetRefAODtrack();
         if (aodt) pidObj = aodt->GetDetPid();
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
         break;
      default:
         AliError(Form("[%s] Wrong type", GetName()));
         return kFALSE;
   }
   leadingID = fEvent->GetLeadingParticleID();
   esdev = fEvent->GetRefESD();
   aodev = fEvent->GetRefAOD();
   
   // if leading index is negative, assume that leading particle was not searched
   // and then searches for it
   if (leadingID < 0) {
      fEvent->SelectLeadingParticle();
      leadingID = fEvent->GetLeadingParticleID();
   }

   if (esdt) status = esdt->GetStatus();
   if (aodt) status = aodt->GetStatus();
   
   // these objects are all types of supports
   // which could be needed for some values
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
   // the method returng kFALSE and sets the computed value to a meaningless number
   switch (fValueType) {
      case kTrackP:
         // single track:
         // total momentum 
         fComputedValue = useMC ? pSim.Vect().Mag() : pRec.Vect().Mag();
         return kTRUE;
      case kTrackPt:
         // single track:
         // transverse momentum
         fComputedValue = useMC ? pSim.Perp() : pRec.Perp();
         return kTRUE;
      case kTrackPtpc:
         // single track:
         // transverse momentum
         if (esdt) {
            if (esdt->GetInnerParam()) {
               fComputedValue = esdt->GetInnerParam()->P();
               return kTRUE;
            } else {
               AliError(Form("[%s] TPC inner param is not initialized", GetName()));
               return kFALSE;
            }
         }
         else if (aodt && pidObj) {
            fComputedValue = pidObj->GetTPCmomentum();
            return kTRUE;
         } else {
            AliError(Form("[%s] Cannot retrieve TPC momentum", GetName()));
            return kFALSE;
         }
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
         // works only if the status is OK
         if ((status & AliESDtrack::kITSin) == 0) {
            AliDebug(AliLog::kDebug + 2, "Rejecting non-ITS track");
            return kFALSE;
         }
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
         // works only if the status is OK
         if ((status & AliESDtrack::kTPCin) == 0) {
            AliDebug(AliLog::kDebug + 2, "Rejecting non-TPC track");
            return kFALSE;
         }
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
         // works only if the status is OK
         if ((status & AliESDtrack::kTOFout) == 0 || (status & AliESDtrack::kTIME) == 0) {
            AliDebug(AliLog::kDebug + 2, "Rejecting non-TOF track");
            return kFALSE;
         }
         if (esdt) {
            if (!esdPID || !esdev) {
               AliError(Form("[%s] Required a correctly initialized AliRsnEvent and AliESDpid support object to compute this value with ESDs", GetName()));
               fComputedValue = fgkVeryBig;
               return kFALSE;
            }
            else {
               esdPID->SetTOFResponse(esdev, AliESDpid::kTOF_T0);
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
      case kTrackTOFbeta:
         // single track:
         // TOF beta (successful only for tracks, for ESD requires an AliESDpid support)
         if (esdt) {
            if (!esdPID) {
               AliError(Form("[%s] Required a correctly initialized AliESDpid support object to compute this value with ESDs", GetName()));
               fComputedValue = fgkVeryBig;
               return kFALSE;
            }
            else if (!esdev) {
               AliError(Form("[%s] Required a correctly initialized AliESDEvent to compute this value with ESDs", GetName()));
               fComputedValue = fgkVeryBig;
               return kFALSE;
            }
            else {
               esdPID->SetTOFResponse(esdev, AliESDpid::kTOF_T0);
               fComputedValue = esdt->GetIntegratedLength();
               time = (esdt->GetTOFsignal() - esdPID->GetTOFResponse().GetStartTime(esdt->P()));
               if (time > 0.0) {
                  fComputedValue /= time;
                  fComputedValue /= 2.99792458E-2;
                  return kTRUE;
               } else {
                  fComputedValue = fgkVeryBig;
                  return kFALSE;
               }
            }
         }
         else {
            AliError(Form("[%s] Length information not available in AODs", GetName()));
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
         fComputedValue = fMother->AngleToLeading(success);
         return success;
      //---------------------------------------------------------------------------------------------------------------------
      case kEventMult:
         // event:
         // multiplicity of tracks
         fComputedValue = (Double_t)fEvent->GetMultiplicityFromTracks();
         return (fComputedValue >= 0);
      case kEventMultMC:
         // event:
         // multiplicity of MC tracks
         fComputedValue = (Double_t)fEvent->GetMultiplicityFromMC();
         return (fComputedValue >= 0);
      case kEventMultESDCuts:
         // event:
         // multiplicity of good quality tracks 
         fComputedValue = fEvent->GetMultiplicityFromESDCuts();
         return (fComputedValue >= 0);
      case kEventMultSPD:
         // event:
         // multiplicity of good quality tracks 
         fComputedValue = fEvent->GetMultiplicityFromSPD();
         return (fComputedValue >= 0);
      case kEventLeadingPt: 
         // event:
         // transverse momentum of leading particle
         if (leadingID >= 0) {
            AliRsnDaughter leadingPart = fEvent->GetDaughter(leadingID);
            AliVParticle *ref = leadingPart.GetRef();
            fComputedValue = ref->Pt();
            return kTRUE;
         } else {
            AliError(Form("[%s] Leading ID has bad value (%d)", GetName(), leadingID));
            return kFALSE;
         }
      case kEventVz:
         // event:
         // Z position of primary vertex
         fComputedValue = fEvent->GetRef()->GetPrimaryVertex()->GetZ();
         return kTRUE;
      case kEventCentralityV0:
         // event:
         // centrality using V0 method
         if (esdev) {
            AliCentrality *centrality = esdev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("V0M");
            return kTRUE;
         } else if (aodev) {
            AliCentrality *centrality = aodev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("V0M");
            return kTRUE;
         } else {
            AliError(Form("[%s] Neither the ESD nor the AOD events are initialized", GetName()));
            return kFALSE;
         }
      case kEventCentralityTrack:
         // event:
         // centrality using tracks method
         if (esdev) {
            AliCentrality *centrality = esdev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("TRK");
            return kTRUE;
         } else if (aodev) {
            AliCentrality *centrality = aodev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("TRK");
            return kTRUE;
         } else {
            AliError(Form("[%s] Neither the ESD nor the AOD events are initialized", GetName()));
            return kFALSE;
         }
      case kEventCentralityCL1:
         // event:
         // centrality using CL1 method
         if (esdev) {
            AliCentrality *centrality = esdev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("CL1");
            return kTRUE;
         } else if (aodev) {
            AliCentrality *centrality = aodev->GetCentrality();
            fComputedValue = centrality->GetCentralityPercentile("CL1");
            return kTRUE;
         } else {
            AliError(Form("[%s] Neither the ESD nor the AOD events are initialized", GetName()));
            return kFALSE;
         }
      default:
         AliError(Form("[%s] Invalid value type for this computation", GetName()));
         return kFALSE;
   }
}

//_____________________________________________________________________________
void AliRsnValueStd::Print(Option_t *option) const
{
//
// Print informations about this object
//

   AliInfo("=== VALUE INFO =================================================");
   AliInfo(Form(" Name                  : %s", GetName()));
   AliInfo(Form(" Type                  : %s", GetValueTypeName()));
   AliInfo(Form(" Current computed value: %f", fComputedValue));
   if (!strcmp(option, "BINS")) {
      Int_t i;
      for (i = 0; i < fBinArray.GetSize(); i++) {
         AliInfo(Form(" Bin limit #%03d        = %f", i, fBinArray[i]));
      }
   }
   AliInfo(Form(" Support object        : %s", (fSupportObject ? fSupportObject->ClassName() : " NO SUPPORT")));
   AliInfo("=== END VALUE INFO =============================================");
}

//_____________________________________________________________________________
RSNTARGET AliRsnValueStd::TargetType(EValueType type)
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
