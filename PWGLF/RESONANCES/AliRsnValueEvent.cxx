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

#include "AliRsnValueEvent.h"

ClassImp(AliRsnValueEvent)

//_____________________________________________________________________________
AliRsnValueEvent::AliRsnValueEvent(const char *name, EType type) :
   AliRsnValue(name, AliRsnTarget::kEvent),
   fType(type)
{
//
// Constructor
//
}

//_____________________________________________________________________________
AliRsnValueEvent::AliRsnValueEvent(const AliRsnValueEvent &copy) :
   AliRsnValue(copy),
   fType(copy.fType)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnValueEvent &AliRsnValueEvent::operator=(const AliRsnValueEvent &copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

   AliRsnValue::operator=(copy);
   if (this == &copy)
      return *this;
   fType = copy.fType;

   return (*this);
}

//_____________________________________________________________________________
const char *AliRsnValueEvent::GetTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (fType) {
      case kLeadingPt:       return "EventLeadingPt";
      case kMult:            return "EventMult";
      case kMultMC:          return "EventMultMC";
      case kMultESDCuts:     return "EventMultESDCuts";
      case kMultSPD:         return "EventMultSPD";
      case kVz:              return "EventVz";
      case kCentralityV0:    return "EventCentralityV0";
      case kCentralityTrack: return "EventCentralityTrack";
      case kCentralityCL1:   return "EventCentralityCL1";
      default:               return "Undefined";
   }
}

//_____________________________________________________________________________
Bool_t AliRsnValueEvent::Eval(TObject *object)
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
   if (!fEvent->GetRef()) {
      AliWarning("NULL ref");
      return kFALSE;
   }

   // declare support variables
   AliCentrality *centrality = fEvent->GetRef()->GetCentrality();

   // compute value depending on types in the enumeration
   // if the type does not match any available choice, or if
   // the computation is not doable due to any problem
   // (not initialized support object, wrong values, risk of floating point errors)
   // the method returng kFALSE and sets the computed value to a meaningless number
   switch (fType) {
      case kMult:
         fComputedValue = (Double_t)fEvent->GetRef()->GetNumberOfTracks();
         return (fComputedValue >= 0);
      case kMultMC:
         fComputedValue = -999.0;
         if (fEvent->GetRefMC()) {
            if (fEvent->IsESD())
               fComputedValue = (Double_t)fEvent->GetRefMC()->GetNumberOfTracks();
            else {
               AliAODEvent *aod = (AliAODEvent *)fEvent->GetRefMC();
               TClonesArray *mcArray = (TClonesArray *)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
               if (mcArray) fComputedValue = (Double_t)mcArray->GetEntries();
            }
         }
         return (fComputedValue >= 0);
      case kMultESDCuts:
         fComputedValue = -999.0;
         if (fEvent->IsESD()) {
            fComputedValue = AliESDtrackCuts::GetReferenceMultiplicity(fEvent->GetRefESD(), kTRUE);
         } else {
            AliWarning("Cannot compute ESD cuts multiplicity in AOD");
            return kFALSE;
         }
         return (fComputedValue >= 0);
      case kMultSPD:
         fComputedValue = -999.0;
         if (fEvent->IsESD()) {
            const AliMultiplicity *mult = fEvent->GetRefESD()->GetMultiplicity();
            Float_t nClusters[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
            for(Int_t ilay = 0; ilay < 6; ilay++) nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
            fComputedValue = AliESDUtils::GetCorrSPD2(nClusters[1], fEvent->GetRef()->GetPrimaryVertex()->GetZ());
         } else {
            AliWarning("Cannot compute SPD multiplicity with AOD");
            return kFALSE;
         }
         return (fComputedValue >= 0);
      case kLeadingPt:
         if (fEvent->GetLeadingIndex() >= 0) {
            AliRsnDaughter leadingPart;
            fEvent->SetLeadingParticle(leadingPart);
            fComputedValue = leadingPart.GetRef()->Pt();
            return kTRUE;
         } else {
            AliError("Not found good leading particle");
            return kFALSE;
         }
      case kVz:
         fComputedValue = fEvent->GetRef()->GetPrimaryVertex()->GetZ();
         return kTRUE;
      case kCentralityV0:
         if (centrality) {
            fComputedValue = centrality->GetCentralityPercentile("V0M");
            return kTRUE;
         } else {
            AliError("Centrality undefined");
            return kFALSE;
         }
      case kCentralityTrack:
         if (centrality) {
            fComputedValue = centrality->GetCentralityPercentile("TRK");
            return kTRUE;
         } else {
            AliError("Centrality undefined");
            return kFALSE;
         }
      case kCentralityCL1:
         if (centrality) {
            fComputedValue = centrality->GetCentralityPercentile("CL1");
            return kTRUE;
         } else {
            AliError("Centrality undefined");
            return kFALSE;
         }
      default:
         AliError(Form("[%s] Invalid value type for this computation", GetName()));
         return kFALSE;
   }
}

//___________________________________________________________________
void AliRsnValueEvent::ApplyCentralityPatchAOD049(TObject *object)
{
  //
  //Apply centrality patch for AOD049 outliers
  //
  fComputedValue = -999.;
  if (!TargetOK(object)) return;
  if  (fEvent->IsESD()) {
    AliWarning(Form("Requested patch for AOD049 for ESD. Setting cent = %5.2f", fComputedValue));
    return;
  }

  if (fType!=AliRsnValueEvent::kCentralityV0){
    AliWarning(Form("Requested patch forAOD049 for wrong value (not centrality from V0). Setting cent = %5.2f", fComputedValue));
    return;
  }
  
  AliAODEvent * aodEvent = (AliAODEvent*)fEvent->GetRefAOD();
  if (!aodEvent) {
    AliWarning(Form("NULL ref to AOD event. Setting cent = %5.2f", fComputedValue));
    return;
  }
  
  // declare support variables
  AliCentrality *centrality = aodEvent->GetCentrality();
  if (!centrality) {
    AliWarning(Form("Cannot get centrality from AOD event. Setting cent = %5.2f", fComputedValue));
    return;
  }

  Float_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
  
  /*
  Bool_t isSelRun = kFALSE;
  Int_t selRun[5] = {138364, 138826, 138828, 138836, 138871};
    if(cent<0){
    Int_t quality = centrality->GetQuality();
    if(quality<=1){
      cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
    } else {
      Int_t runnum=aodEvent->GetRunNumber();
      for(Int_t ir=0;ir<5;ir++){
	if(runnum==selRun[ir]){
	  isSelRun=kTRUE;
	  break;
	}
      }
      if((quality==8||quality==9)&&isSelRun) cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
    }  
  } 
  */
  if(cent>=0.0) {
    Float_t v0 = 0.0;
    AliAODVZERO* aodV0 = (AliAODVZERO*) aodEvent->GetVZEROData();
    v0+=aodV0->GetMTotV0A();
    v0+=aodV0->GetMTotV0C();
    if ( (cent==0) && (v0<19500) ) {
      fComputedValue = -999.;//filtering issue
      AliDebug(3, Form("Filtering issue in centrality -> cent = %5.2f",fComputedValue));
      return;
    }
    Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
    Float_t val = 1.30552 +  0.147931 * v0;

    Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86, 
			     120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654, 
			     92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334, 
			     68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224, 
			     51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255, 
			     37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398, 
			     26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235, 
			     19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504, 
			     12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544, 
			     13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544 };

    if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )  {
      fComputedValue = -999.;//outlier
      AliDebug(3, Form("Outlier event in centrality -> cent = %5.2f",fComputedValue));
      return;
    }
  } else {
    //force it to be -999. whatever the negative value was
    cent = -999.;
  }
  fComputedValue=cent;
  return;
}
