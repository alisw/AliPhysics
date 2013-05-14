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
AliRsnValueDaughter::AliRsnValueDaughter(const AliRsnValueDaughter &copy) :
   AliRsnValue(copy),
   fType(copy.fType)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnValueDaughter &AliRsnValueDaughter::operator=(const AliRsnValueDaughter &copy)
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
const char *AliRsnValueDaughter::GetTypeName() const
{
//
// This method returns a string to give a name to each possible
// computation value.
//

   switch (fType) {
      case kP:            	       return "SingleTrackPtot";
      case kPt:           	       return "SingleTrackPt";
      case kPtpc:         	       return "SingleTrackPtpc";
      case kEta:          	       return "SingleTrackEta";
      case kMass:         	       return "SingleTrackMass";
      case kITSsignal:    	       return "SingleTrackITSsignal";
      case kTPCsignal:    	        return "SingleTrackTPCsignal";
      case kTOFsignal:    	        return "SingleTrackTOFsignal";
      case kTPCnsigmaPi:  	        return "SingleTrackTPCnsigmaPion";
      case kTPCnsigmaK:   	        return "SingleTrackTPCnsigmaKaon";
      case kTPCnsigmaP:   	        return "SingleTrackTPCnsigmaProton";
      case kTOFnsigmaPi:  	        return "SingleTrackTOFnsigmaPion";
      case kTOFnsigmaK:   	        return "SingleTrackTOFnsigmaKaon";
      case kTOFnsigmaP:   	        return "SingleTrackTOFnsigmaProton";
      case kCharge:                     return "SingleTrackCharge";
      case kPhi:                        return "SingleTrackPhi";
      case kPhiOuterTPC:                return "SingleTrackPhiOuterTPC";
      case kV0DCA:        	        return "V0DCAToPrimaryVertex";
      case kDaughterDCA:  	        return "V0DaughterDCA";
      case kCosPointAng:  	        return "V0CosineOfPointingAngle";
      case kLambdaProtonPIDCut:         return "V0LambdaProtonNsigma";
      case kAntiLambdaAntiProtonPIDCut: return "V0AntiLambdaAntiProtonNsigma";
      case kLambdaPionPIDCut:           return "V0LambdaPionNsigma";
      case kAntiLambdaAntiPionPIDCut:   return "V0AntiLambdaPionNsigma";
      default:                          return "Undefined";
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
   AliESDv0       *v0esd = fDaughter->Ref2ESDv0();
   AliAODv0       *v0aod = fDaughter->Ref2AODv0();
   AliESDEvent    *lESDEvent = fEvent->GetRefESD();
   
   Double_t xPrimaryVertex = -999.9;
   Double_t yPrimaryVertex = -999.9;
   Double_t zPrimaryVertex = -999.9;
   
   if(lESDEvent){
   
   xPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetX();
   yPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetY();
   zPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetZ();
   
   }
   
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
      case kMass:
         fComputedValue = (fUseMCInfo ? refMC->M() : ref->M());
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
      case kNITSclusters:
         if (track) {
            AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
            if (trackESD) {
               fComputedValue =  trackESD->GetITSclusters(0);
            } else {
               fComputedValue =  ((AliAODTrack *)track)->GetITSNcls();
            }
            return kTRUE;
         } else {
            AliWarning("Cannot get n ITS clusters for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kNTPCclusters:
         if (track) {
            AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
            if (trackESD) {
               fComputedValue =  trackESD->GetTPCclusters(0);
            } else {
               fComputedValue =  ((AliAODTrack *)track)->GetTPCNcls();
            }
            return kTRUE;
         } else {
            AliWarning("Cannot get n TPC clusters for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kITSchi2:
         if (track) {
            AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
            if (trackESD) {
               UShort_t nClustersITS = trackESD->GetITSclusters(0);
               fComputedValue =  trackESD->GetITSchi2()/Float_t(nClustersITS);
            } else {
               fComputedValue = ((AliAODTrack *)track)->Chi2perNDF();
            }
            return kTRUE;
         } else {
            AliWarning("Cannot get ITS chi^2 for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kTPCchi2:
         if (track) {
            AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
            if (trackESD) {
               UShort_t nClustersTPC = trackESD->GetTPCclusters(0);
               fComputedValue =  trackESD->GetTPCchi2()/Float_t(nClustersTPC);
            } else {
               fComputedValue = ((AliAODTrack *)track)->Chi2perNDF();
            }
            return kTRUE;
         } else {
            AliWarning("Cannot get TPC chi^2 for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kDCAXY:
         if (track) {
            AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
            if (trackESD) {
               Float_t b[2], bCov[3];
               trackESD->GetImpactParameters(b, bCov);
               fComputedValue =  b[0];
            } else {
               Double_t b[2]= {-999,-999}, cov[3];
               AliAODVertex *vertex = fEvent->GetRefAOD()->GetPrimaryVertex();
               if(vertex) {
                  track->PropagateToDCA(vertex, fEvent->GetRefAOD()->GetMagneticField(),kVeryBig, b, cov);
                  fComputedValue = b[0];
               } else {
                  fComputedValue = -999;
               }
            }
            return kTRUE;
         } else {
            AliWarning("Cannot get TPC chi^2 for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }
      case kDCAZ:
         if (track) {
            AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
            if (trackESD) {
               Float_t b[2], bCov[3];
               trackESD->GetImpactParameters(b, bCov);
               fComputedValue =  b[1];
            } else {
               Double_t b[2]= {-999,-999}, cov[3];
               AliAODVertex *vertex = fEvent->GetRefAOD()->GetPrimaryVertex();
               if(vertex) {
                  track->PropagateToDCA(vertex, fEvent->GetRefAOD()->GetMagneticField(),kVeryBig, b, cov);
                  fComputedValue = b[1];
               } else {
                  fComputedValue = -999;
               }

            }
            return kTRUE;
         } else {
            AliWarning("Cannot get TPC chi^2 for non-track object");
            fComputedValue = 0.0;
            return kFALSE;
         }

     case kCharge:
         fComputedValue = (fUseMCInfo ? refMC->Charge() : ref->Charge());
         return kTRUE;
   
      case kPhi:         
        fComputedValue = (fUseMCInfo ? (refMC->Phi()*TMath::RadToDeg()) : (ref->Phi()*TMath::RadToDeg()));      
        return kTRUE;

      case kPhiOuterTPC:    
        if (track) {
          Double_t pos[3]={0.,0.,0.};
          Double_t phiOut = -999.0;
          Double_t radius = 278.;//TPC outer (vessel) = 278 cm, TOF radius (active surf.) = 378 cm;  ref. PPR.1
          AliExternalTrackParam etp; //thanks to Andrea and Cristina
          etp.CopyFromVTrack(track);
          if(etp.GetXYZAt(radius, 5., pos)){
            phiOut=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
            if (phiOut<0) phiOut+= (2*TMath::Pi()*TMath::RadToDeg());
          }
          fComputedValue = phiOut;      
        } else {
          AliWarning("Cannot get phi at outer TPC radius for non-track object");
          fComputedValue = -99.0;
          return kFALSE;
        }
        return kTRUE;

      case kV0DCA:
         if(v0esd) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   fComputedValue = v0ESD->GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
	   return kTRUE;
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   fComputedValue = v0AOD->DcaV0ToPrimVertex();
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }	 	 
      case kDaughterDCA:
         if(v0esd) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   fComputedValue = v0ESD->GetDcaV0Daughters();
	   return kTRUE;
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   fComputedValue = v0AOD->DcaV0Daughters();
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }
      case kCosPointAng:
         if(v0esd) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   fComputedValue = v0ESD->GetV0CosineOfPointingAngle();
	   return kTRUE;	 
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   fComputedValue = TMath::Cos(v0AOD->OpenAngleV0());
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }
      case kLambdaProtonPIDCut:
         if(v0esd && lESDEvent) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   // retrieve the V0 daughters
   	   UInt_t lIdxPos      = (UInt_t) TMath::Abs(v0ESD->GetPindex());
   	   AliESDtrack *pTrack = lESDEvent->GetTrack(lIdxPos);
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, AliPID::kProton));
	   return kTRUE;	 
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   AliAODTrack *pTrack = (AliAODTrack *) (v0AOD->GetSecondaryVtx()->GetDaughter(0));
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, AliPID::kProton));
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }
       case kAntiLambdaAntiProtonPIDCut:
         if(v0esd && lESDEvent) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   // retrieve the V0 daughters
   	   UInt_t lIdxNeg      = (UInt_t) TMath::Abs(v0ESD->GetNindex());
   	   AliESDtrack *nTrack = lESDEvent->GetTrack(lIdxNeg);
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, AliPID::kProton));
	   return kTRUE;	 
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   AliAODTrack *nTrack = (AliAODTrack *) (v0AOD->GetSecondaryVtx()->GetDaughter(1));
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, AliPID::kProton));
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }
      case kLambdaPionPIDCut:
         if(v0esd && lESDEvent) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   // retrieve the V0 daughters
   	   UInt_t lIdxNeg      = (UInt_t) TMath::Abs(v0ESD->GetNindex());
   	   AliESDtrack *nTrack = lESDEvent->GetTrack(lIdxNeg);
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, AliPID::kPion));
	   return kTRUE;	 
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   AliAODTrack *nTrack = (AliAODTrack *) (v0AOD->GetSecondaryVtx()->GetDaughter(1));
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, AliPID::kPion));
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }
       case kAntiLambdaAntiPionPIDCut:
         if(v0esd && lESDEvent) {
	   AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
	   // retrieve the V0 daughters
   	   UInt_t lIdxPos      = (UInt_t) TMath::Abs(v0ESD->GetPindex());
   	   AliESDtrack *pTrack = lESDEvent->GetTrack(lIdxPos);
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, AliPID::kPion));
	   return kTRUE;	 
	 }
	 if(v0aod) {
	   AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	   AliAODTrack *pTrack = (AliAODTrack *) (v0AOD->GetSecondaryVtx()->GetDaughter(0));
	   AliPIDResponse *pid = fEvent->GetPIDResponse(); 
	   fComputedValue = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, AliPID::kPion));
	   return kTRUE;
	 }
	 else {
	  fComputedValue = -999;
	  return kFALSE;
	 }
	 
      
      default:
         AliError(Form("[%s] Invalid value type for this computation", GetName()));
         return kFALSE;
   }
}
