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
// modified: Kunal Garg (kgarg@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>


#include <TFormula.h>
#include <TBits.h>


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
      case kY:                         return "SingleTrackY";
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
      case kTOFdeltaPi:  	        return "SingleTrackTOFdeltaPion";
      case kTOFdeltaK:   	        return "SingleTrackTOFdeltaKaon";
      case kTOFdeltaP:   	        return "SingleTrackTOFdeltaProton";
      case kNITSclusters:               return "SingleTrackNITSclusters";
      case kNTPCclusters:               return "SingleTrackNTPCclusters";
      case kNTPCcrossedRows:            return "SingleTrackNTPCcrossedRows";
      case kNTPCcrossedRowsFclusters:   return "SingleTrackNTPCcrossedRowsFclusters";
      case kITSchi2:                    return "SingleTrackITSchi2";
      case kTPCchi2:                    return "SingleTrackTPCchi2";
      case kDCAXY:                      return "SingleTrackDCAz";
      case kDCAZ:                       return "SingleTrackDCAz";
      case kCharge:                     return "SingleTrackCharge";
      case kPhi:                        return "SingleTrackPhi";
      case kPhiOuterTPC:                return "SingleTrackPhiOuterTPC";
      case kV0DCA:        	        return "V0DCAToPrimaryVertex";
      case kV0Radius:        	        return "V0Radius";
      case kV0Mass:        	        return "V0Mass";
      case kV0P:        	        return "V0Momentum";
      case kV0Pt:        	        return "V0TransverseMomentum";
      case kV0NPt:        	        return "V0NegativeDaughterTransverseMomentum";
      case kV0PPt:        	        return "V0PositiveDaughterTransverseMomentum";
      case kV0DCAXY:                return "V0TracksDCAXY";
      case kV0Lifetime:             return "V0Lifetime";
      case kDaughterDCA:  	        return "V0DaughterDCA";
      case kCosPointAng:  	        return "V0CosineOfPointingAngle";
      case kLambdaProtonPIDCut:         return "V0LambdaProtonNsigma";
      case kAntiLambdaAntiProtonPIDCut: return "V0AntiLambdaAntiProtonNsigma";
      case kLambdaPionPIDCut:           return "V0LambdaPionNsigma";
      case kAntiLambdaAntiPionPIDCut:   return "V0AntiLambdaPionNsigma";
      case kK0SMass:                    return "K0SMass";
      case kLambdaMass:                 return "LambdaMass";
      case kAntiLambdaMass:             return "AntiLambdaMass";
      case kXiMass:                     return "XiMass";
      case kOmegaMass:                  return "OmegaMass";
      case kCascadeP:                   return "CascadeMomentum";
      case kCascadePt:                  return "CascadeTransverseMomentum";
      case kCascadeDCA:                 return "CascadeDCA";
      case kCascadeRadius:              return "CascadeRadius";
      case kCascadeDaughterDCA:         return "CascadeDaughterDCA";
      case kCascadeCosPointAng:         return "CascadeCosPointAng";
      case kCascadeV0CosPointAng:       return "CascadeV0CosPointAng";
      case kCascadeV0Lifetime:          return "CascadeV0Lifetime";
      case kBachelorPt:                 return "CascadeBachelorPt";
      case kBachelorPionTPCnsigma:      return "CascadeBachelorPionTPCnsigma";
      case kBachelorKaonTPCnsigma:      return "CascadeBachelorKaonTPCnsigma";
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
   AliESDcascade  *caesd = fDaughter->Ref2ESDcascade();
   AliAODcascade  *caaod = fDaughter->Ref2AODcascade();
   AliESDEvent    *lESDEvent = fEvent->GetRefESD();
   AliAODEvent    *lAODEvent = fEvent->GetRefAOD();
   
   Double_t xPrimaryVertex = -999.9;
   Double_t yPrimaryVertex = -999.9;
   Double_t zPrimaryVertex = -999.9;
   
   if(lESDEvent){
      xPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetX();
      yPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetY();
      zPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetZ();
   }else if(lAODEvent){
      xPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetX();
      yPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetY();
      zPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetZ();
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

   case kY:
     fComputedValue = (fUseMCInfo ? refMC->Y() : ref->Y());
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

   case kTOFdeltaPi:
     if (track) {
       AliPIDResponse *pid = fEvent->GetPIDResponse();
       pid->GetSignalDelta(AliPIDResponse::kTOF, track, AliPID::kPion, fComputedValue, kFALSE);//(==AliPIDResponse::kDetPidOk);
	 return kTRUE;
     } else {
       AliWarning("Cannot get TOF delta for non-track object");
       fComputedValue = 0.0;
       return kFALSE;
     }

   case kTOFdeltaK:
     if (track) {
       AliPIDResponse *pid = fEvent->GetPIDResponse();
       pid->GetSignalDelta(AliPIDResponse::kTOF, track, AliPID::kKaon, fComputedValue, kFALSE);//(==AliPIDResponse::kDetPidOk);
       return kTRUE;
     } else {
       AliWarning("Cannot get TOF delta for non-track object");
       fComputedValue = 0.0;
       return kFALSE;
     }

   case kTOFdeltaP:
     if (track) {
       AliPIDResponse *pid = fEvent->GetPIDResponse();
       pid->GetSignalDelta(AliPIDResponse::kTOF, track, AliPID::kProton, fComputedValue, kFALSE);//(==AliPIDResponse::kDetPidOk);
       return kTRUE;
     } else {
       AliWarning("Cannot get TOF delta for non-track object");
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

   case kNTPCcrossedRows:
     if (track) {
       AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
       if (trackESD) {
	 fComputedValue =  trackESD->GetTPCCrossedRows();
       } else {
	 fComputedValue =  ((AliAODTrack *)track)->GetTPCNCrossedRows();
       }
       return kTRUE;
     } else {
       AliWarning("Cannot get n TPC crossed rows for non-track object");
       fComputedValue = 0.0;
       return kFALSE;
     }

   case kNTPCcrossedRowsFclusters:
     if (track) {
       AliESDtrack *trackESD = dynamic_cast<AliESDtrack *>(track);
       fComputedValue = 1.0;
       if (trackESD) {
	 if (trackESD->GetTPCNclsF()>0) fComputedValue = trackESD->GetTPCCrossedRows() / trackESD->GetTPCNclsF();
       } else {
	 Float_t nCrossedRows = ((AliAODTrack*) track)->GetTPCNCrossedRows();
	 Float_t nFcls = ((AliAODTrack*) track)->GetTPCNclsF();
	 if (nFcls>0) fComputedValue = nCrossedRows / nFcls;
       }
       return kTRUE;
     } else {
       AliWarning("Cannot get n TPC crossed rows/findable clusters for non-track object");
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

   case kV0Radius:
     if(v0esd) {
       Double_t v0Position[3]; // from $ALICE_ROOT/ANALYSIS/AliESDV0Cuts.cxx
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       v0ESD->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);
       fComputedValue = TMath::Sqrt(TMath::Power(v0Position[0],2) + TMath::Power(v0Position[1],2));
       return kTRUE;
      }
     if(v0aod) {
        AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
	fComputedValue = v0AOD->RadiusV0();
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kV0Mass:
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       fComputedValue = v0ESD->GetEffMass();
       return kTRUE;
     }
     if(v0aod) {//returns the mass for the best hypothesis
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
       Double_t mass[3], dmass[3];
       mass[0] = v0AOD->MassK0Short();
       dmass[0] = TMath::Abs(mass[0] - 0.497611);
       mass[1] = v0AOD->MassLambda();
       dmass[1] = TMath::Abs(mass[1] - 1.115683);
       mass[2] = v0AOD->MassAntiLambda();
       dmass[2] = TMath::Abs(mass[2] - 1.115683);
       Int_t jmin=0;
       Double_t dmass_min = dmass[0];
       for(Int_t j=1; j<3; j++) if(dmass[j] < dmass_min){
         dmass_min = dmass[j];
         jmin = j;
       }
       fComputedValue = mass[jmin];
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }
	 
   case kV0P:
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       fComputedValue = v0ESD->P();
       return kTRUE;
     }
     if(v0aod) {
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
       fComputedValue = TMath::Sqrt(v0AOD->Ptot2V0());
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kV0Pt:
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       fComputedValue = v0ESD->Pt();
       return kTRUE;
     }
     if(v0aod) {
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
       fComputedValue = TMath::Sqrt(v0AOD->Pt2V0());
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kV0NPt:
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       Double_t px, py, pz;
       v0ESD->GetNPxPyPz(px, py, pz);
       fComputedValue = TMath::Sqrt(px*px+py*py);
       return kTRUE;
     }
     if(v0aod) {
       fComputedValue = -999;
       AliAODTrack *nTrack = (AliAODTrack *) (v0aod->GetDaughter(1));
       if(!nTrack) return kFALSE;
       fComputedValue=nTrack->Pt();
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kV0PPt:
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       Double_t px, py, pz;
       v0ESD->GetPPxPyPz(px, py, pz);
       fComputedValue = TMath::Sqrt(px*px+py*py);
       return kTRUE;
     }
     if(v0aod) {
       fComputedValue = -999;
       AliAODTrack *pTrack = (AliAODTrack *) (v0aod->GetDaughter(0));
       if(!pTrack) return kFALSE;
       fComputedValue=pTrack->Pt();
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kV0DCAXY:
     if(v0esd){
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       UInt_t lIdxNeg      = (UInt_t) TMath::Abs(v0ESD->GetNindex());
       AliESDtrack *nTrack = lESDEvent->GetTrack(lIdxNeg);
       if(!nTrack) return kFALSE;

       Float_t b[2], bCov[3];
       nTrack->GetImpactParameters(b, bCov);
       fComputedValue =  b[0];
       return kTRUE;
     }
     if(v0aod){
       fComputedValue = -999;
       AliAODTrack *nTrack = (AliAODTrack *) (v0aod->GetDaughter(1));
       if(!nTrack) return kFALSE;
	   Double_t b[2]= {-999,-999}, cov[3];
	   AliAODVertex *vertex = fEvent->GetRefAOD()->GetPrimaryVertex();
       if(!vertex) return kFALSE;
       nTrack->PropagateToDCA(vertex, fEvent->GetRefAOD()->GetMagneticField(),kVeryBig, b, cov);
       fComputedValue = b[0];
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kV0Lifetime:
     if(v0esd){
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       Double_t lMass = v0ESD->GetEffMass(); //Mass of V0 particle

       Double_t tV0mom[3];
       v0ESD->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
       Double_t lV0TotalMomentum = TMath::Sqrt(
                                        tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );  //Total Momentum of V0 particle

       Double_t v0Position[3];
       v0ESD->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);

       Double_t lLength = TMath::Sqrt(TMath::Power(v0Position[0]- xPrimaryVertex,2) + TMath::Power(v0Position[1] - yPrimaryVertex,2)+ TMath::Power(v0Position[2]- zPrimaryVertex,2));     //Distance of V0 from primary vertex

       fComputedValue = TMath::Abs(lMass*lLength/lV0TotalMomentum);
       return kTRUE;
     }
     if(v0aod) {
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);

       Double_t mass[3], dmass[3];
       mass[0] = v0AOD->MassK0Short();
       dmass[0] = TMath::Abs(mass[0] - 0.497611);
       mass[1] = v0AOD->MassLambda();
       dmass[1] = TMath::Abs(mass[1] - 1.115683);
       mass[2] = v0AOD->MassAntiLambda();
       dmass[2] = TMath::Abs(mass[2] - 1.115683);
       Int_t jmin=0;
       Double_t dmass_min = dmass[0];
       for(Int_t j=1; j<3; j++) if(dmass[j] < dmass_min){
         dmass_min = dmass[j];
         jmin = j;
       }
       Double_t lMass = mass[jmin];

       Double_t lV0TotalMomentum = TMath::Sqrt(v0AOD->Ptot2V0());
       Double_t v0Position[3] = {v0AOD->DecayVertexV0X(), v0AOD->DecayVertexV0Y(), v0AOD->DecayVertexV0Z()};
       Double_t lLength = TMath::Sqrt(TMath::Power(v0Position[0] - xPrimaryVertex,2) + TMath::Power(v0Position[1] - yPrimaryVertex,2)+ TMath::Power(v0Position[2] - zPrimaryVertex,2)); //Distance of V0 from primary vertex

       fComputedValue = TMath::Abs(lMass*lLength/lV0TotalMomentum);
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
       AliAODVertex *vertex = fEvent->GetRefAOD()->GetPrimaryVertex();
       fComputedValue = v0AOD->CosPointingAngle(vertex);
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

   case kK0SMass://V0 mass for K0S hypothesis
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       fComputedValue = v0ESD->GetEffMass(2,2);
       return kTRUE;
     }
     if(v0aod) {
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
       fComputedValue = v0AOD->MassK0Short();
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kLambdaMass://V0 mass for Lambda hypothesis
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       fComputedValue = v0ESD->GetEffMass(4,2);
       return kTRUE;
     }
     if(v0aod) {
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
       fComputedValue = v0AOD->MassLambda();
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kAntiLambdaMass://V0 mass for anti-Lambda hypothesis
     if(v0esd) {
       AliESDv0 *v0ESD = dynamic_cast<AliESDv0 *>(v0esd);
       fComputedValue = v0ESD->GetEffMass(2,4);
       return kTRUE;
     }
     if(v0aod) {
       AliAODv0 *v0AOD = dynamic_cast<AliAODv0 *>(v0aod);
       fComputedValue = v0AOD->MassAntiLambda();
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kXiMass://cascade mass for Xi hypothesis
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Int_t code = caESD->GetPdgCodeXi();
       Double_t v0q=0;
       if(code == kOmegaMinus) caESD->ChangeMassHypothesis(v0q, kXiMinus);
       if(code == kOmegaPlusBar) caESD->ChangeMassHypothesis(v0q, kXiPlusBar);
       fComputedValue = caESD->M();
       if(code == kOmegaMinus) caESD->ChangeMassHypothesis(v0q, kOmegaMinus);
       if(code == kOmegaPlusBar) caESD->ChangeMassHypothesis(v0q, kOmegaPlusBar);
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = caAOD->MassXi();
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kOmegaMass://cascade mass for Omega hypothesis
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Int_t code = caESD->GetPdgCodeXi();
       Double_t v0q=0;
       if(code == kXiMinus) caESD->ChangeMassHypothesis(v0q, kOmegaMinus);
       if(code == kXiPlusBar) caESD->ChangeMassHypothesis(v0q, kOmegaPlusBar);
       fComputedValue = caESD->M();
       if(code == kXiMinus) caESD->ChangeMassHypothesis(v0q, kXiMinus);
       if(code == kXiPlusBar) caESD->ChangeMassHypothesis(v0q, kXiPlusBar);
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = caAOD->MassOmega();
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }
	 
   case kCascadeP:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = caESD->P();
       return kTRUE;
     }
     if(caaod) {
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = TMath::Sqrt(caAOD->Ptot2Xi());
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadePt:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = caESD->Pt();
       return kTRUE;
     }
     if(caaod) {
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = TMath::Sqrt(caAOD->Pt2Xi());
       return kTRUE;
     }
     else {
       fComputedValue = -999;
       return kFALSE;
     }
           
   case kCascadeDCA:
     if(caesd && lESDEvent) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = caESD->GetDcascade(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex);
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = caAOD->DcaXiToPrimVertex();
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadeRadius:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Double_t XiPosition[3];
       caESD->GetXYZcascade(XiPosition[0],XiPosition[1],XiPosition[2]);
       fComputedValue = TMath::Sqrt(TMath::Power(XiPosition[0],2) + TMath::Power(XiPosition[1],2));
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = TMath::Sqrt(TMath::Power(caAOD->DecayVertexXiX(),2) + TMath::Power(caAOD->DecayVertexXiY(),2));
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadeDaughterDCA:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = caESD->GetDcaXiDaughters();
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = caAOD->DcaXiDaughters();
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadeCosPointAng:
     if(caesd && lESDEvent) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = caESD->GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
       return kTRUE;
     } else if(caaod && lAODEvent){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = caAOD->CosPointingAngleXi(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadeV0CosPointAng:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Double_t XiPosition[3];
       caESD->GetXYZcascade(XiPosition[0],XiPosition[1],XiPosition[2]);
       fComputedValue = caESD->GetV0CosineOfPointingAngle(XiPosition[0],XiPosition[1],XiPosition[2]);
       return kTRUE;
     } else if(caaod && lAODEvent){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       Double_t XiPosition[3];
       XiPosition[0] = caAOD->DecayVertexXiX();
       XiPosition[1] = caAOD->DecayVertexXiY();
       XiPosition[2] = caAOD->DecayVertexXiZ();
       fComputedValue = caAOD->CosPointingAngle(XiPosition);
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadeV0Lifetime:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Double_t XiPosition[3], v0Position[3], pmom[3], nmom[3];
       caESD->GetXYZcascade(XiPosition[0],XiPosition[1],XiPosition[2]);
       caESD->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);
       caESD->GetPPxPyPz(pmom[0], pmom[1], pmom[2]);
       caESD->GetNPxPyPz(nmom[0], nmom[1], nmom[2]);
       Double_t length = TMath::Sqrt(TMath::Power(v0Position[0] - XiPosition[0],2) + TMath::Power(v0Position[1] - XiPosition[1],2) + TMath::Power(v0Position[2] - XiPosition[2],2));
       Double_t v0mom = TMath::Sqrt(TMath::Power(pmom[0] + nmom[0],2) + TMath::Power(pmom[1] + nmom[1],2) + TMath::Power(pmom[2] + nmom[2],2));
       fComputedValue = TMath::Abs(1.115683*length/v0mom);
       return kTRUE;
     } else if(caaod && lAODEvent){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       Double_t length = TMath::Sqrt(TMath::Power(caAOD->DecayVertexV0X() - caAOD->DecayVertexXiX(),2) + TMath::Power(caAOD->DecayVertexV0Y() - caAOD->DecayVertexXiY(),2) + TMath::Power(caAOD->DecayVertexV0Z() - caAOD->DecayVertexXiZ(),2));
       Double_t v0mom = TMath::Sqrt(caAOD->Ptot2V0());
       fComputedValue = TMath::Abs(1.115683*length/v0mom);
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kCascadeV0Pt:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Double_t ppx, ppy, ppz, npx, npy, npz;
       caESD->GetPPxPyPz(ppx, ppy, ppz);
       caESD->GetNPxPyPz(npx, npy, npz);
       fComputedValue = TMath::Sqrt(TMath::Power(ppx+npx,2) + TMath::Power(ppy+npy,2));
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = -999;
       AliAODTrack *pTrack = (AliAODTrack *) (caAOD->GetDaughter(0));
       AliAODTrack *nTrack = (AliAODTrack *) (caAOD->GetDaughter(1));
       if(!pTrack || !nTrack) return kFALSE;
       fComputedValue=TMath::Sqrt(TMath::Power(pTrack->Px()+nTrack->Px(),2) + TMath::Power(pTrack->Py()+nTrack->Py(),2));
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kBachelorPt:
     if(caesd) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       Double_t px, py, pz;
       caESD->GetBPxPyPz(px, py, pz);
       fComputedValue = TMath::Sqrt(px*px+py*py);
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = -999;
       AliAODTrack* bTrack = (AliAODTrack *) (caAOD->GetDecayVertexXi()->GetDaughter(0));
       if(!bTrack) return kFALSE;
       fComputedValue=bTrack->Pt();
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kBachelorPionTPCnsigma:
     if(caesd && lESDEvent) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = -999;
       UInt_t lIdxBach     = (UInt_t) TMath::Abs(caESD->GetBindex());
       AliESDtrack* bTrack = lESDEvent->GetTrack(lIdxBach);
       if(!bTrack) return kFALSE;
       AliPIDResponse* pid = fEvent->GetPIDResponse();
       fComputedValue = pid->NumberOfSigmasTPC(bTrack, AliPID::kPion);
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = -999;
       AliAODTrack* bTrack = (AliAODTrack *) (caAOD->GetDecayVertexXi()->GetDaughter(0));
       if(!bTrack) return kFALSE;
       AliPIDResponse* pid = fEvent->GetPIDResponse();
       fComputedValue = pid->NumberOfSigmasTPC(bTrack, AliPID::kPion);
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   case kBachelorKaonTPCnsigma:
     if(caesd && lESDEvent) {
       AliESDcascade* caESD = dynamic_cast<AliESDcascade *>(caesd);
       fComputedValue = -999;
       UInt_t lIdxBach     = (UInt_t) TMath::Abs(caESD->GetBindex());
       AliESDtrack * bTrack = lESDEvent->GetTrack(lIdxBach);
       if(!bTrack) return kFALSE;
       AliPIDResponse* pid = fEvent->GetPIDResponse();
       fComputedValue = pid->NumberOfSigmasTPC(bTrack, AliPID::kKaon);
       return kTRUE;
     } else if(caaod){
       AliAODcascade* caAOD = dynamic_cast<AliAODcascade *>(caaod);
       fComputedValue = -999;
       AliAODTrack* bTrack = (AliAODTrack *) (caAOD->GetDecayVertexXi()->GetDaughter(0));
       if(!bTrack) return kFALSE;
       AliPIDResponse* pid = fEvent->GetPIDResponse();
       fComputedValue = pid->NumberOfSigmasTPC(bTrack, AliPID::kKaon);
       return kTRUE;
     } else {
       fComputedValue = -999;
       return kFALSE;
     }

   default:
     AliError(Form("[%s] Invalid value type for this computation", GetName()));
     return kFALSE;
   }
}
