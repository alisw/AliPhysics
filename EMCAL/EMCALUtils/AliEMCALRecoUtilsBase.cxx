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

// ROOT includes
#include <TGeoManager.h>
#include <TGeoBBox.h>

// STEER includes
#include "AliVCluster.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"

// EMCAL includes
#include "AliEMCALRecoUtilsBase.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliEMCALRecoUtilsBase) ;
/// \endcond

///
/// Constructor.
/// Initialize all constant values which have to be used
/// during Reco algorithm execution
///
//_____________________________________
AliEMCALRecoUtilsBase::AliEMCALRecoUtilsBase()
{
}

//
// Copy constructor.
//
//______________________________________________________________________
AliEMCALRecoUtilsBase::AliEMCALRecoUtilsBase(const AliEMCALRecoUtilsBase & reco) 
: TNamed(reco)
{  
}

///
/// Assignment operator.
///
//______________________________________________________________________
AliEMCALRecoUtilsBase & AliEMCALRecoUtilsBase::operator = (const AliEMCALRecoUtilsBase & reco) 
{  
  if (this == &reco)return *this;
  ((TNamed *)this)->operator=(reco);
  
  return *this;
}

///
/// Extrapolate track to EMCAL surface. 
/// NOTE, on success the call will change the track!
/// Used in AliAnalysisTaskESDfilter and AliEMCALRecoUtilsBase::ExtrapolateTrackToEMCalSurface(AliExternalTrackParam *trkParam, ...)
///
/// \param track: pointer to track
/// \param emcalR: distance
/// \param mass: mass hypothesis
/// \param step: ...
/// \param minpt: minimum track pT for matching
/// \param useMassForTracking: switch to use the mass hypothesis
/// \param useDCA: switch to use different track vertex
///
/// \return bool true if track could be extrapolated
///
//------------------------------------------------------------------------------------
Bool_t AliEMCALRecoUtilsBase::ExtrapolateTrackToEMCalSurface(AliVTrack *track,
                                                             Double_t emcalR, 
                                                             Double_t mass,
                                                             Double_t step, 
                                                             Double_t minpt,
                                                             Bool_t useMassForTracking, 
                                                             Bool_t useDCA,
                                                             Bool_t useOuterParam)
{ 
  track->SetTrackPhiEtaPtOnEMCal(-999, -999, -999);
  
  if ( track->Pt() < minpt )
    return kFALSE;
  
  if ( TMath::Abs(track->Eta()) > 0.9 ) return kFALSE;
  
  // Save some time and memory in case of no DCal present
  AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance();
  
  Int_t       nSupMod = AliEMCALGeoParams::fgkEMCALModules;
  if ( geom ) nSupMod = geom->GetNumberOfSuperModules();
  
  if ( nSupMod < 13 ) // Run1 10 (12, 2 not active but present)
  {
    Double_t phi = track->Phi()*TMath::RadToDeg();
    if ( phi <= 10 || phi >= 250 ) return kFALSE;
  }
  
  AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(track);
  AliAODTrack *aodt = 0;
  if (!esdt) 
  {
    aodt = dynamic_cast<AliAODTrack*>(track);
    
    if (!aodt)
      return kFALSE;
  }
  
  // Select the mass hypothesis
  if ( mass < 0 )
  {
    Bool_t onlyTPC = kFALSE;
    if ( mass == -99 ) onlyTPC=kTRUE;
    
    if (esdt)
    {
      if ( useMassForTracking ) mass = esdt->GetMassForTracking();
      else                      mass = esdt->GetMass(onlyTPC);
    }
    else
    {
      if ( useMassForTracking ) mass = aodt->GetMassForTracking();
      else                      mass = aodt->M();
    }
  }
  
  AliExternalTrackParam *trackParam = 0;
  if (esdt) 
  {
    AliExternalTrackParam *in = 0;

    if ( useOuterParam ) 
      in = const_cast<AliExternalTrackParam*>(esdt->GetOuterParam()) ;
    else                
      in = const_cast<AliExternalTrackParam*>(esdt->GetInnerParam()) ;
    
    if (!in)
      return kFALSE;
    
    trackParam = new AliExternalTrackParam(*in);
  } 
  else 
  {
    Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
    aodt->PxPyPz(pxpypz);
    if (useDCA) {
      // Note: useDCA is by default false in this function only for backwards compatibility
      aodt->GetXYZ(xyz);
    }
    else {
      aodt->XvYvZv(xyz);
    }
    aodt->GetCovarianceXYZPxPyPz(cv);  
    trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
  }
  
  if (!trackParam)
    return kFALSE;
  
  Float_t etaout=-999, phiout=-999, ptout=-999;
  Bool_t ret = ExtrapolateTrackToEMCalSurface(trackParam, 
                                              emcalR,
                                              mass,
                                              step,
                                              etaout, 
                                              phiout,
                                              ptout);
  
  delete trackParam;
  
  if (!ret) return kFALSE;
  
  if ( TMath::Abs(etaout) > 0.75 ) return kFALSE;
  
  // Save some time and memory in case of no DCal present
  if ( nSupMod < 13 ) // Run1 10 (12, 2 not active but present)
  {
    if ( (phiout < 70*TMath::DegToRad()) || (phiout > 190*TMath::DegToRad()) )  return kFALSE;
  }
  
  track->SetTrackPhiEtaPtOnEMCal(phiout, etaout, ptout);
  
  return kTRUE;
}


///
/// Extrapolate track to EMCAL surface.
/// Used in AliEMCALTracker.
///
/// \param trkParam: pointer to external track param
/// \param emcalR: distance
/// \param mass: mass hypothesis
/// \param step: ...
/// \param eta: track extrapolated eta
/// \param phi: track extrapolated phi
/// \param pt: track extrapolated pt
///
/// \return bool true if track could be extrapolated
///
//---------------------------------------------------------------------------------------
Bool_t AliEMCALRecoUtilsBase::ExtrapolateTrackToEMCalSurface(AliExternalTrackParam *trkParam, 
                                                             Double_t emcalR,
                                                             Double_t mass, 
                                                             Double_t step, 
                                                             Float_t &eta, 
                                                             Float_t &phi,
                                                             Float_t &pt)
{
  eta = -999, phi = -999, pt = -999;
  
  if (!trkParam) return kFALSE;
  
  if (!AliTrackerBase::PropagateTrackToBxByBz(trkParam, emcalR, mass, step, kTRUE, 0.8, -1)) return kFALSE;
  
  Double_t trkPos[3] = {0.,0.,0.};
  
  if (!trkParam->GetXYZ(trkPos)) return kFALSE;
  
  TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
  
  eta = trkPosVec.Eta();
  phi = trkPosVec.Phi();
  pt  = trkParam->Pt();
  
  if ( phi < 0 )
    phi += TMath::TwoPi();
  
  return kTRUE;
}

///
/// Return the residual by extrapolating a track param to a global position.
/// Used in AliEMCALTracker
///
/// \param trkParam: pointer to external track param
/// \param clsPos: cluster position array xyz
/// \param mass: mass hypothesis
/// \param step: ...
/// \param tmpEta: residual eta
/// \param tmpPhi: residual phi
///
/// \return bool true if track could be extrapolated
//-----------------------------------------------------------------------------------
Bool_t AliEMCALRecoUtilsBase::ExtrapolateTrackToPosition(AliExternalTrackParam *trkParam, 
                                                         const Float_t *clsPos, 
                                                         Double_t mass, 
                                                         Double_t step, 
                                                         Float_t &tmpEta, 
                                                         Float_t &tmpPhi)
{
  tmpEta = -999;
  tmpPhi = -999;
  
  if (!trkParam) return kFALSE;
  
  Double_t trkPos[3] = {0.,0.,0.};
  TVector3 vec(clsPos[0],clsPos[1],clsPos[2]);
  
  Float_t phi = vec.Phi();
  
  if ( phi < 0 )
    phi += TMath::TwoPi();
  
  // Rotate the cluster to the local extrapolation coordinate system
  Double_t alpha = ((int)(phi*TMath::RadToDeg()/20)+0.5)*20*TMath::DegToRad();
  
  vec.RotateZ(-alpha); 
  
  if (!AliTrackerBase::PropagateTrackToBxByBz(trkParam, vec.X(), mass, step,kTRUE, 0.8, -1)) 
    return kFALSE;
  
  if (!trkParam->GetXYZ(trkPos)) 
    return kFALSE; // Get the extrapolated global position
  
  TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);
  TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
  
  // Track cluster matching
  tmpPhi = clsPosVec.DeltaPhi(trkPosVec);    // tmpPhi is between -pi and pi
  tmpEta = clsPosVec.Eta()-trkPosVec.Eta();
  
  return kTRUE;
}

///
/// Return the residual by extrapolating a track param to a cluster.
/// Used in AliEMCALReconstructor.
///
/// \param trkParam: pointer to external track param
/// \param cluster: pointer to AliVCluster
/// \param mass: mass hypothesis
/// \param step: ...
/// \param tmpEta: residual eta
/// \param tmpPhi: residual phi
///
/// \return bool true if track could be extrapolated
///
//----------------------------------------------------------------------------------
Bool_t AliEMCALRecoUtilsBase::ExtrapolateTrackToCluster(AliExternalTrackParam *trkParam, 
                                                        const AliVCluster *cluster, 
                                                        Double_t mass, 
                                                        Double_t step, 
                                                        Float_t &tmpEta, 
                                                        Float_t &tmpPhi)
{
  tmpEta = -999;
  tmpPhi = -999;
  
  if (!cluster || !trkParam) 
    return kFALSE;
  
  Float_t clsPos[3] = {0.,0.,0.};
  cluster->GetPosition(clsPos);
  
  return ExtrapolateTrackToPosition(trkParam, clsPos, mass, step, tmpEta, tmpPhi);
}

///
/// Calculate shower depth for a given cluster energy and particle type.
/// Needed to calculate cluster position.
///
/// \param energy: cluster energy
/// \param iParticle: particle assumption defined in enum ParticleType
/// \param iSM: supermodule number
///
//_________________________________________________________
Float_t  AliEMCALRecoUtilsBase::GetDepth(Float_t energy, 
                                         Int_t   iParticle, 
                                         Int_t   iSM ) const 
{
  // parameters 
  Float_t x0    = 1.31;
  Float_t ecr   = 8;
  Float_t depth = 0;
  Float_t arg   = energy*1000/ ecr; //Multiply energy by 1000 to transform to MeV
  
  switch ( iParticle )
  {
    case kPhoton:
      if (arg < 1) 
        depth = 0;
      else
        depth = x0 * (TMath::Log(arg) + 0.5); 
      break;
      
    case kElectron:
      if (arg < 1) 
        depth = 0;
      else
        depth = x0 * (TMath::Log(arg) - 0.5); 
      break;
      
    case kHadron:
      // hadron 
      // boxes anc. here
      if (gGeoManager) 
      {
        gGeoManager->cd("ALIC_1/XEN1_1");
        TGeoNode        *geoXEn1    = gGeoManager->GetCurrentNode();
        TGeoNodeMatrix  *geoSM      = dynamic_cast<TGeoNodeMatrix *>(geoXEn1->GetDaughter(iSM));
        if (geoSM) 
        {
          TGeoVolume      *geoSMVol   = geoSM->GetVolume(); 
          TGeoShape       *geoSMShape = geoSMVol->GetShape();
          TGeoBBox        *geoBox     = dynamic_cast<TGeoBBox *>(geoSMShape);
          if (geoBox) depth = 0.5 * geoBox->GetDX()*2 ;
          else AliFatal("Null GEANT box");
        }
        else AliFatal("NULL  GEANT node matrix");
      }
      else
      {//electron
        if (arg < 1) 
          depth = 0;
        else
          depth = x0 * (TMath::Log(arg) - 0.5); 
      }
      
      break;
      
    default://photon
      if (arg < 1) 
        depth = 0;
      else
        depth = x0 * (TMath::Log(arg) + 0.5);
  }  
  
  return depth;
}

