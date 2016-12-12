/*************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliMCParticleContainerToyModel.h"
#include <TMath.h>
#include <TRandom3.h>
/// \cond CLASSIMP
ClassImp(AliMCParticleContainerToyModel);
/// \endcond

/**
 * This is the default constructor, used for ROOT I/O purposes.
 */
AliMCParticleContainerToyModel::AliMCParticleContainerToyModel() :
  AliMCParticleContainer(),
  fTrackScalePt(1),
  fTrackEtaWindow(0.9),
  fRandomizeEtaPhi(0)
 
 

{
}

/**
 * This is the standard named constructor.
 * \param name Name of the particle collection
 */
AliMCParticleContainerToyModel::AliMCParticleContainerToyModel(const char *name) :
  AliMCParticleContainer(name),
  fTrackScalePt(1),
  fTrackEtaWindow(0.9),
  fRandomizeEtaPhi(0) 
{
}

/**
 * Retrieve momentum information of a track and fill a TLorentzVector
 * with it. In case the optional parameter mass is provided, it is used as mass
 * hypothesis, otherwise the mass hypothesis from the particle itself is used.
 * Overrides the AliMCParticleContainer methods by applying a pT scaling factor to the TLorentzVector.
 * @param[out] mom Momentum vector to be filled
 * @param[in] track MCParticle from which the momentum information is obtained.
 * @param[in] mass (Optional) Mass hypothesis
 * @return
 */
Bool_t AliMCParticleContainerToyModel::GetMomentumFromParticle(TLorentzVector &mom, const AliAODMCParticle* track, Double_t mass) const
{
  Bool_t r = AliMCParticleContainer::GetMomentumFromParticle(mom, track, mass);
  ScalePtOfLorentzVector(mom);
  SetRandomEtaPhiOfLorentzVector(mom);
  return r;
}

/**
 * Fills a TLorentzVector with the momentum information of the
 * \f$ i^{th} \f$ particle in the container, using a global
 * mass hypothesis. In case the provided index is out of
 * range, false is returned as return value.
 * Overrides the AliMCParticleContainer methods by applying a pT scaling factor to the TLorentzVector.
 * @param[out] mom Momentum vector of the \f$ i^{th} \f$ particle in the array
 * @param[in] i Index of the particle to check
 * @return True if the request was successful, false otherwise
 */
Bool_t AliMCParticleContainerToyModel::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  Bool_t r = AliMCParticleContainer::GetMomentum(mom, i);
  ScalePtOfLorentzVector(mom);
  SetRandomEtaPhiOfLorentzVector(mom);
  return r;
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * \f$ i^{th} \f$ accepted particle in the container, using a
 * global mass hypothesis. In case the provided index is out of
 * range, or the particle under the index is not accepted, false
 * is returned as return value.
 * Overrides the AliMCParticleContainer methods by applying a pT scaling factor to the TLorentzVector.
 * @param[out] mom Momentum vector of the accepted particle
 * @param[in] i Index to check
 * @return True if the request was successfull, false otherwise
 */
Bool_t AliMCParticleContainerToyModel::GetAcceptMomentum(TLorentzVector &mom, Int_t i) const
{
  Bool_t r = AliMCParticleContainer::GetAcceptMomentum(mom, i);
  ScalePtOfLorentzVector(mom);
  SetRandomEtaPhiOfLorentzVector(mom);
  return r;
}

/**
 * Fills a TLorentzVector with the momentum information of the
 * next particle in the container, using a global mass hypothesis.
 * In case the iterator reached the end of the array, false
 * is returned as return value.
 * Overrides the AliMCParticleContainer methods by applying a pT scaling factor to the TLorentzVector.
 * @deprecated Old style iterator - use all_iterator instead
 * @param[out] mom Momentum vector of the next particle
 * @return True if the request was successful, false otherwise
 */
Bool_t AliMCParticleContainerToyModel::GetNextMomentum(TLorentzVector &mom)
{
  Bool_t r = AliMCParticleContainer::GetNextMomentum(mom);
  ScalePtOfLorentzVector(mom);
  SetRandomEtaPhiOfLorentzVector(mom);
  return r;
}

/**
 * Fills a TLorentzVector with the monentum infomation of the
 * next accepted particle in the container, using a global
 * mass hypothesis. In case the iteration reached the end of
 * the array, false is returned as return value.
 * Overrides the AliMCParticleContainer methods by applying a pT scaling factor to the TLorentzVector.
 * @deprecated Old style iterator - use accept_iterator instead
 * @param[out] mom Momentum vector of the next particle in the array
 * @return True if the request was successfull, false (no more entries) otherwise
 */
Bool_t AliMCParticleContainerToyModel::GetNextAcceptMomentum(TLorentzVector &mom)
{
  Bool_t r = AliMCParticleContainer::GetNextAcceptMomentum(mom);
  ScalePtOfLorentzVector(mom);
  SetRandomEtaPhiOfLorentzVector(mom);
  return r;
}

/**
 * Scales the pt of a TLorentzVector with a constant factor.
 * @param mom TLorentzVector object reference to be scaled.
 */
void AliMCParticleContainerToyModel::ScalePtOfLorentzVector(TLorentzVector &mom) const
{
  if(fTrackScalePt<1.){
   
    Double_t pTscale = fTrackScalePt*mom.Pt();
    Double_t phiscale   = mom.Phi();
    Double_t thetascale = 2.*TMath::ATan(TMath::Exp(-1.*(mom.Eta())));
    Double_t pXscale    = pTscale * TMath::Cos(phiscale);
    Double_t pYscale    = pTscale * TMath::Sin(phiscale);
    Double_t pZscale    = pTscale/TMath::Tan(thetascale);
    Double_t pscale=TMath::Sqrt(pTscale*pTscale+pZscale*pZscale);
    mom.SetPxPyPzE(pXscale, pYscale, pZscale, pscale);
  }
}
/**
 * Assigns random phi,eta to thetracks,keeping their momentum
 * @param mom TLorentzVector object reference to be scaled.
 */
void AliMCParticleContainerToyModel::SetRandomEtaPhiOfLorentzVector(TLorentzVector &mom) const
{
  if(fRandomizeEtaPhi==1){
   
   
      
    Double_t pTscale = mom.Pt();
 
    Double_t etascale = 2.*fTrackEtaWindow * gRandom->Rndm() - fTrackEtaWindow;
    Double_t phiscale = 2.* TMath::Pi() * gRandom->Rndm();

    
    Double_t thetascale = 2.*TMath::ATan(TMath::Exp(-1.*(etascale)));
    Double_t pXscale    = pTscale * TMath::Cos(phiscale);
    Double_t pYscale    = pTscale * TMath::Sin(phiscale);
    Double_t pZscale    = pTscale/TMath::Tan(thetascale);
    Double_t pscale=TMath::Sqrt(pTscale*pTscale+pZscale*pZscale);
    mom.SetPxPyPzE(pXscale, pYscale, pZscale, pscale);
    
  }
}

void AliMCParticleContainerToyModel::ExecOnce()
{
  
    if (gRandom) delete gRandom;
     gRandom = new TRandom3(0);
  }
