/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliCaloTrackParticleCorrelation.h"
#include "AliAODJet.h"

/// \cond CLASSIMP
ClassImp(AliCaloTrackParticleCorrelation)
/// \endcond

//______________________________________________________________________________
///
/// Default Constructor.
///
 AliCaloTrackParticleCorrelation::AliCaloTrackParticleCorrelation() :
   AliCaloTrackParticle(),
   fLeadingDetector(-1), fLeading(), fCorrJet(),  fCorrBkg(), fRefJet(0),
   fListOfObjArrays(0)
{
}

//______________________________________________________________________________
///
/// Constructor.
///
/// \param px particle momentum in x
/// \param py particle momentum in y
/// \param pz particle momentum in z
/// \param e particle energy
///
/// particle: cluster or track
///
AliCaloTrackParticleCorrelation::AliCaloTrackParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e):
  AliCaloTrackParticle(),
  fLeadingDetector(-1),  fLeading(), fCorrJet(),
  fCorrBkg(), fRefJet(0),  fListOfObjArrays(new TList)
{
  SetMomentum(new TLorentzVector(px, py, pz, e));
  
  fListOfObjArrays->SetOwner(kTRUE);
}

//______________________________________________________________________________
///
/// Constructor.
///
/// \param p: TLorentzVector of particle kinematics.
///
/// particle: cluster or track
///
AliCaloTrackParticleCorrelation::AliCaloTrackParticleCorrelation(TLorentzVector & p):
  AliCaloTrackParticle(p),
  fLeadingDetector(-1),  fLeading(), fCorrJet(), fCorrBkg(), fRefJet(0), fListOfObjArrays(new TList)
{
  fListOfObjArrays->SetOwner(kTRUE);
}

//______________________________________________________________________________
///
/// Constructor.
///
/// \param p: AliCaloTrackParticle of particle kinematics and detector other inputs.
///
/// particle: cluster or track
///
AliCaloTrackParticleCorrelation::AliCaloTrackParticleCorrelation(AliCaloTrackParticle & p):
  AliCaloTrackParticle(p),
  fLeadingDetector(-1),  fLeading(), fCorrJet(), fCorrBkg(),fRefJet(0), fListOfObjArrays(new TList)
{
  fListOfObjArrays->SetOwner(kTRUE);
}

//______________________________________________________________________________
///
/// Destructor.
///
AliCaloTrackParticleCorrelation::~AliCaloTrackParticleCorrelation() 
{
  if(fListOfObjArrays)
  {
    fListOfObjArrays->Clear();
    delete   fListOfObjArrays ;
  }
}

//______________________________________________________________________________
///
/// Clear object.
///
void AliCaloTrackParticleCorrelation::Clear(const Option_t* /*opt*/) 
{  
  AliCaloTrackParticle::Clear(""); //delete fMomentum
  
  if(fListOfObjArrays)
  {
    fListOfObjArrays->Clear();
    delete   fListOfObjArrays ;
  }
}


//______________________________________________________________________________
///
/// Copy constructor
///
AliCaloTrackParticleCorrelation::AliCaloTrackParticleCorrelation(const AliCaloTrackParticleCorrelation& part) :
  AliCaloTrackParticle(part), fLeadingDetector(part.fLeadingDetector), fLeading(part.fLeading),  
  fCorrJet(part.fCorrJet), fCorrBkg(part.fCorrBkg), fRefJet(part.fRefJet),   
  fListOfObjArrays(new TList)
{
}

//______________________________________________________________________________
///
/// Assignment operator
///
//AliCaloTrackParticleCorrelation& AliCaloTrackParticleCorrelation::operator=(const AliCaloTrackParticleCorrelation& part)
//{
//  if(this!=&part) 
//  {
//    fRefJet   = part.fRefJet ;
//    fLeading  = part.fLeading;
//    fCorrJet  = part.fCorrJet ;
//    fCorrBkg  = part.fCorrBkg; 
//    fListOfObjArrays = fListOfObjArrays;
//  }
//  return *this;
//}

//______________________________________________________________________________
///
/// Print information of all data members.
///
void AliCaloTrackParticleCorrelation::Print(Option_t* /*option*/) const 
{
  AliCaloTrackParticle::Print("");

  if(GetJet()) GetJet()->Print("");

  printf("Leading Detector : %d\n",fLeadingDetector);
  printf("Leading Particle 4-vector:\n");
  printf("     E  = %13.3f",   fLeading.E() );
  printf("     Px = %13.3f",   fLeading.Px());
  printf("     Py = %13.3f",   fLeading.Py());
  printf("     Pz = %13.3f\n", fLeading.Pz());

  if( fListOfObjArrays)   fListOfObjArrays->Print("");
}
