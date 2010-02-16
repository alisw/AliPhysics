#ifndef ALIDIELECTRONVARMANAGER_H
#define ALIDIELECTRONVARMANAGER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronVarManager                     #
//#         Class for management of available variables       #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################


#include <TNamed.h>

#include <AliVEvent.h>
#include <AliESDEvent.h>

#include <AliVParticle.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliKFParticle.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include <AliVTrack.h>  // ?

#include "AliDielectronPair.h"

class AliVEvent;

//________________________________________________________________
class AliDielectronVarManager : public TNamed {
  
public:
  
  // Particle specific variables
  enum ValueTypes {
    kPx = 0,                 // px
    kPy,                     // py
    kPz,                     // pz
    kPt,                     // transverse momentum
    kP,                      // momentum
    kXv,                     // vertex position in x
    kYv,                     // vertex position in y
    kZv,                     // vertex position in z
    kOneOverPt,              // 1/pt
    kPhi,                    // phi angle
    kTheta,                  // theta angle
    kEta,                    // pseudo-rapidity
    kY,                      // rapidity
    kE,                      // energy
    kM,                      // mass
    kCharge,                 // charge
    kNclsITS,                // number of clusters assigned in the ITS
    kNclsTPC,                // number of clusters assigned in the TPC
    kNFclsTPC,               // number of findable clusters in the TPC
    kTPCsignalN,             // number of points used for dEdx
    kNclsTRD,                // number of clusters assigned in the TRD
    kTRDntracklets,          // number of TRD tracklets used for tracking/PID TODO: correct getter
    kTRDpidQuality,          // number of TRD tracklets used for PID
    kImpactParXY,            // Impact parameter in XY plane
    kImpactParZ,             // Impact parameter in Z
    kTrackLength,            // Track length
    kPdgCode,                // PDG code
    kPIn,                    // momentum at inner wall of TPC (if available), used for PID
    kTPCsignal,              // TPC dE/dx signal
    kParticleMax,             //
    // TODO: kRNClusters ??
  // AliDielectronPair specific variables
    kChi2NDF = kParticleMax, // Chi^2/NDF
    kDecayLength,            // decay length
    kR,                      // distance to the origin
    kOpeningAngle,           // opening angle
    kMerr,                   // error of mass calculation
    kDCA,                    // distance of closest approach TODO: not implemented yet
    kPairType,               // type of the pair, like like sign ++ unlikesign ...
    kPairMax,                 //
  // Event specific variables
    kXvPrim=kPairMax,        // TODO: prim vertex
    kYvPrim,                 // TODO: prim vertex
    kZvPrim,                 // TODO: prim vertex
    kXRes,                   // primary vertex x-resolution (AliESDVertex)
    kYRes,                   // primary vertex y-resolution (AliESDVertex)
    kZRes,                   // primary vertex z-resolution (AliESDVertex)
    kNTrk,                   // number of tracks (or tracklets)
    kTracks,                 // ESD tracks (AliESDEvent)
    kNMaxValues              //
    // TODO: (for A+A) ZDCEnergy, impact parameter, Iflag??
  };
  

  AliDielectronVarManager();
  AliDielectronVarManager(const char* name, const char* title);
  virtual ~AliDielectronVarManager();
  static void Fill(const TObject* particle, Double_t * const values);

  static const char* GetValueName(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i]:""; }
private:

  static const char* fgkParticleNames[kNMaxValues];  //variable names

  static void FillVarVParticle(const AliVParticle *particle,         Double_t * const values);
  static void FillVarESDtrack(const AliESDtrack *particle,           Double_t * const values);
  static void FillVarAODTrack(const AliAODTrack *particle,           Double_t * const values);
  static  void FillVarMCParticle(const AliMCParticle *particle,       Double_t * const values);
  static void FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values);
  static void FillVarDielectronPair(const AliDielectronPair *pair,   Double_t * const values);

  AliDielectronVarManager(const AliDielectronVarManager &c);
  AliDielectronVarManager &operator=(const AliDielectronVarManager &c);
  
  ClassDef(AliDielectronVarManager,1);
};


//Inline functions
inline void AliDielectronVarManager::Fill(const TObject* particle, Double_t * const values)
{
  //
  // Main function to fill all available variables according to the type of particle
  //

  if      (particle->IsA() == AliESDtrack::Class())       FillVarESDtrack(static_cast<const AliESDtrack*>(particle), values);
  else if (particle->IsA() == AliAODTrack::Class())       FillVarAODTrack(static_cast<const AliAODTrack*>(particle), values);
  else if (particle->IsA() == AliMCParticle::Class())     FillVarMCParticle(static_cast<const AliMCParticle*>(particle), values);
  else if (particle->IsA() == AliAODMCParticle::Class())  FillVarAODMCParticle(static_cast<const AliAODMCParticle*>(particle), values);
  else if (particle->IsA() == AliDielectronPair::Class()) FillVarDielectronPair(static_cast<const AliDielectronPair*>(particle), values);
//   else Error("Fill",Form("Type %s is not supported by AliDielectronVarManager!", particle->ClassName())); //TODO: implement without object needed
}

inline void AliDielectronVarManager::FillVarVParticle(const AliVParticle *particle, Double_t * const values)
{
  //
  // Fill track information available in AliVParticle into an array
  //
  values[AliDielectronVarManager::kPx]        = particle->Px();
  values[AliDielectronVarManager::kPy]        = particle->Py();
  values[AliDielectronVarManager::kPz]        = particle->Pz();
  values[AliDielectronVarManager::kPt]        = particle->Pt();
  values[AliDielectronVarManager::kP]         = particle->P();

  values[AliDielectronVarManager::kXv]        = particle->Xv();
  values[AliDielectronVarManager::kYv]        = particle->Yv();
  values[AliDielectronVarManager::kZv]        = particle->Zv();

  values[AliDielectronVarManager::kOneOverPt] = particle->OneOverPt();
  values[AliDielectronVarManager::kPhi]       = particle->Phi();
  values[AliDielectronVarManager::kTheta]     = particle->Theta();
  values[AliDielectronVarManager::kEta]       = particle->Eta();
  values[AliDielectronVarManager::kY]         = particle->Y();

  values[AliDielectronVarManager::kE]         = particle->E();
  values[AliDielectronVarManager::kM]         = particle->M();
  values[AliDielectronVarManager::kCharge]    = particle->Charge();
}

inline void AliDielectronVarManager::FillVarESDtrack(const AliESDtrack *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  // Fill AliESDtrack interface specific information
  values[AliDielectronVarManager::kNclsITS]       = particle->GetNcls(0); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNclsTPC]       = particle->GetNcls(1); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNFclsTPC]      = particle->GetTPCNclsF();
  values[AliDielectronVarManager::kTPCsignalN]    = particle->GetTPCsignalN();
  values[AliDielectronVarManager::kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
  values[AliDielectronVarManager::kTRDpidQuality] = particle->GetTRDpidQuality();

  Float_t impactParXY, impactParZ;
  particle->GetImpactParameters(impactParXY, impactParZ);
  values[AliDielectronVarManager::kImpactParXY]   = impactParXY;
  values[AliDielectronVarManager::kImpactParZ]    = impactParZ;

  values[AliDielectronVarManager::kTrackLength]   = particle->GetIntegratedLength();

  //dEdx information
  Double_t mom = particle->GetP();
  const AliExternalTrackParam *in=particle->GetInnerParam();
  if (in) mom = in->GetP();
  values[AliDielectronVarManager::kPIn]=mom;
  values[AliDielectronVarManager::kTPCsignal]=particle->GetTPCsignal();
  
}

inline void AliDielectronVarManager::FillVarAODTrack(const AliAODTrack *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  // Fill AliAODTrack interface information
  // ...

}

inline void AliDielectronVarManager::FillVarMCParticle(const AliMCParticle *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  // Fill AliMCParticle interface specific information
  values[AliDielectronVarManager::kPdgCode] = particle->PdgCode();
}

inline void AliDielectronVarManager::FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  // Fill AliAODMCParticle interface specific information
  values[AliDielectronVarManager::kPdgCode] = particle->GetPdgCode();
}

inline void AliDielectronVarManager::FillVarDielectronPair(const AliDielectronPair *pair, Double_t * const values)
{
  //
  // Fill pair information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(pair, values);

  // Fill AliDielectronPair specific information
  const AliKFParticle &kfPair = pair->GetKFParticle();
  
  values[AliDielectronVarManager::kChi2NDF]      = kfPair.GetChi2()/kfPair.GetNDF();
  values[AliDielectronVarManager::kDecayLength]  = kfPair.GetDecayLength();
  values[AliDielectronVarManager::kR]            = kfPair.GetR();
  values[AliDielectronVarManager::kOpeningAngle] = pair->OpeningAngle();
  values[AliDielectronVarManager::kMerr]         = kfPair.GetErrMass()>0?kfPair.GetErrMass()/kfPair.GetMass():1000000;
  values[AliDielectronVarManager::kPairType]     = pair->GetType();
}

/*
inline void AliDielectronVarManager::FillValues(const TParticle *particle, Double_t *values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill TParticle interface information
  values[AliDielectronVarManager::kPx]     = particle->Px();
  values[AliDielectronVarManager::kPy]     = particle->Py();
  values[AliDielectronVarManager::kPz]     = particle->Pz();
  values[AliDielectronVarManager::kPt]     = particle->Pt();
  values[AliDielectronVarManager::kP]      = particle->P();

  values[AliDielectronVarManager::kXv]     = particle->Vx();
  values[AliDielectronVarManager::kYv]     = particle->Vy();
  values[AliDielectronVarManager::kZv]     = particle->Vz();

  values[AliDielectronVarManager::kOneOverPt] = 1./particle->Pt();
  values[AliDielectronVarManager::kPhi]    = particle->Phi();
  values[AliDielectronVarManager::kTheta]  = 
  values[AliDielectronVarManager::kEta]    = particle->Eta();
  values[AliDielectronVarManager::kY]      = 

  values[AliDielectronVarManager::kE]      = particle->Energy();
  values[AliDielectronVarManager::kM]      = particle->GetMass();

  values[AliDielectronVarManager::kCharge] = particle->GetPDG()->Charge()/3; // uggly

}*/

#endif
