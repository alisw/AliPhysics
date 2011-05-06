#ifndef ALIDIELECTRONVARCONTAINER_H
#define ALIDIELECTRONVARCONTAINER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronVarContainer                   #
//#         Class for management of available variables       #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Markus    Köhler,   GSI / M.Koehler@gsi.de              #
//#   Frederick Kramer,   Uni Ffm / Frederick.Kramer@cern.ch  #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################


#include <TNamed.h>
#include <TDatabasePDG.h>

#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliESDVertex.h>
#include <AliAODVertex.h>

#include <AliVParticle.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliAODPid.h>
#include <AliKFParticle.h>
#include <AliKFVertex.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include <AliVTrack.h>  // ?

#include <AliExternalTrackParam.h>
#include <AliESDpid.h>
#include <AliCentrality.h>
#include <AliAODpidUtil.h>
#include <AliPID.h>
#include <AliPIDResponse.h>

#include "AliDielectronPair.h"
#include "AliDielectronMC.h"
#include "AliDielectronPID.h"
#include "AliDielectronHelper.h"

class AliVEvent;

//________________________________________________________________
class AliDielectronVarContainer : public TNamed {
  
public:

  enum ValueTypes {
    // Leg specific variables
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
    kNclsTPCiter1,           // number of clusters assigned in the TPC after first iteration
    kNFclsTPC,               // number of findable clusters in the TPC
    kNFclsTPCr,              // number of findable clusters in the TPC with more robust definition
    kNFclsTPCrFrac,          // number of found/findable clusters in the TPC with more robust definition
    kTPCsignalN,             // number of points used for dEdx
    kTPCsignalNfrac,         // fraction of points used for dEdx / cluster used for tracking
    kTPCchi2Cl,              // chi2/cl in TPC
    kTrackStatus,            // track status bits
    kNclsTRD,                // number of clusters assigned in the TRD
    kTRDntracklets,          // number of TRD tracklets used for tracking/PID TODO: correct getter
    kTRDpidQuality,          // number of TRD tracklets used for PID
    kTRDprobEle,             // TRD electron pid probability
    kTRDprobPio,             // TRD electron pid probability
    kImpactParXY,            // Impact parameter in XY plane
    kImpactParZ,             // Impact parameter in Z
    kTrackLength,            // Track length
    kITSsignal,		           // ITS dE/dx signal
    kITSsignalSSD1,	         // SSD1 dE/dx signal
    kITSsignalSSD2,	         // SSD2 dE/dx signal
    kITSsignalSDD1,	         // SDD1 dE/dx signal
    kITSsignalSDD2,	         // SDD2 dE/dx signal
    kITSclusterMap,          // ITS cluster map
    kITSnSigmaEle,           // number of sigmas to the dE/dx electron line in the ITS
    kITSnSigmaPio,           // number of sigmas to the dE/dx pion line in the ITS
    kITSnSigmaMuo,           // number of sigmas to the dE/dx muon line in the ITS
    kITSnSigmaKao,           // number of sigmas to the dE/dx kaon line in the ITS
    kITSnSigmaPro,           // number of sigmas to the dE/dx proton line in the ITS
    kPIn,                    // momentum at inner wall of TPC (if available), used for PID
    kTPCsignal,              // TPC dE/dx signal
  	kTOFsignal,              // TOF signal
	  kTOFbeta,                // TOF beta
    kTPCnSigmaEle,           // number of sigmas to the dE/dx electron line in the TPC
    kTPCnSigmaPio,           // number of sigmas to the dE/dx pion line in the TPC
    kTPCnSigmaMuo,           // number of sigmas to the dE/dx muon line in the TPC
    kTPCnSigmaKao,           // number of sigmas to the dE/dx kaon line in the TPC
    kTPCnSigmaPro,           // number of sigmas to the dE/dx proton line in the TPC
	  kTOFnSigmaEle,           // number of sigmas to the pion line in the TOF
    kTOFnSigmaPio,           // number of sigmas to the pion line in the TOF
    kTOFnSigmaMuo,           // number of sigmas to the muon line in the TOF
    kTOFnSigmaKao,           // number of sigmas to the kaon line in the TOF
    kTOFnSigmaPro,           // number of sigmas to the proton line in the TOF
    kKinkIndex0,             // kink index 0
    kParticleMax,            // TODO: kRNClusters ??

    // Pair specific variables
    kChi2NDF = kParticleMax, // Chi^2/NDF
    kDecayLength,            // decay length
    kR,                      // distance to the origin
    kOpeningAngle,           // opening angle
    // Helicity picture: Z-axis is considered the direction of the mother's 3-momentum vector
    kThetaHE,                // theta in mother's rest frame in the helicity picture 
    kPhiHE,                  // phi in mother's rest frame in the helicity picture
    // Collins-Soper picture: Z-axis is considered the direction of the vectorial difference between 
    // the 3-mom vectors of target and projectile beams
    kThetaCS,                // theta in mother's rest frame in Collins-Soper picture
    kPhiCS,                  // phi in mother's rest frame in Collins-Soper picture
    kLegDist,                // distance of the legs
    kLegDistXY,              // distance of the legs in XY
    kDeltaEta,               // Absolute value of Delta Eta for the legs
    kDeltaPhi,               // Absolute value of Delta Phi for the legs
    kMerr,                   // error of mass calculation
    kDCA,                    // distance of closest approach TODO: not implemented yet
    kPairType,               // type of the pair, like like sign ++ unlikesign ...
    kPseudoProperTime,       // pseudo proper time
    kPairMax,                //

    // Event specific variables
    kXvPrim = kPairMax,      // prim vertex x
    kYvPrim,                 // prim vertex y
    kZvPrim,                 // prim vertex z
    kXRes,                   // primary vertex x-resolution
    kYRes,                   // primary vertex y-resolution
    kZRes,                   // primary vertex z-resolution
    kNContrib,               // Number of vertex contributors
    kBzkG,                   // z componenent of field in kG
    kNTrk,                   // number of tracks
    kNacc,                   // Number of accepted tracks
    kNaccTrcklts,            // Number of accepted tracklets (MUON definition)
    kCentrality,             // event centrality fraction
    kNevents,                // event counter
    kNMaxValues              // TODO: (for A+A) ZDCEnergy, impact parameter, Iflag??
  };


  enum ValueTypesMC {
    // Leg specific variables
    kPx_MC = 0,              // px
    kPy_MC,                  // py
    kPz_MC,                  // pz
    kPt_MC,                  // transverse momentum
    kP_MC,                   // momentum
    kXv_MC,                  // vertex position in x
    kYv_MC,                  // vertex position in y
    kZv_MC,                  // vertex position in z
    kOneOverPt_MC,           // 1/pt
    kPhi_MC,                 // phi angle
    kTheta_MC,               // theta angle
    kEta_MC,                 // pseudo-rapidity
    kY_MC,                   // rapidity
    kE_MC,                   // energy
    kM_MC,                   // mass
    kCharge_MC,              // charge
    kImpactParXY_MC,         // Impact parameter in XY plane
    kImpactParZ_MC,          // Impact parameter in Z
    kPdgCode,                // PDG code
    kPdgCodeMother,          // PDG code of the mother
    kPdgCodeGrandMother,     // PDG code of the grand mother
    kNumberOfDaughters,      // number of daughters
    kHaveSameMother,         // check that particles have the same mother (MC)
    kIsJpsiPrimary,          // check if the particle is primary (MC)
    kParticleMax_MC,         //

// Pair specific variables
    kDecayLengthv_MC = kParticleMax_MC, // decay length
    kR_MC,                      // distance to the origin
    kOpeningAngle_MC,           // opening angle
    // Helicity picture: Z-axis is considered the direction of the mother's 3-momentum vector
    kThetaHE_MC,                // theta in mother's rest frame in the helicity picture 
    kPhiHE_MC,                  // phi in mother's rest frame in the helicity picture
    // Collins-Soper picture: Z-axis is considered the direction of the vectorial difference between 
    // the 3-mom vectors of target and projectile beams
    kThetaCS_MC,                // theta in mother's rest frame in Collins-Soper picture
    kPhiCS_MC,                  // phi in mother's rest frame in Collins-Soper picture
    kLegDist_MC,                // distance of the legs
    kLegDistXY_MC,              // distance of the legs in XY
    kDeltaEta_MC,               // Absolute value of Delta Eta for the legs
    kDeltaPhi_MC,               // Absolute value of Delta Phi for the legs
    kDCA_MC,                    // distance of closest approach TODO: not implemented yet
    kPairType_MC,               // type of the pair, like like sign ++ unlikesign ...
    kPseudoProperTime_MC,       // pseudo proper time
    kPairMax_MC,                //

// Event specific variables
    kXvPrim_MC = kPairMax_MC,   // prim vertex x
    kYvPrim_MC,                 // prim vertex y
    kZvPrim_MC,                 // prim vertex z
    kNch,                       // Number of charged MC tracks
    kCentrality_MC,             // event centrality fraction
    kNevents_MC,                // event counter
    kNMaxValues_MC              // TODO: (for A+A) ZDCEnergy, impact parameter, Iflag??
  };



  AliDielectronVarContainer();
  AliDielectronVarContainer(const char* name, const char* title);
  virtual ~AliDielectronVarContainer();

  static void Fill(const TObject* object);
  static void InitESDpid(Int_t type=0);
  static void InitAODpidUtil(Int_t type=0);

  static void SetESDpid(AliESDpid * const pid)            { fgPIDResponse=pid;                                   }
  static void SetPIDResponse(AliPIDResponse *pidResponse) { fgPIDResponse=pidResponse;                           }
  static void SetEvent(AliVEvent * const event); // TODO: needed?

  static AliESDpid* GetESDpid()                           { return (AliESDpid*)fgPIDResponse;                    }
  static AliAODpidUtil* GetAODpidUtil()                   { return (AliAODpidUtil*)fgPIDResponse;                }
  static const AliKFVertex* GetKFVertex()                 { return fgKFVertex;                                   }
  static const char* GetValueName(Int_t i)                { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i]:""; }
  static const Double_t* GetData()                        { return fgData;                                       }
  static const Double_t* GetDataMC()                      { return fgDataMC;                                     }
  static Double_t GetValue(ValueTypes val)                { return fgData[val];                                  }
  static Double_t GetValueMC(ValueTypesMC val)            { return fgDataMC[val];                                }
  static Bool_t GetDCA(const AliAODTrack *track, Double_t d0z0[2]);

private:

  static Double_t fgData[kNMaxValues];                   //! data
  static Double_t fgDataMC[kNMaxValues_MC];              //! MC data
  static const char* fgkParticleNames[kNMaxValues];      // variable names
  static const char* fgkParticleNamesMC[kNMaxValues_MC]; // MC variable names

  static void FillVarVParticle(const AliVParticle *particle);
  static void FillVarVParticleMC(const AliVParticle *particle);
  static void FillVarESDtrack(const AliESDtrack *particle);
  static void FillVarAODTrack(const AliAODTrack *particle);
  static void FillVarMCParticle(const AliMCParticle *particle);
  static void FillVarAODMCParticle(const AliAODMCParticle *particle);
  static void FillVarDielectronPair(const AliDielectronPair *pair);
  static void FillVarKFParticle(const AliKFParticle *pair);
  
  static void FillVarVEvent(const AliVEvent *event);
  static void FillVarESDEvent(const AliESDEvent *event);
  static void FillVarAODEvent(const AliAODEvent *event);
  static void FillVarMCEvent(const AliMCEvent *event);
  static void ResetArrayData(Int_t to);
  static void ResetArrayDataMC(Int_t to);

  static Double_t GetPseudoProperTime(const AliDielectronPair *pair);
  
  static AliPIDResponse *fgPIDResponse;        // PID response object
  static AliKFVertex     *fgKFVertex;          // KF vertex
  static AliAODVertex    *fgAODVertex;         // AOD vertex

  
  AliDielectronVarContainer(const AliDielectronVarContainer &c);
  AliDielectronVarContainer &operator=(const AliDielectronVarContainer &c);
  
  ClassDef(AliDielectronVarContainer,1);
};


//Inline functions
inline void AliDielectronVarContainer::Fill(const TObject* object)
{
  //
  // Main function to fill all available variables according to the type of particle
  //

  // Protect
  if (!object) return;

  if      (object->IsA() == AliESDtrack::Class())       FillVarESDtrack(static_cast<const AliESDtrack*>(object));
  else if (object->IsA() == AliAODTrack::Class())       FillVarAODTrack(static_cast<const AliAODTrack*>(object));
  else if (object->IsA() == AliMCParticle::Class())     FillVarMCParticle(static_cast<const AliMCParticle*>(object));
  else if (object->IsA() == AliAODMCParticle::Class())  FillVarAODMCParticle(static_cast<const AliAODMCParticle*>(object));
  else if (object->IsA() == AliDielectronPair::Class()) FillVarDielectronPair(static_cast<const AliDielectronPair*>(object));
  else if (object->IsA() == AliKFParticle::Class())     FillVarKFParticle(static_cast<const AliKFParticle*>(object));
  //else if (object->IsA() == TParticle::Class())         FillVarTParticle(static_cast<const TParticle*>(object));

  // Main function to fill all available variables according to the type of event
  else if (object->IsA() == AliVEvent::Class())         FillVarVEvent(static_cast<const AliVEvent*>(object));
  else if (object->IsA() == AliESDEvent::Class())       FillVarESDEvent(static_cast<const AliESDEvent*>(object));
  else if (object->IsA() == AliAODEvent::Class())       FillVarAODEvent(static_cast<const AliAODEvent*>(object));
  else if (object->IsA() == AliMCEvent::Class())        FillVarMCEvent(static_cast<const AliMCEvent*>(object));
//   else printf(Form("AliDielectronVarContainer::Fill: Type %s is not supported by AliDielectronVarContainer!", object->ClassName())); //TODO: implement without object needed
}



inline void AliDielectronVarContainer::ResetArrayData(Int_t to)
{
  // Protect
  if (to >= AliDielectronVarContainer::kNMaxValues) return;
  // Reset
  for (Int_t i=0; i<to; ++i) fgData[i] = 0.;

  fgData[AliDielectronVarContainer::kTPCchi2Cl] = -1;
}


inline void AliDielectronVarContainer::ResetArrayDataMC(Int_t to)
{
  // Protect
  if (to >= AliDielectronVarContainer::kNMaxValues_MC) return;
  // Reset
  //for (Int_t i=0; i<to; ++i) fgDataMC[i] = 0.;

  //fgDataMC[AliDielectronVarContainer::kPdgCode]            = -1.;
  //fgDataMC[AliDielectronVarContainer::kPdgCodeMother]      = -1.;
  //fgDataMC[AliDielectronVarContainer::kPdgCodeGrandMother] = -1.;
  //fgDataMC[AliDielectronVarContainer::kNumberOfDaughters]  = -999.;
}


inline void AliDielectronVarContainer::FillVarVParticle(const AliVParticle *particle)
{
  //
  // Fill track information available in AliVParticle into array
  //

  // Protect
  if (!particle) return;

  // Reset
  ResetArrayData(AliDielectronVarContainer::kPairMax);

  // Set
  fgData[AliDielectronVarContainer::kPx]        = particle->Px();
  fgData[AliDielectronVarContainer::kPy]        = particle->Py();
  fgData[AliDielectronVarContainer::kPz]        = particle->Pz();
  fgData[AliDielectronVarContainer::kPt]        = particle->Pt();
  fgData[AliDielectronVarContainer::kP]         = particle->P();
  fgData[AliDielectronVarContainer::kXv]        = particle->Xv();
  fgData[AliDielectronVarContainer::kYv]        = particle->Yv();
  fgData[AliDielectronVarContainer::kZv]        = particle->Zv();
  fgData[AliDielectronVarContainer::kOneOverPt] = particle->OneOverPt();
  fgData[AliDielectronVarContainer::kPhi]       = particle->Phi();
  fgData[AliDielectronVarContainer::kTheta]     = particle->Theta();
  fgData[AliDielectronVarContainer::kEta]       = particle->Eta();
  fgData[AliDielectronVarContainer::kY]         = particle->Y();
  fgData[AliDielectronVarContainer::kE]         = particle->E();
  fgData[AliDielectronVarContainer::kM]         = particle->M();
  fgData[AliDielectronVarContainer::kCharge]    = particle->Charge();
  fgData[AliDielectronVarContainer::kPdgCode]   = particle->PdgCode();
}


inline void AliDielectronVarContainer::FillVarVParticleMC(const AliVParticle *particle)
{
  //
  // Fill MC track information available in AliVParticle into array
  //

  // Protect
  if (!particle) return;

  // Get the MC interface if available
  AliDielectronMC *mc = AliDielectronMC::Instance();
  if (!mc->HasMC()) return;

  // Reset
  ResetArrayDataMC(AliDielectronVarContainer::kPairMax_MC);

  // If called for a reco track, get the MC track first
  if (particle->IsA() == AliESDtrack::Class()) particle = mc->GetMCTrack((AliESDtrack*)particle);
  if (particle->IsA() == AliAODTrack::Class()) particle = mc->GetMCTrack((AliAODTrack*)particle);
  if (!particle) return;

  // Set the AliDielectronMC specific info
  if (particle->IsA() == AliMCParticle::Class()) {
    AliMCParticle* mcParticle = (AliMCParticle*)particle;
    if (mc->GetMCTrackMother(mcParticle)) {
      fgDataMC[AliDielectronVarContainer::kPdgCodeMother] = mc->GetMCTrackMother(mcParticle)->PdgCode();
      if (mc->GetMCTrackMother(mc->GetMCTrackMother(mcParticle)))
        fgDataMC[AliDielectronVarContainer::kPdgCodeGrandMother] = mc->GetMCTrackMother(mc->GetMCTrackMother(mcParticle))->PdgCode();;
    }
    fgDataMC[AliDielectronVarContainer::kNumberOfDaughters] = mc->NumberOfDaughters(mcParticle);
  } else if (particle->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle* mcParticle = (AliAODMCParticle*)particle;
    if (mc->GetMCTrackMother(mcParticle)) {
      fgDataMC[AliDielectronVarContainer::kPdgCodeMother] = mc->GetMCTrackMother(mcParticle)->PdgCode();
      if (mc->GetMCTrackMother(mc->GetMCTrackMother(mcParticle)))
        fgDataMC[AliDielectronVarContainer::kPdgCodeGrandMother] = mc->GetMCTrackMother(mc->GetMCTrackMother(mcParticle))->PdgCode();;
    }
    fgDataMC[AliDielectronVarContainer::kNumberOfDaughters] = mc->NumberOfDaughters(mcParticle);
  }


  // Set the common info
  fgData[AliDielectronVarContainer::kIsJpsiPrimary]       = mc->IsJpsiPrimary(particle);
  fgDataMC[AliDielectronVarContainer::kPdgCode]           = particle->PdgCode();
  fgDataMC[AliDielectronVarContainer::kPx_MC]             = particle->Px();
  fgDataMC[AliDielectronVarContainer::kPy_MC]             = particle->Py();
  fgDataMC[AliDielectronVarContainer::kPz_MC]             = particle->Pz();
  fgDataMC[AliDielectronVarContainer::kPt_MC]             = particle->Pt();
  fgDataMC[AliDielectronVarContainer::kP_MC]              = particle->P();
  fgDataMC[AliDielectronVarContainer::kXv_MC]             = particle->Xv();
  fgDataMC[AliDielectronVarContainer::kYv_MC]             = particle->Yv();
  fgDataMC[AliDielectronVarContainer::kZv_MC]             = particle->Zv();
  fgDataMC[AliDielectronVarContainer::kOneOverPt_MC]      = particle->OneOverPt();
  fgDataMC[AliDielectronVarContainer::kPhi_MC]            = particle->Phi();
  fgDataMC[AliDielectronVarContainer::kTheta_MC]          = particle->Theta();
  fgDataMC[AliDielectronVarContainer::kEta_MC]            = particle->Eta();
  fgDataMC[AliDielectronVarContainer::kY_MC]              = particle->Y();
  fgDataMC[AliDielectronVarContainer::kE_MC]              = particle->E();
  fgDataMC[AliDielectronVarContainer::kM_MC]              = particle->M();
  fgDataMC[AliDielectronVarContainer::kCharge_MC]         = particle->Charge();
}


inline void AliDielectronVarContainer::FillVarESDtrack(const AliESDtrack *particle)
{
  //
  // Fill AliESDtrack interface specific information
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle);
  // Fill common MC information if available
  FillVarVParticleMC(particle);

  Double_t tpcNcls=particle->GetTPCNcls();
  Double_t tpcSignalN=particle->GetTPCsignalN();
  fgData[AliDielectronVarContainer::kNclsITS]       = particle->GetNcls(0); // TODO: get rid of the plain numbers
  fgData[AliDielectronVarContainer::kNclsTPC]       = tpcNcls; // TODO: get rid of the plain numbers
  fgData[AliDielectronVarContainer::kNclsTPCiter1]  = particle->GetTPCNclsIter1(); // TODO: get rid of the plain numbers
  fgData[AliDielectronVarContainer::kNFclsTPC]      = particle->GetTPCNclsF();
  fgData[AliDielectronVarContainer::kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
  fgData[AliDielectronVarContainer::kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
  fgData[AliDielectronVarContainer::kTPCsignalN]    = tpcSignalN;
  fgData[AliDielectronVarContainer::kTPCsignalNfrac]= tpcNcls>0?tpcSignalN/tpcNcls:0;
  fgData[AliDielectronVarContainer::kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
  fgData[AliDielectronVarContainer::kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
  fgData[AliDielectronVarContainer::kTRDpidQuality] = particle->GetTRDpidQuality();
  fgData[AliDielectronVarContainer::kTrackStatus]   = (Double_t)particle->GetStatus();
  if (tpcNcls>0) fgData[AliDielectronVarContainer::kTPCchi2Cl] = particle->GetTPCchi2() / tpcNcls;

  //TRD pidProbs
  Double_t pidProbs[AliPID::kSPECIES];
  particle->GetTRDpid(pidProbs);
  fgData[AliDielectronVarContainer::kTRDprobEle]    = pidProbs[AliPID::kElectron];
  fgData[AliDielectronVarContainer::kTRDprobPio]    = pidProbs[AliPID::kPion];

  fgData[AliDielectronVarContainer::kKinkIndex0]    = particle->GetKinkIndex(0);
  
  Float_t impactParXY, impactParZ;
  particle->GetImpactParameters(impactParXY, impactParZ);
  fgData[AliDielectronVarContainer::kImpactParXY]   = impactParXY;
  fgData[AliDielectronVarContainer::kImpactParZ]    = impactParZ;

  fgData[AliDielectronVarContainer::kITSsignal]       = particle->GetITSsignal();
  
  Double_t itsdEdx[4];
  particle->GetITSdEdxSamples(itsdEdx);

  fgData[AliDielectronVarContainer::kITSsignalSSD1]   = itsdEdx[0];
  fgData[AliDielectronVarContainer::kITSsignalSSD2]   = itsdEdx[1];
  fgData[AliDielectronVarContainer::kITSsignalSDD1]   = itsdEdx[2];
  fgData[AliDielectronVarContainer::kITSsignalSDD2]   = itsdEdx[3];
  fgData[AliDielectronVarContainer::kITSclusterMap]   = particle->GetITSClusterMap();
  
  fgData[AliDielectronVarContainer::kTrackLength]     = particle->GetIntegratedLength();

  //dEdx information
  Double_t mom = particle->GetP();
  const AliExternalTrackParam *in=particle->GetInnerParam();
  if (in) mom = in->GetP();
  fgData[AliDielectronVarContainer::kPIn]=mom;
  fgData[AliDielectronVarContainer::kTPCsignal]       = particle->GetTPCsignal();
  fgData[AliDielectronVarContainer::kTOFsignal]       = particle->GetTOFsignal();
  
  Double_t l = particle->GetIntegratedLength();  // cm
  Double_t t = particle->GetTOFsignal();
  Double_t t0 = fgPIDResponse->GetTOFResponse().GetTimeZero(); // ps

  if( (l < 360. || l > 800.) || (t <= 0.) || (t0 >999990.0) ) {
	  fgData[AliDielectronVarContainer::kTOFbeta]=0.0;
  } else {
    t -= t0; // subtract the T0
    l *= 0.01;  // cm ->m
    t *= 1e-12; //ps -> s

    Double_t v = l / t;
    Float_t beta = v / TMath::C();
    fgData[AliDielectronVarContainer::kTOFbeta] = beta;
  }

  // nsigma to Electron band
  // TODO: for the moment we set the bethe bloch parameters manually
  //       this should be changed in future!
  
  fgData[AliDielectronVarContainer::kTPCnSigmaEle]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron)-AliDielectronPID::GetCorrVal();
  fgData[AliDielectronVarContainer::kTPCnSigmaPio]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
  fgData[AliDielectronVarContainer::kTPCnSigmaMuo]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
  fgData[AliDielectronVarContainer::kTPCnSigmaKao]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
  fgData[AliDielectronVarContainer::kTPCnSigmaPro]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);

  fgData[AliDielectronVarContainer::kITSnSigmaEle]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron);
  fgData[AliDielectronVarContainer::kITSnSigmaPio]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kPion);
  fgData[AliDielectronVarContainer::kITSnSigmaMuo]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kMuon);
  fgData[AliDielectronVarContainer::kITSnSigmaKao]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kKaon);
  fgData[AliDielectronVarContainer::kITSnSigmaPro]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kProton);

  fgData[AliDielectronVarContainer::kTOFnSigmaEle]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kElectron);
  fgData[AliDielectronVarContainer::kTOFnSigmaPio]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kPion);
  fgData[AliDielectronVarContainer::kTOFnSigmaMuo]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kMuon);
  fgData[AliDielectronVarContainer::kTOFnSigmaKao]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kKaon);
  fgData[AliDielectronVarContainer::kTOFnSigmaPro]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kProton);
}

inline void AliDielectronVarContainer::FillVarAODTrack(const AliAODTrack *particle)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle);
  // Fill common MC information if available
  FillVarVParticleMC(particle);

  Double_t tpcNcls=particle->GetTPCNcls();
  // Reset AliESDtrack interface specific information
  fgData[AliDielectronVarContainer::kNclsTPC]       = tpcNcls;
  fgData[AliDielectronVarContainer::kNclsTPCiter1]  = tpcNcls; // not really available in AOD
  fgData[AliDielectronVarContainer::kTrackStatus]   = (Double_t)particle->GetStatus();
  
  //TODO: set TRD pidProbs correctly
  
  // Fill AliAODTrack interface information
  //
  Double_t d0z0[2];
  GetDCA(particle, d0z0);
  fgData[AliDielectronVarContainer::kImpactParXY]   = d0z0[0];
  fgData[AliDielectronVarContainer::kImpactParZ]    = d0z0[1];
  fgData[AliDielectronVarContainer::kITSclusterMap] = particle->GetITSClusterMap();
  
  AliAODPid *pid=particle->GetDetPid();
  if (pid){
    Double_t mom =pid->GetTPCmomentum();
    Double_t tpcSignalN=pid->GetTPCsignalN();
    fgData[AliDielectronVarContainer::kTPCsignalN] = tpcSignalN;
    fgData[AliDielectronVarContainer::kTPCsignalN] = tpcNcls>0?tpcSignalN/tpcNcls:0;
    Double_t tpcNsigmaEle=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron);
    Double_t tpcNsigmaPio=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
    Double_t tpcNsigmaMuo=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
    Double_t tpcNsigmaKao=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
    Double_t tpcNsigmaPro=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);
    
    fgData[AliDielectronVarContainer::kPIn]=mom;
    fgData[AliDielectronVarContainer::kTPCsignal]=pid->GetTPCsignal();

    fgData[AliDielectronVarContainer::kTPCnSigmaEle]=tpcNsigmaEle;
    fgData[AliDielectronVarContainer::kTPCnSigmaPio]=tpcNsigmaPio;
    fgData[AliDielectronVarContainer::kTPCnSigmaMuo]=tpcNsigmaMuo;
    fgData[AliDielectronVarContainer::kTPCnSigmaKao]=tpcNsigmaKao;
    fgData[AliDielectronVarContainer::kTPCnSigmaPro]=tpcNsigmaPro;
  }
}


inline void AliDielectronVarContainer::FillVarMCParticle(const AliMCParticle *particle)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle);
  // Fill common MC information if available
  FillVarVParticleMC(particle);


  // Fill AliMCParticle interface specific information

}

inline void AliDielectronVarContainer::FillVarAODMCParticle(const AliAODMCParticle *particle)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle);
  // Fill common MC information if available
  FillVarVParticleMC(particle);
 

  // Fill AliAODMCParticle interface specific information

}

inline void AliDielectronVarContainer::FillVarDielectronPair(const AliDielectronPair *pair)
{
  //
  // Fill pair information available for histogramming into an array
  //
  
  // Fill common AliVParticle interface information
  FillVarVParticle(pair);
  // Reset MC array
  ResetArrayDataMC(AliDielectronVarContainer::kPairMax_MC);

  // Fill AliDielectronPair specific information
  const AliKFParticle &kfPair = pair->GetKFParticle();

  Double_t thetaHE=0;
  Double_t phiHE=0;
  Double_t thetaCS=0;
  Double_t phiCS=0;

  pair->GetThetaPhiCM(thetaHE,phiHE,thetaCS,phiCS);
    
  fgData[AliDielectronVarContainer::kChi2NDF]          = kfPair.GetChi2()/kfPair.GetNDF();
  fgData[AliDielectronVarContainer::kDecayLength]      = kfPair.GetDecayLength();
  fgData[AliDielectronVarContainer::kR]                = kfPair.GetR();
  fgData[AliDielectronVarContainer::kOpeningAngle]     = pair->OpeningAngle();
  fgData[AliDielectronVarContainer::kThetaHE]          = thetaHE;
  fgData[AliDielectronVarContainer::kPhiHE]            = phiHE;
  fgData[AliDielectronVarContainer::kThetaCS]          = thetaCS;
  fgData[AliDielectronVarContainer::kPhiCS]            = phiCS;
  fgData[AliDielectronVarContainer::kLegDist]          = pair->DistanceDaughters();
  fgData[AliDielectronVarContainer::kLegDistXY]        = pair->DistanceDaughtersXY();
  fgData[AliDielectronVarContainer::kDeltaEta]         = pair->DeltaEta();
  fgData[AliDielectronVarContainer::kDeltaPhi]         = pair->DeltaPhi();
  fgData[AliDielectronVarContainer::kMerr]             = kfPair.GetErrMass()>1e-30&&kfPair.GetMass()>1e-30?kfPair.GetErrMass()/kfPair.GetMass():1000000;
  fgData[AliDielectronVarContainer::kPairType]         = pair->GetType();
  fgData[AliDielectronVarContainer::kPseudoProperTime] = GetPseudoProperTime(pair);

  
  AliDielectronMC *mc=AliDielectronMC::Instance();
  if (mc->HasMC()){
    Bool_t samemother =  mc->HaveSameMother(pair);
    fgDataMC[AliDielectronVarContainer::kIsJpsiPrimary] = mc->IsJpsiPrimary(pair);
    fgDataMC[AliDielectronVarContainer::kHaveSameMother] = samemother ;
  }

}


inline void AliDielectronVarContainer::FillVarKFParticle(const AliKFParticle *particle)
{
  //
  // Fill track information available in AliKFParticle into an array
  //

  // Reset data array
  ResetArrayData(AliDielectronVarContainer::kPairMax);
  // Reset MC array
  ResetArrayDataMC(AliDielectronVarContainer::kPairMax_MC);


  fgData[AliDielectronVarContainer::kPx]        = particle->GetPx();
  fgData[AliDielectronVarContainer::kPy]        = particle->GetPy();
  fgData[AliDielectronVarContainer::kPz]        = particle->GetPz();
  fgData[AliDielectronVarContainer::kPt]        = particle->GetPt();
  fgData[AliDielectronVarContainer::kP]         = particle->GetP();
  fgData[AliDielectronVarContainer::kXv]        = particle->GetX();
  fgData[AliDielectronVarContainer::kYv]        = particle->GetY();
  fgData[AliDielectronVarContainer::kZv]        = particle->GetZ();
  fgData[AliDielectronVarContainer::kPhi]       = particle->GetPhi();
  fgData[AliDielectronVarContainer::kEta]       = particle->GetEta();
  fgData[AliDielectronVarContainer::kY]         = ((particle->GetE()*particle->GetE()-particle->GetPx()*particle->GetPx()-particle->GetPy()*particle->GetPy()-particle->GetPz()*particle->GetPz())>0.) ? TLorentzVector(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE()).Rapidity() : -1111.;
  fgData[AliDielectronVarContainer::kE]         = particle->GetE();
  fgData[AliDielectronVarContainer::kM]         = particle->GetMass();
  fgData[AliDielectronVarContainer::kCharge]    = particle->GetQ();
}


inline void AliDielectronVarContainer::FillVarVEvent(const AliVEvent *event)
{
  //
  // Fill event information available for histogramming into an array
  //
  
  // Reset data array
  ResetArrayData(AliDielectronVarContainer::kNMaxValues);
  // Reset MC array
  ResetArrayDataMC(AliDielectronVarContainer::kNMaxValues_MC);

  // set the KF vertex
  if (fgKFVertex) delete fgKFVertex;
  fgKFVertex = 0x0;
  if (!event->GetPrimaryVertex()) return;
  fgKFVertex = new AliKFVertex(*event->GetPrimaryVertex());


  fgData[AliDielectronVarContainer::kXvPrim]       = event->GetPrimaryVertex()->GetX();
  fgData[AliDielectronVarContainer::kYvPrim]       = event->GetPrimaryVertex()->GetY();
  fgData[AliDielectronVarContainer::kZvPrim]       = event->GetPrimaryVertex()->GetZ();
  fgData[AliDielectronVarContainer::kNContrib]     = event->GetPrimaryVertex()->GetNContributors();
  //fgData[AliDielectronVarContainer::kChi2NDF]      = event->GetPrimaryVertex()->GetChi2perNDF(); // This is the pair value!
  fgData[AliDielectronVarContainer::kNTrk]         = event->GetNumberOfTracks();
  fgData[AliDielectronVarContainer::kBzkG]         = event->GetMagneticField();
  fgData[AliDielectronVarContainer::kNacc]         = AliDielectronHelper::GetNacc(event);
  fgData[AliDielectronVarContainer::kNaccTrcklts]  = AliDielectronHelper::GetNaccTrcklts(event);
}


inline void AliDielectronVarContainer::FillVarESDEvent(const AliESDEvent *event)
{
  //
  // Fill event information available for histogramming into an array
  // 
  
  // Fill common AliVEvent interface information
  FillVarVEvent(event);

  Double_t centralityF=-1;
  AliCentrality *esdCentrality = const_cast<AliESDEvent*>(event)->GetCentrality();
  if (esdCentrality) centralityF = esdCentrality->GetCentralityPercentile("V0M");
  
  // Fill AliESDEvent interface specific information
  const AliESDVertex *esdVtx = event->GetPrimaryVertex();
  fgData[AliDielectronVarContainer::kXRes]       = esdVtx->GetXRes();
  fgData[AliDielectronVarContainer::kYRes]       = esdVtx->GetYRes();
  fgData[AliDielectronVarContainer::kZRes]       = esdVtx->GetZRes();
  fgData[AliDielectronVarContainer::kCentrality] = centralityF;
}


inline void AliDielectronVarContainer::FillVarAODEvent(const AliAODEvent *event)
{
  //
  // Fill event information available for histogramming into an array
  //   

  // Fill common AliVEvent interface information
  FillVarVEvent(event);

  // Fill AliAODEvent interface specific information
  // set the AOD vertex
  if (fgAODVertex) delete fgAODVertex;
  fgAODVertex = 0x0;
  if (!event->GetPrimaryVertex()) return;
  fgAODVertex = new AliAODVertex(*event->GetPrimaryVertex());

}


inline void AliDielectronVarContainer::FillVarMCEvent(const AliMCEvent *event)
{ 
  //
  // Fill event information available for histogramming into an array
  //   
        
  // Fill common AliVEvent interface information
  FillVarVEvent(event);

  // Fill AliMCEvent interface specific information
  fgDataMC[AliDielectronVarContainer::kNch] = AliDielectronHelper::GetNch(event, 1.6);
} 


inline Double_t AliDielectronVarContainer::GetPseudoProperTime(const AliDielectronPair *pair) 
{
  //
  // Calculate the pseudo proper time
  //

  if(!pair) return 0.;

  Double_t pt  = pair->Pt();
  Double_t dx  = pair->Xv() - fgData[AliDielectronVarContainer::kXvPrim];
  Double_t dy  = pair->Yv() - fgData[AliDielectronVarContainer::kYvPrim];
  Double_t lxy = ((dx * pair->Px()) + (dy * pair->Py()))/pt;
  Double_t ppt = lxy * (TDatabasePDG::Instance()->GetParticle(443)->Mass())/pt;

  return ppt;
}


inline void AliDielectronVarContainer::SetEvent(AliVEvent * const event)
{
  //
  // Set the event
  //

  FillVarVEvent(event);
}


inline void AliDielectronVarContainer::InitESDpid(Int_t type)
{
  //
  // initialize PID parameters
  // type=0 is simulation
  // type=1 is data

  if (!fgPIDResponse) fgPIDResponse=new AliESDpid((Bool_t)(type==0));
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fgPIDResponse->GetTOFResponse().SetTimeResolution(80.);
  
  // data
  if (type==1){    
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fgPIDResponse->GetTOFResponse().SetTimeResolution(130.);
    fgPIDResponse->GetTPCResponse().SetMip(50.);
  }

  fgPIDResponse->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  fgPIDResponse->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}

inline void AliDielectronVarContainer::InitAODpidUtil(Int_t type)
{
  if (!fgPIDResponse) fgPIDResponse=new AliAODpidUtil;
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fgPIDResponse->GetTOFResponse().SetTimeResolution(80.);
  
  // data
  if (type==1){
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fgPIDResponse->GetTOFResponse().SetTimeResolution(130.);
    fgPIDResponse->GetTPCResponse().SetMip(50.);
  }
  
  fgPIDResponse->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  fgPIDResponse->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}


inline Bool_t AliDielectronVarContainer::GetDCA(const AliAODTrack *track, Double_t d0z0[2])
{
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    return kTRUE;
  }
  
  Double_t covd0z0[3];
  AliAODTrack copy(*track);
  AliAODVertex *vtx =(AliAODVertex*)fgAODVertex;
  Bool_t ok = copy.PropagateToDCA(vtx,fgData[AliDielectronVarContainer::kBzkG],kVeryBig,d0z0,covd0z0);
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}

/*
inline void AliDielectronVarContainer::FillVarTParticle(const TParticle *particle)
{
  //
  // Fill TParticle interface information
  //

  fgData[AliDielectronVarContainer::kPx]        = particle->Px();
  fgData[AliDielectronVarContainer::kPy]        = particle->Py();
  fgData[AliDielectronVarContainer::kPz]        = particle->Pz();
  fgData[AliDielectronVarContainer::kPt]        = particle->Pt();
  fgData[AliDielectronVarContainer::kP]         = particle->P();
  fgData[AliDielectronVarContainer::kXv]        = particle->Vx();
  fgData[AliDielectronVarContainer::kYv]        = particle->Vy();
  fgData[AliDielectronVarContainer::kZv]        = particle->Vz();
  fgData[AliDielectronVarContainer::kOneOverPt] = 1./particle->Pt();
  fgData[AliDielectronVarContainer::kPhi]       = particle->Phi();
  fgData[AliDielectronVarContainer::kTheta]     = particle->Theta();
  fgData[AliDielectronVarContainer::kEta]       = particle->Eta();
  fgData[AliDielectronVarContainer::kY]         = particle->Y();
  fgData[AliDielectronVarContainer::kE]         = particle->Energy();
  fgData[AliDielectronVarContainer::kM]         = particle->GetMass();
  fgData[AliDielectronVarContainer::kCharge]    = particle->GetPDG()->Charge()/3;
}
*/




#endif

