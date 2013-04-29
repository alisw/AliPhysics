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
//#   Markus    Köhler,   GSI / M.Koehler@gsi.de              #
//#   Frederick Kramer,   Uni Ffm / Frederick.Kramer@cern.ch  #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <TNamed.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TKey.h>
#include <TBits.h>

#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliVVertex.h>
#include <AliESDVertex.h>
#include <AliAODVertex.h>
#include <AliEventplane.h>

#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

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

#include "AliAnalysisManager.h"
#include "AliVZEROEPSelectionTask.h"

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
    kITSchi2Cl,              // chi2/cl in the ITS
    kNclsTPC,                // number of clusters assigned in the TPC
    kNclsSTPC,                // number of shared clusters assigned in the TPC
    kNclsSFracTPC,           // fraction of shared clusters assigned in the TPC
    kNclsTPCiter1,           // number of clusters assigned in the TPC after first iteration
    kNFclsTPC,               // number of findable clusters in the TPC
    kNFclsTPCr,              // number of findable clusters(crossed rows) in the TPC with more robust definition
    kNFclsTPCrFrac,          // number of found/findable clusters in the TPC with more robust definition
    kNFclsTPCfCross,         // fraction crossed rows/findable clusters in the TPC, as done in AliESDtrackCuts
    kTPCsignalN,             // number of points used for dEdx
    kTPCsignalNfrac,         // fraction of points used for dEdx / cluster used for tracking
    kTPCchi2Cl,              // chi2/cl in TPC
    kTPCclsDiff,             // TPC cluster difference
    kTPCclsSegments,         // TPC cluster segments
    kTrackStatus,            // track status bits
    kFilterBit,              // AOD filter bits

    kNclsTRD,                // number of clusters assigned in the TRD
    kTRDntracklets,          // number of TRD tracklets used for tracking/PID TODO: correct getter
    kTRDpidQuality,          // number of TRD tracklets used for PID
    kTRDchi2,                // chi2 in TRD
    kTRDprobEle,             // TRD electron pid probability
    kTRDprobPio,             // TRD electron pid probability
    kTRDprob2DEle,           // TRD electron pid probability 2D LQ 
    kTRDprob2DPio,           // TRD electron pid probability 2D LQ
    kTRDphi,                 // Phi angle of the track at the entrance of the TRD
    kTRDpidEffLeg,           // TRD pid efficiency from conversion electrons
    kTRDsignal,              // TRD signal
      
    kImpactParXY,            // Impact parameter in XY plane
    kImpactParZ,             // Impact parameter in Z
    kTrackLength,            // Track length


    kPdgCode,                // PDG code
    kPdgCodeMother, 
    kPdgCodeGrandMother,     // PDG code of the grandmother
    kHasCocktailMother,      // true if particle is added via MC generator cocktail (AliDielectronSignal::kDirect)
    kHasCocktailGrandMother, // true if particle is added via MC generator cocktail (AliDielectronSignal::kDirect)
    kNumberOfDaughters,      // number of daughters
    kHaveSameMother,         // check that particles have the same mother (MC)
    kIsJpsiPrimary,          // check if the particle is primary (MC)
    kNumberOfJPsis,          // number of generated inclusive jpsis per event (MC)
    kNumberOfJPsisPrompt,    // number of generated prompt jpsis per event (MC)
    kNumberOfJPsisNPrompt,   // number of generated non-prompt jpsis per event (MC)

    kITSsignal,              // ITS dE/dx signal
    kITSsignalSSD1,          // SSD1 dE/dx signal
    kITSsignalSSD2,          // SSD2 dE/dx signal
    kITSsignalSDD1,          // SDD1 dE/dx signal
    kITSsignalSDD2,          // SDD2 dE/dx signal
    kITSclusterMap,          // ITS cluster map
    kITSLayerFirstCls,       // No of innermost ITS layer with a cluster of a track
    kITSnSigmaEle,           // number of sigmas to the dE/dx electron line in the ITS
    kITSnSigmaPio,           // number of sigmas to the dE/dx pion line in the ITS
    kITSnSigmaMuo,           // number of sigmas to the dE/dx muon line in the ITS
    kITSnSigmaKao,           // number of sigmas to the dE/dx kaon line in the ITS
    kITSnSigmaPro,           // number of sigmas to the dE/dx proton line in the ITS

    kPIn,                    // momentum at inner wall of TPC (if available), used for PID
    kPOut,                   // momentum at outer wall of TPC, used for TRD studies
    kYsignedIn,              // signed local y at inner wall of TPC
    kTPCsignal,              // TPC dE/dx signal
    
    kTOFsignal,              // TOF signal
    kTOFbeta,                // TOF beta
    kTOFPIDBit,              // TOF PID bit (1:set, 0:TOF not available)a
    kTOFmismProb, 	         // and mismatchPorbability as explain in TOF-twiki
	
    kTPCnSigmaEle,           // number of sigmas to the dE/dx electron line in the TPC
    kTPCnSigmaPio,           // number of sigmas to the dE/dx pion line in the TPC
    kTPCnSigmaMuo,           // number of sigmas to the dE/dx muon line in the TPC
    kTPCnSigmaKao,           // number of sigmas to the dE/dx kaon line in the TPC
    kTPCnSigmaPro,           // number of sigmas to the dE/dx proton line in the TPC
      
    kTOFnSigmaEle,           // number of sigmas to the electron line in the TOF
    kTOFnSigmaPio,           // number of sigmas to the pion line in the TOF
    kTOFnSigmaMuo,           // number of sigmas to the muon line in the TOF
    kTOFnSigmaKao,           // number of sigmas to the kaon line in the TOF
    kTOFnSigmaPro,           // number of sigmas to the proton line in the TOF

    kEMCALnSigmaEle,         // number of sigmas to the proton line in the TOF
    kEMCALEoverP,            // E over P from EMCAL
    kEMCALE,                 // E from EMCAL
    kEMCALNCells,            // NCells from EMCAL
    kEMCALM02,               // M02 showershape parameter
    kEMCALM20,               // M20 showershape parameter
    kEMCALDispersion,        // Dispersion paramter
    
    kV0Index0,               // v0 index 0
    kKinkIndex0,             // kink index 0
      
    kParticleMax,             //
    // TODO: kRNClusters ??
    // AliDielectronPair specific variables
    kChi2NDF = kParticleMax, // Chi^2/NDF
    kDecayLength,            // decay length
    kR,                      // distance to the origin
    kOpeningAngle,           // opening angle
    kCosPointingAngle,       // cosine of the pointing angle
    kArmAlpha,               // Armenteros-Podolanski alpha
    kArmPt,                  // Armenteros-Podolanski pt
    // helicity picture: Z-axis is considered the direction of the mother's 3-momentum vector
    kThetaHE,                // theta in mother's rest frame in the helicity picture 
    kPhiHE,                  // phi in mother's rest frame in the helicity picture
    kThetaSqHE,              // squared value of kThetaHE
    kCos2PhiHE,              // Cosine of 2*phi in mother's rest frame in the helicity picture
    kCosTilPhiHE,            // Shifted phi depending on kThetaHE
    // Collins-Soper picture: Z-axis is considered the direction of the vectorial difference between 
    // the 3-mom vectors of target and projectile beams
    kThetaCS,                // theta in mother's rest frame in Collins-Soper picture
    kPhiCS,                  // phi in mother's rest frame in Collins-Soper picture
    kThetaSqCS,              // squared value of kThetaCS
    kPsiPair,                // phi in mother's rest frame in Collins-Soper picture
	kPhivPair,               // angle between ee plane and the magnetic field (can be useful for conversion rejection)
	kPairPlanev0rpH2Angle,   // angle between ee plane and VZERO-C reaction plane           
	kPairPlaneMagAngle,      // angle between ee plane and strong magnetic field  
	kPairPlaneAngle,         // angle between ee plane and strong magnetic field
	kRotPairx,               //ee plane vector
	kRotPairy,               //ee plane vector
	kRotPairz,               //ee plane vector
	kCos2PhiCS,              // Cosine of 2*phi in mother's rest frame in the Collins-Soper picture
    kCosTilPhiCS,            // Shifted phi depending on kThetaCS
    kDeltaPhiV0ArpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-A
    kDeltaPhiV0CrpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-C
    kDeltaPhiV0ACrpH2,       // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-A + V0-C
    kV0ArpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-A
    kV0CrpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-C
    kV0ACrpH2FlowV2,         // v2 coefficient with respect to the 2nd order reaction plane from V0-A + V0-C
    kDeltaPhiv0ArpH2,          // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-A (EPtask)
    kDeltaPhiv0CrpH2,          // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-C
    kDeltaPhiv0ACrpH2,         // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-AC
    kv0ArpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-A (EPtask)
    kv0CrpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-C
    kv0ACrpH2FlowV2,         // v2 coefficient with respect to the 2nd order reaction plane from V0-A + V0-C

    kLegDist,                // distance of the legs
    kLegDistXY,              // distance of the legs in XY
    kDeltaEta,         // Absolute value of Delta Eta for the legs
    kDeltaPhi,           // Absolute value of Delta Phi for the legs
    kMerr,                   // error of mass calculation
    kDCA,                    // distance of closest approach TODO: not implemented yet
    kPairType,               // type of the pair, like like sign ++ unlikesign ...
    kPseudoProperTime,       // pseudo proper time
    kPseudoProperTimeErr,    // pseudo proper time error
    kPseudoProperTimeResolution,     // resolution for pseudo proper decay time (reconstructed - MC truth)
    kPseudoProperTimePull,   // normalizd resolution for pseudo proper time = (reco - MC truth)/dReco
    kTRDpidEffPair,          // TRD pid efficieny from conversion electrons
    kMomAsymDau1,            // momentum fraction of daughter1
    kMomAsymDau2,            // momentum fraction of daughter2
    kPairMax,                 //
  // Event specific variables
    kXvPrim=kPairMax,        // prim vertex
    kYvPrim,                 // prim vertex
    kZvPrim,                 // prim vertex
    kXRes,                   // primary vertex x-resolution
    kYRes,                   // primary vertex y-resolution
    kZRes,                   // primary vertex z-resolution
    kPhiMaxPt,               // phi angle of the track with maximum pt
    kMaxPt,                  // track with maximum pt

    //// v0 reaction plane quantities from AliEPSelectionTaks, angles interval [-pi,+pi]
    kv0ArpH2,                // VZERO-A reaction plane of the Q vector for 2nd harmonic
    kv0CrpH2,                //         reaction plane
    kv0ACrpH2,               // VZERO-AC reaction plane of the Q vector for 2nd harmonic
    kv0AxH2,                 // VZERO-A x-component of the Q vector for 2nd harmonic
    kv0AyH2,                 // VZERO-A y-component of the Q vector for 2nd harmonic
    kv0CxH2,                 // VZERO-C x-component of the Q vector for 2nd harmonic
    kv0CyH2,                 // VZERO-C y-component of the Q vector for 2nd harmonic
    kv0ACxH2,                // VZERO-AC x-component of the Q vector for 2nd harmonic
    kv0ACyH2,                // VZERO-AC y-component of the Q vector for 2nd harmonic
    kv0AmagH2,               // VZERO-A the Q vectors magnitude for 2nd harmonic
    kv0CmagH2,               // VZERO-A the Q vectors magnitude for 2nd harmonic
    kv0ACmagH2,              // VZERO-A the Q vectors magnitude for 2nd harmonic
    kv0A0rpH2,                 // VZERO-A 1st  ring reaction plane of the Q vector for 2nd harmonic
    kv0A3rpH2,                 // VZERO-A last ring reaction plane of the Q vector for 2nd harmonic
    kv0C0rpH2,                 // VZERO-C 1st  ring reaction plane of the Q vector for 2nd harmonic
    kv0C3rpH2,                 // VZERO-C last ring reaction plane of the Q vector for 2nd harmonic
    kv0ATPCDiffH2,             // V0A-TPC reaction plane difference for 2nd harmonic
    kv0CTPCDiffH2,             // V0C-TPC reaction plane difference for 2nd harmonic
    kv0Av0CDiffH2,             // V0A-V0C reaction plane difference for 2nd harmonic
    kv0Av0C0DiffH2,             // V0A-ring 0 ofV0C reaction plane difference for 2nd harmonic
    kv0Av0C3DiffH2,             // V0A-ring 3 ofV0C reaction plane difference for 2nd harmonic
    kv0Cv0A0DiffH2,             // V0C-ring 0 ofV0A reaction plane difference for 2nd harmonic
    kv0Cv0A3DiffH2,             // V0C-ring 3 ofV0A reaction plane difference for 2nd harmonic
    kv0A0v0A3DiffH2,             // V0C-ring 0 ofV0A reaction plane difference for 2nd harmonic
    kv0C0v0C3DiffH2,             // V0C-ring 0 ofV0A reaction plane difference for 2nd harmonic

    kMultV0A,                // VZERO multiplicity and ADC amplitudes
    kMultV0C,
    kMultV0,
    kAdcV0A,
    kAdcV0C,
    kAdcV0,
    kVZEROchMult,
    // VZERO reaction plane quantities
    kV0AxH2=kVZEROchMult+64,   // VZERO-A x-component of the Q vector for 2nd harmonic
    kV0AyH2,                   // VZERO-A y-component of the Q vector for 2nd harmonic
    kV0ArpH2,                  // VZERO-A reaction plane of the Q vector for 2nd harmonic
    kV0CxH2,                   // VZERO-C x-component of the Q vector for 2nd harmonic
    kV0CyH2,                   //         y-component
    kV0CrpH2,                  //         reaction plane
    kV0ACxH2,                  // VZERO-AC x-component of the Q vector for 2nd harmonic
    kV0ACyH2,                  // VZERO-AC y-component of the Q vector for 2nd harmonic
    kV0ACrpH2,                 // VZERO-AC reaction plane of the Q vector for 2nd harmonic
    kV0ArpResH2,               // 2nd harmonic reaction plane resolution for V0A
    kV0CrpResH2,               //                               V0C
    kV0ACrpResH2,              //                             V0A+V0C
    kV0XaXcH2,                 // Correlation quantities to check V0 reaction plane quality
    kV0XaYaH2,
    kV0XaYcH2,
    kV0YaXcH2,
    kV0YaYcH2,
    kV0XcYcH2,
    kV0ATPCDiffH2,             // V0A-TPC reaction plane difference for 2nd harmonic
    kV0CTPCDiffH2,             // V0C-TPC reaction plane difference for 2nd harmonic
    kV0AV0CDiffH2,             // V0A-V0C reaction plane difference for 2nd harmonic
    // TPC reaction plane quantities, angle interval [0,+pi]
    kTPCxH2,                  // TPC x-component of the Q vector for 2nd harmonic (corrected)
    kTPCyH2,                  // TPC y-component of the Q vector for 2nd harmonic (corrected)
    kTPCmagH2,                // TPC reaction plane the Q vectors magnitude for 2nd harmonic (corrected)
    kTPCrpH2,                 // TPC reaction plane angle of the Q vector for 2nd harmonic (corrected)
    kTPCsub1xH2,              // TPC x-component of the Q vector for 2nd harmonic (corrected, sub event 1) 
    kTPCsub1yH2,              // TPC y-component of the Q vector for 2nd harmonic (corrected, sub event 1)
    kTPCsub1rpH2,             // TPC reaction plane of the Q vector for 2nd harmonic (corrected, sub event 1)
    kTPCsub2xH2,              // TPC x-component of the Q vector for 2nd harmonic (corrected, sub event 2)
    kTPCsub2yH2,              // TPC y-component of the Q vector for 2nd harmonic (corrected, sub event 2)
    kTPCsub2rpH2,             // TPC reaction plane of the Q vector for 2nd harmonic (corrected, sub event 2)
    kTPCsub12DiffH2,          // TPC reaction plane difference of sub event 1,2 for 2nd harmonic
    kTPCsub12DiffH2Sin,       // TPC reaction plane difference of sub event 1,2 for 2nd harmonic, sinus term

    kTPCxH2uc,                  // TPC x-component of the Q vector for 2nd harmonic (uncorrected)
    kTPCyH2uc,                  // TPC y-component of the Q vector for 2nd harmonic (uncorrected)
    kTPCmagH2uc,                // TPC reaction plane the Q vectors magnitude for 2nd harmonic (uncorrected)
    kTPCrpH2uc,                 // TPC reaction plane angle of the Q vector for 2nd harmonic (uncorrected)
    kTPCsub1xH2uc,              // TPC x-component of the Q vector for 2nd harmonic (uncorrected, sub event 1) 
    kTPCsub1yH2uc,              // TPC y-component of the Q vector for 2nd harmonic (uncorrected, sub event 1)
    kTPCsub1rpH2uc,             // TPC reaction plane of the Q vector for 2nd harmonic (uncorrected, sub event 1)
    kTPCsub2xH2uc,              // TPC x-component of the Q vector for 2nd harmonic (uncorrected, sub event 2)
    kTPCsub2yH2uc,              // TPC y-component of the Q vector for 2nd harmonic (uncorrected, sub event 2)
    kTPCsub2rpH2uc,             // TPC reaction plane of the Q vector for 2nd harmonic (uncorrected, sub event 2)
    kTPCsub12DiffH2uc,          // TPC reaction plane difference of sub event 1,2 for 2nd harmonic (uncorrected)

    kNTrk,                   // number of tracks (or tracklets) TODO: ambiguous
    kTracks,                 // ESD tracks TODO: ambiguous
    kNVtxContrib,             // number of primary vertex contibutors
    kNVtxContribTPC,         // number of TPC vertex contibutors
    kNacc,                   // Number of accepted tracks
    kMatchEffITSTPC,         // ruff estimate on the ITS-TPC matching efficiceny
    kNaccTrcklts,            // number of accepted SPD tracklets in |eta|<1.6        
    kNaccTrcklts0916,        // number of accepted SPD tracklets in 0.9<|eta|<1.6
    
    kNaccTrckltsEsd05,       // number of accepted SPD tracklets in |eta|<0.5 (AliESDEvent::EstimateMultiplicity())
    kNaccTrckltsEsd10,       // number of accepted SPD tracklets in |eta|<1.0 (AliESDEvent::EstimateMultiplicity())
    kNaccTrckltsEsd16,       // number of accepted SPD tracklets in |eta|<1.6 (AliESDEvent::EstimateMultiplicity())
    kNaccTrckltsEsd05Corr,   //
    kNaccTrckltsEsd10Corr,   //
    kNaccTrckltsEsd16Corr,   //
    kNaccItsTpcEsd05,        // ITS-TPC tracks + ITS SA complementary tracks + tracklets from unassigned tracklets in |eta|<0.5 (AliESDEvent::EstimateMultiplicity())
    kNaccItsTpcEsd10,        // ITS-TPC tracks + ITS SA complementary tracks + tracklets from unassigned tracklets in |eta|<1.0 (AliESDEvent::EstimateMultiplicity())
    kNaccItsTpcEsd16,        // ITS-TPC tracks + ITS SA complementary tracks + tracklets from unassigned tracklets in |eta|<1.6 (AliESDEvent::EstimateMultiplicity())
    kNaccItsTpcEsd05Corr,        // 
    kNaccItsTpcEsd10Corr,        // 
    kNaccItsTpcEsd16Corr,        // 
    
    kNaccItsPureEsd05,       // ITS SA tracks + tracklets from unassigned tracklets in |eta|<0.5 (AliESDEvent::EstimateMultiplicity())
    kNaccItsPureEsd10,       // ITS SA tracks + tracklets from unassigned tracklets in |eta|<1.0 (AliESDEvent::EstimateMultiplicity())
    kNaccItsPureEsd16,       // ITS SA tracks + tracklets from unassigned tracklets in |eta|<1.6 (AliESDEvent::EstimateMultiplicity())
    kNaccItsPureEsd05Corr,   // 
    kNaccItsPureEsd10Corr,   // 
    kNaccItsPureEsd16Corr,   // 
    kRefMult,                // reference multiplicity (only in AODs) should be Ntrk w/o double counts
    kRefMultTPConly,         // TPC only Reference Multiplicty (AliESDtrackCuts::GetReferenceMultiplicity(&esd, kTRUE))
    
    kNch,                    // MC true number of charged particles in |eta|<1.6
    kNch05,                  // MC true number of charged particles in |eta|<0.5
    kNch10,                  // MC true number of charged particles in |eta|<1.0

    kCentrality,             // event centrality fraction
    kCentralitySPD,          // centrality using SPD
    kNevents,                // event counter
    kRunNumber,              // run number
    kMixingBin,
    kNMaxValues              //
    // TODO: (for A+A) ZDCEnergy, impact parameter, Iflag??
  };
  

  AliDielectronVarManager();
  AliDielectronVarManager(const char* name, const char* title);
  virtual ~AliDielectronVarManager();
  static void Fill(const TObject* particle, Double_t * const values);
  static void FillVarMCParticle2(const AliVParticle *p1, const AliVParticle *p2, Double_t * const values);
  static void FillVarVParticle(const AliVParticle *particle,         Double_t * const values);

  static void InitESDpid(Int_t type=0);
  static void InitAODpidUtil(Int_t type=0);
  static void InitEstimatorAvg(const Char_t* filename);
  static void InitTRDpidEffHistograms(const Char_t* filename);
  static void SetVZEROCalibrationFile(const Char_t* filename) {fgVZEROCalibrationFile = filename;}
  
  static void SetVZERORecenteringFile(const Char_t* filename) {fgVZERORecenteringFile = filename;}
  static void SetPIDResponse(AliPIDResponse *pidResponse) {fgPIDResponse=pidResponse;}
  static AliPIDResponse* GetPIDResponse() { return fgPIDResponse; }
  static void SetEvent(AliVEvent * const ev);
  static void SetEventData(const Double_t data[AliDielectronVarManager::kNMaxValues]);
  static Bool_t GetDCA(const AliAODTrack *track, Double_t d0z0[2]);
  static void SetTPCEventPlane(AliEventplane *const evplane);
  static void GetVzeroRP(const AliVEvent* event, Double_t* qvec, Int_t sideOption);      // 0- V0A; 1- V0C; 2- V0A+V0C
  static AliAODVertex* GetVertex(const AliAODEvent *event, AliAODVertex::AODVtx_t vtype);

  static TProfile* GetEstimatorHistogram(Int_t period, Int_t type) {return fgMultEstimatorAvg[period][type];}
  static Double_t GetTRDpidEfficiency(Int_t runNo, Double_t centrality, Double_t eta, Double_t trdPhi, Double_t pout, Double_t& effErr);

  static const AliKFVertex* GetKFVertex() {return fgKFVertex;}
  
  static const char* GetValueName(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i][0]:""; }
  static const char* GetValueLabel(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i][1]:""; }
  static const char* GetValueUnit(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i][2]:""; }
  static UInt_t GetValueType(const char* valname);
  static const Double_t* GetData() {return fgData;}
  static AliVEvent* GetCurrentEvent() {return fgEvent;}

  static Double_t GetValue(ValueTypes var) {return fgData[var];}
  static void SetValue(ValueTypes var, Double_t val) { fgData[var]=val; }
  
private:

  static const char* fgkParticleNames[kNMaxValues][3];  //variable names


  static void FillVarESDtrack(const AliESDtrack *particle,           Double_t * const values);
  static void FillVarAODTrack(const AliAODTrack *particle,           Double_t * const values);
  static void FillVarMCParticle(const AliMCParticle *particle,       Double_t * const values);
  static void FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values);
  static void FillVarDielectronPair(const AliDielectronPair *pair,   Double_t * const values);
  static void FillVarKFParticle(const AliKFParticle *pair,   Double_t * const values);
  
  static void FillVarVEvent(const AliVEvent *event,                  Double_t * const values);
  static void FillVarESDEvent(const AliESDEvent *event,              Double_t * const values);
  static void FillVarAODEvent(const AliAODEvent *event,              Double_t * const values);
  static void FillVarMCEvent(const AliMCEvent *event,                Double_t * const values);
  static void FillVarTPCEventPlane(const AliEventplane *evplane,     Double_t * const values);

  static void InitVZEROCalibrationHistograms(Int_t runNo);
  static void InitVZERORecenteringHistograms(Int_t runNo);
  
  static AliPIDResponse  *fgPIDResponse;        // PID response object
  static AliVEvent       *fgEvent;              // current event pointer
  static AliEventplane   *fgTPCEventPlane;      // current event tpc plane pointer
  static AliKFVertex     *fgKFVertex;           // kf vertex
  static TProfile        *fgMultEstimatorAvg[4][9];  // multiplicity estimator averages (4 periods x 9 estimators)
  static Double_t         fgTRDpidEffCentRanges[10][4];   // centrality ranges for the TRD pid efficiency histograms
  static TH3D            *fgTRDpidEff[10][4];   // TRD pid efficiencies from conversion electrons
  static TString          fgVZEROCalibrationFile;  // file with VZERO channel-by-channel calibrations
  static TString          fgVZERORecenteringFile;  // file with VZERO Q-vector averages needed for event plane recentering
  static TProfile2D      *fgVZEROCalib[64];           // 1 histogram per VZERO channel
  static TProfile2D      *fgVZERORecentering[2][2];   // 2 VZERO sides x 2 Q-vector components
  static Int_t            fgCurrentRun;               // current run number
  
  static Double_t fgData[kNMaxValues];        //! data
  
  AliDielectronVarManager(const AliDielectronVarManager &c);
  AliDielectronVarManager &operator=(const AliDielectronVarManager &c);
  
  ClassDef(AliDielectronVarManager,1);
};


//Inline functions
inline void AliDielectronVarManager::Fill(const TObject* object, Double_t * const values)
{
  //
  // Main function to fill all available variables according to the type of particle
  //
  if (!object) return;
  if      (object->IsA() == AliESDtrack::Class())       FillVarESDtrack(static_cast<const AliESDtrack*>(object), values);
  else if (object->IsA() == AliAODTrack::Class())       FillVarAODTrack(static_cast<const AliAODTrack*>(object), values);
  else if (object->IsA() == AliMCParticle::Class())     FillVarMCParticle(static_cast<const AliMCParticle*>(object), values);
  else if (object->IsA() == AliAODMCParticle::Class())  FillVarAODMCParticle(static_cast<const AliAODMCParticle*>(object), values);
  else if (object->IsA() == AliDielectronPair::Class()) FillVarDielectronPair(static_cast<const AliDielectronPair*>(object), values);
  else if (object->IsA() == AliKFParticle::Class())     FillVarKFParticle(static_cast<const AliKFParticle*>(object),values);
  // Main function to fill all available variables according to the type of event
  
  else if (object->IsA() == AliVEvent::Class())         FillVarVEvent(static_cast<const AliVEvent*>(object), values);
  else if (object->IsA() == AliESDEvent::Class())       FillVarESDEvent(static_cast<const AliESDEvent*>(object), values);
  else if (object->IsA() == AliAODEvent::Class())       FillVarAODEvent(static_cast<const AliAODEvent*>(object), values);
  else if (object->IsA() == AliMCEvent::Class())        FillVarMCEvent(static_cast<const AliMCEvent*>(object), values);
  else if (object->IsA() == AliEventplane::Class())     FillVarTPCEventPlane(static_cast<const AliEventplane*>(object), values);
//   else printf(Form("AliDielectronVarManager::Fill: Type %s is not supported by AliDielectronVarManager!", object->ClassName())); //TODO: implement without object needed
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

  values[AliDielectronVarManager::kOneOverPt] = (particle->Pt()>1.0e-3 ? particle->OneOverPt() : 0.0);
  values[AliDielectronVarManager::kPhi]       = particle->Phi();
  values[AliDielectronVarManager::kTheta]     = particle->Theta();
  values[AliDielectronVarManager::kEta]       = particle->Eta();
  values[AliDielectronVarManager::kY]         = particle->Y();
  
  values[AliDielectronVarManager::kE]         = particle->E();
  values[AliDielectronVarManager::kM]         = particle->M();
  values[AliDielectronVarManager::kCharge]    = particle->Charge();
  
  values[AliDielectronVarManager::kPdgCode]   = particle->PdgCode();
    
//   if ( fgEvent ) AliDielectronVarManager::Fill(fgEvent, values);
  for (Int_t i=AliDielectronVarManager::kPairMax; i<AliDielectronVarManager::kNMaxValues; ++i)
    values[i]=fgData[i];
}

inline void AliDielectronVarManager::FillVarESDtrack(const AliESDtrack *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  AliESDtrack *esdTrack=0x0;
  Double_t origdEdx=particle->GetTPCsignal();
  
  // apply ETa correction, remove once this is in the tender
  esdTrack=const_cast<AliESDtrack*>(particle);
  if (!esdTrack) return;
  esdTrack->SetTPCsignal(origdEdx/AliDielectronPID::GetEtaCorr(esdTrack)/AliDielectronPID::GetCorrValdEdx(),esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());
  
  Double_t pidProbs[AliPID::kSPECIES];
  // Fill AliESDtrack interface specific information
  Double_t tpcNcls=particle->GetTPCNcls();
  Double_t tpcNclsS = particle->GetTPCnclsS(); 
  Double_t itsNcls=particle->GetNcls(0);
  Double_t tpcSignalN=particle->GetTPCsignalN();
  Double_t tpcClusFindable=particle->GetTPCNclsF();
  values[AliDielectronVarManager::kNclsITS]       = itsNcls; // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNclsTPC]       = tpcNcls; // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNclsSTPC]      = tpcNclsS;
  values[AliDielectronVarManager::kNclsSFracTPC]  = tpcNcls>0?tpcNclsS/tpcNcls:0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = particle->GetTPCNclsIter1(); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kNFclsTPC]       = tpcClusFindable;
  values[AliDielectronVarManager::kNFclsTPCr]      = particle->GetTPCClusterInfo(2,1);
  values[AliDielectronVarManager::kNFclsTPCrFrac]  = particle->GetTPCClusterInfo(2);
  values[AliDielectronVarManager::kNFclsTPCfCross]= (tpcClusFindable>0)?(particle->GetTPCClusterInfo(2,1)/tpcClusFindable):0;
  values[AliDielectronVarManager::kTPCsignalN]    = tpcSignalN;
  values[AliDielectronVarManager::kTPCsignalNfrac]= tpcNcls>0?tpcSignalN/tpcNcls:0;
  values[AliDielectronVarManager::kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
  values[AliDielectronVarManager::kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
  values[AliDielectronVarManager::kTRDpidQuality] = particle->GetTRDntrackletsPID();
  values[AliDielectronVarManager::kTRDchi2]       = particle->GetTRDchi2();
  values[AliDielectronVarManager::kTRDsignal]     = particle->GetTRDsignal();
  values[AliDielectronVarManager::kTPCclsDiff]    = tpcSignalN-tpcNcls;
  values[AliDielectronVarManager::kTPCclsSegments] = 0.0;
  const UChar_t threshold = 5;
  TBits tpcClusterMap = particle->GetTPCClusterMap();
  UChar_t n=0; UChar_t j=0;
  for(UChar_t i=0; i<8; ++i) {
    n=0;
    for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliDielectronVarManager::kTPCclsSegments] += 1.0;
  }
  values[AliDielectronVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  values[AliDielectronVarManager::kFilterBit]     = 0;
  
  values[AliDielectronVarManager::kTPCchi2Cl] = -1;
  if (tpcNcls>0) values[AliDielectronVarManager::kTPCchi2Cl] = particle->GetTPCchi2() / tpcNcls;
  values[AliDielectronVarManager::kITSchi2Cl] = -1;
  if (itsNcls>0) values[AliDielectronVarManager::kITSchi2Cl] = particle->GetITSchi2() / itsNcls;
  //TRD pidProbs
  particle->GetTRDpid(pidProbs);
  values[AliDielectronVarManager::kTRDprobEle]    = pidProbs[AliPID::kElectron];
  values[AliDielectronVarManager::kTRDprobPio]    = pidProbs[AliPID::kPion];

  values[AliDielectronVarManager::kV0Index0]      = particle->GetV0Index(0);
  values[AliDielectronVarManager::kKinkIndex0]    = particle->GetKinkIndex(0);
  
  Float_t impactParXY, impactParZ;
  particle->GetImpactParameters(impactParXY, impactParZ);
  values[AliDielectronVarManager::kImpactParXY]   = impactParXY;
  values[AliDielectronVarManager::kImpactParZ]    = impactParZ;


  values[AliDielectronVarManager::kPdgCode]=-1;
  values[AliDielectronVarManager::kPdgCodeMother]=-1;
  values[AliDielectronVarManager::kPdgCodeGrandMother]=-1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;
  
  values[AliDielectronVarManager::kNumberOfDaughters]=-999;
  
  AliDielectronMC *mc=AliDielectronMC::Instance();
  
  if (mc->HasMC()){
    if (mc->GetMCTrack(particle)) {
      values[AliDielectronVarManager::kPdgCode]=mc->GetMCTrack(particle)->PdgCode();
      Int_t trkLbl = mc->GetMCTrack(particle)->GetLabel();
      values[AliDielectronVarManager::kHasCocktailMother]=mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
    }
    AliMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
    if (motherMC){
      values[AliDielectronVarManager::kPdgCodeMother]=motherMC->PdgCode();
      Int_t motherLbl = motherMC->GetLabel();
      values[AliDielectronVarManager::kHasCocktailGrandMother]=mc->CheckParticleSource(motherLbl, AliDielectronSignalMC::kDirect);
      
      motherMC=mc->GetMCTrackMother(motherMC);  //grand motherMC
      if (motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=motherMC->PdgCode();;
    }
      
    values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
  } //if(mc->HasMC())
  


  values[AliDielectronVarManager::kITSsignal]   =   particle->GetITSsignal();
  
  Double_t itsdEdx[4];
  particle->GetITSdEdxSamples(itsdEdx);

  values[AliDielectronVarManager::kITSsignalSSD1]   =   itsdEdx[0];
  values[AliDielectronVarManager::kITSsignalSSD2]   =   itsdEdx[1];
  values[AliDielectronVarManager::kITSsignalSDD1]   =   itsdEdx[2];
  values[AliDielectronVarManager::kITSsignalSDD2]   =   itsdEdx[3];
  values[AliDielectronVarManager::kITSclusterMap]   =   particle->GetITSClusterMap();
  values[AliDielectronVarManager::kITSLayerFirstCls] = -1.;

  for (Int_t iC=0; iC<6; iC++) {
    if (((particle->GetITSClusterMap()) & (1<<(iC))) > 0) {
      values[AliDielectronVarManager::kITSLayerFirstCls] = iC;
      break;
    }
  }

  
  values[AliDielectronVarManager::kTrackLength]   = particle->GetIntegratedLength();
  //dEdx information
  Double_t mom = particle->GetP();
  const AliExternalTrackParam *in=particle->GetInnerParam();
  Double_t ysignedIn=-100;
  if (in) {
    mom = in->GetP();
    ysignedIn=particle->Charge()*in->GetY();
  }
  values[AliDielectronVarManager::kPIn]=mom;
  values[AliDielectronVarManager::kYsignedIn]=ysignedIn;
  const AliExternalTrackParam *out=particle->GetOuterParam();
  if(out) values[AliDielectronVarManager::kPOut] = out->GetP();
  else values[AliDielectronVarManager::kPOut] = mom;
  if(out && fgEvent) {
    Double_t localCoord[3]={0.0};
    Bool_t localCoordGood = out->GetXYZAt(298.0, ((AliESDEvent*)fgEvent)->GetMagneticField(), localCoord);
    values[AliDielectronVarManager::kTRDphi] = (localCoordGood && TMath::Abs(localCoord[0])>1.0e-6 && TMath::Abs(localCoord[1])>1.0e-6 ? TMath::ATan2(localCoord[1], localCoord[0]) : -999.);
  }
  if(mc->HasMC() && fgTRDpidEff[0][0]) {
    Int_t runNo = (fgEvent ? fgEvent->GetRunNumber() : -1);
    Float_t centrality=-1.0;
    AliCentrality *esdCentrality = (fgEvent ? fgEvent->GetCentrality() : 0x0);
    if(esdCentrality) centrality = esdCentrality->GetCentralityPercentile("V0M");
    Double_t effErr=0.0;
    values[kTRDpidEffLeg] = GetTRDpidEfficiency(runNo, centrality, values[AliDielectronVarManager::kEta], 
 						values[AliDielectronVarManager::kTRDphi], 
 						values[AliDielectronVarManager::kPOut], effErr);
  }
  values[AliDielectronVarManager::kTPCsignal]=particle->GetTPCsignal();
  
  values[AliDielectronVarManager::kTOFsignal]=particle->GetTOFsignal();
  
  Double_t l = particle->GetIntegratedLength();  // cm
  Double_t t = particle->GetTOFsignal();
  Double_t t0 = fgPIDResponse->GetTOFResponse().GetTimeZero(); // ps

  if( (l < 360. || l > 800.) || (t <= 0.) || (t0 >999990.0) ) {
	values[AliDielectronVarManager::kTOFbeta]=0.0;
  }
  else {
	t -= t0; // subtract the T0
	l *= 0.01;  // cm ->m
	t *= 1e-12; //ps -> s
    
	Double_t v = l / t;
	Float_t beta = v / TMath::C();
	values[AliDielectronVarManager::kTOFbeta]=beta;
  }
  values[AliDielectronVarManager::kTOFPIDBit]=(particle->GetStatus()&AliESDtrack::kTOFpid? 1: 0);

  values[AliDielectronVarManager::kTOFmismProb] = fgPIDResponse->GetTOFMismatchProbability(particle);
  
  // nsigma to Electron band
  // TODO: for the moment we set the bethe bloch parameters manually
  //       this should be changed in future!
  values[AliDielectronVarManager::kTPCnSigmaEle]=(fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron)-AliDielectronPID::GetCorrVal()-AliDielectronPID::GetCntrdCorr(particle)) / AliDielectronPID::GetWdthCorr(particle);
  values[AliDielectronVarManager::kTPCnSigmaPio]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
  values[AliDielectronVarManager::kTPCnSigmaMuo]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
  values[AliDielectronVarManager::kTPCnSigmaKao]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
  values[AliDielectronVarManager::kTPCnSigmaPro]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);

  values[AliDielectronVarManager::kITSnSigmaEle]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kITSnSigmaPio]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kPion);
  values[AliDielectronVarManager::kITSnSigmaMuo]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kMuon);
  values[AliDielectronVarManager::kITSnSigmaKao]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kKaon);
  values[AliDielectronVarManager::kITSnSigmaPro]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kProton);

  values[AliDielectronVarManager::kTOFnSigmaEle]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kTOFnSigmaPio]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kPion);
  values[AliDielectronVarManager::kTOFnSigmaMuo]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kMuon);
  values[AliDielectronVarManager::kTOFnSigmaKao]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kKaon);
  values[AliDielectronVarManager::kTOFnSigmaPro]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kProton);

  //EMCAL PID information
  Double_t eop=0;
  Double_t showershape[4]={0.,0.,0.,0.};
//   values[AliDielectronVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron,eop,showershape);
  values[AliDielectronVarManager::kEMCALEoverP]     = eop;
  values[AliDielectronVarManager::kEMCALE]          = eop*values[AliDielectronVarManager::kP];
  values[AliDielectronVarManager::kEMCALNCells]     = showershape[0];
  values[AliDielectronVarManager::kEMCALM02]        = showershape[1];
  values[AliDielectronVarManager::kEMCALM20]        = showershape[2];
  values[AliDielectronVarManager::kEMCALDispersion] = showershape[3];
  
  
  //restore TPC signal if it was changed
  if (esdTrack) esdTrack->SetTPCsignal(origdEdx,esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());
}

inline void AliDielectronVarManager::FillVarAODTrack(const AliAODTrack *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  Double_t tpcNcls=particle->GetTPCNcls();

  //GetNclsS not present in AODtrack
  //Replace with method as soon as available
  TBits tpcSharedMap = particle->GetTPCSharedMap();   
  Double_t tpcNclsS=  tpcSharedMap.CountBits(0)-tpcSharedMap.CountBits(159);
  Double_t tpcClusFindable=particle->GetTPCNclsF();

  // Reset AliESDtrack interface specific information
  values[AliDielectronVarManager::kNclsITS]       = particle->GetITSNcls();
  values[AliDielectronVarManager::kITSchi2Cl]     = -1;
  values[AliDielectronVarManager::kNclsTPC]       = tpcNcls;
  values[AliDielectronVarManager::kNclsSTPC]      = tpcNclsS;
  values[AliDielectronVarManager::kNclsSFracTPC]  = tpcNcls>0?tpcNclsS/tpcNcls:0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = tpcNcls; // not really available in AOD
  values[AliDielectronVarManager::kNFclsTPC]      = tpcClusFindable;
  values[AliDielectronVarManager::kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
  values[AliDielectronVarManager::kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
  values[AliDielectronVarManager::kNFclsTPCfCross]= (tpcClusFindable>0)?(particle->GetTPCClusterInfo(2,1)/tpcClusFindable):0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = particle->GetTRDntrackletsPID();
  values[AliDielectronVarManager::kTRDchi2]       = (particle->GetTRDntrackletsPID()!=0.?particle->GetTRDchi2():-1);
  values[AliDielectronVarManager::kTRDsignal]     = particle->GetTRDsignal();
  values[AliDielectronVarManager::kTPCclsSegments] = 0.0;
  const UChar_t threshold = 5;
  TBits tpcClusterMap = particle->GetTPCClusterMap();
  UChar_t n=0; UChar_t j=0;
  for(UChar_t i=0; i<8; ++i) {
    n=0;
    for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliDielectronVarManager::kTPCclsSegments] += 1.0;
  }

  values[AliDielectronVarManager::kTPCchi2Cl]     = (tpcNcls>0)?particle->Chi2perNDF()*(tpcNcls-5)/tpcNcls:-1.;  // it is stored as normalized to tpcNcls-5 (see AliAnalysisTaskESDfilter)
  values[AliDielectronVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  values[AliDielectronVarManager::kFilterBit]     = (Double_t)particle->GetFilterMap();

  //TRD pidProbs
  //TODO: set correctly
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;

  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kTPCsignalNfrac]= 0;

  // Fill AliAODTrack interface information
  //
  Int_t v0Index=-1;
  Int_t kinkIndex=-1;
  if (particle->GetProdVertex()) {
    v0Index   = particle->GetProdVertex()->GetType()==AliAODVertex::kV0   ? 1 : 0;
    kinkIndex = particle->GetProdVertex()->GetType()==AliAODVertex::kKink ? 1 : 0;
  }
  values[AliDielectronVarManager::kV0Index0]      = v0Index;
  values[AliDielectronVarManager::kKinkIndex0]    = kinkIndex;

  Double_t d0z0[2];
  GetDCA(particle, d0z0);
  values[AliDielectronVarManager::kImpactParXY]   = d0z0[0];
  values[AliDielectronVarManager::kImpactParZ]    = d0z0[1];

  values[AliDielectronVarManager::kPIn]            =  0.;
  values[AliDielectronVarManager::kTPCsignal]      =  0.;
  values[AliDielectronVarManager::kTPCsignalN]     = -1.;
  values[AliDielectronVarManager::kTPCsignalNfrac] = -1.;
  values[AliDielectronVarManager::kTPCclsDiff]     = -999.;
  
  values[AliDielectronVarManager::kTOFsignal]=0;
  //values[AliDielectronVarManager::kTOFbeta]=0;

  values[AliDielectronVarManager::kTPCnSigmaEle]=0;
  values[AliDielectronVarManager::kTPCnSigmaPio]=0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]=0;
  values[AliDielectronVarManager::kTPCnSigmaKao]=0;
  values[AliDielectronVarManager::kTPCnSigmaPro]=0;

  values[AliDielectronVarManager::kTOFnSigmaEle]=0;
  values[AliDielectronVarManager::kTOFnSigmaPio]=0;
  values[AliDielectronVarManager::kTOFnSigmaMuo]=0;
  values[AliDielectronVarManager::kTOFnSigmaKao]=0;
  values[AliDielectronVarManager::kTOFnSigmaPro]=0;
  
  values[AliDielectronVarManager::kITSclusterMap]   =   particle->GetITSClusterMap();
  values[AliDielectronVarManager::kITSLayerFirstCls] = -1.;
  for (Int_t iC=0; iC<6; iC++) {
    if (((particle->GetITSClusterMap()) & (1<<(iC))) > 0) {
      values[AliDielectronVarManager::kITSLayerFirstCls] = iC;
      break;
    }
  }

  AliAODPid *pid=const_cast<AliAODPid*>(particle->GetDetPid());
  if (pid) {
    Double_t origdEdx=pid->GetTPCsignal();
    //overwrite signal
    pid->SetTPCsignal(origdEdx/AliDielectronPID::GetEtaCorr(particle)/AliDielectronPID::GetCorrValdEdx());
    
    Double_t mom =pid->GetTPCmomentum();
    Double_t tpcSignalN=pid->GetTPCsignalN();
    values[AliDielectronVarManager::kTPCsignalN]     = tpcSignalN;
    values[AliDielectronVarManager::kTPCsignalNfrac] = tpcNcls>0?tpcSignalN/tpcNcls:0;
    values[AliDielectronVarManager::kTPCclsDiff]     = tpcSignalN-tpcNcls;
    Double_t tpcNsigmaEle=(fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron)-AliDielectronPID::GetCntrdCorr(particle))/AliDielectronPID::GetWdthCorr(particle);
    Double_t tpcNsigmaPio=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
    Double_t tpcNsigmaMuo=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
    Double_t tpcNsigmaKao=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
    Double_t tpcNsigmaPro=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);
    
    values[AliDielectronVarManager::kPIn]=mom;
    values[AliDielectronVarManager::kTPCsignal]=pid->GetTPCsignal();

    values[AliDielectronVarManager::kTPCnSigmaEle]=tpcNsigmaEle;
    values[AliDielectronVarManager::kTPCnSigmaPio]=tpcNsigmaPio;
    values[AliDielectronVarManager::kTPCnSigmaMuo]=tpcNsigmaMuo;
    values[AliDielectronVarManager::kTPCnSigmaKao]=tpcNsigmaKao;
    values[AliDielectronVarManager::kTPCnSigmaPro]=tpcNsigmaPro;
    
    Double_t prob[AliPID::kSPECIES];
    fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob);
    values[AliDielectronVarManager::kTRDprobEle]      = prob[AliPID::kElectron];
    values[AliDielectronVarManager::kTRDprobPio]      = prob[AliPID::kPion];
    //   fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob, AliTRDPIDResponse::kLQ2D);
    values[AliDielectronVarManager::kTRDprob2DEle]    = prob[AliPID::kElectron];
    values[AliDielectronVarManager::kTRDprob2DPio]    = prob[AliPID::kPion];

    values[AliDielectronVarManager::kTOFsignal]=pid->GetTOFsignal();
    Double32_t expt[5];
    particle->GetIntegratedTimes(expt);         // ps
    Double_t l  = TMath::C()* expt[0]*1e-12;    // m
    Double_t t  = pid->GetTOFsignal();          // ps start time subtracted (until v5-02-Rev09)
    AliTOFHeader* tofH=0x0;                     // from v5-02-Rev10 on subtract the start time
    if(fgEvent) tofH = (AliTOFHeader*)fgEvent->GetTOFHeader();
    if(tofH) t -= fgPIDResponse->GetTOFResponse().GetStartTime(particle->P()); // ps

    if( (l < 360.e-2 || l > 800.e-2) || (t <= 0.) ) {
      values[AliDielectronVarManager::kTOFbeta]  =0;
    }
    else {
      t *= 1e-12; //ps -> s
      
      Double_t v = l / t;
      Float_t beta = v / TMath::C();
      values[AliDielectronVarManager::kTOFbeta]=beta;
    }

    Double_t tofNsigmaEle=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kElectron);
    Double_t tofNsigmaPio=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kPion);
    Double_t tofNsigmaMuo=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kMuon);
    Double_t tofNsigmaKao=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kKaon);
    Double_t tofNsigmaPro=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kProton);
    
    values[AliDielectronVarManager::kTOFnSigmaEle]=tofNsigmaEle;
    values[AliDielectronVarManager::kTOFnSigmaPio]=tofNsigmaPio;
    values[AliDielectronVarManager::kTOFnSigmaMuo]=tofNsigmaMuo;
    values[AliDielectronVarManager::kTOFnSigmaKao]=tofNsigmaKao;
    values[AliDielectronVarManager::kTOFnSigmaPro]=tofNsigmaPro;

    values[AliDielectronVarManager::kTOFmismProb] = fgPIDResponse->GetTOFMismatchProbability(particle);
  
    pid->SetTPCsignal(origdEdx);
  }

  //EMCAL PID information
  Double_t eop=0;
  Double_t showershape[4]={0.,0.,0.,0.};
//   values[AliDielectronVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron,eop,showershape);
  values[AliDielectronVarManager::kEMCALEoverP]     = eop;
  values[AliDielectronVarManager::kEMCALE]          = eop*values[AliDielectronVarManager::kP];
  values[AliDielectronVarManager::kEMCALNCells]     = showershape[0];
  values[AliDielectronVarManager::kEMCALM02]        = showershape[1];
  values[AliDielectronVarManager::kEMCALM20]        = showershape[2];
  values[AliDielectronVarManager::kEMCALDispersion] = showershape[3];

  values[AliDielectronVarManager::kPdgCode]=-1;
  values[AliDielectronVarManager::kPdgCodeMother]=-1;
  values[AliDielectronVarManager::kPdgCodeGrandMother]=-1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;

  values[AliDielectronVarManager::kNumberOfDaughters]=-1;
  
  AliDielectronMC *mc=AliDielectronMC::Instance();
  
  if (mc->HasMC()){
    if (mc->GetMCTrack(particle)) {
      values[AliDielectronVarManager::kPdgCode]=mc->GetMCTrack(particle)->PdgCode();
      Int_t trkLbl = mc->GetMCTrack(particle)->GetLabel();
      //      printf("trklbl %d for %p->%p \n",trkLbl,particle,mc->GetMCTrack(particle));
      values[AliDielectronVarManager::kHasCocktailMother]=mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
    }
    AliAODMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
    if (motherMC){
      values[AliDielectronVarManager::kPdgCodeMother]=motherMC->PdgCode();
      Int_t motherLbl = motherMC->GetLabel();
      values[AliDielectronVarManager::kHasCocktailGrandMother]=mc->CheckParticleSource(motherLbl, AliDielectronSignalMC::kDirect);
     
      motherMC=mc->GetMCTrackMother(motherMC);  //grand motherMC
      if (motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=motherMC->PdgCode();;
    }
    
    values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
  } //if(mc->HasMC())
  
  values[AliDielectronVarManager::kTOFPIDBit]=(particle->GetStatus()&AliESDtrack::kTOFpid? 1: 0);
}

inline void AliDielectronVarManager::FillVarMCParticle(const AliMCParticle *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kITSchi2Cl]     = 0;
  values[AliDielectronVarManager::kNclsTPC]       = 0;
  values[AliDielectronVarManager::kNclsSTPC]      = 0;
  values[AliDielectronVarManager::kNclsSFracTPC]  = 0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = 0; 
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNFclsTPCr]     = 0;
  values[AliDielectronVarManager::kNFclsTPCrFrac] = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  values[AliDielectronVarManager::kTPCchi2Cl]     = 0;
  values[AliDielectronVarManager::kTrackStatus]   = 0;
  values[AliDielectronVarManager::kFilterBit]     = 0;
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kTPCclsDiff]    = 0;
  values[AliDielectronVarManager::kTPCsignalNfrac]    = 0;
  values[AliDielectronVarManager::kImpactParXY]   = 0;
  values[AliDielectronVarManager::kImpactParZ]    = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kYsignedIn]     = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTOFsignal]     = 0;
  values[AliDielectronVarManager::kTOFbeta]       = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;
  
  values[AliDielectronVarManager::kPdgCode]       = -1;
  values[AliDielectronVarManager::kPdgCodeMother] = -1;
  values[AliDielectronVarManager::kPdgCodeGrandMother] = -1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  
  AliDielectronMC *mc=AliDielectronMC::Instance();

  // Fill AliMCParticle interface specific information
  values[AliDielectronVarManager::kPdgCode] = particle->PdgCode();
  Int_t trkLbl = particle->GetLabel();
  values[AliDielectronVarManager::kHasCocktailMother]=mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);

  AliMCParticle *motherMC = mc->GetMCTrackMother(particle);
  if (motherMC){
    values[AliDielectronVarManager::kPdgCodeMother]=motherMC->PdgCode();
    Int_t motherLbl = motherMC->GetLabel();
    values[AliDielectronVarManager::kHasCocktailGrandMother]=mc->CheckParticleSource(motherLbl, AliDielectronSignalMC::kDirect);

    motherMC=mc->GetMCTrackMother(motherMC);  //grand mother
    if (motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=motherMC->PdgCode();;
  }
  
  values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);
  values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
}


inline void AliDielectronVarManager::FillVarMCParticle2(const AliVParticle *p1, const AliVParticle *p2, Double_t * const values) {
  //
  // fill 2 track information starting from MC legs
  //

  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kITSchi2Cl]     = -1;
  values[AliDielectronVarManager::kNclsTPC]       = 0;
  values[AliDielectronVarManager::kNclsSTPC]       = 0;
  values[AliDielectronVarManager::kNclsSFracTPC]  = 0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = 0; 
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNFclsTPCr]     = 0;
  values[AliDielectronVarManager::kNFclsTPCrFrac] = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  values[AliDielectronVarManager::kTPCchi2Cl]     = 0;
  values[AliDielectronVarManager::kTrackStatus]   = 0;
  values[AliDielectronVarManager::kFilterBit]     = 0;
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kTPCclsDiff]    = 0;
  values[AliDielectronVarManager::kTPCsignalNfrac]    = 0;
  values[AliDielectronVarManager::kImpactParXY]   = 0;
  values[AliDielectronVarManager::kImpactParZ]    = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kYsignedIn]     = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;
  
  values[AliDielectronVarManager::kPdgCode]       = 0;
  values[AliDielectronVarManager::kPdgCodeMother] = 0;
  values[AliDielectronVarManager::kHasCocktailMother]=0;

  AliDielectronMC *mc=AliDielectronMC::Instance();
  AliVParticle* mother=0x0;
  Int_t mLabel1 = mc->GetMothersLabel(p1->GetLabel());
  Int_t mLabel2 = mc->GetMothersLabel(p2->GetLabel());
  if(mLabel1==mLabel2)
    mother = mc->GetMCTrackFromMCEvent(mLabel1);

  values[AliDielectronVarManager::kPseudoProperTime] = -2e10;
  if(mother) {    // same mother
    FillVarVParticle(mother, values);
    Double_t vtxX, vtxY, vtxZ;
    mc->GetPrimaryVertex(vtxX,vtxY,vtxZ);
    Double_t lxy = ((mother->Xv()- vtxX) * mother->Px() + 
		    (mother->Yv()- vtxY) * mother->Py() )/mother->Pt();
    values[AliDielectronVarManager::kPseudoProperTime] = lxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/mother->Pt();
  }
  // AliVParticle part
  values[AliDielectronVarManager::kPx]        = p1->Px()+p2->Px();
  values[AliDielectronVarManager::kPy]        = p1->Py()+p2->Py();
  values[AliDielectronVarManager::kPz]        = p1->Pz()+p2->Pz();
  values[AliDielectronVarManager::kPt]        = TMath::Sqrt(values[AliDielectronVarManager::kPx]*
						      values[AliDielectronVarManager::kPx]+
						      values[AliDielectronVarManager::kPy]*
						      values[AliDielectronVarManager::kPy]);
  values[AliDielectronVarManager::kP]         = TMath::Sqrt(values[AliDielectronVarManager::kPt]*
						      values[AliDielectronVarManager::kPt]+
						      values[AliDielectronVarManager::kPz]*
						      values[AliDielectronVarManager::kPz]);
    
  values[AliDielectronVarManager::kXv]        = 0;
  values[AliDielectronVarManager::kYv]        = 0;
  values[AliDielectronVarManager::kZv]        = 0;
    
  values[AliDielectronVarManager::kOneOverPt] = (values[AliDielectronVarManager::kPt]>1.0e-6 ? 1.0/values[AliDielectronVarManager::kPt] : 0.0);
  values[AliDielectronVarManager::kPhi]       = TMath::ATan2(values[AliDielectronVarManager::kPy],values[AliDielectronVarManager::kPx]);
  values[AliDielectronVarManager::kTheta]     = TMath::ATan2(values[AliDielectronVarManager::kPt],values[AliDielectronVarManager::kPz]);
  values[AliDielectronVarManager::kEta]       = ((values[AliDielectronVarManager::kP]-values[AliDielectronVarManager::kPz])>1.0e-6 && (values[AliDielectronVarManager::kP]+values[AliDielectronVarManager::kPz])>1.0e-6 ? 0.5*TMath::Log((values[AliDielectronVarManager::kP]+values[AliDielectronVarManager::kPz])/(values[AliDielectronVarManager::kP]-values[AliDielectronVarManager::kPz])) : -9999.);
  values[AliDielectronVarManager::kE]         = p1->E()+p2->E();
  values[AliDielectronVarManager::kY]         = ((values[AliDielectronVarManager::kE]-values[AliDielectronVarManager::kPz])>1.0e-6 && (values[AliDielectronVarManager::kE]+values[AliDielectronVarManager::kPz])>1.0e-6 ? 0.5*TMath::Log((values[AliDielectronVarManager::kE]+values[AliDielectronVarManager::kPz])/(values[AliDielectronVarManager::kE]-values[AliDielectronVarManager::kPz])) : -9999.);
  values[AliDielectronVarManager::kCharge]    = p1->Charge()+p2->Charge();

  values[AliDielectronVarManager::kM]         = p1->M()*p1->M()+p2->M()*p2->M()+
                       2.0*(p1->E()*p2->E()-p1->Px()*p2->Px()-p1->Py()*p2->Py()-p1->Pz()*p2->Pz());
  values[AliDielectronVarManager::kM]         = (values[AliDielectronVarManager::kM]>1.0e-8 ? TMath::Sqrt(values[AliDielectronVarManager::kM]) : -1.0);

  if ( fgEvent ) AliDielectronVarManager::Fill(fgEvent, values);  

  values[AliDielectronVarManager::kThetaHE]   = AliDielectronPair::ThetaPhiCM(p1,p2,kTRUE,  kTRUE);
  values[AliDielectronVarManager::kPhiHE]     = AliDielectronPair::ThetaPhiCM(p1,p2,kTRUE,  kFALSE);
  values[AliDielectronVarManager::kThetaSqHE]  = values[AliDielectronVarManager::kThetaHE] * values[AliDielectronVarManager::kThetaHE];
  values[AliDielectronVarManager::kCos2PhiHE] = TMath::Cos(2*values[AliDielectronVarManager::kPhiHE]);
  values[AliDielectronVarManager::kThetaCS]   = AliDielectronPair::ThetaPhiCM(p1,p2,kFALSE, kTRUE);
  values[AliDielectronVarManager::kPhiCS]     = AliDielectronPair::ThetaPhiCM(p1,p2,kFALSE, kFALSE);
  values[AliDielectronVarManager::kThetaSqCS]  = values[AliDielectronVarManager::kThetaCS] * values[AliDielectronVarManager::kThetaCS];
  values[AliDielectronVarManager::kCos2PhiCS] = TMath::Cos(2*values[AliDielectronVarManager::kPhiCS]);
  values[AliDielectronVarManager::kCosTilPhiHE]  = (values[AliDielectronVarManager::kThetaHE]>0)?(TMath::Cos(values[AliDielectronVarManager::kPhiHE]-TMath::Pi()/4.)):(TMath::Cos(values[AliDielectronVarManager::kPhiHE]-3*TMath::Pi()/4.));
  values[AliDielectronVarManager::kCosTilPhiCS]  = (values[AliDielectronVarManager::kThetaCS]>0)?(TMath::Cos(values[AliDielectronVarManager::kPhiCS]-TMath::Pi()/4.)):(TMath::Cos(values[AliDielectronVarManager::kPhiCS]-3*TMath::Pi()/4.));
}


inline void AliDielectronVarManager::FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kITSchi2Cl]     = -1;
  values[AliDielectronVarManager::kNclsTPC]       = 0;
  values[AliDielectronVarManager::kNclsSTPC]      = 0;
  values[AliDielectronVarManager::kNclsSFracTPC]  = 0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = 0;
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  values[AliDielectronVarManager::kTPCchi2Cl]     = 0;
  values[AliDielectronVarManager::kTrackStatus]   = 0;
  values[AliDielectronVarManager::kFilterBit]     = 0;
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kTPCclsDiff]    = 0;
  values[AliDielectronVarManager::kTPCsignalNfrac]= 0;
  values[AliDielectronVarManager::kImpactParXY]   = 0;
  values[AliDielectronVarManager::kImpactParZ]    = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kYsignedIn]     = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;
  
  values[AliDielectronVarManager::kPdgCode]       = -1;
  values[AliDielectronVarManager::kPdgCodeMother] = -1;
  values[AliDielectronVarManager::kPdgCodeGrandMother] = -1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);
  
  AliDielectronMC *mc=AliDielectronMC::Instance();


  // Fill AliAODMCParticle interface specific information
  values[AliDielectronVarManager::kPdgCode] = particle->PdgCode();
  Int_t trkLbl = particle->GetLabel();
  values[AliDielectronVarManager::kHasCocktailMother]=mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);

  AliAODMCParticle *motherMC = mc->GetMCTrackMother(particle);
  if (motherMC){
    values[AliDielectronVarManager::kPdgCodeMother]=motherMC->PdgCode();
    Int_t motherLbl = motherMC->GetLabel();
    values[AliDielectronVarManager::kHasCocktailGrandMother]=mc->CheckParticleSource(motherLbl, AliDielectronSignalMC::kDirect);

    motherMC=mc->GetMCTrackMother(motherMC);  //grand mother
    if (motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=motherMC->PdgCode();;
  }
  
  
  values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);

  values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
}

inline void AliDielectronVarManager::FillVarDielectronPair(const AliDielectronPair *pair, Double_t * const values)
{
  //
  // Fill pair information available for histogramming into an array
  //
  
  values[AliDielectronVarManager::kPdgCode]=-1;
  values[AliDielectronVarManager::kPdgCodeMother]=-1;
  values[AliDielectronVarManager::kPdgCodeGrandMother]=-1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;
  
  Double_t errPseudoProperTime2 = -1;
  // Fill common AliVParticle interface information
  FillVarVParticle(pair, values);

  // Fill AliDielectronPair specific information
  const AliKFParticle &kfPair = pair->GetKFParticle();

  Double_t thetaHE=0;
  Double_t phiHE=0;
  Double_t thetaCS=0;
  Double_t phiCS=0;

  pair->GetThetaPhiCM(thetaHE,phiHE,thetaCS,phiCS);
  
  values[AliDielectronVarManager::kChi2NDF]      = kfPair.GetChi2()/kfPair.GetNDF();
  values[AliDielectronVarManager::kDecayLength]  = kfPair.GetDecayLength();
  values[AliDielectronVarManager::kR]            = kfPair.GetR();
  values[AliDielectronVarManager::kOpeningAngle] = pair->OpeningAngle();
  values[AliDielectronVarManager::kCosPointingAngle] = fgEvent ? pair->GetCosPointingAngle(fgEvent->GetPrimaryVertex()) : -1;
  values[AliDielectronVarManager::kThetaHE]      = thetaHE;
  values[AliDielectronVarManager::kPhiHE]        = phiHE;
  values[AliDielectronVarManager::kThetaSqHE]    = thetaHE * thetaHE;
  values[AliDielectronVarManager::kCos2PhiHE]    = TMath::Cos(2.0*phiHE);
  values[AliDielectronVarManager::kCosTilPhiHE]  = (thetaHE>0)?(TMath::Cos(phiHE-TMath::Pi()/4.)):(TMath::Cos(phiHE-3*TMath::Pi()/4.));
  values[AliDielectronVarManager::kThetaCS]      = thetaCS;
  values[AliDielectronVarManager::kPhiCS]        = phiCS;
  values[AliDielectronVarManager::kThetaSqCS]    = thetaCS * thetaCS;
  values[AliDielectronVarManager::kCos2PhiCS]    = TMath::Cos(2.0*phiCS);
  values[AliDielectronVarManager::kCosTilPhiCS]  = (thetaCS>0)?(TMath::Cos(phiCS-TMath::Pi()/4.)):(TMath::Cos(phiCS-3*TMath::Pi()/4.));
  values[AliDielectronVarManager::kLegDist]      = pair->DistanceDaughters();
  values[AliDielectronVarManager::kLegDistXY]    = pair->DistanceDaughtersXY();
  values[AliDielectronVarManager::kDeltaEta]     = pair->DeltaEta();
  values[AliDielectronVarManager::kDeltaPhi]     = pair->DeltaPhi();
  values[AliDielectronVarManager::kMerr]         = kfPair.GetErrMass()>1e-30&&kfPair.GetMass()>1e-30?kfPair.GetErrMass()/kfPair.GetMass():1000000;
  values[AliDielectronVarManager::kPairType]     = pair->GetType();
  // Armenteros-Podolanski quantities
  values[AliDielectronVarManager::kArmAlpha]     = pair->GetArmAlpha();
  values[AliDielectronVarManager::kArmPt]        = pair->GetArmPt();

  values[AliDielectronVarManager::kPsiPair]      = fgEvent ? pair->PsiPair(fgEvent->GetMagneticField()) : -5;
  values[AliDielectronVarManager::kPhivPair]      = fgEvent ? pair->PhivPair(fgEvent->GetMagneticField()) : -5;
  values[AliDielectronVarManager::kPseudoProperTime] = fgEvent ? kfPair.GetPseudoProperDecayTime(*(fgEvent->GetPrimaryVertex()), TDatabasePDG::Instance()->GetParticle(443)->Mass(), &errPseudoProperTime2 ) : -1e10;
  // values[AliDielectronVarManager::kPseudoProperTime] = fgEvent ? pair->GetPseudoProperTime(fgEvent->GetPrimaryVertex()): -1e10;
  values[AliDielectronVarManager::kPseudoProperTimeErr] = (errPseudoProperTime2 > 0) ? TMath::Sqrt(errPseudoProperTime2) : -1e10;

 
  if (!(pair->GetKFUsage())) {
	//if KF Pairing is not enabled, overwrite values that can be easily derived from legs
	//use the INDIVIDUAL KF particles as source, which should be a copy of the corresponding properties
	//the ESDtrack, the reference to the ESDtrack is not (always) accessible in Mixing, while KF
	//particles are copied in the Pair-Object
	static const Double_t mElectron = AliPID::ParticleMass(AliPID::kElectron); // MeV

	const AliKFParticle& fD1 = pair->GetKFFirstDaughter();
	const AliKFParticle& fD2 = pair->GetKFSecondDaughter();

	//Define local buffer variables for leg properties
	Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
	Double_t px2=-9999.,py2=-9999.,pz2=-9999.;
	Double_t e1 =-9999.,e2 =-9999.;
	Double_t feta1=-9999.;//,fphi1=-9999.;
	Double_t feta2=-9999.;//,fphi2=-9999.;

	px1 = fD1.GetPx(); 
	py1 = fD1.GetPy(); 
	pz1 = fD1.GetPz(); 
	feta1 = fD1.GetEta();
	//	fphi1 = fD1.GetPhi();

	px2 = fD2.GetPx(); 
	py2 = fD2.GetPy(); 
	pz2 = fD2.GetPz(); 
	feta2 = fD2.GetEta();
	//	fphi2 = fD2.GetPhi();

	//Calculate Energy per particle by hand
	e1 = TMath::Sqrt(mElectron*mElectron+px1*px1+py1*py1+pz1*pz1);
	e2 = TMath::Sqrt(mElectron*mElectron+px2*px2+py2*py2+pz2*pz2);

	//Now Create TLorentzVector:
	TLorentzVector lv1,lv2;
	lv1.SetPxPyPzE(px1,py1,pz1,e1);
	lv2.SetPxPyPzE(px2,py2,pz2,e2);

	values[AliDielectronVarManager::kPx]        = (lv1+lv2).Px();
	values[AliDielectronVarManager::kPy]        = (lv1+lv2).Py();
	values[AliDielectronVarManager::kPz]        = (lv1+lv2).Pz();

	values[AliDielectronVarManager::kPt]        =  (lv1+lv2).Pt();

	values[AliDielectronVarManager::kP]         =  (lv1+lv2).P();

	//Not overwritten, could take event vertex in next iteration
	values[AliDielectronVarManager::kXv]        = (lv1+lv2).X(); 
	values[AliDielectronVarManager::kYv]        = (lv1+lv2).Y();
	values[AliDielectronVarManager::kZv]        = (lv1+lv2).Z();

	values[AliDielectronVarManager::kE]         = (lv1+lv2).E();


	values[AliDielectronVarManager::kM]         = (lv1+lv2).M();

	values[AliDielectronVarManager::kOpeningAngle] =  lv1.Angle(lv2.Vect());

	values[AliDielectronVarManager::kOneOverPt] = (values[AliDielectronVarManager::kPt]>0. ? 1./values[AliDielectronVarManager::kPt] : -9999.);
	values[AliDielectronVarManager::kPhi]       = (lv1+lv2).Phi();
	values[AliDielectronVarManager::kEta]       = (lv1+lv2).Eta();

	values[AliDielectronVarManager::kY]       = (lv1+lv2).Rapidity();

	for (Int_t i=AliDielectronVarManager::kPairMax; i<AliDielectronVarManager::kNMaxValues; ++i)
	  values[i]=fgData[i];

	// Fill AliDielectronPair specific information
	values[AliDielectronVarManager::kDeltaEta]     = TMath::Abs(feta1 -feta2 );
	values[AliDielectronVarManager::kDeltaPhi]     = lv1.DeltaPhi(lv2);
	values[AliDielectronVarManager::kPairType]     = pair->GetType();

	/*
	//Also not overwritten, still coming from KF particle
	//where needed to be replaced by independent determination
	values[AliDielectronVarManager::kCharge]    = 0.;
	values[AliDielectronVarManager::kPdgCode]   = 0.;
	values[AliDielectronVarManager::kChi2NDF]      = 0.;
	values[AliDielectronVarManager::kDecayLength]  = 0.;
	values[AliDielectronVarManager::kR]            = 0.;
	values[AliDielectronVarManager::kCosPointingAngle] = 0.;
	values[AliDielectronVarManager::kThetaHE]      = 0.;
	values[AliDielectronVarManager::kPhiHE]        = 0.;
	values[AliDielectronVarManager::kThetaSqHE]    = 0.;
	values[AliDielectronVarManager::kCos2PhiHE]    = 0.;
	values[AliDielectronVarManager::kCosTilPhiHE]  = 0.;
	values[AliDielectronVarManager::kThetaCS]      = 0.;
	values[AliDielectronVarManager::kPhiCS]        = 0.;
	values[AliDielectronVarManager::kThetaSqCS]    = 0.;
	values[AliDielectronVarManager::kCos2PhiCS]    = 0.;
	values[AliDielectronVarManager::kCosTilPhiCS]  = 0.;
	values[AliDielectronVarManager::kLegDist]      = 0.;
	values[AliDielectronVarManager::kLegDistXY]    = 0.;
	values[AliDielectronVarManager::kMerr]         = 0.;
	values[AliDielectronVarManager::kPseudoProperTime] = 0.;
	values[AliDielectronVarManager::kPseudoProperTimeErr] = 0.;
	//Fill in Taku's PhiV?
	values[AliDielectronVarManager::kPsiPair]      = 0.;

	 */
  }
  //common, regardless of calculation method 
   // Flow quantities
  Double_t delta=0.0;
  // v2 with respect to VZERO-A event plane
  delta = values[AliDielectronVarManager::kPhi] - fgData[AliDielectronVarManager::kV0ArpH2];
  if(delta>TMath::Pi()) delta -= 2.0*TMath::Pi();             // keep the [-pi,+pi] interval
  if(delta<-1.0*TMath::Pi()) delta += 2.0*TMath::Pi();
  values[AliDielectronVarManager::kV0ArpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
  values[AliDielectronVarManager::kDeltaPhiV0ArpH2] = delta;
  // v2 with respect to VZERO-C event plane
  delta = values[AliDielectronVarManager::kPhi] - fgData[AliDielectronVarManager::kV0CrpH2];
  if(delta>TMath::Pi()) delta -= 2.0*TMath::Pi();             // keep the [-pi,+pi] interval
  if(delta<-1.0*TMath::Pi()) delta += 2.0*TMath::Pi();
  values[AliDielectronVarManager::kV0CrpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
  values[AliDielectronVarManager::kDeltaPhiV0CrpH2] = delta;
  // v2 with respect to the combined VZERO-A and VZERO-C event plane
  delta = values[AliDielectronVarManager::kPhi] - fgData[AliDielectronVarManager::kV0ACrpH2];
  if(delta>TMath::Pi()) delta -= 2.0*TMath::Pi();             // keep the [-pi,+pi] interval
  if(delta<-1.0*TMath::Pi()) delta += 2.0*TMath::Pi();
  values[AliDielectronVarManager::kV0ACrpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
  values[AliDielectronVarManager::kDeltaPhiV0ACrpH2] = delta;


  // quantities using the values of  AliEPSelectionTask
  values[AliDielectronVarManager::kDeltaPhiv0ArpH2]  = values[AliDielectronVarManager::kPhi] - values[AliDielectronVarManager::kv0ArpH2];
  values[AliDielectronVarManager::kDeltaPhiv0CrpH2]  = values[AliDielectronVarManager::kPhi] - values[AliDielectronVarManager::kv0CrpH2];
  values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] = values[AliDielectronVarManager::kPhi] - values[AliDielectronVarManager::kv0ACrpH2];
  values[AliDielectronVarManager::kv0ACrpH2FlowV2]   = TMath::Cos( 2*(values[AliDielectronVarManager::kPhi] - values[AliDielectronVarManager::kv0ACrpH2]) );
  values[AliDielectronVarManager::kv0ArpH2FlowV2]    = TMath::Cos( 2*(values[AliDielectronVarManager::kPhi] - values[AliDielectronVarManager::kv0ArpH2]) );
  values[AliDielectronVarManager::kv0CrpH2FlowV2]    = TMath::Cos( 2*(values[AliDielectronVarManager::kPhi] - values[AliDielectronVarManager::kv0CrpH2]) );

  // keep the interval [-pi,+pi]
  if ( values[AliDielectronVarManager::kDeltaPhiv0ArpH2] > TMath::Pi() ) 
    values[AliDielectronVarManager::kDeltaPhiv0ArpH2] -= TMath::TwoPi(); 
  if ( values[AliDielectronVarManager::kDeltaPhiv0CrpH2] > TMath::Pi() ) 
    values[AliDielectronVarManager::kDeltaPhiv0CrpH2] -= TMath::TwoPi(); 
  if ( values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] > TMath::Pi() ) 
    values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] -= TMath::TwoPi(); 

  if ( values[AliDielectronVarManager::kDeltaPhiv0ArpH2] < -1.*TMath::Pi() ) 
    values[AliDielectronVarManager::kDeltaPhiv0ArpH2] += TMath::TwoPi(); 
  if ( values[AliDielectronVarManager::kDeltaPhiv0CrpH2] < -1.*TMath::Pi() ) 
    values[AliDielectronVarManager::kDeltaPhiv0CrpH2] += TMath::TwoPi(); 
  if ( values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] < -1.*TMath::Pi() )
    values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] += TMath::TwoPi(); 

  //angle between ee plane and Mag/Reaction plane
  values[AliDielectronVarManager::kPairPlanev0rpH2Angle] = pair->PairPlanev0rpH2Angle(values[AliDielectronVarManager::kv0CrpH2]);
  values[AliDielectronVarManager::kPairPlaneMagAngle] = pair->PairPlaneMagAngle(values[AliDielectronVarManager::kv0CrpH2]);
  values[AliDielectronVarManager::kPairPlaneAngle] = pair->PairPlaneAngle(values[AliDielectronVarManager::kv0CrpH2]);
  //ee plane vector
  Double_t RotPairx = 0;
  Double_t RotPairy = 0;
  Double_t RotPairz = 0;
  pair->GetRotPair(RotPairx,RotPairy,RotPairz);
  values[AliDielectronVarManager::kRotPairx]=RotPairx;
  values[AliDielectronVarManager::kRotPairy]=RotPairy;
  values[AliDielectronVarManager::kRotPairz]=RotPairz;


  AliDielectronMC *mc=AliDielectronMC::Instance();
  
  if (mc->HasMC()){
    values[AliDielectronVarManager::kPseudoProperTimeResolution] = -10.0e+10;
    Bool_t samemother =  mc->HaveSameMother(pair);
    values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(pair);
    values[AliDielectronVarManager::kHaveSameMother] = samemother ;

    // fill kPseudoProperTimeResolution
    values[AliDielectronVarManager::kPseudoProperTimeResolution] = -1e10;
    // values[AliDielectronVarManager::kPseudoProperTimePull] = -1e10;
    if(samemother && fgEvent) {
      if(pair->GetFirstDaughter()->GetLabel() > 0) {
        const AliVParticle *motherMC = 0x0;
        if(fgEvent->IsA() == AliESDEvent::Class())  motherMC = (AliMCParticle*)mc->GetMCTrackMother((AliESDtrack*)pair->GetFirstDaughter());
        else if(fgEvent->IsA() == AliAODEvent::Class())  motherMC = (AliAODMCParticle*)mc->GetMCTrackMother((AliAODTrack*)pair->GetFirstDaughter());
        Double_t vtxX, vtxY, vtxZ;
	if(motherMC && mc->GetPrimaryVertex(vtxX,vtxY,vtxZ)) {
	  Int_t motherLbl = motherMC->GetLabel();
	  values[AliDielectronVarManager::kHasCocktailMother]=mc->CheckParticleSource(motherLbl, AliDielectronSignalMC::kDirect);
      	  const Double_t lxyMC = ( (motherMC->Xv() - vtxX) * motherMC->Px() +
                                   (motherMC->Yv() - vtxY) * motherMC->Py()   ) / motherMC->Pt();
	  const Double_t pseudoMC = lxyMC * (TDatabasePDG::Instance()->GetParticle(443)->Mass())/motherMC->Pt();
	  values[AliDielectronVarManager::kPseudoProperTimeResolution] = values[AliDielectronVarManager::kPseudoProperTime] - pseudoMC;
          if (errPseudoProperTime2 > 0)
            values[AliDielectronVarManager::kPseudoProperTimePull] = values[AliDielectronVarManager::kPseudoProperTimeResolution]/sqrt(errPseudoProperTime2);
      }
      }
    }
    
	values[AliDielectronVarManager::kTRDpidEffPair] = 0.;
	if (fgTRDpidEff[0][0]){
	  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
	  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
	  AliVParticle* leg1 = pair->GetFirstDaughter();
	  AliVParticle* leg2 = pair->GetSecondDaughter();
	  if (leg1 && leg2){
		Fill(leg1, valuesLeg1);
		Fill(leg2, valuesLeg2);
		values[AliDielectronVarManager::kTRDpidEffPair] = valuesLeg1[AliDielectronVarManager::kTRDpidEffLeg]*valuesLeg2[AliDielectronVarManager::kTRDpidEffLeg];
	  }
	}


  }//if (mc->HasMC())

  AliVParticle* leg1 = pair->GetFirstDaughter();
  AliVParticle* leg2 = pair->GetSecondDaughter();
  if (leg1)
	values[AliDielectronVarManager::kMomAsymDau1] = (values[AliDielectronVarManager::kP] != 0)? leg1->P()  / values[AliDielectronVarManager::kP]: 0;
  else 
	values[AliDielectronVarManager::kMomAsymDau1] = -9999.;
  if (leg2)
	values[AliDielectronVarManager::kMomAsymDau2] = (values[AliDielectronVarManager::kP] != 0)? leg2->P()  / values[AliDielectronVarManager::kP]: 0;
  else 
	values[AliDielectronVarManager::kMomAsymDau2] = -9999.;
}

inline void AliDielectronVarManager::FillVarKFParticle(const AliKFParticle *particle, Double_t * const values)
{
  //
  // Fill track information available in AliVParticle into an array
  //
  values[AliDielectronVarManager::kPx]        = particle->GetPx();
  values[AliDielectronVarManager::kPy]        = particle->GetPy();
  values[AliDielectronVarManager::kPz]        = particle->GetPz();
  values[AliDielectronVarManager::kPt]        = particle->GetPt();
  values[AliDielectronVarManager::kP]         = particle->GetP();
  
  values[AliDielectronVarManager::kXv]        = particle->GetX();
  values[AliDielectronVarManager::kYv]        = particle->GetY();
  values[AliDielectronVarManager::kZv]        = particle->GetZ();
  
  values[AliDielectronVarManager::kOneOverPt] = 0;
  values[AliDielectronVarManager::kPhi]       = particle->GetPhi();
  values[AliDielectronVarManager::kTheta]     = 0.;
  values[AliDielectronVarManager::kEta]       = particle->GetEta();
  values[AliDielectronVarManager::kY]         = ((particle->GetE()*particle->GetE()-particle->GetPx()*particle->GetPx()-particle->GetPy()*particle->GetPy()-particle->GetPz()*particle->GetPz())>0.) ? TLorentzVector(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE()).Rapidity() : -1111.;
  
  values[AliDielectronVarManager::kE]         = particle->GetE();
  values[AliDielectronVarManager::kM]         = particle->GetMass();
  values[AliDielectronVarManager::kCharge]    = particle->GetQ();
  
  values[AliDielectronVarManager::kNclsITS]       = 0;
  values[AliDielectronVarManager::kITSchi2Cl]     = -1;
  values[AliDielectronVarManager::kNclsTPC]       = 0;
  values[AliDielectronVarManager::kNclsSTPC]      = 0;
  values[AliDielectronVarManager::kNclsSFracTPC]  = 0;
  values[AliDielectronVarManager::kNclsTPCiter1]  = 0;
  values[AliDielectronVarManager::kNFclsTPC]      = 0;
  values[AliDielectronVarManager::kNclsTRD]       = 0;
  values[AliDielectronVarManager::kTRDntracklets] = 0;
  values[AliDielectronVarManager::kTRDpidQuality] = 0;
  values[AliDielectronVarManager::kTPCchi2Cl]     = 0;
  values[AliDielectronVarManager::kTrackStatus]   = 0;
  values[AliDielectronVarManager::kFilterBit]     = 0;
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;
  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kTPCclsDiff]    = 0;
  values[AliDielectronVarManager::kTPCsignalNfrac]= 0;
  values[AliDielectronVarManager::kImpactParXY]   = 0;
  values[AliDielectronVarManager::kImpactParZ]    = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kYsignedIn]     = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTOFsignal]     = 0;
  values[AliDielectronVarManager::kTOFbeta]       = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;
  
  values[AliDielectronVarManager::kPdgCode]       = -1;
  values[AliDielectronVarManager::kPdgCodeMother] = -1;
  values[AliDielectronVarManager::kPdgCodeGrandMother] = -1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;
  
//   if ( fgEvent ) AliDielectronVarManager::Fill(fgEvent, values);
  for (Int_t i=AliDielectronVarManager::kPairMax; i<AliDielectronVarManager::kNMaxValues; ++i)
    values[i]=fgData[i];
  
}

inline void AliDielectronVarManager::FillVarVEvent(const AliVEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //
  values[AliDielectronVarManager::kRunNumber]    = event->GetRunNumber();
  if(fgCurrentRun!=event->GetRunNumber()) {
    if(fgVZEROCalibrationFile.Contains(".root")) InitVZEROCalibrationHistograms(event->GetRunNumber());
    if(fgVZERORecenteringFile.Contains(".root")) InitVZERORecenteringHistograms(event->GetRunNumber());
    fgCurrentRun=event->GetRunNumber();
  }
  values[AliDielectronVarManager::kMixingBin]=0;

  const AliVVertex *primVtx = event->GetPrimaryVertex();
  
  values[AliDielectronVarManager::kXvPrim]       = 0;
  values[AliDielectronVarManager::kYvPrim]       = 0;
  values[AliDielectronVarManager::kZvPrim]       = 0;
  values[AliDielectronVarManager::kNVtxContrib]  = 0;
//   values[AliDielectronVarManager::kChi2NDF]      = 0; //This is the pair value!!!

  values[AliDielectronVarManager::kNTrk]            = 0;
  values[AliDielectronVarManager::kNVtxContrib]     = 0;
  values[AliDielectronVarManager::kNacc]            = 0;
  values[AliDielectronVarManager::kNaccTrcklts]     = 0;
  values[AliDielectronVarManager::kNaccTrcklts0916] = 0;
  values[AliDielectronVarManager::kNevents]         = 0; //always fill bin 0;
  values[AliDielectronVarManager::kRefMult]         = 0;
  values[AliDielectronVarManager::kRefMultTPConly]  = 0;
  
  if (primVtx){
    values[AliDielectronVarManager::kXvPrim]       = primVtx->GetX();
    values[AliDielectronVarManager::kYvPrim]       = primVtx->GetY();
    values[AliDielectronVarManager::kZvPrim]       = primVtx->GetZ();
    values[AliDielectronVarManager::kNVtxContrib]  = primVtx->GetNContributors();
  }
  //   values[AliDielectronVarManager::kChi2NDF]      = primVtx->GetChi2perNDF(); //this is the pair value
  
  values[AliDielectronVarManager::kNTrk]            = event->GetNumberOfTracks();
  values[AliDielectronVarManager::kNacc]            = AliDielectronHelper::GetNacc(event);
  values[AliDielectronVarManager::kMatchEffITSTPC]  = AliDielectronHelper::GetITSTPCMatchEff(event);
  values[AliDielectronVarManager::kNaccTrcklts]     = AliDielectronHelper::GetNaccTrcklts(event);      // etaRange = 1.6 (default)
  values[AliDielectronVarManager::kNaccTrcklts0916] = AliDielectronHelper::GetNaccTrcklts(event,1.6)-AliDielectronHelper::GetNaccTrcklts(event,.9);
  //  values[AliDielectronVarManager::kNaccTrcklts05]   = AliDielectronHelper::GetNaccTrcklts(event, 0.5);
  //  values[AliDielectronVarManager::kNaccTrcklts10]   = AliDielectronHelper::GetNaccTrcklts(event, 1.0);
  //  values[AliDielectronVarManager::kNaccTrckltsCorr] = AliDielectronHelper::GetNaccTrckltsCorrected(event, values[AliDielectronVarManager::kNaccTrcklts], values[AliDielectronVarManager::kZvPrim]);

  values[AliDielectronVarManager::kPhiMaxPt]          = 0;
  values[AliDielectronVarManager::kMaxPt]             = 0;
  Double_t ptMaxEv    = -1., phiptMaxEv= -1.;
  for(Int_t itrk=0; itrk<event->GetNumberOfTracks(); itrk++) {
    AliVParticle *part= event->GetTrack(itrk);
    if(part->Pt() > ptMaxEv) {
      ptMaxEv    = part->Pt();
      phiptMaxEv = part->Phi();
    }
  }
  values[AliDielectronVarManager::kPhiMaxPt]          = phiptMaxEv;
  values[AliDielectronVarManager::kMaxPt]             = ptMaxEv;


  // event plane quantities from the AliEPSelectionTask
  for(Int_t ivar=AliDielectronVarManager::kv0ArpH2; ivar<=kv0C0v0C3DiffH2;   ivar++) values[ivar] = 0.0; // v0  variables
  for(Int_t ivar=AliDielectronVarManager::kTPCxH2;  ivar<=kTPCsub12DiffH2uc; ivar++) values[ivar] = 0.0; // tpc variables

  // ep angle interval [todo, fill]
  AliEventplane *ep = const_cast<AliVEvent*>(event)->GetEventplane();
  if(ep) {

    // TPC event plane quantities (uncorrected)
    TVector2 *qstd  = ep->GetQVector();  // This is the "standard" Q-Vector for TPC
    TVector2 *qsub1 = ep->GetQsub1();    // random subevent plane
    TVector2 *qsub2 = ep->GetQsub2();
    if(qstd && qsub1 && qsub2) {
      values[AliDielectronVarManager::kTPCxH2uc]       = qstd->X();
      values[AliDielectronVarManager::kTPCyH2uc]       = qstd->Y();
      values[AliDielectronVarManager::kTPCmagH2uc]     = qstd->Mod();
      values[AliDielectronVarManager::kTPCrpH2uc]      = ((TMath::Abs(qstd->X())>1.0e-10) ? TMath::ATan2(qstd->Y(),qstd->X())/2.0 : 0.0);
      values[AliDielectronVarManager::kTPCsub1xH2uc]   = qsub1->X();
      values[AliDielectronVarManager::kTPCsub1yH2uc]   = qsub1->Y();
      values[AliDielectronVarManager::kTPCsub1rpH2uc]  = ((TMath::Abs(qsub1->X())>1.0e-10) ? TMath::ATan2(qsub1->Y(),qsub1->X())/2.0 : 0.0);
      values[AliDielectronVarManager::kTPCsub2xH2uc]   = qsub2->X();
      values[AliDielectronVarManager::kTPCsub2yH2uc]   = qsub2->Y();
      values[AliDielectronVarManager::kTPCsub2rpH2uc]  = ((TMath::Abs(qsub2->X())>1.0e-10) ? TMath::ATan2(qsub2->Y(),qsub2->X())/2.0 : 0.0);

      values[AliDielectronVarManager::kTPCsub12DiffH2uc] = TMath::Cos( 2.*(values[AliDielectronVarManager::kTPCsub1rpH2uc] -
									   values[AliDielectronVarManager::kTPCsub2rpH2uc]) );
    }

    // VZERO event plane
    TVector2 qvec;
    Double_t qx = 0, qy = 0;
    ep->CalculateVZEROEventPlane(event,10, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ACrpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    values[AliDielectronVarManager::kv0ACxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0ACyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0ACmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 8, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ArpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    values[AliDielectronVarManager::kv0AxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0AyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0AmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 9, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0CrpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    values[AliDielectronVarManager::kv0CxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0CyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0CmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 0, 0, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C0rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    ep->CalculateVZEROEventPlane(event, 3, 3, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C3rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    ep->CalculateVZEROEventPlane(event, 4, 4, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A0rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    ep->CalculateVZEROEventPlane(event, 7, 7, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A3rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
  } //if: eventplane

  // ESD VZERO information
  AliVVZERO* vzeroData = event->GetVZEROData();
  values[AliDielectronVarManager::kMultV0A] = 0.0;
  values[AliDielectronVarManager::kMultV0C] = 0.0;
  values[AliDielectronVarManager::kAdcV0A]  = 0.0;
  values[AliDielectronVarManager::kAdcV0C]  = 0.0;
  for(Int_t i=0; i<32; ++i) {
    values[AliDielectronVarManager::kVZEROchMult+i] = vzeroData->GetMultiplicity(i);
    values[AliDielectronVarManager::kVZEROchMult+32+i] = vzeroData->GetMultiplicity(i+32);
    //values[AliDielectronVarManager::kVZEROchMult+i] = event->GetVZEROEqMultiplicity(i);
    //values[AliDielectronVarManager::kVZEROchMult+32+i] = event->GetVZEROEqMultiplicity(i+32);
    values[AliDielectronVarManager::kMultV0A] += vzeroData->GetMultiplicityV0A(i);
    values[AliDielectronVarManager::kMultV0C] += vzeroData->GetMultiplicityV0C(i);
    //values[AliDielectronVarManager::kAdcV0A] += vzeroData->GetAdcV0A(i);
    //values[AliDielectronVarManager::kAdcV0C] += vzeroData->GetAdcV0C(i);
  }
  values[AliDielectronVarManager::kMultV0] = values[AliDielectronVarManager::kMultV0A] + values[AliDielectronVarManager::kMultV0C];
  values[AliDielectronVarManager::kAdcV0] = values[AliDielectronVarManager::kAdcV0A] + values[AliDielectronVarManager::kAdcV0C];
  // VZERO event plane quantities
  Double_t qvec[3]={0.0};
  GetVzeroRP(event, qvec,0);      // V0-A
  values[AliDielectronVarManager::kV0AxH2] = qvec[0]; values[AliDielectronVarManager::kV0AyH2] = qvec[1]; 
  values[AliDielectronVarManager::kV0ArpH2] = qvec[2];
  qvec[0]=0.0; qvec[1]=0.0; qvec[2]=0.0;
  GetVzeroRP(event, qvec,1);      // V0-C
  values[AliDielectronVarManager::kV0CxH2] = qvec[0]; values[AliDielectronVarManager::kV0CyH2] = qvec[1]; 
  values[AliDielectronVarManager::kV0CrpH2] = qvec[2];
  qvec[0]=0.0; qvec[1]=0.0; qvec[2]=0.0;
  GetVzeroRP(event, qvec,2);      // V0-A and V0-C combined
  values[AliDielectronVarManager::kV0ACxH2] = qvec[0]; values[AliDielectronVarManager::kV0ACyH2] = qvec[1]; 
  values[AliDielectronVarManager::kV0ACrpH2] = qvec[2];
  // VZERO event plane resolution
  values[AliDielectronVarManager::kV0ArpResH2] = 1.0;
  values[AliDielectronVarManager::kV0CrpResH2] = 1.0;
  values[AliDielectronVarManager::kV0ACrpResH2] = 1.0;
  // Q vector components correlations  
  values[AliDielectronVarManager::kV0XaXcH2] = values[AliDielectronVarManager::kV0AxH2]*values[AliDielectronVarManager::kV0CxH2];
  values[AliDielectronVarManager::kV0XaYaH2] = values[AliDielectronVarManager::kV0AxH2]*values[AliDielectronVarManager::kV0AyH2];
  values[AliDielectronVarManager::kV0XaYcH2] = values[AliDielectronVarManager::kV0AxH2]*values[AliDielectronVarManager::kV0CyH2];
  values[AliDielectronVarManager::kV0YaXcH2] = values[AliDielectronVarManager::kV0AyH2]*values[AliDielectronVarManager::kV0CxH2];
  values[AliDielectronVarManager::kV0YaYcH2] = values[AliDielectronVarManager::kV0AyH2]*values[AliDielectronVarManager::kV0CyH2];
  values[AliDielectronVarManager::kV0XcYcH2] = values[AliDielectronVarManager::kV0CxH2]*values[AliDielectronVarManager::kV0CyH2];


  // event plane differences used for EP resolution calculation
  values[AliDielectronVarManager::kV0ATPCDiffH2]   = TMath::Cos( 2.*(values[AliDielectronVarManager::kV0ArpH2] - 
								     values[AliDielectronVarManager::kTPCrpH2]) ); 
  
  values[AliDielectronVarManager::kV0CTPCDiffH2]   = TMath::Cos( 2.*(values[AliDielectronVarManager::kV0CrpH2] - 
								     values[AliDielectronVarManager::kTPCrpH2]) ); 
  
  values[AliDielectronVarManager::kV0AV0CDiffH2]   = TMath::Cos( 2.*(values[AliDielectronVarManager::kV0ArpH2] - 
								     values[AliDielectronVarManager::kV0CrpH2]) ); 

  values[AliDielectronVarManager::kv0ATPCDiffH2]   = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0ArpH2] - 
								     values[AliDielectronVarManager::kTPCrpH2]) ); 
  
  values[AliDielectronVarManager::kv0CTPCDiffH2]   = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0CrpH2] - 
								     values[AliDielectronVarManager::kTPCrpH2]) ); 
  
  values[AliDielectronVarManager::kv0Av0CDiffH2]   = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0ArpH2] - 
								     values[AliDielectronVarManager::kv0CrpH2]) ); 

  values[AliDielectronVarManager::kv0Av0C0DiffH2]  = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0ArpH2] - 
								     values[AliDielectronVarManager::kv0C0rpH2]) ); 

  values[AliDielectronVarManager::kv0Av0C3DiffH2]  = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0ArpH2] - 
								     values[AliDielectronVarManager::kv0C3rpH2]) ); 

  values[AliDielectronVarManager::kv0Cv0A0DiffH2]  = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0CrpH2] - 
								     values[AliDielectronVarManager::kv0A0rpH2]) ); 

  values[AliDielectronVarManager::kv0Cv0A3DiffH2]  = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0CrpH2] - 
								     values[AliDielectronVarManager::kv0A3rpH2]) ); 

  values[AliDielectronVarManager::kv0A0v0A3DiffH2] = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0A0rpH2] - 
								     values[AliDielectronVarManager::kv0A3rpH2]) ); 

  values[AliDielectronVarManager::kv0C0v0C3DiffH2] = TMath::Cos( 2.*(values[AliDielectronVarManager::kv0C0rpH2] - 
								     values[AliDielectronVarManager::kv0C3rpH2]) ); 


}

inline void AliDielectronVarManager::FillVarESDEvent(const AliESDEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  // 
  
  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  Double_t centralityF=-1; Double_t centralitySPD=-1;
  AliCentrality *esdCentrality = const_cast<AliESDEvent*>(event)->GetCentrality();
  if (esdCentrality) centralityF = esdCentrality->GetCentralityPercentile("V0M");
  if (esdCentrality) centralitySPD = esdCentrality->GetCentralityPercentile("CL1");
  
  // Fill AliESDEvent interface specific information
  const AliESDVertex *primVtx = event->GetPrimaryVertex();
  values[AliDielectronVarManager::kXRes]       = primVtx->GetXRes();
  values[AliDielectronVarManager::kYRes]       = primVtx->GetYRes();
  values[AliDielectronVarManager::kZRes]       = primVtx->GetZRes();
  values[AliDielectronVarManager::kCentrality] = centralityF;
  values[AliDielectronVarManager::kCentralitySPD] = centralitySPD;

  const AliESDVertex *vtxTPC = event->GetPrimaryVertexTPC(); 
  values[AliDielectronVarManager::kNVtxContribTPC] = (vtxTPC ? vtxTPC->GetNContributors() : 0);

  // Event multiplicity estimators
  Int_t nTrSPD05=0; Int_t nTrITSTPC05=0; Int_t nTrITSSA05=0;
  event->EstimateMultiplicity(nTrSPD05, nTrITSTPC05, nTrITSSA05, 0.5);
  values[AliDielectronVarManager::kNaccTrckltsEsd05] = nTrSPD05;
  values[AliDielectronVarManager::kNaccItsTpcEsd05] = nTrITSTPC05;
  values[AliDielectronVarManager::kNaccItsPureEsd05] = nTrITSSA05;
  values[AliDielectronVarManager::kNaccTrckltsEsd05Corr] = 
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrSPD05),values[AliDielectronVarManager::kZvPrim],0);
  values[AliDielectronVarManager::kNaccItsTpcEsd05Corr] = 
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrITSTPC05),values[AliDielectronVarManager::kZvPrim],3);
  values[AliDielectronVarManager::kNaccItsPureEsd05Corr] = 
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrITSSA05),values[AliDielectronVarManager::kZvPrim],6);
  
  Int_t nTrSPD10=0; Int_t nTrITSTPC10=0; Int_t nTrITSSA10=0;
  event->EstimateMultiplicity(nTrSPD10, nTrITSTPC10, nTrITSSA10, 1.0);
  values[AliDielectronVarManager::kNaccTrckltsEsd10] = nTrSPD10;
  values[AliDielectronVarManager::kNaccItsTpcEsd10] = nTrITSTPC10;
  values[AliDielectronVarManager::kNaccItsPureEsd10] = nTrITSSA10;
  values[AliDielectronVarManager::kNaccTrckltsEsd10Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrSPD10),values[AliDielectronVarManager::kZvPrim],1);
  values[AliDielectronVarManager::kNaccItsTpcEsd10Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrITSTPC10),values[AliDielectronVarManager::kZvPrim],4);
  values[AliDielectronVarManager::kNaccItsPureEsd10Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrITSSA10),values[AliDielectronVarManager::kZvPrim],7); 

  Int_t nTrSPD16=0; Int_t nTrITSTPC16=0; Int_t nTrITSSA16=0;
  event->EstimateMultiplicity(nTrSPD16, nTrITSTPC16, nTrITSSA16, 1.6);
  values[AliDielectronVarManager::kNaccTrckltsEsd16] = nTrSPD16;
  values[AliDielectronVarManager::kNaccItsTpcEsd16] = nTrITSTPC16;
  values[AliDielectronVarManager::kNaccItsPureEsd16] = nTrITSSA16;
  values[AliDielectronVarManager::kNaccTrckltsEsd16Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrSPD16),values[AliDielectronVarManager::kZvPrim],2);
  values[AliDielectronVarManager::kNaccItsTpcEsd16Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrITSTPC16),values[AliDielectronVarManager::kZvPrim],5);
  values[AliDielectronVarManager::kNaccItsPureEsd16Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event,Double_t(nTrITSSA16),values[AliDielectronVarManager::kZvPrim],8);
 
}

inline void AliDielectronVarManager::FillVarAODEvent(const AliAODEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //

  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  // Fill AliAODEvent interface specific information
  AliAODHeader *header = event->GetHeader();

  Double_t centralityF=-1; Double_t centralitySPD=-1;
  AliCentrality *aodCentrality = header->GetCentralityP();
  if (aodCentrality) centralityF = aodCentrality->GetCentralityPercentile("V0M");
  if (aodCentrality) centralitySPD = aodCentrality->GetCentralityPercentile("CL1");
  values[AliDielectronVarManager::kCentrality] = centralityF;
  values[AliDielectronVarManager::kCentralitySPD] = centralitySPD;

  values[AliDielectronVarManager::kRefMult]        = header->GetRefMultiplicity();        // similar to Ntrk
  values[AliDielectronVarManager::kRefMultTPConly] = header->GetTPConlyRefMultiplicity(); // similar to Nacc

  // nanoAODs (w/o AliCentrality branch) should have the VOM centrality stored in the header
  if(!header->GetCentralityP())
    values[AliDielectronVarManager::kCentrality] = header->GetCentrality();
  // nanoAODs (w/o AliEventPlane branch) should have the tpc event plane angle stored in the header
  if(!header->GetEventplaneP()) {

    //    values[AliDielectronVarManager::kNTrk] = header->GetRefMultiplicity();    // overwritten datamembers in "our" nanoAODs
    //    values[AliDielectronVarManager::kNacc] = header->GetRefMultiplicityPos(); // overwritten datamembers in "our" nanoAODs

    TVector2 qvec;
    // TPC
    qvec.Set(header->GetEventplaneQx(), header->GetEventplaneQy());
    values[AliDielectronVarManager::kTPCxH2uc]   = qvec.X();
    values[AliDielectronVarManager::kTPCyH2uc]   = qvec.Y();
    values[AliDielectronVarManager::kTPCmagH2uc] = qvec.Mod();
    values[AliDielectronVarManager::kTPCrpH2uc]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);

    // VZERO
    AliEventplane ep2;
    // get event plane corrections from the VZERO EP selection task
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliVZEROEPSelectionTask *eptask = dynamic_cast<AliVZEROEPSelectionTask *>(man->GetTask("AliVZEROEPSelectionTask"));
    if(eptask) eptask->SetEventplaneParams(&ep2,centralityF);
    else printf("no VZERO event plane selection task added! \n");

    Double_t qx = 0, qy = 0;
    ep2.CalculateVZEROEventPlane(event,10, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ACrpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    values[AliDielectronVarManager::kv0ACxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0ACyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0ACmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 8, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ArpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    values[AliDielectronVarManager::kv0AxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0AyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0AmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 9, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0CrpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    values[AliDielectronVarManager::kv0CxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0CyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0CmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 0, 0, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C0rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    ep2.CalculateVZEROEventPlane(event, 3, 3, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C3rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    ep2.CalculateVZEROEventPlane(event, 4, 4, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A0rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
    ep2.CalculateVZEROEventPlane(event, 7, 7, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A3rpH2]  = ((TMath::Abs(qvec.X())>1.0e-10) ? TMath::ATan2(qvec.Y(),qvec.X())/2.0 : 0.0);
  }

  const AliAODVertex *vtxtpc = GetVertex(event, AliAODVertex::kMainTPC);
  values[AliDielectronVarManager::kNVtxContribTPC] = (vtxtpc ? vtxtpc->GetNContributors() : 0);

}
  
inline void AliDielectronVarManager::FillVarMCEvent(const AliMCEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //

  // Fill common AliVEvent interface information
  //  FillVarVEvent(event, values);
  const AliVVertex* vtx = event->GetPrimaryVertex();
  values[AliDielectronVarManager::kXvPrim]       = (vtx ? vtx->GetX() : 0.0);
  values[AliDielectronVarManager::kYvPrim]       = (vtx ? vtx->GetY() : 0.0);
  values[AliDielectronVarManager::kZvPrim]       = (vtx ? vtx->GetZ() : 0.0);
  // Fill AliMCEvent interface specific information
  values[AliDielectronVarManager::kNch]   = AliDielectronHelper::GetNch(event, 1.6);
  values[AliDielectronVarManager::kNch05] = AliDielectronHelper::GetNch(event, 0.5);
  values[AliDielectronVarManager::kNch10] = AliDielectronHelper::GetNch(event, 1.0);
  
  values[AliDielectronVarManager::kNumberOfJPsis] = AliDielectronHelper::GetNMothers(event, 0.9, 443, 11);
  values[AliDielectronVarManager::kNumberOfJPsisPrompt]  = AliDielectronHelper::GetNMothers(event, 0.9, 443, 11, 1);
  values[AliDielectronVarManager::kNumberOfJPsisNPrompt] = AliDielectronHelper::GetNMothers(event, 0.9, 443, 11, 0);
}

inline void AliDielectronVarManager::FillVarTPCEventPlane(const AliEventplane *evplane, Double_t * const values)
{
  //
  // Fill TPC event plane information after correction
  //
  if(evplane) {
    TVector2 *qcorr  = const_cast<AliEventplane *>(evplane)->GetQVector();  // This is the "corrected" Q-Vector
    TVector2 *qcsub1 = const_cast<AliEventplane *>(evplane)->GetQsub1();
    TVector2 *qcsub2 = const_cast<AliEventplane *>(evplane)->GetQsub2();
    if(qcorr) {
      values[AliDielectronVarManager::kTPCxH2]   = qcorr->X();
      values[AliDielectronVarManager::kTPCyH2]   = qcorr->Y();
      values[AliDielectronVarManager::kTPCmagH2] = qcorr->Mod();
      values[AliDielectronVarManager::kTPCrpH2]  = ((TMath::Abs(qcorr->X())>1.0e-10) ? TMath::ATan2(qcorr->Y(),qcorr->X())/2.0 : 0.0);
    }
    if(qcsub1 && qcsub2) {
      values[AliDielectronVarManager::kTPCsub1xH2]   = qcsub1->X();
      values[AliDielectronVarManager::kTPCsub1yH2]   = qcsub1->Y();
      values[AliDielectronVarManager::kTPCsub1rpH2]  = ((TMath::Abs(qcsub1->X())>1.0e-10) ? TMath::ATan2(qcsub1->Y(),qcsub1->X())/2.0 : 0.0);

      values[AliDielectronVarManager::kTPCsub2xH2]   = qcsub2->X();
      values[AliDielectronVarManager::kTPCsub2yH2]   = qcsub2->Y();
      values[AliDielectronVarManager::kTPCsub2rpH2]  = ((TMath::Abs(qcsub2->X())>1.0e-10) ? TMath::ATan2(qcsub2->Y(),qcsub2->X())/2.0 : 0.0);

      values[AliDielectronVarManager::kTPCsub12DiffH2] = TMath::Cos( 2.*(values[AliDielectronVarManager::kTPCsub1rpH2] -
									 values[AliDielectronVarManager::kTPCsub2rpH2]) );
      values[AliDielectronVarManager::kTPCsub12DiffH2Sin] = TMath::Sin( 2.*(values[AliDielectronVarManager::kTPCsub1rpH2] -
									    values[AliDielectronVarManager::kTPCsub2rpH2]) );
    }
  }
}

inline void AliDielectronVarManager::InitESDpid(Int_t type)
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

inline void AliDielectronVarManager::InitAODpidUtil(Int_t type)
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


inline void AliDielectronVarManager::InitEstimatorAvg(const Char_t* filename)
{
  //
  // initialize the profile histograms neccessary for the correction of the multiplicity estimators in pp collisions
  //
  
  const Char_t* estimatorNames[9] = {"SPDmult05","SPDmult10","SPDmult16",
				     "ITSTPC05", "ITSTPC10", "ITSTPC16", 
				     "ITSSA05",  "ITSSA10",  "ITSSA16"};
  const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
  TFile* file=TFile::Open(filename);
  if(!file) return;
  
  for(Int_t ip=0; ip<4; ++ip) {
    for(Int_t ie=0; ie<9; ++ie) {
      fgMultEstimatorAvg[ip][ie] = (TProfile*)(file->Get(Form("%s_%s",estimatorNames[ie],periodNames[ip]))->Clone(Form("%s_%s_clone",estimatorNames[ie],periodNames[ip])));
    }
  }
}


inline void AliDielectronVarManager::InitTRDpidEffHistograms(const Char_t* filename)
{
  //
  // initialize the 3D histograms with the TRD pid efficiency histograms
  //
  
  // reset the centrality ranges and the efficiency histograms
  for(Int_t i=0; i<10; ++i) {         // centrality ranges
    for(Int_t j=0; j<4; ++j) fgTRDpidEffCentRanges[i][j] = -1.;
    if(fgTRDpidEff[i][0]) {
      delete fgTRDpidEff[i][0];
      fgTRDpidEff[i][0] = 0x0;
    }
    if(fgTRDpidEff[i][1]) {
      delete fgTRDpidEff[i][1];
      fgTRDpidEff[i][1] = 0x0;
    }
  }
  
  TFile* file=TFile::Open(filename);
  TList* keys=file->GetListOfKeys();
  Int_t idxp=0; Int_t idxn=0;
  for(Int_t i=0; i<keys->GetEntries(); ++i) {
    if(idxp>=10) continue;
    if(idxn>=10) continue;
    TString name=((TKey*)keys->At(i))->ReadObj()->GetName();    
    // Name of histograms should be in the format:
    // TRDeff<field>_cent_<centLow>_<centHigh>
    // <field> is either "BPLUS" or "BMINUS"
    if(!(name.Contains("BPLUS") || name.Contains("BMINUS"))) continue;
    TObjArray* arr = name.Tokenize("_");
    Bool_t isBplus = kTRUE;
    if(name.Contains("BMINUS")) isBplus = kFALSE;
    TString centMinStr = arr->At(2)->GetName();
    TString centMaxStr = arr->At(3)->GetName();
    delete arr;
    if(isBplus) {
      fgTRDpidEffCentRanges[idxp][2] = centMinStr.Atof();
      fgTRDpidEffCentRanges[idxp][3] = centMaxStr.Atof();
      fgTRDpidEff[idxp][1] = (TH3D*)(file->Get(name.Data())->Clone(Form("%s_clone",name.Data())));
      ++idxp;
    }
    else {
      fgTRDpidEffCentRanges[idxn][0] = centMinStr.Atof();
      fgTRDpidEffCentRanges[idxn][1] = centMaxStr.Atof();
      fgTRDpidEff[idxn][0] = (TH3D*)(file->Get(name.Data())->Clone(Form("%s_clone",name.Data())));
      ++idxn;
    }
  }
}


inline void AliDielectronVarManager::InitVZEROCalibrationHistograms(Int_t runNo) {
  //
  // Initialize the VZERO channel-by-channel calibration histograms
  //

  //initialize only once
  if(fgVZEROCalib[0]) return;
  
  for(Int_t i=0; i<64; ++i) 
    if(fgVZEROCalib[i]) {
      delete fgVZEROCalib[i];
      fgVZEROCalib[i] = 0x0;
    }
  
  TFile file(fgVZEROCalibrationFile.Data());
  
  for(Int_t i=0; i<64; ++i){
    fgVZEROCalib[i] = (TProfile2D*)(file.Get(Form("RUN%d_ch%d_VtxCent", runNo, i)));
    if (fgVZEROCalib[i]) fgVZEROCalib[i]->SetDirectory(0x0);
  }
}


inline void AliDielectronVarManager::InitVZERORecenteringHistograms(Int_t runNo) {
  //
  // Initialize the VZERO event plane recentering histograms
  //

  //initialize only once
  if(fgVZERORecentering[0][0]) return;
  
  for(Int_t i=0; i<2; ++i)
    for(Int_t j=0; j<2; ++j)
      if(fgVZERORecentering[i][j]) {
        delete fgVZERORecentering[i][j];
        fgVZERORecentering[i][j] = 0x0;
      }
  
  TFile file(fgVZERORecenteringFile.Data());
  if (!file.IsOpen()) return;
  
  fgVZERORecentering[0][0] = (TProfile2D*)(file.Get(Form("RUN%d_QxA_CentVtx", runNo)));
  fgVZERORecentering[0][1] = (TProfile2D*)(file.Get(Form("RUN%d_QyA_CentVtx", runNo)));
  fgVZERORecentering[1][0] = (TProfile2D*)(file.Get(Form("RUN%d_QxC_CentVtx", runNo)));
  fgVZERORecentering[1][1] = (TProfile2D*)(file.Get(Form("RUN%d_QyC_CentVtx", runNo)));

  if (fgVZERORecentering[0][0]) fgVZERORecentering[0][0]->SetDirectory(0x0);
  if (fgVZERORecentering[0][1]) fgVZERORecentering[0][1]->SetDirectory(0x0);
  if (fgVZERORecentering[1][0]) fgVZERORecentering[1][0]->SetDirectory(0x0);
  if (fgVZERORecentering[1][1]) fgVZERORecentering[1][1]->SetDirectory(0x0);
  
}


inline Double_t AliDielectronVarManager::GetTRDpidEfficiency(Int_t runNo, Double_t centrality, 
 				                             Double_t eta, Double_t trdPhi, Double_t pout,
 				                             Double_t& effErr) {
  //
  // return the efficiency in the given phase space cell
  //
  // LHC10h data----------------------------------------------
  Bool_t isBplus = kTRUE;
  if(runNo<=138275) isBplus = kFALSE;
  // TODO: check magnetic polarity for runs in 2011 data
  // ---------------------------------------------------------
  Int_t centIdx = -1;
  for(Int_t icent=0; icent<10; ++icent) {
    if(isBplus) {
      if(centrality>=fgTRDpidEffCentRanges[icent][2] && centrality<fgTRDpidEffCentRanges[icent][3]) {
 	centIdx = icent;
 	break;
      }
    }
    else {
      if(centrality>=fgTRDpidEffCentRanges[icent][0] && centrality<fgTRDpidEffCentRanges[icent][1]) {
 	centIdx = icent;
 	break;
      }
    }
  }
  //TODO: chek logick
  if (centIdx<0) return 1;
  
  TH3D* effH = fgTRDpidEff[centIdx][(isBplus ? 1 : 0)];
  if(!effH) {effErr=0x0; return 1.0;}
  Int_t etaBin = effH->GetXaxis()->FindBin(eta);
  if(eta<effH->GetXaxis()->GetXmin()) etaBin=1;
  if(eta>effH->GetXaxis()->GetXmax()) etaBin=effH->GetXaxis()->GetNbins();
  Int_t phiBin = effH->GetYaxis()->FindBin(trdPhi);
  if(trdPhi<effH->GetYaxis()->GetXmin()) phiBin=1;
  if(trdPhi>effH->GetYaxis()->GetXmax()) phiBin=effH->GetYaxis()->GetNbins();
  Int_t poutBin = effH->GetZaxis()->FindBin(pout);
  if(pout<effH->GetZaxis()->GetXmin()) poutBin=1;
  if(pout>effH->GetZaxis()->GetXmax()) poutBin=effH->GetZaxis()->GetNbins();
  Double_t eff = effH->GetBinContent(etaBin, phiBin, poutBin);
  effErr = effH->GetBinError(etaBin, phiBin, poutBin);
  if(eff<-0.0001) {
    effErr = 0.0;
    eff = 1.0;
  }
  return eff;
}


inline void AliDielectronVarManager::SetEvent(AliVEvent * const ev)
{
  
  fgEvent = ev;
  if (fgKFVertex) delete fgKFVertex;
  fgKFVertex=0x0;
  if (!ev) return;
  if (ev->GetPrimaryVertex()) fgKFVertex=new AliKFVertex(*ev->GetPrimaryVertex());

  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues;++i) fgData[i]=0.;
  AliDielectronVarManager::Fill(fgEvent, fgData);
}

inline void AliDielectronVarManager::SetEventData(const Double_t data[AliDielectronVarManager::kNMaxValues])
{
  for (Int_t i=0; i<kNMaxValues;++i) fgData[i]=0.;
  for (Int_t i=kPairMax; i<kNMaxValues;++i) fgData[i]=data[i];
}


inline Bool_t AliDielectronVarManager::GetDCA(const AliAODTrack *track, Double_t d0z0[2])
{
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    return kTRUE;
  }
  
  Bool_t ok=kFALSE;
  if(fgEvent) {
    Double_t covd0z0[3];
    //AliAODTrack copy(*track);
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);

    Float_t xstart = etp.GetX();
    if(xstart>3.) {
    d0z0[0]=-999.;
    d0z0[1]=-999.;
    //printf("This method can be used only for propagation inside the beam pipe \n");
    return kFALSE;
    }


    AliAODVertex *vtx =(AliAODVertex*)(fgEvent->GetPrimaryVertex());
    Double_t fBzkG = fgEvent->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
    //ok = copy.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}

inline void AliDielectronVarManager::SetTPCEventPlane(AliEventplane *const evplane)
{
  
  fgTPCEventPlane = evplane;
  FillVarTPCEventPlane(evplane,fgData);
  //  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues;++i) fgData[i]=0.;
  //  AliDielectronVarManager::Fill(fgEvent, fgData);
}


//_________________________________________________________________
inline void AliDielectronVarManager::GetVzeroRP(const AliVEvent* event, Double_t* qvec, Int_t sideOption) {
  //
  // Get the reaction plane from the VZERO detector for a given harmonic
  //
  // sideOption = 0- V0A, 1- V0C, 2-both
  //  Q{x,y} = SUM_i mult(i) * {cos(n*phi_i), sin(n*phi_i)} 
  //  phi_i - phi angle of the VZERO sector i
  //          Each sector covers 45 degrees(8 sectors per ring). Middle of sector 0 is at 45/2
  //        channel 0: 22.5
  //                1: 22.5+45
  //                2: 22.5+45*2
  //               ...
  //        at the next ring continues the same
  //        channel 8: 22.5
  //        channel 9: 22.5 + 45
  //               ... 
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  Int_t phi;
  Float_t mult;
  
  // get centrality and vertex for this event
  Double_t centralitySPD = -1; Double_t vtxZ = -999.;
  if(event->IsA() == AliESDEvent::Class()) {
    const AliESDEvent* esdEv = static_cast<const AliESDEvent*>(event);
    AliCentrality *esdCentrality = const_cast<AliESDEvent*>(esdEv)->GetCentrality();
    if(esdCentrality) centralitySPD = esdCentrality->GetCentralityPercentile("CL1");
  }
  if(event->IsA() == AliAODEvent::Class()) {
    const AliAODEvent* aodEv = static_cast<const AliAODEvent*>(event);
    AliAODHeader *header = aodEv->GetHeader();
    AliCentrality *aodCentrality = header->GetCentralityP();
    if(aodCentrality) centralitySPD = aodCentrality->GetCentralityPercentile("CL1");
  }
  const AliVVertex *primVtx = event->GetPrimaryVertex();
  if(!primVtx) return;
  vtxZ = primVtx->GetZ();
  if(TMath::Abs(vtxZ)>10.) return;
  if(centralitySPD<0. || centralitySPD>80.) return;
  
  Int_t binCent = -1; Int_t binVtx = -1;
  if(fgVZEROCalib[0]) {
    binVtx = fgVZEROCalib[0]->GetXaxis()->FindBin(vtxZ);
    binCent = fgVZEROCalib[0]->GetYaxis()->FindBin(centralitySPD);
  }
  
  AliVVZERO* vzero = event->GetVZEROData();
  Double_t average = 0.0;
  for(Int_t iChannel=0; iChannel<64; ++iChannel) {
    if(iChannel<32 && sideOption==0) continue;
    if(iChannel>=32 && sideOption==1) continue;
    phi=iChannel%8;
    mult = vzero->GetMultiplicity(iChannel);
    if(fgVZEROCalib[iChannel])
      average = fgVZEROCalib[iChannel]->GetBinContent(binVtx, binCent);
    if(average>1.0e-10 && mult>0.5) 
      mult /= average;
    else
      mult = 0.0;
    //  2nd harmonic
    qvec[0] += mult*(2.0*TMath::Power(kX[phi],2.0)-1);
    qvec[1] += mult*(2.0*kX[phi]*kY[phi]);
  }    // end loop over channels 
  
  // do recentering
  if(fgVZERORecentering[0][0]) {
//     printf("vzero: %p\n",fgVZERORecentering[0][0]);
    Int_t binCentRecenter = -1; Int_t binVtxRecenter = -1;
    binCentRecenter = fgVZERORecentering[0][0]->GetXaxis()->FindBin(centralitySPD);
    binVtxRecenter = fgVZERORecentering[0][0]->GetYaxis()->FindBin(vtxZ);
    if(sideOption==0) {  // side A
      qvec[0] -= fgVZERORecentering[0][0]->GetBinContent(binCentRecenter, binVtxRecenter);
      qvec[1] -= fgVZERORecentering[0][1]->GetBinContent(binCentRecenter, binVtxRecenter);
    }
    if(sideOption==1) {  // side C
      qvec[0] -= fgVZERORecentering[1][0]->GetBinContent(binCentRecenter, binVtxRecenter);
      qvec[1] -= fgVZERORecentering[1][1]->GetBinContent(binCentRecenter, binVtxRecenter);
    }
    if(sideOption==2) {  // side A and C together
      qvec[0] -= fgVZERORecentering[0][0]->GetBinContent(binCentRecenter, binVtxRecenter);
      qvec[0] -= fgVZERORecentering[1][0]->GetBinContent(binCentRecenter, binVtxRecenter);
      qvec[1] -= fgVZERORecentering[0][1]->GetBinContent(binCentRecenter, binVtxRecenter);
      qvec[1] -= fgVZERORecentering[1][1]->GetBinContent(binCentRecenter, binVtxRecenter);
    }
  }
  
  // calculate the reaction plane
  if(TMath::Abs(qvec[0])>1.0e-10)
    qvec[2] = TMath::ATan2(qvec[1],qvec[0])/2.0;
}

//______________________________________________________________________________                                                                                                                                                                                     
inline AliAODVertex* AliDielectronVarManager::GetVertex(const AliAODEvent* event, AliAODVertex::AODVtx_t vtype) {
  // Get vertex
  Int_t nVertices=event->GetNumberOfVertices();
  for(Int_t iVert=0; iVert<nVertices; iVert++){
    AliAODVertex *v=event->GetVertex(iVert);
    //    printf(" vtx %d  contrib %d  daughters %d \n ",v->GetType(),v->GetNContributors(), v->GetNDaughters());
    if(v->GetType()==vtype) return v;
  }
  return 0;
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

