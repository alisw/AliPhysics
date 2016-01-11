#ifndef ALICMEVARMANAGER_H
#define ALICMEVARMANAGER_H


/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


// This class is adopted from AliDiElectronVarManager
// main changes:
// dropped: dependencies on other dielectron classes, mixed event handling, pair reconstruction, mc
// added  : few variables, track status flags

//#############################################################
//#                                                           #
//#         Class AliDiElectronVarManager                     #
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
#include <TProfile3D.h>
#include <TH3D.h>
#include <THn.h>
#include <THnBase.h>
#include <TSpline.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TKey.h>
#include <TBits.h>
#include <TRandom3.h>

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
#include <AliMultSelection.h>
#include <AliAODpidUtil.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliESDtrackCuts.h>

//#include "AliDielectronMC.h"
//#include "AliDielectronPID.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVZEROEPSelectionTask.h"

#include "AliAODMCHeader.h"
#include "assert.h"

class AliVEvent;

//________________________________________________________________
class AliCMEVarManager : public TNamed {

public:


  // tracking flags as in AliESDtrack.h (or AliVTrack?)
  // NOTE: when in doubt check consistency with aliroot
  enum TrackingStatus {
    kITSin=0,
    kITSout,
    kITSrefit,
    kITSpid,
    kTPCin,
    kTPCout,
    kTPCrefit,
    kTPCpid,
    kTRDin,
    kTRDout,
    kTRDrefit,
    kTRDpid,
    kTOFin,
    kTOFout,
    kTOFrefit,
    kTOFpid,
    kTOFmismatch,
    kHMPIDout,
    kHMPIDpid,
    kEMCALmatch,
    kPHOSmatch,
    kTRDbackup,
    kTRDStop,
    kESDpid,
    kTIME,
    kGlobalMerge,
    kITSpureSA,
    kMultInV0,
    kMultSec,
    kTRDnPlanes,
    kEMCALNoMatch,
    kNTrackingStatus
  };


  // to be moved to AliReducedTrackInfo
  enum TrackingFlags {
    kTpcEP=0,
    kGammaConv,
    kK0s,
    kLambda,
    kALambda,
    kKink0,
    kKink1,
    kKink3,
    kPureGammaConv,
    kPureK0s,
    kPureLambda,
    kPureALambda,
    kNTrackingFlags
  };



  // Particle specific variables
  enum ValueTypes {
    kPx = 0,                 // px
    kPy,                     // py
    kPz,                     // pz
    kPt,                     // transverse momentum
    kPtSq,                   // transverse momentum squared
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
    kTPCclsIRO,             // TPC clusters inner read out
    kTPCclsORO,             // TPC clusters outer read out
    kTrackStatus,            // track status bits
    kTrackingStatus, // kTrackStatus above is not good because of float precision used for the arrays, now we read all bits individually
    kFilterBit=kTrackingStatus+kNTrackingStatus,      // AOD filter bits

    kNclsTRD,                // number of clusters assigned in the TRD
    kTRDntracklets,          // number of TRD tracklets used for tracking/PID TODO: correct getter
    kTRDpidQuality,          // number of TRD tracklets used for PID
    kTRDchi2,                // chi2 in TRD
    kTRDchi2Trklt,	     // chi2/trklts in TRD
    kTRDprobEle,             // TRD electron pid probability
    kTRDprobPio,             // TRD electron pid probability
    kTRDprob2DEle,           // TRD electron pid probability 2D LQ
    kTRDprob2DPio,           // TRD electron pid probability 2D LQ
    kTRDphi,                 // Phi angle of the track at the entrance of the TRD
    kTRDpidEffLeg,           // TRD pid efficiency from conversion electrons
    kTRDsignal,              // TRD signal

    kImpactParXY,            // Impact parameter in XY plane
    kImpactParZ,             // Impact parameter in Z
    kImpactPar2D,            // Impact parameter in 2d set in analysis class
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

    kITSnSigmaEleRaw,        // raw number of sigmas to the dE/dx electron line in the ITS
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

    kTPCnSigmaEleRaw,        // raw number of sigmas to the dE/dx electron line in the TPC
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

    kLegEff,                 // single electron efficiency
    kOneOverLegEff,          // 1 / single electron efficiency (correction factor)
    kV0Index0,               // v0 index 0
    kKinkIndex0,             // kink index 0

    //TRD online trackinfo
    kTRDonlineA,             //transverse offset from nominal primary vertex
    kTRDonlineLayerMask,     //tracklet map for given TRD online track
    kTRDonlineFirstLayer,    //First layer having a tracklet associated to the track
    kTRDonlinePID,           //TRD online track PID (number between 0-256, 2013 trigger threshols SE: 144, QU: 164)
    kTRDonlinePt,            //TRD online track pt estimate (via ->Pt())
    kTRDonlineStack,         //TRD online track Stack number (0-89), unique due to stack-wise tracking
    kTRDonlineSector,        //TRD online track Sector number
    kTRDonlineTrackInTime,   //TRD online track boolean, if in time window
    kTRDonlineFlagsTiming,   //TRD online timing flags
    kTRDonlineLabel,         //TRD online track label
    kTRDonlineNTracklets,     //TRD online: number of contributing online tracklets


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

    kPairPlaneAngle1A,         // angle between ee decay plane and x'-z reaction plane by using V0-A
    kPairPlaneAngle2A,         // angle between ee decay plane and (p1+p2) rot ez
    kPairPlaneAngle3A,         // angle between ee decay plane and (p1+p2) rot (p1+p2)x'z
    kPairPlaneAngle4A,         // angle between ee decay plane and x'-y' plane
    kPairPlaneAngle1C,         // using v0-C
    kPairPlaneAngle2C,
    kPairPlaneAngle3C,
    kPairPlaneAngle4C,
    kPairPlaneAngle1AC,        // using v0-AC
    kPairPlaneAngle2AC,
    kPairPlaneAngle3AC,
    kPairPlaneAngle4AC,
    kPairPlaneAngle1Ran,       // using random reaction plane
    kPairPlaneAngle2Ran,
    kPairPlaneAngle3Ran,
    kPairPlaneAngle4Ran,
    kRandomRP,                //Random reaction plane
    kDeltaPhiRandomRP,        //delta phi of the pair

    kPairPlaneMagInPro,     // Inner Product of strong magnetic field and ee plane
	kCos2PhiCS,              // Cosine of 2*phi in mother's rest frame in the Collins-Soper picture
    kCosTilPhiCS,            // Shifted phi depending on kThetaCS
    kCosPhiH2,               // cosine of pair phi for 2nd harmonic
    kSinPhiH2,               // sinus  of pair phi for 2nd harmonic
    kDeltaPhiV0ArpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-A
    kDeltaPhiV0CrpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-C
    kDeltaPhiV0ACrpH2,       // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-A + V0-C
    kV0ArpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-A
    kV0CrpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-C
    kV0ACrpH2FlowV2,         // v2 coefficient with respect to the 2nd order reaction plane from V0-A + V0-C
    kDeltaPhiv0ArpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-A (EPtask)
    kDeltaPhiv0CrpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-C
    kDeltaPhiv0ACrpH2,       // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from V0-AC
    kDeltaPhiTPCrpH2,        // Delta phi of the pair with respect to the 2nd order harmonic reaction plane from TPC
    kv0ArpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-A (EPtask)
    kv0CrpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from V0-C
    kv0ACrpH2FlowV2,         // v2 coefficient with respect to the 2nd order reaction plane from V0-A + V0-C
    kTPCrpH2FlowV2,          // v2 coefficient with respect to the 2nd order reaction plane from TPC
    kTPCrpH2FlowV2Sin,       // sinus of v2 coefficient with respect to the 2nd order reaction plane from TPC

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
    kPairEff,                // pair efficiency
    kOneOverPairEff,         // 1 / pair efficiency (correction factor)
    kOneOverPairEffSq,        // 1 / pair efficiency squared (correction factor)
    kRndmPair,               // radomly created number (used to apply special signal reduction cuts)
    kRndm=kRndmPair,
    kLeg1DCAsigXY,            //DCA in sigma for first daughter of the pair in xy-plane
    kLeg1DCAabsXY,            //DCA in cm for first daughter of the pair in xy-plane
    kLeg1DCAresXY,            //resolution from kov matrix for first daughter of the pair in xy-plane
    // pair dinstance of closest approach (dca) variables
    kPairDCAsigXY,             // dca in xy-plane calculated in orders of sigma calculated as sqrt(dcaD1^2 + dcaD2^2)
    kPairDCAsigZ,              // dca in z-plane calculated in orders of sigma calculated as sqrt(dcaD1^2 + dcaD2^2)
    kPairDCAabsXY,             // dca in xy-plane in absolute values (cm) calculated as sqrt(dcaD1^2 + dcaD2^2)
    kPairDCAabsZ,              // dca in z-plane in absolute values (cm) calculated as sqrt(dcaD1^2 + dcaD2^2)
    kPairLinDCAsigXY,          // dca in xy-plane calculated in orders of sigma calculated as (dcaD1 + dcaD2)/2)
    kPairLinDCAsigZ,           // dca in z-plane calculated in orders of sigma calculated as (dcaD1 + dcaD2)/2)
    kPairLinDCAabsXY,          // dca in xy-plane in absolute values (cm) calculated as (dcaD1 + dcaD2)/2)
    kPairLinDCAabsZ,           // dca in z-plane in absolute values (cm) calculated as (dcaD1 + dcaD2)/2)
    kPairs,                  // number of Ev1PM pair candidates after all cuts
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

    //// v0 reaction plane quantities from AliEPSelectionTaks, angles interval [-pi/2,+pi/2]
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
    kEqMultV0A,              // equalized VZERO multiplicity
    kEqMultV0C,
    kEqMultV0,
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
    // TPC reaction plane quantities, angle interval [-pi/2,+pi/2]
    kTPCxH2,                  // TPC x-component of the Q vector for 2nd harmonic (corrected)
    kTPCyH2,                  // TPC y-component of the Q vector for 2nd harmonic (corrected)
    kTPCmagH2,                // TPC reaction plane the Q vectors magnitude for 2nd harmonic (corrected)
    kTPCrpH2,                 // TPC reaction plane angle of the Q vector for 2nd harmonic (corrected)
    kCosTPCrpH2,              // cosine of TPC reaction plane angle of the Q vector for 2nd harmonic (corrected)
    kSinTPCrpH2,              // sinus of TPC reaction plane angle of the Q vector for 2nd harmonic (corrected)
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

    //ZDC reaction plane(v1 plane) quantities

    kZDCArpH1,                  // ZDC-A reaction plane of the Q vector for 1st harmonic
    kZDCCrpH1,                  // ZDC-C reaction plane of the Q vector for 1st harmonic
    kZDCACrpH1,                  // ZDC-AC reaction plane of the Q vector for 1st harmonic
    kZDCrpResH1,                  // 1st harmonic reaction plane resolution for ZDC
    kv0ZDCrpRes,                //ZDC reaction plane for 1st harmonic and VZERO reaction plane for 2nd harmonic correlation


    kNTrk,                   // number of tracks (or tracklets) TODO: ambiguous
    kTracks,                 // track after all cuts
    kNVtxContrib,             // number of primary vertex contibutors
    kNVtxContribTPC,         // number of TPC vertex contibutors
    kNacc,                   // Number of accepted tracks
    kMatchEffITSTPC,         // ruff estimate on the ITS-TPC matching efficiceny
    kNaccTrcklts,            // number of accepted SPD tracklets in |eta|<1.6
    kNaccTrcklts10,          // number of accepted SPD tracklets in |eta|<1.
    kNaccTrcklts0916,        // number of accepted SPD tracklets in 0.9<|eta|<1.6
    kNaccTrckltsCorr,        // number of accepted SPD tracklets in |eta|<1.6 (corrected)
    kNaccTrcklts10Corr,      // number of accepted SPD tracklets in |eta|<1. (corrected)
    kNaccTrcklts0916Corr,    // number of accepted SPD tracklets in 0.9<|eta|<1.6 (corrected)

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

    kCentrality,             // event centrality fraction V0M
    kCentralityV0A,          // event centrality fraction V0A
    kCentralityV0C,          // event centrality fraction V0C
    kCentralityZNA,          // event centrality fraction ZNA
    kCentralitySPD,          // centrality using SPD (from second layer)
    kCentralityV0mSPD,          // difference centrality fraction V0M-SPD
    kTriggerInclONL,         // online trigger bits fired (inclusive)
    kTriggerInclOFF,         // offline trigger bits fired (inclusive)
    kTriggerExclOFF,         // offline only this trigger bit fired (exclusive)
    kNevents,                // event counter
    kRunNumber,              // run number
    kMixingBin,
    kNMaxValues
    // TODO: (for A+A) ZDCEnergy, impact parameter, Iflag??
  };

  AliCMEVarManager();
  AliCMEVarManager(const char* name, const char* title);
  virtual ~AliCMEVarManager();
  static void Fill(const TObject* particle, Float_t * const values);
  static void FillVarMCParticle2(const AliVParticle *p1, const AliVParticle *p2, Float_t * const values);
  static void FillVarVParticle(const AliVParticle *particle,         Float_t * const values);

  static void InitESDpid(Int_t type=0);
  static void InitAODpidUtil(Int_t type=0);
  static void InitEstimatorAvg(const Char_t* filename);
  static void InitEstimatorObjArrayAvg(const TObjArray* array);
  static void InitTRDpidEffHistograms(const Char_t* filename);
  static void SetLegEffMap( TObject *map) { fgLegEffMap=map; }
  static void SetPairEffMap(TObject *map) { fgPairEffMap=map; }
  static void SetFillMap(   TBits   *map) { fgFillMap=map; }
  static void SetVZEROCalibrationFile(const Char_t* filename) {fgVZEROCalibrationFile = filename;}

  static void SetVZERORecenteringFile(const Char_t* filename) {fgVZERORecenteringFile = filename;}
  static void SetZDCRecenteringFile(const Char_t* filename) {fgZDCRecenteringFile = filename;}
  static void SetPIDResponse(AliPIDResponse *pidResponse) {fgPIDResponse=pidResponse;}
  static AliPIDResponse* GetPIDResponse() { return fgPIDResponse; }
  static void SetEvent(AliVEvent * const ev);
  static void SetEventData(const Double_t data[AliCMEVarManager::kNMaxValues]);
  static Bool_t GetDCA(const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0=0);
  static void SetTPCEventPlane(AliEventplane *const evplane);
  static void GetVzeroRP(const AliVEvent* event, Double_t* qvec, Int_t sideOption);      // 0- V0A; 1- V0C; 2- V0A+V0C
  static void GetZDCRP(const AliVEvent* event, Double_t qvec[][2]);
  static AliAODVertex* GetVertex(const AliAODEvent *event, AliAODVertex::AODVtx_t vtype);
  static THnF* CreateTHnF( const Char_t* name, const Char_t* title, Int_t nDimensions,TAxis* binLimits);

  static TProfile* GetEstimatorHistogram(Int_t period, Int_t type) {return fgMultEstimatorAvg[period][type];}
  static Double_t GetTRDpidEfficiency(Int_t runNo, Double_t centrality, Double_t eta, Double_t trdPhi, Double_t pout, Double_t& effErr);
  static Double_t GetSingleLegEff(Float_t * const values);
  static Double_t GetPairEff(Float_t * const values);

  static const AliKFVertex* GetKFVertex() {return fgKFVertex;}

  static const char* GetValueName(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i][0]:""; }
  static const char* GetValueLabel(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i][1]:""; }
  static const char* GetValueUnit(Int_t i) { return (i>=0&&i<kNMaxValues)?fgkParticleNames[i][2]:""; }
  static UInt_t GetValueType(const char* valname);
  static const Float_t* GetData() {return fgData;}
  static AliVEvent* GetCurrentEvent() {return fgEvent;}

  static Double_t GetValue(ValueTypes var) {return fgData[var];}
  static void SetValue(ValueTypes var, Double_t val) { fgData[var]=val; }


  static const char* fgkParticleNames[kNMaxValues][3];  //variable names

  static Bool_t Req(ValueTypes var) { return (fgFillMap ? fgFillMap->TestBitNumber(var) : kTRUE); }
  static void FillVarESDtrack(const AliESDtrack *particle,           Float_t * const values);
  static void FillVarAODTrack(const AliAODTrack *particle,           Float_t * const values);
  static void FillVarVTrdTrack(const AliVParticle *particle,         Float_t * const values);
  static void FillVarMCParticle(const AliMCParticle *particle,       Float_t * const values);
  static void FillVarAODMCParticle(const AliAODMCParticle *particle, Float_t * const values);
  //static void FillVarDielectronPair(const AliDielectronPair *pair,   Float_t * const values);
  static void FillVarKFParticle(const AliKFParticle *pair,           Float_t * const values);

  static void FillVarVEvent(const AliVEvent *event,                  Float_t * const values);
  static void FillVarESDEvent(const AliESDEvent *event,              Float_t * const values);
  static void FillVarAODEvent(const AliAODEvent *event,              Float_t * const values);
  static void FillVarMCEvent(const AliMCEvent *event,                Float_t * const values);
  static void FillVarTPCEventPlane(const AliEventplane *evplane,     Float_t * const values);

private:

  static void InitVZEROCalibrationHistograms(Int_t runNo);
  static void InitVZERORecenteringHistograms(Int_t runNo);
  static void InitZDCRecenteringHistograms(Int_t runNo);

  static AliPIDResponse  *fgPIDResponse;        // PID response object
  static AliVEvent       *fgEvent;              // current event pointer
  static AliEventplane   *fgTPCEventPlane;      // current event tpc plane pointer
  static AliKFVertex     *fgKFVertex;           // kf vertex
  static TProfile        *fgMultEstimatorAvg[6][9];  // multiplicity estimator averages (6 periods x 18 estimators)
  static Double_t         fgTRDpidEffCentRanges[10][4];   // centrality ranges for the TRD pid efficiency histograms
  static TH3D            *fgTRDpidEff[10][4];   // TRD pid efficiencies from conversion electrons
  static TObject         *fgLegEffMap;             // single electron efficiencies
  static TObject         *fgPairEffMap;             // pair efficiencies
  static TBits           *fgFillMap;             // map for requested variable filling
  static TString          fgVZEROCalibrationFile;  // file with VZERO channel-by-channel calibrations
  static TString          fgVZERORecenteringFile;  // file with VZERO Q-vector averages needed for event plane recentering
  static TProfile2D      *fgVZEROCalib[64];           // 1 histogram per VZERO channel
  static TProfile2D      *fgVZERORecentering[2][2];   // 2 VZERO sides x 2 Q-vector components
  static Int_t            fgCurrentRun;               // current run number

  static TString          fgZDCRecenteringFile; // file with ZDC Q-vector averages needed for event plane recentering
  static TProfile3D      *fgZDCRecentering[3][2];   // 2 VZERO sides x 2 Q-vector components

  static Float_t fgData[kNMaxValues];        //! data

  AliCMEVarManager(const AliCMEVarManager &c);
  AliCMEVarManager &operator=(const AliCMEVarManager &c);

  ClassDef(AliCMEVarManager,1);
};


//Inline functions
inline void AliCMEVarManager::Fill(const TObject* object, Float_t * const values)
{
  //
  // Main function to fill all available variables according to the type of particle
  //
  if (!object) return;
  if      (object->IsA() == AliESDtrack::Class())       FillVarESDtrack(static_cast<const AliESDtrack*>(object), values);
  else if (object->IsA() == AliAODTrack::Class())       FillVarAODTrack(static_cast<const AliAODTrack*>(object), values);
  else if (object->IsA() == AliMCParticle::Class())     FillVarMCParticle(static_cast<const AliMCParticle*>(object), values);
  else if (object->IsA() == AliAODMCParticle::Class())  FillVarAODMCParticle(static_cast<const AliAODMCParticle*>(object), values);
  //else if (object->IsA() == AliDielectronPair::Class()) FillVarDielectronPair(static_cast<const AliDielectronPair*>(object), values);
  else if (object->IsA() == AliKFParticle::Class())     FillVarKFParticle(static_cast<const AliKFParticle*>(object),values);
  // Main function to fill all available variables according to the type of event

  else if (object->IsA() == AliVEvent::Class())         FillVarVEvent(static_cast<const AliVEvent*>(object), values);
  else if (object->IsA() == AliESDEvent::Class())       FillVarESDEvent(static_cast<const AliESDEvent*>(object), values);
  else if (object->IsA() == AliAODEvent::Class())       FillVarAODEvent(static_cast<const AliAODEvent*>(object), values);
  else if (object->IsA() == AliMCEvent::Class())        FillVarMCEvent(static_cast<const AliMCEvent*>(object), values);
  else if (object->IsA() == AliEventplane::Class())     FillVarTPCEventPlane(static_cast<const AliEventplane*>(object), values);
//   else printf(Form("AliCMEVarManager::Fill: Type %s is not supported by AliCMEVarManager!", object->ClassName())); //TODO: implement without object needed
}

inline void AliCMEVarManager::FillVarVParticle(const AliVParticle *particle, Float_t * const values)
{
  //
  // Fill track information available in AliVParticle into an array
  //
  values[AliCMEVarManager::kPx]        = particle->Px();
  values[AliCMEVarManager::kPy]        = particle->Py();
  values[AliCMEVarManager::kPz]        = particle->Pz();
  values[AliCMEVarManager::kPt]        = particle->Pt();
  values[AliCMEVarManager::kPtSq]      = particle->Pt()*particle->Pt();
  values[AliCMEVarManager::kP]         = particle->P();

  values[AliCMEVarManager::kXv]        = particle->Xv();
  values[AliCMEVarManager::kYv]        = particle->Yv();
  values[AliCMEVarManager::kZv]        = particle->Zv();

  values[AliCMEVarManager::kOneOverPt] = (particle->Pt()>1.0e-3 ? particle->OneOverPt() : 0.0);
  values[AliCMEVarManager::kPhi]       = TVector2::Phi_0_2pi(particle->Phi());
  values[AliCMEVarManager::kTheta]     = particle->Theta();
  values[AliCMEVarManager::kEta]       = particle->Eta();
  values[AliCMEVarManager::kY]         = particle->Y();

  values[AliCMEVarManager::kE]         = particle->E();
  values[AliCMEVarManager::kM]         = particle->M();
  values[AliCMEVarManager::kCharge]    = particle->Charge();

  values[AliCMEVarManager::kPdgCode]   = particle->PdgCode();

  values[AliCMEVarManager::kRndm]      = gRandom->Rndm();

//   if ( fgEvent ) AliCMEVarManager::Fill(fgEvent, values);
  //for (Int_t i=AliCMEVarManager::kPairMax; i<AliCMEVarManager::kNMaxValues; ++i)
    //values[i]=fgData[i];
}

inline void AliCMEVarManager::FillVarESDtrack(const AliESDtrack *particle, Float_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  AliESDtrack *esdTrack=0x0;
  //Double_t origdEdx=particle->GetTPCsignal();

  // apply ETa correction, remove once this is in the tender
  esdTrack=const_cast<AliESDtrack*>(particle);
  if (!esdTrack) return;
  //esdTrack->SetTPCsignal(origdEdx/AliDielectronPID::GetEtaCorr(esdTrack)/AliDielectronPID::GetCorrValdEdx(),esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());

  Double_t pidProbs[AliPID::kSPECIES];
  // Fill AliESDtrack interface specific information
  Double_t tpcNcls=particle->GetTPCNcls();
  Double_t tpcNclsS = particle->GetTPCnclsS();
  Double_t itsNcls=particle->GetNcls(0);
  Double_t tpcSignalN=particle->GetTPCsignalN();
  Double_t tpcClusFindable=particle->GetTPCNclsF();
  values[AliCMEVarManager::kNclsITS]       = itsNcls; // TODO: get rid of the plain numbers
  values[AliCMEVarManager::kNclsTPC]       = tpcNcls; // TODO: get rid of the plain numbers
  values[AliCMEVarManager::kNclsSTPC]      = tpcNclsS;
  values[AliCMEVarManager::kNclsSFracTPC]  = tpcNcls>0?tpcNclsS/tpcNcls:0;
  values[AliCMEVarManager::kNclsTPCiter1]  = particle->GetTPCNclsIter1(); // TODO: get rid of the plain numbers
  values[AliCMEVarManager::kNFclsTPC]       = tpcClusFindable;
  values[AliCMEVarManager::kNFclsTPCr]      = particle->GetTPCClusterInfo(2,1);
  values[AliCMEVarManager::kNFclsTPCrFrac]  = particle->GetTPCClusterInfo(2);
  values[AliCMEVarManager::kNFclsTPCfCross]= (tpcClusFindable>0)?(particle->GetTPCClusterInfo(2,1)/tpcClusFindable):0;
  values[AliCMEVarManager::kTPCsignalN]    = tpcSignalN;
  values[AliCMEVarManager::kTPCsignalNfrac]= tpcNcls>0?tpcSignalN/tpcNcls:0;
  values[AliCMEVarManager::kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
  values[AliCMEVarManager::kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
  values[AliCMEVarManager::kTRDpidQuality] = particle->GetTRDntrackletsPID();
  values[AliCMEVarManager::kTRDchi2]       = particle->GetTRDchi2();
  values[AliCMEVarManager::kTRDchi2Trklt]  = (particle->GetTRDntrackletsPID() > 0 ? particle->GetTRDchi2() / particle->GetTRDntrackletsPID() : -1.);
  values[AliCMEVarManager::kTRDsignal]     = particle->GetTRDsignal();
  values[AliCMEVarManager::kTPCclsDiff]    = tpcSignalN-tpcNcls;
  values[AliCMEVarManager::kTPCclsSegments] = 0.0;
  UChar_t threshold = 5;
  TBits tpcClusterMap = particle->GetTPCClusterMap();
  UChar_t n=0; UChar_t j=0;
  for(UChar_t i=0; i<8; ++i) {
    n=0;
    for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliCMEVarManager::kTPCclsSegments] += 1.0;
  }

  n=0;
  threshold=0;
  values[AliCMEVarManager::kTPCclsIRO]=0.;
  for(j=0; j<63; ++j) n+=tpcClusterMap.TestBitNumber(j);
  if(n>=threshold) values[AliCMEVarManager::kTPCclsIRO] = n;
  n=0;
  threshold=0;
  values[AliCMEVarManager::kTPCclsORO]=0.;
  for(j=63; j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
  if(n>=threshold) values[AliCMEVarManager::kTPCclsORO] = n;


  values[AliCMEVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  for(Int_t i=0; i<kNTrackingStatus; i++) values[kTrackingStatus+i] = (((particle->GetStatus())&(1<<i)) > 0 ? 1. : 0.);
  values[AliCMEVarManager::kFilterBit]     = 0;

  values[AliCMEVarManager::kTPCchi2Cl] = -1;
  if (tpcNcls>0) values[AliCMEVarManager::kTPCchi2Cl] = particle->GetTPCchi2() / tpcNcls;
  values[AliCMEVarManager::kITSchi2Cl] = -1;
  if (itsNcls>0) values[AliCMEVarManager::kITSchi2Cl] = particle->GetITSchi2() / itsNcls;
  //TRD pidProbs
  particle->GetTRDpid(pidProbs);
  values[AliCMEVarManager::kTRDprobEle]    = pidProbs[AliPID::kElectron];
  values[AliCMEVarManager::kTRDprobPio]    = pidProbs[AliPID::kPion];

  values[AliCMEVarManager::kV0Index0]      = particle->GetV0Index(0);
  values[AliCMEVarManager::kKinkIndex0]    = particle->GetKinkIndex(0);

  Float_t impactParXY, impactParZ;
  particle->GetImpactParameters(impactParXY, impactParZ);
  values[AliCMEVarManager::kImpactParXY]   = impactParXY;
  values[AliCMEVarManager::kImpactParZ]    = impactParZ;


  values[AliCMEVarManager::kPdgCode]=-1;
  values[AliCMEVarManager::kPdgCodeMother]=-1;
  values[AliCMEVarManager::kPdgCodeGrandMother]=-1;
  values[AliCMEVarManager::kHasCocktailMother]=0;
  values[AliCMEVarManager::kHasCocktailGrandMother]=0;

  values[AliCMEVarManager::kNumberOfDaughters]=-999;

  //AliDielectronMC *mc=AliDielectronMC::Instance();
  //if (mc->HasMC()){
  //  if (mc->GetMCTrack(particle)) {
  //    Int_t trkLbl = TMath::Abs(mc->GetMCTrack(particle)->GetLabel());
  //    values[AliCMEVarManager::kPdgCode]           =mc->GetMCTrack(particle)->PdgCode();
  //    values[AliCMEVarManager::kHasCocktailMother] =mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
  //    values[AliCMEVarManager::kPdgCodeMother]     =mc->GetMotherPDG(particle);
  //    AliMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
  //    if(motherMC) values[AliCMEVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);
  //  }
  //  values[AliCMEVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
  //} //if(mc->HasMC())


  values[AliCMEVarManager::kITSsignal]   =   particle->GetITSsignal();

  Double_t itsdEdx[4];
  particle->GetITSdEdxSamples(itsdEdx);

  values[AliCMEVarManager::kITSsignalSSD1]   =   itsdEdx[0];
  values[AliCMEVarManager::kITSsignalSSD2]   =   itsdEdx[1];
  values[AliCMEVarManager::kITSsignalSDD1]   =   itsdEdx[2];
  values[AliCMEVarManager::kITSsignalSDD2]   =   itsdEdx[3];
  values[AliCMEVarManager::kITSclusterMap]   =   particle->GetITSClusterMap();
  values[AliCMEVarManager::kITSLayerFirstCls] = -1.;

  for (Int_t iC=0; iC<6; iC++) {
    if (((particle->GetITSClusterMap()) & (1<<(iC))) > 0) {
      values[AliCMEVarManager::kITSLayerFirstCls] = iC;
      break;
    }
  }


  values[AliCMEVarManager::kTrackLength]   = particle->GetIntegratedLength();
  //dEdx information
  Double_t mom = particle->GetP();
  const AliExternalTrackParam *in=particle->GetInnerParam();
  Double_t ysignedIn=-100;
  if (in) {
    mom = in->GetP();
    ysignedIn=particle->Charge()*in->GetY();
  }
  values[AliCMEVarManager::kPIn]=mom;
  values[AliCMEVarManager::kYsignedIn]=ysignedIn;
  const AliExternalTrackParam *out=particle->GetOuterParam();
  if(out) values[AliCMEVarManager::kPOut] = out->GetP();
  else values[AliCMEVarManager::kPOut] = mom;
  if(out && fgEvent) {
    Double_t localCoord[3]={0.0};
    Bool_t localCoordGood = out->GetXYZAt(298.0, ((AliESDEvent*)fgEvent)->GetMagneticField(), localCoord);
    values[AliCMEVarManager::kTRDphi] = (localCoordGood && TMath::Abs(localCoord[0])>1.0e-6 && TMath::Abs(localCoord[1])>1.0e-6 ? TMath::ATan2(localCoord[1], localCoord[0]) : -999.);
  }
  //if(mc->HasMC() && fgTRDpidEff[0][0]) {
  //  Int_t runNo = (fgEvent ? fgEvent->GetRunNumber() : -1);
  //  Float_t centrality=-1.0;
  //  AliCentrality *esdCentrality = (fgEvent ? fgEvent->GetCentrality() : 0x0);
  //  if(esdCentrality) centrality = esdCentrality->GetCentralityPercentile("V0M");
  //  Double_t effErr=0.0;
  //  values[kTRDpidEffLeg] = GetTRDpidEfficiency(runNo, centrality, values[AliCMEVarManager::kEta],
	//					values[AliCMEVarManager::kTRDphi],
	//					values[AliCMEVarManager::kPOut], effErr);
  //}
  values[AliCMEVarManager::kTPCsignal]=particle->GetTPCsignal();

  values[AliCMEVarManager::kTOFsignal]=particle->GetTOFsignal();

  Double_t l = particle->GetIntegratedLength();  // cm
  Double_t t = particle->GetTOFsignal();
  Double_t t0 = fgPIDResponse->GetTOFResponse().GetTimeZero(); // ps

  if( (l < 360. || l > 800.) || (t <= 0.) || (t0 >999990.0) ) {
	values[AliCMEVarManager::kTOFbeta]=0.0;
  }
  else {
	t -= t0; // subtract the T0
	l *= 0.01;  // cm ->m
	t *= 1e-12; //ps -> s

	Double_t v = l / t;
	Float_t beta = v / TMath::C();
	values[AliCMEVarManager::kTOFbeta]=beta;
  }
  values[AliCMEVarManager::kTOFPIDBit]=(particle->GetStatus()&AliESDtrack::kTOFpid? 1: 0);

  values[AliCMEVarManager::kTOFmismProb] = fgPIDResponse->GetTOFMismatchProbability(particle);

  // nsigma to Electron band
  // TODO: for the moment we set the bethe bloch parameters manually
  //       this should be changed in future!
  values[AliCMEVarManager::kTPCnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron);
  //values[AliCMEVarManager::kTPCnSigmaEle]   =(fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron)
  //                                                   -AliDielectronPID::GetCorrVal()
  //                                                   -AliDielectronPID::GetCntrdCorr(particle)
  //                                                   ) / AliDielectronPID::GetWdthCorr(particle);

  values[AliCMEVarManager::kTPCnSigmaPio]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
  values[AliCMEVarManager::kTPCnSigmaMuo]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
  values[AliCMEVarManager::kTPCnSigmaKao]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
  values[AliCMEVarManager::kTPCnSigmaPro]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);

  values[AliCMEVarManager::kITSnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron);
  //values[AliCMEVarManager::kITSnSigmaEle]   =(fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron)
  //                                                   -AliDielectronPID::GetCntrdCorrITS(particle)
  //                                                   ) / AliDielectronPID::GetWdthCorrITS(particle);

  values[AliCMEVarManager::kITSnSigmaPio]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kPion);
  values[AliCMEVarManager::kITSnSigmaMuo]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kMuon);
  values[AliCMEVarManager::kITSnSigmaKao]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kKaon);
  values[AliCMEVarManager::kITSnSigmaPro]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kProton);

  values[AliCMEVarManager::kTOFnSigmaEle]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kElectron);
  values[AliCMEVarManager::kTOFnSigmaPio]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kPion);
  values[AliCMEVarManager::kTOFnSigmaMuo]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kMuon);
  values[AliCMEVarManager::kTOFnSigmaKao]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kKaon);
  values[AliCMEVarManager::kTOFnSigmaPro]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kProton);

  //EMCAL PID information
  Double_t eop=0;
  Double_t showershape[4]={0.,0.,0.,0.};
//   values[AliCMEVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron);
  //values[AliCMEVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron,eop,showershape);
  values[AliCMEVarManager::kEMCALEoverP]     = eop;
  values[AliCMEVarManager::kEMCALE]          = eop*values[AliCMEVarManager::kP];
  values[AliCMEVarManager::kEMCALNCells]     = showershape[0];
  values[AliCMEVarManager::kEMCALM02]        = showershape[1];
  values[AliCMEVarManager::kEMCALM20]        = showershape[2];
  values[AliCMEVarManager::kEMCALDispersion] = showershape[3];

  values[AliCMEVarManager::kLegEff]        = GetSingleLegEff(values);
  values[AliCMEVarManager::kOneOverLegEff] = (values[AliCMEVarManager::kLegEff]>0.0 ? 1./values[AliCMEVarManager::kLegEff] : 0.0);
  //restore TPC signal if it was changed
  //if (esdTrack) esdTrack->SetTPCsignal(origdEdx,esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());

  //fill info from AliVTrdTrack
  //if(Req(kTRDonlineA)||Req(kTRDonlineLayerMask)||Req(kTRDonlinePID)||Req(kTRDonlinePt)||Req(kTRDonlineStack)||Req(kTRDonlineTrackInTime)||Req(kTRDonlineSector)||Req(kTRDonlineFlagsTiming)||Req(kTRDonlineLabel)||Req(kTRDonlineNTracklets)||Req(kTRDonlineFirstLayer))
    //FillVarVTrdTrack(particle,values);


}

inline void AliCMEVarManager::FillVarAODTrack(const AliAODTrack *particle, Float_t * const values)
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

  // Reset AliESDtrack interface specific information
  if(Req(kNclsITS))      values[AliCMEVarManager::kNclsITS]       = particle->GetITSNcls();
  if(Req(kITSchi2Cl))    values[AliCMEVarManager::kITSchi2Cl]     = -1;
  if(Req(kNclsTPC))      values[AliCMEVarManager::kNclsTPC]       = tpcNcls;
  if(Req(kNclsSTPC))     values[AliCMEVarManager::kNclsSTPC]      = tpcNclsS;
  if(Req(kNclsSFracTPC)) values[AliCMEVarManager::kNclsSFracTPC]  = tpcNcls>0?tpcNclsS/tpcNcls:0;
  if(Req(kNclsTPCiter1)) values[AliCMEVarManager::kNclsTPCiter1]  = tpcNcls; // not really available in AOD
  if(Req(kNFclsTPC)  || Req(kNFclsTPCfCross))  values[AliCMEVarManager::kNFclsTPC]      = particle->GetTPCNclsF();
  if(Req(kNFclsTPCr) || Req(kNFclsTPCfCross))  values[AliCMEVarManager::kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
  if(Req(kNFclsTPCrFrac))  values[AliCMEVarManager::kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
  if(Req(kNFclsTPCfCross)) values[AliCMEVarManager::kNFclsTPCfCross]= (values[kNFclsTPC]>0)?(values[kNFclsTPCr]/values[kNFclsTPC]):0;
  if(Req(kNclsTRD))        values[AliCMEVarManager::kNclsTRD]       = particle->GetNcls(2);
  if(Req(kTRDntracklets))  values[AliCMEVarManager::kTRDntracklets] = 0;
  if(Req(kTRDpidQuality))  values[AliCMEVarManager::kTRDpidQuality] = particle->GetTRDntrackletsPID();
  if(Req(kTRDchi2))        values[AliCMEVarManager::kTRDchi2]       = (particle->GetTRDntrackletsPID()!=0.?particle->GetTRDchi2():-1);
  if(Req(kTRDchi2Trklt))   values[AliCMEVarManager::kTRDchi2Trklt]  = (particle->GetTRDntrackletsPID()>0 ? particle->GetTRDchi2() / particle->GetTRDntrackletsPID() : -1.);
  if(Req(kTRDsignal))      values[AliCMEVarManager::kTRDsignal]     = particle->GetTRDsignal();


  TBits tpcClusterMap = particle->GetTPCClusterMap();
  UChar_t n=0; UChar_t j=0;
  UChar_t threshold = 5;

  values[AliCMEVarManager::kTPCclsSegments] = 0.0;
  if(Req(kTPCclsSegments)) {
    for(UChar_t i=0; i<8; ++i) {
      n=0;
      for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
      if(n>=threshold) values[AliCMEVarManager::kTPCclsSegments] += 1.0;
    }
  }

  values[AliCMEVarManager::kTPCclsIRO]=0.;
  if(Req(kTPCclsIRO)) {
    n=0;
    threshold=0;
    for(j=0; j<63; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliCMEVarManager::kTPCclsIRO] = n;
  }

  values[AliCMEVarManager::kTPCclsORO]=0.;
  if(Req(kTPCclsORO)) {
    n=0;
    threshold=0;
    for(j=63; j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliCMEVarManager::kTPCclsORO] = n;
  }

  // it is stored as normalized to tpcNcls-5 (see AliAnalysisTaskESDfilter)
  if(Req(kTPCchi2Cl))   values[AliCMEVarManager::kTPCchi2Cl]     = (tpcNcls>0)?particle->Chi2perNDF()*(tpcNcls-5)/tpcNcls:-1.;
  //if(Req(kTrackStatus)) values[AliCMEVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  //if(Req(kFilterBit))   values[AliCMEVarManager::kFilterBit]     = (Double_t)particle->GetFilterMap();
  values[AliCMEVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  values[AliCMEVarManager::kFilterBit]     = (Double_t)particle->GetFilterMap();

  //TRD pidProbs
  values[AliCMEVarManager::kTRDprobEle]    = 0;
  values[AliCMEVarManager::kTRDprobPio]    = 0;

  values[AliCMEVarManager::kTPCsignalN]    = 0;
  values[AliCMEVarManager::kTPCsignalNfrac]= 0;


  // Fill AliAODTrack interface information
  //
  Int_t v0Index=-1;
  Int_t kinkIndex=-1;
  if( (Req(kV0Index0) || Req(kKinkIndex0)) && particle->GetProdVertex()) {
    v0Index   = particle->GetProdVertex()->GetType()==AliAODVertex::kV0   ? 1 : 0;
    kinkIndex = particle->GetProdVertex()->GetType()==AliAODVertex::kKink ? 1 : 0;
  }
  values[AliCMEVarManager::kV0Index0]      = v0Index;
  values[AliCMEVarManager::kKinkIndex0]    = kinkIndex;

  Double_t d0z0[2]={-999.0,-999.0};
  if(Req(kImpactParXY) || Req(kImpactParZ)) GetDCA(particle, d0z0);
  values[AliCMEVarManager::kImpactParXY]   = d0z0[0];
  values[AliCMEVarManager::kImpactParZ]    = d0z0[1];

  values[AliCMEVarManager::kPIn]            =  0.;
  values[AliCMEVarManager::kTPCsignal]      =  0.;
  values[AliCMEVarManager::kTPCsignalN]     = -1.;
  values[AliCMEVarManager::kTPCsignalNfrac] = -1.;
  values[AliCMEVarManager::kTPCclsDiff]     = -999.;

  values[AliCMEVarManager::kTOFsignal]=0;
  values[AliCMEVarManager::kTOFbeta]=0;

  values[AliCMEVarManager::kTPCnSigmaEleRaw]=0;
  values[AliCMEVarManager::kTPCnSigmaEle]=0;
  values[AliCMEVarManager::kTPCnSigmaPio]=0;
  values[AliCMEVarManager::kTPCnSigmaMuo]=0;
  values[AliCMEVarManager::kTPCnSigmaKao]=0;
  values[AliCMEVarManager::kTPCnSigmaPro]=0;

  values[AliCMEVarManager::kTOFnSigmaEle]=0;
  values[AliCMEVarManager::kTOFnSigmaPio]=0;
  values[AliCMEVarManager::kTOFnSigmaMuo]=0;
  values[AliCMEVarManager::kTOFnSigmaKao]=0;
  values[AliCMEVarManager::kTOFnSigmaPro]=0;

  if(Req(kITSsignal))        values[AliCMEVarManager::kITSsignal]        =   particle->GetITSsignal();
  if(Req(kITSclusterMap))    values[AliCMEVarManager::kITSclusterMap]    =   particle->GetITSClusterMap();
  if(Req(kITSLayerFirstCls)) values[AliCMEVarManager::kITSLayerFirstCls] = -1.;
  for (Int_t iC=0; iC<6; iC++) {
    if (((particle->GetITSClusterMap()) & (1<<(iC))) > 0) {
      if(Req(kITSLayerFirstCls)) values[AliCMEVarManager::kITSLayerFirstCls] = iC;
      break;
    }
  }

  AliAODPid *pid=const_cast<AliAODPid*>(particle->GetDetPid());
  if (pid) {
    //Double_t origdEdx=pid->GetTPCsignal();
    //overwrite signal
    //pid->SetTPCsignal(origdEdx/AliDielectronPID::GetEtaCorr(particle)/AliDielectronPID::GetCorrValdEdx());

    Double_t tpcSignalN=0.0;
    if(Req(kTPCsignalN) || Req(kTPCsignalNfrac) || Req(kTPCclsDiff)) tpcSignalN = pid->GetTPCsignalN();
    values[AliCMEVarManager::kTPCsignalN]     = tpcSignalN;
    values[AliCMEVarManager::kTPCsignalNfrac] = tpcNcls>0?tpcSignalN/tpcNcls:0;
    values[AliCMEVarManager::kTPCclsDiff]     = tpcSignalN-tpcNcls;

    values[AliCMEVarManager::kPIn]         = pid->GetTPCmomentum();
    if(Req(kTPCsignal))   values[AliCMEVarManager::kTPCsignal]   = pid->GetTPCsignal();
    if(Req(kTOFsignal))   values[AliCMEVarManager::kTOFsignal]   = pid->GetTOFsignal();
    if(Req(kTOFmismProb)) values[AliCMEVarManager::kTOFmismProb] = fgPIDResponse->GetTOFMismatchProbability(particle);

    // TOF beta calculation
    if(Req(kTOFbeta)) {
      Double32_t expt[5];
      particle->GetIntegratedTimes(expt);         // ps
      Double_t l  = TMath::C()* expt[0]*1e-12;    // m
      Double_t t  = pid->GetTOFsignal();          // ps start time subtracted (until v5-02-Rev09)
      AliTOFHeader* tofH=0x0;                     // from v5-02-Rev10 on subtract the start time
      if(fgEvent) tofH = (AliTOFHeader*)fgEvent->GetTOFHeader();
      if(tofH) t -= fgPIDResponse->GetTOFResponse().GetStartTime(particle->P()); // ps

    if( (l < 360.e-2 || l > 800.e-2) || (t <= 0.) ) {
      values[AliCMEVarManager::kTOFbeta]  =0;
    }
    else {
      t *= 1e-12; //ps -> s

	Double_t v = l / t;
	Float_t beta = v / TMath::C();
	values[AliCMEVarManager::kTOFbeta]=beta;
      }
    }

    // nsigma for various detectors
    if(Req(kTPCnSigmaEleRaw)) values[kTPCnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron);
    //if(Req(kTPCnSigmaEle))    values[kTPCnSigmaEle]   =(fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron)
    //                                                    -AliDielectronPID::GetCorrVal()
    //                                                    -AliDielectronPID::GetCntrdCorr(particle)
    //                                                    ) / AliDielectronPID::GetWdthCorr(particle);

    if(Req(kTPCnSigmaPio)) values[kTPCnSigmaPio]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
    if(Req(kTPCnSigmaMuo)) values[kTPCnSigmaMuo]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
    if(Req(kTPCnSigmaKao)) values[kTPCnSigmaKao]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
    if(Req(kTPCnSigmaPro)) values[kTPCnSigmaPro]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);

    if(Req(kITSnSigmaEleRaw)) values[kITSnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron);
    //if(Req(kITSnSigmaEle))    values[kITSnSigmaEle]   =(fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron)
    //                                                    -AliDielectronPID::GetCntrdCorrITS(particle)
    //                                                    ) / AliDielectronPID::GetWdthCorrITS(particle);

    if(Req(kITSnSigmaPio)) values[kITSnSigmaPio]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kPion);
    if(Req(kITSnSigmaMuo)) values[kITSnSigmaMuo]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kMuon);
    if(Req(kITSnSigmaKao)) values[kITSnSigmaKao]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kKaon);
    if(Req(kITSnSigmaPro)) values[kITSnSigmaPro]=fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kProton);

    if(Req(kTOFnSigmaEle)) values[kTOFnSigmaEle]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kElectron);
    if(Req(kTOFnSigmaPio)) values[kTOFnSigmaPio]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kPion);
    if(Req(kTOFnSigmaMuo)) values[kTOFnSigmaMuo]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kMuon);
    if(Req(kTOFnSigmaKao)) values[kTOFnSigmaKao]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kKaon);
    if(Req(kTOFnSigmaPro)) values[kTOFnSigmaPro]=fgPIDResponse->NumberOfSigmasTOF(particle,AliPID::kProton);

    Double_t prob[AliPID::kSPECIES]={0.0};
    // switch computation off since it takes 70% of the CPU time for filling all AODtrack variables
    // TODO: find a solution when this is needed (maybe at fill time in histos, CFcontainer and cut selection)
    if( Req(kTRDprobEle) || Req(kTRDprobPio) ){
      fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob);
      values[AliCMEVarManager::kTRDprobEle]      = prob[AliPID::kElectron];
      values[AliCMEVarManager::kTRDprobPio]      = prob[AliPID::kPion];
    }
    if( Req(kTRDprob2DEle) || Req(kTRDprob2DPio) ){
      fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob, AliTRDPIDResponse::kLQ2D);
      values[AliCMEVarManager::kTRDprob2DEle]    = prob[AliPID::kElectron];
      values[AliCMEVarManager::kTRDprob2DPio]    = prob[AliPID::kPion];
    }
    //restore TPC signal if it was changed
    //pid->SetTPCsignal(origdEdx);
  }

  //EMCAL PID information
  Double_t eop=0;
  Double_t showershape[4]={0.,0.,0.,0.};
//   if(Req()) values[AliCMEVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron);
  if(Req(kEMCALnSigmaEle) || Req(kEMCALE) || Req(kEMCALEoverP) ||
     Req(kEMCALNCells) || Req(kEMCALM02) || Req(kEMCALM20) || Req(kEMCALDispersion))
    values[AliCMEVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron,eop,showershape);
  values[AliCMEVarManager::kEMCALEoverP]     = eop;
  values[AliCMEVarManager::kEMCALE]          = eop*values[AliCMEVarManager::kP];
  values[AliCMEVarManager::kEMCALNCells]     = showershape[0];
  values[AliCMEVarManager::kEMCALM02]        = showershape[1];
  values[AliCMEVarManager::kEMCALM20]        = showershape[2];
  values[AliCMEVarManager::kEMCALDispersion] = showershape[3];

  values[AliCMEVarManager::kPdgCode]=-1;
  values[AliCMEVarManager::kPdgCodeMother]=-1;
  values[AliCMEVarManager::kPdgCodeGrandMother]=-1;
  values[AliCMEVarManager::kHasCocktailMother]=0;
  values[AliCMEVarManager::kHasCocktailGrandMother]=0;

  values[AliCMEVarManager::kNumberOfDaughters]=-1;

  //AliDielectronMC *mc=AliDielectronMC::Instance();
  //if (mc->HasMC()){
  //  if (mc->GetMCTrack(particle)) {
  //    Int_t trkLbl = TMath::Abs(mc->GetMCTrack(particle)->GetLabel());
  //    values[AliCMEVarManager::kPdgCode]           =mc->GetMCTrack(particle)->PdgCode();
  //    values[AliCMEVarManager::kHasCocktailMother] =mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
  //    values[AliCMEVarManager::kPdgCodeMother]     =mc->GetMotherPDG(particle);
  //    AliAODMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
  //    if(motherMC) values[AliCMEVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);
  //  }
  //  values[AliCMEVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
  //} //if(mc->HasMC())

  if(Req(kTOFPIDBit))     values[AliCMEVarManager::kTOFPIDBit]=(particle->GetStatus()&AliESDtrack::kTOFpid? 1: 0);
  values[AliCMEVarManager::kLegEff]=0.0;
  values[AliCMEVarManager::kOneOverLegEff]=0.0;
  if(Req(kLegEff) || Req(kOneOverLegEff)) {
    values[AliCMEVarManager::kLegEff] = GetSingleLegEff(values);
    values[AliCMEVarManager::kOneOverLegEff] = (values[AliCMEVarManager::kLegEff]>0.0 ? 1./values[AliCMEVarManager::kLegEff] : 0.0);
  }

  //fill info from AliVTrdTrack
  if(Req(kTRDonlineA)||Req(kTRDonlineLayerMask)||Req(kTRDonlinePID)||Req(kTRDonlinePt)||Req(kTRDonlineStack)||Req(kTRDonlineSector)||Req(kTRDonlineTrackInTime)||Req(kTRDonlineFlagsTiming)||Req(kTRDonlineLabel)||Req(kTRDonlineNTracklets)||Req(kTRDonlineFirstLayer))
    FillVarVTrdTrack(particle,values);
}

inline void AliCMEVarManager::FillVarVTrdTrack(const AliVParticle *particle, Float_t * const values)
{


  //Initialisation of values
  values[AliCMEVarManager::kTRDonlineLayerMask] = -1.0;
  values[AliCMEVarManager::kTRDonlinePID] = -1.0 ;
  values[AliCMEVarManager::kTRDonlinePt] = 0;
  values[AliCMEVarManager::kTRDonlineStack] = -1.0;
  values[AliCMEVarManager::kTRDonlineSector] = -1.0;
  values[AliCMEVarManager::kTRDonlineTrackInTime] = -1.0;
  values[AliCMEVarManager::kTRDonlineFlagsTiming] = -1.0;
  //	if(Req(kTRDonlineLabel))values[AliCMEVarManager::kTRDonlineLabel] = ; ???
  values[AliCMEVarManager::kTRDonlineNTracklets]= -1.0;
  values[AliCMEVarManager::kTRDonlineFirstLayer] = -1.;


  AliESDtrack *esdtrack = 0x0;
  AliAODTrack *aodtrack = 0x0;
  AliVEvent* ev = 0x0;

  if(particle->IsA() == AliESDtrack::Class()){
    esdtrack=(AliESDtrack*)(particle);
    ev =  (AliVEvent*) (esdtrack->GetESDEvent());
  }
  if(particle->IsA() == AliAODTrack::Class()){
    aodtrack= (AliAODTrack*)(particle);
    ev= (AliVEvent*) (aodtrack->GetAODEvent());
  }

  Int_t ngtutrk=ev->GetNumberOfTrdTracks();
  AliVTrdTrack* gtutrk=0x0;
  //loop over gtu track in order to find right matched track for offline
  for(Int_t i=0;i<ngtutrk;i++){

    gtutrk = ev->GetTrdTrack(i);
    if(gtutrk->GetTrackMatch()==particle){

	  values[AliCMEVarManager::kTRDonlineA] = gtutrk->GetA();

      values[AliCMEVarManager::kTRDonlineLayerMask] = gtutrk->GetLayerMask();
      values[AliCMEVarManager::kTRDonlinePID] = gtutrk->GetPID() ;
      values[AliCMEVarManager::kTRDonlinePt] = gtutrk->Pt();
      values[AliCMEVarManager::kTRDonlineStack] = gtutrk->GetStack();
      values[AliCMEVarManager::kTRDonlineSector] = gtutrk->GetSector();
      values[AliCMEVarManager::kTRDonlineTrackInTime] = gtutrk->GetTrackInTime();
      values[AliCMEVarManager::kTRDonlineFlagsTiming] = gtutrk->GetFlagsTiming();
      values[AliCMEVarManager::kTRDonlineLabel] = gtutrk->GetLabel();
      values[AliCMEVarManager::kTRDonlineNTracklets]= gtutrk->GetNTracklets();

      for (Int_t iC=0; iC<6; iC++) {
	    if (((gtutrk->GetLayerMask()) & (1<<(iC))) > 0) {
	      values[AliCMEVarManager::kTRDonlineFirstLayer] = iC;
	      break;
	    }
	  }
      }//if matching
      //TO DO: what is the initialising value, if no match? -1? is this a problem?, is the PT-signed?
  }//for loop over gtutracks

}

inline void AliCMEVarManager::FillVarMCParticle(const AliMCParticle *particle, Float_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  values[AliCMEVarManager::kNclsITS]       = 0;
  values[AliCMEVarManager::kITSchi2Cl]     = 0;
  values[AliCMEVarManager::kNclsTPC]       = 0;
  values[AliCMEVarManager::kNclsSTPC]      = 0;
  values[AliCMEVarManager::kNclsSFracTPC]  = 0;
  values[AliCMEVarManager::kNclsTPCiter1]  = 0;
  values[AliCMEVarManager::kNFclsTPC]      = 0;
  values[AliCMEVarManager::kNFclsTPCr]     = 0;
  values[AliCMEVarManager::kNFclsTPCrFrac] = 0;
  values[AliCMEVarManager::kNclsTRD]       = 0;
  values[AliCMEVarManager::kTRDntracklets] = 0;
  values[AliCMEVarManager::kTRDpidQuality] = 0;
  values[AliCMEVarManager::kTPCchi2Cl]     = 0;
  values[AliCMEVarManager::kTrackStatus]   = 0;
  values[AliCMEVarManager::kFilterBit]     = 0;
  values[AliCMEVarManager::kTRDprobEle]    = 0;
  values[AliCMEVarManager::kTRDprobPio]    = 0;
  values[AliCMEVarManager::kTPCsignalN]    = 0;
  values[AliCMEVarManager::kTPCclsDiff]    = 0;
  values[AliCMEVarManager::kTPCsignalNfrac]    = 0;
  values[AliCMEVarManager::kImpactParXY]   = 0;
  values[AliCMEVarManager::kImpactParZ]    = 0;
  values[AliCMEVarManager::kImpactPar2D]   = 0;
  values[AliCMEVarManager::kPIn]           = 0;
  values[AliCMEVarManager::kYsignedIn]     = 0;
  values[AliCMEVarManager::kTPCsignal]     = 0;
  values[AliCMEVarManager::kTOFsignal]     = 0;
  values[AliCMEVarManager::kTOFbeta]       = 0;
  values[AliCMEVarManager::kTPCnSigmaEleRaw]  = 0;
  values[AliCMEVarManager::kTPCnSigmaEle]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPio]  = 0;
  values[AliCMEVarManager::kTPCnSigmaMuo]  = 0;
  values[AliCMEVarManager::kTPCnSigmaKao]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPro]  = 0;
  values[AliCMEVarManager::kITSclusterMap] = 0;

  values[AliCMEVarManager::kPdgCode]       = -1;
  values[AliCMEVarManager::kPdgCodeMother] = -1;
  values[AliCMEVarManager::kPdgCodeGrandMother] = -1;
  values[AliCMEVarManager::kHasCocktailMother]=0;
  values[AliCMEVarManager::kHasCocktailGrandMother]=0;

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  // Fill AliMCParticle interface specific information
  //AliDielectronMC *mc=AliDielectronMC::Instance();
  //Int_t trkLbl = TMath::Abs(particle->GetLabel());
  //values[AliCMEVarManager::kPdgCode]           = particle->PdgCode();
  //values[AliCMEVarManager::kHasCocktailMother] = mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
  //values[AliCMEVarManager::kPdgCodeMother]     = mc->GetMotherPDG(particle);
  //AliMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
  //if(motherMC) values[AliCMEVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);


  //values[AliCMEVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);
  //values[AliCMEVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
}


inline void AliCMEVarManager::FillVarMCParticle2(const AliVParticle *p1, const AliVParticle *p2, Float_t * const values) {
  //
  // fill 2 track information starting from MC legs
  //

  values[AliCMEVarManager::kNclsITS]       = 0;
  values[AliCMEVarManager::kITSchi2Cl]     = -1;
  values[AliCMEVarManager::kNclsTPC]       = 0;
  values[AliCMEVarManager::kNclsSTPC]       = 0;
  values[AliCMEVarManager::kNclsSFracTPC]  = 0;
  values[AliCMEVarManager::kNclsTPCiter1]  = 0;
  values[AliCMEVarManager::kNFclsTPC]      = 0;
  values[AliCMEVarManager::kNFclsTPCr]     = 0;
  values[AliCMEVarManager::kNFclsTPCrFrac] = 0;
  values[AliCMEVarManager::kNclsTRD]       = 0;
  values[AliCMEVarManager::kTRDntracklets] = 0;
  values[AliCMEVarManager::kTRDpidQuality] = 0;
  values[AliCMEVarManager::kTPCchi2Cl]     = 0;
  values[AliCMEVarManager::kTrackStatus]   = 0;
  values[AliCMEVarManager::kFilterBit]     = 0;
  values[AliCMEVarManager::kTRDprobEle]    = 0;
  values[AliCMEVarManager::kTRDprobPio]    = 0;
  values[AliCMEVarManager::kTPCsignalN]    = 0;
  values[AliCMEVarManager::kTPCclsDiff]    = 0;
  values[AliCMEVarManager::kTPCsignalNfrac]    = 0;
  values[AliCMEVarManager::kImpactParXY]   = 0;
  values[AliCMEVarManager::kImpactParZ]    = 0;
  values[AliCMEVarManager::kImpactPar2D]   = 0;
  values[AliCMEVarManager::kPIn]           = 0;
  values[AliCMEVarManager::kYsignedIn]     = 0;
  values[AliCMEVarManager::kTPCsignal]     = 0;
  values[AliCMEVarManager::kTPCnSigmaEleRaw]  = 0;
  values[AliCMEVarManager::kTPCnSigmaEle]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPio]  = 0;
  values[AliCMEVarManager::kTPCnSigmaMuo]  = 0;
  values[AliCMEVarManager::kTPCnSigmaKao]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPro]  = 0;
  values[AliCMEVarManager::kITSclusterMap] = 0;

  values[AliCMEVarManager::kPdgCode]       = -1;
  values[AliCMEVarManager::kPdgCodeMother] = -1;
  values[AliCMEVarManager::kHasCocktailMother]=0;

  //AliDielectronMC *mc=AliDielectronMC::Instance();
  //AliVParticle* mother=0x0;
  //Int_t mLabel1 = mc->GetMothersLabel(p1->GetLabel());
  //Int_t mLabel2 = mc->GetMothersLabel(p2->GetLabel());
  //if(mLabel1==mLabel2)
  //  mother = mc->GetMCTrackFromMCEvent(mLabel1);

  //values[AliCMEVarManager::kPseudoProperTime] = -2e10;
  //if(mother) {    // same mother
  //  FillVarVParticle(mother, values);
  //  Double_t vtxX, vtxY, vtxZ;
  //  mc->GetPrimaryVertex(vtxX,vtxY,vtxZ);
  //  Double_t lxy = ((mother->Xv()- vtxX) * mother->Px() +
	//	    (mother->Yv()- vtxY) * mother->Py() )/mother->Pt();
  //  values[AliCMEVarManager::kPseudoProperTime] = lxy*(TDatabasePDG::Instance()->GetParticle(443)->Mass())/mother->Pt();
  //}
  // AliVParticle part
  values[AliCMEVarManager::kPx]        = p1->Px()+p2->Px();
  values[AliCMEVarManager::kPy]        = p1->Py()+p2->Py();
  values[AliCMEVarManager::kPz]        = p1->Pz()+p2->Pz();
  values[AliCMEVarManager::kPt]        = TMath::Sqrt(values[AliCMEVarManager::kPx]*
						      values[AliCMEVarManager::kPx]+
						      values[AliCMEVarManager::kPy]*
						      values[AliCMEVarManager::kPy]);
  values[AliCMEVarManager::kPtSq]      = values[AliCMEVarManager::kPt] * values[AliCMEVarManager::kPt];
  values[AliCMEVarManager::kP]         = TMath::Sqrt(values[AliCMEVarManager::kPt]*
						      values[AliCMEVarManager::kPt]+
						      values[AliCMEVarManager::kPz]*
						      values[AliCMEVarManager::kPz]);

  values[AliCMEVarManager::kXv]        = 0;
  values[AliCMEVarManager::kYv]        = 0;
  values[AliCMEVarManager::kZv]        = 0;

  values[AliCMEVarManager::kOneOverPt] = (values[AliCMEVarManager::kPt]>1.0e-6 ? 1.0/values[AliCMEVarManager::kPt] : 0.0);
  values[AliCMEVarManager::kPhi]       = TVector2::Phi_0_2pi( TMath::ATan2(values[AliCMEVarManager::kPy],values[AliCMEVarManager::kPx]) );
  values[AliCMEVarManager::kTheta]     = TMath::ATan2(values[AliCMEVarManager::kPt],values[AliCMEVarManager::kPz]);
  values[AliCMEVarManager::kEta]       = ((values[AliCMEVarManager::kP]-values[AliCMEVarManager::kPz])>1.0e-6 && (values[AliCMEVarManager::kP]+values[AliCMEVarManager::kPz])>1.0e-6 ? 0.5*TMath::Log((values[AliCMEVarManager::kP]+values[AliCMEVarManager::kPz])/(values[AliCMEVarManager::kP]-values[AliCMEVarManager::kPz])) : -9999.);
  values[AliCMEVarManager::kE]         = p1->E()+p2->E();
  values[AliCMEVarManager::kY]         = ((values[AliCMEVarManager::kE]-values[AliCMEVarManager::kPz])>1.0e-6 && (values[AliCMEVarManager::kE]+values[AliCMEVarManager::kPz])>1.0e-6 ? 0.5*TMath::Log((values[AliCMEVarManager::kE]+values[AliCMEVarManager::kPz])/(values[AliCMEVarManager::kE]-values[AliCMEVarManager::kPz])) : -9999.);
  values[AliCMEVarManager::kCharge]    = p1->Charge()+p2->Charge();

  values[AliCMEVarManager::kM]         = p1->M()*p1->M()+p2->M()*p2->M()+
                       2.0*(p1->E()*p2->E()-p1->Px()*p2->Px()-p1->Py()*p2->Py()-p1->Pz()*p2->Pz());
  values[AliCMEVarManager::kM]         = (values[AliCMEVarManager::kM]>1.0e-8 ? TMath::Sqrt(values[AliCMEVarManager::kM]) : -1.0);

  if ( fgEvent ) AliCMEVarManager::Fill(fgEvent, values);

  //values[AliCMEVarManager::kThetaHE]   = AliDielectronPair::ThetaPhiCM(p1,p2,kTRUE,  kTRUE);
  //values[AliCMEVarManager::kPhiHE]     = AliDielectronPair::ThetaPhiCM(p1,p2,kTRUE,  kFALSE);
  values[AliCMEVarManager::kThetaSqHE]  = values[AliCMEVarManager::kThetaHE] * values[AliCMEVarManager::kThetaHE];
  values[AliCMEVarManager::kCos2PhiHE] = TMath::Cos(2*values[AliCMEVarManager::kPhiHE]);
  //values[AliCMEVarManager::kThetaCS]   = AliDielectronPair::ThetaPhiCM(p1,p2,kFALSE, kTRUE);
  //values[AliCMEVarManager::kPhiCS]     = AliDielectronPair::ThetaPhiCM(p1,p2,kFALSE, kFALSE);
  values[AliCMEVarManager::kThetaSqCS]  = values[AliCMEVarManager::kThetaCS] * values[AliCMEVarManager::kThetaCS];
  values[AliCMEVarManager::kCos2PhiCS] = TMath::Cos(2*values[AliCMEVarManager::kPhiCS]);
  values[AliCMEVarManager::kCosTilPhiHE]  = (values[AliCMEVarManager::kThetaHE]>0)?(TMath::Cos(values[AliCMEVarManager::kPhiHE]-TMath::Pi()/4.)):(TMath::Cos(values[AliCMEVarManager::kPhiHE]-3*TMath::Pi()/4.));
  values[AliCMEVarManager::kCosTilPhiCS]  = (values[AliCMEVarManager::kThetaCS]>0)?(TMath::Cos(values[AliCMEVarManager::kPhiCS]-TMath::Pi()/4.)):(TMath::Cos(values[AliCMEVarManager::kPhiCS]-3*TMath::Pi()/4.));
}


inline void AliCMEVarManager::FillVarAODMCParticle(const AliAODMCParticle *particle, Float_t * const values)
{
  //
  // Fill track information available for histogramming into an array
  //

  values[AliCMEVarManager::kNclsITS]       = 0;
  values[AliCMEVarManager::kITSchi2Cl]     = -1;
  values[AliCMEVarManager::kNclsTPC]       = 0;
  values[AliCMEVarManager::kNclsSTPC]      = 0;
  values[AliCMEVarManager::kNclsSFracTPC]  = 0;
  values[AliCMEVarManager::kNclsTPCiter1]  = 0;
  values[AliCMEVarManager::kNFclsTPC]      = 0;
  values[AliCMEVarManager::kNclsTRD]       = 0;
  values[AliCMEVarManager::kTRDntracklets] = 0;
  values[AliCMEVarManager::kTRDpidQuality] = 0;
  values[AliCMEVarManager::kTPCchi2Cl]     = 0;
  values[AliCMEVarManager::kTrackStatus]   = 0;
  values[AliCMEVarManager::kFilterBit]     = 0;
  values[AliCMEVarManager::kTRDprobEle]    = 0;
  values[AliCMEVarManager::kTRDprobPio]    = 0;
  values[AliCMEVarManager::kTPCsignalN]    = 0;
  values[AliCMEVarManager::kTPCclsDiff]    = 0;
  values[AliCMEVarManager::kTPCsignalNfrac]= 0;
  values[AliCMEVarManager::kImpactParXY]   = 0;
  values[AliCMEVarManager::kImpactParZ]    = 0;
  values[AliCMEVarManager::kImpactPar2D]   = 0;
  values[AliCMEVarManager::kPIn]           = 0;
  values[AliCMEVarManager::kYsignedIn]     = 0;
  values[AliCMEVarManager::kTPCsignal]     = 0;
  values[AliCMEVarManager::kTPCnSigmaEleRaw]  = 0;
  values[AliCMEVarManager::kTPCnSigmaEle]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPio]  = 0;
  values[AliCMEVarManager::kTPCnSigmaMuo]  = 0;
  values[AliCMEVarManager::kTPCnSigmaKao]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPro]  = 0;
  values[AliCMEVarManager::kITSclusterMap] = 0;

  values[AliCMEVarManager::kPdgCode]       = -1;
  values[AliCMEVarManager::kPdgCodeMother] = -1;
  values[AliCMEVarManager::kPdgCodeGrandMother] = -1;
  values[AliCMEVarManager::kHasCocktailMother]=0;
  values[AliCMEVarManager::kHasCocktailGrandMother]=0;

  // Fill common AliVParticle interface information
  FillVarVParticle(particle, values);

  // Fill AliAODMCParticle interface specific information
  //AliDielectronMC *mc=AliDielectronMC::Instance();
  //Int_t trkLbl = TMath::Abs(particle->GetLabel());
  //values[AliCMEVarManager::kPdgCode]           = particle->PdgCode();
  //values[AliCMEVarManager::kHasCocktailMother] = mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
  //values[AliCMEVarManager::kPdgCodeMother]     = mc->GetMotherPDG(particle);
  //AliAODMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
  //if(motherMC) values[AliCMEVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);

  //values[AliCMEVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);
  //values[AliCMEVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);

  //// using AODMCHEader information
  //AliAODMCHeader *mcHeader = (AliAODMCHeader*)fgEvent->FindListObject(AliAODMCHeader::StdBranchName());
  //if(mcHeader) {
  //  values[AliCMEVarManager::kImpactParZ]  = mcHeader->GetVtxZ()-particle->Zv();
  //  values[AliCMEVarManager::kImpactParXY] = TMath::Sqrt(TMath::Power(mcHeader->GetVtxX()-particle->Xv(),2) +
	//							TMath::Power(mcHeader->GetVtxY()-particle->Yv(),2));
  //}

}

//inline void AliCMEVarManager::FillVarDielectronPair(const AliDielectronPair *pair, Float_t * const values)
//{
//  //
//  // Fill pair information available for histogramming into an array
//  //
//
//  values[AliCMEVarManager::kPdgCode]=-1;
//  values[AliCMEVarManager::kPdgCodeMother]=-1;
//  values[AliCMEVarManager::kPdgCodeGrandMother]=-1;
//  values[AliCMEVarManager::kHasCocktailMother]=0;
//  values[AliCMEVarManager::kHasCocktailGrandMother]=0;
//
//  Double_t errPseudoProperTime2 = -1;
//  // Fill common AliVParticle interface information
//  FillVarVParticle(pair, values);
//
//  // Fill AliDielectronPair specific information
//  const AliKFParticle &kfPair = pair->GetKFParticle();
//
//
//  values[AliCMEVarManager::kThetaHE]      = 0.0;
//  values[AliCMEVarManager::kPhiHE]        = 0.0;
//  values[AliCMEVarManager::kThetaSqHE]    = 0.0;
//  values[AliCMEVarManager::kCos2PhiHE]    = 0.0;
//  values[AliCMEVarManager::kCosTilPhiHE]  = 0.0;
//
//  values[AliCMEVarManager::kThetaCS]      = 0.0;
//  values[AliCMEVarManager::kPhiCS]        = 0.0;
//  values[AliCMEVarManager::kThetaSqCS]    = 0.0;
//  values[AliCMEVarManager::kCos2PhiCS]    = 0.0;
//  values[AliCMEVarManager::kCosTilPhiCS]  = 0.0;
//
//  Double_t thetaHE=0;
//  Double_t phiHE=0;
//  Double_t thetaCS=0;
//  Double_t phiCS=0;
//  if(Req(kThetaHE) || Req(kPhiHE) || Req(kThetaCS) || Req(kPhiCS)) {
//    pair->GetThetaPhiCM(thetaHE,phiHE,thetaCS,phiCS);
//
//    values[AliCMEVarManager::kThetaHE]      = thetaHE;
//    values[AliCMEVarManager::kPhiHE]        = phiHE;
//    values[AliCMEVarManager::kThetaSqHE]    = thetaHE * thetaHE;
//    values[AliCMEVarManager::kCos2PhiHE]    = TMath::Cos(2.0*phiHE);
//    values[AliCMEVarManager::kCosTilPhiHE]  = (thetaHE>0)?(TMath::Cos(phiHE-TMath::Pi()/4.)):(TMath::Cos(phiHE-3*TMath::Pi()/4.));
//    values[AliCMEVarManager::kThetaCS]      = thetaCS;
//    values[AliCMEVarManager::kPhiCS]        = phiCS;
//    values[AliCMEVarManager::kThetaSqCS]    = thetaCS * thetaCS;
//    values[AliCMEVarManager::kCos2PhiCS]    = TMath::Cos(2.0*phiCS);
//    values[AliCMEVarManager::kCosTilPhiCS]  = (thetaCS>0)?(TMath::Cos(phiCS-TMath::Pi()/4.)):(TMath::Cos(phiCS-3*TMath::Pi()/4.));
//  }
//
//  if(Req(kChi2NDF))          values[AliCMEVarManager::kChi2NDF]          = kfPair.GetChi2()/kfPair.GetNDF();
//  if(Req(kDecayLength))      values[AliCMEVarManager::kDecayLength]      = kfPair.GetDecayLength();
//  if(Req(kR))                values[AliCMEVarManager::kR]                = kfPair.GetR();
//  if(Req(kOpeningAngle))     values[AliCMEVarManager::kOpeningAngle]     = pair->OpeningAngle();
//  if(Req(kCosPointingAngle)) values[AliCMEVarManager::kCosPointingAngle] = fgEvent ? pair->GetCosPointingAngle(fgEvent->GetPrimaryVertex()) : -1;
//
//  if(Req(kLegDist))   values[AliCMEVarManager::kLegDist]      = pair->DistanceDaughters();
//  if(Req(kLegDistXY)) values[AliCMEVarManager::kLegDistXY]    = pair->DistanceDaughtersXY();
//  if(Req(kDeltaEta))  values[AliCMEVarManager::kDeltaEta]     = pair->DeltaEta();
//  if(Req(kDeltaPhi))  values[AliCMEVarManager::kDeltaPhi]     = pair->DeltaPhi();
//  if(Req(kMerr))      values[AliCMEVarManager::kMerr]         = kfPair.GetErrMass()>1e-30&&kfPair.GetMass()>1e-30?kfPair.GetErrMass()/kfPair.GetMass():1000000;
//
//  values[AliCMEVarManager::kPairType]     = pair->GetType();
//  // Armenteros-Podolanski quantities
//  if(Req(kArmAlpha)) values[AliCMEVarManager::kArmAlpha]     = pair->GetArmAlpha();
//  if(Req(kArmPt))    values[AliCMEVarManager::kArmPt]        = pair->GetArmPt();
//
//  if(Req(kPsiPair))  values[AliCMEVarManager::kPsiPair]      = fgEvent ? pair->PsiPair(fgEvent->GetMagneticField()) : -5;
//  if(Req(kPhivPair)) values[AliCMEVarManager::kPhivPair]      = fgEvent ? pair->PhivPair(fgEvent->GetMagneticField()) : -5;
//  if(Req(kPseudoProperTime) || Req(kPseudoProperTimeErr)) {
//    values[AliCMEVarManager::kPseudoProperTime] =
//      fgEvent ? kfPair.GetPseudoProperDecayTime(*(fgEvent->GetPrimaryVertex()), TDatabasePDG::Instance()->GetParticle(443)->Mass(), &errPseudoProperTime2 ) : -1e10;
//  // values[AliCMEVarManager::kPseudoProperTime] = fgEvent ? pair->GetPseudoProperTime(fgEvent->GetPrimaryVertex()): -1e10;
//    values[AliCMEVarManager::kPseudoProperTimeErr] = (errPseudoProperTime2 > 0) ? TMath::Sqrt(errPseudoProperTime2) : -1e10;
//  }
//
//  // impact parameter
//  Double_t d0z0[2]={-999., -999.};
//  if( (Req(kImpactParXY) || Req(kImpactParZ)) && fgEvent) pair->GetDCA(fgEvent->GetPrimaryVertex(), d0z0);
//  values[AliCMEVarManager::kImpactParXY]   = d0z0[0];
//  values[AliCMEVarManager::kImpactParZ]    = d0z0[1];
//
//
//  //calculate pair dca in sigma and cm
//  values[AliCMEVarManager::kPairDCAsigXY]     = -999.;
//  values[AliCMEVarManager::kPairDCAsigZ]      = -999.;
//  values[AliCMEVarManager::kPairDCAabsXY]     = -999.;
//  values[AliCMEVarManager::kPairDCAabsZ]      = -999.;
//  values[AliCMEVarManager::kPairLinDCAsigXY]  = -999.;
//  values[AliCMEVarManager::kPairLinDCAsigZ]   = -999.;
//  values[AliCMEVarManager::kPairLinDCAabsXY]  = -999.;
//  values[AliCMEVarManager::kPairLinDCAabsZ]   = -999.;
//  values[AliCMEVarManager::kLeg1DCAsigXY]     = -999.;
//  values[AliCMEVarManager::kLeg1DCAabsXY]     = -999.;
//  values[AliCMEVarManager::kLeg1DCAresXY]     = -999.;
//
//  // check if calculation is requested
//  if(Req(kPairDCAsigXY) || Req(kPairDCAsigZ) || Req(kPairDCAabsXY) || Req(kPairDCAabsZ) ||
//     Req(kPairLinDCAsigXY) || Req(kPairLinDCAsigZ) || Req(kPairLinDCAabsXY) || Req(kPairLinDCAabsZ)) {
//
//    // get track references from pair
//    AliVParticle* d1 = pair-> GetFirstDaughterP();
//    AliVParticle* d2 = pair->GetSecondDaughterP();
//
//    if (d1 && d2) {
//      // check for ESD or AOD
//      Bool_t isESD = (d1->IsA() == AliESDtrack::Class());
//
//      if (d1->IsA() == d2->IsA()) { // Don't mix AOD with ESD. Needed because AliAnalysisTaskRandomRejection always creates AliAODTracks (should be fixed).
//
//        ////// first daughter
//        Double_t dca1[2]       = {-999.,-999.};      // xy,z absolute values
//        Double_t dcaSig1[2]    = {-999.,-999.};      // xy,z sigma values
//        Double_t dcaRes1[3]    = {-999.,-999.,-999.};// Covariance matrix
//        //Float_t dcaTPC1[2]    = {-999.,-999.};      // xy,z TPC-only absolute values
//        //Float_t dcaSigTPC1[2] = {-999.,-999.};      // xy,z TPC-only sigma values
//        //Float_t dcaResTPC1[3] = {-999.,-999.,-999.};// Covariance matrix TPC-only
//
//        ////// second daughter
//        Double_t dca2[2]       = {-999.,-999.};      // xy,z absolute values
//        Double_t dcaSig2[2]    = {-999.,-999.};      // xy,z sigma values
//        Double_t dcaRes2[3]    = {-999.,-999.,-999.};// Covariance matrix
//        //Float_t dcaTPC2[2]    = {-999.,-999.};      // xy,z TPC-only absolute values
//        //Float_t dcaSigTPC2[2] = {-999.,-999.};      // xy,z TPC-only sigma values
//        //Float_t dcaResTPC2[3] = {-999.,-999.,-999.};// Covariance matrix TPC-only
//
//        if (isESD) {
//          // 'Float_t' needed for 'virtual void AliESDtrack::GetImpactParameters(Float_t p[2], Float_t cov[3]) const'
//          Float_t dca_tmp[2] = {-999.,-999.};
//          Float_t res_tmp[3] = {-999.,-999.,-999.};
//          static_cast<AliESDtrack*>(d1)->GetImpactParameters(dca_tmp, res_tmp);
//          dca1[0]   =dca_tmp[0]; dca1[1]   =dca_tmp[1];
//          dcaRes1[0]=res_tmp[0]; dcaRes1[1]=res_tmp[1]; dcaRes1[2]=res_tmp[2];
//          //static_cast<AliESDtrack*>(d1)->GetImpactParametersTPC(dcaTPC1, dcaResTPC1);
//
//          dca_tmp[0]=-999.; dca_tmp[1]=-999.;
//          res_tmp[0]=-999.; res_tmp[1]=-999.; res_tmp[2]=-999.;
//          static_cast<AliESDtrack*>(d2)->GetImpactParameters(dca_tmp, res_tmp);
//          dca2[0]   =dca_tmp[0]; dca2[1]   =dca_tmp[1];
//          dcaRes2[0]=res_tmp[0]; dcaRes2[1]=res_tmp[1]; dcaRes2[2]=res_tmp[2];
//          //static_cast<AliESDtrack*>(d2)->GetImpactParametersTPC(dcaTPC2, dcaResTPC2);
//        }
//        else { // AOD
//          GetDCA(static_cast<AliAODTrack*>(d1), dca1, dcaRes1);
//          GetDCA(static_cast<AliAODTrack*>(d2), dca2, dcaRes2);
//        }
//
//        // compute normalized DCAs
//        // neglect the mixed term 'dcaResX[1]'
//        if(dcaRes1[0]>0.) dcaSig1[0] = dca1[0]/TMath::Sqrt(dcaRes1[0]);
//        if(dcaRes1[2]>0.) dcaSig1[1] = dca1[1]/TMath::Sqrt(dcaRes1[2]);
//        //if(dcaResTPC1[0]>0.) dcaSigTPC1[0] = dcaTPC1[0]/TMath::Sqrt(dcaResTPC1[0]);
//        //if(dcaResTPC1[2]>0.) dcaSigTPC1[1] = dcaTPC1[1]/TMath::Sqrt(dcaResTPC1[2]);
//        if(dcaRes2[0]>0.) dcaSig2[0] = dca2[0]/TMath::Sqrt(dcaRes2[0]);
//        if(dcaRes2[2]>0.) dcaSig2[1] = dca2[1]/TMath::Sqrt(dcaRes2[2]);
//        //if(dcaResTPC2[0]>0.) dcaSigTPC2[0] = dcaTPC2[0]/TMath::Sqrt(dcaResTPC2[0]);
//        //if(dcaResTPC2[2]>0.) dcaSigTPC2[1] = dcaTPC2[1]/TMath::Sqrt(dcaResTPC2[2]);
//
//        // set first daughter variables for cross-checks
//        values[AliCMEVarManager::kLeg1DCAabsXY]   = dca1[0];
//        values[AliCMEVarManager::kLeg1DCAsigXY]   = dcaSig1[0];
//        values[AliCMEVarManager::kLeg1DCAresXY]   = dcaRes1[0];
//
//        // set pair dca values
//        // quadratic summation
//        values[AliCMEVarManager::kPairDCAabsXY]       = TMath::Sqrt( (dca1[0]*dca1[0] + dca2[0]*dca2[0]) / 2 );
//        values[AliCMEVarManager::kPairDCAabsZ]        = TMath::Sqrt( (dca1[1]*dca1[1] + dca2[1]*dca2[1]) / 2 );
//        values[AliCMEVarManager::kPairDCAsigXY]       = TMath::Sqrt( (dcaSig1[0]*dcaSig1[0] + dcaSig2[0]*dcaSig2[0]) / 2 );
//        values[AliCMEVarManager::kPairDCAsigZ]        = TMath::Sqrt( (dcaSig1[1]*dcaSig1[1] + dcaSig2[1]*dcaSig2[1]) / 2 );
//        // linear summation
//        values[AliCMEVarManager::kPairLinDCAabsXY]    = (TMath::Abs(dca1[0]) + TMath::Abs(dca2[0])) / 2;
//        values[AliCMEVarManager::kPairLinDCAabsZ]     = (TMath::Abs(dca1[1]) + TMath::Abs(dca2[1])) / 2;
//        values[AliCMEVarManager::kPairLinDCAsigXY]    = (TMath::Abs(dcaSig1[0]) + TMath::Abs(dcaSig2[0])) / 2;
//        values[AliCMEVarManager::kPairLinDCAsigZ]     = (TMath::Abs(dcaSig1[1]) + TMath::Abs(dcaSig2[1])) / 2;
//      }
//    }
//  }
//
//
//  if (!(pair->GetKFUsage())) {
//	//if KF Pairing is not enabled, overwrite values that can be easily derived from legs
//	//use the INDIVIDUAL KF particles as source, which should be a copy of the corresponding properties
//	//the ESDtrack, the reference to the ESDtrack is not (always) accessible in Mixing, while KF
//	//particles are copied in the Pair-Object
//	static const Double_t mElectron = AliPID::ParticleMass(AliPID::kElectron); // MeV
//
//	const AliKFParticle& fD1 = pair->GetKFFirstDaughter();
//	const AliKFParticle& fD2 = pair->GetKFSecondDaughter();
//
//	//Define local buffer variables for leg properties
//	Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
//	Double_t px2=-9999.,py2=-9999.,pz2=-9999.;
//	Double_t e1 =-9999.,e2 =-9999.;
//	Double_t feta1=-9999.;//,fphi1=-9999.;
//	Double_t feta2=-9999.;//,fphi2=-9999.;
//
//	px1 = fD1.GetPx();
//	py1 = fD1.GetPy();
//	pz1 = fD1.GetPz();
//	feta1 = fD1.GetEta();
//	//	fphi1 = fD1.GetPhi();
//
//	px2 = fD2.GetPx();
//	py2 = fD2.GetPy();
//	pz2 = fD2.GetPz();
//	feta2 = fD2.GetEta();
//	//	fphi2 = fD2.GetPhi();
//
//	//Calculate Energy per particle by hand
//	e1 = TMath::Sqrt(mElectron*mElectron+px1*px1+py1*py1+pz1*pz1);
//	e2 = TMath::Sqrt(mElectron*mElectron+px2*px2+py2*py2+pz2*pz2);
//
//	//Now Create TLorentzVector:
//	TLorentzVector lv1,lv2;
//	lv1.SetPxPyPzE(px1,py1,pz1,e1);
//	lv2.SetPxPyPzE(px2,py2,pz2,e2);
//
//	values[AliCMEVarManager::kPx]        = (lv1+lv2).Px();
//	values[AliCMEVarManager::kPy]        = (lv1+lv2).Py();
//	values[AliCMEVarManager::kPz]        = (lv1+lv2).Pz();
//
//	values[AliCMEVarManager::kPt]        =  (lv1+lv2).Pt();
//	values[AliCMEVarManager::kPtSq]      = values[AliCMEVarManager::kPt] * values[AliCMEVarManager::kPt];
//
//	values[AliCMEVarManager::kP]         =  (lv1+lv2).P();
//
//	//Not overwritten, could take event vertex in next iteration
//	values[AliCMEVarManager::kXv]        = (lv1+lv2).X();
//	values[AliCMEVarManager::kYv]        = (lv1+lv2).Y();
//	values[AliCMEVarManager::kZv]        = (lv1+lv2).Z();
//
//	values[AliCMEVarManager::kE]         = (lv1+lv2).E();
//
//
//	values[AliCMEVarManager::kM]         = (lv1+lv2).M();
//
//	values[AliCMEVarManager::kOpeningAngle] =  lv1.Angle(lv2.Vect());
//
//	values[AliCMEVarManager::kOneOverPt] = (values[AliCMEVarManager::kPt]>0. ? 1./values[AliCMEVarManager::kPt] : -9999.);
//	values[AliCMEVarManager::kPhi]       = TVector2::Phi_0_2pi( (lv1+lv2).Phi() );
//	values[AliCMEVarManager::kEta]       = (lv1+lv2).Eta();
//
//	values[AliCMEVarManager::kY]       = (lv1+lv2).Rapidity();
//
//	for (Int_t i=AliCMEVarManager::kPairMax; i<AliCMEVarManager::kNMaxValues; ++i)
//	  values[i]=fgData[i];
//
//	// Fill AliDielectronPair specific information
//	values[AliCMEVarManager::kDeltaEta]     = TMath::Abs(feta1 -feta2 );
//	values[AliCMEVarManager::kDeltaPhi]     = lv1.DeltaPhi(lv2);
//	values[AliCMEVarManager::kPairType]     = pair->GetType();
//
//	/*
//	//Also not overwritten, still coming from KF particle
//	//where needed to be replaced by independent determination
//	values[AliCMEVarManager::kCharge]    = 0.;
//	values[AliCMEVarManager::kPdgCode]   = 0.;
//	values[AliCMEVarManager::kChi2NDF]      = 0.;
//	values[AliCMEVarManager::kDecayLength]  = 0.;
//	values[AliCMEVarManager::kR]            = 0.;
//	values[AliCMEVarManager::kCosPointingAngle] = 0.;
//	values[AliCMEVarManager::kThetaHE]      = 0.;
//	values[AliCMEVarManager::kPhiHE]        = 0.;
//	values[AliCMEVarManager::kThetaSqHE]    = 0.;
//	values[AliCMEVarManager::kCos2PhiHE]    = 0.;
//	values[AliCMEVarManager::kCosTilPhiHE]  = 0.;
//	values[AliCMEVarManager::kThetaCS]      = 0.;
//	values[AliCMEVarManager::kPhiCS]        = 0.;
//	values[AliCMEVarManager::kThetaSqCS]    = 0.;
//	values[AliCMEVarManager::kCos2PhiCS]    = 0.;
//	values[AliCMEVarManager::kCosTilPhiCS]  = 0.;
//	values[AliCMEVarManager::kLegDist]      = 0.;
//	values[AliCMEVarManager::kLegDistXY]    = 0.;
//	values[AliCMEVarManager::kMerr]         = 0.;
//	values[AliCMEVarManager::kPseudoProperTime] = 0.;
//	values[AliCMEVarManager::kPseudoProperTimeErr] = 0.;
//	//Fill in Taku's PhiV?
//	values[AliCMEVarManager::kPsiPair]      = 0.;
//
//	 */
//  }
//  //common, regardless of calculation method
//
//  // Flow quantities
//  Double_t phi=values[AliCMEVarManager::kPhi];
//  if(Req(kCosPhiH2)) values[AliCMEVarManager::kCosPhiH2] = TMath::Cos(2*phi);
//  if(Req(kSinPhiH2)) values[AliCMEVarManager::kSinPhiH2] = TMath::Sin(2*phi);
//  Double_t delta=0.0;
//  // v2 with respect to VZERO-A event plane
//  delta = TVector2::Phi_mpi_pi(phi - fgData[AliCMEVarManager::kV0ArpH2]);
//  if(Req(kV0ArpH2FlowV2))   values[AliCMEVarManager::kV0ArpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
//  if(Req(kDeltaPhiV0ArpH2)) values[AliCMEVarManager::kDeltaPhiV0ArpH2] = delta;
//  // v2 with respect to VZERO-C event plane
//  delta = TVector2::Phi_mpi_pi(phi - fgData[AliCMEVarManager::kV0CrpH2]);
//  if(Req(kV0CrpH2FlowV2))   values[AliCMEVarManager::kV0CrpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
//  if(Req(kDeltaPhiV0CrpH2)) values[AliCMEVarManager::kDeltaPhiV0CrpH2] = delta;
//  // v2 with respect to the combined VZERO-A and VZERO-C event plane
//  delta = TVector2::Phi_mpi_pi(phi - fgData[AliCMEVarManager::kV0ACrpH2]);
//  if(Req(kV0ACrpH2FlowV2))   values[AliCMEVarManager::kV0ACrpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
//  if(Req(kDeltaPhiV0ACrpH2)) values[AliCMEVarManager::kDeltaPhiV0ACrpH2] = delta;
//
//
//  // quantities using the values of  AliEPSelectionTask , interval [-pi,+pi]
//  values[AliCMEVarManager::kDeltaPhiv0ArpH2]  = TVector2::Phi_mpi_pi(phi - values[AliCMEVarManager::kv0ArpH2]);
//  values[AliCMEVarManager::kDeltaPhiv0CrpH2]  = TVector2::Phi_mpi_pi(phi - values[AliCMEVarManager::kv0CrpH2]);
//  values[AliCMEVarManager::kDeltaPhiv0ACrpH2] = TVector2::Phi_mpi_pi(phi - values[AliCMEVarManager::kv0ACrpH2]);
//  values[AliCMEVarManager::kDeltaPhiTPCrpH2]  = TVector2::Phi_mpi_pi(phi - values[AliCMEVarManager::kTPCrpH2]);
//  values[AliCMEVarManager::kv0ACrpH2FlowV2]   = TMath::Cos( 2.*values[AliCMEVarManager::kDeltaPhiv0ACrpH2] );
//  values[AliCMEVarManager::kv0ArpH2FlowV2]    = TMath::Cos( 2.*values[AliCMEVarManager::kDeltaPhiv0ArpH2] );
//  values[AliCMEVarManager::kv0CrpH2FlowV2]    = TMath::Cos( 2.*values[AliCMEVarManager::kDeltaPhiv0CrpH2] );
//  values[AliCMEVarManager::kTPCrpH2FlowV2]    = TMath::Cos( 2.*values[AliCMEVarManager::kDeltaPhiTPCrpH2] );
//  values[AliCMEVarManager::kTPCrpH2FlowV2Sin] = TMath::Sin( 2.*values[AliCMEVarManager::kDeltaPhiTPCrpH2] );
//
//  //calculate inner product of strong Mag and ee plane
//  if(Req(kPairPlaneMagInPro)) values[AliCMEVarManager::kPairPlaneMagInPro] = pair->PairPlaneMagInnerProduct(values[AliCMEVarManager::kZDCACrpH1]);
//
//  //Calculate the angle between electrons decay plane and variables 1-4
//  if(Req(kPairPlaneAngle1A)) values[AliCMEVarManager::kPairPlaneAngle1A] = pair->GetPairPlaneAngle(values[kv0ArpH2],1);
//  if(Req(kPairPlaneAngle2A)) values[AliCMEVarManager::kPairPlaneAngle2A] = pair->GetPairPlaneAngle(values[kv0ArpH2],2);
//  if(Req(kPairPlaneAngle3A)) values[AliCMEVarManager::kPairPlaneAngle3A] = pair->GetPairPlaneAngle(values[kv0ArpH2],3);
//  if(Req(kPairPlaneAngle4A)) values[AliCMEVarManager::kPairPlaneAngle4A] = pair->GetPairPlaneAngle(values[kv0ArpH2],4);
//
//  if(Req(kPairPlaneAngle1C)) values[AliCMEVarManager::kPairPlaneAngle1C] = pair->GetPairPlaneAngle(values[kv0CrpH2],1);
//  if(Req(kPairPlaneAngle2C)) values[AliCMEVarManager::kPairPlaneAngle2C] = pair->GetPairPlaneAngle(values[kv0CrpH2],2);
//  if(Req(kPairPlaneAngle3C)) values[AliCMEVarManager::kPairPlaneAngle3C] = pair->GetPairPlaneAngle(values[kv0CrpH2],3);
//  if(Req(kPairPlaneAngle4C)) values[AliCMEVarManager::kPairPlaneAngle4C] = pair->GetPairPlaneAngle(values[kv0CrpH2],4);
//
//  if(Req(kPairPlaneAngle1AC)) values[AliCMEVarManager::kPairPlaneAngle1AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],1);
//  if(Req(kPairPlaneAngle2AC)) values[AliCMEVarManager::kPairPlaneAngle2AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],2);
//  if(Req(kPairPlaneAngle3AC)) values[AliCMEVarManager::kPairPlaneAngle3AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],3);
//  if(Req(kPairPlaneAngle4AC)) values[AliCMEVarManager::kPairPlaneAngle4AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],4);
//
//  //Random reaction plane
//  values[AliCMEVarManager::kRandomRP] = gRandom->Uniform(-TMath::Pi()/2.0,TMath::Pi()/2.0);
//  //delta phi of pair fron random reaction plane
//  values[AliCMEVarManager::kDeltaPhiRandomRP] = phi - values[kRandomRP];
//  // keep the interval [-pi,+pi]
//  if ( values[AliCMEVarManager::kDeltaPhiRandomRP] > TMath::Pi() )
//    values[AliCMEVarManager::kDeltaPhiRandomRP] -= TMath::TwoPi();
//
//  if(Req(kPairPlaneAngle1Ran)) values[AliCMEVarManager::kPairPlaneAngle1Ran]= pair->GetPairPlaneAngle(values[kRandomRP],1);
//  if(Req(kPairPlaneAngle2Ran)) values[AliCMEVarManager::kPairPlaneAngle2Ran]= pair->GetPairPlaneAngle(values[kRandomRP],2);
//  if(Req(kPairPlaneAngle3Ran)) values[AliCMEVarManager::kPairPlaneAngle3Ran]= pair->GetPairPlaneAngle(values[kRandomRP],3);
//  if(Req(kPairPlaneAngle4Ran)) values[AliCMEVarManager::kPairPlaneAngle4Ran]= pair->GetPairPlaneAngle(values[kRandomRP],4);
//
//
//
//  //AliDielectronMC *mc=AliDielectronMC::Instance();
//
//  //if (mc->HasMC()){
//  //  values[AliCMEVarManager::kPseudoProperTimeResolution] = -10.0e+10;
//  //  Bool_t samemother =  mc->HaveSameMother(pair);
//  //  values[AliCMEVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(pair);
//  //  values[AliCMEVarManager::kHaveSameMother] = samemother ;
//
//  //  // fill kPseudoProperTimeResolution
//  //  values[AliCMEVarManager::kPseudoProperTimeResolution] = -1e10;
//  //  // values[AliCMEVarManager::kPseudoProperTimePull] = -1e10;
//  //  if(samemother && fgEvent) {
//  //    if(pair->GetFirstDaughterP()->GetLabel() > 0) {
//  //      const AliVParticle *motherMC = 0x0;
//  //      if(fgEvent->IsA() == AliESDEvent::Class())  motherMC = (AliMCParticle*)mc->GetMCTrackMother((AliESDtrack*)pair->GetFirstDaughterP());
//  //      else if(fgEvent->IsA() == AliAODEvent::Class())  motherMC = (AliAODMCParticle*)mc->GetMCTrackMother((AliAODTrack*)pair->GetFirstDaughterP());
//  //      Double_t vtxX, vtxY, vtxZ;
//	//if(motherMC && mc->GetPrimaryVertex(vtxX,vtxY,vtxZ)) {
//	//  Int_t motherLbl = motherMC->GetLabel();
//	//  values[AliCMEVarManager::kHasCocktailMother]=mc->CheckParticleSource(motherLbl, AliDielectronSignalMC::kDirect);
//  //    	  const Double_t lxyMC = ( (motherMC->Xv() - vtxX) * motherMC->Px() +
//  //                                 (motherMC->Yv() - vtxY) * motherMC->Py()   ) / motherMC->Pt();
//	//  const Double_t pseudoMC = lxyMC * (TDatabasePDG::Instance()->GetParticle(443)->Mass())/motherMC->Pt();
//	//  values[AliCMEVarManager::kPseudoProperTimeResolution] = values[AliCMEVarManager::kPseudoProperTime] - pseudoMC;
//  //        if (errPseudoProperTime2 > 0)
//  //          values[AliCMEVarManager::kPseudoProperTimePull] = values[AliCMEVarManager::kPseudoProperTimeResolution]/sqrt(errPseudoProperTime2);
//  //    }
//  //    }
//  //  }
//
//	values[AliCMEVarManager::kTRDpidEffPair] = 0.;
//	if (fgTRDpidEff[0][0]){
//	  Double_t valuesLeg1[AliCMEVarManager::kNMaxValues];
//	  Double_t valuesLeg2[AliCMEVarManager::kNMaxValues];
//	  AliVParticle* leg1 = pair->GetFirstDaughterP();
//	  AliVParticle* leg2 = pair->GetSecondDaughterP();
//	  if (leg1 && leg2){
//		Fill(leg1, valuesLeg1);
//		Fill(leg2, valuesLeg2);
//		values[AliCMEVarManager::kTRDpidEffPair] = valuesLeg1[AliCMEVarManager::kTRDpidEffLeg]*valuesLeg2[AliCMEVarManager::kTRDpidEffLeg];
//	  }
//	}
//
//
//  }//if (mc->HasMC())
//
//  AliVParticle* leg1 = pair->GetFirstDaughterP();
//  AliVParticle* leg2 = pair->GetSecondDaughterP();
//  if (leg1)
//    values[AliCMEVarManager::kMomAsymDau1] = (values[AliCMEVarManager::kP] != 0)? leg1->P()  / values[AliCMEVarManager::kP]: 0;
//  else
//    values[AliCMEVarManager::kMomAsymDau1] = -9999.;
//  if (leg2)
//    values[AliCMEVarManager::kMomAsymDau2] = (values[AliCMEVarManager::kP] != 0)? leg2->P()  / values[AliCMEVarManager::kP]: 0;
//  else
//    values[AliCMEVarManager::kMomAsymDau2] = -9999.;
//
//  Double_t valuesLeg1[AliCMEVarManager::kNMaxValues];
//  Double_t valuesLeg2[AliCMEVarManager::kNMaxValues];
//  values[AliCMEVarManager::kPairEff]=0.0;
//  values[AliCMEVarManager::kOneOverPairEff]=0.0;
//  values[AliCMEVarManager::kOneOverPairEffSq]=0.0;
//  if (leg1 && leg2 && fgLegEffMap) {
//    Fill(leg1, valuesLeg1);
//    Fill(leg2, valuesLeg2);
//    values[AliCMEVarManager::kPairEff] = valuesLeg1[AliCMEVarManager::kLegEff] *valuesLeg2[AliCMEVarManager::kLegEff];
//  }
//  else if(fgPairEffMap) {
//    values[AliCMEVarManager::kPairEff] = GetPairEff(values);
//  }
//  if(fgLegEffMap || fgPairEffMap) {
//    values[AliCMEVarManager::kOneOverPairEff] = (values[AliCMEVarManager::kPairEff]>0.0 ? 1./values[AliCMEVarManager::kPairEff] : 1.0);
//    values[AliCMEVarManager::kOneOverPairEffSq] = (values[AliCMEVarManager::kPairEff]>0.0 ? 1./values[AliCMEVarManager::kPairEff]/values[AliCMEVarManager::kPairEff] : 1.0);
//  }
//
//  if(kRndmPair) values[AliCMEVarManager::kRndmPair] = gRandom->Rndm();
//}

inline void AliCMEVarManager::FillVarKFParticle(const AliKFParticle *particle, Float_t * const values)
{
  //
  // Fill track information available in AliVParticle into an array
  //
  values[AliCMEVarManager::kPx]        = particle->GetPx();
  values[AliCMEVarManager::kPy]        = particle->GetPy();
  values[AliCMEVarManager::kPz]        = particle->GetPz();
  values[AliCMEVarManager::kPt]        = particle->GetPt();
  values[AliCMEVarManager::kPtSq]      = particle->GetPt() * particle->GetPt();
  values[AliCMEVarManager::kP]         = particle->GetP();

  values[AliCMEVarManager::kXv]        = particle->GetX();
  values[AliCMEVarManager::kYv]        = particle->GetY();
  values[AliCMEVarManager::kZv]        = particle->GetZ();

  values[AliCMEVarManager::kOneOverPt] = 0;
  values[AliCMEVarManager::kPhi]       = TVector2::Phi_0_2pi(particle->GetPhi()); // interval [0,+2pi]
  values[AliCMEVarManager::kTheta]     = 0.;
  values[AliCMEVarManager::kEta]       = particle->GetEta();
  values[AliCMEVarManager::kY]         = ((particle->GetE()*particle->GetE()-particle->GetPx()*particle->GetPx()-particle->GetPy()*particle->GetPy()-particle->GetPz()*particle->GetPz())>0.) ? TLorentzVector(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE()).Rapidity() : -1111.;

  values[AliCMEVarManager::kE]         = particle->GetE();
  values[AliCMEVarManager::kM]         = particle->GetMass();
  values[AliCMEVarManager::kCharge]    = particle->GetQ();

  values[AliCMEVarManager::kNclsITS]       = 0;
  values[AliCMEVarManager::kITSchi2Cl]     = -1;
  values[AliCMEVarManager::kNclsTPC]       = 0;
  values[AliCMEVarManager::kNclsSTPC]      = 0;
  values[AliCMEVarManager::kNclsSFracTPC]  = 0;
  values[AliCMEVarManager::kNclsTPCiter1]  = 0;
  values[AliCMEVarManager::kNFclsTPC]      = 0;
  values[AliCMEVarManager::kNclsTRD]       = 0;
  values[AliCMEVarManager::kTRDntracklets] = 0;
  values[AliCMEVarManager::kTRDpidQuality] = 0;
  values[AliCMEVarManager::kTPCchi2Cl]     = 0;
  values[AliCMEVarManager::kTrackStatus]   = 0;
  values[AliCMEVarManager::kFilterBit]     = 0;
  values[AliCMEVarManager::kTRDprobEle]    = 0;
  values[AliCMEVarManager::kTRDprobPio]    = 0;
  values[AliCMEVarManager::kTPCsignalN]    = 0;
  values[AliCMEVarManager::kTPCclsDiff]    = 0;
  values[AliCMEVarManager::kTPCsignalNfrac]= 0;
  values[AliCMEVarManager::kImpactParXY]   = 0;
  values[AliCMEVarManager::kImpactParZ]    = 0;
  values[AliCMEVarManager::kImpactPar2D]   = 0;
  values[AliCMEVarManager::kPIn]           = 0;
  values[AliCMEVarManager::kYsignedIn]     = 0;
  values[AliCMEVarManager::kTPCsignal]     = 0;
  values[AliCMEVarManager::kTOFsignal]     = 0;
  values[AliCMEVarManager::kTOFbeta]       = 0;
  values[AliCMEVarManager::kTPCnSigmaEleRaw]  = 0;
  values[AliCMEVarManager::kTPCnSigmaEle]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPio]  = 0;
  values[AliCMEVarManager::kTPCnSigmaMuo]  = 0;
  values[AliCMEVarManager::kTPCnSigmaKao]  = 0;
  values[AliCMEVarManager::kTPCnSigmaPro]  = 0;
  values[AliCMEVarManager::kITSclusterMap] = 0;

  values[AliCMEVarManager::kPdgCode]       = -1;
  values[AliCMEVarManager::kPdgCodeMother] = -1;
  values[AliCMEVarManager::kPdgCodeGrandMother] = -1;
  values[AliCMEVarManager::kHasCocktailMother]=0;
  values[AliCMEVarManager::kHasCocktailGrandMother]=0;

//   if ( fgEvent ) AliCMEVarManager::Fill(fgEvent, values);
  for (Int_t i=AliCMEVarManager::kPairMax; i<AliCMEVarManager::kNMaxValues; ++i)
    values[i]=fgData[i];

}

inline void AliCMEVarManager::FillVarVEvent(const AliVEvent *event, Float_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //
  values[AliCMEVarManager::kRunNumber]    = event->GetRunNumber();
  if(fgCurrentRun!=event->GetRunNumber()) {
    if(fgVZEROCalibrationFile.Contains(".root")) InitVZEROCalibrationHistograms(event->GetRunNumber());
    if(fgVZERORecenteringFile.Contains(".root")) InitVZERORecenteringHistograms(event->GetRunNumber());
    if(fgZDCRecenteringFile.Contains(".root")) InitZDCRecenteringHistograms(event->GetRunNumber());
    fgCurrentRun=event->GetRunNumber();
  }
  values[AliCMEVarManager::kMixingBin]=0;

  const AliVVertex *primVtx = event->GetPrimaryVertex();

  values[AliCMEVarManager::kXvPrim]       = 0;
  values[AliCMEVarManager::kYvPrim]       = 0;
  values[AliCMEVarManager::kZvPrim]       = 0;
  values[AliCMEVarManager::kNVtxContrib]  = 0;
//   values[AliCMEVarManager::kChi2NDF]      = 0; //This is the pair value!!!

  values[AliCMEVarManager::kNTrk]            = 0;
  values[AliCMEVarManager::kNVtxContrib]     = 0;
  values[AliCMEVarManager::kNacc]            = 0;
  values[AliCMEVarManager::kMatchEffITSTPC]  = 0;
  values[AliCMEVarManager::kNaccTrcklts]     = 0;
  values[AliCMEVarManager::kNaccTrcklts10]   = 0;
  values[AliCMEVarManager::kNaccTrcklts0916] = 0;
  values[AliCMEVarManager::kNevents]         = 0; //always fill bin 0;
  values[AliCMEVarManager::kRefMult]         = 0;
  values[AliCMEVarManager::kRefMultTPConly]  = 0;

  if (primVtx){
    //    printf("prim vertex reco: %f \n",primVtx->GetX());
    values[AliCMEVarManager::kXvPrim]       = primVtx->GetX();
    values[AliCMEVarManager::kYvPrim]       = primVtx->GetY();
    values[AliCMEVarManager::kZvPrim]       = primVtx->GetZ();
    values[AliCMEVarManager::kNVtxContrib]  = primVtx->GetNContributors();
  }
  //   values[AliCMEVarManager::kChi2NDF]      = primVtx->GetChi2perNDF(); //this is the pair value

  // online and offline trigger maps
  values[AliCMEVarManager::kTriggerInclONL]     = event->GetTriggerMask();
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  UInt_t maskOff = ((AliInputEventHandler*)man->GetInputEventHandler())->IsEventSelected();
  values[AliCMEVarManager::kTriggerInclOFF]     = maskOff;
  values[AliCMEVarManager::kTriggerExclOFF]        = -1;
  for(Int_t i=0; i<30; i++) { if(maskOff==BIT(i)) values[AliCMEVarManager::kTriggerExclOFF]=i; }

  values[AliCMEVarManager::kNTrk]            = event->GetNumberOfTracks();


  // event plane quantities from the AliEPSelectionTask
  for(Int_t ivar=AliCMEVarManager::kv0ArpH2; ivar<=kv0C0v0C3DiffH2;   ivar++) values[ivar] = 0.0; // v0  variables
  for(Int_t ivar=AliCMEVarManager::kTPCxH2;  ivar<=kTPCsub12DiffH2uc; ivar++) values[ivar] = 0.0; // tpc variables

  // ep angle interval [todo, fill]
  AliEventplane *ep = const_cast<AliVEvent*>(event)->GetEventplane();
  if(ep) {

    // TPC event plane quantities (uncorrected)
    TVector2 *qstd  = ep->GetQVector();  // This is the "standard" Q-Vector for TPC
    TVector2 *qsub1 = ep->GetQsub1();    // random subevent plane
    TVector2 *qsub2 = ep->GetQsub2();
    if(qstd && qsub1 && qsub2) {
      values[AliCMEVarManager::kTPCxH2uc]       = qstd->X();
      values[AliCMEVarManager::kTPCyH2uc]       = qstd->Y();
      values[AliCMEVarManager::kTPCmagH2uc]     = qstd->Mod();
      values[AliCMEVarManager::kTPCrpH2uc]      = TVector2::Phi_mpi_pi(qstd->Phi())/2;
      values[AliCMEVarManager::kTPCsub1xH2uc]   = qsub1->X();
      values[AliCMEVarManager::kTPCsub1yH2uc]   = qsub1->Y();
      values[AliCMEVarManager::kTPCsub1rpH2uc]  = TVector2::Phi_mpi_pi(qsub1->Phi())/2;
      values[AliCMEVarManager::kTPCsub2xH2uc]   = qsub2->X();
      values[AliCMEVarManager::kTPCsub2yH2uc]   = qsub2->Y();
      values[AliCMEVarManager::kTPCsub2rpH2uc]  = TVector2::Phi_mpi_pi(qsub2->Phi())/2;

      values[AliCMEVarManager::kTPCsub12DiffH2uc] = TMath::Cos( 2.*(values[AliCMEVarManager::kTPCsub1rpH2uc] -
									   values[AliCMEVarManager::kTPCsub2rpH2uc]) );
    }

    // VZERO event plane
    TVector2 qvec;
    Double_t qx = 0, qy = 0;
    ep->CalculateVZEROEventPlane(event,10, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0ACrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliCMEVarManager::kv0ACxH2]   = qvec.X();
    values[AliCMEVarManager::kv0ACyH2]   = qvec.Y();
    values[AliCMEVarManager::kv0ACmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 8, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0ArpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliCMEVarManager::kv0AxH2]   = qvec.X();
    values[AliCMEVarManager::kv0AyH2]   = qvec.Y();
    values[AliCMEVarManager::kv0AmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 9, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0CrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliCMEVarManager::kv0CxH2]   = qvec.X();
    values[AliCMEVarManager::kv0CyH2]   = qvec.Y();
    values[AliCMEVarManager::kv0CmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 0, 0, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0C0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep->CalculateVZEROEventPlane(event, 3, 3, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0C3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep->CalculateVZEROEventPlane(event, 4, 4, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0A0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep->CalculateVZEROEventPlane(event, 7, 7, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0A3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
  } //if: eventplane

  // ESD VZERO information
  AliVVZERO* vzeroData = event->GetVZEROData();
  values[AliCMEVarManager::kMultV0A] = 0.0;
  values[AliCMEVarManager::kMultV0C] = 0.0;
  values[AliCMEVarManager::kEqMultV0A] = 0.0;
  values[AliCMEVarManager::kEqMultV0C] = 0.0;
  values[AliCMEVarManager::kAdcV0A]  = 0.0;
  values[AliCMEVarManager::kAdcV0C]  = 0.0;
  for(Int_t i=0; i<32; ++i) {
    values[AliCMEVarManager::kVZEROchMult+i] = vzeroData->GetMultiplicity(i);
    values[AliCMEVarManager::kVZEROchMult+32+i] = vzeroData->GetMultiplicity(i+32);
    //values[AliCMEVarManager::kVZEROchMult+i] = event->GetVZEROEqMultiplicity(i);
    //values[AliCMEVarManager::kVZEROchMult+32+i] = event->GetVZEROEqMultiplicity(i+32);
    values[AliCMEVarManager::kMultV0A] += vzeroData->GetMultiplicityV0A(i);
    values[AliCMEVarManager::kMultV0C] += vzeroData->GetMultiplicityV0C(i);
    values[AliCMEVarManager::kEqMultV0A] += event->GetVZEROEqMultiplicity(i);
    values[AliCMEVarManager::kEqMultV0C] += event->GetVZEROEqMultiplicity(i+32);
    //values[AliCMEVarManager::kAdcV0A] += vzeroData->GetAdcV0A(i);
    //values[AliCMEVarManager::kAdcV0C] += vzeroData->GetAdcV0C(i);
  }
  values[AliCMEVarManager::kMultV0] = values[AliCMEVarManager::kMultV0A] + values[AliCMEVarManager::kMultV0C];
  values[AliCMEVarManager::kEqMultV0] = values[AliCMEVarManager::kEqMultV0A] + values[AliCMEVarManager::kEqMultV0C];
  values[AliCMEVarManager::kAdcV0] = values[AliCMEVarManager::kAdcV0A] + values[AliCMEVarManager::kAdcV0C];
  // VZERO event plane quantities
  Double_t qvec[3]={0.0};
  GetVzeroRP(event, qvec,0);      // V0-A
  values[AliCMEVarManager::kV0AxH2] = qvec[0]; values[AliCMEVarManager::kV0AyH2] = qvec[1];
  values[AliCMEVarManager::kV0ArpH2] = qvec[2];
  qvec[0]=0.0; qvec[1]=0.0; qvec[2]=0.0;
  GetVzeroRP(event, qvec,1);      // V0-C
  values[AliCMEVarManager::kV0CxH2] = qvec[0]; values[AliCMEVarManager::kV0CyH2] = qvec[1];
  values[AliCMEVarManager::kV0CrpH2] = qvec[2];
  qvec[0]=0.0; qvec[1]=0.0; qvec[2]=0.0;
  GetVzeroRP(event, qvec,2);      // V0-A and V0-C combined
  values[AliCMEVarManager::kV0ACxH2] = qvec[0]; values[AliCMEVarManager::kV0ACyH2] = qvec[1];
  values[AliCMEVarManager::kV0ACrpH2] = qvec[2];
  // VZERO event plane resolution
  values[AliCMEVarManager::kV0ArpResH2] = 1.0;
  values[AliCMEVarManager::kV0CrpResH2] = 1.0;
  values[AliCMEVarManager::kV0ACrpResH2] = 1.0;
  // Q vector components correlations
  values[AliCMEVarManager::kV0XaXcH2] = values[AliCMEVarManager::kV0AxH2]*values[AliCMEVarManager::kV0CxH2];
  values[AliCMEVarManager::kV0XaYaH2] = values[AliCMEVarManager::kV0AxH2]*values[AliCMEVarManager::kV0AyH2];
  values[AliCMEVarManager::kV0XaYcH2] = values[AliCMEVarManager::kV0AxH2]*values[AliCMEVarManager::kV0CyH2];
  values[AliCMEVarManager::kV0YaXcH2] = values[AliCMEVarManager::kV0AyH2]*values[AliCMEVarManager::kV0CxH2];
  values[AliCMEVarManager::kV0YaYcH2] = values[AliCMEVarManager::kV0AyH2]*values[AliCMEVarManager::kV0CyH2];
  values[AliCMEVarManager::kV0XcYcH2] = values[AliCMEVarManager::kV0CxH2]*values[AliCMEVarManager::kV0CyH2];


  // event plane differences used for EP resolution calculation
  values[AliCMEVarManager::kV0ATPCDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kV0ArpH2] -
								     values[AliCMEVarManager::kTPCrpH2]) );

  values[AliCMEVarManager::kV0CTPCDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kV0CrpH2] -
								     values[AliCMEVarManager::kTPCrpH2]) );

  values[AliCMEVarManager::kV0AV0CDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kV0ArpH2] -
								     values[AliCMEVarManager::kV0CrpH2]) );

  values[AliCMEVarManager::kv0ATPCDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
								     values[AliCMEVarManager::kTPCrpH2]) );

  values[AliCMEVarManager::kv0CTPCDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kv0CrpH2] -
								     values[AliCMEVarManager::kTPCrpH2]) );

  values[AliCMEVarManager::kv0Av0CDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
								     values[AliCMEVarManager::kv0CrpH2]) );

  values[AliCMEVarManager::kv0Av0C0DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
								     values[AliCMEVarManager::kv0C0rpH2]) );

  values[AliCMEVarManager::kv0Av0C3DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
								     values[AliCMEVarManager::kv0C3rpH2]) );

  values[AliCMEVarManager::kv0Cv0A0DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0CrpH2] -
								     values[AliCMEVarManager::kv0A0rpH2]) );

  values[AliCMEVarManager::kv0Cv0A3DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0CrpH2] -
								     values[AliCMEVarManager::kv0A3rpH2]) );

  values[AliCMEVarManager::kv0A0v0A3DiffH2] = TMath::Cos( 2.*(values[AliCMEVarManager::kv0A0rpH2] -
								     values[AliCMEVarManager::kv0A3rpH2]) );

  values[AliCMEVarManager::kv0C0v0C3DiffH2] = TMath::Cos( 2.*(values[AliCMEVarManager::kv0C0rpH2] -
								     values[AliCMEVarManager::kv0C3rpH2]) );

  Double_t ZDCqvec[3][2] = {{999., 999.}, {999., 999.}, {999., 999.} };
  GetZDCRP(event, ZDCqvec);

  values[AliCMEVarManager::kZDCArpH1] = TMath::ATan2(ZDCqvec[0][1], ZDCqvec[0][0]);
  values[AliCMEVarManager::kZDCCrpH1] = TMath::ATan2(ZDCqvec[1][1], ZDCqvec[1][0]);
  values[AliCMEVarManager::kZDCACrpH1] = TMath::ATan2(ZDCqvec[2][1], ZDCqvec[2][0]);

  if(TMath::Abs(ZDCqvec[0][0] - 999.) < 1e-10 || TMath::Abs(ZDCqvec[0][1] - 999.) < 1e-10 || TMath::Abs(ZDCqvec[1][0] - 999.) < 1e-10 || TMath::Abs(ZDCqvec[1][1] - 999.) < 1e-10){
    values[AliCMEVarManager::kZDCArpH1] = 999;
    values[AliCMEVarManager::kZDCCrpH1] = 999;
    values[AliCMEVarManager::kZDCACrpH1] = 999;
  }



  values[AliCMEVarManager::kv0ZDCrpRes] = cos(2*(values[AliCMEVarManager::kZDCArpH1] - values[AliCMEVarManager::kv0ArpH2]));
  values[AliCMEVarManager::kZDCrpResH1] = cos(values[AliCMEVarManager::kZDCArpH1] - values[AliCMEVarManager::kZDCCrpH1]);


}

inline void AliCMEVarManager::FillVarESDEvent(const AliESDEvent *event, Float_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //

  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  Double_t centralityF=-1;
  Double_t centralitySPD=-1;
  Double_t centralityV0A = -1;
  Double_t centralityV0C = -1;
  Double_t centralityZNA = -1;
  AliCentrality *esdCentrality = const_cast<AliESDEvent*>(event)->GetCentrality();
  //AliMultSelection *MultSelection = (AliMultSelection * ) const_cast<AliESDEvent*>(event)->FindListObject("MultSelection");
  //if (MultSelection) centralityF = MultSelection->GetMultiplicityPercentile("V0M");
  //if (MultSelection) centralitySPD = MultSelection->GetMultiplicityPercentile("CL1");
  if (esdCentrality) centralityF = esdCentrality->GetCentralityPercentile("V0M");
  if (esdCentrality) centralitySPD = esdCentrality->GetCentralityPercentile("CL1");
  if (esdCentrality) centralityV0A = esdCentrality->GetCentralityPercentile("V0A");
  if (esdCentrality) centralityV0C = esdCentrality->GetCentralityPercentile("V0C");
  if (esdCentrality) centralityZNA = esdCentrality->GetCentralityPercentile("ZNA");

  // Fill AliESDEvent interface specific information
  const AliESDVertex *primVtx = event->GetPrimaryVertex();
  values[AliCMEVarManager::kXRes]       = primVtx->GetXRes();
  values[AliCMEVarManager::kYRes]       = primVtx->GetYRes();
  values[AliCMEVarManager::kZRes]       = primVtx->GetZRes();
  values[AliCMEVarManager::kCentrality] = centralityF;
  values[AliCMEVarManager::kCentralitySPD] = centralitySPD;
  values[AliCMEVarManager::kCentralityV0A] = centralityV0A;
  values[AliCMEVarManager::kCentralityV0C] = centralityV0C;
  values[AliCMEVarManager::kCentralityZNA] = centralityZNA;
  values[AliCMEVarManager::kCentralityV0mSPD] = centralityF-centralitySPD;


  const AliESDVertex *vtxTPC = event->GetPrimaryVertexTPC();
  values[AliCMEVarManager::kNVtxContribTPC] = (vtxTPC ? vtxTPC->GetNContributors() : 0);

  // Event multiplicity estimators
  Int_t nTrSPD05=0; Int_t nTrITSTPC05=0; Int_t nTrITSSA05=0;
  nTrSPD05    = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTracklets, 0.5);
  nTrITSTPC05 = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSTPC, 0.5);
  nTrITSSA05  = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSSA, 0.5);
  values[AliCMEVarManager::kNaccTrckltsEsd05] = nTrSPD05;
  values[AliCMEVarManager::kNaccItsTpcEsd05] = nTrITSTPC05;
  values[AliCMEVarManager::kNaccItsPureEsd05] = nTrITSSA05;

  Int_t nTrSPD10=0; Int_t nTrITSTPC10=0; Int_t nTrITSSA10=0;
  nTrSPD10    = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTracklets, 1.0);
  nTrITSTPC10 = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSTPC, 1.0);
  nTrITSSA10  = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSSA, 1.0);
  values[AliCMEVarManager::kNaccTrckltsEsd10] = nTrSPD10;
  values[AliCMEVarManager::kNaccItsTpcEsd10] = nTrITSTPC10;
  values[AliCMEVarManager::kNaccItsPureEsd10] = nTrITSSA10;

  Int_t nTrSPD16=0; Int_t nTrITSTPC16=0; Int_t nTrITSSA16=0;
  nTrSPD16    = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTracklets, 1.6);
  nTrITSTPC16 = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSTPC, 1.6);
  nTrITSSA16  = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSSA, 1.6);
  values[AliCMEVarManager::kNaccTrckltsEsd16] = nTrSPD16;
  values[AliCMEVarManager::kNaccItsTpcEsd16] = nTrITSTPC16;
  values[AliCMEVarManager::kNaccItsPureEsd16] = nTrITSSA16;

}

inline void AliCMEVarManager::FillVarAODEvent(const AliAODEvent *event, Float_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //

  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  // Fill AliAODEvent interface specific information
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  assert(header&&"Not a standard AOD");

  Double_t centralityF=-1;
  Double_t centralitySPD=-1;
  Double_t centralityV0A = -1;
  Double_t centralityV0C = -1;
  Double_t centralityZNA = -1;
  AliCentrality *aodCentrality = header->GetCentralityP();
  if (aodCentrality) centralityF = aodCentrality->GetCentralityPercentile("V0M");
  if (aodCentrality) centralitySPD = aodCentrality->GetCentralityPercentile("CL1");
  if (aodCentrality) centralityV0A = aodCentrality->GetCentralityPercentile("V0A");
  if (aodCentrality) centralityV0C = aodCentrality->GetCentralityPercentile("V0C");
  if (aodCentrality) centralityZNA = aodCentrality->GetCentralityPercentile("ZNA");
  values[AliCMEVarManager::kCentrality] = centralityF;
  values[AliCMEVarManager::kCentralitySPD] = centralitySPD;
  values[AliCMEVarManager::kCentralityV0A] = centralityV0A;
  values[AliCMEVarManager::kCentralityV0C] = centralityV0C;
  values[AliCMEVarManager::kCentralityZNA] = centralityZNA;
  values[AliCMEVarManager::kCentralityV0mSPD] = centralityF-centralitySPD;

  values[AliCMEVarManager::kRefMult]        = header->GetRefMultiplicity();        // similar to Ntrk
  values[AliCMEVarManager::kRefMultTPConly] = header->GetTPConlyRefMultiplicity(); // similar to Nacc

  ///////////////////////////////////////////
  //////////// NANO AODs ////////////////////
  ///////////////////////////////////////////

  // (w/o AliCentrality branch), VOM centrality should be stored in the header
  if(!header->GetCentralityP())
    values[AliCMEVarManager::kCentrality] = header->GetCentrality();
  // (w/o AliEventPlane branch) tpc event plane stuff stored in the header
  if(!header->GetEventplaneP()) {

    //    values[AliCMEVarManager::kNTrk] = header->GetRefMultiplicity();    // overwritten datamembers in "our" nanoAODs
    //    values[AliCMEVarManager::kNacc] = header->GetRefMultiplicityPos(); // overwritten datamembers in "our" nanoAODs

    TVector2 qvec;
    // TPC
    qvec.Set(header->GetEventplaneQx(), header->GetEventplaneQy());
    values[AliCMEVarManager::kTPCxH2uc]   = qvec.X();
    values[AliCMEVarManager::kTPCyH2uc]   = qvec.Y();
    values[AliCMEVarManager::kTPCmagH2uc] = qvec.Mod();
    values[AliCMEVarManager::kTPCrpH2uc]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;

    // VZERO
    AliEventplane ep2;
    // get event plane corrections from the VZERO EP selection task
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliVZEROEPSelectionTask *eptask = dynamic_cast<AliVZEROEPSelectionTask *>(man->GetTask("AliVZEROEPSelectionTask"));
    if(eptask) eptask->SetEventplaneParams(&ep2,centralityF);
    else printf("no VZERO event plane selection task added! \n");

    Double_t qx = 0, qy = 0;
    ep2.CalculateVZEROEventPlane(event,10, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0ACrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliCMEVarManager::kv0ACxH2]   = qvec.X();
    values[AliCMEVarManager::kv0ACyH2]   = qvec.Y();
    values[AliCMEVarManager::kv0ACmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 8, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0ArpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliCMEVarManager::kv0AxH2]   = qvec.X();
    values[AliCMEVarManager::kv0AyH2]   = qvec.Y();
    values[AliCMEVarManager::kv0AmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 9, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0CrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliCMEVarManager::kv0CxH2]   = qvec.X();
    values[AliCMEVarManager::kv0CyH2]   = qvec.Y();
    values[AliCMEVarManager::kv0CmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 0, 0, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0C0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep2.CalculateVZEROEventPlane(event, 3, 3, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0C3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep2.CalculateVZEROEventPlane(event, 4, 4, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0A0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep2.CalculateVZEROEventPlane(event, 7, 7, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliCMEVarManager::kv0A3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;

  }

  const AliAODVertex *vtxtpc = GetVertex(event, AliAODVertex::kMainTPC);
  values[AliCMEVarManager::kNVtxContribTPC] = (vtxtpc ? vtxtpc->GetNContributors() : 0);

}

inline void AliCMEVarManager::FillVarMCEvent(const AliMCEvent *event, Float_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //

  // Fill common AliVEvent interface information
  //  FillVarVEvent(event, values);
  const AliVVertex* vtx = event->GetPrimaryVertex();
  values[AliCMEVarManager::kXvPrim]       = (vtx ? vtx->GetX() : 0.0);
  values[AliCMEVarManager::kYvPrim]       = (vtx ? vtx->GetY() : 0.0);
  values[AliCMEVarManager::kZvPrim]       = (vtx ? vtx->GetZ() : 0.0);
}

inline void AliCMEVarManager::FillVarTPCEventPlane(const AliEventplane *evplane, Float_t * const values)
{
  //
  // Fill TPC event plane information after correction
  //
  if(evplane) {
    TVector2 *qcorr  = const_cast<AliEventplane *>(evplane)->GetQVector();  // This is the "corrected" Q-Vector
    TVector2 *qcsub1 = const_cast<AliEventplane *>(evplane)->GetQsub1();
    TVector2 *qcsub2 = const_cast<AliEventplane *>(evplane)->GetQsub2();
    if(qcorr) {
      values[AliCMEVarManager::kTPCxH2]   = qcorr->X();
      values[AliCMEVarManager::kTPCyH2]   = qcorr->Y();
      values[AliCMEVarManager::kTPCmagH2] = qcorr->Mod();
      values[AliCMEVarManager::kTPCrpH2]  = TVector2::Phi_mpi_pi(qcorr->Phi())/2;
      // detector effects
      values[AliCMEVarManager::kCosTPCrpH2]     = TMath::Cos( 2.* values[AliCMEVarManager::kTPCrpH2] );
      values[AliCMEVarManager::kSinTPCrpH2]     = TMath::Sin( 2.* values[AliCMEVarManager::kTPCrpH2] );

      // correlations for event plane resoultion
      values[AliCMEVarManager::kv0ATPCDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
									 values[AliCMEVarManager::kTPCrpH2]) );
      values[AliCMEVarManager::kv0CTPCDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kv0CrpH2] -
									 values[AliCMEVarManager::kTPCrpH2]) );
      values[AliCMEVarManager::kv0Av0CDiffH2]   = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
									 values[AliCMEVarManager::kv0CrpH2]) );
      values[AliCMEVarManager::kv0Av0C0DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
									 values[AliCMEVarManager::kv0C0rpH2]) );
      values[AliCMEVarManager::kv0Av0C3DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0ArpH2] -
									 values[AliCMEVarManager::kv0C3rpH2]) );
      values[AliCMEVarManager::kv0Cv0A0DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0CrpH2] -
									 values[AliCMEVarManager::kv0A0rpH2]) );
      values[AliCMEVarManager::kv0Cv0A3DiffH2]  = TMath::Cos( 2.*(values[AliCMEVarManager::kv0CrpH2] -
									 values[AliCMEVarManager::kv0A3rpH2]) );
      values[AliCMEVarManager::kv0A0v0A3DiffH2] = TMath::Cos( 2.*(values[AliCMEVarManager::kv0A0rpH2] -
									 values[AliCMEVarManager::kv0A3rpH2]) );
      values[AliCMEVarManager::kv0C0v0C3DiffH2] = TMath::Cos( 2.*(values[AliCMEVarManager::kv0C0rpH2] -
									 values[AliCMEVarManager::kv0C3rpH2]) );
    }
    if(qcsub1 && qcsub2) {
      values[AliCMEVarManager::kTPCsub1xH2]   = qcsub1->X();
      values[AliCMEVarManager::kTPCsub1yH2]   = qcsub1->Y();
      values[AliCMEVarManager::kTPCsub1rpH2]  = TVector2::Phi_mpi_pi(qcsub1->Phi())/2;

      values[AliCMEVarManager::kTPCsub2xH2]   = qcsub2->X();
      values[AliCMEVarManager::kTPCsub2yH2]   = qcsub2->Y();
      values[AliCMEVarManager::kTPCsub2rpH2]  = TVector2::Phi_mpi_pi(qcsub2->Phi())/2;

      values[AliCMEVarManager::kTPCsub12DiffH2] = TMath::Cos( 2.*(values[AliCMEVarManager::kTPCsub1rpH2] -
									 values[AliCMEVarManager::kTPCsub2rpH2]) );
      values[AliCMEVarManager::kTPCsub12DiffH2Sin] = TMath::Sin( 2.*(values[AliCMEVarManager::kTPCsub1rpH2] -
									    values[AliCMEVarManager::kTPCsub2rpH2]) );
    }
  }
}

inline void AliCMEVarManager::InitESDpid(Int_t type)
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

inline void AliCMEVarManager::InitAODpidUtil(Int_t type)
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


inline void AliCMEVarManager::InitEstimatorAvg(const Char_t* filename) //Not Grid compatible
{
  //
  // initialize the profile histograms neccessary for the correction of the multiplicity estimators in pp collisions
  //

  const Char_t* estimatorNames[9] = {"SPDmult05","SPDmult10","SPDmult16",
				     "ITSTPC05", "ITSTPC10", "ITSTPC16",
				     "ITSSA05",  "ITSSA10",  "ITSSA16"};
  const Char_t* periodNames[6] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e", "LHC13b", "LHC13c"};
  TFile* file=TFile::Open(filename);
  if(!file) return;

  for(Int_t ip=0; ip<6; ++ip) {
    for(Int_t ie=0; ie<9; ++ie) {
      fgMultEstimatorAvg[ip][ie] = (TProfile*)(file->Get(Form("%s_%s",estimatorNames[ie],periodNames[ip]))->Clone(Form("%s_%s_clone",estimatorNames[ie],periodNames[ip])));
    }
  }
}


inline void AliCMEVarManager::InitEstimatorObjArrayAvg(const TObjArray* array) //Grid compatible
{
  //
  // initialize the profile histograms neccessary for the correction of the multiplicity estimators in pp collisions
  // SPDmult05 does not exist yet for Pass4 AODs
  // ITS correction maps do not exist yet for Pass4 AODs

  const Char_t* estimatorNames[9] = {"SPDmult05","SPDmult10","SPDmult16",
				     "ITSTPC05", "ITSTPC10", "ITSTPC16",
				     "ITSSA05",  "ITSSA10",  "ITSSA16"};
  const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};

  TString key;

  Int_t ieTotal = 9;

  for(Int_t ip=0; ip<4; ++ip) {
    for(Int_t ie=0; ie<ieTotal; ++ie) {
      key = Form("%s_%s",estimatorNames[ie],periodNames[ip]);
      if(array->FindObject(key.Data())){
	fgMultEstimatorAvg[ip][ie] = (TProfile*)(array->FindObject(key.Data()))->Clone((key+"_clone").Data());
	continue;
      }
      key += Form("_Pass4_AOD");
      if(array->FindObject(key.Data())){
	printf("ip = %d, ie = %d\n estimator = %s, period = %s\n key = %s\n", ip, ie, estimatorNames[ie], periodNames[ip], key.Data());
        fgMultEstimatorAvg[ip][ie] = (TProfile*)(array->FindObject(key.Data()))->Clone((key+"_clone").Data());
      }
    }
  }
}
inline void AliCMEVarManager::InitTRDpidEffHistograms(const Char_t* filename)
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

inline Double_t AliCMEVarManager::GetSingleLegEff(Float_t * const values) {
  //
  // get the single leg efficiency for a given particle
  //
  if(!fgLegEffMap) return -1.;

  if(fgLegEffMap->IsA()== THnBase::Class()) {
    THnBase *eff = static_cast<THnBase*>(fgLegEffMap);
    Int_t dim=eff->GetNdimensions();
    Int_t *idx=new Int_t[dim];
    for(Int_t idim=0; idim<dim; idim++) {
      UInt_t var = GetValueType(eff->GetAxis(idim)->GetName());
      idx[idim] = eff->GetAxis(idim)->FindBin(values[var]);
      if(idx[idim] < 0 || idx[idim]>eff->GetAxis(idim)->GetNbins()) return 0.0;
    }
    //  printf(" bin content %f+-%f \n",eff->GetBinContent(idx), eff->GetBinError(idx));
    const Double_t ret=(eff->GetBinContent(idx));
    delete [] idx;
    return ret;
  }
  return -1.;
}

inline Double_t AliCMEVarManager::GetPairEff(Float_t * const values) {
  //
  // get the pair efficiency for given pair kinematics
  //
  if(!fgPairEffMap) return -1.;

  if(fgPairEffMap->IsA()== THnBase::Class()) {
    THnBase *eff = static_cast<THnBase*>(fgPairEffMap);
    Int_t dim=eff->GetNdimensions();
    Int_t *idx=new Int_t[dim];
    for(Int_t idim=0; idim<dim; idim++) {
      UInt_t var = GetValueType(eff->GetAxis(idim)->GetName());
    idx[idim] = eff->GetAxis(idim)->FindBin(values[var]);
    if(idx[idim] < 0 || idx[idim]>eff->GetAxis(idim)->GetNbins()) return 0.0;
    }
    //  printf(" bin content %f+-%f \n",eff->GetBinContent(idx), eff->GetBinError(idx));
    const Double_t ret=(eff->GetBinContent(idx));
    delete [] idx;
    return ret;
  }
  if(fgPairEffMap->IsA()== TSpline3::Class()) {
    TSpline3 *eff = static_cast<TSpline3*>(fgPairEffMap);
    if(!eff->GetHistogram()) { printf("no histogram added to the spline\n"); return -1.;}
    UInt_t var = GetValueType(eff->GetHistogram()->GetXaxis()->GetName());
    //printf(" bin content %f \n",eff->Eval(values[var]) );
    return (eff->Eval(values[var]));
  }

  return -1.;
}


inline void AliCMEVarManager::InitVZEROCalibrationHistograms(Int_t runNo) {
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


inline void AliCMEVarManager::InitVZERORecenteringHistograms(Int_t runNo) {
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

inline void AliCMEVarManager::InitZDCRecenteringHistograms(Int_t runNo) {

  //initialize only once
  if(fgZDCRecentering[0][0]) return;

  for(Int_t i=0; i<2; ++i)
    for(Int_t j=0; j<2; ++j)
      if(fgZDCRecentering[i][j]) {
        delete fgZDCRecentering[i][j];
        fgZDCRecentering[i][j] = 0x0;
      }

  TFile* file=TFile::Open(fgZDCRecenteringFile.Data());
  if(!file) return;


  fgZDCRecentering[0][0] = (TProfile3D*)file->Get(Form("RUN%06d_QxA_Recent", runNo));
  fgZDCRecentering[0][1] = (TProfile3D*)file->Get(Form("RUN%06d_QyA_Recent", runNo));
  fgZDCRecentering[1][0] = (TProfile3D*)file->Get(Form("RUN%06d_QxC_Recent", runNo));
  fgZDCRecentering[1][1] = (TProfile3D*)file->Get(Form("RUN%06d_QyC_Recent", runNo));
  fgZDCRecentering[2][0] = (TProfile3D*)file->Get(Form("RUN%06d_QxAC_Recent", runNo));
  fgZDCRecentering[2][1] = (TProfile3D*)file->Get(Form("RUN%06d_QyAC_Recent", runNo));


  if (fgZDCRecentering[0][0]) fgZDCRecentering[0][0]->SetDirectory(0x0);
  if (fgZDCRecentering[0][1]) fgZDCRecentering[0][1]->SetDirectory(0x0);
  if (fgZDCRecentering[1][0]) fgZDCRecentering[1][0]->SetDirectory(0x0);
  if (fgZDCRecentering[1][1]) fgZDCRecentering[1][1]->SetDirectory(0x0);
  if (fgZDCRecentering[2][0]) fgZDCRecentering[2][0]->SetDirectory(0x0);
  if (fgZDCRecentering[2][1]) fgZDCRecentering[2][1]->SetDirectory(0x0);

  delete file;

}


inline Double_t AliCMEVarManager::GetTRDpidEfficiency(Int_t runNo, Double_t centrality,
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


inline void AliCMEVarManager::SetEvent(AliVEvent * const ev)
{

  fgEvent = ev;
  if (fgKFVertex) delete fgKFVertex;
  fgKFVertex=0x0;
  if (!ev) return;
  if (ev->GetPrimaryVertex()) fgKFVertex=new AliKFVertex(*ev->GetPrimaryVertex());

  for (Int_t i=0; i<AliCMEVarManager::kNMaxValues;++i) fgData[i]=0.;
  AliCMEVarManager::Fill(fgEvent, fgData);
}

inline void AliCMEVarManager::SetEventData(const Double_t data[AliCMEVarManager::kNMaxValues])
{
  for (Int_t i=0; i<kNMaxValues;++i) fgData[i]=0.;
  for (Int_t i=kPairMax; i<kNMaxValues;++i) fgData[i]=data[i];
}


//______________________________________________________________________________
inline Bool_t AliCMEVarManager::GetDCA(const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0)
{
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    // the covariance matrix is not stored in case of AliAODTrack::kIsDCA
    return kTRUE;
  }

  Bool_t ok=kFALSE;
  if(fgEvent) {
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);

    Float_t xstart = etp.GetX();
    if(xstart>3.) {
      d0z0[0]=-999.;
      d0z0[1]=-999.;
      return kFALSE;
    }

    AliAODVertex *vtx =(AliAODVertex*)(fgEvent->GetPrimaryVertex());
    Double_t fBzkG = fgEvent->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}

inline void AliCMEVarManager::SetTPCEventPlane(AliEventplane *const evplane)
{

  fgTPCEventPlane = evplane;
  FillVarTPCEventPlane(evplane,fgData);
  //  for (Int_t i=0; i<AliCMEVarManager::kNMaxValues;++i) fgData[i]=0.;
  //  AliCMEVarManager::Fill(fgEvent, fgData);
}


//_________________________________________________________________
inline void AliCMEVarManager::GetVzeroRP(const AliVEvent* event, Double_t* qvec, Int_t sideOption) {
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
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(aodEv->GetHeader());
    assert(header&&"Not a standard AOD");
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
inline void AliCMEVarManager::GetZDCRP(const AliVEvent* event, Double_t qvec[][2]) {

  //
  // Get the reaction plane from the ZDC detector for first harmonic
  //
  //  Q{x,y} = SUM{ri(x,y)*Ei} / SUM{Ei}
  //

  const Int_t   nZDCSides  = 2;
  const Int_t   nZDCplanes = 3;
  const Int_t   Aside = 0, Cside = 1, ACside = 2;
  const Int_t   nZDCTowers = 4;// number of ZDCtowers
  const Double_t ZDCTowerCenters[nZDCTowers][2] = { {-1.75, -1.75}, { 1.75, -1.75},
                                                    {-1.75,  1.75}, { 1.75,  1.75} };

  Double_t   *ZDCTEnergy[nZDCSides]; //reco E in 5 ZDC sectors - high gain chain
  Double_t    qvecNUM[nZDCplanes][2];
  Double_t    qvecDEN[nZDCplanes];
  memset(   qvecNUM,    0,     sizeof(qvecNUM));  //format
  memset(qvecDEN,     0,     sizeof(qvecDEN));  //format

  Double_t TPCRefMulti = 999, vtxX = 999, vtxY = 999;
  Int_t multiBin = 0, vtxXBin = 0, vtxYBin = 0;
  Double_t recentdim[3][3] = { { 50, 0, 2500},   //multiplicity nbin, min, max
                               { 20, 0.04, 0.08},   //    vertex x nbin, min, max
                               { 20, 0.25, 0.29} }; //    vertex y nbin, min, max

  if(!event->GetZDCData()) return;
  AliVZDC* aliZDC = event->GetZDCData();
  ZDCTEnergy[Aside] = (Double_t *)aliZDC -> GetZNATowerEnergy();
  ZDCTEnergy[Cside] = (Double_t *)aliZDC -> GetZNCTowerEnergy();

  for(int j = 0;  j < nZDCSides   ; j++){
    for(int k = 0;   k < nZDCTowers ; k++){
      qvecNUM[j][0] += ZDCTowerCenters[k][0]*ZDCTEnergy[j][k+1]; //    zdcQ += xE
      qvecNUM[j][1] += ZDCTowerCenters[k][1]*ZDCTEnergy[j][k+1]; //    zdcQ += yE
      qvecDEN[j]    += ZDCTEnergy[j][k+1];                   // zdcQsum +=  E

    }
    if(j == Aside){
      qvecNUM[j][0] = -qvecNUM[j][0];
    }

    if(j == Cside){
      qvecNUM[j][0] = -qvecNUM[j][0];
      qvecNUM[j][1] = -qvecNUM[j][1];
    }


    qvecNUM[ACside][0] += qvecNUM[j][0];
    qvecNUM[ACside][1] += qvecNUM[j][1];
    qvecDEN[ACside] += qvecDEN[j];

  }

  for(int j = 0; j < nZDCplanes; j++){
    if(qvecDEN[j] != 0){
      qvec[j][0] = (qvecNUM[j][0] / qvecDEN[j]);
      qvec[j][1] = (qvecNUM[j][1] / qvecDEN[j]);
    }
    else if(qvecDEN[j] == 0) {
      qvec[j][0] = 999;
      qvec[j][1] = 999;
    }

  }

  if(fgZDCRecentering[0][0]){
    const AliAODEvent* aodEv = static_cast<const AliAODEvent*>(event);
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(aodEv->GetHeader());
    if(!header) return;
    TPCRefMulti = header -> GetTPConlyRefMultiplicity();

    const AliVVertex *primVtx = event->GetPrimaryVertex();
    if(!primVtx) return;
    vtxX = primVtx->GetX();
    vtxY = primVtx->GetY();

    multiBin = (Int_t)((TPCRefMulti-recentdim[0][1])*recentdim[0][0] / (recentdim[0][2] - recentdim[0][1])) + 1;
    vtxXBin  = (Int_t)((vtxX-recentdim[1][1])*recentdim[1][0] / (recentdim[1][2] - recentdim[1][1])) + 1;
    vtxYBin  = (Int_t)((vtxY-recentdim[2][1])*recentdim[2][0] / (recentdim[2][2] - recentdim[2][1])) + 1;

    for(int j = 0; j < nZDCplanes; j++)
      if(qvecDEN[j] != 0){
        qvec[j][0] -= fgZDCRecentering[j][0] -> GetBinContent(multiBin, vtxXBin, vtxYBin);
        qvec[j][1] -= fgZDCRecentering[j][1] -> GetBinContent(multiBin, vtxXBin, vtxYBin);
      }
  }

}



//______________________________________________________________________________
inline AliAODVertex* AliCMEVarManager::GetVertex(const AliAODEvent* event, AliAODVertex::AODVtx_t vtype) {
  // Get vertex
  Int_t nVertices=event->GetNumberOfVertices();
  for(Int_t iVert=0; iVert<nVertices; iVert++){
    AliAODVertex *v=event->GetVertex(iVert);
    //    printf(" vtx %d  contrib %d  daughters %d \n ",v->GetType(),v->GetNContributors(), v->GetNDaughters());
    if(v->GetType()==vtype) return v;
  }
  return 0;
}



//_________________________________________________________________
inline THnF* AliCMEVarManager::CreateTHnF( const Char_t* name, const Char_t* title, Int_t nDimensions,TAxis* binLimits){
  //
  // create a multi-dimensional histogram THnF with equal or variable bin widths
  //
  if(!binLimits) return 0x0;
  TString hname = name;

  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");


  Double_t* xmin = new Double_t[nDimensions];
  Double_t* xmax = new Double_t[nDimensions];
  Int_t* nBins = new Int_t[nDimensions];
  for(Int_t idim=0;idim<nDimensions;++idim) {
    nBins[idim] = binLimits[idim].GetNbins();
    xmin[idim] = binLimits[idim].GetBinLowEdge(1);
    xmax[idim] = binLimits[idim].GetBinUpEdge(nBins[idim]);
  }

  THnF* h=new THnF(hname.Data(),title,nDimensions,nBins,xmin,xmax);
  for(Int_t idim=0;idim<nDimensions;++idim) {
    TAxis* axis=h->GetAxis(idim);
    *axis=TAxis(binLimits[idim]);
    if(arr->GetEntries()>(idim+1)) axis->SetTitle(arr->At(idim+1)->GetName());
  }

  h->Sumw2();

  delete [] xmin;
  delete [] xmax;
  delete [] nBins;

  return h;
}



/*
inline void AliCMEVarManager::FillValues(const TParticle *particle, Double_t *values)
{
  //
  // Fill track information available for histogramming into an array
  //

  // Fill TParticle interface information
  values[AliCMEVarManager::kPx]     = particle->Px();
  values[AliCMEVarManager::kPy]     = particle->Py();
  values[AliCMEVarManager::kPz]     = particle->Pz();
  values[AliCMEVarManager::kPt]     = particle->Pt();
  values[AliCMEVarManager::kP]      = particle->P();

  values[AliCMEVarManager::kXv]     = particle->Vx();
  values[AliCMEVarManager::kYv]     = particle->Vy();
  values[AliCMEVarManager::kZv]     = particle->Vz();

  values[AliCMEVarManager::kOneOverPt] = 1./particle->Pt();
  values[AliCMEVarManager::kPhi]    = particle->Phi();
  values[AliCMEVarManager::kTheta]  =
  values[AliCMEVarManager::kEta]    = particle->Eta();
  values[AliCMEVarManager::kY]      =

  values[AliCMEVarManager::kE]      = particle->Energy();
  values[AliCMEVarManager::kM]      = particle->GetMass();

  values[AliCMEVarManager::kCharge] = particle->GetPDG()->Charge()/3; // uggly

}*/

#endif
