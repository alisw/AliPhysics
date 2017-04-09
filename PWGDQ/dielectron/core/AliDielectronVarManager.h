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
#include <TProfile3D.h>
#include <TH3D.h>
#include <THnBase.h>
#include <TSpline.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TKey.h>
#include <TBits.h>
#include <TRandom3.h>

#include <AliLog.h>

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
#include "AliMultSelection.h"
#include <AliAODpidUtil.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliESDtrackCuts.h>

#include "AliDielectronPair.h"
#include "AliDielectronMC.h"
#include "AliDielectronPID.h"
#include "AliDielectronHelper.h"
#include "AliDielectronQnEPcorrection.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVZEROEPSelectionTask.h"

#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"

#include "AliAODMCHeader.h"
#include "AliTRDgeometry.h"
#include "assert.h"

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
    kPtMC,                   // MC transverse momentum
    kPtSq,                   // transverse momentum squared
    kP,                      // momentum
    kPMC,                    // MC momentum
    kXv,                     // vertex position in x
    kYv,                     // vertex position in y
    kZv,                     // vertex position in z
    kOneOverPt,              // 1/pt
    kPhi,                    // phi angle
    kPhiMC,                  // MC phi angle
    kTheta,                  // theta angle
    kEta,                    // pseudo-rapidity
    kEtaMC,                  // MC pseudo-rapidity
    kY,                      // rapidity
    kE,                      // energy
    kM,                      // mass
    kMCorr,                  // mass, corrected (for photons)
    kMMC,                    // MC mass
    kCharge,                 // charge
    kNclsITS,                // number of clusters assigned in the ITS
    kITSFakeFlag,            // ITS fake flag
    kITSchi2Cl,              // chi2/cl in the ITS
    kNclsTPC,                // number of clusters assigned in the TPC
    kNclsSTPC,               // number of shared clusters assigned in the TPC
    kNclsSFracTPC,           // fraction of shared clusters assigned in the TPC
    kNclsSITS,                // number of shared clusters assigned in the ITS
    kNclsSFracITS,           // fraction of shared clusters assigned in the ITS
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
    kTPCclsIRO,              // TPC clusters inner read out
    kTPCclsORO,              // TPC clusters outer read out

    kTPCActiveLength,        // TPC length in active volume
    kTPCGeomLength,          // TPC pT dependent geometrical length

    kTrackStatus,            // track status bits
    kFilterBit,              // AOD filter bits

    kNclsTRD,                // number of clusters assigned in the TRD
    kTRDntracklets,          // number of TRD tracklets used for tracking/PID TODO: correct getter
    kTRDpidQuality,          // number of TRD tracklets used for PID
    kTRDchi2,                // chi2 in TRD
    kTRDchi2Trklt,	         // chi2/trklts in TRD
    kTRDprobEle,             // TRD electron pid probability
    kTRDprobPio,             // TRD pion pid probability
    kTRDprob2DEle,           // TRD electron pid probability 2D LQ
    kTRDprob2DPio,           // TRD pion pid probability 2D LQ
    kTRDprob2DPro,           // TRD proton pid probability 2D LQ
    kTRDprob3DEle,           // TRD electron pid probability 3D LQ
    kTRDprob3DPio,           // TRD pion pid probability 3D LQ
    kTRDprob3DPro,           // TRD proton pid probability 3D LQ
    kTRDprob7DEle,           // TRD electron pid probability 7D LQ
    kTRDprob7DPio,           // TRD pion pid probability 7D LQ
    kTRDprob7DPro,           // TRD proton pid probability 7D LQ
    kTRDphi,                 // Phi angle of the track at the entrance of the TRD
    kTRDpidEffLeg,           // TRD pid efficiency from conversion electrons
    kTRDsignal,              // TRD signal
    kTRDeta,                 // eta of the track at the entrance of the TRD
    kInTRDacceptance,        // in TRD acceptance

    kImpactParXY,            // Impact parameter in XY plane
    kImpactParZ,             // Impact parameter in Z
    kImpactParXYsigma,       // Impact parameter in XY plane normalized to resolution
    kImpactParZsigma,        // Impact parameter in Z normalized to resolution
    kTrackLength,            // Track length
    kDistPrimToSecVtxXYMC,      // Distance of secondary vertex to primary vertex  in the XY plane
    kDistPrimToSecVtxZMC,       // Distance of secondary vertex to primary vertex in Z


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
    kTRDonlineNTracklets,    //TRD online: number of contributing online tracklets


    kParticleMax,             //
    // TODO: kRNClusters ??
    // AliDielectronPair specific variables
    kChi2NDF = kParticleMax, // Chi^2/NDF
    kDecayLength,            // decay length
    kR,                      // distance to the origin
    kOpeningAngle,           // opening angle
    kOpeningAngleCorr,        // opening angle, corrected (for photons)
    kOpeningAngleXY,           // opening angle at in XY direction
    kOpeningAngleRZ,           // opening angle at in RZ direction
    kTriangularConversionCut, // triangular cut on opening angle and kPhivPair
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
    kDeltaCotTheta,          // difference of cotangens of theta of daughters

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
    kDeltaPhiChargeOrdered,  // Absolute value of Delta Phi for the legs
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
    kXvPrimMCtruth,          // MC true prim vertex, so that it is available also for reco tracks
    kYvPrimMCtruth,          // MC true prim vertex, so that it is available also for reco tracks
    kZvPrimMCtruth,          // MC true prim vertex, so that it is available also for reco tracks
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
    kTPCrpH2uc,                 // TPC reaction plane angle of the Q vector for 2nd harmonic (uncorrected) -- corrected if the QnCorrections framework est. 2016 is used
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
    // Beginning of Eventplane variables from Qn Framework
    // Eventplanes for 2nd harmonic from QnCorrections framework est. 2016
    kQnTPCrpH2,                // TPC eventplane from QnCorrections framework
    kQnTPCxH2,
    kQnTPCyH2,
    kQnTPCaSiderpH2,                // TPC A-Side eventplane from QnCorrections framework
    kQnTPCaSidexH2,
    kQnTPCaSideyH2,
    kQnTPCcSiderpH2,                // TPC C-Side eventplane from QnCorrections framework
    kQnTPCcSidexH2,
    kQnTPCcSideyH2,

    kQnV0ArpH2,                // V0A eventplane from QnCorrections framework
    kQnV0AxH2,
    kQnV0AyH2,
    kQnV0CrpH2,                // V0C eventplane from QnCorrections framework
    kQnV0CxH2,
    kQnV0CyH2,
    kQnV0rpH2,                // V0 combined eventplane from QnCorrections framework
    kQnV0xH2,
    kQnV0yH2,

    kQnSPDrpH2,                // SPD eventplane from QnCorrections framework
    kQnSPDxH2,
    kQnSPDyH2,
    kQnFMDArpH2,               // FMDA eventplane from QnCorrections framework
    kQnFMDAxH2,
    kQnFMDAyH2,
    kQnFMDCrpH2,               // FMDA eventplane from QnCorrections framework
    kQnFMDCxH2,
    kQnFMDCyH2,
    // Average Eventplane differences for 2nd harmonics from QnCorrections framework est. 2016 - as input for 3 sub-detector method
    // Returns cos(2(psi_DetA-psi_DetB))
    kQnDiffTPC_V0A,
    kQnDiffTPC_V0C,
    kQnDiffTPC_SPD,
    kQnDiffTPC_FMDA,
    kQnDiffTPC_FMDC,
    kQnDiffTPCa_V0,
    kQnDiffTPCa_V0A,
    kQnDiffTPCa_V0C,
    kQnDiffTPCa_TPCc,
    kQnDiffTPCc_V0,
    kQnDiffTPCc_V0A,
    kQnDiffTPCc_V0C,
    kQnDiffV0A_V0C,
    kQnDiffV0A_SPD,
    kQnDiffV0A_FMDA,
    kQnDiffV0A_FMDC,
    kQnDiffV0C_SPD,
    kQnDiffV0C_FMDA,
    kQnDiffV0C_FMDC,
    kQnDiffSPD_FMDA,
    kQnDiffSPD_FMDC,
    kQnDiffFMDA_FMDC,

    // // Average XX YX XY differences for 2nd harmonics from QnCorrections framework est. 2016
    kQnCorrTPCx_V0Ax,
    kQnCorrTPCx_V0Ay,
    kQnCorrTPCy_V0Ax,
    kQnCorrTPCy_V0Ay,
    kQnCorrTPCx_V0Cx,
    kQnCorrTPCx_V0Cy,
    kQnCorrTPCy_V0Cx,
    kQnCorrTPCy_V0Cy,
    kQnCorrTPCx_SPDx,
    kQnCorrTPCx_SPDy,
    kQnCorrTPCy_SPDx,
    kQnCorrTPCy_SPDy,
    kQnCorrTPCx_FMDAx,
    kQnCorrTPCx_FMDAy,
    kQnCorrTPCy_FMDAx,
    kQnCorrTPCy_FMDAy,
    kQnCorrTPCx_FMDCx,
    kQnCorrTPCx_FMDCy,
    kQnCorrTPCy_FMDCx,
    kQnCorrTPCy_FMDCy,
    kQnCorrV0Ax_V0Cx,
    kQnCorrV0Ax_V0Cy,
    kQnCorrV0Ay_V0Cx,
    kQnCorrV0Ay_V0Cy,
    kQnCorrV0Ax_SPDx,
    kQnCorrV0Ax_SPDy,
    kQnCorrV0Ay_SPDx,
    kQnCorrV0Ay_SPDy,
    kQnCorrV0Ax_FMDAx,
    kQnCorrV0Ax_FMDAy,
    kQnCorrV0Ay_FMDAx,
    kQnCorrV0Ay_FMDAy,
    kQnCorrV0Ax_FMDCx,
    kQnCorrV0Ax_FMDCy,
    kQnCorrV0Ay_FMDCx,
    kQnCorrV0Ay_FMDCy,
    kQnCorrV0Cx_FMDAx,
    kQnCorrV0Cx_FMDAy,
    kQnCorrV0Cy_FMDAx,
    kQnCorrV0Cy_FMDAy,
    kQnCorrV0Cx_FMDCx,
    kQnCorrV0Cx_FMDCy,
    kQnCorrV0Cy_FMDCx,
    kQnCorrV0Cy_FMDCy,
    kQnCorrV0Cx_SPDx,
    kQnCorrV0Cx_SPDy,
    kQnCorrV0Cy_SPDx,
    kQnCorrV0Cy_SPDy,
    kQnCorrSPDx_FMDAx,
    kQnCorrSPDx_FMDAy,
    kQnCorrSPDy_FMDAx,
    kQnCorrSPDy_FMDAy,
    kQnCorrSPDx_FMDCx,
    kQnCorrSPDx_FMDCy,
    kQnCorrSPDy_FMDCx,
    kQnCorrSPDy_FMDCy,
    kQnCorrFMDAx_FMDCx,
    kQnCorrFMDAx_FMDCy,
    kQnCorrFMDAy_FMDCx,
    kQnCorrFMDAy_FMDCy,

    // Flow estimators for measured Jpsis

    kQnDeltaPhiTPCrpH2,
    kQnDeltaPhiV0ArpH2,
    kQnDeltaPhiV0CrpH2,
    kQnDeltaPhiV0rpH2,
    kQnDeltaPhiSPDrpH2,
    kQnTPCrpH2FlowV2,
    kQnV0ArpH2FlowV2,
    kQnV0CrpH2FlowV2,
    kQnV0rpH2FlowV2,
    kQnSPDrpH2FlowV2,

    // End of Eventplane variables from Qn Framework

    kNTrk,                   // number of tracks (or tracklets) TODO: ambiguous
    kTracks,                 // track after all cuts
    kNVtxContrib,             // number of primary vertex contibutors
    kNVtxContribTPC,         // number of TPC vertex contibutors
    kNacc,                   // Number of accepted tracks
    kMatchEffITSTPC,         // ruff estimate on the ITS-TPC matching efficiceny
    kNaccTrcklts,            // number of accepted SPD tracklets in |eta|<1.6
    kNaccTrcklts09,          // number of accepted SPD tracklets in |eta|<0.9
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
    kRefMultOvRefMultTPConly,   // ref mult / tpc only ref mult should give a hint on out of bunch pile-up if much higher than factor ~4 (LHC15o)

    kNch,                    // MC true number of charged particles in |eta|<1.6
    kNchJpsiExcl,            // MC true number of charged particles in |eta|<1.6 without J/psi daughter tracks
    kNch05,                  // MC true number of charged particles in |eta|<0.5
    kNch05JpsiExcl,          // MC true number of charged particles in |eta|<0.5 without J/psi daughter tracks
    kNch10,                  // MC true number of charged particles in |eta|<1.0
    kNch10JpsiExcl,          // MC true number of charged particles in |eta|<1.0 without J/psi daughter tracks

    kCentrality,             // event centrality fraction V0M
    kCentralityV0A,          // event centrality fraction V0A
    kCentralityV0C,          // event centrality fraction V0C
    kCentralityZNA,          // event centrality fraction ZNA
    kCentralitySPD,          // centrality using SPD (from second layer)
    //centrality determination for Run 2
    kCentralityNew,          //event centrality V0M
    kCentralityCL0,          //event centrality CL0
    kCentralityCL1,          //event centrality CL1
    kCentralitySPDClusters,          //event centrality SPD
    kCentralitySPDTracklets,          //event centrality SPDTracklets
    kCentralityCL0plus05,    //event centrality VO AP +0.5%
    kCentralityCL0minus05,    //event centrality VO AP -0.5%
    kCentralityCL0plus10,    //event centrality VO AP +1.0%
    kCentralityCL0minus10,    //event centrality VO AP -1.0%



    kTriggerInclONL,         // online trigger bits fired (inclusive)
    kTriggerInclOFF,         // offline trigger bits fired (inclusive)
    kTriggerExclOFF,         // offline only this trigger bit fired (exclusive)
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
  static void SetEventData(const Double_t data[AliDielectronVarManager::kNMaxValues]);
  static Bool_t GetDCA(const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0=0);
  static void SetTPCEventPlane(AliEventplane *const evplane);
  static void SetTPCEventPlaneACremoval(AliDielectronQnEPcorrection *acCuts) {fgQnEPacRemoval = acCuts; fgEventPlaneACremoval = kTRUE;}
  static void SetQnVectorNormalisation(TString qnNorm) {fgQnVectorNorm = qnNorm;}
  static void GetVzeroRP(const AliVEvent* event, Double_t* qvec, Int_t sideOption);      // 0- V0A; 1- V0C; 2- V0A+V0C
  static void GetZDCRP(const AliVEvent* event, Double_t qvec[][2]);
  static AliAODVertex* GetVertex(const AliAODEvent *event, AliAODVertex::AODVtx_t vtype);
  static TProfile* GetEstimatorHistogram(Int_t period, Int_t type) {return fgMultEstimatorAvg[period][type];}
  static Double_t GetTRDpidEfficiency(Int_t runNo, Double_t centrality, Double_t eta, Double_t trdPhi, Double_t pout, Double_t& effErr);
  static Double_t GetSingleLegEff(Double_t * const values);
  static Double_t GetPairEff(Double_t * const values);

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

  static Bool_t Req(ValueTypes var) { return (fgFillMap ? fgFillMap->TestBitNumber(var) : kTRUE); }
  static void FillVarESDtrack(const AliESDtrack *particle,           Double_t * const values);
  static void FillVarAODTrack(const AliAODTrack *particle,           Double_t * const values);
  static void FillVarVTrdTrack(const AliVParticle *particle,         Double_t * const values);
  static void FillVarMCParticle(const AliMCParticle *particle,       Double_t * const values);
  static void FillVarAODMCParticle(const AliAODMCParticle *particle, Double_t * const values);
  static void FillVarDielectronPair(const AliDielectronPair *pair,   Double_t * const values);
  static void FillVarKFParticle(const AliKFParticle *pair,           Double_t * const values);

  static void FillVarVEvent(const AliVEvent *event,                  Double_t * const values);
  static void FillVarESDEvent(const AliESDEvent *event,              Double_t * const values);
  static void FillVarAODEvent(const AliAODEvent *event,              Double_t * const values);
  static void FillVarMCEvent(const AliMCEvent *event,                Double_t * const values);
  static void FillVarTPCEventPlane(const AliEventplane *evplane,     Double_t * const values);
  static void FillQnEventplanes(TList *qnlist,                       Double_t * const values);

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

  static AliDielectronQnEPcorrection *fgQnEPacRemoval; //! filter for auto correlation removal within Qn Framework
  static Bool_t fgEventPlaneACremoval;
  static TString fgQnVectorNorm;                       // String containing the normalisation for the QnVector if the non-default AddTask is used


  static Double_t CalculateEPDiff(Double_t detArp, Double_t detBrp);


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
  ///
  /// Fill track information available in AliVParticle into an array
  /// Also fill event information from local buffer into the array
  ///
  values[AliDielectronVarManager::kPx]        = particle->Px();
  values[AliDielectronVarManager::kPy]        = particle->Py();
  values[AliDielectronVarManager::kPz]        = particle->Pz();
  values[AliDielectronVarManager::kPt]        = particle->Pt();
  values[AliDielectronVarManager::kPtSq]      = particle->Pt()*particle->Pt();
  values[AliDielectronVarManager::kP]         = particle->P();

  values[AliDielectronVarManager::kXv]        = particle->Xv();
  values[AliDielectronVarManager::kYv]        = particle->Yv();
  values[AliDielectronVarManager::kZv]        = particle->Zv();

  values[AliDielectronVarManager::kOneOverPt] = (particle->Pt()>1.0e-3 ? particle->OneOverPt() : 0.0);
  values[AliDielectronVarManager::kPhi]       = TVector2::Phi_0_2pi(particle->Phi());
  values[AliDielectronVarManager::kTheta]     = particle->Theta();
  values[AliDielectronVarManager::kEta]       = particle->Eta();
  values[AliDielectronVarManager::kY]         = particle->Y();

  values[AliDielectronVarManager::kE]         = particle->E();
  values[AliDielectronVarManager::kM]         = particle->M();
  values[AliDielectronVarManager::kCharge]    = particle->Charge();

  values[AliDielectronVarManager::kPdgCode]   = particle->PdgCode();

  values[AliDielectronVarManager::kRndm]      = gRandom->Rndm();

  if(Req(kPtMC)||Req(kPMC)||Req(kPhiMC)||Req(kEtaMC)){
    values[AliDielectronVarManager::kPtMC]      = -999.;
    values[AliDielectronVarManager::kPMC]       = -999.;
    values[AliDielectronVarManager::kPhiMC]     = -999.;
    values[AliDielectronVarManager::kEtaMC]     = -999.;
    AliVParticle *mcTrack(0x0);
    if(AliDielectronMC::Instance()->HasMC())
      mcTrack = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(TMath::Abs(particle->GetLabel()));
    if(mcTrack){
      values[AliDielectronVarManager::kPtMC]   = mcTrack->Pt();
      values[AliDielectronVarManager::kPMC]    = mcTrack->P();
      values[AliDielectronVarManager::kPhiMC]  = TVector2::Phi_0_2pi(mcTrack->Phi());
      values[AliDielectronVarManager::kEtaMC]  = mcTrack->Eta();
    }
  }

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
  values[AliDielectronVarManager::kITSFakeFlag]   = particle->GetITSFakeFlag();
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
  values[AliDielectronVarManager::kTRDchi2Trklt]  = (particle->GetTRDntrackletsPID() > 0 ? particle->GetTRDchi2() / particle->GetTRDntrackletsPID() : -1.);
  values[AliDielectronVarManager::kTRDsignal]     = particle->GetTRDsignal();
  values[AliDielectronVarManager::kTPCclsDiff]    = tpcSignalN-tpcNcls;
  values[AliDielectronVarManager::kTPCclsSegments] = 0.0;

  Double_t itsNclsS = 0.;
  for(int i=0; i<6; i++){
    if( particle->HasSharedPointOnITSLayer(i) ) itsNclsS ++;
  }
  values[AliDielectronVarManager::kNclsITS]     = itsNcls;
  values[AliDielectronVarManager::kNclsSITS]     = itsNclsS;
  values[AliDielectronVarManager::kNclsSFracITS] = itsNcls ? itsNclsS/ itsNcls :0;


  UChar_t threshold = 5;
  TBits tpcClusterMap = particle->GetTPCClusterMap();
  UChar_t n=0; UChar_t j=0;
  for(UChar_t i=0; i<8; ++i) {
    n=0;
    for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliDielectronVarManager::kTPCclsSegments] += 1.0;
  }

  n=0;
  threshold=0;
  values[AliDielectronVarManager::kTPCclsIRO]=0.;
  for(j=0; j<63; ++j) n+=tpcClusterMap.TestBitNumber(j);
  if(n>=threshold) values[AliDielectronVarManager::kTPCclsIRO] = n;
  n=0;
  threshold=0;
  values[AliDielectronVarManager::kTPCclsORO]=0.;
  for(j=63; j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
  if(n>=threshold) values[AliDielectronVarManager::kTPCclsORO] = n;

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
  values[AliDielectronVarManager::kImpactParXYsigma]   = -1.;
  values[AliDielectronVarManager::kImpactParZsigma]    = -1.;

  Float_t dca[2] = {-999.,-999.};
  Float_t dcaRes[3] = {-999.,-999.,-999.};
  esdTrack->GetImpactParameters(dca, dcaRes);
  if(dcaRes[0]>0.) values[AliDielectronVarManager::kImpactParXYsigma] = dca[0]/TMath::Sqrt(dcaRes[0]);
  if(dcaRes[2]>0.) values[AliDielectronVarManager::kImpactParZsigma]  = dca[1]/TMath::Sqrt(dcaRes[2]);


  values[AliDielectronVarManager::kPdgCode]=-1;
  values[AliDielectronVarManager::kPdgCodeMother]=-1;
  values[AliDielectronVarManager::kPdgCodeGrandMother]=-1;
  values[AliDielectronVarManager::kHasCocktailMother]=0;
  values[AliDielectronVarManager::kHasCocktailGrandMother]=0;

  values[AliDielectronVarManager::kNumberOfDaughters]=-999;

  AliDielectronMC *mc=AliDielectronMC::Instance();
  if (mc->HasMC()){
    if (mc->GetMCTrack(particle)) {
      Int_t trkLbl = TMath::Abs(mc->GetMCTrack(particle)->GetLabel());
      values[AliDielectronVarManager::kPdgCode]           =mc->GetMCTrack(particle)->PdgCode();
      values[AliDielectronVarManager::kHasCocktailMother] =mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
      values[AliDielectronVarManager::kPdgCodeMother]     =mc->GetMotherPDG(particle);
      AliMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
      if(motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);
      AliMCParticle *MCpart = mc->GetMCTrack(particle);
      // Fill distance of primary vertex to secondary vertex (as an alternative to the IP)
      // Pure MC variable by intention, no reconstucted value filled.
      if (Req(kDistPrimToSecVtxXYMC) || Req(kDistPrimToSecVtxZMC)) {
        values[AliDielectronVarManager::kDistPrimToSecVtxXYMC] = TMath::Sqrt(  TMath::Power(MCpart->Xv() - values[AliDielectronVarManager::kXvPrimMCtruth],2)
                                                                             + TMath::Power(MCpart->Yv() - values[AliDielectronVarManager::kYvPrimMCtruth],2));
        values[AliDielectronVarManager::kDistPrimToSecVtxZMC] = TMath::Abs(MCpart->Zv() - values[AliDielectronVarManager::kZvPrimMCtruth]);
      }
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
  values[AliDielectronVarManager::kTPCnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kTPCnSigmaEle]   =(fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron) - AliDielectronPID::GetCorrVal() - AliDielectronPID::GetCntrdCorr(particle)) / AliDielectronPID::GetWdthCorr(particle);

  values[AliDielectronVarManager::kTPCnSigmaPio]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
  values[AliDielectronVarManager::kTPCnSigmaMuo]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
  values[AliDielectronVarManager::kTPCnSigmaKao]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
  values[AliDielectronVarManager::kTPCnSigmaPro]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);

  values[AliDielectronVarManager::kITSnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron);
  values[AliDielectronVarManager::kITSnSigmaEle]   =(fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron)
                                                     -AliDielectronPID::GetCntrdCorrITS(particle)
                                                     ) / AliDielectronPID::GetWdthCorrITS(particle);

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

  values[AliDielectronVarManager::kLegEff]        = GetSingleLegEff(values);
  values[AliDielectronVarManager::kOneOverLegEff] = (values[AliDielectronVarManager::kLegEff]>0.0 ? 1./values[AliDielectronVarManager::kLegEff] : 0.0);
  //restore TPC signal if it was changed
  if (esdTrack) esdTrack->SetTPCsignal(origdEdx,esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());

  //fill info from AliVTrdTrack
  if(Req(kTRDonlineA)||Req(kTRDonlineLayerMask)||Req(kTRDonlinePID)||Req(kTRDonlinePt)||Req(kTRDonlineStack)||Req(kTRDonlineTrackInTime)||Req(kTRDonlineSector)||Req(kTRDonlineFlagsTiming)||Req(kTRDonlineLabel)||Req(kTRDonlineNTracklets)||Req(kTRDonlineFirstLayer))
    FillVarVTrdTrack(particle,values);

  if( fgEvent && fgEvent->GetMagneticField() ){
    if(out){
      AliExternalTrackParam out_tmp(*out);
      out_tmp.PropagateTo(AliTRDgeometry::GetXtrdBeg(), fgEvent->GetMagneticField());
      values[AliDielectronVarManager::kTRDeta] = out_tmp.Eta();
    }
    else{
      AliESDtrack particle_tmp(*particle);
      particle_tmp.PropagateTo(AliTRDgeometry::GetXtrdBeg(), fgEvent->GetMagneticField());
      values[AliDielectronVarManager::kTRDeta] = particle_tmp.Eta();
    }
    int mode = particle->GetInnerParam() ? 1:0;
    values[kTPCActiveLength] = particle->GetLengthInActiveZone(mode, 2., 220., fgEvent->GetMagneticField());
    values[kTPCGeomLength] = values[kTPCActiveLength] / ( 130 - TMath::Power( TMath::Abs( particle->GetSigned1Pt() ),1.5 ) );
    values[AliDielectronVarManager::kInTRDacceptance] = TMath::Abs( values[AliDielectronVarManager::kTRDeta] )<0.85 && (  (values[AliDielectronVarManager::kCharge]<0&&(  values[AliDielectronVarManager::kPhi]<1.32 || (values[AliDielectronVarManager::kPhi]>1.98 && values[AliDielectronVarManager::kPhi]<4.10)||  ( values[AliDielectronVarManager::kPhi]>5.12  && values[AliDielectronVarManager::kPhi]<5.48  && TMath::Abs( values[AliDielectronVarManager::kTRDeta] )>0.155 )  || values[AliDielectronVarManager::kPhi]>5.48 )) ||   (values[AliDielectronVarManager::kCharge]>0&&(  values[AliDielectronVarManager::kPhi]<1.52 || (values[AliDielectronVarManager::kPhi]>2.20 && values[AliDielectronVarManager::kPhi]<4.32)||  ( values[AliDielectronVarManager::kPhi]>5.32  && values[AliDielectronVarManager::kPhi]<5.68  && TMath::Abs( values[AliDielectronVarManager::kTRDeta]  )>0.155 )  || values[AliDielectronVarManager::kPhi]>5.68 )) )  ? 1: 0;
  }

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

  // Reset AliESDtrack interface specific information
  if(Req(kNclsITS))      values[AliDielectronVarManager::kNclsITS]       = particle->GetITSNcls();
  if(Req(kITSchi2Cl))    values[AliDielectronVarManager::kITSchi2Cl]     = (particle->GetITSNcls()>0)? particle->GetITSchi2() / particle->GetITSNcls() : 0;
  if(Req(kNclsTPC))      values[AliDielectronVarManager::kNclsTPC]       = tpcNcls;
  if(Req(kNclsSTPC))     values[AliDielectronVarManager::kNclsSTPC]      = tpcNclsS;
  if(Req(kNclsSFracTPC)) values[AliDielectronVarManager::kNclsSFracTPC]  = tpcNcls>0?tpcNclsS/tpcNcls:0;
  if(Req(kNclsTPCiter1)) values[AliDielectronVarManager::kNclsTPCiter1]  = tpcNcls; // not really available in AOD
  if(Req(kNFclsTPC)  || Req(kNFclsTPCfCross))  values[AliDielectronVarManager::kNFclsTPC]      = particle->GetTPCNclsF();
  if(Req(kNFclsTPCr) || Req(kNFclsTPCfCross))  values[AliDielectronVarManager::kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
  if(Req(kNFclsTPCrFrac))  values[AliDielectronVarManager::kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
  if(Req(kNFclsTPCfCross)) values[AliDielectronVarManager::kNFclsTPCfCross]= (values[kNFclsTPC]>0)?(values[kNFclsTPCr]/values[kNFclsTPC]):0;
  if(Req(kNclsTRD))        values[AliDielectronVarManager::kNclsTRD]       = particle->GetNcls(2);
  if(Req(kTRDntracklets))  values[AliDielectronVarManager::kTRDntracklets] = 0;
  if(Req(kTRDpidQuality))  values[AliDielectronVarManager::kTRDpidQuality] = particle->GetTRDntrackletsPID();
  if(Req(kTRDchi2))        values[AliDielectronVarManager::kTRDchi2]       = (particle->GetTRDntrackletsPID()!=0.?particle->GetTRDchi2():-1);
  if(Req(kTRDchi2Trklt))   values[AliDielectronVarManager::kTRDchi2Trklt]  = (particle->GetTRDntrackletsPID()>0 ? particle->GetTRDchi2() / particle->GetTRDntrackletsPID() : -1.);
  if(Req(kTRDsignal))      values[AliDielectronVarManager::kTRDsignal]     = particle->GetTRDsignal();


  TBits tpcClusterMap = particle->GetTPCClusterMap();
  UChar_t n=0; UChar_t j=0;
  UChar_t threshold = 5;

  values[AliDielectronVarManager::kTPCclsSegments] = 0.0;
  if(Req(kTPCclsSegments)) {
    for(UChar_t i=0; i<8; ++i) {
      n=0;
      for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
      if(n>=threshold) values[AliDielectronVarManager::kTPCclsSegments] += 1.0;
    }
  }

  values[AliDielectronVarManager::kTPCclsIRO]=0.;
  if(Req(kTPCclsIRO)) {
    n=0;
    threshold=0;
    for(j=0; j<63; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliDielectronVarManager::kTPCclsIRO] = n;
  }

  values[AliDielectronVarManager::kTPCclsORO]=0.;
  if(Req(kTPCclsORO)) {
    n=0;
    threshold=0;
    for(j=63; j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) values[AliDielectronVarManager::kTPCclsORO] = n;
  }

  // it is stored as normalized to tpcNcls-5 (see AliAnalysisTaskESDfilter)
  if(Req(kTPCchi2Cl))   values[AliDielectronVarManager::kTPCchi2Cl]     = (tpcNcls>0)?particle->Chi2perNDF()*(tpcNcls-5)/tpcNcls:-1.;
  if(Req(kTrackStatus)) values[AliDielectronVarManager::kTrackStatus]   = (Double_t)particle->GetStatus();
  if(Req(kFilterBit))   values[AliDielectronVarManager::kFilterBit]     = (Double_t)particle->GetFilterMap();

  //TRD pidProbs
  values[AliDielectronVarManager::kTRDprobEle]    = 0;
  values[AliDielectronVarManager::kTRDprobPio]    = 0;

  values[AliDielectronVarManager::kTPCsignalN]    = 0;
  values[AliDielectronVarManager::kTPCsignalNfrac]= 0;

  // Fill AliAODTrack interface information
  //
  Int_t v0Index=-1;
  Int_t kinkIndex=-1;
  if( (Req(kV0Index0) || Req(kKinkIndex0)) && particle->GetProdVertex()) {
    v0Index   = particle->GetProdVertex()->GetType()==AliAODVertex::kV0   ? 1 : 0;
    kinkIndex = particle->GetProdVertex()->GetType()==AliAODVertex::kKink ? 1 : 0;
  }
  values[AliDielectronVarManager::kV0Index0]      = v0Index;
  values[AliDielectronVarManager::kKinkIndex0]    = kinkIndex;

  Double_t d0z0[2]={-999.0,-999.0};
  Double_t dcaRes[3] = {-999.,-999.,-999.};
  if(Req(kImpactParXY) || Req(kImpactParZ) || Req(kImpactParXYsigma) || Req(kImpactParZsigma) ) GetDCA(particle, d0z0, dcaRes);
  values[AliDielectronVarManager::kImpactParXY]   = d0z0[0];
  values[AliDielectronVarManager::kImpactParZ]    = d0z0[1];
  values[AliDielectronVarManager::kImpactParXYsigma] = -999.0;
  values[AliDielectronVarManager::kImpactParZsigma] = -999.0;
  if(dcaRes[0]>0.) values[AliDielectronVarManager::kImpactParXYsigma] = d0z0[0]/TMath::Sqrt(dcaRes[0]);
  if(dcaRes[2]>0.) values[AliDielectronVarManager::kImpactParZsigma]  = d0z0[1]/TMath::Sqrt(dcaRes[2]);


  values[AliDielectronVarManager::kPIn]            =  0.;
  values[AliDielectronVarManager::kTPCsignal]      =  0.;
  values[AliDielectronVarManager::kTPCsignalN]     = -1.;
  values[AliDielectronVarManager::kTPCsignalNfrac] = -1.;
  values[AliDielectronVarManager::kTPCclsDiff]     = -999.;

  values[AliDielectronVarManager::kTOFsignal]=0;
  values[AliDielectronVarManager::kTOFbeta]=0;

  values[AliDielectronVarManager::kTPCnSigmaEleRaw]=0;
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

  if(Req(kITSsignal))        values[AliDielectronVarManager::kITSsignal]        =   particle->GetITSsignal();
  if(Req(kITSclusterMap))    values[AliDielectronVarManager::kITSclusterMap]    =   particle->GetITSClusterMap();
  if(Req(kITSLayerFirstCls)) values[AliDielectronVarManager::kITSLayerFirstCls] = -1.;
  for (Int_t iC=0; iC<6; iC++) {
    if (((particle->GetITSClusterMap()) & (1<<(iC))) > 0) {
      if(Req(kITSLayerFirstCls)) values[AliDielectronVarManager::kITSLayerFirstCls] = iC;
      break;
    }
  }

  AliAODPid *pid=const_cast<AliAODPid*>(particle->GetDetPid());
  if (pid) {
    Double_t origdEdx=pid->GetTPCsignal();
    //overwrite signal
    pid->SetTPCsignal(origdEdx/AliDielectronPID::GetEtaCorr(particle)/AliDielectronPID::GetCorrValdEdx());

    Double_t tpcSignalN=0.0;
    if(Req(kTPCsignalN) || Req(kTPCsignalNfrac) || Req(kTPCclsDiff)) tpcSignalN = pid->GetTPCsignalN();
    values[AliDielectronVarManager::kTPCsignalN]     = tpcSignalN;
    values[AliDielectronVarManager::kTPCsignalNfrac] = tpcNcls>0?tpcSignalN/tpcNcls:0;
    values[AliDielectronVarManager::kTPCclsDiff]     = tpcSignalN-tpcNcls;

    values[AliDielectronVarManager::kPIn]         = pid->GetTPCmomentum();
    if(Req(kTPCsignal))   values[AliDielectronVarManager::kTPCsignal]   = pid->GetTPCsignal();
    if(Req(kTOFsignal))   values[AliDielectronVarManager::kTOFsignal]   = pid->GetTOFsignal();
    if(Req(kTOFmismProb)) values[AliDielectronVarManager::kTOFmismProb] = fgPIDResponse->GetTOFMismatchProbability(particle);

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
      values[AliDielectronVarManager::kTOFbeta]  =0;
    }
    else {
      t *= 1e-12; //ps -> s

	Double_t v = l / t;
	Float_t beta = v / TMath::C();
	values[AliDielectronVarManager::kTOFbeta]=beta;
      }
    }

    // nsigma for various detectors
    if(Req(kTPCnSigmaEleRaw)) values[kTPCnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron);
    if(Req(kTPCnSigmaEle))    values[kTPCnSigmaEle]   =(fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kElectron)-AliDielectronPID::GetCorrVal()-AliDielectronPID::GetCntrdCorr(particle)) / AliDielectronPID::GetWdthCorr(particle);

    if(Req(kTPCnSigmaPio)) values[kTPCnSigmaPio]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kPion);
    if(Req(kTPCnSigmaMuo)) values[kTPCnSigmaMuo]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kMuon);
    if(Req(kTPCnSigmaKao)) values[kTPCnSigmaKao]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kKaon);
    if(Req(kTPCnSigmaPro)) values[kTPCnSigmaPro]=fgPIDResponse->NumberOfSigmasTPC(particle,AliPID::kProton);

    if(Req(kITSnSigmaEleRaw)) values[kITSnSigmaEleRaw]= fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron);
    if(Req(kITSnSigmaEle))    values[kITSnSigmaEle]   =(fgPIDResponse->NumberOfSigmasITS(particle,AliPID::kElectron) - AliDielectronPID::GetCntrdCorrITS(particle)) / AliDielectronPID::GetWdthCorrITS(particle);

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
    // 1D TRD PID
    if( Req(kTRDprobEle) || Req(kTRDprobPio) ){
      fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob);
      values[AliDielectronVarManager::kTRDprobEle]      = prob[AliPID::kElectron];
      values[AliDielectronVarManager::kTRDprobPio]      = prob[AliPID::kPion];
    }
    // 2D TRD PID
    if( Req(kTRDprob2DEle) || Req(kTRDprob2DPio) || Req(kTRDprob2DPro) ){
      fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob, AliTRDPIDResponse::kLQ2D);
      values[AliDielectronVarManager::kTRDprob2DEle]    = prob[AliPID::kElectron];
      values[AliDielectronVarManager::kTRDprob2DPio]    = prob[AliPID::kPion];
      values[AliDielectronVarManager::kTRDprob2DPro]    = prob[AliPID::kProton];
    }
    // 3D TRD PID
     if( Req(kTRDprob3DEle) || Req(kTRDprob3DPio) || Req(kTRDprob3DPro) ){
       fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob, AliTRDPIDResponse::kLQ3D);
       values[AliDielectronVarManager::kTRDprob3DEle]    = prob[AliPID::kElectron];
       values[AliDielectronVarManager::kTRDprob3DPio]    = prob[AliPID::kPion];
       values[AliDielectronVarManager::kTRDprob3DPro]    = prob[AliPID::kProton];
     }
    // 7D TRD PID
     if( Req(kTRDprob7DEle) || Req(kTRDprob7DPio) || Req(kTRDprob7DPro) ){
       fgPIDResponse->ComputeTRDProbability(particle,AliPID::kSPECIES,prob, AliTRDPIDResponse::kLQ7D);
       values[AliDielectronVarManager::kTRDprob7DEle]    = prob[AliPID::kElectron];
       values[AliDielectronVarManager::kTRDprob7DPio]    = prob[AliPID::kPion];
       values[AliDielectronVarManager::kTRDprob7DPro]    = prob[AliPID::kProton];
     }


    //restore TPC signal if it was changed
    pid->SetTPCsignal(origdEdx);
  }

  //EMCAL PID information
  Double_t eop=0;
  Double_t showershape[4]={0.,0.,0.,0.};
//   if(Req()) values[AliDielectronVarManager::kEMCALnSigmaEle]  = fgPIDResponse->NumberOfSigmasEMCAL(particle,AliPID::kElectron);
  if(Req(kEMCALnSigmaEle) || Req(kEMCALE) || Req(kEMCALEoverP) ||
     Req(kEMCALNCells) || Req(kEMCALM02) || Req(kEMCALM20) || Req(kEMCALDispersion))
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
      Int_t trkLbl = TMath::Abs(mc->GetMCTrack(particle)->GetLabel());
      values[AliDielectronVarManager::kPdgCode]           =mc->GetMCTrack(particle)->PdgCode();
      values[AliDielectronVarManager::kHasCocktailMother] =mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
      values[AliDielectronVarManager::kPdgCodeMother]     =mc->GetMotherPDG(particle);
      AliAODMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
      if(motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);
    }
    values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);
  } //if(mc->HasMC())

  if(Req(kTOFPIDBit))     values[AliDielectronVarManager::kTOFPIDBit]=(particle->GetStatus()&AliESDtrack::kTOFpid? 1: 0);
  values[AliDielectronVarManager::kLegEff]=0.0;
  values[AliDielectronVarManager::kOneOverLegEff]=0.0;
  if(Req(kLegEff) || Req(kOneOverLegEff)) {
    values[AliDielectronVarManager::kLegEff] = GetSingleLegEff(values);
    values[AliDielectronVarManager::kOneOverLegEff] = (values[AliDielectronVarManager::kLegEff]>0.0 ? 1./values[AliDielectronVarManager::kLegEff] : 0.0);
  }

  //fill info from AliVTrdTrack
  if(Req(kTRDonlineA)||Req(kTRDonlineLayerMask)||Req(kTRDonlinePID)||Req(kTRDonlinePt)||Req(kTRDonlineStack)||Req(kTRDonlineSector)||Req(kTRDonlineTrackInTime)||Req(kTRDonlineFlagsTiming)||Req(kTRDonlineLabel)||Req(kTRDonlineNTracklets)||Req(kTRDonlineFirstLayer))
    FillVarVTrdTrack(particle,values);
}

inline void AliDielectronVarManager::FillVarVTrdTrack(const AliVParticle *particle, Double_t * const values)
{


  //Initialisation of values
  values[AliDielectronVarManager::kTRDonlineLayerMask] = -1.0;
  values[AliDielectronVarManager::kTRDonlinePID] = -1.0 ;
  values[AliDielectronVarManager::kTRDonlinePt] = 0;
  values[AliDielectronVarManager::kTRDonlineStack] = -1.0;
  values[AliDielectronVarManager::kTRDonlineSector] = -1.0;
  values[AliDielectronVarManager::kTRDonlineTrackInTime] = -1.0;
  values[AliDielectronVarManager::kTRDonlineFlagsTiming] = -1.0;
  //	if(Req(kTRDonlineLabel))values[AliDielectronVarManager::kTRDonlineLabel] = ; ???
  values[AliDielectronVarManager::kTRDonlineNTracklets]= -1.0;
  values[AliDielectronVarManager::kTRDonlineFirstLayer] = -1.;


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

	  values[AliDielectronVarManager::kTRDonlineA] = gtutrk->GetA();

      values[AliDielectronVarManager::kTRDonlineLayerMask] = gtutrk->GetLayerMask();
      values[AliDielectronVarManager::kTRDonlinePID] = gtutrk->GetPID() ;
      values[AliDielectronVarManager::kTRDonlinePt] = gtutrk->Pt();
      values[AliDielectronVarManager::kTRDonlineStack] = gtutrk->GetStack();
      values[AliDielectronVarManager::kTRDonlineSector] = gtutrk->GetSector();
      values[AliDielectronVarManager::kTRDonlineTrackInTime] = gtutrk->GetTrackInTime();
      values[AliDielectronVarManager::kTRDonlineFlagsTiming] = gtutrk->GetFlagsTiming();
      values[AliDielectronVarManager::kTRDonlineLabel] = gtutrk->GetLabel();
      values[AliDielectronVarManager::kTRDonlineNTracklets]= gtutrk->GetNTracklets();

      for (Int_t iC=0; iC<6; iC++) {
	    if (((gtutrk->GetLayerMask()) & (1<<(iC))) > 0) {
	      values[AliDielectronVarManager::kTRDonlineFirstLayer] = iC;
	      break;
	    }
	  }
      }//if matching
      //TO DO: what is the initialising value, if no match? -1? is this a problem?, is the PT-signed?
  }//for loop over gtutracks

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
  values[AliDielectronVarManager::kDistPrimToSecVtxXYMC] = 0;
  values[AliDielectronVarManager::kDistPrimToSecVtxZMC] = 0;
  values[AliDielectronVarManager::kPIn]           = 0;
  values[AliDielectronVarManager::kYsignedIn]     = 0;
  values[AliDielectronVarManager::kTPCsignal]     = 0;
  values[AliDielectronVarManager::kTOFsignal]     = 0;
  values[AliDielectronVarManager::kTOFbeta]       = 0;
  values[AliDielectronVarManager::kTPCnSigmaEleRaw]  = 0;
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

  // Fill distance of primary vertex to secondary vertex (as a well-defined alternative to the IP-approximation below)
  if (Req(kDistPrimToSecVtxXYMC) || Req(kDistPrimToSecVtxZMC)) {
    values[AliDielectronVarManager::kDistPrimToSecVtxXYMC] = TMath::Sqrt(  TMath::Power(particle->Xv() - values[AliDielectronVarManager::kXvPrim],2)
                                                                         + TMath::Power(particle->Yv() - values[AliDielectronVarManager::kYvPrim],2));
    values[AliDielectronVarManager::kDistPrimToSecVtxZMC] = TMath::Abs(particle->Zv() - values[AliDielectronVarManager::kZvPrim]);
  }
  //Approximation of the Impact Parameter
  //Get TVectors for primary and secondary vertex as well as particle momentum
  // distance of space point to a straight line
  // d = |b x (p-a)|/|b|
  TVector3 priVtx(values[AliDielectronVarManager::kXvPrim],values[AliDielectronVarManager::kYvPrim],values[AliDielectronVarManager::kZvPrim]);
  TVector3 secVtx(values[AliDielectronVarManager::kXv],values[AliDielectronVarManager::kYv],values[AliDielectronVarManager::kZv]);
  TVector3 momPart(values[AliDielectronVarManager::kPx],values[AliDielectronVarManager::kPy],values[AliDielectronVarManager::kPz]);
  priVtx -= secVtx;
  TVector3 denom = momPart.Cross(priVtx);
  values[AliDielectronVarManager::kImpactParXY]   = TMath::Sqrt(denom.X()*denom.X() + denom.Y()*denom.Y()) / TMath::Sqrt(momPart.X()*momPart.X() + momPart.Y()*momPart.Y() + momPart.Z()*momPart.Z());
  values[AliDielectronVarManager::kImpactParZ]   = TMath::Abs(denom.Z()) / TMath::Sqrt(momPart.X()*momPart.X() + momPart.Y()*momPart.Y() + momPart.Z()*momPart.Z());


  // Fill AliMCParticle interface specific information
  AliDielectronMC *mc=AliDielectronMC::Instance();
  Int_t trkLbl = TMath::Abs(particle->GetLabel());
  values[AliDielectronVarManager::kPdgCode]           = particle->PdgCode();
  values[AliDielectronVarManager::kHasCocktailMother] = mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
  values[AliDielectronVarManager::kPdgCodeMother]     = mc->GetMotherPDG(particle);
  AliMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
  if(motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);


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
  values[AliDielectronVarManager::kTPCnSigmaEleRaw]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaEle]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPio]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaMuo]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaKao]  = 0;
  values[AliDielectronVarManager::kTPCnSigmaPro]  = 0;
  values[AliDielectronVarManager::kITSclusterMap] = 0;

  values[AliDielectronVarManager::kPdgCode]       = -1;
  values[AliDielectronVarManager::kPdgCodeMother] = -1;
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
  values[AliDielectronVarManager::kPtSq]      = values[AliDielectronVarManager::kPt] * values[AliDielectronVarManager::kPt];
  values[AliDielectronVarManager::kP]         = TMath::Sqrt(values[AliDielectronVarManager::kPt]*
						      values[AliDielectronVarManager::kPt]+
						      values[AliDielectronVarManager::kPz]*
						      values[AliDielectronVarManager::kPz]);

  values[AliDielectronVarManager::kXv]        = 0;
  values[AliDielectronVarManager::kYv]        = 0;
  values[AliDielectronVarManager::kZv]        = 0;

  values[AliDielectronVarManager::kOneOverPt] = (values[AliDielectronVarManager::kPt]>1.0e-6 ? 1.0/values[AliDielectronVarManager::kPt] : 0.0);
  values[AliDielectronVarManager::kPhi]       = TVector2::Phi_0_2pi( TMath::ATan2(values[AliDielectronVarManager::kPy],values[AliDielectronVarManager::kPx]) );
  values[AliDielectronVarManager::kTheta]     = TMath::ATan2(values[AliDielectronVarManager::kPt],values[AliDielectronVarManager::kPz]);
  values[AliDielectronVarManager::kEta]       = ((values[AliDielectronVarManager::kP]-values[AliDielectronVarManager::kPz])>1.0e-6 && (values[AliDielectronVarManager::kP]+values[AliDielectronVarManager::kPz])>1.0e-6 ? 0.5*TMath::Log((values[AliDielectronVarManager::kP]+values[AliDielectronVarManager::kPz])/(values[AliDielectronVarManager::kP]-values[AliDielectronVarManager::kPz])) : -9999.);
  values[AliDielectronVarManager::kE]         = p1->E()+p2->E();
  values[AliDielectronVarManager::kY]         = ((values[AliDielectronVarManager::kE]-values[AliDielectronVarManager::kPz])>1.0e-6 && (values[AliDielectronVarManager::kE]+values[AliDielectronVarManager::kPz])>1.0e-6 ? 0.5*TMath::Log((values[AliDielectronVarManager::kE]+values[AliDielectronVarManager::kPz])/(values[AliDielectronVarManager::kE]-values[AliDielectronVarManager::kPz])) : -9999.);
  values[AliDielectronVarManager::kCharge]    = p1->Charge()+p2->Charge();

  values[AliDielectronVarManager::kM]         = p1->M()*p1->M()+p2->M()*p2->M()+
                       2.0*(p1->E()*p2->E()-p1->Px()*p2->Px()-p1->Py()*p2->Py()-p1->Pz()*p2->Pz());
  values[AliDielectronVarManager::kM]         = (values[AliDielectronVarManager::kM]>1.0e-8 ? TMath::Sqrt(values[AliDielectronVarManager::kM]) : -1.0);
  values[AliDielectronVarManager::kMMC] = values[AliDielectronVarManager::kM];
  values[AliDielectronVarManager::kPtMC] = values[AliDielectronVarManager::kPt];

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
  values[AliDielectronVarManager::kTPCnSigmaEleRaw]  = 0;
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

  // Fill AliAODMCParticle interface specific information
  AliDielectronMC *mc=AliDielectronMC::Instance();
  Int_t trkLbl = TMath::Abs(particle->GetLabel());
  values[AliDielectronVarManager::kPdgCode]           = particle->PdgCode();
  values[AliDielectronVarManager::kHasCocktailMother] = mc->CheckParticleSource(trkLbl, AliDielectronSignalMC::kDirect);
  values[AliDielectronVarManager::kPdgCodeMother]     = mc->GetMotherPDG(particle);
  AliAODMCParticle *motherMC=mc->GetMCTrackMother(particle); //mother
  if(motherMC) values[AliDielectronVarManager::kPdgCodeGrandMother]=mc->GetMotherPDG(motherMC);

  values[AliDielectronVarManager::kIsJpsiPrimary] = mc->IsJpsiPrimary(particle);
  values[AliDielectronVarManager::kNumberOfDaughters]=mc->NumberOfDaughters(particle);

  // using AODMCHEader information
  AliAODMCHeader *mcHeader = (AliAODMCHeader*)fgEvent->FindListObject(AliAODMCHeader::StdBranchName());
  if(mcHeader) {
    values[AliDielectronVarManager::kImpactParZ]  = mcHeader->GetVtxZ()-particle->Zv();
    values[AliDielectronVarManager::kImpactParXY] = TMath::Sqrt(TMath::Power(mcHeader->GetVtxX()-particle->Xv(),2) +
								TMath::Power(mcHeader->GetVtxY()-particle->Yv(),2));
  }

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
  FillVarVParticle(pair, values); // this also filles the event information into 'values'.

  // Fill AliDielectronPair specific information
  const AliKFParticle &kfPair = pair->GetKFParticle();


  values[AliDielectronVarManager::kThetaHE]      = 0.0;
  values[AliDielectronVarManager::kPhiHE]        = 0.0;
  values[AliDielectronVarManager::kThetaSqHE]    = 0.0;
  values[AliDielectronVarManager::kCos2PhiHE]    = 0.0;
  values[AliDielectronVarManager::kCosTilPhiHE]  = 0.0;

  values[AliDielectronVarManager::kThetaCS]      = 0.0;
  values[AliDielectronVarManager::kPhiCS]        = 0.0;
  values[AliDielectronVarManager::kThetaSqCS]    = 0.0;
  values[AliDielectronVarManager::kCos2PhiCS]    = 0.0;
  values[AliDielectronVarManager::kCosTilPhiCS]  = 0.0;

  Double_t thetaHE=0;
  Double_t phiHE=0;
  Double_t thetaCS=0;
  Double_t phiCS=0;
  if(Req(kThetaHE) || Req(kPhiHE) || Req(kThetaCS) || Req(kPhiCS)) {
    pair->GetThetaPhiCM(thetaHE,phiHE,thetaCS,phiCS);

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
  }

  if(Req(kChi2NDF))          values[AliDielectronVarManager::kChi2NDF]          = kfPair.GetChi2()/kfPair.GetNDF();
  if(Req(kDecayLength))      values[AliDielectronVarManager::kDecayLength]      = kfPair.GetDecayLength();
  if(Req(kR))                values[AliDielectronVarManager::kR]                = kfPair.GetR();
  if(Req(kOpeningAngle))     values[AliDielectronVarManager::kOpeningAngle]     = pair->OpeningAngle();
  if(Req(kOpeningAngleXY))     values[AliDielectronVarManager::kOpeningAngleXY] = pair->OpeningAngleXY();
  if(Req(kOpeningAngleRZ))     values[AliDielectronVarManager::kOpeningAngleRZ] = pair->OpeningAngleRZ();
  if(Req(kCosPointingAngle)) values[AliDielectronVarManager::kCosPointingAngle] = fgEvent ? pair->GetCosPointingAngle(fgEvent->GetPrimaryVertex()) : -1;

  if(Req(kLegDist))   values[AliDielectronVarManager::kLegDist]      = pair->DistanceDaughters();
  if(Req(kLegDistXY)) values[AliDielectronVarManager::kLegDistXY]    = pair->DistanceDaughtersXY();
  if(Req(kDeltaEta))  values[AliDielectronVarManager::kDeltaEta]     = pair->DeltaEta();
  if(Req(kDeltaPhi))  values[AliDielectronVarManager::kDeltaPhi]     = pair->DeltaPhi();
  if(Req(kMerr))      values[AliDielectronVarManager::kMerr]         = kfPair.GetErrMass()>1e-30&&kfPair.GetMass()>1e-30?kfPair.GetErrMass()/kfPair.GetMass():1000000;

  values[AliDielectronVarManager::kPairType]     = pair->GetType();
  // Armenteros-Podolanski quantities
  if(Req(kArmAlpha)) values[AliDielectronVarManager::kArmAlpha]     = pair->GetArmAlpha();
  if(Req(kArmPt))    values[AliDielectronVarManager::kArmPt]        = pair->GetArmPt();

  if(Req(kPsiPair))  values[AliDielectronVarManager::kPsiPair]      = fgEvent ? pair->PsiPair(fgEvent->GetMagneticField()) : -5;
  if(Req(kPhivPair)) values[AliDielectronVarManager::kPhivPair]      = fgEvent ? pair->PhivPair(fgEvent->GetMagneticField()) : -5;
  if(Req(kDeltaCotTheta)) values[kDeltaCotTheta] =  pair->DeltaCotTheta();
  if(Req(kTriangularConversionCut)) values[AliDielectronVarManager::kTriangularConversionCut] = fgEvent ? pair->PhivPair(fgEvent->GetMagneticField()) - 21. * pair->M() : -999.;
  if(Req(kPseudoProperTime) || Req(kPseudoProperTimeErr)) {
    values[AliDielectronVarManager::kPseudoProperTime] =
      fgEvent ? kfPair.GetPseudoProperDecayTime(*(fgEvent->GetPrimaryVertex()), TDatabasePDG::Instance()->GetParticle(443)->Mass(), &errPseudoProperTime2 ) : -1e10;
  // values[AliDielectronVarManager::kPseudoProperTime] = fgEvent ? pair->GetPseudoProperTime(fgEvent->GetPrimaryVertex()): -1e10;
    values[AliDielectronVarManager::kPseudoProperTimeErr] = (errPseudoProperTime2 > 0) ? TMath::Sqrt(errPseudoProperTime2) : -1e10;
  }

  // impact parameter
  Double_t d0z0[2]={-999., -999.};
  if( (Req(kImpactParXY) || Req(kImpactParZ)) && fgEvent) pair->GetDCA(fgEvent->GetPrimaryVertex(), d0z0);
  values[AliDielectronVarManager::kImpactParXY]   = d0z0[0];
  values[AliDielectronVarManager::kImpactParZ]    = d0z0[1];


  //calculate pair dca in sigma and cm
  values[AliDielectronVarManager::kPairDCAsigXY]     = -999.;
  values[AliDielectronVarManager::kPairDCAsigZ]      = -999.;
  values[AliDielectronVarManager::kPairDCAabsXY]     = -999.;
  values[AliDielectronVarManager::kPairDCAabsZ]      = -999.;
  values[AliDielectronVarManager::kPairLinDCAsigXY]  = -999.;
  values[AliDielectronVarManager::kPairLinDCAsigZ]   = -999.;
  values[AliDielectronVarManager::kPairLinDCAabsXY]  = -999.;
  values[AliDielectronVarManager::kPairLinDCAabsZ]   = -999.;
  values[AliDielectronVarManager::kLeg1DCAsigXY]     = -999.;
  values[AliDielectronVarManager::kLeg1DCAabsXY]     = -999.;
  values[AliDielectronVarManager::kLeg1DCAresXY]     = -999.;

  // check if calculation is requested
  if(Req(kPairDCAsigXY) || Req(kPairDCAsigZ) || Req(kPairDCAabsXY) || Req(kPairDCAabsZ) ||
     Req(kPairLinDCAsigXY) || Req(kPairLinDCAsigZ) || Req(kPairLinDCAabsXY) || Req(kPairLinDCAabsZ)) {

    // get track references from pair
    AliVParticle* d1 = pair-> GetFirstDaughterP();
    AliVParticle* d2 = pair->GetSecondDaughterP();

    if (d1 && d2) {
      // check for ESD or AOD
      Bool_t isESD = (d1->IsA() == AliESDtrack::Class());

      if (d1->IsA() == d2->IsA()) { // Don't mix AOD with ESD. Needed because AliAnalysisTaskRandomRejection always creates AliAODTracks (should be fixed).

        ////// first daughter
        Double_t dca1[2]       = {-999.,-999.};      // xy,z absolute values
        Double_t dcaSig1[2]    = {-999.,-999.};      // xy,z sigma values
        Double_t dcaRes1[3]    = {-999.,-999.,-999.};// Covariance matrix
        //Float_t dcaTPC1[2]    = {-999.,-999.};      // xy,z TPC-only absolute values
        //Float_t dcaSigTPC1[2] = {-999.,-999.};      // xy,z TPC-only sigma values
        //Float_t dcaResTPC1[3] = {-999.,-999.,-999.};// Covariance matrix TPC-only

        ////// second daughter
        Double_t dca2[2]       = {-999.,-999.};      // xy,z absolute values
        Double_t dcaSig2[2]    = {-999.,-999.};      // xy,z sigma values
        Double_t dcaRes2[3]    = {-999.,-999.,-999.};// Covariance matrix
        //Float_t dcaTPC2[2]    = {-999.,-999.};      // xy,z TPC-only absolute values
        //Float_t dcaSigTPC2[2] = {-999.,-999.};      // xy,z TPC-only sigma values
        //Float_t dcaResTPC2[3] = {-999.,-999.,-999.};// Covariance matrix TPC-only

        if (isESD) {
          // 'Float_t' needed for 'virtual void AliESDtrack::GetImpactParameters(Float_t p[2], Float_t cov[3]) const'
          Float_t dca_tmp[2] = {-999.,-999.};
          Float_t res_tmp[3] = {-999.,-999.,-999.};
          static_cast<AliESDtrack*>(d1)->GetImpactParameters(dca_tmp, res_tmp);
          dca1[0]   =dca_tmp[0]; dca1[1]   =dca_tmp[1];
          dcaRes1[0]=res_tmp[0]; dcaRes1[1]=res_tmp[1]; dcaRes1[2]=res_tmp[2];
          //static_cast<AliESDtrack*>(d1)->GetImpactParametersTPC(dcaTPC1, dcaResTPC1);

          dca_tmp[0]=-999.; dca_tmp[1]=-999.;
          res_tmp[0]=-999.; res_tmp[1]=-999.; res_tmp[2]=-999.;
          static_cast<AliESDtrack*>(d2)->GetImpactParameters(dca_tmp, res_tmp);
          dca2[0]   =dca_tmp[0]; dca2[1]   =dca_tmp[1];
          dcaRes2[0]=res_tmp[0]; dcaRes2[1]=res_tmp[1]; dcaRes2[2]=res_tmp[2];
          //static_cast<AliESDtrack*>(d2)->GetImpactParametersTPC(dcaTPC2, dcaResTPC2);
        }
        else { // AOD
          GetDCA(static_cast<AliAODTrack*>(d1), dca1, dcaRes1);
          GetDCA(static_cast<AliAODTrack*>(d2), dca2, dcaRes2);
        }

        // compute normalized DCAs
        // neglect the mixed term 'dcaResX[1]'
        if(dcaRes1[0]>0.) dcaSig1[0] = dca1[0]/TMath::Sqrt(dcaRes1[0]);
        if(dcaRes1[2]>0.) dcaSig1[1] = dca1[1]/TMath::Sqrt(dcaRes1[2]);
        //if(dcaResTPC1[0]>0.) dcaSigTPC1[0] = dcaTPC1[0]/TMath::Sqrt(dcaResTPC1[0]);
        //if(dcaResTPC1[2]>0.) dcaSigTPC1[1] = dcaTPC1[1]/TMath::Sqrt(dcaResTPC1[2]);
        if(dcaRes2[0]>0.) dcaSig2[0] = dca2[0]/TMath::Sqrt(dcaRes2[0]);
        if(dcaRes2[2]>0.) dcaSig2[1] = dca2[1]/TMath::Sqrt(dcaRes2[2]);
        //if(dcaResTPC2[0]>0.) dcaSigTPC2[0] = dcaTPC2[0]/TMath::Sqrt(dcaResTPC2[0]);
        //if(dcaResTPC2[2]>0.) dcaSigTPC2[1] = dcaTPC2[1]/TMath::Sqrt(dcaResTPC2[2]);

        // set first daughter variables for cross-checks
        values[AliDielectronVarManager::kLeg1DCAabsXY]   = dca1[0];
        values[AliDielectronVarManager::kLeg1DCAsigXY]   = dcaSig1[0];
        values[AliDielectronVarManager::kLeg1DCAresXY]   = dcaRes1[0];

        // set pair dca values
        // quadratic summation
        values[AliDielectronVarManager::kPairDCAabsXY]       = TMath::Sqrt( (dca1[0]*dca1[0] + dca2[0]*dca2[0]) / 2 );
        values[AliDielectronVarManager::kPairDCAabsZ]        = TMath::Sqrt( (dca1[1]*dca1[1] + dca2[1]*dca2[1]) / 2 );
        values[AliDielectronVarManager::kPairDCAsigXY]       = TMath::Sqrt( (dcaSig1[0]*dcaSig1[0] + dcaSig2[0]*dcaSig2[0]) / 2 );
        values[AliDielectronVarManager::kPairDCAsigZ]        = TMath::Sqrt( (dcaSig1[1]*dcaSig1[1] + dcaSig2[1]*dcaSig2[1]) / 2 );
        // linear summation
        values[AliDielectronVarManager::kPairLinDCAabsXY]    = (TMath::Abs(dca1[0]) + TMath::Abs(dca2[0])) / 2;
        values[AliDielectronVarManager::kPairLinDCAabsZ]     = (TMath::Abs(dca1[1]) + TMath::Abs(dca2[1])) / 2;
        values[AliDielectronVarManager::kPairLinDCAsigXY]    = (TMath::Abs(dcaSig1[0]) + TMath::Abs(dcaSig2[0])) / 2;
        values[AliDielectronVarManager::kPairLinDCAsigZ]     = (TMath::Abs(dcaSig1[1]) + TMath::Abs(dcaSig2[1])) / 2;
      }
    }
  }


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
	values[AliDielectronVarManager::kPtSq]      = values[AliDielectronVarManager::kPt] * values[AliDielectronVarManager::kPt];

	values[AliDielectronVarManager::kP]         =  (lv1+lv2).P();

	//Not overwritten, could take event vertex in next iteration
	values[AliDielectronVarManager::kXv]        = (lv1+lv2).X();
	values[AliDielectronVarManager::kYv]        = (lv1+lv2).Y();
	values[AliDielectronVarManager::kZv]        = (lv1+lv2).Z();

	values[AliDielectronVarManager::kE]         = (lv1+lv2).E();

	values[AliDielectronVarManager::kM]         = (lv1+lv2).M();

	values[AliDielectronVarManager::kOpeningAngle] =  lv1.Angle(lv2.Vect());

	values[AliDielectronVarManager::kOneOverPt] = (values[AliDielectronVarManager::kPt]>0. ? 1./values[AliDielectronVarManager::kPt] : -9999.);
	values[AliDielectronVarManager::kPhi]       = TVector2::Phi_0_2pi( (lv1+lv2).Phi() );
	values[AliDielectronVarManager::kEta]       = (lv1+lv2).Eta();

	values[AliDielectronVarManager::kY]       = (lv1+lv2).Rapidity();

	// Fill AliDielectronPair specific information
	values[AliDielectronVarManager::kDeltaEta]     = TMath::Abs(feta1 -feta2 );
	values[AliDielectronVarManager::kDeltaPhi]     = lv1.DeltaPhi(lv2);

       if( Req(kDeltaPhiChargeOrdered) && fgEvent ) values[AliDielectronVarManager::kDeltaPhiChargeOrdered] = fD1.GetQ() * fgEvent->GetMagneticField() > 0 ? lv1.Phi() - lv2.Phi() :lv2.Phi() - lv1.Phi() ;
	values[AliDielectronVarManager::kPairType]     = pair->GetType();

        // Calculate pair variables for corresponding generated pair
        if(AliDielectronMC::Instance()->HasMC() && (Req(kMMC)||Req(kPtMC)||Req(kPMC)||Req(kEtaMC)||Req(kPhiMC))){
          values[AliDielectronVarManager::kMMC]   = -999.;
          values[AliDielectronVarManager::kPtMC]  = -999.;
          values[AliDielectronVarManager::kPMC]   = -999.;
          values[AliDielectronVarManager::kPhiMC] = -999.;
          values[AliDielectronVarManager::kEtaMC] = -999.;
          AliVParticle *mcDaughter1 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(TMath::Abs((pair->GetFirstDaughterP() )->GetLabel()));
          AliVParticle *mcDaughter2 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(TMath::Abs((pair->GetSecondDaughterP())->GetLabel()));
          if(mcDaughter1 && mcDaughter2){
            TLorentzVector lv1MC,lv2MC;
            lv1MC.SetPtEtaPhiM(mcDaughter1->Pt(),mcDaughter1->Eta(),mcDaughter1->Phi(),mElectron);
            lv2MC.SetPtEtaPhiM(mcDaughter2->Pt(),mcDaughter2->Eta(),mcDaughter2->Phi(),mElectron);
            values[AliDielectronVarManager::kMMC]   = (lv1MC+lv2MC).M();
            values[AliDielectronVarManager::kPtMC]  = (lv1MC+lv2MC).Pt();
            values[AliDielectronVarManager::kPMC]   = (lv1MC+lv2MC).P();
            values[AliDielectronVarManager::kPhiMC] = TVector2::Phi_0_2pi( (lv1MC+lv2MC).Phi() );
            values[AliDielectronVarManager::kEtaMC] = (lv1MC+lv2MC).Eta();
          }
        }
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

    if(Req(kOpeningAngleCorr)) {
      Float_t a = 1.54e-01;
      values[AliDielectronVarManager::kOpeningAngleCorr]  =
        values[AliDielectronVarManager::kOpeningAngle]
        - a * TMath::Sqrt(  values[AliDielectronVarManager::kPairDCAabsXY] * values[AliDielectronVarManager::kOneOverPt] );
    }

    if(Req(kMCorr)) {
      Float_t a =  7.59e-02;
      values[AliDielectronVarManager::kMCorr]  =
        values[AliDielectronVarManager::kM]
        - a * TMath::Sqrt( values[AliDielectronVarManager::kPairDCAabsXY] * values[AliDielectronVarManager::kPt] );
    }

  }
  //common, regardless of calculation method

  // Flow quantities
  Double_t phi=values[AliDielectronVarManager::kPhi];
  if(Req(kCosPhiH2)) values[AliDielectronVarManager::kCosPhiH2] = TMath::Cos(2*phi);
  if(Req(kSinPhiH2)) values[AliDielectronVarManager::kSinPhiH2] = TMath::Sin(2*phi);
  Double_t delta=0.0;
  // v2 with respect to VZERO-A event plane
  delta = TVector2::Phi_mpi_pi(phi - fgData[AliDielectronVarManager::kV0ArpH2]);
  if(Req(kV0ArpH2FlowV2))   values[AliDielectronVarManager::kV0ArpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
  if(Req(kDeltaPhiV0ArpH2)) values[AliDielectronVarManager::kDeltaPhiV0ArpH2] = delta;
  // v2 with respect to VZERO-C event plane
  delta = TVector2::Phi_mpi_pi(phi - fgData[AliDielectronVarManager::kV0CrpH2]);
  if(Req(kV0CrpH2FlowV2))   values[AliDielectronVarManager::kV0CrpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
  if(Req(kDeltaPhiV0CrpH2)) values[AliDielectronVarManager::kDeltaPhiV0CrpH2] = delta;
  // v2 with respect to the combined VZERO-A and VZERO-C event plane
  delta = TVector2::Phi_mpi_pi(phi - fgData[AliDielectronVarManager::kV0ACrpH2]);
  if(Req(kV0ACrpH2FlowV2))   values[AliDielectronVarManager::kV0ACrpH2FlowV2] = TMath::Cos(2.0*delta);  // 2nd harmonic flow coefficient
  if(Req(kDeltaPhiV0ACrpH2)) values[AliDielectronVarManager::kDeltaPhiV0ACrpH2] = delta;


  // quantities using the values of  AliEPSelectionTask , interval [-pi,+pi]
  values[AliDielectronVarManager::kDeltaPhiv0ArpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kv0ArpH2]);
  values[AliDielectronVarManager::kDeltaPhiv0CrpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kv0CrpH2]);
  values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kv0ACrpH2]);
  values[AliDielectronVarManager::kDeltaPhiTPCrpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kTPCrpH2]);
  values[AliDielectronVarManager::kv0ACrpH2FlowV2]   = TMath::Cos( 2.*values[AliDielectronVarManager::kDeltaPhiv0ACrpH2] );
  values[AliDielectronVarManager::kv0ArpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kDeltaPhiv0ArpH2] );
  values[AliDielectronVarManager::kv0CrpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kDeltaPhiv0CrpH2] );
  values[AliDielectronVarManager::kTPCrpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kDeltaPhiTPCrpH2] );
  values[AliDielectronVarManager::kTPCrpH2FlowV2Sin] = TMath::Sin( 2.*values[AliDielectronVarManager::kDeltaPhiTPCrpH2] );

  //calculate inner product of strong Mag and ee plane
  if(Req(kPairPlaneMagInPro)) values[AliDielectronVarManager::kPairPlaneMagInPro] = pair->PairPlaneMagInnerProduct(values[AliDielectronVarManager::kZDCACrpH1]);

  //Calculate the angle between electrons decay plane and variables 1-4
  if(Req(kPairPlaneAngle1A)) values[AliDielectronVarManager::kPairPlaneAngle1A] = pair->GetPairPlaneAngle(values[kv0ArpH2],1);
  if(Req(kPairPlaneAngle2A)) values[AliDielectronVarManager::kPairPlaneAngle2A] = pair->GetPairPlaneAngle(values[kv0ArpH2],2);
  if(Req(kPairPlaneAngle3A)) values[AliDielectronVarManager::kPairPlaneAngle3A] = pair->GetPairPlaneAngle(values[kv0ArpH2],3);
  if(Req(kPairPlaneAngle4A)) values[AliDielectronVarManager::kPairPlaneAngle4A] = pair->GetPairPlaneAngle(values[kv0ArpH2],4);

  if(Req(kPairPlaneAngle1C)) values[AliDielectronVarManager::kPairPlaneAngle1C] = pair->GetPairPlaneAngle(values[kv0CrpH2],1);
  if(Req(kPairPlaneAngle2C)) values[AliDielectronVarManager::kPairPlaneAngle2C] = pair->GetPairPlaneAngle(values[kv0CrpH2],2);
  if(Req(kPairPlaneAngle3C)) values[AliDielectronVarManager::kPairPlaneAngle3C] = pair->GetPairPlaneAngle(values[kv0CrpH2],3);
  if(Req(kPairPlaneAngle4C)) values[AliDielectronVarManager::kPairPlaneAngle4C] = pair->GetPairPlaneAngle(values[kv0CrpH2],4);

  if(Req(kPairPlaneAngle1AC)) values[AliDielectronVarManager::kPairPlaneAngle1AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],1);
  if(Req(kPairPlaneAngle2AC)) values[AliDielectronVarManager::kPairPlaneAngle2AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],2);
  if(Req(kPairPlaneAngle3AC)) values[AliDielectronVarManager::kPairPlaneAngle3AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],3);
  if(Req(kPairPlaneAngle4AC)) values[AliDielectronVarManager::kPairPlaneAngle4AC] = pair->GetPairPlaneAngle(values[kv0ACrpH2],4);

  //Random reaction plane
  values[AliDielectronVarManager::kRandomRP] = gRandom->Uniform(-TMath::Pi()/2.0,TMath::Pi()/2.0);
  //delta phi of pair fron random reaction plane
  values[AliDielectronVarManager::kDeltaPhiRandomRP] = phi - values[kRandomRP];
  // keep the interval [-pi,+pi]
  if ( values[AliDielectronVarManager::kDeltaPhiRandomRP] > TMath::Pi() )
    values[AliDielectronVarManager::kDeltaPhiRandomRP] -= TMath::TwoPi();

  if(Req(kPairPlaneAngle1Ran)) values[AliDielectronVarManager::kPairPlaneAngle1Ran]= pair->GetPairPlaneAngle(values[kRandomRP],1);
  if(Req(kPairPlaneAngle2Ran)) values[AliDielectronVarManager::kPairPlaneAngle2Ran]= pair->GetPairPlaneAngle(values[kRandomRP],2);
  if(Req(kPairPlaneAngle3Ran)) values[AliDielectronVarManager::kPairPlaneAngle3Ran]= pair->GetPairPlaneAngle(values[kRandomRP],3);
  if(Req(kPairPlaneAngle4Ran)) values[AliDielectronVarManager::kPairPlaneAngle4Ran]= pair->GetPairPlaneAngle(values[kRandomRP],4);

  // Calculate v2 of Jpsi using the EP from the 2016 est. qVecQnFramework
  Double_t qnTPCeventplane = values[AliDielectronVarManager::kQnTPCrpH2];
  if(fgEventPlaneACremoval)
    if(fgQnEPacRemoval->IsSelected(pair)){
      AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
      if( AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*> (man->GetTask("FlowQnVectorCorrections")) ){
        if(flowQnVectorTask != NULL){
          AliQnCorrectionsManager *flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
          TList *qnlist = flowQnVectorMgr->GetQnVectorList();
          if(qnlist != NULL){
            qnTPCeventplane = fgQnEPacRemoval->GetACcorrectedQnTPCEventplane(pair, qnlist); // Remove auto correlations from the eventplane for the given pair
          }
          if(qnTPCeventplane == -999.) qnTPCeventplane = values[AliDielectronVarManager::kQnTPCrpH2];
        }
      }
    }

  if(Req(kQnDeltaPhiTPCrpH2) || Req(kQnTPCrpH2FlowV2))   values[AliDielectronVarManager::kQnDeltaPhiTPCrpH2]  = TVector2::Phi_mpi_pi(phi - qnTPCeventplane);
  if(Req(kQnDeltaPhiV0ArpH2) || Req(kQnV0ArpH2FlowV2))   values[AliDielectronVarManager::kQnDeltaPhiV0ArpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kQnV0ArpH2]);
  if(Req(kQnDeltaPhiV0CrpH2) || Req(kQnV0CrpH2FlowV2))   values[AliDielectronVarManager::kQnDeltaPhiV0CrpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kQnV0CrpH2]);
  if(Req(kQnDeltaPhiV0rpH2) || Req(kQnV0rpH2FlowV2))   values[AliDielectronVarManager::kQnDeltaPhiV0rpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kQnV0rpH2]);
  if(Req(kQnDeltaPhiSPDrpH2) || Req(kQnSPDrpH2FlowV2))   values[AliDielectronVarManager::kQnDeltaPhiSPDrpH2]  = TVector2::Phi_mpi_pi(phi - values[AliDielectronVarManager::kQnSPDrpH2]);
  if(Req(kQnTPCrpH2FlowV2)) values[AliDielectronVarManager::kQnTPCrpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kQnDeltaPhiTPCrpH2] );
  if(Req(kQnV0ArpH2FlowV2)) values[AliDielectronVarManager::kQnV0ArpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kQnDeltaPhiV0ArpH2] );
  if(Req(kQnV0CrpH2FlowV2)) values[AliDielectronVarManager::kQnV0CrpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kQnDeltaPhiV0CrpH2] );
  if(Req(kQnV0rpH2FlowV2)) values[AliDielectronVarManager::kQnV0rpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kQnDeltaPhiV0rpH2] );
  if(Req(kQnSPDrpH2FlowV2)) values[AliDielectronVarManager::kQnSPDrpH2FlowV2]    = TMath::Cos( 2.*values[AliDielectronVarManager::kQnDeltaPhiSPDrpH2] );

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
      if(pair->GetFirstDaughterP()->GetLabel() > 0) {
        const AliVParticle *motherMC = 0x0;
        if(fgEvent->IsA() == AliESDEvent::Class())  motherMC = (AliMCParticle*)mc->GetMCTrackMother((AliESDtrack*)pair->GetFirstDaughterP());
        else if(fgEvent->IsA() == AliAODEvent::Class())  motherMC = (AliAODMCParticle*)mc->GetMCTrackMother((AliAODTrack*)pair->GetFirstDaughterP());
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
	  AliVParticle* leg1 = pair->GetFirstDaughterP();
	  AliVParticle* leg2 = pair->GetSecondDaughterP();
	  if (leg1 && leg2){
		Fill(leg1, valuesLeg1);
		Fill(leg2, valuesLeg2);
		values[AliDielectronVarManager::kTRDpidEffPair] = valuesLeg1[AliDielectronVarManager::kTRDpidEffLeg]*valuesLeg2[AliDielectronVarManager::kTRDpidEffLeg];
	  }
	}


  }//if (mc->HasMC())

  AliVParticle* leg1 = pair->GetFirstDaughterP();
  AliVParticle* leg2 = pair->GetSecondDaughterP();
  if (leg1)
    values[AliDielectronVarManager::kMomAsymDau1] = (values[AliDielectronVarManager::kP] != 0)? leg1->P()  / values[AliDielectronVarManager::kP]: 0;
  else
    values[AliDielectronVarManager::kMomAsymDau1] = -9999.;
  if (leg2)
    values[AliDielectronVarManager::kMomAsymDau2] = (values[AliDielectronVarManager::kP] != 0)? leg2->P()  / values[AliDielectronVarManager::kP]: 0;
  else
    values[AliDielectronVarManager::kMomAsymDau2] = -9999.;

  Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
  Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
  values[AliDielectronVarManager::kPairEff]=0.0;
  values[AliDielectronVarManager::kOneOverPairEff]=0.0;
  values[AliDielectronVarManager::kOneOverPairEffSq]=0.0;
  if (leg1 && leg2 && fgLegEffMap) {
    Fill(leg1, valuesLeg1);
    Fill(leg2, valuesLeg2);
    values[AliDielectronVarManager::kPairEff] = valuesLeg1[AliDielectronVarManager::kLegEff] *valuesLeg2[AliDielectronVarManager::kLegEff];
  }
  else if(fgPairEffMap) {
    values[AliDielectronVarManager::kPairEff] = GetPairEff(values);
  }
  if(fgLegEffMap || fgPairEffMap) {
    values[AliDielectronVarManager::kOneOverPairEff] = (values[AliDielectronVarManager::kPairEff]>0.0 ? 1./values[AliDielectronVarManager::kPairEff] : 1.0);
    values[AliDielectronVarManager::kOneOverPairEffSq] = (values[AliDielectronVarManager::kPairEff]>0.0 ? 1./values[AliDielectronVarManager::kPairEff]/values[AliDielectronVarManager::kPairEff] : 1.0);
  }

  if(kRndmPair) values[AliDielectronVarManager::kRndmPair] = gRandom->Rndm();
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
  values[AliDielectronVarManager::kPtSq]      = particle->GetPt() * particle->GetPt();
  values[AliDielectronVarManager::kP]         = particle->GetP();

  values[AliDielectronVarManager::kXv]        = particle->GetX();
  values[AliDielectronVarManager::kYv]        = particle->GetY();
  values[AliDielectronVarManager::kZv]        = particle->GetZ();

  values[AliDielectronVarManager::kOneOverPt] = 0;
  values[AliDielectronVarManager::kPhi]       = TVector2::Phi_0_2pi(particle->GetPhi()); // interval [0,+2pi]
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
  values[AliDielectronVarManager::kTPCnSigmaEleRaw]  = 0;
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
    if(fgZDCRecenteringFile.Contains(".root")) InitZDCRecenteringHistograms(event->GetRunNumber());
    fgCurrentRun=event->GetRunNumber();
  }
  values[AliDielectronVarManager::kMixingBin]=0;

  values[AliDielectronVarManager::kXvPrim]       = 0;
  values[AliDielectronVarManager::kYvPrim]       = 0;
  values[AliDielectronVarManager::kZvPrim]       = 0;
  values[AliDielectronVarManager::kNVtxContrib]  = 0;
//   values[AliDielectronVarManager::kChi2NDF]      = 0; //This is the pair value!!!

  values[AliDielectronVarManager::kNTrk]            = 0;
  values[AliDielectronVarManager::kNVtxContrib]     = 0;
  values[AliDielectronVarManager::kNacc]            = 0;
  values[AliDielectronVarManager::kMatchEffITSTPC]  = 0;
  values[AliDielectronVarManager::kNaccTrcklts]     = 0;
  values[AliDielectronVarManager::kNaccTrcklts09]   = 0;
  values[AliDielectronVarManager::kNaccTrcklts10]   = 0;
  values[AliDielectronVarManager::kNaccTrcklts0916] = 0;
  values[AliDielectronVarManager::kNevents]         = 0; //always fill bin 0;
  values[AliDielectronVarManager::kRefMult]         = 0;
  values[AliDielectronVarManager::kRefMultTPConly]  = 0;

  //Centrality Run2 - from 2015 on
  values[AliDielectronVarManager::kCentralityNew] = 0.;
  values[AliDielectronVarManager::kCentralityCL0] = 0.;
  values[AliDielectronVarManager::kCentralityCL1] = 0.;
  values[AliDielectronVarManager::kCentralitySPDClusters]  = 0.;
  values[AliDielectronVarManager::kCentralitySPDTracklets] = 0.;
  values[AliDielectronVarManager::kCentralityCL0plus05]    = 0.;
  values[AliDielectronVarManager::kCentralityCL0minus05]   = 0.;
  values[AliDielectronVarManager::kCentralityCL0plus10]    = 0.;
  values[AliDielectronVarManager::kCentralityCL0minus10]   = 0.;

  AliMultSelection *multSelection = (AliMultSelection*)event->FindListObject("MultSelection");
  if(!multSelection){
    // only for debugging purpose
    // printf("Did not find AliMultSelection!!! Mostly defined for run2 data");
  }
  else {
    values[AliDielectronVarManager::kCentralityNew]          = multSelection->GetMultiplicityPercentile("V0M",kFALSE);
    values[AliDielectronVarManager::kCentralityCL0]          = multSelection->GetMultiplicityPercentile("CL0",kFALSE);
    values[AliDielectronVarManager::kCentralityCL1]          = multSelection->GetMultiplicityPercentile("CL1",kFALSE);
    values[AliDielectronVarManager::kCentralitySPDClusters]  = multSelection->GetMultiplicityPercentile("SPDClustersCorr",kFALSE);
    values[AliDielectronVarManager::kCentralitySPDTracklets] = multSelection->GetMultiplicityPercentile("SPDTrackletsCorr",kFALSE);
    values[AliDielectronVarManager::kCentralityCL0plus05]    = multSelection->GetMultiplicityPercentile("CL0plus05",kFALSE);
    values[AliDielectronVarManager::kCentralityCL0minus05]   = multSelection->GetMultiplicityPercentile("CL0minus05",kFALSE);
    values[AliDielectronVarManager::kCentralityCL0plus10]    = multSelection->GetMultiplicityPercentile("CL0plus10",kFALSE);
    values[AliDielectronVarManager::kCentralityCL0minus10]   = multSelection->GetMultiplicityPercentile("CL0minus10",kFALSE);
  }
  const AliVVertex *primVtx = event->GetPrimaryVertex();
  if (primVtx){
    //    printf("prim vertex reco: %f \n",primVtx->GetX());
    values[AliDielectronVarManager::kXvPrim]       = primVtx->GetX();
    values[AliDielectronVarManager::kYvPrim]       = primVtx->GetY();
    values[AliDielectronVarManager::kZvPrim]       = primVtx->GetZ();
    values[AliDielectronVarManager::kNVtxContrib]  = primVtx->GetNContributors();
  }
  //   values[AliDielectronVarManager::kChi2NDF]      = primVtx->GetChi2perNDF(); //this is the pair value

  // online and offline trigger maps
  values[AliDielectronVarManager::kTriggerInclONL]     = event->GetTriggerMask();
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  UInt_t maskOff = ((AliInputEventHandler*)man->GetInputEventHandler())->IsEventSelected();
  values[AliDielectronVarManager::kTriggerInclOFF]     = maskOff;
  values[AliDielectronVarManager::kTriggerExclOFF]        = -1;
  for(Int_t i=0; i<30; i++) { if(maskOff==BIT(i)) values[AliDielectronVarManager::kTriggerExclOFF]=i; }

  values[AliDielectronVarManager::kNTrk]            = event->GetNumberOfTracks();
  if(Req(kNacc))            values[AliDielectronVarManager::kNacc]            = AliDielectronHelper::GetNacc(event);
  if(Req(kMatchEffITSTPC))  values[AliDielectronVarManager::kMatchEffITSTPC]  = AliDielectronHelper::GetITSTPCMatchEff(event);
  if(Req(kNaccTrcklts) || Req(kNaccTrckltsCorr))
    values[AliDielectronVarManager::kNaccTrcklts]     = AliDielectronHelper::GetNaccTrcklts(event,1.6);
  if(Req(kNaccTrcklts09))
      values[AliDielectronVarManager::kNaccTrcklts09]     = AliDielectronHelper::GetNaccTrcklts(event,0.9);
  if(Req(kNaccTrcklts10) || Req(kNaccTrcklts10Corr))
    values[AliDielectronVarManager::kNaccTrcklts10]   = AliDielectronHelper::GetNaccTrcklts(event,1.0);
  if(Req(kNaccTrcklts0916))
    values[AliDielectronVarManager::kNaccTrcklts0916] = AliDielectronHelper::GetNaccTrcklts(event,1.6)-AliDielectronHelper::GetNaccTrcklts(event,.9);

  if(Req(kNaccTrckltsCorr))
  values[AliDielectronVarManager::kNaccTrckltsCorr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event, values[AliDielectronVarManager::kNaccTrcklts],
						 values[AliDielectronVarManager::kZvPrim],2);
  if(Req(kNaccTrcklts10Corr))
  values[AliDielectronVarManager::kNaccTrcklts10Corr] =
    AliDielectronHelper::GetNaccTrckltsCorrected(event, values[AliDielectronVarManager::kNaccTrcklts10],
						 values[AliDielectronVarManager::kZvPrim],1);


  Double_t ptMaxEv    = -1., phiptMaxEv= -1.;
  if(Req(kMaxPt) || Req(kPhiMaxPt)) AliDielectronHelper::GetMaxPtAndPhi(event, ptMaxEv, phiptMaxEv);
  values[AliDielectronVarManager::kPhiMaxPt]          = phiptMaxEv;
  values[AliDielectronVarManager::kMaxPt]             = ptMaxEv;


  // event plane quantities from the AliEPSelectionTask
  for(Int_t ivar=AliDielectronVarManager::kv0ArpH2; ivar<=kv0C0v0C3DiffH2;   ivar++) values[ivar] = 0.0; // v0  variables
  for(Int_t ivar=AliDielectronVarManager::kTPCxH2;  ivar<=kTPCsub12DiffH2uc; ivar++) values[ivar] = 0.0; // tpc variables

  // If QnCorrections framework (est. 2016) task should be used run the AddTask with your train. The following function will overwrite the existing AliEventplane object with the information extracted from the QnCorrections Task. Then the following code can be used as usual. The current implementation uses only the second harmonic but this could be adaptet if needed.
  if( AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
      dynamic_cast<AliAnalysisTaskFlowVectorCorrections*> (man->GetTask("FlowQnVectorCorrections")) )
  if(flowQnVectorTask != NULL){
    AliQnCorrectionsManager *flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    TList *qnlist = flowQnVectorMgr->GetQnVectorList();
    if(qnlist != NULL)  AliDielectronVarManager::FillQnEventplanes(qnlist, values);
  }

  // ep angle interval [todo, fill]
  AliEventplane *ep = const_cast<AliVEvent*>(event)->GetEventplane();
  if(ep) {

    // TPC event plane quantities (uncorrected)
    TVector2 *qstd  = ep->GetQVector();  // This is the "standard" Q-Vector for TPC
    TVector2 *qsub1 = ep->GetQsub1();    // random subevent plane
    TVector2 *qsub2 = ep->GetQsub2();

    if(qstd) {
      values[AliDielectronVarManager::kTPCxH2uc]       = qstd->X();
      values[AliDielectronVarManager::kTPCyH2uc]       = qstd->Y();
      values[AliDielectronVarManager::kTPCmagH2uc]     = qstd->Mod();
      values[AliDielectronVarManager::kTPCrpH2uc]      = TVector2::Phi_mpi_pi(qstd->Phi())/2;
      if(qsub1 && qsub2){
        values[AliDielectronVarManager::kTPCsub1xH2uc]   = qsub1->X();
        values[AliDielectronVarManager::kTPCsub1yH2uc]   = qsub1->Y();
        values[AliDielectronVarManager::kTPCsub1rpH2uc]  = TVector2::Phi_mpi_pi(qsub1->Phi())/2;
        values[AliDielectronVarManager::kTPCsub2xH2uc]   = qsub2->X();
        values[AliDielectronVarManager::kTPCsub2yH2uc]   = qsub2->Y();
        values[AliDielectronVarManager::kTPCsub2rpH2uc]  = TVector2::Phi_mpi_pi(qsub2->Phi())/2;

        values[AliDielectronVarManager::kTPCsub12DiffH2uc] = TMath::Cos( 2.*(values[AliDielectronVarManager::kTPCsub1rpH2uc] -
  									   values[AliDielectronVarManager::kTPCsub2rpH2uc]) );
      }
    }

    // VZERO event plane
    TVector2 qvec;
    TVector2 *qVecQnFramework;
    Double_t qx = 0, qy = 0;

    ep->CalculateVZEROEventPlane(event,10, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ACrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliDielectronVarManager::kv0ACxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0ACyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0ACmagH2] = qvec.Mod();

    qx = qy = 0.;
    ep->CalculateVZEROEventPlane(event, 8, 2, qy, qy);
    qvec.Set(qx,qy);

    values[AliDielectronVarManager::kv0ArpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliDielectronVarManager::kv0AxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0AyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0AmagH2] = qvec.Mod();

    qx = qy = 0.;
    ep->CalculateVZEROEventPlane(event, 9, 2, qx, qy);
    qvec.Set(qx,qy);

    values[AliDielectronVarManager::kv0CrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliDielectronVarManager::kv0CxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0CyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0CmagH2] = qvec.Mod();
    ep->CalculateVZEROEventPlane(event, 0, 0, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep->CalculateVZEROEventPlane(event, 3, 3, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep->CalculateVZEROEventPlane(event, 4, 4, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep->CalculateVZEROEventPlane(event, 7, 7, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
  } //if: eventplane

  // ESD VZERO information
  AliVVZERO* vzeroData = event->GetVZEROData();
  values[AliDielectronVarManager::kMultV0A] = 0.0;
  values[AliDielectronVarManager::kMultV0C] = 0.0;
  values[AliDielectronVarManager::kEqMultV0A] = 0.0;
  values[AliDielectronVarManager::kEqMultV0C] = 0.0;
  values[AliDielectronVarManager::kAdcV0A]  = 0.0;
  values[AliDielectronVarManager::kAdcV0C]  = 0.0;
  for(Int_t i=0; i<32; ++i) {
    values[AliDielectronVarManager::kVZEROchMult+i] = vzeroData->GetMultiplicity(i);
    values[AliDielectronVarManager::kVZEROchMult+32+i] = vzeroData->GetMultiplicity(i+32);
    //values[AliDielectronVarManager::kVZEROchMult+i] = event->GetVZEROEqMultiplicity(i);
    //values[AliDielectronVarManager::kVZEROchMult+32+i] = event->GetVZEROEqMultiplicity(i+32);
    values[AliDielectronVarManager::kMultV0A] += vzeroData->GetMultiplicityV0A(i);
    values[AliDielectronVarManager::kMultV0C] += vzeroData->GetMultiplicityV0C(i);
    values[AliDielectronVarManager::kEqMultV0A] += event->GetVZEROEqMultiplicity(i);
    values[AliDielectronVarManager::kEqMultV0C] += event->GetVZEROEqMultiplicity(i+32);
    //values[AliDielectronVarManager::kAdcV0A] += vzeroData->GetAdcV0A(i);
    //values[AliDielectronVarManager::kAdcV0C] += vzeroData->GetAdcV0C(i);
  }
  values[AliDielectronVarManager::kMultV0] = values[AliDielectronVarManager::kMultV0A] + values[AliDielectronVarManager::kMultV0C];
  values[AliDielectronVarManager::kEqMultV0] = values[AliDielectronVarManager::kEqMultV0A] + values[AliDielectronVarManager::kEqMultV0C];
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

  Double_t ZDCqvec[3][2] = {{999., 999.}, {999., 999.}, {999., 999.} };
  GetZDCRP(event, ZDCqvec);

  values[AliDielectronVarManager::kZDCArpH1] = TMath::ATan2(ZDCqvec[0][1], ZDCqvec[0][0]);
  values[AliDielectronVarManager::kZDCCrpH1] = TMath::ATan2(ZDCqvec[1][1], ZDCqvec[1][0]);
  values[AliDielectronVarManager::kZDCACrpH1] = TMath::ATan2(ZDCqvec[2][1], ZDCqvec[2][0]);

  if(TMath::Abs(ZDCqvec[0][0] - 999.) < 1e-10 || TMath::Abs(ZDCqvec[0][1] - 999.) < 1e-10 || TMath::Abs(ZDCqvec[1][0] - 999.) < 1e-10 || TMath::Abs(ZDCqvec[1][1] - 999.) < 1e-10){
    values[AliDielectronVarManager::kZDCArpH1] = 999;
    values[AliDielectronVarManager::kZDCCrpH1] = 999;
    values[AliDielectronVarManager::kZDCACrpH1] = 999;
  }



  values[AliDielectronVarManager::kv0ZDCrpRes] = cos(2*(values[AliDielectronVarManager::kZDCArpH1] - values[AliDielectronVarManager::kv0ArpH2]));
  values[AliDielectronVarManager::kZDCrpResH1] = cos(values[AliDielectronVarManager::kZDCArpH1] - values[AliDielectronVarManager::kZDCCrpH1]);


}

inline void AliDielectronVarManager::FillVarESDEvent(const AliESDEvent *event, Double_t * const values)
{
  //
  // Fill event information available for histogramming into an array
  //

  // Fill common AliVEvent interface information
  FillVarVEvent(event, values);

  // Centrality Run1
  Double_t centralityF=-1;
  Double_t centralitySPD=-1;
  Double_t centralityV0A = -1;
  Double_t centralityV0C = -1;
  Double_t centralityZNA = -1;

  AliCentrality *esdCentrality = const_cast<AliESDEvent*>(event)->GetCentrality();
  if(esdCentrality) {
    centralityF   = esdCentrality->GetCentralityPercentile("V0M");
    centralitySPD = esdCentrality->GetCentralityPercentile("CL1");
    centralityV0A = esdCentrality->GetCentralityPercentile("V0A");
    centralityV0C = esdCentrality->GetCentralityPercentile("V0C");
    centralityZNA = esdCentrality->GetCentralityPercentile("ZNA");
  }

  // Fill AliESDEvent interface specific information
  const AliESDVertex *primVtx = event->GetPrimaryVertex();
  values[AliDielectronVarManager::kXRes]       = primVtx->GetXRes();
  values[AliDielectronVarManager::kYRes]       = primVtx->GetYRes();
  values[AliDielectronVarManager::kZRes]       = primVtx->GetZRes();
  values[AliDielectronVarManager::kCentrality] = centralityF;
  values[AliDielectronVarManager::kCentralitySPD] = centralitySPD;
  values[AliDielectronVarManager::kCentralityV0A] = centralityV0A;
  values[AliDielectronVarManager::kCentralityV0C] = centralityV0C;
  values[AliDielectronVarManager::kCentralityZNA] = centralityZNA;

  const AliESDVertex *vtxTPC = event->GetPrimaryVertexTPC();
  values[AliDielectronVarManager::kNVtxContribTPC] = (vtxTPC ? vtxTPC->GetNContributors() : 0);

  // The true vertex is needed for the pair DCA analysis (needs DCA of reco track w.r.t. true vertex).
  if (AliDielectronMC::Instance()->HasMC()){
    if (Req(kDistPrimToSecVtxXYMC) || Req(kDistPrimToSecVtxZMC) || Req(kXvPrimMCtruth) || Req(kYvPrimMCtruth) || Req(kZvPrimMCtruth)) {
      AliMCEvent* mcevent = AliDielectronMC::Instance()->GetMCEvent();
      const AliVVertex* mcvtx = mcevent->GetPrimaryVertex();
      values[AliDielectronVarManager::kXvPrimMCtruth] = (mcvtx ? mcvtx->GetX() : 0.0);
      values[AliDielectronVarManager::kYvPrimMCtruth] = (mcvtx ? mcvtx->GetY() : 0.0);
      values[AliDielectronVarManager::kZvPrimMCtruth] = (mcvtx ? mcvtx->GetZ() : 0.0);
    }
  }

  // Event multiplicity estimators
  Int_t nTrSPD05=0; Int_t nTrITSTPC05=0; Int_t nTrITSSA05=0;
  nTrSPD05    = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTracklets, 0.5);
  nTrITSTPC05 = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSTPC, 0.5);
  nTrITSSA05  = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSSA, 0.5);
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
  nTrSPD10    = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTracklets, 1.0);
  nTrITSTPC10 = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSTPC, 1.0);
  nTrITSSA10  = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSSA, 1.0);
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
  nTrSPD16    = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTracklets, 1.6);
  nTrITSTPC16 = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSTPC, 1.6);
  nTrITSSA16  = AliESDtrackCuts::GetReferenceMultiplicity(event, AliESDtrackCuts::kTrackletsITSSA, 1.6);
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
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  assert(header&&"Not a standard AOD");


  // Centrality Run1 / Run2
  Double_t centralityF=-1;
  Double_t centralitySPD=-1;
  Double_t centralityV0A = -1;
  Double_t centralityV0C = -1;
  Double_t centralityZNA = -1;

  if(AliMultSelection *multSelection = (AliMultSelection*) event->FindListObject("MultSelection")){
    centralityF   = multSelection->GetMultiplicityPercentile("V0M",kFALSE);
    centralitySPD = multSelection->GetMultiplicityPercentile("CL1",kFALSE);
    centralityV0A = multSelection->GetMultiplicityPercentile("V0A",kFALSE); // Currently not supported for 15o
    centralityV0C = multSelection->GetMultiplicityPercentile("V0C",kFALSE); // Currently not supported for 15o
    centralityZNA = multSelection->GetMultiplicityPercentile("ZDC",kFALSE);
  }
  else{
    AliCentrality *aodCentrality = header->GetCentralityP();
    if (aodCentrality) centralityF = aodCentrality->GetCentralityPercentile("V0M");
    if (aodCentrality) centralitySPD = aodCentrality->GetCentralityPercentile("CL1");
    if (aodCentrality) centralityV0A = aodCentrality->GetCentralityPercentile("V0A");
    if (aodCentrality) centralityV0C = aodCentrality->GetCentralityPercentile("V0C");
    if (aodCentrality) centralityZNA = aodCentrality->GetCentralityPercentile("ZNA");
  }

  values[AliDielectronVarManager::kCentrality] = centralityF;
  values[AliDielectronVarManager::kCentralitySPD] = centralitySPD;
  values[AliDielectronVarManager::kCentralityV0A] = centralityV0A;
  values[AliDielectronVarManager::kCentralityV0C] = centralityV0C;
  values[AliDielectronVarManager::kCentralityZNA] = centralityZNA;

  values[AliDielectronVarManager::kRefMult]        = header->GetRefMultiplicity();        // similar to Ntrk
  values[AliDielectronVarManager::kRefMultTPConly] = header->GetTPConlyRefMultiplicity(); // similar to Nacc
  values[AliDielectronVarManager::kRefMultOvRefMultTPConly] = (values[AliDielectronVarManager::kRefMultTPConly] > 0. ? (values[AliDielectronVarManager::kRefMult]/values[AliDielectronVarManager::kRefMultTPConly]) : 0.);

  // The true vertex is needed for the pair DCA analysis (needs DCA of reco track w.r.t. true vertex).
  if (AliDielectronMC::Instance()->HasMC()){
    if (Req(kDistPrimToSecVtxXYMC) || Req(kDistPrimToSecVtxZMC) || Req(kXvPrimMCtruth) || Req(kYvPrimMCtruth) || Req(kZvPrimMCtruth)) {
      // @TODO: adopt the code from FillVarESDEvent() for AOD...
      printf("WARNING: filling of MC true vertex not implemented for AOD tracks!\n");
      values[AliDielectronVarManager::kXvPrimMCtruth] = 0.;
      values[AliDielectronVarManager::kYvPrimMCtruth] = 0.;
      values[AliDielectronVarManager::kZvPrimMCtruth] = 0.;
    }
  }

  ///////////////////////////////////////////
  //////////// NANO AODs ////////////////////
  ///////////////////////////////////////////

  // (w/o AliCentrality branch), VOM centrality should be stored in the header
  if((!header->GetCentralityP()) && (centralityF == -1))
    values[AliDielectronVarManager::kCentrality] = header->GetCentrality();
  // (w/o AliEventPlane branch) tpc event plane stuff stored in the header
  if(!header->GetEventplaneP()) {
    //    values[AliDielectronVarManager::kNTrk] = header->GetRefMultiplicity();    // overwritten datamembers in "our" nanoAODs
    //    values[AliDielectronVarManager::kNacc] = header->GetRefMultiplicityPos(); // overwritten datamembers in "our" nanoAODs

    TVector2 qvec;
    // TPC

    TList *qnlist = (TList*) event->FindListObject("qnVectorList");
    if(Req(kQnTPCrpH2) && qnlist ==NULL){
      for (Int_t i = AliDielectronVarManager::kQnTPCrpH2; i <= AliDielectronVarManager::kQnCorrFMDAy_FMDCy; i++) {
        values[i] = -999.;
      }
    }
    if(qnlist != NULL)  AliDielectronVarManager::FillQnEventplanes(qnlist, values);

    qvec.Set(header->GetEventplaneQx(), header->GetEventplaneQy());
    values[AliDielectronVarManager::kTPCxH2uc]   = qvec.X();
    values[AliDielectronVarManager::kTPCyH2uc]   = qvec.Y();
    values[AliDielectronVarManager::kTPCmagH2uc] = qvec.Mod();
    values[AliDielectronVarManager::kTPCrpH2uc]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;

    // VZERO
    AliEventplane ep2;
    // get event plane corrections from the VZERO EP selection task
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliVZEROEPSelectionTask *eptask = dynamic_cast<AliVZEROEPSelectionTask *>(man->GetTask("AliVZEROEPSelectionTask"));
    if(eptask) eptask->SetEventplaneParams(&ep2,centralityF);
    else if(!qnlist) printf("no VZERO event plane selection task added! \n");


    Double_t qx = 0, qy = 0;
    ep2.CalculateVZEROEventPlane(event,10, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ACrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliDielectronVarManager::kv0ACxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0ACyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0ACmagH2] = qvec.Mod();



    ep2.CalculateVZEROEventPlane(event, 8, 2, qx, qy);
    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0ArpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliDielectronVarManager::kv0AxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0AyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0AmagH2] = qvec.Mod();

    ep2.CalculateVZEROEventPlane(event, 9, 2, qx, qy);
    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0CrpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    values[AliDielectronVarManager::kv0CxH2]   = qvec.X();
    values[AliDielectronVarManager::kv0CyH2]   = qvec.Y();
    values[AliDielectronVarManager::kv0CmagH2] = qvec.Mod();
    ep2.CalculateVZEROEventPlane(event, 0, 0, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep2.CalculateVZEROEventPlane(event, 3, 3, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0C3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep2.CalculateVZEROEventPlane(event, 4, 4, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A0rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;
    ep2.CalculateVZEROEventPlane(event, 7, 7, 2, qx, qy);    qvec.Set(qx,qy);
    values[AliDielectronVarManager::kv0A3rpH2]  = TVector2::Phi_mpi_pi(qvec.Phi())/2;

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
  // For MC truth, these variables are identical to the above. (different in FillVarESDEvent() / FillVarAODEvent()).
  values[AliDielectronVarManager::kXvPrimMCtruth]       = values[AliDielectronVarManager::kXvPrim];
  values[AliDielectronVarManager::kYvPrimMCtruth]       = values[AliDielectronVarManager::kYvPrim];
  values[AliDielectronVarManager::kZvPrimMCtruth]       = values[AliDielectronVarManager::kZvPrim];
  // Fill AliMCEvent interface specific information
  values[AliDielectronVarManager::kNch]   = AliDielectronHelper::GetNch(event, 1.6);
  values[AliDielectronVarManager::kNchJpsiExcl]   = AliDielectronHelper::GetNch(event, 1.6, kTRUE);
  values[AliDielectronVarManager::kNch05] = AliDielectronHelper::GetNch(event, 0.5);
  values[AliDielectronVarManager::kNch05JpsiExcl] = AliDielectronHelper::GetNch(event, 0.5, kTRUE);
  values[AliDielectronVarManager::kNch10] = AliDielectronHelper::GetNch(event, 1.0);
  values[AliDielectronVarManager::kNch10JpsiExcl] = AliDielectronHelper::GetNch(event, 1.0, kTRUE);

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
      values[AliDielectronVarManager::kTPCrpH2]  = TVector2::Phi_mpi_pi(qcorr->Phi())/2;
      // detector effects
      values[AliDielectronVarManager::kCosTPCrpH2]     = TMath::Cos( 2.* values[AliDielectronVarManager::kTPCrpH2] );
      values[AliDielectronVarManager::kSinTPCrpH2]     = TMath::Sin( 2.* values[AliDielectronVarManager::kTPCrpH2] );

      // correlations for event plane resoultion
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
    if(qcsub1 && qcsub2) {
      values[AliDielectronVarManager::kTPCsub1xH2]   = qcsub1->X();
      values[AliDielectronVarManager::kTPCsub1yH2]   = qcsub1->Y();
      values[AliDielectronVarManager::kTPCsub1rpH2]  = TVector2::Phi_mpi_pi(qcsub1->Phi())/2;

      values[AliDielectronVarManager::kTPCsub2xH2]   = qcsub2->X();
      values[AliDielectronVarManager::kTPCsub2yH2]   = qcsub2->Y();
      values[AliDielectronVarManager::kTPCsub2rpH2]  = TVector2::Phi_mpi_pi(qcsub2->Phi())/2;

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


inline void AliDielectronVarManager::InitEstimatorAvg(const Char_t* filename) //Not Grid compatible
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


inline void AliDielectronVarManager::InitEstimatorObjArrayAvg(const TObjArray* array) //Grid compatible
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

inline Double_t AliDielectronVarManager::GetSingleLegEff(Double_t * const values) {
  //
  // get the single leg efficiency for a given particle
  //
  if(!fgLegEffMap) return -1.;

  if(fgLegEffMap->InheritsFrom(THnBase::Class())) {
    THnBase *eff = static_cast<THnBase*>(fgLegEffMap);
    Int_t dim=eff->GetNdimensions();
    Int_t idx[dim];
    for(Int_t idim=0; idim<dim; idim++) {
      UInt_t var = GetValueType(eff->GetAxis(idim)->GetName());
      idx[idim] = eff->GetAxis(idim)->FindBin(values[var]);
      if(idx[idim] < 0 || idx[idim]>eff->GetAxis(idim)->GetNbins()) return 0.0;
    }
    const Double_t ret=(eff->GetBinContent(idx));
    return ret;
  }
  return -1.;
}

inline Double_t AliDielectronVarManager::GetPairEff(Double_t * const values) {
  //
  // get the pair efficiency for given pair kinematics
  //
  if(!fgPairEffMap) return -1.;

  if(fgPairEffMap->IsA()== THnBase::Class()) {
    THnBase *eff = static_cast<THnBase*>(fgPairEffMap);
    Int_t dim=eff->GetNdimensions();
    Int_t idx[dim];
    for(Int_t idim=0; idim<dim; idim++) {
      UInt_t var = GetValueType(eff->GetAxis(idim)->GetName());
      idx[idim] = eff->GetAxis(idim)->FindBin(values[var]);
      if(idx[idim] < 0 || idx[idim]>eff->GetAxis(idim)->GetNbins()) return 0.0;
    }
    const Double_t ret=(eff->GetBinContent(idx));
    return ret;
  }
  if(fgPairEffMap->IsA()== TSpline3::Class()) {
    TSpline3 *eff = static_cast<TSpline3*>(fgPairEffMap);
    if(!eff->GetHistogram()) { printf("no histogram added to the spline\n"); return -1.;}
    UInt_t var = GetValueType(eff->GetHistogram()->GetXaxis()->GetName());
    return (eff->Eval(values[var]));
  }

  return -1.;
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

inline void AliDielectronVarManager::InitZDCRecenteringHistograms(Int_t runNo) {

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


//______________________________________________________________________________
inline Bool_t AliDielectronVarManager::GetDCA(const AliAODTrack *track, Double_t* d0z0, Double_t* covd0z0)
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
    //   }
   AliCentrality *esdCentrality = const_cast<AliESDEvent*>(esdEv)->GetCentrality();
    if(esdCentrality) centralitySPD = esdCentrality->GetCentralityPercentile("CL1");
    //2015 cent
    // AliMultSelection *MultSelection = (AliMultSelection*)const_cast<AliESDEvent*>(esdEv)->FindListObject("MultSelection");
    //if(MultSelection){
    // centralitySPD = MultSelection->GetMultiplicityPercentile("CL1",kFALSE);
    // }
  }
  if(event->IsA() == AliAODEvent::Class()) {
    const AliAODEvent* aodEv = static_cast<const AliAODEvent*>(event);
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(aodEv->GetHeader());
    assert(header&&"Not a standard AOD");
    AliCentrality *aodCentrality = header->GetCentralityP();
    // Run1 aodCentrality -- Run2 multSelection
    if(AliMultSelection *multSelection = (AliMultSelection*) aodEv->FindListObject("MultSelection")){
      centralitySPD = multSelection->GetMultiplicityPercentile("CL1",kFALSE);
    }
    else{
      if(aodCentrality) centralitySPD = aodCentrality->GetCentralityPercentile("CL1");
      else printf("GetVzeroRP: No centrality estimation avaible!\n");
    }
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
inline void AliDielectronVarManager::GetZDCRP(const AliVEvent* event, Double_t qvec[][2]) {

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

//________________________________________________________________
inline void AliDielectronVarManager::FillQnEventplanes(TList *qnlist, Double_t * const values){
  Bool_t bTPCqVector(kFALSE), bTPCaSideqVector(kFALSE), bTPCcSideqVector(kFALSE), bV0AqVector(kFALSE), bV0CqVector(kFALSE), bV0qVector(kFALSE),bSPDqVector(kFALSE), bFMDAqVector(kFALSE), bFMDCqVector(kFALSE);
  for (Int_t i = AliDielectronVarManager::kQnTPCrpH2; i <= AliDielectronVarManager::kQnCorrFMDAy_FMDCy; i++) {
    values[i] = -999.;
  }
  TString qnListDetector;
  // TPC Eventplane q-Vector
  qnListDetector = "TPC" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkTPC = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorTPC = new TVector2(-200.,-200.);
  if(qVecQnFrameworkTPC != NULL){
    bTPCqVector = kTRUE;
    qVectorTPC->Set(qVecQnFrameworkTPC->Qx(2),qVecQnFrameworkTPC->Qy(2));
    values[AliDielectronVarManager::kQnTPCrpH2] = TVector2::Phi_mpi_pi(qVectorTPC->Phi())/2;
    values[AliDielectronVarManager::kQnTPCxH2]  = qVecQnFrameworkTPC->Qx(2);
    values[AliDielectronVarManager::kQnTPCyH2]  = qVecQnFrameworkTPC->Qy(2);
  }
  delete qVectorTPC;

  // TPC A-Side/Neg. Eta Eventplane q-Vector
  qnListDetector = "TPCNegEta" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkTPCaSide = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorTPCaSide = new TVector2(-200.,-200.);
  if(qVecQnFrameworkTPCaSide != NULL){
    bTPCaSideqVector = kTRUE;
    qVectorTPCaSide->Set(qVecQnFrameworkTPCaSide->Qx(2),qVecQnFrameworkTPCaSide->Qy(2));
    values[AliDielectronVarManager::kQnTPCaSiderpH2] = TVector2::Phi_mpi_pi(qVectorTPCaSide->Phi())/2;
    values[AliDielectronVarManager::kQnTPCaSidexH2]  = qVecQnFrameworkTPCaSide->Qx(2);
    values[AliDielectronVarManager::kQnTPCaSideyH2]  = qVecQnFrameworkTPCaSide->Qy(2);
  }
  delete qVectorTPCaSide;

  // TPC C-Side/Pos. Eta Eventplane q-Vector
  qnListDetector = "TPCPosEta" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkTPCcSide = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorTPCcSide = new TVector2(-200.,-200.);
  if(qVecQnFrameworkTPCcSide != NULL){
    bTPCcSideqVector = kTRUE;
    qVectorTPCcSide->Set(qVecQnFrameworkTPCcSide->Qx(2),qVecQnFrameworkTPCcSide->Qy(2));
    values[AliDielectronVarManager::kQnTPCcSiderpH2] = TVector2::Phi_mpi_pi(qVectorTPCcSide->Phi())/2;
    values[AliDielectronVarManager::kQnTPCcSidexH2]  = qVecQnFrameworkTPCcSide->Qx(2);
    values[AliDielectronVarManager::kQnTPCcSideyH2]  = qVecQnFrameworkTPCcSide->Qy(2);
  }
  delete qVectorTPCcSide;

  // VZEROA Eventplane q-Vector
  qnListDetector = "VZEROA" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkV0A = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorV0A = new TVector2(-200.,-200.);
  if(qVecQnFrameworkV0A != NULL){
    bV0AqVector = kTRUE;
    qVectorV0A->Set(qVecQnFrameworkV0A->Qx(2),qVecQnFrameworkV0A->Qy(2));
    values[AliDielectronVarManager::kQnV0ArpH2] = TVector2::Phi_mpi_pi(qVectorV0A->Phi())/2;
    values[AliDielectronVarManager::kQnV0AxH2]  = qVecQnFrameworkV0A->Qx(2);
    values[AliDielectronVarManager::kQnV0AyH2]  = qVecQnFrameworkV0A->Qy(2);
  }
  delete qVectorV0A;

  // VZEROC Eventplane q-Vector
  qnListDetector = "VZEROC" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkV0C = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorV0C = new TVector2(-200.,-200.);
  if(qVecQnFrameworkV0C != NULL){
    bV0CqVector = kTRUE;
    qVectorV0C->Set(qVecQnFrameworkV0C->Qx(2),qVecQnFrameworkV0C->Qy(2));
    values[AliDielectronVarManager::kQnV0CrpH2] = TVector2::Phi_mpi_pi(qVectorV0C->Phi())/2;
    values[AliDielectronVarManager::kQnV0CxH2]  = qVecQnFrameworkV0C->Qx(2);
    values[AliDielectronVarManager::kQnV0CyH2]  = qVecQnFrameworkV0C->Qy(2);
  }
  delete qVectorV0C;

  // VZERO Eventplane q-Vector only accessible with NewDetConfig AddTask for QnFramework
  qnListDetector = "VZERO" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkV0 = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorV0 = new TVector2(-200.,-200.);
  if(qVecQnFrameworkV0 != NULL){
    bV0qVector = kTRUE;
    qVectorV0->Set(qVecQnFrameworkV0->Qx(2),qVecQnFrameworkV0->Qy(2));
    values[AliDielectronVarManager::kQnV0rpH2] = TVector2::Phi_mpi_pi(qVectorV0->Phi())/2;
    values[AliDielectronVarManager::kQnV0xH2]  = qVecQnFrameworkV0->Qx(2);
    values[AliDielectronVarManager::kQnV0yH2]  = qVecQnFrameworkV0->Qy(2);
  }
  delete qVectorV0;

  // SPD Eventplane q-Vector
  qnListDetector = "SPD" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkSPD = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorSPD = new TVector2(-200.,-200.);
  if(qVecQnFrameworkSPD != NULL){
    bSPDqVector = kTRUE;
    qVectorSPD->Set(qVecQnFrameworkSPD->Qx(2),qVecQnFrameworkSPD->Qy(2));
    values[AliDielectronVarManager::AliDielectronVarManager::kQnSPDrpH2] = TVector2::Phi_mpi_pi(qVectorSPD->Phi())/2;
    values[AliDielectronVarManager::kQnSPDxH2]  = qVecQnFrameworkSPD->Qx(2);
    values[AliDielectronVarManager::kQnSPDyH2]  = qVecQnFrameworkSPD->Qy(2);
  }
  delete qVectorSPD;

  // FMDA Eventplane q-Vector
  qnListDetector = "FMDA" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkFMDA = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorFMDA = new TVector2(-200.,-200.);
  if(qVecQnFrameworkFMDA != NULL){
    bFMDAqVector = kTRUE;
    qVectorFMDA->Set(qVecQnFrameworkFMDA->Qx(2),qVecQnFrameworkFMDA->Qy(2));
    values[AliDielectronVarManager::kQnFMDArpH2] = TVector2::Phi_mpi_pi(qVectorFMDA->Phi())/2;
    values[AliDielectronVarManager::kQnFMDAxH2]  = qVecQnFrameworkFMDA->Qx(2);
    values[AliDielectronVarManager::kQnFMDAyH2]  = qVecQnFrameworkFMDA->Qy(2);
  }
  delete qVectorFMDA;

  // FMDC Eventplane q-Vector
  qnListDetector = "FMDC" + fgQnVectorNorm;
  const AliQnCorrectionsQnVector *qVecQnFrameworkFMDC = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,qnListDetector.Data(),"latest","latest");
  TVector2 *qVectorFMDC = new TVector2(-200.,-200.);
  if(qVecQnFrameworkFMDC != NULL){
    bFMDCqVector = kTRUE;
    qVectorFMDC->Set(qVecQnFrameworkFMDC->Qx(2),qVecQnFrameworkFMDC->Qy(2));
    values[AliDielectronVarManager::kQnFMDCrpH2] = TVector2::Phi_mpi_pi(qVectorFMDC->Phi())/2;
    values[AliDielectronVarManager::kQnFMDCxH2]  = qVecQnFrameworkFMDC->Qx(2);
    values[AliDielectronVarManager::kQnFMDCyH2]  = qVecQnFrameworkFMDC->Qy(2);
  }
  delete qVectorFMDC;

  // TPC Diff
  if(bTPCqVector){
    if(bV0AqVector){
      values[AliDielectronVarManager::kQnDiffTPC_V0A] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCrpH2],values[AliDielectronVarManager::kQnV0ArpH2]);
      values[kQnCorrTPCx_V0Ax] = values[kQnTPCxH2] * values[kQnV0AxH2];
      values[kQnCorrTPCx_V0Ay] = values[kQnTPCxH2] * values[kQnV0AyH2];
      values[kQnCorrTPCy_V0Ax] = values[kQnTPCyH2] * values[kQnV0AxH2];
      values[kQnCorrTPCy_V0Ay] = values[kQnTPCyH2] * values[kQnV0AyH2];
    }
    if(bV0CqVector){
      values[AliDielectronVarManager::kQnDiffTPC_V0C] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCrpH2],values[AliDielectronVarManager::kQnV0CrpH2]);
      values[kQnCorrTPCx_V0Cx] = values[kQnTPCxH2] * values[kQnV0CxH2];
      values[kQnCorrTPCx_V0Cy] = values[kQnTPCxH2] * values[kQnV0CyH2];
      values[kQnCorrTPCy_V0Cx] = values[kQnTPCyH2] * values[kQnV0CxH2];
      values[kQnCorrTPCy_V0Cy] = values[kQnTPCyH2] * values[kQnV0CyH2];
    }
    if(bSPDqVector){
      values[AliDielectronVarManager::kQnDiffTPC_SPD] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCrpH2],values[AliDielectronVarManager::kQnSPDrpH2]);
      values[kQnCorrTPCx_SPDx] = values[kQnTPCxH2] * values[kQnSPDxH2];
      values[kQnCorrTPCx_SPDy] = values[kQnTPCxH2] * values[kQnSPDyH2];
      values[kQnCorrTPCy_SPDx] = values[kQnTPCyH2] * values[kQnSPDxH2];
      values[kQnCorrTPCy_SPDy] = values[kQnTPCyH2] * values[kQnSPDyH2];
    }
    if(bFMDAqVector){
      values[AliDielectronVarManager::kQnDiffTPC_FMDA] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCrpH2],values[AliDielectronVarManager::kQnFMDArpH2]);
      values[kQnCorrTPCx_FMDAx] = values[kQnTPCxH2] * values[kQnFMDAxH2];
      values[kQnCorrTPCx_FMDAy] = values[kQnTPCxH2] * values[kQnFMDAyH2];
      values[kQnCorrTPCy_FMDAx] = values[kQnTPCyH2] * values[kQnFMDAxH2];
      values[kQnCorrTPCy_FMDAy] = values[kQnTPCyH2] * values[kQnFMDAyH2];
    }
    if(bFMDCqVector){
      values[AliDielectronVarManager::kQnDiffTPC_FMDC] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCrpH2],values[AliDielectronVarManager::kQnFMDCrpH2]);
      values[kQnCorrTPCx_FMDCx] = values[kQnTPCxH2] * values[kQnFMDCxH2];
      values[kQnCorrTPCx_FMDCy] = values[kQnTPCxH2] * values[kQnFMDCyH2];
      values[kQnCorrTPCy_FMDCx] = values[kQnTPCyH2] * values[kQnFMDCxH2];
      values[kQnCorrTPCy_FMDCy] = values[kQnTPCyH2] * values[kQnFMDCyH2];
    }
  }

  // TPC A-Side diff
  if(bTPCaSideqVector){
    if(bTPCcSideqVector){
      values[AliDielectronVarManager::kQnDiffTPCa_TPCc] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCaSiderpH2],values[AliDielectronVarManager::kQnTPCcSiderpH2]);
    }
    if(bV0qVector){
      values[AliDielectronVarManager::kQnDiffTPCa_V0] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCaSiderpH2],values[AliDielectronVarManager::kQnV0rpH2]);
    }
    if(bV0AqVector){
      values[AliDielectronVarManager::kQnDiffTPCa_V0A] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCaSiderpH2],values[AliDielectronVarManager::kQnV0ArpH2]);
    }
    if(bV0CqVector){
      values[AliDielectronVarManager::kQnDiffTPCa_V0C] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCaSiderpH2],values[AliDielectronVarManager::kQnV0CrpH2]);
    }
  }

  // TPC C-Side diff
  if(bTPCcSideqVector){
    if(bV0qVector){
      values[AliDielectronVarManager::kQnDiffTPCc_V0] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCcSiderpH2],values[AliDielectronVarManager::kQnV0rpH2]);
    }
    if(bV0AqVector){
      values[AliDielectronVarManager::kQnDiffTPCc_V0A] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCcSiderpH2],values[AliDielectronVarManager::kQnV0ArpH2]);
    }
    if(bV0CqVector){
      values[AliDielectronVarManager::kQnDiffTPCc_V0C] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnTPCcSiderpH2],values[AliDielectronVarManager::kQnV0CrpH2]);
    }
  }

  // V0A Diff
  if(bV0AqVector){
    if(bV0CqVector){
      values[AliDielectronVarManager::kQnDiffV0A_V0C] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0ArpH2],values[AliDielectronVarManager::kQnV0CrpH2]);
      values[kQnCorrV0Ax_V0Cx] = values[kQnV0AxH2] * values[kQnV0CxH2];
      values[kQnCorrV0Ax_V0Cy] = values[kQnV0AxH2] * values[kQnV0CyH2];
      values[kQnCorrV0Ay_V0Cx] = values[kQnV0AyH2] * values[kQnV0CxH2];
      values[kQnCorrV0Ay_V0Cy] = values[kQnV0AyH2] * values[kQnV0CyH2];
    }
    if(bSPDqVector){
      values[AliDielectronVarManager::kQnDiffV0A_SPD] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0ArpH2],values[AliDielectronVarManager::kQnSPDrpH2]);
      values[kQnCorrV0Ax_SPDx] = values[kQnV0AxH2] * values[kQnSPDxH2];
      values[kQnCorrV0Ax_SPDy] = values[kQnV0AxH2] * values[kQnSPDyH2];
      values[kQnCorrV0Ay_SPDx] = values[kQnV0AyH2] * values[kQnSPDxH2];
      values[kQnCorrV0Ay_SPDy] = values[kQnV0AyH2] * values[kQnSPDyH2];
    }
    if(bFMDAqVector){
      values[AliDielectronVarManager::kQnDiffV0A_FMDA] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0ArpH2],values[AliDielectronVarManager::kQnFMDArpH2]);
      values[kQnCorrV0Ax_FMDAx] = values[kQnV0AxH2] * values[kQnFMDAxH2];
      values[kQnCorrV0Ax_FMDAy] = values[kQnV0AxH2] * values[kQnFMDAyH2];
      values[kQnCorrV0Ay_FMDAx] = values[kQnV0AyH2] * values[kQnFMDAxH2];
      values[kQnCorrV0Ay_FMDAy] = values[kQnV0AyH2] * values[kQnFMDAyH2];
    }
    if(bFMDCqVector){
      values[AliDielectronVarManager::kQnDiffV0A_FMDC] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0ArpH2],values[AliDielectronVarManager::kQnFMDCrpH2]);
      values[kQnCorrV0Ax_FMDCx] = values[kQnV0AxH2] * values[kQnFMDCxH2];
      values[kQnCorrV0Ax_FMDCy] = values[kQnV0AxH2] * values[kQnFMDCyH2];
      values[kQnCorrV0Ay_FMDCx] = values[kQnV0AyH2] * values[kQnFMDCxH2];
      values[kQnCorrV0Ay_FMDCy] = values[kQnV0AyH2] * values[kQnFMDCyH2];
    }
  }
  // V0C Diff
  if(bV0CqVector){
    if(bSPDqVector){
      values[AliDielectronVarManager::kQnDiffV0C_SPD] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0CrpH2],values[AliDielectronVarManager::kQnSPDrpH2]);
      values[kQnCorrV0Cx_SPDx] = values[kQnV0CxH2] * values[kQnSPDxH2];
      values[kQnCorrV0Cx_SPDy] = values[kQnV0CxH2] * values[kQnSPDyH2];
      values[kQnCorrV0Cy_SPDx] = values[kQnV0CyH2] * values[kQnSPDxH2];
      values[kQnCorrV0Cy_SPDy] = values[kQnV0CyH2] * values[kQnSPDyH2];
    }
    if(bFMDAqVector){
      values[AliDielectronVarManager::kQnDiffV0C_FMDA] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0CrpH2],values[AliDielectronVarManager::kQnFMDArpH2]);
      values[kQnCorrV0Cx_FMDAx] = values[kQnV0CxH2] * values[kQnFMDAxH2];
      values[kQnCorrV0Cx_FMDAy] = values[kQnV0CxH2] * values[kQnFMDAyH2];
      values[kQnCorrV0Cy_FMDAx] = values[kQnV0CyH2] * values[kQnFMDAxH2];
      values[kQnCorrV0Cy_FMDAy] = values[kQnV0CyH2] * values[kQnFMDAyH2];
    }
    if(bFMDCqVector){
      values[AliDielectronVarManager::kQnDiffV0C_FMDC] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnV0CrpH2],values[AliDielectronVarManager::kQnFMDCrpH2]);
      values[kQnCorrV0Cx_FMDCx] = values[kQnV0CxH2] * values[kQnFMDCxH2];
      values[kQnCorrV0Cx_FMDCy] = values[kQnV0CxH2] * values[kQnFMDCyH2];
      values[kQnCorrV0Cy_FMDCx] = values[kQnV0CyH2] * values[kQnFMDCxH2];
      values[kQnCorrV0Cy_FMDCy] = values[kQnV0CyH2] * values[kQnFMDCyH2];
    }
  }

  // SPD Diff
  if(bSPDqVector){
    if(bFMDAqVector){
      values[AliDielectronVarManager::kQnDiffSPD_FMDA] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnSPDrpH2],values[AliDielectronVarManager::kQnFMDArpH2]);
      values[kQnCorrSPDx_FMDAx] = values[kQnSPDxH2] * values[kQnFMDAxH2];
      values[kQnCorrSPDx_FMDAy] = values[kQnSPDxH2] * values[kQnFMDAyH2];
      values[kQnCorrSPDy_FMDAx] = values[kQnSPDyH2] * values[kQnFMDAxH2];
      values[kQnCorrSPDy_FMDAy] = values[kQnSPDyH2] * values[kQnFMDAyH2];
    }
    if(bFMDCqVector){
      values[AliDielectronVarManager::kQnDiffSPD_FMDC] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnSPDrpH2],values[AliDielectronVarManager::kQnFMDCrpH2]);
      values[kQnCorrSPDx_FMDCx] = values[kQnSPDxH2] * values[kQnFMDCxH2];
      values[kQnCorrSPDx_FMDCy] = values[kQnSPDxH2] * values[kQnFMDCyH2];
      values[kQnCorrSPDy_FMDCx] = values[kQnSPDyH2] * values[kQnFMDCxH2];
      values[kQnCorrSPDy_FMDCy] = values[kQnSPDyH2] * values[kQnFMDCyH2];
    }
  }

  // FMDA Diff
  if(bFMDAqVector && bFMDCqVector){
    values[AliDielectronVarManager::kQnDiffFMDA_FMDC] = AliDielectronVarManager::CalculateEPDiff(values[AliDielectronVarManager::kQnFMDArpH2],values[AliDielectronVarManager::kQnFMDCrpH2]);
    values[kQnCorrFMDAx_FMDCx] = values[kQnFMDAxH2] * values[kQnFMDCxH2];
    values[kQnCorrFMDAx_FMDCy] = values[kQnFMDAxH2] * values[kQnFMDCyH2];
    values[kQnCorrFMDAy_FMDCx] = values[kQnFMDAyH2] * values[kQnFMDCxH2];
    values[kQnCorrFMDAy_FMDCy] = values[kQnFMDAyH2] * values[kQnFMDCyH2];
  }
}

//________________________________________________________________
inline Double_t AliDielectronVarManager::CalculateEPDiff(Double_t detArp, Double_t detBrp){
  Double_t diff = TMath::Cos( 2.*(detArp - detBrp));
  return diff;
}

#endif
