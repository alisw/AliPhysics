//
// Author: Ionut-Cristian Arsene, 2015/04/05
// email: iarsene@cern.ch
//
// Variable manager
//

#ifndef ALIREDUCEDVARMANAGER_H
#define ALIREDUCEDVARMANAGER_H

#include <TObject.h>
#include <TString.h>
#include <TChain.h>
#include <TH2F.h>
#include <TProfile2D.h>

class AliReducedBaseEvent;
class AliReducedEventInfo;
class AliReducedEventPlaneInfo;
class AliReducedBaseTrack;
class AliReducedTrackInfo;
class AliReducedPairInfo;
class AliReducedCaloClusterInfo;

//_____________________________________________________________________
class AliReducedVarManager : public TObject {

 public:
   
  enum ParticleId {
    kUnknown = -1,
    kElectron = 0,
    kPion,
    kKaon,
    kProton,
    kKaonZero,
    kPhiMeson,
    kNSpecies
  };

     
  enum Detectors {
    kTPC = 0,
    kVZERO,
    kZDC,
    kTZERO,
    kFMD
  };

  // offline triggers as defined in AliVEvent.h
  // NOTE: Check consistency with updates in aliroot!!!
  enum EOfflineTriggerTypes { 
    kMB                = BIT(0), // Minimum bias trigger, i.e. interaction trigger, offline SPD or V0 selection
    kINT1              = BIT( 0), // V0A | V0C | SPD minimum bias trigger
    kINT7              = BIT(1), // V0AND trigger, offline V0 selection
    kMUON              = BIT(2), // Muon trigger, offline SPD or V0 selection
    kHighMult          = BIT(3), // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
    kHighMultSPD    = BIT(3), // offline SPD high multiplicity trigger
    kEMC1              = BIT(4), // EMCAL trigger
    kCINT5             = BIT(5), // Minimum bias trigger without SPD. i.e. interaction trigger, offline V0 selection
    kINT5              = BIT( 5), // V0OR minimum bias trigger
    kCMUS5             = BIT(6), // Muon trigger, offline V0 selection
    kMUSPB             = BIT(6), // idem for PbPb
    kINT7inMUON     = BIT( 6), // INT7 in MUON or MUFAST cluster
    kMuonSingleHighPt7 = BIT( 7), // Single muon high-pt, INT7 suite
    kMUSH7             = BIT(7), // Muon trigger: high pt, single muon, offline V0 selection, CINT7 suite
    kMUSHPB            = BIT(7), // idem for PbPb
    kMuonLikeLowPt7    = BIT( 8), // Like-sign dimuon low-pt, INT7 suite
    kMUL7              = BIT(8), // Muon trigger: like sign dimuon, offline V0 selection, CINT7 suite
    kMuonLikePB        = BIT(8), // idem for PbPb
    kMuonUnlikeLowPt7  = BIT( 9), // Unlike-sign dimuon low-pt, INT7 suite
    kMUU7              = BIT(9), // Muon trigger, unlike sign dimuon, offline V0 selection, CINT7 suite
    kMuonUnlikePB      = BIT(9), // idem for PbPb
    kEMC7              = BIT(10), // EMCAL trigger, INT7 suite
    kEMC8              = BIT(10), // EMCAL/DCAL L0 trigger, INT8 suite
    kMUS7              = BIT(11), // Muon trigger: low pt, single muon, offline V0 selection, CINT7 suite
    kMuonSingleLowPt7  = BIT(11), // Single muon low-pt, INT7 suite 
    kPHI1              = BIT(12), // PHOS trigger, CINT1 suite
    kPHI7              = BIT(13), // PHOS trigger, CINT7 suite
    kPHI8              = BIT(13), // PHOS trigger, INT8 suite
    kPHOSPb            = BIT(13), // idem for PbPb
    kEMCEJE            = BIT(14), // EMCAL jet patch trigger
    kEMCEGA            = BIT(15), // EMCAL gamma trigger
    kCentral           = BIT(16), // PbPb central collision trigger
    kHighMultV0    = BIT(16), // offline V0 high multiplicity trigger
    kSemiCentral       = BIT(17), // PbPb semicentral collision trigger
    kDG                = BIT(18), // Double gap diffractive
    kDG5               = BIT(18), // Double gap diffractive
    kZED               = BIT(19), // ZDC electromagnetic dissociation
    kSPI7              = BIT(20), // Power interaction trigger
    kSPI               = BIT(20), // Power interaction trigger
    kINT8              = BIT(21), // CINT8 trigger: 0TVX (T0 vertex) triger
    kMuonSingleLowPt8  = BIT(22), // Muon trigger : single muon, low pt, T0 selection, CINT8 suite
    kMuonSingleHighPt8 = BIT(23), // Muon trigger : single muon, high pt, T0 selection, CINT8 suite
    kMuonLikeLowPt8    = BIT(24), // Muon trigger : like sign muon, low pt, T0 selection, CINT8 suite
    kMuonUnlikeLowPt8  = BIT(25), // Muon trigger : unlike sign muon, low pt, T0 selection, CINT8 suite
    kMuonUnlikeLowPt0  = BIT(26), // Unlike-sign dimuon low-pt, no additional L0 requirement
    kINT6              = BIT(26),
    kUserDefined       = BIT(27), // Set when custom trigger classes are set in AliPhysicsSelection, offline SPD or V0 selection
    kTRD               = BIT(28), // TRD trigger  
    // Bits 28 and above are reserved for FLAGS
    kFastOnly          = BIT(30), // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB
    kAny               = 0xffffffff, // to accept any trigger
    kAnyINT            = kMB | kINT7 | kCINT5, // to accept any interaction (aka minimum bias) trigger
    kNTriggers         = 47
  };
  
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

  
  static const Float_t fgkParticleMass[kNSpecies];
  
  enum Variables {
    kNothing = -1,
    // Run wise variables (LHC and ALICE GRP information)
    kTotalLuminosity = 0,   // total run delivered lumi
    kBeamIntensity0,    // beam 0 intensity
    kBeamIntensity1,    // beam 1 intensity
    kLHCFillNumber,     // LHC fill number
    kDipolePolarity,    // dipole magnet polarity
    kL3Polarity,        // L3 magnet polarity
    kRunTimeStart,   // run start time
    kRunTimeEnd,     // run end time
    kNRunWiseVariables,
    
    // Event wise variables
    kEventTag = kNRunWiseVariables,      // event tag
    kEventNumberInFile, // event number in file
    kL0TriggerInput,    // L0 trigger input
    kL1TriggerInput,    // L1 trigger input
    kL2TriggerInput,    // L2 trigger input
    kRunNo,             // run number         
    kRunID,             // variable for easy filling of histograms vs. run number, without empty bins
    kBeamEnergy,        // LHC beam energy
    kDetectorMask,      // detector mask
    kNumberOfDetectors, // number of active detectors
    kBC,                // bunch crossing     
    kTimeStamp,         // time stamp of the event
    kTimeRelativeSOR,   // time relative to the start of run, in minutes
    kTimeRelativeSORfraction,   // time relative to the start of runs, expressed as fraction of the whole run duration 
    kEventType,         // event type
    kTriggerMask,       // trigger mask       
    kOnlineTriggersFired,// 1 if fired, 0 if not fired, for each trigger
    kOnlineTrigger=kOnlineTriggersFired+kNTriggers,  // online trigger    
    kOnlineTriggerFired,  // online trigger fired
    kOnlineTriggerFired2,  // online trigger if fired, -1 if not fired
    kIsPhysicsSelection,    // physics selection 
    kIsSPDPileup,          // whether is SPD pileup
    kIsSPDPileup5,          // whether is SPD pileup (5 vertex contributors)
    kIsSPDPileupMultBins,  // whether is SPD pileup in multiplicity bins
    kNSPDpileups,         // number of pileup events from SPD
    kNTrackPileups,       // number of pileup events from tracks
    kIRIntClosestIntMap,  // map of closest out of bunch interactions; [0]-Int1, [1]-Int2
    kNPMDtracks=kIRIntClosestIntMap+2,   // number of PMD tracks
    kNTRDtracks,          // number of TRD tracks
    kNTRDtracklets,          // number of TRD tracklets
    kNVtxContributors,    // number of vertex contributors
    kNVtxTPCContributors,  // number of TPC vertex contributors
    kNVtxSPDContributors,  // number of SPD vertex contributors
    kVtxX,              // vtx X                      
    kVtxY,              // vtx Y                      
    kVtxZ,              // vtx Z 
    kVtxXtpc,           // vtx X from tpc
    kVtxYtpc,           // vtx Y from tpc
    kVtxZtpc,           // vtx Z from tpc
    kDeltaVtxZ,         // vtxZ - vtxZtpc
    kVtxXspd,           // vtx X from spd
    kVtxYspd,           // vtx Y from spd
    kVtxZspd,           // vtx Z from spd
    kDeltaVtxZspd,         // vtxZ - vtxZspd
    kNTracksPerTrackingStatus,  // number of tracks with a given tracking flag
    kNTracksTPCoutVsITSout=kNTracksPerTrackingStatus+kNTrackingStatus,   //  TPCout/ITSout
    kNTracksTRDoutVsITSout,                              //  TRDout/ITSout
    kNTracksTOFoutVsITSout,                              //  TOFout/ITSout
    kNTracksTRDoutVsTPCout,                              //  TRDout/TPCout
    kNTracksTOFoutVsTPCout,                              //  TOFout/TPCout
    kNTracksTOFoutVsTRDout,                              //  TOFout/TRDout
    kNTracksITSoutVsSPDtracklets,                        //  ITSout/SPDtracklets
    kNTracksTPCoutVsSPDtracklets,                        //  TPCout/SPDtracklets
    kNTracksTRDoutVsSPDtracklets,                        //  TRDout/SPDtracklets
    kNTracksTOFoutVsSPDtracklets,                        //  TOFout/SPDtracklets
    kNTracksTPCoutFromPileup,                       // number of tracks from (kNTracksPerTrackingStatus+kTPCout) minus the no-pileup expectation
    kNTracksTPCoutVsVZEROTotalMult,      // number of kTPCout tracks / VZERO multiplicity
    kCentVZERO,         // centrality from VZERO
    kCentSPD,           // centrality from SPD
    kCentSPDcorr,       // corrected centrality from SPD
    kCentTPC,           // centrality from TPC  
    kCentZDC,           // centrality from ZDC  
    kCentVZEROA,        // centrality from VZERO-A
    kCentVZEROC,        // centrality from VZERO-C
    kCentZNA,           // centrality from ZNA
    kCentQuality,       // centrality quality   
    kNV0total,          // total number of V0s in the esd      
    kNV0selected,       // number of V0s selected              
    kNpairsSelected,    // number of selected pairs per event  
    kEvAverageTPCchi2,   // average TPC chi2 for the tracks in a given event
    kNDplusToK0sPiplusSelected,       // D+           -> K0s pi+
    kNDplusToK0sKplusSelected,        // D+           -> K0s K+
    kNDplusToPhiPiplusSelected,       // D+           -> phi pi+
    kNDminusToK0sPiminusSelected,     // D-           -> K0s pi-
    kNDminusToK0sKminusSelected,      // D-           -> K0s K-
    kNDminusToPhiPiminusSelected,     // D-           -> phi pi-
    kNDzeroToKminusPiplusSelected,    // D0           -> K- pi+
    kNADzeroToKplusPiminusSelected,   // anti-D0      -> K+ pi-
    kNDsplusToK0sKplusSelected,       // Ds+          -> K0s K+
    kNDsminusToK0sKminusSelected,     // Ds-          -> K0s K-    
    kNtracksTotal,      // total number of tracks               
    kNtracksSelected,   // number of selected tracks            
    kNtracksPosAnalyzed,// number of positive tracks passing analysis cuts      
    kNtracksNegAnalyzed,// number of negative tracks passing analysis cuts      
    kNtracksPiPlusAnalyzed,     // number of pi plus selected tracks
    kNtracksPiMinusAnalyzed,    // number of pi minus selected tracks
    kNtracksKPlusAnalyzed,      // number of K plus selected tracks
    kNtracksKMinusAnalyzed,     // number of K minus selected tracks
    kNK0sAnalyzed,              // number of K0s candidates selected
    kNPhiAnalyzed,              // number of phi candidates selected
    kNtracksAnalyzed,   // number of positive+negative tracks passing analysis cuts
    kNtracksAnalyzedInPhiBins,   // number of positive+negative tracks passing analysis cuts in 36 phi bins (18 for each side of the TPC)
    kNtracksSubEvLeft=kNtracksAnalyzedInPhiBins+36,  // number of tracks in the left sub-event (negative pseudo-rapidity)            
    kNtracksSubEvRight, // number of tracks in the left sub-event (positive pseudo-rapidity)            
    kNtracksEventPlane, // number of tracks used for event plane                
    kNCaloClusters,     // number of calorimeter clusters
    kNTPCclusters,    // number of TPC clusters
    kSPDntracklets,     // SPD number of tracklets in |eta|<1.0                 
    kSPDntracklets08,     // SPD number of tracklets in |eta|<0.8
    kSPDntracklets16,     // SPD number of tracklets in |eta|<1.6
    kSPDntrackletsCorr, // SPD number of tracklets in |eta|<1.0, corrected for detector effects
    kSPDntrackletsCorrSmear, // SPD number of tracklets in |eta|<1.0, corrected for detector effects, Poisson smeared
    kSPDntrackletsOuter,     // SPD number of tracklets in outer eta region
    kSPDntrackletsOuterCorr,     // SPD number of tracklets in outer eta region, corrected for detector effects
    kSPDntrackletsOuterCorrSmear,     // SPD number of tracklets in outer eta region, corrected for detector effects, Poisson smeared
    kSPDntrackletsEta,  // SPD number of tracklets in -1.6+0.1*i < eta < -1.6+0.1*(i+1)
    kSPDFiredChips=kSPDntrackletsEta+32,   // SPD fired chips in first and second layer
    kITSnClusters=kSPDFiredChips+2,        // number of ITS clusters in each layer
    kSPDnSingleClusters=kITSnClusters+6,   // number of clusters in SPD layer 1 not mached to tracklets from layer 2
    kEventMixingId,     // Id of the event mixing category 
    // VZERO event plane related variables
    kVZEROATotalMult,   // total multiplicity of VZEROA                         
    kVZEROCTotalMult,   // total multiplicity of VZEROC                         
    kVZEROTotalMult,    // total multiplicity of VZERO                          
    kVZEROTotalMultCorr,    // total multiplicity of VZERO, corrected for detector effects
    kVZEROTotalMultCorrSmear,    // total multiplicity of VZERO, corrected for detector effects, Poisson smeared
    kVZEROAemptyChannels,  // Number of empty VZERO channels in A side          
    kVZEROCemptyChannels,  // Number of empty VZERO channels in C side          
    kVZEROChannelMult,                        // VZERO multiplicity per channel           
    kVZEROChannelEta = kVZEROChannelMult+64,  // pseudo-rapidity of a VZERO channel       
    kVZEROQvecX      = kVZEROChannelEta+64,   // Q-vector components for harmonics 1-6 and  
    kVZEROQvecY      = kVZEROQvecX+6*3,        //  6- n-harmonics; 3- A,C and A&C options   
    kVZEROQvecMag    = kVZEROQvecY+6*3,       // magnitude of the Q vector
    kVZERORP         = kVZEROQvecMag+6*3,      // VZERO reaction plane from A,C and A&C sides (harmonics 1-6)  
    kVZERORPres      = kVZERORP+6*3,     // VZERO reaction plane resolution (sqrt(n*(RPa-RPc)))        
    kVZEROXaXc       = kVZERORPres+6,           // correlations for the components of the Q vector     
    kVZEROXaYa       = kVZEROXaXc+6,                            
    kVZEROXaYc       = kVZEROXaYa+6,                            
    kVZEROYaXc       = kVZEROXaYc+6,                            
    kVZEROYaYc       = kVZEROYaXc+6,                            
    kVZEROXcYc       = kVZEROYaYc+6,                            
    kVZEROdeltaRPac  = kVZEROXcYc+6,         // Psi_VZEROA-Psi_VZEROC
    kVZEROflowV2TPC  = kVZEROdeltaRPac+6,     // vzero v2 using TPC event plane        
    kVZEROQaQcSP     = kVZEROflowV2TPC+64,     // scalar product for VZERO-A Q  times  VZERO-C Q  (just the cosine term)
    kVZEROQaQcSPsine = kVZEROQaQcSP + 6,     // sine term from the scalar product
    // TPC event plane variables
    kTPCQvecX = kVZEROQaQcSPsine+6,   // TPC Q-vector components for harmonics 1-6     
    kTPCQvecY = kTPCQvecX+6,                                                           
    kTPCRP    = kTPCQvecY+6,                // Event plane using TPC                    
    kTPCRPres = kTPCRP+6,                // Event plane resolution variables sqrt(n*(RPtpc-RPvzeroa)),sqrt(n*(RPtpc-RPvzeroc))
    // Correlations between TPC and VZERO event planes
    kRPXtpcXvzeroa    = kTPCRPres+6*2,          
    kRPXtpcXvzeroc    = kRPXtpcXvzeroa+6,       
    kRPYtpcYvzeroa    = kRPXtpcXvzeroc+6,       
    kRPYtpcYvzeroc    = kRPYtpcYvzeroa+6,       
    kRPXtpcYvzeroa    = kRPYtpcYvzeroc+6,       
    kRPXtpcYvzeroc    = kRPXtpcYvzeroa+6,       
    kRPYtpcXvzeroa    = kRPXtpcYvzeroc+6,       
    kRPYtpcXvzeroc    = kRPYtpcXvzeroa+6,       
    kRPdeltaVZEROAtpc = kRPYtpcXvzeroc+6,       
    kRPdeltaVZEROCtpc = kRPdeltaVZEROAtpc+6,    
    // TPC event plane using sub-intervals in pseudo-rapidity
    kTPCQvecXleft   = kRPdeltaVZEROCtpc+6,      
    kTPCQvecYleft   = kTPCQvecXleft+6,          
    kTPCRPleft      = kTPCQvecYleft+6,          
    kTPCQvecXright  = kTPCRPleft+6,             
    kTPCQvecYright  = kTPCQvecXright+6,         
    kTPCRPright     = kTPCQvecYright+6,         
    kTPCQvecXtotal  = kTPCRPright+6,            
    kTPCQvecYtotal  = kTPCQvecXtotal+6,         
    kTPCRPtotal     = kTPCQvecYtotal+6,         
    kTPCsubResCos   = kTPCRPtotal+6,            
    // ZDC variables
    kZDCnEnergyCh   = kTPCsubResCos+6,         // ZDCn energy in each channel
    kZDCpEnergyCh   = kZDCnEnergyCh+10,        // ZDCp energy in each channel
    // TZERO variables
    kTZEROAmplitudeCh = kZDCpEnergyCh+10,      // TZERO aplitudes in all channels
    kTZEROTOF         = kTZEROAmplitudeCh+26,  // TZERO TOF start times
    kTZEROTOFbest     = kTZEROTOF+3,           // TZERO TOF best start times
    kTZEROzVtx        = kTZEROTOFbest+3,       // TZERO event z vertex
    kTZEROstartTime,                           // TZERO event start time
    kTZEROpileup,                              // TZERO pileup flag
    kTZEROsatellite,                           // TZERO satellite flag
    // Multiplicity estimators
    kMultEstimatorOnlineV0M,
    kMultEstimatorOnlineV0A,
    kMultEstimatorOnlineV0C,
    kMultEstimatorADM,
    kMultEstimatorADA,
    kMultEstimatorADC,
    kMultEstimatorSPDClusters,
    kMultEstimatorSPDTracklets,
    kMultEstimatorRefMult05,
    kMultEstimatorRefMult08,
    kMultEstimatorPercentileOnlineV0M,
    kMultEstimatorPercentileOnlineV0A,
    kMultEstimatorPercentileOnlineV0C,
    kMultEstimatorPercentileADM,
    kMultEstimatorPercentileADA,
    kMultEstimatorPercentileADC,
    kMultEstimatorPercentileSPDClusters,
    kMultEstimatorPercentileSPDTracklets,
    kMultEstimatorPercentileRefMult05,
    kMultEstimatorPercentileRefMult08,
    kINT7Triggered,
    kHighMultV0Triggered,
    kNEventVars,                               // number of event variables  
    // Particle variables --------------------------------------
    // Common pair/track variables
    kPt=kNEventVars,
    kPtMC,
    kP,      
    kPMC,
    kPx,   
    kPxMC,
    kPy,     
    kPyMC,
    kPz,     
    kPzMC,
    kTheta,
    kThetaMC,
    kEta,
    kEtaMC,
    kPhi,     
    kPhiMC,
    kCosNPhi,   
    kSinNPhi = kCosNPhi+6,
    kPtSquared = kSinNPhi+6,
    kOneOverSqrtPt,                   // one over square root of pT
    kMass,
    kMassMC,
    kRap,
    kRapMC,
    kPdgMC,
    kCharge = kPdgMC+4,
    kVZEROFlowVn,                     // v_n using VZERO RP
    kTPCFlowVn=kVZEROFlowVn+6*3,      // v_n using TPC RP
    kVZEROFlowSine=kTPCFlowVn+6,      // sin(n*(phi-Psi)) using VZERO RP
    kTPCFlowSine=kVZEROFlowSine+6*3,  // sin(n*(phi-Psi)) using TPC RP
    kVZEROuQ = kTPCFlowSine+6,        // cosine term from the u*Q products from VZERO (harmonics 1-6; VZERO-A and VZERO-C)
    kVZEROuQsine = kVZEROuQ+6*2,      // sine terms from the u*Q products from VZERO (harmonics 1-6; VZERO-A and VZERO-C)
    kTPCuQ=kVZEROuQsine+6*2,          // cosine terms from the u*Q products from TPC (harmonics 1-6)
    kTPCuQsine=kTPCuQ+6,              // sine terms from the u*Q products from TPC (harmonics 1-6)
    // Pair-only variables
    kCandidateId=kTPCuQsine+6,
    kPairType,                  // 0 ++; 1 +-; 2 --    
    kMassV0,                    // masses for all 4 V0 assumptions (0-K0s, 1-Lambda, 2-ALambda, 3-Gamma)
    kPairChisquare=kMassV0+4,     
    kPairLxy,           
    kPairOpeningAngle,  
    kPairPointingAngle, 
    kPairThetaCS,                // cos (theta*) in Collins-Soper frame       
    kPairPhiCS,                    // phi* in Collins-Soper frame
    kPairThetaHE,                // cos (theta*) in helicity frame       
    kPairPhiHE,                    // phi* in helicity frame
    kPairQualityFlag,
    kDMA,                        // Distance of minimal approach
    kPairPhiV,                   // angle between pair plane and magnetic field
    kPairDca,                    // pair DCA: sqare root of quadratic sum of daughter DCAs
    kPairDcaXY,
    kPairDcaZ,
    kPairDcaSqrt,                // square root of pair DCA
    kPairDcaXYSqrt,
    kPairDcaZSqrt,
    kMassDcaPtCorr,             // invariant mass, corrected for DCA and pT effects
    kOpAngDcaPtCorr,            // opening angle, corrected for DCA and pT effects
    kPairLegITSchi2,              // the ITS chi2 for the pair legs, used in correlations between pair legs
    kPairLegTPCchi2=kPairLegITSchi2+2,              // the TPC chi2 for the pair legs, used in correlations between pair legs
    // Track-only variables -------------------------------------
    kPtTPC=kPairLegTPCchi2+2,     
    kPhiTPC,    
    kEtaTPC,    
    kDcaXYTPC,    
    kDcaZTPC,    
    kPin,       
    kDcaXY,     
    kDcaZ,              
    kTrackLength,       // track length
    kChi2TPCConstrainedVsGlobal,
    kMassUsedForTracking,
    kITSncls,
    kNclsSFracITS,
    kITSchi2,
    kITSnclsShared,
    kITSlayerHit,       
    kITSsignal,         
    kITSnSig,
    kTPCncls=kITSnSig+4,    
    kTPCchi2,
    kTPCclusBitFired,   
    kTPCNclusBitsFired, 
    kTPCclustersPerBit, 
    kTPCcrossedRows,    
    kTPCnclsF,  
    kTPCnclsShared,
    kTPCnclsSharedRatio,
    kTPCnclsRatio,       // TPCncls / TPCnclsF          
    kTPCnclsRatio2,      // TPCncls / TPCCrossedRows
    //TODO: TPC number of crossed rows over findable clusters has at the moment 2 variables assigned: kTPCcrossedRowsOverFindableClusters and kTPCnclsRatio3
    kTPCcrossedRowsOverFindableClusters,
    kTPCnclsRatio3,      // TPCCrossedRows/TPCnclsF
    kTPCsignal,         
    kTPCsignalN,
    kTPCnSig,  
    kTPCnSigCorrected=kTPCnSig+4,
    kTOFbeta=kTPCnSigCorrected+4,
    kTOFtime,
    kTOFdx,
    kTOFdz,
    kTOFmismatchProbability,
    kTOFchi2,
    kTOFdeltaBC,
    kTOFnSig,                   
    kBayes=kTOFnSig+4,
    kTRDntracklets=kBayes+4,  
    kTRDntrackletsPID,          
    kTRDpidProbabilitiesLQ1D,   
    kTRDpidProbabilitiesLQ2D=kTRDpidProbabilitiesLQ1D+2,
    kEMCALmatchedEnergy=kTRDpidProbabilitiesLQ2D+2,         
    kEMCALmatchedClusterId,
    kEMCALmatchedEOverP,        
    // Calorimeter cluster variables --------------------------------------
    kEMCALclusterEnergy,        
    kEMCALclusterDx,            
    kEMCALclusterDz,            
    kEMCALdetector,         // 0 - EMCAL; 1 - PHOS     
    kEMCALm20,
    kEMCALm02,
    kEMCALdispersion,
    // Track flags -----------------------------------------------------
    kTrackingFlag,
    kTrackQualityFlag,
    // Correlation variables ----------------------------------------------
    kDeltaPhi,          
    kDeltaTheta,        
    kDeltaEta,          
    kTrackingFlags,
    kTrackingStatus=kTrackingFlags+kNTrackingFlags,
    kNVars=kTrackingStatus+kNTrackingStatus,     
  };
  
  static TString fgVariableNames[kNVars];         // variable names
  static TString fgVariableUnits[kNVars];         // variable units  
  static const Char_t* fgkTrackingStatusNames[kNTrackingStatus];  // tracking flags name
  
  static const Char_t* fgkOfflineTriggerNames[64];
  
  static const Double_t fgkVZEROChannelRadii[64]; // radii of VZERO channels centers (in cm)
  static const Double_t fgkVZEROAz;      // z-position for VZERO-A
  static const Double_t fgkVZEROCz;    // z-position for VZERO-C
  static const Double_t fgkVZEROminMult;   // minimum VZERO channel multiplicity
  static const Float_t fgkTPCQvecRapGap;    // symmetric interval in the middle of the TPC excluded from EP calculation
    
  AliReducedVarManager();
  AliReducedVarManager(const Char_t* name);
  virtual ~AliReducedVarManager();
  
  static void SetBeamMomentum(Float_t beamMom) {fgBeamMomentum = beamMom;}
  static Float_t GetBeamMomentum() {return fgBeamMomentum;}
  
  static void SetEvent(AliReducedBaseEvent* const ev) {fgEvent = ev;};
  static void SetEventPlane(AliReducedEventPlaneInfo* const ev) {fgEventPlane = ev;};
  static void SetUseVariable(Variables var) {fgUsedVars[var] = kTRUE; SetVariableDependencies();}
  static void SetUseVars(Bool_t* usedVars) {
    for(Int_t i=0;i<kNVars;++i) {
      if(usedVars[i]) fgUsedVars[i]=kTRUE;    // overwrite only the variables that are being used since there are more channels to modify the used variables array, independently
    }
    SetVariableDependencies();
  }
  static Bool_t GetUsedVar(Variables var) {return fgUsedVars[var];}
  
  static void FillEventInfo(Float_t* values);
  static void FillEventInfo(AliReducedBaseEvent* event, Float_t* values, AliReducedEventPlaneInfo* eventPlane=0x0);
  static void FillEventOnlineTriggers(AliReducedEventInfo* event, Float_t* values);
  static void FillEventOnlineTrigger(UShort_t triggerBit, Float_t* values);
  static void FillEventTagInput(AliReducedBaseEvent* event, Int_t input, Float_t* values);
  static void FillL0TriggerInputs(AliReducedEventInfo* event, Int_t input, Float_t* values);
  static void FillL1TriggerInputs(AliReducedEventInfo* event, Int_t input, Float_t* values);
  static void FillL2TriggerInputs(AliReducedEventInfo* event, Int_t input, Float_t* values);
  static void FillTrackingFlag(AliReducedTrackInfo* track, UInt_t flag, Float_t* values);
  static void FillTrackQualityFlag(AliReducedBaseTrack* track, UShort_t flag, Float_t* values);
  static void FillPairQualityFlag(AliReducedPairInfo* p, UShort_t flag, Float_t* values);
  static void FillTrackInfo(AliReducedBaseTrack* p, Float_t* values);
  static void FillITSlayerFlag(AliReducedTrackInfo* track, Int_t layer, Float_t* values);
  static void FillTPCclusterBitFlag(AliReducedTrackInfo* track, Int_t bit, Float_t* values);
  static void FillPairInfo(AliReducedPairInfo* p, Float_t* values);
  static void FillPairInfo(AliReducedBaseTrack* t1, AliReducedBaseTrack* t2, Int_t type, Float_t* values);
  static void FillPairInfo(AliReducedPairInfo* leg1, AliReducedBaseTrack* leg2, Int_t type, Float_t* values);
  static void FillPairInfoME(AliReducedBaseTrack* t1, AliReducedBaseTrack* t2, Int_t type, Float_t* values);
  static void FillCorrelationInfo(AliReducedPairInfo* p, AliReducedBaseTrack* t, Float_t* values);
  static void FillCaloClusterInfo(AliReducedCaloClusterInfo* cl, Float_t* values);
  static void FillTrackingStatus(AliReducedTrackInfo* p, Float_t* values);
  static void FillTrackingFlags(AliReducedTrackInfo* p, Float_t* values);
  static void FillMCTruthInfo(AliReducedTrackInfo* p, Float_t* values, AliReducedTrackInfo* leg1 = 0x0, AliReducedTrackInfo* leg2 = 0x0);
  
  static void PrintTrackFlags(AliReducedTrackInfo* track);
  static void PrintBits(ULong_t mask, Int_t maxBit=64);
  static void PrintBits(UInt_t mask, Int_t maxBit=32);
  static void SetDefaultVarNames();
  //static TChain* GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
  //                        TChain* friendChain=0x0, const Char_t* friendChainFile=0x0);
  static TChain* GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                       TChain* friendChain=0x0, const Char_t* friendDir=0x0, const Char_t* friendFileName=0x0);

  static TString GetVarName(Int_t var) {return fgVariableNames[var];} 
  static TString GetVarUnit(Int_t var) {return fgVariableUnits[var];} 

  static void SetTPCelectronCorrectionMaps(TH2F* centroidMap, TH2F* widthMap, Variables xVarDep, Variables yVarDep);
  static void SetLHCDataInfo(TH1F* totalLumi, TH1F* totalInt0, TH1F* totalInt1, TH1I* fillNumber);
  static void SetGRPDataInfo(TH1I* dipolePolarity, TH1I* l3Polarity, TH1I* timeStart, TH1I* timeStop);
  static void SetRunNumbers( TString runNumbers );
  static void SetMultiplicityProfile( TH1* profile, Int_t estimator = 0 );
  static void SetVZEROCalibrationPath(const Char_t* path);
  static void SetCalibrateVZEROqVector(Bool_t option);
  static void SetRecenterVZEROqVector(Bool_t option);
  
 private:
  static Int_t     fgCurrentRunNumber;               // current run number
  static Float_t fgBeamMomentum;                  // beam energy (needed when calculating polarization angles) 
  static AliReducedBaseEvent* fgEvent;            // pointer to the current event
  static AliReducedEventPlaneInfo* fgEventPlane;  // pointer to the current event plane
  static Bool_t fgUsedVars[kNVars];              // array of flags toggled when the corresponding variable is required (e.g., in the histogram manager, in cuts, mixing handler, etc.) 
                                                 //   when a variable is used
  static void SetVariableDependencies();       // toggle those variables on which other used variables might depend 
  

  static Double_t DeltaPhi(Double_t phi1, Double_t phi2);  
  static void GetThetaPhiCM(AliReducedBaseTrack* leg1, AliReducedBaseTrack* leg2,
                            Float_t &thetaHE, Float_t &phiHE, 
			    Float_t &thetaCS, Float_t &phiCS,
			    Float_t leg1Mass=fgkParticleMass[kElectron], Float_t leg2Mass=fgkParticleMass[kElectron]);
  static void GetLegMassAssumption(Int_t id, Float_t& m1, Float_t& m2);

  static TH2F* fgTPCelectronCentroidMap;    // TPC electron centroid 2D map
  static TH2F* fgTPCelectronWidthMap;       // TPC electron width 2D map
  static Variables fgVarDependencyX;        // varX in the 2-D electron correction maps
  static Variables fgVarDependencyY;        // varY in the 2-D electron correction maps
  
  static TH1F* fgRunTotalLuminosity;      // total luminosity, GRP/GRP/LHCData::GetLumiAliceSBDelivered()
  static TH1F* fgRunTotalIntensity0;        // total intensity beam 1, GRP/GRP/LHCData::GetTotalIntensity(0)
  static TH1F* fgRunTotalIntensity1;         // total intensity beam 2, GRP/GRP/LHCData::GetTotalIntensity(1)
  static TH1I* fgRunLHCFillNumber;         // LHC fill number, GRP/GRP/LHCData::GetFillNumber()
  static TH1I* fgRunDipolePolarity;          // dipole magnet polarity, GRP/GRP/Data::GetDipolePolarity()
  static TH1I* fgRunL3Polarity;                // L3 magnet polarity, GRP/GRP/Data::GetL3Polarity()
  static TH1I* fgRunTimeStart;                // run start time, GRP/GRP/Data::GetTimeStart()
  static TH1I* fgRunTimeEnd;                  // run stop time, GRP/GRP/Data::GetTimeEnd()
  static std::vector<Int_t> fgRunNumbers;     // vector with run numbers (for histograms vs. run number)
  static Int_t fgRunID;                       // run ID
  static TH1* fgAvgMultVsVertex[3];           // average multiplicity vs. z-vertex position
  static Double_t fgRefMult[3];                  // reference multiplicity for z-vertex correction
  static TString fgVZEROCalibrationPath;       // path to the VZERO calibration histograms
  static TProfile2D* fgAvgVZEROChannelMult[64];       // average multiplicity in VZERO channels vs (vtxZ,centSPD)
  static TProfile2D* fgVZEROqVecRecentering[4];       // (vtxZ,centSPD) maps of the VZERO A and C recentering Qvector offsets
  static Bool_t fgOptionCalibrateVZEROqVec;
  static Bool_t fgOptionRecenterVZEROqVec;
  
  AliReducedVarManager(AliReducedVarManager const&);
  AliReducedVarManager& operator=(AliReducedVarManager const&);  
  
  ClassDef(AliReducedVarManager,2);
};

#endif
