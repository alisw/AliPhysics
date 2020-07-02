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
#include <TH3F.h>
#include <TProfile2D.h>
#include <TGraphErrors.h>
#include <THn.h>

#include <AliReducedPairInfo.h>

class AliReducedBaseEvent;
class AliReducedEventInfo;
class AliReducedEventPlaneInfo;
class AliReducedBaseTrack;
class AliReducedTrackInfo;
class AliReducedCaloClusterInfo;
class AliReducedCaloClusterTrackMatcher;
class AliKFParticle;

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
    kINT1              = BIT(0), // V0A | V0C | SPD minimum bias trigger
    kINT7              = BIT(1), // V0AND trigger, offline V0 selection
    kMUON              = BIT(2), // Muon trigger, offline SPD or V0 selection
    kHighMult          = BIT(3), // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
    kHighMultSPD       = BIT(3), // offline SPD high multiplicity trigger
    kEMC1              = BIT(4), // EMCAL trigger
    kCINT5             = BIT(5), // Minimum bias trigger without SPD. i.e. interaction trigger, offline V0 selection
    kINT5              = BIT(5), // V0OR minimum bias trigger
    kCMUS5             = BIT(6), // Muon trigger, offline V0 selection
    kMUSPB             = BIT(6), // idem for PbPb
    kINT7inMUON        = BIT(6), // INT7 in MUON or MUFAST cluster
    kMuonSingleHighPt7 = BIT(7), // Single muon high-pt, INT7 suite
    kMUSH7             = BIT(7), // Muon trigger: high pt, single muon, offline V0 selection, CINT7 suite
    kMUSHPB            = BIT(7), // idem for PbPb
    kMuonLikeLowPt7    = BIT(8), // Like-sign dimuon low-pt, INT7 suite
    kMUL7              = BIT(8), // Muon trigger: like sign dimuon, offline V0 selection, CINT7 suite
    kMuonLikePB        = BIT(8), // idem for PbPb
    kMuonUnlikeLowPt7  = BIT(9), // Unlike-sign dimuon low-pt, INT7 suite
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
    kHighMultV0        = BIT(16), // offline V0 high multiplicity trigger
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
    kMuonCalo          = BIT(29), // Muon-calo triggers
    kCaloOnly          = BIT(29), // MB, EMCAL and PHOS triggers in CALO or CALOFAST cluster
    // Bits 28 and above are reserved for FLAGS
    kFastOnly          = BIT(30), // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB
    kAny               = 0xffffffff, // to accept any trigger
    kAnyINT            = kMB | kINT7 | kCINT5 | kINT8 | kSPI7, // to accept any interaction (aka minimum bias) trigger
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

  enum Corrections {
    kVertexCorrectionGlobal=0,
    kVertexCorrectionRunwise,
    kVertexCorrectionGlobalGainLoss,
    kVertexCorrectionRunwiseGainLoss,
    kVertexCorrection2D,
    kGainLossCorrection,
    kNCorrections
  };

  enum ReferenceMultiplicities {
   kMaximumMultiplicity=0,
   kMinimumMultiplicity,
   kMeanMultiplicity,
   kNReferenceMultiplicities
  };

  enum SmearingMethods {
   kNoSmearing=0,
   kPoissonSmearing,
   kNSmearingMethods
  };

  
  static const Float_t fgkParticleMass[kNSpecies];
  static const Float_t fgkPairMass[AliReducedPairInfo::kNMaxCandidateTypes];
  
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
    kL0TriggerInput2,    // L0 trigger input, used for correlations between inputs
    kL1TriggerInput,    // L1 trigger input
    kL1TriggerInput2,    // L1 trigger input, used for correlations between inputs
    kL2TriggerInput,    // L2 trigger input
    kL2TriggerInput2,    // L2 trigger input, used for correlations between inputs
    kRunNo,             // run number         
    kRunID,             // variable for easy filling of histograms vs. run number, without empty bins
    kBeamEnergy,        // LHC beam energy
    kInstLumi,           // instantaneous interaction rate
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
    kIsPileupMV,            // pileup from multi vertexer
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
    kTPCpileupZAC,      // TPC pileup event Z from A&C sides  
    kTPCpileupZA,       // TPC pileup event Z from A side
    kTPCpileupZC,       // TPC pileup event Z from C side
    kTPCpileupContributorsAC,    // TPC pileup event contributors from A&C sides
    kTPCpileupContributorsA,     // TPC pileup event contributors from A side
    kTPCpileupContributorsC,     // TPC pileup event contributors from C side
    kTPCpileupZAC2,      // TPC pileup event Z with larger DCAz selection from A&C sides  
    kTPCpileupZA2,       // TPC pileup event Z with larger DCAz selection from A side
    kTPCpileupZC2,       // TPC pileup event Z with larger DCAz selection from C side
    kTPCpileupContributorsAC2,    // TPC pileup event contributors with larger DCAz selection from A&C sides
    kTPCpileupContributorsA2,     // TPC pileup event contributors with larger DCAz selection from A side
    kTPCpileupContributorsC2,     // TPC pileup event contributors with larger DCAz selection from C side
    kNTracksPerTrackingStatus,  // number of tracks with a given tracking flag
    kNTracksTPCoutBeforeClean=kNTracksPerTrackingStatus+kNTrackingStatus,      // TPCout tracks before ESD cleaning
    kNTracksTPCoutVsITSout,                              //  TPCout/ITSout
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
    kNTPCclustersFromPileup,            // number of TPC clusters minus the expected TPC clusters if no pileup is present
    kNTPCclustersFromPileupRelative,    // number of TPC clusters minus the expected TPC clusters w/o pileup relative to the TPC clusters w/o pileup 
    kMultiplicity,
    kSPDntracklets = kMultiplicity,
    kSPDntracklets08,
    kSPDntracklets16,
    kSPDntrackletsOuterEta,
    kSPDntrackletsEtaBin,
    kSPDnTracklets10EtaVtxCorr = kSPDntrackletsEtaBin + 32,
    kVZEROTotalMult,
    kVZEROATotalMult,
    kVZEROCTotalMult,
    kVZEROTotalMultFromChannels,
    kVZEROATotalMultFromChannels,
    kVZEROCTotalMultFromChannels,
    kVZEROTPCoutDiff,
    kVZEROACTotalMult,
    kCorrectedMultiplicity,
    kNMultiplicityEstimators = kCorrectedMultiplicity - kMultiplicity,
    kSPDFiredChips = kCorrectedMultiplicity + kNMultiplicityEstimators * ( 1 + kNCorrections * kNReferenceMultiplicities * kNSmearingMethods), // SPD fired chips in first and second layer
    kITSnClusters=kSPDFiredChips+2,        // number of ITS clusters in each layer
    kSPDnSingleClusters=kITSnClusters+6,   // number of clusters in SPD layer 1 not mached to tracklets from layer 2
    kSDDandSSDclusters,                    // number of clusters in the SDD and SSD layers
    kEventMixingId,     // Id of the event mixing category 
    // VZERO event plane related variables
    kVZEROCurrentChannel,         // current VZERO channel
    kVZEROCurrentChannelMult,     // current VZERO channel multiplicity
    kVZEROCurrentChannelMultCalib,  // current VZERO channel calibrated multiplicity
    kVZEROAemptyChannels,  // Number of empty VZERO channels in A side          
    kVZEROCemptyChannels,  // Number of empty VZERO channels in C side          
    kVZEROChannelMult,                        // VZERO multiplicity per channel 
    kVZEROChannelMultCalib=kVZEROChannelMult+64,                   // VZERO multiplicity per channel calibrated
    kVZEROChannelEta = kVZEROChannelMultCalib+64,  // pseudo-rapidity of a VZERO channel       
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
    kVZEROARPres=kTPCRPres+6,           //event plane resolution using V0A as reference detector 
    kVZEROCRPres=kVZEROARPres+6,       //event plane resolution using V0C as reference detector
    kVZEROTPCRPres=kVZEROCRPres+6,    //event plane resolution using tpc as reference detector
    
    // Correlations between TPC and VZERO event planes
    kRPXtpcXvzeroa    = kVZEROTPCRPres+6*2,          
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
    // TPC event plane obtained from the precomputed Q vector in the trees
    kTPCQvecXtree   = kTPCsubResCos+6,
    kTPCQvecYtree   = kTPCQvecXtree+6,
    kTPCRPtree      = kTPCQvecYtree+6,
    kTPCQvecXptWeightsTree = kTPCRPtree+6,
    kTPCQvecYptWeightsTree = kTPCQvecXptWeightsTree+6,
    kTPCRPptWeightsTree    = kTPCQvecYptWeightsTree+6,
    kTPCQvecXposTree   = kTPCRPptWeightsTree+6,
    kTPCQvecYposTree   = kTPCQvecXposTree+6,
    kTPCRPposTree      = kTPCQvecYposTree+6,
    kTPCQvecXnegTree   = kTPCRPposTree+6,
    kTPCQvecYnegTree   = kTPCQvecXnegTree+6,
    kTPCRPnegTree      = kTPCQvecYnegTree+6,
    // ZDC variables
    kZDCnEnergyCh   = kTPCRPnegTree+6,         // ZDCn energy in each channel
    kZDCpEnergyCh   = kZDCnEnergyCh+10,        // ZDCp energy in each channel
    // TZERO variables
    kTZEROAmplitudeCh = kZDCpEnergyCh+10,      // TZERO aplitudes in all channels
    kTZEROTOF         = kTZEROAmplitudeCh+26,  // TZERO TOF start times
    kTZEROTOFbest     = kTZEROTOF+3,           // TZERO TOF best start times
    kTZEROzVtx        = kTZEROTOFbest+3,       // TZERO event z vertex
    kTZEROstartTime,                           // TZERO event start time
    kTZEROpileup,                              // TZERO pileup flag
    kTZEROsatellite,                           // TZERO satellite flag
    // External Multiplicity estimators
    kMultEstimatorV0M,
    kMultEstimatorV0A,
    kMultEstimatorV0C,
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
    kMultEstimatorPercentileV0M,
    kMultEstimatorPercentileV0A,
    kMultEstimatorPercentileV0C,
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
    kEMCEGATriggered,
    kEMCEGAHighTriggered,
    kEtaBinForSPDtracklets,
    kMCNch,                                  // number of primary charged particles in the MC, in |eta|<1
    kMCNchNegSide,                     // number of primary charged particles in the MC, in -1<eta<0
    kMCNchPosSide,                     // number of primary charged particles in the MC, in 0<eta<1
    kMCNchSPDacc,                       // number of primary charged particles in the MC, in |eta|<1 but limited to the SPD acceptance
    kDiffNchSPDtrklts,
    kDiffNchSPDaccSPDtrklts,
    kRelDiffNchSPDtrklts,
    kRelDiffNchSPDaccSPDtrklts,
    kRelDiff2NchSPDtrklts,
    kRelDiff2NchSPDaccSPDtrklts,
    kSPDntrackletsInCurrentEtaBin,
    kNEventVars,                               // number of event variables  
    // Particle variables --------------------------------------
    // Common pair/track variables
    kPt=kNEventVars,
    kPtMC,
    kPt_weight,
    kPtMCfromLegs,             // MC truth pt computed using the decay leg kinematics
    kP,      
    kPMC,
    kPMCfromLegs,
    kPx,   
    kPxMC,
    kPxMCfromLegs,
    kPy,     
    kPyMC,
    kPyMCfromLegs,
    kPz,     
    kPzMC,
    kPzMCfromLegs,
    kTheta,
    kThetaMC,
    kThetaMCfromLegs,
    kEta,
    kEtaMC,
    kEtaMCfromLegs,
    kPhi,     
    kPhiMC,
    kPhiMCfromLegs,
    kCosNPhi,   
    kSinNPhi = kCosNPhi+6,
    kPtSquared = kSinNPhi+6,
    kOneOverSqrtPt,                   // one over square root of pT
    kMass,
    kMassMC,
    kMassMCfromLegs,
    kRap,
    kRapAbs,
    kRapMC,
    kRapMCAbs,
    kRapMCfromLegs,
    kPdgMC,
    kCharge = kPdgMC+4,
    kVZEROFlowVn,                     // v_n using VZERO RP
    kVZERODeltaPhiPsiN = kVZEROFlowVn+6*3,   // delta phi = phi - Psi  for VZERO event plane
    kTPCFlowVn=kVZERODeltaPhiPsiN+6*3,      // v_n using TPC RP
    kTPCDeltaPhiPsiN=kTPCFlowVn+6,         // delta phi = phi - Psi  for TPC event plane
    kVZEROFlowSine=kTPCDeltaPhiPsiN+6,      // sin(n*(phi-Psi)) using VZERO RP
    kTPCFlowSine=kVZEROFlowSine+6*3,  // sin(n*(phi-Psi)) using TPC RP
    kVZEROuQ = kTPCFlowSine+6,        // cosine term from the u*Q products from VZERO (harmonics 1-6; VZERO-A and VZERO-C)
    kVZEROuQsine = kVZEROuQ+6*2,      // sine terms from the u*Q products from VZERO (harmonics 1-6; VZERO-A and VZERO-C)
    kTPCuQ=kVZEROuQsine+6*2,          // cosine terms from the u*Q products from TPC (harmonics 1-6)
    kTPCuQsine=kTPCuQ+6,              // sine terms from the u*Q products from TPC (harmonics 1-6)
    // Pair-only variables
    kCandidateId=kTPCuQsine+6,
    kPairType,                  // 0 ++; 1 +-; 2 --    
    kPairTypeSPD,               // 2 (both); 1 (one) 0 (none) of the legs has an hit in the first SPD layer;     
    kMassV0,                    // masses for all 4 V0 assumptions (0-K0s, 1-Lambda, 2-ALambda, 3-Gamma)
    kPairChisquare=kMassV0+4,     
    kPairLxy,           
    kPseudoProperDecayTime,
    kPseudoProperDecayTimeMC,
    kPairOpeningAngle,  
    kPairPointingAngle, 
    kPairThetaCS,                // cos (theta*) in Collins-Soper frame       
    kPairPhiCS,                    // phi* in Collins-Soper frame
    kPairThetaHE,                // cos (theta*) in helicity frame       
    kPairPhiHE,                    // phi* in helicity frame
    kPairVZEROFlowNom,            // Nominator of combinatorial pair flow for VZERO, 6 harmonics, 3 VZERO sides (A, C and A&C) 
    kPairVZEROFlowDenom=kPairVZEROFlowNom+6*3,  // Denominator of combinatorial pair flow for VZERO, 6 harmonics, 3 VZERO sides (A, C and A&C) 
    kPairTPCFlowNom=kPairVZEROFlowDenom+6*3,
    kPairTPCFlowDenom=kPairTPCFlowNom+6,
    kPairVZEROFlowSPNom=kPairTPCFlowDenom+6,     // Nominator of combinatorial Scalar Product pair flow for VZERO, 6 harmonics, 3 VZERO sides (A, C and A&C) 
    kPairVZEROFlowSPDenom=kPairVZEROFlowSPNom+6*3,  // Denominator of combinatorial pair flow for VZERO, 6 harmonics, 3 VZERO sides (A, C and A&C) 
    kPairTPCFlowSPNom=kPairVZEROFlowSPDenom+6*3,
    kPairTPCFlowSPDenom=kPairTPCFlowSPNom+6,
    kPairQualityFlag=kPairTPCFlowSPDenom+6,
    kPairQualityFlag2,
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
    kPairEff,                     // pair efficiency
    kOneOverPairEff,             // 1 / pair efficiency (correction factor) 
    kOneOverPairEffSq,             // 1 / pair efficiency squared (correction factor)
    kPairLegITSchi2,              // the ITS chi2 for the pair legs, used in correlations between pair legs
    kPairLegTPCchi2=kPairLegITSchi2+2,              // the TPC chi2 for the pair legs, used in correlations between pair legs
    kPairLegPt=kPairLegTPCchi2+2,                 // pair leg pt
    kPairLegPtSum=kPairLegPt+2,                   // sum of the pt of the two legs
    kPairLegPtMC,                                // MC truth pair leg pt
    kPairLegPtMCSum=kPairLegPtMC+2,               // sum of the MC truth leg pt's
    kPairLegEMCALmatchedEnergy,                   // pair leg EMCal cluster energy

    // Track-only variables -------------------------------------
    kPtTPC=kPairLegEMCALmatchedEnergy+2,
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
    kITSlayerShared,
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
    kTPCActiveLength,
    kTPCGeomLength,
    kTPCsignal,         
    kTPCsignalN,
    kTPCdEdxQmax,                 // dEdx info from Qmax (IROC, medium OROC, long OROC, all OROC)
    kTPCdEdxQtot=kTPCdEdxQmax+4,  // dEdx info from Qtot (IROC, medium OROC, long OROC, all OROC)
    kTPCdEdxQmaxOverQtot=kTPCdEdxQtot+4,    // Qmax / Qtot
    kTPCnSig=kTPCdEdxQmaxOverQtot+4,  
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
    kEMCALmatchedM02,
    kEMCALmatchedM20,
    kEMCALmatchedNCells,
    kEMCALmatchedNMatchedTracks,
    kEMCALmatchedDeltaPhi,
    kEMCALmatchedDeltaEta,
    kEMCALmatchedDistance,
    kEMCALmatchedNSigmaElectron,
    kNTrackVars,            // variable to mark end of track vars, introduce new tracks vars before this one
    // Calorimeter cluster variables --------------------------------------
    kEMCALclusterEnergy,        
    kEMCALclusterDx,            
    kEMCALclusterDz,            
    kEMCALdetector,         // 0 - EMCAL; 1 - PHOS     
    kEMCALm20,
    kEMCALm02,
    kEMCALdispersion,
    kEMCALnCells,
    kEMCALnMatchedTracks,
    kEMCALclusterPhi,
    kEMCALclusterEta,
    kNEMCALvars,            // variable to mark end of EMCal vars, introduce new EMCal vars before this one
    // Track flags -----------------------------------------------------
    kTrackingFlag,
    kTrackQualityFlag,
    kTrackQualityFlag2,
    kTrackMCFlag,
    kTrackMCFlag2,
    // Correlation variables ----------------------------------------------
    kDeltaPhi,        // shifted to [-pi/2, 3/2 * pi]
    kDeltaPhiBoosted, // after boost of associated track to trigger rest fram
    kDeltaPhiSym,     // shifted to [0, pi]
    kDeltaPhiSymBoosted,
    kDeltaTheta,
    kDeltaThetaBoosted,
    kDeltaEta,
    kDeltaEtaBoosted,
    kDeltaEtaAbs,
    kDeltaEtaAbsBoosted,
    kTriggerPt,       // pt of J/psi candidate
    kTriggerRap,      // rapidity of J/psi candidate
    kTriggerRapAbs,   // absolute rapidity of J/psi candidate
    kTriggerPseudoProperDecayTime,  // pseudo-proper decay length of J/psi candidate
    kTriggerPairTypeSPD,            // SPD pair type of J/psi candidate
    kAssociatedPt,          // pt of associated track
    kAssociatedPtBoosted,   // pt of associated track, after boost to trigger rest frame
    kAssociatedPtOverTriggerGammaT, // pt of associated track / transverse gamma of J/psi candidate
    kTriggerGammaT,                 // transverse gamma of J/psi candidate
    kAssociatedEta,         // eta of associated track
    kAssociatedEtaBoosted,
    kAssociatedPhi,         // phi of associated track
    kAssociatedPhiBoosted,
    kTriggerEff,                            // J/psi candidate efficiency
    kOneOverTriggerEff,                     // 1 / J/psi candidate efficiency
    kAssocHadronEff,                        // associated hadron efficiency
    kOneOverAssocHadronEff,                 // 1 / associated hadron efficiency (correction factor)
    kTriggerEffTimesAssocHadronEff,         // J/psi candidate efficiency x associated hadron efficiency
    kOneOverTriggerEffTimesAssocHadronEff,  // 1 / (J/psi candidate efficiency x associated hadron efficiency)
    // vars related to psiprime decay into J/psi + pi+pi- channel
    kAssociated2Pt,
    kAssociated2Eta,
    kAssociated2Phi,

    kOpAngleMotherPosPion,
    kDeltaRPosPi,
    kPPosPi,
    kPtPosPi,
    kDeltaRPosPiJPsi,

    kOpAngleMotherNegPion,
    kDeltaRNegPi,
    kPNegPi,
    kPtNegPi,
    kDeltaRNegPiJPsi,
    
    kOpAngleMotherJPsi,
    kDeltaRJPsi,
    kPJPsi,
    kPtJPsi,
    kEtaJPsi,
    
    kMassPionPair,
    kP_PionPair,
    kPt_PionPair,
    kPhi_PionPair,
    kEta_PionPair,
    kOpAngleMotherDiPion,
    kDeltaRDiPion,
    
    kMassPsiPrime,
    kPt_PsiPrime, 
    kPhi_PsiPrime,
    kEta_PsiPrime,

    kJPsiPosPionOpeningAngle,
    kJPsiNegPionOpeningAngle,
    kPionsOpeningAngle,
    kJPsiDiPionOpeningAngle,
    
    kQValue,

    kMassElecPairPosPion,
    kPt_ElecPairPosPion,
    kEta_ElecPairPosPion,
    kPhi_ElecPairPosPion,
    
    kMassPsiPrime_II,
    kPt_PsiPrime_II,
    kEta_PsiPrime_II,
    kJPsiPosPion_NegPionOpeningAngle,

    // TRD GTU online tracks
    kTRDGTUtracklets,   // TRD online track #tracklets
    kTRDGTUlayermask,   // TRD online track hit in layer0 yes/no
    kTRDGTUpt,          // TRD online track pT
    kTRDGTUsagitta,     // TRD online track sagitta
    kTRDGTUPID,         // TRD online track pid
    kTrackingFlags,
    kTRDTriggeredType,
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
  
  static const Double_t fgkSPDEtaCutsVsVtxZ[20][2];      // eta interval coverage for the SPDntracklets estimator as a function of vtx 
    
  AliReducedVarManager();
  AliReducedVarManager(const Char_t* name);
  virtual ~AliReducedVarManager();
  
  static void SetBeamMomentum(Float_t beamMom) {fgBeamMomentum = beamMom;}
  static Float_t GetBeamMomentum() {return fgBeamMomentum;}
  
  static void SetEvent(AliReducedBaseEvent* const ev) {fgEvent = ev;};
  static void SetEventPlane(AliReducedEventPlaneInfo* const ev) {fgEventPlane = ev;};
  static void SetUseVariable(Int_t var) {fgUsedVars[var] = kTRUE; SetVariableDependencies();}
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
  static void FillEventOnlineTrigger(UShort_t triggerBit, Float_t* values, UShort_t triggerBit2=999);
  static void FillEventTagInput(AliReducedBaseEvent* event, Int_t input, Float_t* values);
  static void FillL0TriggerInputs(AliReducedEventInfo* event, Int_t input, Float_t* values, Int_t input2=999);
  static void FillL1TriggerInputs(AliReducedEventInfo* event, Int_t input, Float_t* values, Int_t input2=999);
  static void FillL2TriggerInputs(AliReducedEventInfo* event, Int_t input, Float_t* values, Int_t input2=999);
  static void FillV0Channel(Int_t ich, Float_t* values);
  static void FillTrackingFlag(AliReducedTrackInfo* track, UInt_t flag, Float_t* values);
  static void FillTrackQualityFlag(AliReducedBaseTrack* track, UShort_t flag, Float_t* values, UShort_t flag2=999);
  static void FillTrackMCFlag(AliReducedBaseTrack* track, UShort_t flag, Float_t* values, UShort_t flag2=999);
  static void FillPairQualityFlag(AliReducedPairInfo* p, UShort_t flag, Float_t* values, UShort_t flag2=999);
  static void FillTrackInfo(AliReducedBaseTrack* p, Float_t* values);
  static void FillClusterMatchedTrackInfo(AliReducedBaseTrack* p, Float_t* values, TList* clusterList=0x0, AliReducedCaloClusterTrackMatcher* matcher=0x0);
  static void FillITSlayerFlag(AliReducedTrackInfo* track, Int_t layer, Float_t* values);
  static void FillITSsharedLayerFlag(AliReducedTrackInfo* track, Int_t layer, Float_t* values);
  static void FillTPCclusterBitFlag(AliReducedTrackInfo* track, Int_t bit, Float_t* values);
  static void FillPairInfo(AliReducedPairInfo* p, Float_t* values);
  static void FillPairInfo(AliReducedBaseTrack* t1, AliReducedBaseTrack* t2, Int_t type, Float_t* values);
  static void FillPairInfo(AliReducedPairInfo* leg1, AliReducedBaseTrack* leg2, Int_t type, Float_t* values);
  static void FillPairInfoME(AliReducedBaseTrack* t1, AliReducedBaseTrack* t2, Int_t type, Float_t* values);
  static void FillPairMEflow(AliReducedBaseTrack* t1, AliReducedBaseTrack* t2, Float_t* values/*, Int_t idx=0*/);
  static void FillCorrelationInfo(AliReducedBaseTrack* p, AliReducedBaseTrack* t, Float_t* values);
  static void FillCaloClusterInfo(AliReducedCaloClusterInfo* cl, Float_t* values);
  static void FillPsiPrimeInfo(AliReducedBaseTrack* p, AliReducedBaseTrack* t1, AliReducedBaseTrack* t2, Float_t* values);//tariq
  static void FillTrackingStatus(AliReducedTrackInfo* p, Float_t* values);
 // static void FillTrackingFlags(AliReducedTrackInfo* p, Float_t* values);
  static void FillMCTruthInfo(AliReducedTrackInfo* p, Float_t* values, AliReducedTrackInfo* leg1 = 0x0, AliReducedTrackInfo* leg2 = 0x0);
  static void FillMCTruthInfo(AliReducedTrackInfo* leg1, AliReducedTrackInfo* leg2, Float_t* values);
  static void FillMCEventInfo(AliReducedEventInfo* event, Float_t* values);
  
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
  static void SetTPCpidCalibMaps(Int_t pid, THnF* centroidMap, THnF* widthMap, THnI* statusMap);
  static void SetTPCpidCalibDepVars(Variables vars[]);
  static void SetPairEfficiencyMap(TH1* map, Variables varX, Variables varY=kNothing, Variables varZ=kNothing);
  static void SetPairEfficiencyMapDependeciesCorrelation(Variables varX, Variables varY=kNothing, Variables varZ=kNothing);
  static void SetAssociatedHadronEfficiencyMap(TH1* map, Variables varX, Variables varY=kNothing, Variables varZ=kNothing);
  static void SetLHCDataInfo(TH1F* totalLumi, TH1F* totalInt0, TH1F* totalInt1, TH1I* fillNumber);
  static void SetGRPDataInfo(TH1I* dipolePolarity, TH1I* l3Polarity, TH1I* timeStart, TH1I* timeStop);
  static void SetupGRPinformation(const Char_t* filename);
  static void SetRunNumbers( TString runNumbers );
  static void SetMultiplicityProfile( TH2* profile, Int_t estimator );
  static void SetVZEROCalibrationPath(const Char_t* path);
  static void SetCalibrateVZEROqVector(Bool_t option);
  static void SetRecenterVZEROqVector(Bool_t option);
  static void SetRecenterTPCqVector(Bool_t option);
  static void SetEventResolution(Bool_t option);
  static Int_t GetCorrectedMultiplicity( Int_t estimator = kMultiplicity, Int_t correction = 0, Int_t reference = 0, Int_t smearing = 0 );
  static void SetWeightSpectrum(TH1F *gReweightMCpt);
  static Double_t CalculateWeightFactor(Double_t Mcpt, Double_t Centrality);
  
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
  static AliKFParticle BuildKFcandidate(AliReducedTrackInfo* track1, Float_t mh1, AliReducedTrackInfo* track2, Float_t mh2);
  static AliKFParticle BuildKFvertex( AliReducedEventInfo * event );
  
  static TH2F* fgTPCelectronCentroidMap;    // TPC electron centroid 2D map
  static TH2F* fgTPCelectronWidthMap;       // TPC electron width 2D map
  static Variables fgVarDependencyX;        // varX in the 2-D electron correction maps
  static Variables fgVarDependencyY;        // varY in the 2-D electron correction maps
  static THnF* fgTPCpidCalibCentroid[3];    // TPC calib centroid 4D maps; [0] - electron, [1] - pion, [2] - proton
  static THnF* fgTPCpidCalibWidth[3];       // TPC calib width 4D map
  static THnI* fgTPCpidCalibStatus[3];      // TPC calib status 4D map
  static Variables fgTPCpidCalibVars[4];    // variables used for TPC pid 4D calibration
  static TH1*      fgPairEffMap;                  // pair efficiency map
  static Variables fgEffMapVarDependencyX;        // varX in the pair eff maps
  static Variables fgEffMapVarDependencyY;        // varY in the pair eff maps
  static Variables fgEffMapVarDependencyZ;        // varZ in the pair eff maps
  static Variables fgEffMapVarDependencyXCorr;        // varX in the pair eff maps, used for correlation
  static Variables fgEffMapVarDependencyYCorr;        // varY in the pair eff maps, used for correlation
  static Variables fgEffMapVarDependencyZCorr;        // varZ in the pair eff maps, used for correlation
  static TH1*      fgAssocHadronEffMap;               // assoc hadron efficiency map
  static Variables fgAssocHadronEffMapVarDependencyX; // varX in assoc hadron eff map
  static Variables fgAssocHadronEffMapVarDependencyY; // varY in assoc hadron eff map
  static Variables fgAssocHadronEffMapVarDependencyZ; // varZ in assoc hadron eff map
  
  static TH1F* fgRunTotalLuminosity;      // total luminosity, GRP/GRP/LHCData::GetLumiAliceSBDelivered()
  static TH1F* fgRunTotalIntensity0;        // total intensity beam 1, GRP/GRP/LHCData::GetTotalIntensity(0)
  static TH1F* fgRunTotalIntensity1;         // total intensity beam 2, GRP/GRP/LHCData::GetTotalIntensity(1)
  static TH1I* fgRunLHCFillNumber;         // LHC fill number, GRP/GRP/LHCData::GetFillNumber()
  static TH1I* fgRunDipolePolarity;          // dipole magnet polarity, GRP/GRP/Data::GetDipolePolarity()
  static TH1I* fgRunL3Polarity;                // L3 magnet polarity, GRP/GRP/Data::GetL3Polarity()
  static TH1I* fgRunTimeStart;                // run start time, GRP/GRP/Data::GetTimeStart()
  static TH1I* fgRunTimeEnd;                  // run stop time, GRP/GRP/Data::GetTimeEnd()
  static TFile* fgGRPfile;                    // file containing GRP information
  static TGraphErrors* fgRunInstLumi;         // time dependence of the instantaneous lumi for the current run, AliLumiTools::GetLumiFromCTP(run)
  static std::vector<Int_t> fgRunNumbers;     // vector with run numbers (for histograms vs. run number)
  static Int_t fgRunID;                       // run ID
  static TH1* fgAvgMultVsVtxGlobal      [kNMultiplicityEstimators];        // average multiplicity vs. z-vertex position (global)
  static TH1* fgAvgMultVsVtxRunwise     [kNMultiplicityEstimators];        // average multiplicity vs. z-vertex position (run-by-run)
  static TH1* fgAvgMultVsRun            [kNMultiplicityEstimators];           // average multiplicity vs. run number
  static TH2* fgAvgMultVsVtxAndRun      [kNMultiplicityEstimators];  // 2D : average multiplicity vs. run number and z-vertex position
  static Double_t fgRefMultVsVtxGlobal  [kNMultiplicityEstimators] [kNReferenceMultiplicities];  // reference multiplicity for z-vertex correction (global)
  static Double_t fgRefMultVsVtxRunwise [kNMultiplicityEstimators] [kNReferenceMultiplicities];  // reference multiplicity for z-vertex correction (run-by-run)
  static Double_t fgRefMultVsRun        [kNMultiplicityEstimators] [kNReferenceMultiplicities];  // reference multiplicity for run correction
  static Double_t fgRefMultVsVtxAndRun  [kNMultiplicityEstimators] [kNReferenceMultiplicities];  // reference multiplicity for run, vtx correction
  static TString fgVZEROCalibrationPath;       // path to the VZERO calibration histograms
  static TProfile2D* fgAvgVZEROChannelMult[64];       // average multiplicity in VZERO channels vs (vtxZ,centSPD)
  static TProfile2D* fgVZEROqVecRecentering[4];       // (vtxZ,centSPD) maps of the VZERO A and C recentering Qvector offsets
  static TProfile2D* fgTPCqVecRecentering[2];       // (vtxZ,centV0) maps of the TPC recentering Qvector offsets
  static Bool_t fgOptionCalibrateVZEROqVec;         //option to calibrate V0
  static Bool_t fgOptionRecenterVZEROqVec;         //option to do Q vector recentering for V0
  static Bool_t fgOptionRecenterTPCqVec;           //option to do Q vector recentering for TPC
  static Bool_t fgOptionEventRes;                 //option to divide by resolution
  static TH1F* fgReweightMCpt;  // ratio between nature pt shape and gernareted pT shape
  
  AliReducedVarManager(AliReducedVarManager const&);
  AliReducedVarManager& operator=(AliReducedVarManager const&);  
  
  ClassDef(AliReducedVarManager, 17);
};

#endif
