#ifndef ALIALIQNCORRECTIONS_VARMANAGER_H
#define ALIALIQNCORRECTIONS_VARMANAGER_H

/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/
/*  Based on work of Ionut-Cristian Arsene
 ************************************************************/
#include "TNamed.h"
#include "Rtypes.h"



class AliVEvent;

class AliQnCorrectionsVarManagerTask : public AliAnalysisTaskSE {

public:
  AliQnCorrectionsVarManagerTask();
  AliQnCorrectionsVarManagerTask(const char *name);
  virtual ~AliQnCorrectionsVarManagerTask();

  enum Detector {
    kVZERO=0,
    kTPC,    // 1
    kZDC,    // 2
    kTZERO,   // 3
    kFMD,    // 4
    kFMDraw,    // 5
    kSPD,    // 6
    kNdetectors
  };

  enum ParticleId {
    kUnknown = -1,
    kElectron = 0,
    kPion,
    kKaon,
    kProton
  };

  enum Variables {
    kNothing = -1,
    // Event wise variables
    kRandom1 = 0,   // slot for in-macro assignment
    kRandom2,       // slot for in-macro assignment
    kRunNo,         // run number         
    kLHCFillNumber,     // LHC fill number
    kBeamEnergy,        // LHC beam energy
    kDetectorMask,      // detector mask
    kNumberOfDetectors, // number of active detectors
    kDipolePolarity,    // dipole magnet polarity
    kL3Polarity,        // L3 magnet polarity
    kBeamIntensity0,    // beam 0 intensity
    kBeamIntensity1,    // beam 1 intensity
    kBC,                // bunch crossing     
    kTimeStamp,         // time stamp of the event
    kEventType,         // event type
    kTriggerMask,       // trigger mask       
    kOfflineTrigger,    // offline trigger    
    kOfflineTriggerFired,  // offline trigger fired
    kOfflineTriggerFired2,  // offline trigger if fired, -1 if not fired
    kIsPhysicsSelection,    // physics selection 
    kIsSPDPileup,          // whether is SPD pileup
    kNSPDpileups,         // number of pileup events from SPD
    kNTrackPileups,       // number of pileup events from tracks
    kNPMDtracks,          // number of PMD tracks
    kNTRDtracks,          // number of TRD tracks
    kNTRDtracklets,          // number of TRD tracklets
    kNVtxContributors,    // number of vertex contributors
    kNVtxTPCContributors,  // number of TPC vertex contributors
    kVtxX,              // vtx X                      
    kVtxY,              // vtx Y                      
    kVtxZ,              // vtx Z 
    kVtxXtpc,           // vtx X from tpc
    kVtxYtpc,           // vtx Y from tpc
    kVtxZtpc,           // vtx Z from tpc
    kDeltaVtxZ,         // vtxZ - vtxZtpc
    kNTracksPerTrackingFlag,  // number of tracks with a given tracking flag
    kNTracksTPCoutVsITSout=kNTracksPerTrackingFlag+32,   //  TPCout/ITSout
    kNTracksTRDoutVsITSout,                              //  TRDout/ITSout
    kNTracksTOFoutVsITSout,                              //  TOFout/ITSout
    kNTracksTRDoutVsTPCout,                              //  TRDout/TPCout
    kNTracksTOFoutVsTPCout,                              //  TOFout/TPCout
    kNTracksTOFoutVsTRDout,                              //  TOFout/TRDout
    kNTracksITSoutVsSPDtracklets,                        //  ITSout/SPDtracklets
    kNTracksTPCoutVsSPDtracklets,                        //  TPCout/SPDtracklets
    kNTracksTRDoutVsSPDtracklets,                        //  TRDout/SPDtracklets
    kNTracksTOFoutVsSPDtracklets,                        //  TOFout/SPDtracklets
    kVZEROMultPercentile,  // multiplicity percentile for VZERO (PbPb2015)
    kCentVZERO,         // centrality from VZERO
    kCentSPD,           // centrality from SPD  
    kCentSPDcorr,       // corrected centrality from SPD
    kCentTPC,           // centrality from TPC  
    kCentZDC,           // centrality from ZDC  
    kCentQuality,       // centrality quality   
    kNV0total,          // total number of V0s in the esd       
    kNV0selected,       // number of V0s selected               
    kNdielectrons,      // number of dielectron pairs           
    kNpairsSelected,    // number of selected dielectron pairs per event        
    kNtracksTotal,      // total number of tracks               
    kNtracksSelected,   // number of selected tracks            
    kNtracksPosAnalyzed,// number of positive tracks passing analysis cuts      
    kNtracksNegAnalyzed,// number of negative tracks passing analysis cuts      
    kNtracksAnalyzed,   // number of positive+negative tracks passing analysis cuts
    kNtracksSubEvLeft,  // number of tracks in the left sub-event (negative pseudo-rapidity)            
    kNtracksSubEvRight, // number of tracks in the left sub-event (positive pseudo-rapidity)            
    kNtracksEventPlane, // number of tracks used for event plane                
    kSPDnSingleClusters,// SPD number of single clusters
    kSPDntracklets,     // SPD number of tracklets in |eta|<1.0                 
    kSPDntrackletsCorr, // SPD number of tracklets in |eta|<1.0                 
    kSPDntrackletsEta,  // SPD number of tracklets in -1.6+0.2*i < eta < -1.6+0.2*(i+1)
    kSPDtrackletEta,    // SPD tracklet eta
    kSPDtrackletPhi,    // SPD tracklet phi
    kNFMD1channels,      // FMD number of channels with nonzero multiplicity
    kNFMD2Ichannels,      // FMD number of channels with nonzero multiplicity
    kNFMD2Ochannels,      // FMD number of channels with nonzero multiplicity
    kNFMD3Ichannels,      // FMD number of channels with nonzero multiplicity
    kNFMD3Ochannels,      // FMD number of channels with nonzero multiplicity
    kFMD1TotalMult,       // FMD multiplicity of FMD1
    kFMD2ITotalMult,       // FMD multiplicity of FMD2
    kFMD2OTotalMult,       // FMD multiplicity of FMD2
    kFMD3ITotalMult,       // FMD multiplicity of FMD2
    kFMD3OTotalMult,       // FMD multiplicity of FMD3
    kEventMixingId=kSPDntrackletsEta+16,     // Id of the event mixing category                      
    // VZERO event plane related variables
    kVZEROATotalMult,   // total multiplicity of VZEROA                         
    kVZEROCTotalMult,   // total multiplicity of VZEROC                         
    kVZEROTotalMult,    // total multiplicity of VZERO                          
    kVZEROAemptyChannels,  // Number of empty VZERO channels in A side          
    kVZEROCemptyChannels,  // Number of empty VZERO channels in C side          
    kVZEROChannel,      // For filling histograms at correct bin
    kVZEROChannel1,      // For filling histograms at correct bin
    kVZEROChannel2,      // For filling histograms at correct bin
    kVZEROChannelMulttmp,      // For filling histograms with correct weight
    kChannelMult,      // For filling histograms with correct weight
    kVZEROChannelMult,                        // VZERO multiplicity per channel                 
    kVZEROChannelEta = kVZEROChannelMult+64,  // pseudo-rapidity of a VZERO channel             
    // ZDC event plane related variables
    kZDCATotalEnergy = kVZEROChannelEta+64,  // total energy of ZDCA                         
    kZDCCTotalEnergy,   // total energy of ZDCC                         
    kZDCTotalEnergy,    // total energy of ZDC                          
    kZDCAemptyTowers,  // Number of empty ZDC towers in A side          
    kZDCCemptyTowers,  // Number of empty ZDC towers in C side          
    kZDCTower,      // For filling histograms at correct bin
    kZDCTowerEnergytmp,      // For filling histograms with correct weight
    kZDCTowerEnergy,                        // ZDC energy per tower                 
    kZDCTowerEta = kZDCTowerEnergy+10,  // pseudo-rapidity of a ZDC tower             
    // TZERO event plane variables
    kTZEROChannelMult = kZDCTowerEta+10,                       // TZERO multiplicity per channel                 
    kTZEROTotalMult=kTZEROChannelMult+24,                        // TZERO multiplicity                 
    kTZEROATotalMult,                        // TZERO multiplicity A-side                 
    kTZEROCTotalMult,                        // TZERO multiplicity C-side                 
    kTZEROAemptyChannels,  // Number of empty TZERO channels in A side          
    kTZEROCemptyChannels, 			 // Number of empty TZERO channels in C side          
    kTZEROChannelMulttmp,      // For filling histograms with correct weight
    kTZEROChannel,      // For filling histograms at correct bin
    kNEventVars,           // number of event variables
    // Track variables -------------------------------------
    kPt,               
    kP,         
    kPx,        
    kPy,        
    kPz,        
    kCharge,
    kTheta,     
    kPhi,       
    kCosNPhi,   
    kSinNPhi = kCosNPhi+6,      
    kCos2NPhi = kSinNPhi+6,
    kSin2NPhi = kCos2NPhi+6,
    kEta = kSin2NPhi+6,          
    kRap,       
    kPtTPC,     
    kPhiTPC,    
    kEtaTPC,    
    kPin,       
    kDcaXY,     
    kDcaZ,              
    kITSncls,           
    kITSlayerHit,       
    kITSsignal,         
    kITSnSig,
    kTPCncls=kITSnSig+4,    
    kTPCclusBitFired,   
    kTPCNclusBitsFired, 
    kTPCclustersPerBit, 
    kTPCcrossedRows,    
    kTPCnclsIter1,      
    kTPCnclsF,          
    kTPCnclsRatio,       // TPCncls / TPCnclsF          
    kTPCnclsRatio2,      // TPCncls / TPCCrossedRows    
    kTPCsignal,         
    kTPCnSig,           
    kTPCchi2=kTPCnSig+4,           
    kTPCchi2Iter1,
    kTOFbeta,
    kTOFnSig,                   
    kTRDntracklets=kTOFnSig+4,  
    kTRDntrackletsPID,          
    kTRDpidProbabilities,       
    kEMCALmatchedEnergy=kTRDpidProbabilities+2,         
    kEMCALmatchedEOverP,        
    // Calorimeter cluster variables --------------------------------------
    kEMCALclusterEnergy,        
    kEMCALclusterDx,            
    kEMCALclusterDz,            
    kEMCALdetector,         // 0 - EMCAL; 1 - PHOS      
    // FMD channel variables
    kFMDId,
    kFMDMultiplicity,
    kFMDEqualizedMultiplicity,
    kFMDEta,
    kFMDPhi,
    kFMDdetector,
    // Tracking flags -----------------------------------------------------
    kTrackingFlag,              
    // Correlation variables ----------------------------------------------
    kDeltaPhi,          
    kDeltaTheta,        
    kDeltaEta,          
    kFilterBit,         
    kFilterBitMask768=kFilterBit+9,
    kDetector,
    kQvectorX,
    kQvectorY,
    kNVars
  };


  //-----------Variable names and units to be used as default when not specified via AddHistogram() -------------

  // tracking flags as in AliESDtrack.h
  // NOTE: check consistency with aliroot
  enum TrackingFlags {
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
    kNTrackingFlags
  };

  enum FilterBits { 
    kFilterBit0 = 0, 
    kFilterBit1, 
    kFilterBit2, 
    kFilterBit3, 
    kFilterBit4, 
    kFilterBit5, 
    kFilterBit6, 
    kFilterBit7, 
    kFilterBit8        
  };

  // offline triggers as defined in AliVEvent.h
  // NOTE: Check consistency with updates in aliroot!!!
  enum EOfflineTriggerTypes { 
    kMB                = BIT(0), // Minimum bias trigger, i.e. interaction trigger, offline SPD or V0 selection
    kINT7              = BIT(1), // V0AND trigger, offline V0 selection
    kMUON              = BIT(2), // Muon trigger, offline SPD or V0 selection
    kHighMult          = BIT(3), // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
    kEMC1              = BIT(4), // EMCAL trigger
    kCINT5             = BIT(5), // Minimum bias trigger without SPD. i.e. interaction trigger, offline V0 selection
    kCMUS5             = BIT(6), // Muon trigger, offline V0 selection
    kMUSPB             = BIT(6), // idem for PbPb
    kMUSH7             = BIT(7), // Muon trigger: high pt, single muon, offline V0 selection, CINT7 suite
    kMUSHPB            = BIT(7), // idem for PbPb
    kMUL7              = BIT(8), // Muon trigger: like sign dimuon, offline V0 selection, CINT7 suite
    kMuonLikePB        = BIT(8), // idem for PbPb
    kMUU7              = BIT(9), // Muon trigger, unlike sign dimuon, offline V0 selection, CINT7 suite
    kMuonUnlikePB      = BIT(9), // idem for PbPb
    kEMC7              = BIT(10), // EMCAL trigger, CINT7 suite
    kMUS7              = BIT(11), // Muon trigger: low pt, single muon, offline V0 selection, CINT7 suite
    kPHI1              = BIT(12), // PHOS trigger, CINT1 suite
    kPHI7              = BIT(13), // PHOS trigger, CINT7 suite
    kPHOSPb            = BIT(13), // idem for PbPb
    kEMCEJE            = BIT(14), // EMCAL jet patch trigger
    kEMCEGA            = BIT(15), // EMCAL gamma trigger
    kCentral           = BIT(16), // PbPb central collision trigger
    kSemiCentral       = BIT(17), // PbPb semicentral collision trigger
    kDG5               = BIT(18), // Double gap diffractive
    kZED               = BIT(19), // ZDC electromagnetic dissociation
    kSPI7              = BIT(20), // Power interaction trigger
    kINT8              = BIT(21), // CINT8 trigger: 0TVX (T0 vertex) triger
    kMuonSingleLowPt8  = BIT(22), // Muon trigger : single muon, low pt, T0 selection, CINT8 suite
    kMuonSingleHighPt8 = BIT(23), // Muon trigger : single muon, high pt, T0 selection, CINT8 suite
    kMuonLikeLowPt8    = BIT(24), // Muon trigger : like sign muon, low pt, T0 selection, CINT8 suite
    kMuonUnlikeLowPt8  = BIT(25), // Muon trigger : unlike sign muon, low pt, T0 selection, CINT8 suite
    kINT6              = BIT(26),
    kUserDefined       = BIT(27), // Set when custom trigger classes are set in AliPhysicsSelection, offline SPD or V0 selection
    // Bits 28 and above reserved for FLAGS
    kFastOnly          = BIT(30), // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB
    kAny               = 0xffffffff, // to accept any trigger
    kAnyINT            = kMB | kINT7 | kCINT5 // to accept any interaction (aka minimum bias) trigger
  };

  enum ITSLayerMap {
    kITSfirst  =  1,
    kITSsecond =  2,
    kITSthird  =  4,
    kITSfourth =  8,
    kITSfifth  = 16,
    kITSsixth  = 32
  };

public:
  virtual void UserExec(Option_t *) = 0;
  virtual void UserCreateOutputObjects() = 0;
  virtual void FinishTaskOutput() = 0;

  const Char_t* VarName(Int_t var) const;
  const Char_t* VarUnits(Int_t var) const;
protected:
  void SetDefaultVarNames();

public:
  static const Char_t* fTrackingFlagNames[kNTrackingFlags];
  static const Char_t* fOfflineTriggerNames[64];
  const Char_t* fVariableNames[kNVars][2];      //!<! The variable names. Transient!

private:

  AliQnCorrectionsVarManagerTask(const AliQnCorrectionsVarManagerTask &c);
  AliQnCorrectionsVarManagerTask & operator= (const AliQnCorrectionsVarManagerTask &c);

  ClassDef(AliQnCorrectionsVarManagerTask, 1);
};  

inline const Char_t* AliQnCorrectionsVarManagerTask::VarName(Int_t var) const {
  if ((!(var < 0)) && (var < kNVars))
    return fVariableNames[var][0];
  else
    return "";
}

inline const Char_t* AliQnCorrectionsVarManagerTask::VarUnits(Int_t var) const {
  if ((!(var < 0)) && (var < kNVars))
    return fVariableNames[var][1];
  else
    return "";
}


#endif    
