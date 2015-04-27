/*
***********************************************************
    Variable definitions for event plane correction framework
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
    Based on work of Ionut-Cristian Arsene
***********************************************************
*/
#ifndef ALIEVENTPLANEVARMANAGER_H
#define ALIEVENTPLANEVARMANAGER_H

#include <iostream>
#include <fstream>


//#include <AliCentrality.h>
//#include <AliVEvent.h>
//#include <AliESDEvent.h>
//#include <AliESDHeader.h>
//#include <AliTriggerAnalysis.h>
//#include <AliMultiplicity.h>
//#include <AliESDtrack.h>
//#include "AliReducedEvent.h"
#include "AliEventPlaneQvector.h"
//#include "AliEventPlaneESDVarManager.h"

//#ifdef ALIEVENTPLANEESDVARMANAGER_CXX
//#define FRIEND AliEventPlaneESDVarManager
//#endif

//class AliReducedEvent;
class AliEventPlaneQvector;

#ifdef ALIREDUCEDEVENT_H
#define EVENT AliReducedEvent
//#define FRIEND AliReducedEventFriend
#define TRACK AliReducedTrack
#define FMD AliReducedFMD
#define CLUSTER AliReducedCaloCluster
#endif

class AliVEvent;

//_________________________________________________________
class AliEventPlaneVarManager : public TNamed {
//namespace AliEventPlaneVarManager{

  public:
  AliEventPlaneVarManager();
  virtual ~AliEventPlaneVarManager();


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
    kCentVZERO,         // centrality from VZERO
    kCentSPD,           // centrality from SPD  
    kCentSPDcorr,       // corrected centrality from SPD
    kCentTPC,           // centrality from TPC  
    kCentZDC,           // centrality from ZDC  
    kCentQuality,       // centrality quality   
    kNV0total,          // total number of V0s in the esd       //18
    kNV0selected,       // number of V0s selected               //19
    kNdielectrons,      // number of dielectron pairs           //20
    kNpairsSelected,    // number of selected dielectron pairs per event        //21
    kNtracksTotal,      // total number of tracks               //22
    kNtracksSelected,   // number of selected tracks            //23
    kNtracksPosAnalyzed,// number of positive tracks passing analysis cuts      //24
    kNtracksNegAnalyzed,// number of negative tracks passing analysis cuts      //25
    kNtracksAnalyzed,   // number of positive+negative tracks passing analysis cuts
    kNtracksSubEvLeft,  // number of tracks in the left sub-event (negative pseudo-rapidity)            //26
    kNtracksSubEvRight, // number of tracks in the left sub-event (positive pseudo-rapidity)            //27
    kNtracksEventPlane, // number of tracks used for event plane                //28
    kSPDntracklets,     // SPD number of tracklets in |eta|<1.0                 //29
    kSPDntrackletsCorr, // SPD number of tracklets in |eta|<1.0                 //29
    kSPDntrackletsEta,  // SPD number of tracklets in -1.6+0.2*i < eta < -1.6+0.2*(i+1)
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
    kEventMixingId=kSPDntrackletsEta+16,     // Id of the event mixing category                      //30
    // VZERO event plane related variables
    kVZEROATotalMult,   // total multiplicity of VZEROA                         //31
    kVZEROCTotalMult,   // total multiplicity of VZEROC                         //32
    kVZEROTotalMult,    // total multiplicity of VZERO                          //33
    kVZEROAemptyChannels,  // Number of empty VZERO channels in A side          //34
    kVZEROCemptyChannels,  // Number of empty VZERO channels in C side          //35
    kVZEROChannel,      // For filling histograms at correct bin
    kVZEROChannel1,      // For filling histograms at correct bin
    kVZEROChannel2,      // For filling histograms at correct bin
    kVZEROChannelMulttmp,      // For filling histograms with correct weight
    kChannelMult,      // For filling histograms with correct weight
    kVZEROChannelMult,                        // VZERO multiplicity per channel                 //36
    kVZEROChannelEta = kVZEROChannelMult+64,  // pseudo-rapidity of a VZERO channel             //100
    kVZEROQvecX      = kVZEROChannelEta+64,   // Q-vector components for harmonics 1-6 and      //164
    kVZEROQvecY      = kVZEROQvecX+6*11,        //  6- n-harmonics; 3- A,C and A&C options       //182
    kVZEROQvecMag    = kVZEROQvecY+6*11,       // magnitude of the Q vector
    kVZERORP         = kVZEROQvecMag+6*11,      // VZERO reaction plane from A,C and A&C sides (harmonics 1-6)        //200
    kVZERORPres      = kVZERORP+6*11,     // VZERO reaction plane resolution (sqrt(n*(RPa-RPc)))         //218
    kRPres           = kVZERORPres+6,     // Allocate 6 slots for user defined RP resolution using 2 detectors of choice
    kVZEROXaXc       = kRPres+6,           // correlations for the components of the Q vector      //224
    kVZEROXaYa       = kVZEROXaXc+6,                            //230
    kVZEROXaYc       = kVZEROXaYa+6,                            //236
    kVZEROYaXc       = kVZEROXaYc+6,                            //242
    kVZEROYaYc       = kVZEROYaXc+6,                            //248
    kVZEROXcYc       = kVZEROYaYc+6,                            //254
    kVZEROdeltaRPac  = kVZEROXcYc+6,         // Psi_VZEROA-Psi_VZEROC           //260
    kVZEROdeltaRPa   = kVZEROdeltaRPac+6,      // Psi_VZEROA6-Psi_VZEROA5, Psi_VZEROA7-Psi_VZEROA5, Psi_VZEROA6-Psi_VZEROA5 //266
    kVZEROdeltaRPc   = kVZEROdeltaRPa+6*2,     // Psi_VZEROA2-Psi_VZEROA1, Psi_VZEROA3-Psi_VZEROA1, Psi_VZEROA4-Psi_VZEROA1 //278
    kVZEROflowV2TPC  = kVZEROdeltaRPc+6*2,     // vzero v2 using TPC event plane        //290
        // ZDC event plane related variables
    kZDCATotalEnergy = kVZEROflowV2TPC+64,  // total energy of ZDCA                         //31
    kZDCCTotalEnergy,   // total energy of ZDCC                         //32
    kZDCTotalEnergy,    // total energy of ZDC                          //33
    kZDCAemptyTowers,  // Number of empty ZDC towers in A side          //34
    kZDCCemptyTowers,  // Number of empty ZDC towers in C side          //35
    kZDCTower,      // For filling histograms at correct bin
    kZDCTowerEnergytmp,      // For filling histograms with correct weight
    kZDCTowerEnergy,                        // ZDC energy per tower                 //36
    kZDCTowerEta = kZDCTowerEnergy+10,  // pseudo-rapidity of a ZDC tower             //100
    kZDCQvecX      = kZDCTowerEta+10,   // Q-vector components for harmonics 1-6 and      //164
    kZDCQvecY      = kZDCQvecX+6*3,        //  6- n-harmonics; 3- A,C and A&C options       //182
    kZDCQvecMag    = kZDCQvecY+6*3,       // magnitude of the Q vector
    kZDCRP         = kZDCQvecMag+6*3,      // ZDC reaction plane from A,C and A&C sides (harmonics 1-6)        //200
    kZDCRPres      = kZDCRP+6*3,     // ZDC reaction plane resolution (sqrt(n*(RPa-RPc)))         //218
    kZDCXaXc       = kZDCRPres+6,           // correlations for the components of the Q vector      //224
    kZDCXaYa       = kZDCXaXc+6,                            //230
    kZDCXaYc       = kZDCXaYa+6,                            //236
    kZDCYaXc       = kZDCXaYc+6,                            //242
    kZDCYaYc       = kZDCYaXc+6,                            //248
    kZDCXcYc       = kZDCYaYc+6,                            //254
    kZDCdeltaRPac  = kZDCXcYc+6,         // Psi_ZDCA-Psi_ZDCC           //260
    kZDCdeltaRPa   = kZDCdeltaRPac+6,      // Psi_ZDCA6-Psi_ZDCA5, Psi_ZDCA7-Psi_ZDCA5, Psi_ZDCA6-Psi_ZDCA5 //266
    kZDCdeltaRPc   = kZDCdeltaRPa+6*2,     // Psi_ZDCA2-Psi_ZDCA1, Psi_ZDCA3-Psi_ZDCA1, Psi_ZDCA4-Psi_ZDCA1 //278
    kZDCflowV2TPC  = kZDCdeltaRPc+6*2,     // vzero v2 using TPC event plane        //290
    // TPC event plane variables
    kTPCQvecX = kZDCflowV2TPC+8,   // TPC Q-vector components for harmonics 1-6      //354
    kTPCQvecY = kTPCQvecX+6,                                                            //360
    kTPCRP    = kTPCQvecY+6,                // Event plane using TPC                    //366
    kTPCRPres = kTPCRP+6,                // Event plane resolution variables sqrt(n*(RPtpc-RPvzeroa)),sqrt(n*(RPtpc-RPvzeroc))  //372
    // TZERO event plane variables
    kTZEROChannelMult= kTPCRPres+6*2,                        // TZERO multiplicity per channel                 //36
    kTZEROTotalMult=kTZEROChannelMult+24,                        // TZERO multiplicity                 //36
    kTZEROATotalMult,                        // TZERO multiplicity A-side                 //36
    kTZEROCTotalMult,                        // TZERO multiplicity C-side                 //36
    kTZEROAemptyChannels,  // Number of empty TZERO channels in A side          //34
    kTZEROCemptyChannels, 			 // Number of empty TZERO channels in C side          //35
    kTZEROChannelMulttmp,      // For filling histograms with correct weight
    kTZEROChannel,      // For filling histograms at correct bin
    kTZERORP         ,          // TZERO A- and C-side, plus slot for total
    kTZEROQvecX      = kTZERORP+6*3,
    kTZEROQvecY      = kTZEROQvecX+6*3,
    // FMD event plane variables
    kFMDRP         = kTZEROQvecY+6*3,// FMD{1,2I,2O,3I,3O} plus slot for total    kFMDQvecX      = kFMDRP+6*4,
    kFMDQvecX      = kFMDRP+6*6,
    kFMDQvecY      = kFMDQvecX+6*6,
    // Correlations between TPC and VZERO event planes
    kRPXtpcXvzeroa    = kFMDQvecY+6*6,
    kRPXtpcXvzeroc    = kRPXtpcXvzeroa+6,       //390
    kRPYtpcYvzeroa    = kRPXtpcXvzeroc+6,       //396
    kRPYtpcYvzeroc    = kRPYtpcYvzeroa+6,       //402
    kRPXtpcYvzeroa    = kRPYtpcYvzeroc+6,       //408
    kRPXtpcYvzeroc    = kRPXtpcYvzeroa+6,       //414
    kRPYtpcXvzeroa    = kRPXtpcYvzeroc+6,       //420
    kRPYtpcXvzeroc    = kRPYtpcXvzeroa+6,       //426
    kRPdeltaVZEROAtpc = kRPYtpcXvzeroc+6,       //432
    kRPdeltaVZEROCtpc = kRPdeltaVZEROAtpc+6,    //438
    // TPC event plane using sub-intervals in pseudo-rapidity
    kTPCtmp   = kRPdeltaVZEROCtpc+6,      //444
    kTPCQvecXleft,
    kTPCQvecYleft   = kTPCQvecXleft+6,          //450
    kTPCRPleft      = kTPCQvecYleft+6,          //456
    kTPCQvecXright  = kTPCRPleft+6,             //462
    kTPCQvecYright  = kTPCQvecXright+6,         //468
    kTPCRPright     = kTPCQvecYright+6,         //474
    kTPCQvecXneg   = kTPCRPright+6,
    kTPCQvecYneg   = kTPCQvecXneg+6,
    kTPCRPneg      = kTPCQvecYneg+6,
    kTPCQvecXpos  = kTPCRPneg+6,
    kTPCQvecYpos  = kTPCQvecXpos+6,
    kTPCRPpos     = kTPCQvecYpos+6,
    kTPCQvecXran1  = kTPCRPpos+6,
    kTPCQvecYran1  = kTPCQvecXran1+6,
    kTPCRPran1     = kTPCQvecYran1+6,
    kTPCQvecXran2  = kTPCRPran1+6,
    kTPCQvecYran2  = kTPCQvecXran2+6,
    kTPCRPran2     = kTPCQvecYran2+6,
    kTPCQvecXtotal  = kTPCRPran2+6,
    kTPCQvecYtotal  = kTPCQvecXtotal+6,         //486
    kTPCRPtotal     = kTPCQvecYtotal+6,         //492
    kTPCsubResCos   = kTPCRPtotal+6,            //498
    // Event plane variables for macro defined detectors
    kDetAQvecX      = kTPCsubResCos+6,
    kDetAQvecY      = kDetAQvecX+6,
    kDetBQvecX      = kDetAQvecY+6,
    kDetBQvecY      = kDetBQvecX+6,
    kDetCQvecX      = kDetBQvecY+6,
    kDetCQvecY      = kDetCQvecX+6,
    kDetARP         = kDetCQvecY+6,
    kDetBRP         = kDetARP+6,
    kDetCRP         = kDetBRP+6,
    kDetAmult         = kDetCRP+6,
    kDetBmult         = kDetAmult+6,
    kDetCmult         = kDetBmult+6,
    kCosDetARP         = kDetCmult+6,
    kCosDetBRP         = kCosDetARP+6,
    kCosDetCRP         = kCosDetBRP+6,
    kSinDetARP         = kCosDetCRP+6,
    kSinDetBRP         = kSinDetARP+6,
    kSinDetCRP         = kSinDetBRP+6,
    kDetADetBDeltaRP      = kSinDetCRP+6,
    kDetBDetCDeltaRP      = kDetADetBDeltaRP+6,
    kDetADetCDeltaRP      = kDetBDetCDeltaRP+6,
    kDetADetBRPres      = kDetADetCDeltaRP+6,
    kDetBDetCRPres      = kDetADetBRPres+6,
    kDetADetCRPres      = kDetBDetCRPres+6,
    kDetADetBRPcosplus      = kDetADetCRPres+6,
    kDetBDetCRPcosplus      = kDetADetBRPcosplus+6,
    kDetADetCRPcosplus      = kDetBDetCRPcosplus+6,
    kDetADetBRPsinplus      = kDetADetCRPcosplus+6,
    kDetBDetCRPsinplus      = kDetADetBRPsinplus+6,
    kDetADetCRPsinplus      = kDetBDetCRPsinplus+6,
    kDetADetBRPsin      = kDetADetCRPsinplus+6,
    kDetBDetCRPsin      = kDetADetBRPsin+6,
    kDetADetCRPsin      = kDetBDetCRPsin+6,
    kDetAQxDetBQx    = kDetADetCRPsin+6,
    kDetAQxDetCQx    = kDetAQxDetBQx+6,
    kDetBQxDetCQx    = kDetAQxDetCQx+6,
    kDetAQxDetBQy    = kDetBQxDetCQx+6,
    kDetAQxDetCQy    = kDetAQxDetBQy+6,
    kDetBQxDetCQy    = kDetAQxDetCQy+6,
    kDetAQyDetBQx    = kDetBQxDetCQy+6,
    kDetAQyDetCQx    = kDetAQyDetBQx+6,
    kDetBQyDetCQx    = kDetAQyDetCQx+6,
    kDetAQyDetBQy    = kDetBQyDetCQx+6,
    kDetAQyDetCQy    = kDetAQyDetBQy+6,
    kDetBQyDetCQy    = kDetAQyDetCQy+6,

    kDetAQxDetBQxh    = kDetBQyDetCQy+6,
    kDetAQxhDetBQx    = kDetAQxDetBQxh+6,
    kDetAQxDetBQyh    = kDetAQxhDetBQx+6,
    kDetAQxhDetBQy    = kDetAQxDetBQyh+6,
    kDetAQyDetBQxh    = kDetAQxhDetBQy+6,
    kDetAQyhDetBQx    = kDetAQyDetBQxh+6,
    kDetAQyDetBQyh    = kDetAQyhDetBQx+6,
    kDetAQyhDetBQy    = kDetAQyDetBQyh+6,
                        
    kDetAQxDetCQxh    = kDetAQyhDetBQy+6,
    kDetAQxhDetCQx    = kDetAQxDetCQxh+6,
    kDetAQxDetCQyh    = kDetAQxhDetCQx+6,
    kDetAQxhDetCQy    = kDetAQxDetCQyh+6,
    kDetAQyDetCQxh    = kDetAQxhDetCQy+6,
    kDetAQyhDetCQx    = kDetAQyDetCQxh+6,
    kDetAQyDetCQyh    = kDetAQyhDetCQx+6,
    kDetAQyhDetCQy    = kDetAQyDetCQyh+6,

    kDetBQxDetCQxh    = kDetAQyhDetCQy+6,
    kDetBQxhDetCQx    = kDetBQxDetCQxh+6,
    kDetBQxDetCQyh    = kDetBQxhDetCQx+6,
    kDetBQxhDetCQy    = kDetBQxDetCQyh+6,
    kDetBQyDetCQxh    = kDetBQxhDetCQy+6,
    kDetBQyhDetCQx    = kDetBQyDetCQxh+6,
    kDetBQyDetCQyh    = kDetBQyhDetCQx+6,
    kDetBQyhDetCQy    = kDetBQyDetCQyh+6,

    kXaXbXc    = kDetBQyhDetCQy+6,
    kXaYbYc,
    kYaXbYc,
    kYaYbXc,
    kXbXcXa,
    kXbYcYa,
    kYbXcYa,
    kYbYcXa,
    kXaXcXb,
    kXaYcYb,
    kYaXcYb,
    kYaYcXb,
    kDetABCRPcos, 
    kDetABCRPsin,
    kDetACBRPcos,
    kDetACBRPsin,
    kDetBCARPcos,
    kDetBCARPsin,

    kNEventVars,           // number of event variables
    // Track variables -------------------------------------
    kPt,               //605
    kP,         //606
    kPx,        //607
    kPy,        //608
    kPz,        //609
    kCharge,
    kTheta,     //610
    kPhi,       //611
    kCosNPhi,   //612
    kSinNPhi = kCosNPhi+6,      //618
    kCos2NPhi = kSinNPhi+6,
    kSin2NPhi = kCos2NPhi+6,
    kEta = kSin2NPhi+6,          //624
    kRap,       //625
    kPtTPC,     //626
    kPhiTPC,    //627
    kEtaTPC,    //628
    kPin,       //629
    kTrackVZERORPdeltaPhi,      //630
    kTrackVZEROFlowVn = kTrackVZERORPdeltaPhi+6*3,          //631
    kTrackZDCRPdeltaPhi = kTrackVZEROFlowVn+6*3,      //630
    kTrackZDCFlowVn = kTrackZDCRPdeltaPhi+6*3,          //631
    kDcaXY = kTrackZDCFlowVn+6*3,     //643
    kDcaZ,              //644
    kITSncls,           //645
    kITSlayerHit,       //646
    kITSsignal,         //647
    kITSnSig,
    kTPCncls=kITSnSig+4,    //648
    kTPCclusBitFired,   //649
    kTPCNclusBitsFired, //650
    kTPCclustersPerBit, //651
    kTPCcrossedRows,    //652
    kTPCnclsIter1,      //653
    kTPCnclsF,          //654
    kTPCnclsRatio,       // TPCncls / TPCnclsF          //655
    kTPCnclsRatio2,      // TPCncls / TPCCrossedRows    //656
    kTPCsignal,         //657
    kTPCnSig,           //658
    kTPCchi2=kTPCnSig+4,           //658
    kTPCchi2Iter1,
    kTOFbeta,
    kTOFnSig,                   //663
    kTRDntracklets=kTOFnSig+4,  //667
    kTRDntrackletsPID,          //668
    kTRDpidProbabilities,       //669
    kEMCALmatchedEnergy=kTRDpidProbabilities+2,         //671
    kEMCALmatchedEOverP,        //672
    // Calorimeter cluster variables --------------------------------------
    kEMCALclusterEnergy,        //673
    kEMCALclusterDx,            //674
    kEMCALclusterDz,            //675
    kEMCALdetector,         // 0 - EMCAL; 1 - PHOS      //676
    // FMD channel variables
    kFMDId,
    kFMDMultiplicity,
    kFMDEqualizedMultiplicity,
    kFMDEta,
    kFMDPhi,
    kFMDdetector,
    // Tracking flags -----------------------------------------------------
    kTrackingFlag,              //677
    // Correlation variables ----------------------------------------------
    kDeltaPhi,          //678
    kDeltaTheta,        //679
    kDeltaEta,          //680
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


  //void GetUsedVars(Float_t* fUsedVars);

  // Function prototypes
  //static void FillEventInfo(EVENT* event, Float_t* values);
  //FRIEND static void FillEventInfo(TObject* event, Float_t* values);
  static void FillRPcorrelationInfo(Float_t* values, AliEventPlaneQvector* QvecA, AliEventPlaneQvector* QvecB);
  static void SetRPcorrelationVars(Float_t* values, Int_t ih, Bool_t enable);
  //static void FillEventOfflineTriggers(UShort_t triggerBit, Float_t* values);
  //static void FillEventOfflineTriggers(UShort_t triggerBit, AliReducedEvent*event, Float_t* values);
  //static void FillTrackingFlag(TRACK* track, UShort_t flag, Float_t* values);
  //static void FillTrackInfo(EVENT* event, TRACK* p, Float_t* values);
  //FRIEND static void FillTrackInfo(TObject* particle, Float_t* values);
  //static void FillITSlayerFlag(TRACK* track, Int_t layer, Float_t* values);
  //static void FillTPCclusterBitFlag(TRACK* track, Int_t bit, Float_t* values);
  //static void FillCaloClusterInfo(CLUSTER* cl, Float_t* values);
  //static void FillFMDInfo(FMD* fmd, Float_t* values);
  static Double_t DeltaPhi(Double_t phi1, Double_t phi2);  
  //static void GetThetaPhiCM(TRACK* leg1, TRACK* leg2,
  //Float_t &thetaHE, Float_t &phiHE, 
  //Float_t &thetaCS, Float_t &phiCS);
  //static void PrintTrackFlags(TRACK* track);
  static void PrintBits(ULong_t mask);
  static void SetDefaultVarNames();
  static void UnsetDefaultVarNames();
  static void UnsetUsedVars();
  static void SetUsedVar(Int_t var);
  static void UnsetUsedVar(Int_t var);
  static const Char_t* VarName(Int_t var) {SetDefaultVarNames(); const Char_t* ch = fVariableNames[var][0]; UnsetDefaultVarNames(); return ch;}


  static const Char_t* fTrackingFlagNames[kNTrackingFlags];
  static const Char_t* fOfflineTriggerNames[64];
  static const Char_t* fVariableNames[kNVars][2];
  static Bool_t fUsedVars[kNVars];
  static Bool_t fUseDefaultVariablesName;

private:

  AliEventPlaneVarManager(const AliEventPlaneVarManager &c);
  AliEventPlaneVarManager& operator= (const AliEventPlaneVarManager &c);

  ClassDef(AliEventPlaneVarManager, 1);

};  

////__________________________________________________________________
//inline void AliEventPlaneVarManager::FillEventInfo(EVENT* event, Float_t* values) {
//  //
//  // fill event wise info
//  //
//
//  const Double_t VZEROChannelRadii[64] = {6.0567, 6.0567, 6.0567, 6.0567, 6.0567, 6.0567, 6.0567, 6.0567,
//					    9.6977, 9.6977, 9.6977, 9.6977, 9.6977, 9.6977, 9.6977, 9.6977,
//					   15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504,
//					   26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031,
//					    5.9347, 5.9347, 5.9347, 5.9347, 5.9347, 5.9347, 5.9347, 5.9347,
//					   10.685, 10.685, 10.685, 10.685, 10.685, 10.685, 10.685, 10.685,
//					   18.116, 18.116, 18.116, 18.116, 18.116, 18.116, 18.116, 18.116,
//					   31.84, 31.84, 31.84, 31.84, 31.84, 31.84, 31.84, 31.84};
//
//
//  const Double_t VZEROAz = 340.0;   // cm
//  const Double_t VZEROCz = 90.0;    // cm
//
//  values[kRunNo]       = event->RunNo();
//  values[kBC]          = event->BC();
//#ifdef ALIREDUCEDEVENTV3_H
//  values[kTimeStamp]   = event->TimeStamp();
//  values[kEventType]   = event->EventType();
//#endif
//  values[kTriggerMask] = event->TriggerMask();
//  values[kIsPhysicsSelection]  = (event->IsPhysicsSelection() ? 1.0 : 0.0);
//  values[kNVtxContributors]    = event->VertexNContributors();
//  values[kNVtxTPCContributors] = event->VertexTPCContributors();
//  values[kVtxX]        = event->Vertex(0);
//  values[kVtxY]        = event->Vertex(1);
//  values[kVtxZ]        = event->Vertex(2);
//  values[kVtxXtpc]     = event->VertexTPC(0);
//  values[kVtxYtpc]     = event->VertexTPC(1);
//  values[kVtxZtpc]     = event->VertexTPC(2);
//  values[kDeltaVtxZ]   = values[kVtxZ] - values[kVtxZtpc];
//  values[kCentVZERO]   = event->CentralityVZERO();
//  values[kCentSPD]     = event->CentralitySPD();
//  values[kCentTPC]     = event->CentralityTPC();
//  values[kCentZDC]     = event->CentralityZEMvsZDC();
//  values[kCentQuality] = event->CentralityQuality();
//  values[kNV0total]        = event->NV0CandidatesTotal();
//  values[kNV0selected]     = event->NV0Candidates();
//  values[kNdielectrons]    = event->NDielectrons();
//  values[kNtracksTotal]    = event->NTracksTotal();
//  values[kNtracksSelected] = event->NTracks();
//  values[kSPDntracklets]   = event->SPDntracklets();
//  //values[kNFMD1channels]    =  event->NFMDchannels(0);
//  //values[kNFMD2Ichannels]   = event->NFMDchannels(1);
//  //values[kNFMD2Ochannels]   = event->NFMDchannels(2);
//  //values[kNFMD3Ichannels]   = event->NFMDchannels(3);
//  //values[kNFMD3Ochannels]   = event->NFMDchannels(4);
//  //values[kFMD1TotalMult]   =  event->FMDtotalMult(0);
//  //values[kFMD2ITotalMult]   = event->FMDtotalMult(1);
//  //values[kFMD2OTotalMult]   = event->FMDtotalMult(2);
//  //values[kFMD3ITotalMult]   = event->FMDtotalMult(3);
//  //values[kFMD3OTotalMult]   = event->FMDtotalMult(4);
//  values[kNPMDtracks]      = event->NPMDtracks();
//  
//  //VZERO detector information
//  values[kVZEROATotalMult] = event->MultVZEROA();
//  values[kVZEROCTotalMult] = event->MultVZEROC();
//  values[kVZEROTotalMult] = event->MultVZERO();
//  values[kVZEROAemptyChannels] = 0;
//  values[kVZEROCemptyChannels] = 0;
//  for(Int_t ich=0;ich<64;++ich) fUsedVars[kVZEROChannelMult+ich] = kTRUE; 
//  Float_t theta=0.0;
//  for(Int_t ich=0;ich<64;++ich) {
//      values[kVZEROChannelMult+ich] = event->MultChannelVZERO(ich);
//      if(values[kVZEROChannelMult+ich]<0.5) {
//        fUsedVars[kVZEROChannelMult+ich] = kFALSE;   // will not be filled in histograms by the histogram manager
//        if(ich<32) values[kVZEROCemptyChannels] += 1;
//        else values[kVZEROAemptyChannels] += 1;
//    }
//    if(fUsedVars[kVZEROChannelEta+ich]) {
//      if(ich<32) theta = TMath::ATan(VZEROChannelRadii[ich]/(VZEROCz-values[kVtxZ]));
//      else theta = TMath::Pi()-TMath::ATan(VZEROChannelRadii[ich]/(VZEROAz-values[kVtxZ]));
//      values[kVZEROChannelEta+ich] = -1.0*TMath::Log(TMath::Tan(theta/2.0));
//    }
//  }
//
//  
//  //ZDC detector information
//  values[kZDCATotalEnergy] = event->EnergyZDCA();
//  values[kZDCCTotalEnergy] = event->EnergyZDCC();
//  values[kZDCTotalEnergy] = event->EnergyZDC();
//  values[kZDCAemptyTowers] = 0;
//  values[kZDCCemptyTowers] = 0;
//  for(Int_t ich=0;ich<10;++ich) fUsedVars[kZDCTowerEnergy+ich] = kTRUE; 
//  for(Int_t ich=0;ich<10;++ich) {
//      values[kZDCTowerEnergy+ich] = event->EnergyZDCn(ich);
//      if(values[kZDCTowerEnergy+ich]<0.01||!(values[kZDCTowerEnergy+ich]<1e6)) {
//        fUsedVars[kZDCTowerEnergy+ich] = kFALSE;   // will not be filled in histograms by the histogram manager
//        if(ich<5&&ich!=0) values[kZDCCemptyTowers] += 1;
//        else if(ich!=5) values[kZDCAemptyTowers] += 1;
//      }
//    }
//    
//
//  
//  //TZERO information
//  values[kTZEROATotalMult] = event->AmplitudeTZEROA();
//  values[kTZEROCTotalMult] = event->AmplitudeTZEROC();
//  values[kTZEROTotalMult] = event->AmplitudeTZERO();
//
//      // TZERO channels multiplicities
//  values[kTZEROAemptyChannels] = 0;
//  values[kTZEROCemptyChannels] = 0;
//  for(Int_t ich=0;ich<24;++ich) fUsedVars[kTZEROChannelMult+ich] = kTRUE; 
//  theta=0.0;
//  for(Int_t ich=0;ich<24;++ich) {
//      values[kTZEROChannelMult+ich] = event->AmplitudeTZEROch(ich);
//      if(values[kTZEROChannelMult+ich]<1.e-2) {
//        fUsedVars[kTZEROChannelMult+ich] = kFALSE;   // will not be filled in histograms by the histogram manager
//        if(ich<12) values[kTZEROCemptyChannels] += 1;
//        else values[kTZEROAemptyChannels] += 1;
//      }
//    }
//
//    
//}


////__________________________________________________________________
//inline void AliEventPlaneVarManager::FillEventInfo(AliVEvent* event, Float_t* values) {
//  //
//  // fill event wise info
//  //
//  values[kRunNo]       = event->GetRunNumber();
//  //values[kTriggerMask] = event->TriggerMask();
//  //values[kIsPhysicsSelection]  = (event->IsPhysicsSelection() ? 1.0 : 0.0);
//  //values[kNVtxTPCContributors] = event->VertexTPCContributors();
//  values[kVtxX]        = -999.;
//  values[kVtxY]        = -999.;
//  values[kVtxZ]        = -999.;
//  const AliVVertex *primVtx = event->GetPrimaryVertex();
//  if (primVtx){
//    values[kVtxX]        = primVtx->GetX();
//    values[kVtxY]        = primVtx->GetY();
//    values[kVtxZ]        = primVtx->GetZ();
//    values[kNVtxContributors]    = primVtx->GetNContributors();
//  }
//
//  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(event);
//  AliCentrality* cent = esdEvent->GetCentrality();
//  if(cent){
//    values[kCentVZERO]   = cent->GetCentralityPercentile("V0M");
//    values[kCentSPD]     = cent->GetCentralityPercentile("CL1");
//    values[kCentTPC]     = cent->GetCentralityPercentile("TRK");
//    values[kCentQuality] = cent->GetQuality();
//  }
//    
//}




inline void AliEventPlaneVarManager::FillRPcorrelationInfo(Float_t* values, AliEventPlaneQvector* QvecA, AliEventPlaneQvector* QvecB) {

      values[kDetAmult] = QvecA->Multiplicity();
      values[kDetBmult] = QvecB->Multiplicity();

    for(Int_t ih=0; ih<6; ++ih) {


      // Macro defined detector event plane
      values[kDetARP   +ih] = QvecA->EventPlane(ih+1);
      values[kDetBRP   +ih] = QvecB->EventPlane(ih+1);
      values[kCosDetARP   +ih] = TMath::Cos(Double_t(ih+1)*QvecA->EventPlane(ih+1));
      values[kCosDetBRP   +ih] = TMath::Cos(Double_t(ih+1)*QvecB->EventPlane(ih+1));
      values[kSinDetARP   +ih] = TMath::Sin(Double_t(ih+1)*QvecA->EventPlane(ih+1));
      values[kSinDetBRP   +ih] = TMath::Sin(Double_t(ih+1)*QvecB->EventPlane(ih+1));
      values[kDetAQvecX+ih] = QvecA->Qx(ih+1);
      values[kDetAQvecY+ih] = QvecA->Qy(ih+1);
      values[kDetBQvecX+ih] = QvecB->Qx(ih+1);
      values[kDetBQvecY+ih] = QvecB->Qy(ih+1);

      // Qx,Qy correlations for macro detectors

      values[kDetAQxDetBQx+ih] = values[kDetAQvecX+ih]*values[kDetBQvecX+ih];
      values[kDetAQxDetBQy+ih] = values[kDetAQvecX+ih]*values[kDetBQvecY+ih];
      values[kDetAQyDetBQx+ih] = values[kDetAQvecY+ih]*values[kDetBQvecX+ih];
      values[kDetAQyDetBQy+ih] = values[kDetAQvecY+ih]*values[kDetBQvecY+ih];

      values[kDetAQxhDetBQx+ih] = values[kDetAQvecX+ih+1]*values[kDetBQvecX+ih];
      values[kDetAQxDetBQxh+ih] = values[kDetAQvecX+ih]*values[kDetBQvecX+ih+1];
      values[kDetAQxhDetBQy+ih] = values[kDetAQvecX+ih+1]*values[kDetBQvecY+ih];
      values[kDetAQxDetBQyh+ih] = values[kDetAQvecX+ih]*values[kDetBQvecY+ih+1];
      values[kDetAQyhDetBQx+ih] = values[kDetAQvecY+ih+1]*values[kDetBQvecX+ih];
      values[kDetAQyDetBQxh+ih] = values[kDetAQvecY+ih]*values[kDetBQvecX+ih+1];
      values[kDetAQyhDetBQy+ih] = values[kDetAQvecY+ih+1]*values[kDetBQvecY+ih];
      values[kDetAQyDetBQyh+ih] = values[kDetAQvecY+ih]*values[kDetBQvecY+ih+1];


      // cos/sin (n*(psi_A-psi_B))
      values[kDetADetBRPres+ih] = DeltaPhi(QvecA->EventPlane(ih+1), QvecB->EventPlane(ih+1));
      values[kDetADetBRPsin+ih]  = TMath::Sin(Double_t(ih+1)*values[kDetADetBRPres+ih]);
      values[kDetADetBRPres+ih]  = TMath::Cos(Double_t(ih+1)*values[kDetADetBRPres+ih]);


      values[kDetADetBRPcosplus+ih] = QvecA->EventPlane(ih+1)+QvecB->EventPlane(ih+1);
      values[kDetADetBRPcosplus+ih]  = TMath::Cos(Double_t(ih+1)*values[kDetADetBRPcosplus+ih]);

      values[kDetADetBRPsinplus+ih] = QvecA->EventPlane(ih+1)+QvecB->EventPlane(ih+1);
      values[kDetADetBRPsinplus+ih]  = TMath::Sin(Double_t(ih+1)*values[kDetADetBRPsinplus+ih]);

    }


}



inline void AliEventPlaneVarManager::SetRPcorrelationVars(Float_t* values, Int_t ih, Bool_t enable) {

      fUsedVars[kDetARP+ih] = enable;
      fUsedVars[kDetBRP+ih] = enable;
      fUsedVars[kCosDetARP+ih] = enable;
      fUsedVars[kCosDetBRP+ih] = enable;
      fUsedVars[kSinDetARP+ih] = enable;
      fUsedVars[kSinDetBRP+ih] = enable;
      fUsedVars[kDetAQvecX+ih] = enable;
      fUsedVars[kDetAQvecY+ih] = enable;
      fUsedVars[kDetBQvecX+ih] = enable;
      fUsedVars[kDetBQvecY+ih] = enable;
      fUsedVars[kDetAQxDetBQx+ih] = enable;
      fUsedVars[kDetAQxDetBQy+ih] = enable;
      fUsedVars[kDetAQyDetBQx+ih] = enable;
      fUsedVars[kDetAQyDetBQy+ih] = enable;
      fUsedVars[kDetAQxhDetBQx+ih] = enable;
      fUsedVars[kDetAQxDetBQxh+ih] = enable;
      fUsedVars[kDetAQxhDetBQy+ih] = enable;
      fUsedVars[kDetAQxDetBQyh+ih] = enable;
      fUsedVars[kDetAQyhDetBQx+ih] = enable;
      fUsedVars[kDetAQyDetBQxh+ih] = enable;
      fUsedVars[kDetAQyhDetBQy+ih] = enable;
      fUsedVars[kDetAQyDetBQyh+ih] = enable;
      fUsedVars[kDetADetBRPres+ih]= enable;
      fUsedVars[kDetADetBRPsin+ih]= enable;
      fUsedVars[kDetADetBRPres+ih]= enable;
      fUsedVars[kDetADetBRPcosplus+ih]= enable;
      fUsedVars[kDetADetBRPcosplus+ih]= enable;
      fUsedVars[kDetADetBRPsinplus+ih]= enable;
      fUsedVars[kDetADetBRPsinplus+ih]= enable;


}






////_________________________________________________________________
//inline void AliEventPlaneVarManager::FillITSlayerFlag(TRACK* track, Int_t layer, Float_t* values) {
//  //
//  // fill the ITS layer hit
//  //
//  values[kITSlayerHit] = -1.0*(layer+1);
//  if(fUsedVars[kITSlayerHit] && track->ITSLayerHit(layer)) values[kITSlayerHit] = layer+1;
//}
//
//
////_________________________________________________________________
//inline void AliEventPlaneVarManager::FillTPCclusterBitFlag(TRACK* track, Int_t bit, Float_t* values) {
//  //
//  // fill the TPC cluster map
//  //
//  values[kTPCclusBitFired] = -1;
//  if(fUsedVars[kTPCclusBitFired] && track->TPCClusterMapBitFired(bit)) values[kTPCclusBitFired] = bit;
//}
//
//
////_________________________________________________________________
//inline void AliEventPlaneVarManager::FillTrackingFlag(TRACK* track, UShort_t flag, Float_t* values) {
//  //
//  // fill the tracking flag
//  //
//  values[kTrackingFlag] = -1;
//  if(track->CheckTrackStatus(flag)) values[kTrackingFlag] = flag;
//}


////_________________________________________________________________
//inline void AliEventPlaneVarManager::FillEventOfflineTriggers(UShort_t triggerBit, AliReducedEvent*event, Float_t* values) {
//  //
//  // fill the trigger bit input
//  //
//  if(triggerBit>=64) return;
//  ULong64_t trigger = BIT(0);
//  values[kOfflineTrigger] = triggerBit;
//  values[kOfflineTriggerFired] = (event->TriggerMask()&(trigger<<triggerBit) ? 1.0 : 0.0);
//  values[kOfflineTriggerFired2] = (values[kOfflineTriggerFired]>0.01 ? triggerBit : -1.0); 
//}


////_________________________________________________________________
//inline void AliEventPlaneVarManager::FillTrackInfo(EVENT* event, TRACK* p, Float_t* values) {
//  //
//  // fill track information
//  //
//                        values[kPt]     = p->Pt();
//                        values[kPtTPC]  = p->PtTPC();
//  if(fUsedVars[kP])     values[kP]      = p->P();
//  if(fUsedVars[kPx])    values[kPx]      = p->Px();
//  if(fUsedVars[kPy])    values[kPy]      = p->Py();
//  if(fUsedVars[kPz])    values[kPz]      = p->Pz();
//  values[kCharge]  = p->Charge();
//  //if(fUsedVars[kCharge]) values[kCharge]  = p->Charge();
//  if(fUsedVars[kTheta]) values[kTheta]  = p->Theta();
//                        values[kPhi]    = p->Phi();
//  for(Int_t ih=1; ih<=6; ++ih) {
//    if(fUsedVars[kCosNPhi+ih-1]) values[kCosNPhi+ih-1] = TMath::Cos(p->Phi()*ih);
//    if(fUsedVars[kSinNPhi+ih-1]) values[kSinNPhi+ih-1] = TMath::Sin(p->Phi()*ih);
//  }
//  for(Int_t ih=1; ih<=6; ++ih) {
//    if(fUsedVars[kCos2NPhi+ih-1]) values[kCos2NPhi+ih-1] = TMath::Cos(p->Phi()*ih*2.);
//    if(fUsedVars[kSin2NPhi+ih-1]) values[kSin2NPhi+ih-1] = TMath::Sin(p->Phi()*ih*2.);
//  }
//                        values[kPhiTPC] = p->PhiTPC();
//                        values[kEta]    = p->Eta();
//                        values[kEtaTPC] = p->EtaTPC();
//                        values[kPin]    = p->Pin();
//                        values[kDcaXY]  = p->DCAxy();
//                        values[kDcaZ]   = p->DCAz();
//  
//  if(fUsedVars[kITSncls]) values[kITSncls] = p->ITSncls();
//                          values[kITSsignal] = p->ITSsignal();
//  
//  
//  values[kTPCncls] = p->TPCncls();
//  if(fUsedVars[kTPCnclsRatio]) 
//    values[kTPCnclsRatio] = (p->TPCFindableNcls()>0 ? Float_t(p->TPCncls())/Float_t(p->TPCFindableNcls()) : 0.0);
//  if(fUsedVars[kTPCnclsRatio2]) 
//    values[kTPCnclsRatio2] = (p->TPCCrossedRows()>0 ? Float_t(p->TPCncls())/Float_t(p->TPCCrossedRows()) : 0.0);
//  values[kTPCnclsIter1]   = p->TPCnclsIter1();
//  values[kTPCnclsF]       = p->TPCFindableNcls();
//  values[kTPCcrossedRows] = p->TPCCrossedRows();
//  values[kTPCsignal]      = p->TPCsignal();
//  if(fUsedVars[kTPCNclusBitsFired]) values[kTPCNclusBitsFired] = p->TPCClusterMapBitsFired();
//  if(fUsedVars[kTPCclustersPerBit]) {
//    Int_t nbits = p->TPCClusterMapBitsFired();
//    values[kTPCclustersPerBit] = (nbits>0 ? values[kTPCncls]/Float_t(nbits) : 0.0);
//  }
//  
//  values[kTOFbeta] = p->TOFbeta();
//  for(Int_t specie=kElectron; specie<=kProton; ++specie) {
//    values[kITSnSig+specie] = p->ITSnSig(specie);
//    values[kTPCnSig+specie] = p->TPCnSig(specie);
//    values[kTOFnSig+specie] = p->TOFnSig(specie);
//  }
//  values[kTRDpidProbabilities]   = p->TRDpid(0);
//  values[kTRDpidProbabilities+1] = p->TRDpid(1);
//  
//  values[kTRDntracklets]    = p->TRDntracklets(0);
//  values[kTRDntrackletsPID] = p->TRDntracklets(1);
//  
//  if(fUsedVars[kEMCALmatchedEnergy] || fUsedVars[kEMCALmatchedEOverP]) {
//    CLUSTER* cluster = event->GetCaloCluster(p->CaloClusterId());
//    values[kEMCALmatchedEnergy] = (cluster ? cluster->Energy() : -999.0);
//    Float_t mom = p->P();
//    values[kEMCALmatchedEnergy] = (TMath::Abs(mom)>1.e-8 && cluster ? values[kEMCALmatchedEOverP]/mom : -999.0);
//  }
//  
//  // Fill track flow variables
//  for(Int_t iVZEROside=0; iVZEROside<3; ++iVZEROside) {
//    for(Int_t ih=0; ih<6; ++ih) {
//      if(fUsedVars[kTrackVZEROFlowVn+iVZEROside*6+ih]) {
//        //values[kTrackVZEROFlowVn+iVZEROside*6+ih]         = TMath::Cos(DeltaPhi(values[kPhi],values[kVZERORP+iVZEROside*6+ih])*(ih+1));
//        values[kTrackVZEROFlowVn+iVZEROside*6+ih] = TMath::Cos((values[kPhi]-values[kVZERORP+iVZEROside*6+ih])*(ih+1));        
//      }
//      if(fUsedVars[kTrackVZERORPdeltaPhi+iVZEROside*6+ih]) {
//        //values[kTrackVZERORPdeltaPhi+iVZEROside*6+ih] = values[kPhi]/Float_t(ih+1) - values[kVZERORP+iVZEROside*6+ih];
//        values[kTrackVZERORPdeltaPhi+iVZEROside*6+ih] = TMath::ACos(TMath::Cos((values[kPhi]-values[kVZERORP+iVZEROside*6+ih])*(ih+1)))/Float_t(ih+1);
//      }
//    }
//  }
//    for(Int_t iZDCside=0; iZDCside<3; ++iZDCside) {
//    for(Int_t ih=0; ih<6; ++ih) {
//      if(fUsedVars[kTrackZDCFlowVn+iZDCside*6+ih]) {
//        //values[kTrackZDCFlowVn+iZDCside*6+ih]         = TMath::Cos(DeltaPhi(values[kPhi],values[kZDCRP+iZDCside*6+ih])*(ih+1));
//        values[kTrackZDCFlowVn+iZDCside*6+ih] = TMath::Cos((values[kPhi]-values[kZDCRP+iZDCside*6+ih])*(ih+1));        
//      }
//      if(fUsedVars[kTrackZDCRPdeltaPhi+iZDCside*6+ih]) {
//        //values[kTrackZDCRPdeltaPhi+iZDCside*6+ih] = values[kPhi]/Float_t(ih+1) - values[kZDCRP+iZDCside*6+ih];
//        values[kTrackZDCRPdeltaPhi+iZDCside*6+ih] = TMath::ACos(TMath::Cos((values[kPhi]-values[kZDCRP+iZDCside*6+ih])*(ih+1)))/Float_t(ih+1);
//      }
//    }
//  }
//}


////__________________________________________________________________
//inline void AliEventPlaneVarManager::FillTrackInfo(AliESDtrack* particle, Float_t* values) {
//
//  Float_t dcaxy=0.0;
//  Float_t dcaz=0.0;
//  particle->GetImpactParameters(dcaxy,dcaz);
//
//  values[kPx]        = particle->Px();
//  values[kPy]        = particle->Py();
//  values[kPz]        = particle->Pz();
//  values[kPt]        = particle->Pt();
//  values[kP]         = particle->P();
//  values[kPhi]       = particle->Phi();
//  values[kTheta]     = particle->Theta();
//  values[kEta]       = particle->Eta();
//  values[kCharge]    = particle->Charge();
//  values[kDcaXY]     = dcaxy;
//  values[kDcaZ]      = dcaz;
//
//  values[kITSncls]       = particle->GetNcls(0); 
//  values[kTPCncls]       = particle->GetTPCNcls();
//  values[kTPCnclsIter1]  = particle->GetTPCNclsIter1();
//  values[kTPCchi2]       = particle->GetTPCchi2()/values[kTPCncls];
//  values[kTPCchi2Iter1]  = particle->GetTPCchi2Iter1()/values[kTPCnclsIter1];
//  values[kTPCsignal]     = particle->GetTPCsignal();
//  ////values[kNclsTPCiter1]  = particle->GetTPCNclsIter1(); // TODO: get rid of the plain numbers
//  //values[kNFclsTPC]      = particle->GetTPCNclsF();
//  //values[kNFclsTPCr]     = particle->GetTPCClusterInfo(2,1);
//  //values[kNFclsTPCrFrac] = particle->GetTPCClusterInfo(2);
//  //values[kTPCsignalNfrac]= tpcNcls>0?tpcSignalN/tpcNcls:0;
//  //values[kNclsTRD]       = particle->GetNcls(2); // TODO: get rid of the plain numbers
//  //values[kTRDntracklets] = particle->GetTRDntracklets(); // TODO: GetTRDtracklets/GetTRDntracklets?
//  //values[kTRDpidQuality] = particle->GetTRDpidQuality();
//  //values[kTrackStatus]   = (Double_t)particle->GetStatus();
//  //if (tpcNcls>0) values[AliDielectronVarContainer::kTPCchi2Cl] = particle->GetTPCchi2() / tpcNcls;
//
//}



////_________________________________________________________________
//inline void AliEventPlaneVarManager::FillCaloClusterInfo(CLUSTER* cl, Float_t* values) {
//  //
//  // Fill calorimeter cluster information
//  //
//  values[kEMCALclusterEnergy] = cl->Energy();
//  values[kEMCALclusterDx] = cl->Dx();
//  values[kEMCALclusterDz] = cl->Dz();
//  values[kEMCALdetector] = (cl->IsEMCAL() ? CLUSTER::kEMCAL : CLUSTER::kPHOS);
//}

//__________________________________________________________________
inline void AliEventPlaneVarManager::SetUsedVar(Int_t var) {
  (var<0||var>=kNVars ? fUsedVars[var] = kFALSE : fUsedVars[var] = kTRUE);
  return;
}



//__________________________________________________________________
inline void AliEventPlaneVarManager::UnsetUsedVar(Int_t var) {
//  var<0||var>=kNVars ? fUsedVars[var] = kFALSE : 0;
//  return;
}

#endif    
