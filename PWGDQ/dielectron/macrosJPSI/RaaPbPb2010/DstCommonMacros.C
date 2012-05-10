//
// Common macros used to analyze dst trees
// Author: Ionut-Cristian Arsene, 2012/03/04
// email: i.c.arsene@gsi.de
//
#include <iostream>
#include <fstream>
using namespace std;

#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TFile.h>
#include <TDirectory.h>
#include <THashList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TIterator.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TAxis.h>

#ifndef ALICORRELATIONREDUCEDEVENT_H
#include "AliCorrelationReducedEvent.h"
#endif

namespace DstCommonMacros {

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
    kRunNo = 0,         // run number
    kBC,                // bunch crossing
    kTriggerMask,       // trigger mask
    kOfflineTrigger,    // offline trigger
    kOfflineTriggerFired,  // offline trigger fired
    kOfflineTriggerFired2,  // offline trigger if fired, -1 if not fired
    kIsPhysicsSelection,    // physics selection
    kVtxX,              // vtx X
    kVtxY,              // vtx Y
    kVtxZ,              // vtx Z
    kVtxXtpc,           // vtx X from tpc
    kVtxYtpc,           // vtx Y from tpc
    kVtxZtpc,           // vtx Z from tpc
    kCentVZERO,         // centrality from VZERO
    kCentSPD,           // centrality from SPD
    kCentTPC,           // centrality from TPC
    kCentZDC,           // centrality from ZDC
    kCentQuality,       // centrality quality
    kNV0total,          // total number of V0s in the esd
    kNV0selected,       // number of V0s selected
    kNdielectrons,      // number of dielectron pairs
    kNpairsSelected,    // number of selected dielectron pairs per event
    kNtracksTotal,      // total number of tracks
    kNtracksSelected,   // number of selected tracks
    kSPDntracklets,     // SPD number of tracklets in |eta|<1.0
    kEventMixingId,     // Id of the event mixing category
    // VZERO event plane related variables
    kVZEROAemptyChannels,  // Number of empty VZERO channels in A side
    kVZEROCemptyChannels,  // Number of empty VZERO channels in C side
    kVZEROChannelMult,                        // VZERO multiplicity per channel
    kVZEROChannelEta = kVZEROChannelMult+64,  // pseudo-rapidity of a VZERO channel
    kVZEROQvecX      = kVZEROChannelEta+64,   // Q-vector components for harmonics 1-6 and 
    kVZEROQvecY      = kVZEROQvecX+6*3,        //  6- n-harmonics; 3- A,C and A&C options
    kVZERORP         = kVZEROQvecY+6*3,           // VZERO reaction plane from A,C and A&C sides (harmonics 1-6)
    kVZERORPres      = kVZERORP+6*3,     // VZERO reaction plane resolution (sqrt(n*(RPa-RPc)))
    kVZEROXaXc       = kVZERORPres+6,           // correlations for the components of the Q vector
    kVZEROXaYa       = kVZEROXaXc+6,
    kVZEROXaYc       = kVZEROXaYa+6,
    kVZEROYaXc       = kVZEROXaYc+6,
    kVZEROYaYc       = kVZEROYaXc+6,
    kVZEROXcYc       = kVZEROYaYc+6,
    kVZEROdeltaRPac  = kVZEROXcYc+6,         // Psi_VZEROA-Psi_VZEROC
    kVZEROdeltaRPa   = kVZEROdeltaRPac+6,      // Psi_VZEROA6-Psi_VZEROA5, Psi_VZEROA7-Psi_VZEROA5, Psi_VZEROA6-Psi_VZEROA5 
    kVZEROdeltaRPc   = kVZEROdeltaRPa+6*2,     // Psi_VZEROA2-Psi_VZEROA1, Psi_VZEROA3-Psi_VZEROA1, Psi_VZEROA4-Psi_VZEROA1
    kVZEROflowV2TPC  = kVZEROdeltaRPc+6*2,     // vzero v2 using TPC event plane
    // TPC event plane variables
    kTPCQvecX = kVZEROflowV2TPC+64,   // TPC Q-vector components for harmonics 1-6
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
    kNEventVars       = kRPdeltaVZEROCtpc+6,           // number of event variables
    // Pair variables --------------------------------------
    kMass=kNEventVars,
    kCandidateId,
    kPairType,                          // 0 ++; 1 +-; 2 --
    kPairPt,
    kPairPx,
    kPairPy,
    kPairPz,
    kPairP,
    kPairRap,
    kPairEta,
    kPairTheta,
    kPairPhi,
    kPairLxy,
    kPairOpeningAngle,
    kPairCosNPhi,                                                // cos (n*phi)
    kPairSinNPhi             = kPairCosNPhi+6,                   // sin (n*phi)
    kPairDeltaPhiVZEROFlowVn = kPairSinNPhi+6,                   // phi - Psi_{VZERO}
    kPairDeltaPhiTPCFlowVn   = kPairDeltaPhiVZEROFlowVn+6*3,     // phi - Psi_{TPC}
    kPairVZEROFlowVn         = kPairDeltaPhiTPCFlowVn+6,         // vn{Psi_{n,VZERO}}
    kPairTPCFlowVn           = kPairVZEROFlowVn+6*3,             // vn{Psi_{n,TPC}}
    kPairVZEROFlowV3Psi2     = kPairTPCFlowVn+6,                 // v3{Psi_{2,VZERO}}
    kPairTPCFlowV3Psi2       = kPairVZEROFlowV3Psi2+3,           // v3{Psi_{2,TPC}}
    // Track variables -------------------------------------
    kPt=kPairTPCFlowV3Psi2+1,
    kP,
    kTheta,
    kPhi,
    kEta,
    kRap,
    kPtTPC,
    kPhiTPC,
    kEtaTPC,
    kPin,
    kTrackDeltaPhiVZEROFlowVn,
    kTrackVZEROFlowVn = kTrackDeltaPhiVZEROFlowVn+6*2,
    kDcaXY            = kTrackVZEROFlowVn+6*2,
    kDcaZ,
    kITSncls,
    kITSsignal,
    kTPCncls,
    kTPCcrossedRows,
    kTPCnclsIter1,
    kTPCnclsF,
    kTPCnclsRatio,
    kTPCsignal,
    kTPCnSig,
    kTOFbeta=kTPCnSig+4,
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
    // Tracking flags -----------------------------------------------------
    kTrackingFlag,
    // Correlation variables ----------------------------------------------
    kDeltaPhi,
    kDeltaTheta,
    kDeltaEta,
    kNVars
  };


  Bool_t gUsedVars[kNVars] = {kFALSE};

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

  const Char_t* gkTrackingFlagNames[kNTrackingFlags] = {
    "kITSin", "kITSout", "kITSrefit", "kITSpid",
    "kTPCin", "kTPCout", "kTPCrefit", "kTPCpid",
    "kTRDin", "kTRDout", "kTRDrefit", "kTRDpid",
    "kTOFin", "kTOFout", "kTOFrefit", "kTOFpid", "kTOFmismatch",
    "kHMPIDout", "kHMPIDpid", 
    "kEMCALmatch", "kPHOSmatch", 
    "kTRDbackup", "kTRDStop",
    "kESDpid", "kTIME", "kGlobalMerge",
    "kITSpureSA", 
    "kMultInV0",
    "kMultSec",
    "kTRDnPlanes",
    "kEMCALNoMatch"
  };

  // offline triggers as defined in AliVEvent.h
  // NOTE: Check consistency with updates in aliroot!!!
  enum EOfflineTriggerTypes { 
    kMB           = BIT(0), // Minimum bias trigger, i.e. interaction trigger, offline SPD or V0 selection
    kINT7         = BIT(1), // V0AND trigger, offline V0 selection
    kMUON         = BIT(2), // Muon trigger, offline SPD or V0 selection
    kHighMult     = BIT(3), // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
    kEMC1         = BIT(4), // EMCAL trigger
    kCINT5        = BIT(5), // Minimum bias trigger without SPD. i.e. interaction trigger, offline V0 selection
    kCMUS5        = BIT(6), // Muon trigger, offline V0 selection
    kMUSPB        = BIT(6), // idem for PbPb
    kMUSH7        = BIT(7), // Muon trigger: high pt, single muon, offline V0 selection, CINT7 suite
    kMUSHPB       = BIT(7), // idem for PbPb
    kMUL7         = BIT(8), // Muon trigger: like sign dimuon, offline V0 selection, CINT7 suite
    kMuonLikePB   = BIT(8), // idem for PbPb
    kMUU7         = BIT(9), // Muon trigger, unlike sign dimuon, offline V0 selection, CINT7 suite
    kMuonUnlikePB = BIT(9), // idem for PbPb
    kEMC7         = BIT(10), // EMCAL trigger, CINT7 suite
    kMUS7         = BIT(11), // Muon trigger: low pt, single muon, offline V0 selection, CINT7 suite
    kPHI1         = BIT(12), // PHOS trigger, CINT1 suite
    kPHI7         = BIT(13), // PHOS trigger, CINT7 suite
    kPHOSPb       = BIT(13), // idem for PbPb
    kEMCEJE       = BIT(14), // EMCAL jet patch trigger
    kEMCEGA       = BIT(15), // EMCAL gamma trigger
    kCentral      = BIT(16), // PbPb central collision trigger
    kSemiCentral  = BIT(17), // PbPb semicentral collision trigger
    kDG5          = BIT(18), // Double gap diffractive
    kZED          = BIT(19), // ZDC electromagnetic dissociation
    kUserDefined  = BIT(27), // Set when custom trigger classes are set in AliPhysicsSelection, offline SPD or V0 selection
    // Bits 28 and above are reserved for FLAGS
    kFastOnly     = BIT(30), // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB
    kAny          = 0xffffffff, // to accept any trigger
    kAnyINT       = kMB | kINT7 | kCINT5 // to accept any interaction (aka minimum bias) trigger
  };

  const Char_t* gkOfflineTriggerNames[64] = {
    "MB",              "INT7",              "MUON", "HighMult",    "EMC1", "CINT5",       "CMUS5/MUSPB", "MUSH7/MUSHPB",
    "MUL7/MuonLikePB", "MUU7/MuonUnlikePB", "EMC7", "MUS7",        "PHI1", "PHI7/PHOSPb", "EMCEJE",      "EMCEGA",
    "Central",         "SemiCentral",       "DG5",  "ZED",         "N/A",  "N/A",         "N/A",         "N/A",  
    "N/A",             "N/A",               "N/A",  "UserDefined", "N/A",  "N/A",         "FastOnly",    "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",         "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",         "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",         "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",         "N/A"
  };

  enum ITSLayerMap {
    kITSfirst  =  1,
    kITSsecond =  2,
    kITSthird  =  4,
    kITSfourth =  8,
    kITSfifth  = 16,
    kITSsixth  = 32
  };

  // radii of VZERO channels centers (in cm)
  const Double_t gkVZEROChannelRadii[64] = {6.0567, 6.0567, 6.0567, 6.0567, 6.0567, 6.0567, 6.0567, 6.0567,
					    9.6977, 9.6977, 9.6977, 9.6977, 9.6977, 9.6977, 9.6977, 9.6977,
					   15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504,
					   26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031,
					    5.9347, 5.9347, 5.9347, 5.9347, 5.9347, 5.9347, 5.9347, 5.9347,
					   10.685, 10.685, 10.685, 10.685, 10.685, 10.685, 10.685, 10.685,
					   18.116, 18.116, 18.116, 18.116, 18.116, 18.116, 18.116, 18.116,
					   31.84, 31.84, 31.84, 31.84, 31.84, 31.84, 31.84, 31.84};
  const Double_t gkVZEROAz = 340.0;   // cm
  const Double_t gkVZEROCz = 90.0;    // cm
  const Double_t gkVZEROminMult = 0.5;   // minimum VZERO channel multiplicity

  // Pointer to the current event
  AliCorrelationReducedEvent*       gCurrentEvent       = 0x0;
  AliCorrelationReducedEventFriend* gCurrentEventFriend = 0x0;

  // Event mixing variables
  Int_t gEMCategories=0;
  TString* gEMCategoryNames;

  TObjArray* gHistLists=0x0;   // main histogram list for the current running process
  TDirectoryFile* gHistListsOld=0x0;  // main directory for a standard tree analysis output (used for calibration, plotting etc.)
  TFile* gHistFile = 0x0;      // pointer to a TFile opened for reading

  TProfile2D* gVzeroAvMult[64] = {0x0};   // pointer to the array of average VZERO multiplicity per channel
  TProfile2D* gQvecCentering[AliCorrelationReducedEventFriend::kNdetectors][fgkNMaxHarmonics][2] = {{{0x0}}};  // pointer to the array of Qvec centering histograms

  // pt range for J/psi's
  const Float_t gkJpsiPtCut[2] = {0.0, 20.0};

  // Function prototypes
  void WriteOutput(TFile* saveFile);
  //void DefineHistograms(const Char_t* histClasses);
  TChain* GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries, TChain* friendChain=0x0, const Char_t* friendChainFile=0x0);
  void FillEventInfo(AliCorrelationReducedEvent* event, Float_t* values, AliCorrelationReducedEventFriend* eventF=0x0);
  void FillEventOfflineTriggers(UShort_t triggerBit, Float_t* values);
  void FillTrackingFlag(AliCorrelationReducedTrack* track, UShort_t flag, Float_t* values);
  void FillTrackInfo(AliCorrelationReducedTrack* p, Float_t* values);
  void FillPairInfo(AliCorrelationReducedPair* p, Float_t* values);
  void FillPairInfo(AliCorrelationReducedTrack* t1, AliCorrelationReducedTrack* t2, Int_t type, Float_t* values);
  void FillCorrelationInfo(AliCorrelationReducedPair* p, AliCorrelationReducedTrack* t, Float_t* values);
  void FillCaloClusterInfo(AliCorrelationReducedCaloCluster* cl, Float_t* values);
  void FillHistClass(const Char_t* className, Float_t* values);
  void DoEventMixing(TList* list1, TList* list2, Float_t* values, Int_t mixingType, const Char_t* histClass);
  void EventMixingPairTracks(TList* pairs, TList* tracks, Float_t* values);
  void EventMixingResonanceLegs(TList* posLegs, TList* negLegs, Float_t* values, Int_t type, const Char_t* histClass);
  Double_t DeltaPhi(Double_t phi1, Double_t phi2);  // calculate delta phi in the (-pi,+pi) interval
  Bool_t IsPairSelectedEM(Float_t* values);        // pair selection used in the mixed event
  void AddHistClass(const Char_t* histClass);
  Int_t ValidateHistogramName(THashList* hList, const Char_t* name);
  void AddHistogram(const Char_t* histClass,
		    const Char_t* name, const Char_t* title, Bool_t isProfile,
                    Int_t nXbins, Double_t xmin, Double_t xmax, Int_t varX,
		    Int_t nYbins=0, Double_t ymin=0, Double_t ymax=0, Int_t varY=kNothing,
		    Int_t nZbins=0, Double_t zmin=0, Double_t zmax=0, Int_t varZ=kNothing,
		    const Char_t* xLabels="", const Char_t* yLabels="", const Char_t* zLabels="");
  void AddHistogram(const Char_t* histClass,
		    const Char_t* name, const Char_t* title, Bool_t isProfile,
                    Int_t nXbins, Double_t* xbins, Int_t varX,
		    Int_t nYbins=0, Double_t* ybins=0x0, Int_t varY=kNothing,
		    Int_t nZbins=0, Double_t* zbins=0x0, Int_t varZ=kNothing,
  		    const Char_t* xLabels="", const Char_t* yLabels="", const Char_t* zLabels="");
  void MakeAxisLabels(TAxis* ax, const Char_t* labels);
  void InitFile(const Char_t* filename);    // open an old output filename
  void CloseFile();
  TObject* GetHistogram(const Char_t* listname, const Char_t* hname);  // get a histogram from an old output

}  // end declarations for namespace DstCommonMacros

//__________________________________________________________________
void DstCommonMacros::FillEventInfo(AliCorrelationReducedEvent* event, Float_t* values, AliCorrelationReducedEventFriend* eventF/*=0x0*/) {
  //
  // fill event wise info
  //
  values[kRunNo]       = event->RunNo();
  values[kBC]          = event->BC();
  values[kTriggerMask] = event->TriggerMask();
  values[kIsPhysicsSelection] = (event->IsPhysicsSelection() ? 1.0 : 0.0);
  values[kVtxX]        = event->Vertex(0);
  values[kVtxY]        = event->Vertex(1);
  values[kVtxZ]        = event->Vertex(2);
  values[kVtxXtpc]     = event->VertexTPC(0);
  values[kVtxYtpc]     = event->VertexTPC(1);
  values[kVtxZtpc]     = event->VertexTPC(2);
  values[kCentVZERO]   = event->CentralityVZERO();
  values[kCentSPD]     = event->CentralitySPD();
  values[kCentTPC]     = event->CentralityTPC();
  values[kCentZDC]     = event->CentralityZEMvsZDC();
  values[kCentQuality] = event->CentralityQuality();
  values[kNV0total]        = event->NV0CandidatesTotal();
  values[kNV0selected]     = event->NV0Candidates();
  values[kNdielectrons]    = event->NDielectrons();
  values[kNtracksTotal]    = event->NTracksTotal();
  values[kNtracksSelected] = event->NTracks();
  values[kSPDntracklets]   = event->SPDntracklets();
  values[kVZEROAemptyChannels] = 0;
  values[kVZEROCemptyChannels] = 0;
  for(Int_t ich=0;ich<64;++ich) gUsedVars[kVZEROChannelMult+ich] = kTRUE; 
  Float_t theta=0.0;
  for(Int_t ich=0;ich<64;++ich) {
    if(gUsedVars[kVZEROChannelMult+ich]) {
      values[kVZEROChannelMult+ich] = event->MultChannelVZERO(ich);
      if(!gVzeroAvMult[0] && values[kVZEROChannelMult+ich]<gkVZEROminMult) {
        gUsedVars[kVZEROChannelMult+ich] = kFALSE;   // will not be filled in histograms by the histogram manager
        if(ich<32) values[kVZEROCemptyChannels] += 1;
        else values[kVZEROAemptyChannels] += 1;
      }
    }
    if(gUsedVars[kVZEROChannelEta+ich]) {
      if(ich<32) theta = TMath::ATan(gkVZEROChannelRadii[ich]/(gkVZEROCz-values[kVtxZ]));
      else theta = TMath::Pi()-TMath::ATan(gkVZEROChannelRadii[ich]/(gkVZEROAz-values[kVtxZ]));
      values[kVZEROChannelEta+ich] = -1.0*TMath::Log(TMath::Tan(theta/2.0));
    }
  }
  
  if(eventF) {
    for(Int_t ih=0; ih<6; ++ih) {
      // VZERO event plane variables
      values[kVZEROQvecX+2*6+ih] = 0.0;
      values[kVZEROQvecY+2*6+ih] = 0.0;
      values[kVZERORP   +2*6+ih] = 0.0;
      for(Int_t iVZEROside=0; iVZEROside<2; ++iVZEROside) {
        values[kVZEROQvecX+iVZEROside*6+ih] = eventF->Qx(AliCorrelationReducedEventFriend::kVZEROA+iVZEROside, ih+1);
        values[kVZEROQvecY+iVZEROside*6+ih] = eventF->Qy(AliCorrelationReducedEventFriend::kVZEROA+iVZEROside, ih+1);
        values[kVZERORP   +iVZEROside*6+ih] = eventF->EventPlane(AliCorrelationReducedEventFriend::kVZEROA+iVZEROside, ih+1);
	values[kVZEROQvecX+2*6         +ih] += values[kVZEROQvecX+iVZEROside*6+ih];
	values[kVZEROQvecY+2*6         +ih] += values[kVZEROQvecY+iVZEROside*6+ih];
	// cos(n(EPtpc-EPvzero A/C))	
        values[kTPCRPres+iVZEROside*6+ih] = DeltaPhi(eventF->EventPlane(AliCorrelationReducedEventFriend::kTPC, ih+1), eventF->EventPlane(AliCorrelationReducedEventFriend::kVZEROA+iVZEROside, ih+1));
        values[kTPCRPres+iVZEROside*6+ih] = TMath::Cos(values[kTPCRPres+iVZEROside*6+ih]*(ih+1));
      }
      values[kVZERORP   +2*6+ih] = TMath::ATan2(values[kVZEROQvecY+2*6+ih],values[kVZEROQvecX+2*6+ih]);
      // cos (n*(psi_A-psi_C))
      values[kVZERORPres + ih] = DeltaPhi(eventF->EventPlane(AliCorrelationReducedEventFriend::kVZEROA, ih+1), 
					  eventF->EventPlane(AliCorrelationReducedEventFriend::kVZEROC, ih+1));
      values[kVZERORPres + ih] = TMath::Cos(values[kVZERORPres + ih]*(ih+1));
      // Qx,Qy correlations for VZERO
      values[kVZEROXaXc+ih] = eventF->Qx(AliCorrelationReducedEventFriend::kVZEROA, ih+1)*eventF->Qx(AliCorrelationReducedEventFriend::kVZEROC, ih+1);
      values[kVZEROXaYa+ih] = eventF->Qx(AliCorrelationReducedEventFriend::kVZEROA, ih+1)*eventF->Qy(AliCorrelationReducedEventFriend::kVZEROA, ih+1);
      values[kVZEROXaYc+ih] = eventF->Qx(AliCorrelationReducedEventFriend::kVZEROA, ih+1)*eventF->Qy(AliCorrelationReducedEventFriend::kVZEROC, ih+1);
      values[kVZEROYaXc+ih] = eventF->Qy(AliCorrelationReducedEventFriend::kVZEROA, ih+1)*eventF->Qx(AliCorrelationReducedEventFriend::kVZEROC, ih+1);
      values[kVZEROYaYc+ih] = eventF->Qy(AliCorrelationReducedEventFriend::kVZEROA, ih+1)*eventF->Qy(AliCorrelationReducedEventFriend::kVZEROC, ih+1);
      values[kVZEROXcYc+ih] = eventF->Qx(AliCorrelationReducedEventFriend::kVZEROC, ih+1)*eventF->Qy(AliCorrelationReducedEventFriend::kVZEROC, ih+1);
      // Psi_A - Psi_C
      values[kVZEROdeltaRPac+ih] = DeltaPhi(eventF->EventPlane(AliCorrelationReducedEventFriend::kVZEROA, ih+1), 
					    eventF->EventPlane(AliCorrelationReducedEventFriend::kVZEROC, ih+1));
      
      // TPC event plane
      values[kTPCQvecX+ih] = eventF->Qx(AliCorrelationReducedEventFriend::kTPC, ih+1);
      values[kTPCQvecY+ih] = eventF->Qy(AliCorrelationReducedEventFriend::kTPC, ih+1);
      values[kTPCRP   +ih] = eventF->EventPlane(AliCorrelationReducedEventFriend::kTPC, ih+1);
      // TPC VZERO Q-vector correlations
      values[kRPXtpcXvzeroa+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecX+ih];
      values[kRPXtpcXvzeroc+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecX+6+ih];
      values[kRPYtpcYvzeroa+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecY+ih];
      values[kRPYtpcYvzeroc+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecY+6+ih];
      values[kRPXtpcYvzeroa+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecY+ih];
      values[kRPXtpcYvzeroc+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecY+6+ih];
      values[kRPYtpcXvzeroa+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecX+ih];
      values[kRPYtpcXvzeroc+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecX+6+ih];
      // Psi_TPC - Psi_VZERO A/C      
      values[kRPdeltaVZEROAtpc+ih] = DeltaPhi(values[kVZERORP+0*6+ih], values[kTPCRP+ih]);
      values[kRPdeltaVZEROCtpc+ih] = DeltaPhi(values[kVZERORP+1*6+ih], values[kTPCRP+ih]);
    }  // end loop over harmonics
    
    // VZERO v2 using TPC event plane
    Double_t vzeroChannelPhi[8] = {0.3927, 1.1781, 1.9635, 2.7489, -2.7489, -1.9635, -1.1781, -0.3927};
    for(Int_t ich=0; ich<64; ++ich) {
      if(gUsedVars[kVZEROflowV2TPC]) {
        values[kVZEROflowV2TPC+ich] = values[kVZEROChannelMult+ich]*TMath::Cos(2.0*DeltaPhi(vzeroChannelPhi[ich%8],values[kTPCRP+1]));
      } 
    }
  }  // end if (eventF)  
}


//_________________________________________________________________
void DstCommonMacros::FillTrackingFlag(AliCorrelationReducedTrack* track, UShort_t flag, Float_t* values) {
  //
  // fill the tracking flag
  //
  values[kTrackingFlag] = -1;
  if(track->CheckTrackStatus(flag)) values[kTrackingFlag] = flag;
}


//_________________________________________________________________
void DstCommonMacros::FillEventOfflineTriggers(UShort_t triggerBit, Float_t* values) {
  //
  // fill the trigger bit input
  //
  if(triggerBit>=64) return;
  if(!gCurrentEvent) return;
  ULong64_t trigger = BIT(0);
  values[kOfflineTrigger] = triggerBit;
  values[kOfflineTriggerFired] = (gCurrentEvent->TriggerMask()&(trigger<<triggerBit) ? 1.0 : 0.0);
  values[kOfflineTriggerFired2] = (values[kOfflineTriggerFired]>0.01 ? triggerBit : -1.0); 
}


//_________________________________________________________________
void DstCommonMacros::FillTrackInfo(AliCorrelationReducedTrack* p, Float_t* values) {
  //
  // fill track information
  //
                        values[kPt]     = p->Pt();
                        values[kPtTPC]  = p->PtTPC();
  if(gUsedVars[kP])     values[kP]      = p->P();
  if(gUsedVars[kTheta]) values[kTheta]  = p->Theta();
                        values[kPhi]    = p->Phi();
                        values[kPhiTPC] = p->PhiTPC();
                        values[kEta]    = p->Eta();
                        values[kEtaTPC] = p->EtaTPC();
                        values[kPin]    = p->Pin();
                        values[kDcaXY]  = p->DCAxy();
                        values[kDcaZ]   = p->DCAz();
  
  if(gUsedVars[kITSncls]) values[kITSncls] = p->ITSncls();
                          values[kITSsignal] = p->ITSsignal();
  
  values[kTPCncls] = p->TPCncls();
  if(gUsedVars[kTPCnclsRatio]) 
    values[kTPCnclsRatio] = (p->TPCFindableNcls()>0 ? Float_t(p->TPCncls())/Float_t(p->TPCFindableNcls()) : 0.0);
  values[kTPCnclsIter1]   = p->TPCnclsIter1();
  values[kTPCnclsF]       = p->TPCFindableNcls();
  values[kTPCsignal]      = p->TPCsignal();
  
  values[kTOFbeta] = p->TOFbeta();
  for(Int_t specie=kElectron; specie<=kProton; ++specie) {
    values[kTPCnSig+specie] = p->TPCnSig(specie);
    values[kTOFnSig+specie] = p->TOFnSig(specie);
  }
  values[kTRDpidProbabilities]   = p->TRDpid(0);
  values[kTRDpidProbabilities+1] = p->TRDpid(1);
  
  values[kTRDntracklets]    = p->TRDntracklets(0);
  values[kTRDntrackletsPID] = p->TRDntracklets(1);
  
  if(gUsedVars[kEMCALmatchedEnergy] || gUsedVars[kEMCALmatchedEOverP]) {
    AliCorrelationReducedCaloCluster* cluster = gCurrentEvent->GetCaloCluster(p->CaloClusterId());
    values[kEMCALmatchedEnergy] = (cluster ? cluster->Energy() : -999.0);
    Float_t mom = p->P();
    values[kEMCALmatchedEnergy] = (TMath::Abs(mom)>1.e-8 && cluster ? values[kEMCALmatchedEOverP]/mom : -999.0);
  }
  
  // Fill track flow variables
  for(Int_t iVZEROside=0; iVZEROside<2; ++iVZEROside) {
    for(Int_t ih=0; ih<6; ++ih) {
      if(gUsedVars[kTrackVZEROFlowVn+iVZEROside*6+ih] || gUsedVars[kTrackDeltaPhiVZEROFlowVn+iVZEROside*6+ih]) {
        values[kTrackVZEROFlowVn+iVZEROside*6+ih]         = TMath::Cos(DeltaPhi(values[kPhi],values[kVZERORP+iVZEROside*6+ih])*(ih+1));
        values[kTrackDeltaPhiVZEROFlowVn+iVZEROside*6+ih] = DeltaPhi(values[kPhi],values[kVZERORP+iVZEROside*6+ih]);
      }
    }
  }
}


//_________________________________________________________________
void DstCommonMacros::FillCaloClusterInfo(AliCorrelationReducedCaloCluster* cl, Float_t* values) {
  //
  // Fill calorimeter cluster information
  //
  values[kEMCALclusterEnergy] = cl->Energy();
  values[kEMCALclusterDx] = cl->Dx();
  values[kEMCALclusterDz] = cl->Dz();
  values[kEMCALdetector] = (cl->IsEMCAL() ? AliCorrelationReducedCaloCluster::kEMCAL : AliCorrelationReducedCaloCluster::kPHOS);
}


//_________________________________________________________________
void DstCommonMacros::FillPairInfo(AliCorrelationReducedPair* p, Float_t* values) {
  //
  // fill pair information
  //
  
  values[kCandidateId] = p->CandidateId();
  values[kPairType]    = p->PairType();
  if(gUsedVars[kMass]) {
    values[kMass] = p->Mass();
    if(p->CandidateId()==AliCorrelationReducedPair::kLambda0ToPPi) values[kMass] = p->Mass(1);
    if(p->CandidateId()==AliCorrelationReducedPair::kALambda0ToPPi) values[kMass] = p->Mass(2);
  }
                            values[kPairPt]           = p->Pt();
  if(gUsedVars[kPairP])     values[kPairP]            = p->P();
  if(gUsedVars[kPairPx])    values[kPairPx]           = p->Px();
  if(gUsedVars[kPairPy])    values[kPairPy]           = p->Py();
  if(gUsedVars[kPairPz])    values[kPairPz]           = p->Pz();
                            values[kPairEta]          = p->Eta();
  if(gUsedVars[kPairRap])   values[kPairRap]          = p->Rapidity();
                            values[kPairPhi]          = p->Phi();
                            values[kPairLxy]          = p->Lxy();
                            values[kPairOpeningAngle] = p->OpeningAngle();
  if(gUsedVars[kPairTheta]) values[kPairTheta]        = p->Theta();
  
  // Flow variables
  // cos(n*phi), sin(n*phi)
  for(Int_t ih=0; ih<6; ++ih) {
    if(gUsedVars[kPairCosNPhi+ih]) values[kPairCosNPhi+ih] = TMath::Cos(values[kPairPhi]*(1.0+ih));
    if(gUsedVars[kPairSinNPhi+ih]) values[kPairSinNPhi+ih] = TMath::Sin(values[kPairPhi]*(1.0+ih));
  }
  // VZERO  
  for(Int_t iVZEROside=0; iVZEROside<3; ++iVZEROside) {
    // v3 using VZERO Psi_2
    if(gUsedVars[kPairVZEROFlowV3Psi2+iVZEROside]) 
      values[kPairVZEROFlowV3Psi2+iVZEROside] = TMath::Cos(3.0*DeltaPhi(values[kPairPhi],values[kVZERORP+iVZEROside*6+1]));
    for(Int_t ih=0; ih<6; ++ih) {
      // vn using VZERO Psi_n
      if(gUsedVars[kPairVZEROFlowVn+iVZEROside*6+ih]) 
	values[kPairVZEROFlowVn+iVZEROside*6+ih] = TMath::Cos(DeltaPhi(values[kPairPhi], values[kVZERORP+iVZEROside*6+ih])*(ih+1));
      // phi - Psi_n
      if(gUsedVars[kPairDeltaPhiVZEROFlowVn+iVZEROside*6+ih])
        values[kPairDeltaPhiVZEROFlowVn+iVZEROside*6+ih] = DeltaPhi(values[kPairPhi], values[kVZERORP+iVZEROside*6+ih]);
    }
  }
  // TPC
  // v3 using Psi_2
  if(gUsedVars[kPairTPCFlowV3Psi2]) 
    values[kPairTPCFlowV3Psi2] = TMath::Cos(3.0*DeltaPhi(values[kPairPhi],values[kTPCRP+1]));
  for(Int_t ih=0; ih<6; ++ih) {
    // vn using Psi_n
    if(gUsedVars[kPairTPCFlowVn+ih]) 
      values[kPairTPCFlowVn+ih] = TMath::Cos(DeltaPhi(values[kPairPhi],values[kTPCRP+ih])*(ih+1));
    if(gUsedVars[kPairDeltaPhiTPCFlowVn+ih])
      values[kPairDeltaPhiTPCFlowVn+ih] = DeltaPhi(values[kPairPhi], values[kTPCRP+ih]);
  }
}


//_________________________________________________________________
void DstCommonMacros::FillPairInfo(AliCorrelationReducedTrack* t1, AliCorrelationReducedTrack* t2, Int_t type, Float_t* values) {
  //
  // fill pair information from 2 tracks
  //
  // type - Parameter encoding the resonance type 
  //        This is needed for making a mass assumption on the legs
  //
  if(gUsedVars[kPairType])  {
    if(t1->Charge()*t2->Charge()<0) values[kPairType]  = 1;
    else if(t1->Charge()>0)         values[kPairType]  = 0;
    else                            values[kPairType]  = 2;
  }
  Float_t kMass1 = 0.0;
  Float_t kMass2 = 0.0;
  switch (type) {
    case AliCorrelationReducedPair::kK0sToPiPi :
      kMass1 = 0.13957; kMass2 = 0.13957;
      break;
    case AliCorrelationReducedPair::kPhiToKK :
      kMass1 = 0.493677; kMass2 = 0.493677;
      break;
    case AliCorrelationReducedPair::kLambda0ToPPi :
      kMass1 = 0.938272; kMass2 = 0.13957;
      break;
    case AliCorrelationReducedPair::kALambda0ToPPi :
      kMass1 = 0.13957; kMass2 = 0.938272;
      break;
    case AliCorrelationReducedPair::kJpsiToEE :
      kMass1 = 0.000511; kMass2 = 0.000511;
      break;
    default :
      break;
  }
  
  if(gUsedVars[kMass]) {     
    values[kMass]      = kMass1*kMass1+kMass2*kMass2 + 
                         2.0*(TMath::Sqrt(kMass1*kMass1+t1->P()*t1->P())*TMath::Sqrt(kMass2*kMass2+t2->P()*t2->P()) - 
                         t1->Px()*t2->Px() - t1->Py()*t2->Py() - t1->Pz()*t2->Pz());
    if(values[kMass]<0.0) {
      cout << "FillPairInfo(track, track, type, values): Warning: Very small squared mass found. "
           << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
      cout << "   mass2: " << values[kMass] << endl;
      cout << "p1(p,x,y,z): " << t1->P() << ", " << t1->Px() << ", " << t1->Py() << ", " << t1->Pz() << endl;
      cout << "p2(p,x,y,z): " << t2->P() << ", " << t2->Px() << ", " << t2->Py() << ", " << t2->Pz() << endl;
      values[kMass] = 0.0;
    }
    else
      values[kMass]      = TMath::Sqrt(values[kMass]);
  }
  
  if(gUsedVars[kPairPt])    values[kPairPt]    = TMath::Sqrt((t1->Px()+t2->Px())*(t1->Px()+t2->Px()) +
                                                             (t1->Py()+t2->Py())*(t1->Py()+t2->Py()));
  if(gUsedVars[kPairP])     values[kPairP]     = TMath::Sqrt((t1->Px()+t2->Px())*(t1->Px()+t2->Px()) +
                                                             (t1->Py()+t2->Py())*(t1->Py()+t2->Py()) +
                                                             (t1->Pz()+t2->Pz())*(t1->Pz()+t2->Pz()));
  if(gUsedVars[kPairEta])   {
    Float_t p = TMath::Sqrt((t1->Px()+t2->Px())*(t1->Px()+t2->Px()) +
                            (t1->Py()+t2->Py())*(t1->Py()+t2->Py()) +
                            (t1->Pz()+t2->Pz())*(t1->Pz()+t2->Pz()));
    values[kPairEta] = p-t1->Pz()-t2->Pz();
    values[kPairEta] = (TMath::Abs(values[kPairEta])>1.0e-8 ? (p+t1->Pz()+t2->Pz())/values[kPairEta] : 0.0);
    values[kPairEta]   = (values[kPairEta]>1.0e-8 ? 0.5*TMath::Log(values[kPairEta]) : -999.);
  }
  if(gUsedVars[kPairRap])   {
    Float_t mass = kMass1*kMass1+kMass2*kMass2 +
                   2.0*(TMath::Sqrt(kMass1*kMass1+t1->P()*t1->P())*TMath::Sqrt(kMass2*kMass2+t2->P()*t2->P()) - 
                   t1->Px()*t2->Px() - t1->Py()*t2->Py() - t1->Pz()*t2->Pz());
    if(mass<0.0) {
      cout << "Negative squared mass (Float_t resolution). Setting mass to zero" << endl;
      mass = 0.0;
    }
    else mass = TMath::Sqrt(mass);
    Float_t e = TMath::Sqrt(mass*mass+
                            (t1->Px()+t2->Px())*(t1->Px()+t2->Px()) +
                            (t1->Py()+t2->Py())*(t1->Py()+t2->Py()) +
                            (t1->Pz()+t2->Pz())*(t1->Pz()+t2->Pz()));
    values[kPairRap]   = 0.5*TMath::Log((e+t1->Pz()+t2->Pz())/(e-t1->Pz()-t2->Pz()));
  }
  if(gUsedVars[kPairPhi])   {
    values[kPairPhi]   = TMath::ATan2(t1->Py()+t2->Py(),t1->Px()+t2->Px());
    if(values[kPairPhi]<0.0) values[kPairPhi] = 2.0*TMath::Pi() + values[kPairPhi];
  }
  if(gUsedVars[kPairTheta]) values[kPairTheta] = TMath::ACos((t1->Pz()+t2->Pz())/
                                                             TMath::Sqrt((t1->Px()+t2->Px())*(t1->Px()+t2->Px()) +
                                                                         (t1->Py()+t2->Py())*(t1->Py()+t2->Py()) +
                                                                         (t1->Pz()+t2->Pz())*(t1->Pz()+t2->Pz())));
  // Flow variables
  // cos(n*phi), sin(n*phi)
  for(Int_t ih=0; ih<6; ++ih) {
    if(gUsedVars[kPairCosNPhi+ih]) values[kPairCosNPhi+ih] = TMath::Cos(values[kPairPhi]*(1.0+ih));
    if(gUsedVars[kPairSinNPhi+ih]) values[kPairSinNPhi+ih] = TMath::Sin(values[kPairPhi]*(1.0+ih));
  }
  // VZERO  
  for(Int_t iVZEROside=0; iVZEROside<3; ++iVZEROside) {
    // v3 using VZERO Psi_2
    if(gUsedVars[kPairVZEROFlowV3Psi2+iVZEROside]) 
      values[kPairVZEROFlowV3Psi2+iVZEROside] = TMath::Cos(3.0*DeltaPhi(values[kPairPhi],values[kVZERORP+iVZEROside*6+1]));
    for(Int_t ih=0; ih<6; ++ih) {
      // vn using VZERO Psi_n
      if(gUsedVars[kPairVZEROFlowVn+iVZEROside*6+ih]) 
	values[kPairVZEROFlowVn+iVZEROside*6+ih] = TMath::Cos(DeltaPhi(values[kPairPhi], values[kVZERORP+iVZEROside*6+ih])*(ih+1));
      // phi - Psi_n
      if(gUsedVars[kPairDeltaPhiVZEROFlowVn+iVZEROside*6+ih])
        values[kPairDeltaPhiVZEROFlowVn+iVZEROside*6+ih] = DeltaPhi(values[kPairPhi], values[kVZERORP+iVZEROside*6+ih]);
    }
  }
  // TPC
  // v3 using Psi_2
  if(gUsedVars[kPairTPCFlowV3Psi2]) 
    values[kPairTPCFlowV3Psi2] = TMath::Cos(3.0*DeltaPhi(values[kPairPhi],values[kTPCRP+1]));
  for(Int_t ih=0; ih<6; ++ih) {
    // vn using Psi_n
    if(gUsedVars[kPairTPCFlowVn+ih]) 
      values[kPairTPCFlowVn+ih] = TMath::Cos(DeltaPhi(values[kPairPhi],values[kTPCRP+ih])*(ih+1));
    if(gUsedVars[kPairDeltaPhiTPCFlowVn+ih])
      values[kPairDeltaPhiTPCFlowVn+ih] = DeltaPhi(values[kPairPhi], values[kTPCRP+ih]);
  }
}


//__________________________________________________________________
void DstCommonMacros::FillCorrelationInfo(AliCorrelationReducedPair* p, AliCorrelationReducedTrack* t, Float_t* values) {
  //
  // fill pair-track correlation information
  //
  if(gUsedVars[kDeltaPhi]) 
    values[kDeltaPhi] = DeltaPhi(p->Phi(), t->Phi());
  
  if(gUsedVars[kDeltaTheta]) 
    values[kDeltaTheta] = p->Theta() - t->Theta();
  
  if(gUsedVars[kDeltaEta]) values[kDeltaEta] = p->Eta() - t->Eta();
}


//__________________________________________________________________
void DstCommonMacros::DoEventMixing(TList* list1, TList* list2, Float_t* values, Int_t type, const Char_t* histClass) {
  //
  // Do the event mixing
  //
  if(type==0) return;
    
  cout << "mixing ..." << endl; 
  switch(type) {
    case -1:
      EventMixingPairTracks(list1, list2, values);    // list1 contains pairs; list2 contains tracks
      break;
    case AliCorrelationReducedPair::kJpsiToEE:
      // type is needed to encode the type of resonance which is needed for the mass assumption of legs when calculating the invariant mass
      EventMixingResonanceLegs(list1, list2, values, type, histClass);  // list1 contains positive legs; list2 contains negative legs
      break;
    default:
      break;
  }
  
  for(Int_t i=list1->GetEntries()-1; i>=0; --i) ((TList*)list1->At(i))->RemoveAll();
  for(Int_t j=list2->GetEntries()-1; j>=0; --j) ((TList*)list2->At(j))->RemoveAll();
  list1->Clear();
  list2->Clear();
}


//__________________________________________________________________
void DstCommonMacros::EventMixingPairTracks(TList* pairLists, TList* trackLists, Float_t* values) {
  //
  //  Make event mixing for pair - track correlations
  //
  cout << "Mixing category " << gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data() << endl;
  Int_t entries = pairLists->GetEntries();   // we should have the same number of entries in both pair and track lists
  if(entries<2) return;
  
  TList* pairs = 0x0;
  TList* tracks = 0x0;
  AliCorrelationReducedPair* pair = 0x0;
  AliCorrelationReducedTrack* track = 0x0;
  
  TIter nextPairList(pairLists);
  for(Int_t ie1=0; ie1<entries; ++ie1) {
    pairs = (TList*)nextPairList();
    
    TIter nextTrackList(trackLists); 
    for(Int_t ie2=0; ie2<entries; ++ie2) {
      tracks = (TList*)nextTrackList();
      if(ie1==ie2) continue;
      
      TIter nextPair(pairs);
      for(Int_t ip=0; ip<pairs->GetEntries(); ++ip) {
        pair = (AliCorrelationReducedPair*)nextPair();
        FillPairInfo(pair, values);

        TIter nextTrack(tracks);
	Int_t leadingIdx = -1;
	Float_t leadingPt = 0.0;
        for(Int_t it=0; it<tracks->GetEntries(); ++it) {
          track = (AliCorrelationReducedTrack*)nextTrack();
          FillTrackInfo(track, values);

          FillCorrelationInfo(pair, track, values);
          if(pair->PairType()==1) 
            FillHistClass(Form("Correlation_ME_US%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()),values);
          if(pair->PairType()==0) {
            FillHistClass(Form("Correlation_ME_LS%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()), values); 
            FillHistClass(Form("Correlation_ME_LSpp%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()), values);
          }
          if(pair->PairType()==2) {
            FillHistClass(Form("Correlation_ME_LS%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()), values); 
            FillHistClass(Form("Correlation_ME_LSmm%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()), values);
          }
          if(track->Pt()>leadingPt) {
	    leadingPt = track->Pt();
	    leadingIdx = it;
	  }
        }   // end loop over tracks
        if(leadingIdx!=-1) {
	  track = (AliCorrelationReducedTrack*)tracks->At(leadingIdx);
	  FillTrackInfo(track, values);
	  FillCorrelationInfo(pair, track, values);
          if(pair->PairType()==1) 
	    FillHistClass(Form("Correlation_LeadingPt_ME_US%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()),values);
          if(pair->PairType()==0) {
            FillHistClass(Form("Correlation_LeadingPt_ME_LS%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()),values);
            FillHistClass(Form("Correlation_LeadingPt_ME_LSpp%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()),values);
          }
          if(pair->PairType()==2) {
            FillHistClass(Form("Correlation_LeadingPt_ME_LS%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()),values);
            FillHistClass(Form("Correlation_LeadingPt_ME_LSmm%s", gEMCategoryNames[TMath::Nint(values[kEventMixingId])].Data()),values);
	  }
	}
      }   // end loop over pairs      
    }  // end loop over second event
  }  // end loop over first event
}


//__________________________________________________________________
void DstCommonMacros::EventMixingResonanceLegs(TList* posLists, TList* negLists, Float_t* values, Int_t type, const Char_t* histClass) {
  //
  // Do event mixing with the legs of a resonance to obtain the background
  //
  cout << "Mixing event class " << histClass << ", (cent/vtx): " << values[kCentVZERO] << "/" << values[kVtxZ] << endl;
  Int_t entries = posLists->GetEntries();
  cout << "mixing positives: " << entries << endl;
  cout << "mixing negatives: " << negLists->GetEntries() << endl;
  if(entries<2) return;
  
  TList* positives1 = 0x0;   // list of tracks in the first event
  //TList* negatives1 = 0x0;
  //TList* positives2 = 0x0;   // list of tracks in the second event
  TList* negatives2 = 0x0;
  AliCorrelationReducedTrack* posTrack1 = 0x0;
  //AliCorrelationReducedTrack* negTrack1 = 0x0;
  //AliCorrelationReducedTrack* posTrack2 = 0x0;
  AliCorrelationReducedTrack* negTrack2 = 0x0;
  
  TIter nextPosList1(posLists);
  TIter nextNegList1(negLists);
  for(Int_t ie1=0; ie1<entries; ++ie1) {
    positives1 = (TList*)nextPosList1();
    //negatives1 = (TList*)nextNegList1();
    
    TIter nextPosList2(posLists);
    TIter nextNegList2(negLists);
    for(Int_t ie2=0; ie2<entries; ++ie2) {
      //positives2 = (TList*)nextPosList2();
      negatives2 = (TList*)nextNegList2();
      if(ie1==ie2) continue;    // no SE mixing
      
      TIter nextPosLeg1(positives1);
      while((posTrack1=(AliCorrelationReducedTrack*)nextPosLeg1())) {
        /*TIter nextPosLeg2(positives2);
        while((posTrack2=(AliCorrelationReducedTrack*)nextPosLeg2())) {
          FillPairInfo(posTrack1, posTrack2, type, values);             // ++ combinations
          FillHistClass(Form("%s_LSpp", histClass),values);
        }    // end while over pos legs in event 2
        */
        TIter nextNegLeg2(negatives2);
        while((negTrack2=(AliCorrelationReducedTrack*)nextNegLeg2())) {
          FillPairInfo(posTrack1, negTrack2, type, values);             // +- combinations
	  if(IsPairSelectedEM(values))
            FillHistClass(histClass, values);
        }    // end while over neg legs in event 2
      }    // end while over pos legs in event 1
      /*TIter nextNegLeg1(negatives1);
      while((negTrack1=(AliCorrelationReducedTrack*)nextNegLeg1())) {
        TIter nextNegLeg2(negatives2);
        while((negTrack2=(AliCorrelationReducedTrack*)nextNegLeg2())) {
          FillPairInfo(negTrack1, negTrack2, type, values);             // -- combinations
          FillHistClass(Form("%s_LSmm", histClass),values);
        }    // end while over neg legs in event 2
      }    // end while over neg legs in event 1
      */
    }    // end for over event 2
  }    // end for over event 1
}




//________________________________________________________________
Double_t DstCommonMacros::DeltaPhi(Double_t phi1, Double_t phi2) {
  //
  // compute the delta of two angles defined in the (-pi,+pi) interval
  //
  Double_t delta = phi1-phi2;
  if(delta>2.0*TMath::Pi()) delta -= 2.0*TMath::Pi();
  if(delta<0.0) delta += 2.0*TMath::Pi();
  /*Double_t delta = phi2;
  if(phi2<0.0) delta += 2.0*TMath::Pi();
  delta = phi1-delta;
  if(delta>TMath::Pi()) delta = delta - 2.0*TMath::Pi();
  if(delta<-1.*TMath::Pi()) delta = 2.0*TMath::Pi() + delta;
  */
  return delta;
}


//__________________________________________________________________
void DstCommonMacros::FillHistClass(const Char_t* className, Float_t* values) {
  //
  //  fill a class of histograms
  //
  THashList* hList = (THashList*)gHistLists->FindObject(className);
  if(!hList) {
    //cout << "Warning in FillHistClass(): Histogram class " << className << " not found" << endl;
    return;
  }
    
  TIter next(hList);
  TObject* h=0x0;
  while((h=next())) {
    UInt_t id = h->GetUniqueID();
    UInt_t varX = id%kNVars;
    UInt_t varY = (id/kNVars)%kNVars;
    UInt_t varZ = (id/kNVars/kNVars)%kNVars;
    if(id<kNVars) { 
      ((TH1F*)h)->Fill(values[varX]);
    } else if(id<kNVars*kNVars) {
      if(((TH1*)h)->GetDimension()==1) ((TProfile*)h)->Fill(values[varX],values[varY]);
      if(((TH1*)h)->GetDimension()==2) ((TH2F*)h)->Fill(values[varX],values[varY]);
    } else {
      if(((TH1*)h)->GetDimension()==2) ((TProfile2D*)h)->Fill(values[varX],values[varY],values[varZ]);
      if(((TH1*)h)->GetDimension()==3) ((TH3F*)h)->Fill(values[varX],values[varY],values[varZ]);
    };
  }
}


//__________________________________________________________________
void DstCommonMacros::AddHistClass(const Char_t* histClass) {
  //
  // Add a new histogram list
  //
  if(!gHistLists)
    gHistLists = new TObjArray();
  gHistLists->SetOwner();
  gHistLists->SetName("histos");
  
  if(gHistLists->FindObject(histClass)) {
    cout << "Warning in AddHistClass: Cannot add histogram class " << histClass
         << " because it already exists." << endl;
    return;
  }
  THashList* hList=new THashList;
  hList->SetOwner(kTRUE);
  hList->SetName(histClass);
  gHistLists->Add(hList);
}


//_________________________________________________________________
Int_t DstCommonMacros::ValidateHistogramName(THashList* hList, const Char_t* name) {
  //
  // check whether a histogram with this name already exist, and how many of them
  //
  for(Int_t i=0; i<hList->GetEntries(); ++i) {
    TString hname = hList->At(i)->GetName();
    TObjArray* arr=hname.Tokenize("#");
    TString nameRoot = arr->At(0)->GetName();
    if(nameRoot.CompareTo(name)==0) {
      cout << "Warning in AddHistogram(): Histogram " << name << " already exists in this list (" << nameRoot.Data() << ")" << endl;
      return -1;
    }
  }
  Int_t nnames = 0;
  for(Int_t il=0; il<gHistLists->GetEntries(); ++il) {
    THashList* tlist = (THashList*)gHistLists->At(il);
    for(Int_t ih=0; ih<tlist->GetEntries(); ++ih) {
      TString hname = tlist->At(ih)->GetName();
      TObjArray* arr=hname.Tokenize("#");
      TString nameRoot = arr->At(0)->GetName();
      if(nameRoot.CompareTo(name)==0) ++nnames;
    }
  }
  return nnames;
}


//_________________________________________________________________
void DstCommonMacros::AddHistogram(const Char_t* histClass,
		                   const Char_t* name, const Char_t* title, Bool_t isProfile,
                                   Int_t nXbins, Double_t xmin, Double_t xmax, Int_t varX,
		                   Int_t nYbins, Double_t ymin, Double_t ymax, Int_t varY,
		                   Int_t nZbins, Double_t zmin, Double_t zmax, Int_t varZ,
		                   const Char_t* xLabels, const Char_t* yLabels, const Char_t* zLabels) {
  //
  // add a histogram
  //
  THashList* hList = (THashList*)gHistLists->FindObject(histClass);
  if(hList->FindObject(name)) {
    cout << "Warning in AddHistogram(): Histogram " << name << " already exists" << endl;
    return;
  }
  /*Int_t nnames = ValidateHistogramName(hList, name);
  if(nnames<0) return;
  TString hname = Form("%s#v%d", name, nnames);
  cout << "Assigned name " << hname.Data() << endl;*/
  TString hname = name;
  
  Int_t dimension = 1;
  if(varY!=kNothing) dimension = 2;
  if(varZ!=kNothing) dimension = 3;
  
  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");
  
  TH1* h=0x0;
  switch(dimension) {
    case 1:
      h=new TH1F(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax);
      //cout << "HISTOGRAM " << hname.Data() << ", ID: " << UInt_t(varX) << endl;
      h->SetUniqueID(UInt_t(varX));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      gUsedVars[varX] = kTRUE;
      hList->Add(h);
      break;
    case 2:
      if(isProfile)
	h=new TProfile(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax);
      else
	h=new TH2F(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax);
      //cout << "HISTOGRAM " << hname.Data() << ", ID: " << UInt_t(varX+kNVars*varY) << endl;
      h->SetUniqueID(UInt_t(varX+kNVars*varY));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      gUsedVars[varX] = kTRUE;
      gUsedVars[varY] = kTRUE;
      hList->Add(h);
      break;
    case 3:
      if(isProfile)
	h=new TProfile2D(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax);
      else
	h=new TH3F(hname.Data(),arr->At(0)->GetName(),nXbins,xmin,xmax,nYbins,ymin,ymax,nZbins,zmin,zmax);
      //cout << "HISTOGRAM " << hname.Data() << ", ID: " << UInt_t(varX+kNVars*varY+kNVars*kNVars*varZ) << endl;
      h->SetUniqueID(UInt_t(varX+kNVars*varY+kNVars*kNVars*varZ));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      if(arr->At(3)) h->GetZaxis()->SetTitle(arr->At(3)->GetName());
      if(zLabels[0]!='\0') MakeAxisLabels(h->GetZaxis(), zLabels);
      gUsedVars[varX] = kTRUE;
      gUsedVars[varY] = kTRUE;
      gUsedVars[varZ] = kTRUE;
      hList->Add(h);
      break;
  }
}

//_________________________________________________________________
void DstCommonMacros::AddHistogram(const Char_t* histClass,
		                   const Char_t* name, const Char_t* title, Bool_t isProfile,
                                   Int_t nXbins, Double_t* xbins, Int_t varX,
		                   Int_t nYbins, Double_t* ybins, Int_t varY,
		                   Int_t nZbins, Double_t* zbins, Int_t varZ,
		                   const Char_t* xLabels, const Char_t* yLabels, const Char_t* zLabels) {
  //
  // add a histogram
  //
  THashList* hList = (THashList*)gHistLists->FindObject(histClass);
  if(hList->FindObject(name)) {
    cout << "Warning in AddHistogram(): Histogram " << name << " already exists" << endl;
    return;
  }
  //Int_t nnames = ValidateHistogramName(hList, name);
  //if(nnames<0) return;
  //TString hname = Form("%s#v%d", name, nnames);
  TString hname = name;
  
  Int_t dimension = 1;
  if(varY!=kNothing) dimension = 2;
  if(varZ!=kNothing) dimension = 3;
  
  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");
  
  TH1* h=0x0;
  switch(dimension) {
    case 1:
      h=new TH1F(hname.Data(),arr->At(0)->GetName(),nXbins,xbins);
      //cout << "HISTOGRAM " << hname.Data() << ", ID: " << UInt_t(varX) << endl;
      h->SetUniqueID(UInt_t(varX));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      gUsedVars[varX] = kTRUE;
      hList->Add(h);
      break;
    case 2:
      if(isProfile)
	h=new TProfile(hname.Data(),arr->At(0)->GetName(),nXbins,xbins);
      else
	h=new TH2F(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins);
      //cout << "HISTOGRAM " << hname.Data() << ", ID: " << UInt_t(varX+kNVars*varY) << endl;
      h->SetUniqueID(UInt_t(varX+kNVars*varY));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      gUsedVars[varX] = kTRUE;
      gUsedVars[varY] = kTRUE;
      hList->Add(h);
      break;
    case 3:
      if(isProfile)
	h=new TProfile2D(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins);
      else
	h=new TH3F(hname.Data(),arr->At(0)->GetName(),nXbins,xbins,nYbins,ybins,nZbins,zbins);
      //cout << "HISTOGRAM " << hname.Data() << ", ID: " << UInt_t(varX+kNVars*varY+kNVars*kNVars*varZ) << endl;
      h->SetUniqueID(UInt_t(varX+kNVars*varY+kNVars*kNVars*varZ));
      if(arr->At(1)) h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      if(xLabels[0]!='\0') MakeAxisLabels(h->GetXaxis(), xLabels);
      if(arr->At(2)) h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      if(yLabels[0]!='\0') MakeAxisLabels(h->GetYaxis(), yLabels);
      if(arr->At(3)) h->GetZaxis()->SetTitle(arr->At(3)->GetName());
      if(zLabels[0]!='\0') MakeAxisLabels(h->GetZaxis(), zLabels);
      gUsedVars[varX] = kTRUE;
      gUsedVars[varY] = kTRUE;
      gUsedVars[varZ] = kTRUE;
      hList->Add(h);
      break;
  }
}


//____________________________________________________________________________________
void DstCommonMacros::MakeAxisLabels(TAxis* ax, const Char_t* labels) {
  //
  // add bin labels to an axis
  //
  TString labelsStr(labels);
  TObjArray* arr=labelsStr.Tokenize(";");
  for(Int_t ib=1; ib<=ax->GetNbins(); ++ib) {
    if(ib>=arr->GetEntries()+1) break;
    ax->SetBinLabel(ib, arr->At(ib-1)->GetName());
  }
}


//____________________________________________________________________________________
void DstCommonMacros::InitFile(const Char_t* filename) {
  //
  // Open an existing ROOT file containing lists of histograms and initialize the global list pointer
  //
  //if(gHistFile) gHistFile->Close();
  TString histfilename="";
  if(gHistFile) histfilename = gHistFile->GetName();
  if(!histfilename.Contains(filename)) 
    gHistFile = new TFile(filename);    // open file only if not already open
  
  if(!gHistFile) {
    cout << "GlobalMacros.C::GetHistogram() : File " << filename << " not opened!!" << endl;
    return;
  }
  if(gHistFile->IsZombie()) {
    cout << "GlobalMacros.C::GetHistogram() : File " << filename << " not opened!!" << endl;
    return;
  }
  TList* list1 = gHistFile->GetListOfKeys();
  TKey* key1 = (TKey*)list1->At(0);
  gHistListsOld = (TDirectoryFile*)key1->ReadObj();
}


//____________________________________________________________________________________
void DstCommonMacros::CloseFile() {
  //
  // Close the opened file
  //
  gHistListsOld = 0x0;
  if(gHistFile && gHistFile->IsOpen()) gHistFile->Close();
}


//____________________________________________________________________________________
TObject* DstCommonMacros::GetHistogram(const Char_t* listname, const Char_t* hname) {
  //
  // Retrieve a histogram from the list hlist
  //
  if(!gHistListsOld) {
    cout << "GlobalMacros.C::GetHistogram() : The main list was not initialized." << endl;
    cout << "                   A ROOT file must pe initialized first!!" << endl;
    return 0x0;
  }
  TKey* listKey = gHistListsOld->FindKey(listname);
  TDirectoryFile* hlist = (TDirectoryFile*)listKey->ReadObj();
  TKey* key = hlist->FindKey(hname);
  return key->ReadObj();
}


//_________________________________________________________________
TChain* DstCommonMacros::GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                  TChain* friendChain/*=0x0*/, const Char_t* friendChainFile/*=0x0*/) {
  //
  // read an ascii file containing a list of root files with reduced trees
  // and build a TChain
  //
  cout << "Creating the data chain from " << filename << " ..." << flush; 
  TChain* chain = new TChain("DstTree");
  ifstream inBuf;
  inBuf.open(filename);
  Int_t index = 0;
  while(inBuf.good()) {
    Char_t str[512];
    inBuf.getline(str,512,'\n');
    
    if(index<offset) {++index; continue;}
    if(index>=offset+howMany) break;
    
    TString strstr = str;
    if(!strstr.Contains(".root")) continue;
    cout << endl << "Adding file " << str << endl;
    chain->Add(str);
    if(friendChain) {
      TObjArray* arr = strstr.Tokenize("/");
      strstr.ReplaceAll(arr->At(arr->GetEntries()-1)->GetName(),friendChainFile);
      friendChain->Add(strstr.Data());
    }    
    ++index;
  }
  inBuf.close();
  entries = chain->GetEntries();
  Long64_t entriesFriend = (friendChain ? friendChain->GetEntries() : 0);
  cout << "DstCommonMacros::GetChain() Chain entries = " << entries << endl;
  if(friendChain)
    cout << "DstCommonMacros::GetChain() Friend chain entries = " << entriesFriend << endl;
  if(friendChain && (entries!=entriesFriend)) {
    cout << "DstCommonMacros::GetChain() The friend chain does not have the same number of entries as the main chain!!!" << endl;
    cout << "                            Check it out and retry !!" << endl;
    return 0x0;
  }
  cout << " done" << endl;
  return chain;
}


//__________________________________________________________________
void DstCommonMacros::WriteOutput(TFile* save) {
  //
  // Write the histogram lists in the output file
  //
  cout << "Writing the output to " << save->GetName() << " ... " << flush;
  //TFile* save=new TFile(filename,"RECREATE");
  TDirectory* mainDir = save->mkdir(gHistLists->GetName());
  mainDir->cd();
  for(Int_t i=0; i<gHistLists->GetEntries(); ++i) {
    THashList* list = (THashList*)gHistLists->At(i);
    TDirectory* dir = mainDir->mkdir(list->GetName());
    dir->cd();
    list->Write();
    mainDir->cd();
  }
  save->Close();
  cout << "done" << endl;
}


//__________________________________________________________________
Bool_t DstCommonMacros::IsPairSelectedEM(Float_t* values) {
  //
  // pair selection on the pair resulting from the event mixing
  //
  if(values[kPairPt]<gkJpsiPtCut[0]) return kFALSE;
  if(values[kPairPt]>gkJpsiPtCut[1]) return kFALSE;
  return kTRUE;
}
