/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections ALICE glue                                                *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
/***********************************************************
    Variable definitions for event plane correction framework
    Based on work of Ionut-Cristian Arsene
 ***********************************************************/

#include <TMath.h>

#include "AliAnalysisTaskSE.h"
#include "AliQnCorrectionsVarManagerTask.h"

ClassImp(AliQnCorrectionsVarManagerTask)

const Char_t* AliQnCorrectionsVarManagerTask::fTrackingFlagNames[] = {
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


const Char_t* AliQnCorrectionsVarManagerTask::fOfflineTriggerNames[64] = {
    "MB",              "INT7",              "MUON", "HighMult",    "EMC1", "CINT5",       "CMUS5/MUSPB",      "MUSH7/MUSHPB",
    "MUL7/MuonLikePB", "MUU7/MuonUnlikePB", "EMC7", "MUS7",        "PHI1", "PHI7/PHOSPb", "EMCEJE",           "EMCEGA",
    "Central",         "SemiCentral",       "DG5",  "ZED",         "SPI7", "CINT8",       "MuonSingleLowPt8", "MuonSingleHighPt8",
    "MuonLikeLowPt8",  "MuonUnlikeLowPt8",  "N/A",  "UserDefined", "N/A",  "N/A",         "FastOnly",         "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",              "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",              "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",              "N/A",
    "N/A",             "N/A",               "N/A",  "N/A",         "N/A",  "N/A",         "N/A",              "N/A"
};

//__________________________________________________________________
void AliQnCorrectionsVarManagerTask::SetDefaultVarNames() {
  fVariableNames[kRandom1][0]              = "User";                            fVariableNames[kRandom1][1] = "";
  fVariableNames[kRunNo][0]                = "Run number";                      fVariableNames[kRunNo][1] = "";
  fVariableNames[kLHCFillNumber][0]        = "LHC fill number";                 fVariableNames[kLHCFillNumber][1] = ""; 
  fVariableNames[kBeamEnergy][0]           = "Beam energy";                     fVariableNames[kBeamEnergy][1] = "GeV";
  fVariableNames[kDetectorMask][0]         = "Detector mask";                   fVariableNames[kDetectorMask][1] = "";
  fVariableNames[kNumberOfDetectors][0]    = "Number of active detectors";      fVariableNames[kNumberOfDetectors][1] = "";  
  fVariableNames[kDipolePolarity][0]       = "Dipole magnet polarity";          fVariableNames[kDipolePolarity][1] = "";
  fVariableNames[kL3Polarity][0]           = "L3 magnet polarity";              fVariableNames[kL3Polarity][1] = "";
  fVariableNames[kBeamIntensity0][0]       = "Beam 0 intensity";                fVariableNames[kBeamIntensity0][1] = "";
  fVariableNames[kBeamIntensity1][0]       = "Beam 1 intensity";                fVariableNames[kBeamIntensity1][1] = "";
  fVariableNames[kBC][0]                   = "Bunch crossing";                  fVariableNames[kBC][1] = "";
  fVariableNames[kTimeStamp][0]            = "Time stamp";                      fVariableNames[kTimeStamp][1] = "";
  fVariableNames[kEventType][0]            = "Event type";                      fVariableNames[kEventType][1] = "";
  fVariableNames[kTriggerMask][0]          = "Trigger mask";                    fVariableNames[kTriggerMask][1] = "";
  fVariableNames[kOfflineTrigger][0]       = "Offline trigger";                 fVariableNames[kOfflineTrigger][1] = "";
  fVariableNames[kOfflineTriggerFired][0]  = "Offline trigger fired";           fVariableNames[kOfflineTriggerFired][1] = "";
  fVariableNames[kOfflineTriggerFired2][0] = "Offline trigger fired2";          fVariableNames[kOfflineTriggerFired2][1] = "";
  fVariableNames[kIsPhysicsSelection][0]   = "Physics selection ON";            fVariableNames[kIsPhysicsSelection][1] = "";
  fVariableNames[kIsSPDPileup][0]          = "SPD pileup ON";                   fVariableNames[kIsSPDPileup][1] = "";
  fVariableNames[kNSPDpileups][0]          = "Number of SPD pileup events";     fVariableNames[kNSPDpileups][1] = "";
  fVariableNames[kNTrackPileups][0]        = "Number of track pileup events";   fVariableNames[kNTrackPileups][1] = "";
  fVariableNames[kNPMDtracks][0]           = "Number of PMD tracks";            fVariableNames[kNPMDtracks][1] = "";
  fVariableNames[kNTRDtracks][0]           = "Number of TRD tracks";            fVariableNames[kNTRDtracks][1] = "";
  fVariableNames[kNTRDtracklets][0]        = "Number of TRD tracklets";         fVariableNames[kNTRDtracklets][1] = "";
  fVariableNames[kNVtxContributors][0]     = "Number of vtx. contributors";     fVariableNames[kNVtxContributors][1] = "";
  fVariableNames[kNVtxTPCContributors][0]  = "Number of TPC vtx. contributors"; fVariableNames[kNVtxTPCContributors][1] = "";
  fVariableNames[kVtxX][0]                 = "Vtx X";                           fVariableNames[kVtxX][1] = "cm";
  fVariableNames[kVtxY][0]                 = "Vtx Y";                           fVariableNames[kVtxY][1] = "cm";
  fVariableNames[kVtxZ][0]                 = "Vtx Z";                           fVariableNames[kVtxZ][1] = "cm";
  fVariableNames[kVtxXtpc][0]              = "Vtx X TPC";                       fVariableNames[kVtxXtpc][1] = "cm";
  fVariableNames[kVtxYtpc][0]              = "Vtx Y TPC";                       fVariableNames[kVtxXtpc][1] = "cm";
  fVariableNames[kVtxZtpc][0]              = "Vtx Z TPC";                       fVariableNames[kVtxXtpc][1] = "cm";
  fVariableNames[kDeltaVtxZ][0]            = "#Delta Z";                        fVariableNames[kDeltaVtxZ][1] = "cm";
  for(Int_t iflag=0; iflag<kNTrackingFlags; ++iflag) {
    fVariableNames[kNTracksPerTrackingFlag+iflag][0] = Form("Tracks with %s on",fTrackingFlagNames[iflag]); 
    fVariableNames[kNTracksPerTrackingFlag+iflag][1] = ""; 
  }
  fVariableNames[kNTracksTPCoutVsITSout][0]       = "TPCout/ITSout";                   fVariableNames[kNTracksTPCoutVsITSout][1] = "";
  fVariableNames[kNTracksTRDoutVsITSout][0]       = "TRDout/ITSout";                   fVariableNames[kNTracksTRDoutVsITSout][1] = "";
  fVariableNames[kNTracksTOFoutVsITSout][0]       = "TOFout/ITSout";                   fVariableNames[kNTracksTOFoutVsITSout][1] = "";
  fVariableNames[kNTracksTRDoutVsTPCout][0]       = "TRDout/TPCout";                   fVariableNames[kNTracksTRDoutVsTPCout][1] = "";
  fVariableNames[kNTracksTOFoutVsTPCout][0]       = "TOFout/TPCout";                   fVariableNames[kNTracksTOFoutVsTPCout][1] = "";
  fVariableNames[kNTracksTOFoutVsTRDout][0]       = "TOFout/TRDout";                   fVariableNames[kNTracksTOFoutVsTRDout][1] = "";
  fVariableNames[kNTracksITSoutVsSPDtracklets][0] = "ITSout/SPDtracklets";             fVariableNames[kNTracksITSoutVsSPDtracklets][1] = "";
  fVariableNames[kNTracksTPCoutVsSPDtracklets][0] = "TPCout/SPDtracklets";             fVariableNames[kNTracksTPCoutVsSPDtracklets][1] = "";
  fVariableNames[kNTracksTRDoutVsSPDtracklets][0] = "TRDout/SPDtracklets";             fVariableNames[kNTracksTRDoutVsSPDtracklets][1] = "";
  fVariableNames[kNTracksTOFoutVsSPDtracklets][0] = "TOFout/SPDtracklets";             fVariableNames[kNTracksTOFoutVsSPDtracklets][1] = "";
  fVariableNames[kCentVZERO][0]                   = "VZERO centrality";                fVariableNames[kCentVZERO][1] = "%";
  fVariableNames[kCentSPD][0]                     = "CL1 centrality";                  fVariableNames[kCentSPD][1] = "%";
  fVariableNames[kCentSPDcorr][0]                 = "SPD trklts centrality";           fVariableNames[kCentSPDcorr][1] = "%";
  fVariableNames[kCentTPC][0]                     = "TPC centrality";                  fVariableNames[kCentTPC][1] = "%";
  fVariableNames[kCentZDC][0]                     = "ZDC centrality";                  fVariableNames[kCentZDC][1] = "%";
  fVariableNames[kCentQuality][0]                 = "Centrality quality";              fVariableNames[kCentQuality][1] = "";
  fVariableNames[kNV0total][0]                    = "Total number of V0s";             fVariableNames[kNV0total][1] = "";  
  fVariableNames[kNV0selected][0]                 = "Number of selected V0s";          fVariableNames[kNV0selected][1] = "";  
  fVariableNames[kNdielectrons][0]                = "Number of dielectrons";           fVariableNames[kNdielectrons][1] = "";  
  fVariableNames[kNpairsSelected][0]              = "Number of pairs per event";       fVariableNames[kNpairsSelected][1] = "";    
  fVariableNames[kNtracksTotal][0]                = "No. of tracks in original event"; fVariableNames[kNtracksTotal][1] = "";
  fVariableNames[kNtracksSelected][0]             = "No. of selected tracks";          fVariableNames[kNtracksSelected][1] = "";
  fVariableNames[kNtracksPosAnalyzed][0]          = "No.selected positive tracks";     fVariableNames[kNtracksPosAnalyzed][1] = "";  
  fVariableNames[kNtracksNegAnalyzed][0]          = "No.selected negative tracks";     fVariableNames[kNtracksNegAnalyzed][1] = "";  
  fVariableNames[kNtracksAnalyzed][0]             = "No.selected tracks";              fVariableNames[kNtracksAnalyzed][1] = "";  
  fVariableNames[kNtracksSubEvLeft][0]            = "No.tracks sub-event left";        fVariableNames[kNtracksSubEvLeft][1] = "";  
  fVariableNames[kNtracksSubEvRight][0]           = "No.tracks sub-event right";       fVariableNames[kNtracksSubEvRight][1] = "";  
  fVariableNames[kNtracksEventPlane][0]           = "No.tracks";                       fVariableNames[kNtracksEventPlane][1] = "";  
  fVariableNames[kSPDnSingleClusters][0]          = "No.SPD single clusters";          fVariableNames[kSPDnSingleClusters][1] = "";  
  fVariableNames[kSPDntracklets][0]               = "No.SPD tracklets";                fVariableNames[kSPDntracklets][1] = "";  
  fVariableNames[kSPDntrackletsCorr][0]           = "No.corrected SPD tracklets";      fVariableNames[kSPDntrackletsCorr][1] = "";
  for(Int_t ieta=0;ieta<16;++ieta) {
    fVariableNames[kSPDntrackletsEta+ieta][0] = Form("No.SPD tracklets in %.1f<#eta<%.1f", -1.6+0.2*ieta, -1.6+0.2*(ieta+1));
    fVariableNames[kSPDntrackletsEta+ieta][1] = "";
  }  
  fVariableNames[kSPDtrackletEta][0]      = "SPD tracklet #eta";                fVariableNames[kSPDtrackletEta][1] = "";  
  fVariableNames[kSPDtrackletPhi][0]      = "SPD tracklet #phi";                fVariableNames[kSPDtrackletPhi][1] = "";  
  fVariableNames[kEventMixingId][0]       = "Event mixing id";        fVariableNames[kEventMixingId][1] = "";  
  fVariableNames[kVZEROATotalMult][0]     = "Multiplicity VZERO-A";   fVariableNames[kVZEROATotalMult][1] = "";
  fVariableNames[kVZEROCTotalMult][0]     = "Multiplicity VZERO-C";   fVariableNames[kVZEROCTotalMult][1] = "";
  fVariableNames[kVZEROTotalMult][0]      = "Multiplicity VZERO";     fVariableNames[kVZEROTotalMult][1] = "";
  fVariableNames[kVZEROMultPercentile][0] = "Multiplicity % VZERO";   fVariableNames[kVZEROMultPercentile][1] = "";
  fVariableNames[kTZEROATotalMult][0]     = "Multiplicity TZERO-A";   fVariableNames[kTZEROATotalMult][1] = "";
  fVariableNames[kTZEROCTotalMult][0]     = "Multiplicity TZERO-C";   fVariableNames[kTZEROCTotalMult][1] = "";
  fVariableNames[kTZEROTotalMult][0]      = "Multiplicity TZERO";     fVariableNames[kTZEROTotalMult][1] = "";
  fVariableNames[kFMD1TotalMult][0]      = "Multiplicity FMD1";     fVariableNames[kFMD1TotalMult][1] = "";
  fVariableNames[kFMD2ITotalMult][0]      = "Multiplicity FMD2I";     fVariableNames[kFMD2ITotalMult][1] = "";
  fVariableNames[kFMD2OTotalMult][0]      = "Multiplicity FMD2O";     fVariableNames[kFMD2OTotalMult][1] = "";
  fVariableNames[kFMD3ITotalMult][0]      = "Multiplicity FMD3I";     fVariableNames[kFMD3ITotalMult][1] = "";
  fVariableNames[kFMD3OTotalMult][0]      = "Multiplicity FMD3O";     fVariableNames[kFMD3OTotalMult][1] = "";
  fVariableNames[kVZEROAemptyChannels][0] = "VZERO-A empty channels"; fVariableNames[kVZEROAemptyChannels][1] = "";
  fVariableNames[kVZEROCemptyChannels][0] = "VZERO-C empty channels"; fVariableNames[kVZEROCemptyChannels][1] = "";
  for(Int_t ich=0;ich<64;++ich) {
    fVariableNames[kVZEROChannelMult+ich][0] = Form("Multiplicity VZERO ch.%d", ich);
    fVariableNames[kVZEROChannelMult+ich][1] = "";
    fVariableNames[kVZEROChannelEta+ich][0] = Form("#eta for VZERO ch.%d", ich);
    fVariableNames[kVZEROChannelEta+ich][1] = "";
  }  
  TString vzeroSideNames[3] = {"A","C","AC"};
  for(Int_t iHarmonic=0;iHarmonic<12;++iHarmonic) {
    fVariableNames[kCosNPhi+iHarmonic][0] = Form("cos(%d#varphi)",iHarmonic+1); fVariableNames[kCosNPhi+iHarmonic][1] = "";
    fVariableNames[kSinNPhi+iHarmonic][0] = Form("sin(%d#varphi)",iHarmonic+1); fVariableNames[kSinNPhi+iHarmonic][1] = "";
  }

  fVariableNames[kPt][0] = "p_{T}"; fVariableNames[kPt][1] = "GeV/c";
  fVariableNames[kP][0] = "p"; fVariableNames[kP][1] = "GeV/c";
  fVariableNames[kPx][0] = "p_{x}"; fVariableNames[kPx][1] = "GeV/c";
  fVariableNames[kPy][0] = "p_{y}"; fVariableNames[kPy][1] = "GeV/c";
  fVariableNames[kPz][0] = "p_{z}"; fVariableNames[kPz][1] = "GeV/c";
  fVariableNames[kCharge][0] = "Charge"; fVariableNames[kPz][1] = "e";
  fVariableNames[kTheta][0] = "#theta"; fVariableNames[kTheta][1] = "rad.";
  fVariableNames[kPhi][0] = "#varphi"; fVariableNames[kPhi][1] = "rad.";
  fVariableNames[kEta][0] = "#eta"; fVariableNames[kEta][1] = "";
  fVariableNames[kRap][0] = "y"; fVariableNames[kRap][1] = "";
  fVariableNames[kPtTPC][0] = "p_{T}^{TPC}"; fVariableNames[kPtTPC][1] = "GeV/c";
  fVariableNames[kPhiTPC][0] = "#varphi^{TPC}"; fVariableNames[kPhiTPC][1] = "rad.";
  fVariableNames[kEtaTPC][0] = "#eta^{TPC}"; fVariableNames[kEtaTPC][1] = "";
  fVariableNames[kPin][0] = "p_{IN}"; fVariableNames[kPin][1] = "GeV/c";
  fVariableNames[kDcaXY][0] = "DCA_{xy}"; fVariableNames[kDcaXY][1] = "cm.";  
  fVariableNames[kDcaZ][0] = "DCA_{z}"; fVariableNames[kDcaZ][1] = "cm.";  
  fVariableNames[kITSncls][0] = "No.ITS clusters"; fVariableNames[kITSncls][1] = "";    
  fVariableNames[kITSlayerHit][0] = "ITS layer"; fVariableNames[kITSlayerHit][1] = "";
  fVariableNames[kITSsignal][0] = "ITS dE/dx"; fVariableNames[kITSsignal][1] = "";    
  fVariableNames[kITSnSig][0] = "ITS n_{#sigma}^{e}"; fVariableNames[kITSnSig][1] = "#sigma";
  fVariableNames[kITSnSig+1][0] = "ITS n_{#sigma}^{#pi}"; fVariableNames[kITSnSig+1][1] = "#sigma";
  fVariableNames[kITSnSig+2][0] = "ITS n_{#sigma}^{K}"; fVariableNames[kITSnSig+2][1] = "#sigma";
  fVariableNames[kITSnSig+3][0] = "ITS n_{#sigma}^{p}"; fVariableNames[kITSnSig+3][1] = "#sigma";
  fVariableNames[kTPCncls][0] = "No.TPC clusters"; fVariableNames[kTPCncls][1] = "";
  fVariableNames[kTPCclusBitFired][0] = "TPC segment"; fVariableNames[kTPCclusBitFired][1] = "";
  fVariableNames[kTPCNclusBitsFired][0] = "No.TPC segments"; fVariableNames[kTPCNclusBitsFired][1] = "";
  fVariableNames[kTPCclustersPerBit][0] = "No.TPC clusters/segment"; fVariableNames[kTPCclustersPerBit][1] = "";
  fVariableNames[kTPCcrossedRows][0] = "No.TPC crossed rows"; fVariableNames[kTPCcrossedRows][1] = "";
  fVariableNames[kTPCnclsIter1][0] = "No.TPC clusters iter.1"; fVariableNames[kTPCnclsIter1][1] = "";
  fVariableNames[kTPCnclsF][0] = "No.TPC findable clusters"; fVariableNames[kTPCnclsF][1] = "";
  fVariableNames[kTPCnclsRatio][0] = "No.TPC clusters/findable"; fVariableNames[kTPCnclsRatio][1] = "";
  fVariableNames[kTPCnclsRatio2][0] = "No.TPC clusters/crossed rows"; fVariableNames[kTPCnclsRatio2][1] = "";
  fVariableNames[kTPCsignal][0] = "TPC dE/dx"; fVariableNames[kTPCsignal][1] = "";  
  fVariableNames[kTPCnSig][0] = "TPC n_{#sigma}^{e}"; fVariableNames[kTPCnSig][1] = "#sigma";  
  fVariableNames[kTPCnSig+1][0] = "TPC n_{#sigma}^{#pi}"; fVariableNames[kTPCnSig+1][1] = "#sigma";
  fVariableNames[kTPCnSig+2][0] = "TPC n_{#sigma}^{K}"; fVariableNames[kTPCnSig+2][1] = "#sigma";
  fVariableNames[kTPCnSig+3][0] = "TPC n_{#sigma}^{p}"; fVariableNames[kTPCnSig+3][1] = "#sigma";
  fVariableNames[kTOFbeta][0] = "TOF #beta"; fVariableNames[kTOFbeta][1] = "";
  fVariableNames[kTOFnSig][0] = "TOF n_{#sigma}^{e}"; fVariableNames[kTOFnSig][1] = "#sigma";  
  fVariableNames[kTOFnSig+1][0] = "TOF n_{#sigma}^{#pi}"; fVariableNames[kTOFnSig+1][1] = "#sigma";
  fVariableNames[kTOFnSig+2][0] = "TOF n_{#sigma}^{K}"; fVariableNames[kTOFnSig+2][1] = "#sigma";
  fVariableNames[kTOFnSig+3][0] = "TOF n_{#sigma}^{p}"; fVariableNames[kTOFnSig+3][1] = "#sigma";
  fVariableNames[kTRDntracklets][0] = "No.TRD tracklets"; fVariableNames[kTRDntracklets][1] = "#sigma";
  fVariableNames[kTRDntrackletsPID][0] = "No.TRD PID tracklets"; fVariableNames[kTRDntrackletsPID][1] = "#sigma";
  fVariableNames[kTRDpidProbabilities][0] = "TRD electron probability"; fVariableNames[kTRDpidProbabilities][1] = "";  
  fVariableNames[kTRDpidProbabilities+1][0] = "TRD pion probability"; fVariableNames[kTRDpidProbabilities+1][1] = "";
  fVariableNames[kEMCALmatchedEnergy][0] = "EMCAL energy"; fVariableNames[kEMCALmatchedEnergy][1] = "GeV";
  fVariableNames[kEMCALmatchedEOverP][0] = "EMCAL E/p"; fVariableNames[kEMCALmatchedEOverP][1] = "";  
  fVariableNames[kEMCALclusterEnergy][0] = "EMCAL cls. energy"; fVariableNames[kEMCALclusterEnergy][1] = "GeV";    
  fVariableNames[kEMCALclusterDx][0] = "EMCAL cls. dx"; fVariableNames[kEMCALclusterDx][1] = "";  
  fVariableNames[kEMCALclusterDz][0] = "EMCAL cls. dz"; fVariableNames[kEMCALclusterDz][1] = "";  
  fVariableNames[kEMCALdetector][0] = "EMCAL detector"; fVariableNames[kEMCALdetector][1] = "";  
  fVariableNames[kTrackingFlag][0] = "Tracking flag"; fVariableNames[kTrackingFlag][1] = "";  
  fVariableNames[kDeltaPhi][0] = "#Delta #varphi"; fVariableNames[kDeltaPhi][1] = "rad.";  
  fVariableNames[kDeltaTheta][0] = "#Delta #theta"; fVariableNames[kDeltaTheta][1] = "rad.";
  fVariableNames[kDeltaEta][0] = "#Delta #eta"; fVariableNames[kDeltaEta][1] = "";
  for(Int_t ibit=0;ibit<9;++ibit) { fVariableNames[kFilterBit+ibit][0] = Form("filter bit %d", ibit); fVariableNames[kFilterBit+ibit][1] = "";}
  fVariableNames[kFilterBitMask768][0] = "filter bit 768"; fVariableNames[kFilterBitMask768][1] = "";

}




AliQnCorrectionsVarManagerTask::AliQnCorrectionsVarManagerTask():
    AliAnalysisTaskSE()
{
  //
  // Default constructor
  //
  SetDefaultVarNames();
}

AliQnCorrectionsVarManagerTask::AliQnCorrectionsVarManagerTask(const char *name):
    AliAnalysisTaskSE(name)
{
  //
  // Default constructor
  //
  SetDefaultVarNames();
}

AliQnCorrectionsVarManagerTask::~AliQnCorrectionsVarManagerTask()
{
  //
  // Default destructor
  //
}

