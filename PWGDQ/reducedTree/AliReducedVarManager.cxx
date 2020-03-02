/*
***********************************************************
  Implementation of AliReducedVarManager
  Contact: iarsene@cern.ch
  2015/04/16
  *********************************************************
*/

#ifndef ALIREDUCEDVARMANAGER_H
#include "AliReducedVarManager.h"
#endif

#include <iostream>
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
#include <fstream>

#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TProfile2D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <THashList.h>

#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedEventPlaneInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedCaloClusterTrackMatcher.h"
#include "AliKFParticle.h"

#define BASEEVENT AliReducedBaseEvent
#define EVENT AliReducedEventInfo
#define TRACK AliReducedTrackInfo
#define PAIR  AliReducedPairInfo
#define EVENTPLANE AliReducedEventPlaneInfo
#define BASETRACK AliReducedBaseTrack
#define CLUSTER AliReducedCaloClusterInfo

ClassImp(AliReducedVarManager)

const Float_t AliReducedVarManager::fgkParticleMass[AliReducedVarManager::kNSpecies] = {
    0.000511,     // electron
    0.13957,      // pion+-
    0.493677,     // kaon+-
    0.938272,     // proton
    0.497614,     // K0
    1.019455      // phi(1020) meson
};

const Float_t AliReducedVarManager::fgkPairMass[AliReducedPairInfo::kNMaxCandidateTypes] = 
{
   0.0, // gamma
   0.497614, //K0s 
   1.11568, //Lambda
   1.11568, //ALambda
   1.019455, // Phi
   3.09691599, //Jpsi
   3.686097,   // psi(2S)
   9.460300, //Upsilon
   1.86962, // D+-
   1.86962, // D+-
   1.86962, // D+-
   1.86962, // D+-
   1.86962, // D+-
   1.86962, // D+-
   1.86484, // D0
   1.86962, // D+-
   1.96850, // Ds
};

const Char_t* AliReducedVarManager::fgkTrackingStatusNames[AliReducedVarManager::kNTrackingStatus] = {
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

const Char_t* AliReducedVarManager::fgkOfflineTriggerNames[64] = {
    "MB/INT1",                                             "INT7",                                                         "MUON",                                    "HighMult/HighMultSPD",                      
    "EMC1",                                                 "CINT5/INT5",                                               "CMUS5/MUSPB/INT7inMUON", "MuonSingleHighPt7/MUSH7/MUSHPB", 
    "MuonLikeLowPt7/MUL7/MuonLikePB", "MuonUnlikeLowPt7/MUU7/MuonUnlikePB", "EMC7/EMC8",                           "MUS7/MuonSingleLowPt7",    
    "PHI1",                                                  "PHI7/PHI8/PHOSPb",                                    "EMCEJE",                                  "EMCEGA",           
    "Central/HighMultV0",                           "SemiCentral",                                              "DG/DG5",                                "ZED",         
    "SPI7/SPI",                                             "INT8",                                                          "MuonSingleLowPt8",                "MuonSingleHighPt8", 
    "MuonLikeLowPt8",                               "MuonUnlikeLowPt8",                                    "MuonUnlikeLowPt0/INT6",       "UserDefined",      
    "TRD",                                                   "MuonCalo/CaloOnly",                                                            "FastOnly",                               "N/A",         
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A",              
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A",
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A",               
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A",      
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A",     
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A",   
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A", 
    "N/A",                                                    "N/A",                                                            "N/A",                                       "N/A"
};


// radii of VZERO channels centers (in cm)
const Double_t AliReducedVarManager::fgkVZEROChannelRadii[64] = {
  6.0567,  6.0567,  6.0567,  6.0567,  6.0567,  6.0567,  6.0567,  6.0567,
  9.6977,  9.6977,  9.6977,  9.6977,  9.6977,  9.6977,  9.6977,  9.6977,
 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504, 15.9504,
 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031, 26.4031,
  5.9347,  5.9347,  5.9347,  5.9347,  5.9347,  5.9347,  5.9347,  5.9347,
 10.685,  10.685,  10.685,  10.685,  10.685,  10.685,  10.685,  10.685,
 18.116,  18.116,  18.116,  18.116,  18.116,  18.116,  18.116,  18.116,
 31.84,   31.84,   31.84,   31.84,   31.84,   31.84,   31.84,   31.84 
};
const Double_t AliReducedVarManager::fgkVZEROAz = 340.0;   // cm
const Double_t AliReducedVarManager::fgkVZEROCz = 90.0;    // cm
const Double_t AliReducedVarManager::fgkVZEROminMult = -0.01;   // minimum VZERO channel multiplicity
const Float_t  AliReducedVarManager::fgkTPCQvecRapGap = 0.8;    // symmetric interval in the middle of the TPC excluded from EP calculation
      Float_t  AliReducedVarManager::fgBeamMomentum = 1380.;   // beam momentum in GeV/c
     
const Double_t AliReducedVarManager::fgkSPDEtaCutsVsVtxZ[20][2] = {
   {-0.5, 1.0}, {-0.6, 1.0}, {-0.8, 1.0}, {-0.9, 1.0}, {-1.0, 1.0},
   {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
   {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
   {-1.0, 1.0}, {-1.0, 0.9}, {-1.0, 0.8}, {-1.0, 0.7}, {-1.0, 0.6}
};
     
Int_t                           AliReducedVarManager::fgCurrentRunNumber = -1;
TString                         AliReducedVarManager::fgVariableNames[AliReducedVarManager::kNVars] = {""};
TString                         AliReducedVarManager::fgVariableUnits[AliReducedVarManager::kNVars] = {""};
AliReducedBaseEvent*            AliReducedVarManager::fgEvent = 0x0;
AliReducedEventPlaneInfo*       AliReducedVarManager::fgEventPlane = 0x0;
Bool_t                          AliReducedVarManager::fgUsedVars[AliReducedVarManager::kNVars] = {kFALSE};
TH2F*                           AliReducedVarManager::fgTPCelectronCentroidMap = 0x0;
TH2F*                           AliReducedVarManager::fgTPCelectronWidthMap = 0x0;
AliReducedVarManager::Variables AliReducedVarManager::fgVarDependencyX = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgVarDependencyY = kNothing;
THnF*                           AliReducedVarManager::fgTPCpidCalibCentroid[3] = {0x0};
THnF*                           AliReducedVarManager::fgTPCpidCalibWidth[3] = {0x0};
THnI*                           AliReducedVarManager::fgTPCpidCalibStatus[3] = {0x0};
AliReducedVarManager::Variables AliReducedVarManager::fgTPCpidCalibVars[4] = {kNothing, kNothing, kNothing, kNothing};
TH1*                            AliReducedVarManager::fgPairEffMap = 0x0;
AliReducedVarManager::Variables AliReducedVarManager::fgEffMapVarDependencyX = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgEffMapVarDependencyY = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgEffMapVarDependencyZ = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgEffMapVarDependencyXCorr = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgEffMapVarDependencyYCorr = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgEffMapVarDependencyZCorr = kNothing;
TH1*                            AliReducedVarManager::fgAssocHadronEffMap = 0x0;
AliReducedVarManager::Variables AliReducedVarManager::fgAssocHadronEffMapVarDependencyX = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgAssocHadronEffMapVarDependencyY = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgAssocHadronEffMapVarDependencyZ = kNothing;
TH1F*                           AliReducedVarManager::fgRunTotalLuminosity = 0x0;
TH1F*                           AliReducedVarManager::fgRunTotalIntensity0 = 0x0;
TH1F*                           AliReducedVarManager::fgRunTotalIntensity1 = 0x0;
TH1I*                           AliReducedVarManager::fgRunLHCFillNumber = 0x0;
TH1I*                           AliReducedVarManager::fgRunDipolePolarity = 0x0;
TH1I*                           AliReducedVarManager::fgRunL3Polarity = 0x0;
TH1I*                           AliReducedVarManager::fgRunTimeStart = 0x0;
TH1I*                           AliReducedVarManager::fgRunTimeEnd = 0x0;
TFile*                          AliReducedVarManager::fgGRPfile = 0x0;
TGraphErrors*                   AliReducedVarManager::fgRunInstLumi = 0x0;
std::vector<Int_t>              AliReducedVarManager::fgRunNumbers;
Int_t                           AliReducedVarManager::fgRunID = -1;
TH1*                            AliReducedVarManager::fgAvgMultVsVtxGlobal[kNMultiplicityEstimators] = {0x0};
TH1*                            AliReducedVarManager::fgAvgMultVsVtxRunwise[kNMultiplicityEstimators] = {0x0};
TH1*                            AliReducedVarManager::fgAvgMultVsRun[kNMultiplicityEstimators] = {0x0};
TH2*                            AliReducedVarManager::fgAvgMultVsVtxAndRun[kNMultiplicityEstimators] = {0x0};
Double_t                        AliReducedVarManager::fgRefMultVsVtxGlobal[kNMultiplicityEstimators] [kNReferenceMultiplicities] = {0.};
Double_t                        AliReducedVarManager::fgRefMultVsVtxRunwise[kNMultiplicityEstimators] [kNReferenceMultiplicities] = {0.};
Double_t                        AliReducedVarManager::fgRefMultVsRun[kNMultiplicityEstimators] [kNReferenceMultiplicities] = {0.};
Double_t                        AliReducedVarManager::fgRefMultVsVtxAndRun[kNMultiplicityEstimators] [kNReferenceMultiplicities] = {0.};

TString                         AliReducedVarManager::fgVZEROCalibrationPath = "";
TProfile2D*                     AliReducedVarManager::fgAvgVZEROChannelMult[64] = {0x0};
TProfile2D*                     AliReducedVarManager::fgVZEROqVecRecentering[4] = {0x0};
TProfile2D*                     AliReducedVarManager::fgTPCqVecRecentering[2] = {0x0};
Bool_t                          AliReducedVarManager::fgOptionCalibrateVZEROqVec = kFALSE;
Bool_t                          AliReducedVarManager::fgOptionRecenterVZEROqVec = kFALSE;
Bool_t                          AliReducedVarManager::fgOptionRecenterTPCqVec = kFALSE;
Bool_t                          AliReducedVarManager::fgOptionEventRes = kFALSE;
TH1F*                           AliReducedVarManager::fgpT_reweight_centrality_0_10=0x0;
//__________________________________________________________________
AliReducedVarManager::AliReducedVarManager() :
  TObject()
{
  //
  // constructor
  //
  SetDefaultVarNames();
}

//__________________________________________________________________
AliReducedVarManager::AliReducedVarManager(const Char_t* name) :
  TObject()
{
  //
  // named constructor
  //
  SetDefaultVarNames();
}

//__________________________________________________________________
AliReducedVarManager::~AliReducedVarManager() {
  //
  // destructor
  //
}

//__________________________________________________________________
void AliReducedVarManager::SetVariableDependencies() {
  //
  // Set as used those variables on which other variables calculation depends
  //
  if(fgUsedVars[kDeltaVtxZ]) {
    fgUsedVars[kVtxZ] = kTRUE;
    fgUsedVars[kVtxZtpc] = kTRUE;
  }
  if(fgUsedVars[kRap] || fgUsedVars[kRapAbs]) {
    fgUsedVars[kMass] = kTRUE;
    fgUsedVars[kP] = kTRUE;
    fgUsedVars[kEta] = kTRUE;
  }
  if(fgUsedVars[kTriggerRap] || fgUsedVars[kTriggerRapAbs]) {
	  fgUsedVars[kMass] = kTRUE;
	  fgUsedVars[kP] = kTRUE;
	  fgUsedVars[kEta] = kTRUE;
  }

  if(fgUsedVars[kEta]) fgUsedVars[kP] = kTRUE;
  
  for(Int_t ih=0; ih<6; ++ih) {
    if(fgUsedVars[kVZEROQvecX+2*6+ih]) {
      fgUsedVars[kVZEROQvecX+0*6+ih] = kTRUE;
      fgUsedVars[kVZEROQvecX+1*6+ih] = kTRUE;
    }
    if(fgUsedVars[kVZEROQvecY+2*6+ih]) {
      fgUsedVars[kVZEROQvecY+0*6+ih] = kTRUE;
      fgUsedVars[kVZEROQvecY+1*6+ih] = kTRUE;
    }
    if(fgUsedVars[kVZERORP+2*6+ih]) {
      fgUsedVars[kVZEROQvecX+2*6+ih] = kTRUE; fgUsedVars[kVZEROQvecY+2*6+ih] = kTRUE;
      fgUsedVars[kVZEROQvecX+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecX+1*6+ih] = kTRUE;
      fgUsedVars[kVZEROQvecY+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecY+1*6+ih] = kTRUE;
    }
    if(fgUsedVars[kVZEROQaQcSP+ih] || fgUsedVars[kVZEROQaQcSPsine+ih]) {
      fgUsedVars[kVZERORP+0*6+ih]    = kTRUE; fgUsedVars[kVZERORP+1*6+ih]    = kTRUE;
      fgUsedVars[kVZEROQvecX+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecX+1*6+ih] = kTRUE;
      fgUsedVars[kVZEROQvecY+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecY+1*6+ih] = kTRUE;
    }
    if(fgUsedVars[kRPXtpcXvzeroa+ih]) {
      fgUsedVars[kTPCQvecX+ih] = kTRUE; fgUsedVars[kVZEROQvecX+ih] = kTRUE;
    }
    if(fgUsedVars[kRPXtpcXvzeroc+ih]) {
      fgUsedVars[kTPCQvecX+ih] = kTRUE; fgUsedVars[kVZEROQvecX+6+ih] = kTRUE;
    }
    if(fgUsedVars[kRPYtpcYvzeroa+ih]) {
      fgUsedVars[kTPCQvecY+ih] = kTRUE; fgUsedVars[kVZEROQvecY+ih] = kTRUE;
    }
    if(fgUsedVars[kRPYtpcYvzeroc+ih]) {
      fgUsedVars[kTPCQvecY+ih] = kTRUE; fgUsedVars[kVZEROQvecY+6+ih] = kTRUE;
    }  
    if(fgUsedVars[kRPXtpcYvzeroa+ih]) {
      fgUsedVars[kTPCQvecX+ih] = kTRUE; fgUsedVars[kVZEROQvecY+ih] = kTRUE;
    }  
    if(fgUsedVars[kRPXtpcYvzeroc+ih]) {
      fgUsedVars[kTPCQvecX+ih] = kTRUE; fgUsedVars[kVZEROQvecY+6+ih] = kTRUE;
    }  
    if(fgUsedVars[kRPYtpcXvzeroa+ih]) {
      fgUsedVars[kTPCQvecY+ih] = kTRUE; fgUsedVars[kVZEROQvecX+ih] = kTRUE;
    }  
    if(fgUsedVars[kRPYtpcXvzeroc+ih]) {
      fgUsedVars[kTPCQvecY+ih] = kTRUE; fgUsedVars[kVZEROQvecX+6+ih] = kTRUE;
    }  
    if(fgUsedVars[kRPdeltaVZEROAtpc+ih]) {
      fgUsedVars[kVZERORP+0*6+ih] = kTRUE; fgUsedVars[kTPCRP+ih] = kTRUE;
    }
    if(fgUsedVars[kRPdeltaVZEROCtpc+ih]) {
      fgUsedVars[kVZERORP+1*6+ih] = kTRUE; fgUsedVars[kTPCRP+ih] = kTRUE;
    }
    if(fgUsedVars[kTPCsubResCos+ih]) {
      fgUsedVars[kTPCRPleft+ih] = kTRUE; fgUsedVars[kTPCRPright+ih] = kTRUE;
    }
    for(Int_t iVZEROside=0; iVZEROside<3; ++iVZEROside) {
      if(fgUsedVars[kVZEROFlowVn+iVZEROside*6+ih] || fgUsedVars[kVZEROFlowSine+iVZEROside*6+ih] ||
	 fgUsedVars[kVZEROuQ+iVZEROside*6+ih] || fgUsedVars[kVZEROuQsine+iVZEROside*6+ih]) {
	fgUsedVars[kPhi] = kTRUE; fgUsedVars[kVZERORP+iVZEROside*6+ih] = kTRUE;
	if(iVZEROside<2 && (fgUsedVars[kVZEROuQ+iVZEROside*6+ih] || fgUsedVars[kVZEROuQsine+iVZEROside*6+ih])) {
	  fgUsedVars[kVZEROQvecX+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecX+1*6+ih] = kTRUE;
          fgUsedVars[kVZEROQvecY+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecY+1*6+ih] = kTRUE;
	}
        if(iVZEROside==2) {
	  fgUsedVars[kVZEROQvecX+2*6+ih] = kTRUE; fgUsedVars[kVZEROQvecY+2*6+ih] = kTRUE;
          fgUsedVars[kVZEROQvecX+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecX+1*6+ih] = kTRUE;
          fgUsedVars[kVZEROQvecY+0*6+ih] = kTRUE; fgUsedVars[kVZEROQvecY+1*6+ih] = kTRUE;
	}
      }
    }
    if(fgUsedVars[kTPCFlowVn+ih] || fgUsedVars[kTPCFlowSine+ih] || fgUsedVars[kTPCuQ+ih] || fgUsedVars[kTPCuQsine+ih]) {
      fgUsedVars[kPhi] = kTRUE;
      fgUsedVars[kTPCQvecXtotal+ih] = kTRUE;
      fgUsedVars[kTPCQvecYtotal+ih] = kTRUE;
    }
   
  } // end loop over harmonics
  for(Int_t ich=0; ich<64; ++ich) {
    if(fgUsedVars[kVZEROflowV2TPC+ich]) {
      fgUsedVars[kVZEROChannelMult+ich] = kTRUE; fgUsedVars[kTPCRP+1] = kTRUE;
    }
  }
  if(fgUsedVars[kPtSquared]) fgUsedVars[kPt]=kTRUE;  
  if(fgUsedVars[kTPCnSigCorrected+kElectron]) fgUsedVars[kTPCnSig+kElectron] = kTRUE; 
  if(fgUsedVars[kTPCnSigCorrected+kPion])     fgUsedVars[kTPCnSig+kPion] = kTRUE; 
  if(fgUsedVars[kTPCnSigCorrected+kProton])   fgUsedVars[kTPCnSig+kProton] = kTRUE;
  if(fgUsedVars[kTPCnSigCorrected+kElectron] || fgUsedVars[kTPCnSigCorrected+kPion] || fgUsedVars[kTPCnSigCorrected+kProton]) {
     fgUsedVars[fgVarDependencyX] = kTRUE; 
     fgUsedVars[fgVarDependencyY] = kTRUE;
     for(Int_t i=0;i<4;++i) fgUsedVars[fgTPCpidCalibVars[i]] = kTRUE;
  }
  
  
  if (fgUsedVars[kTriggerEffTimesAssocHadronEff]) {
    fgUsedVars[kTriggerEff]             = kTRUE;
    fgUsedVars[kAssocHadronEff]         = kTRUE;
  }
  if (fgUsedVars[kOneOverTriggerEffTimesAssocHadronEff]) {
    fgUsedVars[kOneOverTriggerEff]      = kTRUE;
    fgUsedVars[kOneOverAssocHadronEff]  = kTRUE;
  }
  if(fgUsedVars[kPairEff] || fgUsedVars[kOneOverPairEff] || fgUsedVars[kOneOverPairEffSq]){
    fgUsedVars[fgEffMapVarDependencyX] = kTRUE;
    fgUsedVars[fgEffMapVarDependencyY] = kTRUE;
    fgUsedVars[fgEffMapVarDependencyZ] = kTRUE;
  }
  if (fgUsedVars[kTriggerEff] || fgUsedVars[kOneOverTriggerEff]) {
    fgUsedVars[fgEffMapVarDependencyXCorr] = kTRUE;
    fgUsedVars[fgEffMapVarDependencyYCorr] = kTRUE;
    fgUsedVars[fgEffMapVarDependencyZCorr] = kTRUE;
  }
  if(fgUsedVars[kAssocHadronEff] || fgUsedVars[kOneOverAssocHadronEff]){
    fgUsedVars[fgAssocHadronEffMapVarDependencyX] = kTRUE;
    fgUsedVars[fgAssocHadronEffMapVarDependencyY] = kTRUE;
    fgUsedVars[fgAssocHadronEffMapVarDependencyZ] = kTRUE;
  }
  if(fgUsedVars[kNTracksITSoutVsSPDtracklets] || fgUsedVars[kNTracksTPCoutVsSPDtracklets] ||
     fgUsedVars[kNTracksTOFoutVsSPDtracklets] || fgUsedVars[kNTracksTRDoutVsSPDtracklets])
     fgUsedVars[kSPDntracklets] = kTRUE;
  
  if(fgUsedVars[kRapMC] || fgUsedVars[kRapMCAbs]) fgUsedVars[kMassMC] = kTRUE;

  if(fgUsedVars[kPairPhiV]){
    fgUsedVars[kL3Polarity] = kTRUE;
  }
  if(fgUsedVars[kMassDcaPtCorr] ) {
    fgUsedVars[kMass]          = kTRUE;
    fgUsedVars[kPt]            = kTRUE;
    fgUsedVars[kPairDcaXYSqrt] = kTRUE;
  }
  if(fgUsedVars[kOpAngDcaPtCorr] ) {
    fgUsedVars[kPairOpeningAngle] = kTRUE;
    fgUsedVars[kOneOverSqrtPt]    = kTRUE;
    fgUsedVars[kPt]               = kTRUE;
    fgUsedVars[kPairDcaXYSqrt]    = kTRUE;
  }
  if(fgUsedVars[kNTPCclustersFromPileupRelative]) {
    fgUsedVars[kNTPCclusters] = kTRUE;
    fgUsedVars[kNTPCclustersFromPileup] = kTRUE;
  }
  if(fgUsedVars[kNTracksTPCoutFromPileup]) {
    fgUsedVars[kNTracksTPCoutBeforeClean] = kTRUE;
    fgUsedVars[kVZEROTotalMultFromChannels] = kTRUE;
  }
}

//__________________________________________________________________
void AliReducedVarManager::FillEventInfo(Float_t* values) {
  //
  // Fill event information
  //
  FillEventInfo(fgEvent, values, fgEventPlane);
}

void AliReducedVarManager::FillMCEventInfo(AliReducedEventInfo* event, Float_t* values) {
   //
   // fill MC event information
   //
   
}

//__________________________________________________________________
void AliReducedVarManager::FillEventInfo(BASEEVENT* baseEvent, Float_t* values, EVENTPLANE* eventF/*=0x0*/) {
  //
  // fill event wise info
  //
  // Basic event information

  values[kVtxX]                      = baseEvent->Vertex(0);
  values[kVtxY]                       = baseEvent->Vertex(1);
  values[kVtxZ]                      = baseEvent->Vertex(2);
  values[kNVtxContributors]= baseEvent->VertexNContributors(); 
  
  values[kCentVZERO]         = baseEvent->CentralityVZERO();
  values[kCentSPD]              = baseEvent->CentralitySPD();
  values[kCentTPC]              = baseEvent->CentralityTPC();
  values[kCentZDC]             = baseEvent->CentralityZEMvsZDC();
  values[kCentVZEROA]      = baseEvent->CentralityVZEROA();
  values[kCentVZEROC]      = baseEvent->CentralityVZEROC();
  values[kCentZNA]             = baseEvent->CentralityZNA();
  values[kCentQuality]        = baseEvent->CentralityQuality();
  
  values[kNV0total]             = baseEvent->NV0CandidatesTotal();
  values[kNV0selected]       = baseEvent->NV0Candidates();
  values[kNtracksTotal]       = baseEvent->NTracksTotal();
  values[kNtracksSelected] = baseEvent->NTracks();
  
  if(baseEvent->IsA()!=EVENT::Class()) return;
  
  EVENT* event = (EVENT*)baseEvent;
  
  // Update run wise information if available (needed for the first event filled and whenever the run changes)
  if(fgCurrentRunNumber!=baseEvent->RunNo()) {
    fgCurrentRunNumber = baseEvent->RunNo();
    // GRP and LHC information
    if(fgRunTotalLuminosity) values[kTotalLuminosity] = fgRunTotalLuminosity->GetBinContent(fgRunTotalLuminosity->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunTotalIntensity0) values[kBeamIntensity0] = fgRunTotalIntensity0->GetBinContent(fgRunTotalIntensity0->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunTotalIntensity1) values[kBeamIntensity1] = fgRunTotalIntensity1->GetBinContent(fgRunTotalIntensity1->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunLHCFillNumber) values[kLHCFillNumber] = fgRunLHCFillNumber->GetBinContent(fgRunLHCFillNumber->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunDipolePolarity) values[kDipolePolarity] = fgRunDipolePolarity->GetBinContent(fgRunDipolePolarity->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunL3Polarity) values[kL3Polarity] = fgRunL3Polarity->GetBinContent(fgRunL3Polarity->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunTimeStart) values[kRunTimeStart] = fgRunTimeStart->GetBinContent(fgRunTimeStart->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgRunTimeEnd) values[kRunTimeEnd] = fgRunTimeEnd->GetBinContent(fgRunTimeEnd->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber)));
    if(fgGRPfile) {
       fgRunInstLumi = (TGraphErrors*)fgGRPfile->Get(Form("InstLumi_Run%d", fgCurrentRunNumber));
       //cout << "Lumi hist loaded " << fgRunInstLumi << endl;
    }
    
    // VZERO calibration
    if(fgVZEROCalibrationPath.Data()[0]!='\0') {
       cout << "AliReducedVarManager::Info  Attempting to load VZERO calibration and/or recentering histograms from path: " << endl << fgVZEROCalibrationPath.Data() << endl;
      TFile* calibFile = TFile::Open(Form("%s/000%d/dstAnalysisHistograms.root", fgVZEROCalibrationPath.Data(), fgCurrentRunNumber));
      THashList* mainList = (THashList*)calibFile->Get("jpsi2eeHistos");
      THashList* calibList = (THashList*)mainList->FindObject("Event_AfterCuts");
      if(!calibList) {
         cout << "AliReducedVarManager::Info  Cannot open calibration file for run " << fgCurrentRunNumber << endl;
         cout << "                        Will run uncalibrated and not-recentered!" << endl;
         fgOptionCalibrateVZEROqVec = kFALSE;
         fgOptionRecenterVZEROqVec = kFALSE;
         fgOptionRecenterTPCqVec = kFALSE;
      }
      cout << "AliReducedVarManager::Info  Loading VZERO calibration and/or recentering parameters for run " << fgCurrentRunNumber << endl;
      if(fgOptionCalibrateVZEROqVec) {
        for(Int_t iCh=0; iCh<64; ++iCh) {
           
           fgAvgVZEROChannelMult[iCh] = (TProfile2D*)calibList->FindObject(Form("VZEROmult_ch%d_VtxCent_prof", iCh))->Clone(Form("run%d_ch%d", fgCurrentRunNumber, iCh));
           fgAvgVZEROChannelMult[iCh]->SetDirectory(0x0);
        }
      }
      if(fgOptionRecenterVZEROqVec){
         fgVZEROqVecRecentering[0] = (TProfile2D*)calibList->FindObject(Form("QvecX_sideA_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecX_VZEROA", fgCurrentRunNumber));
         fgVZEROqVecRecentering[0]->SetDirectory(0x0);
         fgVZEROqVecRecentering[1] = (TProfile2D*)calibList->FindObject(Form("QvecY_sideA_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecY_VZEROA", fgCurrentRunNumber));
         fgVZEROqVecRecentering[1]->SetDirectory(0x0);
         fgVZEROqVecRecentering[2] = (TProfile2D*)calibList->FindObject(Form("QvecX_sideC_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecX_VZEROC", fgCurrentRunNumber));
         fgVZEROqVecRecentering[2]->SetDirectory(0x0);
         fgVZEROqVecRecentering[3] = (TProfile2D*)calibList->FindObject(Form("QvecY_sideC_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecY_VZEROC", fgCurrentRunNumber));
         fgVZEROqVecRecentering[3]->SetDirectory(0x0);
        
      }
      
      if(fgOptionRecenterTPCqVec){
         fgTPCqVecRecentering[0] = (TProfile2D*)calibList->FindObject(Form("QvecX_TPC_h2_CentV0VtxZ_prof"))->Clone(Form("run%d_QvecX_TPC", fgCurrentRunNumber));
         fgTPCqVecRecentering[0]->SetDirectory(0x0);
         fgTPCqVecRecentering[1] = (TProfile2D*)calibList->FindObject(Form("QvecY_TPC_h2_CentV0VtxZ_prof"))->Clone(Form("run%d_QvecY_TPC", fgCurrentRunNumber));
         fgTPCqVecRecentering[1]->SetDirectory(0x0);
         
      }
      calibFile->Close();
    }

    if(fgUsedVars[kRunID] && fgRunNumbers.size() && fgRunID < 0  ){
      for( fgRunID = 0; fgRunNumbers[ fgRunID ] != fgCurrentRunNumber && fgRunID< (Int_t) fgRunNumbers.size() ; ++fgRunID );
    }
    for( int iEstimator =0 ; iEstimator < kNMultiplicityEstimators ; ++iEstimator ){
      if( fgAvgMultVsVtxAndRun[iEstimator] ){
        Bool_t fillGlobal = !fgAvgMultVsVtxGlobal[iEstimator];
        fgAvgMultVsVtxRunwise  [iEstimator] = fgAvgMultVsVtxAndRun[iEstimator]->ProjectionY( Form("AvgMultVsVtxRunwise%d",iEstimator ), fgRunID, fgRunID );
        if( fillGlobal ){
          fgAvgMultVsVtxGlobal [iEstimator] = fgAvgMultVsVtxAndRun[iEstimator]->ProjectionY( Form("AvgMultVsVtxGlobal%d", iEstimator) );
          fgAvgMultVsVtxGlobal [iEstimator] -> Scale(1. / fgAvgMultVsVtxAndRun[iEstimator]->GetXaxis()->GetNbins());
          fgAvgMultVsRun       [iEstimator] = fgAvgMultVsVtxAndRun[iEstimator]->ProjectionX( Form("AvgMultVsRun%d", iEstimator)  );
          fgAvgMultVsRun       [iEstimator] -> Scale(1. / fgAvgMultVsVtxAndRun[iEstimator]->GetYaxis()->GetNbins());
        }
        for( int iReference = 0; iReference < kNReferenceMultiplicities; ++ iReference  ){
          Double_t refVsVtx, refVsVtxGlobal, refVsRun, refVsVtxAndRun;
          switch ( iReference ){
            case kMaximumMultiplicity :
              refVsVtx = fgAvgMultVsVtxRunwise[iEstimator]->GetMaximum();
              if( fillGlobal ){
                refVsVtxGlobal = fgAvgMultVsVtxGlobal[iEstimator]->GetMaximum();
                refVsRun       = fgAvgMultVsRun[iEstimator]->GetMaximum();
                refVsVtxAndRun = fgAvgMultVsVtxAndRun[iEstimator]->GetMaximum();
              }
              break;
            case kMinimumMultiplicity :
              refVsVtx = fgAvgMultVsVtxRunwise[iEstimator]->GetMinimum();
              if( fillGlobal ){
                refVsVtxGlobal = fgAvgMultVsVtxGlobal[iEstimator]->GetMinimum();
                refVsRun       = fgAvgMultVsRun[iEstimator]->GetMinimum();
                refVsVtxAndRun = fgAvgMultVsVtxAndRun[iEstimator]->GetMinimum();
              }
              break;
            case kMeanMultiplicity :
              refVsVtx = 0.5 * ( fgAvgMultVsVtxRunwise[iEstimator]->GetMaximum() +  fgAvgMultVsVtxRunwise[iEstimator]->GetMinimum() );
              if( fillGlobal ){
                refVsVtxGlobal = 0.5 * ( fgAvgMultVsVtxGlobal[iEstimator]->GetMaximum() + fgAvgMultVsVtxGlobal[iEstimator]->GetMinimum() ) ;
                refVsRun       = 0.5 * ( fgAvgMultVsRun[iEstimator]->GetMaximum() + fgAvgMultVsRun[iEstimator]->GetMinimum() );
                refVsVtxAndRun = 0.5 * ( fgAvgMultVsVtxAndRun[iEstimator]->GetMaximum() + fgAvgMultVsVtxAndRun[iEstimator]->GetMinimum() );
              }
              break;
          }
          fgRefMultVsVtxRunwise  [iEstimator][iReference] = refVsVtx;
          if(fillGlobal){
            fgRefMultVsVtxGlobal [iEstimator][iReference] = refVsVtxGlobal;
            fgRefMultVsRun       [iEstimator][iReference] = refVsRun;
            fgRefMultVsVtxAndRun [iEstimator][iReference] = refVsVtxAndRun;
          }
        }
      }
    }
  }

  values[kRunNo] = fgCurrentRunNumber;
  values[kRunID] = fgRunID;
  
  values[kEventNumberInFile]    = event->EventNumberInFile();
  values[kBC]                   = event->BC();
  values[kTimeStamp]            = event->TimeStamp();
  Int_t timeSOR = 0; Int_t timeEOR = 0;
  if(fgUsedVars[kTimeRelativeSOR] || fgUsedVars[kTimeRelativeSORfraction]) {
     timeSOR = (fgRunTimeStart ? fgRunTimeStart->GetBinContent(fgRunTimeStart->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber))) : 0);
     timeEOR = (fgRunTimeEnd ? fgRunTimeEnd->GetBinContent(fgRunTimeEnd->GetXaxis()->FindBin(Form("%d",fgCurrentRunNumber))) : 0);
  }
  if(fgUsedVars[kTimeRelativeSOR]) {
     values[kTimeRelativeSOR] = Double_t(event->TimeStamp() - timeSOR) / 60.;     // in minutes
  }
  if(fgUsedVars[kTimeRelativeSORfraction] && 
     (values[kRunTimeEnd]-values[kRunTimeStart])>1.)   // the run should be longer than 1 second ... 
     values[kTimeRelativeSORfraction] = (event->TimeStamp() - timeSOR) / (timeEOR - timeSOR);
  if(fgRunInstLumi) {
     values[kInstLumi]          = fgRunInstLumi->Eval(event->TimeStamp()); 
  }
  values[kEventType]            = event->EventType();
  values[kTriggerMask]          = event->TriggerMask();
  values[kINT7Triggered]        = event->TriggerMask() & kINT7 ?1:0;
  values[kTRDTriggeredType]     = event->TRDfired();
  values[kHighMultV0Triggered]  = event->TriggerMask() & kHighMultV0 ?1:0;
  values[kEMCEGATriggered]      = event->TriggerMask() & kEMCEGA ?1:0;
  values[kEMCEGAHighTriggered]  = 0;
  if (values[kEMCEGATriggered]) {
    TString trgClasses = event->TriggerClass();
    if (trgClasses.Contains("EG1") || trgClasses.Contains("DG1")) values[kEMCEGAHighTriggered] = 1;
  }
  values[kIsPhysicsSelection]   = (event->IsPhysicsSelection() ? 1.0 : 0.0);
  values[kIsSPDPileup]          = event->IsSPDPileup();
  values[kIsSPDPileup5]         = event->EventTag(11);
  values[kIsPileupMV]           = event->EventTag(1);
  values[kIsSPDPileupMultBins]  = event->IsSPDPileupMultBins();
  values[kNSPDpileups]          = event->NpileupSPD();
  values[kNTrackPileups]        = event->NpileupTracks();
  values[kIRIntClosestIntMap]   = event->IRIntClosestIntMap(0);
  values[kIRIntClosestIntMap+1] = event->IRIntClosestIntMap(1);
  values[kNPMDtracks]           = event->NPMDtracks();
  values[kNTRDtracks]           = event->NTRDtracks();
  values[kNTRDtracklets]        = event->NTRDtracklets();
  values[kNVtxTPCContributors]  = event->VertexTPCContributors();
  values[kVtxXtpc]              = event->VertexTPC(0);
  values[kVtxYtpc]              = event->VertexTPC(1);
  values[kVtxZtpc]              = event->VertexTPC(2);
  values[kNVtxTPCContributors]  = event->VertexTPCContributors();
  values[kVtxXspd]              = event->VertexSPD(0);
  values[kVtxYspd]              = event->VertexSPD(1);
  values[kVtxZspd]              = event->VertexSPD(2);
  values[kNVtxSPDContributors]  = event->VertexSPDContributors();
  if(fgUsedVars[kDeltaVtxZ]) values[kDeltaVtxZ] = values[kVtxZ] - values[kVtxZtpc];
  if(fgUsedVars[kDeltaVtxZspd]) values[kDeltaVtxZspd] = values[kVtxZ] - values[kVtxZspd];
  values[kTPCpileupZAC]         = event->TPCpileupZ();
  values[kTPCpileupZA]          = event->TPCpileupZ(1);
  values[kTPCpileupZC]          = event->TPCpileupZ(2);
  values[kTPCpileupContributorsAC] = event->TPCpileupContributors();
  values[kTPCpileupContributorsA]  = event->TPCpileupContributors(1);
  values[kTPCpileupContributorsC]  = event->TPCpileupContributors(2);
  values[kTPCpileupZAC2]         = event->TPCpileupZ2();
  values[kTPCpileupZA2]          = event->TPCpileupZ2(1);
  values[kTPCpileupZC2]          = event->TPCpileupZ2(2);
  values[kTPCpileupContributorsAC2] = event->TPCpileupContributors2();
  values[kTPCpileupContributorsA2]  = event->TPCpileupContributors2(1);
  values[kTPCpileupContributorsC2]  = event->TPCpileupContributors2(2);  
  
  for(Int_t iflag=0;iflag<32;++iflag) 
    values[kNTracksPerTrackingStatus+iflag] = event->TracksPerTrackingFlag(iflag);
  values[kNTracksTPCoutBeforeClean] = event->TracksWithTPCout();
  
  // set the fgUsedVars to true as these might have been set to false in the previous event
  fgUsedVars[kNTracksTPCoutVsITSout] = kTRUE;
  fgUsedVars[kNTracksTRDoutVsITSout] = kTRUE;
  fgUsedVars[kNTracksTOFoutVsITSout] = kTRUE;
  fgUsedVars[kNTracksTRDoutVsTPCout] = kTRUE;
  fgUsedVars[kNTracksTOFoutVsTPCout] = kTRUE;
  fgUsedVars[kNTracksTOFoutVsTRDout] = kTRUE;
  if(TMath::Abs(values[kNTracksPerTrackingStatus+kITSout])>0.01) {
    values[kNTracksTPCoutVsITSout] = values[kNTracksPerTrackingStatus+kTPCout]/values[kNTracksPerTrackingStatus+kITSout];
    values[kNTracksTRDoutVsITSout] = values[kNTracksPerTrackingStatus+kTRDout]/values[kNTracksPerTrackingStatus+kITSout]; 
    values[kNTracksTOFoutVsITSout] = values[kNTracksPerTrackingStatus+kTOFout]/values[kNTracksPerTrackingStatus+kITSout];
  }
  else {
     // if these values are undefined, set fgUsedVars as false such that the values are not filled in histograms
     fgUsedVars[kNTracksTPCoutVsITSout] = kFALSE; fgUsedVars[kNTracksTRDoutVsITSout] = kFALSE; fgUsedVars[kNTracksTOFoutVsITSout] = kFALSE;
  }
  
  if(TMath::Abs(values[kNTracksPerTrackingStatus+kTPCout])>0.01) {
    values[kNTracksTRDoutVsTPCout] = values[kNTracksPerTrackingStatus+kTRDout]/values[kNTracksPerTrackingStatus+kTPCout];
    values[kNTracksTOFoutVsTPCout] = values[kNTracksPerTrackingStatus+kTOFout]/values[kNTracksPerTrackingStatus+kTPCout];
  }
  else {
     fgUsedVars[kNTracksTRDoutVsTPCout] = kFALSE; fgUsedVars[kNTracksTOFoutVsTPCout] = kFALSE; 
  }
  
  if(TMath::Abs(values[kNTracksPerTrackingStatus+kTRDout])>0.01)
    values[kNTracksTOFoutVsTRDout] = values[kNTracksPerTrackingStatus+kTOFout]/values[kNTracksPerTrackingStatus+kTRDout];
  else
     fgUsedVars[kNTracksTOFoutVsTRDout] = kFALSE;

  // Multiplicity estimators

  values[ kVZEROATotalMult ] = event->MultVZEROA();
  values[ kVZEROCTotalMult ] = event->MultVZEROC();
  values[ kVZEROTotalMult  ] = event->MultVZERO();
  
  values[ kVZEROATotalMultFromChannels ] = event->MultVZEROA(kTRUE);
  values[ kVZEROCTotalMultFromChannels ] = event->MultVZEROC(kTRUE);
  values[ kVZEROTotalMultFromChannels  ] = event->MultVZERO(kTRUE);
  if(TMath::Abs(values[kVZEROTotalMultFromChannels])>1.0e-6) 
     values[kVZEROTPCoutDiff] = (values[kVZEROTotalMultFromChannels]-values[kNTracksPerTrackingStatus+kTPCout]) / (values[kVZEROTotalMultFromChannels]);
  else
     values[kVZEROTPCoutDiff] = 0.0;
  
  values[ kVZEROACTotalMult ] = event->MultVZEROA() + event->MultVZEROC();

  values[ kSPDntracklets ]   = event->SPDntracklets();
  values[ kSPDntracklets08 ] = 0.;
  values[ kSPDntracklets16 ] = 0.;
  values[ kSPDntrackletsOuterEta ] = 0.;
  values[ kSPDnTracklets10EtaVtxCorr ] = 0.;
  
  
  
  
  for(Int_t ieta=0;ieta<32;++ieta) {
    values[ kSPDntrackletsEtaBin+ieta ] = event->SPDntracklets(ieta);
    if( ieta > 7 && ieta < 24 ) values[ kSPDntracklets08 ] += event->SPDntracklets(ieta);
    if( ieta < 8 || ieta > 23 ) values[ kSPDntrackletsOuterEta ] += event->SPDntracklets(ieta);
  }

  for( Int_t iEstimator = 0; iEstimator < kNMultiplicityEstimators; ++iEstimator){
    Int_t estimator = kMultiplicity + iEstimator;
    if( fgAvgMultVsVtxAndRun[iEstimator] ){
      Int_t vtxBin = fgAvgMultVsVtxAndRun[iEstimator]->GetYaxis()->FindBin( values[kVtxZ] );
      Int_t runBin = fgAvgMultVsVtxAndRun[iEstimator]->GetXaxis()->FindBin( values[kRunID] );
      Double_t multRaw = values[ estimator ];
      for( Int_t iCorrection = 0; iCorrection < kNCorrections; ++iCorrection  ){
        for(Int_t iReference = 0 ; iReference <  kNReferenceMultiplicities; ++iReference ){
          Int_t indexNotSmeared = GetCorrectedMultiplicity( estimator, iCorrection, iReference, kNoSmearing );
          Int_t indexSmeared    = GetCorrectedMultiplicity( estimator, iCorrection, iReference, kPoissonSmearing );
          Double_t multCorr        = multRaw;
          Double_t multCorrSmeared = multRaw;
    // apply vertex and gain loss correction simultaneously
          if( iCorrection == kVertexCorrection2D  ){
            Double_t localAvg = fgAvgMultVsVtxAndRun[iEstimator]->GetBinContent( runBin, vtxBin );
            Double_t refMult  = fgRefMultVsVtxAndRun[iEstimator][iReference];
            multCorr *=  localAvg ?  refMult / localAvg : 1.;
            Double_t deltaM =  localAvg ?  multRaw * ( refMult/localAvg - 1) : 0.;
            multCorrSmeared += (deltaM>0 ? 1. : -1.) * gRandom->Poisson(TMath::Abs(deltaM));
          }
          else{
    // first apply vertex correction
            Double_t localAvgVsVtx , refMultVsVtx  ;
            switch( iCorrection ){
              case kVertexCorrectionGlobal:
              case kVertexCorrectionGlobalGainLoss:
                localAvgVsVtx = fgAvgMultVsVtxGlobal[iEstimator]->GetBinContent( vtxBin );
                refMultVsVtx  = fgRefMultVsVtxGlobal[iEstimator][iReference];
                break;
              case kVertexCorrectionRunwise:
              case kVertexCorrectionRunwiseGainLoss:
                localAvgVsVtx = fgAvgMultVsVtxRunwise[iEstimator]->GetBinContent( vtxBin );
                refMultVsVtx  = fgRefMultVsVtxRunwise[iEstimator][iReference];
                break;
              default:
                localAvgVsVtx = 0.;
                refMultVsVtx  = 0.;
                break;  
            }
            multCorr        *= localAvgVsVtx > 0. ? refMultVsVtx / localAvgVsVtx : 1.;
            Double_t deltaM  = localAvgVsVtx > 0. ? multRaw  * ( refMultVsVtx/localAvgVsVtx - 1) : 0.;
            multCorrSmeared += (deltaM>0 ? 1. : -1.) * gRandom->Poisson(TMath::Abs(deltaM));
    // then apply gain loss correction
            if( iCorrection == kVertexCorrectionGlobalGainLoss || 
                iCorrection == kVertexCorrectionRunwiseGainLoss || 
                iCorrection == kGainLossCorrection   ){
              Double_t localAvgVsRun = fgAvgMultVsRun[iEstimator]->GetBinContent( runBin );
              Double_t refMultVsRun  = fgRefMultVsRun[iEstimator][iReference];
              multCorr        *= localAvgVsRun ? refMultVsRun / localAvgVsRun : 1.;
              deltaM           = localAvgVsRun ? multCorrSmeared  * ( refMultVsRun/localAvgVsRun - 1) : 0;
              multCorrSmeared += (deltaM>0 ? 1. : -1.) * gRandom->Poisson(TMath::Abs(deltaM));
            }
          }
          values[ indexNotSmeared ] = multCorr;
          values[ indexSmeared ]    = multCorrSmeared;
          fgUsedVars [indexNotSmeared] = kTRUE;
          fgUsedVars [indexSmeared] = kTRUE;
        }
      }
    }
    else if( ( estimator == kVZEROACTotalMult && fgAvgMultVsVtxAndRun[kVZEROATotalMult-kMultiplicity] && fgAvgMultVsVtxAndRun[kVZEROCTotalMult-kMultiplicity]  )  
      || (estimator == kSPDnTracklets10EtaVtxCorr && fgAvgMultVsVtxAndRun[kSPDntrackletsEtaBin-kMultiplicity] ) ){
      for( Int_t iCorrection = 0; iCorrection < kNCorrections; ++iCorrection  ){
        for(Int_t iReference = 0 ; iReference <  kNReferenceMultiplicities; ++iReference ){
          Int_t indexNotSmeared = GetCorrectedMultiplicity( estimator, iCorrection, iReference, kNoSmearing );
          Int_t indexSmeared    = GetCorrectedMultiplicity( estimator, iCorrection, iReference, kPoissonSmearing );
          values[indexNotSmeared] = 0.;
          values[indexSmeared] = 0.;
          if( estimator == kSPDnTracklets10EtaVtxCorr ){
            for( Int_t ieta=6; ieta<26; ++ieta ) {
              Int_t vtxBin = fgAvgMultVsVtxAndRun[kSPDntrackletsEtaBin+ieta-kMultiplicity]->GetYaxis()->FindBin( values[kVtxZ] );
              if( fgAvgMultVsVtxGlobal[kSPDntrackletsEtaBin+ieta-kMultiplicity]->GetBinContent( vtxBin ) > .3 ){
                Int_t indexBinNotSmeared = GetCorrectedMultiplicity( kSPDntrackletsEtaBin+ieta, iCorrection, iReference, kNoSmearing );
                Int_t indexBinSmeared    = GetCorrectedMultiplicity( kSPDntrackletsEtaBin+ieta, iCorrection, iReference, kPoissonSmearing );
                if( fgUsedVars[indexBinNotSmeared]) values[ indexNotSmeared ] += values[ indexBinNotSmeared ];
                if( fgUsedVars[indexBinSmeared]) values[ indexSmeared ] += values[ indexBinSmeared ];
              }
            }
          }
          else{
            Int_t indexAnotSmeared = GetCorrectedMultiplicity( kVZEROATotalMult, iCorrection, iReference, kNoSmearing );
            Int_t indexCnotSmeared = GetCorrectedMultiplicity( kVZEROCTotalMult, iCorrection, iReference, kNoSmearing );

            Int_t indexAsmeared = GetCorrectedMultiplicity( kVZEROATotalMult, iCorrection, iReference, kPoissonSmearing );
            Int_t indexCsmeared = GetCorrectedMultiplicity( kVZEROCTotalMult, iCorrection, iReference, kPoissonSmearing );

            values[ indexNotSmeared ] = values[ indexAnotSmeared ] + values[ indexCnotSmeared ];
            values[ indexSmeared ]    = values[ indexAsmeared ]    + values[ indexCsmeared ];
          }
        }
      }
    }
  }

  fgUsedVars[kNTracksITSoutVsSPDtracklets] = kTRUE;  
  fgUsedVars[kNTracksTPCoutVsSPDtracklets] = kTRUE;
  fgUsedVars[kNTracksTRDoutVsSPDtracklets] = kTRUE;
  fgUsedVars[kNTracksTOFoutVsSPDtracklets] = kTRUE;
  if(values[kSPDntracklets]>0.01) {
    values[kNTracksITSoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kITSout] / values[kSPDntracklets];
    values[kNTracksTPCoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kTPCout] / values[kSPDntracklets];
    values[kNTracksTRDoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kTRDout] / values[kSPDntracklets];
    values[kNTracksTOFoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kTOFout] / values[kSPDntracklets];
  }
  else {
     fgUsedVars[kNTracksITSoutVsSPDtracklets] = kFALSE;  
     fgUsedVars[kNTracksTPCoutVsSPDtracklets] = kFALSE;
     fgUsedVars[kNTracksTRDoutVsSPDtracklets] = kFALSE;
     fgUsedVars[kNTracksTOFoutVsSPDtracklets] = kFALSE;
  }
    
  values[kNCaloClusters]   = event->GetNCaloClusters();
  values[kNTPCclusters]    = event->NTPCClusters();

  
  
  for(Int_t i=0;i<2;++i) values[kSPDFiredChips+i] = event->SPDFiredChips(i+1);
  for(Int_t i=0;i<6;++i) values[kITSnClusters+i] = event->ITSClusters(i+1);
  values[kSPDnSingleClusters] = event->SPDnSingleClusters();
  if(fgUsedVars[kSDDandSSDclusters]) {
     values[kSDDandSSDclusters] = 0.0;
     for(Int_t i=2;i<6;++i) values[kSDDandSSDclusters] += event->ITSClusters(i+1);  
  }
  
  //VZERO detector information
  fgUsedVars[kNTracksTPCoutVsVZEROTotalMult] = kTRUE;
  if(values[kVZEROTotalMult]>1.0e-5)
     values[kNTracksTPCoutVsVZEROTotalMult] = values[kNTracksPerTrackingStatus+kTPCout] / values[kVZEROTotalMult];
  else
     fgUsedVars[kNTracksTPCoutVsVZEROTotalMult] = kFALSE;
  
  values[kVZEROAemptyChannels] = 0;
  values[kVZEROCemptyChannels] = 0;
  for(Int_t ich=0;ich<64;++ich) fgUsedVars[kVZEROChannelMult+ich] = kTRUE; 
  Float_t theta=0.0;
  for(Int_t ich=0;ich<64;++ich) {
    if(fgUsedVars[kVZEROChannelMult+ich]) {
      values[kVZEROChannelMult+ich] = event->MultChannelVZERO(ich);
      if(values[kVZEROChannelMult+ich]<fgkVZEROminMult) {
        fgUsedVars[kVZEROChannelMult+ich] = kFALSE;   // will not be filled in histograms by the histogram manager
        if(ich<32) values[kVZEROCemptyChannels] += 1;
        else values[kVZEROAemptyChannels] += 1;
      }
    }
    if(fgUsedVars[kVZEROChannelEta+ich]) {
      if(ich<32) theta = TMath::ATan(fgkVZEROChannelRadii[ich]/(fgkVZEROCz-values[kVtxZ]));
      else theta = TMath::Pi()-TMath::ATan(fgkVZEROChannelRadii[ich]/(fgkVZEROAz-values[kVtxZ]));
      values[kVZEROChannelEta+ich] = -1.0*TMath::Log(TMath::Tan(theta/2.0));
    }
  }
  
  if(values[kVZEROTotalMultFromChannels]>0.0) 
     values[kNTracksTPCoutFromPileup] = values[kNTracksTPCoutBeforeClean] - (-3.2+TMath::Sqrt(3.2*3.2+4.0*1.6e-5*values[kVZEROTotalMultFromChannels]))/(2.0*1.6e-5);
  
  Float_t tpcClustersExpectationWOpileup = 0.001;
  if(fgUsedVars[kNTPCclustersFromPileup] && values[kVZEROTotalMultFromChannels]>0.0) {
     tpcClustersExpectationWOpileup = (-0.0132+TMath::Sqrt(0.0132*0.0132-4.0*(-200.0-values[kVZEROTotalMultFromChannels])*1.4e-9))/2.0/1.4e-9;
     values[kNTPCclustersFromPileup] = values[kNTPCclusters] - tpcClustersExpectationWOpileup;
  }
  else 
     values[kNTPCclustersFromPileup] = 0.0;
  
  if(fgUsedVars[kNTPCclustersFromPileupRelative] && values[kNTPCclusters]>0.0)
     values[kNTPCclustersFromPileupRelative] = values[kNTPCclustersFromPileup] / tpcClustersExpectationWOpileup;
  
  if(fgUsedVars[kVZEROQvecX+0*6+1] || fgUsedVars[kVZEROQvecY+0*6+1] || fgUsedVars[kVZERORP+0*6+1]) {
    Double_t qvecVZEROA[EVENTPLANE::fgkNMaxHarmonics][2] = {{0.0}};
    Double_t qvecVZEROC[EVENTPLANE::fgkNMaxHarmonics][2] = {{0.0}};
    if(fgOptionCalibrateVZEROqVec && fgAvgVZEROChannelMult[0]) {
       Float_t calibVZEROMult[64] = {0.};
       Float_t refMult=0;
      for(Int_t ich=0;ich<64;++ich) fgUsedVars[kVZEROChannelMultCalib+ich] = kTRUE; 
       
      for(Int_t iCh=0; iCh<64; ++iCh) {
         if(event->MultChannelVZERO(iCh)>=fgkVZEROminMult) {
                         
            Float_t avMult = fgAvgVZEROChannelMult[iCh]->GetBinContent(fgAvgVZEROChannelMult[iCh]->FindBin(event->Vertex(2), event->CentralitySPD()));
            Int_t refCh = iCh-(iCh%8);
            if(iCh==refCh)
                refMult=fgAvgVZEROChannelMult[iCh]->GetBinContent(fgAvgVZEROChannelMult[iCh]->GetXaxis()->FindBin(0.0),fgAvgVZEROChannelMult[iCh]->GetYaxis()->FindBin(event->CentralitySPD()));
            
            calibVZEROMult[iCh] = event->MultChannelVZERO(iCh) / (avMult>1.0e-6 ? avMult : 1.0)*refMult;
            values[kVZEROChannelMultCalib+iCh]=calibVZEROMult[iCh];
            
          //  cout<<"V0 channel"<<" "<<iCh<<" "<<"Reference Channel "<<refCh<<" "<<"Reference multiplicity"<<refMult<<"Avg mult"<<avMult<<" "<<values[kVZEROChannelMult+iCh]<<" "<<values[kVZEROChannelMultCalib+iCh]<<endl;
            
         }
         else
             fgUsedVars[kVZEROChannelMultCalib+iCh] = kFALSE; // will not be filled in histograms by the histogram manager
      }
      event->GetVZEROQvector(qvecVZEROA, EVENTPLANE::kVZEROA, calibVZEROMult);
      event->GetVZEROQvector(qvecVZEROC, EVENTPLANE::kVZEROC, calibVZEROMult);
      
     
    }
    else {
      event->GetVZEROQvector(qvecVZEROA, EVENTPLANE::kVZEROA);
      event->GetVZEROQvector(qvecVZEROC, EVENTPLANE::kVZEROC);
    }
    if(fgOptionRecenterVZEROqVec && fgVZEROqVecRecentering[0]) {
         Float_t recenterOffset = fgVZEROqVecRecentering[0]->GetBinContent(fgVZEROqVecRecentering[0]->FindBin(event->CentralitySPD(), event->Vertex(2)));
         Float_t widthEqVZERO = fgVZEROqVecRecentering[0]->GetBinError(fgVZEROqVecRecentering[0]->FindBin(event->CentralitySPD(), event->Vertex(2)));

         qvecVZEROA[1][0] -= recenterOffset;
           if(widthEqVZERO >0.0)
             qvecVZEROA[1][0] /=widthEqVZERO;
           else
             qvecVZEROA[1][0]=0;
                
        recenterOffset = fgVZEROqVecRecentering[1]->GetBinContent(fgVZEROqVecRecentering[1]->FindBin(event->CentralitySPD(), event->Vertex(2)));
        widthEqVZERO = fgVZEROqVecRecentering[1]->GetBinError(fgVZEROqVecRecentering[1]->FindBin(event->CentralitySPD(), event->Vertex(2)));
        qvecVZEROA[1][1] -= recenterOffset;
           if(widthEqVZERO >0.0)
              qvecVZEROA[1][1] /=widthEqVZERO;
           else
            qvecVZEROA[1][1]=0;
          
        recenterOffset = fgVZEROqVecRecentering[2]->GetBinContent(fgVZEROqVecRecentering[2]->FindBin(event->CentralitySPD(), event->Vertex(2)));
        widthEqVZERO = fgVZEROqVecRecentering[2]->GetBinError(fgVZEROqVecRecentering[2]->FindBin(event->CentralitySPD(), event->Vertex(2)));
        qvecVZEROC[1][0] -= recenterOffset;
           if(widthEqVZERO >0.0)
              qvecVZEROC[1][0] /=widthEqVZERO;
           else
              qvecVZEROC[1][0]=0;
       
        recenterOffset = fgVZEROqVecRecentering[3]->GetBinContent(fgVZEROqVecRecentering[3]->FindBin(event->CentralitySPD(), event->Vertex(2)));
        widthEqVZERO = fgVZEROqVecRecentering[3]->GetBinError(fgVZEROqVecRecentering[3]->FindBin(event->CentralitySPD(), event->Vertex(2)));
        qvecVZEROC[1][1] -= recenterOffset;
           if(widthEqVZERO >0.0)
              qvecVZEROC[1][1] /=widthEqVZERO;
           else
              qvecVZEROC[1][1]=0;
     
    }
    for(Int_t ih=1; ih<2; ++ih) {
       // VZERO event plane variables
       values[kVZEROQvecX+0*6+ih] = qvecVZEROA[ih][0];
       values[kVZEROQvecY+0*6+ih] = qvecVZEROA[ih][1];
       values[kVZEROQvecX+1*6+ih] = qvecVZEROC[ih][0];
       values[kVZEROQvecY+1*6+ih] = qvecVZEROC[ih][1];
       values[kVZERORP+0*6+ih] = TMath::ATan2(qvecVZEROA[ih][1], qvecVZEROA[ih][0])/Double_t(ih+1);
       values[kVZERORP+1*6+ih] = TMath::ATan2(qvecVZEROC[ih][1], qvecVZEROC[ih][0])/Double_t(ih+1);
       values[kVZEROQvecX+2*6+ih] = qvecVZEROA[ih][0] + qvecVZEROC[ih][0];
       values[kVZEROQvecY+2*6+ih] = qvecVZEROA[ih][1] + qvecVZEROC[ih][1];
       values[kVZERORP   +2*6+ih] = TMath::ATan2(values[kVZEROQvecY+2*6+ih], values[kVZEROQvecX+2*6+ih])/Double_t(ih+1);
     
       if(fgUsedVars[kVZEROQaQcSP+ih]) {
          values[kVZEROQaQcSP+ih] = TMath::Cos((ih+1)*(values[kVZERORP+0*6+ih]-values[kVZERORP+1*6+ih]));
          values[kVZEROQaQcSP+ih] *= TMath::Sqrt(values[kVZEROQvecX+0*6+ih]*values[kVZEROQvecX+0*6+ih]+
          values[kVZEROQvecY+0*6+ih]*values[kVZEROQvecY+0*6+ih]);
          values[kVZEROQaQcSP+ih] *= TMath::Sqrt(values[kVZEROQvecX+1*6+ih]*values[kVZEROQvecX+1*6+ih]+
          values[kVZEROQvecY+1*6+ih]*values[kVZEROQvecY+1*6+ih]);
       }
       values[kVZEROQaQcSPsine+ih] = TMath::Sin((ih+1)*(values[kVZERORP+0*6+ih]-values[kVZERORP+1*6+ih]));
       values[kVZEROQaQcSPsine+ih] *= TMath::Sqrt(values[kVZEROQvecX+0*6+ih]*values[kVZEROQvecX+0*6+ih]+
       values[kVZEROQvecY+0*6+ih]*values[kVZEROQvecY+0*6+ih]);
       values[kVZEROQaQcSPsine+ih] *= TMath::Sqrt(values[kVZEROQvecX+1*6+ih]*values[kVZEROQvecX+1*6+ih]+
       values[kVZEROQvecY+1*6+ih]*values[kVZEROQvecY+1*6+ih]);
       values[kVZERORP   +2*6+ih] = TMath::ATan2(values[kVZEROQvecY+2*6+ih],values[kVZEROQvecX+2*6+ih])/Double_t(ih+1);
       // cos (n*(psi_A-psi_C))
       if(fgUsedVars[kVZERORPres + ih]) {
          values[kVZERORPres + ih] = DeltaPhi(values[kVZERORP+0*6+ih], values[kVZERORP+1*6+ih]);
          values[kVZERORPres + ih] = TMath::Cos(values[kVZERORPres + ih]*(ih+1));
       }
       // Qx,Qy correlations for VZERO
       if(fgUsedVars[kVZEROXaXc+ih]) 
          values[kVZEROXaXc+ih] = qvecVZEROA[ih][0]*qvecVZEROC[ih][0];
       if(fgUsedVars[kVZEROXaYa+ih]) 
          values[kVZEROXaYa+ih] = qvecVZEROA[ih][0]*qvecVZEROA[ih][1];
       if(fgUsedVars[kVZEROXaYc+ih]) 
          values[kVZEROXaYc+ih] = qvecVZEROA[ih][0]*qvecVZEROC[ih][1];
       if(fgUsedVars[kVZEROYaXc+ih]) 
          values[kVZEROYaXc+ih] = qvecVZEROA[ih][1]*qvecVZEROC[ih][0];
       if(fgUsedVars[kVZEROYaYc+ih]) 
          values[kVZEROYaYc+ih] = qvecVZEROA[ih][1]*qvecVZEROC[ih][1];
       if(fgUsedVars[kVZEROXcYc+ih]) 
          values[kVZEROXcYc+ih] = qvecVZEROC[ih][0]*qvecVZEROC[ih][1];
       // Psi_A - Psi_C
       if(fgUsedVars[kVZEROdeltaRPac+ih])
          values[kVZEROdeltaRPac+ih] = DeltaPhi(values[kVZERORP+0*6+ih], values[kVZERORP+1*6+ih]);
    }    // end loop over harmonics
  }
  
  // Get the TPC event plane in case it was written in the trees 
  for(Int_t ih=0; ih<3; ++ih) {
     if(event->GetEventPlaneStatus(EVENTPLANE::kTPC,ih+1)!=EVENTPLANE::kUnset) {
        values[kTPCQvecXtree+ih] = event->GetQx(EVENTPLANE::kTPC,ih+1);
        values[kTPCQvecYtree+ih] = event->GetQy(EVENTPLANE::kTPC,ih+1);
        values[kTPCRPtree+ih] = event->GetEventPlane(EVENTPLANE::kTPC,ih+1);
        
          //TPC Q vector recentering       
          if(fgOptionRecenterTPCqVec && fgTPCqVecRecentering[0] && ih==1) {
            Float_t recenterOffsetTPC = fgTPCqVecRecentering[0]->GetBinContent(fgTPCqVecRecentering[0]->FindBin(event->CentralityVZERO(), event->Vertex(2)));
            Double_t widthEqTPC = fgTPCqVecRecentering[0]->GetBinError(fgTPCqVecRecentering[0]->FindBin(event->CentralityVZERO(), event->Vertex(2)));
            values[kTPCQvecXtree+1] -= recenterOffsetTPC;

                if(widthEqTPC >0.0)
                   values[kTPCQvecXtree+1] /=widthEqTPC;
                else
                   values[kTPCQvecXtree+1]=0;
                
            recenterOffsetTPC = fgTPCqVecRecentering[1]->GetBinContent(fgTPCqVecRecentering[1]->FindBin(event->CentralityVZERO(), event->Vertex(2)));
            widthEqTPC = fgTPCqVecRecentering[1]->GetBinError(fgTPCqVecRecentering[1]->FindBin(event->CentralityVZERO(), event->Vertex(2)));
            values[kTPCQvecYtree+1] -= recenterOffsetTPC;
                
                if(widthEqTPC >0.0)
                   values[kTPCQvecYtree+1] /=widthEqTPC;
                else
                   values[kTPCQvecYtree+1]=0;
                 
            values[kTPCRPtree+ih]=TMath::ATan2(values[kTPCQvecYtree+1],values[kTPCQvecXtree+1])/Double_t(ih+1);
      
        }
 
     }
     if(event->GetEventPlaneStatus(EVENTPLANE::kTPCptWeights,ih+1)!=EVENTPLANE::kUnset) {
        values[kTPCQvecXptWeightsTree+ih] = event->GetQx(EVENTPLANE::kTPCptWeights,ih+1);
        values[kTPCQvecYptWeightsTree+ih] = event->GetQy(EVENTPLANE::kTPCptWeights,ih+1);
        values[kTPCRPptWeightsTree+ih] = event->GetEventPlane(EVENTPLANE::kTPCptWeights,ih+1);
     }
     if(event->GetEventPlaneStatus(EVENTPLANE::kTPCpos,ih+1)!=EVENTPLANE::kUnset) {
        values[kTPCQvecXposTree+ih] = event->GetQx(EVENTPLANE::kTPCpos,ih+1);
        values[kTPCQvecYposTree+ih] = event->GetQy(EVENTPLANE::kTPCpos,ih+1);
        values[kTPCRPposTree+ih] = event->GetEventPlane(EVENTPLANE::kTPCpos,ih+1);
     }
     if(event->GetEventPlaneStatus(EVENTPLANE::kTPCneg,ih+1)!=EVENTPLANE::kUnset) {
        values[kTPCQvecXnegTree+ih] = event->GetQx(EVENTPLANE::kTPCneg,ih+1);
        values[kTPCQvecYnegTree+ih] = event->GetQy(EVENTPLANE::kTPCneg,ih+1);
        values[kTPCRPnegTree+ih] = event->GetEventPlane(EVENTPLANE::kTPCneg,ih+1);
     }
     
      // TPC VZERO Q-vector correlations
     if(fgUsedVars[kRPXtpcXvzeroa+ih]) 
	values[kRPXtpcXvzeroa+ih] = values[kTPCQvecXtree+ih]*values[kVZEROQvecX+ih];
     if(fgUsedVars[kRPXtpcXvzeroc+ih]) 
	values[kRPXtpcXvzeroc+ih] = values[kTPCQvecXtree+ih]*values[kVZEROQvecX+6+ih];
     if(fgUsedVars[kRPYtpcYvzeroa+ih]) 
	values[kRPYtpcYvzeroa+ih] = values[kTPCQvecYtree+ih]*values[kVZEROQvecY+ih];
     if(fgUsedVars[kRPYtpcYvzeroc+ih]) 
	values[kRPYtpcYvzeroc+ih] = values[kTPCQvecYtree+ih]*values[kVZEROQvecY+6+ih];
     if(fgUsedVars[kRPXtpcYvzeroa+ih]) 
	values[kRPXtpcYvzeroa+ih] = values[kTPCQvecXtree+ih]*values[kVZEROQvecY+ih];
     if(fgUsedVars[kRPXtpcYvzeroc+ih]) 
	values[kRPXtpcYvzeroc+ih] = values[kTPCQvecXtree+ih]*values[kVZEROQvecY+6+ih];
     if(fgUsedVars[kRPYtpcXvzeroa+ih]) 
	values[kRPYtpcXvzeroa+ih] = values[kTPCQvecYtree+ih]*values[kVZEROQvecX+ih];
     if(fgUsedVars[kRPYtpcXvzeroc+ih]) 
	values[kRPYtpcXvzeroc+ih] = values[kTPCQvecYtree+ih]*values[kVZEROQvecX+6+ih];
      // Psi_TPC - Psi_VZERO A/C      
     if(fgUsedVars[kRPdeltaVZEROAtpc+ih]) 
	values[kRPdeltaVZEROAtpc+ih] = DeltaPhi(values[kVZERORP+0*6+ih], values[kTPCRPtree+ih]);
     if(fgUsedVars[kRPdeltaVZEROCtpc+ih])
        values[kRPdeltaVZEROCtpc+ih] = DeltaPhi(values[kVZERORP+1*6+ih], values[kTPCRPtree+ih]);
     
     
     // cos(n(EPtpc-EPvzero A/C))
     for(Int_t iVZEROside=0; iVZEROside<2; ++iVZEROside) {
          if(fgUsedVars[kTPCRPres+iVZEROside*6+ih]) {
	  values[kTPCRPres+iVZEROside*6+ih] = DeltaPhi(values[kTPCRPtree+ih], values[kVZERORP+iVZEROside*6+ih]);
          values[kTPCRPres+iVZEROside*6+ih] = TMath::Cos(values[kTPCRPres+iVZEROside*6+ih]*(ih+1));
	}
      }
     // cout<<values[kCentVZERO]<<" "<<values[kCentSPD]<<" "<<values[kVtxZ]<<" "<<values[kTPCRPres+1*6+ih]<<" "<<values[kTPCRPres+0*6+ih]<<" "<<values[kVZERORPres+ih]<<endl;
//      cout<<values[kTPCRPres+0*6+ih]<<endl;
//      cout<<values[kTPCRPres+1*6+ih]<<endl;
//      cout<<values[kVZERORPres+ih]<<endl;
      //resolution of V0A, V0C or TPC as reference detector
      if(fgOptionEventRes && (fgUsedVars[kVZEROARPres+ih]||fgUsedVars[kVZEROCRPres+ih]||fgUsedVars[kVZEROTPCRPres+ih])){
         
         if(values[kTPCRPres+1*6+ih]>1.0e-7 && values[kTPCRPres+0*6+ih]>1.0e-7 && values[kVZERORPres+ih]>1.0e-7){
    
          values[kVZEROARPres+ih] = TMath::Sqrt(values[kTPCRPres+1*6+ih]/(values[kVZERORPres + ih]*values[kTPCRPres+0*6+ih]));
          values[kVZEROCRPres+ih] = TMath::Sqrt(values[kTPCRPres+0*6+ih]/(values[kVZERORPres + ih]*values[kTPCRPres+1*6+ih]));
          values[kVZEROTPCRPres+ih] = TMath::Sqrt(values[kVZERORPres+ih]/(values[kTPCRPres+0*6 + ih]*values[kTPCRPres+1*6+ih]));
        }
        else{
          values[kVZEROARPres+ih]=0;
          values[kVZEROCRPres+ih]=0;
          values[kVZEROTPCRPres+ih]=0;
        }
    }//end if fgOptionEventRes
     
  }// end loop over harmonics
  
 
     
    
   
  
  if(eventF) {
   for(Int_t ih=0; ih<6; ++ih) {
     // VZERO event plane variables
     values[kVZEROQvecX+2*6+ih] = 0.0;
     values[kVZEROQvecY+2*6+ih] = 0.0;
     values[kVZERORP   +2*6+ih] = 0.0;
     for(Int_t iVZEROside=0; iVZEROside<2; ++iVZEROside) {
       values[kVZEROQvecX+iVZEROside*6+ih] = eventF->Qx(EVENTPLANE::kVZEROA+iVZEROside, ih+1);
       values[kVZEROQvecY+iVZEROside*6+ih] = eventF->Qy(EVENTPLANE::kVZEROA+iVZEROside, ih+1);
       if(fgUsedVars[kVZERORP+iVZEROside*6+ih]) 
        values[kVZERORP+iVZEROside*6+ih] = eventF->EventPlane(EVENTPLANE::kVZEROA+iVZEROside, ih+1);
	if(fgUsedVars[kVZEROQvecX+2*6+ih])
	  values[kVZEROQvecX+2*6+ih] += values[kVZEROQvecX+iVZEROside*6+ih];
	if(fgUsedVars[kVZEROQvecY+2*6+ih])
	  values[kVZEROQvecY+2*6+ih] += values[kVZEROQvecY+iVZEROside*6+ih];
	// cos(n(EPtpc-EPvzero A/C))	
        if(fgUsedVars[kTPCRPres+iVZEROside*6+ih]) {
	  values[kTPCRPres+iVZEROside*6+ih] = DeltaPhi(eventF->EventPlane(EVENTPLANE::kTPC, ih+1), eventF->EventPlane(EVENTPLANE::kVZEROA+iVZEROside, ih+1));
          values[kTPCRPres+iVZEROside*6+ih] = TMath::Cos(values[kTPCRPres+iVZEROside*6+ih]*(ih+1));
	}
      }
      
      if(fgUsedVars[kVZEROQaQcSP+ih]) {
        values[kVZEROQaQcSP+ih] = TMath::Cos((ih+1)*(values[kVZERORP+0*6+ih]-values[kVZERORP+1*6+ih]));
        values[kVZEROQaQcSP+ih] *= TMath::Sqrt(values[kVZEROQvecX+0*6+ih]*values[kVZEROQvecX+0*6+ih]+
                                               values[kVZEROQvecY+0*6+ih]*values[kVZEROQvecY+0*6+ih]);
        values[kVZEROQaQcSP+ih] *= TMath::Sqrt(values[kVZEROQvecX+1*6+ih]*values[kVZEROQvecX+1*6+ih]+
                                               values[kVZEROQvecY+1*6+ih]*values[kVZEROQvecY+1*6+ih]);
      }
      values[kVZEROQaQcSPsine+ih] = TMath::Sin((ih+1)*(values[kVZERORP+0*6+ih]-values[kVZERORP+1*6+ih]));
      values[kVZEROQaQcSPsine+ih] *= TMath::Sqrt(values[kVZEROQvecX+0*6+ih]*values[kVZEROQvecX+0*6+ih]+
                                             values[kVZEROQvecY+0*6+ih]*values[kVZEROQvecY+0*6+ih]);
      values[kVZEROQaQcSPsine+ih] *= TMath::Sqrt(values[kVZEROQvecX+1*6+ih]*values[kVZEROQvecX+1*6+ih]+
                                             values[kVZEROQvecY+1*6+ih]*values[kVZEROQvecY+1*6+ih]);
      values[kVZERORP   +2*6+ih] = TMath::ATan2(values[kVZEROQvecY+2*6+ih],values[kVZEROQvecX+2*6+ih])/Double_t(ih+1);
      // cos (n*(psi_A-psi_C))
      if(fgUsedVars[kVZERORPres + ih]) {
	values[kVZERORPres + ih] = DeltaPhi(eventF->EventPlane(EVENTPLANE::kVZEROA, ih+1), 
					    eventF->EventPlane(EVENTPLANE::kVZEROC, ih+1));
        values[kVZERORPres + ih] = TMath::Cos(values[kVZERORPres + ih]*(ih+1));
      }
      // Qx,Qy correlations for VZERO
      if(fgUsedVars[kVZEROXaXc+ih]) 
	values[kVZEROXaXc+ih] = eventF->Qx(EVENTPLANE::kVZEROA, ih+1)*eventF->Qx(EVENTPLANE::kVZEROC, ih+1);
      if(fgUsedVars[kVZEROXaYa+ih]) 
	values[kVZEROXaYa+ih] = eventF->Qx(EVENTPLANE::kVZEROA, ih+1)*eventF->Qy(EVENTPLANE::kVZEROA, ih+1);
      if(fgUsedVars[kVZEROXaYc+ih]) 
	values[kVZEROXaYc+ih] = eventF->Qx(EVENTPLANE::kVZEROA, ih+1)*eventF->Qy(EVENTPLANE::kVZEROC, ih+1);
      if(fgUsedVars[kVZEROYaXc+ih]) 
	values[kVZEROYaXc+ih] = eventF->Qy(EVENTPLANE::kVZEROA, ih+1)*eventF->Qx(EVENTPLANE::kVZEROC, ih+1);
      if(fgUsedVars[kVZEROYaYc+ih]) 
	values[kVZEROYaYc+ih] = eventF->Qy(EVENTPLANE::kVZEROA, ih+1)*eventF->Qy(EVENTPLANE::kVZEROC, ih+1);
      if(fgUsedVars[kVZEROXcYc+ih]) 
	values[kVZEROXcYc+ih] = eventF->Qx(EVENTPLANE::kVZEROC, ih+1)*eventF->Qy(EVENTPLANE::kVZEROC, ih+1);
      // Psi_A - Psi_C
      if(fgUsedVars[kVZEROdeltaRPac+ih])
        values[kVZEROdeltaRPac+ih] = DeltaPhi(eventF->EventPlane(EVENTPLANE::kVZEROA, ih+1), 
	  				      eventF->EventPlane(EVENTPLANE::kVZEROC, ih+1));
      
      // TPC event plane
      values[kTPCQvecX+ih] = eventF->Qx(EVENTPLANE::kTPC, ih+1);
      values[kTPCQvecY+ih] = eventF->Qy(EVENTPLANE::kTPC, ih+1);
      if(fgUsedVars[kTPCRP+ih]) 
	values[kTPCRP+ih] = eventF->EventPlane(EVENTPLANE::kTPC, ih+1);
      // TPC VZERO Q-vector correlations
      if(fgUsedVars[kRPXtpcXvzeroa+ih]) 
	values[kRPXtpcXvzeroa+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecX+ih];
      if(fgUsedVars[kRPXtpcXvzeroc+ih]) 
	values[kRPXtpcXvzeroc+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecX+6+ih];
      if(fgUsedVars[kRPYtpcYvzeroa+ih]) 
	values[kRPYtpcYvzeroa+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecY+ih];
      if(fgUsedVars[kRPYtpcYvzeroc+ih]) 
	values[kRPYtpcYvzeroc+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecY+6+ih];
      if(fgUsedVars[kRPXtpcYvzeroa+ih]) 
	values[kRPXtpcYvzeroa+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecY+ih];
      if(fgUsedVars[kRPXtpcYvzeroc+ih]) 
	values[kRPXtpcYvzeroc+ih] = values[kTPCQvecX+ih]*values[kVZEROQvecY+6+ih];
      if(fgUsedVars[kRPYtpcXvzeroa+ih]) 
	values[kRPYtpcXvzeroa+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecX+ih];
      if(fgUsedVars[kRPYtpcXvzeroc+ih]) 
	values[kRPYtpcXvzeroc+ih] = values[kTPCQvecY+ih]*values[kVZEROQvecX+6+ih];
      // Psi_TPC - Psi_VZERO A/C      
      if(fgUsedVars[kRPdeltaVZEROAtpc+ih]) 
	values[kRPdeltaVZEROAtpc+ih] = DeltaPhi(values[kVZERORP+0*6+ih], values[kTPCRP+ih]);
      if(fgUsedVars[kRPdeltaVZEROCtpc+ih])
        values[kRPdeltaVZEROCtpc+ih] = DeltaPhi(values[kVZERORP+1*6+ih], values[kTPCRP+ih]);
      // TPC event planes with sub-event method
      values[kTPCQvecXleft+ih] = eventF->Qx(EVENTPLANE::kTPCneg, ih+1);
      values[kTPCQvecYleft+ih] = eventF->Qy(EVENTPLANE::kTPCneg, ih+1);
      if(fgUsedVars[kTPCRPleft+ih])
	values[kTPCRPleft+ih] = eventF->EventPlane(EVENTPLANE::kTPCneg, ih+1);
      values[kTPCQvecXright+ih] = eventF->Qx(EVENTPLANE::kTPCpos, ih+1);
      values[kTPCQvecYright+ih] = eventF->Qy(EVENTPLANE::kTPCpos, ih+1);
      if(fgUsedVars[kTPCRPright+ih])
        values[kTPCRPright+ih] = eventF->EventPlane(EVENTPLANE::kTPCpos, ih+1); 
      if(fgUsedVars[kTPCsubResCos+ih]) 
	values[kTPCsubResCos+ih] = TMath::Cos(Double_t(ih+1)*(values[kTPCRPleft+ih]-values[kTPCRPright+ih]));
    }  // end loop over harmonics
    
    // VZERO v2 using TPC event plane
    Double_t vzeroChannelPhi[8] = {0.3927, 1.1781, 1.9635, 2.7489, -2.7489, -1.9635, -1.1781, -0.3927};
    
    for(Int_t ich=0; ich<64; ++ich) {
      if(fgUsedVars[kVZEROflowV2TPC+ich])
	values[kVZEROflowV2TPC+ich] = values[kVZEROChannelMult+ich]*
                                      TMath::Cos(2.0*DeltaPhi(vzeroChannelPhi[ich%8],values[kTPCRP+1]));
    } 
  }  // end if (eventF)
    
  for(Int_t izdc=0;izdc<10;++izdc) values[kZDCnEnergyCh+izdc] = event->EnergyZDCnTree(izdc);
  for(Int_t izdc=0;izdc<10;++izdc) values[kZDCpEnergyCh+izdc] = event->EnergyZDCpTree(izdc);
  for(Int_t itzero=0;itzero<26;++itzero) values[kTZEROAmplitudeCh+itzero] = event->AmplitudeTZEROch(itzero);
  for(Int_t itzero=0;itzero<3;++itzero) values[kTZEROTOF+itzero] = event->EventTZEROStartTimeTOFfirst(itzero);
  for(Int_t itzero=0;itzero<3;++itzero) values[kTZEROTOFbest+itzero] = event->EventTZEROStartTimeTOFbest(itzero);
  values[kTZEROzVtx] = event->VertexTZERO();
  values[kTZEROstartTime] = event->EventTZEROStartTime();
  values[kTZEROpileup] = event->IsPileupTZERO();
  values[kTZEROsatellite] = event->IsSatteliteCollisionTZERO();  

  values[kMultEstimatorV0M]          = event->MultEstimatorV0M();
  values[kMultEstimatorV0A]          = event->MultEstimatorV0A();
  values[kMultEstimatorV0C]          = event->MultEstimatorV0C();
  values[kMultEstimatorOnlineV0M]    = event->MultEstimatorOnlineV0M();
  values[kMultEstimatorOnlineV0A]    = event->MultEstimatorOnlineV0A();
  values[kMultEstimatorOnlineV0C]    = event->MultEstimatorOnlineV0C();
  values[kMultEstimatorADM]          = event->MultEstimatorADM();
  values[kMultEstimatorADA]          = event->MultEstimatorADA();
  values[kMultEstimatorADC]          = event->MultEstimatorADC();
  values[kMultEstimatorSPDClusters]  = event->MultEstimatorSPDClusters();
  values[kMultEstimatorSPDTracklets] = event->MultEstimatorSPDTracklets();
  values[kMultEstimatorRefMult05]    = event->MultEstimatorRefMult05();
  values[kMultEstimatorRefMult08]    = event->MultEstimatorRefMult08();

  values[kMultEstimatorPercentileV0M]    = event->MultEstimatorPercentileV0M();
  values[kMultEstimatorPercentileV0A]    = event->MultEstimatorPercentileV0A();
  values[kMultEstimatorPercentileV0C]    = event->MultEstimatorPercentileV0C();
  values[kMultEstimatorPercentileOnlineV0M]    = event->MultEstimatorPercentileOnlineV0M();
  values[kMultEstimatorPercentileOnlineV0A]    = event->MultEstimatorPercentileOnlineV0A();
  values[kMultEstimatorPercentileOnlineV0C]    = event->MultEstimatorPercentileOnlineV0C();
  values[kMultEstimatorPercentileADM]          = event->MultEstimatorPercentileADM();
  values[kMultEstimatorPercentileADA]          = event->MultEstimatorPercentileADA();
  values[kMultEstimatorPercentileADC]          = event->MultEstimatorPercentileADC();
  values[kMultEstimatorPercentileSPDClusters]  = event->MultEstimatorPercentileSPDClusters();
  values[kMultEstimatorPercentileSPDTracklets] = event->MultEstimatorPercentileSPDTracklets();
  values[kMultEstimatorPercentileRefMult05]    = event->MultEstimatorPercentileRefMult05();
  values[kMultEstimatorPercentileRefMult08]    = event->MultEstimatorPercentileRefMult08();

}

//_________________________________________________________________
void AliReducedVarManager::FillITSlayerFlag(TRACK* track, Int_t layer, Float_t* values) {
  //
  // fill the ITS layer hit
  //
  values[kITSlayerHit] = -1.0*(layer+1);
  if(fgUsedVars[kITSlayerHit] && track->ITSLayerHit(layer)) values[kITSlayerHit] = layer+1;
}

//_________________________________________________________________
void AliReducedVarManager::FillITSsharedLayerFlag(TRACK* track, Int_t layer, Float_t* values) {
   //
   // fill the ITS layer having shared cluster
   //
   values[kITSlayerShared] = -1.0*(layer+1);
   if(fgUsedVars[kITSlayerShared] && track->ITSLayerHit(layer) && track->ITSClusterIsShared(layer)) values[kITSlayerShared] = layer+1;
}

//_________________________________________________________________
void AliReducedVarManager::FillL0TriggerInputs(EVENT* event, Int_t input, Float_t* values, Int_t input2 /*=999*/) {
  //
  // fill the L0 trigger inputs
  //
  values[kL0TriggerInput] = -1.0;
  if(fgUsedVars[kL0TriggerInput] && event->L0TriggerInput(input)) values[kL0TriggerInput] = input;
  values[kL0TriggerInput2] = -1.0;
  if(fgUsedVars[kL0TriggerInput2] && event->L0TriggerInput(input2)) values[kL0TriggerInput2] = input2;
}


//_________________________________________________________________
void AliReducedVarManager::FillL1TriggerInputs(EVENT* event, Int_t input, Float_t* values, Int_t input2 /*=999*/) {
  //
  // fill the L1 trigger inputs
  //
  values[kL1TriggerInput] = -1.0;
  if(fgUsedVars[kL1TriggerInput] && event->L1TriggerInput(input)) values[kL1TriggerInput] = input;
  values[kL1TriggerInput2] = -1.0;
  if(fgUsedVars[kL1TriggerInput2] && event->L1TriggerInput(input2)) values[kL1TriggerInput2] = input2;
}

//_________________________________________________________________
void AliReducedVarManager::FillL2TriggerInputs(EVENT* event, Int_t input, Float_t* values, Int_t input2 /*=999*/) {
  //
  // fill the L2 trigger inputs
  //
  values[kL2TriggerInput] = -1.0;
  if(fgUsedVars[kL2TriggerInput] && event->L2TriggerInput(input)) values[kL2TriggerInput] = input;
  values[kL2TriggerInput2] = -1.0;
  if(fgUsedVars[kL2TriggerInput2] && event->L2TriggerInput(input2)) values[kL2TriggerInput2] = input2;
}

//_________________________________________________________________
void AliReducedVarManager::FillEventTagInput(BASEEVENT* event, Int_t input, Float_t* values) {
  //
  // fill the event tag inputs
  //
  values[kEventTag] = -1.0;
  if(fgUsedVars[kEventTag] && event->EventTag(input)) values[kEventTag] = input;
}

//_________________________________________________________________
void AliReducedVarManager::FillTPCclusterBitFlag(TRACK* track, Int_t bit, Float_t* values) {
  //
  // fill the TPC cluster map
  //
  values[kTPCclusBitFired] = -1;
  if(fgUsedVars[kTPCclusBitFired] && track->TPCClusterMapBitFired(bit)) values[kTPCclusBitFired] = bit;
}


//_________________________________________________________________
void AliReducedVarManager::FillTrackingStatus(TRACK* p, Float_t* values) {
  //
  // Fill tracking flags
  //
  for(Int_t i=0; i<kNTrackingStatus; i++)
    values[kTrackingStatus+i] = p->CheckTrackStatus(i)+10e-6;
}


//_________________________________________________________________
void AliReducedVarManager::FillTrackingFlag(TRACK* track, UInt_t flag, Float_t* values) {
  //
  // fill the tracking flag
  //
  values[kTrackingFlag] = -1;
  if(track->CheckTrackStatus(flag)) values[kTrackingFlag] = flag;
}

//_________________________________________________________________
void AliReducedVarManager::FillTrackQualityFlag(BASETRACK* track, UShort_t flag, Float_t* values, UShort_t flag2 /*=999*/) {
  //
  // fill the track quality flag
  //
  values[kTrackQualityFlag] = -1;
  if(track->TestQualityFlag(flag)) values[kTrackQualityFlag] = flag;
  values[kTrackQualityFlag2] = -1;
  if(track->TestQualityFlag(flag2)) values[kTrackQualityFlag2] = flag2;
}

//_________________________________________________________________
void AliReducedVarManager::FillTrackMCFlag(BASETRACK* track, UShort_t flag, Float_t* values, UShort_t flag2 /*=999*/) {
   //
   // fill the track MC flag
   //
   values[kTrackMCFlag] = -1;
   if(flag<32 && track->TestMCFlag(flag)) values[kTrackMCFlag] = flag;
   values[kTrackMCFlag2] = -1;
   if(flag2<32 && track->TestMCFlag(flag2)) values[kTrackMCFlag2] = flag2;
}

//_________________________________________________________________
void AliReducedVarManager::FillPairQualityFlag(PAIR* p, UShort_t flag, Float_t* values, UShort_t flag2 /*=999*/) {
  //
  // fill the pair quality flag
  //
  values[kPairQualityFlag] = -1;
  if(p->TestQualityFlag(flag)) values[kPairQualityFlag] = flag;
  values[kPairQualityFlag2] = -1;
  if(p->TestQualityFlag(flag2)) values[kPairQualityFlag2] = flag2;
}

//_________________________________________________________________
void AliReducedVarManager::FillEventOnlineTriggers(AliReducedEventInfo* event, Float_t* values){
  //
  // fill the trigger bit input
  //
  for(UShort_t i=0; i<kNTriggers; i++){
    values[kOnlineTriggersFired+i] = (event->TriggerMask()&(ULong_t(1)<<i) ? 1.0 : 0.0);
  }
}

//_________________________________________________________________
void AliReducedVarManager::FillEventOnlineTrigger(UShort_t triggerBit, Float_t* values, UShort_t triggerBit2 /*=999*/) {
  //
  // fill the trigger bit input
  //  The second trigger bit (triggerBit2) is used for correlation histograms between the different trigger inputs
  //
  if(triggerBit>=64) return;
  if(!fgEvent) return;
  values[kOnlineTrigger] = triggerBit;
  values[kOnlineTriggerFired] = (((AliReducedEventInfo*)fgEvent)->TriggerMask()&(ULong_t(1)<<triggerBit) ? triggerBit : -1.0);
  values[kOnlineTriggerFired2] = 0.0;
  if(triggerBit<64)
     values[kOnlineTriggerFired2] = (((AliReducedEventInfo*)fgEvent)->TriggerMask()&(ULong_t(1)<<triggerBit2) ? triggerBit2 : -1.0);
}

//________________________________________________________________
void AliReducedVarManager::FillV0Channel(Int_t ich, Float_t* values) {
   //
   // Update the V0 channel number
   //
   if(ich>=0 && ich<64) {
      values[kVZEROCurrentChannel] = ich;
      values[kVZEROCurrentChannelMult] = values[kVZEROChannelMult+ich];
      values[kVZEROCurrentChannelMultCalib] = values[kVZEROChannelMultCalib+ich];
   }
   else {
      values[kVZEROCurrentChannel] = -1;
      values[kVZEROCurrentChannelMult] = -999.;
      values[kVZEROCurrentChannelMultCalib] = -999.;
   }
}

//_________________________________________________________________
void AliReducedVarManager::FillMCTruthInfo(TRACK* p, Float_t* values, TRACK* leg1 /* = 0x0 */, TRACK* leg2 /* = 0x0 */) {
   //
   //  Fill pure MC truth information
   //
   if(fgUsedVars[kPtMC]) values[kPtMC] = p->PtMC();
   if(fgUsedVars[kPMC]) values[kPMC] = p->PMC();
   values[kPxMC] = p->MCmom(0);
   values[kPyMC] = p->MCmom(1);
   values[kPzMC] = p->MCmom(2);
   if(fgUsedVars[kPt_weight]) values[kPt_weight] =  calculate_weight_factor(values[kPtMC],values[kCentVZERO]);
   
   if(fgUsedVars[kThetaMC]) values[kThetaMC] = p->ThetaMC();
   if(fgUsedVars[kEtaMC]) values[kEtaMC] = p->EtaMC();
   if(fgUsedVars[kPhiMC]) values[kPhiMC] = p->PhiMC();
   if(fgUsedVars[kMassMC]) {
      if(TMath::Abs(p->MCPdg(0))==443)
         values[kMassMC] = fgkPairMass[AliReducedPairInfo::kJpsiToEE];  
      if(TMath::Abs(p->MCPdg(0))==100443)
         values[kMassMC] = fgkPairMass[AliReducedPairInfo::kPsi2SToEE];  
   }
   if(fgUsedVars[kRapMC]) {
      if(TMath::Abs(p->MCPdg(0))==443)
         values[kRapMC] = p->RapidityMC(fgkPairMass[AliReducedPairInfo::kJpsiToEE]); 
      if(TMath::Abs(p->MCPdg(0))==100443)
         values[kRapMC] = p->RapidityMC(fgkPairMass[AliReducedPairInfo::kPsi2SToEE]); 
   }
  if(fgUsedVars[kRapMCAbs]) {
    if(TMath::Abs(p->MCPdg(0))==443)
      values[kRapMCAbs] = TMath::Abs(p->RapidityMC(fgkPairMass[AliReducedPairInfo::kJpsiToEE]));
    if(TMath::Abs(p->MCPdg(0))==100443)
       values[kRapMCAbs] = TMath::Abs(p->RapidityMC(fgkPairMass[AliReducedPairInfo::kPsi2SToEE]));
  }

  if(fgUsedVars[kPseudoProperDecayTimeMC]){
     if(fgEvent->IsA()==EVENT::Class()){
     EVENT* eventInfo = (EVENT*)fgEvent;
     Double_t lxyMC = ( (p->MCFreezeout(0) - eventInfo->VertexMC(0)) * p->MCmom(0) + (p->MCFreezeout(1) - eventInfo->VertexMC(1)) * p->MCmom(1) ) / p->PtMC();
     values[kPseudoProperDecayTimeMC] = lxyMC * (fgkPairMass[AliReducedPairInfo::kJpsiToEE])/p->PtMC();
     }
  }
   // compute MC truth variables from decay legs, e.g. from the 2 electrons of a J/psi decay
   // NOTE: this may be different from the kinematics of the mother, if not all decay legs are considered / tracked
   Bool_t requestMCfromLegs = kFALSE;
   if(fgUsedVars[kPtMCfromLegs] || fgUsedVars[kPMCfromLegs] || 
      fgUsedVars[kPxMCfromLegs] || fgUsedVars[kPyMCfromLegs] || fgUsedVars[kPzMCfromLegs] ||
      fgUsedVars[kThetaMCfromLegs] || fgUsedVars[kEtaMCfromLegs] || fgUsedVars[kPhiMCfromLegs] ||
      fgUsedVars[kMassMCfromLegs] || fgUsedVars[kRapMCfromLegs] ||
      fgUsedVars[kPairLegPtMC] || fgUsedVars[kPairLegPtMC+1] || fgUsedVars[kPairLegPtMCSum]) 
      requestMCfromLegs = kTRUE;
   
   if(leg1 && leg2 && requestMCfromLegs) {
      values[kPxMCfromLegs] = leg1->MCmom(0) + leg2->MCmom(0);
      values[kPyMCfromLegs] = leg1->MCmom(1) + leg2->MCmom(1);
      values[kPzMCfromLegs] = leg1->MCmom(2) + leg2->MCmom(2);
      values[kPairLegPtMC+0] = TMath::Sqrt(leg1->MCmom(0)*leg1->MCmom(0)+leg1->MCmom(1)*leg1->MCmom(1));
      values[kPairLegPtMC+1] = TMath::Sqrt(leg2->MCmom(0)*leg2->MCmom(0)+leg2->MCmom(1)*leg2->MCmom(1));
      values[kPairLegPtMCSum] = values[kPairLegPtMC]+values[kPairLegPtMC+1];
      values[kPtMCfromLegs] = TMath::Sqrt(values[kPxMCfromLegs]*values[kPxMCfromLegs]+values[kPyMCfromLegs]*values[kPyMCfromLegs]);
      values[kPMCfromLegs] = TMath::Sqrt(values[kPtMCfromLegs]*values[kPtMCfromLegs]+values[kPzMCfromLegs]*values[kPzMCfromLegs]);
      values[kThetaMCfromLegs] = (values[kPMCfromLegs]>=1.0e-6 ? TMath::ACos(values[kPzMCfromLegs]/values[kPMCfromLegs]) : 0.0);
      values[kEtaMCfromLegs] = TMath::Tan(0.5*values[kThetaMCfromLegs]);
      values[kEtaMCfromLegs] = (values[kEtaMCfromLegs]>1.0e-6 ? -1.0*TMath::Log(values[kEtaMCfromLegs]) : 0.0);
      values[kPhiMCfromLegs] = TMath::ATan2(values[kPyMCfromLegs],values[kPxMCfromLegs]);
      values[kPhiMCfromLegs] = (values[kPhiMCfromLegs]<0.0 ? (TMath::TwoPi()+values[kPhiMCfromLegs]) : values[kPhiMCfromLegs]);
      if(TMath::Abs(p->MCPdg(0))==443 || TMath::Abs(p->MCPdg(0))==100443) {
         Float_t m1 = fgkParticleMass[kElectron];
         Float_t m2 = fgkParticleMass[kElectron];
         values[kMassMCfromLegs] = m1*m1+m2*m2 + 
         2.0*(TMath::Sqrt(m1*m1+leg1->P()*leg1->P())*TMath::Sqrt(m2*m2+leg2->P()*leg2->P()) - 
         leg1->Px()*leg2->Px() - leg1->Py()*leg2->Py() - leg1->Pz()*leg2->Pz());
         if(values[kMassMCfromLegs]<0.0) {
            cout << "FillMCTruthInfo(mother, values, track1, track2): Warning: Very small squared mass found. "
            << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
            cout << "   mass^2: " << values[kMassMCfromLegs] << endl;
            cout << "p1(p,x,y,z): " << leg1->P() << ", " << leg1->Px() << ", " << leg1->Py() << ", " << leg1->Pz() << endl;
            cout << "p2(p,x,y,z): " << leg2->P() << ", " << leg2->Px() << ", " << leg2->Py() << ", " << leg2->Pz() << endl;
            values[kMassMCfromLegs] = 0.0;
         }
         else
            values[kMassMCfromLegs] = TMath::Sqrt(values[kMassMCfromLegs]);
         
         Float_t e = TMath::Sqrt(values[kPMCfromLegs]*values[kPMCfromLegs] + values[kMassMCfromLegs] * values[kMassMCfromLegs]);
         Float_t factor = e - values[kPzMCfromLegs];
         values[kRapMCfromLegs] = (TMath::Abs(factor)>1.0e-6 ? (e+values[kPzMCfromLegs])/factor : 0.0);
         values[kRapMCfromLegs] = (values[kRapMCfromLegs]>1.0e-6 ? 0.5*TMath::Log(values[kRapMCfromLegs]) : 0.0);
      }
   }
   
   values[kPdgMC] = p->MCPdg(0);
   values[kPdgMC+1] = p->MCPdg(1);
   values[kPdgMC+2] = p->MCPdg(2);
   values[kPdgMC+3] = p->MCPdg(3);
   
   // polarization variables
   Bool_t usePolarization=kFALSE;
   if(leg1 && leg2 && (fgUsedVars[kPairThetaCS] || fgUsedVars[kPairThetaHE] || fgUsedVars[kPairPhiCS] || fgUsedVars[kPairPhiHE]))
      usePolarization = kTRUE;
   if(usePolarization)
      GetThetaPhiCM(leg1, leg2, values[kPairThetaHE], values[kPairPhiHE], values[kPairThetaCS], values[kPairPhiCS]);
}


//_________________________________________________________________
void AliReducedVarManager::FillMCTruthInfo(TRACK* leg1, TRACK* leg2, Float_t* values) {
   //
   // Compute MC truth variables from decay legs only, no mother assumption used
   // NOTE: This function is done specifically for pairs of MC particles which do not have any mother in the stack, e.g. electrons from gamma-gamma processes produced with Starlight
   //      All quantities are computed just from the leg kinematics
   //
   
   values[kPxMCfromLegs] = leg1->MCmom(0) + leg2->MCmom(0);
   values[kPyMCfromLegs] = leg1->MCmom(1) + leg2->MCmom(1);
   values[kPzMCfromLegs] = leg1->MCmom(2) + leg2->MCmom(2);
   values[kPtMCfromLegs] = TMath::Sqrt(values[kPxMCfromLegs]*values[kPxMCfromLegs]+values[kPyMCfromLegs]*values[kPyMCfromLegs]);
   values[kPMCfromLegs] = TMath::Sqrt(values[kPtMCfromLegs]*values[kPtMCfromLegs]+values[kPzMCfromLegs]*values[kPzMCfromLegs]);
   values[kThetaMCfromLegs] = (values[kPMCfromLegs]>=1.0e-6 ? TMath::ACos(values[kPzMCfromLegs]/values[kPMCfromLegs]) : 0.0);
   values[kEtaMCfromLegs] = TMath::Tan(0.5*values[kThetaMCfromLegs]);
   values[kEtaMCfromLegs] = (values[kEtaMCfromLegs]>1.0e-6 ? -1.0*TMath::Log(values[kEtaMCfromLegs]) : 0.0);
   values[kPhiMCfromLegs] = TMath::ATan2(values[kPyMCfromLegs],values[kPxMCfromLegs]);
   values[kPhiMCfromLegs] = (values[kPhiMCfromLegs]<0.0 ? (TMath::TwoPi()+values[kPhiMCfromLegs]) : values[kPhiMCfromLegs]);
   values[kPairLegPtMC+0] = TMath::Sqrt(leg1->MCmom(0)*leg1->MCmom(0)+leg1->MCmom(1)*leg1->MCmom(1));
   values[kPairLegPtMC+1] = TMath::Sqrt(leg2->MCmom(0)*leg2->MCmom(0)+leg2->MCmom(1)*leg2->MCmom(1));
   values[kPairLegPtMCSum] = values[kPairLegPtMC]+values[kPairLegPtMC+1];
   
   // NOTE: for other particle species, check the PDG code of the particles
   Float_t m1 = fgkParticleMass[kElectron];
   Float_t m2 = fgkParticleMass[kElectron];
   values[kMassMCfromLegs] = m1*m1+m2*m2 + 
      2.0*(TMath::Sqrt(m1*m1+leg1->P()*leg1->P())*TMath::Sqrt(m2*m2+leg2->P()*leg2->P()) - 
            leg1->Px()*leg2->Px() - leg1->Py()*leg2->Py() - leg1->Pz()*leg2->Pz());
   if(values[kMassMCfromLegs]<0.0) {
      cout << "FillMCTruthInfo(track1, track2, values): Warning: Very small squared mass found. "
      << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
      cout << "   mass^2: " << values[kMassMCfromLegs] << endl;
      cout << "p1(p,x,y,z): " << leg1->P() << ", " << leg1->Px() << ", " << leg1->Py() << ", " << leg1->Pz() << endl;
      cout << "p2(p,x,y,z): " << leg2->P() << ", " << leg2->Px() << ", " << leg2->Py() << ", " << leg2->Pz() << endl;
      values[kMassMCfromLegs] = 0.0;
   }
   else
      values[kMassMCfromLegs] = TMath::Sqrt(values[kMassMCfromLegs]);
   values[kMassMC] = values[kMassMCfromLegs];
   
   Float_t e = TMath::Sqrt(values[kPMCfromLegs]*values[kPMCfromLegs] + values[kMassMCfromLegs] * values[kMassMCfromLegs]);
   Float_t factor = e - values[kPzMCfromLegs];
   values[kRapMCfromLegs] = (TMath::Abs(factor)>1.0e-6 ? (e+values[kPzMCfromLegs])/factor : 0.0);
   values[kRapMCfromLegs] = (values[kRapMCfromLegs]>1.0e-6 ? 0.5*TMath::Log(values[kRapMCfromLegs]) : 0.0);
   
   values[kPtMC] = values[kPtMCfromLegs];
   values[kPMC] = values[kPMCfromLegs];
   values[kPxMC] = values[kPxMCfromLegs];
   values[kPyMC] = values[kPyMCfromLegs];
   values[kPzMC] = values[kPzMCfromLegs];
   values[kThetaMC] = values[kThetaMCfromLegs];
   values[kEtaMC] = values[kEtaMCfromLegs];
   values[kPhiMC] = values[kPhiMCfromLegs];
   values[kMassMC] = values[kMassMCfromLegs];
   values[kRapMC] = values[kRapMCfromLegs];
   values[kRapMCAbs] = TMath::Abs(values[kRapMCfromLegs]); 
      
   // polarization variables
   Bool_t usePolarization=kFALSE;
   if(leg1 && leg2 && (fgUsedVars[kPairThetaCS] || fgUsedVars[kPairThetaHE] || fgUsedVars[kPairPhiCS] || fgUsedVars[kPairPhiHE]))
      usePolarization = kTRUE;
   if(usePolarization)
      GetThetaPhiCM(leg1, leg2, values[kPairThetaHE], values[kPairPhiHE], values[kPairThetaCS], values[kPairPhiCS]);
}


//_________________________________________________________________
void AliReducedVarManager::FillTrackInfo(BASETRACK* p, Float_t* values) {
  //
  // fill track information
  //
  
  // Fill base track information
  if(fgUsedVars[kPt])        values[kPt]        = p->Pt();
  if(fgUsedVars[kPtSquared]) values[kPtSquared] = values[kPt]*values[kPt];
  if(fgUsedVars[kOneOverSqrtPt]) {
    values[kOneOverSqrtPt] = values[kPt] > 0. ? 1./TMath::Sqrt(values[kPt]) : 999.;
  }
  if(fgUsedVars[kP])         values[kP]         = p->P();
  if(fgUsedVars[kPx])        values[kPx]        = p->Px();
  if(fgUsedVars[kPy])        values[kPy]        = p->Py();
  if(fgUsedVars[kPz])        values[kPz]        = p->Pz();
  if(fgUsedVars[kTheta])     values[kTheta]     = p->Theta();
  if(fgUsedVars[kPhi])       values[kPhi]       = p->Phi();
  if(fgUsedVars[kEta])       values[kEta]       = p->Eta();
  for(Int_t ih=1; ih<=6; ++ih) {
     if(fgUsedVars[kCosNPhi+ih-1]) values[kCosNPhi+ih-1] = TMath::Cos(p->Phi()*ih);
     if(fgUsedVars[kSinNPhi+ih-1]) values[kSinNPhi+ih-1] = TMath::Sin(p->Phi()*ih);
  }
  values[kCharge] = p->Charge();
  
  //pair efficiency variables
  if((fgUsedVars[kPairEff] || fgUsedVars[kOneOverPairEff] || fgUsedVars[kOneOverPairEffSq]) && fgPairEffMap) {
    Int_t binX = 0;
    if (fgEffMapVarDependencyX!=kNothing) {
      binX = fgPairEffMap->GetXaxis()->FindBin(values[fgEffMapVarDependencyX]);
      if(binX==0) binX = 1;
      if(binX==fgPairEffMap->GetXaxis()->GetNbins()+1) binX -= 1;
    }
    Int_t binY = 0;
    if (fgEffMapVarDependencyY!=kNothing) {
      binY = fgPairEffMap->GetYaxis()->FindBin(values[fgEffMapVarDependencyY]);
      if(binY==0) binY = 1;
      if(binY==fgPairEffMap->GetXaxis()->GetNbins()+1) binY -= 1;
    }
    Int_t binZ = 0;
    if (fgEffMapVarDependencyZ!=kNothing) {
      binZ = fgPairEffMap->GetZaxis()->FindBin(values[fgEffMapVarDependencyZ]);
      if(binZ==0) binZ = 1;
      if(binZ==fgPairEffMap->GetZaxis()->GetNbins()+1) binZ -= 1;
    }

    Float_t                   pairEff = 1.;
    if (binX && binY && binZ) pairEff = fgPairEffMap->GetBinContent(binX, binY, binZ);
    else if (binX && binY)    pairEff = fgPairEffMap->GetBinContent(binX, binY);
    else if (binX)            pairEff = fgPairEffMap->GetBinContent(binX);

    Float_t               oneOverPairEff = 1.;
    if (pairEff > 1.0e-6) oneOverPairEff = 1/pairEff;

    values[kPairEff]          = pairEff;
    values[kOneOverPairEff]   = oneOverPairEff;
    values[kOneOverPairEffSq] = oneOverPairEff*oneOverPairEff;
  }

  // Fill VZERO flow variables
  for(Int_t iVZEROside=0; iVZEROside<3; ++iVZEROside) {
     for(Int_t ih=0; ih<6; ++ih) {
        if(fgUsedVars[kVZEROFlowVn+iVZEROside*6+ih])
           values[kVZEROFlowVn+iVZEROside*6+ih] = TMath::Cos((values[kPhi]-values[kVZERORP+iVZEROside*6+ih])*(ih+1));
        if(fgUsedVars[kVZEROFlowSine+iVZEROside*6+ih])
           values[kVZEROFlowSine+iVZEROside*6+ih] = TMath::Sin((values[kPhi]-values[kVZERORP+iVZEROside*6+ih])*(ih+1));
        if(iVZEROside<2) {
           if(fgUsedVars[kVZEROuQ+iVZEROside*6+ih]) {
              values[kVZEROuQ+iVZEROside*6+ih] = TMath::Cos((values[kPhi]-values[kVZERORP+iVZEROside*6+ih])*(ih+1));
              values[kVZEROuQ+iVZEROside*6+ih] *= TMath::Sqrt(values[kVZEROQvecX+iVZEROside*6+ih]*values[kVZEROQvecX+iVZEROside*6+ih] +
              values[kVZEROQvecY+iVZEROside*6+ih]*values[kVZEROQvecY+iVZEROside*6+ih]); 
           }
           if(fgUsedVars[kVZEROuQsine+iVZEROside*6+ih]) {
              values[kVZEROuQsine+iVZEROside*6+ih] = TMath::Sin((values[kPhi]-values[kVZERORP+iVZEROside*6+ih])*(ih+1));
              values[kVZEROuQsine+iVZEROside*6+ih] *= TMath::Sqrt(values[kVZEROQvecX+iVZEROside*6+ih]*values[kVZEROQvecX+iVZEROside*6+ih] +
              values[kVZEROQvecY+iVZEROside*6+ih]*values[kVZEROQvecY+iVZEROside*6+ih]); 
           }	    
        }
     }  // end loop over harmonics
  }  // end loop over VZERO sides
  
  // Fill TPC flow variables
  // Subtract the q vector of the track or of the pair legs from the event q-vector 
  Bool_t tpcEPUsed = kFALSE;
  for(Int_t ih=0; ih<6; ++ih) {
     if(fgUsedVars[kTPCFlowVn+ih]) {tpcEPUsed = kTRUE; break;}
     if(fgUsedVars[kTPCFlowSine+ih]) {tpcEPUsed = kTRUE; break;}
     if(fgUsedVars[kTPCuQ+ih]) {tpcEPUsed = kTRUE; break;}
     if(fgUsedVars[kTPCuQsine+ih]) {tpcEPUsed = kTRUE; break;}
  }

  if(tpcEPUsed) {
//      Float_t tpcEPsubtracted[6] = {0.0};
//      Double_t qVec[6][2] = {{0.0}};
//      for(Int_t ih=0; ih<6; ++ih) {qVec[ih][0]=values[kTPCQvecXtotal+ih]; qVec[ih][1]=values[kTPCQvecYtotal+ih];}
//      EVENT* eventInfo = NULL;
//      if(fgEvent->IsA()==EVENT::Class()) eventInfo = (EVENT*)fgEvent;
//      if((p->IsA() == AliReducedTrackInfo::Class()) && eventInfo) {
//         eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)p,qVec,EVENTPLANE::kTPC,-0.8,-0.5*fgkTPCQvecRapGap);
//         eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)p,qVec,EVENTPLANE::kTPC,0.5*fgkTPCQvecRapGap,0.8);
//      }

     // TODO: Make sure the pair legs are properly subtracted from the TPC event plane calculation
     //              For the moment this part of the code is commented out
     /* else if((p->IsA() == AliReducedPairInfo::Class()) && eventInfo) {
        cout<<"id  "<<((AliReducedPairInfo*)p)->LegId(1)<<endl;
        eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)eventInfo->GetTrack(((AliReducedPairInfo*)p)->LegId(0)),qVec,EVENTPLANE::kTPC,-0.8,-0.5*fgkTPCQvecRapGap);
        eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)eventInfo->GetTrack(((AliReducedPairInfo*)p)->LegId(0)),qVec,EVENTPLANE::kTPC,0.5*fgkTPCQvecRapGap,0.8);
        eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)eventInfo->GetTrack(((AliReducedPairInfo*)p)->LegId(1)),qVec,EVENTPLANE::kTPC,-0.8,-0.5*fgkTPCQvecRapGap);
        eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)eventInfo->GetTrack(((AliReducedPairInfo*)p)->LegId(1)),qVec,EVENTPLANE::kTPC,0.5*fgkTPCQvecRapGap,0.8);
     }  */
     // recalculate the TPC event plane
//      for(Int_t ih=0; ih<6;++ih) 
//         tpcEPsubtracted[ih] = TMath::ATan2(qVec[ih][1], qVec[ih][0])/Double_t(ih+1);
//      for(Int_t ih=0; ih<6; ++ih) {
//         // vn using Psi_n
//         if(fgUsedVars[kTPCFlowVn+ih])
//            values[kTPCFlowVn+ih] = TMath::Cos(DeltaPhi(values[kPhi],tpcEPsubtracted[ih])*(ih+1));
//         if(fgUsedVars[kTPCFlowSine+ih]) 
//            values[kTPCFlowSine+ih] = TMath::Sin(DeltaPhi(values[kPhi],tpcEPsubtracted[ih])*(ih+1));
//         if(fgUsedVars[kTPCuQ+ih]) {
//            values[kTPCuQ+ih] = TMath::Cos((values[kPhi]-tpcEPsubtracted[ih])*(ih+1));
//            values[kTPCuQ+ih] *= TMath::Sqrt(qVec[ih][0]*qVec[ih][0] + qVec[ih][1]*qVec[ih][1]);
//         }
//         if(fgUsedVars[kTPCuQsine+ih]) {
//            values[kTPCuQsine+ih] = TMath::Sin((values[kPhi]-tpcEPsubtracted[ih])*(ih+1));
//            values[kTPCuQsine+ih] *= TMath::Sqrt(qVec[ih][0]*qVec[ih][0] + qVec[ih][1]*qVec[ih][1]);
//         }
//      }

        for(Int_t ih=0; ih<6; ++ih) {
        // vn using Psi_n
        if(fgUsedVars[kTPCFlowVn+ih])
           values[kTPCFlowVn+ih] = TMath::Cos(DeltaPhi(values[kPhi],values[kTPCRPtree+ih])*(ih+1));
        if(fgUsedVars[kTPCFlowSine+ih]) 
           values[kTPCFlowSine+ih] = TMath::Sin(DeltaPhi(values[kPhi],values[kTPCRPtree+ih])*(ih+1));
        if(fgUsedVars[kTPCuQ+ih]) {
           values[kTPCuQ+ih] = TMath::Cos((values[kPhi]-values[kTPCRPtree+ih])*(ih+1));
           values[kTPCuQ+ih] *= TMath::Sqrt(values[kVZEROQvecX+ih]*values[kVZEROQvecX+ih] + values[kVZEROQvecY+ih]*values[kVZEROQvecY+ih]);
        }
        if(fgUsedVars[kTPCuQsine+ih]) {
           values[kTPCuQsine+ih] = TMath::Sin((values[kPhi]-values[kTPCRPtree+ih])*(ih+1));
           values[kTPCuQsine+ih] *= TMath::Sqrt(values[kVZEROQvecX+ih]*values[kVZEROQvecX+ih] + values[kVZEROQvecY+ih]*values[kVZEROQvecY+ih]);
        }
     }
  }

  if(p->IsA()!=TRACK::Class()) return;
  TRACK* pinfo = (TRACK*)p;

  values[kPtTPC]       = pinfo->PtTPC();
  values[kTrackLength] = pinfo->TrackLength();
  values[kChi2TPCConstrainedVsGlobal] = pinfo->Chi2TPCConstrainedVsGlobal();
  values[kMassUsedForTracking] = pinfo->MassForTracking();
  values[kPhiTPC]      = pinfo->PhiTPC();
  values[kEtaTPC]      = pinfo->EtaTPC();
  values[kPin]         = pinfo->Pin();
  values[kDcaXY]       = pinfo->DCAxy();
  values[kDcaZ]        = pinfo->DCAz();
  values[kDcaXYTPC]    = pinfo->DCAxyTPC();
  values[kDcaZTPC]     = pinfo->DCAzTPC();

  if(fgUsedVars[kITSncls]) values[kITSncls] = pinfo->ITSncls();
  values[kITSsignal] = pinfo->ITSsignal();
  values[kITSchi2] = pinfo->ITSchi2();

  if(fgUsedVars[kITSnclsShared]) values[kITSnclsShared] = pinfo->ITSnSharedCls();
  values[kTPCncls] = pinfo->TPCncls();

  if(fgUsedVars[kNclsSFracITS])
  values[kNclsSFracITS] = (pinfo-> ITSncls()>0 ? Float_t (pinfo->ITSnSharedCls())/Float_t(pinfo->ITSncls()) :0.0) ;
  if(fgUsedVars[kTPCnclsRatio]) 
    values[kTPCnclsRatio] = (pinfo->TPCFindableNcls()>0 ? Float_t(pinfo->TPCncls())/Float_t(pinfo->TPCFindableNcls()) : 0.0);
  if(fgUsedVars[kTPCnclsRatio2]) 
    values[kTPCnclsRatio2] = (pinfo->TPCCrossedRows()>0 ? Float_t(pinfo->TPCncls())/Float_t(pinfo->TPCCrossedRows()) : 0.0);

  if(fgUsedVars[kTPCcrossedRowsOverFindableClusters]) { 
     if(pinfo->TPCFindableNcls()>0)
       values[kTPCcrossedRowsOverFindableClusters] = Float_t(pinfo->TPCCrossedRows()) / Float_t(pinfo->TPCFindableNcls());
     else 
        values[kTPCcrossedRowsOverFindableClusters] = 0.0;
  }
  if(fgUsedVars[kTPCnclsSharedRatio]) {
     if(pinfo->TPCncls()>0) 
        values[kTPCnclsSharedRatio] = Float_t(pinfo->TPCnclsShared()) / Float_t(pinfo->TPCncls());
     else
        values[kTPCnclsSharedRatio] = 0.0;
  }

  if(fgUsedVars[kTPCnclsRatio3])
    values[kTPCnclsRatio3] = (pinfo->TPCFindableNcls()>0 ? Float_t(pinfo->TPCCrossedRows())/Float_t(pinfo->TPCFindableNcls()) : 0.0);

  values[kTPCnclsF]       = pinfo->TPCFindableNcls();
  values[kTPCnclsShared]  = pinfo->TPCnclsShared();
  values[kTPCcrossedRows] = pinfo->TPCCrossedRows();
  values[kTPCsignal]      = pinfo->TPCsignal();
  values[kTPCsignalN]     = pinfo->TPCsignalN();
  values[kTPCActiveLength] = pinfo->TPCActiveLength();
  values[kTPCGeomLength] = pinfo->TPCGeomLength();
  for(Int_t i=0; i<4; ++i) {
     values[kTPCdEdxQmax+i] = pinfo->TPCdEdxInfoQmax(i);
     values[kTPCdEdxQtot+i] = pinfo->TPCdEdxInfoQtot(i);
     values[kTPCdEdxQmaxOverQtot+i] = ( values[kTPCdEdxQtot+i]>1.0e-7 ? values[kTPCdEdxQmax+i] / values[kTPCdEdxQtot+i] : -999. );
  }
  values[kTPCchi2] = pinfo->TPCchi2();
  if(fgUsedVars[kTPCNclusBitsFired]) values[kTPCNclusBitsFired] = pinfo->TPCClusterMapBitsFired();
  if(fgUsedVars[kTPCclustersPerBit]) {
    Int_t nbits = pinfo->TPCClusterMapBitsFired();
    values[kTPCclustersPerBit] = (nbits>0 ? values[kTPCncls]/Float_t(nbits) : 0.0);
  }

  values[kTOFbeta] = pinfo->TOFbeta();
  values[kTOFdeltaBC] = pinfo->TOFdeltaBC();
  values[kTOFtime] = pinfo->TOFtime();
  values[kTOFdx] = pinfo->TOFdx();
  values[kTOFdz] = pinfo->TOFdz();
  values[kTOFmismatchProbability] = pinfo->TOFmismatchProbab();
  values[kTOFchi2] = pinfo->TOFchi2();

  for(Int_t specie=kElectron; specie<=kProton; ++specie) {
    values[kITSnSig+specie] = pinfo->ITSnSig(specie);
    values[kTPCnSig+specie] = pinfo->TPCnSig(specie);
    values[kTOFnSig+specie] = pinfo->TOFnSig(specie);
    values[kBayes+specie]   = pinfo->GetBayesProb(specie);
  }
  if(fgUsedVars[kTPCnSigCorrected+kElectron] && fgTPCelectronCentroidMap && fgTPCelectronWidthMap) {
     Int_t binX = fgTPCelectronCentroidMap->GetXaxis()->FindBin(values[fgVarDependencyX]);
     if(binX==0) binX = 1;
     if(binX==fgTPCelectronCentroidMap->GetXaxis()->GetNbins()+1) binX -= 1;
     Int_t binY = fgTPCelectronCentroidMap->GetYaxis()->FindBin(values[fgVarDependencyY]);
     if(binY==0) binY=1;
     if(binY==fgTPCelectronCentroidMap->GetYaxis()->GetNbins()+1) binY -= 1;
     Float_t centroid = fgTPCelectronCentroidMap->GetBinContent(binX, binY);
     Float_t width = fgTPCelectronWidthMap->GetBinContent(binX, binY);
     if(TMath::Abs(width)<1.0e-6) width = 1.;
     values[kTPCnSigCorrected+kElectron] = (values[kTPCnSig+kElectron] - centroid)/width;
  }
  
  
  if(fgUsedVars[kTPCnSigCorrected+kElectron] && fgTPCpidCalibCentroid[0] && fgTPCpidCalibWidth[0] && fgTPCpidCalibStatus[0]) {
     Int_t bin[4];
     for(Int_t i=0; i<4; i++) {
        bin[i] = fgTPCpidCalibCentroid[0]->GetAxis(i)->FindBin(values[fgTPCpidCalibVars[i]]);
        if(bin[i]==0) bin[i] = 1;
        if(bin[i]==fgTPCpidCalibCentroid[0]->GetAxis(i)->GetNbins()+1) bin[i] -= 1;
     }
     Int_t status = fgTPCpidCalibStatus[0]->GetBinContent(bin);
     if(status<0) {
        values[kTPCnSigCorrected+kElectron] = -999.0;
     }
     else {
        Float_t centroid = fgTPCpidCalibCentroid[0]->GetBinContent(bin);
        Float_t width = fgTPCpidCalibWidth[0]->GetBinContent(bin);
        values[kTPCnSigCorrected+kElectron] = (values[kTPCnSig+kElectron] - centroid)/width;
     }
  }
  if(fgUsedVars[kTPCnSigCorrected+kPion] && fgTPCpidCalibCentroid[1] && fgTPCpidCalibWidth[1] && fgTPCpidCalibStatus[1]) {
     Int_t bin[4];
     for(Int_t i=0; i<4; i++) {
        bin[i] = fgTPCpidCalibCentroid[1]->GetAxis(i)->FindBin(values[fgTPCpidCalibVars[i]]);
        if(bin[i]==0) bin[i] = 1;
        if(bin[i]==fgTPCpidCalibCentroid[1]->GetAxis(i)->GetNbins()+1) bin[i] -= 1;
     }
     Int_t status = fgTPCpidCalibStatus[1]->GetBinContent(bin);
     if(status<0) {
        values[kTPCnSigCorrected+kPion] = -999.0;
     }
     else {
        Float_t centroid = fgTPCpidCalibCentroid[1]->GetBinContent(bin);
        Float_t width = fgTPCpidCalibWidth[1]->GetBinContent(bin);
        values[kTPCnSigCorrected+kPion] = (values[kTPCnSig+kPion] - centroid)/width;
     }
  }
  if(fgUsedVars[kTPCnSigCorrected+kProton] && fgTPCpidCalibCentroid[2] && fgTPCpidCalibWidth[2] && fgTPCpidCalibStatus[2]) {
     Int_t bin[4];
     for(Int_t i=0; i<4; i++) {
        bin[i] = fgTPCpidCalibCentroid[2]->GetAxis(i)->FindBin(values[fgTPCpidCalibVars[i]]);
        if(bin[i]==0) bin[i] = 1;
        if(bin[i]==fgTPCpidCalibCentroid[2]->GetAxis(i)->GetNbins()+1) bin[i] -= 1;
     }
     Int_t status = fgTPCpidCalibStatus[2]->GetBinContent(bin);
     if(status<0) {
        values[kTPCnSigCorrected+kProton] = -999.0;
     }
     else {
        Float_t centroid = fgTPCpidCalibCentroid[2]->GetBinContent(bin);
        Float_t width = fgTPCpidCalibWidth[2]->GetBinContent(bin);
        values[kTPCnSigCorrected+kProton] = (values[kTPCnSig+kProton] - centroid)/width;
     }        
  }

  values[kTRDpidProbabilitiesLQ1D]   = pinfo->TRDpidLQ1D(0);
  values[kTRDpidProbabilitiesLQ1D+1] = pinfo->TRDpidLQ1D(1);
  values[kTRDpidProbabilitiesLQ2D]   = pinfo->TRDpidLQ2D(0);
  values[kTRDpidProbabilitiesLQ2D+1] = pinfo->TRDpidLQ2D(1);
  values[kTRDntracklets]    = pinfo->TRDntracklets(0);
  values[kTRDntrackletsPID] = pinfo->TRDntracklets(1);

  // TRD GTU online tracks
  values[kTRDGTUtracklets]   = pinfo->TRDGTUtracklets();
  values[kTRDGTUlayermask]   = pinfo->TRDGTUlayermask();
  values[kTRDGTUpt]          = pinfo->TRDGTUpt();
  values[kTRDGTUsagitta]     = pinfo->TRDGTUsagitta();
  values[kTRDGTUPID]         = pinfo->TRDGTUPID();

  FillTrackingStatus(pinfo,values);
  //FillTrackingFlags(pinfo,values);

  if(fgUsedVars[kPtMC]) values[kPtMC] = pinfo->PtMC();
  if(fgUsedVars[kPMC]) values[kPMC] = pinfo->PMC();
  values[kPxMC] = pinfo->MCmom(0);
  values[kPyMC] = pinfo->MCmom(1);
  values[kPzMC] = pinfo->MCmom(2);
  if(fgUsedVars[kThetaMC]) values[kThetaMC] = pinfo->ThetaMC();
  if(fgUsedVars[kEtaMC]) values[kEtaMC] = pinfo->EtaMC();
  if(fgUsedVars[kPhiMC]) values[kPhiMC] = pinfo->PhiMC();
  //TODO: add also the massMC and RapMC   
  values[kPdgMC] = pinfo->MCPdg(0);
  values[kPdgMC+1] = pinfo->MCPdg(1);
  values[kPdgMC+2] = pinfo->MCPdg(2);
  values[kPdgMC+3] = pinfo->MCPdg(3);
  
  if(fgUsedVars[kRap] && pinfo->IsMCKineParticle())  {
     if(pinfo->MCPdg(0)==443) values[kRap] = p->Rapidity(fgkPairMass[AliReducedPairInfo::kJpsiToEE]);
     if(pinfo->MCPdg(0)==100443) values[kRap] = p->Rapidity(fgkPairMass[AliReducedPairInfo::kPsi2SToEE]);
     if(TMath::Abs(pinfo->MCPdg(0))==11) values[kRap] = p->Rapidity(fgkParticleMass[AliReducedVarManager::kElectron]);
  }
  if(fgUsedVars[kRapAbs] && pinfo->IsMCKineParticle())  {
    if(pinfo->MCPdg(0)==443) values[kRapAbs] = TMath::Abs(p->Rapidity(fgkPairMass[AliReducedPairInfo::kJpsiToEE]));
    if(pinfo->MCPdg(0)==100443) values[kRapAbs] = TMath::Abs(p->Rapidity(fgkPairMass[AliReducedPairInfo::kPsi2SToEE]));
    if(TMath::Abs(pinfo->MCPdg(0))==11) values[kRapAbs] = TMath::Abs(p->Rapidity(fgkParticleMass[AliReducedVarManager::kElectron]));
  }
}

//_________________________________________________________________
void AliReducedVarManager::FillClusterMatchedTrackInfo(AliReducedBaseTrack* p, Float_t* values, TList* clusterList/*=0x0*/, AliReducedCaloClusterTrackMatcher* matcher/*=0x0*/) {
  //
  // fill calorimeter cluster matched track info
  //
  if(p->IsA()!=TRACK::Class()) return;
  TRACK* pinfo = (TRACK*)p;

  if (!fgUsedVars[kEMCALmatchedEnergy] &&
      !fgUsedVars[kEMCALmatchedEOverP] &&
      !fgUsedVars[kEMCALmatchedM02] &&
      !fgUsedVars[kEMCALmatchedM20] &&
      !fgUsedVars[kEMCALmatchedClusterId]) return;

  if (fgEvent && (fgEvent->IsA()==EVENT::Class())) {
    CLUSTER* cluster = NULL;
    if (clusterList) {
      for (Int_t i=0; i<clusterList->GetEntries(); ++i) {
        cluster = (AliReducedCaloClusterInfo*)clusterList->At(i);
        if (cluster->ClusterID()==pinfo->CaloClusterId()) break;
        cluster = NULL;
      }
    } else {
      cluster = ((EVENT*)fgEvent)->GetCaloClusterFromID(pinfo->CaloClusterId());
    }
    // track matcher:
    Float_t deltaPhi  = -9999.;
    Float_t deltaEta  = -9999.;
    Float_t dist      = -9999.;
    if (matcher) {
      if (!matcher->IsClusterMatchedToTrack(pinfo, cluster, deltaPhi, deltaEta, dist)) cluster = NULL;
    }
    if (fgUsedVars[kEMCALmatchedClusterId])       values[kEMCALmatchedClusterId]      = (cluster ? pinfo->CaloClusterId() : -9999.);
    if (fgUsedVars[kEMCALmatchedEnergy])          values[kEMCALmatchedEnergy]         = (cluster ? cluster->Energy() : -9999.);
    if (fgUsedVars[kEMCALmatchedM02])             values[kEMCALmatchedM02]            = (cluster ? cluster->M02() : -9999.);
    if (fgUsedVars[kEMCALmatchedM20])             values[kEMCALmatchedM20]            = (cluster ? cluster->M20() : -9999.);
    if (fgUsedVars[kEMCALmatchedNCells])          values[kEMCALmatchedNCells]         = (cluster ? cluster->NCells() : -9999.);
    if (fgUsedVars[kEMCALmatchedNMatchedTracks])  values[kEMCALmatchedNMatchedTracks] = (cluster ? cluster->NMatchedTracks() : -9999.);
    if (fgUsedVars[kEMCALmatchedNSigmaElectron])  values[kEMCALmatchedNSigmaElectron] = (cluster ? pinfo->EMCALnSigEle() : -9999.);
    if (fgUsedVars[kEMCALmatchedDeltaPhi])        values[kEMCALmatchedDeltaPhi]       = (cluster ? deltaPhi : -9999.);
    if (fgUsedVars[kEMCALmatchedDeltaEta])        values[kEMCALmatchedDeltaEta]       = (cluster ? deltaEta : -9999.);
    if (fgUsedVars[kEMCALmatchedDistance])        values[kEMCALmatchedDistance]       = (cluster ? dist : -9999.);
    if (fgUsedVars[kEMCALmatchedEOverP]) {
      Float_t               mom = 0.0;
      if (pinfo->PonCalo()) mom = pinfo->PonCalo();
      else                  mom = pinfo->P();
      values[kEMCALmatchedEOverP] = (TMath::Abs(mom)>1.e-8 && cluster ? values[kEMCALmatchedEnergy]/mom : -9999.);
    }
  }
}

//_________________________________________________________________
void AliReducedVarManager::FillCaloClusterInfo(CLUSTER* cl, Float_t* values) {
  //
  // Fill calorimeter cluster information
  //
  values[kEMCALclusterEnergy] = cl->Energy();
  values[kEMCALclusterDx] = cl->Dx();
  values[kEMCALclusterDz] = cl->Dz();
  values[kEMCALdetector] = (cl->IsEMCAL() ? CLUSTER::kEMCAL : CLUSTER::kPHOS);
  values[kEMCALm02] = cl->M02();
  values[kEMCALm20] = cl->M20();
  values[kEMCALdispersion] = cl->Dispersion();
  values[kEMCALnCells] = cl->NCells();
  values[kEMCALnMatchedTracks] = cl->NMatchedTracks();
  TVector3 clusterVector(cl->X(), cl->Y(), cl->Z());
  Float_t phiCluster = clusterVector.Phi();
  if (phiCluster<0) phiCluster += 2*TMath::Pi();
  values[kEMCALclusterPhi] = phiCluster;
  values[kEMCALclusterEta] = clusterVector.Eta();
}

//_________________________________________________________________
void AliReducedVarManager::GetLegMassAssumption(Int_t id, Float_t& m1, Float_t& m2) {
  //
  // Get the mass assumption for the pair legs depending on its candidate ID
  //
  switch (id) {
    case PAIR::kK0sToPiPi :
      m1 = fgkParticleMass[kPion]; m2 = fgkParticleMass[kPion];
      break;
    case PAIR::kPhiToKK :
      m1 = fgkParticleMass[kKaon]; m2 = fgkParticleMass[kKaon];
      break;
    case PAIR::kLambda0ToPPi :
      m1 = fgkParticleMass[kProton]; m2 = fgkParticleMass[kPion];
      break;
    case PAIR::kALambda0ToPPi :
      m1 = fgkParticleMass[kPion]; m2 = fgkParticleMass[kProton];
      break;
    case PAIR::kJpsiToEE :
      m1 = fgkParticleMass[kElectron]; m2 = fgkParticleMass[kElectron];
      break;        
    case PAIR::kUpsilon :
      m1 = fgkParticleMass[kElectron]; m2 = fgkParticleMass[kElectron];
      break;
    case PAIR::kGammaConv :
      m1 = fgkParticleMass[kElectron]; m2 = fgkParticleMass[kElectron];
      break;
    case PAIR::kDplusToK0sPiplus :
      m1 = fgkParticleMass[kKaonZero]; m2 = fgkParticleMass[kPion];
      break;
    case PAIR::kDplusToK0sKplus :
      m1 = fgkParticleMass[kKaonZero]; m2 = fgkParticleMass[kKaon];
      break;
    case PAIR::kDplusToPhiPiplus :
      m1 = fgkParticleMass[kPhiMeson]; m2 = fgkParticleMass[kPion];
      break;	
    case PAIR::kDminusToK0sPiminus :
      m1 = fgkParticleMass[kKaonZero]; m2 = fgkParticleMass[kPion];
      break;	
    case PAIR::kDminusToK0sKminus :
      m1 = fgkParticleMass[kKaonZero]; m2 = fgkParticleMass[kKaon];
      break;	
    case PAIR::kDminusToPhiPiminus :
      m1 = fgkParticleMass[kPhiMeson]; m2 = fgkParticleMass[kPion];
      break;	
    case PAIR::kDzeroToKminusPiplus :
      m1 = fgkParticleMass[kKaon]; m2 = fgkParticleMass[kPion];
      break;	
    case PAIR::kADzeroToKplusPiminus :
      m1 = fgkParticleMass[kKaon]; m2 = fgkParticleMass[kPion];
      break;	
    case PAIR::kDsplusToK0sKplus :
      m1 = fgkParticleMass[kKaonZero]; m2 = fgkParticleMass[kKaon];
      break;	
    case PAIR::kDsminusToK0sKminus :
      m1 = fgkParticleMass[kKaonZero]; m2 = fgkParticleMass[kKaon];
      break;		
      default :
        break;
  }
}

//_________________________________________________________________
void AliReducedVarManager::FillPairInfo(PAIR* p, Float_t* values) {
  //
  // fill pair information
  //
  FillTrackInfo(p, values);
  
  values[kCandidateId]   = p->CandidateId();
  values[kPairType]      = p->PairType();
  values[kPairTypeSPD]      = p->PairTypeSPD();
  values[kPairChisquare] = p->Chi2();
  if(fgUsedVars[kMass]) {
    values[kMass] = p->Mass();
    if(p->CandidateId()==PAIR::kLambda0ToPPi)  values[kMass] = p->Mass(1);
    if(p->CandidateId()==PAIR::kALambda0ToPPi) values[kMass] = p->Mass(2);
    if(p->CandidateId()==PAIR::kGammaConv)     values[kMass] = p->Mass(3);    
  }
  
  Float_t m1 = 0.0; Float_t m2 = 0.0;
  GetLegMassAssumption(p->CandidateId(),m1,m2); 
  
  values[kMassV0]   = p->Mass(0);
  values[kMassV0+1] = p->Mass(1);
  values[kMassV0+2] = p->Mass(2);
  values[kMassV0+3] = p->Mass(3);
  
  if(fgUsedVars[kRap])    values[kRap]              = p->Rapidity();
  if(fgUsedVars[kRapAbs]) values[kRapAbs]           = TMath::Abs(p->Rapidity());
                          values[kPairLxy]          = p->Lxy();
                          values[kPairPointingAngle]= p->PointingAngle();

  // polarization variables
  Bool_t usePolarization=kFALSE;
  if(fgUsedVars[kPairThetaCS] || fgUsedVars[kPairThetaHE] || fgUsedVars[kPairPhiCS] || fgUsedVars[kPairPhiHE])
    usePolarization = kTRUE;
  if(usePolarization)
    GetThetaPhiCM(fgEvent->GetTrack(((AliReducedPairInfo*)p)->LegId(0)), 
		  fgEvent->GetTrack(((AliReducedPairInfo*)p)->LegId(1)), 
		  values[kPairThetaHE], values[kPairPhiHE], values[kPairThetaCS], values[kPairPhiCS], m1, m2);
}


//_________________________________________________________________
void AliReducedVarManager::FillPairInfo(BASETRACK* t1, BASETRACK* t2, Int_t type, Float_t* values) {
  //
  // fill pair information from 2 tracks
  //
  // type - Parameter encoding the resonance type 
  //        This is needed for making a mass assumption on the legs
  //
  PAIR p;
  p.PxPyPz(t1->Px()+t2->Px(), t1->Py()+t2->Py(), t1->Pz()+t2->Pz());
  p.CandidateId(type);
  //p.SetLegIds(t1->TrackId(), t2->TrackId());  
  
  if(t1->Charge()*t2->Charge()<0) p.PairType(1);
  else if(t1->Charge()>0)         p.PairType(0);
  else                            p.PairType(2);
  values[kPairType] = p.PairType();
  values[kPairTypeSPD] = -1.;
  if(t1->IsA()==TRACK::Class() && t2->IsA()==TRACK::Class() ){
   TRACK* ti1=(TRACK*)t1; TRACK* ti2=(TRACK*)t2;
   values[kPairTypeSPD] = ti1->ITSLayerHit(0)+ti2->ITSLayerHit(0);
  }
  values[kCandidateId] = type;
  values[kPairChisquare] = -999.;
  
  Float_t m1 = 0.0; Float_t m2 = 0.0;
  GetLegMassAssumption(type,m1,m2); 
    
  if(fgUsedVars[kMass]) {     
    values[kMass] = m1*m1+m2*m2 + 
                    2.0*(TMath::Sqrt(m1*m1+t1->P()*t1->P())*TMath::Sqrt(m2*m2+t2->P()*t2->P()) - 
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
      values[kMass] = TMath::Sqrt(values[kMass]);
    p.SetMass(values[kMass]);
  }

  if(fgUsedVars[kRap])    values[kRap]    = p.Rapidity();
  if(fgUsedVars[kRapAbs]) values[kRapAbs] = TMath::Abs(p.Rapidity());
  values[kPairLegPt+0] = t1->Pt();
  values[kPairLegPt+1] = t2->Pt();
  values[kPairLegPtSum] = t1->Pt()+t2->Pt();
  
  values[kMassV0]   = -999.0;
  values[kMassV0+1] = -999.0;
  values[kMassV0+2] = -999.0;
  values[kMassV0+3] = -999.0;

  FillTrackInfo(&p, values);
  
  // polarization variables
  Bool_t usePolarization=kFALSE;
  if(fgUsedVars[kPairThetaCS] || fgUsedVars[kPairThetaHE] || fgUsedVars[kPairPhiCS] || fgUsedVars[kPairPhiHE])
    usePolarization = kTRUE;
  if(usePolarization)
    GetThetaPhiCM(t1, t2, values[kPairThetaHE], values[kPairPhiHE], values[kPairThetaCS], values[kPairPhiCS]);
  
  if(fgUsedVars[kDMA] && (t1->IsA()==TRACK::Class()) && (t2->IsA()==TRACK::Class())) {
     TRACK* ti1=(TRACK*)t1; TRACK* ti2=(TRACK*)t2;
     values[kDMA]=TMath::Sqrt((ti1->HelixX()-ti2->HelixX())*(ti1->HelixX()-ti2->HelixX())+(ti1->HelixY()-ti2->HelixY())*(ti1->HelixY()-ti2->HelixY()))-ti1->HelixR()-ti2->HelixR();   
  }
  
  if((fgUsedVars[kPairLegTPCchi2] || fgUsedVars[kPairLegTPCchi2+1]) && (t1->IsA()==TRACK::Class()) && (t2->IsA()==TRACK::Class())) {
     TRACK* ti1=(TRACK*)t1; TRACK* ti2=(TRACK*)t2;
     values[kPairLegTPCchi2] = ti1->TPCchi2();
     values[kPairLegTPCchi2+1] = ti2->TPCchi2();
  }
  if((fgUsedVars[kPairLegITSchi2] || fgUsedVars[kPairLegITSchi2+1]) && (t1->IsA()==TRACK::Class()) && (t2->IsA()==TRACK::Class())) {
     TRACK* ti1=(TRACK*)t1; TRACK* ti2=(TRACK*)t2;
    values[kPairLegITSchi2] = ti1->ITSchi2();
    values[kPairLegITSchi2+1] = ti2->ITSchi2();
  }
  
  if((fgUsedVars[kPseudoProperDecayTime] || fgUsedVars[kPairLxy]) &&  
     (t1->IsA()==TRACK::Class()) && (t2->IsA()==TRACK::Class()) && 
     (fgEvent->IsA()==EVENT::Class())) {
     TRACK* ti1=(TRACK*)t1; 
     TRACK* ti2=(TRACK*)t2;
     AliKFParticle pairKF = BuildKFcandidate(ti1,m1,ti2,m2);
     Double_t errPseudoProperTime2;
     EVENT* eventInfo = (EVENT*)fgEvent;
     AliKFParticle primVtx = BuildKFvertex(eventInfo);
     if(fgUsedVars[kPseudoProperDecayTime]) 
        values[kPseudoProperDecayTime] = pairKF.GetPseudoProperDecayTime(primVtx, fgkPairMass[type], &errPseudoProperTime2);
     if(fgUsedVars[kPairLxy]) values[kPairLxy] =  ( (pairKF.X() - primVtx.X())*p.Px() + (pairKF.Y() - primVtx.Y())*p.Py() )/p.Pt(); // = values[kPseudoProperDecayTime]*(p.Pt()/PAIR::fgkPairMass[type]);
  }
  
  // fill MC information
  if(p.PairType()==1) {
     TRACK* pinfo1 = 0x0;
     if(t1->IsA()==TRACK::Class()) pinfo1 = (TRACK*)t1;
     TRACK* pinfo2 = 0x0;
     if(t2->IsA()==TRACK::Class()) pinfo2 = (TRACK*)t2;
     
     PAIR pMC;
     if(pinfo1 && pinfo2 && !pinfo1->IsMCTruth() && !pinfo2->IsMCTruth())
        pMC.PxPyPz(pinfo1->MCmom(0)+pinfo2->MCmom(0), pinfo1->MCmom(1)+pinfo2->MCmom(1), pinfo1->MCmom(2)+pinfo2->MCmom(2));
     else
        pMC.PxPyPz(t1->Px()+t2->Px(), t1->Py()+t2->Py(), t1->Pz()+t2->Pz());
     pMC.CandidateId(type);
     
     if(fgUsedVars[kPtMC]) values[kPtMC] = pMC.Pt();
     if(fgUsedVars[kPMC]) values[kPMC] = pMC.P();
     values[kPxMC] = pMC.Px();
     values[kPyMC] = pMC.Py();
     values[kPzMC] = pMC.Pz();
     values[kPt_weight]=0.0;
     if(fgUsedVars[kPt_weight]) values[kPt_weight] =  calculate_weight_factor(pMC.Pt(),,values[kCentVZERO]);
     
     if(fgUsedVars[kThetaMC]) values[kThetaMC] = pMC.Theta();
     if(fgUsedVars[kEtaMC]) values[kEtaMC] = pMC.Eta();
     if(fgUsedVars[kPhiMC]) values[kPhiMC] = pMC.Phi();
     if(fgUsedVars[kMassMC]) {
        if(pinfo1 && pinfo2 && !pinfo1->IsMCTruth() && !pinfo2->IsMCTruth())
           values[kMassMC] = m1*m1+m2*m2 + 
              2.0*(TMath::Sqrt(m1*m1+pinfo1->PMC()*pinfo1->PMC())*TMath::Sqrt(m2*m2+pinfo2->PMC()*pinfo2->PMC()) - 
              pinfo1->MCmom(0)*pinfo2->MCmom(0) - pinfo1->MCmom(1)*pinfo2->MCmom(1) - pinfo1->MCmom(2)*pinfo2->MCmom(2));
         else
            values[kMassMC] = m1*m1+m2*m2 + 
               2.0*(TMath::Sqrt(m1*m1+t1->P()*t1->P())*TMath::Sqrt(m2*m2+t2->P()*t2->P()) - 
               t1->Px()*t2->Px() - t1->Py()*t2->Py() - t1->Pz()*t2->Pz());
         if(values[kMassMC]<0.0) {
            cout << "FillPairInfo(track, track, type, values): Warning: Very small squared mass found. "
            << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
            cout << "   massMC2: " << values[kMassMC] << endl;
            values[kMassMC] = 0.0;
         }
         else
         values[kMassMC] = TMath::Sqrt(values[kMassMC]);
     }
     
     // TODO: think about whether to use the PDG mass or the calculated mass from the legs for rapidity
     if(fgUsedVars[kRapMC]) {
       pMC.SetMass(values[kMassMC]);
       values[kRapMC] = pMC.Rapidity();   
     }
     if(fgUsedVars[kRapMCAbs]) values[kRapMCAbs] = TMath::Abs(pMC.Rapidity());
  }

   if( fgUsedVars[kPairPhiV] ){
    // implementation taken from AliDielectronPair.cxx
    Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
    Double_t px2=-9999.,py2=-9999.,pz2=-9999.;

    if (t1->Charge()*t2->Charge() > 0.) { // Like Sign
      if(values[kL3Polarity]<0){ // inverted behaviour
        if(t1->Charge()>0){
          px1 = t1->Px();   py1 = t1->Py();   pz1 = t1->Pz();
          px2 = t2->Px();   py2 = t2->Py();   pz2 = t2->Pz();
        }else{
          px1 = t2->Px();   py1 = t2->Py();   pz1 = t2->Pz();
          px2 = t1->Px();   py2 = t1->Py();   pz2 = t1->Pz();
        }
      }else{
        if(t1->Charge()>0){
          px1 = t2->Px();   py1 = t2->Py();   pz1 = t2->Pz();
          px2 = t1->Px();   py2 = t1->Py();   pz2 = t1->Pz();
        }else{
          px1 = t1->Px();   py1 = t1->Py();   pz1 = t1->Pz();
          px2 = t2->Px();   py2 = t2->Py();   pz2 = t2->Pz();
        }
      }
    }
    else { // Unlike Sign
      if(values[kL3Polarity]>0){ // regular behaviour
        if(t1->Charge()>0){
          px1 = t1->Px();
          py1 = t1->Py();
          pz1 = t1->Pz();

          px2 = t2->Px();
          py2 = t2->Py();
          pz2 = t2->Pz();
        }else{
          px1 = t2->Px();
          py1 = t2->Py();
          pz1 = t2->Pz();

          px2 = t1->Px();
          py2 = t1->Py();
          pz2 = t1->Pz();
        }
      }else{
        if(t1->Charge()>0){
          px1 = t2->Px();
          py1 = t2->Py();
          pz1 = t2->Pz();

          px2 = t1->Px();
          py2 = t1->Py();
          pz2 = t1->Pz();
        }else{
          px1 = t1->Px();
          py1 = t1->Py();
          pz1 = t1->Pz();

          px2 = t2->Px();
          py2 = t2->Py();
          pz2 = t2->Pz();
        }
      }
    }

    Double_t px = px1+px2;
    Double_t py = py1+py2;
    Double_t pz = pz1+pz2;
    Double_t dppair = TMath::Sqrt(px*px+py*py+pz*pz);

    //unit vector of (pep+pem)
    Double_t pl = dppair;
    Double_t ux = px/pl;
    Double_t uy = py/pl;
    Double_t uz = pz/pl;
    Double_t ax = uy/TMath::Sqrt(ux*ux+uy*uy);
    Double_t ay = -ux/TMath::Sqrt(ux*ux+uy*uy);

    //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
    //Double_t ptep = iep->Px()*ax + iep->Py()*ay;
    //Double_t ptem = iem->Px()*ax + iem->Py()*ay;

    Double_t pxep = px1;
    Double_t pyep = py1;
    Double_t pzep = pz1;
    Double_t pxem = px2;
    Double_t pyem = py2;
    Double_t pzem = pz2;

    //vector product of pep X pem
    Double_t vpx = pyep*pzem - pzep*pyem;
    Double_t vpy = pzep*pxem - pxep*pzem;
    Double_t vpz = pxep*pyem - pyep*pxem;
    Double_t vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz);

    //unit vector of pep X pem
    Double_t vx = vpx/vp;
    Double_t vy = vpy/vp;
    Double_t vz = vpz/vp;

    //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
    Double_t wx = uy*vz - uz*vy;
    Double_t wy = uz*vx - ux*vz;
    // by construction, (wx,wy,wz) must be a unit vector.
    // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them
    // should be small if the pair is conversion
    // this function then returns values close to pi!
    Double_t cosPhiV = wx*ax + wy*ay;
    Double_t phiv = TMath::ACos(cosPhiV);
    values[kPairPhiV] = phiv;
  }

  if( fgUsedVars[kPairOpeningAngle] ){
    TVector3 v1(t1->Px(), t1->Py(), t1->Pz());
    TVector3 v2(t2->Px(), t2->Py(), t2->Pz());
    values[kPairOpeningAngle] = v1.Angle(v2);
  }
  
  if(  t1->IsA()==TRACK::Class() && t2->IsA()==TRACK::Class()     ) {
    TRACK* ti1=(TRACK*)t1;
    TRACK* ti2=(TRACK*)t2;

    if( fgUsedVars[kPairDca]   ) values[kPairDca]   = TMath::Sqrt(ti1->DCAxy() * ti1->DCAxy() + ti2->DCAxy() * ti2->DCAxy() + ti1->DCAz() * ti1->DCAz() + ti2->DCAz() * ti2->DCAz());
    if( fgUsedVars[kPairDcaXY] ) values[kPairDcaXY] = TMath::Sqrt( ti1->DCAxy() * ti1->DCAxy() + ti2->DCAxy() * ti2->DCAxy() );
    if( fgUsedVars[kPairDcaZ]  ) values[kPairDcaZ]  = TMath::Sqrt(ti1->DCAz() * ti1->DCAz() + ti2->DCAz() * ti2->DCAz() );

    if( fgUsedVars[kPairDcaSqrt]   ) values[kPairDcaSqrt]   = TMath::Power(ti1->DCAxy() * ti1->DCAxy() + ti2->DCAxy() * ti2->DCAxy() + ti1->DCAz() * ti1->DCAz() + ti2->DCAz() * ti2->DCAz(), 0.25);
    if( fgUsedVars[kPairDcaXYSqrt] ) values[kPairDcaXYSqrt] = TMath::Power( ti1->DCAxy() * ti1->DCAxy() + ti2->DCAxy() * ti2->DCAxy(), 0.25);
    if( fgUsedVars[kPairDcaZSqrt]  ) values[kPairDcaZSqrt]  = TMath::Power(ti1->DCAz() * ti1->DCAz() + ti2->DCAz() * ti2->DCAz(), 0.25);

    if( fgUsedVars[kOpAngDcaPtCorr] ) {
      Float_t a = -1.56316e-03;
      Float_t b =  1.22515e-02;
      Float_t c =  3.39455e-03;
      Float_t d =  1.00681e-01;
      values[kOpAngDcaPtCorr] = values[kPairOpeningAngle] - a - b * values[kPairDcaXYSqrt] - c * values[kOneOverSqrtPt] -  d * values[kPairDcaXYSqrt]  * values[kOneOverSqrtPt];
    }
    if( fgUsedVars[kMassDcaPtCorr] ) {
      Float_t a =  1.87774e-03;
      Float_t b =  4.53156e-02;
      Float_t c = -9.72947e-05;
      Float_t d =  1.83003e-02;
      values[kMassDcaPtCorr] = values[kMass] - a - b * values[kPairDcaXYSqrt] - c * values[kPt] -  d * values[kPairDcaXYSqrt]  * values[kPt];
    }
  }
}


//_________________________________________________________________
void AliReducedVarManager::FillPairInfoME(BASETRACK* t1, BASETRACK* t2, Int_t type, Float_t* values) {
  //
  // Lightweight fill pair information from 2 base track objects.
  // NOTE: Mostly intended for making pairing during event mixing.
  //       FillBaseTrackInfo() is not used to avoid calculation of not needed variables (like flow)
  // type - Parameter encoding the resonance type 
  //        This is needed for making a mass assumption on the legs
  //
  PAIR p;
  p.PxPyPz(t1->Px()+t2->Px(), t1->Py()+t2->Py(), t1->Pz()+t2->Pz());
  p.CandidateId(type);
 
  values[kPairTypeSPD] = -1.;
   if(t1->IsA()==TRACK::Class() && t2->IsA()==TRACK::Class() ){
   TRACK* ti1=(TRACK*)t1; TRACK* ti2=(TRACK*)t2;
   values[kPairTypeSPD] = ti1->ITSLayerHit(0)+ti2->ITSLayerHit(0);
   }
   
  if(t1->Charge()*t2->Charge()<0) p.PairType(1);
  else if(t1->Charge()>0)         p.PairType(0);
  else                            p.PairType(2);
  values[kPairType] = p.PairType();
  values[kCandidateId] = type;
  values[kPairChisquare] = -999.;
  
  Float_t m1 = 0.0; Float_t m2 = 0.0;
  GetLegMassAssumption(type,m1,m2); 
    
  if(fgUsedVars[kMass]) {     
    values[kMass] = m1*m1+m2*m2 + 
                    2.0*(TMath::Sqrt(m1*m1+t1->P()*t1->P())*TMath::Sqrt(m2*m2+t2->P()*t2->P()) - 
                    t1->Px()*t2->Px() - t1->Py()*t2->Py() - t1->Pz()*t2->Pz());
    if(values[kMass]<0.0) {
      cout << "FillPairInfoME(track, track, type, values): Warning: Very small squared mass found. "
           << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
      cout << "   mass2: " << values[kMass] << endl;
      cout << "p1(p,x,y,z): " << t1->P() << ", " << t1->Px() << ", " << t1->Py() << ", " << t1->Pz() << endl;
      cout << "p2(p,x,y,z): " << t2->P() << ", " << t2->Px() << ", " << t2->Py() << ", " << t2->Pz() << endl;
      values[kMass] = 0.0;
    }
    else
      values[kMass] = TMath::Sqrt(values[kMass]);
    p.SetMass(values[kMass]);
  }
  
  values[kPx] = p.Px();
  values[kPy] = p.Py();
  values[kPz] = p.Pz();
  if(fgUsedVars[kPt] || fgUsedVars[kPtSquared]) {
    values[kPt] = p.Pt();
    if(fgUsedVars[kPtSquared]) values[kPtSquared] = values[kPt]*values[kPt];
  }
  values[kPairLegPt] = t1->Pt();
  values[kPairLegPt+1] = t2->Pt();
  values[kPairLegPtSum] = t1->Pt() + t2->Pt();
  if(fgUsedVars[kP])      values[kP]      = p.P();
  if(fgUsedVars[kEta])    values[kEta]    = p.Eta();
  if(fgUsedVars[kRap])    values[kRap]    = p.Rapidity();
  if(fgUsedVars[kRapAbs]) values[kRapAbs] = TMath::Abs(p.Rapidity());
  if(fgUsedVars[kPhi])    values[kPhi]    = p.Phi();
  if(fgUsedVars[kTheta])  values[kTheta]  = p.Theta();
  
  if((fgUsedVars[kPairEff] || fgUsedVars[kOneOverPairEff] || fgUsedVars[kOneOverPairEffSq]) && fgPairEffMap) {
    Int_t binX = 0;
    if (fgEffMapVarDependencyX!=kNothing) {
      binX = fgPairEffMap->GetXaxis()->FindBin(values[fgEffMapVarDependencyX]); //make sure the values[XVar] are filled for EM
      if(binX==0) binX = 1;
      if(binX==fgPairEffMap->GetXaxis()->GetNbins()+1) binX -= 1;
    }
    Int_t binY = 0;
    if (fgEffMapVarDependencyY!=kNothing) {
      binY = fgPairEffMap->GetYaxis()->FindBin(values[fgEffMapVarDependencyY]); //make sure the values[XVar] are filled for EM
      if(binY==0) binY = 1;
      if(binY==fgPairEffMap->GetXaxis()->GetNbins()+1) binY -= 1;
    }
    Int_t binZ = 0;
    if (fgEffMapVarDependencyZ!=kNothing) {
      binZ = fgPairEffMap->GetZaxis()->FindBin(values[fgEffMapVarDependencyZ]); //make sure the values[ZVar] are filled for EM
      if(binZ==0) binZ = 1;
      if(binZ==fgPairEffMap->GetZaxis()->GetNbins()+1) binZ -= 1;
    }

    Float_t                   pairEff = 1.;
    if (binX && binY && binZ) pairEff = fgPairEffMap->GetBinContent(binX, binY, binZ);
    else if (binX && binY)    pairEff = fgPairEffMap->GetBinContent(binX, binY);
    else if (binX)            pairEff = fgPairEffMap->GetBinContent(binX);

    Float_t               oneOverPairEff = 1.;
    if (pairEff > 1.0e-6) oneOverPairEff = 1/pairEff;

    values[kPairEff]          = pairEff;
    values[kOneOverPairEff]   = oneOverPairEff;
    values[kOneOverPairEffSq] = oneOverPairEff*oneOverPairEff;
  }
}


//_________________________________________________________________
void AliReducedVarManager::FillPairInfo(PAIR* t1, BASETRACK* t2, Int_t type, Float_t* values) {
  //
  // fill pair information for a pair with one leg being an AliReducedPair and 
  // the other leg being an AliReducedBaseTrack
  //
  // type - Parameter encoding the resonance type 
  //        This is needed for making a mass assumption on the legs
  //
  PAIR p;
  p.PxPyPz(t1->Px()+t2->Px(), t1->Py()+t2->Py(), t1->Pz()+t2->Pz());
  p.CandidateId(type);
  
  values[kPairType]  = 1;
  values[kCandidateId] = type;
  values[kPairChisquare] = -999.;
  
  Float_t m1 = 0.0; Float_t m2 = 0.0;
  GetLegMassAssumption(type,m1,m2); 
  
  if(fgUsedVars[kMass]) {     
    values[kMass] = m1*m1+m2*m2 + 
                    2.0*(TMath::Sqrt(m1*m1+t1->P()*t1->P())*TMath::Sqrt(m2*m2+t2->P()*t2->P()) - 
                    t1->Px()*t2->Px() - t1->Py()*t2->Py() - t1->Pz()*t2->Pz());
    if(values[kMass]<0.0) {
      cout << "FillPairInfo(pair, track, type, values): Warning: Very small squared mass found. "
           << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
      cout << "   mass2: " << values[kMass] << endl;
      cout << "p1(p,x,y,z): " << t1->P() << ", " << t1->Px() << ", " << t1->Py() << ", " << t1->Pz() << endl;
      cout << "p2(p,x,y,z): " << t2->P() << ", " << t2->Px() << ", " << t2->Py() << ", " << t2->Pz() << endl;
      values[kMass] = 0.0;
    }
    else
      values[kMass] = TMath::Sqrt(values[kMass]);
    p.SetMass(values[kMass]);
  }  
  values[kMassV0] = -1.0;
  values[kMassV0+1] = -1.0;
  values[kMassV0+2] = -1.0;
  values[kMassV0+3] = -1.0;
  
  FillTrackInfo(&p, values);
}


//__________________________________________________________________
void AliReducedVarManager::FillCorrelationInfo(BASETRACK* trig, BASETRACK* assoc, Float_t* values) {
  //
  // fill pair-track correlation information
  // NOTE:  Add here only NEEDED information because this function is called during event mixing in the innermost loop
  //
  if(fgUsedVars[kTriggerPt]) values[kTriggerPt] = trig->Pt();
  if(fgUsedVars[kTriggerRap] && (trig->IsA()==PAIR::Class())) 	  values[kTriggerRap]     = ((PAIR*)trig)->Rapidity();
  if(fgUsedVars[kTriggerRapAbs] && (trig->IsA()==PAIR::Class()))  values[kTriggerRapAbs]  = TMath::Abs(((PAIR*)trig)->Rapidity());
  if(fgUsedVars[kTriggerPseudoProperDecayTime] && (trig->IsA()==PAIR::Class())) values[kTriggerPseudoProperDecayTime] = ((PAIR*)trig)->PsProper();
  if(fgUsedVars[kTriggerPairTypeSPD] && (trig->IsA()==PAIR::Class()))           values[kTriggerPairTypeSPD]           = ((PAIR*)trig)->PairTypeSPD();
  if(fgUsedVars[kAssociatedPt]) values[kAssociatedPt] = assoc->Pt();
  if(fgUsedVars[kAssociatedEta]) values[kAssociatedEta] = assoc->Eta();
  if(fgUsedVars[kAssociatedPhi]) values[kAssociatedPhi] = assoc->Phi();

  // associated pT / transverse trigger gamma
  if (trig->IsA()==PAIR::Class() &&
      (fgUsedVars[kAssociatedPtOverTriggerGammaT] || fgUsedVars[kTriggerGammaT])) {

    // NOTE:  only interested in transverse beta (gamma) -> eta is set to zero for beta vector calculation
    //        gives same result as 'manual' calculation, i.e.:  betaT  = pT / (m^2 + pT^2)
    //                                                          gammaT = 1 / sqrt(1-betaT^2) = sqrt(1 + pT^2/m^2)
    TLorentzVector trigVec;
    trigVec.SetPtEtaPhiM(trig->Pt(), 0.0, trig->Phi(), 3.096916); //NOTE: J/psi mass from PDG
    TVector3 betaVec = trigVec.BoostVector();

    Float_t betaT   = betaVec.Mag();
    Float_t gammaT  = 1./TMath::Sqrt(1-betaT*betaT);

    if (fgUsedVars[kTriggerGammaT])                 values[kTriggerGammaT]                  = gammaT;
    if (fgUsedVars[kAssociatedPtOverTriggerGammaT]) values[kAssociatedPtOverTriggerGammaT]  = assoc->Pt()/gammaT;
  }

  // values after boost of hadrons to pair rest frame
  if (trig->IsA()==PAIR::Class() &&
      (fgUsedVars[kDeltaPhiBoosted] || fgUsedVars[kDeltaPhiSymBoosted] || fgUsedVars[kDeltaThetaBoosted] || fgUsedVars[kDeltaEtaBoosted] ||
       fgUsedVars[kDeltaEtaAbsBoosted] || fgUsedVars[kAssociatedPtBoosted] || fgUsedVars[kAssociatedEtaBoosted] || fgUsedVars[kAssociatedPhiBoosted])) {

    // get boost vector
    TLorentzVector trigVec;
    trigVec.SetPtEtaPhiM(trig->Pt(), trig->Eta(), trig->Phi(), 3.096916); //NOTE: J/psi mass from PDG
    TVector3 boostVec = trigVec.BoostVector();

    // fill TLorentzVector for associated track
    TLorentzVector assocVec;
    assocVec.SetPtEtaPhiM(assoc->Pt(), assoc->Eta(), assoc->Phi(), 0.13957061); // NOTE: pion mass from PDG
    assocVec.Boost(-boostVec);

    if(fgUsedVars[kAssociatedPtBoosted]) values[kAssociatedPtBoosted] = assocVec.Pt();
    if(fgUsedVars[kAssociatedEtaBoosted]) values[kAssociatedEtaBoosted] = assocVec.Eta();
    if(fgUsedVars[kAssociatedPhiBoosted]) values[kAssociatedPhiBoosted] = assocVec.Phi();

    if(fgUsedVars[kDeltaPhiBoosted]) {
      Double_t delta = trig->Phi() - assocVec.Phi();
      if(delta>3.0/2.0*TMath::Pi()) delta -= 2.0*TMath::Pi();
      if(delta<-0.5*TMath::Pi()) delta += 2.0*TMath::Pi();
      values[kDeltaPhiBoosted] = delta;
    }
    if(fgUsedVars[kDeltaPhiSymBoosted]) {
      Double_t delta = TMath::Abs(trig->Phi() - assocVec.Phi());
      if(delta>TMath::Pi()) delta = 2*TMath::Pi()-delta;
      values[kDeltaPhiSymBoosted] = delta;
    }

    if(fgUsedVars[kDeltaThetaBoosted]) values[kDeltaThetaBoosted] = trig->Theta() - assocVec.Theta();

    if(fgUsedVars[kDeltaEtaBoosted])     values[kDeltaEtaBoosted]     = trig->Eta() - assocVec.Eta();
    if(fgUsedVars[kDeltaEtaAbsBoosted])  values[kDeltaEtaAbsBoosted]  = TMath::Abs(trig->Eta() - assocVec.Eta());
  }

  if(fgUsedVars[kDeltaPhi]) {
    Double_t delta = trig->Phi() - assoc->Phi();
    if(delta>3.0/2.0*TMath::Pi()) delta -= 2.0*TMath::Pi();
    if(delta<-0.5*TMath::Pi()) delta += 2.0*TMath::Pi();
    values[kDeltaPhi] = delta;
  }
  if(fgUsedVars[kDeltaPhiSym]) {
    Double_t delta = TMath::Abs(trig->Phi() - assoc->Phi());
    if(delta>TMath::Pi()) delta = 2*TMath::Pi()-delta;
    values[kDeltaPhiSym] = delta;
  }

  if(fgUsedVars[kDeltaTheta]) values[kDeltaTheta] = trig->Theta() - assoc->Theta();
  
  if(fgUsedVars[kDeltaEta])     values[kDeltaEta]     = trig->Eta() - assoc->Eta();
  if(fgUsedVars[kDeltaEtaAbs])  values[kDeltaEtaAbs]  = TMath::Abs(trig->Eta() - assoc->Eta());
  if(fgUsedVars[kMass] && (trig->IsA()==PAIR::Class())) values[kMass] = ((PAIR*)trig)->Mass();

  // J/psi efficiency variables
  if ((fgUsedVars[kTriggerEff] || fgUsedVars[kOneOverTriggerEff]) && fgPairEffMap) {
    Int_t binX = 0;
    if (fgEffMapVarDependencyXCorr!=kNothing) {
      binX = fgPairEffMap->GetXaxis()->FindBin(values[fgEffMapVarDependencyXCorr]);
      if(binX==0) binX = 1;
      if(binX==fgPairEffMap->GetXaxis()->GetNbins()+1) binX -= 1;
    }
    Int_t binY = 0;
    if (fgEffMapVarDependencyYCorr!=kNothing) {
      binY = fgPairEffMap->GetYaxis()->FindBin(values[fgEffMapVarDependencyYCorr]);
      if(binY==0) binY = 1;
      if(binY==fgPairEffMap->GetXaxis()->GetNbins()+1) binY -= 1;
    }
    Int_t binZ = 0;
    if (fgEffMapVarDependencyZCorr!=kNothing) {
      binZ = fgPairEffMap->GetZaxis()->FindBin(values[fgEffMapVarDependencyZCorr]);
      if(binZ==0) binZ = 1;
      if(binZ==fgPairEffMap->GetZaxis()->GetNbins()+1) binZ -= 1;
    }

    Float_t                   pairEff = 1.;
    if (binX && binY && binZ) pairEff = fgPairEffMap->GetBinContent(binX, binY, binZ);
    else if (binX && binY)    pairEff = fgPairEffMap->GetBinContent(binX, binY);
    else if (binX)            pairEff = fgPairEffMap->GetBinContent(binX);

    Float_t               oneOverPairEff = 1.;
    if (pairEff > 1.0e-6) oneOverPairEff = 1/pairEff;

    values[kTriggerEff]         = pairEff;
    values[kOneOverTriggerEff]  = oneOverPairEff;
  }

  // hadron efficiency variables
  if ((fgUsedVars[kAssocHadronEff] || fgUsedVars[kOneOverAssocHadronEff]) && fgAssocHadronEffMap) {
    Int_t binX = 0;
    if (fgAssocHadronEffMapVarDependencyX!=kNothing) {
      binX = fgAssocHadronEffMap->GetXaxis()->FindBin(values[fgAssocHadronEffMapVarDependencyX]);
      if(binX==0) binX = 1;
      if(binX==fgAssocHadronEffMap->GetXaxis()->GetNbins()+1) binX -= 1;
    }
    Int_t binY = 0;
    if (fgAssocHadronEffMapVarDependencyY!=kNothing) {
      binY = fgAssocHadronEffMap->GetYaxis()->FindBin(values[fgAssocHadronEffMapVarDependencyY]);
      if(binY==0) binY = 1;
      if(binY==fgAssocHadronEffMap->GetXaxis()->GetNbins()+1) binY -= 1;
    }
    Int_t binZ = 0;
    if (fgAssocHadronEffMapVarDependencyZ!=kNothing) {
      binZ = fgAssocHadronEffMap->GetZaxis()->FindBin(values[fgAssocHadronEffMapVarDependencyZ]);
      if(binZ==0) binZ = 1;
      if(binZ==fgAssocHadronEffMap->GetZaxis()->GetNbins()+1) binZ -= 1;
    }

    Float_t                   hadronEff = 1.;
    if (binX && binY && binZ) hadronEff = fgAssocHadronEffMap->GetBinContent(binX, binY, binZ);
    else if (binX && binY)    hadronEff = fgAssocHadronEffMap->GetBinContent(binX, binY);
    else if (binX)            hadronEff = fgAssocHadronEffMap->GetBinContent(binX);

    Float_t                 oneOverHadronEff = 1.;
    if (hadronEff > 1.0e-6) oneOverHadronEff = 1./hadronEff;

    values[kAssocHadronEff]         = hadronEff;
    values[kOneOverAssocHadronEff]  = oneOverHadronEff;
  }

  // J/psi x hadron efficiency variables
  if ((fgUsedVars[kTriggerEffTimesAssocHadronEff] || fgUsedVars[kOneOverTriggerEffTimesAssocHadronEff]) &&
      fgPairEffMap && fgAssocHadronEffMap) {
    values[kTriggerEffTimesAssocHadronEff]        = values[kTriggerEff]*values[kAssocHadronEff];
    values[kOneOverTriggerEffTimesAssocHadronEff] = values[kOneOverTriggerEff]*values[kOneOverAssocHadronEff];
  }
}


//________________________________________________________________
Double_t AliReducedVarManager::DeltaPhi(Double_t phi1, Double_t phi2) {
  //
  // compute the delta of two angles defined in the (-pi,+pi) interval
  //
  Double_t delta = phi1-phi2;
  //if(delta>2.0*TMath::Pi()) delta -= 2.0*TMath::Pi();
  //if(delta<0.0) delta += 2.0*TMath::Pi();
  /*Double_t delta = phi2;
  if(phi2<0.0) delta += 2.0*TMath::Pi();
  delta = phi1-delta;
  if(delta>TMath::Pi()) delta = delta - 2.0*TMath::Pi();
  if(delta<-1.*TMath::Pi()) delta = 2.0*TMath::Pi() + delta;
  */
  return delta;
}


//____________________________________________________________________________________
void AliReducedVarManager::GetThetaPhiCM(BASETRACK* leg1, BASETRACK* leg2,
                                    Float_t &thetaHE, Float_t &phiHE, 
                                    Float_t &thetaCS, Float_t &phiCS,
				    Float_t leg1Mass /*=gkParticleMass[kElectron]*/, Float_t leg2Mass /*=gkParticleMass[kElectron]*/)
{
  //
  // Calculate theta and phi in helicity and Collins-Soper coordinate frame
  //
  if(!leg1||!leg2) {cout<<"AliReducedVarManager::GetThetaPhiCM:  base leg doesn't exist"<<endl; return;}
  Double_t pxyz1[3]={leg1->Px(),leg1->Py(),leg1->Pz()};
  Double_t pxyz2[3]={leg2->Px(),leg2->Py(),leg2->Pz()};
    
  TLorentzVector projMom(0.,0.,-fgBeamMomentum,TMath::Sqrt(fgBeamMomentum*fgBeamMomentum+fgkParticleMass[kProton]*fgkParticleMass[kProton]));
  TLorentzVector targMom(0.,0., fgBeamMomentum,TMath::Sqrt(fgBeamMomentum*fgBeamMomentum+fgkParticleMass[kProton]*fgkParticleMass[kProton]));
  
  // first & second daughter 4-mom
  TLorentzVector p1Mom(pxyz1[0],pxyz1[1],pxyz1[2],
                       TMath::Sqrt(pxyz1[0]*pxyz1[0]+pxyz1[1]*pxyz1[1]+pxyz1[2]*pxyz1[2]+leg1Mass*leg1Mass));
  TLorentzVector p2Mom(pxyz2[0],pxyz2[1],pxyz2[2],
                       TMath::Sqrt(pxyz2[0]*pxyz2[0]+pxyz2[1]*pxyz2[1]+pxyz2[2]*pxyz2[2]+leg2Mass*leg2Mass));
  // J/Psi 4-momentum vector
  TLorentzVector motherMom=p1Mom+p2Mom;
  
  // boost all the 4-mom vectors to the mother rest frame
  TVector3 beta = (-1.0/motherMom.E())*motherMom.Vect();
  p1Mom.Boost(beta);
  p2Mom.Boost(beta);
  projMom.Boost(beta);
  targMom.Boost(beta);
  
  // x,y,z axes
  TVector3 zAxisHE = (motherMom.Vect()).Unit();
  TVector3 zAxisCS = ((projMom.Vect()).Unit()-(targMom.Vect()).Unit()).Unit();
  TVector3 yAxis = ((projMom.Vect()).Cross(targMom.Vect())).Unit();
  TVector3 xAxisHE = (yAxis.Cross(zAxisHE)).Unit();
  TVector3 xAxisCS = (yAxis.Cross(zAxisCS)).Unit();
  
  // fill theta and phi
  if(leg1->Charge()>0){
    thetaHE = zAxisHE.Dot((p1Mom.Vect()).Unit());
    thetaCS = zAxisCS.Dot((p1Mom.Vect()).Unit());
    phiHE   = TMath::ATan2((p1Mom.Vect()).Dot(yAxis), (p1Mom.Vect()).Dot(xAxisHE));
    phiCS   = TMath::ATan2((p1Mom.Vect()).Dot(yAxis), (p1Mom.Vect()).Dot(xAxisCS));
  } else {
    thetaHE = zAxisHE.Dot((p2Mom.Vect()).Unit());
    thetaCS = zAxisCS.Dot((p2Mom.Vect()).Unit());
    phiHE   = TMath::ATan2((p2Mom.Vect()).Dot(yAxis), (p2Mom.Vect()).Dot(xAxisHE));
    phiCS   = TMath::ATan2((p2Mom.Vect()).Dot(yAxis), (p2Mom.Vect()).Dot(xAxisCS));
  }
}

//__________________________________________________________________
void AliReducedVarManager::PrintTrackFlags(TRACK* track) {
  //
  // Print the track flags
  //
  for(UShort_t i=0;i<64;++i)
    cout << (track->TestFlag(i) ? 1 : 0) << flush;
}

//__________________________________________________________________
void AliReducedVarManager::PrintBits(ULong_t mask, Int_t maxBit) {
  //
  // Bit-wise print of mask
  //
  for(UShort_t i=0;i<maxBit;++i)
    cout << (mask&(ULong_t(1)<<i) ? 1 : 0) << flush;
}

//__________________________________________________________________
void AliReducedVarManager::PrintBits(UInt_t mask, Int_t maxBit) {
  //
  // Bit-wise print of mask
  //
  for(UShort_t i=0;i<maxBit;++i)
    cout << (mask&(ULong_t(1)<<i) ? 1 : 0) << flush;
}

//__________________________________________________________________
void AliReducedVarManager::SetDefaultVarNames() {
  //
  // Set default variable names
  //
  for(Int_t ivar=0; ivar<kNVars; ++ivar) {
    fgVariableNames[ivar] = "DEFAULT NOT DEFINED"; fgVariableUnits[ivar] = "n/a";
  }
  
  fgVariableNames[kEventTag]             = "Event tag";                       fgVariableUnits[kEventTag]             = "";
  fgVariableNames[kEventNumberInFile]    = "Event no. in ESD file";           fgVariableUnits[kEventNumberInFile]    = "";
  fgVariableNames[kL0TriggerInput]       = "  ";               fgVariableUnits[kL0TriggerInput]       = "";
  fgVariableNames[kL1TriggerInput]       = "  ";               fgVariableUnits[kL1TriggerInput]       = "";
  fgVariableNames[kL2TriggerInput]       = "  ";               fgVariableUnits[kL2TriggerInput]       = "";
  fgVariableNames[kL0TriggerInput2]       = "  ";               fgVariableUnits[kL0TriggerInput2]       = "";
  fgVariableNames[kL1TriggerInput2]       = "  ";               fgVariableUnits[kL1TriggerInput2]       = "";
  fgVariableNames[kL2TriggerInput2]       = "  ";               fgVariableUnits[kL2TriggerInput2]       = "";
  fgVariableNames[kRunNo]                = "Run number";                      fgVariableUnits[kRunNo]                = "";
  fgVariableNames[kRunID]                = "Run number";                      fgVariableUnits[kRunID]                = "";
  fgVariableNames[kLHCFillNumber]        = "LHC fill number";                 fgVariableUnits[kLHCFillNumber]        = ""; 
  fgVariableNames[kBeamEnergy]           = "Beam energy";                     fgVariableUnits[kBeamEnergy]           = "GeV";
  fgVariableNames[kInstLumi]             = "Instantaneous luminosity";        fgVariableUnits[kInstLumi]             = "Hz/mb";   
  fgVariableNames[kDetectorMask]         = "Detector mask";                   fgVariableUnits[kDetectorMask]         = "";
  fgVariableNames[kNumberOfDetectors]    = "Number of active detectors";      fgVariableUnits[kNumberOfDetectors]    = "";  
  fgVariableNames[kDipolePolarity]       = "Dipole magnet polarity";          fgVariableUnits[kDipolePolarity]       = "";
  fgVariableNames[kL3Polarity]           = "L3 magnet polarity";              fgVariableUnits[kL3Polarity]           = "";
  fgVariableNames[kTotalLuminosity]  = "Run luminosity";                    fgVariableUnits[kTotalLuminosity]  = "";
  fgVariableNames[kBeamIntensity0]       = "Beam 0 intensity";                fgVariableUnits[kBeamIntensity0]       = "";
  fgVariableNames[kBeamIntensity1]       = "Beam 1 intensity";                fgVariableUnits[kBeamIntensity1]       = "";
  fgVariableNames[kRunTimeStart]         = "Run start time";                fgVariableUnits[kRunTimeStart] = "sec";
  fgVariableNames[kRunTimeEnd]         = "Run end time";                fgVariableUnits[kRunTimeEnd] = "sec";
  fgVariableNames[kBC]                   = "Bunch crossing";                  fgVariableUnits[kBC]                   = "";
  fgVariableNames[kTimeStamp]            = "Time stamp";                      fgVariableUnits[kTimeStamp]            = "";
  fgVariableNames[kTimeRelativeSOR]   = "Event time from SOR";      fgVariableUnits[kTimeRelativeSOR]  = "min";
  fgVariableNames[kTimeRelativeSORfraction] = "Event time from SOR";   fgVariableUnits[kTimeRelativeSORfraction] = "fraction of total run duration";
  fgVariableNames[kEventType]            = "Event type";                      fgVariableUnits[kEventType]            = "";
  fgVariableNames[kTriggerMask]          = "Trigger mask";                    fgVariableUnits[kTriggerMask]          = "";
  fgVariableNames[kOnlineTrigger]        = "Online trigger";                  fgVariableUnits[kOnlineTrigger]        = "";
  fgVariableNames[kOnlineTriggerFired]   = "  ";            fgVariableUnits[kOnlineTriggerFired]   = "";
  fgVariableNames[kOnlineTriggerFired2]  = "  ";           fgVariableUnits[kOnlineTriggerFired2]  = "";
  fgVariableNames[kIsPhysicsSelection]   = "Physics selection ON";            fgVariableUnits[kIsPhysicsSelection]   = "";
  fgVariableNames[kIsSPDPileup]          = "SPD pileup ON";                   fgVariableUnits[kIsSPDPileup]          = "";
  fgVariableNames[kIsSPDPileup5]          = "SPD pileup (5 contributors) ON";                   fgVariableUnits[kIsSPDPileup5]          = "";
  fgVariableNames[kIsPileupMV]          = "MV pileup ON";                   fgVariableUnits[kIsPileupMV]          = "";
  fgVariableNames[kIsSPDPileupMultBins]  = "SPD pileup multiplicity bins ON"; fgVariableUnits[kIsSPDPileupMultBins]  = "";
  fgVariableNames[kNSPDpileups]          = "Number of SPD pileup events";     fgVariableUnits[kNSPDpileups]          = "";
  fgVariableNames[kNTrackPileups]        = "Number of track pileup events";   fgVariableUnits[kNTrackPileups]        = "";
  fgVariableNames[kIRIntClosestIntMap]   = "Closest out of bunch int. IRInt1";fgVariableUnits[kIRIntClosestIntMap]   = "";
  fgVariableNames[kIRIntClosestIntMap+1] = "Closest out of bunch int. IRInt2";fgVariableUnits[kIRIntClosestIntMap+1] = "";
  fgVariableNames[kNPMDtracks]           = "Number of PMD tracks";            fgVariableUnits[kNPMDtracks]           = "";
  fgVariableNames[kNTRDtracks]           = "Number of TRD tracks";            fgVariableUnits[kNTRDtracks]           = "";
  fgVariableNames[kNTRDtracklets]        = "Number of TRD tracklets";         fgVariableUnits[kNTRDtracklets]        = "";
  fgVariableNames[kNVtxContributors]     = "Number of vtx. contributors";     fgVariableUnits[kNVtxContributors]     = "";
  fgVariableNames[kNVtxTPCContributors]  = "Number of TPC vtx. contributors"; fgVariableUnits[kNVtxTPCContributors]  = "";
  fgVariableNames[kNVtxSPDContributors]  = "Number of SPD vtx. contributors"; fgVariableUnits[kNVtxSPDContributors]  = "";
  fgVariableNames[kVtxX]                 = "Vtx X";                           fgVariableUnits[kVtxX]                 = "cm";
  fgVariableNames[kVtxY]                 = "Vtx Y";                           fgVariableUnits[kVtxY]                 = "cm";
  fgVariableNames[kVtxZ]                 = "Vtx Z";                           fgVariableUnits[kVtxZ]                 = "cm";
  fgVariableNames[kVtxXtpc]              = "Vtx X TPC";                       fgVariableUnits[kVtxXtpc]              = "cm";
  fgVariableNames[kVtxYtpc]              = "Vtx Y TPC";                       fgVariableUnits[kVtxYtpc]              = "cm";
  fgVariableNames[kVtxZtpc]              = "Vtx Z TPC";                       fgVariableUnits[kVtxZtpc]              = "cm";
  fgVariableNames[kDeltaVtxZ]            = "#Delta Z";                        fgVariableUnits[kDeltaVtxZ]            = "cm";
  fgVariableNames[kVtxXspd]              = "Vtx X SPD";                       fgVariableUnits[kVtxXspd]              = "cm";
  fgVariableNames[kVtxYspd]              = "Vtx Y SPD";                       fgVariableUnits[kVtxYspd]              = "cm";
  fgVariableNames[kVtxZspd]              = "Vtx Z SPD";                       fgVariableUnits[kVtxZspd]              = "cm";
  fgVariableNames[kDeltaVtxZspd]         = "#Delta Z (global-SPD)";           fgVariableUnits[kDeltaVtxZspd]         = "cm";
  fgVariableNames[kTPCpileupZAC]         = "TPC pileup Z from A&C sides";     fgVariableUnits[kTPCpileupZAC]         = "cm";
  fgVariableNames[kTPCpileupZA]          = "TPC pileup Z from A side";        fgVariableUnits[kTPCpileupZA]          = "cm";
  fgVariableNames[kTPCpileupZC]          = "TPC pileup Z from C side";        fgVariableUnits[kTPCpileupZC]          = "cm";
  fgVariableNames[kTPCpileupContributorsAC] = "TPC pileup n-contributors from A&C sides"; fgVariableUnits[kTPCpileupContributorsAC] = "";
  fgVariableNames[kTPCpileupContributorsA]  = "TPC pileup n-contributors from A side"; fgVariableUnits[kTPCpileupContributorsA] = "";
  fgVariableNames[kTPCpileupContributorsC]  = "TPC pileup n-contributors from C side"; fgVariableUnits[kTPCpileupContributorsC] = "";
  fgVariableNames[kTPCpileupZAC2]         = "TPC pileup Z (2) from A&C sides";     fgVariableUnits[kTPCpileupZAC2]         = "cm";
  fgVariableNames[kTPCpileupZA2]          = "TPC pileup Z (2) from A side";        fgVariableUnits[kTPCpileupZA2]          = "cm";
  fgVariableNames[kTPCpileupZC2]          = "TPC pileup Z (2) from C side";        fgVariableUnits[kTPCpileupZC2]          = "cm";
  fgVariableNames[kTPCpileupContributorsAC2] = "TPC pileup n-contributors (2) from A&C sides"; fgVariableUnits[kTPCpileupContributorsAC2] = "";
  fgVariableNames[kTPCpileupContributorsA2]  = "TPC pileup n-contributors (2) from A side"; fgVariableUnits[kTPCpileupContributorsA2] = "";
  fgVariableNames[kTPCpileupContributorsC2]  = "TPC pileup n-contributors (2) from C side"; fgVariableUnits[kTPCpileupContributorsC2] = "";
  for(Int_t iflag=0; iflag<kNTrackingStatus; ++iflag) {
    fgVariableNames[kNTracksPerTrackingStatus+iflag] = Form("Tracks with %s on", fgkTrackingStatusNames[iflag]); 
    fgVariableUnits[kNTracksPerTrackingStatus+iflag] = ""; 
  }
  fgVariableNames[kNTracksTPCoutBeforeClean] = "Number of kTPCout before ESD cleanup"; fgVariableUnits[kNTracksTPCoutBeforeClean] = "";
  fgVariableNames[kNTracksTPCoutVsITSout]       = "TPCout/ITSout";                   fgVariableUnits[kNTracksTPCoutVsITSout] = "";
  fgVariableNames[kNTracksTRDoutVsITSout]       = "TRDout/ITSout";                   fgVariableUnits[kNTracksTRDoutVsITSout] = "";
  fgVariableNames[kNTracksTOFoutVsITSout]       = "TOFout/ITSout";                   fgVariableUnits[kNTracksTOFoutVsITSout] = "";
  fgVariableNames[kNTracksTRDoutVsTPCout]       = "TRDout/TPCout";                   fgVariableUnits[kNTracksTRDoutVsTPCout] = "";
  fgVariableNames[kNTracksTOFoutVsTPCout]       = "TOFout/TPCout";                   fgVariableUnits[kNTracksTOFoutVsTPCout] = "";
  fgVariableNames[kNTracksTOFoutVsTRDout]       = "TOFout/TRDout";                   fgVariableUnits[kNTracksTOFoutVsTRDout] = "";
  fgVariableNames[kNTracksITSoutVsSPDtracklets] = "ITSout/SPDtracklets";             fgVariableUnits[kNTracksITSoutVsSPDtracklets] = "";
  fgVariableNames[kNTracksTPCoutVsSPDtracklets] = "TPCout/SPDtracklets";             fgVariableUnits[kNTracksTPCoutVsSPDtracklets] = "";
  fgVariableNames[kNTracksTRDoutVsSPDtracklets] = "TRDout/SPDtracklets";             fgVariableUnits[kNTracksTRDoutVsSPDtracklets] = "";
  fgVariableNames[kNTracksTOFoutVsSPDtracklets] = "TOFout/SPDtracklets";             fgVariableUnits[kNTracksTOFoutVsSPDtracklets] = "";
  fgVariableNames[kNTracksTPCoutFromPileup] = "# kTPCout tracks - expectation";   fgVariableUnits[kNTracksTPCoutFromPileup] = "";
  fgVariableNames[kNTracksTPCoutVsVZEROTotalMult] = "TPCout / VZERO multiplicity"; fgVariableUnits[kNTracksTPCoutVsVZEROTotalMult] = "";
  fgVariableNames[kCentVZERO]                   = "VZERO centrality";                fgVariableUnits[kCentVZERO]      = "%";
  fgVariableNames[kCentSPD]                     = "CL1 centrality";                  fgVariableUnits[kCentSPD]        = "%";
  fgVariableNames[kCentSPDcorr]                 = "SPD trklts centrality";           fgVariableUnits[kCentSPDcorr]    = "%";
  fgVariableNames[kCentTPC]                     = "TPC centrality";                  fgVariableUnits[kCentTPC]        = "%";
  fgVariableNames[kCentZDC]                     = "ZDC centrality";                  fgVariableUnits[kCentZDC]        = "%";
  fgVariableNames[kCentVZEROA]                  = "VZERO-A centrality";              fgVariableUnits[kCentVZEROA]     = "%";
  fgVariableNames[kCentVZEROC]                  = "VZERO-C centrality";              fgVariableUnits[kCentVZEROC]     = "%";
  fgVariableNames[kCentZNA]                     = "ZNA centrality";                  fgVariableUnits[kCentZNA]        = "%";
  fgVariableNames[kCentQuality]                 = "Centrality quality";              fgVariableUnits[kCentQuality]    = "";
  fgVariableNames[kNV0total]                    = "Total number of V0s";             fgVariableUnits[kNV0total]       = "";  
  fgVariableNames[kNV0selected]                 = "Number of selected V0s";          fgVariableUnits[kNV0selected]    = "";  
  fgVariableNames[kNpairsSelected]              = "Number of pairs per event";       fgVariableUnits[kNpairsSelected] = "";    
  fgVariableNames[kEvAverageTPCchi2]        = "Event TPC <chi2>";            fgVariableUnits[kEvAverageTPCchi2] = "";
  fgVariableNames[kNDplusToK0sPiplusSelected]     = "Number of D+ ->K0s pi+ pairs per event";     fgVariableUnits[kNDplusToK0sPiplusSelected]     = "";
  fgVariableNames[kNDplusToK0sKplusSelected]      = "Number of D+ ->K0s K+ pairs per event";      fgVariableUnits[kNDplusToK0sKplusSelected]      = "";
  fgVariableNames[kNDplusToPhiPiplusSelected]     = "Number of D+ ->phi pi+ pairs per event";     fgVariableUnits[kNDplusToPhiPiplusSelected]     = "";
  fgVariableNames[kNDminusToK0sPiminusSelected]   = "Number of D- ->K0s pi- pairs per event";     fgVariableUnits[kNDminusToK0sPiminusSelected]   = "";
  fgVariableNames[kNDminusToK0sKminusSelected]    = "Number of D- ->K0s K- pairs per event";      fgVariableUnits[kNDminusToK0sKminusSelected]    = "";
  fgVariableNames[kNDminusToPhiPiminusSelected]   = "Number of D- ->phi pi- pairs per event";     fgVariableUnits[kNDminusToPhiPiminusSelected]   = "";  
  fgVariableNames[kNDzeroToKminusPiplusSelected]  = "Number of D0 ->K- pi+ pairs per event";      fgVariableUnits[kNDzeroToKminusPiplusSelected]  = "";    
  fgVariableNames[kNADzeroToKplusPiminusSelected] = "Number of anti-D0 ->K+ pi- pairs per event"; fgVariableUnits[kNADzeroToKplusPiminusSelected] = "";
  fgVariableNames[kNDsplusToK0sKplusSelected]     = "Number of Ds+ ->K0s K+ pairs per event";     fgVariableUnits[kNDsplusToK0sKplusSelected]     = "";
  fgVariableNames[kNDsminusToK0sKminusSelected]   = "Number of Ds- ->K0s K- pairs per event";     fgVariableUnits[kNDsminusToK0sKminusSelected]   = "";  
  fgVariableNames[kNtracksTotal]                = "No. of tracks in original event"; fgVariableUnits[kNtracksTotal]           = "";
  fgVariableNames[kNtracksSelected]             = "No. of selected tracks";          fgVariableUnits[kNtracksSelected]        = "";
  fgVariableNames[kNtracksPosAnalyzed]          = "No.selected positive tracks";     fgVariableUnits[kNtracksPosAnalyzed]     = "";  
  fgVariableNames[kNtracksNegAnalyzed]          = "No.selected negative tracks";     fgVariableUnits[kNtracksNegAnalyzed]     = ""; 
  fgVariableNames[kNtracksPiPlusAnalyzed]       = "No.selected pos. pions";          fgVariableUnits[kNtracksPiPlusAnalyzed]  = "";  
  fgVariableNames[kNtracksPiMinusAnalyzed]      = "No.selected neg. pions";          fgVariableUnits[kNtracksPiMinusAnalyzed] = "";  
  fgVariableNames[kNtracksKPlusAnalyzed]        = "No.selected pos. kaons";          fgVariableUnits[kNtracksKPlusAnalyzed]   = "";  
  fgVariableNames[kNtracksKMinusAnalyzed]       = "No.selected neg. kaons";          fgVariableUnits[kNtracksKMinusAnalyzed]  = "";  
  fgVariableNames[kNK0sAnalyzed]                = "No.selected K0s candidates";      fgVariableUnits[kNK0sAnalyzed]           = "";  
  fgVariableNames[kNPhiAnalyzed]                = "No.selected phi candidates";      fgVariableUnits[kNPhiAnalyzed]           = "";  
  fgVariableNames[kNtracksAnalyzed]             = "No.selected tracks";              fgVariableUnits[kNtracksAnalyzed]        = "";  
  for(Int_t i=0; i<18; ++i) fgVariableNames[kNtracksAnalyzedInPhiBins+i] = Form("# selected tracks in #varphi sector %d and #eta<0.", i);
  for(Int_t i=0; i<18; ++i) fgVariableNames[kNtracksAnalyzedInPhiBins+18+i] = Form("# selected tracks in #varphi sector %d and #eta>0.",i);
  fgVariableNames[kNtracksSubEvLeft]            = "No.tracks sub-event left";        fgVariableUnits[kNtracksSubEvLeft]       = "";  
  fgVariableNames[kNtracksSubEvRight]           = "No.tracks sub-event right";       fgVariableUnits[kNtracksSubEvRight]      = "";  
  fgVariableNames[kNtracksEventPlane]           = "No.tracks";                       fgVariableUnits[kNtracksEventPlane]      = "";  
  fgVariableNames[kNCaloClusters]               = "No.calorimeter clusters";         fgVariableUnits[kNCaloClusters]          = "";
  fgVariableNames[kNTPCclusters]                = "No. TPC clusters";                fgVariableUnits[kNTPCclusters]  = "";
  fgVariableNames[kNTPCclustersFromPileup]      = "No. TPC clusters above expectation w/o pileup"; fgVariableUnits[kNTPCclustersFromPileup] = "";
  fgVariableNames[kNTPCclustersFromPileupRelative] = "Relative no. TPC clusters above expectation w/o pileup"; fgVariableUnits[kNTPCclustersFromPileupRelative] = "";
  
  TString multEstimators[kNMultiplicityEstimators] = {
    "SPDntracklets10",
    "SPDntracklets08",
    "SPDntracklets16",
    "SPDntrackletsOuterEta",
    "SPDntracklets -1.6<#eta<-1.5",
    "SPDntracklets -1.5<#eta<-1.4",
    "SPDntracklets -1.4<#eta<-1.3",
    "SPDntracklets -1.3<#eta<-1.2",
    "SPDntracklets -1.2<#eta<-1.1",
    "SPDntracklets -1.1<#eta<-1.0",
    "SPDntracklets -1.0<#eta<-0.9",
    "SPDntracklets -0.9<#eta<-0.8",
    "SPDntracklets -0.8<#eta<-0.7",
    "SPDntracklets -0.7<#eta<-0.6",
    "SPDntracklets -0.6<#eta<-0.5",
    "SPDntracklets -0.5<#eta<-0.4",
    "SPDntracklets -0.4<#eta<-0.3",
    "SPDntracklets -0.3<#eta<-0.2",
    "SPDntracklets -0.2<#eta<-0.1",
    "SPDntracklets -0.1<#eta<0.0",
    "SPDntracklets 0.0<#eta<0.1",
    "SPDntracklets 0.1<#eta<0.2",
    "SPDntracklets 0.2<#eta<0.3",
    "SPDntracklets 0.3<#eta<0.4",
    "SPDntracklets 0.4<#eta<0.5",
    "SPDntracklets 0.5<#eta<0.6",
    "SPDntracklets 0.6<#eta<0.7",
    "SPDntracklets 0.7<#eta<0.8",
    "SPDntracklets 0.8<#eta<0.9",
    "SPDntracklets 0.9<#eta<1.0",
    "SPDntracklets 1.0<#eta<1.1",
    "SPDntracklets 1.1<#eta<1.2",
    "SPDntracklets 1.2<#eta<1.3",
    "SPDntracklets 1.3<#eta<1.4",
    "SPDntracklets 1.4<#eta<1.5",
    "SPDntracklets 1.5<#eta<1.6",
    "SPDntrackletsVtxEta",
    "VZEROTotalMult",
    "VZEROATotalMult",
    "VZEROCTotalMult",
    "VZEROAplusCTotalMult"
  };


  
  
  
  TString corrections[kNCorrections] = {
    ", Vtx. corr. (global)",
    ", Vtx. corr. (run-wise)",
    ", Vtx. corr. (global) + gain loss corr.",
    ", Vtx. corr. (run-wise) + gain loss corr.",
    ", 2D Vertex + gain loss corr.",
    ", gain loss corr."
  };

  TString referenceMultiplicities[kNReferenceMultiplicities] = {
    ", (max mult)",
    ", (min mult)",
    ", (mean mult)"
  };

  TString smearingMethods[kNSmearingMethods] = {
    "",
    ", Poisson smeared"
  };

  for( int iEstimator=0; iEstimator < kNMultiplicityEstimators; ++iEstimator){
    Int_t estimator = kMultiplicity + iEstimator;
    fgVariableNames[estimator] = multEstimators[iEstimator];
    fgVariableUnits[estimator]  = "";
    for( int iCorrection = 0; iCorrection < kNCorrections; ++iCorrection){
      for( int iReference=0; iReference < kNReferenceMultiplicities; ++iReference ){
        for( int iSmearing=0; iSmearing < kNSmearingMethods; ++iSmearing){
          Int_t index = GetCorrectedMultiplicity( estimator, iCorrection, iReference, iSmearing );
          fgVariableNames[index] = Form("%s%s%s%s",
                                        multEstimators[iEstimator].Data(),
                                        corrections[iCorrection].Data(),
                                        referenceMultiplicities[iReference].Data(),
                                        smearingMethods[iSmearing].Data()
                                   );
          fgVariableUnits[index] = "";
        }
      }
    }
  }

  fgVariableNames[kVZEROATotalMultFromChannels] = "VZERO-A total multiplicity"; fgVariableUnits[kVZEROATotalMultFromChannels] = "";
  fgVariableNames[kVZEROCTotalMultFromChannels] = "VZERO-C total multiplicity"; fgVariableUnits[kVZEROCTotalMultFromChannels] = "";
  fgVariableNames[kVZEROTotalMultFromChannels] = "VZERO total multiplicity"; fgVariableUnits[kVZEROTotalMultFromChannels] = "";

  for(Int_t il=0;il<2;++il) {
    fgVariableNames[kSPDFiredChips+il] = Form("Fired chips in SPD layer %d", il+1); 
    fgVariableUnits[kSPDFiredChips+il] = "";
  }
  for(Int_t il=0;il<6;++il) {
    fgVariableNames[kITSnClusters+il] = Form("No. clusters in ITS layer %d", il+1); 
    fgVariableUnits[kITSnClusters+il] = "";
  }
  fgVariableNames[kSPDnSingleClusters]  = "SPD single clusters";    fgVariableUnits[kSPDnSingleClusters]  = "";  
  fgVariableNames[kSDDandSSDclusters] = "SDD+SSD clusters";         fgVariableUnits[kSDDandSSDclusters] = "";
  fgVariableNames[kEventMixingId]       = "Event mixing id";        fgVariableUnits[kEventMixingId]       = "";  
  fgVariableNames[kVZEROCurrentChannel] = "VZERO channel";          fgVariableUnits[kVZEROCurrentChannel] = "";
  fgVariableNames[kVZEROCurrentChannelMult] = "VZERO channel multiplicity";   fgVariableUnits[kVZEROCurrentChannelMult] = "";
  fgVariableNames[kVZEROCurrentChannelMultCalib] = "VZERO channel calibrated multiplicity";   fgVariableUnits[kVZEROCurrentChannelMultCalib] = "";
  fgVariableNames[kVZEROAemptyChannels] = "VZERO-A empty channels"; fgVariableUnits[kVZEROAemptyChannels] = "";
  fgVariableNames[kVZEROCemptyChannels] = "VZERO-C empty channels"; fgVariableUnits[kVZEROCemptyChannels] = "";
  for(Int_t ich=0;ich<64;++ich) {
    fgVariableNames[kVZEROChannelMult+ich] = Form("Multiplicity VZERO ch.%d", ich);
    fgVariableUnits[kVZEROChannelMult+ich] = "";
    fgVariableNames[kVZEROChannelMultCalib+ich] = Form("Multiplicity VZERO ch.%d calibrated", ich);
    fgVariableUnits[kVZEROChannelMultCalib+ich] = "";
    fgVariableNames[kVZEROChannelEta+ich] = Form("#eta for VZERO ch.%d", ich);
    fgVariableUnits[kVZEROChannelEta+ich] = "";
    fgVariableNames[kVZEROflowV2TPC+ich] = Form("v_{2}^{VZERO-ch%d}{EP,TPC}", ich);
    fgVariableUnits[kVZEROflowV2TPC+ich] = "";
  }
  for(Int_t ich=0;ich<10;++ich) {
    fgVariableNames[kZDCnEnergyCh+ich] = Form("Energy ZDCn ch.%d", ich);
    fgVariableUnits[kZDCnEnergyCh+ich] = "GeV";
    fgVariableNames[kZDCpEnergyCh+ich] = Form("Energy ZDCp ch.%d", ich);
    fgVariableUnits[kZDCpEnergyCh+ich] = "GeV";
  }
  for(Int_t ich=0;ich<26;++ich) {
    fgVariableNames[kTZEROAmplitudeCh+ich] = Form("TZERO amplitude ch.%d",ich);
    fgVariableUnits[kTZEROAmplitudeCh+ich] = "";
  }
  fgVariableNames[kTZEROTOF]       = "TOF TZERO A&C first";
  fgVariableNames[kTZEROTOF+1]     = "TOF TZERO A first";
  fgVariableNames[kTZEROTOF+2]     = "TOF TZERO C first";
  fgVariableNames[kTZEROTOFbest]   = "TOF TZERO A&C best";
  fgVariableNames[kTZEROTOFbest+1] = "TOF TZERO A best";
  fgVariableNames[kTZEROTOFbest+2] = "TOF TZERO C best";
  for(Int_t i=0; i<3; ++i) {fgVariableUnits[kTZEROTOF+i] = "ps"; fgVariableUnits[kTZEROTOFbest+i] = "ps";}
  fgVariableNames[kTZEROzVtx]      = "TZERO Vtx. Z";                 fgVariableUnits[kTZEROzVtx]      = "cm.";
  fgVariableNames[kTZEROstartTime] = "TZERO start time";             fgVariableUnits[kTZEROstartTime] = "ps";
  fgVariableNames[kTZEROpileup]    = "TZERO pileup ON";              fgVariableUnits[kTZEROpileup]    = "";
  fgVariableNames[kTZEROsatellite] = "TZERO satellite collision BC"; fgVariableUnits[kTZEROsatellite] = "";

  fgVariableNames[kMultEstimatorV0M]    = "Muliplicity estimator V0M";
  fgVariableNames[kMultEstimatorV0A]    = "Muliplicity estimator V0A";
  fgVariableNames[kMultEstimatorV0C]    = "Muliplicity estimator V0C";
  fgVariableNames[kMultEstimatorOnlineV0M]    = "Muliplicity estimator OnlineV0M";
  fgVariableNames[kMultEstimatorOnlineV0A]    = "Muliplicity estimator OnlineV0A";
  fgVariableNames[kMultEstimatorOnlineV0C]    = "Muliplicity estimator OnlineV0C";
  fgVariableNames[kMultEstimatorADM]          = "Muliplicity estimator ADM ";
  fgVariableNames[kMultEstimatorADA]          = "Muliplicity estimator ADA";
  fgVariableNames[kMultEstimatorADC]          = "Muliplicity estimator ADC";
  fgVariableNames[kMultEstimatorSPDClusters]  = "Muliplicity estimator SPDclusters";
  fgVariableNames[kMultEstimatorSPDTracklets] = "Muliplicity estimator SPDtracklets";
  fgVariableNames[kMultEstimatorRefMult05]    = "Muliplicity estimator RefMult05";
  fgVariableNames[kMultEstimatorRefMult08]    = "Muliplicity estimator RefMult08";

  fgVariableUnits[kMultEstimatorV0M]    = "";
  fgVariableUnits[kMultEstimatorV0A]    = "";
  fgVariableUnits[kMultEstimatorV0C]    = "";
  fgVariableUnits[kMultEstimatorOnlineV0M]    = "";
  fgVariableUnits[kMultEstimatorOnlineV0A]    = "";
  fgVariableUnits[kMultEstimatorOnlineV0C]    = "";
  fgVariableUnits[kMultEstimatorADM]          = "";
  fgVariableUnits[kMultEstimatorADA]          = "";
  fgVariableUnits[kMultEstimatorADC]          = "";
  fgVariableUnits[kMultEstimatorSPDClusters]  = "";
  fgVariableUnits[kMultEstimatorSPDTracklets] = "";
  fgVariableUnits[kMultEstimatorRefMult05]    = "";
  fgVariableUnits[kMultEstimatorRefMult08]    = "";


  fgVariableNames[kMultEstimatorPercentileV0M]    = "Muliplicity estimator percentile V0M";
  fgVariableNames[kMultEstimatorPercentileV0A]    = "Muliplicity estimator percentile V0A";
  fgVariableNames[kMultEstimatorPercentileV0C]    = "Muliplicity estimator percentile V0C";
  fgVariableNames[kMultEstimatorPercentileOnlineV0M]    = "Muliplicity estimator percentile OnlineV0M";
  fgVariableNames[kMultEstimatorPercentileOnlineV0A]    = "Muliplicity estimator percentile OnlineV0A";
  fgVariableNames[kMultEstimatorPercentileOnlineV0C]    = "Muliplicity estimator percentile OnlineV0C";
  fgVariableNames[kMultEstimatorPercentileADM]          = "Muliplicity estimator percentile ADM ";
  fgVariableNames[kMultEstimatorPercentileADA]          = "Muliplicity estimator percentile ADA";
  fgVariableNames[kMultEstimatorPercentileADC]          = "Muliplicity estimator percentile ADC";
  fgVariableNames[kMultEstimatorPercentileSPDClusters]  = "Muliplicity estimator percentile SPDclusters";
  fgVariableNames[kMultEstimatorPercentileSPDTracklets] = "Muliplicity estimator percentile SPDtracklets";
  fgVariableNames[kMultEstimatorPercentileRefMult05]    = "Muliplicity estimator percentile RefMult05";
  fgVariableNames[kMultEstimatorPercentileRefMult08]    = "Muliplicity estimator percentile RefMult08";

  fgVariableUnits[kMultEstimatorPercentileV0M]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileV0A]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileV0C]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileOnlineV0M]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileOnlineV0A]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileOnlineV0C]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileADM]          = "\%";
  fgVariableUnits[kMultEstimatorPercentileADA]          = "\%";
  fgVariableUnits[kMultEstimatorPercentileADC]          = "\%";
  fgVariableUnits[kMultEstimatorPercentileSPDClusters]  = "\%";
  fgVariableUnits[kMultEstimatorPercentileSPDTracklets] = "\%";
  fgVariableUnits[kMultEstimatorPercentileRefMult05]    = "\%";
  fgVariableUnits[kMultEstimatorPercentileRefMult08]    = "\%";


  fgVariableNames[kINT7Triggered]       = "event was triggered with INT7";       fgVariableUnits[kINT7Triggered]       = "";
  fgVariableNames[kTRDTriggeredType]    = "event was triggered by TRD ele trigger"; fgVariableUnits[kTRDTriggeredType]       = "";
  fgVariableNames[kHighMultV0Triggered] = "event was triggered with HighMultV0"; fgVariableUnits[kHighMultV0Triggered] = "";
  fgVariableNames[kEMCEGATriggered]     = "event was triggered with EMCEGA";              fgVariableUnits[kEMCEGATriggered] = "";
  fgVariableNames[kEMCEGAHighTriggered] = "event was triggered with EMCEGA high thresh."; fgVariableUnits[kEMCEGAHighTriggered] = "";

  TString vzeroSideNames[3] = {"A","C","AC"};
  for(Int_t iHarmonic=0;iHarmonic<6;++iHarmonic) {
    for(Int_t iSide=0;iSide<3;++iSide) {
      fgVariableNames[kVZEROQvecX+iSide*6+iHarmonic]    = Form("Q_{x,%d}^{VZERO-%s}",iHarmonic+1,vzeroSideNames[iSide].Data());
      fgVariableUnits[kVZEROQvecX+iSide*6+iHarmonic]    = "";
      fgVariableNames[kVZEROQvecY+iSide*6+iHarmonic]    = Form("Q_{y,%d}^{VZERO-%s}",iHarmonic+1,vzeroSideNames[iSide].Data());
      fgVariableUnits[kVZEROQvecY+iSide*6+iHarmonic]    = "";
      fgVariableNames[kVZEROQvecMag+iSide*6+iHarmonic]  = Form("|Q_{%d}^{VZERO-%s}|",iHarmonic+1,vzeroSideNames[iSide].Data());
      fgVariableUnits[kVZEROQvecMag+iSide*6+iHarmonic]  = "";
      fgVariableNames[kVZERORP+iSide*6+iHarmonic]       = Form("#Psi_{%d}^{VZERO-%s}",iHarmonic+1,vzeroSideNames[iSide].Data());
      fgVariableUnits[kVZERORP+iSide*6+iHarmonic]       = "rad.";
      fgVariableNames[kVZEROFlowVn+iSide*6+iHarmonic]   = Form("v_{%d}{EP,VZERO-%s}",iHarmonic+1,vzeroSideNames[iSide].Data());
      fgVariableUnits[kVZEROFlowVn+iSide*6+iHarmonic]   = "";
      fgVariableNames[kVZEROFlowSine+iSide*6+iHarmonic] = Form("sin(%d(#varphi-#Psi_{%d}^{VZERO-%s}))",
								   iHarmonic+1,iHarmonic+1,vzeroSideNames[iSide].Data());
      fgVariableUnits[kVZEROFlowSine+iSide*6+iHarmonic] = "";
      if(iSide<2) {
	fgVariableNames[kVZEROuQ+iSide*6+iHarmonic] = Form("cos(%d(#phi - #Psi_{%d}^{VZERO-%s})|Q^{VZERO-%s}_{%d}|",
							       iHarmonic+1,iHarmonic+1,vzeroSideNames[iSide].Data(),vzeroSideNames[iSide].Data(),iHarmonic+1);
        fgVariableUnits[kVZEROuQ+iSide*6+iHarmonic] = "";
	fgVariableNames[kVZEROuQsine+iSide*6+iHarmonic] = Form("sin(%d(#phi - #Psi_{%d}^{VZERO-%s})|Q^{VZERO-%s}_{%d}|",
							           iHarmonic+1,iHarmonic+1,vzeroSideNames[iSide].Data(),vzeroSideNames[iSide].Data(),iHarmonic+1);
        fgVariableUnits[kVZEROuQsine+iSide*6+iHarmonic] = "";
      }
    }
    fgVariableNames[kVZERORPres+iHarmonic] = Form("#sqrt{%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{VZERO-C})}",
						    iHarmonic+1,iHarmonic+1,iHarmonic+1);
    fgVariableUnits[kVZERORPres+iHarmonic] = "rad.^{1/2}";
    fgVariableNames[kVZEROARPres+iHarmonic] = Form("#sqrt{cos(%d(#Psi_{%d}^{VZERO-C}-#Psi_{%d}^{TPC}))/cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{VZERO-C}))*cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{TPC}))}",
						    iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1);
    fgVariableUnits[kVZEROARPres+iHarmonic] = "rad.^{1/2}";
    fgVariableNames[kVZEROCRPres+iHarmonic] = Form("#sqrt{cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{TPC}))/cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{VZERO-C}))*cos(%d(#Psi_{%d}^{VZERO-C}-#Psi_{%d}^{TPC}))}",
						    iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1);
    fgVariableUnits[kVZEROCRPres+iHarmonic] = "rad.^{1/2}";
    fgVariableNames[kVZEROXaXc+iHarmonic] = Form("Q_{x,%d}^{VZERO-A} #times Q_{x,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROXaXc+iHarmonic] = "";
    fgVariableNames[kVZEROXaYa+iHarmonic] = Form("Q_{x,%d}^{VZERO-A} #times Q_{y,%d}^{VZERO-A}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROXaYa+iHarmonic] = "";
    fgVariableNames[kVZEROXaYc+iHarmonic] = Form("Q_{x,%d}^{VZERO-A} #times Q_{y,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROXaYc+iHarmonic] = "";
    fgVariableNames[kVZEROYaXc+iHarmonic] = Form("Q_{y,%d}^{VZERO-A} #times Q_{x,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROYaXc+iHarmonic] = "";
    fgVariableNames[kVZEROYaYc+iHarmonic] = Form("Q_{y,%d}^{VZERO-A} #times Q_{y,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROYaYc+iHarmonic] = "";
    fgVariableNames[kVZEROXcYc+iHarmonic] = Form("Q_{x,%d}^{VZERO-C} #times Q_{y,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROXcYc+iHarmonic] = "";
    fgVariableNames[kVZEROdeltaRPac+iHarmonic]  = Form("#Psi_{%d}^{VZERO-A} - #Psi_{%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROdeltaRPac+iHarmonic]  = "rad.";
    fgVariableNames[kVZEROQaQcSP+iHarmonic] = Form("cos(%d(#Psi_{%d}^{VZERO-A} - #Psi_{%d}^{VZERO-C})|Q^{VZERO-A}_{%d}||Q^{VZERO-C}_{%d}|", 
						       iHarmonic+1, iHarmonic+1, iHarmonic+1, iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROQaQcSP+iHarmonic] = "";
    fgVariableNames[kVZEROQaQcSPsine+iHarmonic] = Form("sin(%d(#Psi_{%d}^{VZERO-A} - #Psi_{%d}^{VZERO-C})|Q^{VZERO-A}_{%d}||Q^{VZERO-C}_{%d}|", 
							   iHarmonic+1, iHarmonic+1, iHarmonic+1, iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kVZEROQaQcSPsine+iHarmonic] = "";
    fgVariableUnits[kVZEROdeltaRPac+iHarmonic]  = "rad.";
    fgVariableNames[kTPCQvecX+iHarmonic] = Form("Q_{x,%d}^{TPC}", iHarmonic+1);
    fgVariableUnits[kTPCQvecX+iHarmonic] = "";
    fgVariableNames[kTPCQvecY+iHarmonic] = Form("Q_{y,%d}^{TPC}", iHarmonic+1);
    fgVariableUnits[kTPCQvecY+iHarmonic] = "";
    fgVariableNames[kTPCRP+iHarmonic]    = Form("#Psi_{%d}^{TPC}", iHarmonic+1);
    fgVariableUnits[kTPCRP+iHarmonic]    = "rad.";
    fgVariableNames[kTPCRPres+iHarmonic] = Form("#sqrt{%d(#Psi_{%d}^{TPC}-#Psi_{%d}^{VZERO-A})}", 
						  iHarmonic+1, iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kTPCRPres+iHarmonic] = "rad^{1/2}";
    fgVariableNames[kVZEROTPCRPres+6+iHarmonic] = Form("#sqrt{cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{VZERO-C}))/cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{TPC}))*cos(%d(#Psi_{%d}^{VZERO-A}-#Psi_{%d}^{TPC}))}",
						    iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1,iHarmonic+1);
    fgVariableUnits[kVZEROTPCRPres+6+iHarmonic] = "rad^{1/2}";
    fgVariableNames[kRPXtpcXvzeroa+iHarmonic] = Form("Q_{x,%d}^{TPC} #times Q_{x,%d}^{VZERO-A}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPXtpcXvzeroa+iHarmonic] = "";
    fgVariableNames[kRPXtpcXvzeroc+iHarmonic] = Form("Q_{x,%d}^{TPC} #times Q_{x,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPXtpcXvzeroc+iHarmonic] = "";
    fgVariableNames[kRPYtpcYvzeroa+iHarmonic] = Form("Q_{y,%d}^{TPC} #times Q_{y,%d}^{VZERO-A}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPYtpcYvzeroa+iHarmonic] = "";
    fgVariableNames[kRPYtpcYvzeroc+iHarmonic] = Form("Q_{y,%d}^{TPC} #times Q_{y,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPYtpcYvzeroc+iHarmonic] = "";
    fgVariableNames[kRPXtpcYvzeroa+iHarmonic] = Form("Q_{x,%d}^{TPC} #times Q_{y,%d}^{VZERO-A}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPXtpcYvzeroa+iHarmonic] = "";
    fgVariableNames[kRPXtpcYvzeroc+iHarmonic] = Form("Q_{x,%d}^{TPC} #times Q_{y,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPXtpcYvzeroc+iHarmonic] = "";
    fgVariableNames[kRPYtpcXvzeroa+iHarmonic] = Form("Q_{y,%d}^{TPC} #times Q_{x,%d}^{VZERO-A}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPYtpcXvzeroa+iHarmonic] = "";
    fgVariableNames[kRPYtpcXvzeroc+iHarmonic] = Form("Q_{y,%d}^{TPC} #times Q_{x,%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPYtpcXvzeroc+iHarmonic] = "";
    fgVariableNames[kRPdeltaVZEROAtpc+iHarmonic] = Form("#Psi_{%d}^{TPC} - #Psi_{%d}^{VZERO-A}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPdeltaVZEROAtpc+iHarmonic] = "rad.";
    fgVariableNames[kRPdeltaVZEROCtpc+iHarmonic] = Form("#Psi_{%d}^{TPC} - #Psi_{%d}^{VZERO-C}", iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kRPdeltaVZEROCtpc+iHarmonic] = "rad.";
    fgVariableNames[kTPCQvecXleft+iHarmonic] = Form("Q_{x,%d}^{TPC,left}", iHarmonic+1);
    fgVariableUnits[kTPCQvecXleft+iHarmonic] = "";
    fgVariableNames[kTPCQvecYleft+iHarmonic] = Form("Q_{y,%d}^{TPC,left}", iHarmonic+1);
    fgVariableUnits[kTPCQvecYleft+iHarmonic] = "";
    fgVariableNames[kTPCRPleft+iHarmonic] = Form("#Psi_{%d}^{TPC,left}", iHarmonic+1);
    fgVariableUnits[kTPCRPleft+iHarmonic] = "rad.";
    fgVariableNames[kTPCQvecXright+iHarmonic] = Form("Q_{x,%d}^{TPC,right}", iHarmonic+1);
    fgVariableUnits[kTPCQvecXright+iHarmonic] = "";
    fgVariableNames[kTPCQvecYright+iHarmonic] = Form("Q_{y,%d}^{TPC,right}", iHarmonic+1);
    fgVariableUnits[kTPCQvecYright+iHarmonic] = "";
    fgVariableNames[kTPCRPright+iHarmonic] = Form("#Psi_{%d}^{TPC,right}", iHarmonic+1);
    fgVariableUnits[kTPCRPright+iHarmonic] = "rad.";
    fgVariableNames[kTPCQvecXtotal+iHarmonic] = Form("Q_{x,%d}^{TPC,total}", iHarmonic+1);
    fgVariableUnits[kTPCQvecXtotal+iHarmonic] = "";
    fgVariableNames[kTPCQvecYtotal+iHarmonic] = Form("Q_{y,%d}^{TPC,total}", iHarmonic+1);
    fgVariableUnits[kTPCQvecYtotal+iHarmonic] = "";
    fgVariableNames[kTPCRPtotal+iHarmonic] = Form("#Psi_{%d}^{TPC,total}", iHarmonic+1);
    fgVariableUnits[kTPCRPtotal+iHarmonic] = "rad.";
    fgVariableNames[kTPCsubResCos+iHarmonic] = Form("cos(%d(#Psi_{%d}^{TPC,left}-#Psi_{%d}^{TPC,right}))", 
						      iHarmonic+1, iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kTPCsubResCos+iHarmonic] = "";
    fgVariableNames[kTPCFlowVn+iHarmonic] = Form("v_{%d}{EP,TPC}", iHarmonic+1); 
    fgVariableUnits[kTPCFlowVn] = "";
    fgVariableNames[kTPCFlowSine+iHarmonic] = Form("sin(%d(#varphi-#Psi_{%d}^{TPC}))", iHarmonic+1, iHarmonic+1); 
    fgVariableUnits[kTPCFlowSine] = "";
    fgVariableNames[kTPCuQ+iHarmonic] = Form("cos(%d(#phi - #Psi_{%d}^{TPC})|Q^{TPC}_{%d}|", iHarmonic+1,iHarmonic+1,iHarmonic+1);
    fgVariableUnits[kTPCuQ+iHarmonic] = "";
    fgVariableNames[kTPCuQsine+iHarmonic] = Form("sin(%d(#phi - #Psi_{%d}^{TPC})|Q^{TPC}_{%d}|", iHarmonic+1,iHarmonic+1,iHarmonic+1);
    fgVariableUnits[kTPCuQsine+iHarmonic] = "";
    fgVariableNames[kTPCQvecXtree+iHarmonic] = Form("Q_{x,%d}^{TPC}", iHarmonic+1);
    fgVariableUnits[kTPCQvecXtree+iHarmonic] = "";
    fgVariableNames[kTPCQvecYtree+iHarmonic] = Form("Q_{y,%d}^{TPC}", iHarmonic+1);
    fgVariableUnits[kTPCQvecYtree+iHarmonic] = "";
    fgVariableNames[kTPCRPtree+iHarmonic]    = Form("#Psi_{%d}^{TPC}", iHarmonic+1);
    fgVariableUnits[kTPCRPtree+iHarmonic]    = "rad.";
    fgVariableNames[kTPCQvecXptWeightsTree+iHarmonic] = Form("Q_{x,%d}^{TPC}, pt weights", iHarmonic+1);
    fgVariableUnits[kTPCQvecXptWeightsTree+iHarmonic] = "";
    fgVariableNames[kTPCQvecYptWeightsTree+iHarmonic] = Form("Q_{y,%d}^{TPC}, pt weights", iHarmonic+1);
    fgVariableUnits[kTPCQvecYptWeightsTree+iHarmonic] = "";
    fgVariableNames[kTPCRPptWeightsTree+iHarmonic]    = Form("#Psi_{%d}^{TPC}, pt weights", iHarmonic+1);
    fgVariableUnits[kTPCRPptWeightsTree+iHarmonic]    = "rad.";
    fgVariableNames[kTPCQvecXposTree+iHarmonic] = Form("Q_{x,%d}^{TPC}, positives", iHarmonic+1);
    fgVariableUnits[kTPCQvecXposTree+iHarmonic] = "";
    fgVariableNames[kTPCQvecYposTree+iHarmonic] = Form("Q_{y,%d}^{TPC}, positives", iHarmonic+1);
    fgVariableUnits[kTPCQvecYposTree+iHarmonic] = "";
    fgVariableNames[kTPCRPposTree+iHarmonic]    = Form("#Psi_{%d}^{TPC}, positives", iHarmonic+1);
    fgVariableUnits[kTPCRPposTree+iHarmonic]    = "rad.";
    fgVariableNames[kTPCQvecXnegTree+iHarmonic] = Form("Q_{x,%d}^{TPC}, negatives", iHarmonic+1);
    fgVariableUnits[kTPCQvecXnegTree+iHarmonic] = "";
    fgVariableNames[kTPCQvecYnegTree+iHarmonic] = Form("Q_{y,%d}^{TPC}, negatives", iHarmonic+1);
    fgVariableUnits[kTPCQvecYnegTree+iHarmonic] = "";
    fgVariableNames[kTPCRPnegTree+iHarmonic]    = Form("#Psi_{%d}^{TPC}, negatives", iHarmonic+1);
    fgVariableUnits[kTPCRPnegTree+iHarmonic]    = "rad.";
  }  // end loop over harmonics 
  
  fgVariableNames[kMCNch] = "N_{ch} in |#eta|<1"; fgVariableUnits[kMCNch] = ""; 
  fgVariableNames[kMCNchSPDacc] = "N_{ch} in SPD acceptance"; fgVariableUnits[kMCNchSPDacc] = ""; 
  fgVariableNames[kEtaBinForSPDtracklets] = "#eta"; fgVariableUnits[kEtaBinForSPDtracklets] = "";
  fgVariableNames[kMCNchNegSide] = "N_{ch} in -1<#eta<0"; fgVariableUnits[kMCNchNegSide] = "";
  fgVariableNames[kMCNchPosSide] = "N_{ch} in 0<#eta<1"; fgVariableUnits[kMCNchPosSide] = ""; 
  fgVariableNames[kDiffNchSPDtrklts] = "N_{ch}^{|#eta|<1} - SPD #trklts"; fgVariableUnits[kDiffNchSPDtrklts] = "";
  fgVariableNames[kDiffNchSPDaccSPDtrklts] = "N_{ch}^{SPD} - SPD #trklts"; fgVariableUnits[kDiffNchSPDaccSPDtrklts] = "";
  fgVariableNames[kRelDiffNchSPDtrklts] = "(N_{ch}^{|#eta|<1} - SPD #trklts) / N_{ch}^{|#eta|<1}"; fgVariableUnits[kRelDiffNchSPDtrklts] = "";
  fgVariableNames[kRelDiffNchSPDaccSPDtrklts] = "(N_{ch}^{SPD} - SPD #trklts) / N_{ch}^{SPD}"; fgVariableUnits[kRelDiffNchSPDaccSPDtrklts] = "";
  fgVariableNames[kRelDiff2NchSPDtrklts] = "2*(N_{ch}^{|#eta|<1} - SPD #trklts) / (N_{ch}^{|#eta|<1} + SPD #trklts)"; fgVariableUnits[kRelDiff2NchSPDtrklts] = "";
  fgVariableNames[kRelDiff2NchSPDaccSPDtrklts] = "2*(N_{ch}^{SPD} - SPD #trklts) / (N_{ch}^{SPD} + SPD #trklts)"; fgVariableUnits[kRelDiff2NchSPDaccSPDtrklts] = "";
  
  fgVariableNames[kPt]    = "p_{T}";   fgVariableUnits[kPt]    = "GeV/c";
  fgVariableNames[kPtMC]    = "p^{MC}_{T}";   fgVariableUnits[kPtMC]    = "GeV/c";
  fgVariableNames[kPtMCfromLegs]    = "p^{MC-legs}_{T}";   fgVariableUnits[kPtMCfromLegs]    = "GeV/c";
  fgVariableNames[kP]     = "p";       fgVariableUnits[kP]     = "GeV/c";
  fgVariableNames[kPMC]     = "p^{MC}";       fgVariableUnits[kPMC]     = "GeV/c";
  fgVariableNames[kPMCfromLegs]     = "p^{MC-legs}";       fgVariableUnits[kPMCfromLegs]     = "GeV/c";
  fgVariableNames[kPx]    = "p_{x}";   fgVariableUnits[kPx]    = "GeV/c";
  fgVariableNames[kPxMC]    = "p^{MC}_{x}";   fgVariableUnits[kPxMC]    = "GeV/c";
  fgVariableNames[kPxMCfromLegs]    = "p^{MC-legs}_{x}";   fgVariableUnits[kPxMCfromLegs]    = "GeV/c";
  fgVariableNames[kPy]    = "p_{y}";   fgVariableUnits[kPy]    = "GeV/c";
  fgVariableNames[kPyMC]    = "p^{MC}_{y}";   fgVariableUnits[kPyMC]    = "GeV/c";
  fgVariableNames[kPyMCfromLegs]    = "p^{MC-legs}_{y}";   fgVariableUnits[kPyMCfromLegs]    = "GeV/c";
  fgVariableNames[kPz]    = "p_{z}";   fgVariableUnits[kPz]    = "GeV/c";
  fgVariableNames[kPzMC]    = "p^{MC}_{z}";   fgVariableUnits[kPzMC]    = "GeV/c";
  fgVariableNames[kPzMCfromLegs]    = "p^{MC-legs}_{z}";   fgVariableUnits[kPzMCfromLegs]    = "GeV/c";
  fgVariableNames[kTheta] = "#theta";  fgVariableUnits[kTheta] = "rad.";
  fgVariableNames[kThetaMC] = "#theta^{MC}";  fgVariableUnits[kThetaMC] = "rad.";
  fgVariableNames[kThetaMCfromLegs] = "#theta^{MC-legs}";  fgVariableUnits[kThetaMCfromLegs] = "rad.";
  fgVariableNames[kEta]   = "#eta";    fgVariableUnits[kEta]   = "";
  fgVariableNames[kEtaMC]   = "#eta^{MC}";    fgVariableUnits[kEtaMC]   = "";
  fgVariableNames[kEtaMCfromLegs]   = "#eta^{MC-legs}";    fgVariableUnits[kEtaMCfromLegs]   = "";
  fgVariableNames[kPhi]   = "#varphi"; fgVariableUnits[kPhi]   = "rad.";
  fgVariableNames[kPhiMC]   = "#varphi^{MC}"; fgVariableUnits[kPhiMC]   = "rad.";
  fgVariableNames[kPhiMCfromLegs]   = "#varphi^{MC-legs}"; fgVariableUnits[kPhiMCfromLegs]   = "rad.";
  for(Int_t iHarmonic=0;iHarmonic<6;++iHarmonic) {
    fgVariableNames[kCosNPhi+iHarmonic] = Form("cos(%d#varphi)",iHarmonic+1); fgVariableUnits[kCosNPhi+iHarmonic] = "";
    fgVariableNames[kSinNPhi+iHarmonic] = Form("sin(%d#varphi)",iHarmonic+1); fgVariableUnits[kSinNPhi+iHarmonic] = "";
  }
  fgVariableNames[kPtSquared]     = "p_{T}^{2}";     fgVariableUnits[kPtSquared]     = "GeV^{2}/c^{2}";
  fgVariableNames[kOneOverSqrtPt] = "1/sqrt(p_{T})"; fgVariableUnits[kOneOverSqrtPt] = "GeV^{-1/2}";
  fgVariableNames[kMass]      = "m";         fgVariableUnits[kMass]        = "GeV/c^{2}";
  fgVariableNames[kMassMC]      = "m^{MC}";         fgVariableUnits[kMassMC]        = "GeV/c^{2}";
  fgVariableNames[kMassMCfromLegs]      = "m^{MC-legs}";         fgVariableUnits[kMassMCfromLegs]        = "GeV/c^{2}";
  fgVariableNames[kRap]       = "y";          fgVariableUnits[kRap]         = "";
  fgVariableNames[kRapAbs]    = "|y|";        fgVariableUnits[kRapAbs]      = "";
  fgVariableNames[kRapMC]     = "y^{MC}";     fgVariableUnits[kRapMC]       = "";
  fgVariableNames[kRapMCAbs]  = "|y^{MC}|";   fgVariableUnits[kRapMCAbs]    = "";
  fgVariableNames[kRapMCfromLegs]       = "y^{MC-legs}";         fgVariableUnits[kRapMCfromLegs]         = "";
  fgVariableNames[kPdgMC]       = "PDG code";          fgVariableUnits[kPdgMC]         = "";
  fgVariableNames[kPdgMC+1]  = "mother's PDG code";     fgVariableUnits[kPdgMC+1]         = "";
  fgVariableNames[kPdgMC+2]  = "grand-mother's PDG code";     fgVariableUnits[kPdgMC+2]         = "";
  fgVariableNames[kPdgMC+3]  = "grand-grand-mother's PDG code";     fgVariableUnits[kPdgMC+3]         = "";
  
  fgVariableNames[kCandidateId]       = "pair id.";              fgVariableUnits[kCandidateId]       = "";
  fgVariableNames[kPairType]          = "pair type";             fgVariableUnits[kPairType]          = "";
  fgVariableNames[kPairTypeSPD]          = "pair type spd legs";             fgVariableUnits[kPairTypeSPD]          = "";
  fgVariableNames[kMassV0]            = "m_{K^{0}_{S}}";         fgVariableUnits[kMassV0]            = "GeV/c^{2}";
  fgVariableNames[kMassV0+1]          = "m_{#Lambda^{0}}";       fgVariableUnits[kMassV0+1]          = "GeV/c^{2}";
  fgVariableNames[kMassV0+2]          = "m_{#bar{#Lambda^{0}}}"; fgVariableUnits[kMassV0+2]          = "GeV/c^{2}";
  fgVariableNames[kMassV0+3]          = "m_{#gamma}";            fgVariableUnits[kMassV0+3]          = "GeV/c^{2}";
  fgVariableNames[kPairChisquare]     = "pair #chi^{2}";         fgVariableUnits[kPairChisquare]     = "";
  fgVariableNames[kPairLxy]           = "L_{xy}";                fgVariableUnits[kPairLxy]           = "cm.";
  fgVariableNames[kPseudoProperDecayTime]  = "t";                fgVariableUnits[kPseudoProperDecayTime]  = "cm./c";
  fgVariableNames[kPseudoProperDecayTimeMC]  = "t_{MC}";              fgVariableUnits[kPseudoProperDecayTimeMC]  = "cm./c";
  fgVariableNames[kPairOpeningAngle]  = "pair opening angle";    fgVariableUnits[kPairOpeningAngle]  = "rad.";    
  fgVariableNames[kPairPointingAngle] = "#theta_{pointing}";     fgVariableUnits[kPairPointingAngle] = "rad.";
  fgVariableNames[kPairThetaCS]       = "cos(#theta^{*}_{CS})";  fgVariableUnits[kPairThetaCS]       = "";  
  fgVariableNames[kPairPhiCS]         = "#varphi^{*}_{CS}";      fgVariableUnits[kPairPhiCS]         = "rad.";  
  fgVariableNames[kPairThetaHE]       = "cos(#theta^{*}_{HE})";  fgVariableUnits[kPairThetaHE]       = "";  
  fgVariableNames[kPairPhiHE]         = "#varphi^{*}_{HE}";      fgVariableUnits[kPairPhiHE]         = "rad.";
  fgVariableNames[kPairPhiV]          = "#varphi^{*}_{v}";       fgVariableUnits[kPairPhiV]          = "rad.";
  fgVariableNames[kPairEff]           = "pair eff.";             fgVariableUnits[kPairEff]           = "";
  fgVariableNames[kOneOverPairEff]    = "1/pair eff.";           fgVariableUnits[kOneOverPairEff]    = "";
  fgVariableNames[kOneOverPairEffSq]  = "1/pair eff. squared";   fgVariableUnits[kOneOverPairEffSq]  = "";
  for(Int_t i=0;i<2;++i) {
     fgVariableNames[kPairLegTPCchi2+i] = Form("TPC #chi^{2}, leg %d", i+1);
     fgVariableUnits[kPairLegTPCchi2+i] = "";
     fgVariableNames[kPairLegITSchi2+i] = Form("ITS #chi^{2}, leg %d", i+1);
     fgVariableUnits[kPairLegITSchi2+i] = "";
     fgVariableNames[kPairLegPt+i] = Form("Leg%d p_{T}", i+1);
     fgVariableNames[kPairLegPtMC+i] = Form("Leg%d p^{MC}_{T}", i+1);
     fgVariableUnits[kPairLegPt+i] = "GeV/c"; fgVariableUnits[kPairLegPtMC+i] = "GeV/c";
  }
  fgVariableNames[kPairLegPtSum] = "Pair leg p_{T,1} + p_{T,2}"; fgVariableUnits[kPairLegPtSum] = "GeV/c";
  fgVariableNames[kPairLegPtMCSum] = "Pair leg p^{MC}_{T,1} + p^{MC}_{T,2}"; fgVariableUnits[kPairLegPtMCSum] = "GeV/c";
  
  
  fgVariableNames[kPtTPC]             = "p_{T}^{TPC}";                  fgVariableUnits[kPtTPC] = "GeV/c";
  fgVariableNames[kPhiTPC]            = "#varphi^{TPC}";                fgVariableUnits[kPhiTPC] = "rad.";
  fgVariableNames[kEtaTPC]            = "#eta^{TPC}";                   fgVariableUnits[kEtaTPC] = "";
  fgVariableNames[kPin]               = "p_{IN}";                       fgVariableUnits[kPin] = "GeV/c";
  fgVariableNames[kDcaXY]             = "DCA_{xy}";                     fgVariableUnits[kDcaXY] = "cm.";  
  fgVariableNames[kDcaZ]              = "DCA_{z}";                      fgVariableUnits[kDcaZ] = "cm.";  
  fgVariableNames[kPairDca]           = "DCA_{pair}";                   fgVariableUnits[kPairDca] = "cm";
  fgVariableNames[kPairDcaXY]         = "DCA_{xy,pair}";                fgVariableUnits[kPairDcaXY] = "cm";
  fgVariableNames[kPairDcaZ]          = "DCA_{z,pair}";                 fgVariableUnits[kPairDcaZ] = "cm";
  fgVariableNames[kPairDcaSqrt]       = "DCA^{1/2}_{pair}";             fgVariableUnits[kPairDcaSqrt]   = "cm^{1/2}";
  fgVariableNames[kPairDcaXYSqrt]     = "DCA^{1/2}_{xy,pair}";          fgVariableUnits[kPairDcaXYSqrt] = "cm^{1/2}";
  fgVariableNames[kPairDcaZSqrt]      = "DCA^{1/2}_{z,pair}";           fgVariableUnits[kPairDcaZSqrt]   = "cm^{1/2}";
  fgVariableNames[kMassDcaPtCorr]     = "DCA-, p_{T}-corrected mass";           fgVariableUnits[kMassDcaPtCorr]  = "GeV/c^{2}";
  fgVariableNames[kOpAngDcaPtCorr]    = "DCA-, p_{T}-corrected opening angle";  fgVariableUnits[kOpAngDcaPtCorr] = "rad.";
  fgVariableNames[kTrackLength]       = "Track length";                 fgVariableUnits[kTrackLength] = "cm.";
  fgVariableNames[kChi2TPCConstrainedVsGlobal] = "#chi^{2} TPC constrained vs global"; fgVariableUnits[kChi2TPCConstrainedVsGlobal] = "";
  fgVariableNames[kMassUsedForTracking] = "Mass used for tracking"; fgVariableUnits[kMassUsedForTracking] = "GeV/c^{2}";
  fgVariableNames[kITSncls]           = "No.ITS clusters";              fgVariableUnits[kITSncls] = "";
  fgVariableNames[kITSchi2]           = "ITS #chi^{2}";                 fgVariableUnits[kITSchi2] = "";
  fgVariableNames[kITSnclsShared]     = "No.of shared ITS clusters";              fgVariableUnits[kITSnclsShared] = "";
  fgVariableNames[kNclsSFracITS]      = "Fraction of shared ITS clusters/ITS clusters";fgVariableUnits[kNclsSFracITS] = "";
  fgVariableNames[kITSlayerHit]       = "ITS layer";                    fgVariableUnits[kITSlayerHit] = "";
  fgVariableNames[kITSlayerShared]    = "ITS layer";                    fgVariableUnits[kITSlayerShared] = "";
  fgVariableNames[kITSsignal]         = "ITS dE/dx";                    fgVariableUnits[kITSsignal] = "";    
  fgVariableNames[kITSnSig]           = "ITS n_{#sigma}^{e}";           fgVariableUnits[kITSnSig] = "#sigma";
  fgVariableNames[kITSnSig+1]         = "ITS n_{#sigma}^{#pi}";         fgVariableUnits[kITSnSig+1] = "#sigma";
  fgVariableNames[kITSnSig+2]         = "ITS n_{#sigma}^{K}";           fgVariableUnits[kITSnSig+2] = "#sigma";
  fgVariableNames[kITSnSig+3]         = "ITS n_{#sigma}^{p}";           fgVariableUnits[kITSnSig+3] = "#sigma";
  fgVariableNames[kTPCncls]           = "No.TPC clusters";              fgVariableUnits[kTPCncls]           = "";
  fgVariableNames[kTPCchi2]           = "TPC #chi^{2}";                 fgVariableUnits[kTPCchi2]           = "";
  fgVariableNames[kTPCclusBitFired]   = "TPC segment";                  fgVariableUnits[kTPCclusBitFired]   = "";
  fgVariableNames[kTPCNclusBitsFired] = "No.TPC segments";              fgVariableUnits[kTPCNclusBitsFired] = "";
  fgVariableNames[kTPCclustersPerBit] = "No.TPC clusters/segment";      fgVariableUnits[kTPCclustersPerBit] = "";
  fgVariableNames[kTPCcrossedRows]    = "No.TPC crossed rows";          fgVariableUnits[kTPCcrossedRows] = "";
  fgVariableNames[kTPCnclsF]          = "No.TPC findable clusters";     fgVariableUnits[kTPCnclsF] = "";
  fgVariableNames[kTPCnclsShared]     = "No.TPC shared clusters";       fgVariableUnits[kTPCnclsShared] = "";
  fgVariableNames[kTPCnclsSharedRatio] = "# TPC shared clusters / all TPC clusters"; fgVariableUnits[kTPCnclsSharedRatio] = "fraction";
  fgVariableNames[kTPCnclsRatio]      = "No.TPC clusters/findable";     fgVariableUnits[kTPCnclsRatio] = "";
  fgVariableNames[kTPCnclsRatio2]     = "No.TPC clusters/crossed rows"; fgVariableUnits[kTPCnclsRatio2] = "";
  fgVariableNames[kTPCcrossedRowsOverFindableClusters] = "Crossed rows / findable clusters"; fgVariableUnits[kTPCcrossedRowsOverFindableClusters] = "";
  fgVariableNames[kTPCnclsRatio3]     = "No.TPC crossed rows/findable clusters"; fgVariableUnits[kTPCnclsRatio3] = "";
  fgVariableNames[kTPCActiveLength]   = "Track length in active TPC";   fgVariableUnits[kTPCActiveLength] = "";
  fgVariableNames[kTPCGeomLength]     = "Track relative geometric length in active TPC";   fgVariableUnits[kTPCGeomLength] = "";
  fgVariableNames[kTPCsignal]         = "TPC dE/dx";                    fgVariableUnits[kTPCsignal] = "";  
  fgVariableNames[kTPCsignalN]        = "No. TPC clusters PID";         fgVariableUnits[kTPCsignalN] = "";  
  fgVariableNames[kTPCdEdxQmax+0]     = "TPC dE/dx Qmax from IROC";        fgVariableUnits[kTPCdEdxQmax+0] = "";
  fgVariableNames[kTPCdEdxQmax+1]     = "TPC dE/dx Qmax from medium OROC";  fgVariableUnits[kTPCdEdxQmax+1] = "";
  fgVariableNames[kTPCdEdxQmax+2]     = "TPC dE/dx Qmax from long OROC"; fgVariableUnits[kTPCdEdxQmax+2] = "";
  fgVariableNames[kTPCdEdxQmax+3]     = "TPC dE/dx Qmax from all OROC";  fgVariableUnits[kTPCdEdxQmax+3] = "";
  fgVariableNames[kTPCdEdxQtot+0]     = "TPC dE/dx Qtot from IROC";        fgVariableUnits[kTPCdEdxQtot+0] = "";
  fgVariableNames[kTPCdEdxQtot+1]     = "TPC dE/dx Qtot from medium OROC";  fgVariableUnits[kTPCdEdxQtot+1] = "";
  fgVariableNames[kTPCdEdxQtot+2]     = "TPC dE/dx Qtot from long OROC"; fgVariableUnits[kTPCdEdxQtot+2] = "";
  fgVariableNames[kTPCdEdxQtot+3]     = "TPC dE/dx Qtot from all OROC";  fgVariableUnits[kTPCdEdxQtot+3] = "";
  fgVariableNames[kTPCdEdxQmaxOverQtot+0] = "TPC dE/dx Qmax/Qtot from IROC";        fgVariableUnits[kTPCdEdxQmaxOverQtot+0] = "";
  fgVariableNames[kTPCdEdxQmaxOverQtot+1] = "TPC dE/dx Qmax/Qtot from medium OROC"; fgVariableUnits[kTPCdEdxQmaxOverQtot+1] = "";
  fgVariableNames[kTPCdEdxQmaxOverQtot+2] = "TPC dE/dx Qmax/Qtot from long OROC";   fgVariableUnits[kTPCdEdxQmaxOverQtot+2] = "";
  fgVariableNames[kTPCdEdxQmaxOverQtot+3] = "TPC dE/dx Qmax/Qtot from all OROC";    fgVariableUnits[kTPCdEdxQmaxOverQtot+3] = "";
  fgVariableNames[kTPCnSig]           = "TPC n_{#sigma}^{e}";           fgVariableUnits[kTPCnSig] = "#sigma";  
  fgVariableNames[kTPCnSig+1]         = "TPC n_{#sigma}^{#pi}";         fgVariableUnits[kTPCnSig+1] = "#sigma";
  fgVariableNames[kTPCnSig+2]         = "TPC n_{#sigma}^{K}";           fgVariableUnits[kTPCnSig+2] = "#sigma";
  fgVariableNames[kTPCnSig+3]         = "TPC n_{#sigma}^{p}";           fgVariableUnits[kTPCnSig+3] = "#sigma";
  fgVariableNames[kTPCnSigCorrected]         = "TPC n_{#sigma}^{e} - corrected";           fgVariableUnits[kTPCnSigCorrected] = "#sigma";  
  fgVariableNames[kTPCnSigCorrected+1]    = "TPC n_{#sigma}^{#pi} - corrected";         fgVariableUnits[kTPCnSigCorrected+1] = "#sigma";
  fgVariableNames[kTPCnSigCorrected+2]    = "TPC n_{#sigma}^{K} - corrected";           fgVariableUnits[kTPCnSigCorrected+2] = "#sigma";
  fgVariableNames[kTPCnSigCorrected+3]    = "TPC n_{#sigma}^{p} - corrected";           fgVariableUnits[kTPCnSigCorrected+3] = "#sigma";
  fgVariableNames[kTOFbeta]           = "TOF #beta";                    fgVariableUnits[kTOFbeta] = "";
  fgVariableNames[kTOFtime]           = "TOF time";                     fgVariableUnits[kTOFtime] = "ps";
  fgVariableNames[kTOFdx]             = "TOF dx";                       fgVariableUnits[kTOFdx] = "mm";
  fgVariableNames[kTOFdz]             = "TOF dz";                       fgVariableUnits[kTOFdz] = "mm";
  fgVariableNames[kTOFmismatchProbability] = "TOF mismatch probab";     fgVariableUnits[kTOFmismatchProbability] = "";
  fgVariableNames[kTOFchi2]           = "TOF #chi^{2}";                 fgVariableUnits[kTOFchi2] = "";
  fgVariableNames[kTOFdeltaBC]        = "TOF delta BC";                 fgVariableUnits[kTOFdeltaBC] = "";
  fgVariableNames[kTOFnSig]           = "TOF n_{#sigma}^{e}";           fgVariableUnits[kTOFnSig] = "#sigma";  
  fgVariableNames[kTOFnSig+1]         = "TOF n_{#sigma}^{#pi}";         fgVariableUnits[kTOFnSig+1] = "#sigma";
  fgVariableNames[kTOFnSig+2]         = "TOF n_{#sigma}^{K}";           fgVariableUnits[kTOFnSig+2] = "#sigma";
  fgVariableNames[kTOFnSig+3]         = "TOF n_{#sigma}^{p}";           fgVariableUnits[kTOFnSig+3] = "#sigma";
  fgVariableNames[kTRDntracklets]             = "No.TRD tracklets";              fgVariableUnits[kTRDntracklets] = "";
  fgVariableNames[kTRDntrackletsPID]          = "No.TRD PID tracklets";          fgVariableUnits[kTRDntrackletsPID] = "";
  fgVariableNames[kTRDpidProbabilitiesLQ1D]   = "TRD LQ1D electron probability"; fgVariableUnits[kTRDpidProbabilitiesLQ1D] = "";  
  fgVariableNames[kTRDpidProbabilitiesLQ1D+1] = "TRD LQ1D pion probability";     fgVariableUnits[kTRDpidProbabilitiesLQ1D+1] = "";
  fgVariableNames[kTRDpidProbabilitiesLQ2D]   = "TRD LQ2D electron probability"; fgVariableUnits[kTRDpidProbabilitiesLQ2D] = "";  
  fgVariableNames[kTRDpidProbabilitiesLQ2D+1] = "TRD LQ2D pion probability";     fgVariableUnits[kTRDpidProbabilitiesLQ2D+1] = "";
  fgVariableNames[kTRDGTUtracklets]             = "No. TRD GTU tracklets";           fgVariableUnits[kTRDGTUtracklets]        = "";
  fgVariableNames[kTRDGTUlayermask]             = "No. TRD GTU layer0 cond";         fgVariableUnits[kTRDGTUlayermask]        = "";
  fgVariableNames[kTRDGTUpt]                    = "No. TRD GTU pt";                  fgVariableUnits[kTRDGTUpt]               = "";
  fgVariableNames[kTRDGTUsagitta]               = "No. TRD GTU sagitta";             fgVariableUnits[kTRDGTUsagitta]          = "";
  fgVariableNames[kTRDGTUPID]                   = "No. TRD GTU PID";                 fgVariableUnits[kTRDGTUPID]              = "";

  fgVariableNames[kEMCALmatchedEnergy]    = "Calo energy";             fgVariableUnits[kEMCALmatchedEnergy] = "GeV";
  fgVariableNames[kEMCALmatchedClusterId] = "matched Calo cluster id"; fgVariableUnits[kEMCALmatchedClusterId] = "";
  fgVariableNames[kEMCALmatchedEOverP]    = "Calo E/p";                fgVariableUnits[kEMCALmatchedEOverP] = "";  
  fgVariableNames[kEMCALmatchedM02]       = "Calo M02";                fgVariableUnits[kEMCALmatchedM02] = "";
  fgVariableNames[kEMCALmatchedM20]       = "Calo M20";                fgVariableUnits[kEMCALmatchedM20] = "";
  fgVariableNames[kEMCALmatchedNCells]          = "Calo No. cells";           fgVariableUnits[kEMCALmatchedNCells] = "";
  fgVariableNames[kEMCALmatchedNMatchedTracks]  = "Calo No. matched tracks";  fgVariableUnits[kEMCALmatchedNMatchedTracks] = "";
  fgVariableNames[kEMCALmatchedDeltaPhi]        = "#Delta#varphi track-cluster";  fgVariableUnits[kEMCALmatchedDeltaPhi] = "";
  fgVariableNames[kEMCALmatchedDeltaEta]        = "#Delta#eta track-cluster";     fgVariableUnits[kEMCALmatchedDeltaEta] = "";
  fgVariableNames[kEMCALmatchedDistance]        = "track-cluster distance";       fgVariableUnits[kEMCALmatchedDistance] = "";
  fgVariableNames[kEMCALmatchedNSigmaElectron]  = "EMCal n_{#sigma}^{e}";         fgVariableUnits[kEMCALmatchedNSigmaElectron] = "#sigma";
  fgVariableNames[kEMCALclusterEnergy]    = "Calo cls. energy";        fgVariableUnits[kEMCALclusterEnergy] = "GeV";
  fgVariableNames[kEMCALclusterDx]        = "Calo cls. dx";            fgVariableUnits[kEMCALclusterDx] = "";  
  fgVariableNames[kEMCALclusterDz]        = "Calo cls. dz";            fgVariableUnits[kEMCALclusterDz] = "";  
  fgVariableNames[kEMCALdetector]         = "Calo detector";           fgVariableUnits[kEMCALdetector] = "";  
  fgVariableNames[kEMCALm20]              = "Cluster short axis M20";  fgVariableUnits[kEMCALm20] = "";  
  fgVariableNames[kEMCALm02]              = "Cluster short axis M02";  fgVariableUnits[kEMCALm02] = "";  
  fgVariableNames[kEMCALdispersion]       = "Cluster dispersion";      fgVariableUnits[kEMCALdispersion] = "";
  fgVariableNames[kEMCALnCells]           = "Cluster No. cells";       fgVariableUnits[kEMCALnCells] = "";
  fgVariableNames[kEMCALnMatchedTracks]   = "Cluster No. matched tracks"; fgVariableUnits[kEMCALnMatchedTracks] = "";
  fgVariableNames[kEMCALclusterPhi]       = "Cluster #varphi";         fgVariableUnits[kEMCALclusterPhi] = "rad.";
  fgVariableNames[kEMCALclusterEta]       = "Cluster #eta";            fgVariableUnits[kEMCALclusterEta] = "";
  fgVariableNames[kTrackingFlag] = "Tracking flag";  fgVariableUnits[kTrackingFlag] = "";
  fgVariableNames[kTrackingStatus] = "Tracking status";  fgVariableUnits[kTrackingStatus] = "";  
  fgVariableNames[kDeltaPhi]              = "#Delta #varphi";             fgVariableUnits[kDeltaPhi]              = "rad.";
  fgVariableNames[kDeltaPhiBoosted]       = "#Delta #varphi";             fgVariableUnits[kDeltaPhiBoosted]       = "rad.";
  fgVariableNames[kDeltaPhiSym]           = "#Delta #varphi";             fgVariableUnits[kDeltaPhiSym]           = "rad.";
  fgVariableNames[kDeltaPhiSymBoosted]    = "#Delta #varphi";             fgVariableUnits[kDeltaPhiSymBoosted]    = "rad.";
  fgVariableNames[kDeltaTheta]            = "#Delta #theta";              fgVariableUnits[kDeltaTheta]            = "rad.";
  fgVariableNames[kDeltaThetaBoosted]     = "#Delta #theta";              fgVariableUnits[kDeltaThetaBoosted]     = "rad.";
  fgVariableNames[kDeltaEta]              = "#Delta #eta";                fgVariableUnits[kDeltaEta]              = "";
  fgVariableNames[kDeltaEtaBoosted]       = "#Delta #eta";                fgVariableUnits[kDeltaEtaBoosted]       = "";
  fgVariableNames[kDeltaEtaAbs]           = "|#Delta #eta|";              fgVariableUnits[kDeltaEtaAbs]           = "";
  fgVariableNames[kDeltaEtaAbsBoosted]    = "|#Delta #eta|";              fgVariableUnits[kDeltaEtaAbsBoosted]    = "";
  fgVariableNames[kTriggerPt]             = "p_{T} trigger particle";     fgVariableUnits[kTriggerPt]             = "GeV/c";
  fgVariableNames[kTriggerRap]            = "#it{y} trigger particle";    fgVariableUnits[kTriggerRap]            = "";
  fgVariableNames[kTriggerRapAbs]         = "|#it{y}| trigger particle";  fgVariableUnits[kTriggerRapAbs]         = "";
  fgVariableNames[kTriggerPseudoProperDecayTime]  = "t";                      fgVariableUnits[kTriggerPseudoProperDecayTime]  = "cm./c";
  fgVariableNames[kTriggerPairTypeSPD]            = "trigger type spd legs";  fgVariableUnits[kTriggerPairTypeSPD]            = "";
  fgVariableNames[kAssociatedPt]          = "p_{T} associated particle";  fgVariableUnits[kAssociatedPt]          = "GeV/c";
  fgVariableNames[kAssociatedPtBoosted]   = "p_{T} associated particle (boosted)"; fgVariableUnits[kAssociatedPtBoosted] = "GeV/c";
  fgVariableNames[kAssociatedPtOverTriggerGammaT] = "p_{T} associated particle / #gamma_{T} trigger particle"; fgVariableUnits[kAssociatedPtOverTriggerGammaT] = "GeV/c";
  fgVariableNames[kTriggerGammaT]         = "#gamma_{T} trigger particle"; fgVariableUnits[kTriggerGammaT]        = "";
  fgVariableNames[kAssociatedEta]         = "#eta associated particle";   fgVariableUnits[kAssociatedEta]         = "";
  fgVariableNames[kAssociatedEtaBoosted]  = "#eta associated particle";   fgVariableUnits[kAssociatedEtaBoosted]  = "";
  fgVariableNames[kAssociatedPhi]         = "#varphi associated particle";fgVariableUnits[kAssociatedPhi]         = "rad.";
  fgVariableNames[kAssociatedPhiBoosted]  = "#varphi associated particle";fgVariableUnits[kAssociatedPhiBoosted]  = "rad.";
  fgVariableNames[kTriggerEff]                            = "#epsilon_{trigger}";                           fgVariableUnits[kTriggerEff]                            = "";
  fgVariableNames[kOneOverTriggerEff]                     = "1/#epsilon_{trigger}";                         fgVariableUnits[kOneOverTriggerEff]                     = "";
  fgVariableNames[kAssocHadronEff]                        = "#epsilon_{assoc}";                             fgVariableUnits[kAssocHadronEff]                        = "";
  fgVariableNames[kOneOverAssocHadronEff]                 = "1/#epsilon_{assoc}";                           fgVariableUnits[kOneOverAssocHadronEff]                 = "";
  fgVariableNames[kTriggerEffTimesAssocHadronEff]         = "#epsilon_{trigger} #times #epsilon_{assoc}";   fgVariableUnits[kTriggerEffTimesAssocHadronEff]         = "";
  fgVariableNames[kOneOverTriggerEffTimesAssocHadronEff]  = "1/#epsilon_{trigger} #times #epsilon_{assoc}"; fgVariableUnits[kOneOverTriggerEffTimesAssocHadronEff]  = "";
}


//_________________________________________________________________
TChain* AliReducedVarManager::GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                       TChain* friendChain/*=0x0*/, const Char_t* friendDir/*=0x0*/, const Char_t* friendFileName/*=0x0*/) {
  //
  // read an ascii file containing a list of root files with reduced trees
  // and build a TChain
  //
  cout << "Creating the data chain from " << filename << " ..." << endl; 
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
    TString frstr = friendDir;

    if(!strstr.Contains(".root")) continue;
    cout << endl << "Adding file " << str << endl;
    chain->Add(str);
    if(friendChain) {
      TObjArray* arr = strstr.Tokenize("/");
      for(Int_t is=arr->GetEntries()-3; is<arr->GetEntries()-1; is++){
        frstr+="/";
        frstr+=arr->At(is)->GetName();
      }
      frstr+=Form("/%s",friendFileName);
      cout << endl << "Adding friend file " << frstr << endl;
      friendChain->Add(frstr.Data());
    }
    ++index;
  }
  inBuf.close();
  entries = chain->GetEntries();
  Long64_t entriesFriend = (friendChain ? friendChain->GetEntries() : 0);
  cout << "AliReducedVarManager::GetChain() Chain entries = " << entries << endl;
  if(friendChain)
    cout << "AliReducedVarManager::GetChain() Friend chain entries = " << entriesFriend << endl;
  if(friendChain && (entries!=entriesFriend)) {
    cout << "AliReducedVarManager::GetChain() The friend chain does not have the same number of entries as the main chain!!!" << endl;
    cout << "                            Check it out and retry !!" << endl;
    return 0x0;
  }
  cout << " done" << endl;
  return chain;
}

//____________________________________________________________________________________
void AliReducedVarManager::SetTPCelectronCorrectionMaps(TH2F* centroidMap, TH2F* widthMap, AliReducedVarManager::Variables varX, AliReducedVarManager::Variables varY) {
   //
   // initialize the electron TPC pid correction maps
   //
   if(varX>kNVars || varX<=kNothing) {
      cout << "AliReducedVarManager::SetTPCelectronCorrectionMaps() The X-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Correction maps not used! Check it out!" << endl;
      return;
   }
   if(varY>kNVars || varY<=kNothing) {
      cout << "AliReducedVarManager::SetTPCelectronCorrectionMaps() The Y-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Correction maps not used! Check it out!" << endl;
      return;
   }
   fgVarDependencyX = varX;
   fgVarDependencyY = varY;
   if(centroidMap) {
     fgTPCelectronCentroidMap = (TH2F*)centroidMap->Clone(Form("AliReducedVarManager_TPCelectronCentroidMap"));
     fgTPCelectronCentroidMap->SetDirectory(0x0);
   }
   if(widthMap) {
     fgTPCelectronWidthMap = (TH2F*)widthMap->Clone(Form("AliReducedVarManager_TPCelectronWidthMap"));
     fgTPCelectronWidthMap->SetDirectory(0x0);
   }
}

//____________________________________________________________________________________
void AliReducedVarManager::SetTPCpidCalibMaps(Int_t pid, THnF* centroidMap, THnF* widthMap, THnI* statusMap) {
   //
   // initialize the 4D calibration histograms for a given particle (pid = 0-electron, 1-pion, 2-proton)
   //
   if(pid<0 || pid>2) {
      cout << "AliReducedVarManager::SetTPCpidCalibMaps(): Invalid pid code used" << endl;
      return;
   }
   fgTPCpidCalibCentroid[pid] = centroidMap;
   fgTPCpidCalibWidth[pid] = widthMap;
   fgTPCpidCalibStatus[pid] = statusMap;
}

//____________________________________________________________________________________
void AliReducedVarManager::SetTPCpidCalibDepVars(Variables vars[]) {
   //
   //
   //
   for(Int_t i=0; i<4; ++i) {
      if(vars[i]>kNVars || vars[i]<=kNothing) {
         cout << "AliReducedVarManager::SetTPCpidCalibDepVars(): Invalid dependence variable used" << endl;
         return;
      }
      fgTPCpidCalibVars[i] = vars[i];
      fgUsedVars[vars[i]] = kTRUE;
      if(vars[i]==kTPCpileupZA || vars[i]==kTPCpileupZC) {
         fgUsedVars[kTPCpileupZA] = kTRUE;
         fgUsedVars[kTPCpileupZC] = kTRUE;
      }
      if(vars[i]==kTPCpileupContributorsA || vars[i]==kTPCpileupContributorsC) {
         fgUsedVars[kTPCpileupContributorsA] = kTRUE;
         fgUsedVars[kTPCpileupContributorsC] = kTRUE;
      }
   }
}

//____________________________________________________________________________________
void AliReducedVarManager::SetPairEfficiencyMap(TH1* map, AliReducedVarManager::Variables varX, AliReducedVarManager::Variables varY/*=kNothing*/, AliReducedVarManager::Variables varZ/*=kNothing*/) {
  //
  // initialize the pair efficiency map
  //
  if (!map) {
    cout << "AliReducedVarManager::SetPairEfficiencyMap() No efficiency map provided!" << endl;
    return;
  }
  if(varX>kNVars || varX<=kNothing) {
    cout << "AliReducedVarManager::SetPairEfficiencyMap() The X-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
    cout << "                           Efficiency map not used! Check it out!" << endl;
    return;
  }
  if (map->IsA()==TH3F::Class() || map->IsA()==TH3D::Class()) {
    if (varY>kNVars || varY<=kNothing) {
      cout << "AliReducedVarManager::SetPairEfficiencyMap() The Y-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Efficiency map not used! Check it out!" << endl;
      return;
    }
    if (varZ>kNVars || varZ<=kNothing) {
      cout << "AliReducedVarManager::SetPairEfficiencyMap() The Z-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Efficiency map not used! Check it out!" << endl;
      return;
    }
    fgEffMapVarDependencyX  = varX;
    fgEffMapVarDependencyY  = varY;
    fgEffMapVarDependencyZ  = varZ;
    fgPairEffMap            = (TH3F*)map->Clone(Form("AliReducedVarManager_PairEffMap"));
  } else if (map->IsA()==TH2F::Class() || map->IsA()==TH2D::Class()) {
    if (varY>kNVars || varY<=kNothing) {
      cout << "AliReducedVarManager::SetPairEfficiencyMap() The Y-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Efficiency map not used! Check it out!" << endl;
      return;
    }
    fgEffMapVarDependencyX  = varX;
    fgEffMapVarDependencyY  = varY;
    fgPairEffMap            = (TH2F*)map->Clone(Form("AliReducedVarManager_PairEffMap"));
  } else {
    fgEffMapVarDependencyX  = varX;
    fgPairEffMap            = (TH1F*)map->Clone(Form("AliReducedVarManager_PairEffMap"));
  }
  fgPairEffMap->SetDirectory(0x0);
}

//____________________________________________________________________________________
void AliReducedVarManager::SetPairEfficiencyMapDependeciesCorrelation(AliReducedVarManager::Variables varX, AliReducedVarManager::Variables varY/*=kNothing*/, AliReducedVarManager::Variables varZ/*=kNothing*/) {
  //
  // initialize variable dependencies for pair efficiency map used in correlation
  //
  if (!fgPairEffMap) {
    cout << "AliReducedVarManager::SetPairEfficiencyMapDependeciesCorrelation() No efficiency map provided!" << endl;
    return;
  }
  if(varX>kNVars || varX<=kNothing) {
    cout << "AliReducedVarManager::SetPairEfficiencyMapDependeciesCorrelation() The X-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
    cout << "                           Efficiency map not used in correlation! Check it out!" << endl;
    return;
  }
  fgEffMapVarDependencyXCorr = varX;
  if (varY>kNothing && varY<=kNVars) fgEffMapVarDependencyYCorr = varY;
  if (varZ>kNothing && varZ<=kNVars) fgEffMapVarDependencyZCorr = varZ;
}


//____________________________________________________________________________________
void AliReducedVarManager::SetAssociatedHadronEfficiencyMap(TH1* map, AliReducedVarManager::Variables varX, AliReducedVarManager::Variables varY/*=kNothing*/, AliReducedVarManager::Variables varZ/*=kNothing*/) {
  //
  // initialize the associated hadron efficiency map, used for correlation analysis
  //
  if (!map) {
    cout << "AliReducedVarManager::SetAssociatedHadronEfficiencyMap() No efficiency map provided!" << endl;
    return;
  }
  if (fgAssocHadronEffMap) {
    cout << "AliReducedVarManager::SetAssociatedHadronEfficiencyMap() Efficiency map already defined!" << endl;
    return;
  }
  if (varX>kNVars || varX<=kNothing) {
    cout << "AliReducedVarManager::SetAssociatedHadronEfficiencyMap() The X-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
    cout << "                           Efficiency map not used! Check it out!" << endl;
    return;
  }
  if (map->IsA()==TH3F::Class() || map->IsA()==TH3D::Class()) {
    if (varY>kNVars || varY<=kNothing) {
      cout << "AliReducedVarManager::SetAssociatedHadronEfficiencyMap() The Y-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Efficiency map not used! Check it out!" << endl;
      return;
    }
    if (varZ>kNVars || varZ<=kNothing) {
      cout << "AliReducedVarManager::SetAssociatedHadronEfficiencyMap() The Z-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Efficiency map not used! Check it out!" << endl;
      return;
    }
    fgAssocHadronEffMapVarDependencyX   = varX;
    fgAssocHadronEffMapVarDependencyY   = varY;
    fgAssocHadronEffMapVarDependencyZ   = varZ;
    fgAssocHadronEffMap                 = (TH3F*)map->Clone(Form("AliReducedVarManager_AssocHadronEffMap"));
  } else if (map->IsA()==TH2F::Class() || map->IsA()==TH2D::Class()) {
    if (varY>kNVars || varY<=kNothing) {
      cout << "AliReducedVarManager::SetAssociatedHadronEfficiencyMap() The Y-dependency variable is not a valid variable defined in AliReducedVarManager" << endl;
      cout << "                           Efficiency map not used! Check it out!" << endl;
      return;
    }
    fgAssocHadronEffMapVarDependencyX   = varX;
    fgAssocHadronEffMapVarDependencyY   = varY;
    fgAssocHadronEffMap                 = (TH2F*)map->Clone(Form("AliReducedVarManager_AssocHadronEffMap"));
  } else {
    fgAssocHadronEffMapVarDependencyX   = varX;
    fgAssocHadronEffMap                 = (TH1F*)map->Clone(Form("AliReducedVarManager_AssocHadronEffMap"));
  }
  fgAssocHadronEffMap->SetDirectory(0x0);
}

//____________________________________________________________________________________
void AliReducedVarManager::SetLHCDataInfo(TH1F* totalLumi, TH1F* totalInt0, TH1F* totalInt1, TH1I* fillNumber) {
   //
   // initialize the LHC data histograms
   //
   if(totalLumi) {
      fgRunTotalLuminosity = (TH1F*)totalLumi->Clone("AliReducedVarManager_TotalLuminosity");
      fgRunTotalLuminosity->SetDirectory(0x0);
   }
   if(totalInt0) {
      fgRunTotalIntensity0 = (TH1F*)totalInt0->Clone("AliReducedVarManager_TotalIntensity0");
      fgRunTotalIntensity0->SetDirectory(0x0);
   }
   if(totalInt1) {
      fgRunTotalIntensity1 = (TH1F*)totalInt1->Clone("AliReducedVarManager_TotalIntensity1");
      fgRunTotalIntensity1->SetDirectory(0x0);
   }
   if(fillNumber) {
      fgRunLHCFillNumber = (TH1I*)fillNumber->Clone("AliReducedVarManager_LHCFillNumber");
      fgRunLHCFillNumber->SetDirectory(0x0);
   }
}

//____________________________________________________________________________________
void AliReducedVarManager::SetGRPDataInfo(TH1I* dipolePolarity, TH1I* l3Polarity, TH1I* timeStart, TH1I* timeStop) {
   //
   // initialize the GRP data histograms
   //
   if(dipolePolarity) {
      fgRunDipolePolarity = (TH1I*)dipolePolarity->Clone("AliReducedVarManager_DipolePolarity");
      fgRunDipolePolarity->SetDirectory(0x0);
   }
   if(l3Polarity) {
      fgRunL3Polarity = (TH1I*)l3Polarity->Clone("AliReducedVarManager_L3Polarity");
      fgRunL3Polarity->SetDirectory(0x0);
   }
   if(timeStart) {
      fgRunTimeStart = (TH1I*)timeStart->Clone("AliReducedVarManager_TimeStart");
      fgRunTimeStart->SetDirectory(0x0);
   }
   if(timeStop) {
      fgRunTimeEnd = (TH1I*)timeStop->Clone("AliReducedVarManager_TimeEnd");
      fgRunTimeEnd->SetDirectory(0x0);
   }
}

//____________________________________________________________________________________
void AliReducedVarManager::SetupGRPinformation(const Char_t* filename) {
   //
   // open the root file containing GRP information and initialize needed objects
   //
   fgGRPfile = new TFile(filename);
   fgRunTotalLuminosity = (TH1F*)fgGRPfile->Get("lumiTotal")->Clone("AliReducedVarManager_TotalLuminosity");
   fgRunTotalLuminosity->SetDirectory(0x0);
   fgRunTotalIntensity0 = (TH1F*)fgGRPfile->Get("intTotal0")->Clone("AliReducedVarManager_TotalIntensity0");
   fgRunTotalIntensity0->SetDirectory(0x0);
   fgRunTotalIntensity1 = (TH1F*)fgGRPfile->Get("intTotal1")->Clone("AliReducedVarManager_TotalIntensity1");
   fgRunTotalIntensity1->SetDirectory(0x0);
   fgRunLHCFillNumber = (TH1I*)fgGRPfile->Get("fillNumber")->Clone("AliReducedVarManager_LHCFillNumber");
   fgRunLHCFillNumber->SetDirectory(0x0);
   fgRunDipolePolarity = (TH1I*)fgGRPfile->Get("dipolePolarity")->Clone("AliReducedVarManager_DipolePolarity");
   fgRunDipolePolarity->SetDirectory(0x0);
   fgRunL3Polarity = (TH1I*)fgGRPfile->Get("l3Polarity")->Clone("AliReducedVarManager_L3Polarity");
   fgRunL3Polarity->SetDirectory(0x0);
   fgRunTimeStart = (TH1I*)fgGRPfile->Get("timeStart")->Clone("AliReducedVarManager_TimeStart");
   fgRunTimeStart->SetDirectory(0x0);
   fgRunTimeEnd = (TH1I*)fgGRPfile->Get("timeEnd")->Clone("AliReducedVarManager_TimeEnd");
   fgRunTimeEnd->SetDirectory(0x0);
}

//____________________________________________________________________________________
void AliReducedVarManager::SetRunNumbers( TString runNumbers ){
  TObjArray* runNumbersArr = runNumbers.Tokenize(";");
  runNumbersArr->SetOwner();
  for( Int_t iRun=0; iRun < runNumbersArr->GetEntries(); ++iRun){
    TString runNumberString = runNumbersArr->At(iRun)->GetName();
    fgRunNumbers.push_back( runNumberString.Atoi() );
  }
}
//___________________________________________________________________________________
void AliReducedVarManager::Set_pT_re_weight(TH1F *pT_centrality_0_10)
{
  if(pT_centrality_0_10)
    fgpT_reweight_centrality_0_10=pT_centrality_0_10;
  else
    {
      cout<< "There is no pT reweight spectra, please double check !!!"<<endl;
    }
}
//____________________________________________________________________________________
Double_t AliReducedVarManager::calculate_weight_factor( Double_t McpT, Double_t Centrality)
{
  // currently only one pT shape was used, this can be used for pp, more centralitys bins for PbPb
  if(!fgpT_reweight_centrality_0_10){
    cout<< "No pT reweight spectra !!"<<endl;
    return 0.0;
  }

  if(McpT<=fgpT_reweight_centrality_0_10->GetBinCenter(fgpT_reweight_centrality_0_10->GetNbinsX()))
    {
      Int_t bin=fgpT_reweight_centrality_0_10->FindBin(McpT);
      Double_t weight= fgpT_reweight_centrality_0_10->GetBinContent(bin);
      return weight;
    }
  else
    {
      cout<< " pt out the range of the reweight "<<endl;
      return 0.0;
    }
}
//____________________________________________________________________________________
void AliReducedVarManager::SetMultiplicityProfile(TH2* profile, Int_t estimator) {
   //
   // initialize the profile for the z-vertex equalization of the multiplicity estimator
   //
  Int_t iEstimator = estimator - kMultiplicity;
  if( iEstimator >= kNMultiplicityEstimators ){
    cout << "Multiplcity estimator " << estimator << " not defined!" <<endl;
    return;
  }
  if(!profile){
    cout <<"AliReducedVarManager::SetMultiplicityProfile : Profile null!"  << endl;
    return;
  }
  fgAvgMultVsVtxAndRun[iEstimator] = (TH2*)profile->Clone( Form("profile_%d", estimator  ));
  fgAvgMultVsVtxAndRun[iEstimator]->SetDirectory(0x0);
}

//____________________________________________________________________________________
void AliReducedVarManager::SetVZEROCalibrationPath(const Char_t* path) {
   //
   // initialize the path to the VZERO calibration histograms
   //
   fgVZEROCalibrationPath = path;
}

//____________________________________________________________________________________
void AliReducedVarManager::SetCalibrateVZEROqVector(Bool_t option) {
   //
   // set the option whether to calibrate the VZERO event plane
   //
   fgOptionCalibrateVZEROqVec = option;
}

//____________________________________________________________________________________
void AliReducedVarManager::SetRecenterVZEROqVector(Bool_t option) {
   //
   // set the option whether to recenter the VZERO event plane
   //   
   fgOptionRecenterVZEROqVec = option;
   //if(fgOptionRecenterVZEROqVec) fgOptionCalibrateVZEROqVec = kTRUE;
}

void AliReducedVarManager::SetRecenterTPCqVector(Bool_t option) {
   //
   // set the option whether to recenter the TPC event plane
   //   
   fgOptionRecenterTPCqVec = option;
  
}

void AliReducedVarManager::SetEventResolution(Bool_t option) {
   //
   // set the option whether to divide by resolution the flow coefficients
   //   
   fgOptionEventRes = option;
  
}

//____________________________________________________________________________________
Int_t AliReducedVarManager::GetCorrectedMultiplicity( Int_t estimator,  Int_t correction, Int_t reference, Int_t smearing ){
  //
  // Return the index of the multiplicity estimator, for given
  // - estimator
  // - correction ( correction wrt. vertex, vertex and run number (1D or 2D) )
  // - reference bin (bin with maximum or minimum efficiency)
  // - smearing method
  //
  Int_t iEstimator = estimator - kMultiplicity;
  Int_t nPerEstimator = kNCorrections * kNReferenceMultiplicities * kNSmearingMethods;
  Int_t ret = kCorrectedMultiplicity + iEstimator * nPerEstimator;
  ret += correction * kNReferenceMultiplicities * kNSmearingMethods;
  ret += reference * kNSmearingMethods;
  ret += smearing;
  return ret;
}


//____________________________________________________________________________________
AliKFParticle AliReducedVarManager::BuildKFcandidate(TRACK* track1, Float_t mh1, TRACK* track2, Float_t mh2) {
   //
   // build a KF pair from the 2 legs
   //
 Double_t p[6]; Double_t cov[21];
 for(int i=0; i<6; i++) p[i] = track1->TrackParam(i);
 for(int i=0; i<21; i++) cov[i] = track1->CovMatrix(i);
 AliKFParticle kfPart1; 
 kfPart1.Initialize();
 kfPart1.Create(p,cov,track1->Charge(),mh1);
 //
 for(int i=0; i<6; i++) p[i] = track2->TrackParam(i);
 for(int i=0; i<21; i++) cov[i] = track2->CovMatrix(i);
 AliKFParticle kfPart2;  
 kfPart2.Initialize();
 kfPart2.Create(p,cov,track2->Charge(),mh2);
 //
 AliKFParticle pairKF; 
 pairKF.Initialize();
 pairKF.AddDaughter(kfPart1);
 pairKF.AddDaughter(kfPart2);
 
 return pairKF;
}


//____________________________________________________________________________________
AliKFParticle AliReducedVarManager::BuildKFvertex( AliReducedEventInfo * event )
{
   // 
   // build an AliKF vertex using the cov matrix 
   //
   AliKFParticle kVtx; kVtx.Initialize();
   for(int i=0; i<3; i++) kVtx.Parameter(i) = event->Vertex(i);
   for(int i=0; i<6; i++) kVtx.Covariance(i) = event->VertexCovMatrix(i);
   kVtx.Chi2() = 1.; // FIX THIS! to be added in AliReducedEventInfo
   kVtx.NDF() = 2*event->VertexNContributors() - 3;
   // printf("ndf %d %d %f \n",kVtx.GetNDF(),kVtx.GetQ(),kVtx.GetParameter(0)); getchar();
 
   return kVtx;
}
