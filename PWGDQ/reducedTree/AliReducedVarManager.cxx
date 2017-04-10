/*
***********************************************************
  Implementation of AliVarManager
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
#include <TFile.h>
#include <THashList.h>

#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedEventPlaneInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliReducedCaloClusterInfo.h"

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
    "TRD",                                                   "N/A",                                                            "FastOnly",                               "N/A",         
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
const Double_t AliReducedVarManager::fgkVZEROminMult = 0.5;   // minimum VZERO channel multiplicity
const Float_t  AliReducedVarManager::fgkTPCQvecRapGap = 0.8;    // symmetric interval in the middle of the TPC excluded from EP calculation
      Float_t  AliReducedVarManager::fgBeamMomentum = 1380.;   // beam momentum in GeV/c
     
Int_t      AliReducedVarManager::fgCurrentRunNumber = -1;
TString AliReducedVarManager::fgVariableNames[AliReducedVarManager::kNVars] = {""};
TString AliReducedVarManager::fgVariableUnits[AliReducedVarManager::kNVars] = {""};
AliReducedBaseEvent* AliReducedVarManager::fgEvent = 0x0;
AliReducedEventPlaneInfo* AliReducedVarManager::fgEventPlane = 0x0;
Bool_t AliReducedVarManager::fgUsedVars[AliReducedVarManager::kNVars] = {kFALSE};
TH2F* AliReducedVarManager::fgTPCelectronCentroidMap = 0x0;
TH2F* AliReducedVarManager::fgTPCelectronWidthMap = 0x0;
AliReducedVarManager::Variables AliReducedVarManager::fgVarDependencyX = kNothing;
AliReducedVarManager::Variables AliReducedVarManager::fgVarDependencyY = kNothing;
TH1F* AliReducedVarManager::fgRunTotalLuminosity = 0x0;
TH1F* AliReducedVarManager::fgRunTotalIntensity0 = 0x0;
TH1F* AliReducedVarManager::fgRunTotalIntensity1 = 0x0;
TH1I* AliReducedVarManager::fgRunLHCFillNumber = 0x0;
TH1I* AliReducedVarManager::fgRunDipolePolarity = 0x0;
TH1I* AliReducedVarManager::fgRunL3Polarity = 0x0;
TH1I* AliReducedVarManager::fgRunTimeStart = 0x0;
TH1I* AliReducedVarManager::fgRunTimeEnd = 0x0;
std::vector<Int_t>  AliReducedVarManager::fgRunNumbers;
Int_t AliReducedVarManager::fgRunID = -1;
TH1* AliReducedVarManager::fgAvgMultVsVertex[3] = {0x0, 0x0, 0x0};
Double_t AliReducedVarManager::fgRefMult[3] = {0., 0., 0.};
TString AliReducedVarManager::fgVZEROCalibrationPath = "";
TProfile2D* AliReducedVarManager::fgAvgVZEROChannelMult[64] = {0x0};
TProfile2D* AliReducedVarManager::fgVZEROqVecRecentering[4] = {0x0};
Bool_t AliReducedVarManager::fgOptionCalibrateVZEROqVec = kFALSE;
Bool_t AliReducedVarManager::fgOptionRecenterVZEROqVec = kFALSE;

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
  if(fgUsedVars[kRap]) {
    fgUsedVars[kMass] = kTRUE;
    fgUsedVars[kP] = kTRUE;
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
  if(fgUsedVars[kTPCnSigCorrected+kElectron]) {
     fgUsedVars[kTPCnSig+kElectron] = kTRUE; 
     fgUsedVars[fgVarDependencyX] = kTRUE; 
     fgUsedVars[fgVarDependencyY] = kTRUE;
  }
  if(fgUsedVars[kNTracksITSoutVsSPDtracklets] || fgUsedVars[kNTracksTPCoutVsSPDtracklets] ||
     fgUsedVars[kNTracksTOFoutVsSPDtracklets] || fgUsedVars[kNTracksTRDoutVsSPDtracklets])
     fgUsedVars[kSPDntracklets] = kTRUE;
  
  if(fgUsedVars[kRapMC]) fgUsedVars[kMassMC] = kTRUE;

  if(fgUsedVars[kPairPhiV]){
    fgUsedVars[kL3Polarity] = kTRUE;
  }
  if(fgUsedVars[kSPDntrackletsCorr]){
    fgUsedVars[kSPDntracklets] = kTRUE;
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
}

//__________________________________________________________________
void AliReducedVarManager::FillEventInfo(Float_t* values) {
  //
  // Fill event information
  //
  FillEventInfo(fgEvent, values, fgEventPlane);
}

//__________________________________________________________________
void AliReducedVarManager::FillEventInfo(BASEEVENT* baseEvent, Float_t* values, EVENTPLANE* eventF/*=0x0*/) {
  //
  // fill event wise info
  //
  // Basic event information
  values[kRunNo]                   = baseEvent->RunNo();
  if(fgUsedVars[kRunID]){
    if( fgRunID < 0 ){
      for( fgRunID = 0; fgRunNumbers[ fgRunID ] != values[kRunNo] && fgRunID< fgRunNumbers.size() ; ++fgRunID );
    }
    values[kRunID] = fgRunID;
  }
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
      }
      cout << "AliReducedVarManager::Info  Loading VZERO calibration and/or recentering parameters for run " << fgCurrentRunNumber << endl;
      if(fgOptionCalibrateVZEROqVec) {
        for(Int_t iCh=0; iCh<64; ++iCh) {
           
           fgAvgVZEROChannelMult[iCh] = (TProfile2D*)calibList->FindObject(Form("VZEROmult_ch%d_VtxCent_prof", iCh))->Clone(Form("run%d_ch%d", fgCurrentRunNumber, iCh));
           fgAvgVZEROChannelMult[iCh]->SetDirectory(0x0);
        }
      }
      if(fgOptionRecenterVZEROqVec) {
         fgVZEROqVecRecentering[0] = (TProfile2D*)calibList->FindObject(Form("QvecX_sideA_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecX_VZEROA", fgCurrentRunNumber));
         fgVZEROqVecRecentering[0]->SetDirectory(0x0);
         fgVZEROqVecRecentering[1] = (TProfile2D*)calibList->FindObject(Form("QvecY_sideA_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecY_VZEROA", fgCurrentRunNumber));
         fgVZEROqVecRecentering[1]->SetDirectory(0x0);
         fgVZEROqVecRecentering[2] = (TProfile2D*)calibList->FindObject(Form("QvecX_sideC_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecX_VZEROC", fgCurrentRunNumber));
         fgVZEROqVecRecentering[2]->SetDirectory(0x0);
         fgVZEROqVecRecentering[3] = (TProfile2D*)calibList->FindObject(Form("QvecY_sideC_h2_CentSPDVtxZ_prof"))->Clone(Form("run%d_QvecY_VZEROC", fgCurrentRunNumber));
         fgVZEROqVecRecentering[3]->SetDirectory(0x0);
      }
      calibFile->Close();
    }
  }
  
  values[kEventNumberInFile]    = event->EventNumberInFile();
  values[kBC]                   = event->BC();
  values[kTimeStamp]            = event->TimeStamp();
  if(fgUsedVars[kTimeRelativeSOR]) values[kTimeRelativeSOR] = (event->TimeStamp() - values[kRunTimeStart]) / 60.;
  if(fgUsedVars[kTimeRelativeSORfraction] && 
     (values[kRunTimeEnd]-values[kRunTimeStart])>1.)   // the run should be longer than 1 second ... 
    values[kTimeRelativeSORfraction] = (event->TimeStamp() - values[kRunTimeStart]) / (values[kRunTimeEnd] - values[kRunTimeStart]);
  values[kEventType]            = event->EventType();
  values[kTriggerMask]          = event->TriggerMask();
  values[kINT7Triggered]        = event->TriggerMask() & kINT7 ?1:0;
  values[kHighMultV0Triggered]  = event->TriggerMask() & kHighMultV0 ?1:0;
  values[kIsPhysicsSelection]   = (event->IsPhysicsSelection() ? 1.0 : 0.0);
  values[kIsSPDPileup]          = event->IsSPDPileup();
  values[kIsSPDPileup5]         = event->EventTag(11);
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
  
  for(Int_t iflag=0;iflag<32;++iflag) 
    values[kNTracksPerTrackingStatus+iflag] = event->TracksPerTrackingFlag(iflag);
  
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
  
  values[kSPDntracklets]   = event->SPDntracklets();
  if( fgUsedVars[kSPDntrackletsCorr] || fgUsedVars[kSPDntrackletsCorrSmear] ){
      values[kSPDntrackletsCorr] = values[kSPDntracklets];
      values[kSPDntrackletsCorrSmear] = values[kSPDntracklets];
      if( fgAvgMultVsVertex[0] && TMath::Abs(baseEvent->Vertex(2)) < 10. ){
        if( !fgRefMult[0] ) fgRefMult[0] = fgAvgMultVsVertex[0]->GetMaximum();
        Double_t localAvg = fgAvgMultVsVertex[0]->GetBinContent( fgAvgMultVsVertex[0]->FindBin(baseEvent->Vertex(2)) );
        Double_t deltaM = values[kSPDntracklets] * (fgRefMult[0]/localAvg - 1);
        values[kSPDntrackletsCorrSmear] += (deltaM>0 ? 1. : -1.) * gRandom->Poisson(TMath::Abs(deltaM));
        values[kSPDntrackletsCorr]      *= fgRefMult[0]/localAvg;
      }
  }
  fgUsedVars[kNTracksITSoutVsSPDtracklets] = kTRUE;  
  fgUsedVars[kNTracksTPCoutVsSPDtracklets] = kTRUE;
  fgUsedVars[kNTracksTRDoutVsSPDtracklets] = kTRUE;
  fgUsedVars[kNTracksTOFoutVsSPDtracklets] = kTRUE;
  if(values[kSPDntracklets]>0.01) {
    values[kNTracksITSoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kITSout]/values[kSPDntracklets];
    values[kNTracksTPCoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kTPCout]/values[kSPDntracklets];
    values[kNTracksTRDoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kTRDout]/values[kSPDntracklets];
    values[kNTracksTOFoutVsSPDtracklets] = values[kNTracksPerTrackingStatus+kTOFout]/values[kSPDntracklets];
  }
  else {
     fgUsedVars[kNTracksITSoutVsSPDtracklets] = kFALSE;  
     fgUsedVars[kNTracksTPCoutVsSPDtracklets] = kFALSE;
     fgUsedVars[kNTracksTRDoutVsSPDtracklets] = kFALSE;
     fgUsedVars[kNTracksTOFoutVsSPDtracklets] = kFALSE;
  }
    
  values[kNCaloClusters]   = event->GetNCaloClusters();
  values[kNTPCclusters]    = event->NTPCClusters();
  values[kSPDntracklets08] = 0.;
  values[kSPDntracklets16] = 0.;
  values[kSPDntrackletsOuter] = 0.;
  for(Int_t ieta=0;ieta<32;++ieta) {
    values[kSPDntrackletsEta+ieta] = event->SPDntracklets(ieta);
    values[kSPDntracklets16] += event->SPDntracklets(ieta);
    if(ieta >7 && ieta <24 ) values[kSPDntracklets08] += event->SPDntracklets(ieta);
    if( ieta <7 || ieta > 24 ) values[kSPDntrackletsOuter] += event->SPDntracklets(ieta);
  }
  if( fgUsedVars[kSPDntrackletsOuterCorr] || fgUsedVars[kSPDntrackletsOuterCorrSmear] ){
      values[kSPDntrackletsOuterCorr] = values[kSPDntrackletsOuter];
      values[kSPDntrackletsOuterCorrSmear] = values[kSPDntrackletsOuter];
      if( fgAvgMultVsVertex[1] && TMath::Abs(baseEvent->Vertex(2)) < 10. ){
        if( !fgRefMult[1] ) fgRefMult[1] = fgAvgMultVsVertex[1]->GetMaximum();
        Double_t localAvg = fgAvgMultVsVertex[1]->GetBinContent( fgAvgMultVsVertex[1]->FindBin(baseEvent->Vertex(2)) );
        Double_t deltaM = values[kSPDntrackletsOuter] * (fgRefMult[1]/localAvg - 1);
        values[kSPDntrackletsOuterCorrSmear] += (deltaM>0 ? 1. : -1.) * gRandom->Poisson(TMath::Abs(deltaM));
        values[kSPDntrackletsOuterCorr]      *= fgRefMult[1]/localAvg;
      }
  }
  
  
  for(Int_t i=0;i<2;++i) values[kSPDFiredChips+i] = event->SPDFiredChips(i+1);
  for(Int_t i=0;i<6;++i) values[kITSnClusters+i] = event->ITSClusters(i+1);
  values[kSPDnSingleClusters] = event->SPDnSingleClusters();

  //VZERO detector information
  values[kVZEROATotalMult] = event->MultVZEROA();
  values[kVZEROCTotalMult] = event->MultVZEROC();
  values[kVZEROTotalMult]  = event->MultVZERO();
  
    if( fgUsedVars[kVZEROTotalMultCorr] || fgUsedVars[kVZEROTotalMultCorrSmear] ){
      values[kVZEROTotalMultCorr] = values[kVZEROTotalMult];
      values[kVZEROTotalMultCorrSmear] = values[kVZEROTotalMult];
      if( fgAvgMultVsVertex[1] && TMath::Abs(baseEvent->Vertex(2)) < 10. ){
        if( !fgRefMult[2] ) fgRefMult[2] = fgAvgMultVsVertex[2]->GetMaximum();
        Double_t localAvg = fgAvgMultVsVertex[2]->GetBinContent( fgAvgMultVsVertex[2]->FindBin(baseEvent->Vertex(2)) );
        Double_t deltaM = values[kVZEROTotalMult] * (fgRefMult[2]/localAvg - 1);
        values[kVZEROTotalMultCorrSmear] += (deltaM>0 ? 1. : -1.) * gRandom->Poisson(TMath::Abs(deltaM));
        values[kVZEROTotalMultCorr]      *= fgRefMult[2]/localAvg;
      }
  }
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
  
  fgUsedVars[kNTracksTPCoutFromPileup] = kTRUE;
  if(values[kVZEROTotalMult]>0.0)
     values[kNTracksTPCoutFromPileup] = values[kNTracksPerTrackingStatus+kTPCout] - (-2.55+TMath::Sqrt(2.55*2.55+4.0e-5*values[kVZEROTotalMult]))/2.0e-5;
  else fgUsedVars[kNTracksTPCoutFromPileup] = kFALSE;
  
  if(!eventF && (fgUsedVars[kVZEROQvecX+0*6+1] || fgUsedVars[kVZEROQvecY+0*6+1] || fgUsedVars[kVZERORP+0*6+1])) {
    Double_t qvecVZEROA[EVENTPLANE::fgkNMaxHarmonics][2] = {{0.0}};
    Double_t qvecVZEROC[EVENTPLANE::fgkNMaxHarmonics][2] = {{0.0}};
    if(fgOptionCalibrateVZEROqVec && fgAvgVZEROChannelMult[0]) {
      Float_t calibVZEROMult[64] = {0.};
      for(Int_t iCh=0; iCh<64; ++iCh) {
         if(event->MultChannelVZERO(iCh)>=fgkVZEROminMult) {
            Float_t avMult = fgAvgVZEROChannelMult[iCh]->GetBinContent(fgAvgVZEROChannelMult[iCh]->FindBin(event->Vertex(2), event->CentralitySPD()));
            calibVZEROMult[iCh] = event->MultChannelVZERO(iCh) / (avMult>1.0e-6 ? avMult : 1.0);
         }
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
       qvecVZEROA[1][0] -= recenterOffset;
       recenterOffset = fgVZEROqVecRecentering[1]->GetBinContent(fgVZEROqVecRecentering[1]->FindBin(event->CentralitySPD(), event->Vertex(2)));
       qvecVZEROA[1][1] -= recenterOffset;
       recenterOffset = fgVZEROqVecRecentering[2]->GetBinContent(fgVZEROqVecRecentering[2]->FindBin(event->CentralitySPD(), event->Vertex(2)));
       qvecVZEROC[1][0] -= recenterOffset;
       recenterOffset = fgVZEROqVecRecentering[3]->GetBinContent(fgVZEROqVecRecentering[3]->FindBin(event->CentralitySPD(), event->Vertex(2)));
       qvecVZEROC[1][1] -= recenterOffset;
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
void AliReducedVarManager::FillL0TriggerInputs(EVENT* event, Int_t input, Float_t* values) {
  //
  // fill the L0 trigger inputs
  //
  values[kL0TriggerInput] = -1.0;
  if(fgUsedVars[kL0TriggerInput] && event->L0TriggerInput(input)) values[kL0TriggerInput] = input;
}


//_________________________________________________________________
void AliReducedVarManager::FillL1TriggerInputs(EVENT* event, Int_t input, Float_t* values) {
  //
  // fill the L1 trigger inputs
  //
  values[kL1TriggerInput] = -1.0;
  if(fgUsedVars[kL1TriggerInput] && event->L1TriggerInput(input)) values[kL1TriggerInput] = input;
}

//_________________________________________________________________
void AliReducedVarManager::FillL2TriggerInputs(EVENT* event, Int_t input, Float_t* values) {
  //
  // fill the L2 trigger inputs
  //
  values[kL2TriggerInput] = -1.0;
  if(fgUsedVars[kL2TriggerInput] && event->L2TriggerInput(input)) values[kL2TriggerInput] = input;
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
void AliReducedVarManager::FillTrackingFlags(TRACK* p, Float_t* values) {
  //
  // Fill tracking flags
  //
  for(Int_t i=0; i<kNTrackingFlags; i++){
    values[kTrackingFlags+i] = p->TestFlag(i)+10e-6;
  }
}



//_________________________________________________________________
void AliReducedVarManager::FillTrackingFlag(TRACK* track, UInt_t flag, Float_t* values) {
  //
  // fill the tracking flag
  //
   //cout << "AliReducedVarManager::FillTrackingFlag track/flag/status :: " << track << "/" << flag << "/" << track->CheckTrackStatus(flag) << endl;
  values[kTrackingFlag] = -1;
  if(track->CheckTrackStatus(flag)) values[kTrackingFlag] = flag;
  //cout << "AliReducedVarManager::FillTrackingFlag values[kTrackingFlag] :: " << values[kTrackingFlag] << endl;
}

//_________________________________________________________________
void AliReducedVarManager::FillTrackQualityFlag(BASETRACK* track, UShort_t flag, Float_t* values) {
  //
  // fill the track quality flag
  //
  values[kTrackQualityFlag] = -1;
  if(track->TestQualityFlag(flag)) values[kTrackQualityFlag] = flag;
}

//_________________________________________________________________
void AliReducedVarManager::FillPairQualityFlag(PAIR* p, UShort_t flag, Float_t* values) {
  //
  // fill the pair quality flag
  //
  values[kPairQualityFlag] = -1;
  if(p->TestQualityFlag(flag)) values[kPairQualityFlag] = flag;
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
void AliReducedVarManager::FillEventOnlineTrigger(UShort_t triggerBit, Float_t* values) {
  //
  // fill the trigger bit input
  //
  if(triggerBit>=64) return;
  if(!fgEvent) return;
  ULong64_t trigger = 1;
  values[kOnlineTrigger] = triggerBit;
  values[kOnlineTriggerFired] = (((AliReducedEventInfo*)fgEvent)->TriggerMask()&(trigger<<triggerBit) ? 1.0 : 0.0);
  values[kOnlineTriggerFired2] = (values[kOnlineTriggerFired]>0.01 ? triggerBit : -1.0); 
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
   if(fgUsedVars[kThetaMC]) values[kThetaMC] = p->ThetaMC();
   if(fgUsedVars[kEtaMC]) values[kEtaMC] = p->EtaMC();
   if(fgUsedVars[kPhiMC]) values[kPhiMC] = p->PhiMC();
   if(fgUsedVars[kMassMC]) {
      if(TMath::Abs(p->MCPdg(0))==443)
      values[kMassMC] = 3.1;   // TODO: use correct PDG mass
   }
   if(fgUsedVars[kRapMC]) {
      if(TMath::Abs(p->MCPdg(0))==443)
         values[kRapMC] = p->RapidityMC(3.1);   // TODO: use correct PDG mass
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
     Float_t tpcEPsubtracted[6] = {0.0};
     Double_t qVec[6][2] = {{0.0}};
     for(Int_t ih=0; ih<6; ++ih) {qVec[ih][0]=values[kTPCQvecXtotal+ih]; qVec[ih][1]=values[kTPCQvecYtotal+ih];}
     EVENT* eventInfo = NULL;
     if(fgEvent->IsA()==EVENT::Class()) eventInfo = (EVENT*)fgEvent;
     if((p->IsA() == AliReducedTrackInfo::Class()) && eventInfo) {
        eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)p,qVec,EVENTPLANE::kTPC,-0.8,-0.5*fgkTPCQvecRapGap);
        eventInfo->SubtractParticleFromQvector((AliReducedTrackInfo*)p,qVec,EVENTPLANE::kTPC,0.5*fgkTPCQvecRapGap,0.8);
     }
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
     for(Int_t ih=0; ih<6;++ih) 
        tpcEPsubtracted[ih] = TMath::ATan2(qVec[ih][1], qVec[ih][0])/Double_t(ih+1);
     
     for(Int_t ih=0; ih<6; ++ih) {
        // vn using Psi_n
        if(fgUsedVars[kTPCFlowVn+ih])
           values[kTPCFlowVn+ih] = TMath::Cos(DeltaPhi(values[kPhi],tpcEPsubtracted[ih])*(ih+1));
        if(fgUsedVars[kTPCFlowSine+ih]) 
           values[kTPCFlowSine+ih] = TMath::Sin(DeltaPhi(values[kPhi],tpcEPsubtracted[ih])*(ih+1));
        if(fgUsedVars[kTPCuQ+ih]) {
           values[kTPCuQ+ih] = TMath::Cos((values[kPhi]-tpcEPsubtracted[ih])*(ih+1));
           values[kTPCuQ+ih] *= TMath::Sqrt(qVec[ih][0]*qVec[ih][0] + qVec[ih][1]*qVec[ih][1]);
        }
        if(fgUsedVars[kTPCuQsine+ih]) {
           values[kTPCuQsine+ih] = TMath::Sin((values[kPhi]-tpcEPsubtracted[ih])*(ih+1));
           values[kTPCuQsine+ih] *= TMath::Sqrt(qVec[ih][0]*qVec[ih][0] + qVec[ih][1]*qVec[ih][1]);
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
  values[kCharge]      = pinfo->Charge();
  
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

  values[kTRDpidProbabilitiesLQ1D]   = pinfo->TRDpidLQ1D(0);
  values[kTRDpidProbabilitiesLQ1D+1] = pinfo->TRDpidLQ1D(1);
  values[kTRDpidProbabilitiesLQ2D]   = pinfo->TRDpidLQ2D(0);
  values[kTRDpidProbabilitiesLQ2D+1] = pinfo->TRDpidLQ2D(1);
  values[kTRDntracklets]    = pinfo->TRDntracklets(0);
  values[kTRDntrackletsPID] = pinfo->TRDntracklets(1);
  
  if(fgUsedVars[kEMCALmatchedEnergy] || fgUsedVars[kEMCALmatchedEOverP]) {
    values[kEMCALmatchedClusterId] = pinfo->CaloClusterId();
    if(fgEvent && (fgEvent->IsA()==EVENT::Class())){
      CLUSTER* cluster = ((EVENT*)fgEvent)->GetCaloCluster(pinfo->CaloClusterId());
      values[kEMCALmatchedEnergy] = (cluster ? cluster->Energy() : -999.0);
      Float_t mom = pinfo->P();
      values[kEMCALmatchedEOverP] = (TMath::Abs(mom)>1.e-8 && cluster ? values[kEMCALmatchedEnergy]/mom : -999.0);
    }
  }  

  FillTrackingStatus(pinfo,values);
  FillTrackingFlags(pinfo,values);
  
  if(pinfo->HasMCTruthInfo()) {
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
  
  if(fgUsedVars[kRap]) values[kRap]              = p->Rapidity();
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

  if(fgUsedVars[kRap]) values[kRap] = p.Rapidity();

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
  
  if(p.PairType()==1 && t1->HasMCTruthInfo() && t2->HasMCTruthInfo()) {
     TRACK* pinfo1 = 0x0;
     if(t1->IsA()==TRACK::Class()) pinfo1 = (TRACK*)t1;
     TRACK* pinfo2 = 0x0;
     if(t2->IsA()==TRACK::Class()) pinfo2 = (TRACK*)t2;
     
     PAIR pMC;
     pMC.PxPyPz(pinfo1->MCmom(0)+pinfo2->MCmom(0), pinfo1->MCmom(1)+pinfo2->MCmom(1), pinfo1->MCmom(2)+pinfo2->MCmom(2));
     pMC.CandidateId(type);
     if(fgUsedVars[kPtMC]) values[kPtMC] = pMC.Pt();
     if(fgUsedVars[kPMC]) values[kPMC] = pMC.P();
     values[kPxMC] = pMC.Px();
     values[kPyMC] = pMC.Py();
     values[kPzMC] = pMC.Pz();
     if(fgUsedVars[kThetaMC]) values[kThetaMC] = pMC.Theta();
     if(fgUsedVars[kEtaMC]) values[kEtaMC] = pMC.Eta();
     if(fgUsedVars[kPhiMC]) values[kPhiMC] = pMC.Phi();
     if(fgUsedVars[kMassMC]) {
        values[kMassMC] = m1*m1+m2*m2 + 
            2.0*(TMath::Sqrt(m1*m1+pinfo1->PMC()*pinfo1->PMC())*TMath::Sqrt(m2*m2+pinfo2->PMC()*pinfo2->PMC()) - 
            pinfo1->MCmom(0)*pinfo2->MCmom(0) - pinfo1->MCmom(1)*pinfo2->MCmom(1) - pinfo1->MCmom(2)*pinfo2->MCmom(2));
         if(values[kMassMC]<0.0) {
            cout << "FillPairInfo(track, track, type, values): Warning: Very small squared mass found. "
            << "   Could be negative due to resolution of Float_t so it will be set to a small positive value." << endl; 
            cout << "   massMC2: " << values[kMassMC] << endl;
            cout << "p1(p,x,y,z): " << pinfo1->PMC() << ", " << pinfo1->MCmom(0) << ", " << pinfo1->MCmom(1) << ", " << pinfo1->MCmom(2) << endl;
            cout << "p2(p,x,y,z): " << pinfo2->PMC() << ", " << pinfo2->MCmom(0) << ", " << pinfo2->MCmom(1) << ", " << pinfo2->MCmom(2) << endl;
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
  if(fgUsedVars[kP]) values[kP] = p.P();
  if(fgUsedVars[kEta]) values[kEta] = p.Eta();
  if(fgUsedVars[kRap]) values[kRap] = p.Rapidity();
  if(fgUsedVars[kPhi]) values[kPhi] = p.Phi();
  if(fgUsedVars[kTheta]) values[kTheta] = p.Theta();
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
void AliReducedVarManager::FillCorrelationInfo(PAIR* p, BASETRACK* t, Float_t* values) {
  //
  // fill pair-track correlation information
  //
  if(fgUsedVars[kDeltaPhi]) {
    Double_t delta = p->Phi() - t->Phi();
    if(delta>3.0/2.0*TMath::Pi()) delta -= 2.0*TMath::Pi();
    if(delta<-0.5*TMath::Pi()) delta += 2.0*TMath::Pi();
    values[kDeltaPhi] = delta;
  }
  
  if(fgUsedVars[kDeltaTheta]) 
    values[kDeltaTheta] = p->Theta() - t->Theta();
  
  if(fgUsedVars[kDeltaEta]) values[kDeltaEta] = p->Eta() - t->Eta();
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
  fgVariableNames[kL0TriggerInput]       = "L0 trigger inputs";               fgVariableUnits[kL0TriggerInput]       = "";
  fgVariableNames[kL1TriggerInput]       = "L1 trigger inputs";               fgVariableUnits[kL1TriggerInput]       = "";
  fgVariableNames[kL2TriggerInput]       = "L2 trigger inputs";               fgVariableUnits[kL2TriggerInput]       = "";
  fgVariableNames[kRunNo]                = "Run number";                      fgVariableUnits[kRunNo]                = "";
  fgVariableNames[kRunID]                = "Run number";                      fgVariableUnits[kRunID]                = "";
  fgVariableNames[kLHCFillNumber]        = "LHC fill number";                 fgVariableUnits[kLHCFillNumber]        = ""; 
  fgVariableNames[kBeamEnergy]           = "Beam energy";                     fgVariableUnits[kBeamEnergy]           = "GeV";
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
  fgVariableNames[kOnlineTriggerFired]   = "Online trigger fired";            fgVariableUnits[kOnlineTriggerFired]   = "";
  fgVariableNames[kOnlineTriggerFired2]  = "Online trigger fired2";           fgVariableUnits[kOnlineTriggerFired2]  = "";
  fgVariableNames[kIsPhysicsSelection]   = "Physics selection ON";            fgVariableUnits[kIsPhysicsSelection]   = "";
  fgVariableNames[kIsSPDPileup]          = "SPD pileup ON";                   fgVariableUnits[kIsSPDPileup]          = "";
  fgVariableNames[kIsSPDPileup5]          = "SPD pileup (5 contributors) ON";                   fgVariableUnits[kIsSPDPileup5]          = "";
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
  fgVariableNames[kDeltaVtxZspd]            = "#Delta Z (global-SPD)";                        fgVariableUnits[kDeltaVtxZspd]            = "cm";
  for(Int_t iflag=0; iflag<kNTrackingStatus; ++iflag) {
    fgVariableNames[kNTracksPerTrackingStatus+iflag] = Form("Tracks with %s on", fgkTrackingStatusNames[iflag]); 
    fgVariableUnits[kNTracksPerTrackingStatus+iflag] = ""; 
  }
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
  fgVariableNames[kNTPCclusters]                = "No. TPC clusters";                  fgVariableUnits[kNTPCclusters]  = "";
  fgVariableNames[kSPDntracklets]               = "No.SPD tracklets in |#eta| < 1.0";                fgVariableUnits[kSPDntracklets]          = "";
  fgVariableNames[kSPDntracklets08]             = "No.SPD tracklets in |#eta| < 0.8";                fgVariableUnits[kSPDntracklets08]        = "";
  fgVariableNames[kSPDntracklets16]             = "No.SPD tracklets in |#eta| < 1.6";                fgVariableUnits[kSPDntracklets16]        = "";
  fgVariableNames[kSPDntrackletsCorr]           = "Corrected no. SPD tracklets in |#eta| < 1.0";     fgVariableUnits[kSPDntrackletsCorr]      = "";
  fgVariableNames[kSPDntrackletsCorrSmear]      = "Corrected no. SPD tracklets in |#eta| < 1.0, Poisson smeared";    
  fgVariableNames[kSPDntrackletsOuterCorr]      = "Corrected no. SPD tracklets in outer eta region";   fgVariableUnits[kSPDntrackletsOuterCorr]     = "";
  fgVariableNames[kSPDntrackletsOuterCorrSmear]          = "Corrected no.SPD tracklets in outer eta region, Poisson smeared";                fgVariableUnits[kSPDntrackletsOuterCorrSmear]     = ""; fgVariableUnits[kSPDntrackletsCorrSmear]      = "";
  for(Int_t ieta=0;ieta<32;++ieta) {
    fgVariableNames[kSPDntrackletsEta+ieta] = Form("No.SPD tracklets in %.1f<#eta<%.1f", -1.6+0.1*ieta, -1.6+0.1*(ieta+1));
    fgVariableUnits[kSPDntrackletsEta+ieta] = "";
  }
  for(Int_t il=0;il<2;++il) {
    fgVariableNames[kSPDFiredChips+il] = Form("Fired chips in SPD layer %d", il+1); 
    fgVariableUnits[kSPDFiredChips+il] = "";
  }
  for(Int_t il=0;il<6;++il) {
    fgVariableNames[kITSnClusters+il] = Form("No. clusters in ITS layer %d", il+1); 
    fgVariableUnits[kITSnClusters+il] = "";
  }
  fgVariableNames[kSPDnSingleClusters]  = "SPD single clusters";    fgVariableUnits[kSPDnSingleClusters]  = "";  
  fgVariableNames[kEventMixingId]       = "Event mixing id";        fgVariableUnits[kEventMixingId]       = "";  
  fgVariableNames[kVZEROATotalMult]     = "Multiplicity VZERO-A";   fgVariableUnits[kVZEROATotalMult]     = "";
  fgVariableNames[kVZEROCTotalMult]     = "Multiplicity VZERO-C";   fgVariableUnits[kVZEROCTotalMult]     = "";
  fgVariableNames[kVZEROTotalMult]      = "Multiplicity VZERO";     fgVariableUnits[kVZEROTotalMult]      = "";
  fgVariableNames[kVZEROTotalMultCorr]      = "Corrected Multiplicity VZERO";     fgVariableUnits[kVZEROTotalMult]      = "";
  fgVariableNames[kVZEROTotalMultCorrSmear]      = "Corrected Multiplicity VZERO, Poisson smeared";     fgVariableUnits[kVZEROTotalMult]      = "";
  fgVariableNames[kVZEROAemptyChannels] = "VZERO-A empty channels"; fgVariableUnits[kVZEROAemptyChannels] = "";
  fgVariableNames[kVZEROCemptyChannels] = "VZERO-C empty channels"; fgVariableUnits[kVZEROCemptyChannels] = "";
  for(Int_t ich=0;ich<64;++ich) {
    fgVariableNames[kVZEROChannelMult+ich] = Form("Multiplicity VZERO ch.%d", ich);
    fgVariableUnits[kVZEROChannelMult+ich] = "";
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

  fgVariableNames[kMultEstimatorOnlineV0M]    = "Muliplicity estimator V0M";
  fgVariableNames[kMultEstimatorOnlineV0A]    = "Muliplicity estimator V0A";
  fgVariableNames[kMultEstimatorOnlineV0C]    = "Muliplicity estimator V0C";
  fgVariableNames[kMultEstimatorADM]          = "Muliplicity estimator ADM ";
  fgVariableNames[kMultEstimatorADA]          = "Muliplicity estimator ADA";
  fgVariableNames[kMultEstimatorADC]          = "Muliplicity estimator ADC";
  fgVariableNames[kMultEstimatorSPDClusters]  = "Muliplicity estimator SPDclusters";
  fgVariableNames[kMultEstimatorSPDTracklets] = "Muliplicity estimator SPDtracklets";
  fgVariableNames[kMultEstimatorRefMult05]    = "Muliplicity estimator RefMult05";
  fgVariableNames[kMultEstimatorRefMult08]    = "Muliplicity estimator RefMult08";

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


  fgVariableNames[kMultEstimatorPercentileOnlineV0M]    = "Muliplicity estimator percentile V0M";
  fgVariableNames[kMultEstimatorPercentileOnlineV0A]    = "Muliplicity estimator percentile V0A";
  fgVariableNames[kMultEstimatorPercentileOnlineV0C]    = "Muliplicity estimator percentile V0C";
  fgVariableNames[kMultEstimatorPercentileADM]          = "Muliplicity estimator percentile ADM ";
  fgVariableNames[kMultEstimatorPercentileADA]          = "Muliplicity estimator percentile ADA";
  fgVariableNames[kMultEstimatorPercentileADC]          = "Muliplicity estimator percentile ADC";
  fgVariableNames[kMultEstimatorPercentileSPDClusters]  = "Muliplicity estimator percentile SPDclusters";
  fgVariableNames[kMultEstimatorPercentileSPDTracklets] = "Muliplicity estimator percentile SPDtracklets";
  fgVariableNames[kMultEstimatorPercentileRefMult05]    = "Muliplicity estimator percentile RefMult05";
  fgVariableNames[kMultEstimatorPercentileRefMult08]    = "Muliplicity estimator percentile RefMult08";

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
  fgVariableNames[kHighMultV0Triggered] = "event was triggered with HighMultV0"; fgVariableUnits[kHighMultV0Triggered] = "";

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
    fgVariableNames[kTPCRPres+6+iHarmonic] = Form("#sqrt{%d(#Psi_{%d}^{TPC}-#Psi_{%d}^{VZERO-C})}", 
						  iHarmonic+1, iHarmonic+1, iHarmonic+1);
    fgVariableUnits[kTPCRPres+6+iHarmonic] = "rad^{1/2}";
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
  }  // end loop over harmonics 
  
  fgVariableNames[kPt]    = "p_{T}";   fgVariableUnits[kPt]    = "GeV/c";
  fgVariableNames[kPtMC]    = "p^{MC}_{T}";   fgVariableUnits[kPtMC]    = "GeV/c";
  fgVariableNames[kP]     = "p";       fgVariableUnits[kP]     = "GeV/c";
  fgVariableNames[kPMC]     = "p^{MC}";       fgVariableUnits[kPMC]     = "GeV/c";
  fgVariableNames[kPx]    = "p_{x}";   fgVariableUnits[kPx]    = "GeV/c";
  fgVariableNames[kPxMC]    = "p^{MC}_{x}";   fgVariableUnits[kPxMC]    = "GeV/c";
  fgVariableNames[kPy]    = "p_{y}";   fgVariableUnits[kPy]    = "GeV/c";
  fgVariableNames[kPyMC]    = "p^{MC}_{y}";   fgVariableUnits[kPyMC]    = "GeV/c";
  fgVariableNames[kPz]    = "p_{z}";   fgVariableUnits[kPz]    = "GeV/c";
  fgVariableNames[kPzMC]    = "p^{MC}_{z}";   fgVariableUnits[kPzMC]    = "GeV/c";
  fgVariableNames[kTheta] = "#theta";  fgVariableUnits[kTheta] = "rad.";
  fgVariableNames[kThetaMC] = "#theta^{MC}";  fgVariableUnits[kThetaMC] = "rad.";
  fgVariableNames[kEta]   = "#eta";    fgVariableUnits[kEta]   = "";
  fgVariableNames[kEtaMC]   = "#eta^{MC}";    fgVariableUnits[kEtaMC]   = "";
  fgVariableNames[kPhi]   = "#varphi"; fgVariableUnits[kPhi]   = "rad.";
  fgVariableNames[kPhiMC]   = "#varphi^{MC}"; fgVariableUnits[kPhiMC]   = "rad.";
  for(Int_t iHarmonic=0;iHarmonic<6;++iHarmonic) {
    fgVariableNames[kCosNPhi+iHarmonic] = Form("cos(%d#varphi)",iHarmonic+1); fgVariableUnits[kCosNPhi+iHarmonic] = "";
    fgVariableNames[kSinNPhi+iHarmonic] = Form("sin(%d#varphi)",iHarmonic+1); fgVariableUnits[kSinNPhi+iHarmonic] = "";
  }
  fgVariableNames[kPtSquared]     = "p_{T}^{2}";     fgVariableUnits[kPtSquared]     = "GeV^{2}/c^{2}";
  fgVariableNames[kOneOverSqrtPt] = "1/sqrt(p_{T})"; fgVariableUnits[kOneOverSqrtPt] = "GeV^{-1/2}";
  fgVariableNames[kMass]      = "m";         fgVariableUnits[kMass]        = "GeV/c^{2}";
  fgVariableNames[kMassMC]      = "m^{MC}";         fgVariableUnits[kMassMC]        = "GeV/c^{2}";
  fgVariableNames[kRap]       = "y";         fgVariableUnits[kRap]         = "";  
  fgVariableNames[kRapMC]       = "y^{MC}";         fgVariableUnits[kRapMC]         = "";  
  fgVariableNames[kPdgMC]       = "PDG code";          fgVariableUnits[kPdgMC]         = "";
  fgVariableNames[kPdgMC+1]  = "mother's PDG code";     fgVariableUnits[kPdgMC+1]         = "";
  fgVariableNames[kPdgMC+2]  = "grand-mother's PDG code";     fgVariableUnits[kPdgMC+2]         = "";
  fgVariableNames[kPdgMC+3]  = "grand-grand-mother's PDG code";     fgVariableUnits[kPdgMC+3]         = "";
  
  fgVariableNames[kCandidateId]       = "pair id.";              fgVariableUnits[kCandidateId]       = "";
  fgVariableNames[kPairType]          = "pair type";             fgVariableUnits[kPairType]          = "";
  fgVariableNames[kMassV0]            = "m_{K^{0}_{S}}";         fgVariableUnits[kMassV0]            = "GeV/c^{2}";
  fgVariableNames[kMassV0+1]          = "m_{#Lambda^{0}}";       fgVariableUnits[kMassV0+1]          = "GeV/c^{2}";
  fgVariableNames[kMassV0+2]          = "m_{#bar{#Lambda^{0}}}"; fgVariableUnits[kMassV0+2]          = "GeV/c^{2}";
  fgVariableNames[kMassV0+3]          = "m_{#gamma}";            fgVariableUnits[kMassV0+3]          = "GeV/c^{2}";
  fgVariableNames[kPairChisquare]     = "pair #chi^{2}";         fgVariableUnits[kPairChisquare]     = "";
  fgVariableNames[kPairLxy]           = "L_{xy}";                fgVariableUnits[kPairLxy]           = "cm.";
  fgVariableNames[kPairOpeningAngle]  = "pair opening angle";    fgVariableUnits[kPairOpeningAngle]  = "rad.";    
  fgVariableNames[kPairPointingAngle] = "#theta_{pointing}";     fgVariableUnits[kPairPointingAngle] = "rad.";
  fgVariableNames[kPairThetaCS]       = "cos(#theta^{*}_{CS})";  fgVariableUnits[kPairThetaCS]       = "";  
  fgVariableNames[kPairPhiCS]         = "#varphi^{*}_{CS}";      fgVariableUnits[kPairPhiCS]         = "rad.";  
  fgVariableNames[kPairThetaHE]       = "cos(#theta^{*}_{HE})";  fgVariableUnits[kPairThetaHE]       = "";  
  fgVariableNames[kPairPhiHE]         = "#varphi^{*}_{HE}";      fgVariableUnits[kPairPhiHE]         = "rad.";
  fgVariableNames[kPairPhiV]         = "#varphi^{*}_{v}";        fgVariableUnits[kPairPhiV]         = "rad.";
  for(Int_t i=0;i<2;++i) {
     fgVariableNames[kPairLegTPCchi2+i] = Form("TPC #chi^{2}, leg %d", i+1);
     fgVariableUnits[kPairLegTPCchi2+i] = "";
     fgVariableNames[kPairLegITSchi2+i] = Form("ITS #chi^{2}, leg %d", i+1);
     fgVariableUnits[kPairLegITSchi2+i] = "";
  }
  
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
  fgVariableNames[kTPCcrossedRowsOverFindableClusters] = "Crossed rows / findable clusters"; fgVariableUnits[kTPCcrossedRowsOverFindableClusters] = "fraction";
  fgVariableNames[kTPCnclsRatio3]     = "No.TPC crossed rows/findable clusters"; fgVariableUnits[kTPCnclsRatio3] = "";
  fgVariableNames[kTPCsignal]         = "TPC dE/dx";                    fgVariableUnits[kTPCsignal] = "";  
  fgVariableNames[kTPCsignalN]        = "No. TPC clusters PID";         fgVariableUnits[kTPCsignalN] = "";  
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
  fgVariableNames[kEMCALmatchedEnergy]    = "Calo energy";             fgVariableUnits[kEMCALmatchedEnergy] = "GeV";
  fgVariableNames[kEMCALmatchedClusterId] = "matched Calo cluster id"; fgVariableUnits[kEMCALmatchedClusterId] = "";
  fgVariableNames[kEMCALmatchedEOverP]    = "Calo E/p";                fgVariableUnits[kEMCALmatchedEOverP] = "";  
  fgVariableNames[kEMCALclusterEnergy]    = "Calo cls. energy";        fgVariableUnits[kEMCALclusterEnergy] = "GeV";    
  fgVariableNames[kEMCALclusterDx]        = "Calo cls. dx";            fgVariableUnits[kEMCALclusterDx] = "";  
  fgVariableNames[kEMCALclusterDz]        = "Calo cls. dz";            fgVariableUnits[kEMCALclusterDz] = "";  
  fgVariableNames[kEMCALdetector]         = "Calo detector";           fgVariableUnits[kEMCALdetector] = "";  
  fgVariableNames[kEMCALm20]              = "Cluster short axis M20";  fgVariableUnits[kEMCALm20] = "";  
  fgVariableNames[kEMCALm02]              = "Cluster short axis M02";  fgVariableUnits[kEMCALm02] = "";  
  fgVariableNames[kEMCALdispersion]       = "Cluster dispersion";      fgVariableUnits[kEMCALdispersion] = "";  
  fgVariableNames[kTrackingFlag] = "Tracking flag";  fgVariableUnits[kTrackingFlag] = "";  
  fgVariableNames[kTrackingStatus] = "Tracking status";  fgVariableUnits[kTrackingStatus] = "";  
  fgVariableNames[kDeltaPhi]     = "#Delta #varphi"; fgVariableUnits[kDeltaPhi] = "rad.";  
  fgVariableNames[kDeltaTheta]   = "#Delta #theta";  fgVariableUnits[kDeltaTheta] = "rad.";
  fgVariableNames[kDeltaEta]     = "#Delta #eta";    fgVariableUnits[kDeltaEta] = "";
  
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
void AliReducedVarManager::SetRunNumbers( TString runNumbers ){
  TObjArray* runNumbersArr = runNumbers.Tokenize(";");
  runNumbersArr->SetOwner();
  for( Int_t iRun=0; iRun < runNumbersArr->GetEntries(); ++iRun){
    TString runNumberString = runNumbersArr->At(iRun)->GetName();
    fgRunNumbers.push_back( runNumberString.Atoi() );
  }
}

//____________________________________________________________________________________
void AliReducedVarManager::SetMultiplicityProfile(TH1* profile, Int_t estimator) {
   //
   // initialize the profile for the z-vertex equalization of the multiplicity estimator
   //
   if(profile) {
     fgAvgMultVsVertex[estimator] = (TH1*)profile->Clone(Form("AliReducedVarManager_AverageMultiplicityVsVertex"));
     fgAvgMultVsVertex[estimator]->SetDirectory(0x0);
   }
   else{
     cout <<"AliReducedVarManager::SetMultiplicityProfile : Profile empty. " << estimator << endl;
   }
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
