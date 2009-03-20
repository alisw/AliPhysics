/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//The design of this class is based on the AliFMDParameters class. Its purpose
//is to hold parameters for the analysis such as background correction and 
//fit functions.
//
//Author: Hans Hjersing Dalsgaard, NBI, hans.dalsgaard@cern.ch
//

#include "AliFMDDebug.h"		   // ALILOG_H
#include "AliFMDAnaParameters.h"	   // ALIFMDPARAMETERS_H
#include <AliCDBManager.h>         // ALICDBMANAGER_H
#include <AliCDBEntry.h>           // ALICDBMANAGER_H
#include <AliLog.h>
#include <Riostream.h>
#include <sstream>
#include <TSystem.h>
#include <TH2D.h>
#include <TF1.h>

//====================================================================
ClassImp(AliFMDAnaParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

const char* AliFMDAnaParameters::fgkBackgroundCorrection  = "FMD/Correction/Background";
const char* AliFMDAnaParameters::fgkEnergyDists  = "FMD/Correction/EnergyDistribution";
//____________________________________________________________________
AliFMDAnaParameters* AliFMDAnaParameters::fgInstance = 0;

//____________________________________________________________________

AliFMDAnaParameters* 
AliFMDAnaParameters::Instance() 
{
  // Get static instance 
  if (!fgInstance) fgInstance = new AliFMDAnaParameters;
  return fgInstance;
}

//____________________________________________________________________
AliFMDAnaParameters::AliFMDAnaParameters() :
  fIsInit(kFALSE),
  fBackground(0),
  fEnergyDistribution(0)
{
  
  // Default constructor 
}
//____________________________________________________________________
void AliFMDAnaParameters::Init(Bool_t forceReInit, UInt_t what)
{
  // Initialize the parameters manager.  We need to get stuff from the
  // CDB here. 
  if (forceReInit) fIsInit = kFALSE;
  if (fIsInit) return;
  if (what & kBackgroundCorrection)       InitBackground();
  if (what & kEnergyDistributions)        InitEnergyDists();
  
  fIsInit = kTRUE;
}
//____________________________________________________________________

void AliFMDAnaParameters::InitBackground() {
  
  AliCDBEntry*   background = GetEntry(fgkBackgroundCorrection);
  if (!background) return;
  
  fBackground = dynamic_cast<AliFMDAnaCalibBackgroundCorrection*>(background->GetObject());
  if (!fBackground) AliFatal("Invalid background object from CDB");
  
}
//____________________________________________________________________

void AliFMDAnaParameters::InitEnergyDists() {
  
  AliCDBEntry*   edist = GetEntry(fgkEnergyDists);
  if (!edist) return;
  
  fEnergyDistribution = dynamic_cast<AliFMDAnaCalibEnergyDistribution*>(edist->GetObject());
  
  if (!fEnergyDistribution) AliFatal("Invalid background object from CDB");
  
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetVtxCutZ() {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  return fBackground->GetVtxCutZ();
}

//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNvtxBins() {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  return fBackground->GetNvtxBins();
}
//____________________________________________________________________
TH1F* AliFMDAnaParameters::GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta) {
  
  return fEnergyDistribution->GetEnergyDistribution(det, ring, eta);
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetSigma(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist       = GetEnergyDistribution(det,ring, eta);
  TF1*  fitFunc           = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) {
    AliWarning(Form("No function for FMD%d%c, eta %f",det,ring,eta));
    return 1024;
  }
  Float_t sigma           = fitFunc->GetParameter(2);
  return sigma;
}


//____________________________________________________________________
Float_t AliFMDAnaParameters::GetMPV(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) {
    AliWarning(Form("No function for FMD%d%c, eta %f",det,ring,eta));
    return 1024;
  }
    
  Float_t MPV           = fitFunc->GetParameter(1);
  return MPV;
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::Get2MIPWeight(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) return 0;
  Float_t TwoMIPweight    = fitFunc->GetParameter(3);
  
  
  
  if(TwoMIPweight < 1e-05)
    TwoMIPweight = 0;
  
  return TwoMIPweight;
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::Get3MIPWeight(Int_t det, Char_t ring, Float_t eta) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TH1F* hEnergyDist     = GetEnergyDistribution(det,ring,eta);
  TF1*  fitFunc         = hEnergyDist->GetFunction("FMDfitFunc");
  if(!fitFunc) return 0;
  Float_t ThreeMIPweight    = fitFunc->GetParameter(4);
  
  if(ThreeMIPweight < 1e-05)
    ThreeMIPweight = 0;
  
  Float_t TwoMIPweight    = fitFunc->GetParameter(3);
  
  if(TwoMIPweight < 1e-05)
    ThreeMIPweight = 0;
    
  return ThreeMIPweight;
}
//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNetaBins() {
  return GetBackgroundCorrection(1,'I',0)->GetNbinsX();
  
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaMin() {

  return GetBackgroundCorrection(1,'I',0)->GetXaxis()->GetXmin();
} 
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetEtaMax() {

return GetBackgroundCorrection(1,'I',0)->GetXaxis()->GetXmax();

}
//____________________________________________________________________

TH2F* AliFMDAnaParameters::GetBackgroundCorrection(Int_t det, 
						   Char_t ring, 
						   Int_t vtxbin) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  
  
  if(vtxbin > fBackground->GetNvtxBins()) {
    AliWarning(Form("No background object for vertex bin %d", vtxbin));
    return 0;
  } 
  
  return fBackground->GetBgCorrection(det,ring,vtxbin);
}

//____________________________________________________________________
AliCDBEntry* AliFMDAnaParameters::GetEntry(const char* path, Bool_t fatal) const
{
  // Get an entry from the CDB or via preprocessor 
  AliCDBEntry* entry = 0;
  AliCDBManager* cdb = AliCDBManager::Instance();
  entry              = cdb->Get(path);
  
  if (!entry) { 
    TString msg(Form("No %s found in CDB, perhaps you need to "
		     "use AliFMDCalibFaker?", path));
    if (fatal) { AliFatal(msg.Data()); }
    else       AliLog::Message(AliLog::kWarning, msg.Data(), "FMD", 
			       "AliFMDParameters", "GetEntry", __FILE__, 
			       __LINE__);
    return 0;
  }
  return entry;
}


//____________________________________________________________________
//
// EOF
//
