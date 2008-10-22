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

//====================================================================
ClassImp(AliFMDAnaParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

const char* AliFMDAnaParameters::fgkBackgroundCorrection  = "FMD/AnalysisCalib/Background";
const char* AliFMDAnaParameters::fgkEnergyDists  = "FMD/AnalysisCalib/EnergyDistribution";
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
  fBackgroundArray(0),
  fEdistArray(0)
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
  
  fBackgroundArray = dynamic_cast<TObjArray*>(background->GetObject());
  if (!fBackgroundArray) AliFatal("Invalid background object from CDB");
  
}
//____________________________________________________________________

void AliFMDAnaParameters::InitEnergyDists() {
  
  AliCDBEntry*   edist = GetEntry(fgkEnergyDists);
  if (!edist) return;
  
  fEdistArray = dynamic_cast<TObjArray*>(edist->GetObject());
  
  if (!fEdistArray) AliFatal("Invalid background object from CDB");
  
}
//____________________________________________________________________
Float_t AliFMDAnaParameters::GetVtxCutZ() {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  TAxis* refAxis = GetRefAxis();
  
  return refAxis->GetXmax();
}

//____________________________________________________________________
Int_t AliFMDAnaParameters::GetNvtxBins() {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return -1;
  }
  
  TAxis* refAxis = GetRefAxis();
  
  return refAxis->GetNbins();
}
//____________________________________________________________________
TH1F* AliFMDAnaParameters::GetEnergyDistribution(Int_t det, Char_t ring) {
  
  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  
  TObjArray* detArray   = (TObjArray*)fEdistArray->At(det);
  Int_t ringNumber      = (ring == 'I' ? 0 : 1);
  TH1F* hEnergyDist     = (TH1F*)detArray->At(ringNumber);  
  return hEnergyDist;
}
//____________________________________________________________________
TH2F* AliFMDAnaParameters::GetBackgroundCorrection(Int_t det, 
						   Char_t ring, 
						   Int_t vtxbin) {

  if(!fIsInit) {
    AliWarning("Not initialized yet. Call Init() to remedy");
    return 0;
  }
  if(vtxbin > GetNvtxBins()) {
    AliWarning(Form("No background object for vertex bin %d", vtxbin));
    return 0;
  } 
  
  TObjArray* correction = GetBackgroundArray();
  
  TObjArray* detArray   = (TObjArray*)correction->At(det);
  Int_t ringNumber      = (ring == 'I' ? 0 : 1);
  TObjArray* ringArray  = (TObjArray*)detArray->At(ringNumber);
  TH2F* bgHist          = (TH2F*)ringArray->At(vtxbin);
  
  return bgHist;
}
//____________________________________________________________________
TAxis* AliFMDAnaParameters::GetRefAxis() {
  
  TAxis* refAxis = (TAxis*)fBackgroundArray->At(0);
  
  return refAxis;
}
//____________________________________________________________________
TObjArray* AliFMDAnaParameters::GetBackgroundArray() {
  
  TObjArray* correction = (TObjArray*)fBackgroundArray->At(1);
  
  return correction;
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
