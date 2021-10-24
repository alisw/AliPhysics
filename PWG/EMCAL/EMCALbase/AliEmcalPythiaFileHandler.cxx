/************************************************************************************
 * Copyright (C) 2021, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include "AliEmcalPythiaFileHandler.h"
#include "AliLog.h"
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <memory>
#include <iostream>

using namespace PWG::EMCAL;

AliEmcalPythiaFileHandler *AliEmcalPythiaFileHandler::fgInstance = nullptr;

AliEmcalPythiaFileHandler::AliEmcalPythiaFileHandler():
  fCrossSectionNrials(),
  fCurrentFile(),
  fInitialized(false)
{
}

AliEmcalPythiaFileHandler::~AliEmcalPythiaFileHandler() {
  if(fgInstance) delete fgInstance;
}

AliEmcalPythiaFileHandler *AliEmcalPythiaFileHandler::Instance() {
  if(!fgInstance) fgInstance = new AliEmcalPythiaFileHandler;
  return fgInstance;
}

AliEmcalPythiaCrossSectionData AliEmcalPythiaFileHandler::GetCrossSectionAndNTrials(const char *filename) {
  if(fCurrentFile != std::string(filename)) {
    UpdateCache(filename);
    fCurrentFile = filename;
  } 
  AliDebugGeneralStream("AliEmcalPythiaFileHandler::GetCrossSectionAndNTrials", 2) << "Returning cached value for cross section and number of trials" << std::endl;
  return fCrossSectionNrials;
}

void AliEmcalPythiaFileHandler::UpdateCache(const char *filename){
  AliDebugGeneralStream("AliEmcalPythiaFileHandler::GetCrossSectionAndNTrials", 1) << "New file - updating cache" << std::endl;
  fCrossSectionNrials.fCrossSection = 0;
  fCrossSectionNrials.fNTrials = 1;
  fInitialized = false;

  TString file(filename);
  bool isAOD = file.Contains("AliAOD.root"),
       isESD = file.Contains("AliESDs.root");
  // Determine archive type
  TString archivetype;
  std::unique_ptr<TObjArray> walk(file.Tokenize("/"));
  for(auto t : *walk){
    TString &tok = static_cast<TObjString *>(t)->String();
    if(tok.Contains(".zip")){
      archivetype = tok;
      Int_t pos = archivetype.Index(".zip");
      archivetype.Replace(pos, archivetype.Length() - pos, "");
    }
  }
  if(archivetype.Length()){
    Ssiz_t pos1 = file.Index(archivetype,archivetype.Length(),0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }

  if(isAOD) {
    UpdateFromXsecHistFile(Form("%s%s",file.Data(),"pyxsec_hists.root"));
  } else if(isESD){
    UpdateFromXsecFile(Form("%s%s",file.Data(),"pyxsec.root"));
  } else {
    throw FileNotFoundException(file);
  }
  fInitialized = true;
}

void AliEmcalPythiaFileHandler::UpdateFromXsecFile(const char *pyxsecfile) {
  AliDebugGeneralStream("AliEmcalPythiaFileHandler::UpdateFromXsecFile", 1) << "Reading cross section from tree file " << pyxsecfile << " (ESD mode)" << std::endl;
  std::unique_ptr<TFile> xsecfile(TFile::Open(pyxsecfile));
  if(!xsecfile || xsecfile->IsZombie()) throw 1;
  auto xtree = (TTree*)xsecfile->Get("Xsection");
  if (!xtree) throw 2;
  UInt_t   ntrials  = 0;
  Double_t  xsection  = 0;
  xtree->SetBranchAddress("xsection",&xsection);
  xtree->SetBranchAddress("ntrials",&ntrials);
  xtree->GetEntry(0);
  fCrossSectionNrials.fNTrials = ntrials;
  fCrossSectionNrials.fCrossSection = xsection;
}

void AliEmcalPythiaFileHandler::UpdateFromXsecHistFile(const char *pyxsechistfile) {
  AliDebugGeneralStream("AliEmcalPythiaFileHandler::UpdateFromXsecHistFile", 1) << "Reading cross section from hist file " << pyxsechistfile << " (AOD mode)" << std::endl;
  std::unique_ptr<TFile> xsecfile(TFile::Open(pyxsechistfile));
  if(!xsecfile || xsecfile->IsZombie()) throw FileNotFoundException(pyxsechistfile);
  // find the tlist we want to be independtent of the name so use the Tkey
  auto key = static_cast<TKey *>(xsecfile->GetListOfKeys()->At(0)); 
  if (!key) throw FileContentException(pyxsechistfile);
  auto list = dynamic_cast<TList *>(key->ReadObj());
  if (!list) throw FileContentException(pyxsechistfile);
  auto xSecHist = static_cast<TProfile*>(list->FindObject("h1Xsec"));
  // check for failure
  if(xSecHist->GetEntries()) fCrossSectionNrials.fCrossSection = xSecHist->GetBinContent(1);
  fCrossSectionNrials.fNTrials  = static_cast<TH1 *>(list->FindObject("h1Trials"))->GetBinContent(1);
}
