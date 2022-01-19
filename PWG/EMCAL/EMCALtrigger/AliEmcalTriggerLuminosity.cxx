/************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                 *
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
#include <array>
#include <iostream>
#include <memory>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TList.h>
#include "AliEmcalTriggerLuminosity.h"

ClassImp(PWG::EMCAL::AliEmcalTriggerLuminosity)

using namespace PWG::EMCAL;

AliEmcalTriggerLuminosity::AliEmcalTriggerLuminosity():
  TObject(),
  fCollisionType(CollisionType_t::kUnknown),
  fYear(-1),
  fLuminosityHist(nullptr),
  fClusterCounters(),
  fLuminosities()
{

}

AliEmcalTriggerLuminosity::AliEmcalTriggerLuminosity(const char *fname, const char *directory):
  TObject(),
  fCollisionType(),
  fYear(-1),
  fLuminosityHist(nullptr),
  fClusterCounters(),
  fLuminosities()
{
  InitFromFile(fname, directory);
}

AliEmcalTriggerLuminosity::~AliEmcalTriggerLuminosity() {
  delete fLuminosityHist; 
  for(auto clustercounter : fClusterCounters) delete clustercounter.second;
}

void AliEmcalTriggerLuminosity::LoadHistograms(const TList &histlist) {
  fLuminosityHist = static_cast<TH2 *>(histlist.FindObject("hTriggerLuminosity"));
  if(fLuminosityHist) fLuminosityHist->SetDirectory(nullptr);
  for(auto en : histlist) {
    std::string histname = en->GetName();
    if(histname.find("hClusterCounter") != std::string::npos) {
      std::string triggerclass = histname.substr(15);
      auto clustercounter = static_cast<TH1 *>(en);
      clustercounter->SetDirectory(nullptr);
      fClusterCounters[triggerclass] = clustercounter;
    }
  }

  auto yearhist = dynamic_cast<TH1 *>(histlist.FindObject("hYear"));
  if(yearhist) {
    fYear = determineYear(yearhist);
    std::cout << "Auto-detected year " << fYear << std::endl; 
  } else {
    std::cout << "Output does not contain the histogram for years, the year must be set manually" << std::endl;
  }
  auto collhist = dynamic_cast<TH1 *>(histlist.FindObject("hCollisionSystem"));
  if(collhist) {
    fCollisionType = determineCollisionType(collhist);
    std::cout << "Auto-detected collsion system " << getCollisionLabel(fCollisionType) << std::endl; 
  } else { 
    std::cout << "Output does not contain the histogram for collision systems, the collision system must be set manually" << std::endl;
  }
}

void AliEmcalTriggerLuminosity::InitFromFile(const char *filename, const char *dirname) {
  std::unique_ptr<TFile> reader(TFile::Open(filename, "READ"));
  if(!reader || reader->IsZombie()) throw InputDataException(filename, dirname);
  auto dirkey = reader->GetListOfKeys()->FindObject(dirname);
  if(!dirkey) throw InputDataException(filename, dirname);
  if(dirkey->IsA() != TDirectoryFile::Class()) throw InputDataException(filename, dirname);
  reader->cd(dirname);
  auto histlist = static_cast<TList *>(static_cast<TKey *>(gDirectory->GetListOfKeys()->At(0))->ReadObj());
  if(!histlist) throw InputDataException(filename, dirname);
  LoadHistograms(*histlist);
}

void AliEmcalTriggerLuminosity::Evaluate() {
  fLuminosities.clear();
  std::vector<std::string> errorfields;
  if(fYear == 0) errorfields.push_back("year");
  if(fCollisionType == CollisionType_t::kUnknown) errorfields.push_back("collision type");
  if(!fLuminosityHist) errorfields.push_back("lumihist");
  if(!fClusterCounters.size()) errorfields.push_back("clustercounter");
  if(errorfields.size()) throw UninitException(errorfields);
  switch (fCollisionType)
  {
  case CollisionType_t::kPP13TeV:
    evaluatePP13TeV();
    break;

  case CollisionType_t::kPPB8TeV:
    evaluatePPB8TeV();
    break;
  
  default:
    std::cout << "Evaluation for collision system " << getCollisionLabel(fCollisionType) << " not (yet) implemented" << std::endl;
    break;
  }
}

double AliEmcalTriggerLuminosity::GetLuminosityForTrigger(const char *trigger) {
  auto found = fLuminosities.find(trigger);
  if(found == fLuminosities.end()) throw TriggerNotFoundException(trigger);
  return found->second;
}

void AliEmcalTriggerLuminosity::evaluatePP13TeV() {
  TH1 *counterEG1 = fClusterCounters["EG1"];
  std::vector<std::string> errorfields;
  if(!counterEG1) errorfields.push_back("counter EG1");
  if(errorfields.size()) throw UninitException(errorfields);

	double corrCENTNOTRD = getTriggerClusterCounts(counterEG1, "CENTNOTRD") / getTriggerClusterCounts(counterEG1, "CENTNOPMD");

	std::cout << "Correction CENTNOTRD: " << corrCENTNOTRD << std::endl;

	// normalize by reference cross section
  std::map<int, double> xsecsMB {{2016, 58.44}, {2017, 58.10}, {2018, 57.52}};
  auto found = xsecsMB.find(fYear);
  if(found == xsecsMB.end()) throw AmbiguityException("Crosssection");
	double xrefMB =  found->second,
	       barnToPB = 1e9,
	       xrefPB = xrefMB * barnToPB;
  for(int ib = 0; ib < fLuminosityHist->GetXaxis()->GetNbins(); ib++) {
    std::string triggerclass = fLuminosityHist->GetXaxis()->GetBinLabel(ib+1);
    std::unique_ptr<TH1> sliceTrigger(fLuminosityHist->ProjectionY("sliceTrigger", ib+1, ib+1));
    double rawcounts = sliceTrigger->Integral();
    double clustercorrection = 1.;
    if(triggerclass.find("G1") != std::string::npos || triggerclass.find("J1") != std::string::npos) {
      clustercorrection = corrCENTNOTRD;
    }
    double correctedCounts = rawcounts * clustercorrection;
    double luminosity = correctedCounts / xrefPB;
    std::cout << "Trigger class " << triggerclass << ": raw counts " << rawcounts << ", corrected counts " << correctedCounts << ", luminosity " << luminosity << " pb-1" << std::endl; 
    fLuminosities[triggerclass] = luminosity;
  }
}

void AliEmcalTriggerLuminosity::evaluatePPB8TeV() {
  TH1 *counterEMC7 = fClusterCounters["EMC7"],
      *counterEG1 = fClusterCounters["EG1"];
  std::vector<std::string> errorfields;
  if(!counterEMC7) errorfields.push_back("counter EMC7");
  if(!counterEG1) errorfields.push_back("counter EG1");
  if(errorfields.size()) throw UninitException(errorfields);

	double corrCENTNOPMD = getTriggerClusterCounts(counterEMC7, "CENTNOPMD") / getTriggerClusterCounts(counterEMC7, "CENT");
	double corrCENTNOTRD = getTriggerClusterCounts(counterEG1, "CENTNOTRD") / getTriggerClusterCounts(counterEG1, "CENTNOPMD");
  double corrCombined = corrCENTNOPMD * corrCENTNOTRD;

	std::cout << "Correction CENTNOPMD: " << corrCENTNOPMD << std::endl;
	std::cout << "Correction CENTNOTRD: " << corrCENTNOTRD << std::endl;
  std::cout << "Combined Correction:  " << corrCombined << std::endl;

	// normalize by reference cross section
	double xrefBarn = 2.09,
	       barnToNB = 1e9,
	       xrefNB = xrefBarn * barnToNB;
  for(int ib = 0; ib < fLuminosityHist->GetXaxis()->GetNbins(); ib++) {
    std::string triggerclass = fLuminosityHist->GetXaxis()->GetBinLabel(ib+1);
    std::unique_ptr<TH1> sliceTrigger(fLuminosityHist->ProjectionY("sliceTrigger", ib+1, ib+1));
    double rawcounts = sliceTrigger->Integral();
    double clustercorrection = 1.;
    if(triggerclass.find("MC7") != std::string::npos || triggerclass.find("G2") != std::string::npos || triggerclass.find("J2") != std::string::npos) {
      clustercorrection = corrCENTNOPMD;
    } else if(triggerclass.find("G1") != std::string::npos || triggerclass.find("J1") != std::string::npos) {
      clustercorrection = corrCombined;
    }
    double correctedCounts = rawcounts * clustercorrection;
    double luminosity = correctedCounts / xrefNB;
    std::cout << "Trigger class " << triggerclass << ": raw counts " << rawcounts << ", corrected counts " << correctedCounts << ", luminosity " << luminosity << " nb-1" << std::endl; 
    fLuminosities[triggerclass] = luminosity;
  }
}

double AliEmcalTriggerLuminosity::getTriggerClusterCounts(TH1 * clustercounter, const std::string &clustername) const {
  int binID = clustercounter->GetXaxis()->FindBin(clustername.data());
  if(binID <= 0 || binID > clustercounter->GetXaxis()->GetNbins()) return 0;
  return clustercounter->GetBinContent(binID);
}

int AliEmcalTriggerLuminosity::determineYear(const TH1 * const yearhist) const {
  std::vector<int> binsNonZero;
  for(int ib = 0; ib < yearhist->GetXaxis()->GetNbins(); ib++) {
    if(yearhist->GetBinContent(ib+1)) {
      binsNonZero.push_back(ib);
    }
  }
  if(binsNonZero.size() > 1) throw AmbiguityException("year");
  return static_cast<int>(yearhist->GetXaxis()->GetBinCenter(binsNonZero[0]+1));
}

AliEmcalTriggerLuminosity::CollisionType_t AliEmcalTriggerLuminosity::determineCollisionType(const TH1 *const collisionhist) const {
  std::vector<int> binsNonZero;
  for(int ib = 0; ib < collisionhist->GetXaxis()->GetNbins(); ib++) {
    if(collisionhist->GetBinContent(ib+1)) {
      binsNonZero.push_back(ib);
    }
  }
  if(binsNonZero.size() > 1) throw AmbiguityException("year");
  std::string bintitle = collisionhist->GetXaxis()->GetBinLabel(binsNonZero[0]+1);  
  return getCollisionType(bintitle);
}

std::string AliEmcalTriggerLuminosity::getCollisionLabel(AliEmcalTriggerLuminosity::CollisionType_t coltype) const {
  switch(coltype) {
    case CollisionType_t::kPP13TeV: return "pp, 13 TeV";
    case CollisionType_t::kPP5TeV: return "pp, 5.02 TeV"; 
    case CollisionType_t::kPPB5TeV: return "p-Pb, 5.02 TeV";
    case CollisionType_t::kPPB8TeV: return "p-Pb, 8.16 TeV";
    case CollisionType_t::kPBPB5TeV: return "Pb-Pb, 5.02 TeV";
    default: return "Unknown collision type";
  }
}

AliEmcalTriggerLuminosity::CollisionType_t AliEmcalTriggerLuminosity::getCollisionType(const std::string &collabel) const {
  std::map<std::string, CollisionType_t> collisionsystems = {{"pp, 13 TeV", CollisionType_t::kPP13TeV}, 
                                                             {"pp, 5.02 TeV", CollisionType_t::kPP5TeV}, 
                                                             {"p-Pb, 5.02 TeV", CollisionType_t::kPPB5TeV}, 
                                                             {"p-Pb, 8.16 TeV", CollisionType_t::kPPB8TeV}, 
                                                             {"Pb-Pb, 5.02 TeV", CollisionType_t::kPBPB5TeV}};
  auto found = collisionsystems.find(collabel);
  if(found == collisionsystems.end()) return CollisionType_t::kUnknown;
  return found->second;
}

void AliEmcalTriggerLuminosity::UninitException::buildErrorMessage(){
  std::stringstream msgbuilder;
  msgbuilder << "Fields uninitialized: ";
  bool first = true;
  // Thanks ALICE offline for forcing us to go the hard way
  for(const auto &field : fFields){
    if(first) {
      first = false;
    } else {
      msgbuilder << ", ";
    }
    msgbuilder << field;
  }
  fMessage = msgbuilder.str();

}