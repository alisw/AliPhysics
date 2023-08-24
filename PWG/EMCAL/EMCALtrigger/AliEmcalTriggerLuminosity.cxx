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
  fLuminosities(),
  fCountsCorr(),
  fCorrCENTNOTRD(0)
{

}

AliEmcalTriggerLuminosity::AliEmcalTriggerLuminosity(const char *fname, const char *directory):
  TObject(),
  fCollisionType(),
  fYear(-1),
  fLuminosityHist(nullptr),
  fClusterCounters(),
  fLuminosities(),
  fCountsCorr(),
  fCorrCENTNOTRD(0)
{
  InitFromFile(fname, directory);
}

AliEmcalTriggerLuminosity::AliEmcalTriggerLuminosity(TList *luminosityHistograms):
  TObject(),
  fCollisionType(),
  fYear(-1),
  fLuminosityHist(nullptr),
  fClusterCounters(),
  fLuminosities(),
  fCountsCorr(),
  fCorrCENTNOTRD(0)
{
  InitFromList(luminosityHistograms);
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
    } else if(histname.find("hTriggerCorrelation") != std::string::npos) {
      std::string triggercluster = histname.substr(19);
      auto corrhist = static_cast<TH2 *>(en);
      corrhist->SetDirectory(nullptr);
      fCorrelationHistograms[triggercluster] = corrhist;
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

void AliEmcalTriggerLuminosity::InitFromList(TList *luminosityHistograms) {
  LoadHistograms(*luminosityHistograms);
}

void AliEmcalTriggerLuminosity::Evaluate() {
  fLuminosities.clear();
  fCountsCorr.clear();
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

    case CollisionType_t::kPP5TeV:
      evaluatePP5TeV();
      break;

    case CollisionType_t::kPBPB5TeV:
      evaluatePBPB5TeV();
      break;

  case CollisionType_t::kPPB8TeV:
    evaluatePPB8TeV();
    break;
  
  default:
    std::cout << "Evaluation for collision system " << getCollisionLabel(fCollisionType) << " not (yet) implemented" << std::endl;
    break;
  }
}

double AliEmcalTriggerLuminosity::GetLuminosityForTrigger(const char *trigger, LuminosityUnit_t unit) const {
  auto found = fLuminosities.find(trigger);
  if(found == fLuminosities.end()) throw TriggerNotFoundException(trigger);
  return convertLuminosity(found->second, unit);
}

double AliEmcalTriggerLuminosity::GetCorrectedCountsForTrigger(const char *trigger) const {
  auto found = fCountsCorr.find(trigger);
  if(found == fCountsCorr.end()) throw TriggerNotFoundException(trigger);
  return found->second;
}

double AliEmcalTriggerLuminosity::GetEffectiveLuminosityForTrigger(const char *trigger, LuminosityUnit_t unit) const {
  auto found = fEffectiveLuminosity.find(trigger);
  if(found == fEffectiveLuminosity.end()) throw TriggerNotFoundException(trigger);
  return convertLuminosity(found->second, unit);
}

double AliEmcalTriggerLuminosity::GetEffectiveDownscalingForTrigger(const char *trigger) const {
  auto found = fEffectiveDownscaling.find(trigger);
  if(found == fEffectiveDownscaling.end()) throw TriggerNotFoundException(trigger);
  return found->second;
}

double AliEmcalTriggerLuminosity::GetLuminosityUncertaintyForTrigger(const char *trigger) const {
  auto foundExpected = fLuminosities.find(trigger);
  auto foundEffective = fEffectiveLuminosity.find(trigger);
  if(foundExpected == fLuminosities.end() || foundEffective == fEffectiveLuminosity.end()) throw TriggerNotFoundException(trigger);
  auto expectedLumi = foundExpected->second,
       effectiveLumi = foundEffective->second;
  return TMath::Abs(expectedLumi - effectiveLumi)/expectedLumi;
}

void AliEmcalTriggerLuminosity::evaluatePP13TeV() {
  TH1 *counterEG1 = fClusterCounters["EG1"];
  std::vector<std::string> errorfields;
  if(!counterEG1) errorfields.push_back("counter EG1");
  if(errorfields.size()) throw UninitException(errorfields);

	double corrCENTNOTRD = getTriggerClusterCounts(counterEG1, "CENTNOTRD") / getTriggerClusterCounts(counterEG1, "CENT");
  fCorrCENTNOTRD =  corrCENTNOTRD;

	std::cout << "Correction CENTNOTRD: " << corrCENTNOTRD << std::endl;

	// normalize by reference cross section
  std::map<int, double> xsecsMB {{2016, 58.44}, {2017, 58.10}, {2018, 57.52}}; // From https://twiki.cern.ch/twiki/pub/ALICE/EventNormalization/LumiUncCombinationExample.pdf
  auto found = xsecsMB.find(fYear);
  if(found == xsecsMB.end()) throw AmbiguityException("Crosssection");
	double xrefMB =  found->second,
	       xrefPB = xrefMB * getConversionToPB(LuminosityUnit_t::kMb); // convert cross section from mb to pb
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
    fCountsCorr  [triggerclass] = correctedCounts;
  }

  // extract effective downscalings and effective luminosities based on the trigger correlation of the CENT cluster
  auto foundCorrelationhist = fCorrelationHistograms.find("CENT");
  if(foundCorrelationhist != fCorrelationHistograms.end()) {
    auto corrhist = foundCorrelationhist->second;
    fEffectiveDownscaling["INT7"] = getEffectiveDownscaling(corrhist, "INT7", "EG1");
    fEffectiveDownscaling["EMC7"] = getEffectiveDownscaling(corrhist, "EMC7", "EG1");
    fEffectiveDownscaling["EG1"] = getEffectiveDownscaling(corrhist, "EG1", "EG1");
    fEffectiveDownscaling["EG2"] = getEffectiveDownscaling(corrhist, "EG2", "EG1");
    fEffectiveDownscaling["EJ1"] = getEffectiveDownscaling(corrhist, "EJ1", "EJ1");
    fEffectiveDownscaling["EJ2"] = getEffectiveDownscaling(corrhist, "EJ2", "EJ1");
    fEffectiveDownscaling["DMC7"] = getEffectiveDownscaling(corrhist, "DMC7", "DG1");
    fEffectiveDownscaling["DG1"] = getEffectiveDownscaling(corrhist, "DG1", "DG1");
    fEffectiveDownscaling["DG2"] = getEffectiveDownscaling(corrhist, "DG2", "DG1");
    fEffectiveDownscaling["DJ1"] = getEffectiveDownscaling(corrhist, "DJ1", "DJ1");
    fEffectiveDownscaling["DJ2"] = getEffectiveDownscaling(corrhist, "DJ2", "DJ2");

    TH1 *counterINT7 = fClusterCounters["INT7"];
    auto countsINT7 = getTriggerClusterCounts(counterINT7, "CENT"),
         refLuminosity = countsINT7 / fEffectiveDownscaling["INT7"];
    fEffectiveLuminosity["INT7"] = countsINT7 / xrefPB;
    for(auto &trgdownscale : fEffectiveDownscaling) {
      auto &triggerclass = trgdownscale.first;
      if(triggerclass == "INT7") continue;
      auto &downscaling = trgdownscale.second;
      double clustercorrection = 1.;
      if(triggerclass.find("G1") != std::string::npos || triggerclass.find("J1") != std::string::npos) {
        clustercorrection = corrCENTNOTRD;
      }
      fEffectiveLuminosity[triggerclass] = refLuminosity * downscaling * clustercorrection / xrefPB;
    }
  }
}

void AliEmcalTriggerLuminosity::evaluatePP5TeV() {
  TH1 *counterEG2 = fClusterCounters["EG2"];
  std::vector<std::string> errorfields;
  if(!counterEG2) errorfields.push_back("counter EG2");
  if(errorfields.size()) throw UninitException(errorfields);

  double corrCENTNOTRD = getTriggerClusterCounts(counterEG2, "CENTNOTRD") / getTriggerClusterCounts(counterEG2, "CENT");
  fCorrCENTNOTRD =  corrCENTNOTRD;

  std::cout << "Correction CENTNOTRD: " << corrCENTNOTRD << std::endl;

  // normalize by reference cross section
  //double xrefMB = 51.2, // From: https://cds.cern.ch/record/2202638/files/vdmNote.pdf, LHC15n
  double xrefMB = 50.87, // V0, From: http://cds.cern.ch/record/2648933/files/vdmNote.pdf, LHC17p+q, unc 2.1%
         xrefPB = xrefMB * getConversionToPB(LuminosityUnit_t::kMb); // convert cross section from mb to pb
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
    fCountsCorr  [triggerclass] = correctedCounts;
  }

  // extract effective downscalings and effective luminosities based on the trigger correlation of the CENT cluster
  auto foundCorrelationhist = fCorrelationHistograms.find("CENT");
  if(foundCorrelationhist != fCorrelationHistograms.end()) {
    auto corrhist = foundCorrelationhist->second;
    fEffectiveDownscaling["INT7"] = getEffectiveDownscaling(corrhist, "INT7", "EG2");
    fEffectiveDownscaling["EMC7"] = getEffectiveDownscaling(corrhist, "EMC7", "EG2");
    fEffectiveDownscaling["EG1"] = getEffectiveDownscaling(corrhist, "EG1", "EG2");
    fEffectiveDownscaling["EG2"] = getEffectiveDownscaling(corrhist, "EG2", "EG2");
    fEffectiveDownscaling["EJ1"] = getEffectiveDownscaling(corrhist, "EJ1", "EJ2");
    fEffectiveDownscaling["EJ2"] = getEffectiveDownscaling(corrhist, "EJ2", "EJ2");
    fEffectiveDownscaling["DMC7"] = getEffectiveDownscaling(corrhist, "DMC7", "DG2");
    fEffectiveDownscaling["DG1"] = getEffectiveDownscaling(corrhist, "DG1", "DG2");
    fEffectiveDownscaling["DG2"] = getEffectiveDownscaling(corrhist, "DG2", "DG2");
    fEffectiveDownscaling["DJ1"] = getEffectiveDownscaling(corrhist, "DJ1", "DJ2");
    fEffectiveDownscaling["DJ2"] = getEffectiveDownscaling(corrhist, "DJ2", "DJ2");

    TH1 *counterINT7 = fClusterCounters["INT7"];
    auto countsINT7 = getTriggerClusterCounts(counterINT7, "CENT"),
         refLuminosity = countsINT7 / fEffectiveDownscaling["INT7"];
    fEffectiveLuminosity["INT7"] = countsINT7 / xrefPB;
    for(auto &trgdownscale : fEffectiveDownscaling) {
      auto &triggerclass = trgdownscale.first;
      if(triggerclass == "INT7") continue;
      auto &downscaling = trgdownscale.second;
      double clustercorrection = 1.;
      if(triggerclass.find("G2") != std::string::npos || triggerclass.find("J2") != std::string::npos) {
        clustercorrection = corrCENTNOTRD;
      }
      fEffectiveLuminosity[triggerclass] = refLuminosity * downscaling * clustercorrection / xrefPB;
    }
  }
}

void AliEmcalTriggerLuminosity::evaluatePBPB5TeV() {
  TH1 *counterEG1 = fClusterCounters["EG1"];
  std::vector<std::string> errorfields;
  if(!counterEG1) errorfields.push_back("counter EG1");
  if(errorfields.size()) throw UninitException(errorfields);

  double corrCENTNOTRD = getTriggerClusterCounts(counterEG1, "CENTNOTRD") / getTriggerClusterCounts(counterEG1, "CENT");
  fCorrCENTNOTRD =  corrCENTNOTRD;

  std::cout << "Correction CENTNOTRD: " << corrCENTNOTRD << std::endl;

  // normalize by reference cross section
  double xrefMB =  67.6, //From: https://cds.cern.ch/record/2636623/files/centrality%20determination%20note.pdf
         xrefPB = xrefMB * getConversionToPB(LuminosityUnit_t::kMb); // convert cross section from mb to pb
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
    fCountsCorr  [triggerclass] = correctedCounts;
  }

  // extract effective downscalings and effective luminosities based on the trigger correlation of the CENT cluster
  auto foundCorrelationhist = fCorrelationHistograms.find("CENT");
  if(foundCorrelationhist != fCorrelationHistograms.end()) {
    auto corrhist = foundCorrelationhist->second;
    fEffectiveDownscaling["INT7"] = getEffectiveDownscaling(corrhist, "INT7", "EG1");
    fEffectiveDownscaling["EMC7"] = getEffectiveDownscaling(corrhist, "EMC7", "EG1");
    fEffectiveDownscaling["EG1"] = getEffectiveDownscaling(corrhist, "EG1", "EG1");
    fEffectiveDownscaling["EG2"] = getEffectiveDownscaling(corrhist, "EG2", "EG1");
    fEffectiveDownscaling["EJ1"] = getEffectiveDownscaling(corrhist, "EJ1", "EJ1");
    fEffectiveDownscaling["EJ2"] = getEffectiveDownscaling(corrhist, "EJ2", "EJ1");
    fEffectiveDownscaling["DMC7"] = getEffectiveDownscaling(corrhist, "DMC7", "DG1");
    fEffectiveDownscaling["DG1"] = getEffectiveDownscaling(corrhist, "DG1", "DG1");
    fEffectiveDownscaling["DG2"] = getEffectiveDownscaling(corrhist, "DG2", "DG1");
    fEffectiveDownscaling["DJ1"] = getEffectiveDownscaling(corrhist, "DJ1", "DJ1");
    fEffectiveDownscaling["DJ2"] = getEffectiveDownscaling(corrhist, "DJ2", "DJ2");

    TH1 *counterINT7 = fClusterCounters["INT7"];
    auto countsINT7 = getTriggerClusterCounts(counterINT7, "CENT"),
         refLuminosity = countsINT7 / fEffectiveDownscaling["INT7"];
    fEffectiveLuminosity["INT7"] = countsINT7 / xrefPB;
    for(auto &trgdownscale : fEffectiveDownscaling) {
      auto &triggerclass = trgdownscale.first;
      if(triggerclass == "INT7") continue;
      auto &downscaling = trgdownscale.second;
      double clustercorrection = 1.;
      if(triggerclass.find("G1") != std::string::npos || triggerclass.find("J1") != std::string::npos) {
        clustercorrection = corrCENTNOTRD;
      }
      fEffectiveLuminosity[triggerclass] = refLuminosity * downscaling * clustercorrection / xrefPB;
    }
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
  fCorrCENTNOTRD =  corrCENTNOTRD;

	std::cout << "Correction CENTNOPMD: " << corrCENTNOPMD << std::endl;
	std::cout << "Correction CENTNOTRD: " << corrCENTNOTRD << std::endl;
  std::cout << "Combined Correction:  " << corrCombined << std::endl;

	// normalize by reference cross section
	double xrefBarn = 2.09,
	       xrefPB = xrefBarn * getConversionToPB(LuminosityUnit_t::kB); // convert cross section from mb to pb
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
    double luminosity = correctedCounts / xrefPB;
    std::cout << "Trigger class " << triggerclass << ": raw counts " << rawcounts << ", corrected counts " << correctedCounts << ", luminosity " << luminosity << " pb-1" << std::endl; 
    fLuminosities[triggerclass] = luminosity;
    fCountsCorr  [triggerclass] = correctedCounts;
  }

  // extract effective downscalings and effective luminosities based on the trigger correlation of the CENT cluster
  auto foundCorrelationhist = fCorrelationHistograms.find("CENT");
  if(foundCorrelationhist != fCorrelationHistograms.end()) {
    auto corrhist = foundCorrelationhist->second;
    fEffectiveDownscaling["INT7"] = getEffectiveDownscaling(corrhist, "INT7", "EG1");
    fEffectiveDownscaling["EMC7"] = getEffectiveDownscaling(corrhist, "EMC7", "EG1");
    fEffectiveDownscaling["EG1"] = getEffectiveDownscaling(corrhist, "EG1", "EG1");
    fEffectiveDownscaling["EG2"] = getEffectiveDownscaling(corrhist, "EG2", "EG1");
    fEffectiveDownscaling["EJ1"] = getEffectiveDownscaling(corrhist, "EJ1", "EJ1");
    fEffectiveDownscaling["EJ2"] = getEffectiveDownscaling(corrhist, "EJ2", "EJ1");
    fEffectiveDownscaling["DMC7"] = getEffectiveDownscaling(corrhist, "DMC7", "DG1");
    fEffectiveDownscaling["DG1"] = getEffectiveDownscaling(corrhist, "DG1", "DG1");
    fEffectiveDownscaling["DG2"] = getEffectiveDownscaling(corrhist, "DG2", "DG1");
    fEffectiveDownscaling["DJ1"] = getEffectiveDownscaling(corrhist, "DJ1", "DJ1");
    fEffectiveDownscaling["DJ2"] = getEffectiveDownscaling(corrhist, "DJ2", "DJ2");

    TH1 *counterINT7 = fClusterCounters["INT7"];
    auto countsINT7 = getTriggerClusterCounts(counterINT7, "CENT"),
         refLuminosity = countsINT7 / fEffectiveDownscaling["INT7"];
    fEffectiveLuminosity["INT7"] = countsINT7 / xrefPB;
    for(auto &trgdownscale : fEffectiveDownscaling) {
      auto &triggerclass = trgdownscale.first;
      if(triggerclass == "INT7") continue;
      auto &downscaling = trgdownscale.second;
      double clustercorrection = 1.;
      if(triggerclass.find("MC7") != std::string::npos || triggerclass.find("G2") != std::string::npos || triggerclass.find("J2") != std::string::npos) {
        clustercorrection = corrCENTNOPMD;
      } else if(triggerclass.find("G1") != std::string::npos || triggerclass.find("J1") != std::string::npos) {
        clustercorrection = corrCombined;
      }
      fEffectiveLuminosity[triggerclass] = refLuminosity * downscaling * clustercorrection / xrefPB;
    }
  }
}

double AliEmcalTriggerLuminosity::getTriggerClusterCounts(TH1 * clustercounter, const std::string &clustername) const {
  int binID = clustercounter->GetXaxis()->FindBin(clustername.data());
  if(binID <= 0 || binID > clustercounter->GetXaxis()->GetNbins()) return 0;
  return clustercounter->GetBinContent(binID);
}

double AliEmcalTriggerLuminosity::getEffectiveDownscaling(TH2 *corrhist, const std::string &trigger, const std::string &reftrigger) const {
  int binID = corrhist->GetXaxis()->FindBin(reftrigger.data());
  std::unique_ptr<TH1> refhist(corrhist->ProjectionY("refhist", binID, binID));
  int bintrigger = refhist->GetXaxis()->FindBin(trigger.data()),
      binref = refhist->GetXaxis()->FindBin(reftrigger.data());
  double countstrigger = refhist->GetBinContent(bintrigger),
         countsref = refhist->GetBinContent(binref);
  return countstrigger / countsref;
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

double AliEmcalTriggerLuminosity::convertLuminosity(double luminosityPB, LuminosityUnit_t unit) const {
  return luminosityPB * getConversionToPB(unit);
}

double AliEmcalTriggerLuminosity::getConversionToPB(LuminosityUnit_t unit) const {
  double conversion = 1.;
  switch(unit){
    case LuminosityUnit_t::kPb: conversion = 1.; break;
    case LuminosityUnit_t::kNb: conversion = 1e3; break;
    case LuminosityUnit_t::kMub: conversion = 1e6; break; 
    case LuminosityUnit_t::kMb: conversion = 1e9; break;
    case LuminosityUnit_t::kB: conversion = 1e12; break;
  };
  return conversion;
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
