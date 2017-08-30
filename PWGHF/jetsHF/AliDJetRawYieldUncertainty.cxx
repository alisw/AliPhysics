/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <iostream>
#include <cstring>
#include <set>

#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TRandom2.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <Riostream.h>
#include <TLine.h>

class TNtuple;
#include "AliHFMultiTrials.h"

#include "AliDJetVReader.h"

#include "AliDJetRawYieldUncertainty.h"

/// \cond CLASSIMP
ClassImp(AliDJetRawYieldUncertainty);
/// \endcond

const Int_t AliDJetRawYieldUncertainty::fgkNSigmaVar;
const Int_t AliDJetRawYieldUncertainty::fgkNBkgVar;

/**
 * Default constructor.
 */
AliDJetRawYieldUncertainty::AliDJetRawYieldUncertainty():
  TObject(),
  fDJetReader(nullptr),
  fSaveInvMassFitCanvases(kFALSE),
  fDmesonSpecie(kUnknownMeson),
  fDmesonLabel(),
  fYieldApproach(kUnknownMethod),
  fMethodLabel(),
  fpTmin(0),
  fpTmax(0),
  fzmin(0),
  fzmax(0),
  fnDbins(0),
  fDbinpTedges(nullptr),
  fnJetPtbins(0),
  fJetPtBinEdges(nullptr),
  fnJetzbins(0),
  fJetzBinEdges(nullptr),
  fDEffValues(nullptr),
  fnSigmaSignReg(0),
  fnSigmaSideBandLeft1(0),
  fnSigmaSideBandLeft2(0),
  fnSigmaSideBandRight1(0),
  fnSigmaSideBandRight2(0),
  fSigmaToFixDPtBins(nullptr),
  fSigmaToFixJetPtBins(nullptr),
  fnRebinSteps(0),
  fRebinSteps(nullptr),
  fnMinMassSteps(0),
  fMinMassSteps(nullptr),
  fnMaxMassSteps(0),
  fMaxMassSteps(nullptr),
  fnSigmaBC(0),
  fSigmaBC(nullptr),
  fnMask(0),
  fMask(nullptr),
  fChi2Cut(0),
  fnMaxTrials(0),
  fAllowRepetitions(kFALSE),
  fFitRefl(kFALSE),
  fReflFilenameInput(),
  fSigMCFilenameInput(),
  fReflHistoName(),
  fSigMCHistoName(),
  fFixRiflOverS(0),
  fReflRangeL(0),
  fReflRangeR(0),
  fCoherentChoice(kFALSE),
  fUseBkgInBinEdges(kTRUE),
  fDebug(0),
  fMassPlot(nullptr),
  fMassVsJetPtPlot(nullptr),
  fMassVsJetzPlot(nullptr),
  fJetPtYieldCentral(nullptr),
  fJetPtYieldUnc(nullptr),
  fJetPtSpectrSBVars(nullptr),
  fJetPtSpectrSBDef(nullptr),
  fJetPtBinYieldDistribution(nullptr),
  fJetzYieldCentral(nullptr),
  fJetzYieldUnc(nullptr),
  fJetzSpectrSBVars(nullptr),
  fJetzSpectrSBDef(nullptr),
  fJetzBinYieldDistribution(nullptr),
  fSuccess(kFALSE),
  fCanvases()
{
  memset(fMeanSigmaVar, 0, sizeof(Bool_t)*fgkNSigmaVar);
  memset(fBkgVar, 0, sizeof(Bool_t)*fgkNBkgVar);
}

/**
 * Copy constructor.
 * @param[in] source Const reference to an object to copy from
 */
AliDJetRawYieldUncertainty::AliDJetRawYieldUncertainty(const AliDJetRawYieldUncertainty &source):
  TObject(),
  fDJetReader(nullptr),
  fSaveInvMassFitCanvases(source.fSaveInvMassFitCanvases),
  fDmesonSpecie(source.fDmesonSpecie),
  fDmesonLabel(source.fDmesonLabel),
  fYieldApproach(source.fYieldApproach),
  fMethodLabel(source.fMethodLabel),
  fpTmin(source.fpTmin),
  fpTmax(source.fpTmax),
  fzmin(source.fzmin),
  fzmax(source.fzmax),
  fnDbins(0),
  fDbinpTedges(nullptr),
  fnJetPtbins(0),
  fJetPtBinEdges(nullptr),
  fnJetzbins(0),
  fJetzBinEdges(nullptr),
  fDEffValues(nullptr),
  fnSigmaSignReg(source.fnSigmaSignReg),
  fnSigmaSideBandLeft1(source.fnSigmaSideBandLeft1),
  fnSigmaSideBandLeft2(source.fnSigmaSideBandLeft2),
  fnSigmaSideBandRight1(source.fnSigmaSideBandRight1),
  fnSigmaSideBandRight2(source.fnSigmaSideBandRight2),
  fSigmaToFixDPtBins(nullptr),
  fSigmaToFixJetPtBins(nullptr),
  fnRebinSteps(0),
  fRebinSteps(nullptr),
  fnMinMassSteps(0),
  fMinMassSteps(nullptr),
  fnMaxMassSteps(0),
  fMaxMassSteps(nullptr),
  fnSigmaBC(0),
  fSigmaBC(nullptr),
  fnMask(0),
  fMask(nullptr),
  fChi2Cut(source.fChi2Cut),
  fnMaxTrials(source.fnMaxTrials),
  fAllowRepetitions(source.fAllowRepetitions),
  fFitRefl(source.fFitRefl),
  fReflFilenameInput(source.fReflFilenameInput),
  fSigMCFilenameInput(source.fSigMCFilenameInput),
  fReflHistoName(source.fReflHistoName),
  fSigMCHistoName(source.fSigMCHistoName),
  fFixRiflOverS(source.fFixRiflOverS),
  fReflRangeL(source.fReflRangeL),
  fReflRangeR(source.fReflRangeR),
  fCoherentChoice(source.fCoherentChoice),
  fUseBkgInBinEdges(source.fUseBkgInBinEdges),
  fDebug(source.fDebug),
  fMassPlot(nullptr),
  fMassVsJetPtPlot(nullptr),
  fMassVsJetzPlot(nullptr),
  fJetPtYieldCentral(nullptr),
  fJetPtYieldUnc(nullptr),
  fJetPtSpectrSBVars(nullptr),
  fJetPtSpectrSBDef(nullptr),
  fJetPtBinYieldDistribution(nullptr),
  fJetzYieldCentral(nullptr),
  fJetzYieldUnc(nullptr),
  fJetzSpectrSBVars(nullptr),
  fJetzSpectrSBDef(nullptr),
  fJetzBinYieldDistribution(nullptr),
  fSuccess(kFALSE),
  fCanvases()
{
  if (source.fnDbins > 0) {
    fnDbins = source.fnDbins;
    fDbinpTedges = new Double_t[fnDbins+1];
    memcpy(fDbinpTedges, source.fDbinpTedges, sizeof(Double_t)*(fnDbins+1));
    fDEffValues = new Double_t[fnDbins+1];
    memcpy(fDEffValues, source.fDEffValues, sizeof(Double_t)*(fnDbins+1));
    fSigmaToFixDPtBins = new Double_t[fnDbins+1];
    memcpy(fSigmaToFixDPtBins, source.fSigmaToFixDPtBins, sizeof(Double_t)*(fnDbins+1));
  }
  if (source.fnJetPtbins > 0) {
    fnJetPtbins = source.fnJetPtbins;
    fJetPtBinEdges = new Double_t[fnJetPtbins+1];
    memcpy(fJetPtBinEdges, source.fJetPtBinEdges, sizeof(Double_t)*(fnJetPtbins+1));
    fSigmaToFixJetPtBins = new Double_t[fnJetPtbins+1];
    memcpy(fSigmaToFixJetPtBins, source.fSigmaToFixJetPtBins, sizeof(Double_t)*(fnJetPtbins+1));
  }
  if (source.fnRebinSteps > 0) {
    fnRebinSteps = source.fnRebinSteps;
    fRebinSteps = new Int_t[fnRebinSteps];
    memcpy(fRebinSteps, source.fRebinSteps, sizeof(Int_t)*fnRebinSteps);
  }
  if (source.fnMinMassSteps > 0) {
    fnMinMassSteps = source.fnMinMassSteps;
    fMinMassSteps = new Double_t[fnMinMassSteps];
    memcpy(fMinMassSteps, source.fMinMassSteps, sizeof(Double_t)*fnMinMassSteps);
  }
  if (source.fnMaxMassSteps > 0) {
    fnMaxMassSteps = source.fnMaxMassSteps;
    fMaxMassSteps = new Double_t[fnMaxMassSteps];
    memcpy(fMaxMassSteps, source.fMaxMassSteps, sizeof(Double_t)*fnMaxMassSteps);
  }
  if (source.fnSigmaBC > 0) {
    fnSigmaBC = source.fnSigmaBC;
    fSigmaBC = new Double_t[fnSigmaBC];
    memcpy(fSigmaBC, source.fSigmaBC, sizeof(Double_t)*fnSigmaBC);
  }
  if (source.fnMask > 0) {
    fnMask = source.fnMask;
    fMask = new Bool_t[fnMask];
    memcpy(fMask, source.fMask, sizeof(Bool_t)*fnMask);
  }

  memcpy(fMeanSigmaVar, source.fMeanSigmaVar, sizeof(Bool_t)*fgkNSigmaVar);
  memcpy(fBkgVar, source.fBkgVar, sizeof(Bool_t)*fgkNBkgVar);
}

/**
 * Destructor
 */
AliDJetRawYieldUncertainty::~AliDJetRawYieldUncertainty()
{
  ClearObjects();

  if (fMassPlot) delete fMassPlot;
  if (fMassVsJetPtPlot) delete fMassVsJetPtPlot;
  if (fMassVsJetzPlot) delete fMassVsJetzPlot;
  if (fJetPtYieldCentral) delete fJetPtYieldCentral;
  if (fJetPtYieldUnc) delete fJetPtYieldUnc;
  if (fJetPtSpectrSBVars) {
    for (int i = 0; i < fnMaxTrials; i++) delete fJetPtSpectrSBVars[i];
    delete[] fJetPtSpectrSBVars;
  }
  if (fJetPtSpectrSBDef) delete fJetPtSpectrSBDef;
  if (fJetPtBinYieldDistribution) {
    for (int i = 0; i < fnJetPtbins; i++) delete fJetPtBinYieldDistribution[i];
    delete[] fJetPtBinYieldDistribution;
  }
  if (fJetzYieldCentral) delete fJetzYieldCentral;
  if (fJetzYieldUnc) delete fJetzYieldUnc;
  if (fJetzSpectrSBVars) {
    for (int i = 0; i < fnMaxTrials; i++) delete fJetzSpectrSBVars[i];
    delete[] fJetzSpectrSBVars;
  }
  if (fJetzSpectrSBDef) delete fJetzSpectrSBDef;
  if (fJetzBinYieldDistribution) {
    for (int i = 0; i < fnJetzbins; i++) delete fJetzBinYieldDistribution[i];
    delete[] fJetzBinYieldDistribution;
  }
  for (auto c : fCanvases) delete c;
}

/**
 * Set the D meson specie (D0 or D*)
 * @param[in] k D meson specie
 * @return kTRUE on success, kFALSE otherwise
 */
Bool_t AliDJetRawYieldUncertainty::SetDmesonSpecie(EDMesonSpecies_t k)
{
  switch(k){
  case kD0toKpi:
    fDmesonLabel = "Dzero";
    break;
  case kDStarD0pi:
    fDmesonLabel = "Dstar";
    break;
  default:
    Printf("Error! D meson specie not correctly set!\n");
    return kFALSE;
  }

  fDmesonSpecie = k;
  return kTRUE;
}

/**
 * Set the yield extraction method (side-band or eff.scaled)
 * @param[in] meth Yield extraction method
 * @return kTRUE on success, kFALSE otherwise
 */
Bool_t AliDJetRawYieldUncertainty::SetYieldMethod(EYieldMethod_t meth)
{
  fYieldApproach = meth;
  switch (fYieldApproach){
  case kEffScale:
    fMethodLabel = "InvMassFit";
    break;
  case kSideband:
    fMethodLabel = "SideBand";
    break;
  default:
    Printf("Error! Method not correctly set!");
    return kFALSE;
  }

  return kTRUE;
}

/**
 * Extract the input mass plots.
 */
Bool_t AliDJetRawYieldUncertainty::ExtractInputMassPlot()
{
  TH1::AddDirectory(kFALSE);

  std::cout << "Configuration:\nD meson: " << fDmesonLabel << "\nMethod: " << fMethodLabel << std::endl;

  fDJetReader->SetPtBinEdgesForMassPlot(fpTmin, fpTmax);
  fDJetReader->SetDmesonPtBins(fnDbins, fDbinpTedges);
  fDJetReader->SetJetPtBins(fnJetPtbins, fJetPtBinEdges);
  fDJetReader->SetJetzBins(fnJetzbins, fJetzBinEdges);
  fDJetReader->SetDmesonEfficiency(fDEffValues);

  Bool_t success = kFALSE;
  switch (fYieldApproach) {
  case kEffScale:
    success = fDJetReader->ExtractInputMassPlotEffScale();
    break;
  case kSideband:
    success = fDJetReader->ExtractInputMassPlotSideband();
    break;
  default:
    Printf("Unknown method '%d'! Exiting...", fYieldApproach);
    return kFALSE;
  }

  if (success) {
    fMassPlot = fDJetReader->GetMassPlot();
    fMassVsJetPtPlot = fDJetReader->GetMassVsJetPtPlot();
    fMassVsJetzPlot = fDJetReader->GetMassVsJetzPlot();
    std::cout << "Extracted mass spectrum for fit variations" << std::endl;
    std::cout << "Mass spectrum entries: " << fMassPlot->GetEntries() << std::endl;
  }

  return success;
}

/**
 * Run the multi-trial analysis.
 */
AliHFMultiTrials* AliDJetRawYieldUncertainty::RunMultiTrial()
{
  TH1::AddDirectory(kFALSE);

  std::cout << "Running MultiTrial on pT bin" << fpTmin << " to " << fpTmax << std::endl;

  Double_t massD = 0.0;
  switch (fDmesonSpecie) {
  case kD0toKpi:
    massD = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    break;
  case kDStarD0pi:
    massD = TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass();
    break;
  default:
    break;
  }

  if(fDebug) std::cout << "D-meson mass: " << massD << std::endl;

  TString outfilnam = TString::Format("RawYieldVariations_%s_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fMethodLabel.Data(), fpTmin, fpTmax);

  Double_t sigmaToFix = 0;
  switch (fYieldApproach) {
  case kEffScale:
    for (Int_t i = 0; i < fnJetPtbins; i++){
      if (fpTmin == fJetPtBinEdges[i]) {
        sigmaToFix = fSigmaToFixJetPtBins[i];
        Printf("InvMassFit, Jet pt Bin %d (%.2f, %.2f), sigma = %.5f", i, fpTmin, fpTmax, sigmaToFix);
        break;
      }
    }
    break;
  case kSideband:
    for (Int_t i = 0; i < fnDbins; i++){
      if (fpTmin == fDbinpTedges[i]) {
        sigmaToFix = fSigmaToFixDPtBins[i];
        Printf("SB approach, D pt Bin %d (%.2f, %.2f), sigma = %.5f", i, fpTmin, fpTmax, sigmaToFix);
        break;
      }
    }
    break;
  default:
    break;
  }

  if (sigmaToFix == 0){
    Printf("**** ERROR **** sigmaToFix = 0!!");
    return nullptr;
  }
  AliHFMultiTrials* mt = new AliHFMultiTrials();
  mt->SetSuffixForHistoNames("");
  mt->SetMass(massD);
  mt->SetSigmaGaussMC(sigmaToFix);
  mt->SetUseFixSigFreeMean(fMeanSigmaVar[0]);
  mt->SetUseFixSigUpFreeMean(fMeanSigmaVar[1]);
  mt->SetUseFixSigDownFreeMean(fMeanSigmaVar[2]);
  mt->SetUseFreeS(fMeanSigmaVar[3]);
  mt->SetUseFixedMeanFreeS(fMeanSigmaVar[4]);
  mt->SetUseFixSigFixMean(fMeanSigmaVar[5]);
  mt->SetUseExpoBackground(fBkgVar[0]);
  mt->SetUseLinBackground(fBkgVar[1]);
  mt->SetUsePol2Background(fBkgVar[2]);
  mt->SetUsePol3Background(fBkgVar[3]);
  mt->SetUsePol4Background(fBkgVar[4]);
  mt->SetUsePol5Background(fBkgVar[5]);
  mt->SetUsePowerLawBackground(fBkgVar[6]);
  mt->SetUsePowerLawTimesExpoBackground(fBkgVar[7]);
  mt->ConfigureRebinSteps(fnRebinSteps,fRebinSteps);
  mt->ConfigureLowLimFitSteps(fnMinMassSteps,fMinMassSteps);
  mt->ConfigureUpLimFitSteps(fnMaxMassSteps,fMaxMassSteps);
  if (fSaveInvMassFitCanvases) {
    mt->AddInvMassFitSaveAsFormat("pdf");
    mt->AddInvMassFitSaveAsFormat("root");
  }
  Printf("sigmaBC=%d", fnSigmaBC);

  if (fnSigmaBC) mt->ConfigurenSigmaBinCSteps(fnSigmaBC, fSigmaBC);
  mt->SetSaveBkgValue(kTRUE, fnSigmaSignReg);
  mt->SetDrawIndividualFits(kFALSE);
  if (fDebug >= 2) mt->SetDrawIndividualFits(kTRUE);

  TString cname;
  cname = TString::Format("MassFit_%s_%s_%1.1fto%1.1f", fDmesonLabel.Data(), fMethodLabel.Data(), fpTmin, fpTmax);
  TCanvas* c0 = new TCanvas(cname,cname);
  fCanvases.push_back(c0);

  if (fFitRefl) { //reflection treatment
    std::cout << " Reflection template file: " << fReflFilenameInput.Data() << " - Reflection histo name:  " << fReflHistoName.Data() << std::endl;
    TFile fRefl(fReflFilenameInput.Data(), "read");
    if (fRefl.IsZombie()) {
      std::cout << " Reflection file not found! Exiting... " << std::endl;
      delete mt;
      return nullptr;
    }
    TH1 *hTemplRefl_temp = dynamic_cast<TH1*>(fRefl.Get(fReflHistoName.Data()));
    if (!hTemplRefl_temp) {
      std::cout << " Reflection histo not found! Exiting... " << std::endl;
      fRefl.Close();
      delete mt;
      return nullptr;
    }
    Printf("Copying histogram...");
    TH1F *hTemplRefl = new TH1F("temp", "temp", 1, 0, 1); // needed to convert other types of TH1 (TH1D in particular)
    hTemplRefl_temp->Copy(*hTemplRefl);
    Printf("Histogram '%s' copied successfully", hTemplRefl->GetName());
    fRefl.Close();
    Printf("File closed!");
    hTemplRefl_temp = nullptr;

    if (fFixRiflOverS < 0) {
      Double_t factor = -fFixRiflOverS;
      std::cout << " MC signal file: " << fSigMCFilenameInput.Data() << " - MC signal histo name:  " << fSigMCHistoName.Data() << std::endl;
      TFile fSigMC(fSigMCFilenameInput.Data(),"read"); //needed if refl/sigMC ratio is not fixed
      if (fSigMC.IsZombie()) {
        std::cout << " MC signal file not found! Exiting... " << std::endl;
        delete mt;
        return nullptr;
      }
      TH1 *hSignMC = dynamic_cast<TH1*>(fSigMC.Get(fSigMCHistoName.Data()));
      if (!hSignMC) {
        std::cout << " MC signal histo not found! Exiting... " << std::endl;
        fSigMC.Close();
        delete mt;
        return nullptr;
      }
      fFixRiflOverS = factor * hTemplRefl->Integral(hTemplRefl->FindBin(fReflRangeL + 0.00001), hTemplRefl->FindBin(fReflRangeR - 0.00001)) / hSignMC->Integral(hSignMC->FindBin(fReflRangeL + 0.00001), hSignMC->FindBin(fReflRangeR - 0.00001));
      Printf("Sign over Refl bin %1.1f-%1.1f = %f", fpTmin, fpTmax, fFixRiflOverS);
      fSigMC.Close();
    }

    mt->SetTemplateRefl(hTemplRefl);
    mt->SetFixRefoS(fFixRiflOverS);
  }

  Bool_t isOK = mt->DoMultiTrials(fMassPlot, c0);

  cname = TString::Format("AllTrials_%s_%s_%1.1fto%1.1f", fDmesonLabel.Data(), fMethodLabel.Data(), fpTmin, fpTmax);
  TCanvas* cOut = new TCanvas(cname, cname);
  fCanvases.push_back(cOut);
  std::cout << "Has MultiTrial suceeded? " << isOK << std::endl;

  if (isOK) {
    mt->DrawHistos(cOut);
    mt->SaveToRoot(outfilnam.Data(), "recreate");
  }
  else {
    delete mt;
    return nullptr;
  }

  fSuccess = CombineMultiTrialOutcomes();

  return mt;
}

/**
 * Combines outcomes of the multi-trial run.
 */
Bool_t AliDJetRawYieldUncertainty::CombineMultiTrialOutcomes()
{
  TH1::AddDirectory(kFALSE);

  TString infilnam = TString::Format("RawYieldVariations_%s_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fMethodLabel.Data(), fpTmin, fpTmax);
  TFile fil(infilnam.Data());

  TString confCase[fgkNSigmaVar] = { "", "", "", "", "", "" };
  TString bkgFunc[fgkNBkgVar] = { "", "", "", "", "", "", "", "" };

  Int_t nConfigCases = 0;
  for (int i = 0; i < fgkNSigmaVar; i++) {
    if (i == 0) confCase[nConfigCases] = "FixedS";
    else if (i == 1) confCase[nConfigCases] = "FixedSp20";
    else if (i == 2) confCase[nConfigCases] = "FixedSm20";
    else if (i == 3) confCase[nConfigCases] = "FreeS";
    else if (i == 4) confCase[nConfigCases] = "FixedMeanFreeS";
    else if (i == 5) confCase[nConfigCases] = "FixedMeanFixedS";
    if (fMeanSigmaVar[i]) nConfigCases++;
  }

  Int_t nBackFuncCases = 0;
  for(int i = 0; i < fgkNBkgVar; i++) {
    if (i == 0) bkgFunc[nBackFuncCases] = "Expo";
    if (i == 1) bkgFunc[nBackFuncCases] = "Lin";
    if (i == 2) bkgFunc[nBackFuncCases] = "Pol2";
    if (i == 3) bkgFunc[nBackFuncCases] = "Pol3";
    if (i == 4) bkgFunc[nBackFuncCases] = "Pol4";
    if (i == 5) bkgFunc[nBackFuncCases] = "Pol5";
    if (i == 6) bkgFunc[nBackFuncCases] = "PowLaw";
    if (i == 7) bkgFunc[nBackFuncCases] = "PowLawExpo";
    if (fBkgVar[i]) nBackFuncCases++;
  }

  Int_t totCases = nBackFuncCases * nConfigCases;
  if (totCases != fnMask) {
    std::cout << "Error in the configuration of the mask! Mismatch with the number of active sigma/mean and bkg types settings!" << std::endl;
    return kFALSE;
  }

  TH1F* histo[totCases];
  for (int i = 0; i < totCases; i++) histo[i] = nullptr;
  std::cout << " Total cases (sigma/mean * bkg configs): " << totCases << std::endl;

  Int_t jh = 0;
  for (Int_t iConf = 0; iConf < nConfigCases; iConf++){
    for( Int_t iType = 0; iType < nBackFuncCases; iType++){
      histo[jh++] = dynamic_cast<TH1F*>(fil.Get(Form("hRawYieldTrial%s%s", bkgFunc[iType].Data(), confCase[iConf].Data())));
      if (fDebug) std::cout << "Loading histo: " << Form("hRawYieldTrial%s%s",bkgFunc[iType].Data(),confCase[iConf].Data()) << std::endl;
    }
  }

  Int_t totTrials = 0;
  Int_t totHistos = 0;
  Int_t first[totCases];
  for (int i = 0; i < totCases; i++) first[i] = 0;
  for (Int_t j = 0; j < totCases; j++) {
    if (fMask[j]) {
      if (fDebug) std::cout << " case (with fMask active): " << j << std::endl;
      if (histo[j]) {
        first[j] = totTrials;
        totTrials += histo[j]->GetNbinsX();
        totHistos++;
      }
      else {
        fMask[j] = 0;
      }
    }
  }

  TLine* vlines[totCases];
  for (int i = 0; i < totCases; i++) vlines[i] = nullptr;
  if (fDebug) Printf("Histos merged = %d, totTrials = %d", totHistos, totTrials);
  for (Int_t j = 0; j < totCases; j++) {
    if (fMask[j]) {
      if (fDebug) Printf("  %d) %s  -- %d",j,histo[j]->GetName(), first[j]);
      vlines[j] = new TLine(first[j], 0., first[j], 50000.);
      vlines[j]->SetLineColor(kMagenta + 2);
      vlines[j]->SetLineStyle(2);
    }
  }

  TH1F* hRawYieldAll = new TH1F("hRawYieldAll", " ; Trial # ; Raw Yield", totTrials, -0.5, totTrials - 0.5);
  TH1F* hMeanAll = new TH1F("hMeanAll", " ; Trial # ; Gaussian mean", totTrials, -0.5, totTrials - 0.5);
  TH1F* hSigmaAll = new TH1F("hSigmaAll", " ; Trial # ; Gaussian #sigma", totTrials, -0.5, totTrials - 0.5);
  TH1F* hChi2All = new TH1F("hChi2All", " ; Trial # ; #chi^{2}", totTrials, -0.5, totTrials - 0.5);
  TH1F* hBkgAll = new TH1F("hBkgAll", " ; Trial # ; Background under peak", totTrials, -0.5, totTrials - 0.5);
  TH1F* hRawYieldDistAll = new TH1F("hRawYieldDistAll", "  ; Raw Yield", 25000, 0., 25000.);
  TH1F* hStatErrDistAll = new TH1F("hStatErrDistAll", "  ; Stat Unc on Yield", 1000, 0., 20000.);
  TH1F* hRelStatErrDistAll = new TH1F("hRelStatErrDistAll", "  ; Rel Stat Unc on Yield", 100, 0., 1.);

  Double_t minYield = 999999.;
  Double_t maxYield = 0.;
  Double_t sumy[4] = {0};
  Double_t sumwei[4] ={0};
  Double_t sumerr[4] = {0};
  Double_t counts = 0.;
  Double_t wei[4] = {0};

  if (fDebug) Printf("Building overall variation plot");
  for (Int_t j = 0; j < totCases; j++) {
    if (fMask[j]){
      if (fDebug) Printf("-> Case %d", j);
      TString hmeanname = histo[j]->GetName();
      hmeanname.ReplaceAll("RawYield", "Mean");
      TH1F* hmeant = dynamic_cast<TH1F*>(fil.Get(hmeanname.Data()));

      TString hsigmaname = histo[j]->GetName();
      hsigmaname.ReplaceAll("RawYield", "Sigma");
      TH1F* hsigmat = dynamic_cast<TH1F*>(fil.Get(hsigmaname.Data()));

      TString hchi2name = histo[j]->GetName();
      hchi2name.ReplaceAll("RawYield","Chi2");
      TH1F* hchi2t = dynamic_cast<TH1F*>(fil.Get(hchi2name.Data()));

      TString hbkgname = histo[j]->GetName();
      if (fUseBkgInBinEdges) {
        hbkgname.ReplaceAll("RawYield","BkgInBinEdges");
      }
      else {
        hbkgname.ReplaceAll("RawYield","Bkg");
      }
      TH1F* hbkg = dynamic_cast<TH1F*>(fil.Get(hbkgname.Data()));

      for (Int_t ib = 1; ib <= histo[j]->GetNbinsX(); ib++) {
        Double_t ry = histo[j]->GetBinContent(ib);//rawyield
        Double_t ery = histo[j]->GetBinError(ib);

        Double_t pos = hmeant->GetBinContent(ib);
        Double_t epos = hmeant->GetBinError(ib);

        Double_t sig = hsigmat->GetBinContent(ib);
        Double_t esig = hsigmat->GetBinError(ib);

        Double_t chi2 = hchi2t->GetBinContent(ib);
        Double_t bkg = hbkg->GetBinContent(ib);
        Double_t ebkg = hbkg->GetBinError(ib);

        if (fDebug) std::cout << " ry " << ry << " ery " << ery << " chi2 " << chi2 << std::endl;
        if(ry > 0.001 && ery > (0.01 * ry) && ery < (0.5 * ry) && chi2 < fChi2Cut) {
          hRawYieldDistAll->Fill(ry);
          hStatErrDistAll->Fill(ery);
          hRelStatErrDistAll->Fill(ery/ry);
          hRawYieldAll->SetBinContent(first[j]+ib,ry);
          hRawYieldAll->SetBinError(first[j]+ib,ery);
          if (ry < minYield) minYield = ry;
          if (ry > maxYield) maxYield = ry;
          wei[0] = 1.;
          wei[1] = 1. / (ery * ery);
          wei[2] = 1. / (ery * ery / (ry * ry));
          wei[3] = 1. / (ery * ery / ry);
          for (Int_t kw = 0; kw < 4; kw++){
            sumy[kw] += wei[kw] * ry;
            sumerr[kw] += wei[kw] * wei[kw] * ery * ery;
            sumwei[kw] += wei[kw];
          }
          counts += 1.;
          hSigmaAll->SetBinContent(first[j] + ib, sig);
          hSigmaAll->SetBinError(first[j] + ib, esig);
          hMeanAll->SetBinContent(first[j] + ib, pos);
          hMeanAll->SetBinError(first[j] + ib, epos);
          hChi2All->SetBinContent(first[j] + ib, chi2);
          hChi2All->SetBinError(first[j] + ib, 0.000001);
          hBkgAll->SetBinContent(first[j] + ib, bkg);
          hBkgAll->SetBinError(first[j] + ib, ebkg);
        }
      }
    }
  }

  Double_t weiav[4] = {0};
  Double_t eweiav[4] = {0};
  for (Int_t kw = 0; kw < 4; kw++){
    if(sumwei[kw] > 0.) {
      weiav[kw] = sumy[kw] / sumwei[kw];
      eweiav[kw] = TMath::Sqrt(sumerr[kw]) / sumwei[kw];
    }
  }

  hRawYieldAll->SetStats(0);
  hMeanAll->SetStats(0);
  hSigmaAll->SetStats(0);
  hChi2All->SetStats(0);
  hChi2All->SetMarkerStyle(7);
  hBkgAll->SetStats(0);
  hSigmaAll->SetMinimum(0.);

  TString cname = TString::Format("All_%s_%s_%1.1fto%1.1f", fDmesonLabel.Data(), fMethodLabel.Data(), fpTmin, fpTmax);
  TCanvas* call = new TCanvas(cname, cname, 1400 ,800);
  fCanvases.push_back(call);
  call->Divide(4,2);
  call->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hSigmaAll->GetYaxis()->SetTitleOffset(1.7);
  hSigmaAll->Draw();
  call->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hMeanAll->GetYaxis()->SetTitleOffset(1.7);
  hMeanAll->Draw();
  call->cd(3);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hChi2All->GetYaxis()->SetTitleOffset(1.7);
  hChi2All->Draw();
  call->cd(4);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  hBkgAll->GetYaxis()->SetTitleOffset(1.7);
  hBkgAll->Draw();
  call->cd(5);
  hRawYieldAll->SetTitle(Form("%s",infilnam.Data()));
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.06);
  Double_t newmax = 1.25 * (hRawYieldAll->GetMaximum() + hRawYieldAll->GetBinError(1));
  hRawYieldAll->GetYaxis()->SetTitleOffset(1.7);
  hRawYieldAll->Draw();
  TLatex* tweimean[4] = {0};
  for (Int_t kw = 0; kw < 4; kw++){
    tweimean[kw] = new TLatex(0.16, 0.84 - 0.06 * kw, Form("<Yield>_{wei%d} = %.1f #pm %.1f", kw, weiav[kw], eweiav[kw] * TMath::Sqrt(counts)));
    tweimean[kw]->SetNDC();
    tweimean[kw]->SetTextColor(4);
  }

  for(Int_t j = 1; j < totCases; j++){
    if (fMask[j]){
      vlines[j]->SetY2(newmax);
      vlines[j]->Draw("same");
    }
  }

  call->cd(6);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  hRawYieldDistAll->SetTitle(Form("%s", infilnam.Data()));
  hRawYieldDistAll->Draw();
  hRawYieldDistAll->GetXaxis()->SetRangeUser(minYield * 0.8, maxYield * 1.2);
  Double_t perc[3] = {0.15, 0.5, 0.85}; // quantiles for +-1 sigma
  Double_t lim70[3] = {0};
  hRawYieldDistAll->GetQuantiles(3, lim70, perc);
  call->cd(7);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  Double_t aver = hRawYieldDistAll->GetMean();
  TLatex* tmean = new TLatex(0.15, 0.94, Form("mean=%.1f", aver));
  tmean->SetNDC();
  tmean->Draw();
  TLatex* tmedian = new TLatex(0.15, 0.86, Form("median=%.1f", lim70[1]));
  tmedian->SetNDC();
  tmedian->Draw();
  Double_t val = hRawYieldDistAll->GetRMS();
  TLatex* thrms = new TLatex(0.15, 0.78, Form("rms=%.1f  (%.1f%%)", val, val / aver * 100.));
  Double_t rms=val/aver*100.;
  thrms->SetNDC();
  thrms->Draw();
  TLatex* tmin = new TLatex(0.15, 0.64, Form("min=%.1f      max=%.1f", minYield, maxYield));
  tmin->SetNDC();
  tmin->Draw();
  val = (maxYield - minYield) / TMath::Sqrt(12);
  TLatex* trms = new TLatex(0.15, 0.56, Form("(max-min) / sqrt(12)=%.1f  (%.1f%%)", val, val / aver * 100.));
  trms->SetNDC();
  trms->Draw();
  val = (maxYield - aver) / TMath::Sqrt(3);
  TLatex* tup = new TLatex(0.15, 0.48, Form("(max-mean)/sqrt(3)=%.1f  (%.1f%%)", val, val / aver * 100.));
  tup->SetNDC();
  tup->Draw();
  val=(aver-minYield)/sqrt(3);
  TLatex* tdw = new TLatex(0.15, 0.40, Form("(mean-min)/sqrt(3)=%.1f  (%.1f%%)", val, val / aver * 100.));
  tdw->SetNDC();
  tdw->Draw();
  TLatex* tl15 = new TLatex(0.15, 0.26, Form("15 percentile=%.1f", lim70[0]));
  tl15->SetNDC();
  tl15->Draw();
  TLatex* tl85 = new TLatex(0.15, 0.18, Form("85 percentile=%.1f", lim70[2]));
  tl85->SetNDC();
  tl85->Draw();
  val = (lim70[2] - lim70[0]) / 2.;
  TLatex* t1s = new TLatex(0.15, 0.1, Form("70%% range =%.1f  (%.1f%%)", val, val / aver * 100.));
  t1s->SetNDC();
  t1s->Draw();

  for (Int_t kw = 0; kw < 4; kw++) {
    if (fDebug) Printf("Weight %d: %.1f +- %.1f(stat) +- %.1f (syst)\n", kw, weiav[kw], eweiav[kw] * TMath::Sqrt(counts), (maxYield - minYield) / TMath::Sqrt(12));
  }

  TString outfilnam;
  TFile* outfile = nullptr;

  outfilnam = TString::Format("RawYieldSyst_%s_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fMethodLabel.Data(), fpTmin, fpTmax);
  outfile = TFile::Open(outfilnam, "recreate");
  if (!outfile || outfile->IsZombie()) {
    Printf("Could not create file '%s'!", outfilnam.Data());
  }
  else {
    Printf("Saving results in '%s'", outfilnam.Data());
    outfile->cd();
    call->Write();
    if (fMassVsJetPtPlot) fMassVsJetPtPlot->Write();
    if (fMassVsJetzPlot) fMassVsJetzPlot->Write();
    outfile->Close();
    delete outfile;
    outfile = nullptr;
  }

  TH1D *hUnc = new TH1D(Form("hUnc_%1.1fto%1.1f", fpTmin, fpTmax), "hUnc", 1, 0., 1.);
  hUnc->SetBinContent(1, aver);
  hUnc->SetBinError(1, rms * aver / 100.);

  outfilnam = TString::Format("Hist_%s", outfilnam.Data());
  outfile = TFile::Open(outfilnam, "recreate");
  if (!outfile || outfile->IsZombie()) {
    Printf("Could not create file '%s'!", outfilnam.Data());
  }
  else {
    Printf("Saving results in '%s'", outfilnam.Data());
    outfile->cd();
    hUnc->Write();
    outfile->Close();
    delete outfile;
    outfile = nullptr;
  }

  fil.Close();

  Printf("Combine multi-trials done!");
  return kTRUE;
}

/**
 * Evaluate final uncertainty.
 */
Bool_t AliDJetRawYieldUncertainty::EvaluateUncertainty()
{
  TH1::AddDirectory(kFALSE);

  Bool_t success = kTRUE;
  switch (fYieldApproach) {
  case kEffScale:
    success = EvaluateUncertaintyEffScale();
    break;
  case kSideband:
    success = EvaluateUncertaintySideband();
    break;
  default:
    break;
  }

  if (success) std::cout << "Evaluated raw yield uncertainty" << std::endl;

  return success;
}

/**
 * Evaluate final uncertainty for the efficiency scaled method.
 */
Bool_t AliDJetRawYieldUncertainty::EvaluateUncertaintyEffScale()
{
  std::cout << "Jet spectrum pT bin edges: ";
  for (int i = 0; i < fnJetPtbins; i++) std::cout << fJetPtBinEdges[i] << " - ";
  std::cout << fJetPtBinEdges[fnJetPtbins] << std::endl;

  fJetPtYieldUnc = new TH1D("JetPtRawYieldUncert", "Raw yield uncertainty on jet pt spectrum - Dzero - Eff.scaling", fnJetPtbins, fJetPtBinEdges);
  fJetPtYieldCentral = new TH1D("JetPtRawYieldCentral", "Jet pt spectrum central values + syst yield uncertainty - Dzero - Eff.scaling", fnJetPtbins, fJetPtBinEdges);

  TH1D *hUnc = 0;

  // loop over the jet pT bins already extracted
  for (int ibin = 0; ibin < fnJetPtbins; ibin++) {
    TString fname = TString::Format("Hist_RawYieldSyst_%s_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fMethodLabel.Data(), fJetPtBinEdges[ibin], fJetPtBinEdges[ibin+1]);
    if (fDebug) std::cout << fname.Data() << std::endl;
    TFile f(fname, "read");
    if (f.IsZombie()) {
      std::cout << "Uncertainty file for bin " << fJetPtBinEdges[ibin] << " - " << fJetPtBinEdges[ibin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return kFALSE;
    }
    else {
      if (fDebug) std::cout << "Building uncertainty for bin " << fJetPtBinEdges[ibin] << " - " << fJetPtBinEdges[ibin+1] << std::endl;
    }

    hUnc = dynamic_cast<TH1D*>(f.Get(Form("hUnc_%1.1fto%1.1f", fJetPtBinEdges[ibin], fJetPtBinEdges[ibin+1])));
    if (!hUnc) {
      std::cout << "Histogram with uncertainty from mass plot not found! Returning..." << std::endl;
      f.Close();
      return kFALSE;
    }

    Double_t centrYield = hUnc->GetBinContent(1);
    Double_t rmsPct = hUnc->GetBinError(1);

    fJetPtYieldUnc->SetBinContent(ibin+1, rmsPct);
    fJetPtYieldCentral->SetBinContent(ibin+1, centrYield);
    fJetPtYieldCentral->SetBinError(ibin+1, rmsPct);
    f.Close();
  }

  fJetPtYieldUnc->SetStats(kFALSE);
  fJetPtYieldUnc->Draw();
  fJetPtYieldUnc->SaveAs(Form("FinalRawYieldUncertaintyJetPt_%s_%s.root", fDmesonLabel.Data(), fMethodLabel.Data()));
  fJetPtYieldCentral->SetStats(kFALSE);
  fJetPtYieldCentral->Draw();
  fJetPtYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertaintyJetPt_%s_%s.root", fDmesonLabel.Data(), fMethodLabel.Data()));

  // print distribution of yields for each variation
  for (int ibin = 0; ibin < fnJetPtbins; ibin++) {
    TFile f2(Form("RawYieldSyst_%s_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fMethodLabel.Data(), fJetPtBinEdges[ibin], fJetPtBinEdges[ibin+1]), "read");
    TString cname = TString::Format("All_%s_%s_%1.1fto%1.1f", fDmesonLabel.Data(), fMethodLabel.Data(), fJetPtBinEdges[ibin], fJetPtBinEdges[ibin+1]);
    TCanvas *c = dynamic_cast<TCanvas*>(f2.Get(cname));
    TH1F *hDist = dynamic_cast<TH1F*>(c->FindObject("hRawYieldDistAll"));
    hDist->SetStats(kTRUE);
    hDist->SaveAs(Form("YieldDistribution_%s_%s_%1.1fto%1.1f.root", fDmesonLabel.Data(), fMethodLabel.Data(), fJetPtBinEdges[ibin], fJetPtBinEdges[ibin+1]));
    f2.Close();
  }

  return kTRUE;
}

/**
 * Generate the jet pt spectrum for a given D-meson pt bin (side-band method).
 * @param[in] hInvMassJetObs Valid pointer to a 2-dimensional histogram with x = inv.mass and y = jet observable (pt, z, etc.)
 * @param[in] mean Mean of the invariant mass fit in the D meson pt bin
 * @param[in] sigma Sigma of the invariant mass fit in the D meson pt bin
 * @param[in] bkg Background of the invariant mass fit in the D meson pt bin
 * @param[in] iDbin D meson pt bin
 * @param[out] hjetobs Valid pointer to a histogram which will be filled with the subtracted jet spectrum
 * @param[out] hjetobs_s Valid pointer to a histogram which will be filled with the side-band jet spectrum
 * @param[out] hjetobs_s1 Valid pointer to a histogram which will be filled with the side-band jet spectrum (left)
 * @param[out] hjetobs_s2 Valid pointer to a histogram which will be filled with the side-band jet spectrum (right)
 * @return kTRUE if successful
 */
Bool_t AliDJetRawYieldUncertainty::GenerateJetSpectrum(TH2* hInvMassJetObs, Double_t mean, Double_t sigma, Double_t bkg, Int_t iDbin, TH1* hjetobs, TH1* hjetobs_s, TH1* hjetobs_s1, TH1* hjetobs_s2)
{
  Double_t jetmin = hjetobs->GetXaxis()->GetXmin();
  Double_t jetmax = hjetobs->GetXaxis()->GetXmax();

  Float_t signal_l_min = mean - fnSigmaSideBandLeft1 * sigma;
  Float_t signal_l_max = mean - fnSigmaSideBandLeft2 * sigma;
  Float_t signal_u_min = mean + fnSigmaSideBandRight1 * sigma;
  Float_t signal_u_max = mean + fnSigmaSideBandRight2 * sigma;
  Float_t signal_c_min = mean - fnSigmaSignReg * sigma;
  Float_t signal_c_max = mean + fnSigmaSignReg * sigma;

  //extract signal and sideband region spectra
  TH1* tmphjet = hInvMassJetObs->ProjectionY(Form("tmphjetobs%d",iDbin), hInvMassJetObs->GetXaxis()->FindBin(signal_c_min), hInvMassJetObs->GetXaxis()->FindBin(signal_c_max));
  TH1* tmphjet_s1 = hInvMassJetObs->ProjectionY(Form("tmphjetobs_s1%d",iDbin), hInvMassJetObs->GetXaxis()->FindBin(signal_l_min), hInvMassJetObs->GetXaxis()->FindBin(signal_l_max));
  TH1* tmphjet_s2 = hInvMassJetObs->ProjectionY(Form("tmphjetobs_s2%d",iDbin), hInvMassJetObs->GetXaxis()->FindBin(signal_u_min), hInvMassJetObs->GetXaxis()->FindBin(signal_u_max));
  TH1* tmphjet_s = static_cast<TH1*>(tmphjet_s1->Clone(Form("tmphjetobs_s%d",iDbin)));
  tmphjet_s->Add(tmphjet_s2);

  // scale background from side bands to the background under the peak
  if (tmphjet_s->Integral() == 0) {
    std::cout << "Error! At least one variation with no entries! Exiting..." << std::endl;
    return kFALSE;
  }
  Double_t scaling = bkg / tmphjet_s->Integral(tmphjet_s->FindBin(jetmin+0.0001), tmphjet_s->FindBin(jetmax-0.0001)); //integral btw jetmin and jetmax (where you get the bkg from the mass plot)
  Printf("Background scaling factor = %.6f", scaling);
  for (int j = 1; j <= tmphjet->GetNbinsX(); j++) {
    Double_t centerbin = tmphjet->GetBinCenter(j);
    hjetobs->Fill(centerbin, tmphjet->GetBinContent(j));
    hjetobs_s1->Fill(centerbin, tmphjet_s1->GetBinContent(j));
    hjetobs_s2->Fill(centerbin, tmphjet_s2->GetBinContent(j));
    hjetobs_s->Fill(centerbin, tmphjet_s->GetBinContent(j));
  }
  for (int j = 1; j <= hjetobs->GetNbinsX(); j++) {
    hjetobs->SetBinError(j, TMath::Sqrt(hjetobs->GetBinContent(j)));
    hjetobs_s1->SetBinError(j, TMath::Sqrt(hjetobs_s1->GetBinContent(j)));
    hjetobs_s2->SetBinError(j, TMath::Sqrt(hjetobs_s2->GetBinContent(j)));
    hjetobs_s->SetBinError(j, TMath::Sqrt(hjetobs_s->GetBinContent(j)));
  }

  hjetobs_s->Scale(scaling);

  // subtract background from signal jet
  hjetobs->Add(hjetobs_s, -1);

  // correct for D* efficiency
  hjetobs->Scale(1. / fDEffValues[iDbin]); // D efficiency
  hjetobs->SetMarkerColor(kBlue + 3);
  hjetobs->SetLineColor(kBlue + 3);

  //Normalize to full range the signal range (from fnSigmaSignReg range)
  Double_t normNsigma = 1.0;
  if (!fUseBkgInBinEdges) {
    normNsigma -= TMath::Erfc(fnSigmaSignReg / TMath::Sqrt2());
  }
  else {
    Double_t effSigma1 = (mean - hInvMassJetObs->GetXaxis()->GetBinLowEdge(hInvMassJetObs->GetXaxis()->FindBin(mean - fnSigmaSignReg * sigma))) / sigma;
    Double_t effSigma2 = (hInvMassJetObs->GetXaxis()->GetBinUpEdge(hInvMassJetObs->GetXaxis()->FindBin(mean + fnSigmaSignReg * sigma)) - mean) / sigma;
    normNsigma -= TMath::Erfc(effSigma1 / TMath::Sqrt2()) / 2 + TMath::Erfc(effSigma2 / TMath::Sqrt2()) / 2;
    Printf("The left effective sigma is %.3f. The right effective sigma is %.3f.", effSigma1, effSigma2);
  }
  hjetobs->Scale(1.0 / normNsigma);

  return kTRUE;
}

/**
 * Evaluate final uncertainty for the side-band method.
 * @return kTRUE if successful
 */
Bool_t AliDJetRawYieldUncertainty::EvaluateUncertaintySideband()
{
  Bool_t s = kTRUE;

  SBResults jetPtRes = EvaluateUncertaintySideband("Pt", fnJetPtbins, fJetPtBinEdges);
  fJetPtYieldCentral = jetPtRes.fJetYieldCentral;
  fJetPtYieldUnc = jetPtRes.fJetYieldUnc;
  fJetPtSpectrSBVars = jetPtRes.fJetSpectrSBVars;
  fJetPtSpectrSBDef = jetPtRes.fJetSpectrSBDef;
  fJetPtBinYieldDistribution = jetPtRes.fJetBinYieldDistribution;
  s = s && jetPtRes.fSuccess;

  SBResults jetzRes = EvaluateUncertaintySideband("z", fnJetzbins, fJetzBinEdges);
  fJetzYieldCentral = jetzRes.fJetYieldCentral;
  fJetzYieldUnc = jetzRes.fJetYieldUnc;
  fJetzSpectrSBVars = jetzRes.fJetSpectrSBVars;
  fJetzSpectrSBDef = jetzRes.fJetSpectrSBDef;
  fJetzBinYieldDistribution = jetzRes.fJetBinYieldDistribution;
  s = s && jetzRes.fSuccess;

  return s;
}

/**
 * Evaluate final uncertainty for the side-band method.
 * @return kTRUE if successful
 */
AliDJetRawYieldUncertainty::SBResults AliDJetRawYieldUncertainty::EvaluateUncertaintySideband(TString obs, Int_t nJetBins, Double_t* jetBinEdges)
{
  TRandom2 gen;
  gen.SetSeed(0);

  TH1F* jetSpectrSBDef = 0;
  TH1D* jetYieldUnc = 0;
  TH1D* jetYieldCentral = 0;
  TH1F** jetBinYieldDistribution = 0;

  //define list of histograms, one per variation
  TH1F** jetSpectrSBVars = new TH1F*[fnMaxTrials];
  for(Int_t k=0; k<fnMaxTrials; k++) jetSpectrSBVars[k] = nullptr;

  //for debug and thorough studies
  Double_t arrYldBinPerBin[fnDbins][nJetBins][fnMaxTrials];
  for (Int_t i = 0; i < fnDbins; i++) for(Int_t j=0; j<nJetBins; j++) for(Int_t k=0; k<fnMaxTrials; k++) arrYldBinPerBin[i][j][k] = 0;

  TString fname;
  TString cname;
  for (int iDbin = 0; iDbin < fnDbins; iDbin++) {
    if (fDebug) std::cout << "Running bin pT(D) " << iDbin << std::endl;

    TH1* hjet = new TH1F(Form("hjet%d",iDbin), "hjet_signReg_Rebinned", nJetBins, jetBinEdges);
    TH1* hjet_s1 = new TH1F(Form("hjet_s1%d",iDbin), "hjet_sb1_Rebinned", nJetBins, jetBinEdges);
    TH1* hjet_s2 = new TH1F(Form("hjet_s2%d",iDbin), "hjet_sb2_Rebinned", nJetBins, jetBinEdges);
    TH1* hjet_s = new TH1F(Form("hjet_s%d",iDbin), "hjet_sb2_Rebinned", nJetBins, jetBinEdges);
    hjet->Sumw2();
    hjet_s1->Sumw2();
    hjet_s2->Sumw2();
    hjet_s->Sumw2();

    //open file with summary of variations from MultiTrial - get histos of variations
    fname = TString::Format("RawYieldSyst_%s_%s_%1.1fto%1.1f.root",fDmesonLabel.Data(),fMethodLabel.Data(),fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]);
    TFile fileMult(fname, "read");
    if (fileMult.IsZombie()) {
      std::cout << "Uncertainty file for bin " << fDbinpTedges[iDbin] << " - " << fDbinpTedges[iDbin+1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
    }
    if (fDebug) Printf("File '%s' open successfully.", fname.Data());
    TString hname = TString::Format("hInvMassJet%s", obs.Data());
    TH2* hInvMassJet = dynamic_cast<TH2*>(fileMult.Get(hname));
    if (!hInvMassJet) {
      Printf("Could not find histogram %s!", hname.Data());
      fileMult.Close();
      return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
    }
    cname = Form("All_%s_%s_%1.1fto%1.1f", fDmesonLabel.Data(), fMethodLabel.Data(),fDbinpTedges[iDbin],fDbinpTedges[iDbin+1]);
    TCanvas *c = dynamic_cast<TCanvas*>(fileMult.Get(cname));
    if (!c) {
      Printf("Could not find canvas '%s'!", cname.Data());
      fileMult.Close();
      return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
    }
    TH1F *hMean = dynamic_cast<TH1F*>(c->FindObject("hMeanAll"));
    TH1F *hSigma = dynamic_cast<TH1F*>(c->FindObject("hSigmaAll"));
    TH1F *hBkg = dynamic_cast<TH1F*>(c->FindObject("hBkgAll"));
    fileMult.Close();

    if (fDebug) Printf("File '%s' closed successfully.", fname.Data());

    fname = TString::Format("RawYieldVariations_Dzero_SideBand_%1.1fto%1.1f.root", fDbinpTedges[iDbin], fDbinpTedges[iDbin + 1]);
    TFile fileMultVar(fname, "read");
    if (fileMultVar.IsZombie()) {
      std::cout << "Uncertainty file for bin " << fDbinpTedges[iDbin] << " - " << fDbinpTedges[iDbin + 1] << " cannot be opened! Did you already evaluate it?" << std::endl;
      return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
    }
    if (fDebug) Printf("File '%s' open successfully.", fname.Data());
    if (fDebug) std::cout << "Default trial" << " TrialExpoFreeS" << std::endl;
    TH1F *hMeanDef = static_cast<TH1F*>(fileMultVar.Get("hMeanTrialExpoFreeS"));
    TH1F *hSigmaDef = static_cast<TH1F*>(fileMultVar.Get("hSigmaTrialExpoFreeS"));
    TH1F *hBkgDef = 0;
    if (fUseBkgInBinEdges) {
      hBkgDef = static_cast<TH1F*>(fileMultVar.Get("hBkgInBinEdgesTrialExpoFreeS"));
    }
    else {
      hBkgDef = static_cast<TH1F*>(fileMultVar.Get("hBkgTrialExpoFreeS"));
    }

    TH1F *hRawYieldDef = static_cast<TH1F*>(fileMultVar.Get("hRawYieldTrialExpoFreeS"));
    fileMultVar.Close();
    if (fDebug) Printf("File '%s' closed successfully.", fname.Data());

    Double_t meanDef = hMeanDef->GetBinContent(1);
    Double_t sigmaDef = hSigmaDef->GetBinContent(1);
    Double_t bkgDef = hBkgDef->GetBinContent(1);
    Double_t rawYieldDef = hRawYieldDef->GetBinContent(1);

    std::cout << "Mean " << meanDef << ", sigma " << sigmaDef << ", bkg " << bkgDef << ", raw yield " << rawYieldDef << std::endl;
    if (sigmaDef == 0) {
      std::cout << "Error while generating spectrum for TrialExpoFreeS variation: trial did not converge!" << std::endl;
    }

    Bool_t resDefSpectrum = GenerateJetSpectrum(hInvMassJet, meanDef, sigmaDef, bkgDef, iDbin, hjet, hjet_s, hjet_s1, hjet_s2);
    if (!resDefSpectrum) {
      std::cout << "Error while generating spectrum for TrialExpoFreeS variation" << std::endl;
      return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
    }

    // add 'iDbin' pT(D) bin to total spectrum for variation 'iTrial'
    if (!iDbin) jetSpectrSBDef = static_cast<TH1F*>(hjet->Clone(TString::Format("fJet%sSpectrSBDef", obs.Data())));
    else jetSpectrSBDef->Add(hjet);

    hjet->SaveAs(Form("TrialExpoFreeS_Jet%s_%s_%s_%d.root", obs.Data(), fDmesonLabel.Data(),fMethodLabel.Data(), iDbin));

    hjet->Reset();
    hjet_s1->Reset();
    hjet_s2->Reset();
    hjet_s->Reset();

    if (fnMaxTrials > 0) {
      if (!fAllowRepetitions && fnMaxTrials > hMean->GetNbinsX()) {
        std::cout << "Error! you set more set spectrum total variations than those done for pT(D) bin" << fDbinpTedges[iDbin] << " - " << fDbinpTedges[iDbin+1] << "! ";
        std::cout << "Impossible to do without allowing repetitions! Exiting..." << std::endl;
        return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
      }

      if (fDebug) std::cout << "Running bin pT(D) " << iDbin << std::endl;

      std::set<int> extracted;
      for (int iTrial = 0; iTrial < fnMaxTrials; iTrial++) {
        Bool_t extractOk = kFALSE;
        Int_t jTrial = -1;

        if (fCoherentChoice) {
          jTrial = iTrial;
        }
        else {
          do {  //just one time if fAllowRepetitions==kTRUE, repeat extraction till new number is obtained if fAllowRepetitions==kFALSE
            jTrial = gen.Integer(hMean->GetNbinsX()) + 1;

            //avoid 'empty' cases
            if (hSigma->GetBinContent(jTrial) > 0) extractOk = kTRUE;

            //check if already extracted for this pT(D) bin
            if (!fAllowRepetitions && extracted.find(jTrial) != extracted.end()) extractOk = kFALSE;
          } while (extractOk == kFALSE);
        }
        extracted.insert(jTrial);

        Double_t mean = hMean->GetBinContent(jTrial);
        Double_t sigma = hSigma->GetBinContent(jTrial);
        Double_t bkg = hBkg->GetBinContent(jTrial);

        std::cout << "Mean " << mean << ", sigma " << sigma << ", bkg " << bkg << std::endl;

        Bool_t resSpectrum = GenerateJetSpectrum(hInvMassJet, mean, sigma, bkg, iDbin, hjet, hjet_s, hjet_s1, hjet_s2);
        if (!resSpectrum) {
          std::cout << "Error while generating spectrum for one of the variations" << std::endl;
          return {kFALSE, nullptr, nullptr, nullptr, nullptr, nullptr};
        }

        // for every trial of every pT(D) bin, save the value of the yield, after eff scaling, in each pT(jet) bin (to study pT(D)->pT(jet) yield correlations)
        for (int l = 0; l < hjet->GetNbinsX(); l++) arrYldBinPerBin[iDbin][l][iTrial] = hjet->GetBinContent(l + 1);

        // add 'iDbin' pT(D) bin to total spectrum for variation 'iTrial'
        if( !iDbin) jetSpectrSBVars[iTrial] = static_cast<TH1F*>(hjet->Clone(Form("Jet%sRawYieldUncert_%d", obs.Data(), iTrial)));
        else jetSpectrSBVars[iTrial]->Add(hjet);

        hjet->Reset();
        hjet_s1->Reset();
        hjet_s2->Reset();
        hjet_s->Reset();
      } //end loop on trials for sideband approach
    }
    delete hjet;
    delete hjet_s1;
    delete hjet_s2;
    delete hjet_s;
  } //end loop on pT(D) bins

  jetSpectrSBDef->SetStats(kFALSE);
  jetSpectrSBDef->Draw();
  jetSpectrSBDef->SaveAs(Form("TrialExpoFreeS_Jet%s_%s_%s.root", obs.Data(), fDmesonLabel.Data(),fMethodLabel.Data()));
  if (fnMaxTrials > 0) {
    //Now evaluate central value + rms in each pT(jet) bin to build the uncertainty
    Double_t arrYld[nJetBins][fnMaxTrials];
    for (Int_t i = 0; i < nJetBins; i++) for (Int_t j = 0; j < fnMaxTrials; j++) arrYld[i][j] = 0;

    jetBinYieldDistribution = new TH1F*[nJetBins];
    for (int i = 0; i < nJetBins; i++) jetBinYieldDistribution[i] = nullptr;

    jetYieldUnc = static_cast<TH1D*>(jetSpectrSBVars[0]->Clone(TString::Format("Jet%sRawYieldUncert", obs.Data())));
    jetYieldUnc->Reset();
    jetYieldUnc->SetTitle("Raw yield uncertainty on jet spectrum - Dstar - Sideband subtraction");
    jetYieldCentral =  static_cast<TH1D*>(jetSpectrSBVars[0]->Clone(TString::Format("Jet%sRawYieldCentral", obs.Data())));
    jetYieldCentral->Reset();
    jetYieldCentral->SetTitle("Jet spectrum central values + syst yield uncertainty - Dstar - Sideband subtraction");

    for (Int_t iJetbin = 0; iJetbin < nJetBins; iJetbin++) { //loop on jet spectrum pT bins

      jetBinYieldDistribution[iJetbin] = new TH1F(Form("fJet%sBinYieldDistribution_Bin%d", obs.Data(), iJetbin), "  ; Yield distribution", 50000, 0., 50000.);

      for (Int_t iTrial = 0; iTrial < fnMaxTrials; iTrial++) { //loop on trials and build array of variations for a given pT(jet) bin
        arrYld[iJetbin][iTrial] = jetSpectrSBVars[iTrial]->GetBinContent(iJetbin + 1);
        jetBinYieldDistribution[iJetbin]->Fill(arrYld[iJetbin][iTrial]);
      }

      Double_t mean = TMath::Mean(fnMaxTrials, arrYld[iJetbin]);
      Double_t rms = TMath::RMS(fnMaxTrials, arrYld[iJetbin]);
      if (fDebug) {
        std::cout << "Jet bin " << iJetbin << " (" << jetSpectrSBVars[0]->GetXaxis()->GetBinLowEdge(iJetbin+1) << "-" << jetSpectrSBVars[0]->GetXaxis()->GetBinUpEdge(iJetbin+1) << ")";
        std::cout << ": Mean = " << mean << ", RMS = " << rms << std::endl;
      }

      jetYieldUnc->SetBinContent(iJetbin + 1, rms);
      jetYieldCentral->SetBinContent(iJetbin + 1, mean);
      jetYieldCentral->SetBinError(iJetbin + 1, rms);

      jetBinYieldDistribution[iJetbin]->SaveAs(Form("YieldDistributionJet%s_%s_%s_%1.1fto%1.1f.root", obs.Data(), fDmesonLabel.Data(), fMethodLabel.Data(), jetBinEdges[iJetbin], jetBinEdges[iJetbin + 1]));
    }

    jetYieldUnc->SetStats(kFALSE);
    jetYieldUnc->Draw();
    jetYieldUnc->SaveAs(Form("FinalRawYieldUncertainty_Jet%s_%s_%s.root", obs.Data(), fDmesonLabel.Data(), fMethodLabel.Data()));
    jetYieldCentral->SetStats(kFALSE);
    jetYieldCentral->Draw();
    jetYieldCentral->SaveAs(Form("FinalRawYieldCentralPlusSystUncertainty_Jet%s_%s_%s.root", obs.Data(), fDmesonLabel.Data(),fMethodLabel.Data()));

    if (fDebug) {
      //ADVANCED - save distribution of final jet yields (summing all pT(D) bins) in a single plot
      cname = Form("cDistr_Jet%s_%s_%s", obs.Data(), fDmesonLabel.Data(), fMethodLabel.Data());
      TCanvas *cDistr = new TCanvas(cname, cname, 900, 600);
      fCanvases.push_back(cDistr);
      for (Int_t iTrial = 0; iTrial < fnMaxTrials; iTrial++) {
        for (int l=0; l < jetSpectrSBVars[iTrial]->GetNbinsX();l++) jetSpectrSBVars[iTrial]->SetBinError(l + 1, 0.0001);
        jetSpectrSBVars[iTrial]->SetMarkerColor(iTrial + 1);
        jetSpectrSBVars[iTrial]->SetLineColor(iTrial + 1);
        if (!iTrial) jetSpectrSBVars[iTrial]->Draw();
        else jetSpectrSBVars[iTrial]->Draw("same");
      }
      cDistr->SaveAs(Form("DistributionOfFinalYields_SBApproach_Jet%s_%s_AfterDbinSum.root", obs.Data(), fDmesonLabel.Data()));

      //ADVANCED - save distribution of final jet yields from each single pT(D) bin in a single plot (one per each pT(D) bin)
      for (int iDbin = 0; iDbin < fnDbins; iDbin++) {
        cname = Form("cDistr_%s_%s_%d", fDmesonLabel.Data(), fMethodLabel.Data(), iDbin);
        TCanvas *cDistr1 = new TCanvas(cname, cname, 900, 600);
        fCanvases.push_back(cDistr1);
        TH1F** hJetSpectrFromSingleDbin = new TH1F*[fnMaxTrials];
        for(Int_t iTrial = 0; iTrial < fnMaxTrials; iTrial++) {
          hJetSpectrFromSingleDbin[iTrial] = static_cast<TH1F*>(jetSpectrSBVars[0]->Clone(Form("Jet%sRawYieldDistr_Dbin%d",obs.Data(),iDbin)));
          for (int l = 0; l < hJetSpectrFromSingleDbin[iTrial]->GetNbinsX(); l++) {
            hJetSpectrFromSingleDbin[iTrial]->SetBinContent(l + 1, arrYldBinPerBin[iDbin][l][iTrial]);
            hJetSpectrFromSingleDbin[iTrial]->SetBinError(l + 1, 0.0001);
          }
          hJetSpectrFromSingleDbin[iTrial]->SetMarkerColor(iTrial + 1);
          hJetSpectrFromSingleDbin[iTrial]->SetLineColor(iTrial + 1);
          if (!iTrial) hJetSpectrFromSingleDbin[iTrial]->Draw();
          else hJetSpectrFromSingleDbin[iTrial]->Draw("same");
        }
        cDistr1->SaveAs(Form("DistributionOfFinalYields_SBApproach_Jet%s_%s_Bin%d.root", obs.Data(), fDmesonLabel.Data(), iDbin));
      }

      //ADVANCED - save averages of final jet yields from each single pT(D) bin, with their RMS, without summing them, in a single plot
      cname = Form("cDistrAllAvgs_Jet%s_%s_%s", obs.Data(), fDmesonLabel.Data(), fMethodLabel.Data());
      TCanvas *cDistr2 = new TCanvas(cname, cname, 900, 600);
      fCanvases.push_back(cDistr2);

      jetYieldCentral->SetLineWidth(3);
      jetYieldCentral->Draw();

      TH1F** hJetSpectrFromSingleDbin_Avg = new TH1F*[fnDbins];

      for (Int_t iDbin = 0; iDbin < fnDbins; iDbin++) {
        hJetSpectrFromSingleDbin_Avg[iDbin] = static_cast<TH1F*>(jetSpectrSBVars[0]->Clone(Form("Jet%sRawYieldAvgDistr_%d",obs.Data(),iDbin)));

        for(Int_t iJetbin = 0; iJetbin < nJetBins; iJetbin++) { //loop on jet spectrum pT bins

          Double_t mean = TMath::Mean(fnMaxTrials, arrYldBinPerBin[iDbin][iJetbin]);
          Double_t rms = TMath::RMS(fnMaxTrials, arrYldBinPerBin[iDbin][iJetbin]);

          jetYieldUnc->SetBinContent(iJetbin + 1, rms);
          hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinContent(iJetbin + 1, mean);
          hJetSpectrFromSingleDbin_Avg[iDbin]->SetBinError(iJetbin + 1, rms);
        }

        hJetSpectrFromSingleDbin_Avg[iDbin]->SetMarkerColor(iDbin + 1);
        hJetSpectrFromSingleDbin_Avg[iDbin]->SetLineColor(iDbin + 1);
        hJetSpectrFromSingleDbin_Avg[iDbin]->Draw("same");
      }

      TLegend* leg = new TLegend(0.1, 0.7, 0.48, 0.9);
      leg->AddEntry(jetYieldCentral, "Average after pT(D) bin sum", "pl");
      for (Int_t iDbin = 0; iDbin < fnDbins; iDbin++) leg->AddEntry(hJetSpectrFromSingleDbin_Avg[iDbin], Form("pt(D) %1.1f - %1.1f", fDbinpTedges[iDbin], fDbinpTedges[iDbin+1]), "pl");
      leg->Draw();

      cDistr2->SaveAs(Form("AverageOfFinalYields_SBApproach_Jet%s_%s_AllDBins.root", obs.Data(), fDmesonLabel.Data()));

    } //end of advanced plots
  }
  SBResults result = {kTRUE, jetYieldCentral, jetYieldUnc, jetSpectrSBVars, jetSpectrSBDef, jetBinYieldDistribution};
  return result;
}

/**
 * Fit the reflection distribution (D0 analysis)
 * @param[in] nPtBins Number of pt bins
 * @param[in] inputfile Input file name
 * @param[in] fitType Type of the fit
 */
void AliDJetRawYieldUncertainty::FitReflDistr(Int_t nPtBins, TString inputfile, TString fitType)
{
  TH1::AddDirectory(kFALSE);

  Printf("Fitting reflections.");
  TFile fReflections(inputfile.Data());
  TString inputfileNoExt = inputfile.ReplaceAll(".root", "");
  TFile fFitReflection(Form("%s_fitted_%s.root", inputfileNoExt.Data(), fitType.Data()), "recreate");

  for (Int_t i = 0; i < nPtBins; i++) {
    TH1 *hSignMC = dynamic_cast<TH1*>(fReflections.Get(Form("histSgn_%d", i)));
    if (hSignMC) {
      fFitReflection.cd();
      hSignMC->Write(Form("histSgn_%d", i));
    }
  }

  TCanvas *cy = new TCanvas(Form("%s_%s_fitCanv", fitType.Data(), inputfileNoExt.Data()), Form("%s_%s_fitCanv", fitType.Data(), inputfileNoExt.Data()));
  cy->Divide(4,3);
  TCanvas *cy2 = new TCanvas(Form("%s_%s_fitCanv2", fitType.Data(), inputfileNoExt.Data()), Form("%s_%s_fitCanv2", fitType.Data(), inputfileNoExt.Data()));
  cy2->Divide(4,3);
  TCanvas *cyRatio = new TCanvas(Form("%s_%s_fitCanvRatio", fitType.Data(), inputfileNoExt.Data()), Form("%s_%s_fitCanvRatio", fitType.Data(), inputfileNoExt.Data()));
  cyRatio->Divide(4,3);

  for (Int_t iBin = 0; iBin < nPtBins; iBin++) {
    Printf("Bin %d", iBin);
    TH1 *hfitRefl= dynamic_cast<TH1*>(fReflections.Get(Form("histRfl_%d", iBin)));
    hfitRefl->SetName(Form("histoRfl_%d",iBin));
    hfitRefl->SetMarkerStyle(1);
    hfitRefl->SetLineStyle(1);
    hfitRefl->Rebin(2);
    cy->cd(iBin+1);
    hfitRefl->Draw();
    hfitRefl->Sumw2();

    TF1 *finput = nullptr;
    if (fitType == "DoubleGaus"){//DoubleGaus
      finput= new TF1 ("finput","[0]/(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]/( TMath::Sqrt(2.*TMath::Pi())*[5])*TMath::Exp(-(x-[4])*(x-[4])/(2.*[5]*[5]))", hfitRefl->GetXaxis()->GetXmin(), hfitRefl->GetXaxis()->GetXmax());
      finput->SetParameter(0,30);
      finput->SetParameter(1,1.85);
      finput->SetParameter(2,0.2);
      finput->SetParameter(3,0.0001);
      finput->SetParameter(4,1.9);
      finput->SetParameter(5,0.05);

    }
    else if (fitType == "pol3"){
      finput= new TF1 ("finput", "pol3");
      finput->SetParameter(0,-2.5);
      finput->SetParameter(1,1.5);
      finput->SetParameter(2,0.5);
      finput->SetParameter(3,-0.4);
    }
    else if (fitType == "pol6"){
      finput= new TF1 ("finput", "pol6");
      finput->SetParameter(0,-2.5);
      finput->SetParameter(1,1.5);
      finput->SetParameter(2,0.5);
      finput->SetParameter(3,-0.4);
      finput->SetParameter(4,0.02);
      finput->SetParameter(5,0.005);
      finput->SetParameter(6,-0.005);
    }
    else if (fitType == "gaus"){
      finput= new TF1 ("finput", "gaus");
      finput->SetParameter(0,30);
      finput->SetParameter(1,1.85);
      finput->SetParameter(2,0.2);
    }

    hfitRefl->Fit("finput", "MLFI");

    Printf("Fit done");

    TF1 *fFitRefl = hfitRefl->GetFunction("finput");

    TH1 *hFitReflNewTemp = static_cast<TH1*>(hfitRefl->Clone(Form("histRflFitted%s_ptBin%d", fitType.Data(), iBin)));
    TH1 *ratio = static_cast<TH1*>(hfitRefl->Clone(Form("ratioRelDistr_%s_bin%d", fitType.Data(), iBin)));
    for (Int_t iBin2 = 1; iBin2 <= hfitRefl->GetNbinsX(); iBin2++){
      hFitReflNewTemp->SetBinContent(iBin2, 0.);
      ratio->SetBinContent(iBin2, 0.);
      hFitReflNewTemp->SetBinContent(iBin2, fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2)));
      ratio->SetBinContent(iBin2, (hfitRefl->GetBinContent(iBin2) / fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2))));
      ratio->SetBinError(iBin2, (hfitRefl->GetBinError(iBin2) / fFitRefl->Eval(hfitRefl->GetBinCenter(iBin2))));
    }

    cy2->cd(iBin + 1);
    hFitReflNewTemp->Draw();

    cyRatio->cd(iBin+1);
    ratio->GetYaxis()->SetRangeUser(-1.5, 3.);
    ratio->SetMarkerStyle(20);
    ratio->Fit("pol0", "FM");
    ratio->Draw("p");
    gPad->Update();

    fFitReflection.cd();
    hfitRefl->Write();
    hFitReflNewTemp->Write();
    ratio->Write();
  }
  fReflections.Close();
  fFitReflection.Close();
}

/**
 * Set the D meson pt bins
 * @param[in] nbins Number of pt bins
 * @param[in] ptedges Edges of the pt bins
 */
void AliDJetRawYieldUncertainty::SetDmesonPtBins(Int_t nbins, Double_t* ptedges)
{
  fnDbins = nbins;
  if (fDbinpTedges) {
    delete[] fDbinpTedges;
    fDbinpTedges = nullptr;
  }
  if (nbins == 0) return;
  fDbinpTedges = new Double_t[fnDbins + 1];
  memcpy(fDbinpTedges, ptedges, sizeof(Double_t) * (fnDbins + 1));
}

/**
 * Set the jet pt bins
 * @param[in] nbins Number of pt bins
 * @param[in] ptedges Edges of the pt bins
 */
void AliDJetRawYieldUncertainty::SetJetPtBins(Int_t nbins, Double_t* ptedges)
{
  fnJetPtbins = nbins;
  if (fJetPtBinEdges) {
    delete[] fJetPtBinEdges;
    fJetPtBinEdges = nullptr;
  }
  if (nbins == 0) return;
  fJetPtBinEdges = new Double_t[fnJetPtbins + 1];
  memcpy(fJetPtBinEdges, ptedges, sizeof(Double_t) * (fnJetPtbins + 1));
}

/**
 * Set the jet z bins
 * @param[in] nbins Number of z bins
 * @param[in] zedges Edges of the z bins
 */
void AliDJetRawYieldUncertainty::SetJetzBins(Int_t nbins, Double_t* zedges)
{
  fnJetzbins = nbins;
  if (fJetzBinEdges) {
    delete[] fJetzBinEdges;
    fJetzBinEdges = nullptr;
  }
  if (nbins == 0) return;
  fJetzBinEdges = new Double_t[fnJetzbins + 1];
  memcpy(fJetzBinEdges, zedges, sizeof(Double_t) * (fnJetzbins + 1));
}

/**
 * Set the fixed sigma values in bins of D meson pt
 * @param[in] sigmafix Values of the fixed sigmas
 */
void AliDJetRawYieldUncertainty::SetSigmaToFixDPtBins(Double_t* sigmafix)
{
  if (fSigmaToFixDPtBins) {
    delete[] fSigmaToFixDPtBins;
    fSigmaToFixDPtBins = nullptr;
  }
  if (fnDbins == 0) return;
  fSigmaToFixDPtBins = new Double_t[fnDbins];
  memcpy(fSigmaToFixDPtBins, sigmafix, sizeof(Double_t) * fnDbins);
}

/**
 * Set the fixed sigma values in bins of jet pt
 * @param[in] sigmafix Values of the fixed sigmas
 */
void AliDJetRawYieldUncertainty::SetSigmaToFixJetPtBins(Double_t* sigmafix)
{
  if (fSigmaToFixJetPtBins) {
    delete[] fSigmaToFixJetPtBins;
    fSigmaToFixJetPtBins = nullptr;
  }
  if (fnJetPtbins == 0) return;
  fSigmaToFixJetPtBins = new Double_t[fnJetPtbins];
  memcpy(fSigmaToFixJetPtBins, sigmafix, sizeof(Double_t) * fnJetPtbins);
}

/**
 * Set the efficiency values in bins of D meson pt
 * @param[in] sigmafix Values of the efficiency
 */
void AliDJetRawYieldUncertainty::SetDmesonEfficiency(Double_t* effvalues)
{
  if (fDEffValues) {
    delete[] fDEffValues;
    fDEffValues = nullptr;
  }
  if (fnDbins == 0) return;
  fDEffValues = new Double_t[fnDbins];
  memcpy(fDEffValues, effvalues, sizeof(Double_t) * fnDbins);
}

/**
 * Enable/disable the various mean/sigma cases
 * @param[in] cases True/False for each case
 */
void AliDJetRawYieldUncertainty::SetMeanSigmaVariations(Bool_t* cases)
{
  memcpy(fMeanSigmaVar, cases, sizeof(Bool_t) * fgkNSigmaVar);
}

/**
 * Enable/disable the various background function cases
 * @param[in] cases True/False for each case
 */
void AliDJetRawYieldUncertainty::SetBkgVariations(Bool_t* cases)
{
  memcpy(fBkgVar, cases, sizeof(Bool_t) * fgkNBkgVar);
}

/**
 * Set the rebin steps
 * @param[in] nsteps Number of steps
 * @param[in] cases Rebin values for each step
 */
void AliDJetRawYieldUncertainty::SetRebinSteps(Int_t nsteps, Int_t* cases)
{
  fnRebinSteps = nsteps;
  if (fRebinSteps) {
    delete[] fRebinSteps;
    fRebinSteps = nullptr;
  }
  if (fnRebinSteps == 0) return;
  fRebinSteps = new Int_t[fnRebinSteps];
  memcpy(fRebinSteps, cases, sizeof(Int_t) * fnRebinSteps);
}

/**
 * Set the min mass steps
 * @param[in] nsteps Number of steps
 * @param[in] cases Min mass values for each step
 */
void AliDJetRawYieldUncertainty::SetMinMassSteps(Int_t nsteps, Double_t* cases)
{
  fnMinMassSteps = nsteps;
  if (fMinMassSteps) {
    delete[] fMinMassSteps;
    fMinMassSteps = nullptr;
  }
  if (fnMinMassSteps == 0) return;
  fMinMassSteps = new Double_t[fnMinMassSteps];
  memcpy(fMinMassSteps, cases, sizeof(Double_t) * fnMinMassSteps);
}

/**
 * Set the max mass steps
 * @param[in] nsteps Number of steps
 * @param[in] cases Max mass values for each step
 */
void AliDJetRawYieldUncertainty::SetMaxMassSteps(Int_t nsteps, Double_t* cases)
{
  fnMaxMassSteps = nsteps;
  if (fMaxMassSteps) {
    delete[] fMaxMassSteps;
    fMaxMassSteps = nullptr;
  }
  if (fnMaxMassSteps == 0) return;
  fMaxMassSteps = new Double_t[fnMaxMassSteps];
  memcpy(fMaxMassSteps, cases, sizeof(Double_t) * fnMaxMassSteps);
}

/**
 * Set the bin counting steps
 * @param[in] nsteps Number of steps
 * @param[in] cases Bin counting for each step
 */
void AliDJetRawYieldUncertainty::SetSigmaBinCounting(Int_t nsteps, Double_t* cases)
{
  fnSigmaBC = nsteps;
  if (fSigmaBC) {
    delete[] fSigmaBC;
    fSigmaBC = nullptr;
  }
  if (fnSigmaBC == 0) return;
  fSigmaBC = new Double_t[fnSigmaBC];
  memcpy(fSigmaBC, cases, sizeof(Double_t) * fnSigmaBC);
}

/**
 * Set the variation mask
 * @param[in] ncases Number of cases
 * @param[in] cases Enable/Disable each case
 */
void AliDJetRawYieldUncertainty::SetMaskOfVariations(Int_t ncases, Bool_t* cases)
{
  fnMask = ncases;
  if (fMask) {
    delete[] fMask;
    fMask = nullptr;
  }
  if (fnMask == 0) return;
  fMask = new Bool_t[fnMask];
  memcpy(fMask, cases, sizeof(Bool_t) * fnMask);
}

/**
 * Release memory.
 */
void AliDJetRawYieldUncertainty::ClearObjects()
{
  if (fDbinpTedges) {
    delete[] fDbinpTedges;
    fDbinpTedges = nullptr;
  }
  if (fJetPtBinEdges) {
    delete[] fJetPtBinEdges;
    fJetPtBinEdges = nullptr;
  }
  if (fJetzBinEdges) {
    delete[] fJetzBinEdges;
    fJetzBinEdges = nullptr;
  }
  if (fDEffValues) {
    delete[] fDEffValues;
    fDEffValues = nullptr;
  }
  if (fSigmaToFixDPtBins) {
    delete[] fSigmaToFixDPtBins;
    fSigmaToFixDPtBins = nullptr;
  }
  if (fSigmaToFixJetPtBins) {
    delete[] fSigmaToFixJetPtBins;
    fSigmaToFixJetPtBins = nullptr;
  }
  if (fRebinSteps) {
    delete[] fRebinSteps;
    fRebinSteps = nullptr;
  }
  if (fMinMassSteps) {
    delete[] fMinMassSteps;
    fMinMassSteps = nullptr;
  }
  if (fMaxMassSteps) {
    delete[] fMaxMassSteps;
    fMaxMassSteps = nullptr;
  }
  if (fSigmaBC) {
    delete[] fSigmaBC;
    fSigmaBC = nullptr;
  }
  if (fMask) {
    delete[] fMask;
    fMask = nullptr;
  }

  fnDbins = 0;
  fnJetPtbins = 0;
  fnJetzbins = 0;
  fnRebinSteps = 0;
  fnMinMassSteps = 0;
  fnMaxMassSteps = 0;
  fnSigmaBC = 0;
  fnMask = 0;
}
