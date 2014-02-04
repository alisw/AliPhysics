#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include "AliCFContainer.h"
#include "AliCFEffGrid.h"
#include "AliCFDataGrid.h"
#include "AliPID.h"

#include <iostream>

#include "THnSparseDefinitions.h"

enum type { kTrackPt = 0, kZ = 1, kXi = 2, kNtypes = 3 };


const Int_t nPtBinsType2 = 44;
const Double_t binsPtType2[nPtBinsType2+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
          0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4,
          1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 3.8,
          4.5, 5.5, 6.5, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0,
          28.0, 32.0, 36.0, 45.0, 50.0 };


//___________________________________________________________________
const Double_t* GetBins(Int_t type, Int_t& nBins)
{
  if (type == 2) {
    nBins = nPtBinsType2;
    
    return binsPtType2;
  }
  
  return 0x0;
}


//___________________________________________________________________
Double_t trackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1. - 0.129758 * TMath::Exp(-pTmc * 0.679612));
}


//___________________________________________________________________
Double_t trackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093 * pTmc), 1.);
}


//___________________________________________________________________
Bool_t geantFlukaCorrection(AliCFContainer* data, Int_t genStepToDownscale)
{
  // Issue: GEANT/FLUKA correction factor is for MC_pT. 
  // Finally the effeciency should be DIVIDED by this correction factor.
  // To include resolution effects, it is therefore the best way to just
  // multiply the generated step with the correction factor.
  
  if (!data) {
    printf("No CFContainer for GEANT/FLUKA correction!\n");
    return kFALSE;
  }
  
  const Int_t iPt     = data->GetVar("P_{T} (GeV/c)");
  const Int_t iMCid   = data->GetVar("MC ID");
  const Int_t iCharge = data->GetVar("Charge (e_{0})");
  
  if (iPt < 0 || iMCid < 0 || iCharge < 0) {
    printf("Data axis for GEANT/FLUKA correction not found!\n");
    return kFALSE;
  }
  
  if (!data->GetGrid(genStepToDownscale)) {
    printf("Step for downscaling (GEANT/FLUKA) not found!\n");
    return kFALSE;
  }
  
  const Int_t nDim = data->GetNVar();
  Int_t coord[nDim];
  Double_t binCenterCoord[nDim];
  
  Long64_t nBinsGrid = data->GetGrid(genStepToDownscale)->GetGrid()->GetNbins();
  
  for (Long64_t iBin = 0; iBin < nBinsGrid; iBin++) {
    Double_t binContent = data->GetGrid(genStepToDownscale)->GetGrid()->GetBinContent(iBin, coord);
    Double_t binError  = data->GetGrid(genStepToDownscale)->GetGrid()->GetBinError(iBin);
    
    for (Int_t iDim = 0; iDim < nDim; iDim++) 
      binCenterCoord[iDim] = data->GetBinCenter(iDim, coord[iDim]);

    if (binCenterCoord[iCharge] < 0) {
      Double_t corrFactor = 1.;
      
      if (binCenterCoord[iMCid] - 0.5 == AliPID::kProton) 
        corrFactor = trackingPtGeantFlukaCorrectionPrMinus(binCenterCoord[iPt]);
      else if (binCenterCoord[iMCid] - 0.5 == AliPID::kKaon)
        corrFactor = trackingPtGeantFlukaCorrectionKaMinus(binCenterCoord[iPt]);
      else
        continue;
      
      data->GetGrid(genStepToDownscale)->GetGrid()->SetBinContent(iBin, binContent * corrFactor);
      data->GetGrid(genStepToDownscale)->GetGrid()->SetBinError(iBin, binError * corrFactor);
    }
  }
  
  return kTRUE;
}



//___________________________________________________________________
void setupHist(TH1* h, TString histName, TString histTitle, TString xAxisTitle, TString yAxisTitle, Int_t color)
{
  if (histName != "")
    h->SetName(histName.Data());
  h->SetTitle(histTitle.Data());
  
  if (xAxisTitle != "")
    h->GetXaxis()->SetTitle(xAxisTitle.Data());
  if (yAxisTitle != "")
    h->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  h->SetMarkerStyle(24);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  
  h->SetStats(kFALSE);
}


//____________________________________________________________________________________________________________________
void undoJetNormalisationHist2D(TH2* hData, TH2* hNumJets, const Int_t lowerCentralityBinLimit, const Int_t upperCentralityBinLimit)
{
  // Undo normalisation to 1/numJets. NOTE: jetPt binning of hData and hNumJets assumed to be the same!
  
  if (!hData || !hNumJets)
    return;
  
  for (Int_t binJetPt = 0; binJetPt <= hData->GetNbinsY() + 1; binJetPt++) {
    const Double_t numJets = hNumJets->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, binJetPt, binJetPt);
    Bool_t noJets = numJets < 1e-13;
    
    for (Int_t binObs = 0; binObs <= hData->GetNbinsX() + 1; binObs++) {
      if (noJets) 
        continue;

      hData->SetBinContent(binObs, binJetPt, hData->GetBinContent(binObs, binJetPt) * numJets);
      hData->SetBinError(binObs, binJetPt, hData->GetBinError(binObs, binJetPt) * numJets);
    }
  }
}


//___________________________________________________________________
void normaliseHist(TH1* h, Double_t scale = 1.0)
{
  if (h->GetSumw2N() <= 0)
    h->Sumw2();
  
  h->Scale(scale);
  
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t normFactor = h->GetBinWidth(i);
    h->SetBinContent(i, h->GetBinContent(i) / normFactor);
    h->SetBinError(i, h->GetBinError(i) / normFactor);
  }
}


//___________________________________________________________________
void multiplyHistsDifferentBinning(TH1* h1, const TH1* h2, Double_t c1, Double_t c2, Bool_t ignoreErrorOfSecondHist = kFALSE)
{
   Int_t nbinsx = h1->GetNbinsX();
   Int_t nbinsy = h1->GetNbinsY();
   Int_t nbinsz = h1->GetNbinsZ();
   
   if (h1->GetDimension() < 2)
     nbinsy -= 1;
   if (h1->GetDimension() < 3)
     nbinsz -= 1;
   
   if (h1->GetSumw2N() == 0)
     h1->Sumw2();
   
   Int_t bin, bin2, binx, biny, binz;
   Double_t b1,b2,w,d1,d2;
   d1 = c1*c1;
   d2 = c2*c2;
   
   // Loop over bins (including underflows/overflows)
   for (binz = 0; binz <= nbinsz + 1; binz++) {
      for (biny = 0; biny <= nbinsy + 1; biny++) {
         for (binx = 0; binx <= nbinsx + 1;binx++) {
            bin = binx + (nbinsx + 2) * (biny + (nbinsy + 2) * binz);
            bin2 = h2->FindFixBin(h1->GetXaxis()->GetBinCenter(binx), h1->GetYaxis()->GetBinCenter(biny),
                                  h1->GetZaxis()->GetBinCenter(binz));
            b1  = h1->GetBinContent(bin);
            b2  = h2->GetBinContent(bin2);
            w   = (c1*b1)*(c2*b2);
            h1->SetBinContent(bin, w);
            Double_t e1 = h1->GetBinError(bin);
            Double_t e2 = ignoreErrorOfSecondHist ? 0. : h2->GetBinError(bin2);
            h1->SetBinError(bin, TMath::Sqrt(d1*d2*(e1*e1*b2*b2 + e2*e2*b1*b1)));
         }
      }
   }
}


//___________________________________________________________________
void divideHistsDifferentBinning(TH1* h1, const TH1* h2, Double_t c1, Double_t c2, Bool_t ignoreErrorOfSecondHist = kFALSE)
{
   Int_t nbinsx = h1->GetNbinsX();
   Int_t nbinsy = h1->GetNbinsY();
   Int_t nbinsz = h1->GetNbinsZ();
   
   if (h1->GetDimension() < 2)
     nbinsy -= 1;
   if (h1->GetDimension() < 3)
     nbinsz -= 1;
   
   if (h1->GetSumw2N() == 0)
     h1->Sumw2();
   
   Int_t bin, bin2, binx, biny, binz;
   Double_t b1,b2,w,d1,d2;
   d1 = c1*c1;
   d2 = c2*c2;
   
   // Loop over bins (including underflows/overflows)
   for (binz = 0; binz <= nbinsz + 1; binz++) {
      for (biny = 0; biny <= nbinsy + 1; biny++) {
         for (binx = 0; binx <= nbinsx + 1;binx++) {
            bin = binx + (nbinsx + 2) * (biny + (nbinsy + 2) * binz);
            bin2 = h2->FindFixBin(h1->GetXaxis()->GetBinCenter(binx), h1->GetYaxis()->GetBinCenter(biny),
                                  h1->GetZaxis()->GetBinCenter(binz));
            b1  = h1->GetBinContent(bin);
            b2  = h2->GetBinContent(bin2);
            if (b2)
              w = (c1*b1)/(c2*b2);
            else 
              w = 0;
            h1->SetBinContent(bin, w);
            
            if (!b2) {
              h1->SetBinError(bin, 0);
              continue;
            }
            
            Double_t b22 = b2*b2*d2;
            Double_t e1 = h1->GetBinError(bin);
            Double_t e2 = ignoreErrorOfSecondHist ? 0. : h2->GetBinError(bin2);
            h1->SetBinError(bin, TMath::Sqrt(d1*d2*(e1*e1*b2*b2 + e2*e2*b1*b1)/(b22*b22)));
         }
      }
   }
}


//___________________________________________________________________
// Efficiency for inclusive spectra vs. pT (or also jets, but without using the generation)
// E.g. a 'calcEfficiency.C+("finalCuts/MC_pp/7TeV/LHC10f6a/corrected/finalisedSplines/analytical/efficiency_noClCut_preliminary/bhess_PID_Jets_Inclusive_PureGauss_efficiency.root", "finalCuts/pp/7TeV/LHC10e.pass2/corrected/finalisedSplines/finalMapsAndTail/Jets/noCutOn_ncl_or_liav/outputSystematicsTotal_SummedSystematicErrors__2013_10_21.root", kTRUE, kFALSE, 0, -2, -2, -1, -1, 8, 2)' -b -q
Int_t calcEfficiency(TString pathNameEfficiency, TString pathNameData, Bool_t correctGeantFluka, Bool_t scaleStrangeness,
                     Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                     Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                     Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/,
                     Double_t constantCorrectionAbovePtThreshold,
                     Int_t rebinEfficiency)
{
  TString pathData = pathNameData;
  pathData.Replace(pathData.Last('/'), pathData.Length(), "");
  
  TFile* fileEff = TFile::Open(pathNameEfficiency.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", pathNameEfficiency.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
  }
  
  const Int_t iMCid   = data->GetVar("MC ID");
  const Int_t iPt     = data->GetVar("P_{T} (GeV/c)");
  const Int_t iEta    = data->GetVar("#eta");
  const Int_t iCharge = data->GetVar("Charge (e_{0})");
  const Int_t iMult   = data->GetVar("Centrality Percentile");
  const Int_t iJetPt  = data->GetVar("P_{T}^{jet} (GeV/c)");
  //const Int_t iZ      = data->GetVar("z = P_{T}^{track} / P_{T}^{jet}");
  //const Int_t iXi     = data->GetVar("#xi = ln(P_{T}^{jet} / P_{T}^{track})");
  
  TFile* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TH1D* hYield[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldSysError[AliPID::kSPECIES] = { 0x0, };
  TGraphAsymmErrors* gYieldSysError[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimYield[AliPID::kSPECIES] = { 0x0, };
  Int_t numMCgenPrimYieldHistsFound = 0;
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    TString speciesName = AliPID::ParticleName(species);
    TString firstLetter = speciesName(0);
    firstLetter.ToUpper();
    speciesName.Replace(0, 1, firstLetter.Data());
    TString histName = Form("hYield%ss", speciesName.Data());
    hYield[species] = (TH1D*)fileData->Get(histName.Data());
    if (!hYield[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    hYield[species]->SetFillStyle(0);
    
    TString graphName = Form("systematicErrorYields_%s", AliPID::ParticleName(species));
    gYieldSysError[species] = (TGraphAsymmErrors*)fileData->Get(graphName.Data());
    
    // In case of MC also retrieve the MC truth generated yields
    TString histNameMCgenYields = Form("hMCgenYieldsPrimSpecies_%s", AliPID::ParticleShortName(species));
    hMCgenPrimYield[species] = (TH1D*)fileData->Get(histNameMCgenYields.Data());
    if (hMCgenPrimYield[species])
      numMCgenPrimYieldHistsFound++;
  }
  
  if (numMCgenPrimYieldHistsFound > 0 && numMCgenPrimYieldHistsFound != AliPID::kSPECIES) {
    printf("Error: Unable to retrieve all MC generated prim yield histos! Got %d.\n", numMCgenPrimYieldHistsFound);
    return -1;
  }
  
  TCanvas* cFractions = (TCanvas*)fileData->Get("cFractionsWithTotalSystematicError");
  if (!cFractions)
    cFractions = (TCanvas*)fileData->Get("cFractions");
  
  
  // Convert graphs with systematic errors into histograms (assume symmetric error, which is currently true)
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (gYieldSysError[species]) {
      hYieldSysError[species] = new TH1D(*hYield[species]);
      hYieldSysError[species]->SetName(Form("%s_sysError", hYieldSysError[species]->GetName()));
      hYieldSysError[species]->SetFillStyle(0);
      
      for (Int_t binX = 1; binX <= hYieldSysError[species]->GetNbinsX(); binX++) {
        hYieldSysError[species]->SetBinContent(binX, gYieldSysError[species]->GetY()[binX - 1]);
        hYieldSysError[species]->SetBinError(binX, gYieldSysError[species]->GetErrorY(binX - 1));
      }
    }
  }
  
  // Take pT binning from pion yield (binning for all species the same) and create a new AliCFContainer with this new binning for pT
  TH1D hDummy(*hYield[AliPID::kPion]);
  hDummy.SetName("hDummy");
  TAxis* axis = 0x0;
  
  if (rebinEfficiency > 1) {
    Int_t nBinsNew = 0;
    const Double_t* newBins = GetBins(rebinEfficiency, nBinsNew);
    
    hDummy.SetBins(nBinsNew, newBins);
    
    axis = hDummy.GetXaxis();
  
    //axis = hDummy.Rebin(rebinEfficiency, "", 0)->GetXaxis();
  }
  else
    axis = hDummy.GetXaxis();

  const TArrayD* binsPtCurrent = axis->GetXbins();
  TArrayD* binsPtNew = new TArrayD(*binsPtCurrent);
  
  const Int_t nEffDims = data->GetNVar();
  Int_t nEffBins[nEffDims];
  
  for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
    if (iDim == iPt)
      nEffBins[iDim] = axis->GetNbins();
    else 
      nEffBins[iDim] = data->GetNBins(iDim);
  }
  
  
  // Just make one large pT bin above some threshold, if desired
  if (binsPtNew->fN != 0 && constantCorrectionAbovePtThreshold > 0) {
    for (Int_t iBin = 0; iBin < nEffBins[iPt]; iBin++) {
      // Find the first bin edged really larger than the threshold.
      // If the bin edge before equals the threshold, just set the
      // current bin edge to the right end of the spectrum -> Done.
      // If the bin edge before is different, set the bin edge to the
      // threshold
      if (binsPtNew->fArray[iBin] > constantCorrectionAbovePtThreshold) {
        if (binsPtNew->fArray[iBin - 1] == constantCorrectionAbovePtThreshold) {
          binsPtNew->fArray[iBin] = binsPtNew->fArray[nEffBins[iPt]];
          nEffBins[iPt] = iBin;
          break;
        }
        else {
          binsPtNew->fArray[iBin] = constantCorrectionAbovePtThreshold;
        }
      }
    }
  }
  
  
  AliCFContainer *dataRebinned = new AliCFContainer(Form("%s_rebinned", data->GetName()), Form("%s (rebinned)", data->GetTitle()),
                                                    data->GetNStep(), nEffDims, nEffBins);
  
  for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
    dataRebinned->SetVarTitle(iDim, data->GetVarTitle(iDim));
    
    if (iDim == iPt) {
      if (binsPtNew->fN == 0)
        dataRebinned->SetBinLimits(iDim, axis->GetXmin(), axis->GetXmax());
      else
        dataRebinned->SetBinLimits(iDim, binsPtNew->fArray);
    }
    else {
      dataRebinned->SetBinLimits(iDim, data->GetBinLimits(iDim));
    }
  }
  
  for (Int_t iStep = 0; iStep < data->GetNStep(); iStep++)
    dataRebinned->SetStepTitle(iStep, data->GetStepTitle(iStep));
  
  Int_t coord[nEffDims];
  Double_t binCenterCoord[nEffDims];
  
  // Fill content from old grid into the new grid with proper binning
  for (Int_t iStep = 0; iStep < data->GetNStep(); iStep++) {
    Long64_t nBinsGrid = data->GetGrid(iStep)->GetGrid()->GetNbins();
    
    for (Long64_t iBin = 0; iBin < nBinsGrid; iBin++) {
      Double_t binContent = data->GetGrid(iStep)->GetGrid()->GetBinContent(iBin, coord);
      Double_t binError2  = data->GetGrid(iStep)->GetGrid()->GetBinError2(iBin);
      
      for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
        binCenterCoord[iDim] = data->GetBinCenter(iDim, coord[iDim]);
      }

      Long64_t iBinRebinned = dataRebinned->GetGrid(iStep)->GetGrid()->GetBin(binCenterCoord);
      dataRebinned->GetGrid(iStep)->GetGrid()->AddBinContent(iBinRebinned, binContent);
      dataRebinned->GetGrid(iStep)->GetGrid()->AddBinError2(iBinRebinned, binError2);
    }
  }
  
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = dataRebinned->GetAxis(iMult, 0)->FindBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = dataRebinned->GetAxis(iMult, 0)->FindBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= dataRebinned->GetAxis(iMult, 0)->GetNbins()) {
      actualLowerCentrality = dataRebinned->GetAxis(iMult, 0)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = dataRebinned->GetAxis(iMult, 0)->GetBinUpEdge(upperCentralityBinLimit);

      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis) {
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
    dataRebinned->SetRangeUser(iMult, lowerCentralityBinLimit, upperCentralityBinLimit, kTRUE);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -1;
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 &&
        upperJetPtBinLimit <= dataRebinned->GetAxis(iJetPt, 0)->GetNbins()) {
      actualLowerJetPt = dataRebinned->GetAxis(iJetPt, 0)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = dataRebinned->GetAxis(iJetPt, 0)->GetBinUpEdge(upperJetPtBinLimit);

      restrictJetPtAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested jet pT range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "jet pT: ";
  if (restrictJetPtAxis) {
    std::cout << actualLowerJetPt << " - " << actualUpperJetPt << std::endl;
    dataRebinned->SetRangeUser(iJetPt, lowerJetPtBinLimit, upperJetPtBinLimit, kTRUE);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  // If desired, restrict charge axis
  std::cout << "Charge selection (efficiency): ";
  if (chargeMode == kAllCharged)
    std::cout << "All charged particles" << std::endl;
  else if (chargeMode == kNegCharge)
    std::cout << "Negative particles only" << std::endl;
  else if (chargeMode == kPosCharge)
    std::cout << "Positive particles only" << std::endl;
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  const Bool_t restrictCharge = (chargeMode != kAllCharged);
  
  Int_t lowerChargeBinLimit = -1;
  Int_t upperChargeBinLimit = -2;
  Double_t actualLowerCharge = -999;
  Double_t actualUpperCharge = -999;
  
  if (restrictCharge) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(-1. + 0.001);
      upperChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(0. + 0.001);
      upperChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimit <= upperChargeBinLimit && lowerChargeBinLimit >= 1
        && upperChargeBinLimit <= dataRebinned->GetAxis(iCharge, 0)->GetNbins()) {
      actualLowerCharge = dataRebinned->GetAxis(iCharge, 0)->GetBinLowEdge(lowerChargeBinLimit);
      actualUpperCharge = dataRebinned->GetAxis(iCharge, 0)->GetBinUpEdge(upperChargeBinLimit);
      
      std::cout << "Charge range (efficiency): " << actualLowerCharge << " - " << actualUpperCharge << std::endl;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range (efficiency) out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    dataRebinned->SetRangeUser(iCharge, lowerChargeBinLimit, upperChargeBinLimit, kTRUE);
    data->SetRangeUser(iCharge, lowerChargeBinLimit, upperChargeBinLimit, kTRUE);
  }
  
  std::cout << std::endl;
  
  
  // If jet axis is restricted (i.e. jet input), also need num jet histos for proper normalisation.
  // Load numJet histos from bhess_PID*.root file related to efficiency file.
  TH2D* hNjetsGen = 0x0;
  TH2D* hNjetsRec = 0x0;

  if (restrictJetPtAxis) {
    TString pathNameDataMC = pathNameEfficiency;
    pathNameDataMC.ReplaceAll("_efficiency", "");
    
    TFile* fDataMC = TFile::Open(pathNameDataMC.Data());
    if (!fDataMC)  {
      std::cout << std::endl;
      std::cout << "Failed to open file \"" << pathNameData.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
      
      return -1;
    }
    
    TString listName = pathNameDataMC;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
      
    TObjArray* histList = (TObjArray*)(fDataMC->Get(listName.Data()));
    if (!histList) {
      std::cout << std::endl;
      std::cout << "Failed to load list \"" << listName.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
      return -1;
    }
    
    hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
    hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
    
    if (!hNjetsRec || !hNjetsGen) {
      std::cout << "Failed to load number of jets histos!" << std::endl;
      
      
      // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
      // period is considered; also: No multiplicity information)
      TString pathEfficiency = pathNameEfficiency;
      pathEfficiency.Replace(pathEfficiency.Last('/'), pathEfficiency.Length(), "");
      TString pathBackward = Form("%s/AnalysisResults.root", pathEfficiency.Data());
      TFile* fBackward = TFile::Open(pathBackward.Data());
      
      TString dirDataInFile = "";
      TDirectory* dirData = fBackward ? (TDirectory*)fBackward->Get(fBackward->GetListOfKeys()->At(0)->GetName()) : 0x0;
    
      TList* list = dirData ? (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName()) : 0x0;

      TH1D* hFFJetPtRec = list ? (TH1D*)list->FindObject("fh1FFJetPtRecCuts") : 0x0;
      TH1D* hFFJetPtGen = list ? (TH1D*)list->FindObject("fh1FFJetPtGen") : 0x0;
      
      if (hFFJetPtRec && hFFJetPtGen) {
        printf("***WARNING: For backward compatibility, using file \"%s\" to get number of jets. BUT: Might be wrong period and has no mult info!***\n",
          pathBackward.Data());
        
        hNjetsRec = new TH2D("fh2FFJetPtRec", "", 1, -1, 1, dataRebinned->GetAxis(iJetPt, 0)->GetNbins(),
                            dataRebinned->GetAxis(iJetPt, 0)->GetXbins()->GetArray());
        
        for (Int_t iJet = 1; iJet <= hNjetsRec->GetNbinsY(); iJet++) {
          Int_t lowerBin = hFFJetPtRec->FindBin(hNjetsRec->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
          Int_t upperBin = hFFJetPtRec->FindBin(hNjetsRec->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
          hNjetsRec->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
        }
        
        hNjetsGen = new TH2D("fh2FFJetPtGen", "", 1, -1, 1,  dataRebinned->GetAxis(iJetPt, 0)->GetNbins(),
                            dataRebinned->GetAxis(iJetPt, 0)->GetXbins()->GetArray());
        
        for (Int_t iJet = 1; iJet <= hNjetsGen->GetNbinsY(); iJet++) {
          Int_t lowerBin = hFFJetPtGen->FindBin(hNjetsGen->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
          Int_t upperBin = hFFJetPtGen->FindBin(hNjetsGen->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
          hNjetsGen->SetBinContent(1, iJet, hFFJetPtGen->Integral(lowerBin, upperBin));
        }
      }
      
      if (!hNjetsRec || ! hNjetsGen)
        return -1;
    }
  }

  // For normalisation to number of jets
  // NOTE: These numbers are for the efficiency only! The data will be normalised to its own number!!!
  const Double_t nJetsGen = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                            upperJetPtBinLimit) : 1.;
  const Double_t nJetsRec = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                            upperJetPtBinLimit) : 1.;
  
  // Secondary correction
  AliCFEffGrid* sec = new AliCFEffGrid("sec", "Secondary Contamination", *dataRebinned);
  AliCFEffGrid* secStrangeScale = new AliCFEffGrid("secStrangeScale", "Secondary Contamination with Strangeness Scaling",
                                                   *dataRebinned);
  
  // Either one can take kStepRecWithGenCutsMeasuredObs or, what I prefer, one can take
  // kStepRecWithRecCutsMeasuredObsPrimaries => The difference is only the eta cut, which is on the rec level
  // in the latter case, i.e. one corrects for eta resolution (although the effect should be very small)
  // => TESTED: There is NO difference (only 1 bin vs. pt shows a deviation by one entry), so this cut has,
  // in principle, more or less no effect
  
  // For testing with data set w/o strangeness secStrangeScale->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObs);
  secStrangeScale->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObsStrangenessScaled);
  sec->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObs);
  
  // QA plots for secondaries
  TCanvas* cSec = new TCanvas("cSec", "Secondary Contamination", 100, 10, 1200, 800);
  cSec->Divide(2, 1);
  cSec->GetPad(1)->SetLogx(kTRUE);
  cSec->cd(1);
  TH1* hSecPt = sec->Project(iPt); 
  hSecPt->SetName("hSecPt");
  hSecPt->SetStats(kFALSE);
  hSecPt->GetXaxis()->SetRangeUser(0.15, 50.);
  hSecPt->GetXaxis()->SetMoreLogLabels(kTRUE);
  hSecPt->GetXaxis()->SetNoExponent(kTRUE);
  hSecPt->GetYaxis()->SetTitle("Primary Fraction");
  hSecPt->Draw("E1");
  cSec->cd(2);
  TH1* hEtaSec = sec->Project(iEta);
  hEtaSec->SetName("hEtaSec");
  hEtaSec->SetStats(kFALSE);
  hEtaSec->GetYaxis()->SetTitle("Primary Fraction");
  hEtaSec->Draw("E1");
  TH2D* hSecID2Pt = (TH2D*)sec->Project(iPt, iMCid);
  hSecID2Pt->SetName("hSecID2Pt");
  
  // Get the secondary contamination vs. pT for each species
  TH1D* hSec[AliPID::kSPECIES];
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSec[species] = hSecID2Pt->ProjectionX(Form("hSec_%s", AliPID::ParticleShortName(species)), species + 1, species + 1, "e");
    hSec[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hSec[species]->SetLineColor(hYield[species]->GetLineColor());
    hSec[species]->SetMarkerColor(hYield[species]->GetLineColor());
    hSec[species]->SetLineWidth(2.);
    hSec[species]->GetXaxis()->SetRangeUser(0.15, 50);
    hSec[species]->GetYaxis()->SetRangeUser(0., 1.01);
    hSec[species]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hSec[species]->GetXaxis()->SetNoExponent(kTRUE);
    hSec[species]->SetStats(kFALSE);
    hSec[species]->GetYaxis()->SetTitle("Primary Fraction");
  }
  
  TCanvas* cSec2 = new TCanvas("cSec2", "Primary fraction for different species", 100, 10, 1200, 800);
  cSec2->SetGridx(1);
  cSec2->SetGridy(1);
  cSec2->SetLogx(1);
  
  hSec[0]->Draw("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    hSec[i]->Draw("E1 same");
  }
  cSec2->BuildLegend()->SetFillColor(kWhite);
  
  ClearTitleFromHistoInCanvas(cSec2);
  
  // QA plots for secondaries with strangeness scaling
  TCanvas* cSecSS = new TCanvas("cSecSS", "Secondary Contamination (strangeness scaled)", 100, 10, 1200, 800);
  cSecSS->Divide(2, 1);
  cSecSS->GetPad(1)->SetLogx(kTRUE);
  cSecSS->cd(1);
  TH1* hSecSSPt = secStrangeScale->Project(iPt); 
  hSecSSPt->SetName("hSecSSPt");
  hSecSSPt->SetStats(kFALSE);
  hSecSSPt->GetXaxis()->SetRangeUser(0.15, 50.);
  hSecSSPt->GetXaxis()->SetMoreLogLabels(kTRUE);
  hSecSSPt->GetXaxis()->SetNoExponent(kTRUE);
  hSecSSPt->GetYaxis()->SetTitle("Primary Fraction");
  hSecSSPt->Draw("E1");
  cSecSS->cd(2);
  TH1* hEtaSecSS = secStrangeScale->Project(iEta);
  hEtaSecSS->SetName("hEtaSecSS");
  hEtaSecSS->SetStats(kFALSE);
  hEtaSecSS->GetYaxis()->SetTitle("Primary Fraction");
  hEtaSecSS->Draw("E1");
  TH2D* hSecSSID2Pt = (TH2D*)secStrangeScale->Project(iPt, iMCid);
  hSecSSID2Pt->SetName("hSecSSID2Pt");
  
  // Get the secondary contamination vs. pT for each species
  TH1D* hSecSS[AliPID::kSPECIES];
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSecSS[species] = hSecSSID2Pt->ProjectionX(Form("hSecSS_%s", AliPID::ParticleShortName(species)), species + 1, species + 1, "e");
    hSecSS[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hSecSS[species]->SetLineColor(hYield[species]->GetLineColor());
    hSecSS[species]->SetMarkerColor(hYield[species]->GetLineColor());
    hSecSS[species]->SetLineWidth(2.);
    hSecSS[species]->GetXaxis()->SetRangeUser(0.15, 50);
    hSecSS[species]->GetYaxis()->SetRangeUser(0., 1.01);
    hSecSS[species]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hSecSS[species]->GetXaxis()->SetNoExponent(kTRUE);
    hSecSS[species]->SetStats(kFALSE);
    hSecSS[species]->GetYaxis()->SetTitle("Primary Fraction");
  }
  
  TCanvas* cSecSS2 = new TCanvas("cSecSS2", "Primary fraction for different species (strangeness scaled)", 100, 10, 1200, 800);
  cSecSS2->SetGridx(1);
  cSecSS2->SetGridy(1);
  cSecSS2->SetLogx(1);
  
  hSecSS[0]->Draw("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    hSecSS[i]->Draw("E1 same");
  }
  cSecSS2->BuildLegend()->SetFillColor(kWhite);
  
  ClearTitleFromHistoInCanvas(cSecSS2);
  
  
  // Secondary correction for to-pi-ratios
  TH1D* hSecToPiRatio[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecSSToPiRatio[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kPion)
      continue; // Do not consider pion-to-pion ratio
      
    hSecToPiRatio[species] = new TH1D(*hSec[species]);
    hSecToPiRatio[species]->Reset();
    hSecToPiRatio[species]->SetName(Form("hSecSSToPionRatio_%s", AliPID::ParticleShortName(species)));
    hSecToPiRatio[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
    hSecToPiRatio[species]->SetLineColor(getLineColorAliPID(species));
    hSecToPiRatio[species]->SetMarkerColor(getLineColorAliPID(species));
    hSecToPiRatio[species]->SetMarkerStyle(24);
    hSecToPiRatio[species]->SetLineWidth(2.);
    hSecToPiRatio[species]->SetStats(kFALSE);
    hSecToPiRatio[species]->GetYaxis()->SetTitle("Primary Fraction of Ratio");
    
    hSecSSToPiRatio[species] = new TH1D(*hSecToPiRatio[species]);
    hSecSSToPiRatio[species]->SetName(Form("hSecSSToPionRatio_%s", AliPID::ParticleShortName(species)));
    
    // Samples for different species are independent, so just divide correction factors
    hSecToPiRatio[species]->Divide(hSec[species], hSec[AliPID::kPion], 1., 1., ""); 
    hSecToPiRatio[species]->GetYaxis()->SetRangeUser(0., 1.1);
    
    hSecSSToPiRatio[species]->Divide(hSecSS[species], hSecSS[AliPID::kPion], 1., 1., ""); 
    hSecSSToPiRatio[species]->GetYaxis()->SetRangeUser(0., 1.1);
  }
  
  TCanvas* cSecToPiRatio = new TCanvas("cSecToPiRatio", "Primary fraction of to-#pi-ratio for different species", 100, 10, 1200, 800);
  cSecToPiRatio->SetGridx(1);
  cSecToPiRatio->SetGridy(1);
  cSecToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    hSecToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  cSecToPiRatio->BuildLegend()->SetFillColor(kWhite);
  
  ClearTitleFromHistoInCanvas(cSecToPiRatio);
  
  
  
  TCanvas* cSecSSToPiRatio = new TCanvas("cSecSSToPiRatio",
                                         "Primary fraction of to-#pi-ratio for different species (strangeness scaled)",
                                         100, 10, 1200, 800);
  cSecSSToPiRatio->SetGridx(1);
  cSecSSToPiRatio->SetGridy(1);
  cSecSSToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    hSecSSToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  cSecSSToPiRatio->BuildLegend()->SetFillColor(kWhite);
  
  ClearTitleFromHistoInCanvas(cSecSSToPiRatio);
  
  

  
  // If desired, apply the GEANT/FLUKA correction first
  // NOTE: This will change dataRebinned! If anything else should also be calculated, it must be done before.
  // Otherwise, the GEANT/FLUKA correction factor for the efficiency will also affect it!
  const Int_t genStepEff = kStepGenWithGenCuts;
  if (correctGeantFluka) {
    printf("Applying GEANT/FLUKA correction...\n");
    if (!geantFlukaCorrection(dataRebinned, genStepEff)) {
      printf("GEANT/FLUKA correction could not be applied!\n");
      return kFALSE;
    }
  }
  
  // Construct the efficiency grid from the data container 
  AliCFEffGrid* eff = new AliCFEffGrid("eff", "Efficiency x Acceptance x pT Resolution", *dataRebinned);
  
  // Either one can take kStepRecWithGenCutsMeasuredObs or, what I prefer, one can take
  // kStepRecWithRecCutsMeasuredObsPrimaries => The difference is only the eta cut, which is on the rec level
  // in the latter case, i.e. one corrects for eta resolution (although the effect should be very small)
  eff->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, genStepEff);
  
  // If the jet axis is restricted (i.e. jet input), scale with the corresponding number of jets.
  // Note: Since this is supposed to be a real scaling, set the error of the scale factor to zero
  // (second element of factor array)
  if (restrictJetPtAxis) {
    Double_t factor_Numerator[2] = { nJetsRec > 0 ? 1. / nJetsRec : 0., 0.  };
    Double_t factor_Denominator[2] = { nJetsGen > 0 ? 1. / nJetsGen : 0., 0.  };
    eff->GetNum()->Scale(factor_Numerator);
    eff->GetDen()->Scale(factor_Denominator);
  }
  
  //The efficiency along pt and vertex, and 2-D projection
  TCanvas* cEff = new TCanvas("cEff", "Efficiency x Acceptance x pT Resolution", 100, 10, 1200, 800);
  cEff->Divide(2, 1);
  cEff->GetPad(1)->SetLogx(kTRUE);
  cEff->cd(1);
  TH1* hEffPt = eff->Project(iPt); //the efficiency vs pt
  hEffPt->SetStats(kFALSE);
  hEffPt->GetXaxis()->SetRangeUser(0.15, 50.);
  hEffPt->GetXaxis()->SetMoreLogLabels(kTRUE);
  hEffPt->GetXaxis()->SetNoExponent(kTRUE);
  hEffPt->GetYaxis()->SetTitle("Efficiency x Acceptance x pT Resolution");
  hEffPt->Draw("E1");
  cEff->cd(2);
  TH1* hEffEta = eff->Project(iEta); //the efficiency vs eta
  hEffEta->SetStats(kFALSE);
  hEffEta->GetYaxis()->SetTitle("Efficiency x Acceptance x pT Resolution");
  hEffEta->GetYaxis()->SetTitleSize(0.05);
  hEffEta->Draw("E1");
  TH2D* hEffID2Pt = (TH2D*)eff->Project(iPt, iMCid);

  // Get the efficiencies vs. pT for each species
  TH1D* hEfficiency[AliPID::kSPECIES];
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hEfficiency[species] = hEffID2Pt->ProjectionX(Form("hEfficiency_%s", AliPID::ParticleShortName(species)), species + 1, species + 1, "e");
    hEfficiency[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hEfficiency[species]->SetLineColor(hYield[species]->GetLineColor());
    hEfficiency[species]->SetMarkerColor(hYield[species]->GetLineColor());
    hEfficiency[species]->SetLineWidth(2.);
    hEfficiency[species]->GetXaxis()->SetRangeUser(0.15, 50);
    hEfficiency[species]->GetYaxis()->SetRangeUser(0., 1.01);
    hEfficiency[species]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hEfficiency[species]->GetXaxis()->SetNoExponent(kTRUE);
    hEfficiency[species]->SetStats(kFALSE);
    hEfficiency[species]->GetYaxis()->SetTitle("Efficiency x Acceptance x pT Resolution");
    hEfficiency[species]->GetYaxis()->SetTitleSize(0.05);
  }
  
  TCanvas* cEff2 = new TCanvas("cEff2", "Efficiency x Acceptance x pT Resolution for different species", 0, 300, 900, 900);
  cEff2->SetGridx(1);
  cEff2->SetGridy(1);
  cEff2->SetLogx(1);
  
  hEfficiency[0]->Draw("E1");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    hEfficiency[i]->Draw("E1 same");
  }
  cEff2->BuildLegend()->SetFillColor(kWhite);
  
  ClearTitleFromHistoInCanvas(cEff2);
  
  
  
  
  // Efficiency correction for to-pi-ratios
  TH1D* hEfficiencyToPiRatio[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (species == AliPID::kPion)
      continue; // Do not consider pion-to-pion ratio
      
    hEfficiencyToPiRatio[species] = new TH1D(*hEfficiency[species]);
    hEfficiencyToPiRatio[species]->Reset();
    hEfficiencyToPiRatio[species]->SetName(Form("hEfficiencyToPionRatio_%s", AliPID::ParticleShortName(species)));
    hEfficiencyToPiRatio[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
    
    // Samples for different species are independent, so just divide correction factors
    hEfficiencyToPiRatio[species]->Divide(hEfficiency[species], hEfficiency[AliPID::kPion], 1., 1., ""); 
    hEfficiencyToPiRatio[species]->GetYaxis()->SetRangeUser(0., 2.0);
  }
  
  TCanvas* cEffToPiRatio = new TCanvas("cEffToPiRatio", "Efficiency x Acceptance x pT Resolution of to-#pi-ratio for different species",
                                       100, 10, 1200, 800);
  cEffToPiRatio->SetGridx(1);
  cEffToPiRatio->SetGridy(1);
  cEffToPiRatio->SetLogx(1);
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (i == AliPID::kPion)
      continue;
    
    hEfficiencyToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
  }
  
  cEffToPiRatio->BuildLegend()->SetFillColor(kWhite);
  
  ClearTitleFromHistoInCanvas(cEffToPiRatio);
  
  
  // Correct the yields with the efficiencies and primary fractions
  TH1D* hYieldCorrected[AliPID::kSPECIES] = { 0x0, };
  TH1D* hYieldCorrectedSysError[AliPID::kSPECIES] = { 0x0, };
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hYieldCorrected[species] = new TH1D(*hYield[species]);
    hYieldCorrected[species]->SetName(Form("%s_corrected", hYield[species]->GetName()));
    hYieldCorrected[species]->SetTitle(Form("%s", hYield[species]->GetTitle()));
    //hYieldCorrected[species]->SetTitle(Form("%s, secondary and efficiency x acceptance x pT resolution corrected",
    //                                        hYield[species]->GetTitle()));
    
    // Correction factor histos can have different binning than data histos (constant factor above some threshold)
    // -> Need special functions to multiply and divide such histos
    multiplyHistsDifferentBinning(hYieldCorrected[species], scaleStrangeness ? hSecSS[species] : hSec[species], 1., 1.);
    divideHistsDifferentBinning(hYieldCorrected[species], hEfficiency[species], 1., 1.);

    if (hYieldSysError[species]) {
      hYieldCorrectedSysError[species] = new TH1D(*hYieldSysError[species]);
      hYieldCorrectedSysError[species]->SetName(Form("%s_corrected", hYieldSysError[species]->GetName()));
      hYieldCorrectedSysError[species]->SetTitle(Form("%s", hYieldSysError[species]->GetTitle()));
      
      multiplyHistsDifferentBinning(hYieldCorrectedSysError[species], scaleStrangeness ? hSecSS[species] : hSec[species], 1., 1.,
                                    kTRUE);
      divideHistsDifferentBinning(hYieldCorrectedSysError[species], hEfficiency[species], 1., 1., kTRUE);
    }
  }
  
  // Calculate the total corrected yield. The original total yield had no error (just the number of detected tracks in a pT bin),
  // but due to the correction there is some error for the total yield. Also the error of the fractions introduces uncertainties
  // for the yields of individual species
  TH1D* hYieldCorrectedTotal = new TH1D(*hYieldCorrected[0]);
  hYieldCorrectedTotal->SetLineColor(kBlack);
  hYieldCorrectedTotal->SetMarkerColor(kBlack);
  hYieldCorrectedTotal->SetName("hYieldCorrectedTotal");
  hYieldCorrectedTotal->SetTitle("Total");
  //hYieldCorrectedTotal->SetTitle("Total yield, secondary and efficiency x acceptance x pT resolution corrected");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++)
    hYieldCorrectedTotal->Add(hYieldCorrected[i], 1.);
  
  // Calculate the corrected fractions
  TH1D* hFractionCorrected[AliPID::kSPECIES] = { 0x0, };
  TH1D* hFractionCorrectedSysError[AliPID::kSPECIES] = { 0x0, };
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hFractionCorrected[species] = new TH1D(*hYield[species]);
    TString oldName = hYield[species]->GetName();
    TString newName = oldName.ReplaceAll("Yield", "Fraction");
    TString oldTitle = hYield[species]->GetTitle();
    TString newTitle = oldTitle.ReplaceAll("Yield", "Fraction");
    hFractionCorrected[species]->SetName(Form("%s_corrected", newName.Data()));
    hFractionCorrected[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    hFractionCorrected[species]->GetYaxis()->SetTitle("Corrected Fraction");
    hFractionCorrected[species]->GetXaxis()->SetMoreLogLabels(kTRUE);
    hFractionCorrected[species]->GetXaxis()->SetNoExponent(kTRUE);

    // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
    // (numerator is a "selected" subset of the denominator). It doesn't matter that the error of the histos is not just "sqrt(content)"
    // because the error formula also works for weighted histograms (which means that the error can be more or less anything)
    hFractionCorrected[species]->Divide(hYieldCorrected[species], hYieldCorrectedTotal, 1., 1., "B"); 
    

    //  The systematic errors just need to be scaled in the same way as the fractions.
    // So, just take the ratio to the uncorrected fraction and scale the sys. error accordingly
    // or, in this case, just divide by the same total yield as for yield -> fractions
    if (hYieldCorrectedSysError[species]) {
      hFractionCorrectedSysError[species] = new TH1D(*hFractionCorrected[species]);
      hFractionCorrectedSysError[species]->SetName(Form("%s_sysError", hFractionCorrected[species]->GetName()));
      hFractionCorrectedSysError[species]->SetTitle(Form("%s (sys. error)", hFractionCorrected[species]->GetTitle()));
      
      for (Int_t binX = 1; binX <= hFractionCorrectedSysError[species]->GetNbinsX(); binX++) {
        const Double_t corrTotalYield = hYieldCorrectedTotal->GetBinContent(binX);
        const Double_t scaleFactor = corrTotalYield > 0 ? 1.0 / corrTotalYield : 1.;
        hFractionCorrectedSysError[species]->SetBinError(binX,   hYieldCorrectedSysError[species]->GetBinError(binX) * scaleFactor);
      }
    }
  }
  
  // If MC is available, calculate the generated fractions
  TH1D* hMCgenPrimYieldTotal = 0x0;
  TH1D* hMCgenPrimFraction[AliPID::kSPECIES];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    hMCgenPrimFraction[i] = 0x0;
  
  if (numMCgenPrimYieldHistsFound > 0) {
    hMCgenPrimYieldTotal = new TH1D(*hMCgenPrimYield[0]);
    hMCgenPrimYieldTotal->SetLineColor(kBlack);
    hMCgenPrimYieldTotal->SetMarkerColor(kBlack);
    hMCgenPrimYieldTotal->SetName("hMCgenPrimYieldTotal");
    hMCgenPrimYieldTotal->SetTitle("Total (MC truth)");
    //hMCgenPrimYieldTotal->SetTitle("Total generated primary yield (MC truth)");
    
    for (Int_t i = 1; i < AliPID::kSPECIES; i++)
      hMCgenPrimYieldTotal->Add(hMCgenPrimYield[i], 1.);
    
    // Calculate the MC fractions
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hMCgenPrimYield[species]->SetTitle(Form("%s (MC truth)", AliPID::ParticleLatexName(species)));
      
      hMCgenPrimFraction[species] = new TH1D(*hMCgenPrimYield[species]);
      TString oldName = hMCgenPrimYield[species]->GetName();
      TString newName = oldName.ReplaceAll("Yield", "Fraction");
      TString oldTitle = hMCgenPrimYield[species]->GetTitle();
      TString newTitle = oldTitle.ReplaceAll("yield", "fraction");
      hMCgenPrimFraction[species]->SetName(newName.Data());
      hMCgenPrimFraction[species]->SetTitle(newTitle.Data());

      // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
      // (numerator is a "selected" subset of the denominator).
      hMCgenPrimFraction[species]->Divide(hMCgenPrimFraction[species], hMCgenPrimYieldTotal, 1., 1., "B"); 
    }
  }
  
  TCanvas* cCorrData = new TCanvas("cCorrData", "Corrected data", 0, 300, 900, 900);
  cCorrData->Divide(2, 1);//, 0., 0.01);
  cCorrData->GetPad(1)->SetLogx(1);
  cCorrData->GetPad(1)->SetLogy(1);
  cCorrData->GetPad(2)->SetLogx(1);
  cCorrData->GetPad(2)->SetLogy(1);
  
  cCorrData->GetPad(1)->SetRightMargin(0.001);
  cCorrData->GetPad(2)->SetRightMargin(0.001);
  
  cCorrData->GetPad(1)->SetLeftMargin(0.2);
  cCorrData->GetPad(2)->SetLeftMargin(0.2);
  
  cCorrData->cd(1); // uncorrected
  hYield[AliPID::kPion]->GetYaxis()->SetTitleOffset(1.4);
  hYield[AliPID::kPion]->Draw("E1");
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldSysError[i])
      hYieldSysError[i]->Draw("E2 same");
    
    if (i == AliPID::kPion) continue;
    hYield[i]->Draw("E1 same");
  }
  
  ClearTitleFromHistoInCanvas(cCorrData, 1);
  
  cCorrData->cd(2); // corrected
  hYieldCorrectedTotal->GetYaxis()->SetTitleOffset(1.4);
  hYieldCorrectedTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
  hYieldCorrectedTotal->GetXaxis()->SetNoExponent(kTRUE);
  hYieldCorrectedTotal->GetYaxis()->SetRangeUser(hYieldCorrected[AliPID::kMuon]->GetBinContent(hYieldCorrected[AliPID::kMuon]->FindLastBinAbove(0.)) * 0.1,
                                                 hYieldCorrectedTotal->GetBinContent(hYieldCorrectedTotal->GetMaximumBin()) * 10.);
  
  if (hMCgenPrimYieldTotal) {
    hMCgenPrimYieldTotal->GetYaxis()->SetTitleOffset(1.4);
    hMCgenPrimYieldTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
    hMCgenPrimYieldTotal->GetXaxis()->SetNoExponent(kTRUE);
    hMCgenPrimYieldTotal->GetXaxis()->SetTitle(hYieldCorrectedTotal->GetXaxis()->GetTitle());
    hMCgenPrimYieldTotal->GetYaxis()->SetRangeUser(hYieldCorrected[AliPID::kMuon]->GetBinContent(hYieldCorrected[AliPID::kMuon]->FindLastBinAbove(0.)) * 0.1,
                                                   hYieldCorrectedTotal->GetBinContent(hYieldCorrectedTotal->GetMaximumBin()) * 10.);
  }
  
  if (hMCgenPrimYieldTotal) {
    hMCgenPrimYieldTotal->Draw("E1");
    hYieldCorrectedTotal->Draw("E1 same");
  }
  else
    hYieldCorrectedTotal->Draw("E1");
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenPrimYield[i])
      hMCgenPrimYield[i]->Draw("E1 same");
    hYieldCorrected[i]->Draw("E1 same");
  }
  
  TLegend* legTemp = cCorrData->cd(2)->BuildLegend(0.25, 0.16, 0.65, 0.51);
  legTemp->SetNColumns(2);
  legTemp->SetFillColor(kWhite);
  
  // Do not include in legend
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrectedSysError[i])
      hYieldCorrectedSysError[i]->Draw("E2 same");
  }
  
  ClearTitleFromHistoInCanvas(cCorrData, 2);
  
  TCanvas* cCorrYieldsRatio = 0x0;
  
  TH1D* hYieldCorrectedRatioToMC[AliPID::kSPECIES];
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) 
    hYieldCorrectedRatioToMC[i] = 0x0;

  TH1D* hYieldCorrectedTotalRatioToMC = 0x0;
  
  if (numMCgenPrimYieldHistsFound > 0) {
    // Compare with MC truth
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hYieldCorrectedRatioToMC[species] = new TH1D(*hYieldCorrected[species]);
      hYieldCorrectedRatioToMC[species]->SetName(Form("hYieldCorrectedRatioToMC_%s", AliPID::ParticleShortName(species)));
      hYieldCorrectedRatioToMC[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
      hYieldCorrectedRatioToMC[species]->GetYaxis()->SetTitle("Corrected Yield: Fit / MC Truth");
      hYieldCorrectedRatioToMC[species]->Divide(hMCgenPrimYield[species]);
    }
    
    hYieldCorrectedTotalRatioToMC = new TH1D(*hYieldCorrectedTotal);
    hYieldCorrectedTotalRatioToMC->SetName("hYieldCorrectedTotalRatioToMC");
    hYieldCorrectedTotalRatioToMC->SetTitle("Total");
    hYieldCorrectedTotalRatioToMC->GetYaxis()->SetTitle("Corrected Yield: Fit / MC Truth");
    hYieldCorrectedTotalRatioToMC->Divide(hMCgenPrimYieldTotal);
    
    cCorrYieldsRatio = new TCanvas("cCorrYieldsRatio", "Corrected Yields Comparison to MC", 0, 300, 900, 900);
    cCorrYieldsRatio->SetGridx(1);
    cCorrYieldsRatio->SetGridy(1);
    cCorrYieldsRatio->SetLogx(1);
    
    hYieldCorrectedTotalRatioToMC->GetYaxis()->SetRangeUser(0.6, 1.6);
    hYieldCorrectedTotalRatioToMC->GetYaxis()->SetTitleOffset(0.85);
    hYieldCorrectedTotalRatioToMC->Draw("E1");
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++)
      hYieldCorrectedRatioToMC[species]->Draw("E1 same");
    
    cCorrYieldsRatio->BuildLegend()->SetFillColor(kWhite);
    
    ClearTitleFromHistoInCanvas(cCorrYieldsRatio);
  }
  
  
  if (cFractions)
    cFractions->Draw();
  
  TCanvas* cCorrFractions = new TCanvas("cCorrFractions", "Corrected particleFractions", 0, 300, 900, 900);
  cCorrFractions->SetLogx(1);
  hFractionCorrected[0]->GetYaxis()->SetRangeUser(0., 1.);
  hFractionCorrected[0]->Draw("E1");
  if (hMCgenPrimFraction[0])
    hMCgenPrimFraction[0]->Draw("E1 same");
  
  for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
    hFractionCorrected[i]->Draw("E1 same");
    if (hMCgenPrimFraction[i])
      hMCgenPrimFraction[i]->Draw("E1 same");
  }
  
  cCorrFractions->BuildLegend()->SetFillColor(kWhite);
  
  // Do not include in legend!!
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hFractionCorrectedSysError[i])
      hFractionCorrectedSysError[i]->Draw("E2 same");
  }
  
  
  ClearTitleFromHistoInCanvas(cCorrFractions);
  
  // Save results to file
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  TString saveFileName = pathNameData;
  saveFileName.ReplaceAll(Form("%s/", pathData.Data()), "");
  saveFileName.Prepend("output_EfficiencyCorrection_");
  saveFileName.ReplaceAll(".root", Form("%s.root", chargeString.Data()));
  
  TString saveFilePathName = Form("%s/%s", pathData.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  saveFile->cd();
  
  if (cSec)
    cSec->Write();
  
  if (cSecSS)
    cSecSS->Write();
  
  if (cSec2)
    cSec2->Write();
  
  if (cSecSS2)
    cSecSS2->Write();
  
  if (cSecToPiRatio)
    cSecToPiRatio->Write();
  
  if (cSecSSToPiRatio)
    cSecSSToPiRatio->Write();
  
  if (cEff)
    cEff->Write();
  
  if (cEff2)
    cEff2->Write();
  
  if (cEffToPiRatio)
    cEffToPiRatio->Write();
  
  if (cCorrData)
    cCorrData->Write();
  
  if (cCorrYieldsRatio)
    cCorrYieldsRatio->Write();
  
  if (cFractions)
    cFractions->Write();
  
  if (cCorrFractions)
    cCorrFractions->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSec[i])
      hSec[i]->Write();
    
    if (hSecSS[i])
      hSecSS[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecToPiRatio[i])
      hSecToPiRatio[i]->Write();
    
    if (hSecSSToPiRatio[i])
      hSecSSToPiRatio[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiency[i])
      hEfficiency[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyToPiRatio[i])
      hEfficiencyToPiRatio[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYield[i])
      hYield[i]->Write();
    
    if (hYieldSysError[i])
      hYieldSysError[i]->Write();
  }
  
  if (hYieldCorrectedTotal)
    hYieldCorrectedTotal->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrected[i])
      hYieldCorrected[i]->Write();
    
    if (hYieldCorrectedSysError[i])
      hYieldCorrectedSysError[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenPrimYield[i])
      hMCgenPrimYield[i]->Write();
  }
  
  if (hMCgenPrimYieldTotal)
      hMCgenPrimYieldTotal->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hYieldCorrectedRatioToMC[i])
      hYieldCorrectedRatioToMC[i]->Write();
  }
  
  if (hYieldCorrectedTotalRatioToMC)
      hYieldCorrectedTotalRatioToMC->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hFractionCorrected[i])
      hFractionCorrected[i]->Write();
    
    if (hFractionCorrectedSysError[i])
      hFractionCorrectedSysError[i]->Write();
  }
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hMCgenPrimFraction[i])
      hMCgenPrimFraction[i]->Write();
  }
  
  TNamed* settings = new TNamed(
      Form("Settings: Efficiency file \"%s\", Data file \"%s\", lowerCentrality %.3f, upperCentrality %.3f, lowerJetPt %.1f, upperJetPt %.1f, constantCorrectionAbovePtThreshold %.3f\n",
           pathNameEfficiency.Data(), pathNameData.Data(), lowerCentrality, upperCentrality, lowerJetPt,
           upperJetPt, constantCorrectionAbovePtThreshold), "");
  settings->Write();
  
  saveFile->Close();
  
  return 0;
}


//___________________________________________________________________
// Efficiency for jet spectra (pT, z and xi)
// E.g. a 'calcEfficiency.C+("finalCuts/MC_pp/7TeV/LHC11a1/corrected/allJetPtMergedWeighted/bhess_PID_Jets_PureGauss_efficiency.root", "finalCuts/pp/7TeV/LHC10e.pass2/corrected/finalisedSplines/finalMapsAndTail/Jets/noCutOn_ncl_or_liav/output_extractedFFs_bhess_PID_Jets_centrality_all_highN_rougherXi.root", kTRUE, kFALSE, 0, -2, -2, -2, -2, 40, 80, -100, -100, 2, 1, 2)' -b -q
Int_t calcEfficiency(TString pathNameEfficiency, TString pathNameData, Bool_t correctGeantFluka, Bool_t scaleStrangeness,
                     Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                     Double_t lowerCentralityData /*= -2*/, Double_t upperCentralityData /*= -2*/,
                     Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                     Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/,
                     Double_t constantCorrectionAbovePtThreshold,
                     Double_t constantCorrectionAboveXiThreshold,
                     Int_t rebinEfficiencyPt,
                     Int_t rebinEfficiencyZ,
                     Int_t rebinEfficiencyXi)
{
  TString titles[kNtypes];
  titles[kTrackPt] = "trackPt";
  titles[kZ] = "z";
  titles[kXi] = "xi";
  
  TString pathData = pathNameData;
  pathData.Replace(pathData.Last('/'), pathData.Length(), "");
  
  TFile* fileEff = TFile::Open(pathNameEfficiency.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", pathNameEfficiency.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
  }
  
  const Int_t iMCid   = data->GetVar("MC ID");
  const Int_t iPt     = data->GetVar("P_{T} (GeV/c)");
  //const Int_t iEta    = data->GetVar("#eta");
  const Int_t iCharge = data->GetVar("Charge (e_{0})");
  const Int_t iMult   = data->GetVar("Centrality Percentile");
  const Int_t iJetPt  = data->GetVar("P_{T}^{jet} (GeV/c)");
  const Int_t iZ      = data->GetVar("z = P_{T}^{track} / P_{T}^{jet}");
  const Int_t iXi     = data->GetVar("#xi = ln(P_{T}^{jet} / P_{T}^{track})");
  
  TH2D* hNjetsGen = 0x0;
  TH2D* hNjetsRec = 0x0;
  
  TH2D* hNjetsGenData = 0x0;
  TH2D* hNjetsRecData = 0x0;
  
  TH1D* hYieldPt[AliPID::kSPECIES] = {0x0, };
  TH1D* hYieldPtSysError[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCgenPrimYieldPt[AliPID::kSPECIES] = {0x0, };
  
  TH1D* hYieldZ[AliPID::kSPECIES] = {0x0, };
  TH1D* hYieldZSysError[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCgenPrimYieldZ[AliPID::kSPECIES] = {0x0, };
  
  TH1D* hYieldXi[AliPID::kSPECIES] = {0x0, };
  TH1D* hYieldXiSysError[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCgenPrimYieldXi[AliPID::kSPECIES] = {0x0, };
  
  
  TH2D* hTemp = 0x0;
  
  Int_t numMCgenPrimYieldHistsFound = 0;
  
  Int_t lowerJetPtBinYields = -1;
  Int_t upperJetPtBinYields = -1;
  
  Double_t actualLowerJetPtYields = -10;
  Double_t actualUpperJetPtYields = -10;
  
  TFile* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  
  hNjetsGenData = (TH2D*)fileData->Get("fh2FFJetPtGen");
  
  hNjetsRecData = (TH2D*)fileData->Get("fh2FFJetPtRec");
  
  if (!hNjetsRecData) {
    printf("Failed to load numJet histo for data!\n");
    return -1;
  }
  
  const Bool_t restrictCentralityData = ((lowerCentralityData >= -1) && (upperCentralityData >= -1));
  // Integral(lowerCentBinLimitData, uppCentBinLimitData) will not be restricted if these values are kept
  const Int_t lowerCentralityBinLimitData = restrictCentralityData ? hNjetsRecData->GetXaxis()->FindBin(lowerCentralityData + 0.001) 
                                                                   : -1;
  const Int_t upperCentralityBinLimitData = restrictCentralityData ? hNjetsRecData->GetXaxis()->FindBin(upperCentralityData - 0.001) 
                                                                   : -2;
  
  // Integral(lowerJetBinLimitData, uppJetBinLimitData) will not be restricted if these values are kept
  const Int_t lowerJetPtBinLimitData = ((lowerJetPt >= 0) && (upperJetPt >= 0)) 
                                       ? hNjetsRecData->GetYaxis()->FindBin(lowerJetPt + 0.001) : -1;
  const Int_t upperJetPtBinLimitData  = ((lowerJetPt >= 0) && (upperJetPt >= 0)) 
                                       ? hNjetsRecData->GetYaxis()->FindBin(upperJetPt - 0.001) : -2;
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    TString histName = Form("hIDFFz_%s", AliPID::ParticleShortName(species));
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp && (lowerJetPtBinYields < 0 || upperJetPtBinYields < 0)) {
      // If the jetPt range is to be restricted
      if (lowerJetPt >= 0 && upperJetPt >= 0) {
        lowerJetPtBinYields = hTemp->GetYaxis()->FindBin(lowerJetPt + 0.001);
        upperJetPtBinYields = hTemp->GetYaxis()->FindBin(upperJetPt - 0.001);
        
        actualLowerJetPtYields = hTemp->GetYaxis()->GetBinLowEdge(lowerJetPtBinYields);
        actualUpperJetPtYields = hTemp->GetYaxis()->GetBinUpEdge(upperJetPtBinYields);
      }
      else {
        lowerJetPtBinYields = 0;
        upperJetPtBinYields = -1;
      }
    }
    
    if ((lowerJetPtBinYields < 0 || upperJetPtBinYields < 0) && (lowerJetPt >= 0 && upperJetPt >= 0)) {
      printf("Failed to determine jet pt bin limits!\n");
      return -1;
    }
    
    // NOTE: The yields have already been normalised to nJets in extractFFs. However, the binning in jetPt might be different now.
    // Thus, in order to get the correct result, the old normalisation (just the same binnnig as in the nJet histo) has to be 
    // undone and then it has to be normalised to the number of jets in the new binning again (see below).
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsRecData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hYieldZ[species] = hTemp->ProjectionX(histName.ReplaceAll("hIDFF", "hYieldIDFF").Data(), lowerJetPtBinYields, upperJetPtBinYields, 
                                            "e");
    }
      
    if (!hYieldZ[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    
    
    histName = Form("hIDFFxi_%s", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsRecData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hYieldXi[species] = hTemp->ProjectionX(histName.ReplaceAll("hIDFF", "hYieldIDFF").Data(),
                                             lowerJetPtBinYields, upperJetPtBinYields, "e");
    }
    
    if (!hYieldXi[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    histName = Form("hIDFFtrackPt_%s", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsRecData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hYieldPt[species] = hTemp->ProjectionX(histName.ReplaceAll("hIDFF", "hYieldIDFF").Data(), 
                                             lowerJetPtBinYields, upperJetPtBinYields, "e");
    }
    
    if (!hYieldPt[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    
    
    histName = Form("hIDFFz_%s_sysError", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsRecData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hYieldZSysError[species] = hTemp->ProjectionX(histName.ReplaceAll("hIDFF", "hYieldIDFF").Data(), lowerJetPtBinYields, 
                                                    upperJetPtBinYields, "e");
    }
    
    if (!hYieldZSysError[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    histName = Form("hIDFFxi_%s_sysError", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsRecData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hYieldXiSysError[species] = hTemp->ProjectionX(histName.ReplaceAll("hIDFF", "hYieldIDFF").Data(), lowerJetPtBinYields, 
                                                     upperJetPtBinYields, "e");
    }
    
    if (!hYieldXiSysError[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    histName = Form("hIDFFtrackPt_%s_sysError", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsRecData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hYieldPtSysError[species] = hTemp->ProjectionX(histName.ReplaceAll("hIDFF", "hYieldIDFF").Data(), lowerJetPtBinYields, 
                                                     upperJetPtBinYields, "e");
    }
    
    if (!hYieldPtSysError[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    // In case of MC also retrieve the MC truth generated yields
    histName = Form("hMCgenYieldsPrimZ_%s", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsGenData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hMCgenPrimYieldZ[species] = hTemp->ProjectionX(histName.ReplaceAll("hMCgen", "hYieldMCgen").Data(), lowerJetPtBinYields,
                                                     upperJetPtBinYields, "e");
    }
    
    
    histName = Form("hMCgenYieldsPrimXi_%s", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsGenData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hMCgenPrimYieldXi[species] = hTemp->ProjectionX(histName.ReplaceAll("hMCgen", "hYieldMCgen").Data(), lowerJetPtBinYields,
                                                      upperJetPtBinYields, "e");
    }
    
    
    histName = Form("hMCgenYieldsPrimPt_%s", AliPID::ParticleShortName(species));
    
    hTemp = (TH2D*)fileData->Get(histName.Data());
    
    if (hTemp) {
      undoJetNormalisationHist2D(hTemp, hNjetsGenData, lowerCentralityBinLimitData, upperCentralityBinLimitData);
      hMCgenPrimYieldPt[species] = hTemp->ProjectionX(histName.ReplaceAll("hMCgen", "hYieldMCgen").Data(), lowerJetPtBinYields,
                                                      upperJetPtBinYields, "e");
    }
    
    if (hMCgenPrimYieldPt[species] && hMCgenPrimYieldZ[species] && hMCgenPrimYieldXi[species])
      numMCgenPrimYieldHistsFound++;
  }
  
  if (numMCgenPrimYieldHistsFound > 0 && numMCgenPrimYieldHistsFound != AliPID::kSPECIES) {
    printf("Error: Unable to retrieve all MC generated prim yield histos! Got %d.\n", numMCgenPrimYieldHistsFound);
    return -1;
  }
  
  
  /*OLD: Load from AnalysisResults.root
  
  File* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TString dirDataInFile = "";
  TDirectory* dirData = (TDirectory*)fileData->Get(fileData->GetListOfKeys()->At(0)->GetName());
  if (!dirData) {
    printf("No directory found inside data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TList* list = (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName());
  if (!list) {
    printf("No list found in directory \"%s\" inside data file \"%s\"\n",
           dirData->GetListOfKeys()->At(0)->GetName(),pathNameData.Data());
    return -1;
  }
  
  hNjetsGen = (TH1D*)list->FindObject("fh1FFJetPtGen");
  hNjetsRec = (TH1D*)list->FindObject("fh1FFJetPtRecCuts");
  
  if (!hNjetsRec) {
    printf("Failed to load number of jets (rec) histo!\n");
    return -1;
  }
  
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    TString histName = Form("fh2FFZRecCuts_%s", AliPID::ParticleName(species));
    hTemp = (TH2D*)list->FindObject(histName.Data());
    
    if (hTemp && (lowerJetPtBinYields < 0 || upperJetPtBinYields < 0)) {
      // If the jetPt range is to be restricted
      if (lowerJetPt >= 0 && upperJetPt >= 0) {
        lowerJetPtBinYields = hTemp->GetXaxis()->FindBin(lowerJetPt + 0.001);
        upperJetPtBinYields = hTemp->GetXaxis()->FindBin(upperJetPt - 0.001);
        
        actualLowerJetPtYields = hTemp->GetXaxis()->GetBinLowEdge(lowerJetPtBinYields);
        actualUpperJetPtYields = hTemp->GetXaxis()->GetBinUpEdge(upperJetPtBinYields);
      }
      else {
        lowerJetPtBinYields = 0;
        upperJetPtBinYields = -1;
      }
    }
    
    if ((lowerJetPtBinYields < 0 || upperJetPtBinYields < 0) && (lowerJetPt >= 0 && upperJetPt >= 0)) {
      printf("Failed to determine jet pt bin limits!\n");
      return -1;
    }
    
    if (hTemp)
      hYieldZ[species] = hTemp->ProjectionY(histName.ReplaceAll("fh2", "fh").Data(), lowerJetPtBinYields, upperJetPtBinYields, "e");
      
    if (!hYieldZ[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    
    
    histName = Form("fh2FFXiRecCuts_%s", AliPID::ParticleName(species));
    
    hTemp = (TH2D*)list->FindObject(histName.Data());
    
    if (hTemp)
      hYieldXi[species] = hTemp->ProjectionY(histName.ReplaceAll("fh2", "fh").Data(), lowerJetPtBinYields, upperJetPtBinYields, "e");
    
    if (!hYieldXi[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    histName = Form("fh2FFTrackPtRecCuts_%s", AliPID::ParticleName(species));
    
    hTemp = (TH2D*)list->FindObject(histName.Data());
    
    if (hTemp)
      hYieldPt[species] = hTemp->ProjectionY(histName.ReplaceAll("fh2", "fh").Data(), lowerJetPtBinYields, upperJetPtBinYields, "e");
    
    if (!hYieldPt[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    
    // In case of MC also retrieve the MC truth generated yields
    histName = Form("fh2FFZGen_%s", AliPID::ParticleName(species));
    
    hTemp = (TH2D*)list->FindObject(histName.Data());
    if (hTemp)
      hMCgenPrimYieldZ[species] = hTemp->ProjectionY(histName.ReplaceAll("fh2", "fh").Data(), lowerJetPtBinYields, upperJetPtBinYields, 
                                                     "e");
    
    
    histName = Form("fh2FFXiGen_%s", AliPID::ParticleName(species));
    
    hTemp = (TH2D*)list->FindObject(histName.Data());
    if (hTemp)
      hMCgenPrimYieldXi[species] = hTemp->ProjectionY(histName.ReplaceAll("fh2", "fh").Data(), lowerJetPtBinYields, upperJetPtBinYields, 
                                                      "e");
    
    
    histName = Form("fh2FFTrackPtGen_%s", AliPID::ParticleName(species));
    
    hTemp = (TH2D*)list->FindObject(histName.Data());
    if (hTemp)
      hMCgenPrimYieldPt[species] = hTemp->ProjectionY(histName.ReplaceAll("fh2", "fh").Data(), lowerJetPtBinYields, upperJetPtBinYields, 
                                                      "e");
    
    if (hMCgenPrimYieldPt[species] && hMCgenPrimYieldZ[species] && hMCgenPrimYieldXi[species])
      numMCgenPrimYieldHistsFound++;
  }
  
  if (numMCgenPrimYieldHistsFound > 0 && numMCgenPrimYieldHistsFound != AliPID::kSPECIES) {
    printf("Error: Unable to retrieve all MC generated prim yield histos! Got %d.\n", numMCgenPrimYieldHistsFound);
    return -1;
  }
  
  
  const Double_t nJetsGen = hNjetsGen ? hNjetsGen->Integral(lowerJetPtBinYields, upperJetPtBinYields) : -1;
  const Double_t nJetsRec = hNjetsRec->Integral(lowerJetPtBinYields, upperJetPtBinYields);
  */
  
  // Take binning from pion yield (binning for all species the same) and create a new AliCFContainer with this new binning
  TH1D hDummyPt(*hYieldPt[AliPID::kPion]);
  hDummyPt.SetName("hDummyPt");
  TAxis* axisPt = 0x0;
  
  if (rebinEfficiencyPt > 1) {
    Int_t nBinsNew = 0;
    const Double_t* newBins = GetBins(rebinEfficiencyPt, nBinsNew);
    
    hDummyPt.SetBins(nBinsNew, newBins);
    
    axisPt = hDummyPt.GetXaxis();
  
    //axisPt = hDummyPt.Rebin(rebinEfficiencyPt, "", 0)->GetXaxis();
  }
  else
    axisPt = hDummyPt.GetXaxis();
  
  const TArrayD* binsPtCurrent = axisPt->GetXbins();
  TArrayD* binsPtNew = new TArrayD(*binsPtCurrent);
  
  TH1D hDummyZ(*hYieldZ[AliPID::kPion]);
  hDummyZ.SetName("hDummyZ");
  TAxis* axisZ = 0x0;
  
  if (rebinEfficiencyZ > 1)
    axisZ = hDummyZ.Rebin(rebinEfficiencyZ, "", 0)->GetXaxis();
  else
    axisZ = hDummyZ.GetXaxis();

  const TArrayD* binsZCurrent = axisZ->GetXbins();
  TArrayD* binsZNew = new TArrayD(*binsZCurrent);
  
  TH1D hDummyXi(*hYieldXi[AliPID::kPion]);
  hDummyXi.SetName("hDummyXi");
  TAxis* axisXi = 0x0;
  
  if (rebinEfficiencyXi > 1)
    axisXi = hDummyXi.Rebin(rebinEfficiencyXi, "", 0)->GetXaxis();
  else
    axisXi = hDummyXi.GetXaxis();

  const TArrayD* binsXiCurrent = axisXi->GetXbins();
  TArrayD* binsXiNew = new TArrayD(*binsXiCurrent);
  
  
  const Int_t nEffDims = data->GetNVar();
  Int_t nEffBins[nEffDims];
  
  for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
    if (iDim == iPt)
      nEffBins[iDim] = axisPt->GetNbins();
    else if (iDim == iZ)
      nEffBins[iDim] = axisZ->GetNbins();
    else if (iDim == iXi)
      nEffBins[iDim] = axisXi->GetNbins();
    else 
      nEffBins[iDim] = data->GetNBins(iDim);
  }
  
 
  // Just make one large pT bin above some threshold, if desired
  if (binsPtNew->fN != 0 && constantCorrectionAbovePtThreshold > 0) {
    for (Int_t iBin = 0; iBin <= nEffBins[iPt] + 1; iBin++) {
      // Find the first bin edged really larger than the threshold.
      // If the bin edge before equals the threshold, just set the
      // current bin edge to the right end of the spectrum -> Done.
      // If the bin edge before is different, set the bin edge to the
      // threshold
      if (binsPtNew->fArray[iBin] > constantCorrectionAbovePtThreshold) {
        if (binsPtNew->fArray[iBin - 1] == constantCorrectionAbovePtThreshold) {
          binsPtNew->fArray[iBin] = binsPtNew->fArray[nEffBins[iPt]];
          nEffBins[iPt] = iBin;
          break;
        }
        else {
          binsPtNew->fArray[iBin] = constantCorrectionAbovePtThreshold;
        }
      }
    }
  }
  
  // Just make one large Xi bin above some threshold, if desired
  if (binsXiNew->fN != 0 && constantCorrectionAboveXiThreshold > 0) {
    for (Int_t iBin = 0; iBin <= nEffBins[iXi] + 1; iBin++) {
      // Find the first bin edged really larger than the threshold.
      // If the bin edge before equals the threshold, just set the
      // current bin edge to the right end of the spectrum -> Done.
      // If the bin edge before is different, set the bin edge to the
      // threshold
      if (binsXiNew->fArray[iBin] > constantCorrectionAboveXiThreshold) {
        if (binsXiNew->fArray[iBin - 1] == constantCorrectionAboveXiThreshold) {
          binsXiNew->fArray[iBin] = binsXiNew->fArray[nEffBins[iXi]];
          nEffBins[iXi] = iBin;
          break;
        }
        else {
          binsXiNew->fArray[iBin] = constantCorrectionAboveXiThreshold;
        }
      }
    }
  }
  
  AliCFContainer *dataRebinned = new AliCFContainer(Form("%s_rebinned", data->GetName()), Form("%s (rebinned)", data->GetTitle()),
                                                    data->GetNStep(), nEffDims, nEffBins);
  
  for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
    dataRebinned->SetVarTitle(iDim, data->GetVarTitle(iDim));
    
    if (iDim == iPt) {
      if (binsPtNew->fN == 0)
        dataRebinned->SetBinLimits(iDim, axisPt->GetXmin(), axisPt->GetXmax());
      else
        dataRebinned->SetBinLimits(iDim, binsPtNew->fArray);
    }
    else if (iDim == iZ) {
      if (binsZNew->fN == 0)
        dataRebinned->SetBinLimits(iDim, axisZ->GetXmin(), axisZ->GetXmax());
      else
        dataRebinned->SetBinLimits(iDim, binsZNew->fArray);
    }
    else if (iDim == iXi) {
      if (binsXiNew->fN == 0)
        dataRebinned->SetBinLimits(iDim, axisXi->GetXmin(), axisXi->GetXmax());
      else
        dataRebinned->SetBinLimits(iDim, binsXiNew->fArray);
    }
    else {
      dataRebinned->SetBinLimits(iDim, data->GetBinLimits(iDim));
    }
  }
  
  for (Int_t iStep = 0; iStep < data->GetNStep(); iStep++)
    dataRebinned->SetStepTitle(iStep, data->GetStepTitle(iStep));
  
  Int_t coord[nEffDims];
  Double_t binCenterCoord[nEffDims];
  
  // Fill content from old grid into the new grid with proper binning
  for (Int_t iStep = 0; iStep < data->GetNStep(); iStep++) {
    Long64_t nBinsGrid = data->GetGrid(iStep)->GetGrid()->GetNbins();
    
    for (Long64_t iBin = 0; iBin <= nBinsGrid + 1; iBin++) {
      Double_t binContent = data->GetGrid(iStep)->GetGrid()->GetBinContent(iBin, coord);
      Double_t binError2  = data->GetGrid(iStep)->GetGrid()->GetBinError2(iBin);
      
      for (Int_t iDim = 0; iDim < nEffDims; iDim++) {
        binCenterCoord[iDim] = data->GetBinCenter(iDim, coord[iDim]);
      }

      Long64_t iBinRebinned = dataRebinned->GetGrid(iStep)->GetGrid()->GetBin(binCenterCoord);
      dataRebinned->GetGrid(iStep)->GetGrid()->AddBinContent(iBinRebinned, binContent);
      dataRebinned->GetGrid(iStep)->GetGrid()->AddBinError2(iBinRebinned, binError2);
    }
  }
  
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = dataRebinned->GetAxis(iMult, 0)->FindBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = dataRebinned->GetAxis(iMult, 0)->FindBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= dataRebinned->GetAxis(iMult, 0)->GetNbins()) {
      actualLowerCentrality = dataRebinned->GetAxis(iMult, 0)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = dataRebinned->GetAxis(iMult, 0)->GetBinUpEdge(upperCentralityBinLimit);

      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis) {
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
    dataRebinned->SetRangeUser(iMult, lowerCentralityBinLimit, upperCentralityBinLimit, kTRUE);
    data->SetRangeUser(iMult, lowerCentralityBinLimit, upperCentralityBinLimit, kTRUE);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -2; // Integral(lowerJetBinLimit, uppJetBinLimit) will not be restricted if these values are kept
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 &&
        upperJetPtBinLimit <= dataRebinned->GetAxis(iJetPt, 0)->GetNbins()) {
      actualLowerJetPt = dataRebinned->GetAxis(iJetPt, 0)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = dataRebinned->GetAxis(iJetPt, 0)->GetBinUpEdge(upperJetPtBinLimit);

      restrictJetPtAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested jet pT range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "jet pT: ";
  if (restrictJetPtAxis) {
    std::cout << actualLowerJetPt << " - " << actualUpperJetPt << std::endl;
    dataRebinned->SetRangeUser(iJetPt, lowerJetPtBinLimit, upperJetPtBinLimit, kTRUE);
    data->SetRangeUser(iJetPt, lowerJetPtBinLimit, upperJetPtBinLimit, kTRUE);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  if (restrictJetPtAxis && 
      (TMath::Abs(actualLowerJetPt - actualLowerJetPtYields) > 1e-3 || TMath::Abs(actualUpperJetPt - actualUpperJetPtYields) > 1e-3)) {
    std::cout << "ERROR: Disagreement between jetPt bounds of data and efficiency!" << std::endl;
    return -1;
  }
  
  // If desired, restrict charge axis
  std::cout << "Charge selection (efficiency): ";
  if (chargeMode == kAllCharged)
    std::cout << "All charged particles" << std::endl;
  else if (chargeMode == kNegCharge)
    std::cout << "Negative particles only" << std::endl;
  else if (chargeMode == kPosCharge)
    std::cout << "Positive particles only" << std::endl;
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  const Bool_t restrictCharge = (chargeMode != kAllCharged);
  
  Int_t lowerChargeBinLimit = -1;
  Int_t upperChargeBinLimit = -2;
  Double_t actualLowerCharge = -999;
  Double_t actualUpperCharge = -999;
  
  if (restrictCharge) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(-1. + 0.001);
      upperChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(0. + 0.001);
      upperChargeBinLimit = dataRebinned->GetAxis(iCharge, 0)->FindBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimit <= upperChargeBinLimit && lowerChargeBinLimit >= 1
        && upperChargeBinLimit <= dataRebinned->GetAxis(iCharge, 0)->GetNbins()) {
      actualLowerCharge = dataRebinned->GetAxis(iCharge, 0)->GetBinLowEdge(lowerChargeBinLimit);
      actualUpperCharge = dataRebinned->GetAxis(iCharge, 0)->GetBinUpEdge(upperChargeBinLimit);
      
      std::cout << "Charge range (efficiency): " << actualLowerCharge << " - " << actualUpperCharge << std::endl;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range (efficiency) out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    dataRebinned->SetRangeUser(iCharge, lowerChargeBinLimit, upperChargeBinLimit, kTRUE);
    data->SetRangeUser(iCharge, lowerChargeBinLimit, upperChargeBinLimit, kTRUE);
  }
  
  std::cout << std::endl;
  
  
  // Load numJet histos from bhess_PID*.root file related to efficiency file.
  TString pathNameDataMC = pathNameEfficiency;
  pathNameDataMC.ReplaceAll("_efficiency", "");
  
  TFile* fDataMC = TFile::Open(pathNameDataMC.Data());
  if (!fDataMC)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameData.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
    
    return -1;
  }
  
  TString listName = pathNameDataMC;
  listName.Replace(0, listName.Last('/') + 1, "");
  listName.ReplaceAll(".root", "");
    
  TObjArray* histList = (TObjArray*)(fDataMC->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
    return -1;
  }
  
  hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
  hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
  
  if (!hNjetsRec || !hNjetsGen) {
    std::cout << "Failed to load number of jets histos!" << std::endl;
    
    
    // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
    // period is considered; also: No multiplicity information)
    TString pathEfficiency = pathNameEfficiency;
    pathEfficiency.Replace(pathEfficiency.Last('/'), pathEfficiency.Length(), "");
    TString pathBackward = Form("%s/AnalysisResults.root", pathEfficiency.Data());
    TFile* fBackward = TFile::Open(pathBackward.Data());
    
    TString dirDataInFile = "";
    TDirectory* dirData = fBackward ? (TDirectory*)fBackward->Get(fBackward->GetListOfKeys()->At(0)->GetName()) : 0x0;
  
    TList* list = dirData ? (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName()) : 0x0;

    TH1D* hFFJetPtRec = list ? (TH1D*)list->FindObject("fh1FFJetPtRecCutsInc") : 0x0;
    TH1D* hFFJetPtGen = list ? (TH1D*)list->FindObject("fh1FFJetPtGenInc") : 0x0;
    
    if (hFFJetPtRec && hFFJetPtGen) {
      printf("***WARNING: For backward compatibility, using file \"%s\" to get number of jets. BUT: Might be wrong period and has no mult info!***\n",
        pathBackward.Data());
      printf("ALSO: Using Njets for inclusive jets!!!!\n");
      
      hNjetsRec = new TH2D("fh2FFJetPtRec", "", 1, -1, 1, dataRebinned->GetAxis(iJetPt, 0)->GetNbins(),
                          dataRebinned->GetAxis(iJetPt, 0)->GetXbins()->GetArray());
      
      for (Int_t iJet = 1; iJet <= hNjetsRec->GetNbinsY(); iJet++) {
        Int_t lowerBin = hFFJetPtRec->FindBin(hNjetsRec->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
        Int_t upperBin = hFFJetPtRec->FindBin(hNjetsRec->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
        hNjetsRec->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
      }
      
      hNjetsGen = new TH2D("fh2FFJetPtGen", "", 1, -1, 1,  dataRebinned->GetAxis(iJetPt, 0)->GetNbins(),
                          dataRebinned->GetAxis(iJetPt, 0)->GetXbins()->GetArray());
      
      for (Int_t iJet = 1; iJet <= hNjetsGen->GetNbinsY(); iJet++) {
        Int_t lowerBin = hFFJetPtGen->FindBin(hNjetsGen->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
        Int_t upperBin = hFFJetPtGen->FindBin(hNjetsGen->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
        hNjetsGen->SetBinContent(1, iJet, hFFJetPtGen->Integral(lowerBin, upperBin));
      }
    }
    
    if (!hNjetsRec || ! hNjetsGen)
      return -1;
  }

  // For normalisation to number of jets
  // NOTE: These numbers are for the efficiency only! The data will be normalised to its own number!!!
  const Double_t nJetsGen = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                            upperJetPtBinLimit) : 1.;
  const Double_t nJetsRec = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, lowerJetPtBinLimit,
                                                            upperJetPtBinLimit) : 1.;
  
  
  
  // Secondary correction
  AliCFEffGrid* sec = new AliCFEffGrid("sec", "Secondary Contamination", *dataRebinned);
  sec->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObs);
  
  TH1 *hSecAllPt = 0x0, *hSecID2Pt = 0x0;
  TH1D* hSecPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecToPiRatioPt[AliPID::kSPECIES] = { 0x0, };
  TCanvas *cSecPt = 0x0;
  TCanvas *cSecToPiRatioPt = 0x0;
  
  TH1 *hSecAllZ = 0x0, *hSecID2Z = 0x0;
  TH1D* hSecZ[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecToPiRatioZ[AliPID::kSPECIES] = { 0x0, };
  TCanvas *cSecZ = 0x0;
  TCanvas *cSecToPiRatioZ = 0x0;
  
  TH1 *hSecAllXi = 0x0, *hSecID2Xi = 0x0;
  TH1D* hSecXi[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecToPiRatioXi[AliPID::kSPECIES] = { 0x0, };
  TCanvas *cSecXi = 0x0;
  TCanvas *cSecToPiRatioXi = 0x0;
  
  // Secondary correction with strangeness scaling
  AliCFEffGrid* secStrangeScale = new AliCFEffGrid("secStrangeScale", "Secondary Contamination with strangeness scaling",
                                                   *dataRebinned);
  secStrangeScale->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepRecWithRecCutsMeasuredObsStrangenessScaled);
  
  TH1 *hSecSSAllPt = 0x0, *hSecSSID2Pt = 0x0;
  TH1D* hSecSSPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecSSToPiRatioPt[AliPID::kSPECIES] = { 0x0, };
  TCanvas *cSecSSPt = 0x0;
  TCanvas *cSecSSToPiRatioPt = 0x0;
  
  TH1 *hSecSSAllZ = 0x0, *hSecSSID2Z = 0x0;
  TH1D* hSecSSZ[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecSSToPiRatioZ[AliPID::kSPECIES] = { 0x0, };
  TCanvas *cSecSSZ = 0x0;
  TCanvas *cSecSSToPiRatioZ = 0x0;
 
  TH1 *hSecSSAllXi = 0x0, *hSecSSID2Xi = 0x0;
  TH1D* hSecSSXi[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSecSSToPiRatioXi[AliPID::kSPECIES] = { 0x0, };
  TCanvas *cSecSSXi = 0x0;
  TCanvas *cSecSSToPiRatioXi = 0x0;
  
  
  const Double_t upperTrackPt = restrictJetPtAxis ? actualUpperJetPt : 50.;
  
  TString strangenessString[2] = { "", "SS" };
  
  // Get the secondary contamination vs. pT,z,xi for each species w/ and w/o strangeness scaling
  
  for (Int_t type = 0; type < kNtypes; type++) {
    for (Int_t strangenessScaling = 0; strangenessScaling < 2; strangenessScaling++) {
      TH1 **hSecAll, **hSecID2;
      TH1D **hSec, **hSecToPiRatio;
      TCanvas **cSec, **cSecToPiRatio;
      
      Int_t iDimProj = -1;
      
      AliCFEffGrid* secCurr = strangenessScaling ? secStrangeScale : sec;
    
      if (type == kTrackPt) {
        hSecAll = strangenessScaling ? &hSecSSAllPt : &hSecAllPt;
        hSecID2 = strangenessScaling ? &hSecSSID2Pt : &hSecID2Pt;
        hSec = strangenessScaling ? &hSecSSPt[0] : &hSecPt[0];
        hSecToPiRatio = strangenessScaling ? &hSecSSToPiRatioPt[0] : &hSecToPiRatioPt[0];
        cSec = strangenessScaling ? &cSecSSPt : &cSecPt;
        cSecToPiRatio = strangenessScaling ? &cSecSSToPiRatioPt : &cSecToPiRatioPt;
        iDimProj = iPt;
      }
      else if (type == kZ) {
        hSecAll = strangenessScaling ? &hSecSSAllZ : &hSecAllZ;
        hSecID2 = strangenessScaling ? &hSecSSID2Z : &hSecID2Z;
        hSec = strangenessScaling ? &hSecSSZ[0] : &hSecZ[0];
        hSecToPiRatio = strangenessScaling ? &hSecSSToPiRatioZ[0] : &hSecToPiRatioZ[0];
        cSec = strangenessScaling ? &cSecSSZ : &cSecZ;
        cSecToPiRatio = strangenessScaling ? &cSecSSToPiRatioZ : &cSecToPiRatioZ;
        iDimProj = iZ;
      }
      else if (type == kXi) {
        hSecAll = strangenessScaling ? &hSecSSAllXi : &hSecAllXi;
        hSecID2 = strangenessScaling ? &hSecSSID2Xi : &hSecID2Xi;
        hSec = strangenessScaling ? &hSecSSXi[0] : &hSecXi[0];
        hSecToPiRatio = strangenessScaling ? &hSecSSToPiRatioXi[0] : &hSecToPiRatioXi[0];
        cSec = strangenessScaling ? &cSecSSXi : &cSecXi;
        cSecToPiRatio = strangenessScaling ? &cSecSSToPiRatioXi : &cSecToPiRatioXi;
        iDimProj = iXi;
      }
      else
        continue;
      
      *hSecAll = secCurr->Project(iDimProj); 
      (*hSecAll)->SetName(Form("hSec%s_%s_all", strangenessString[strangenessScaling].Data(), titles[type].Data()));
      (*hSecAll)->SetTitle("All");
      (*hSecAll)->SetMarkerStyle(24);
      (*hSecAll)->SetLineWidth(2.);
      (*hSecAll)->SetStats(kFALSE);
      (*hSecAll)->GetYaxis()->SetTitle("Primary Fraction");
      (*hSecID2) = (TH2D*)secCurr->Project(iDimProj, iMCid);
      (*hSecID2)->SetName(Form("hSecID2_%s", titles[type].Data()));
      
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hSec[species] = ((TH2D*)*hSecID2)->ProjectionX(Form("hSec%s_%s_%s", strangenessString[strangenessScaling].Data(),
                                                            titles[type].Data(), AliPID::ParticleShortName(species)),
                                                      species + 1, species + 1, "e");
        hSec[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
        hSec[species]->SetLineColor(getLineColorAliPID(species));
        hSec[species]->SetMarkerColor(getLineColorAliPID(species));
        hSec[species]->SetMarkerStyle(24);
        hSec[species]->SetLineWidth(2.);
        hSec[species]->GetYaxis()->SetRangeUser(0., 1.01);
        hSec[species]->SetStats(kFALSE);
        hSec[species]->GetYaxis()->SetTitle("Primary Fraction");
      }
      
      *cSec = new TCanvas(Form("cSec%s_%s", strangenessString[strangenessScaling].Data(), titles[type].Data()),
                          Form("Primary fraction for different species%s", strangenessScaling ? " (strangeness scaled)" : ""),
                          100, 10, 1200, 800);
      (*cSec)->cd();
      (*cSec)->SetGridx(1);
      (*cSec)->SetGridy(1);
      
      if (type == kTrackPt)
        (*hSecAll)->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
      
      (*hSecAll)->Draw("E1");
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (type == kTrackPt)
          hSec[i]->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
        
        hSec[i]->Draw("E1 same");
      }
      (*cSec)->BuildLegend()->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(*cSec);
      
      // To-pi ratios
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (species == AliPID::kPion)
          continue; // Do not consider pion-to-pion ratio
          
        hSecToPiRatio[species] = new TH1D(*hSec[species]);
        hSecToPiRatio[species]->Reset();
        hSecToPiRatio[species]->SetName(Form("hSec%sToPionRatio_%s_%s", strangenessString[strangenessScaling].Data(),
                                             titles[type].Data(), AliPID::ParticleShortName(species)));
        hSecToPiRatio[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
        hSecToPiRatio[species]->SetLineColor(getLineColorAliPID(species));
        hSecToPiRatio[species]->SetMarkerColor(getLineColorAliPID(species));
        hSecToPiRatio[species]->SetMarkerStyle(24);
        hSecToPiRatio[species]->SetLineWidth(2.);
        hSecToPiRatio[species]->SetStats(kFALSE);
        hSecToPiRatio[species]->GetYaxis()->SetTitle("Primary Fraction of Ratio");
        
        // Samples for different species are independent, so just divide correction factors
        hSecToPiRatio[species]->Divide(hSec[species], hSec[AliPID::kPion], 1., 1., ""); 
        hSecToPiRatio[species]->GetYaxis()->SetRangeUser(0., 1.1);
      }
      
      *cSecToPiRatio = new TCanvas(Form("cSec%sToPiRatio_%s", strangenessString[strangenessScaling].Data(), titles[type].Data()),
                                   Form("Primary fraction of to-#pi-ratio for different species%s",
                                        strangenessScaling ? " (strangeness scaled)" : ""),
                                   100, 10, 1200, 800);
      (*cSecToPiRatio)->cd();
      (*cSecToPiRatio)->SetGridx(1);
      (*cSecToPiRatio)->SetGridy(1);
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (i == AliPID::kPion)
          continue;
        
        if (type == kTrackPt)
          hSecToPiRatio[i]->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
        
        hSecToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
      }
      (*cSecToPiRatio)->BuildLegend()->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(*cSecToPiRatio);
    }
  }
  
  
  // Construct the efficiency grid from the data container 
  // NOTE: Here, the scale factors 1/nJet are different. This effects the error of the ratio.
  // Thus, the histograms must be extracted, scaled and THEN divided to arrive at the correction
  // factor and the standard "calculateEfficiency" function CANNOT be used (multiplying the ratio
  // of the scale factor afterwards would result in wrong errors)
  
  
  
  // If desired, apply the GEANT/FLUKA correction first
  // NOTE: This will change dataRebinned! If anything else should also be calculated, it must be done before.
  // Otherwise, the GEANT/FLUKA correction factor for the efficiency will also affect it!
  const Int_t genStepEff = kStepGenWithGenCuts;
  if (correctGeantFluka) {
    printf("Applying GEANT/FLUKA correction...\n");
    if (!geantFlukaCorrection(dataRebinned, genStepEff)) {
      printf("GEANT/FLUKA correction could not be applied!\n");
      return kFALSE;
    }
  }
  
  
  // Either one can take kStepRecWithGenCutsMeasuredObs or, what I prefer, one can take
  // kStepRecWithRecCutsMeasuredObsPrimaries => The difference is only the eta cut, which is on the rec level
  // in the latter case, i.e. one corrects for eta resolution (although the effect should be very small)
  const Int_t recStepEff = kStepRecWithRecCutsMeasuredObsPrimaries;
  
  // The efficiency vs. pT, z, xi
  TH1* hSpecPt_RecPrimaries = 0x0;
  TH1* hSpecPt_GenPrimaries = 0x0;
  TH1D* hEfficiencyAllPt = 0x0;
  
  TH2* hSpecID2Pt_RecPrimaries = 0x0;
  TH2* hSpecID2Pt_GenPrimaries = 0x0;
  TH2D* hEffID2Pt = 0x0;
  
  TH1D* hSpecIDPt_GenPrimaries[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSpecIDPt_RecPrimaries[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatioPt[AliPID::kSPECIES] = { 0x0, };
  
  TCanvas *cSpecPt = 0x0, *cEffPt = 0x0, *cEffToPiRatioPt = 0x0;
  
  
  TH1* hSpecZ_RecPrimaries = 0x0;
  TH1* hSpecZ_GenPrimaries = 0x0;
  TH1D* hEfficiencyAllZ = 0x0;
  
  TH2* hSpecID2Z_RecPrimaries = 0x0;
  TH2* hSpecID2Z_GenPrimaries = 0x0;
  TH2D* hEffID2Z = 0x0;
  
  TH1D* hSpecIDZ_GenPrimaries[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSpecIDZ_RecPrimaries[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyZ[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatioZ[AliPID::kSPECIES] = { 0x0, };
  
  TCanvas *cSpecZ = 0x0, *cEffZ = 0x0, *cEffToPiRatioZ = 0x0;
  
  
  TH1* hSpecXi_RecPrimaries = 0x0;
  TH1* hSpecXi_GenPrimaries = 0x0;
  TH1D* hEfficiencyAllXi = 0x0;
  
  TH2* hSpecID2Xi_RecPrimaries = 0x0;
  TH2* hSpecID2Xi_GenPrimaries = 0x0;
  TH2D* hEffID2Xi = 0x0;
  
  TH1D* hSpecIDXi_GenPrimaries[AliPID::kSPECIES] = { 0x0, };
  TH1D* hSpecIDXi_RecPrimaries[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyXi[AliPID::kSPECIES] = { 0x0, };
  TH1D* hEfficiencyToPiRatioXi[AliPID::kSPECIES] = { 0x0, };

  TCanvas *cSpecXi = 0x0, *cEffXi = 0x0, *cEffToPiRatioXi = 0x0;
  
  TString titleYaxis[kNtypes];
  titleYaxis[kTrackPt] = "1/N_{Jets} dN/dP_{T} (GeV/c)^{-1}";
  titleYaxis[kZ] = "1/N_{Jets} dN/dz";
  titleYaxis[kXi] = "1/N_{Jets} dN/d#xi";
  
  for (Int_t type = 0; type < kNtypes; type++) {
    TH1** hSpec_RecPrimaries;
    TH1** hSpec_GenPrimaries;
    TH1D** hEfficiencyAll;
    
    TH2** hSpecID2_RecPrimaries;
    TH2** hSpecID2_GenPrimaries;
    TH2D** hEffID2;
    
    TH1D** hSpecID_GenPrimaries;
    TH1D** hSpecID_RecPrimaries;
    TH1D** hEfficiency;
    TH1D** hEfficiencyToPiRatio;
    
    TCanvas **cSpec, **cEff, **cEffToPiRatio;
    
    Int_t iDimProj = -1;
    
  
    if (type == kTrackPt) {
      hSpec_RecPrimaries = &hSpecPt_RecPrimaries;
      hSpec_GenPrimaries = &hSpecPt_GenPrimaries;
      hEfficiencyAll = &hEfficiencyAllPt;
      
      hSpecID2_RecPrimaries = &hSpecID2Pt_RecPrimaries;
      hSpecID2_GenPrimaries = &hSpecID2Pt_GenPrimaries;
      hEffID2 = &hEffID2Pt;
      
      hSpecID_RecPrimaries = &hSpecIDPt_RecPrimaries[0];
      hSpecID_GenPrimaries = &hSpecIDPt_GenPrimaries[0];
      hEfficiency = &hEfficiencyPt[0];
      hEfficiencyToPiRatio = &hEfficiencyToPiRatioPt[0];
      
      cSpec = &cSpecPt;
      cEff = &cEffPt;
      cEffToPiRatio = &cEffToPiRatioPt;
      
      iDimProj = iPt;
    }
    else if (type == kZ) {
      hSpec_RecPrimaries = &hSpecZ_RecPrimaries;
      hSpec_GenPrimaries = &hSpecZ_GenPrimaries;
      hEfficiencyAll = &hEfficiencyAllZ;
      
      hSpecID2_RecPrimaries = &hSpecID2Z_RecPrimaries;
      hSpecID2_GenPrimaries = &hSpecID2Z_GenPrimaries;
      hEffID2 = &hEffID2Z;
      
      hSpecID_RecPrimaries = &hSpecIDZ_RecPrimaries[0];
      hSpecID_GenPrimaries = &hSpecIDZ_GenPrimaries[0];
      hEfficiency = &hEfficiencyZ[0];
      hEfficiencyToPiRatio = &hEfficiencyToPiRatioZ[0];
      
      cSpec = &cSpecZ;
      cEff = &cEffZ;
      cEffToPiRatio = &cEffToPiRatioZ;
      
      iDimProj = iZ;
    }
    else if (type == kXi) {
      hSpec_RecPrimaries = &hSpecXi_RecPrimaries;
      hSpec_GenPrimaries = &hSpecXi_GenPrimaries;
      hEfficiencyAll = &hEfficiencyAllXi;
      
      hSpecID2_RecPrimaries = &hSpecID2Xi_RecPrimaries;
      hSpecID2_GenPrimaries = &hSpecID2Xi_GenPrimaries;
      hEffID2 = &hEffID2Xi;
      
      hSpecID_RecPrimaries = &hSpecIDXi_RecPrimaries[0];
      hSpecID_GenPrimaries = &hSpecIDXi_GenPrimaries[0];
      hEfficiency = &hEfficiencyXi[0];
      hEfficiencyToPiRatio = &hEfficiencyToPiRatioXi[0];
      
      cSpec = &cSpecXi;
      cEff = &cEffXi;
      cEffToPiRatio = &cEffToPiRatioXi;
      
      iDimProj = iXi;
    }
    else
      continue;
    
    
    // Divide steps: recStepEff (=kStepRecWithRecCutsMeasuredObsPrimaries) / genStepEff (=kStepGenWithGenCuts)
    *hSpec_RecPrimaries = dataRebinned->GetGrid(recStepEff)->Project(iDimProj);  
    setupHist(*hSpec_RecPrimaries, Form("hSpec_%s_RecPrimaries", titles[type].Data()), "All, rec", "",
              titleYaxis[type].Data(), kBlack);
    normaliseHist(*hSpec_RecPrimaries, nJetsRec > 0 ? 1. / nJetsRec : 0.);
    (*hSpec_RecPrimaries)->SetLineWidth(2.);
    (*hSpec_RecPrimaries)->SetMarkerStyle(20);
    
    *hSpec_GenPrimaries = dataRebinned->GetGrid(genStepEff)->Project(iDimProj);
    setupHist(*hSpec_GenPrimaries, Form("hSpec_%s_GenPrimaries", titles[type].Data()), "All, gen", "",
              titleYaxis[type].Data(), kBlack);
    normaliseHist(*hSpec_GenPrimaries, nJetsGen > 0 ? 1. / nJetsGen : 0.);
    (*hSpec_GenPrimaries)->SetLineWidth(2.);
    (*hSpec_GenPrimaries)->SetLineStyle(2);
    
    *hEfficiencyAll = new TH1D(*(TH1D*)(*hSpec_RecPrimaries));
    (*hEfficiencyAll)->SetName(Form("hEfficiencyAll_%s", titles[type].Data()));
    (*hEfficiencyAll)->GetYaxis()->SetTitle("");
    (*hEfficiencyAll)->SetTitle("All");
    (*hEfficiencyAll)->SetLineWidth(2.0);
    (*hEfficiencyAll)->SetStats(kFALSE);
    (*hEfficiencyAll)->GetYaxis()->SetTitle("Efficiency x Acceptance x pT Resolution");
    (*hEfficiencyAll)->GetYaxis()->SetTitleSize(0.05);
    (*hEfficiencyAll)->Divide(*hSpec_RecPrimaries, *hSpec_GenPrimaries, 1., 1., "B");
    (*hEfficiencyAll)->GetYaxis()->SetRangeUser(0., 2.0);
    
    *hSpecID2_RecPrimaries = (TH2*)dataRebinned->GetGrid(recStepEff)->Project(iDimProj, iMCid);  
    (*hSpecID2_RecPrimaries)->SetName(Form("hSpecID2_RecPrimaries_%s", titles[type].Data()));
    (*hSpecID2_RecPrimaries)->Scale(nJetsRec > 0 ? 1. / nJetsRec : 0.);
    
    *hSpecID2_GenPrimaries = (TH2*)dataRebinned->GetGrid(genStepEff)->Project(iDimProj, iMCid); 
    (*hSpecID2_GenPrimaries)->SetName(Form("hSpecID2_GenPrimaries_%s", titles[type].Data()));
    (*hSpecID2_GenPrimaries)->Scale(nJetsGen > 0 ? 1. / nJetsGen : 0.);
    
    *hEffID2 = new TH2D(*(TH2D*)(*hSpecID2_RecPrimaries));
    (*hEffID2)->SetName(Form("hEff_%s", titles[type].Data()));
    (*hEffID2)->SetTitle("All");
    (*hEffID2)->Divide(*hSpecID2_RecPrimaries, *hSpecID2_GenPrimaries, 1., 1., "B");

    // Get the efficiencies vs. pT,z,xi for each species
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hEfficiency[species] = ((TH2D*)*hEffID2)->ProjectionX(Form("hEfficiency_%s_%s", titles[type].Data(),
                                                                 AliPID::ParticleShortName(species)), 
                                                            species + 1, species + 1, "e");
      setupHist(hEfficiency[species], "", Form("%s", AliPID::ParticleLatexName(species)), "", "",
                getLineColorAliPID(species));
      hEfficiency[species]->SetLineWidth(2.);
      hEfficiency[species]->GetYaxis()->SetRangeUser(0., 2.0);
      hEfficiency[species]->SetStats(kFALSE);
      hEfficiency[species]->GetYaxis()->SetTitle("Efficiency x Acceptance x pT Resolution");
      hEfficiency[species]->GetYaxis()->SetTitleSize(0.05);
      
      hSpecID_GenPrimaries[species] = ((TH2D*)*hSpecID2_GenPrimaries)->ProjectionX(Form("hSpecID_%s_gen_%s",
                                                                                        titles[type].Data(),
                                                                                        AliPID::ParticleShortName(species)),
                                                                                        species + 1, species + 1, "e");
      setupHist(hSpecID_GenPrimaries[species], "", Form("%s, gen", AliPID::ParticleLatexName(species)), "",
                titleYaxis[type].Data(), getLineColorAliPID(species));
      normaliseHist(hSpecID_GenPrimaries[species], 1.);
      hSpecID_GenPrimaries[species]->SetLineWidth(2.);
      hSpecID_GenPrimaries[species]->SetLineStyle(2.);
      hSpecID_GenPrimaries[species]->SetStats(kFALSE);
      
      
      hSpecID_RecPrimaries[species] = ((TH2D*)*hSpecID2_RecPrimaries)->ProjectionX(Form("hSpecID_%s_rec_%s",
                                                                                        titles[type].Data(),
                                                                                        AliPID::ParticleShortName(species)),
                                                                                        species + 1, species + 1, "e");
      setupHist(hSpecID_RecPrimaries[species], "", Form("%s, rec", AliPID::ParticleLatexName(species)), "",
                titleYaxis[type].Data(), getLineColorAliPID(species));
      normaliseHist(hSpecID_RecPrimaries[species], 1.);
      hSpecID_RecPrimaries[species]->SetLineWidth(2.);
      hSpecID_RecPrimaries[species]->SetMarkerStyle(20);
      hSpecID_RecPrimaries[species]->SetStats(kFALSE);
    }

    *cSpec = new TCanvas(Form("cSpec_%s", titles[type].Data()), "Spectra for different species", 0, 300, 900, 900);
    (*cSpec)->cd();
    (*cSpec)->SetGridx(1);
    (*cSpec)->SetGridy(1);
    (*cSpec)->SetLogy(1);
    
    if (type == kTrackPt)
      (*hSpec_GenPrimaries)->GetYaxis()->SetRangeUser(1e-8, 1e2);
    else if (type == kZ)
      (*hSpec_GenPrimaries)->GetYaxis()->SetRangeUser(1e-5, 1e3);
    else if (type == kXi)
      (*hSpec_GenPrimaries)->GetYaxis()->SetRangeUser(1e-5, 1e3);
    
    
    if (type == kTrackPt)
      (*hSpec_GenPrimaries)->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
    (*hSpec_GenPrimaries)->Draw("E1");
    
    if (type == kTrackPt)
      (*hSpec_RecPrimaries)->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
    (*hSpec_RecPrimaries)->Draw("E1 same");
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (type == kTrackPt) {
        hSpecID_GenPrimaries[i]->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
        hSpecID_RecPrimaries[i]->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
      }
      
      hSpecIDPt_GenPrimaries[i]->Draw("E1 same");
      hSpecIDPt_RecPrimaries[i]->Draw("E1 same");
    }
    TLegend* leg = (*cSpec)->BuildLegend();
    leg->SetFillColor(kWhite);
    leg->SetNColumns(2);
    
    ClearTitleFromHistoInCanvas(*cSpec);
    
    
    *cEff = new TCanvas(Form("cEff_%s", titles[type].Data()), "Efficiency x Acceptance x pT Resolution for different species",
                        0, 300, 900, 900);
    (*cEff)->cd();
    (*cEff)->SetGridx(1);
    (*cEff)->SetGridy(1);
    
    if (type == kTrackPt)
      (*hEfficiencyAll)->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
    (*hEfficiencyAll)->Draw("E1");
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (type == kTrackPt)
        hEfficiency[i]->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
      hEfficiency[i]->Draw("E1 same");
    }
    (*cEff)->BuildLegend()->SetFillColor(kWhite);
    
    ClearTitleFromHistoInCanvas(*cEff);
    
    
    // To-pi ratios
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (species == AliPID::kPion)
        continue; // Do not consider pion-to-pion ratio
        
      hEfficiencyToPiRatio[species] = new TH1D(*hEfficiency[species]);
      hEfficiencyToPiRatio[species]->Reset();
      hEfficiencyToPiRatio[species]->SetName(Form("hEfficiency_ToPionRatio_%s_%s", titles[type].Data(),
                                                  AliPID::ParticleShortName(species)));
      
      setupHist(hEfficiencyToPiRatio[species], "", Form("%s/#pi", AliPID::ParticleLatexName(species)), "", "",
                getLineColorAliPID(species));
      hEfficiencyToPiRatio[species]->SetLineWidth(2.);
      hEfficiencyToPiRatio[species]->SetStats(kFALSE);
      hEfficiencyToPiRatio[species]->GetYaxis()->SetTitle("Efficiency x Acceptance x pT Resolution");
      hEfficiencyToPiRatio[species]->GetYaxis()->SetTitleSize(0.05);
      
      // Samples for different species are independent, so just divide correction factors
      hEfficiencyToPiRatio[species]->Divide(hEfficiency[species], hEfficiency[AliPID::kPion], 1., 1., ""); 
      hEfficiencyToPiRatio[species]->GetYaxis()->SetRangeUser(0., 2.0);
    }
    
    *cEffToPiRatio = new TCanvas(Form("cEffToPiRatio_%s", titles[type].Data()),
                                 "Efficiency x Acceptance x pT Resolution for different to-#pi-ratios",
                                 100, 10, 1200, 800);
    (*cEffToPiRatio)->cd();
    (*cEffToPiRatio)->SetGridx(1);
    (*cEffToPiRatio)->SetGridy(1);
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i == AliPID::kPion)
        continue;
      
      if (type == kTrackPt)
        hEfficiencyToPiRatio[i]->GetXaxis()->SetRangeUser(0.15, upperTrackPt);
      
      hEfficiencyToPiRatio[i]->Draw(Form("E1%s", i == 0 ? "" : " same"));
    }
    (*cEffToPiRatio)->BuildLegend()->SetFillColor(kWhite);
    
    ClearTitleFromHistoInCanvas(*cEffToPiRatio);
  }
  
  // Save results to file
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  TString saveFileName = pathNameData;
  saveFileName.ReplaceAll(Form("%s/", pathData.Data()), "");
  saveFileName.Prepend("output_EfficiencyCorrection_");
  saveFileName.ReplaceAll(".root", Form("_jetPt%.1f_%.1f%s.root", actualLowerJetPt, actualUpperJetPt, chargeString.Data()));
  
  TString saveFilePathName = Form("%s/%s", pathData.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  
  if (!saveFile) {
    printf("Could not save results!\n");
    
    return -1;
  }
  
  saveFile->cd();
  
  if (cSecPt)
    cSecPt->Write();
  
  if (cSecSSPt)
    cSecSSPt->Write();
  
  if (hSecAllPt)
      hSecAllPt->Write();
  
  if (hSecSSAllPt)
      hSecSSAllPt->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecPt[i])
      hSecPt[i]->Write();
    
    if (hSecSSPt[i])
      hSecSSPt[i]->Write();
  }
  
  if (cSecToPiRatioPt)
    cSecToPiRatioPt->Write();
  
  if (cSecSSToPiRatioPt)
    cSecSSToPiRatioPt->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecToPiRatioPt[i])
      hSecToPiRatioPt[i]->Write();
    
    if (hSecSSToPiRatioPt[i])
      hSecSSToPiRatioPt[i]->Write();
  }

  
  
  if (cSecZ)
    cSecZ->Write();
  
  if (cSecSSZ)
    cSecSSZ->Write();
  
  if (hSecAllZ)
      hSecAllZ->Write();
  
  if (hSecSSAllZ)
      hSecSSAllZ->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecZ[i])
      hSecZ[i]->Write();
    
    if (hSecSSZ[i])
      hSecSSZ[i]->Write();
  }
  
  if (cSecToPiRatioZ)
    cSecToPiRatioZ->Write();
  
  if (cSecSSToPiRatioZ)
    cSecSSToPiRatioZ->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecToPiRatioZ[i])
      hSecToPiRatioZ[i]->Write();
    
    if (hSecSSToPiRatioZ[i])
      hSecSSToPiRatioZ[i]->Write();
  }
  
  
  
  if (cSecXi)
    cSecXi->Write();
  
  if (cSecSSXi)
    cSecSSXi->Write();
  
  if (hSecAllXi)
      hSecAllXi->Write();
  
  if (hSecSSAllXi)
      hSecSSAllXi->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecXi[i])
      hSecXi[i]->Write();
    
    if (hSecSSXi[i])
      hSecSSXi[i]->Write();
  }
  
  if (cSecToPiRatioXi)
    cSecToPiRatioXi->Write();
  
  if (cSecSSToPiRatioXi)
    cSecSSToPiRatioXi->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSecToPiRatioXi[i])
      hSecToPiRatioXi[i]->Write();
    
    if (hSecSSToPiRatioXi[i])
      hSecSSToPiRatioXi[i]->Write();
  }
  
  
  
  if (cSpecPt)
    cSpecPt->Write();
  
  if (hSpecPt_RecPrimaries)
    hSpecPt_RecPrimaries->Write();
  
  if (hSpecPt_GenPrimaries)
    hSpecPt_GenPrimaries->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSpecIDPt_GenPrimaries[i])
      hSpecIDPt_GenPrimaries[i]->Write();
    
    if (hSpecIDPt_RecPrimaries[i])
      hSpecIDPt_RecPrimaries[i]->Write();
  }
  
  
  if (cEffPt)
    cEffPt->Write();
  
  if (hEfficiencyAllPt)
      hEfficiencyAllPt->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyPt[i])
      hEfficiencyPt[i]->Write();
  }
  
  if (cEffToPiRatioPt)
    cEffToPiRatioPt->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyToPiRatioPt[i])
      hEfficiencyToPiRatioPt[i]->Write();
  }
  
  
  
  if (cSpecZ)
    cSpecZ->Write();
  
  if (hSpecZ_RecPrimaries)
    hSpecZ_RecPrimaries->Write();
  
  if (hSpecZ_GenPrimaries)
    hSpecZ_GenPrimaries->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSpecIDZ_GenPrimaries[i])
      hSpecIDZ_GenPrimaries[i]->Write();
    
    if (hSpecIDZ_RecPrimaries[i])
      hSpecIDZ_RecPrimaries[i]->Write();
  }
  
  
  if (cEffZ)
    cEffZ->Write();
  
  if (hEfficiencyAllZ)
      hEfficiencyAllZ->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyZ[i])
      hEfficiencyZ[i]->Write();
  }
  
  if (cEffToPiRatioZ)
    cEffToPiRatioZ->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyToPiRatioZ[i])
      hEfficiencyToPiRatioZ[i]->Write();
  }
  
  
  
  if (cSpecXi)
    cSpecXi->Write();
  
  if (hSpecXi_RecPrimaries)
    hSpecXi_RecPrimaries->Write();
  
  if (hSpecXi_GenPrimaries)
    hSpecXi_GenPrimaries->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hSpecIDXi_GenPrimaries[i])
      hSpecIDXi_GenPrimaries[i]->Write();
    
    if (hSpecIDXi_RecPrimaries[i])
      hSpecIDXi_RecPrimaries[i]->Write();
  }
  
  
  if (cEffXi)
    cEffXi->Write();
  
  if (hEfficiencyAllXi)
      hEfficiencyAllXi->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyXi[i])
      hEfficiencyXi[i]->Write();
  }
  
  if (cEffToPiRatioXi)
    cEffToPiRatioXi->Write();
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hEfficiencyToPiRatioXi[i])
      hEfficiencyToPiRatioXi[i]->Write();
  }
  
  
  
  // Correct the yields with the efficiencies and primary fractions
  for (Int_t type = 0; type < kNtypes; type++) {
    TH1D* hYieldCorrected[AliPID::kSPECIES];
    TH1D* hYieldCorrectedSysError[AliPID::kSPECIES];
    
    TH1D** hYield;
    TH1D** hYieldSysError;
    TH1D** hMCgenPrimYield;
    TH1D** hSecSS;
    TH1D** hSec;
    TH1D** hEfficiency;
    
    if (type == kTrackPt) {
      hYield = &hYieldPt[0];
      hYieldSysError = &hYieldPtSysError[0];
      hMCgenPrimYield = &hMCgenPrimYieldPt[0];
      hSecSS = &hSecSSPt[0];
      hSec = &hSecPt[0];
      hEfficiency = &hEfficiencyPt[0];
    }
    else if (type == kZ) {
      hYield = &hYieldZ[0];
      hYieldSysError = &hYieldZSysError[0];
      hMCgenPrimYield = &hMCgenPrimYieldZ[0];
      hSecSS = &hSecSSZ[0];
      hSec = &hSecZ[0];
      hEfficiency = &hEfficiencyZ[0];
    }
    else if (type == kXi) {
      hYield = &hYieldXi[0];
      hYieldSysError = &hYieldXiSysError[0];
      hMCgenPrimYield = &hMCgenPrimYieldXi[0];
      hSecSS = &hSecSSXi[0];
      hSec = &hSecXi[0];
      hEfficiency = &hEfficiencyXi[0];
    }
    else
      continue;
    
    
    // NOTE: The yields have already been normalised to nJets in extractFFs. However, the binning in jetPt might be different now.
    // Thus, in order to get the correct result, the old normalisation (just the same binnnig as in the nJet histo) has to be 
    // undone (was done above) and then it has to be normalised to the number of jets in the new binning again (now).
    // Also NOTE: Histos are already normalised to bin width via extractFFs!
    const Double_t nJetsGenData = hNjetsGenData ? hNjetsGenData->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData, 
                                                                          lowerJetPtBinLimitData, upperJetPtBinLimitData) : 1.;
    const Double_t nJetsRecData = hNjetsRecData ? hNjetsRecData->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData,                 
                                                                          lowerJetPtBinLimitData, upperJetPtBinLimitData) : 1.;
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hYield[species]->Scale((nJetsRecData > 0) ? 1. / nJetsRecData : 0.);
      hYieldSysError[species]->Scale((nJetsRecData > 0) ? 1. / nJetsRecData : 0.);
        
      hYield[species]->SetStats(kFALSE);
      hYieldSysError[species]->SetStats(kFALSE);
      
      if (hMCgenPrimYield[species]) {
        hMCgenPrimYield[species]->Scale((nJetsGenData > 0) ? 1. / nJetsGenData : 0.);
        hMCgenPrimYield[species]->SetStats(kFALSE);
      }
      
      hYield[species]->GetYaxis()->SetTitle(titleYaxis[type].Data());
      if (hMCgenPrimYield[species]) {
        hMCgenPrimYield[species]->GetYaxis()->SetTitle(titleYaxis[type].Data());
        hMCgenPrimYield[species]->SetTitle(Form("%s (MC)", AliPID::ParticleLatexName(species)));
        hMCgenPrimYield[species]->SetMarkerStyle(24);
        hMCgenPrimYield[species]->SetLineColor(getLineColorAliPID(species));
        hMCgenPrimYield[species]->SetMarkerColor(getLineColorAliPID(species));
      }

      hYield[species]->SetLineColor(getLineColorAliPID(species));
      hYield[species]->SetMarkerColor(getLineColorAliPID(species));
      hYield[species]->SetMarkerStyle(20);
      
      hYieldSysError[species]->SetLineColor(getLineColorAliPID(species));
      hYieldSysError[species]->SetMarkerColor(getLineColorAliPID(species));
      hYieldSysError[species]->SetMarkerStyle(20);
      hYieldSysError[species]->SetFillStyle(0);
      
      hYieldCorrected[species] = new TH1D(*hYield[species]);
      hYieldCorrected[species]->SetName(Form("%s_corrected", hYield[species]->GetName()));
      hYieldCorrected[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
      
      hYieldCorrectedSysError[species] = new TH1D(*hYieldSysError[species]);
      hYieldCorrectedSysError[species]->SetName(Form("%s_corrected", hYieldSysError[species]->GetName()));
      hYieldCorrectedSysError[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
    
      // Correction factor histos can have different binning than data histos (constant factor above some threshold)
      // -> Need special functions to multiply and divide such histos
      multiplyHistsDifferentBinning(hYieldCorrected[species], scaleStrangeness ? hSecSS[species] : hSec[species], 1., 1.);
      divideHistsDifferentBinning(hYieldCorrected[species], hEfficiency[species], 1., 1.);
      
      multiplyHistsDifferentBinning(hYieldCorrectedSysError[species], scaleStrangeness ? hSecSS[species] : hSec[species], 1., 1.,
                                    kTRUE);
      divideHistsDifferentBinning(hYieldCorrectedSysError[species], hEfficiency[species], 1., 1., kTRUE);

      //hYieldCorrected[species]->Multiply(hYieldCorrected[species], scaleStrangeness ? hSecSS[species] : hSec[species], 1., 1.,
      //"");
      //hYieldCorrected[species]->Divide(hYieldCorrected[species], hEfficiency[species], 1., 1., ""); 
    }
    
    // Calculate the total corrected yield. The original total yield had no error (just the number of detected tracks in a pT bin),
    // but due to the correction there is some error for the total yield. Also the error of the fractions introduces uncertainties
    // for the yields of individual species
    TH1D* hYieldCorrectedTotal = new TH1D(*hYieldCorrected[0]);
    hYieldCorrectedTotal->SetLineColor(kBlack);
    hYieldCorrectedTotal->SetMarkerColor(kBlack);
    hYieldCorrectedTotal->SetMarkerStyle(20);
    hYieldCorrectedTotal->SetName(Form("hYieldCorrectedTotal_%s", titles[type].Data()));
    hYieldCorrectedTotal->SetTitle("Total");
    //hYieldCorrectedTotal->SetTitle("Total, secondary and efficiency x acceptance x pT resolution corrected");
    
    //TODO doesn't it require the covariance matrix here to get the error of the corrected total yield?
    // Because before the total count was taken into account just by the fit. Now the fractions and their errors
    // influence the MC correction (and its errors), so that there is some error on the corrected total yield.
    // However, the errors of the contributing yields are correlated!
    for (Int_t i = 1; i < AliPID::kSPECIES; i++)
      hYieldCorrectedTotal->Add(hYieldCorrected[i], 1.);
    
    // Calculate the corrected fractions
    TH1D* hFractionCorrected[AliPID::kSPECIES];
    TH1D* hFractionCorrectedSysError[AliPID::kSPECIES];
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      hFractionCorrected[species] = new TH1D(*hYield[species]);
      TString oldName = hYield[species]->GetName();
      TString newName = oldName.ReplaceAll("Yield", "Fraction");
      hFractionCorrected[species]->SetName(Form("%s_corrected", newName.Data()));
      hFractionCorrected[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
      hFractionCorrected[species]->GetYaxis()->SetTitle("Corrected Fraction");
      hFractionCorrected[species]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hFractionCorrected[species]->GetXaxis()->SetNoExponent(kTRUE);
      hFractionCorrected[species]->SetFillStyle(0);

      
      // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
      // (numerator is a "selected" subset of the denominator). It doesn't matter that the error of the histos is not just 
      // "sqrt(content)" because the error formula also works for weighted histograms (which means that the error can be more
      // or less anything)
      hFractionCorrected[species]->Divide(hYieldCorrected[species], hYieldCorrectedTotal, 1., 1., "B"); 
      
      
      
      //  The systematic errors just need to be scaled in the same way as the fractions.
      // So, just take the ratio to the uncorrected fraction and scale the sys. error accordingly
      // or, in this case, just divide by the same total yield as for yield -> fractions
      hFractionCorrectedSysError[species] = new TH1D(*hFractionCorrected[species]);
      hFractionCorrectedSysError[species]->SetName(Form("%s_sysError", hFractionCorrected[species]->GetName()));
      hFractionCorrectedSysError[species]->SetTitle(Form("%s (sys. error)", hFractionCorrected[species]->GetTitle()));
      
      for (Int_t binX = 1; binX <= hFractionCorrectedSysError[species]->GetNbinsX(); binX++) {
        const Double_t corrTotalYield = hYieldCorrectedTotal->GetBinContent(binX);
        const Double_t scaleFactor = corrTotalYield > 0 ? 1.0 / corrTotalYield : 1.;
        hFractionCorrectedSysError[species]->SetBinError(binX,   hYieldCorrectedSysError[species]->GetBinError(binX) * scaleFactor);
      }
    }
    
    // If MC is available, calculate the generated fractions
    TH1D* hMCgenPrimYieldTotal = 0x0;
    TH1D* hMCgenPrimFraction[AliPID::kSPECIES];
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      hMCgenPrimFraction[i] = 0x0;
    
    if (numMCgenPrimYieldHistsFound > 0) {
      hMCgenPrimYieldTotal = new TH1D(*hMCgenPrimYield[0]);
      hMCgenPrimYieldTotal->SetLineColor(kBlack);
      hMCgenPrimYieldTotal->SetMarkerColor(kBlack);
      hMCgenPrimYieldTotal->SetMarkerStyle(24);
      hMCgenPrimYieldTotal->SetName(Form("hMCgenPrimYieldTotal%s", titles[type].Data()));
      hMCgenPrimYieldTotal->SetTitle("Total, MC truth");
      //hMCgenPrimYieldTotal->SetTitle("Total generated primary yield (MC truth)");
      
      for (Int_t i = 1; i < AliPID::kSPECIES; i++)
        hMCgenPrimYieldTotal->Add(hMCgenPrimYield[i], 1.);
      
      // Calculate the MC fractions
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hMCgenPrimYield[species]->SetLineColor(getLineColorAliPID(species));
        hMCgenPrimYield[species]->SetMarkerColor(getLineColorAliPID(species));
        hMCgenPrimYield[species]->SetMarkerStyle(24);
      
        hMCgenPrimYield[species]->SetTitle(Form("%s, MC truth", AliPID::ParticleLatexName(species)));
        
        hMCgenPrimFraction[species] = new TH1D(*hMCgenPrimYield[species]);
        TString oldName = hMCgenPrimYield[species]->GetName();
        TString newName = oldName.ReplaceAll("Yield", "Fraction");
        hMCgenPrimFraction[species]->SetName(newName.Data());

        // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
        // (numerator is a "selected" subset of the denominator).
        hMCgenPrimFraction[species]->Divide(hMCgenPrimFraction[species], hMCgenPrimYieldTotal, 1., 1., "B"); 
      }
    }
    
    TCanvas* cCorrData = new TCanvas(Form("cCorrData_%s", titles[type].Data()), "Corrected data", 0, 300, 900, 900);
    cCorrData->Divide(2, 1, 0., 0.01);
    cCorrData->GetPad(1)->SetLogx(1);
    cCorrData->GetPad(1)->SetLogy(1);
    cCorrData->GetPad(2)->SetLogx(1);
    cCorrData->GetPad(2)->SetLogy(1);
    
    cCorrData->GetPad(1)->SetRightMargin(0.001);
    cCorrData->GetPad(2)->SetRightMargin(0.001);
    
    cCorrData->GetPad(1)->SetLeftMargin(0.2);
    cCorrData->GetPad(2)->SetLeftMargin(0.2);
    
    cCorrData->cd(1); // uncorrected
    hYield[AliPID::kPion]->GetYaxis()->SetTitleOffset(1.4);
    hYield[AliPID::kPion]->Draw("E1");
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      hYieldSysError[i]->Draw("E2 same");
      
      if (i == AliPID::kPion) continue;
      hYield[i]->Draw("E1 same");
    }
    
    
    ClearTitleFromHistoInCanvas(cCorrData, 1);
    
    Double_t maximum = hYieldCorrectedTotal->GetBinContent(hYieldCorrectedTotal->GetMaximumBin()) * 10.;
    Double_t minimum = maximum;
    for (Int_t spec = 0; spec < AliPID::kSPECIES; spec++) {
      Double_t temp = hYieldCorrected[spec]->GetBinContent(hYieldCorrected[spec]->FindLastBinAbove(0.)) * 0.1;
      if (temp > 0)
        minimum = temp;
    }
    
    
    cCorrData->cd(2); // corrected
    hYieldCorrectedTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
    hYieldCorrectedTotal->GetXaxis()->SetNoExponent(kTRUE);
    hYieldCorrectedTotal->GetYaxis()->SetRangeUser(minimum, maximum);
    
    if (hMCgenPrimYieldTotal) {
      hMCgenPrimYieldTotal->GetYaxis()->SetTitleOffset(1.4);
      hMCgenPrimYieldTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
      hMCgenPrimYieldTotal->GetXaxis()->SetNoExponent(kTRUE);
      hMCgenPrimYieldTotal->GetXaxis()->SetTitle(hYieldCorrectedTotal->GetXaxis()->GetTitle());
      hMCgenPrimYieldTotal->GetYaxis()->SetRangeUser(minimum, maximum);
    }
    
    if (hMCgenPrimYieldTotal) {
      hMCgenPrimYieldTotal->Draw("E1");
      hYieldCorrectedTotal->Draw("E1 same");
    }
    else
      hYieldCorrectedTotal->Draw("E1");
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hMCgenPrimYield[i])
        hMCgenPrimYield[i]->Draw("E1 same");
      hYieldCorrected[i]->Draw("E1 same");
    }
    
    TLegend* legTemp = cCorrData->cd(2)->BuildLegend();//0.25, 0.16, 0.65, 0.51);
    legTemp->SetNColumns(2);
    legTemp->SetFillColor(kWhite);
    
    // Do not include into legend
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      hYieldCorrectedSysError[i]->Draw("E2 same");
    
    ClearTitleFromHistoInCanvas(cCorrData, 2);
    
    TCanvas* cCorrYieldsRatio = 0x0;
    
    TH1D* hYieldCorrectedRatioToMC[AliPID::kSPECIES];
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) 
      hYieldCorrectedRatioToMC[i] = 0x0;

    TH1D* hYieldCorrectedTotalRatioToMC = 0x0;
    
    if (numMCgenPrimYieldHistsFound > 0) {
      // Compare with MC truth
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hYieldCorrectedRatioToMC[species] = new TH1D(*hYieldCorrected[species]);
        hYieldCorrectedRatioToMC[species]->SetName(Form("hYieldCorrectedRatioToMC_%s_%s", titles[type].Data(), 
                                                        AliPID::ParticleShortName(species)));
        hYieldCorrectedRatioToMC[species]->SetTitle(Form("%s", AliPID::ParticleLatexName(species)));
        hYieldCorrectedRatioToMC[species]->GetYaxis()->SetTitle("Corrected Yield: Fit / MC Truth");
        hYieldCorrectedRatioToMC[species]->Divide(hMCgenPrimYield[species]);
      }
      
      hYieldCorrectedTotalRatioToMC = new TH1D(*hYieldCorrectedTotal);
      hYieldCorrectedTotalRatioToMC->SetName(Form("hYieldCorrectedTotalRatioToMC_%s", titles[type].Data()));
      hYieldCorrectedTotalRatioToMC->SetTitle("Total yield");
      hYieldCorrectedTotalRatioToMC->GetYaxis()->SetTitle("Corrected Yield: Fit / MC Truth");
      hYieldCorrectedTotalRatioToMC->Divide(hMCgenPrimYieldTotal);
      
      cCorrYieldsRatio = new TCanvas(Form("cCorrYieldsRatio_%s", titles[type].Data()), "Corrected Yields Comparison to MC",
                                     0, 300, 900, 900);
      cCorrYieldsRatio->SetGridx(1);
      cCorrYieldsRatio->SetGridy(1);
      cCorrYieldsRatio->SetLogx(1);
      
      hYieldCorrectedTotalRatioToMC->GetYaxis()->SetRangeUser(0.6, 1.6);
      hYieldCorrectedTotalRatioToMC->GetYaxis()->SetTitleOffset(0.85);
      hYieldCorrectedTotalRatioToMC->Draw("E1");
      
      for (Int_t species = 0; species < AliPID::kSPECIES; species++)
        hYieldCorrectedRatioToMC[species]->Draw("E1 same");
      
      cCorrYieldsRatio->BuildLegend()->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cCorrYieldsRatio);
    }
    
    TCanvas* cCorrFractions = new TCanvas(Form("cCorrFractions_%s", titles[type].Data()), "Corrected particleFractions",
                                          0, 300, 900, 900);
    cCorrFractions->SetLogx(1);
    hFractionCorrected[0]->GetYaxis()->SetRangeUser(0., 1.);
    hFractionCorrected[0]->Draw("E1");
    if (hMCgenPrimFraction[0])
      hMCgenPrimFraction[0]->Draw("E1 same");
    
    for (Int_t i = 1; i < AliPID::kSPECIES; i++) {
      hFractionCorrected[i]->Draw("E1 same");
      if (hMCgenPrimFraction[i])
        hMCgenPrimFraction[i]->Draw("E1 same");
    }
    
    cCorrFractions->BuildLegend()->SetFillColor(kWhite);
    
    // Do not include into legend
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      hFractionCorrectedSysError[i]->Draw("E2 same");
    
    ClearTitleFromHistoInCanvas(cCorrFractions);
    
    // Save results
    saveFile->cd();
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hYield[i])
        hYield[i]->Write();
    }
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hYieldCorrected[i])
        hYieldCorrected[i]->Write();
      
      if (hYieldCorrectedSysError[i])
        hYieldCorrectedSysError[i]->Write();
    }
    
    if (hYieldCorrectedTotal)
      hYieldCorrectedTotal->Write();
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hMCgenPrimYield[i])
        hMCgenPrimYield[i]->Write();
    }
    
    if (hMCgenPrimYieldTotal)
        hMCgenPrimYieldTotal->Write();
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hYieldCorrectedRatioToMC[i])
        hYieldCorrectedRatioToMC[i]->Write();
    }
    
    if (hYieldCorrectedTotalRatioToMC)
        hYieldCorrectedTotalRatioToMC->Write();
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hFractionCorrected[i])
        hFractionCorrected[i]->Write();
      
      if (hFractionCorrectedSysError[i])
        hFractionCorrectedSysError[i]->Write();
    }
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (hMCgenPrimFraction[i])
        hMCgenPrimFraction[i]->Write();
    }
    
    if (cCorrData)
      cCorrData->Write();
    
    if (cCorrYieldsRatio)
      cCorrYieldsRatio->Write();
    
    if (cCorrFractions)
      cCorrFractions->Write();
  }
  
  saveFile->cd();
  
  TNamed* settings = new TNamed(
      Form("Settings: Efficiency file \"%s\", numJetsForEfficiency file \"%s\", Data file \"%s\", lowerCentrality %.3f, upperCentrality %.3f, lowerJetPt %.1f, upperJetPt %.1f, constantCorrectionAbovePtThreshold %.3f\n",
           pathNameEfficiency.Data(), pathNameDataMC.Data(), pathNameData.Data(), lowerCentrality, upperCentrality, lowerJetPt,
           upperJetPt, constantCorrectionAbovePtThreshold), "");
  settings->Write();
  
  saveFile->Close();
  
  
  return 0;
}