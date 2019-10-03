#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"

#include "AliPID.h"

#include <iostream>

#include "THnSparseDefinitions.h"

enum type { kTrackPt = 0, kZ = 1, kXi = 2, kNtypes = 3 };


//___________________________________________________________________
THnSparse* rebinPtOfTHnSparse(THnSparse* h, Int_t iDimPt)
{
  // Create a new THnSparse with pT binning as in header file
  
  if (!h)
    return 0x0;
  
  TH1D hDummyPt("hDummyPt", "", nPtBins, binsPt);
  TAxis* axisPt = hDummyPt.GetXaxis();
  
  const Int_t nDimRebinned = h->GetNdimensions();
  
  Int_t binsRebinned[nDimRebinned];
  Double_t xminRebinned[nDimRebinned];
  Double_t xmaxRebinned[nDimRebinned];
  
  for (Int_t iDim = 0; iDim < nDimRebinned; iDim++) {
    if (iDim == iDimPt) {
      binsRebinned[iDim] = axisPt->GetNbins();
      xminRebinned[iDim] = axisPt->GetBinLowEdge(1);
      xmaxRebinned[iDim] = axisPt->GetBinUpEdge(axisPt->GetNbins());
    }
    else {
      binsRebinned[iDim] = h->GetAxis(iDim)->GetNbins();
      xminRebinned[iDim] = h->GetAxis(iDim)->GetBinLowEdge(1);
      xmaxRebinned[iDim] = h->GetAxis(iDim)->GetBinUpEdge(h->GetAxis(iDim)->GetNbins());
    }
  }
  
  THnSparse* hRebinned = new THnSparseD(Form("%s_rebinned", h->GetName()), h->GetTitle(),
                                        nDimRebinned, binsRebinned, xminRebinned, xmaxRebinned);
  hRebinned->Sumw2();
  
  for (Int_t iDim = 0; iDim < nDimRebinned; iDim++) {
    if (iDim == iDimPt) {
      hRebinned->SetBinEdges(iDim, axisPt->GetXbins()->fArray);
    }
    else {
      if (h->GetAxis(iDim)->GetXbins()->fN != 0)
        hRebinned->SetBinEdges(iDim, h->GetAxis(iDim)->GetXbins()->fArray);
    }
    
    hRebinned->GetAxis(iDim)->SetTitle(h->GetAxis(iDim)->GetTitle());
  }
  
  // Now just fill the THnSparse with the original data. RebinnedAdd already takes into account the different binning properly
  // (was also tested explicitely!).
  hRebinned->RebinnedAdd(h, 1.);
  
  return hRebinned;
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
Int_t extractGenMCtruth(TString pathNameData, TString listName,
                        Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                        Double_t lowerCentrality /*= -2*/, Double_t upperCentrality /*= -2*/,
                        Double_t lowerJetPt /*= -1*/ , Double_t upperJetPt/* = -1*/)
{
  if (listName == "") {
    listName = pathNameData;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
  }
  
  TString titles[kNtypes];
  titles[kTrackPt] = "trackPt";
  titles[kZ] = "z";
  titles[kXi] = "xi";
  
  TString pathData = pathNameData;
  pathData.Replace(pathData.Last('/'), pathData.Length(), "");
  
  
  TH2D* hNjetsGenData = 0x0;
  TH2D* hNjetsRecData = 0x0;
  
  TH1D* hMCgenPrimYieldPt[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCgenPrimYieldZ[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCgenPrimYieldXi[AliPID::kSPECIES] = {0x0, };

  TFile* fileData = TFile::Open(pathNameData.Data());
  if (!fileData) {
    printf("Failed to open data file \"%s\"\n", pathNameData.Data());
    return -1;
  }
  
  TObjArray* histList = (TObjArray*)(fileData->Get(listName.Data()));
  THnSparse* hMCgeneratedYieldsPrimaries = (THnSparse*)histList->FindObject("fhMCgeneratedYieldsPrimaries");
  
  if (!hMCgeneratedYieldsPrimaries) {
    printf("Failed to load generated primary yields!\n");
    return -1;
  }
  
  // Set proper errors, if not yet calculated
  if (!hMCgeneratedYieldsPrimaries->GetCalculateErrors()) {
    std::cout << "Re-calculating errors of " << hMCgeneratedYieldsPrimaries->GetName() << "..." << std::endl;
    
    hMCgeneratedYieldsPrimaries->Sumw2();
    
    Long64_t nBinsTHnSparseGenYield = hMCgeneratedYieldsPrimaries->GetNbins();
    Double_t binContentGenYield = 0;
    for (Long64_t bin = 0; bin < nBinsTHnSparseGenYield; bin++) {
      binContentGenYield = hMCgeneratedYieldsPrimaries->GetBinContent(bin);
      hMCgeneratedYieldsPrimaries->SetBinError(bin, TMath::Sqrt(binContentGenYield));
    }
  }
  
  // Set the pointer to new THnSparse with rebinned pT (keep the name)
  hMCgeneratedYieldsPrimaries = rebinPtOfTHnSparse(hMCgeneratedYieldsPrimaries, kPidGenYieldPt);
  
  hNjetsGenData = (TH2D*)histList->FindObject("fh2FFJetPtGen");
  
  hNjetsRecData = (TH2D*)histList->FindObject("fh2FFJetPtRec");
  
  Bool_t restrictJetPtAxis = (lowerJetPt >= 0 && upperJetPt >= 0);
  
  if (restrictJetPtAxis && (!hNjetsRecData || !hNjetsGenData)) {
    printf("Failed to load numJet histo for data!\n");
    
    // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
    // period is considered; also: No multiplicity information)
    
    TString pathBackward = Form("%s/AnalysisResults.root", pathData.Data());
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
      
      hNjetsRecData = new TH2D("fh2FFJetPtRec", "", 1, -1, 1, hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetNbins(),
                               hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetXbins()->GetArray());
      
      for (Int_t iJet = 1; iJet <= hNjetsRecData->GetNbinsY(); iJet++) {
        Int_t lowerBin = hFFJetPtRec->FindBin(hNjetsRecData->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
        Int_t upperBin = hFFJetPtRec->FindBin(hNjetsRecData->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
        hNjetsRecData->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
      }
      
      hNjetsGenData = new TH2D("fh2FFJetPtGen", "", 1, -1, 1,  hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetNbins(),
                               hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetXbins()->GetArray());
      
      for (Int_t iJet = 1; iJet <= hNjetsGenData->GetNbinsY(); iJet++) {
        Int_t lowerBin = hFFJetPtGen->FindBin(hNjetsGenData->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
        Int_t upperBin = hFFJetPtGen->FindBin(hNjetsGenData->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
        hNjetsGenData->SetBinContent(1, iJet, hFFJetPtGen->Integral(lowerBin, upperBin));
      }
    }
    
    if (!hNjetsRecData || ! hNjetsGenData)
      return -1;
  }
  
  const Bool_t restrictCentralityData = ((lowerCentrality >= -1) && (upperCentrality >= -1));
  // Integral(lowerCentBinLimitData, uppCentBinLimitData) will not be restricted if these values are kept
  const Int_t lowerCentralityBinLimitData = restrictCentralityData ? hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->FindBin(lowerCentrality + 0.001) 
                                                                   : -1;
  const Int_t upperCentralityBinLimitData = restrictCentralityData ? hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->FindBin(upperCentrality - 0.001) 
                                                                   : -2;
  
  const Double_t actualLowerCentralityData = hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->GetBinLowEdge(lowerCentralityBinLimitData);
  
  const Double_t actualUpperCentralityData = hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->GetBinLowEdge(upperCentralityBinLimitData);
  
  // Integral(lowerJetBinLimitData, uppJetBinLimitData) will not be restricted if these values are kept
  const Int_t lowerJetPtBinLimitData = restrictJetPtAxis 
                                       ? hNjetsRecData->GetYaxis()->FindBin(lowerJetPt + 0.001) : -1;
  const Int_t upperJetPtBinLimitData  = restrictJetPtAxis 
                                       ? hNjetsRecData->GetYaxis()->FindBin(upperJetPt - 0.001) : -2;
  
  const Double_t actualLowerJetPtData = restrictJetPtAxis ? hNjetsRecData->GetYaxis()->GetBinLowEdge(lowerJetPtBinLimitData) : -1;
  const Double_t actualUpperJetPtData = restrictJetPtAxis ? hNjetsRecData->GetYaxis()->GetBinUpEdge(upperJetPtBinLimitData)  : -1;
  
  if (restrictCentralityData)
    hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->SetRange(lowerCentralityBinLimitData, upperCentralityBinLimitData);
  
  const Bool_t restrictCharge = (chargeMode != kAllCharged);
  
  Int_t lowerChargeBinLimitGenYield = -1;
  Int_t upperChargeBinLimitGenYield = -2;
    
  if (restrictCharge) {
    const Int_t indexChargeAxisGenYield = GetAxisByTitle(hMCgeneratedYieldsPrimaries, "Charge (e_{0})");
    if (indexChargeAxisGenYield < 0) {
      std::cout << "Error: Charge axis not found for gen yield histogram!" << std::endl;
      return -1;
    }

    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(-1. + 0.001);
      upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(0. + 0.001);
      upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimitGenYield <= upperChargeBinLimitGenYield && lowerChargeBinLimitGenYield >= 1
        && upperChargeBinLimitGenYield <= hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetNbins()) {
      // OK
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range (gen yield) out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->SetRange(lowerChargeBinLimitGenYield, upperChargeBinLimitGenYield);
  }
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -1;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (restrictJetPtAxis) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->FindBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->FindBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 &&
        upperJetPtBinLimit <= hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetNbins()) {
      actualLowerJetPt = hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->GetBinUpEdge(upperJetPtBinLimit);

      restrictJetPtAxis = kTRUE;
    
      if (TMath::Abs(actualLowerJetPt - actualLowerJetPtData) > 1e-3 || TMath::Abs(actualUpperJetPt - actualUpperJetPtData) > 1e-3) {
        std::cout << std::endl;
        std::cout << "Mismatch between jet pT range of gen histo and num jet histo!" << std::endl;
        return -1;
      }
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
    hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  
  
  for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
    hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(MCid + 1, MCid + 1);
    
    hMCgenPrimYieldPt[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldPt, "e");
    hMCgenPrimYieldPt[MCid]->SetName(Form("hMCgenYieldsPrimPt_%s", AliPID::ParticleShortName(MCid)));
    hMCgenPrimYieldPt[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
    hMCgenPrimYieldPt[MCid]->GetYaxis()->SetTitle("1/N_{Jets} dN/dp_{T} (GeV/c)^{-1}");
    hMCgenPrimYieldPt[MCid]->SetStats(kFALSE);
    
    if (hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldZ)) { // For backward compatibility
      hMCgenPrimYieldZ[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldZ, "e");
      hMCgenPrimYieldZ[MCid]->SetName(Form("hMCgenYieldsPrimZ_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimYieldZ[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      hMCgenPrimYieldZ[MCid]->GetYaxis()->SetTitle("1/N_{Jets} dN/dz");
      hMCgenPrimYieldZ[MCid]->SetStats(kFALSE);
    }

    if (hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldXi)) {
      hMCgenPrimYieldXi[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldXi, "e");
      hMCgenPrimYieldXi[MCid]->SetName(Form("hMCgenYieldsPrimXi_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimYieldXi[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      hMCgenPrimYieldXi[MCid]->GetYaxis()->SetTitle("1/N_{Jets} dN/d#xi");
      hMCgenPrimYieldXi[MCid]->SetStats(kFALSE);
    }
    
    hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(0, -1);
  }
  
  // Obtain uncorrected fractions
  THnSparse* hPIDdata = dynamic_cast<THnSparse*>(histList->FindObject("hPIDdataAll"));
  if (!hPIDdata) {
    std::cout << std::endl;
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  
  // Set proper errors, if not yet calculated
  if (!hPIDdata->GetCalculateErrors()) {
    std::cout << "Re-calculating errors of " << hPIDdata->GetName() << "..." << std::endl;
    hPIDdata->Sumw2();
    Long64_t nBinsTHnSparse = hPIDdata->GetNbins();
    Double_t binContent = 0;
    
    for (Long64_t bin = 0; bin < nBinsTHnSparse; bin++) {
      binContent = hPIDdata->GetBinContent(bin);
      hPIDdata->SetBinError(bin, TMath::Sqrt(binContent));
    }
  }
  
  // Set the pointer to new THnSparse with rebinned pT (keep the name)
  hPIDdata = rebinPtOfTHnSparse(hPIDdata, kPidPt);
  
  // Bin limits are the same as in gen yield histo
  if (restrictCentralityData) 
    hPIDdata->GetAxis(kPidCentrality)->SetRange(lowerCentralityBinLimitData, upperCentralityBinLimitData);
  
  if (restrictJetPtAxis)
    hPIDdata->GetAxis(kPidJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
  
  const Int_t indexChargeAxisData = GetAxisByTitle(hPIDdata, "Charge (e_{0})");
  if (indexChargeAxisData < 0 && restrictCharge) {
    std::cout << "Error: Charge axis not found for data histogram!" << std::endl;
    return -1;
  }
  
  if (restrictCharge)
    hPIDdata->GetAxis(indexChargeAxisData)->SetRange(lowerChargeBinLimitGenYield, upperChargeBinLimitGenYield);
  
  
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1); // Do not count each particle more than once
  
  
  TH1D* hMCuncorrYieldPt[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCuncorrYieldZ[AliPID::kSPECIES] = {0x0, };
  TH1D* hMCuncorrYieldXi[AliPID::kSPECIES] = {0x0, };
  
  for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
    // Unfortunately, different MC bins than AliPID
    Int_t MCbinInHisto = 0;
    
    switch (MCid) {
      case AliPID::kElectron:
        MCbinInHisto = 1;
        break;
      case AliPID::kKaon:
        MCbinInHisto = 2;
        break;
      case AliPID::kMuon:
        MCbinInHisto = 3;
        break;
      case AliPID::kPion:
        MCbinInHisto = 4;
        break;
      case AliPID::kProton:
        MCbinInHisto = 5;
        break;
      default:
        break;
    }
    
    hPIDdata->GetAxis(kPidMCpid)->SetRange(MCbinInHisto, MCbinInHisto);
    
    hMCuncorrYieldPt[MCid] = hPIDdata->Projection(kPidPt, "e");
    hMCuncorrYieldPt[MCid]->SetName(Form("hMCuncorrYieldsPrimPt_%s", AliPID::ParticleShortName(MCid)));
    hMCuncorrYieldPt[MCid]->SetTitle(Form("MC truth uncorrected yield, %s", AliPID::ParticleName(MCid)));
    hMCuncorrYieldPt[MCid]->GetYaxis()->SetTitle("1/N_{Jets} dN/dp_{T} (GeV/c)^{-1}");
    hMCuncorrYieldPt[MCid]->SetStats(kFALSE);
    
    if (hPIDdata->GetAxis(kPidZ)) { // For backward compatibility
      hMCuncorrYieldZ[MCid] = hPIDdata->Projection(kPidZ, "e");
      hMCuncorrYieldZ[MCid]->SetName(Form("hMCuncorrYieldsPrimZ_%s", AliPID::ParticleShortName(MCid)));
      hMCuncorrYieldZ[MCid]->SetTitle(Form("MC truth uncorrected yield, %s", AliPID::ParticleName(MCid)));
      hMCuncorrYieldZ[MCid]->GetYaxis()->SetTitle("1/N_{Jets} dN/dz");
      hMCuncorrYieldZ[MCid]->SetStats(kFALSE);
    }

    if (hPIDdata->GetAxis(kPidXi)) {
      hMCuncorrYieldXi[MCid] = hPIDdata->Projection(kPidXi, "e");
      hMCuncorrYieldXi[MCid]->SetName(Form("hMCuncorrYieldsPrimXi_%s", AliPID::ParticleShortName(MCid)));
      hMCuncorrYieldXi[MCid]->SetTitle(Form("MC truth uncorrected yield, %s", AliPID::ParticleName(MCid)));
      hMCuncorrYieldXi[MCid]->GetYaxis()->SetTitle("1/N_{Jets} dN/d#xi");
      hMCuncorrYieldXi[MCid]->SetStats(kFALSE);
    }
    
    hPIDdata->GetAxis(kPidMCpid)->SetRange(0, -1);
  }
  
  // Save results to file
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  TString saveFileName = pathNameData;
  saveFileName.Replace(0, pathNameData.Last('/') + 1, "");
  
  TString savePath = pathNameData;
  savePath.ReplaceAll(Form("/%s", saveFileName.Data()), "");
  
  saveFileName.Prepend("output_extractedMC_");
  TString centralityString = restrictCentralityData ? Form("_centrality_%.0f_%.0f", actualLowerCentralityData,
                                                           actualUpperCentralityData)
                                                    : "_centrality_all";
  TString jetPtString = restrictJetPtAxis ? Form("_jetPt_%.0f_%.0f", actualLowerJetPt, actualUpperJetPt)
                                          : "";  
  saveFileName.ReplaceAll(".root", Form("%s%s%s.root", centralityString.Data(), jetPtString.Data(), chargeString.Data()));
  
  TString saveFilePathName = Form("%s/%s", savePath.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  
  if (!saveFile) {
    printf("Failed to save results to file \"%s\"!\n", saveFilePathName.Data());
    return -1;
  }
  
  saveFile->cd();
  
  
  
  // Calculate the fractions and do the normalisation
  
  const Double_t nJetsGenData = hNjetsGenData ? hNjetsGenData->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData, 
                                                                        lowerJetPtBinLimitData, upperJetPtBinLimitData) : 1.;
  const Double_t nJetsRecData = hNjetsRecData ? hNjetsRecData->Integral(lowerCentralityBinLimitData, upperCentralityBinLimitData,                 
                                                                        lowerJetPtBinLimitData, upperJetPtBinLimitData) : 1.;
  
  for (Int_t type = 0; type < kNtypes; type++) {
    TH1D** hMCgenPrimYield;
    TH1D** hMCuncorrYield;
    
    if (type == kTrackPt) {
      hMCgenPrimYield = &hMCgenPrimYieldPt[0];
      hMCuncorrYield = &hMCuncorrYieldPt[0];
    }
    else if (type == kZ) {
      hMCgenPrimYield = &hMCgenPrimYieldZ[0];
      hMCuncorrYield = &hMCuncorrYieldZ[0];
    }
    else if (type == kXi) {
      hMCgenPrimYield = &hMCgenPrimYieldXi[0];
      hMCuncorrYield = &hMCuncorrYieldXi[0];
    }
    else
      continue;
    
    // Gen yields
    if (hMCgenPrimYield[0]) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        normaliseHist(hMCgenPrimYield[species], (nJetsGenData > 0) ? 1. / nJetsGenData : 0.);
        hMCgenPrimYield[species]->SetStats(kFALSE);

        hMCgenPrimYield[species]->SetTitle(Form("%s (MC)", AliPID::ParticleLatexName(species)));
        hMCgenPrimYield[species]->SetMarkerStyle(24);
        hMCgenPrimYield[species]->SetLineColor(getLineColorAliPID(species));
        hMCgenPrimYield[species]->SetMarkerColor(getLineColorAliPID(species));
      }
      
      TH1D* hMCgenPrimYieldTotal = 0x0;
      TH1D* hMCgenPrimFraction[AliPID::kSPECIES];
      for (Int_t i = 0; i < AliPID::kSPECIES; i++)
        hMCgenPrimFraction[i] = 0x0;
      
      hMCgenPrimYieldTotal = new TH1D(*hMCgenPrimYield[0]);
      hMCgenPrimYieldTotal->SetLineColor(kBlack);
      hMCgenPrimYieldTotal->SetMarkerColor(kBlack);
      hMCgenPrimYieldTotal->SetMarkerStyle(24);
      hMCgenPrimYieldTotal->SetName(Form("hMCgenPrimYieldTotal%s", titles[type].Data()));
      hMCgenPrimYieldTotal->SetTitle("Total (MC)");
      //hMCgenPrimYieldTotal->SetTitle("Total generated primary yield (MC truth)");
      
      for (Int_t i = 1; i < AliPID::kSPECIES; i++)
        hMCgenPrimYieldTotal->Add(hMCgenPrimYield[i], 1.);
      
      // Calculate the MC fractions
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hMCgenPrimYield[species]->SetLineColor(getLineColorAliPID(species));
        hMCgenPrimYield[species]->SetMarkerColor(getLineColorAliPID(species));
        hMCgenPrimYield[species]->SetMarkerStyle(24);
      
        hMCgenPrimYield[species]->SetTitle(Form("%s (MC)", AliPID::ParticleLatexName(species)));
        
        hMCgenPrimFraction[species] = new TH1D(*hMCgenPrimYield[species]);
        TString oldName = hMCgenPrimYield[species]->GetName();
        TString newName = oldName.ReplaceAll("Yield", "Fraction");
        hMCgenPrimFraction[species]->SetName(newName.Data());
        hMCgenPrimFraction[species]->GetYaxis()->SetTitle("Particle Fraction");

        // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
        // (numerator is a "selected" subset of the denominator).
        hMCgenPrimFraction[species]->Divide(hMCgenPrimFraction[species], hMCgenPrimYieldTotal, 1., 1., "B"); 
        hMCgenPrimFraction[species]->GetYaxis()->SetRangeUser(0., 1.);
      }
      
      if (hMCgenPrimYieldTotal) {
        hMCgenPrimYieldTotal->GetYaxis()->SetTitleOffset(1.4);
        hMCgenPrimYieldTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
        hMCgenPrimYieldTotal->GetXaxis()->SetNoExponent(kTRUE);
      }
      
      TCanvas* cGenYield = new TCanvas(Form("cGenYield_%s", titles[type].Data()), "Generated Yield", 0, 300, 900, 900);
      if (type == kTrackPt)
        cGenYield->SetLogx(1);
      
      hMCgenPrimYieldTotal->Draw("E1");
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) 
        hMCgenPrimYield[i]->Draw("E1 same");
      
      TLegend* legTemp = cGenYield->BuildLegend();//0.25, 0.16, 0.65, 0.51);
      legTemp->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cGenYield);
      

      TCanvas* cGenFractions = new TCanvas(Form("cGenFractions_%s", titles[type].Data()), "Generated Fractions", 0, 300, 900, 900);
      if (type == kTrackPt)
        cGenFractions->SetLogx(1);
      
      hMCgenPrimFraction[0]->Draw("E1");
      
      for (Int_t i = 1; i < AliPID::kSPECIES; i++)
        hMCgenPrimFraction[i]->Draw("E1 same");
      
      cGenFractions->BuildLegend()->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cGenFractions);
      

      // Calculate the generated MC to-pi-ratios
      TH1D* hMCgenPrimRatioToPi[AliPID::kSPECIES];
      for (Int_t i = 0; i < AliPID::kSPECIES; i++)
        hMCgenPrimRatioToPi[i] = 0x0;
      
      
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (species == AliPID::kPion)
          continue;
        
        hMCgenPrimRatioToPi[species] = new TH1D(*hMCgenPrimYield[species]);
        TString oldName = hMCgenPrimYield[species]->GetName();
        TString newName = oldName.ReplaceAll("Yield", "RatioToPi");
        hMCgenPrimRatioToPi[species]->SetName(newName.Data());
        hMCgenPrimRatioToPi[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
        hMCgenPrimRatioToPi[species]->GetYaxis()->SetTitle("Ratio");

        // Statistically independent data sets -> Just divide
        hMCgenPrimRatioToPi[species]->Divide(hMCgenPrimYield[species], hMCgenPrimYield[AliPID::kPion]); 
        hMCgenPrimRatioToPi[species]->GetYaxis()->SetRangeUser(0., 1.);
      }
      
      TCanvas* cGenRatioToPi = new TCanvas(Form("cGenRatioToPi_%s", titles[type].Data()),
                                                "Generated To-Pi-Ratio",
                                                0, 300, 900, 900);
      if (type == kTrackPt)
        cGenRatioToPi->SetLogx(1);
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (i == AliPID::kPion)
          continue;
        hMCgenPrimRatioToPi[i]->Draw(i == 0 ? "E1" :"E1 same");
      }
      
      legTemp = cGenRatioToPi->BuildLegend();//0.25, 0.16, 0.65, 0.51);
      legTemp->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cGenRatioToPi);
      
      
      // Save results
      saveFile->cd();
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (hMCgenPrimYield[i])
          hMCgenPrimYield[i]->Write();
      }
      
      if (hMCgenPrimYieldTotal)
          hMCgenPrimYieldTotal->Write();
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (hMCgenPrimFraction[i])
          hMCgenPrimFraction[i]->Write();
      }
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (hMCgenPrimRatioToPi[i])
          hMCgenPrimRatioToPi[i]->Write();
      }
      
      if (cGenYield)
        cGenYield->Write();
      
      if (cGenFractions)
        cGenFractions->Write();
      
      if (cGenRatioToPi)
        cGenRatioToPi->Write();
    }
    
    
    // Uncorr rec MC yields
    if (hMCuncorrYield[0]) {
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        // NOTE: These are RECONSTRUCTED yields (just used MC label for PID) -> Must be normalised to number of reconstructed jets!
        normaliseHist(hMCuncorrYield[species], (nJetsRecData > 0) ? 1. / nJetsRecData : 0.);
        hMCuncorrYield[species]->SetStats(kFALSE);

        hMCuncorrYield[species]->SetTitle(Form("%s (MC)", AliPID::ParticleLatexName(species)));
        hMCuncorrYield[species]->SetMarkerStyle(24);
        hMCuncorrYield[species]->SetLineColor(getLineColorAliPID(species));
        hMCuncorrYield[species]->SetMarkerColor(getLineColorAliPID(species));
      }
      
      TH1D* hMCuncorrYieldTotal = 0x0;
      TH1D* hMCuncorrFraction[AliPID::kSPECIES];
      for (Int_t i = 0; i < AliPID::kSPECIES; i++)
        hMCuncorrFraction[i] = 0x0;
      
      hMCuncorrYieldTotal = new TH1D(*hMCuncorrYield[0]);
      hMCuncorrYieldTotal->SetLineColor(kBlack);
      hMCuncorrYieldTotal->SetMarkerColor(kBlack);
      hMCuncorrYieldTotal->SetMarkerStyle(24);
      hMCuncorrYieldTotal->SetName(Form("hMCuncorrYieldTotal%s", titles[type].Data()));
      hMCuncorrYieldTotal->SetTitle("Total (MC)");
      //hMCuncorrYieldTotal->SetTitle("Total uncorrected yield (MC truth)");
      
      for (Int_t i = 1; i < AliPID::kSPECIES; i++)
        hMCuncorrYieldTotal->Add(hMCuncorrYield[i], 1.);
      
      // Calculate the MC fractions
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        hMCuncorrYield[species]->SetLineColor(getLineColorAliPID(species));
        hMCuncorrYield[species]->SetMarkerColor(getLineColorAliPID(species));
        hMCuncorrYield[species]->SetMarkerStyle(24);
      
        hMCuncorrYield[species]->SetTitle(Form("%s (MC)", AliPID::ParticleLatexName(species)));
        
        hMCuncorrFraction[species] = new TH1D(*hMCuncorrYield[species]);
        TString oldName = hMCuncorrYield[species]->GetName();
        TString newName = oldName.ReplaceAll("Yield", "Fraction");
        hMCuncorrFraction[species]->SetName(newName.Data());
        hMCuncorrFraction[species]->GetYaxis()->SetTitle("Particle Fraction");

        // Binomial error as for efficiencies (numerator and denominator are statistically not independent) for correct error calculation
        // (numerator is a "selected" subset of the denominator).
        hMCuncorrFraction[species]->Divide(hMCuncorrFraction[species], hMCuncorrYieldTotal, 1., 1., "B"); 
        hMCuncorrFraction[species]->GetYaxis()->SetRangeUser(0., 1.);
      }
      
      hMCuncorrYieldTotal->GetYaxis()->SetTitleOffset(1.4);
      hMCuncorrYieldTotal->GetXaxis()->SetMoreLogLabels(kTRUE);
      hMCuncorrYieldTotal->GetXaxis()->SetNoExponent(kTRUE);
      
      TCanvas* cMCRecUncorrYield = new TCanvas(Form("cMCRecUncorrYield_%s", titles[type].Data()), "MC Reconstructed Uncorrected Yield",
                                               0, 300, 900, 900);
      if (type == kTrackPt)
        cMCRecUncorrYield->SetLogx(1);
      
      hMCuncorrYieldTotal->Draw("E1");
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) 
        hMCuncorrYield[i]->Draw("E1 same");
      
      TLegend* legTemp = cMCRecUncorrYield->BuildLegend();//0.25, 0.16, 0.65, 0.51);
      legTemp->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cMCRecUncorrYield);
      

      TCanvas* cMCRecUncorrFractions = new TCanvas(Form("cMCRecUncorrFractions_%s", titles[type].Data()),
                                                   "MC Reconstructed Uncorrected Fractions",
                                                   0, 300, 900, 900);
      if (type == kTrackPt)
        cMCRecUncorrFractions->SetLogx(1);
      
      hMCuncorrFraction[0]->Draw("E1");
      
      for (Int_t i = 1; i < AliPID::kSPECIES; i++)
        hMCuncorrFraction[i]->Draw("E1 same");
      
      cMCRecUncorrFractions->BuildLegend()->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cMCRecUncorrFractions);
      
      // Calculate the MC to-pi-ratios
      TH1D* hMCuncorrRatioToPi[AliPID::kSPECIES];
      for (Int_t i = 0; i < AliPID::kSPECIES; i++)
        hMCuncorrRatioToPi[i] = 0x0;
      
      
      for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
        if (species == AliPID::kPion)
          continue;
        
        hMCuncorrRatioToPi[species] = new TH1D(*hMCuncorrYield[species]);
        TString oldName = hMCuncorrYield[species]->GetName();
        TString newName = oldName.ReplaceAll("Yield", "RatioToPi");
        hMCuncorrRatioToPi[species]->SetName(newName.Data());
        hMCuncorrRatioToPi[species]->SetTitle(Form("%s/#pi", AliPID::ParticleLatexName(species)));
        hMCuncorrRatioToPi[species]->GetYaxis()->SetTitle("Ratio");

        // Statistically independent data sets -> Just divide
        hMCuncorrRatioToPi[species]->Divide(hMCuncorrYield[species], hMCuncorrYield[AliPID::kPion]); 
        hMCuncorrRatioToPi[species]->GetYaxis()->SetRangeUser(0., 1.);
      }
      
      TCanvas* cMCRecUncorrRatioToPi = new TCanvas(Form("cMCRecUncorrRatioToPi_%s", titles[type].Data()),
                                                   "MC Reconstructed Uncorrected To-Pi-Ratio",
                                                   0, 300, 900, 900);
      if (type == kTrackPt)
        cMCRecUncorrRatioToPi->SetLogx(1);
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (i == AliPID::kPion)
          continue;
        hMCuncorrRatioToPi[i]->Draw(i == 0 ? "E1" :"E1 same");
      }
      
      legTemp = cMCRecUncorrRatioToPi->BuildLegend();//0.25, 0.16, 0.65, 0.51);
      legTemp->SetFillColor(kWhite);
      
      ClearTitleFromHistoInCanvas(cMCRecUncorrRatioToPi);
      
      // Save results
      saveFile->cd();
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (hMCuncorrYield[i])
          hMCuncorrYield[i]->Write();
      }
      
      if (hMCuncorrYieldTotal)
          hMCuncorrYieldTotal->Write();
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (hMCuncorrFraction[i])
          hMCuncorrFraction[i]->Write();
      }
      
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (hMCuncorrRatioToPi[i])
          hMCuncorrRatioToPi[i]->Write();
      }
      
      if (cMCRecUncorrYield)
        cMCRecUncorrYield->Write();
      
      if (cMCRecUncorrFractions)
        cMCRecUncorrFractions->Write();
      
      if (cMCRecUncorrRatioToPi)
        cMCRecUncorrRatioToPi->Write();
    }
  }
  
  saveFile->cd();
  
  TNamed* settings = new TNamed(
      Form("Settings: Data file \"%s\", lowerCentrality %.3f, upperCentrality %.3f, lowerJetPt %.1f, upperJetPt %.1f\n",
           pathNameData.Data(), lowerCentrality, upperCentrality, lowerJetPt, upperJetPt), "");
  settings->Write();
  
  saveFile->Close();
  
  
  return 0;
}