#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliPID.h"

#include <iostream>

#include "THnSparseDefinitions.h"

#include "../../UserTasks/AliAnalysisTaskPID.h"
#include "SystematicErrorUtils.h"

enum axisDataProj { kPidPtProj = 0, kPidJetPtProj = 1, kPidZProj = 2, kPidXiProj = 3, kNDimsProj = 4 };



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
void normaliseYieldHist2D(TH2* hData, TH2* hNumJets, const Int_t lowerCentralityBinLimit, const Int_t upperCentralityBinLimit)
{
  // Normalise to 1/numJets and bin width. NOTE: jetPt binning of hData and hNumJets assumed to be the same!
  
  if (!hData)
    return;
  
  for (Int_t binJetPt = 0; binJetPt <= hData->GetNbinsY() + 1; binJetPt++) {
    // No normalisation to numJets, if no histo is provided
    const Double_t numJets = hNumJets ? hNumJets->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, binJetPt, binJetPt) : 1.;
    Bool_t noJets = numJets < 1e-13;
    
    for (Int_t binObs = 0; binObs <= hData->GetNbinsX() + 1; binObs++) {
      if (noJets) {
        if (hData->GetBinContent(binObs, binJetPt) > 0.) {
          printf("Error: No jets for jetPt ~ %f, but found content %f at y-coord %f!\n", hData->GetYaxis()->GetBinCenter(binJetPt),
                 hData->GetBinContent(binObs, binJetPt),  hData->GetXaxis()->GetBinCenter(binObs));
        }
        continue;
      }
      const Double_t dObservable = hData->GetXaxis()->GetBinWidth(binObs);
      const Double_t normFactor = 1. / (numJets * dObservable);
      
      hData->SetBinContent(binObs, binJetPt, hData->GetBinContent(binObs, binJetPt) * normFactor);
      hData->SetBinError(binObs, binJetPt, hData->GetBinError(binObs, binJetPt) * normFactor);
    }
  }
}

//____________________________________________________________________________________________________________________
void normaliseYieldHist(TH1* h, Double_t numJets)
{
  // Normalise to 1/numJets and bin width
  
  if (numJets <= 0) // Do not normalise
    numJets = 1.;
  
  for (Int_t bin = 0; bin <= h->GetNbinsX() + 1; bin++) {
    const Double_t dObservable = h->GetBinWidth(bin);
    const Double_t normFactor = 1. / (numJets * dObservable);
    h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
    h->SetBinError(bin, h->GetBinError(bin) * normFactor);
  }
}


//___________________________________________________________________
void GenerateParticleYields(THnSparse* hDataProj, AliAnalysisTaskPID* pidTask, const Double_t centrality,
                            TH2D** hIDFFtrackPt, TH2D** hIDFFz, TH2D** hIDFFxi,
                            const Bool_t setMean, const Bool_t addErrorsQuadratically, const Bool_t smearByError,
                            const Bool_t uniformSystematicError, const Bool_t takeIntoAccountSysError, Int_t nGenerations)
{
  if (!hDataProj || !pidTask || !hIDFFtrackPt || !hIDFFz || !hIDFFxi) {
    printf("Cannot generate particle fractions - missing input!\n");
    return;
  }
  
  if (takeIntoAccountSysError && smearByError) {
    printf("It doesn't make sense to correlate statistical and systematic errors!\n");
    return;
  }
  
  const Int_t nDimProj = hDataProj->GetNdimensions();
  const Long64_t nBinsTHnSparseProj = hDataProj->GetNbins();
  Double_t binContent = 0., binError2 = 0.;
  Int_t binCoord[nDimProj];
  Double_t binCentre[nDimProj];
  Double_t prob[AliPID::kSPECIES];
  
  Bool_t success = kTRUE;
  
  // NOTE: In the following, the bin error is divided according to the fraction. The idea is to divide bin content (and error)
  // into independent sub-samples according to the fractions. But then, adding up the errors of the sub-samples quadratically
  // should recover the original error (think of 100 tracks, error sqrt(100), divided into 10 samples a 10 tracks => error of
  // samples should be sqrt(10) then).
  // This means, that one needs to add error^2 * fraction and NOT error^2 * fraction^2!!
  // However, the errors are ignored anyway at the moment....
  
  if (!smearByError && !takeIntoAccountSysError) {
    nGenerations = 1;
    // Just calculate fractions/spectra w/o any smearing
    const Int_t smearSpeciesByError = -1;
    const Int_t takeIntoAccountSysErrorOfSpecies = -1;
  
    for (Long64_t bin = 0; bin < nBinsTHnSparseProj; bin++) {
      binContent = hDataProj->GetBinContent(bin, binCoord);
      binError2  = hDataProj->GetBinError2(bin);
      
      // If the bin is empty, do not compute the particle fractions -> This bin contributes nothing anyway
      if (binContent < 2. * std::numeric_limits<double>::min())
        continue;
      
      for (Int_t dim = 0; dim < nDimProj; dim++) 
        binCentre[dim] = hDataProj->GetAxis(dim)->GetBinCenter(binCoord[dim]);
      
      // Since there is no smearing, the fraction will stay the same and there is no need to save the result for one iteration in a
      // histogram (cf. case smearByError >= 0).
      success = pidTask->GetParticleFractions(binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, prob,
                                              smearSpeciesByError, takeIntoAccountSysErrorOfSpecies, uniformSystematicError);
      
      if (success) {
        for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
          // Track pT
          hIDFFtrackPt[species]->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                              hIDFFtrackPt[species]->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj])
                                              + binContent * prob[species]);
          Double_t tempError = hIDFFtrackPt[species]->GetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
          hIDFFtrackPt[species]->SetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                            TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
          
          // z
          hIDFFz[species]->SetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj], 
                                        hIDFFz[species]->GetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj])
                                        + binContent * prob[species]);
          tempError = hIDFFz[species]->GetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj]);
          hIDFFz[species]->SetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj],
                                      TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
          
          // xi
          hIDFFxi[species]->SetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj], 
                                          hIDFFxi[species]->GetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj])
                                          + binContent * prob[species]);
          tempError = hIDFFxi[species]->GetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj]);
          hIDFFxi[species]->SetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj],
                                        TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
        }
      }
      else 
        printf("Problem obtaining fractions for: trackPt %f, jetPt %f, cent %f\n", binCentre[kPidPtProj], binCentre[kPidJetPtProj],
              centrality);
    }
  }
  else {
    // Calculate fractions/spectra w/ smearing "nGenerations" times for every species and obtain the statistical or systematic error 
    // from the spread of the results.
    // Note: For a given species, only the spread of the variation of exactly that species is considered, i.e. one is not
    // looking at the change of the fraction if another species is varied (this would bias the result towards smaller errors
    // by construction since the fraction of the considered species changes weighted with its statistical/systematic error)
    
    TH2D* hIDFFtrackPtVaried[AliPID::kSPECIES][nGenerations];
    TH2D* hIDFFzVaried[AliPID::kSPECIES][nGenerations];
    TH2D* hIDFFxiVaried[AliPID::kSPECIES][nGenerations];
    
    // Create this histo and set all entries to -1 (= not calculated yet) in every loop
    TH2D* hPartFractionsSmeared = hDataProj->Projection(kPidJetPtProj, kPidPtProj);
    hPartFractionsSmeared->SetName("hPartFractionsSmeared");
    
    // Generate results with varying fractions
    for (Int_t smearSpeciesByError = 0; smearSpeciesByError < AliPID::kSPECIES; smearSpeciesByError++) {
      for (Int_t i = 0; i < nGenerations; i++) {
        hIDFFtrackPtVaried[smearSpeciesByError][i] = new TH2D(*(TH2D*)hIDFFtrackPt[smearSpeciesByError]);
        hIDFFtrackPtVaried[smearSpeciesByError][i]->SetName(Form("hIDFFtrackPtVaried_%d_%d", smearSpeciesByError, i));
        hIDFFtrackPtVaried[smearSpeciesByError][i]->Reset();
        
        hIDFFzVaried[smearSpeciesByError][i] = new TH2D(*(TH2D*)hIDFFz[smearSpeciesByError]);
        hIDFFzVaried[smearSpeciesByError][i]->SetName(Form("hIDFFzVaried_%d_%d", smearSpeciesByError, i));
        hIDFFzVaried[smearSpeciesByError][i]->Reset();
        
        hIDFFxiVaried[smearSpeciesByError][i] = new TH2D(*(TH2D*)hIDFFxi[smearSpeciesByError]);
        hIDFFxiVaried[smearSpeciesByError][i]->SetName(Form("hIDFFxiVaried_%d_%d", smearSpeciesByError, i));
        hIDFFxiVaried[smearSpeciesByError][i]->Reset();
        
        // In each iteration (moving over all bins), one want only ONE fixed particle fraction per bin in the spectra map.
        // Otherwise, one would assing different fractions to e.g. the same trackPt bin, if only z and xi is different
        // (but jetPt bin is still the same). This would then cause the fluctuations to cancel partially.
        // Thus: Calculate the smeared fraction for each bin in the spectra map only once and store it in a histo. 
        // If it is requested again (i.e. histoEntry >= 0), use the value from the histo.
        
        // Set all entries to -1 (= not calculated yet)
        for (Int_t iX = 0; iX <= hPartFractionsSmeared->GetNbinsX(); iX++) {
          for (Int_t iY = 0; iY <= hPartFractionsSmeared->GetNbinsY(); iY++) {
            hPartFractionsSmeared->SetBinContent(iX, iY, -1);
          }
        }
        
        for (Long64_t bin = 0; bin < nBinsTHnSparseProj; bin++) {
          binContent = hDataProj->GetBinContent(bin, binCoord);
          binError2  = hDataProj->GetBinError2(bin);
          
          // If the bin is empty, do not compute the particle fractions -> This bin contributes nothing anyway
          if (binContent < 2. * std::numeric_limits<double>::min())
            continue;
          
          for (Int_t iS = 0; iS < AliPID::kSPECIES; iS++)
            prob[iS] = 0.;
          
          for (Int_t dim = 0; dim < nDimProj; dim++) 
            binCentre[dim] = hDataProj->GetAxis(dim)->GetBinCenter(binCoord[dim]);
          
          if (hPartFractionsSmeared->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj]) >= 0) {
            success = kTRUE;
            prob[smearSpeciesByError] = hPartFractionsSmeared->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
          }
          else {
            const Int_t smearSpeciesByStatisticalError = smearByError ? smearSpeciesByError : -1;
            const Int_t takeIntoAccountSysErrorOfSpecies = takeIntoAccountSysError ? smearSpeciesByError : -1;
            
            success = pidTask->GetParticleFractions(binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, prob, 
                                                    smearSpeciesByStatisticalError, takeIntoAccountSysErrorOfSpecies, 
                                                    uniformSystematicError);
            
            if (success)
              hPartFractionsSmeared->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], prob[smearSpeciesByError]);
          }
          
          // NOTE: Since only the bin content of the current species is stored in hPartFractionsSmeared,
          // only prob[smearSpeciesByError] can be used in the following!!!!
          
          if (success) {
            // To make things readable...
            TH2D* hIDFFtrackPtVariedCurr = hIDFFtrackPtVaried[smearSpeciesByError][i];
            TH2D* hIDFFzVariedCurr       = hIDFFzVaried[smearSpeciesByError][i];
            TH2D* hIDFFxiVariedCurr      = hIDFFxiVaried[smearSpeciesByError][i];
            
            // Track pT
            hIDFFtrackPtVariedCurr->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                  hIDFFtrackPtVariedCurr->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj])
                                                  + binContent * prob[smearSpeciesByError]);
            Double_t tempError = hIDFFtrackPtVariedCurr->GetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
            hIDFFtrackPtVariedCurr->SetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                TMath::Sqrt(tempError * tempError + binError2 * prob[smearSpeciesByError]));
            
            // z
            hIDFFzVariedCurr->SetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj], 
                                            hIDFFzVariedCurr->GetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj])
                                            + binContent * prob[smearSpeciesByError]);
            tempError = hIDFFzVariedCurr->GetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj]);
            hIDFFzVariedCurr->SetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj],
                                          TMath::Sqrt(tempError * tempError + binError2 * prob[smearSpeciesByError]));
            
            // xi
            hIDFFxiVariedCurr->SetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj], 
                                              hIDFFxiVariedCurr->GetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj])
                                              + binContent * prob[smearSpeciesByError]);
            tempError = hIDFFxiVariedCurr->GetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj]);
            hIDFFxiVariedCurr->SetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj],
                                            TMath::Sqrt(tempError * tempError + binError2 * prob[smearSpeciesByError]));
          }
          else 
            printf("Problem obtaining fractions for: trackPt %f, jetPt %f, cent %f, smearSpeciesByError %d\n", 
                   binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, smearSpeciesByError);
        }
      }
    }
    
    delete hPartFractionsSmeared;
    
    const Int_t nHistos = nGenerations;
    
    
    /*OLD error for each species by variation of ALL species (will bias towards smaller errors!)
    // TODO Still does not store all the values in histogram to avoid different fractions for same map bin in same iteration.
    // => If this is build in, one needs a histogram for EACH species!!!

    // Calculate fractions/spectra w/ smearing "nGenerations" times for every species and obtain the statistical error from
    // the spread of the results
    
    TH2D* hIDFFtrackPtVaried[AliPID::kSPECIES][AliPID::kSPECIES * nGenerations];
    TH2D* hIDFFzVaried[AliPID::kSPECIES][AliPID::kSPECIES * nGenerations];
    TH2D* hIDFFxiVaried[AliPID::kSPECIES][AliPID::kSPECIES * nGenerations];
    
    // Generate results with varying fractions
    for (Int_t smearSpeciesByError = 0; smearSpeciesByError < AliPID::kSPECIES; smearSpeciesByError++) {
      for (Int_t i = 0; i < nGenerations; i++) {
        for (Int_t consideredSpecies = 0; consideredSpecies < AliPID::kSPECIES; consideredSpecies++) {
          hIDFFtrackPtVaried[consideredSpecies][smearSpeciesByError * nGenerations + i] = 
            new TH2D(*(TH2D*)hIDFFtrackPt[consideredSpecies]);
          hIDFFtrackPtVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->Reset();
          hIDFFtrackPtVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->SetName(Form("hIDFFtrackPtVaried_%d_%d",
                                                                                          consideredSpecies,
                                                                                          smearSpeciesByError * nGenerations + i));
          hIDFFzVaried[consideredSpecies][smearSpeciesByError * nGenerations + i] = 
            new TH2D(*(TH2D*)hIDFFz[consideredSpecies]);
          hIDFFzVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->Reset();
          hIDFFzVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->SetName(Form("hIDFFzVaried_%d_%d",
                                                                                          consideredSpecies,
                                                                                          smearSpeciesByError * nGenerations + i));
          hIDFFxiVaried[consideredSpecies][smearSpeciesByError * nGenerations + i] = 
            new TH2D(*(TH2D*)hIDFFxi[consideredSpecies]);
          hIDFFxiVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->Reset();
          hIDFFxiVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->SetName(Form("hIDFFxiVaried_%d_%d",
                                                                                          consideredSpecies,
                                                                                          smearSpeciesByError * nGenerations + i));
        }
        
        for (Long64_t bin = 0; bin < nBinsTHnSparseProj; bin++) {
          binContent = hDataProj->GetBinContent(bin, binCoord);
          binError2  = hDataProj->GetBinError2(bin);
          
          // If the bin is empty, do not compute the particle fractions -> This bin contributes nothing anyway
          if (binContent < 2. * std::numeric_limits<double>::min())
            continue;
          
          for (Int_t dim = 0; dim < nDimProj; dim++) 
            binCentre[dim] = hDataProj->GetAxis(dim)->GetBinCenter(binCoord[dim]);
          
          success = pidTask->GetParticleFractions(binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, prob, 
                                                  smearSpeciesByError);
          
          if (success) {
            for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
              // To make things readable...
              TH2D* hIDFFtrackPtVariedCurr = hIDFFtrackPtVaried[species][smearSpeciesByError * nGenerations + i];
              TH2D* hIDFFzVariedCurr       = hIDFFzVaried[species][smearSpeciesByError * nGenerations + i];
              TH2D* hIDFFxiVariedCurr      = hIDFFxiVaried[species][smearSpeciesByError * nGenerations + i];
              
              // Track pT
              hIDFFtrackPtVariedCurr->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                    hIDFFtrackPtVariedCurr->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj])
                                                    + binContent * prob[species]);
              Double_t tempError = hIDFFtrackPtVariedCurr->GetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
              hIDFFtrackPtVariedCurr->SetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                  TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
              
              // z
              hIDFFzVariedCurr->SetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj], 
                                              hIDFFzVariedCurr->GetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj])
                                              + binContent * prob[species]);
              tempError = hIDFFzVariedCurr->GetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj]);
              hIDFFzVariedCurr->SetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj],
                                            TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
              
              // xi
              hIDFFxiVariedCurr->SetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj], 
                                               hIDFFxiVariedCurr->GetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj])
                                               + binContent * prob[species]);
              tempError = hIDFFxiVariedCurr->GetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj]);
              hIDFFxiVariedCurr->SetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj],
                                             TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
            }
          }
          else 
            printf("Problem obtaining fractions for: trackPt %f, jetPt %f, cent %f, smearSpeciesByError %d\n", 
                   binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, smearSpeciesByError);
        }
      }
    }
    
    const Int_t nHistos = AliPID::kSPECIES * nGenerations;
    */
    
    // Compare results to obtain error
    const Bool_t ignoreSigmaErrors = kTRUE;
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!extractSystematicError(nHistos, hIDFFtrackPtVaried[species], hIDFFtrackPt[species], setMean, addErrorsQuadratically,
                                  ignoreSigmaErrors))
        printf("Failed to determine systematic error for trackPt, species %d\n", species);
      
      if (!extractSystematicError(nHistos, hIDFFzVaried[species], hIDFFz[species], setMean, addErrorsQuadratically, ignoreSigmaErrors))
        printf("Failed to determine systematic error for z, species %d\n", species);
      
      if (!extractSystematicError(nHistos, hIDFFxiVaried[species], hIDFFxi[species], setMean, addErrorsQuadratically, 
                                  ignoreSigmaErrors))
        printf("Failed to determine systematic error for xi, species %d\n", species);
      
      for (Int_t i = 0; i < nHistos; i++) {
        delete hIDFFtrackPtVaried[species][i];
        hIDFFtrackPtVaried[species][i] = 0x0;
        
        delete hIDFFzVaried[species][i];
        hIDFFzVaried[species][i] = 0x0;
        
        delete hIDFFxiVaried[species][i];
        hIDFFxiVaried[species][i] = 0x0;
      }
    }
  }
}
  
  
//___________________________________________________________________
Int_t extractFFs(TString particleFractionPackagePathName, TString pathNameData, TString listName /*= ""*/,
                 Bool_t uniformSystematicError, Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                 Double_t lowerCentrality = -2, Double_t upperCentrality = -2, Int_t rebinPt = 1, Int_t rebinZ = 1, Int_t rebinXi = 1)
{
  TObjArray* histList = 0x0;
  
  if (listName == "") {
    listName = pathNameData;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
  }
  
  TFile* f = TFile::Open(pathNameData.Data());
  if (!f)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameData.Data() << "\"!" << std::endl;
    return -1;
  }
  
  histList = (TObjArray*)(f->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data histogram
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
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= hPIDdata->GetAxis(kPidCentrality)->GetNbins()) {
      actualLowerCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinUpEdge(upperCentralityBinLimit);

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
  }
  else {
    std::cout << "All" << std::endl;
  }
    
  if (restrictCentralityAxis) {
    hPIDdata->GetAxis(kPidCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  }
  
  // If desired, restrict charge axis
  std::cout << "Charge selection: ";
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
  
  const Int_t indexChargeAxisData = GetAxisByTitle(hPIDdata, "Charge (e_{0})");
  if (indexChargeAxisData < 0 && restrictCharge) {
    std::cout << "Error: Charge axis not found for data histogram!" << std::endl;
    return -1;
  }
  Int_t lowerChargeBinLimitData = -1;
  Int_t upperChargeBinLimitData = -2;
  Double_t actualLowerChargeData = -999;
  Double_t actualUpperChargeData = -999;
  
  if (restrictCharge) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(-1. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(0. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimitData <= upperChargeBinLimitData && lowerChargeBinLimitData >= 1
        && upperChargeBinLimitData <= hPIDdata->GetAxis(indexChargeAxisData)->GetNbins()) {
      actualLowerChargeData = hPIDdata->GetAxis(indexChargeAxisData)->GetBinLowEdge(lowerChargeBinLimitData);
      actualUpperChargeData = hPIDdata->GetAxis(indexChargeAxisData)->GetBinUpEdge(upperChargeBinLimitData);
      
      std::cout << "Charge range data: " << actualLowerChargeData << " - " << actualUpperChargeData << std::endl;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    hPIDdata->GetAxis(indexChargeAxisData)->SetRange(lowerChargeBinLimitData, upperChargeBinLimitData);
  }
  
  // Just take one arbitrary selectSpecies to avoid multiple counting
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1);
  
  // Get projection on dimensions relevant for FFs
  const Int_t nDimProj = kNDimsProj;
  Int_t dimProj[nDimProj];
  dimProj[kPidPtProj] = kPidPt;
  dimProj[kPidJetPtProj] = kPidJetPt;
  dimProj[kPidZProj] = kPidZ;
  dimProj[kPidXiProj] =kPidXi;
  
  THnSparse* hDataProj = hPIDdata->Projection(nDimProj, dimProj, "e");
  
  // If desired, rebin axes
  if (rebinPt != 1 || rebinZ != 1 || rebinXi != 1) {
    Int_t group[hDataProj->GetNdimensions()];
    
    for (Int_t i = 0; i < hDataProj->GetNdimensions(); i++) {
      if (i == kPidPtProj)
        group[i] = rebinPt;
      else if (i == kPidZProj)
        group[i] = rebinZ;
      else if (i == kPidXiProj)
        group[i] = rebinXi;
      else
        group[i] = 1;
    }
    
    THnSparse* hTemp = hDataProj->Rebin(group);
    hDataProj->SetName("temp");
    delete hDataProj;

    hDataProj = hTemp;
  }
  
  /* OLD - Now normalisation to nJets
  Double_t numEvents = -1;
  TH1* hNumEvents = dynamic_cast<TH1*>(histList->FindObject("fhEventsProcessed"));
  if (!hNumEvents) {
    std::cout << std::endl;
    std::cout << "Histo with number of processed events not found! Yields will NOT be normalised to this number!" << std::endl 
              << std::endl;
  }
  else {
    numEvents = hNumEvents->Integral(lowerCentralityBinLimit, upperCentralityBinLimit);
    
    if (numEvents <= 0) {
      numEvents = -1;
      std::cout << std::endl;
      std::cout << "Number of processed events < 1 in selected range! Yields will NOT be normalised to this number!"
                << std::endl << std::endl;
    }
  }*/
  
  
  TH2D* hNjetsGen = 0x0;
  TH2D* hNjetsRec = 0x0;
  
  TH2D* hMCgenPrimYieldPt[AliPID::kSPECIES];
  TH2D* hMCgenPrimYieldZ[AliPID::kSPECIES];
  TH2D* hMCgenPrimYieldXi[AliPID::kSPECIES];
  
  hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
  hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
  
  if (!hNjetsRec) {
    printf("Failed to load number of jets (rec) histo!\n");
    
    // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
    // period is considered; also: No multiplicity information)
    TFile* fBackward = TFile::Open("finalCuts/pp/7TeV/LHC10e.pass2/corrected/finalisedSplines/finalMapsAndTail/Jets/noCutOn_ncl_or_liav/AnalysisResults.root");
    
    TString dirDataInFile = "";
    TDirectory* dirData = fBackward ? (TDirectory*)fBackward->Get(fBackward->GetListOfKeys()->At(0)->GetName()) : 0x0;
  
    TList* list = dirData ? (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName()) : 0x0;

    TH1D* hFFJetPtRec = list ? (TH1D*)list->FindObject("fh1FFJetPtRecCutsInc") : 0x0;
    
    if (hFFJetPtRec) {
      printf("***WARNING: For backward compatibility, using file \"finalCuts/pp/7TeV/LHC10e.pass2/corrected/finalisedSplines/finalMapsAndTail/Jets/noCutOn_ncl_or_liav/AnalysisResults.root\" to get number of jets. BUT: Might be wrong period and has no mult info!***\n");
      printf("ALSO: Using Njets for inclusive jets!!!!\n");
      
      hNjetsRec = new TH2D("fh2FFJetPtRec", "", 1, -1, 1,  hDataProj->GetAxis(kPidJetPtProj)->GetNbins(),
                          hDataProj->GetAxis(kPidJetPtProj)->GetXbins()->GetArray());
      
      for (Int_t iJet = 1; iJet <= hNjetsRec->GetNbinsY(); iJet++) {
        Int_t lowerBin = hFFJetPtRec->FindBin(hNjetsRec->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
        Int_t upperBin = hFFJetPtRec->FindBin(hNjetsRec->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
        hNjetsRec->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
      }
    }
    
    if (!hNjetsRec)
      return -1;
  }
  
  THnSparse* hMCgeneratedYieldsPrimaries = (THnSparse*)histList->FindObject("fhMCgeneratedYieldsPrimaries");
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    hMCgenPrimYieldPt[i] = 0x0;
    hMCgenPrimYieldZ[i] = 0x0;
    hMCgenPrimYieldXi[i] = 0x0;
  }
  
  if (hMCgeneratedYieldsPrimaries && hMCgeneratedYieldsPrimaries->GetEntries() > 0) {
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
    
    // If desired, rebin axes
    if (rebinPt != 1 || rebinZ != 1 || rebinXi != 1) {
      Int_t group[hMCgeneratedYieldsPrimaries->GetNdimensions()];
      
      for (Int_t k = 0; k < hMCgeneratedYieldsPrimaries->GetNdimensions(); k++) {
        if (k == kPidGenYieldPt)
          group[k] = rebinPt;
        else if (k == kPidGenYieldZ)
          group[k] = rebinZ;
        else if (k == kPidGenYieldXi)
          group[k] = rebinXi;
        else
          group[k] = 1;
      }
      
      THnSparse* hTemp = hMCgeneratedYieldsPrimaries->Rebin(group);
      hMCgeneratedYieldsPrimaries->SetName("temp");
      delete hMCgeneratedYieldsPrimaries;
      
      hMCgeneratedYieldsPrimaries = hTemp;
    }


    if (restrictCentralityAxis)
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    
    if (restrictCharge) {
      const Int_t indexChargeAxisGenYield = GetAxisByTitle(hMCgeneratedYieldsPrimaries, "Charge (e_{0})");
      if (indexChargeAxisGenYield < 0) {
        std::cout << "Error: Charge axis not found for gen yield histogram!" << std::endl;
        return -1;
      }
  
      Int_t lowerChargeBinLimitGenYield = -1;
      Int_t upperChargeBinLimitGenYield = -2;
      Double_t actualLowerChargeGenYield = -999;
      Double_t actualUpperChargeGenYield = -999;
  
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
        actualLowerChargeGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetBinLowEdge(lowerChargeBinLimitGenYield);
        actualUpperChargeGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetBinUpEdge(upperChargeBinLimitGenYield);
        
        if (TMath::Abs(actualLowerChargeGenYield - actualLowerChargeData) > 1e-4 ||
            TMath::Abs(actualUpperChargeGenYield - actualUpperChargeData) > 1e-4) {
          std::cout << std::endl;
          std::cout << "Error: Charge range gen yield: " << actualLowerChargeGenYield << " - " << actualUpperChargeGenYield
                    << std::endl << "differs from that of data: " << actualLowerChargeData << " - " << actualUpperChargeData
                    << std::endl;
          return -1;
        }
      }
      else {
        std::cout << std::endl;
        std::cout << "Requested charge range (gen yield) out of limits or upper and lower limit are switched!" << std::endl;
        return -1;
      }
      
      hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->SetRange(lowerChargeBinLimitGenYield, 
                                                                              upperChargeBinLimitGenYield);
    }
  
    for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(MCid + 1, MCid + 1);
      
      hMCgenPrimYieldPt[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldJetPt, kPidGenYieldPt, "e");
      hMCgenPrimYieldPt[MCid]->SetName(Form("hMCgenYieldsPrimPt_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimYieldPt[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      hMCgenPrimYieldPt[MCid]->GetZaxis()->SetTitle("1/N_{Jets} dN/dP_{T} (GeV/c)^{-1}");
      hMCgenPrimYieldPt[MCid]->SetStats(kFALSE);
      
      hMCgenPrimYieldZ[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldJetPt, kPidGenYieldZ, "e");
      hMCgenPrimYieldZ[MCid]->SetName(Form("hMCgenYieldsPrimZ_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimYieldZ[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      hMCgenPrimYieldZ[MCid]->GetZaxis()->SetTitle("1/N_{Jets} dN/dz");
      hMCgenPrimYieldZ[MCid]->SetStats(kFALSE);
      
      hMCgenPrimYieldXi[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldJetPt, kPidGenYieldXi, "e");
      hMCgenPrimYieldXi[MCid]->SetName(Form("hMCgenYieldsPrimXi_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimYieldXi[MCid]->SetTitle(Form("MC truth generated primary yield, %s", AliPID::ParticleName(MCid)));
      hMCgenPrimYieldXi[MCid]->GetZaxis()->SetTitle("1/N_{Jets} dN/d#xi");
      hMCgenPrimYieldXi[MCid]->SetStats(kFALSE);
      
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(0, -1);
    }
  }
  
  
  AliAnalysisTaskPID *pidTask = new AliAnalysisTaskPID("spectrumExtractorTask");
  
  if (!pidTask->SetParticleFractionHistosFromFile(particleFractionPackagePathName, kFALSE)) {
    printf("Failed to load particle fraction package from file \"%s\"!\n", particleFractionPackagePathName.Data());
    return -1;
  }
  
  if (!pidTask->SetParticleFractionHistosFromFile(particleFractionPackagePathName, kTRUE)) {
    printf("Failed to load particle fraction sys error package from file \"%s\"!\n", particleFractionPackagePathName.Data());
    return -1;
  }
  
  printf("Loaded particle fraction package from file \"%s\"!\n", particleFractionPackagePathName.Data());
  
  TH2D* hIDFFtrackPt[AliPID::kSPECIES];
  TH2D* hIDFFz[AliPID::kSPECIES];
  TH2D* hIDFFxi[AliPID::kSPECIES];
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    hIDFFtrackPt[i] = hDataProj->Projection(kPidJetPtProj, kPidPtProj);
    hIDFFtrackPt[i]->Reset();
    if (hIDFFtrackPt[i]->GetSumw2N() <= 0)
      hIDFFtrackPt[i]->Sumw2();
    setupHist(hIDFFtrackPt[i], Form("hIDFFtrackPt_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
              hDataProj->GetAxis(kPidPtProj)->GetTitle(), hDataProj->GetAxis(kPidJetPtProj)->GetTitle(), kBlack);
    
    hIDFFz[i] = hDataProj->Projection(kPidJetPtProj, kPidZProj);
    hIDFFz[i]->Reset();
    if (hIDFFz[i]->GetSumw2N() <= 0)
      hIDFFz[i]->Sumw2();
    setupHist(hIDFFz[i], Form("hIDFFz_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
              hDataProj->GetAxis(kPidZProj)->GetTitle(), hDataProj->GetAxis(kPidJetPtProj)->GetTitle(), kBlack);
    
    hIDFFxi[i] = hDataProj->Projection(kPidJetPtProj, kPidXiProj);
    hIDFFxi[i]->Reset();
    if (hIDFFxi[i]->GetSumw2N() <= 0)
      hIDFFxi[i]->Sumw2();    
    setupHist(hIDFFxi[i], Form("hIDFFxi_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
              hDataProj->GetAxis(kPidXiProj)->GetTitle(), hDataProj->GetAxis(kPidJetPtProj)->GetTitle(), kBlack);
  }
  
  const Double_t centrality = (actualLowerCentrality + actualUpperCentrality) / 2.;
  
  // First iteration: Just take the default fractions to obtain the "mean" of the yields
  Bool_t setMean = kTRUE;
  Bool_t addErrorsQuadratically = kFALSE;
  Bool_t smearByError = kFALSE;
  Bool_t takeIntoAccountSysError = kFALSE;
  // setMean and all the follwoing parameters are anyway irrelevant, if smearByError = kFALSE and takeIntoAccountSysError = kFALSE
  GenerateParticleYields(hDataProj, pidTask, centrality, hIDFFtrackPt, hIDFFz, hIDFFxi, setMean, addErrorsQuadratically,
                         smearByError, uniformSystematicError, takeIntoAccountSysError, 1);
  

  // Next iteration: Vary the PID map within statistical errors several times and calculate the resulting statistical error
  // of the spectra from this. But since the mean of the fraction is some kind of "best estimate" of the truth, leave the means
  // as they have been set during the last step.
  
  
  // NOTE: Do NOT add the statistical errors from the last step, since the yields/fractions (rel. error is the same) already
  // contain the stat. uncertainty of the yield of this species in the corresponding bin! I.e. one just takes some kind of weighted
  // mean of the fractions (mean and error) (equivalent of taking the yields with corresponding error).
  setMean = kFALSE;
  addErrorsQuadratically = kFALSE;
  smearByError = kTRUE;
  takeIntoAccountSysError = kFALSE;
  const Int_t nGenerations = 5000; //TODO set to 5000 in the end, after all testing was successful
  GenerateParticleYields(hDataProj, pidTask, centrality, hIDFFtrackPt, hIDFFz, hIDFFxi, setMean, addErrorsQuadratically,
                         smearByError, uniformSystematicError, takeIntoAccountSysError, nGenerations);
  
  
  // Next iteration: Vary the PID map within systematic errors several times and calculate the resulting systematic error
  // of the spectra from this. But since the mean of the fraction is some kind of "best estimate" of the truth, leave the means
  // as they have been set during the last step.
  setMean = kFALSE;
  addErrorsQuadratically = kFALSE;
  smearByError = kFALSE;
  takeIntoAccountSysError = kTRUE;
  // Clone histograms with final statistical errors -> Only set systematic errors, but leave mean as it is
  TH2D* hIDFFtrackPtSysError[AliPID::kSPECIES];
  TH2D* hIDFFzSysError[AliPID::kSPECIES];
  TH2D* hIDFFxiSysError[AliPID::kSPECIES];
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    hIDFFtrackPtSysError[i] = new TH2D(*hIDFFtrackPt[i]);
    hIDFFtrackPtSysError[i]->SetName(Form("%s_sysError", hIDFFtrackPt[i]->GetName()));
    hIDFFtrackPtSysError[i]->SetFillStyle(0);
    
    hIDFFzSysError[i] = new TH2D(*hIDFFz[i]);
    hIDFFzSysError[i]->SetName(Form("%s_sysError", hIDFFz[i]->GetName()));
    hIDFFzSysError[i]->SetFillStyle(0);
    
    hIDFFxiSysError[i] = new TH2D(*hIDFFxi[i]);
    hIDFFxiSysError[i]->SetName(Form("%s_sysError", hIDFFxi[i]->GetName()));
    hIDFFxiSysError[i]->SetFillStyle(0);
  }
  
  GenerateParticleYields(hDataProj, pidTask, centrality, hIDFFtrackPtSysError, hIDFFzSysError, hIDFFxiSysError, setMean,
                         addErrorsQuadratically, smearByError, uniformSystematicError, takeIntoAccountSysError, nGenerations);
  
  delete pidTask;
  
  
  // Normalise properly
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    normaliseYieldHist2D(hIDFFtrackPt[species], hNjetsRec, lowerCentralityBinLimit, upperCentralityBinLimit);
    normaliseYieldHist2D(hIDFFz[species], hNjetsRec, lowerCentralityBinLimit, upperCentralityBinLimit);
    normaliseYieldHist2D(hIDFFxi[species], hNjetsRec, lowerCentralityBinLimit, upperCentralityBinLimit);
    
    normaliseYieldHist2D(hIDFFtrackPtSysError[species], hNjetsRec, lowerCentralityBinLimit, upperCentralityBinLimit);
    normaliseYieldHist2D(hIDFFzSysError[species], hNjetsRec, lowerCentralityBinLimit, upperCentralityBinLimit);
    normaliseYieldHist2D(hIDFFxiSysError[species], hNjetsRec, lowerCentralityBinLimit, upperCentralityBinLimit);
    
    normaliseYieldHist2D(hMCgenPrimYieldPt[species], hNjetsGen, lowerCentralityBinLimit, upperCentralityBinLimit);
    normaliseYieldHist2D(hMCgenPrimYieldZ[species], hNjetsGen, lowerCentralityBinLimit, upperCentralityBinLimit);
    normaliseYieldHist2D(hMCgenPrimYieldXi[species], hNjetsGen, lowerCentralityBinLimit, upperCentralityBinLimit);
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
  
  saveFileName.Prepend("output_extractedFFs_");
  saveFileName.ReplaceAll(".root", Form("_centrality_%s%s.root",
                                        restrictCentralityAxis ? Form("%.0f_%.0f.root", actualLowerCentrality, actualUpperCentrality)
                                                               : "all",
                                        chargeString.Data()));
  
  TString saveFilePathName = Form("%s/%s", savePath.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  
  if (!saveFile) {
    printf("Failed to save results to file \"%s\"!\n", saveFilePathName.Data());
    return -1;
  }
  
  saveFile->cd();
  
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hIDFFtrackPt[i])
      hIDFFtrackPt[i]->Write();
    
    if (hIDFFz[i])
      hIDFFz[i]->Write();
    
    if (hIDFFxi[i])
      hIDFFxi[i]->Write();
    
    
    if (hIDFFtrackPtSysError[i])
      hIDFFtrackPtSysError[i]->Write();
    
    if (hIDFFzSysError[i])
      hIDFFzSysError[i]->Write();
    
    if (hIDFFxiSysError[i])
      hIDFFxiSysError[i]->Write();
    
    
    if (hMCgenPrimYieldPt[i])
      hMCgenPrimYieldPt[i]->Write();
    
    if (hMCgenPrimYieldZ[i])
      hMCgenPrimYieldZ[i]->Write();
    
    if (hMCgenPrimYieldXi[i])
      hMCgenPrimYieldXi[i]->Write();
  }
  
  if (hNjetsGen)
    hNjetsGen->Write();
  
  if (hNjetsRec)
    hNjetsRec->Write();
  

  TNamed* settings = new TNamed(
      Form("Settings: Fraction package \"%s\", Data file \"%s\", lowerCentrality %.0f, upperCentrality %.0f, chargeMode %d, uniformSystematicError %d, rebinPt %d, rebinZ %d, rebinXi %d\n",
           particleFractionPackagePathName.Data(), pathNameData.Data(), lowerCentrality, upperCentrality, chargeMode, 
           uniformSystematicError, rebinPt, rebinZ, rebinXi), "");
  settings->Write();
  
  saveFile->Close();
  
  
  
  /*
  Int_t t1 = hIDFFtrackPt[AliPID::kPion]->GetYaxis()->FindBin(5.1);
  Int_t t2 = hIDFFtrackPt[AliPID::kPion]->GetYaxis()->FindBin(9.9);
  
  printf("t1 %d, t2 %d\n", t1, t2);
  new TCanvas();
  TH1D* hTemp = hIDFFtrackPt[AliPID::kPion]->ProjectionX("_pfy1", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFtrackPtSysError[AliPID::kPion]->ProjectionX("sys_pfy1", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldPt[AliPID::kPion]) {
    hTemp = hMCgenPrimYieldPt[AliPID::kPion]->ProjectionX("MC_pfy1", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFz[AliPID::kPion]->ProjectionX("_pfy2", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFzSysError[AliPID::kPion]->ProjectionX("sys_pfy2", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldZ[AliPID::kPion]) {
    hTemp = hMCgenPrimYieldZ[AliPID::kPion]->ProjectionX("MC_pfy2", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFxi[AliPID::kPion]->ProjectionX("_pfy3", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFxiSysError[AliPID::kPion]->ProjectionX("sys_pfy3", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldXi[AliPID::kPion]) {
    hTemp = hMCgenPrimYieldXi[AliPID::kPion]->ProjectionX("MC_pfy3", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  
  
  new TCanvas();
  hTemp = hIDFFtrackPt[AliPID::kProton]->ProjectionX("_pfy4", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFtrackPtSysError[AliPID::kProton]->ProjectionX("sys_pfy4", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldPt[AliPID::kProton]) {
    hTemp = hMCgenPrimYieldPt[AliPID::kProton]->ProjectionX("MC_pfy4", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFz[AliPID::kProton]->ProjectionX("_pfy5", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFzSysError[AliPID::kProton]->ProjectionX("sys_pfy5", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldZ[AliPID::kProton]) {
    hTemp = hMCgenPrimYieldZ[AliPID::kProton]->ProjectionX("MC_pfy5", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFxi[AliPID::kProton]->ProjectionX("_pfy", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFxiSysError[AliPID::kProton]->ProjectionX("sys_pfy", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldXi[AliPID::kProton]) {
    hTemp = hMCgenPrimYieldXi[AliPID::kProton]->ProjectionX("MC_pfy", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }*/
  
  return 0;
}