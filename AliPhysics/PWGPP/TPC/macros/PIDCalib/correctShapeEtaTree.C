#include "TTree.h"
#include "TH2D.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TDatime.h"
#include "TSpline.h"
#include "TF1.h"

#include "AliPID.h"

#include <iostream>

#include "ProgressBar.h"

Double_t getEtaCorrection(TH2D* hMap, Double_t tanTheta, Double_t dEdx)
{
  if (!hMap) 
    return 1.;
    
  if (dEdx < 1.)
    return 1.;
  
  
  //TODO WARNING: Here we can use tanTheta directly but in data we have to extract tanTheta from thetaGlobal!!!
    
  /*
  // For ESD tracks, the local tanTheta could be used (esdTrack->GetInnerParam()->GetTgl()).
  // However, this value is not available for AODs and, thus, not for AliVTrack.
  // Fortunately, the following formula allows to approximate the local tanTheta with the 
  // global theta angle -> This is for by far most of the tracks the same, but gives at
  // maybe the percent level differences within +- 0.2 in tanTheta -> Which is still ok.
  Double_t tanThetaLocal = TMath::Tan(-thetaGlobal + TMath::Pi() / 2.0);
  */
  Int_t binX = hMap->GetXaxis()->FindBin(tanTheta);
  Int_t binY = hMap->GetYaxis()->FindBin(1. / dEdx);
  
  if (binX == 0) 
    binX = 1;
  if (binX > hMap->GetXaxis()->GetNbins())
    binX = hMap->GetXaxis()->GetNbins();
  
  if (binY == 0)
    binY = 1;
  if (binY > hMap->GetYaxis()->GetNbins())
    binY = hMap->GetYaxis()->GetNbins();
  
  return hMap->GetBinContent(binX, binY);
}

Int_t correctShapeEtaTree(Bool_t correctData, TString pathMap, TString fileNameMap, TString mapSuffix,
                          Bool_t recalculateExpecteddEdx, TString pathNameSplinesFile, TString prSplinesName,
                          TString pathTree,
                          Bool_t hasMultiplicity,
                          Bool_t correctMult,
                          TString fileNameTree = "bhess_PIDetaTree.root", TString treeName = "fTree") 
{ 
  const Double_t massProton = AliPID::ParticleMass(AliPID::kProton);
  
  if (!correctData && !recalculateExpecteddEdx) {
    std::cout << "Nothing to be done: Correction AND recalculation of expected dEdx are disabled!" << std::endl;
    return -1;
  }
  // Extract the splines, if desired
  TSpline3* splPr = 0x0;
  if (recalculateExpecteddEdx) {
    std::cout << "Loading splines to recalculate expected dEdx!" << std::endl << std::endl;
    
    TFile* fSpl = TFile::Open(pathNameSplinesFile.Data());
    if (!fSpl) {
      std::cout << "Failed to open spline file \"" << pathNameSplinesFile.Data() << "\"!" << std::endl;
      return 0x0;
    }
    
    TObjArray* TPCPIDResponse = (TObjArray*)fSpl->Get("TPCPIDResponse");
    if (!TPCPIDResponse) {
      splPr = (TSpline3*)fSpl->Get(prSplinesName.Data());
      
      // If splines are in file directly, without TPCPIDResponse object, try to load them
      if (!splPr) {
        std::cout << "Failed to load object array from spline file \"" << pathNameSplinesFile.Data() << "\"!" << std::endl;
        return 0x0;
      }
    }
    else {
      splPr = (TSpline3*)TPCPIDResponse->FindObject(prSplinesName.Data());
      
      if (!splPr) {
        std::cout << "Failed to load splines from file \"" << pathNameSplinesFile.Data() << "\"!" << std::endl;
        return 0x0;
      }
    }
  }
  else
    std::cout << "Taking dEdxExpected from Tree..." << std::endl << std::endl;
  
  if (correctMult)
    std::cout << "Correcting multiplicity dependence..." << std::endl << std::endl;
  else
    std::cout << "NOT correcting multiplicity dependence..." << std::endl << std::endl;
  
  // Extract the correction map, if desired
  TH2D* hMap = 0x0;
  if (correctData) {
    std::cout << "Loading map to correct dEdx (NOT dEdx_expected in this special case!)!" << std::endl << std::endl;
    
    TFile* fMap = TFile::Open(Form("%s/%s", pathMap.Data(), fileNameMap.Data()));
    if (!fMap)  {
      std::cout << "Failed to open file \"" << Form("%s/%s", pathMap.Data(), fileNameMap.Data()) << "\"!" << std::endl;
      return -1;
    }
    
    hMap = dynamic_cast<TH2D*>(fMap->Get(Form("hRefined%s", mapSuffix.Data())));
    if (!hMap) {
      std::cout << "Failed to load correction map!" << std::endl;
      return -1;
    }
  }
  else
    std::cout << "NOT correcting eta dependence..." << std::endl << std::endl;
  
  // Extract the data Tree		       
  TFile* fTree = TFile::Open(Form("%s/%s", pathTree.Data(), fileNameTree.Data()));
  if (!fTree)  {
    std::cout << "Failed to open file \"" << Form("%s/%s", pathTree.Data(), fileNameTree.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  TTree* tree = dynamic_cast<TTree*>(fTree->Get(treeName.Data()));
  if (!tree) {
    std::cout << "Failed to load data tree!" << std::endl << std::endl;
    return -1;
  }

  Long64_t nTreeEntries = tree->GetEntriesFast();
  
  Double_t pTPC = 0.; // Momentum at TPC inner wall
  //Double_t pT = 0.;
  Double_t dEdx = 0.; // Measured dE/dx
  Double_t dEdxExpected = 0.; // Expected dE/dx according to parametrisation
  Double_t tanTheta = 0.; // Tangens of (local) theta at TPC inner wall
  //Double_t sinAlpha = 0.; // Sine of (local) phi at TPC inner wall
  //Double_t y = 0.; // Local Y
  Double_t phiPrime = 0; // Phi prime
  UShort_t tpcSignalN = 0; // Number of TPC clusters for PID
  Int_t multiplicity = 0;
  UChar_t  pidType = 0; // Type of identification (TPC dEdx, V0, ...)
  
  tree->SetBranchAddress("pTPC", &pTPC);
  //tree->SetBranchAddress("pT", &Pt);
  tree->SetBranchAddress("dEdx", &dEdx);
  if (!recalculateExpecteddEdx)
    tree->SetBranchAddress("dEdxExpected", &dEdxExpected);
  tree->SetBranchAddress("tanTheta", &tanTheta);
  //tree->SetBranchAddress("sinAlpha", &sinAlpha);
  //tree->SetBranchAddress("y", &y);
  //TODO not needed for the moment tree->SetBranchAddress("phiPrime", &phiPrime);
  tree->SetBranchAddress("tpcSignalN", &tpcSignalN);
  if (hasMultiplicity)
    tree->SetBranchAddress("fMultiplicity", &multiplicity);
  tree->SetBranchAddress("pidType", &pidType);
  
  // Output file
  TDatime daTime;  
  TString savefileName = "";
  if (correctData) {
    savefileName = Form("%s_%sCorrectedWithMap_%s___%04d_%02d_%02d__%02d_%02d.root", fileNameTree.ReplaceAll(".root", "").Data(),
                        correctMult ? "multiplicityCorrected_" : "",
                        fileNameMap.ReplaceAll(".root", "").Data(), daTime.GetYear(), 
                        daTime.GetMonth(), daTime.GetDay(), daTime.GetHour(), daTime.GetMinute());
  }
  else if (recalculateExpecteddEdx) {
    savefileName = Form("%s_%sNewSplines___%04d_%02d_%02d__%02d_%02d.root", fileNameTree.ReplaceAll(".root", "").Data(),
                        correctMult ? "multiplicityCorrected_" : "",
                        daTime.GetYear(), daTime.GetMonth(), daTime.GetDay(), daTime.GetHour(), daTime.GetMinute());
  }
  else
    return -1;
    
    
  TFile* fSave = TFile::Open(Form("%s/%s", pathTree.Data(), savefileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", pathTree.Data(), savefileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  
  fSave->cd();
  TTree* treeCorrected = new TTree("fTree", "Tree for analysis of #eta dependence of TPC signal; #eta corrected");
  treeCorrected->Write(0, TObject::kOverwrite);
  
  treeCorrected->Branch("pTPC", &pTPC);
  //treeCorrected->Branch("pT", &Pt);
  treeCorrected->Branch("dEdx", &dEdx);
  treeCorrected->Branch("dEdxExpected", &dEdxExpected);
  treeCorrected->Branch("tpcSignalN", &tpcSignalN);
  
  treeCorrected->Branch("tanTheta", &tanTheta);
  //treeCorrected->Branch("sinAlpha", &sinAlpha);
  //treeCorrected->Branch("y", &y);
  treeCorrected->Branch("phiPrime", &phiPrime);
  if (hasMultiplicity)
    treeCorrected->Branch("fMultiplicity", &multiplicity);
  treeCorrected->Branch("pidType", &pidType);
  
  Double_t corrFactor = 0;
  
  
  TF1 corrFuncMult("corrFuncMult", "[0] + [1]*TMath::Max([4], TMath::Min(x, [3])) + [2] * TMath::Power(TMath::Max([4], TMath::Min(x, [3])), 2)",
                   0., 0.2);
  TF1 corrFuncMultTanTheta("corrFuncMultTanTheta", "[0] * (x -[2]) + [1] * (x * x - [2] * [2])", -1.5, 1.5);
  //OLD TF1 corrFuncMult("corrFuncMult", "[0] + [1]*TMath::Min(x, [3]) + [2] * TMath::Power(TMath::Min(x, [3]), 2)", 0., 0.1);
  //OLD acceptable try, but with coarser tanTheta binning and absTanTheta TF1 corrFuncMultTanTheta("corrFuncMultTanTheta", "[0] * (TMath::Abs(x) -TMath::Abs([2])) + [1] * (x * x - [2] * [2])", -1.5, 1.5);
  
  // LHC13b.pass2
  if (correctMult)
    printf("Using corr Parameters for 13b.pass2\n!");
  
  corrFuncMult.SetParameter(0, -6.27187e-06);
  corrFuncMult.SetParameter(1, -4.60649e-04);
  corrFuncMult.SetParameter(2, -4.26450e-02);
  corrFuncMult.SetParameter(3, 2.40590e-02);
  corrFuncMult.SetParameter(4, 0);
  
  corrFuncMultTanTheta.SetParameter(0, -5.338e-06);
  corrFuncMultTanTheta.SetParameter(1,  1.220e-05);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  
  /*
  // LHC11a10a
  if (correctMult)
    printf("Using corr Parameters for 11a10a\n!");
  
  corrFuncMult.SetParameter(0, 6.90133e-06);
  corrFuncMult.SetParameter(1, -1.22123e-03);
  corrFuncMult.SetParameter(2, 1.80220e-02);
  corrFuncMult.SetParameter(3, 0.1);
  corrFuncMult.SetParameter(4, 6.45306e-03);
  
  corrFuncMultTanTheta.SetParameter(0, -2.85505e-07);
  corrFuncMultTanTheta.SetParameter(1, -1.31911e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  /* OLD very good try, but with fewer pBins for the fitting
  corrFuncMult.SetParameter(0, 6.88365e-06);
  corrFuncMult.SetParameter(1, -1.22324e-03);
  corrFuncMult.SetParameter(2, 1.81625e-02);
  corrFuncMult.SetParameter(3, 0.1);
  corrFuncMult.SetParameter(4, 6.36890e-03);
  
  corrFuncMultTanTheta.SetParameter(0, -2.85505e-07);
  corrFuncMultTanTheta.SetParameter(1, -1.31911e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  */
  /*OLD good try
  corrFuncMult.SetParameter(0, 7.50321e-06);
  corrFuncMult.SetParameter(1, -1.25250e-03);
  corrFuncMult.SetParameter(2, 1.85437e-02);
  corrFuncMult.SetParameter(3, 0.1);
  corrFuncMult.SetParameter(4, 6.21192e-03);
  
  corrFuncMultTanTheta.SetParameter(0, -1.43112e-07);
  corrFuncMultTanTheta.SetParameter(1, -1.53e-06);
  corrFuncMultTanTheta.SetParameter(2, 0.3);
  */
  /* OLD acceptable try, but with coarser tanTheta binning and absTanTheta
  corrFuncMult.SetParameter(0, 6.78255e-6);
  corrFuncMult.SetParameter(1, -0.00117312);
  corrFuncMult.SetParameter(2, 0.0162423);
  corrFuncMult.SetParameter(3, 0.0563968);
  corrFuncMult.SetParameter(4, 0.00663576);
  
  corrFuncMultTanTheta.SetParameter(0, -1.85779e-6);
  corrFuncMultTanTheta.SetParameter(1, 5.40642e-7);
  corrFuncMultTanTheta.SetParameter(2, 0.35);
  */
  
  /*OLD
  corrFuncMult.SetParameter(0, 6.798e-6);
  corrFuncMult.SetParameter(1, -0.001176);
  corrFuncMult.SetParameter(2, 0.01603);
  corrFuncMult.SetParameter(3, 0.1955);
  */
  
  
  /*
  // LHC10h.pass2
  if (correctMult)
    printf("Using corr Parameters for 10h.pass2\n!");
  
  corrFuncMult.SetParameter(0, 3.21636e-07);
  corrFuncMult.SetParameter(1, -6.65876e-04);
  corrFuncMult.SetParameter(2, 1.28786e-03);
  corrFuncMult.SetParameter(3, 1.47677e-02);
  corrFuncMult.SetParameter(4, 0.);
  
  corrFuncMultTanTheta.SetParameter(0, 7.23591e-08);
  corrFuncMultTanTheta.SetParameter(1, 2.7469e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  */
  /*OLD bad try
  corrFuncMult.SetParameter(0, 2.71514e-07);
  corrFuncMult.SetParameter(1, -6.92031e-04);
  corrFuncMult.SetParameter(2, 3.56042e-03);
  corrFuncMult.SetParameter(3, 1.47497e-02);
  corrFuncMult.SetParameter(4, 0.);
  
  corrFuncMultTanTheta.SetParameter(0, 8.53204e-08);
  corrFuncMultTanTheta.SetParameter(1, 2.85591e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  */
  
  /*OLD OLD acceptable try, but with coarser tanTheta binning and absTanTheta
  corrFuncMult.SetParameter(0, 1.167e-6);
  corrFuncMult.SetParameter(1, -0.0009747);
  corrFuncMult.SetParameter(2, 0.02117);
  corrFuncMult.SetParameter(3, 0.01778);
  corrFuncMult.SetParameter(4, 0.);
  
  corrFuncMultTanTheta.SetParameter(0, 7.036e-7);
  corrFuncMultTanTheta.SetParameter(1, 1.868e-6);
  corrFuncMultTanTheta.SetParameter(2, 0.5);
  */
  
  
  /* OLD definition of multiplicity with nContriutorsToPrimVertex
  TF1 corrFuncMult("corrFuncMult", "pol2", 0, 0.1);
  corrFuncMult.SetParameter(0, 5.5e-6);
  corrFuncMult.SetParameter(1, -0.00436);
  corrFuncMult.SetParameter(2, 0.103);
  */

  progressbar(0.);
  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);

    if (recalculateExpecteddEdx) {
      //Double_t old = dEdxExpected;
      dEdxExpected = 50. * splPr->Eval(pTPC / massProton); //WARNING: What, if MIP is different from 50.? Seems not to be used (tested for pp, MC_pp, PbPb and MC_PbPb), but can in principle happen
      /*if (TMath::Abs(dEdxExpected - old) > 1e-4) {
       p *rintedSomething = kTRUE;
       printf("%.2f - ", dEdxExpected);
       printf("%.2f\n", old);
      }*/
    }
    
    if (correctData) {
      corrFactor = getEtaCorrection(hMap, tanTheta, dEdxExpected);
      
      if (corrFactor <= 0)  {
        printf("Error: Bad correction factor (%f)\n", corrFactor);
        printedSomething = kTRUE;
        continue;
      }
      
      // Correct the track dEdx and leave the expected dEdx as it is, when creating the sigma map!
      // NOTE: This is due to the method the maps are created. The track dEdx (not the expected one!)
      // is corrected to uniquely relate a momemtum bin with an expected dEdx, where the expected dEdx
      // equals the track dEdx for all eta (thanks to the correction and the re-fit of the splines).
      // Consequently, looking up the uncorrected expected dEdx at a given tanTheta yields the correct
      // sigma parameter!
      
      // In principle, during map creation, one could also calculate for each tanTheta-p bin pair the expected dEdx
      // and then fill the map at the corresponding dEdxExpected-tanTheta bin. Then, one would also need to correct
      // the dEdxExpected to create the sigma maps.
      
      // In summary, both correction methods are valid in this case, but one has to be consistent! For the different methods,
      // the maps will be slightly distorted, but overall give the same results.
      
      
      
      //Scale eta dependence -> Doesn't seem to work!!!
      //corrFactor = (corrFactor - 1.)*(2e-6*30*multiplicity + 1.0) + 1.0; 
      
      dEdx /= corrFactor;
    }
    
    if (correctMult) {
      // Multiplicity depends on pure dEdx. Therefore, correction factor depends indirectly on eta
      // => Use eta correction factor to scale dEdxExpected accordingly
      Double_t dEdxExpectedInv = 1. / (dEdxExpected * (correctData ? corrFactor : 1.));
      Double_t relSlope = corrFuncMult.Eval(dEdxExpectedInv);
      
      //Correct eta dependence of slope
      relSlope += corrFuncMultTanTheta(tanTheta);
      
      Double_t corrFactorMult = relSlope * multiplicity;
      dEdx /= 1. + corrFactorMult;
    }
    
    treeCorrected->Fill();
    
    if (i % 10000 == 0)
      progressbar(100. * (((Double_t)i) / nTreeEntries));
  }
  
  progressbar(100.);
  
  fSave->cd();
  treeCorrected->Write(0, TObject::kOverwrite);
  fSave->Close();
  
  return 0;
}