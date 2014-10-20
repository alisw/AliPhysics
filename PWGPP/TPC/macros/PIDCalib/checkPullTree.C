#include "TTree.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TSpline.h"

#include "TF1.h" //TODO NOW

#include <iostream>
#include <iomanip>

#include "THnSparseDefinitions.h"

//__________________________________________________________________________________________
void normaliseHisto(TH2D* h)
{
  h->Sumw2();
  //TODO NOW: Disabled normalisation (doesn't make sense here)
  /*
  for (Int_t x = 1; x <= h->GetXaxis()->GetNbins(); x++) {
    Double_t binWidth = h->GetXaxis()->GetBinWidth(x);
    for (Int_t y = 1; y <= h->GetYaxis()->GetNbins(); y++) {
      Double_t binError = h->GetBinError(x,y);
      h->SetBinContent(x, y, h->GetBinContent(x, y) / binWidth);
      h->SetBinError(x, y, binError / binWidth);
    }
  }*/
}


//_________________________________________________________________________________________
Int_t getBinX(TH2D* h, Double_t tanTheta)
{
  Int_t binX = h->GetXaxis()->FindBin(tanTheta);
  
  if (binX <= 0)
    binX = 1;
  if (binX > h->GetXaxis()->GetNbins())
    binX = h->GetXaxis()->GetNbins();

    return binX;
}


//_________________________________________________________________________________________
Int_t getBinY(TH2D* h, Double_t dEdxInv)
{
  Int_t binY = h->GetYaxis()->FindBin(dEdxInv);
  
  if (binY <= 0)
    binY = 1;
  if (binY > h->GetYaxis()->GetNbins())
    binY = h->GetYaxis()->GetNbins();

  return binY;
}


//_________________________________________________________________________________________
Int_t checkPullTree(TString pathTree,  TString pathNameThetaMap, TString pathNameSigmaMap,
                    TString mapSuffix, const Int_t collType /*0: pp, 1: pPb, 2: PbPb*/, const Bool_t plotPull = kTRUE,
                    const Double_t downScaleFactor = 1,
                    TString pathNameSplinesFile = "", TString prSplinesName = "",
                    TString fileNameTree = "bhess_PIDetaTree.root", TString treeName = "fTree")
{
  const Bool_t isNonPP = collType != 0;
  const Double_t massProton = AliPID::ParticleMass(AliPID::kProton);
  
  Bool_t recalculateExpecteddEdx = pathNameSplinesFile != "";
  
  TFile* f = 0x0;
	
  f = TFile::Open(Form("%s/%s", pathTree.Data(), fileNameTree.Data()));
  if (!f)  {
    std::cout << "Failed to open tree file \"" << Form("%s/%s", pathTree.Data(), fileNameTree.Data()) << "\"!" << std::endl;
    return -1;
  }
      
  // Extract the data Tree
  TTree* tree = dynamic_cast<TTree*>(f->Get(treeName.Data()));
  if (!tree) {
    std::cout << "Failed to load data tree!" << std::endl;
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

  // Extract the correction maps
  TFile* fMap = TFile::Open(pathNameThetaMap.Data());
  if (!fMap)  {
    std::cout << "Failed to open thetaMap file \"" << pathNameThetaMap.Data() << "\"! Will not additionally correct data...." << std::endl;
  }

  TH2D* hMap = 0x0;
  
  if (fMap) {
    hMap = dynamic_cast<TH2D*>(fMap->Get(Form("hRefined%s", mapSuffix.Data())));
    if (!hMap) {
      std::cout << "Failed to load theta map!" << std::endl;
      return -1;
    }
  }

  TFile* fSigmaMap = TFile::Open(pathNameSigmaMap.Data());
  if (!fSigmaMap)  {
    std::cout << "Failed to open simgaMap file \"" << pathNameSigmaMap.Data() << "\"!" << std::endl;
    return -1;
  }

  TH2D* hThetaMapSigmaPar1 = dynamic_cast<TH2D*>(fSigmaMap->Get("hThetaMapSigmaPar1"));
  if (!hThetaMapSigmaPar1) {
    std::cout << "Failed to load sigma map for par 1!" << std::endl;
    return -1;
  }

  Double_t c0 = -1;
  TNamed* c0Info = dynamic_cast<TNamed*>(fSigmaMap->Get("c0"));
  if (!c0Info) {
    std::cout << "Failed to extract c0 from file with sigma map!" << std::endl;
    return -1;
  }

  TString c0String = c0Info->GetTitle();
  c0 = c0String.Atof();
  printf("Loaded parameter 0 for sigma: %f\n\n", c0);

  if (plotPull)
    std::cout << "Plotting pull..." << std::endl << std::endl;
  else
    std::cout << "Plotting delta'..." << std::endl << std::endl;

  Long64_t nTreeEntries = tree->GetEntriesFast();

  Double_t dEdx = 0.; // Measured dE/dx
  Double_t dEdxExpected = 0.; // Expected dE/dx according to parametrisation
  Double_t tanTheta = 0.; // Tangens of (local) theta at TPC inner wall
  Double_t pTPC = 0.; // Momentum at TPC inner wall
  UShort_t tpcSignalN = 0; // Number of clusters used for dEdx
  UChar_t  pidType = 0;
  Int_t    fMultiplicity = 0;
  //Double_t phiPrime = 0;

  // Only activate the branches of interest to save processing time
  tree->SetBranchStatus("*", 0); // Disable all branches
  tree->SetBranchStatus("pTPC", 1);
  tree->SetBranchStatus("dEdx", 1);
  tree->SetBranchStatus("dEdxExpected", 1);
  tree->SetBranchStatus("tanTheta", 1);
  tree->SetBranchStatus("tpcSignalN", 1);
  tree->SetBranchStatus("pidType", 1);
  //tree->SetBranchStatus("phiPrime", 1);
  if (isNonPP)
    tree->SetBranchStatus("fMultiplicity", 1);

  
  tree->SetBranchAddress("dEdx", &dEdx);
  tree->SetBranchAddress("dEdxExpected", &dEdxExpected);
  tree->SetBranchAddress("tanTheta", &tanTheta);
  tree->SetBranchAddress("tpcSignalN", &tpcSignalN);
  tree->SetBranchAddress("pTPC", &pTPC);
  tree->SetBranchAddress("pidType", &pidType);
  //tree->SetBranchAddress("phiPrime", &phiPrime);
  if (isNonPP)
    tree->SetBranchAddress("fMultiplicity", &fMultiplicity);

  
  // Output file
  TDatime daTime;
  TString savefileName = Form("%s%s_checkPullSigma_%04d_%02d_%02d__%02d_%02d.root", fileNameTree.ReplaceAll(".root", "").Data(),
                              recalculateExpecteddEdx ? "_recalcdEdx" : "",
                              daTime.GetYear(), daTime.GetMonth(), daTime.GetDay(), daTime.GetHour(), daTime.GetMinute());

  TFile* fSave = TFile::Open(Form("%s/%s", pathTree.Data(), savefileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", pathTree.Data(), savefileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  const Double_t pBoundLow = 0.1;
  const Double_t pBoundUp = 5;

  const Int_t nBins1 = TMath::Ceil(180 / downScaleFactor);
  const Int_t nBins2 = TMath::Ceil(100 / downScaleFactor);
  const Int_t nBins3 = TMath::Ceil(60 / downScaleFactor);
  
  const Int_t nPbinsForMap = nBins1 + nBins2 + nBins3;
  Double_t binsPforMap[nPbinsForMap + 1];
  
  Double_t binWidth1 = (1.0 - pBoundLow) / nBins1;
  Double_t binWidth2 = (2.0 - 1.0 ) / nBins2;
  Double_t binWidth3 = (pBoundUp - 2.0) / nBins3;
  
  for (Int_t i = 0; i < nBins1; i++)  {
    binsPforMap[i] = pBoundLow + i * binWidth1;
  }
  for (Int_t i = nBins1, j = 0; i < nBins1 + nBins2; i++, j++)  {
    binsPforMap[i] = 1.0 + j * binWidth2;
  }
  for (Int_t i = nBins1 + nBins2, j = 0; i < nBins1 + nBins2 + nBins3; i++, j++)  {
    binsPforMap[i] = 2.0 + j * binWidth3;
  }
  binsPforMap[nPbinsForMap] = pBoundUp;

  TH2D* hPull = new TH2D("hPull", "Pull vs. p_{TPC} integrated over tan(#Theta);p_{TPC} (GeV/c);Pull", nPbinsForMap, binsPforMap, 
                         plotPull ? 120 : 240, plotPull ? -6 : -0.6, plotPull ? 6 : 0.6);
  TH2D* hPullAdditionalCorr = (TH2D*)hPull->Clone("hPullAdditionalCorr");
  hPullAdditionalCorr->SetTitle("Pull vs. p_{TPC} integrated over tan(#Theta) with additional dEdx correction w.r.t. tan(#Theta)");
  /*
  const Int_t nThetaHistos = 3;
  TH2D* hPullTheta[nThetaHistos];
  TH2D* hPullAdditionalCorrTheta[nThetaHistos];
  Double_t tThetaLow[nThetaHistos] = { 0.0, 0.4, 0.9 };
  Double_t tThetaHigh[nThetaHistos] = { 0.1, 0.5, 1.0 };
  */
  const Int_t nThetaHistos = 10;
  TH2D* hPullTheta[nThetaHistos];
  TH2D* hPullAdditionalCorrTheta[nThetaHistos];
  Double_t tThetaLow[nThetaHistos] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
  Double_t tThetaHigh[nThetaHistos] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  

  for (Int_t i = 0; i < nThetaHistos; i++)    {
    hPullTheta[i] = new TH2D(Form("hPullTheta_%d", i),
                             Form("Pull vs. p_{TPC} for %.2f <= |tan(#Theta)| < %.2f;p_{TPC} (GeV/c);Pull", tThetaLow[i], tThetaHigh[i]),
                             nPbinsForMap, binsPforMap, plotPull ? 120 : 240, plotPull ? -6 : -0.6, plotPull ? 6 : 0.6);

    hPullAdditionalCorrTheta[i] =
      new TH2D(Form("hPullAdditionalCorrTheta_%d", i),
               Form("Pull vs. p_{TPC} for %.2f <= |tan(#Theta)| < %.2f with additional dEdx correction w.r.t. tan(#Theta);p_{TPC} (GeV/c);Pull",
                    tThetaLow[i], tThetaHigh[i]),
               nPbinsForMap, binsPforMap, plotPull ? 120 : 240, plotPull ? -6 : -0.6, plotPull ? 6 : 0.6);
  }

  
  
  
  
  
  TF1 corrFuncMult("corrFuncMult", "[0] + [1]*TMath::Max([4], TMath::Min(x, [3])) + [2] * TMath::Power(TMath::Max([4], TMath::Min(x, [3])), 2)",
                   0., 0.2);
  TF1 corrFuncMultTanTheta("corrFuncMultTanTheta", "[0] * (x -[2]) + [1] * (x * x - [2] * [2])", -1.5, 1.5);
  TF1 corrFuncSigmaMult("corrFuncSigmaMul", "TMath::Max(0, [0] + [1]*TMath::Min(x, [3]) + [2] * TMath::Power(TMath::Min(x, [3]), 2))", 0., 0.2);
  
  
  // LHC13b.pass2
  if (isNonPP)
    printf("Using corr Parameters for 13b.pass2\n!");
  
  corrFuncMult.SetParameter(0, -5.906e-06);
  corrFuncMult.SetParameter(1, -5.064e-04);
  corrFuncMult.SetParameter(2, -3.521e-02);
  corrFuncMult.SetParameter(3,  2.469e-02);
  corrFuncMult.SetParameter(4, 0);
  
  corrFuncMultTanTheta.SetParameter(0, -5.32e-06);
  corrFuncMultTanTheta.SetParameter(1,  1.177e-05);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  corrFuncSigmaMult.SetParameter(0, 0.);
  corrFuncSigmaMult.SetParameter(1, 0.);
  corrFuncSigmaMult.SetParameter(2, 0.);
  corrFuncSigmaMult.SetParameter(3, 0.);
  
  
  /* OK, but PID task was not very satisfying
  corrFuncMult.SetParameter(0, -6.27187e-06);
  corrFuncMult.SetParameter(1, -4.60649e-04);
  corrFuncMult.SetParameter(2, -4.26450e-02);
  corrFuncMult.SetParameter(3, 2.40590e-02);
  corrFuncMult.SetParameter(4, 0);
  
  corrFuncMultTanTheta.SetParameter(0, -5.338e-06);
  corrFuncMultTanTheta.SetParameter(1,  1.220e-05);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  corrFuncSigmaMult.SetParameter(0, 7.89237e-05);
  corrFuncSigmaMult.SetParameter(1, -1.30662e-02);
  corrFuncSigmaMult.SetParameter(2, 8.91548e-01);
  corrFuncSigmaMult.SetParameter(3, 1.47931e-02);
  */
  
  
  /*
  // LHC11a10a
  if (isNonPP)
    printf("Using corr Parameters for 11a10a\n!");
  
  corrFuncMult.SetParameter(0, 6.90133e-06);
  corrFuncMult.SetParameter(1, -1.22123e-03);
  corrFuncMult.SetParameter(2, 1.80220e-02);
  corrFuncMult.SetParameter(3, 0.1);
  corrFuncMult.SetParameter(4, 6.45306e-03);
  
  corrFuncMultTanTheta.SetParameter(0, -2.85505e-07);
  corrFuncMultTanTheta.SetParameter(1, -1.31911e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  corrFuncSigmaMult.SetParameter(0, -4.29665e-05);
  corrFuncSigmaMult.SetParameter(1, 1.37023e-02);
  corrFuncSigmaMult.SetParameter(2, -6.36337e-01);
  corrFuncSigmaMult.SetParameter(3, 1.13479e-02);
  */
  
  /* OLD without saturation and large error for negative slopes
  corrFuncSigmaMult.SetParameter(0, -4.79684e-05);
  corrFuncSigmaMult.SetParameter(1, 1.49938e-02);
  corrFuncSigmaMult.SetParameter(2, -7.15269e-01);
  corrFuncSigmaMult.SetParameter(3, 1.06855e-02);
  */
  
  /* OLD very good try, but with fewer pBins for the fitting
  corrFuncMult.SetParameter(0, 6.88365e-06);
  corrFuncMult.SetParameter(1, -1.22324e-03);
  corrFuncMult.SetParameter(2, 1.81625e-02);
  corrFuncMult.SetParameter(3, 0.1);
  corrFuncMult.SetParameter(4, 6.36890e-03);
  
  corrFuncMultTanTheta.SetParameter(0, -2.85505e-07);
  corrFuncMultTanTheta.SetParameter(1, -1.31911e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  corrFuncSigmaMult.SetParameter(0, -4.28401e-05);
  corrFuncSigmaMult.SetParameter(1, 1.24812e-02);
  corrFuncSigmaMult.SetParameter(2, -5.28531e-01);
  corrFuncSigmaMult.SetParameter(3, 1.25147e-02);
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
  
  corrFuncSigmaMult.SetParameter(0, -2.54019e-05);
  corrFuncSigmaMult.SetParameter(1, 8.68883e-03);
  corrFuncSigmaMult.SetParameter(2, -3.36176e-01);
  corrFuncSigmaMult.SetParameter(3, 1.29230e-02);
  */
  
  /*
  // LHC10h.pass2
  if (isNonPP)
    printf("Using corr Parameters for 10h.pass2\n!");
  
  corrFuncMult.SetParameter(0, 3.21636e-07);
  corrFuncMult.SetParameter(1, -6.65876e-04);
  corrFuncMult.SetParameter(2, 1.28786e-03);
  corrFuncMult.SetParameter(3, 1.47677e-02);
  corrFuncMult.SetParameter(4, 0.);
  
  corrFuncMultTanTheta.SetParameter(0, 7.23591e-08);
  corrFuncMultTanTheta.SetParameter(1, 2.7469e-06);
  corrFuncMultTanTheta.SetParameter(2, -0.5);
  
  corrFuncSigmaMult.SetParameter(0, -1.22590e-05);
  corrFuncSigmaMult.SetParameter(1, 6.88888e-03);
  corrFuncSigmaMult.SetParameter(2, -3.20788e-01);
  corrFuncSigmaMult.SetParameter(3, 1.07345e-02);
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
  
  corrFuncSigmaMult.SetParameter(0, -6.82477e-06);
  corrFuncSigmaMult.SetParameter(1, 4.97051e-03);
  corrFuncSigmaMult.SetParameter(2, -1.64954e-01);
  corrFuncSigmaMult.SetParameter(3, 9.21061e-03);
  */
  

  //TODO NOW
  TF1* fShapeSmallP = new TF1("fShapeSmallP", "pol5", -0.4, 0.4);
  fShapeSmallP->SetParameters(1.01712, -0.0202725, -0.260692, 0.261623, 0.671854, -1.14014);
    
  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);

    if (dEdx <= 0 || dEdxExpected <= 0 || tpcSignalN <= 10)
      continue;
    /*
    Double_t pT = pTPC*TMath::Sin(-TMath::ATan(tanTheta)+TMath::Pi()/2.0);
    if ((phiPrime > 0.072/pT+TMath::Pi()/18.0-0.035 && phiPrime < 0.07/pT/pT+0.1/pT+TMath::Pi()/18.0+0.035)) 
      continue;
    */
      
    if (pidType != kMCid) {
      if (pidType == kTPCid && pTPC > 0.6)
        continue;
      if (pidType == kTPCandTOFid && (pTPC < 0.6 || pTPC > 2.0))
        continue;
      if ((collType == 2) && pidType == kTPCandTOFid && pTPC > 1.0)
        continue;// Only V0's in case of PbPb above 1.0 GeV/c
      if (pidType == kV0idPlusTOFrejected) //TODO NOW NEW
        continue;
    }
    
    if (recalculateExpecteddEdx) {
      dEdxExpected = 50. * splPr->Eval(pTPC / massProton); //WARNING: What, if MIP is different from 50.? Seems not to be used (tested for pp, MC_pp, PbPb and MC_PbPb), but can in principle happen
    }
      
    //TODO NOW
    /*
    if (TMath::Abs(tanTheta) <= 0.4) {
      Double_t p0 = fShapeSmallP->Eval(tanTheta) - 1.0; // Strength of the correction
      Double_t p1 = -9.0; // How fast the correction is turned off
      Double_t p2 = -0.209; // Turn off correction around 0.2 GeV/c
      Double_t p3 = 1.0; // Delta' for large p should be 1

      Double_t corrFactor = TMath::Erf((pTPC + p2) * p1) * p0 + p3 + p0; // Add p0 to have 1 for p3 = 1 and large pTPC
      dEdxExpected *= corrFactor;
    }*/
    
     /*TODO old unsuccessful try 
    Double_t thetaGlobalTPC = -TMath::ATan(tanTheta) + TMath::Pi() / 2.;
    Double_t pTtpc = pTPC * TMath::Sin(thetaGlobalTPC);
    Double_t pTtpcInv = (pTtpc > 0) ? 1. / pTtpc : 0;
    Double_t p0 = 1.0;
    Double_t p1 = 1./ 0.5;//TODO 2.0;
    Double_t p2 = -0.2;//TODO 0.1
    Double_t pTcorrFactor = p0 + (pTtpcInv > p1) * p2 * (pTtpcInv - p1);
    
    dEdxExpected *= pTcorrFactor;
    */
    
      
    // From the momentum (via dEdxExpected) and the tanTheta of the track, the expected dEdx can be calculated (correctedDeDxExpected).
    // If the splines are correct, this should give in average the same value as dEdx. 
    // Now valid: Maps created from corrected data with splines adopted to corrected data, so lookup should be for dEdxExpected=dEdxSplines (no further
    // eta correction) or the corrected dEdx from the track (which should ideally be = dEdxSplines)
    
    // Tested with corrected data for LHC10d.pass2: using dEdx for the lookup (which is the corrected value and should ideally be = dEdxSplines):
    // Results almost the same. Maybe slightly better for dEdxExpected.
    
    // No longer valid: Note that the maps take always the uncorrected dEdx w.r.t.
    // tanTheta, so that correctedDeDxExpected is needed here normally. However, the information for the correction will be lost at some point.
    // Therefore, dEdxExpected can be used instead and should provide a good approximation.
    Double_t c1FromSigmaMap = hThetaMapSigmaPar1->GetBinContent(getBinX(hThetaMapSigmaPar1, tanTheta), getBinY(hThetaMapSigmaPar1, 1./dEdxExpected));
    
    Double_t expectedSigma = dEdxExpected * TMath::Sqrt( c0 * c0 + (c1FromSigmaMap * c1FromSigmaMap) / tpcSignalN);
    Double_t pull = (dEdx - dEdxExpected) / (plotPull ? expectedSigma: dEdxExpected);

    // Fill pull histo
    hPull->Fill(pTPC, pull);

    Double_t tanThetaAbs = TMath::Abs(tanTheta);
    
    for (Int_t j = 0; j < nThetaHistos; j++)    {
      if (tanThetaAbs  >= tThetaLow[j] && tanThetaAbs < tThetaHigh[j])  {
        hPullTheta[j]->Fill(pTPC, pull);
      }
    }

    if (!hMap)
      continue;

    Double_t correctionFactor = 1.;
    
    if (isNonPP) {
      // 1. Correct eta dependence
      correctionFactor = hMap->GetBinContent(getBinX(hMap, tanTheta), getBinY(hMap, 1./dEdxExpected));
      
      // 2. Correct for multiplicity dependence:
      Double_t multCorrectionFactor = 1.;
      
      if (fMultiplicity > 0) {
        Double_t relSlope = corrFuncMult.Eval(1. / (dEdxExpected * correctionFactor));
        relSlope += corrFuncMultTanTheta.Eval(tanTheta);

        multCorrectionFactor = 1. + relSlope * fMultiplicity;
      }

      c1FromSigmaMap = hThetaMapSigmaPar1->GetBinContent(getBinX(hThetaMapSigmaPar1, tanTheta), getBinY(hThetaMapSigmaPar1, 1./dEdxExpected));
      
      // Multiplicity dependence of sigma depends on the real dEdx at zero multiplicity, i.e. the eta (only) corrected dEdxExpected value has to be used
      // since all maps etc. have been created for ~zero multiplicity
      Double_t relSigmaSlope = corrFuncSigmaMult.Eval(1. / (dEdxExpected * correctionFactor));
      Double_t multSigmaCorrectionFactor = 1. + relSigmaSlope * fMultiplicity;
      
      dEdxExpected *= correctionFactor * multCorrectionFactor; 
      
      expectedSigma = dEdxExpected * TMath::Sqrt( c0 * c0 + (c1FromSigmaMap * c1FromSigmaMap) / tpcSignalN);
      expectedSigma *= multSigmaCorrectionFactor;
      
      pull = (dEdx - dEdxExpected) / (plotPull ? expectedSigma: dEdxExpected);
    }
    else {
      correctionFactor = hMap->GetBinContent(getBinX(hMap, tanTheta), getBinY(hMap, 1./dEdxExpected));
      
      c1FromSigmaMap = hThetaMapSigmaPar1->GetBinContent(getBinX(hThetaMapSigmaPar1, tanTheta), getBinY(hThetaMapSigmaPar1, 1./dEdxExpected));
   
      dEdxExpected *= correctionFactor; // If data is not corrected, but the sigma map is for corrected data, re-do analysis with corrected dEdx
      
      expectedSigma = dEdxExpected * TMath::Sqrt( c0 * c0 + (c1FromSigmaMap * c1FromSigmaMap) / tpcSignalN);
      pull = (dEdx - dEdxExpected) / (plotPull ? expectedSigma: dEdxExpected);
    }

    pull = (dEdx - dEdxExpected) / (plotPull ? expectedSigma: dEdxExpected);

    hPullAdditionalCorr->Fill(pTPC, pull);

    for (Int_t j = 0; j < nThetaHistos; j++)    {
      if (tanThetaAbs  >= tThetaLow[j] && tanThetaAbs < tThetaHigh[j])  {
        hPullAdditionalCorrTheta[j]->Fill(pTPC, pull);
      }
    }
  }
/*
  // Mean, Sigma, chi^2/NDF of pull of different theta bins and all in one plot
  TCanvas* canvPullMean = new TCanvas("canvPullMean", "canvPullMean", 100,10,1380,800);
  canvPullMean->SetLogx(kTRUE);
  canvPullMean->SetGridx(kTRUE);
  canvPullMean->SetGridy(kTRUE);
  TCanvas* canvPullSigma = new TCanvas("canvPullSigma", "canvPullSigma", 100,10,1380,800);
  canvPullSigma->SetLogx(kTRUE);
  canvPullSigma->SetGridx(kTRUE);
  canvPullSigma->SetGridy(kTRUE);
  TCanvas* canvPullChi2 = new TCanvas("canvPullChi2", "canvPullChi2", 100,10,1380,800);
  canvPullChi2->SetLogx(kTRUE);
  canvPullChi2->SetGridx(kTRUE);
  canvPullChi2->SetGridy(kTRUE);
  

  TCanvas* canvPull[nThetaHistos + 1];
  for (Int_t i = 0, j = nThetaHistos; i < nThetaHistos + 1; i++, j--)  {
    canvPull[i] = new TCanvas(Form("canvPull_%d", i), "canvPull", 100,10,1380,800);
    canvPull[i]->cd();
    canvPull[i]->SetLogx(kTRUE);
    canvPull[i]->SetLogz(kTRUE);
    canvPull[i]->SetGrid(kTRUE, kTRUE);

    TH2D* hTemp = 0x0;
    TString thetaString = "";
    if (i == nThetaHistos)  {
      hTemp = hPull;
      thetaString = "tan(#Theta) integrated";
    }
    else {
      hTemp = hPullTheta[i];
      thetaString = Form("%.2f #leq |tan(#Theta)| < %.2f", tThetaLow[i], tThetaHigh[i]);
    }
    
    normaliseHisto(hTemp);
    hTemp->FitSlicesY();
    hTemp->GetYaxis()->SetNdivisions(12);
    hTemp->GetXaxis()->SetMoreLogLabels(kTRUE);
    TH1D* hTempMean = (TH1D*)gDirectory->Get(Form("%s_1", hTemp->GetName()));
    hTempMean->SetTitle(Form("mean(pull), %s", thetaString.Data()));
    hTempMean->GetXaxis()->SetMoreLogLabels(kTRUE);
    hTempMean->SetLineWidth(2);
    hTempMean->SetMarkerStyle(20);
    TH1D* hTempSigma = (TH1D*)gDirectory->Get(Form("%s_2", hTemp->GetName()));
    hTempSigma->SetTitle(Form("#sigma(pull), %s", thetaString.Data()));
    hTempSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
    hTempSigma->SetLineColor(kMagenta);
    hTempSigma->SetMarkerStyle(20);
    hTempSigma->SetMarkerColor(kMagenta);
    hTempSigma->SetLineWidth(2);
    TH1D* hTempChi2 = (TH1D*)gDirectory->Get(Form("%s_chi2", hTemp->GetName()));
    hTempChi2->SetTitle(Form("#chi^{2} / NDF (pull), %s", thetaString.Data()));
    hTempChi2->GetXaxis()->SetMoreLogLabels(kTRUE);
    hTempChi2->SetLineColor(kMagenta + 2);
    hTempChi2->SetMarkerStyle(20);
    hTempChi2->SetMarkerColor(kMagenta + 2);
    hTempChi2->SetLineWidth(2);

    hTemp->DrawCopy("colz");
    hTempMean->DrawCopy("same");
    hTempSigma->DrawCopy("same");
    hTempChi2->Scale(-1./10.);
    hTempChi2->DrawCopy("same");
    hTempChi2->Scale(-10.);

    canvPullMean->cd();
    hTempMean->SetLineColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
    hTempMean->SetMarkerColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
    hTempMean->DrawCopy((i == 0 ? "" : "same"));

    canvPullSigma->cd();
    hTempSigma->SetLineColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
    hTempSigma->SetMarkerColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
    hTempSigma->DrawCopy((i == 0 ? "" : "same"));

    canvPullChi2->cd();
    hTempChi2->SetLineColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
    hTempChi2->SetMarkerColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
    hTempChi2->DrawCopy((i == 0 ? "" : "same"));
  }

  canvPullMean->BuildLegend();
  canvPullSigma->BuildLegend();
  canvPullChi2->BuildLegend();
*/
  // Histograms with additional correction
  TCanvas* canvPullMeanCorr = 0x0;
  TCanvas* canvPullSigmaCorr = 0x0;
  TCanvas* canvPullChi2Corr = 0x0;
  TCanvas* canvPullCorr[nThetaHistos + 1];
  for (Int_t i = 0; i < nThetaHistos + 1; i++) 
    canvPullCorr[i] = 0x0;
  
  if (hMap) {
    // Mean, Sigma, chi^2/NDF of pull of different theta bins and all in one plot
    canvPullMeanCorr = new TCanvas("canvPullMeanCorr", "canvPullMeanCorr", 100,10,1380,800);
    canvPullMeanCorr->SetLogx(kTRUE);
    canvPullMeanCorr->SetGridx(kTRUE);
    canvPullMeanCorr->SetGridy(kTRUE);
    canvPullSigmaCorr = new TCanvas("canvPullSigmaCorr", "canvPullSigmaCorr", 100,10,1380,800);
    canvPullSigmaCorr->SetLogx(kTRUE);
    canvPullSigmaCorr->SetGridx(kTRUE);
    canvPullSigmaCorr->SetGridy(kTRUE);
    canvPullChi2Corr = new TCanvas("canvPullChi2Corr", "canvPullChi2Corr", 100,10,1380,800);
    canvPullChi2Corr->SetLogx(kTRUE);
    canvPullChi2Corr->SetGridx(kTRUE);
    canvPullChi2Corr->SetGridy(kTRUE);
    
    for (Int_t i = 0, j = nThetaHistos; i < nThetaHistos + 1; i++, j--)  {
      canvPullCorr[i] = new TCanvas(Form("canvPullCorr_%d", i), "canvPullCorr", 100,10,1380,800);
      canvPullCorr[i]->cd();
      canvPullCorr[i]->SetLogx(kTRUE);
      canvPullCorr[i]->SetLogz(kTRUE);
      canvPullCorr[i]->SetGrid(kTRUE, kTRUE);

      TH2D* hTemp = 0x0;
      TString thetaString = "";
      
      if (i == nThetaHistos)  {
        hTemp = hPullAdditionalCorr;
        thetaString = "tan(#Theta) integrated";
      }
      else    {
        hTemp = hPullAdditionalCorrTheta[i];
        thetaString = Form("%.2f #leq |tan(#Theta)| < %.2f", tThetaLow[i], tThetaHigh[i]);
      }

      normaliseHisto(hTemp);
      hTemp->FitSlicesY();
      hTemp->GetYaxis()->SetNdivisions(12);
      hTemp->GetXaxis()->SetMoreLogLabels(kTRUE);
      TH1D* hTempMean = (TH1D*)gDirectory->Get(Form("%s_1", hTemp->GetName()));
      hTempMean->SetTitle(Form("mean(pull), %s", thetaString.Data()));
      hTempMean->GetXaxis()->SetMoreLogLabels(kTRUE);
      hTempMean->SetLineWidth(2);
      hTempMean->SetMarkerStyle(20);
      TH1D* hTempSigma = (TH1D*)gDirectory->Get(Form("%s_2", hTemp->GetName()));
      hTempSigma->SetTitle(Form("#sigma(pull), %s", thetaString.Data()));
      hTempSigma->GetXaxis()->SetMoreLogLabels(kTRUE);
      hTempSigma->SetLineColor(kMagenta);
      hTempSigma->SetMarkerStyle(20);
      hTempSigma->SetMarkerColor(kMagenta);
      hTempSigma->SetLineWidth(2);
      TH1D* hTempChi2 = (TH1D*)gDirectory->Get(Form("%s_chi2", hTemp->GetName()));
      hTempChi2->SetTitle(Form("#chi^{2} / NDF (pull), %s", thetaString.Data()));
      hTempChi2->GetXaxis()->SetMoreLogLabels(kTRUE);
      hTempChi2->SetLineColor(kMagenta + 2);
      hTempChi2->SetMarkerStyle(20);
      hTempChi2->SetMarkerColor(kMagenta + 2);
      hTempChi2->SetLineWidth(2);

      hTemp->DrawCopy("colz");
      hTempMean->DrawCopy("same");
      hTempSigma->DrawCopy("same");
      hTempChi2->Scale(-1./10.);
      hTempChi2->DrawCopy("same");
      hTempChi2->Scale(-10.);
  
      canvPullMeanCorr->cd();
      hTempMean->SetLineColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
      hTempMean->SetMarkerColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
      hTempMean->DrawCopy((i == 0 ? "" : "same"));
      
      canvPullSigmaCorr->cd();
      hTempSigma->SetLineColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
      hTempSigma->SetMarkerColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
      hTempSigma->DrawCopy((i == 0 ? "" : "same"));
      
      canvPullChi2Corr->cd();
      hTempChi2->SetLineColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
      hTempChi2->SetMarkerColor(1 + ((j >= 9) ? (39 + 2 * (j - 9)) : j));
      hTempChi2->DrawCopy((i == 0 ? "" : "same"));
    }
    
    canvPullMeanCorr->BuildLegend();
    canvPullSigmaCorr->BuildLegend();
    canvPullChi2Corr->BuildLegend();
  }
  
  
  
  
  
  fSave->cd();
  /*canvPullMean->Write();
  canvPullSigma->Write();
  canvPullChi2->Write();
  
  for (Int_t  i = 0; i < nThetaHistos + 1; i++) {
    canvPull[i]->Write();
  }*/
  
  canvPullMeanCorr->Write();
  canvPullSigmaCorr->Write();
  canvPullChi2Corr->Write();
  
  for (Int_t  i = 0; i < nThetaHistos + 1; i++) {
    canvPullCorr[i]->Write();
  }

  TNamed* info = new TNamed(Form("Theta map: %s\n\nSigma map: %s\n\nSplines file: %s\n\nSplines name: %s", pathNameThetaMap.Data(), 
                                 pathNameSigmaMap.Data(), pathNameSplinesFile.Data(), prSplinesName.Data()),
                            "info");
  info->Write();
  fSave->Close();
  
  return 0;
}
