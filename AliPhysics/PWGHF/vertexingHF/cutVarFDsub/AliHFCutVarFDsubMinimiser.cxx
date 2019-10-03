#include "AliHFCutVarFDsubMinimiser.h"

#include "TH1F.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "Riostream.h"
#include "TCanvas.h"

/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubMinimiser);
/// \endcond


AliHFCutVarFDsubMinimiser::AliHFCutVarFDsubMinimiser()
  : TObject()
  , fhRawYields(0x0)
  , fhEffPrompt(0x0)
  , fhEffFD(0x0)
  , fNiterations(0)
  , fUseWeights(kFALSE)
  , fRelSystEffErr(0.)
  , fPromptYield(-1.)
  , fPromptYieldErr(-1.)
  , fFDYield(-1.)
  , fFDYieldErr(-1.)
  , fhResiduals(0x0)
  , fhPulls(0x0)
  , fhfPrompt(0x0)
  , fhfPromptRaw(0x0)
  , fMinimised(kFALSE)
  , fnSim(0)
  , fhIncDistError(0x0)
{
  /// Default constructor
}


AliHFCutVarFDsubMinimiser::AliHFCutVarFDsubMinimiser(TH1F* hRawYields, TH1F* hEffPrompt, TH1F* hEffFD,
                                                     UInt_t method/*=0*/, UInt_t nIterations/*=10*/,
                                                     Bool_t useWeights/*=kTRUE*/,
                                                     Double_t relSystEffErr/*=0.*/,
                                                     Int_t nSim/*= 1000*/)
  : TObject()
  , fhRawYields(hRawYields)
  , fhEffPrompt(hEffPrompt)
  , fhEffFD(hEffFD)
  , fNiterations(nIterations)
  , fUseWeights(useWeights)
  , fRelSystEffErr(relSystEffErr)
  , fPromptYield(-1.)
  , fPromptYieldErr(-1.)
  , fFDYield(-1.)
  , fFDYieldErr(-1.)
  , fhResiduals(0x0)
  , fhPulls(0x0)
  , fhfPrompt(0x0)
  , fhfPromptRaw(0x0)
  , fMinimised(kFALSE)
  , fnSim(nSim)
  , fhIncDistError(0x0)
{
  /// Constructor

  Int_t nCutSets = fhRawYields->GetNbinsX();

  fhResiduals  = new TH1F("hResiduals",  ";Cut Set; Residuals (a.u.)",     nCutSets, 0., (Double_t)nCutSets+1.);
  fhPulls      = new TH1F("hPulls",      ";Cut Set; Pulls (a.u.)",         nCutSets, 0., (Double_t)nCutSets+1.);
  fhfPrompt    = new TH1F("hfPrompt",    ";Cut Set; f_{Prompt} (a.u.)",    nCutSets, 0., (Double_t)nCutSets+1.);
  fhfPromptRaw = new TH1F("hfPromptRaw", ";Cut Set; f_{PromptRaw} (a.u.)", nCutSets, 0., (Double_t)nCutSets+1.);

  switch (method) {
  case 0:
    fMinimised = MinimiseDefault();
    break;
  case 1:
    fMinimised = MinimiseInCentre();
    break;
  default:
    fMinimised = kFALSE;
  }
}


AliHFCutVarFDsubMinimiser::~AliHFCutVarFDsubMinimiser() {
  delete fhResiduals;
  fhResiduals = 0x0;
  delete fhPulls;
  fhPulls = 0x0;
  delete fhPulls;
  fhPulls = 0x0;
  delete fhPulls;
  fhPulls = 0x0;
  delete fhIncDistError;
  fhIncDistError = 0x0;
}


Bool_t AliHFCutVarFDsubMinimiser::MinimiseDefault() {
  /// Minimiser default method

  Int_t nCutSets = fhRawYields->GetNbinsX();
  TMatrixD mRawYield   = TMatrixD(nCutSets,        1);
  TMatrixD mWeights    = TMatrixD(nCutSets, nCutSets);
  TMatrixD mEff        = TMatrixD(nCutSets,        2);
  TMatrixD mCorrYield  = TMatrixD(       2,        1);
  TMatrixD mCovariance = TMatrixD(       2,        2);

  for(UInt_t iter=0; iter<fNiterations; ++iter) {
    for(Int_t iCutSet=0; iCutSet<nCutSets;++iCutSet) {
      // check the raw yields
      if(fhRawYields->GetBinContent(iCutSet+1) < 1.e-9 ||
         fhEffPrompt->GetBinContent(iCutSet+1) < 1.e-9 ||
         fhEffFD->GetBinContent(iCutSet+1)     < 1.e-9) {
        Printf("Zero or negative yield or efficiency, cannot proceed (cut set %d)", iCutSet);
        return kFALSE;
      }

      mRawYield(iCutSet, 0) = fhRawYields->GetBinContent(iCutSet+1); // fill raw yields

      // caculate weights
      if(!fUseWeights) mWeights(iCutSet, iCutSet) = 1.; // not using weights
      else {
        if(iter==0) {
          mWeights(iCutSet,iCutSet) = 1./(fhRawYields->GetBinError(iCutSet+1)*fhRawYields->GetBinError(iCutSet+1));
        }
        else {
          Double_t errEffPrompt = (fhEffPrompt->GetBinError(iCutSet+1)
                                   +fRelSystEffErr*fhEffPrompt->GetBinContent(iCutSet+1)) * mCorrYield(0,0);
          Double_t errEffFD     = (fhEffFD->GetBinError(iCutSet+1)
                                   +fRelSystEffErr*fhEffFD->GetBinContent(iCutSet+1)) * mCorrYield(1,0);
          Double_t errRawYield  = fhRawYields->GetBinError(iCutSet+1);
          mWeights(iCutSet, iCutSet) = 1./(errRawYield*errRawYield+errEffPrompt*errEffPrompt+errEffFD*errEffFD);
        }
      }
      mEff(iCutSet, 0) = fhEffPrompt->GetBinContent(iCutSet+1);
      mEff(iCutSet, 1) = fhEffFD->GetBinContent(iCutSet+1);
    }

    // minimisation
    TMatrixD mEffT = TMatrixD(mEff);
    mEffT.T();

    mCovariance = (mEffT*mWeights)*mEff;

    mCovariance.InvertFast(); // using closed formula

    TMatrixD mR = (mEffT*mWeights)*mRawYield;
    if(mCovariance(0,0)<1.e-36 || mCovariance(1,1)<1.e-36){
      Printf("Close to zero or negative diagonal element in covariance matrix: likely inversion failed, cannot proceed!");
      return kFALSE;
    }
    mCorrYield = mCovariance*mR;
  }
  if(mCorrYield(0,0)<1.e-24 || mCorrYield(0,0)>1.e+99) {
    Printf("Meaningless yield values");
    return kFALSE;
  }
  else {
    // Resulting yields (+their errors)
    fPromptYield    = mCorrYield(0,0);
    fPromptYieldErr = TMath::Sqrt(mCovariance(0,0));
    fFDYield        =  mCorrYield(1,0);
    fFDYieldErr     = TMath::Sqrt(mCovariance(1,1));

    // Residuals, pulls, fPrompt, fPrompt*
    for(Int_t iCutSet=0;iCutSet<nCutSets;++iCutSet) {
      fhResiduals->SetBinContent(iCutSet+1,mRawYield(iCutSet,0)-mEff(iCutSet,0)*mCorrYield(0,0)-mEff(iCutSet,1)*mCorrYield(1,0));
      fhResiduals->SetBinError(iCutSet+1,fhRawYields->GetBinError(iCutSet+1));// not real error, there is a correlation neglected
      //      fhResiduals->SetBinError(iCutSet+1,TMath::Sqrt(fhRawYields->GetBinError(iCutSet+1)*fhRawYields->GetBinError(iCutSet+1)+mEff(iCutSet,0)*mEff(iCutSet,0)*mCovariance(0,0)+mEff(iCutSet,1)*mEff(iCutSet,1)*mCovariance(1,1)));// NOT CORRECT... NEED TO THINK ABOUT IT

      fhPulls->SetBinContent(iCutSet+1,fhResiduals->GetBinContent(iCutSet+1)/fhResiduals->GetBinError(iCutSet+1));

      fhfPromptRaw->SetBinContent(iCutSet+1,mEff(iCutSet,0)*mCorrYield(0,0)/mRawYield(iCutSet,0));
      fhfPromptRaw->SetBinError(  iCutSet+1,mEff(iCutSet,0)*TMath::Sqrt(mCovariance(0,0))/mRawYield(iCutSet,0));

      fhfPrompt->SetBinContent(iCutSet+1,mEff(iCutSet,0)*mCorrYield(0,0)/(mEff(iCutSet,0)*mCorrYield(0,0)+mEff(iCutSet,1)*mCorrYield(1,0)));
      fhfPrompt->SetBinError(  iCutSet+1,mEff(iCutSet,0)*TMath::Sqrt(mCovariance(0,0))/(mEff(iCutSet,0)*mCorrYield(0,0)+mEff(iCutSet,1)*mCorrYield(1,0)));// ERROR NOT PROPERLY CALCULATED
    }
  }
  return kTRUE;
}

Bool_t AliHFCutVarFDsubMinimiser::MinimiseInCentre() {

  const Int_t nCutSets = fhRawYields->GetNbinsX();

  Bool_t incentre = kFALSE;

  Double_t effPrompt[nCutSets];
  Double_t effFD[nCutSets];
  Double_t rawYield[nCutSets];

  Double_t Ncorr[2] = {0,0};
  Double_t fprompt[nCutSets];
  Double_t fpromptraw[nCutSets];
  Double_t residuals[nCutSets];
  
  for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
    effPrompt[iCutSet]=fhEffPrompt->GetBinContent(iCutSet+1);
    effFD[iCutSet]=fhEffFD->GetBinContent(iCutSet+1);
    rawYield[iCutSet]=fhRawYields->GetBinContent(iCutSet+1);
    fprompt[iCutSet]=0;
    fpromptraw[iCutSet]=0;
    residuals[iCutSet]=0;
  }

  incentre = InCentre(effPrompt,effFD,rawYield,Ncorr,fprompt,fpromptraw,residuals);

  if(!incentre)
    return kFALSE;

  fPromptYield = Ncorr[0];
  fFDYield = Ncorr[1];

  //evaluate error with a toy MC
  Double_t rndmEffPrompt[nCutSets];
  Double_t rndmEffFD[nCutSets];
  Double_t rndmRawYield[nCutSets];

  Int_t nBins = 200;
  Double_t FDbins[2] = {fFDYield-4*TMath::Abs(fFDYield),fFDYield+4*TMath::Abs(fFDYield)};
  Double_t Promptbins[2] = {fPromptYield-4*TMath::Abs(fPromptYield),fPromptYield+4*TMath::Abs(fPromptYield)};

  fhIncDistError = new TH2F("hIncDistError",";Cut Set; Distribution for incentre error",nBins,FDbins[0],FDbins[1],nBins,Promptbins[0],Promptbins[1]);
  //xaxis->FD, yaxis->Prompt
  fhIncDistError->GetXaxis()->SetTitle("N_{FD}");
  fhIncDistError->GetYaxis()->SetTitle("N_{Prompt}");

  Double_t rndmfprompt[nCutSets];
  Double_t rndmfpromptraw[nCutSets];
  Double_t rndmresiduals[nCutSets];

  TH1F** hFpromptDistError = new TH1F*[nCutSets]; 
  TH1F** hFpromptRawDistError = new TH1F*[nCutSets];
  TH1F** hResidualsDistError = new TH1F*[nCutSets];

  Double_t resmin[nCutSets];
  Double_t resmax[nCutSets];
  
  for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
    rndmfprompt[iCutSet] = 0;
    rndmfpromptraw[iCutSet] = 0;
    rndmresiduals[iCutSet] = 0;
    
    resmin[iCutSet] = residuals[iCutSet]-10*TMath::Abs(residuals[iCutSet]);
    resmax[iCutSet] = residuals[iCutSet]+10*TMath::Abs(residuals[iCutSet]);

    hFpromptDistError[iCutSet] = new TH1F("hFpromptDistError","",200,0,2);
    hFpromptRawDistError[iCutSet] = new TH1F("hFpromptRawDistError","",200,0,2);
    hResidualsDistError[iCutSet] = new TH1F("hResidualsDistError","",200,resmin[iCutSet],resmax[iCutSet]);
  }
  
  for(Int_t iSim=0; iSim<fnSim; iSim++) {
    for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
      rndmEffPrompt[iCutSet]=gRandom->Gaus(fhEffPrompt->GetBinContent(iCutSet+1), fhEffPrompt->GetBinError(iCutSet+1));
      rndmEffFD[iCutSet]=gRandom->Gaus(fhEffFD->GetBinContent(iCutSet+1),fhEffFD->GetBinError(iCutSet+1));
      rndmRawYield[iCutSet]=gRandom->Gaus(fhRawYields->GetBinContent(iCutSet+1),fhRawYields->GetBinError(iCutSet+1));
    }
    
    incentre = kFALSE;
    incentre = InCentre(rndmEffPrompt,rndmEffFD,rndmRawYield,Ncorr,rndmfprompt,rndmfpromptraw,rndmresiduals);
    
    if(incentre) {
      fhIncDistError->Fill(Ncorr[1],Ncorr[0]);
      for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
        hFpromptDistError[iCutSet]->Fill(rndmfprompt[iCutSet]);
        hFpromptRawDistError[iCutSet]->Fill(rndmfpromptraw[iCutSet]);
        hResidualsDistError[iCutSet]->Fill(rndmresiduals[iCutSet]);
      }
    }
  }
  
  TH1F *hProjX = (TH1F*)fhIncDistError->ProjectionX();
  TH1F *hProjY = (TH1F*)fhIncDistError->ProjectionY();

  fPromptYieldErr = hProjY->GetRMS();
  fFDYieldErr = hProjX->GetRMS();

  Double_t errfprompt[nCutSets];
  Double_t errfpromptraw[nCutSets];
  Double_t errresiduals[nCutSets];
  
  //fPrompt, fPrompt* - toy MC errors
  for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
                        
    errfprompt[iCutSet] = hFpromptDistError[iCutSet]->GetRMS();
    errfpromptraw[iCutSet] = hFpromptRawDistError[iCutSet]->GetRMS();
    errresiduals[iCutSet] = hResidualsDistError[iCutSet]->GetRMS();
      
    fhfPromptRaw->SetBinContent(iCutSet+1,fpromptraw[iCutSet]);
    fhfPromptRaw->SetBinError(iCutSet+1,errfpromptraw[iCutSet]);

    fhfPrompt->SetBinContent(iCutSet+1,fprompt[iCutSet]);
    fhfPrompt->SetBinError(iCutSet+1,errfprompt[iCutSet]);

    fhResiduals->SetBinContent(iCutSet+1,residuals[iCutSet]);
    fhResiduals->SetBinError(iCutSet+1,errresiduals[iCutSet]);
    
    fhPulls->SetBinContent(iCutSet+1,residuals[iCutSet]/hResidualsDistError[iCutSet]->GetRMS());
    
    delete hFpromptDistError[iCutSet];
    delete hFpromptRawDistError[iCutSet];
    delete hResidualsDistError[iCutSet];
  }

  delete[] hFpromptDistError;
  delete[] hFpromptRawDistError;
  delete[] hResidualsDistError;
    
  delete hProjX;
  delete hProjY;

  return kTRUE;
  
}

Bool_t AliHFCutVarFDsubMinimiser::InCentre(Double_t* effPrompt,
                                           Double_t* effFD,
                                           Double_t* rawYield,
                                           Double_t* Ncorr,
                                           Double_t* fprompt,
                                           Double_t* fpromptraw,
                                           Double_t* residuals) {

  Int_t nCutSets = fhRawYields->GetNbinsX();

  TMatrixD mEff(2,2);
  TMatrixD mInt(2,1);
  TMatrixD mRaw(2,1);

  Double_t det[1];

  if(nCutSets == 2)//determinated system
  {
    for(Int_t iCutSet=0; iCutSet<nCutSets; iCutSet++) {
      mEff(iCutSet,0) = effPrompt[iCutSet];
      mEff(iCutSet,1) = effFD[iCutSet];
      mRaw(iCutSet,0) = rawYield[iCutSet];
    }

    mEff.InvertFast(det);

    if(det[0]==0) {
      Printf("The system has no solution");
      return kFALSE;
    }

    mInt = mEff*mRaw;

    //calculate the incentre and relatives quantities
    Ncorr[0] = mInt(0,0);
    Ncorr[1] = mInt(1,0);

    for(Int_t iCutSet=0; iCutSet<nCutSets; ++iCutSet) {
      fprompt[iCutSet] = Ncorr[0]*effPrompt[iCutSet]/(Ncorr[0]*effPrompt[iCutSet]+Ncorr[1]*effFD[iCutSet]);
      fpromptraw[iCutSet] = Ncorr[0]*effPrompt[iCutSet]/rawYield[iCutSet];
      residuals[iCutSet] = Ncorr[0]*effPrompt[iCutSet]+Ncorr[1]*effFD[iCutSet]-rawYield[iCutSet];
    }
    
    return kTRUE;
  }

  else if(nCutSets == 3) {
    //calculate 3 intersections
    //intersection between line 1 and 2
    Double_t NFDInt12 = 0;
    Double_t NPromptInt12 = 0;

    mEff(0,0) = effPrompt[0];
    mEff(0,1) = effFD[0];
    mRaw(0,0) = rawYield[0];

    mEff(1,0) = effPrompt[1];
    mEff(1,1) = effFD[1];
    mRaw(1,0) = rawYield[1];

    mEff.InvertFast(det);

    if(det[0]==0) {
      Printf("Line 1 and line 2 have no intersection");
      return kFALSE;
    }

    mInt = mEff*mRaw;
    NPromptInt12 = mInt(0,0);
    NFDInt12 = mInt(1,0);

    //intersection between line 1 and 3
    Double_t NFDInt13 = 0;
    Double_t NPromptInt13 = 0;

    mEff.InvertFast();
    mEff(1,0) = effPrompt[2];
    mEff(1,1) = effFD[2];
    mRaw(1,0) = rawYield[2];

    mEff.InvertFast(det);

    if(det[0]==0) {
      Printf("Line 1 and line 3 have no intersection");
      return kFALSE;
    }

    mInt = mEff*mRaw;
    NPromptInt13 = mInt(0,0);
    NFDInt13 = mInt(1,0);

    //intersection between line 2 and 3
    Double_t NFDInt23 = 0;
    Double_t NPromptInt23 = 0;

    mEff.InvertFast();
    mEff(0,0) = effPrompt[1];
    mEff(0,1) = effFD[1];
    mRaw(0,0) = rawYield[1];

    mEff.InvertFast(det);

    if(det[0]==0) {
      Printf("Line 2 and line 3 have no intersection");
      return kFALSE;
    }

    mInt = mEff*mRaw;
    NPromptInt23 = mInt(0,0);
    NFDInt23 = mInt(1,0);

    //calculate the sides of the triangle
    Double_t side1 = TMath::Sqrt((NFDInt12-NFDInt13)*(NFDInt12-NFDInt13)+(NPromptInt12-NPromptInt13)*(NPromptInt12-NPromptInt13));
    Double_t side2 = TMath::Sqrt((NFDInt23-NFDInt12)*(NFDInt23-NFDInt12)+(NPromptInt23-NPromptInt12)*(NPromptInt23-NPromptInt12));
    Double_t side3 = TMath::Sqrt((NFDInt23-NFDInt13)*(NFDInt23-NFDInt13)+(NPromptInt23-NPromptInt13)*(NPromptInt23-NPromptInt13));

    //calculate the incentre and relatives quantities
    Ncorr[0] = (NPromptInt12*side3+NPromptInt23*side1+NPromptInt13*side2)/(side1+side2+side3); //0 -> NPROMPT
    Ncorr[1] = (NFDInt12*side3+NFDInt23*side1+NFDInt13*side2)/(side1+side2+side3); //1 -> NFD;

    for(Int_t iCutSet=0; iCutSet<nCutSets; ++iCutSet) {
      fprompt[iCutSet] = Ncorr[0]*effPrompt[iCutSet]/(Ncorr[0]*effPrompt[iCutSet]+Ncorr[1]*effFD[iCutSet]);
      fpromptraw[iCutSet] = Ncorr[0]*effPrompt[iCutSet]/rawYield[iCutSet];
      residuals[iCutSet] = Ncorr[0]*effPrompt[iCutSet]+Ncorr[1]*effFD[iCutSet]-rawYield[iCutSet];
    }
    
    return kTRUE;
  }
  else {
    Printf("The number of sets it's different from 2 or 3, it's impossible to evaluate the incentre");
    return kFALSE;
  }
}

AliHFCutVarFDsubMinimiser::AliHFCutVarFDsubMinimiser(const AliHFCutVarFDsubMinimiser& m)
  : TObject()
  , fhRawYields(m.fhRawYields)
  , fhEffPrompt(m.fhEffPrompt)
  , fhEffFD(m.fhEffFD)
  , fNiterations(m.fNiterations)
  , fUseWeights(m.fUseWeights)
  , fRelSystEffErr(m.fRelSystEffErr)
  , fPromptYield(m.fPromptYield)
  , fPromptYieldErr(m.fPromptYieldErr)
  , fFDYield(m.fFDYield)
  , fFDYieldErr(m.fFDYieldErr)
  , fnSim(m.fnSim)
  , fhResiduals((TH1F*)m.fhResiduals->Clone())
  , fhPulls((TH1F*)m.fhPulls->Clone())
  , fhfPrompt((TH1F*)m.fhfPrompt->Clone())
  , fhfPromptRaw((TH1F*)m.fhfPromptRaw->Clone())
  , fhIncDistError((TH2F*)m.fhIncDistError->Clone())
{
  /// Copy constructor
}


AliHFCutVarFDsubMinimiser AliHFCutVarFDsubMinimiser::operator=(const AliHFCutVarFDsubMinimiser& m)
{
  /// Assignment operator
  if (this != &m) {
    fhRawYields = m.fhRawYields;
    fhEffPrompt = m.fhEffPrompt;
    fhEffFD = m.fhEffFD;
    fNiterations = m.fNiterations;
    fUseWeights = m.fUseWeights;
    fRelSystEffErr = m.fRelSystEffErr;
    fPromptYield = m.fPromptYield;
    fPromptYieldErr = m.fPromptYieldErr;
    fFDYield = m.fFDYield;
    fFDYieldErr = m.fFDYieldErr;
    fnSim = m.fnSim;
    fhResiduals = (TH1F*)m.fhResiduals->Clone();
    fhPulls = (TH1F*)m.fhPulls->Clone();
    fhfPrompt = (TH1F*)m.fhfPrompt->Clone();
    fhfPromptRaw = (TH1F*)m.fhfPromptRaw->Clone();
    fhIncDistError = (TH2F*)m.fhIncDistError->Clone();
  }
  return *this;
}
