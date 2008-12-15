#include "TObject.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "AliLog.h"
#include "Cal/AliTRDCalPID.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliTRDpidUtil.h"

ClassImp(AliTRDpidUtil)


Float_t AliTRDpidUtil::fEleEffi = 0.9;

//________________________________________________________________________
AliTRDpidUtil::AliTRDpidUtil() 
  : TObject()
  ,fCalcEleEffi(0.)
  ,fPionEffi(-1.)
  ,fError(-1.)
  ,fThreshold(-1.)
{
  //
  // Default constructor
  //
}



//________________________________________________________________________
Bool_t  AliTRDpidUtil::CalculatePionEffi(TH1* histo1, TH1* histo2)
// Double_t  AliTRDpidUtil::GetError()
{
  //
  // Calculates the pion efficiency
  //

  histo1 -> SetLineColor(kRed);
  histo2 -> SetLineColor(kBlue); 
  if(!histo1 -> GetEntries() || !histo2 -> GetEntries()){
    AliWarning("Histo with no entries !");
    return kFALSE;
  }
  histo1 -> Scale(1./histo1->GetEntries());
  histo2 -> Scale(1./histo2->GetEntries());

  Int_t abinE, bbinE, cbinE = -1;                    
  Double_t aBinSumE, bBinSumE;                  // content of a single bin
  Bool_t bFirst = 1;                            // checks if threshold is crossed for the first time
  Double_t SumElecsE[kBins+2], SumPionsE[kBins+2];  // array of the integrated sum in each bin
  memset(SumElecsE, 0, (kBins+2)*sizeof(Double_t));
  memset(SumPionsE, 0, (kBins+2)*sizeof(Double_t));


  // calculate electron efficiency of each bin
  for (abinE = (histo1 -> GetNbinsX()); abinE >= 0; abinE--){  
   aBinSumE = 0;
    aBinSumE = histo1 -> GetBinContent(abinE);
    
    SumElecsE[abinE] = SumElecsE[abinE+1] + aBinSumE;

    if((SumElecsE[abinE] >= fEleEffi) && (bFirst == 1)){
      bFirst = 0;
      cbinE = abinE;
      fCalcEleEffi = (SumElecsE[cbinE]); 
    }
  }
  
  fThreshold = histo1 -> GetBinCenter(cbinE);

  // calculate pion efficiency of each bin
  for (bbinE = (histo2 -> GetNbinsX()); bbinE >= abinE; bbinE--){	
    bBinSumE = 0;
    bBinSumE = histo2 -> GetBinContent(bbinE);

    SumPionsE[bbinE] = SumPionsE[bbinE+1] + bBinSumE;
    if(bbinE == cbinE){
      fPionEffi = (SumPionsE[cbinE]);
    }
  }
  

  // pion efficiency vs electron efficiency
  TGraph gEffis(kBins, SumElecsE, SumPionsE);

  // use fit function to get derivate of the TGraph for error calculation
  TF1 f1("f1","[0]*x*x+[1]*x+[2]", fEleEffi-.05, fEleEffi+.05);
  gEffis.Fit(&f1, "Q", "", fEleEffi-.05, fEleEffi+.05);
  
  // return the error of the pion efficiency
  if(((1.-fPionEffi) < 0) || ((1.-fCalcEleEffi) < 0)){
    AliWarning(" EleEffi or PioEffi > 1. Error can not be calculated. Please increase statistics or check your simulation!");
    return kFALSE;
  }
  fError = sqrt(((1/histo2 -> GetEntries())*fPionEffi*(1-fPionEffi))+((f1.Derivative(fEleEffi))*(f1.Derivative(fEleEffi))*(1/histo1 -> GetEntries())*fCalcEleEffi*(1-fCalcEleEffi)));

//   AliInfo(Form("Pion Effi at [%f] : [%f +/- %f], Threshold[%f]", fCalcEleEffi, fPionEffi, fError, fThreshold));
//   AliInfo(Form("Derivative at %4.2f : %f\n", fEleEffi, f1.Derivative(fEleEffi)));
  return kTRUE;
}


//__________________________________________________________________________
Int_t AliTRDpidUtil::GetMomentumBin(Double_t p)
{
  //
  // returns the momentum bin for a given momentum
  //

  Int_t pBin1 = -1;                                    // check bin1
  Int_t pBin2 = -1;                                    // check bin2

  if(p < 0) return -1;                                 // return -1 if momentum < 0
  if(p < AliTRDCalPID::GetMomentum(0)) return 0;                                      // smallest momentum bin
  if(p >= AliTRDCalPID::GetMomentum(AliTRDCalPID::kNMom-1)) return AliTRDCalPID::kNMom-1; // largest momentum bin


  // calculate momentum bin for non extremal momenta
  for(Int_t iMomBin = 1; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
    if(p < AliTRDCalPID::GetMomentum(iMomBin)){
      pBin1 = iMomBin - 1;
      pBin2 = iMomBin;
    }
    else
      continue;

    if(p - AliTRDCalPID::GetMomentum(pBin1) >= AliTRDCalPID::GetMomentum(pBin2) - p){
       return pBin2;
    }
    else{
      return pBin1;
    }
  }

  return -1;
}


//__________________________________________________________________________
Bool_t AliTRDpidUtil::IsElectron(const AliESDtrack *track, PIDmethod_t method){
  //
  // Do PID decision for the TRD based on 90% Electron efficiency threshold
  //
  // Author: Markus Fasel (M.Fasel@gsi.de)
  //
  if(method == kESD) method = kNN;
  TString histname[2] = {"fHistThreshLQ", "fHistThreshNN"};
  AliCDBManager *cdb = AliCDBManager::Instance(); 
  AliCDBEntry *cdb_thresholds = cdb->Get("TRD/Calib/PIDThresholds");
  TObjArray *histos = dynamic_cast<TObjArray *>(cdb_thresholds->GetObject());
  TH1 * threshold_hist = dynamic_cast<TH1F *>(histos->FindObject(histname[method].Data()));
  Double_t threshold = threshold_hist->GetBinContent(GetMomentumBin(track->P()) + 1);
  
  // Do Decision
  Double_t pid_probs[5];
  track->GetTRDpid(pid_probs);
  if(pid_probs[AliPID::kElectron] >= threshold) return kTRUE;
  return kFALSE; 
}

//__________________________________________________________________________
Double_t AliTRDpidUtil::GetSystematicError(const AliESDtrack *track, PIDmethod_t method){
  //
  // Returns the pion efficiency at 90% electron efficiency from the OCDB
  //
  // Author: Markus Fasel (M.Fasel@gsi.de)
  //
  if(method == kESD) method = kNN;
  TString histname[2] = {"fHistPionEffLQ", "fHistPionEffNN"};
  AliCDBManager *cdb = AliCDBManager::Instance(); 
  AliCDBEntry *cdb_thresholds = cdb->Get("TRD/Calib/PIDThresholds");
  TObjArray *histos = dynamic_cast<TObjArray *>(cdb_thresholds->GetObject());
  TH1 * threshold_hist = dynamic_cast<TH1F *>(histos->FindObject(histname[method].Data()));
  return threshold_hist->GetBinContent(GetMomentumBin(track->P()) + 1);
}
