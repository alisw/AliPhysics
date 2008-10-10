#include "TObject.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "AliLog.h"
#include "Cal/AliTRDCalPID.h"
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
void  AliTRDpidUtil::CalculatePionEffi(TH1* histo1, TH1* histo2)
// Double_t  AliTRDpidUtil::GetError()
{
  //
  // Calculates the pion efficiency
  //

  histo1 -> SetLineColor(kRed);
  histo2 -> SetLineColor(kBlue); 
  AliInfo(Form("Histo1[%d] Histo2[%d]", (Int_t)histo1 -> GetEntries(), (Int_t)histo2 -> GetEntries()));
  if(!(histo1 -> GetEntries() > 0 && histo2 -> GetEntries() > 0)){
    AliWarning("Histo has no Entries !");
    return;
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
  TGraph *gEffis = new TGraph(kBins, SumElecsE, SumPionsE);

  // use fit function to get derivate of the TGraph for error calculation
  TF1 *f1 = new TF1("f1","[0]*x*x+[1]*x+[2]", fEleEffi-.05, fEleEffi+.05);
  gEffis -> Fit("f1","Q","",fEleEffi-.05, fEleEffi+.05);
  
  // return the error of the pion efficiency
  if(((1.-fPionEffi) < 0) || ((1.-fCalcEleEffi) < 0)){
    AliWarning(" EleEffi or PioEffi > 1. Error can not be calculated. Please increase statistics or check your simulation!");
    return;
  }
  fError = sqrt(((1/histo2 -> GetEntries())*fPionEffi*(1-fPionEffi))+((f1 -> Derivative(fEleEffi))*(f1 -> Derivative(fEleEffi))*(1/histo1 -> GetEntries())*fCalcEleEffi*(1-fCalcEleEffi)));

  AliInfo(Form("Pion Effi at [%f] : [%f +/- %f], Threshold[%f]", fCalcEleEffi, fPionEffi, fError, fThreshold));
  AliInfo(Form("Derivative at %4.2f : %f\n", fEleEffi, f1 -> Derivative(fEleEffi)));

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

