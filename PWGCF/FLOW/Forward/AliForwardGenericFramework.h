//
// Helper class to calculate Q cumulant in forward & central regions
//
#ifndef AliForwardGenericFramework_H
#define AliForwardGenericFramework_H
/**
 * @file AliForwardGenericFramework.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "TObject.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TString.h"
#include "TNtuple.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom.h"
#include "THn.h"
#include "TString.h"
#include "TComplex.h"

#include "AliForwardSettings.h"
#include "AliForwardNUATask.h"
#include <iostream>

/**
 * Class to handle cumulant calculations.
 */
class AliForwardGenericFramework
{
public:
  /*
  * Constructor
  */
  AliForwardGenericFramework(Int_t refbins);

  /**
   * Destructor
   */
  virtual ~AliForwardGenericFramework(){}


  AliForwardSettings fSettings;

    // Utility class for filling histograms
  /**
   * Do cumulants calculations for current event with
   * centrality cent
   *
   * @param cent Event centrality
   */
  void CumulantsAccumulate(TH2D*& dNdetadphi, double cent,double vertexpos,Bool_t useFMD,Bool_t doRefFlow, Bool_t doDiffFlow);

  void saveEvent(double cent, double vertexpos,UInt_t r, Int_t ptn);

  /**
   * Constants
   */
  enum {
    ktpcOnly = 128,        // TPC only tracks
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };

  /**
   * Reset histograms
   */
  void reset();

  TH1D* fAutoDiff;//!     // Accumulated reference particles
  TH1D* fAutoRef;//!     // Accumulated reference particles
  THnD* fQvector;//!     // Accumulated reference particles
  THnD* fpvector;//!    // Accumulated differential particles
  THnD* fqvector;//!    // Accumulated differential particles

  TComplex Q(Int_t n, Int_t p, Int_t etaBin);
  TComplex p(Int_t n, Int_t p, Int_t etaBin);
  TComplex q(Int_t n, Int_t p, Int_t etaBin);

  TComplex Two(Int_t n1, Int_t n2, Int_t eta1, Int_t eta2);
  TComplex TwoDiff(Int_t n1, Int_t n2, Int_t refetabin, Int_t diffetabin);
  TComplex TwoTwoDiff(Int_t n1, Int_t n2, Int_t diffetabin1, Int_t diffetabin2);
  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4,Int_t eta1, Int_t eta2);
  TComplex FourDiff(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t refetabinA, Int_t refetabinB, Int_t diffetabin,Int_t qetabin);
  TComplex FourDiff_SC(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t diffetabinA, Int_t refetabinA, Int_t diffetabinB,Int_t refetabinB);


  void fill(THnD*& cumu, Int_t n, Int_t ptn, Double_t sample,
                        Double_t zvertex,  
                        Double_t eta,     
                        Double_t cent,    
                        Double_t value)
  { 
    if (TMath::IsNaN(value)) return; 
    if (n < 0){
      Double_t values[5] = {sample, zvertex, eta, cent};
      cumu->Fill(values,value);
      return;
    }
    else{
      Double_t values[6] = {Double_t(n-2),sample, zvertex, eta, cent};
      cumu->Fill(values,value);
      return;
    }
  }

  Double_t applyNUAcentral(Double_t eta, Double_t phi, Double_t zvertex, Double_t weight){
    if (!fSettings.use_primaries_cen) {
      Int_t nuaeta = fSettings.nuacentral->GetXaxis()->FindBin(eta);
      Int_t nuaphi = fSettings.nuacentral->GetYaxis()->FindBin(phi);
      Int_t nuavtz = fSettings.nuacentral->GetZaxis()->FindBin(zvertex);
      Double_t factor = fSettings.nuacentral->GetBinContent(nuaeta,nuaphi,nuavtz+10*fSettings.nua_runnumber);
      if (!(TMath::IsNaN(factor)) & (factor > 0.)) weight = weight*factor;
      else weight = 0.;
    }
    return weight;
  }


  Double_t applyNUAforward(TH2D*& dNdetadphi, Int_t etaBin, Int_t phiBin, Double_t eta, Double_t phi, Double_t zvertex, Double_t weight){

    if (fSettings.nua_mode & fSettings.kInterpolate) {
      weight = AliForwardNUATask::InterpolateWeight(*dNdetadphi,phiBin,etaBin,weight);
    }
    if (!fSettings.use_primaries_fwd) {
      Int_t nuaeta = fSettings.nuaforward->GetXaxis()->FindBin(eta);
      Int_t nuaphi = fSettings.nuaforward->GetYaxis()->FindBin(phi);
      Int_t nuavtz = fSettings.nuaforward->GetZaxis()->FindBin(zvertex);
      Double_t factor = fSettings.nuaforward->GetBinContent(nuaeta,nuaphi,nuavtz+10*fSettings.nua_runnumber);
      if (!(TMath::IsNaN(factor)) & (factor > 0.)) weight = weight*factor;
      else weight = 0.;
    }
    return weight;
  }

  Double_t applyInterpolateCorr(Int_t n, Double_t cent, Int_t etaBin, Double_t weight){
    Int_t ccbin = 0;

    if ((cent > 0.) & (cent < 5.)) ccbin = 0;
    if ((cent > 5.) & (cent < 10.)) ccbin = 1;
    else{
      for (Int_t cbin = 1; cbin < 7; cbin++){
        if ((cent > (cbin)*10 ) & (cent < (cbin+1)*10 )){
          ccbin = cbin+1;
          break;
        }
      }
    }
 
    if (n == 2){
      Double_t hole1_e1[] = {1.002, 1.003, 1.009, 1.017, 1.023, 1.023, 1.018, 1.023, 1.023};
      Double_t hole1_e2[] = {1.002, 1.004, 1.011, 1.023, 1.031, 1.029, 1.023, 1.018, 1.018};
      Double_t hole1_e3[] = {1.001, 1.004, 1.01 , 1.019, 1.024, 1.023, 1.021, 1.02 , 1.02 };
      Double_t hole2_e1[] = {1.003, 1.006, 1.01 , 1.014, 1.018, 1.017, 1.015, 1.011, 1.011};
      Double_t hole2_e2[] = {1.006, 1.008, 1.014, 1.019, 1.023, 1.022, 1.022, 1.011, 1.011};

      if (etaBin == 24) return weight*hole1_e1[ccbin];
      if (etaBin == 25) return weight*hole1_e2[ccbin];
      if (etaBin == 26) return weight*hole1_e3[ccbin];
      if (etaBin == 34) return weight*hole2_e1[ccbin];
      if (etaBin == 35) return weight*hole2_e2[ccbin];
    }

    if (n == 3){
      Double_t hole1_e1[] = {1.002, 1.004, 1.009, 1.023, 1.034, 1.039, 1.036, 1.036, 1.036};
      Double_t hole1_e2[] = {1.001, 1.005, 1.014, 1.031, 1.047, 1.045, 1.059, 1.059, 1.059};
      Double_t hole1_e3[] = {1.002, 1.004, 1.012, 1.024, 1.037, 1.04 , 1.051, 1.051, 1.051};
      Double_t hole2_e1[] = {1.003, 1.006, 1.01 , 1.018, 1.023, 1.015, 1.052, 1.052, 1.052};
      Double_t hole2_e2[] = {1.007, 1.008, 1.016, 1.02 , 1.028, 1.024, 1.05 , 1.05 , 1.05 };
      if (etaBin == 24) return weight*hole1_e1[ccbin];
      if (etaBin == 25) return weight*hole1_e2[ccbin];
      if (etaBin == 26) return weight*hole1_e3[ccbin];
      if (etaBin == 34) return weight*hole2_e1[ccbin];
      if (etaBin == 35) return weight*hole2_e2[ccbin];      
    }
    if (n == 4){
      Double_t hole1_e1[] = {1.001, 1.004, 1.009, 1.026, 1.054, 1.061, 1.061, 1.061, 1.061};
      Double_t hole1_e2[] = {1.001, 1.005, 1.012, 1.029, 1.068, 1.089, 1.089, 1.089, 1.089};
      Double_t hole1_e3[] = {1.002, 1.003, 1.011, 1.031, 1.054, 1.104, 1.104, 1.104, 1.104};
      Double_t hole2_e1[] = {1.01 , 1.006, 1.01 , 1.025, 1.02 , 1.058, 1.058, 1.058, 1.058};
      Double_t hole2_e2[] = {1.003, 1.01 , 1.038, 1.038, 1.015, 1.041, 1.041, 1.041, 1.041};
      if (etaBin == 24) return weight*hole1_e1[ccbin];
      if (etaBin == 25) return weight*hole1_e2[ccbin];
      if (etaBin == 26) return weight*hole1_e3[ccbin];
      if (etaBin == 34) return weight*hole2_e1[ccbin];
      if (etaBin == 35) return weight*hole2_e2[ccbin];      
    }
  return weight;
  }



  Double_t applySecondaryCorr(Int_t n, Double_t eta,Double_t zvertex,Double_t cent, Double_t weight){
    if (fSettings.seccorr_fwd){
      Int_t seceta = fSettings.seccorr_fwd->GetZaxis()->FindBin(eta);
      Int_t secvtz = fSettings.seccorr_fwd->GetYaxis()->FindBin(zvertex);
      Double_t factor = fSettings.seccorr_fwd->GetBinContent(n,secvtz,seceta);
      if (!(TMath::IsNaN(factor)) & (factor > 0.)) weight = weight*factor;
      else weight = 0.;
    }

    if (fSettings.seccorr_cent){
      Int_t seceta = fSettings.seccorr_cent->GetYaxis()->FindBin(eta);
      Int_t seccent = 1;
      if (fSettings.XeXe){
        if (cent < 10) seccent = 2;
        if (cent > 10) seccent = fSettings.seccorr_cent->GetZaxis()->FindBin(cent+10);
      }
      if (!fSettings.XeXe){
        seccent = fSettings.seccorr_cent->GetZaxis()->FindBin(cent);
      }
      Double_t factor = fSettings.seccorr_cent->GetBinContent(n,seceta,seccent);
      if (!(TMath::IsNaN(factor)) & (factor > 0.)) weight = weight*factor;
      else weight = 0.; 
    }
    
    return weight;
  }

  // reference
  THnD* cumu_rW2         ;//!
  THnD* cumu_rW2Two      ;//!
  THnD* cumu_rW4         ;//!
  THnD* cumu_rW4Four     ;//!

  // normal
  THnD* cumu_dW2B        ;//!
  THnD* cumu_dW2TwoB     ;//!
  THnD* cumu_dW4         ;//!
  THnD* cumu_dW4Four     ;//!

  // SC 
  THnD* cumu_dW2TwoTwoN  ;//!
  THnD* cumu_dW2TwoTwoD  ;//!
  THnD* cumu_dW4ThreeTwo ;//!
  THnD* cumu_dW4FourTwo  ;//!
  THnD* cumu_wSC  ;//!

  // decorrelation
  THnD* cumu_dW2A        ;//!
  THnD* cumu_dW2TwoA     ;//!
  THnD* cumu_dW22TwoTwoN ;//! // Numerator of R_{n,n; 2}
  THnD* cumu_dW22TwoTwoD ;//! // Denominator of R_{n,n; 2}



  ClassDef(AliForwardGenericFramework, 1); // object for eta dependent cumulant ananlysis
};

#endif
