#ifndef ALITRDPROBDIST_H
#define ALITRDPROBDIST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*-----------------------------------------------------------------
   Class for dE/dx and Time Bin of Max. Cluster for Electrons and 
   pions in TRD. 
   It is instantiated in class AliTRDpidESD for particle identification
   in TRD
   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>
   -----------------------------------------------------------------*/

#include "TObject.h"

const Int_t kNo_Mom=7; 
const Int_t kNo_EnBins=250;
const Int_t kNo_TBins=20;

class AliTRDprobdist : public TObject {
  public:
    AliTRDprobdist(Int_t multiplicity=1);   // multiplicity can take 
                         // values 1, 2000, 4000, 6000, 8000
    //    ~AliTRDprobdist(); 
    void FillData();      // Fills Data 
    //    void FillData2000();  // Fills Data of multiplicity 2000
    //    void FillData4000();  // Fills Data of multiplicity 4000
    //    void FillData6000();  // Fills Data of multiplicity 6000
    //    void FillData8000();  // Fills Data of multiplicity 8000
    Double_t GetBM(Int_t ip) const { return fTrackMomentum[ip];}
    Double_t GetMeanPI(Int_t ip) const;  // Gets mean of de/dx dist. of pi
    Double_t GetMeanEL(Int_t ip) const;  // Gets mean of de/dx dist. of e
    Double_t GetNormalizationPI(Int_t ip) const;  // Gets Norm. of de/dx dist. of pi
    Double_t GetNormalizationEL(Int_t ip) const;  // Gets Norm. of de/dx dist. of e
    Double_t GetProbability(Int_t k, Double_t mom, Double_t dedx) const;
                          // Gets the Probability of having dedx
    Double_t GetProbabilityT(Int_t k, Double_t mom, Int_t timbin) const;
                        // Gets the Probability of having timbin
    void GetData(Int_t ip, Double_t *ebin, Double_t *ppi, Double_t *pel) const { 
      for(Int_t ie=0; ie<fNEbins; ie++){
        ebin[ie]=fEnergyLoss[ip][ie];
        ppi[ie]=fProbPI[ip][ie];
        pel[ie]=fProbEL[ip][ie];
      }
    }
    void GetDataT(Int_t ip, Int_t *tbin, Double_t *ppi, Double_t *pel) const { 
      for(Int_t it=0; it<fNTbins; it++){
        tbin[it]=fTimBin[ip][it];
        ppi[it]=fProbPIT[ip][it];
        pel[it]=fProbELT[ip][it];
      }
    }
  protected:
    Double_t fADCNorm;               // Ratio of de/dx from Det. to prob. dist.
    Int_t fNMom;                     // Number of momenta  
    Double_t fTrackMomentum[kNo_Mom];           // Track Momentum 
    Int_t fNEbins;                              // Number of Energy bins
    Double_t fEnBinSize;                        // Size of Energy bin
    Double_t fEnergyLoss[kNo_Mom][kNo_EnBins];  // dE/dx 
    Double_t fProbPI[kNo_Mom][kNo_EnBins];      // Prob. of dEdx  for pi
    Double_t fProbEL[kNo_Mom][kNo_EnBins];      // Prob. of dEdx  for e

    Int_t fNTbins;                             // Number of Tim bins=20      
    Double_t fTBinSize;                        // Size of Time Bin =1 
    Int_t fTimBin[kNo_Mom][kNo_TBins];      // Time Bin  
    Double_t fProbPIT[kNo_Mom][kNo_TBins];      // Prob. of dEdx for pi
    Double_t fProbELT[kNo_Mom][kNo_TBins];    // Prob. of dEdx for e
    ClassDef(AliTRDprobdist,1)
};


#endif

