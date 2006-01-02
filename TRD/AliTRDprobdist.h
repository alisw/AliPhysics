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

#include <TNamed.h>
#include <TFile.h>
#include <TH1F.h>

const Int_t kNMom=11; 

class AliTRDprobdist : public TNamed {
  public:
    AliTRDprobdist(Char_t *responseFile="$ALICE_ROOT/TRD/TRDdEdxHistogramsV1.root");
                       // Main Constructor
    AliTRDprobdist(const AliTRDprobdist& pd);  // Copy Constructor
    virtual ~AliTRDprobdist();               // Destructor

    Double_t GetADCNorm() const {return fADCNorm;}  // Get ADC Normalization

    Double_t GetMomentum(Int_t ip) const {return fTrackMomentum[ip];}
                      // Gets the momentum for given histogram number ip
    Double_t GetMean(Int_t iType, Int_t ip) const;        
                      // Particle type is iType and histogram number is ip         
    Double_t GetNormalization(Int_t iType, Int_t ip) const;

    TH1F* GetHistogram(Int_t iType, Int_t ip) const;

    Double_t GetProbability(Int_t iType, Double_t mom, Double_t dedx) const;
                      // Gets the Probability of having dedx
    Double_t GetProbabilityT(Int_t iType, Double_t mom, Int_t timbin) const;
                      // Gets the Probability of having timbin

    Int_t GetNbins() const {return fNbins;}         // Number of Energy bins
    Double_t GetBinSize() const {return fBinSize;}  // Size of Energy bin

    Bool_t UpdateData(Char_t *responseFile) {return ReadData(responseFile);}
                      // Update the histograms from responseFile
    void SetADCNorm(Double_t norm){fADCNorm=norm;}  // Set ADC Normalization

  protected:
    Bool_t ReadData(Char_t *responseFile);      // Reads Data
    Char_t *fpartName[5];             // Names of particle species
    Double_t fADCNorm;          // Ratio of mean charge from real Det. to prob. dist.
    Int_t fNMom;                     // Number of momenta  
    Double_t fTrackMomentum[kNMom];          // Track Momentum 

    Int_t fNbins;                   // Number of Energy bins
    Double_t fBinSize;              // Size of Energy bin

    TH1F *fh1dEdxEL[kNMom];      // Prob. of dEdx  for e
    TH1F *fh1dEdxPI[kNMom];      // Prob. of dEdx  for pi
    TH1F *fh1dEdxMU[kNMom];      // Prob. of dEdx  for muon
    TH1F *fh1dEdxKA[kNMom];      // Prob. of dEdx  for Kaon
    TH1F *fh1dEdxPR[kNMom];      // Prob. of dEdx  for proton
    TH1F *fh1MaxTimBinEL[kNMom];    // Prob. of max Time Bin for e
    TH1F *fh1MaxTimBinPI[kNMom];    // Prob. of max Time Bin for pi
    ClassDef(AliTRDprobdist,1)
};


#endif

