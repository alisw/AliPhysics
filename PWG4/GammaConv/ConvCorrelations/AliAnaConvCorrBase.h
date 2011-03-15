/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnaConvCorrPhoton.h
/// @author Svein Lindal
/// @brief  Base class for analysis of correlations between conversion particles and charged tracks

#ifndef ALIANACONVCORRBASE_CXX
#define ALIANACONVCORRBASE_CXX

#include "Rtypes.h"
#include "TMath.h"
#include "TList.h"
#include "TH1.h"

class TH1F;
class TH1I;
class TH2F;
class AliAODConversionParticle;
class TClonesArray;
class TNtuple;
class TString;


class AliAnaConvCorrBase : public TObject {

public:

  

  //Constructor / desctructor
  AliAnaConvCorrBase(TString name); 
  virtual ~AliAnaConvCorrBase();
  
  //Set and get min pt for triggers
  void SetTriggerPt(Float_t pt) { fTriggerPt = pt; }
  inline Float_t GetTriggerPt() const {return fTriggerPt; }


  //Set and get min pt for correlation particles
  void SetCorrelatedPt(Float_t pt) { fCorrelatedPt = pt; }
  inline Float_t GetCorrelatedPt() const {return fCorrelatedPt; }

  //CreateHistograms
  void CreateHistograms();
  
  //Get list of histograms
  TList * GetHistograms() const { return fHistograms;}

  //Add histogram to list
  void AddHistogram(TH1 * histogram) { fHistograms->Add(dynamic_cast<TObject*>(histogram));}

  //Set and get number of bins in phi direction
  void SetNPhiBins(Int_t bins) { fNPhiBins = bins; }
  Int_t GetNPhiBins() const { return fNPhiBins;}

  ///Get the distance in phi between trigger particle and correlated particle
  Float_t GetDPhi(Float_t dPhi) { 
    if ( dPhi < 3*TMath::PiOver2() && dPhi > - TMath::PiOver2() ) return dPhi;
    else return ( (dPhi>0)? dPhi - TMath::TwoPi() : dPhi + TMath::TwoPi() ); 
  }

  ///Get bin limits of various pt bins
  Float_t GetLowerBinLimit(const Int_t bin, const Float_t * const bins) const;
  Float_t GetUpperBinLimit(const Int_t bin, const Float_t * const bins) const;

  //Get trigger bin for particle carrying given pt
  Int_t GetTriggerBin(Float_t pt) const;
  //Get correlation particle pt bin for particle carrying given pt
  Int_t GetCorrBin(Float_t pt) const;


  //Print statistics for histograms
  void PrintStatistics();


protected:

  //Fill histograms
  void FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated);
  //Fill trigger counter histograms
  void FillTriggerCounters(Float_t tPt, Bool_t isolated);

 private:
  
  TString fName; //name of analysis

  static const Int_t fNTBins = 3; //Number of trigger pt bins
  static const Int_t fNCBins = 3; //Number of corr particles pt bins

  TList * fHistograms; //List of histograms

  TH2F * fHEtaPhiPt[2]; //2D eta phi correlations histogram
  TH1F * fHdEta[2]; //Eta correlations histogram
  TH1F * fHdPhi[2]; //Phi correlations histogram

  TH1F * fHdPhiBins[2][fNTBins][fNCBins]; //dPhi correlations histograms in several bins
  TH1I * fHNTriggers[2]; //Histograms containing number of triggers in various bins
  
  Float_t fTBins[fNTBins]; ///Array of trigger bin limits
  Float_t fCBins[fNCBins];///Array of corr particle pt bin limits
  
  Float_t fTriggerPt; //Min trigger particle pt
  Float_t fCorrelatedPt; // Min correlation particle pt

  Int_t fNPhiBins; //Number of bins in phi 
  
  //Default constructor prohibited
  AliAnaConvCorrBase(); //not implemented
  AliAnaConvCorrBase(const AliAnaConvCorrBase&); // not implemented
  AliAnaConvCorrBase& operator=(const AliAnaConvCorrBase&); // not implemented
  ClassDef(AliAnaConvCorrBase, 2); // example of analysis
};

#endif
