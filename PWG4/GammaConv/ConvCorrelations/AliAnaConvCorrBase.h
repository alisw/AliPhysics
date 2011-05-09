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
class TH3F;
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
  // void SetTriggerPt(Float_t pt) { fTriggerPt = pt; }
  // inline Float_t GetTriggerPt() const {return fTriggerPt; }


  // //Set and get min pt for correlation particles
  // void SetCorrelatedPt(Float_t pt) { fCorrelatedPt = pt; }
  // inline Float_t GetCorrelatedPt() const {return fCorrelatedPt; }

  //CreateHistograms
  void CreateBaseHistograms();
  //To be overrriden by children. Should always call CreateBaseHistograms()
  virtual void CreateHistograms();
  
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


  TArrayD * GetTriggerBins() const { return fPtBins; }

  //Print statistics for histograms
  void PrintStatistics();


protected:

  //Fill histograms
  void FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated);
  //Fill trigger counter histograms
  void FillTriggerCounters(Float_t tPt, Bool_t isolated);

 private:
  
  TString fName; //name of analysis
  TList * fHistograms; //List of histograms

  Int_t fNPhiBins;  //Nbins in phi direction
  TArrayD * fdPhiBins; //!transient phi bin array
  TArrayD * fPtBins; //!Array of trigger bins
  

  TH3F * fHdPhi[2]; //dPhi pt histogram
  TH1F * fHNTriggers[2]; //Histograms containing number of triggers in various bins


  //Default constructor prohibited
  AliAnaConvCorrBase(); //not implemented
  AliAnaConvCorrBase(const AliAnaConvCorrBase&); // not implemented
  AliAnaConvCorrBase& operator=(const AliAnaConvCorrBase&); // not implemented
  ClassDef(AliAnaConvCorrBase, 3); // example of analysis
};

#endif
