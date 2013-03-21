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
#include <THnSparse.h>

class TH1F;
class TH3F;
class TH2F;
class AliAODConversionParticle;
class TClonesArray;
class TNtuple;
class TString;

class AliAnaConvCorrBase : public TNamed {

public:

  

  //Constructor / desctructor
  AliAnaConvCorrBase(TString name, TString title); 
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

  void AddAxis(TAxis * axis) { fAxesList.Add(axis); }

  ///Get the distance in phi between trigger particle and correlated particle
  Float_t GetDPhi(Float_t dPhi) { 
    if ( dPhi < 3*TMath::PiOver2() && dPhi > - TMath::PiOver2() ) return dPhi;
    else return ( (dPhi>0)? dPhi - TMath::TwoPi() : dPhi + TMath::TwoPi() ); 
  }

  void PrintStatistics();

  void CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray * tracks, const Int_t tIDs[4], Int_t isolated);
  //virtual void FillTriggerCounters(const AliAODConversionParticle * particle, Int_t leading);

  TAxis& GetAxistPt()       { return fAxistPt;   }
  TAxis& GetAxiscPt()       { return fAxiscPt;   }
  TAxis& GetAxisdEta()      { return fAxisdEta;  }
  TAxis& GetAxisdPhi()      { return fAxisdPhi;  }
  TAxis& GetAxisIso()       { return fAxisIso;   }

  TAxis& GetAxisMEEta()     { return fAxisMEEta; }
  TAxis& GetAxisMEPhi()     { return fAxisMEPhi; }

  TList& GetAxisList()      { return fAxesList;  }
  TList& GetTrackAxisList() { return fTrackAxisList; }
  TList& GetTrigAxisList()  { return fTrigAxisList; }


protected:

  //Fill histograms
  //void FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated);
  //Fill trigger counter histograms
  //virtual void FillTriggerCounters(Float_t tPt, Bool_t isolated) = NULL;
  THnSparseF * CreateSparse(TString name, TString title, TList * axes);

private:

  void SetUpDefaultBins();
  
  //TString fName; //name of analysis
  TList * fHistograms; //List of histograms
  TList fAxesList;  //List over axes to be used in sparse
  TList fTrigAxisList; //list
  TList fTrackAxisList; //list

  TAxis fAxistPt; //Pt axis
  TAxis fAxiscPt; //correlated particle pt axis
  TAxis fAxisdEta; //delta eta axis
  TAxis fAxisdPhi; //delta phi axis
  TAxis fAxisIso; //Isolated particle axis
  
  TAxis fAxisMEEta ; //Eta axis for ME
  TAxis fAxisMEPhi ; //Phi axis for ME

  THnSparseF * fCorrSparse; // Sparse for corr
  THnSparseF * fTrigSparse; // ME Sparse
  THnSparseF * fTrackSparse; //Track Sparse

  //Default constructor prohibited
  AliAnaConvCorrBase(); //not implemented
  AliAnaConvCorrBase(const AliAnaConvCorrBase&); // not implemented
  AliAnaConvCorrBase& operator=(const AliAnaConvCorrBase&); // not implemented

  ClassDef(AliAnaConvCorrBase, 5); // example of analysis

};

#endif
