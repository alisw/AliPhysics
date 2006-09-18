#ifndef AliD0toKpiAnalysis_H
#define AliD0toKpiAnalysis_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Class AliD0toKpiAnalysis
//             Reconstruction and analysis D0 -> K^- pi^+
//      
//         Origin: A. Dainese    andrea.dainese@pd.infn.it                  
//-------------------------------------------------------------------------

#include <TString.h>
#include <TNamed.h>
#include "AliESD.h"

//-----------------------------------------------------------------------------
class AliD0toKpiAnalysis : public TNamed {
 public:
  //
  AliD0toKpiAnalysis();
  virtual ~AliD0toKpiAnalysis();

  void ApplySelection(const Char_t *inName="AliD0toKpi.root",
		      const Char_t *outName="AliD0toKpi_sele.root") const;
  //  void FindCandidates(Int_t evFirst=0,Int_t evLast=0,
  //		      const Char_t *outName="AliD0toKpi.root");
  void FindCandidatesESD(Int_t evFirst=0,Int_t evLast=0,
  			 const Char_t *outName="AliD0toKpi.root");
  void PrintStatus() const;
  void SetVertexOnTheFly() { fVertexOnTheFly=kTRUE; }
  void SetSimulation() { fSim=kTRUE; }
  void SetOnlySignal() { fOnlySignal=kTRUE; }
  void SetPtCut(Double_t pt=0.) { fPtCut=pt; }
  void Setd0Cut(Double_t d0=0.) { fd0Cut=d0; } 
  void SetMassCut(Double_t deltaM=1000.) { fMassCut=deltaM; }
  void SetD0Cuts(Double_t cut0=1000.,Double_t cut1=100000.,
		 Double_t cut2=1.1,Double_t cut3=0.,Double_t cut4=0.,
		 Double_t cut5=100000.,Double_t cut6=100000.,
		 Double_t cut7=100000000.,Double_t cut8=-1.1); 
  void SetD0Cuts(const Double_t cuts[9]); 
  void SetPID(const Char_t * pid="TOFparam_PbPb") { fPID=pid; }
  //
 private:
  //
  Bool_t   fVertexOnTheFly; // flag for primary vertex reco on the fly
  Bool_t   fSim;            // flag for the analysis of simulated events
  Bool_t   fOnlySignal;     // write to file only signal candidates (for sim)
  TString  fPID;           // PID scheme

  Double_t fV1[3]; // primary vertex position (in cm)
  Double_t fPtCut;   // minimum track pt (in GeV/c)
  Double_t fd0Cut;   // minimum track |rphi impact parameter| (in micron) 
  Double_t fMassCut; // maximum of |InvMass-MD0| (in GeV)
  Double_t fD0Cuts[9]; // cuts on D0 candidates (see SetD0Cuts())
                       // (to be passed to function AliD0toKpi::Select())
                       // 0 = inv. mass half width [GeV]   
                       // 1 = dca [micron]
                       // 2 = cosThetaStar 
                       // 3 = pTK [GeV/c]
                       // 4 = pTPi [GeV/c]
                       // 5 = d0K [micron]   upper limit!
                       // 6 = d0Pi [micron]  upper limit!
                       // 7 = d0d0 [micron^2]
                       // 8 = cosThetaPoint

  //
  Double_t CalculateTOFmass(Double_t mom,Double_t length,Double_t time) const;
  //void     MakeTracksRefFile(Int_t evFirst=0,Int_t evLast=0) const;
  void     MakeTracksRefFileESD() const;
  Bool_t   SelectInvMass(const Double_t p[6]) const;
  //void     SelectTracks(TTree &trkTree,
  //			TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
  //			TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const;
  void     SelectTracksESD(AliESD &event,
			   TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
			   TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const;
  void     SelectTracksESDvtx(AliESD &event,TTree *trkTree,
			      TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
			      TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const;
  void     SetVertex1(Double_t x=0.,Double_t y=0.,Double_t z=0.) 
    { fV1[0]=x;fV1[1]=y;fV1[2]=z; }
  void     SimulationInfo(TTree *treeD0in,TTree *treeD0out) const;
  Bool_t   SingleTrkCuts(const AliESDtrack& trk, Double_t b) const;
  //
  ClassDef(AliD0toKpiAnalysis,2)  // Reconstruction of D0 candidates class
};


#endif








