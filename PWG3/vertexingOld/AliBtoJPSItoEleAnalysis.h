#ifndef AliBtoJPSItoEleAnalysis_H
#define AliBtoJPSItoEleAnalysis_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Class AliBtoJPSItoEleAnalysis
//             Reconstruction and analysis B -> J/\psi + X 
//						  |_> e+ e- 
//      
//         Origin: G.E Bruno    giuseppe.bruno@ba.infn.it                
//  based on Class for charm golden channel (D0->Kpi)
//-------------------------------------------------------------------------

#include <TString.h>
#include <TNamed.h>
#include "AliESDEvent.h"
#include "AliRun.h"

//-----------------------------------------------------------------------------
class AliBtoJPSItoEleAnalysis : public TNamed {
 public:
  //
  AliBtoJPSItoEleAnalysis();
  virtual ~AliBtoJPSItoEleAnalysis();

  void ApplySelection(const Char_t *inName="AliBtoJPSItoEle.root",
		      const Char_t *outName="AliBtoJPSItoEle_sele.root") const;
  void FindCandidates(Int_t evFirst=0,Int_t evLast=0,
		      const Char_t *outName="AliBtoJPSItoEle.root");
  void MakeTracksRefFile(AliRun *mygAlice,Int_t evFirst=0,Int_t evLast=0) const;
  void PrintStatus() const;
  void SetVertexOnTheFly() { fVertexOnTheFly=kTRUE; }
  void SetSimulation() { fSim=kTRUE; }
  void SetOnlySignal() { fOnlySignal=kTRUE; }
  void SetOnlyPrimaryJpsi() { fOnlyPrimaryJpsi=kTRUE; }
  void SetOnlySignalAndPrimaryJpsi() { fOnlySignal=kTRUE; fOnlyPrimaryJpsi=kTRUE; }
  void SetPtCut(Double_t pt=0.) { fPtCut=pt; }
  void Setd0Cut(Double_t d0=0.) { fd0Cut=d0; } 
  void SetPidCut(Double_t pid=0.) { fPidCut=pid; }
  void SetMassCut(Double_t deltaM=1000.) { fMassCut=deltaM; }
  void SetBCuts(Double_t cut0=1000.,Double_t cut1=100000.,
		 Double_t cut2=1.1,Double_t cut3=0.,Double_t cut4=0.,
		 Double_t cut5=100000.,Double_t cut6=100000.,
		 Double_t cut7=100000000.,Double_t cut8=-1.1); 
  void SetBCuts(const Double_t cuts[9]); 
  void SetPID(const Char_t * pid="TRDTPCparam") { fPID=pid; }
//  void SetKFPrimVertex() { fKFPrimVertex=kTRUE; }          //new setter
   void SetKFSecondVertex() {fKFSecondVertex=kTRUE;}         //new setter
   void UnSetKFSecondVertex() {fKFSecondVertex=kFALSE;}      //new setter
//  void SetKFTopConstr() { fKFTopConstr=kTRUE; }            //new setter
  //
 private:
  //
  Bool_t   fVertexOnTheFly; // flag for primary vertex reco on the fly
  Bool_t   fSim;            // flag for the analysis of simulated events
  Bool_t   fOnlySignal;     // write to file only signal candidates (for sim)
  Bool_t   fOnlyPrimaryJpsi;// write to file only primary Jpsi candidates (for sim)
  TString  fPID;           // PID scheme

  Double_t fV1[3]; // primary vertex position (in cm)
  Double_t fPtCut;   // minimum track pt (in GeV/c)
  Double_t fd0Cut;   // minimum track |rphi impact parameter| (in micron) 
  Double_t fMassCut; // maximum of |InvMass-M(J/Psi)| (in GeV)
  Double_t fPidCut; // min. pid probability as an electron
  Bool_t   fKFSecondVertex;  // flag for Kalmann Filter reco of secondary vertex
//  Bool_t   fKFTopConstr;     // flag for Kalmann Filter topological constraint in primary vtx reco
//  Bool_t   fKFPrimVertex;    // flag for Kalmann Filter reco of primary vertex
  Double_t fBCuts[9]; // cuts on b candidates (see SetBCuts())
                       // (to be passed to function AliBtoJPSItoEle::Select())
                       // 0 = inv. mass half width [GeV]   
                       // 1 = dca [micron]
                       // 2 = cosThetaStar 
                       // 3 = pTP [GeV/c] (positron)
                       // 4 = pTN [GeV/c] (electron)
                       // 5 = d0P [micron]   upper limit!
                       // 6 = d0N [micron]  upper limit!
                       // 7 = d0d0 [micron^2]
                       // 8 = cosThetaPoint

  //
  Double_t CalculateTOFmass(Double_t mom,Double_t length,Double_t time) const;
  Bool_t   SelectInvMass(const Double_t p[6]) const;
  void     SelectTracks(AliESDEvent *event,
			TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
			TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const;
  void     SetVertex1(Double_t x=0.,Double_t y=0.,Double_t z=0.) 
    { fV1[0]=x;fV1[1]=y;fV1[2]=z; }
  void     SimulationInfo(TTree *treeBin,TTree *treeBout) const;
  Bool_t   SingleTrkCuts(const AliESDtrack& trk, Double_t b) const;
  //
  ClassDef(AliBtoJPSItoEleAnalysis,2)  // Reconstruction of B->JPSI-> e+e- candidates class
};


#endif


