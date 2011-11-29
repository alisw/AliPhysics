/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// class to extract omega(782)->pi0+gamma->3gamma
//
//-- Author: Renzhuo Wan (IOPP-Wuhan, China)
//_________________________________________________________________________
#ifndef ALIANAOMEGATOPI0GAMMA_H
#define ALIANAOMEGATOPI0GAMMA_H
//Root
class TList;
class TH2F ;
class TLorentzVector;
//Analysis
#include "AliAnaPartCorrBaseClass.h"
class TParticle;

class AliAnaOmegaToPi0Gamma : public AliAnaPartCorrBaseClass {
  
  public: 
  
  AliAnaOmegaToPi0Gamma() ; // default ctor
  AliAnaOmegaToPi0Gamma(const char *name) ; // default ctor
  virtual ~AliAnaOmegaToPi0Gamma() ;//virtual dtor
  
  TList * GetCreateOutputObjects(); 
  void Print(const Option_t * opt) const;
  
  void InitParameters();
  void MakeAnalysisFillHistograms();
  void Terminate(TList * outList);

  TString GetInputAODPhotonName()  const { return fInputAODGammaName;}
  void    SetInputAODPhotonName(TString & name) { fInputAODGammaName = name; }
  Bool_t  IsBadRun(Int_t /*iRun*/) const { return kFALSE;} //Tests if this run bad according to private list

  void SetNCentBin(Int_t nbin){fNCentBin = nbin;}
  void SetNPID(Int_t pid) {fNpid=pid;} //N pid cut 
  void SetNVtxZ(Int_t vtx){fNVtxZBin=vtx;} //N vertex Z cut
  void SetNDistToBadChannel(Int_t ndist){fNBadChDistBin = ndist;}
  void SetPi0MassPeakWidthCut(Double_t win){fPi0MassWindow=win;} 

  void SetPi0OverOmegaPtCut(Double_t cut){fPi0OverOmegaPtCut=cut;}
  void SetGammaOverOmegaPtCut(Double_t cut){fGammaOverOmegaPtCut=cut;}
  void SetEOverlapCluster(Double_t e){fEOverlapCluster=e;}
  void ReadHistograms(TList * outputList);

  private:

  TClonesArray * fInputAODPi0;   //Input AOD pi0 array
  TString fInputAODGammaName;    //Input AOD gamma name
  TList ** fEventsList;          //event list for mixing 
 
  Int_t fNVtxZBin;               //Number of vertex z cut
  Int_t fNCentBin;               //Number of centrality cut
  Int_t fNRpBin;                 //Number of reaction plane cut
  Int_t fNBadChDistBin;          //Number of bad channel dist cut
  Int_t fNpid;                   //Number of PID cut

  Double_t *fVtxZCut;            //![fNVtxZBin] vtertx z cut
  Double_t *fCent;               //![fNCentBin] centrality cut
  Double_t *fRp;                 //![fNRpBin] reaction plane cut
  
  Double_t fPi0Mass;             //nominal pi0 mass
  Double_t fPi0MassWindow;       //pi0 mass windows
  Double_t fPi0OverOmegaPtCut;   //pi0 Pt over omega pt cut
  Double_t fGammaOverOmegaPtCut; //gamma pt over omega pt cut
  Double_t fEOverlapCluster;     //the pt when the two photons overlapped

  TH2F * fhEtalon;               //an etalon of 3D histograms
  TH2F **fRealOmega0;            //real omega IVM(asy, pt, m), with Asy_pi0<1 
  TH2F **fMixAOmega0;            //mixA omega IVM(asy, pt, m) 
  TH2F **fMixBOmega0;            //mixB omega IVM(asy, pt, m) 
  TH2F **fMixCOmega0;            //mixC omega IVM(asy, pt, m) 
  TH2F **fRealOmega1;            //real omega IVM(asy, pt, m), with Asy_pi0<0.7
  TH2F **fMixAOmega1;            //mixA omega IVM(asy, pt, m)
  TH2F **fMixBOmega1;            //mixB omega IVM(asy, pt, m)
  TH2F **fMixCOmega1;            //mixC omega IVM(asy, pt, m)
  TH2F **fRealOmega2;            //real omega IVM(asy, pt, m), with Asy_pi0<0.8
  TH2F **fMixAOmega2;            //mixA omega IVM(asy, pt, m)
  TH2F **fMixBOmega2;            //mixB omega IVM(asy, pt, m)
  TH2F **fMixCOmega2;            //mixC omega IVM(asy, pt, m)

  TH2F **fhFakeOmega;            //high pt clusters assumed as pi0 + another gamma 

  TH1F *fhOmegaPriPt;            //MC primary omega pt in 2pi and |y|<0.5
  
  AliAnaOmegaToPi0Gamma(const AliAnaOmegaToPi0Gamma & ex) ;              // cpy ctor
  AliAnaOmegaToPi0Gamma & operator = (const AliAnaOmegaToPi0Gamma & ex) ;// cpy assignment
  
  ClassDef(AliAnaOmegaToPi0Gamma,2)

} ;

#endif 



