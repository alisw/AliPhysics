#ifndef ALITPCTRACKERPARAM_H
#define ALITPCTRACKERPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice                               */

/* $Id$ */


//-----------------------------------------------------------------------------
//                    TPC Tracking Parameterization Class
//
//   Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it
//-----------------------------------------------------------------------------

//----- Root headers ---------
#include <TMatrixD.h>
//---- AliRoot headers -------
#include "alles.h"
#include "AliGausCorr.h"
#include "AliMagF.h"
#include "AliTPCkineGrid.h"
#include "AliTPCtrack.h"
#include "AliTrackReference.h"
//----------------------------

class AliTPCtrackerParam {
  /////////////////////////////////////////////////////////////////////////
  //                                                                        
  // This class builds AliTPCtrack objects from generated tracks to feed    
  // ITS tracking (V2). The AliTPCtrack is built from its first hit in      
  // the TPC. The track is assigned a Kalman-like covariance matrix         
  // depending on its pT and pseudorapidity and track parameters are        
  // smeared according to this covariance matrix.                           
  // Output file contains sorted tracks, ready for matching with ITS.        
  //                                                                        
  // See implementation file for more details.  
  //                                  
  //                                                                        
  //  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it     
  //                                                                        
  /////////////////////////////////////////////////////////////////////////
 public:
  AliTPCtrackerParam(const Int_t coll=0,const Double_t Bz=0.4,const Int_t n=1);
  virtual ~AliTPCtrackerParam();

  // this function performs the parameterized tracking
  Int_t BuildTPCtracks(const TFile *inp, TFile *out);

  // these functions are used to create a DB of cov. matrices,
  // including regularization, efficiencies and dE/dx
  void  AllGeantTracks() { fSelAndSmear=kFALSE; return; }
  void  AnalyzedEdx(const Char_t *outName,Int_t pdg);
  void  AnalyzePulls(const Char_t *outName);
  void  AnalyzeResolutions(Int_t pdg);
  void  CompareTPCtracks(const Char_t *galiceName="galice.root",
			 const Char_t *trkGeaName="AliTPCtracksGeant.root",
			 const Char_t *trkKalName="AliTPCtracksSorted.root",
			 const Char_t *covmatName="CovMatrix.root",
			 const Char_t *tpceffasciiName="TPCeff.dat",
			 const Char_t *tpceffrootName="TPCeff.root");
  void  DrawEffs(const Char_t *inName,Int_t pdg=211);
  void  DrawPulls(const Char_t *inName,Int_t pdg=211,Int_t par=0);
  void  MakeDataBase();
  void  MergeEvents(Int_t evFirst=1,Int_t evLast=1);
  void  RegularizeCovMatrix(const Char_t *outName,Int_t pdg);


  //********* Internal class definition *******
  class AliTPCtrackParam : public AliTPCtrack {
  public:
    AliTPCtrackParam():AliTPCtrack(){}
    AliTPCtrackParam(const AliTPCtrack &t):AliTPCtrack(t){}

    void AssignMass(Double_t mass) {SetMass(mass); return;}
    
  private:

  };
  //********* end of internal class ***********

  //********* Internal class definition *******
  class AliTPCseedGeant : public TObject {
  public:
    AliTPCseedGeant(Double_t x,Double_t y,Double_t z,
		    Double_t px,Double_t py,Double_t pz,
		    Int_t lab) {
      fXg = x;
      fYg = y;
      fZg = z;
      fPx = px;
      fPy = py;
      fPz = pz;
      fLabel = lab;
      Double_t a = TMath::ATan2(y,x)*180./TMath::Pi();
      if(a<0) a += 360.;
      fSector = (Int_t)(a/20.);
      fAlpha = 10.+20.*fSector;
      fAlpha /= 180.;
      fAlpha *= TMath::Pi();
    }
    Int_t    GetLabel() { return fLabel; }
    Double_t GetAlpha() { return fAlpha; }      
    Double_t GetXL() { return fXg*TMath::Cos(fAlpha)+fYg*TMath::Sin(fAlpha); }
    Double_t GetYL() { return -fXg*TMath::Sin(fAlpha)+fYg*TMath::Cos(fAlpha); }
    Double_t GetZL() { return fZg; }
    Double_t GetPx() { return fPx; } 
    Double_t GetPy() { return fPy; } 
    Double_t GetPz() { return fPz; } 
    Double_t GetPt() { return TMath::Sqrt(fPx*fPx+fPy*fPy); }
    Double_t GetEta() { return -TMath::Log(TMath::Tan(0.25*TMath::Pi()-0.5*TMath::ATan(fPz/GetPt()))); }
    void     SetLabel(Int_t lab) { fLabel=lab; return; }
    Bool_t   InTPCAcceptance() {
      if(TMath::Abs(GetZL()+(244.-GetXL())*fPz/GetPt())>252.) return kFALSE;
      return kTRUE;
    }

  private:
    Double_t fXg;
    Double_t fYg;
    Double_t fZg;
    Double_t fPx;
    Double_t fPy;
    Double_t fPz;
    Double_t fAlpha;
    Int_t    fLabel;
    Int_t    fSector;
  };
  //******* end of internal class ****************
  
 private:
  Int_t           fNevents;     // number of events in the file to be processed
  Double_t        fBz;          // value of the z component of L3 field (Tesla)
  Int_t           fColl;        // collision code (0: PbPb6000; 1: pp)
  Bool_t          fSelAndSmear; // if kFALSE returns GEANT tracks 
                                // at TPC first hit 
  TString         fDBfileName;  // DataBase file name
  
  AliTPCtrack     fTrack;    // current track

  TTree          *fCovTree;  // tree with regularized cov matrices 
                             // for the current track

  AliTPCkineGrid *fDBgrid;   // grid for the cov matrix look-up table  
  AliTPCkineGrid  fDBgridPi; //               "                  for pions  
  AliTPCkineGrid  fDBgridKa; //               "                  for kaons
  AliTPCkineGrid  fDBgridPr; //               "                  for protons
  AliTPCkineGrid  fDBgridEl; //               "                  for electrons
  AliTPCkineGrid  fDBgridMu; //               "                  for muons

  AliTPCkineGrid *fEff;   // TPC efficiencies for the current track
  AliTPCkineGrid  fEffPi; //           "        pions 
  AliTPCkineGrid  fEffKa; //           "        kaons 
  AliTPCkineGrid  fEffPr; //           "        protons 
  AliTPCkineGrid  fEffEl; //           "        electrons 
  AliTPCkineGrid  fEffMu; //           "        muons 

  AliTPCkineGrid *fPulls;      // pulls for the current track
  AliTPCkineGrid  fPullsPi[5]; //        "      pions
  AliTPCkineGrid  fPullsKa[5]; //        "      muons
  AliTPCkineGrid  fPullsPr[5]; //        "      protons
  AliTPCkineGrid  fPullsEl[5]; //        "      electrons
  AliTPCkineGrid  fPullsMu[5]; //        "      muons

  TMatrixD       *fRegPar;     // regularization parameters for the curr. track
  TMatrixD        fRegParPi;   //                  "        for pions         
  TMatrixD        fRegParKa;   //                  "        for kaons
  TMatrixD        fRegParPr;   //                  "        for protons
  TMatrixD        fRegParEl;   //                  "        for electrons
  TMatrixD        fRegParMu;   //                  "        for muons

  AliTPCkineGrid *fdEdxMean;   // dEdx mean for the current track
  AliTPCkineGrid  fdEdxMeanPi; //                "     pions
  AliTPCkineGrid  fdEdxMeanKa; //                "     kaons
  AliTPCkineGrid  fdEdxMeanPr; //                "     protons    
  AliTPCkineGrid  fdEdxMeanEl; //                "     electrons
  AliTPCkineGrid  fdEdxMeanMu; //                "     muons

  AliTPCkineGrid *fdEdxRMS;    // dEdx RMS for the current track
  AliTPCkineGrid  fdEdxRMSPi;  //                "     pions
  AliTPCkineGrid  fdEdxRMSKa;  //                "     kaons
  AliTPCkineGrid  fdEdxRMSPr;  //                "     protons    
  AliTPCkineGrid  fdEdxRMSEl;  //                "     electrons
  AliTPCkineGrid  fdEdxRMSMu;  //                "     muons


  void     BuildTrack(AliTPCseedGeant *s,Int_t ch);
  Int_t    CheckLabel(AliTPCseedGeant *s,Int_t nPart,
		      Double_t *ptkine,Double_t *pzkine) const;
  void     CookdEdx(Double_t pt,Double_t eta);
  void     CookTrack(Double_t pt,Double_t eta);
  TMatrixD GetSmearingMatrix(Double_t *cc, Double_t pt,Double_t eta) const;
  void     InitializeKineGrid(Option_t *which);    
  void     MakeSeedsFromHits(AliTPC *TPC,TTree *TH,TObjArray &seedArray) const;
  void     MakeSeedsFromRefs(TTree *TTR,
			     TObjArray &seedArray) const;
  Int_t    ReadAllData(const Char_t *inName);
  Int_t    ReadDBgrid(const Char_t *inName);
  Int_t    ReaddEdx(const Char_t *inName,Int_t pdg);
  Int_t    ReadEffs(const Char_t *inName);
  Int_t    ReadPulls(const Char_t *inName);
  Int_t    ReadRegParams(const Char_t *inName,Int_t pdg);
  Bool_t   SelectedTrack(Double_t pt, Double_t eta) const;
  void     SetParticle(Int_t pdg);  
  void     SmearTrack(Double_t *xx,Double_t *xxsm,TMatrixD cov) const;
  Int_t    WritedEdx(const Char_t *outName,Int_t pdg);
  Int_t    WriteEffs(const Char_t *outName);
  Int_t    WritePulls(const Char_t *outName);
  Int_t    WriteRegParams(const Char_t *outName,Int_t pdg);
  
  
  ClassDef(AliTPCtrackerParam,1) // TPC tracking parameterization class
};

#endif













