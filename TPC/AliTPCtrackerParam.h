#ifndef ALITPCTRACKERPARAM_H
#define ALITPCTRACKERPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice                               */

/* $Id$ */


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
//  adapted to ESD output by: Marcello Lunardon, Padova
/////////////////////////////////////////////////////////////////////////


//----- Root headers ---------
class TTree;
#include <TMatrixD.h>
//---- AliRoot headers -------
#include "AliConfig.h"
#include "AliTPCkineGrid.h"
#include "AliTPCtrack.h"
#include "AliESD.h"
//----------------------------

class AliTPC;

class AliTPCtrackerParam:
  public TObject
{   
 public:
  
  AliTPCtrackerParam(Int_t coll=0, Double_t Bz=0.4,
		     const char* evfoldname = AliConfig::GetDefaultEventFolderName());
  virtual ~AliTPCtrackerParam();

  // this function performs the parameterized tracking
  //
  AliTPCtrackerParam(const AliTPCtrackerParam& p);
  //
  Int_t Init();
  Int_t BuildTPCtracks(AliESD* event);
  // these functions are used to create a DB of cov. matrices,
  // including regularization, efficiencies and dE/dx
  void  AllGeantTracks() { fSelAndSmear=kFALSE; return; }
  void  AnalyzedEdx(const Char_t *outName,Int_t pdg);
  void  AnalyzePulls(const Char_t *outName);
  void  AnalyzeResolutions(Int_t pdg);
  void  CompareTPCtracks(Int_t nEvents=1,
			 const Char_t *galiceName="galice.root",
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
    AliTPCseedGeant():TObject(),
      fXg(0.),
      fYg(0.),
      fZg(0.),
      fPx(0.),
      fPy(0.),
      fPz(0.),
      fAlpha(0.),
      fLabel(0),
      fSector(0){}

    AliTPCseedGeant(Double_t x=0.,Double_t y=0.,Double_t z=0.,
		    Double_t px=0.,Double_t py=0.,Double_t pz=0.,
		    Int_t lab=0);
    Int_t    GetLabel() const { return fLabel; }
    Double_t GetAlpha() const { return fAlpha; }      
    Double_t GetXG() const { return fXg; }
    Double_t GetYG() const { return fYg; }
    Double_t GetZG() const { return fZg; }
    Double_t GetXL() const 
      { return fXg*TMath::Cos(fAlpha)+fYg*TMath::Sin(fAlpha); }
    Double_t GetYL() const 
      { return -fXg*TMath::Sin(fAlpha)+fYg*TMath::Cos(fAlpha); }
    Double_t GetZL() const { return fZg; }
    Double_t GetPx() const { return fPx; } 
    Double_t GetPy() const { return fPy; } 
    Double_t GetPz() const { return fPz; } 
    Double_t GetPt() const { return TMath::Sqrt(fPx*fPx+fPy*fPy); }
    Double_t GetEta() const 
      { return -TMath::Log(TMath::Tan(0.25*TMath::Pi()-0.5*TMath::ATan(fPz/GetPt()))); }
    void     SetLabel(Int_t lab) { fLabel=lab; return; }
    Bool_t   InTPCAcceptance() const {
      if(TMath::Abs(GetZL()+(244.-GetXL())*fPz/GetPt())>252.) return kFALSE;
      return kTRUE;
    }

  private:
    Double_t fXg;     // global x of seed 
    Double_t fYg;     // global y of seed
    Double_t fZg;     // global z of seed
    Double_t fPx;     // global px
    Double_t fPy;     // global py
    Double_t fPz;     // global pz
    Double_t fAlpha;  // alpha angle
    Int_t    fLabel;  // track label
    Int_t    fSector; // TPC sector
  };
  //******* end of internal class ****************
  
 private:
  AliTPCtrackerParam & operator=(const AliTPCtrackerParam & );

  TString fEvFolderName;//! name of data folder

  Double_t        fBz;          // value of the z component of L3 field (Tesla)
  Int_t           fColl;        // collision code (0: PbPb6000; 1: pp)
  Bool_t          fSelAndSmear; // if kFALSE returns GEANT tracks 
                                // at TPC first hit 
  TString         fDBfileName;  // DataBase file name
  
  AliTPCtrack     fTrack;    // current track

  TTree          *fCovTree;  // tree with regularized cov matrices 
                             // for the current track
  TTree           *fCovTreePi[30];
  TTree           *fCovTreeKa[30];
  TTree           *fCovTreePr[30];
  TTree           *fCovTreeEl[30];
  TTree           *fCovTreeMu[30];

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
  void     CookdEdx(Double_t pt,Double_t eta);
  void     CookTrack(Double_t pt,Double_t eta);
  TMatrixD GetSmearingMatrix(Double_t *cc, Double_t pt,Double_t eta) const;
  void     InitializeKineGrid(Option_t *which);    
  void     MakeSeedsFromHits(AliTPC *tpc,TTree *th,TObjArray &seedArray) const;
  void     MakeSeedsFromRefs(TTree *ttr,
			     TObjArray &seedArray) const;
  Int_t    ReadAllData(const Char_t *inName);
  Int_t    ReadDBgrid(const Char_t *inName);
  Int_t    ReaddEdx(const Char_t *inName,Int_t pdg);
  Int_t    ReadEffs(const Char_t *inName);
  Int_t    ReadPulls(const Char_t *inName);
  Int_t    ReadRegParams(const Char_t *inName,Int_t pdg);
  Int_t    ReadCovTrees(const Char_t* inName); 
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













