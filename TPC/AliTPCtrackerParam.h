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
#include <TMatrixD.h>
#include "alles.h"
#include "AliGausCorr.h"
#include "AliMagF.h"
#include "AliTPCkineGrid.h"
#include "AliTPCtrack.h"

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
  // For details:                                                           
  // http://www.pd.infn.it/alipd/talks/soft/adIII02/TPCtrackingParam.htm    
  //                                                                        
  // Test macro is: AliBarrelRec_TPCparam.C                                    
  //                                                                        
  //  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it     
  //                                                                        
  /////////////////////////////////////////////////////////////////////////
 public:
  AliTPCtrackerParam(const Int_t coll=0,const Double_t Bz=0.4);
  virtual ~AliTPCtrackerParam();

  // this function performs the parameterized tracking
  Int_t BuildTPCtracks(const TFile *inp, TFile *out,Int_t n=1);

  // these functions are used to create a DB of cov. matrices,
  // including regularization, efficiencies and dE/dx
  void  AllGeantTracks() { fSelAndSmear=kFALSE; return; }
  void  AnalyzedEdx(const Char_t *outName,Int_t pdg);
  void  AnalyzePulls(const Char_t *outName);
  void  CompareTPCtracks(const Char_t *galiceName="galice.root",
			 const Char_t *trkGeaName="AliTPCtracksGeant.root",
			 const Char_t *trkKalName="AliTPCtracksSorted.root",
			 const Char_t *covmatName="CovMatrix.root",
			 const Char_t *tpceffName="TPCeff.dat") const;
  void  DrawEffs(const Char_t *inName,Int_t pdg=211);
  void  DrawPulls(const Char_t *inName,Int_t pdg=211,Int_t par=0);
  void  MakeDataBase();
  void  MergeEvents(Int_t evFirst=1,Int_t evLast=1);
  void  RegularizeCovMatrix(const Char_t *outName,Int_t pdg);

  class AliTPCtrackParam : public AliTPCtrack {
  public:
    AliTPCtrackParam():AliTPCtrack(){}
    AliTPCtrackParam(const AliTPCtrack &t):AliTPCtrack(t){}

    void AssignMass(Double_t mass) {SetMass(mass); return;}
    
  private:
  
  };

 private:
  Double_t        fBz;          // value of the z component of L3 field (Tesla)
  Int_t           fColl;        // collision code (0: PbPb6000)
  Bool_t          fSelAndSmear; // if kFALSE returns GEANT tracks 
                                // at TPC first hit 
  TString         fDBfileName;  // DataBase file name
  
  AliTPCtrack     fTrack;    // current track

  TTree          *fCovTree;  // tree with regularized cov matrices 
                             // for the current track

  AliTPCkineGrid *fDBgrid;   // grid for the cov matrix look-up table  
  AliTPCkineGrid  fDBgridPi; //               "                  for pions  
  AliTPCkineGrid  fDBgridKa; //               "                  for kaons
  AliTPCkineGrid  fDBgridEl; //               "                  for electrons

  AliTPCkineGrid *fEff;   // TPC efficiencies for the current track
  AliTPCkineGrid  fEffPi; //           "        pions 
  AliTPCkineGrid  fEffKa; //           "        kaons 
  AliTPCkineGrid  fEffPr; //           "        protons 
  AliTPCkineGrid  fEffEl; //           "        electrons 
  AliTPCkineGrid  fEffMu; //           "        muons 

  AliTPCkineGrid *fPulls;      // pulls for the current track
  AliTPCkineGrid  fPullsPi[5]; //        "      pions
  AliTPCkineGrid  fPullsKa[5]; //        "      muons
  AliTPCkineGrid  fPullsEl[5]; //        "      electrons

  TMatrixD       *fRegPar;     // regularization parameters for the curr. track
  TMatrixD        fRegParPi;   //                  "        for pions         
  TMatrixD        fRegParKa;   //                  "        for kaons
  TMatrixD        fRegParEl;   //                  "        for electrons

  AliTPCkineGrid *fdEdxMean;   // dEdx mean for the current track
  AliTPCkineGrid  fdEdxMeanPi; //                "     pions
  AliTPCkineGrid  fdEdxMeanKa; //                "     kaons
  AliTPCkineGrid  fdEdxMeanPr; //                "     protons    
  AliTPCkineGrid  fdEdxMeanEl; //                "     electrons

  AliTPCkineGrid *fdEdxRMS;    // dEdx RMS for the current track
  AliTPCkineGrid  fdEdxRMSPi;  //                "     pions
  AliTPCkineGrid  fdEdxRMSKa;  //                "     kaons
  AliTPCkineGrid  fdEdxRMSPr;  //                "     protons    
  AliTPCkineGrid  fdEdxRMSEl;  //                "     electrons


  void     BuildTrack(Double_t alpha,Double_t x,Double_t y,Double_t z,
		      Double_t px,Double_t py,Double_t pz,Double_t pt,
		      Int_t ch);
  void     CookdEdx(Double_t pt,Double_t eta);
  void     CookTrack(Double_t pt,Double_t eta);
  Int_t    GetBin(Double_t pt,Double_t eta) const;
  TMatrixD GetSmearingMatrix(Double_t *cc, Double_t pt,Double_t eta) const;
  void     InitializeKineGrid(Option_t *which,Option_t *how);    
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













