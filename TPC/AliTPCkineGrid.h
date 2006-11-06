#ifndef ALITPCKINEGRID_H
#define ALITPCKINEGRID_H
/* Copyright(c) 2001-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-- Root headers ---
#include <TNamed.h>
#include <TMatrixD.h>
#include <TArrayD.h>
//-------------------

class AliTPCkineGrid : public TNamed {
  ////////////////////////////////////////////////////////////////////////
  // Class used by TPC tracking parameterization to handle the tracking
  // parameters (efficiencies, etc...) on a kinematic grid [pt,eta].
  // User has to provide the grid steps and the values of the parameter
  // in the points of the grid. The function GetValueAt(pt,eta) returns
  // the result of a linear interpolation at the point [pt,eta].
  //
  //  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it
  ////////////////////////////////////////////////////////////////////////
 public:
  
  AliTPCkineGrid();
  AliTPCkineGrid(Int_t npt,Int_t neta,Double_t* pt,Double_t* eta);
  AliTPCkineGrid(const AliTPCkineGrid& grid);
  AliTPCkineGrid &operator = (const AliTPCkineGrid & param);
  virtual ~AliTPCkineGrid();
  void     GetArrayEta(Double_t* eta) const;
  void     GetArrayPt(Double_t* pt) const;
  Int_t    GetBin(Double_t pt,Double_t eta) const;
  Int_t    GetBinsEta() const {return fNeta+1;}
  Int_t    GetBinsPt() const {return fNpt+1;}
  Double_t GetParam(Int_t i) const;
  Double_t GetParam(Int_t ipt,Int_t ieta) const {return (*fParams)(ipt,ieta);}
  Int_t    GetPointsEta() const {return fNeta;}
  Int_t    GetPointsPt() const {return fNpt;}
  Int_t    GetTotBins() const {return (fNeta+1)*(fNpt+1);}
  Int_t    GetTotPoints() const {return fNeta*fNpt;}
  Double_t GetValueAt(Double_t pt,Double_t eta) const;
  void     SetParam(Int_t i,Double_t par);
  void     SetParam(Int_t ipt,Int_t ieta,Double_t par) {
                                 (*fParams)(ipt,ieta)=par; return; }
  

 private:

  Int_t     fNpt;    // number of points in the grid:   pt
  Int_t     fNeta;   //               "                 eta
  TArrayD*  fPt;     // grid points in pt
  TArrayD*  fEta;    // grid points in eta
  TMatrixD* fParams; // matrix of parameters in the grid points  


  ClassDef(AliTPCkineGrid,1) // Parameters used by AliTPCtrackerParam 
};


#endif




