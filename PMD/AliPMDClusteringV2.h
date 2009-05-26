#ifndef ALIPMDCLUSTERINGV2_H
#define ALIPMDCLUSTERINGV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDClusteringV2.h,                   //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//
// Author      : S.C. Phatak
// Modified by : B.K. Nandi, Ajay Dash
//
#include "Rtypes.h"
#include "AliPMDClustering.h"

class TObjArray;
class TArrayI;
class AliPMDcluster;
class AliPMDcludata;
class AliPMDClusteringV2 : public AliPMDClustering
{
  
 public:
  AliPMDClusteringV2();
  AliPMDClusteringV2(const AliPMDClusteringV2 &pmdclv2);
  AliPMDClusteringV2 &operator=(const AliPMDClusteringV2 &pmdclv2);
  virtual ~AliPMDClusteringV2();
  
  void     DoClust(Int_t idet, Int_t ismn, Int_t celltrack[][96],
		   Int_t cellpid[][96], Double_t celladc[][96],
		   TObjArray *pmdisocell, TObjArray *pmdcont);
  Int_t    CrClust(Double_t ave, Double_t cutoff, Int_t nmx1,
		   Int_t iord1[], Double_t edepcell[]);
  void     RefClust(Int_t incr, Double_t edepcell[]);
	
  void     ClustDetails(Int_t ncell, Int_t nclust, Double_t x[],
			Double_t y[], Double_t z[], Double_t xc[],
			Double_t yc[], Double_t zc[],
			Double_t rcl[], Double_t rcs[], Double_t cells[],
			TArrayI &testncl, TArrayI &testindex);
  Double_t Distance(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
  void     SetEdepCut(Float_t decut);
  
 protected:
  
  TObjArray *fPMDclucont;
  
  static const Double_t fgkSqroot3by2;  // fgkSqroot3by2 = sqrt(3.)/2.
  enum {
    kNMX    = 11424, // no. of cells in a module
    kNDIMX  = 119,   // max no. of cells along x direction
    kNDIMY  = 96     // max no. of cells along axis at 60 deg with x axis
  };
  Int_t    fInfocl[2][kNDIMX][kNDIMY]; // cellwise information on the 
                                       // cluster to which the cell
  Int_t    fInfcl[3][kNMX];            // cluster information [0][i]
                                       // -- cluster number
  Double_t fCoord[2][kNDIMX][kNDIMY];

  Float_t fCutoff; // Energy(ADC) cutoff per cell before clustering
  
  ClassDef(AliPMDClusteringV2,4) // Does clustering for PMD
};
#endif
    
