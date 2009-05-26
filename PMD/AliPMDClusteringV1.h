#ifndef ALIPMDCLUSTERINGV1_H
#define ALIPMDCLUSTERINGV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDClusteringV1.h, Version 00        //
//                                                     //
//  Date   : September 26 2002                         //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//
// -- Author   : S.C. Phatak
// -- Modified : B.K. Nandi, Ajay Dash
//               T. Nayak, N. Sharma
//
#include "Rtypes.h"
#include "AliPMDClustering.h"

class TNtuple;
class TObjArray;
class AliPMDcluster;
class AliPMDcludata;

class AliPMDClusteringV1: public AliPMDClustering
{
 public:
  AliPMDClusteringV1();
  AliPMDClusteringV1(const AliPMDClusteringV1 &pmdclv1);
  AliPMDClusteringV1 &operator=(const AliPMDClusteringV1 &pmdclv1);
  virtual ~AliPMDClusteringV1();

  void     DoClust(Int_t idet, Int_t ismn, Int_t celltrack[][96],
		   Int_t cellpid[][96], Double_t celladc[][96],
		   TObjArray *pmdisocell, TObjArray *pmdcont);
  Int_t    CrClust(Double_t ave, Double_t cutoff, Int_t nmx1,
		   Int_t iord1[], Double_t edepcell[]);
  void     RefClust(Int_t incr, Double_t edepcell[]);
  Double_t Distance(Double_t x1, Double_t y1,
		    Double_t x2, Double_t y2);
  void     CalculateIsoCell(Int_t idet, Int_t ism,
			    Double_t celladc[][96], TObjArray *pmdisocell);
  void     SetEdepCut(Float_t decut);
  
 protected:
  
  TObjArray *fPMDclucont;    // carry cluster informations
  
  static const Double_t fgkSqroot3by2;  // fgkSqroot3by2 = sqrt(3.)/2.
  
  enum {
    kNMX    = 11424,     // no. of cells in a module
    kNDIMX  = 119,       // max no. of cells along x direction
    kNDIMY  = 96         // max no. of cells along axis at 60 deg with x axis
  };

  //Variables for association
  Int_t    fInfocl[2][kNDIMX][kNDIMY]; // cellwise information on the 
                                       // cluster to which the cell
  Int_t    fInfcl[3][kNMX];            // cluster information [0][i]
                                       // -- cluster number
  Double_t fCoord[2][kNDIMX][kNDIMY];

  Float_t  fCutoff; // Energy(ADC) cutoff per cell before clustering

  ClassDef(AliPMDClusteringV1,6) // Does clustering for PMD
};
#endif
