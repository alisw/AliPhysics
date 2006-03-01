#ifndef ALIPMDCLUSTERINGV1_H
#define ALIPMDCLUSTERINGV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDClustering.h, Version 00          //
//                                                     //
//  Date   : September 26 2002                         //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//
/* --------------------------------------------------------------------
   Code developed by S. C. Phatak, Institute of Physics,
   Bhubaneswar 751 005 ( phatak@iopb.res.in ) Given the energy deposited
   ( or ADC value ) in each cell of supermodule ( pmd or cpv ), the code
   builds up superclusters and breaks them into clusters. The input is 
   in array d[ndimx][ndimy] and cluster information is in array
   clusters[5][5000]. integer clno gives total number of clusters in the
   supermodule.
   d, clno  and clusters are the only global ( public ) variables. Others
   are local ( private ) to the code.
   At the moment, the data is read for whole detector ( all supermodules
   and pmd as well as cpv. This will have to be modify later )
   LAST UPDATE  :  October 23, 2002
-----------------------------------------------------------------------*/
#include "Rtypes.h"
#include "AliPMDClustering.h"

class TNtuple;
class TObjArray;
class AliPMDcluster;

class AliPMDClusteringV1: public AliPMDClustering
{

 public:
  AliPMDClusteringV1();
  virtual ~AliPMDClusteringV1();

  void     DoClust(Int_t idet, Int_t ismn, Double_t celladc[][96],
		   TObjArray *pmdcont);
  void     Order();
  
  Int_t    CrClust(Double_t ave, Double_t cutoff, Int_t nmx1);
  void     RefClust(Int_t incr);
  void     GaussFit(Int_t ncell, Int_t nclust, Double_t &x,
		    Double_t &y, Double_t &z, Double_t &xc,
		    Double_t &yc, Double_t &zc, Double_t &rc);
  Double_t Distance(Double_t x1, Double_t y1,
		    Double_t x2, Double_t y2);
  Double_t Ranmar() const;
  void     SetEdepCut(Float_t decut);
  
 protected:

  static const Double_t fgkSqroot3by2;  // fgkSqroot3by2 = sqrt(3.)/2.
  /*enum {
    kNMX   = 4608,
    kNDIMX = 48,
    kNDIMY = 96
  };*/
  /*
    Proposed changes inNMX, kNDIMX and kNDIMY by S. C. Phatak to account
    for rectangular ( vs rhomboid ) geometry.
    To keep the clustering functional, we define a rhomboid which
    superscribes the rectangle. So we need to pad up dummy cells in x
    direction. The number of these cells is 96/2-1=47 in each row ( value
    of x ). For first two rows, all dummy cells are to the left. For
    every two rows add one cell to right and subtract one from left.
    So previous (i,j) values go over to ( i',j) i'=i+(96-j)/2-1
    Note we use C++ convention so i and j run from 0 to 47 or 95.
  */

  enum {
    kNMX    = 9120,
    kNDIMX  = 95,
    kNDIMY  = 96,
    kNDIMXr = 48,
    kNDIMYr = 96
  };
  /*
    kNMX   : # of cells in a supermodule
    kNDIMX : maximum number of cells along x direction (origin at one corner)
    kNDIMY : maximum number of cells along axis at 60 degrees with x axis
  */

  Double_t fEdepCell[kNDIMX][kNDIMY]; //energy(ADC) in each cell of the supermodule
  Double_t fClusters[5][5000]; // Cluster informations
  Int_t    fClno;   // number of clusters in a supermodule

  /*
    clusters[0][i] --- x position of the cluster center
    clusters[1][i] --- y position of the cluster center
    clusters[2][i] --- total energy in the cluster
    clusters[3][i] --- number of cells forming the cluster
                       ( possibly fractional )
    clusters[4][i] --- cluster radius
  */

  Int_t    fIord[2][kNMX]; // ordered list of i and j according to decreasing energy dep.
  Int_t    fInfocl[2][kNDIMX][kNDIMY]; // cellwise information on the cluster to which the cell
  Int_t    fInfcl[3][kNMX]; // cluster information [0][i] -- cluster number
  Double_t fCoord[2][kNDIMX][kNDIMY];

  /*
    fIord --- ordered list of i and j according to decreasing energy dep.
    fInfocl --- cellwise information on the cluster to which the cell
    belongs and whether it has largest energy dep. or not
    ( now redundant - probably )
    fInfcl ---  cluster information [0][i] -- cluster number
    [1][i] -- i of the cell
    [2][i] -- j of the cell
    coord --- x and y coordinates of center of each cell
  */

  Float_t fCutoff; // Energy(ADC) cutoff per cell before clustering

  ClassDef(AliPMDClusteringV1,1) // Does clustering for PMD
};
#endif
