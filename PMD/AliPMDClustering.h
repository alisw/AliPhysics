#ifndef PMDClustering_H
#define PMDClustering_H
//-----------------------------------------------------//
//                                                     //
//  Header File : PMDClustering.h, Version 00          //
//                                                     //
//  Date   : September 26 2002                         //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//
/* 
   --------------------------------------------------------------------
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
-----------------------------------------------------------------------
*/


#include <Riostream.h> // define cout stream
#include <stdlib.h>   // defines exit() functions
#include <time.h> // for time function
#include <math.h> // for mathematical functions
#include "Rtypes.h"

class TNtuple;
class TObjArray;
class AliPMDcluster;
class AliPMDClustering
{
  
 protected:

  static const double pi;
  static const double sqrth;  // sqrth = sqrt(3.)/2.
  enum {
    nmx=5184,
    ndimx=72,
    ndimy=72
  };

  /*
    nmx : # of cells in a supermodule
    ndimx : maximum number of cells along x direction (origin at one corner)
    ndimy : maximum number of cells along axis at 60 degrees with x axis
  */

  double d[ndimx][ndimy], clusters[5][5000]; 
  int clno;

  /*
    d ---- energy deposited ( or ADC ) in each cell of the supermodule
    clno --- number of clusters in a supermodule
    A cell is defined in terms of two integers (i,j) giving the its location
    clusters[0][i] --- x position of the cluster center
    clusters[1][i] --- y position of the cluster center
    clusters[2][i] --- total energy in the cluster
    clusters[3][i] --- number of cells forming the cluster 
                       ( possibly fractional )
    clusters[4][i] --- cluster radius
    One corner of the supermodule is chosen as the origin
  */


  int iord[2][nmx], infocl[2][ndimx][ndimy], infcl[3][nmx];
  double coord[2][ndimx][ndimy];

  /* 
     iord --- ordered list of i and j according to decreasing energy dep.
     infocl --- cellwise information on the cluster to which the cell
     belongs and whether it has largest energy dep. or not
     ( now redundant - probably )
     infcl ---  cluster information [0][i] -- cluster number
     [1][i] -- i of the cell
     [2][i] -- j of the cell
     coord --- x and y coordinates of center of each cell
  */

  Int_t fMessage;

 public:
  AliPMDClustering();
  virtual ~AliPMDClustering();
  
  void DoClust(int, int, double [][72], TObjArray *);
  int crclust(double, double, int, int);
  void refclust(int, int, int);
  double ranmar();
  void order(int);
  double Dist(double, double, double, double);
  void gaussfit(int, int, double &, double &, double &, double &, double &, 
		double &, double &);
  void SetMessage(Int_t);

  ClassDef(AliPMDClustering,1)
};
#endif
