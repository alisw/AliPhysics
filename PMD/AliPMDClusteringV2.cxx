/***************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------//
//                                                     //
//  Source File : PMDClusteringV2.cxx                  //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//

/* --------------------------------------------------------------------
   Code developed by S. C. Phatak, Institute of Physics,
   Bhubaneswar 751 005 ( phatak@iopb.res.in ) Given the energy deposited
   ( or ADC value ) in each cell of supermodule ( pmd or cpv ), the code
   builds up superclusters and breaks them into clusters. The input is
   in array fEdepCell[kNDIMX][kNDIMY] and cluster information is in array
   fClusters[5][5000]. integer fClno gives total number of clusters in the
   supermodule.

   fEdepCell, fClno  and fClusters are the only global ( public ) variables.
   Others are local ( private ) to the code.
   At the moment, the data is read for whole detector ( all supermodules
   and pmd as well as cpv. This will have to be modify later )
   LAST UPDATE  :  October 23, 2002
-----------------------------------------------------------------------*/

#include <Riostream.h>
#include <TMath.h>
#include <TObjArray.h>
#include <stdio.h>

#include "AliPMDcluster.h"
#include "AliPMDClustering.h"
#include "AliPMDClusteringV2.h"
#include "AliLog.h"

ClassImp(AliPMDClusteringV2)

const Double_t AliPMDClusteringV2::fgkSqroot3by2=0.8660254;  // sqrt(3.)/2.

AliPMDClusteringV2::AliPMDClusteringV2():
  fClno(0),
  fCutoff(0.0)
{
  for(int i = 0; i < kNDIMX; i++)
    {
      for(int j = 0; j < kNDIMY; j++)
	{
	  fCoord[0][i][j] = i+j/2.;
	  fCoord[1][i][j] = fgkSqroot3by2*j;
	  fEdepCell[i][j] = 0;
	}
    }
}
// ------------------------------------------------------------------------ //
AliPMDClusteringV2::~AliPMDClusteringV2()
{

}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::DoClust(Int_t idet, Int_t ismn, Double_t celladc[48][96], TObjArray *pmdcont)
{
  // main function to call other necessary functions to do clustering
  //
  AliPMDcluster *pmdcl = 0;

  Int_t    i, i1, i2, j, nmx1, incr, id, jd;
  Int_t    celldataX[15], celldataY[15];
  Float_t  clusdata[6];
  Double_t cutoff, ave;

  const float ktwobysqrt3 = 1.1547; // 2./sqrt(3.)

  Int_t ndimXr =0;
  Int_t ndimYr =0;

  if (ismn < 12)
    {
      ndimXr = 96;
      ndimYr = 48;
    }
  else if (ismn >= 12 && ismn <= 23)
    {
      ndimXr = 48;
      ndimYr = 96;
    }

  for (Int_t i =0; i < kNDIMX; i++)
    {
      for (Int_t j =0; j < kNDIMY; j++)
	{
	  fEdepCell[i][j] = 0;
	}
    }


  for (id = 0; id < ndimXr; id++)
    {
      for (jd = 0; jd < ndimYr; jd++)
	{
	  j=jd;
	  i=id+(ndimYr/2-1)-(jd/2);
	  
	  if (ismn < 12)
	    {
	      fEdepCell[i][j] = celladc[jd][id];
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      fEdepCell[i][j] = celladc[id][jd];
	    }

	}
    }

  Order();          // order the data
  cutoff = fCutoff; // cutoff used to discard cells having ener. dep.
  ave=0.;
  nmx1=-1;

  for(j=0;j<kNMX; j++)
    {
      i1 = fIord[0][j];
      i2 = fIord[1][j];
      if (fEdepCell[i1][i2] > 0.) {ave = ave + fEdepCell[i1][i2];}
      if (fEdepCell[i1][i2] > cutoff ) nmx1 = nmx1 + 1;
    }
  // nmx1 --- number of cells having ener dep >= cutoff

  AliDebug(1,Form("Number of cells having energy >= %f are %d",cutoff,nmx1));

  if (nmx1 == 0) nmx1 = 1;
  ave=ave/nmx1;

  AliDebug(1,Form("Number of cells in a SuperM = %d and Average = %f",
		  kNMX,ave));
	   
  incr = CrClust(ave, cutoff, nmx1);
  RefClust(incr);

  AliDebug(1,Form("Detector Plane = %d  Serial Module No = %d Number of clusters = %d",idet, ismn, fClno));
  
  for(i1=0; i1<=fClno; i1++)
    {
      Float_t cluXC    = (Float_t) fClusters[0][i1];
      Float_t cluYC    = (Float_t) fClusters[1][i1];
      Float_t cluADC   = (Float_t) fClusters[2][i1];
      Float_t cluCELLS = (Float_t) fClusters[3][i1];
      Float_t sigmaX   = (Float_t) fClusters[4][i1];
      Float_t sigmaY   = (Float_t) fClusters[5][i1];
      Float_t cluY0    = ktwobysqrt3*cluYC;
      Float_t cluX0    = cluXC - cluY0/2.;
      // 
      // Cluster X centroid is back transformed
      //
      if (ismn < 12)
	{
	  clusdata[0] = cluX0 - (24-1) + cluY0/2.;
	}
      else if (ismn >= 12 && ismn <= 23)
	{
	  clusdata[0] = cluX0 - (48-1) + cluY0/2.;
	}	  

      clusdata[1]      = cluY0;
      clusdata[2]      = cluADC;
      clusdata[3]      = cluCELLS;
      clusdata[4]      = sigmaX;
      clusdata[5]      = sigmaY;

      //
      // Cells associated with a cluster
      //
      for (Int_t ihit = 0; ihit < 15; ihit++)
	{
	  celldataX[ihit] = 1;  // dummy nos. -- will be changed
	  celldataY[ihit] = 1;  // dummy nos. -- will be changed
	}

      pmdcl = new AliPMDcluster(idet, ismn, clusdata, celldataX, celldataY);
      pmdcont->Add(pmdcl);
    }
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::Order()
{
  // Sorting algorithm
  // sorts the ADC values from higher to lower
  //
  double dd[kNMX];
  // matrix fEdepCell converted into
  // one dimensional array dd. adum a place holder for double
  int i, j, i1, i2, iord1[kNMX];
  // information of
  // ordering is stored in iord1, original array not ordered
  //
  // define arrays dd and iord1
  for(i1=0; i1 < kNDIMX; i1++)
    {
      for(i2=0; i2 < kNDIMY; i2++)
	{
	  i        = i1 + i2*kNDIMX;
	  iord1[i] = i;
	  dd[i]    = fEdepCell[i1][i2];
	}
    }
  // sort and store sorting information in iord1

  TMath::Sort(kNMX,dd,iord1);

  // store the sorted information in fIord for later use
  for(i=0; i<kNMX; i++)
    {
      j  = iord1[i];
      i2 = j/kNDIMX;
      i1 = j-i2*kNDIMX;
      fIord[0][i]=i1;
      fIord[1][i]=i2;
    }
}
// ------------------------------------------------------------------------ //
Int_t AliPMDClusteringV2::CrClust(Double_t ave, Double_t cutoff, Int_t nmx1)
{
  // Does crude clustering
  // Finds out only the big patch by just searching the
  // connected cells
  //

  int i,j,k,id1,id2,icl, numcell;
  int jd1,jd2, icell, cellcount;
  int clust[2][5000];
  static int neibx[6]={1,0,-1,-1,0,1}, neiby[6]={0,1,1,0,-1,-1};

  // neibx and neiby define ( incremental ) (i,j) for the neighbours of a
  // cell. There are six neighbours.
  // cellcount --- total number of cells having nonzero ener dep
  // numcell --- number of cells in a given supercluster
  // ofstream ofl0("cells_loc",ios::out);
  // initialize fInfocl[2][kNDIMX][kNDIMY]

  AliDebug(1,Form("kNMX = %d nmx1 = %d kNDIMX = %d kNDIMY = %d ave = %f cutoff = %f",kNMX,nmx1,kNDIMX,kNDIMY,ave,cutoff));
  
  for (j=0; j < kNDIMX; j++){
    for(k=0; k < kNDIMY; k++){
      fInfocl[0][j][k] = 0;
      fInfocl[1][j][k] = 0;
    }
  }
  for(i=0; i < kNMX; i++){
    fInfcl[0][i] = -1;
    id1=fIord[0][i];
    id2=fIord[1][i];
    if(fEdepCell[id1][id2] <= cutoff){fInfocl[0][id1][id2]=-1;}
  }
  // ---------------------------------------------------------------
  // crude clustering begins. Start with cell having largest adc
  // count and loop over the cells in descending order of adc count
  // ---------------------------------------------------------------
  icl=-1;
  cellcount=-1;
  for(icell=0; icell <= nmx1; icell++){
    id1=fIord[0][icell];
    id2=fIord[1][icell];
    if(fInfocl[0][id1][id2] == 0 ){
      // ---------------------------------------------------------------
      // icl -- cluster #, numcell -- # of cells in it, clust -- stores
      // coordinates of the cells in a cluster, fInfocl[0][i1][i2] is 1 for
      // primary and 2 for secondary cells,
      // fInfocl[1][i1][i2] stores cluster #
      // ---------------------------------------------------------------
      icl=icl+1;
      numcell=0;
      cellcount = cellcount + 1;
      fInfocl[0][id1][id2]=1;
      fInfocl[1][id1][id2]=icl;
      fInfcl[0][cellcount]=icl;
      fInfcl[1][cellcount]=id1;
      fInfcl[2][cellcount]=id2;

      clust[0][numcell]=id1;
      clust[1][numcell]=id2;
      for(i=1; i<5000; i++)clust[0][i] = -1;
      // ---------------------------------------------------------------
      // check for adc count in neib. cells. If ne 0 put it in this clust
      // ---------------------------------------------------------------
      for(i=0; i<6; i++){
	jd1=id1+neibx[i];
	jd2=id2+neiby[i];
	if( (jd1 >= 0 && jd1 < kNDIMX) && (jd2 >= 0 && jd2 < kNDIMY) &&
	    fInfocl[0][jd1][jd2] == 0){
	  numcell=numcell+1;
	  fInfocl[0][jd1][jd2]=2;
	  fInfocl[1][jd1][jd2]=icl;
	  clust[0][numcell]=jd1;
	  clust[1][numcell]=jd2;
	  cellcount=cellcount+1;
	  fInfcl[0][cellcount]=icl;
	  fInfcl[1][cellcount]=jd1;
	  fInfcl[2][cellcount]=jd2;
	}
      }
      // ---------------------------------------------------------------
      // check adc count for neighbour's neighbours recursively and
      // if nonzero, add these to the cluster.
      // ---------------------------------------------------------------
      for(i=1;i < 5000;i++){
	if(clust[0][i] != -1){
	  id1=clust[0][i];
	  id2=clust[1][i];
	  for(j=0; j<6 ; j++){
	    jd1=id1+neibx[j];
	    jd2=id2+neiby[j];
	    if( (jd1 >= 0 && jd1 < kNDIMX) && (jd2 >= 0 && jd2 < kNDIMY) &&
		fInfocl[0][jd1][jd2] == 0 ){
	      fInfocl[0][jd1][jd2] = 2;
	      fInfocl[1][jd1][jd2] = icl;
	      numcell              = numcell + 1;
	      clust[0][numcell]    = jd1;
	      clust[1][numcell]    = jd2;
	      cellcount            = cellcount+1;
	      fInfcl[0][cellcount] = icl;
	      fInfcl[1][cellcount] = jd1;
	      fInfcl[2][cellcount] = jd2;
	    }
	  }
	}
      }
    }
  }
  //  for(icell=0; icell<=cellcount; icell++){
  //    ofl0 << fInfcl[0][icell] << " " << fInfcl[1][icell] << " " <<
  //      fInfcl[2][icell] << endl;
  //  }
  return cellcount;
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::RefClust(Int_t incr)
{
  // Does the refining of clusters
  // Takes the big patch and does gaussian fitting and
  // finds out the more refined clusters
  //

  const Int_t kndim = 4500;

  int i, j, k, i1, i2, id, icl, itest;
  int ihld;
  int ig, nsupcl;
  int ncl[kndim], iord[kndim];

  double x1, y1, z1, x2, y2, z2;
  double rr;

  double x[kndim], y[kndim], z[kndim];
  double xc[kndim], yc[kndim], zc[kndim], cells[kndim];
  double rcl[kndim], rcs[kndim];

  // fClno counts the final clusters
  // nsupcl =  # of superclusters; ncl[i]= # of cells in supercluster i
  // x, y and z store (x,y) coordinates of and energy deposited in a cell
  // xc, yc store (x,y) coordinates of the cluster center
  // zc stores the energy deposited in a cluster
  // rc is cluster radius
  // finally the cluster information is put in 2-dimensional array clusters
  // ofstream ofl1("checking.5",ios::app);

  fClno  = -1;
  nsupcl = -1;
  for(i=0; i<4500; i++){ncl[i]=-1;}
  for(i=0; i<incr; i++){
    if(fInfcl[0][i] != nsupcl){ nsupcl=nsupcl+1; }
    if (nsupcl > 4500) {
      AliWarning("RefClust: Too many superclusters!");
      nsupcl = 4500;
      break;
    }
    ncl[nsupcl]=ncl[nsupcl]+1;
  }

  AliDebug(1,Form("Number of cells = %d Number of Superclusters = %d",
		  incr+1,nsupcl+1));

  id=-1;
  icl=-1;
  for(i=0; i<nsupcl; i++){
    if(ncl[i] == 0){
      id++;
      icl++;
      // one  cell super-clusters --> single cluster
      // cluster center at the centyer of the cell
      // cluster radius = half cell dimension
      if (fClno >= 5000) {
	AliWarning("RefClust: Too many clusters! more than 5000");
	return;
      }
      fClno++;
      i1 = fInfcl[1][id];
      i2 = fInfcl[2][id];
      fClusters[0][fClno] = fCoord[0][i1][i2];
      fClusters[1][fClno] = fCoord[1][i1][i2];
      fClusters[2][fClno] = fEdepCell[i1][i2];
      fClusters[3][fClno] = 1.;
      fClusters[4][fClno] = 0.0;
      fClusters[5][fClno] = 0.0;
      //ofl1 << icl << " " << fCoord[0][i1][i2] << " " << fCoord[1][i1][i2] <<
      //" " << fEdepCell[i1][i2] << " " << fClusters[3][fClno] <<endl;
    }else if(ncl[i] == 1){
      // two cell super-cluster --> single cluster
      // cluster center is at ener. dep.-weighted mean of two cells
      // cluster radius == half cell dimension
      id++;
      icl++;
      if (fClno >= 5000) {
	AliWarning("RefClust: Too many clusters! more than 5000");
	return;
      }
      fClno++;
      i1   = fInfcl[1][id];
      i2   = fInfcl[2][id];
      x1   = fCoord[0][i1][i2];
      y1   = fCoord[1][i1][i2];
      z1   = fEdepCell[i1][i2];

      id++;
      i1   = fInfcl[1][id];
      i2   = fInfcl[2][id];
      x2   = fCoord[0][i1][i2];
      y2   = fCoord[1][i1][i2];
      z2   = fEdepCell[i1][i2];

      fClusters[0][fClno] = (x1*z1+x2*z2)/(z1+z2);
      fClusters[1][fClno] = (y1*z1+y2*z2)/(z1+z2);
      fClusters[2][fClno] = z1+z2;
      fClusters[3][fClno] = 2.;
      fClusters[4][fClno] = sqrt(z1*z2)/(z1+z2);
      fClusters[5][fClno] = 0;  // sigma large nonzero, sigma small zero

      //ofl1 << icl << " " << fClusters[0][fClno] << " " << fClusters[1][fClno]
      //   << " " << fClusters[2][fClno] << " " <<fClusters[3][fClno] <<endl;
    }
    else{
      id      = id + 1;
      iord[0] = 0;
      // super-cluster of more than two cells - broken up into smaller
      // clusters gaussian centers computed. (peaks separated by > 1 cell)
      // Begin from cell having largest energy deposited This is first
      // cluster center
      // *****************************************************************
      // NOTE --- POSSIBLE MODIFICATION: ONE MAY NOT BREAKING SUPERCLUSTERS
      // IF NO. OF CELLS IS NOT TOO LARGE ( SAY 5 OR 6 )
      // SINCE WE EXPECT THE SUPERCLUSTER 
      // TO BE A SINGLE CLUSTER
      //*******************************************************************

      i1      = fInfcl[1][id];
      i2      = fInfcl[2][id];
      x[0]    = fCoord[0][i1][i2];
      y[0]    = fCoord[1][i1][i2];
      z[0]    = fEdepCell[i1][i2];
      iord[0] = 0;
      for(j=1;j<=ncl[i];j++){

	id      = id + 1;
	i1      = fInfcl[1][id];
	i2      = fInfcl[2][id];
	iord[j] = j;
	x[j]    = fCoord[0][i1][i2];
	y[j]    = fCoord[1][i1][i2];
	z[j]    = fEdepCell[i1][i2];
      }
      // arranging cells within supercluster in decreasing order
      for(j=1;j<=ncl[i];j++)
	{
	  itest = 0;
	  ihld  = iord[j];
	  for(i1=0; i1<j; i1++)
	    {
	      if(itest == 0 && z[iord[i1]] < z[ihld])
		{
		  itest = 1;
		  for(i2=j-1;i2>=i1;i2--)
		    {
		      iord[i2+1] = iord[i2];
		    }
		  iord[i1] = ihld;
		}
	    }
	}

      // compute the number of clusters and their centers ( first
      // guess )
      // centers must be separated by cells having smaller ener. dep.
      // neighbouring centers should be either strong or well-separated
      ig     = 0;
      xc[ig] = x[iord[0]];
      yc[ig] = y[iord[0]];
      zc[ig] = z[iord[0]];
      for(j=1;j<=ncl[i];j++){
	itest = -1;
	x1    = x[iord[j]];
	y1    = y[iord[j]];
	for(k=0;k<=ig;k++){
	  x2 = xc[k];
	  y2 = yc[k];
	  rr = Distance(x1,y1,x2,y2);
	  //***************************************************************
	  // finetuning cluster splitting
	  // the numbers zc/4 and zc/10 may need to be changed. 
	  // Also one may need to add one more layer because our 
	  // cells are smaller in absolute scale
	  //****************************************************************


	  if( rr >= 1.1 && rr < 1.8 && z[iord[j]] > zc[k]/4.)
	    itest++;
	  if( rr >= 1.8 && rr < 2.1 && z[iord[j]] > zc[k]/10.)
	    itest++;
	  if( rr >= 2.1)itest++;
	}
	if(itest == ig){
	  ig++;
	  xc[ig] = x1;
	  yc[ig] = y1;
	  zc[ig] = z[iord[j]];
	}
      }

      ClustDetails(ncl[i], ig, x[0], y[0] ,z[0], xc[0], yc[0], zc[0],
		   rcl[0], rcs[0], cells[0]);

      icl = icl + ig + 1;

      for(j=0; j<=ig; j++)
	{
	  if (fClno >= 5000)
	    {
	      AliWarning("RefClust: Too many clusters! more than 5000");
	      return;
	    }
	  fClno++;
	  fClusters[0][fClno] = xc[j];
	  fClusters[1][fClno] = yc[j];
	  fClusters[2][fClno] = zc[j];
	  fClusters[4][fClno] = rcl[j];
	  fClusters[5][fClno] = rcs[j];
	  if(ig == 0)
	    {
	      fClusters[3][fClno] = ncl[i];
	    }
	  else
	    {
	      fClusters[3][fClno] = cells[j];
	    }
	}


    }
  }
}


// ------------------------------------------------------------------------ //

void AliPMDClusteringV2::ClustDetails(Int_t ncell, Int_t nclust,
				      Double_t &x, Double_t &y, Double_t &z,
				      Double_t &xc, Double_t &yc, Double_t &zc,
				      Double_t &rcl, Double_t &rcs,
				      Double_t &cells)
{
  // function begins
  //
  
  const Int_t kndim1 = 4500;
  const Int_t kndim2 = 10;
  const Int_t kndim3 = 100;

  int i, j, k, i1, i2;
  int cluster[kndim1][kndim2];
  
  double x1, y1, x2, y2, rr;
  double sumx, sumy, sumxy, sumxx;
  double sum, sum1, sumyy;
  double b, c, r1, r2;

  double xx[kndim1], yy[kndim1], zz[kndim1];
  double xxc[kndim1], yyc[kndim1];

  double str[kndim1];

  double str1[kndim1];
  double xcl[kndim1], ycl[kndim1], cln[kndim1];
  double clustcell[kndim1][kndim3];

  for(i=0; i<=nclust; i++){
   xxc[i]=*(&xc+i); 
   yyc[i]=*(&yc+i); 
   str[i]=0.; 
   str1[i]=0.;
  }
  for(i=0; i<=ncell; i++){
    xx[i]=*(&x+i); 
    yy[i]=*(&y+i); 
    zz[i]=*(&z+i);
  }
  // INITIALIZE 
  for(i=0; i<4500; i++){
    for(j=0; j<100; j++){
      clustcell[i][j]=0.;
    }
  }

  // INITIALIZE
  for(i=0;i<4500;i++){
    for(j=0;j<10;j++){
      cluster[i][j]=0;
    }
  }


  if(nclust > 0){
    // more than one cluster
    // checking cells shared between several  clusters.
    // First check if the cell is within
    // one cell unit ( nearest neighbour). Else, 
    // if it is within 1.74 cell units ( next nearest )
    // Else if it is upto 2 cell units etc.

    for (i=0; i<=ncell; i++){
      x1            = xx[i];
      y1            = yy[i];
      cluster[i][0] = 0;
      // distance <= 1 cell unit
      for(j=0; j<=nclust; j++)
	{
	  x2 = xxc[j];
	  y2 = yyc[j];
	  rr = Distance(x1, y1, x2, y2);
	  if(rr <= 1.)
	    {
	      cluster[i][0]++;
	      i1             = cluster[i][0];
	      cluster[i][i1] = j;
	    }
	}
      // next nearest neighbour
      if(cluster[i][0] == 0)
	{
	  for(j=0; j<=nclust; j++)
	    {
	      x2 = xxc[j];
	      y2 = yyc[j];
	      rr = Distance(x1, y1, x2, y2);
	      if(rr <= sqrt(3.))
		{
		  cluster[i][0]++;
		  i1             = cluster[i][0];
		  cluster[i][i1] = j;
		}
	    }
	}
      // next-to-next nearest neighbour
      if(cluster[i][0] == 0)
	{
	  for(j=0; j<=nclust; j++)
	    {
	      x2 = xxc[j];
	      y2 = yyc[j];
	      rr = Distance(x1, y1, x2, y2);
	      if(rr <= 2.)
		{
		  cluster[i][0]++;
		  i1 = cluster[i][0];
		  cluster[i][i1] = j;
		}
	    }
	}
      // one more
      if(cluster[i][0] == 0)
	{
	  for(j=0; j<=nclust; j++)
	    {
	      x2 = xxc[j];
	      y2 = yyc[j];
	      rr = Distance(x1, y1, x2, y2);
	      if(rr <= 2.7)
		{
		  cluster[i][0]++;
		  i1 = cluster[i][0];
		  cluster[i][i1] = j;
		}
	    }
	}
    }


    // computing cluster strength. Some cells are shared.
    for(i=0; i<=ncell; i++){
      if(cluster[i][0] != 0){
	i1 = cluster[i][0];
	for(j=1; j<=i1; j++){
	  i2      = cluster[i][j];
	  str[i2] = str[i2]+zz[i]/i1;
	}
      }
    }

    for(k=0; k<5; k++)
      {
	for(i=0; i<=ncell; i++)
	  {
	    if(cluster[i][0] != 0)
	      {
		i1=cluster[i][0];
		sum=0.;
		for(j=1; j<=i1; j++)
		  {
		    sum=sum+str[cluster[i][j]];
		  }

		for(j=1; j<=i1; j++)
		  {
		    i2 = cluster[i][j]; 
		    str1[i2]         = str1[i2] + zz[i]*str[i2]/sum;
		    clustcell[i2][i] = zz[i]*str[i2]/sum;
		  }
	      }
	  }


	  for(j=0; j<=nclust; j++)
	    {
	      str[j]=str1[j];
	      str1[j]=0.;
	    }
      }

    for(i=0; i<=nclust; i++){
      sumx = 0.;
      sumy = 0.;
      sum  = 0.;
      sum1 = 0.;
      for(j=0; j<=ncell; j++){
	if(clustcell[i][j] != 0){
	  sumx = sumx+clustcell[i][j]*xx[j];
	  sumy = sumy+clustcell[i][j]*yy[j];
	  sum  = sum+clustcell[i][j];
	  sum1 = sum1+clustcell[i][j]/zz[j];
	}
      }
      //***** xcl and ycl are cluster centroid positions ( center of gravity )

      xcl[i] = sumx/sum;
      ycl[i] = sumy/sum;
      cln[i] = sum1;
      sumxx = 0.;
      sumyy = 0.;
      sumxy = 0.;
      for(j=0; j<=ncell; j++){
	sumxx = sumxx+clustcell[i][j]*(xx[j]-xcl[i])*(xx[j]-xcl[i])/sum;
	sumyy = sumyy+clustcell[i][j]*(yy[j]-ycl[i])*(yy[j]-ycl[i])/sum;
	sumxy = sumxy+clustcell[i][j]*(xx[j]-xcl[i])*(yy[j]-ycl[i])/sum;
      }
      b = sumxx+sumyy;
      c = sumxx*sumyy-sumxy*sumxy;
      // ******************r1 and r2 are major and minor axes ( r1 > r2 ). 
      r1 = b/2.+sqrt(b*b/4.-c);
      r2 = b/2.-sqrt(b*b/4.-c);
      // final assignments to proper external variables
      *(&xc + i) = xcl[i];
      *(&yc + i) = ycl[i];
      *(&zc + i) = str[i];
      *(&cells + i) = cln[i];
      *(&rcl+i) = r1;
      *(&rcs+i) = r2;
    }
  }else{
    sumx = 0.;
    sumy = 0.;
    sum  = 0.;
    sum1 = 0.;
    i    = 0;
    for(j=0; j<=ncell; j++){
      sumx = sumx+zz[j]*xx[j];
      sumy = sumy+zz[j]*yy[j];
      sum  = sum+zz[j];
      sum1 = sum1+1.;
    }
    xcl[i] = sumx/sum;
    ycl[i] = sumy/sum;
    cln[i] = sum1;
    sumxx  = 0.;
    sumyy  = 0.;
    sumxy  = 0.;
    for(j=0; j<=ncell; j++){
      sumxx = sumxx+clustcell[i][j]*(xx[j]-xcl[i])*(xx[j]-xcl[i])/sum;
      sumyy = sumyy+clustcell[i][j]*(yy[j]-ycl[i])*(yy[j]-ycl[i])/sum;
      sumxy = sumxy+clustcell[i][j]*(xx[j]-xcl[i])*(yy[j]-ycl[i])/sum;
    }
    b  = sumxx+sumyy;
    c  = sumxx*sumyy-sumxy*sumxy;
    r1 = b/2.+sqrt(b*b/4.-c);
    r2 = b/2.-sqrt(b*b/4.-c);
    // final assignments
    *(&xc + i)    = xcl[i];
    *(&yc + i)    = ycl[i];
    *(&zc + i)    = str[i];
    *(&cells + i) = cln[i];
    *(&rcl+i)     = r1;
    *(&rcs+i)     = r2;
  }
}

// ------------------------------------------------------------------------ //
Double_t AliPMDClusteringV2::Distance(Double_t x1, Double_t y1,
				      Double_t x2, Double_t y2)
{
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::SetEdepCut(Float_t decut)
{
  fCutoff = decut;
}
// ------------------------------------------------------------------------ //
