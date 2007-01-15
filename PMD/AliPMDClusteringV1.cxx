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
//  Source File : PMDClusteringV1.cxx, Version 00      //
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
#include <TNtuple.h>
#include <TObjArray.h>
#include <stdio.h>

#include "AliPMDcluster.h"
#include "AliPMDClustering.h"
#include "AliPMDClusteringV1.h"
#include "AliLog.h"

ClassImp(AliPMDClusteringV1)

const Double_t AliPMDClusteringV1::fgkSqroot3by2=0.8660254;  // sqrt(3.)/2.

AliPMDClusteringV1::AliPMDClusteringV1():
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
AliPMDClusteringV1::~AliPMDClusteringV1()
{

}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::DoClust(Int_t idet, Int_t ismn, Double_t celladc[48][96], TObjArray *pmdcont)
{
  // main function to call other necessary functions to do clustering
  //
  AliPMDcluster *pmdcl = 0;
  /*
    int id and jd defined to read the input data.
    It is assumed that for data we have 0 <= id <= 48
    and 0 <= jd <=96
  */

  int i, i1, i2, j, nmx1, incr, id, jd;
  Int_t   celldataX[15], celldataY[15];
  Float_t clusdata[6];

  Double_t  cutoff, ave;

  const float ktwobysqrt3 = 1.1547; // 2./sqrt(3.)

  // ndimXr and ndimYr are different because of different module size

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
	  fCellTrNo[i][j] = -1;
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
	      fCellTrNo[i][j] = jd*10000+id; /* for association */
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      fEdepCell[i][j] = celladc[id][jd];
	      fCellTrNo[i][j] = id*10000+jd; /* for association */
	    }

	}
    }
  Order(); // order the data
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
      Float_t cluRAD   = (Float_t) fClusters[4][i1];
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

      clusdata[1]     = cluY0;
      clusdata[2]     = cluADC;
      clusdata[3]     = cluCELLS;
      clusdata[4]     = cluRAD;
      clusdata[5]     = 0.;

      //
      // Cells associated with a cluster
      //
      for (Int_t ihit = 0; ihit < 15; ihit++)
	{

	  if (ismn < 12)
	    {
	      celldataX[ihit] = fClTr[ihit][i1]%10000;
	      celldataY[ihit] = fClTr[ihit][i1]/10000;
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      celldataX[ihit] = fClTr[ihit][i1]/10000;
	      celldataY[ihit] = fClTr[ihit][i1]%10000;
	    }
	}
      //printf("%d %f %f\n",idet,cluXC,cluYC );
      pmdcl = new AliPMDcluster(idet, ismn, clusdata, celldataX, celldataY);
      pmdcont->Add(pmdcl);
    }
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::Order()
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
  
  TMath::Sort(kNMX,dd,iord1); //PH Using much better algorithm...
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
int AliPMDClusteringV1::CrClust(double ave, double cutoff, int nmx1)
{
  // Does crude clustering
  // Finds out only the big patch by just searching the
  // connected cells
  //
  int i,j,k,id1,id2,icl, numcell, clust[2][5000];
  int jd1,jd2, icell, cellcount;
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
      for(i=1; i<5000; i++)clust[0][i]=0;
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
	if(clust[0][i] != 0){
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
  //}

  return cellcount;
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::RefClust(int incr)
{
  // Does the refining of clusters
  // Takes the big patch and does gaussian fitting and
  // finds out the more refined clusters
  //
  int i, j, k, i1, i2, id, icl, ncl[4500], iord[4500], itest;
  int ihld;
  int ig, nsupcl, lev1[20], lev2[20];
  double x[4500], y[4500], z[4500], x1, y1, z1, x2, y2, z2, dist;
  double xc[4500], yc[4500], zc[4500], cells[4500], sum, rc[4500], rr;
  
  
  //asso
  Int_t t[4500],cellCount[4500];
  for(i=0; i<4500; i++)
    {
      t[i]=-1;
      cellCount[i]=0;
    }
  
  
  // fClno counts the final clusters
  // nsupcl =  # of superclusters; ncl[i]= # of cells in supercluster i
  // x, y and z store (x,y) coordinates of and energy deposited in a cell
  // xc, yc store (x,y) coordinates of the cluster center
  // zc stores the energy deposited in a cluster
  // rc is cluster radius
  // finally the cluster information is put in 2-dimensional array clusters
  //ofstream ofl1("checking.5",ios::app);
  fClno  = -1;
  nsupcl = -1;
  for(i=0; i<4500; i++){ncl[i]=-1;}
  for(i=0; i<= incr; i++){
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
  for(i=0; i<=nsupcl; i++) {
    if(ncl[i] == 0){
      id=id+1;
      icl=icl+1;
      // one  cell super-clusters --> single cluster
      // cluster center at the centyer of the cell
      // cluster radius = half cell dimension
      if (fClno >= 5000) {
	AliWarning("RefClust: Too many clusters! more than 5000");
	return;
      }
      fClno = fClno + 1;
      i1 = fInfcl[1][id];
      i2 = fInfcl[2][id];
      fClusters[0][fClno] = fCoord[0][i1][i2];
      fClusters[1][fClno] = fCoord[1][i1][i2];
      fClusters[2][fClno] = fEdepCell[i1][i2];
      fClusters[3][fClno] = 1.;
      fClusters[4][fClno] = 0.5;

      //association

      fClTr[0][fClno]=fCellTrNo[i1][i2];
      for(Int_t icltr=1;icltr<14;icltr++)
	{
	  fClTr[icltr][fClno]=-1;
	}
      
      //ofl1 << icl << " " << fCoord[0][i1][i2] << " " << fCoord[1][i1][i2] <<
      //" " << fEdepCell[i1][i2] << " " << fClusters[3][fClno] <<endl;
      
    }else if(ncl[i] == 1){
      // two cell super-cluster --> single cluster
      // cluster center is at ener. dep.-weighted mean of two cells
      // cluster radius == half cell dimension
      id   = id + 1;
      icl  = icl+1;
      if (fClno >= 5000) {
	AliWarning("RefClust: Too many clusters! more than 5000");
	return;
      }
      fClno = fClno+1;
      i1   = fInfcl[1][id];
      i2   = fInfcl[2][id];
      x1   = fCoord[0][i1][i2];
      y1   = fCoord[1][i1][i2];
      z1   = fEdepCell[i1][i2];

      //asso
      fClTr[0][fClno]=fCellTrNo[i1][i2];
      //

      id   = id+1;
      i1   = fInfcl[1][id];
      i2   = fInfcl[2][id];
      x2   = fCoord[0][i1][i2];
      y2   = fCoord[1][i1][i2];
      z2   = fEdepCell[i1][i2];

      //asso

      fClTr[1][fClno]=fCellTrNo[i1][i2];
      for(Int_t icltr=2;icltr<14;icltr++)
	{
	  fClTr[icltr][fClno] = -1;
	}
      //

      fClusters[0][fClno] = (x1*z1+x2*z2)/(z1+z2);
      fClusters[1][fClno] = (y1*z1+y2*z2)/(z1+z2);
      fClusters[2][fClno] = z1+z2;
      fClusters[3][fClno] = 2.;
      fClusters[4][fClno] = 0.5;


      //ofl1 << icl << " " << fClusters[0][fClno] << " " << fClusters[1][fClno]
      //  << " " << fClusters[2][fClno] << " " <<fClusters[3][fClno] <<endl;
    }
    else{
      
      //asso
      for(Int_t icg=0;icg<4500;icg++)
	{
	  cellCount[icg]=0;
	}
      //

      id      = id + 1;
      iord[0] = 0;
      // super-cluster of more than two cells - broken up into smaller
      // clusters gaussian centers computed. (peaks separated by > 1 cell)
      // Begin from cell having largest energy deposited This is first
      // cluster center
      i1      = fInfcl[1][id];
      i2      = fInfcl[2][id];
      x[0]    = fCoord[0][i1][i2];
      y[0]    = fCoord[1][i1][i2];
      z[0]    = fEdepCell[i1][i2];
      
      //asso
      t[0]=fCellTrNo[i1][i2];
      //

      iord[0] = 0;
      for(j=1;j<=ncl[i];j++){

	id      = id + 1;
	i1      = fInfcl[1][id];
	i2      = fInfcl[2][id];
	iord[j] = j;
	x[j]    = fCoord[0][i1][i2];
	y[j]    = fCoord[1][i1][i2];
	z[j]    = fEdepCell[i1][i2];

	//asso
	t[j]=fCellTrNo[i1][i2];
	//


      }
      // arranging cells within supercluster in decreasing order
      for(j=1;j<=ncl[i];j++){
	itest=0;
	ihld=iord[j];
	for(i1=0;i1<j;i1++){
	  if(itest == 0 && z[iord[i1]] < z[ihld]){
	    itest=1;
	    for(i2=j-1;i2>=i1;i2--){
	      iord[i2+1]=iord[i2];
	    }
	    iord[i1]=ihld;
	  }
	}
      }

      // compute the number of Gaussians and their centers ( first
      // guess )
      // centers must be separated by cells having smaller ener. dep.
      // neighbouring centers should be either strong or well-separated
      ig=0;
      xc[ig]=x[iord[0]];
      yc[ig]=y[iord[0]];
      zc[ig]=z[iord[0]];
      for(j=1;j<=ncl[i];j++){
	itest=-1;
	x1=x[iord[j]];
	y1=y[iord[j]];
	for(k=0;k<=ig;k++){
	  x2=xc[k]; y2=yc[k];
	  rr=Distance(x1,y1,x2,y2);
	  if( rr >= 1.1 && rr < 1.8 && z[iord[j]] > zc[k]/4.)
	    itest=itest+1;
	  if( rr >= 1.8 && rr < 2.1 && z[iord[j]] > zc[k]/10.)
	    itest=itest+1;
	  if( rr >= 2.1)itest=itest+1;
	}
	if(itest == ig){
	  ig=ig+1;
	  xc[ig]=x1;
	  yc[ig]=y1;
	  zc[ig]=z[iord[j]];
	}
      }
      // for(j=0; j<=ig; j++){
      //ofl1 << icl+j+1 << " " << xc[j] << " " <<yc[j] <<" "<<zc[j]<<endl;
      //}
      // GaussFit to adjust cluster parameters to minimize
      GaussFit(ncl[i], ig, x[0], y[0] ,z[0], xc[0], yc[0], zc[0], rc[0]);
      icl=icl+ig+1;
      // compute the number of cells belonging to each cluster.
      // cell is shared between several clusters ( if they are equidistant
      // from it ) in the ratio of cluster energy deposition
      for(j=0; j<=ig; j++){
	cells[j]=0.;
      }
      if(ig > 0){
	for(j=0; j<=ncl[i]; j++){
	  lev1[j]=0;
	  lev2[j]=0;
	  for(k=0; k<=ig; k++){
	    dist=Distance(x[j], y[j], xc[k], yc[k]);
	    if(dist < sqrt(3.) ){

	      //asso
	      fClTr[cellCount[k]][fClno+k+1]=t[j];
	      cellCount[k]++;
	      //

	      lev1[0]++;
	      i1=lev1[0];
	      lev1[i1]=k;
	    }else{
	      if(dist < 2.1){
		lev2[0]++;
		i1=lev2[0];
		lev2[i1]=k;
	      }
	    }
	  }
	  if(lev1[0] != 0){
	    if(lev1[0] == 1){cells[lev1[1]]=cells[lev1[1]]+1.;}
	    else{
	      sum=0.;
	      for(k=1; k<=lev1[0]; k++){
		sum=sum+zc[lev1[k]];
	      }
	      for(k=1; k<=lev1[0]; k++){
		cells[lev1[k]]=cells[lev1[k]]+zc[lev1[k]]/sum;
	      }
	    }
	  }else{
	    if(lev2[0] == 0){cells[lev2[1]]=cells[lev2[1]]+1.;}
	    else{
	      sum=0.;
	      for(k=1; k<=lev2[0]; k++){
		sum=sum+zc[lev2[k]];
	      }
	      for(k=1; k<=lev2[0]; k++){
		cells[lev2[k]]=cells[lev2[k]]+zc[lev2[k]]/sum;
	      }
	    }
	  }
	}
      }

      // zero rest of the cell array
      //asso
      for(k=0; k<=ig; k++)
	{
	  for(Int_t icltr=cellCount[k];icltr<14;icltr++)
	    {
	      fClTr[icltr][fClno]=-1;
	    }
	}
      //



      for(j=0; j<=ig; j++){
	if (fClno >= 5000) {
	  AliWarning("RefClust: Too many clusters! more than 5000");
	  return;
	}
	fClno               = fClno + 1;
	fClusters[0][fClno] = xc[j];
	fClusters[1][fClno] = yc[j];
	fClusters[2][fClno] = zc[j];
	fClusters[4][fClno] = rc[j];
	if(ig == 0){
	  fClusters[3][fClno] = ncl[i];
	}else{
	  fClusters[3][fClno] = cells[j];
	}
      }
    }
  }
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::GaussFit(Int_t ncell, Int_t nclust, Double_t &x, Double_t &y ,Double_t &z, Double_t &xc, Double_t &yc, Double_t &zc, Double_t &rc)
{
  // Does gaussian fitting
  //
  int i, j, i1, i2, novar, idd, jj;
  double xx[4500], yy[4500], zz[4500], xxc[4500], yyc[4500];
  double a[4500], b[4500], c[4500], d[4500], ha[4500], hb[4500];
  double hc[4500], hd[4500], zzc[4500], rrc[4500];
  int neib[4500][50];
  double sum, dx, dy, str, str1, aint, sum1, rr, dum;
  double x1, x2, y1, y2;
  str   = 0.;
  str1  = 0.;
  rr    = 0.3;
  novar = 0;
  j = 0;  // Just put not to see the compiler warning, BKN

  for(i=0; i<=ncell; i++)
    {
      xx[i] = *(&x+i);
      yy[i] = *(&y+i);
      zz[i] = *(&z+i);
      str   = str + zz[i];
    }
  for(i=0; i<=nclust; i++)
    {
      xxc[i] = *(&xc+i);
      yyc[i] = *(&yc+i);
      zzc[i] = *(&zc+i);
      str1   = str1 + zzc[i];
      rrc[i] = 0.5;
    }
  for(i=0; i<=nclust; i++)
    {
      zzc[i] = str/str1*zzc[i];
      ha[i]  = xxc[i];
      hb[i]  = yyc[i];
      hc[i]  = zzc[i];
      hd[i]  = rrc[i];
      x1     = xxc[i];
      y1     = yyc[i];
    }
  for(i=0; i<=ncell; i++){
    idd=0;
    x1=xx[i];
    y1=yy[i];
    for(j=0; j<=nclust; j++){
      x2=xxc[j];
      y2=yyc[j];
      if(Distance(x1,y1,x2,y2) <= 3.){ idd=idd+1; neib[i][idd]=j; }
    }
    neib[i][0]=idd;
  }
  sum=0.;
  for(i1=0; i1<=ncell; i1++){
    aint=0.;
    idd=neib[i1][0];
    for(i2=1; i2<=idd; i2++){
      jj=neib[i1][i2];
      dx=xx[i1]-xxc[jj];
      dy=yy[i1]-yyc[jj];
      dum=rrc[j]*rrc[jj]+rr*rr;
      aint=aint+exp(-(dx*dx+dy*dy)/dum)*zzc[idd]*rr*rr/dum;
    }
    sum=sum+(aint-zz[i1])*(aint-zz[i1])/str;
  }
//   jmax=nclust*1000;
//   if(nclust > 20)jmax=20000;
//   for(j=0; j<jmax; j++){
    str1=0.;
    for(i=0; i<=nclust; i++){
      a[i]=xxc[i]+0.6*(Ranmar()-0.5);
      b[i]=yyc[i]+0.6*(Ranmar()-0.5);
      c[i]=zzc[i]*(1.+(Ranmar()-0.5)*0.2);
      str1=str1+zzc[i];
      d[i]=rrc[i]*(1.+(Ranmar()-0.5)*0.1);
      if(d[i] < 0.25)d[i]=0.25;
    }
    for(i=0; i<=nclust; i++){ c[i]=c[i]*str/str1; }
    sum1=0.;
    for(i1=0; i1<=ncell; i1++){
      aint=0.;
      idd=neib[i1][0];
      for(i2=1; i2<=idd; i2++){
	jj=neib[i1][i2];
	dx=xx[i1]-a[jj];
	dy=yy[i1]-b[jj];
	dum=d[jj]*d[jj]+rr*rr;
	aint=aint+exp(-(dx*dx+dy*dy)/dum)*c[i2]*rr*rr/dum;
      }
      sum1=sum1+(aint-zz[i1])*(aint-zz[i1])/str;
    }

    if(sum1 < sum){
      for(i2=0; i2<=nclust; i2++){
	xxc[i2]=a[i2];
	yyc[i2]=b[i2];
	zzc[i2]=c[i2];
	rrc[i2]=d[i2];
	sum=sum1;
      }
    }
//   }
  for(j=0; j<=nclust; j++){
    *(&xc+j)=xxc[j];
    *(&yc+j)=yyc[j];
    *(&zc+j)=zzc[j];
    *(&rc+j)=rrc[j];
  }
}
// ------------------------------------------------------------------------ //
double AliPMDClusteringV1::Distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
// ------------------------------------------------------------------------ //
double AliPMDClusteringV1::Ranmar() const
{
  //  Universal random number generator proposed by Marsaglia and Zaman
  //  in report FSU-SCRI-87-50

  //  clock_t start;
  int ii, jj;
  static int i=96, j=32, itest=0, i1, i2, i3, i4, i5;
  static double u[97], c, cd, cm, s, t;
  static double uni;
  int count1,count2,idum;
  /*    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
  if (itest == 0) {
    //*******************************************************
    // following three lines if the seed to be provided by computer
    // start = time(NULL);
    // ii=start;
    // jj=start;
    //*******************************************************
    //following two lines for fixed seed ( during testing only. Else
    //use preceeing three lines
    ii=8263;
    jj=5726;
    if(ii > 31328 ) ii = ii - ( ii / 31328 ) * 31328;
    if(jj > 30081 ) jj = jj - ( jj / 30081 ) * 30081;
    itest=itest+1;
    if((( ii > 0 ) &&  ( ii <= 31328 )) && (( jj > 0 ) &&
					    ( jj <= 30081 ))){
      i1=ii/177+2; i2=ii-(i1-2)*177+2; i3=jj/169+1; i4=jj-(i3-1)*169;
      i4 = jj - (i3-1)*169;
      count1=0;
      while ( count1 < 97 ){
	s=0.;
	t=0.5;
	count2=0;
	while( count2 < 24 ){
	  idum=i1*i2/179;
	  idum=( i1*i2 - (i1*i2/179)*179 ) * i3;
	  i5=idum-(idum/179)*179;
	  i1=i2; i2=i3; i3=i5; idum=53*i4+1; i4=idum-(idum/169)*169;
	  if( i4*i5-((i4*i5)/64)*64 >= 32 ) s=s+t;
	  t=0.5*t;
	  count2=count2+1;
	}
	u[count1] = s;
	count1 = count1 +1;
      }
      c = 362436./16777216.;  cd = 7654321./16777216.;
      cm = 16777213./16777216.;
    }
    else{
      AliWarning("Wrong initialization");
    }
  }
  else{
    uni = u[i] - u[j];
    if( uni < 0.) uni = uni + 1;
    u[i] = uni;
    i = i -1;
    if( i < 0 ) i = 96;
    j = j - 1;
    if ( j < 0 ) j = 96;
    c = c - cd;
    if( c < 0. ) c = c+cm;
    uni = uni-c ;
    if( uni < 0. )uni = uni+1.;
  }
  return uni;
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::SetEdepCut(Float_t decut)
{
  fCutoff = decut;
}
// ------------------------------------------------------------------------ //
