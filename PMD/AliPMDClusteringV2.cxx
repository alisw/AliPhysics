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
   in TObjarray  and cluster information is in TObjArray.
   integer clno gives total number of clusters in the  supermodule.
   fClusters is the  global ( public ) variables.
   Others are local ( private ) to the code.
   At the moment, the data is read for whole detector ( all supermodules
   and pmd as well as cpv. This will have to be modify later )
   LAST UPDATE  :  October 23, 2002
-----------------------------------------------------------------------*/

#include <Riostream.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TArrayI.h>

#include "AliPMDcludata.h"
#include "AliPMDcluster.h"
#include "AliPMDisocell.h"
#include "AliPMDClustering.h"
#include "AliPMDClusteringV2.h"
#include "AliLog.h"

ClassImp(AliPMDClusteringV2)

const Double_t AliPMDClusteringV2::fgkSqroot3by2=0.8660254;  // sqrt(3.)/2.

AliPMDClusteringV2::AliPMDClusteringV2():
  fPMDclucont(new TObjArray()),
  fCutoff(0.0)
{
  for(int i = 0; i < kNDIMX; i++)
    {
      for(int j = 0; j < kNDIMY; j++)
	{
	  fCoord[0][i][j] = i+j/2.;
	  fCoord[1][i][j] = fgkSqroot3by2*j;
	}
    }
}
// ------------------------------------------------------------------------ //


AliPMDClusteringV2::AliPMDClusteringV2(const AliPMDClusteringV2& pmdclv2):
  AliPMDClustering(pmdclv2),
  fPMDclucont(0),
  fCutoff(0)
{
  // copy constructor
  AliError("Copy constructor not allowed ");
  
}
// ------------------------------------------------------------------------ //
AliPMDClusteringV2 &AliPMDClusteringV2::operator=(const AliPMDClusteringV2& /*pmdclv2*/)
{
  // copy constructor
  AliError("Assignment operator not allowed ");
  return *this;
}
// ------------------------------------------------------------------------ //
AliPMDClusteringV2::~AliPMDClusteringV2()
{
  delete fPMDclucont;
}
// ------------------------------------------------------------------------ //

void AliPMDClusteringV2::DoClust(Int_t idet, Int_t ismn, 
				 Int_t celltrack[48][96],
				 Int_t cellpid[48][96],
				 Double_t celladc[48][96],
				 TObjArray *pmdisocell, TObjArray *pmdcont)
{
  // main function to call other necessary functions to do clustering
  //
  AliPMDcluster *pmdcl = 0;

  const Float_t ktwobysqrt3 = 1.1547; // 2./sqrt(3.)
  const Int_t   kNmaxCell   = 19;     // # of cells surrounding a cluster center
  Int_t    i, j, nmx1, incr, id, jd;
  Int_t    ndimXr = 0;
  Int_t    ndimYr = 0;
  Int_t    celldataX[kNmaxCell], celldataY[kNmaxCell];
  Int_t    celldataTr[kNmaxCell], celldataPid[kNmaxCell];
  Float_t  celldataAdc[kNmaxCell];
  Float_t  clusdata[6];  
  Double_t cutoff, ave;
  Double_t edepcell[kNMX];


  // call the isolated cell search method

  CalculateIsoCell(idet, ismn, celladc, pmdisocell);



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
  
  for (i =0; i < kNMX; i++)
    {
     edepcell[i] = 0.;
    }
    
  for (id = 0; id < ndimXr; id++)
    {
      for (jd = 0; jd < ndimYr; jd++)
	{
	  j = jd;
	  i = id + (ndimYr/2-1) - (jd/2);
	  Int_t ij = i + j*kNDIMX;
	  if (ismn < 12)
	    {
	      edepcell[ij]    = celladc[jd][id];
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	     edepcell[ij]    = celladc[id][jd];
	    }

	}
    }

  Int_t iord1[kNMX];
  TMath::Sort((Int_t)kNMX,edepcell,iord1);// order the data
  cutoff = fCutoff; // cutoff used to discard cells having ener. dep.
  ave  = 0.;
  nmx1 = -1;

  for(i = 0;i < kNMX; i++)
    {
      if(edepcell[i] > 0.) 
	{
	  ave += edepcell[i];
	}
      if(edepcell[i] > cutoff )
	{
	  nmx1++;
	}
    }
  
  AliDebug(1,Form("Number of cells having energy >= %f are %d",cutoff,nmx1));
  
  if (nmx1 == 0) 
    {
      nmx1 = 1;
    }
  ave = ave/nmx1;
  
  AliDebug(1,Form("Number of cells in a SuperM = %d and Average = %f",
		  kNMX,ave));
  
  incr = CrClust(ave, cutoff, nmx1,iord1, edepcell);
  RefClust(incr,edepcell );
  
  Int_t nentries1 = fPMDclucont->GetEntries();
  AliDebug(1,Form("Detector Plane = %d  Serial Module No = %d Number of clusters = %d",idet, ismn, nentries1));
  AliDebug(1,Form("Total number of clusters/module = %d",nentries1));
  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
    {
      AliPMDcludata *pmdcludata = 
	(AliPMDcludata*)fPMDclucont->UncheckedAt(ient1);
      Float_t cluXC    = pmdcludata->GetClusX();
      Float_t cluYC    = pmdcludata->GetClusY();
      Float_t cluADC   = pmdcludata->GetClusADC();
      Float_t cluCELLS = pmdcludata->GetClusCells();
      Float_t cluSIGX  = pmdcludata->GetClusSigmaX();
      Float_t cluSIGY  = pmdcludata->GetClusSigmaY();
      
      Float_t cluY0    = ktwobysqrt3*cluYC;
      Float_t cluX0    = cluXC - cluY0/2.;
      
      // 
      // Cluster X centroid is back transformed
      //
      if (ismn < 12)
	{
	  clusdata[0] = cluX0 - (24-1) + cluY0/2.;
	}
      else if (ismn  >= 12 && ismn <= 23)
	{
	  clusdata[0] = cluX0 - (48-1) + cluY0/2.;
	}	  

      clusdata[1]     = cluY0;
      clusdata[2]     = cluADC;
      clusdata[3]     = cluCELLS;
      clusdata[4]     = cluSIGX;
      clusdata[5]     = cluSIGY;
      //
      // Cells associated with a cluster
      //
      for (Int_t ihit = 0; ihit < kNmaxCell; ihit++)
	{
	  Int_t dummyXY = pmdcludata->GetCellXY(ihit);
	 
	  Int_t celldumY   = dummyXY%10000;
	  Int_t celldumX   = dummyXY/10000;
          Float_t cellY    = (Float_t) celldumY/10;
	  Float_t cellX    = (Float_t) celldumX/10;

	  // 
	  // Cell X centroid is back transformed
	  //
	  if (ismn < 12)
	    {
	      celldataX[ihit] = (Int_t) ((cellX - (24-1) + cellY/2.) + 0.5);
	    }
	  else if (ismn  >= 12 && ismn <= 23)
	    {
	      celldataX[ihit] = (Int_t) ((cellX - (48-1) + cellY/2.) + 0.5 );
	    }	  
	  celldataY[ihit]   = (Int_t) (cellY + 0.5);

	  Int_t irow = celldataX[ihit];
	  Int_t icol = celldataY[ihit];

	  if ((irow >= 0 && irow < 48) && (icol >= 0 && icol < 96))
	    {
	      celldataTr[ihit]  = celltrack[irow][icol];
	      celldataPid[ihit] = cellpid[irow][icol];
	      celldataAdc[ihit] = (Float_t) celladc[irow][icol];
	    }
	  else
	    {
	      celldataTr[ihit]  = -1;
	      celldataPid[ihit] = -1;
	      celldataAdc[ihit] = -1;
	    }

	}

      pmdcl = new AliPMDcluster(idet, ismn, clusdata, celldataX, celldataY,
				celldataTr, celldataPid, celldataAdc);
      pmdcont->Add(pmdcl);
    }
  fPMDclucont->Delete();
}
// ------------------------------------------------------------------------ //
Int_t AliPMDClusteringV2::CrClust(Double_t ave, Double_t cutoff, Int_t nmx1,
				  Int_t iord1[], Double_t edepcell[])
{
  // Does crude clustering
  // Finds out only the big patch by just searching the
  // connected cells
  //

  Int_t i,j,k,id1,id2,icl, numcell;
  Int_t jd1,jd2, icell, cellcount;
  Int_t clust[2][5000];
  static Int_t neibx[6] = {1,0,-1,-1,0,1}, neiby[6] = {0,1,1,0,-1,-1};

  // neibx and neiby define ( incremental ) (i,j) for the neighbours of a
  // cell. There are six neighbours.
  // cellcount --- total number of cells having nonzero ener dep
  // numcell --- number of cells in a given supercluster
  
  AliDebug(1,Form("kNMX = %d nmx1 = %d kNDIMX = %d kNDIMY = %d ave = %f cutoff = %f",kNMX,nmx1,kNDIMX,kNDIMY,ave,cutoff));
  
  for (j=0; j < kNDIMX; j++)
    {
      for(k=0; k < kNDIMY; k++)
	{
	  fInfocl[0][j][k] = 0;
	  fInfocl[1][j][k] = 0;
	}
    }
 
  for(i=0; i < kNMX; i++)
    {
      fInfcl[0][i] = -1;
      
      j  = iord1[i];
      id2 = j/kNDIMX;
      id1 = j-id2*kNDIMX;
      
      if(edepcell[j] <= cutoff)
	{
	  fInfocl[0][id1][id2] = -1;
	}
    }
  // ---------------------------------------------------------------
  // crude clustering begins. Start with cell having largest adc
  // count and loop over the cells in descending order of adc count
  // ---------------------------------------------------------------
  icl       = -1;
  cellcount = -1;
  for(icell=0; icell <= nmx1; icell++)
    {
      j  = iord1[icell];
      id2 = j/kNDIMX;
      id1 = j-id2*kNDIMX;
      if(fInfocl[0][id1][id2] == 0 )
	{
	  // ---------------------------------------------------------------
	  // icl -- cluster #, numcell -- # of cells in it, clust -- stores
	  // coordinates of the cells in a cluster, fInfocl[0][i1][i2] is 1 for
	  // primary and 2 for secondary cells,
	  // fInfocl[1][i1][i2] stores cluster #
	  // ---------------------------------------------------------------
	  icl++;
	  numcell = 0;
	  cellcount++;
	  fInfocl[0][id1][id2]  = 1;
	  fInfocl[1][id1][id2]  = icl;
	  fInfcl[0][cellcount]  = icl;
	  fInfcl[1][cellcount]  = id1;
	  fInfcl[2][cellcount]  = id2;

	  clust[0][numcell]     = id1;
	  clust[1][numcell]     = id2;
	  for(i = 1; i < 5000; i++)
	    {
	      clust[0][i] = -1;
	    }
	  // ---------------------------------------------------------------
	  // check for adc count in neib. cells. If ne 0 put it in this clust
	  // ---------------------------------------------------------------
	  for(i = 0; i < 6; i++)
	    {
	    jd1 = id1 + neibx[i];
	    jd2 = id2 + neiby[i];
	    if( (jd1 >= 0 && jd1 < kNDIMX) && (jd2 >= 0 && jd2 < kNDIMY) &&
		fInfocl[0][jd1][jd2] == 0)
	      {
		numcell++;
		fInfocl[0][jd1][jd2] = 2;
		fInfocl[1][jd1][jd2] = icl;
		clust[0][numcell]    = jd1;
		clust[1][numcell]    = jd2;
		cellcount++;
		fInfcl[0][cellcount] = icl;
		fInfcl[1][cellcount] = jd1;
		fInfcl[2][cellcount] = jd2;
	      }
	    }
	  // ---------------------------------------------------------------
	  // check adc count for neighbour's neighbours recursively and
	  // if nonzero, add these to the cluster.
	  // ---------------------------------------------------------------
	  for(i = 1;i < 5000; i++)
	    { 
	      if(clust[0][i] != -1)
		{
		  id1 = clust[0][i];
		  id2 = clust[1][i];
		  for(j = 0; j < 6 ; j++)
		    {
		      jd1 = id1 + neibx[j];
		      jd2 = id2 + neiby[j];
		      if( (jd1 >= 0 && jd1 < kNDIMX) && 
			  (jd2 >= 0 && jd2 < kNDIMY) 
			  && fInfocl[0][jd1][jd2] == 0 )
			{
			  fInfocl[0][jd1][jd2] = 2;
			  fInfocl[1][jd1][jd2] = icl;
			  numcell++;
			  clust[0][numcell]    = jd1;
			  clust[1][numcell]    = jd2;
			  cellcount++;
			  fInfcl[0][cellcount] = icl;
			  fInfcl[1][cellcount] = jd1;
			  fInfcl[2][cellcount] = jd2;
			}
		    }
		}
	    }
	}
    }
  return cellcount;
}
// ------------------------------------------------------------------------ //
  void AliPMDClusteringV2::RefClust(Int_t incr, Double_t edepcell[])
{
  // Does the refining of clusters
  // Takes the big patch and does gaussian fitting and
  // finds out the more refined clusters

  const Float_t ktwobysqrt3 = 1.1547;
  const Int_t   kNmaxCell   = 19;

  AliPMDcludata *pmdcludata = 0;

  Int_t i12;
  Int_t    i, j, k, i1, i2, id, icl, itest, ihld;
  Int_t    ig, nsupcl, clno, clX,clY;
  Int_t    clxy[kNmaxCell];

  Float_t  clusdata[6];
  Double_t x1, y1, z1, x2, y2, z2, rr;

  Int_t kndim = incr + 1;

  TArrayI testncl;
  TArrayI testindex;

  Int_t    *ncl, *iord;

  Double_t *x, *y, *z, *xc, *yc, *zc, *cells, *rcl, *rcs;

  ncl   = new Int_t [kndim];
  iord  = new Int_t [kndim];
  x     = new Double_t [kndim];
  y     = new Double_t [kndim];
  z     = new Double_t [kndim];
  xc    = new Double_t [kndim];
  yc    = new Double_t [kndim];
  zc    = new Double_t [kndim];
  cells = new Double_t [kndim];
  rcl   = new Double_t [kndim];
  rcs   = new Double_t [kndim];
  
  for(Int_t kk = 0; kk < 15; kk++)
    {
      if( kk < 6 )clusdata[kk] = 0.;
    }
   
  // nsupcl =  # of superclusters; ncl[i]= # of cells in supercluster i
  // x, y and z store (x,y) coordinates of and energy deposited in a cell
  // xc, yc store (x,y) coordinates of the cluster center
  // zc stores the energy deposited in a cluster, rc is cluster radius

  clno   = -1;
  nsupcl = -1;

  for(i = 0; i < kndim; i++)
    {
      ncl[i] = -1;
    }
  for(i = 0; i <= incr; i++)
    {
      if(fInfcl[0][i] != nsupcl)
	{
	  nsupcl++;
	}
      if (nsupcl > 4500) 
	{
	  AliWarning("RefClust: Too many superclusters!");
	  nsupcl = 4500;
	  break;
	}
      ncl[nsupcl]++;
    }
  
  AliDebug(1,Form("Number of cells = %d Number of Superclusters = %d",
		  incr+1,nsupcl+1));
  
  id  = -1;
  icl = -1;
  for(i = 0; i <= nsupcl; i++)
    {
      if(ncl[i] == 0)
	{
	  id++;
	  icl++;
	  // one  cell super-clusters --> single cluster
	  // cluster center at the centyer of the cell
	  // cluster radius = half cell dimension
	  if (clno >= 5000) 
	    {
	      AliWarning("RefClust: Too many clusters! more than 5000");
	      return;
	    }
	  clno++;
	  i1          = fInfcl[1][id];
	  i2          = fInfcl[2][id];
	  i12         = i1 + i2*kNDIMX;
	  clusdata[0] = fCoord[0][i1][i2];
	  clusdata[1] = fCoord[1][i1][i2];
	  clusdata[2] = edepcell[i12];
	  clusdata[3] = 1.;
	  clusdata[4] = 0.0;
	  clusdata[5] = 0.0;
	  
	  //cell information
	  
	  clY = (Int_t)((ktwobysqrt3*fCoord[1][i1][i2])*10);
	  clX = (Int_t)((fCoord[0][i1][i2] - clY/20.)*10);
	  clxy[0] = clX*10000 + clY ;

	  for(Int_t icltr = 1; icltr < kNmaxCell; icltr++)
	    {
	      clxy[icltr] = -1;
	    }
	  pmdcludata  = new AliPMDcludata(clusdata,clxy);
	  fPMDclucont->Add(pmdcludata);
	  
	  
	}
      else if(ncl[i] == 1)
	{
	  // two cell super-cluster --> single cluster
	  // cluster center is at ener. dep.-weighted mean of two cells
	  // cluster radius == half cell dimension
	  id++;
	  icl++;
	  if (clno >= 5000) 
	    {
	      AliWarning("RefClust: Too many clusters! more than 5000");
	      return;
	    }
	  clno++;
	  i1   = fInfcl[1][id];
	  i2   = fInfcl[2][id];
	  i12  = i1 + i2*kNDIMX;
	  
	  x1   = fCoord[0][i1][i2];
	  y1   = fCoord[1][i1][i2];
	  z1   = edepcell[i12];
	  
	  id++;
	  i1   = fInfcl[1][id];
	  i2   = fInfcl[2][id];
	  i12  = i1 + i2*kNDIMX;
	  
	  x2   = fCoord[0][i1][i2];
	  y2   = fCoord[1][i1][i2];
	  z2   = edepcell[i12];
	  
	  clusdata[0] = (x1*z1+x2*z2)/(z1+z2);
	  clusdata[1] = (y1*z1+y2*z2)/(z1+z2);
	  clusdata[2] = z1+z2;
	  clusdata[3] = 2.;
	  clusdata[4] = (TMath::Sqrt(z1*z2))/(z1+z2);
	  clusdata[5] = 0.0;

          clY = (Int_t)((ktwobysqrt3*y1)*10);
	  clX = (Int_t)((x1 - clY/20.)*10);
	  clxy[0] = clX*10000 + clY ;

	  clY = (Int_t)((ktwobysqrt3*y2)*10);
	  clX = (Int_t)((x2 - clY/20.)*10);
	  clxy[1] = clX*10000 + clY ;

	  for(Int_t icltr = 2; icltr < kNmaxCell; icltr++)
	    {
	      clxy[icltr] = -1;
	    }
	  pmdcludata  = new AliPMDcludata(clusdata, clxy);
	  fPMDclucont->Add(pmdcludata);
	}
      else{
	id++;
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
	i12     = i1 + i2*kNDIMX;
	
	x[0]    = fCoord[0][i1][i2];
	y[0]    = fCoord[1][i1][i2];
	z[0]    = edepcell[i12];
	
	iord[0] = 0;
	for(j = 1; j <= ncl[i]; j++)
	  {
	    
	    id++;
	    i1      = fInfcl[1][id];
	    i2      = fInfcl[2][id];
	    i12     = i1 + i2*kNDIMX;
	    iord[j] = j;
	    x[j]    = fCoord[0][i1][i2];
	    y[j]    = fCoord[1][i1][i2];
	    z[j]    = edepcell[i12];
	  }
	
	// arranging cells within supercluster in decreasing order
	for(j = 1; j <= ncl[i];j++)
	  {
	    itest = 0;
	    ihld  = iord[j];
	    for(i1 = 0; i1 < j; i1++)
	      {
		if(itest == 0 && z[iord[i1]] < z[ihld])
		  {
		    itest = 1;
		    for(i2 = j-1;i2 >= i1;i2--)
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
	for(j = 1; j <= ncl[i]; j++)
	  {
	    itest = -1;
	    x1    = x[iord[j]];
	    y1    = y[iord[j]];
	    for(k = 0; k <= ig; k++)
	      {
		x2 = xc[k];
		y2 = yc[k];
		rr = Distance(x1,y1,x2,y2);
		//************************************************************
		// finetuning cluster splitting
		// the numbers zc/4 and zc/10 may need to be changed. 
		// Also one may need to add one more layer because our 
		// cells are smaller in absolute scale
		//************************************************************
		
		
		if( rr >= 1.1 && rr < 1.8 && z[iord[j]] > zc[k]/4.) itest++;
		if( rr >= 1.8 && rr < 2.1 && z[iord[j]] > zc[k]/10.) itest++;
		if( rr >= 2.1)itest++;
	      }
	    
	    if(itest == ig)
	      {
		ig++;
		xc[ig] = x1;
		yc[ig] = y1;
		zc[ig] = z[iord[j]];
	      }
	  }
	ClustDetails(ncl[i], ig, x, y ,z, xc, yc, zc, rcl, rcs, cells, 
		     testncl, testindex);
	
	Int_t pp = 0;
	for(j = 0; j <= ig; j++)
	  { 
	    clno++;
	    if (clno >= 5000)
	      {
		AliWarning("RefClust: Too many clusters! more than 5000");
		return;
	      }
	    clusdata[0] = xc[j];
	    clusdata[1] = yc[j];
	    clusdata[2] = zc[j];
	    clusdata[4] = rcl[j];
	    clusdata[5] = rcs[j];
	    if(ig == 0)
	      {
		clusdata[3] = ncl[i] + 1;
	      }
	    else
	      {
		clusdata[3] = cells[j];
	      }
	    // cell information
	    Int_t ncellcls =  testncl[j];
	    if( ncellcls < kNmaxCell )
	      {
		for(Int_t kk = 1; kk <= ncellcls; kk++)
		  {
		    Int_t ll =  testindex[pp];
                    clY = (Int_t)((ktwobysqrt3*y[ll])*10);
	            clX = (Int_t)((x[ll] - clY/20.)*10);
		    clxy[kk-1] = clX*10000 + clY ;

		    pp++;
		  }
		for(Int_t icltr = ncellcls ; icltr < kNmaxCell; icltr++)
		  {
		    clxy[icltr] = -1;
		  }
	      }
	    pmdcludata = new AliPMDcludata(clusdata, clxy);
	    fPMDclucont->Add(pmdcludata);
	  }
	testncl.Set(0);
	testindex.Set(0);
      }
    }
  delete [] ncl;
  delete [] iord;
  delete [] x;
  delete [] y;
  delete [] z;
  delete [] xc;
  delete [] yc;
  delete [] zc;
  delete [] cells;
  delete [] rcl;
  delete [] rcs;
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::ClustDetails(Int_t ncell, Int_t nclust, Double_t x[], 
				      Double_t y[], Double_t z[],Double_t xc[],
				      Double_t yc[], Double_t zc[],
				      Double_t rcl[], Double_t rcs[], 
				      Double_t cells[], TArrayI &testncl,
				      TArrayI &testindex)
{
  // function begins
  //

  Int_t kndim1 = ncell + 1;//ncell
  Int_t kndim2 = 20;
  Int_t kndim3 = nclust + 1;//nclust

  Int_t    i, j, k, i1, i2;
  Double_t x1, y1, x2, y2, rr, b, c, r1, r2;
  Double_t sumx, sumy, sumxy, sumxx, sum, sum1, sumyy;

  Double_t  *str, *str1, *xcl, *ycl, *cln; 
  Int_t    **cell;
  Int_t    ** cluster;
  Double_t **clustcell;
  str  = new Double_t [kndim3];
  str1 = new Double_t [kndim3];
  xcl  = new Double_t [kndim3];
  ycl  = new Double_t [kndim3];
  cln  = new Double_t [kndim3];

  clustcell = new Double_t *[kndim3];
  cell      = new Int_t    *[kndim3];
  cluster   = new Int_t    *[kndim1];
  for(i = 0; i < kndim1; i++)
    {
      cluster[i] = new Int_t [kndim2];
    }
  
  for(i = 0; i < kndim3; i++)
    {
      str[i]  = 0;
      str1[i] = 0;
      xcl[i]  = 0;
      ycl[i]  = 0;
      cln[i]  = 0;
      
      cell[i]    = new Int_t [kndim2];
      clustcell[i] = new Double_t [kndim1];	  
      for(j = 0; j < kndim1; j++)
	{
	  clustcell[i][j] = 0;
	}
      for(j = 0; j < kndim2; j++)
	{
	  cluster[i][j] = 0;
	  cell[i][j] = 0;
	}
    }
  
  if(nclust > 0)
    {
      // more than one cluster
      // checking cells shared between several  clusters.
      // First check if the cell is within
      // one cell unit ( nearest neighbour). Else, 
      // if it is within 1.74 cell units ( next nearest )
      // Else if it is upto 2 cell units etc.
      
      for (i = 0; i <= ncell; i++)
	{
	  x1            = x[i];
	  y1            = y[i];
	  cluster[i][0] = 0;

	  // distance <= 1 cell unit

	  for(j = 0; j <= nclust; j++)
	    {
	      x2 = xc[j];
	      y2 = yc[j];
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
		  x2 = xc[j];
		  y2 = yc[j];
		  rr = Distance(x1, y1, x2, y2);
		  if(rr <= TMath::Sqrt(3.))
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
		  x2 = xc[j];
		  y2 = yc[j];
		  rr = Distance(x1, y1, x2, y2);
		  if(rr <= 2.)
		    {
		      cluster[i][0]++;
		      i1             = cluster[i][0];
		      cluster[i][i1] = j;
		    }
		}
	    }
	  // one more
	  if(cluster[i][0] == 0)
	    {
	      for(j = 0; j <= nclust; j++)
		{
		  x2 = xc[j];
		  y2 = yc[j];
		  rr = Distance(x1, y1, x2, y2);
		  if(rr <= 2.7)
		    {
		      cluster[i][0]++;
		      i1             = cluster[i][0];
		      cluster[i][i1] = j;
		    }
		}
	    }
	}
      
      // computing cluster strength. Some cells are shared.
      for(i = 0; i <= ncell; i++)
	{
	  if(cluster[i][0] != 0)
	    {
	      i1 = cluster[i][0];
	      for(j = 1; j <= i1; j++)
		{
		  i2       = cluster[i][j];
		  str[i2] += z[i]/i1;
		}
	    }
	}
      
      for(k = 0; k < 5; k++)
	{
	  for(i = 0; i <= ncell; i++)
	    {
	      if(cluster[i][0] != 0)
		{
		  i1=cluster[i][0];
		  sum=0.;
		  for(j = 1; j <= i1; j++)
		    {
		      sum += str[cluster[i][j]];
		    }
		  
		  for(j = 1; j <= i1; j++)
		    {
		      i2               = cluster[i][j]; 
		      str1[i2]        +=  z[i]*str[i2]/sum;
		      clustcell[i2][i] = z[i]*str[i2]/sum;
		    }
		}
	    }
	  
	  
	  for(j = 0; j <= nclust; j++)
	    {
	      str[j]  = str1[j];
	      str1[j] = 0.;
	    }
	}
      
      for(i = 0; i <= nclust; i++)
	{
	  sumx = 0.;
	  sumy = 0.;
	  sum  = 0.;
	  sum1 = 0.;
	  for(j = 0; j <= ncell; j++)
	    {
	      if(clustcell[i][j] != 0)
		{
		  sumx  +=  clustcell[i][j]*x[j];
		  sumy  +=  clustcell[i][j]*y[j];
		  sum   +=  clustcell[i][j];
		  sum1  +=  clustcell[i][j]/z[j];
		}
	    }
	  //** xcl and ycl are cluster centroid positions ( center of gravity )
	  
	  xcl[i] = sumx/sum;
	  ycl[i] = sumy/sum;
	  cln[i] = sum1;
	  sumxx = 0.;
	  sumyy = 0.;
	  sumxy = 0.;
	  for(j = 0; j <= ncell; j++)
	    {
	      sumxx += clustcell[i][j]*(x[j]-xcl[i])*(x[j]-xcl[i])/sum;
	      sumyy += clustcell[i][j]*(y[j]-ycl[i])*(y[j]-ycl[i])/sum;
	      sumxy += clustcell[i][j]*(x[j]-xcl[i])*(y[j]-ycl[i])/sum;
	    }
	  b = sumxx+sumyy;
	  c = sumxx*sumyy-sumxy*sumxy;
	  // ******************r1 and r2 are major and minor axes ( r1 > r2 ). 
	  r1 = b/2.+TMath::Sqrt(b*b/4.-c);
	  r2 = b/2.-TMath::Sqrt(b*b/4.-c);
	  // final assignments to proper external variables
 	  xc[i]    = xcl[i];
	  yc[i]    = ycl[i];
	  zc[i]    = str[i];
	  cells[i] = cln[i];
	  rcl[i]   = r1;
	  rcs[i]   = r2;

	}
      
      //To get the cell position in a cluster
      
      for(Int_t ii=0; ii<= ncell; ii++)
	{
	  Int_t jj = cluster[ii][0]; 
	  for(Int_t kk=1; kk<= jj; kk++)
	    {
	      Int_t ll = cluster[ii][kk];
	      cell[ll][0]++;
	      cell[ll][cell[ll][0]] = ii;
	    }
	}
      
      testncl.Set(nclust+1);
      Int_t counter = 0;
      
      for(Int_t ii=0; ii <= nclust; ii++)
	{
	  testncl[ii] =  cell[ii][0];
	  counter += testncl[ii];
	}
      testindex.Set(counter);
      Int_t ll = 0;
      for(Int_t ii=0; ii<= nclust; ii++)
	{
	  for(Int_t jj = 1; jj<= testncl[ii]; jj++)
	    {
	      Int_t kk = cell[ii][jj];
	      testindex[ll] = kk;
	      ll++;
	    }
	}
      
    }
  else if(nclust == 0)
    {
      sumx = 0.;
      sumy = 0.;
      sum  = 0.;
      sum1 = 0.;
      i    = 0;
      for(j = 0; j <= ncell; j++)
	{
	  sumx += z[j]*x[j];
	  sumy += z[j]*y[j];
	  sum  += z[j];
	  sum1++;
	}
      xcl[i] = sumx/sum;
      ycl[i] = sumy/sum;
      cln[i] = sum1;
      sumxx  = 0.;
      sumyy  = 0.;
      sumxy  = 0.;
      for(j = 0; j <= ncell; j++)
	{
	  sumxx += clustcell[i][j]*(x[j]-xcl[i])*(x[j]-xcl[i])/sum;
	  sumyy += clustcell[i][j]*(y[j]-ycl[i])*(y[j]-ycl[i])/sum;
	  sumxy += clustcell[i][j]*(x[j]-xcl[i])*(y[j]-ycl[i])/sum;
	}
      b  = sumxx+sumyy;
      c  = sumxx*sumyy-sumxy*sumxy;
      r1 = b/2.+ TMath::Sqrt(b*b/4.-c);
      r2 = b/2.- TMath::Sqrt(b*b/4.-c);
      
      // To get the cell position in a cluster
      testncl.Set(nclust+1);
      testindex.Set(ncell+1);
      cell[0][0] = ncell + 1;
      testncl[0] = cell[0][0];
      Int_t ll   = 0;
      for(Int_t ii = 1; ii <= ncell; ii++)
	{
	  cell[0][ii]=ii;
	  Int_t kk = cell[0][ii];
	  testindex[ll] = kk;
	  ll++;
	}
      // final assignments
      xc[i]    = xcl[i];
      yc[i]    = ycl[i];
      zc[i]    = sum;
      cells[i] = cln[i];
      rcl[i]   = r1;
      rcs[i]   = r2;
    }
  for(i = 0; i < kndim3; i++)
    {
      delete [] clustcell[i];
      delete [] cell[i];
    }
  delete [] clustcell;
  delete [] cell;
  for(i = 0; i <kndim1 ; i++)
    {
      delete [] cluster[i];
    }
  delete [] cluster;
  delete [] str;
  delete [] str1;
  delete [] xcl;
  delete [] ycl;
  delete [] cln;
}

// ------------------------------------------------------------------------ //
Double_t AliPMDClusteringV2::Distance(Double_t x1, Double_t y1,
				      Double_t x2, Double_t y2)
{
  return TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::CalculateIsoCell(Int_t idet, Int_t ismn, Double_t celladc[][96], TObjArray *pmdisocell)
{
  // Does isolated cell search for offline calibration

  AliPMDisocell *isocell = 0;

  const Int_t kMaxRow = 48;
  const Int_t kMaxCol = 96;
  const Int_t kCellNeighbour = 6;

  Int_t id1, jd1;

  Int_t neibx[6] = {1,0,-1,-1,0,1};
  Int_t neiby[6] = {0,1,1,0,-1,-1};


  for(Int_t irow = 0; irow < kMaxRow; irow++)
    {
      for(Int_t icol = 0; icol < kMaxCol; icol++)
	{
	  if(celladc[irow][icol] > 0)
	    {
	      Int_t isocount = 0;
	      for(Int_t ii = 0; ii < kCellNeighbour; ii++)
		{
		  id1 = irow + neibx[ii];
		  jd1 = icol + neiby[ii];
		  Float_t adc = (Float_t) celladc[id1][jd1];
		  if(adc == 0.)
		    {
		      isocount++;
		      if(isocount == kCellNeighbour)
			{
			  Float_t cadc = (Float_t) celladc[irow][icol];

			  isocell = new AliPMDisocell(idet,ismn,irow,icol,cadc);
			  pmdisocell->Add(isocell);
			  
			}
		    }
		}  // neigh cell cond.
	    }
	}
    }


}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV2::SetEdepCut(Float_t decut)
{
  fCutoff = decut;
}
// ------------------------------------------------------------------------ //
