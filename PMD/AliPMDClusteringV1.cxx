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
   in array edepcell[kNMX] and cluster information is in a
   TObjarray. Integer clno gives total number of clusters in the
   supermodule.

   fClusters is the only global ( public ) variables.
   Others are local ( private ) to the code.
   At the moment, the data is read for whole detector ( all supermodules
   and pmd as well as cpv. This will have to be modify later )
   LAST UPDATE  :  October 23, 2002
-----------------------------------------------------------------------*/

#include <Riostream.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include "TRandom.h"
#include <stdio.h>

#include "AliPMDcludata.h"
#include "AliPMDcluster.h"
#include "AliPMDisocell.h"
#include "AliPMDClustering.h"
#include "AliPMDClusteringV1.h"
#include "AliLog.h"


ClassImp(AliPMDClusteringV1)

const Double_t AliPMDClusteringV1::fgkSqroot3by2=0.8660254;  // sqrt(3.)/2.

AliPMDClusteringV1::AliPMDClusteringV1():
  fPMDclucont(new TObjArray()),
  fCutoff(0.0)
{
  for(Int_t i = 0; i < kNDIMX; i++)
    {
      for(Int_t j = 0; j < kNDIMY; j++)
	{
	  fCoord[0][i][j] = i+j/2.;
	  fCoord[1][i][j] = fgkSqroot3by2*j;
	}
    }
}
// ------------------------------------------------------------------------ //
AliPMDClusteringV1::AliPMDClusteringV1(const AliPMDClusteringV1& pmdclv1):
  AliPMDClustering(pmdclv1),
  fPMDclucont(0),
  fCutoff(0)
{
  // copy constructor
  AliError("Copy constructor not allowed ");
  
}
// ------------------------------------------------------------------------ //
AliPMDClusteringV1 &AliPMDClusteringV1::operator=(const AliPMDClusteringV1& /*pmdclv1*/)
{
  // copy constructor
  AliError("Assignment operator not allowed ");
  return *this;
}
// ------------------------------------------------------------------------ //
AliPMDClusteringV1::~AliPMDClusteringV1()
{
  delete fPMDclucont;
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::DoClust(Int_t idet, Int_t ismn,
				 Int_t celltrack[48][96],
				 Int_t cellpid[48][96],
				 Double_t celladc[48][96],
				 TObjArray *pmdisocell, TObjArray *pmdcont)
{
  // main function to call other necessary functions to do clustering
  //

  AliPMDcluster *pmdcl = 0;

  const float ktwobysqrt3 = 1.1547; // 2./sqrt(3.)
  const Int_t kNmaxCell   = 19;     // # of cells surrounding a cluster center

  Int_t    i,  j, nmx1, incr, id, jd;
  Int_t    celldataX[kNmaxCell], celldataY[kNmaxCell];
  Int_t    celldataTr[kNmaxCell], celldataPid[kNmaxCell];
  Float_t  celldataAdc[kNmaxCell];
  Float_t  clusdata[6];
  Double_t cutoff, ave;
  Double_t edepcell[kNMX];
  
  
  Double_t cellenergy[11424];
  

  // call the isolated cell search method

  FindIsoCell(idet, ismn, celladc, pmdisocell);

  // ndimXr and ndimYr are different because of different module size

  Int_t ndimXr = 0;
  Int_t ndimYr = 0;

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

  for (i =0; i < 11424; i++)
  {
      cellenergy[i] = 0.;
  }

  Int_t kk = 0;
  for (i = 0; i < kNDIMX; i++)
    {
      for (j = 0; j < kNDIMY; j++)
	{
	  edepcell[kk] = 0.;
	  kk++;
	}
    }

  for (id = 0; id < ndimXr; id++)
    {
      for (jd = 0; jd < ndimYr; jd++)
	{
	  j = jd;
	  i = id+(ndimYr/2-1)-(jd/2);

	  Int_t ij = i + j*kNDIMX;
	  
	  if (ismn < 12)
	    {
	      cellenergy[ij]    = celladc[jd][id];
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      cellenergy[ij]    = celladc[id][jd];
	    }
	}
    }
  
  for (i = 0; i < kNMX; i++)
    {
      edepcell[i] = cellenergy[i];
    }

  Int_t iord1[kNMX];
  TMath::Sort((Int_t)kNMX,edepcell,iord1);// order the data
  cutoff = fCutoff;                       // cutoff to discard cells
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

  if (nmx1 == 0) nmx1 = 1;
  ave = ave/nmx1;
  AliDebug(1,Form("Number of cells in a SuperM = %d and Average = %f",
		  kNMX,ave));
  
  incr = CrClust(ave, cutoff, nmx1,iord1, edepcell );
  RefClust(incr,edepcell);
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
      else if ( ismn >= 12 && ismn <= 23)
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
	  Int_t cellrow = pmdcludata->GetCellXY(ihit)/10000;
	  Int_t cellcol = pmdcludata->GetCellXY(ihit)%10000;

	  if (ismn < 12)
	    {
	      celldataX[ihit] = cellrow - (24-1) + int(cellcol/2.);
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      celldataX[ihit] = cellrow - (48-1) + int(cellcol/2.);
	    }
	  
	  celldataY[ihit] = cellcol;
	  
	  Int_t irow = celldataX[ihit];
	  Int_t icol = celldataY[ihit];

	  if (ismn < 12)
	    {
	      if ((irow >= 0 && irow < 96) && (icol >= 0 && icol < 48))
		{
		  celldataTr[ihit]  = celltrack[icol][irow];
		  celldataPid[ihit] = cellpid[icol][irow];
		  celldataAdc[ihit] = (Float_t) celladc[icol][irow];
		}
	      else
		{
		  celldataTr[ihit]  = -1;
		  celldataPid[ihit] = -1;
		  celldataAdc[ihit] = -1;
		}
	    }
	  else if (ismn >= 12 && ismn < 24)
	    {
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
	  
	}
      
      pmdcl = new AliPMDcluster(idet, ismn, clusdata, celldataX, celldataY,
				celldataTr, celldataPid, celldataAdc);
      pmdcont->Add(pmdcl);
    }
  
  fPMDclucont->Delete();
}
// ------------------------------------------------------------------------ //
Int_t AliPMDClusteringV1::CrClust(Double_t ave, Double_t cutoff, Int_t nmx1,
				  Int_t iord1[], Double_t edepcell[])
{
  // Does crude clustering 
  // Finds out only the big patch by just searching the
  // connected cells
  //
  const Int_t kndim = 4609;
  Int_t i,j,k,id1,id2,icl, numcell, clust[2][kndim];
  Int_t jd1,jd2, icell, cellcount;
  static Int_t neibx[6]={1,0,-1,-1,0,1}, neiby[6]={0,1,1,0,-1,-1};

  AliDebug(1,Form("kNMX = %d nmx1 = %d kNDIMX = %d kNDIMY = %d ave = %f cutoff = %f",kNMX,nmx1,kNDIMX,kNDIMY,ave,cutoff));
  
  for (j = 0; j < kNDIMX; j++)
    {
      for(k = 0; k < kNDIMY; k++)
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

  for(icell = 0; icell <= nmx1; icell++)
    {
      j  = iord1[icell];
      id2 = j/kNDIMX;
      id1 = j-id2*kNDIMX;

      if(fInfocl[0][id1][id2] == 0 )
	{
	  icl++;
	  numcell = 0;
	  cellcount++; 
	  fInfocl[0][id1][id2] = 1;
	  fInfocl[1][id1][id2] = icl;
	  fInfcl[0][cellcount] = icl;
	  fInfcl[1][cellcount] = id1;
	  fInfcl[2][cellcount] = id2;

	  clust[0][numcell] = id1;
	  clust[1][numcell] = id2;
	  
	  for(i = 1; i < kndim; i++)
	    {
	      clust[0][i]=0;
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
	  for(i = 1; i < kndim;i++)
	    {
	      if(clust[0][i] != 0)
		{
		  id1 = clust[0][i];
		  id2 = clust[1][i];
		  for(j = 0; j < 6 ; j++)
		    {
		      jd1 = id1 + neibx[j];
		      jd2 = id2 + neiby[j];
		      if( (jd1 >= 0 && jd1 < kNDIMX) && 
			  (jd2 >= 0 && jd2 < kNDIMY) &&
			  fInfocl[0][jd1][jd2] == 0 )
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
void AliPMDClusteringV1::RefClust(Int_t incr, Double_t edepcell[])
{
  // Does the refining of clusters
  // Takes the big patch and does gaussian fitting and
  // finds out the more refined clusters
  //
  
  AliPMDcludata *pmdcludata = 0;
  
  const Int_t kNmaxCell   = 19;    // # of cells surrounding a cluster center
  
  Int_t ndim = incr + 1;
  
  Int_t    *ncl  = 0x0;
  Int_t    *clxy = 0x0;  
  Int_t    i12, i22;
  Int_t    i, j, k, i1, i2, id, icl,  itest,ihld, ig, nsupcl,clno, t1, t2;
  Float_t  clusdata[6];
  Double_t x1, y1, z1, x2, y2, z2, rr;
  
  ncl   = new Int_t [ndim];
  clxy  = new Int_t [kNmaxCell];
  
  // Initialisation  
  for(i = 0; i<ndim; i++)
    {
      ncl[i]  = -1; 
      if (i < 6) clusdata[i] = 0.;
      if (i < kNmaxCell) clxy[i]    = 0;
    }

  // clno counts the final clusters
  // nsupcl =  # of superclusters; ncl[i]= # of cells in supercluster i
  // x, y and z store (x,y) coordinates of and energy deposited in a cell
  // xc, yc store (x,y) coordinates of the cluster center
  // zc stores the energy deposited in a cluster
  // rc is cluster radius
  
  clno  = -1;
  nsupcl = -1;

  for(i = 0; i <= incr; i++)
    {
      if(fInfcl[0][i] != nsupcl)
	{
	  nsupcl++;
	}
      if (nsupcl > ndim) 
	{
	  AliWarning("RefClust: Too many superclusters!");
	  nsupcl = ndim;
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
	  if (clno >= 4608) 
	    {
	      AliWarning("RefClust: Too many clusters! more than 4608");
	      return;
	    }
	  clno++;
	  i1 = fInfcl[1][id];
	  i2 = fInfcl[2][id];
	  
	  i12 = i1 + i2*kNDIMX;
	  
	  clusdata[0] = fCoord[0][i1][i2];
	  clusdata[1] = fCoord[1][i1][i2];
	  clusdata[2] = edepcell[i12];
	  clusdata[3] = 1.;
	  clusdata[4] = 0.5;
	  clusdata[5] = 0.0;

	  clxy[0] = i1*10000 + i2;
	  
	  for(Int_t icltr = 1; icltr < kNmaxCell; icltr++)
	    {
	      clxy[icltr] = -1;
	    }

	  pmdcludata  = new AliPMDcludata(clusdata,clxy);
	  fPMDclucont->Add(pmdcludata);
	}
      else if(ncl[i] == 1) 
	{
	  id++;
	  icl++;
	  if (clno >= 4608) 
	    {
	      AliWarning("RefClust: Too many clusters! more than 4608");
	      return;
	    }
	  clno++;
	  i1   = fInfcl[1][id];
	  i2   = fInfcl[2][id];
	  i12  = i1 + i2*kNDIMX;

	  x1   = fCoord[0][i1][i2];
	  y1   = fCoord[1][i1][i2];
	  z1   = edepcell[i12];

	  clxy[0] = i1*10000 + i2;
	  
	  id++;
	  i1   = fInfcl[1][id];
	  i2   = fInfcl[2][id];

	  i22  = i1 + i2*kNDIMX;
	  x2   = fCoord[0][i1][i2];
	  y2   = fCoord[1][i1][i2];
	  z2   = edepcell[i22];

	  clxy[1] = i1*10000 + i2;
	  

	  for(Int_t icltr = 2; icltr < kNmaxCell; icltr++)
	    {
	      clxy[icltr] = -1;
	    }
	  
	  clusdata[0] = (x1*z1+x2*z2)/(z1+z2);
	  clusdata[1] = (y1*z1+y2*z2)/(z1+z2);
	  clusdata[2] = z1+z2;
	  clusdata[3] = 2.;
	  clusdata[4] = 0.5;
	  clusdata[5] = 0.0;
	  pmdcludata  = new AliPMDcludata(clusdata,clxy);
	  fPMDclucont->Add(pmdcludata);
	}
      else
	{
	  
	  Int_t    *iord, *tc, *t;
	  Double_t *x, *y, *z, *xc, *yc, *zc;

	  iord = new Int_t [ncl[i]+1];
	  tc   = new Int_t [ncl[i]+1];
	  t    = new Int_t [ncl[i]+1];
	  
	  x    = new Double_t [ncl[i]+1];
	  y    = new Double_t [ncl[i]+1];
	  z    = new Double_t [ncl[i]+1];
	  xc   = new Double_t [ncl[i]+1];
	  yc   = new Double_t [ncl[i]+1];
	  zc   = new Double_t [ncl[i]+1];
	  
	  for( k = 0; k < ncl[i]+1; k++)
	    {
	      iord[k] = -1;
	      t[k]    = -1;
	      tc[k]   = -1;
	      x[k]    = -1;
	      y[k]    = -1;
	      z[k]    = -1;
	      xc[k]   = -1;
	      yc[k]   = -1;
	      zc[k]   = -1;
	    }
	  id++;
	  // super-cluster of more than two cells - broken up into smaller
	  // clusters gaussian centers computed. (peaks separated by > 1 cell)
	  // Begin from cell having largest energy deposited This is first
	  // cluster center
	  i1      = fInfcl[1][id];
	  i2      = fInfcl[2][id];
	  i12     = i1 + i2*kNDIMX;
	  
	  x[0]    = fCoord[0][i1][i2];
	  y[0]    = fCoord[1][i1][i2];
	  z[0]    = edepcell[i12];
	  t[0]    = i1*10000 + i2;
	  

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
	      t[j]    = i1*10000 + i2;

	    }
	  
	  // arranging cells within supercluster in decreasing order
	  
	  for(j = 1;j <= ncl[i]; j++)
	    {
	      itest = 0;
	      ihld  = iord[j];
	      for(i1 = 0; i1 < j; i1++)
		{
		  if(itest == 0 && z[iord[i1]] < z[ihld])
		    {
		      itest = 1;
		      for(i2 = j-1; i2 >= i1; i2--)
			{
			  iord[i2+1] = iord[i2];
			}
		      iord[i1] = ihld;
		    }
		}
	    }
	  /* MODIFICATION PART STARTS (Tapan July 2008)
	     iord[0] is the cell with highest ADC in the crude-cluster
	     ig is the number of local maxima in the crude-cluster
	     For the higest peak we make ig=0 which means first local maximum.
	     Next we go down in terms of the ADC sequence and find out if any
	     more of the cells form local maxima. The definition of local 
	     maxima is that all its neighbours are of less ADC compared to it. 
	  */
	  ig = 0;
	  xc[ig] = x[iord[0]];
	  yc[ig] = y[iord[0]];
	  zc[ig] = z[iord[0]];
	  tc[ig] = t[iord[0]];
	  Int_t ivalid = 0,  icount = 0;
	  
	  for(j=1;j<=ncl[i];j++)
	    {
	      x1  = x[iord[j]];
	      y1  = y[iord[j]]; 
	      z1  = z[iord[j]];
	      t1  = t[iord[j]];
	      rr=Distance(x1,y1,xc[ig],yc[ig]);
	      
	      // Check the cells which are outside the neighbours (rr>1.2)
	      if(rr>1.2 ) 
		{
		  ivalid=0;
		  icount=0;
		  for(Int_t j1=1;j1<j;j1++)
		    {
		      icount++;
		      Float_t rr1=Distance(x1,y1,x[iord[j1]],y[iord[j1]]);		
		      if(rr1>1.2) ivalid++;
		    }
		  if(ivalid == icount && z1>0.5*zc[ig])
		    {
		      ig++;
		      xc[ig]=x1; 
		      yc[ig]=y1; 
		      zc[ig]=z1;
		      tc[ig]=t1;
		    }
		}	  
	    }
	  
	  icl=icl+ig+1;
	  
	  //  We use simple Gaussian weighting. (Tapan Jan 2005)
	  // compute the number of cells belonging to each cluster.
	  // cell can be shared between several clusters  
	  // in the ratio of cluster energy deposition
	  // To calculate: 
	  // (1) number of cells belonging to a cluster (ig) and 
	  // (2) total ADC of the cluster (ig) 
	  // (3) x and y positions of the cluster
	  
	  
	  Int_t *cellCount;
	  Int_t **cellXY;
	  
	  Int_t    *status;
	  Double_t *totaladc, *totaladc2, *ncell,*weight;
	  Double_t *xclust, *yclust, *sigxclust, *sigyclust;
	  Double_t *ax, *ay, *ax2, *ay2;
	  
	  
	  status    = new Int_t [ncl[i]+1];
	  cellXY    = new Int_t *[ncl[i]+1];
	  
	  cellCount = new Int_t [ig+1];
	  totaladc  = new Double_t [ig+1];
	  totaladc2 = new Double_t [ig+1];
	  ncell     = new Double_t [ig+1];
	  weight    = new Double_t [ig+1];
	  xclust    = new Double_t [ig+1];
	  yclust    = new Double_t [ig+1];
	  sigxclust = new Double_t [ig+1];
	  sigyclust = new Double_t [ig+1];
	  ax        = new Double_t [ig+1];
	  ay        = new Double_t [ig+1];
	  ax2       = new Double_t [ig+1];
	  ay2       = new Double_t [ig+1];
	  
	  for(j = 0; j < ncl[i]+1; j++) 
	    {
	      status[j]     = 0;
	      cellXY[j] = new Int_t[ig+1];
	    }
	  //initialization
	  for(Int_t kcl = 0; kcl < ig+1; kcl++)
	    {
	      cellCount[kcl] = 0;
	      totaladc[kcl]  = 0.;
	      totaladc2[kcl] = 0.;
	      ncell[kcl]     = 0.;
	      weight[kcl]    = 0.;	
	      xclust[kcl]    = 0.; 
	      yclust[kcl]    = 0.;
	      sigxclust[kcl] = 0.; 
	      sigyclust[kcl] = 0.;
	      ax[kcl]        = 0.;      
	      ay[kcl]        = 0.;      
	      ax2[kcl]       = 0.;      
	      ay2[kcl]       = 0.;    
	      for(j = 0; j < ncl[i]+1; j++)
		{
		  cellXY[j][kcl] = 0;
		}
	    }
	  Double_t sumweight, gweight; 
	  
	  for(j = 0;j <= ncl[i]; j++)
	    {
	      x1 = x[iord[j]];
	      y1 = y[iord[j]];
	      z1 = z[iord[j]];
	      t1 = t[iord[j]];
	      
	      for(Int_t kcl=0; kcl<=ig; kcl++)
		{
		  x2 = xc[kcl];
		  y2 = yc[kcl];
		  rr = Distance(x1,y1,x2,y2);
		  t2 = tc[kcl];		  
		  
		  if(rr==0)
		    {
		      ncell[kcl]     = 1.;
		      totaladc[kcl]  = z1;
		      totaladc2[kcl] = z1*z1;
		      ax[kcl]        = x1 * z1;
		      ay[kcl]        = y1 * z1;
		      ax2[kcl]       = 0.;
		      ay2[kcl]       = 0.;
		      status[j]      = 1;
		    }
		}
	    }
	  
	  for(j = 0; j <= ncl[i]; j++)
	    {
	      Int_t   maxweight = 0;
	      Double_t     max  = 0.;
	      
	      if(status[j] == 0)
		{ 
		  x1 = x[iord[j]]; 
		  y1 = y[iord[j]];
		  z1 = z[iord[j]];
		  t1 = t[iord[j]];
		  sumweight = 0.;

		  for(Int_t kcl = 0; kcl <= ig; kcl++)
		    {
		      x2 = xc[kcl]; 
		      y2 = yc[kcl]; 
		      rr = Distance(x1,y1,x2,y2);
		      gweight     = exp(-(rr*rr)/(2*(1.2*1.2)));
		      weight[kcl] = zc[kcl] * gweight;
		      sumweight   = sumweight + weight[kcl];
		      
		      if(weight[kcl] > max)
			{
			  max       =  weight[kcl];
			  maxweight =  kcl;
			}
		    }
		  
		  cellXY[cellCount[maxweight]][maxweight] = iord[j];
		  		  
		  cellCount[maxweight]++;
		  
		  x2 = xc[maxweight];
		  y2 = yc[maxweight];
		  totaladc[maxweight]  +=  z1;
		  ax[maxweight]        +=  x1*z1;
		  ay[maxweight]        +=  y1*z1;
		  totaladc2[maxweight] +=  z1*z1;
		  ax2[maxweight]       +=  z1*(x1-x2)*(x1-x2);
		  ay2[maxweight]       +=  z1*(y1-y2)*(y1-y2);
		  ncell[maxweight]++;

		}
	    }
	  
	  for(Int_t kcl = 0; kcl <= ig; kcl++)
	    {
	      
	      if(totaladc[kcl] > 0.)
		{
		  xclust[kcl] = (ax[kcl])/ totaladc[kcl];
		  yclust[kcl] = (ay[kcl])/ totaladc[kcl];
		  
		  //natasha
		  Float_t sqtotadc = totaladc[kcl]*totaladc[kcl];
		  if(totaladc2[kcl] >= sqtotadc)
		    {
		      sigxclust[kcl] = 0.25;
		      sigyclust[kcl] = 0.25;
		    }
		  else
		    {
		      sigxclust[kcl] = (totaladc[kcl]/(sqtotadc-totaladc2[kcl]))*ax2[kcl];
		      sigyclust[kcl] = (totaladc[kcl]/(sqtotadc-totaladc2[kcl]))*ay2[kcl];
		    }	
		}
	      
	      for(j = 0; j < cellCount[kcl]; j++) clno++; 
	      
	      if (clno >= 4608) 
		{
		  AliWarning("RefClust: Too many clusters! more than 4608");
		  return;
		}
	      clusdata[0] = xclust[kcl];
	      clusdata[1] = yclust[kcl];
	      clusdata[2] = totaladc[kcl];
	      clusdata[3] = ncell[kcl];
	      
	      
	      if(sigxclust[kcl] > sigyclust[kcl]) 
		{
		  clusdata[4] = TMath::Sqrt(sigxclust[kcl]);
		  clusdata[5] = TMath::Sqrt(sigyclust[kcl]);
		}
	      else
		{
		  clusdata[4] = TMath::Sqrt(sigyclust[kcl]);
		  clusdata[5] = TMath::Sqrt(sigxclust[kcl]);
		}
	      
	      clxy[0] = tc[kcl];
	      
	      Int_t Ncell=1;
	      for (Int_t ii = 0; ii < cellCount[kcl]; ii++)
		{
		  if(ii<18) 
		    {	
		      clxy[Ncell] = t[cellXY[ii][kcl]];
		      Ncell++;
		    }
		} 
	      
	      pmdcludata = new AliPMDcludata(clusdata,clxy);
	      fPMDclucont->Add(pmdcludata);
	    }
	  
	  delete [] iord;
	  delete [] tc;	  
	  delete [] t;
	  delete [] x;
	  delete [] y;
	  delete [] z;
	  delete [] xc;
	  delete [] yc;
	  delete [] zc;
	  
	  delete [] cellCount;
	  for(Int_t jj = 0; jj < ncl[i]+1; jj++) delete [] cellXY[jj];
	  
	  delete [] status;
	  delete [] totaladc;
	  delete [] totaladc2;
	  delete [] ncell;
	  delete [] xclust;
	  delete [] yclust;
	  delete [] sigxclust;
	  delete [] sigyclust;
	  delete [] ax;
	  delete [] ay;
	  delete [] ax2;
	  delete [] ay2;
	  delete [] weight;
	}
    }
  delete [] ncl;
  delete [] clxy;
}
// ------------------------------------------------------------------------ //
Double_t AliPMDClusteringV1::Distance(Double_t x1, Double_t y1, 
				      Double_t x2, Double_t y2)
{
  return TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::FindIsoCell(Int_t idet, Int_t ismn, Double_t celladc[][96], TObjArray *pmdisocell)
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
		  if (id1 < 0) id1 = 0;
		  if (id1 > kMaxRow-1) id1 = kMaxRow - 1;
		  if (jd1 < 0) jd1 = 0;
		  if (jd1 > kMaxCol-1) jd1 = kMaxCol - 1;
		  Float_t adc = (Float_t) celladc[id1][jd1];
		  if(adc < 1.)
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
void AliPMDClusteringV1::SetEdepCut(Float_t decut)
{
  fCutoff = decut;
}
// ------------------------------------------------------------------------ //
