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
				 Double_t celladc[48][96], TObjArray *pmdcont)
{
  // main function to call other necessary functions to do clustering
  //

  AliPMDcluster *pmdcl = 0;

  Int_t    i,  j, nmx1, incr, id, jd;
  Int_t    celldataX[15], celldataY[15];
  Float_t  clusdata[6];
  Double_t cutoff, ave;
  Double_t edepcell[kNMX];
  
  
  Double_t *cellenergy = new Double_t [11424];// Ajay
  
  const float ktwobysqrt3 = 1.1547; // 2./sqrt(3.)

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
	  fCellTrNo[i][j] = -1;
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
	      cellenergy[ij]    = celladc[jd][id];//Ajay
	      fCellTrNo[i][j]   = jd*10000+id;    // for association 
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      cellenergy[ij]    = celladc[id][jd];//Ajay
	      fCellTrNo[i][j] = id*10000+jd;      // for association 
	    }
	}
    }
  
  for (i = 0; i < kNMX; i++)
  {
    edepcell[i] = cellenergy[i];
  }

  delete [] cellenergy;

  Int_t iord1[kNMX];
  TMath::Sort((Int_t)kNMX,edepcell,iord1);// order the data
  cutoff = fCutoff;                // cutoff to discard cells
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

      for (Int_t ihit = 0; ihit < 15; ihit++)
	{
	  if (ismn < 12)
	    {
	      celldataX[ihit] = pmdcludata->GetCellXY(ihit)%10000;
	      celldataY[ihit] = pmdcludata->GetCellXY(ihit)/10000;
	    }
	  else if (ismn >= 12 && ismn <= 23)
	    {
	      celldataX[ihit] = pmdcludata->GetCellXY(ihit)/10000;
	      celldataY[ihit] = pmdcludata->GetCellXY(ihit)%10000;
	    }
	}
      pmdcl = new AliPMDcluster(idet, ismn, clusdata, celldataX, celldataY);
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
  Int_t i,j,k,id1,id2,icl, numcell, clust[2][5000];
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
	  
	  for(i = 1; i < 5000; i++)
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
	  for(i = 1; i < 5000;i++)
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

  Int_t *cellCount = 0x0;
  Int_t **cellXY = 0x0;
  const Int_t kdim = 4609;

  Int_t    i12;
  Int_t    i, j, k, i1, i2, id, icl,  itest;
//  Int_t    ihld;
  Int_t    ig, nsupcl,clno;
  Int_t    t[kdim];
  Int_t    ncl[kdim], iord[kdim], lev1[kdim], lev2[kdim];
  Int_t    clxy[15];
  Float_t  clusdata[6];
  Double_t x1, y1, z1, x2, y2, z2, dist,rr,sum;
  Double_t x[kdim], y[kdim], z[kdim];
  Double_t xc[kdim], yc[kdim], zc[kdim], cells[kdim], rc[kdim];

  // Initialisation  
  for(i = 0; i<kdim; i++)
    { 
      t[i]         = -1;
      ncl[i]       = -1;
      if (i < 6) clusdata[i] = 0.;
      if (i < 15) clxy[i] = 0;
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
      if (nsupcl > kdim) 
	{
	  AliWarning("RefClust: Too many superclusters!");
	  nsupcl = kdim;
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
	  if (clno >= 5000) 
	    {
	      AliWarning("RefClust: Too many clusters! more than 5000");
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

	  clxy[0] = fCellTrNo[i1][i2];      	  //association
	  for(Int_t icltr = 1; icltr < 15; icltr++)
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
	  clxy[0] = fCellTrNo[i1][i2];	  //asso
	  id++;
	  i1   = fInfcl[1][id];
	  i2   = fInfcl[2][id];

	  Int_t i22 = i1 + i2*kNDIMX;
	  x2   = fCoord[0][i1][i2];
	  y2   = fCoord[1][i1][i2];
	  z2   = edepcell[i22];
	  clxy[1] = fCellTrNo[i1][i2];  	  //asso
	  for(Int_t icltr = 2; icltr < 15; icltr++)
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
	  id++;
	  iord[0] = 0;
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

	  t[0] = fCellTrNo[i1][i2];	  //asso
	  
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

	      t[j]    = fCellTrNo[i1][i2];	      //asso
	    }
	  
	  // arranging cells within supercluster in decreasing order

/*	  
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
*/

	  Int_t imaxdim = ncl[i] + 1;
	  TMath::Sort(imaxdim,z,iord);// order the data


	  // compute the number of Gaussians and their centers ( first
	  // guess )
	  // centers must be separated by cells having smaller ener. dep.
	  // neighbouring centers should be either strong or well-separated
	  ig = 0;
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
		  if(rr >= 1.1 && rr < 1.8 && z[iord[j]] > zc[k]/4.)itest++;
		  if(rr >= 1.8 && rr < 2.1 && z[iord[j]] > zc[k]/10.)itest++;
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
	  GaussFit(ncl[i], ig, x[0], y[0] ,z[0], xc[0], yc[0], zc[0], rc[0]);
	  icl += ig+1;
	  // compute the number of cells belonging to each cluster.
	  // cell is shared between several clusters ( if they are equidistant
	  // from it ) in the ratio of cluster energy deposition

	  Int_t jj = 15;
	  cellCount = new Int_t [ig+1];
	  cellXY = new Int_t *[jj];
	  for(Int_t ij = 0; ij < 15; ij++) cellXY[ij] = new Int_t [ig+1];

	  for(j = 0; j <= ig; j++)
	    {
	      cellCount[j] = 0;
	      cells[j]     = 0.;
	    }
	  
	  if(ig > 0)
	    {
	      for(j = 0; j <= ncl[i]; j++)
		{
		  lev1[j] = 0;
		  lev2[j] = 0;
		  for(k = 0; k <= ig; k++)
		    {
		      dist = Distance(x[j], y[j], xc[k], yc[k]);
		      if(dist < TMath::Sqrt(3.) )
			{
			  //asso
			  if (cellCount[k] < 15)
			    {
			      cellXY[cellCount[k]][k] = t[j];
			    }
			  cellCount[k]++;
			  //
			  lev1[0]++;
			  i1       = lev1[0];
			  lev1[i1] = k;
			}
		      else
			{
			  if(dist < 2.1)
			    {
			      lev2[0]++;
			      i1       = lev2[0];
			      lev2[i1] = k;
			    }
			}
		    }
		  if(lev1[0] != 0)
		    {
		      if(lev1[0] == 1)
			{
			  cells[lev1[1]]++;
			} 
		      else 
			{
			  sum=0.;
			  for(k = 1; k <= lev1[0]; k++)
			    {
			      sum  += zc[lev1[k]];
			    }
			  for(k = 1; k <= lev1[0]; k++)
			    {
			      cells[lev1[k]] += zc[lev1[k]]/sum;
			    }
			}
		    }
		  else
		    {
		      if(lev2[0] == 0)
			{
			  cells[lev2[1]]++;
			}
		      else
			{
			  sum=0.;
			  for( k = 1; k <= lev2[0]; k++)
			    {
			      sum += zc[lev2[k]];
			    }
			  for(k = 1; k <= lev2[0]; k++)
			    {
			      cells[lev2[k]] +=  zc[lev2[k]]/sum;
			    }
			}
		    }
		}
	    }
	  
	  // zero rest of the cell array
	  //asso
	  for( k = 0; k <= ig; k++)
	    {
	      for(Int_t icltr = cellCount[k]; icltr < 15; icltr++)
		{
		  cellXY[icltr][k] = -1;
		}
	    }
	  //
	  
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
	      clusdata[4] = rc[j];
	      clusdata[5] = 0.0;
	      if(ig == 0)
		{
		  clusdata[3] = ncl[i];
		}
	      else
		{
		  clusdata[3] = cells[j];
		}


	      for (Int_t ii=0; ii < 15; ii++)
		{
		  clxy[ii] = cellXY[ii][j];
		}  
	      pmdcludata = new AliPMDcludata(clusdata,clxy);
	      fPMDclucont->Add(pmdcludata);
	    }
	  delete [] cellCount;
	  for(jj = 0; jj < 15; jj++)delete [] cellXY[jj];
	  delete [] cellXY;
	}
    }
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::GaussFit(Int_t ncell, Int_t nclust, Double_t &x, 
				  Double_t &y ,Double_t &z, Double_t &xc, 
				  Double_t &yc, Double_t &zc, Double_t &rc)
{
  // Does gaussian fitting
  //

  const Int_t kdim = 4609;
  Int_t i, j, i1, i2, novar, idd, jj;
  Int_t neib[kdim][50];

  Double_t sum, dx, dy, str, str1, aint, sum1, rr, dum;
  Double_t x1, x2, y1, y2;
  Double_t xx[kdim], yy[kdim], zz[kdim], xxc[kdim], yyc[kdim];
  Double_t a[kdim], b[kdim], c[kdim], d[kdim], ha[kdim], hb[kdim];
  Double_t hc[kdim], hd[kdim], zzc[kdim], rrc[kdim];
  
  TRandom rnd;
  
  str   = 0.;
  str1  = 0.;
  rr    = 0.3;
  novar = 0;
  j     = 0;  

  for(i = 0; i <= ncell; i++)
    {
      xx[i] = *(&x+i);
      yy[i] = *(&y+i);
      zz[i] = *(&z+i);
      str  += zz[i];
    }
  for(i=0; i<=nclust; i++)
    {
      xxc[i] = *(&xc+i);
      yyc[i] = *(&yc+i);
      zzc[i] = *(&zc+i);
      str1  += zzc[i];
      rrc[i] = 0.5;
    }
  for(i = 0; i <= nclust; i++)
    {
      zzc[i] = str/str1*zzc[i];
      ha[i]  = xxc[i];
      hb[i]  = yyc[i];
      hc[i]  = zzc[i];
      hd[i]  = rrc[i];
      x1     = xxc[i];
      y1     = yyc[i];
    }
 
  for(i = 0; i <= ncell; i++)
    {
      idd = 0;
      x1  = xx[i];
      y1  = yy[i];
      for(j = 0; j <= nclust; j++)
	{
	  x2 = xxc[j];
	  y2 = yyc[j];
	  if(Distance(x1,y1,x2,y2) <= 3.)
	    { 
	      idd++;
	      neib[i][idd] = j; 
	    }
	}
      neib[i][0] = idd;
    }
  sum = 0.;
  for(i1 = 0; i1 <= ncell; i1++)
    {
      aint = 0.;
      idd = neib[i1][0];
      for(i2 = 1; i2 <= idd; i2++)
	{
	  jj    = neib[i1][i2];
	  dx    = xx[i1] - xxc[jj];
	  dy    = yy[i1] - yyc[jj];
	  dum   = rrc[j]*rrc[jj] + rr*rr;
	  aint += exp(-(dx*dx+dy*dy)/dum)*zzc[idd]*rr*rr/dum;
	}
      sum += (aint - zz[i1])*(aint - zz[i1])/str;
    } 
  str1 = 0.;
 
  for(i = 0; i <= nclust; i++)
    {
      a[i]  = xxc[i] + 0.6*(rnd.Uniform() - 0.5);
      b[i]  = yyc[i] + 0.6*(rnd.Uniform() - 0.5);
      c[i]  = zzc[i]*(1.+ ( rnd.Uniform() - 0.5)*0.2);
      str1 += zzc[i];
      d[i]  = rrc[i]*(1.+ ( rnd.Uniform() - 0.5)*0.1);
      
      if(d[i] < 0.25)
	{
	  d[i]=0.25;
	}
    }
  for(i = 0; i <= nclust; i++)
    {
      c[i] = c[i]*str/str1; 
    }
  sum1=0.;
  
  for(i1 = 0; i1 <= ncell; i1++)
    {
      aint = 0.;
      idd = neib[i1][0];
      for(i2 = 1; i2 <= idd; i2++)
	{
	  jj    = neib[i1][i2];
	  dx    = xx[i1] - a[jj];
	  dy    = yy[i1] - b[jj];
	  dum   = d[jj]*d[jj]+rr*rr;
	  aint += exp(-(dx*dx+dy*dy)/dum)*c[i2]*rr*rr/dum;
	}
      sum1 += (aint - zz[i1])*(aint - zz[i1])/str;
    }

    if(sum1 < sum)
      {
	for(i2 = 0; i2 <= nclust; i2++)
	{
	  xxc[i2] = a[i2];
	  yyc[i2] = b[i2];
	  zzc[i2] = c[i2];
	  rrc[i2] = d[i2];
	  sum     = sum1;
	}
      }
    for(j = 0; j <= nclust; j++)
      {
	*(&xc+j) = xxc[j];
	*(&yc+j) = yyc[j];
	*(&zc+j) = zzc[j];
	*(&rc+j) = rrc[j];
      }
}
// ------------------------------------------------------------------------ //
Double_t AliPMDClusteringV1::Distance(Double_t x1, Double_t y1, 
				      Double_t x2, Double_t y2)
{
  return TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
// ------------------------------------------------------------------------ //
void AliPMDClusteringV1::SetEdepCut(Float_t decut)
{
  fCutoff = decut;
}
// ------------------------------------------------------------------------ //
