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

//-----------------------------------------------------//
//                                                     //
//           Date   : August 05 2003                   //
//  This reads the file PMD.digits.root(TreeD),        //
//  calls the Clustering algorithm and stores the      //
//  clustering output in PMD.RecPoints.root(TreeR)     // 
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <TMath.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliPMDcluster.h"
#include "AliPMDclupid.h"
#include "AliPMDDiscriminator.h"

ClassImp(AliPMDDiscriminator)

AliPMDDiscriminator::AliPMDDiscriminator() :
  fDiscrim(0)
{
//
// Default Constructor
//

}
// -----------------------------------------------------------------------
AliPMDDiscriminator::~AliPMDDiscriminator()
{
  // Destructor
}

// -----------------------------------------------------------------------

void AliPMDDiscriminator::Discrimination(TObjArray *pmdcontin, TObjArray *pmdcontout)
{
  // Does Photon/Hadron discrimination

  if(fDiscrim == 0)
    {
      EmpDiscrimination(pmdcontin, pmdcontout);
    }
  else if(fDiscrim == 1)
    {
      NNDiscrimination();
    }
}
// -----------------------------------------------------------------------

void AliPMDDiscriminator::EmpDiscrimination(TObjArray *pmdcontin, TObjArray *pmdcontout)
{
  // Does Photon/Hadron discrimination
  // matching the clusters of CPV and PREshower plane
  //
  const  Int_t kumperdet = 24;
  static Int_t neibx[6]={1,0,-1,-1,0,1}, neiby[6]={0,1,1,0,-1,-1}; 
  Int_t   det,smn;
  Int_t   iprecount[24], icpvcount[24];
  Float_t xpos,ypos;
  Float_t adc, ncell, rad;
  Float_t clusdata[6];

  for(Int_t i = 0; i < kumperdet; i++)
    {
      iprecount[i] = 0;
      icpvcount[i] = 0;
    }
  AliPMDcluster  *pmdcl    = 0;
  AliPMDclupid   *pmdclout = 0;

  Int_t nentries1 = pmdcontin->GetEntries();

  AliDebug(1,Form("Number of total clusters from CPV PRE = %d",nentries1));
  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
    {
      pmdcl = (AliPMDcluster*)pmdcontin->UncheckedAt(ient1);

      det   = pmdcl->GetDetector();
      smn   = pmdcl->GetSMN();
      if(det == 0) iprecount[smn]++;
      if(det == 1) icpvcount[smn]++;
    } // Entries of TObjArray loop

  Int_t   idet, ismn;
  Float_t edepcpv[48][96];
  Int_t   statuscpv[48][96];

  for(Int_t i = 0; i < kumperdet; i++) // unit module
    {
      // Initialisation of the ADC of CPV (1UM)
      for (Int_t ix = 0; ix < 48;ix++)
	{
	  for (Int_t iy = 0; iy < 96;iy++)
	    {
	      edepcpv[ix][iy]   = 0.;
	      statuscpv[ix][iy] = 0;
	    }
	}
      Int_t precounter   = iprecount[i];
      Int_t cpvcounter   = icpvcount[i];

      Float_t *xpadpre   = new Float_t[precounter];
      Float_t *ypadpre   = new Float_t[precounter];
      Float_t *adcpre    = new Float_t[precounter];
      Float_t *ncellpre  = new Float_t[precounter];
      Float_t *radpre    = new Float_t[precounter];
      Int_t   *sortcoord = new Int_t[precounter];
      Int_t   *clupidpre = new Int_t[precounter];

      Float_t *xpadcpv   = new Float_t[cpvcounter];
      Float_t *ypadcpv   = new Float_t[cpvcounter];
      Float_t *adccpv    = new Float_t[cpvcounter];
      Float_t *ncellcpv  = new Float_t[cpvcounter];
      Float_t *radcpv    = new Float_t[cpvcounter];

      Int_t ii = 0;
      Int_t ij = 0;
      for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	{
	  pmdcl = (AliPMDcluster*)pmdcontin->UncheckedAt(ient1);
	  
	  det   = pmdcl->GetDetector();
	  smn   = pmdcl->GetSMN();
	  xpos  = pmdcl->GetClusX();
	  ypos  = pmdcl->GetClusY();
	  adc   = pmdcl->GetClusADC();
	  ncell = pmdcl->GetClusCells();
	  rad   = pmdcl->GetClusRadius();
	  
	  if(det == 0 && smn == i)
	    {
	      xpadpre[ii]  = xpos;
	      ypadpre[ii]  = ypos;
	      adcpre[ii]   = adc;
	      ncellpre[ii] = ncell;
	      radpre[ii]   = rad;
	      ii++;
	    }
	  if(det == 1 && smn == i)
	    {
	      Int_t ix = (Int_t) (xpos+0.5);
	      Int_t iy = (Int_t) (ypos+0.5);
	      if(ix > 47) ix = 47;
	      if(iy > 95) iy = 95;
	      edepcpv[ix][iy] = adc;
	      xpadcpv[ij]  = xpos;
	      ypadcpv[ij]  = ypos;
	      adccpv[ij]   = adc;
	      ncellcpv[ij] = ncell;
	      radcpv[ij]   = rad;
	      ij++;
	    }
	} // Entries of TObjArray loop
      // sorts from lowest ADC to highest ADC
      // and returns the coordinates
      Bool_t jsort = false;
      TMath::Sort(precounter,adcpre,sortcoord,jsort);

      Int_t jjsort = 0;
      for(Int_t jj=0; jj<precounter; jj++)
	{
	  // PRE information
	  // PIDs for PRE clusters are 0(photon) and 1(hadron)

	  jjsort = sortcoord[jj];

	  Int_t ix = (Int_t) (xpadpre[jjsort]+0.5);
	  Int_t iy = (Int_t) (ypadpre[jjsort]+0.5);
	  if(ix > 47) ix = 47;
	  if(iy > 95) iy = 95;

	  for(Int_t jk=0; jk<6; jk++)
	    {
	      Int_t jd1 = ix + neibx[jk]; 
	      Int_t jd2 = iy + neiby[jk];
	      if(jd1 <0 ) jd1 = 0;
	      if(jd1 >47) jd1 = 47;
	      if(jd2 <0 ) jd2 = 0;
	      if(jd2 >47) jd2 = 47;
	      if(edepcpv[jd1][jd2] > 0.0 && statuscpv[jd1][jd2] == 0)
		{
		  statuscpv[jd1][jd2] = 1;
		  clupidpre[jjsort]   = 1;
		  break;
		}
	    }

	  idet        = 0;
	  ismn        = i;
	  clusdata[0] = xpadpre[jjsort];
	  clusdata[1] = ypadpre[jjsort];
	  clusdata[2] = adcpre[jjsort];
	  clusdata[3] = ncellpre[jjsort];
	  clusdata[4] = radpre[jjsort];
	  //PH	  clusdata[5] = (Float_t) clupidpre[jjsort];

	  // Temporary the cluster PID is set to 1 if the
	  // adc > 3MIP units which will be changed later on.

	  if (adcpre[jjsort] > 4500.)
	    {
	      clusdata[5] = 1.0;
	    }
	  else
	    {
	      clusdata[5] = 0.0;
	    }
	  pmdclout = new AliPMDclupid(idet,ismn,clusdata);
	  pmdcontout->Add(pmdclout);
	} // jj loop

      for(Int_t jj=0; jj<cpvcounter; jj++)
	{
	  // CPV information
	  // PID for CPV clusters is 1

	  idet        = 1;
	  ismn        = i;
	  clusdata[0] = xpadcpv[jj];
	  clusdata[1] = ypadcpv[jj];
	  clusdata[2] = adccpv[jj];
	  clusdata[3] = ncellcpv[jj];
	  clusdata[4] = radcpv[jj];
	  clusdata[5] = 0.;

	  pmdclout = new AliPMDclupid(idet,ismn,clusdata);
	  pmdcontout->Add(pmdclout);
	}

      // delete all the pointers
      delete [] xpadpre;
      delete [] ypadpre;
      delete [] adcpre;
      delete [] ncellpre;
      delete [] radpre;
      delete [] clupidpre;
      delete [] xpadcpv;
      delete [] ypadcpv;
      delete [] adccpv;
      delete [] ncellcpv;
      delete [] radcpv;
    } // i loop

}
// -----------------------------------------------------------------------
void AliPMDDiscriminator::NNDiscrimination()
{
  // This method does discrimination using Neural Network technique

}
// -----------------------------------------------------------------------
void AliPMDDiscriminator::SetDiscrimination(Int_t idiscrim)
{
  fDiscrim = idiscrim;
}
// -----------------------------------------------------------------------

