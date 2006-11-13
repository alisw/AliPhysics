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
//           Date   : March 25 2004                    //
//  This reads the file PMD.RecPoints.root(TreeR),     //
//  calls the Clustering algorithm and stores the      //
//  clustering output in PMD.RecPoints.root(TreeR)     // 
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <TMath.h>
#include <TBRIK.h>
#include <TNode.h>
#include <TTree.h>
#include <TGeometry.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TBranch.h>
#include <TNtuple.h>
#include <TParticle.h>

#include "AliPMDcluster.h"
#include "AliPMDclupid.h"
#include "AliPMDrecpoint1.h"
#include "AliPMDUtility.h"
#include "AliPMDDiscriminator.h"
#include "AliPMDEmpDiscriminator.h"
#include "AliPMDtracker.h"

#include "AliESDPmdTrack.h"
#include "AliESD.h"
#include "AliLog.h"

ClassImp(AliPMDtracker)

AliPMDtracker::AliPMDtracker():
  fTreeR(0),
  fRecpoints(new TClonesArray("AliPMDrecpoint1", 1000)),
  fPMDcontin(new TObjArray()),
  fPMDcontout(new TObjArray()),
  fPMDutil(new AliPMDUtility()),
  fPMDrecpoint(0),
  fPMDclin(0),
  fPMDclout(0),
  fXvertex(0.),
  fYvertex(0.),
  fZvertex(0.),
  fSigmaX(0.),
  fSigmaY(0.),
  fSigmaZ(0.)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------//
AliPMDtracker:: AliPMDtracker(const AliPMDtracker & /* tracker */):
  TObject(/* tracker */),
  fTreeR(0),
  fRecpoints(NULL),
  fPMDcontin(NULL),
  fPMDcontout(NULL),
  fPMDutil(NULL),
  fPMDrecpoint(0),
  fPMDclin(0),
  fPMDclout(0),
  fXvertex(0.),
  fYvertex(0.),
  fZvertex(0.),
  fSigmaX(0.),
  fSigmaY(0.),
  fSigmaZ(0.)
{
  // copy constructor
  AliError("Copy constructor not allowed");
}

//--------------------------------------------------------------------//
AliPMDtracker& AliPMDtracker::operator=(const AliPMDtracker & /* tracker */)
{
 // assignment operator
  AliError("Assignment operator not allowed");
  return *this;
}

//--------------------------------------------------------------------//
AliPMDtracker::~AliPMDtracker()
{
  // Destructor
  if (fRecpoints)
    {
      fRecpoints->Delete();
      delete fRecpoints;
      fRecpoints=0;
    }
  if (fPMDcontin)
    {
      fPMDcontin->Delete();
      delete fPMDcontin;
      fPMDcontin=0;
    }
  if (fPMDcontout)
    {
      fPMDcontout->Delete();
      delete fPMDcontout;
      fPMDcontout=0;
    }
}
//--------------------------------------------------------------------//
void AliPMDtracker::LoadClusters(TTree *treein)
{
  // Load the Reconstructed tree
  fTreeR = treein;
}
//--------------------------------------------------------------------//
void AliPMDtracker::Clusters2Tracks(AliESD *event)
{
  // Converts digits to recpoints after running clustering
  // algorithm on CPV plane and PREshower plane
  //

  Int_t   idet;
  Int_t   ismn;
  Float_t clusdata[6];

  TBranch *branch = fTreeR->GetBranch("PMDRecpoint");
  if (!branch)
    {
      AliError("PMDRecpoint branch not found");
      return;
    }
  branch->SetAddress(&fRecpoints);  
  
  Int_t   nmodules = (Int_t) branch->GetEntries();
  
  AliDebug(1,Form("Number of modules filled in treeR = %d",nmodules));
  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      branch->GetEntry(imodule); 
      Int_t nentries = fRecpoints->GetLast();
      AliDebug(2,Form("Number of clusters per modules filled in treeR = %d"
		      ,nentries));
      for(Int_t ient = 0; ient < nentries+1; ient++)
	{
	  fPMDrecpoint = (AliPMDrecpoint1*)fRecpoints->UncheckedAt(ient);
	  idet        = fPMDrecpoint->GetDetector();
	  ismn        = fPMDrecpoint->GetSMNumber();
	  clusdata[0] = fPMDrecpoint->GetClusX();
	  clusdata[1] = fPMDrecpoint->GetClusY();
	  clusdata[2] = fPMDrecpoint->GetClusADC();
	  clusdata[3] = fPMDrecpoint->GetClusCells();
	  clusdata[4] = fPMDrecpoint->GetClusSigmaX();
	  clusdata[5] = fPMDrecpoint->GetClusSigmaY();

	  fPMDclin = new AliPMDrecpoint1(idet,ismn,clusdata);
	  fPMDcontin->Add(fPMDclin);
	}
    }

  AliPMDDiscriminator *pmddiscriminator = new AliPMDEmpDiscriminator();
  pmddiscriminator->Discrimination(fPMDcontin,fPMDcontout);

  const Float_t kzpos = 361.5;    // middle of the PMD

  Int_t   det,smn;
  Float_t xpos,ypos;
  Float_t adc, ncell, rad;
  Float_t xglobal = 0., yglobal = 0., zglobal = 0;
  Float_t pid;


  Int_t nentries2 = fPMDcontout->GetEntries();
  AliDebug(1,Form("Number of clusters coming after discrimination = %d"
		  ,nentries2));
  for (Int_t ient1 = 0; ient1 < nentries2; ient1++)
    {
      fPMDclout = (AliPMDclupid*)fPMDcontout->UncheckedAt(ient1);
      
      det   = fPMDclout->GetDetector();
      smn   = fPMDclout->GetSMN();
      xpos  = fPMDclout->GetClusX();
      ypos  = fPMDclout->GetClusY();
      adc   = fPMDclout->GetClusADC();
      ncell = fPMDclout->GetClusCells();
      rad   = fPMDclout->GetClusRadius();
      pid   = fPMDclout->GetClusPID();
      
      //
      /**********************************************************************
       *    det   : Detector, 0: PRE & 1:CPV                                *
       *    smn   : Serial Module Number 0 to 23 for each plane             *
       *    xpos  : x-position of the cluster                               *
       *    ypos  : y-position of the cluster                               *
       *            THESE xpos & ypos are not the true xpos and ypos        *
       *            for some of the unit modules. They are rotated.         *
       *    adc   : ADC contained in the cluster                            *
       *    ncell : Number of cells contained in the cluster                *
       *    rad   : radius of the cluster (1d fit)                          *
       **********************************************************************/
      //

      fPMDutil->RectGeomCellPos(smn,xpos,ypos,xglobal,yglobal);

      if (det == 0)
	{
	  zglobal = kzpos + 1.6; // PREshower plane
	}
      else if (det == 1)
	{
	  zglobal = kzpos - 1.7; // CPV plane
	}

      // Fill ESD

      AliESDPmdTrack *esdpmdtr = new  AliESDPmdTrack();

      esdpmdtr->SetDetector(det);
      esdpmdtr->SetClusterX(xglobal);
      esdpmdtr->SetClusterY(yglobal);
      esdpmdtr->SetClusterZ(zglobal);
      esdpmdtr->SetClusterADC(adc);
      esdpmdtr->SetClusterCells(ncell);
      esdpmdtr->SetClusterPID(pid);

      event->AddPmdTrack(esdpmdtr);
    }
}
//--------------------------------------------------------------------//
void AliPMDtracker::SetVertex(Double_t vtx[3], Double_t evtx[3])
{
  fXvertex = vtx[0];
  fYvertex = vtx[1];
  fZvertex = vtx[2];
  fSigmaX  = evtx[0];
  fSigmaY  = evtx[1];
  fSigmaZ  = evtx[2];
}
//--------------------------------------------------------------------//
void AliPMDtracker::ResetClusters()
{
  if (fRecpoints) fRecpoints->Clear();
}
//--------------------------------------------------------------------//
