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
#include <TBRIK.h>
#include <TNode.h>
#include <TTree.h>
#include <TGeometry.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliPMD.h"
#include "AliDetector.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliHeader.h"

#include "AliPMDdigit.h"
#include "AliPMDClusterFinder.h"
#include "AliPMDClustering.h"
#include "AliPMDcluster.h"
#include "AliPMDrecpoint1.h"


ClassImp(AliPMDClusterFinder)

AliPMDClusterFinder::AliPMDClusterFinder()
{
//
// Default Constructor
//
  if (!fRecpoints) fRecpoints = new TClonesArray("AliPMDrecpoint1", 1000);  
  fNpoint = 0;

  fDebug = 0;
  fEcut  = 0.;

}
AliPMDClusterFinder::~AliPMDClusterFinder()
{
  // Destructor
  delete fRecpoints;
}

void AliPMDClusterFinder::OpengAliceFile(Char_t *file, Option_t *option)
{
  // Loads galice.root file and corresponding header, kinematics
  // hits and sdigits or digits depending on the option
  //
  fRunLoader = AliRunLoader::Open(file,AliConfig::fgkDefaultEventFolderName,
				  "UPDATE");
  
  if (!fRunLoader)
   {
     Error("Open","Can not open session for file %s.",file);
   }
  
  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice)
    {
      printf("<AliPMDdigitizer::Open> ");
      printf("AliRun object found on file.\n");
    }
  else
    {
      printf("<AliPMDdigitizer::Open> ");
      printf("Could not find AliRun object.\n");
    }
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  if (fPMDLoader == 0x0)
    {
      cerr<<"OpengAlice : Can not find PMD or PMDLoader\n";
    }

  const char *cDR = strstr(option,"DR");

  if (cDR)
    {
      fPMDLoader->LoadDigits("READ");
      fPMDLoader->LoadRecPoints("recreate");
    }
}

void AliPMDClusterFinder::Digits2RecPoints(Int_t ievt)
{
  // Converts digits to recpoints after running clustering
  // algorithm on CPV plane and PREshower plane
  //
  Int_t    det = 0,smn = 0;
  Int_t    xpos,ypos;
  Float_t  adc;
  Int_t    isup;
  Int_t    idet;
  Float_t  clusdata[7];

  TObjArray *pmdcont = new TObjArray();
  AliPMDcluster  *pmdcl  = new AliPMDcluster;
  AliPMDClustering *pmdclust = new AliPMDClustering();
  pmdclust->SetDebug(fDebug);
  pmdclust->SetEdepCut(fEcut);

  fRunLoader->GetEvent(ievt);
  //cout << " ***** Beginning::Digits2RecPoints *****" << endl;
  fTreeD = fPMDLoader->TreeD();
  if (fTreeD == 0x0)
    {
      cout << " Can not get TreeD" << endl;
    }
  AliPMDdigit  *pmddigit;
  TBranch *branch = fTreeD->GetBranch("PMDDigit");
  branch->SetAddress(&fDigits);

  ResetRecpoint();
  fTreeR = fPMDLoader->TreeR();
  if (fTreeR == 0x0)
    {
      fPMDLoader->MakeTree("R");
      fTreeR = fPMDLoader->TreeR();
    }

  Int_t bufsize = 16000;
  fTreeR->Branch("PMDRecpoint", &fRecpoints, bufsize); 

  Int_t nmodules = (Int_t) fTreeD->GetEntries();
  
  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      ResetCellADC();
      fTreeD->GetEntry(imodule); 
      Int_t nentries = fDigits->GetLast();
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmddigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
	  
	  det    = pmddigit->GetDetector();
	  smn    = pmddigit->GetSMNumber();
	  xpos   = pmddigit->GetRow();
	  ypos   = pmddigit->GetColumn();
	  adc    = pmddigit->GetADC();
	  //Int_t trno   = pmddigit->GetTrackNumber();

	  fCellADC[xpos][ypos] = (Double_t) adc;
	}

      idet = det;
      isup = smn;
      pmdclust->DoClust(fCellADC,pmdcont);
      
      Int_t nentries1 = pmdcont->GetEntries();
      cout << " nentries1 = " << nentries1 << endl;
      for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
	{
	  clusdata[0] = (Float_t) idet;
	  clusdata[1] = (Float_t) isup;
	      
	  pmdcl = (AliPMDcluster*)pmdcont->UncheckedAt(ient1);
	      
	  clusdata[2] = pmdcl->GetClusX();
	  clusdata[3] = pmdcl->GetClusY();
	  clusdata[4] = pmdcl->GetClusADC();
	  clusdata[5] = pmdcl->GetClusCells();
	  clusdata[6] = pmdcl->GetClusRadius();
	  
	  AddRecPoint(clusdata);
	}
      pmdcont->Clear();
      
      fTreeR->Fill();
      ResetRecpoint();

    } // modules

  ResetCellADC();
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");  
  fPMDLoader->WriteRecPoints("OVERWRITE");

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
    
  //  cout << " ***** End::Digits2RecPoints *****" << endl;
}

void AliPMDClusterFinder::SetCellEdepCut(Float_t ecut)
{
  fEcut = ecut;
}
void AliPMDClusterFinder::SetDebug(Int_t idebug)
{
  fDebug = idebug;
}

void AliPMDClusterFinder::AddRecPoint(Float_t *clusdata)
{
  // Add Reconstructed points
  //
  TClonesArray &lrecpoints = *fRecpoints;
  AliPMDrecpoint1 *newrecpoint;
  newrecpoint = new AliPMDrecpoint1(clusdata);
  new(lrecpoints[fNpoint++]) AliPMDrecpoint1(newrecpoint);
  delete newrecpoint;
}
void AliPMDClusterFinder::ResetCellADC()
{
  // Reset the individual cell ADC value to zero
  //
  for(Int_t irow = 0; irow < fgkRow; irow++)
    {
      for(Int_t icol = 0; icol < fgkCol; icol++)
	{
	  fCellADC[irow][icol] = 0.;
	}
    }
}

void AliPMDClusterFinder::ResetRecpoint()
{
  // Clear the list of reconstructed points
  fNpoint = 0;
  if (fRecpoints) fRecpoints->Clear();
}
void AliPMDClusterFinder::UnLoad(Option_t *option)
{
  // Unload all the *.root files
  //
  const char *cR = strstr(option,"R");

  fRunLoader->UnloadgAlice();

  if (cR)
    {
      fPMDLoader->UnloadDigits();
      fPMDLoader->UnloadRecPoints();
    }
}
