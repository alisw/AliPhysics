//-----------------------------------------------------//
//                                                     //
//           Date   : August 05 2003                   //
//  This reads the file PMD.digits.root(TreeD),        //
//  calls the Clustering algorithm and stores the      //
//  clustering output in PMD.RecPoints.root(TreeR)     // 
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
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
#include "AliPMDrecpoint.h"


ClassImp(AliPMDClusterFinder)
//
// Constructor
//
AliPMDClusterFinder::AliPMDClusterFinder()
{
  if (!fRecpoints) fRecpoints = new TClonesArray("AliPMDrecpoint", 1000);  
  fNpoint = 0;

  for (Int_t i = 0; i < fTotSM; i++)
    {
      for (Int_t j = 0; j < fNCell; j++)
	{
	  for (Int_t k = 0; k < fNCell; k++)
	    {
	      fCPV[i][j][k] = 0.; 
	      fPMD[i][j][k] = 0.; 
	    }
	}
    }

}
AliPMDClusterFinder::~AliPMDClusterFinder()
{
  delete fRecpoints;
}
//
// Member functions
//
void AliPMDClusterFinder::OpengAliceFile(Char_t *file, Option_t *option)
{

  fRunLoader = AliRunLoader::Open(file,AliConfig::fgkDefaultEventFolderName,
				  "UPDATE");
  
  if (!fRunLoader)
   {
     Error("Open","Can not open session for file %s.",file);
   }
  
  fRunLoader->LoadgAlice();
  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();

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
  PMD  = (AliPMD*)gAlice->GetDetector("PMD");
  pmdloader = fRunLoader->GetLoader("PMDLoader");
  if (pmdloader == 0x0)
    {
      cerr<<"OpengAlice : Can not find PMD or PMDLoader\n";
    }

  const char *cDR = strstr(option,"DR");

  if (cDR)
    {
      pmdloader->LoadDigits("READ");
      pmdloader->LoadRecPoints("recreate");
    }
}

void AliPMDClusterFinder::Digits2RecPoints(Int_t ievt)
{
  Int_t    det,smn;
  Int_t    cellno;
  Int_t    xpos,ypos;
  Float_t  adc;
  Int_t    isup, ix, iy;
  Int_t    idet;
  Double_t d[72][72];
  Float_t  clusdata[7];
  Int_t    fMessage = 1;

  fRunLoader->GetEvent(ievt);
  //cout << " ***** Beginning::Digits2RecPoints *****" << endl;
  treeD = pmdloader->TreeD();
  if (treeD == 0x0)
    {
      cout << " Can not get TreeD" << endl;
    }
  AliPMDdigit  *pmddigit;
  TBranch *branch = treeD->GetBranch("PMDDigit");
  branch->SetAddress(&fDigits);

  ResetRecpoint();
  treeR = pmdloader->TreeR();
  if (treeR == 0x0)
    {
      pmdloader->MakeTree("R");
      treeR = pmdloader->TreeR();
    }

  Int_t bufsize = 16000;
  treeR->Branch("PMDRecpoint", &fRecpoints, bufsize); 

  Int_t nmodules = (Int_t) treeD->GetEntries();
  
  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      treeD->GetEntry(imodule); 
      Int_t nentries = fDigits->GetLast();
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmddigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
	  
	  det    = pmddigit->GetDetector();
	  smn    = pmddigit->GetSMNumber();
	  cellno = pmddigit->GetCellNumber();
	  adc    = pmddigit->GetADC();

	  ypos = cellno/fNCell;
	  xpos = cellno - ypos*fNCell;

	  if (det == 0)
	    {
	      fPMD[smn][xpos][ypos] = adc;
	    }
	  else if(det == 1)
	    {
	      fCPV[smn][xpos][ypos] = adc;
	    }
	}
    } // modules

  //
  // Clustering started
  //

  TObjArray *pmdcont = new TObjArray();

  AliPMDcluster  *pmdcl  = new AliPMDcluster;

  AliPMDClustering *pmdclust = new AliPMDClustering();
  pmdclust->SetMessage(fMessage);
  
  for (idet = 0; idet < 2; idet++)
    {
      for (isup = 0; isup < fTotSM; isup++)
	{
	  for (ix = 0; ix < fNCell; ix++)
	    {
	      for (iy = 0; iy < fNCell; iy++)
		{
		  if (idet == 0)
		    {
		      d[ix][iy] = (Double_t) fPMD[isup][ix][iy];
		    }
		  else if (idet == 1)
		    {
		      d[ix][iy] = (Double_t) fCPV[isup][ix][iy];
		    }
		    
		}
	    }
	  pmdclust->DoClust(idet,isup,d,pmdcont);
	  
	  Int_t nentries = pmdcont->GetEntries();
	  cout << " nentries = " << nentries << endl;
	  for (Int_t ient = 0; ient < nentries; ient++)
	    {
	      clusdata[0] = (Float_t) idet;
	      clusdata[1] = (Float_t) isup;
	      
	      pmdcl = (AliPMDcluster*)pmdcont->UncheckedAt(ient);
	      
	      clusdata[2] = pmdcl->GetClusX();
	      clusdata[3] = pmdcl->GetClusY();
	      clusdata[4] = pmdcl->GetClusADC();
	      clusdata[5] = pmdcl->GetClusCells();
	      clusdata[6] = pmdcl->GetClusRadius();
	      
	      AddRecPoint(clusdata);
	    }
	  pmdcont->Clear();
	  
	  treeR->Fill();
	  ResetRecpoint();
	}  // SuperModule
    }  // Detector
  
  ResetCellADC();
  
  pmdloader->WriteRecPoints("OVERWRITE");

  //   delete the pointers
  delete pmdclust;
  delete pmdcont;
    
  //  cout << " ***** End::Digits2RecPoints *****" << endl;
}


void AliPMDClusterFinder::AddRecPoint(Float_t *clusdata)
{
  TClonesArray &lrecpoints = *fRecpoints;
  AliPMDrecpoint *newrecpoint;
  newrecpoint = new AliPMDrecpoint(clusdata);
  new(lrecpoints[fNpoint++]) AliPMDrecpoint(newrecpoint);
  delete newrecpoint;
}
void AliPMDClusterFinder::ResetCellADC()
{
  for (Int_t i = 0; i < fTotSM; i++)
    {
      for (Int_t j = 0; j < fNCell; j++)
	{
	  for (Int_t k = 0; k < fNCell; k++)
	    {
	      fCPV[i][j][k] = 0.; 
	      fPMD[i][j][k] = 0.; 
	    }
	}
    }
}

void AliPMDClusterFinder::ResetRecpoint()
{
  fNpoint = 0;
  if (fRecpoints) fRecpoints->Clear();
}
void AliPMDClusterFinder::UnLoad(Option_t *option)
{
  const char *cR = strstr(option,"R");

  fRunLoader->UnloadgAlice();
  fRunLoader->UnloadHeader();
  fRunLoader->UnloadKinematics();

  if (cR)
    {
      pmdloader->UnloadDigits();
    }
}
