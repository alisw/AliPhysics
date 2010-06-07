// ----------------------------------------------------//
//                                                     //
//       This macro reads the PMD clusters which       //
//       are stored in the file "PMD.RecPoints.root"   //
//                                                     //
// ----------------------------------------------------//

#include <Riostream.h>
#include "TBranch.h"
#include "TStopwatch.h"

extern AliRun *gAlice;

Int_t AliPMDRecpointRead(Int_t nevent = 1)
{
  if (gAlice)
    { 
      delete AliRunLoader::Instance();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }
  AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root","Event","update");
  if (!fRunLoader)
    {
      cerr<<"Can't load RunLoader"<<endl;
      return 1;
    }
  AliLoader *pmdloader = fRunLoader->GetLoader("PMDLoader");
  Int_t nevent = fRunLoader->GetNumberOfEvents();
  cout << " * *********** nevent = " << nevent << endl;
  
  if (pmdloader == 0x0)
    {
      cerr<<" ===> Can not find PMD or PMDLoader <===\n";
      delete fRunLoader;
      return 2;
    }
  
  pmdloader->LoadRecPoints("READ");


  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  
  TClonesArray *fRecpoints;
  AliPMDUtility *cc = new AliPMDUtility();

  cc->ApplyAlignment();

  TH2F *h2 = new TH2F("h2"," ",100,-100.,100.,100,-100.,100.);

  FILE *fpw = fopen("junk_rec.dat","w");
  
  for (Int_t ievt = 0; ievt < nevent; ievt++)
    {
      fRunLoader->GetEvent(ievt);
      TTree *treeR = pmdloader->TreeR();
      if (treeR == 0x0)
	{
	  cout << " Can not get TreeR" << endl;
	  return 3;
	}

      AliPMDrecpoint1  *pmdrecpoint;
      TBranch *branch1 = treeR->GetBranch("PMDRecpoint");
      branch1->SetAddress(&fRecpoints);
      /**********************************************************************
       *    det   : Detector, 0: PRE & 1:CPV                                *
       *    smn   : Serial Module Number from 0 to 23 for both detector     *
       *    xpos  : x-position of the cluster                               *
       *    ypos  : y-position of the cluster                               *
       *            THESE xpos & ypos are not the true xpos and ypos        *
       *            for some of the unit modules. They are rotated.         *
       *    adc   : ADC contained in the cluster                            *
       *    ncell : Number of cells contained in the cluster                *
       *    rad   : radius of the cluster (1d fit)                          *
       *    xpad  : TRUE x-position of the cluster                          *
       *    ypad  : TRUE y-position of the cluster                          *
       **********************************************************************/

      Int_t   det,smn;
      Float_t xpos,ypos, xpad, ypad;
      Float_t adc, ncell, sigx, sigy;
      Float_t xx, yy;
      Int_t   nmodules = branch1->GetEntries();
      cout << " nmodules = " << nmodules << endl;
      for (Int_t imodule = 0; imodule < nmodules; imodule++)
	{
	  branch1->GetEntry(imodule); 
	  Int_t nentries = fRecpoints->GetLast();
	  for(Int_t ient = 0; ient < nentries+1; ient++)
	    {
	      pmdrecpoint = (AliPMDrecpoint1*)fRecpoints->UncheckedAt(ient);
	      det   = (Int_t) pmdrecpoint->GetDetector();
	      smn   = (Int_t) pmdrecpoint->GetSMNumber();
	      xpos  = pmdrecpoint->GetClusX();
	      ypos  = pmdrecpoint->GetClusY();
	      adc   = pmdrecpoint->GetClusADC();
	      ncell = pmdrecpoint->GetClusCells();
	      sigx  = pmdrecpoint->GetClusSigmaX();
	      sigy  = pmdrecpoint->GetClusSigmaY();

	      //
	      // User has to plug in his analysis code here
	      //

	      fprintf(fpw,"%d %d %d %d\n",
		      det,smn,xpos,ypos);
	      //
	      // Plot the cluster centroid to see the PMD geometry
	      // using the PMD Utility class
	      //
	      if (det == 1)
		{
		  // Draw only for PRE plane
		  cc->RectGeomCellPos(smn,xpos,ypos,xx,yy);
		  h2->Fill(xx,yy);
		}

	      //
	      // End of the User code
	      //
	    }
	}
      
    }

  h2->Draw();

  fclose(fpw);
  return 0;
}

