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
      delete gAlice->GetRunLoader();
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
  TClonesArray *fRecpoints;
  AliPMDUtility *cc = new AliPMDUtility();
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
       *    smn   : Serial Module Number from which Super Module Number     *
       *            and Unit Module Numbers are extracted                   *
       *    xpos  : x-position of the cluster                               *
       *    ypos  : y-position of the cluster                               *
       *            THESE xpos & ypos are not the true xpos and ypos        *
       *            for some of the unit modules. They are rotated.         *
       *    adc   : ADC contained in the cluster                            *
       *    ncell : Number of cells contained in the cluster                *
       *    rad   : radius of the cluster (1d fit)                          *
       *    ism   : Supermodule number extracted from smn                   *
       *    ium   : Unit module number extracted from smn                   *
       *    xpad  : TRUE x-position of the cluster                          *
       *    ypad  : TRUE y-position of the cluster                          *
       **********************************************************************/
      Int_t   ism, ium;
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
	      // Now change the xpos and ypos to its original values
	      // for the unit modules which are earlier changed.
	      // xpad and ypad are the real positions.
	      //
	      if(det == 0 || det == 1)
		{
		  if(smn < 12)
		    {
		      ism  = smn/6;
		      ium  = smn - ism*6;
		      xpad = ypos;
		      ypad = xpos;
		    }
		  else if( smn >= 12 && smn < 24)
		    {
		      ism  = smn/6;
		      ium  = smn - ism*6;
		      xpad = xpos;
		      ypad = ypos;
		    }
		}
	      //
	      // User has to plug in his analysis code here
	      //

	      fprintf(fpw,"%d %d %d %d %f %f %f %f %f %f\n",
		      det,smn,ism,ium,xpad,ypad,adc,ncell,sigx,sigy);
	      //
	      // Plot the cluster centroid to see the PMD geometry
	      // using the PMD Utility class
	      //
	      if (det == 1)
		{
		  // Draw only for PRE plane
		  cc->RectGeomCellPos(ism,ium,xpad,ypad,xx,yy);
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

