// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TParticle.h"
#include "TTree.h"
//
// // STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
//
// // MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
// //
// 


void rawddl4 (Int_t evNumber1 = 0, Int_t evNumber2 = 0)
{

  long int mappheader2, mapp2[123456];
  FILE *fb = fopen ("lutraw41.dat", "r");
  for (Int_t n = 0; n < 48448; n++)
    {
      fscanf (fb, "%ld", &mappheader2);
      fscanf (fb, "%ld", &mapp2[mappheader2]);
    }
  fclose (fb);
  long int mappheader4, mapp4[123456];
  FILE *fc = fopen ("lutraw42.dat", "r");
  for (Int_t n = 0; n < 48448; n++)
    {
      fscanf (fc, "%ld", &mappheader4);
      fscanf (fc, "%ld", &mapp4[mappheader4]);
    }
  fclose (fc);
  FILE *fp = fopen ("ddl13.dat", "w+");
  FILE *fq = fopen ("ddl14.dat", "w+");
  FILE *fr = fopen ("ddl15.dat", "w+");
  FILE *fs = fopen ("ddl16.dat", "w+");

  AliRunLoader *RunLoader =
    AliRunLoader::Open ("galice.root", "MUONFolder", "READ");
  if (RunLoader == 0x0)
    {
      printf (">>> Error : Error Opening %s file \n", "galice.root");
      return;
    }
  // Loading MUON subsystem
  AliLoader *MUONLoader = RunLoader->GetLoader ("MUONLoader");
  MUONLoader->LoadDigits ("READ");

  // Creating MUON data container
  AliMUONData muondata (MUONLoader, "MUON", "MUON");

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents ();
  printf (">>> No. of Event %d \n", nevents);
  AliMUONDigit *mDigit;
  
  Int_t dspheader[1000][4];
  Int_t dsphead[1000][100];

  //   Start loop over events 
  for (int ievent = evNumber1; ievent <= evNumber2; ievent++)	// event start 
    {
      printf ("Event:%d\n", ievent + 1);
      for (Int_t ii = 0; ii < 1000; ii++)
	{
	  for (Int_t ij = 0; ij < 100; ij++)
	    {
	      dsphead[ii][ij] = 0;
	    }
	}
      for (Int_t ik = 0; ik < 1000; ik++)
	{
	  for (Int_t il = 0; il < 4; il++)
	    {
	      dspheader[ik][il] = 0;
	    }
	}
//       printf(">>> Event %d \n",ievent);
      RunLoader->GetEvent (ievent);

      muondata.SetTreeAddress ("D");

      Int_t ncathodes = 2;
      for (Int_t icathode = 0; icathode < ncathodes; icathode++)
	{
	  muondata.GetCathode (icathode);

	  Int_t ichamber;
	  for (ichamber = 6; ichamber < 8; ichamber++)
	    {
	      Int_t idigit, ndigits;
	      ndigits = (Int_t) muondata.Digits (ichamber)->GetEntriesFast ();	//7 or 8

	      for (idigit = 0; idigit < ndigits; idigit++)
		{
		  mDigit =
		    static_cast <
		    AliMUONDigit * >(muondata.Digits (ichamber)->At (idigit));
		  {		// pads start  
		    Int_t dsp = 0, xydata = 0, manu = 0, adc = 0;
		    Int_t index = 0, counter = 0;
		    Int_t iqpad = mDigit->Signal ();	// charge per pad
		    if (iqpad >= 4096)
		      iqpad = 4095;
		    Int_t iPx = mDigit->PadX ();	// pad number on X
		    Int_t iPy = mDigit->PadY ();	// pad number on Y
//printf("root:%f\t%f\n",x,y);
		    Int_t modx = abs (iPx);
		    if (icathode == 0)
		      {
			xydata = 1024 * modx + iPy;
			manu = ((mapp2[xydata] >> 6) & 1023);
			adc = (((mapp2[xydata] << 11) >> 11) & 63);
			dsp = iPy / 80;
		      }
		    if (icathode == 1)
		      {
			xydata = 128 * modx + iPy;
			manu = ((mapp4[xydata] >> 6) & 1023);
			if (manu > 0)
			  {
			    adc = (((mapp4[xydata] << 11) >> 11) & 63);
			    if (manu <= 509)
			      dsp = 0;
			    if (manu > 509 && manu <= 525)
			      dsp = 1;
			    if (manu > 525 && manu <= 546)
			      dsp = 2;
			    if (manu > 546 && manu <= 579)
			      dsp = 3;
			    if (manu > 579 && manu <= 612)
			      dsp = 4;
			    if (manu > 612 && manu <= 645)
			      dsp = 5;
			    if (manu > 635 && manu <= 678)
			      dsp = 6;
			    if (manu > 668 && manu <= 711)
			      dsp = 7;
			    if (manu > 709 && manu <= 732)
			      dsp = 8;
			    if (manu > 722 && manu <= 748)
			      dsp = 9;
			    if (manu > 738 && manu <= 757)
			      dsp = 10;
			  }
		      }
		    if (iPx < 0)
		      {
			if (ichamber == 6)
			  {
			    dsp += 501;
			    index = 0;
			  }
			if (ichamber == 7)
			  {
			    dsp += 541;
			    index = 2;
			  }
			dspheader[dsp][index] += 1;
			counter = dspheader[dsp][index];
			dsphead[dsp][counter] =
			  4096 * (64 * (manu) + adc) + iqpad;
		      }
		    else
		      {
			if (ichamber == 6)
			  {
			    dsp += 521;
			    index = 1;
			  }
			if (ichamber == 7)
			  {
			    dsp += 561;
			    index = 3;
			  }
			dspheader[dsp][index] += 1;
			counter = dspheader[dsp][index];
			dsphead[dsp][counter] =
			  4096 * (64 * (manu) + adc) + iqpad;
		      }
		  }		//pad
		}		// digit
	    }			// chamber
	}			// cathode
      for (Int_t dsp = 501; dsp < 580; dsp++)
	{
	  for (Int_t aa = 0; aa < 4; aa++)
	    {
	      Int_t dspcount = dspheader[dsp][aa];
	      if (aa == 0)
		{
		  fwrite (&dsp, 2, 1, fp);
		  fwrite (&dspcount, 2, 1, fp);
		  for (Int_t ax = 1; ax <= dspcount; ax++)
		    fwrite (&dsphead[dsp][ax], 4, 1, fp);
		}
	      if (aa == 1)
		{
		  fwrite (&dsp, 2, 1, fq);
		  fwrite (&dspcount, 2, 1, fq);
		  for (Int_t ax = 1; ax <= dspcount; ax++)
		    fwrite (&dsphead[dsp][ax], 4, 1, fq);
		}
	      if (aa == 2)
		{
		  fwrite (&dsp, 2, 1, fr);
		  fwrite (&dspcount, 2, 1, fr);
		  for (Int_t ax = 1; ax <= dspcount; ax++)
		    fwrite (&dsphead[dsp][ax], 4, 1, fr);
		}
	      if (aa == 3)
		{
		  fwrite (&dsp, 2, 1, fs);
		  fwrite (&dspcount, 2, 1, fs);
		  for (Int_t ax = 1; ax <= dspcount; ax++)
		    fwrite (&dsphead[dsp][ax], 4, 1, fs);
		}
	    }
	}

    }				//event
  fclose (fp);
  fclose (fq);
  fclose (fr);
  fclose (fs);
}				// end of funtion     
