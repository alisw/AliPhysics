// ----------------------------------------------------//
//                                                     //
//       This macro does digits to Raw Data            //
//                                                     //
// ----------------------------------------------------//

#include <Riostream.h>
#include "TBranch.h"
#include "TStopwatch.h"

extern AliRun *gAlice;

Int_t AliPMDDigits2RawData(Int_t nevent = 1)
{
  if (gAlice)
    { 
      delete gAlice->GetRunLoader();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }

  AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");
  if (!fRunLoader)
    {
      cerr<<"Can't load RunLoader"<<endl;
      return 1;
    }
  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();

  AliLoader *pmdloader = fRunLoader->GetLoader("PMDLoader");
  //  Int_t nevent = fRunLoader->GetNumberOfEvents();
  cout << " ************ nevent = " << nevent << endl;
  
  if (pmdloader == 0x0)
    {
      cerr<<" ===> Can not find PMD or PMDLoader <===\n";
      delete fRunLoader;
      return 2;
    }
  pmdloader->LoadDigits("READ");


  for (Int_t ievt = 0; ievt < nevent; ievt++)
    {

      fRunLoader->GetEvent(ievt);
      TTree *treeD = pmdloader->TreeD();
      if (treeD == 0x0)
	{
	  cout << " Can not get TreeD" << endl;
	  return 3;
	}

      AliPMDDDLRawData rawdata;
      rawdata.WritePMDRawData(treeD, ievt);


    }

  return 0;
}

