// ----------------------------------------------------//
//                                                     //
//       This macro reads the PMD clusters,            //
//       does photon/hadron discrimination             //
//       and stores in the ESD                         //
//                                                     //
// ----------------------------------------------------//

#include <Riostream.h>
#include "TBranch.h"
#include "TStopwatch.h"

extern AliRun *gAlice;

Int_t AliESDPmdTest(Int_t nevent = 1)
{
  if (gAlice)
    { 
      delete AliRunLoader::Instance();
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
  pmdloader->LoadRecPoints("READ");

  TFile *bf=TFile::Open("AliESDcheck.root","RECREATE");
  if (!bf || !bf->IsOpen()) {
    cerr<<"Can't open AliESDcheck.root !\n"; return 1;
  }

  Char_t ename[100]; 

  for (Int_t ievt = 0; ievt < nevent; ievt++)
    {

      AliESD *event=new AliESD(); 
      Int_t run = 0;
      event->SetRunNumber(run);
      event->SetEventNumber(ievt);

      fRunLoader->GetEvent(ievt);
      TTree *treeR = pmdloader->TreeR();
      if (treeR == 0x0)
	{
	  cout << " Can not get TreeR" << endl;
	  return 3;
	}
      AliPMDtracker *pmdtracker = new AliPMDtracker();
      pmdtracker->LoadClusters(treeR);
      pmdtracker->Clusters2Tracks(event);

      bf->cd();
      sprintf(ename,"in%d",ievt);
      event->Write(ename); bf->Flush();
      
      delete pmdtracker;  
      delete event;
    }
  bf->Close();
  return 0;
}

