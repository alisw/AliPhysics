#if !defined(__CINT__) || defined(__MAKECINT__)
#include<Riostream.h>
#include<TROOT.h>
#include<TClassTable.h>
#include<TClonesArray.h>
#include<TGeoManager.h>
#include <TInterpreter.h>
#include<TTree.h>
#include "AliRun.h"
#include "AliITSgeom.h"
#include "AliGeomManager.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSdigit.h"
#include "AliITSdigitSSD.h"
#include "AliITShit.h"
#include "AliITSmodule.h" 
#include "AliITSsegmentation.h"
#include "AliITSsegmentationSPD.h" 
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliHeader.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif
void AliITSPrintRecPoints(Int_t outtype=1,TString rfn="galice.root",
                          Int_t mod=-1,Int_t evnt=-1){
  // Macro to print out the recpoints for all or a specific module

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
  else {
    if(gAlice){
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice=0;
    }
  }

  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number 0\n");
    man->SetDefaultStorage("local://$ALICE_ROOT");
    man->SetRun(0);
  }
  else {
    printf("Using deafult storage \n");
  }

  AliRunLoader* rl = AliRunLoader::Open(rfn.Data());
  if (rl == 0x0){
    cerr<<"AliITSPrintRecPoints.C : Can not open session RL=NULL"<< endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"AliITSPrintRecPoints.C : LoadgAlice returned error"<<endl;
    return;
  }
  gAlice=rl->GetAliRun();

  retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSPrintRecPoints.C : LoadHeader returned error"<<endl;
    return;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");

  if(!ITSloader){
    cerr<<"AliITSPrintRecPoints.C :  ITS loader not found"<<endl;
    return;
  }
  cout <<"ITSloader ok"<<endl;

  if(!gGeoManager){
    AliGeomManager::LoadGeometry("geometry.root");
  } 
  ITSloader->LoadHits("read");
  ITSloader->LoadDigits("read");
  ITSloader->LoadRecPoints("read");
  cout << "loaded hits, digits, and RecPoints"<< endl;
 
  AliITSgeom *gm=0;
  gm = ITSloader->GetITSgeom();

  Int_t evNumber1 = 0;
  Int_t evNumber2 = AliRunLoader::GetNumberOfEvents();
  if(evnt>=0){
    evNumber1 = evnt;
    evNumber2 = evnt+1;
  } // end if evnt>=0
  Int_t mod1 = 0;
  Int_t mod2 = gm->GetIndexMax();
  if(mod>=0){
    mod1 = mod;
    mod2 = mod+1;
  } // end if mod>=0
  TClonesArray *rpa;
  AliITSRecPoint *rp = 0;
  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetITSgeom(gm);
  rec->SetDefaults();

  Int_t event,m,i,i2;
  Float_t xyz[3];
  for(event = evNumber1; event < evNumber2; event++){
    rl->GetEvent(event);
    //    rec->SetTreeAddress();
    for(m=mod1;m<mod2;m++){
      rec->ResetRecPoints();
      TTree *TR = ITSloader->TreeR();
      rec->SetTreeAddressR(TR);
      TR->GetEvent(m);
      rpa = rec->RecPoints();
      i2 = rpa->GetEntriesFast();
      cout <<  "Event=" << event << " module=" << m <<
	" Number of Recpoints=" << i2 <<endl;
      for(i=0;i<i2;i++){
          rp = (AliITSRecPoint*)(rpa->At(i));
          switch(outtype){
          case 1:
              rp->GetGlobalXYZ(xyz);
              cout << i << " lx=" << rp->GetDetLocalX() 
                        << " lz=" << rp->GetDetLocalZ() 
                   << " x=" << rp->GetX() 
                   << " y=" << rp->GetY()<< " z=" << rp->GetZ() 
                   <<" gx=" << xyz[0] << " gy="<< xyz[1] <<" gz="<<xyz[2]
                   << endl;
              break;
          default:
              cout << i << " ";
              rp->Print();
              cout << endl;
              break;
          } // end switch
      } // end for i
    } // end for m
  } // end for event

}
