#if 0
  #include "/afs/cern.ch/user/s/skowron/aliroot/my/TPC/alles.h"
  #include "AliHBTReader.h"
  #include "AliHBTReaderKineTree.h"
  #include "AliHBTReaderITSv2.h"
  #include "AliHBTReaderITSv1.h"
  #include "AliHBTReaderTPC.h"
  #include "AliHBTParticleCut.h"
  #include "AliHBTEvent.h"
  #include "AliHBTPairCut.h"
  #include "AliHBTQResolutionFctns.h"
  #include "AliHBTTwoTrackEffFctn.h"
  #include "AliHBTCorrelFctn.h"
  #include "TSystem.h"
  #include "TObjString.h"
  #include "TString.h"
  #include "AliPDG.h"
#endif


void AliHBTWriteInternFormat(Option_t* datatype, Option_t* processopt="TracksAndParticles", 
                Int_t first = -1,Int_t last = -1,
                char *outfile = "data.root")
 {
//HBT Anlysis Macro
//Anlyzes TPC recontructed tracks and simulated particles that corresponds to them

//datatype defines type of data to be read
//  Kine  - analyzes Kine Tree: simulated particles
//  TPC   - analyzes TPC   tracking + particles corresponding to these tracks
//  ITSv1 - analyzes ITSv1 ----------------------//--------------------------
//  ITSv2 - analyzes ITSv2 ----------------------//--------------------------

//processopt defines option passed to AliHBTAnlysis::Process method
// default: TracksAndParticles - process both recontructed tracks and sim. particles corresponding to them 
//          Tracks - process only recontructed tracks
//          Particles - process only simulated particles

//Reads data from diroctories from first to last(including)
// For examples if first=3 and last=5 it reads from
//  ./3/
//  ./4/
//  ./5/
//if first or last is negative (or both), it reads from current directory
//
//these names I use when analysis is done directly from CASTOR files via RFIO
  const char* basedir="rfio:/castor/cern.ch/user/s/skowron";
  const char* serie="10k.12fm";
  const char* field = "0.2";
 
  // Dynamically link some shared libs                    
  
//  gROOT->LoadMacro("loadlibs.C");
//  loadlibs();
  
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libHBTAN");
  

  /***********************************************************/
  
  AliHBTReader* reader;
  Int_t kine = strcmp(datatype,"Kine");
  Int_t TPC = strcmp(datatype,"TPC");
  Int_t ITSv1 = strcmp(datatype,"ITSv1");
  Int_t ITSv2 = strcmp(datatype,"ITSv2");

  if(!kine)
   {
    reader = new AliHBTReaderKineTree();
    processopt="Particles"; //this reader by definition reads only simulated particles
   }
  else if(!TPC)
   {
    reader = new AliHBTReaderTPC();
   }
  else if(!ITSv1)
   {
    reader = new AliHBTReaderITSv1();
   }
  else if(!ITSv2)
   {
    reader = new AliHBTReaderITSv2();
   }
  else
   {
    cerr<<"Option "<<datatype<<"  not recognized. Exiting"<<endl;
    return;
   }

  TObjArray* dirs=0;
  if ( (first >= 0) && (last>=0) && ( (last-first)>=0 ) )
   {//read from many dirs dirs
     char buff[50];
     dirs = new TObjArray(last-first+1);
     for (Int_t i = first; i<=last; i++)
      { 
        sprintf(buff,"%s/%s/%s/%d",basedir,field,serie,i);
        TObjString *odir= new TObjString(buff);
        dirs->Add(odir);
      }
    }
   reader->SetDirs(dirs);
   AliHBTReaderInternal::Write(reader,outfile);
   
 }

