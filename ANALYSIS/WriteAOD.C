#if 0
  #include "$(ALICE_ROOT)/TPC/alles.h"
  #include "AliReader.h"
  #include "AliReaderKineTree.h"
  #include "AliAODParticleCut.h"
  #include "AliAOD.h"
  #include "AliAODPairCut.h"
  #include "TSystem.h"
  #include "TObjString.h"
  #include "TString.h"
  #include "AliPDG.h"
#endif


void WriteAOD(Option_t* datatype, Option_t* processopt="TracksAndParticles",
                Int_t first = -1,Int_t last = -1,
                char *outfile = "AOD.root")
 {
//datatype defines type of data to be read
//  Kine  - analyzes Kine Tree: simulated particles
//  ESD
//  AOD

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

  const char* basedir=".";
  const char* serie="";
  const char* field = "";
  cout<<"WriteAOD.C: datatype is "<<datatype<<" dir is basedir"<<endl;
  // Dynamically link some shared libs                    
  
  cout<<"WriteAOD.C: Loading  ANALYSIS .....\n";
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSIS");
  cout<<"WriteAOD.C: ..... Loaded\n";
  
  Bool_t multcheck = kTRUE;
  /***********************************************************/
   
  AliReader* reader;
  Int_t kine = strcmp(datatype,"Kine");
  Int_t ESD = strcmp(datatype,"ESD");
  Int_t intern = strcmp(datatype,"AOD");

  if(!kine)
   {
    reader = new AliReaderKineTree();
    processopt="Particles"; //this reader by definition reads only simulated particles
    multcheck = kFALSE;
   }
  else if(!ESD)
   {
    AliReaderESD* esdreader = new AliReaderESD();
    esdreader->ReadSimulatedData(kTRUE);
    reader = esdreader;
    multcheck = kTRUE;
   }

  else if(!intern)
   {
    reader = new AliHBTReaderAOD("AOD.root");
    multcheck = kTRUE;
   }
  else
   {
    cerr<<"Option "<<datatype<<"  not recognized. Exiting"<<endl;
    return;
   }

  TObjArray* dirs=0;
  if ( (first >= 0) && (last>=0) && ( (last-first)>=0 ) )
   {//read from many dirs dirs
     cout<<"WriteAOD.C: ..... Setting dirs first="<<first<<" last="<<last<<"\n";
     char buff[50];
     dirs = new TObjArray(last-first+1);
     dirs->SetOwner();
     for (Int_t i = first; i<=last; i++)
      { 
        sprintf(buff,"%s/%s/%s/%d",basedir,field,serie,i);
        TObjString *odir= new TObjString(buff);
        dirs->Add(odir);
      }
    }

   reader->SetDirs(dirs);
    
   AliAODParticleCut* readerpartcut= new AliAODParticleCut();
   readerpartcut->SetPtRange(0.4,1.2);
   readerpartcut->SetPID(kPiPlus);
   AliAODPIDCut* pidcut = new AliAODPIDCut(kPiPlus,0.5);
   readerpartcut->AddBasePartCut(pidcut);
   
   reader->AddParticleCut(readerpartcut);//read this particle type with this cut

   cout<<"WriteAOD.C:   P R O C S E S S I N G .....\n\n";
   AliReaderAOD::WriteAOD(reader,outfile,"AliAODParticle",multcheck);
   cout<<"\n\nWriteAOD.C:   F I N I S H E D\n";
   
   if (dirs) delete dirs;
   delete reader;
 }

