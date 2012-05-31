#if 0
  #include "$(ALICE_ROOT)/TPC/alles.h"
  #include "AliReader.h"
  #include "AliReaderKineTree.h"
  #include "AliReaderESDTree.h"
  #include "AliAODParticleCut.h"
  #include "AliAOD.h"
  #include "AliAODPairCut.h"
  #include "TSystem.h"
  #include "TObjString.h"
  #include "TString.h"
  #include "AliPDG.h"
#endif


void WriteAOD(Option_t* datatype, Int_t first = -1,Int_t last = -1,
                Option_t* processopt="TracksAndParticles",
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

  Bool_t multcheck = kTRUE;

  const char* basedir=".";
  const char* serie="";
  const char* field = "";
  cout<<"WriteAOD.C: datatype is "<<datatype<<" dir is basedir"<<endl;
  // Dynamically link some shared libs                    
  
  cout<<"WriteAOD.C: Loading  ANALYSIS .....\n";
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSIS");
  cout<<"WriteAOD.C: ..... Loaded\n";
  
  Int_t PID[11];

  PID[0]=kProton;
  PID[1]=kProtonBar;
  PID[2]=kKPlus;
  PID[3]=kKMinus;
  PID[4]=kPiPlus;
  PID[5]=kPiMinus;
  PID[6]=kElectron;
  PID[7]=kPositron;
  PID[8]=kMuonMinus;
  PID[9]=kMuonPlus;
  PID[10]=0;//Last must be 0!!!!!!!!!!!!!!!!!!
  
  Float_t PIDprob[11];
  PIDprob[0] = 0.5;
  PIDprob[1] = 0.5;
  PIDprob[2] = 0.5;
  PIDprob[3] = 0.5;
  PIDprob[4] = 0.5;
  PIDprob[5] = 0.5;
  PIDprob[6] = 0.5;
  PIDprob[7] = 0.5;
  PIDprob[8] = 0.5;
  PIDprob[9] = 0.5;
  PIDprob[10] = 0.5;
  /***********************************************************/
   
  AliReader* reader;
  Int_t kine = strcmp(datatype,"Kine");
  Int_t ESD = strcmp(datatype,"ESD");
  Int_t ESDMuon = strcmp(datatype,"ESDMuon");
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
    esdreader->SetNumberOfTrackPoints(5,30.);//5 points every 30 cm
    esdreader->SetITSTrackPoints(kTRUE);
    esdreader->SetClusterMap(kTRUE);
    
    reader = esdreader;
    multcheck = kTRUE;
   }

  else if(!intern)
   {
    AliReaderAOD* aodreader = new AliReaderAOD("AOD.root");
    if (strstr(processopt,"Particles"))
      aodreader->ReadSimulatedData(kTRUE);
    else
      aodreader->ReadSimulatedData(kFALSE);
    
    reader = aodreader;
    
    multcheck = kTRUE;
   }

  else if (!ESDMuon)
   {
     // set reader for ESD
     AliReaderESDTree* muonreader = new AliReaderESDTree("AliESDs.root");
     // active muon ESD reader
     muonreader->SetReadMuon(kTRUE);
     // disable central barrel (default = kTRUE)
     muonreader->SetReadCentralBarrel(kFALSE);
     // disable simulated data (not implemented yet)
     muonreader->ReadSimulatedData(kFALSE);

     reader = muonreader;
     multcheck = kFALSE;
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

  Int_t l = 0;
  while (PID[l] != 0)
   {
     cout<<"WriteAOD.C: Adding PID  = "<<PID[l]<<" l = "<<l<<endl;
     readerpartcut->SetPID(PID[l]);
     AliAODParticleCut * pcut = (AliAODParticleCut*)readerpartcut->Clone();
     AliAODPIDCut* pidcut = new AliAODPIDCut(PID[l],PIDprob[l]);
     pcut->AddBasePartCut(pidcut);
     reader->AddParticleCut(pcut);//read this particle type with this cut
     delete pcut;
     l++;
   }
  
//   readerpartcut->SetPtRange(0.0,1.2);

   cout<<"WriteAOD.C:   P R O C S E S S I N G .....\n\n";
   AliReaderAOD::WriteAOD(reader,outfile,"AliAODParticle",multcheck);
   cout<<"\n\nWriteAOD.C:   F I N I S H E D\n";
   
   if (dirs) delete dirs;
   delete reader;
 }

