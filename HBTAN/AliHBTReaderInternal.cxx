#include "AliHBTReaderInternal.h"

#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TError.h>
#include <AliRun.h>
#include <AliMagF.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"


ClassImp(AliHBTReaderInternal)
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal()
{
//Defalut constructor
  fParticles = 0x0; 
  fTracks = 0x0;
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal(const char *filename):fFileName(filename)
{ 
//constructor 
//filename - name of file to open
 fParticles = new AliHBTRun();
 fTracks    = new AliHBTRun();
 fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal(TObjArray* dirs, const char *filename):
  AliHBTReader(dirs),fFileName(filename)
{ 
//ctor
//dirs contains strings with directories to look data in
//filename - name of file to open
 fParticles = new AliHBTRun();
 fTracks    = new AliHBTRun();
 fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderInternal::~AliHBTReaderInternal()
 {
 //desctructor
   delete fParticles;
   delete fTracks;
 }
/********************************************************************/

AliHBTEvent* AliHBTReaderInternal::GetParticleEvent(Int_t n)
 {
 //returns nth event with simulated particles
   if (!fIsRead) 
    if(Read(fParticles,fTracks))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }
   return fParticles->GetEvent(n);
 }
/********************************************************************/

AliHBTEvent* AliHBTReaderInternal::GetTrackEvent(Int_t n)
 {
 //returns nth event with reconstructed tracks
   if (!fIsRead) 
    if(Read(fParticles,fTracks))
     {
       Error("GetTrackEvent","Error in reading");
       return 0x0;
     }
   return fTracks->GetEvent(n);
 }
/********************************************************************/

Int_t AliHBTReaderInternal::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead) 
    if ( Read(fParticles,fTracks))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
   return fParticles->GetNumberOfEvents();
 }

/********************************************************************/
Int_t AliHBTReaderInternal::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (!fIsRead)
    if(Read(fParticles,fTracks))
     {
       Error("GetNumberOfTrackEvents","Error in reading");
       return 0;
     }
  return fTracks->GetNumberOfEvents();
 }
/********************************************************************/

Int_t AliHBTReaderInternal::Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("Read","");
  Int_t i; //iterator and some temprary values
  Int_t totalnevents = 0; //total number of read events
  Int_t nevents = 0;
  Int_t currentdir = 0;
  Int_t Ndirs;
  Int_t counter;
  TFile *aFile;//file with tracks
  AliHBTParticle* tpart = 0x0, *ttrack = 0x0;
  TDatabasePDG* pdgdb = TDatabasePDG::Instance();  
  if (pdgdb == 0x0)
   {
     Error("Read","Can not get PDG Particles Data Base");
     return 1;
   }
  if (!particles) //check if an object is instatiated
   {
     Error("Read"," particles object must instatiated before passing it to the reader");
     return 1;
   }
  if (!tracks)  //check if an object is instatiated
   {
     Error("Read"," tracks object must instatiated before passing it to the reader");
     return 1;
   }
  particles->Reset();//clear runs == delete all old events
  tracks->Reset();
 
  if (fDirs) //if array with directories is supplied by user
   {
     Ndirs = fDirs->GetEntries(); //get the number if directories
   }
  else
   {
     Ndirs = 0; //if the array is not supplied read only from current directory
   }

  TClonesArray* pbuffer = new TClonesArray("AliHBTParticle",15000);
  TClonesArray* tbuffer = new TClonesArray("AliHBTParticle",15000);
//  pbuffer->BypassStreamer(kFALSE);
//  tbuffer->BypassStreamer(kFALSE);
  
  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
   {
    
    if( (i=OpenFile(aFile, currentdir)) )
     {
       Error("Read","Skippimg directory due to problems with opening files. Errorcode %d",i);
       currentdir++;
       continue;
     }
   /***************************/
   /***************************/
   /***************************/
    
     TTree* tree = (TTree*)aFile->Get("data");
     if (tree == 0x0)
      {
       Error("Read","Can not get the tree");
       currentdir++;
       continue;
      }
     
     TBranch *trackbranch=tree->GetBranch("tracks");//get the branch with tracks
     if (trackbranch == 0x0) ////check if we got the branch
       {//if not return with error
         Warning("Read","Can't find a branch with tracks !\n"); 
       }
     else
      {
        trackbranch->SetAddress(&tbuffer);
      }

     TBranch *partbranch=tree->GetBranch("particles");//get the branch with particles
     if (partbranch == 0x0) ////check if we got the branch
       {//if not return with error
         Warning("Read","Can't find a branch with particles !\n"); 
       }
     else
      {
        partbranch->SetAddress(&pbuffer);
      }
     
     nevents = (Int_t)tree->GetEntries();
     Info("Read","________________________________________________________");
     Info("Read","Found %d event(s) in directory %s",nevents,GetDirName(currentdir).Data());
     
     for (Int_t currentEvent =0; currentEvent<nevents;currentEvent++)
      {
       Info("Read","Event %d",currentEvent);
       tree->GetEvent(currentEvent);
       
       counter = 0;  
       if (partbranch && trackbranch)
        {
           for(i = 0; i < pbuffer->GetEntries(); i++)
             {
               tpart = dynamic_cast<AliHBTParticle*>(pbuffer->At(i));
               ttrack =  dynamic_cast<AliHBTParticle*>(tbuffer->At(i));

               if( tpart == 0x0 ) continue; //if returned pointer is NULL
               
               if (ttrack->GetUID() != tpart->GetUID())
                 {
                   Error("Read","Sth. is wrong: Track and Particle has different UID.");
                   Error("Read","They probobly do not correspond to each other.");
                 }
               
               for (Int_t s = 0; s < tpart->GetNumberOfPids(); s++)
                {
                  if( pdgdb->GetParticle(tpart->GetNthPid(s)) == 0x0 ) continue; //if particle has crazy PDG code (not known to our database)
                  if( Pass(tpart->GetNthPid(s)) ) continue; //check if we are intersted with particles of this type
                                              //if not take next partilce
                  AliHBTParticle* part = new AliHBTParticle(*tpart);
                  part->SetPdgCode(tpart->GetNthPid(s),tpart->GetNthPidProb(s));
                  if( Pass(part) )
                    {
	  delete part;
	  continue; 
 	}
                  AliHBTParticle* track = new AliHBTParticle(*ttrack);
                  
                  particles->AddParticle(totalnevents,part);//put track and particle on the run
                  tracks->AddParticle(totalnevents,track);
                  counter++;
                }
             }
            Info("Read","   Read: %d particles and tracks.",counter);
        }
       else
        {  
         if (partbranch)
          {
            Info("Read","Found %d particles in total.",pbuffer->GetEntries());	
            for(i = 0; i < pbuffer->GetEntries(); i++)
             { 
               tpart = dynamic_cast<AliHBTParticle*>(pbuffer->At(i));
               if(tpart == 0x0) continue; //if returned pointer is NULL
               
               for (Int_t s = 0; s < tpart->GetNumberOfPids(); s++)
                {
                  if( pdgdb->GetParticle(tpart->GetNthPid(s)) == 0x0 ) continue; //if particle has crazy PDG code (not known to our database)
                  if( Pass(tpart->GetNthPid(s)) ) continue; //check if we are intersted with particles of this type
               
                  AliHBTParticle* part = new AliHBTParticle(*tpart);
                  part->SetPdgCode(tpart->GetNthPid(s),tpart->GetNthPidProb(s));
                  if( Pass(part) )
                    {
	  delete part;
	  continue; 
 	}
                  particles->AddParticle(totalnevents,part);//put track and particle on the run
                  counter++;
                }
             }
            Info("Read","   Read: %d particles.",counter);
          }
         else Info("Read","   Read: 0 particles.");
         
         if (trackbranch)
          {
            for(i = 0; i < tbuffer->GetEntries(); i++)
             {
               tpart = dynamic_cast<AliHBTParticle*>(tbuffer->At(i));
               if(tpart == 0x0) continue; //if returned pointer is NULL
               for (Int_t s = 0; s < tpart->GetNumberOfPids(); s++)
                {
                  if( pdgdb->GetParticle(tpart->GetNthPid(s)) == 0x0 ) continue; //if particle has crazy PDG code (not known to our database)
                  if( Pass(tpart->GetNthPid(s)) ) continue; //check if we are intersted with particles of this type
               
                  AliHBTParticle* part = new AliHBTParticle(*tpart);
                  part->SetPdgCode(tpart->GetNthPid(s),tpart->GetNthPidProb(s));
                  if( Pass(part) )
                    {
	  delete part;
	  continue; 
 	}
                  tracks->AddParticle(totalnevents,part);//put track and particle on the run
                  counter++;
                }
             }
            Info("Read","   Read: %d tracks",counter);
          }
         else Info("Read","   Read: 0 tracks.");
        }
        totalnevents++;
       }

    /***************************/
    /***************************/
    /***************************/
    currentdir++;
    delete tree;
    aFile->Close();
    delete aFile;
    aFile = 0x0;

   }while(currentdir < Ndirs);

  delete pbuffer;
  delete tbuffer;
  fIsRead = kTRUE;
  return 0;
}
/********************************************************************/

Int_t AliHBTReaderInternal::OpenFile(TFile*& aFile,Int_t event)
{

   const TString& dirname = GetDirName(event); 
   if (dirname == "")
    {
      Error("OpenFile","Can not get directory name");
      return 4;
    }
   
   TString filename = dirname +"/"+ fFileName;
   aFile = TFile::Open(filename.Data());
   if ( aFile  == 0x0 ) 
     {
       Error("OpenFiles","Can't open file with tacks named %s",filename.Data());
       return 1;
     }
   if (!aFile->IsOpen())
     {
       Error("OpenFiles","Can't open file with tacks named %s",filename.Data());
       return 1;
     }
   return 0; 
}
/********************************************************************/

Int_t AliHBTReaderInternal::Write(AliHBTReader* reader,const char* outfile, Bool_t multcheck)
 {
  //reads tracks from reader and writes runs to file
  //reader - provides data for writing in internal format
  //name of output file
  //multcheck - switches of checking if particle was stored with other incarnation
  // usefull e.g. when using kine data, where all particles have 100% pid prob.and saves a lot of time
  
  Int_t i,j;
  
  ::Info("AliHBTReaderInternal::Write","________________________________________________________");
  ::Info("AliHBTReaderInternal::Write","________________________________________________________");
  ::Info("AliHBTReaderInternal::Write","________________________________________________________");

  TFile *histoOutput = TFile::Open(outfile,"recreate");
  
  if (!histoOutput->IsOpen())
   {
     ::Error("AliHBTReaderInternal::Write","File is not opened");
     return 1;
   }
    
  TTree *tracktree = new TTree("data","Tree with tracks");

  TClonesArray* pbuffer = new TClonesArray("AliHBTParticle",15000);
  TClonesArray* tbuffer = new TClonesArray("AliHBTParticle",15000);
//  pbuffer->BypassStreamer(kFALSE);
//  tbuffer->BypassStreamer(kFALSE);

  TClonesArray &particles = *pbuffer;
  TClonesArray &tracks = *tbuffer;
    
  TString name("Tracks");
  
  Int_t nt = reader->GetNumberOfTrackEvents();
  Int_t np = reader->GetNumberOfPartEvents();
  
  if (AliHBTParticle::GetDebug() > 0)
   ::Info("Write","Reader has %d track events and %d particles events.",nt,np);
   
  Bool_t trck = (nt > 0) ? kTRUE : kFALSE;
  Bool_t part = (np > 0) ? kTRUE : kFALSE;

  TBranch *trackbranch = 0x0, *partbranch = 0x0;
  
  if (trck) trackbranch = tracktree->Branch("tracks",&tbuffer,32000,0);
  if (part) partbranch = tracktree->Branch("particles",&pbuffer,32000,0);
  
  if ( (trck) && (part) && (np != nt))
   {
     ::Warning("AliHBTReaderInternal::Write","Number of track and particle events is different");
   }
  
  Int_t n;
  if (nt >= np ) n = nt; else n = np;
  
  if (AliHBTParticle::GetDebug() > 0)
   ::Info("Write","Will loop over %d events",n);

  for ( i =0;i< n; i++)
    {
      ::Info("AliHBTReaderInternal::Write","Event %d",i+1);
      Int_t counter = 0;
      if (trck && (i<=nt))
       { 
         AliHBTEvent* trackev = reader->GetTrackEvent(i);
         for ( j = 0; j< trackev->GetNumberOfParticles();j++)
          {
            const AliHBTParticle& t = *(trackev->GetParticle(j));
            if (multcheck)
             {
              if (FindIndex(tbuffer,t.GetUID())) 
               {
                 if (AliHBTParticle::GetDebug()>4)
                  { 
                   ::Info("Write","Track with Event UID %d already stored",t.GetUID());
                  }
                 continue; //not to write the same particles with other incarnations
               }
             }
            new (tracks[counter++]) AliHBTParticle(t);
          }
         ::Info("AliHBTReaderInternal::Write","    Tracks: %d",tracks.GetEntries());
       }else ::Info("AliHBTReaderInternal::Write","NO TRACKS");
      
      counter = 0;
      if (part && (i<=np))
       {
//        ::Warning("AliHBTReaderInternal::Write","Find index switched off!!!");

        AliHBTEvent* partev = reader->GetParticleEvent(i);
        for ( j = 0; j< partev->GetNumberOfParticles();j++)
         {
           const AliHBTParticle& part= *(partev->GetParticle(j));
            if (multcheck)
             {
              if (FindIndex(pbuffer,part.GetUID())) 
               {
                 if (AliHBTParticle::GetDebug()>4)
                  { 
                   ::Info("Write","Particle with Event UID %d already stored",part.GetUID());
                  }
                 continue; //not to write the same particles with other incarnations
               }
             } 
           new (particles[counter++]) AliHBTParticle(part);
         }
         ::Info("AliHBTReaderInternal::Write","    Particles: %d",particles.GetEntries());
       }else ::Info("AliHBTReaderInternal::Write","NO PARTICLES");

      histoOutput->cd();
      tracktree->Fill();
      tracktree->AutoSave();
      tbuffer->Delete();
      pbuffer->Delete();
     }

  histoOutput->cd();
  tracktree->Write(0,TObject::kOverwrite);
  delete tracktree;

  tbuffer->SetOwner();
  pbuffer->SetOwner();
  delete pbuffer;
  delete tbuffer;

  histoOutput->Close();
  return 0;
 }
/********************************************************************/

Bool_t AliHBTReaderInternal::FindIndex(TClonesArray* arr,Int_t idx)
{
//Checks if in the array exists already partilce with Unique ID idx
  if (arr == 0x0)
   {
     ::Error("FindIndex","Array is 0x0");
     return kTRUE;
   }
  TIter next(arr);
  AliHBTParticle* p;
  while (( p = (AliHBTParticle*)next()))
   {
     if (p->GetUID() == idx) return kTRUE;
   }
  return kFALSE;
}
