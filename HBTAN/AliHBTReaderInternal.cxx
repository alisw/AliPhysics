#include "AliHBTReaderInternal.h"

#include <iostream.h>
//#include <fstream.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

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
  fParticles = 0x0; 
  fTracks = 0x0;
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderInternal::AliHBTReaderInternal(const char *filename):fFileName(filename)
{ 
 fParticles = new AliHBTRun();
 fTracks    = new AliHBTRun();
 fIsRead = kFALSE;
}
/********************************************************************/
AliHBTReaderInternal::AliHBTReaderInternal(TObjArray* dirs, const char *filename):
  AliHBTReader(dirs),fFileName(filename)
{ 
 fParticles = new AliHBTRun();
 fTracks    = new AliHBTRun();
 fIsRead = kFALSE;
}
AliHBTReaderInternal::~AliHBTReaderInternal()
 {
 //desctructor
   delete fParticles;
   delete fTracks;
 }
/********************************************************************/
AliHBTEvent* AliHBTReaderInternal::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
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
 //returns Nth event with reconstructed tracks
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
  cout<<"AliHBTReaderInternal::Read()"<<endl;
  Int_t i; //iterator and some temprary values
  Int_t totalNevents = 0; //total number of read events
  Int_t Nevents = 0;
  Int_t currentdir = 0;
  Int_t Ndirs;
  
  TFile *aFile;//file with tracks
  AliHBTParticle* tpart = 0x0, *ttrack = 0x0;
  
  if (!particles) //check if an object is instatiated
   {
     Error("Read"," particles object must instatiated before passing it to the reader");
   }
  if (!tracks)  //check if an object is instatiated
   {
     Error("Read"," tracks object must instatiated before passing it to the reader");
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

  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
   {
    
    if( (i=OpenFile(aFile, currentdir)) )
     {
       Error("Read","Exiting due to problems with opening files. Errorcode %d",i);
       return i;
     }
   /***************************/
   /***************************/
   /***************************/
    
     TTree* tree = (TTree*)aFile->Get("data");
     if (tree == 0x0)
      {
       Error("Read","Can not get the tree");
       return 1;
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
     
     Nevents = (Int_t)tree->GetEntries();
     cout<<"________________________________________________________\n";
     cout<<"Found "<<Nevents<<" event(s) in directory "<<GetDirName(currentdir)<<endl;
     
     for (Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)
      {
       cout<<"Event "<<currentEvent<<"\n";
       tree->GetEvent(currentEvent);
         
       if (partbranch && trackbranch)
        {
           for(i = 0; i < pbuffer->GetEntries(); i++)
             {
               cout<<i<<"\r";
               tpart = dynamic_cast<AliHBTParticle*>(pbuffer->At(i));
               ttrack =  dynamic_cast<AliHBTParticle*>(tbuffer->At(i));

               if( (tpart == 0x0) || (ttrack == 0x0) ) continue; //if returned pointer is NULL
               if( (tpart->GetPDG()==0x0)||(ttrack->GetPDG()==0x0) ) continue; //if particle has crezy PDG code (not known to our database)
               if((Pass(tpart))||(Pass(ttrack))) continue; //check if we are intersted with particles of this type 
                                                          //if not take next partilce
               AliHBTParticle* part = new AliHBTParticle(*tpart);
               AliHBTParticle* track = new AliHBTParticle(*ttrack);
               particles->AddParticle(totalNevents,part);//put track and particle on the run
               tracks->AddParticle(totalNevents,track);
             }
            cout<<"Read: "<<particles->GetNumberOfParticlesInEvent(totalNevents)<<" particles and tracks \n";
        }
       else
        { 
         if (partbranch)
          {
            for(i = 0; i < pbuffer->GetEntries(); i++)
             {
               cout<<i<<"\r";
               tpart = dynamic_cast<AliHBTParticle*>(pbuffer->At(i));
               if(tpart == 0x0) continue; //if returned pointer is NULL
               if(tpart->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)
               if(Pass(tpart)) continue; //check if we are intersted with particles of this type 
                                         //if not take next partilce
 
               AliHBTParticle* part = new AliHBTParticle(*tpart);
               particles->AddParticle(totalNevents,part);//put track and particle on the run
             }
            cout<<"\nRead: "<<particles->GetNumberOfParticlesInEvent(totalNevents)<<" particles  \n";
          }
         else cout<<"Read: 0 particles  \n";
         
         if (trackbranch)
          {
            for(i = 0; i < tbuffer->GetEntries(); i++)
             {
               cout<<i<<"\r";
               tpart = dynamic_cast<AliHBTParticle*>(tbuffer->At(i));
               if(tpart == 0x0) continue; //if returned pointer is NULL
               if(tpart->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)
               if(Pass(tpart)) continue; //check if we are intersted with particles of this type 
                                                   //if not take next partilce
               AliHBTParticle* part = new AliHBTParticle(*tpart);
               tracks->AddParticle(totalNevents,part);//put track and particle on the run
             }
            cout<<"\nRead: "<<tracks->GetNumberOfParticlesInEvent(totalNevents)<<" tracks"<<endl;
          }
         else cout<<" 0 tracks"<<endl;
        }
        totalNevents++;
         
       }

   /***************************/
   /***************************/
   /***************************/
   currentdir++;
   aFile->Close();
   aFile = 0x0;

   }while(currentdir < Ndirs);

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




Int_t AliHBTReaderInternal::Write(AliHBTReader* reader,const char* outfile)
 {
  //reads tracks from runs and writes them to file
    Int_t i,j;
  
  
  cout<<"________________________________________________________\n";
  cout<<"________________________________________________________\n";
  cout<<"________________________________________________________\n";

  TFile *histoOutput = TFile::Open(outfile,"recreate");
  
  if (!histoOutput->IsOpen())
   {
     cout<<"File is not opened"<<endl;
     return 1;
   }
    

  TTree *tracktree = new TTree("data","Tree with tracks");

  TClonesArray* pbuffer = new TClonesArray("AliHBTParticle",15000);
  TClonesArray* tbuffer = new TClonesArray("AliHBTParticle",15000);

  TClonesArray &particles = *pbuffer;
  TClonesArray &tracks = *tbuffer;
    
  TString name("Tracks");
 
  
  Int_t NT = reader->GetNumberOfTrackEvents();
  Int_t NP = reader->GetNumberOfPartEvents();
  
  Bool_t trck = (NT > 0) ? kTRUE : kFALSE;
  Bool_t part = (NP > 0) ? kTRUE : kFALSE;

  TBranch *trackbranch = 0x0, *partbranch = 0x0;
  
  
  if (trck) trackbranch = tracktree->Branch("tracks","TClonesArray",&tbuffer);
  if (part) partbranch = tracktree->Branch("particles","TClonesArray",&pbuffer);


  
  if ( (trck) && (part) && (NP != NT))
   {
     cerr<<"Warning number of track and particle events is different"<<endl;
   }
  
  Int_t N;
  if (NT >= NP ) N = NT; else N = NP;

  for ( i =0;i< N; i++)
    {
      cout<<"Event "<<i+1<<"\n";
      if (trck && (i<=NT))
       {
         cout<<"Tracks: \n";
         AliHBTEvent* trackev = reader->GetTrackEvent(i);
         for ( j = 0; j< trackev->GetNumberOfParticles();j++)
          {
            cout<<j<<"\r";
            const AliHBTParticle& t= *(trackev->GetParticle(j));
            new (tracks[j]) AliHBTParticle(t);
          }

       }else cout<<"NO TRACKS";
      
      cout<<endl;
      
      if (part && (i<=NP))
       {
        cout<<"Particles: \n";
        AliHBTEvent* partev = reader->GetParticleEvent(i);
        for ( j = 0; j< partev->GetNumberOfParticles();j++)
         {
           const AliHBTParticle& part= *(partev->GetParticle(j));
           cout<<j<<"\r";
           new (particles[j]) AliHBTParticle(part);
         }
       }else cout<<"NO PARTICLES";
      cout<<endl;

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
