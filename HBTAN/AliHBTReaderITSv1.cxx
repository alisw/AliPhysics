//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
//  class AliHBTReaderITSv1                                         //
//                                                                  //
//  Reader for ITSv1 tracks. Not maintained since v1 is not         //
//  supposed to be used                                             //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliHBTReaderITSv1.h"
#include "AliHBTEvent.h"
#include "AliHBTRun.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"

#include <Riostream.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TString.h>
#include <TObjString.h>

#include <AliRun.h>
#include <AliStack.h>
#include <AliMagF.h>
#include <AliKalmanTrack.h>
#include <AliITSIOTrack.h>

ClassImp(AliHBTReaderITSv1)
/********************************************************************/

AliHBTReaderITSv1::
AliHBTReaderITSv1(const Char_t* tracksfilename,const Char_t* galicefilename):
                 fITSTracksFileName(tracksfilename),fGAliceFileName(galicefilename)
 {
     fParticles = new AliHBTRun();
     fTracks    = new AliHBTRun();
     fIsRead = kFALSE;
 }
/********************************************************************/

AliHBTReaderITSv1::
AliHBTReaderITSv1(TObjArray* dirs, const Char_t* tracksfilename,const Char_t* galicefilename):
                 AliHBTReader(dirs),
                 fITSTracksFileName(tracksfilename),fGAliceFileName(galicefilename)
 {
   fParticles = new AliHBTRun();
   fTracks    = new AliHBTRun();
   fIsRead    = kFALSE;
 }
/********************************************************************/

AliHBTReaderITSv1::~AliHBTReaderITSv1()
{
   delete fParticles;
   delete fTracks;
}
/********************************************************************/

AliHBTEvent* AliHBTReaderITSv1::GetParticleEvent(Int_t n)
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

AliHBTEvent* AliHBTReaderITSv1::GetTrackEvent(Int_t n)
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

Int_t AliHBTReaderITSv1::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead)
    if(Read(fParticles,fTracks))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
   return fParticles->GetNumberOfEvents();
 }

/********************************************************************/
Int_t AliHBTReaderITSv1::GetNumberOfTrackEvents()
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

Int_t AliHBTReaderITSv1::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
 Int_t Nevents = 0;
 AliITSIOTrack *iotrack=new AliITSIOTrack;
 Int_t currentdir = 0;
 Int_t Ndirs;
 Int_t totalNevents = 0;
 
 if (fDirs)
  {
    Ndirs = fDirs->GetEntries();
  }
 else
  {
    Ndirs = 0;
  }
 
 do //do while is good even if 
  {  
   TFile* gAliceFile = OpenGAliceFile(currentdir);
   if(gAliceFile == 0x0)
    {
       Error("Read","Can not open the file with gAlice");
       delete iotrack;
       return 1;
    }
   if (gAlice->TreeE())//check if tree E exists
     {
      Nevents = (Int_t)gAlice->TreeE()->GetEntries();//if yes get number of events in gAlice
      cout<<"________________________________________________________\n";
      cout<<"Found "<<Nevents<<" event(s) in directory "<<GetDirName(currentdir)<<endl;
      cout<<"Setting Magnetic Field. Factor is "<<gAlice->Field()->Factor()<<endl;
      AliKalmanTrack::SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());
     }
    else
     {//if not return an error
       Error("Read","Can not find Header tree (TreeE) in gAlice");
       delete iotrack;
       return 4;
     }

   TFile *file = OpenTrackFile(currentdir);
   if(file == 0x0)
    {
       Error("Read","Can not open the file with ITS tracks V1");
       delete iotrack;
       return 2;
    }
    
   Int_t naccepted = 0;
   char tname[30];
   
   for (Int_t currentEvent = 0; currentEvent < Nevents; currentEvent++)
    { 
      cout<<"Reading Event "<<currentEvent;
      
      sprintf(tname,"TreeT%d",currentEvent);
      file->cd(); 
      TTree *tracktree=(TTree*)file->Get(tname);
      TBranch *tbranch=tracktree->GetBranch("ITStracks");
      tbranch->SetAddress(&iotrack);
      
      gAliceFile->cd();
      gAlice->GetEvent(currentEvent);
      gAlice->Stack()->Particles();

      Int_t nentr=(Int_t)tracktree->GetEntries();
      
      cout<<".  Found "<<nentr<<" tracks.";
      fflush(0);
      
      for (Int_t i=0; i<nentr; i++) 
       {

        tracktree->GetEvent(i);
        if(!iotrack) continue;       
        Int_t label = iotrack->GetLabel();
        if (label < 0) 
         {
           continue;
         }

        TParticle *p = (TParticle*)gAlice->Stack()->Particle(label);
        if(!p)
         {
           Warning("Read","Can not get particle with label &d",label);
           continue;
         }
        if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type
                                           //if not take next partilce

        AliHBTParticle* part = new AliHBTParticle(*p,i);
        if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                //if it does not delete it and take next good track
        
        Double_t px=iotrack->GetPx();
        Double_t py=iotrack->GetPy();
        Double_t pz=iotrack->GetPz();
        Double_t mass = p->GetMass();
        Double_t tEtot = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);//total energy of the track
 
        Double_t x= iotrack->GetX();
        Double_t y= iotrack->GetY();
        Double_t z= iotrack->GetZ();
        
        AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(), i, px, py , pz, tEtot, x, y, z, 0.);
        if(Pass(track)) { delete  track;continue;}//check if meets all criteria of any of our cuts
                                                  //if it does not delete it and take next good track

        particles->AddParticle(totalNevents,part);//put track and particle on the run
        tracks->AddParticle(totalNevents,track);
        naccepted++;
       }//end loop over tracks in the event

       totalNevents++;
       cout<<"  Accepted "<<naccepted<<" tracks"<<endl;
     }//end of loop over events in current directory
    
    gAliceFile->Close();
    delete gAliceFile;
    gAliceFile = 0;
    
    file->Close(); 
    delete file;
    file = 0;
    currentdir++;
   }while(currentdir < Ndirs);//end of loop over directories specified in fDirs Obj Array


  delete iotrack;
  fIsRead = kTRUE;
  return 0;
 
 }
/********************************************************************/

TFile* AliHBTReaderITSv1::OpenTrackFile(Int_t ndir)
{
//opens files to be read for given directoru nomber in fDirs Array
   const TString& dirname = GetDirName(ndir); 
   if (dirname == "")
    {
      Error("OpenGAliceFile","Can not get directory name");
      return 0x0;
    }
   TString filename = dirname + "/" + fITSTracksFileName;

   TFile *file = TFile::Open(filename.Data());   
   if (!file)
    {
      Error("Read","Can not open file %s",filename.Data());
      return 0x0;
    }
   if (!file->IsOpen())
    {
      Error("Read","Can not open file %s",filename.Data());
      return 0x0;
    }
   
   return file;
}


/********************************************************************/
TFile* AliHBTReaderITSv1::OpenGAliceFile(Int_t ndir)
{
  const TString& dirname = GetDirName(ndir); 
   if (dirname == "")
    {
      Error("OpenGAliceFile","Can not get directory name");
      return 0x0;
    }
  
  TString filename = dirname + "/" + fGAliceFileName;

  TFile* gAliceFile = TFile::Open(filename.Data());
  if ( gAliceFile== 0x0)
   {
     Error("OpenFiles","Can't open file named %s",filename.Data());
     return 0x0;
   }
  if (!gAliceFile->IsOpen())
   {
     Error("OpenFiles","Can't open file named %s",filename.Data());
     return 0x0;
   }

  if (!(gAlice=(AliRun*)gAliceFile->Get("gAlice")))
   {
     Error("OpenFiles","gAlice have not been found on %s !\n",filename.Data());
     gAliceFile->Close();
     delete gAliceFile;
     return 0x0;
   }
  
  return gAliceFile;
}

/********************************************************************/
/********************************************************************/


