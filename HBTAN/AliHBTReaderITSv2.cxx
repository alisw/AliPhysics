
#include "AliHBTReaderITSv2.h"

#include <Riostream.h>
#include <Riostream.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliLoader.h>
#include <AliMagF.h>
#include <AliITS.h>
#include <AliITStrackV2.h>
#include <AliITStrackerV2.h>
#include <AliITSgeom.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"


ClassImp(AliHBTReaderITSv2)

AliHBTReaderITSv2::AliHBTReaderITSv2():fFileName("galice.root")
{
  //constructor, 
  //Defaults:
  //  galicefilename = "galice.root"
  fParticles = 0x0;
  fTracks    = 0x0;
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderITSv2::AliHBTReaderITSv2(const Char_t* galicefilename):fFileName(galicefilename)
{
  //constructor, 
  //Defaults:
  //  galicefilename = "galice.root"
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderITSv2::AliHBTReaderITSv2(TObjArray* dirs, const Char_t* galicefilename): 
       AliHBTReader(dirs), fFileName(galicefilename)
{
  //constructor, 
  //Defaults:
  //  galicefilename = "galice.root"
  
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderITSv2::~AliHBTReaderITSv2()
 {
   if (fParticles) delete fParticles;
   if (fTracks) delete fTracks;
 }
/********************************************************************/
/********************************************************************/

AliHBTEvent* AliHBTReaderITSv2::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
 if (!fIsRead) 
  { 
    if (fParticles == 0x0) fParticles = new AliHBTRun();
    if (fTracks == 0x0) fTracks    = new AliHBTRun();
    if(Read(fParticles,fTracks))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }
  }
 return (fParticles)?fParticles->GetEvent(n):0x0;
}
/********************************************************************/

AliHBTEvent* AliHBTReaderITSv2::GetTrackEvent(Int_t n)
 {
 //returns Nth event with reconstructed tracks
 if (!fIsRead) 
  { 
    if (fParticles == 0x0) fParticles = new AliHBTRun();
    if (fTracks == 0x0) fTracks    = new AliHBTRun();
    if(Read(fParticles,fTracks))
     {
       Error("GetTrackEvent","Error in reading");
       return 0x0;
     }
   }
  return (fTracks)?fTracks->GetEvent(n):0x0;
 }
/********************************************************************/

Int_t AliHBTReaderITSv2::GetNumberOfPartEvents()
 {
 //returns number of events of particles
  if (!fIsRead)
   {
    if (fParticles == 0x0) fParticles = new AliHBTRun();
    if (fTracks == 0x0) fTracks    = new AliHBTRun();
    if(Read(fParticles,fTracks))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
   }
  return (fParticles)?fParticles->GetNumberOfEvents():0;
 }

/********************************************************************/
Int_t AliHBTReaderITSv2::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (!fIsRead) 
   {
    if (fParticles == 0x0) fParticles = new AliHBTRun();
    if (fTracks == 0x0) fTracks    = new AliHBTRun();
    if(Read(fParticles,fTracks))
     {
       Error("GetNumberOfTrackEvents","Error in reading");
       return 0;
     }
   }
  return (fTracks)?fTracks->GetNumberOfEvents():0;
 }
/********************************************************************/
/********************************************************************/

 
Int_t AliHBTReaderITSv2::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
//reads data
 Int_t Nevents = 0; //number of events found in given directory
 Int_t Ndirs; //number of the directories to be read
 Int_t Ntracks; //number of tracks in current event
 Int_t currentdir = 0; //number of events in the current directory 
 Int_t totalNevents = 0; //total number of events read from all directories up to now
 register Int_t i = 0; //iterator
 
// AliITStrackerV2 *tracker; // ITS tracker - used for cooking labels
 TTree *tracktree; // tree for tracks
 
 Double_t xk;
 Double_t par[5]; //Kalman track parameters
 Float_t phi, lam, pt;//angles and transverse momentum
 Int_t label; //label of the current track

 AliITStrackV2 *iotrack = 0x0; //buffer track for reading data from tree

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
 
// cout<<"Found "<<Ndirs<<" directory entries"<<endl;
 
 do //do while is good even if Ndirs==0 (than read from current directory)
   {
    TString filename = GetDirName(currentdir);
    if (filename.IsNull())
     {
       Error("Read","Can not get directory name");
       currentdir++;
       continue;
     }
    filename = filename +"/"+ fFileName;
    AliRunLoader* rl = AliRunLoader::Open(filename);
    if( rl == 0x0)
     {
       Error("Read","Exiting due to problems with opening files.");
       currentdir++;
       continue;
     }
    
    rl->LoadHeader();
    rl->LoadKinematics();
    rl->LoadgAlice();
    gAlice = rl->GetAliRun();
    AliITS* its = (AliITS*)gAlice->GetModule("ITS");
    
    AliLoader* itsl = rl->GetLoader("ITSLoader");
    
    if ((its == 0x0) || ( itsl== 0x0))
     {
       Error("Read","Can not found ITS in this run");
       delete rl;
       rl = 0x0;
       currentdir++;
       continue;
     }
    Nevents = rl->GetNumberOfEvents();
 
    if (Nevents > 0)//check if tree E exists
     {
      Info("Read","________________________________________________________");
      Info("Read","Found %d event(s) in directory %s",Nevents,GetDirName(currentdir).Data());
      Float_t mf;
      if (fUseMagFFromRun)
       {
         mf = gAlice->Field()->SolenoidField();
         Info("Read","Setting Magnetic Field from run: B=%fT",mf/10.);
       }
      else
       {
         Info("Read","Setting Own Magnetic Field: B=%fT",fMagneticField);
         mf = fMagneticField*10.;
       }
      AliKalmanTrack::SetConvConst(1000/0.299792458/mf);
      if (iotrack == 0x0) iotrack = new AliITStrackV2();
     }
    else
     {//if not return an error
       Error("Read","No events in this run");
       delete rl;
       rl = 0x0;
       currentdir++;
       continue;
     }
    
    AliITSgeom *geom= its->GetITSgeom();
    if (!geom) 
     { 
       Error("Read","Can't get the ITS geometry!"); 
       delete rl;
       rl = 0x0;
       currentdir++;
       continue;
     }

    itsl->LoadTracks();

    for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
     {
       cout<<"Reading Event "<<currentEvent<<endl;
       rl->GetEvent(currentEvent);
       tracktree=itsl->TreeT();
       
       if (!tracktree) 
         {
           Error("Read","Can't get a tree with ITS tracks"); 
           continue;
         }
       TBranch *tbranch=tracktree->GetBranch("tracks");
       Ntracks=(Int_t)tracktree->GetEntries();

       Int_t accepted = 0;
       Int_t tpcfault = 0;
       Int_t itsfault = 0;
       for (i=0; i<Ntracks; i++) //loop over all tpc tracks
        { 
          if(i%100 == 0)cout<<"all: "<<i<<"   accepted: "<<accepted<<"   tpc faults: "<<tpcfault<<"\r";
          
          tbranch->SetAddress(&iotrack);
          tracktree->GetEvent(i);

          label=iotrack->GetLabel();
          if (label < 0) 
           {
             tpcfault++;
             continue;
           }

          TParticle *p = (TParticle*)gAlice->Particle(label);
          if(p == 0x0) continue; //if returned pointer is NULL
          if(p->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)

          if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                              //if not take next partilce
            
          AliHBTParticle* part = new AliHBTParticle(*p,i);
          if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                  //if it does not delete it and take next good track

          iotrack->PropagateTo(3.,0.0028,65.19);
          iotrack->PropagateToVertex();
 
          iotrack->GetExternalParameters(xk,par);     //get properties of the track
          phi=TMath::ASin(par[2]) + iotrack->GetAlpha(); 
          if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
          if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
          lam=par[3]; 
          pt=1.0/TMath::Abs(par[4]);
            
          Double_t tpx = pt * TMath::Cos(phi); //track x coordinate of momentum
          Double_t tpy = pt * TMath::Sin(phi); //track y coordinate of momentum
          Double_t tpz = pt * lam; //track z coordinate of momentum
           
          Double_t mass = p->GetMass();
          Double_t tEtot = TMath::Sqrt( tpx*tpx + tpy*tpy + tpz*tpz + mass*mass);//total energy of the track
            
          AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(), i, tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
          if(Pass(track))//check if meets all criteria of any of our cuts
                         //if it does not delete it and take next good track
           { 
            delete track;
            delete part;
            continue;
           }
          particles->AddParticle(totalNevents,part);//put track and particle on the run
          tracks->AddParticle(totalNevents,track);
          accepted++;
        }//end of loop over tracks in the event
       
       totalNevents++;
       cout<<"all: "<<i<<"   accepted: "<<accepted<<"   tpc faults: "<<tpcfault<<"   its faults: "<<itsfault<<endl;
     
     }//end of loop over events in current directory
    delete rl;
    currentdir++;
   }while(currentdir < Ndirs);//end of loop over directories specified in fDirs Obj Array

 delete iotrack;
 fIsRead = kTRUE;
 return 0;
}

/********************************************************************/
/********************************************************************/


