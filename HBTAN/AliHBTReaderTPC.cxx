#include "AliHBTReaderTPC.h"

#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliLoader.h>
#include <AliStack.h>
#include <AliMagF.h>
#include <AliTPCtrack.h>
#include <AliTPCParam.h>
#include <AliTPCtracker.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"


ClassImp(AliHBTReaderTPC)
//______________________________________________
//
// class AliHBTReaderTPC
//
//reader for TPC tracking
//needs galice.root, AliTPCtracks.root, AliTPCclusters.root
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//Piotr.Skowronski@cern.ch
AliHBTReaderTPC::AliHBTReaderTPC():fFileName("galice.root")
{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untouched
  
  fParticles = 0x0;
  fTracks    = 0x0;
  fIsRead = kFALSE;
}

AliHBTReaderTPC::AliHBTReaderTPC(const Char_t* galicefilename):fFileName(galicefilename)
{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untouched
  
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();
  fIsRead = kFALSE;
}
/********************************************************************/
AliHBTReaderTPC::AliHBTReaderTPC(TObjArray* dirs, const Char_t* galicefilename):
                  AliHBTReader(dirs), fFileName(galicefilename)

{
  //constructor, 
  //Defaults:
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untached
  
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderTPC::~AliHBTReaderTPC()
 {
 //desctructor
   delete fParticles;
   delete fTracks;
 }
/********************************************************************/

AliHBTEvent* AliHBTReaderTPC::GetParticleEvent(Int_t n)
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

AliHBTEvent* AliHBTReaderTPC::GetTrackEvent(Int_t n)
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

Int_t AliHBTReaderTPC::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead) 
    {
      if (fParticles == 0x0) fParticles = new AliHBTRun();
      if (fTracks == 0x0) fTracks    = new AliHBTRun();
      if ( Read(fParticles,fTracks))
       {
         Error("GetNumberOfPartEvents","Error in reading");
         return 0;
       }
    }
   return (fParticles)?fParticles->GetNumberOfEvents():0;
 }

/********************************************************************/
Int_t AliHBTReaderTPC::GetNumberOfTrackEvents()
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


Int_t AliHBTReaderTPC::Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("Read","");
  Int_t Nevents = 0;
  Int_t totalNevents = 0;

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
 
  TObjArray *tarray = new TObjArray(5000); //cotainer for tpc tracks
  tarray->SetOwner(); //set the ownership of the objects it contains
                      //when array is is deleted or cleared all objects 
                      //that it contains are deleted
  Int_t currentdir = 0;

  Int_t Ndirs;
  if (fDirs) //if array with directories is supplied by user
   {
     Ndirs = fDirs->GetEntries(); //get the number if directories
   }
  else
   {
     Ndirs = 0; //if the array is not supplied read only from current directory
   }

  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
   {
    TString filename = GetDirName(currentdir);
    if (filename.IsNull())
     {
       Error("Read","Can not get directory name");
       return 4;
     }
    filename = filename +"/"+ fFileName;
    AliRunLoader* rl = AliRunLoader::Open(filename,AliConfig::fgkDefaultEventFolderName);
    if( rl == 0x0)
     {
       Error("Read","Can not open session.");
       currentdir++;
       continue;
     }
    
    rl->LoadHeader();
    rl->LoadKinematics();
    AliLoader* tpcl = rl->GetLoader("TPCLoader");
    
    if ( tpcl== 0x0)
     {
       Error("Read","Exiting due to problems with opening files.");
       currentdir++;
       continue;
     }
    Nevents = rl->GetNumberOfEvents();
 
    if (Nevents > 0)//check if tree E exists
     {
      Info("Read","________________________________________________________");
      Info("Read","Found %d event(s) in directory %s",Nevents,GetDirName(currentdir).Data());
      rl->LoadgAlice();
      Info("Read","Setting Magnetic Field: B=%fT",rl->GetAliRun()->Field()->SolenoidField());
      AliKalmanTrack::SetConvConst(1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField());
      rl->UnloadgAlice();
     }
    else
     {//if not return an error
       Error("Read","Can not find Header tree (TreeE) in gAlice");
       currentdir++;
       continue;
     }
    
   rl->CdGAFile();
   AliTPCParam *TPCParam= (AliTPCParam*)gDirectory->Get("75x40_100x60");
   
   if (!TPCParam) 
    {
     TPCParam=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
     if (!TPCParam) 
      { 
        Error("Read","TPC parameters have not been found !\n");
        delete rl;
        rl = 0x0;
        currentdir++;
        continue;
      }
    }

    tpcl->LoadTracks();
  
    for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
     {
       Info("Read","Reading Event %d",currentEvent);
       /**************************************/
        /**************************************/
         /**************************************/ 
         rl->GetEvent(currentEvent);
         TTree *tracktree = tpcl->TreeT();//get the tree 
         if (!tracktree) //check if we got the tree
          {//if not return with error
            Error("Read","Can't get a tree with TPC tracks !\n"); 
            continue;
          }
   
         TBranch *trackbranch=tracktree->GetBranch("tracks");//get the branch with tracks
         if (!trackbranch) ////check if we got the branch
          {//if not return with error
            Error("Read","Can't get a branch with TPC tracks !\n"); 
            continue;
          }
         Int_t NTPCtracks=(Int_t)tracktree->GetEntries();//get number of TPC tracks 
         Info("Read","Found %d TPC tracks.",NTPCtracks);
         //Copy tracks to array
         
         AliTPCtrack *iotrack=0;
         
printf("This method is not converted to the NewIO !\n"); //I.B.
         //AliTPCtracker *tracker = new AliTPCtracker(TPCParam,currentEvent,AliConfig::fgkDefaultEventFolderName);//create the tacker for this event
         AliTPCtracker *tracker = new AliTPCtracker(TPCParam); //I.B.
         if (!tracker) //check if it has created succeffuly
          {//if not return with error
            Error("Read","Can't get a tracker !\n"); 
            continue;
          }
         tracker->LoadClusters(0);//I.Belikov, "0" must be a pointer to a tree
   
         for (Int_t i=0; i<NTPCtracks; i++) //loop over all tpc tracks
          {
            iotrack=new AliTPCtrack;   //create new tracks
            trackbranch->SetAddress(&iotrack); //tell the branch ehere to put track data from tree(file)
            tracktree->GetEvent(i); //stream track i to the iotrack
            tracker->CookLabel(iotrack,0.1); //calculate (cook) the label of the tpc track
                                             //which is the label of corresponding simulated particle 
            tarray->AddLast(iotrack); //put the track in the array
          }
         
         delete tracker; //delete tracker
         
         tracker = 0x0;
         trackbranch = 0x0;
         tracktree = 0x0;
   
         Double_t xk;
         Double_t par[5];
         Float_t phi, lam, pt;//angles and transverse momentum
         Int_t label; //label of the current track

         rl->Stack()->Particles();
         
         for (Int_t i=0; i<NTPCtracks; i++) //loop over all good tracks
          { 
            iotrack = (AliTPCtrack*)tarray->At(i);
            label = iotrack->GetLabel();

            if (label < 0) continue;
            
            TParticle *p = (TParticle*)rl->Stack()->Particle(label);
             
            if(p == 0x0) continue; //if returned pointer is NULL
            if(p->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)
           
            if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                        //if not take next partilce
            
            AliHBTParticle* part = new AliHBTParticle(*p,i);
            if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                    //if it does not delete it and take next good track
         
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
            
            AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(),i, tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
            if(Pass(track))//check if meets all criteria of any of our cuts
                         //if it does not delete it and take next good track
             { 
               delete track;
               delete part;
               continue;
             }
            particles->AddParticle(totalNevents,part);//put track and particle on the run
            tracks->AddParticle(totalNevents,track);

          }
         tarray->Clear(); //clear the array
         
        /**************************************/
       /**************************************/
      /**************************************/  
     totalNevents++;
    }
  
    //save environment (resouces) --
    //clean your place after the work
    delete rl;
    currentdir++;
   }while(currentdir < Ndirs);

  delete tarray;
  fIsRead = kTRUE;
  
  return 0;
 }


/********************************************************************/
/********************************************************************/
/********************************************************************/

