#include "AliHBTReaderTPC.h"

#include <iostream.h>
//#include <fstream.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliMagF.h>
#include <AliTPCtrack.h>
#include <AliTPCParam.h>
#include <AliTPCtracker.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"


ClassImp(AliHBTReaderTPC)
//reader for TPC tracking
//needs galice.root, AliTPCtracks.root, AliTPCclusters.root, good_tracks_tpc 
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//Piotr.Skowronski@cern.ch

AliHBTReaderTPC::
 AliHBTReaderTPC(const Char_t* trackfilename,const Char_t* clusterfilename,
                 const Char_t* galicefilename):
                 fTrackFileName(trackfilename),fClusterFileName(clusterfilename),
                 fGAliceFileName(galicefilename)
{
  //constructor, 
  //Defaults:
  //  trackfilename = "AliTPCtracks.root"
  //  clusterfilename = "AliTPCclusters.root"
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untached
  
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();
  fIsRead = kFALSE;
}
/********************************************************************/
AliHBTReaderTPC::
AliHBTReaderTPC(TObjArray* dirs,
                  const Char_t* trackfilename, const Char_t* clusterfilename,
                  const Char_t* galicefilename):
                  AliHBTReader(dirs), fTrackFileName(trackfilename),
                  fClusterFileName(clusterfilename),fGAliceFileName(galicefilename)

{
  //constructor, 
  //Defaults:
  //  trackfilename = "AliTPCtracks.root"
  //  clusterfilename = "AliTPCclusters.root"
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
    if(Read(fParticles,fTracks))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }
   return fParticles->GetEvent(n);
 }
/********************************************************************/
AliHBTEvent* AliHBTReaderTPC::GetTrackEvent(Int_t n)
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

Int_t AliHBTReaderTPC::GetNumberOfPartEvents()
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
Int_t AliHBTReaderTPC::GetNumberOfTrackEvents()
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


Int_t AliHBTReaderTPC::Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  cout<<"AliHBTReaderTPC::Read()"<<endl;
  Int_t i; //iterator and some temprary values
  Int_t Nevents = 0;
  Int_t totalNevents = 0;
  TFile *aTracksFile;//file with tracks
  TFile *aClustersFile;//file with clusters
  TFile *aGAliceFile;//!ile name with galice

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
    
    if( (i=OpenFiles(aTracksFile,aClustersFile,aGAliceFile,currentdir)) )
     {
       Error("Read","Exiting due to problems with opening files. Errorcode %d",i);
       currentdir++;
       continue;
     }
  
    
    if (gAlice->TreeE())//check if tree E exists
     {
      Nevents = (Int_t)gAlice->TreeE()->GetEntries();//if yes get number of events in gAlice
      cout<<"________________________________________________________\n";
      cout<<"Found "<<Nevents<<" event(s) in directory "<<GetDirName(currentdir)<<endl;
      cout<<"Setting Magnetic Field: B="<<gAlice->Field()->SolenoidField()<<"T"<<endl;
      AliKalmanTrack::SetConvConst(1000/0.299792458/gAlice->Field()->SolenoidField());
     }
    else
     {//if not return an error
       Error("Read","Can not find Header tree (TreeE) in gAlice");
       currentdir++;
       continue;
     }
  
    aClustersFile->cd();//set cluster file active 
    AliTPCParam *TPCParam= (AliTPCParam*)aClustersFile->Get("75x40_100x60");
    if (!TPCParam) 
      { 
       Error("Read","TPC parameters have not been found !\n");
       currentdir++;
       continue;
      }

  
    for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
     {
       cout<<"Reading Event "<<currentEvent<<endl;
       /**************************************/
        /**************************************/
         /**************************************/ 
         
         aTracksFile->cd();//set track file active
          
         Char_t  treename[100];
         sprintf(treename,"TreeT_TPC_%d",currentEvent);//prepare name of the tree
   
         TTree *tracktree=0;
         
         tracktree=(TTree*)aTracksFile->Get(treename);//get the tree 
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
         cout<<"Found "<<NTPCtracks<<" TPC tracks.\n";
         //Copy tracks to array
         
         AliTPCtrack *iotrack=0;
         
         aClustersFile->cd();//set cluster file active 
         AliTPCtracker *tracker = new AliTPCtracker(TPCParam,currentEvent);//create the tacker for this event
         if (!tracker) //check if it has created succeffuly
          {//if not return with error
            Error("Read","Can't get a tracker !\n"); 
            continue;
          }
         tracker->LoadInnerSectors();
         tracker->LoadOuterSectors();
   
         for (i=0; i<NTPCtracks; i++) //loop over all tpc tracks
          {
            iotrack=new AliTPCtrack;   //create new tracks
            trackbranch->SetAddress(&iotrack); //tell the branch ehere to put track data from tree(file)
            tracktree->GetEvent(i); //stream track i to the iotrack
            tracker->CookLabel(iotrack,0.1); //calculate (cook) the label of the tpc track
                                             //which is the label of corresponding simulated particle 
            tarray->AddLast(iotrack); //put the track in the array
          }
         
         aTracksFile->Delete(treename);//delete tree from memmory (and leave untached on disk)- we do not need it any more
         aTracksFile->Delete("tracks");//delete branch from memmory
         delete tracker; //delete tracker
         
         tracker = 0x0;
         trackbranch = 0x0;
         tracktree = 0x0;
   
         Double_t xk;
         Double_t par[5];
         Float_t phi, lam, pt;//angles and transverse momentum
         Int_t label; //label of the current track
         
         aGAliceFile->cd();
         gAlice->GetEvent(currentEvent); 

         gAlice->Particles();
         
         for (i=0; i<NTPCtracks; i++) //loop over all good tracks
          { 
            iotrack = (AliTPCtrack*)tarray->At(i);
            label = iotrack->GetLabel();

            if (label < 0) continue;
            
            TParticle *p = (TParticle*)gAlice->Particle(label);
            
            if(p == 0x0) continue; //if returned pointer is NULL
            if(p->GetPDG() == 0x0) continue; //if particle has crezy PDG code (not known to our database)
	    
            if(Pass(p->GetPdgCode())) continue; //check if we are intersted with particles of this type 
                                        //if not take next partilce
            
            AliHBTParticle* part = new AliHBTParticle(*p);
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
            
            AliHBTParticle* track = new AliHBTParticle(p->GetPdgCode(), tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
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
    CloseFiles(aTracksFile,aClustersFile,aGAliceFile); 
    currentdir++;
   }while(currentdir < Ndirs);


  delete tarray;
  fIsRead = kTRUE;
  return 0;
 }

/********************************************************************/
Int_t AliHBTReaderTPC::OpenFiles
(TFile*& aTracksFile, TFile*& aClustersFile, TFile*& agAliceFile,Int_t event)
{
 //opens all the files
   
   
   const TString& dirname = GetDirName(event); 
   if (dirname == "")
    {
      Error("OpenFiles","Can not get directory name");
      return 4;
    }
   
   TString filename = dirname +"/"+ fTrackFileName;
   aTracksFile = TFile::Open(filename.Data());
   if ( aTracksFile  == 0x0 ) 
     {
       Error("OpenFiles","Can't open file with tacks named %s",filename.Data());
       return 1;
     }
   if (!aTracksFile->IsOpen())
     {
       Error("OpenFiles","Can't open file with tacks named %s",filename.Data());
       return 1;
     }
  
   filename = dirname +"/"+ fClusterFileName;
   aClustersFile = TFile::Open(filename.Data());
   if ( aClustersFile == 0x0 )
    {
      Error("OpenFiles","Can't open file with TPC clusters named %s",filename.Data());
      return 2;
    }
   if (!aClustersFile->IsOpen())
    {
      Error("OpenFiles","Can't open file with TPC clusters named %s",filename.Data());
      return 2;
    }

   filename = dirname +"/"+ fGAliceFileName;
   agAliceFile = TFile::Open(filename.Data());
   if ( agAliceFile== 0x0)
    {
      Error("OpenFiles","Can't open file with TPC clusters named %s",filename.Data());
      return 3;
    }
   if (!agAliceFile->IsOpen())
    {
      Error("OpenFiles","Can't open file with TPC clusters named %s",filename.Data());
      return 3;
    } 
   
   if (!(gAlice=(AliRun*)agAliceFile->Get("gAlice"))) 
    {
      Error("OpenFiles","gAlice have not been found on %s !\n",filename.Data());
      return 5;
    }

   return 0; 
}
/********************************************************************/

/********************************************************************/
  
void AliHBTReaderTPC::CloseFiles(TFile*& tracksFile, TFile*& clustersFile, TFile*& gAliceFile)
{
  //closes the files
  tracksFile->Close();
  delete tracksFile;
  tracksFile = 0x0;
  clustersFile->Close();
  delete clustersFile;
  clustersFile = 0x0;

  delete gAlice;
  gAlice = 0;

  if (gAliceFile) 
   {
     gAliceFile->Close();
     delete gAliceFile;
     gAliceFile = 0x0;
   }
}

/********************************************************************/

/********************************************************************/
/********************************************************************/
/********************************************************************/

