
#include "AliHBTReaderITSv2.h"

#include <Riostream.h>
#include <Riostream.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include <AliRun.h>
#include <AliMagF.h>
#include <AliITStrackV2.h>
//#include <AliITSParam.h>
#include <AliITStrackerV2.h>
#include <AliITSgeom.h>

#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"


ClassImp(AliHBTReaderITSv2)

AliHBTReaderITSv2::
 AliHBTReaderITSv2(const Char_t* trackfilename, const Char_t* clusterfilename,
	const Char_t* galicefilename)
                  :fTrackFileName(trackfilename),fClusterFileName(clusterfilename),
                    fGAliceFileName(galicefilename)
{
  //constructor, 
  //Defaults:
  //  trackfilename = "AliITStracksV2.root"
  //  clusterfilename = "AliITSclustersV2.root"
  //  galicefilename = "galice.root"
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();
  fIsRead = kFALSE;
}
/********************************************************************/

AliHBTReaderITSv2::
 AliHBTReaderITSv2(TObjArray* dirs, const Char_t* trackfilename, 
                   const Char_t* clusterfilename, const Char_t* galicefilename)
                  : AliHBTReader(dirs),
                    fTrackFileName(trackfilename),fClusterFileName(clusterfilename),
                    fGAliceFileName(galicefilename)
{
  //constructor, 
  //Defaults:
  //  trackfilename = "AliITStracksV2.root"
  //  clusterfilename = "AliITSclustersV2.root"
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
    if(Read(fParticles,fTracks))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }

   return fParticles->GetEvent(n);
 }
/********************************************************************/

AliHBTEvent* AliHBTReaderITSv2::GetTrackEvent(Int_t n)
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

Int_t AliHBTReaderITSv2::GetNumberOfPartEvents()
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
Int_t AliHBTReaderITSv2::GetNumberOfTrackEvents()
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
/********************************************************************/
Int_t AliHBTReaderITSv2::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
 Int_t Nevents = 0; //number of events found in given directory
 Int_t Ndirs; //number of the directories to be read
 Int_t Ntracks; //number of tracks in current event
 Int_t currentdir = 0; //number of events in the current directory 
 Int_t totalNevents = 0; //total number of events read from all directories up to now
 register Int_t i; //iterator
 
 TFile *aTracksFile;//file with tracks
 TFile *aClustersFile;//file with clusters
 TFile *aGAliceFile;//file name with galice

// AliITStrackerV2 *tracker; // ITS tracker - used for cooking labels
 TTree *tracktree; // tree for tracks
 
 Double_t xk;
 Double_t par[5]; //Kalman track parameters
 Float_t phi, lam, pt;//angles and transverse momentum
 Int_t label; //label of the current track

 char tname[100]; //buffer for tree name
 AliITStrackV2 *iotrack= new AliITStrackV2(); //buffer track for reading data from tree

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
    
    AliITSgeom *geom=(AliITSgeom*)aClustersFile->Get("AliITSgeom");
    if (!geom) 
     { 
       Error("Read","Can't get the ITS geometry!"); 
       currentdir++;
       continue;
     }

    for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
     {
       cout<<"Reading Event "<<currentEvent<<endl;
       
       aGAliceFile->cd();
       gAlice->GetEvent(currentEvent);

       aClustersFile->cd();
       sprintf(tname,"TreeT_ITS_%d",currentEvent);
       
       tracktree=(TTree*)aTracksFile->Get(tname);
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
            
          AliHBTParticle* part = new AliHBTParticle(*p);
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
          accepted++;
        }//end of loop over tracks in the event
        
       aTracksFile->Delete(tname);
       aTracksFile->Delete("tracks");
//       delete tracker;
       
       totalNevents++;
       cout<<"all: "<<i<<"   accepted: "<<accepted<<"   tpc faults: "<<tpcfault<<"   its faults: "<<itsfault<<endl;
     
     }//end of loop over events in current directory
    CloseFiles(aTracksFile,aClustersFile,aGAliceFile);     
    currentdir++;
   }while(currentdir < Ndirs);//end of loop over directories specified in fDirs Obj Array

 delete iotrack;
 fIsRead = kTRUE;
 return 0;
}

/********************************************************************/
Int_t AliHBTReaderITSv2::OpenFiles
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
  
void AliHBTReaderITSv2::CloseFiles(TFile*& tracksFile, TFile*& clustersFile, TFile*& gAliceFile)
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


