#include "AliHBTReaderTPC.h"

#include <iostream.h>
#include <fstream.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>

#include <AliTPCtrack.h>
#include <AliTPCParam.h>
#include <AliTPCtracker.h>

#include "AliRun.h"
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
                 const Char_t* goodtracksfilename,const Char_t* galicefilename):
                 fTrackFileName(trackfilename),fClusterFileName(clusterfilename),
                 fGAliceFileName(galicefilename),
                 fGoodTPCTracksFileName(goodtracksfilename)
{
  //constructor, only file names are set
  //Defaults:
  //  trackfilename = "AliTPCtracks.root"
  //  clusterfilename = "AliTPCclusters.root"
  //  goodtracksfilename = "good_tracks_tpc"
  //  galicefilename = ""  - this means: Do not open gAlice file - 
  //                         just leave the global pointer untached
  
  fParticles = new AliHBTRun();
  fTracks    = new AliHBTRun();

  fTracksFile   = 0x0;  //files are opened during reading only
  fClustersFile = 0x0;
  
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
   if (!fIsRead) Read(fParticles,fTracks);
   return fParticles->GetEvent(n);
 }
/********************************************************************/
AliHBTEvent* AliHBTReaderTPC::GetTrackEvent(Int_t n)
 {
 //returns Nth event with reconstructed tracks
   if (!fIsRead) Read(fParticles,fTracks);
   return fTracks->GetEvent(n);
 }
/********************************************************************/

Int_t AliHBTReaderTPC::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead) Read(fParticles,fTracks);
   return fParticles->GetNumberOfEvents();
 }

/********************************************************************/
Int_t AliHBTReaderTPC::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (!fIsRead) Read(fParticles,fTracks);
  return fTracks->GetNumberOfEvents();
 }
/********************************************************************/


Int_t AliHBTReaderTPC::Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Int_t i; //iterator and some temprary values
  Int_t Nevents;
  if (!particles) //check if an object is instatiated
   {
     Error("AliHBTReaderTPC::Read"," particles object must instatiated before passing it to the reader");
   }
  if (!tracks)  //check if an object is instatiated
   {
     Error("AliHBTReaderTPC::Read"," tracks object must instatiated before passing it to the reader");
   }
  particles->Reset();//clear runs == delete all old events
  tracks->Reset();
    
  if( (i=OpenFiles()) )
   {
     Error("AliHBTReaderTPC::Read","Exiting due to problems with opening files. Errorcode %d",i);
     return i;
   }
  
  AliGoodTracks *goodTPCTracks = new AliGoodTracks(fGoodTPCTracksFileName);
  if (!goodTPCTracks)
   {
     Error("AliHBTReaderTPC::Read","Exiting due to problems with opening files. Errorcode %d",i);
     return 1;
   }
  
    
  if (gAlice->TreeE())//check if tree E exists
   {
    Nevents = (Int_t)gAlice->TreeE()->GetEntries();//if yes get number of events in gAlice
    cout<<"Found "<<Nevents<<endl;
   }
  else
   {//if not return an error
     Error("AliHBTReaderPPprod::Read","Can not find Header tree (TreeE) in gAlice");
     return 1;
   }
  
  fClustersFile->cd();//set cluster file active 
  AliTPCParam *TPCParam= (AliTPCParam*)fClustersFile->Get("75x40_100x60");
  if (!TPCParam) 
    { 
     Error("AliHBTReaderTPC::Read","TPC parameters have not been found !\n");
     return 1;
    }

  TObjArray *tarray = new TObjArray(5000); //cotainer for tpc tracks
  tarray->SetOwner(); //set the ownership of the objects it contains
                      //when array is is deleted or cleared all objects 
                      //that it contains are deleted
  
  for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)//loop over all events
   {
     cout<<"Reading Event "<<currentEvent<<endl;
     /**************************************/
      /**************************************/
       /**************************************/ 
         fTracksFile->cd();//set track file active
         
         Char_t  treename[100];
         sprintf(treename,"TreeT_TPC_%d",currentEvent);//prepare name of the tree
   
         TTree *tracktree=0;
         
         tracktree=(TTree*)fTracksFile->Get(treename);//get the tree 
         if (!tracktree) //check if we got the tree
          {//if not return with error
            Error("AliHBTReaderTPC::Read","Can't get a tree with TPC tracks !\n"); 
            return 1;
          }
   
         TBranch *trackbranch=tracktree->GetBranch("tracks");//get the branch with tracks
         if (!trackbranch) ////check if we got the branch
          {//if not return with error
            Error("AliHBTReaderTPC::Read","Can't get a branch with TPC tracks !\n"); 
            return 2;
          }
         Int_t NTPCtracks=(Int_t)tracktree->GetEntries();//get number of TPC tracks 
         cout<<"Found "<<NTPCtracks<<" TPC tracks.\n";
         //Copy tracks to array
         
         AliTPCtrack *iotrack=0;
         
         fClustersFile->cd();//set cluster file active 
         AliTPCtracker *tracker = new AliTPCtracker(TPCParam,currentEvent);//create the tacker for this event
         if (!tracker) //check if it has created succeffuly
          {//if not return with error
            Error("AliHBTReaderTPC::Read","Can't get a tracker !\n"); 
            return 3;
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
         
         fTracksFile->Delete(treename);//delete tree from memmory (and leave untached on disk)- we do not need it any more
         fTracksFile->Delete("tracks");//delete branch from memmory
         delete tracker; //delete tracker
         
         tracker = 0x0;
         trackbranch = 0x0;
         tracktree = 0x0;

         Int_t & ngood = goodTPCTracks->fGoodInEvent[currentEvent]; //number of good tracks in the current event
   
         Double_t xk;
         Double_t par[5];
         Float_t phi, lam, pt;//angles and transverse momentum
         Int_t label; //label of the current track
         Bool_t found; //flag indicated wether we managed to match good_tpc_track with track
   
         for (i=0; i<ngood; i++) //loop over all good tracks
          { 
            const struct GoodTrack & gt = goodTPCTracks->GetTrack(currentEvent,i); //get ith goog track
            
            if(Pass(gt.code)) continue; //check if we are intersted with particles of this type 
                                        //if not take next partilce
            
            label = gt.lab;
            found = kFALSE; //guard in case we don't find track with such a label
            for (Int_t j=0;j<NTPCtracks;j++)//lopp over all tpc tracks
              {
                iotrack = (AliTPCtrack*)tarray->At(j);
                if (iotrack->GetLabel() == label) //if the label is the same 
                  {
                    found = kTRUE; //we found the track
                    break;
                  }
              }  
            if(!found) //check if we found the track
              {
                Warning("Read",
                "Sth is going wrong with tracks - there is no TPC track corresponding to goodtrack.\nGood tack label %d",label);
                continue; //put comunicate on the screen and continue loop
              }
        
            Double_t mass = TDatabasePDG::Instance()->GetParticle(gt.code)->Mass();//CMS mass of this particle 
            Double_t pEtot = TMath::Sqrt(gt.px*gt.px + gt.py*gt.py + gt.pz*gt.pz + mass*mass); //particle total energy
            
            AliHBTParticle* part = new AliHBTParticle(gt.code, gt.px, gt.py, gt.pz, pEtot, gt.x, gt.y, gt.z, 0.0);
            if(Pass(part)) { delete part; continue;}//check if meets all criteria of any of our cuts
                                                    //if it does not delete it and take next good track
         
            iotrack->PropagateTo(gt.x);
            iotrack->GetExternalParameters(xk,par);     //get properties of the track
            phi=TMath::ASin(par[2]) + iotrack->GetAlpha(); 
            if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
            if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
            lam=par[3]; 
            pt=1.0/TMath::Abs(par[4]);
            
            Double_t tpx = pt * TMath::Cos(phi); //track x coordinate of momentum
            Double_t tpy = pt * TMath::Sin(phi); //track y coordinate of momentum
            Double_t tpz = pt * lam; //track z coordinate of momentum
            
            Double_t tEtot = TMath::Sqrt( tpx*tpx + tpy*tpy + tpz*tpz + mass*mass);//total energy of the track
            
            AliHBTParticle* track = new AliHBTParticle(gt.code, tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
            if(Pass(track)) { delete  track;continue;}//check if meets all criteria of any of our cuts
                                                      //if it does not delete it and take next good track
            
            particles->AddParticle(currentEvent,part);//put track and particle on the run
            tracks->AddParticle(currentEvent,track);

          }
         tarray->Clear(); //clear the array
    
       /**************************************/
      /**************************************/
     /**************************************/  
   }
  
  //save environment (resouces) --
  //clean your place after the work
  CloseFiles(); 
  delete tarray;
  delete goodTPCTracks;
  fIsRead = kTRUE;
  return 0;
 }

/********************************************************************/
Int_t AliHBTReaderTPC::OpenFiles()
{
 //opens all the files
   fTracksFile = 0;
   fTracksFile=TFile::Open(fTrackFileName.Data());
   if (!fTracksFile->IsOpen()) 
     {
       Error("AliHBTReaderTPC::OpenFiles","Can't open file with tacks named ",fTrackFileName.Data());
       return 1;
     }
   
   fClustersFile = 0;
   
   fClustersFile=TFile::Open(fClusterFileName.Data());
   if (!fClustersFile->IsOpen()) 
    {
      Error("AliHBTReaderTPC::OpenFiles","Can't open file with TPC clusters named ",fClusterFileName.Data());
      return 2;
    }

 return 0; 
}
  


/********************************************************************/
  
void AliHBTReaderTPC::CloseFiles()
{
  //closes the files
  fTracksFile->Close();
  fClustersFile->Close();
}

/********************************************************************/

/********************************************************************/
/********************************************************************/
/********************************************************************/


AliGoodTracks::~AliGoodTracks()
{
//destructor
 delete [] fGoodInEvent;
 for (Int_t i = 0;i<fNevents;i++)
   delete [] fData[i];
 delete [] fData;
}
/********************************************************************/
AliGoodTracks::AliGoodTracks(const TString& infilename)
{

  cout<<"AliGoodTracks::AliGoodTracks()  ....\n";
  if(!gAlice) 
    {
      cerr<<"There is no gAlice"<<endl;
      delete this;
      return;
    }
  
  if (!gAlice->TreeE())
   {
     cerr<<"Can not find Header tree (TreeE) in gAlice"<<endl;
     delete this;
     return;
   }
   
  fNevents = (Int_t)gAlice->TreeE()->GetEntries();
  //fNevents = 100;
  cout<<fNevents<<" FOUND"<<endl;
  ifstream in(infilename.Data());

  if(!in)
    {
      cerr<<"Can not open file with Good TPC Tracks named:"<<infilename.Data()<<endl;
      delete this;
      return;
    }

  
  fGoodInEvent = new Int_t[fNevents];
  fData = new struct GoodTrack* [fNevents];

  Int_t i;
  for( i = 0;i<fNevents;i++)
   {
    fGoodInEvent[i] =0;
    fData[i] = new struct GoodTrack[50000];
   }

  Int_t evno;
  while(in>>evno)
   {
    if(fGoodInEvent[evno]>=50000)
     {
      cerr<<"AliGoodTracks::AliGoodTracks() : Not enough place in the array\n";
      continue;
     }
    in>>fData[evno][fGoodInEvent[evno]].lab;
    in>>fData[evno][fGoodInEvent[evno]].code;
    in>>fData[evno][fGoodInEvent[evno]].px;
    in>>fData[evno][fGoodInEvent[evno]].py;
    in>>fData[evno][fGoodInEvent[evno]].pz;
    in>>fData[evno][fGoodInEvent[evno]].x;
    in>>fData[evno][fGoodInEvent[evno]].y;
    in>>fData[evno][fGoodInEvent[evno]].z;
    
 /* cout<<evno<<" ";
  cout<<fData[evno][fGoodInEvent[evno]].lab;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].code;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].px;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].py;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].pz;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].x;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].y;
  cout<<" ";cout<<fData[evno][fGoodInEvent[evno]].z;
  cout<<"\n";
 */ 
  fGoodInEvent[evno]++;
 }
 in.close();
 cout<<"AliGoodTracks::AliGoodTracks()  ....  Done\n";
}



const GoodTrack& AliGoodTracks::GetTrack(Int_t event, Int_t n) const
 {
  
  if( (event>fNevents) || (event<0))
   {
     gAlice->Fatal("AliGoodTracks::GetTrack","No such Event %d",event);
   }
  if( (n>fGoodInEvent[event]) || (n<0))
   {
     gAlice->Fatal("AliGoodTracks::GetTrack","No such Good TPC Track %d",n);
   }
  return fData[event][n];

 }
