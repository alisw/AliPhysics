#include "AliHBTReaderPPprod.h"

#include <Riostream.h>
#include <Riostream.h>
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


ClassImp(AliHBTReaderPPprod)

AliHBTReaderPPprod::
 AliHBTReaderPPprod(const Char_t* trackfilename,const Char_t* clusterfilename,
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

AliHBTReaderPPprod::~AliHBTReaderPPprod()
 {
   delete fParticles;
   delete fTracks;
 }
/********************************************************************/

AliHBTEvent* AliHBTReaderPPprod::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
   if (!fIsRead) Read(fParticles,fTracks);
   return fParticles->GetEvent(n);
 }
/********************************************************************/
AliHBTEvent* AliHBTReaderPPprod::GetTrackEvent(Int_t n)
 {
 //returns Nth event with reconstructed tracks
   if (!fIsRead) Read(fParticles,fTracks);
   return fTracks->GetEvent(n);
 }
/********************************************************************/

Int_t AliHBTReaderPPprod::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead) Read(fParticles,fTracks);
   return fParticles->GetNumberOfEvents();
 }

/********************************************************************/
Int_t AliHBTReaderPPprod::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (!fIsRead) Read(fParticles,fTracks);
  return fTracks->GetNumberOfEvents();
 }
/********************************************************************/


Int_t AliHBTReaderPPprod::Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Int_t i; //iterator and some temprary values
  Int_t Nevents;
  if (!particles) //check if an object is instatiated
   {
     Error("AliHBTReaderPPprod::Read"," particles object must instatiated before passing it to the reader");
   }
  if (!tracks)  //check if an object is instatiated
   {
     Error("AliHBTReaderPPprod::Read"," tracks object must instatiated before passing it to the reader");
   }
  particles->Reset();//clear runs == delete all old events
  tracks->Reset();
    
  if( (i=OpenFiles()) )
   {
     Error("AliHBTReaderPPprod::Read","Exiting due to problems with opening files. Errorcode %d",i);
     return i;
   }
  
  AliGoodTracksPP *goodTPCTracks = new AliGoodTracksPP(fGoodTPCTracksFileName);
  if (!goodTPCTracks)
   {
     Error("AliHBTReaderPPprod::Read","Exiting due to problems with opening files. Errorcode %d",i);
     return 1;
   }

  Nevents = 100;

  fClustersFile->cd();
  AliTPCParam *TPCParam= (AliTPCParam*)fClustersFile->Get("75x40_100x60");
  if (!TPCParam) 
    { 
     Error("AliHBTReaderPPprod::Read","TPC parameters have not been found !\n");
     return 1;
    }

  TObjArray *tarray = new TObjArray(5000);
  tarray->SetOwner();//  this causes memory leak, but in some cases deleting is infinite loop
  
 
  for(Int_t currentEvent =0; currentEvent<Nevents;currentEvent++)
   {
     cout<<"Reading Event "<<currentEvent<<endl;
     /**************************************/
      /**************************************/
       /**************************************/ 
         fTracksFile->cd();
         Char_t  treename[100];
         sprintf(treename,"TreeT_TPC_%d",currentEvent);
   
         TTree *tracktree=0;
         
         tracktree=(TTree*)fTracksFile->Get(treename);
         if (!tracktree) 
          {
            Error("AliHBTReaderPPprod::Read","Can't get a tree with TPC tracks !\n"); 
            return 1;
          }
   
         TBranch *trackbranch=tracktree->GetBranch("tracks");
         if (!trackbranch) 
          {
            Error("AliHBTReaderPPprod::Read","Can't get a branch with TPC tracks !\n"); 
            return 2;
          }
         Int_t NTPCtracks=(Int_t)tracktree->GetEntries();
         cout<<"Found "<<NTPCtracks<<" TPC tracks.\n";
         //Copy tracks to array
         
         AliTPCtrack *iotrack=0;
         
         fClustersFile->cd();
       //AliTPCtracker *tracker=new AliTPCtracker(TPCParam,currentEvent);//I.B.
       AliTPCtracker *tracker=new AliTPCtracker(TPCParam);               //I.B.
       tracker->SetEventNumber(currentEvent);                            //I.B.
         if (!tracker) 
          {
            Error("AliHBTReaderPPprod::Read","Can't get a tracker !\n"); 
            return 3;
          }
         //tracker->LoadInnerSectors(); //I.B.
         //tracker->LoadOuterSectors(); //I.B.
         tracker->LoadClusters();     //I.B.
   
         for (i=0; i<NTPCtracks; i++)
          {
            iotrack=new AliTPCtrack;
            trackbranch->SetAddress(&iotrack);
            tracktree->GetEvent(i);
            tracker->CookLabel(iotrack,0.1);
            tarray->AddLast(iotrack);
          }
        
         
         fTracksFile->Delete(treename);//delete tree from memmory (and leave untached on disk)- we do not need it any more
         fTracksFile->Delete("tracks");
         
         delete tracker;
         
         tracker = 0x0;
         trackbranch = 0x0;
         tracktree = 0x0;

         Int_t & ngood = goodTPCTracks->fGoodInEvent[currentEvent];
   
         Double_t xk;
         Double_t par[5];
         Float_t phi, lam, pt;
         Int_t label;
         Bool_t found;
   
         for (i=0; i<ngood; i++)
          { 
            const struct GoodTrack & gt = goodTPCTracks->GetTrack(currentEvent,i);
            
            if(Pass(gt.code)) continue;
            
            label = gt.lab;
            found = kFALSE; //guard in case we don't find track with such a label
            for (Int_t j=0;j<NTPCtracks;j++)
              {
                iotrack = (AliTPCtrack*)tarray->At(j);
                if (iotrack->GetLabel() == label) 
                  {
                    found = kTRUE;
                    break;
                  }
              }  
            if(!found) 
              {
                Warning("Read",
                "Sth is going wrong with tracks - there is no TPC track corresponding to goodtrack.\nGood tack label %d",label);
                continue; //put comunicate on the screen and continue loop
              }
        
            Double_t mass = TDatabasePDG::Instance()->GetParticle(gt.code)->Mass();
            Double_t pEtot = TMath::Sqrt(gt.px*gt.px + gt.py*gt.py + gt.pz*gt.pz + mass*mass);
            
            AliHBTParticle* part = new AliHBTParticle(gt.code, gt.px, gt.py, gt.pz, pEtot, gt.x, gt.y, gt.z, 0.0);
            if(Pass(part)) continue;
            
         
            iotrack->PropagateTo(gt.x);
            iotrack->GetExternalParameters(xk,par);
            phi=TMath::ASin(par[2]) + iotrack->GetAlpha();
            if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
            if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
            lam=par[3]; 
            pt=1.0/TMath::Abs(par[4]);
            
            Double_t tpx = pt * TMath::Cos(phi); //track x coordinate of momentum
            Double_t tpy = pt * TMath::Sin(phi); //track x coordinate of momentum
            Double_t tpz = pt * lam;
            
            Double_t tEtot = TMath::Sqrt( tpx*tpx + tpy*tpy + tpz*tpz + mass*mass);
            
            AliHBTParticle* track = new AliHBTParticle(gt.code, tpx, tpy , tpz, tEtot, 0., 0., 0., 0.);
            if(Pass(track)) continue;
            
            particles->AddParticle(currentEvent,part);
            tracks->AddParticle(currentEvent,track);

          }
         tarray->Clear();
       /**************************************/
      /**************************************/
     /**************************************/  
   }
  
  CloseFiles();
  delete tarray;
  delete goodTPCTracks;
  fIsRead = kTRUE;
  return 0;
 }

/********************************************************************/
Int_t AliHBTReaderPPprod::OpenFiles()
{
 
   fTracksFile = 0;
   fTracksFile=TFile::Open("AliTPCtracks.root");
   if (!fTracksFile->IsOpen()) 
     {
       Error("AliHBTReaderPPprod::OpenFiles","Can't open AliTPCtracks.root");
       return 1;
     }
   
   fClustersFile = 0;
   
   fClustersFile=TFile::Open("AliTPCclusters.root");
   if (!fClustersFile->IsOpen()) 
    {
      Error("AliHBTReaderPPprod::OpenFiles","Can't open AliTPCclusters.root");
      return 2;
    }

    

 return 0; 
}
  


/********************************************************************/
  
void AliHBTReaderPPprod::CloseFiles()
{
  fTracksFile->Close();
  fClustersFile->Close();
}

/********************************************************************/

/********************************************************************/
/********************************************************************/
/********************************************************************/


AliGoodTracksPP::~AliGoodTracksPP()
{
 delete [] fGoodInEvent;
 for (Int_t i = 0;i<fNevents;i++)
   delete [] fData[i];
 delete [] fData;
}
/********************************************************************/
AliGoodTracksPP::AliGoodTracksPP(const TString& infilename)
{

  fNevents = 100;
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
    fData[i] = new struct GoodTrack[500];
   }
  Float_t tmp;
  Int_t evno;
  while(in>>evno)
   {
    if(fGoodInEvent[evno]>=500)
     {
      cerr<<"AliGoodTracksPP::AliGoodTracksPP() : Not enough place in the array\n";
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
    
    in>>tmp;
    in>>tmp;
    in>>tmp;
    
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
 cout<<"AliGoodTracksPP::AliGoodTracksPP()  ....  Done\n";
}



const GoodTrack& AliGoodTracksPP::GetTrack(Int_t event, Int_t n) const
 {
  
  if( (event>fNevents) || (event<0))
   {
     gAlice->Fatal("AliGoodTracksPP::GetTrack","No such Event %d",event);
   }
  if( (n>fGoodInEvent[event]) || (n<0))
   {
     gAlice->Fatal("AliGoodTracksPP::GetTrack","No such Good TPC Track %d",n);
   }
  return fData[event][n];

 }
