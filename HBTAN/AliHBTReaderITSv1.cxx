
#include "AliHBTReader.h"

ClassImp(AliHBTReaderITSv1)

AliHBTReaderITSv1::AliHBTReaderITSv1(const Char_t* goodtracksfilename):
                 fGoodITSTracksFileName(goodtracksfilename)

 {
     fParticles = new AliHBTRun();
     fTracks    = new AliHBTRun();
     fIsRead = kFALSE;
 }
AliHBTReaderITSv1::AliHBTReaderITSv1()
{
   delete fParticles;
   delete fTracks;
}
AliHBTEvent* AliHBTReaderITSv1::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
   if (!fIsRead) Read(fParticles,fTracks);
   return fParticles->GetEvent(n);
 }
/********************************************************************/
AliHBTEvent* AliHBTReaderITSv1::GetTrackEvent(Int_t n)
 {
 //returns Nth event with reconstructed tracks
   if (!fIsRead) Read(fParticles,fTracks);
   return fTracks->GetEvent(n);
 }
/********************************************************************/

Int_t AliHBTReaderITSv1::GetNumberOfPartEvents()
 {
 //returns number of events of particles
   if (!fIsRead) Read(fParticles,fTracks);
   return fParticles->GetNumberOfEvents();
 }

/********************************************************************/
Int_t AliHBTReaderITSv1::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (!fIsRead) Read(fParticles,fTracks);
  return fTracks->GetNumberOfEvents();
 }
/********************************************************************/

Int_t AliHBTReaderITSv1::Read(AliHBTRun* particles, AliHBTRun *tracks)
 {
   AliGoodTracksITSv1 *goodITStracks = new AliGoodTracksITSv1(fGoodITSTracksFileName);
   
   if (!goodITStracks)
    {
     Error("AliHBTReaderITSv1::Read","Exiting due to problems with opening files. Errorcode %d",i);
     return 1;
    }
   for (Int_t currentEvent = 0; currentEvent < goodITStracks->fNevents; currentEvent++)
    {
      for(Int_t i =0; i<goodITStracks->fGoodInEvent[currentEvent]; i++)
       {
         const struct GoodTrackITSv1 & gt = goodTPCTracks->GetTrack(currentEvent,i);
         
         if(Pass(gt.code)) continue;
         
         Double_t mass = TDatabasePDG::Instance()->GetParticle(gt.code)->Mass();
         Double_t pEtot = TMath::Sqrt(gt.px*gt.px + gt.py*gt.py + gt.pz*gt.pz + mass*mass);
         
         AliHBTParticle* part = new AliHBTParticle(gt.code, gt.px, gt.py, gt.pz, pEtot, gt.x, gt.y, gt.z, 0.0);
         if(Pass(part)) { delete part; continue;}
         
         Double_t tEtot = TMath::Sqrt( gt.pxg*gt.pxg + gt.pyg*gt.pyg + gt.pzg*gt.pzg + mass*mass);
         
         AliHBTParticle* track = new AliHBTParticle(gt.code, gt.pxg, gt.pyg , gt.pzg, tEtot, 0., 0., 0., 0.);
         if(Pass(track)) { delete  track;continue;}
            
         particles->AddParticle(currentEvent,part);
         tracks->AddParticle(currentEvent,track);
       }
    }
  delete goodITStracks;
  fIsRead = kTRUE;
  return 0;
 
 }
/********************************************************************/

/********************************************************************/
/********************************************************************/
/********************************************************************/


AliGoodTracksITSv1::~AliGoodTracksITSv1()
{
 delete [] fGoodInEvent;
 for (Int_t i = 0;i<fNevents;i++)
   delete [] fData[i];
 delete [] fData;
}
/********************************************************************/
AliGoodTracksITSv1::AliGoodTracksITSv1(const TString& infilename)
{

  cout<<"AliGoodTracksITSv1::AliGoodTracksITSv1()  ....\n";


  fNevents = 0;
  Int_t maxevents = 500;
  
  ifstream in(infilename.Data());

  if(!in)
    {
      cerr<<"Can not open file with Good ITSv1 Tracks named:"<<infilename.Data()<<endl;
      delete this;
      return;
    }

  
  fGoodInEvent = new Int_t[maxevents];
  fData = new struct GoodTrack* [maxevents];

  Int_t i;
  Int_t lastevent = 0;
  fGoodInEvent[0] =0;
  fData[0] = new struct GoodTrack[50000];

  Int_t evno;
  while(in>>evno)
   {
    if(lastevent>evno)
     {
       for(i = lastevent;i<=evno;i++)
        {
            fGoodInEvent[i] =0;
            fData[i] = new struct GoodTrack[50000];
        }
      lastevent = evno;
     }
    if(fGoodInEvent[evno]>=50000)
     {
      cerr<<"AliGoodTracksITSv1::AliGoodTracksITSv1() : Not enough place in the array\n";
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
 fNevents = evno+1;
 in.close();
 cout<<"AliGoodTracksITSv1::AliGoodTracksITSv1()  ....  Done\n";
}



const GoodTrack& AliGoodTracksITSv1::GetTrack(Int_t event, Int_t n) const
 {
  
  if( (event>fNevents) || (event<0))
   {
     gROOT->Fatal("AliGoodTracksITSv1::GetTrack","No such Event %d",event);
   }
  if( (n>fGoodInEvent[event]) || (n<0))
   {
     gROOT->Fatal("AliGoodTracksITSv1::GetTrack","No such Good TPC Track %d",n);
   }
  return fData[event][n];

 }
