// $Id$

//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliJetParticlesReader
//
// loizides@ikf.uni-frankfurt.de
///////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TClass.h>
#include <TString.h>
#include <TObjString.h>

#include "AliJetEventParticles.h"
#include "AliJetParticlesReader.h"

 
ClassImp(AliJetParticlesReader)

AliJetParticlesReader::AliJetParticlesReader()
  : TNamed(),
    fEventParticles(0),
    fOwner(kTRUE),
    fDirs(0),
    fCurrentEvent(0),
    fCurrentDir(0),
    fNEventsRead(0),
    fFirst(0),
    fLast(0),
    fPtMin(0),fPtMax(1000),
    fEtaMin(-1),fEtaMax(1),
    fPhiMin(0),fPhiMax(2*TMath::Pi())
{
}

AliJetParticlesReader::AliJetParticlesReader(TObjArray *dirs)
  : TNamed(),
    fEventParticles(0),
    fOwner(kTRUE),
    fDirs(dirs),
    fCurrentEvent(0),
    fCurrentDir(0),
    fNEventsRead(0),
    fFirst(0),
    fLast(0),
    fPtMin(0),fPtMax(1000),
    fEtaMin(-1),fEtaMax(1),
    fPhiMin(0),fPhiMax(2*TMath::Pi())
{
}

AliJetParticlesReader::~AliJetParticlesReader()
{
  if((fOwner) && (fEventParticles)) delete fEventParticles;
}

Int_t AliJetParticlesReader::Next()
{
  //moves to next event

  //if asked to read up to event nb. fLast, 
  //and it is overcome, report no more events
  if ((fNEventsRead > fLast) && (fLast > 0) ) return kFALSE;
  
  do //if asked to read from event fFirst, rewind to it
    {
      if ( ReadNext() == kFALSE) 
	return kFALSE; //if no more evets, return it
    } while (fNEventsRead < fFirst);
   
  //here we have event

  return kTRUE;
}

TString& AliJetParticlesReader::GetDirName(Int_t entry)
{
  //returns directory name of entry to read

  TString* retval;//return value
  if (fDirs == 0)
    {
     retval = new TString(".");
     return *retval;
   }

  if ((entry>fDirs->GetEntries()) || (entry<0))
                 //if out of bounds return empty string
                 //note that entry==0 is accepted even if array is empty (size=0)
    {
      Error("GetDirName","Entry out of bounds");
      retval = new TString();
      return *retval;
    }

  if (fDirs->GetEntries() == 0)
    { 
      retval = new TString(".");
      return *retval;
    }

  TClass *objclass = fDirs->At(entry)->IsA();
  TClass *stringclass = TObjString::Class();

  TObjString *dir = (TObjString*)objclass->DynamicCast(stringclass,fDirs->At(entry));
  if(dir == 0)
   {
     Error("GetDirName","Object in TObjArray is not a TObjString");
     retval = new TString();
     return *retval;
   }

  //Info("GetDirName","Returned ok %s",dir->String().Data());
  return dir->String();
}


#if 0
/*************************************************************************************/

void AliHBTReader::AddParticleCut(AliHBTParticleCut* cut)
{
 //sets the new cut 
 
  if (!cut) //if cut is NULL return with error
   {
    Error("AddParticleType","NULL pointers are not accepted any more.\nIf You want to accept all particles of this type, set an empty cut ");
    return;
   }
  AliHBTParticleCut *c = (AliHBTParticleCut*)cut->Clone();
  fCuts->Add(c);
}
/********************************************************************/

AliHBTEvent* AliHBTReader::GetParticleEvent(Int_t n)
 {
 //returns Nth event with simulated particles
  if (ReadsParticles() == kFALSE)
   {
     Error("GetParticleEvent","This reader is not able to provide simulated particles.");
     return 0;
   } 
   
  if (!fIsRead) 
   { 
    if (ReadsParticles() && (fParticles == 0x0)) fParticles = new AliHBTRun();
    if (ReadsTracks() && (fTracks == 0x0)) fTracks = new AliHBTRun();
    
    if (Read(fParticles,fTracks))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }
    else fIsRead = kTRUE;
   }
  return fParticles->GetEvent(n);
 }
/********************************************************************/

AliHBTEvent* AliHBTReader::GetTrackEvent(Int_t n)
 {
 //returns Nth event with reconstructed tracks
  if (ReadsTracks() == kFALSE)
   {
     Error("GetTrackEvent","This reader is not able to provide recosntructed tracks.");
     return 0;
   } 
  if (!fIsRead) 
   {
    if (ReadsParticles() && (fParticles == 0x0)) fParticles = new AliHBTRun();
    if (ReadsTracks() && (fTracks == 0x0)) fTracks = new AliHBTRun();
    
    if(Read(fParticles,fTracks))
     {
       Error("GetTrackEvent","Error in reading");
       return 0x0;
     }
    else fIsRead = kTRUE;
   }
  return fTracks->GetEvent(n);
 }
/********************************************************************/

Int_t AliHBTReader::GetNumberOfPartEvents()
 {
 //returns number of events of particles
  if (ReadsParticles() == kFALSE)
   {
     Error("GetNumberOfPartEvents","This reader is not able to provide simulated particles.");
     return 0;
   } 
   
  if (!fIsRead) 
   {
    if (ReadsParticles() && (fParticles == 0x0)) fParticles = new AliHBTRun();
    if (ReadsTracks() && (fTracks == 0x0)) fTracks = new AliHBTRun();
    
    if (Read(fParticles,fTracks))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
    else fIsRead = kTRUE;
   }
   return fParticles->GetNumberOfEvents();
 }
/********************************************************************/
 
Int_t AliHBTReader::GetNumberOfTrackEvents()
 {
 //returns number of events of tracks
  if (ReadsTracks() == kFALSE)
   {
     Error("GetNumberOfTrackEvents","This reader is not able to provide recosntructed tracks.");
     return 0;
   } 
  if (!fIsRead)
   {
     if (ReadsParticles() && (fParticles == 0x0)) fParticles = new AliHBTRun();
     if (ReadsTracks() && (fTracks == 0x0)) fTracks = new AliHBTRun();
     
     if(Read(fParticles,fTracks))
      {
        Error("GetNumberOfTrackEvents","Error in reading");
        return 0;
      }
     else fIsRead = kTRUE;
   }
  return fTracks->GetNumberOfEvents();
 }
/********************************************************************/

Int_t AliHBTReader::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("Read","");
  
  if ( ReadsParticles() && (particles == 0x0) ) //check if an object is instatiated
   {
     Error("Read"," particles object must be instatiated before passing it to the reader");
     return 1;
   }
  if ( ReadsTracks() && (tracks == 0x0) )  //check if an object is instatiated
   {
     Error("Read"," tracks object must be instatiated before passing it to the reader");
     return 1;
   }
   
  if (ReadsParticles()) particles->Reset();//clear runs == delete all old events
  if (ReadsTracks()) tracks->Reset();
  
  Rewind();
  
  Int_t i = 0;
  while(Next() == kFALSE)
   {
     if (ReadsTracks()) tracks->SetEvent(i,fTracksEvent);
     if (ReadsParticles()) particles->SetEvent(i,fParticlesEvent);
     i++;
   }
  return 0;
}      
/*************************************************************************************/

Bool_t AliHBTReader::Pass(AliHBTParticle* p)
{
 //Method examines whether particle meets all cut and particle type criteria
  
   if(p==0x0)//of corse we not pass NULL pointers
    {
     Warning("Pass()","No Pasaran! We never accept NULL pointers");
     return kTRUE;
    }
   //if no particle is specified, we pass all particles
   //excluding NULL pointers, of course
  if ( fCuts->GetEntriesFast() == 0 ) return kFALSE; //if no cut specified accept all particles
  for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliHBTParticleCut &cut = *((AliHBTParticleCut*)fCuts->At(i));
     if(!cut.Pass(p)) return kFALSE;  //accepted
   }
   
  return kTRUE;//not accepted
}
/*************************************************************************************/

Bool_t  AliHBTReader::Pass(Int_t pid)
{
//this method checks if any of existing cuts accepts this pid particles
//or any cuts accepts all particles

 if(pid == 0)
  return kTRUE;

 if ( fCuts->GetEntriesFast() == 0 ) return kFALSE; //if no cut specified accept all particles
  
 for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliHBTParticleCut &cut = *((AliHBTParticleCut*)fCuts->At(i));
     //if some of cuts accepts all particles or some accepts particles of this type, accept
     if ( (cut.GetPID() == 0) || (cut.GetPID() == pid) ) return kFALSE; 
   }
 return kTRUE;
}

#endif
