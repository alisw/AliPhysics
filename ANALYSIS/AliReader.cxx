/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliReader
//
// Reader Base class 
// Reads particles and tracks and
// puts them to the AliAOD objects and eventually, if needed, buffers AliAODs in AliAODRun(s)
//
// User loops over events calling method Next. In case of success this method returns 0.
// In case of error or if there is no more events to read, non-0 value is returned
//
// Reading can be rewound to the beginning using method Rewind.
//
// Tracks are read to the fEventRec (contains reconstructed tracks) 
// and fEventSim (corresponding MC simulated data) data members,
// that are of the type AliAOD. 
//
// If a given reader has ability of reading both, reconstructed and simulated data, 
// these are structured in AODs so a "n'th" simulated particle 
// (the one stored in the fEventSim at slot n) 
// corresponds to the n'th reconstructed track (the one stored in the fEventRec at slot n).
//
// The same reconstructed track can be present more than ones in the AOD,
// but with a different PID. In this case
// pointer to the corresponding MC simulated particles is also present more than ones.
// This situation happens if you want to read all particles 
// with PID probability of being , e.g.,  pion higher than 60%
// and being kaon higher than 40%. Than, if a given track has probability Ppid(pi)=52% and Ppid(K)=48% 
// than it is read twise.
//
// Provides functionality for both buffering and non-buffering reading
// This can be switched on/off via method SetEventBuffering(bool)
// The main method that inheriting classes need to implement is ReadNext()
// that read next event in queue.
//
// The others are:
// Bool_t  ReadsRec() const; specifies if reader is able to read simulated particles
// Bool_t  ReadsSim() const; specifies if reader is able to read reconstructed tracks
// void    Rewind(); rewind reading to the beginning
//
// This class provides full functionality for reading from many sources
// User can provide TObjArray of TObjStrings (SetDirs method or via parameter 
// in the constructor) which desribes paths of directories to search data in.
// If none specified current directory is searched.
// 
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TGliteXmlEventlist.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRandom.h>
#include <TString.h>

#include "AliAOD.h"
#include "AliAODParticleCut.h"
#include "AliAODRun.h"
#include "AliLog.h"
#include "AliReader.h"
 
ClassImp(AliReader)
//pure virtual
    
/*************************************************************************************/

AliReader::AliReader():
 fEventList(0x0),
 fCuts(new TObjArray()),
 fDirs(0x0),
 fCurrentEvent(0),
 fCurrentDir(0),
 fNEventsRead(0),
 fEventRec(0x0),
 fEventSim(0x0),
 fRunSim(0x0),
 fRunRec(0x0),
 fIsRead(kFALSE),
 fBufferEvents(kFALSE),
 fBlend(kFALSE),
 fFirst(0),
 fLast(0),
 fTrackCounter(0x0)
{
//constructor
}
/*************************************************************************************/

AliReader::AliReader(TObjArray* dirs):
 fEventList(0x0),
 fCuts(new TObjArray()),
 fDirs(dirs),
 fCurrentEvent(0),
 fCurrentDir(0),
 fNEventsRead(0),
 fEventRec(0x0),
 fEventSim(0x0),
 fRunSim(0x0),
 fRunRec(0x0),
 fIsRead(kFALSE),
 fBufferEvents(kFALSE),
 fBlend(kFALSE),
 fFirst(0),
 fLast(0),
 fTrackCounter(0x0)
{
//ctor with array of directories to read as parameter
}
/*************************************************************************************/
AliReader::AliReader(const AliReader& in):
 TNamed(in),
 fEventList((in.fEventList)?(TGliteXmlEventlist*)in.fEventList->Clone():0x0),
 fCuts((in.fCuts)?(TObjArray*)in.fCuts->Clone():0x0),
 fDirs((in.fDirs)?(TObjArray*)in.fDirs->Clone():0x0),
 fCurrentEvent(0),
 fCurrentDir(0),
 fNEventsRead(0),
 fEventRec(0x0),
 fEventSim(0x0),
 fRunSim(0x0),
 fRunRec(0x0),
 fIsRead(kFALSE),
 fBufferEvents(in.fBufferEvents),
 fBlend(in.fBlend),
 fFirst(in.fFirst),
 fLast(in.fLast),
 fTrackCounter(0x0)
{
 //cpy constructor
}

AliReader::~AliReader()
{
//destructor
 if(fCuts)
  {
   fCuts->SetOwner();
   delete fCuts;
  }
 delete fEventSim;
 delete fEventRec;
 delete fTrackCounter;
 delete fEventList;
}
/*************************************************************************************/

AliReader& AliReader::operator=(const AliReader& in)
{
  //Assigment operator
 if (this == &in) return *this;  
 TNamed::operator=( (const TNamed&)in );
  
 fCuts = (in.fCuts)?(TObjArray*)in.fCuts->Clone():0x0;
 fDirs = (in.fDirs)?(TObjArray*)in.fDirs->Clone():0x0;
 fCurrentEvent = 0;
 fCurrentDir = 0;
 fNEventsRead = 0;
 fEventRec = 0x0;
 fEventSim = 0x0;
 fRunSim = 0x0;
 fRunRec = 0x0;
 fIsRead = kFALSE;
 fBufferEvents = in.fBufferEvents;
 fBlend = in.fBlend;
 fFirst = in.fFirst;
 fLast = in.fLast;
 fTrackCounter = 0x0;
 return *this;  
}
/*************************************************************************************/

Int_t AliReader::Next()
{
//moves to next event

  //if asked to read up to event nb. fLast, and it is overcome, report no more events
  if ((fNEventsRead >= fLast) && (fLast > 0) ) return kTRUE;
  
  if (fTrackCounter == 0x0)//create Track Counter
   {
     fTrackCounter = new TH1I("trackcounter","Track Counter",20000,0,20000);
     fTrackCounter->SetDirectory(0x0);
   }
  
  do //if asked to read from event fFirst, rewind to it
   {
    if ( ReadNext() == kTRUE) //if no more evets, return it
      return kTRUE;
   }while (fNEventsRead < fFirst);
   
  //here we have event
  
  if (fBlend) Blend();//Mix particles order 
  
  if (fBufferEvents)//store events if buffering is on
   {
     if ( ReadsRec() && fEventRec) 
       fRunRec->SetEvent(fNEventsRead-1-fFirst,fEventRec);
     if ( ReadsSim() && fEventSim)
       fRunSim->SetEvent(fNEventsRead-1-fFirst,fEventSim);
   }
  return kFALSE;
}
/*************************************************************************************/

void AliReader::AddParticleCut(AliAODParticleCut* cut)
{
 //sets the new cut. MAKES A COPY OF THE CUT !!!!
 
  if (!cut) //if cut is NULL return with error
   {
    Error("AddParticleType","NULL pointers are not accepted any more.\nIf You want to accept all particles of this type, set an empty cut ");
    return;
   }
  AliAODParticleCut *c = (AliAODParticleCut*)cut->Clone();
  fCuts->Add(c);
}
/********************************************************************/

AliAOD* AliReader::GetEventSim(Int_t n)
 {
 //returns Nth event with simulated particles
  if (ReadsSim() == kFALSE)
   {
     Error("GetParticleEvent","This reader is not able to provide simulated particles.");
     return 0;
   } 
   
  if (!fIsRead) 
   { 
    if (ReadsSim() && (fRunSim == 0x0)) fRunSim = new AliAODRun();
    if (ReadsRec() && (fRunRec == 0x0)) fRunRec = new AliAODRun();
    
    if (Read(fRunSim,fRunRec))
     {
       Error("GetParticleEvent","Error in reading");
       return 0x0;
     }
    else fIsRead = kTRUE;
   }
  return fRunSim->GetEvent(n);
 }
/********************************************************************/

AliAOD* AliReader::GetEventRec(Int_t n)
 {
 //returns Nth event with reconstructed tracks
  if (ReadsRec() == kFALSE)
   {
     Error("GetTrackEvent","This reader is not able to provide recosntructed tracks.");
     return 0;
   } 
  if (!fIsRead) 
   {
    if (ReadsSim() && (fRunSim == 0x0)) fRunSim = new AliAODRun();
    if (ReadsRec() && (fRunRec == 0x0)) fRunRec = new AliAODRun();
    
    if(Read(fRunSim,fRunRec))
     {
       Error("GetTrackEvent","Error in reading");
       return 0x0;
     }
    else fIsRead = kTRUE;
   }
  return fRunRec->GetEvent(n);
 }
/********************************************************************/

Int_t AliReader::GetNumberOfSimEvents()
 {
 //returns number of events of particles
  if (ReadsSim() == kFALSE)
   {
     Error("GetNumberOfPartEvents","This reader is not able to provide simulated particles.");
     return 0;
   } 
   
  if (!fIsRead) 
   {
    if (ReadsSim() && (fRunSim == 0x0)) fRunSim = new AliAODRun();
    if (ReadsRec() && (fRunRec == 0x0)) fRunRec = new AliAODRun();
    
    if (Read(fRunSim,fRunRec))
     {
       Error("GetNumberOfPartEvents","Error in reading");
       return 0;
     }
    else fIsRead = kTRUE;
   }
   return fRunSim->GetNumberOfEvents();
 }
/********************************************************************/
 
Int_t AliReader::GetNumberOfRecEvents()
 {
 //returns number of events of tracks
  if (ReadsRec() == kFALSE)
   {
     Error("GetNumberOfTrackEvents","This reader is not able to provide recosntructed tracks.");
     return 0;
   } 
  if (!fIsRead)
   {
     if (ReadsSim() && (fRunSim == 0x0)) fRunSim = new AliAODRun();
     if (ReadsRec() && (fRunRec == 0x0)) fRunRec = new AliAODRun();
     
     if(Read(fRunSim,fRunRec))
      {
        Error("GetNumberOfTrackEvents","Error in reading");
        return 0;
      }
     else fIsRead = kTRUE;
   }
  return fRunRec->GetNumberOfEvents();
 }
/********************************************************************/

Int_t AliReader::Read(AliAODRun* particles, AliAODRun *tracks)
{
 //reads data and puts put to the particles and tracks objects
 //reurns 0 if everything is OK
 //
  Info("Read","");
  
  if ( ReadsSim() && (particles == 0x0) ) //check if an object is instatiated
   {
     Error("Read"," particles object must be instatiated before passing it to the reader");
     return 1;
   }
  if ( ReadsRec() && (tracks == 0x0) )  //check if an object is instatiated
   {
     Error("Read"," tracks object must be instatiated before passing it to the reader");
     return 1;
   }
   
  if (ReadsSim()) particles->Reset();//clear runs == delete all old events
  if (ReadsRec()) tracks->Reset();
  
  Rewind();
  
  Int_t i = 0;
  while(Next() == kFALSE)
   {
     if (ReadsRec()) tracks->SetEvent(i,fEventRec);
     if (ReadsSim()) particles->SetEvent(i,fEventSim);
     i++;
   }
  return 0;
}      
/*************************************************************************************/

Bool_t AliReader::Rejected(AliVAODParticle* p)
{
 //Method examines whether particle meets all cut and particle type criteria
  
   if(p==0x0)//of corse we not pass NULL pointers
    {
     Warning("Rejected()","No Pasaran! We never accept NULL pointers");
     return kTRUE;
    }
   //if no particle is specified, we pass all particles
   //excluding NULL pointers, of course
  if ( fCuts->GetEntriesFast() == 0 ) return kFALSE; //if no cut specified accept all particles
  for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliAODParticleCut &cut = *((AliAODParticleCut*)fCuts->At(i));
     if(!cut.Rejected(p)) return kFALSE;  //accepted
   }
   
  return kTRUE;//not accepted
}
/*************************************************************************************/

Bool_t  AliReader::Rejected(Int_t pid)
{
//this method checks if any of existing cuts accepts this pid particles
//or any cuts accepts all particles

 if(pid == 0)
  return kTRUE;

 if ( fCuts->GetEntriesFast() == 0 ) return kFALSE; //if no cut specified accept all particles
  
 for(Int_t i=0; i<fCuts->GetEntriesFast(); i++)   
   {
     AliAODParticleCut &cut = *((AliAODParticleCut*)fCuts->At(i));
     //if some of cuts accepts all particles or some accepts particles of this type, accept
     if ( (cut.GetPID() == 0) || (cut.GetPID() == pid) ) return kFALSE; 
   }
 return kTRUE;
}
/*************************************************************************************/

TString AliReader::GetDirName(Int_t entry)
{
//returns directory name of next one to read
  TString  retval;//return value
  if (fDirs ==  0x0)
   { 
     if (entry == 0)
      {
       retval = ".";
       return retval;
      }
     else
      {
       return retval;
      }  
   }
  
  
  if ( (entry >= fDirs->GetEntries()) || (entry < 0))//if out of bounds return empty string
   {                                            //note that entry==0 is accepted even if array is empty (size=0)
    if ( (fDirs->GetEntries() == 0) && (entry == 0) )
      { 
        retval = ".";
        return retval;
      }
     AliDebug(1,Form("Index %d out of bounds",entry));

     return retval;
   }


  TClass *objclass = fDirs->At(entry)->IsA();
  TClass *stringclass = TObjString::Class();

  TObjString *dir = (TObjString*)objclass->DynamicCast(stringclass,fDirs->At(entry));

  if(dir == 0x0)
   {
     Error("GetDirName","Object in TObjArray is not a TObjString or its descendant");
     return retval;
   }
  AliDebug(1,Form("Returned ok %s",dir->String().Data()));
  retval = dir->String();
  return retval;
}
/*************************************************************************************/

void AliReader::Blend()
{
  //randomly change positions of the particles after reading
  //is used to check if some distributions (of many particle properties) 
  //depend on the order of particles
  //(tracking gives particles Pt sorted)
  Int_t npart = 0;

  if (fEventSim ) 
   {
     npart = fEventSim->GetNumberOfParticles();
    } 
  else 
    if (fEventRec ) 
     {
        npart = fEventRec->GetNumberOfParticles();
     }
    else
     {
       return;
     }
  for (Int_t i = 2;  i < npart; i++)
   {
     Int_t with = gRandom->Integer(i);
//     Info("Blend","%d %d",i, with);
     if (fEventSim) fEventSim->SwapParticles(i,with);
     if (fEventRec) fEventRec->SwapParticles(i,with);
   }
}
/*************************************************************************************/

void AliReader::WriteTrackCounter() const
{
 //writes the counter histogram
 
 if (fTrackCounter) fTrackCounter->Write(0,TObject::kOverwrite);
 else 
  {
    Warning("WriteTrackCounter","Counter is NULL");
  }
}
