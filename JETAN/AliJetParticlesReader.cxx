// $Id$

//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//                                                                       //
// class AliJetParticlesReader                                           //
//                                                                       //
// This reader reads tracks from Event Summary Data                      //
// taken from Piotr.Skowronski@cern.ch                                   //
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html  //
//                                                                       //
// loizides@ikf.uni-frankfurt.de                                         //
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TObjString.h>
#include <TTree.h>
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
    fPhiMin(0),fPhiMax(2*TMath::Pi()),
    fNewTree(0),
    fTree(0)
{
  //Constructor
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
    fPhiMin(0),fPhiMax(2*TMath::Pi()),
    fNewTree(0),
    fTree(0)
{
}

AliJetParticlesReader::~AliJetParticlesReader()
{
  //Constructor
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
  if(fTree && fEventParticles){
      if(fNewTree){
	fTree->Branch("particles","AliJetEventParticles",&fEventParticles,32000,99);
	fNewTree=0;
      }
      fEventParticles->SetEventNr(fCurrentDir*1000+fCurrentEvent);
      fTree->Fill();
  }
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

