#ifndef ALIJETPARTICLESREADER_H
#define ALIJETPARTICLESREADER_H

/* $Id$ */

//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliJetReader
//
// loizides@ikf.uni-frankfurt.de
///////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TObjArray.h>
class TClonesArray;
class TString;
class AliJetEventParticles;

class AliJetParticlesReader: public TNamed
{
  public:
  AliJetParticlesReader();
  AliJetParticlesReader(TObjArray* dirs);
  virtual ~AliJetParticlesReader();

  void SetDirs(TObjArray* dirs){fDirs = dirs;} //sets array directories names
  void ReadEventsFromTo(Int_t first,Int_t last){fFirst = first; fLast = last;}

  virtual Int_t Next(); //call this if you want the next event
  virtual void  Rewind() = 0;

  void SetPtCut(Float_t ptmin=0, Float_t ptmax=1000)
    {fPtMin=ptmin;fPtMax=ptmax;}
  void SetPhiCut(Float_t phi=2*TMath::Pi()){SetPhiCut(0,phi);}
  void SetPhiCut(Float_t phimin, Float_t phimax)
    {fPhiMin=phimin;fPhiMax=phimax;}
  void SetEtaCut(Float_t e=1){SetEtaCut(-e,e);}
  void SetEtaCut(Float_t emin, Float_t emax)
    {fEtaMin=emin;fEtaMax=emax;}

  const AliJetEventParticles* GetEventParticles() const {return fEventParticles;}
  AliJetEventParticles* GetEventParticles(Bool_t o) 
    {fOwner=o; return fEventParticles;} //return particles and set ownership

  Int_t GetNumberOfDirs()   const {return (fDirs)?fDirs->GetEntries():0;}
  Int_t GetCurEventNumber() const {return fCurrentEvent;}
  Int_t GetCurDirNumber()   const {return fCurrentDir;}
  Int_t GetTotEventsRead()  const {return fNEventsRead;} 
  Int_t GetFirstEvent()     const {return fFirst;}
  Int_t GetLastEvent()      const {return fLast;}

  protected:

  virtual Int_t ReadNext() = 0; //this methods reads next event and 
                                 //put result in fParticles
  TString& GetDirName(Int_t entry);

  AliJetEventParticles* fEventParticles; //array with read particles
  Bool_t fOwner;              //ownership of particles
  TObjArray*   fDirs;         //array with directories to read data from
  Int_t        fCurrentEvent; //number of current event in current file
  Int_t        fCurrentDir;   //number of current directory in array
  Int_t        fNEventsRead;  //total number of processed events 
  Int_t        fFirst;        //first event to return (all before are skipped)
  Int_t        fLast;         //last event to return (relative to total number)

  Float_t fPtMin;   //min pt cut
  Float_t fPtMax;   //max pt cut
  Float_t fEtaMin;  //min eta cut
  Float_t fEtaMax;  //max eta cut
  Float_t fPhiMin;  //min phi cut
  Float_t fPhiMax;  //max phi cut

  ClassDef(AliJetParticlesReader,1) // Basic AliJetParticles Reader class
};

#endif
