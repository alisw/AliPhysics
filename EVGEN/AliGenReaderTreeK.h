#ifndef ALIGENREADERTreeK_H
#define ALIGENREADERTREEK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//
// Realisation of AliGenReader to be used with AliGenExtFile
// It reads events from a kinematics TreeK.
// Author: andreas.morsch@cern.ch
//
#include "AliGenReader.h"

class TFile;
class AliRunLoader;
class AliStack;
class TString;
class TObjArray;

class AliGenReaderTreeK : public AliGenReader
{
 public:
    AliGenReaderTreeK();
    AliGenReaderTreeK(const AliGenReaderTreeK &reader);
    virtual ~AliGenReaderTreeK();
    // Initialise 
    virtual void Init();
    // Read
    virtual Int_t NextEvent();
    virtual TParticle*  NextParticle();
    virtual void RewindEvent();
    virtual void SetOnlyPrimaries(Bool_t flag){fOnlyPrimaries = flag;}
    AliGenReaderTreeK & operator=(const AliGenReaderTreeK & rhs);
    void SetDirs(TObjArray* dirs){fDirs = dirs;} //sets array directories names
    void AddDir(const char* dirname);
    AliRunLoader * GetRunLoader() const {return fInRunLoader;}

 protected:
    Int_t             fNcurrent;          // points to the next entry
    Int_t             fNparticle;         // Next particle in list
    Int_t             fNp;                // number of particles
    AliRunLoader     *fInRunLoader;       //!Run Loader of the input event
    TFile            *fBaseFile;          //! pointer to base file
    AliStack         *fStack;             //! Particle stack
    Bool_t            fOnlyPrimaries;     // Flag indicating wether only primaries are read
    TObjArray        *fDirs;              //arry with directories to read data from
    Int_t             fCurrentDir;        //Number of current directory
    static const TString fgkEventFolderName;//!name of folder where event to read is mounted
    
    TString&   GetDirName(Int_t entry);
    TParticle* GetParticle(Int_t i);
 private:
    void Copy(TObject&) const;    
    ClassDef(AliGenReaderTreeK,1) // Read particles from TreeK
};

inline 
TParticle* AliGenReaderTreeK::GetParticle(Int_t i)
 {
  if (fStack && i<fNp) return fStack->Particle(i);
  return 0x0;
 }

#endif






