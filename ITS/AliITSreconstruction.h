#ifndef ALIITSRECONSTRUCTION_H
#define ALIITSRECONSTRUCTION_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
 
/*
  $Id$
 */

#include <TTask.h>

class AliRunLoader;
class TString;
class AliITS;

class AliITSreconstruction : public TTask{
 public:
    AliITSreconstruction(); // default constructor
    AliITSreconstruction(const char *filename); // standard constructor
    AliITSreconstruction(AliRunLoader *rl); // standard constructor
    virtual ~AliITSreconstruction();//Destructor
    virtual Bool_t Init();
    virtual void Exec(const Option_t *opt="ALL");
    virtual void SetOutputFile(TString filename);
 private:
    AliITSreconstruction(const AliITSreconstruction &);
    AliITSreconstruction & operator=(const AliITSreconstruction &);
    Bool_t InitRec();  // Standard Reconstrution initilization.
 private:
    Bool_t  fDet[3];   //! logical specifing which detectors to reconstruct.
    Bool_t  fInit;     //! True if Init was sucessfull, else false.
    Int_t   fEnt;      //! Number of events to processevent index.
    Int_t   fEnt0;     //! first event to process, default 0.
    AliITS  *fITS;     //! Local pointer to ITS class.
    Bool_t  fDfArp;    //! if True then delete fRunLoader in destructor.

    AliITSLoader *fLoader; //! ITS loader
    AliRunLoader* fRunLoader;//!Run Loader

    ClassDef(AliITSreconstruction,2) // Task to Reconstruct ITS from Digits.

};
#endif
