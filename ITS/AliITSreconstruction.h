#ifndef ALIITSRECONSTRUCTION_H
#define ALIITSRECONSTRUCTION_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
 
/*
  $Id$
 */
/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Class for ITS RecPoint reconstruction                               //
//                                                                     //
////////////////////////////////////////////////////////////////////////

#include <TTask.h>

class AliRunLoader;
class TString;

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
    Bool_t InitRec();  // Standard Reconstrution initilization.
 private:

    AliITSreconstruction(const AliITSreconstruction& rec);
    AliITSreconstruction& operator=(const AliITSreconstruction &source);

    Bool_t  fDet[3];   //! logical specifing which detectors to reconstruct.
    Bool_t  fInit;     //! True if Init was sucessfull, else false.
    Int_t   fEnt;      //! Number of events to processevent index.
    Int_t   fEnt0;     //! first event to process, default 0.
    AliITSDetTypeRec *fDetTypeRec; //!ITS obj. for reconstruction
    Bool_t  fDfArp;    //! if True then delete fRunLoader in destructor.
    AliITSgeom*   fITSgeom;//! ITS geometry
    AliITSLoader *fLoader; //! ITS loader
    AliRunLoader* fRunLoader;//!Run Loader
 
    ClassDef(AliITSreconstruction,3) // Task to Reconstruct ITS from Digits.

};
#endif
