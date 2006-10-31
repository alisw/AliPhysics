#ifndef ALIITSSDIGITIZE_H
#define ALIITSSDIGITIZE_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */
 
/*
  $Id$
 */
/////////////////////////////////////////////////////////
//                                                     //
//                                                     //
//                                                     //
///////////////////////////////////////////////////////// 
 
#include <TTask.h>
 
class TString;
class AliITS;
 
class AliITSsDigitize : public TTask{
 public:
    AliITSsDigitize(); // default constructor
    AliITSsDigitize(const char *filename); // standard constructor
    virtual ~AliITSsDigitize();//Destructor
    virtual Bool_t Init();
    virtual void Exec(const Option_t *opt="ALL");
 private:
    Bool_t InitSDig() const;  // Standard SDigitization initilization.
 private:
//    TFile   *fFile;    //! pointer to the file contatining the hits and
//                       // and will contain the SDigits
    AliITSsDigitize(const AliITSsDigitize &source); // copy  constructor
    AliITSsDigitize& operator=(const AliITSsDigitize &source); // operator=

    
    Bool_t  fDet[3];   //! logical specifing which detectors to reconstruct.
    Bool_t  fInit;     //! True if Init was sucessfull, else false.
    Int_t   fEnt;      //! Number of events to processevent index.
    Int_t   fEnt0;     //! first event to process, default 0.
    AliITS  *fITS;     //! Local pointer to ITS class.
    AliRunLoader* fRunLoader;//!Run Loader 
    ClassDef(AliITSsDigitize,1) // Task to SDigitize ITS from Hits.
 
};
#endif
