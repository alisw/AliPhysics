#ifndef ALIMCQA_H
#define ALIMCQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Quality assurance services for MC                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
class TObjArray;
class TBrowser;

class AliMCQA : public TObject
{
public:
  AliMCQA() {}
  AliMCQA(Int_t ndets);
  virtual ~AliMCQA() {delete fQAList;fQAList=0;}
  Bool_t  IsFolder() const {return kTRUE;}
  virtual  void  Browse(TBrowser *b);
  virtual  void  PreTrack();
  virtual  TObjArray *GetQAList() const {return fQAList;}
  
  // QA step manager
  virtual void StepManager(Int_t id);

protected:
  Int_t      fNdets;       // Number of detectors
  TObjArray *fQAList;      // QA histograms
  Int_t      fOldId;       //! ID of the current module
  Int_t     *fDetDone;     //! Detector done flag 

private:
  AliMCQA(const AliMCQA &) {}
  AliMCQA & operator=(const AliMCQA &) {return (*this);}

  ClassDef(AliMCQA,1)  //Quality Assurance class for the MC
};

#endif 

