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
class TPaveLabel;

class AliMCQA : public TObject
{
public:
  AliMCQA();
  AliMCQA(Int_t ndets);
  virtual ~AliMCQA();
  Bool_t  IsFolder() const {return kTRUE;}
  virtual  void  Browse(TBrowser *b);
  virtual  void  PreTrack();
  virtual  TObjArray *GetQAList() const {return fQAList;}
  void DrawModuleName();
  void AddModuleName();
  void DrawVolumeName();
  void AddVolumeName();


  // QA step manager
  virtual void StepManager(Int_t id);

protected:
  Int_t       fNdets;       // Number of detectors
  Int_t       fNvolumes;    // Number of volumes
  TObjArray  *fQAList;      // QA lists of histograms
  Int_t       fOldId;       //! ID of the current module
  Int_t      *fDetDone;     //! Detector done flag 
  TObjArray  *fQAHist;      // Global QA histograms
  TObjArray  *fVolNames;    // Volume names
  TObjArray  *fModNames;    // Module names
  TPaveLabel *fMPaveLabel;  //! PaveLabel for the Modules
  TPaveLabel *fVPaveLabel;  //! PaveLabel for the Volumes

private:
  AliMCQA(const AliMCQA &) {}
  AliMCQA & operator=(const AliMCQA &) {return (*this);}
  void DrawPaveLabel(TPaveLabel *&pv);
  Int_t GetHBin(const char* hname);

  ClassDef(AliMCQA,1)  //Quality Assurance class for the MC
};

#endif 

