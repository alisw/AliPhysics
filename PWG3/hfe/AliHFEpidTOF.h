#ifndef __ALIHFEPIDTOF_H__
#define __ALIHFEPIDTOF_H__

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice   */   

/************************************************************************
 *                                                                      *
 * Class for TOF PID                                                    *
 * Implements the abstract base class AliHFEpidBase                     *
 * IsInitialized() does the PID decision                                *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel  <M.Fasel@gsi.de>                                     *
 *   Matus Kalisky <matus.kalisky@cern.ch>  (contact)                   *
 ************************************************************************/
 
#ifndef __ALIHFEPIDBASE_H__
#include "AliHFEpidBase.h"
#endif

class TList;
class TH2F;

class AliVParticle;

class AliHFEpidTOF : public AliHFEpidBase{
  typedef enum{
    kHistTOFpidFlags = 0,
      kHistTOFpid_beta_v_P = 1,
      kHistTOFsignal = 2,
      kHistTOFlength =3,
      kHistTOFpid_0 = 4,
      kHistTOFpid_1 = 5,
      kHistTOFpid_2 = 6,
      kHistTOFpid_3 = 7,
      kHistTOFpid_4 = 8
      
      } QAHist_t;
 public:
  AliHFEpidTOF(const Char_t *name);
  virtual ~AliHFEpidTOF();
  AliHFEpidTOF(const AliHFEpidTOF &c);
  AliHFEpidTOF &operator=(const AliHFEpidTOF &c);
  
  virtual Bool_t    InitializePID();
  virtual Int_t     IsSelected(AliVParticle *track);
  virtual Bool_t    HasQAhistos() const { return kTRUE; };
  
  
 protected:
  void Copy(TObject &ref) const;
  void AddQAhistograms(TList *qaHist);
  
 private:
  
  AliPID *fPID;           //! PID Object
  TList *fQAList;         //! QA histograms
  ClassDef(AliHFEpidTOF, 1)
};

#endif
