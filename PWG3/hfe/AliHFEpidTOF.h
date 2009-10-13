#ifndef ALIHFEPIDTOF_H
#define ALIHFEPIDTOF_H

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
 
#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class TList;
class TH2F;

class AliAODTrack;
class AliAODMCParticle;
class AliESDtrack;
class AliMCParticle;

class AliHFEpidTOF : public AliHFEpidBase{
  typedef enum{
    kHistTOFpidFlags = 0,
      kHistTOFpidBetavP = 1,
      kHistTOFsignal = 2,
      kHistTOFlength =3,
      kHistTOFpid0 = 4,
      kHistTOFpid1 = 5,
      kHistTOFpid2 = 6,
      kHistTOFpid3 = 7,
      kHistTOFpid4 = 8
      
      } QAHist_t;
 public:
  AliHFEpidTOF(const Char_t *name);
  virtual ~AliHFEpidTOF();
  AliHFEpidTOF(const AliHFEpidTOF &c);
  AliHFEpidTOF &operator=(const AliHFEpidTOF &c);
  
  virtual Bool_t    InitializePID();
  virtual Int_t     IsSelected(AliHFEpidObject *track);
  virtual Bool_t    HasQAhistos() const { return kTRUE; };
  
  
 protected:
  void Copy(TObject &ref) const;
  void AddQAhistograms(TList *qaHist);
  Int_t MakePIDesd(AliESDtrack *esdTrack, AliMCParticle *mcTrack);
  Int_t MakePIDaod(AliAODTrack *aodTrack, AliAODMCParticle *mcTrack);
  
 private:
  
  AliPID *fPID;           //! PID Object
  TList *fQAList;         //! QA histograms
  ClassDef(AliHFEpidTOF, 1)
};

#endif
