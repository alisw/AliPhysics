#ifndef ALIJETCORRELMAKER_H
#define ALIJETCORRELMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________
// Class that uses info from the AliJetCorrelSelector object to set up the
// two-particle correlations to be run in one instance of the analysis module.
//-- Author: Paul Constantin

#include "AliJetCorrelSelector.h"

class AliJetCorrelMaker : public TObject {
  
 public:
  AliJetCorrelMaker();
  ~AliJetCorrelMaker();
  
  Bool_t Init(UInt_t s, UInt_t * const v);
  Bool_t Check() const;
  void Show()const;
  
  UInt_t NoOfCorrel() const {return fNumCorrel;}
  UInt_t NoOfTrigg()  const {return fNumTrigg;}
  UInt_t NoOfAssoc()  const {return fNumAssoc;}
  UInt_t IdxTrigg(UInt_t k) const;
  UInt_t IdxAssoc(UInt_t k) const;
  cPartType_t TriggType(UInt_t k) const;
  cPartType_t AssocType(UInt_t k) const;
  TString Descriptor(UInt_t k) const;
  Bool_t RecoTrigger(UInt_t k) const;
  Bool_t RecoTrigger() const;
  
 private: 
  UInt_t fNumCorrel, fNumTrigg, fNumAssoc; //! counters
  UInt_t* fCorrelType;     //! array of correlation types
  TString* fCorrelStr;     //! array of correlation string descriptors
  cPartType_t* fTriggType;  //! array of trigger particle types
  cPartType_t* fAssocType;  //! array of associated particle types
  UInt_t* fIdxTrigg;       //! array with trigger indices
  UInt_t* fIdxAssoc;       //! array with associated indices
  
  // disable (make private) copy constructor, and assignment operator:
  AliJetCorrelMaker(const AliJetCorrelMaker&);
  AliJetCorrelMaker& operator=(const AliJetCorrelMaker&);
  
  ClassDef(AliJetCorrelMaker, 1);
};

#endif 
