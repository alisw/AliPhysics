#ifndef ALIMUONREJECTLIST_H
#define ALIMUONREJECTLIST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONRejectList
/// \brief Object to hold the list of elements we want to reject
/// from the reconstruction
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVStore;

class AliMUONRejectList : public TObject
{
public:
  AliMUONRejectList(TRootIOCtor* ioCtor);
  AliMUONRejectList();
  AliMUONRejectList(const AliMUONRejectList& rl);
  AliMUONRejectList& operator=(const AliMUONRejectList& rl);
  virtual ~AliMUONRejectList();

  /// Clone this object
  virtual TObject* Clone(const char* /*name*/="") const { return new AliMUONRejectList(*this); }
  
  /// Whether or not all our probabilities are 0.0 or 1.0
  Bool_t IsBinary() const { return fIsBinary; }
  
  Float_t DetectionElementProbability(Int_t detElemId) const;
  Float_t BusPatchProbability(Int_t busPatchId) const;
  Float_t ManuProbability(Int_t detElemId, Int_t manuId) const;
  Float_t ChannelProbability(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  void SetDetectionElementProbability(Int_t detElemId, Float_t proba=1.0);
  void SetBusPatchProbability(Int_t busPatchId, Float_t proba=1.0);
  void SetManuProbability(Int_t detElemId, Int_t manuId, Float_t proba=1.0);
  void SetChannelProbability(Int_t detElemId, Int_t manuId, Int_t manuChannel, Float_t proba=1.0);
  
  void SetPCBProbability(Int_t detElemId, Int_t pcbNumber, Float_t proba=1.0);
  void SetHVProbability(const char* dcsName, Float_t proba=1.0);
  
  void Print(Option_t* opt="") const;
  
private:
  void ZeroOrOne(Float_t proba);
  
private:
  
  Bool_t fIsBinary; ///< Whether or not we only store zeros and ones for probabilities  

  UInt_t fMaxNofDEs; ///< max number of detection elements (for allocation purposes)  
  UInt_t fMaxNofBPs; ///< max number of bus patches (for allocation purposes)
  UInt_t fMaxNofManus; ///< max number of manus (for allocation purposes)
  
  UInt_t fNofDEs; ///< actual number of detection elements for which we have probabilities
  UInt_t fNofBPs; ///< actual number of bus patches for which we have probabilities
  UInt_t fNofManus; ///< actual number of manus for which we have probabilities

  /// array of detection element ids
  UInt_t* fDEIds; //[fMaxNofDEs] 

  /// array of probabilities of DEs
  Float_t* fDEProbas; //[fMaxNofDEs] 

  /// array of bus patch ids
  UInt_t* fBPIds; //[fMaxNofBPs] 
  
  /// array of proba for bus patches
  Float_t* fBPProbas; //[fMaxNofBPs]

  /// array of manu ids
  UInt_t* fManuIds; //[fMaxNofManus]
  
  /// array of proba for manus
  Float_t* fManuProbas; //[fMaxNofManus]

  AliMUONVStore* fChannels; ///< probabilities for all channels
  
  ClassDef(AliMUONRejectList,1) // (probabilistic) Reject list for MUON Tracker
};

#endif
