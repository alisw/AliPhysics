#ifndef AliMFTHalfDiskSegmentation_H
#define AliMFTHalfDiskSegmentation_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class for the description of the virtual structure for an half-disk of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TClonesArray.h"
#include "AliLog.h"
#include "AliMFTLadderSegmentation.h"
#include "TXMLEngine.h"
#include "AliMFTVSegmentation.h"


//====================================================================================================================================================

class AliMFTHalfDiskSegmentation : public AliMFTVSegmentation {

public:

  enum {kFront, kBack, kBoth, kProfile};

  AliMFTHalfDiskSegmentation();
  AliMFTHalfDiskSegmentation(UInt_t uniqueID);
  AliMFTHalfDiskSegmentation(const AliMFTHalfDiskSegmentation& pt);
  AliMFTHalfDiskSegmentation& operator=(const AliMFTHalfDiskSegmentation &source);
  
  virtual ~AliMFTHalfDiskSegmentation();

  virtual void Clear(const Option_t* /*opt*/);
  
  virtual void Print(Option_t* opt="");
  
  void CreateLadders(TXMLEngine* xml, XMLNodePointer_t node);
  
  Int_t    GetNLaddersBuild()  const {return fLadders->GetEntriesFast();};
  Int_t    GetNLadders()  const {return fNLadders;};
  void    SetNLadders(Int_t val)   {fNLadders = val;};

  AliMFTLadderSegmentation* GetLadder(Int_t iLadder) { return (AliMFTLadderSegmentation*) fLadders->At(iLadder); }
  
  Double_t GetZ() const {const Double_t *pos = GetTransformation()->GetTranslation(); return pos[2];};


  Int_t GetNChips();
  
private:
  
  Int_t fNLadders;

  TClonesArray *fLadders;
  
  ClassDef(AliMFTHalfDiskSegmentation, 1)

};

//====================================================================================================================================================
	
#endif

