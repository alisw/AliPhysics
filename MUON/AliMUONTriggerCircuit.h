#ifndef ALIMUONTRIGGERCIRCUIT_H
#define ALIMUONTRIGGERCIRCUIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONTriggerCircuit
/// \brief MUON Trigger circuit
///
/// \author Philippe Crochet (LPCCFd)

#include <TObject.h>
#include <TObjArray.h>
#include "AliMpPad.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMpVSegmentation.h"

class TObjArray;
class AliMUONTriggerCrateStore;
class AliMUONGeometryTransformer;

class AliMUONTriggerCircuit : public TObject 
{
public: 
  AliMUONTriggerCircuit();  
  virtual ~AliMUONTriggerCircuit();
     // copy constructor
  AliMUONTriggerCircuit(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 
  // assignment operator
  AliMUONTriggerCircuit& operator=(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 

  // initializations
  void Init(Int_t iCircuit, const AliMUONTriggerCrateStore& crates);    
  
  //--- methods which return member data related info
  Float_t GetY11Pos(Int_t istrip) const;
  Float_t GetY21Pos(Int_t istrip) const;
  Float_t GetX11Pos(Int_t istrip) const;
  Int_t   DetElemId(Int_t ichamber, char side, Int_t iline);
  void DecodeBoardName(const char* boardName, char& side,
                       Int_t& iLine,
                       Int_t& iCol);
  Int_t DetElemId(Int_t ichamber, const char* localBoardName);
  
  //  void Print(Option_t* opt="") const;
  //  void dump(const char* what, const Float_t* array, Int_t size);
  //  void dump(const char* what, const Int_t* array, Int_t size);
  
  void  SetTransformer(const AliMUONGeometryTransformer* transformer) {fTransformer = transformer;}

private:

  void LoadYPos(const AliMUONTriggerCrateStore& crates);
  void LoadXPos(const AliMUONTriggerCrateStore& crates);
  Int_t FirstStrip(const char* localBoardName);

  void FillXstrips(const AliMpVSegmentation* seg,
		   const Int_t detElemId, const Int_t icol, 
                   const Int_t iFirstStrip, const Int_t iLastStrip,
                   Int_t liStripCircuit, Float_t *tab);
  
  void FillYstrips(const AliMpVSegmentation* seg,
		   const Int_t detElemId, const Int_t iFirstStrip,
                   const Int_t iLastStrip, Int_t liStripCircuit,
                   const Bool_t doubling);

  void XYGlobal(Int_t detElemId, const AliMpPad& pad,
                Double_t xyGlobal[4]);    
  

private:    
  Int_t fILocalBoard;          ///< local board number
  Float_t fXpos11[16];         ///< X position of Y strips in MC11
  Float_t fYpos11[31];         ///< Y position of X strips in MC11
  Float_t fYpos21[63];         ///< Y position of X strips in MC21

  const AliMUONGeometryTransformer* fTransformer; //!< pointer to transformation

  ClassDef(AliMUONTriggerCircuit,1) // Trigger Circuit class
};
#endif
