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
//  Author Philippe Crochet (LPCCFd)

#include <TObject.h>
#include <TArrayF.h>

class AliMpLocalBoard;
class AliMUONGeometryTransformer;
class AliMpPad;
class AliMpVSegmentation;

class AliMUONTriggerCircuit : public TObject 
{
public: 
  AliMUONTriggerCircuit(const AliMUONGeometryTransformer* transformer);  
  virtual ~AliMUONTriggerCircuit();
     // copy constructor
  AliMUONTriggerCircuit(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 
  // assignment operator
  AliMUONTriggerCircuit& operator=(const AliMUONTriggerCircuit& AliMUONTriggerCircuit); 

  //--- methods which return member data related info
  Float_t GetX11Pos(Int_t localBoardId, Int_t istrip) const;
  Float_t GetY11Pos(Int_t localBoardId, Int_t istrip) const;
  Float_t GetY21Pos(Int_t localBoardId, Int_t istrip) const;
  Float_t GetZ11Pos(Int_t localBoardId, Int_t istrip) const;
  Float_t GetZ21Pos(Int_t localBoardId, Int_t istrip) const;
  Float_t GetX11Width(Int_t localBoardId, Int_t istrip) const;
  Float_t GetY11Width(Int_t localBoardId, Int_t istrip) const;
  Float_t GetY21Width(Int_t localBoardId, Int_t istrip) const;

  //  void Print(Option_t* opt="") const;
  //  void dump(const char* what, const Float_t* array, Int_t size);
  //  void dump(const char* what, const Int_t* array, Int_t size);
  
  /// Set pointer to transformations
  void  SetTransformer(const AliMUONGeometryTransformer* transformer) {fkTransformer = transformer;}
  /// Get pointer to transformations
  const AliMUONGeometryTransformer* GetTransformer() const {return fkTransformer;}
  Float_t PtCal(Int_t localBoardId, Int_t istripX, Int_t idev, Int_t istripY) const;
  
private:

  void LoadYPos(AliMpLocalBoard* const localBoard);
  void LoadXPos(AliMpLocalBoard* const localBoard);

  Int_t FirstStrip(AliMpLocalBoard* localBoard);

  void FillXstrips(const Int_t icol, 
                   const Int_t iFirstStrip, const Int_t iLastStrip,
                   Int_t liStripCircuit, const Bool_t is11);
  
  void FillYstrips(const Int_t iFirstStrip,
                   const Int_t iLastStrip, Int_t liStripCircuit,
                   const Bool_t doubling);

  void XYGlobal(const AliMpPad& pad,
                Double_t* xyGlobal);    
  

private:    
  TArrayF fXpos11[235];         ///< X position of Y strips in MC11
  TArrayF fYpos11[235];         ///< Y position of X strips in MC11
  TArrayF fYpos21[235];         ///< Y position of X strips in MC21
  TArrayF fZpos11[235];         ///< Z position of X strips in MC11
  TArrayF fZpos21[235];         ///< Z position of X strips in MC21
  TArrayF fXwidth11[235];         ///< width of Y strips in MC11
  TArrayF fYwidth11[235];         ///< width of X strips in MC11
  TArrayF fYwidth21[235];         ///< width of X strips in MC21

  const AliMUONGeometryTransformer* fkTransformer; //!<! pointer to transformation
  const AliMpVSegmentation* fkCurrentSeg;          //!<! current segmentation

  Int_t fCurrentDetElem;                          //!<! current detection elt id
  Int_t fCurrentLocalBoard;                       //!<! current local board id

  ClassDef(AliMUONTriggerCircuit,2) // Trigger Circuit class
};
#endif
