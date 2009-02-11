#ifndef ALIMUONTRIGGERGUIBOARD_H
#define ALIMUONTRIGGERGUIBOARD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup evaluation
/// \class AliMUONTriggerGUIboard
/// \brief Trigger GUI utility class: single board object
//  Author Bogdan Vulpescu, LPC Clermont-Ferrand

#include "AliMpPad.h"

#include <TString.h>
#include <TObject.h>

class TClonesArray;
class TBox;

class AliMUONTriggerGUIboard : public TObject
{

public:

  AliMUONTriggerGUIboard();
  virtual ~AliMUONTriggerGUIboard();

  /// get the standard name of this board
  Char_t  *GetBoardName()   const { return (Char_t*)(fName->Data()); };
  /// get the name of the crate containing this board
  Char_t  *GetCrateName()   const { return (Char_t*)(fCrateName->Data()); };
  /// get the working status of this board
  UShort_t GetStatus() const { return fStatus; };
  /// get the number of this board
  Int_t    GetNumber() const { return fID; };
  /// get the id of the detector element
  Int_t    GetDetElemId() const { return fDetElemId; };
  /// get the id of the circuit
  Int_t    GetIdCircuit() const { return fIdCircuit; };
  /// get detector side (Left=0 , Right=1)
  Int_t GetSide() const;
  /// get line
  Int_t GetLine() const;
  /// get column
  Int_t GetCol() const;

  /// set the working status of this board
  void SetStatus(UShort_t s) { fStatus = s; };
  /// set the standard name of this board
  void SetBoardName(const Char_t *name) { fName = new TString(name); };
  /// set the name of the crate containing this board
  void SetCrateName(const Char_t *name) { fCrateName = new TString(name); };
  /// set the number of the detector element containing this board
  void SetDetElemId(Int_t id) { fDetElemId = id; };
  /// set the number of this board
  void SetNumber(Int_t id) { fID = id; }

  /// add a mapping x-pad
  void AddPadX(const AliMpPad &pad, Int_t ich) 
  { 
    new ((*fPadsX[ich])[fNPadsX[ich]++]) AliMpPad(pad); 
  }
  /// add a mapping y-pad
  void AddPadY(const AliMpPad &pad, Int_t ich) 
  { 
    new ((*fPadsY[ich])[fNPadsY[ich]++]) AliMpPad(pad); 
  }
  /// create the display geometry from the mapping pads
  void MakeGeometry();

  /// set an x-strip digit in a chamber with amplitude = amp
  void SetDigitX(Int_t imt, Int_t is, Int_t amp) { 
    fXDig[imt][is] = (UChar_t)amp; }; 
  /// set a  y-strip digit in a chamber with amplitude = amp
  void SetDigitY(Int_t imt, Int_t is, Int_t amp) { 
    fYDig[imt][is] = (UChar_t)amp; }; 
  /// get neighbouring boards with common y strips
  UChar_t GetYOver() const  { return fYOver; };
  /// get the board position inside the detector element in y direction
  UChar_t GetPosition() const { return fPosition; };
  /// get the digit amplitude for an x-strip in a given chamber
  Int_t GetXDig(Int_t imt, Int_t is) const { return fXDig[imt][is]; };
  /// get the digit amplitude for a  y-strip in a given chamber
  Int_t GetYDig(Int_t imt, Int_t is) const { return fYDig[imt][is]; };

  /// set x-strip box for display
  void SetXDigBox(Int_t imt, Int_t is, Double_t x1, Double_t y1, Double_t x2, Double_t y2);
  /// set y-strip box for display
  void SetYDigBox(Int_t imt, Int_t is, Double_t x1, Double_t y1, Double_t x2, Double_t y2);

  /// get x-strip box for display
  TBox *GetXDigBox(Int_t imt, Int_t is) const { return fXDigBox[imt][is]; };
  /// get y-strip box for display
  TBox *GetYDigBox(Int_t imt, Int_t is) const { return fYDigBox[imt][is]; };

  /// get x-center of the board in chamber imt
  Float_t GetXCenter(Int_t imt) const { return fXCenter[imt]; };
  /// get y-center of the board in chamber imt
  Float_t GetYCenter(Int_t imt) const { return fYCenter[imt]; };
  /// get z-center of the board in chamber imt
  Float_t GetZCenter(Int_t imt) const { return fZCenter[imt]; };
  /// get x-width of the board in chamber imt
  Float_t GetXWidth(Int_t imt)  const { return fXWidth[imt]; };
  /// get y-width of the board in chamber imt
  Float_t GetYWidth(Int_t imt)  const { return fYWidth[imt]; };

  /// get x-index in detector element for an x-strip
  Int_t GetXSix()  const { return fXSix;  };
  /// get first y-index in detector element for an x-strip
  Int_t GetXSiy1() const { return fXSiy1; };
  /// get last  y-index in detector element for an x-strip
  Int_t GetXSiy2() const { return fXSiy2; };
  /// get first x-index in detector element for a  y-strip
  Int_t GetYSix1() const { return fYSix1; };
  /// get last  x-index in detector element for a  y-strip
  Int_t GetYSix2() const { return fYSix2; };
  /// get y-index in detector element for a y-strip
  Int_t GetYSiy()  const { return fYSiy;  };
  /// get number of x strips
  Int_t GetNStripX() const { return GetXSiy2() - GetXSiy1() + 1; };
  /// get number of y strips
  Int_t GetNStripY() const { return GetYSix2() - GetYSix1() + 1; };

  /// set true if this board has a gui active
  void   SetOpen(Bool_t open) { fIsOpen = open; };
  /// true if this board has a gui active
  Bool_t IsOpen() const       { return fIsOpen; };

  /// delete the set x-digits
  void  ClearXDigits();
  /// delete the set y-digits
  void  ClearYDigits();

  /// print information on this board
  void PrintBoard() const;

private:

  enum { kNMT = 4, kNS = 16 };     ///< constants

  /// Not implemented
  AliMUONTriggerGUIboard (const AliMUONTriggerGUIboard& board);
  /// Not implemented
  AliMUONTriggerGUIboard& operator=(const AliMUONTriggerGUIboard& board);

  TString       *fName;            ///< Board name LCxLxBx or RCxLxBx
  TString       *fCrateName;       ///< Crate name
  Int_t          fID;              ///< Board array number
  UShort_t       fStatus;          ///< Board status
  UChar_t        fPosition;        ///< Y-boards position
  UChar_t        fYOver;           ///< Y-boards with common y-strips

  Float_t        fXCenter[kNMT];   ///< X-center of the board
  Float_t        fYCenter[kNMT];   ///< Y-center of the board
  Float_t        fZCenter[kNMT];   ///< Z-center of the board
  Float_t        fXWidth[kNMT];    ///< X-width  of the board
  Float_t        fYWidth[kNMT];    ///< Y-width  of the board

  Int_t          fXSix;            ///< X-strips ix index in the board
  Int_t          fXSiy1;           ///< X-strips first iy index in the board
  Int_t          fXSiy2;           ///< X-strips last  iy index in the board

  Int_t          fYSix1;           ///< Y-strips first ix index in the board
  Int_t          fYSix2;           ///< Y-strips last  ix index in the board
  Int_t          fYSiy;            ///< Y-strips iy index in the board

  Int_t          fDetElemId;       ///< Detector element ID (modulo 100)

  Int_t          fIdCircuit;       ///< Circuit number

  UChar_t        fXDig[kNMT][kNS]; ///< X-digits amplitude, set by GUI
  UChar_t        fYDig[kNMT][kNS]; ///< Y-digits amplitude, set by GUI

  TBox          *fXDigBox[kNMT][kNS]; ///< X-digits boxes
  TBox          *fYDigBox[kNMT][kNS]; ///< Y-digits boxes

  Bool_t         fIsOpen;          ///< Selection flag for the digits map

  /// adding pads from mapping to calculate the board geometry
  Int_t fNPadsX[kNMT];             ///< nr of added mapping pads in x
  Int_t fNPadsY[kNMT];             ///< nr of added mapping pads in y
  TClonesArray *fPadsX[kNMT];      ///< array of mapping pads in x
  TClonesArray *fPadsY[kNMT];      ///< array of mapping pads in y

  ClassDef(AliMUONTriggerGUIboard,2) //Trigger GUI utility class: single board object

};

#endif
