/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpManuGeo.h,v 1.5 2006/05/24 13:58:16 ivana Exp $

/// \ingroup mptrigger
/// \class AliMpLocalBoard
/// \brief Class that manages the properties of the local board
///
/// \author Ch. Finck; Subatech Nantes

#ifndef ALI_MP_LOCAL_BOARD_H
#define ALI_MP_LOCAL_BOARD_H

#include <TNamed.h>

#include "AliMpArrayI.h"
#include "AliMpEncodePair.h"

class TString;

class AliMpLocalBoard : public TNamed
{

 public:
    AliMpLocalBoard(Int_t id, const Char_t* name, Int_t slot);
    AliMpLocalBoard(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpLocalBoard();
    
    // detElem
    Bool_t AddDE(Int_t detElemId);
    Int_t  GetNofDEs() const;
    Int_t  GetDEId(Int_t index) const;
    Int_t  GetDEIdByChamber(Int_t chamberId) const;
    Bool_t HasDEId(Int_t detElemId) const;

    // get methods
    //
           /// Return the identifier (unique)
    Int_t  GetId()   const {return fId;}
           /// Return the slot Identifier in the given crate
    Int_t  GetSlot() const {return fSlot;}

    /// set slot
    void   SetSlot(Int_t slot) {fSlot = slot;}

    // Switches
    //
            /// Get switch bit wise (return a inteter for backware compatibility)
    Int_t  GetSwitch(Int_t index) const;
            /// Set switch in a compact way
    void   SetSwitch(UInt_t swit);
            /// Get switch in a compact way
    UInt_t  GetSwitch() const {return fSwitch;}

    // switch enum for local board (see PRR, chpt: 2.4.4)
    enum {kX2d, kX2m, kX2u, ///< (1) indicate a change of strip pitch in Y circuit
	  kOR0, kOR1,  ///< taking into account the different segmentation in Y from MT1 to MT2
	  kENY,        ///< (0) enable communication in Y to n+/-1 board via tranverse connector, (1) disable
	  kZeroAllYLSB,///< (1) reset the LSB for special configuration of board RC2L5B4 & RC2L6B1
	  kZeroDown,   ///< (0) information is expected from n-1 board for X input, (1) not
	  kZeroMiddle, ///< (0) always, not used
	  kZeroUp };   ///< (0) information is expected from n+1 board for X input, (1) not

    // Transverse connector
    //     
             /// Set transverse connector
    void     SetTC(Bool_t connect) {fTC = connect;}
             /// Return transverse connector
    Bool_t   GetTC() const {return fTC;}

    // Crate name
    //
             /// Set crate name
    void     SetCrate(TString name) {fCrate = name;}
             /// Return crate name
    TString  GetCrate() const {return fCrate;}

    // Notify
    //
             /// Set notified flag (not copy card)
    void     SetNotified(Bool_t notify) {fNotified = notify;}
             /// Return notified flag (not copy card) 
    Bool_t   IsNotified() const {return fNotified;}

    
    /// given position (line, col)
    MpPair_t GetPosition() const;

    // Id to be copy to or from
    
    /// Get Id from where the X input are copied
    Int_t GetInputXfrom() const {return fInputXfrom;}
    /// Get Id to where the X input are copied
    Int_t GetInputXto() const   {return fInputXto;}
    /// Get Id from where the Y input are copied
    Int_t GetInputYfrom() const {return fInputYfrom;}
    /// Get Id to where the Y input are copied
    Int_t GetInputYto() const  {return fInputYto;}

    /// Set Id from where the X input are copied 
    void SetInputXfrom(Int_t id) {fInputXfrom = id;}
    /// Set Id to where the X input are copied    
    void SetInputXto(Int_t id)   {fInputXto   = id;}
    /// Set Id from where the Y input are copied 
    void SetInputYfrom(Int_t id) {fInputYfrom = id;}
    /// Set Id to where the Y input are copied 
    void SetInputYto(Int_t id)   {fInputYto   = id;}

 private:
  /// Not implemented
   AliMpLocalBoard();
  /// Not implemented
   AliMpLocalBoard(const AliMpLocalBoard& src);
   /// Not implemented
   AliMpLocalBoard& operator = (const AliMpLocalBoard& src) ;

   Int_t GetIndex(Int_t chamberId) const;
   
   Int_t       fId;       ///< Identifier (unique)
   Int_t       fSlot;     ///< Slot Identifier in the given crate 
   Bool_t      fTC;       ///< Transverse connector
   TString     fCrate;    ///< Crate name
   UInt_t      fSwitch;   ///< switches in compact way
   Bool_t      fNotified; ///< notified flag (not copy card)
   AliMpArrayI fDEId;     ///< list of Detection element to which this local board is connected
   Int_t       fInputXfrom;///< local id of x3-4 inputs copied from (zero: not copied)
   Int_t       fInputXto;  ///< local id of x3-4 inputs copied to (zero: not copied)
   Int_t       fInputYfrom;///< local id of y1-4 inputs copied from (zero: not copied)
   Int_t       fInputYto;  ///< local id of y1-4 inputs copied to (zero: not copied)

  ClassDef(AliMpLocalBoard,3) //utility class for the motif type
};


#endif 
