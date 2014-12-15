/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpPad.h,v 1.11 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpPad
/// \brief Class which encapsuate all information about a pad
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_PAD_H
#define ALI_MP_PAD_H

#include "AliMpEncodePair.h"

#include <TObject.h>

#include <TClonesArray.h>

class AliMpPad : public TObject
{
 public:
  AliMpPad(Int_t manuId, Int_t channel,
           Int_t ix, Int_t iy,
           Double_t x,  Double_t y, 
           Double_t dx,  Double_t dy,
	   Bool_t validity = true);
  AliMpPad(Int_t manuId, Int_t channel,
           MpPair_t indices,
           Double_t positionX,  Double_t positionY, 
           Double_t dx,  Double_t dy,
	   Bool_t validity = true);

  AliMpPad();
  AliMpPad(const AliMpPad& src);
  ~AliMpPad();

  //
  // operators  
  //
  Bool_t operator == (const AliMpPad& pos2) const;
  Bool_t operator != (const AliMpPad& pos2) const;
  AliMpPad& operator = (const AliMpPad& src) ;
  
  //
  // methods
  //
  void PrintOn(ostream& out) const;
  void Print(const char* /*option*/ = "") const;

  //
  // static get methods
  //
               /// Return invalid pad
  static AliMpPad Invalid() {return AliMpPad();}

  //
  // set methods
  //
  Bool_t  AddLocation(Int_t localBoardId, Int_t localBoardChannel, 
                      Bool_t warn = true);

  //
  // get methods
  //
               /// Return pad location as encoded pair (manuId, manuChannel)
  MpPair_t     GetLocation() const { return fLLocation; }
  Int_t        GetManuId() const;
  Int_t        GetManuChannel() const;
  
               /// Return pad indices as encoded pair (ix, iy)
  MpPair_t     GetIndices()  const { return fLIndices; }
  Int_t        GetIx() const;
  Int_t        GetIy() const;
  
               /// Return the pad x position (in cm)
  Double_t     GetPositionX() const { return fPositionX; }
               /// Return the pad x position (in cm)
  Double_t     GetPositionY() const { return fPositionY; }
  
               /// Return the x pad dimension - half length (in cm)
  Double_t     GetDimensionX()  const {return fDimensionX;}
               /// Return the y pad dimension - half length (in cm)
  Double_t     GetDimensionY()  const {return fDimensionY;}

               /// Return validity
  Bool_t       IsValid()     const {return fValidity  ;}
  
  Int_t        GetNofLocations() const;
  MpPair_t     GetLocation(Int_t i) const;  
  Int_t        GetLocalBoardId(Int_t i) const;
  Int_t        GetLocalBoardChannel(Int_t i) const;

  Bool_t       HasLocation(Int_t localBoardId, Int_t localBoardChannel) const; 

 private:

  // static data members
  static const Int_t  fgkMaxNofLocations; ///< \brief maximum number of pad locations
                                          /// in the collection
  // data members
  UInt_t          fNofLocations;   ///<  number of locations in fLocations
  /// Collection of pad locations - encoded pair (localBoardId, localBoardChannel) 
  MpPair_t*       fLLocations;     //[fNofLocations]
  MpPair_t        fLLocation;      ///<  pad location as encoded pair (manuId, manuChannel) 
  MpPair_t        fLIndices;       ///<  pad indices as encoded pair (ix, iy)  
  Double_t        fPositionX;      ///<  the pad x position (in cm)
  Double_t        fPositionY;      ///<  the pad y position (in cm)
  Double_t        fDimensionX;     ///<  the pad x dimension - half length (in cm)
  Double_t        fDimensionY;     ///<  the pad y dimension - half length(in cm)
  Bool_t          fValidity;       ///<  validity

  ClassDef(AliMpPad,4) //utility class for the motif type
};

ostream& operator << (ostream &out, const AliMpPad& op);
Bool_t operator < (const AliMpPad& left, const AliMpPad& right);

#endif //ALI_MP_PAD_H
