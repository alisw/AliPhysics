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

#include <TVector2.h>
#include <TClonesArray.h>

class AliMpPad : public TObject
{
 public:
  AliMpPad(Int_t manuId, Int_t channel,
           Int_t ix, Int_t iy,
           const TVector2& position, const TVector2& dimensions,
	   Bool_t validity = true);
  AliMpPad(Int_t manuId, Int_t channel,
           MpPair_t indices,
           const TVector2& position, const TVector2& dimensions,
	   Bool_t validity = true);

  AliMpPad();
  AliMpPad(const AliMpPad& src);
  virtual ~AliMpPad();

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
  virtual void Print(const char* /*option*/ = "") const;

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
  
               /// Return the pad position (in cm)
  TVector2     Position()    const {return fPosition  ;}
               /// Return the pad dimensions (in cm)
  TVector2     Dimensions()  const {return fDimensions;}
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
  MpPair_t*       fLLocations;     ///<  collection of pad locations - encoded pair (localBoardId, localBoardChannel) 
  UInt_t          fNofLocations;   ///<  number of locations in fLocations
  MpPair_t        fLLocation;      ///<  pad location as encoded pair (manuId, manuChannel) 
  MpPair_t        fLIndices;       ///<  pad indices as encoded pair (ix, iy)  
  TVector2        fPosition;       ///<  the pad position (in cm)
  TVector2        fDimensions;     ///<  the pad dimensions (in cm)
  Bool_t          fValidity;       ///<  validity

  ClassDef(AliMpPad,3) //utility class for the motif type
};

ostream& operator << (ostream &out, const AliMpPad& op);
Bool_t operator < (const AliMpPad& left, const AliMpPad& right);

#endif //ALI_MP_PAD_H
