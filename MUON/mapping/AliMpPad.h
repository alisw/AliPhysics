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

#include <TObject.h>

#include "AliMpIntPair.h"

#include <TVector2.h>
#include <TClonesArray.h>

class AliMpPad : public TObject
{
 public:
  AliMpPad(const AliMpIntPair& location, const AliMpIntPair& indices,
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
  Bool_t  AddLocation(const AliMpIntPair& location, Bool_t warn = true);

  //
  // get methods
  //
               /// Return pad location
  AliMpIntPair GetLocation() const {return fLocation;}
               /// Return pad indices
  AliMpIntPair GetIndices()  const {return fIndices;}
               /// Return the pad position (in cm)
  TVector2     Position()    const {return fPosition  ;}
               /// Return the pad dimensions (in cm)
  TVector2     Dimensions()  const {return fDimensions;}
               /// Return validity
  Bool_t       IsValid()     const {return fValidity  ;}
  
  Int_t        GetNofLocations() const;
  AliMpIntPair GetLocation(Int_t i) const;
  Bool_t       HasLocation(const AliMpIntPair& location) const; 

 private:
  // static data members
  static const Int_t  fgkMaxNofLocations; ///< \brief maximum number of pad locations
                                          /// in the collection
  // data members
  AliMpIntPair*   fLocations;      ///<  collection of pad locations 
  UInt_t          fNofLocations;   ///<  number of locations in fLocations
  AliMpIntPair    fLocation;       ///<  pad location
  AliMpIntPair    fIndices;        ///<  pad indices
  TVector2        fPosition;       ///<  the pad position (in cm)
  TVector2        fDimensions;     ///<  the pad dimensions (in cm)
  Bool_t          fValidity;       ///<  validity

  ClassDef(AliMpPad,2) //utility class for the motif type
};

ostream& operator << (ostream &out, const AliMpPad& op);
Bool_t operator < (const AliMpPad& left, const AliMpPad& right);

#endif //ALI_MP_PAD_H
