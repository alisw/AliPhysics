#ifndef ALI_MUON_ST1_ELECTRONIC_ELEMENT_H
#define ALI_MUON_ST1_ELECTRONIC_ELEMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1ElectronicElement
// ---------------------------------
// Describes a set of pads either by defining
// a range of indices, or
// a range of (centimeters) positions or
// a range of electronic channel numbers or
// a range of MANU numbers or, finally,
// a range of gassiplex/MANAS numbers, in a given range of MANU addresses

#include <TObject.h>

class AliMpPad;

class AliMUONSt1ElectronicElement : public TObject 
{
  public:
    enum TDescription {kNone, kIJ, kXY, kMGC, kMG, kM};

  public:
    AliMUONSt1ElectronicElement();
    AliMUONSt1ElectronicElement(TDescription descr);
    virtual ~AliMUONSt1ElectronicElement();
    
    // methods
    Bool_t Contains(const AliMpPad& pad) const;
    void   SetRange(Int_t numVar,Int_t i1,Int_t i2);
    void   SetRange(Int_t numVar,Double_t x1,Double_t x2);
    Bool_t IsInRange(Int_t numVar,Int_t i) const;
    Bool_t IsInRange(Int_t numVar,Double_t x) const;

  private:
    typedef union {Int_t i; Double_t x;}  TData;
    
    TDescription fDescription; // how the pad range is described
    TData fRanges[2][2];       // range of the 2 variables

  ClassDef(AliMUONSt1ElectronicElement,1) //range of electronic elements
};
#endif //ALI_MUON_ST1_ELECTRONIC_ELEMENT_H

