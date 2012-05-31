#ifndef ALIBITPACKING_H
#define ALIBITPACKING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a service class for packing and unpacking bits in a 32 bit word.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>


class AliBitPacking: public TObject {
  public :
    static Bool_t  PackWord(UInt_t data, UInt_t &word, 
			    Int_t startBit, Int_t stopBit);
    static UInt_t  UnpackWord(UInt_t word, Int_t startBit, Int_t stopBit);

    ClassDef(AliBitPacking, 0) // class for packing and unpacking bits
};

#endif
