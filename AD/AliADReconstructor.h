#ifndef ALIADRECONSTRUCTOR_H
#define ALIADRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.*/
/* See cxx source for full Copyright notice                              */
/* $Id: AliADReconstructor.h 20956 2007-09-26 14:22:18Z cvetan $  */

///////////////////////////////////////////////////////////////////////////
///                                                                      //
/// class for AD reconstruction                                       //
///                                                                      //
///////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliLog.h"

class AliESDAD;
class AliESDEvent;

class AliADReconstructor: public AliReconstructor {
public:
  AliADReconstructor();
  virtual ~AliADReconstructor();
  virtual void   Init();
  
  virtual void   Reconstruct(AliRawReader* /*rawReader*/, 
		             TTree* /*clustersTree*/) const {
    AliError("Method not implemented"); return;};
  virtual void   Reconstruct(TTree*, TTree*) const {return;};
  
  virtual void   FillESD(TTree* digitsTree, TTree* /*clustersTree*/, 
			 AliESDEvent* esd) const;

  virtual void   FillESD(AliRawReader* /*rawReader*/, TTree* /*clustersTree*/, 
			 AliESDEvent* /*esd*/) const { 
    AliError("Method not implemented"); return;};
  
  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  virtual void ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;


protected:

  AliESDAD*        fESDAD;      // AD ESD object  

private:
  AliADReconstructor(const AliADReconstructor&); //Not implemented
  AliADReconstructor& operator = (const AliADReconstructor&); //Not implemented
  

  mutable TClonesArray *fDigitsArray;  // clones-array for ConvertDigits() and FillESD()

  ClassDef(AliADReconstructor, 1)  // class for the AD reconstruction
};

#endif
