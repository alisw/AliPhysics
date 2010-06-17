#ifndef ALITRACKLETER_H
#define ALITRACKLETER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//   class Alitrackleter
//   An abstract interface for tracklet reconstruction
//-------------------------------------------------------------------------

#include <TObject.h>
class TTree;
class AliESDEvent;
class AliMultiplicity;

class AliTrackleter : public TObject {
public:
 AliTrackleter() :fMult(0) {}
  virtual ~AliTrackleter() { delete fMult; }
  virtual void Reconstruct(AliESDEvent* esd, TTree* treeRP) = 0;
  virtual AliMultiplicity* GetMultiplicity() const {return fMult;}
  //
protected:
 AliTrackleter(const AliTrackleter &src) : TObject(src) {}
  AliTrackleter & operator=(const AliTrackleter &src) {if (&src!=this) TObject::operator=(src); return *this;}
  
 protected:
  AliMultiplicity* fMult;   // multiplicity object

  ClassDef(AliTrackleter,1) //base trackleter
};

#endif
