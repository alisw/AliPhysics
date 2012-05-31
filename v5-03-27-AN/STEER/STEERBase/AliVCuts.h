#ifndef AliVCuts_H
#define AliVCuts_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Event cuts base class
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>

class AliVCuts : public TNamed {

 public :
  AliVCuts(); 
  AliVCuts(const char* name, const char* title); 
  virtual ~AliVCuts() { };
  AliVCuts(const AliVCuts& evt); 
  AliVCuts& operator=(const AliVCuts& evt);
  virtual Bool_t IsSelected(TObject* /* obj  */)  = 0;
  virtual UInt_t GetSelectionMask(const TObject* /* obj */) { return 0; }
  virtual TObject *GetStatistics(Option_t *) const {return 0;}
  ClassDef(AliVCuts,1);
};

#endif
