#ifndef ALIOBJECTLOADER_H
#define ALIOBJECTLOADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////
//                                        //
//  class AliObjectLoader                 //
//                                        //
//                                        //
////////////////////////////////////////////

class TObject;
class TFolder;
class TString;
class AliDataLoader;

#include "AliBaseLoader.h"

class AliObjectLoader: public AliBaseLoader
{
 public:
     AliObjectLoader(){};
     AliObjectLoader(const TString& name, AliDataLoader* dl, Bool_t storeontop = kFALSE);
     virtual          ~AliObjectLoader(){};
     TObject*          Get() const;

 protected:
     TFolder*          GetFolder() const;
     Int_t             AddToBoard(TObject* obj);
     void              RemoveFromBoard(TObject* obj);

 private:
     AliObjectLoader(const AliObjectLoader&);            //Not implemented
     AliObjectLoader& operator=(const AliObjectLoader&); //Not implemented


 ClassDef(AliObjectLoader,1)    
  
};
#endif


