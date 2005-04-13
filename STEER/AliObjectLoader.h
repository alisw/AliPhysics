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

#include "AliBaseLoader.h"

class AliObjectLoader: public AliBaseLoader
 {
   public:
     AliObjectLoader(){};
     AliObjectLoader(const TString& name, AliDataLoader* dl, Bool_t storeontop = kFALSE);
     AliObjectLoader(const AliObjectLoader& source);
     AliObjectLoader& operator=(const AliObjectLoader& source);
     virtual          ~AliObjectLoader(){};
     TObject*          Get() const;

   protected:
     TFolder*          GetFolder() const;
     Int_t             AddToBoard(TObject* obj);
     void              RemoveFromBoard(TObject* obj);

 ClassDef(AliObjectLoader,1)    
  
};
#endif


