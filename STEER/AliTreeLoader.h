#ifndef ALITREELOADER_H
#define ALITREELOADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////
//                                        //
//  class AliTreeLoader                   //
//                                        //
//  Loader responsible for one data type  //
//  i.e. Hits, Kine, etc.                 //
//  many objects type can be assciated    //
//  with one data type: storing object    //
//  (usually tree), task producing it,    //
//  Quality Assurance(QA), QA Task, and   //
//  others.                               //
//                                        //
//                                        //
////////////////////////////////////////////

class TString;
class TTree;

#include "AliObjectLoader.h"

class AliTreeLoader: public AliObjectLoader
 {
   public:
     AliTreeLoader(){};
     AliTreeLoader(const TString& name, AliDataLoader* dl, Bool_t storeontop = kFALSE);
     AliTreeLoader(const AliTreeLoader& source);
     AliTreeLoader& operator=(const AliTreeLoader& source);
     virtual ~AliTreeLoader(){};
     
     virtual TTree*     Tree() const {return dynamic_cast<TTree*>(Get());}
     virtual void       MakeTree();
     virtual Int_t      WriteData(Option_t* opt="");

   ClassDef(AliTreeLoader,1)    
 };

#endif


