#ifndef ALIXMLCOLLECTION_H
#define ALIXMLCOLLECTION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliXMLCollection
//   This is the class that creates XML collections after querying the tags
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>
#include <Riostream.h>
#include <TString.h>

class TEventList;

//___________________________________________________________________________
class AliXMLCollection : public TObject {
 
 public:
  AliXMLCollection();
  ~AliXMLCollection();
  
  //____________________________________________________//
  Bool_t WriteHeader();
  Bool_t WriteBody(Int_t counter, const char* guid, const char *turl, TEventList *fEventList);
  Bool_t Export();

  void SetCollectionName(const char* name) {fCollectionName = name;}
  
  //____________________________________________________//
  const char* GetCollectionName() {return fCollectionName.Data();}

  //____________________________________________________//
 protected:
  TString fCollectionName;   //the name of the xml file
  ofstream fout; // The output stream
  
  ClassDef(AliXMLCollection,0)  //(ClassName, ClassVersion)
};
//___________________________________________________________________________

#endif
