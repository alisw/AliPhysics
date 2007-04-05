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
#include "TGridCollection.h"
#include <Riostream.h>
#include <TString.h>

class TMap;
class TIter;
class TEntryList;

//___________________________________________________________________________
class AliXMLCollection : public TGridCollection {
 
 public:
  AliXMLCollection();
  AliXMLCollection(const char *localCollectionFile);
  AliXMLCollection(const AliXMLCollection& collection);
  
  virtual ~AliXMLCollection();
  
  //____________________________________________________//
  Bool_t WriteHeader();
  Bool_t WriteBody(Int_t counter, const char* guid, const char *lfn, const char *turl, TEntryList *fEntryList);
  Bool_t Export();

  void SetCollectionName(const char* name) {fCollectionName = name;}
  
  //____________________________________________________//
  const char* GetCollectionName() {return fCollectionName.Data();}

  //____________________________________________________//
  void        Reset();
  TMap       *Next();
  Bool_t      Remove(TMap *map);
  const char *GetTURL(const char *name) const;
  const char *GetLFN(const char *name) const;
  const char *GetGUID(const char *name) const;
  TEntryList *GetEventList(const char *filename) const;
  Bool_t      OverlapCollection(AliXMLCollection * comparator);

  static AliXMLCollection *Open(const char *localcollectionfile);

  //____________________________________________________//
 protected:
  virtual void ParseXML();

  TString  fXmlFile;        // collection XML file
  TList   *fEventList;      // list with event file maps
  TIter   *fEventListIter;  // event file list iterator
  TMap    *fCurrent;        // current event file map
  TString  fCollectionName;   //the name of the xml file
  ofstream fout; // The output stream
  
  AliXMLCollection & operator=(const AliXMLCollection & ) {return *this;}

  ClassDef(AliXMLCollection,0)  //(ClassName, ClassVersion)
};
//___________________________________________________________________________

#endif
