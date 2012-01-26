/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id: AliMergeableCollection.cxx 50593 2011-07-14 17:42:28Z martinez $

///
/// A mergeable object container. 
///
/// For each tuple (key1,key2,..,keyN) a (hash)list of mergeable objects is associated.
/// Note that key1, key2 (optional), ..., keyN (optional) are strings. 
/// Those strings should not contain "/" themselves.
///
/// More helper functions might be added in the future (e.g. Project, etc...)

#include "AliMergeableCollection.h"

ClassImp(AliMergeableCollection)

#include "AliLog.h"
#include "Riostream.h"
#include "TError.h"
#include "THashList.h"
#include "TKey.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
//#include "TH2.h"

//_____________________________________________________________________________
AliMergeableCollection::AliMergeableCollection(const char* name, const char* title) 
: TNamed(name,title), fMap(0x0), fMustShowEmptyObject(0), fMapVersion(0), fMessages()
{
  /// Ctor
}

//_____________________________________________________________________________
AliMergeableCollection::~AliMergeableCollection()
{
  /// dtor. Note that the map is owner
  if ( fMap ) fMap->DeleteAll();
  delete fMap;
}

//_____________________________________________________________________________
Bool_t 
AliMergeableCollection::Adopt(TObject* obj)
{
  /// Adopt a given object at top level (i.e. no key)
  return InternalAdopt("",obj);
}

//_____________________________________________________________________________
Bool_t 
AliMergeableCollection::Adopt(const char* identifier, TObject* obj)
{
  /// Adopt a given object, and associate it with pair key
  TString sidentifier(identifier);
  if ( ! sidentifier.IsNull() ){
    if ( ! sidentifier.EndsWith("/") ) sidentifier.Append("/");
    if ( ! sidentifier.BeginsWith("/") ) sidentifier.Prepend("/");
  }
  return InternalAdopt(sidentifier.Data(),obj);
}

//_____________________________________________________________________________
void AliMergeableCollection::ClearMessages()
{
  /// clear pending messages
  fMessages.clear();
}

//_____________________________________________________________________________
TIterator*
AliMergeableCollection::CreateIterator(Bool_t direction) const
{
  /// Create an iterator (must be deleted by the client)
  return fMap ? new AliMergeableCollectionIterator(this,direction) : 0x0;
}

//_____________________________________________________________________________
AliMergeableCollection*
AliMergeableCollection::Clone(const char* name) const
{
  /// Clone this collection.
  /// We loose the messages.
  
  AliMergeableCollection* newone = new AliMergeableCollection(name,GetTitle());
  
  newone->fMap = static_cast<TMap*>(fMap->Clone());
  newone->fMustShowEmptyObject = fMustShowEmptyObject;
  newone->fMapVersion = fMapVersion;  
  
  return newone;
}

//_____________________________________________________________________________
void 
AliMergeableCollection::Delete(Option_t*)
{
  /// Delete all the objects
  fMap->DeleteAll();
  delete fMap;
  fMap=0x0;
}

//_____________________________________________________________________________
TObject* 
AliMergeableCollection::FindObject(const char* fullIdentifier) const
{
  /// Find an object by its full identifier.
  
  return GetObject(fullIdentifier);
}

//_____________________________________________________________________________
TObject* 
AliMergeableCollection::FindObject(const TObject *object) const
{
  /// Find an object 
  AliWarning("This method is awfully inefficient. Please improve it or use FindObject(const char*)");
  TIter next(CreateIterator());
  TObject* obj;
  while ( ( obj=next() ) )
  {
    if ( obj->IsEqual(object) ) return obj;
  }
  return 0x0;
}


//_____________________________________________________________________________
TList*
AliMergeableCollection::CreateListOfKeys(Int_t index) const
{
  /// Create the list of keys at level index
  
  TList* list = new TList;
  list->SetOwner(kTRUE);
  
  TObjArray* ids = SortAllIdentifiers();
  TIter next(ids);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString oneid = GetKey(str->String().Data(),index,kFALSE);
    if (oneid.Length()>0 && !list->Contains(oneid))
    {
      list->Add(new TObjString(oneid));
    }
  }
  
  delete ids;
  return list;
}

//_____________________________________________________________________________
TList* 
AliMergeableCollection::CreateListOfObjectNames(const char* identifier) const
{
  /// Create list of object names for /key1/key2/key...
  /// Returned list must be deleted by client
  
  TList* listOfNames = new TList;
  listOfNames->SetOwner(kTRUE);
  
  TIter next(Map());
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString currIdentifier = str->String();
    if ( currIdentifier.CompareTo(identifier) ) continue;
    
    THashList* list = static_cast<THashList*>(Map()->GetValue(identifier));
    
    TIter nextObject(list);
    TObject* obj;
    
    while ( ( obj = nextObject() ) )
    {
      listOfNames->Add(new TObjString(obj->GetName()));
    }    
  }
  
  return listOfNames;
}


//_____________________________________________________________________________
TString
AliMergeableCollection::GetIdentifier(const char* fullIdentifier) const
{
  /// Extract the identifier from the fullIdentifier
  TString sfullIdentifier(fullIdentifier);
  TObjArray* arr = sfullIdentifier.Tokenize("/");
  TString identifier = "";
  for ( Int_t istr=0; istr<arr->GetLast(); istr++ ) {
    identifier += "/" + GetKey(fullIdentifier, istr, kTRUE);
  }
  delete arr;
  identifier.Append("/");
  return identifier;
}

//_____________________________________________________________________________
TString
AliMergeableCollection::GetKey(const char* identifier, Int_t index, Bool_t idContainsObjName) const
{
  /// Extract the index element of the key pair from the fullIdentifier
  TString sidentifier(identifier);
  if ( ! idContainsObjName ) sidentifier.Append("/dummy");
  return InternalDecode(sidentifier.Data(),index);
}

//_____________________________________________________________________________
TString
AliMergeableCollection::GetObjectName(const char* fullIdentifier) const
{
  /// Extract the object name from an identifier
  
  return InternalDecode(fullIdentifier,-1);  
}


//_____________________________________________________________________________
TObject* 
AliMergeableCollection::GetObject(const char* fullIdentifier) const
{
  /// Get object key1/key2/.../objectName:action
  
  TObjArray* arr = TString(fullIdentifier).Tokenize(":");
  
  TString fullIdWithoutAction(static_cast<TObjString*>(arr->At(0))->String());
  TString action;
  
  if ( arr->GetLast() > 0 ) 
  {
    action = static_cast<TObjString*>(arr->At(1))->String();
    action.ToUpper();
  }
  
  delete arr;
  
  return GetObject(GetIdentifier(fullIdWithoutAction).Data(), GetObjectName(fullIdWithoutAction));
  
  
  //  if (obj)
  //  {
  //    TH2* h2(0x0);
  //    
  //    if ( action == "PX" && ( (h2 = dynamic_cast<TH2*>(obj)) ) ) 
  //    {
  //      return h2->ProjectionX(NormalizeName(identifier.Data(),action.Data()).Data());
  //    }
  //    else if ( action == "PY" && ( (h2 = dynamic_cast<TH2*>(obj)) ) ) 
  //    {
  //      return h2->ProjectionY(NormalizeName(identifier.Data(),action.Data()).Data());
  //    }
  //    else if ( action == "PFX" && ( (h2 = dynamic_cast<TH2*>(obj)) ) ) 
  //    {
  //      return h2->ProfileX(NormalizeName(identifier.Data(),action.Data()).Data());
  //    }
  //    else if ( action == "PFY" && ( (h2 = dynamic_cast<TH2*>(obj)) ) ) 
  //    {
  //      return h2->ProfileY(NormalizeName(identifier.Data(),action.Data()).Data());
  //    }
  //    
  //  }
  //  else
  //  {
  //    AliDebug(1,Form("Object %s not found",sidentifier));
  //  }
  //  return obj;
}


//_____________________________________________________________________________
TObject* 
AliMergeableCollection::GetObject(const char* identifier, 
                                  const char* objectName) const
{
  /// Get object for (identifier,objectName) triplet
  
  TString sidentifier(identifier);
  if ( ! sidentifier.IsNull() ) {
    if ( ! sidentifier.BeginsWith("/") ) sidentifier.Prepend("/");
    if ( ! sidentifier.EndsWith("/") ) sidentifier.Append("/");
  }
  return InternalObject(sidentifier.Data(),objectName);
}

//_____________________________________________________________________________
TObject* AliMergeableCollection::GetSum(const char* idPattern)
{
  /// Sum objects
  /// The pattern must be in the form:
  /// /key1_pattern1 key1_pattern2,key1_pattern,.../key2_pattern1,key2_pattern,.../.../objectName_pattern1,objectName_pattern...
  /// The logical or between patterns separated by commas is taken
  
  TObject* sumObject = 0x0;
  
  // Build array of lists of pattern
  TString idPatternString(idPattern);
  TObjArray* keyList = idPatternString.Tokenize("/");
  TObjArray keyMatrix(keyList->GetEntries());
  keyMatrix.SetOwner();
  for ( Int_t ikey=0; ikey<keyList->GetEntries(); ikey++ ) {
    TObjArray* subKeyList = ((TObjString*)keyList->At(ikey))->GetString().Tokenize(",");
    keyMatrix.AddAt(subKeyList, ikey);
  }
  delete keyList;
  
  TString debugMsg = "Adding objects:";
  
  TIter next(Map());
  TObjString* str;
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString identifier = str->String();
    
    Bool_t listMatchPattern = kTRUE;
    for ( Int_t ikey=0; ikey<keyMatrix.GetEntries()-1; ikey++ ) {
      TString currKey = GetKey(identifier, ikey, kFALSE);
      Bool_t matchKey = kFALSE;
      TObjArray* subKeyList = static_cast<TObjArray*> ( keyMatrix.At(ikey) );
      for ( Int_t isub=0; isub<subKeyList->GetEntries(); isub++ ) {
        TString subKeyString = static_cast<TObjString*> (subKeyList->At(isub))->GetString();
        if ( currKey.Contains(subKeyString.Data()) ) {
          matchKey = kTRUE;
          break;
        }
      } // loop on the list of patterns of each key
      if ( ! matchKey ) {
        listMatchPattern = kFALSE;
        break;
      }
    } // loop on keys in the idPattern
    if ( ! listMatchPattern ) continue;
    
    THashList* list = static_cast<THashList*>(Map()->GetValue(identifier.Data()));
    
    TIter nextObj(list);
    TObject* obj;
    
    while ( ( obj = nextObj()) )
    {
      TString currKey = obj->GetName();
      Bool_t matchKey = kFALSE;
      TObjArray* subKeyList = static_cast<TObjArray*> ( keyMatrix.Last() );
      for ( Int_t isub=0; isub<subKeyList->GetEntries(); isub++ ) {
        TString subKeyString = static_cast<TObjString*> (subKeyList->At(isub))->GetString();
        if ( currKey.Contains(subKeyString.Data()) ) {
          matchKey = kTRUE;
          break;
        }
      }
      if ( ! matchKey ) continue;
      if ( ! sumObject ) sumObject = obj->Clone();
      else MergeObject(sumObject, obj);
      debugMsg += Form(" %s%s",identifier.Data(),obj->GetName());
    } // loop on objects in list
  } // loop on identifiers in map
  
  AliDebug(1,debugMsg.Data());
  
  return sumObject;
}

//_____________________________________________________________________________
Bool_t AliMergeableCollection::InternalAdopt(const char* identifier, TObject* obj)
{
  /// Adopt an obj
  
  if (!obj)
  {
    Error("Adopt","Cannot adopt a null object");
    return kFALSE;
  }
  
  if ( ! obj->IsA()->InheritsFrom(TObject::Class()) ||
        ! obj->IsA()->GetMethodWithPrototype("Merge", "TCollection*") ) {
    Error("Adopt","Cannot adopt an object which is not mergeable!"); 
  }
  
  THashList* hlist = 0x0;
  
  hlist = static_cast<THashList*>(Map()->GetValue(identifier));
  
  if (!hlist)
  {
    hlist = new THashList;
    hlist->SetOwner(kTRUE);
    Map()->Add(new TObjString(identifier),hlist);
    hlist->SetName(identifier);
  }
  
  TObject* existingObj = hlist->FindObject(obj->GetName());
  
  if (existingObj)
  {
    AliError(Form("Cannot adopt an already existing object : %s -> %s",identifier,existingObj->GetName()));
    return kFALSE;
  }
  
  if ( obj->IsA()->InheritsFrom(TH1::Class()) ) (static_cast<TH1*> ( obj ))->SetDirectory(0);  
  
  hlist->AddLast(obj);
  
  return kTRUE;
  
}

//_____________________________________________________________________________
TString
AliMergeableCollection::InternalDecode(const char* identifier, Int_t index) const
{
  /// Extract the index-th element of the identifier (/key1/key2/.../keyN/objectName)
  /// object is index=-1 (i.e. last)
  
  if ( strlen(identifier) > 0 && identifier[0] != '/' ) 
  {    
    AliError(Form("identifier %s is malformed (should start with /)",identifier));
    return "";
  }
  
  TObjArray* array = TString(identifier).Tokenize("/");

  if ( index >= array->GetLast() ) 
  {
    AliError(Form("Requiring index %i of identifier %s which only have %i",index, identifier, array->GetLast()));
    delete array;
    return "";
  }

  TString value("");
  
  if ( index < 0 ) 
  {
    value = static_cast<TObjString*>(array->Last())->String();    
  }
  else if ( index <= array->GetLast() ) 
  {
    value = static_cast<TObjString*>(array->At(index))->String();
  }
  
  delete array;
  
  return value;
}

//_____________________________________________________________________________
TObject* 
AliMergeableCollection::InternalObject(const char* identifier,
                                       const char* objectName) const
{
  /// Get object for (identifier,objectName) 
  
  if (!fMap) 
  {
    return 0x0;
  }
  
  THashList* hlist = static_cast<THashList*>(Map()->GetValue(identifier));
  if (!hlist) 
  {
    TString msg(Form("Did not find hashlist for identifier=%s dir=%s",identifier,gDirectory ? gDirectory->GetName() : "" ));
    fMessages[msg.Data()]++;
    return 0x0;
  }
  
  TObject* obj = hlist->FindObject(objectName);
  if (!obj)
  {
    TString msg(Form("Did not find objectName=%s in %s",objectName,identifier));
    fMessages[msg.Data()]++;
  }
  return obj;
}


//_____________________________________________________________________________
Bool_t AliMergeableCollection::IsEmptyObject(TObject* obj) const
{
  /// Check if object is empty
  /// (done only for TH1, so far)
    
  if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
    TH1* histo = static_cast<TH1*> (obj);
    if ( histo->GetEntries() == 0 ) return kTRUE;
  }

  return kFALSE;

}


//_____________________________________________________________________________
TMap* AliMergeableCollection::Map() const
{
  /// Wrapper to insure proper key formats (i.e. new vs old)
  
  if (!fMap)
  {
    fMap = new TMap;
    fMap->SetOwnerKeyValue(kTRUE,kTRUE);
    fMapVersion = 1;
  }
  else
  {
    if ( fMapVersion < 1 ) 
    {
      AliInfo("Remapping");
      // change the keys
      TIter next(fMap);
      TObjString* str;
      
      while ( ( str = static_cast<TObjString*>(next()) ) )
      {
        if ( str->String().Contains("./") )
        {
          TString newkey(str->String());
          
          newkey.ReplaceAll("./","");
          
          TObject* o = fMap->GetValue(str);
          
          TPair* p = fMap->RemoveEntry(str);
          if (!p)
          {
            AliError("oups oups oups");
            return 0x0;
          }
          
          fMap->Add(new TObjString(newkey.Data()),o);
          
          delete p;
        }
      }
      
      fMapVersion = 1;
    }
  }
  
  return fMap;
}

//_____________________________________________________________________________
Long64_t
AliMergeableCollection::Merge(TCollection* list)
{
  // Merge a list of AliMergeableCollection objects with this
  // Returns the number of merged objects (including this).
  
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  TIter next(list);
  TObject* currObj;
  TList mapList;
  Int_t count(0);
  
  while ( ( currObj = next() ) )
  {
    AliMergeableCollection* mergeCol = dynamic_cast<AliMergeableCollection*>(currObj);
    if (!mergeCol) {
      AliFatal(Form("object named \"%s\" is a %s instead of an AliMergeableCollection!", currObj->GetName(), currObj->ClassName()));
      continue;
    }
    
    ++count;
    
    if ( mergeCol->fMap ) mergeCol->Map(); // to insure keys in the new format
    
    TIter nextIdentifier(mergeCol->fMap);
    TObjString* identifier;

    while ( ( identifier = static_cast<TObjString*>(nextIdentifier()) ) )
    {
      THashList* otherList = static_cast<THashList*>(mergeCol->fMap->GetValue(identifier->String().Data()));

      TIter nextObject(otherList);
      TObject* obj;
      
      while ( ( obj = nextObject() ) )
      {
        TString newid(Form("%s%s",identifier->String().Data(),obj->GetName()));
        
        TObject* thisObject = GetObject(newid.Data());
        
        if (!thisObject)
        {
          AliDebug(1,Form("Adopting a new object = %s%s",identifier->String().Data(),obj->GetName()));
          
          Bool_t ok = Adopt(identifier->String(), obj->Clone());
          
          if (!ok)
          {
            AliError(Form("Adoption of object %s failed",obj->GetName()));
          }
        }
        else
        {
          // add it...
          AliDebug(1,Form("Merging object = %s%s",
                          identifier->String().Data(),
                          obj->GetName()));
          
          MergeObject(thisObject, obj);
        }
      } // loop on objects in map
    } // loop on identifiers
  } // loop on collections in list
         
  AliDebug(1,Form("count=%d",count));
  
  return count+1;
}

//_____________________________________________________________________________
Bool_t AliMergeableCollection::MergeObject(TObject* baseObject, TObject* objToAdd)
{
  /// Add objToAdd to baseObject
  
  if ( baseObject->IsA()->Class() != objToAdd->IsA()->Class() ) {
    printf("MergeObject: Cannot add %s to %s", objToAdd->ClassName(), baseObject->ClassName());
    return kFALSE;
  }
  if ( ! baseObject->IsA()->InheritsFrom(TObject::Class()) ||
      ! baseObject->IsA()->GetMethodWithPrototype("Merge", "TCollection*") ) {
    printf("MergeObject: Objects are not mergeable!");
    return kFALSE;
  }  
  
  TList list;
  list.Add(objToAdd);
  
  TString listArgs = Form("((TCollection*)0x%lx)", (ULong_t)&list);
  Int_t error = 0;
  baseObject->Execute("Merge", listArgs.Data(), &error);
  return kTRUE;
}

//_____________________________________________________________________________
TString AliMergeableCollection::NormalizeName(const char* identifier,const char* action) const
{
  /// Replace / by _ to build a root-compliant histo name
  TString name(GetName());
  
  name += "_";
  name += identifier;
  name += "_";
  name += action;
  name.ReplaceAll("/","_");
  name.ReplaceAll("-","_");
  return name;
}

//_____________________________________________________________________________
Int_t 
AliMergeableCollection::NumberOfObjects() const
{
  /// Get the number of objects we hold
  TIter next(CreateIterator());
  Int_t count(0);
  while ( next() ) ++count;
  return count;
}

//_____________________________________________________________________________
Int_t 
AliMergeableCollection::NumberOfKeys() const
{
  /// Get the number of keys we have
  return fMap ? fMap->GetSize() : 0;
}

//_____________________________________________________________________________
void 
AliMergeableCollection::Print(Option_t* option) const
{
  /// Print all the objects we hold, in a hopefully visually pleasing
  /// way.
  ///
  /// Option can be used to select given part only, using the schema :
  /// /*/*/*/*/*
  /// Where the stars are wilcards for /key1/key2/.../objectName
  ///
  /// if * is used it is assumed to be a wildcard for objectName
  ///
  /// For other selections the full syntax /*/*/*/*/* must be used.
  /// 
  /// Use "-" as objectName to disable object's name output
  
  cout << Form("AliMergeableCollection(%s,%s) : %d keys and %d objects",
               GetName(),GetTitle(),
               NumberOfKeys(), NumberOfObjects()) << endl;
  
  if (!strlen(option)) return;
    
  TObjArray* select = TString(option).Tokenize("/");
  
  TString sreObjectName(select->Last()->GetName());
  TRegexp reObjectName(sreObjectName.Data(),kTRUE);
  
  TObjArray* identifiers = SortAllIdentifiers();
  
  printf("identifiers entries %i\n", identifiers->GetEntries());
    
  TIter nextIdentifier(identifiers);
  
  TObjString* sid(0x0);
  
  while ( ( sid = static_cast<TObjString*>(nextIdentifier()) ) )
  {
    Bool_t identifierPrinted(kFALSE);

    TString identifier(sid->String());
    
    Bool_t matchPattern = kTRUE;
    for ( Int_t isel=0; isel<select->GetLast(); isel++ ) {
      if ( ! GetKey(identifier.Data(), isel, kFALSE).Contains(TRegexp(select->At(isel)->GetName(),kTRUE)) ) {
        matchPattern = kFALSE;
        break;
      }
    }
    if ( ! matchPattern ) continue;
    
    if ( sreObjectName == "*" ) {
      identifierPrinted = kTRUE;
      cout << identifier.Data() << endl;
    }
      
    THashList * list = static_cast<THashList*>(Map()->GetValue(sid->String().Data()));      
    TObjArray names;
    names.SetOwner(kTRUE);
    TIter nextUnsortedObj(list);
    TObject* obj;
    while ( ( obj = nextUnsortedObj() ) )
    {
      names.Add(new TObjString(obj->GetName()));
    }
    names.Sort();
    TIter nextObjName(&names);
    TObjString* oname;
    while ( ( oname = static_cast<TObjString*>(nextObjName()) ) )
    {
      TString objName(oname->String());
      if (objName.Contains(reObjectName) )
      {
        obj = list->FindObject(objName.Data());
        if ( IsEmptyObject(obj) && ! fMustShowEmptyObject ) continue;
        if (!identifierPrinted)
        {
          cout << identifier.Data() << endl;
          identifierPrinted = kTRUE;
        }
        cout << Form("    (%s) %s", obj->ClassName(), obj->GetName());
        if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
          TH1* histo = static_cast<TH1*> (obj);
          cout << Form(" %s Entries=%d Sum=%g",histo->GetTitle(),Int_t(histo->GetEntries()),histo->GetSumOfWeights());
        }
        cout << endl;
      }
    }
    if (!identifierPrinted && sreObjectName=="-" )
    { 
      // to handle the case where we used objectName="-" to disable showing the objectNames,
      // but we still want to see the matching keys maybe...
      cout << identifier.Data() << endl;
    }
  }
  
  delete select;
  
  delete identifiers;
}

//_____________________________________________________________________________
void 
AliMergeableCollection::PrintMessages(const char* prefix) const
{
  /// Print pending messages
  
  std::map<std::string,int>::const_iterator it;
  
  for ( it = fMessages.begin(); it != fMessages.end(); ++it ) 
  {
    cout << Form("%s : message %s appeared %5d times",prefix,it->first.c_str(),it->second) << endl;
  }
}


//_____________________________________________________________________________
UInt_t 
AliMergeableCollection::EstimateSize(Bool_t show) const
{
  /// Estimate the memory (in kilobytes) used by some objects

//  For TH1:
//  sizeof(TH1) + (nbins+2)*(nbytes_per_bin) +name+title_sizes 
//  if you have errors add (nbins+2)*8 
    
  TIter next(CreateIterator());
  
  TObject* obj;
  UInt_t size(0);
  
  while ( ( obj = next() ) )
  {
    UInt_t thissize=0;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      TH1* histo = static_cast<TH1*> (obj);
      Int_t nbins = (histo->GetNbinsX()+2);
    
      if (histo->GetNbinsY()>1)
      {
        nbins *= (histo->GetNbinsY()+2);
      }
    
      if (histo->GetNbinsZ()>1)
      {
        nbins *= (histo->GetNbinsZ()+2);
      }
      
      Bool_t hasErrors = ( histo->GetSumw2N() > 0 );
    
      TString cname(histo->ClassName());
    
      Int_t nbytesPerBin(0);
    
      if (cname.Contains(TRegexp("C$")) ) nbytesPerBin = sizeof(Char_t);
      if (cname.Contains(TRegexp("S$")) ) nbytesPerBin = sizeof(Short_t);
      if (cname.Contains(TRegexp("I$")) ) nbytesPerBin = sizeof(Int_t);
      if (cname.Contains(TRegexp("F$")) ) nbytesPerBin = sizeof(Float_t);
      if (cname.Contains(TRegexp("D$")) ) nbytesPerBin = sizeof(Double_t);
        
      if (!nbytesPerBin)
      {
        AliError(Form("Could not get the number of bytes per bin for histo %s of class %s. Thus the size estimate will be wrong !",
                      histo->GetName(),histo->ClassName()));
        continue;
      }
    
      thissize = sizeof(histo) + nbins*(nbytesPerBin) + strlen(histo->GetName())
      + strlen(histo->GetTitle());
      
      if ( hasErrors) thissize += nbins*8;
    }
    else {
      AliWarning(Form("Cannot estimate size of %s\n", obj->ClassName()));
      continue;
    }

    size += thissize;
    
    if ( show ) 
    {
      AliInfo(Form("Size of %30s is %20d bytes",obj->GetName(),thissize));
    }
  } // loop on objects

  return size;
}

//_____________________________________________________________________________
void AliMergeableCollection::PruneEmptyObjects()
{
  /// Delete the empty objects
  /// (Implemented for TH1 only)
  TIter next(Map());
  TObjString* key;
  
  TList toBeRemoved;
  toBeRemoved.SetOwner(kTRUE);
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    TString identifier(key->String());
    THashList* hlist = static_cast<THashList*>(Map()->GetValue(identifier.Data()));
    TIter nextObject(hlist);
    TObject* obj;
    while ( ( obj = nextObject() ) )
    {
      if ( IsEmptyObject(obj) ) toBeRemoved.Add(new TObjString(Form("%s%s",identifier.Data(),obj->GetName())));
    }
  }
  
  TIter nextTBR(&toBeRemoved);
  while ( ( key = static_cast<TObjString*>(nextTBR()) ) )
  {
    Remove(key->GetString().Data());
    AliDebug(2,Form("Removing %s", key->GetString().Data()));
  }
}

//_____________________________________________________________________________
AliMergeableCollection* 
AliMergeableCollection::Project(const char* identifier) const
{
  /// To be implemented : would create a new collection starting at /key1/key2/...
  
  if (!fMap) return 0x0;
  
  AliMergeableCollection* mergCol = new AliMergeableCollection(Form("%s %s",GetName(),identifier),
                                                               GetTitle());
  
  TIter next(Map());
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString currIdentifier = str->String();
    if ( ! currIdentifier.Contains(identifier) ) continue;
    
    THashList* list = static_cast<THashList*>(Map()->GetValue(identifier));
    
    TIter nextObj(list);
    TObject* obj;
    
    while ( ( obj = nextObj()) )
    {
      TObject* clone = obj->Clone();

      TString newkey(currIdentifier.Data());
      newkey.ReplaceAll(identifier,"");

      if (newkey=="/") newkey="";
      
      mergCol->InternalAdopt(newkey.Data(),clone);
    }    
  }

  return mergCol;
}

//_____________________________________________________________________________
TObject* 
AliMergeableCollection::Remove(const char* fullIdentifier)
{
  ///
  /// Remove a given object (given its fullIdentifier=/key1/key2/.../objectName)
  ///
  /// Note that we do *not* remove the /key1/key2/... entry even if there's no
  /// more object for this triplet.
  ///
  /// Not very efficient. Could be improved ?
  ///
  
  TString identifier = GetIdentifier(fullIdentifier);
  
  THashList* hlist = dynamic_cast<THashList*>(Map()->GetValue(identifier.Data()));
  
  if (!hlist)
  {
    AliWarning(Form("Could not get hlist for key=%s",identifier.Data()));
    return 0x0;
  }
    
  TObject* obj = GetObject(fullIdentifier);
  if (!obj)
  {
    AliError(Form("Could not find object %s",fullIdentifier));
    return 0x0;
  }
  
  TObject* rmObj = hlist->Remove(obj);
  if (!rmObj)
  {
    AliError("Remove failed");
    return 0x0;
  }
  
  return rmObj;
}

//_____________________________________________________________________________
TObjArray*
AliMergeableCollection::SortAllIdentifiers() const
{
  /// Sort our internal identifiers. Returned array must be deleted.
  TObjArray* identifiers = new TObjArray;
  identifiers->SetOwner(kFALSE); 
  TIter next(Map());
  TObjString* sid;
  
  while ( ( sid = static_cast<TObjString*>(next()) ) )
  {
    if ( !identifiers->FindObject(sid->String().Data()) )
    {
      identifiers->Add(sid);      
    }
  }
  identifiers->Sort();
  return identifiers;
}


///////////////////////////////////////////////////////////////////////////////
//
// AliMergeableCollectionIterator
//
///////////////////////////////////////////////////////////////////////////////

class AliMergeableCollectionIterator;

//_____________________________________________________________________________
AliMergeableCollectionIterator::AliMergeableCollectionIterator(const AliMergeableCollection* mcol, Bool_t dir)
: fkMergeableCollection(mcol), fMapIterator(0x0), fHashListIterator(0x0), fDirection(dir)
{
  /// Default ctor
}

//_____________________________________________________________________________
AliMergeableCollectionIterator&
AliMergeableCollectionIterator::operator=(const TIterator&)
{
  /// Overriden operator= (imposed by Root's declaration of TIterator ?)
  Fatal("TIterator::operator=","Not implementeable"); // because there's no clone in TIterator :-(
  return *this;
}

//_____________________________________________________________________________
AliMergeableCollectionIterator::~AliMergeableCollectionIterator()
{
  /// dtor
  Reset();
}

//_____________________________________________________________________________
TObject* AliMergeableCollectionIterator::Next()
{
  /// Advance to next object in the collection
  
  if (!fHashListIterator)
  {
    if ( !fMapIterator ) 
    {
      fMapIterator = fkMergeableCollection->fMap->MakeIterator(fDirection);
    }
    TObjString* key = static_cast<TObjString*>(fMapIterator->Next());
    if (!key)
    {
      // we are done
      return 0x0;
    }      
    THashList* list = static_cast<THashList*>(fkMergeableCollection->Map()->GetValue(key->String().Data()));
    if (!list) return 0x0;
    fHashListIterator = list->MakeIterator(fDirection);
  }

  TObject* obj = fHashListIterator->Next();
  
  if (!obj) 
  {
    delete fHashListIterator;
    fHashListIterator = 0x0;
    return Next();
  }
  
  return obj;
}

//_____________________________________________________________________________
void AliMergeableCollectionIterator::Reset()
{
  /// Reset the iterator
  delete fHashListIterator;
  delete fMapIterator;
}
