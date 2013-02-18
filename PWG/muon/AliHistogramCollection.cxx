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

// $Id$

///
/// A (fake) "4D" histogram container. 
///
/// For each tuple (keyA,keyB,keyC,keyD) a (hash)list of histogram is associated.
/// Note that keyA, keyB (optional), keyC (optional) and keyD (optional) are strings. 
/// Those strings should not contain "/" themselves.
///
/// More helper functions might be added in the future (e.g. Project, etc...)

#include "AliHistogramCollection.h"

ClassImp(AliHistogramCollection)

#include "AliLog.h"
#include "AliMergeableCollection.h"
#include "Riostream.h"
#include "TError.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TKey.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TSystem.h"

using std::cout;
using std::endl;
ClassImp(AliHistogramCollection)

//_____________________________________________________________________________
AliHistogramCollection::AliHistogramCollection(const char* name, const char* title) 
: TNamed(name,title), fMap(0x0), fMustShowEmptyHistogram(kFALSE), fMapVersion(0), fMessages()
{
  /// Ctor
}

//_____________________________________________________________________________
AliHistogramCollection::~AliHistogramCollection()
{
  /// dtor. Note that the map is owner
  if ( fMap ) fMap->DeleteAll();
  delete fMap;
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(TH1* histo)
{
  /// Adopt a given histogram at top level (i.e. no key)
  return InternalAdopt("",histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* key, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA)
  return InternalAdopt(Form("/%s/",key),histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* keyA, const char* keyB, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA,keyB)
  return InternalAdopt(Form("/%s/%s/",keyA,keyB),histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* keyA, const char* keyB, const char* keyC, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA,keyB,keyC)
  return InternalAdopt(Form("/%s/%s/%s/",keyA,keyB,keyC),histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* keyA, const char* keyB, const char* keyC, const char* keyD, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA,keyB,keyC,keyD)
  return InternalAdopt(Form("/%s/%s/%s/%s/",keyA,keyB,keyC,keyD),histo);
}

//_____________________________________________________________________________
void AliHistogramCollection::ClearMessages()
{
  /// clear pending messages
  fMessages.clear();
}

//_____________________________________________________________________________
TIterator*
AliHistogramCollection::CreateIterator(Bool_t direction) const
{
  /// Create an iterator (must be deleted by the client)
  return fMap ? new AliHistogramCollectionIterator(this,direction) : 0x0;
}

//_____________________________________________________________________________
AliHistogramCollection*
AliHistogramCollection::Clone(const char* name) const
{
  /// Clone this collection.
  /// We loose the messages.
  
  AliHistogramCollection* newone = new AliHistogramCollection(name,GetTitle());
  
  newone->fMap = static_cast<TMap*>(fMap->Clone());
  newone->fMustShowEmptyHistogram = fMustShowEmptyHistogram;
  newone->fMapVersion = fMapVersion;  
  
  return newone;
}

//_____________________________________________________________________________
void 
AliHistogramCollection::Delete(Option_t*)
{
  /// Delete all the histograms
  fMap->DeleteAll();
  delete fMap;
  fMap=0x0;
}

//_____________________________________________________________________________
TObject* 
AliHistogramCollection::FindObject(const char* identifier) const
{
  /// Find an object by its full identifier.
  
  return Histo(KeyA(identifier),
               KeyB(identifier),
               KeyC(identifier),
               KeyD(identifier),
               HistoName(identifier));
}

//_____________________________________________________________________________
TObject* 
AliHistogramCollection::FindObject(const TObject *key) const
{
  /// Find an object 
  AliWarning("This method is awfully inefficient. Please improve it or use FindObject(const char*)");
  TIter next(CreateIterator());
  TObject* o;
  while ( ( o=next() ) )
  {
    if ( o->IsEqual(key) ) return o;
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMergeableCollection* AliHistogramCollection::Convert() const
{
  /// Convert this into a mergeable collection so it can receive objects
  /// that are not histograms
  
  AliMergeableCollection* mc = new AliMergeableCollection(GetName(),GetTitle());
  
  TObjArray* ids = SortAllIdentifiers();
  TIter next(ids);
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    THashList* list = static_cast<THashList*>(fMap->GetValue(key->String().Data()));
    TIter nextHisto(list);
    TH1* histo;
    
    while ( ( histo = static_cast<TH1*>(nextHisto())))
    {
      TH1* h = static_cast<TH1*>(histo->Clone());
      mc->Adopt(key->String().Data(),h);
    }
  }
  
  delete ids;
  
  return mc;
}


//_____________________________________________________________________________
TList*
AliHistogramCollection::CreateListOfKeysA() const
{
  /// Create list of keys at level "A" (/A/B/C/D/histoname)
  /// Must be delete by client
  return CreateListOfKeys(0);
}

//_____________________________________________________________________________
TList*
AliHistogramCollection::CreateListOfKeysB() const
{
  /// Create list of keys at level "B" (/A/B/C/D/histoname)
  /// Must be delete by client
  return CreateListOfKeys(1);
}

//_____________________________________________________________________________
TList*
AliHistogramCollection::CreateListOfKeysC() const
{
  /// Create list of keys at level "C" (/A/B/C/D/histoname)
  /// Must be delete by client
  return CreateListOfKeys(2);
}

//_____________________________________________________________________________
TList*
AliHistogramCollection::CreateListOfKeysD() const
{
  /// Create list of keys at level "D" (/A/B/C/D/histoname)
  /// Must be delete by client
  return CreateListOfKeys(3);
}

//_____________________________________________________________________________
TList* 
AliHistogramCollection::CreateListOfHistogramNames(const char* keyA, const char* keyB, const char* keyC, const char* keyD) const
{
  /// Create list of histogram names for /keyA/keyB/keyC/KeyD
  /// Returned list must be deleted by client
  
  TList* listOfNames = new TList;
  listOfNames->SetOwner(kTRUE);
  
  TIter next(Map());
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString identifier = str->String();
    
    Bool_t copy= ( KeyA(identifier) == keyA &&
                  ( strlen(keyB)==0 || KeyB(identifier)==keyB ) && 
                  ( strlen(keyC)==0 || KeyC(identifier)==keyC ) && 
                  ( strlen(keyD)==0 || KeyD(identifier)==keyD ));
    
    if ( !copy ) continue;
    
    THashList* list = static_cast<THashList*>(Map()->GetValue(identifier.Data()));
    
    TIter nextHisto(list);
    TH1* h;
    
    while ( ( h = static_cast<TH1*>(nextHisto()) ) )
    {
      listOfNames->Add(new TObjString(h->GetName()));
    }    
  }
  
  return listOfNames;
}

//_____________________________________________________________________________
TList*
AliHistogramCollection::CreateListOfKeys(Int_t index) const
{
  /// Create the list of keys at level index
  
  TList* list = new TList;
  list->SetOwner(kTRUE);
  
  TObjArray* ids = SortAllIdentifiers();
  TIter next(ids);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString oneid = InternalDecode(str->String().Data(),index);
    if (oneid.Length()>0 && !list->Contains(oneid))
    {
      list->Add(new TObjString(oneid));
    }
  }
  
  delete ids;
  return list;
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* keyA, 
                              const char* histoname) const
{
  /// Get histo for (keyA,histoname) triplet
  
  return InternalHisto(Form("/%s/",keyA),histoname);
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* keyA, const char* keyB, 
                              const char* histoname) const
{
  /// Get histo for (keyA,keyB,histoname) triplet

  return InternalHisto(Form("/%s/%s/",keyA,keyB),histoname);
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* keyA, const char* keyB, const char* keyC,
                              const char* histoname) const
{
  /// Get histo for (keyA,keyB,keyC,histoname) quad
  
  return InternalHisto(Form("/%s/%s/%s/",keyA,keyB,keyC),histoname);
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* keyA, const char* keyB, 
                              const char* keyC, const char* keyD,
                              const char* histoname) const
{
  /// Get histo for (keyA,keyB,keyC,histoname) quad
  
  return InternalHisto(Form("/%s/%s/%s/%s/",keyA,keyB,keyC,keyD),histoname);
}

//_____________________________________________________________________________
TString
AliHistogramCollection::HistoName(const char* identifier) const
{
  /// Extract the histogram name from an identifier
  
  return InternalDecode(identifier,-1);  
}

//_____________________________________________________________________________
Bool_t AliHistogramCollection::InternalAdopt(const char* identifier, TH1* histo)
{
  /// Adopt an histogram
  
  if (!histo)
  {
    Error("Adopt","Cannot adopt a null histogram");
    return kFALSE;
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
  
  TH1* h = static_cast<TH1*>(hlist->FindObject(histo->GetName()));
  
  if (h)
  {
    AliError(Form("Cannot adopt an already existing histogram : %s -> %s",identifier,h->GetName()));
    return kFALSE;
  }
  
  
  histo->SetDirectory(0);  
  
  hlist->AddLast(histo);
  
  return kTRUE;
  
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* sidentifier) const
{
  /// Get histogram keyA/keyB/keyC/keyD/histoname:action
  
  TObjArray* a = TString(sidentifier).Tokenize(":");
  
  TString identifier(static_cast<TObjString*>(a->At(0))->String());
  TString action;
  
  if ( a->GetLast() > 0 ) 
  {
    action = static_cast<TObjString*>(a->At(1))->String();
    action.ToUpper();
  }

  delete a;

  Int_t nslashes = identifier.CountChar('/');

  TH1* h(0x0);
  
  switch (nslashes)
  {
    case 0:
    case 1:
      // no slash : the identifier is just the histogram name
      h = InternalHisto("",identifier);
      break;      
    case 2:
      h = Histo(InternalDecode(identifier,0),
                InternalDecode(identifier,-1));
      break;
    case 3:
      h = Histo(InternalDecode(identifier,0),
                InternalDecode(identifier,1),
                InternalDecode(identifier,-1));
      break;
    case 4:
      h = Histo(InternalDecode(identifier,0),
                InternalDecode(identifier,1),
                InternalDecode(identifier,2),
                InternalDecode(identifier,-1));
      break;
    case 5:
      h = Histo(InternalDecode(identifier,0),
                InternalDecode(identifier,1),
                InternalDecode(identifier,2),
                InternalDecode(identifier,3),
                InternalDecode(identifier,-1));
      break;
    default:
      AliError(Form("Invalid identifier %s",identifier.Data()));
      break;
  }
  
  if (h)
  {
    TH2* h2(0x0);
    
    if ( action == "PX" && ( (h2 = dynamic_cast<TH2*>(h)) ) ) 
    {
      return h2->ProjectionX(NormalizeName(identifier.Data(),action.Data()).Data());
    }
    else if ( action == "PY" && ( (h2 = dynamic_cast<TH2*>(h)) ) ) 
    {
      return h2->ProjectionY(NormalizeName(identifier.Data(),action.Data()).Data());
    }
    else if ( action == "PFX" && ( (h2 = dynamic_cast<TH2*>(h)) ) ) 
    {
      return h2->ProfileX(NormalizeName(identifier.Data(),action.Data()).Data());
    }
    else if ( action == "PFY" && ( (h2 = dynamic_cast<TH2*>(h)) ) ) 
    {
      return h2->ProfileY(NormalizeName(identifier.Data(),action.Data()).Data());
    }
    
  }
  else
  {
    AliDebug(1,Form("histogram %s not found",sidentifier));
  }
  return h;
}

//_____________________________________________________________________________
TString
AliHistogramCollection::InternalDecode(const char* identifier, Int_t index) const
{
  /// Extract the index-th element of the identifier (/keyA/keyB/keyC/keyD/histoname)
  /// keyA is index=0
  /// keyB is index=1
  /// keyC is index=2
  /// keyD is index=3
  /// histo is index=-1 (i.e. last)
  
  if ( strlen(identifier) > 0 && identifier[0] != '/' ) 
  {    
    AliError(Form("identifier %s is malformed (should start with /)",identifier));
    return "";
  }
  
  TObjArray* array = TString(identifier).Tokenize("/");

  if ( array->GetLast()>5 ) 
  {
    AliError(Form("identifier %s is malformed (more than 5 /)",identifier));
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
TH1* 
AliHistogramCollection::InternalHisto(const char* identifier,
                                      const char* histoname) const
{
  /// Get histo for (identifier,histoname) 
  
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
  
  TH1* h = static_cast<TH1*>(hlist->FindObject(histoname));  
  if (!h)
  {
    TString msg(Form("Did not find histoname=%s in %s",histoname,identifier));
    fMessages[msg.Data()]++;
  }
  return h;
}

//_____________________________________________________________________________
TString
AliHistogramCollection::KeyA(const char* identifier) const
{
  /// Extract the first element of the key pair from an identifier
  
  return InternalDecode(identifier,0);
}

//_____________________________________________________________________________
TString
AliHistogramCollection::KeyB(const char* identifier) const
{
  /// Extract the second element (if present) 
  return InternalDecode(identifier,1);
}

//_____________________________________________________________________________
TString
AliHistogramCollection::KeyC(const char* identifier) const
{
  /// Extract the 3rd element (if present) 
  return InternalDecode(identifier,2);
}

//_____________________________________________________________________________
TString
AliHistogramCollection::KeyD(const char* identifier) const
{
  /// Extract the 4th element (if present) 
  return InternalDecode(identifier,3);
}

//_____________________________________________________________________________
TMap* AliHistogramCollection::Map() const
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
AliHistogramCollection::Merge(TCollection* list)
{
  // Merge a list of AliHistogramCollection objects with this
  // Returns the number of merged objects (including this).
  
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  TIter next(list);
  TObject* o;
  TList mapList;
  Int_t count(0);
  
  while ( ( o = next() ) )
  {
    AliHistogramCollection* hcol = dynamic_cast<AliHistogramCollection*>(o);
    if (!hcol) {
      AliFatal(Form("object named \"%s\" is a %s instead of an AliHistogramCollection!", o->GetName(), o->ClassName()));
      continue;
    }
    
    ++count;
    
    if ( hcol->fMap ) hcol->Map(); // to insure keys in the new format
    
    TIter nextIdentifier(hcol->fMap);
    TObjString* identifier;

    while ( ( identifier = static_cast<TObjString*>(nextIdentifier()) ) )
    {
      THashList* otherList = static_cast<THashList*>(hcol->fMap->GetValue(identifier->String().Data()));

      TIter nextHisto(otherList);
      TH1* h;
      
      while ( ( h = static_cast<TH1*>(nextHisto()) ) )
      {
        TString newid(Form("%s%s",identifier->String().Data(),h->GetName()));
        
        TH1* thisHisto = Histo(newid.Data());
        
        if (!thisHisto)
        {
          // this is an histogram we don't have yet. Let's add it
          
          Int_t nslashes = TString(newid).CountChar('/');
          
          Bool_t ok(kFALSE);
          
          switch (nslashes)
          {
            case 0:
              ok = Adopt(static_cast<TH1*>(h->Clone()));
              break;
            case 2:
              ok = Adopt(KeyA(identifier->String()),
                        static_cast<TH1*>(h->Clone()));
              break;
            case 3:
              ok = Adopt(KeyA(identifier->String()),
                        KeyB(identifier->String()),
                        static_cast<TH1*>(h->Clone()));
              break;
            case 4:
              ok = Adopt(KeyA(identifier->String()),
                        KeyB(identifier->String()),
                        KeyC(identifier->String()),
                        static_cast<TH1*>(h->Clone()));
              break;
            case 5:
              ok = Adopt(KeyA(identifier->String()),
                        KeyB(identifier->String()),
                        KeyC(identifier->String()),
                        KeyD(identifier->String()),
                        static_cast<TH1*>(h->Clone()));
              break;
            default:
              AliError(Form("Invalid identifier : %s",identifier->String().Data()));
              break;
          }
          
          if (!ok)
          {
            AliError(Form("Adoption of histogram %s failed",h->GetName()));
          }
        }
        else
        {
          // add it...
          if ( HistoSameAxis(h,thisHisto) )
          {
            thisHisto->Add(h);
          }
          else
          {
            TList l;
            l.Add(h);
          
            thisHisto->Merge(&l);
          }
        }
      }
    }
  }
  
  return count+1;
}

//_____________________________________________________________________________
TString AliHistogramCollection::NormalizeName(const char* identifier,const char* action) const
{
  // Replace / by _ to build a root-compliant histo name
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
AliHistogramCollection::NumberOfHistograms() const
{
  /// Get the number of histograms we hold
  TIter next(CreateIterator(this));
  Int_t n(0);
  while ( next() ) ++n;
  return n;
}

//_____________________________________________________________________________
Int_t 
AliHistogramCollection::NumberOfKeys() const
{
  /// Get the number of keys we have
  return fMap ? fMap->GetSize() : 0;
}

//_____________________________________________________________________________
void 
AliHistogramCollection::Print(Option_t* option) const
{
  /// Print all the histograms we hold, in a hopefully visually pleasing
  /// way.
  ///
  /// Option can be used to select given part only, using the schema :
  /// /*/*/*/*/*
  /// Where the stars are wilcards for /keyA/keyB/keyC/KeyD/histoname
  ///
  /// if * is used it is assumed to be a wildcard for histoname
  ///
  /// For other selections the full syntax /*/*/*/*/* must be used.
  /// 
  /// Use "-" as histoname to disable histogram's name output
  
  cout << Form("AliHistogramCollection(%s,%s) : %d keys and %d histos",
               GetName(),GetTitle(),
               NumberOfKeys(), NumberOfHistograms()) << endl;
  
  if (!strlen(option)) return;
  
  TString sreKeyA("*");
  TString sreKeyB("*");
  TString sreKeyC("*");
  TString sreKeyD("*");
  TString sreHistoname("*");
  
  TObjArray* select = TString(option).Tokenize("/");
  Int_t n = select->GetLast();
  
  if (n>=0)
  {
    sreHistoname = static_cast<TObjString*>(select->At(n))->String();
  }
  if (n>=1)
  {
    sreKeyD = static_cast<TObjString*>(select->At(n-1))->String();    
  }
  if (n>=2)
  {
    sreKeyC = static_cast<TObjString*>(select->At(n-2))->String();    
  }
  if (n>=3)
  {
    sreKeyB = static_cast<TObjString*>(select->At(n-3))->String();    
  }
  if (n>=4)
  {
    sreKeyA = static_cast<TObjString*>(select->At(n-3))->String();    
  }
  
  TRegexp reKeyA(sreKeyA,kTRUE);
  TRegexp reKeyB(sreKeyB,kTRUE);
  TRegexp reKeyC(sreKeyC,kTRUE);
  TRegexp reKeyD(sreKeyD,kTRUE);
  TRegexp reHistoname(sreHistoname,kTRUE);

  delete select;
  
  TObjArray* identifiers = SortAllIdentifiers();
    
  TIter nextIdentifier(identifiers);
  
  TObjString* sid(0x0);
  
  while ( ( sid = static_cast<TObjString*>(nextIdentifier()) ) )
  {
    Bool_t identifierPrinted(kFALSE);

    TString identifier(sid->String());
    
    if ( ( InternalDecode(identifier,0).Contains(reKeyA) &&
          InternalDecode(identifier,1).Contains(reKeyB) &&
          InternalDecode(identifier,2).Contains(reKeyC) &&
          InternalDecode(identifier,3).Contains(reKeyD) )
        )
    {
      if ( sreHistoname == "*" ) 
      {
        identifierPrinted = kTRUE;
        cout << identifier.Data() << endl;
      }
      
      THashList * list = static_cast<THashList*>(Map()->GetValue(sid->String().Data()));      
      TObjArray names;
      names.SetOwner(kTRUE);
      TIter nextUnsortedHisto(list);
      TH1* h;
      while ( ( h = static_cast<TH1*>(nextUnsortedHisto()) ) )
      {
        names.Add(new TObjString(h->GetName()));
      }
      names.Sort();
      TIter nextHistoName(&names);
      TObjString* hname;
      while ( ( hname = static_cast<TObjString*>(nextHistoName()) ) )
      {
        TString histoName(hname->String());
        if (histoName.Contains(reHistoname) )
        {
          h = static_cast<TH1*>(list->FindObject(histoName.Data()));
          if ( h->GetEntries()==0 && !fMustShowEmptyHistogram ) continue;
          if (!identifierPrinted)
          {
            cout << identifier.Data() << endl;
            identifierPrinted = kTRUE;
          }
          cout << Form("    %s %s Entries=%d Sum=%g",h->GetName(),h->GetTitle(),Int_t(h->GetEntries()),h->GetSumOfWeights()) << endl;
        }
      }
      if (!identifierPrinted && sreHistoname=="-" )
      { 
        // to handle the case where we used histoname="-" to disable showing the histonames,
        // but we still want to see the matching keys maybe...
        cout << identifier.Data() << endl;
      }
    }
  }
  
  delete identifiers;
}

//_____________________________________________________________________________
void 
AliHistogramCollection::PrintMessages(const char* prefix) const
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
AliHistogramCollection::EstimateSize(Bool_t show) const
{
  /// Estimate the memory (in bytes) used by our histograms
  
//  sizeof(TH1) + (nbins+2)*(nbytes_per_bin) +name+title_sizes 
//  if you have errors add (nbins+2)*8 
  
  TMap sizeMap;
  sizeMap.SetOwnerKeyValue(kTRUE,kTRUE);
  
  TIter next(CreateIterator());
  TH1* h;
  UInt_t n(0);
  
  while ( ( h = static_cast<TH1*>(next()) ) )
  {
    Int_t nbins = (h->GetNbinsX()+2);
    
    if (h->GetNbinsY()>1)
    {
      nbins *= (h->GetNbinsY()+2);
    }
    
    if (h->GetNbinsZ()>1)
    {
      nbins *= (h->GetNbinsZ()+2);
    }
      
    Bool_t hasErrors = ( h->GetSumw2N() > 0 );
    
    TString cname(h->ClassName());
    
    Int_t nbytesPerBin(0);
    
    if (cname.Contains(TRegexp("C$")) ) nbytesPerBin = sizeof(Char_t);
    if (cname.Contains(TRegexp("S$")) ) nbytesPerBin = sizeof(Short_t);
    if (cname.Contains(TRegexp("I$")) ) nbytesPerBin = sizeof(Int_t);
    if (cname.Contains(TRegexp("F$")) ) nbytesPerBin = sizeof(Float_t);
    if (cname.Contains(TRegexp("D$")) ) nbytesPerBin = sizeof(Double_t);
        
    if (!nbytesPerBin)
    {
      AliError(Form("Could not get the number of bytes per bin for histo %s of class %s. Thus the size estimate will be wrong !",
                    h->GetName(),h->ClassName()));
      continue;
    }
    
    UInt_t thissize = sizeof(h) + nbins*(nbytesPerBin) + strlen(h->GetName())
    + strlen(h->GetTitle());
    
    TObjString* m = static_cast<TObjString*>(sizeMap.GetValue(h->GetName()));
    
    if (!m)
    {
      m = new TObjString(Form("%u %u",thissize,1));
      sizeMap.Add(new TObjString(h->GetName()),m);
    }
    else
    {
      UInt_t s;
      UInt_t d;
      
      sscanf(m->String().Data(),"%u %u",&s,&d);
      
      s += thissize;
      ++d;
      
      m->String().Form("%u %u",s,d);
    }
    
    if ( hasErrors) thissize += nbins*8;

    n += thissize;
    
  }

  if ( show )
  {
    std::multimap<UInt_t,std::string> sorted;
    
    TIter nextMapSize(&sizeMap);
    TObjString* m;

    UInt_t s;
    UInt_t d;

    while ( ( m = static_cast<TObjString*>(nextMapSize()) ) )
    {

      TObjString* sizeMsg = static_cast<TObjString*>(sizeMap.GetValue(m->String()));
      
      sscanf(sizeMsg->String().Data(),"%u %u",&s,&d);

      sorted.insert(std::pair<UInt_t,std::string>(s,m->String().Data()));
    }
    
    std::multimap<UInt_t,std::string>::const_iterator it;
    
    for ( it = sorted.begin(); it != sorted.end(); ++it )
    {
      TObjString* sizeMsg = static_cast<TObjString*>(sizeMap.GetValue(it->second.c_str()));
      
      sscanf(sizeMsg->String().Data(),"%u %u",&s,&d);

      AliInfo(Form("%40s : size %10u (%7u per object, %6d objects)",
                   it->second.c_str(),s,s/d,d));
                   
    }
  }

  return n;
}

//_____________________________________________________________________________
void AliHistogramCollection::PruneEmptyHistograms()
{
  /// Delete the empty histograms
  TIter next(Map());
  TObjString* key;
  
  TList toBeRemoved;
  toBeRemoved.SetOwner(kTRUE);
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    TString identifier(key->String());
    THashList* hlist = static_cast<THashList*>(Map()->GetValue(identifier.Data()));
    TIter nextHisto(hlist);
    TH1* h;
    while ( ( h = static_cast<TH1*>(nextHisto())))
    {
      if ( h->GetEntries()==0)
      {
        toBeRemoved.Add(new TObjString(Form("%s%s",identifier.Data(),h->GetName())));
      }
    }
  }
  
  TIter nextTBR(&toBeRemoved);
  while ( ( key = static_cast<TObjString*>(nextTBR()) ) )
  {
    Remove(key);
  }
}

//_____________________________________________________________________________
AliHistogramCollection* 
AliHistogramCollection::Project(const char* keyA, const char* keyB, const char* keyC, const char* keyD) const
{
  /// Create a new collection starting at keyA/keyB/keyC/keyD
  /// Histograms are *copied*
  
  if (!fMap) return 0x0;
  
  AliHistogramCollection* hc = new AliHistogramCollection(Form("%s %s/%s/%s/%s",GetName(),keyA,keyB,keyC,keyD),
                                                          GetTitle());
  
  TIter next(Map());
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TString identifier = str->String();
    
    Bool_t copy= ( KeyA(identifier) == keyA &&
                  ( strlen(keyB)==0 || KeyB(identifier)==keyB ) && 
                  ( strlen(keyC)==0 || KeyC(identifier)==keyC ) && 
                  ( strlen(keyD)==0 || KeyD(identifier)==keyD ));

    if ( !copy ) continue;
    
    THashList* list = static_cast<THashList*>(Map()->GetValue(identifier.Data()));
    
    TIter nextHisto(list);
    TH1* h;
    
    while ( ( h = static_cast<TH1*>(nextHisto()) ) )
    {
      TH1* hclone = static_cast<TH1*>(h->Clone());

      TString newkey(identifier.Data());
      
      if ( strlen(keyD) > 0 ) newkey.ReplaceAll(Form("/%s",keyD),"");
      if ( strlen(keyC) > 0 ) newkey.ReplaceAll(Form("/%s",keyC),"");
      if ( strlen(keyB) > 0 ) newkey.ReplaceAll(Form("/%s",keyB),"");
      if ( strlen(keyA) > 0 ) newkey.ReplaceAll(Form("/%s",keyA),"");

      if (newkey=="/") newkey="";
      
      hc->InternalAdopt(newkey.Data(),hclone);
    }    
  }

  return hc;
}

//_____________________________________________________________________________
TObject* 
AliHistogramCollection::Remove(TObject* key)
{
  ///
  /// Remove a given histogram (given its key=full identifier=/keyA/keyB/keyC/keyD/histoname)
  ///
  /// Note that we do *not* remove the /keyA/keyB/keyC/keyD entry even if there's no
  /// more histogram for this triplet.
  ///
  /// Not very efficient. Could be improved ?
  ///
  
  TObjString* str = dynamic_cast<TObjString*>(key);
  
  if (!str)
  {
    AliError(Form("key is not of the expected TObjString type, but of %s",key->ClassName()));
    return 0x0;
  }
  
  TString identifier(str->String());
    
  Int_t nslashes = TString(identifier).CountChar('/');

  TString skey;
  
  switch (nslashes )
  {
    case 2:
      skey = Form("/%s/",
                  KeyA(identifier).Data());
      break;
    case 3:
      skey = Form("/%s/%s/",
                  KeyA(identifier).Data(),
                  KeyB(identifier).Data());
      break;
    case 4:
      skey = Form("/%s/%s/%s/",
                  KeyA(identifier).Data(),
                  KeyB(identifier).Data(),
                  KeyC(identifier).Data());
      break;
    case 5:
      skey = Form("/%s/%s/%s/%s/",
                  KeyA(identifier).Data(),
                  KeyB(identifier).Data(),
                  KeyC(identifier).Data(),
                  KeyD(identifier).Data()
                  );
      break;
    default:
      AliError("Oups");
      break;
  };
  
  THashList* hlist = dynamic_cast<THashList*>(Map()->GetValue(skey.Data()));
  
  if (!hlist)
  {
    AliWarning(Form("Could not get hlist for key=%s",skey.Data()));
    return 0x0;
  }
    
  TH1* h = InternalHisto(skey,HistoName(identifier.Data()));
  if (!h)
  {
    AliError(Form("Could not find histo %s",identifier.Data()));
    return 0x0;
  }
  
  TObject* o = hlist->Remove(h);
  if (!o)
  {
    AliError("Remove failed");
    return 0x0;
  }
  
  return o;
}

//______________________________________________________________________________
Bool_t AliHistogramCollection::HistoSameAxis(TH1 *h0, TH1 *h1) const
{
  // shameless copy from TProofPlayerRemote::HistoSameAxis
  //
  // Return kTRUE is the histograms 'h0' and 'h1' have the same binning and ranges
  // on the axis (i.e. if they can be just Add-ed for merging).
  
  Bool_t rc = kFALSE;
  if (!h0 || !h1) return rc;
  
  TAxis *a0 = 0, *a1 = 0;
  
  // Check X
  a0 = h0->GetXaxis();
  a1 = h1->GetXaxis();
  if (a0->GetNbins() == a1->GetNbins())
    if (TMath::Abs(a0->GetXmax() - a1->GetXmax()) < 1.e-9)
      if (TMath::Abs(a0->GetXmin() - a1->GetXmin()) < 1.e-9) rc = kTRUE;
  
  // Check Y, if needed
  if (h0->GetDimension() > 1) {
    rc = kFALSE;
    a0 = h0->GetYaxis();
    a1 = h1->GetYaxis();
    if (a0->GetNbins() == a1->GetNbins())
      if (TMath::Abs(a0->GetXmax() - a1->GetXmax()) < 1.e-9)
        if (TMath::Abs(a0->GetXmin() - a1->GetXmin()) < 1.e-9) rc = kTRUE;
  }
  
  // Check Z, if needed
  if (h0->GetDimension() > 2) {
    rc = kFALSE;
    a0 = h0->GetZaxis();
    a1 = h1->GetZaxis();
    if (a0->GetNbins() == a1->GetNbins())
      if (TMath::Abs(a0->GetXmax() - a1->GetXmax()) < 1.e-9)
        if (TMath::Abs(a0->GetXmin() - a1->GetXmin()) < 1.e-9) rc = kTRUE;
  }
  
  // Done
  return rc;
}

//_____________________________________________________________________________
TObjArray*
AliHistogramCollection::SortAllIdentifiers() const
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
// AliHistogramCollectionIterator
//
///////////////////////////////////////////////////////////////////////////////

class AliHistogramCollectionIterator;

//_____________________________________________________________________________
AliHistogramCollectionIterator::AliHistogramCollectionIterator(const AliHistogramCollection* hcol, Bool_t dir)
: fkHistogramCollection(hcol), fMapIterator(0x0), fHashListIterator(0x0), fDirection(dir)
{
  /// Default ctor
}

//_____________________________________________________________________________
AliHistogramCollectionIterator&
AliHistogramCollectionIterator::operator=(const TIterator&)
{
  /// Overriden operator= (imposed by Root's declaration of TIterator ?)
  Fatal("TIterator::operator=","Not implementeable"); // because there's no clone in TIterator :-(
  return *this;
}

//_____________________________________________________________________________
AliHistogramCollectionIterator::~AliHistogramCollectionIterator()
{
  /// dtor
  Reset();
}

//_____________________________________________________________________________
TObject* AliHistogramCollectionIterator::Next()
{
  /// Advance to next object in the collection
  
  if (!fHashListIterator)
  {
    if ( !fMapIterator ) 
    {
      fMapIterator = fkHistogramCollection->Map()->MakeIterator(fDirection);
    }
    TObjString* key = static_cast<TObjString*>(fMapIterator->Next());
    if (!key)
    {
      // we are done
      return 0x0;
    }      
    THashList* list = static_cast<THashList*>(fkHistogramCollection->Map()->GetValue(key->String().Data()));
    if (!list) return 0x0;
    fHashListIterator = list->MakeIterator(fDirection);
  }

  TObject* o = fHashListIterator->Next();
  
  if (!o) 
  {
    delete fHashListIterator;
    fHashListIterator = 0x0;
    return Next();
  }
  
  return o;
}

//_____________________________________________________________________________
void AliHistogramCollectionIterator::Reset()
{
  /// Reset the iterator
  delete fHashListIterator;
  delete fMapIterator;
}
