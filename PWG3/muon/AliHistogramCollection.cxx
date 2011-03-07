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

#include "AliLog.h"
#include "Riostream.h"
#include "TError.h"
#include "TH1.h"
#include "THashList.h"
#include "TKey.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TROOT.h"
#include "TSystem.h"

ClassImp(AliHistogramCollection)

//_____________________________________________________________________________
AliHistogramCollection::AliHistogramCollection(const char* name, const char* title) 
: TNamed(name,title), fMap(0x0), fMustShowEmptyHistogram(kFALSE)
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
AliHistogramCollection::Adopt(const char* key, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA)
  return InternalAdopt(Form("/%s/./././",key),histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* keyA, const char* keyB, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA,keyB)
  return InternalAdopt(Form("/%s/%s/././",keyA,keyB),histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* keyA, const char* keyB, const char* keyC, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA,keyB,keyC)
  return InternalAdopt(Form("/%s/%s/%s/./",keyA,keyB,keyC),histo);
}

//_____________________________________________________________________________
Bool_t 
AliHistogramCollection::Adopt(const char* keyA, const char* keyB, const char* keyC, const char* keyD, TH1* histo)
{
  /// Adopt a given histogram, and associate it with pair (keyA,keyB,keyC,keyD)
  return InternalAdopt(Form("/%s/%s/%s/%s/",keyA,keyB,keyC,keyD),histo);
}

//_____________________________________________________________________________
TIterator*
AliHistogramCollection::CreateIterator(Bool_t direction) const
{
  /// Create an iterator (must be deleted by the client)
  return fMap ? new AliHistogramCollectionIterator(this,direction) : 0x0;
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
TH1* 
AliHistogramCollection::Histo(const char* keyA, 
                              const char* histoname) const
{
  /// Get histo for (keyA,histoname) triplet
  
  return InternalHisto(Form("/%s/./././",keyA),histoname);
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* keyA, const char* keyB, 
                              const char* histoname) const
{
  /// Get histo for (keyA,keyB,histoname) triplet

  return InternalHisto(Form("/%s/%s/././",keyA,keyB),histoname);
}

//_____________________________________________________________________________
TH1* 
AliHistogramCollection::Histo(const char* keyA, const char* keyB, const char* keyC,
                              const char* histoname) const
{
  /// Get histo for (keyA,keyB,keyC,histoname) quad
  
  return InternalHisto(Form("/%s/%s/%s/./",keyA,keyB,keyC),histoname);
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
  
  return InternalDecode(identifier,4);  
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
  
  if ( !fMap ) 
  {
    fMap = new TMap;
    fMap->SetOwner(kTRUE);
  }
  
  hlist = static_cast<THashList*>(fMap->GetValue(identifier));
  
  if (!hlist)
  {
    hlist = new THashList;
    hlist->SetOwner(kTRUE);
    fMap->Add(new TObjString(identifier),hlist);
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
AliHistogramCollection::Histo(const char* identifier) const
{
  /// Get histogram keyA/keyB/keyC/keyD/histoname
  return Histo(InternalDecode(identifier,0),
               InternalDecode(identifier,1),
               InternalDecode(identifier,2),
               InternalDecode(identifier,3),
               InternalDecode(identifier,4));
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
  /// histo is index=4
  
  if ( identifier[0] != '/' ) 
  {    
    AliError(Form("identifier %s is malformed.",identifier));
    return "";
  }
  
  TObjArray* array = TString(identifier).Tokenize("/");

  if ( array->GetLast()>5 ) 
  {
    AliError(Form("identifier %s is malformed.",identifier));
    delete array;
    return "";    
  }

  TString value("");
  
  if ( index <= array->GetLast() ) 
  {
    value = static_cast<TObjString*>(array->At(index))->String();
  }
  
  delete array;
  
  return value;
  
// "Custom" implementation below is even slower that Tokenize... ??
// or at least did not change the timing result enough to be worthwhile (indicating
// the cpu time is wasted elsewhere ?)
//
//  Int_t slashPos[6] = {0};
//  Int_t nslashes(0);
//  
//  for ( Int_t i = 0; i < identifier.Length(); ++i ) 
//  {
//    if ( identifier[i] == '/' ) 
//    {
//      slashPos[nslashes++]=i;      
//    }
//    if ( nslashes > 5 ) 
//    {
//      AliError(Form("identifier %s is malformed.",identifier.Data()));
//      return "";
//    }
//  }
//  
//  slashPos[nslashes++]=identifier.Length();
//  --nslashes;
//
//  if ( index < nslashes )
//  {
//    return TString(identifier(slashPos[index]+1,slashPos[index+1]-slashPos[index]-1));
//  }
//  return "";
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
  
  THashList* hlist = static_cast<THashList*>(fMap->GetValue(identifier));
  if (!hlist) 
  {
    return 0x0;
  }
  
  return static_cast<TH1*>(hlist->FindObject(histoname));  
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
    
    TIter nextIdentifier(hcol->fMap);
    TObjString* identifier;

    while ( ( identifier = static_cast<TObjString*>(nextIdentifier()) ) )
    {
      THashList* otherList = static_cast<THashList*>(hcol->fMap->GetValue(identifier->String().Data()));

      TIter nextHisto(otherList);
      TH1* h;
      
      while ( ( h = static_cast<TH1*>(nextHisto()) ) )
      {
        TH1* thisHisto = Histo(KeyA(identifier->String()).Data(),
                               KeyB(identifier->String()).Data(),
                               KeyC(identifier->String()).Data(),
                               KeyD(identifier->String()).Data(),
                               h->GetName());
        
        if (!thisHisto)
        {
          AliDebug(1,Form("Adopting a new histo = %s/%s",identifier->String().Data(),h->GetName()));
          
          // this is an histogram we don't have yet. Let's add it
          Adopt(KeyA(identifier->String()),
                KeyB(identifier->String()),
                KeyC(identifier->String()),
                KeyD(identifier->String()),
                static_cast<TH1*>(h->Clone()));
        }
        else
        {
          // add it...
          AliDebug(1,Form("Merging histo = %s/%s (%g vs %g)",
                          identifier->String().Data(),
                          h->GetName(),
                          h->GetSumOfWeights(),
                          thisHisto->GetSumOfWeights()));
          TList l;
          l.Add(h);
          
          thisHisto->Merge(&l);
        }
      }
    }
  }
         
  AliDebug(1,Form("count=%d",count));
  
  return count+1;
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
  
  cout << Form("AliHistogramCollection : %d keys and %d histos",
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
    
    if ( InternalDecode(identifier,0).Contains(reKeyA) &&
        InternalDecode(identifier,1).Contains(reKeyB) &&
        InternalDecode(identifier,2).Contains(reKeyC) &&
        InternalDecode(identifier,3).Contains(reKeyD)
        )
    {
      if ( sreHistoname == "*" ) 
      {
        identifierPrinted = kTRUE;
        cout << identifier.Data() << endl;
      }
      
      THashList * list = static_cast<THashList*>(fMap->GetValue(sid->String().Data()));
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
UInt_t 
AliHistogramCollection::EstimateSize(Bool_t show) const
{
  /// Estimate the memory (in kilobytes) used by our histograms
  
//  sizeof(TH1) + (nbins+2)*(nbytes_per_bin) +name+title_sizes 
//  if you have errors add (nbins+2)*8 
    
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
    
    if ( hasErrors) thissize += nbins*8;

    n += thissize;
    
    if ( show ) 
    {
      AliInfo(Form("Size of %30s is %20d bytes",h->GetName(),thissize));
    }
  }

  return n;
}

//_____________________________________________________________________________
void AliHistogramCollection::PruneEmptyHistograms()
{
  /// Delete the empty histograms
  TIter next(fMap);
  TObjString* key;
  
  TList toBeRemoved;
  toBeRemoved.SetOwner(kTRUE);
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    TString identifier(key->String());
    THashList* hlist = static_cast<THashList*>(fMap->GetValue(identifier.Data()));
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
AliHistogramCollection::Project(const char* /*keyA*/, const char* /*keyB*/) const
{
  /// To be implemented : would create a new collection starting at keyA/keyB
  AliError("Implement me !");
  return 0x0;
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
    
  TString skey = Form("/%s/%s/%s/%s/",
                      KeyA(identifier).Data(),
                      KeyB(identifier).Data(),
                      KeyC(identifier).Data(),
                      KeyD(identifier).Data());
  
  THashList* hlist = dynamic_cast<THashList*>(fMap->GetValue(skey.Data()));
  
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

//  if ( hlist->IsEmpty() ) 
//  {
//    // we should remove the key as well
//    TObject* k = fMap->Remove(key);
//    if (!k)
//    {
//      AliError("Removal of the key failed");
//    }
//  }
  
  return o;
}

//_____________________________________________________________________________
TObjArray*
AliHistogramCollection::SortAllIdentifiers() const
{
  /// Sort our internal identifiers. Returned array must be deleted.
  TObjArray* identifiers = new TObjArray;
  identifiers->SetOwner(kFALSE); 
  TIter next(fMap);
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
      fMapIterator = fkHistogramCollection->fMap->MakeIterator(fDirection);
    }
    TObjString* key = static_cast<TObjString*>(fMapIterator->Next());
    if (!key)
    {
      // we are done
      return 0x0;
    }      
    THashList* list = static_cast<THashList*>(fkHistogramCollection->fMap->GetValue(key->String().Data()));
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
