#include "AliAnalysisMuMuBinning.h"

//
// AliAnalysisMuMuBinning.h : a class to hold bins of various sizes
//
// The idea behind this class is to store it together with
// the histograms (corresponding to the bins)
// so we can then loop easily on all bins afterwards.
// 
// author: L. Aphecetche (Subatech)
//

#include "AliLog.h"
#include "Riostream.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <cassert>

ClassImp(AliAnalysisMuMuBinning::Range)
ClassImp(AliAnalysisMuMuBinning)

//______________________________________________________________________________
AliAnalysisMuMuBinning::AliAnalysisMuMuBinning(const char* name, const char* title)
: TNamed(name,title), fBins(0x0)
{
  // default ctor
}

//______________________________________________________________________________
AliAnalysisMuMuBinning::AliAnalysisMuMuBinning(const AliAnalysisMuMuBinning& rhs)
: TNamed(), fBins(0x0)
{
  // copy ctor
  TObjArray* bins = rhs.CreateBinObjArray();
  TIter next(bins);
  AliAnalysisMuMuBinning::Range* b;
  
  while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    AddBin(*b);
  }
}

//______________________________________________________________________________
AliAnalysisMuMuBinning& AliAnalysisMuMuBinning::operator=(const AliAnalysisMuMuBinning& rhs)
{
  // assignment  operator
  
  if ( this != &rhs )
  {
    delete fBins;
    fBins = 0x0;
    TObjArray* bins = rhs.CreateBinObjArray();
    TIter next(bins);
    AliAnalysisMuMuBinning::Range* b;
  
    while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
    {
      AddBin(*b);
    }
  }
  return *this;
}

//______________________________________________________________________________
AliAnalysisMuMuBinning::~AliAnalysisMuMuBinning()
{
  // dtor
  delete fBins;
}

//______________________________________________________________________________
void AliAnalysisMuMuBinning::AddBin(const AliAnalysisMuMuBinning::Range& bin)
{
  /// add one bin
  AddBin(bin.What().Data(),bin.Quantity().Data(),
         bin.Xmin(),bin.Xmax(),
         bin.Ymin(),bin.Ymax(),bin.Flavour());
}

//______________________________________________________________________________
void AliAnalysisMuMuBinning::AddBin(const char* what, const char* quantity,
                                    Double_t xmin, Double_t xmax,
                                    Double_t ymin, Double_t ymax,
                                    const char* flavour)
{
  /// Add a bin
  /// Note that what and Quantity are not case sensitive.
  if (!fBins)
  {
    fBins = new TMap;
    fBins->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  
  TString swhat(what);
  swhat.ToUpper();
  
  TString sQuantity(quantity);
  sQuantity.ToUpper();
  
  TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(swhat.Data()));
  if (!b)
  {
    b = new TObjArray;
    b->SetOwner(kTRUE);
    fBins->Add(new TObjString(swhat),b);
  }

  Range* r = new Range(swhat.Data(),sQuantity.Data(),xmin,xmax,ymin,ymax,flavour);
  
  if ( b->FindObject(r) )
  {
    AliDebug(1,Form("Trying to add an already existing bin : %s. Not doing it.",r->AsString().Data()));
    delete r;
  }
  else
  {
    b->Add(r);
    
    TString bQuantity(Form("%s-%s",swhat.Data(),sQuantity.Data()));
    
    TString name(GetName());
    
    if ( !name.Contains(bQuantity) )
    {
      if (name.Length()>0)
      {
        name += " ";
      }
        
      name += bQuantity;
      SetName(name);
    }
  }
  
  b->Sort();
}

//______________________________________________________________________________
Double_t* AliAnalysisMuMuBinning::CreateBinArray() const
{
  /// Create a (variable) bin array suitable for TH1 constructor
  /// The returned array must be deleted by the user
  /// (using delete[] )
  
  TObjArray* binArray = CreateBinObjArray();
  if (!binArray) return 0x0;
  Double_t* bins = new Double_t[binArray->GetEntries()+1];
  TIter next(binArray);
  AliAnalysisMuMuBinning::Range* b;
  Int_t i(0);
  while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    bins[i] = b->Xmin();
    ++i;
  }

  b = static_cast<AliAnalysisMuMuBinning::Range*>(binArray->At(binArray->GetEntries()-1));
  
  bins[i] = b->Xmax();
  
  delete binArray;
  return bins;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateBinObjArray() const
{
  /// Get the list of all the bins
  /// The returned array must be deleted by the user
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);

  TIter nextwhat(fBins);
  TObjString* what;
  
  while ( ( what = static_cast<TObjString*>(nextwhat()) ) )
  {
    TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(what->String().Data()));
    TIter next(b);
    Range* r;
  
    while ( ( r = static_cast<Range*>(next()) ) )
    {
      a->Add(r->Clone());
    }
  }
  
  if ( a->GetLast() < 0 )
  {
    delete a;
    a = 0x0;
  }
  return a;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateBinObjArray(const char* what) const
{
  /// Get the list of bins for a given what 
  /// The returned array must be deleted by the user
  
  if (!fBins) return 0x0;
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  TString swhat(what);
  swhat.ToUpper();
  
  TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(swhat.Data()));
  if (!b) return 0x0;
  
  TIter next(b);
  Range* r;
  
  while ( ( r = static_cast<Range*>(next()) ) )
  {
    a->Add(r->Clone());
  }
  
  if ( a->GetLast() < 0 )
  {
    delete a;
    a = 0x0;
  }
  return a;
}


//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateBinObjArray(const char* what, const char* quantity, const char* flavour) const
{
  /// Get the list of bins for a given what and given Quantity
  /// The returned array must be deleted by the user
  /// Quantity can be a single Quantity or several Quantities separated by comma
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  TString swhat(what);
  swhat.ToUpper();

  TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(swhat.Data()));
  if (!b) return 0x0;
  
  TIter next(b);
  Range* r;

  TString sQuantity(quantity);
  sQuantity.ToUpper();

  TString sflavour(flavour);
  
  TObjArray* Quantitys = sQuantity.Tokenize(",");
  TObjString* oneQuantity;
  TIter nextQuantity(Quantitys);
  
  while ( ( r = static_cast<Range*>(next()) ) )
  {
    nextQuantity.Reset();
    while ( ( oneQuantity = static_cast<TObjString*>(nextQuantity()) ) )
    {
      if ( r->Quantity() == oneQuantity->String() &&
          ( ( sflavour.Length() > 0 && r->Flavour() == sflavour.Data() ) || sflavour.Length()==0 ) )
      {
        a->Add(r->Clone());
      }
    }
  }
  
  if ( a->GetLast() < 0 )
  {
    delete a;
    a = 0x0;
  }
  
  delete Quantitys;
  return a;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateBinStrArray() const
{
  /// Get the list of all the bin names
  /// The returned array must be deleted by the user

  TObjArray* a = CreateBinObjArray();
  TObjArray* s = new TObjArray;
  s->SetOwner(kTRUE);
  
  Range* r;
  
  TIter next(a);
  while ( ( r = static_cast<Range*>(next())))
  {
    s->Add(new TObjString(r->AsString()));
  }
  
  delete a;
  
  return s;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateBinStrArray(const char* what) const
{
  /// Get the list of bin names for a given what
  /// The returned array must be deleted by the user

  TObjArray* a = CreateBinObjArray(what);
  TObjArray* s = new TObjArray;
  s->SetOwner(kTRUE);
  
  Range* r;
  
  TIter next(a);
  while ( ( r = static_cast<Range*>(next())))
  {
    s->Add(new TObjString(r->AsString()));
  }
  
  delete a;
  
  return s;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateBinStrArray(const char* what, const char* quantity, const char* flavour) const
{
  /// Get the list of bin names, for a given what and given quantity
  /// The returned array must be deleted by the user
  /// Quantity can be a single Quantity or several Quantities separated by comma

  TObjArray* a = CreateBinObjArray(what,quantity,flavour);
  TObjArray* s = new TObjArray;
  s->SetOwner(kTRUE);
  
  Range* r;
  
  TIter next(a);
  while ( ( r = static_cast<Range*>(next())))
  {
    s->Add(new TObjString(r->AsString()));
  }
  
  delete a;
  
  return s;
}

//______________________________________________________________________________
Double_t* AliAnalysisMuMuBinning::CreateBinArrayY() const
{
  /// Create a TObjArray with 2 (variable) bin array with x and y binning, suitable for TH1
  /// The returned array must be deleted by the user
  /// (using delete[] )
  
  TObjArray* binArray = CreateBinObjArray();
  if (!binArray) return 0x0;
  
  AliAnalysisMuMuBinning::Range* firstBin = static_cast<AliAnalysisMuMuBinning::Range*>(binArray->First());
  
  Double_t* binsY = new Double_t[binArray->GetEntries()+1];
  
  Double_t minY = firstBin->Ymin();
  Double_t maxY = firstBin->Ymax();
  
  if ( !(minY < maxY))
  {
    std::cout << "No 2D binning" << std::endl;
    delete[] binsY;
    delete binArray;
    return 0x0;
  }
  
  TIter next(binArray);
  AliAnalysisMuMuBinning::Range* b;
  Int_t i(0);
  while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if ( (i != 0) && (minY == b->Ymin()) ) break;
    
    binsY[i] = b->Ymin();
    ++i;
  }
  
  b = static_cast<AliAnalysisMuMuBinning::Range*>(binArray->At(i-1));
  
  binsY[i] = b->Ymax();
  
  Double_t* bins = new Double_t[i + 1];
  
  for (Int_t j = 0 ; j < (i+1) ; j++)
  {
    bins[j] = binsY[j];
  }
  
  delete binArray;
  delete[] binsY;
  return bins;
}

//______________________________________________________________________________
Int_t AliAnalysisMuMuBinning::GetNBinsX() const
{
  /// Gets the number of x bins, suitable for TH1
  
  
  TObjArray* binArray = CreateBinObjArray();
  if (!binArray) return 0;
  
  Double_t binsX(0);
  
  
  TIter next(binArray);
  AliAnalysisMuMuBinning::Range* b;
  Int_t i(0);
  while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if ( (i != 0) && (binsX == b->Xmin()) ) continue;
    
    binsX = b->Xmin();
    ++i;
  }
  
  delete binArray;
  return i;
}

//______________________________________________________________________________
Int_t AliAnalysisMuMuBinning::GetNBinsY() const
{
  /// Gets the number of x bins, suitable for TH1
  
  
  TObjArray* binArray = CreateBinObjArray();
  if (!binArray) return 0x0;
  
  Double_t binsY(0);
  
  AliAnalysisMuMuBinning::Range* firstBin = static_cast<AliAnalysisMuMuBinning::Range*>(binArray->First());
  
  Double_t minY = firstBin->Ymin();
  Double_t maxY = firstBin->Ymax();
  
  if ( !(minY < maxY))
  {
    std::cout << "No 2D binning" << std::endl;
    delete binArray;
    return 0;
  }
  
  TIter next(binArray);
  AliAnalysisMuMuBinning::Range* b;
  Int_t i(0);
  while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if ( (i != 0) && (minY == b->Ymin()) ) break;
    
    binsY = b->Ymin();
    ++i;
  }
  
  delete binArray;
  return i;
}

//______________________________________________________________________________
Double_t* AliAnalysisMuMuBinning::CreateBinArrayX() const
{
  /// Create a TObjArray with 2 (variable) bin array with x and y binning , suitable for TH1
  /// The returned array must be deleted by the user
  /// (using delete[] )
  
  TObjArray* binArray = CreateBinObjArray();
  if (!binArray) return 0x0;
  
  Double_t* binsX = new Double_t[binArray->GetEntries()+1];
  
  
  TIter next(binArray);
  AliAnalysisMuMuBinning::Range* b;
  Int_t i(0);
  while ( ( b = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if ( (i != 0) && (binsX[i-1] == b->Xmin()) ) continue;
    
    binsX[i] = b->Xmin();
    ++i;
  }
  
  b = static_cast<AliAnalysisMuMuBinning::Range*>(binArray->At(binArray->GetEntries()-1));
  
  binsX[i] = b->Xmax();
  
  Double_t* bins = new Double_t[i + 1];
  
  for (Int_t j = 0 ; j < (i+1) ; j++)
  {
    bins[j] = binsX[j];
  }
  
  delete binArray;
  delete[] binsX;
  return bins;
}

//______________________________________________________________________________
void AliAnalysisMuMuBinning::CreateMesh(const char* what,
                                        const char* quantity1, const char* quantity2,
                                        const char* flavour,
                                        Bool_t remove12)
{
  /// Create 2D bins from existing 1d ones of Quantity1 and Quantity2
  TObjArray* a1 = CreateBinObjArray(what,quantity1,flavour);
  if (!a1)
  {
    AliError(Form("No bin for Quantity %s. Done nothing.",quantity1));
    return;
  }
  TObjArray* a2 = CreateBinObjArray(what,quantity2,flavour);
  if (!a2)
  {
    AliError(Form("No bin for Quantity %s. Done nothing.",quantity2));
    return;
  }
  
  TString meshQuantity(Form("%s VS %s - %s",quantity1,quantity2,flavour));
  
  for ( Int_t i1 = 0; i1 <= a1->GetLast(); ++i1 )
  {
    Range* r1 = static_cast<Range*>(a1->At(i1));

    for ( Int_t i2 = 0; i2 <= a2->GetLast(); ++i2 )
    {
      Range* r2 = static_cast<Range*>(a2->At(i2));
      
      AddBin(what,meshQuantity,r2->Xmin(),r2->Xmax(),r1->Xmin(),r1->Xmax(),Form("%s VS %s",r1->Flavour().Data(),r2->Flavour().Data()));
    }
  }
  
  delete a1;
  delete a2;

  if ( remove12 )
  {
    TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(what));
    TIter next(b);
    Range* r;
    TString sQuantity1(quantity1);
    TString sQuantity2(quantity2);

    sQuantity1.ToUpper();
    sQuantity2.ToUpper();
    
    while ( ( r = static_cast<Range*>(next())) )
    {
      if ( r->Quantity() == quantity1 ||
           r->Quantity() == quantity2 )
      {
        b->Remove(r);
      }
    }
  }
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateWhatArray() const
{
  /// Create a TObjString array with the names of the what we're holding results for
  /// Returned array must be delete by user
  
  TObjArray* whatArray(0x0);
  
  TIter nextwhat(fBins);
  TObjString* what;

  while ( ( what = static_cast<TObjString*>(nextwhat()) ) )
  {
    if (!whatArray)
    {
      whatArray = new TObjArray;
      whatArray->SetOwner(kTRUE);
    }
    whatArray->Add(new TObjString(*what));
  }
  return whatArray;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateQuantityArray() const
{
  /// Create a TObjString array with the names of the binning Quantitys
  /// Returned array must be delete by user

  TObjArray* QuantityArray(0x0);
  
  TIter nextwhat(fBins);
  TObjString* what;
  
  while ( ( what = static_cast<TObjString*>(nextwhat()) ) )
  {
    TObjArray* whats = static_cast<TObjArray*>(fBins->GetValue(what->String()));
    
    TIter next(whats);
    Range* r;
    
    while ( ( r = static_cast<Range*>(next())) )
    {
      if (!QuantityArray)
      {
        QuantityArray = new TObjArray;
        QuantityArray->SetOwner(kTRUE);
      }
      if ( !QuantityArray->FindObject(r->Quantity()) )
      {
        QuantityArray->Add(new TObjString(r->Quantity()));
      }
    }
  }
  return QuantityArray;
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuBinning::IsEqual(const TObject* obj) const
{
  /// Return kTRUE if obj is an AliAnalysisMuMuBinning object and is
  /// equal to *this
  
  if ( obj->IsA() == AliAnalysisMuMuBinning::Class() )
  {
    const AliAnalysisMuMuBinning* other = static_cast<const AliAnalysisMuMuBinning*>(obj);
    
    TIter nextOther(other->fBins);
    TObjString* str;
    
    while ( ( str = static_cast<TObjString*>(nextOther()) ) )
    {
      TObject* o = fBins->GetValue(str->String());
      if (!o) return kFALSE;
      if (o->IsA() != TObjArray::Class()) return kFALSE;

      TObjArray* thisArray = static_cast<TObjArray*>(o);

      o = other->fBins->GetValue(str->String());
      if (!o) return kFALSE;
      if (o->IsA() != TObjArray::Class()) return kFALSE;
      
      TObjArray* otherArray = static_cast<TObjArray*>(o);
      
      Int_t n = thisArray->GetEntries();
      
      if ( n != otherArray->GetEntries() ) return kFALSE;
      
      for ( Int_t i = 0; i < n; ++i )
      {
        Range* thisRange = static_cast<Range*>(thisArray->At(i));
        Range* otherRange = static_cast<Range*>(otherArray->At(i));
        
        if ( !thisRange->IsEqual(otherRange) ) return kFALSE;
      }
    }
    return kTRUE;
  }
  
  return kFALSE;
}


//______________________________________________________________________________
Long64_t AliAnalysisMuMuBinning::Merge(TCollection* list)
{
  /// Merge method
  
  // Merge a list of AliAnalysisMuMuBinning objects with this
  // Returns the number of merged objects (including this).
  
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  TIter next(list);
  TObject* currObj;
  Int_t count(0);
  
  while ( ( currObj = next() ) )
  {
    AliAnalysisMuMuBinning* binning = dynamic_cast<AliAnalysisMuMuBinning*>(currObj);
    if (!binning)
    {
      AliFatal(Form("object named \"%s\" is a %s instead of an AliAnalysisMuMuBinning!", currObj->GetName(), currObj->ClassName()));
      continue;
    }
    
    if ( IsEqual(binning) )
    {
      // nothing to do if we have the same binning already ;-)
    }
    else
    {
      AliWarning("Implement me!");
      std::cout << ">>>> this=" << std::endl;
      Print();
      std::cout << ">>>> other=" << std::endl;
      binning->Print();
    }
    ++count;
  }
  
  return count+1;
}

//______________________________________________________________________________
AliAnalysisMuMuBinning*
AliAnalysisMuMuBinning::Project(const char* what, const char* quantity, const char* flavour) const
{
  /// Create a sub-binning object with only the bins pertaining to (what,Quantity)
  
  TObjArray* bins = CreateBinObjArray(what,quantity,flavour);
  if (!bins) return 0x0;
  AliAnalysisMuMuBinning* p = new AliAnalysisMuMuBinning;
  TIter next(bins);
  AliAnalysisMuMuBinning::Range* bin;
  TString sQuantity(quantity);
  sQuantity.ToUpper();
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    if  (bin->Quantity()!=sQuantity || bin->Flavour()!=flavour)
    {
      AliDebug(1,Form("sQuantity=%s flavour=%s bin=%s => skip",sQuantity.Data(),flavour,bin->AsString().Data()));
      continue;
    }
    {
      p->AddBin(what,bin->Quantity(),bin->Xmin(),bin->Xmax(),bin->Ymin(),bin->Ymax(),bin->Flavour().Data());
    }
  }
  
  delete bins;
  
  return p;
}

//______________________________________________________________________________
void AliAnalysisMuMuBinning::Print(Option_t* /*opt*/) const
{
  /// Print the bins
  
  if (!fBins)
  {
    std::cout << "Empty object. No bin defined. " << std::endl;
    return;
  }
  
  TIter nextwhat(fBins);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(nextwhat()) ) )
  {
    std::cout << "what : " << str->String().Data() << std::endl;
    TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(str->String()));
    TIter next(b);
    Range* r(0x0);
    next.Reset();
    Int_t i(1);
    while ( ( r = static_cast<Range*>(next()) ) )
    {
      std::cout << Form("BIN %3d ",i++);
      r->Print();
    }
  }
}

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

//______________________________________________________________________________
AliAnalysisMuMuBinning::Range::Range(const char* what, const char* quantity,
                                     Double_t xmin, Double_t xmax,
                                     Double_t ymin, Double_t ymax,
                                     const char* flavour)
: TObject(), fWhat(what), fQuantity(quantity),
  fXmin(xmin), fXmax(xmax), fYmin(ymin), fYmax(ymax),
  fFlavour(flavour)
{
  /// ctor
  fWhat.ToUpper();
  fQuantity.ToUpper();
}

//______________________________________________________________________________
TString AliAnalysisMuMuBinning::Range::AsString() const
{
  /// Return a string representation of this range
  
  if ( IsIntegrated()) return Quantity().Data();
  
  TString s;
  
  if ( fFlavour.Length() > 0 )
  {
    s.Form("%s_%s_%05.2f_%05.2f",Quantity().Data(),Flavour().Data(),Xmin(),Xmax());
  }
  else
  {
    s.Form("%s_%05.2f_%05.2f",Quantity().Data(),Xmin(),Xmax());
  }
  
  if (Is2D())
  {
    s += TString::Format("_%06.2f_%06.2f",Ymin(),Ymax());
  }

  s.ReplaceAll(" ","");
  s.ReplaceAll("+","");
  s.ReplaceAll("-","m");

  return s;
}

//______________________________________________________________________________
Int_t	AliAnalysisMuMuBinning::Range::Compare(const TObject* obj) const
{
  //  Must return -1 if this is smaller
  // than obj, 0 if objects are equal and 1 if this is larger than obj.
  const Range* other = static_cast<const Range*>(obj);
  
  int s = strcmp(What().Data(),other->What().Data());
  
  if ( s ) return s;
  
  s = strcmp(Quantity().Data(),other->Quantity().Data());
  
  if (s) return s;
  
  s = strcmp(Flavour().Data(),other->Flavour().Data());
  
  if (s) return s;
  
//  if ( IsIntegrated() && other->IsIntegrated() ) return 0;
  
  if ( Xmin() < other->Xmin() )
  {
    return -1;
  }
  else if ( Xmin() > other->Xmin() )
  {
    return 1;
  }
  else
  {
    if  ( Xmax() < other->Xmax() )
    {
      return -1;
    }
    else if ( Xmax() > other->Xmax() )
    {
      return 1;
    }
    else {
      if ( Ymin() < other->Ymin() )
      {
        return -1;
      }
      else if ( Ymin() > other->Ymin() )
      {
        return 1;
      }
      else
      {
        if ( Ymax() < other->Ymax() )
        {
          return -1;
        }
        else if ( Ymax() > other->Ymax() )
        {
          return 1;
        }
      }
    }
  }
  return 0;
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuBinning::Range::IsInRange(Double_t x, Double_t y) const
{
  /// If Range is 1D, returns true if x is in range
  /// If Range is 2D, returns true if (x,y) is in range
  
  if ( IsIntegrated() )
  {
    return kTRUE;
  }
  else
  {
    if ( Is2D() )
    {
      return ( x >= Xmin() && x < Xmax() && y >= Ymin() && y < Ymax());
    }
    else
    {
      return ( x >= Xmin() && x < Xmax() );
    }
  }
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuBinning::Range::IsIntegrated() const
{
  /// Whether we're a null object (aka integrated) or not
  return
  Xmin() >= TMath::Limits<Double_t>::Max() &&
  Ymin() >= TMath::Limits<Double_t>::Max() &&
  Xmax() >= TMath::Limits<Double_t>::Max() &&
  Ymax() >= TMath::Limits<Double_t>::Max() ;
}


//______________________________________________________________________________
void AliAnalysisMuMuBinning::Range::Print(Option_t* /*opt*/) const
{
  /// Output to stdout
  
  if (IsIntegrated())
  {
    std::cout << Form("%s : %s : INTEGRATED",What().Data(),Quantity().Data()) << std::endl;
  
    return;
  }
  
  std::cout << Form("%s : %s : %5.2f : %5.2f",What().Data(),Quantity().Data(),Xmin(),Xmax());
  
  if (Is2D())
  {
    std::cout << Form(" ; %5.2f : %5.2f",Ymin(),Ymax());
  }
  
  if (Flavour().Length()>0)
  {
    std::cout << " - " << Flavour().Data();
  }
  
  std::cout << "->" << AsString().Data() << std::endl;
  
}


