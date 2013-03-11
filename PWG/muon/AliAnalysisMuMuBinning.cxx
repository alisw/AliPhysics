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
  AddBin(bin.Particle().Data(),bin.Type().Data(),
         bin.Xmin(),bin.Xmax(),
         bin.Ymin(),bin.Ymax(),bin.Flavour());
}

//______________________________________________________________________________
void AliAnalysisMuMuBinning::AddBin(const char* particle, const char* type,
                                    Double_t xmin, Double_t xmax,
                                    Double_t ymin, Double_t ymax,
                                    const char* flavour)
{
  /// Add a bin
  /// Note that particle and type are not case sensitive.
  if (!fBins)
  {
    fBins = new TMap;
    fBins->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  
  TString sparticle(particle);
  sparticle.ToUpper();
  
  TString stype(type);
  stype.ToUpper();
  
  TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(sparticle.Data()));
  if (!b)
  {
    b = new TObjArray;
    b->SetOwner(kTRUE);
    fBins->Add(new TObjString(sparticle),b);
  }

  Range* r = new Range(sparticle.Data(),stype.Data(),xmin,xmax,ymin,ymax,flavour);
  
  if ( b->FindObject(r) )
  {
    AliDebug(1,"Trying to add an already existing bin. Not doing it.");
    delete r;
  }
  else
  {
    b->Add(r);
    
    TString btype(Form("%s-%s",sparticle.Data(),stype.Data()));
    
    TString name(GetName());
    
    if ( !name.Contains(btype) )
    {
      if (name.Length()>0)
      {
        name += " ";
      }
        
      name += btype;
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

  TIter nextParticle(fBins);
  TObjString* particle;
  
  while ( ( particle = static_cast<TObjString*>(nextParticle()) ) )
  {
    TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(particle->String().Data()));
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
TObjArray* AliAnalysisMuMuBinning::CreateBinObjArray(const char* particle) const
{
  /// Get the list of bins for a given particle 
  /// The returned array must be deleted by the user
  
  if (!fBins) return 0x0;
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  TString sparticle(particle);
  sparticle.ToUpper();
  
  TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(sparticle.Data()));
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
TObjArray* AliAnalysisMuMuBinning::CreateBinObjArray(const char* particle, const char* type, const char* flavour) const
{
  /// Get the list of bins for a given particle and given type
  /// The returned array must be deleted by the user
  /// Type can be a single type or several types separated by comma
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  TString sparticle(particle);
  sparticle.ToUpper();

  TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(sparticle.Data()));
  if (!b) return 0x0;
  
  TIter next(b);
  Range* r;

  TString stype(type);
  stype.ToUpper();

  TString sflavour(flavour);
  
  TObjArray* types = stype.Tokenize(",");
  TObjString* onetype;
  TIter nextType(types);
  
  while ( ( r = static_cast<Range*>(next()) ) )
  {
    nextType.Reset();
    while ( ( onetype = static_cast<TObjString*>(nextType()) ) )
    {
      if ( r->Type() == onetype->String() &&
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
  
  delete types;
  return a;
}

//______________________________________________________________________________
void AliAnalysisMuMuBinning::CreateMesh(const char* particle,
                                        const char* type1, const char* type2,
                                        const char* flavour,
                                        Bool_t remove12)
{
  /// Create 2D bins from existing 1d ones of type1 and type2
  TObjArray* a1 = CreateBinObjArray(particle,type1,flavour);
  if (!a1)
  {
    AliError(Form("No bin for type %s. Done nothing.",type1));
    return;
  }
  TObjArray* a2 = CreateBinObjArray(particle,type2,flavour);
  if (!a2)
  {
    AliError(Form("No bin for type %s. Done nothing.",type2));
    return;
  }
  
  TString meshType(Form("%s VS %s - %s",type1,type2,flavour));
  
  for ( Int_t i1 = 0; i1 <= a1->GetLast(); ++i1 )
  {
    Range* r1 = static_cast<Range*>(a1->At(i1));

    for ( Int_t i2 = 0; i2 <= a2->GetLast(); ++i2 )
    {
      Range* r2 = static_cast<Range*>(a2->At(i2));
      
      AddBin(particle,meshType,r2->Xmin(),r2->Xmax(),r1->Xmin(),r1->Xmax(),Form("%s VS %s",r1->Flavour().Data(),r2->Flavour().Data()));
    }
  }
  
  delete a1;
  delete a2;

  if ( remove12 )
  {
    TObjArray* b = static_cast<TObjArray*>(fBins->GetValue(particle));
    TIter next(b);
    Range* r;
    TString stype1(type1);
    TString stype2(type2);

    stype1.ToUpper();
    stype2.ToUpper();
    
    while ( ( r = static_cast<Range*>(next())) )
    {
      if ( r->Type() == type1 ||
           r->Type() == type2 )
      {
        b->Remove(r);
      }
    }
  }
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateParticleArray() const
{
  /// Create a TObjString array with the names of the particle we're holding results for
  /// Returned array must be delete by user
  
  TObjArray* particleArray(0x0);
  
  TIter nextParticle(fBins);
  TObjString* particle;

  while ( ( particle = static_cast<TObjString*>(nextParticle()) ) )
  {
    if (!particleArray)
    {
      particleArray = new TObjArray;
      particleArray->SetOwner(kTRUE);
    }
    particleArray->Add(new TObjString(*particle));
  }
  return particleArray;
}

//______________________________________________________________________________
TObjArray* AliAnalysisMuMuBinning::CreateTypeArray() const
{
  /// Create a TObjString array with the names of the binning types
  /// Returned array must be delete by user

  TObjArray* typeArray(0x0);
  
  TIter nextParticle(fBins);
  TObjString* particle;
  
  while ( ( particle = static_cast<TObjString*>(nextParticle()) ) )
  {
    TObjArray* particles = static_cast<TObjArray*>(fBins->GetValue(particle->String()));
    
    TIter next(particles);
    Range* r;
    
    while ( ( r = static_cast<Range*>(next())) )
    {
      if (!typeArray)
      {
        typeArray = new TObjArray;
        typeArray->SetOwner(kTRUE);
      }
      if ( !typeArray->FindObject(r->Type()) )
      {
        typeArray->Add(new TObjString(r->Type()));
      }
    }
  }
  return typeArray;
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
AliAnalysisMuMuBinning::Project(const char* particle, const char* type, const char* flavour) const
{
  /// Create a sub-binning object with only the bins pertaining to (particle,type)
  
  TObjArray* bins = CreateBinObjArray(particle,type,flavour);
  if (!bins) return 0x0;
  AliAnalysisMuMuBinning* p = new AliAnalysisMuMuBinning;
  TIter next(bins);
  AliAnalysisMuMuBinning::Range* bin;
  TString stype(type);
  stype.ToUpper();
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    assert  (bin->Type()==stype && bin->Flavour()==flavour);
    {
      p->AddBin(particle,bin->Type(),bin->Xmin(),bin->Xmax(),bin->Ymin(),bin->Ymax(),bin->Flavour().Data());
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
  
  TIter nextParticle(fBins);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(nextParticle()) ) )
  {
    std::cout << "Particle : " << str->String().Data() << std::endl;
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
AliAnalysisMuMuBinning::Range::Range(const char* particle, const char* type,
                                     Double_t xmin, Double_t xmax,
                                     Double_t ymin, Double_t ymax,
                                     const char* flavour)
: TObject(), fParticle(particle), fType(type),
  fXmin(xmin), fXmax(xmax), fYmin(ymin), fYmax(ymax),
  fFlavour(flavour)
{
  /// ctor
  fParticle.ToUpper();
  fType.ToUpper();
}

//______________________________________________________________________________
TString AliAnalysisMuMuBinning::Range::AsString() const
{
  /// Return a string representation of this range
  
  if ( IsNullObject()) return "";
  
  TString s;
  
  if ( fFlavour.Length() > 0 )
  {
    s.Form("%s_%s_%05.2f_%05.2f",Type().Data(),Flavour().Data(),Xmin(),Xmax());
  }
  else
  {
    s.Form("%s_%05.2f_%05.2f",Type().Data(),Xmin(),Xmax());
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
  
  int s = strcmp(Particle().Data(),other->Particle().Data());
  
  if ( s ) return s;
  
  s = strcmp(Type().Data(),other->Type().Data());
  
  if (s) return s;
  
  s = strcmp(Flavour().Data(),other->Flavour().Data());
  
  if (s) return s;
  
  if ( IsNullObject() ) return 0;
  
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
  
  if ( IsNullObject() )
  {
    return kTRUE;
  }
  else
  {
    if ( Is2D() )
    {
      return ( x > Xmin() && x < Xmax() && y > Ymin() && y < Ymax());
    }
    else
    {
      return ( x > Xmin() && x < Xmax() );
    }
  }
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuBinning::Range::IsNullObject() const
{
  /// Whether we're a null object or not
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
  
  if (IsNullObject())
  {
    std::cout << Form("%s : %s : NullObject",Particle().Data(),Type().Data()) << std::endl;
  
    return;
  }
  
  std::cout << Form("%s : %s : %5.2f : %5.2f",Particle().Data(),Type().Data(),Xmin(),Xmax());
  
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


