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

///
/// Base class to hold a set of results for the same quantity,
/// computed using various methods, each with their errors
///
/// author: Laurent Aphecetche (Subatech)
///

#include "AliAnalysisMuMuResult.h"

ClassImp(AliAnalysisMuMuResult)

#include "THashList.h"
#include "TLine.h"
#include "TList.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TParameter.h"
#include "AliLog.h"
#include <map>

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(const char* name, const char* title) :
TNamed(name,title),
fSubResults(0x0),
fMap(0x0),
fMother(0x0),
fKeys(0x0),
fWeight(0.0),
fAlias(),
fSubResultsToBeIncluded(0x0)
{
  /// default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::AliAnalysisMuMuResult(const AliAnalysisMuMuResult& rhs)
:
TNamed(rhs),
fSubResults(0x0),
fMap(0x0),
fMother(0x0),
fKeys(0x0),
fWeight(rhs.fWeight),
fAlias(),
fSubResultsToBeIncluded(0x0)
{
  /// copy ctor
  /// Note that the mother is lost
  /// fKeys remains 0x0 so it will be recomputed if need be

  if (rhs.fSubResults)
  {
    fSubResults = static_cast<TObjArray*>(rhs.fSubResults->Clone());
  }
  
  if ( rhs.fMap )
  {
    fMap = static_cast<TMap*>(rhs.fMap->Clone());
  }
  
  if ( rhs.fAlias.Length() > 0 )
  {
    fAlias = rhs.fAlias;
  }
  
  if ( rhs.fSubResultsToBeIncluded )
  {
    fSubResultsToBeIncluded = static_cast<TList*>(rhs.fSubResultsToBeIncluded->Clone());
  }

}

//_____________________________________________________________________________
AliAnalysisMuMuResult& AliAnalysisMuMuResult::operator=(const AliAnalysisMuMuResult& rhs)
{
  /// Assignment operator
  
  if (this!=&rhs)
  {
    delete fMap;
    delete fSubResults;
    delete fSubResultsToBeIncluded;
    
    fMap = 0x0;
    fSubResults = 0x0;
    fKeys = 0x0;
    fSubResultsToBeIncluded = 0x0;
    
    if (rhs.fSubResults)
    {
      fSubResults = static_cast<TObjArray*>(rhs.fSubResults->Clone());
    }
    
    if ( rhs.fMap )
    {
      fMap = static_cast<TMap*>(rhs.fMap->Clone());
    }

    if ( rhs.fSubResultsToBeIncluded )
    {
      fSubResultsToBeIncluded = static_cast<TList*>(rhs.fSubResultsToBeIncluded->Clone());
    }

    static_cast<TNamed&>(*this)=rhs;

    fWeight = rhs.fWeight;
    fAlias="";
    
    if ( rhs.fAlias.Length() > 0 )
    {
      fAlias = rhs.fAlias;
    }
  }
  
  return *this;
}

//_____________________________________________________________________________
AliAnalysisMuMuResult::~AliAnalysisMuMuResult()
{
  // dtor
  delete fMap;
  delete fSubResults;
  delete fKeys;
  delete fSubResultsToBeIncluded;
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::AdoptSubResult(AliAnalysisMuMuResult* r)
{
  /// Adopt (i.e. becomes owner) of a subresult
  if (!fSubResults)
  {
    fSubResults = new TObjArray;
    fSubResults->SetOwner(kTRUE);
  }

  fSubResults->Add(r);
  
  SubResultsToBeIncluded()->Add(new TObjString(r->Alias()));
}

//_____________________________________________________________________________
TObject* AliAnalysisMuMuResult::Clone(const char* /*newname*/) const
{
  /// Clone this result
  return new AliAnalysisMuMuResult(*this);
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::ErrorAB(Double_t a, Double_t aerr, Double_t b, Double_t berr)
{
  /// Compute the quadratic sum of 2 errors

  Double_t e(0.0);
  
  if ( TMath::Abs(a) > 1E-12 )
  {
    e += (aerr*aerr)/(a*a);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (berr*berr)/(b*b);
  }
  
  return TMath::Sqrt(e);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::ErrorABC(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror)
{
  /// Compute the quadratic sum of 3 errors
  
  Double_t e(0.0);
  
  if ( TMath::Abs(a) > 1E-12 )
  {
    e += (aerr*aerr)/(a*a);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (berr*berr)/(b*b);
  }
  
  if ( TMath::Abs(c) > 1E-12 )
  {
    e += (cerror*cerror)/(c*c);
  }
  
  return TMath::Sqrt(e);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::ErrorABCD(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror, Double_t d, Double_t derror)
{
  /// Compute the quadratic sum of 4 errors
  
  Double_t e(0.0);
  
  if ( TMath::Abs(a) > 1E-12 )
  {
    e += (aerr*aerr)/(a*a);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (berr*berr)/(b*b);
  }
  
  if ( TMath::Abs(c) > 1E-12 )
  {
    e += (cerror*cerror)/(c*c);
  }

  if ( TMath::Abs(d) > 1E-12 )
  {
    e += (derror*derror)/(d*d);
  }

  return TMath::Sqrt(e);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::ErrorABCDE(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror, Double_t d, Double_t derror, Double_t ee, Double_t eeerror)
{
  /// Compute the quadratic sum of 4 errors
  
  Double_t e(0.0);
  
  if ( TMath::Abs(a) > 1E-12 )
  {
    e += (aerr*aerr)/(a*a);
  }
  
  if ( TMath::Abs(b) > 1E-12 )
  {
    e += (berr*berr)/(b*b);
  }
  
  if ( TMath::Abs(c) > 1E-12 )
  {
    e += (cerror*cerror)/(c*c);
  }
  
  if ( TMath::Abs(d) > 1E-12 )
  {
    e += (derror*derror)/(d*d);
  }

  if ( TMath::Abs(e) > 1E-12 )
  {
    e += (eeerror*eeerror)/(ee*ee);
  }

  return TMath::Sqrt(e);
}


//_____________________________________________________________________________
void AliAnalysisMuMuResult::Exclude(const char* subResultList)
{
  // exclude some subresult names from the list of subresult
  // to be used when computing the mean 
  
  TString slist(subResultList);
  
  TString tobeincluded = GetSubResultNameList();
  
  if ( slist == "*" )
  {
    Exclude(tobeincluded);
    return;
  }
  
  if ( fSubResultsToBeIncluded )
  {
    TObjArray* a = slist.Tokenize(",");
    TIter nextA(a);
    TObjString* s;

    while ( ( s = static_cast<TObjString*>(nextA())) )
    {
      TObject* o = fSubResultsToBeIncluded->FindObject(s->String());
      fSubResultsToBeIncluded->Remove(o);
    }
    
    delete a;
  }
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::GetErrorStat(const char* name, const char* subResultName) const
{
  // compute the mean value from all subresults that are included
  
  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    return sub->GetErrorStat(name);
  }
  
  if ( fMap )
  {
    TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
    if (p)
    {
      TParameter<double>* val = static_cast<TParameter<double>*>(p->At(kErrorStat));
      return val->GetVal();
    }
  }
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  Int_t n(0);
  Double_t v1(0);
  Double_t v2(0);
  Double_t d(0);
  
  Double_t mean = GetValue(name);
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {    
    if ( IsIncluded(r->Alias()) && r->HasValue(name) )
    {
      Double_t error = r->GetErrorStat(name);
      if (error)
      {
        v1 += 1.0/(error*error);
        v2 += 1.0/(error*error*error*error);
       
        d += ( (r->GetValue(name) - mean)*(r->GetValue(name)-mean) / (error*error));
      }
      ++n;
    }
  }
  
  if ( n<1 ) return 0.0;
  
  if ( n == 1 )
  {
    next.Reset();
    while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
    {
      if ( IsIncluded(r->Alias()) && r->HasValue(name) )
      {
        return r->GetErrorStat(name);
      }
    }
    return 0.0;
  }
  
  Double_t variance= (1.0/v1)*(1.0/(n-1))*d;
    // variance corrected by over/under-estimation of errors
    // i.e. scaled by chisquare per ndf
  
  return TMath::Sqrt(variance);
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::GetRMS(const char* name, const char* subResultName) const
{
  // compute the rms of the subresults
  // returns zero if no subresults

  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    return sub->GetRMS(name);
  }
  
  if ( fMap )
  {
    TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
    if (p)
    {
      TParameter<double>* val = static_cast<TParameter<double>*>(p->At(kRMS));
      return val ? val->GetVal() : 0.0; // val can be null for old results which did not have the rms set
    }
  }
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  Double_t v1(0);
  Double_t v2(0);
  Double_t sum(0);
  Int_t n(0);
  
  Double_t xmean = GetValue(name);
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( IsIncluded(r->Alias()) && r->HasValue(name) )
    {
      if ( r->GetErrorStat(name) > 0 )
      {
        Double_t ei = r->GetErrorStat(name);
        Double_t wi = 1.0/(ei*ei);
        v1 += wi;
        v2 += wi*wi;
        Double_t diff = r->GetValue(name) - xmean;
        sum += wi*diff*diff;
        ++n;
      }
    }
  }
  
  if ( n < 1 ) return 0.0;
  
  if ( n == 1 )
  {
    next.Reset();
    while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
    {
      if ( IsIncluded(r->Alias()) && r->HasValue(name) )
      {
        return r->GetRMS(name);
      }
    }
  }
  
  Double_t unbiased = TMath::Sqrt( (v1/(v1*v1-v2)) * sum);
  
  Double_t biased = TMath::Sqrt( sum/v1 );
  
  AliDebug(1,Form("v1 %e v1*v1 %e v2 %e -> biased %e unbiased %e (ratio %e)",v1,v1*v1,v2,biased,unbiased,unbiased/biased));
  
  return unbiased;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuResult::GetSubResultNameList() const
{
  // get a comma separated list of our subresult aliases
  TString tobeincluded;
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next())) )
  {
    if (tobeincluded.Length()>0) tobeincluded+=",";
    tobeincluded += r->Alias();
  }
  return tobeincluded;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuResult::GetValue(const char* name, const char* subResultName) const
{
  // get a value (either directly or by computing the mean of the subresults)
  
  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return TMath::Limits<Double_t>::Max();
    }
    return sub->GetValue(name);
  }
  
  if (fMap)
  {
    TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
    if (p)
    {
      TParameter<double>* val = static_cast<TParameter<double>*>(p->At(kValue));
      return val->GetVal();
    }
  }
  
  // compute the mean value from all subresults
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  Double_t mean(0);
  Double_t errorSum(0.0);
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( IsIncluded(r->Alias()) && r->HasValue(name) )
    {
      Double_t e = r->GetErrorStat(name);
      Double_t e2 = e*e;
      if ( e != 0.0 )
      {
        mean += r->GetValue(name)/e2;
        errorSum += 1.0/e2;
      }
    }
  }
  if ( errorSum != 0.0 )
  {
    return mean/errorSum;
  }
  else
  {
    return TMath::Limits<Double_t>::Max();
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuResult::HasValue(const char* name, const char* subResultName) const
{
  /// Whether this result (or subresult if subResultName is provided) has a property
  /// named "name"
  
  if ( strlen(subResultName) > 0 )
  {
    if ( !fSubResults)
    {
      AliError(Form("No subresult from which I could get the %s one...",subResultName));
      return kFALSE;
    }
    AliAnalysisMuMuResult* sub = static_cast<AliAnalysisMuMuResult*>(fSubResults->FindObject(subResultName));
    if (!sub)
    {
      AliError(Form("Could not get subresult named %s",subResultName));
      return kFALSE;
    }
    return sub->HasValue(name);
  }

  if ( fMap && ( fMap->GetValue(name) != 0x0 ) )
  {
    return kTRUE;
  }
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;

  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( r->HasValue(name) ) return kTRUE;
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::Include(const char* subResultList)
{
  // (re)include some subresult names
  
  TString slist(subResultList);

  if ( slist.Length()==0 )
  {
    Exclude("*");
    return;
  }
  
  if ( slist == "*" )
  {
    slist = GetSubResultNameList();
  }

  TObjArray* a = slist.Tokenize(",");
  a->SetOwner(kFALSE);
  TIter next(a);
  TObjString* s;
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    if (!fSubResultsToBeIncluded )
    {
      fSubResultsToBeIncluded = new TList;
      fSubResultsToBeIncluded->SetOwner(kTRUE);
    }
    if (!IsIncluded(s->String()))
    {
      fSubResultsToBeIncluded->Add(s);      
    }
  }
  
  delete a;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuResult::IsIncluded(const TString& alias) const
{
  // whether that subresult alias should be included when computing means, etc...
  
  if (!fSubResultsToBeIncluded) return kTRUE;
  
  return ( fSubResultsToBeIncluded->FindObject(alias) != 0x0 );
}

//_____________________________________________________________________________
THashList* AliAnalysisMuMuResult::Keys() const
{
  /// Return the complete list of keys we're using
  if (!fKeys)
  {
    fKeys = new THashList;
    fKeys->SetOwner(kTRUE);
    TIter next(fMap);
    TObjString* key;
    
    while ( ( key = static_cast<TObjString*>(next()) ) )
    {
      if ( !fKeys->FindObject(key->String()) )
      {
        fKeys->Add(new TObjString(key->String()));
      }
    }
    
    AliAnalysisMuMuResult* r;
    TIter nextResult(fSubResults);
    
    while ( ( r = static_cast<AliAnalysisMuMuResult*>(nextResult())) )
    {
      TIter nextHL(r->Keys());
      TObjString* s;
      
      while ( ( s = static_cast<TObjString*>(nextHL())) )
      {
        if ( !fKeys->FindObject(s->String()) )
        {
          fKeys->Add(new TObjString(s->String()));
        }
      }
    }

  }
  return fKeys;
}

//_____________________________________________________________________________
Long64_t AliAnalysisMuMuResult::Merge(TCollection* list)
{
  /// Merge method
  ///
  /// Merge a list of AliAnalysisMuMuResult objects with this
  /// Returns the number of merged objects (including this).
  ///
  /// Note that the merging is to be understood here as a weighed mean operation
  
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  TIter next(list);
  TObject* currObj;

  Double_t sumw(Weight()); // sum of weights
  
  while ( ( currObj = next() ) )
  {
    AliAnalysisMuMuResult* result = dynamic_cast<AliAnalysisMuMuResult*>(currObj);
    if (!result)
    {
      AliFatal(Form("object named \"%s\" is a %s instead of an AliAnalysisMuMuResult!", currObj->GetName(), currObj->ClassName()));
      continue;
    }
    
    sumw += result->Weight();
  }
  
  TIter nextKey(Keys());
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(nextKey())) )
  {
    Double_t value = GetValue(key->String())*Weight()/sumw;
    Double_t e = GetErrorStat(key->String());
    Double_t e2 = e*e*Weight()*Weight()/sumw/sumw;

    next.Reset();
    
    while ( ( currObj = next() ) )
    {
      AliAnalysisMuMuResult* result = dynamic_cast<AliAnalysisMuMuResult*>(currObj);
    
      if (!result)
      { 
        continue;
      }
      
      if (!result->HasValue(key->String()))
      {
        AliError(Form("Did not find key %s in of the result to merge",key->String().Data()));
        continue;
      }
      
      // can only merge under the condition we have the same bin
    
      Double_t w = result->Weight()/sumw;
      
      Double_t w2 = w*w;
      
      value += result->GetValue(key->String())*w;
      e2 += result->GetErrorStat(key->String())*result->GetErrorStat(key->String())*w2;
      
    }
    
    Set(key->String(),
        value,
        TMath::Sqrt(e2));
  }
  
  TIter nextSubresult(fSubResults);
  AliAnalysisMuMuResult* r;
  
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(nextSubresult())) )
  {
    TList sublist;

    next.Reset();
    
    while ( ( currObj = next() ) )
    {
      sublist.Add(currObj);
    }

    r->Merge(&sublist);
  }
  
  fWeight = sumw;
  
  return list->GetEntries()+1;
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuResult::NofIncludedSubResults(const char* name) const
{
  // Return the number of subresults which have key name and are included
  
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  Int_t n(0);
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if ( IsIncluded(r->Alias()) && r->HasValue(name) )
    {
      ++n;
    }
  }
  return n;
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::Print(Option_t* opt) const
{
  /// printout
  
  TString sopt(opt);
  sopt.ToUpper();
  
  for ( Int_t i = 0; i < 9; ++i )
  {
    sopt.ReplaceAll(Form("%d",i),"");
  }
  
  TString pot(sopt);
  pot.ReplaceAll("ALL","");
  pot.ReplaceAll("FULL","");

  std::cout << pot.Data();
  
    if ( fAlias.Length() > 0 )
    {
      std::cout << Form("%s - ",fAlias.Data());
    }
  
    std::cout << Form("%s %s %s",
                      GetName(),GetTitle(),fWeight > 0.0 ? Form(" WEIGHT %e",fWeight) : "");
  
  if ( fSubResults && fSubResults->GetEntries()>1 )
  {
    std::cout << " (" << fSubResults->GetEntries() << " subresults)";
  }

  std::cout << std::endl;
  
  TIter next(Keys());
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(next())) )
  {
    PrintValue(key->String().Data(),pot.Data(),
               GetValue(key->String()),
               GetErrorStat(key->String()),
               GetRMS(key->String()));
  }

  if ( fSubResults /* && fSubResults->GetEntries() > 1 */ && ( sopt.Contains("ALL") || sopt.Contains("FULL") ) )
  {
    std::cout << pot.Data() << "\t===== sub results ===== " << std::endl;
    
    sopt += "\t\t";
    
    TIter nextSubresult(fSubResults);
    AliAnalysisMuMuResult* r;
    
    while ( ( r = static_cast<AliAnalysisMuMuResult*>(nextSubresult()) ) )
    {
      if ( !IsIncluded(r->Alias()) )
      {
        std::cout << " [EXCLUDED]";
      }
      r->Print(sopt.Data());
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::PrintValue(const char* key, const char* opt,
                                       Double_t value, Double_t errorStat, Double_t rms) const
{
  // print one value and its associated error
  
  if ( TString(key).Contains("AccEff") )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f %% (%5.2f %%)",key,value*100,errorStat*100,
                             value != 0.0 ? errorStat*100.0/value : 0.0 );
    
    if ( rms )
    {
      std::cout << Form(" RMS %9.2f (%5.2f %%)",rms,100.0*rms/value);
    }

    std::cout << std::endl;
  }
  else if ( TString(key).BeginsWith("Sigma") ||  TString(key).BeginsWith("Mean") )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f (%5.2f %%) MeV/c^2",key,value*1E3,1E3*errorStat,
                             value != 0.0 ? errorStat*100.0/value : 0.0);
    
    if ( rms )
    {
      std::cout << Form(" RMS %9.2f (%5.2f %%)",rms,100.0*rms/value);
    }

    std::cout << std::endl;
  }
  else if ( TString(key).Contains("Nof") || ( TString(key).Contains("Fnorm") && !TString(key).Contains("persion") ) )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f (%5.2f %%)",key,value,errorStat,
                             value != 0.0 ? errorStat*100.0/value : 0.0);
    
    if ( rms )
    {
      std::cout << Form(" RMS %9.2f (%5.2f %%)",rms,100.0*rms/value);
    }
    std::cout << std::endl;
  }
  else if ( value > 1E-3 && value < 1E3 )
  {
    std::cout << opt << Form("\t\t%20s %9.2f +- %5.2f (%5.2f %%)",key,value,errorStat,
                             value != 0.0 ? errorStat*100.0/value : 0.0);
    if ( rms )
    {
      std::cout << Form(" RMS %9.2f (%5.2f %%)",rms,100.0*rms/value);
    }
    std::cout << std::endl;
  }
  else
  {
    std::cout << opt << Form("\t\t%20s %9.2e +- %9.2e (%5.2f %%)",key,value,errorStat,
                             value != 0.0 ? errorStat*100.0/value : 0.0);
    if ( rms )
    {
      std::cout << Form(" RMS %9.2e (%5.2f %%)",rms,100.0*rms/value);
    }
    std::cout << std::endl;
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::Scale(Double_t w)
{
  /// Scale all our internal values by value
  
  TIter next(Keys());
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(next())) )
  {
    Double_t value = GetValue(key->String());
    Double_t error = GetErrorStat(key->String());
    Double_t rms = GetRMS(key->String());
    
    Set(key->String(),value*w,error*w,rms*w);
  }

}

//_____________________________________________________________________________
void AliAnalysisMuMuResult::Set(const char* name, Double_t value, Double_t errorStat, Double_t rms)
{
  /// Set a (value,error) pair with a given name
  
  if ( !fMap )
  {
    fMap = new TMap;
    fMap->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  
  TObjArray* p = static_cast<TObjArray*>(fMap->GetValue(name));
  if (!p)
  {
    p = new TObjArray(4);
    
    p->SetOwner(kTRUE);
    
    p->AddAt(new TParameter<Double_t>(name,value),kValue);
    p->AddAt(new TParameter<Double_t>(name,errorStat),kErrorStat);
    p->AddAt(new TParameter<Double_t>(name,rms),kRMS);
    
    fMap->Add(new TObjString(name),p);
  }
  else
  {
    static_cast<TParameter<double>*>(p->At(kValue))->SetVal(value);
    static_cast<TParameter<double>*>(p->At(kErrorStat))->SetVal(errorStat);
    static_cast<TParameter<double>*>(p->At(kRMS))->SetVal(rms);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuResult*
AliAnalysisMuMuResult::SubResult(const char* subResultName) const
{
  /// get a given subresult
  if (!fSubResults)
  {
    return 0x0;
  }
  TIter next(fSubResults);
  AliAnalysisMuMuResult* r;
  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next())) )
  {
    if ( r->Alias() == subResultName )
    {
      return r;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
TList* AliAnalysisMuMuResult::SubResultsToBeIncluded() const
{
  if (!fSubResultsToBeIncluded)
  {
    fSubResultsToBeIncluded = new TList;
    fSubResultsToBeIncluded->SetOwner(kTRUE);
  }
  return fSubResultsToBeIncluded;
}

