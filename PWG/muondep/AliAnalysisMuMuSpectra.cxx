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


// AliAnalysisMuMuSpectra : a class to encapsulate results from MuMu analysis
//
// Spectra can be merged and converted into histograms
//
// author: L. Aphecetche (Subatech)
//

#include "AliAnalysisMuMuSpectra.h"

#include "AliLog.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuResult.h"
#include "Riostream.h"
#include "TH1.h"
#include "TList.h"
#include "TObjArray.h"

ClassImp(AliAnalysisMuMuSpectra)

//______________________________________________________________________________
AliAnalysisMuMuSpectra::AliAnalysisMuMuSpectra(const char* name, const char* title) :
TNamed(name,title),
fBinning(0x0),
fBins(0x0)
{
  // default ctor
}

//______________________________________________________________________________
AliAnalysisMuMuSpectra::AliAnalysisMuMuSpectra(const AliAnalysisMuMuSpectra& rhs)
: TNamed(rhs.GetName(),rhs.GetTitle()),
fBinning(0x0),
fBins(0x0)
{
  // copy ctor

  if ( rhs.fBinning )
  {
    fBinning = new AliAnalysisMuMuBinning(*rhs.fBinning);
  }

  TIter next(rhs.fBins);
  AliAnalysisMuMuBinning::Range* bin;
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if (!fBins)
    {
      fBins = new TObjArray;
      fBins->SetOwner(kTRUE);
    }
    fBins->Add(bin);
  }
}

//______________________________________________________________________________
AliAnalysisMuMuSpectra&
AliAnalysisMuMuSpectra::operator=(const AliAnalysisMuMuSpectra& rhs)
{
  // assignment operator
  
  if (this==&rhs) return *this;
  
  delete fBinning;
  fBinning = 0x0;
  delete fBins;
  fBins = 0x0;
  
  if ( rhs.fBinning )
  {
    fBinning = new AliAnalysisMuMuBinning(*rhs.fBinning);
  }
  
  TIter next(rhs.fBins);
  AliAnalysisMuMuBinning::Range* bin;
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if (!fBins)
    {
      fBins = new TObjArray;
      fBins->SetOwner(kTRUE);
    }
    fBins->Add(bin);
  }
  
  return *this;
}

//______________________________________________________________________________
AliAnalysisMuMuSpectra::~AliAnalysisMuMuSpectra()
{
  // dtor
  delete fBinning;
  delete fBins;
}

//______________________________________________________________________________
void AliAnalysisMuMuSpectra::AdoptResult(const AliAnalysisMuMuBinning::Range& bin,
                                         AliAnalysisMuMuResult* result)
{
  // adopt (i.e. we are becoming the owner) a result for a given bin
  if (!fBinning)
  {
    fBinning = new AliAnalysisMuMuBinning;
    fBins = new TObjArray;
    fBins->SetOwner(kTRUE);
  }
  fBinning->AddBin(bin);
  fBins->Add(result);
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuSpectra::Correct(const AliAnalysisMuMuSpectra& accEff, const char* particle, const char* subResultName)
{
  /// Correct this spectra by acceff one
  if (IsEmpty())
  {
    AliError("Cannot correct an empty spectra !");
    return kFALSE;
  }
  
  // check we have the same binning first
  AliAnalysisMuMuBinning* accEffBins = accEff.Binning();
  
  if ( !fBinning->IsEqual(accEffBins) )
  {
    AliError("Cannot correct with a spectra which does not have the same binning");
    return kFALSE;
  }
  
  TObjArray* particles = fBinning->CreateParticleArray();
  TObjArray* types = fBinning->CreateTypeArray();
  
  if (particles->GetEntries()!=1 || types->GetEntries()!=1 )
  {
    delete particles;
    delete types;
    return kFALSE;
  }
  
  TObjArray* bins = accEff.Bins();

  
  for ( Int_t i = 0; i < bins->GetEntries(); ++i )
  {
    AliAnalysisMuMuResult* thisResult = static_cast<AliAnalysisMuMuResult*>(fBins->At(i));
    AliAnalysisMuMuResult* accResult = static_cast<AliAnalysisMuMuResult*>(bins->At(i));
//    AliInfoClass(Form("i=%d",i ));
//    StdoutToAliInfoClass(thisResult->Print("full");
//                         std::cout << "----" << std::endl;
//                         accResult->Print("full"));
    
    thisResult->Correct(*accResult,particle,subResultName);
  }
  
  delete particles;
  delete types;
  return kTRUE;
}

//______________________________________________________________________________
AliAnalysisMuMuResult*
AliAnalysisMuMuSpectra::GetResultForBin(const AliAnalysisMuMuBinning::Range& bin) const
{
  /// Get result for a given bin
  /// Warning: this requires a loop on bins
  
  if ( IsEmpty() ) return 0x0;
  
  TObjArray* bins = fBinning->CreateBinObjArray();
  
  Int_t j(-1);
  
  StdoutToAliDebug(1,std::cout << "searching for "; bin.Print());
  
  for ( Int_t i = 0; i <= bins->GetLast() && j < 0 ; ++i )
  {
    AliAnalysisMuMuBinning::Range* b = static_cast<AliAnalysisMuMuBinning::Range*>(bins->At(i));

    StdoutToAliDebug(1,b->Print(););
    
    if ( bin == *b )
    {
      j = i;
    }
  }
  
  delete bins;
  
  if (j>=0)
  {
    return static_cast<AliAnalysisMuMuResult*>(fBins->At(j));
  }
  else
  {
    StdoutToAliDebug(1,std::cout << "Could not find result for bin:" << std::endl; bin.Print(););
  }
  return 0x0;
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuSpectra::HasValue(const char* what) const
{
    // whether or not our result(s) has a given property
  if ( IsEmpty() ) return kFALSE;
  
  AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(fBins->First());
  
  return r->HasValue(what);
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuSpectra::IsEmpty() const
{
  // whether this spectra is empty or not
  return ( fBins==0x0 || fBins->GetEntries()<=0 );
}

//______________________________________________________________________________
Long64_t AliAnalysisMuMuSpectra::Merge(TCollection* list)
{
  /// Merge method
  
  // Merge a list of AliAnalysisMuMuSpectra objects with this
  // Returns the number of merged objects (including this).
  
  if (!list) return 0;
  
  if (list->IsEmpty()) return 1;
  
  TIter next(list);
  TObject* currObj;
  Int_t count(0);
  
  TList binningList;
  
  for ( Int_t i = 0; i <= fBins->GetLast(); ++i )
  {
    next.Reset();
    
    TList binList;
    
    while ( ( currObj = next() ) )
    {
      AliAnalysisMuMuSpectra* spectra = dynamic_cast<AliAnalysisMuMuSpectra*>(currObj);
      if (!spectra)
      {
        AliFatal(Form("object named \"%s\" is a %s instead of an AliAnalysisMuMuSpectra!", currObj->GetName(), currObj->ClassName()));
        continue;
      }
    
      if (i==0)
      {
        binningList.Add(spectra->Binning());
  
        if ( !fBinning->IsEqual(spectra->Binning()) || spectra->Bins()->GetLast() != Bins()->GetLast() )
        {
          AliError("Cannot merge spectra with different binning");
          continue;
        }
        
        ++count;
      }
      
      binList.Add(fBins->At(i));
    }
    
    AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(fBins->At(i));
    r->Merge(&binList);
    
  }
  
  fBinning->Merge(&binningList);
  
  return count+1;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuSpectra::Plot(const char* what, const char* subresult, Bool_t divideByBinWidth) const
{
  // Convert the spectra into an histogram
  
  TString swhat(what);
  swhat.ToUpper();
  
  Double_t* bins = fBinning->CreateBinArray();
  
  TObjArray* binArray = fBinning->CreateBinObjArray();
  
  TH1* h(0x0);
  
  AliDebug(1,Form("nbins=%d nresults=%d",binArray->GetEntries(),fBins->GetEntries()));
  
  for ( Int_t j = 0; j < TMath::Min(binArray->GetEntries(),fBins->GetEntries()); ++j )
  {
    AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(fBins->At(j));
    
    if ( strlen(subresult) > 0 && r->SubResults() )
    {
      TString sub(subresult);
      sub.ToUpper();
      r = r->SubResult(sub.Data());
      if (!r) continue;
    }
    
    const AliAnalysisMuMuBinning::Range& b = r->Bin();
    
    if (!h)
    {
      h = new TH1F(r->GetName(),r->GetName(),binArray->GetEntries(),bins);
      h->SetDirectory(0);
    }
    
    Double_t y = r->GetValue(what);
    Double_t yerr = r->GetErrorStat(what);
  
    if ( divideByBinWidth && b.WidthX()>0 )
    {
      y /= (b.WidthX());
      yerr /= (b.WidthX());
    }
    
    std::cout << b.AsString();
    AliAnalysisMuMuResult::PrintValue(swhat.Data(),"",y,yerr);
    
    h->SetBinContent(j+1,y);
    h->SetBinError(j+1,yerr);
    
  }
  
  delete binArray;
  delete[] bins;
  
  return h;
}


//______________________________________________________________________________
void AliAnalysisMuMuSpectra::Print(Option_t* opt) const
{
  // printout
  
  if (!IsEmpty())
  {
    TString sopt(opt);
    Int_t nmax = sopt.Atoi();
    if ( nmax <= 0 ) nmax = fBins->GetEntries();
    for ( Int_t i = 0; i < nmax; ++i )
    {
      AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(fBins->At(i));
      if (r) r->Print(opt);
    }
  }
}
