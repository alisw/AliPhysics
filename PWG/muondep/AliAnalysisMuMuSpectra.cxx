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

#include "AliAnalysisMuMuSpectra.h"

#include "AliLog.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "Riostream.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TObjArray.h"

ClassImp(AliAnalysisMuMuSpectra)

//______________________________________________________________________________
AliAnalysisMuMuSpectra::AliAnalysisMuMuSpectra(const char* name, const char* title) :
TNamed(name,title),
fBinning(0x0),
fBins(0x0),
fWeight(1.0)
{
  /// default ctor
}

//______________________________________________________________________________
AliAnalysisMuMuSpectra::AliAnalysisMuMuSpectra(const AliAnalysisMuMuSpectra& rhs)
: TNamed(rhs.GetName(),rhs.GetTitle()),
fBinning(0x0),
fBins(0x0),
fWeight(rhs.Weight())
{
  /// copy ctor

  if ( rhs.fBinning )
  {
    fBinning = new AliAnalysisMuMuBinning(*rhs.fBinning);
  }

  TIter next(rhs.fBins);
  AliAnalysisMuMuResult* bin;

  while ( ( bin = static_cast<AliAnalysisMuMuResult*>(next()) ) )
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
  /// assignment operator

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
  AliAnalysisMuMuResult* bin;

  while ( ( bin = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    if (!fBins)
    {
      fBins = new TObjArray;
      fBins->SetOwner(kTRUE);
    }
    fBins->Add(bin);
  }

  fWeight = rhs.Weight();

  return *this;
}

//______________________________________________________________________________
AliAnalysisMuMuSpectra::~AliAnalysisMuMuSpectra()
{
  /// dtor
  delete fBinning;
  delete fBins;
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuSpectra::AdoptResult(const AliAnalysisMuMuBinning::Range& bin,
                                         AliAnalysisMuMuResult* result)
{
  /// adopt (i.e. we are becoming the owner) a result for a given bin
  if ( !result )
  {
    AliError("Cannot adopt a null result");
    return kFALSE
    ;
  }
  if (!fBinning)
  {
    fBinning = new AliAnalysisMuMuBinning;
    fBins = new TObjArray;
    fBins->SetOwner(kTRUE);
  }
  fBinning->AddBin(bin); // Add a bin to the spectra

  Int_t sizeBeforeAdd = fBins->GetEntriesFast();
  fBins->Add(result); // Add the result to the recently added bin
  Int_t sizeAfterAdd = fBins->GetEntriesFast();
  if ( sizeBeforeAdd >= sizeAfterAdd )
  {
   AliError(Form("Error adopting result %s to spectra %s",result->GetName(),GetName()));
    return kFALSE;
  }
  else return kTRUE;
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

  TObjArray* particles = fBinning->CreateWhatArray();
  TObjArray* types = fBinning->CreateQuantityArray();

  // FIXME : cross-check what happens here with a binning containing
  // centralities, etc...

  if (particles->GetEntries()!=1 || types->GetEntries()!=1 )
  {
    delete particles;
    delete types;
    return kFALSE;
  }

  TObjArray* bins = accEff.BinContentArray();


  for ( Int_t i = 0; i < bins->GetEntries(); ++i )
  {
    AliAnalysisMuMuJpsiResult* thisResult = static_cast<AliAnalysisMuMuJpsiResult*>(fBins->At(i));
    AliAnalysisMuMuJpsiResult* accResult = static_cast<AliAnalysisMuMuJpsiResult*>(bins->At(i));
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
AliAnalysisMuMuResult*
AliAnalysisMuMuSpectra::GetResultForBin(const char* binName) const
{
  /// Get result for a given bin
  /// Warning: this requires a loop on bins

  if ( IsEmpty() ) return 0x0;

  TObjArray* bins = fBinning->CreateBinObjArray();

  Int_t j(-1);

  AliDebug(1,Form("searching for bin %s",binName));

  for ( Int_t i = 0; i <= bins->GetLast() && j < 0 ; ++i )
  {
    AliAnalysisMuMuBinning::Range* b = static_cast<AliAnalysisMuMuBinning::Range*>(bins->At(i));

    if ( b->AsString() == binName )
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
    AliDebug(1,Form("Could not find result for bin %s",binName));
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
  /// whether this spectra is empty or not
  return ( fBins==0x0 || fBins->GetEntries()<=0 );
}

//______________________________________________________________________________
Long64_t AliAnalysisMuMuSpectra::Merge(TCollection* list)
{
  /// Merge method
  ///
  /// Merge a list of AliAnalysisMuMuSpectra objects with this
  /// Returns the number of merged objects (including this).

  if (!list) return 0;

  if (list->IsEmpty()) return 1;

  TIter next(list);
  TObject* currObj;
  Int_t count(0);

  TList binningList;

  // for each bin must do a list of results, and merge that list

  TObjArray* bins = fBinning->CreateBinObjArray();
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* bin;

  Int_t i(0);

  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    next.Reset();

    TList binList;

    while ( ( currObj = next() ) )
    {
      AliAnalysisMuMuSpectra* spectra = static_cast<AliAnalysisMuMuSpectra*>(currObj);

      if (i==0)
      {
        binningList.Add(spectra->Binning());

        if ( !fBinning->IsEqual(spectra->Binning()) || spectra->BinContentArray()->GetLast() != BinContentArray()->GetLast() )
        {
          AliError("Cannot merge spectra with different binning");
          continue;
        }
        ++count;
      }

      binList.Add(spectra->GetResultForBin(*bin));
    }

    ++i;

    AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(GetResultForBin(*bin));
    r->Merge(&binList);

  }

  delete bins;

  return count+1;
}

//______________________________________________________________________________
Long64_t AliAnalysisMuMuSpectra::Merge(AliAnalysisMuMuSpectra* spectraToAdd)
{
  /// Merge method
  ///
  /// Merge a list of AliAnalysisMuMuSpectra objects with this
  /// Returns the number of merged objects (including this).

  if (!spectraToAdd) return 0;

  Int_t count(0);

  TList binningList;

  // for each bin must do a list of results, and merge that list

  TObjArray* bins = fBinning->CreateBinObjArray();
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* bin;

  Int_t i(0);

  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {

    TList binList;

    AliAnalysisMuMuSpectra* spectra = static_cast<AliAnalysisMuMuSpectra*>(spectraToAdd);

    if (i==0)
    {
      binningList.Add(spectra->Binning());

      if ( !fBinning->IsEqual(spectra->Binning()) || spectra->BinContentArray()->GetLast() != BinContentArray()->GetLast() )
      {
        AliError("Cannot merge spectra with different binning");
        continue;
      }
      ++count;
    }

    binList.Add(spectra->GetResultForBin(*bin));

    ++i;

    AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(GetResultForBin(*bin));
    r->Merge(&binList);

  }

  delete bins;

  return count+1;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuSpectra::Plot(const char* what, const char* subresult, Bool_t divideByBinWidth) const
{
  /// Convert the spectra into an histogram

  TString swhat(what);
  TString swhatT(what);
  swhat.ToUpper();

  Double_t* bins  = fBinning->CreateBinArray();
  Double_t* binsX = fBinning->CreateBinArrayX();
  Double_t* binsY = fBinning->CreateBinArrayY();

  Int_t nbinsX    = fBinning->GetNBinsX();
  Int_t nbinsY    = fBinning->GetNBinsY();

  TObjArray* binArray = fBinning->CreateBinObjArray();

  TH1* h(0x0);

  AliDebug(1,Form("nbins=%d nresults=%d",binArray->GetEntries(),fBins->GetEntries()));

  for ( Int_t j = 0; j < TMath::Min(binArray->GetEntries(),fBins->GetEntries()); ++j )
  {
    AliAnalysisMuMuJpsiResult* r = static_cast<AliAnalysisMuMuJpsiResult*>(fBins->At(j));

    if ( strlen(subresult) > 0 && r->SubResults() )
    {
      TString sub(subresult);
      r = static_cast<AliAnalysisMuMuJpsiResult*>(r->SubResult(sub.Data()));
      if (!r) continue;
    }

    const AliAnalysisMuMuBinning::Range& b = r->Bin();

    if ( b.Is2D() )
    {
      if (!h)
      {
        h = new TH2F(r->GetName(),r->GetName(),nbinsX,binsX,nbinsY,binsY);
        h->SetDirectory(0);
      }

      Double_t z = r->GetValue(swhatT.Data());
      Double_t zerr = r->GetErrorStat(swhatT.Data());

      if ( swhat.Contains("AccEff",TString::kIgnoreCase) )
      {
        if ( z < 0. || z > 5.)
        {
          z = -1.;
          zerr = 0.;
        }
      }

      if ( divideByBinWidth && b.WidthX()>0 )
      {
        z /= (b.WidthX()*b.WidthY());
        zerr /= (b.WidthX()*b.WidthY());

      }

      std::cout << b.AsString();
      r->PrintValue(swhat.Data(),"",z,zerr);

      Int_t x = j/nbinsY + 1;
      Int_t y = (j % nbinsY) + 1;

      // std::cout << "x="  << x << ";" << "y=" << y << std::endl;

      h->SetBinContent(x,y,z);
      h->SetBinError(x,y,zerr);

    } else {

      if (!h)
      {
        if(b.IsIntegrated()) h = new TH1F(r->GetName(),r->GetName(),1,-1,1.);
        else h = new TH1F(r->GetName(),r->GetName(),binArray->GetEntries(),bins);
        h->SetDirectory(0);
      }

      Double_t y = r->GetValue(swhatT.Data());
      Double_t yerr = r->GetErrorStat(swhatT.Data());

      if ( divideByBinWidth && b.WidthX()>0 )
      {
        y /= (b.WidthX());
        yerr /= (b.WidthX());
      }

      if (!TMath::Finite(y)) y = 0.0;
      if (!TMath::Finite(yerr)) yerr = 0.0;

      std::cout << b.AsString();
      r->PrintValue(swhat.Data(),"",y,yerr);

      //__________Config. Histo
      h->SetBinContent(j+1,y);
      h->SetBinError(j+1,yerr);
      h->GetXaxis()->SetBinLabel(j+1,b.AsString().Data());
      h->GetYaxis()->SetTitle(swhat.Data());
      //__________
    }
  }

  delete binArray;
  delete[] binsX;
  delete[] binsY;
  delete[] bins;

  return h;
}


//______________________________________________________________________________
void AliAnalysisMuMuSpectra::Print(Option_t* opt) const
{
  /// printout

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

//______________________________________________________________________________
void AliAnalysisMuMuSpectra::Scale(Double_t value)
{
  /// scale all bins by value

  TIter next(fBins);
  AliAnalysisMuMuResult* r;

  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    r->Scale(value);
  }
}

//______________________________________________________________________________
void AliAnalysisMuMuSpectra::SetWeight(Double_t w)
{
  /// Set the weight of this spectra
  fWeight = w;
  TIter next(fBins);
  AliAnalysisMuMuResult* r;

  while ( ( r = static_cast<AliAnalysisMuMuResult*>(next()) ) )
  {
    r->SetWeight(Weight());
  }
}
