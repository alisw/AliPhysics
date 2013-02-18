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
  
}

//______________________________________________________________________________
AliAnalysisMuMuSpectra::~AliAnalysisMuMuSpectra()
{
  delete fBinning;
  delete fBins;
}

//______________________________________________________________________________
void AliAnalysisMuMuSpectra::AdoptResult(const AliAnalysisMuMuBinning::Range& bin,
                                         AliAnalysisMuMuResult* result)
{
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
    AliInfoClass(Form("i=%d",i ));
    StdoutToAliInfoClass(thisResult->Print("full");
                         std::cout << "----" << std::endl;
                         accResult->Print("full"));
    
    thisResult->Correct(*accResult,particle,subResultName);
  }
  
  delete particles;
  delete types;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisMuMuSpectra::IsEmpty() const
{
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
TH1* AliAnalysisMuMuSpectra::Plot(const char* what, const char* subresult) const
{
  TString swhat(what);
  swhat.ToUpper();
  
  Double_t* bins = fBinning->CreateBinArray();
  
  TObjArray* binArray = fBinning->CreateBinObjArray();
  
  TH1* h(0x0);
  
  AliDebugClass(1,Form("nbins=%d nresults=%d",binArray->GetEntries(),fBins->GetEntries()));
  
  for ( Int_t j = 0; j < TMath::Min(binArray->GetEntries(),fBins->GetEntries()); ++j )
  {
    AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(fBins->At(j));
    
    if ( strlen(subresult) > 0 )
    {
      TObjArray* sr = r->SubResults();
      if (!sr) continue;
      TString sub(subresult);
      sub.ToUpper();
      r = static_cast<AliAnalysisMuMuResult*>(sr->FindObject(sub.Data()));
      if (!r) continue;
    }
    
    const AliAnalysisMuMuBinning::Range& b = r->Bin();
    
    r->Print();
    b.Print();
    
    if (!h)
    {
      h = new TH1F(r->GetName(),r->GetName(),binArray->GetEntries(),bins);
      h->SetDirectory(0);
    }
    
    Double_t y = r->GetValue(what);
    Double_t yerr = r->GetErrorStat(what);
    
//    if ( !swhat.BeginsWith("ACC") )
      if ( swhat.Contains("NOF") )
    {
      y /= (b.WidthX());
      yerr /= (b.WidthX());
    }
    
    h->SetBinContent(j+1,y);
    h->SetBinError(j+1,yerr);
    
    AliInfoClass(Form("%e +- %e",y,yerr));
  }
  
  delete binArray;
  delete[] bins;
  
  return h;
}


//______________________________________________________________________________
void AliAnalysisMuMuSpectra::Print(Option_t* opt) const
{
  if (!IsEmpty())
  {
    TString sopt(opt);
    sopt.ToUpper();
    if ( sopt.Contains("BINNING") )
    {
      fBinning->Print(opt);
    }
    Int_t nmax = sopt.Atoi();
    if ( nmax <= 0 ) nmax = fBins->GetEntries();
    for ( Int_t i = 0; i < nmax; ++i )
    {
      AliAnalysisMuMuBinning::Range* r = static_cast<AliAnalysisMuMuBinning::Range*>(fBins->At(i));
      if (r) r->Print(opt);
    }
    
  }
}
