//
// Class AliRsnFunction
//
// This class defines a base classe to implement a function
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from 
// this class and override the virtual functions defined in it, which 
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TString.h>

#include "AliLog.h"

#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnCut.h"

#include "AliRsnFunction.h"

ClassImp(AliRsnFunction)

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction() :
  fFcnType(kFcnTypes),
  fUseBins(kFALSE),
  fSkipOutsideInterval(kFALSE),
  fBins(0),
  fBinningCut(),
  fBinningCutType(AliRsnCut::kLastCutType),
  fHistoDef(0x0)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //
  
  Int_t i;
  for (i = 0; i < 100; i++) {
    fHisto[i] = 0x0;
  }
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction
(EFcnType type, AliRsnHistoDef *hd, Bool_t skipOut) :
  fFcnType(type),
  fUseBins(kFALSE),
  fSkipOutsideInterval(skipOut),
  fBins(0),
  fBinningCut(),
  fBinningCutType(AliRsnCut::kLastCutType),
  fHistoDef(hd)
{
  //
  // Constructor.
  // The histogram data member cannot be passed externally,
  // its initialization MUST be defined inside the Init() method,
  // which must be overridden in any derivate implementation.
  //
  
  Int_t i;
  for (i = 0; i < 100; i++) {
    fHisto[i] = 0x0;
  }
}

//________________________________________________________________________________________
AliRsnFunction::AliRsnFunction(const AliRsnFunction &copy) :
  TObject(copy),
  fFcnType(copy.fFcnType),
  fUseBins(copy.fUseBins),
  fSkipOutsideInterval(copy.fSkipOutsideInterval),
  fBins(0),
  fBinningCut(),
  fBinningCutType(AliRsnCut::kLastCutType),
  fHistoDef(copy.fHistoDef)
{
  //
  // Copy constructor.
  // Calls the function to define binning.
  //

  Int_t i, n = 100;
  for (i = 0; i < n; i++) {
    fHisto[i] = 0x0;
  }
    
  if (fUseBins) {
    n = copy.fBins.GetSize();
    Double_t *array = new Double_t[n];
    for (i = 0; i < n; i++) array[i] = copy.fBins[i];
    SetBinningCut(copy.fBinningCutType, copy.fBins.GetSize(), array);
    delete [] array;
  }
}
//________________________________________________________________________________________
const AliRsnFunction& AliRsnFunction::operator=(const AliRsnFunction& /*copy*/)
{
  //
  // Assignment operator.
  // Behaves like copy constructor.
  // Also in this case, the histogram is not copied, and,
  // if it was present, it is destroyed and will need to be recreated.
  //

  return (*this);
}
//________________________________________________________________________________________
void AliRsnFunction::Clear(Option_t* /*option*/)
{
  //
  // Clear arrays and histogram.
  // For the sake of security, all pointers are also set explicitly to NULL.
  //

  Int_t i;
  for (i = 0; i < 100; i++) {
    delete fHisto[i];
    fHisto[i] = 0x0;
  }
}

//________________________________________________________________________________________
TList* AliRsnFunction::Init(const char *histoName, const char *histoTitle)
{
  //
  // Initialization function.
  // By default, it initializes the owned histogram using the method
  // from AliRsnHistoDef class, giving the same name and title of this.
  // A user can override this behaviour, if necessary.
  // Before creating, the HistoDef is checked for proper initialization.
  //

  Clear();
    
  Int_t i, ibin, nbins = fHistoDef->GetNBins();
  Double_t min = fHistoDef->GetMin(), max = fHistoDef->GetMax();
    
  // list is created and named after the general
  // settings used for the contained histograms
  TList *histos = new TList;
  histos->SetName(Form("%s", GetFcnName().Data()));
    
  // a general histogram is always added, 
  // which overrides the binning and collects everything
  fHisto[0] = new TH1D(histoName, histoTitle, nbins, min, max);
  histos->AddLast(fHisto[0]);
    
  // if requested a binning w.r. to some cut variable, histograms are added 
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  Char_t hTitle[255];
  if (fUseBins) {
    for (ibin = 0, i = 1; ibin < fBins.GetSize() - 1; ibin++, i++) {
      sprintf(hName, "%s[%.2f-%.2f]", histoName, fBins[ibin], fBins[ibin+1]);
      sprintf(hTitle, "%s [%.2f-%.2f]", histoTitle, fBins[ibin], fBins[ibin+1]);
      fHisto[i] = new TH1D(hName, hTitle, nbins, min, max);
      histos->AddLast(fHisto[i]);
    }
  }
    
  // returns the full list at the end
  return histos;
}

//________________________________________________________________________________________
void AliRsnFunction::Init(const char *histoName, const char *histoTitle, TList *histos)
{
  //
  // Initialization function.
  // By default, it initializes the owned histogram using the method
  // from AliRsnHistoDef class, giving the same name and title of this.
  // A user can override this behaviour, if necessary.
  // Before creating, the HistoDef is checked for proper initialization.
  //

  Clear();
    
  Int_t i, ibin, nbins = fHistoDef->GetNBins();
  Double_t min = fHistoDef->GetMin(), max = fHistoDef->GetMax();
    
  // list is created and named after the general
  // settings used for the contained histograms
  if (!histos) {
    AliError("NULL target list!");
    return;
  }
    
  // a general histogram is always added, 
  // which overrides the binning and collects everything
  fHisto[0] = new TH1D(histoName, histoTitle, nbins, min, max);
  histos->AddLast(fHisto[0]);
    
  // if requested a binning w.r. to some cut variable, histograms are added 
  // for that in this part of the method (one per each bin)
  Char_t hName[255];
  Char_t hTitle[255];
  if (fUseBins) {
    for (ibin = 0, i = 1; ibin < fBins.GetSize() - 1; ibin++, i++) {
      sprintf(hName, "%s[%.2f-%.2f]", histoName, fBins[ibin], fBins[ibin+1]);
      sprintf(hTitle, "%s [%.2f-%.2f]", histoTitle, fBins[ibin], fBins[ibin+1]);
      fHisto[i] = new TH1D(hName, hTitle, nbins, min, max);
      histos->AddLast(fHisto[i]);
    }
  }
}

//________________________________________________________________________________________
void AliRsnFunction::SetBinningCut
(AliRsnCut::EType type, Double_t min, Double_t max, Double_t step)
{
  //
  // Set fixed bins
  //
  
  fUseBins = kTRUE;

  Int_t i, nBins = (Int_t)((max - min) / step) + 1;
  fBinningCutType = type;
  fBins.Set(nBins);
  for (i = 0; i < nBins; i++) {
    fBins[i] = min + (Double_t)i * step;
  }
}

//________________________________________________________________________________________
void AliRsnFunction::SetBinningCut
(AliRsnCut::EType type, Int_t nbins, Double_t *bins)
{
  //
  // Set variable bins
  //

  fUseBins = kTRUE;

  Int_t i;
  fBinningCutType = type;
  fBins.Set(nbins);
  for (i = 0; i < nbins; i++) {
    fBins[i] = bins[i];
  }
}

//________________________________________________________________________________________
TString AliRsnFunction::GetFcnName()
{
  //
  // Return a string which names the function type
  //

  TString text("Undef");

  switch (fFcnType) {
  case kInvMass:
    text = "IM";
    break;
  case kInvMassMC:
    text = "IM_MC";
    break;
  case kResolution:
    text = "RES";
    break;
  case kPtSpectrum:
    text = "PT";
  default:
    AliError("Type not defined");
  }
  
  return text;
}

//________________________________________________________________________________________
TString AliRsnFunction::GetFcnTitle()
{
  //
  // Return a string which names the function type
  //
    
  TString text("Undef");

  switch (fFcnType) {
  case kInvMass:
    text = "Invariant mass";
    break;
  case kInvMassMC:
    text = "Invariant mass (MC)";
    break;
  case kResolution:
    text = "Resolution";
    break;
  case kPtSpectrum:
    text = "p_{#perp} distribution";
  default:
    AliError("Type not defined");
  }
  
  return text;
}

//________________________________________________________________________________________
Bool_t AliRsnFunction::Fill(AliRsnPairParticle *pair, AliRsnPairDef *ref, Double_t weight)
{
  //
  // Fillse the histogram with data contained in a defined pair.
  // This method must be overidden by an appropriate definition in each inheriting class.
  //

  Double_t value = FcnValue(pair, ref);
  if (fSkipOutsideInterval)
    {
      if (value < fHistoDef->GetMin()) return kFALSE;
      if (value > fHistoDef->GetMax()) return kFALSE;
    }
    
  // fill global histogram
  if (weight == 0.0) fHisto[0]->Fill(value);
  else fHisto[0]->Fill(value, weight);
    
  // if bins are allocated, find right one and fill it
  if (fUseBins) {
    Int_t i, ibin;
    for (ibin = 0, i = 1; ibin < fBins.GetSize() - 1; ibin++, i++) {
      if (!fHisto[i]) continue;
      fBinningCut.SetCutValues(fBinningCutType, fBins[ibin], fBins[ibin+1]);
      if (fBinningCut.IsSelected(AliRsnCut::kPair, pair)) {
        if (weight == 0.0) fHisto[i]->Fill(value);
        else fHisto[i]->Fill(value, weight);
        break;
      }
    }
  }
    
  return kTRUE;
}

//________________________________________________________________________________________
Double_t AliRsnFunction::FcnValue(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
  //
  // This method must be overridden in all inheritin functions.
  // It computes the value which must be used to fill the histogram.
  //

  switch (fFcnType) {
  case kInvMass:
    return pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
  case kInvMassMC:
    return pair->GetInvMassMC(ref->GetMass(0), ref->GetMass(1));
  case kResolution:
    return FcnResolution(pair, ref);
  case kPtSpectrum:
    return pair->GetPt();
  default:
    AliError("Type not defined");
  }
  
  return 0.0;
}

//________________________________________________________________________________________
inline Double_t AliRsnFunction::FcnResolution(AliRsnPairParticle *pair, AliRsnPairDef *ref)
{
  //
  // Invariant mass resolution (compared between reconstructed and montecarlo)
  //

  Double_t recInvMass = pair->GetInvMass(ref->GetMass(0), ref->GetMass(1));
  Double_t simInvMass = pair->GetInvMassMC(ref->GetMass(0), ref->GetMass(1));
    
  return (simInvMass - recInvMass) / simInvMass;
}
