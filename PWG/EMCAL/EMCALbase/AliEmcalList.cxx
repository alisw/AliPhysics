/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake                                                       *
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

#include "TList.h"
#include "TH1.h"
#include "THnBase.h"
#include "AliLog.h"
#include "AliEmcalList.h"

/// \cond CLASSIMP
ClassImp(AliEmcalList)
/// \endcond

//________________________________________________________________________
AliEmcalList::AliEmcalList() : TList(), fUseScaling(kFALSE)
{
  // constructor
}

/// Overridden ::Merge function
//________________________________________________________________________
Long64_t AliEmcalList::Merge(TCollection *hlist)
{
  if(!hlist)
    return 0;

  AliInfo(Form("Scaled merging for list %s is %sactivated.", hlist->GetName(),(fUseScaling) ? "" : "not "));
  if(!fUseScaling)
    return TList::Merge(hlist);

  // #### Retrieve xsection and ntrials from histograms in this list
  // NOTE: they must be directly added to the AliEmcalList, not nested in sublists!
  TH1* xsection = static_cast<TH1*>(FindObject("fHistXsection"));
  TH1* ntrials  = static_cast<TH1*>(FindObject("fHistTrials"));

  // #### If xsection or ntrials are not given, go to fatal
  if(!(xsection && ntrials))
    AliFatal("Scaled merging is active but AliEmcalList does not contain fHistXsection or fHistTrials. Do not activate scaling for those lists or include these histograms.");

  // #### Do the scaling only on the last level when we mix several pt hard bins
  // This is easy to find out checking the std histos in hlist
  Bool_t  isLastLevel = IsLastMergeLevel(hlist);
 
  // #### On last level, do the scaling
  if(isLastLevel)
  {
    AliInfo(Form("===== LAST LEVEL OF MERGING ====="));

    // Scale all histograms in this list
    Double_t scalingFactor = GetScalingFactor(xsection, ntrials);
    ScaleAllHistograms(this, scalingFactor);

    // Scale all histograms in the lists to be merged
    TIter listIterator(hlist);
    while (AliEmcalList* tmpList = static_cast<AliEmcalList*>(listIterator()))
    {
      xsection = static_cast<TH1*>(tmpList->FindObject("fHistXsection"));
      ntrials  = static_cast<TH1*>(tmpList->FindObject("fHistTrials"));
      scalingFactor = GetScalingFactor(xsection, ntrials);
      ScaleAllHistograms(tmpList, scalingFactor);
    }
  }

  AliInfo("Merge() done.");

  TList::Merge(hlist);
  return hlist->GetEntries() + 1;
}

/// Function that does the scaling of all histograms in hlist recursively
//________________________________________________________________________
void AliEmcalList::ScaleAllHistograms(TCollection *hlist, Double_t scalingFactor)
{
  TIter listIterator(hlist);
  while (TObject* listObject = listIterator())
  {
    // Recurse into nested all lists
    TCollection* sublist = dynamic_cast<TCollection*>(listObject);
    if(sublist)
    {
      ScaleAllHistograms(sublist, scalingFactor);
      continue;
    }

    // Otherwise, scale TH1-derived / THnBase-derived histograms
    if (!(listObject->InheritsFrom(TH1::Class()) || listObject->InheritsFrom(THnBase::Class())))
      continue;

    // Don't scale profiles and histograms used for scaling
    TString histogram_class (listObject->ClassName());
    if (!strcmp(listObject->GetName(), "fHistXsection") || !strcmp(listObject->GetName(), "fHistTrials"))
    {
      AliInfo(Form("Histogram %s will not be scaled, because a scaling histogram", listObject->GetName()));
      continue;
    }
    if (histogram_class.Contains("TProfile"))
    {
      AliInfo(Form("Histogram %s will not be scaled, because it is a TProfile", listObject->GetName()));
      continue;
    }

    TH1 *histogram = dynamic_cast<TH1 *>(listObject);
    if(histogram) {
      // Handle TH1/TH2/TH3
      histogram->Sumw2();
      histogram->Scale(scalingFactor);
    } else {
      // Handle THn
      THnBase *histogramND = dynamic_cast<THnBase *>(listObject);
      histogramND->Sumw2();
      histogramND->Scale(scalingFactor);
    }
    AliInfo(Form("Histogram %s (%s) was scaled...", listObject->GetName(), histogram_class.Data()));

  }
}

/// Helper function scaling factor
//________________________________________________________________________
Double_t AliEmcalList::GetScalingFactor(TH1* xsection, TH1* ntrials)
{
  Int_t binNumber = GetFilledBinNumber(xsection);
  if(!binNumber)
  {
    AliInfo("List already scaled or scaling invalid. Scaling factor = 1.");
    return 1.0;
  }

  Double_t valNTRIALS = ntrials->GetBinContent(binNumber);
  Double_t valXSEC = xsection->GetBinContent(binNumber);
  Double_t scalingFactor = 0;
  if(valNTRIALS)
    scalingFactor = valXSEC/valNTRIALS;

  AliInfo(Form("## Bin %i: trials=%f, xsec=%f -> scaling=%f", binNumber, valNTRIALS, valXSEC, scalingFactor));
  return scalingFactor;
}

/// Helper function to determine whether we are in last merge step
/// \param collection Collection of AliEmcalList objects
//________________________________________________________________________
Bool_t AliEmcalList::IsLastMergeLevel(TCollection* collection)
{
  // Get the pt hard bin number that is filled in this object
  TH1* xsection = static_cast<TH1*>(FindObject("fHistXsection"));
  Int_t binNumberThis = GetFilledBinNumber(xsection);
  
  TIter listIterator(collection);
  while (AliEmcalList* tmpList = static_cast<AliEmcalList*>(listIterator()))
  {
    // For every list in the collection, test if they are from different pt hard bins
    xsection = static_cast<TH1*>(tmpList->FindObject("fHistXsection"));
    if(binNumberThis != GetFilledBinNumber(xsection))
      return kTRUE;
  }
  return kFALSE;
}

/// Helper function that returns the bin in a TH1 that is filled
/// \param hist  Histogram
/// \return bin number that is filled. If no or more than one bin is filled, 0.
//________________________________________________________________________
Int_t AliEmcalList::GetFilledBinNumber(TH1* hist)
{
  AliInfo(Form("%s: nbinsX=%i", hist->GetName(), hist->GetNbinsX()));

  Int_t binFound = 0;
  for(Int_t i=1; i<=hist->GetNbinsX(); i++)
  {
    AliInfo(Form("%s: bin=%i, val=%f", hist->GetName(), i, hist->GetBinContent(i)));
    if(hist->GetBinContent(i))
    {
      if(!binFound)
        binFound = i;
      else
      {
        return 0; // more than one bin filled (e.g. when merging is partly done)
      }
    }
  }

  if(!binFound)
    AliError("No bin filled in scaling histogram.");

  // 0 if no bin found or more than one filled, otherwise returns bin number
  return binFound;
}
