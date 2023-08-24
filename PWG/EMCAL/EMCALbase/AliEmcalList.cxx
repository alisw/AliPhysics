/************************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <iostream>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "THnBase.h"
#include "AliLog.h"
#include "AliEmcalList.h"

#ifdef WITH_ROOUNFOLD
#include "RooUnfoldResponse.h"
#endif

ClassImp(AliEmcalList)

//________________________________________________________________________
AliEmcalList::AliEmcalList() : TList(), 
  fUseScaling(kFALSE),
  fNameXsec("fHistXsection"),
  fNameNTrials("fHistTrials")
{
}

//________________________________________________________________________
Long64_t AliEmcalList::Merge(TCollection *hlist)
{
  if(!hlist)
    return 0;

  AliInfoStream() << "Scaled merging for list " <<hlist->GetName() << " is " << (fUseScaling ? "" : "not ") << "activated." << std::endl;
  if(!fUseScaling)
    return TList::Merge(hlist);

  // #### Retrieve xsection and ntrials from histograms in this list
  // NOTE: they must be directly added to the AliEmcalList, not nested in sublists!
  TH1* xsection = static_cast<TH1*>(FindObject(fNameXsec.Data()));
  TH1* ntrials  = static_cast<TH1*>(FindObject(fNameNTrials.Data()));

  // #### If xsection or ntrials are not given, go to fatal
  if(!(xsection && ntrials))
    AliFatal("Scaled merging is active but AliEmcalList does not contain fHistXsection or fHistTrials. Do not activate scaling for those lists or include these histograms.");

  // #### Do the scaling only on the last level when we mix several pt hard bins
  // This is easy to find out checking the std histos in hlist
  Bool_t  isLastLevel = IsLastMergeLevel(hlist);
 
  // #### On last level, do the scaling
  if(isLastLevel)
  {
    AliInfoStream() << "===== LAST LEVEL OF MERGING =====" << std::endl;

    // Scale all histograms in this list
    Double_t scalingFactor = GetScalingFactor(xsection, ntrials);
    ScaleAllHistograms(this, scalingFactor);

    // Scale all histograms in the lists to be merged
    TIter listIterator(hlist);
    while (AliEmcalList* tmpList = static_cast<AliEmcalList*>(listIterator()))
    {
      xsection = static_cast<TH1*>(tmpList->FindObject(fNameXsec.Data()));
      ntrials  = static_cast<TH1*>(tmpList->FindObject(fNameNTrials.Data()));
      scalingFactor = GetScalingFactor(xsection, ntrials);
      ScaleAllHistograms(tmpList, scalingFactor);
    }
  }

  AliInfoStream() << "Merge() done." << std::endl;

  TList::Merge(hlist);
  return hlist->GetEntries() + 1;
}

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

    // Otherwise, scale object if scaling is supported
    if (!IsScalingSupported(listObject)){
      AliInfoStream() << "Scaling of objects of type " << listObject->IsA()->GetName() << " unsupported - object " << listObject->GetName() << " will not be scaled" << std::endl;
      continue;
    }

    // Don't scale profiles and histograms used for scaling
    TString histogram_class (listObject->ClassName());
    TString histogram_name (listObject->GetName());
    if ((histogram_name.Contains("fHistXsection") || histogram_name.Contains("fHistTrials") || histogram_name.Contains("fHistEvents") || histogram_name.Contains("fHistWeights")) && (!histogram_name.Contains("PtHard")))
    {
      AliInfoStream() << "Histogram " << listObject->GetName() << " will not be scaled, because a scaling histogram" << std::endl;
      continue;
    }
    if (histogram_class.Contains("TProfile"))
    {
      AliInfoStream() << "Histogram " << listObject->GetName() << " will not be scaled, because it is a TProfile" << std::endl;
      continue;
    }

    TH1 *histogram = dynamic_cast<TH1 *>(listObject);
    if(histogram) {
      // Handle TH1/TH2/TH3
      histogram->Sumw2();
      histogram->Scale(scalingFactor);
    } else if(listObject->InheritsFrom(THnBase::Class())){
      // Handle THn
      THnBase *histogramND = dynamic_cast<THnBase *>(listObject);
      histogramND->Sumw2();
      histogramND->Scale(scalingFactor);
    } 
#ifdef WITH_ROOUNFOLD
    else if(listObject->InheritsFrom(RooUnfoldResponse::Class())){
      RooUnfoldResponse *response = dynamic_cast<RooUnfoldResponse *>(listObject);
      if(auto measured = response->Hmeasured()) measured->Scale(scalingFactor);
      if(auto fakes = response->Hfakes()) fakes->Scale(scalingFactor);
      if(auto truth = response->Htruth()) truth->Scale(scalingFactor);
      if(auto responseND = response->Hresponse()) responseND->Scale(scalingFactor);
    }
#endif
    AliInfoStream() << "Histogram " << listObject->GetName() << " (" << histogram_class.Data() << ") was scaled..." << std::endl;

  }
}

//________________________________________________________________________
Double_t AliEmcalList::GetScalingFactor(const TH1* xsection, const TH1* ntrials) const
{
  Int_t binNumber = GetFilledBinNumber(xsection);
  if(!binNumber)
  {
    AliInfoStream() << "List already scaled or scaling invalid. Scaling factor = 1." << std::endl;
    return 1.0;
  }

  Double_t valNTRIALS = ntrials->GetBinContent(binNumber);
  Double_t valXSEC = xsection->GetBinContent(binNumber);
  Double_t scalingFactor = 0;
  if(valNTRIALS)
    scalingFactor = valXSEC/valNTRIALS;

  AliInfoStream() << "## Bin " << binNumber << ": trials=" << valNTRIALS << ", xsec=" << valXSEC << " -> scaling=" << scalingFactor << std::endl;
  return scalingFactor;
}

//________________________________________________________________________
Bool_t AliEmcalList::IsLastMergeLevel(const TCollection* collection) const 
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

//________________________________________________________________________
bool AliEmcalList::IsScalingSupported(const TObject *scaleobject) const {
  if(scaleobject->InheritsFrom(TH1::Class()) || scaleobject->InheritsFrom(THnBase::Class())) return true;
#ifdef WITH_ROOUNFOLD
  if(scaleobject->InheritsFrom(RooUnfoldResponse::Class())) return true;
#endif
  return false;
}

//________________________________________________________________________
Int_t AliEmcalList::GetFilledBinNumber(const TH1* hist) const
{
  AliInfoStream() << hist->GetName() << ": nbinsX=" << hist->GetNbinsX() << std::endl;;

  Int_t binFound = 0;
  for(Int_t i=1; i<=hist->GetNbinsX(); i++)
  {
    AliInfoStream() << hist->GetName() << ": bin=" << i << ", val=" << hist->GetBinContent(i) << std::endl;
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
    AliErrorStream() << "No bin filled in scaling histogram." << std::endl;

  // 0 if no bin found or more than one filled, otherwise returns bin number
  return binFound;
}
