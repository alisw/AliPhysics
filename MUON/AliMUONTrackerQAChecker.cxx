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

#include "AliMUONTrackerQAChecker.h"

/// \class AliMUONTrackerQAChecker
///
/// Implementation of AliQACheckerBase for MCH and MTR
///
/// For the moment we only implement the checking of raw data QA for the tracker
/// by looking at the occupancy at the bus patch level.
///
/// \author Laurent Aphecetche, Subatech

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONQAIndices.h"
#include "AliMUONRecoParam.h"
#include "AliMpBusPatch.h"
#include "AliMpDDLStore.h"
#include "AliQAv1.h"
#include "AliQAv1.h"
#include "Riostream.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TVirtualPad.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackerQAChecker)
/// \endcond

namespace {
  
  //___________________________________________________________________________  
  int trim(Int_t n, 
           Double_t* x,
           Double_t alpha,
           Double_t& tmean,
           Double_t& tvar,
           Double_t& min,
           Double_t& max)
  {
    //
    // Calculates the trimmed (tmean) mean
    // of a sample (x) and estimates the variance (tvar)
    // of that mean.
    //
    
    // First check input parameters
    
    // number of observations
    if ( n < 2 )
    {
      return -1;
    }
    
    if ( alpha < 0 || alpha >= 0.5 )
      // proportion of observations
      // to be trimmed at each end of the sorted sample
    {
      return -2;
    }
    
    // Input parameters are good. Let's move on.
    
    // Insure we use a sample sorted into ascending order.
    
    Int_t* indices = new Int_t[n];
    
    TMath::Sort(n,x,indices,kFALSE);
    
    Double_t* sx = new Double_t[n];
    
    for ( Int_t i = 0; i < n; ++i )
    {
      sx[i] = x[indices[i]];
    }
    delete[] indices;
    
    
    // Number of observations trimmed at each end.
    
    Int_t k = TMath::FloorNint(alpha * n);
    
    double sum = 0.0;
    
    for ( Int_t i = k; i < n - k ; ++i )
    {
      sum += sx[i];
    }
    
    tmean = sum / ( n - 2 * k );
    
    double t2 = 0.0;
    
    for ( Int_t i = k; i < n - k; ++i )
    {
      t2 += (sx[i] - tmean) * (sx[i] - tmean);
    }
    
    tvar = (
            t2 +
            k * (sx[k] - tmean) * (sx[k] - tmean) +
            k * (sx[n - k - 1] - tmean) * (sx[n - k - 1] - tmean)
            ) / (n * n);
    
    // get the min and max for the non-rejected values
    min = DBL_MAX;
    max = 0.0;
    
    for ( Int_t i = k; i < n-k; ++i ) 
    {
      min = TMath::Min(min,sx[i]);
      max = TMath::Max(max,sx[i]);
    }
    
    delete[] sx;
    
    return 0;
  }
}

//__________________________________________________________________
AliMUONTrackerQAChecker::AliMUONTrackerQAChecker() : AliMUONVQAChecker()
{
	/// ctor
}          

//__________________________________________________________________
AliMUONTrackerQAChecker::~AliMUONTrackerQAChecker() 
{
	/// dtor
}

//__________________________________________________________________
AliMUONTrackerQAChecker::AliMUONTrackerQAChecker(const AliMUONTrackerQAChecker& qac) : 
    AliMUONVQAChecker(qac) 
{
	/// copy ctor 
}   

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode*
AliMUONTrackerQAChecker::CheckRecPoints(TObjArray ** list, const AliMUONRecoParam* /*recoParam*/)
{
  /// Check rec points
  /// Very binary check for the moment. 
  
  AliCodeTimerAuto("",0);
  
  AliMUONVQAChecker::ECheckCode * rv = new AliMUONVQAChecker::ECheckCode[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] =  AliMUONVQAChecker::kInfo; 
  
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
  {
    TH1* h = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerNumberOfClustersPerDE,AliRecoParam::ConvertIndex(specie));

    if ( !h ) rv[specie] = AliMUONVQAChecker::kWarning; // only a warning if histo not found, in order not to kill anything because QA failed...
  
    else if ( h->GetMean() == 0.0 ) rv[specie] =  MarkHisto(*h,AliMUONVQAChecker::kFatal);
  }
  return rv;
}

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode 
AliMUONTrackerQAChecker::MarkHisto(TH1& histo, AliMUONVQAChecker::ECheckCode value) const
{
  /// Mark histo as originator of some QA error/warning
  
  if ( value != AliMUONVQAChecker::kInfo )
  {
    histo.SetBit(AliQAv1::GetQABit());
  }
  
  return value;
}

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode*
AliMUONTrackerQAChecker::CheckESD(TObjArray ** list, const AliMUONRecoParam* /*recoParam*/)
{
  /// Check ESD
  
  AliCodeTimerAuto("",0);
  
  AliMUONVQAChecker::ECheckCode * rv = new AliMUONVQAChecker::ECheckCode[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = AliMUONVQAChecker::kInfo; 
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    
    TH1* h = AliQAv1::GetData(list,AliMUONQAIndices::kESDnTracks,AliRecoParam::ConvertIndex(specie));
  
    if (!h) rv[specie] = AliMUONVQAChecker::kWarning;
  
    else if ( h->GetMean() == 0.0 ) rv[specie] =  MarkHisto(*h,AliMUONVQAChecker::kFatal); // no track -> fatal
  
    h = AliQAv1::GetData(list,AliMUONQAIndices::kESDMatchTrig,AliRecoParam::ConvertIndex(specie));
  
    if (!h) rv[specie] = AliMUONVQAChecker::kWarning;
  
    else if (h->GetMean() == 0.0 ) rv[specie] = MarkHisto(*h,AliMUONVQAChecker::kError); // no trigger matching -> error
  }
  return rv;
}

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode*
AliMUONTrackerQAChecker::CheckRaws(TObjArray ** list, const AliMUONRecoParam* recoParam)
{
  /// Check raws

  AliCodeTimerAuto("",0);
  
  if (!recoParam) return 0x0;
  
  AliMUONVQAChecker::ECheckCode * rv = new AliMUONVQAChecker::ECheckCode[AliRecoParam::kNSpecies] ; 

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
  {
    rv[specie] = AliMUONVQAChecker::kInfo; 
  }
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
  {
    TH1* hbp = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerBusPatchOccupancy,AliRecoParam::ConvertIndex(specie));    

    TH1* hnpads = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerBusPatchNofPads,AliRecoParam::ConvertIndex(specie));

    TH1* hbpconfig = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerBusPatchConfig,AliRecoParam::ConvertIndex(specie));

    TH1* hnevents = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerNofRawEventSeen,AliRecoParam::ConvertIndex(specie));

    TH1* hddl = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerDDLOccupancy,AliRecoParam::ConvertIndex(specie));
    
    if ( !hbp || !hnpads || !hnevents || !hddl ) 
    {
      continue;
    }
    
    Int_t nevents = TMath::Nint(hnevents->GetBinContent(1));
        
    rv[specie] = BeautifyTrackerBusPatchOccupancy(*hddl,*hbp,hbpconfig,*hnpads,nevents,*recoParam);    
  }
  
  return rv;
}

//____________________________________________________________________________ 
AliMUONVQAChecker::ECheckCode
AliMUONTrackerQAChecker::BeautifyTrackerBusPatchOccupancy(TH1& hddl,
                                                          TH1& hbp, 
                                                          const TH1* hbuspatchconfig, 
                                                          const TH1& hnpads, 
                                                          Int_t nevents,
                                                          const AliMUONRecoParam& recoParam)
{
  /// Put labels, limits and so on on the TrackerBusPatchOccupancy histogram
  /// hbuspatchconfig and hbp must have the same bin definitions
  
  if ( hbuspatchconfig )
  {
    if ( hbp.GetNbinsX() != hbuspatchconfig->GetNbinsX() ||
        hbp.GetXaxis()->GetXmin() != hbuspatchconfig->GetXaxis()->GetXmin() ||
        hbp.GetXaxis()->GetXmax() != hbuspatchconfig->GetXaxis()->GetXmax() )
    {
      AliError("hbp and hbuspatchconfig histograms are not compatible !");
      return AliMUONVQAChecker::kFatal;
    }
  }
  
  hbp.SetXTitle("Absolute Bus Patch Id");
  hbp.SetYTitle("Occupancy (percent)");
  hbp.SetStats(kFALSE);
  
  Double_t xmin = hbp.GetXaxis()->GetXmin();
  Double_t xmax = hbp.GetXaxis()->GetXmax();
  
  Double_t occMax(0.1); // 0.1% y-limit for the plot
  Double_t maxToleratedOccupancy(recoParam.BuspatchOccupancyHighLimit()*100.0); 
  Double_t minToleratedOccupancy(recoParam.BuspatchOccupancyLowLimit()*100.0);   
  TLine* line1 = new TLine(xmin,maxToleratedOccupancy,xmax,maxToleratedOccupancy);
  line1->SetLineColor(1);
  line1->SetLineWidth(1);

  TLine* line2 = new TLine(xmin,minToleratedOccupancy,xmax,minToleratedOccupancy);
  line2->SetLineColor(1);
  line2->SetLineWidth(1);
  
  hbp.GetListOfFunctions()->Add(line1);
  hbp.GetListOfFunctions()->Add(line2);
    
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp(0x0);

  Int_t nMissingPads(0);
  Int_t nPads(0);
  Int_t nBusPatches(0);
  Int_t nMissingBusPatches(0);
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
  {
    Int_t bin = hbp.FindBin(bp->GetId());
    Int_t n = TMath::Nint(hnpads.GetBinContent(bin));
    
    ++nBusPatches;
    
    nPads += n;
    
    if ( hbp.GetBinContent(bin) <= 0 ) 
    {
      nMissingPads += n;
      ++nMissingBusPatches;
    }
  }
  
  next.Reset();
  
  Int_t ok(-1);
  Int_t n(0);
  Int_t nBusPatchesAboveLimit(0);
  Int_t nBusPatchesBelowLimit(0);
  Double_t alpha(0.1); // trim 10% of data
  Double_t tmean(0.0),tvar(0.0);
  Double_t ymin(0.0),ymax(0.0);
  AliMUONVQAChecker::ECheckCode rv(AliMUONVQAChecker::kFatal); // default value = serious problem
  
  if ( nBusPatches ) 
  {
    Double_t* x = new Double_t[nBusPatches];
    
    while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
    {
      Int_t bin = hbp.FindBin(bp->GetId());
      if ( hbp.GetBinContent(bin) > 0 )
      {
        x[n] = hbp.GetBinContent(bin);
        ++n;
      }
      if ( hbp.GetBinContent(bin) > maxToleratedOccupancy )
      {
        ++nBusPatchesAboveLimit;
      }
      if ( hbp.GetBinContent(bin) < minToleratedOccupancy )
      {
        // check whether this buspatch has a reason to be absent (only valid
        // if we got the config, otherwise we cannot do the test)
        if ( hbuspatchconfig && hbuspatchconfig->GetBinContent(bin) > 0 )
        {
          // should be there, so it's an error
          ++nBusPatchesBelowLimit;
        }
      }
    }
    
    // computed the truncated mean of the occupancy values, in order to get a 
    // reasonable y-range for the histogram (without giant peaks to the roof 
    // for misbehaving buspatches).
    ok = trim(n,x,alpha,tmean,tvar,ymin,ymax);
    
    delete[] x;
  }
  
  if ( ok < 0 ) 
  {
    ymax = occMax;
  }
  else
  {
    ymax = TMath::Max(ymax,occMax);
  }
  
  hbp.SetMaximum(ymax*1.4);
  
  TPaveText* text = new TPaveText(0.30,0.50,0.99,0.99,"NDC");
  
  text->AddText(Form("MCH RUN %d - %d events",AliCDBManager::Instance()->GetRun(),nevents));
  
  if ( ok < 0 ) 
  {
    text->AddText("Could not compute truncated mean. Not enough events ?");
    text->AddText(Form("nBusPatches=%d n=%d",nBusPatches,n));
  }
  else if (!nPads || !nBusPatches)
  {
    text->AddText(Form("Could not get the total number of pads (%d) or total number of buspatches (%d). ERROR !!!",
                       nPads,nBusPatches));
  }
  else
  {
    Float_t missingPadFraction = nMissingPads*100.0/nPads;
    Float_t missingBusPatchFraction = nMissingBusPatches*100.0/nBusPatches;
    Float_t aboveLimitFraction = nBusPatchesAboveLimit*100.0/nBusPatches;
    Float_t belowLimitFraction = nBusPatchesBelowLimit*100.0/nBusPatches;
    
    text->AddText(Form("%5.2f %% of missing buspatches (%d out of %d)",missingBusPatchFraction,nMissingBusPatches,nBusPatches));
    text->AddText(Form("%5.2f %% of missing pads (%d out of %d)",missingPadFraction,nMissingPads,nPads));
    text->AddText(Form("%5.2f %% bus patches above the %5.2f %% limit",aboveLimitFraction,maxToleratedOccupancy));
    text->AddText(Form("%5.2f %% bus patches below the %e %% limit",belowLimitFraction,minToleratedOccupancy));
    text->AddText(Form("Bus patch mean occupancy (truncated at %2d %%) is %7.2f %%",(Int_t)(alpha*100),tmean));
    
    if ( missingPadFraction >= 100.0 ) 
    {
      rv = AliMUONVQAChecker::kFatal;       
    }
    
    else if ( missingPadFraction > recoParam.MissingPadFractionLimit()*100.0 || 
             aboveLimitFraction > recoParam.FractionOfBuspatchOutsideOccupancyLimit()*100.0 || 
             belowLimitFraction > recoParam.FractionOfBuspatchOutsideOccupancyLimit()*100.0 ) 
    {
      rv = AliMUONVQAChecker::kError;
    }
    else
    {
      rv = AliMUONVQAChecker::kInfo;
    }
  }
  
  hbp.GetListOfFunctions()->Add(text);
  
  if ( rv == AliMUONVQAChecker::kInfo ) 
  {
    text->SetFillColor(3); // green = INFO
  }
  else if ( rv == AliMUONVQAChecker::kWarning ) 
  {
    text->SetFillColor(5); // yellow = WARNING
  }
  else if ( rv ==  AliMUONVQAChecker::kFatal) 
  {
    text->SetFillColor(2); // red = FATAL
  }
  else 
  {
    text->SetFillColor(6); // pink = ERROR
  }
  
  /// Make as well a version for DDL occupancy, that'll be used by the shifter
  
  hddl.GetListOfFunctions()->Add(text->Clone());
  
  Bool_t aboveOnePercent(kFALSE);
  Bool_t aboveTwoPercent(kFALSE);
  
  for ( Int_t i = 1; i <= hddl.GetXaxis()->GetNbins(); ++i )
  {
    Double_t b = hddl.GetBinContent(i);
    if ( b > 1.0 ) aboveOnePercent = kTRUE;
    if ( b > 2.0 ) aboveTwoPercent = kTRUE;
    
  }
  
  hddl.SetMaximum(2);
  hddl.SetFillStyle(0);
  if ( aboveOnePercent ) 
  {
    hddl.SetFillStyle(1001);
    hddl.SetFillColor(kOrange);    
  }
  if ( aboveTwoPercent ) 
  {
    hddl.SetFillStyle(1001);
    hddl.SetFillColor(kRed);
  }
  hddl.SetLineWidth(3);
  hddl.SetStats(kFALSE);
  
  return rv;
}
