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
#include "Riostream.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TLatex.h"

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

  //___________________________________________________________________________  
  Int_t GetColorFromCheckCode(AliMUONVQAChecker::ECheckCode code)
  {
    if ( code == AliMUONVQAChecker::kInfo ) return AliMUONVQAChecker::kInfoColor;
    else if ( code == AliMUONVQAChecker::kWarning ) return AliMUONVQAChecker::kWarningColor;
    else if ( code ==  AliMUONVQAChecker::kFatal) return AliMUONVQAChecker::kFatalColor;
    else return AliMUONVQAChecker::kErrorColor;
  }
  
  const char* NOTENOUGHEVENTMESSAGE = "Not enough event to judge. Please wait a bit";
  const char* NOTIFYEXPERTMESSAGE = "PLEASE NOTIFY EXPERT !";
  const char* ALLISFINEMESSAGE = "All is fine. Just enjoy.";
  
  const int NOTENOUGHEVENTLIMIT = 50;
  
  //___________________________________________________________________________  
  void SetupHisto(Int_t neventsseen, Int_t neventsused, const TObjArray& messages, TH1& histo, AliMUONVQAChecker::ECheckCode code, const char* extraopt="")
  {
    Bool_t allIsFine(kFALSE);
    
    if ( code == AliMUONVQAChecker::kInfo ) 
    {
      allIsFine = kTRUE;
    }    
    
    Double_t y1 = 0.99 - (messages.GetLast()+(allIsFine?1:0)+2)*0.075;
    
    TPaveText* text = new TPaveText(0.5,y1,0.99,0.99,"NDC");
    
    text->AddText(Form("MCH RUN %d - %d events seen - %d events used",AliCDBManager::Instance()->GetRun(),neventsseen,neventsused));
    
    TIter next(&messages);
    TObjString* str;
    
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      text->AddText(str->String());
    }
    
    if ( allIsFine )
    {
      text->AddText(ALLISFINEMESSAGE);
    }
                    
    text->SetFillColor(GetColorFromCheckCode(code));
                      
    Int_t color = GetColorFromCheckCode(code);
    
    histo.SetFillStyle(1001);
    histo.SetFillColor(color);    
    TString sopt("BAR");
    sopt += extraopt;
    histo.SetOption(sopt.Data());
    
    histo.SetTitle(kFALSE);
    
    TList* lstF = histo.GetListOfFunctions();
    TObject* title = lstF->FindObject("title");
    if (title) delete lstF->Remove(title);
    // 
    lstF->Add(text);  
    
  }

  //___________________________________________________________________________  
  void ShowTrueValue(TH1& hrostatusnorm, Int_t v)
  {
    // For a bar > 100% we put the text inside the bar (as TEXT45 option is putting above which
    // would go offscreen)
    
    Int_t bin = hrostatusnorm.FindBin(1.0*v);
    
    Double_t value = hrostatusnorm.GetBinContent(bin);
    
    if ( value > 100.0 ) 
    {
      TLatex* l = new TLatex(hrostatusnorm.GetBinCenter(bin),50,Form("%g",value));
      l->SetTextFont(gStyle->GetTextFont());
      l->SetTextAlign(22);
      l->SetTextAngle(45);
      l->SetTextSize(0.02*hrostatusnorm.GetMarkerSize());
      hrostatusnorm.GetListOfFunctions()->Add(l);
    }
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
    TH1* hneventsseen = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerNofPhysicsEventsSeen,AliRecoParam::ConvertIndex(specie));
    TH1* hneventsused = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerNofGoodPhysicsEventsUsed,AliRecoParam::ConvertIndex(specie));

    if (!hneventsseen || !hneventsused ) continue;
    
    Int_t neventsseen = TMath::Nint(hneventsseen->GetBinContent(1));
    Int_t neventsused = TMath::Nint(hneventsused->GetBinContent(1));
        
    AliMUONVQAChecker::ECheckCode c1 = AliMUONVQAChecker::kInfo;
    AliMUONVQAChecker::ECheckCode c2 = AliMUONVQAChecker::kInfo;
    AliMUONVQAChecker::ECheckCode c3 = AliMUONVQAChecker::kInfo;
    
    TH1* hbp = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerBusPatchOccupancy,AliRecoParam::ConvertIndex(specie));    
    TH1* hbpconfig = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerBusPatchConfig,AliRecoParam::ConvertIndex(specie));
    TH1* hddl = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerDDLOccupancy,AliRecoParam::ConvertIndex(specie));    

    if ( hbp && hbpconfig && hddl ) 
    {
      c1 = BeautifyOccupancyHistograms(*hddl,*hbp,hbpconfig,neventsseen,neventsused,*recoParam);  
    }
    else
    {
      AliError(Form("Could not BeautifyOccupancyHistograms : hddl=%p hbpconfig=%p hbp=%p",hddl,hbpconfig,hbp));
    }
    
    TH1* hbuspatchtokenerrors = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerBusPatchTokenLostErrors,AliRecoParam::ConvertIndex(specie));
    TH1* hrostatus = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerReadoutStatus,AliRecoParam::ConvertIndex(specie));    
    TH1* hrostatusnorm = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerReadoutStatusPerEvent,AliRecoParam::ConvertIndex(specie));
    
    if ( hrostatus && hrostatusnorm && hbuspatchtokenerrors )
    {
      c2 = BeautifyReadoutHistograms(*hrostatus,*hrostatusnorm,*hbuspatchtokenerrors,
                                     neventsseen,neventsused,*recoParam);
    }
    else
    {
      AliError(Form("Could not BeautifyReadoutHistograms : hrostatus=%p hrostatusnorm=%p hbuspatchtokenerrors=%p",
                    hrostatus,hrostatusnorm,hbuspatchtokenerrors));
    }
    
    
    TH1* heventsize = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerDDLEventSize,AliRecoParam::ConvertIndex(specie));    
    TH1* heventsizeperevent = AliQAv1::GetData(list,AliMUONQAIndices::kTrackerDDLEventSizePerEvent,AliRecoParam::ConvertIndex(specie));
    
    if ( heventsize && heventsizeperevent ) 
    {
      c3 = BeautifyEventsizeHistograms(*heventsize,*heventsizeperevent,
                                       neventsseen,neventsused,*recoParam);
    }
    else
    {
      AliError(Form("Could not BeautifyEventsizeHistograms heventsize=%p heventsizeperevent=%p",heventsize,heventsizeperevent));
    }
    
    if ( c1 == AliMUONVQAChecker::kFatal || c2 == AliMUONVQAChecker::kFatal || c3 == AliMUONVQAChecker::kFatal )
    {
      rv[specie] = AliMUONVQAChecker::kFatal;
    }
    else if ( c1 == AliMUONVQAChecker::kError || c2 == AliMUONVQAChecker::kError || c3 == AliMUONVQAChecker::kError )
    {
      rv[specie] = AliMUONVQAChecker::kError;
    }
    else if ( c1 == AliMUONVQAChecker::kWarning || c2 == AliMUONVQAChecker::kWarning || c3 == AliMUONVQAChecker::kWarning )
    {
      rv[specie] = AliMUONVQAChecker::kWarning;
    }
    else
    {
      rv[specie] = AliMUONVQAChecker::kInfo;
    }
  }
  
  return rv;
}

//____________________________________________________________________________ 
AliMUONVQAChecker::ECheckCode
AliMUONTrackerQAChecker::BeautifyOccupancyHistograms(TH1& hddl,
                                                     TH1& hbp, 
                                                     const TH1* hbuspatchconfig, 
                                                     Int_t neventsseen,
                                                     Int_t neventsused,
                                                     const AliMUONRecoParam& recoParam)
{
  /// Put labels, limits and so on on the TrackerBusPatchOccupancy histograms
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
  
  hddl.SetStats(kFALSE);
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
  
  TList* lstF = hbp.GetListOfFunctions();
  if (lstF) {
    TObject *stats = lstF->FindObject("stats");
    lstF->Remove(stats);
    TObject *obj;
    while ((obj = lstF->First())) { while(lstF->Remove(obj)) { } delete obj; }
    if (stats) lstF->Add(stats);
  } 
  
  hbp.GetListOfFunctions()->Add(line1);
  hbp.GetListOfFunctions()->Add(line2);
    
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp(0x0);

  Int_t nBusPatches(0);
  Int_t nMissingBusPatches(0);
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
  {
    ++nBusPatches;
    
    Int_t bin = hbp.FindBin(bp->GetId());
    
    if ( hbp.GetBinContent(bin) <= 0.0 ) 
    {
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
  
  AliMUONVQAChecker::ECheckCode rv(AliMUONVQAChecker::kInfo);

  TObjArray messages;
  messages.SetOwner(kTRUE);
  
  if ( neventsseen < NOTENOUGHEVENTLIMIT ) 
  {
    messages.Add(new TObjString(NOTENOUGHEVENTMESSAGE));
  }
  else
  {
    if ( ok < 0 ) 
    {
      messages.Add(new TObjString("Could not compute truncated mean. Not enough events ?"));
      
      if ( neventsused < TMath::Nint(0.1*neventsseen) )
      {
        messages.Add(new TObjString(Form("We've actually seen %d events, but used only %d !",neventsseen,neventsused)));
        messages.Add(new TObjString("For a normal physics run, this is highly suspect !"));
        messages.Add(new TObjString(NOTIFYEXPERTMESSAGE));
        rv = AliMUONVQAChecker::kFatal;
      }
    }
    else if (!nBusPatches)
    {
      messages.Add(new TObjString(Form("Could not get the total number of buspatches (%d). ERROR !!!",
                         nBusPatches)));
      messages.Add(new TObjString("Please check with expert !"));
      rv = AliMUONVQAChecker::kFatal;
    }
    else
    {
      Float_t missingBusPatchFraction = nMissingBusPatches*100.0/nBusPatches;
      Float_t aboveLimitFraction = nBusPatchesAboveLimit*100.0/nBusPatches;
      Float_t belowLimitFraction = nBusPatchesBelowLimit*100.0/nBusPatches;
      
      messages.Add(new TObjString(Form("%5.2f %% of missing buspatches (%d out of %d)",missingBusPatchFraction,nMissingBusPatches,nBusPatches)));
      messages.Add(new TObjString(Form("%5.2f %% bus patches above the %5.2f %% limit",aboveLimitFraction,maxToleratedOccupancy)));
      messages.Add(new TObjString(Form("%5.2f %% bus patches below the %e %% limit",belowLimitFraction,minToleratedOccupancy)));
      messages.Add(new TObjString(Form("Bus patch mean occupancy (truncated at %2d %%) is %7.2f %%",(Int_t)(alpha*100),tmean)));
      
      if ( aboveLimitFraction > recoParam.FractionOfBuspatchOutsideOccupancyLimit()*100.0 || 
          belowLimitFraction > recoParam.FractionOfBuspatchOutsideOccupancyLimit()*100.0 ) 
      {
        rv = AliMUONVQAChecker::kError;
      }
      else
      {
        rv = AliMUONVQAChecker::kInfo;
      }
    }
  }
  
  SetupHisto(neventsseen,neventsused,messages,hbp,rv);
  
  /// Make as well a version for DDL occupancy, that'll be used by the shifter
  SetupHisto(neventsseen,neventsused,messages,hddl,rv);
  
  
  hddl.SetStats(kFALSE);
  
  return rv;
}

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode 
AliMUONTrackerQAChecker::BeautifyReadoutHistograms(TH1& hrostatus,
                                                   TH1& hrostatusnorm,
                                                   const TH1& hbuspatchtokenerrors,
                                                   Int_t neventsseen,
                                                   Int_t neventsused,
                                                   const AliMUONRecoParam& recoParam)
{
  /// Normalize and put some text on the readout error histogram
  /// Note in particular the treatment of tokenlost errors !
  
  hrostatusnorm.Reset();
  
  AliMUONVQAChecker::ECheckCode rv(AliMUONVQAChecker::kInfo);
  
  TObjArray messages;
  messages.SetOwner(kTRUE);
  
  if ( neventsseen < NOTENOUGHEVENTLIMIT )
  {
    messages.Add(new TObjString(NOTENOUGHEVENTMESSAGE));
  }
  else
  {
    hrostatusnorm.Add(&hrostatus,100.0/neventsseen);
    hrostatusnorm.SetOption("TEXT45"); // note : cannot use HIST option, otherwise the associated function (i.e. the tpavetext) is not drawn...
    hrostatusnorm.SetMinimum(0.0);
    hrostatusnorm.SetMaximum(100.0);
    
    Double_t tokenLost = hrostatusnorm.GetBinContent(hrostatusnorm.FindBin(1.0*AliMUONQAIndices::kTrackerRawNofTokenLostErrors));
    
    if ( tokenLost > recoParam.TokenLostLimit() )
    {
      rv = AliMUONVQAChecker::kError;
      
      messages.Add(new TObjString("There are some token lost errors !"));
      messages.Add(new TObjString("PLEASE CHECK THE BUSY TIME FOR MUON !"));
      messages.Add(new TObjString("If above 5 ms please have the MUON expert"));
      messages.Add(new TObjString("check the following bus patches :"));
      
      for ( Int_t i = 1; i <= hbuspatchtokenerrors.GetNbinsX(); ++i ) 
      {
        if ( hbuspatchtokenerrors.GetBinContent(i) > 0 ) 
        {
          messages.Add(new TObjString(Form("BP %4d",i)));
        }
      }
    }
    
    if ( hrostatusnorm.GetBinContent(hrostatusnorm.FindBin(AliMUONQAIndices::kTrackerRawNofEmptyEvents)) > 25.0 )
    {
      messages.Add(new TObjString("Too many empty events !"));
      messages.Add(new TObjString(NOTIFYEXPERTMESSAGE));
      rv = AliMUONVQAChecker::kFatal;      
    }    
  }
   
  SetupHisto(neventsseen,neventsused,messages,hrostatusnorm,rv,"text45");
  
  ShowTrueValue(hrostatusnorm,AliMUONQAIndices::kTrackerRawNofTokenLostErrors);
  ShowTrueValue(hrostatusnorm,AliMUONQAIndices::kTrackerRawNofParityErrors);
  ShowTrueValue(hrostatusnorm,AliMUONQAIndices::kTrackerRawNofPaddingErrors);
  ShowTrueValue(hrostatusnorm,AliMUONQAIndices::kTrackerRawNofGlitchErrors);
  
  return rv;
}

//______________________________________________________________________________
AliMUONVQAChecker::ECheckCode 
AliMUONTrackerQAChecker::BeautifyEventsizeHistograms(TH1& heventsize,
                                                     TH1& heventsizeperevent,
                                                     Int_t neventsseen,
                                                     Int_t neventsused,
                                                     const AliMUONRecoParam& recoParam)
{
  /// Normalize and put some text on the event size histogram

  AliMUONVQAChecker::ECheckCode rv(AliMUONVQAChecker::kInfo);

  heventsizeperevent.Reset();

  TObjArray messages;
  messages.SetOwner(kTRUE);
  
  if ( neventsseen < NOTENOUGHEVENTLIMIT )
  {
    messages.Add(new TObjString(NOTENOUGHEVENTMESSAGE));
  }
  else
  {
    heventsizeperevent.Add(&heventsize,1.0/neventsseen/1024.0);
    heventsizeperevent.SetMinimum(0);
    
    Double_t totalEventSizePerEvent = heventsizeperevent.Integral();
    
    TString msg;
    TString action;
    
    if ( totalEventSizePerEvent > recoParam.EventSizeHardLimit() ) 
    {
      rv = AliMUONVQAChecker::kError;
      msg = "That is really too high.";
      action = NOTIFYEXPERTMESSAGE;
    }
    else if ( totalEventSizePerEvent > recoParam.EventSizeSoftLimit() ) 
    {
      msg = "That is a bit high.";
      action = "Please keep an eye on it.";
      rv = AliMUONVQAChecker::kWarning;
    }
    else 
    {
      rv = AliMUONVQAChecker::kInfo;
    }
    
    messages.Add(new TObjString(Form("<MCH event size> %5.1f KB/event\n",totalEventSizePerEvent)));
    if (msg.Length()>0)
    {
      messages.Add(new TObjString(msg));
      messages.Add(new TObjString(action));
    }
  }
  
  SetupHisto(neventsseen,neventsused,messages,heventsizeperevent,rv);
  
  return rv;
}



