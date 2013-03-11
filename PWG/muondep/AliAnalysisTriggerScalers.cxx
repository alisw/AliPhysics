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

//
// AliAnalysisTriggerScalers : utility class to play with OCDB scalers.
//
// Can use it to e.g. :
//
// - get the integrated luminosity for a given period
//   (see IntegratedLuminosity method)
//
// - plot the trigger rates (see PlotTriggerEvolution)
//
// Please note that this class is doing an OCDB "scan" (for as many run
// as you put in the SetRunList method), so it can be really long
// if you're using the raw:// OCDB. If you're planning on using this
// class for anything more than a brief test, you'd better think
// or making a local copy of the OCDB objects required by this class :
//
// GRP/CTP/Config
// GRP/CTP/Scalers
//
// and also (to get the fill and period ranges, mainly for drawing)
//
// GRP/GRP/Data
// GRP/GRP/LHCData
//
//
// Please note also that this class evolved from a set a quick macros to
// follow the luminosity and trigger rates during data taking to this stage.
//
// Now (Feb 2013) would probably be a good time to regroup a bit, think about it and
// make it more general, more robust and just less "had-hoc"...
//
// Ideas for improvement :
//
// - make sure the "what" parameter mean the same thing in all methods
//   and so can take the same values in all methods
//
// - get the lumi trigger name and cross-section from e.g. OADB instead of
//   hard-coding them ?
//
// - organize the class so that the CTP Scalers and Config can be fetched
//   from elsewhere (e.g. from a run based container available in some
//   distant future in the ESDs/AODs ?)
//
// - robustify the PAC fraction computation
//   
//
// L. Aphecetche (Subatech)
//

#include "AliAnalysisTriggerScalers.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGRPObject.h"
#include "AliLog.h"
#include "AliTimeStamp.h"
#include "AliTriggerBCMask.h"
#include "AliTriggerClass.h"
#include "AliTriggerCluster.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerDescriptor.h"
#include "AliTriggerInput.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersESD.h"
#include "AliTriggerScalersRecord.h"
#include "Riostream.h"
#include "TDatime.h"
#include "TError.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <set>
#include <map>
#include "TStyle.h"
#include "AliLHCData.h"
#include "TMath.h"
#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TText.h"

ClassImp(AliAnalysisTriggerScalers)

namespace {
  
  //______________________________________________________________________________
  UChar_t GetIndex(ULong64_t mask)
  {
    for ( Int_t i = 0; i < 64; ++i )
    {
      if ( mask & ( 1ull << i ) ) return i+1;
    }
    return 0;
  }
  
  
  //______________________________________________________________________________
  Double_t Mu(Double_t L0B, Double_t Nb)
  {
    /// L0B = trigger rate before any veto
    /// Nb = number of crossing bunches
    
    Double_t p0 = 1-L0B/(Nb*11245.0); // proba to get *no* collision per bunch crossing
    
    return -TMath::Log(p0);
  }
  
  //______________________________________________________________________________
  TCanvas* NewCanvas(const char* name)
  {
    TCanvas* c = new TCanvas(name,name);
    c->SetLeftMargin(0.15);
    return c;
  }
  
  //______________________________________________________________________________
  void TimeAxis(TGraph* g)
  {
//    g->GetXaxis()->SetTimeDisplay(1);
//    g->GetXaxis()->SetTimeFormat("%d/%m %H:%M%F2010-12-31 24:00:00");
//    g->GetXaxis()->SetNdivisions(505);
    g->GetYaxis()->SetTitleOffset(1.5);
    
    g->GetXaxis()->SetTimeDisplay(1);
    //  g->GetXaxis()->SetTimeFormat("%d/%m %H:%M%F2010-12-31 24:00:00");
    g->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
    g->GetXaxis()->SetTimeOffset(0,"gmt");
    g->GetXaxis()->SetNdivisions(505);

  }
  
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::AliAnalysisTriggerScalers(Int_t runNumber, const char* ocdbPath) : fRunList(), fOCDBPath(ocdbPath), fShouldCorrectForPileUp(kTRUE), fCrossSectionUnit("UB")
{
  // ctor using a single run number
  SetRunList(runNumber);
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::AliAnalysisTriggerScalers(const std::vector<int>& runs, const char* ocdbPath) : fRunList(), fOCDBPath(ocdbPath), fShouldCorrectForPileUp(kTRUE), fCrossSectionUnit("UB")
{
  // ctor using a vector of run numbers
  SetRunList(runs);
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::AliAnalysisTriggerScalers(const std::set<int>& runs, const char* ocdbPath) : fRunList(), fOCDBPath(ocdbPath), fShouldCorrectForPileUp(kTRUE), fCrossSectionUnit("UB")
{
  // ctor using a set of run numbers
  SetRunList(runs);
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::AliAnalysisTriggerScalers(const char* runlist, const char* ocdbPath) : fRunList(), fOCDBPath(ocdbPath), fShouldCorrectForPileUp(kTRUE), fCrossSectionUnit("UB")
{
  // ctor from a run list (txt file)
  SetRunList(runlist);
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::~AliAnalysisTriggerScalers()
{
  // dtor
}

//______________________________________________________________________________
Bool_t AliAnalysisTriggerScalers::CheckRecord(const AliTriggerScalersRecord& record,
                                              UInt_t index,
                                              UInt_t refb,
                                              UInt_t refa,
                                              UInt_t timelapse) const
{
  /// Check if a trigger scaler record should be kept for our
  /// luminosity computations
  
  const AliTriggerScalers* scaler = record.GetTriggerScalersForClass(index);
  
  UInt_t a = scaler->GetLOCA() - refa;
  UInt_t b = scaler->GetLOCB() - refb;
  
  Bool_t ok(kTRUE);
  
  if ( b <= 2 || ( a <= 2 ) || timelapse <= 9 ) // empirical cuts
  {
    ok = kFALSE;
  }
  
  return ok;
}


//______________________________________________________________________________
void AliAnalysisTriggerScalers::DrawFill(Int_t run1, Int_t run2, double ymin, double ymax, const char* label)
{
  // Draw a yellow box to indicate a run range
  
  AliDebugClass(1,Form("RUN1 %09d RUN2 %09d YMIN %e YMAX %e %s",
                       run1,run2,ymin,ymax,label));
  TBox* b = new TBox(run1*1.0,ymin,run2*1.0,ymax);
  b->SetFillColor(5);
  b->Draw();
  TText* text = new TText((run1+run2)/2.0,ymin + (ymax-ymin)*0.85,label);
  text->SetTextSize(0.025);
  text->SetTextFont(42);
  text->SetTextAlign(23);
  text->SetTextAngle(45);
  text->Draw();
}

//_____________________________________________________________________________
void AliAnalysisTriggerScalers::DrawFills(Double_t ymin, Double_t ymax)
{
  /// Draw the fill ranges corresponding to the list of runs
  /// Note that this method will scan the OCDB to get the run -> fill number relationship,
  /// so it's better in this case to use a local copy of the OCDB if you have one. Otherwise
  /// it will be long.
  ///
  std::map<int, std::pair<int,int> > fills;
  
  GetFillBoundaries(fills);
  
  std::map<int, std::pair<int,int> >::const_iterator it;
  
  for ( it = fills.begin(); it != fills.end(); ++it )
  {
    const std::pair<int,int>& p = it->second;
    TString fillnumber;
    fillnumber.Form("%d",it->first);
    DrawFill(p.first,p.second,ymin,ymax,fillnumber.Data());
  }
}

//_____________________________________________________________________________
void AliAnalysisTriggerScalers::GetCTPObjects(Int_t runNumber,
                                              AliTriggerConfiguration*& tc,
                                              AliTriggerRunScalers*& trs,
                                              AliLHCData*& lhc) const
{
  /// For a given run, get the CTP objects we need to do our work
  ///
  /// FIXME: this is the method that should probably be virtualized so
  /// we can get those objects from either OCDB or run based container in AODs/ESDs
  
  tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
  trs = static_cast<AliTriggerRunScalers*>(GetOCDBObject("GRP/CTP/Scalers",runNumber));
  lhc = static_cast<AliLHCData*>(GetOCDBObject("GRP/GRP/LHCData",runNumber));
}

//_____________________________________________________________________________
void AliAnalysisTriggerScalers::GetFillBoundaries(std::map<int, std::pair<int,int> >& fills)
{
  /// Given a list of runs get the fills they are in
  
  fills.clear();
  
  for ( std::vector<int>::const_iterator it = fRunList.begin(); it != fRunList.end(); ++it )
  {
    int runnumber = *it;
    
    int fill = GetFillNumberFromRunNumber(runnumber);
    
    if (fills.count(fill))
    {
      std::pair<int,int>& p = fills[fill];
      p.first = TMath::Min(runnumber,p.first);
      p.second = TMath::Max(runnumber,p.second);
    }
    else
    {
      fills[fill] = make_pair<int,int>(runnumber,runnumber);
    }
  }
}


//______________________________________________________________________________
Int_t AliAnalysisTriggerScalers::GetFillNumberFromRunNumber(Int_t runNumber)
{
  /// Get the fill number of a run
  
  AliLHCData* lhcData = static_cast<AliLHCData*>(GetOCDBObject("GRP/GRP/LHCData",runNumber));
  if (lhcData)
  {
    Int_t fillNumber = lhcData->GetFillNumber();
    if ( fillNumber == 0 && runNumber == 189694)
    {
      // manual hack because GRP info incorrect for this run ?
      fillNumber = 3135;
    }

    return fillNumber;
  }
  return -1;
}

//______________________________________________________________________________
TObject* AliAnalysisTriggerScalers::GetOCDBObject(const char* path, Int_t runNumber) const
{
  // Helper method to get an object from the OCDB
  
  if ( !AliCDBManager::Instance()->IsDefaultStorageSet() )
  {
    AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  }
  
  AliCDBManager::Instance()->SetRun(runNumber);
  
  gErrorIgnoreLevel=kWarning; // to avoid all the TAlienFile::Open messages...
  
  AliCDBEntry* e = AliCDBManager::Instance()->Get(path);
  return e ? e->GetObject() : 0x0;
}


//______________________________________________________________________________
TString AliAnalysisTriggerScalers::GetLHCPeriodFromRunNumber(Int_t runNumber) const
{
  /// Get a list of LHC periods from a list of run numbers
  
  AliGRPObject* grp = static_cast<AliGRPObject*>(GetOCDBObject("GRP/GRP/Data",runNumber));

  if (!grp) return "";
  
  return grp->GetLHCPeriod();
}

//_____________________________________________________________________________
void
AliAnalysisTriggerScalers::GetLuminosityTriggerAndCrossSection(Int_t runNumber,
                                                               TString& lumiTriggerClassName,
                                                               Double_t& lumiTriggerCrossSection,
                                                               Double_t& lumiTriggerCrossSectionError) const
{
  /// For one given run, get the trigger class name to be used as the luminometer,
  /// and its corresponding cross-section
  ///
  /// FIXME: all is hard-coded for now... (use OADB for this ?)
  
  lumiTriggerClassName="";
  lumiTriggerCrossSection=0.0;
  lumiTriggerCrossSectionError=0.0; // FIXME
  
  TString period = GetLHCPeriodFromRunNumber(runNumber);
  
  if ( period.BeginsWith("LHC11") )
  {
    AliError("Not implemented yet");
  }
  else if ( period.BeginsWith("LHC12") )
  {
     if ( period == "LHC12h" || period == "LHC12i" )
     {
       // pp 8 TeV main-satellites
       lumiTriggerClassName = "C0TVX-S-NOPF-ALLNOTRD";
       lumiTriggerCrossSection = 28.0; // FIXME: not from a vdM yet       
     }
    else
    {
      AliError("Not implemented yet");
    }
  }
  else if ( period.BeginsWith("LHC13") )
  {
    if ( period == "LHC13g" )
    {
      // pp 2.76 TeV
      lumiTriggerClassName = "C0TVX-B-NOPF-ALLNOTRD";
      lumiTriggerCrossSection = 18.0; // FIXME: not from a vdM yet
    }
    else
    {
      // p-Pb or Pb-p 5.02 TeV
      lumiTriggerClassName = "C0TVX-B-NOPF-ALLNOTRD";
      lumiTriggerCrossSection = 0.755*2000; // FIXME: not from a vdM yet
    }
  }
  else
  {
    AliError("Not implemented yet");
  }
  
  float csMult(1.0);
  
  if ( fCrossSectionUnit=="PB" )
  {
    csMult=1E9;
  }
  else if (fCrossSectionUnit=="NB")
  {
    csMult=1E6;
  }
  else if ( fCrossSectionUnit=="UB" )
  {
    csMult=1E3;
  }
  else if ( fCrossSectionUnit=="MB" )
  {
    csMult=1.0;
  }

  lumiTriggerCrossSection *= csMult;
  lumiTriggerCrossSectionError *= csMult;
}

//_____________________________________________________________________________
void AliAnalysisTriggerScalers::GetLHCPeriodBoundaries(std::map<std::string, std::pair<int,int> >& periods)
{
  /// Given a list of runs get the fills they are in
  
  periods.clear();
  
  for ( std::vector<int>::const_iterator it = fRunList.begin(); it != fRunList.end(); ++it )
  {
    int runnumber = *it;
    
    std::string period = GetLHCPeriodFromRunNumber(runnumber).Data();
    
    if (periods.count(period))
    {
      std::pair<int,int>& p = periods[period];
      p.first = TMath::Min(runnumber,p.first);
      p.second = TMath::Max(runnumber,p.second);
    }
    else
    {
      periods[period] = make_pair<int,int>(runnumber,runnumber);
    }
  }
}


//______________________________________________________________________________
Int_t AliAnalysisTriggerScalers::GetTriggerInput(Int_t runNumber, const char* inputname)
{
  /// Get one trigger input for a given run
  AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
  if (!tc) return -1;
  
  const TObjArray& inputs = tc->GetInputs();
  
  AliTriggerInput* ti = static_cast<AliTriggerInput*>(inputs.FindObject(inputname));
  
  if (!ti) return -1;
  
  return ti->GetSignature();
}

//______________________________________________________________________________
Float_t
AliAnalysisTriggerScalers::GetPauseAndConfigCorrection(Int_t runNumber, const char* triggerClassName)
{
  /// Tries to estimate the duration of the Pause and Configure phase(s) in a given run
  /// Probably not really accurate for the moment though.
  
  ULong64_t total(0);
  ULong64_t nopac(0);
  
  AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));

  AliTriggerRunScalers* trs = static_cast<AliTriggerRunScalers*>(GetOCDBObject("GRP/CTP/Scalers",runNumber));
  const TObjArray* scalers = trs->GetScalersRecords();
  
  TIter next(scalers);
  AliTriggerScalersRecord* record;
  UInt_t reft(0);
  UInt_t refl0b(0);
  UInt_t refl0a(0);
  Bool_t first(kTRUE);

  Int_t classIndex = tc->GetClassIndexFromName(triggerClassName);
  
  while ( ( record = static_cast<AliTriggerScalersRecord*>(next()) ) )
  {
    const AliTriggerScalers* scaler = record->GetTriggerScalersForClass(classIndex);
        
    const AliTimeStamp* ats = record->GetTimeStamp();
    
    UInt_t seconds = ats->GetSeconds();// - TTimeStamp::GetZoneOffset();
    
    TTimeStamp ts(seconds,ats->GetMicroSecs());
    
    UInt_t l0b = scaler->GetLOCB() - refl0b;
    UInt_t l0a = scaler->GetLOCA() - refl0a;
    UInt_t timelapse = seconds - reft;

    if ( l0b <=2 || timelapse < 9 ) continue;

    reft = seconds;
    refl0b = scaler->GetLOCB();
    refl0a = scaler->GetLOCA();

    if ( first )
    {
      first = kFALSE;
      continue;
    }
    
    total += l0b;
    
    if ( l0a > 0 )
    {
      nopac += l0b;
    }
  }

  return total > 0 ? nopac*1.0/total : 0.0;
}

//______________________________________________________________________________
AliAnalysisTriggerScalerItem* 
AliAnalysisTriggerScalers::GetTriggerScaler(Int_t runNumber, const char* level, const char* triggerClassName)
{
  // Get the scaler for a given trigger class and a given trigger level.
  // Returned object must be deleted by the user.
  
  AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
  AliTriggerRunScalers* trs = static_cast<AliTriggerRunScalers*>(GetOCDBObject("GRP/CTP/Scalers",runNumber));
  AliGRPObject* grp = static_cast<AliGRPObject*>(GetOCDBObject("GRP/GRP/Data",runNumber));

  if (!tc || !trs || !grp) return 0x0;
  
  TString diCurrent(Form("L3:%5.0f;DIP:%5.0f [L3:%d;DIP:%d]",
                         grp->GetL3Current((AliGRPObject::Stats)0),
                         grp->GetDipoleCurrent((AliGRPObject::Stats)0),
                         grp->GetL3Polarity(),
                         grp->GetDipolePolarity()));
  
  const TObjArray& trClasses = tc->GetClasses();
  
  const TObjArray* scalers = trs->GetScalersRecords();  
  const AliTriggerScalersRecord* begin = (AliTriggerScalersRecord*)(scalers->First());
  const AliTriggerScalersRecord* end = (AliTriggerScalersRecord*)(scalers->Last());

  time_t duration = TMath::Nint((end->GetTimeStamp()->GetBunchCross() - begin->GetTimeStamp()->GetBunchCross())*AliTimeStamp::fNanosecPerBC*1E-9);
  
  AliTriggerClass* triggerClass = static_cast<AliTriggerClass*>(trClasses.FindObject(triggerClassName));
  if (!triggerClass)
  {
    return 0x0;
  }
  UChar_t index = GetIndex(triggerClass->GetMask());
  
  ULong64_t value(0);
  
  const AliTriggerScalers* sbegin = begin->GetTriggerScalersForClass(index);
  const AliTriggerScalers* send = end->GetTriggerScalersForClass(index);
  
  TString swhat(level);
  swhat.ToUpper();
  
  if ( swhat.BeginsWith("L0A") )
  {
    value = send->GetLOCA() - sbegin->GetLOCA();
  }
  else if ( swhat.BeginsWith("L0B") )
  {
    value = send->GetLOCB() - sbegin->GetLOCB();
  }
  else if ( swhat.BeginsWith("L1A") )
  {
    value = send->GetL1CA() - sbegin->GetL1CA();
  }
  else if ( swhat.BeginsWith("L1B") )
  {
    value = send->GetL1CB() - sbegin->GetL1CB();
  }
  else if ( swhat.BeginsWith("L2A") )
  {
    value = send->GetL2CA() - sbegin->GetL2CA();
  }
  else if ( swhat.BeginsWith("L2B") )
  {
    value = send->GetL2CB() - sbegin->GetL2CB();
  }
  Int_t ds(1);
  
  // FIXME: get the downscaling factor here for all cases (only BC1 supported here so far)
  if ( TString(triggerClassName).Contains("_B1") )
  {
    ds = GetTriggerInput(runNumber,"BC1");
    if (!ds) ds=1;
  }
  
  swhat.ReplaceAll("RATE","");
  
  return new AliAnalysisTriggerScalerItem(runNumber,swhat.Data(),diCurrent.Data(),triggerClass->GetName(),value,triggerClass->GetBCMask(),ds,duration);
}

//______________________________________________________________________________
TGraph*
AliAnalysisTriggerScalers::IntegratedLuminosityGraph(const char* triggerClassName, const char* triggerClassNameForPACEstimate)
{
  // Get the integrated luminosity graph for one trigger
  
  if ( fRunList.empty() )
  {
    AliError("Cannot work without a run list");
    return 0x0;
  }
  
  Double_t offset(0.0);
  
  std::vector<double> vx;
  std::vector<double> vy;
  
  Double_t x,y;
  
  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    TGraph* g = IntegratedLuminosityGraph(fRunList[i],triggerClassName,triggerClassNameForPACEstimate);
    
    if (!g) continue;
    
    for ( Int_t j = 0 ; j < g->GetN(); ++j )
    {
      g->GetPoint(j,x,y);
      vx.push_back(x);
      vy.push_back(y+offset);
    }
    
    offset += y;
    
    delete g;
  }
  
  TGraph* g = new TGraph(vx.size(),&vx[0],&vy[0]);
  
  TimeAxis(g);
  
  return g;
}

//______________________________________________________________________________
TGraph*
AliAnalysisTriggerScalers::IntegratedLuminosityGraph(Int_t runNumber, const char* triggerClassName, const char* triggerClassNameForPACEstimate)
{
  // Get the integrated luminosity graph for one trigger for the given run
  //
  // the triggerClassNameForPACEstimate is only used if triggerClassName is,
  // for the given run, the same as the luminometer trigger class name,
  // and should be the trigger class with the higher and non-zero L0A rate (except during PACs)
  //
  
  AliTriggerConfiguration* tc(0x0);
  AliTriggerRunScalers* trs(0x0);
  AliLHCData* lhc(0x0);
  
  GetCTPObjects(runNumber,tc,trs,lhc);
  
  if (!tc || !trs || !lhc)
  {
    return 0x0;
  }
  
  const TObjArray* scalerRecords= trs->GetScalersRecords();
  const TObjArray& triggerClasses = tc->GetClasses();
  
  TString lumiTriggerClassName;
  Double_t lumiTriggerCrossSection(0.0);
  Double_t lumiTriggerCrossSectionError(0.0);
  
  GetLuminosityTriggerAndCrossSection(runNumber,lumiTriggerClassName,lumiTriggerCrossSection,lumiTriggerCrossSectionError);

  AliTriggerClass* lumiTriggerClass = static_cast<AliTriggerClass*>(triggerClasses.FindObject(lumiTriggerClassName));

  if (!lumiTriggerClass)
  {
    AliError(Form("Could not find lumiTriggerClass=%s",lumiTriggerClassName.Data()));
    return 0x0;
  }
  
  if (lumiTriggerCrossSection==0.0)
  {
    AliError(Form("Could not get cross-section for lumiTriggerClass=%s",lumiTriggerClassName.Data()));
    return 0x0;
  }
  
  AliTriggerClass* triggerClass = static_cast<AliTriggerClass*>(triggerClasses.FindObject(triggerClassName));

  if (!triggerClass)
  {
    AliError(Form("Could not find triggerClass=%s",triggerClassName));
    return 0x0;
  }

  Int_t nbcx = NumberOfInteractingBunches(*lhc);
  
  if (nbcx <= 0 && ShouldCorrectForPileUp())
  {
    AliError("Could not get the number of bunches, so won't be able to correct for pile-up");
  }
  
  TIter nextScaler(scalerRecords);
  AliTriggerScalersRecord* record;
  UInt_t reft(0);
  
  UInt_t refl0b[] = { 0, 0 };
  UInt_t refl0a[] = { 0, 0 };

  UInt_t l0b[] = { 0, 0 };
  UInt_t l0a[] = { 0, 0 };
  
  Bool_t first(kTRUE);
  
  UChar_t classIndices[2];
  
  // class 0 = class for luminomiter
  // class 1 = class for livetime estimate
  //           or for PAC estimate if triggerClassNameForPACEstimate
  //              is given and triggerClass is the same as the lumi class
  
  classIndices[0] = GetIndex(lumiTriggerClass->GetMask());
  classIndices[1] = GetIndex(triggerClass->GetMask());

  Bool_t sameClass = ( classIndices[0] == classIndices[1] );
  Bool_t pacRemoval(kFALSE);
  
  if ( sameClass && strlen(triggerClassNameForPACEstimate) > 0 )
  {
    AliTriggerClass* triggerClassForPACEstimate = static_cast<AliTriggerClass*>(triggerClasses.FindObject(triggerClassNameForPACEstimate));
    if (!triggerClassForPACEstimate)
    {
      AliError(Form("Could not find triggerClassForPACEstimate=%s. Will not correct for PAC durations",triggerClassNameForPACEstimate));
      return 0x0;
    }
    classIndices[1] = GetIndex(triggerClassForPACEstimate->GetMask());
    sameClass = ( classIndices[0] == classIndices[1] );
    if (!sameClass)
    {
      pacRemoval=kTRUE;
    }
  }
  
  ULong64_t counter(0);
  
  std::vector<Double_t> vseconds;
  std::vector<Double_t> vlivetime;
  std::vector<Double_t> vlumi;
  
  while ( ( record = static_cast<AliTriggerScalersRecord*>(nextScaler()) ) )
  {
    const AliTimeStamp* ats = record->GetTimeStamp();

    UInt_t seconds = ats->GetSeconds();// - TTimeStamp::GetZoneOffset();
    
    TTimeStamp ts(seconds,ats->GetMicroSecs());

    UInt_t timelapse = seconds - reft;

    Bool_t ok = sameClass || CheckRecord(*record,classIndices[1],refl0b[1],refl0a[1],timelapse);
    
    if (ok) reft = seconds;
    
    for ( Int_t i = 0; i < 2; ++i )
    {
      const AliTriggerScalers* scaler = record->GetTriggerScalersForClass(classIndices[i]);
      
      if (ok)
      {
        l0b[i] = scaler->GetLOCB() - refl0b[i];
      }
      else
      {
        l0b[i] = 0;
      }
      l0a[i] = scaler->GetLOCA() - refl0a[i];

      refl0b[i] = scaler->GetLOCB();
      refl0a[i] = scaler->GetLOCA();
    }
    
    if ( first )
    {
      first = kFALSE;
      continue;
    }
    
    Double_t factor(1.0);
    
    if ( ShouldCorrectForPileUp() && nbcx > 0 && l0b[0] > 0 )
    {
      Double_t mu = Mu(l0b[0]/timelapse,nbcx);
      factor = mu/(1-TMath::Exp(-mu));
    }
    
    vseconds.push_back(seconds);
    
    Double_t lt(1.0);
    
    if ( ok && !sameClass && !pacRemoval )
    {
      lt = (l0a[1]*1.0)/l0b[1];
    }

    counter += TMath::Nint(factor*l0b[0]*lt);

    vlumi.push_back( counter / lumiTriggerCrossSection );
    
  }
  
  if ( vseconds.empty() ) return 0x0;
  
  TGraph* g = new TGraph(vseconds.size(),&vseconds[0],&vlumi[0]);
  
  TimeAxis(g);

  return g;
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::IntegratedLuminosity(const char* triggerList,
                                                     const char* lumiTrigger,
                                                     Double_t lumiCrossSection,
                                                     const char* csvOutputFile,
                                                     const char sep,
                                                     const char* csUnit)
{
  // Compute the luminosity for a set of triggers
  
  // for T0 based lumi (end of pp 2012), use lumiTrigger = C0TVX-S-NOPF-ALLNOTRD and lumiCrossSection = 28 mb (and csUnit="nb")
  //
  // for T0 based lumi (pPb 2013), use lumiTrigger = C0TVX-B-NOPF-ALLNOTRD
  // and lumiCrossSection = 0.755*2000 mb (and csUnit="mub")
  
  // for T0 based lumi (pp 2.76 TeV 2013), use lumiTrigger = C0TVX-B-NOPF-ALLNOTRD
  // and lumiCrossSection = 18 mb (and csUnit="nb")
  
  double intLumi(0);
  
  std::map<std::string,float> lumiPerTrigger;
  std::map<int, std::map<std::string,float> > lumiPerFillPerTrigger;

  std::map<std::string,time_t> durationPerTrigger;
  
  TString sTriggerList(triggerList);
  
  if ( sTriggerList.Length()==0 )
  {
    if ( lumiCrossSection < 30 && lumiCrossSection  > 20 )
    {
      // assuming it's pp 2012
      sTriggerList = "CMUL8-S-NOPF-ALLNOTRD,CMUL7-S-NOPF-ALLNOTRD,CMUL8-S-NOPF-MUON,CMUL7-S-NOPF-MUON";
    
    // for C0MUL-SC-NOPF-MUON must use C0TVX-SC as lumiTrigger (with same cross-section as C0TVX=28mb)
    }
    else if ( lumiCrossSection < 30 )
    {
      // assuming it's pp 2013
      sTriggerList = "CMSL7-B-NOPF-MUON,CMSH7-B-NOPF-MUON,CMUL7-B-NOPF-MUON,CMSL7-B-NOPF-ALLNOTRD,CMSH7-B-NOPF-ALLNOTRD,CMUL7-B-NOPF-ALLNOTRD";
    }
    else
    {
      sTriggerList = "CMSL7-B-NOPF-MUON,CMSH7-B-NOPF-MUON,CMUL7-B-NOPF-MUON,CMSL7-B-NOPF-ALLNOTRD,CMSH7-B-NOPF-ALLNOTRD,CMUL7-B-NOPF-ALLNOTRD";
    }
  }
  
  TObjArray* triggerArray = sTriggerList.Tokenize(",");
  triggerArray->SetOwner(kTRUE);

  std::ofstream* out(0x0);

  if ( TString(csvOutputFile).Length() > 0 )
  {
    out = new std::ofstream(gSystem->ExpandPathName(csvOutputFile));
    if (!out || out->bad())
    {
      delete out;
      out = 0x0;
    }
  }
  
  TIter nextTrigger(triggerArray);
  TObjString* trigger;
  
  if (out)
  {
    (*out) << "Fill number" << sep;
    
    nextTrigger.Reset();
    
    std::map<std::string,float>::const_iterator it;

    while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
    {
      lumiPerTrigger[trigger->String().Data()] = 0;
    }
    
    for ( it = lumiPerTrigger.begin(); it != lumiPerTrigger.end(); ++it )
    {
      (*out) << "lumi from " << it->first.c_str() << sep;
    }

    (*out) << "comments" << sep;
    
    for ( it = lumiPerTrigger.begin(); it != lumiPerTrigger.end(); ++it )
    {
      (*out) << "recorded " <<  it->first.c_str() << " (integrated)" << sep;
    }

    (*out) << Form("LHC delivered (%s^-1) per fill ",csUnit)
      << sep << Form("LHC delivered (%s^-1 integrated)",csUnit) << sep;
    (*out) << "lumi tot muon" << sep << "eff (%)" << sep;
    (*out) << std::endl;
    
    nextTrigger.Reset();
  }
  
  Int_t currentFillNumber(-1);
  Int_t fillNumber(0);
  
  float csMult(1.0);
  TString scsUnit(csUnit);
  scsUnit.ToUpper();
  
  if ( scsUnit=="PB" )
  {
    csMult=1E9;
  }
  else if (scsUnit=="NB")
  {
    csMult=1E6;
  }
  else if ( scsUnit=="UB" )
  {
    csMult=1E3;
  }
  else if ( scsUnit=="MB" )
  {
    csMult=1.0;
  }
  else
  {
    AliError(Form("Don't know how to deal with cross-section unit=%s",csUnit));
    return;
  }
  

  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];
    Bool_t atLeastOneTriggerFound(kFALSE);
    
    // find out which trigger classes to use
    
    AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
    const TObjArray& trClasses = tc->GetClasses();
        
    if (out)
    {
      fillNumber = GetFillNumberFromRunNumber(runNumber);
      if ( fillNumber == 0 )
      {
        AliError(Form("Got fillNumber = 0 for run %09d !",runNumber));
      }
      
      if ( fillNumber != currentFillNumber )
      {
        std::map<std::string,float>::const_iterator it;
        
        for ( it = lumiPerTrigger.begin(); it != lumiPerTrigger.end(); ++it )
        {
          lumiPerFillPerTrigger[fillNumber][it->first.c_str()] = 0;
        }
        currentFillNumber = fillNumber;
      }
    }

    AliTriggerRunScalers* trs = static_cast<AliTriggerRunScalers*>(GetOCDBObject("GRP/CTP/Scalers",runNumber));
    const TObjArray* scalers = trs->GetScalersRecords();
    const AliTriggerScalersRecord* begin = (AliTriggerScalersRecord*)(scalers->First());
    const AliTriggerScalersRecord* end = (AliTriggerScalersRecord*)(scalers->Last());
    
    time_t duration = TMath::Nint((end->GetTimeStamp()->GetBunchCross() - begin->GetTimeStamp()->GetBunchCross())*AliTimeStamp::fNanosecPerBC*1E-9);

    nextTrigger.Reset();

    TString lumiTriggerClassName(lumiTrigger);
    float lumiSigma = lumiCrossSection*csMult; //in csUnit
    AliAnalysisTriggerScalerItem* lumiB = GetTriggerScaler(runNumber,"L0B",lumiTriggerClassName.Data());

    if (!lumiB)
    {
      AliError(Form("Did not find lumiTrigger %s for run %09d",lumiTriggerClassName.Data(),runNumber));
      continue;
    }
        
    Float_t pacCorrection(1.0);
    
    while ( ( trigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      TString muTriggerClassName(trigger->String());
      Int_t n(0);
      
      if ( !trClasses.FindObject(muTriggerClassName.Data() ) )
      {
        continue;
      }
      
      if ( muTriggerClassName.Contains("CMUL8") ) ++n;
      if ( muTriggerClassName.Contains("CMUL7") ) ++n;
      
      if ( n>1 )
      {
        AliError(Form("More than 1 relevant trigger class found for run %09d ! Check that !",runNumber));
        trClasses.Print();
        continue;
      }
      
      AliAnalysisTriggerScalerItem* muonA = GetTriggerScaler(runNumber,"L0A",muTriggerClassName.Data());
      AliAnalysisTriggerScalerItem* muonB = GetTriggerScaler(runNumber,"L0B",muTriggerClassName.Data());
      
      if (!muonA || !muonB) continue;
      
      atLeastOneTriggerFound = kTRUE;
      
      Float_t ratio(0);
      
      if ( muonB->ValueCorrectedForDownscale() > 0 )
      {
        ratio = muonA->ValueCorrectedForDownscale()/muonB->ValueCorrectedForDownscale();
      }
      
      ratio *= lumiB->ValueCorrectedForDownscale()/lumiSigma;
      
      if ( muTriggerClassName.BeginsWith("CMUL"))
      {
        intLumi += ratio;
      }
      
      if ( muTriggerClassName.Contains("CMSH") )
      {
        pacCorrection = GetPauseAndConfigCorrection(runNumber,muTriggerClassName.Data());
      }
      
      lumiPerTrigger[muTriggerClassName.Data()] += ratio;
      durationPerTrigger[muTriggerClassName.Data()] += duration;
      lumiPerFillPerTrigger[currentFillNumber][muTriggerClassName.Data()] += ratio;
      
    }
    
    lumiPerTrigger[lumiTriggerClassName.Data()] += lumiB->ValueCorrectedForDownscale()/lumiSigma;
    durationPerTrigger[lumiTriggerClassName.Data()] += duration;
    lumiPerFillPerTrigger[currentFillNumber][lumiTriggerClassName.Data()] += lumiB->ValueCorrectedForDownscale()/lumiSigma;

    TString lumiPACCorrected(Form("%s(-PAC)",lumiTriggerClassName.Data()));
    
    lumiPerTrigger[lumiPACCorrected.Data()] += pacCorrection*lumiB->ValueCorrectedForDownscale()/lumiSigma;
    durationPerTrigger[lumiPACCorrected.Data()] += duration;
    lumiPerFillPerTrigger[currentFillNumber][lumiPACCorrected.Data()] += pacCorrection*lumiB->ValueCorrectedForDownscale()/lumiSigma;

    if (!atLeastOneTriggerFound && sTriggerList.Contains("CMUL") )
    {
      AliError(Form("Found no relevant trigger for run %09d",runNumber));
    }
  }
  
  AliInfo(Form("Integrated lumi %7.4f %s^-1",intLumi,csUnit));
  
  std::map<std::string,float>::const_iterator it;
  
  for ( it = lumiPerTrigger.begin(); it != lumiPerTrigger.end(); ++it )
  {
    AliInfo(Form("Trigger %30s Lumi %10.4f %s^-1 duration %10ld s",it->first.c_str(),it->second,csUnit,durationPerTrigger[it->first]));
  }
  
  if (out)
  {
    std::map<int, std::map<std::string, float> >::const_iterator fit;
    
    lumiPerTrigger.clear();
    
    for ( fit = lumiPerFillPerTrigger.begin(); fit != lumiPerFillPerTrigger.end(); ++fit )
    {
      int fill = fit->first;
      std::map<std::string,float>::const_iterator tit;
      (*out) << fill << sep;
      
      for ( tit = fit->second.begin(); tit != fit->second.end(); ++tit )
      {
        (*out) << Form("%e",tit->second) << sep;
      }
      
      (*out) << sep; // comment (empty)
      
      for ( tit = fit->second.begin(); tit != fit->second.end(); ++tit )
      {
        lumiPerTrigger[tit->first] += tit->second;        
        
        (*out) <<  Form("%e",lumiPerTrigger[tit->first]) << sep;
      }
      (*out) << sep << "0" << sep << "0" << sep << "0" << sep << "0" << sep << std::endl; // LHC per fill, LHC integrated, lumi tot muon , efficiency
    }
  }
  //
  delete triggerArray;
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::MakeGraph(const std::vector<int>& vx,
                                             const std::vector<int>& vex,
                                             const std::vector<double>& vy,
                                             const std::vector<double>& vey)
{
  /// Build a graph from a set of stl vectors
  
  if ( ! ( vx.size() == vex.size() && vx.size() == vy.size() && vx.size() == vey.size() ) )
  {
    std::cerr << "incompatible sizes" << std::endl;
    return 0x0;
  }
  
  double* x = new double[vx.size()];
  double* ex = new double[vx.size()];
  double* y = new double[vx.size()];
  double* ey = new double[vx.size()];
  
  for ( size_t i = 0; i < vx.size(); ++i )
  {
    x[i] = vx[i];
    ex[i] = vex[i];
    y[i] = vy[i];
    ey[i] = vey[i];
  }
  
  TGraph* g = new TGraphErrors(vx.size(),x,y,ex,ey);
  
  TimeAxis(g);
  
  delete[] x;
  delete[] y;
  delete[] ex;
  delete[] ey;
  
  return g;
}

//______________________________________________________________________________
Int_t AliAnalysisTriggerScalers::NumberOfInteractingBunches(const AliLHCData& lhc) const
{
  /// Extract the number of colliding bunches from the LHC data
  
  Int_t numberOfInteractingBunches(0);
  Int_t numberOfInteractingBunchesMeasured(0);
  Int_t nIBM2(0);

  int beam1(0);
  int beam2(1);

  AliLHCDipValI* val = lhc.GetBunchConfigDeclared(beam1,0);

  for ( Int_t i = 0; i < val->GetSizeTotal(); ++i )
  {
    if ( val->GetValue(i) < 0 ) ++numberOfInteractingBunches;
  }

  AliLHCDipValI* valm = lhc.GetBunchConfigMeasured(beam1,0);

  for ( Int_t i = 0; i < valm->GetSizeTotal(); ++i )
  {
    if ( valm->GetValue(i) < 0 ) ++numberOfInteractingBunchesMeasured;
  }

  valm = lhc.GetBunchConfigMeasured(beam2,0);

  for ( Int_t i = 0; i < valm->GetSizeTotal(); ++i )
  {
    if ( valm->GetValue(i) <= 0 ) ++nIBM2;
  }

  if ( numberOfInteractingBunches != numberOfInteractingBunchesMeasured ||
       numberOfInteractingBunches != nIBM2 )
  {
    AliWarning(Form("Got some different number of interacting bunches here ! NumberOfInteractingBunches=%3d NumberOfInteractingBunchesMeasured=%3d NIBM2=%3d",
                  numberOfInteractingBunches,numberOfInteractingBunchesMeasured,nIBM2));
  }
  
  return numberOfInteractingBunches;
}


//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTrigger(const char* triggerClassName,
                                               const char* what)
{
  // Plot one of the scalers (L0A,L0B,L0AOVERB,etc...) of a given triggerClass
  // Get one value per run.
  
  std::vector<Float_t> x;
  std::vector<Float_t> y;

  double integral(0);
  double mean(0);
  
  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];
    
    AliAnalysisTriggerScalerItem* s = GetTriggerScaler(runNumber,what,triggerClassName);
    
    if (!s) continue;
    
    x.push_back(runNumber);

    Double_t value = s->ValueCorrectedForDownscale();
    
    if ( TString(what).Contains("RATE") )
    {
      value = s->Rate();
    }
  
    integral += value;
    mean += value;
    
    y.push_back(value);
  }
  
  if ( fRunList.size() ) {
    mean /= fRunList.size();
  }
  
  AliInfo(Form("Integral %e mean %e",integral,mean));
  
  return new TGraph(x.size(),&x[0],&y[0]);
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTriggerEvolution(const char* triggerClassName,
                                                        const char* what,
                                                        bool draw,
                                                        double* mean,
                                                        bool removeZeros)
{
  /// Make a graph of "what" for a given trigger class.
  ///
  /// What can be :
  ///
  /// - L0B (level 0 before any veto applied)
  /// - L0A (level 0 after all vetoes)
  /// - L0AOVERB (L0A/L0B)
  /// - mu ( = -TMath::Log( 1 - P(0) ) where P(0) is the proba to have zero collisions in a bunch crossing)
  /// - pileupfactor = mu/(1-exp(-mu)) : the factor to apply to correct the cint1b count rate
  /// - vsnb = L0B/(NumberOfInteractingBunches*11245)
  /// - NBCX = NumberOfInteractingBunches
  
  TString swhat(what);
  swhat.ToUpper();
  
  if ( swhat.Contains(";"))
  {
    swhat.ReplaceAll(";",",");
    AliWarningClass("; is not a valid separator, replaced it with ,");
  }
  
  int color(1);
  int marker(20);
  
  TObjArray* a = swhat.Tokenize(",");
  if (a->GetEntries()>1)
  {
    TObjArray graphs;
    TIter next(a);
    TObjString* str(0x0);
    Double_t ymin(TMath::Limits<Double_t>::Max());
    Double_t ymax(0);
    Double_t xmin(TMath::Limits<Double_t>::Max());
    Double_t xmax(0);
    TGraph* g(0x0);
    
    while ( ( str = static_cast<TObjString*>(next())) )
    {
      g = PlotTriggerEvolution(triggerClassName,str->String().Data(),false,0x0,removeZeros);
      graphs.Add(g);
      for ( Int_t i = 0; i < g->GetN(); ++i )
      {
        ymin = TMath::Min(ymin,g->GetY()[i]);
        ymax = TMath::Max(ymax,g->GetY()[i]);
      }
      xmin = TMath::Min(xmin,g->GetX()[0]);
      xmax = TMath::Max(xmax,g->GetX()[g->GetN()-1]);
    }
    
    gStyle->SetOptTitle(0);
    
    AliInfoClass(Form("x %e ; %e y %e ; %e",xmin,xmax,ymin,ymax));
    TH2* h = new TH2F("h",triggerClassName,100,xmin,xmax,100,ymin,ymax);
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTimeDisplay(1);
    h->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
    h->GetXaxis()->SetTimeOffset(0,"gmt");
    h->GetXaxis()->SetNdivisions(505);
    h->Draw();
    
    TIter nextGraph(&graphs);
    
    while ( ( g = static_cast<TGraph*>(nextGraph())) )
    {
      AliInfoClass(g->GetTitle());
      g->Draw("lp");
      g->SetLineColor(color);
      g->SetMarkerColor(color);
      g->SetMarkerStyle(marker);
      ++color;
      ++marker;
    }
    delete a;
    return 0x0;
  }
  delete a;
  
  std::vector<int> vx;
  std::vector<int> vex;
  std::vector<double> vy;
  std::vector<double> vey;
  
  
  if (mean) *mean=0;
  double nvalues(0);
  
  for ( std::vector<int>::size_type iRun = 0; iRun < fRunList.size(); ++iRun )
  {
    Int_t runNumber = fRunList[iRun];
    
    AliTriggerConfiguration* tc(0x0);
    AliTriggerRunScalers* trs(0x0);
    AliLHCData* lhc(0x0);
    
    GetCTPObjects(runNumber,tc,trs,lhc);
    
    Int_t numberOfInteractingBunches = NumberOfInteractingBunches(*lhc);
    
    const TObjArray* scalers = trs->GetScalersRecords();
    
    const TObjArray& trClasses = tc->GetClasses();
    TIter next(&trClasses);
    AliTriggerClass* triggerClass;
    
    while ( ( triggerClass = static_cast<AliTriggerClass*>(next()) ) )
    {
      UChar_t index = GetIndex(triggerClass->GetMask());
      
      if ( !TString(triggerClass->GetName()).Contains(triggerClassName) ) continue;
      
      TIter nextScaler(scalers);
      AliTriggerScalersRecord* record;
      UInt_t reft(0);
      UInt_t refl0b(0);
      UInt_t refl1b(0);
      UInt_t refl2b(0);
      UInt_t refl0a(0);
      UInt_t refl1a(0);
      UInt_t refl2a(0);
      Bool_t first(kTRUE);
      
      while ( ( record = static_cast<AliTriggerScalersRecord*>(nextScaler()) ) )
      {
        const AliTriggerScalers* scaler = record->GetTriggerScalersForClass(index);
        
        
        const AliTimeStamp* ats = record->GetTimeStamp();
        
        UInt_t seconds = ats->GetSeconds();// - TTimeStamp::GetZoneOffset();
        
        TTimeStamp ts(seconds,ats->GetMicroSecs());
        
        UInt_t l0b = scaler->GetLOCB() - refl0b;
        UInt_t l0a = scaler->GetLOCA() - refl0a;
        UInt_t l1b = scaler->GetL1CB() - refl1b;
        UInt_t l1a = scaler->GetL1CA() - refl1a;
        UInt_t l2b = scaler->GetL2CB() - refl2b;
        UInt_t l2a = scaler->GetL2CA() - refl2a;
        UInt_t timelapse = seconds - reft;
        
        if ( removeZeros && ( l0b <= 2 || ( l0a <= 2 && l0a != 0 ) || timelapse <= 9 ) )
        {
          AliInfo(Form("Skipping point for %s l0b %d l0a %d timelapse %d ts %s",
                       triggerClassName,l0b,l0a,timelapse,ts.AsString()));
          continue;
        }
        
        reft = seconds;
        refl0b = scaler->GetLOCB();
        refl0a = scaler->GetLOCA();
        refl1b = scaler->GetL1CB();
        refl1a = scaler->GetL1CA();
        refl2b = scaler->GetL2CB();
        refl2a = scaler->GetL2CA();
        
        if ( first )
        {
          first = kFALSE;
          continue;
        }
        
        double value(1.0);
        double error(0.0);
        
        if (swhat.Contains("L0AOVERB") )
        {
          value = l0a*1.0/l0b;
          if ( l0a > 0 )
          {
            error = value*TMath::Sqrt(1.0/l0b+1.0/l0a);
          }
          else
          {
            error = 0.0;
          }
        }
        else if ( swhat.Contains("L0B") )
        {
          value = l0b;
        }
        else if (swhat.Contains("L0A") )
        {
          value = l0a;
        }
        else if ( swhat.Contains("L1B") )
        {
          value = l1b;
        }
        else if (swhat.Contains("L1A") )
        {
          value = l1a;
        }
        else if ( swhat.Contains("L2B") )
        {
          value = l2b;
        }
        else if (swhat.Contains("L2A") )
        {
          value = l2a;
        }
        else if ( swhat.Contains("MU") )
        {
          value = Mu(l0b/timelapse,numberOfInteractingBunches);
          error = 0.0; // FIXME
        }
        else if ( swhat.Contains("PILEUPFACTOR") )
        {
          Double_t mu = Mu(l0b/timelapse,numberOfInteractingBunches);
          value = mu/(1-TMath::Exp(-mu));
          error = -1.0; // FIXME
        }
        else if ( swhat.Contains("VSNB") )
        {
          value = l0b/(11245.0*numberOfInteractingBunches);
          error = -1.0; // FIXME
        }
        else if ( swhat.Contains("NBCX"))
        {
          value = numberOfInteractingBunches;
          error = 1.0; // FIXME          
        }
        else
        {
          value = timelapse;
          AliInfo(Form("Unknown what %s",what));
        }
        
        if ( ! swhat.Contains("OVER") && ! swhat.Contains("RATIO") &&
            ! swhat.Contains("MU") && ! swhat.Contains("PILEUPFACTOR") &&
            ! swhat.Contains("RAW") & ! swhat.Contains("NBCX") )
        {
          value /= timelapse;
        }
        
        if ( !TMath::Finite(value) ) continue;
        
        if ( !swhat.Contains("L0AOVERB") && error > 0 )
        {
          error = value >0 ? 1.0/TMath::Sqrt(value) : 0.0;
        }
        
        if (removeZeros && TMath::Abs(value) < 1E-6 )
        {
          continue;
        }
        if (mean)
        {
          (*mean) += value;
          nvalues += 1.0;
        }
        
        vx.push_back(seconds);
        vex.push_back(1);
        
        vy.push_back(value);
        
        vey.push_back(error);
        
      }
      
    }
    
  }
  
  if (mean && nvalues)
  {
    (*mean) /= nvalues;
  }
  
  if ( vx.empty() ) return 0;
  
  TGraph* g = MakeGraph(vx,vex,vy,vey);
  TString title(Form("TriggerEvolution%s-%s",triggerClassName,swhat.Data()));
  
  g->SetName(title.Data());
  g->SetTitle(title.Data());
  
  g->GetYaxis()->SetTitle(title.Data());
  
  if (draw)
  {
    TCanvas* c = NewCanvas(g->GetName());
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(marker);
    g->Draw("ALP");
    c->SaveAs(Form("%s.png",title.Data()));
  }
  
  
  return g;
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTriggerRatio(const char* triggerClassName1,
                                                    const char* what1,
                                                    const char* triggerClassName2,
                                                    const char* what2)
{
  // Plot the ratio of two scalers.
  // Get one value per run.

  std::vector<Float_t> x;
  std::vector<Float_t> y;
  
  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];
    
    AliAnalysisTriggerScalerItem* s1 = GetTriggerScaler(runNumber,what1,triggerClassName1);
    AliAnalysisTriggerScalerItem* s2 = GetTriggerScaler(runNumber,what2,triggerClassName2);
    
    if (!s1 || !s2) continue;
    
    x.push_back(runNumber);
    
    Float_t ratio(0);
    
    if ( s2->ValueCorrectedForDownscale() > 0 )
    {
      ratio = s1->ValueCorrectedForDownscale()/s2->ValueCorrectedForDownscale();
    }
    
    y.push_back(ratio);
    
      AliDebug(1,Form("RUN %09d %20s (%s) %12llu (%5d) %20s (%s) %12llu (%5d) R %7.2f",
                   runNumber,
                   triggerClassName1,what1,
                   s1->Value(),s1->DownscalingFactor(),
                   triggerClassName2,what2,
                   s2->Value(),s2->DownscalingFactor(),
                   ratio));
  }
  
  return new TGraph(x.size(),&x[0],&y[0]);
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTriggerRatioEvolution(const char* triggerClassName1,
                                                             const char* what1,
                                                             const char* triggerClassName2,
                                                             const char* what2)
{
  /// Plots the evolution of one trigger ratio
  
  Bool_t draw(kFALSE);
  Bool_t removeZeros(kFALSE);
  
  TGraph* g1 = PlotTriggerEvolution(triggerClassName1,what1,draw,0x0,removeZeros);
  TGraph* g2 = PlotTriggerEvolution(triggerClassName2,what2,draw,0x0,removeZeros);
  
  if (!g1 || !g2) return 0x0;
  
  if ( g1->GetN() != g2->GetN() )
  {
    AliError("Oups. Did not get the same number of points for the 2 graphs ?!");
    return 0x0;
  }

  Double_t x1,x2;
  Double_t y1,y2;
  Double_t x1err,x2err;
  Double_t y1err,y2err;
  
  TGraphErrors* g = new TGraphErrors(g1->GetN());
  Int_t j(0);
  
  for ( Int_t i = 0; i < g1->GetN(); ++i )
  {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    
    x1err = g1->GetErrorX(i);
    x2err = g2->GetErrorX(i);

    y1err = g1->GetErrorY(i);
    y2err = g2->GetErrorY(i);
    
    if  (x1!=x2)
    {
      AliError(Form("Points at %d don't have the same x ! ",i));
      continue;
    }
    
    Double_t y(0.0);
    
    if ( TMath::Abs(y2) > 1E-12 )
    {
      y = y1/y2;
    }
    
    Double_t yerr(0.0);
    
    if ( TMath::Abs(x1) > 1E-12)
    {
      yerr += TMath::Sqrt(x1err*x1err/(x1*x1));
    }

    if ( TMath::Abs(x2) > 1E-12)
    {
      yerr += TMath::Sqrt(x2err*x2err/(x2*x2));
    }

    yerr *= y;
      
    g->SetPoint(j,x1,y);
    g->SetPointError(j,x1err,yerr);
    
    ++j;
  }
  
  delete g1;
  delete g2;
  
  TimeAxis(g);
  
  return g;
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::Print(Option_t* /* opt */) const
{
  /// print our runlist
  AliAnalysisTriggerScalers::PrintIntegers(fRunList,',');
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::PrintIntegers(const std::vector<int>& integers,
                                              const char sep,
                                              std::ostream& out)
{
  /// print a list of integers
  for ( std::vector<int>::size_type i = 0; i < integers.size(); ++i )
  {
    out << integers[i] << sep;
  }
  out << std::endl;
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::ReadIntegers(const char* filename,
                                             std::vector<int>& integers,
                                             Bool_t resetVector)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
    std::ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
  std::set<int> runset;
  
  if (!resetVector)
  {
    for ( std::vector<int>::size_type j = 0; j < integers.size(); ++ j )
    {
      runset.insert(integers[j]);
    }
  }
  
  char line[10000];
  
  in.getline(line,10000,'\n');
  
  TString sline(line);
  
  if (sline.Contains(","))
  {
    TObjArray* a = sline.Tokenize(",");
    TIter next(a);
    TObjString* s;
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }
    delete a;
  }
  else
  {
    runset.insert(sline.Atoi());
    
    while ( in >> i )
    {
      runset.insert(i);
    }
  }
  
  integers.clear();
  
  for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it ) 
  {
    integers.push_back((*it)); 
  }
  
  std::sort(integers.begin(),integers.end());
}


//______________________________________________________________________________
void AliAnalysisTriggerScalers::SetRunList(Int_t runNumber)
{
  // Make the runlist be a single run
  fRunList.clear();
  fRunList.push_back(runNumber);
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::SetRunList(const std::vector<int>& runs)
{
  // Make the runlist be a single run
  fRunList = runs;
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::SetRunList(const std::set<int>& runs)
{
  // Make the runlist be a single run
  fRunList.clear();
  for ( std::set<int>::const_iterator it = runs.begin(); it != runs.end(); ++it )
  {
    fRunList.push_back(*it);
  }
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::SetRunList(const char* runlist)
{
  // Read the runlist from an ASCII file or a comma separated list
  // or a space separated list
  
  fRunList.clear();
  
  if ( TString(runlist).Contains(",") || TString(runlist).Contains(" ") )
  {
    TObjArray* runs = 0x0;
    if ( TString(runlist).Contains(",") ) 
    {
      runs = TString(runlist).Tokenize(",");
    }
    else
    {
      runs = TString(runlist).Tokenize(" ");
    }
    TIter next(runs);
    TObjString* s;
    std::set<int> runset;
    
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }

    for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it ) 
    {
      fRunList.push_back((*it)); 
    }
    
    std::sort(fRunList.begin(),fRunList.end());

    delete runs;
  }
  else
  {
    ReadIntegers(runlist,fRunList);
  }
}

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

ClassImp(AliAnalysisTriggerScalerItem)

//______________________________________________________________________________
const char* 
AliAnalysisTriggerScalerItem::BCMaskName() const 
{
  // return bc mask name
  if ( BCMask() ) return BCMask()->GetName(); else return ""; 
}

//______________________________________________________________________________
Int_t AliAnalysisTriggerScalerItem::Compare(const TObject* obj) const
{
  /// compare two scaler items (by means of their run number only)
  const AliAnalysisTriggerScalerItem* s = static_cast<const AliAnalysisTriggerScalerItem*>(obj);
  
  if ( s->RunNumber() < RunNumber() ) 
  {
    return -1;
  }
  else if ( s->RunNumber() > RunNumber() )
  {
    return 1;
  }
  return 0;
}

//______________________________________________________________________________
void AliAnalysisTriggerScalerItem::Print(Option_t* opt) const
{
  /// Printout
  TString sopt(opt);
  
  if ( RunNumber() > 0 )
  {
    sopt.ToUpper();
    
    if ( sopt.Contains("FULL") )
    {
    }
    else
    {
      std::cout << Form("RUN %6d CLASS %24s (%5s %4d BCs) SCALER %s %12llu DS %6d DURATION %ld",
                   RunNumber(),TriggerClassName(),
                   BCMaskName(),
                   BCMask() ? BCMask()->GetNUnmaskedBCs() : 0,
                   Level(),
                   Value(),DownscalingFactor(),Duration()) << std::endl;
    }
  }
  else
  {
    std::cout << Form("CLASS %24s SCALER %20llu NRUNS %d",
                 TriggerClassName(),Value(),NofRuns()) << std::endl;    
  }
}
