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
#include "TObjArray.h"
#include "TObjString.h"
#include <set>
#include <map>
#include "TStyle.h"
#include "AliLHCData.h"
#include "TMath.h"
#include "TAxis.h"
#include "TCanvas.h"

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
AliAnalysisTriggerScalers::AliAnalysisTriggerScalers(Int_t runNumber, const char* ocdbPath) : fRunList(), fVerbose(0), fOCDBPath(ocdbPath)
{
  // ctor using a single run number
  SetRunList(runNumber);
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::AliAnalysisTriggerScalers(const char* runlist, const char* ocdbPath) : fRunList(), fVerbose(0), fOCDBPath(ocdbPath)
{
  // ctor from a run list (txt file)
  SetRunList(runlist);
}

//_____________________________________________________________________________
AliAnalysisTriggerScalers::~AliAnalysisTriggerScalers()
{
}

//______________________________________________________________________________
TObject* AliAnalysisTriggerScalers::GetOCDBObject(const char* path, Int_t runNumber)
{
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
Int_t AliAnalysisTriggerScalers::GetTriggerInput(Int_t runNumber, const char* inputname)
{
  AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
  if (!tc) return -1;
  
  const TObjArray& inputs = tc->GetInputs();
  
  AliTriggerInput* ti = static_cast<AliTriggerInput*>(inputs.FindObject(inputname));
  
  if (!ti) return -1;
  
  return ti->GetSignature();
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
  
  TString diCurrent(Form("L3:%5.0f;DIP:%5.0f [L3:%d;DIP:%d]",
                         grp->GetL3Current((AliGRPObject::Stats)0),
                         grp->GetDipoleCurrent((AliGRPObject::Stats)0),
                         grp->GetL3Polarity(),
                         grp->GetDipolePolarity()));
  
  if (!tc || !trs || !grp) return 0x0;
  
  time_t duration = grp->GetTimeEnd() - grp->GetTimeStart();

  const TObjArray& trClasses = tc->GetClasses();
  
  const TObjArray* scalers = trs->GetScalersRecords();  
  const AliTriggerScalersRecord* begin = (AliTriggerScalersRecord*)(scalers->First());
  const AliTriggerScalersRecord* end = (AliTriggerScalersRecord*)(scalers->Last());
  
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
void AliAnalysisTriggerScalers::IntegratedLuminosity()
{
    const char* T0LUMITRIGGER = "C0TVX-S-NOPF-ALLNOTRD";
//    const char* V0LUMITRIGGER = "CVBAND-S-NOPF-ALLNOTRD";
    
    float T0LUMI=26.0; //mb
//    float V0LUMI=55.0; //mb
  
  double intLumi(0);
  
    std::map<std::string,float> lumiPerTrigger;
    std::map<std::string,time_t> durationPerTrigger;
    
  TObjArray triggerList;
  triggerList.SetOwner(kTRUE);
  
  triggerList.Add(new TObjString("CMUL8-S-NOPF-ALLNOTRD"));
  triggerList.Add(new TObjString("CMUL7-S-NOPF-ALLNOTRD"));
  triggerList.Add(new TObjString("CMUL8-S-NOPF-MUON"));
  triggerList.Add(new TObjString("CMUL7-S-NOPF-MUON"));
  triggerList.Add(new TObjString("C0MUL-SC-NOPF-MUON"));
//  triggerList.Add(new TObjString("C0MUL-SA-NOPF-MUON"));
  
    for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
    {
        Int_t runNumber = fRunList[i];
        
        // find out which trigger classes to use
        
        AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
        const TObjArray& trClasses = tc->GetClasses();

        AliGRPObject* grp = static_cast<AliGRPObject*>(GetOCDBObject("GRP/GRP/Data",runNumber));
        time_t duration = grp->GetTimeEnd() - grp->GetTimeStart();

      TIter nextTrigger(&triggerList);
      TObjString* trigger;
      
      while ( ( trigger = static_cast<TObjString*>(nextTrigger()) ) )
      {
        TString lumiTriggerClassName(T0LUMITRIGGER);
        TString muTriggerClassName(trigger->String());
        Int_t n(0);
        float lumiSigma = T0LUMI*1E6; //nb
        
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
        
        AliAnalysisTriggerScalerItem* lumiB = GetTriggerScaler(runNumber,"L0B",lumiTriggerClassName.Data());
        AliAnalysisTriggerScalerItem* muonA = GetTriggerScaler(runNumber,"L0A",muTriggerClassName.Data());
        AliAnalysisTriggerScalerItem* muonB = GetTriggerScaler(runNumber,"L0B",muTriggerClassName.Data());
        
        if (!lumiB || !muonA || !muonB) continue;
        
        Float_t ratio(0);
        
        if ( muonB->ValueCorrectedForDownscale() > 0 )
        {
            ratio = muonA->ValueCorrectedForDownscale()/muonB->ValueCorrectedForDownscale();
        }
        
//        muonA->Print();
//        muonB->Print();
        
        ratio *= lumiB->ValueCorrectedForDownscale()/lumiSigma;

        if ( muTriggerClassName.BeginsWith("CMUL")) {
          intLumi += ratio;
        }
        
        lumiPerTrigger[muTriggerClassName.Data()] += ratio;
        durationPerTrigger[muTriggerClassName.Data()] += duration;
        
      }
    }
  
    AliInfo(Form("Integrated lumi %7.4f nb^-1",intLumi));
    
    std::map<std::string,float>::const_iterator it;
    
    for ( it = lumiPerTrigger.begin(); it != lumiPerTrigger.end(); ++it )
    {
        AliInfo(Form("Trigger %30s Lumi %7.4f nb^-1 duration %10ld s",it->first.c_str(),it->second,durationPerTrigger[it->first]));
    }
  
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTriggerEvolution(const char* triggerClassName,
                                                        const char* what,
                                                        bool draw,
                                                        double* mean)
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
  
  std::vector<int> vx;
  std::vector<int> vex;
  std::vector<double> vy;
  std::vector<double> vey;
  
  TString swhat(what);
  swhat.ToUpper();
  
  if (mean) *mean=0;
  double nvalues(0);
  

  for ( std::vector<int>::size_type iRun = 0; iRun < fRunList.size(); ++iRun )
  {
    Int_t runNumber = fRunList[iRun];
    
    AliTriggerConfiguration* tc = static_cast<AliTriggerConfiguration*>(GetOCDBObject("GRP/CTP/Config",runNumber));
    AliTriggerRunScalers* trs = static_cast<AliTriggerRunScalers*>(GetOCDBObject("GRP/CTP/Scalers",runNumber));
    
    AliLHCData* lhc = static_cast<AliLHCData*>(GetOCDBObject("GRP/GRP/LHCData",runNumber));
    
    Int_t NumberOfInteractingBunches(0);
    Int_t NumberOfInteractingBunchesMeasured(0);
    Int_t NIBM2(0);
    
    int beam1(0);
    int beam2(1);
    
    AliLHCDipValI* val = lhc->GetBunchConfigDeclared(beam1,0);
    
    for ( Int_t i = 0; i < val->GetSizeTotal(); ++i )
    {
      if ( val->GetValue(i) < 0 ) ++NumberOfInteractingBunches;
    }
    
    AliLHCDipValI* valm = lhc->GetBunchConfigMeasured(beam1,0);
    
    for ( Int_t i = 0; i < valm->GetSizeTotal(); ++i )
    {
      if ( valm->GetValue(i) < 0 ) ++NumberOfInteractingBunchesMeasured;
    }
    
    valm = lhc->GetBunchConfigMeasured(beam2,0);
    
    for ( Int_t i = 0; i < valm->GetSizeTotal(); ++i )
    {
      if ( valm->GetValue(i) <= 0 ) ++NIBM2;
    }
    
    AliInfo(Form("RUN %09d NumberOfInteractingBunches=%3d NumberOfInteractingBunchesMeasured=%3d NIBM2=%3d",
                 runNumber,NumberOfInteractingBunches,NumberOfInteractingBunchesMeasured,NIBM2));
    
//    if ( NumberOfInteractingBunches <= 0 )
//    {
//      AliInfo(Form("RUN %09d NumberOfInteractingBunches=%3d !!!!",runNumber,NumberOfInteractingBunches));
//      return 0x0;
//    }
    
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
        
        if ( l0b <= 2 || ( l0a <= 2 && l0a != 0 ) || timelapse <= 9 ) continue;
        
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
        
        if (swhat.Contains("L0AOVERB") )
        {
          value = l0a*1.0/l0b;
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
          value = Mu(l0b/timelapse,NumberOfInteractingBunches);
        }
        else if ( swhat.Contains("PILEUPFACTOR") )
        {
          Double_t mu = Mu(l0b/timelapse,NumberOfInteractingBunches);
          value = mu/(1-TMath::Exp(-mu));
        }
        else if ( swhat.Contains("VSNB") )
        {
          value = l0b/(11245.0*NumberOfInteractingBunches);
        }
        else
        {
          value = timelapse;
          AliInfo(Form("Unknown what %s",what));
        }
        
        if ( ! swhat.Contains("OVER") && ! swhat.Contains("RATIO") &&
            ! swhat.Contains("MU") && ! swhat.Contains("PILEUPFACTOR") )
        {
          value /= timelapse;
        }
        
        if ( !TMath::Finite(value) ) continue;
        
        if (mean)
        {
          (*mean) += value;
          ++nvalues;
        }
        
        vx.push_back(seconds);
        vex.push_back(1);
        
        vy.push_back(value);
        
        vey.push_back(1);
        
      }
    }
    
    if (mean && nvalues)
    {
      (*mean) /= nvalues;
    }
    
    if ( vx.empty() ) return 0;
    
  }
  
  TGraph* g = MakeGraph(vx,vex,vy,vey);
  TString title(Form("TriggerEvolution%s-%s",triggerClassName,swhat.Data()));
  
  g->SetName(title.Data());
  g->SetTitle(title.Data());
  
  g->GetYaxis()->SetTitle(title.Data());
  
  if (draw)
  {
    NewCanvas(g->GetName());
    g->Draw("AL");
    g->SaveAs(Form("%s.pdf",title.Data()));
  }
  
  
  return g;
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::MakeGraph(const std::vector<int>& vx,
                                             const std::vector<int>& vex,
                                             const std::vector<double>& vy,
                                             const std::vector<double>& vey)
{
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
  
  TGraph* g = new TGraph(vx.size(),x,y);
  
  TimeAxis(g);
  
  return g;
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTrigger(const char* triggerClassName,
                                               const char* what)
{
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
    
    y.push_back(value);
  }
  
  if ( fRunList.size() ) {
    mean /= fRunList.size();
  }
  
  AliInfo(Form("Integral %e mean %e",integral,mean));
  
  return new TGraph(x.size(),&x[0],&y[0]);
}

//______________________________________________________________________________
TGraph* AliAnalysisTriggerScalers::PlotTriggerRatio(const char* triggerClassName1,
                                                    const char* what1,
                                                    const char* triggerClassName2,
                                                    const char* what2)
{
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
    
    if ( fVerbose ) 
    {
      AliInfo(Form("RUN %09d %20s (%s) %12llu (%5d) %20s (%s) %12llu (%5d) R %7.2f",
                   runNumber,
                   triggerClassName1,what1,
                   s1->Value(),s1->DownscalingFactor(),
                   triggerClassName2,what2,
                   s2->Value(),s2->DownscalingFactor(),
                   ratio));
    }
  }
  
  return new TGraph(x.size(),&x[0],&y[0]);
}

//______________________________________________________________________________
void AliAnalysisTriggerScalers::PrintIntegers(std::vector<int>& integers,
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
  if ( BCMask() ) return BCMask()->GetName(); else return ""; 
}

//______________________________________________________________________________
Int_t AliAnalysisTriggerScalerItem::Compare(const TObject* obj) const
{
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

