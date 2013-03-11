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

#include "AliAnalysisMuMu.h"

#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliAnalysisTriggerScalers.h"
#include "AliCounterCollection.h"
#include "AliHistogramCollection.h"
#include "AliLog.h"
#include "AliMergeableCollection.h"
#include "Riostream.h"
#include "TArrayL64.h"
#include "TASImage.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGrid.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cassert>
#include <map>
#include <set>
#include <string>

ClassImp(AliAnalysisMuMu)

TString AliAnalysisMuMu::fgOCDBPath("raw://");

TString AliAnalysisMuMu::fgDefaultDimuonTriggers("CMUL7-B-NOPF-MUON");

//,CMUL7-S-NOPF-ALLNOTRD,CMUL7-S-NOPF-MUON,CMUL8-S-NOPF-MUON,CMUL7-B-NOPF-ALLNOTRD,CMUU7-B-NOPF-ALLNOTRD,CMUU7-B-NOPF-MUON,CPBI1MUL-B-NOPF-MUON,CMULLO-B-NOPF-MUON");

TString AliAnalysisMuMu::fgDefaultMuonTriggers("CMSL7-S-NOPF-MUON");

//,CMSL7-S-NOPF-ALLNOTRD,CMSL8-S-NOPF-MUON,CMSL8-S-NOPF-ALLNOTRD,CMSL7-B-NOPF-MUON,CMUS1-B-NOPF-MUON,CMUS7-B-NOPF-MUON,CMSNGL-B-NOPF-MUON");

TString AliAnalysisMuMu::fgDefaultMinbiasTriggers("CINT7-B-NOPF-ALLNOTRD");

//,CINT7-S-NOPF-ALLNOTRD,CINT8-B-NOPF-ALLNOTRD,CINT8-S-NOPF-ALLNOTRD,CINT1-B-NOPF-ALLNOTRD,CPBI2_B1-B-NOPF-ALLNOTRD");

TString AliAnalysisMuMu::fgDefaultEventSelectionList("PSALL"); // for real data, for simulation see AliAnalysisMuMu ctor

TString AliAnalysisMuMu::fgDefaultPairSelectionList("pMATCHLOWRABSBOTH");

TString AliAnalysisMuMu::fgDefaultCentralitySelectionList("PP");

//TString AliAnalysisMuMu::fgDefaultFitTypeList("PSILOW:2,PSILOWalphaLow0.984nLow5.839alphaUp1.972nUp3.444:2,PSILOWMCTAILS:2");
TString AliAnalysisMuMu::fgDefaultFitTypeList("PSILOWalphaLow0.984nLow5.839alphaUp1.972nUp3.444:2,PSILOWMCTAILS:2");

TString AliAnalysisMuMu::fgDefaultEventSelectionForSimulations("ALL");
TString AliAnalysisMuMu::fgDefaultDimuonTriggerForSimulations("CMULLO-B-NOPF-MUON");

Bool_t AliAnalysisMuMu::fgIsCompactGraphs(kFALSE);

//_____________________________________________________________________________
TString First(const TString s)
{
  TString rv;
  
  TObjArray* tokens = s.Tokenize(",");
  
  if (!tokens) return rv;
  
  rv = static_cast<TObjString*>(tokens->First())->String();
  
  delete tokens;
  
  return rv;
}

//_____________________________________________________________________________
TString FindTrigger(const AliMergeableCollection& mc,
                    const char* base,
                    const char* selection,
                    const char* paircut,
                    const char* centrality)
{
  /// find the trigger containing the MinvPt histograms
  
  std::vector<std::string> trigger2test;
  
  //  trigger2test.push_back(Form("%s5-B-NOPF-ALLNOTRD",base));
  //  trigger2test.push_back(Form("%s1-B-NOPF-ALLNOTRD",base));
  //  trigger2test.push_back(Form("%s1B-ABCE-NOPF-MUON",base));
  if ( TString(base).Contains("||") || TString(base).Contains("-") )
  {
    trigger2test.push_back(base);
  }
  else
  {
    trigger2test.push_back(Form("%s-B-NOPF-ALLNOTRD",base));
    trigger2test.push_back(Form("%s-B-NOPF-MUON",base));
    trigger2test.push_back(Form("%s-S-NOPF-ALLNOTRD",base));
    trigger2test.push_back(Form("%s-S-NOPF-MUON",base));
  }
  trigger2test.push_back("ANY");
  
  for ( std::vector<std::string>::size_type i = 0; i < trigger2test.size(); ++i )
  {
    std::string trigger = trigger2test[i];
    
    if ( mc.GetObject(Form("/%s/%s/%s/%s",selection,trigger.c_str(),centrality,paircut),"MinvUS") ||
        mc.GetObject(Form("/%s/%s/%s/%s",selection,trigger.c_str(),centrality,paircut),"MinvUSPt")
        )
    {
      return trigger.c_str();
    }
  }
  
  //  AliWarningGeneral("FindTrigger",Form("DID NOT FIND TRIGGER base=%s selection=%s paircut=%s centrality=%s",
  //                  base,selection,paircut,centrality));
  //  for ( std::vector<std::string>::size_type i = 0; i < trigger2test.size(); ++i )
  //  {
  //    AliWarningGeneral("FindTrigger",Form("tested trigger = %s",trigger2test[i].c_str()));
  //  }
  return "";
}

//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu(const char* filename, const char* associatedSimFileName) : TObject(),
fFilename(filename),
fCounterCollection(0x0),
fDimuonTriggers(fgDefaultDimuonTriggers),
fMuonTriggers(fgDefaultMuonTriggers),
fMinbiasTriggers(fgDefaultMinbiasTriggers),
fEventSelectionList(fgDefaultEventSelectionList),
fPairSelectionList(fgDefaultPairSelectionList),
fCentralitySelectionList(fgDefaultCentralitySelectionList),
fFitTypeList(fgDefaultFitTypeList),
fBinning(0x0),
fMergeableCollection(0x0),
fRunNumbers(),
fCorrectionPerRun(0x0),
fAssociatedSimulation(0x0)
{
  // ctor
  
  GetCollections(fFilename,fMergeableCollection,fCounterCollection,fBinning,fRunNumbers);
  
  if (IsSimulation())
  {
    SetEventSelectionList("ALL");
    SetDimuonTriggerList("CMULLO-B-NOPF-MUON");
    SetFitTypeList("PSI1:1,COUNTJPSI:1");
  }
  
  if ( strlen(associatedSimFileName) )
  {
    fAssociatedSimulation = new AliAnalysisMuMu(associatedSimFileName);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMu::~AliAnalysisMuMu()
{
  // dtor
  delete fCounterCollection;
  delete fBinning;
  delete fMergeableCollection;
  delete fCorrectionPerRun;
  delete fAssociatedSimulation;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::BasicCounts(Bool_t detailTriggers,
                                  ULong64_t* totalNmb,
                                  ULong64_t* totalNmsl,
                                  ULong64_t* totalNmul)
{
  // Report of some basic numbers, like number of MB and MUON triggers, 
  // both before and after physics selection, and comparison with 
  // the total number of such triggers (taken from the OCDB scalers)
  // if requested.
  //
  // filename is assumed to be a root filecontaining a list containing
  //    an AliCounterCollection (or directly an AliCounterCollection)
  //
  // if detailTriggers is kTRUE, each kind of (MB,MUL,MSL) is counted separately
  //
  
  if (!fMergeableCollection || !fCounterCollection) return;
  
  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);

  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);

  TObjArray* events = fCounterCollection->GetKeyWords("event").Tokenize(",");

  Bool_t doPS = (events->FindObject("PSALL") != 0x0);
  
  TObjString* srun;
  TObjString* strigger;

  ULong64_t localNmb(0);
  ULong64_t localNmsl(0);
  ULong64_t localNmul(0);
  
  if ( totalNmb) *totalNmb = 0;
  if ( totalNmsl) *totalNmsl = 0;
  if ( totalNmul ) *totalNmul = 0;

  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    std::cout << Form("RUN %09d ",srun->String().Atoi());
    
    TString details;
    ULong64_t nmb(0);
    ULong64_t nmsl(0);
    ULong64_t nmul(0);
    
    nextTrigger.Reset();
    
    Int_t nofPS(0);
    
    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      
      if ( !fgDefaultMinbiasTriggers.Contains(strigger->String().Data()) &&
           !fgDefaultMuonTriggers.Contains(strigger->String().Data()) &&
           !fgDefaultDimuonTriggers.Contains(strigger->String().Data()) ) continue;
          
      ULong64_t n = fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                    strigger->String().Data(),"ALL",srun->String().Atoi()));

      details += TString::Format("\n%50s %10lld",strigger->String().Data(),n);
      

      ULong64_t nps = fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                      strigger->String().Data(),"PSALL",srun->String().Atoi()));

      if ( doPS )
      {
        details += TString::Format(" PS %5.1f %%",nps*100.0/n);
      }

      if (nps)
      {
        ++nofPS;
      }
      
      if ( fMinbiasTriggers.Contains(strigger->String()) )
      {
        nmb += n;
        if ( totalNmb) (*totalNmb) += n;
        localNmb += n;
      }
      else if ( fMuonTriggers.Contains(strigger->String()) )
      {
        nmsl += n;
        if ( totalNmsl) (*totalNmsl) += n;
        localNmsl += n;
      }
      else if ( fDimuonTriggers.Contains(strigger->String()) )
      {
        nmul += n;
        if ( totalNmul ) (*totalNmul) += n;
        localNmul += n;
      }      
    }
    
    std::cout << Form("MB %10lld MSL %10lld MUL %10lld %s",
                 nmb,nmsl,nmul,(nofPS == 0 ? "(NO PS AVAIL)": ""));
    
    if ( detailTriggers )
    {
      std::cout << details.Data();
    }
    std::cout << std::endl;
  }

  if ( !totalNmul && !totalNmsl && !totalNmb )
  {
    std::cout << std::endl << Form("%13s MB %10lld MSL %10lld MUL %10lld ","TOTAL",
                                   localNmb,localNmsl,localNmul) << std::endl;
  }

  delete runs;
  delete triggers;
  delete events;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::BasicCountsEvolution(const char* filelist, Bool_t detailTriggers)
{
  // Report of some basic numbers, like number of MB and MUON triggers,
  // both before and after physics selection, and comparison with
  // the total number of such triggers (taken from the OCDB scalers)
  // if requested.
  //
  // if detailTriggers is kTRUE, each kind of (MB,MUL,MSL) is counted separately
  //
  // To change the list of (single muon, dimuon, MB) triggers, use
  // the SetDefault*TriggerList methods prior to call this one
  //
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return;
  
  TIter next(files);
  TObjString* str;
  
  ULong64_t totalNmb(0);
  ULong64_t totalNmsl(0);
  ULong64_t totalNmul(0);

  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliAnalysisMuMu m(str->String().Data());
    
    ULong64_t nmb(0);
    ULong64_t nmsl(0);
    ULong64_t nmul(0);
    
    m.BasicCounts(detailTriggers,&nmb,&nmsl,&nmul);
    
    totalNmb += nmb;
    totalNmsl += nmsl;
    totalNmul += nmul;
  }
  
  std::cout << std::endl << Form("%13s MB %10lld MSL %10lld MUL %10lld ","TOTAL",
                                 totalNmb,totalNmsl,totalNmul) << std::endl;

}

//_____________________________________________________________________________
void AliAnalysisMuMu::CentralityCheck(const char* filelist)
{
  // Check if we get correctly filled centrality
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return;
  
  TIter next(files);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliMergeableCollection* mc(0x0);
    AliCounterCollection* cc(0x0);
    AliAnalysisMuMuBinning* bin(0x0);
    std::set<int> runnumbers;
    
    if (!GetCollections(str->String().Data(),mc,cc,bin,runnumbers)) continue;
    
    int run = RunNumberFromFileName(str->String().Data());
    
    TH1* h = mc->Histo("/ALL/CPBI1MUL-B-NOPF-MUON/Centrality");
    
    float percent(0);
    
    if (h)
    {
      percent = 100*h->Integral(1,1)/h->Integral();
    }
    
    std::cout << Form("RUN %09d PERCENT %7.2f",run,percent) << std::endl;
    
    delete mc;
  }
  
  gROOT->CloseFiles();
  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::CleanAllSpectra()
{
  /// Delete all the spectra we may have

  MC()->RemoveByType("AliAnalysisMuMuSpectra");
  Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Compact(TGraph& g)
{
  /// Compact (i.e. get the equivalent of 1 bin = 1 run number for an histogram)
  /// the graph.
  /// Only works if the x content of this graph represents run numbers. Otherwise
  /// result is unpredictable ;-)
  
  if ( !g.GetN() ) return;
  
  TGraph* newgraph = static_cast<TGraph*>(g.Clone());
  Double_t x,xerr,y,yerr;
  TGraphErrors* ge = dynamic_cast<TGraphErrors*>(newgraph);
  
  TAxis* axis = g.GetXaxis();
  
  for ( Int_t i = 0; i < newgraph->GetN(); ++i )
  {
    g.GetPoint(i,x,y);
    if (ge)
    {
      xerr = ge->GetErrorX(i);
      yerr = ge->GetErrorY(i);
    }
    
    g.SetPoint(i,i+0.5,y);
    if (ge)
    {
      static_cast<TGraphErrors&>(g).SetPointError(i,0.5,yerr);
    }
    
    axis->SetBinLabel(i,Form("%d",TMath::Nint(x)));
  }
  
  
}


//_____________________________________________________________________________
TObjArray* AliAnalysisMuMu::CompareJpsiPerCMUUWithBackground(const char* jpsiresults,
                                                                   const char* backgroundresults)
{
  TFile* fjpsi = FileOpen(jpsiresults);
  TFile* fbck = FileOpen(backgroundresults);
  
  if (!fjpsi || !fbck) return 0x0;
  
  TGraph* gjpsi = static_cast<TGraph*>(fjpsi->Get("jpsipercmuu"));
    
  std::vector<std::string> checks;

  checks.push_back("muminus-CMUU7-B-NOPF-ALLNOTRD");
  checks.push_back("muplus-CMUU7-B-NOPF-ALLNOTRD");
  checks.push_back("muminus-CMUSH7-B-NOPF-MUON");
  checks.push_back("muplus-CMUSH7-B-NOPF-MUON");
  
  if (!gjpsi) return 0x0;

  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  for ( std::vector<std::string>::size_type j = 0; j < checks.size(); ++j )
  {
    
    TGraph* gback = static_cast<TGraph*>(fbck->Get(checks[j].c_str()));
    
    if (!gback) continue;

    if ( gjpsi->GetN() != gback->GetN() )
    {
      AliErrorClass("graphs have different number of points !");
      continue;
    }
    
    TGraphErrors* g = new TGraphErrors(gjpsi->GetN());
    
    for ( int i = 0; i < gjpsi->GetN(); ++i ) 
    {
      double r1,r2,y1,y2;
      
      gjpsi->GetPoint(i,r1,y1);
      gback->GetPoint(i,r2,y2);
      
      if ( r1 != r2 ) 
      {
        AliWarningClass(Form("run[%d]=%d vs %d",i,(int)r1,(int)r2));
        continue;
      }
      
      g->SetPoint(i,y2,y1);
      //    g->SetPointError(i,gjpsi->GetErrorY(i),gback->GetErrorY(i));
    }
    
    g->SetMarkerStyle(25+j);
    g->SetMarkerSize(1.2);
    if (j==0)
    {
      g->Draw("ap");
    }
    else
    {
      g->Draw("p");
    }
    g->SetLineColor(j+1);
    g->SetMarkerColor(j+1);
    g->SetName(checks[j].c_str());
    a->AddLast(g);
  }
  
  return a;
}

//_____________________________________________________________________________
TGraph* AliAnalysisMuMu::CompareJpsiPerCMUUWithSimu(const char* realjpsiresults,
                                                             const char* simjpsiresults)
{
  TFile* freal = FileOpen(realjpsiresults);
  TFile* fsim = FileOpen(simjpsiresults);
  
  if (!freal || !fsim) return 0x0;
  
  TGraph* greal = static_cast<TGraph*>(freal->Get("jpsipercmuu"));
  TGraph* gsim = static_cast<TGraph*>(fsim->Get("jpsipercmuu"));
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  if ( greal->GetN() != gsim->GetN() )
  {
    AliErrorClass("graphs have different number of points !");
    return 0x0;
  }
    
  TGraphErrors* g = new TGraphErrors(greal->GetN());
  TGraphErrors* gratio = new TGraphErrors(greal->GetN());
    
  for ( int i = 0; i < greal->GetN(); ++i ) 
  {
    double r1,r2,y1,y2;
    
    greal->GetPoint(i,r1,y1);
    gsim->GetPoint(i,r2,y2);
    
    if ( r1 != r2 ) 
    {
      AliWarningClass(Form("run[%d]=%d vs %d",i,(int)r1,(int)r2));
      continue;
    }
    
    double ratio(0.0);
    
    if ( TMath::Abs(y1)<1E-6 || TMath::Abs(y2)<1E-6)
    {
      g->SetPoint(i,0,0);
      g->SetPointError(i,0,0);
    }
    else
    {    
      g->SetPoint(i,y2,y1);
      g->SetPointError(i,greal->GetErrorY(i),gsim ->GetErrorY(i));
      ratio = y2/y1;
    }
    gratio->SetPoint(i,r1,ratio);
  }
    
  g->SetMarkerStyle(25);
  g->SetMarkerSize(1.2);

  new TCanvas;
  
  g->Draw("ap");

  g->SetLineColor(1);
  g->SetMarkerColor(1);
  g->SetName("jpsipercmuurealvssim");

  new TCanvas;
  
  greal->Draw("alp");
  gsim->SetLineColor(4);
  
  gsim->Draw("lp");

  new TCanvas;
  gratio->Draw("alp");
  
  return g;
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMu::ComputeBackgroundEvolution(const char* filelist, 
                                                       const char* triggerList,
                                                       Double_t ptmin,
                                                       const char* outputFile,
                                                       const char* outputMode)
{
  // triggerList is a list of complete trigger names, separated by space
  // of the triggers to consider : only the first one found in the list 
  // is used for each run.
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return 0x0;
  
  TIter next(files);
  TObjString* str;
  
  const char* ps = "PSALL";
  const char* centrality = "PP";
  const char* ts1 = "sMATCHLOWRABS";
  const char* ts2 = "sMATCHLOWRABSDCA";
  
  std::map<std::string, std::vector<float> > runs;
  std::map<std::string, std::vector<float> > errruns;
  std::map<std::string, std::vector<float> > yplus,erryplus;
  std::map<std::string, std::vector<float> > yminus,erryminus;
  
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  
  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);

  Bool_t bothSigns(kFALSE);
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    AliInfoClass(str->String().Data());
    
    AliMergeableCollection* mc(0x0);
    AliCounterCollection* cc(0x0);
    AliAnalysisMuMuBinning* bin(0x0);
    std::set<int> runnumbers;
    
    if (!GetCollections(str->String().Data(),mc,cc,bin,runnumbers)) continue;
    
    TIter nextObject(mc->CreateIterator());
    TObject* o;
    int nplus(0), nminus(0);
    
    while ( (  o = nextObject() ) )
    {
      if ( o->InheritsFrom("TH1") )
      {
        continue;
      }
      
      TH1* h = static_cast<TH1*>(o);
      
      if ( TString(h->GetName()).EndsWith("Plus") )
      {
        nplus++;
      }
      if ( TString(h->GetName()).EndsWith("Minus") )
      {
        nminus++;
      }
    }
    
    if  (nminus==nplus && nplus>0 ) 
    {
      bothSigns = kTRUE;
    }

    AliInfoClass(Form("Both signs = %d",bothSigns));
    
    TIter nextTrigger(triggers);
    TObjString* trigger;
    TH1* h1p(0x0);
    TH1* h1m(0x0);
    TH1* h2p(0x0);
    TH1* h2m(0x0);
    TString format;
    
    while ( ( trigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      if  (bothSigns)
      {
        format = "/%s/%s/%s/%s/PtEtaMuPlus:py";
      }
      else
      {
        format = "/%s/%s/%s/%s/PtEtaMu:py";
      }
      
      TString hname(Form(format.Data(),ps,trigger->String().Data(),centrality,ts1));
      
      h1p = mc->Histo(hname.Data());
      
      if (!h1p)
      {
        AliInfoClass(Form("Could not get %s",hname.Data()));
        continue;
      }
      
      AliInfoClass(Form("Will use trigger %s",trigger->String().Data()));
      
      h2p = mc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts2));
      
      if ( bothSigns )
      {
        format.ReplaceAll("Plus","Minus");
        h1m = mc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts1));
        h2m = mc->Histo(Form(format.Data(),ps,trigger->String().Data(),centrality,ts2));
      }
      else
      {
        h2m=h2p;
        h1m=h1p;
      }
      
      if (h1m && h2m && h1p && h2p)
      {
        Int_t bin1 = h1m->GetXaxis()->FindBin(ptmin);
        Int_t bin2 = h1m->GetXaxis()->GetNbins();

        runs[trigger->String().Data()].push_back(RunNumberFromFileName(str->String().Data()));
        errruns[trigger->String().Data()].push_back(0.5);

        double e1,e2;
        double v1 = h2m->IntegralAndError(bin1,bin2,e1);
        double v2 = h1m->IntegralAndError(bin1,bin2,e2);
        double value = 100*(1.0-v1/v2);
        e1/=v1;
        e2/=v2;
        yminus[trigger->String().Data()].push_back(value);
//        double e1 = 1.0/TMath::Sqrt(h1m->GetEntries());
//        double e2 = 1.0/TMath::Sqrt(h2m->GetEntries());
        erryminus[trigger->String().Data()].push_back(TMath::Sqrt(e1*e1+e2*e2)*value);
        
        v1=h2p->IntegralAndError(bin1,bin2,e1);
        v2=h1p->IntegralAndError(bin1,bin2,e1);
        value = 100*(1.0-v1/v2);
        e1/=v1;
        e2/=v2;
        yplus[trigger->String().Data()].push_back(value);
//        e1 = 1.0/TMath::Sqrt(h1p->GetEntries());
//        e2 = 1.0/TMath::Sqrt(h2p->GetEntries());
        erryplus[trigger->String().Data()].push_back(TMath::Sqrt(e1*e1+e2*e2)*value);
      }
      else
      {
        std::cout << Form("Error : h1m %p h2m %p h1p %p h2p %p",h1m,h2m,h1p,h2p) << std::endl;
      }
    }
    
    delete mc;
    delete cc;
    TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(str->String().Data()));    
    delete f;
  }
  
  delete triggers;
  
  TFile* f = new TFile(outputFile,outputMode);
  
  std::map<std::string, std::vector<float> >::const_iterator it;
  
  for ( it = runs.begin(); it != runs.end(); ++it )
  {
    std::string triggerName = it->first;
    
  TGraphErrors* gp = new TGraphErrors(runs[triggerName].size(),&runs[triggerName][0],&yplus[triggerName][0],&errruns[triggerName][0],&erryplus[triggerName][0]);
  TGraphErrors* gm(0x0);
  
  if ( bothSigns ) 
  {
    gm = new TGraphErrors(runs[triggerName].size(),&runs[triggerName][0],&yminus[triggerName][0],&errruns[triggerName][0],&erryminus[triggerName][0]);
  }
  
  if ( bothSigns ) 
  {
    gp->Write(Form("muplus_%s",triggerName.c_str()),TObject::kOverwrite);
    gm->Write(Form("muminus_%s",triggerName.c_str()),TObject::kOverwrite);
  }
  else
  {
    gp->Write(Form("mu_%s",triggerName.c_str()),TObject::kOverwrite);
  }
  
  }
  
  delete f;
  
  return a;
}

//_____________________________________________________________________________
TMap*
AliAnalysisMuMu::ComputeJpsiEvolution(const char* filelist, const char* triggerList,
                                      const char* outputFile)
{
  /// Compute some jpsi information for a list of files / trigger combinations
  
  TObjArray* files = ReadFileList(filelist);
  
  if (!files || files->IsEmpty() ) return 0x0;
  
  TMap results; // one TObjString->TObjArray per file
  results.SetOwnerKeyValue(kTRUE,kTRUE);
  
  TIter nextFile(files);
  TObjString* str;
  TString fitType;
  
//  while ( ( str = static_cast<TObjString*>(nextFile()) ) )
//  {
//    std::cout << str->String().Data() << std::endl;
//    
//    AliAnalysisMuMu m(str->String().Data());
//    
//    m.SetDimuonTriggerList(triggerList);
//        
//    TMap* map = m.Jpsi();
//    
//    if (!map)
//    {
//      AliWarningClass(Form("Got no jpsi for %s",str->String().Data()));
//    }
//    
//    results.Add(new TObjString(str->String()), map);
//  }
//  
//  if (!results.GetSize()) return 0x0;
  
  // compute the total over all files
  
  TMap* total = new TMap;
  total->SetOwnerKeyValue(kTRUE,kTRUE);
  
  nextFile.Reset();
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  
  TIter nextTrigger(triggers);
  TObjString* trigger(0x0);
  
  while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
  {
    nextFile.Reset();
    
    TList l;
    AliAnalysisMuMuResult* ref(0x0);
    
    while ( ( str = static_cast<TObjString*>(nextFile()) ) )
    {
//      TObjArray* a = static_cast<TObjArray*>(results.GetValue(str->String().Data()));
      
      AliInfoClass("FIXME: write the merging of AliAnalysisMuMuResult objects !");
      
//      AliAnalysisMuMuResult* r(0x0);
//      
//      if (a)
//      {
//        r = static_cast<AliAnalysisMuMuResult*>(a->FindObject(trigger->String().Data()));
//        
//        if (r)
//        {
//          if (!ref) ref = r;
//
//          if ( !hminv )
//          {
//            TH1* htmp = static_cast<TH1*>(r->Minv()->Clone());
//            if (!htmp)
//            {
//              continue;
//            }
//            hminv = htmp;
//          }
//          else
//          {
//            l.Add(r->Minv());
//          }
//          
//          n += r->NofTriggers();
//          ++nruns;
//        }
//      }
    }
    
//    sum->Merge(&l);
    
    if (!ref) continue;
    

//    AliAnalysisMuMuResult* sum = new AliAnalysisMuMuResult(*hminv,
//                                                           ref->TriggerClass(),
//                                                           ref->EventSelection(),
//                                                           ref->PairSelection(),
//                                                           ref->CentralitySelection(),
//                                                           AliAnalysisMuMuBinning::Range());
//
//    sum->SetNofTriggers(n);
//    
//    sum->SetNofRuns(nruns);
//    
//    sum->Fit(1);
//    
//    total->Add(new TObjString(trigger->String().Data()),sum);
    
  }

  AliInfoClass("--------------------------------------");
  StdoutToAliInfoClass(total->Print(););

  AliInfoClass("---------------------------Going to write file");

  TFile* fout = new TFile(outputFile,"RECREATE");
  
  results.Write("rbr",TObject::kSingleKey);
  
  total->Write("total",TObject::kSingleKey);

  fout->Close();
  
  delete fout;
  
  AliInfoClass(Form("%d files analyzed",files->GetEntries()));
  
  return total;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::DecodeFileName(const char* filename,
                                             TString& period,
                                             int& esdpass,
                                             int& aodtrain,
                                             int& runnumber)
{
  esdpass=aodtrain=runnumber=-1;
  period="";
  
  TString sfile(gSystem->BaseName(filename));
  
  if (!sfile.BeginsWith("LHC") && !sfile.BeginsWith("SIM") ) 
  {
    std::cerr << Form("filename %s does not start with LHC or SIM",filename) << std::endl;
    return kFALSE;
  }
  
  int year;
  char p;
  Bool_t ok(kFALSE);
  
  if ( sfile.BeginsWith("LHC") ) 
  {
    if ( sfile.Contains("pass") && sfile.Contains("AODMUON") )
    {
      int pass;
      sscanf(sfile.Data(),"LHC%2d%c_pass%d_AODMUON%03d_%09d",&year,&p,&pass,&aodtrain,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("pass") && sfile.Contains("_muon_") && sfile.Contains("AOD000") )
    {
      // LHC11c_pass2_muon_AOD000_000152087.saf.root
      sscanf(sfile.Data(),"LHC%2d%c_pass%d_muon_AOD000_%09d",&year,&p,&esdpass,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;      
    }
    else if ( sfile.Contains("_muon_calo_") && sfile.Contains("AODMUON000") )
    {
      //      LHC12h_muon_calo_AODMUON000_000190112.saf.root
      sscanf(sfile.Data(),"LHC%2d%c_muon_calo_AODMUON000_%09d",&year,&p,&runnumber);
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("_muon_calo_") && sfile.Contains("AOD000") )
    {
      //      LHC12h_muon_calo_AOD000_000190112.saf.root
      sscanf(sfile.Data(),"LHC%2d%c_muon_calo_AOD000_%09d",&year,&p,&runnumber);
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("AODMUON" ) )
    {    
      sscanf(sfile.Data(),"LHC%2d%c_AODMUON%03d_%09d",&year,&p,&aodtrain,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("AOD000") ) 
    {
      sscanf(sfile.Data(),"LHC%2d%c_muons_AOD000_%09d",&year,&p,&runnumber); 
      period = TString::Format("LHC%2d%c",year,p);
      ok=kTRUE;
    }
    else if ( sfile.Contains("ESD_OUTER000"))
    {
      sscanf(sfile.Data(),"LHC%2d%c_cpass1_ESD_OUTER000_%09d",&year,&p,&runnumber);
      ok=kTRUE;
    }
  }
  else if ( sfile.BeginsWith("SIM_JPSI3" ) )
  {
    sscanf(sfile.Data(),"SIM_JPSI3_%09d",&runnumber);
    ok = kTRUE;
  }
  else if ( sfile.BeginsWith("SIM_UPSILON" ) )
  {
    sscanf(sfile.Data(),"SIM_UPSILON_%09d",&runnumber);
    ok = kTRUE;
  }
  else if ( sfile.BeginsWith("SIM_JPSI" ) )
  {
    sscanf(sfile.Data(),"SIM_JPSI_LHC%2d%c_%09d",&year,&p,&runnumber);
    period = TString::Format("LHC%2d%c",year,p);
    ok = kTRUE;
  }
  
  if (!ok)
  {
    std::cerr << Form("Can not decode %s",filename) << std::endl;
    return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DrawMinv(const char* type,
                               const char* particle,
                               const char* trigger,
                               const char* eventType,
                               const char* pairCut,
                               const char* centrality,
                               const char* subresultname,
                               const char* flavour) const
{
  /// Draw minv spectra for binning of given type
  
  if (!MC() || !BIN()) return;
  
  TObjArray* bins = BIN()->CreateBinObjArray(particle,type,flavour);
  if (!bins)
  {
    AliError(Form("Could not get %s bins",type));
    return;
  }
  
  Double_t xmin(-1);
  Double_t xmax(-1);
  
  TString sparticle(particle);
  if ( sparticle=="PSI" )
  {
    xmin = 2;
    xmax = 6;
  }
  
  Int_t nx(1);
  Int_t ny(1);
  
  Int_t n = bins->GetEntries();
  
  if ( n == 2 )
  {
    nx = 2;
  }
  else if ( n > 2 )
  {
    ny = TMath::Nint(TMath::Sqrt(n));
    nx = n/ny;
  }
  
  TString stype(type);
  stype.ToUpper();
  
  TString spectraName(Form("/%s/%s/%s/%s/%s-%s-%s",eventType,trigger,centrality,pairCut,particle,stype.Data(),flavour));
  
  AliAnalysisMuMuSpectra* spectra = static_cast<AliAnalysisMuMuSpectra*>(MC()->GetObject(spectraName.Data()));
  
  AliDebug(1,Form("spectraName=%s spectra=%p",spectraName.Data(),spectra));

  TObjArray* spectraBins(0x0);
  if ( spectra )
  {
    spectraBins = spectra->Bins();
  }
  
  TCanvas* c = new TCanvas;
  c->Divide(nx,ny);
  c->Draw();
  gStyle->SetOptFit(1112);
  
  c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "AliAnalysisMuMu",
              (void*)this, "ExecuteCanvasEvent(Int_t,Int_t,Int_t,TObject*)");

  
  TIter next(bins);
  AliAnalysisMuMuBinning::Range* r;
  Int_t ci(0);
  
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    TString name(Form("/%s/%s/%s/%s/MinvUS%s",eventType,trigger,centrality,pairCut,r->AsString().Data()));

    AliDebug(1,name.Data());
    
    AliAnalysisMuMuResult* spectraBin(0x0);
    
    if ( spectraBins )
    {
      AliAnalysisMuMuResult* sr = static_cast<AliAnalysisMuMuResult*>(spectraBins->At(ci));
      
      spectraBin = sr->SubResult(subresultname);
      
      AliDebug(1,Form("spectraBin(%s)=%p",subresultname,spectraBin));
    }
    
    TH1* h = MC()->Histo(name.Data());
    
    if ( spectraBin )
    {
      h = spectraBin->Minv();
    }
    
    if (h)
    {
      ++ci;
      c->cd(ci);
      gPad->SetLogy();
      if (xmin>0)
      {
        h->GetXaxis()->SetRangeUser(xmin,xmax);
      }
      h->Draw("histes");
      
      TObject* f1 = h->GetListOfFunctions()->FindObject("fitTotal");
      if (f1)
      {
        f1->Draw("same");
      }
      
      gPad->Modified();
      gPad->Update();
      
      TObject* stats = h->FindObject("stats");
      if (stats)
      {
        stats->Draw("same");
      }
    }
  }
  
  delete bins;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DrawMinv(const char* type, const char* particle, const char* flavour, const char* subresultname) const
{
  /// Draw minv spectra for binning of given type

  DrawMinv(type,particle,           
           First(DimuonTriggerList()).Data(),
           First(EventSelectionList()).Data(),
           First(PairSelectionList()).Data(),
           First(CentralitySelectionList()).Data(),
           subresultname,
           flavour);
}

//___________________________________________________________________
void AliAnalysisMuMu::ExecuteCanvasEvent(Int_t event, Int_t /*px*/, Int_t /*py*/, TObject *sel)
{
  // Actions in reponse to mouse button events.
  
  TCanvas* c = static_cast<TCanvas*>(gTQSender);
  TPad* pad = static_cast<TPad*>(c->GetSelectedPad());
  if (!pad) return;
  
//  if ((event == kButton1Down) ||
  if (event == kButton1Double) 
  {
    
//    Float_t x = pad->AbsPixeltoX(px);
//    Float_t y = pad->AbsPixeltoY(py);
//    x = pad->PadtoX(x);
//    y = pad->PadtoY(y);

//    std::cout << "event=" << event << " px=" << px << " py=" << py << " ";
    
    if ( sel && sel->InheritsFrom("TH1") )
    {
      TCanvas* clocal = new TCanvas;
      clocal->SetLogy();
      clocal->Draw();
      sel->Draw();
    }
    else
    {
      TList* list = pad->GetListOfPrimitives();
      TIter next(list);
      TObject* h;
      
      while ( ( h = next() ) )
      {
        if ( h->InheritsFrom("TH1") )
        {
          TCanvas* clocal = new TCanvas;
          clocal->SetLogy();
          clocal->Draw();
          h->Draw();
          break;
        }
      }
      
    }

//      std::cout  << std::endl;

      pad->Modified();
  }
  
}

//_____________________________________________________________________________
TString 
AliAnalysisMuMu::ExpandPathName(const char* file)
{
  // An expand method that lives alien URL as they are
  TString sfile;
  
  if ( !sfile.BeginsWith("alien://") )
  {
    return gSystem->ExpandPathName(file);
  }
  else
  {
    if (!gGrid) TGrid::Connect("alien://");
    if (!gGrid) return "";    
  }
  
  return file;
}

//_____________________________________________________________________________
TFile* 
AliAnalysisMuMu::FileOpen(const char* file)
{
  // Open a file after expansion of its name
  
  return TFile::Open(ExpandPathName(file).Data());
}

//_____________________________________________________________________________
TString AliAnalysisMuMu::First(const TString& list) const
{
  TObjArray* a = list.Tokenize(",");
  if ( a->GetLast() < 0 ) return "";
  
  TString rv = static_cast<TObjString*>(a->First())->String();
  
  delete a;
  
  return rv;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra*
AliAnalysisMuMu::FitParticle(const char* particle,
                             const char* trigger,
                             const char* eventType,
                             const char* pairCut,
                             const char* centrality,
                             const AliAnalysisMuMuBinning& binning)
{
  // Fit the minv spectra to find the given particle
  // Returns an array of AliAnalysisMuMuResult objects
  
  static int n(0);
  
  TObjArray* bins = binning.CreateBinObjArray(particle);
  if (!bins)
  {
    AliError(Form("Did not get any bin for particle %s",particle));
    return 0x0;
  }
  
  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  if ( !triggers->FindObject(trigger) )
  {
    AliDebug(1,Form("Did not find trigger %s",trigger));
    delete bins;
    delete triggers;
    return 0x0;
  }
  delete triggers;
  
  TObjArray* events = fCounterCollection->GetKeyWords("event").Tokenize(",");
  if ( !events->FindObject(eventType) )
  {
    AliError(Form("Did not find eventType %s",eventType));
    delete bins;
    delete events;
    return 0x0;
  }
  delete events;

  Int_t ntrigger = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s",trigger,eventType)));
  
  if  (ntrigger<=0)
  {
    AliError(Form("No trigger for trigger:%s/event:%s",trigger,eventType));
    delete bins;
    return 0x0;
  }
  
//  binning.Print();
  
  AliAnalysisMuMuSpectra* spectra(0x0);
  
  AliAnalysisMuMuBinning::Range* bin;
  TIter next(bins);
  
  TObjArray* fitTypeArray = fFitTypeList.Tokenize(",");
  TIter nextFitType(fitTypeArray);
  TObjString* fitType;
  TString flavour;
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    TString hname(Form("MinvUS%s",bin->AsString().Data()));
    
    TH1* hminv = fMergeableCollection->Histo(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut),hname.Data());
    
    if (!hminv)
    {
      if (!fBinning && bin->IsNullObject() )
      {
        // old file, we only had MinvUSPt
        hminv = fMergeableCollection->Histo(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut),"MinvUSPt:py");
      }
      
      if (!hminv)
      {
        AliDebug(1,Form("Could not find histo %s",hname.Data()));
        continue;
      }
    }
    
    hminv = static_cast<TH1*>(hminv->Clone(Form("minv%d",n++)));
    
    
    AliAnalysisMuMuResult* r = new AliAnalysisMuMuResult(*hminv,
                                                         trigger,
                                                         eventType,
                                                         pairCut,
                                                         centrality,
                                                         *bin);
    
    r->SetNofTriggers(ntrigger);
    
    nextFitType.Reset();
    
    while ( ( fitType = static_cast<TObjString*>(nextFitType())) )
    {
      AliDebug(1,Form("<<<<<< fitType=%s bin=%s",fitType->String().Data(),bin->Flavour().Data()));
      
      if ( fitType->String().BeginsWith("PSILOWMCTAILS") )
      {
        std::vector<Double_t> par;
        par = GetMCCB2Tails(*bin);
        if (!par.empty())
        {
          r->AddFit(fitType->String().Data(),par.size(),&par[0]);
        }
      }
      else
      {
        r->AddFit(fitType->String());
      }
    }
  
    flavour = bin->Flavour();
    
    if (!spectra)
    {
      TString spectraName(binning.GetName());
      if ( flavour.Length() > 0 )
      {
        spectraName += "-";
        spectraName += flavour;
      }
      spectra = new AliAnalysisMuMuSpectra(spectraName.Data());
    }
    
    spectra->AdoptResult(*bin,r);
    
    if ( IsSimulation() )
    {
      SetNofInputParticles(*r);
    }
  
    
  } // loop on bins
  
  delete fitTypeArray;
  delete bins;
  
  return spectra;
}

//_____________________________________________________________________________
std::vector<Double_t>
AliAnalysisMuMu::GetMCCB2Tails(const AliAnalysisMuMuBinning::Range& bin) const
{
  /// Get the tails from the associated simulation
  
  std::vector<Double_t> par;
  
  if (!SIM())
  {
    AliError("Cannot get MC tails without an associated simulation file !");
    return par;
  }


  AliAnalysisMuMuSpectra* s = static_cast<AliAnalysisMuMuSpectra*>(SIM()->GetSpectra(bin.Type().Data(),bin.Flavour().Data()));
  
  if (!s)
  {
    AliError(Form("Could not find spectra %s,%s for associated simulation",bin.Type().Data(),bin.Flavour().Data()));
    fAssociatedSimulation->MC()->Print("*:Ali*");
    return par;
  }
  else
  {
    AliDebug(1,Form("AliAnalysisMuMuSpectra* s = reinterpret_cast<AliAnalysisMuMuSpectra*>(%p)",s));
  }
  
  AliAnalysisMuMuResult* r = s->GetResultForBin(bin);
  
  if ( r )
  {
    AliAnalysisMuMuResult* r1 = r->SubResult("JPSI:1");
    if  (r1)
    {
      TF1* func = static_cast<TF1*>(r1->Minv()->GetListOfFunctions()->FindObject("fitTotal"));
      if (func)
      {
        par.push_back(func->GetParameter("alphaLow"));
        par.push_back(func->GetParameter("nLow"));
        par.push_back(func->GetParameter("alphaUp"));
        par.push_back(func->GetParameter("nUp"));
      }
    }
  }
  
  return par;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisMuMu::GetTriggerScalerCount(const char* triggerList, Int_t runNumber)
{
  // Get the expected (from OCDB scalers) trigger count
  
  AliAnalysisTriggerScalers ts(runNumber,fgOCDBPath.Data());
  
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  TObjString* trigger;
  TIter next(triggers);
  ULong64_t n(0);
  
  while ( ( trigger = static_cast<TObjString*>(next()) ) )
  {
    AliAnalysisTriggerScalerItem* item = ts.GetTriggerScaler(runNumber,"L2A",trigger->String().Data());
    if (item)
    {
      n += item->Value();
    }
    delete item;
  }
  delete triggers;
  
  return n;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::GetSpectra(const char* what, const char* flavour) const
{
  /// get a given spectra
  
  TString swhat(what);
  TString sflavour(flavour);
  swhat.ToUpper();
  sflavour.ToUpper();
  
  TString spectraName(Form("/PSALL/%s/PP/%s/PSI-%s",
                           First(fgDefaultDimuonTriggers).Data(),
                           First(fgDefaultPairSelectionList).Data(),
                           swhat.Data()));

  if (sflavour.Length()>0)
  {
    spectraName += "-";
    spectraName += sflavour.Data();
  }

  if (IsSimulation())
  {
    spectraName.ReplaceAll("PSALL",fgDefaultEventSelectionForSimulations.Data());
    spectraName.ReplaceAll(First(fgDefaultDimuonTriggers).Data(),fgDefaultDimuonTriggerForSimulations.Data());
  }

  return SPECTRA(spectraName.Data());
}

//_____________________________________________________________________________
UInt_t AliAnalysisMuMu::GetSum(AliCounterCollection& cc, const char* triggerList,
                               const char* eventSelection, Int_t runNumber)
{
  TObjArray* ktrigger = cc.GetKeyWords("trigger").Tokenize(",");
  TObjArray* kevent = cc.GetKeyWords("event").Tokenize(",");
  TObjArray* a = TString(triggerList).Tokenize(" ");
  TIter next(a);
  TObjString* str;
  
  UInt_t n(0);
  
  TString sEventSelection(eventSelection);
  sEventSelection.ToUpper();
  
  if ( kevent->FindObject(sEventSelection.Data()) ) 
  {
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      if ( ktrigger->FindObject(str->String().Data()) )
      {
        if ( runNumber < 0 ) 
        {
          n +=  static_cast<UInt_t>(cc.GetSum(Form("trigger:%s/event:%s",str->String().Data(),eventSelection)));              
        }
        else
        {
          n +=  static_cast<UInt_t>(cc.GetSum(Form("trigger:%s/event:%s/run:%d",str->String().Data(),eventSelection,runNumber)));                        
        }
      }
    }
  }
  
  delete a;
  delete ktrigger;
  delete kevent;
  return n;
}

//_____________________________________________________________________________
Bool_t
AliAnalysisMuMu::GetCollections(const char* rootfile,
                                AliMergeableCollection*& mc,
                                AliCounterCollection*& cc,
                                AliAnalysisMuMuBinning*& bin,
                                std::set<int>& runnumbers)
{
  mc = 0x0;
  cc = 0x0;
  bin = 0x0;
  
  TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(rootfile));
  if (!f)
  {
    f = TFile::Open(rootfile);
  }
  
  if ( !f || f->IsZombie() )
  {
    return kFALSE;
  }

  f->GetObject("MC",mc);
  f->GetObject("CC",cc);
  
  TIter next(f->GetListOfKeys());
  TKey* key;
  
  while ( ( key = static_cast<TKey*>(next())) && !bin )
  {
    if ( strcmp(key->GetClassName(),"AliAnalysisMuMuBinning")==0 )
    {
      bin = dynamic_cast<AliAnalysisMuMuBinning*>(key->ReadObj());
    }
  }
  
  if (!mc || !cc)
  {
    AliErrorClass("Old file. Please upgrade it!");
    
    return kFALSE;
  }
  
  // get run list
  TObjArray* runs = cc->GetKeyWords("run").Tokenize(",");
  runs->Sort();
  TIter nextRun(runs);
  TObjString* srun;

  runnumbers.clear();
  
  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    runnumbers.insert(srun->String().Atoi());
  }
  
  delete runs;
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::IsSimulation() const
{
  // whether or not we have MC information
  
  return ( fMergeableCollection->Histo(Form("/INPUT/%s/MinvUS",fgDefaultEventSelectionForSimulations.Data())) != 0x0 );
}

//_____________________________________________________________________________
Int_t
AliAnalysisMuMu::Jpsi(const char* what, const char* binningFlavour)
{
  // Fit the J/psi (and psiprime) peaks for the triggers in fDimuonTriggers list
  // what="integrated" => fit only fully integrated MinvUS
  // what="pt" => fit MinvUS in pt bins
  // what="y" => fit MinvUS in y bins
  // what="pt,y" => fit MinvUS in (pt,y) bins
  
  TStopwatch timer;
  
  if (!fMergeableCollection)
  {
    AliError("No mergeable collection. Consider Upgrade()");
    return 0;
  }
  
  Int_t nfits(0);
  
  TObjArray* triggerArray = fDimuonTriggers.Tokenize(",");
  TObjArray* eventTypeArray = fEventSelectionList.Tokenize(",");
  TObjArray* pairCutArray = fPairSelectionList.Tokenize(",");
  TObjArray* whatArray = TString(what).Tokenize(",");
  
  TIter nextTrigger(triggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextWhat(whatArray);
  
  TObjString* trigger;
  TObjString* eventType;
  TObjString* pairCut;
  TObjString* swhat;
  
  while ( ( swhat = static_cast<TObjString*>(nextWhat()) ) )
  {    
    AliAnalysisMuMuBinning* binning(0x0);
    
    if ( fBinning && swhat->String().Length() > 0 )
    {
      binning = fBinning->Project("psi",swhat->String().Data(),binningFlavour);
    }
    else
    {
      binning = new AliAnalysisMuMuBinning;
      binning->AddBin("psi",swhat->String().Data());
    }
    
    std::cout << "++++++++++++ swhat=" << swhat->String().Data() << std::endl;
    
    if (!binning)
    {
      AliError("oups. binning is NULL");
      continue;
    }
    binning->Print();
    
    nextTrigger.Reset();
    
    while ( ( trigger = static_cast<TObjString*>(nextTrigger())) )
    {
      AliDebug(1,Form("TRIGGER %s",trigger->String().Data()));
      
      nextEventType.Reset();
      
      while ( ( eventType = static_cast<TObjString*>(nextEventType())) )
      {
        AliDebug(1,Form("--EVENTTYPE %s",eventType->String().Data()));
        
        nextPairCut.Reset();
        
        while ( ( pairCut = static_cast<TObjString*>(nextPairCut())) )
        {
          AliDebug(1,Form("----PAIRCUT %s",pairCut->String().Data()));
          
          AliDebug(1,"----Fitting...");
          
          AliAnalysisMuMuSpectra* spectra = FitParticle("psi",
                                                  trigger->String().Data(),
                                                  eventType->String().Data(),
                                                  pairCut->String().Data(),
                                                  "PP",
                                                  *binning);
          
          AliDebug(1,Form("----fitting done spectra = %p",spectra));
          
          if ( spectra )
          {
            ++nfits;
            
            TString id(Form("/%s/%s/PP/%s",eventType->String().Data(),
                            trigger->String().Data(),
                            pairCut->String().Data()));
            
            TObject* o = fMergeableCollection->GetObject(id.Data(),spectra->GetName());
          
            AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));
            
            if (o)
            {
              AliWarning(Form("Replacing %s/%s",id.Data(),spectra->GetName()));
              fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectra->GetName()));
            }
          
            fMergeableCollection->Adopt(id.Data(),spectra);
            
            spectra->Print();
          }
        }
      }
    }
  }
  
  delete whatArray;
  delete triggerArray;
  delete eventTypeArray;
  delete pairCutArray;

  timer.Print();

  if (nfits)
  {
    Update();
//    ReOpen(fFilename,"UPDATE");
//    fMergeableCollection->Write("MC",TObjArray::kOverwrite);// | TObjArray::kSingleKey);
//    ReOpen(fFilename,"READ");
  }
  
  
  return nfits;
  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::PlotBackgroundEvolution(const char* gfile, const char* triggerList, Double_t ymax, Bool_t fillBoundaries)
{
  // plot the graphs found in the file (and generated using the ComputeBackgroundEvolution() method)
  
  TFile* f = TFile::Open(ExpandPathName(gfile).Data());
  
  if ( !f || !f->IsOpen() )
  {
    return;
  }
  
  SetColorScheme();
  
  
  TCanvas* c = new TCanvas("background-evolution","background-evolution");
  
  c->Draw();
  
  TLegend* l = new TLegend(0.4,0.6,0.97,0.97);
  l->SetFillColor(0);
  l->SetTextColor(AliAnalysisMuMu::kBlue);
  l->SetLineColor(AliAnalysisMuMu::kBlue);
  
  TObjArray* triggers = TString(triggerList).Tokenize(",");
  
  gStyle->SetOptTitle(0);
  
  TObjString* str(0x0);
  TIter next(triggers);
  Int_t i(0);
  Int_t run1(99999999);
  Int_t run2(0);
  
  std::set<int> runnumbers;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TGraph* g = static_cast<TGraph*>(f->Get(Form("mu_%s",str->String().Data())));
    if (!g) continue;
    for ( Int_t ir = 0; ir < g->GetN(); ++ir )
    {
      Int_t run = TMath::Nint(g->GetX()[ir]);
      runnumbers.insert(run);
      run1 = TMath::Min(run1,run);
      run2 = TMath::Max(run2,run);
    }
  }
  
  AliInfoClass(Form("run1 %d run2 %d",run1,run2));
  
  Double_t ymin(0);
  
  TH2* hframe = new TH2F("hframe","hframe",run2-run1+1,run1,run2,100,ymin,ymax);
  hframe->Draw();
  hframe->GetXaxis()->SetNoExponent();
  hframe->GetYaxis()->SetTitle("Background percentage");
  
  if (fillBoundaries)
  {
    AliAnalysisTriggerScalers ts(runnumbers,fgOCDBPath.Data());
    ts.DrawFills(ymin,ymax);
  }
  
  next.Reset();
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TGraph* g = static_cast<TGraph*>(f->Get(Form("mu_%s",str->String().Data())));
    if (!g)
    {
      AliErrorClass(Form("Graph mu_%s not found",str->String().Data()));
      continue;
    }
    
    Int_t color(i+1);
    
    if (i==0) color = AliAnalysisMuMu::kBlue;
    if (i==1) color = AliAnalysisMuMu::kOrange;
    
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20+i);
    
    g->Draw("LP");
    
    TLegendEntry* le = l->AddEntry(g,str->String().Data(),"lp");
    le->SetTextColor(color);
    
    g->GetYaxis()->SetTitleColor(AliAnalysisMuMu::kBlue);
    g->GetXaxis()->SetTitleColor(AliAnalysisMuMu::kBlue);
    //    g->Print();
    
    ++i;
  }
  
  hframe->Draw("sameaxis");
  
  l->Draw();
}

//_____________________________________________________________________________
void
AliAnalysisMuMu::PlotJpsiEvolution(const char* resultFile, const char* triggerList, Bool_t fillBoundaries,
                                   const char* efficiencyFile, Bool_t simulation)
{
  /// Will plot the Jpsi rate (or AccxEff if simulation=kTRUE) evolution.
  /// (JpsiRate is the number of Jpsi divided by the number of triggers)
  
  std::map<int, float> efficiencies;
  
  if ( efficiencyFile && strlen(efficiencyFile) > 0 )
  {
    std::ifstream in(gSystem->ExpandPathName(efficiencyFile));
    if (!in.bad())
    {
      char line[1024];
      int run;
      float eff;
      float dummy,errorl,errorh;
      int idummy;
      while ( in.getline(line,1023,'\n') )
      {
        sscanf(line,"%d, x[%d]=%f, y[%d]=%f, exl[%d]=%f, exh[%d]=%f, eyl[%d]=%f, eyh[%d]=%f",
               &run,&idummy,&dummy,&idummy,&eff,&idummy,&dummy,&idummy,&dummy,&idummy,&errorl,&idummy,&errorh);

        AliDebugClass(1,Form("%09d %8.6f +- %8.6f ",run,eff,errorl+errorh));

        efficiencies[run] = eff;
      }
    }
  }
  
  TFile* f = TFile::Open(gSystem->ExpandPathName(resultFile));
  
  std::set<int> runnumbers;
  
  TMap* m = static_cast<TMap*>(f->Get("rbr"));

  TIter next(m);
  TObjString* str;
  
  TObjArray files;
  files.SetOwner(kTRUE);
  
  while ( ( str = static_cast<TObjString*>(next())) )
  {
    files.Add(new TObjString(str->String()));
  }
  
  files.Sort();
  
  std::map<std::string, std::vector<float> > x_jpsirate;
  std::map<std::string, std::vector<float> > y_jpsirate;
  std::map<std::string, std::vector<float> > xerr_jpsirate;
  std::map<std::string, std::vector<float> > yerr_jpsirate;
  
  TIter nextTrigger(TString(triggerList).Tokenize(","));
  TObjString* trigger(0x0);
  
  int runMin(100000000);
  int runMax(0);

  TIter nextFile(&files);

  Double_t ymin(TMath::Limits<double>::Max());
  Double_t ymax(TMath::Limits<double>::Min());
  

  while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
  {
    Double_t sumw(0);
    Double_t n(0);
    
    TString triggerClass(trigger->String());
    
    nextFile.Reset();
        
    while ( ( str = static_cast<TObjString*>(nextFile())) )
    {
      TObjArray* a = static_cast<TObjArray*>(m->GetValue(str->String().Data()));
      if (!a) continue;
      AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(a->FindObject(triggerClass.Data()));
      if (!r) continue;

      TString period;
      int aodtrain,esdpass,runnumber;

      if ( DecodeFileName(str->String().Data(),period,esdpass,aodtrain,runnumber) )
      {
        runnumbers.insert(runnumber);
        
        runMin = TMath::Min(runMin,runnumber);
        runMax = TMath::Max(runMax,runnumber);
        
        x_jpsirate[triggerClass.Data()].push_back(runnumber);
        xerr_jpsirate[triggerClass.Data()].push_back(0.5);
        
        Double_t y(0.0);
        Double_t yerr(0.0);
        TString what("RateJpsi");
        if ( simulation )
        {
          what = "AccEffJpsi";
        }
        
        if ( TMath::Finite(r->GetValue("SigmaJpsi")) && r->NofTriggers() > 10 )
        {
          y = 100*r->GetValue(what.Data());
          yerr = 100*r->GetErrorStat(what.Data());
          
          if  (!efficiencies.empty() )
          {
            if (efficiencies.count(runnumber))
            {
              y /= ( efficiencies[runnumber] );
            }
            else
            {
              continue;
            }
          }
          
          ymin = TMath::Min(ymin,y);
          ymax = TMath::Max(ymax,y);
          
          sumw += y*r->NofTriggers();
          n += r->NofTriggers();
        }
        
        y_jpsirate[triggerClass.Data()].push_back(y);
        yerr_jpsirate[triggerClass.Data()].push_back(yerr);
      }
    }
    
    AliInfoClass(Form("Trigger %30s ponderated mean is %7.2f",trigger->String().Data(),sumw/n));
  }

  delete f;
  
  TString canvasName("cJpsiRateEvolution");
  
  if ( !efficiencies.empty() )
  {
    canvasName += "Corr";
    
  }
  TCanvas* c = new TCanvas(canvasName.Data(),canvasName.Data());
  
  c->Draw();
  
  Int_t nbins = runnumbers.size();
  Int_t xmin(runMin);
  Int_t xmax(runMax);
  
  if ( CompactGraphs() )
  {
    xmin = 0;
    xmax = nbins-1;
  }
  
  TH2* h = new TH2F("h",Form("h;RunNumber;%s",(simulation ? "AccxEff (%)":"J/#psi per CMUL (%)")),
                    nbins,xmin,xmax,100,ymin,ymax*1.2);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  if (!CompactGraphs())
  {
    h->GetXaxis()->SetNoExponent();
  }
  else
  {
    std::set<int>::const_iterator it;
    int i(0);
    
    for ( it = runnumbers.begin(); it != runnumbers.end(); ++it )
    {
      h->GetXaxis()->SetBinLabel(i,Form("%d",*it));
      ++i;
    }
    h->GetXaxis()->SetNdivisions(1,kFALSE);
    
  }
  
  h->Draw();

  if (fillBoundaries)
  {
    AliAnalysisTriggerScalers ts(runnumbers,fgOCDBPath);
    ts.DrawFills(ymin,ymax);
  }

  h->Draw("sameaxis");
  
  //c->RedrawAxis("g");

  nextTrigger.Reset();
  
  int i(0);
  int color[] = { 2,1,4,5,6 };
  int marker[] = { 20,23,25,21,22 };
  
  while ( ( trigger = static_cast<TObjString*>(nextTrigger())))
  {
    std::vector<float>& x = x_jpsirate[trigger->String().Data()];
    std::vector<float>& y = y_jpsirate[trigger->String().Data()];
    std::vector<float>& xerr = xerr_jpsirate[trigger->String().Data()];
    std::vector<float>& yerr = yerr_jpsirate[trigger->String().Data()];
    
    TGraphErrors* g = new TGraphErrors(x.size(),&x[0],&y[0],&xerr[0],&yerr[0]);
    
    g->SetLineColor(1);//color[i]);
    g->SetMarkerColor(color[i]);
    g->SetMarkerStyle(marker[i]);
    g->SetMarkerSize(0.7);
    g->GetXaxis()->SetNoExponent();
    
    if ( CompactGraphs() )
    {
      Compact(*g);
    }
    
    g->Draw("P");
    TString gname(trigger->String());
    gname.ReplaceAll("-","_");
    g->SetName(gname.Data());
//    g->Print();

    Double_t m2 = g->GetMean(2);
    
    TLine* line = new TLine(runMin,m2,runMax,m2);
    line->SetLineColor(color[i]);
    line->Draw();
    
    AliInfoClass(Form("TRIGGER %s MEAN %7.2f",trigger->String().Data(),m2));
    ++i;
  }
  
  
}

//_____________________________________________________________________________
TGraph* AliAnalysisMuMu::PlotEventSelectionEvolution(const char* trigger1, const char* event1,
                                                     const char* trigger2, const char* event2,
                                                     Bool_t drawFills,
                                                     Bool_t asRejection) const
{
  if (!CC()) return 0x0;
  
  const std::set<int>& runnumbers = RunNumbers();
  
  TGraphErrors* g = new TGraphErrors(runnumbers.size());
  
  std::set<int>::const_iterator it;
  Int_t i(0);

  Double_t ymin(TMath::Limits<double>::Max());
  Double_t ymax(TMath::Limits<double>::Min());

  for ( it = runnumbers.begin(); it != runnumbers.end(); ++it )
  {
    Int_t runNumber = *it;
    Double_t n = CC()->GetSum(Form("trigger:%s/event:%s/run:%d",trigger1,event1,runNumber));
    Double_t d = CC()->GetSum(Form("trigger:%s/event:%s/run:%d",trigger2,event2,runNumber));
    if (n>0 && d>0)
    {
      Double_t y = n/d;
      
      if ( fCorrectionPerRun )
      {
        Double_t xcorr,ycorr;
        fCorrectionPerRun->GetPoint(i,xcorr,ycorr); // note that the fact that xcorr==runNumber has been checked by the SetCorrectionPerRun method
        y *= ycorr;
        // FIXME: should get the correction error here
      }
      
      if ( asRejection ) y = 100*(1.0 - y);
      ymin = TMath::Min(ymin,y);
      ymax = TMath::Max(ymax,y);
      Double_t yerr = y*AliAnalysisMuMuResult::ErrorAB(n,TMath::Sqrt(n),d,TMath::Sqrt(d));
      g->SetPoint(i,runNumber,y);
      g->SetPointError(i,0.5,yerr);
      
      ++i;
    }
    
  }

  TH2* hframe = new TH2F(Form("%s %s-%s / %s-%s",(asRejection ? "1 - ":""),trigger1,event1,trigger2,event2),
                         Form("%s %s-%s / %s-%s",(asRejection ? "1 - ":""),trigger1,event1,trigger2,event2),
                         runnumbers.size()+50,
                         *(runnumbers.begin())-25,
                         *(runnumbers.rbegin())+25,100,0,ymax*1.3);
  
  gStyle->SetOptStat(0);
  
  hframe->Draw();
  
  hframe->GetXaxis()->SetNoExponent();
           
  hframe->GetYaxis()->SetTitle(asRejection ? "Rejection (%)" : "Ratio");
  
  g->Set(i);
  g->SetTitle(Form("%s %s-%s / %s-%s",(asRejection ? "1 - ":""),trigger1,event1,trigger2,event2));
  g->GetXaxis()->SetNoExponent();
  g->Draw("lp");

  AliAnalysisTriggerScalers ts(RunNumbers(),fgOCDBPath.Data());

  if ( drawFills )
  {
    ts.DrawFills(ymin,ymax);
    g->Draw("lp");
  }
  
  
  std::map<std::string, std::pair<int,int> > periods;
  
  ts.GetLHCPeriodBoundaries(periods);
  
  TLegend* legend = new TLegend(0.15,0.82,0.90,0.92);
  legend->SetFillColor(0);
  Int_t n(0);
  

  for ( std::map<std::string, std::pair<int,int> >::const_iterator pit = periods.begin(); pit != periods.end(); ++pit )
  {
    std::string period = pit->first;
    int run1 = (pit->second).first;
    int run2 = (pit->second).second;
    int nruns(0);
    for ( std::set<int>::const_iterator rit = RunNumbers().begin(); rit != RunNumbers().end(); ++ rit )
    {
      if ( (*rit) >= run1 && (*rit) <= run2 )
      {
        ++nruns;
      }
    }
    AliInfo(Form("Period %s runs %6d-%6d ; %d actual runs",period.c_str(),run1,run2,nruns));
    
    g->Fit("pol0","+Q","",run1,run2);
    TF1* func = static_cast<TF1*>(g->GetListOfFunctions()->Last());
    if (func)
    {
      func->SetLineColor(2+n);
      legend->AddEntry(func,Form("%s %5.2f #pm %5.2f %s (rel. error %5.2f %%)",period.c_str(),func->GetParameter(0),func->GetParError(0),
                                 (asRejection ? "%":""),100*func->GetParError(0)/func->GetParameter(0)));
      ++n;
    }
  }

  legend->SetNColumns(3);

  Double_t mean = TMath::Mean(g->GetN(),g->GetY());
  Double_t rms = TMath::RMS(g->GetN(),g->GetY());
  
  legend->AddEntry("",Form("Mean %5.2f RMS %5.2f (%5.2f %%)",mean,rms,(mean) ? 100.0*rms/mean : 0.0),"");
  
  legend->Draw();
  
  return g;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Print(Option_t* opt) const
{
    /// printout
  std::cout << "Reading from file : " << fFilename.Data() << std::endl;
  std::cout << "List of dimuon triggers to consider : " << DimuonTriggerList() << std::endl;
  std::cout << "List of   muon triggers to consider : " << MuonTriggerList() << std::endl;
  std::cout << "List of     MB triggers to consider : " << MinbiasTriggerList() << std::endl;
  std::cout << "Event selection list : " << EventSelectionList() << std::endl;
  std::cout << "Pair  selection list : " << PairSelectionList() << std::endl;
  
  std::cout << RunNumbers().size() << " runs";
  if ( fCorrectionPerRun )
  {
    std::cout << " with correction factors";
  }
  std::cout << std::endl;
  Int_t i(0);
  for ( std::set<int>::const_iterator it = RunNumbers().begin(); it != RunNumbers().end(); ++it )
  {
    std::cout << (*it);
    if ( fCorrectionPerRun )
    {
      std::cout << Form("(%e)",fCorrectionPerRun->GetY()[i]);
    }
    std::cout << ",";
    ++i;
  }
  std::cout << std::endl;
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( sopt.Contains("BIN") && BIN() )
  {
    std::cout << "Binning : " << std::endl;
    TString topt(sopt);
    topt.ReplaceAll("BIN","");
    BIN()->Print(topt.Data());
  }
  if ( sopt.Contains("MC") && MC() )
  {
    TString topt(sopt);
    topt.ReplaceAll("MC","");
    MC()->Print(topt.Data());
  }
  if ( sopt.Contains("CC") && CC() )
  {
    CC()->Print("trigger/event");
  }
  
  if ( sopt.Contains("SIZE") )
  {
    TFile* f = ReOpen(fFilename,"READ");
    TIter next(f->GetListOfKeys());
    TKey* key;
    
    while ( ( key = static_cast<TKey*>(next()) ) )
    {
      std::cout << key->GetName() << " " << key->GetNbytes() << " " << key->GetObjlen() << std::endl;
    }
  }
}

//_____________________________________________________________________________
TObjArray*
AliAnalysisMuMu::ReadFileList(const char* filelist)
{
  //
  // read the filelist and try to order it by runnumber
  //
  // filelist can either be a real filelist (i.e. a text file containing
  // root filenames) or a root file itself.
  //
  
  char line[1024];
  
  TObjArray* files = new TObjArray;
  files->SetOwner(kTRUE);
  
  TString sfilelist(ExpandPathName(filelist));
  
  if ( sfilelist.EndsWith(".root") )
  {
    files->Add(new TObjString(sfilelist.Data()));
    return files;
  }
  
  std::set<int> runnumbers;
  std::map<int,std::string> filemap;
  
  std::ifstream in(sfilelist.Data());
  
  TString period;
  int aodtrain,esdpass,runnumber;
  
  while ( in.getline(line,1022,'\n') )
  {
    DecodeFileName(line,period,esdpass,aodtrain,runnumber);
    
    AliDebugClass(1,Form("line %s => period %s esdpass %d aodtrain %d runnumber %09d",
                         line,period.Data(),esdpass,aodtrain,runnumber));
    
    filemap.insert(std::make_pair<int,std::string>(runnumber,line));
    runnumbers.insert(runnumber);
  }
  
  in.close();
  
  std::set<int>::const_iterator it;
  
  for ( it = runnumbers.begin(); it != runnumbers.end(); ++it )
  {
    files->Add(new TObjString(filemap[*it].c_str()));
  }
  
  return files;
}

//_____________________________________________________________________________
TFile* AliAnalysisMuMu::ReOpen(const char* filename, const char* mode) const
{
  /// Tries to reopen the file with a new mode
  
  TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(filename));
  
  if (f)
  {
    delete f;
  }
  
  f = TFile::Open(filename,mode);
  
  if ( !f || !f->IsOpen() )
  {
    AliError(Form("Cannot open file %s in mode %s",filename,mode));
    return 0x0;
  }
  
  return f;
}

//_____________________________________________________________________________
Int_t
AliAnalysisMuMu::RunNumberFromFileName(const char* filename)
{
  TString period;
  int esdpass,aodtrain,runnumber;
  Bool_t ok = DecodeFileName(filename,period,esdpass,aodtrain,runnumber);
  if ( ok ) return runnumber;
  return -1;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SetColorScheme()
{
  new TColor(AliAnalysisMuMu::kBlue,4/255.0,44/255.0,87/255.0,"my blue");
  new TColor(AliAnalysisMuMu::kOrange,255/255.0,83/255.0,8/255.0,"my orange");
  new TColor(AliAnalysisMuMu::kGreen,152/255.0,202/255.0,52/255.0,"my green");
  
  gStyle->SetGridColor(AliAnalysisMuMu::kBlue);
  
  gStyle->SetFrameLineColor(AliAnalysisMuMu::kBlue);
  gStyle->SetAxisColor(AliAnalysisMuMu::kBlue,"xyz");
  gStyle->SetLabelColor(AliAnalysisMuMu::kBlue,"xyz");
  
  gStyle->SetTitleColor(AliAnalysisMuMu::kBlue);
  gStyle->SetTitleTextColor(AliAnalysisMuMu::kBlue);
  gStyle->SetLabelColor(AliAnalysisMuMu::kBlue);
  gStyle->SetStatTextColor(AliAnalysisMuMu::kBlue);
  
  gStyle->SetOptStat(0);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::SetCorrectionPerRun(const TGraph& corr, const char* formula)
{
    /// Sets the graph used to correct values per run
  delete fCorrectionPerRun;
  fCorrectionPerRun=0x0;
  
  // check that corr has the same runs as we do
  
  Int_t i(0);
  
  for ( std::set<int>::const_iterator it = RunNumbers().begin(); it != RunNumbers().end(); ++it )
  {
    Int_t corrRun = TMath::Nint(corr.GetX()[i]);
    if (corrRun != *it)
    {
      AliError(Form("%d-th run mistmatch %d vs %d",i,corrRun,*it));
      
      return kFALSE;
    }
    ++i;
  }
  
  fCorrectionPerRun = new TGraphErrors(corr.GetN());

  TFormula* tformula(0x0);
  if ( strlen(formula) > 0 )
  {
    tformula = new TFormula("SetCorrectionPerRunFormula",formula);
  }

  i = 0;
  
  for ( std::set<int>::const_iterator it = RunNumbers().begin(); it != RunNumbers().end(); ++it )
  {
    Double_t y = corr.GetY()[i];
    
    if ( tformula )
    {
      y = tformula->Eval(y);
    }
    fCorrectionPerRun->SetPoint(i,corr.GetX()[i],y);
    ++i;
  }

  delete formula;
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SetNofInputParticles(AliAnalysisMuMuResult& r)
{
  /// Set the "NofInput" variable(s) of one result
  
  TString hname(Form("MinvUS%s",r.Bin().AsString().Data()));

  TH1* hinput = fMergeableCollection->Histo("/INPUT/INYRANGE",hname.Data());

  if (!hinput)
  {
    AliError(Form("Got a simulation file where I did not find histogram /INPUT/INYRANGE/%s",hname.Data()));

  }
  else
  {
    r.SetNofInputParticles(*hinput);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::SPECTRA(const char* fullpath) const
{
  /// Shortcut method to get to a spectra
  if (!MC()) return 0x0;
  
  return static_cast<AliAnalysisMuMuSpectra*>(MC()->GetObject(fullpath));
}

//_____________________________________________________________________________
void AliAnalysisMuMu::TriggerCountCoverage(const char* triggerList,
                                           Bool_t compact,
                                           Bool_t orderByTriggerCount)
{
  // Give the fraction of triggers (in triggerList) relative 
  // to what is expected in the scalers
  
  TGrid::Connect("alien://"); // to insure the "Trying to connect to server... message does not pollute our output later on...
  
  AliLog::EType_t oldLevel = static_cast<AliLog::EType_t>(AliLog::GetGlobalLogLevel());
  
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  
  if (!fMergeableCollection || !fCounterCollection) return;
  
  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);
  
  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);
  
  TObjString* srun;
  TObjString* strigger;
  
  TString striggerList(triggerList);
  
  ULong64_t total(0);
  ULong64_t totalExpected(0);
  TString msg;
  std::multimap<ULong64_t,std::string> messages;
  
  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    msg.Form("RUN %09d ",srun->String().Atoi());
    
    if (!compact)
    {
        msg += "\n";
    }

    ULong64_t nmax(0);

    nextTrigger.Reset();
    
    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      if ( !striggerList.Contains(strigger->String().Data()) ) 
      {
        continue;
      }
      ULong64_t n = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                                strigger->String().Data(),"ALL",srun->String().Atoi())));
   
      ULong64_t expected = GetTriggerScalerCount(strigger->String().Data(),srun->String().Atoi());
    
      
      nmax = TMath::Max(n,nmax);
      
      total += n;
      totalExpected += expected;
      
      msg += TString::Format("%30s %9lld expected %9lld [%s] ",strigger->String().Data(),n,expected,
                             (n>expected ? "!" : " "));
      
      if ( expected > 0 ) {
        msg += TString::Format("fraction %5.1f %%",n*100.0/expected);
      }

      if (!compact)
      {
        msg += "\n";
      }
    }
    if (nmax>0)
    {
      if (!orderByTriggerCount)
      {
        std::cout << msg.Data() << std::endl;
      }
      else
      {
        messages.insert(std::make_pair<ULong64_t,std::string>(nmax,msg.Data()));
      }
    }
  }
  
  std::multimap<ULong64_t,std::string>::const_reverse_iterator it;
  
  ULong64_t current(0);
  Int_t n(0);
  
  for ( it = messages.rbegin(); it != messages.rend(); ++it )
  {
    ++n;
    current += it->first;
    Double_t percent = ( total > 0.0 ? current*100.0/total : 0.0);
    std::cout << Form("%10lld",it->first) << " " << it->second << " percentage of total = " << Form("%7.2f %% %3d",percent,n ) << std::endl;
  }

  std::cout << Form("--- TOTAL %lld expected %lld fraction %5.1f %%",
                    total,totalExpected,totalExpected ? total*100.0/totalExpected : 0.0) << std::endl;
  

  AliLog::SetGlobalLogLevel(oldLevel);
  delete triggers;
  delete runs;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::UnsetCorrectionPerRun()
{
    // drop the correction factors
  delete fCorrectionPerRun;
  fCorrectionPerRun=0x0;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Update()
{
  /// update the current file with memory
 
  ReOpen(fFilename,"UPDATE");

  MC()->Write("MC",TObject::kSingleKey);

  ReOpen(fFilename,"READ");
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::Upgrade(const char* filename)
{
  /// Upgrade a file
  AliAnalysisMuMu m(filename);
  
  return m.Upgrade();
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::Upgrade()
{
  /// Upgrade the current file
  /// - from single list to one key per object, if needed
  /// - from histogramCollection to mergeableCollection, if needed

  TFile* f = ReOpen(fFilename,"UPDATE");
  
  TList* list = static_cast<TList*>(f->Get("chist"));
  
  if (list)
  {
    // really old file where everything was in a single list
  
    AliHistogramCollection* hc = static_cast<AliHistogramCollection*>(list->At(0));
    AliCounterCollection* cc = static_cast<AliCounterCollection*>(list->At(1));
    
    AliMergeableCollection* mc = hc->Convert();
    
    f->cd();
    
    mc->Write("MC",TObject::kSingleKey);
    cc->Write("CC",TObject::kSingleKey);
    
    f->Delete("chist;*");
    
    f->Write();
    
  }
  else
  {
    AliHistogramCollection* hc = static_cast<AliHistogramCollection*>(f->Get("HC"));

    if ( hc )
    {
      // old file with histogram collection instead of mergeable collection
      
      AliMergeableCollection* mc = hc->Convert();

      f->cd();

      mc->Write("MC",TObject::kSingleKey);

      f->Delete("HC;*");
      
      f->Write();
    }
  }

  delete f;
  
  return kTRUE;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::CorrectSpectra(const char* type, const char* flavour)
{
  /// Correct one spectra
  
  if (!SIM())
  {
    AliError("Cannot compute corrected yield without associated MC file !");
    return 0x0;
  }

  const char* accEffSubResultName="COUNTJPSI:1";
  
  AliAnalysisMuMuSpectra* realSpectra = GetSpectra(type,flavour);
  AliAnalysisMuMuSpectra* simSpectra = SIM()->GetSpectra(type,flavour);
  
  if ( !realSpectra )
  {
    AliError("could not get real spectra");
    return 0x0;
  }
  
  if ( !simSpectra)
  {
    AliError("could not get sim spectra");
    return 0x0;
  }
  
  realSpectra->Correct(*simSpectra,"Jpsi",accEffSubResultName);

  Update();
  
  return realSpectra;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::ComputeYield(const char* type, const char* flavour)
{
  if (!SIM())
  {
    AliError("Cannot compute corrected yield without associated MC file !");
    return 0x0;
  }
  
  const char* accEffSubResultName="COUNTJPSI:1";
  
  AliAnalysisMuMuSpectra* realSpectra = GetSpectra(type,flavour);
  
  if ( !realSpectra )
  {
    AliError("could not get real spectra");
    return 0x0;
  }
  
  if (!realSpectra->HasValue("CoffNofJpsi"))
  {
    if (!CorrectSpectra(type,flavour))
    {
      AliError("Could not get corrected spectra");
      return 0x0;
    }
  }
  
  AliAnalysisMuMuSpectra* simSpectra = SIM()->GetSpectra(type,flavour);
  
  if ( !simSpectra)
  {
    AliErrorClass("could not get sim spectra");
    return 0x0;
  }
  
  Double_t nofCMUL7 = CC()->GetSum(Form("trigger:CMUL7-B-NOPF-MUON/event:PSALL"));
  Double_t nofCINT7 = CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD/event:PSALL"));
  Double_t nofCINT7w0MUL = CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD&0MUL/event:PSALL"));
  
  AliAnalysisMuMuBinning* binning = realSpectra->Binning();
  TObjArray* bins = binning->CreateBinObjArray();
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* bin;
  Int_t i(0);
  AliAnalysisMuMuResult* r;
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    r = static_cast<AliAnalysisMuMuResult*>(realSpectra->Bins()->At(i));
   
    StdoutToAliDebug(1,std::cout << "bin=";r->Print(););
    
    AliAnalysisMuMuResult* rsim = static_cast<AliAnalysisMuMuResult*>(simSpectra->Bins()->At(i));
    
    Double_t mbeq = nofCINT7w0MUL / ( nofCINT7 * nofCMUL7);
    Double_t mbeqError = mbeq * AliAnalysisMuMuResult::ErrorABC( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
                                                                nofCINT7,TMath::Sqrt(nofCINT7),
                                                                nofCMUL7,TMath::Sqrt(nofCMUL7));
    
    r->Set("MBR",nofCINT7/nofCINT7w0MUL,(nofCINT7/nofCINT7w0MUL)*AliAnalysisMuMuResult::ErrorAB( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
                                                                                                nofCINT7,TMath::Sqrt(nofCINT7)));
    
    Double_t yield =  r->GetValue("CorrNofJpsi") * mbeq;
    
    Double_t yieldError = yield * AliAnalysisMuMuResult::ErrorAB( r->GetValue("CorrNofJpsi"), r->GetErrorStat("CorrNofJpsi"),
                                                                 mbeq,mbeqError);
    
    r->Set("YJpsi",yield,yieldError);
    
    r->Set("NofInputJpsi",rsim->GetValue("NofInputJpsi",accEffSubResultName),rsim->GetErrorStat("NofInputJpsi",accEffSubResultName));
    r->Set("AccEffJpsi",rsim->GetValue("AccEffJpsi",accEffSubResultName),rsim->GetErrorStat("AccEffJpsi",accEffSubResultName));
    
    ++i;
  }
  
  delete bins;
  
  Update();
  
  return realSpectra;
}

////_____________________________________________________________________________
//AliAnalysisMuMuSpectra* AliAnalysisMuMu::ComputeYield(const char* realFile, const char* simFile,
//                                                      const  char* type)
//{
//  const char* accEffSubResultName="COUNTJPSI-1";
//
//  AliAnalysisMuMu real(realFile);
//  AliAnalysisMuMu sim(simFile);
//  
//  
//  AliAnalysisMuMuSpectra* realSpectra = static_cast<AliAnalysisMuMuSpectra*>(real.MC()->GetObject(Form("/PSALL/CMUL7-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-%s",type)));
//  AliAnalysisMuMuSpectra* simSpectra = static_cast<AliAnalysisMuMuSpectra*>(sim.MC()->GetObject(Form("/ALL/CMULLO-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-%s",type)));
//  
//  if ( !realSpectra )
//  {
//    AliErrorClass("could not get real spectra");
//    return 0x0;
//  }
//  
//  if ( !simSpectra)
//  {
//    AliErrorClass("could not get sim spectra");
//    return 0x0;
//  }
//  
//  AliAnalysisMuMuSpectra* corrSpectra = static_cast<AliAnalysisMuMuSpectra*>(realSpectra->Clone());
//  corrSpectra->Correct(*simSpectra,"Jpsi",accEffSubResultName);
//  
//  Double_t nofCMUL7 = real.CC()->GetSum(Form("trigger:CMUL7-B-NOPF-MUON/event:PSALL"));
//  Double_t nofCINT7 = real.CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD/event:PSALL"));
//  Double_t nofCINT7w0MUL = real.CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD&0MUL/event:PSALL"));
//  
//  AliAnalysisMuMuBinning* binning = corrSpectra->Binning();
//  TObjArray* bins = binning->CreateBinObjArray();
//  TIter nextBin(bins);
//  AliAnalysisMuMuBinning::Range* bin;
//  Int_t i(0);
//  AliAnalysisMuMuResult* r;
//  
//  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
//  {
//    r = static_cast<AliAnalysisMuMuResult*>(corrSpectra->Bins()->At(i));
//
//    AliAnalysisMuMuResult* rsim = static_cast<AliAnalysisMuMuResult*>(simSpectra->Bins()->At(i));
//    
//    Double_t mbeq = nofCINT7w0MUL / ( nofCINT7 * nofCMUL7);
//    Double_t mbeqError = mbeq * AliAnalysisMuMuResult::ErrorABC( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
//                                                                nofCINT7,TMath::Sqrt(nofCINT7),
//                                                                nofCMUL7,TMath::Sqrt(nofCMUL7));
//    
//    r->Set("MBR",nofCINT7/nofCINT7w0MUL,(nofCINT7/nofCINT7w0MUL)*AliAnalysisMuMuResult::ErrorAB( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
//                                                                                                nofCINT7,TMath::Sqrt(nofCINT7)));
//    
//    Double_t yield =  r->GetValue("CorrNofJpsi") * mbeq;
//    
//    Double_t yieldError = yield * AliAnalysisMuMuResult::ErrorAB( r->GetValue("CorrNofJpsi"), r->GetErrorStat("CorrNofJpsi"),
//                                                                 mbeq,mbeqError);
//    
//    r->Set("YJpsi",yield,yieldError);
//        
//    r->Set("NofInputJpsi",rsim->GetValue("NofInputJpsi",accEffSubResultName),rsim->GetErrorStat("NofInputJpsi",accEffSubResultName));
//    r->Set("AccEffJpsi",rsim->GetValue("AccEffJpsi",accEffSubResultName),rsim->GetErrorStat("AccEffJpsi",accEffSubResultName));
//    
//    ++i;
//  }
//
//  delete bins;
//  
//  return corrSpectra;
//}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::RABy(const char* realFile, const char* simFile, const char* type,
                                               const char* direction)
{
  /// Compute the RAB...
  Double_t rapidityShift = 0.465;// 0.5*TMath::Log(208.0/82.0);
  const Double_t sqrts=5.023;
  const Double_t ymax=TMath::Log(sqrts*1000.0/3.096916);
  const Double_t tab = 0.093e-6; // nb^-1
  const Double_t tabError = 0.0035E-6; // nb^-1
  const char* accEffSubResultName="COUNTJPSI:1";
  
  TF1 ydist("ydist","[0]*TMath::Exp(-(x*x)/(2.0*0.39*0.39))",0.,0.5);
  ydist.SetParameter(0,1.);

  //Normalization to the values presented by Zaida and Rosana on January 11th 2013 https://indico.cern.ch/conferenceDisplay.py?confId=224985 slide 22
  // Normalization is done in the rapidity range 2.75<y<3.25 where Rosanas values is 230.8+212.1
  Double_t y1_norma= 2.75/ymax;
  Double_t y2_norma= 3.25/ymax;
  Double_t normalization = 0.25*(230.8+212.1)/ydist.Integral(y1_norma, y2_norma);
  ydist.SetParameter(0,normalization);
//  AliInfoClass(Form("ymax=%e normalization=%f",ymax,ydist.Integral(y1_norma, y2_norma)));
  
  AliAnalysisMuMu real(realFile);
  AliAnalysisMuMu sim(simFile);
  
  
  AliAnalysisMuMuSpectra* realSpectra = static_cast<AliAnalysisMuMuSpectra*>(real.MC()->GetObject(Form("/PSALL/CMUL7-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-%s",type)));
  AliAnalysisMuMuSpectra* simSpectra = static_cast<AliAnalysisMuMuSpectra*>(sim.MC()->GetObject(Form("/ALL/CMULLO-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-%s",type)));
  
  if ( !realSpectra )
  {
    AliErrorClass("could not get real spectra");
    return 0x0;
  }
  
  if ( !simSpectra)
  {
    AliErrorClass("could not get sim spectra");
    return 0x0;
  }
  
  AliAnalysisMuMuSpectra* corrSpectra = static_cast<AliAnalysisMuMuSpectra*>(realSpectra->Clone());
  corrSpectra->Correct(*simSpectra,"Jpsi",accEffSubResultName);
  
  Double_t nofCMUL7 = real.CC()->GetSum(Form("trigger:CMUL7-B-NOPF-MUON/event:PSALL"));
  Double_t nofCINT7 = real.CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD/event:PSALL"));
  Double_t nofCINT7w0MUL = real.CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD&0MUL/event:PSALL"));
  
  AliAnalysisMuMuBinning* binning = realSpectra->Binning();
  TObjArray* bins = binning->CreateBinObjArray();
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* bin;
  Int_t i(0);
  AliAnalysisMuMuResult* r;
  
  Int_t n = bins->GetLast();
  
  TObjArray finalBins(n+1);
  finalBins.SetOwner(kTRUE);
  
  TObjArray finalResults(n+1);
  finalResults.SetOwner(kFALSE);
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    Double_t ylowlab = bin->Xmin();
    Double_t yhighlab = bin->Xmax();

    Double_t ylowcms, yhighcms;
    Double_t ylownorm, yhighnorm;
    
    if ( bin->IsNullObject() )
    {
      ylowlab = -4;
      yhighlab = -2.5;
    }
    
    if ( strcmp(direction,"pPb")==0 )
    {
      ylowcms = TMath::Abs(yhighlab) -  rapidityShift;
      yhighcms = TMath::Abs(ylowlab) - rapidityShift;
      ylownorm = ylowcms/ymax;
      yhighnorm = yhighcms/ymax;
    }
    else
    {
      ylowcms = ylowlab - rapidityShift;
      yhighcms = yhighlab - rapidityShift;
      ylownorm = -yhighcms/ymax;
      yhighnorm = -ylowcms/ymax;
    }
    
    
    Double_t brsigmapp = ydist.Integral(ylownorm,yhighnorm);
    Double_t brsigmappError = 0.0; // FIXME
    
    AliInfoClass(Form("y range : LAB %f ; %f CMS %f ; %f -> ynorm : %f ; %f -> BR x sigmapp = %f",
                      ylowlab,yhighlab,ylowcms,yhighcms,ylownorm,yhighnorm,brsigmapp));
    
    r = static_cast<AliAnalysisMuMuResult*>(corrSpectra->Bins()->At(i)->Clone());

    AliAnalysisMuMuResult* rsim = static_cast<AliAnalysisMuMuResult*>(simSpectra->Bins()->At(i));
    
    Double_t mbeq = nofCINT7w0MUL / ( nofCINT7 * nofCMUL7);
    Double_t mbeqError = mbeq * AliAnalysisMuMuResult::ErrorABC( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
                                         nofCINT7,TMath::Sqrt(nofCINT7),
                                         nofCMUL7,TMath::Sqrt(nofCMUL7));
    
    r->Set("MBR",nofCINT7/nofCINT7w0MUL,(nofCINT7/nofCINT7w0MUL)*AliAnalysisMuMuResult::ErrorAB( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
                                                                      nofCINT7,TMath::Sqrt(nofCINT7)));
    
    Double_t yield =  r->GetValue("CorrNofJpsi") * mbeq;

    Double_t yieldError = yield * AliAnalysisMuMuResult::ErrorAB( r->GetValue("CorrNofJpsi"), r->GetErrorStat("CorrNofJpsi"),
                                          mbeq,mbeqError);
    
    r->Set(Form("Y%sJpsi",direction),yield,yieldError);

    Double_t raa = yield/(tab*brsigmapp);
    Double_t raaError = AliAnalysisMuMuResult::ErrorABC(yield,yieldError,
                                                        tab,tabError,
                                                        brsigmapp,brsigmappError);
    r->Set(Form("R%sJpsi",direction),raa,raaError);

    r->Set("NofInputJpsi",rsim->GetValue("NofInputJpsi",accEffSubResultName),rsim->GetErrorStat("NofInputJpsi",accEffSubResultName));
    r->Set("AccEffJpsi",rsim->GetValue("AccEffJpsi",accEffSubResultName),rsim->GetErrorStat("AccEffJpsi",accEffSubResultName));
    
    AliAnalysisMuMuBinning::Range* bincm = new AliAnalysisMuMuBinning::Range(bin->Particle(),bin->Type(),ylowcms,yhighcms);
    
    r->SetBin(*bincm);
        
    finalBins.Add(bincm);
    finalResults.Add(r);
    
    ++i;
  }
  
  delete bins;
  
  AliAnalysisMuMuSpectra* spectra = new AliAnalysisMuMuSpectra(type,direction);
  
  for ( i = 0; i <= n; ++i )
  {
    Int_t j(i);
    if ( strcmp(direction,"pPb")==0 )
    {
      j = n-i;
    }
    
    r = static_cast<AliAnalysisMuMuResult*>(finalResults.At(j));

    bin = static_cast<AliAnalysisMuMuBinning::Range*>(finalBins.At(j));
    
    spectra->AdoptResult(*bin,r);
  }
  

  delete corrSpectra;
  
  return spectra;
}

//_____________________________________________________________________________
TGraph* AliAnalysisMuMu::ResultEvolution(const char* runlist, const char* period, const char* what, Bool_t forceRecomputation)
{
  std::vector<int> runs;
  AliAnalysisTriggerScalers::ReadIntegers(runlist,runs);

  TGraphErrors* g = new TGraphErrors(runs.size());
  
  TString direction("Pbp");
  
  if (TString(period) == "LHC13b" ||
    TString(period) == "LHC13c" ||
      TString(period) == "LHC13d" ||
      TString(period) == "LHC13e"
      )
  {
    direction = "pPb";
  }
  
  AliInfoClass(Form("period %s direction %s",period,direction.Data()));
  
  Double_t weightedMean(0.0);
  Double_t sumOfWeights(0.0);

  Double_t mean(0.0);
  TString subResultName("");
  TString swhat(what);
    
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i )
  {
    Int_t runNumber = runs[i];
    
    AliInfoClass(Form("RUN %09d",runNumber));
    
    TString realFile(Form("RUNBYRUN/%s_muon_calo_AODMUON000_%09d.saf.root",period,runNumber));
    
    TString simFileName(Form("RUNBYRUN/SIM_JPSI_%s_pp503_newalign_%09d.saf.root",period,runNumber));
    if ( direction == "pPb" )
    {
      simFileName = Form("RUNBYRUN/SIM_JPSI_%s_pp503_%09d.saf.root",period,runNumber);
    }

    TString simFile(simFileName);

    TString resultName(Form("%s%sJpsi",what,direction.Data()));
    
    if ( swhat == "MBR")
    {
      resultName = "MBR";
    }
    else if ( swhat.Contains("Acc") )
    {
      resultName.ReplaceAll(direction,"");
    }
    
    AliAnalysisMuMu mreal(realFile);
    
    AliAnalysisMuMuSpectra* real = static_cast<AliAnalysisMuMuSpectra*>(mreal.MC()->GetObject("/PSALL/CMUL7-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-"));

    if (!real || forceRecomputation)
    {
      mreal.Jpsi();
      real = static_cast<AliAnalysisMuMuSpectra*>(mreal.MC()->GetObject("/PSALL/CMUL7-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-"));
      if (!real)
      {
        AliErrorClass(Form("Could not get real spectra for run %d",runNumber));
        return 0x0;
      }
    }
    
    AliAnalysisMuMu msim(simFile);
    
    AliAnalysisMuMuSpectra* sim = static_cast<AliAnalysisMuMuSpectra*>(msim.MC()->GetObject("/ALL/CMULLO-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-"));

    if (!sim || forceRecomputation)
    {
      msim.SetEventSelectionList("ALL");
      msim.Jpsi();
      sim = static_cast<AliAnalysisMuMuSpectra*>(msim.MC()->GetObject("/ALL/CMULLO-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-"));
      if (!sim)
      {
        AliErrorClass(Form("Could not get sim spectra for run %d",runNumber));
        return 0x0;
      }
    }
    
    AliAnalysisMuMuSpectra* corrected = AliAnalysisMuMu::RABy(realFile.Data(),simFile.Data(),"",direction.Data());

    AliAnalysisMuMuResult* result = static_cast<AliAnalysisMuMuResult*>(corrected->Bins()->First());

    result->Print();

    Double_t value = result->GetValue(resultName.Data());
    Double_t error = result->GetErrorStat(resultName.Data());
    
    g->SetPoint(i,runNumber,value);
    g->SetPointError(i,1,error);
    
    Double_t n = mreal.CC()->GetSum(Form("trigger:CMUL7-B-NOPF-MUON/event:PSALL/run:%d",runNumber));
    
    weightedMean += n*result->GetValue(resultName.Data());
    sumOfWeights += n;
    
    mean += result->GetValue(resultName.Data());
    
//    std::cout << result->SubResults() << " " << result->GetError(resultName.Data()) << std::endl;

  }
  
  gStyle->SetOptFit(1);
  g->Draw("alp");
  g->Fit("pol0");
  g->SetTitle("");
  g->GetXaxis()->SetTitle("Run number");
  g->GetXaxis()->SetNoExponent();
  if ( TString(what) ==  "Y" )
  {
    g->GetYaxis()->SetTitle("J/#psi yield");
  }
  else if ( TString(what) == "R" )
  {
    g->GetYaxis()->SetTitle(Form("R_{%s}^{J/#psi}",direction.Data()));
  }
  else if ( TString(what).Contains("Acc") )
  {
    g->GetYaxis()->SetTitle("Acc#timesEff_{J/#psi}");    
  }
  else if ( TString(what).Contains("MBR") )
  {
    g->GetYaxis()->SetTitle("CINT7 / CINT7&0MUL");
  }
  
  if (CompactGraphs())
  {
    Compact(*g);
  }
  
  mean /= runs.size();
  weightedMean /= sumOfWeights;
  
  AliInfoClass(Form("Mean %e Weighted Mean %e",mean,weightedMean));
  
  return g;
}

//______________________________________________________________________________
void AliAnalysisMuMu::GetMBR(Int_t runNumber, const char* eventSelection, Double_t& value, Double_t& error) const
{
   // Get the scaling factor to go from CMUL to CINT7 for a given event selection
  value = 0.0;
  error = 0.0;
  if ( strlen(eventSelection) > 0 )
  {
    ULong64_t a = fCounterCollection->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD/event:%s/run:%d",
                                                  eventSelection,runNumber));
    
    ULong64_t b = fCounterCollection->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD&0MUL/event:%s/run:%d",
                                                  eventSelection,runNumber));
    
    value = b > 0 ? a/b : 0;
    error = value*AliAnalysisMuMuResult::ErrorAB(a,TMath::Sqrt(a),b,TMath::Sqrt(b));
  }
}

//______________________________________________________________________________
TGraph* AliAnalysisMuMu::MBREvolution(const char* eventSelection1, const char* eventSelection2) const
{
  if (!fCounterCollection) return 0x0;
  
  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  runs->Sort();
  TIter nextRun(runs);
  TObjString* srun;
  
  std::set<int> runnumbers;
  
  TGraphErrors* g = new TGraphErrors(runs->GetEntries());
  Int_t i(0);
  
  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    Int_t runNumber = srun->String().Atoi();
    
    runnumbers.insert(runNumber);
    
    Double_t mbr1,mbrError1;
    Double_t mbr2,mbrError2;
    
    GetMBR(runNumber,eventSelection1,mbr1,mbrError1);
    
    GetMBR(runNumber,eventSelection2,mbr2,mbrError2);

    Double_t mbr = mbr1;
    
    if ( mbr2 > 0 ) mbr /= mbr2;

    Double_t mbrError = mbr*AliAnalysisMuMuResult::ErrorAB(mbr1,mbrError1,mbr2,mbrError2);

    g->SetPoint(i,runNumber,mbr);
    g->SetPointError(i,0.5,mbrError);
    
    ++i;
  }
  
  g->GetXaxis()->SetNoExponent();
  
  AliAnalysisTriggerScalers ts(RunNumbers(),fgOCDBPath.Data());
  ts.DrawFills(0,1000);

  return g;
}




