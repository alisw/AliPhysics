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
#include "AliAnalysisMuMuBase.h" //Just to have avaliable the MCInputPrefix() method

#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuConfig.h"
#include "AliAnalysisMuMuFnorm.h"
#include "AliAnalysisMuMuGraphUtil.h"
#include "AliAnalysisMuMuJpsiResult.h"
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
#include "TProfile.h"
#include "THashList.h"
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
#include "TPRegexp.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cassert>
#include <map>
#include <set>
#include <string>
#include "TLatex.h"

using std::cout;
using std::endl;
using std::ifstream;

ClassImp(AliAnalysisMuMu)

//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu(const char* filename, const char* associatedSimFileName, const char* associatedSimFileName2, const char* beamYear) : TObject(),
fFilename(filename),
fCounterCollection(0x0),
fBinning(0x0),
fMergeableCollection(0x0),
fRunNumbers(),
fCorrectionPerRun(0x0),
fAssociatedSimulation(0x0),
fAssociatedSimulation2(0x0),
fParticleName(""),
fConfig(0x0)
{
  // ctor
  
  GetCollections(fFilename,fMergeableCollection,fCounterCollection,fBinning,fRunNumbers);
  
  if ( IsSimulation() )
  {
    TString sFilename(filename);
    if ( sFilename.Contains("JPSI",TString::kIgnoreCase) ) SetParticleName("JPsi"); //FIXME: DO a method for this
    else if ( sFilename.Contains("PSIP",TString::kIgnoreCase) ) SetParticleName("PsiP");
    else AliError("Unknown Particle Name in simulation");
  }
  
  if ( fCounterCollection )
  {
    if ( strlen(associatedSimFileName) )
    {
      fAssociatedSimulation = new AliAnalysisMuMu(associatedSimFileName);
      
      TString sAssociatedSimFileName(associatedSimFileName);
      if ( sAssociatedSimFileName.Contains("JPSI",TString::kIgnoreCase) ) fAssociatedSimulation->SetParticleName("JPsi");
      else if ( sAssociatedSimFileName.Contains("PSIP",TString::kIgnoreCase) ) fAssociatedSimulation->SetParticleName("PsiP");
      else AliError("Unknown Particle Name in associated simulation");     
    }
    
    if ( strlen(associatedSimFileName2) )
    {
      fAssociatedSimulation2 = new AliAnalysisMuMu(associatedSimFileName2);
      
      TString sAssociatedSimFileName2(associatedSimFileName2);
      if ( sAssociatedSimFileName2.Contains("JPSI",TString::kIgnoreCase) ) fAssociatedSimulation2->SetParticleName("JPsi");
      else if ( sAssociatedSimFileName2.Contains("PSIP",TString::kIgnoreCase) ) fAssociatedSimulation2->SetParticleName("PsiP");
      else AliError("Unknown Particle Name in associated simulation 2");
    }
    
    fConfig = new AliAnalysisMuMuConfig(beamYear);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMu::~AliAnalysisMuMu()
{
  // dtor
  
  if ( fAssociatedSimulation )
  {
    fAssociatedSimulation->Update();
  }
  if ( fAssociatedSimulation2 )
  {
    fAssociatedSimulation2->Update();
  }
  
  Update();
  
  delete fCounterCollection;
  delete fBinning;
  delete fMergeableCollection;
  delete fCorrectionPerRun;
  delete fAssociatedSimulation;
  delete fAssociatedSimulation2;
  delete fConfig;
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
      
      if ( !Config()->GetList(AliAnalysisMuMuConfig::kMinbiasTriggerList,IsSimulation()).Contains(strigger->String().Data()) &&
           !Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,IsSimulation()).Contains(strigger->String().Data()) &&
           !Config()->GetList(AliAnalysisMuMuConfig::kMuonTriggerList,IsSimulation()).Contains(strigger->String().Data()) ) continue;
          
      ULong64_t n = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                    strigger->String().Data(),"ALL",srun->String().Atoi())));

      details += TString::Format("\n%50s %10lld",strigger->String().Data(),n);
      

      ULong64_t nps = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                      strigger->String().Data(),"PSALL",srun->String().Atoi())));

      if ( doPS )
      {
        details += TString::Format(" PS %5.1f %%",nps*100.0/n);
      }

      if (nps)
      {
        ++nofPS;
      }
      
      if ( Config()->GetList(AliAnalysisMuMuConfig::kMinbiasTriggerList,IsSimulation()).Contains(strigger->String()) )
      {
        nmb += n;
        if ( totalNmb) (*totalNmb) += n;
        localNmb += n;
      }
      else if ( Config()->GetList(AliAnalysisMuMuConfig::kMuonTriggerList,IsSimulation()).Contains(strigger->String()) )
      {
        nmsl += n;
        if ( totalNmsl) (*totalNmsl) += n;
        localNmsl += n;
      }
      else if ( Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,IsSimulation()).Contains(strigger->String()) )
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
void AliAnalysisMuMu::CleanAllSpectra()
{
  /// Delete all the spectra we may have

  OC()->RemoveByType("AliAnalysisMuMuSpectra");
  Update();
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
AliAnalysisMuMuConfig* AliAnalysisMuMu::Config()
{
  /// Return the configuration
  return fConfig;
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
  
  if (!OC() || !BIN()) return;
  
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
  
  TString spectraName(Form("/%s/%s/%s/%s/%s-%s",eventType,trigger,centrality,pairCut,particle,stype.Data()));
  
  if ( strlen(flavour))
  {
    spectraName += "-";
    spectraName += flavour;
  }
  
  AliAnalysisMuMuSpectra* spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraName.Data()));
  
  AliDebug(1,Form("spectraName=%s spectra=%p",spectraName.Data(),spectra));

  TObjArray* spectraBins(0x0);
  if ( spectra )
  {
    spectraBins = spectra->BinContentArray();
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
    
    AliAnalysisMuMuJpsiResult* spectraBin(0x0);
    
    if ( spectraBins )
    {
      AliAnalysisMuMuResult* sr = static_cast<AliAnalysisMuMuResult*>(spectraBins->At(ci));
      
      spectraBin = static_cast<AliAnalysisMuMuJpsiResult*>(sr->SubResult(subresultname));
      
      AliDebug(1,Form("spectraBin(%s)=%p",subresultname,spectraBin));
    }
    
    TH1* h = OC()->Histo(name.Data());
    
    if ( spectraBin )
    {
      h = spectraBin->Histo();
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

//  AliWarning("Reimplement me!");
  
  if (!fConfig)
  {
    AliError("No configuration available yet. Don't know what to draw");
    return;
  }
  
  DrawMinv(type,particle,
           First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,IsSimulation())),
           First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,IsSimulation())),
           First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,IsSimulation())),
           First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,IsSimulation())),
           subresultname,
           flavour);
//           First(Config()->kEventSelectionList()).Data(),
//           First(Config()->kPairSelectionList()).Data(),
//           First(Config()->kCentralitySelectionList()).Data(),
//           subresultname,
//           flavour);
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
void AliAnalysisMuMu::TwikiOutputFnorm(const char* series) const
{
  // make a twiki-compatible output of the Fnorm factor(s)
  TObjArray* what = TString(series).Tokenize(",");
  TObjString* s;
  TObjArray graphs;
  TIter next(what);

  std::cout << "| *Run* |";
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TGraph* g = static_cast<TGraph*>(OC()->GetObject(Form("/FNORM/GRAPHS/%s",s->String().Data())));
    if (!g)
    {
      AliError(Form("Could not find graph for %s",s->String().Data()));
      continue;
    }
    std::cout << " *" << s->String().Data();
    if ( s->String().BeginsWith("RelDif") ) std::cout << " %";
    std::cout << "*|";
    graphs.Add(g);
  }
  
  std::cout << std::endl;
  
  TGraphErrors* g0 = static_cast<TGraphErrors*>(graphs.First());
  if (!g0) return;
  
  for ( Int_t i = 0; i < g0->GetN(); ++i )
  {
    TString msg;
    
    msg.Form("|%6d|",TMath::Nint(g0->GetX()[i]));
    
    for ( Int_t j = 0; j < graphs.GetEntries(); ++j )
    {
      TGraphErrors* g = static_cast<TGraphErrors*>(graphs.At(j));
      
      msg += TString::Format(" %6.2f +- %6.2f |",g->GetY()[i],g->GetEY()[i]);
    }
    
    std::cout << msg.Data() << std::endl;
  }
  
  next.Reset();
  
  std::cout << "|*Weigthed mean (*)*|";

  AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/Fnorm"));
  
  if (!r)
  {
    AliError("Could not find Fnorm result !");
    return;
  }

  
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TString var("Fnorm");
    TString unit;
    
    if ( s->String().BeginsWith("Fnorm") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/Fnorm"));
    }
    else if ( s->String().BeginsWith("RelDif") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/RelDif"));
      unit = "%";
    }
      
    r->Exclude("*");
    r->Include(s->String().Data());

    std::cout << Form(" * %5.2f +- %5.2f %s * |",
                      r->GetValue(var.Data()),
                      r->GetErrorStat(var.Data()),
                      unit.Data());
  }
  
  next.Reset();
  
  std::cout << std::endl;

  std::cout << "|*RMS*|";

  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TString var("Fnorm");
    
    if ( s->String().BeginsWith("Fnorm") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/Fnorm"));
    }
    else if ( s->String().BeginsWith("RelDif") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/RelDif"));
    }
    
    r->Exclude("*");
    r->Include(s->String().Data());
    
    Double_t d = 100.0*r->GetRMS(var.Data())/r->GetValue(var.Data());
    
    std::cout << Form(" * %5.2f (%5.2f %%) * |",
                      r->GetRMS(var.Data()),d);
  }
  
  std::cout << std::endl;
  std::cout << "(*) weight is the number of CMUL7-B-NOPF-MUON triggers (physics-selected and pile-up corrected) in each run" << std::endl;
  
  delete what;
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
                             const AliAnalysisMuMuBinning& binning,
                             const char* spectraType,
                             Bool_t corrected)
{
  // Fit the minv/mpt spectra to find the given particle
  // Returns an array of AliAnalysisMuMuResult objects
  
  TProfile::Approximate(); //To avoid bins with error=0 due to low statstics
  
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
    AliError(Form("Did not find trigger %s",trigger));
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

  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  Int_t nruns = runs->GetEntries();
  delete runs;
  
  
  TString id(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut));
  
//  binning.Print();
  
  AliAnalysisMuMuSpectra* spectra(0x0);
  
  AliAnalysisMuMuBinning::Range* bin;
  TIter next(bins);
  
  TObjArray* fitTypeArray = Config()->GetList(AliAnalysisMuMuConfig::kFitTypeList,IsSimulation()).Tokenize(",");
  TIter nextFitType(fitTypeArray);
  TObjString* fitType;
  TString flavour;
  TString sSpectraType(spectraType);
  
  TString spectraName(binning.GetName());
  if ( flavour.Length() > 0 )
  {
    spectraName += "-";
    spectraName += flavour;
  }
  if ( corrected )
  {
    spectraName += "-";
    spectraName += "AccEffCorr";
  }
  
//  Int_t binN(0);
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    TString hname;
    if (!sSpectraType.CompareTo("minv"))
    {
      hname = corrected ? Form("MinvUS%s_AccEffCorr",bin->AsString().Data()) : Form("MinvUS%s",bin->AsString().Data());
    }
    else if (!sSpectraType.CompareTo("mpt"))
    {
      hname = corrected ? Form("MeanPtVsMinvUS%s_AccEffCorr",bin->AsString().Data()) : Form("MeanPtVsMinvUS%s",bin->AsString().Data());
    }
    else
    {
      AliError("Wrong spectra type choice: Posibilities are: 'minv' or 'mpt' ");
      return 0x0;
    }
    
    TString isCorr(corrected ? " AccEffCorr " : " ");
    std::cout << "---------------------------------//---------------------------------" << std::endl;
    std::cout << "Fitting" << isCorr.Data() << sSpectraType.Data() << " spectra in " << id.Data() << std::endl;
    
    TH1* histo = OC()->Histo(id.Data(),hname.Data());

    if (!histo)
    {
//      if (!fBinning && bin->IsIntegrated() )
//      {
//        // old file, we only had MinvUSPt
//        hminv = fMergeableCollection->Histo(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut),"MinvUSPt:py");
//      }
//      
//      if (!hminv)
//      {
        AliError(Form("Could not find histo %s",hname.Data()));
        continue;
//      }
    }
    
    histo = static_cast<TH1*>(histo->Clone(Form("%s%d",sSpectraType.Data(),n++)));
    
    const char* particleTmp = IsSimulation() ? GetParticleName() : "JPsi"; //At some point particleTmp should become particle (but for now particle is always = "psi")
  
    TString sparticleTmp(particleTmp);
    
    AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(particleTmp,
                                                                 *histo, // Result for current bin
                                                                 trigger,
                                                                 eventType,
                                                                 pairCut,
                                                                 centrality,
                                                                 *bin);
    
    r->SetNofTriggers(ntrigger);
    r->SetNofRuns(nruns);
    
    nextFitType.Reset();

    Int_t added(0);    
    
    while ( ( fitType = static_cast<TObjString*>(nextFitType())) )
    {
      // In this loop we create a Subresult for each fit inside the Result for current bin (AddFit will do)
      
      TString sFitType(fitType->String());
      
      if ( !sFitType.Contains(sSpectraType.Data()) ) continue;
      
      AliDebug(1,Form("<<<<<< fitType=%s bin=%s",sFitType.Data(),bin->Flavour().Data()));
      
      std::cout << "" << std::endl;
      std::cout << "---------------" << "Fit " << added + 1 << "------------------" << std::endl;
      std::cout << "Fitting " << hname.Data() << " with " << sFitType.Data() << std::endl;
      std::cout << "" << std::endl;
      
      if ( sFitType.Contains("mctails",TString::kIgnoreCase) ) //FIXME: Find a univoque way to determine the correctly the fit type
      {
        TString sbin = bin->AsString();
        TString spectraMCName = spectraName;
        AliAnalysisMuMuBinning::Range* binMC = bin;
        
        if ((sbin.Contains("MULT") || sbin.Contains("NCH") || sbin.Contains("DNCHDETA") || sbin.Contains("V0A") || sbin.Contains("V0ACENT") || sbin.Contains("NTRCORR")|| sbin.Contains("RELNTRCORR")) && !sbin.Contains("NTRCORRPT") && !sbin.Contains("NTRCORRY"))
        {
          //-------has to have a better way to do it
          AliAnalysisMuMuBinning* b = new AliAnalysisMuMuBinning;
          b->AddBin("psi","INTEGRATED");
          
          binMC = static_cast<AliAnalysisMuMuBinning::Range*>(b->CreateBinObjArray()->At(0));
          
          spectraMCName = b->GetName();
          delete b;
          
//          if ( flavour.Length() > 0 ) //Commented cause we set no flavour in "INTEGRATED" bin in the analysis task
//          {
//            spectraMCName += "-";
//            spectraMCName += flavour;
//          }
          if ( corrected )
          {
            spectraMCName += "-";
            spectraMCName += "AccEffCorr";
          }
          //-----------

        }
//        if( sbin.Contains("NTRCORRPT") || !sbin.Contains("NTRCORRY") )
//        {
//          spectraMCName.Remove(4,7);
//          
//          if ( spectraMCName.Contains("PT") )
//          {
//            AliAnalysisMuMuBinning* b = SIM()->BIN();
//            
//            TObjArray* binsPt = b->CreateBinObjArray(particle,"PT","");
//            
//            Int_t nEntrBinMC = binsPt->GetEntries();
//            
////            binsPt->Print();
//            
//            binMC = static_cast<AliAnalysisMuMuBinning::Range*>(binsPt->At(binN));
//            
//            if ( binN == nEntrBinMC - 1 ) binN = -1;
//            
//          }
//          
//        }
        
        //par = GetCB2Tails(*binInt,"MCTAILS",eventType,trigger,pairCut,corrected);Why I was taking the tails from the fitted data spectra?
        GetParametersFromMC(sFitType,Form("/%s/%s",centrality,pairCut),spectraMCName.Data(),binMC);
        
        if (sFitType.Length()>0)
        {
          added += ( r->AddFit(sFitType.Data()) == kTRUE );
        }
      }
      
      else if ( sFitType.Contains("mpt",TString::kIgnoreCase) && !sFitType.Contains("minv",TString::kIgnoreCase) )
      {
        std::cout << "++The Minv parameters will be taken from " << spectraName.Data() << std::endl;
        std::cout << "" << std::endl;
        
        AliAnalysisMuMuSpectra* minvSpectra = dynamic_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(id.Data(),spectraName.Data()));
       
        if ( !minvSpectra )
        {
          AliError(Form("Cannot fit mean pt: could not get the minv spectra for %s",id.Data()));
          continue;//return 0x0;
        }
        
        AliAnalysisMuMuJpsiResult* minvResult = static_cast<AliAnalysisMuMuJpsiResult*>(minvSpectra->GetResultForBin(*bin));
        
        if ( !minvResult )
        {
          AliError(Form("Cannot fit mean pt: could not get the minv result for bin %s in %s",bin->AsString().Data(),id.Data()));
          continue; //return 0x0;
        }
        
        TObjArray* minvSubResults = minvResult->SubResults();
        TIter nextSubResult(minvSubResults);
        AliAnalysisMuMuJpsiResult* fitMinv;
        TString subResultName;
        
        Int_t nSubFit(0);
        while ( ( fitMinv = static_cast<AliAnalysisMuMuJpsiResult*>(nextSubResult())) )
        {
          std::cout << "" << std::endl;
          std::cout <<  "      /-- SubFit " << nSubFit + 1 << " --/ " << std::endl;
          std::cout << "" << std::endl;
          
          TString fitMinvName(fitMinv->GetName());
          fitMinvName.Remove(fitMinvName.First("_"),fitMinvName.Sizeof()-fitMinvName.First("_"));
          
//          std::cout << fitMinvName.Data() << " ; " << fitMinv->GetName() << std::endl;
//          std::cout << "" << std::endl;
          
          if ( !sFitType.Contains(fitMinvName) ) continue; //FIXME: Ambiguous, i.e. NA60NEWPOL2EXP & NA60NEWPOL2 (now its ok cause only VWG and POL2EXP are used, but care)
          
          TString sMinvFitType(sFitType);
          
          GetParametersFromResult(sMinvFitType,fitMinv);//FIXME: Think about if this is necessary
          
          added += ( r->AddFit(sMinvFitType.Data()) == kTRUE );
          
          nSubFit++;
        }
      }
      
      else if ( sFitType.Contains("minv&mpt",TString::kIgnoreCase) )
      {
        AliWarning("Implement here the method to do the combined minv mpt fits");
        //FIXME: Shall we use the fitType or spectraType to choose to perform combined fits? Cause we have to check what kind of object is returned by the combined fit in order to decide if we put it in a different spectra(spectraType would be the flag,and here we should update the spectraName) or as a subresult(fitType in this case)
      }
      
      else
      {
        if ( sFitType.Contains("PSICB2",TString::kIgnoreCase) || sFitType.Contains("PSINA60NEW",TString::kIgnoreCase)) std::cout << "+Free tails fit... " << std::endl;
        else if ( sFitType.Contains("PSICOUNT",TString::kIgnoreCase) )  std::cout << Form("+Just counting %s...",GetParticleName()) << std::endl;
        else std::cout << "+Using predefined tails... " << std::endl; 
        
        if ( sFitType.Contains("minvJPsi") && !sparticleTmp.Contains("JPsi") )
        {
          std::cout << "This fitting funtion is set to fit JPsi: Skipping fit..." << std::endl;
          continue;
        }
        if ( sFitType.Contains("minvPsiP") && !sparticleTmp.Contains("PsiP") )
        {
          std::cout << "This fitting funtion is set to fit PsiP: Skipping fit..." << std::endl;
          continue;
        }
        
        added += ( r->AddFit(sFitType.Data()) == kTRUE );
      }
      
      std::cout << "-------------------------------------" << std::endl;
      std::cout << "" << std::endl;
    }
  
    if ( !added ) continue;
    
//    TObjArray* a = r->SubResults(); // TEST
//    a->Print();
    
    flavour = bin->Flavour();
    
    if (!spectra)
    {
      TString spectraSaveName = spectraName;
      if ( !sSpectraType.CompareTo("mpt") )
      {
        spectraSaveName += "-";
        spectraSaveName += "MeanPtVsMinvUS";
      }
      
      spectra = new AliAnalysisMuMuSpectra(spectraSaveName.Data());
    }
    
    Bool_t adoptOk = spectra->AdoptResult(*bin,r); // We adopt the Result for current bin into the spectra
    
    if ( adoptOk ) std::cout << "Result " << r->GetName() << " adopted in spectra " << spectra->GetName() << std::endl;
    else AliError(Form("Error adopting result %s in spectra %s",r->GetName(),spectra->GetName()));
      
    
    if ( IsSimulation() )
    {
      SetNofInputParticles(*r);
    }
  
//    binN++;
  } // loop on bins
  
  delete fitTypeArray;
  delete bins;
  
  return spectra;
}


////_____________________________________________________________________________
//AliAnalysisMuMuSpectra*
//AliAnalysisMuMu::FitMeanPt(const char* particle,
//                           const char* trigger,
//                           const char* eventType,
//                           const char* pairCut,
//                           const char* centrality,
//                           const AliAnalysisMuMuBinning& binning,const AliAnalysisMuMuSpectra& spectraMinv,Bool_t corrected)
//{
//  // Fit the dimuon mean pt spectra to find the particle mean pt
//  // Returns an array of AliAnalysisMuMuResult objects
//  
//  TObjArray* bins = binning.CreateBinObjArray(particle);
//  if (!bins)
//  {
//    AliError(Form("Did not get any bin for particle %s",particle));
//    return 0x0;
//  }
//  
//  //  AliAnalysisMuMuBinning* binningMinv(0x0);
//  //  TObjArray* binsMinv(0x0);
//  //  AliAnalysisMuMuBinning::Range* binInt(0x0);
//  //  AliAnalysisMuMuBinning::Range* testbin = static_cast<AliAnalysisMuMuBinning::Range*>(bins->At(0));
//  //  TString sbins = testbin->AsString();
//  
//  //  if ( sbins.Contains("MULT")) //For the multiplicity bins we will use the tails from the integrated Minv spectra
//  //  {
//  //    binningMinv = spectraMinv.Binning();
//  //    binsMinv = binningMinv->CreateBinObjArray(particle);
//  //    binInt = static_cast<AliAnalysisMuMuBinning::Range*>(binsMinv->At(0));
//  //  }
//  
//  
//  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
//  if ( !triggers->FindObject(trigger) )
//  {
//    AliDebug(1,Form("Did not find trigger %s",trigger));
//    delete bins;
//    delete triggers;
//    return 0x0;
//  }
//  delete triggers;
//  
//  TObjArray* events = fCounterCollection->GetKeyWords("event").Tokenize(",");
//  if ( !events->FindObject(eventType) )
//  {
//    AliError(Form("Did not find eventType %s",eventType));
//    delete bins;
//    delete events;
//    return 0x0;
//  }
//  delete events;
//  
//  Int_t ntrigger = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s",trigger,eventType)));
//  
//  if  (ntrigger<=0)
//  {
//    AliError(Form("No trigger for trigger:%s/event:%s",trigger,eventType));
//    delete bins;
//    return 0x0;
//  }
//  
//  AliAnalysisMuMuSpectra* sMinv = static_cast<AliAnalysisMuMuSpectra*>(spectraMinv.Clone());
//  if(!sMinv)
//  {
//    AliError("Did not find inv mass spectra");
//  }
//  
//  
//  //  binning.Print();
//  
//  AliAnalysisMuMuSpectra* spectra(0x0);
//  
//  AliAnalysisMuMuBinning::Range* bin;
//  //  AliAnalysisMuMuBinning::Range* binMinv;
//  TIter next(bins);
//  
//  //  TObjArray* fitTypeArray = fFitTypeList.Tokenize(",");
//  //  TIter nextFitType(fitTypeArray);
//  //  TObjString* fitType;
//  TString flavour;
//  
//  
//  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
//  {
//    std::cout << "---------------------------------" << std::endl;
//    std::cout << "Fitting mean Pt in " << bin->AsString().Data() << " " << "for" << " " << eventType << "/" << trigger << "/" << centrality << "/" << pairCut << std::endl;
//    std::cout << "---------------------------------" << std::endl;
//    
//    
//    TString pname;
//    if ( corrected ) pname = Form("MeanPtVsMinvUS%s_Corr",bin->AsString().Data());
//    else pname = Form("MeanPtVsMinvUS%s",bin->AsString().Data());
//    
//    TProfile* hmPt = static_cast<TProfile*>(fMergeableCollection->GetObject(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut),pname.Data()));
//    
//    if (!hmPt)
//    {
//      //      if (!fBinning && bin->IsNullObject() )
//      //      {
//      //        // old file, we only had MinvUSPt
//      //        hminv = fMergeableCollection->Histo(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut),"MinvUSPt:py");
//      //      }
//      
//      //      if (!hmPt)
//      //      {
//      AliDebug(1,Form("Could not find histo profile %s",pname.Data()));
//      continue;
//      //      }
//    }
//    hmPt->Approximate();
//    
//    
//    //    hmPt = static_cast<TH1*>(hmPt->Clone(Form("mPtv%d",n++))); //Reuse this
//    
//    //    if ( sbins.Contains("MULT") )
//    //    {
//    //      binMinv = binInt;
//    //    }
//    //    else
//    //    {
//    //      binMinv = bin;
//    //    }
//    
//    AliAnalysisMuMuJpsiResult* rMinv = static_cast<AliAnalysisMuMuJpsiResult*>(sMinv->GetResultForBin(*bin/*Minv*/));
//    
//    if( !rMinv )
//    {
//      AliError(Form("Did not find inv mass result inside spectra for bin %s",bin/*Minv*/->Quantity().Data()));
//    }
//    
//    TObjArray* fits = rMinv->SubResults();
//    AliAnalysisMuMuResult* fitMinv;
//    TIter nextFit(fits);
//    TString fitName;
//    
//    AliAnalysisMuMuJpsiResult* r = new AliAnalysisMuMuJpsiResult(*hmPt,
//                                                         trigger,
//                                                         eventType,
//                                                         pairCut,
//                                                         centrality,
//                                                         *bin);
//    
//    r->SetNofTriggers(ntrigger);
//    
//    nextFit.Reset();
//    
//    //    new TCanvas;
//    Int_t i = 1;
//    
//    while ( ( fitMinv = static_cast<AliAnalysisMuMuJpsiResult*>(nextFit())) )
//    {
//      std::cout << "" << std::endl;
//      std::cout << "---------------" << "Fit " << i << "------------------" << std::endl;
//      i++;
//      
//      fitName = fitMinv->Alias();
//      std::cout << "" << std::endl;
//      std::cout << "Getting signal parameters from " << eventType << "/" << trigger << "/" << centrality << "/" << pairCut << "/" << spectraMinv.GetName()  << std::endl;
//      std::cout << "" << std::endl;
//      
//      if(fitName.BeginsWith("JPSI2CB2VWG") || fitName.BeginsWith("MCTAILS"))
//      {
//        Double_t par[14] = {fitMinv->GetValue("kVWG"),fitMinv->GetValue("mVWG"),fitMinv->GetValue("sVWG1"),
//          fitMinv->GetValue("sVWG2"),fitMinv->GetValue("kPsi"),fitMinv->GetValue("MeanJpsi"),fitMinv->GetValue("SigmaJpsi"),
//          fitMinv->GetValue("alPsi"),fitMinv->GetValue("nlPsi"),fitMinv->GetValue("auPsi"),fitMinv->GetValue("nuPsi"),fitMinv->GetValue("kPsi'"),fitMinv->GetValue("NofJpsi"),fitMinv->GetErrorStat("NofJpsi")}; //Create a method to get the parameters //The last 2 parameters are not used in the fit
//        
//        if(fitName.BeginsWith("JPSI2CB2VWG"))
//        {
//          std::cout << "PRE-SET TAILS" << std::endl;
//          std::cout << "" << std::endl;
//          r->AddMeanPtFit("JPSI2CB2VWG:2",&par[0]); // FIXME: :2 should be in the default fit types but for meanpt fit naming is not there, change it
//        }
//        else
//        {
//          std::cout << "MC TAILS" << std::endl;
//          std::cout << "" << std::endl;
//          r->AddMeanPtFit("MCTAILS-JPSI2CB2VWG:2",&par[0]);
//        }
//        
//        
//      }
//      
//      //Make a part for NA48 and the rest of fitting functions.
//      
//    }
//    
//    flavour = bin->Flavour();
//    
//    if (!spectra)
//    {
//      TString spectraName(Form("MeanPtVsMinvUS-%s",binning.GetName()));
//      if ( flavour.Length() > 0 )
//      {
//        spectraName += "-";
//        spectraName += flavour;
//      }
//      if ( corrected )
//      {
//        spectraName += "-";
//        spectraName += "Corr";
//      }
//      
//      spectra = new AliAnalysisMuMuSpectra(spectraName.Data());
//    }
//    
//    spectra->AdoptResult(*bin,r);
//    
//    //    if ( IsSimulation() )
//    //    {
//    //      SetNofInputParticles(*r);
//    //    }
//    
//    
//  } // loop on bins
//  
//  //  delete fitTypeArray;
//  delete sMinv;
//  delete bins;
//  
//  return spectra;
//}

//_____________________________________________________________________________
void
AliAnalysisMuMu::GetParametersFromMC(TString& fitType, const char* pathCentrPairCut, const char* spectraName, AliAnalysisMuMuBinning::Range* bin) const
{
  /// Get the parameters from the associated simulation. Is intended to be used in minv fits, where we need the tails from JPsi (and Psi')
  // We can choose between the 2 associated simulations (a JPsi and Psi' ones) by setting the sim variable to "sim1" (fAssociatedSimulation by default) or "sim2" (fAssociatedSimulation2)
  
  if ( !SIM() && !SIM2() )
  {
    AliError("Cannot get MC tails without associated simulation(s) file(s) !");
    fitType = "";
    return;
  }
  
  TString subResultName("");
  if ( fitType.Contains("NA60NEW",TString::kIgnoreCase) ) subResultName = "PSINA60NEW";//_1.5_4.2
  else if ( fitType.Contains("CB2",TString::kIgnoreCase) ) subResultName = "PSICB2";//_2.2_3.9
  else
  {
    AliError("I don't know from which MC subresult take the tails");
    return;
  }
    
  TObjArray* simArr = new TObjArray;
  if ( SIM() ) simArr->Add(SIM());
  if ( SIM2() && fitType.Contains("INDEPTAILS",TString::kIgnoreCase) ) simArr->Add(SIM2()); // If we have set the key to get the JPsi ans PsiP tails
  
  TIter nextSim(simArr);
  AliAnalysisMuMu* currentSIM;
  
  TString spath(pathCentrPairCut);
  
  spath.Prepend(Form("/%s",First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kTRUE)).Data()));//FIXME: Care with this when there is more than one selection in the list
  spath.Prepend(Form("/%s",First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kTRUE)).Data()));
  
  
  while ( (currentSIM = static_cast<AliAnalysisMuMu*>(nextSim())) )
  {
    AliAnalysisMuMuSpectra* minvMCSpectra = currentSIM->SPECTRA(Form("%s/%s",spath.Data(),spectraName));
    if (!minvMCSpectra)
    {
      AliError(Form("Could not find spectra %s/%s for associated simulation",spath.Data(),spectraName));
      currentSIM->OC()->Print("*:Ali*");
      fitType = "";
      continue;
    }
    
    AliAnalysisMuMuJpsiResult* minvMCResult = static_cast<AliAnalysisMuMuJpsiResult*>(minvMCSpectra->GetResultForBin(*bin));
    if ( !minvMCResult )
    {
      AliError(Form("Cannot get MC tails cause the minv result for bin %s in %s/%s was not found",bin->AsString().Data(),spath.Data(),spectraName));
      fitType = "";
      continue;
    }
    
    AliAnalysisMuMuJpsiResult* r = dynamic_cast<AliAnalysisMuMuJpsiResult*>(minvMCResult->SubResult(subResultName.Data())); //FIXME: Find an independet way of naming results
    if  ( r && subResultName.Contains("PSICB2") )
    {
      fitType += Form(":al%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("al%s",currentSIM->GetParticleName())));
      fitType += Form(":nl%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("nl%s",currentSIM->GetParticleName())));
      fitType += Form(":au%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("au%s",currentSIM->GetParticleName())));
      fitType += Form(":nu%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("nu%s",currentSIM->GetParticleName())));
      
      std::cout << " Using MC " << currentSIM->GetParticleName() << " CB2 tails... " << std::endl;
      std::cout << std::endl;
    }
    else if ( r && subResultName.Contains("PSINA60NEW") )
    {
      fitType += Form(":p1L%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p1L%s",currentSIM->GetParticleName())));
      fitType += Form(":p2L%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p2L%s",currentSIM->GetParticleName())));
      fitType += Form(":p3L%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p3L%s",currentSIM->GetParticleName())));
      fitType += Form(":p1R%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p1R%s",currentSIM->GetParticleName())));
      fitType += Form(":p2R%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p2R%s",currentSIM->GetParticleName())));
      fitType += Form(":p3R%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p3R%s",currentSIM->GetParticleName())));
      
      fitType += Form(":aL%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("aL%s",currentSIM->GetParticleName())));
      fitType += Form(":aR%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("aR%s",currentSIM->GetParticleName())));
      
      std::cout << " Using MC " << currentSIM->GetParticleName() << " NA60New tails... " << std::endl;
      std::cout << std::endl;
    }
    else
    {
      AliError(Form("Cannot get MC tails. MC Subresult %s not found",minvMCResult->GetName()));
      fitType = "";
    }
  }
  delete simArr;
}

//_____________________________________________________________________________
void
AliAnalysisMuMu::GetParametersFromResult(TString& fitType, AliAnalysisMuMuJpsiResult* minvResult) const
{
  // Gets the parameters from a result, is intended to be used for mean pt fits where we need the signal and backgroud parameters
  
  AliWarning("Re-implement me !!!"); //FIXME: The parameters to get will depend on the fit function and also in this way is not suitable for other particles (ie Upsilon)(Find a way to get the particle(s) name)
  
  TString msg("");
  if ( minvResult )
  {
    // We take the usual parameters (the ones from JPsi and the normalization of the Psi')
    fitType += Form(":kJPsi=%f",minvResult->GetValue("kJPsi")); //FIXME: Names are not correct
    fitType += Form(":mJPsi=%f",minvResult->GetValue("mJPsi"));
    fitType += Form(":sJPsi=%f",minvResult->GetValue("sJPsi"));
    
    fitType += Form(":NofJPsi=%f",minvResult->GetValue("NofJPsi"));
    fitType += Form(":ErrStatNofJPsi=%f",minvResult->GetErrorStat("NofJPsi"));
    
    fitType += Form(":kPsiP=%f",minvResult->GetValue("kPsiP"));
    
    TString minvName(minvResult->GetName());
    
    TString minvRangeParam = minvName;
    minvRangeParam.Remove(0,minvRangeParam.First("_") + 1);
    fitType += Form(":MinvRS=%s",minvRangeParam.Data());
    
    fitType += Form(":FSigmaPsiP=%f",minvResult->GetValue("FSigmaPsiP"));
    
    if ( fitType.Contains("CB2",TString::kIgnoreCase) )
    {
      fitType += Form(":alJPsi=%f",minvResult->GetValue("alJPsi"));
      fitType += Form(":nlJPsi=%f",minvResult->GetValue("nlJPsi"));
      fitType += Form(":auJPsi=%f",minvResult->GetValue("auJPsi"));
      fitType += Form(":nuJPsi=%f",minvResult->GetValue("nuJPsi"));
      
      msg += "JPsi CB2 signal parameters";
      //    fitType += Form(":NofPsiP=%f",minvResult->GetValue("NofPsiP"));
      //    fitType += Form(":ErrStatNofPsiP=%f",minvResult->GetErrorStat("NofPsiP"));
      
      if ( fitType.Contains("INDEPTAILS") )
      {
//        minvName = minvResult->GetName();
        if ( minvName.Contains("INDEPTAILS") )
        {
          // In case we use independent parameters tails for JPsi and Psi' we take also the Psi' ones
          fitType += Form(":alPsiP=%f",minvResult->GetValue("alPsiP"));
          fitType += Form(":nlPsiP=%f",minvResult->GetValue("nlPsiP"));
          fitType += Form(":auPsiP=%f",minvResult->GetValue("auPsiP"));
          fitType += Form(":nuPsiP=%f",minvResult->GetValue("nuPsiP"));
          fitType += Form(":mPsiP=%f",minvResult->GetValue("mPsiP"));
          fitType += Form(":sPsiP=%f",minvResult->GetValue("sPsiP"));
          
          msg += " + PsiP CB2 signal parameters";
        }
        else
        {
          AliError(Form("Cannot get PsiP tails from result. Result %s does not contain PsiP tails info => Fit will fail ",minvResult->GetName()));
          fitType = "";
          return;
        }
      }
    }
    else if ( fitType.Contains("NA60NEW",TString::kIgnoreCase) )
    {
      fitType += Form(":p1LJPsi=%f",minvResult->GetValue("p1LJPsi"));
      fitType += Form(":p2LJPsi=%f",minvResult->GetValue("p2LJPsi"));
      fitType += Form(":p3LJPsi=%f",minvResult->GetValue("p3LJPsi"));
      fitType += Form(":p1RJPsi=%f",minvResult->GetValue("p1RJPsi"));
      fitType += Form(":p2RJPsi=%f",minvResult->GetValue("p2RJPsi"));
      fitType += Form(":p3RJPsi=%f",minvResult->GetValue("p3RJPsi"));
      
      fitType += Form(":aLJPsi=%f",minvResult->GetValue("aLJPsi"));
      fitType += Form(":aRJPsi=%f",minvResult->GetValue("aRJPsi"));
      
      msg += "JPsi NA60New signal parameters";
      
      if ( fitType.Contains("INDEPTAILS") )
      {
//        TString minvName(minvResult->GetName());
        if ( minvName.Contains("INDEPTAILS") )
        {
          // In case we use independent parameters tails for JPsi and Psi' we take also the Psi' ones
          fitType += Form(":p1LPsiP=%f",minvResult->GetValue("p1LPsiP"));
          fitType += Form(":p2LPsiP=%f",minvResult->GetValue("p2LPsiP"));
          fitType += Form(":p3LPsiP=%f",minvResult->GetValue("p3LPsiP"));
          fitType += Form(":p1RPsiP=%f",minvResult->GetValue("p1RPsiP"));
          fitType += Form(":p2RPsiP=%f",minvResult->GetValue("p2RPsiP"));
          fitType += Form(":p3RPsiP=%f",minvResult->GetValue("p3RPsiP"));
          
          fitType += Form(":aLPsiP=%f",minvResult->GetValue("aLPsiP"));
          fitType += Form(":aRPsiP=%f",minvResult->GetValue("aRPsiP"));
          
          msg += " + PsiP NA60New signal parameters";
          
        }
        else
        {
          AliError(Form("Cannot get PsiP tails from result. Result %s does not contain PsiP tails info => Fit will fail ",minvResult->GetName()));
          fitType = "";
          return;
        }
      }
    }
    else
    {
      AliError(Form("Cannot get the parameters from %s",minvResult->GetName()));
      fitType = "";
      return;
    }
    // Now we take the background parameters
    if ( fitType.Contains("VWG_") || fitType.Contains("VWGINDEPTAILS") ) //FIXME: Check that cannot be misunderstood(like Exp x Pol2..). In fact it can be misunderstood since the meanpt function name has also the name of the function to fit the bkg (free parameters). Also add the rest of the BKG functions
    {
      fitType += Form(":kVWG=%f",minvResult->GetValue("kVWG"));
      fitType += Form(":mVWG=%f",minvResult->GetValue("mVWG"));
      fitType += Form(":sVWG1=%f",minvResult->GetValue("sVWG1"));
      fitType += Form(":sVWG2=%f",minvResult->GetValue("sVWG2"));
      
      msg += " + VWG Bkg parameters";
    }
    else if ( fitType.Contains("POL2EXP_") || fitType.Contains("POL2EXPINDEPTAILS") )
    {
      fitType += Form(":kPol2Exp=%f",minvResult->GetValue("kPol2Exp"));
      fitType += Form(":pol0=%f",minvResult->GetValue("pol0"));
      fitType += Form(":pol1=%f",minvResult->GetValue("pol1"));
      fitType += Form(":pol2=%f",minvResult->GetValue("pol2"));
      fitType += Form(":exp=%f",minvResult->GetValue("exp"));
      
      msg += " + Pol2xExp Bkg parameters";
    }
    else if ( fitType.Contains("POL4EXP_") || fitType.Contains("POL4EXPINDEPTAILS") )
    {
      fitType += Form(":kPol4Exp=%f",minvResult->GetValue("kPol4Exp"));
      fitType += Form(":pol0=%f",minvResult->GetValue("pol0"));
      fitType += Form(":pol1=%f",minvResult->GetValue("pol1"));
      fitType += Form(":pol2=%f",minvResult->GetValue("pol2"));
      fitType += Form(":pol3=%f",minvResult->GetValue("pol3"));
      fitType += Form(":pol4=%f",minvResult->GetValue("pol4"));
      fitType += Form(":exp=%f",minvResult->GetValue("exp"));
      
      msg += " + Pol4xExp Bkg parameters";
    }
    std::cout << "Using " << msg.Data() << " from " << minvResult->GetName() <<  " inv mass result" << std::endl;
    std::cout << "" << std::endl;
  }
  else
  {
    AliError(Form("Cannot get tails from result. Result %s not found",minvResult->GetName()));
    fitType = "";
    return;
  }
}

//_____________________________________________________________________________
ULong64_t AliAnalysisMuMu::GetTriggerScalerCount(const char* triggerList, Int_t runNumber)
{
  // Get the expected (from OCDB scalers) trigger count
  
  AliAnalysisTriggerScalers ts(runNumber,Config()->OCDBPath());
  
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
  
  TString spectraName(Form("/%s/%s/PP/%s/PSI-%s",
                           First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,IsSimulation())).Data(),
                           First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,IsSimulation())).Data(),
                           First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,IsSimulation())).Data(),
                           swhat.Data()));

  if (sflavour.Length()>0)
  {
    spectraName += "-";
    spectraName += sflavour.Data();
  }

  return SPECTRA(spectraName.Data());
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::PlotAccEfficiency(const char* whatever)
{
  //FIXME::Make it general
  if ( !IsSimulation() )
  {
    AliError("Could not get AccxEff histo: Not simulation file");
    return 0x0;
  }
  
  TString path(Form("/%s/%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kTRUE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kTRUE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kTRUE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,kTRUE)).Data()));
  
  AliAnalysisMuMuSpectra* s = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("%s/%s",path.Data(),whatever)));
  if ( !s )
  {
    AliError(Form("No AccEff spectra %s found in %s",whatever,path.Data()));
    return 0x0;
  }
  
  return s->Plot(Form("AccEff%s",GetParticleName()),"PSICOUNT",kFALSE);//_2.2_3.9
                                                                                         
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::PlotSystematicsTestsRelative(const char* quantity,const char* flavour,const char* value2Test)
{
  /// what,quantity and flavour defines de binning to test (output will be 1 plot per bin)
  /// value2test is the observable we want to test ( i.e. NJpsi(bin)/integratedNJpsi, <pt>(bin)/integrated<pt>... )
  /// --value2test == yield -> NJpsi(bin)/integratedNJpsi
  /// --value2test == mpt -> <pt>(bin)/integrated<pt>

  TString path(Form("%s/%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,kFALSE)).Data()));
  
  TString svalue2Test(value2Test);
  TString shName("");
  TString sObsInt("");
  TString sObsBin("");
  TString sflavour(flavour);
  TString squantity(quantity);
  squantity.ToUpper();
  
  if ( !svalue2Test.CompareTo("yield",TString::kIgnoreCase) )
  {
    sObsInt = "PSI-INTEGRATED-AccEffCorr";
    sObsBin = Form("PSI-%s-AccEffCorr",squantity.Data());
    svalue2Test = "NofJPsi";
    shName = "N^{J/#psi}_{bin}/N^{J/#psi}_{int}";
  }
  else if ( !svalue2Test.CompareTo("mpt",TString::kIgnoreCase) )
  {
    sObsInt = "PSI-INTEGRATED-AccEffCorr-MeanPtVsMinvUS";
    sObsBin = Form("PSI-%s-AccEffCorr-MeanPtVsMinvUS",squantity.Data());
    svalue2Test = "MeanPtJPsi";
    shName = "<p_{t}>^{bin}/<p_{t}>^{int}";
  }
  else
  {
    AliError("unrecognized value to test");
    return  0x0;
  }
  
  TString id(Form("/TESTSYST/%s",path.Data()));
  
  //--Get the integrated results
  AliAnalysisMuMuSpectra* sInt = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/%s/%s",path.Data(),sObsInt.Data())));
  if ( !sInt )
  {
    AliError(Form("No integrated spectra %s found in %s",sObsInt.Data(),path.Data()));
    return 0x0;
  }
  
  TObjArray* bin = BIN()->CreateBinObjArray("psi","integrated","");
  if ( !bin )
  {
    AliError("No integrated bin found");
    return 0x0;
  }
  AliAnalysisMuMuBinning::Range* r = static_cast<AliAnalysisMuMuBinning::Range*>(bin->At(0));
  
  AliAnalysisMuMuJpsiResult* resInt =  static_cast<AliAnalysisMuMuJpsiResult*>(sInt->GetResultForBin(*r));
  if ( !resInt )
  {
    AliError(Form("No integrated result found in spectra %s at %s",sObsInt.Data(),path.Data()));
    return 0x0;
  }
  TObjArray* sresIntArray = resInt->SubResults();
//  Int_t nsresInt = sresIntArray->GetEntriesFast();

  bin = BIN()->CreateBinObjArray("psi",squantity.Data(),sflavour.Data());//memory leak?
  if ( !bin )
  {
    AliError(Form("%s-%s-%s binning does not exist","psi",squantity.Data(),sflavour.Data()));
    return 0x0;
  }
  
  Int_t nbin = bin->GetEntries();
  Int_t nsres = sresIntArray->GetEntries();
  TObjArray* sResultNameArray= new TObjArray();
  std::vector<std::vector<double> > valuesArr;
  valuesArr.resize(nbin+1, std::vector<double>(nsres,0));
  std::vector<std::vector<double> > valuesErrorArr;
  valuesErrorArr.resize(nbin+1, std::vector<double>(nsres,0));
  
  TIter nextIntSubResult(sresIntArray);
  AliAnalysisMuMuResult* sresInt(0x0);
  Int_t nBkgParametrizations(0);
  Int_t j(0);
  while ( ( sresInt = static_cast<AliAnalysisMuMuResult*>(nextIntSubResult()) ) ) //Integrated SubResult loop
  {
    TString sresIntName(sresInt->GetName());
    
    valuesArr[0][j] = sresInt->GetValue(svalue2Test.Data());
    valuesErrorArr[0][j] = sresInt->GetErrorStat(svalue2Test.Data());
    sResultNameArray->Add(new TObjString(sresIntName.Data()));
        
    if ( sresIntName.Contains("CB2") )
    {
      nBkgParametrizations++;
    }
    j++;
  }
  
  //--Get the bin per bin results
  AliAnalysisMuMuSpectra* sBin = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/%s/%s",path.Data(),sObsBin.Data())));
  if ( !sBin )
  {
    AliError(Form("No integrated spectra %s found in %s",sObsBin.Data(),path.Data()));
    return 0x0;
    delete bin;
  }
  
  TIter nextBin(bin);
  AliAnalysisMuMuJpsiResult* sresBin(0x0);
  Int_t i(1);
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ) //Bin loop
  {
    sresBin =  static_cast<AliAnalysisMuMuJpsiResult*>(sBin->GetResultForBin(*r));
    if ( !sresBin )
    {
      AliError(Form("No result found in spectra %s at %s for bin %s",sObsInt.Data(),path.Data(),r->AsString().Data()));
      return 0x0;
      delete bin;
    }
    
    TObjArray* sresBinArray = sresBin->SubResults();
    
    TIter nextSubResult(sresBinArray);
    j = 0;
    while ( (sresBin = static_cast<AliAnalysisMuMuJpsiResult*>(nextSubResult())) ) // Subresults loop
    {
      valuesArr[i][j] = sresBin->GetValue(svalue2Test.Data());
      valuesErrorArr[i][j] = sresBin->GetErrorStat(svalue2Test.Data());
      
      j++;

    } //End Subresults loop
    i++;

  } //End bin loop
  
  TH1* hsyst = new TH1F(Form("%s_Systematics",value2Test),Form("%s Systematics results",value2Test),nbin,0,nbin);
  
  TString binName("");
  for ( Int_t b = 1 ; b < i ; b++ ) //Bin loop
  {
    binName = static_cast<AliAnalysisMuMuBinning::Range*>(bin->At(b-1))->AsString().Data();
    TString savePath(Form("%s/%s",id.Data(),binName.Data()));
    
    Int_t binUCTotal(1);
    
    TH1* hratiosBin = new TH1D(Form("SystTests_%s_Bkg_%s",sObsBin.Data(),binName.Data()),
                          Form("%s Systematics tests for %s",binName.Data(),shName.Data()),j*nBkgParametrizations,0,
                          j*nBkgParametrizations);
    
    for ( Int_t k = 0 ; k <= j ; k++ )
    {
      TString hName("");
      
      if ( k == 0 ) hName = "UC";
      else hName = sResultNameArray->At(k-1)->GetName();
     
      
      TH1* hratios = new TH1D(Form("SystTests_%s_%s",sObsBin.Data(),hName.Data()),
                              Form("%s Systematics tests for %s(%s)",binName.Data(),shName.Data(),hName.Data()),j,0,j);
      hratios->SetNdivisions(j+1);
      
      TH1* hratiosUC(0x0);
      if ( k != 0 )
      {
        hratiosUC = new TH1D(Form("SystTests_%s_%s_Bkg",sObsBin.Data(),hName.Data()),
                                         Form("%s Systematics tests for %s(%s bkg)",binName.Data(),shName.Data(),hName.Data()),nBkgParametrizations,0,nBkgParametrizations);
        hratiosUC->SetNdivisions(nBkgParametrizations+1);
      }
      
      TString signalName(hName.Data());
      Int_t sizeName = signalName.Sizeof();
      signalName.Remove(2,sizeName-3);
      Int_t binUC(1);
      
      for ( Int_t l = 0; l < j ; l++) //Subresults loop
      {
        TString binSignalName(sResultNameArray->At(l)->GetName());
        
        Double_t ratio,ratioError;
        if ( k == 0)
        {
          ratio = valuesArr[b][l] / valuesArr[0][l];
          ratioError = TMath::Sqrt( TMath::Power(valuesErrorArr[b][l] / valuesArr[0][l],2.) + TMath::Power(valuesArr[b][l]*valuesErrorArr[0][l] / TMath::Power(valuesArr[0][l],2.),2.) );
          
          hratios->GetXaxis()->SetBinLabel(l+1,sResultNameArray->At(l)->GetName());
        }
        else
        {
          ratio = valuesArr[b][l] / valuesArr[0][k-1];
          ratioError = TMath::Sqrt( TMath::Power(valuesErrorArr[b][l] / valuesArr[0][k-1],2.) + TMath::Power(valuesArr[b][l]*valuesErrorArr[0][k-1] / TMath::Power(valuesArr[0][k-1],2.),2.) );
          
          hratios->GetXaxis()->SetBinLabel(l+1,Form("%s/%s",sResultNameArray->At(l)->GetName(),sResultNameArray->At(k-1)->GetName()));
          
          if ( binSignalName.Contains(signalName.Data()) )
          {
            hratiosUC->GetXaxis()->SetBinLabel(binUC,Form("%s/%s",sResultNameArray->At(l)->GetName(),sResultNameArray->At(k-1)->GetName()));
            
            hratiosUC->SetBinContent(binUC,ratio);
            hratiosUC->SetBinError(binUC,ratioError);
            
            hratiosBin->GetXaxis()->SetBinLabel(binUCTotal,Form("%s/%s",sResultNameArray->At(l)->GetName(),sResultNameArray->At(k-1)->GetName()));
            
            hratiosBin->SetBinContent(binUCTotal,ratio);
            hratiosBin->SetBinError(binUCTotal,ratioError);
            
            binUC++;
            binUCTotal++;
          }
        }
        
        hratios->SetBinContent(l+1,ratio);
        hratios->SetBinError(l+1,ratioError);
      }
      
      
      
      
//      TH1* o = OC()->Histo(Form("%s",savePath.Data()),hratios->GetName());
//      
//      if (o)
//      {
//        AliWarning(Form("Replacing %s/%s",savePath.Data(),hratios->GetName()));
//        OC()->Remove(Form("%s/%s",savePath.Data(),hratios->GetName()));
//      }
//      
//      Bool_t adoptOK = OC()->Adopt(savePath.Data(),hratios);
//      
//      if ( adoptOK ) std::cout << "+++syst histo " << hratios->GetName() << " adopted" << std::endl;
//      else AliError(Form("Could not adopt syst histo %s",hratios->GetName()));
//      
//      if ( hratiosUC )
//      {
//        o = OC()->Histo(Form("%s",savePath.Data()),hratiosUC->GetName());
//        
//        if (o)
//        {
//          AliWarning(Form("Replacing %s/%s",savePath.Data(),hratiosUC->GetName()));
//          OC()->Remove(Form("%s/%s",savePath.Data(),hratiosUC->GetName()));
//        }
//        
//        adoptOK = OC()->Adopt(savePath.Data(),hratiosUC);
//        
//        if ( adoptOK ) std::cout << "+++syst histo " << hratiosUC->GetName() << " adopted" << std::endl;
//        else AliError(Form("Could not adopt syst histo %s",hratiosUC->GetName()));
//      }
    }
    
    
    //Syst computation for bin
    Double_t num(0.),deno(0.);
    for ( Int_t m = 1 ; m <= hratiosBin->GetNbinsX() ; m++ )
    {
      Double_t value = hratiosBin->GetBinContent(m);
      Double_t error = hratiosBin->GetBinError(m);
      
      num += value/TMath::Power(error,2.);
      deno += 1./TMath::Power(error,2.);
    }
    
    Double_t mean = num/deno;
    Double_t v1(0.),v2(0.),sum(0.);
    for ( Int_t l = 1 ; l <= hratiosBin->GetNbinsX() ; l++ )
    {
      Double_t value = hratiosBin->GetBinContent(l);
      Double_t error = hratiosBin->GetBinError(l);
      
      Double_t wi = 1./TMath::Power(error,2.);
      v1 += wi;
      v2 += wi*wi;
      Double_t diff = value - mean;
      sum += wi*diff*diff;
      
    }

    Double_t syst = TMath::Sqrt( (v1/(v1*v1-v2)) * sum);
    
    hsyst->GetXaxis()->SetBinLabel(b,binName.Data());
    hsyst->SetBinContent(b,(syst*100.)/mean);
    
    
    //___
    TF1* meanF = new TF1("mean","[0]",0,j*nBkgParametrizations);
    meanF->SetParameter(0,mean);
    
    TF1* meanFPS = new TF1("meanPS","[0]",0,j*nBkgParametrizations);
    meanFPS->SetParameter(0,mean+syst);
    meanFPS->SetLineStyle(2);
    
    TF1* meanFMS = new TF1("meanMS","[0]",0,j*nBkgParametrizations);
    meanFMS->SetParameter(0,mean-syst);
    meanFMS->SetLineStyle(2);

    hratiosBin->GetListOfFunctions()->Add(meanF);
    hratiosBin->GetListOfFunctions()->Add(meanFPS);
    hratiosBin->GetListOfFunctions()->Add(meanFMS);
    
    TH1* o = OC()->Histo(Form("%s",savePath.Data()),hratiosBin->GetName());
    
    if (o)
    {
      AliWarning(Form("Replacing %s/%s",savePath.Data(),hratiosBin->GetName()));
      OC()->Remove(Form("%s/%s",savePath.Data(),hratiosBin->GetName()));
    }
    
    Bool_t adoptOK = OC()->Adopt(savePath.Data(),hratiosBin);
    
    if ( adoptOK ) std::cout << "+++syst histo " << hratiosBin->GetName() << " adopted" << std::endl;
    else AliError(Form("Could not adopt syst histo %s",hratiosBin->GetName()));
  }
  
  
  TH1* o = OC()->Histo(Form("%s",id.Data()),hsyst->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hsyst->GetName()));
    OC()->Remove(Form("%s/%s",id.Data(),hsyst->GetName()));
  }
  
  Bool_t adoptOK = OC()->Adopt(id.Data(),hsyst);
  
  if ( adoptOK ) std::cout << "+++syst histo " << hsyst->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt syst histo %s",hsyst->GetName()));
  
  
  delete bin;
  delete sResultNameArray;
  
  return 0x0;
  
}


//_____________________________________________________________________________
TH1* AliAnalysisMuMu::PlotJpsiYield(const char* whatever)
{
  
  //FIXME::Make it general
  if ( IsSimulation() )
  {
    AliError("Cannot compute J/Psi yield: Is a simulation file");
    return 0x0;
  }
  
  TString path(Form("/%s/%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,kFALSE)).Data()));
  
  AliAnalysisMuMuSpectra* s = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("%s/%s",path.Data(),whatever)));
  if ( !s )
  {
    AliError(Form("No spectra %s found in %s",whatever,path.Data()));
    return 0x0;
  }
  
  std::cout << "Number of J/Psi:" << std::endl;
  TH1* hry = s->Plot("NofJPsi","PSIPSIPRIMECB2VWGINDEPTAILS",kFALSE); //Number of Jpsi
  std::cout << "" << std::endl;
  
  std::cout << "Equivalent number of MB events:" << std::endl;
  TH1* hN = ComputeEquNofMB();
  std::cout << "" << std::endl;
  
  TH1* hy = static_cast<TH1*>(hry->Clone("CorrJPsiYields"));
  Double_t bR = 0.0593; // BR(JPsi->mu+mu-)
  Double_t bRerror = 0.0006 ;
  
  for (Int_t i = 1 ; i <= hy->GetNbinsX() ; i++)
  {
    Double_t yield = hry->GetBinContent(i)/(hN->GetBinContent(i)*bR);
    Double_t yieldError = TMath::Sqrt(TMath::Power(hry->GetBinError(i)/(hN->GetBinContent(i)*bR),2.) +
                                      TMath::Power(hN->GetBinError(i)*bR/TMath::Power(hN->GetBinContent(i)*bR,2.),2.) +
                                      TMath::Power(hry->GetBinContent(i)*hN->GetBinContent(i)*bRerror/TMath::Power(hN->GetBinContent(i)*bR,2.),2.));
    
    std::cout << yield << " +- " << yieldError << std::endl;
    
    hy->SetBinContent(i,yield);
    hy->SetBinError(i,yieldError);
  }
  
  delete hry;
  delete hN;
  
  return hy;
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
                                AliMergeableCollection*& oc,
                                AliCounterCollection*& cc,
                                AliAnalysisMuMuBinning*& bin,
                                std::set<int>& runnumbers)
{
  oc = 0x0;
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

  f->GetObject("OC",oc);
  if (!oc)
  {
    f->GetObject("MC",oc);
  }
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
  
  if (!oc || !cc)
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
  
//  return kFALSE;
  
  if (!fMergeableCollection) return kFALSE;
  
  TList* list = fMergeableCollection->CreateListOfKeys(0);
  TIter next(list);
  TObjString* str;
  Bool_t ok(kFALSE);
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    if ( str->String().Contains(AliAnalysisMuMuBase::MCInputPrefix()) ) ok = kTRUE;
  }
  delete list;
  
  return ok;
}

//_____________________________________________________________________________
Int_t
AliAnalysisMuMu::Jpsi(const char* what, const char* binningFlavour, Bool_t fitmPt)
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
  
  TObjArray* triggerArray = Config()->GetListElements(AliAnalysisMuMuConfig::kDimuonTriggerList,IsSimulation());
  TObjArray* eventTypeArray = Config()->GetListElements(AliAnalysisMuMuConfig::kEventSelectionList,IsSimulation());
  TObjArray* pairCutArray = Config()->GetListElements(AliAnalysisMuMuConfig::kPairSelectionList,IsSimulation());
  TObjArray* whatArray = TString(what).Tokenize(",");
  TObjArray* centralityArray = Config()->GetListElements(AliAnalysisMuMuConfig::kCentralitySelectionList,IsSimulation());
  
  TIter nextTrigger(triggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextWhat(whatArray);
  TIter nextCentrality(centralityArray);
  
  TObjString* trigger;
  TObjString* eventType;
  TObjString* pairCut;
  TObjString* swhat;
  TObjString* centrality;
  
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
    
    StdoutToAliDebug(1,std::cout << "++++++++++++ swhat=" << swhat->String().Data() << std::endl;);
    
    std::cout << "" << std::endl;
    std::cout << "++++++++++++++++++" << "NEW BIN TYPE" << "+++++++++++++++++++" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "++++++++++++ swhat=" << swhat->String().Data() << "++++++++++++++++++++" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    
    if (!binning)
    {
      AliError("oups. binning is NULL");
      continue;
    }
    
    StdoutToAliDebug(1,binning->Print(););
    
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
          
          nextCentrality.Reset();
          
          while ( ( centrality = static_cast<TObjString*>(nextCentrality()) ) )
          {
            AliDebug(1,"----Fitting...");
            
            AliAnalysisMuMuSpectra* spectra = FitParticle("psi",
                                                          trigger->String().Data(), //Uncomment
                                                          eventType->String().Data(),
                                                          pairCut->String().Data(),
                                                          centrality->String().Data(),
                                                          *binning);
            
            AliDebug(1,Form("----fitting done spectra = %p",spectra));
            
            TString id(Form("/%s/%s/%s/%s",eventType->String().Data(),
                            trigger->String().Data(),
                            centrality->String().Data(),
                            pairCut->String().Data()));
            TObject* o;
            
            if ( spectra )
            {
              ++nfits;
              
              o = fMergeableCollection->GetObject(id.Data(),spectra->GetName());
              
              AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));
              
              if (o)
              {
                AliWarning(Form("Replacing %s/%s",id.Data(),spectra->GetName()));
                fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectra->GetName()));
              }
              
              Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),spectra);
              
              if ( adoptOK ) std::cout << "+++Spectra " << spectra->GetName() << " adopted" << std::endl;
              else AliError(Form("Could not adopt spectra %s",spectra->GetName()));
           
              StdoutToAliDebug(1,spectra->Print(););
            }
            else AliError("Error creating spectra"); //Uncomment
            
            AliDebug(1,"----Fitting corrected spectra...");
            
            AliAnalysisMuMuSpectra* spectraCorr = FitParticle("psi",
                                                              trigger->String().Data(),
                                                              eventType->String().Data(),
                                                              pairCut->String().Data(),
                                                              centrality->String().Data(),
                                                              *binning,"minv",kTRUE);
            
            AliDebug(1,Form("----fitting done corrected spectra = %p",spectraCorr));
            
            o = 0x0;
            if ( spectraCorr )
            {
              ++nfits;
              
              o = fMergeableCollection->GetObject(id.Data(),spectraCorr->GetName());
              
              AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));
              
              if (o)
              {
                AliWarning(Form("Replacing %s/%s",id.Data(),spectraCorr->GetName()));
                fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectraCorr->GetName()));
              }
              
              Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),spectraCorr);
              
              if ( adoptOK ) std::cout << "+++Spectra " << spectraCorr->GetName() << " adopted" << std::endl;
              else AliError(Form("Could not adopt spectra %s",spectraCorr->GetName()));
              
              StdoutToAliDebug(1,spectraCorr->Print(););
            }
            else AliError("Error creating spectra");
            
            
            if (fitmPt)
            {
              AliDebug(1,"----Fitting mean pt...");
              
              std::cout << "" << std::endl;
              std::cout << "" << std::endl;
              std::cout << "++++++++++++ Fitting mean Pt for " << swhat->String().Data() << " " << "slices" << std::endl; //Uncomment
              
              if ( spectra )
              {
                AliAnalysisMuMuSpectra* spectraMeanPt = FitParticle("psi",
                                                                  trigger->String().Data(),
                                                                  eventType->String().Data(),
                                                                  pairCut->String().Data(),
                                                                  centrality->String().Data(),
                                                                  *binning,"mpt"/*,*spectra*/);
                
                
                
                AliDebug(1,Form("----fitting done spectra = %p",spectraMeanPt));
                o = 0x0;
                
                if ( spectraMeanPt )
                {
                  ++nfits; //Review this
                  
                  o = fMergeableCollection->GetObject(id.Data(),spectraMeanPt->GetName());
                  
                  AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));
                  
                  if (o)
                  {
                    AliWarning(Form("Replacing %s/%s",id.Data(),spectraMeanPt->GetName()));
                    fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectraMeanPt->GetName()));
                  }
                  
                  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),spectraMeanPt);
                  //                  spectraMeanPt->Print();
                  if ( adoptOK ) std::cout << "+++Spectra " << spectraMeanPt->GetName() << " adopted" << std::endl;
                  else AliError(Form("Could not adopt spectra %s",spectraMeanPt->GetName()));
                }
                else AliError("Error creating spectra");

              }
              else std::cout << "Mean pt fit failed: No inv mass spectra for " << swhat->String().Data() << " " << "slices" << std::endl; //Uncomment
              
              
              std::cout << "++++++++++++ Fitting corrected mean Pt for" << " " << swhat->String().Data() << " " << "slices" << std::endl;
              
              if ( spectraCorr )
              {
                AliAnalysisMuMuSpectra* spectraMeanPtCorr  = FitParticle("psi",
                                                                       trigger->String().Data(),
                                                                       eventType->String().Data(),
                                                                       pairCut->String().Data(),
                                                                       centrality->String().Data(),
                                                                       *binning,"mpt"/*,*spectraCorr*/,kTRUE);
                
                
                
                AliDebug(1,Form("----fitting done spectra = %p",spectraMeanPtCorr));
                
                o = 0x0;
                
                if ( spectraMeanPtCorr )
                {
                  ++nfits; //Review this
                  
                  o = fMergeableCollection->GetObject(id.Data(),spectraMeanPtCorr->GetName());
                  
                  AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));
                  
                  if (o)
                  {
                    AliWarning(Form("Replacing %s/%s",id.Data(),spectraMeanPtCorr->GetName()));
                    fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectraMeanPtCorr->GetName()));
                  }
                  
                  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),spectraMeanPtCorr);
                  //                  spectraMeanPtCorr->Print();
                  if ( adoptOK ) std::cout << "+++Spectra " << spectraMeanPtCorr->GetName() << " adopted" << std::endl;
                  else AliError(Form("Could not adopt spectra %s",spectraMeanPtCorr->GetName()));
                  
                }
                else AliError("Error creating spectra");
                
              }
              
              else std::cout << "Corrected mean pt fit failed: No corrected inv mass spectra for " << swhat->String().Data() << " " << "slices" << std::endl;
              
            }
          }
        }
      }
    }
  }
  
  delete whatArray;
  delete triggerArray;
  delete eventTypeArray;
  delete pairCutArray;
  delete centralityArray;
  
  StdoutToAliDebug(1,timer.Print(););

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

  AliAnalysisTriggerScalers ts(RunNumbers(),Config()->OCDBPath());

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
void AliAnalysisMuMu::ShowList(const char* title, const TString& list, const char separator) const
{
  /// Show a list of strings
  TObjArray* parts = list.Tokenize(separator);
  
  std::cout << title << " (" << parts->GetEntries() << ") : " << std::endl;
  
  TIter next(parts);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    std::cout << "    " << str->String().Data() << std::endl;
  }
  
  if ( parts->GetEntries()==0) std::cout << endl;
  
  delete parts;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Print(Option_t* opt) const
{
    /// printout
  std::cout << "Reading from file : " << fFilename.Data() << std::endl;
  
  TString copt(opt);
  
  if (IsSimulation() || SIM() )
  {
    copt += "SIM";
  }

  if ( !IsSimulation() )
  {
    copt += "REAL";
  }
  
  Config()->Print(copt.Data());

  if ( RunNumbers().size() > 1 )
  {
    std::cout << RunNumbers().size() << " runs";
  }
  else
  {
    std::cout << RunNumbers().size() << " run";
  }
  
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
  if ( sopt.Contains("MC") && OC() )
  {
    TString topt(sopt);
    topt.ReplaceAll("MC","");
    OC()->Print(topt.Data());
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
void AliAnalysisMuMu::SetCentralitySelectionList(const char* centralitySelectionList)
{
  /// Set centralities to be used during fitting
  /// centralitySelectionList is a regular expression.

  TObjArray* centralities = BIN()->CreateBinObjArray("centrality");
  TIter next(centralities);
  AliAnalysisMuMuBinning::Range* r;

  TString csl;

  TPRegexp re(centralitySelectionList);
  
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    AliDebug(1,Form("r=%s",r->AsString().Data()));
    if ( re.MatchB(r->AsString()) )
    {
      csl += r->AsString();
      csl += ",";
    }
  }
  
  Config()->SetList(AliAnalysisMuMuConfig::kCentralitySelectionList,IsSimulation(),csl);
  
  delete centralities;
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
void AliAnalysisMuMu::SetNofInputParticles(AliAnalysisMuMuJpsiResult& r)
{
  /// Set the "NofInput" variable(s) of one result
  
  TString hname(Form("MinvUS%s",r.Bin().AsString().Data()));

  TH1* hinput = fMergeableCollection->Histo(Form("/%s/ALL/ANY/V0A/INYRANGE",AliAnalysisMuMuBase::MCInputPrefix()),hname.Data());

  if (!hinput)
  {
    AliError(Form("Got a simulation file where I did not find histogram /%s/ALL/EVERYTHING/ANY/INYRANGE/%s",AliAnalysisMuMuBase::MCInputPrefix(),hname.Data()));

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
  if (!OC()) return 0x0;
  
  return static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(fullpath));
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SelectRunByTrigger(const char* triggerList)
{
  if (!fMergeableCollection || !fCounterCollection) return;
  
  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);
  
  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);
  
  TObjString* srun;
  TObjString* strigger;
  
  TString striggerList(triggerList);
  
  TList* runList = new TList();
  
  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    
    nextTrigger.Reset();
    
    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      if ( !striggerList.Contains(strigger->String().Data()) )
      {
        continue;
      }
      
      ULong64_t n = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                                strigger->String().Data(),"PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00",srun->String().Atoi())));
      if ( n > 0 ) runList->Add(srun);
  
    }
  }
    runList->Sort();
    TIter nextRunOK(runList);
    while ( ( srun = static_cast<TObjString*>(nextRunOK()) ) )
    {
      std::cout << srun->String().Atoi() << std::endl;
    }
    delete runList;

    delete triggers;
    delete runs;
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
        messages.insert(std::make_pair(nmax,static_cast<std::string>(msg.Data())));
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
 
  if (!CC() || !OC()) return;
  
  ReOpen(fFilename,"UPDATE");

  if (OC())
  {
    OC()->Write("OC",TObject::kSingleKey|TObject::kOverwrite);
  }

  ReOpen(fFilename,"READ");
  
  GetCollections(fFilename,fMergeableCollection,fCounterCollection,fBinning,fRunNumbers);
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

  
  AliWarning("Out of date method");
  
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
TH2* AliAnalysisMuMu::ComputeSPDCorrection(const char* type, const char* eventSel, const char* triggerSel, Bool_t bkgReject)
{
  TString stype(type);
  TString evtype(eventSel);
  TString trigtype(triggerSel);
  
  TH2* hNch = static_cast<TH2*>(OC()->Histo(Form("/MCINPUT/%s/%s/V0A/NchVsZVertexVsEta",evtype.Data(),
                                                 trigtype.Data()))); // Input Nch // //"/INPUT/QASPDZPSALL/NchVSEtaVSZVertMC"
  if ( !hNch )
  {
    AliError("No Nch histo found");
    return 0x0;
  }
  TH2* hNtr = static_cast<TH2*>(OC()->Histo(Form("/%s/%s/V0A/TrackletsVsZVertexVsEta",evtype.Data(),
                                                trigtype.Data()))); // Reco tracklets //  //"/RECO/QASPDZPSALL/MB1/NtrVSEtaVSZVertMC"
  if ( !hNtr )
  {
    AliError("No tracklets histo found");
    return 0x0;
  }
  
  TH2* hNtrBkg = static_cast<TH2*>(OC()->Histo(Form("/MCINPUT/%s/%s/V0A/NBkgTrackletsVsZVertexVsEta",evtype.Data(),
                                                 trigtype.Data()))); // Reco tracklets //  //"/RECO/QASPDZPSALL/MB1/NtrVSEtaVSZVertMC"
  if ( !hNtrBkg )
  {
    AliError("No background tracklets histo found");
    return 0x0;
  }

  
  TH2D* hSPDCorr = static_cast<TH2D*>(hNtr->Clone("SPDCorr"));
  TString title("");\
  if ( stype.Contains("oneOverAccEff")) hSPDCorr->SetTitle("SPD 1/AccxEff correction");
  else if ( stype.Contains("AccEffOnly")) hSPDCorr->SetTitle("SPD AccxEff correction");
  else if ( stype.Contains("statOneOverAccEff")) hSPDCorr->SetTitle("SPD 1/AccxEff correction stat. unc.");
  
  for (Int_t i = 1 ; i < hNch->GetNbinsX() ; i++)
  {
    for (Int_t j = 1 ; j < hNch->GetNbinsY() ; j++)
    {
      Int_t n = hNch->GetBin(i,j);
      Double_t nch = hNch->GetBinContent(n);
      Double_t ntr = hNtr->GetBinContent(n);
      Double_t nBkgtr(0.);
      if ( bkgReject ) nBkgtr = hNtrBkg->GetBinContent(n);
      
      Double_t corr(0.),corrErr(0.);
      if ( nch != 0. )
      {
        corr = (ntr - nBkgtr)/nch;
        corrErr = TMath::Max( 1./nch,TMath::Sqrt( corr*(1.-corr)/nch ) );
      }
      
      if ( stype.Contains("oneOverAccEff"))
      {
        if ( corr > 0. )
        {
          hSPDCorr->SetBinContent(n,1./corr);
          hSPDCorr->SetBinError(n,corrErr/TMath::Power(corr,2.));
        }
        else
        {
          hSPDCorr->SetBinContent(n,0.);
          hSPDCorr->SetBinError(n,1.);
        }
        
      }
      else if ( stype.Contains("AccEffOnly"))
      {
        hSPDCorr->SetBinContent(n,corr);
        hSPDCorr->SetBinError(n,corrErr);
      }
      else if ( stype.Contains("statOneOverAccEff"))
      {
        if ( corr != 0. )
        {
          hSPDCorr->SetBinContent(n,(corrErr/TMath::Power(corr,2.)*100)/(1./corr));
        }
        else
        {
          hSPDCorr->SetBinContent(n,-1);
        }

      }
    }
  }
  
  return hSPDCorr;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeFnorm()
{
  /// Compute the CMUL to CINT ratio(s)
  
  if (!CC()) return;
  
  OC()->Prune("/FNORM");
  
  AliAnalysisMuMuFnorm computer(*(CC()),AliAnalysisMuMuFnorm::kMUL,Config()->OCDBPath(),Config()->CompactGraphs());
  
  computer.ComputeFnorm();

  AliMergeableCollection* fnorm = computer.DetachMC();
  
  OC()->Attach(fnorm,"/FNORM/");
  
  Update();
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::ComputeDiffFnormFromHistos(const char* what,const char* quantity,const char* flavour,Bool_t printout)
{
  /// Compute the CMUL to CINT ratio(s) from the histos stored in the OC()
  //FIXME: This is just a patch...
    
  AliAnalysisMuMuBinning* binning = BIN()->Project(what,quantity,flavour);
  if ( !binning )
  {
    AliError(Form("%s-%s-%s binning does not exist",what,quantity,flavour));
    return 0x0;
  }
  TObjArray* dNchdEtas = binning->CreateBinObjArray();
  
  Double_t* binArray = binning->CreateBinArray();
  
  TIter next(dNchdEtas);
  AliAnalysisMuMuBinning::Range* r;
  
  Double_t FNorm(0.);
  Double_t FNormError(0.);
  
  TH1* hFNorm = new TH1F("hFNorm","'Global' normalization factor vs dN_{ch}/d#eta;dN_{ch}/d#eta;FNorm",dNchdEtas->GetEntries(),binArray);
  
  Int_t bin(0);
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    
    TH1* hCMSL = OC()->Histo(Form("/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/CMSL7-B-NOPF-MUON/V0A/%s",
                             Form("EventsIn%s",r->AsString().Data())));
    if ( !hCMSL )
    {
      AliError(Form("No event histo in bin %s found for CMSL7-B-NOPF-MUON",r->AsString().Data()));
      delete binning;
      delete dNchdEtas;
      delete binArray;
      delete hFNorm;
      return 0x0;
    }
    
    TH1* hCMSLandOMUL = OC()->Histo(Form("/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/CMSL7-B-NOPF-MUON&0MUL/V0A/%s",
                                  Form("EventsIn%s",r->AsString().Data())));
    if ( !hCMSLandOMUL )
    {
      AliError(Form("No event histo in bin %s found for CMSL7-B-NOPF-MUON & 0MUL",r->AsString().Data()));
      delete binning;
      delete dNchdEtas;
      delete binArray;
      delete hFNorm;
      return 0x0;
    }
    
    TH1* hCINT = OC()->Histo(Form("/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/CINT7-B-NOPF-ALLNOTRD/V0A/%s",
                                  Form("EventsIn%s",r->AsString().Data())));
    if ( !hCINT )
    {
      AliError(Form("No event histo in bin %s found for CINT7-B-NOPF-ALLNOTRD",r->AsString().Data()));
      delete binning;
      delete dNchdEtas;
      delete binArray;
      delete hFNorm;
      return 0x0;
    }
    
    TH1* hCINTandOMSL = OC()->Histo(Form("/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/CINT7-B-NOPF-ALLNOTRD&0MSL/V0A/%s",
                                  Form("EventsIn%s",r->AsString().Data())));
    if ( !hCINTandOMSL )
    {
      AliError(Form("No event histo in bin %s found for CINT7-B-NOPF-ALLNOTRD & 0MSL",r->AsString().Data()));
      delete binning;
      delete dNchdEtas;
      delete binArray;
      delete hFNorm;
      return 0x0;
    }
  
    FNorm = (hCMSL->GetBinContent(1)/hCMSLandOMUL->GetBinContent(1))*(hCINT->GetBinContent(1)/hCINTandOMSL->GetBinContent(1));
    FNormError = ErrorPropagationAxBoverCxD(hCMSL->GetBinContent(1),hCINT->GetBinContent(1),hCMSLandOMUL->GetBinContent(1),hCINTandOMSL->GetBinContent(1));
    
    if ( printout ) std::cout << r->AsString().Data() << " : " << FNorm << " +- " << FNormError << std::endl;
    
    hFNorm->SetBinContent(++bin,FNorm);
    hFNorm->SetBinError(bin,FNormError);
  }
  
  delete binning;
  delete dNchdEtas;
  delete[] binArray;
  
  return hFNorm;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeDiffFnormFromInt(const char* triggerCluster, const char* eventSelection, AliMergeableCollection* mc, const char* what,const char* quantity,const char* flavour,Bool_t printout)
{
  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }

  TString seventSelection(eventSelection);
  
  TString id(Form("/FNORM-%s/%s/V0A",striggerCluster.Data(),seventSelection.Data()));
  
  AliAnalysisMuMuBinning* binning = BIN()->Project(what,quantity,flavour);
  if ( !binning )
  {
    AliError(Form("%s-%s-%s binning does not exist",what,quantity,flavour));
    return;
  }
  
  TString path(Form("%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data()));
  if ( !mc )
  {
    AliError("Error: No mergeable collection to get Nch histo");
    delete binning;
    return;
  }
  TH1* hNtr = static_cast<TH1*>(mc->Histo(Form("/%s/Nch",path.Data())));
  if ( !hNtr )
  {
    AliError(Form("Error: No /%s/Nch histo in mergeable collection",path.Data()));
    delete binning;
    return;
  }
  Int_t nTrackletsCorrTot = hNtr->Integral();
  
  TObjArray* bin = binning->CreateBinObjArray(what,quantity,flavour);
  Int_t nEntries = bin->GetEntries();
  Double_t* binArray = binning->CreateBinArray();
  Double_t FNormTot(0.);
  Double_t FNormTotError(0.);
  
  TH1* hNorm = OC()->Histo(Form("%s/hFNormInt",id.Data()));
  
  TH1* hFNormTot = new TH1F("hFNormVSdNchdEtaFromInt","Normalization factor vs dN_{ch}/d#eta;dN_{ch}/d#eta;FNorm",nEntries,binArray);
  
  Double_t FNormGlobal = hNorm->GetBinContent(1);
  Double_t FNormGlobalError = hNorm->GetBinError(1);
  
  if ( printout ) std::cout << "Global FNorm = " << FNormGlobal << " + - " << FNormGlobalError << std::endl;
  
  TIter nextBin(bin);
  AliAnalysisMuMuBinning::Range* r;
  Int_t i(1);
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ) //Bin loop
  {
    Double_t nTrklsCorrBin = hNtr->Integral(r->Xmin(),r->Xmax());
    Double_t nTrklsCorrBinFrac = nTrklsCorrBin / nTrackletsCorrTot;
    Double_t nTrklsCorrBinFracError = TMath::Sqrt( TMath::Power(TMath::Sqrt(nTrklsCorrBin) / nTrackletsCorrTot,2.) +
                                                  TMath::Power(TMath::Sqrt(nTrackletsCorrTot)*nTrklsCorrBin / TMath::Power(nTrackletsCorrTot,2.) ,2.) );
    
    FNormTot = FNormGlobal*nTrklsCorrBinFrac;
    FNormTotError = TMath::Sqrt( TMath::Power(FNormGlobalError*nTrklsCorrBinFrac,2.) + TMath::Power(FNormGlobal*nTrklsCorrBinFracError,2.) );

    hFNormTot->SetBinContent(i,FNormTot);
    hFNormTot->SetBinError(i,FNormTotError);
    i++;
    
    if ( printout ) std::cout << "Bin: " << r->AsString().Data() << " ; " << " FNorm = " << FNormTot << " +- " << FNormTotError << std::endl;
  }
  
  TH1* o = fMergeableCollection->Histo(id.Data(),hFNormTot->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hFNormTot->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hFNormTot->GetName()));
  }
  
  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),hFNormTot);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hFNormTot->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hFNormTot->GetName()));
  
  delete binning;
  delete bin;
  delete[] binArray;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeDiffFnormFromCounters(const char* triggerCluster, const char* eventSelection, const char* filePileUpCorr, const char* what,const char* quantity,const char* flavour,Bool_t printout)
{
  /// Compute the CMUL to CINT ratio(s) in 2 steps from the CC(), in bins
  
  TString colType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( colType.Contains("-B-") ) colType = "B";
  else if ( colType.Contains("-S-") ) colType = "S";
  else
  {
    AliError("Unknown collision type");
    return;
  }
  
  TString triggerType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( triggerType.Contains("7-") ) triggerType = "7";
  else if ( triggerType.Contains("8-") ) triggerType = "8";
  else
  {
    AliError("Unknown trigger type");
    return;
  }
  
  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }

  std::cout << striggerCluster.Data() << std::endl;
  
  Bool_t corrPU(kFALSE);
  TObjArray* pUCorr = new TObjArray();
  if ( strlen(filePileUpCorr) > 0 )
  {
    //    std::cout << "Extracting Pile-Up correction factors from " << filePileUpCorr << std::endl;
    char line[1024];
    ifstream in(filePileUpCorr);
    
    while ( in.getline(line,1024,'\n'))
    {
      TString lrun(line);
      TString lvalue(line);
      
      lrun.Remove(0,4);
      lrun.Remove(6,67);
      
      lvalue.Remove(0,57);//71
      
      //      std::cout << "RUN: " << lrun.Data() << " PUFactor = " << lvalue.Data() << std::endl;
      
      pUCorr->Add(new TParameter<Double_t>(lrun.Data(),lvalue.Atof()));
    }
    corrPU = kTRUE;
  }

  TString seventSelection(eventSelection);
  TString sruns = CC()->GetKeyWords("run");
  TObjArray* runs = sruns.Tokenize(",");
  Double_t NofRuns = runs->GetEntries();
  
  TIter nextRun(runs);
  TObjString* s;
  
  AliAnalysisMuMuBinning* binning = BIN()->Project(what,quantity,flavour);
  if ( !binning )
  {
    AliError(Form("%s-%s-%s binning does not exist",what,quantity,flavour));
    return;
  }
  TObjArray* bin = binning->CreateBinObjArray(what,quantity,flavour);
  Double_t* binArray = binning->CreateBinArray();
  Int_t nEntries = bin->GetEntries();
  
  TH1* h;
  TH1* hNofEqMB = new TH1F("hNofEqMBVSdNchdEta","Equivalent MB events per CMUL for vs dN_{ch}/d#eta",bin->GetEntries(),binArray);
  TH1* hFNormTot = new TH1F("hFNormVSdNchdEta","Normalization factor vs dN_{ch}/d#eta;dN_{ch}/d#eta;FNorm",bin->GetEntries(),binArray);

  Double_t* FNormTot = new Double_t[nEntries];
  Double_t* FNormTotError = new Double_t[nEntries];
  
  TString id(Form("/FNORM-%s/%s/V0A",striggerCluster.Data(),seventSelection.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00

  TList* lRun2Reject = new TList();
  lRun2Reject->SetOwner(kTRUE);
  
  Int_t i(0); //dNchdEta bin number
  TObjArray* aCluster = striggerCluster.Tokenize("-");
  TIter nextCluster(aCluster);
  TObjString* striggerClusterS;
  
  TIter nextBin(bin);
  AliAnalysisMuMuBinning::Range* r;
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ) //Bin loop
  {
    FNormTot[i] = 0;
    FNormTotError[i] = 0;
    if ( printout )
    {
      std::cout << "______________________________" << std::endl;
      std::cout << "Bin: " << r->AsString().Data() << std::endl;
    }
    
    h = new TH1F(Form("hFNormVSrun_%s",r->AsString().Data()),Form("Normalization factor vs run for %s ;run;FNorm",r->AsString().Data()),NofRuns,1,NofRuns);
    //Set the run labels
    
    Double_t nCMULBin = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/bin:%s",
                                          seventSelection.Data(),triggerType.Data(),colType.Data(),r->AsString().Data())); //Nof CMUL7/8 events in Bin summed over runs
   //HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
    
    Int_t j(1); //Run label index
    nextRun.Reset();
    while ( ( s = static_cast<TObjString*>(nextRun())) ) //Run loop
    {
      Double_t nCMSL(0.),nCMSLandOMUL(0.);
      nextCluster.Reset();
      while ( (striggerClusterS = static_cast<TObjString*>(nextCluster())) && nCMSL == 0. )
      {
        nCMSL = CC()->GetSum(Form("/event:%s/trigger:CMSL%s-%s-NOPF-%s/centrality:V0A/run:%s/bin:%s",
                                  seventSelection.Data(),triggerType.Data(),colType.Data(),striggerClusterS->GetName(),s->GetName(),r->AsString().Data()));
        
        nCMSLandOMUL = CC()->GetSum(Form("/event:%s/trigger:CMSL%s-%s-NOPF-%s&0MUL/centrality:V0A/run:%s/bin:%s",
                                         seventSelection.Data(),triggerType.Data(),colType.Data(),striggerClusterS->GetName(),s->GetName(),r->AsString().Data()));
      }
      Double_t nCMUL = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/run:%s/bin:%s",
                                seventSelection.Data(),triggerType.Data(),colType.Data(),s->GetName(),r->AsString().Data()));
      
      Double_t nCINT = CC()->GetSum(Form("/event:%s/trigger:CINT%s-%s-NOPF-ALLNOTRD/centrality:V0A/run:%s/bin:%s",
                                         seventSelection.Data(),triggerType.Data(),colType.Data(),s->GetName(),r->AsString().Data()));
      
      Double_t nCINTandOMSL = CC()->GetSum(Form("/event:%s/trigger:CINT%s-%s-NOPF-ALLNOTRD&0MSL/centrality:V0A/run:%s/bin:%s",
                                                seventSelection.Data(),triggerType.Data(),colType.Data(),s->GetName(),r->AsString().Data()));
      
      Double_t FNorm(0.);
      Double_t FNormError(0.);
      Double_t FNormError2(0.);
      Double_t pUfactor = 1.;
      if ( nCMSLandOMUL != 0. && nCINTandOMSL !=0. && nCMSL != 0. && nCINT !=0. )
      {
        if (corrPU)
        {
          TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(pUCorr->FindObject(s->GetName()));
          if ( p ) pUfactor = p->GetVal();
          else
          {
            AliError(Form("Run %s not found in pile-up correction list",s->GetName()));
          }
        }
        
        FNorm = (nCMSL*nCINT)*pUfactor/(nCMSLandOMUL*nCINTandOMSL);
        FNormError = ErrorPropagationAxBoverCxD(nCMSL,nCINT,nCMSLandOMUL,nCINTandOMSL)*pUfactor;
        FNormError2 = AliAnalysisMuMuResult::ErrorABCD(nCMSL, TMath::Sqrt(nCMSL), nCINT, TMath::Sqrt(nCINT), nCMSLandOMUL, TMath::Sqrt(nCMSLandOMUL), nCINTandOMSL, TMath::Sqrt(nCINTandOMSL));
      }
      else
      {
        if ( nCINT == 0 ) std::cout << " Warning: Run " << s->GetName() << " has no MB trigger in this bin" << std::endl;
        
        lRun2Reject->Add(new TObjString(s->GetName()));
        if ( printout ) std::cout << "Run " << s->GetName() << " not used for FNorm cause lack of stats" << std::endl;
        continue;
      }
      FNormTot[i] += FNorm*nCMUL; // This is the sum of equivalent Nof MB per CMUL run by run. NOTE: This sum is NOT always the total equivalent Nof MB per CMUL because in pp 2012 if just one cluster is used at a time this sum is not the sum for all runs
      FNormTotError[i] += TMath::Power(nCMUL*FNormError,2.) + TMath::Power(FNorm*TMath::Sqrt(nCMUL),2.);
      
      if ( printout ) std::cout << "Run " << s->GetName() << " FNorm = " << FNorm << " +- " << FNormError << " (" << FNormError2 << ")" << " ; PUFactor =" << pUfactor << " ; " << "Nof CMUL = " << nCMUL << std::endl;
      
      h->GetXaxis()->SetBinLabel(j,s->GetName());
      h->SetBinContent(j,FNorm);
      h->SetBinError(j++,FNormError);

    }
    
//    std::cout << "NofCMUL in " << i << " = " << nCMULBin << std::endl;
    
    TIter nextRejectRun(lRun2Reject);
    TObjString* run2Rej;
    Double_t nCMULBinRej(0.);
    while ( (run2Rej = static_cast<TObjString*>(nextRejectRun())) )
    {
      nCMULBinRej += CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/bin:%s/run:%s",
                                       seventSelection.Data(),triggerType.Data(),colType.Data(),r->AsString().Data(),
                                       run2Rej->GetName())); //Sum of CMUL7 events from rejected runs
    }
    
    nCMULBin = nCMULBin - nCMULBinRej;
    lRun2Reject->Clear();
    
    FNormTotError[i] =  TMath::Sqrt(TMath::Power(TMath::Sqrt(FNormTotError[i])/nCMULBin,2.) + TMath::Power(FNormTot[i]*TMath::Sqrt(nCMULBin)/TMath::Power(nCMULBin,2.),2.));
    FNormTot[i] = FNormTot[i]/nCMULBin;
    
    std::cout << "Mean FNorm in Bin = " << FNormTot[i]  << " +- " << FNormTotError[i] <<std::endl;
    
    hFNormTot->SetBinContent(i+1,FNormTot[i]);
    hFNormTot->SetBinError(i+1,FNormTotError[i]);
    
    //____
    nCMULBin = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/bin:%s",
                                 seventSelection.Data(),triggerType.Data(),colType.Data(),r->AsString().Data())); //Nof CMUL7/8 events in Bin summed over runs
    
    Double_t nofEqMB = FNormTot[i]*nCMULBin;
    Double_t nofEqMBError = TMath::Sqrt( TMath::Power(FNormTotError[i]*nCMULBin,2.) + TMath::Power(FNormTot[i]*TMath::Sqrt(nCMULBin),2.) );
    
    std::cout << "EqMB in Bin  = " << nofEqMB << " +- " << nofEqMBError << " ; nCMUL = " << nCMULBin << std::endl;
    
    hNofEqMB->SetBinContent(i+1,nofEqMB);
    hNofEqMB->SetBinError(i+1,nofEqMBError);
    //____
    
    TH1* o = fMergeableCollection->Histo(id.Data(),h->GetName());
    
    if (o)
    {
      AliWarning(Form("Replacing %s/%s",id.Data(),h->GetName()));
      fMergeableCollection->Remove(Form("%s/%s",id.Data(),h->GetName()));
    }
    
    Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),h);
    
    if ( adoptOK ) std::cout << "+++FNorm histo " << h->GetName() << " adopted" << std::endl;
    else AliError(Form("Could not adopt FNorm histo %s",h->GetName()));
    
    i++;
    lRun2Reject->Clear();
  }
  
  TH1* o = fMergeableCollection->Histo(id.Data(),hFNormTot->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hFNormTot->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hFNormTot->GetName()));
  }
  
  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),hFNormTot);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hFNormTot->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hFNormTot->GetName()));
  
  
  o = fMergeableCollection->Histo(id.Data(),hNofEqMB->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hNofEqMB->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hNofEqMB->GetName()));
  }
  
  adoptOK = fMergeableCollection->Adopt(id.Data(),hNofEqMB);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hNofEqMB->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hNofEqMB->GetName()));
  
  delete binning;
  delete runs;
  delete aCluster;
  delete bin;
  delete[] binArray;
  delete[] FNormTot;
  delete[] FNormTotError;
  delete lRun2Reject;
  
  return;
 
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeDiffFnormFromGlobal(const char* triggerCluster, const char* eventSelection, const char* what,const char* quantity,
                                                 const char* flavour, Bool_t printout)
{
  TString colType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( colType.Contains("-B-") ) colType = "B";
  else if ( colType.Contains("-S-") ) colType = "S";
  else
  {
    AliError("Unknown collision type");
    return;
  }
  
  TString triggerType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( triggerType.Contains("7-") ) triggerType = "7";
  else if ( triggerType.Contains("8-") ) triggerType = "8";
  else
  {
    AliError("Unknown trigger type");
    return;
  }
  
  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }
  
  TString seventSelection(eventSelection);
  TString id(Form("/FNORM-%s/%s/V0A",striggerCluster.Data(),seventSelection.Data()));
  
  TH1* hFnormGlobal = OC()->Histo(id.Data(),"hFNormVSdNchdEta");
  if( !hFnormGlobal)
  {
    AliError("hFNormVSdNchdEta not found");
    return;
  }
  
  TH1* hFnormGlobalInt = OC()->Histo(id.Data(),"hFNormInt");
  if( !hFnormGlobalInt)
  {
    AliError("hFNormInt not found");
    return;
  }
  Double_t FNormGlobal = hFnormGlobalInt->GetBinContent(1);
  Double_t FNormGlobalError = hFnormGlobalInt->GetBinError(1);
  
  AliAnalysisMuMuBinning* binning = BIN()->Project(what,quantity,flavour);
  if ( !binning )
  {
    AliError(Form("%s-%s-%s binning does not exist",what,quantity,flavour));
    return;
  }
  TObjArray* bin = binning->CreateBinObjArray(what,quantity,flavour);
  
  TH1* hFNormVSNtr = static_cast<TH1*>(hFnormGlobal->Clone());
  hFNormVSNtr->SetTitle("Normalization factor vs dN_{ch}/d#eta;dN_{ch}/d#eta;FNorm");
  hFNormVSNtr->SetName("hFNormVSdNchdEtaFromGlobal");
  
  TH1* hNMBVSNtr = static_cast<TH1*>(hFnormGlobal->Clone());
  hNMBVSNtr->SetTitle("Equivalent MB events per CMUL for vs dN_{ch}/d#eta;NMB");
  hNMBVSNtr->SetName("hNofEqMBVSdNchdEtaFromGlobal");

  
  Double_t nCMULTot = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A",
                                        seventSelection.Data(),triggerType.Data(),colType.Data()));
  
  Double_t nCINTTot = CC()->GetSum(Form("/event:%s/trigger:CINT%s-%s-NOPF-ALLNOTRD/centrality:V0A",
                                     seventSelection.Data(),triggerType.Data(),colType.Data()));
  TIter nextBin(bin);
  AliAnalysisMuMuBinning::Range* r;
  Int_t i(1);
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ) //Bin loop
  {
    Double_t nCMUL = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/bin:%s",
                                       seventSelection.Data(),triggerType.Data(),colType.Data(),r->AsString().Data()));
    
    Double_t nCINT = CC()->GetSum(Form("/event:%s/trigger:CINT%s-%s-NOPF-ALLNOTRD/centrality:V0A/bin:%s",
                                       seventSelection.Data(),triggerType.Data(),colType.Data(),r->AsString().Data()));
    
    Double_t f = nCMUL/nCMULTot;
    Double_t fError = TMath::Sqrt( TMath::Power(TMath::Sqrt(nCMUL)/nCMULTot,2.) + TMath::Power(nCMUL*TMath::Sqrt(nCMULTot)/TMath::Power(nCMULTot,2.),2.) );
    
    Double_t g = nCINT/nCINTTot;
    Double_t gError = TMath::Sqrt( TMath::Power(TMath::Sqrt(nCINT)/nCINTTot,2.) + TMath::Power(nCINT*TMath::Sqrt(nCINTTot)/TMath::Power(nCINTTot,2.),2.) );
    
    Double_t value = FNormGlobal*(g/f);
    Double_t error = TMath::Sqrt( TMath::Power(FNormGlobalError*(g/f),2.) + TMath::Power(FNormGlobal*(gError/f),2.) + TMath::Power(FNormGlobal*g*fError/TMath::Power(f,2.),2.) );
    
    hFNormVSNtr->SetBinContent(i,value);
    hFNormVSNtr->SetBinError(i,error);
    
    if (printout) std::cout << value << " +- " << error << " ; nCMUL = " << nCMUL << std::endl;
    
    hNMBVSNtr->SetBinContent(i,value*nCMUL);
    hNMBVSNtr->SetBinError(i,TMath::Sqrt( TMath::Power(error*nCMUL,2.) + TMath::Power(value*TMath::Sqrt(nCMUL),2.) ));
    
    i++;
  }
  
  TH1* o = fMergeableCollection->Histo(id.Data(),hFNormVSNtr->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hFNormVSNtr->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hFNormVSNtr->GetName()));
  }
  
  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),hFNormVSNtr);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hFNormVSNtr->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hFNormVSNtr->GetName()));
  
  o = fMergeableCollection->Histo(id.Data(),hNMBVSNtr->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hNMBVSNtr->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hNMBVSNtr->GetName()));
  }
  
  adoptOK = fMergeableCollection->Adopt(id.Data(),hNMBVSNtr);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hNMBVSNtr->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hNMBVSNtr->GetName()));

  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeMeanFnorm(const char* triggerCluster, const char* eventSelection, const char* what,const char* quantity,const char* flavour, Bool_t printout)
{
  /// Compute the mean Fnorm and mean NMB from the offline and "rescaled global" methods
  
  TString seventSelection(eventSelection);
  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }
  
  TString triggerType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( triggerType.Contains("7-") ) triggerType = "7";
  else if ( triggerType.Contains("8-") ) triggerType = "8";
  else
  {
    AliError("Unknown trigger type");
    return;
  }
  
  TString colType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( colType.Contains("-B-") ) colType = "B";
  else if ( colType.Contains("-S-") ) colType = "S";
  else
  {
    AliError("Unknown collision type");
    return;
  }

  TH1* hMB = OC()->Histo(Form("/FNORM-%s/%s/V0A/hNofEqMBVSdNchdEta",striggerCluster.Data(),seventSelection.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
  if ( !hMB )
  {
    AliError("Histo hNofEqMBVSdNchdEta not found");
    return;
  }
  
  TH1* hMBG = OC()->Histo(Form("/FNORM-%s/%s/V0A/hNofEqMBVSdNchdEtaFromGlobal",striggerCluster.Data(),seventSelection.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
  if ( !hMBG )
  {
    AliError("Histo hNofEqMBVSdNchdEtaFromGlobal not found");
    return;
  }
  
  AliAnalysisMuMuBinning* binning = BIN()->Project(what,quantity,flavour);
  if ( !binning )
  {
    AliError(Form("%s-%s-%s binning does not exist",what,quantity,flavour));
    return;
  }
  TObjArray* bin = binning->CreateBinObjArray(what,quantity,flavour);
  
  
  TString id(Form("/FNORM-%s/%s/V0A",striggerCluster.Data(),seventSelection.Data()));
  
  TH1* hMBMean = static_cast<TH1*>(hMBG->Clone());
  hMBMean->SetName("hNofEqMBVSdNchdEtaFromMean");
  
  TH1* hFnormMean = static_cast<TH1*>(hMBG->Clone());
  hFnormMean->SetName("hFNormVSdNchdEtaFromMean");
  
  for ( Int_t i = 1 ; i <= hMB->GetNbinsX() ; i++ )
  {
    Double_t Fn = hMB->GetBinContent(i);
    Double_t Fng = hMBG->GetBinContent(i);
    
    Double_t FnE = hMB->GetBinError(i);
    Double_t FngE = hMBG->GetBinError(i);
    
    Double_t meanBin = (Fn + Fng) / 2.;
    Double_t meanBinError = TMath::Sqrt( TMath::Power(FnE/2.,2.) + TMath::Power(FngE/2.,2.) );
    Double_t meanBinSys = TMath::Abs( meanBin - Fn );

//    Double_t meanBin = (Fn/TMath::Power(FnE,2.) + Fng/TMath::Power(FngE,2.)) / ( 1./TMath::Power(FnE,2.) + 1./TMath::Power(FngE,2.) );
//    Double_t meanBinError = 1. / TMath::Sqrt( 1./TMath::Power(FnE,2.) + 1./TMath::Power(FngE,2.) );
//    Double_t meanBinSys = TMath::Sqrt( TMath::Power(Fn - meanBin,2.)/TMath::Power(FnE,2.) + TMath::Power(Fng - meanBin,2.)/TMath::Power(FngE,2.) );
    
    
    hMBMean->SetBinContent(i,meanBin);
    hMBMean->SetBinError(i,meanBinError);
    
    Double_t nCMULBin = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/bin:%s",
                                          seventSelection.Data(),triggerType.Data(),colType.Data(),
                                          static_cast<AliAnalysisMuMuBinning::Range*>(bin->At(i-1))->AsString().Data()));
    
    if (printout) std::cout << meanBinSys/nCMULBin << std::endl;
    
    Double_t meanFnBin = meanBin/nCMULBin;
    Double_t meanFnBinError = TMath::Sqrt( TMath::Power(meanBinError/nCMULBin,2) + TMath::Power(meanBin/TMath::Power(nCMULBin,2.),2) );
    
    if (printout) std::cout << meanBinSys/nCMULBin/meanFnBin << std::endl;

    if (printout) std::cout << meanFnBin << " +- " << meanFnBinError << std::endl;
    
    hFnormMean->SetBinContent(i,meanFnBin);
    hFnormMean->SetBinError(i,meanFnBinError);
  }

  TH1* o = fMergeableCollection->Histo(id.Data(),hMBMean->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hMBMean->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hMBMean->GetName()));
  }
  
  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),hMBMean);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hMBMean->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hMBMean->GetName()));
  
  
  o = fMergeableCollection->Histo(id.Data(),hFnormMean->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hFnormMean->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hFnormMean->GetName()));
  }
  
  adoptOK = fMergeableCollection->Adopt(id.Data(),hFnormMean);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hFnormMean->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hFnormMean->GetName()));


}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeIntFnormFromCounters(const char* triggerCluster, const char* eventSelection, const char* filePileUpCorr, Bool_t printout)
{
  /// Compute the CMUL to CINT ratio(s) in 2 steps from the CC(), integrated
  
  TString colType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( colType.Contains("-B-") ) colType = "B";
  else if ( colType.Contains("-S-") ) colType = "S";
  else
  {
    AliError("Unknown collision type");
    return;
  }
  
  TString triggerType(First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data());
  if ( triggerType.Contains("7-") ) triggerType = "7";
  else if ( triggerType.Contains("8-") ) triggerType = "8";
  else
  {
    AliError("Unknown trigger type");
    return;
  }

  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }

  TString seventSelection(eventSelection);
  TString sruns = CC()->GetKeyWords("run");
  TObjArray* runs = sruns.Tokenize(",");
  Double_t NofRuns = runs->GetEntries();
  
  Bool_t corrPU(kFALSE);
  TObjArray* pUCorr = new TObjArray();
  if ( strlen(filePileUpCorr) > 0 )
  {
//    std::cout << "Extracting Pile-Up correction factors from " << filePileUpCorr << std::endl;
    char line[1024];
    ifstream in(filePileUpCorr);
    
    while ( in.getline(line,1024,'\n'))
    {
      TString lrun(line);
      TString lvalue(line);
      
      lrun.Remove(0,4);
      lrun.Remove(6,67);
      
      lvalue.Remove(0,57);//71
      
//      std::cout << "RUN: " << lrun.Data() << " PUFactor = " << lvalue.Data() << std::endl;
      
      pUCorr->Add(new TParameter<Double_t>(lrun.Data(),lvalue.Atof()));
    }
    corrPU = kTRUE;
  }
  
  TIter nextRun(runs);
  TObjString* s;
  
  TH1* h;
  
  Double_t FNormTot(0.);
  Double_t FNormTotError(0.);
  
  TString id(Form("/FNORM-%s/%s/V0A",striggerCluster.Data(),seventSelection.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
  
  h = new TH1F("hFNormIntVSrun","Integrated Normalization factor vs run;run;FNorm",NofRuns,1.,NofRuns);
  
  Double_t nCMULTot = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A",
                                        seventSelection.Data(),triggerType.Data(),colType.Data())); //Total Nof CMUL7 events
  
  TObjArray* aCluster = striggerCluster.Tokenize("-");
  TIter nextCluster(aCluster);
  TObjString* striggerClusterS;
  
  TList* lRun2Reject = new TList();
  lRun2Reject->SetOwner(kTRUE);
  
  nextRun.Reset();
  Int_t i(0);//Run label index
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Computing FNorm" << std::endl;
  while ( ( s = static_cast<TObjString*>(nextRun())) ) //Run loop
  {
    Double_t nCMSL(0.),nCMSLandOMUL(0.);
    nextCluster.Reset();
    while ( (striggerClusterS = static_cast<TObjString*>(nextCluster())) && nCMSL == 0. )
    {
      nCMSL = CC()->GetSum(Form("/event:%s/trigger:CMSL%s-%s-NOPF-%s/centrality:V0A/run:%s",
                                seventSelection.Data(),triggerType.Data(),colType.Data(),striggerClusterS->GetName(),s->GetName()));
      
      nCMSLandOMUL = CC()->GetSum(Form("/event:%s/trigger:CMSL%s-%s-NOPF-%s&0MUL/centrality:V0A/run:%s",
                                       seventSelection.Data(),triggerType.Data(),colType.Data(),striggerClusterS->GetName(),s->GetName()));
    }
    
    Double_t nCINT = CC()->GetSum(Form("/event:%s/trigger:CINT%s-%s-NOPF-ALLNOTRD/centrality:V0A/run:%s",
                                       seventSelection.Data(),triggerType.Data(),colType.Data(),s->GetName()));
    
    Double_t nCINTandOMSL = CC()->GetSum(Form("/event:%s/trigger:CINT%s-%s-NOPF-ALLNOTRD&0MSL/centrality:V0A/run:%s",
                                              seventSelection.Data(),triggerType.Data(),colType.Data(),s->GetName()));
    
    Double_t nCMUL = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/run:%s",
                                       seventSelection.Data(),triggerType.Data(),colType.Data(),s->GetName()));
    
    Double_t FNorm(0.),FNormError(0.);
    Double_t pUfactor = 1.;
    if ( nCMSLandOMUL != 0. && nCINTandOMSL !=0. && nCMSL != 0. && nCINT !=0. )
    {
      if (corrPU)
      {
        TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(pUCorr->FindObject(s->GetName()));
        if ( p ) pUfactor = p->GetVal();
        else
        {
          AliError(Form("Run %s not found in pile-up correction list",s->GetName()));
        }
      }
      FNorm = (nCMSL*nCINT)*pUfactor/(nCMSLandOMUL*nCINTandOMSL);
      FNormError = ErrorPropagationAxBoverCxD(nCMSL,nCINT,nCMSLandOMUL,nCINTandOMSL)*pUfactor;
    }
    else
    {
      if ( nCINT == 0 ) std::cout << " Warning: Bad run " << s->GetName() << " has no MB trigger in this bin. Remove from analysis" << std::endl;
      
      lRun2Reject->Add(new TObjString(s->GetName()));
      if ( printout ) std::cout << "Run " << s->GetName() << " not used for FNorm cause lack of stats" << std::endl;
      continue;
    }
    
    FNormTot += FNorm*nCMUL; // This is the sum of equivalent Nof MB per CMUL run by run. NOTE: This sum is NOT always the total equivalent Nof MB per CMUL because in pp 2012 if just one cluster is used at a time this sum is not the sum for all runs
    FNormTotError += TMath::Power(nCMUL*FNormError,2.) + TMath::Power(FNorm*TMath::Sqrt(nCMUL),2.);
    
    if ( printout ) std::cout << "Run " << s->GetName() << " FNorm = " << FNorm << " +- " << FNormError << " ; PUFactor =" << pUfactor << " ; " << "Nof CMUL = " << nCMUL << std::endl;
    
    h->GetXaxis()->SetBinLabel(++i,s->GetName());
    h->SetBinContent(i,FNorm);
    h->SetBinError(i,FNormError);
    
  }
  
  TIter nextRejectRun(lRun2Reject);
  TObjString* run2Rej;
  Double_t nCMULTotRej(0.);
  while ( (run2Rej = static_cast<TObjString*>(nextRejectRun())) )
  {
    nCMULTotRej += CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A/run:%s",
                                     seventSelection.Data(),triggerType.Data(),colType.Data(),
                                     run2Rej->GetName())); //Sum of CMUL7 events from rejected runs
  }
  
  nCMULTot = nCMULTot - nCMULTotRej;
  
  TH1* o = fMergeableCollection->Histo(id.Data(),h->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),h->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),h->GetName()));
  }
  
  Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),h);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << h->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",h->GetName()));
  
  //___
  
  FNormTotError =  TMath::Sqrt(TMath::Power(TMath::Sqrt(FNormTotError)/nCMULTot,2.) + TMath::Power(FNormTot*TMath::Sqrt(nCMULTot)/TMath::Power(nCMULTot,2.),2.));
  
  FNormTot = FNormTot/nCMULTot; // nCMULTot is here nCMULTot - nCMULTotRej
  
  std::cout << "FNorm = " << FNormTot << " +- " << FNormTotError << std::endl;
  
  TH1* hFNormTot = new TH1F("hFNormInt","Global Normalization factor",1,0.,1.);
  hFNormTot->SetBinContent(1,FNormTot);
  hFNormTot->SetBinError(1,FNormTotError);
  
  o = fMergeableCollection->Histo(id.Data(),hFNormTot->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hFNormTot->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hFNormTot->GetName()));
  }
  
  adoptOK = fMergeableCollection->Adopt(id.Data(),hFNormTot);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hFNormTot->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hFNormTot->GetName()));
  
  //___
  TH1* hNEqMB = new TH1F("hNEqMB","Equivalent number of MB events per CMUL",1,0.,1.);
  
  nCMULTot = CC()->GetSum(Form("/event:%s/trigger:CMUL%s-%s-NOPF-MUON/centrality:V0A",
                               seventSelection.Data(),triggerType.Data(),colType.Data())); //Total Nof CMUL7 events
  
  Double_t nofEqMB = FNormTot*nCMULTot;
  Double_t nofEqMBError = TMath::Sqrt( TMath::Power(FNormTotError*nCMULTot,2.) + TMath::Power(FNormTot*TMath::Sqrt(nCMULTot),2.) );
  
  std::cout << "EqMB = " << nofEqMB << " +- " << TMath::Sqrt(nofEqMBError) << std::endl;
  
  hNEqMB->SetBinContent(1,nofEqMB);
  hNEqMB->SetBinError(1,nofEqMBError);
  
  o = fMergeableCollection->Histo(id.Data(),hNEqMB->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s",id.Data(),hNEqMB->GetName()));
    fMergeableCollection->Remove(Form("%s/%s",id.Data(),hNEqMB->GetName()));
  }
  
  adoptOK = fMergeableCollection->Adopt(id.Data(),hNEqMB);
  
  if ( adoptOK ) std::cout << "+++FNorm histo " << hNEqMB->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt FNorm histo %s",hNEqMB->GetName()));
  
  //___
  
  delete runs;
  delete lRun2Reject;
  delete aCluster;
  
  return;
  
}

//_____________________________________________________________________________
void AliAnalysisMuMu::PlotYiedWSyst(const char* triggerCluster)
{
  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }
  
  TString path(Form("%s/%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,kFALSE)).Data()));
  
  TH1* hY = OC()->Histo(Form("/RESULTS-%s/%s",striggerCluster.Data(),path.Data()),"hJPsiYieldVSdNchdEtaRelative");
  if ( !hY )
  {
    AliError("No yield found");
    return;
  }

  TString id(Form("/TESTSYST/%s",path.Data()));
  
  TH1* hYSyst = static_cast<TH1*>(hY->Clone("RelYieldSyst"));
  if ( !hYSyst )
  {
    AliError("No systematic found");
    return;
  }

  TH1* hS = OC()->Histo(id.Data(),"yield_Systematics");
  
  for ( Int_t i = 1 ; i <= hY->GetNbinsX() ; i++ )
  {
    hYSyst->SetBinError(i,hS->GetBinContent(i)*hY->GetBinContent(i)/100.);
  }
    
  hY->Draw();
  hYSyst->Draw("same");
}


//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeJpsiYield(AliMergeableCollection* oc, Bool_t relative, const char* fNormType, const char* triggerCluster,
                                       const char* whatever, const char* sResName, AliMergeableCollection* ocMBTrigger, Double_t mNTrCorrection)
{
  // This method is suppossed to be used from the file with the counters, oc is the AliMergeableCollection of the file with the histograms (if separated, which is better since we do not need the minv,mean pt... analysis in CINT&0MUL... triggers)
  // ocMBTrigger is the mergeableCollection with the MB trigger dNchdEta plot (migth be the same as oc, in which case we set ocMBTrigger=0x0)
  //FIXME::Make it general
 
  TString sfNormType(fNormType);
  TString swhat("");
  TString sres("");
  TString swhatever(whatever);
  if ( swhatever.Contains("DNCHDETA"))
  {
    swhat = "dNchdEta";
    if ( strlen(sResName) > 0 ) sres = sResName; //"PSIPSIPRIMECB2VWGINDEPTAILS";
  }
  else if ( swhatever.Contains("NTRCORR") )
  {
    swhat = "Nch";
    if ( strlen(sResName) > 0 ) sres = sResName; //"PSIPSIPRIMECB2VWG_2.0_5.0";
  }
  
  if ( IsSimulation() )
  {
    AliError("Cannot compute J/Psi yield: Is a simulation file");
    return;
  }
  
  TString striggerCluster(triggerCluster);
  if ( striggerCluster.Contains("MUON") && !striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON";
  else if ( striggerCluster.Contains("ALLNOTRD") && !striggerCluster.Contains("MUON") ) striggerCluster = "ALLNOTRD";
  else if ( striggerCluster.Contains("MUON") && striggerCluster.Contains("ALLNOTRD") ) striggerCluster = "MUON-ALLNOTRD";
  else
  {
    AliError("Unknown trigger cluster");
    return;
  }
  
  TString path(Form("%s/%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,kFALSE)).Data()));
  
  AliMergeableCollection* mc;
  if ( !oc ) mc = OC();
  else mc = oc;
  
  Double_t bR = 0.0593; // BR(JPsi->mu+mu-)
  Double_t bRerror = 0.0006 ;
  
  //_________Integrated yield
  AliAnalysisMuMuSpectra* sInt = static_cast<AliAnalysisMuMuSpectra*>(mc->GetObject(Form("/%s/%s",path.Data(),"PSI-INTEGRATED-AccEffCorr")));
  if ( !sInt )
  {
    AliError(Form("No spectra %s found in %s","PSI-INTEGRATED-AccEffCorr",path.Data()));
    return;
  }
  
  AliAnalysisMuMuBinning* b = new AliAnalysisMuMuBinning;
  b->AddBin("psi","INTEGRATED");
  
  AliAnalysisMuMuBinning::Range* bin = static_cast<AliAnalysisMuMuBinning::Range*>(b->CreateBinObjArray()->At(0));
  
  AliAnalysisMuMuResult* result = sInt->GetResultForBin(*bin);
  if ( !result )
  {
    AliError(Form("No result for bin %s found in %s",bin->AsString().Data(),"PSI-INTEGRATED-AccEffCorr"));
    return;
  }
  
//  if ( strlen(sResName) > 0/*sResName.Sizeof() > 0*/ )
//  {
//    result = result->SubResult(sres.Data());//INDEPTAILS
//    if ( !result )
//    {
//      AliError(Form("No subresult %s found in %s",sres.Data(),path.Data()));
//      return;
//    }
//  }
  
  Double_t NofJPsiTot = result->GetValue("NofJPsi");
  Double_t NofJPsiTotError = result->GetErrorStat("NofJPsi");
  
  TH1* hMBTot = OC()->Histo(Form("/FNORM-%s/PSALL/V0A/hNEqMB",striggerCluster.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
  if ( !hMBTot )
  {
    AliError(Form("No eq Nof MB events found in %s",Form("/FNORM-%s/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/V0A/hNEqMB",striggerCluster.Data())));
    return;
  }
  
  Double_t nEqMBTot = hMBTot->GetBinContent(1);
  Double_t nEqMBTotError = hMBTot->GetBinError(1);
  
  Double_t yieldInt = NofJPsiTot/(nEqMBTot*bR);
  Double_t yieldIntError = TMath::Sqrt(TMath::Power(NofJPsiTotError/(nEqMBTot*bR),2.) +
                                       TMath::Power(nEqMBTotError*NofJPsiTot*bR/TMath::Power(nEqMBTot*bR,2.),2.) +
                                       TMath::Power(NofJPsiTot*nEqMBTot*bRerror/TMath::Power(nEqMBTot*bR,2.),2.));
  
  std::cout << "Integrated yield = " << yieldInt << " +- " << yieldIntError << std::endl;
  
  TH1* hYint = new TH1F("hJPsiYieldInt","Integrated J/#psi yield",1,0.,1.);
  hYint->SetBinContent(1,yieldInt);
  hYint->SetBinError(1,yieldIntError);
  
  TH1* o = mc->Histo(Form("/RESULTS-%s/%s",striggerCluster.Data(),path.Data()),hYint->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing /RESULTS-%s/%s/%s",striggerCluster.Data(),path.Data(),hYint->GetName()));
    mc->Remove(Form("/RESULTS-%s/%s/%s",striggerCluster.Data(),path.Data(),hYint->GetName()));
  }
  
  Bool_t adoptOK = mc->Adopt(Form("/RESULTS-%s/%s",striggerCluster.Data(),path.Data()),hYint);
  
  if ( adoptOK ) std::cout << "+++Yield histo " << hYint->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt Yield histo %s",hYint->GetName()));

  delete b;
  
  //_____Differential yield
  
  AliAnalysisMuMuSpectra* s = static_cast<AliAnalysisMuMuSpectra*>(mc->GetObject(Form("/%s/%s",path.Data(),whatever)));
  if ( !s )
  {
    AliError(Form("No spectra %s found in %s",whatever,path.Data()));
    return;
  }
  
  std::cout << "Number of J/Psi:" << std::endl;
  TH1* hry = s->Plot("NofJPsi",sres.Data(),kFALSE);//INDEPTAILS //Number of Jpsi
  
  std::cout << "" << std::endl;
  
//  std::cout << "Equivalent number of MB events:" << std::endl;
  TH1* hMB(0x0);
  if ( sfNormType.Contains("offline") )
  {
    hMB = OC()->Histo(Form("/FNORM-%s/PSALL/V0A/hNofEqMBVSdNchdEta",striggerCluster.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
    if ( !hMB )
    {
      AliError("Histo hNofEqMBVSdNchdEta not found");
      return;
    }
  }
  else if ( sfNormType.Contains("global") )
  {
    hMB = OC()->Histo(Form("/FNORM-%s/PSALL/V0A/hNofEqMBVSdNchdEtaFromGlobal",striggerCluster.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
    if ( !hMB )
    {
      AliError("Histo hNofEqMBVSdNchdEtaFromGlobal not found");
      return;
    }
  }
  else if ( sfNormType.Contains("mean") )
  {
    hMB = OC()->Histo(Form("/FNORM-%s/PSALL/V0A/hNofEqMBVSdNchdEtaFromMean",striggerCluster.Data()));//HASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00
    if ( !hMB )
    {
      AliError("Histo hNofEqMBVSdNchdEtaFromMean not found");
      return;
    }
  }
  else
  {
    AliError("Dont know what Fnorm use");
    return;
  }
  
  TH1* hy;
  if ( relative )
  {
    TString path2(Form("/%s/%s/%s",
                      First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                      First(Config()->GetList(AliAnalysisMuMuConfig::kMinbiasTriggerList,kFALSE)).Data(),
                      First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data()));
    
    TH1* hdNch;
    if ( ocMBTrigger ) hdNch = ocMBTrigger->Histo(path2.Data(),swhat.Data());//dNchdEta
    else hdNch = mc->Histo(path2.Data(),swhat.Data());//dNchdEta
    
    const TArrayD* binArray = hry->GetXaxis()->GetXbins();
    Int_t size = binArray->GetSize();
    Double_t* axis = new Double_t[size];
    for ( Int_t k = 0 ; k < size ; k++ )
    {
      axis[k] = binArray->At(k)/(hdNch->GetMean()*(1 - mNTrCorrection));
    }
    
    hy = new TH1D("hJPsiYieldVSdNchdEtaRelative","Relative J/#psi yield vs dN_{ch}/d#eta;dN_{ch}/d#eta/<dN_{ch}/d#eta>;Y^{J/#psi}/Y^{J/#psi}_{int}",size-1,axis);
    delete axis;
  }
  else
  {
    hy = static_cast<TH1D*>(hry->Clone("hJPsiYieldVSdNchdEta"));
    hy->SetTitle("J/#psi yield vs dN_{ch}/d#eta");
    hy->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    hy->GetYaxis()->SetTitle("Y^{J/#psi}");
  }
                                    // AccxEff(from rel diff or paper)  // Signal extraction
//  Double_t systNofJpsiBin[9] = {TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.01,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.008,2.) ),TMath::Sqrt( TMath::Power(0.022,2.) + TMath::Power(0.007,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.008,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.007,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.009,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.008,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.016 ,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.033,2.) )}; //FIXME: find a way to give this as input
//  Double_t systFNorm[9] = {0.003,0.001,0.002,0.003,0.002,0.004,0.011,0.012,0.071};
//  Double_t systPU[9] = {0.00,0.01,0.012,0.014,0.014,0.019,0.020,0.021,0.040}; //_______pPb
       // AccxEff(from paper)
//  Double_t systNofJpsiTot = 0.015;
  
//  Double_t systNofJpsiBin[9] = {TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.007,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.006,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.005,2.) ),TMath::Sqrt( TMath::Power(0.028,2.) + TMath::Power(0.006,2.) ),TMath::Sqrt( TMath::Power(0.016,2.) + TMath::Power(0.004,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.004,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.006,2.) ),TMath::Sqrt( TMath::Power(0.024,2.) + TMath::Power(0.005 ,2.) ),TMath::Sqrt( TMath::Power(0.015,2.) + TMath::Power(0.016,2.) )}; //FIXME: find a way to give this as input
//  Double_t systFNorm[9] = {0.005,0.004,0.004,0.004,0.003,0.002,0.002,0.04,0.04};
//  Double_t systPU[9] = {0.00,0.007,0.015,0.011,0.014,0.018,0.014,0.011,0.020}; //______Pbp
//  Double_t systNofJpsiTot = 0.015;
  
//  Double_t systNofJpsiBin[9] = {TMath::Sqrt( TMath::Power(0.034,2.) + TMath::Power(0.005,2.) ),TMath::Sqrt( TMath::Power(0.017,2.) + TMath::Power(0.005,2.) ),TMath::Sqrt( TMath::Power(0.017,2.) + TMath::Power(0.004,2.) ),TMath::Sqrt( TMath::Power(0.017,2.) + TMath::Power(0.005,2.) ),TMath::Sqrt( TMath::Power(0.042,2.) + TMath::Power(0.002,2.) ),TMath::Sqrt( TMath::Power(0.063,2.) + TMath::Power(0.014,2.) ),TMath::Sqrt( TMath::Power(0.094,2.) + TMath::Power(0.009,2.) ),TMath::Sqrt( TMath::Power(0.00,2.) + TMath::Power(0.00 ,2.) ),TMath::Sqrt( TMath::Power(0.00,2.) + TMath::Power(0.00,2.) )}; //FIXME: find a way to give this as input
//  Double_t systFNorm[9] = {0.004,0.019,0.002,0.012,0.048,0.063,0.082,0.000,0.000};
//  Double_t systPU[9] = {0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00}; //______pp |eta|<0.5
//  Double_t systNofJpsiTot = 0.017;
  
  Double_t systNofJpsiBin[9] = {TMath::Sqrt( TMath::Power(0.037,2.) + TMath::Power(0.002,2.) ),TMath::Sqrt( TMath::Power(0.021,2.) + TMath::Power(0.002,2.) ),TMath::Sqrt( TMath::Power(0.022,2.) + TMath::Power(0.002,2.) ),TMath::Sqrt( TMath::Power(0.017,2.) + TMath::Power(0.002,2.) ),TMath::Sqrt( TMath::Power(0.019,2.) + TMath::Power(0.001,2.) ),TMath::Sqrt( TMath::Power(0.036,2.) + TMath::Power(0.002,2.) ),TMath::Sqrt( TMath::Power(0.042,2.) + TMath::Power(0.001,2.) ),TMath::Sqrt( TMath::Power(0.039,2.) + TMath::Power(0.012 ,2.) ),TMath::Sqrt( TMath::Power(0.000,2.) + TMath::Power(0.000,2.) )}; //FIXME: find a way to give this as input
  Double_t systFNorm[9] = {0.026,0.002,0.015,0.019,0.012,0.030,0.015,0.119,0.000};
  Double_t systPU[9] = {0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00}; //______pp |eta|<1
  Double_t systNofJpsiTot = 0.017;
  
  for ( Int_t i = 1 ; i <= hy->GetNbinsX() ; i++ )
  {
    Double_t yield = hry->GetBinContent(i)/(hMB->GetBinContent(i)*bR);
    Double_t yieldError = TMath::Sqrt(TMath::Power(hry->GetBinError(i)/(hMB->GetBinContent(i)*bR),2.) +
                                      TMath::Power(hMB->GetBinError(i)*hry->GetBinContent(i)*bR/TMath::Power(hMB->GetBinContent(i)*bR,2.),2.) +
                                      TMath::Power(hry->GetBinContent(i)*hMB->GetBinContent(i)*bRerror/TMath::Power(hMB->GetBinContent(i)*bR,2.),2.));
    
//    std::cout << "Differential yield bin " << i << " = " << yield << " +- " << yieldError << std::endl;
    
    if ( relative )
    {
      yieldError = TMath::Sqrt(TMath::Power(yieldError/yieldInt,2.) + TMath::Power((yield*yieldIntError)/TMath::Power(yieldInt,2.),2.));
      yield /= yieldInt;
      
//      std::cout << "relative yield bin " << i << " = " << yield << " +- " << yieldError << std::endl;
      Double_t sNJpsiBin = hry->GetBinContent(i)*systNofJpsiBin[i-1];
      Double_t sNJpsiTot = NofJPsiTot*systNofJpsiTot;
      Double_t sMBBin = hMB->GetBinContent(i)*systFNorm[i-1];
      Double_t sMBTot = nEqMBTot*0.01;
      
      Double_t syst = TMath::Sqrt( TMath::Power((sNJpsiBin/NofJPsiTot)*(nEqMBTot/hMB->GetBinContent(i)),2.) + TMath::Power((hry->GetBinContent(i)*sNJpsiTot/TMath::Power(NofJPsiTot,2.))*(nEqMBTot/hMB->GetBinContent(i)),2.) + TMath::Power((hry->GetBinContent(i)/NofJPsiTot)*(sMBTot/hMB->GetBinContent(i)),2.) + TMath::Power((hry->GetBinContent(i)/NofJPsiTot)*(sMBBin*nEqMBTot/TMath::Power(hMB->GetBinContent(i),2.)),2.) );
      
      std::cout << "sys" << syst/yield << " w/pu = " << TMath::Sqrt( TMath::Power(syst/yield,2.) + TMath::Power(systPU[i-1],2.)) << std::endl;
      std::cout << yield << " +- " << yieldError << std::endl;
    }

    hy->SetBinContent(i,yield);
    hy->SetBinError(i,yieldError);
  }

  o = mc->Histo(Form("/RESULTS-%s/%s",striggerCluster.Data(),path.Data()),hy->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing %s/%s","/RESULTS-%s/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/V0A",hy->GetName()));
    mc->Remove(Form("/RESULTS-%s/%s/%s",striggerCluster.Data(),path.Data(),hy->GetName()));
  }
  
  adoptOK = mc->Adopt(Form("/RESULTS-%s/%s",striggerCluster.Data(),path.Data()),hy);
  
  if ( adoptOK ) std::cout << "+++Yield histo " << hy->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt Yield histo %s",hy->GetName()));

  
  
  delete hry;

  
  return;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeJpsiMPt(Bool_t relative, const char* whatever, const char* sResName, AliMergeableCollection* ocMBTrigger, Double_t mNTrCorrection)
{
  // ocMBTrigger is the mergeableCollection with the MB trigger dNchdEta plot (migth be the same as oc, in which case we set ocMBTrigger=0x0)
  //FIXME::Make it general
  
  
  TString swhat("");
  TString sres("");
  TString swhatever(whatever);
//  if ( swhatever.Contains("DNCHDETA"))
//  {
//    swhat = "dNchdEta";
//    sres = "MPT2CB2VWGPOL2INDEPTAILS";
//  }
//  else
    if ( swhatever.Contains("NTRCORR") )
  {
    swhat = "Nch";
    if ( strlen(sResName) > 0 ) sres = sResName; //sres = "MPTPSIPSIPRIMECB2VWG_BKGMPTPOL2";
  }

  if ( IsSimulation() )
  {
    AliError("Cannot compute J/Psi yield: Is a simulation file");
    return;
  }
  
  TString path(Form("%s/%s/%s/%s",
                    First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kDimuonTriggerList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data(),
                    First(Config()->GetList(AliAnalysisMuMuConfig::kPairSelectionList,kFALSE)).Data()));
  
  //_________Integrated mean pt
  AliAnalysisMuMuSpectra* sInt = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/%s/%s",path.Data(),"PSI-INTEGRATED-AccEffCorr-MeanPtVsMinvUS")));
  if ( !sInt )
  {
    AliError(Form("No spectra %s found in %s","PSI-INTEGRATED-AccEffCorr-MeanPtVsMinvUS",path.Data()));
    return;
  }
  
  AliAnalysisMuMuBinning* b = new AliAnalysisMuMuBinning;
  b->AddBin("psi","INTEGRATED");
  
  AliAnalysisMuMuBinning::Range* bin = static_cast<AliAnalysisMuMuBinning::Range*>(b->CreateBinObjArray()->At(0));
  
  AliAnalysisMuMuResult* result = sInt->GetResultForBin(*bin);
  if ( !result )
  {
    AliError(Form("No result for bin %s found in spectra %s",bin->AsString().Data(),sInt->GetName()));
    return;
  }
  
//  if ( sres.Sizeof() > 0 )
//  {
//    result = result->SubResult(sres.Data());
//    if ( !result )
//    {
//      AliError(Form("No subresult %s found in result",result->GetName()));
//      return;
//    }
//  AliAnalysisMuMuResult* subresult = result->SubResult(sres.Data());//"MPT2CB2VWGPOL2INDEPTAILS"
//  if ( !subresult )
//  {
//    AliError(Form("No subresult MPT2CB2VWGPOL2 found in result %s",result->GetName()));
//    return;
//  }
    
//  }
//  Double_t JPsiMPtTot = subresult->GetValue("MeanPtJPsi");
//  Double_t JPsiMPtTotError = subresult->GetErrorStat("MeanPtJPsi");
  
  Double_t JPsiMPtTot = result->GetValue("MeanPtJPsi");
  Double_t JPsiMPtTotError = result->GetErrorStat("MeanPtJPsi");
 
  TH1* hMPtint = new TH1F("hJPsiMPtInt","Integrated J/#psi mean p_{T}",1,0.,1.);
  hMPtint->SetBinContent(1,JPsiMPtTot);
  hMPtint->SetBinError(1,JPsiMPtTotError);
  
  TH1* o = OC()->Histo(Form("/RESULTS/%s",path.Data()),hMPtint->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing /RESULTS/%s/%s",path.Data(),hMPtint->GetName()));
    OC()->Remove(Form("/RESULTS/%s/%s",path.Data(),hMPtint->GetName()));
  }
  
  Bool_t adoptOK = OC()->Adopt(Form("/RESULTS/%s",path.Data()),hMPtint);
  
  if ( adoptOK ) std::cout << "+++Mean Pt histo " << hMPtint->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt Mean Pt histo %s",hMPtint->GetName()));
  
  delete b;

   //_____Differential mean pt
  
  AliAnalysisMuMuSpectra* s = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/%s/%s",path.Data(),whatever)));
  if ( !s )
  {
    AliError(Form("No spectra %s found in %s",whatever,path.Data()));
    return;
  }
  
  std::cout << "Mean pt of J/Psi:" << std::endl;
  TH1* hrmPt = s->Plot("MeanPtJPsi",sres.Data(),kFALSE); //MPT2CB2VWGPOL2INDEPTAILS//mean pt of Jpsi
  std::cout << "" << std::endl;
  
  Double_t ptInt,ptIntError;
  TH1* hmPt;
  if ( relative )
  {
    TString path2(Form("/%s/%s/%s",
                       First(Config()->GetList(AliAnalysisMuMuConfig::kEventSelectionList,kFALSE)).Data(),
                       First(Config()->GetList(AliAnalysisMuMuConfig::kMinbiasTriggerList,kFALSE)).Data(),
                       First(Config()->GetList(AliAnalysisMuMuConfig::kCentralitySelectionList,kFALSE)).Data()));
    
    TH1* hdNch;
    if ( ocMBTrigger ) hdNch = ocMBTrigger->Histo(path2.Data(),swhat.Data());
    else hdNch = OC()->Histo(path2.Data(),swhat.Data());
    
    const TArrayD* binArray = hrmPt->GetXaxis()->GetXbins();
    Int_t size = binArray->GetSize();
    Double_t* axis = new Double_t[size];
    for ( Int_t k = 0 ; k < size ; k++ )
    {
      axis[k] = binArray->At(k)/(hdNch->GetMean()*(1 - mNTrCorrection));
    }
    
    hmPt = new TH1D("hJPsiMeanPtVSdNchdEtaRelative","Relative J/#psi mean p_{T} vs dN_{ch}/d#eta/<dN_{ch}/d#eta>;dN_{ch}/d#eta/<dN_{ch}/d#eta>;<p_{T}^{J/#psi}>/<p_{T}^{J/#psi}_{int}>",size-1,axis);
    delete axis;
    
    ptInt = result->GetValue("MeanPtJPsi",sres.Data());
    ptIntError = result->GetErrorStat("MeanPtJPsi",sres.Data());

//    delete b;
  }
  else
  {
    hmPt = static_cast<TH1D*>(hrmPt->Clone("hJPsiMeanPtVSdNchdEta"));
    hmPt->SetTitle("J/#psi mean p_{T} vs dN_{ch}/d#eta");
    hmPt->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    hmPt->GetYaxis()->SetTitle("<p_{T}^{J/#psi}>");
  }
  
  Double_t systMptInt[9] = {0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014}; //FIXME: find a way to give this as input
  Double_t systMptBin[9] = {0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014};
  
  Double_t systMptRel[9] = {0.002,0.001,0.001,0.002,0.002,0.002,0.002,0.004,0.004}; //signal extraction pPb
  
//  Double_t systMptRel[9] = {0.002,0.001,0.012,0.001,0.002,0.002,0.001,0.004,0.003}; //signal extraction Pbp
  
//  Double_t systMptRel[9] = {0.001,0.002,0.001,0.002,0.002,0.003,0.005,0.000,0.000}; //signal extraction pp|eta|<05
  
//  Double_t systMptRel[9] = {0.002,0.002,0.002,0.002,0.001,0.002,0.001,0.012,0.000}; //signal extraction pp|eta|<1
  
  for ( Int_t i = 1 ; i <= hrmPt->GetNbinsX() ; i++ )
  {
    Double_t pt = hrmPt->GetBinContent(i);
    Double_t ptError = hrmPt->GetBinError(i);
    
    if ( relative )
    {
      ptError = TMath::Sqrt(TMath::Power(ptError/ptInt,2.) + TMath::Power((pt*ptIntError)/TMath::Power(ptInt,2.),2.));
      
      Double_t sMptInt = ptInt*systMptInt[i-1];
      Double_t sMptBin = pt*systMptBin[i-1];
      Double_t sysMptRel = TMath::Sqrt( TMath::Power(sMptBin/ptInt,2) + TMath::Power(pt*sMptInt/TMath::Power(ptInt,2.),2.) );
      
      pt /= ptInt;
      
      std::cout << TMath::Sqrt( TMath::Power(sysMptRel/pt,2.) +TMath::Power(systMptRel[i-1],2.) ) << std::endl;
      
      std::cout << pt << " +- " << ptError << std::endl;

    }
    
    hmPt->SetBinContent(i,pt);
    hmPt->SetBinError(i,ptError);
  }
  
  o = fMergeableCollection->Histo(Form("/RESULTS/%s",path.Data()),hmPt->GetName());
  
  if (o)
  {
    AliWarning(Form("Replacing /RESULTS/%s/%s",path.Data(),hmPt->GetName()));
    fMergeableCollection->Remove(Form("/RESULTS/%s/%s",path.Data(),hmPt->GetName()));
  }
  
  adoptOK = fMergeableCollection->Adopt(Form("/RESULTS/%s",path.Data()),hmPt);
  
  if ( adoptOK ) std::cout << "+++Mean Pt histo " << hmPt->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt mean pt histo %s",hmPt->GetName()));
  
  
  
  delete hrmPt;
  
  
  return;

  
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMu::ErrorPropagationAxBoverCxD(Double_t a,Double_t b,Double_t c,Double_t d)
{
  //Just valid for counts
  Double_t error2 = TMath::Power(b/(c*d),2.)*a + TMath::Power(a/(c*d),2.)*b + TMath::Power(a*b*d,2.)*(c/TMath::Power(c*d,4.)) + TMath::Power(a*b*c,2.)*(d/TMath::Power(c*d,4.));
  
  return TMath::Sqrt(error2);
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::ComputeEquNofMB(const char* what,const char* quantity,const char* flavour,Bool_t printout)
{
  
  AliAnalysisMuMuBinning* binning = BIN()->Project(what,quantity,flavour);
  TObjArray* dNchdEtas = binning->CreateBinObjArray();
  
  Double_t* binArray = binning->CreateBinArray();
  
  TIter next(dNchdEtas);
  AliAnalysisMuMuBinning::Range* r;
  
  TH1* hFNorm = ComputeDiffFnormFromHistos(what,quantity,flavour,kFALSE);
  
  TH1* hNMB = new TH1F("hNofEqMB","Equivalent number of MB triggers vs dN_{ch}/d#eta;dN_{ch}/d#eta;FNorm",dNchdEtas->GetEntries(),binArray);
  
  Int_t bin(0);
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    
    TH1* hCMUL = OC()->Histo(Form("/PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00/CMUL7-B-NOPF-MUON/V0A/%s",
                                  Form("EventsIn%s",r->AsString().Data())));
    if ( !hCMUL )
    {
      AliError(Form("No event histo in bin %s found for CMUL7-B-NOPF-MUON",r->AsString().Data()));
      return 0x0;
    }
    
    Double_t NMB = hCMUL->GetBinContent(1)*hFNorm->GetBinContent(++bin);
    Double_t NMBError = TMath::Sqrt(TMath::Power(hCMUL->GetBinContent(1)*hFNorm->GetBinError(bin),2.) + TMath::Power(TMath::Sqrt(hCMUL->GetBinContent(1))*hFNorm->GetBinContent(bin),2));
    
    if ( printout ) std::cout << r->AsString().Data() << " : " << NMB << " +- " << NMBError << std::endl;
    
    hNMB->SetBinContent(bin,NMB);
    hNMB->SetBinError(bin,NMBError);
  }
  
  delete dNchdEtas;
  delete[] binArray;
  
  return hNMB;
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

  const char* accEffSubResultName="PSICOUNT:1";
  
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
  
  const char* accEffSubResultName="PSICOUNT:1";
  
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
  AliAnalysisMuMuJpsiResult* r;
  
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
  {
    r = static_cast<AliAnalysisMuMuJpsiResult*>(realSpectra->BinContentArray()->At(i));
   
    StdoutToAliDebug(1,std::cout << "bin=";r->Print(););
    
    AliAnalysisMuMuJpsiResult* rsim = static_cast<AliAnalysisMuMuJpsiResult*>(simSpectra->BinContentArray()->At(i));
    
    Double_t mbeq = nofCINT7w0MUL / ( nofCINT7 * nofCMUL7);
    Double_t mbeqError = mbeq * AliAnalysisMuMuResult::ErrorABC( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
                                                                nofCINT7,TMath::Sqrt(nofCINT7),
                                                                nofCMUL7,TMath::Sqrt(nofCMUL7));
    
    r->Set("Fnorm",nofCINT7/nofCINT7w0MUL,(nofCINT7/nofCINT7w0MUL)*AliAnalysisMuMuResult::ErrorAB( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
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

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::RABy(const char* type, const char* direction)
{
  /// Compute the RAB...
  
  if (!SIM()) return 0x0;
  
  Double_t rapidityShift = 0.465;// 0.5*TMath::Log(208.0/82.0);
  const Double_t sqrts=5.023;
  const Double_t ymax=TMath::Log(sqrts*1000.0/3.096916);
  const Double_t tab = 0.093e-6; // nb^-1
  const Double_t tabError = 0.0035E-6; // nb^-1
  const char* accEffSubResultName="PSICOUNT:1";
  
  TF1 ydist("ydist","[0]*TMath::Exp(-(x*x)/(2.0*0.39*0.39))",0.,0.5);
  ydist.SetParameter(0,1.);

  //Normalization to the values presented by Zaida and Rosana on January 11th 2013 https://indico.cern.ch/conferenceDisplay.py?confId=224985 slide 22
  // Normalization is done in the rapidity range 2.75<y<3.25 where Rosanas values is 230.8+212.1
  Double_t y1_norma= 2.75/ymax;
  Double_t y2_norma= 3.25/ymax;
  Double_t normalization = 0.25*(230.8+212.1)/ydist.Integral(y1_norma, y2_norma);
  ydist.SetParameter(0,normalization);
//  AliInfoClass(Form("ymax=%e normalization=%f",ymax,ydist.Integral(y1_norma, y2_norma)));
  
  AliAnalysisMuMuSpectra* realSpectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/PSALL/CMUL7-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-%s",type)));
  AliAnalysisMuMuSpectra* simSpectra = static_cast<AliAnalysisMuMuSpectra*>(SIM()->OC()->GetObject(Form("/ALL/CMULLO-B-NOPF-MUON/PP/pMATCHLOWRABSBOTH/PSI-%s",type)));
  
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
  
  Double_t nofCMUL7 = CC()->GetSum(Form("trigger:CMUL7-B-NOPF-MUON/event:PSALL"));
  Double_t nofCINT7 = CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD/event:PSALL"));
  Double_t nofCINT7w0MUL = CC()->GetSum(Form("trigger:CINT7-B-NOPF-ALLNOTRD&0MUL/event:PSALL"));
  
  AliAnalysisMuMuBinning* binning = realSpectra->Binning();
  TObjArray* bins = binning->CreateBinObjArray();
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* bin;
  Int_t i(0);
  AliAnalysisMuMuJpsiResult* r;
  
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
    
    if ( bin->IsIntegrated() )
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
    
    r = static_cast<AliAnalysisMuMuJpsiResult*>(corrSpectra->BinContentArray()->At(i)->Clone());

    AliAnalysisMuMuJpsiResult* rsim = static_cast<AliAnalysisMuMuJpsiResult*>(simSpectra->BinContentArray()->At(i));
    
    Double_t mbeq = nofCINT7w0MUL / ( nofCINT7 * nofCMUL7);
    Double_t mbeqError = mbeq * AliAnalysisMuMuResult::ErrorABC( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
                                         nofCINT7,TMath::Sqrt(nofCINT7),
                                         nofCMUL7,TMath::Sqrt(nofCMUL7));
    
    r->Set("Fnorm",nofCINT7/nofCINT7w0MUL,(nofCINT7/nofCINT7w0MUL)*AliAnalysisMuMuResult::ErrorAB( nofCINT7w0MUL, TMath::Sqrt(nofCINT7w0MUL),
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
    
    AliAnalysisMuMuBinning::Range* bincm = new AliAnalysisMuMuBinning::Range(bin->What(),bin->Quantity(),ylowcms,yhighcms);
    
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
    
    r = static_cast<AliAnalysisMuMuJpsiResult*>(finalResults.At(j));

    bin = static_cast<AliAnalysisMuMuBinning::Range*>(finalBins.At(j));
    
    spectra->AdoptResult(*bin,r);
  }
  

  delete corrSpectra;
  
  return spectra;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SetConfig(const AliAnalysisMuMuConfig& config)
{
  /// (re)set the config
  delete fConfig;
  fConfig = new AliAnalysisMuMuConfig(config);
}

