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

/* $Id$ */

/// \ingroup macros
/// \file MUONStatusMap.C
/// \brief Macro to check/test pad status and pad status map makers
///
/// \author Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliMpCDB.h"
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "Riostream.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraph.h"
#include "TBox.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TText.h"
#include <vector>
#endif

namespace
{
  Int_t NTOTALNUMBEROFPADS(1064008);
}

//______________________________________________________________________________
void ReadIntegers(const char* filename, std::vector<int>& integers)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
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
      integers.push_back(s->String().Atoi());
    }
  }
  else
  {
    integers.push_back(sline.Atoi());
    
    while ( in >> i )
    {
      integers.push_back(i);
    }
  }
  
  std::sort(integers.begin(),integers.end());
}


//______________________________________________________________________________
void MUONStatusMap(AliMUONVStore*& vstatus,
                   AliMUONVStore*& vstatusMap,
                   const char* cdbStorage = "alien://folder=/alice/data/2011/OCDB",
                   Int_t runNumber=145292)
{  

  AliCDBManager::Instance()->SetDefaultStorage(cdbStorage);
  AliCDBManager::Instance()->SetRun(runNumber);
  
  AliMUONCDB::LoadMapping();
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  
  AliMUONCalibrationData cd(runNumber);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
  statusMaker.SetLimits(*recoParam);
  
  UInt_t mask = recoParam->PadGoodnessMask();

  statusMaker.Report(mask);
  
  vstatus = static_cast<AliMUONVStore*>(statusMaker.StatusStore()->Clone());
  
  const Bool_t deferredInitialization = kFALSE;
  
  AliMUONPadStatusMapMaker statusMapMaker(cd,mask,deferredInitialization);
    
  vstatusMap = static_cast<AliMUONVStore*>(statusMapMaker.StatusMap()->Clone());
}

//______________________________________________________________________________
Int_t GetBadChannels(Int_t runNumber,
                     Int_t& nbadped,
                     Int_t& nbadhv,
                     Int_t& nbadgain,
                     Int_t& nbadocc,
                     Int_t& nmissing,
                     Int_t& nreco)
{
  if (!AliCDBManager::Instance()->IsDefaultStorageSet())
  {
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }
  
  AliCDBManager::Instance()->SetRun(runNumber);
  
  AliMpCDB::LoadDDLStore();
  
  AliMUONCalibrationData cd(runNumber,true);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();

  statusMaker.SetLimits(*recoParam);
  
  AliMpManuIterator it;
  Int_t detElemId, manuId;
  
  Int_t pedCheck = ( 
                     AliMUONPadStatusMaker::kPedMeanZero |
                     AliMUONPadStatusMaker::kPedMeanTooLow |
                     AliMUONPadStatusMaker::kPedMeanTooHigh |
                     AliMUONPadStatusMaker::kPedSigmaTooLow |
                     AliMUONPadStatusMaker::kPedSigmaTooHigh );
     
  Int_t hvCheck = ( 
                   AliMUONPadStatusMaker::kHVError |
                   AliMUONPadStatusMaker::kHVTooLow |
                   AliMUONPadStatusMaker::kHVTooHigh |
                   AliMUONPadStatusMaker::kHVChannelOFF |
                   AliMUONPadStatusMaker::kHVSwitchOFF );

  
  Int_t occCheck = (
                    AliMUONPadStatusMaker::kManuOccupancyTooHigh 
                   );
                    
  Int_t ntotal(0);
  Int_t nbad(0);
  nbadped=0;
  nbadocc=0;
  nbadgain=0;
  nbadhv=0;
  nmissing=0;
  nreco=0;
  
  while ( it.Next(detElemId,manuId) )
  {
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
    {
      if ( de->IsConnectedChannel(manuId,manuChannel) )
      {
        ++ntotal;
        
        UInt_t status = statusMaker.PadStatus(detElemId, manuId, manuChannel);
        
        if (!status) continue;
        
        bool bad(false);
        
        if ( status & AliMUONPadStatusMaker::BuildStatus(pedCheck,0,0,0) ) 
        {
          ++nbadped;
          bad=true;
        }
        
        if ( status & AliMUONPadStatusMaker::BuildStatus(0,hvCheck,0,0) ) 
        {
          ++nbadhv;
          bad=true;
        }
        
        if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,0,occCheck) ) 
        {
          ++nbadocc;
          bad=true;
        }
        
        if ( status & recoParam->PadGoodnessMask() ) 
        {
          ++nreco;
        }
        
        if ( status & AliMUONPadStatusMaker::BuildStatus(AliMUONPadStatusMaker::kMissing,0,0,AliMUONPadStatusMaker::kMissing) )
        {
          bad=true;
          ++nmissing;
        }
        
        if (bad) ++nbad;
      }
    }
  }
  
  if (ntotal!=NTOTALNUMBEROFPADS)
  {
    cerr << Form("ERROR ! NOT THE EXPECTED NUMBER OF CHANNELS (%d vs 1064008) FOR RUN %09d",
                 ntotal,runNumber) << endl;
  }
  
  AliCDBManager::Instance()->ClearCache(); 
 
  return nbad;
}

//______________________________________________________________________________
void Draw(TFile* f, const char* gname, TLegend* l, Bool_t normalized)
{
  if (!f) return;
  
  TGraph* g = static_cast<TGraph*>(f->Get(gname));
  
  if (!g) return;
  
  if ( normalized ) 
  {
    g = static_cast<TGraph*>(g->Clone());
    for ( Int_t i = 0; i < g->GetN(); ++i ) 
    {
      Double_t y = g->GetY()[i];
      g->SetPoint(i,g->GetX()[i],y/NTOTALNUMBEROFPADS);
    }
  }
  
  g->Draw("lp");
  g->GetXaxis()->SetNdivisions(505);
  g->GetXaxis()->SetNoExponent();
  
  if (l) l->AddEntry(g,gname,"LP");
}

//______________________________________________________________________________
void DrawPeriod(double run1, double run2, double ymin, double ymax, const char* label)
{
  TBox* b = new TBox(run1,ymin,run2,ymax);
  b->SetFillColor(5);
  b->Draw();
  TText* text = new TText((run1+run2)/2.0,ymax*0.6,label);
  text->SetTextSize(0.05);
  text->Draw();
}

//______________________________________________________________________________
void DrawEvolution(const char* file, bool normalized=true)
{

  TFile* f = TFile::Open(file);
  
  if (!f) return;
  
  TCanvas* c = new TCanvas("mch-status-evolution","mch-status-evolution");
  
  c->SetGridy();
  c->SetTicky();
  
  c->Draw();
  
  TLegend* l = new TLegend(0.1,0.7,0.3,0.95,"ch evolution");

  TGraph* g = static_cast<TGraph*>(f->Get("nbad"));
  if (!g) return;
  
  int runmin = TMath::Nint(g->GetX()[0]);
  int runmax = TMath::Nint(g->GetX()[g->GetN()-1]);
  
  double ymax(0.4);
  
  TH2* h = new TH2F("hframe","hframe;Run number;Fraction of dead channels",100,runmin-200,runmax+200,100,0,ymax);

  gStyle->SetOptStat(kFALSE);
  h->Draw();
  h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetNdivisions(505);

  gStyle->SetOptTitle(kFALSE);

  DrawPeriod(115881,117222,0,ymax,"10b"); 

  DrawPeriod(119159,120824,0,ymax,"10c");

  DrawPeriod(122374,126424,0,ymax,"10d");

  DrawPeriod(127724,130850,0,ymax,"10e");

  DrawPeriod(133005,134929,0,ymax,"10f");

  DrawPeriod(135658,136376,0,ymax,"10g");

  DrawPeriod(137133,139513,0,ymax,"10h");
  
  DrawPeriod(143856,146860,0,ymax,"11a");

  DrawPeriod(148370,150702,0,ymax,"11b");

  DrawPeriod(151566,154583,0,ymax,"11c");

  DrawPeriod(158084,159606,0,ymax,"11d");

  DrawPeriod(160677,162000,0,ymax,"11e");

  Draw(f,"nbad",l,normalized);
  Draw(f,"nbadped",l,normalized);
  Draw(f,"nbadocc",l,normalized);
  Draw(f,"nbadhv",l,normalized);
  Draw(f,"nmissing",l,normalized);
  Draw(f,"nreco",l,normalized);
  
  h->Draw("same");

  c->RedrawAxis("g");
  
  l->Draw();
}

//______________________________________________________________________________
void MUONStatusMapEvolution(const char* runlist, const char* outfile)
{
  // Compute the number of bad pads (because of bad ped, bad hv, bad occupancy
  // or missing in configuration)
  //
  // output a root file with the different graphs.
  //
  // output can be then plotted using the DrawEvolution function
  //
  // Note that the output of different runlists can then be merged simply using
  // the hadd program, and so then the DrawEvolution can be used over
  // a huge period, e.g. a full year, while this method is better restricted
  // to a period or even less (depending on your success of accessing the OCDB)
  //
  
  std::vector<int> runs;

  ReadIntegers(runlist,runs);
  
  if ( runs.empty() ) 
  {
    cout << "No runs to process..." << endl;
    return;    
  }
  
  int year(2011);
  
  if ( runs[0] <= 139699 ) year=2010;
  
  AliCDBManager::Instance()->SetDefaultStorage(Form("alien://folder=/alice/data/%d/OCDB?cacheFold=/local/cdb",year));
  
  TList glist;
  
  glist.SetOwner(kTRUE);
  
  TGraph* gnbad = new TGraph(runs.size());
  gnbad->SetName("nbad");
  glist.Add(gnbad);
  
  TGraph* gnbadped = new TGraph(runs.size());
  gnbadped->SetName("nbadped");
  glist.Add(gnbadped);
  
  TGraph* gnbadocc = new TGraph(runs.size());
  gnbadocc->SetName("nbadocc");
  glist.Add(gnbadocc);
  
  TGraph* gnbadhv = new TGraph(runs.size());
  gnbadhv->SetName("nbadhv");
  glist.Add(gnbadhv);
  
  TGraph* gnmissing = new TGraph(runs.size());
  gnmissing->SetName("nmissing");
  glist.Add(gnmissing);

  TGraph* gnreco = new TGraph(runs.size());
  gnreco->SetName("nreco");
  glist.Add(gnreco);
  
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i ) 
  {
    Int_t runNumber = runs[i];
    Int_t nbadped;
    Int_t nbadhv;
    Int_t nbadgain;
    Int_t nbadocc;
    Int_t nmissing;
    Int_t nreco;
    
    Int_t nbad = GetBadChannels(runNumber,nbadped,nbadhv,nbadgain,nbadocc,nmissing,nreco);
    
    gnbad->SetPoint(i,runNumber,nbad);
    gnbadped->SetPoint(i,runNumber,nbadped);
    gnbadhv->SetPoint(i,runNumber,nbadhv);
    gnbadocc->SetPoint(i,runNumber,nbadocc);
    gnmissing->SetPoint(i,runNumber,nmissing);
    gnreco->SetPoint(i,runNumber,nreco);
  }
  
  TIter next(&glist);
  TGraph* g;
  Int_t index(0);

  TFile f(outfile,"recreate");
  Int_t  color[] = {  1 ,  2 ,  3 ,  4,  6, 8 };
  Int_t marker[] = { 28 , 24 , 23 , 26, 30, 5 };

  while ( ( g = static_cast<TGraph*>(next() ) ) )
  {
    g->GetXaxis()->SetNdivisions(505);
    g->GetXaxis()->SetNoExponent();
    g->SetMinimum(0);
    g->GetYaxis()->SetNoExponent();
    g->SetLineColor(color[index]);
    g->SetMarkerStyle(marker[index]);
    g->SetMarkerColor(color[index]);
    g->SetMarkerSize(1.0);
    ++index;
    g->Write();
  }
    
  f.Close();
  
  
}
