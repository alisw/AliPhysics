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

#include "AliMUONTrackerHV.h"

#include <algorithm>
#include <map>
#include <set>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliDCSValue.h"
#include "AliGRPObject.h"
#include "AliMpDCSNamer.h"
#include "AliMpDEStore.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONCDB.h"
#include "AliLog.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLine.h"
#include "TMap.h"
#include "TMultiGraph.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStyle.h"
#include "Riostream.h"

//
// Class to inspect the MUON TRACKER HV values
//
// With this class you can :
//
// a) get a list of trips (method ReportTrips)
// b) print the values for some (or all) HV channels (method Print)
// c) plot the values for some (or all) HV channels (method Plot)
// d) get a list of HV channels that are "OFF" (methods Scan and HVoff)
//
// Note that in this class, all the output (either text or canvas) or the
// channel *names* used are the same as in the DCS UI at Pt2
// Specifically the chamber ids start at 1, the slat numbers at 1 and
// the quad and sect number at 1 also. And not at zero like for the
// DCS *aliases*. On the contraty, the internal map, coming from the OCDB,
// only contains aliases, not names. Confusing ? It is.
//

ClassImp(AliMUONTrackerHV)

//______________________________________________________________________________
AliMUONTrackerHV::AliMUONTrackerHV(const char* runlist, const char* ocdbPath)
: TObject(), fRunList(), fOCDBPath(ocdbPath), fDCSNamer(0x0)
{
  // ctor from a runlist (txt file)
  SetRunList(runlist);
}

//______________________________________________________________________________
AliMUONTrackerHV::AliMUONTrackerHV(Int_t runNumber, const char* ocdbPath)
: TObject(), fRunList(), fOCDBPath(ocdbPath), fDCSNamer(0x0)
{
  // ctor for a single run
  SetRunList(runNumber);
}

//______________________________________________________________________________
AliMUONTrackerHV::~AliMUONTrackerHV()
{
  // dtor
  delete fDCSNamer;
}

//______________________________________________________________________________
void AliMUONTrackerHV::ReadIntegers(const char* filename, std::vector<int>& integers)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  std::ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
  std::set<int> runset;
  
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
  
  for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
  {
    integers.push_back((*it));
  }
  
  std::sort(integers.begin(),integers.end());
}

//______________________________________________________________________________
AliMpDCSNamer*
AliMUONTrackerHV::DCSNamer() const
{
  // return the dcs namer
  if (!fDCSNamer)
  {
    if (!AliMpDEStore::Instance(false))
    {
      AliMUONCDB::LoadMapping();
    }
    fDCSNamer = new AliMpDCSNamer("TRACKER");
  }
  return fDCSNamer;
}

//______________________________________________________________________________
void AliMUONTrackerHV::SetRunList(Int_t runNumber)
{
  // Make the runlist be a single run
  fRunList.clear();
  fRunList.push_back(runNumber);
}

//______________________________________________________________________________
void
AliMUONTrackerHV::SetRunList(const char* runlist)
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
TGraph*
AliMUONTrackerHV::GraphValues(TMap* m, const char* dcsname)
{
  // make a graph of HV channels' voltage values for a given dcs name (name, not
  // alias)
  
  if ( TString(dcsname).Contains("sw") )
  {
    // do not graph switches
    return 0x0;
  }

  
  AliInfo(dcsname);
  
  TPair* p = static_cast<TPair*>(m->FindObject(DCSNamer()->DCSAliasFromName(dcsname).Data()));
  
  if (!p) return 0x0;
  
  TObjArray* a = static_cast<TObjArray*>(p->Value());
  TIter n2(a);
  AliDCSValue* val;
  Int_t i(0);

  TGraph* g = new TGraph(a->GetEntries());
  while ( ( val = static_cast<AliDCSValue*>(n2()) ) )
  {
    g->SetPoint(i,val->GetTimeStamp(),val->GetFloat());
    ++i;
  }
  return g;
}

//______________________________________________________________________________
void
AliMUONTrackerHV::Scan(Int_t verbose)
{
  /// Retrieve HV values from OCDB for a given run list, and check whether
  /// we have some issues with them...
  /// If you pipe the results of this into a text file, you can then
  /// feed it to the HVoff method for further investigations.
  ///
  
  if ( fRunList.empty() )
  {
    AliError("No runs to process...");
    return;    
  }
    
  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  
  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    AliMUONCDB::CheckHV(fRunList[i],verbose);
  }
}

//______________________________________________________________________________
void AliMUONTrackerHV::HVoff(const char* logfile, const char* outputBaseName)
{
  /// Check the number of HV which have problem
  /// the input is the output of e.g.
  /// .L MUONTrackerHV.C+
  /// ScanHV("lhc11de.list");>  lhc11de.log
  ///
  
  gStyle->SetOptStat(0);
  
  char line[1024];
  
  std::ifstream in(logfile);
  int run(-1),a,b,c,d,e,f,g,h,z,other;
  std::map<int,std::string> results;
  
  std::string message;
  const char* testProblem = "I-AliMUONCDB::CheckHV::CheckHV:      Problem at ";
  
  while ( in.getline(line,1023,'\n') )
  {
    TString sline(line);
    if (sline.Contains("SUMMARY"))
    {
      AliInfo(line);
      int r;
      sscanf(line,"I-AliMUONCDB::CheckHV::CheckHV: RUN %09d HVchannel SUMMARY : # of cases A(%3d) B(%3d) C(%3d) D(%3d) E(%3d) F(%3d) G(%3d) H(%3d) Z(%3d) OTHER(%3d)",
             &r,&a,&b,&c,&d,&e,&f,&g,&h,&z,&other);
      if ( r != run )
      {
        if ( run == -1 )
        {
          AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
          AliCDBManager::Instance()->SetRun(r);
          AliMUONCDB::LoadMapping();
        }
        
        if ( run > 0 )
        {
          results.insert(std::make_pair<int,std::string>(run,message));
          
        }
        message = "";
        run = r;
      }          
    }
    else if ( sline.Contains(testProblem) )
    {
      message += "|";
      message += sline(strlen(testProblem),sline.Length()-1).Data();
    }
  }
  
  results.insert(std::make_pair<int,std::string>(run,message));
  
  TH2* hvoff = new TH2I(outputBaseName,outputBaseName,1,0,1,1,0,1);
  
  std::map<int,std::string>::const_iterator it;
  
  for ( it = results.begin(); it != results.end(); ++it )
  {
    AliInfo(Form("%d -> %s",it->first,it->second.c_str()));
    TObjArray* split = TString(it->second.c_str()).Tokenize("|");
    TIter next(split);
    TObjString* str;
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      TString s(str->String());
      TObjArray* parts = s.Tokenize(":");
      TString alias = (static_cast<TObjString*>(parts->At(0)))->String();
      TString channel = DCSNamer()->DCSNameFromAlias(alias.Data());
      channel.ReplaceAll(".actual.vMon","");
      hvoff->Fill(Form("%6d",it->first),channel.Data(),1.0);
      delete parts;
    }
    delete split;
  }
  
  hvoff->LabelsDeflate("x");
  hvoff->LabelsDeflate("y");
  hvoff->LabelsOption("x","<");
  hvoff->LabelsOption("y","<");
  
  TCanvas* c1 = new TCanvas;
  c1->SetLeftMargin(0.35);
  hvoff->Draw("text");
  c1->Print(Form("%s.pdf",outputBaseName));
  TCanvas* c2 = new TCanvas;
  TH1* hx = hvoff->ProjectionX("hvoffperrun");
  hx->Draw();
  c2->Print(Form("%s-perrun.pdf",outputBaseName));
  TCanvas* c3 = new TCanvas;
  c3->SetBottomMargin(0.5);
  TH1* perchannel = hvoff->ProjectionY("hvoffperchannel");
  perchannel->GetXaxis()->SetBit(TAxis::kLabelsVert);
  perchannel->GetXaxis()->LabelsOption(">");
  perchannel->Draw();
  c3->Print(Form("%s-perchannel.pdf",outputBaseName));
}

//______________________________________________________________________________
void AliMUONTrackerHV::TimeAxis(TMultiGraph* g)
{
  g->GetXaxis()->SetTimeDisplay(1);
//  g->GetXaxis()->SetTimeFormat("%d/%m %H:%M%F2010-12-31 24:00:00");
  g->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
  g->GetXaxis()->SetTimeOffset(0,"gmt");
  g->GetXaxis()->SetNdivisions(505);
}

//______________________________________________________________________________
TMultiGraph*
AliMUONTrackerHV::GraphHV(TMap* m, const char* dcsname)
{
  // Make a graph of the values matching dcsname
  TIter next(m);
  TObjString* s;
  
  TMultiGraph* mg = new TMultiGraph;

  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TString name(DCSNamer()->DCSNameFromAlias(s->String()));
    
    if ( dcsname && !name.Contains(dcsname)) continue;
    
    TGraph* g = GraphValues(m,name);
    
    if ( g ) 
    {
      g->SetMarkerSize(1.5);
      g->SetMarkerStyle(2);
      g->SetLineStyle(2);
      mg->Add(g,"lp");
      g->SetTitle(name.Data());
    }
  }  

  return mg;
}

//______________________________________________________________________________
void
AliMUONTrackerHV::Print(Option_t* dcsname) const
{
  /// Print HV values for a given dcs name (or all if dcsname=0)
  
  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  TList messages;
  messages.SetOwner(kTRUE);
  
  for ( std::vector<int>::size_type iRun = 0; iRun < fRunList.size(); ++iRun )
  {
    Int_t runNumber = fRunList[iRun];
    
    AliInfo("---------------------");
    AliInfo(Form("RUN %09d",runNumber));
    
    messages.Delete();
    
    AliCDBManager::Instance()->SetRun(runNumber);
    
    TMap* m = AliMUONCalibrationData::CreateHV(runNumber,0x0,kFALSE,&messages,kTRUE);
    
    TIter next(m);
    TObjString* s;
    
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {      
      TString name(DCSNamer()->DCSNameFromAlias(s->String()));
      
      if ( dcsname && !name.Contains(dcsname)) continue;
      
      TPair* p = static_cast<TPair*>(m->FindObject(DCSNamer()->DCSAliasFromName(dcsname).Data()));
      
      if (!p) continue;
      
      TObjArray* a = static_cast<TObjArray*>(p->Value());
      TIter n2(a);
      AliDCSValue* val;
      Int_t i(0);
      
      while ( ( val = static_cast<AliDCSValue*>(n2()) ) )
      {
        std::cout << Form("i=%5d ",i) << std::endl;
        val->Print("");
        ++i;
      }
    }
  }
}

//______________________________________________________________________________
void
AliMUONTrackerHV::Plot(const char* dcsname, Bool_t withPatch)
{
  /// Show HV values for a given dcs name (or all if dcsname=0)
  /// Each canvas for each run will go to a separate PDF file
  
  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  TList messages;
  messages.SetOwner(kTRUE);
  
  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];
  
    messages.Delete();
    
    AliCDBManager::Instance()->SetRun(runNumber);
    
    TMap* m = AliMUONCalibrationData::CreateHV(runNumber,0x0,withPatch,&messages,kTRUE);
    
    TMultiGraph* mg = GraphHV(m,dcsname);
    
    if ( !mg ) continue;
    
    TString cname(Form("MCH_HV_RUN%09d",runNumber));
    
    if ( strlen(dcsname) > 0 )
    {
      TString s(dcsname);
      s.ReplaceAll("/","_");
      cname += Form("_dcsname_%s",s.Data());
    }
    
    AliCDBEntry* e = AliCDBManager::Instance()->Get("GRP/GRP/Data",runNumber);
    
    TLine* startRunLine(0);
    TLine* endRunLine(0);
    time_t start(0);
    time_t end(0);
    
    if ( e )
    {
      AliGRPObject* grp = static_cast<AliGRPObject*>(e->GetObject());
      if (grp)
      {
        start = grp->GetTimeStart();
        end = grp->GetTimeEnd();
      }
    }
    
    if ( end )
    {
      TGraph* g = new TGraph(1);
      g->SetPoint(0,end,0);
      mg->Add(g,"");
    }
    
    TCanvas* c = new TCanvas(cname.Data(),cname.Data());
    
    c->Draw();
    
    mg->SetTitle(cname.Data());
    
    mg->Draw("AL");
    
    TimeAxis(mg);
    
    if ( start )
    {
      startRunLine = new TLine(start,mg->GetYaxis()->GetXmin(),start,mg->GetYaxis()->GetXmax());
      startRunLine->SetLineColor(2);
      startRunLine->SetLineWidth(4);
    }
    if  ( end )
    {
      endRunLine = new TLine(end,mg->GetYaxis()->GetXmin(),end,mg->GetYaxis()->GetXmax());
      endRunLine->SetLineColor(2);
      endRunLine->SetLineWidth(4);
    }
    
    if ( startRunLine ) startRunLine->Draw();
    if ( endRunLine ) endRunLine->Draw();
    
    c->SaveAs(Form("%s.pdf",cname.Data()));
  }
}

//______________________________________________________________________________
void
AliMUONTrackerHV::ReportTrips(Bool_t includeLowOnes)
{
  /// Report trips
  /// if includeLowOnes is kTRUE we'll report also the trips which starts from non-operational voltage values
  
  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  
  TList messages;
  messages.SetOwner(kTRUE);
  TObjString* msg(0);

  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];
    
    AliInfo("---------------------");
    
    Int_t ntrips(0);
    
    messages.Delete();
    
    AliCDBManager::Instance()->SetRun(runNumber);
    
    AliMUONCalibrationData::CreateHV(runNumber,0x0,kTRUE,&messages,kTRUE);
    
    if (!AliMpDEStore::Instance(false))
    {
      AliMUONCDB::LoadMapping();
    }
    
    TIter next(&messages);

    while ( ( msg = static_cast<TObjString*>(next())) )
    {
      if ( msg->String().Contains("TRIP") && ( includeLowOnes || !msg->String().Contains("LOWTRIP") ) )
      {
        ++ntrips;
      }
    }

    AliInfo(Form("RUN %09d - %d trip%c",runNumber,ntrips,(ntrips>1 ? 's':' ')));
    
    next.Reset();
    std::map<int,std::string> report;
    
    while ( ( msg = static_cast<TObjString*>(next())) )
    {
      if ( msg->String().Contains("TRIP") )
      {
        TObjArray* parts = msg->String().Tokenize(" ");
        TString channelName(static_cast<TObjString*>(parts->At(0))->String());
        
        for ( Int_t ip = 0; ip <= parts->GetLast(); ++ip)
        {
          TString p(static_cast<TObjString*>(parts->At(ip))->String());
          
          if ( p.Contains("TRIP") )
          {
            if ( includeLowOnes || !p.Contains("LOWTRIP") )
            {
              TString ts(static_cast<TObjString*>(parts->At(ip+2))->String());
          
              ip += 3;
          
              Int_t index = ts.Index("TS:");
          
              UInt_t timeStamp = TString(ts(index+strlen("TS:"),ts.Length()-index)).Atoi();
          
              TString tmp(msg->String());
              tmp.ReplaceAll(channelName.Data(),DCSNamer()->DCSNameFromAlias(channelName.Data()));
              report[timeStamp] = tmp.Data();
            }
          }
        }
        delete parts;
      }
    }

    for ( std::map<int,std::string>::const_iterator it = report.begin(); it != report.end(); ++it )
    {
      AliInfo(Form("%s %s",TTimeStamp(it->first).AsString("s"),it->second.c_str()));
    }
  }
}

