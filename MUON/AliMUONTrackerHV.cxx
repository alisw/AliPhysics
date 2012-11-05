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
#include "Riostream.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliDCSValue.h"
#include "TMap.h"
#include <map>
#include "AliMpDCSNamer.h"
#include "TH2.h"
#include "TStyle.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TLine.h"
#include <set>
#include "AliLog.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "AliMUONCDB.h"
#include "TCanvas.h"
#include "AliMUONCalibrationData.h"
#include "AliGRPObject.h"

ClassImp(AliMUONTrackerHV)

//______________________________________________________________________________
AliMUONTrackerHV::AliMUONTrackerHV(const char* runlist, const char* ocdbPath) : TObject(), fRunList(), fOCDBPath(ocdbPath)
{
  // ctor from a runlist (txt file)
  SetRunList(runlist);
}

//______________________________________________________________________________
AliMUONTrackerHV::AliMUONTrackerHV(Int_t runNumber, const char* ocdbPath) : TObject(), fRunList(), fOCDBPath(ocdbPath)
{
  // ctor for a single run
  SetRunList(runNumber);
}

//______________________________________________________________________________
AliMUONTrackerHV::~AliMUONTrackerHV()
{
  // dtor
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
AliMUONTrackerHV::ShowValues(TMap* m, const char* name)
{
  // make a graph of HV channels' voltage values for a given dcs alias (name)
  
  TGraph* g(0x0);
  
  AliInfo(name);
  
  TPair* p = static_cast<TPair*>(m->FindObject(name));
  TObjArray* a = static_cast<TObjArray*>(p->Value());
  TIter n2(a);
  AliDCSValue* val;
  Int_t i(0);
  
  while ( ( val = static_cast<AliDCSValue*>(n2()) ) )
  {
    StdoutToAliInfo(std::cout << Form("i=%5d ",i);
                    val->Print(""););
    ++i;
  }
  
  if ( TString(name).Contains("sw") ) 
  {
    // do not graph switches
    return 0x0;    
  }
  
  n2.Reset();
  g = new TGraph(a->GetEntries());
  i = 0;
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
  
  if ( fRunList.empty() )
  {
    std::cout << "No runs to process..." << std::endl;
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
  
  AliMpDCSNamer hvNamer("TRACKER");
  
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
      TString channel = hvNamer.DCSChannelNameFromAlias(alias.Data());
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
AliMUONTrackerHV::ShowHV(TMap* m, const char* dcsname)
{
  TIter next(m);
  TObjString* s;
  AliMpDCSNamer hvNamer("TRACKER");
  TMultiGraph* mg = new TMultiGraph;

  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TString name(s->String());
    
    if ( dcsname && !name.Contains(dcsname)) continue;
    
    TGraph* g = ShowValues(m,name);
    
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
AliMUONTrackerHV::Plot(const char* dcsname, Bool_t withPatch)
{
  /// Show HV values for a given dcs alias (or all if dcsname=0)
  /// Each canvas for each run will go to a separate PDF file
  
  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  TList messages;
  messages.SetOwner(kTRUE);
  
  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];
  
    messages.Delete();
    
    AliCDBManager::Instance()->SetRun(runNumber);
    
    TMap* m = AliMUONCalibrationData::CreateHV(runNumber,0x0,withPatch,&messages);
    
    TMultiGraph* mg = ShowHV(m,dcsname);
    
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
AliMUONTrackerHV::ReportTrips()
{
  /// Report trips
  
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
    
    AliMUONCalibrationData::CreateHV(runNumber,0x0,kTRUE,&messages);
    
    TIter next(&messages);

    while ( ( msg = static_cast<TObjString*>(next())) )
    {
      if ( msg->String().Contains("TRIP") )
      {
        ++ntrips;
      }
    }

    AliInfo(Form("RUN %09d - %d trip%c",runNumber,ntrips,(ntrips>1 ? 's':' ')));
    
    next.Reset();
    
    while ( ( msg = static_cast<TObjString*>(next())) )
    {
      if ( msg->String().Contains("TRIP") )
      {
        AliInfo(msg->String().Data());        
      }
    }
  }
}

