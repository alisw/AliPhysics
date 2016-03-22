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

#include "AliMUONTrackerVoltages.h"

#include <algorithm>
#include <map>
#include <set>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliDCSValue.h"
#include "AliGRPObject.h"
#include "AliMpArrayI.h"
#include "AliMpConstants.h"
#include "AliMpDCSNamer.h"
#include "AliMpDEStore.h"
#include "AliMpDetElement.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONCDB.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerDataWrapper.h"
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
// Base class to inspect the MUON TRACKER HV and LV values
//
// With this class you can :
//
// - print the values for some (or all) HV channels or LV groups (method Print)
// - plot the values for some (or all) HV channels or LV groups (method Plot)
//
// Note that in this class, all the output (either text or canvas) or the
// channel *names* used are the same as in the DCS UI at Pt2
// Specifically the chamber ids start at 1, the slat numbers at 1 and
// the quad and sect number at 1 also. And not at zero like for the
// DCS *aliases*. On the contraty, the internal map, coming from the OCDB,
// only contains aliases, not names. Confusing ? It is.
//

ClassImp(AliMUONTrackerVoltages)

//______________________________________________________________________________
AliMUONTrackerVoltages::AliMUONTrackerVoltages(const char* runlist, const char* ocdbPath)
: TObject(), fRunList(), fOCDBPath(ocdbPath), fDCSNamer(0x0)
{
  // ctor from a runlist (txt file)
  SetRunList(runlist);
}

//______________________________________________________________________________
AliMUONTrackerVoltages::AliMUONTrackerVoltages(Int_t runNumber, const char* ocdbPath)
: TObject(), fRunList(), fOCDBPath(ocdbPath), fDCSNamer(0x0)
{
  // ctor for a single run
  SetRunList(runNumber);
}

//______________________________________________________________________________
AliMUONTrackerVoltages::~AliMUONTrackerVoltages()
{
  // dtor
  delete fDCSNamer;
}

//____________________________________________________________________________
TMultiGraph* AliMUONTrackerVoltages::CombineMulti(TObjArray& graphs)
{
  // combine multigraphs

  TMultiGraph* rv = new TMultiGraph;

  TIter next(&graphs);
  TMultiGraph* mg;
  TMultiGraph* ref = static_cast<TMultiGraph*>(next());

  Int_t dref = ref->GetListOfGraphs()->GetEntries();

  while ( ( mg = static_cast<TMultiGraph*>(next())) )
  {
    TList* list = mg->GetListOfGraphs();
    Int_t d1 = list->GetEntries();

    if (  d1 != dref )
    {
      AliError(Form("%d vs %d",d1,dref));
      return 0x0;
    }
  }

  for ( Int_t i = 0; i < dref; ++i )
  {
    TObjArray graph;
    next.Reset();
    while ( ( mg = static_cast<TMultiGraph*>(next())) )
    {
      graph.Add(mg->GetListOfGraphs()->At(i));
      TGraph* g = Combine(graph);
      rv->Add(g);
    }
  }
  return rv;
}

//____________________________________________________________________________
TGraph* AliMUONTrackerVoltages::Combine(TObjArray& graphs)
{
  // make one graph out of several
  // x axis is supposed to be time and will end up ordered in the
  // returned graph

  std::map<int, std::vector<double> > values;
  std::map<int, std::vector<double> >::const_iterator it;

  TIter next(&graphs);
  TGraph* g;

  while ( ( g = static_cast<TGraph*>(next())) )
  {
    for ( Int_t i = 0; i < g->GetN(); ++i )
    {
      std::vector<double> pair;

      pair.push_back(g->GetX()[i]);
      pair.push_back(g->GetY()[i]);

      values.insert( std::make_pair(g->GetX()[i],pair));
    }
  }

  TGraph* rv(0x0);

  if ( values.size() )
  {
    std::vector<double> vx;
    std::vector<double> vy;

    for ( it = values.begin(); it != values.end(); ++it )
    {
      const std::vector<double>& q = it->second;

      vx.push_back(q[0]);
      vy.push_back(q[1]);
    }

    rv = new TGraph(values.size(),&vx[0],&vy[0]);
    rv->GetXaxis()->SetNoExponent();

    g = static_cast<TGraph*>(graphs.At(0));

    rv->SetName(g->GetName());
    rv->SetTitle(g->GetTitle());
  }

  return rv;
}

//______________________________________________________________________________
TMap*
AliMUONTrackerVoltages::CreateMap(Int_t runNumber, Bool_t /*patched*/) const
{
  /// Create the map object
  return dynamic_cast<TMap*>(AliMUONCalibrationData::CreateObject(runNumber,fOCDBObjectPath));
}

//______________________________________________________________________________
bool AliMUONTrackerVoltages::IsAlmostEqualRelative(float A, float B, float maxRelDiff) const
{
    // Calculate the difference.
    float diff = fabs(A - B);
    A = fabs(A);
    B = fabs(B);
    // Find the largest
    float largest = (B > A) ? B : A;

    if (diff <= largest * maxRelDiff)
      return true;
    return false;
}

//______________________________________________________________________________
TMultiGraph*
AliMUONTrackerVoltages::Map2Graph(TMap* m, const char* dcsname)
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
void AliMUONTrackerVoltages::ReadIntegers(const char* filename, std::vector<int>& integers)
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
AliMUONTrackerVoltages::DCSNamer() const
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
void AliMUONTrackerVoltages::SetRunList(Int_t runNumber)
{
  // Make the runlist be a single run
  fRunList.clear();
  fRunList.push_back(runNumber);
}

//______________________________________________________________________________
void
AliMUONTrackerVoltages::SetRunList(const char* runlist)
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
AliMUONTrackerVoltages::GraphValues(TMap* m, const char* dcsname)
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

  std::cout << "a for " << p->Key() << " : " << std::endl;

  a->Print();

  TIter n2(a);
  AliDCSValue* val;
  Int_t i(0);

  TGraph* g = new TGraph(a->GetEntries());
  while ( ( val = static_cast<AliDCSValue*>(n2()) ) )
  {
    g->SetPoint(i,val->GetTimeStamp(),val->GetFloat());
    ++i;
  }
  g->SetName(dcsname);
  return g;
}

//______________________________________________________________________________
void AliMUONTrackerVoltages::TimeAxis(TMultiGraph* g)
{
  g->GetXaxis()->SetTimeDisplay(1);
//  g->GetXaxis()->SetTimeFormat("%d/%m %H:%M%F2010-12-31 24:00:00");
  g->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
  g->GetXaxis()->SetTimeOffset(0,"gmt");
  g->GetXaxis()->SetNdivisions(505);
}

//______________________________________________________________________________
void
AliMUONTrackerVoltages::Print(Option_t* dcsname) const
{
  /// Print HV values for a given dcs name (or all if dcsname=0)

  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());

  for ( std::vector<int>::size_type iRun = 0; iRun < fRunList.size(); ++iRun )
  {
    Int_t runNumber = fRunList[iRun];

    AliInfo("---------------------");
    AliInfo(Form("RUN %09d",runNumber));

    AliCDBManager::Instance()->SetRun(runNumber);

    TMap* m = CreateMap(runNumber);
    TIter next(m);
    TObjString* s;

    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      TString name(DCSNamer()->DCSNameFromAlias(s->String()));

      if ( dcsname && !name.Contains(dcsname)) continue;

      TPair* p = static_cast<TPair*>(m->FindObject(s->String()));

      if (!p) continue;

      std::cout << "==== " << s->String() << std::endl;

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
AliMUONTrackerVoltages::Plot(const char* dcsname, Bool_t withPatch, Bool_t plotIntermediate)
{
  /// Show HV values for a given dcs name (or all if dcsname=0)
  /// Each canvas for each run will go to a separate PDF file

  Ssiz_t ix = fOCDBObjectPath.Last('/');
  TString what = fOCDBObjectPath(ix+1,fOCDBObjectPath.Length()-ix);

  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  TObjArray graphs;

  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];

    AliCDBManager::Instance()->SetRun(runNumber);

    TMap* m = CreateMap(runNumber,withPatch);

    TMultiGraph* mg = Map2Graph(m,dcsname);

    if ( !mg ) continue;

    graphs.Add(mg);

    TString cname(Form("MCH_%s_RUN%09d",what.Data(),runNumber));

    if ( dcsname && strlen(dcsname) > 0 )
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
        TDatime dstart(start);
        TDatime dend(end);

        dstart.Print();
        dend.Print();
      }
    }

    if ( end )
    {
      TGraph* g = new TGraph(1);
      g->SetPoint(0,end,0);
      mg->Add(g,"");
    }

    if ( plotIntermediate )
    {
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

  new TCanvas;

  TMultiGraph* g = CombineMulti(graphs);

  TIter next(g->GetListOfGraphs());
  TGraph* gi;

  while ( ( gi = static_cast<TGraph*>(next())))
  {
    gi->SetMarkerStyle(kPlus);
  }
  g->Draw("alp");
  TimeAxis(g);
}
