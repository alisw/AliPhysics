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

///\cond CLASSIMP
ClassImp(AliMUONTrackerHV)
///\endcond

//______________________________________________________________________________
AliMUONTrackerHV::AliMUONTrackerHV(const char* runlist, const char* ocdbPath)
 : AliMUONTrackerVoltages(runlist,ocdbPath)
{
  // ctor from a runlist (txt file)
  fOCDBObjectPath = "MUON/Calib/HV";
}

//______________________________________________________________________________
AliMUONTrackerHV::AliMUONTrackerHV(Int_t runNumber, const char* ocdbPath)
: AliMUONTrackerVoltages(runNumber,ocdbPath)
{
  // ctor for a single run
  fOCDBObjectPath = "MUON/Calib/HV";
}

//______________________________________________________________________________
AliMUONTrackerHV::~AliMUONTrackerHV()
{
  // dtor
}

//______________________________________________________________________________
TMap*
AliMUONTrackerHV::CreateMap(Int_t runNumber, Bool_t patched) const
{
  if (!patched) {
    return AliMUONTrackerVoltages::CreateMap(runNumber,kFALSE);
  }

  return AliMUONCalibrationData::CreateHV(runNumber,0x0,kTRUE,0x0,kTRUE);
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
          results.insert(std::make_pair(run,message));

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

  results.insert(std::make_pair(run,message));

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
      channel += Form("(%4d)",DCSNamer()->DetElemIdFromDCSAlias(alias.Data()));
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
  c3->SetBottomMargin(0.55);
  TH1* perchannel = hvoff->ProjectionY("hvoffperchannel");
  perchannel->GetXaxis()->SetBit(TAxis::kLabelsVert);
  perchannel->GetXaxis()->LabelsOption(">");
  perchannel->Draw("texthist");
  c3->Print(Form("%s-perchannel.pdf",outputBaseName));
}

//______________________________________________________________________________
void
AliMUONTrackerHV::Plot(const char* dcsname, Bool_t withPatch, Bool_t plotIntermediate)
{
  /// Show HV values for a given dcs name (or all if dcsname=0)
  /// Each canvas for each run will go to a separate PDF file

  AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
  TList messages;
  messages.SetOwner(kTRUE);
  TObjArray graphs;

  for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
  {
    Int_t runNumber = fRunList[i];

    messages.Delete();

    AliCDBManager::Instance()->SetRun(runNumber);

    TMap* m = AliMUONCalibrationData::CreateHV(runNumber,0x0,withPatch,&messages,kTRUE);

    TMultiGraph* mg = Map2Graph(m,dcsname);

    if ( !mg ) continue;

    graphs.Add(mg);

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

  std::map<std::string,int> channels;

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
              channels[channelName.Data()]++;
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

  AliInfo("--------------------------------------------------------------------");
  AliInfo("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

  int totalTrips(0);
  AliMUON2DMap tripMap(kTRUE);
  Int_t nofChannels(AliMpConstants::ManuNofChannels());

  for ( std::map<std::string,int>::const_iterator it = channels.begin(); it != channels.end(); ++it )
  {
    AliInfo(Form("%40s %3d",DCSNamer()->DCSNameFromAlias(it->first.c_str()).Data(),it->second));
    totalTrips += it->second;

    Int_t detElemId = DCSNamer()->DetElemIdFromDCSAlias(it->first.c_str());

    AliMpDetElement* de = AliMpDEStore::Instance()->GetDetElement(detElemId);

    // build the list of manuIds for this channel
    AliMpArrayI manuArray;

    manuArray.SetSize(300);

    Int_t index = DCSNamer()->DCSIndexFromDCSAlias(it->first.c_str());
    Int_t firstIndex(index);
    Int_t lastIndex(index);

    if ( index < 0 )
    {
      // it's a slat, must loop over PCBs
      firstIndex = 0;
      lastIndex = DCSNamer()->NumberOfPCBs(detElemId)-1;
    }

    for ( int i = firstIndex; i <= lastIndex ; ++i )
    {
      const AliMpArrayI* ma = de->ManusForHV(i);
      if (!ma)
      {
        AliError(Form("Could not get ma for de %d index %d",detElemId,i));
        continue;
      }
      for ( int j = 0; j < ma->GetSize(); ++j )
      {
        manuArray.Add(ma->GetValue(j),kFALSE);
      }
    }

    for ( Int_t iManu = 0; iManu < manuArray.GetSize(); ++iManu )
    {
      Int_t manuId = manuArray.GetValue(iManu);

      AliMUONVCalibParam* tripRate = new AliMUONCalibParamND(1,nofChannels,detElemId,manuId,0);

      tripMap.Add(tripRate);

      for ( Int_t j = 0 ; j < nofChannels; ++j )
      {
        tripRate->SetValueAsDouble(j,0,it->second*1.0);
      }
    }
  }

  AliInfo("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
  AliInfo(Form("Total of %3d trips for %4ld runs",totalTrips,fRunList.size()));

  AliMUONTrackerData* data = new AliMUONTrackerData("tripcount","Number of trips",1);
  data->Add(tripMap);
  data->SetDimensionName(0,"ntrips");

  AliMUONVTrackerDataMaker* dw = new AliMUONTrackerDataWrapper(data);

  AliMUONPainterDataRegistry::Instance()->Register(dw);

}

//______________________________________________________________________________
Int_t AliMUONTrackerHV::CompareMaps(const TMap& hv1, const TMap& hv2, Bool_t verbose) const
{
  /// Compare two HV maps (only HV voltages for the moment)
  /// Return the number of HV channels for which there is a difference

  Int_t ndiff(0);
  TIter next(&hv1);
  TObjString* hvChannelName;

  while ( ( hvChannelName = static_cast<TObjString*>(next()) ) )
  {
    TString name(hvChannelName->String());

    if ( name.Contains("sw") ) continue; // skip switches for the moment
    if ( name.Contains("iMon") ) continue; // skip HV currents for the moment

    Bool_t st1Check(kFALSE);

    if ( name.Contains("Chamber00Left") )
    {
      if (name.Contains("Quad1Sect0")) st1Check=kTRUE;
      if (name.Contains("Quad1Sect1")) st1Check=kTRUE;
      if (name.Contains("Quad1Sect2")) st1Check=kTRUE;
      if (name.Contains("Quad2Sect2")) st1Check=kTRUE;
      if (name.Contains("Quad2Sect1")) st1Check=kTRUE;
      if (name.Contains("Quad2Sect0")) st1Check=kTRUE;
    }
    else if ( name.Contains("Chamber01Left"))
    {
      if (name.Contains("Quad2Sect2")) st1Check=kTRUE;
      if (name.Contains("Quad2Sect0")) st1Check=kTRUE;
    }

    TPair* hvPair1 = static_cast<TPair*>(hv1.FindObject(name.Data()));
    TObjArray* values1 = static_cast<TObjArray*>(hvPair1->Value());

    TPair* hvPair2 = static_cast<TPair*>(hv2.FindObject(name.Data()));
    TObjArray* values2 = static_cast<TObjArray*>(hvPair2->Value());

    Bool_t same(kTRUE);

    if ( values1->GetEntries() != values2->GetEntries() )
    {
      same = kFALSE;
    }
    else
    {
      Int_t n = values1->GetEntries();
      for ( Int_t i = 0; i < n && same == kTRUE; ++i )
      {
        AliDCSValue* v1 = dynamic_cast<AliDCSValue*>(values1->At(i));
        AliDCSValue* v2 = dynamic_cast<AliDCSValue*>(values2->At(i));

        if ( v1->Compare(v2) != 0 || v1->GetType() != v2->GetType() || !IsAlmostEqualRelative(v1->GetFloat(),v2->GetFloat() ))
        {
          same = kFALSE;
        }
      }
    }

    if (!same)
    {
      ++ndiff;
    }

    if ( verbose && !same && st1Check )
    {
      std::cout << name << std::endl;
      values1->Print();
      values2->Print();
    }
  }

  return ndiff;
}
