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

#include "AliAnalysisMuMuConfig.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TColor.h"
#include "TStyle.h"
#include "TSystem.h"
#include "AliLog.h"
#include <set>
#include <cassert>
#include "TMap.h"
#include "THashList.h"
#include "TROOT.h"

ClassImp(AliAnalysisMuMuConfig)

namespace {

    Bool_t GetKeyValue(const TString& str, const char separator, TString& key, TString& value)
    {
    /// Get a key value pair, separated by separator character
    key=value="";
    if ( !str.CountChar(separator) ) return kFALSE;
    Int_t index = str.Index(separator);
    key = str(0,index);
    key.Remove(TString::kBoth,' ');
    value = str(index+1,str.Length()-index-1);
    return kTRUE;
    }

    void DecodeRunRanges(const TString& runranges, std::set<int>& runs)
    {
    /// From a string of form A-D,F,H,W-Y,Z return a set containing
    /// A,B,C,D,F,H,W,X,Y,Z
    TObjArray* ranges = runranges.Tokenize(",");
    TObjString* s;
    TIter next(ranges);

    while ( ( s = static_cast<TObjString*>(next())))
        {
        int first = s->String().Atoi();
        int last = first;

        if ( s->String().Contains("-") )
            {
            TString sfirst,slast;
            GetKeyValue(s->String(),'-',sfirst,slast);
            last = slast.Atoi();
            assert(first==sfirst.Atoi());
            }
        for (int i = first; i <= last; ++i )
            {
            runs.insert(i);
            }
        }
    delete ranges;
    }

    void PrintRunInfo(const std::map<AliAnalysisMuMuConfig::EPerRunInfo,std::string>& map)
    {
    if ( map.count(AliAnalysisMuMuConfig::kMBTriggerClassName) )
        {
        std::cout << "MB=" << map.find(AliAnalysisMuMuConfig::kMBTriggerClassName)->second << " ";
        }
    if ( map.count(AliAnalysisMuMuConfig::kMIXTriggerClassName) )
        {
        std::cout << "MIX=" << map.find(AliAnalysisMuMuConfig::kMIXTriggerClassName)->second << " ";
        }
    if ( map.count(AliAnalysisMuMuConfig::kMULTriggerClassName) )
        {
        std::cout << "MUL=" << map.find(AliAnalysisMuMuConfig::kMULTriggerClassName)->second << " ";
        }
    if ( map.count(AliAnalysisMuMuConfig::kMSLTriggerClassName) )
        {
        std::cout << "MSL=" << map.find(AliAnalysisMuMuConfig::kMSLTriggerClassName)->second << " ";
        }
    if ( map.count(AliAnalysisMuMuConfig::kMSHTriggerClassName) )
        {
        std::cout << "MSH=" << map.find(AliAnalysisMuMuConfig::kMSHTriggerClassName)->second << " ";
        }
    if ( map.count(AliAnalysisMuMuConfig::kCentralityName) )
        {
        std::cout << "CENT=" << map.find(AliAnalysisMuMuConfig::kCentralityName)->second << " ";
        }
    std::cout << std::endl;
    }

}

const char* AliAnalysisMuMuConfig::RefMixTriggerKey()       const { return "RefMixTrigger"; }
const char* AliAnalysisMuMuConfig::RefMixEventSelectionKey()const { return "RefMixEventSelection"; }
const char* AliAnalysisMuMuConfig::DimuonTriggerKey()       const { return "DimuonTrigger"; }
const char* AliAnalysisMuMuConfig::MuonTriggerKey()         const { return "MuonTrigger"; }
const char* AliAnalysisMuMuConfig::MinbiasTriggerKey()      const { return "MinbiasTrigger"; }
const char* AliAnalysisMuMuConfig::MixTriggerKey()          const { return "MixTrigger"; }
const char* AliAnalysisMuMuConfig::EventSelectionKey()      const { return "EventSelection"; }
const char* AliAnalysisMuMuConfig::EventSelectionMixKey()   const { return "EventSelectionMix"; }
const char* AliAnalysisMuMuConfig::PairSelectionKey()       const { return "PairSelection"; }
const char* AliAnalysisMuMuConfig::CentralitySelectionKey() const { return "Centrality"; }
const char* AliAnalysisMuMuConfig::FitTypeKey()             const { return "FitType"; }
const char* AliAnalysisMuMuConfig::FitSingleKey()           const { return "FitSingle"; }
const char* AliAnalysisMuMuConfig::CompactGraphKey()        const { return "CompactGraphs"; }
const char* AliAnalysisMuMuConfig::OCDBPathKey()            const { return "OCDBPath"; }

//_____________________________________________________________________________
AliAnalysisMuMuConfig::AliAnalysisMuMuConfig() : TObject(),
fMap(0x0),
fOCDBPath("raw://"),
fIsCompactGraphs(kFALSE),
fPerRunInfo()
{
    // ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuConfig::~AliAnalysisMuMuConfig()
{
    // dtor
    delete fMap;
}

//_____________________________________________________________________________
AliAnalysisMuMuConfig::AliAnalysisMuMuConfig(const AliAnalysisMuMuConfig& other)
{
    if ( &other != this )
    {
    fOCDBPath = other.fOCDBPath;
    fIsCompactGraphs = other.fIsCompactGraphs;
    if (other.fMap)
        {
        fMap = static_cast<TMap*>(other.fMap->Clone());
        }
    fPerRunInfo = other.fPerRunInfo;
    }
}

//_____________________________________________________________________________
AliAnalysisMuMuConfig& AliAnalysisMuMuConfig::operator=(const AliAnalysisMuMuConfig& other)
{
    if ( &other != this )
    {
    fOCDBPath = other.fOCDBPath;
    fIsCompactGraphs = other.fIsCompactGraphs;
    delete fMap;
    fMap = 0x0;
    if (other.fMap)
        {
        fMap = static_cast<TMap*>(other.fMap->Clone());
        }
    fPerRunInfo = other.fPerRunInfo;
    }
    return *this;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::Add(const char* key, const TString& line)
{
    /// Add a configuration element, in the form :
    /// key : line
    /// where line is "actual_value [sim] [real]

    UInt_t dataType(0);

    if ( !line.Contains("SIM",TString::kIgnoreCase) &&
        !line.Contains("REAL",TString::kIgnoreCase) )
    {
    // no REAL or SIM keyword = REAL
    dataType |= kReal;
    }
    else if ( line.Contains("SIM",TString::kIgnoreCase) )
    {
    dataType |= kSim;
    if ( line.Contains("REAL",TString::kIgnoreCase) )
        {
        dataType |= kReal;
        }
    }
    else if ( line.Contains("REAL",TString::kIgnoreCase) )
    {
    dataType |= kReal;
    }

    TObjArray* a = line.Tokenize(" ");
    TString es = static_cast<TObjString*>(a->First())->String();
    delete a;

    THashList* list = static_cast<THashList*>(Map()->GetValue(key));
    if (!list)
    {
    list = new THashList;
    list->SetOwner(kTRUE);
    Map()->Add(new TObjString(key),list);

    }

    TObjString* s = static_cast<TObjString*>(list->FindObject(es));
    if ( !s )
    {
    s = new TObjString(es);
    s->SetUniqueID(dataType);
    list->Add(s);
    }
    else
    {
    s->SetUniqueID(dataType);
    s->String() = es;
    }
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::Clear(Option_t*)
{
    /// Clear the internal map
    if (fMap) fMap->DeleteAll();
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::First(const char* key, Bool_t simulation) const
{
    /// Get the first element of a key

    THashList* list = static_cast<THashList*>(Map()->GetValue(key));
    TIter next(list);
    TObjString* s;

    UInt_t test(kReal);

    if ( simulation )
    {
    test = kSim;
    }

    while ( ( s = static_cast<TObjString*>(next())) )
    {
    if ( s->GetUniqueID() & test )
        {
        return s->String();
        }
    }
    return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetMBTriggerClassName(Int_t runNumber) const
{
    if ( fPerRunInfo.count(runNumber) )
    {
    return (fPerRunInfo.find(runNumber)->second).find(kMBTriggerClassName)->second;
    }
    return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetMIXTriggerClassName(Int_t runNumber) const
{
    if ( fPerRunInfo.count(runNumber) )
    {
    return (fPerRunInfo.find(runNumber)->second).find(kMIXTriggerClassName)->second;
    }
    return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetMULTriggerClassName(Int_t runNumber) const
{
    if ( fPerRunInfo.count(runNumber) )
    {
    return (fPerRunInfo.find(runNumber)->second).find(kMULTriggerClassName)->second;
    }
    return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetMSLTriggerClassName(Int_t runNumber) const
{
    if ( fPerRunInfo.count(runNumber) )
    {
    return (fPerRunInfo.find(runNumber)->second).find(kMSLTriggerClassName)->second;
    }
    return "";

}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetMSHTriggerClassName(Int_t runNumber) const
{
    if ( fPerRunInfo.count(runNumber) )
    {
    return (fPerRunInfo.find(runNumber)->second).find(kMSHTriggerClassName)->second;
    }
    return "";

}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetCentralityName(Int_t runNumber) const
{
    if ( fPerRunInfo.count(runNumber) )
    {
    return (fPerRunInfo.find(runNumber)->second).find(kCentralityName)->second;
    }
    return "";

}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMuConfig::GetListElements(const char* key, Bool_t simulation) const
{
    /// Get list as an array (to be deleted by the user)

    TObjArray* a = new TObjArray;
    a->SetOwner(kTRUE);

    THashList* list = static_cast<THashList*>(Map()->GetValue(key));

    TIter next(list);
    TObjString* s;
    TString rv;
    UInt_t test(kReal);

    if ( simulation )
    {
    test = kSim;
    }

    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
    if ( s->GetUniqueID() & test )
        {
        a->Add(new TObjString(*s));
        }
    }
    return a;
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMuConfig::GetTriggersList(Bool_t simulation) const
{
    /// Get list as an array (to be deleted by the user)

    TObjArray* a = new TObjArray;
    a->SetOwner(kTRUE);

    TObjArray* list = new TObjArray();
    TObjArray* list1 = static_cast<TObjArray*>(Map()->GetValue(RefMixTriggerKey()));
    TObjArray* list2 = static_cast<TObjArray*>(Map()->GetValue(DimuonTriggerKey()));
    TObjArray* list3 = static_cast<TObjArray*>(Map()->GetValue(MuonTriggerKey()));
    TObjArray* list4 = static_cast<TObjArray*>(Map()->GetValue(MinbiasTriggerKey()));
    TObjArray* list5 = static_cast<TObjArray*>(Map()->GetValue(MixTriggerKey()));

    if(list1){
        for (int i = 0; i < list1->GetEntries(); ++i){
            if ( !list->Contains(list1->At(i)->GetName()) )
                list->Add(list1->At(i)->Clone());
        }
    }

    if(list2){
        for (int i = 0; i < list2->GetEntries(); ++i){
            if ( !list->Contains(list2->At(i)->GetName()) )
                list->Add(list2->At(i)->Clone());
        }
    }

    if(list3){
        for (int i = 0; i < list3->GetEntries(); ++i){
            if ( !list->Contains(list3->At(i)->GetName()) )
                list->Add(list3->At(i)->Clone());
        }
    }

    if(list4){
        for (int i = 0; i < list4->GetEntries(); ++i){
            if ( !list->Contains(list4->At(i)->GetName()) )
                list->Add(list4->At(i)->Clone());
        }
    }

    if(list5){
        for (int i = 0; i < list5->GetEntries(); ++i){
            if ( !list->Contains(list5->At(i)->GetName()) )
                list->Add(list5->At(i)->Clone());
        }
    }

    TIter next(list);
    TObjString* s;
    TString rv;
    UInt_t test(kReal);

    if ( simulation ) test = kSim;

    while ( ( s = static_cast<TObjString*>(next()) ) )
        if ( s->GetUniqueID() & test )
            a->Add(new TObjString(*s));

    return a;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetTriggerClassName(ETriggerType tt, Int_t runNumber) const
{
    // get the triggerclass to for a given trigger type and run number

    if ( tt == kMB )
    {
    return GetMBTriggerClassName(runNumber);
    }
    if ( tt == kMIX )
    {
    return GetMIXTriggerClassName(runNumber);
    }
    else if ( tt == kMUL )
    {
    return GetMULTriggerClassName(runNumber);
    }
    else if ( tt == kMSL)
    {
    return GetMSLTriggerClassName(runNumber);
    }
    else if ( tt == kMSH)
    {
    return GetMSHTriggerClassName(runNumber);
    }
    else
    {
    AliError(Form("Unknown trigger type %d ???",tt));
    }
    return "";
}

//_____________________________________________________________________________
TString AliAnalysisMuMuConfig::GetTriggerTypeName(ETriggerType tt)
{
    // get the name of the trigger type
    if ( tt == kMB )
    {
    return "MB";
    }
    if ( tt == kMIX )
    {
    return "MIX";
    }
    else if ( tt == kMUL )
    {
    return "MUL";
    }
    else if ( tt == kMSL)
    {
    return "MSL";
    }
    else if ( tt == kMSH)
    {
    return "MSH";
    }
    return "";
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuConfig::Has(const char* key, const char* value, Bool_t simulation) const
{
    THashList* list = static_cast<THashList*>(Map()->GetValue(key));
    TIter next(list);
    TObjString* s;

    UInt_t test(kReal);

    if ( simulation )
    {
    test = kSim;
    }

    while ( ( s = static_cast<TObjString*>(next())) )
    {
    if ( ( s->GetUniqueID() & test ) && ( s->String() == value ) )
        {
        return kTRUE;
        }
    }
    return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuConfig::HasRunInformation(std::set<int>& runs, Bool_t show) const
{
    /// Return true if we have run information for all the runs in the set.
    /// if show=true, show that information for all those runs

    std::set<int>::size_type n(0);

    std::set<int>::const_iterator it;
    std::map<int,std::map<AliAnalysisMuMuConfig::EPerRunInfo,std::string> >::const_iterator mit;

    for ( it = runs.begin(); it != runs.end(); ++it )
    {
    int runNumber = *it;

    mit = fPerRunInfo.find(runNumber);

    if ( mit != fPerRunInfo.end() )
        {
        ++n;
        if ( show )
            {
            std::cout << Form("RUN %09d ",runNumber);
            PrintRunInfo(mit->second);
            }
        }
    else
        {
        if ( show )
            {
            std::cout << "Missing information for run " << runNumber << std::endl;
            }
        }
    }

    return n == runs.size();
}

//_____________________________________________________________________________
TMap* AliAnalysisMuMuConfig::Map() const
{
    /// Return (and initialize, if needed) the internal map
    if (!fMap)
    {
    fMap = new TMap;
    fMap->SetOwnerKeyValue(kTRUE,kTRUE);
    }
    return fMap;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::ShowList(const char* key, Bool_t simulation) const
{
    /// Show the list for a given key

    TObjArray* a = GetListElements(key,simulation);
    TObjString* s;

    TIter next(a);

    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
    std::cout << "    " << s->String() << std::endl;
    }

    delete a;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::Print(Option_t* opt) const
{
    /// printout
    /// Use opt = "REAL" to show only things relevant to real data
    /// Use opt = "SIM" to show only things relevant to simulation
    /// Use opt = "REAL SIM" or "" to show everything

    std::cout << "OCDBPath : " << fOCDBPath << std::endl;
    std::cout << "Use compact graphs : " << fIsCompactGraphs << std::endl;
    TString sopt(opt);
    sopt.ToUpper();

    if ( sopt == "RAW")
    {
    TIter nextKey(Map());
    TObjString* key;

    while ( ( key = static_cast<TObjString*>(nextKey()) ) )
        {
        THashList* list = static_cast<THashList*>(Map()->GetValue(key->String()));
        TIter next(list);
        TObjString* s;
        while ( ( s = static_cast<TObjString*>(next())))
            {
            std::cout << Form("%s [%u]",s->String().Data(),s->GetUniqueID()) << std::endl;
            }
        }
    return;
    }

    if (sopt.Length()==0)
    {
    sopt = "REAL SIM";
    }

    TIter next(Map());
    TObjString* key;

    while ( ( key = static_cast<TObjString*>(next()) ) )
    {
    std::cout << key->String() << ":" << std::endl;

    if ( sopt.Contains("REAL",TString::kIgnoreCase) )
        {
        std::cout << "  Real:" << std::endl;
        ShowList(key->String(),kFALSE);
        }
    if ( sopt.Contains("SIM",TString::kIgnoreCase) )
        {
        std::cout << "  Sim:" << std::endl;
        ShowList(key->String(),kTRUE);
        }
    }
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::SetRunInfo(const TString& runranges, const TString& runinfo)
{
    /// Associate the run info (specified in string runinfo) to each of the runs
    /// defined by the string runranges

    std::set<int> runs;

    DecodeRunRanges(runranges,runs);

    std::set<int>::const_iterator it;

    std::map<AliAnalysisMuMuConfig::EPerRunInfo,std::string> mapForOneRun;

    TObjArray* info = runinfo.Tokenize(" ");
    TIter next(info);
    TObjString*s ;
    while ( ( s = static_cast<TObjString*>(next())))
    {
    TString key,value;
    GetKeyValue(s->String(),'=',key,value);
    if (!key.CompareTo("MB",TString::kIgnoreCase))
        {
        mapForOneRun[kMBTriggerClassName]=value;
        }
    if (!key.CompareTo("MIX",TString::kIgnoreCase))
        {
        mapForOneRun[kMIXTriggerClassName]=value;
        }
    if (!key.CompareTo("MUL",TString::kIgnoreCase))
        {
        mapForOneRun[kMULTriggerClassName]=value;
        }
    if (!key.CompareTo("MSL",TString::kIgnoreCase))
        {
        mapForOneRun[kMSLTriggerClassName]=value;
        }
    if (!key.CompareTo("MSH",TString::kIgnoreCase))
        {
        mapForOneRun[kMSHTriggerClassName]=value;
        }
    if (!key.CompareTo("CENT",TString::kIgnoreCase))
        {
        mapForOneRun[kCentralityName]=value;
        }
    }

    for ( it = runs.begin(); it != runs.end(); ++it )
    {
    int runNumber = *it;
    fPerRunInfo[runNumber] = mapForOneRun;
    }

    delete info;
}

//_____________________________________________________________________________
void AliAnalysisMuMuConfig::ReadFromFile(const char* inputfile)
{
    /// Read configuration from external file.

    TString filename = gSystem->ExpandPathName(inputfile);
    std::ifstream in(filename.Data());
    if (in.bad())
    {
    AliError(Form("Cannot read input file %s",filename.Data()));
    return;
    }

    std::string line;

    std::vector<TString> selectionKeys;

    selectionKeys.push_back(RefMixTriggerKey());
    selectionKeys.push_back(RefMixEventSelectionKey());
    selectionKeys.push_back(DimuonTriggerKey());
    selectionKeys.push_back(MuonTriggerKey());
    selectionKeys.push_back(MinbiasTriggerKey());
    selectionKeys.push_back(MixTriggerKey());
    selectionKeys.push_back(EventSelectionKey());
    selectionKeys.push_back(EventSelectionMixKey());
    selectionKeys.push_back(PairSelectionKey());
    selectionKeys.push_back(CentralitySelectionKey());
    selectionKeys.push_back(FitTypeKey());
    selectionKeys.push_back(FitSingleKey());

    // read per run info
    while (std::getline(in,line))
    {
    TString sline(line.c_str());
    if ( sline.BeginsWith("#") ) continue;
    if ( sline.BeginsWith("//") ) continue;
    if (!sline.Contains(":")) continue;

    TString left, right; // relative to ":"

    if (!GetKeyValue(sline, ':',left,right)) continue;

    if ( left.Atoi() )
        {
        // we assume left is a run range, so we define runinfo
        SetRunInfo(left,right);
        }
    else if ( !left.CompareTo(OCDBPathKey(),TString::kIgnoreCase) )
        {
        fOCDBPath = right;
        fOCDBPath.Remove(TString::kBoth,' ');
        }
    else if ( !left.CompareTo(CompactGraphKey(),TString::kIgnoreCase) )
        {
        right.Remove(TString::kBoth,' ');
        if ( !right.CompareTo("yes",TString::kIgnoreCase) ||
            !right.CompareTo("on",TString::kIgnoreCase) ||
            !right.CompareTo("1",TString::kIgnoreCase) )
            {
            fIsCompactGraphs = kTRUE;
            }
        else
            {
            fIsCompactGraphs = kFALSE;
            }
        }
    else {

        Bool_t found(kFALSE);

        for ( std::vector<TString>::size_type i = 0; i < selectionKeys.size(); ++i )
        {
        if ( ! left.CompareTo(selectionKeys[i],TString::kIgnoreCase) )
            {
            Add(selectionKeys[i],right);
            found = kTRUE;
            }
        }

        if (!found)
        {
        std::cerr << "Unable to decode line : " << std::endl;
        std::cerr << line << std::endl;
        }
    }
    }
}


//_____________________________________________________________________________
void AliAnalysisMuMuConfig::SetColorScheme()
{
    /// Set a few custom colors

    new TColor(AliAnalysisMuMuConfig::kBlue,4/255.0,44/255.0,87/255.0,"my blue");
    new TColor(AliAnalysisMuMuConfig::kOrange,255/255.0,83/255.0,8/255.0,"my orange");
    new TColor(AliAnalysisMuMuConfig::kGreen,152/255.0,202/255.0,52/255.0,"my green");

    gStyle->SetGridColor(AliAnalysisMuMuConfig::kBlue);

    gStyle->SetFrameLineColor(AliAnalysisMuMuConfig::kBlue);
    gStyle->SetAxisColor(AliAnalysisMuMuConfig::kBlue,"xyz");
    gStyle->SetLabelColor(AliAnalysisMuMuConfig::kBlue,"xyz");

    gStyle->SetTitleColor(AliAnalysisMuMuConfig::kBlue);
    gStyle->SetTitleTextColor(AliAnalysisMuMuConfig::kBlue);
    gStyle->SetLabelColor(AliAnalysisMuMuConfig::kBlue);
    gStyle->SetStatTextColor(AliAnalysisMuMuConfig::kBlue);

    gStyle->SetOptStat(0);
}


// _____________________________________________________________________________
void  AliAnalysisMuMuConfig::LoadAliceStyles(){
  int font = 42;
  gROOT->SetStyle("Plain");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(1.1,"xy");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetMarkerSize(1.3);
  gStyle->SetPalette(1,0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(10);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetEndErrorSize(0);

}


