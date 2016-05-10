/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliMTRChEffAnalysis.h"

// ROOT includes
#include <Riostream.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "TList.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPRegexp.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TKey.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFileMerger.h"
#include "TFitResultPtr.h"

#include "AliLog.h"
#include "AliMergeableCollection.h"
#include "AliCounterCollection.h"

#include "AliAnalysisMuonUtility.h"
#include "AliTrigChEffOutput.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliMUONCDB.h"
#include "AliMUONTriggerEfficiencyCells.h"

/// \cond CLASSIMP
ClassImp(AliMTRChEffAnalysis) // Class implementation in ROOT context
/// \endcond

using std::cout;
using std::endl;
using std::cin;
using std::ofstream;
using std::ifstream;

//________________________________________________________________________
AliMTRChEffAnalysis::AliMTRChEffAnalysis() :
  TObject(),
  fOutputs(0x0),
  fConditions(0x0),
  fMergedOutputs(0x0)
{
  /// Default Ctor.

}

//________________________________________________________________________
AliMTRChEffAnalysis::AliMTRChEffAnalysis ( const char *localFileList, const char* outputName ) :
  TObject(),
  fOutputs(0x0),
  fConditions(0x0),
  fMergedOutputs(0x0)
{
  /// Ctor.
  SetOutList(localFileList,outputName);
}

//________________________________________________________________________
AliMTRChEffAnalysis::~AliMTRChEffAnalysis()
{
  //
  /// Destructor
  //
  delete fOutputs;
  delete fConditions;
  delete fMergedOutputs;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::AddSystematicCondition ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod )
{
  /// Add conditions for the systematic uncertainty estimation
  return SetCondition(physSel, trigClassName, centrality, itrackSel, imatch, imethod, kFALSE);
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::AddToList ( const char *filename, const char *outputName )
{
  /// Add to output list
  AliTrigChEffOutput* trigOut = new AliTrigChEffOutput(filename,outputName);
  if ( ! trigOut->GetMergeableCollection() ) {
    AliError(Form("Cannot find %s in %s\n",outputName,filename));
    return kFALSE;
  }

  TObjArray* condition = static_cast<TObjArray*>(fConditions->At(0));

//  // Delete counter collection to save memory
//  trigOut->RemoveFromList(trigOut->GetCounterCollection());
  // Just keep counter and mergeable collections
  TList* effHistoList = GetEffHistoList(trigOut,condition);
  effHistoList->SetName(filename);
  TString currRun = AliAnalysisMuonUtility::GetRunNumberAsString(filename);
  effHistoList->SetUniqueID(currRun.Atoi());

  if ( ! fOutputs ) {
    fOutputs = new TObjArray();
    fOutputs->SetOwner();
    fOutputs->SetName(outputName);
  }
  fOutputs->Add(effHistoList);
  delete trigOut;
  return kTRUE;
}

//________________________________________________________________________
TArrayI AliMTRChEffAnalysis::BoardsInRPC ( Int_t irpc ) const
{
  /// Return boards contained in the specified RPC

  // FIXME: ugly and hardcoded, but avoid loading the mapping

  TArrayI boards;
  if ( irpc == 0 || irpc == 9 ) {
    Int_t arr[] = {26,27,28,29,48,49,50,51,68,69,84,85,100,101,113};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 9 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 1 || irpc == 8 ) {
    Int_t arr[] = {9,10,11,30,31,32,33,52,53,54,55,70,71,86,87,102,103,114};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 8 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 2 || irpc == 7 ) {
    Int_t arr[] = {12,13,34,35,56,57,72,73,88,89,104,105,115};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 7 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 3 || irpc == 6 ) {
    Int_t arr[] = {14,15,36,37,58,59,74,75,90,91,106,107,116};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 6 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 4 || irpc == 5 ) {
    Int_t arr[] = {16,38,60,76,92,108,117};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 5 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 17 || irpc == 10 ) {
    Int_t arr[] = {6,7,8,22,23,24,25,44,45,46,47,66,67,82,83,98,99,112};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 10 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 16 || irpc == 11 ) {
    Int_t arr[] = {4,5,20,21,42,43,64,65,80,81,96,97,111};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 11 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 15 || irpc == 12 ) {
    Int_t arr[] = {2,3,18,19,40,41,62,63,78,79,94,95,110};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 12 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }
  else if ( irpc == 14 || irpc == 13 ) {
    Int_t arr[] = {1,17,39,61,77,93,109};
    Int_t nBoards = sizeof(arr)/sizeof(arr[0]);
    boards.Set(nBoards,arr);
    if ( irpc == 13 )
      for ( Int_t iboard=0; iboard<nBoards; iboard++ ) boards[iboard] += 117;
  }

  return boards;
}

//________________________________________________________________________
Int_t AliMTRChEffAnalysis::CompareEfficiencies ( const char* sources, const char* titles, const char* opt, const char* canvasNameSuffix ) const
{
  /// Compare efficiency objects
  TString srcs(sources);
  if ( srcs.Contains("raw://") ) {
    AliError("The method assumes that the specified storage is a SpecificStorage. Hence, please replace raw:// with the actual path in alien, e.g.: alien://folder=/alice/data/<year>/OCDB");
    return -2;
  }
  TObjArray* sourceList = srcs.Tokenize(",");
  TObjArray effMapList;
  effMapList.SetOwner();

  TString trigEffCDBdir = "MUON/Calib/TriggerEfficiency";

  TIter next(sourceList);
  TObjString* src = 0x0;
  while ( (src = static_cast<TObjString*>(next())) ) {
    TString currSrc = src->String();
    if ( currSrc.BeginsWith("alien") && ! gGrid ) TGrid::Connect("alien://");

    AliMUONTriggerEfficiencyCells* effMap = 0x0;
    if ( currSrc.EndsWith(".root") ) effMap = new AliMUONTriggerEfficiencyCells(currSrc.Data());
    else {
      TObjArray* dirRun = currSrc.Tokenize("?");
      TString cdbPath = dirRun->UncheckedAt(0)->GetName();
      TString runNum = ( dirRun->GetEntriesFast() > 1 ) ? dirRun->UncheckedAt(1)->GetName() : "";
      AliCDBManager* mgr = AliCDBManager::Instance();
      if ( ! mgr->GetDefaultStorage() ) mgr->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      if ( mgr->GetSpecificStorage(trigEffCDBdir.Data()) ) mgr->UnloadFromCache(trigEffCDBdir.Data());
      mgr->SetSpecificStorage(trigEffCDBdir.Data(),cdbPath.Data());
      Int_t runNumber = 0;
      if ( runNum.IsNull() ) {
        cout << "Please enter run number: " << endl;
        cin >> runNumber;
      }
      else runNumber = runNum.Atoi();
      mgr->SetRun(runNumber);
      AliCDBEntry* cdbEntry = mgr->Get(trigEffCDBdir.Data());

      // The CDB manager replace the current effMap with another in case two specific storage are provided
      // To avoid double deletion, the object must be cloned
      effMap = static_cast<AliMUONTriggerEfficiencyCells*>(cdbEntry->GetObject()->Clone());
    }
    if ( ! effMap ) continue;
    effMapList.Add(effMap);
  }

  return CompareEfficiencies(&effMapList, titles, opt, canvasNameSuffix);
}

//________________________________________________________________________
Int_t AliMTRChEffAnalysis::CompareEfficiencies ( TObjArray* effMapList, const char* titles, const char* opt, const char* canvasNameSuffix ) const
{
  /// Compare efficiency objects

  TString sTitles(titles);
  TObjArray* titleList = sTitles.Tokenize(",");

  Int_t nDiffs = 0;

  Int_t nLists = effMapList->GetEntriesFast();

  TString sCanvasNameSuffix(canvasNameSuffix);
  if ( ! sCanvasNameSuffix.IsNull() && ! sCanvasNameSuffix.BeginsWith("_") ) sCanvasNameSuffix.Prepend("_");

  Double_t xpt, ypt, xref, yref;
  enum {kEff, kDiff, kPull};
  TString sOpt(opt);
  sOpt.ToLower();
  Int_t iopt = kEff;
  TString yTitle = "Eff.";
  if ( sOpt.Contains("diff") ) {
    iopt = kDiff;
    yTitle = "Eff. - (ref.Eff)";
  }
  else if ( sOpt.Contains("pull") ) {
    iopt = kPull;
    yTitle = "(Eff - (ref.Eff)) / err";
  }
  else nDiffs = -1;

  Bool_t needsLegend = ( nLists > 1 );
//  if ( iopt != kEff ) needsLegend = nLists > 2;

  Int_t colors[] = {kBlack, kRed, kSpring, kTeal, kBlue, kViolet, kMagenta, kOrange, kGray};
  Int_t ncolors = sizeof(colors)/sizeof(colors[0]);

  Int_t hTypes[2] = {AliMUONTriggerEfficiencyCells::kHboardCount,AliMUONTriggerEfficiencyCells::kHslatCount};
  Int_t countTypes[3] = {AliMUONTriggerEfficiencyCells::kBendingEff,AliMUONTriggerEfficiencyCells::kNonBendingEff,AliMUONTriggerEfficiencyCells::kBothPlanesEff};

  TString currName = "";
  Int_t ican = 0;
  for ( Int_t itype=0; itype<2; itype++ ) {
    TString xTitle = ( hTypes[itype] == AliMUONTriggerEfficiencyCells::kHslatCount ) ? "RPC" : "Board";
    for ( Int_t icount=0; icount<3; icount++ ) {
      TCanvas* can = 0x0;
      for ( Int_t ich=0; ich<4; ich++ ) {
        TLegend* leg = 0x0;
//        Int_t hrefIdx = AliTrigChEffOutput::kNcounts + AliTrigChEffOutput::kNcounts*4*(itype-1) + 4*AliTrigChEffOutput::kAllTracks + ich;
//        Int_t hIdx = AliTrigChEffOutput::kNcounts + AliTrigChEffOutput::kNcounts*4*(itype-1) + 4*icount + ich;
        TGraphAsymmErrors* refGraph = 0x0;
        for ( Int_t ilist=0; ilist<nLists; ilist++ ) {
          TString currTitle = titleList->At(ilist)->GetName();
          AliMUONTriggerEfficiencyCells* effMap = static_cast<AliMUONTriggerEfficiencyCells*>(effMapList->UncheckedAt(ilist));
          currName = effMap->GetHistoName(hTypes[itype],AliMUONTriggerEfficiencyCells::kAllTracks,ich);
          TH1* histoDen = static_cast<TH1*>(effMap->GetHistoList()->FindObject(currName.Data()));
          currName = effMap->GetHistoName(hTypes[itype],countTypes[icount],ich);
          TH1* histoNum = static_cast<TH1*>(effMap->GetHistoList()->FindObject(currName.Data()));
          TGraphAsymmErrors* graph = new TGraphAsymmErrors(histoNum,histoDen,"e0");
          currName = histoNum->GetName();
          currName.ReplaceAll("Count","Eff");
          currName.Append(Form("_%s",currTitle.Data()));
          graph->SetName(currName.Data());

          if ( iopt != kEff ) {
            if ( refGraph ) {
              for ( Int_t ipt=0; ipt<graph->GetN(); ipt++ ) {
                refGraph->GetPoint(ipt,xref,yref);
                graph->GetPoint(ipt,xpt,ypt);
                Double_t diff = ypt - yref;
                if ( TMath::Abs(diff) > 1.e-4 ) nDiffs++;
                if ( iopt == kDiff ) graph->SetPoint(ipt,xpt,diff);
                else if ( iopt == kPull ) {
                  Double_t err = GetError(graph->GetErrorYlow(ipt),graph->GetErrorYhigh(ipt));
                  Double_t pull = ( err > 0. ) ? diff/err : 0.;
                  graph->SetPoint(ipt,xpt,pull);
                }
              } // loop on points
            }
            else {
              refGraph = graph;
              continue;
            }
          }

          if ( ! can ) {
            currName = graph->GetName();
            currName.Remove(currName.Length()-currTitle.Length()-5);
            currName += sCanvasNameSuffix;
            can = new TCanvas(currName.Data(),currName.Data(),20*ican,20*ican,600,600);
            can->Divide(2,2);
            ican++;
          }
          if ( needsLegend && ! leg ) {
            leg = new TLegend(0.35, 0.15, 0.65, 0.45);
            if ( refGraph ) leg->SetHeader(Form("Ref.: %s",titleList->At(0)->GetName()));
          }
          can->cd(ich+1);
          gPad->SetGridy();
          Int_t icolor = ( ilist < ncolors ) ? colors[ilist] : 20+ilist;
          graph->SetLineColor(icolor);
          graph->SetMarkerColor(icolor);
          graph->SetMarkerStyle(20+ilist);
          graph->SetMarkerSize(0.5);
          graph->GetXaxis()->SetTitle(xTitle.Data());
          graph->GetYaxis()->SetTitle(yTitle.Data());
          graph->Draw((gPad->GetListOfPrimitives()->GetEntries()==0)?"ap":"p");
          if ( leg ) leg->AddEntry(graph,titleList->At(ilist)->GetName(),"lp");
        } // loop on lists
        if (leg ) leg->Draw();
        delete refGraph;
      } // loop on chambers
    } // loop on counts
  } // loop on types

  delete titleList;

  return nDiffs;
}

//________________________________________________________________________
Int_t AliMTRChEffAnalysis::CompareEfficiencyMethods ( const char* source, const char* opt, const char* canvasNameSuffix ) const
{
  /// Compare efficiency methods
  if ( ! fConditions ) {
    AliWarning("No condition found! Please specify the default efficiency condition with SetEffConditions and then add additional tests with AddSystematicCondition");
    return -1;
  }

  AliTrigChEffOutput trigOut(source);

  Int_t nConditions = fConditions->GetEntriesFast();

  TObjArray condTitle(nConditions);
  condTitle.SetOwner();
  GetShortConditionTitles(&trigOut,condTitle);

  TString titles = "";
  TObjArray effMapList;
  effMapList.SetOwner();

  for ( Int_t icond=0; icond<nConditions; icond++ ) {
    TObjArray* condition = static_cast<TObjArray*>(fConditions->UncheckedAt(icond));

    TList* effList = GetEffHistoList(&trigOut, condition);
    AliMUONTriggerEfficiencyCells* effMap = new AliMUONTriggerEfficiencyCells(effList);
    effMapList.Add(effMap);
    titles += Form("%s,",condTitle.At(icond)->GetName());
  }

  titles.Remove(TString::kTrailing,',');
  return CompareEfficiencies(&effMapList, titles, opt, canvasNameSuffix);
}


//________________________________________________________________________
void AliMTRChEffAnalysis::CompareMergedEfficiencies ( const char* opt ) const
{
  if ( ! HasMergedResults() ) return;

  TObjArray* condition = static_cast<TObjArray*>(fConditions->At(0));
  TString titles = "";

  TObjArray effMapList;
  effMapList.SetOwner();

  AliTrigChEffOutput* trigOut = 0x0;
  TIter next(fMergedOutputs);
  while ( (trigOut = static_cast<AliTrigChEffOutput*>(next())) ) {
    TList* effList = GetEffHistoList(trigOut, condition);
    AliMUONTriggerEfficiencyCells* effMap = new AliMUONTriggerEfficiencyCells(effList);
    effMapList.Add(effMap);
    titles += Form("%s,",trigOut->GetName());
  }
  titles.Remove(TString::kTrailing,',');

  CompareEfficiencies(&effMapList, titles.Data(), opt, "MergedComp");
}


//________________________________________________________________________
//TArrayI AliMTRChEffAnalysis::CheckHomogeneity ( const TGraphAsymmErrors* trendGraph, Double_t nSigmas, Bool_t revert ) const
//{
//  /// Check efficiency homogeneity
//  TArrayI ranges(trendGraph->GetN());
//  Double_t xpt=0., ypt=0.;
//
//  Double_t sumEffOverSigma2 = -10.;
//  Double_t sumInvSimga2 = 0.001;
//  Int_t irange = 0;
//  Int_t npt = trendGraph->GetN();
//  for ( Int_t ipt=0; ipt<npt; ipt++ ) {
//    Int_t currPt = ( revert ) ? npt-1-ipt : ipt;
//    trendGraph->GetPoint(currPt,xpt,ypt);
//    Double_t err = GetError(trendGraph->GetErrorYlow(currPt),trendGraph->GetErrorYhigh(currPt));
//    Double_t invSigma2 = 1./(err*err);
//    Double_t effOverSigma2 = ypt*invSigma2;
//
//    Double_t meanEff = sumEffOverSigma2 / sumInvSimga2;
//    Double_t diff = ypt - meanEff;
//    Double_t diffErr = TMath::Sqrt(err*err+1./sumInvSimga2);
////    printf("Eff %g +- %g Mean eff %g +- %g  => %g/%g => %g\n",ypt,err,meanEff,1./sumInvSimga2,diff,diffErr,TMath::Abs(diff)/diffErr);
//    if ( TMath::Abs(diff)/diffErr > nSigmas ) {
//      TString run = trendGraph->GetXaxis()->GetBinLabel(currPt+1);
//      ranges[irange++] = run.Atoi();
//      sumEffOverSigma2 = 0.;
//      sumInvSimga2 = 0.;
//    }
//    sumEffOverSigma2 += effOverSigma2;
//    sumInvSimga2 += invSigma2;
//  }
//  TString run = trendGraph->GetXaxis()->GetBinLabel(trendGraph->GetN());
//  ranges[irange++] = run.Atoi();
//  ranges.Set(irange);
//  return ranges;
//}

//________________________________________________________________________
void AliMTRChEffAnalysis::CopyDir(TDirectory *source) const
{
  ///copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(source->GetName());
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write(obj->GetName(),TObject::kSingleKey);
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}


//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::CopyLocally ( const char* runList, const char* path, const char* pattern, const char* localFileList, const char* outDir, const char* directory ) const
{
  /// Copy files from QA repository or from grid
  TString sPattern(pattern);
  TString sOutDir(outDir);
  TString sPath(path);
  Bool_t isGrid = (! sPattern.IsNull());

  if ( sOutDir.IsNull() ) {
    if ( isGrid ) {
      TString data = sPath(TRegexp("/data/"));
      TString year = sPath(TRegexp("/20[0-9][0-9]/"));
      TString period = sPath(TRegexp("/LHC[0-9][0-9][a-z]"));
      period.Append("/");
      //    TPRegexp re("/(muon|pass).*?/");
      TString pass = AliAnalysisMuonUtility::GetPassName(path);
      if ( pass.IsNull() ) pass = AliAnalysisMuonUtility::GetPassName(pattern);
      sOutDir = data+year+period+pass;
      sOutDir.ReplaceAll("//","/");
      sOutDir.Remove(TString::kTrailing,'/');
      sOutDir.Remove(TString::kLeading,'/');
    }
    else sOutDir = sPath;
    sOutDir.Remove(TString::kTrailing,'/');
    sOutDir.Remove(TString::kLeading,'/');
  }

  TGridResult* res = 0x0;
  if ( isGrid ) {
    if ( ! gGrid ) TGrid::Connect("alien://");
    res = gGrid->Query(path,pattern);
  }

  TList* rl = GetRunList(runList);

  ofstream outFile(localFileList);
  Bool_t allOk = kTRUE;
  Bool_t prompt = kTRUE;
  Bool_t overwrite = kFALSE;
  for ( Int_t irun=0; irun<rl->GetEntries(); irun++ ) {
    TString run = static_cast<TObjString*>(rl->At(irun))->String();
    TString dest = Form("%s/%09i/QAresults.root",sOutDir.Data(),run.Atoi());
    TString destDir = gSystem->DirName(dest.Data());
    if ( gSystem->AccessPathName(destDir.Data()) ) ExecCommand(Form("mkdir -p %s",destDir.Data()),kFALSE);
    Bool_t copyFile = kTRUE;
    Bool_t isFileOk = kTRUE;
    if ( gSystem->AccessPathName(dest) == 0 ) {
      if ( prompt ) {
        TString decision = "";
        cout << "Local file " << dest.Data() << " already exist: overwrite? [y/n/ya/na (a=use decision for all)]" << endl;
        cin >> decision;
        if ( decision.EndsWith("a") ) prompt = kFALSE;
        if ( decision.BeginsWith("y") ) overwrite = kTRUE;
        else overwrite = kFALSE;
      }
      copyFile = overwrite;
    }
    if ( copyFile ) {
      TString src = "";
      if ( isGrid ) {
        for ( Int_t ifile=0; ifile<res->GetEntries(); ifile++ ) {
          TString inFilename = res->GetKey(ifile,"turl");
          if ( inFilename.Contains(run.Data()) ) {
            src = inFilename;
            break;
          }
        }
        if ( src.IsNull() ) {
          AliWarning(Form("Cannot find output for run %s",run.Data()));
          isFileOk = kFALSE;
        }
        else {
          TFile* outFile = TFile::Open(dest,"RECREATE");
          TFile* inFile = TFile::Open(src);
          inFile->cd();
          TDirectory* dir = static_cast<TDirectory*>(inFile->Get(directory));
          if ( dir ) {
            outFile->cd();
            CopyDir(dir);
          }
          else isFileOk = kFALSE;
          delete inFile;
          delete outFile;
        }
      }
      else {
        src = Form("http://aliqamu.web.cern.ch/aliqamu/%s",dest.Data());
        isFileOk = TFile::Cp(src.Data(),dest.Data());
      }
    }

    if ( isFileOk ) outFile << gSystem->pwd() << "/" << dest.Data() << endl;
    else allOk = kFALSE;
  }
  delete rl;
  outFile.close();
  return allOk;
}

//________________________________________________________________________
void AliMTRChEffAnalysis::DrawEffTrend ( Int_t itype, Int_t irpc, Double_t maxNsigmaOutliers, Double_t minEff, Double_t maxEff ) const
{
  /// Draw trenidng
  TString baseNames[3] = {"Chamber","RPC","Board"};
  TString base = baseNames[itype] + "Eff";
  if ( itype == AliTrigChEffOutput::kHboardEff ) {
    if ( irpc < 0 ) {
      AliWarning("Please specify RPC");
      return;
    }
    base += Form("InRPC%i",irpc);
  }

  Int_t nColumns = 2;
  Int_t nRows = 2;
  Int_t width = 600;
  Int_t height = 600;
  Int_t nDetEl = 4;
  Int_t nCh = 1;
  if ( itype != AliTrigChEffOutput::kHchamberEff ) {
    nColumns = 6;
    nRows = 3;
    width = 1200;
    height = 800;
    nDetEl = 18;
    nCh = 4;
  }
  TArrayI boards = BoardsInRPC(irpc);
  if ( itype == AliTrigChEffOutput::kHboardEff ) nDetEl = boards.GetSize();

  for ( Int_t ich=0; ich<nCh; ich++ ) {
    TString canName = base;
    if ( itype != AliTrigChEffOutput::kHchamberEff ) canName += Form("Ch%i",11+ich);
    TCanvas* can = new TCanvas(canName.Data(),canName.Data(),25*ich,25*ich,width,height);
    can->Divide(nColumns,nRows,0,0);
    //    can->SetTopMargin(0.03);
    //    can->SetBottomMargin(0.1);
    for ( Int_t idetelem=0; idetelem<nDetEl; idetelem++ ) {
      can->cd(idetelem+1);
      gPad->SetTicks(1,1);
      gPad->SetGridy();
      Int_t detElemId = idetelem;
      if ( itype == AliTrigChEffOutput::kHchamberEff ) detElemId = 11+idetelem;
      else if ( itype == AliTrigChEffOutput::kHboardEff ) detElemId = boards[idetelem];

      TString title = Form("%s %i",baseNames[itype].Data(),detElemId);
      TLegend* leg = new TLegend(0.2,0.15,0.8,0.4);
      leg->SetHeader(title.Data());
      for ( Int_t icount=0; icount<2; icount++ ) {
        TGraphAsymmErrors* gr = GetTrendEff(itype, icount, ich, detElemId);
        gr->SetLineColor(icount+1);
        gr->SetMarkerColor(icount+1);
        gr->SetMarkerStyle(24+2*icount);
        gr->GetYaxis()->SetRangeUser(minEff,maxEff);
        //        gr->GetYaxis()->SetLabelSize(0.025*nRows);
        //        gr->GetXaxis()->SetLabelSize(0.025*nColumns);
        gr->SetTitle("");
        gr->Draw(icount==0?"ap":"p");
        TString legTitle = ( icount==0 ) ? "bending plane" : "non-bending plane";
        leg->AddEntry(gr,legTitle.Data(),"lp");
        if ( maxNsigmaOutliers > 0. ) {
          TGraphAsymmErrors* outliers = GetOutliers(gr,maxNsigmaOutliers);
          outliers->SetLineColor(6+icount);
          outliers->SetMarkerColor(6+icount);
          outliers->SetMarkerStyle(20+2*icount);
          outliers->SetLineWidth(2);
          outliers->Draw("p");
          legTitle.ReplaceAll("plane","outliers");
          leg->AddEntry(outliers,legTitle.Data(),"lp");
        }
      }
      leg->Draw();
    }
    if ( itype == AliTrigChEffOutput::kHchamberEff ) break;
  }
}

//________________________________________________________________________
void AliMTRChEffAnalysis::DrawStatContribution ( Int_t itype, Int_t irpc, Double_t maxNsigmaOutliers, Double_t minY, Double_t maxY ) const
{
  /// Draw statistical contribution of the element efficiency
  if ( itype == AliTrigChEffOutput::kHchamberEff ) {
    AliWarning("This function is valid only for itype 1 and 2");
    return;
  }

  TString baseNames[3] = {"Chamber","RPC","Board"};
  TString base = baseNames[itype] + "Stat";
  if ( itype == AliTrigChEffOutput::kHboardEff ) {
    if ( irpc < 0 ) {
      AliWarning("Please specify RPC");
      return;
    }
    base += Form("InRPC%i",irpc);
  }

  Int_t nColumns = 6;
  Int_t nRows = 3;
  Int_t width = 1200;
  Int_t height = 800;
  Int_t nDetEl = 18;
  Int_t nCh = 1;

  TArrayI boards = BoardsInRPC(irpc);
  if ( itype == AliTrigChEffOutput::kHboardEff ) nDetEl = boards.GetSize();

  for ( Int_t ich=0; ich<nCh; ich++ ) {
    TH1* histos[nDetEl];
    TH1* sumHistos = 0x0;
    for ( Int_t idetelem=0; idetelem<nDetEl; idetelem++ ) {
      Int_t detElemId = idetelem;
      if ( itype == AliTrigChEffOutput::kHboardEff ) detElemId = boards[idetelem];
      histos[idetelem] = GetTrend(itype, AliTrigChEffOutput::kAllTracks, ich, detElemId);
      histos[idetelem]->SetName(Form("%s_stat",histos[idetelem]->GetName()));
      histos[idetelem]->SetStats(0);
      if ( sumHistos ) sumHistos->Add(histos[idetelem]);
      else sumHistos = static_cast<TH1*>(histos[idetelem]->Clone("sumHistos"));
    }

    TString canName = base; //Form("%sCh%i",base.Data(),11+ich);
    TCanvas* can = new TCanvas(canName.Data(),canName.Data(),25*ich,25*ich,width,height);
    can->Divide(nColumns,nRows,0,0);
    for ( Int_t idetelem=0; idetelem<nDetEl; idetelem++ ) {
      can->cd(idetelem+1);
      gPad->SetTicks(1,1);
      gPad->SetGridy();
      Int_t detElemId = idetelem;
      if ( itype == AliTrigChEffOutput::kHboardEff ) detElemId = boards[idetelem];
      TString title = Form("%s %i",baseNames[itype].Data(),detElemId);
      TLegend* leg = new TLegend(0.2,0.65,0.8,0.9);
      leg->SetHeader(title.Data());
      TGraphAsymmErrors* gr = new TGraphAsymmErrors(histos[idetelem],sumHistos,"e0");
      gr->SetHistogram(static_cast<TH1F*>(histos[idetelem]));

      gr->SetMarkerStyle(24);
      gr->SetMarkerSize(0.5);
      gr->GetYaxis()->SetRangeUser(minY,maxY);
      gr->GetYaxis()->SetTitle(Form("Tracks in %s / Sum of tracks of %ss in %s",baseNames[itype].Data(),baseNames[itype].Data(),baseNames[itype-1].Data()));
      gr->SetTitle("");
      gr->Draw("ap");
      TString legTitle = "stat";
//      leg->AddEntry(gr,legTitle.Data(),"lp");
      if ( maxNsigmaOutliers > 0. ) {
        TGraphAsymmErrors* outliers = GetOutliers(gr,maxNsigmaOutliers);
        outliers->SetLineColor(6);
        outliers->SetMarkerColor(6);
        outliers->SetMarkerStyle(20);
        outliers->SetLineWidth(2);
        outliers->Draw("p");
        legTitle = "outliers";
        leg->AddEntry(outliers,legTitle.Data(),"lp");
      }
//      }
      leg->Draw();
    }
    delete sumHistos;
  }
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::DrawSystematicEnvelope ( Bool_t perRPC ) const
{
  /// Get systematic envelop for merged efficiencies
  if ( ! HasMergedResults() ) return kFALSE;

  Int_t colors[] = {kBlack, kRed, kSpring, kTeal, kBlue, kViolet, kMagenta, kOrange, kGray};
  Int_t ncolors = sizeof(colors)/sizeof(colors[0]);

  Int_t itype = ( perRPC ) ? AliMUONTriggerEfficiencyCells::kHslatCount : AliMUONTriggerEfficiencyCells::kHboardCount;
  Int_t countTypes[2] = {AliMUONTriggerEfficiencyCells::kBendingEff,AliMUONTriggerEfficiencyCells::kNonBendingEff};

  Int_t nConditions = fConditions->GetEntriesFast();

  TObjArray condTitle(nConditions);
  condTitle.SetOwner();

  Double_t xpt, ypt, xref, yref;
  TArrayD eff(4), effErr(4);

  AliTrigChEffOutput* trigOut = 0x0;
  TIter next(fMergedOutputs);
  Int_t imerged = -1;
  while ( (trigOut = static_cast<AliTrigChEffOutput*>(next())) ) {
    imerged++;

    // Get meaningful short titles for each systematic
    if ( imerged == 0 ) GetShortConditionTitles(trigOut,condTitle);

    TArrayI isEmpty(nConditions);

    // First get the needed graphs
    Int_t nDE = 0;
    TObjArray effGraphs[8];
    for ( Int_t iplane=0; iplane<8; iplane++ ) {
      effGraphs[iplane].SetOwner();
    }
    TString titles = "";
    TObjArray effMapList;
    effMapList.SetOwner();
    for ( Int_t icond=0; icond<nConditions; icond++ ) {
      TObjArray* condition = static_cast<TObjArray*>(fConditions->UncheckedAt(icond));

      TList* effList = GetEffHistoList(trigOut, condition);
      AliMUONTriggerEfficiencyCells* effMap = new AliMUONTriggerEfficiencyCells(effList);
      effMapList.Add(effMap);
      titles += Form("%s,",condTitle.At(icond)->GetName());

      for ( Int_t ich=0; ich<4; ich++ ) {
        TString currName = effMap->GetHistoName(itype,AliMUONTriggerEfficiencyCells::kAllTracks,ich);
        TH1* histoDen = static_cast<TH1*>(effMap->GetHistoList()->FindObject(currName.Data()));
        for ( Int_t icount=0; icount<2; icount++ ) {
          Int_t iplane = 4*icount+ich;

          if ( histoDen->GetEntries() > 0 ) {
            currName = effMap->GetHistoName(itype,countTypes[icount],ich);
            TH1* histoNum = static_cast<TH1*>(effMap->GetHistoList()->FindObject(currName.Data()));
            TGraphAsymmErrors* gr = new TGraphAsymmErrors(histoNum,histoDen,"e0");
            nDE = gr->GetN();
            effGraphs[iplane].AddAtAndExpand(gr,icond);
          }
          else {
            isEmpty[icond] = 1;
            AliWarning(Form("No entries in count %i and ch %i for %s\n",icount,ich,condTitle.At(icond)->GetName()));
          }
        }
      }
    }

    titles.Remove(TString::kTrailing,',');
    CompareEfficiencies(&effMapList, titles, "diff", trigOut->GetName());

    // Draw average dispersion per plane
    TString canName = Form("EffSyst_%s",trigOut->GetName());
    Int_t pos = 25*(imerged+1);
    TCanvas* can = new TCanvas(canName.Data(),canName.Data(),pos,pos,1200,800);
    can->Divide(4,2,0,0);

    for ( Int_t icount=0; icount<2; icount++ ) {
      for ( Int_t ich=0; ich<4; ich++ ) {
        Int_t iplane = 4*icount+ich;
        can->cd(iplane+1);
        gPad->SetTicks(1,1);
        gPad->SetLogy();
        TLegend* leg = new TLegend(0.15,0.7,0.9,0.9);
        leg->SetHeader(trigOut->GetHistoName(-1,icount,ich,-1,-1,-1));
        TH1* sumHisto = 0x0;
        for ( Int_t icond=1; icond<nConditions; icond++ ) {
          if ( isEmpty[icond] == 1 ) continue;
          TH1* histo = new TH1D(Form("syst_%s_%s_plane%i_ch%i",trigOut->GetName(),condTitle[icond]->GetName(),icount,11+ich),"",200,-0.1,0.1);
          histo->GetXaxis()->SetTitle("Eff.-(ref.Eff.)");
          histo->GetYaxis()->SetTitle("1/#sigma^{2}");

          TGraphAsymmErrors* gr = static_cast<TGraphAsymmErrors*>(effGraphs[iplane].UncheckedAt(icond));
          TGraphAsymmErrors* grRef = static_cast<TGraphAsymmErrors*>(effGraphs[iplane].UncheckedAt(0));
          for ( Int_t ipt=0; ipt<gr->GetN(); ipt++ ) {
            gr->GetPoint(ipt,xpt,ypt);
            Double_t err = GetError(gr->GetErrorYlow(ipt),gr->GetErrorYhigh(ipt));
            Double_t invErr2 = ( err > 0. ) ? 1./(err*err) : 0.;
            grRef->GetPoint(ipt,xref,yref);

            Double_t diff = ypt-yref;
            histo->Fill(diff,invErr2);
          }

          if ( ! sumHisto ) {
            sumHisto = static_cast<TH1*>(histo->Clone(Form("syst_%s_plane%i_ch%i",trigOut->GetName(),icount,11+ich)));
            sumHisto->SetLineColor(1);
            sumHisto->Draw();
            leg->AddEntry(sumHisto,"All systematics","l");
          }
          else sumHisto->Add(histo);
          Int_t icolor = ( icond < ncolors ) ? colors[icond] : 20+icond;
          histo->SetLineColor(icolor);
          histo->Draw("same");
          leg->AddEntry(histo,condTitle[icond]->GetName(),"l");
        } // loop on conditions
        leg->Draw();
      } // loop on chambers
    } // loop on counts


    canName = Form("TriggerEff_3outOf4_syst_%s",trigOut->GetName());
    pos += 50;
    TCanvas* canSyst = new TCanvas(canName.Data(),canName.Data(),pos,pos,600,600);
    canSyst->SetLogy();
    TLegend* leg = new TLegend(0.15,0.7,0.9,0.4);
//    leg->SetHeader(trigOut->GetHistoName(-1,icount,ich,-1,-1,-1));
    TH1* histo[nConditions];
    for ( Int_t icond=0; icond<nConditions; icond++ ) {
      if ( isEmpty[icond] == 1 ) continue;
      histo[icond] = new TH1D(Form("TriggerEff_syst_%s_%s",condTitle[icond]->GetName(),trigOut->GetName()),"Dispersion of trigger probability (3/4)",200,-0.1,0.1);
      histo[icond]->GetXaxis()->SetTitle("Trig. prob. - (ref. trig. prob)");
      histo[icond]->GetYaxis()->SetTitle("1/#sigma^{2}");
    }

    for ( Int_t ipt=0; ipt<nDE; ipt++ ) {
      Double_t refTrigProb = 0.; // refTrigProbErr = 0.;
      for ( Int_t icond=0; icond<nConditions; icond++ ) {
        if ( isEmpty[icond] == 1 ) continue;
        Double_t trigProb = 1., trigProbErr2 = 0.;
        for ( Int_t icount=0; icount<2; icount++ ) {
          for ( Int_t ich=0; ich<4; ich++ ) {
            Int_t iplane = 4*icount+ich;
            TGraphAsymmErrors* gr = static_cast<TGraphAsymmErrors*>(effGraphs[iplane].UncheckedAt(icond));
            gr->GetPoint(ipt,xpt,ypt);
            eff[ich] = ypt;
            effErr[ich] = GetError(gr->GetErrorYlow(ipt),gr->GetErrorYhigh(ipt));
          } // loop on chambers
          Double_t effErr34 = 0.;
          Double_t eff34 = GetThreeOfFour(eff,effErr,effErr34);
          trigProb *= eff34;
          trigProbErr2 += effErr34*effErr34;
        } // loop on counts

        if ( icond == 0 ) {
          refTrigProb = trigProb;
//          refTrigProbErr = trigProbErr;
        }
        else {
          Double_t invErr2 = ( trigProbErr2>0. ) ? 1./trigProbErr2 : 0.;
          histo[icond]->Fill(trigProb-refTrigProb,invErr2);
        }
      } // loop on conditions
    } // loop on points

    for ( Int_t icond=0; icond<nConditions; icond++ ) {
      if ( isEmpty[icond] == 1 ) continue;
      TString title = ( icond == 0 ) ? "All systematics" : condTitle[icond]->GetName();
      Int_t icolor = ( icond < ncolors ) ? colors[icond] : 20+icond;
      histo[icond]->SetLineColor(icolor);
      histo[icond]->Draw((icond == 0)?"":"same");
      leg->AddEntry(histo[icond],title.Data(),"l");
      if ( icond>0 ) histo[0]->Add(histo[icond]);
    }
    leg->Draw();
  } // loop on merged output
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::ExecCommand ( TString command, Bool_t prompt ) const
{
  /// Execute command
  TString decision = "y";

  if ( gROOT->IsBatch() ) prompt = kFALSE; // To run with crontab

  if ( prompt ) {
    cout << command.Data() << " ? [y/n]" << endl;
    cin >> decision;
  }

  decision.ToLower();
  if ( decision == "y" ) {
    cout << "Executing: " << command.Data() << endl;
    gSystem->Exec(command.Data());
    return kTRUE;
  }

  return kFALSE;
}

//________________________________________________________________________
Double_t AliMTRChEffAnalysis::FitRangesFunc ( Double_t* x, Double_t* par )
{
  /// Function with multiple constant fixes
  /// Parameters are:
  /// [0] = number of break points dividing 2 sub ranges
  /// [1] = value of efficiency in the first (or only) subrange
  /// [2*ipar] (for ipar>=1) = position in x where efficiency value can change
  /// [2*ipar+1] = value of efficiency in the subrange above [2*ipar]

  Double_t xx = x[0];
  Double_t val = par[1];
  Int_t nChanges = par[0];
  Double_t matchDiff = 123456789.;
  Int_t matchChange = -1;
  for ( Int_t iknot=0; iknot<nChanges; iknot++ ) {
    Int_t iparChange = 2*(iknot+1);
    Double_t diff = xx - par[iparChange];
    if ( diff >= 0. && diff < matchDiff ) {
      matchDiff = diff;
      matchChange = iparChange;
    }
  }
  if ( matchChange >= 0 ) val = par[matchChange+1];

  return val;
}

//________________________________________________________________________
Double_t AliMTRChEffAnalysis::GetAverageStat ( Int_t firstRun, Int_t lastRun, Int_t itype,Bool_t excludePeriphericBoard ) const
{
  TH1* statHisto = 0x0;

  UInt_t uFirstRun = (UInt_t)firstRun;
  UInt_t uLastRun = (UInt_t)lastRun;

  TIter next(fOutputs);
  TList* effHistoList = 0x0;
  while ( (effHistoList = static_cast<TList*>(next()) ) ) {
    UInt_t run = effHistoList->GetUniqueID();
    if ( run < uFirstRun || run > uLastRun ) continue;
    TH1* histo = GetHisto(effHistoList,itype,AliTrigChEffOutput::kAllTracks,0);
//    if ( ! histo ) continue;
    if ( statHisto ) statHisto->Add(histo);
    else statHisto = static_cast<TH1*>(histo->Clone("tmpStatHisto"));
//    delete histo;
  }
  if ( ! statHisto ) return 0.;

  Double_t statPerDE = 0.;
  Double_t nDe = 0.;
  if ( itype == AliTrigChEffOutput::kHboardEff && excludePeriphericBoard ) {
    Int_t excludeBoardsHalf[] = {1, 17, 39, 61, 77, 93,109,
                                16, 38, 60, 76, 92, 108, 117,
                                110, 111, 112, 113, 114, 115, 116};
    Int_t nExcludedHalf = sizeof(excludeBoardsHalf)/sizeof(excludeBoardsHalf[0]);
    for ( Int_t ibin=1; ibin<=statHisto->GetNbinsX(); ibin++ ) {
      Bool_t skip = kFALSE;
      for ( Int_t iexcl=0; iexcl<nExcludedHalf; iexcl++ ) {
        if ( ibin == excludeBoardsHalf[iexcl] || ibin == excludeBoardsHalf[iexcl]+117 ) {
          skip = kTRUE;
          break;
        }
      }
      if ( skip ) continue;
      statPerDE += statHisto->GetBinContent(ibin);
      nDe += 1.;
    }
  }
  else {
    statPerDE = statHisto->Integral();
    nDe = (Double_t)statHisto->GetNbinsX();
  }
  statPerDE = nDe > 0. ? statPerDE/nDe : -1.;

  delete statHisto;
  return statPerDE;
}

//________________________________________________________________________
TList* AliMTRChEffAnalysis::GetEffHistoList ( AliTrigChEffOutput* trigOut, TObjArray* condition ) const
{
  /// Get the list of objects for the efficiency calculation
  Int_t itrackSel = static_cast<TObjString*>(condition->At(3))->String().Atoi();
  Int_t imatch = static_cast<TObjString*>(condition->At(4))->String().Atoi();
  Int_t imethod = static_cast<TObjString*>(condition->At(5))->String().Atoi();

  return trigOut->GetEffHistoList(condition->At(0)->GetName(),condition->At(1)->GetName(),condition->At(2)->GetName(),itrackSel,imatch,imethod);
}

//________________________________________________________________________
Double_t AliMTRChEffAnalysis::GetError ( Double_t errLow, Double_t errHigh ) const
{
  /// Get error from the asymmetric error on efficiency
  return TMath::Max(errLow,errHigh);
}

//________________________________________________________________________
TH1* AliMTRChEffAnalysis::GetHisto ( TList* effHistoList, Int_t itype, Int_t icount, Int_t ichamber ) const
{
  /// Get histogram
  Int_t ihisto = ( itype == AliTrigChEffOutput::kHchamberEff ) ? icount : AliTrigChEffOutput::kNcounts + 4*AliTrigChEffOutput::kNcounts*(itype-1) + 4*icount + ichamber;
  return static_cast<TH1*>(effHistoList->At(ihisto));
}

//________________________________________________________________________
TArrayI AliMTRChEffAnalysis::GetHomogeneousRanges ( Double_t chi2Cut, Int_t maxNRanges, Double_t minEffVariation, Bool_t perRPC, TArrayI* forcedChanges, Double_t minEff, Double_t maxEff )
{
  /// Get run ranges with an efficiency compatible with constant

  Int_t nRuns = fOutputs->GetEntriesFast();

  TH1F* hRunChangeCount = new TH1F("runChangeCount","Number of RPCs changing efficiency per run",nRuns,-0.5,-0.5+(Double_t)nRuns);
  hRunChangeCount->SetXTitle("Run num.");
  hRunChangeCount->SetYTitle("Num. of RPCs with change in eff.");
  for ( Int_t irun=0; irun<nRuns; irun++ ) {
    hRunChangeCount->GetXaxis()->SetBinLabel(irun+1,Form("%i",GetRunNumber(irun)));
  }

  Int_t itype = perRPC ? AliTrigChEffOutput::kHslatEff : AliTrigChEffOutput::kHboardEff;

  Int_t nCanvas = perRPC ? 4 : 18;
  TObjArray canList(nCanvas);

  for ( Int_t irpc=0; irpc<18; irpc++ ) {
    Int_t icount = AliTrigChEffOutput::kBothPlanesEff;
    TArrayI boards = BoardsInRPC(irpc);
    Int_t firstDetEl = perRPC ? irpc : 0;
    Int_t lastDetEl = perRPC ? irpc : boards.GetSize()-1;

    for ( Int_t ich=0; ich<4; ich++ ) {
      Int_t ican = perRPC ? ich : irpc;
      TCanvas* can = static_cast<TCanvas*>(canList.At(ican));
      if ( ! can ) {
        TString canName = perRPC ? Form("testRanges_ch%i",11+ich) : Form("testRanges_RPC%i",irpc);
        can = new TCanvas(canName.Data(),canName.Data(),10*ich,10*ich,1200,800);
        can->Divide(6,3,0,0);
        can->SetMargin(0.,0.,0.,0.);
        canList.AddAt(can,ican);
      }

      for ( Int_t idetel=firstDetEl; idetel<=lastDetEl; idetel++ ) {
        Int_t currDE = ( perRPC ) ? idetel : boards[idetel];
//      for ( Int_t icount=0; icount<2; icount++ ) {
        TGraphAsymmErrors* trendGraph = GetTrendEff(itype,icount,ich,currDE);
        TArrayI range = GetHomogeneousRanges(trendGraph,chi2Cut,maxNRanges,minEffVariation,forcedChanges, kTRUE);
        trendGraph->GetYaxis()->SetRangeUser(minEff,maxEff);

        can->cd(idetel+1);
        gPad->SetTicks(1,1);
        gPad->SetMargin(0.08,0.,0.08,0.);
        TString drawOpt = ( gPad->GetListOfPrimitives()->GetEntries() == 0 ) ? "ap" : "p";

        TString legendName = Form("%s_%i",can->GetName(),currDE);
        TLegend* leg = static_cast<TLegend*>(gPad->GetListOfPrimitives()->FindObject(legendName.Data()));
        if ( ! leg ) {
          leg = new TLegend(0.2,0.15,0.8,0.4);
          leg->SetHeader(Form("%s %i",perRPC?"RPC":"Board",currDE));
          leg->SetName(legendName.Data());
        }

        if ( ! perRPC ) {
          Int_t icolor = ich+1;
          trendGraph->SetLineColor(icolor);
          trendGraph->SetMarkerColor(icolor);
          trendGraph->SetMarkerStyle(24+ich);
          TF1* func = static_cast<TF1*>(trendGraph->GetListOfFunctions()->At(0));
          if ( func ) {
            func->SetLineWidth(2);
            func->SetLineColor(icolor);
          }
        }
        trendGraph->GetXaxis()->SetLabelSize(0.07);
        trendGraph->SetTitle("");

        trendGraph->Draw(drawOpt.Data());
        leg->AddEntry(trendGraph,Form("Chamber %i",11+ich),"lp");
        leg->Draw();
        for ( Int_t ichange=2; ichange<range.GetSize(); ichange++ ) {
          // Store only the run when the change applies
          if ( ichange%2 == 1 ) continue;
          Int_t runIdx = range[ichange];
//          if ( ichange != 0 ) {
            TLine* line = new TLine(runIdx,minEff,runIdx,maxEff);
            line->SetLineStyle(2);
            line->Draw("same");
            TLatex text;
            text.DrawLatex((Double_t)(runIdx-3),maxEff-0.1*(maxEff-minEff)*(Double_t)(ichange/2),Form("%i",GetRunNumber(runIdx)));
//          }
          hRunChangeCount->Fill(runIdx);
          if ( hRunChangeCount->GetBinContent(runIdx+1) == 1 ) {
            TString infoMsg = Form("Efficiency change in %i triggered by ch %i  RPC %i",GetRunNumber(runIdx),11+ich,irpc);
            if ( ! perRPC ) infoMsg += Form(" Board %i",currDE);
            AliInfo(infoMsg.Data());
          }
        } // loop on change
      } // loop on detection element
//      }
    } // loop on chambers
  } // loop on RPC

  // Clusterize contiguous runs
  TArrayI runChangeClust(nRuns);
  Double_t sumWgtRun = 0.;
  Double_t sumWgt = 0;
  for ( Int_t irun=0; irun<=nRuns; irun++ ) {
    if ( irun == nRuns || hRunChangeCount->GetBinContent(irun+1) == 0 ) {
      if ( sumWgt > 0. ) {
        Int_t averageRun = TMath::Nint(sumWgtRun / sumWgt);
        AliDebug(2,Form("Average run: %i => %i\n",averageRun,GetRunNumber(averageRun)));
        runChangeClust[averageRun]++;
        sumWgtRun = 0.;
        sumWgt = 0.;
      }
    }
    if ( irun == nRuns ) break;

    AliDebug(2,Form("irun %i => %i: wgt %g",irun,GetRunNumber(irun),hRunChangeCount->GetBinContent(irun+1)));

//    Double_t wgt = (Double_t)runChangeCount[irun];
    Double_t wgt = hRunChangeCount->GetBinContent(irun+1);
    if ( forcedChanges ) {
      for ( Int_t ichange=0; ichange<forcedChanges->GetSize(); ichange++ ) {
        if ( GetRunNumber(irun) == forcedChanges->At(ichange) ) wgt *= 10.;
      }
    }
    sumWgtRun += wgt*(Double_t)irun;
    sumWgt += wgt;
  }

  TCanvas* summaryCan = new TCanvas("effChangeSummary","effChangeSummary",50,50,600,600);
  summaryCan->SetLogy();
  hRunChangeCount->GetXaxis()->LabelsOption("v");
  hRunChangeCount->Draw();

  Int_t ientry = 0;
  TArrayI runRanges(nRuns);
  runRanges[ientry++] = GetRunNumber(0);
  for ( Int_t irun=1; irun<nRuns; irun++ ) {
    if ( runChangeClust[irun] > 0 ) {
      runRanges[ientry++] = GetRunNumber(irun-1);
      runRanges[ientry++] = GetRunNumber(irun);
    }
  }
  runRanges[ientry++] = GetRunNumber(nRuns-1);
  runRanges.Set(ientry);
  return runRanges;
}

//________________________________________________________________________
TArrayI AliMTRChEffAnalysis::GetHomogeneousRanges ( TGraphAsymmErrors* trendGraph, Double_t chi2Cut, Int_t maxNRanges, Double_t minEffVariation, TArrayI* forcedChanges, Bool_t returnIndex )
{
  /// Get run ranges with an efficiency compatible with constant

  TArrayI runRanges;
  TF1* func = 0x0;
  TArrayD forcedChangesBin;
  Int_t nForced = 0;
  if ( forcedChanges ) {
    forcedChangesBin.Set(forcedChanges->GetSize());
    for ( Int_t ichange=0; ichange<forcedChanges->GetSize(); ichange++ ) {
      Int_t idx = GetIndexFromRun(forcedChanges->At(ichange));
      if ( idx >= 0 ) forcedChangesBin[nForced++] = (Double_t)idx;
      else AliWarning(Form("Cannot find run %i\n",forcedChanges->At(ichange)));
    }
  }

  Double_t minNormChi2 = 123456789.;
  Int_t minNormChi2Step = -1;
  TString fitOpt = "NQ0";
//  Double_t fakeVal = -999.;
//  TArrayD params(2*maxNRanges);
//  params.Reset(fakeVal);

  for ( Int_t istep=0; istep<maxNRanges; istep++ ) {
    Int_t nPars = 2*(istep+1);
    Double_t xMin = trendGraph->GetXaxis()->GetXmin();
    Double_t xMax = trendGraph->GetXaxis()->GetXmax();
    func = new TF1("rangeFunc",this,&AliMTRChEffAnalysis::FitRangesFunc,xMin,xMax,nPars,"AliMTRChEffAnalysis","FitRanges");
    func->FixParameter(0,istep);
    for ( Int_t ipar=1; ipar<nPars; ipar++ ) {
//      if ( ipar >= 2*istep ) {
//      if ( TMath::Abs(params[ipar]-fakeVal) < 0.01 ) {
//        if ( ipar%2 == 1 ) params[ipar] = 0.95;
//        else params[ipar] = (xMax-xMin) * (Double_t)(ipar/2)/((Double_t)(istep+1));
//      }
//      func->SetParameter(ipar,params[ipar]);
      Double_t val = ( ipar%2 == 1 ) ? 0.95 : (xMax-xMin) * (Double_t)(ipar/2)/((Double_t)(istep+1));
      func->SetParameter(ipar,val);
      if ( ipar%2 == 0 ) func->SetParLimits(ipar,xMin,xMax);
    }

    TFitResultPtr fitResult = trendGraph->Fit(func,fitOpt.Data());
    // If fit converges close to a point where a break is forced
    // fix parameters and redo the fit
    if ( forcedChanges ) {
      Bool_t hasFixedPars = kFALSE;
      for ( Int_t iforced=0; iforced<nForced; iforced++ ) {
        for ( Int_t jstep=0; jstep<istep; jstep++ ) {
          Int_t ipar = 2*(jstep+1);
          if ( TMath::Abs(forcedChangesBin[iforced]-func->GetParameter(ipar)) > 2. ) continue;
          func->FixParameter(ipar,forcedChangesBin[iforced]);
          hasFixedPars = kTRUE;
        }
      }
      if ( hasFixedPars ) fitResult = trendGraph->Fit(func,fitOpt.Data());
    }

//    for ( Int_t ipar=1; ipar<nPars; ipar++ ) {
//      params[ipar] = func->GetParameter(ipar);
//    }

    Double_t normChi2 = func->GetChisquare() / ((Double_t)func->GetNDF());
    if ( normChi2 < minNormChi2 ) {
      minNormChi2 = normChi2;
      minNormChi2Step = istep;
    }
//    if ( static_cast<int>(fitResult) == 0 && normChi2 < chi2Cut ) break;
    if ( normChi2 < chi2Cut ) break; // REMEMBER TO CHECK
    delete func;
    func = 0x0;
  }

  if ( func ) {
    trendGraph->GetListOfFunctions()->Add(func->Clone());
    Int_t nSteps = (Int_t)func->GetParameter(0);

    // The runs for which the efficiency changes could be unsorted
    // (when forced values are requires)
    // Copy the parameters in arrays to sort them
    Int_t nPoints = nSteps+1;
    TArrayD parRunIdx(nPoints);
    TArrayD parEff(nPoints);
    for ( Int_t ipar=0; ipar<func->GetNpar(); ipar++ ) {
      Int_t istep = ipar/2;
      if ( ipar%2 == 0 ) parRunIdx[istep] = func->GetParameter(ipar);
      else parEff[istep] = func->GetParameter(ipar);
    }
    parRunIdx[0] = 0.;
    TArrayI sortIdx(nPoints);
    TMath::Sort(nPoints,parRunIdx.GetArray(),sortIdx.GetArray(),kFALSE);

    runRanges.Set(2*nPoints);
    Int_t irun = 0;
    runRanges[irun++] = returnIndex ? 0 : GetRunNumber(0);
    for ( Int_t ipoint=1; ipoint<nPoints; ipoint++ ) {
      Double_t deltaEff = TMath::Abs(parEff[sortIdx[ipoint]]-parEff[sortIdx[ipoint-1]]);
//      if ( ipoint>=2 ) deltaEff = TMath::Max(deltaEff,TMath::Abs(parEff[sortIdx[ipoint]]-parEff[sortIdx[ipoint-2]]));
      if ( deltaEff < minEffVariation ) {
        AliWarning(Form("Efficiency variation for %s is %g => consider uniform",trendGraph->GetName(),deltaEff));
        continue;
      }
      Int_t runChangeIdx = TMath::Nint(parRunIdx[sortIdx[ipoint]]);
      AliDebug(1,Form("Change run: %s => %g => %i %i",trendGraph->GetName(),parRunIdx[sortIdx[ipoint]],runChangeIdx,GetRunNumber(runChangeIdx)));
      runRanges[irun++] = returnIndex ? runChangeIdx-1 : GetRunNumber(runChangeIdx-1);
      runRanges[irun++] = returnIndex ? runChangeIdx : GetRunNumber(runChangeIdx);
    }
    Int_t lastPt = trendGraph->GetN()-1;
    runRanges[irun++] = returnIndex ? lastPt : GetRunNumber(lastPt);
    runRanges.Set(irun);
  }
  else {
    AliWarning(Form("Fit did not converge for %s (minimum chi2 %g for step %i)",trendGraph->GetName(),minNormChi2,minNormChi2Step));
  }
  return runRanges;
}

//________________________________________________________________________
TString AliMTRChEffAnalysis::GetIdentifier ( AliTrigChEffOutput* trigOut, TObjArray* condition, Int_t itype, Int_t icount, Int_t ichamber ) const
{
  /// Get histogram identifier
//  AliTrigChEffOutput* trigOut = static_cast<AliTrigChEffOutput*>(fOutputs->At(0));
//  if ( ! trigOut ) {
//    AliWarning("No output object found!");
//    delete arr;
//    return "";
//  }
  TString identifier = Form("/%s/%s/%s/",condition->At(0)->GetName(),condition->At(1)->GetName(),condition->At(2)->GetName());
  Int_t itrackSel = static_cast<TObjString*>(condition->At(3))->String().Atoi();
  Int_t imatch = static_cast<TObjString*>(condition->At(4))->String().Atoi();
  Int_t imethod = static_cast<TObjString*>(condition->At(5))->String().Atoi();
  TString objName = trigOut->GetHistoName(itype, icount, ichamber, itrackSel, imatch, imethod);
  identifier += objName;
  return identifier;
}

//________________________________________________________________________
Int_t AliMTRChEffAnalysis::GetIndexFromRun ( UInt_t runNumber ) const
{
  /// Get index from run number
  for ( Int_t ipt=0; ipt<fOutputs->GetEntriesFast(); ipt++ ) {
    UInt_t run = fOutputs->UncheckedAt(ipt)->GetUniqueID();
    if ( run == runNumber ) return ipt;
  }
  return -1;
}

//________________________________________________________________________
TGraphAsymmErrors* AliMTRChEffAnalysis::GetOutliers ( TGraphAsymmErrors* graph, Double_t maxNsigmas ) const
{
  /// Get outliers
//  TGraphAsymmErrors* outliers = new TGraphAsymmErrors();
  TGraphAsymmErrors* outliers = new TGraphAsymmErrors(*graph);
  outliers->SetHistogram(static_cast<TH1F*>(graph->GetHistogram()->Clone(Form("%s_outliers",graph->GetHistogram()->GetName()))));
//  outliers->SetLineColor(graph->GetLineColor()+1);
  if ( graph->GetListOfFunctions()->GetEntries() == 0 ) graph->Fit("pol0","Q0");
  TF1* func = static_cast<TF1*>(graph->GetListOfFunctions()->At(0));

  Double_t xpt, ypt;
  Int_t nremoved = 0;
  for ( Int_t ipt=0; ipt<graph->GetN(); ipt++ ) {
    graph->GetPoint(ipt,xpt,ypt);
    Double_t diff = ypt - func->Eval(xpt);
    Double_t err = ( diff > 0. ) ? graph->GetErrorYlow(ipt) : graph->GetErrorYhigh(ipt);
    if ( err < 0. || TMath::Abs(diff)/err > maxNsigmas ) continue;
    outliers->RemovePoint(ipt-nremoved);
    nremoved++;
//    outliers->SetPoint(iopt,xpt,ypt);
  }
  return outliers;
}

//________________________________________________________________________
Int_t AliMTRChEffAnalysis::GetRunNumber ( Int_t ipt ) const
{
  /// Get run number from graph
  if ( ipt < 0 || ipt >= fOutputs->GetEntriesFast() ) return -1;
  return fOutputs->UncheckedAt(ipt)->GetUniqueID();
}

//________________________________________________________________________
TList* AliMTRChEffAnalysis::GetRunList ( const char* runList ) const
{
  /// Get run number in string
  TList* rl = new TList;
  rl->SetOwner();
  TString sRunList = gSystem->ExpandPathName(runList);
  if ( gSystem->AccessPathName(sRunList) || sRunList.EndsWith(".root") ) {
    sRunList.ReplaceAll(","," ");
    if ( sRunList.IsDigit() ) {
      TObjArray* arr = sRunList.Tokenize(" ");
      for ( Int_t iarr=0; iarr<arr->GetEntries(); iarr++ ) {
        rl->Add(new TObjString(arr->At(iarr)->GetName()));
      }
      delete arr;
    }
  }
  else {
    ifstream inFile(sRunList.Data());
    TString currLine = "";
    while ( ! inFile.eof() ) {
      currLine.ReadLine(inFile);
      TString currRun = AliAnalysisMuonUtility::GetRunNumberAsString(currLine);
      if ( ! currRun.IsNull() ) rl->Add(new TObjString(currRun));
    }
    inFile.close();
  }
  rl->Sort();
  return rl;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::GetShortConditionTitles ( AliTrigChEffOutput* trigOut, TObjArray& condTitles ) const
{
  /// Get short condition titles

  Int_t nConditions = fConditions->GetEntriesFast();
  TObjArray* refCondition = 0x0;

  for ( Int_t icond=0; icond<nConditions; icond++ ) {
    TObjArray* condition = static_cast<TObjArray*>(fConditions->UncheckedAt(icond));
    TString title = "";
    if ( icond == 0 ) {
      refCondition = condition;
      title = condition->GetName();
    }
    for ( Int_t ic=0; ic<condition->GetEntriesFast(); ic++ ) {
      TString currCond = static_cast<TObjString*>(condition->UncheckedAt(ic))->String();
      TString refCond = static_cast<TObjString*>(refCondition->UncheckedAt(ic))->String();
      if ( currCond == refCond ) continue;
      TString add = currCond;
      if ( ic == 3 ) add = trigOut->GetHistoName(-1,-1,-1,currCond.Atoi(),-1,-1);
      else if ( ic == 4 ) add = trigOut->GetHistoName(-1,-1,-1,-1,currCond.Atoi(),-1);
      else if ( ic == 5 ) add = trigOut->GetHistoName(-1,-1,-1,-1,-1,currCond.Atoi());
      title += Form("_%s",add.Data());
    }
    title.Remove(TString::kLeading,'_');
    title.ReplaceAll(",","|");
    condTitles.AddAt(new TObjString(title),icond);
  }

  return kTRUE;
}

//________________________________________________________________________
TH1* AliMTRChEffAnalysis::GetSum ( AliTrigChEffOutput* trigOut, TObjArray* condition, Int_t itype, Int_t icount, Int_t ichamber ) const
{
  /// Get sum histogram
  Int_t itrackSel = static_cast<TObjString*>(condition->At(3))->String().Atoi();
  Int_t imatch = static_cast<TObjString*>(condition->At(4))->String().Atoi();
  Int_t imethod = static_cast<TObjString*>(condition->At(5))->String().Atoi();
  TString objName = trigOut->GetHistoName(itype, icount, ichamber, itrackSel, imatch, imethod);
  return static_cast<TH1*>(trigOut->GetSum(condition->At(0)->GetName(),condition->At(1)->GetName(),condition->At(2)->GetName(),objName));
}

//________________________________________________________________________
TH1* AliMTRChEffAnalysis::GetTrend ( Int_t itype, Int_t icount, Int_t ichamber, Int_t idetelem ) const
{
  /// Get trending histogram
  if ( itype == AliTrigChEffOutput::kHchamberEff ) {
    if ( idetelem < 0 && ichamber >=0 ) idetelem = 11+ichamber;
  }
  TH1* outHisto = 0x0;

  TIter next(fOutputs);
  TList* effHistoList = 0x0;
  Int_t ibin = 0;
  while ( (effHistoList = static_cast<TList*>(next()) ) ) {
    ibin++;
//    TString identifier = GetIdentifier(trigOut,condition,itype,icount,ichamber);
    if ( ! outHisto ) {
//      TString outName = identifier;
      TString outName = Form("histo_type%i_count%i_ch%i_",itype,icount,11+ichamber);
//      outName.ReplaceAll("/","_");
      outName += Form("%i_trend",idetelem);
      outHisto = new TH1D(outName.Data(),outName.Data(),fOutputs->GetEntriesFast(),0.,(Double_t)fOutputs->GetEntriesFast());
      outHisto->SetDirectory(0);
      outHisto->GetXaxis()->SetTitle("Run num.");
    }
    UInt_t run = effHistoList->GetUniqueID();
    outHisto->GetXaxis()->SetBinLabel(ibin,Form("%u",run));
    TH1* histo = GetHisto(effHistoList,itype,icount,ichamber);
    Int_t currBin = histo->GetXaxis()->FindBin(idetelem);
    outHisto->SetBinContent(ibin,histo->GetBinContent(currBin));
    outHisto->SetBinError(ibin,histo->GetBinError(currBin));
  }
  if ( outHisto ) outHisto->GetXaxis()->LabelsOption("v");
  return outHisto;
}

//________________________________________________________________________
TGraphAsymmErrors* AliMTRChEffAnalysis::GetTrendEff ( Int_t itype, Int_t icount, Int_t ichamber, Int_t idetelem ) const
{
  /// Get trending histogram
  if ( icount == AliTrigChEffOutput::kAllTracks ) {
    AliWarning("Chose either bending plane, non-bending plane or both planes");
    return NULL;
  }
  TH1* histoNum = GetTrend(itype,icount,ichamber,idetelem);
  TH1* histoDen = GetTrend(itype,AliTrigChEffOutput::kAllTracks,ichamber,idetelem);
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(histoNum,histoDen,"e0");
  histoNum->Reset();
  histoNum->SetStats(kFALSE);
  graph->SetHistogram(static_cast<TH1F*>(histoNum));
  graph->GetYaxis()->SetTitle("Efficiency");
  graph->SetMarkerSize(0.5);
//  for ( Int_t ibin=1; ibin<=histoNum->GetXaxis()->GetNbins(); ibin++ ) {
//    graph->GetXaxis()->SetBinLabel(ibin,histoNum->GetXaxis()->GetBinLabel(ibin));
//  }
//  delete histoNum;
  delete histoDen;
  return graph;
}

//________________________________________________________________________
Double_t AliMTRChEffAnalysis::GetThreeOfFour ( TArrayD eff, TArrayD effErr, Double_t &probErr ) const
{
  /// Get probability of firing 3 chambers out of 4
  Double_t binomialEff = 0.;
  Double_t sumErr2 = 0.;
  for ( Int_t jch=-1; jch<4; jch++ ) {
    Double_t prodEff = 1.;
    Double_t prodErr2 = 0.;
    for ( Int_t ich=0; ich<4; ich++ ) {
      Double_t currEff = ( ich == jch ) ? 1.-eff[ich] : eff[ich];
      prodEff *= currEff;
      Double_t relErr = ( currEff>0. ) ? effErr[ich]/currEff : 0.;
      prodErr2 += relErr*relErr;
    }
    binomialEff += prodEff;
    sumErr2 += prodEff*prodEff*prodErr2;
  }
  probErr = TMath::Sqrt(sumErr2);
  return binomialEff;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::HasMergedResults () const
{
  /// Check if merged results are present
  if ( ! fMergedOutputs || fMergedOutputs->IsEmpty() ) {
    AliError("You first need to merge efficiency objects with MergeOutput");
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::MergeOutput ( TArrayI runRanges, Double_t averageStatError, Bool_t isIndex )
{
  /// Merge efficiency output

  if ( runRanges.GetSize()%2 == 1 ) {
    AliError("Run ranges expected in the form: start_period1,end_period1,start_period2,end_period2... => even number expected");
    return kFALSE;
  }

  TArrayI mergedRanges = MergeRangesForStat(runRanges, averageStatError);

  Int_t nRanges = mergedRanges.GetSize()/2;


  if ( fMergedOutputs ) fMergedOutputs->Delete();
  else {
    fMergedOutputs = new TObjArray(nRanges);
    fMergedOutputs->SetOwner();
  }

  TList* effHistoList = 0x0;
  for ( Int_t irange=0; irange<nRanges; irange++ ) {
    Int_t firstRun = mergedRanges[2*irange];
    Int_t lastRun = mergedRanges[2*irange+1];
    if ( isIndex ) {
      firstRun = GetRunNumber(firstRun);
      lastRun = GetRunNumber(lastRun);
    }
    UInt_t uFirstRun = (UInt_t)firstRun;
    UInt_t uLastRun = (UInt_t)lastRun;

    TString filename = "";

    TIter next(fOutputs);
    TFileMerger fileMerger;
    while ( (effHistoList = static_cast<TList*>(next()) ) ) {
      UInt_t run = effHistoList->GetUniqueID();
      if ( run < uFirstRun || run > uLastRun ) continue;
      filename = effHistoList->GetName();
      if ( firstRun == lastRun ) continue;
      fileMerger.AddFile(filename.Data(),kFALSE);
    }

    TString sRange = Form("%i_%i",firstRun,lastRun);
    TString mergedFilename = Form("mergedTrigEff_runs_%s.root",sRange.Data());

    if ( firstRun == lastRun ) TFile::Cp(filename.Data(),mergedFilename.Data(),kFALSE);
    else {
      fileMerger.OutputFile(mergedFilename.Data());
      fileMerger.Merge();
    }
    AliTrigChEffOutput* trigOut = new AliTrigChEffOutput(mergedFilename);
    trigOut->SetName(sRange.Data());
    fMergedOutputs->Add(trigOut);
  }
  return kTRUE;
}

//________________________________________________________________________
TArrayI AliMTRChEffAnalysis::MergeRangesForStat ( TArrayI runRanges, Double_t averageStatError, Bool_t excludePeriphericBoards ) const
{
  if ( averageStatError <= 0. || averageStatError >= 1. ) return runRanges;

  // FIXME: statstical error for binomial depends on the value of epsilon
  // for the moment let's assume an average efficiency of 0.9 (underestimated)
  // in the future one could maybe actually calculate it
  // (still we're working with average values...so the precision is limited)
  Double_t effForStatCalc = 0.9;

  Double_t averageStatNeeded = effForStatCalc*(1.-effForStatCalc)/(averageStatError*averageStatError);

  AliInfo(Form("Average statistics needed to reach precision of %g : %g",averageStatError,averageStatNeeded));

  Int_t nRanges = runRanges.GetSize()/2;

  TArrayD averageStat(nRanges);
  Double_t fullStat = 0.;
  for ( Int_t irange=0; irange<nRanges; irange++ ) {
    averageStat[irange] = GetAverageStat(runRanges[2*irange],runRanges[2*irange+1],AliTrigChEffOutput::kHboardEff,excludePeriphericBoards);
    fullStat += averageStat[irange];
  }

  TArrayI mergedRanges(runRanges.GetSize());
  Int_t imerged = 0;
  mergedRanges[imerged++] = runRanges[0];
  Double_t mergedAverageStat = 0., remainingStat = fullStat;
  for ( Int_t irange=0; irange<nRanges; irange++ ) {
    Int_t istart = 2*irange;
    Int_t iend = istart+1;
    mergedAverageStat += averageStat[irange];
    remainingStat -= averageStat[irange];

    AliInfo(Form("%i - %i => stat %g",runRanges[2*irange],runRanges[2*irange+1],averageStat[irange]));

    if ( ( mergedAverageStat >= averageStatNeeded && remainingStat >= averageStatNeeded ) || iend == runRanges.GetSize()-1 ) {
      mergedRanges[imerged++] = runRanges[iend];
      AliInfo(Form("  merged range %i - %i => stat %g",mergedRanges[imerged-2],mergedRanges[imerged-1],mergedAverageStat));
      mergedAverageStat = 0.;
      Int_t nextRun = iend+1;
      if ( nextRun < runRanges.GetSize() ) mergedRanges[imerged++] = runRanges[nextRun];
    }
  }
  mergedRanges.Set(imerged);
  return mergedRanges;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::SetCondition ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod, Bool_t isBasic )
{
  /// Select the conditions for the efficiency estimation
  if ( ! fConditions ) {
    fConditions = new TObjArray();
    fConditions->SetOwner();
  }

  TString sCondition = Form("%s;%s;%s;%i;%i;%i",physSel,trigClassName,centrality,itrackSel,imatch,imethod);
  TObjArray* foundCondition = static_cast<TObjArray*>(fConditions->FindObject(sCondition.Data()));
  TObjArray* basicCondition = static_cast<TObjArray*>(fConditions->At(0));

  TObjArray* addCondition = 0x0;

  if ( foundCondition ) {
    if ( isBasic ) {
      if ( foundCondition == basicCondition ) return kFALSE;
      else {
        fConditions->Remove(foundCondition);
        addCondition = foundCondition;
      }
    }
    else {
      AliInfo("Systematic condition already added");
      return kFALSE;
    }
  }
  else {
    addCondition = sCondition.Tokenize(";");
    addCondition->SetName(sCondition.Data());
  }


  if ( isBasic ) {
    if ( basicCondition ) {
      AliInfo(Form("Changing current eff. condition: %s\n",basicCondition->GetName()));
      fConditions->Remove(basicCondition);
      delete basicCondition;
    }
    fConditions->AddAt(addCondition,0);
    fConditions->Compress();

    // Update outputs
    if ( fOutputs ) {
      TObjArray* tmpOutputs = static_cast<TObjArray*>(fOutputs->Clone());
      fOutputs->Delete();
      TIter next(tmpOutputs);
      TList* effHistiList = 0x0;
      while ( (effHistiList = static_cast<TList*>(next())) ) {
        AddToList(effHistiList->GetName(), fOutputs->GetName());
      }
      delete tmpOutputs;
    }
  }
  else fConditions->Add(addCondition);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::SetDefaultEffConditions ()
{
  /// Set default conditions for efficiency estimation
  return SetEffConditions("PhysSelPass","ANY","-5_105",AliTrigChEffOutput::kSelectTrack,AliTrigChEffOutput::kMatchApt,AliTrigChEffOutput::kEffFromTrack);
}


//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::SetEffConditions ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod )
{
  /// Select the conditions for the efficiency estimation
  return SetCondition(physSel, trigClassName, centrality, itrackSel, imatch, imethod, kTRUE);
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::SetOutList ( const char *localFileList, const char *outputName )
{
  /// Initialize output list

  if ( ! fConditions ) SetDefaultEffConditions();

  TString filename(localFileList);
  gSystem->ExpandPathName(filename);
  if ( gSystem->AccessPathName(filename.Data()) ) {
    AliWarning(Form("Cannot find %s",filename.Data()));
    return kFALSE;
  }
  if ( filename.EndsWith(".root") ) return AddToList(filename.Data(),outputName);

  Bool_t isOk = kTRUE;
  ifstream inFile(filename.Data());
  TString currLine = "";
  while ( ! inFile.eof() ) {
    currLine.ReadLine(inFile);
    if ( currLine.IsNull() ) continue;
    if ( gSystem->AccessPathName(currLine.Data()) ) continue;
    if ( ! AddToList(currLine.Data(), outputName) ) isOk = kFALSE;
  }
  inFile.close();

  return isOk;
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::SetResultsFromGrid ( const char *runList, const char *path, const char *pattern, const char* localFileList, const char* outDir, const char *directory, const char* outputName )
{
  /// Search for results on grid
  CopyLocally(runList,path,pattern,localFileList,outDir,directory);
  return SetOutList(localFileList,outputName);
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::SetResultsFromWeb ( const char *runList, const char *path, const char* localFileList, const char* outDir, const char *directory, const char* outputName )
{
  /// Search for results on grid
  CopyLocally(runList,path,"",localFileList,outDir,directory);
  return SetOutList(localFileList,outputName);
}

//________________________________________________________________________
Bool_t AliMTRChEffAnalysis::WriteMergedToOCDB ( const char* outputCDB ) const
{
  /// Create the OCDB objects
  if ( ! HasMergedResults() ) return kFALSE;
  TString outCDBdir(outputCDB);
  if ( outCDBdir.Contains("://") && ! outCDBdir.BeginsWith("local") ) {
    AliError("A local directory is required to store the information");
    return kFALSE;
  }
  outCDBdir.ReplaceAll("local://","");
  if ( gSystem->AccessPathName(outCDBdir.Data()) ) {
    if ( ! ExecCommand(Form("mkdir -p %s",outCDBdir.Data()), kTRUE) ) {
      AliError(Form("Failed to create directory %s",outCDBdir.Data()));
      return kFALSE;
    }
  }
  TString outCDB = outCDBdir;
  outCDB.Prepend("local://");

  AliCDBManager* mgr = AliCDBManager::Instance();
  if ( ! mgr->GetDefaultStorage() ) mgr->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  TString trigEffCDBdir = "MUON/Calib/TriggerEfficiency";

  mgr->SetSpecificStorage(trigEffCDBdir.Data(),outCDB.Data());

  AliCDBStorage* specificStorage = mgr->GetSpecificStorage(trigEffCDBdir.Data());

  TObjArray* condition = static_cast<TObjArray*>(fConditions->At(0));

  AliTrigChEffOutput* trigOut = 0x0;
  TIter next(fMergedOutputs);
  while ( (trigOut = static_cast<AliTrigChEffOutput*>(next())) ) {
    TList* effList = GetEffHistoList(trigOut,condition);
    TString runRange = trigOut->GetName();

//    // Write full merged output in file
//    TFile* file = TFile::Open(Form("mergedTrigEff_runs_%s.root",runRange.Data()),"RECREATE");
//    effList->Write("triggerChamberEff",TObject::kSingleKey);
//    trigOut->GetOutput()->Write(trigOut->GetOutput()->GetName(),TObject::kSingleKey);
//    file->Close();
//    delete file;


    // Write OCDB object
    TObjArray* currRange = runRange.Tokenize("_");
    Int_t firstRun = static_cast<TObjString*>(currRange->UncheckedAt(0))->String().Atoi();
    Int_t lastRun = static_cast<TObjString*>(currRange->UncheckedAt(1))->String().Atoi();
    delete currRange;

    // If an object is already there, ask to remove it or keep it
    for ( Int_t irun=0; irun<2; irun++ ) {
      Int_t runnr = ( irun == 0 ) ? firstRun : lastRun;
      specificStorage->QueryCDB(runnr);
      TObjArray* allIdsForRun = specificStorage->GetQueryCDBList();
      TIter nextId(allIdsForRun);
      AliCDBId* id = 0x0;
      while ((id = dynamic_cast<AliCDBId*>(nextId()))) {
        TString path(id->GetPath());
        Int_t foundFirst = id->GetFirstRun();
        Int_t foundLast = id->GetLastRun();
        Int_t version = id->GetVersion();
        Int_t subversion = id->GetSubVersion();
        TString locPath = Form("%s/%s/Run%d_%d_v%d_s%d.root",outCDBdir.Data(),path.Data(),foundFirst,foundLast,version,subversion);
        if ( ! gSystem->AccessPathName(locPath.Data()) ) {
          ExecCommand(Form("rm %s",locPath.Data()), kTRUE);
        }
      }
    }

    // Save the CDB object in the specific storage
    AliMUONTriggerEfficiencyCells* effMap = new AliMUONTriggerEfficiencyCells(effList);
    AliMUONCDB::WriteToCDB(effMap, "MUON/Calib/TriggerEfficiency", firstRun, lastRun, "Measured efficiencies");
    delete effMap; // CAVEAT: effMap is owner of efflist
  }
  return kTRUE;
}

//___________________________________________________________________________
//Bool_t AliMTRChEffAnalysis::WriteMerged() const
//{
//  /// Write merged efficiencies on file
//  if ( ! HasMergedResults ) return kFALSE;
//
//  TObjArray* condition = static_cast<TObjArray*>(fConditions->At(0));
//
//  AliTrigChEffOutput* trigOut = 0x0;
//  TIter next(fMergedOutputs);
//  while ( (trigOut = static_cast<AliTrigChEffOutput*>(next())) ) {
//    TFile* file = TFile::Open(Form("mergedTrigEff_runs_%s.root",trigOut->GetName()),"RECREATE");
//    TList* effList = GetEffHistoList(trigOut,condition);
//    effList->Write("triggerChamberEff",TObject::kSingleKey);
//    trigOut->GetOutput()->Write(trigOut->GetOutput()->GetName(),TObject::kSingleKey);
//    file->Close();
//    delete file;
//    delete effList;
//  }
//  return kTRUE;
//}
