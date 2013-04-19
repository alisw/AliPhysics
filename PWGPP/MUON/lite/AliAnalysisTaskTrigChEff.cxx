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

/* $Id$ */

#include "AliAnalysisTaskTrigChEff.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TFile.h"

// STEER includes
#include "AliVParticle.h"
#include "AliAODTrack.h"
#include "AliESDMuonTrack.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// PWG3 includes
#include "AliVAnalysisMuon.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskTrigChEff) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskTrigChEff::AliAnalysisTaskTrigChEff() :
  AliVAnalysisMuon(),
  fTrackSelKeys(0x0),
  fCountTypeKeys(0x0),
  fHistoTypeKeys(0x0),
  fEffMethodKeys(0x0),
  fUseGhosts(kFALSE),
  fList(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskTrigChEff::AliAnalysisTaskTrigChEff(const char *name, const AliMuonTrackCuts& cuts) :
  AliVAnalysisMuon(name, cuts),
  fTrackSelKeys(0x0),
  fCountTypeKeys(0x0),
  fHistoTypeKeys(0x0),
  fEffMethodKeys(0x0),
  fUseGhosts(kFALSE),
  fList(0x0)
{
  //
  /// Constructor.
  //

  InitLocalKeys();
  
  DefineOutput(2, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskTrigChEff::~AliAnalysisTaskTrigChEff()
{
  //
  /// Destructor
  //
  delete fTrackSelKeys;
  delete fCountTypeKeys;
  delete fHistoTypeKeys;
  delete fEffMethodKeys;
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fList;
  }
}

//___________________________________________________________________________
TList* AliAnalysisTaskTrigChEff::GetEffHistoList(TString physSel, TString trigClassNames, TString centrality, TString trackSelection)
{
  /// Get the list of efficiency objects by merging the 
  // results from the histogram collection
  
  TList* outList = new TList();
  outList->SetOwner();
  FillEffHistoList(physSel, trigClassNames, centrality, trackSelection, outList);
  return outList;
}

//___________________________________________________________________________
Bool_t AliAnalysisTaskTrigChEff::FillEffHistoList(TString physSel, TString trigClassNames, TString centrality, TString trackSelection, TList* outList)
{
  /// Fill the list of objects for the efficiency calculation
  /// merging the splitted output of the fHistogramCollection
  /// The obtained list can be converted in the efficiency map used in simulations
  /// in a backward compatible way
  
  if ( ! fMergeableCollection ) return kFALSE;
  if ( ! outList ) outList = fList;
  TString histoName = "";
  TString histoPattern = "";
  TH1* histo = 0x0;
  Bool_t isOk = kTRUE;
  for ( Int_t icount=0; icount<kNcounts; ++icount ) {
    histoName = GetHistoName(kHchamberEff, icount, -1, -1, -1);
    histoPattern = Form("%s%s", histoName.Data(), trackSelection.Data());
    histo = (TH1*)GetSum(physSel, trigClassNames, centrality, histoPattern);
    if ( histo ) {
      histo->SetName(histoName.Data());
      histo->SetTitle(histoName.Data());
      histo->SetDirectory(0);
      outList->Add(histo);
    }
    else isOk = kFALSE;
  }
  for ( Int_t icount=0; icount<kNcounts; ++icount ) {
    for ( Int_t ich=0; ich<4; ++ich ) {
      histoName = GetHistoName(kHslatEff, icount, ich, -1, -1);
      histoPattern = Form("%s%s", histoName.Data(), trackSelection.Data());
      histo = (TH1*)GetSum(physSel, trigClassNames, centrality, histoPattern);
      if ( histo ) {
        histo->SetName(histoName.Data());
        histo->SetTitle(histoName.Data());
        histo->SetDirectory(0);
        outList->Add(histo);
      }
      else isOk = kFALSE;
    }
  }
  for ( Int_t icount=0; icount<kNcounts; ++icount ) {
    for ( Int_t ich=0; ich<4; ++ich ) {
      histoName = GetHistoName(kHboardEff, icount, ich, -1, -1);
      histoPattern = Form("%s%s", histoName.Data(), trackSelection.Data());
      histo = (TH1*)GetSum(physSel, trigClassNames, centrality, histoPattern);
      if ( histo ) {
        histo->SetName(histoName.Data());
        histo->SetTitle(histoName.Data());
        histo->SetDirectory(0);
        outList->Add(histo);
      }
      else isOk = kFALSE;
    }
  }
  
  histoName = GetHistoName(kHcheckBoard, -1, -1, -1, -1);
  histoPattern = Form("%s%s", histoName.Data(), trackSelection.Data());
  histo = (TH1*)GetSum(physSel, trigClassNames, centrality, histoPattern);
  if ( histo ) {
    histo->SetName(histoName.Data());
    histo->SetTitle(histoName.Data());
    histo->SetDirectory(0);
    outList->Add(histo);
  }
  else isOk = kFALSE;
  
  return isOk;
}

//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::FinishTaskOutput()
{
  //
  /// Merge Apt, Lpt and Hpt
  /// Fill the final efficiency object (for backward compatibility)
  //

  TString histoName = "";
  for ( Int_t isel=0; isel<kNselections; ++isel ) {
    for ( Int_t itrig=0; itrig<GetAllSelectedTrigClasses()->GetEntries(); ++itrig ) {
      for ( Int_t icent=1; icent<=GetCentralityClasses()->GetNbins(); ++icent ) {
        for ( Int_t imethod=0; imethod<kNeffMethods; ++imethod ) {
          for ( Int_t itype=0; itype<kNhistoTypes; ++itype ) {
            for ( Int_t icount=-1; icount<kNcounts; ++icount ) {
              for ( Int_t ich=-1; ich<4; ++ich ) {
                for ( Int_t imatch=kMatchApt; imatch<kMatchHpt; ++imatch ) {
                  TH1* histo = 0x0;
                  for ( Int_t jmatch=imatch+1; jmatch<=kMatchHpt; ++jmatch ) {
                    histoName = GetHistoName(itype, icount, ich, jmatch, imethod);
                    TH1* histoAdd = (TH1*)fMergeableCollection->GetObject(Form("/%s/%s/%s/",fPhysSelKeys->At(isel)->GetName(), GetAllSelectedTrigClasses()->At(itrig)->GetName(), GetCentralityClasses()->GetBinLabel(icent)), histoName);
                    if ( ! histoAdd ) continue;
                    histoName = GetHistoName(itype, icount, ich, imatch, imethod);
                    if ( ! histo ) histo = (TH1*)GetMergeableObject(fPhysSelKeys->At(isel)->GetName(), GetAllSelectedTrigClasses()->At(itrig)->GetName(), GetCentralityClasses()->GetBinLabel(icent), histoName);
                    AliDebug(2,Form("Adding %s (%g) to %s (%g)", histoAdd->GetName(), histoAdd->Integral(), histo->GetName(), histo->Integral()));
                    histo->Add(histoAdd);
                  } // loop on higher pt matching
                } // loop on match trigger
              } // loop on chamber
            } // loop on count type
          } // loop on histo type
        } // loop on eff method
      } // loop on centrality
    } // loop on trigger classes
  } // loop on physics selection

  TString physSel = fPhysSelKeys->At(kPhysSelPass)->GetName();
  TString trigClass = "ANY";
  TString centrality = "all";
  TString histoPattern = Form("%s%s",fTrackSelKeys->At(kMatchApt)->GetName(),fEffMethodKeys->At(kEffFromTrack)->GetName());
  
  FillEffHistoList(physSel, trigClass, centrality, histoPattern, fList);

  AliVAnalysisMuon::FinishTaskOutput();
}


//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::InitLocalKeys()
{
  //
  /// Initialyze objects
  //
  
  TString trackSelNames = "MatchNoTrackSelCut MatchApt MatchLpt MatchHpt"; 
  fTrackSelKeys = trackSelNames.Tokenize(" ");
  
  TString countTypeNames = "bendPlaneCount nonBendPlaneCount bothPlanesCount allTracksCount";
  fCountTypeKeys = countTypeNames.Tokenize(" ");
  
  TString histoTypeKeys = "Chamber Slat Board checkRejectedBoard";
  fHistoTypeKeys = histoTypeKeys.Tokenize(" ");
  
  TString effMethodKeys = "FromTrk FromTrg";
  fEffMethodKeys = effMethodKeys.Tokenize(" ");
}

//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::MyUserCreateOutputObjects()
{
  //
  /// Create prototype objects for mergeable collection
  //
  
  const Char_t* yAxisTitle = "counts";
  
  Int_t nChamberBins = 4;
  Float_t chamberLow = 11.-0.5, chamberHigh = (Float_t)(nChamberBins)+11.-0.5;
  const Char_t* chamberName = "chamber";
  
  Int_t nSlatBins = 18;
  Float_t slatLow = 0.-0.5, slatHigh = (Float_t)nSlatBins-0.5;
  const Char_t* slatName = "slat";
  
  Int_t nBoardBins = 234;
  Float_t boardLow = 1.-0.5, boardHigh = (Float_t)nBoardBins+1.-0.5;
  const Char_t* boardName = "board";
  
  TString baseName, histoName, histoTitle;
  
  TH1* histo;
  TH2F* histo2D;
  
  for ( Int_t imethod=0; imethod<kNeffMethods; ++imethod ) {
    for ( Int_t itrackSel = 0; itrackSel<kNtrackSel; ++itrackSel ) {
      for ( Int_t icount=0; icount<kNcounts; ++icount ) {
        histoName = GetHistoName(kHchamberEff, icount, -1, itrackSel, imethod);
        histo = new TH1F(histoName, histoName,
                         nChamberBins, chamberLow, chamberHigh);
        histo->GetXaxis()->SetTitle(chamberName);
        histo->GetYaxis()->SetTitle(yAxisTitle);
        AddObjectToCollection(histo);
      } // loop on counts
      
      for ( Int_t icount=0; icount<kNcounts; ++icount ) {
        for ( Int_t ich=0; ich<4; ++ich ) {
          histoName = GetHistoName(kHslatEff, icount, ich, itrackSel, imethod);
          histo = new TH1F(histoName, histoName,
                           nSlatBins, slatLow, slatHigh);
          histo->GetXaxis()->SetTitle(slatName);
          histo->GetYaxis()->SetTitle(yAxisTitle);
          AddObjectToCollection(histo);
        } // loop on chamber
      } // loop on counts
      
      for ( Int_t icount=0; icount<kNcounts; ++icount ) {
        for ( Int_t ich=0; ich<4; ++ich ) {
          histoName = GetHistoName(kHboardEff, icount, ich, itrackSel, imethod);
          histo = new TH1F(histoName, histoName,
                           nBoardBins, boardLow, boardHigh);
          histo->GetXaxis()->SetTitle(boardName);
          histo->GetYaxis()->SetTitle(yAxisTitle);
          AddObjectToCollection(histo);
        } // loop on chamber
      } // loop on counts
      
      histoName = GetHistoName(kHcheckBoard, -1, -1, itrackSel, imethod);
      histo2D = new TH2F(histoName.Data(), "Rejected tracks motivation", 
                         5, 20.5, 25.5, nBoardBins, boardLow, boardHigh);
      histo2D->GetXaxis()->SetBinLabel(1,"Many pads");
      histo2D->GetXaxis()->SetBinLabel(2,"Few pads");
      histo2D->GetXaxis()->SetBinLabel(3,"Outside geom");
      histo2D->GetXaxis()->SetBinLabel(4,"Tracker track");
      histo2D->GetXaxis()->SetBinLabel(5,"Masked board");
      histo2D->GetYaxis()->SetTitle(boardName);
      AddObjectToCollection(histo2D);
    } // loop on track selection
  } // loop on eff method

  fMuonTrackCuts->Print();
  
  fList = new TList();
  fList->SetOwner();
  PostData(2, fList);

}

//___________________________________________________________________________
TString AliAnalysisTaskTrigChEff::GetHistoName(Int_t itype, Int_t icount, Int_t ichamber, Int_t itrackSel, Int_t imethod)
{
  /// Get histogram index
  TString histoName = "";
  if ( itype < kHcheckBoard && icount >= 0 ) histoName += ((TObjString*)fCountTypeKeys->At(icount))->GetString();
  histoName += ((TObjString*)fHistoTypeKeys->At(itype))->GetString();
  if ( ichamber >= 0 ) histoName += Form("Ch%i", 11+ichamber);
  if ( itrackSel >= 0 ) histoName += ((TObjString*)fTrackSelKeys->At(itrackSel))->GetString();
  if ( imethod >= 0 ) histoName += ((TObjString*)fEffMethodKeys->At(imethod))->GetString();
  return histoName;
}

//________________________________________________________________________
void AliAnalysisTaskTrigChEff::ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality)
{
  //
  /// Fill histogram
  //

  Int_t slat = 0, board = 0;
  UInt_t pattern = 0;
  TString histoName = "";

  TArrayI othersEfficient(4);

  AliVParticle* track = 0x0;

  Int_t itrackSel = -1;

  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(InputEvent());
  for ( Int_t itrack = 0; itrack < nTracks; ++itrack ) {
    track = AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    Bool_t matchTracker = AliAnalysisMuonUtility::IsMuonTrack(track);
    if ( ! matchTracker && ! fUseGhosts ) continue;

    Int_t matchTrig = AliAnalysisMuonUtility::GetMatchTrigger(track);
    itrackSel = matchTrig;
    UInt_t selection = fMuonTrackCuts->GetSelectionMask(track);
    Bool_t isSelected = ( ( selection & fMuonTrackCuts->GetFilterMask() ) == fMuonTrackCuts->GetFilterMask() );
    if ( matchTrig == 1 && ( ( selection & AliMuonTrackCuts::kMuMatchApt ) == 0 ) ) isSelected = kFALSE;
    if ( matchTrig == 2 && ( ( selection & AliMuonTrackCuts::kMuMatchLpt ) == 0 ) ) isSelected = kFALSE;
    if ( matchTrig == 3 && ( ( selection & AliMuonTrackCuts::kMuMatchHpt ) == 0 ) ) isSelected = kFALSE;
    if ( ! isSelected ) itrackSel = kNoSelCutApt;

    for ( Int_t imethod=0; imethod<kNeffMethods; ++imethod ) {
      if ( imethod == kEffFromTrack ) {
        if ( track->P() < 10. ) continue;
        pattern = AliAnalysisMuonUtility::GetMUONTrigHitsMapTrk(track);
        board = AliESDMuonTrack::GetCrossedBoard(pattern);
      }
      else {
        pattern = AliAnalysisMuonUtility::GetMUONTrigHitsMapTrg(track);
        board = ( AliAnalysisMuonUtility::IsAODEvent(InputEvent()) ) ? AliESDMuonTrack::GetCrossedBoard(pattern) : ((AliESDMuonTrack*)track)->LoCircuit();
      }
      
      Int_t effFlag = AliESDMuonTrack::GetEffFlag(pattern);
      
      if ( effFlag < AliESDMuonTrack::kChEff ) {
        for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
          TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
          
          histoName = GetHistoName(kHcheckBoard, -1, -1, itrackSel, imethod);
          ((TH2F*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(AliESDMuonTrack::GetSlatOrInfo(pattern), board);
        }
        continue; // Track not good for efficiency calculation
      }
      
      othersEfficient.Reset(1);
      for ( Int_t cath=0; cath<2; ++cath ) {
        for ( Int_t ich=0; ich<4; ++ich ) {
          if ( ! AliESDMuonTrack::IsChamberHit(pattern, cath, ich) ) {
            for ( Int_t jch=0; jch<4; jch++ ) {
              if ( jch != ich ) othersEfficient[jch] = 0;
            } // loop on other chambers
            break;
          } // if chamber not efficient
        } // loop on chambers
      } // loop on cathodes
      
      Bool_t rejectTrack = kTRUE;
      for ( Int_t ich=0; ich<4; ++ich ) {
        if ( othersEfficient[ich] > 0 ) {
          rejectTrack = kFALSE;
          break;
        }
      }
      
      if ( rejectTrack ) continue;
      
      slat = AliESDMuonTrack::GetSlatOrInfo(pattern);
      
      for ( Int_t ich=0; ich<4; ++ich ) {
        if ( ! othersEfficient[ich] )
          continue; // Reject track if the info of the chamber under study 
        // is necessary to create the track itself
        
        Int_t iChamber = 11 + ich;
        
        Bool_t hitsBend = AliESDMuonTrack::IsChamberHit(pattern, 0, ich);
        Bool_t hitsNonBend = AliESDMuonTrack::IsChamberHit(pattern, 1, ich);
        
        Bool_t fillHisto[kNcounts] = {
          hitsBend,
          hitsNonBend,
          ( hitsBend && hitsNonBend ),
          kTRUE
        };
        
        for (Int_t icount=0; icount<kNcounts; ++icount){
          if ( ! fillHisto[icount] ) continue;
          for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
            TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
            
            histoName = GetHistoName(kHchamberEff, icount, -1, itrackSel, imethod);
            ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(iChamber);
            
            if ( effFlag < AliESDMuonTrack::kSlatEff ) continue; // Track crossed different slats
            histoName = GetHistoName(kHslatEff, icount, ich, itrackSel, imethod);
            ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(slat);
            
            if ( effFlag < AliESDMuonTrack::kBoardEff ) continue; // Track crossed different boards
            histoName = GetHistoName(kHboardEff, icount, ich, itrackSel, imethod);
            ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(board);
          } // loop on trigger classes
        } // loop on chambers
      } // loop on count types
    } // loop on tracks
  } // loop on eff methods
  
  PostData(2, fList);
}


//________________________________________________________________________
void AliAnalysisTaskTrigChEff::Terminate(Option_t *)
{
  //
  /// Draw some histogram at the end.
  //

  AliVAnalysisMuon::Terminate("");
  
  if ( ! fMergeableCollection ) return;
  

  Int_t xshift = 100;
  Int_t yshift = 20;
  Int_t igroup1 = -1;
  Int_t igroup2 = 0;

  TObjArray* physSel =  ((TObjString*)fTerminateOptions->At(0))->GetString().Tokenize(" ");
  physSel->SetOwner();
  TObjArray* trigClasses = ((TObjString*)fTerminateOptions->At(1))->GetString().Tokenize(" ");
  trigClasses->SetOwner();
  TObjArray* centrality = ((TObjString*)fTerminateOptions->At(2))->GetString().Tokenize(" ");
  centrality->SetOwner();
  TString furtherOpt = ((TObjString*)fTerminateOptions->At(3))->GetString();

  TString currName = "";
  TObjArray* optArr = furtherOpt.Tokenize(" ");
  TObjArray trackSel, methodSel;
  trackSel.SetOwner();
  methodSel.SetOwner();
  TString outFileOpt = "";
  for ( Int_t iopt=0; iopt<optArr->GetEntries(); iopt++ ) {
    currName = optArr->At(iopt)->GetName();
    if ( currName.Contains(".root") ) outFileOpt = currName;
    else if ( currName.Contains("Match") ) trackSel.Add(new TObjString(currName));
    else if ( currName.Contains("From") ) methodSel.Add(new TObjString(currName));
  }
  delete optArr;
  
  if ( trackSel.GetEntries() == 0 ) trackSel.Add(new TObjString(fTrackSelKeys->At(kMatchApt)->GetName()));
  if ( methodSel.GetEntries() == 0 ) methodSel.Add(new TObjString(fEffMethodKeys->At(kEffFromTrack)->GetName()));

  furtherOpt.ToUpper();
  
  Int_t chosenType = ( furtherOpt.Contains("BOARD") ) ? kHboardEff : kHslatEff;

  igroup1++;
  igroup2 = 0;
  TString histoName = "", yAxisTitle = "";
  
  TH1 *num = 0x0;
  TH1 *den = 0x0;
  TGraphAsymmErrors* effGraph = 0x0;
  
  ////////////////
  // Show tests //
  ////////////////
  
  for ( Int_t icount=0; icount<kNcounts-1; ++icount ) {
    currName = Form("%s%s_can", fHistoTypeKeys->At(chosenType)->GetName(), fCountTypeKeys->At(icount)->GetName());
    TCanvas* can = new TCanvas(currName.Data(), currName.Data(), igroup1*xshift,igroup2*yshift,600,600);
    can->Divide(2,2);
    TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->SetBorderSize(1);
    for ( Int_t ich=0; ich<4; ++ich ) {
      TGraph* refGraph = 0x0;
      can->cd(ich+1);
      gPad->SetRightMargin(0.03);
      Int_t icolor = 1;
      Int_t istyle = 0;
      TString drawOpt = "AP";
      for ( Int_t isel=0; isel<physSel->GetEntries(); ++isel ) {
        for ( Int_t itrig=0; itrig<trigClasses->GetEntries(); ++itrig ) {
          for ( Int_t icent=0; icent<centrality->GetEntries(); ++icent ) {
            for ( Int_t imethodSel=0; imethodSel<methodSel.GetEntries(); ++imethodSel ) {
              for ( Int_t itrackSel=0; itrackSel<trackSel.GetEntries(); ++itrackSel ) {
                histoName = GetHistoName(chosenType, kAllTracks, ich, -1, -1); // partial name
                histoName += Form("%s%s",trackSel.At(itrackSel)->GetName(),methodSel.At(imethodSel)->GetName());
                den = (TH1*)GetSum(physSel->At(isel)->GetName(), trigClasses->At(itrig)->GetName(), centrality->At(icent)->GetName(), histoName.Data());
                if ( ! den ) continue;
                histoName = GetHistoName(chosenType, icount, ich, -1, -1); // partial name
                histoName += Form("%s%s",trackSel.At(itrackSel)->GetName(),methodSel.At(imethodSel)->GetName());
                num = (TH1*)GetSum(physSel->At(isel)->GetName(), trigClasses->At(itrig)->GetName(), centrality->At(icent)->GetName(), histoName.Data());
                if ( ! num ) continue;
                effGraph = new TGraphAsymmErrors(num, den, "e0");
                currName = Form("%s_%s_%s_%s_%s", physSel->At(isel)->GetName(), trigClasses->At(itrig)->GetName(), centrality->At(icent)->GetName(), trackSel.At(itrackSel)->GetName(), methodSel.At(imethodSel)->GetName());
                effGraph->SetTitle(currName.Data());

                Double_t ymin = 0.;
                Double_t ymax = 1.1;
                yAxisTitle = "Efficiency";

                if ( furtherOpt.Contains("DIFF") || furtherOpt.Contains("PULL") ) {
                  if ( ! refGraph ) {
                    refGraph = effGraph;
                    continue;
                  }
                  Double_t currX, currY, baseX, baseY, errYlow = 0., errYhigh = 0., newY = 0.;
                  Double_t refVal = 1., errY = 0.; 
                  for ( Int_t ipoint=0; ipoint<effGraph->GetN(); ipoint++ ) {
                    refGraph->GetPoint(ipoint, baseX, baseY);
                    effGraph->GetPoint(ipoint, currX, currY);
                    Double_t errX = effGraph->GetErrorXlow(ipoint);
                    if ( furtherOpt.Contains("DIFF") ) {
                      refVal = ( baseY > 0. ) ? baseY : 1.;
                      newY = ( currY - baseY ) / refVal;
                      Double_t errYlow1 = effGraph->GetErrorYlow(ipoint);
                      Double_t errYlow2 = refGraph->GetErrorYlow(ipoint);
                      Double_t errYhigh1 = effGraph->GetErrorYhigh(ipoint);
                      Double_t errYhigh2 = refGraph->GetErrorYhigh(ipoint);
                      errYlow = TMath::Sqrt(errYlow1*errYlow1 + errYlow2*errYlow2) / refVal;
                      errYhigh = TMath::Sqrt(errYhigh1*errYhigh1 + errYhigh2*errYhigh2) / refVal;
                      //yAxisTitle = Form("(%s - %s) / %s", effGraph->GetTitle(), refGraph->GetTitle(), refGraph->GetTitle());
                      yAxisTitle = "(eff - ref ) / ref";
                      effGraph->SetTitle(Form("Rel. diff. w.r.t. %s", refGraph->GetTitle()));
                      ymin = -0.1;
                      ymax = 0.1;
                    }
                    else if ( furtherOpt.Contains("PULL") ) {
                      errY = 0.5 * ( effGraph->GetErrorYlow(ipoint) + effGraph->GetErrorYhigh(ipoint));
                      newY = ( errY > 0. ) ? ( currY - baseY ) / errY : 0.;
                      errYlow = 1.;
                      errYhigh = 1.;
                      //yAxisTitle = Form("( %s - %s ) / err", effGraph->GetTitle(), refGraph->GetTitle());
                      yAxisTitle = "(eff - ref ) / err";
                      effGraph->SetTitle(Form("Pull w.r.t. %s", refGraph->GetTitle()));
                      ymin = -4.;
                      ymax = 4.;
                    }
                    effGraph->SetPoint(ipoint, currX, newY);
                    effGraph->SetPointError(ipoint, errX, errX, errYlow, errYhigh);
                  } // loop on points
                }
                effGraph->GetYaxis()->SetRangeUser(ymin, ymax);
                effGraph->GetYaxis()->SetTitle(yAxisTitle.Data());
                effGraph->SetLineColor(icolor);
                effGraph->SetMarkerColor(icolor);
                effGraph->SetMarkerStyle(20 + istyle);
                //effGraph->SetMarkerSize(0.3);
                icolor++;
                if ( icolor == 5 || icolor == 10 ) icolor++;
                istyle++;
                effGraph->Draw(drawOpt.Data());
                drawOpt = "P";
                if ( ich < 3 ) continue;
                leg->AddEntry(effGraph, currName.Data(), "lp");
              } // loop on match trigger
            } // loop on eff method
          } // loop on centrality
        } // loop on trigger classes
      } // loop on physics selection
    } // loop on chamber
    leg->Draw("same");
  } // loop on count type
  //} // loop on detection element type

  delete physSel;
  delete trigClasses;
  delete centrality;
  
   
  fList = dynamic_cast<TList*>(GetOutputData(2));
  if ( fList ) {
  
    ///////////////////////////
    // Show final efficiency //
    ///////////////////////////
    TString baseName[3] = {"Chamber", "RPC", "Board"};
    Int_t baseIndex[3] = {kHchamberEff, kHslatEff, kHboardEff};
    TString effName[kNcounts-1] = {"BendPlane", "NonBendPlane", "BothPlanes"};
    for ( Int_t itype=0; itype<3; itype++ ) {
      for ( Int_t icount=0; icount<kNcounts-1; icount++ ){
        TString canName = Form("efficiencyPer%s_%s",baseName[itype].Data(),effName[icount].Data());
        TCanvas* can = new TCanvas(canName.Data(),canName.Data(),10*(1+kNcounts*itype+icount),10*(1+kNcounts*itype+icount),310,310);
        can->SetFillColor(10); can->SetHighLightColor(10);
        can->SetLeftMargin(0.15); can->SetBottomMargin(0.15);  
        if ( itype > 0 )
          can->Divide(2,2);
        
        for ( Int_t ich=-1; ich<4; ich++ ) {
          histoName = GetHistoName(baseIndex[itype], icount, ich, -1, -1);
          num = (TH1*)fList->FindObject(histoName.Data());
          histoName = GetHistoName(baseIndex[itype], kNcounts-1, ich, -1, -1);
          den = (TH1*)fList->FindObject(histoName.Data());
          if ( ! num || ! den ) continue;
          effGraph = new TGraphAsymmErrors(num, den, "e0");
          effGraph->GetYaxis()->SetRangeUser(0., 1.1);
          effGraph->GetYaxis()->SetTitle("Efficiency");
          effGraph->GetXaxis()->SetTitle(baseName[itype].Data());
          can->cd(ich+1);
          effGraph->Draw("AP");
          if ( itype == 0 ) break;
        } // loop on chamber
      } // loop on count types
    } // loop on histo
  }

  
  if ( ! outFileOpt.IsNull() ) {
    TObjArray* outFileOptList = outFileOpt.Tokenize("?");
    AliInfo(Form("Creating file %s", outFileOptList->At(0)->GetName()));
    TList* effList = GetEffHistoList(outFileOptList->At(1)->GetName(), outFileOptList->At(2)->GetName(), outFileOptList->At(3)->GetName(), outFileOptList->At(4)->GetName());
    effList->SetName(GetOutputSlot(2)->GetContainer()->GetName());
    TFile* outFile = TFile::Open(outFileOptList->At(0)->GetName(), "RECREATE");
    effList->Write(effList->GetName(),TObject::kSingleKey);
    outFile->Close();
    delete effList;
    delete outFileOptList;
  }
}
