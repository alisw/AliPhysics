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
#include "AliTrigChEffOutput.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskTrigChEff) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskTrigChEff::AliAnalysisTaskTrigChEff() :
  AliVAnalysisMuon(),
  fAnalysisOutput(0x0),
  fList(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskTrigChEff::AliAnalysisTaskTrigChEff(const char *name, const AliMuonTrackCuts& cuts) :
  AliVAnalysisMuon(name, cuts),
  fAnalysisOutput(0x0),
  fList(0x0)
{
  //
  /// Constructor.
  //

  DefineOutput(2, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskTrigChEff::~AliAnalysisTaskTrigChEff()
{
  //
  /// Destructor
  //
  delete fAnalysisOutput;
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fList;
  }
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
        for ( Int_t itrackSel=0; itrackSel<AliTrigChEffOutput::kNtrackSel; ++itrackSel ) {
          for ( Int_t imethod=0; imethod<AliTrigChEffOutput::kNeffMethods; ++imethod ) {
            for ( Int_t itype=0; itype<AliTrigChEffOutput::kNhistoTypes; ++itype ) {
              for ( Int_t icount=-1; icount<AliTrigChEffOutput::kNcounts; ++icount ) {
                for ( Int_t ich=-1; ich<4; ++ich ) {
                  for ( Int_t imatch=AliTrigChEffOutput::kMatchApt; imatch<AliTrigChEffOutput::kMatchHpt; ++imatch ) {
                    TH1* histo = 0x0;
                    for ( Int_t jmatch=imatch+1; jmatch<=AliTrigChEffOutput::kMatchHpt; ++jmatch ) {
                      histoName = fAnalysisOutput->GetHistoName(itype, icount, ich, itrackSel, jmatch, imethod);
                      TH1* histoAdd = static_cast<TH1*>(fMergeableCollection->GetObject(Form("/%s/%s/%s/",fPhysSelKeys->At(isel)->GetName(), GetAllSelectedTrigClasses()->At(itrig)->GetName(), GetCentralityClasses()->GetBinLabel(icent)), histoName));
                      if ( ! histoAdd ) continue;
                      histoName = fAnalysisOutput->GetHistoName(itype, icount, ich, itrackSel, imatch, imethod);
                      if ( ! histo ) histo = static_cast<TH1*>(GetMergeableObject(fPhysSelKeys->At(isel)->GetName(), GetAllSelectedTrigClasses()->At(itrig)->GetName(), GetCentralityClasses()->GetBinLabel(icent), histoName));
                      AliDebug(2,Form("Adding %s (%g) to %s (%g)", histoAdd->GetName(), histoAdd->Integral(), histo->GetName(), histo->Integral()));
                      histo->Add(histoAdd);
                    } // loop on higher pt matching
                  } // loop on match trigger
                } // loop on chamber
              } // loop on count type
            } // loop on histo type
          } // loop on eff method
        } // loop on track selection
      } // loop on centrality
    } // loop on trigger classes
  } // loop on physics selection

  TString physSel = fPhysSelKeys->At(kPhysSelPass)->GetName();
  TString trigClass = "ANY";
  TString centrality = "-100_200";
  
  TList* outList = fAnalysisOutput->GetEffHistoList(physSel, trigClass, centrality, AliTrigChEffOutput::kSelectTrack,AliTrigChEffOutput::kMatchApt,AliTrigChEffOutput::kEffFromTrack);
  outList->SetOwner(kFALSE);
  TIter next(outList);
  TObject* obj;
  while ( (obj = next()) ) fList->Add(obj);
  delete outList;

  AliVAnalysisMuon::FinishTaskOutput();
}


//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::MyUserCreateOutputObjects()
{
  //
  /// Create prototype objects for mergeable collection
  //

  TString histoName = "";

  fAnalysisOutput = new AliTrigChEffOutput(fOutputList);

  for ( Int_t itrackSel = 0; itrackSel<AliTrigChEffOutput::kNtrackSel; ++itrackSel ) {
    for ( Int_t imethod=0; imethod<AliTrigChEffOutput::kNeffMethods; ++imethod ) {
      for ( Int_t imatch = 0; imatch<AliTrigChEffOutput::kNtrigMatch; ++imatch ) {
        for ( Int_t icount=0; icount<AliTrigChEffOutput::kNcounts; ++icount ) {
          AddObjectToCollection(fAnalysisOutput->GetCountHisto(AliTrigChEffOutput::kHchamberEff, icount, -1, itrackSel, imatch, imethod));
          for ( Int_t ich=0; ich<4; ++ich ) {
            AddObjectToCollection(fAnalysisOutput->GetCountHisto(AliTrigChEffOutput::kHslatEff, icount, ich, itrackSel, imatch, imethod));
            AddObjectToCollection(fAnalysisOutput->GetCountHisto(AliTrigChEffOutput::kHboardEff, icount, ich, itrackSel, imatch, imethod));
          }
        } // loop on counts

        AddObjectToCollection(fAnalysisOutput->GetCountHisto(AliTrigChEffOutput::kHcheckBoard, -1, -1, itrackSel, imatch, imethod));
      } // loop on trig match
    } // loop on eff method
  } // loop on track selection

  fMuonTrackCuts->Print();

  fList = new TList();
  fList->SetOwner();
  PostData(2, fList);
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
  
  UInt_t addMask[4] = {0, AliMuonTrackCuts::kMuMatchApt, AliMuonTrackCuts::kMuMatchApt|AliMuonTrackCuts::kMuMatchLpt, AliMuonTrackCuts::kMuMatchApt|AliMuonTrackCuts::kMuMatchLpt|AliMuonTrackCuts::kMuMatchHpt};

  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(InputEvent());
  for ( Int_t itrack = 0; itrack < nTracks; ++itrack ) {
    track = AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    Bool_t matchTracker = AliAnalysisMuonUtility::IsMuonTrack(track);

    Int_t matchTrig = AliAnalysisMuonUtility::GetMatchTrigger(track);
    UInt_t selection = fMuonTrackCuts->GetSelectionMask(track);
    
    // Apply the sharp pt cut according to the matched pt level of the track
    UInt_t filterMask = fMuonTrackCuts->GetFilterMask() | addMask[matchTrig];
    Bool_t isSelected = ( ( selection & filterMask ) == filterMask );

    for ( Int_t imethod=0; imethod<AliTrigChEffOutput::kNeffMethods; ++imethod ) {
      if ( imethod == AliTrigChEffOutput::kEffFromTrack ) {
        if ( ! matchTracker || track->P() < 10. ) continue;
        pattern = AliAnalysisMuonUtility::GetMUONTrigHitsMapTrk(track);
        board = AliESDMuonTrack::GetCrossedBoard(pattern);
      }
      else {
        if ( matchTrig == 0 ) continue;
        pattern = AliAnalysisMuonUtility::GetMUONTrigHitsMapTrg(track);
        board = ( AliAnalysisMuonUtility::IsAODEvent(InputEvent()) ) ? AliESDMuonTrack::GetCrossedBoard(pattern) : ((AliESDMuonTrack*)track)->LoCircuit();
      }
      
      Int_t effFlag = AliESDMuonTrack::GetEffFlag(pattern);
      
      if ( effFlag < AliESDMuonTrack::kChEff ) {
        for ( Int_t itrackSel=0; itrackSel<AliTrigChEffOutput::kNtrackSel; ++itrackSel ) {
          if ( itrackSel == AliTrigChEffOutput::kSelectTrack && ! isSelected ) continue;
          for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
            TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
          
            histoName = fAnalysisOutput->GetHistoName(AliTrigChEffOutput::kHcheckBoard, -1, -1, itrackSel, matchTrig, imethod);
            ((TH2F*)GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(AliESDMuonTrack::GetSlatOrInfo(pattern), board);
          }
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
        
        Bool_t fillHisto[AliTrigChEffOutput::kNcounts] = {
          hitsBend,
          hitsNonBend,
          ( hitsBend && hitsNonBend ),
          kTRUE
        };
        
        for ( Int_t itrackSel=0; itrackSel<AliTrigChEffOutput::kNtrackSel; ++itrackSel ) {
          if ( itrackSel == AliTrigChEffOutput::kSelectTrack && ! isSelected ) continue;
          for (Int_t icount=0; icount<AliTrigChEffOutput::kNcounts; ++icount){
            if ( ! fillHisto[icount] ) continue;
            for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
              TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
              
              histoName = fAnalysisOutput->GetHistoName(AliTrigChEffOutput::kHchamberEff, icount, -1, itrackSel, matchTrig, imethod);
              static_cast<TH1*>(GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(iChamber);
            
              if ( effFlag < AliESDMuonTrack::kSlatEff ) continue; // Track crossed different slats
              histoName = fAnalysisOutput->GetHistoName(AliTrigChEffOutput::kHslatEff, icount, ich, itrackSel, matchTrig, imethod);
              static_cast<TH1*>(GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(slat);
            
              if ( effFlag < AliESDMuonTrack::kBoardEff ) continue; // Track crossed different boards
              histoName = fAnalysisOutput->GetHistoName(AliTrigChEffOutput::kHboardEff, icount, ich, itrackSel, matchTrig, imethod);
              static_cast<TH1*>(GetMergeableObject(physSel, trigClassName, centrality, histoName))->Fill(board);
            } // loop on trigger classes
          } // loop on count types
        } // loop on track selection
      } // loop on chambers
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

  delete fAnalysisOutput;
  fAnalysisOutput = new AliTrigChEffOutput(fOutputList);

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
  
  if ( trackSel.GetEntries() == 0 ) trackSel.Add(new TObjString(fAnalysisOutput->GetHistoName(-1,-1,-1,AliTrigChEffOutput::kSelectTrack,AliTrigChEffOutput::kMatchApt,-1)));
  if ( methodSel.GetEntries() == 0 ) methodSel.Add(new TObjString(fAnalysisOutput->GetHistoName(-1, -1, -1, -1, -1, AliTrigChEffOutput::kEffFromTrack)));

  furtherOpt.ToUpper();
  
  Int_t chosenType = ( furtherOpt.Contains("BOARD") ) ? AliTrigChEffOutput::kHboardEff : AliTrigChEffOutput::kHslatEff;

  igroup1++;
  igroup2 = 0;
  TString histoName = "", yAxisTitle = "";
  
  TH1 *num = 0x0;
  TH1 *den = 0x0;
  TGraphAsymmErrors* effGraph = 0x0;
  
  ////////////////
  // Show tests //
  ////////////////
  
  for ( Int_t icount=0; icount<AliTrigChEffOutput::kNcounts-1; ++icount ) {
    currName = fAnalysisOutput->GetHistoName(chosenType, icount, -1, -1, -1, -1);
    currName += "_can";
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
                histoName = fAnalysisOutput->GetHistoName(chosenType, AliTrigChEffOutput::kAllTracks, ich, -1, -1, -1); // partial name
                histoName += Form("%s%s",trackSel.At(itrackSel)->GetName(),methodSel.At(imethodSel)->GetName());
                den = static_cast<TH1*>(fAnalysisOutput->GetSum(physSel->At(isel)->GetName(), trigClasses->At(itrig)->GetName(), centrality->At(icent)->GetName(), histoName.Data()));
                if ( ! den ) {
                  printf("Warning: cannot find %s\n", histoName.Data());
                  continue;
                }
                histoName = fAnalysisOutput->GetHistoName(chosenType, icount, ich, -1, -1, -1); // partial name
                histoName += Form("%s%s",trackSel.At(itrackSel)->GetName(),methodSel.At(imethodSel)->GetName());
                num = static_cast<TH1*>(fAnalysisOutput->GetSum(physSel->At(isel)->GetName(), trigClasses->At(itrig)->GetName(), centrality->At(icent)->GetName(), histoName.Data()));
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
    Int_t baseIndex[3] = {AliTrigChEffOutput::kHchamberEff, AliTrigChEffOutput::kHslatEff, AliTrigChEffOutput::kHboardEff};
    TString effName[AliTrigChEffOutput::kNcounts-1] = {"BendPlane", "NonBendPlane", "BothPlanes"};
    for ( Int_t itype=0; itype<3; itype++ ) {
      for ( Int_t icount=0; icount<AliTrigChEffOutput::kNcounts-1; icount++ ){
        TString canName = Form("efficiencyPer%s_%s",baseName[itype].Data(),effName[icount].Data());
        TCanvas* can = new TCanvas(canName.Data(),canName.Data(),10*(1+AliTrigChEffOutput::kNcounts*itype+icount),10*(1+AliTrigChEffOutput::kNcounts*itype+icount),310,310);
        can->SetFillColor(10); can->SetHighLightColor(10);
        can->SetLeftMargin(0.15); can->SetBottomMargin(0.15);  
        if ( itype > 0 )
          can->Divide(2,2);
        
        for ( Int_t ich=-1; ich<4; ich++ ) {
          histoName = fAnalysisOutput->GetHistoName(baseIndex[itype], icount, ich, -1, -1, -1);
          num = static_cast<TH1*>(fList->FindObject(histoName.Data()));
          histoName = fAnalysisOutput->GetHistoName(baseIndex[itype], AliTrigChEffOutput::kNcounts-1, ich, -1, -1, -1);
          den = static_cast<TH1*>(fList->FindObject(histoName.Data()));
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
    TString histoPattern = outFileOptList->At(4)->GetName();
    Int_t itrackSel = histoPattern.Contains("NoSel") ? AliTrigChEffOutput::kNoSelectTrack : AliTrigChEffOutput::kSelectTrack;
    Int_t imatch = AliTrigChEffOutput::kNoMatch;
    if ( histoPattern.Contains("MatchApt") ) imatch = AliTrigChEffOutput::kMatchApt;
    else if ( histoPattern.Contains("MatchLpt") ) imatch = AliTrigChEffOutput::kMatchLpt;
    else if ( histoPattern.Contains("MatchHpt") ) imatch = AliTrigChEffOutput::kMatchHpt;
    Int_t imethod = histoPattern.Contains("FromTrg") ? AliTrigChEffOutput::kEffFromTrig : AliTrigChEffOutput::kEffFromTrack;
    TList* effList = fAnalysisOutput->GetEffHistoList(outFileOptList->At(1)->GetName(), outFileOptList->At(2)->GetName(), outFileOptList->At(3)->GetName(), itrackSel, imatch, imethod);
    effList->SetName(GetOutputSlot(2)->GetContainer()->GetName());
    if ( effList ->GetEntries() == 0 ) {
      printf("\nWarning: no histograms satisfying the requested conditions.\n(%s %s %s %s).\nOutput %s not created.\n", outFileOptList->At(1)->GetName(), outFileOptList->At(2)->GetName(), outFileOptList->At(3)->GetName(), outFileOptList->At(4)->GetName(),outFileOptList->At(0)->GetName());
    }
    else {
      TFile* outFile = TFile::Open(outFileOptList->At(0)->GetName(), "RECREATE");
      effList->Write(effList->GetName(),TObject::kSingleKey);
      outFile->Close();
    }
    delete effList;
    delete outFileOptList;
  }
}
