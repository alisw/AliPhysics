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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCreatePeriodNUA.h"
#include "AliMultSelection.h"
#include <vector>

class AliAnalysisTaskCreatePeriodNUA; // your analysis class

using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCreatePeriodNUA) // classimp: necessary for root

    AliAnalysisTaskCreatePeriodNUA::AliAnalysisTaskCreatePeriodNUA() : AliAnalysisTaskSE(),
                                                           fAOD(0), fOutputList(0), fIsMC(false), fPeriod("LHC17"), Last_RunNumer(0),
                                                           Last_Position(-1), fGFWSelection(NULL), fMinPt(0.2), fMaxPt(3.0), hEventCount(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
  weights = nullptr;
}
//_____________________________________________________________________________
AliAnalysisTaskCreatePeriodNUA::AliAnalysisTaskCreatePeriodNUA(const char *name) : AliAnalysisTaskSE(name),
                                                                       fAOD(0), fOutputList(0), fIsMC(false), fPeriod("LHC17"), Last_RunNumer(0),
                                                                       Last_Position(-1), fGFWSelection(NULL), fMinPt(0.2), fMaxPt(3.0), hEventCount(0)
{
  weights = nullptr;
  // constructor
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                   // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
                                   // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskCreatePeriodNUA::~AliAnalysisTaskCreatePeriodNUA()
{
  // destructor
  if (fOutputList)
  {
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
  }
  if (fGFWSelection)
    delete fGFWSelection;
}
//_____________________________________________________________________________
void AliAnalysisTaskCreatePeriodNUA::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //
  fOutputList = new TList(); // this is a list which will contain all of your histograms
                             // at the end of the analysis, the contents of this list are written
                             // to the output file
  fOutputList->SetName("WeightList");
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all objects it contains and will delete them
                                // if requested (dont worry about this now)
  fGFWSelection = new AliGFWCuts();
  fGFWSelection->PrintSetup();

  //SetRunNumber
  if (fPeriod == "LHC16")
  {
    std::vector<TString> temp = {"averaged","LHC16de","LHC16ghi","LHC16j","LHC16k","LHC16l","LHC16o","LHC16p"};
    Period_LHC16 = temp;

    if (!Period_LHC16.empty())
    {
      for (auto RunIter = Period_LHC16.begin(); RunIter != Period_LHC16.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", (*RunIter).Data()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", (*RunIter).Data(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }

      }
    }
    else
    {
      printf("Period list is empty!\n");
    }
  }
  else if (fPeriod == "LHC17")
  {
    std::vector<TString> temp = {"averaged","LHC17ce","LHC17f","LHC17h","LHC17i","LHC17j","LHC17k","LHC17l","LHC17m","LHC17o","LHC17r"};
    Period_LHC17=temp;

    if (!Period_LHC17.empty())
    {
      for (auto RunIter = Period_LHC17.begin(); RunIter != Period_LHC17.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", (*RunIter).Data()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", (*RunIter).Data(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
    else
    {
      printf("Period_LHC17 list is empty!\n");
    }
  }
  else if (fPeriod == "LHC18")
  {
    std::vector<TString> temp = {"averaged","LHC18b","LHC18d","LHC18bd","LHC18e","LHC18f","LHC18hjk","LHC18ghijk","LHC18l","LHC18m","LHC18mn","LHC18o","LHC18p"};
    Period_LHC18=temp;

    if (!Period_LHC18.empty())
    {
      for (auto RunIter = Period_LHC18.begin(); RunIter != Period_LHC18.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", (*RunIter).Data()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", (*RunIter).Data(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
    else
    {
      printf("Period_LHC18 list is empty!\n");
    }
  }

  //must add in the last
  //avoid changing the index of weight
  hEventCount = new TH1D("hEventCount", "; centrality;;", 5, 0, 5);
	fOutputList->Add(hEventCount);

  //fOutputList->Add(weights);
  // don't forget to add it to the list! the list will be written to file, so if you want
  // your histogram in the output file, add it to the list!

  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the
                            // fOutputList object. the manager will in the end take care of writing your output to file
                            // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskCreatePeriodNUA::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you
  // have access to the current event.
  // once you return from the UserExec function, the manager will retrieve the next event from the chain
  hEventCount->GetXaxis()->SetBinLabel(1,"Loop Number");
  hEventCount->Fill(0.5);
  
  AliMultSelection *lMultSel = (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  Double_t cent = lMultSel->GetMultiplicityPercentile("V0M");


  fAOD = dynamic_cast<AliAODEvent *>(InputEvent()); // get an event (called fAOD) from the input file
                                                    // there's another event format (ESD) which works in a similar wya
                                                    // but is more cpu/memory unfriendly. for now, we'll stick with aod's
  if (!fAOD)
    return; // if the pointer to the event is empty (getting it failed) skip this event
            // example part: i'll show how to loop over the tracks in an event
            // and extract some information from them which we'll store in a histogram
  hEventCount->GetXaxis()->SetBinLabel(2,"AOD OK");
  hEventCount->Fill(1.5);
  //Standard AliEvent Cuts
  if(fUseHM)
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  //Standard AliEventCuts for events
  if (!fEventCuts.AcceptEvent(fAOD))
  { // automatic event selection for Run2
    PostData(1, fOutputList);
    return;
  }
  hEventCount->GetXaxis()->SetBinLabel(3,"After fEventCuts");
  hEventCount->Fill(2.5);
  //AliGFWCuts for Sysmatics
  fGFWSelection->ResetCuts();
  fGFWSelection->SetupCuts(fCurrSystFlag);
  if (!fGFWSelection->AcceptVertex(fAOD))
  {
    PostData(1, fOutputList);
    return;
  }
  hEventCount->GetXaxis()->SetBinLabel(4,"After AliGFWCuts");
  hEventCount->Fill(3.5);
  Int_t iTracks(fAOD->GetNumberOfTracks()); // see how many tracks there are in the event
  Int_t index = ReturnPosi_WeightList(fAOD->GetRunNumber());
  
  //..for DCA
  Double_t pos[3], vz, vx, vy;
  vz = fAOD->GetPrimaryVertex()->GetZ();
  vx = fAOD->GetPrimaryVertex()->GetX();
  vy = fAOD->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};

  if (index > -1)
  {
    for (Int_t i(0); i < iTracks; i++)
    {                                                                     // loop ove rall these tracks
      AliAODTrack *track = static_cast<AliAODTrack *>(fAOD->GetTrack(i)); // get a track (type AliAODTrack) from the event
      if (!track)
      {
        delete track;
        continue;
      } // if we failed, skip this track
      //((AliGFWWeights *)fOutputList->At(674))->Fill(track->Phi(), track->Eta(), vz, track->Pt(), cent, 0);
      track->GetXYZ(pos);
      if (!AcceptAODTrack(track, pos, vtxp))
        continue;
      ((AliGFWWeights *)fOutputList->At(index))->Fill(track->Phi(), track->Eta(), vz, track->Pt(), cent, 0);
      //weights->Fill(track->Phi(),track->Eta(),vz,track->Pt(),cent,0);
    } // continue until all the tracks are processed
  }
  else
  {
    printf("!!!!!! Can't Find RunNumber %d in %s !!!!!!!!\n", fAOD->GetRunNumber(), fPeriod.Data());
  }
  hEventCount->GetXaxis()->SetBinLabel(5,"Final pass");
  hEventCount->Fill(4.5);
  PostData(1, fOutputList); // stream the results the analysis of this event to
                            // the output manager which will take care of writing
                            // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskCreatePeriodNUA::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

//Return the position of w in WeightList
const Int_t AliAnalysisTaskCreatePeriodNUA::ReturnPosi_WeightList(const Int_t RunNumber)
{
  TString Period=ReturnPPperiod(RunNumber);
  if (fPeriod == "LHC16")
  {
    if (!Period_LHC16.empty())
    {
      auto it = find(Period_LHC16.begin(), Period_LHC16.end(), Period);
      if (it == Period_LHC16.end())
      {
        return -1;
      } //can't find the RunNumber
      int index = &*it - &Period_LHC16[0];
      return index;

    }
  }
  else if (fPeriod == "LHC17")
  {
    if (!Period_LHC17.empty())
    {
      auto it = find(Period_LHC17.begin(), Period_LHC17.end(), Period);
      if (it == Period_LHC17.end())
      {
        return -1;
      } //can't find the RunNumber
      int index = &*it - &Period_LHC17[0];
      return index;
    }
  }
  else if (fPeriod == "LHC18")
  {
    if (!Period_LHC18.empty())
    {
      auto it = find(Period_LHC18.begin(), Period_LHC18.end(), Period);
      if (it == Period_LHC18.end())
      {
        return -1;
      } //can't find the RunNumber
      int index = &*it - &Period_LHC18[0];
      return index;
    }
  }
  else
  {
    return -1;
  }
  
}

Bool_t AliAnalysisTaskCreatePeriodNUA::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp)
{
  // Pt cut
  if (mtr->Pt() < fMinPt)
    return kFALSE;
  if (mtr->Pt() > fMaxPt)
    return kFALSE;

  // DCA cut
  if (ltrackXYZ && vtxp)
  {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0] - vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1] - vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2] - vtxp[2];
  }
  else
    return kFALSE; //DCA cut is a must for now

  return fGFWSelection->AcceptTrack(mtr, ltrackXYZ, 0, kFALSE);
}

//============================================================================
const char *AliAnalysisTaskCreatePeriodNUA::ReturnPPperiod(const Int_t runNumber) const
{
  //根据RunNumber返回不同的Period字符串
  Bool_t isHM = kFALSE;
  if (fUseHM)
    isHM = kTRUE;

  if (runNumber >= 252235 && runNumber <= 264347)
  { // LHC16
    if (!isHM && runNumber >= 252235 && runNumber <= 252375)
      return "LHC16de"; //d
    if (!isHM && runNumber >= 253437 && runNumber <= 253591)
      return "LHC16de"; //e
    if (runNumber >= 254128 && runNumber <= 254332)
      return "LHC16ghi"; //g
    if (runNumber >= 254604 && runNumber <= 255467)
      return "LHC16ghi"; //h
    if (runNumber >= 255539 && runNumber <= 255618)
      return "LHC16ghi"; //i
    if (runNumber >= 256219 && runNumber <= 256418)
      return "LHC16j";
    if (runNumber >= 256941 && runNumber <= 258537)
      return "LHC16k";
    if (runNumber >= 258962 && runNumber <= 259888)
      return "LHC16l";
    if (runNumber >= 262424 && runNumber <= 264035)
      return "LHC16o";
    if (runNumber >= 264076 && runNumber <= 264347)
      return "LHC16p";
  }

  if (runNumber >= 270581 && runNumber <= 282704)
  { // LHC17
    if (!isHM && runNumber >= 270581 && runNumber <= 270667)
      return "LHC17ce";
    if (runNumber >= 270822 && runNumber <= 270830)
    {
      if (isHM)
        return "averaged";
      else
        return "LHC17ce";
    }
    if (runNumber >= 270854 && runNumber <= 270865)
    {
      if (isHM)
        return "averaged";
      else
        return "LHC17f";
    }
    if (runNumber >= 271870 && runNumber <= 273103)
      return "LHC17h";
    if (runNumber >= 273591 && runNumber <= 274442)
      return "LHC17i";
    if (!isHM && runNumber >= 274593 && runNumber <= 274671)
      return "LHC17j";
    if (runNumber >= 274690 && runNumber <= 276508)
      return "LHC17k";
    if (runNumber >= 276551 && runNumber <= 278216)
      return "LHC17l";
    if (runNumber >= 278914 && runNumber <= 280140)
      return "LHC17m";
    if (runNumber >= 280282 && runNumber <= 281961)
      return "LHC17o";
    if (runNumber >= 282528 && runNumber <= 282704)
      return "LHC17r";
  }

  if (runNumber >= 285009 && runNumber <= 294925)
  { // LHC18
    if (runNumber >= 285009 && runNumber <= 285396)
    {
      if (isHM)
        return "LHC18bd";
      else
        return "LHC18b";
    }
    if (runNumber >= 285978 && runNumber <= 286350)
    {
      if (isHM)
        return "LHC18bd";
      else
        return "LHC18d";
    }
    if (runNumber >= 286380 && runNumber <= 286937)
      return "LHC18e";
    if (runNumber >= 287000 && runNumber <= 287658)
      return "LHC18f";
    if (runNumber >= 288804 && runNumber <= 288806)
    {
      if (isHM)
        return "LHC18hjk";
      else
        return "LHC18ghijk";
    }
    if (runNumber == 288943)
    {
      if (isHM)
        return "LHC18hjk";
      else
        return "LHC18ghijk";
    }
    if (runNumber >= 289165 && runNumber <= 289201)
    {
      if (isHM)
        return "LHC18hjk";
      else
        return "LHC18ghijk";
    }
    if (!isHM && runNumber >= 288619 && runNumber <= 288750)
      return "LHC18ghijk"; //g, no HM event, only MB
    if (!isHM && runNumber >= 288861 && runNumber <= 288909)
      return "LHC18ghijk"; //i, no HM event, only MB
    if (runNumber >= 289240 && runNumber <= 289971)
      return "LHC18l";
    if (runNumber >= 290323 && runNumber <= 292839)
    {
      if (isHM)
        return "LHC18m";
      else
        return "LHC18mn";
    }
    if (!isHM && runNumber >= 293357 && runNumber <= 293359)
      return "LHC18mn"; //n, no HM event, only MB
    if (runNumber >= 293475 && runNumber <= 293898)
      return "LHC18o";
    if (runNumber >= 294009 && runNumber <= 294925)
      return "LHC18p";
  }

  AliWarning("Unknown period! Returning averaged weights");
  return "averaged";
}