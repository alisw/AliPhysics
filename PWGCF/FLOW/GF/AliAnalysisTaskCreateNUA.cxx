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
#include "AliAnalysisTaskCreateNUA.h"
#include "AliMultSelection.h"
#include <vector>

class AliAnalysisTaskCreateNUA; // your analysis class

using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCreateNUA) // classimp: necessary for root

    AliAnalysisTaskCreateNUA::AliAnalysisTaskCreateNUA() : AliAnalysisTaskSE(),
                                                           fAOD(0), fOutputList(0), fIsMC(false),
  fPeriod("LHC17"),
  Last_RunNumer(0),
  Last_Position(-1),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fMinPt(0.2), fMaxPt(3.0), hEventCount(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
  weights = nullptr;
}
//_____________________________________________________________________________
AliAnalysisTaskCreateNUA::AliAnalysisTaskCreateNUA(const char *name) : AliAnalysisTaskSE(name),
                                                                       fAOD(0), fOutputList(0), fIsMC(false), fPeriod("LHC17"),
  Last_RunNumer(0),
  Last_Position(-1),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fMinPt(0.2),
  fMaxPt(3.0),
  hEventCount(0)
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
AliAnalysisTaskCreateNUA::~AliAnalysisTaskCreateNUA()
{
  // destructor
  if (fOutputList)
  {
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
  }
  if (fGFWSelection)
    delete fGFWSelection;
  if (fGFWSelection15o)
    delete fGFWSelection15o;
}
//_____________________________________________________________________________
void AliAnalysisTaskCreateNUA::UserCreateOutputObjects()
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
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1 and LHC17n
    fGFWSelection15o = new AliGFWNFCuts();
    fGFWSelection15o->PrintSetup();
  } else {
    fGFWSelection = new AliGFWMCuts();
    fGFWSelection->PrintSetup();
  }

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
  else if (fPeriod == "LHC15o")
  {
    vector<int> tempPbPb = {244917, 244918, 244975, 244980, 244982, 244983, 245061, 245064, 245066, 245068, 245145, 245146, 245148, 245151, 245152, 245231, 245232, 245233, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245410, 245411, 245439, 245441, 245446, 245450, 245452, 245453, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554, 245683, 245692, 245700, 245702, 245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793, 245829, 245831, 245833, 245923, 245949, 245952, 245954, 245963, 246001, 246003, 246012, 246036, 246037, 246042, 246048, 246049, 246052, 246053, 246087, 246089, 246113, 246115, 246148, 246151, 246152, 246153, 246178, 246180, 246181, 246182, 246185, 246217, 246222, 246225, 246271, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428, 246431, 246434, 246487, 246488, 246493, 246495, 246540, 246543, 246553, 246567, 246568, 246575, 246583, 246648, 246671, 246675, 246676, 246750, 246751, 246757, 246758, 246759, 246760, 246763, 246765, 246766, 246804, 246805, 246807, 246808, 246809, 246810, 246844, 246845, 246846, 246847, 246851, 246855, 246858, 246859, 246864, 246865, 246867, 246870, 246871, 246928, 246930, 246937, 246942, 246945, 246948, 246949, 246980, 246982, 246984, 246989, 246991, 246994};
    RunNumber_LHC15o = tempPbPb;

    if (!RunNumber_LHC15o.empty())
    {
      for (auto RunIter = RunNumber_LHC15o.begin(); RunIter != RunNumber_LHC15o.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
  } else if (fPeriod == "LHC15o_pass2") {
    vector<int> tempPbPb = {246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246434, 246431, 246424, 246392, 246391, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245793, 245785, 245775, 245766, 245759, 245752, 245731, 245729, 245705, 245702, 245692, 245683, 245554, 245545, 245544, 245543, 245542, 245540, 245535, 245507, 245505, 245504, 245501, 245497, 245496, 245454, 245453, 245450, 245446, 245441, 245411, 245410, 245409, 245407, 245401, 245397, 245396, 245353, 245349, 245347, 245346, 245345, 245343, 245259, 245233, 245232, 245231, 245152, 245151, 245146, 245145, 245068, 245066, 245064, 244983, 244982, 244980, 244975, 244918, 244917};
    RunNumber_LHC15opass2 = tempPbPb;

    if (!RunNumber_LHC15opass2.empty())
    {
      for (auto RunIter = RunNumber_LHC15opass2.begin(); RunIter != RunNumber_LHC15opass2.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
  } else if (fPeriod == "LHC18qr_pass3") {
    vector<int> tempPbPb = {296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551, 296550, 296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420, 296419, 296415, 296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304, 296303, 296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241, 296240, 296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142, 296135, 296134, 296133, 296132, 296123, 296074, 296066, 296065, 296063, 296062, 296060, 296016, 295942, 295941, 295937, 295936, 295913, 295910, 295909, 295861, 295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825, 295822, 295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 295759, 295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718, 295717, 295714, 295712, 295676, 295675, 295673, 295668, 295667, 295666, 295615, 295612, 295611, 295610, 295589, 295588, 295586, 295585,
                            297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380, 297379, 297372, 297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297311, 297310, 297278, 297222, 297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 297117, 297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 296899, 296894, 296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794, 296793, 296790, 296787, 296786, 296785, 296784, 296781, 296752, 296694, 296693, 296691, 296690
    };
    RunNumber_LHC18qrpass3 = tempPbPb;

    if (!RunNumber_LHC18qrpass3.empty())
      {
        for (auto RunIter = RunNumber_LHC18qrpass3.begin(); RunIter != RunNumber_LHC18qrpass3.end(); RunIter++)
          {
            fOutputList->Add(new AliGFWWeights());
            AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
            if (fCurrSystFlag == 0)
              {
                fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
                fWeights->Init(!fIsMC, fIsMC);
              }
            else
              {
                fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
                fWeights->Init(!fIsMC, fIsMC);
              }
          }
      }
  }
  else if (fPeriod == "LHC16qt") {
     vector<int> tempPbPb = {
       265309, 265332, 265334, 265335, 265336, 265338, 265339, 265342, 265343, 265344, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525,
       267163, 267164, 267165, 267166
     };
    RunNumber_LHC16qt = tempPbPb;

    if (!RunNumber_LHC16qt.empty())
      {
        for (auto RunIter = RunNumber_LHC16qt.begin(); RunIter != RunNumber_LHC16qt.end(); RunIter++)
          {
            fOutputList->Add(new AliGFWWeights());
            AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
            if (fCurrSystFlag == 0)
              {
                fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
                fWeights->Init(!fIsMC, fIsMC);
              }
            else
              {
                fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
                fWeights->Init(!fIsMC, fIsMC);
              }
          }
      }
  } else {
    printf("RunNumber list is empty!\n");
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
void AliAnalysisTaskCreateNUA::UserExec(Option_t *)
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
  if(fUseHM){
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  }
  else{
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7, true);
  }
  
  //Standard AliEventCuts for events
  if (!fEventCuts.AcceptEvent(fAOD))
  { // automatic event selection for Run2
    PostData(1, fOutputList);
    return;
  }
  hEventCount->GetXaxis()->SetBinLabel(3,"After fEventCuts");
  hEventCount->Fill(2.5);
  //AliGFWCuts for Sysmatics

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
    fGFWSelection15o->ResetCuts();
    fGFWSelection15o->SetupCuts(fCurrSystFlag);
    if (!fGFWSelection15o->AcceptVertex(fAOD))
    {
      PostData(1, fOutputList);
      return;
    }
  } else {
    fGFWSelection->ResetCuts();
    fGFWSelection->SetupCuts(fCurrSystFlag);
    if (!fGFWSelection->AcceptVertex(fAOD))
      {
        PostData(1, fOutputList);
        return;
    }
  }
  hEventCount->GetXaxis()->SetBinLabel(4,"After AliGFWCuts");
  hEventCount->Fill(3.5);
  Int_t iTracks(fAOD->GetNumberOfTracks()); // see how many tracks there are in the event
  Int_t index = ReturnPosi_WeightList(fAOD->GetRunNumber(), fCurrSystFlag);
  
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
void AliAnalysisTaskCreateNUA::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

//Return the position of w in WeightList
const Int_t AliAnalysisTaskCreateNUA::ReturnPosi_WeightList(const Int_t RunNumber, const Int_t sysflag)
{
  if (fPeriod == "LHC16") {
    TString Period=ReturnPPperiod(RunNumber);
    if (!Period_LHC16.empty()) {
      auto it = find(Period_LHC16.begin(), Period_LHC16.end(), Period);
      if (it == Period_LHC16.end()) {
        return -1;
      } //can't find the RunNumber
      int index = it - Period_LHC16.begin();
      return index;
    }
  } else if (fPeriod == "LHC17") {
    if (!Period_LHC17.empty()) {
      TString Period=ReturnPPperiod(RunNumber);
      auto it = find(Period_LHC17.begin(), Period_LHC17.end(), Period);
      if (it == Period_LHC17.end()) {
        return -1;
      } //can't find the RunNumber
      int index = it - Period_LHC17.begin();
      return index;
    }
  } else if (fPeriod == "LHC18") {
    if (!Period_LHC18.empty()) {
      TString Period=ReturnPPperiod(RunNumber);
      auto it = find(Period_LHC18.begin(), Period_LHC18.end(), Period);
      if (it == Period_LHC18.end()) {
        return -1;
      } //can't find the RunNumber
      int index = it - Period_LHC18.begin();
      return index;
    }
  } else if (fPeriod == "LHC15o") {
    if (!RunNumber_LHC15o.empty())
      {
        auto it = find(RunNumber_LHC15o.begin(), RunNumber_LHC15o.end(), RunNumber);
        if (it == RunNumber_LHC15o.end())
          {
            return -1;
          } //can't find the RunNumber
        int index = it - RunNumber_LHC15o.begin();
        return index;
      }
  } else if (fPeriod == "LHC15o_pass2") {
    if (!RunNumber_LHC15opass2.empty())
      {
        auto it = find(RunNumber_LHC15opass2.begin(), RunNumber_LHC15opass2.end(), RunNumber);
        if (it == RunNumber_LHC15opass2.end())
          {
            return -1;
          } //can't find the RunNumber
        int index = it - RunNumber_LHC15opass2.begin();
        return index;
      }
  } else if (fPeriod == "LHC18qr_pass3") {
    if (!RunNumber_LHC18qrpass3.empty())
      {
        auto it = find(RunNumber_LHC18qrpass3.begin(), RunNumber_LHC18qrpass3.end(), RunNumber);
        if (it == RunNumber_LHC18qrpass3.end())
          {
            return -1;
          } //can't find the RunNumber
        int index = it - RunNumber_LHC18qrpass3.begin();
        return index;
      }
  } else if (fPeriod == "LHC16qt") {
    if (!RunNumber_LHC16qt.empty())
      {
        auto it = find(RunNumber_LHC16qt.begin(), RunNumber_LHC16qt.end(), RunNumber);
        if (it == RunNumber_LHC16qt.end())
          {
            return -1;
          } //can't find the RunNumber
        int index = it - RunNumber_LHC16qt.begin();
        return index;
      }
  } else {
    return -1;
  }

  return -1;
}

Bool_t AliAnalysisTaskCreateNUA::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp)
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
    ltrackXYZ[2] = abs(ltrackXYZ[2] - vtxp[2]);
  }
  else
    return kFALSE; //DCA cut is a must for now

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
    return fGFWSelection15o->AcceptTrack(mtr, ltrackXYZ, 0, kFALSE);
  } else {
    return fGFWSelection->AcceptTrack(mtr, ltrackXYZ, 0, kFALSE);
  }
}

//============================================================================
const char *AliAnalysisTaskCreateNUA::ReturnPPperiod(const Int_t runNumber) const
{
  //根据RunNumber返回不同的Period字符串
  Bool_t isHM = kFALSE;
  if(fUseHM)isHM = kTRUE;

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


/*
  // Run by Run NUA in pp will not be used
  //SetRunNumber
  if (fPeriod == "LHC16")
  {
    vector<int> temp = {252330, 252326, 252325, 252322, 252319, 252317, 252310, 252271, 252248, 252235, 253591, 253589, 253563, 253530, 253529, 253517, 253488, 253482, 253481, 253478, 253437, 254332, 254331, 254330, 254304, 254303, 254302, 254293, 254205, 254204, 254199, 254193, 254178, 254175, 254174, 254149, 254147, 254128, 255467, 255466, 255465, 255463, 255447, 255442, 255440, 255421, 255420, 255419, 255418, 255415, 255407, 255402, 255398, 255352, 255351, 255350, 255283, 255280, 255276, 255275, 255256, 255255, 255253, 255252, 255251, 255249, 255248, 255247, 255242, 255240, 255182, 255181, 255180, 255177, 255176, 255174, 255173, 255171, 255167, 255162, 255159, 255154, 255111, 255091, 255086, 255085, 255082, 255079, 254984, 254983, 254654, 254653, 254652, 254651, 254649, 254648, 254646, 254644, 254640, 254632, 254630, 254629, 254621, 254606, 254604, 255618, 255617, 255616, 255615, 255614, 255591, 255583, 255582, 255577, 255543, 255542, 255541, 255540, 255539, 256418, 256417, 256415, 256373, 256372, 256371, 256368, 256366, 256365, 256364, 256363, 256362, 256361, 256356, 256311, 256309, 256307, 256302, 256299, 256297, 256295, 256292, 256290, 256289, 256287, 256284, 256283, 256282, 256281, 256231, 256228, 256227, 256223, 256219, 258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962, 264035, 264033, 263985, 263984, 263981, 263978, 263977, 263923, 263920, 263917, 263916, 263905, 263866, 263863, 263810, 263803, 263793, 263792, 263790, 263787, 263786, 263785, 263784, 263744, 263743, 263741, 263739, 263738, 263737, 263691, 263690, 263682, 263663, 263662, 263657, 263654, 263652, 263647, 263529, 263497, 263496, 263490, 263487, 263332, 263331, 262858, 262855, 262853, 262849, 262847, 262844, 262842, 262841, 262778, 262777, 262776, 262768, 262760, 262727, 262725, 262723, 262719, 262717, 262713, 262708, 262706, 262705, 262428, 262426, 262425, 262424, 264347, 264346, 264345, 264341, 264336, 264312, 264306, 264305, 264281, 264279, 264277, 264273, 264267, 264266, 264265, 264264, 264262, 264261, 264260, 264259, 264238, 264235, 264233, 264232, 264198, 264197, 264194, 264190, 264188, 264168, 264164, 264139, 264138, 264137, 264129, 264110, 264109, 264086, 264085, 264082, 264078, 264076};
    RunNumber_LHC16 = temp;

    if (!RunNumber_LHC16.empty())
    {
      for (auto RunIter = RunNumber_LHC16.begin(); RunIter != RunNumber_LHC16.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
    else
    {
      printf("RunNumber list is empty!\n");
    }
  }
  else if (fPeriod == "LHC17")
  {
    vector<int> temp = {270830, 270828, 270827, 270824, 270822, 274442, 274390, 274389, 274388, 274387, 274386, 274385, 274364, 274363, 274360, 274352, 274329, 274283, 274281, 274280, 274278, 274276, 274271, 274270, 274269, 274268, 274266, 274264, 274263, 274259, 274258, 274232, 274212, 274174, 274148, 274147, 274125, 274094, 274092, 274058, 273986, 273985, 273946, 273943, 273942, 273918, 273889, 273887, 273886, 273885, 273825, 273824, 273654, 273653, 273593, 273592, 273591, 274671, 274669, 274667, 274657, 274653, 274601, 274596, 274595, 274594, 274593, 276508, 276507, 276506, 276462, 276439, 276438, 276437, 276435, 276351, 276348, 276302, 276297, 276294, 276292, 276290, 276259, 276257, 276230, 276205, 276178, 276177, 276170, 276169, 276166, 276145, 276140, 276135, 276104, 276102, 276099, 276098, 276097, 275847, 275664, 275661, 275650, 275648, 275647, 275624, 275623, 275622, 275621, 275617, 275612, 275559, 275558, 275515, 275472, 275471, 275467, 275459, 275457, 275456, 275453, 275452, 275448, 275443, 275406, 275404, 275401, 275372, 275369, 275361, 275360, 275333, 275332, 275328, 275326, 275324, 275322, 275314, 275283, 275247, 275246, 275245, 275239, 275188, 275184, 275180, 275177, 275174, 275173, 275151, 275150, 275149, 275076, 275075, 275073, 275068, 275067, 274979, 274978, 274886, 274882, 274878, 274877, 274822, 274821, 274815, 274806, 274803, 274802, 274801, 274708, 274690, 278216, 278215, 278191, 278189, 278167, 278166, 278165, 278164, 278158, 278127, 278126, 278123, 278122, 278121, 277996, 277991, 277989, 277987, 277952, 277930, 277907, 277904, 277903, 277900, 277899, 277898, 277897, 277876, 277870, 277848, 277847, 277845, 277842, 277841, 277836, 277834, 277805, 277802, 277801, 277800, 277799, 277795, 277794, 277749, 277747, 277746, 277745, 277725, 277723, 277722, 277721, 277577, 277576, 277575, 277574, 277537, 277536, 277534, 277531, 277530, 277479, 277478, 277477, 277476, 277473, 277472, 277418, 277417, 277416, 277389, 277386, 277385, 277384, 277383, 277360, 277314, 277312, 277310, 277293, 277262, 277257, 277256, 277197, 277196, 277194, 277193, 277189, 277188, 277184, 277183, 277182, 277181, 277180, 277155, 277121, 277117, 277091, 277087, 277082, 277079, 277076, 277073, 277037, 277017, 277016, 277015, 276972, 276971, 276970, 276969, 276967, 276920, 276917, 276916, 276762, 276675, 276674, 276672, 276671, 276670, 276644, 276608, 276557, 276556, 276553, 276552, 276551, 280140, 280135, 280134, 280131, 280126, 280118, 280114, 280111, 280108, 280107, 280066, 280052, 280051, 279879, 279855, 279854, 279853, 279830, 279827, 279826, 279773, 279749, 279747, 279719, 279718, 279715, 279689, 279688, 279687, 279684, 279683, 279682, 279679, 279677, 279676, 279642, 279641, 279632, 279630, 279559, 279550, 279491, 279488, 279487, 279483, 279441, 279439, 279435, 279410, 279391, 279355, 279354, 279349, 279348, 279344, 279342, 279312, 279310, 279309, 279274, 279273, 279270, 279268, 279267, 279265, 279264, 279242, 279238, 279235, 279234, 279232, 279208, 279207, 279201, 279199, 279157, 279155, 279130, 279123, 279122, 279118, 279117, 279107, 279106, 279075, 279074, 279073, 279069, 279068, 279044, 279043, 279041, 279036, 279035, 279008, 279007, 279005, 279000, 278999, 278964, 278963, 278960, 278959, 278941, 278939, 278936, 278915, 278914, 281961, 281956, 281953, 281940, 281939, 281932, 281931, 281928, 281920, 281918, 281916, 281915, 281895, 281894, 281893, 281892, 281633, 281592, 281583, 281574, 281569, 281568, 281563, 281562, 281557, 281511, 281509, 281477, 281475, 281450, 281449, 281446, 281444, 281443, 281441, 281415, 281321, 281301, 281277, 281275, 281273, 281271, 281244, 281243, 281242, 281241, 281240, 281213, 281212, 281191, 281190, 281189, 281181, 281180, 281179, 281081, 281080, 281062, 281061, 281060, 281036, 281035, 281033, 281032, 280999, 280998, 280997, 280996, 280994, 280990, 280947, 280943, 280940, 280936, 280897, 280880, 280856, 280854, 280849, 280848, 280847, 280844, 280842, 280793, 280792, 280787, 280786, 280768, 280767, 280766, 280765, 280764, 280763, 280762, 280761, 280757, 280756, 280755, 280754, 280753, 280729, 280706, 280705, 280681, 280679, 280671, 280647, 280645, 280639, 280637, 280636, 280634, 280613, 280583, 280581, 280574, 280551, 280550, 280547, 280546, 280519, 280518, 280499, 280490, 280448, 280447, 280446, 280445, 280443, 280419, 280415, 280412, 280406, 280405, 280403, 280375, 280374, 280351, 280350, 280349, 280348, 280312, 280310, 280290, 280286, 280285, 280284, 280282, 282704, 282703, 282702, 282700, 282677, 282676, 282673, 282671, 282670, 282667, 282666, 282651, 282629, 282622, 282620, 282618, 282609, 282608, 282607, 282606, 282580, 282579, 282575, 282573, 282546, 282545, 282544, 282528, 273103, 273100, 273099, 273077, 273010, 273009, 272985, 272983, 272976, 272949, 272947, 272939, 272935, 272934, 272933, 272932, 272905, 272903, 272880, 272873, 272871, 272870, 272836, 272834, 272833, 272829, 272828, 272784, 272783, 272782, 272764, 272763, 272760, 272749, 272747, 272712, 272691, 272690, 272620, 272610, 272608, 272607, 272585, 272577, 272575, 272574, 272521, 272468, 272466, 272463, 272462, 272461, 272413, 272411, 272400, 272399, 272395, 272394, 272389, 272388, 272360, 272359, 272340, 272335, 272194, 272156, 272155, 272154, 272153, 272152, 272151, 272123, 272101, 272100, 272076, 272042, 272040, 272039, 272038, 272036, 272020, 272018, 271886, 271880, 271874, 271873, 271871, 271870, 270667, 270665, 270663, 270661, 270581, 270865, 270861, 270856, 270855, 270854};
    RunNumber_LHC17 = temp;

    if (!RunNumber_LHC17.empty())
    {
      for (auto RunIter = RunNumber_LHC17.begin(); RunIter != RunNumber_LHC17.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
    else
    {
      printf("RunNumber list is empty!\n");
    }
  }
  else if (fPeriod == "LHC18")
  {
    vector<int> temp = {285396, 285365, 285364, 285347, 285328, 285327, 285224, 285222, 285203, 285202, 285200, 285165, 285127, 285125, 285108, 285106, 285066, 285065, 285064, 285015, 285014, 285013, 285012, 285011, 285009, 286350, 286349, 286348, 286345, 286341, 286340, 286337, 286336, 286314, 286313, 286312, 286311, 286310, 286309, 286308, 286289, 286288, 286287, 286284, 286282, 286263, 286261, 286258, 286257, 286254, 286231, 286230, 286229, 286203, 286202, 286201, 286199, 286198, 286159, 286130, 286129, 286127, 286124, 286064, 286025, 286014, 285980, 285979, 285978, 286937, 286936, 286933, 286932, 286931, 286930, 286911, 286910, 286907, 286877, 286876, 286874, 286852, 286850, 286846, 286809, 286805, 286801, 286799, 286731, 286695, 286661, 286653, 286633, 286592, 286591, 286569, 286568, 286567, 286566, 286511, 286509, 286508, 286502, 286482, 286455, 286454, 286428, 286427, 286426, 286380, 287658, 287657, 287656, 287654, 287578, 287575, 287524, 287521, 287518, 287517, 287516, 287513, 287486, 287484, 287481, 287480, 287451, 287413, 287389, 287388, 287387, 287385, 287381, 287380, 287360, 287356, 287355, 287353, 287349, 287347, 287346, 287344, 287343, 287325, 287324, 287323, 287283, 287254, 287251, 287250, 287249, 287248, 287209, 287208, 287204, 287203, 287202, 287201, 287185, 287155, 287137, 287077, 287072, 287071, 287066, 287064, 287063, 287021, 287000, 288619, 288640, 288642, 288644, 288650, 288687, 288689, 288690, 288743, 288748, 288750, 288804, 288806, 288861, 288862, 288863, 288864, 288868, 288902, 288903, 288908, 288909, 288943, 289971, 289966, 289965, 289943, 289941, 289940, 289935, 289931, 289928, 289884, 289880, 289879, 289857, 289856, 289855, 289854, 289852, 289849, 289830, 289818, 289817, 289816, 289815, 289814, 289811, 289808, 289775, 289757, 289732, 289731, 289729, 289724, 289723, 289721, 289547, 289521, 289494, 289493, 289468, 289466, 289465, 289463, 289462, 289444, 289426, 289374, 289373, 289370, 289369, 289368, 289367, 289366, 289365, 289356, 289355, 289354, 289353, 289309, 289308, 289306, 289303, 289300, 289281, 289280, 289278, 289277, 289276, 289275, 289254, 289253, 289249, 289247, 289243, 289242, 289241, 289240, 292839, 292836, 292834, 292832, 292831, 292811, 292810, 292809, 292804, 292803, 292752, 292750, 292748, 292747, 292744, 292739, 292737, 292704, 292701, 292698, 292696, 292695, 292693, 292586, 292584, 292563, 292560, 292559, 292557, 292554, 292553, 292526, 292524, 292523, 292521, 292500, 292497, 292496, 292495, 292461, 292460, 292457, 292456, 292434, 292432, 292430, 292429, 292428, 292406, 292405, 292398, 292397, 292298, 292273, 292265, 292242, 292241, 292240, 292218, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114, 292109, 292108, 292107, 292106, 292081, 292080, 292077, 292075, 292067, 292062, 292061, 292060, 292040, 292012, 291982, 291977, 291976, 291953, 291948, 291946, 291945, 291944, 291943, 291942, 291803, 291796, 291795, 291769, 291768, 291766, 291762, 291760, 291756, 291755, 291729, 291706, 291698, 291697, 291690, 291665, 291661, 291657, 291626, 291624, 291622, 291618, 291615, 291614, 291590, 291485, 291484, 291482, 291481, 291457, 291456, 291453, 291451, 291447, 291424, 291420, 291417, 291416, 291402, 291400, 291399, 291397, 291377, 291375, 291363, 291362, 291361, 291360, 291286, 291285, 291284, 291282, 291266, 291265, 291263, 291262, 291257, 291240, 291209, 291188, 291143, 291116, 291111, 291110, 291101, 291100, 291093, 291069, 291066, 291065, 291041, 291037, 291035, 291006, 291005, 291004, 291003, 291002, 290980, 290979, 290976, 290975, 290974, 290948, 290944, 290943, 290941, 290935, 290932, 290895, 290894, 290888, 290887, 290886, 290862, 290860, 290853, 290848, 290846, 290843, 290841, 290790, 290787, 290766, 290689, 290687, 290665, 290660, 290645, 290632, 290627, 290615, 290614, 290613, 290612, 290590, 290588, 290553, 290550, 290549, 290544, 290540, 290539, 290538, 290501, 290500, 290499, 290469, 290467, 290459, 290458, 290456, 290427, 290426, 290425, 290423, 290412, 290411, 290404, 290401, 290399, 290376, 290375, 290374, 290350, 290327, 290323, 293357, 293359, 289165, 289166, 289167, 289169, 289172, 289175, 289176, 289177, 289198, 289199, 289200, 289201, 293898, 293896, 293893, 293891, 293886, 293856, 293831, 293830, 293829, 293809, 293807, 293806, 293805, 293802, 293776, 293774, 293773, 293770, 293741, 293740, 293698, 293696, 293695, 293692, 293691, 293588, 293587, 293583, 293582, 293579, 293578, 293573, 293571, 293570, 293475, 294925, 294916, 294884, 294883, 294880, 294875, 294852, 294818, 294817, 294816, 294815, 294813, 294809, 294805, 294775, 294774, 294772, 294769, 294749, 294747, 294746, 294745, 294744, 294742, 294741, 294722, 294718, 294715, 294710, 294703, 294653, 294636, 294633, 294632, 294593, 294591, 294590, 294587, 294586, 294563, 294562, 294558, 294556, 294553, 294531, 294530, 294529, 294527, 294526, 294525, 294524, 294310, 294308, 294307, 294242, 294241, 294212, 294210, 294208, 294205, 294201, 294200, 294199, 294156, 294155, 294154, 294152, 294131, 294013, 294012, 294011, 294010, 294009};
    RunNumber_LHC18 = temp;

    if (!RunNumber_LHC18.empty())
    {
      for (auto RunIter = RunNumber_LHC18.begin(); RunIter != RunNumber_LHC18.end(); RunIter++)
      {
        fOutputList->Add(new AliGFWWeights());
        AliGFWWeights *fWeights = (AliGFWWeights *)fOutputList->Last();
        if (fCurrSystFlag == 0)
        {
          fWeights->SetName(Form("w%s", to_string(*RunIter).c_str()));
          fWeights->Init(!fIsMC, fIsMC);
        }
        else
        {
          fWeights->SetName(Form("w%s_SystFlag%d_", to_string(*RunIter).c_str(), fCurrSystFlag));
          fWeights->Init(!fIsMC, fIsMC);
        }
      }
    }
    else
    {
      printf("RunNumber list is empty!\n");
    }
  }
*/
