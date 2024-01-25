/* AliAnaysisTaskPtN
 * Maintainer: Yifan Zhang
 * calculating the mean Pt and fluctuations with respect to Nch
 */

#include "TChain.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisPtN.h"
#include <cmath>
#include <TComplex.h>
# include <TRandom3.h>
#include <iostream>
//#include "CorrelationCalculator.h"

class AliAnalysisPtN;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisPtN) // classimp: necessary for root

AliAnalysisPtN::AliAnalysisPtN() : AliAnalysisTaskSE(), 
    fAOD(nullptr), 
    fOutputList(nullptr), 
    fBstList(nullptr),
    fWeightsListNUE(nullptr),
    fWeightNUE(nullptr),
    //correlator(),
    fTestNonWeight(nullptr),
    fPtNch(nullptr),
    dPtNch(nullptr),
    dPt2Nch(nullptr),
    dPt3Nch(nullptr),
    dPtNch_0(nullptr),
    dPtNch_1(nullptr),
    dPtNch_2(nullptr),
    dPtNch_3(nullptr),
    dPtNch_4(nullptr),
    dPtNch_5(nullptr),
    dPtNch_6(nullptr),
    dPtNch_7(nullptr),
    dPtNch_8(nullptr),
    dPtNch_9(nullptr),
    dPt2Nch_0(nullptr),
    dPt2Nch_1(nullptr),
    dPt2Nch_2(nullptr),
    dPt2Nch_3(nullptr),
    dPt2Nch_4(nullptr),
    dPt2Nch_5(nullptr),
    dPt2Nch_6(nullptr),
    dPt2Nch_7(nullptr),
    dPt2Nch_8(nullptr),
    dPt2Nch_9(nullptr),
    dPt3Nch_0(nullptr),
    dPt3Nch_1(nullptr),
    dPt3Nch_2(nullptr),
    dPt3Nch_3(nullptr),
    dPt3Nch_4(nullptr),
    dPt3Nch_5(nullptr),
    dPt3Nch_6(nullptr),
    dPt3Nch_7(nullptr),
    dPt3Nch_8(nullptr),
    dPt3Nch_9(nullptr),
    fPtNchUCC(nullptr),
    TestPtCtr(nullptr),
    TestPt2Ctr(nullptr),
    TestPt3Ctr(nullptr),
    TestNchCtr(nullptr),
    TestNchSelectedCtr(nullptr),
    fEventCuts(0),
    multSelection(nullptr),
    radm(32213)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisPtN::AliAnalysisPtN(const char* name) : AliAnalysisTaskSE(name),
    fAOD(nullptr), 
    fOutputList(nullptr), 
    fBstList(nullptr),
    fWeightsListNUE(nullptr),
    fWeightNUE(nullptr),
    //correlator(),
    fTestNonWeight(nullptr),
    fPtNch(nullptr),
    dPtNch(nullptr),
    dPt2Nch(nullptr),
    dPt3Nch(nullptr),
    dPtNch_0(nullptr),
    dPtNch_1(nullptr),
    dPtNch_2(nullptr),
    dPtNch_3(nullptr),
    dPtNch_4(nullptr),
    dPtNch_5(nullptr),
    dPtNch_6(nullptr),
    dPtNch_7(nullptr),
    dPtNch_8(nullptr),
    dPtNch_9(nullptr),
    dPt2Nch_0(nullptr),
    dPt2Nch_1(nullptr),
    dPt2Nch_2(nullptr),
    dPt2Nch_3(nullptr),
    dPt2Nch_4(nullptr),
    dPt2Nch_5(nullptr),
    dPt2Nch_6(nullptr),
    dPt2Nch_7(nullptr),
    dPt2Nch_8(nullptr),
    dPt2Nch_9(nullptr),
    dPt3Nch_0(nullptr),
    dPt3Nch_1(nullptr),
    dPt3Nch_2(nullptr),
    dPt3Nch_3(nullptr),
    dPt3Nch_4(nullptr),
    dPt3Nch_5(nullptr),
    dPt3Nch_6(nullptr),
    dPt3Nch_7(nullptr),
    dPt3Nch_8(nullptr),
    dPt3Nch_9(nullptr),
    fPtNchUCC(nullptr),
    TestPtCtr(nullptr),
    TestPt2Ctr(nullptr),
    TestPt3Ctr(nullptr),
    TestNchCtr(nullptr),
    TestNchSelectedCtr(nullptr),
    fEventCuts(0),
    multSelection(nullptr),
    radm(32213)
{
    // constructor
    
    DefineInput(0, TChain::Class()); 
    DefineInput(1, TList::Class()); 
      // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisPtN::~AliAnalysisPtN()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisPtN::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
    fEventCuts.AddQAplotsToList(fOutputList);
    
    fTestNonWeight = new TH1F("fHisPhi", "fHisPhi", 100, 0, 6.28);
    fPtNch = new TH2F("fPtNch", "fPtNch", 900, 0, 4500, 500, 0, 3);
    dPtNch = new TProfile("dPtNch", "dPtNch", 90, 0, 4500);
    dPt2Nch = new TProfile("dPt2Nch", "dPt2Nch", 90, 0, 4500);
    dPt3Nch = new TProfile("dPt3Nch", "dPt3Nch", 90, 0, 4500);
    fPtNchUCC = new TH2F("fPtNchUCC", "fPtNchUCC", 600, 1800, 3000, 500, 0.2, 5.0);
    TestPtCtr = new TProfile("TestPtCtr", "TestPtCtr", 100, 0, 100);
    TestPt2Ctr = new TProfile("TestPt2Ctr", "TestPt2Ctr", 100, 0, 100);
    TestPt3Ctr = new TProfile("TestPt3Ctr", "TestPt3Ctr", 100, 0, 100);
    TestNchCtr = new TProfile("TestNchCtr", "TestNchCtr", 100, 0, 100);
    TestNchSelectedCtr = new TProfile("TestNchSelectedCtr", "TestNchSelectedCtr", 100, 0, 100);

    fBstList = new TList();
    dPtNch_0 = new TProfile("dPtNch_0", "dPtNch_0", 90, 0, 4500);
    dPtNch_1 = new TProfile("dPtNch_1", "dPtNch_1", 90, 0, 4500);
    dPtNch_2 = new TProfile("dPtNch_2", "dPtNch_2", 90, 0, 4500);
    dPtNch_3 = new TProfile("dPtNch_3", "dPtNch_3", 90, 0, 4500);
    dPtNch_4 = new TProfile("dPtNch_4", "dPtNch_4", 90, 0, 4500);
    dPtNch_5 = new TProfile("dPtNch_5", "dPtNch_5", 90, 0, 4500);
    dPtNch_6 = new TProfile("dPtNch_6", "dPtNch_6", 90, 0, 4500);
    dPtNch_7 = new TProfile("dPtNch_7", "dPtNch_7", 90, 0, 4500);
    dPtNch_8 = new TProfile("dPtNch_8", "dPtNch_8", 90, 0, 4500);
    dPtNch_9 = new TProfile("dPtNch_9", "dPtNch_9", 90, 0, 4500);
    
    dPt2Nch_0 = new TProfile("dPt2Nch_0", "dPt2Nch_0", 90, 0, 4500);
    dPt2Nch_1 = new TProfile("dPt2Nch_1", "dPt2Nch_1", 90, 0, 4500);
    dPt2Nch_2 = new TProfile("dPt2Nch_2", "dPt2Nch_2", 90, 0, 4500);
    dPt2Nch_3 = new TProfile("dPt2Nch_3", "dPt2Nch_3", 90, 0, 4500);
    dPt2Nch_4 = new TProfile("dPt2Nch_4", "dPt2Nch_4", 90, 0, 4500);
    dPt2Nch_5 = new TProfile("dPt2Nch_5", "dPt2Nch_5", 90, 0, 4500);
    dPt2Nch_6 = new TProfile("dPt2Nch_6", "dPt2Nch_6", 90, 0, 4500);
    dPt2Nch_7 = new TProfile("dPt2Nch_7", "dPt2Nch_7", 90, 0, 4500);
    dPt2Nch_8 = new TProfile("dPt2Nch_8", "dPt2Nch_8", 90, 0, 4500);
    dPt2Nch_9 = new TProfile("dPt2Nch_9", "dPt2Nch_9", 90, 0, 4500);
  
    dPt3Nch_0 = new TProfile("dPt3Nch_0", "dPt3Nch_0", 90, 0, 4500);
    dPt3Nch_1 = new TProfile("dPt3Nch_1", "dPt3Nch_1", 90, 0, 4500);
    dPt3Nch_2 = new TProfile("dPt3Nch_2", "dPt3Nch_2", 90, 0, 4500);
    dPt3Nch_3 = new TProfile("dPt3Nch_3", "dPt3Nch_3", 90, 0, 4500);
    dPt3Nch_4 = new TProfile("dPt3Nch_4", "dPt3Nch_4", 90, 0, 4500);
    dPt3Nch_5 = new TProfile("dPt3Nch_5", "dPt3Nch_5", 90, 0, 4500);
    dPt3Nch_6 = new TProfile("dPt3Nch_6", "dPt3Nch_6", 90, 0, 4500);
    dPt3Nch_7 = new TProfile("dPt3Nch_7", "dPt3Nch_7", 90, 0, 4500);
    dPt3Nch_8 = new TProfile("dPt3Nch_8", "dPt3Nch_8", 90, 0, 4500);
    dPt3Nch_9 = new TProfile("dPt3Nch_9", "dPt3Nch_9", 90, 0, 4500);
     // create your histogram
    fBstList->Add(dPtNch_0);
    fBstList->Add(dPtNch_1);
    fBstList->Add(dPtNch_2);
    fBstList->Add(dPtNch_3);
    fBstList->Add(dPtNch_4);
    fBstList->Add(dPtNch_5);
    fBstList->Add(dPtNch_6);
    fBstList->Add(dPtNch_7);
    fBstList->Add(dPtNch_8);
    fBstList->Add(dPtNch_9);
    fBstList->Add(dPt2Nch_0);
    fBstList->Add(dPt2Nch_1);
    fBstList->Add(dPt2Nch_2);
    fBstList->Add(dPt2Nch_3);
    fBstList->Add(dPt2Nch_4);
    fBstList->Add(dPt2Nch_5);
    fBstList->Add(dPt2Nch_6);
    fBstList->Add(dPt2Nch_7);
    fBstList->Add(dPt2Nch_8);
    fBstList->Add(dPt2Nch_9);
    fBstList->Add(dPt3Nch_0);
    fBstList->Add(dPt3Nch_1);
    fBstList->Add(dPt3Nch_2);
    fBstList->Add(dPt3Nch_3);
    fBstList->Add(dPt3Nch_4);
    fBstList->Add(dPt3Nch_5);
    fBstList->Add(dPt3Nch_6);
    fBstList->Add(dPt3Nch_7);
    fBstList->Add(dPt3Nch_8);
    fBstList->Add(dPt3Nch_9);

    fOutputList->Add(fBstList);
    fOutputList->Add(fTestNonWeight);
    fOutputList->Add(fPtNch);
    fOutputList->Add(dPtNch);
    fOutputList->Add(dPt2Nch);
    fOutputList->Add(dPt3Nch);
    fOutputList->Add(fPtNchUCC);
    fOutputList->Add(TestPtCtr);
    fOutputList->Add(TestPt2Ctr);
    fOutputList->Add(TestPt3Ctr);
    fOutputList->Add(TestNchCtr);
    fOutputList->Add(TestNchSelectedCtr);
         // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
    
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisPtN::UserExec(Option_t *)
{
    
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar way
    // cout<< "event inputted"<<endl;                                                    // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;
    // cout<< "event is AOD"<<endl;
    

    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isTrigselected = false;
    isTrigselected = fSelectMask&(AliVEvent::kINT7+AliVEvent::kMB);
    if(isTrigselected == false) return;
    
    if(!fEventCuts.AcceptEvent(fAOD)) return;                                  // if the pointer to the event is empty (getting it failed) skip this event
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram
    // cout<< "event selected"<<endl;
    Float_t vertexZ = fAOD->GetPrimaryVertex()->GetZ();
    if(vertexZ<=(-10) || vertexZ>=10) return;             //filter Vz range

    Float_t centrality(0);
    AliMultSelection* multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");
    //if(centrality>90) return;

    //initilize the weight list
    fWeightsListNUE = (TList*)GetInputData(1);
    fWeightNUE = (TH1D*)fWeightsListNUE->FindObject(Form("EffRescaled_Cent0"));

    Float_t wtE;

    Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    Double_t sump = 0, sumw = 0, sump2 = 0, sumw2 = 0, sump3 = 0, sumw3 = 0;
    Int_t nTrackSelected = 0;
                                            // how many tracks are selected
    for(Int_t i(0); i < iTracks; i++) {
                                                            // loop ove rall these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track || !track->TestFilterBit(96)) continue; // filterbit if we failed, skip this track
        if(track->Pt()<=0.2 || track->Pt()>=3.0) continue; //filter Pt range                           
        if(track->Eta()<=(-0.8) || track->Eta()>=0.8) continue;//filter Eta range
        ++nTrackSelected;
        wtE = GetWeightNUE(track->Pt());
        fTestNonWeight->Fill(track->Phi());


          sump += wtE*track->Pt();
          sump2 += wtE*wtE*track->Pt()*track->Pt();
          sump3 += wtE*wtE*wtE*track->Pt()*track->Pt()*track->Pt();

          sumw += wtE;
          sumw2 += wtE*wtE;
          sumw3 += wtE*wtE*wtE;
    
    }
    
    //if(M<100) return; //we need enough tracks in one event
    if(sumw!=0 && (sumw*sumw-sumw2)!=0 && (sumw*sumw*sumw-3*sumw2*sumw+2*sumw3)!=0) {
      Float_t pt = sump/sumw;
      fPtNch->Fill(nTrackSelected,pt);
      dPtNch->Fill(nTrackSelected,pt);

      Float_t pt2 = (sump*sump-sump2)/(sumw*sumw-sumw2);
      dPt2Nch->Fill(nTrackSelected,pt2);

      Float_t pt3 = (sump*sump*sump-3*sump2*sump+2*sump3)/(sumw*sumw*sumw-3*sumw2*sumw+2*sumw3);
      dPt3Nch->Fill(nTrackSelected,pt3);


      if(centrality<1) fPtNchUCC->Fill(nTrackSelected,pt);
      TestPtCtr->Fill(centrality,pt);
      TestPt2Ctr->Fill(centrality,pt2);
      TestPt3Ctr->Fill(centrality,pt3);

    //here we fill the boostrap profiles
      Int_t rd = int(floor(radm.Rndm()*10));
      switch (rd){
        case 0:
          dPtNch_0->Fill(nTrackSelected,pt);
          dPt2Nch_0->Fill(nTrackSelected,pt2);
          dPt3Nch_0->Fill(nTrackSelected,pt3);
          break;
        case 1:
          dPtNch_1->Fill(nTrackSelected,pt);
          dPt2Nch_1->Fill(nTrackSelected,pt2);
          dPt3Nch_1->Fill(nTrackSelected,pt3);
          break;
        case 2:
          dPtNch_2->Fill(nTrackSelected,pt);
          dPt2Nch_2->Fill(nTrackSelected,pt2);
          dPt3Nch_2->Fill(nTrackSelected,pt3);
          break;
        case 3:
          dPtNch_3->Fill(nTrackSelected,pt);
          dPt2Nch_3->Fill(nTrackSelected,pt2);
          dPt3Nch_3->Fill(nTrackSelected,pt3);
          break;
        case 4:
          dPtNch_4->Fill(nTrackSelected,pt);
          dPt2Nch_4->Fill(nTrackSelected,pt2);
          dPt3Nch_4->Fill(nTrackSelected,pt3);
          break;
        case 5:
          dPtNch_5->Fill(nTrackSelected,pt);
          dPt2Nch_5->Fill(nTrackSelected,pt2);
          dPt3Nch_5->Fill(nTrackSelected,pt3);
          break;
        case 6:
          dPtNch_6->Fill(nTrackSelected,pt);
          dPt2Nch_6->Fill(nTrackSelected,pt2);
          dPt3Nch_6->Fill(nTrackSelected,pt3);
          break;
        case 7:
          dPtNch_7->Fill(nTrackSelected,pt);
          dPt2Nch_7->Fill(nTrackSelected,pt2);
          dPt3Nch_7->Fill(nTrackSelected,pt3);
          break;
        case 8:
          dPtNch_8->Fill(nTrackSelected,pt);
          dPt2Nch_8->Fill(nTrackSelected,pt2);
          dPt3Nch_8->Fill(nTrackSelected,pt3);
          break;
        case 9:
          dPtNch_9->Fill(nTrackSelected,pt);
          dPt2Nch_9->Fill(nTrackSelected,pt2);
          dPt3Nch_9->Fill(nTrackSelected,pt3);
          break;
      }
    }
    
    
    
    TestNchCtr->Fill(centrality,iTracks);
    TestNchSelectedCtr->Fill(centrality,nTrackSelected);
    

    
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}
//____________________________________________________________________________
// double AliAnalysisPtN::GetWeightNUA(double phi, double eta, double vz) {
//   double weight = fWeightNUA->GetBinContent(fWeightNUA->GetXaxis()->FindBin(phi),
//       fWeightNUA->GetYaxis()->FindBin(eta),
//       fWeightNUA->GetZaxis()->FindBin(vz));
//   return weight;
// }
//_____________________________________________________________________________
double AliAnalysisPtN::GetWeightNUE(double pt)
{
  double binPt = fWeightNUE->GetXaxis()->FindBin(pt);
  double eff = fWeightNUE->GetBinContent(binPt);
  double error = fWeightNUE->GetBinError(binPt);
  double weight = 1;
  //..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error


  if((eff < 0.03) || ((error/eff) > 0.1)) error = 0.00001;
  if((eff < 0.03)) return 1;

  TRandom3 r(0);
  double efficiency = 0;
  efficiency = r.Gaus(eff, error);
  weight = 1./efficiency; 
  return weight;
}
//_____________________________________________________________________________
void AliAnalysisPtN::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
