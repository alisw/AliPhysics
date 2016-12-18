#include "AliAnalysisTaskEventCutsValidation.h"

// ROOT includes
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>

// ALIROOT includes
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskEventCutsValidation);
///\endcond

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskEventCutsValidation::AliAnalysisTaskEventCutsValidation(bool storeCuts, TString taskname) :
  AliAnalysisTaskSE(taskname.Data()),
  fEventCut(false),
  fFB32cuts(),
  fTPConlyCuts(),
  fMultTOFLowCut(0x0),
  fMultTOFHighCut(0x0),
  fMultCentLowCut(0x0),
  fList(nullptr),
  fStoreCuts(storeCuts)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

  AliESDtrackCuts* tmpFB32 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  AliESDtrackCuts* tmpTPConly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fFB32cuts    = *tmpFB32;
  fTPConlyCuts = *tmpTPConly;
  delete tmpFB32;
  delete tmpTPConly;

  if (storeCuts) DefineOutput(2, AliNuclexEventCuts::Class());
}

/// Standard destructor
///
AliAnalysisTaskEventCutsValidation::~AliAnalysisTaskEventCutsValidation() {
  if (fList) delete fList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskEventCutsValidation::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);
  fEventCut.AddQAplotsToList(fList);


  fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000),
                 fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000),
                 fMultCentLowCut = new TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100),

                 fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
  fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
  fMultCentLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);

  string label[2] = {"raw","selected"};
  for (int iS = 0; iS < 2; ++iS) {
    fTOFvsFB32[iS] = new TH2F(Form("fTOFvsFB32%s",label[iS].data()),Form("%s;Multiplicity FB32;Multiplicity FB32+TOF",label[iS].data()),500,0.,3000.,500,0.,2000.);
    fTPCvsAll[iS] = new TH2F(Form("fTPCvsAll%s",label[iS].data()),Form("%s;Multiplicity TPC;Multiplicity ESD - a_{0} #times Multiplicity TPC",label[iS].data()),300,0.,3000.,3000,-200.,28800.);
    fMultvsV0M[iS] = new TH2F(Form("fMultvsV0M%s",label[iS].data()),Form("%s;Centrality (V0M);FB32",label[iS].data()),100,0.,100.,500,0.,3000.);
    fList->Add(fTOFvsFB32[iS]);
    fList->Add(fTPCvsAll[iS]);
    fList->Add(fMultvsV0M[iS]);
  }

  PostData(1,fList);
  if (fStoreCuts) PostData(2,&fEventCut);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskEventCutsValidation::UserExec(Option_t *) {

  AliVEvent* ev = InputEvent();
  if (fEventCut.AcceptEvent(ev)) {
    float multESDTPCDif = float(fEventCut.fContainer.fMultESD) - 3.38f * fEventCut.fContainer.fMultTrkTPC;
    fTOFvsFB32[0]->Fill(fEventCut.fContainer.fMultTrkFB32,fEventCut.fContainer.fMultTrkFB32TOF);
    fTPCvsAll[0]->Fill(fEventCut.fContainer.fMultTrkTPC,multESDTPCDif);
    fMultvsV0M[0]->Fill(fEventCut.GetCentrality(),fEventCut.fContainer.fMultTrkFB32Acc);

    if (float(fEventCut.fContainer.fMultTrkFB32TOF) > fMultTOFLowCut->Eval(float(fEventCut.fContainer.fMultTrkFB32)) &&
        float(fEventCut.fContainer.fMultTrkFB32TOF) < fMultTOFHighCut->Eval(float(fEventCut.fContainer.fMultTrkFB32)) &&
        multESDTPCDif < 15000.f &&
        float(fEventCut.fContainer.fMultTrkFB32Acc) > fMultCentLowCut->Eval(fEventCut.GetCentrality())) {
      fTOFvsFB32[1]->Fill(fEventCut.fContainer.fMultTrkFB32,fEventCut.fContainer.fMultTrkFB32TOF);
      fTPCvsAll[1]->Fill(fEventCut.fContainer.fMultTrkTPC,multESDTPCDif);
      fMultvsV0M[1]->Fill(fEventCut.GetCentrality(),fEventCut.fContainer.fMultTrkFB32Acc);
    }
  }

  PostData(1,fList);
  if (fStoreCuts) PostData(2,&fEventCut);
  return;
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskEventCutsValidation::Terminate(Option_t *) {
  return;
}

