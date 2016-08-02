///
/// \file AliFemtoCutMonitorEventMult.cxx
///

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult():
  fEvMult(NULL),
  fNormEvMult(NULL),
  fSPDMult(NULL),
  fMultSumPt(NULL),
  freadMC(kFALSE),
  faddhists(kFALSE),
  fEstimateITSTPC(NULL),
  fEstimateTracklets(NULL),
  fEstimateITSPure(NULL),
  fEst1Est2(NULL),
  fEst1Est3(NULL),
  fEst2Est3(NULL),
  fEst1Norm(NULL),
  fEst2Norm(NULL),
  fEst3Norm(NULL),
  fPsiVZERO(NULL)
{
  // Default constructor
  fEvMult = new TH1D("EvMult", "Event Multiplicity", 5001, -0.5, 5000.5);
  fMultSumPt = new TH2D("EvMultSumPt",
                        "Event Multiplicity vs Total pT",
                        5001, -0.5, 5000.5,
                        1000, 0.0, 100.0);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const char *aName, int nBins, double multMax):
  AliFemtoCutMonitor(),
  fEvMult(NULL),
  fNormEvMult(NULL),
  fSPDMult(NULL),
  fMultSumPt(NULL),
  freadMC(kFALSE),
  faddhists(kFALSE),
  fEstimateITSTPC(NULL),
  fEstimateTracklets(NULL),
  fEstimateITSPure(NULL),
  fEst1Est2(NULL),
  fEst1Est3(NULL),
  fEst2Est3(NULL),
  fEst1Norm(NULL),
  fEst2Norm(NULL),
  fEst3Norm(NULL),
  fPsiVZERO(NULL)
{
  TString name(aName);

  // Normal constructor
  fEvMult = new TH1D("EvMult" + name,
                     "Event Multiplicity",
                     nBins+1, -0.5, multMax);

  fNormEvMult = new TH1D("NormEvMult" + name,
                         "Normalized Event Multiplicity",
                         nBins+1, -0.5, multMax);

  if (!freadMC) {
    fSPDMult = new TH1D("SPDEvMult" + name,
                        "SPD Tracklet Multiplicity",
                        nBins+1, -0.5,  multMax);
  }

  fMultSumPt = new TH2D("EvMultTotPt" + name,
                        "Event Multiplicity vs Total pT",
                        501, -0.5, 500.5,
                        1000, 0.0, 100.0);



  if (faddhists) {
    fEstimateITSTPC = new TH1D("EvMultEstITSTPC" + name,
                               "ITS+TPC Multiplicity Estimate",
                               5001, -0.5, 5000.5);

    fEstimateTracklets = new TH1D("EvMultEstTracklets" + name,
                                  "Tracklets Multiplicity Estimate",
                                  5001, -0.5, 5000.5);

    fEstimateITSPure = new TH1D("EvMultEstITSPure" + name,
                                "ITS Pure Multiplicity Estimate",
                                8001, -0.5, 8000.5);

    fEst1Est2 = new TH2D("EstITSTPCEstTracklet" + name,
                         "ITS+TPC vs Tracklets",
                         501, -0.5, 5000.5,
                         501, -0.5, 500.5);

    fEst1Est3 = new TH2D("EstITSTPCEstITSPure" + name,
                         "ITS+TPC vs ITS Pure",
                         501, -0.5, 5000.5,
                         801, -0.5, 8000.5);

    fEst2Est3 = new TH2D("EstTrackletEstITSPure" + name,
                         "Tracklets vs ITS Pure",
                         501, -0.5, 5000.5,
                         801, -0.5, 8000.5);

    fEst1Norm = new TH2D("EstITSTPCNormMult" + name,
                         "ITS+TPC vs Normalized Mult",
                         501, -0.5, 5000.5,
                         501, -0.5, 5000.5);

    fEst2Norm = new TH2D("EstTrackletsNormMult" + name,
                         "Tracklets vs Normalized Mult",
                         501, -0.5, 5000.5,
                         501, -0.5, 5000.5);

    fEst3Norm = new TH2D("EstITSPureNormMult" + name,
                         "ITS Pure vs Normalized Mult",
                         501, -0.5, 5000.5,
                         501, -0.5, 5000.5);
  }

  fPsiVZERO = new TH1D("PsiEPVZERO" + name,
                       "event plane angle from vzero",
                       157, -1.575, 1.565);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const AliFemtoCutMonitorEventMult &aCut):
  AliFemtoCutMonitor(aCut),
  fEvMult(NULL),
  fNormEvMult(NULL),
  fSPDMult(NULL),
  fMultSumPt(NULL),
  freadMC(aCut.freadMC),
  faddhists(aCut.faddhists),
  fEstimateITSTPC(NULL),
  fEstimateTracklets(NULL),
  fEstimateITSPure(NULL),
  fEst1Est2(NULL),
  fEst1Est3(NULL),
  fEst2Est3(NULL),
  fEst1Norm(NULL),
  fEst2Norm(NULL),
  fEst3Norm(NULL),
  fPsiVZERO(NULL)
{
  // copy constructor
  fEvMult = new TH1D(*aCut.fEvMult);
  fNormEvMult = new TH1D(*aCut.fNormEvMult);

  if(!freadMC) {
    fSPDMult = new TH1D(*aCut.fSPDMult);
  }

  fMultSumPt = new TH2D(*aCut.fMultSumPt);

  if (faddhists) {
    fEstimateITSTPC = new TH1D(*aCut.fEstimateITSTPC);
    fEstimateTracklets = new TH1D(*aCut.fEstimateTracklets);
    fEstimateITSPure = new TH1D(*aCut.fEstimateITSPure);
    fEst1Est2 = new TH2D(*aCut.fEst1Est2);
    fEst1Est3 = new TH2D(*aCut.fEst1Est3);
    fEst2Est3 = new TH2D(*aCut.fEst2Est3);
    fEst1Norm = new TH2D(*aCut.fEst1Norm);
    fEst2Norm = new TH2D(*aCut.fEst2Norm);
    fEst3Norm = new TH2D(*aCut.fEst3Norm);
  }

  fPsiVZERO = new TH1D(*aCut.fPsiVZERO);
}

AliFemtoCutMonitorEventMult::~AliFemtoCutMonitorEventMult()
{
  // Destructor
  delete fEvMult;
  delete fNormEvMult;

  if (!freadMC) {
    delete fSPDMult;
  }

  delete fMultSumPt;

  if (faddhists) {
    delete fEstimateITSTPC;
    delete fEstimateTracklets;
    delete fEstimateITSPure;
    delete fEst1Est2;
    delete fEst1Est3;
    delete fEst2Est3;
    delete fEst1Norm;
    delete fEst2Norm;
    delete fEst3Norm;
  }

  delete fPsiVZERO;
}

AliFemtoCutMonitorEventMult& AliFemtoCutMonitorEventMult::operator=(const AliFemtoCutMonitorEventMult& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  if (fEvMult) delete fEvMult;
  fEvMult = new TH1D(*aCut.fEvMult);

  if (fNormEvMult) delete fNormEvMult;
  fNormEvMult = new TH1D(*aCut.fNormEvMult);

  if (fPsiVZERO) delete fPsiVZERO;
  fPsiVZERO = new TH1D(*aCut.fPsiVZERO);

  if (!freadMC) {
    if (fSPDMult) delete fSPDMult;
    fSPDMult = new TH1D(*aCut.fSPDMult);
  }

  if (fMultSumPt) delete fMultSumPt;
  fMultSumPt = new TH2D(*aCut.fMultSumPt);


  if (faddhists) {
    if (fEstimateITSTPC) delete fEstimateITSTPC;
    fEstimateITSTPC = new TH1D(*aCut.fEstimateITSTPC);

    if (fEstimateTracklets) delete fEstimateTracklets;
    fEstimateTracklets = new TH1D(*aCut.fEstimateTracklets);

    if (fEstimateITSPure) delete fEstimateITSPure;
    fEstimateITSPure = new TH1D(*aCut.fEstimateITSPure);

    if (fEst1Est2) delete fEst1Est2;
    fEst1Est2 = new TH2D(*aCut.fEst1Est2);

    if (fEst1Est3) delete fEst1Est3;
    fEst1Est3 = new TH2D(*aCut.fEst1Est3);

    if (fEst2Est3) delete fEst2Est3;
    fEst2Est3 = new TH2D(*aCut.fEst2Est3);

    if (fEst1Norm) delete fEst1Norm;
    fEst1Norm = new TH2D(*aCut.fEst1Norm);

    if (fEst2Norm) delete fEst2Norm;
    fEst2Norm = new TH2D(*aCut.fEst2Norm);

    if (fEst3Norm) delete fEst3Norm;
    fEst3Norm = new TH2D(*aCut.fEst3Norm);
  }

  return *this;
}

AliFemtoString AliFemtoCutMonitorEventMult::Report()
{
  // Prepare report from the execution
  TString report = "*** AliFemtoCutMonitorEventMult report";
  return AliFemtoString(report);
}

void AliFemtoCutMonitorEventMult::Fill(const AliFemtoEvent* aEvent)
{
  // Fill in the monitor histograms with the values from the current track
  fEvMult->Fill(aEvent->NumberOfTracks());
  fNormEvMult->Fill(aEvent->UncorrectedNumberOfPrimaries());

  double epvzero = aEvent->ReactionPlaneAngle();

  fPsiVZERO->Fill(epvzero);


  // if(!freadMC){
  //   fSPDMult->Fill(aEvent->SPDMultiplicity());
  // }
  // fMultSumPt->Fill(aEvent->UncorrectedNumberOfPrimaries(), aEvent->ZDCEMEnergy());

  // if(faddhists)
  //   {
  //     fEstimateITSTPC->Fill(aEvent->MultiplicityEstimateITSTPC());
  //     fEstimateTracklets->Fill(aEvent->MultiplicityEstimateTracklets());
  //     fEstimateITSPure->Fill(aEvent->MultiplicityEstimateITSPure());
  //     fEst1Est2->Fill(aEvent->MultiplicityEstimateITSTPC(),aEvent->MultiplicityEstimateTracklets());
  //     fEst1Est3->Fill(aEvent->MultiplicityEstimateITSTPC(),aEvent->MultiplicityEstimateITSPure());
  //     fEst2Est3->Fill(aEvent->MultiplicityEstimateTracklets(),aEvent->MultiplicityEstimateITSPure());
  //     fEst1Norm->Fill(aEvent->MultiplicityEstimateITSTPC(),aEvent->UncorrectedNumberOfPrimaries());
  //     fEst2Norm->Fill(aEvent->MultiplicityEstimateTracklets(),aEvent->UncorrectedNumberOfPrimaries());
  //     fEst3Norm->Fill(aEvent->MultiplicityEstimateITSPure(),aEvent->UncorrectedNumberOfPrimaries());
  //   }

}

void AliFemtoCutMonitorEventMult::Write()
{
  // Write out the relevant histograms
  fEvMult->Write();
  fNormEvMult->Write();
  fPsiVZERO->Write();

  // if(!freadMC){
  //   fSPDMult->Write();
  // }
  // fMultSumPt->Write();

  // if(faddhists)
  //   {
  //     fEstimateITSTPC->Write();
  //     fEstimateTracklets->Write();
  //     fEstimateITSPure->Write();
  //     fEst1Est2->Write();
  //     fEst1Est3->Write();
  //     fEst2Est3->Write();
  //     fEst1Norm->Write();
  //     fEst2Norm->Write();
  //     fEst3Norm->Write();
  //   }

}

TList *AliFemtoCutMonitorEventMult::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fEvMult);
  tOutputList->Add(fNormEvMult);
  tOutputList->Add(fPsiVZERO);
  // tOutputList->Add(fSPDMult);
  // tOutputList->Add(fMultSumPt);

  // if(faddhists)
  //   {
  //     tOutputList->Add(fEstimateITSTPC);
  //     tOutputList->Add(fEstimateTracklets);
  //     tOutputList->Add(fEstimateITSPure);
  //     tOutputList->Add(fEst1Est2);
  //     tOutputList->Add(fEst1Est3);
  //     tOutputList->Add(fEst2Est3);
  //     tOutputList->Add(fEst1Norm);
  //     tOutputList->Add(fEst2Norm);
  //     tOutputList->Add(fEst3Norm);
  //   }
  return tOutputList;
}

void AliFemtoCutMonitorEventMult::SetReadMC(Bool_t mc)
{
  freadMC = mc;
}

void AliFemtoCutMonitorEventMult::AdditionalMultHistsOn(Bool_t addhists)
{
  faddhists = addhists;
}
