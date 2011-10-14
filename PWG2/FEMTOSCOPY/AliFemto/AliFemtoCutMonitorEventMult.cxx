////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorEventMult - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult():
  fEvMult(0),
  fNormEvMult(0),
  fSPDMult(0),
  fMultSumPt(0),
  fEstimateITSTPC(0),
  fEstimateTracklets(0),
  fEstimateITSPure(0),
  fEst1Est2(0),
  fEst1Est3(0),
  fEst2Est3(0),
  fEst1Norm(0),
  fEst2Norm(0),
  fEst3Norm(0),
  freadMC(kFALSE),
  faddhists(kFALSE)
{
  // Default constructor
  fEvMult = new TH1D("EvMult", "Event Multiplicity", 5001, -0.5, 5000.5);
  fMultSumPt = new TH2D("EvMultSumPt","Event Multiplicity vs Total pT",5001,-0.5,5000.5,1000,0.0,100.0);
}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const char *aName):
  AliFemtoCutMonitor(),
  fEvMult(0),
  fNormEvMult(0),
  fSPDMult(0),
  fMultSumPt(0),
  fEstimateITSTPC(0),
  fEstimateTracklets(0),
  fEstimateITSPure(0),
  fEst1Est2(0),
  fEst1Est3(0),
  fEst2Est3(0),
  fEst1Norm(0),
  fEst2Norm(0),
  fEst3Norm(0),
  freadMC(kFALSE),
  faddhists(kFALSE)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "EvMult%s", aName);
  fEvMult = new TH1D(name, "Event Multiplicity", 5001, -0.5, 5000.5);

  snprintf(name, 200, "NormEvMult%s", aName);
  fNormEvMult = new TH1D(name, "Normalized Event Multiplicity", 5001, -0.5, 5000.5);

  if(!freadMC) {
    snprintf(name, 200, "SPDEvMult%s", aName);
    fSPDMult = new TH1D(name, "SPD Tracklet Multiplicity", 5001, -0.5, 5000.5);
  }

  snprintf(name, 200, "EvMultTotPt%s", aName);
  fMultSumPt = new TH2D(name,"Event Multiplicity vs Total pT",501,-0.5,500.5,1000,0.0,100.0);

  if(faddhists)
    {
      snprintf(name, 200, "EvMultEstITSTPC%s", aName);
      fEstimateITSTPC = new TH1D(name, "ITS+TPC Multiplicity Estimate", 5001, -0.5, 5000.5);

      snprintf(name, 200, "EvMultEstTracklets%s", aName);
      fEstimateTracklets = new TH1D(name, "Tracklets Multiplicity Estimate", 5001, -0.5, 5000.5);

      snprintf(name, 200, "EvMultEstITSPure%s", aName);
      fEstimateITSPure = new TH1D(name, "ITS Pure Multiplicity Estimate", 8001, -0.5, 8000.5);

      snprintf(name, 200, "EstITSTPCEstTracklet%s", aName);
      fEst1Est2 = new TH2D(name,"ITS+TPC vs Tracklets",501,-0.5,5000.5,501,-0.5,500.5);

      snprintf(name, 200, "EstITSTPCEstITSPure%s", aName);
      fEst1Est3 = new TH2D(name,"ITS+TPC vs ITS Pure",501,-0.5,5000.5,801,-0.5,8000.5);

      snprintf(name, 200, "EstTrackletEstITSPure%s", aName);
      fEst2Est3 = new TH2D(name,"Tracklets vs ITS Pure",501,-0.5,5000.5,801,-0.5,8000.5);

      snprintf(name, 200, "EstITSTPCNormMult%s", aName);
      fEst1Norm = new TH2D(name,"ITS+TPC vs Normalized Mult",501,-0.5,5000.5,501,-0.5,5000.5);

      snprintf(name, 200, "EstTrackletsNormMult%s", aName);
      fEst2Norm = new TH2D(name,"Tracklets vs Normalized Mult",501,-0.5,5000.5,501,-0.5,5000.5);

      snprintf(name, 200, "EstITSPureNormMult%s", aName);
      fEst3Norm = new TH2D(name,"ITS Pure vs Normalized Mult",501,-0.5,5000.5,501,-0.5,5000.5);
    }

}

AliFemtoCutMonitorEventMult::AliFemtoCutMonitorEventMult(const AliFemtoCutMonitorEventMult &aCut):
  AliFemtoCutMonitor(),
  fEvMult(0),
  fNormEvMult(0),
  fSPDMult(0),
  fMultSumPt(0),
  fEstimateITSTPC(0),
  fEstimateTracklets(0),
  fEstimateITSPure(0),
  fEst1Est2(0),
  fEst1Est3(0),
  fEst2Est3(0),
  fEst1Norm(0),
  fEst2Norm(0),
  fEst3Norm(0),
  freadMC(kFALSE),
  faddhists(kFALSE)
{
  // copy constructor
  if (fEvMult) delete fEvMult;
  fEvMult = new TH1D(*aCut.fEvMult);

  if (fNormEvMult) delete fNormEvMult;
  fNormEvMult = new TH1D(*aCut.fNormEvMult);

  
  if(!freadMC){
    if (fSPDMult) delete fSPDMult;
    fSPDMult = new TH1D(*aCut.fSPDMult);
  }

  if (fMultSumPt) delete fMultSumPt;
  fMultSumPt = new TH2D(*aCut.fMultSumPt);


  if(faddhists)
    {
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
}

AliFemtoCutMonitorEventMult::~AliFemtoCutMonitorEventMult()
{
  // Destructor
  delete fEvMult;
  delete fNormEvMult;
  if(!freadMC){
    delete fSPDMult;
  }
  delete fMultSumPt;

  if(faddhists)
    {
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
  
  if(!freadMC){
    if (fSPDMult) delete fSPDMult;
    fSPDMult = new TH1D(*aCut.fSPDMult);
  }
  
  if (fMultSumPt) delete fMultSumPt;
  fMultSumPt = new TH2D(*aCut.fMultSumPt);


  if(faddhists)
    {
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

AliFemtoString AliFemtoCutMonitorEventMult::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorEventMult report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorEventMult::Fill(const AliFemtoEvent* aEvent)
{
  // Fill in the monitor histograms with the values from the current track
  fEvMult->Fill(aEvent->NumberOfTracks());
  fNormEvMult->Fill(aEvent->UncorrectedNumberOfPrimaries());
  if(!freadMC){
    fSPDMult->Fill(aEvent->SPDMultiplicity());
  }
  fMultSumPt->Fill(aEvent->UncorrectedNumberOfPrimaries(), aEvent->ZDCEMEnergy());

  if(faddhists)
    {
      fEstimateITSTPC->Fill(aEvent->MultiplicityEstimateITSTPC());
      fEstimateTracklets->Fill(aEvent->MultiplicityEstimateTracklets());
      fEstimateITSPure->Fill(aEvent->MultiplicityEstimateITSPure());
      fEst1Est2->Fill(aEvent->MultiplicityEstimateITSTPC(),aEvent->MultiplicityEstimateTracklets());
      fEst1Est3->Fill(aEvent->MultiplicityEstimateITSTPC(),aEvent->MultiplicityEstimateITSPure());
      fEst2Est3->Fill(aEvent->MultiplicityEstimateTracklets(),aEvent->MultiplicityEstimateITSPure());
      fEst1Norm->Fill(aEvent->MultiplicityEstimateITSTPC(),aEvent->UncorrectedNumberOfPrimaries());
      fEst2Norm->Fill(aEvent->MultiplicityEstimateTracklets(),aEvent->UncorrectedNumberOfPrimaries());
      fEst3Norm->Fill(aEvent->MultiplicityEstimateITSPure(),aEvent->UncorrectedNumberOfPrimaries());
    }

}

void AliFemtoCutMonitorEventMult::Write()
{
  // Write out the relevant histograms
  fEvMult->Write();
  fNormEvMult->Write();
  if(!freadMC){
    fSPDMult->Write();
  }
  fMultSumPt->Write();

  if(faddhists)
    {
      fEstimateITSTPC->Write();
      fEstimateTracklets->Write();
      fEstimateITSPure->Write();
      fEst1Est2->Write();
      fEst1Est3->Write();
      fEst2Est3->Write();
      fEst1Norm->Write();
      fEst2Norm->Write();
      fEst3Norm->Write();
    }

}

TList *AliFemtoCutMonitorEventMult::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fEvMult);
  tOutputList->Add(fNormEvMult);
  tOutputList->Add(fSPDMult);
  tOutputList->Add(fMultSumPt);

  if(faddhists)
    {
      tOutputList->Add(fEstimateITSTPC);
      tOutputList->Add(fEstimateTracklets);
      tOutputList->Add(fEstimateITSPure);
      tOutputList->Add(fEst1Est2);
      tOutputList->Add(fEst1Est3);
      tOutputList->Add(fEst2Est3);
      tOutputList->Add(fEst1Norm);
      tOutputList->Add(fEst2Norm);
      tOutputList->Add(fEst3Norm);
    }
  
  return tOutputList;
}

void AliFemtoCutMonitorEventMult::SetReadMC(Bool_t mc)
{
  freadMC=mc;
}

void AliFemtoCutMonitorEventMult::AdditionalMultHistsOn(Bool_t addhists)
{
  faddhists=addhists;
}

