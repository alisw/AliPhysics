// $Id$
//
// Scale task.
//
// Author: R.Reed, M.Connors

#include "AliAnalysisTaskScale.h"

#include <TClonesArray.h>
#include <TF1.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include "AliEMCALGeometry.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

ClassImp(AliAnalysisTaskScale)

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskScale", kTRUE), 
  fScaleFunction(0),
  fEmcalArea(1),
  fTpcArea(1),
  fHistPtTPCvsCent(0), 
  fHistPtEMCALvsCent(0), 
  fHistEtvsCent(0),  
  fHistScalevsCent(0),  
  fHistDeltaScalevsCent(0), 
  fHistScaleEmcalvsCent(0),      
  fHistScale2EmcalvsCent(0),
  fHistDeltaScale2EmcalvsCent(0),
  fHistScale3EmcalvsCent(0),
  fHistDeltaScale3EmcalvsCent(0),     
  fHistChScalevsCent(0),          
  fHistChScale2EmcalvsCent(0),
  fHistChScale3EmcalvsCent(0),      
  fHistPtTPCvsNtrack(0), 
  fHistPtEMCALvsNtrack(0), 
  fHistEtvsNtrack(0),  
  fHistScalevsNtrack(0),  
  fHistDeltaScalevsNtrack(0),
  fHistScaleEmcalvsNtrack(0),      
  fHistScale2EmcalvsNtrack(0),
  fHistScale3EmcalvsNtrack(0),     
  fHistChScalevsNtrack(0),          
  fHistChScale2EmcalvsNtrack(0),
  fHistChScale3EmcalvsNtrack(0),   
  fHistTrackPtvsCent(0), 
  fHistClusterPtvsCent(0),
  fHistTrackEtaPhi(0), 
  fHistClusterEtaPhi(0),
  fHistScalevsScale2Emcal(0),
  fHistScalevsScale3Emcal(0),      
  fHistScalevsScaleEmcal(0),       
  fHistScaleEmcalvsScale2Emcal(0),
  fHistScaleEmcalvsScale3Emcal(0),
  fHistScaleShift1EmcalvsCent(0),
  fHistScaleShift2EmcalvsCent(0),
  fHistScaleShiftMeanEmcalvsCent(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskScale::AliAnalysisTaskScale(const char *name) :
  AliAnalysisTaskEmcal(name, kTRUE), 
  fScaleFunction(0),
  fEmcalArea(1),
  fTpcArea(1),
  fHistPtTPCvsCent(0), 
  fHistPtEMCALvsCent(0), 
  fHistEtvsCent(0),  
  fHistScalevsCent(0),  
  fHistDeltaScalevsCent(0), 
  fHistScaleEmcalvsCent(0),      
  fHistScale2EmcalvsCent(0),
  fHistDeltaScale2EmcalvsCent(0),
  fHistScale3EmcalvsCent(0),
  fHistDeltaScale3EmcalvsCent(0),  
  fHistChScalevsCent(0),          
  fHistChScale2EmcalvsCent(0),
  fHistChScale3EmcalvsCent(0),      
  fHistPtTPCvsNtrack(0), 
  fHistPtEMCALvsNtrack(0), 
  fHistEtvsNtrack(0),  
  fHistScalevsNtrack(0),  
  fHistDeltaScalevsNtrack(0),
  fHistScaleEmcalvsNtrack(0),      
  fHistScale2EmcalvsNtrack(0),
  fHistScale3EmcalvsNtrack(0),          
  fHistChScalevsNtrack(0),          
  fHistChScale2EmcalvsNtrack(0),
  fHistChScale3EmcalvsNtrack(0),      
  fHistTrackPtvsCent(0), 
  fHistClusterPtvsCent(0),
  fHistTrackEtaPhi(0), 
  fHistClusterEtaPhi(0),
  fHistScalevsScale2Emcal(0),
  fHistScalevsScale3Emcal(0),            
  fHistScalevsScaleEmcal(0),       
  fHistScaleEmcalvsScale2Emcal(0),
  fHistScaleEmcalvsScale3Emcal(0),
  fHistScaleShift1EmcalvsCent(0),
  fHistScaleShift2EmcalvsCent(0),
  fHistScaleShiftMeanEmcalvsCent(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskScale::UserCreateOutputObjects()
{
  // Create my user objects.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  //Get track and particle container
  fTracksCont       = GetParticleContainer(0);
  fCaloClustersCont = GetClusterContainer(0);
  if(fTracksCont)       fTracksCont->SetClassName("AliVTrack");
  if(fCaloClustersCont) fCaloClustersCont->SetClassName("AliVCluster");

  //Create histos
  fHistPtTPCvsCent = new TH2F("fHistPtTPCvsCent", "fHistPtTPCvsCent", 101, -1, 100, 750, 0, 1500);
  fHistPtTPCvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistPtTPCvsCent->GetYaxis()->SetTitle("#sum p_{T,track}^{TPC} GeV/c");
  fHistPtTPCvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistPtTPCvsCent);

  fHistPtEMCALvsCent = new TH2F("fHistPtEMCALvsCent", "fHistPtEMCALvsCent", 101, -1, 100, 250, 0, 500);
  fHistPtEMCALvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistPtEMCALvsCent->GetYaxis()->SetTitle("#sum p_{T,track}^{EMCal} GeV/c");
  fHistPtEMCALvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistPtEMCALvsCent);

  fHistEtvsCent = new TH2F("fHistEtvsCent", "fHistEtvsCent", 101, -1, 100, 250, 0, 500);
  fHistEtvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistEtvsCent->GetYaxis()->SetTitle("#sum E_{T,cluster} GeV");
  fHistEtvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistEtvsCent);

  fHistScalevsCent = new TH2F("fHistScalevsCent", "fHistScalevsCent", 101, -1, 100, 500, 0, 5);
  fHistScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScalevsCent->GetYaxis()->SetTitle("s_{TPC} = (#sum E_{T,cluster} + #sum p_{T,track}^{TPC}) / #sum p_{T,track}^{TPC}");
  fHistScalevsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScalevsCent);

  fHistDeltaScalevsCent = new TH2F("fHistDeltaScalevsCent", "fHistDeltaScalevsCent", 101, -1, 100, 500, -2.5, 2.5);
  fHistDeltaScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistDeltaScalevsCent->GetYaxis()->SetTitle("s_{TPC}-s^{old}");
  fHistDeltaScalevsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaScalevsCent);

  fHistScaleEmcalvsCent= new TH2F("fHistScaleEmcalvsCent", "fHistScaleEmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistScaleEmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScaleEmcalvsCent->GetYaxis()->SetTitle("s_{EMC}");
  fHistScaleEmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleEmcalvsCent);

  fHistScale2EmcalvsCent = new TH2F("fHistScale2EmcalvsCent", "fHistScale2EmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistScale2EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScale2EmcalvsCent->GetYaxis()->SetTitle("s_{2 #times EMC}");
  fHistScale2EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScale2EmcalvsCent);

  fHistDeltaScale2EmcalvsCent = new TH2F("fHistDeltaScale2EmcalvsCent", "fHistDeltaScale2EmcalvsCent", 101, -1, 100, 500, -2.5, 2.5);
  fHistDeltaScale2EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistDeltaScale2EmcalvsCent->GetYaxis()->SetTitle("s_{2 #times EMC}-s^{old}");
  fHistDeltaScale2EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaScale2EmcalvsCent);

  fHistScale3EmcalvsCent = new TH2F("fHistScale3EmcalvsCent", "fHistScale3EmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistScale3EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScale3EmcalvsCent->GetYaxis()->SetTitle("s_{3 #times EMC}");
  fHistScale3EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScale3EmcalvsCent);

  fHistScaleShift1EmcalvsCent = new TH2F("fHistScaleShift1EmcalvsCent", "fHistScaleShift1EmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistScaleShift1EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScaleShift1EmcalvsCent->GetYaxis()->SetTitle("s_{EMC}^{1xShift}");
  fHistScaleShift1EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleShift1EmcalvsCent);

  fHistScaleShift2EmcalvsCent = new TH2F("fHistScaleShift2EmcalvsCent", "fHistScaleShift2EmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistScaleShift2EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScaleShift2EmcalvsCent->GetYaxis()->SetTitle("s_{EMC}^{2xShift}");
  fHistScaleShift2EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleShift2EmcalvsCent);

  fHistScaleShiftMeanEmcalvsCent = new TH2F("fHistScaleShiftMeanEmcalvsCent", "fHistScaleShiftMeanEmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistScaleShiftMeanEmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistScaleShiftMeanEmcalvsCent->GetYaxis()->SetTitle("s_{EMC}^{<Shift>}");
  fHistScaleShiftMeanEmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleShiftMeanEmcalvsCent);

  fHistDeltaScale3EmcalvsCent = new TH2F("fHistDeltaScale3EmcalvsCent", "fHistDeltaScale3EmcalvsCent", 101, -1, 100, 500, -2.5, 2.5);
  fHistDeltaScale3EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistDeltaScale3EmcalvsCent->GetYaxis()->SetTitle("s_{3 #times EMC}-s^{old}");
  fHistDeltaScale3EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaScale3EmcalvsCent);

  fHistChScalevsCent = new TH2F("fHistChScalevsCent", "fHistChScalevsCent", 101, -1, 100, 500, 0, 5);
  fHistChScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistChScalevsCent->GetYaxis()->SetTitle("s_{TPC}^{ch}");
  fHistChScalevsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistChScalevsCent);

  fHistChScale2EmcalvsCent = new TH2F("fHistChScale2EmcalvsCent", "fHistChScale2EmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistChScale2EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistChScale2EmcalvsCent->GetYaxis()->SetTitle("s_{2 #times EMC}^{ch}");
  fHistChScale2EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistChScale2EmcalvsCent);

  fHistChScale3EmcalvsCent = new TH2F("fHistChScale3EmcalvsCent", "fHistChScale3EmcalvsCent", 101, -1, 100, 500, 0, 5);
  fHistChScale3EmcalvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistChScale3EmcalvsCent->GetYaxis()->SetTitle("s_{3 #times EMC}^{ch}");
  fHistChScale3EmcalvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistChScale3EmcalvsCent);

  fHistPtTPCvsNtrack = new TH2F("fHistPtTPCvsNtrack", "fHistPtTPCvsNtrack", 800, 0, 4000,  750, 0, 1500);
  fHistPtTPCvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistPtTPCvsNtrack->GetYaxis()->SetTitle("#sum p_{T,track}^{TPC}");
  fHistPtTPCvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistPtTPCvsNtrack);

  fHistPtEMCALvsNtrack = new TH2F("fHistPtEMCALvsNtrack", "fHistPtEMCALvsNtrack", 800, 0, 4000,  500, 0, 1000);
  fHistPtEMCALvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistPtEMCALvsNtrack->GetYaxis()->SetTitle("#sum p_{T,track}^{EMCal}");
  fHistPtEMCALvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistPtEMCALvsNtrack);

  fHistEtvsNtrack = new TH2F("fHistEtvsNtrack", "fHistEtvsNtrack", 800,  0, 4000, 500, 0, 1000);
  fHistEtvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistEtvsNtrack->GetYaxis()->SetTitle("#sum E_{T,cluster}");
  fHistEtvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistEtvsNtrack);

  fHistScalevsNtrack = new TH2F("fHistScalevsNtrack", "fHistScalevsNtrack", 800, 0, 4000,  500, 0, 5);
  fHistScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistScalevsNtrack->GetYaxis()->SetTitle("s_{TPC}");
  fHistScalevsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScalevsNtrack);

  fHistDeltaScalevsNtrack = new TH2F("fHistDeltaScalevsNtrack", "fHistDeltaScalevsNtrack", 800, 0, 4000, 500, -2.5, 2.5);
  fHistDeltaScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistDeltaScalevsNtrack->GetYaxis()->SetTitle("s_{TPC}-s^{old}");
  fHistDeltaScalevsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaScalevsNtrack);

  fHistScaleEmcalvsNtrack = new TH2F("fHistScaleEmcalvsNtrack", "fHistScaleEmcalvsNtrack", 800, 0, 4000, 500, 0, 5);
  fHistScaleEmcalvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistScaleEmcalvsNtrack->GetYaxis()->SetTitle("s_{EMC}");
  fHistScaleEmcalvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleEmcalvsNtrack);

  fHistScale2EmcalvsNtrack = new TH2F("fHistScale2EmcalvsNtrack","fHistScale2EmcalvsNtrack", 800, 0, 4000, 500, 0, 5);
  fHistScale2EmcalvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistScale2EmcalvsNtrack->GetYaxis()->SetTitle("s_{2 #times EMC}");
  fHistScale2EmcalvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScale2EmcalvsNtrack);

  fHistScale3EmcalvsNtrack = new TH2F("fHistScale3EmcalvsNtrack","fHistScale3EmcalvsNtrack", 800, 0, 4000, 500, 0, 5);
  fHistScale3EmcalvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistScale3EmcalvsNtrack->GetYaxis()->SetTitle("s_{3 #times EMC}");
  fHistScale3EmcalvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScale3EmcalvsNtrack);

  fHistChScalevsNtrack = new TH2F("fHistChScalevsNtrack", "fHistChScalevsNtrack", 800, 0, 4000, 500, 0, 5);
  fHistChScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistChScalevsNtrack->GetYaxis()->SetTitle("s_{TPC}^{ch}");
  fHistChScalevsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistChScalevsNtrack);

  fHistChScale2EmcalvsNtrack = new TH2F("fHistChScale2EmcalvsNtrack", "fHistChScale2EmcalvsNtrack", 800,  0, 4000, 500, 0, 5);
  fHistChScale2EmcalvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistChScale2EmcalvsNtrack->GetYaxis()->SetTitle("s_{2 #times EMC}^{ch}");
  fHistChScale2EmcalvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistChScale2EmcalvsNtrack);

  fHistChScale3EmcalvsNtrack = new TH2F("fHistChScale3EmcalvsNtrack", "fHistChScale3EmcalvsNtrack", 800,  0, 4000, 500, 0, 5);
  fHistChScale3EmcalvsNtrack->GetXaxis()->SetTitle("No. of tracks");
  fHistChScale3EmcalvsNtrack->GetYaxis()->SetTitle("s_{3 #times EMC}^{ch}");
  fHistChScale3EmcalvsNtrack->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistChScale3EmcalvsNtrack);

  fHistTrackPtvsCent = new TH2F("fHistTrackPtvsCent", "fHistTrackPtvsCent", 101, -1, 100, 500, 0, 100);
  fHistTrackPtvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistTrackPtvsCent->GetYaxis()->SetTitle("p_{T,track} GeV/c");
  fHistTrackPtvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistTrackPtvsCent);

  fHistClusterPtvsCent = new TH2F("fHistClusterPtvsCent", "fHistClusterPtvsCent", 101, -1, 100, 500, 0, 100);
  fHistClusterPtvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistClusterPtvsCent->GetYaxis()->SetTitle("E_{T,cluster} GeV");
  fHistClusterPtvsCent->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistClusterPtvsCent);

  fHistTrackEtaPhi = new TH2F("fHistTrackEtaPhi", "fHistTrackEtaPhi", 100, -1.0, 1.0, 101, 0, 2.02*TMath::Pi());
  fHistTrackEtaPhi->GetXaxis()->SetTitle("#eta");
  fHistTrackEtaPhi->GetYaxis()->SetTitle("#phi");
  fHistTrackEtaPhi->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistTrackEtaPhi);

  fHistClusterEtaPhi = new TH2F("fHistClusterEtaPhi", "fHistClusterEtaPhi", 100, -1.0, 1.0, 101, 0, 2.02*TMath::Pi());
  fHistClusterEtaPhi->GetXaxis()->SetTitle("#eta");
  fHistClusterEtaPhi->GetYaxis()->SetTitle("#phi");
  fHistClusterEtaPhi->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistClusterEtaPhi);

  fHistScalevsScale3Emcal = new TH2F("fHistScalevsScale3Emcal", "fHistScalevsScale3Emcal",500, 0, 5, 500, 0, 5);
  fHistScalevsScale3Emcal->GetXaxis()->SetTitle("s_{TPC}");
  fHistScalevsScale3Emcal->GetYaxis()->SetTitle("s_{3 #times EMC}");
  fHistScalevsScale3Emcal->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScalevsScale3Emcal);

  fHistScalevsScale2Emcal = new TH2F("fHistScalevsScale2Emcal", "fHistScalevsScale2Emcal",500, 0, 5, 500, 0, 5);
  fHistScalevsScale2Emcal->GetXaxis()->SetTitle("s_{TPC}");
  fHistScalevsScale2Emcal->GetYaxis()->SetTitle("s_{2 #times EMC}");
  fHistScalevsScale2Emcal->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScalevsScale2Emcal);

  fHistScalevsScaleEmcal = new TH2F("fHistScalevsScaleEmcal", "fHistScalevsScaleEmcal", 500, 0, 5, 500, 0, 5);
  fHistScalevsScaleEmcal->GetXaxis()->SetTitle("s_{TPC}");
  fHistScalevsScaleEmcal->GetYaxis()->SetTitle("s_{EMC}");
  fHistScalevsScaleEmcal->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScalevsScaleEmcal);

  fHistScaleEmcalvsScale3Emcal = new TH2F("fHistScaleEmcalvsScale3Emcal", "fHistScaleEmcalvsScale3Emcal", 500, 0, 5, 500, 0, 5);
  fHistScaleEmcalvsScale3Emcal->GetXaxis()->SetTitle("s_{EMC}");
  fHistScaleEmcalvsScale3Emcal->GetYaxis()->SetTitle("s_{3 #times EMC}");
  fHistScaleEmcalvsScale3Emcal->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleEmcalvsScale3Emcal);

  fHistScaleEmcalvsScale2Emcal = new TH2F("fHistScaleEmcalvsScale2Emcal", "fHistScaleEmcalvsScale2Emcal", 500, 0, 5, 500, 0, 5);
  fHistScaleEmcalvsScale2Emcal->GetXaxis()->SetTitle("s_{EMC}");
  fHistScaleEmcalvsScale2Emcal->GetYaxis()->SetTitle("s_{2 #times EMC}");
  fHistScaleEmcalvsScale2Emcal->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistScaleEmcalvsScale2Emcal);

  PostData(1, fOutput);
}

//________________________________________________________________________
Double_t AliAnalysisTaskScale::GetScaleFactor(Double_t cent)
{
  // Get scale function.

  Double_t scale = -1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskScale::FillHistograms() 
{
  // Execute on each event.

  const Double_t EmcalMinEta = fGeom->GetArm1EtaMin();
  const Double_t EmcalMaxEta = fGeom->GetArm1EtaMax();

  Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
  Double_t EmcalMaxPhi = fGeom->GetEMCALPhiMax() * TMath::DegToRad();
  if(InputEvent()->GetRunNumber()>=177295 && InputEvent()->GetRunNumber()<=197470) { //small SM masked in 2012 and 2013
    EmcalMinPhi = 1.405;
    EmcalMaxPhi = 3.135;
  }
  const Double_t EmcalWidth = (EmcalMaxPhi-EmcalMinPhi)/2.0;

  Double_t ptTPC   = 0;
  Double_t ptEMCAL = 0;
  Double_t ptEMCAL2 = 0;
  Double_t ptEMCAL3 = 0;

  Double_t ptEMCALshift1 = 0;
  Double_t ptEMCALshift2 = 0;
  Double_t ptEMCALshiftMean = 0;

  const Int_t Ntracks = fTracksCont->GetNAcceptedParticles();
  if (fTracksCont) {
    fTracksCont->ResetCurrentID();
    AliVTrack *track = NULL;
    while((track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))) {
      if(track->Eta() < EmcalMinEta || track->Eta() > EmcalMaxEta) continue;
      fHistTrackPtvsCent->Fill(fCent,track->Pt());
      fHistTrackEtaPhi->Fill(track->Eta(),track->Phi());
      ptTPC += track->Pt();

      // Shift track acceptance by pi/4 and take mean
      for(int a = 0; a < 7; a++){
        if((track->Phi() < (EmcalMaxPhi+(a*TMath::PiOver4())) && track->Phi() > (EmcalMinPhi+(a*TMath::PiOver4())))) ptEMCALshiftMean += track->Pt();
      }
      ptEMCALshiftMean /= 7;

      // Shift track acceptance by EMCal width
      if ((track->Phi() < (EmcalMaxPhi+(2*EmcalWidth))) && (track->Phi() > (EmcalMinPhi+(2*EmcalWidth)))) ptEMCALshift1 += track->Pt(); // shift by width of emcal
      if ((track->Phi() < (EmcalMaxPhi+(4*EmcalWidth))) && (track->Phi() > (EmcalMinPhi+(4*EmcalWidth)))) ptEMCALshift2 += track->Pt(); // shift by 2x width of emcal


      if ((track->Phi() > (EmcalMaxPhi+(2.0*EmcalWidth))) || (track->Phi() < (EmcalMinPhi-(2.0*EmcalWidth)))) continue;
      ptEMCAL3 += track->Pt();
      if ((track->Phi() > (EmcalMaxPhi+EmcalWidth)) || (track->Phi() < (EmcalMinPhi-EmcalWidth))) continue;
      ptEMCAL2 += track->Pt();
      if ((track->Phi() > EmcalMaxPhi) || (track->Phi() < EmcalMinPhi)) continue;
      ptEMCAL += track->Pt();
    }
  }

  if (ptTPC == 0) 
    return kFALSE;

  Double_t Et = 0;  
  if (fCaloClustersCont) {
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *c = NULL;
    while((c = fCaloClustersCont->GetNextAcceptCluster())) {
      TLorentzVector nPart;
      c->GetMomentum(nPart, fVertex);

      fHistClusterPtvsCent->Fill(fCent, nPart.Pt());
      fHistClusterEtaPhi->Fill(nPart.Eta(), nPart.Phi());

      Et += nPart.Pt();
    }
  }

  if (Et == 0) 
    return kFALSE;
 
  Double_t scalecalc         = -1;
  if (ptEMCAL > 0 && Et > 0 && ptTPC > 0)
    scalecalc         =  ((Et + ptEMCAL) / fEmcalArea) * (fTpcArea / ptTPC);
  const Double_t scale             = GetScaleFactor(fCent);
  Double_t scalecalcemcal          = -1;
  if (ptEMCAL > 0)
    scalecalcemcal                 = (Et+ptEMCAL)/ptEMCAL;
  Double_t scalecalcemcal2         = -1;
  Double_t Chscalecalcemcal2       = -1;
  if (ptEMCAL2 > 0){
    scalecalcemcal2                = 2*(Et+ptEMCAL)/ptEMCAL2;
    Chscalecalcemcal2              = 2*ptEMCAL/ptEMCAL2;}
  const Double_t Chscalecalcemcal  = ((ptEMCAL) / fEmcalArea) * (fTpcArea / ptTPC);
  Double_t scalecalcemcal3         = -1;
  Double_t Chscalecalcemcal3       = -1;
  if (ptEMCAL3 > 0){
    scalecalcemcal3                = 3*(Et+ptEMCAL)/ptEMCAL3;
    Chscalecalcemcal3              = 3*ptEMCAL/ptEMCAL3;}

  // Calculations for shifted EMCal acceptance
  Double_t scalecalcemcalshift1    = -1;
  if (ptEMCALshift1 > 0)
    scalecalcemcalshift1           = (Et+ptEMCALshift1)/ptEMCALshift1;
  Double_t scalecalcemcalshift2    = -1;
  if (ptEMCALshift2 > 0)
    scalecalcemcalshift2           = (Et+ptEMCALshift2)/ptEMCALshift2;
  Double_t scalecalcemcalshiftmean    = -1;
  if (ptEMCALshiftMean > 0)
    scalecalcemcalshiftmean           = (Et+ptEMCALshiftMean)/ptEMCALshiftMean;

  fHistScaleEmcalvsCent->Fill(fCent,scalecalcemcal);      
  fHistScale2EmcalvsCent->Fill(fCent,scalecalcemcal2);    
  fHistScale3EmcalvsCent->Fill(fCent,scalecalcemcal3);
  fHistScaleShift1EmcalvsCent->Fill(fCent,scalecalcemcalshift1);
  fHistScaleShift2EmcalvsCent->Fill(fCent,scalecalcemcalshift2);  
  fHistScaleShiftMeanEmcalvsCent->Fill(fCent,scalecalcemcalshiftmean);    
  fHistChScalevsCent->Fill(fCent,Chscalecalcemcal);    
  fHistChScale2EmcalvsCent->Fill(fCent,Chscalecalcemcal2);   
  fHistChScale3EmcalvsCent->Fill(fCent,Chscalecalcemcal3);   
  fHistScaleEmcalvsNtrack->Fill(Ntracks,scalecalcemcal);      
  fHistScale2EmcalvsNtrack->Fill(Ntracks,scalecalcemcal2);
  fHistScale3EmcalvsNtrack->Fill(Ntracks,scalecalcemcal3);        
  fHistChScalevsNtrack->Fill(Ntracks,Chscalecalcemcal);    
  fHistChScale2EmcalvsNtrack->Fill(Ntracks,Chscalecalcemcal2); 
  fHistChScale3EmcalvsNtrack->Fill(Ntracks,Chscalecalcemcal3);     
  fHistPtTPCvsCent->Fill(fCent, ptTPC);
  fHistPtEMCALvsCent->Fill(fCent, ptEMCAL);
  fHistEtvsCent->Fill(fCent, Et);
  fHistScalevsCent->Fill(fCent, scalecalc);
  fHistDeltaScalevsCent->Fill(fCent, scalecalc - scale);
  fHistDeltaScale2EmcalvsCent->Fill(fCent, scalecalcemcal2 - scale);
  fHistDeltaScale3EmcalvsCent->Fill(fCent, scalecalcemcal3 - scale);
  fHistPtTPCvsNtrack->Fill(Ntracks, ptTPC);
  fHistPtEMCALvsNtrack->Fill(Ntracks, ptEMCAL);
  fHistEtvsNtrack->Fill(Ntracks, Et);
  fHistScalevsNtrack->Fill(Ntracks, scalecalc);
  fHistDeltaScalevsNtrack->Fill(Ntracks, scalecalc - scale);
  fHistScalevsScale2Emcal->Fill(scalecalc,scalecalcemcal2);  
  fHistScalevsScale3Emcal->Fill(scalecalc,scalecalcemcal3);        
  fHistScalevsScaleEmcal->Fill(scalecalc,scalecalcemcal); 
  fHistScaleEmcalvsScale2Emcal->Fill(scalecalcemcal,scalecalcemcal2);
  fHistScaleEmcalvsScale3Emcal->Fill(scalecalcemcal,scalecalcemcal3);


  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskScale::ExecOnce() 
{
  AliAnalysisTaskEmcal::ExecOnce();

  const Double_t EmcalMinEta = fGeom->GetArm1EtaMin();
  const Double_t EmcalMaxEta = fGeom->GetArm1EtaMax();

  Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
  Double_t EmcalMaxPhi = fGeom->GetEMCALPhiMax() * TMath::DegToRad();

  if(InputEvent()->GetRunNumber()>=177295 && InputEvent()->GetRunNumber()<=197470) { //small SM masked in 2012 and 2013
    EmcalMinPhi = 1.405;
    EmcalMaxPhi = 3.135;
  }

  fEmcalArea  = (EmcalMaxPhi - EmcalMinPhi) * (EmcalMinEta - EmcalMaxEta);

  AliParticleContainer *partCont = GetParticleContainer(0);
  if (!partCont) {
    AliError(Form("%s: No particle container found! Assuming tpc area = 1...",GetName()));
    fTpcArea = 1;
    return;
  }

  Float_t TpcMaxPhi = partCont->GetParticlePhiMax();
  Float_t TpcMinPhi = partCont->GetParticlePhiMin();
  
  if (TpcMaxPhi > TMath::Pi()*2) TpcMaxPhi = TMath::Pi()*2;
  if (TpcMinPhi < 0) TpcMinPhi = 0;

  fTpcArea = (TpcMaxPhi - TpcMinPhi) * (EmcalMinEta - EmcalMaxEta);

  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;
}
