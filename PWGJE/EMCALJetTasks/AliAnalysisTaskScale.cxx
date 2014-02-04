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
#include "AliParticleContainer.h"

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
  fHistChScalevsCent(0),          
  fHistChScale2EmcalvsCent(0),   
  fHistPtTPCvsNtrack(0), 
  fHistPtEMCALvsNtrack(0), 
  fHistEtvsNtrack(0),  
  fHistScalevsNtrack(0),  
  fHistDeltaScalevsNtrack(0),
  fHistScaleEmcalvsNtrack(0),      
  fHistScale2EmcalvsNtrack(0),     
  fHistChScalevsNtrack(0),          
  fHistChScale2EmcalvsNtrack(0),   
  fHistTrackPtvsCent(0), 
  fHistClusterPtvsCent(0),
  fHistTrackEtaPhi(0), 
  fHistClusterEtaPhi(0),
  fHistScalevsScale2Emcal(0),      
  fHistScalevsScaleEmcal(0),       
  fHistScaleEmcalvsScale2Emcal(0)
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
  fHistChScalevsCent(0),          
  fHistChScale2EmcalvsCent(0),   
  fHistPtTPCvsNtrack(0), 
  fHistPtEMCALvsNtrack(0), 
  fHistEtvsNtrack(0),  
  fHistScalevsNtrack(0),  
  fHistDeltaScalevsNtrack(0),
  fHistScaleEmcalvsNtrack(0),      
  fHistScale2EmcalvsNtrack(0),     
  fHistChScalevsNtrack(0),          
  fHistChScale2EmcalvsNtrack(0),   
  fHistTrackPtvsCent(0), 
  fHistClusterPtvsCent(0),
  fHistTrackEtaPhi(0), 
  fHistClusterEtaPhi(0),
  fHistScalevsScale2Emcal(0),      
  fHistScalevsScaleEmcal(0),       
  fHistScaleEmcalvsScale2Emcal(0)
{
  // Constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskScale::UserCreateOutputObjects()
{
  // Create my user objects.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fHistPtTPCvsCent             = new TH2F("PtTPCvsCent","rho vs cent",            101, -1, 100,   500,   0, 1000);
  fHistPtEMCALvsCent           = new TH2F("PtEMCALvsCent","rho vs cent",          101, -1, 100,   500,   0, 1000);
  fHistEtvsCent                = new TH2F("EtvsCent","rho vs cent",               101, -1, 100,   500,   0, 1000);
  fHistScalevsCent             = new TH2F("ScalevsCent","rho vs cent",            101, -1, 100,   500,   0, 5);
  fHistDeltaScalevsCent        = new TH2F("DeltaScalevsCent","rho vs cent",       101, -1, 100,   500,  -2.5, 2.5);
  fHistScaleEmcalvsCent        = new TH2F("ScaleEmcalvsCent","",                  101, -1, 100,   500,   0, 5);
  fHistScale2EmcalvsCent       = new TH2F("Scale2EmcalvsCent","",                 101, -1, 100,   500,   0, 5);
  fHistChScalevsCent           = new TH2F("ChScalevsCent","",                     101, -1, 100,   500,   0, 5);
  fHistChScale2EmcalvsCent     = new TH2F("ChScale2EmcalvsCent","",               101, -1, 100,   500,   0, 5);
  fHistPtTPCvsNtrack           = new TH2F("PtTPCvsNtrack","rho vs cent",          500,  0, 2500,  500,   0, 1000);
  fHistPtEMCALvsNtrack         = new TH2F("PtEMCALvsNtrack","rho vs cent",        500,  0, 2500,  500,   0, 1000);
  fHistEtvsNtrack              = new TH2F("EtvsNtrack","rho vs cent",             500,  0, 2500,  500,   0, 1000);
  fHistScalevsNtrack           = new TH2F("ScalevsNtrack","rho vs cent",          500,  0, 2500,  500,   0, 5);
  fHistDeltaScalevsNtrack      = new TH2F("DeltaScalevsNtrack","rho vs cent",     500,  0, 2500,  500,  -2.5, 2.5);
  fHistScaleEmcalvsNtrack      = new TH2F("ScaleEmcalvsNtrack","",                500,  0, 2500,  500,   0, 5);
  fHistScale2EmcalvsNtrack     = new TH2F("Scale2EmcalvsNtrack","",               500,  0, 2500,  500,   0, 5);
  fHistChScalevsNtrack         = new TH2F("ChScalevsNtrack","",                   500,  0, 2500,  500,   0, 5);
  fHistChScale2EmcalvsNtrack   = new TH2F("ChScale2EmcalvsNtrack","",             500,  0, 2500,  500,   0, 5);
  fHistTrackPtvsCent           = new TH2F("TrackPtvsCent","Track pt vs cent",     101, -1, 100,   500,   0, 100);
  fHistClusterPtvsCent         = new TH2F("ClusterPtvsCent","Cluster pt vs cent", 101, -1, 100,   500,   0, 100);
  fHistTrackEtaPhi             = new TH2F("TrackEtaPhi","Track eta phi",          100, -1.0, 1.0, 101,   0, 2.02*TMath::Pi());
  fHistClusterEtaPhi           = new TH2F("ClusterEtaPhi","Cluster eta phi",      100, -1.0, 1.0, 101,   0, 2.02*TMath::Pi());
  fHistScalevsScale2Emcal      = new TH2F("ScalevsScale2Emcal","",                500,  0, 5,     500,   0, 5);
  fHistScalevsScaleEmcal       = new TH2F("ScalevsScaleEmcal","",                 500,  0, 5,     500,   0, 5);
  fHistScaleEmcalvsScale2Emcal = new TH2F("ScaleEmcalvsScale2Emcal","",           500,  0, 5,     500,   0, 5);

  fOutput->Add(fHistPtTPCvsCent);
  fOutput->Add(fHistPtEMCALvsCent);
  fOutput->Add(fHistEtvsCent);
  fOutput->Add(fHistScalevsCent);
  fOutput->Add(fHistDeltaScalevsCent);
  fOutput->Add(fHistScaleEmcalvsCent);      
  fOutput->Add(fHistScale2EmcalvsCent);     
  fOutput->Add(fHistChScalevsCent);    
  fOutput->Add(fHistChScale2EmcalvsCent);   
  fOutput->Add(fHistPtTPCvsNtrack);
  fOutput->Add(fHistPtEMCALvsNtrack);
  fOutput->Add(fHistEtvsNtrack);
  fOutput->Add(fHistScalevsNtrack);
  fOutput->Add(fHistDeltaScalevsNtrack);
  fOutput->Add(fHistScaleEmcalvsNtrack);      
  fOutput->Add(fHistScale2EmcalvsNtrack);     
  fOutput->Add(fHistChScalevsNtrack);    
  fOutput->Add(fHistChScale2EmcalvsNtrack);   
  fOutput->Add(fHistTrackPtvsCent);
  fOutput->Add(fHistClusterPtvsCent);
  fOutput->Add(fHistTrackEtaPhi);
  fOutput->Add(fHistClusterEtaPhi);
  fOutput->Add(fHistScalevsScale2Emcal);      
  fOutput->Add(fHistScalevsScaleEmcal);       
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

  const Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
  const Double_t EmcalMaxPhi = fGeom->GetArm1PhiMax() * TMath::DegToRad();
  const Double_t EmcalWidth = (EmcalMaxPhi-EmcalMinPhi)/2.0;

  Double_t ptTPC   = 0;
  Double_t ptEMCAL = 0;
  Double_t ptEMCAL2 = 0;

  const Int_t Ntracks = fTracks->GetEntries();
  for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
    AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(iTracks));

    if (!track)
      continue;

    if (!AcceptTrack(track))
      continue;

    if (TMath::Abs(track->Eta()) > 0.7)   // only accept tracks in the EMCal eta range
      continue;

    fHistTrackPtvsCent->Fill(fCent,track->Pt());
    fHistTrackEtaPhi->Fill(track->Eta(),track->Phi());
    ptTPC += track->Pt();
    if ((track->Phi() > (EmcalMaxPhi+EmcalWidth)) || (track->Phi() < (EmcalMinPhi-EmcalWidth))) continue;
    ptEMCAL2 += track->Pt();
    if ((track->Phi() > EmcalMaxPhi) || (track->Phi() < EmcalMinPhi)) continue;
    ptEMCAL += track->Pt();
  }

  if (ptTPC == 0) 
    return kFALSE;
  
  Double_t Et = 0;
  const Int_t Nclus = fCaloClusters->GetEntries();
  for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
    AliVCluster *c = static_cast<AliVCluster*>(fCaloClusters->At(iClus));
    if (!c)
      continue;

    if (!AcceptCluster(c))
      continue;

    TLorentzVector nPart;
    c->GetMomentum(nPart, fVertex);

    fHistClusterPtvsCent->Fill(fCent, nPart.Pt());
    fHistClusterEtaPhi->Fill(nPart.Eta(), nPart.Phi());

    Et += nPart.Pt();
  }
 
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

  fHistScaleEmcalvsCent->Fill(fCent,scalecalcemcal);      
  fHistScale2EmcalvsCent->Fill(fCent,scalecalcemcal2);     
  fHistChScalevsCent->Fill(fCent,Chscalecalcemcal);    
  fHistChScale2EmcalvsCent->Fill(fCent,Chscalecalcemcal2);   
  fHistScaleEmcalvsNtrack->Fill(Ntracks,scalecalcemcal);      
  fHistScale2EmcalvsNtrack->Fill(Ntracks,scalecalcemcal2);     
  fHistChScalevsNtrack->Fill(Ntracks,Chscalecalcemcal);    
  fHistChScale2EmcalvsNtrack->Fill(Ntracks,Chscalecalcemcal2);   
  fHistPtTPCvsCent->Fill(fCent, ptTPC);
  fHistPtEMCALvsCent->Fill(fCent, ptEMCAL);
  fHistEtvsCent->Fill(fCent, Et);
  fHistScalevsCent->Fill(fCent, scalecalc);
  fHistDeltaScalevsCent->Fill(fCent, scalecalc - scale);
  fHistPtTPCvsNtrack->Fill(Ntracks, ptTPC);
  fHistPtEMCALvsNtrack->Fill(Ntracks, ptEMCAL);
  fHistEtvsNtrack->Fill(Ntracks, Et);
  fHistScalevsNtrack->Fill(Ntracks, scalecalc);
  fHistDeltaScalevsNtrack->Fill(Ntracks, scalecalc - scale);
  fHistScalevsScale2Emcal->Fill(scalecalc,scalecalcemcal2);      
  fHistScalevsScaleEmcal->Fill(scalecalc,scalecalcemcal); 
  fHistScaleEmcalvsScale2Emcal->Fill(scalecalcemcal,scalecalcemcal2);

  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskScale::ExecOnce() 
{
  AliAnalysisTaskEmcal::ExecOnce();

  const Double_t EmcalMinEta = fGeom->GetArm1EtaMin();
  const Double_t EmcalMaxEta = fGeom->GetArm1EtaMax();
  const Double_t EmcalMinPhi = fGeom->GetArm1PhiMin() * TMath::DegToRad();
  const Double_t EmcalMaxPhi = fGeom->GetArm1PhiMax() * TMath::DegToRad();

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
}
