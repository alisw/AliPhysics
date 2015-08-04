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

// Class that can be plugged in the simulation to monitor transport timing per 
// particle for each geometry volume.
//
//  andrei.gheata@cern.ch 

#include "AliTransportMonitor.h"

#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "THashList.h"
#include "TH2F.h"
#include "TGeoManager.h"
#include "AliPDG.h"
#include "TVirtualMC.h"
#include <RVersion.h>

ClassImp(AliTransportMonitor)
ClassImp(AliTransportMonitor::AliTransportMonitorVol)
ClassImp(AliTransportMonitor::AliTransportMonitorVol::AliPMonData)

typedef AliTransportMonitor::AliTransportMonitorVol::AliPMonData PMonData;

//______________________________________________________________________________
AliTransportMonitor::AliTransportMonitorVol::AliTransportMonitorVol()
                    :TNamed(),
                     fNtypes(0),
                     fTotalTime(0),
                     fNSteps(0),
                     fPData(0),
                     fTimeRZ(0),
                     fParticles()
{
// Default constructor
}

//______________________________________________________________________________
AliTransportMonitor::AliTransportMonitorVol::~AliTransportMonitorVol()
{
// Destructor
  delete [] fPData;
  delete fTimeRZ;
}

//______________________________________________________________________________
void AliTransportMonitor::AliTransportMonitorVol::StepInfo(
                                    Int_t pdg,
                                    Double_t energy, 
                                    Double_t dt,
                                    Double_t x, Double_t y, Double_t z)
{
// This method is called at each N steps to store timing info.
  PMonData &data = GetPMonData(pdg);
  data.fEdt += energy*dt;
  data.fTime += dt;
  fTotalTime += dt;
  fNSteps += 1.;
  if (!fTimeRZ) {
    Bool_t status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    fTimeRZ = new TH2F("h2rz", "", 100, -5000., 5000, 100, 0., 1000.);
    TH1::AddDirectory(status);
  }
  Double_t r = TMath::Sqrt(x*x+y*y);
  fTimeRZ->Fill(z,r,dt);
}

//______________________________________________________________________________
PMonData &AliTransportMonitor::AliTransportMonitorVol::GetPMonData(Int_t pdg)
{
// Retrieve stored monitoring object for a given pdg type. If not existing 
// create one.

  // The object could have been retrieved from file, in which case we have to 
  // build the map.
  //
  // unknown heavy fragment ?
  //  TParticlePDG* pdgP = (TDatabasePDG::Instance())->GetParticle(pdg);
  Int_t apdg = TMath::Abs(pdg);
  if ((apdg > 10000) 
      && (apdg != 1000010020)
      && (apdg != 1000010030)
      && (apdg != 1000020030)
      && (apdg != 1000020040)
      && (apdg != 50000050)
      && (apdg != 50000051)
      ) 
    pdg = 1111111111; 

  PMonData *data;
  if (fNtypes) {
    if (fParticles.empty()) {
      for (Int_t i=0; i<fNtypes; i++) {
        data = &fPData[i];
        fParticles[i] = data->fPDG;
      }
    }
    ParticleMapIt_t it = fParticles.find(pdg);
    if (it == fParticles.end()) {
      data = &fPData[fNtypes]; 
      data->fPDG = pdg;
      fParticles[pdg] = fNtypes++;
    } else {
      data = &fPData[it->second];
    }  
  } else {
     TDatabasePDG *pdgDB = TDatabasePDG::Instance();
     if (!pdgDB->ParticleList()) AliPDG::AddParticlesToPdgDataBase();
     Int_t size = pdgDB->ParticleList()->GetSize();
     // account for heavy fragments coded as "1111111111"
     //
     fPData = new PMonData[size+10];
     data = &fPData[fNtypes];
     data->fPDG = pdg;
     fParticles[pdg] = fNtypes++;
  }
  return *data;
}

void AliTransportMonitor::AliTransportMonitorVol::Merge(AliTransportMonitorVol* volM) 
{
  //
  // Merging
  //
  fTotalTime = (fTotalTime + volM->GetTotalTime());
  if (fTimeRZ && volM->GetHistogram()) {
    fTimeRZ->Add(volM->GetHistogram()); 
  } else if (volM->GetHistogram()) {
    fTimeRZ = (TH2F*)(volM->GetHistogram()->Clone());
  }

  Int_t ntypes = volM->GetNtypes();
  for (Int_t i = 0; i < ntypes; i++) {
    Int_t pdg = volM->GetPDG(i);
     PMonData &data  = GetPMonData(pdg);
     data.fEdt  += (volM->GetEmed(i) * volM->GetTotalTime());
     data.fTime += (volM->GetTime(i));
  }
}

ClassImp(AliTransportMonitor)

//______________________________________________________________________________
AliTransportMonitor::AliTransportMonitor()
                    :TObject(),
                     fTotalTime(0),
                     fTimer(),
                     fVolumeMon(0)
{
// Default constructor
}

//______________________________________________________________________________
AliTransportMonitor::AliTransportMonitor(Int_t nvolumes)
                    :TObject(),
                     fTotalTime(0),
                     fTimer(),
                     fVolumeMon(0)
{
// Default constructor
  fVolumeMon = new TObjArray(nvolumes);
  fVolumeMon->SetOwner();
  for (Int_t i=0; i<nvolumes; i++) {
    AliTransportMonitorVol *volMon = new AliTransportMonitorVol();
    if (TVirtualMC::GetMC()) volMon->SetName(TVirtualMC::GetMC()->VolName(i));
    fVolumeMon->Add(volMon);
  }   
}

//______________________________________________________________________________
AliTransportMonitor::~AliTransportMonitor()
{
// Destructor
  delete fVolumeMon;
}

//______________________________________________________________________________
void AliTransportMonitor::Print(Option_t *volName) const
{
// Inspect the timing statistics for a single volume or for all the setup
  Int_t uid = -1;
  Int_t ntotal = 0;
  if (!fVolumeMon || !(ntotal=fVolumeMon->GetEntriesFast())) {
     Info("Inspect", "Transport monitor is empty !");
     return;
  }   
  if (strlen(volName)) {
    TString svname = volName;
    Int_t i = 0;
    for (i=1; i<ntotal; i++) if (svname == fVolumeMon->At(i)->GetName()) break;
    if (i==ntotal) {
       Error("Inspect", "No monitoring info stored for volume %s", volName);
       return;
    }
    uid = i;
    AliTransportMonitorVol *volMon = (AliTransportMonitorVol*)fVolumeMon->At(uid);
    Int_t ntypes = volMon->GetNtypes();
    if (!ntypes) {
      Info("Inspect", "No particles crossed volume %s", volName);
      return;
    }  
    Double_t *timeperpart = new Double_t[ntypes];
    Int_t *isort = new Int_t[ntypes];
    Double_t timepervol = 0.;
    for (i=0; i<ntypes; i++) {
      timeperpart[i] = volMon->GetTime(i); 
      timepervol += timeperpart[i];
    }
    printf("Volume %s: Transport time: %g%% of %g %g [s]\n", volMon->GetName(), 100.*timepervol/fTotalTime, fTotalTime,
	   volMon->GetTotalTime());
    TMath::Sort(ntypes, timeperpart, isort, kTRUE);
    TString particle;
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();    
    if (!pdgDB->ParticleList()) AliPDG::AddParticlesToPdgDataBase();
    for (i=0; i<ntypes; i++)  {
       timeperpart[i] /=  timepervol;
       TParticlePDG* pdgP =  pdgDB->GetParticle(volMon->GetPDG(isort[i]));
       if (pdgP) {
	 particle = pdgDB->GetParticle(volMon->GetPDG(isort[i]))->GetName();
       } else {
	 particle = Form("pdg code not in DB: %d", volMon->GetPDG(isort[i]));
       }
       printf("   %s: %g%%  mean energy: %g\n", particle.Data(), 100.*timeperpart[i], volMon->GetEmed(i));
    }
    if (volMon->GetHistogram()) {
      TCanvas *c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("crz");
      if (!c1) c1 = new TCanvas("crz");
      c1->cd();
      volMon->GetHistogram()->GetXaxis()->SetTitle("z [cm]");
      volMon->GetHistogram()->GetYaxis()->SetTitle("r [cm]");
      volMon->GetHistogram()->SetTitle(Form("RZ plot weighted by time spent in %s",volMon->GetName()));
      volMon->GetHistogram()->Draw();
    }  
    return;
  }
  // General view
  TIter next(fVolumeMon);
  AliTransportMonitorVol *volMon;
  Int_t ncrossed = 0;
  TH1F *hnames = new TH1F("volume_timing", "relative volume timing", 3,0,3);
  hnames->SetStats(0);
  hnames->SetFillColor(38);
  #if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  hnames->SetBit(TH1::kCanRebin);
  #endif
  while ((volMon=(AliTransportMonitorVol*)next())) {
    if (volMon->GetNtypes()) {
      hnames->Fill(volMon->GetName(), volMon->GetTotalTime());
      ncrossed++;
    }
  }
  
  hnames->LabelsDeflate();
  hnames->GetXaxis()->LabelsOption(">");
  
  TCanvas *c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cvol_timing");
  if (!c) c = new TCanvas("cvol_timing");
  c->cd();
  c->SetLogy();
  hnames->Draw();
  
  printf("=============================================================================\n");
  printf("Effective transport time:  %6.2f minutes\n", fTotalTime/60.);
  printf("Number of crossed volumes: %d from %d\n", ncrossed, fVolumeMon->GetEntriesFast());
  printf("=============================================================================\n");  
}     

//______________________________________________________________________________
void AliTransportMonitor::DummyStep()
{
// Reset timer for zero-length steps
   fTimer.Stop();
   fTimer.Start(kTRUE);
}   

//______________________________________________________________________________
void AliTransportMonitor::StepInfo( Int_t volId,
                                    Int_t pdg,
                                    Double_t energy, 
                                    Double_t x, Double_t y, Double_t z)
{
// This method is called at each N steps to store timing info.
  fTimer.Stop();
  Double_t dt = fTimer.RealTime();
  fTotalTime += dt;
  AliTransportMonitorVol *volMon = (AliTransportMonitorVol*)fVolumeMon->At(volId);
  volMon->StepInfo(pdg,energy,dt,x,y,z);
  fTimer.Start(kTRUE);
}

//______________________________________________________________________________
void AliTransportMonitor::Start()
{
// Start collecting timing information
  if (fTotalTime > 0) {
    Info("Start", "Cannot start twice");
    return;
  }
  fTimer.Start(kTRUE);
}

//______________________________________________________________________________
void AliTransportMonitor::Stop()
{
// Stop the global timer
  fTimer.Stop();
}

//______________________________________________________________________________
void AliTransportMonitor::Export(const char *fname)
{
// Export information to file
  TFile *file = TFile::Open(fname, "RECREATE");
  Write();
  file->Write();
  file->Close();
}  
//______________________________________________________________________________
void AliTransportMonitor::Merge(AliTransportMonitor* mergeMon)
{

  // merge with monitor 
  if (!fVolumeMon) 
    {
      TObjArray* arr = mergeMon->GetVolumes();
      Int_t nvol = arr->GetEntriesFast();
      fVolumeMon = new TObjArray(nvol);
      fVolumeMon->SetOwner();
      for (Int_t i = 0; i < nvol; i++) {
	AliTransportMonitorVol *volMon = new AliTransportMonitorVol();
	volMon->SetName(arr->At(i)->GetName());
	fVolumeMon->Add(volMon);
      }
    } // first time


  Int_t n = fVolumeMon->GetEntriesFast();
  TObjArray* mergeVols = mergeMon->GetVolumes();
  fTotalTime = 0;
  for (Int_t i = 0; i < n; i++)
    {
      AliTransportMonitorVol *volMon1 = (AliTransportMonitorVol*)fVolumeMon->At(i);      
      AliTransportMonitorVol *volMon2 = (AliTransportMonitorVol*)mergeVols->At(i);      
      volMon1->Merge(volMon2);
      fTotalTime += (volMon1->GetTotalTime());
    }
}
//______________________________________________________________________________
AliTransportMonitor *AliTransportMonitor::Import(const char *fname)
{
// Import information from a file
  TFile *file = TFile::Open(fname);
  if (!file) {
    ::Error("Import", "File %s could not be opened", fname);
    return 0;
  }
  AliTransportMonitor *mon = (AliTransportMonitor *)file->Get("AliTransportMonitor");
  if (!mon) {
    ::Error("Import", "No AliTransportMonitor object found n file %s", fname);
    return 0;
  }
  AliPDG::AddParticlesToPdgDataBase();
  return mon;
}

    
