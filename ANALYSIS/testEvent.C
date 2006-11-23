// #####    TO RUN THIS MACRO:
// bash$ aliroot         (due to AliInfo, AliError, ...)
// root[0] gSystem->Load("libANALYSIS_NEW");
// IN case you do not have the include path set to AliRoot includes:
// root[1] gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");
// root[2] .L testEvent.C+;
// root[3] generate();
// root[4] filter_reco();

#include "TClonesArray.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "testEvent.h"   

//============= First step: generate events
void generate()
{
// Simple event generation
   AliAnalysisManager *mgr = new AliAnalysisManager();
   TaskGenerate *task = new TaskGenerate("gener");
   mgr->AddTask(task);

   if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      mgr->ExecAnalysis();
   }   
   delete mgr;
}   

//============= Second step: filter gammas; use TSelector functionality
void filter_reco()
{
// Filter the input events having more than 100 gammas. Reconstruct pi0's
// From gammas coming from the same vertex.
   // Get the input data as a chain
   TChain *chain = new TChain("T");
   chain->Add("event02000.root");
   chain->Add("event04000.root");
   chain->Add("event06000.root");
   chain->Add("event08000.root");
   chain->Add("event10000.root");
   // Create an analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager();
   // Create a filter task and register it
   TaskFilter *task1 = new TaskFilter("TaskFilter");
   mgr->AddTask(task1);
   // Create a reco task and register it
   TaskRecoPi0 *task2 = new TaskRecoPi0("TaskRecoPi0");
   mgr->AddTask(task2);
   // Create containers for input/output
   AliAnalysisDataContainer *cinput = mgr->CreateContainer("input0", 
                      TTree::Class(), AliAnalysisManager::kInputContainer);
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("output0", 
                      TTree::Class(), AliAnalysisManager::kOutputContainer);
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("output1", 
                      TList::Class(), AliAnalysisManager::kOutputContainer);
   AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", 
                      TH1::Class(), AliAnalysisManager::kOutputContainer);
   // Connect containers to the task input/outputs
   mgr->ConnectInput(task1,0,cinput);
   mgr->ConnectOutput(task1,0,coutput1);
   mgr->ConnectOutput(task1,1,coutput2);
   mgr->ConnectInput(task2,0,coutput1);
   mgr->ConnectOutput(task2,0,coutput);
   // Connect input data
   cinput->SetData(chain);
   // Init analysis and start event loop
   if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      chain->Process(mgr);
   }
//   delete mgr;   
}
   
ClassImp(TaskGenerate)

//______________________________________________________________________________
void TaskGenerate::Exec(Option_t *)
{
// Generate 10k events in 5 files
   AnaEvent::CreateEvents(2000, "event02000.root");
   AnaEvent::CreateEvents(2000, "event04000.root");
   AnaEvent::CreateEvents(2000, "event06000.root");
   AnaEvent::CreateEvents(2000, "event08000.root");
   AnaEvent::CreateEvents(2000, "event10000.root");
}

ClassImp(TaskFilter)

//______________________________________________________________________________
TaskFilter::TaskFilter(const char *name) 
           :AliAnalysisTask(name,""), fEvent(0), fOutput(0), fList(0), fHist1(0), fHist2(0)
{
// Constructor
   // Input slot #0 works with a tree of events
   DefineInput(0, TTree::Class());

   // Output slot #0 writes into a tree, #1 produces a hisogram
   DefineOutput(0, TTree::Class());
   DefineOutput(1, TList::Class());
}

//______________________________________________________________________________
void TaskFilter::Init(Option_t *)
{
// Initialize branches.
   printf("   Init %s\n", GetName());
   if (!fEvent) {
      // One should first check if the branch address was taken by some other task
      char ** address = (char **)GetBranchAddress(0, "event");
      if (address) fEvent = (AnaEvent*)(*address);
      if (!fEvent) {
         fEvent = new AnaEvent();
         SetBranchAddress(0, "event", &fEvent);
      }
      // The output tree will be written to gammas.root
      TDirectory *dirsav = gDirectory;
      // Open a file for output #0
      OpenFile(0, "gammas.root", "RECREATE");
      fOutput = new TTree("TGAM", "gammas");
      TBranch *branch = fOutput->Branch("event", &fEvent, 18000,1);
      branch->SetAutoDelete(kFALSE);
      if (dirsav) dirsav->cd();
   } 
   if (!fList) {
      fList = new TList();
      fHist1 = new TH1I("ntracks", "Number of tracks per event", 100, 0, 1000);
      fHist1->SetLineColor(kRed);
      fHist2 = new TH1I("ngammas", "Number of gammas per event", 100, 0, 1000);
      fHist2->SetLineColor(kBlue);
      fList->Add(fHist1);
      fList->Add(fHist2);
   }   
}

//______________________________________________________________________________
void TaskFilter::Exec(Option_t *)
{
// Filtering.
   TTree *tinput = (TTree*)GetInputData(0);
   Long64_t ientry = tinput->GetReadEntry();
   // In this case fEvent address is already connected to the input
   if (!fEvent) return;
   // First check multiplicity
   Int_t ntracks = fEvent->GetNtracks();
   // Loop tracks and get rid of non-gammas
   AnaTrack *track;
   Int_t igamma = 0;
   for (Int_t i=0; i<ntracks; i++) {
      track = fEvent->GetTrack(i);
      if (track->GetMass() < 1.e-3) igamma++;
   }
   if (ientry%100 == 0) printf("TaskFilter -> Event %lld: %i tracks filtered %i gammas\n", ientry, ntracks, igamma);
   fHist1->Fill(ntracks);
   fHist2->Fill(igamma);
   if (igamma > 100) {
      fOutput->Fill();
      PostData(0, fOutput);
   }   
   PostData(1, fList);
}

//______________________________________________________________________________
void TaskFilter::Terminate(Option_t *)
{
// Draw some histogram at the end.
   if (!gROOT->IsBatch()) {
      fHist1->SetMaximum(2500);
      fHist1->Draw();
      fHist2->Draw("SAME");
   }   
}

ClassImp(TaskRecoPi0)

//______________________________________________________________________________
TaskRecoPi0::TaskRecoPi0(const char *name) 
           :AliAnalysisTask(name,""), fEvent(0), fGammas(0), fPions(0), fHist(0)
{
// Constructor
   // Input slot #0 works with a tree of events
   DefineInput(0, TTree::Class());

   // Output slot #1 produces a hisogram
   DefineOutput(0, TH1::Class());
}

//______________________________________________________________________________
TaskRecoPi0::~TaskRecoPi0()
{
// Dtor.
   if (fEvent) delete fEvent;
   if (fGammas) delete fGammas;
   if (fPions) delete fPions;
}   

//______________________________________________________________________________
void TaskRecoPi0::Init(Option_t *)
{
// Initialize branches.
   printf("   Init %s\n", GetName());
   if (!fEvent) {
      // One should first check if the branch address was taken by some other task
      char ** address = (char **)GetBranchAddress(0, "event");
      if (address) fEvent = (AnaEvent*)(*address);
      if (!fEvent) {
         fEvent = new AnaEvent();
         SetBranchAddress(0, "event", &fEvent);
      }
      fGammas = new TObjArray();
      fPions  = new TObjArray();
   } 
   if (!fHist) {
      fHist = new TH1F("Pt_pi0", "Pt distribution for pi0's", 100, 0., 10.);
      fHist->SetLineColor(kRed);
   }   
}

//______________________________________________________________________________
void TaskRecoPi0::Exec(Option_t *)
{
// Reconstruct Pi0's for one event
   AnaTrack *track = 0;
   Int_t ntracks = fEvent->GetNtracks();
   Int_t ngamma = 0;
   Int_t i,j;
   // Clear containers
   fGammas->Clear();
   fPions->Delete();
   fPions->Clear();
   // Loop tracks and move gammas to gamma container
   for (i=0; i<ntracks; i++) {
      track = fEvent->GetTrack(i);
      if (track->GetMass() < 1.e-3) {
         fGammas->Add(track);
         ngamma++;
      }
   }
   printf("TaskRecoPi0 -> Tracks %i \n", ntracks);

   // Loop gammas and check vertex position
   Double_t v1[3], v2[3];
   AnaTrack *tracko = 0;
   Double_t cutoff = 0.001;
   for (i=0; i<ngamma-1; i++) {
      track = (AnaTrack*)fGammas->At(i);
      v1[0] = track->GetVertex(0);
      v1[1] = track->GetVertex(1);
      v1[2] = track->GetVertex(2);
      for (j=i+1; j<ngamma; j++) {
         tracko = (AnaTrack*)fGammas->At(j);
         v2[0] = tracko->GetVertex(0);
         v2[1] = tracko->GetVertex(1);
         v2[2] = tracko->GetVertex(2);
         Double_t dist2 = (v2[0]-v1[0])*(v2[0]-v1[0])+
                          (v2[1]-v1[1])*(v2[1]-v1[1])+
                          (v2[2]-v1[2])*(v2[2]-v1[2]);
         if (dist2>cutoff*cutoff) continue;
         // We have found a pair candidate
         Double_t px = track->GetPx()+tracko->GetPx();
         Double_t py = track->GetPy()+tracko->GetPy();
         Double_t pz = track->GetPz()+tracko->GetPz();
         track = new AnaTrack(px,py,pz,0.135,0,v1[0],v1[1],v1[2]);
         fPions->Add(track);
         fHist->Fill(track->GetPt());
         break;
      }   
   }   
   PostData(0,fHist);         
}

//______________________________________________________________________________
void TaskRecoPi0::Terminate(Option_t *)
{
// Draw some histogram at the end.
   if (!gROOT->IsBatch()) {
      new TCanvas("pi0", "Pt for pi0's", 800,600);
      fHist->Draw();
   }   
}


ClassImp(AnaTrack)

//______________________________________________________________________________
AnaTrack::AnaTrack(Double_t random, Double_t *vertex) 
{
// Constructor
   Int_t itype;
   if (random<0.3) itype = 1;      // pi+
   else if (random<0.6) itype = 2; // pi-
   else if (random<0.9) itype = 3; // p
   else if (random<0.95) itype = 4; // pi0
   else itype = 5;                 // gamma
   gRandom->Rannor(fPx, fPy);
   Double_t vert_width = 0.1;
   
   switch (itype) {
      case 1:
         fMass = 0.13957 + gRandom->Gaus(0.,0.001);
         fCharge = 1;
         break;
      case 2:
         fMass = 0.13957 + gRandom->Gaus(0.,0.001);
         fCharge = -1;
         break;
      case 3:
         fMass = 0.938 + gRandom->Gaus(0.,0.002);
         fCharge = 1;
         fPx *= 0.15;
         fPy *= 0.15;
         break;
      case 4:
         fMass = 0.135 + gRandom->Gaus(0.,0.001);
         fCharge = 0;
         fPx *= 0.8;
         fPy *= 0.8;
         vert_width = 10.;
         break;
      case 5:
         fMass = 0.;
         fCharge = 0;
         fPx *= 0.5;
         fPy *= 1.5;
         break;
   };
   fPz = gRandom->Gaus(4., 2.);
   if (vertex) {
      fVertex[0] = vertex[0];
      fVertex[1] = vertex[1];
      fVertex[2] = vertex[2];
   } else {   
      fVertex[0] = gRandom->Gaus(0,vert_width);
      fVertex[1] = gRandom->Gaus(0,vert_width);
      fVertex[2] = gRandom->Gaus(0,0.01);
   }   
}

//______________________________________________________________________________
Bool_t AnaTrack::Decay(Double_t &px1, Double_t &py1, Double_t &pz1, 
                       Double_t &px2, Double_t &py2, Double_t &pz2)
{
// Decay a pi0 in 2 gammas.
   if (fCharge != 0) return kFALSE;
   if (fMass<0.132 || fMass>0.138) return kFALSE;
   Double_t phi1 = 2.*TMath::Pi()*gRandom->Rndm(); // [0,2*pi]
   Double_t phi2 = phi1 + TMath::Pi();
   if (phi2 > 2.*TMath::Pi()) phi2 -= 2.*TMath::Pi();
   Double_t r2 = gRandom->Rndm();
   Double_t theta1 = TMath::ACos(1.-2.*r2);
   Double_t p0 = GetP();
   Double_t m0 = 0.135;
   Double_t p1 = 0.5*(p0*(1-2*r2)+TMath::Sqrt(p0*p0*(1-2*r2)*(1-2*r2)+2*m0*m0));
   Double_t p2 = TMath::Sqrt(p0*p0+m0*m0-p1*p1);
   Double_t theta2 = TMath::ACos((p0-p1*(1-2*r2))/p2);
   // Px, Py and Pz in the reference frame of the pion
   px1 = p1 * TMath::Sin(theta1)*TMath::Cos(phi1);
   px2 = p2 * TMath::Sin(theta2)*TMath::Cos(phi2);
   py1 = p1 * TMath::Sin(theta1)*TMath::Sin(phi1);
   py2 = p2 * TMath::Sin(theta2)*TMath::Sin(phi2);
   pz1 = p1 * TMath::Cos(theta1);
   pz2 = p2 * TMath::Cos(theta2);
   Double_t phi = TMath::ATan2(fPy, fPx) * TMath::RadToDeg();
   Double_t theta = TMath::ACos(GetPt()/p0) * TMath::RadToDeg();
   TGeoRotation r("rot", phi+90., theta, 0.);
   Double_t loc[3], vect[3];
   p0 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
   loc[0] = px1/p0;
   loc[1] = py1/p0;
   loc[2] = pz1/p0;
   r.LocalToMasterVect(loc, vect);
   px1 = vect[0]*p0;
   py1 = vect[1]*p0;
   pz1 = vect[2]*p0;
//   t1 = new AnaTrack(1., fVertex);
   p0 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);
   loc[0] = px2/p0;
   loc[1] = py2/p0;
   loc[2] = pz2/p0;
   r.LocalToMasterVect(loc, vect);
   px2 = vect[0]*p0;
   py2 = vect[1]*p0;
   pz2 = vect[2]*p0;
//   t2 = new AnaTrack(1., fVertex);
   return kTRUE;
}

//______________________________________________________________________________
ClassImp(AnaEvent)

TClonesArray *AnaEvent::fgTracks = 0;

//______________________________________________________________________________
AnaEvent::AnaEvent()
{
// Ctor
   if (!fgTracks) fgTracks = new TClonesArray("AnaTrack", 2000);
   fTracks = fgTracks;
   fEventNumber = 0;
   fNtracks = 0;
}

//______________________________________________________________________________
AnaTrack *AnaEvent::AddTrack(Double_t rnd, Double_t *vert)
{
// Add a random track
//   printf("track %d\n", fNtracks);
   TClonesArray &tracks = *fTracks;
   AnaTrack *track = new(tracks[fNtracks++]) AnaTrack(rnd, vert);
   return track;
}   

//______________________________________________________________________________
void AnaEvent::Clear(Option_t *)
{
// Clears current event.
   fTracks->Clear();
   fNtracks = 0;
   fEventNumber = 0;
}   

//______________________________________________________________________________
Int_t AnaEvent::Build(Int_t ev)
{
// Create a random event
   Clear();
   fEventNumber = ev;
   Int_t ntracks = Int_t(gRandom->Gaus(500., 100.));
   if (ntracks < 1) ntracks = 1;
   if (ntracks>1000) ntracks = 1000;
   AnaTrack *track, *track0;
   for (Int_t i=0; i<ntracks; i++) {
      Double_t rnd = gRandom->Rndm();
      if (rnd>=0.90 && rnd<0.95) {
      // Pi0 -> decay the track in 2 gammas
         track0 = new AnaTrack(rnd);
         Double_t vert[3];
         vert[0] = track0->GetVertex(0);
         vert[1] = track0->GetVertex(1);
         vert[2] = track0->GetVertex(2);
         Double_t px1,py1,pz1,px2,py2,pz2;
         if (track0->Decay(px1,py1,pz1,px2,py2,pz2)) {
            track = AddTrack(1.,vert);
            track->SetPx(px1);
            track->SetPy(py1);
            track->SetPz(pz1);
            track = AddTrack(1.,vert);
            track->SetPx(px2);
            track->SetPy(py2);
            track->SetPz(pz2);
         }
         delete track0;
      } else {   
         track = AddTrack(rnd);
      }   
   }
   return fNtracks;
}   

//______________________________________________________________________________
void AnaEvent::CreateEvents(Int_t nevents, const char *filename)
{
// Create nevents in one tree and put them in filename.
   TFile *hfile = new TFile(filename, "RECREATE", "Some AnaEvents...");
   TTree *tree  = new TTree("T", "Tree of AnaEvents");
   tree->SetAutoSave(1000000000);  // autosave when 1 Gbyte written
   Int_t bufsize = 16000;
   AnaEvent *event = new AnaEvent();
   TBranch *branch = tree->Branch("event", &event, bufsize,1);
   branch->SetAutoDelete(kFALSE);
   
   for (Int_t ev=0; ev<nevents; ev++) {
      Int_t ntracks = event->Build(ev);
      if (ev%100 == 0) printf("event: %d  ntracks=%d\n", ev, ntracks);
      tree->Fill();
   }   
   hfile->Write();
   tree->Print();
   hfile->Close();
}  
   
