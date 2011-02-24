#define AliAnalysisTaskMuonAODfromGeneral_cxx

/* $Id$ */

// 19 Nov 2007
// Class implementation for the specific muon AOD generation
// Extracts only muon tracks from a general AOD and builds dimuons
// Livio Bianchi, Universita' di Torino


#include "TTree.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom.h"

#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskMuonAODfromGeneral.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODEventInfo.h"
#include "AliAODDimuon.h"

ClassImp(AliAnalysisTaskMuonAODfromGeneral)

//________________________________________________________________________
AliAnalysisTaskMuonAODfromGeneral::AliAnalysisTaskMuonAODfromGeneral(const char *name, Double_t BeamEnergy):
	AliAnalysisTask(name, "AnalysisTaskMuonAODfromGeneral"),
	fInfos(0),
	fDimuons(0),
        fChain(0),
        fOrgAOD(0), 
	fNewAOD(0), 
	ft(0),
	fBeamEnergy(0)
{  // Constructor.
  // Input slot #0 works with an Ntuple
  SetBeamEnergy(BeamEnergy);
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TTree container
  DefineOutput(0, TTree::Class());
}

//________________________________________________________________________
void AliAnalysisTaskMuonAODfromGeneral::ConnectInputData(Option_t *) {
  printf("   ConnectInputData %s\n", GetName());
  // New connection to the AOD
  fChain = (TChain*)GetInputData(0);
  fOrgAOD = new AliAODEvent();
  fOrgAOD->ReadFromTree(fChain); // This checks also if the branch address is already set
}

//________________________________________________________________________
void AliAnalysisTaskMuonAODfromGeneral::CreateOutputObjects() {
  printf("Creating output objects\n");
  if(!fNewAOD)fNewAOD=new AliAODEvent();
  fNewAOD->CreateStdContent();
 
  if (!ft) {
    ft = new TTree("AOD","Muon AOD tree",0);
    ft->BranchRef();
    ft->Branch(fNewAOD->GetList());
    //ft->Print();
 }
 
  //Add new Objects: MuonInfos & Dimuons
   fInfos = new AliAODEventInfo();
   fNewAOD->AddObject(fInfos);
   const char *nameInfos = "MuonInfos";
   fInfos->SetName(nameInfos);
   ft->Branch(nameInfos, &fInfos);
   fDimuons = new TClonesArray("AliAODDimuon",0);
   fNewAOD->AddObject(fDimuons);
   const char *nameDimuons = "Dimuons";
   fDimuons->SetName(nameDimuons);
   ft->Branch(nameDimuons, &fDimuons);
      
}

//________________________________________________________________________
void AliAnalysisTaskMuonAODfromGeneral::Exec(Option_t *) {
  static int ncall=0;
  // Task making a Muon AOD
  // Get input data
  TChain *chain = (TChain*)GetInputData(0);
  Long64_t ientry = chain->GetReadEntry();
  ientry=ientry;

  if (!fOrgAOD) return;
  
  Int_t nTracks=fOrgAOD->GetNumberOfTracks();
  Int_t nMuTracks=0;
  Int_t nPosTracks = 0;
  Int_t mutrNumb[10]; for (Int_t i=0; i<10; mutrNumb[i++]=0) {}
  for (Int_t iTrack=0; iTrack<nTracks; iTrack++){
       const Double_t *trackpid=fOrgAOD->GetTrack(iTrack)->PID();
       if (trackpid[AliAODTrack::kMuon]==1.) {
       		mutrNumb[nMuTracks]=iTrack;
		nMuTracks++;
		}
       if (fOrgAOD->GetTrack(iTrack)->Charge()> 0) nPosTracks++;
       }
  
  Bool_t ExistMuon=kFALSE;
  if (nMuTracks>0) ExistMuon=kTRUE;
  
  //--------------------------------------------------------------------
  // set arrays and pointers
  Double_t pos[3];
  Double_t covVtx[6];
  
    // Access to the header
    AliAODHeader *header = fNewAOD->GetHeader();

  // fill the header
  header->SetRunNumber       (fOrgAOD->GetRunNumber()       );
  header->SetBunchCrossNumber(fOrgAOD->GetBunchCrossNumber());
  header->SetOrbitNumber     (fOrgAOD->GetOrbitNumber()     );
  header->SetPeriodNumber    (fOrgAOD->GetPeriodNumber()    );
  header->SetTriggerMask     (fOrgAOD->GetTriggerMask()     ); 
  header->SetTriggerCluster  (fOrgAOD->GetTriggerCluster()  );
  header->SetEventType       (fOrgAOD->GetEventType()       );
  header->SetMagneticField   (fOrgAOD->GetMagneticField()   );
  header->SetZDCN1Energy     (fOrgAOD->GetZDCN1Energy()     );
  header->SetZDCP1Energy     (fOrgAOD->GetZDCP1Energy()     );
  header->SetZDCN2Energy     (fOrgAOD->GetZDCN2Energy()     );
  header->SetZDCP2Energy     (fOrgAOD->GetZDCP2Energy()     );
  header->SetZDCEMEnergy     (0,fOrgAOD->GetZDCEMEnergy(0)  );
  header->SetZDCEMEnergy     (1,fOrgAOD->GetZDCEMEnergy(1)  );
  header->SetRefMultiplicity   (nTracks);
  header->SetRefMultiplicityPos(nPosTracks);
  header->SetRefMultiplicityNeg(nTracks - nPosTracks);
  header->SetMuonMagFieldScale(-999.); // FIXME
  header->SetCentrality(0);            // FIXME

  const AliAODVertex *vtx = fOrgAOD->GetPrimaryVertex();
  if ( !vtx) {
    // CHECK Gines
    AliError("Primary vertex is not defined");
    return;
  }  

  Int_t nDimuons=nMuTracks*(nMuTracks-1)/2;
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(fNewAOD->GetVertices());
  Int_t jVertices=0;

  // Access to the AOD container of tracks
   TClonesArray &tracks = *(fNewAOD->GetTracks());
   tracks.Delete();
   Int_t jTracks=0; 


  vtx->GetCovMatrix(covVtx); //covariance matrix

  AliAODVertex * primary = new(vertices[jVertices++])
      AliAODVertex(pos, covVtx, vtx->GetChi2perNDF(), NULL, AliAODVertex::kPrimary);

  static int ncal=0;
  static int numtracks=0;
  for (Int_t iTrack=0; iTrack<nMuTracks; iTrack++) {
	primary->AddDaughter(new(tracks[jTracks++])
	AliAODTrack((*(fOrgAOD->GetTrack(mutrNumb[iTrack])))));
        ncal++;
  }
  numtracks+=jTracks;

  fDimuons->Delete();
  fDimuons->Expand(nDimuons);
  TClonesArray &rDimuons = *fDimuons;
    
  Int_t jDimuons=0;
  for(Int_t i=0; i<nMuTracks; i++){
     for(Int_t j=i+1; j<nMuTracks; j++){
       new(rDimuons[jDimuons++]) AliAODDimuon(tracks[i],tracks[j]);
     } 
  }
  
  fInfos->SetBeamEnergy(fBeamEnergy);
  fInfos->SetEv(fNewAOD);
  fInfos->SetHe(header);
  fInfos->SetTr(fNewAOD->GetTracks());
  fInfos->SetDi(fDimuons);
  fInfos->SelectTriggerBits(0,1,2,3,4,5);
  if(ExistMuon) ft->Fill();
  ncall++;

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, ft);
}      

//________________________________________________________________________
void AliAnalysisTaskMuonAODfromGeneral::Terminate(Option_t *) {
  // Draw some histogram at the end.
  ft->GetUserInfo()->Add(fNewAOD);
  ft->Write("",TObject::kOverwrite);
  ft->ResetBranchAddresses();
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Pt",10,10,310,310);
    c1->SetFillColor(10);
    c1->SetHighLightColor(10);
    
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);  
    c1->cd(1)->SetLogy();
    ft->Draw("tracks.Pt()");
    //ft->Draw("dimuons.M()");
  }
}
