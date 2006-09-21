// ////////////////////////////////////////////////////////////////////////////
//   Macro to get TRD signal distribution from all six layers of TRD
//   This function takes the data from following directories 
//   "rmom.6", "rmom.8", "rmom1", "rmom1.5", "rmom2", "rmom3", "rmom4", "rmom5",//   "rmom6", "rmom8", "rmom10"
//    corresponding to momenta
//   0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0
//   Running this macro will generate the root file TRDdEdxHistogramsV1.root 
//   containing histograms for all momenta.
//   Prashant Shukla (shukla@physi.uni-heidelberg.de)
// ////////////////////////////////////////////////////////////////////////////

#if !defined( __CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TTree.h>

#include "AliTRDtrack.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliESD.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#else
#endif


//_____________________________________________________________________________

Bool_t AliTRDhistTRDpid(Int_t oFile=1)
{
  //
  //

  const Int_t kNMom=11;

  TH1F *h1dEdxEL[kNMom];
  TH1F *h1dEdxPI[kNMom];
  TH1F *h1dEdxMU[kNMom];
  TH1F *h1dEdxKA[kNMom];
  TH1F *h1dEdxPR[kNMom];
  TH1F *h1MaxTimBinEL[kNMom];
  TH1F *h1MaxTimBinPI[kNMom];
  Char_t text[200];
 
  // Histograms for dedx of e and pi track in TRD 
  
  for (Int_t imom = 0; imom < kNMom; imom++) {
    sprintf(text,"h1dEdxEL%01d",imom+1);
    h1dEdxEL[imom] = new TH1F(text,"dE/dx  Dist.",200,0,4000);
    h1dEdxEL[imom]->GetXaxis()->SetTitle("dE/dx_{TRD}(arbitrary units)");
    h1dEdxEL[imom]->GetYaxis()->SetTitle("Entries");
    h1dEdxEL[imom]->SetLineColor(kRed);
    
    sprintf(text,"h1dEdxPI%01d",imom+1);  
    h1dEdxPI[imom] = new TH1F(text,"dE/dx  Dist.",200,0,4000);
    //  TH1F *h1dEdxPI[imom] = new TH1F("","",250,0,2500);
    h1dEdxPI[imom]->GetXaxis()->SetTitle("dE/dx_{TRD} (arbitrary units)");
    h1dEdxPI[imom]->GetYaxis()->SetTitle("Entries");
    h1dEdxPI[imom]->SetLineColor(kBlue);

    sprintf(text,"h1dEdxMU%01d",imom+1);  
    h1dEdxMU[imom] = new TH1F(text,"dE/dx  Dist.",200,0,4000);
    h1dEdxMU[imom]->GetXaxis()->SetTitle("dE/dx_{TRD}(arbitrary units)");
    h1dEdxMU[imom]->GetYaxis()->SetTitle("Entries");
    h1dEdxMU[imom]->SetLineColor(kGreen);
    
    sprintf(text,"h1dEdxKA%01d",imom+1);  
    h1dEdxKA[imom] = new TH1F(text,"dE/dx  Dist.",200,0,4000);
    h1dEdxKA[imom]->GetXaxis()->SetTitle("dE/dx_{TRD}(arbitrary units)");
    h1dEdxKA[imom]->GetYaxis()->SetTitle("Entries");
    h1dEdxKA[imom]->SetLineColor(7);
    
    sprintf(text,"h1dEdxPR%01d",imom+1); 
    h1dEdxPR[imom] = new TH1F(text,"dE/dx  Dist.",200,0,4000);
    h1dEdxPR[imom]->GetXaxis()->SetTitle("dE/dx_{TRD}(arbitrary units)");
    h1dEdxPR[imom]->GetYaxis()->SetTitle("Entries");
    h1dEdxPR[imom]->SetLineColor(6);
    
    sprintf(text,"h1MaxTimBinEL%01d",imom+1);  
    h1MaxTimBinEL[imom] = new TH1F(text,"Dist. Time Bin of max. Cluster for e(Red) and pi(Blue)",30,0,30);
    h1MaxTimBinEL[imom]->GetXaxis()->SetTitle("Time Bin of Maximum Cluster");
    h1MaxTimBinEL[imom]->GetYaxis()->SetTitle("Entries");
    h1MaxTimBinEL[imom]->SetLineColor(2);
    
    sprintf(text,"h1MaxTimBinPI%01d",imom+1);  
    h1MaxTimBinPI[imom] = new TH1F(text,"Time Bin of max. Cluster Pion",30,0,30);
    h1MaxTimBinPI[imom]->SetLineColor(4);
  }
  
  // PID
  Int_t partCode[5] = 
    {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  
  // Files

  Char_t *FileDir[11] = {"rmom.6", "rmom.8", "rmom1", "rmom1.5", "rmom2", "rmom3", "rmom4", "rmom5", "rmom6", "rmom8", "rmom10"}; 

  Float_t mome[11]={0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

  Char_t gAliceFileName[200];
  Char_t esdFileName[200];
  for (Int_t imom = 0; imom < kNMom; imom++) {

    sprintf(gAliceFileName,"~/work/data/%s/galice.root", FileDir[imom]);
    sprintf(esdFileName,"~/work/data/%s/AliESDs.root", FileDir[imom]);

  // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open(gAliceFileName);
  if (!runLoader) {
    Error("CheckESD", "getting run loader from file %s failed",gAliceFileName);
    return kFALSE;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("CheckESD", "no galice object found");
    return kFALSE;
  }
  runLoader->LoadKinematics();
  runLoader->LoadHeader();
  
  //  AliKalmanTrack::SetConvConst(1000/0.299792458/gAlice->Field()->SolenoidField());
  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckESD", "opening ESD file %s failed",esdFileName);
    return kFALSE;
  }
  AliESD* esd = new AliESD;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  tree->SetBranchAddress("ESD", &esd);

  // loop over events
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    AliStack* stack = gAlice->Stack();
    Int_t nParticles = stack->GetNtrack();
    TArrayF vertex(3);
    runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
    Int_t nGene=0, nGenpi=0;
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = stack->Particle(iParticle);
      if (!particle) continue;
      if (particle->Pt() < 0.001) continue;
      //      if (TMath::Abs(particle->Eta()) > 0.3) continue;
      TVector3 dVertex(particle->Vx() - vertex[0], 
        	       particle->Vy() - vertex[1],
        	       particle->Vz() - vertex[2]);
      if (dVertex.Mag() > 0.0001) continue;
      if (TMath::Abs(particle->GetPdgCode()) == partCode[0]) nGene++;
      if (TMath::Abs(particle->GetPdgCode()) == partCode[2]) nGenpi++;
    }  // end loop over particles
    printf("Gen Prim. el= %d, Gen Prim. pi= %d in Event No.= %d in File = %s \n", nGene, nGenpi, iEvent, FileDir[imom]);
    
    // get the event summary data
    tree->GetEvent(iEvent);
    if (!esd) {
      Error("CheckESD", "no ESD object found for event %d", iEvent);
      return kFALSE;
    }

    ///////////////////////////////////////////
    //    AliTRDpidESD::MakePID(esd);
    //    AliESDpid::MakePID(esd);
    ///////////////////////////////////////////

    // Initializing number of electrons and pions
    Int_t nne=0, nnpi=0;    
    // loop over tracks
    Int_t nTracks = esd->GetNumberOfTracks();
    printf("Found %d tracks. \n",nTracks);

    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);
      // select tracks of selected particles
      Int_t label = TMath::Abs(track->GetLabel());
      if (label > stack->GetNtrack()) continue;     // background
      TParticle* particle = stack->Particle(label);
      //      if (TMath::Abs(particle->Eta()) > 0.3) continue;
      TVector3 dVertex(particle->Vx() - vertex[0], 
        	       particle->Vy() - vertex[1],
        	       particle->Vz() - vertex[2]);
      if (dVertex.Mag() > 0.0001) continue;
      if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;
      if (track->GetConstrainedChi2() > 1e9) continue;
      Double_t p[3];
      track->GetConstrainedPxPyPz(p);
      TVector3 pTrack(p);
      //      Double_t mom = track->GetP();

      // PID/////////////////////////////////////////////////////
      if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) continue;
      Int_t iGen = 5;
      for (Int_t i = 0; i < 5; i++) {
        if (TMath::Abs(particle->GetPdgCode()) == partCode[i]) iGen = i;
      }

      // Filling muon kan and proton//////////////////////////////
      Int_t Itot=0;
      if(track->GetTRDsignal()) Itot=1;
      if(!Itot) continue;

      // dE/dx and TEST TRD
      Double_t TRDsignal = track->GetTRDsignal();
      //      if (iGen == 0) h1dEdxEL[imom]->Fill(TRDsignal); 
      //      if (iGen == 2) h1dEdxPI[imom]->Fill(TRDsignal); 
      //      if (iGen == 1) h1dEdxMU[imom]->Fill(TRDsignal);
      //      if (iGen == 3) h1dEdxKA[imom]->Fill(TRDsignal);
      //      if (iGen == 4) h1dEdxPR[imom]->Fill(TRDsignal);

      //Filling elcron and pions////////////////////////////////
      const Int_t kNPlane = 6;
      // Track Selection 
      Double_t Ise=1.;
      for (Int_t iPlane = 0; iPlane < kNPlane; iPlane++) {
        Ise*=track->GetTRDsignals(iPlane);
      }
      if(Ise==0) continue;
      
      for (Int_t iPlane = 0; iPlane < kNPlane; iPlane++) {
        Double_t TRDsignal = track->GetTRDsignals(iPlane);
        if (iGen == 0) h1dEdxEL[imom]->Fill(TRDsignal);
        if (iGen == 2) h1dEdxPI[imom]->Fill(TRDsignal);
        if (iGen == 1) h1dEdxMU[imom]->Fill(TRDsignal);
        if (iGen == 3) h1dEdxKA[imom]->Fill(TRDsignal);
        if (iGen == 4) h1dEdxPR[imom]->Fill(TRDsignal);

        if (iGen == 0) h1MaxTimBinEL[imom]->Fill(track->GetTRDTimBin(iPlane));
        if (iGen == 2) h1MaxTimBinPI[imom]->Fill(track->GetTRDTimBin(iPlane));
      }
      
      //      Int_t TRDncls = track->GetTRDncls();
      if (iGen == 0) nne++;
      if (iGen == 2) nnpi++;
    } // end of loop over tracks
    printf("Electron Tracks = %d, Pion Tracks =  %d \n", nne, nnpi);
  }// end of loop over events
  delete esd;
  esdFile->Close();
  delete esdFile;
  
  runLoader->UnloadHeader();
  runLoader->UnloadKinematics();
  delete runLoader;
  } // Loop over files
  
    
  Char_t fileprint[200];
  
  // draw the histograms if not in batch mode

  new TCanvas;
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  //  gStyle->SetStatColor(0);
  h1dEdxPI[0]->DrawCopy();
  h1dEdxEL[0]->DrawCopy("SAME");
  //  h1dEdxMU[0]->DrawCopy("SAME");
  //  h1dEdxKA[0]->DrawCopy("SAME");
  //  h1dEdxPR[0]->DrawCopy("SAME");

  new TCanvas;
  gPad->SetBorderMode(0);
  gPad->SetFillColor(0);
  h1MaxTimBinEL[0]->DrawCopy();
  h1MaxTimBinPI[0]->DrawCopy("SAME");


  // write the output histograms to a file
  
  Char_t fileresponse[200];
  sprintf(fileresponse,"TRDdEdxHistogramsV%d.root",oFile);
  
  TFile* outputFile = TFile::Open(fileresponse, "recreate");
  if (!outputFile || !outputFile->IsOpen()) {
    Error("CheckESD", "opening output file check.root failed");
    return kFALSE;
  }

  for (Int_t imom = 0; imom < kNMom; imom++) {
    h1dEdxPI[imom]->Write();
    h1dEdxEL[imom]->Write();
    h1dEdxMU[imom]->Write();
    h1dEdxKA[imom]->Write();
    h1dEdxPR[imom]->Write();
    h1MaxTimBinEL[imom]->Write();
    h1MaxTimBinPI[imom]->Write();
  }

  delete [] h1dEdxPI;  
  delete [] h1dEdxEL;
  delete [] h1dEdxMU;
  delete [] h1dEdxKA;
  delete [] h1dEdxPR;
  delete [] h1MaxTimBinEL;
  delete [] h1MaxTimBinPI;
  delete outputFile;
  
  // result of check
  Info("GetPIDsignals", "GetPIDsignals was successfull");
  return kTRUE;
  
} // end of main program

