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

////////////////////////////////////////////////////////////
//                                                        //
//  Task for Non-HFE reconstruction efficiency            //
//  Non-Photonic Electron identified with Invariant mass  //
//                                                        //
//  Author: Deepa Thomas (University of Texas at Austin)  //
//          Vivek Kumar Singh (VECC)                      //
////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisUtils.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliVertexerTracks.h"


#include "AliAnalysisTaskIPResolBeautyppCal.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"

#include "AliGenHijingEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliCentrality.h"
#include "AliMagF.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

ClassImp(AliAnalysisTaskIPResolBeautyppCal)
  //________________________________________________________________________
  AliAnalysisTaskIPResolBeautyppCal::AliAnalysisTaskIPResolBeautyppCal(const char *name)
: AliAnalysisTaskSE(name),
  fVevent(0),
  fAOD(0),
  fpVtx(0),
  fpidResponse(0),
  fMultSelection(0),
  fMCHeader(0),
  fMCArray(0),
  fCentrality(-1),
  fCentralityMin(0),
  fCentralityMax(20),
  fMultiplicity(-1),
  fNTotMCpart(0),
  fNpureMC(0),
  fNembMCpi0(0),
  fNembMCeta(0),
  fTPCnSigma(-999.0),
  ftype(-1),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fCentralityNoPass(0),
  fCentralityPass(0),
  fMultiplicityNoPass(0),
  fMultiplicityPass(0),
  fCentMultiplicityNoPass(0),
  fCentMultiplicityPass(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCnsig(0),
  fRecalIP(kTRUE),
fTrkPt_Ele(0),
fTrkPt_HFEle(0),
fTrkPt_NHFEle(0),
fTrkPt_GammaE(0),
fTrkPt_DalitzE(0),
fImpParSprs_All(0),
fImpParSprs_AllE(0),
fImpParSprs_HFEle(0),
fImpParSprs_NHFEle(0),
fImpParSprs_GammaE(0),
fImpParSprs_DalitzE(0)
{
  //Named constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskIPResolBeautyppCal::AliAnalysisTaskIPResolBeautyppCal()
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskIPResolBeautyppCal"),
fVevent(0),
fAOD(0),
fpVtx(0),
fpidResponse(0),
fMultSelection(0),
fMCHeader(0),
fMCArray(0),
fCentrality(-1),
fCentralityMin(0),
fCentralityMax(20),
fMultiplicity(-1),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
fTPCnSigma(-999.0),
ftype(-1),
fOutputList(0),
fNevents(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fCentralityNoPass(0),
fCentralityPass(0),
fMultiplicityNoPass(0),
fMultiplicityPass(0),
fCentMultiplicityNoPass(0),
fCentMultiplicityPass(0),
fNegTrkIDPt(0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
fdEdx(0),
fTPCnsig(0),
fRecalIP(kTRUE),
fTrkPt_Ele(0),
fTrkPt_HFEle(0),
fTrkPt_NHFEle(0),
fTrkPt_GammaE(0),
fTrkPt_DalitzE(0),
fImpParSprs_All(0),
fImpParSprs_AllE(0),
fImpParSprs_HFEle(0),
fImpParSprs_NHFEle(0),
fImpParSprs_GammaE(0),
fImpParSprs_DalitzE(0)
{
  //Default constructor
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(3, TTree::Class());
}
//_________________________________________
AliAnalysisTaskIPResolBeautyppCal::~AliAnalysisTaskIPResolBeautyppCal()
{
  //Destructor
  delete fOutputList;
}
//_________________________________________
void AliAnalysisTaskIPResolBeautyppCal::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  AliDebug(3, "Creating Output Objects");

  Double_t pi = TMath::Pi();
    
  ///////////////
  //Output list//
  ///////////////

  fOutputList = new TList();
  fOutputList->SetOwner();

  fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
  fOutputList->Add(fNevents);
  fNevents->GetYaxis()->SetTitle("counts");
  fNevents->GetXaxis()->SetBinLabel(1,"All");
  fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
  fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");

  fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

  fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
  fOutputList->Add(fVtxY);

  fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
  fOutputList->Add(fVtxX);

  fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass;centrality;counts", 101, -1, 100);
  fOutputList->Add(fCentralityPass);

  fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass;centrality;counts", 101, -1, 100);
  fOutputList->Add(fCentralityNoPass);

  fMultiplicityPass = new TH1F("fMultiplicityPass", "Multiplicity Pass;multiplicity;counts", 2500, 0, 25000);
  fOutputList->Add(fMultiplicityPass);

  fMultiplicityNoPass = new TH1F("fMultiplicityNoPass", "Multiplicity No Pass;multiplicity;counts", 2500, 0, 25000);
  fOutputList->Add(fMultiplicityNoPass);

  fCentMultiplicityNoPass = new TH2F("fCentMultiplicityNoPass", "Multiplicity vs Centrality, No Pass;multiplicity;centrality", 2500, 0, 25000, 101, -1, 100);
  fOutputList->Add(fCentMultiplicityNoPass);

  fCentMultiplicityPass = new TH2F("fCentMultiplicityPass", "Multiplicity vs Centrality, Pass;multiplicity;centrality", 2500, 0, 25000, 101, -1, 100);
  fOutputList->Add(fCentMultiplicityPass);

  fNegTrkIDPt = new TH1F("fNegTrkIDPt", "p_{T} distribution of tracks with negative track id;p_{T} (GeV/c);counts", 500, 0.0, 50.0);
  fOutputList->Add(fNegTrkIDPt);

  fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",700,0.1,35);
  fOutputList->Add(fTrkPt);

  fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fTrketa);

  fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,2*pi);
  fOutputList->Add(fTrkphi);

  fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
  fOutputList->Add(fdEdx);

  fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fTPCnsig);
    
    fTrkPt_Ele = new TH1F("fTrkPt_Ele","p_{T} distribution of all electrons;p_{T} (GeV/c);counts",700,0.1,35);
    fOutputList->Add(fTrkPt_Ele);
    
    fTrkPt_HFEle = new TH1F("fTrkPt_HFEle","p_{T} distribution of HF electrons;p_{T} (GeV/c);counts",700,0.1,35);
    fOutputList->Add(fTrkPt_HFEle);
    
    fTrkPt_NHFEle = new TH1F("fTrkPt_NHFEle","p_{T} distribution of Non-HF electrons;p_{T} (GeV/c);counts",700,0.1,35);
    fOutputList->Add(fTrkPt_NHFEle);
    
    fTrkPt_GammaE = new TH1F("fTrkPt_GammaE","p_{T} distribution of gamma electrons;p_{T} (GeV/c);counts",700,0.1,35);
    fOutputList->Add(fTrkPt_GammaE);
    
    fTrkPt_DalitzE = new TH1F("fTrkPt_DalitzE","p_{T} distribution of Dalitz electrons;p_{T} (GeV/c);counts",700,0.1,35);
    fOutputList->Add(fTrkPt_DalitzE);
    
    Int_t nbins[5] =  {3000, 70, 4, 2, 10};
    Double_t limitLow[5] = {-1500., 0.1, 0., 0., -10.};
    Double_t limitUp[5] =  {1500., 35., 4., 2., +10.};
    TString axTitle[7]={"imp. par. (#mum)", "#it{p}_{T} (GeV/c)", "#phi", "#eta", "{z}_{vtx} (cm)"};
    
    fImpParSprs_All=new THnSparseF("fImpParSprs_All","Sparse for all charged particles",5,nbins,limitLow,limitUp);
    fOutputList->Add(fImpParSprs_All);
    
    fImpParSprs_AllE=new THnSparseF("fImpParSprs_AllE","Sparse for all electrons",5,nbins,limitLow,limitUp);
    fOutputList->Add(fImpParSprs_AllE);
    
    fImpParSprs_HFEle=new THnSparseF("fImpParSprs_HFEle","Sparse for HF electrons",5,nbins,limitLow,limitUp);
    fOutputList->Add(fImpParSprs_HFEle);
    
    fImpParSprs_NHFEle=new THnSparseF("fImpParSprs_NHFEle","Sparse for Non-HF electrons",5,nbins,limitLow,limitUp);
    fOutputList->Add(fImpParSprs_NHFEle);
    
    fImpParSprs_GammaE=new THnSparseF("fImpParSprs_GammaE","Sparse for gamma electrons",5,nbins,limitLow,limitUp);
    fOutputList->Add(fImpParSprs_GammaE);
    
    fImpParSprs_DalitzE=new THnSparseF("fImpParSprs_DalitzE","Sparse for Dalitz electrons",5,nbins,limitLow,limitUp);
    fOutputList->Add(fImpParSprs_DalitzE);
    
    for(Int_t iax=0; iax<5; iax++) {
        fImpParSprs_All->GetAxis(iax)->SetTitle(axTitle[iax].Data());
        fImpParSprs_AllE->GetAxis(iax)->SetTitle(axTitle[iax].Data());
        fImpParSprs_HFEle->GetAxis(iax)->SetTitle(axTitle[iax].Data());
        fImpParSprs_NHFEle->GetAxis(iax)->SetTitle(axTitle[iax].Data());
        fImpParSprs_GammaE->GetAxis(iax)->SetTitle(axTitle[iax].Data());
        fImpParSprs_DalitzE->GetAxis(iax)->SetTitle(axTitle[iax].Data());
    }
    
  PostData(1,fOutputList);
}
//_________________________________________
void AliAnalysisTaskIPResolBeautyppCal::UserExec(Option_t*)
{
  // Main loop
  // Called for each event
  // Post output data.

  UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent) {
    printf("ERROR: fVEvent not available\n");
    return;
  }
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fpVtx = fVevent->GetPrimaryVertex();

  ///////////////////
  //PID initialised//
  ///////////////////
  fpidResponse = fInputHandler->GetPIDResponse();

  ////////////////////
  //event selection///
  ////////////////////
  if(!PassEventSelect(fVevent)) return;

  /////////////////
  // Centrality ///
  /////////////////
/*  Bool_t pass = kFALSE;
  if(fCentralityMin > -0.5){
    CheckCentrality(fAOD,pass);
    if(!pass)return;
  }
*/
  //////////////////
  //Get MC headers//
  //////////////////

  fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!fMCArray){
    AliError("Array of MC particles not found");
    return;
  }
  fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  if (!fMCHeader) {
    AliError("Could not find MC Header in AOD");
    return;
  }

  ////////////////////////////////
  //Get number of Gen particles //
  ////////////////////////////////
  GetNMCPartProduced();

  ///////////////
  //Track loop///
  ///////////////

  Int_t ntracks = fVevent->GetNumberOfTracks();
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVParticle* Vtrack = 0x0;
    Vtrack  = fVevent->GetTrack(iTracks);
    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);

    ////////////////////
    //Apply track cuts//
    ////////////////////
    if(!PassTrackCuts(atrack)) continue;

    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
    TrkPhi = track->Phi();
    TrkPt = track->Pt();
    TrkEta = track->Eta();
    TrkP = track->P();
    Double_t zvtx = fpVtx->GetZ();

      Double_t beampiperadius=3.;
      Double_t d0z0[2], cov[3];
      atrack->PropagateToDCA(fpVtx, fVevent->GetMagneticField(), beampiperadius, d0z0, cov);
      
      Int_t phiBin = PhiBin(TrkPhi);
      
      Int_t etaBin=0.0;
      if(TrkEta < 0) etaBin = 0.;
      if(TrkEta > 0) etaBin = 1.;

      Double_t fsparse[5];
      fsparse[0] = 10000.*d0z0[0];
      fsparse[1] = TrkPt;
      fsparse[2] = phiBin;
      fsparse[3] = etaBin;
      fsparse[4] = zvtx;

      fImpParSprs_All->Fill(fsparse);
      
      Int_t iTrklabel = TMath::Abs(track->GetLabel());
      if(iTrklabel == 0) continue;
      AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
      //if(!MCPart->IsPhysicalPrimary()) continue;
      
      if(TMath::Abs(MCPart->GetPdgCode())!=11) continue;
      
      fTrkPt_Ele->Fill(TrkPt);
      fImpParSprs_AllE->Fill(fsparse);

      Int_t iMCmom=-999, MomPDG = -999;
      iMCmom = MCPart->GetMother();
      AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
      MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
      
      //cout << "ele mom : " << MomPDG <<endl;

      if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
          fTrkPt_NHFEle->Fill(TrkPt);
          fImpParSprs_NHFEle->Fill(fsparse);
          
          if(MomPDG == 22) {
              fTrkPt_GammaE->Fill(TrkPt);
              fImpParSprs_GammaE->Fill(fsparse);
          }
          if((MomPDG == 111) || (MomPDG == 221))
          {
              fTrkPt_DalitzE->Fill(TrkPt);
              fImpParSprs_DalitzE->Fill(fsparse);
          }
      }
      
      if((MomPDG > 400) && (MomPDG < 600)){
          fTrkPt_HFEle->Fill(TrkPt);
          fImpParSprs_HFEle->Fill(fsparse);
      }

      ftype = -1;
      Bool_t fFromHijing = kTRUE;
      Double_t MomPt =-999;
      //Bool_t fNonHFE = IsNonHFE(MCPart, fFromHijing, ftype, iMCmom, MomPDG, MomPt);
      
  }//track loop

  PostData(1, fOutputList);
}
//___________________________________________
Bool_t  AliAnalysisTaskIPResolBeautyppCal::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
{
  //Is electron from pi0, eta and gamma

  iMCmom = MCPart->GetMother();
  AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
  MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
  MomPt = MCPartMom->Pt();

  if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
    if(iMCmom >= fNpureMC)fFromHijing = kFALSE;
    type = GetPi0EtaType(MCPartMom);
    return kTRUE;
  }
  else return kFALSE;
}
//_________________________________________
Bool_t AliAnalysisTaskIPResolBeautyppCal::GetNMCPartProduced()
{
  //Get number of MC particles produced by generators.

  TList *lh = fMCHeader->GetCocktailHeaders();
  fNTotMCpart = 0;
  fNembMCpi0 = 0;
  fNembMCeta = 0;
  fNpureMC = 0;
  TString MCgen;
  TString embpi0("pi");
  TString embeta("eta");

  if(!lh){
    AliError("no MC header");
    return (0);
  }

  for(int igene=0; igene<lh->GetEntries(); igene++)
  {
    AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
    if(!gh) continue;

    MCgen =  gh->GetName();
    //   cout << "Gen name, N produced = " << gh->GetName() << ", " << gh->NProduced() << endl;
    if(igene==0) fNpureMC = gh->NProduced();  // generated by HIJING

    //   if(MCgen.Contains(embpi0))cout << MCgen << endl;
    //   if(MCgen.Contains(embeta))cout << MCgen << endl;

    if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
    if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
    fNTotMCpart += gh->NProduced();
  }
  //  cout << "fNpureMC, fNembMCpi0, fNembMCeta, fNTotMCpart : " <<fNpureMC << ", " << fNembMCpi0 << ", " << fNembMCeta << ", " << fNTotMCpart << endl;

  return kTRUE;
}
//_________________________________________
Int_t AliAnalysisTaskIPResolBeautyppCal::GetPi0EtaType(AliAODMCParticle *part)
{
  // Return the type of particle

  // IsPrimary
  Bool_t primMC = part->IsPrimary();
  if(!primMC) return kNotIsPrimary;

  // Mother
  Int_t motherlabel = part->GetMother();
  if(motherlabel<0) return kNoMother;

  else {
    AliAODMCParticle *mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
    Int_t motherpdg = TMath::Abs(mother->GetPdgCode());

    if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;

    if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
    if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
      
    return kNoFeedDown;
  }
}
//_________________________________________
Int_t AliAnalysisTaskIPResolBeautyppCal::GetPrimary(Int_t id)
{
  // Return the label of the primary that generated the given track

  int current, parent;
  parent=id;
  while (1) {
    current=parent;
    AliAODMCParticle *part = (AliAODMCParticle*)fMCArray->At(current);
    parent = part->GetMother();
    if(parent<0) return current;
  }
}
//_________________________________________
Bool_t AliAnalysisTaskIPResolBeautyppCal::PassEventSelect(AliVEvent *fVevent)
{
  //event selection cuts

  Int_t ntracks = -999;
  ntracks = fVevent->GetNumberOfTracks();
  if(ntracks < 1) printf("There are %d tracks in this event\n",ntracks);

  fNevents->Fill(0); //all events
  Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;

  Double_t NcontV = fpVtx->GetNContributors();
  if(NcontV<2)return kFALSE;

  Bool_t isPileupfromSPDmulbins=fAOD->IsPileupFromSPDInMultBins(); //This function checks if there was a pile up reconstructed with SPD
  if(isPileupfromSPDmulbins)return kFALSE;

  Int_t minContributors=5;    //minimum contributors to the pilup vertices, multi-vertex
  Double_t minChi2=5.; 
  Double_t minWeiZDiff=15;   //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
  Bool_t checkPlpFromDifferentBC=kFALSE;
  
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(minContributors); //Multi Vertex pileup selection
  utils.SetMaxPlpChi2MV(minChi2);   //max value of Chi2perNDF of the pileup vertex, multi-vertex
  utils.SetMinWDistMV(minWeiZDiff);
  utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC); //SPD Pileup slection
  Bool_t isPileupFromMV = utils.IsPileUpMV(fAOD);      //check for multi-vertexer pile-up
  
  if(isPileupFromMV) return kFALSE;
  fNevents->Fill(1); //events after vettex cuts

  Zvertex = fpVtx->GetZ();
  Yvertex = fpVtx->GetY();
  Xvertex = fpVtx->GetX();
  fVtxZ->Fill(Zvertex);
  fVtxX->Fill(Xvertex);
  fVtxY->Fill(Yvertex);
  if(TMath::Abs(Zvertex)>10.0) return kFALSE;
  fNevents->Fill(2); //events after z vtx cut

  return kTRUE;
}
//_________________________________________
/*void AliAnalysisTaskIPResolBeautyppCal::CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass)
{
  //check centrality, Run 2

  if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
  if(!fMultSelection) {
    //If you get this warning (and lPercentiles 300) check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
    fCentrality = fMultSelection->GetMultiplicityPercentile("V0M", false);
  }

  AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  if(!header) AliFatal("Not a standard AOD");
  fMultiplicity = header->GetRefMultiplicity();

  if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
  {
    fCentralityNoPass->Fill(fCentrality);
    fMultiplicityNoPass->Fill(fMultiplicity);
    fCentMultiplicityNoPass->Fill(fMultiplicity,fCentrality);
    //  cout << "--------------Fill no pass-------------------------"<<endl;
    centralitypass = kFALSE;
  }else
  {
    fCentralityPass->Fill(fCentrality);
    fMultiplicityPass->Fill(fMultiplicity);
    fCentMultiplicityPass->Fill(fMultiplicity,fCentrality);
    //  cout << "--------------Fill pass-------------------------"<<endl;
    centralitypass = kTRUE;
  }
}
*/
//___________________________________________
Bool_t AliAnalysisTaskIPResolBeautyppCal::PassTrackCuts(AliAODTrack *atrack)
{
  //apply track cuts

  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 1.0, DCAzCut = 2.0;
  Double_t dEdx =-999;
  Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    
    if(TMath::Abs(atrack->Eta()) > 0.9) return kFALSE;

  //kink daughters
  Int_t numberofvertices = 100;
  numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
  for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      listofmotherkink[numberofmotherkink] = idmother;
      numberofmotherkink++;
    }
  }
  if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; //mimimum cuts
  //reject kink
  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
    if(atrack->GetID() == listofmotherkink[kinkmother]) {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;

    
  //other cuts
  if(atrack->Pt()  < 0.5 ) return kFALSE;
  if (TMath::Abs(atrack->Eta()) > 0.7 ) return kFALSE;                     //Eta cut
  if(atrack->GetITSNcls() < 3) return kFALSE;

  Double_t TPCNClsF =  atrack->GetTPCNclsF();
  Double_t TPCNCrossedRows = atrack->GetTPCNCrossedRows();
  Double_t RatioCrossedRowsOverFindableClusters = -999;

  if( TPCNCrossedRows < 70) return kFALSE; //TPC N clusters
  if(TPCNClsF > 0){ RatioCrossedRowsOverFindableClusters = TPCNCrossedRows/TPCNClsF; }
  if(RatioCrossedRowsOverFindableClusters < 0.8) return kFALSE;

  if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;

  if(fRecalIP) RecalImpactParam(atrack, d0z0, cov);  

  if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
  if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;

  Double_t chi2ndf = atrack->Chi2perNDF();
  if(chi2ndf>4.0) return kFALSE;


  ////////////////////
  //Track properties//
  ////////////////////
  dEdx = atrack->GetTPCsignal();
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(atrack, AliPID::kElectron);
  TrkPhi = atrack->Phi();
  TrkPt = atrack->Pt();
  TrkEta = atrack->Eta();
  TrkP = atrack->P();

  if(atrack->GetID()<0) fNegTrkIDPt->Fill(TrkPt);
  fTrkPt->Fill(TrkPt);
  fTrketa->Fill(TrkEta);
  fTrkphi->Fill(TrkPhi);
  fdEdx->Fill(TrkP,dEdx);
  fTPCnsig->Fill(TrkP,fTPCnSigma);

  return kTRUE;
}
//----------------------------------------------------------------------------
void AliAnalysisTaskIPResolBeautyppCal::RecalImpactParam(const AliAODTrack * const track, Double_t dcaD[2], Double_t covD[3])
{
    //Recalculate impact parameter by recalculating primary vertex
    
    const Double_t kBeampiperadius=3.0;
    Bool_t isRecalcVertex = kFALSE;

    AliAODVertex *vtxAODSkip  = fAOD->GetPrimaryVertex();
    if(!vtxAODSkip) return;
    
    Double_t fMagField = fAOD->GetMagneticField();

    const AliAODTrack *tmptrack = dynamic_cast<const AliAODTrack *>(track);
    if(tmptrack){
        if(vtxAODSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
            
            vtxAODSkip = RemoveDaughtersFromPrimaryVtx(track);
            isRecalcVertex = kTRUE;
        }
        
        if(vtxAODSkip){
            AliAODTrack aodtrack(*tmptrack);
            AliExternalTrackParam etp;
            etp.CopyFromVTrack(&aodtrack);
            
            etp.PropagateToDCA(vtxAODSkip, fMagField, kBeampiperadius, dcaD, covD);
            
            if(isRecalcVertex) delete vtxAODSkip;
        }
    }
}
//________________________________________________________________________
AliAODVertex* AliAnalysisTaskIPResolBeautyppCal::RemoveDaughtersFromPrimaryVtx(const AliAODTrack * const track)
{
    // This method returns a primary vertex without the daughter tracks of the
    // candidate and it recalculates the impact parameters and errors for AOD tracks.
    
    AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
    if(!vtxAOD) return 0;
    TString title=vtxAOD->GetTitle();
    if(!title.Contains("VertexerTracks")) return 0;

    AliVertexerTracks vertexer(fAOD->GetMagneticField());
    
    vertexer.SetITSMode();
    vertexer.SetMinClusters(3);
    vertexer.SetConstraintOff();
    
    if(title.Contains("WithConstraint")) {
        Float_t diamondcovxy[3];
        fAOD->GetDiamondCovXY(diamondcovxy);
        Double_t pos[3]={fAOD->GetDiamondX(),fAOD->GetDiamondY(),0.};
        Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
        AliESDVertex diamond(pos,cov,1.,1);
        vertexer.SetVtxStart(&diamond);
    }
    Int_t skipped[2]; for(Int_t i=0;i<2;i++) skipped[i]=-1;
    Int_t id = (Int_t)track->GetID();
    if(!(id<0)) skipped[0] = id;
    
    vertexer.SetSkipTracks(1,skipped);
    AliESDVertex *vtxESDNew = vertexer.FindPrimaryVertex(fAOD);
    
    if(!vtxESDNew) return 0;
    if(vtxESDNew->GetNContributors()<=0) {
        delete vtxESDNew; vtxESDNew=NULL;
        return 0;
    }
    
    // convert to AliAODVertex
    Double_t pos[3],cov[6],chi2perNDF;
    vtxESDNew->GetXYZ(pos); // position
    vtxESDNew->GetCovMatrix(cov); //covariance matrix
    chi2perNDF = vtxESDNew->GetChi2toNDF();
    delete vtxESDNew; vtxESDNew=NULL;
    
    AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);
    
    return vtxAODNew;
}

//___________________________________________
Int_t AliAnalysisTaskIPResolBeautyppCal::PhiBin(Double_t phi) const {
    //phi bin
    
    Double_t pi=TMath::Pi();
    if(phi>2.*pi || phi<0.) return -1;
    if((phi<=(pi/4.)) || (phi>7.*(pi/4.))) return 0;
    if((phi>(pi/4.)) && (phi<=3.*(pi/4.))) return 1;
    if((phi>3.*(pi/4.)) && (phi<=5.*(pi/4.))) return 2;
    if((phi>(5.*pi/4.)) && (phi<=7.*(pi/4.))) return 3;
    return -1;
}
//___________________________________________
void AliAnalysisTaskIPResolBeautyppCal::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}



