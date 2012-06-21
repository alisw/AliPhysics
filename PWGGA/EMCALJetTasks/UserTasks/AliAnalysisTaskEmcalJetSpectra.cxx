// $Id$
//
// Jet spectrum task.
//
// Author: R.Reed, M.Connors

#include "AliAnalysisTaskEmcalJetSpectra.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEmcalJet.h"
#include "AliVCluster.h"

ClassImp(AliAnalysisTaskEmcalJetSpectra)

//________________________________________________________________________
AliAnalysisTaskEmcalJetSpectra::AliAnalysisTaskEmcalJetSpectra() : 
  AliAnalysisTaskSE(), 
  fTracksName("tracks"),
  fJetsName("jets"),
  fClustersName("clusters"),
  fRhos1Name("fArrRhos"),
  fRhos2Name("fArrRhos"),
  fRhos3Name("fArrRhos"),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fESD(0), 
  fOutputList(0), 
  fHistCentrality(0),  
  fHistDeltaRho12vsCent(0), 
  fHistDeltaRho13vsCent(0), 
  fHistDeltaRho23vsCent(0), 
  fHistDeltaJetPt12vsCent(0), 
  fHistDeltaJetPt13vsCent(0), 
  fHistDeltaJetPt23vsCent(0), 
  fHistRho1vsCent(0), 
  fHistRho2vsCent(0), 
  fHistRho3vsCent(0)
{
  // Default constructor.

  for (Int_t i = 0;i<6;++i) {
    fHistRawJetPt[i]       = 0;
    fHistAreavsRawPt[i]    = 0;
    for (Int_t j = 0;j<4;++j) {
      fHistNEFvsPt[i][j]   = 0;
       fHistZvsPt[i][j]     = 0;
       fHistZchvsPt[i][j]   = 0;
       fHistZemvsPt[i][j]   = 0;
       fHistJetPt[i][j]     = 0;
       fHistNconsvsPt[i][j] = 0;
       fHistJetPt5[i][j]    = 0;
       fHistJetPt6[i][j]    = 0;
       fHistJetPt7[i][j]    = 0;
       fHistJetPt8[i][j]    = 0;
    }
  }
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetSpectra::AliAnalysisTaskEmcalJetSpectra(const char *name) :
  AliAnalysisTaskSE(name), 
  fTracksName("tracks"),
  fJetsName("jets"),
  fClustersName("clusters"),
  fRhos1Name("fArrRhos"),
  fRhos2Name("fArrRhos"),
  fRhos3Name("fArrRhos"),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fESD(0), 
  fOutputList(0), 
  fHistCentrality(0),  
  fHistDeltaRho12vsCent(0), 
  fHistDeltaRho13vsCent(0), 
  fHistDeltaRho23vsCent(0), 
  fHistDeltaJetPt12vsCent(0), 
  fHistDeltaJetPt13vsCent(0), 
  fHistDeltaJetPt23vsCent(0), 
  fHistRho1vsCent(0), 
  fHistRho2vsCent(0), 
  fHistRho3vsCent(0)
{
  // Constructor

  for (Int_t i = 0;i<6;++i) {
    fHistRawJetPt[i]       = 0;
    fHistAreavsRawPt[i]    = 0;
    for (Int_t j = 0;j<4;++j) {
      fHistNEFvsPt[i][j]   = 0;
      fHistZvsPt[i][j]     = 0;
      fHistZchvsPt[i][j]   = 0;
      fHistZemvsPt[i][j]   = 0;
      fHistJetPt[i][j]     = 0;
      fHistNconsvsPt[i][j] = 0;
      fHistJetPt5[i][j]    = 0;
      fHistJetPt6[i][j]    = 0;
      fHistJetPt7[i][j]    = 0;
      fHistJetPt8[i][j]    = 0;
    }
  }
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSpectra::UserCreateOutputObjects()
{
  // Called once at the beginning of the analysis.

  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  fHistCentrality            = new TH1F("Centrality",             "Centrality",            101, -1,  100);
  fHistDeltaRho12vsCent      = new TH2F("DeltaRho12vsCent",       "DeltaRho12vsCent",      100, 0.0, 100.0, 500,-250,250);
  fHistDeltaRho13vsCent      = new TH2F("DeltaRho13vsCent",       "DeltaRho13vsCent",      100, 0.0, 100.0, 500,-250,250);
  fHistDeltaRho23vsCent      = new TH2F("DeltaRho23vsCent",       "DeltaRho23vsCent",      100, 0.0, 100.0, 500,-250,250);
  fHistDeltaJetPt12vsCent    = new TH2F("DeltaJetPt12vsCent",     "DeltaJetPt12vsCent",    100, 0.0, 100.0, 500,-250,250);
  fHistDeltaJetPt13vsCent    = new TH2F("DeltaJetPt13vsCent",     "DeltaJetPt13vsCent",    100, 0.0, 100.0, 500,-250,250);
  fHistDeltaJetPt23vsCent    = new TH2F("DeltaJetPt23vsCent",     "DeltaJetPt23vsCent",    100, 0.0, 100.0, 500,-250,250);
  fHistRho1vsCent            = new TH2F("Rho1vsCent",             "Rho1vsCent",            100, 0.0, 100.0, 500, 0, 500);
  fHistRho2vsCent            = new TH2F("Rho2vsCent",             "Rho2vsCent",            100, 0.0, 100.0, 500, 0, 500);
  fHistRho3vsCent            = new TH2F("Rho3vsCent",             "Rho3vsCent",            100, 0.0, 100.0, 500, 0, 500);

  for (Int_t i = 0;i<6;++i) {
    TString name00(Form("fHistRawJetPt_%i",i));
    fHistRawJetPt[i] = new TH1F(name00,name00,250,0,500);
    fOutputList->Add(fHistRawJetPt[i]);
    TString name01(Form("fHistAreavsRawPt_%i",i));
    fHistAreavsRawPt[i] = new TH2F(name01,name01,250,0,500,100,0,1);
    fOutputList->Add(fHistAreavsRawPt[i]);
    for (Int_t j = 0;j<4;j++) {
      TString name0(Form("fHistNEFvsPt_%i_%i",i,j));
      fHistNEFvsPt[i][j] = new TH2F(name0,name0,250,-250,250,200,0,2);
      fOutputList->Add(fHistNEFvsPt[i][j]);
      TString name1(Form("fHistZvsPt_%i_%i",i,j));
      fHistZvsPt[i][j] = new TH2F(name1,name1,250,-250,250,200,0,2);
      fOutputList->Add(fHistZvsPt[i][j]);
      TString name2(Form("fHistZchvsPt_%i_%i",i,j));
      fHistZchvsPt[i][j] = new TH2F(name2,name2,250,-250,250,200,0,2);
      fOutputList->Add(fHistZchvsPt[i][j]);
      TString name3(Form("fHistZemvsPt_%i_%i",i,j));
      fHistZemvsPt[i][j] = new TH2F(name3,name3,250,-250,250,200,0,2);
      fOutputList->Add(fHistZemvsPt[i][j]);
      TString name5(Form("fHistJetPt_%i_%i",i,j));
      fHistJetPt[i][j] = new TH1F(name5,name5,250,-250,250);
      fOutputList->Add(fHistJetPt[i][j]);
      TString name6(Form("fHistNconsvsPt_%i_%i",i,j));
      fHistNconsvsPt[i][j] = new TH2F(name6,name6,250,-250,250,500,0,500);
      fOutputList->Add(fHistNconsvsPt[i][j]);
      TString name7(Form("fHistJetPt5_%i_%i",i,j));
      fHistJetPt5[i][j] = new TH1F(name7,name7,250,-250,250);
      fOutputList->Add(fHistJetPt5[i][j]);
      TString name8(Form("fHistJetPt6_%i_%i",i,j));
      fHistJetPt6[i][j] = new TH1F(name8,name8,250,-250,250);
      fOutputList->Add(fHistJetPt6[i][j]);
      TString name9(Form("fHistJetPt7_%i_%i",i,j));
      fHistJetPt7[i][j] = new TH1F(name9,name9,250,-250,250);
      fOutputList->Add(fHistJetPt7[i][j]);
      TString name10(Form("fHistJetPt8_%i_%i",i,j));
      fHistJetPt8[i][j] = new TH1F(name10,name10,250,-250,250);
      fOutputList->Add(fHistJetPt8[i][j]);
    }
  }
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistDeltaRho12vsCent);
  fOutputList->Add(fHistDeltaRho13vsCent);
  fOutputList->Add(fHistDeltaRho23vsCent);
  fOutputList->Add(fHistDeltaJetPt12vsCent);
  fOutputList->Add(fHistDeltaJetPt13vsCent);
  fOutputList->Add(fHistDeltaJetPt23vsCent);
  fOutputList->Add(fHistRho1vsCent);
  fOutputList->Add(fHistRho2vsCent);
  fOutputList->Add(fHistRho3vsCent);
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetSpectra::GetCentBin(Double_t cent) const 
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0;
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSpectra::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // esd or aod mode
  Bool_t esdMode = kTRUE;
  if (dynamic_cast<AliAODEvent*>(InputEvent()))
    esdMode = kFALSE;

  if (esdMode) {
    // optimization in case autobranch loading is off
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    if (fTracksName == "Tracks")
      am->LoadBranch("Tracks");
  }

  // get centrality 
  Double_t fCent = -1; 
  TList *l = InputEvent()->GetList();
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  if (centrality)
    fCent = centrality->GetCentralityPercentile("V0M");
  else
    fCent=99; // probably pp data
  if (fCent<0) {
    AliError(Form("Centrality negative: %f", fCent));
    return;
  }

  Int_t centbin = GetCentBin(fCent);

  //don't want to analyze events 90-100
  if (centbin < 0)
    return;

  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;

  tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
    
  jets = dynamic_cast<TClonesArray*>(l->FindObject(fJetsName));
  if (!jets) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }

  Double_t rho1 = -1;
  TParameter<Double_t> *Rho1Param = dynamic_cast<TParameter<Double_t>*>(InputEvent()->FindListObject(fRhos1Name));
  if (Rho1Param)
    rho1 = Rho1Param->GetVal();

  Double_t rho2 = -1;
  TParameter<Double_t> *Rho2Param = dynamic_cast<TParameter<Double_t>*>(InputEvent()->FindListObject(fRhos2Name));
  if (Rho2Param)
    rho2 = Rho2Param->GetVal();

  Double_t rho3 = -1;
  TParameter<Double_t> *Rho3Param = dynamic_cast<TParameter<Double_t>*>(InputEvent()->FindListObject(fRhos3Name));
  if (Rho3Param)
    rho3 = Rho3Param->GetVal();

  if (rho1>0)
    fHistRho1vsCent->Fill(fCent,rho1);
  if (rho2>0)
    fHistRho2vsCent->Fill(fCent,rho2);
  if (rho3>0)
    fHistRho3vsCent->Fill(fCent,rho3);

  if (( rho1>0 ) && ( rho2>0 )) 
    fHistDeltaRho12vsCent->Fill(fCent,rho1-rho2);
  if (( rho1>0 ) && ( rho3>0 )) 
    fHistDeltaRho13vsCent->Fill(fCent,rho1-rho3);
  if (( rho2>0 ) && ( rho3>0 )) 
    fHistDeltaRho23vsCent->Fill(fCent,rho2-rho3);
  fHistCentrality->Fill(fCent);

  Double_t lJetPt   = -500;
  Double_t lJetPtun = -500;
  const Int_t Njets = jets->GetEntries();
  Int_t NjetAcc = 0;
 
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(iJets));
    if (!jet)
      continue; 
    if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
      continue;
    if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
      continue;
    fHistAreavsRawPt[centbin]->Fill(jet->Pt(),jet->Area());
    if (jet->Area()<fAreacut)
      continue;
    //prevents 0 area jets from sneaking by when area cut == 0
    if (jet->Area()==0)
      continue;
    fHistRawJetPt[centbin]->Fill(jet->Pt());
    
    NjetAcc++;
    Double_t Z; 
    Double_t jetPt1 = -500;
    Double_t jetPt2 = -500;
    Double_t jetPt3 = -500;
    if (rho1 > 0)
      jetPt1 = jet->Pt()-jet->Area()*rho1;
    if (rho2 > 0)
      jetPt2 = jet->Pt()-jet->Area()*rho2;
    if (rho3 > 0)
      jetPt3 = jet->Pt()-jet->Area()*rho3;
    if (( rho1>0 ) && ( rho2>0 )) 
      fHistDeltaJetPt12vsCent->Fill(fCent,jetPt1-jetPt2);
    if (( rho1>0 ) && ( rho3>0 )) 
      fHistDeltaJetPt13vsCent->Fill(fCent,jetPt1-jetPt3);
    if (( rho2>0 ) && ( rho3>0 )) 
      fHistDeltaJetPt23vsCent->Fill(fCent,jetPt2-jetPt3);

    if (lJetPt < jetPt2) {
      lJetPt= jetPt2;
      lJetPtun = jet->Pt();
    }

    if (jet->MaxTrackPt()>jet->MaxClusterPt()) {
      Z = jet->MaxTrackPt();
    }
    else{
      Z = jet->MaxClusterPt();
    }

    fHistNEFvsPt[centbin][0]->Fill(jet->Pt(),jet->NEF());
    fHistZvsPt[centbin][0]->Fill(jet->Pt(),Z/jet->Pt());
    fHistZchvsPt[centbin][0]->Fill(jet->Pt(),jet->MaxTrackPt()/jet->Pt());
    fHistZemvsPt[centbin][0]->Fill(jet->Pt(),jet->MaxClusterPt()/jet->Pt());
    fHistJetPt[centbin][0]->Fill(jet->Pt());
    fHistNconsvsPt[centbin][0]->Fill(jet->Pt(),jet->N());
    if ((jet->MaxTrackPt()>5) || (jet->MaxClusterPt()>5))
      fHistJetPt5[centbin][0]->Fill(jet->Pt());
    if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6))
      fHistJetPt6[centbin][0]->Fill(jet->Pt());
    if ((jet->MaxTrackPt()>7) || (jet->MaxClusterPt()>7))
      fHistJetPt7[centbin][0]->Fill(jet->Pt());
    if ((jet->MaxTrackPt()>8) || (jet->MaxClusterPt()>8))
      fHistJetPt8[centbin][0]->Fill(jet->Pt());

    fHistNEFvsPt[centbin][1]->Fill(jetPt1,jet->NEF());
    fHistZvsPt[centbin][1]->Fill(jetPt1,Z/jetPt1);
    fHistZchvsPt[centbin][1]->Fill(jetPt1,jet->MaxTrackPt()/jetPt1);
    fHistZemvsPt[centbin][1]->Fill(jetPt1,jet->MaxClusterPt()/jetPt1);
    fHistJetPt[centbin][1]->Fill(jetPt1);
    fHistNconsvsPt[centbin][1]->Fill(jetPt1,jet->N());
    if ((jet->MaxTrackPt()>5) || (jet->MaxClusterPt()>5))
      fHistJetPt5[centbin][1]->Fill(jetPt1);
    if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6))
      fHistJetPt6[centbin][1]->Fill(jetPt1);
    if ((jet->MaxTrackPt()>7) || (jet->MaxClusterPt()>7))
      fHistJetPt7[centbin][1]->Fill(jetPt1);
    if ((jet->MaxTrackPt()>8) || (jet->MaxClusterPt()>8))
      fHistJetPt8[centbin][1]->Fill(jetPt1);

    fHistNEFvsPt[centbin][2]->Fill(jetPt2,jet->NEF());
    fHistZvsPt[centbin][2]->Fill(jetPt2,Z/jetPt2);
    fHistZchvsPt[centbin][2]->Fill(jetPt2,jet->MaxTrackPt()/jetPt2);
    fHistZemvsPt[centbin][2]->Fill(jetPt2,jet->MaxClusterPt()/jetPt2);
    fHistJetPt[centbin][2]->Fill(jetPt2);
    fHistNconsvsPt[centbin][2]->Fill(jetPt2,jet->N());
    if ((jet->MaxTrackPt()>5) || (jet->MaxClusterPt()>5))
      fHistJetPt5[centbin][2]->Fill(jetPt2);
    if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6))
      fHistJetPt6[centbin][2]->Fill(jetPt2);
    if ((jet->MaxTrackPt()>7) || (jet->MaxClusterPt()>7))
      fHistJetPt7[centbin][2]->Fill(jetPt2);
    if ((jet->MaxTrackPt()>8) || (jet->MaxClusterPt()>8))
      fHistJetPt8[centbin][2]->Fill(jetPt2);

    fHistNEFvsPt[centbin][3]->Fill(jetPt3,jet->NEF());
    fHistZvsPt[centbin][3]->Fill(jetPt3,Z/jetPt3);
    fHistZchvsPt[centbin][3]->Fill(jetPt3,jet->MaxTrackPt()/jetPt3);
    fHistZemvsPt[centbin][3]->Fill(jetPt3,jet->MaxClusterPt()/jetPt3);
    fHistJetPt[centbin][3]->Fill(jetPt3);
    fHistNconsvsPt[centbin][3]->Fill(jetPt3,jet->N());
    if ((jet->MaxTrackPt()>5) || (jet->MaxClusterPt()>5))
      fHistJetPt5[centbin][3]->Fill(jetPt3);
    if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6))
      fHistJetPt6[centbin][3]->Fill(jetPt3);
    if ((jet->MaxTrackPt()>7) || (jet->MaxClusterPt()>7))
      fHistJetPt7[centbin][3]->Fill(jetPt3);
    if ((jet->MaxTrackPt()>8) || (jet->MaxClusterPt()>8))
      fHistJetPt8[centbin][3]->Fill(jetPt3);
  }

  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSpectra::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}


