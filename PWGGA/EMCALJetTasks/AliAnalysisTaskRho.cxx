// $Id: $
//
// Calculation of rho 
//
// Authors: R.Reed, S.Aiola


#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TF1.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliVCluster.h"

#include "AliAnalysisTaskRho.h"

ClassImp(AliAnalysisTaskRho)

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho() : 
  AliAnalysisTaskRhoBase(),
  fTracksName("tracks"),
  fJetsName("KtJets"),
  fClustersName("caloClusters"),
  fRhoScaledName(""),
  fPhiMin(0),
  fPhiMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fAreaCut(0),
  fScaleFunction(0),
  fOutputList(0),
  fHistCentrality(0),
  fHistJetPt(0),
  fHistJetArea(0),
  fHistRhovsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistNjetvsCent(0),
  fHistRhovsNtrack(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fRhoScaled(0),
  fNewRhoFunction(0)
{
  // Constructor
}
//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho(const char *name) :
  AliAnalysisTaskRhoBase(name),
  fTracksName("tracks"),
  fJetsName("KtJets"),
  fClustersName("caloClusters"),
  fRhoScaledName(""),
  fPhiMin(0),
  fPhiMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fAreaCut(0),
  fScaleFunction(0),
  fOutputList(0),
  fHistCentrality(0),
  fHistJetPt(0),
  fHistJetArea(0),
  fHistRhovsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistNjetvsCent(0),
  fHistRhovsNtrack(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fRhoScaled(0),
  fNewRhoFunction(0)
{
  // Constructor

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskRho::UserCreateOutputObjects()
{
  AliAnalysisTaskRhoBase::UserCreateOutputObjects();

  fRhoScaledName = fRhoName;
  fRhoScaledName += "_Scaled";
  fRhoScaled = new TParameter<Double_t>(fRhoScaledName, 0);

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  fHistCentrality             = new TH1F("Centrality",            "Centrality",            101, -1,  100);
  fHistRhovsCent              = new TH2F("RhovsCent",             "RhovsCent",             101, -1,  100,   500,  0,   500);
  fHistDeltaRhovsCent         = new TH2F("DeltaRhovsCent",        "DetlaRhovsCent",        101, -1,  100,   500, -250, 250);
  fHistDeltaRhoScalevsCent    = new TH2F("DeltaRhoScalevsCent",   "DeltaRhoScalevsCent",   101, -1,  100,   500, -250, 250);
  fHistJetPtvsCent            = new TH2F("JetPtvsCent",           "JetPtvsCent",           101, -1,  100,   200,  0,   500);
  fHistJetAreavsCent          = new TH2F("JetAreavsCent",         "JetAreavsCent",         101, -1,  100,   100,  0,   1.0);
  fHistNjetvsCent             = new TH2F("NjetvsCent",            "NjetvsCent",            101, -1,  100,   100,  0,   100);

  fHistRhovsNtrack            = new TH2F("RhovsNtrack",           "RhovsNtrack",           500,  0,  2500,  500,  0,   500);
  fHistDeltaRhovsNtrack       = new TH2F("DeltaRhovsNtrack",      "DeltaRhovsNtrack",      500,  0,  2500,  500, -250, 250);
  fHistDeltaRhoScalevsNtrack  = new TH2F("DeltaRhoScalevsNtrack", "DeltaRhoScalevsNtrack", 500,  0,  2500,  500, -250, 250);
  fHistJetPtvsNtrack          = new TH2F("JetPtvsNtrack",         "JetPtvsNtrack",         500,  0,  2500,  200,  0,   500);
  fHistJetAreavsNtrack        = new TH2F("JetAreavsNtrack",       "JetAreavsNtrack",       500,  0,  2500,  100,  0,   1.0);
  fHistNjetvsNtrack           = new TH2F("NjetvsNtrack",          "rNjetvsNtrack",         500,  0,  2500,  100,  0,   100);

  fHistJetPt                  = new TH1F("JetPt",                  "Jet Pt",               100,   0,    250);
  fHistJetArea                = new TH1F("JetArea",                "Jet Area",             100, 0.0,    1.0);

  fNewRhoFunction             = new TF1("rfunc","pol4", -1, 100);
  
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistRhovsCent);
  fOutputList->Add(fHistJetPt);
  fOutputList->Add(fHistJetArea);
  fOutputList->Add(fHistDeltaRhovsCent);
  fOutputList->Add(fHistDeltaRhoScalevsCent);
  fOutputList->Add(fHistJetPtvsCent);
  fOutputList->Add(fHistJetAreavsCent);
  fOutputList->Add(fHistNjetvsCent);
  
  fOutputList->Add(fHistRhovsNtrack);
  fOutputList->Add(fHistDeltaRhovsNtrack);
  fOutputList->Add(fHistDeltaRhoScalevsNtrack);
  fOutputList->Add(fHistJetPtvsNtrack);
  fOutputList->Add(fHistJetAreavsNtrack);
  fOutputList->Add(fHistNjetvsNtrack);

  //fOutputList->Add(fNewRhoFunction);
  
  PostData(1, fOutputList);
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskRho::GetScaleFactor(Double_t cent)
{
  Double_t scale = 1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}

//________________________________________________________________________
void AliAnalysisTaskRho::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  
  AliAnalysisTaskRhoBase::UserExec("");

  // add rho to event if not yet there
  if (!(InputEvent()->FindListObject(fRhoScaledName))) {
    new(fRhoScaled) TParameter<Double_t>(fRhoScaledName, 0);
    InputEvent()->AddObject(fRhoScaled);
  }

  // optimization in case autobranch loading is off
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");

  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;

  tracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  if (!tracks) {
  AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
   return;
  }
    
  jets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
  if (!jets) {
    AliError(Form("Pointer to jets %s == 0", fTracksName.Data() ));
    return;
  }

  fHistCentrality->Fill(fCent);

  const Int_t Ntracks = tracks->GetEntries();
  const Int_t Njets = jets->GetEntries();
  Int_t NjetAcc = 0;

  vector<Double_t> rhovec;
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
    //push all jets within selected acceptance into stack
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(iJets));
    if (!jet)
      continue; 
    if (jet->Area() < fAreaCut)
      continue;
    if ((jet->Phi() < fPhiMin) || (jet->Phi() > fPhiMax))
      continue;
    if ((jet->Eta() < fEtaMin) || (jet->Eta() > fEtaMax))
      continue;
    if (jet->Area() == 0)
      continue;
    NjetAcc++;
    rhovec.push_back(jet->Pt() / jet->Area());
    fHistJetPt->Fill(jet->Pt());
    fHistJetArea->Fill(jet->Area());
    fHistJetPtvsCent->Fill(fCent, jet->Pt());
    fHistJetPtvsNtrack->Fill(Ntracks, jet->Pt());
    fHistJetAreavsCent->Fill(fCent, jet->Area());
    fHistJetAreavsNtrack->Fill(Ntracks, jet->Area());
  }
  
  fHistNjetvsCent->Fill(fCent,NjetAcc);
  fHistNjetvsNtrack->Fill(Ntracks,NjetAcc);

  Double_t scale = GetScaleFactor(fCent);
  Double_t rhochem = GetRhoFactor(fCent);

  Double_t rho0 = -1;
  
  if (rhovec.size() > 0){
    //find median value
    Sort(rhovec);
    rho0 = GetMedian(rhovec,0);
    fHistRhovsCent->Fill(fCent, rho0);
    fHistDeltaRhovsCent->Fill(fCent, rho0 - rhochem);
    fHistDeltaRhoScalevsCent->Fill(fCent, rho0 * scale - rhochem);
    fHistRhovsNtrack->Fill(Ntracks, rho0);
    fHistDeltaRhovsNtrack->Fill(Ntracks, rho0 - rhochem);
    fHistDeltaRhoScalevsNtrack->Fill(Ntracks, rho0 * scale - rhochem);
  }

  fRho->SetVal(rho0);

  Double_t rho_scaled = rho0 * scale;

  fRhoScaled->SetVal(rho_scaled);
  
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskRho::Sort(vector<Double_t>& v)
{
  vector<Double_t> temp;
  temp.push_back(v[0]);
  for (UInt_t i = 1; i < v.size(); i++) {
    Bool_t insert = kFALSE;
    for (vector<Double_t>::iterator j = temp.begin(); j < temp.end(); j++){
      if (v[i]> * j) {
	temp.insert(j, v[i]);
	insert = kTRUE;
	j = temp.end();
      }
    }
    if (!insert)
      temp.insert(temp.end(), v[i]);
  }
  v = temp;
  return;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRho::GetMedian(vector<Double_t> v, Int_t c)
{
  if (v.size() == 0)
    return -1;
  if ((v.size() < 2) && (c == 1))
    return -1;
  if ((v.size() < 3) && ( c== 2))
    return -1;
  
  if (c == 1)
    v.erase(v.begin());
  
  if (c == 2) {
    v.erase(v.begin());
    v.erase(v.begin());
  }

  Float_t middle = (Float_t)v.size() / 2.0;
  if (middle != ceil(middle)) {
    //odd number
    return v[floor(middle)];
  }
  else{
    //even
    return (v[middle] + v[middle-1]) / 2.0;
  } 
}

//________________________________________________________________________
void AliAnalysisTaskRho::Terminate(Option_t *) 
{
  /*
  fHistRhovsCent->Fit(fNewRhoFunction, "NO");
  PostData(1, fOutputList);
  */
}
