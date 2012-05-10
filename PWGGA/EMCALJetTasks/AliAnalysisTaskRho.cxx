#include <TChain.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TVector3.h>

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
  AliAnalysisTaskSE(),
  fTracksName("tracks"),
  fJetsName("jets"),
  fRhosName("fArrRhos"),
  fClustersName("clusters"),
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
  fArrRhos(0),
  fPhiMin(0),
  fPhiMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fAreaCut(0),
  fCswitch(0)
{
  // Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho(const char *name) :
  AliAnalysisTaskSE(name),
  fTracksName("tracks"),
  fJetsName("jets"),
  fRhosName("fArrRhos"),
  fClustersName("clusters"),
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
  fArrRhos(0),
  fPhiMin(0),
  fPhiMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fAreaCut(0),
  fCswitch(0)

{
  // Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskRho::UserCreateOutputObjects()
{

  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }

  fArrRhos = new TClonesArray("TVector3");
  fArrRhos->SetName(fRhosName);

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
  
  PostData(1, fOutputList);
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskRho::GetScaleFactor(Double_t cent)
{
  Double_t scale = 0.000066*cent*cent-0.0015*cent+1.5;
  return scale;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRho::GetRhoFactor(Double_t cent)
{
  Double_t RhoChEmArray[100] = {201.32,  191.28,  183.72,  176.70,  170.10,  163.63,  157.76,  152.13,  146.62,  141.28,  136.41,  131.69,  126.80,  121.95,  117.85,  113.44,  109.19,  105.35,  101.24,  97.23,  93.49,  90.15,  86.77,  83.26,  79.92,  76.77,  73.75,  70.92,  67.93,  65.09,  62.39,  59.57,  57.26,  54.65,  52.24,  50.03,  47.70,  45.34,  43.29,  41.32,  39.35,  37.35,  35.46,  33.64,  32.15,  30.43,  28.81,  27.28,  25.85,  24.36,  23.09,  21.83,  20.60,  19.42,  18.45,  17.44,  16.62,  15.40,  14.06,  13.31,  12.72,  11.91,  11.01,  10.54,  9.89,  9.23,  8.73,  8.15,  7.73,  7.17,  7.00,  6.55,  6.07,  5.78,  5.42,  5.29,  5.13,  4.69,  4.74,  4.36,  4.40,  4.21,  4.23,  3.77,  3.83,  4.23,  4.17,  4.53,  4.06,  3.63,  3.90,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00};
  Double_t rho = RhoChEmArray[(Int_t)floor(cent)];

  return rho;
}

//________________________________________________________________________
void AliAnalysisTaskRho::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // add rhos to event if not yet there
  if (!(InputEvent()->FindListObject(fRhosName)))
    InputEvent()->AddObject(fArrRhos);

  // optimization in case autobranch loading is off
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");

  // get centrality 
  Double_t fCent = -1; 
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  TList *l = InputEvent()->GetList();
  if (centrality)
    fCent = centrality->GetCentralityPercentile("V0M");
  else
    fCent=99; // probably pp data
  if (fCent<0) {
    AliError(Form("Centrality negative: %f", fCent));
    return;
  }

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

  fHistCentrality->Fill(fCent);

  Double_t scale = GetScaleFactor(fCent);
  Double_t rhochem = GetRhoFactor(fCent);

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
    rhovec.push_back(jet->Pt()/jet->Area());
    fHistJetPt->Fill(jet->Pt());
    fHistJetArea->Fill(jet->Area());
    fHistJetPtvsCent->Fill(fCent,jet->Pt());
    fHistJetPtvsNtrack->Fill(Ntracks,jet->Pt());
    fHistJetAreavsCent->Fill(fCent,jet->Area());
    fHistJetAreavsNtrack->Fill(Ntracks,jet->Area());
  }
  
  fHistNjetvsCent->Fill(fCent,NjetAcc);
  fHistNjetvsNtrack->Fill(Ntracks,NjetAcc);
  Double_t rho0 = -1;
  
  if (rhovec.size()>0){
    //find median value
    Sort(rhovec);
    rho0 = GetMedian(rhovec,0);
    fHistRhovsCent->Fill(fCent,rho0);
    fHistDeltaRhovsCent->Fill(fCent,rho0-rhochem);
    fHistDeltaRhoScalevsCent->Fill(fCent,rho0*scale-rhochem);
    fHistRhovsNtrack->Fill(Ntracks,rho0);
    fHistDeltaRhovsNtrack->Fill(Ntracks,rho0-rhochem);
    fHistDeltaRhoScalevsNtrack->Fill(Ntracks,rho0*scale-rhochem);
  }
  
  //Scale rho for method 2, task does not know whether jet collection
  //is charged only

  Double_t rho2=scale*rho0;
    
  Int_t irx = 0;
  TVector3 *ArrRhos = new ((*fArrRhos)[irx]) TVector3;
  ArrRhos->SetXYZ(rho0, rho2, rhochem);
  irx++;
  
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskRho::Terminate(Option_t *) 
{

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
