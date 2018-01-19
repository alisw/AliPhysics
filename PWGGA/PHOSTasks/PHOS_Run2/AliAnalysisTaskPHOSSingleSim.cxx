#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "THnSparse.h"
#include "THashList.h"
#include "TMath.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"

#include "AliESDtrackCuts.h"
#include "AliESDHeader.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"

#include "AliMultSelection.h"
#include "AliEventplane.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliOADBContainer.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAODMCHeader.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"


#include "AliPHOSEventCuts.h"
#include "AliPHOSClusterCuts.h"
#include "AliPHOSJetJetMC.h"
#include "AliPHOSTriggerHelper.h"
#include "AliAnalysisTaskPHOSPi0EtaToGammaGamma.h"

#include "AliAnalysisTaskPHOSSingleSim.h"

// Author: Daiki Sekihata (Hiroshima University)

using namespace std;

ClassImp(AliAnalysisTaskPHOSSingleSim)

//________________________________________________________________________
AliAnalysisTaskPHOSSingleSim::AliAnalysisTaskPHOSSingleSim(const char *name):
  AliAnalysisTaskPHOSPi0EtaToGammaGamma(name),
  fParticleName(""),
  fMCArray(0x0)
{
  // Constructor


  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());

}
//________________________________________________________________________
AliAnalysisTaskPHOSSingleSim::~AliAnalysisTaskPHOSSingleSim()
{

}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  AliAnalysisTaskPHOSPi0EtaToGammaGamma::UserCreateOutputObjects();

//  const Int_t NpTgg = 101;
//  Double_t pTgg[NpTgg]={};
//  for(Int_t i=0;i<50;i++)     pTgg[i] = 0.1 * i;            //every 0.1 GeV/c, up to 5 GeV/c
//  for(Int_t i=50;i<60;i++)    pTgg[i] = 0.5 * (i-50) + 5.0; //every 0.5 GeV/c, up to 10 GeV/c
//  for(Int_t i=60;i<NpTgg;i++) pTgg[i] = 1.0 * (i-60) + 10.0;//every 1.0 GeV/c, up to 50 GeV/c
//
//  const Int_t Npar = 3;
//  const TString parname[Npar] = {"Pi0","Eta","Gamma"};
//  for(Int_t ipar=0;ipar<Npar;ipar++){
//    TH1F *h1Pt = new TH1F(Form("hGenEmbedded%sPt",parname[ipar].Data()        ),Form("generated %s pT",parname[ipar].Data()        ),NpTgg-1,pTgg);
//    h1Pt->Sumw2();
//    fOutputContainer->Add(h1Pt);
//
//    TH2F *h2EtaPhi = new TH2F(Form("hGenEmbedded%sEtaPhi",parname[ipar].Data()),Form("generated %s eta vs phi",parname[ipar].Data()),200,-1,1,60,0,TMath::TwoPi());
//    h2EtaPhi->Sumw2();
//    fOutputContainer->Add(h2EtaPhi);
//
//    TH2F *h2EtaPt = new TH2F(Form("hGenEmbedded%sEtaPt",parname[ipar].Data()  ),Form("generated %s eta vs pT",parname[ipar].Data() ),200,-1,1,NpTgg-1,pTgg);
//    h2EtaPt->Sumw2();
//    fOutputContainer->Add(h2EtaPt);
//
//  }//end of particle loop


  PostData(1,fOutputContainer);

}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::UserExec(Option_t *option) 
{
  // Main loop
  // Called for each event

  if(!fPHOSEventCuts){
    AliError("fPHOSEventCuts is not set! return");
    return;
  }
  if(!fPHOSClusterCuts){
    AliError("fPHOSClusterCuts is not set! return");
    return;
  }

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }

  fESDEvent = dynamic_cast<AliESDEvent*>(fEvent);
  fAODEvent = dynamic_cast<AliAODEvent*>(fEvent);

  FillHistogramTH1(fOutputContainer,"hEventSummary",1);//all

  const AliVVertex *vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();

  FillHistogramTH1(fOutputContainer,"hVertexZ" ,fVertex[2]);
  fZvtx = (Int_t)((fVertex[2]+10.)/2.);//it should be 0-9.
  if(fZvtx < 0) fZvtx = 0;//protection to avoid fZvtx = -1.
  if(fZvtx > 9) fZvtx = 9;//protection to avoid fZvtx = 10.

  fPHOSClusterArray = (TClonesArray*)fEvent->FindListObject("PHOSClusterArray");
  if(!fPHOSClusterArray){
    AliWarning("fPHOSClusterArray object not found!");
    return;
  }

  GetMCInfo();

  if(!fMCArray){
    AliError("Could not retrieve AOD event!");
    return;
  }

  fCentralityMain = 0.5;
  fEPBin = 0;
  //<- end of event selection
  //-> start physics analysis
  if(fIsPHOSTriggerAnalysis && !(fPHOSTriggerHelper->IsPHI7(fEvent,fPHOSClusterCuts,fEmin))){
    AliInfo("event is rejected. IsPHI7 = kFALSE.");
    return;
  }

  FillHistogramTH1(fOutputContainer,"hVertexZSelectEvent" ,fVertex[2]);
  FillHistogramTH1(fOutputContainer,"hEventSummary",2);//selected event

  if(fRunNumber != fEvent->GetRunNumber()){ // Check run number
    fRunNumber = fEvent->GetRunNumber();
    fPHOSGeo = GetPHOSGeometry();
  }

  if(!fPHOSEvents[fZvtx][fEPBin]) fPHOSEvents[fZvtx][fEPBin] = new TList();
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  ProcessMC();
  SetWeightToClusters();

  AliAnalysisTaskPHOSPi0EtaToGammaGamma::ClusterQA();
  FillPhoton();
  if(fParticleName.Contains("Pi0") || fParticleName.Contains("Eta")){
    FillMgg();
    FillMixMgg();
  }
  AliAnalysisTaskPHOSPi0EtaToGammaGamma::EstimatePIDCutEfficiency();

  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSClusterArray->GetEntriesFast() > 0){
    TClonesArray *clone = new TClonesArray(*fPHOSClusterArray);

    //prevPHOS->AddFirst(fPHOSClusterArray);
    prevPHOS->AddFirst(clone);
    //fPHOSClusterArray=0;
    //clone = 0;

    if(prevPHOS->GetSize() > fNMixed){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last());
      prevPHOS->RemoveLast();
      delete tmp;
      tmp = NULL;
    }
  }

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::Terminate(Option_t *option) 
{
  //Called once at the end of the query
  //In principle, this function is not needed...

  AliInfo(Form("%s is done.",GetName()));

}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::ProcessMC()
{
  //This is for analyzing general purpose MC such as pure PYTHIA, HIJING, DPMJET, PHOJET and so on.
  //get MC information
  TF1 *f1weight = 0x0;
  if(fParticleName.Contains("Pi0"))        f1weight = GetAdditionalPi0PtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Eta"))   f1weight = GetAdditionalEtaPtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Gamma")) f1weight = GetAdditionalGammaPtWeightFunction(fCentralityMain);

  AliAODMCParticle *p_origin = (AliAODMCParticle*)fMCArray->At(0);//0 is always generated particle by AliGenBox.
  Double_t pT_origin = p_origin->Pt();

  //printf("pT_orgin = %f GeV/c.\n",pT_origin);
  //cout <<"fParticleName = " << fParticleName << endl;
  //cout <<"f1weight = " << f1weight << endl;

  Double_t weight = f1weight->Eval(pT_origin) * pT_origin;

  //Int_t genID = -1;
  Double_t pT=0, rapidity=0, phi=0;
  Int_t pdg = 0;
  TString parname = "";
  TString genname = "";

  const Int_t Ntrack = fMCArray->GetEntriesFast();

  for(Int_t i=0;i<Ntrack;i++){
    AliAODMCParticle *p = (AliAODMCParticle*)fMCArray->At(i);
    //genID = p->GetGeneratorIndex();
    pT = p->Pt();
    rapidity = p->Y();
    phi = p->Phi();
    pdg = p->PdgCode();

    //rapidity is Y(), but, pseudo-rapidity is Eta();

    if(pT < 1e-3) continue;//reject below 1 MeV
    if(TMath::Abs(rapidity) > 0.5) continue;

    //printf("pdg = %d , Rho = %e cm\n",pdg,Rho(p));
    //if(RhoEMB(p) > 1.0) continue;

    if(pdg==111){//pi0
      parname = "Pi0";
    }
    else if(pdg==221){//eta
      parname = "Eta";

    }
    else if(pdg==22){//gamma
      parname = "Gamma";
    }
    else{
      continue;
    }

    //Double_t eta = p->Eta();
    //Printf("particle %d is generated at eta = %f , phi = %f and pT = %f.",pdg,eta,phi,pT);

    FillHistogramTH1(fOutputContainer,Form("hGen%sPt"    ,parname.Data()),pT          ,weight);
    FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPhi",parname.Data()),phi,rapidity,weight);
    FillHistogramTH2(fOutputContainer,Form("hGen%sEtaPt" ,parname.Data()),rapidity,pT ,weight);

  }//end of generated particle loop

}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::SetWeightToClusters()
{
  TF1 *f1weight = 0x0;
  if(fParticleName.Contains("Pi0"))        f1weight = GetAdditionalPi0PtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Eta"))   f1weight = GetAdditionalEtaPtWeightFunction(fCentralityMain);
  else if(fParticleName.Contains("Gamma")) f1weight = GetAdditionalGammaPtWeightFunction(fCentralityMain);

  AliAODMCParticle *p_origin = (AliAODMCParticle*)fMCArray->At(0);//0 is always generated particle by AliGenBox.
  Double_t pT_origin = p_origin->Pt();
  Double_t weight = f1weight->Eval(pT_origin) * pT_origin;

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  for(Int_t i=0;i<multClust;i++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(i);
    ph->SetWeight(weight);

  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskPHOSSingleSim::FindCommonParent(Int_t iPart, Int_t jPart)
{
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart

  Int_t ntrack = fMCArray->GetEntriesFast();
  if(iPart==-1 || iPart>=ntrack || jPart==-1 || jPart>=ntrack) return -1;

  Int_t iprim1 = iPart;

  while(iprim1>-1){
    Int_t iprim2=jPart;

    while(iprim2>-1){
      if(iprim1==iprim2) return iprim1;
      //iprim2 = GetParticle(iprim2)->GetMother();
      iprim2 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iprim2))->GetMother();
    }

    //iprim1 = GetParticle(iprim1)->GetMother();
    //iprim1 = (AliAODMCParticle*)(fMCArray->At(iprim1))->GetMother();
    iprim1 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iprim1))->GetMother();
  }

  return -1;
}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::GetMCInfo()
{
  fMCArray = 0x0;
  fMCArray = (TClonesArray*)GetMCInfoAOD();

  if(!fMCArray){
    AliError("Could not retrieve AOD event!");
    return;
  }
}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::FillPhoton() 
{
  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  Double_t pT=0,energy=0;
  Double_t phi = -999, dphi = -999.;
  Double_t eff=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t trgeff=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();
  Double_t value[2] = {};
  Double_t sp1 = -999;

  Double_t weight = 1.;

  for(Int_t iph=0;iph<multClust;iph++){
    AliCaloPhoton *ph = (AliCaloPhoton*)fPHOSClusterArray->At(iph);
    if(!fPHOSClusterCuts->AcceptPhoton(ph)) continue;

     if(fIsPHOSTriggerAnalysis){
      if(!fPHOSTriggerHelper->IsOnActiveTRUChannel(ph)) continue;
      if(!fIsMC && !ph->IsTrig()) continue;//it is meaningless to focus on photon without fired trigger in PHOS triggered data.
    }

    weight = 1.;
    if(fIsMC){
      weight = ph->GetWeight();
    }

    pT = ph->Pt();
    energy = ph->Energy();
    phi = ph->Phi();

    if(fUseCoreEnergy){
      pT = (ph->GetMomV2())->Pt();
      energy = (ph->GetMomV2())->Energy();
      phi = (ph->GetMomV2())->Phi();
    }

    eff = f1tof->Eval(energy);

    if(!fIsMC && fIsPHOSTriggerAnalysis){
      trgeff  = f1trg->Eval(energy);
    }

    //0 < photon phi < 2pi
    if(phi < 0) phi += TMath::TwoPi();
    TVector2 vg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

    if(fIsFlowTask){
      dphi = DeltaPhiIn0Pi(phi - fEventPlane);
      sp1 = vg * fQVector1;
      if(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[1] = TMath::Cos(fHarmonics * dphi);
      else if(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[1] = sp1;
      else                                                value[1] = 0;
    }
    else{
      dphi = phi;
      sp1 = 0;
      value[1] = 0;
    }

    value[0] = pT;

    FillSparse(fOutputContainer,"hSparsePhoton",value,weight * 1/trgeff);

    if(ph->IsTOFOK()){
      FillSparse(fOutputContainer,"hSparsePhoton_TOF",value,1/eff * weight * 1/trgeff);
    }

  }//end of cluster loop

}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::FillMgg() 
{

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();
  TLorentzVector p12, p12core;

  Double_t m12=0,pt12=0,asym=0;
  Double_t e1=0,e2=0;
  Double_t phi = -999, dphi = -999.;

  Double_t eff1=1, eff2=1, eff12=1;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();

  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  Double_t value[4] = {};
  Double_t sp1 = -999;

  Double_t weight = 1., w1 = 1., w2 = 1.;

  for(Int_t i1=0;i1<multClust-1;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;
    if(fIsPHOSTriggerAnalysis && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph1)) continue;

    for(Int_t i2=i1+1;i2<multClust;i2++){
      AliCaloPhoton *ph2 = (AliCaloPhoton*)fPHOSClusterArray->At(i2);
      if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;
      if(fIsPHOSTriggerAnalysis && !fPHOSTriggerHelper->IsOnActiveTRUChannel(ph2)) continue;

      if(!fIsMC && fIsPHOSTriggerAnalysis && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.

      e1 = ph1->Energy();
      e2 = ph2->Energy();

      p12  = *ph1 + *ph2;
      m12  = p12.M();
      pt12 = p12.Pt();
      phi  = p12.Phi();
      asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));//always full energy

      if(fUseCoreEnergy){
        p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        m12     = p12core.M();
        pt12    = p12core.Pt();
        phi     = p12core.Phi();

        e1 = (ph1->GetMomV2())->Energy();
        e2 = (ph2->GetMomV2())->Energy();
        asym = TMath::Abs(e1 - e2) / (e1 + e2);
      }

      eff1 = f1tof->Eval(e1);
      eff2 = f1tof->Eval(e2);
      eff12 = eff1 * eff2;

      if(!fIsMC && fIsPHOSTriggerAnalysis){
        trgeff1  = f1trg->Eval(e1);
        trgeff2  = f1trg->Eval(e2);
        trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR
      }

      weight = 1.;
      if(fIsMC){
        w1 = ph1->GetWeight();
        w2 = ph2->GetWeight();
        weight = w1;//common weighting to all generated particles in embedding.

      }//end of if fIsMC

      if(phi < 0) phi += TMath::TwoPi();

      TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

      if(fIsFlowTask){
        dphi = DeltaPhiIn0Pi(phi - fEventPlane);
        sp1 = vgg * fQVector1;
        if(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[3] = TMath::Cos(fHarmonics * dphi);
        else if(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[3] = sp1;
        else                                                value[3] = 0;
      }
      else{
        dphi = phi;
        sp1 = 0;
        value[3] = 0;
      }

      value[0] = m12;
      value[1] = pt12;
      value[2] = asym;

      if(m12 > 0.96) continue;//reduce entry in THnSparse

      if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,weight * 1/trgeff12);
      FillSparse(fOutputContainer,"hSparseMgg",value,weight * 1/trgeff12);

      if(ph1->IsTOFOK() && ph2->IsTOFOK()){

        FillSparse(fOutputContainer,"hSparseMgg_TOF",value,1/eff12 * weight * 1/trgeff12);

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12 * weight * 1/trgeff12);


      }//end of TOF cut

    }//end of ph2

  }//end of ph1

}
//________________________________________________________________________
void AliAnalysisTaskPHOSSingleSim::FillMixMgg() 
{
  TList *prevPHOS = fPHOSEvents[fZvtx][fEPBin];

  const Int_t multClust = fPHOSClusterArray->GetEntriesFast();

  TLorentzVector p12, p12core;
  Double_t m12=0,pt12=0,asym=0;
  Double_t phi = -999, dphi = -999.;
  Double_t weight = 1., w1 = 1., w2 = 1.;

  Double_t eff1=1, eff2=1, eff12=1;
  Double_t e1=0,e2=0;
  TF1 *f1tof = GetTOFCutEfficiencyFunction();
  Double_t trgeff1=1;
  Double_t trgeff2=1;
  Double_t trgeff12=1;
  TF1 *f1trg = GetTriggerEfficiencyFunction();

  Double_t value[4] = {};
  Double_t sp1 = -999;

  for(Int_t i1=0;i1<multClust;i1++){
    AliCaloPhoton *ph1 = (AliCaloPhoton*)fPHOSClusterArray->At(i1);
    if(!fPHOSClusterCuts->AcceptPhoton(ph1)) continue;

    for(Int_t ev=0;ev<prevPHOS->GetSize();ev++){
      TClonesArray *mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2=0;i2<mixPHOS->GetEntriesFast();i2++){
        AliCaloPhoton *ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if(!fPHOSClusterCuts->AcceptPhoton(ph2)) continue;

        if(!fIsMC && fIsPHOSTriggerAnalysis && (!ph1->IsTrig() && !ph2->IsTrig())) continue;//it is meaningless to reconstruct invariant mass with FALSE-FALSE combination in PHOS triggered data.

        e1 = ph1->Energy();
        e2 = ph2->Energy();

        p12  = *ph1 + *ph2;
        m12  = p12.M();
        pt12 = p12.Pt();
        phi  = p12.Phi();
        asym = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

        if(fUseCoreEnergy){
          p12core = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
          m12     = p12core.M();
          pt12    = p12core.Pt();
          phi     = p12core.Phi();

          e1 = (ph1->GetMomV2())->Energy();
          e2 = (ph2->GetMomV2())->Energy();
        }

        eff1 = f1tof->Eval(e1);
        eff2 = f1tof->Eval(e2);
        eff12 = eff1 * eff2;
        weight = 1.;

        if(!fIsMC && fIsPHOSTriggerAnalysis){
          trgeff1  = f1trg->Eval(e1);
          trgeff2  = f1trg->Eval(e2);
          trgeff12 = trgeff1 + trgeff2 - (trgeff1 * trgeff2);//logical OR
        }

        if(fIsMC){
          w1= ph1->GetWeight();
          w2 = ph2->GetWeight();

          weight = w1*w2;

        }//end of if fIsMC

        if(phi < 0) phi += TMath::TwoPi();
        TVector2 vgg(TMath::Cos(fHarmonics * phi),TMath::Sin(fHarmonics * phi));

        if(fIsFlowTask){
          dphi = DeltaPhiIn0Pi(phi - fEventPlane);
          sp1 = vgg * fQVector1;
          if(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kEP)      value[3] = TMath::Cos(fHarmonics * dphi);
          else if(AliAnalysisTaskPHOSPi0EtaToGammaGamma::kSP) value[3] = sp1;
          else                                                value[3] = 0;
        }
        else{
          dphi = phi;
          sp1 = 0;
          value[3] = 0;
        }

        value[0] = m12;
        value[1] = pt12;
        value[2] = asym;
        if(m12 > 0.96) continue;//reduce entry in THnSparse

        if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12);
        FillSparse(fOutputContainer,"hSparseMixMgg",value,weight * 1/trgeff12);

        //FillHistogramTH3(fOutputContainer,"hMixMgg",m12,pt12,TMath::Cos(fHarmonics * dphi),weight);
        //if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMixMgg_asym08",m12,pt12,TMath::Cos(fHarmonics * dphi),weight);

        if(ph1->IsTOFOK() && ph2->IsTOFOK()){
          FillSparse(fOutputContainer,"hSparseMixMgg_TOF",value,1/eff12 * weight * 1/trgeff12);
          if(TMath::Abs(ph1->Module()-ph2->Module()) < 2) FillHistogramTH2(fOutputContainer,Form("hMixMgg_M%d%d_TOF",TMath::Min(ph1->Module(),ph2->Module()), TMath::Max(ph1->Module(),ph2->Module())),m12,pt12,1/eff12);

          //FillHistogramTH3(fOutputContainer,"hMixMgg_TOF",m12,pt12,TMath::Cos(fHarmonics * dphi),1/eff12 * weight);
          //if(asym < 0.8) FillHistogramTH3(fOutputContainer,"hMixMgg_TOF_asym08",m12,pt12,TMath::Cos(fHarmonics * dphi),1/eff12 * weight);

        }//end of TOF cut

      }//end of ph2

    }//end of mix

  }//end of ph1

}
//________________________________________________________________________

