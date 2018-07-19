/*
 * AliAnalysisTaskMCFemtoDream.cxx
 *
 *  Created on: Mar 15, 2018
 *      Author: hohlweger
 */

#include <iostream>
#include "AliAnalysisTaskMCFemtoDream.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"

ClassImp(AliAnalysisTaskMCFemtoDream);

AliAnalysisTaskMCFemtoDream::AliAnalysisTaskMCFemtoDream()
:AliAnalysisTaskSE()
,foutput()
,fPart()
,fcollConfig()
,fPartCollAll()
,fPartCollAccept()
,fCounter()
,fProtonEta()
,fProtonPt()
,fProtonResonances()
,fAntiProtonEta()
,fAntiProtonPt()
,fAntiProtonResonances()
,fLambdaEta()
,fLambdaPt()
,fLambdaResonances()
,fAntiLambdaEta()
,fAntiLambdaPt()
,fAntiLambdaResonances()
,fXiEta()
,fXiPt()
,fXiResonances()
,fAntiXiEta()
,fAntiXiPt()
,fAntiXiResonances()
{

}

AliAnalysisTaskMCFemtoDream::AliAnalysisTaskMCFemtoDream(const char *name)
:AliAnalysisTaskSE(name)
,foutput()
,fPart(new AliFemtoDreamBasePart())
,fcollConfig()
,fPartCollAll()
,fPartCollAccept()
,fCounter()
,fProtonEta()
,fProtonPt()
,fProtonResonances()
,fAntiProtonEta()
,fAntiProtonPt()
,fAntiProtonResonances()
,fLambdaEta()
,fLambdaPt()
,fLambdaResonances()
,fAntiLambdaEta()
,fAntiLambdaPt()
,fAntiLambdaResonances()
,fXiEta()
,fXiPt()
,fXiResonances()
,fAntiXiEta()
,fAntiXiPt()
,fAntiXiResonances()
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskMCFemtoDream::~AliAnalysisTaskMCFemtoDream() {
  // TODO Auto-generated destructor stub
}

void AliAnalysisTaskMCFemtoDream::UserCreateOutputObjects() {
  foutput=new TList();
  foutput->SetName("OutputList");
  foutput->SetOwner();

  fcollConfig=new AliFemtoDreamCollConfig("AllPart","AllPart");
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);
  PDGParticles.push_back(2212);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3122);
  PDGParticles.push_back(3312);
  PDGParticles.push_back(3312);
  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(10);
  std::vector<int> NBins;
  NBins.push_back(750);
  NBins.push_back(750);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(750);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  NBins.push_back(150);
  //std::vector<double> kMin= {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  std::vector<float> kMin;
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  kMin.push_back(0.);
  //std::vector<double> kMax= {3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.};
  std::vector<float> kMax;
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  kMax.push_back(3.);
  std::vector<int> MultBins;
  MultBins.push_back(0);
  MultBins.push_back(4);
  MultBins.push_back(8);
  MultBins.push_back(12);
  MultBins.push_back(16);
  MultBins.push_back(20);
  MultBins.push_back(24);
  MultBins.push_back(28);
  MultBins.push_back(32);
  MultBins.push_back(36);
  MultBins.push_back(40);
  MultBins.push_back(60);
  MultBins.push_back(80);

  fcollConfig->SetMultBins(MultBins);
  fcollConfig->SetMultBinning(true);
  fcollConfig->SetZBins(ZVtxBins);
  fcollConfig->SetPDGCodes(PDGParticles);
  fcollConfig->SetNBinsHist(NBins);
  fcollConfig->SetMinKRel(kMin);
  fcollConfig->SetMaxKRel(kMax);
  fcollConfig->SetMixingDepth(10);
  fcollConfig->SetSECommonAncestor(true);

  fPartCollAll=new AliFemtoDreamPartCollection(fcollConfig,false);
  fPartCollAccept=new AliFemtoDreamPartCollection(fcollConfig,false);

  TList *collAll=new TList();
  collAll->SetName("CollAll");
  collAll->SetOwner();
  collAll->Add(fPartCollAll->GetHistList());
  collAll->Add(fPartCollAll->GetQAList());

  TList *collAccept=new TList();
  collAccept->SetName("CollAccept");
  collAccept->SetOwner();
  collAccept->Add(fPartCollAccept->GetHistList());
  collAccept->Add(fPartCollAccept->GetQAList());

  foutput->Add(collAll);
  foutput->Add(collAccept);
  fCounter=new TH1F("Counter","Counter",16,0.5,16.5);
  fCounter->GetXaxis()->SetBinLabel(1,"Proton All");
  fCounter->GetXaxis()->SetBinLabel(2,"Proton Accepted");
  fCounter->GetXaxis()->SetBinLabel(3,"AntiProton All");
  fCounter->GetXaxis()->SetBinLabel(4,"AntiProton Accepted");
  fCounter->GetXaxis()->SetBinLabel(5,"Lambda All");
  fCounter->GetXaxis()->SetBinLabel(6,"Lambda Accepted");
  fCounter->GetXaxis()->SetBinLabel(7,"AntiLambda All");
  fCounter->GetXaxis()->SetBinLabel(8,"AntiLambda Accepted");
  fCounter->GetXaxis()->SetBinLabel(9,"Xi All");
  fCounter->GetXaxis()->SetBinLabel(10,"Xi Accepted");
  fCounter->GetXaxis()->SetBinLabel(11,"AntiXi All");
  fCounter->GetXaxis()->SetBinLabel(12,"AntiXi Accepted");
  fCounter->Sumw2();

  foutput->Add(fCounter);

  fProtonEta= new TH1F("fProtonEta","fProtonEta",1000,-10,10);
  fProtonEta->Sumw2();
  foutput->Add(fProtonEta);
  fProtonPt= new TH1F("fProtonPt","fProtonPt",1000,0,100);
  fProtonPt->Sumw2();
  foutput->Add(fProtonPt);
  fProtonResonances= new TH1F("fProtonResonances","fProtonResonances",10000,0.5,10000);
  fProtonResonances->Sumw2();
  foutput->Add(fProtonResonances);

  fAntiProtonEta= new TH1F("fAntiProtonEta","fAntiProtonEta",1000,-10,10);
  fAntiProtonEta->Sumw2();
  foutput->Add(fAntiProtonEta);
  fAntiProtonPt= new TH1F("fAntiProtonPt","fAntiProtonPt",1000,0,100);
  fAntiProtonPt->Sumw2();
  foutput->Add(  fAntiProtonPt);
  fAntiProtonResonances= new TH1F("fAntiProtonResonances","fAntiProtonResonances",10000,0.5,10000);
  fAntiProtonResonances->Sumw2();
  foutput->Add(fAntiProtonResonances);

  fLambdaEta= new TH1F("fLambdaEta","fLambdaEta",1000,-10,10);
  fLambdaEta->Sumw2();
  foutput->Add(  fLambdaEta);
  fLambdaPt= new TH1F("fLambdaPt","fLambdaPt",1000,0,100);
  fLambdaPt->Sumw2();
  foutput->Add(  fLambdaPt);
  fLambdaResonances= new TH1F("fLambdaResonances","fLambdaResonances",10000,0.5,10000);
  fLambdaResonances->Sumw2();
  foutput->Add(fLambdaResonances);

  fAntiLambdaEta= new TH1F("fAntiLambdaEta","fAntiLambdaEta",1000,-10,10);
  fAntiLambdaEta->Sumw2();
  foutput->Add(fAntiLambdaEta);
  fAntiLambdaPt= new TH1F("fAntiLambdaPt","fAntiLambdaPt",1000,0,100);
  fAntiLambdaPt->Sumw2();
  foutput->Add(fAntiLambdaPt);
  fAntiLambdaResonances= new TH1F("fAntiLambdaResonances","fAntiLambdaResonances",10000,0.5,10000);
  fAntiLambdaResonances->Sumw2();
  foutput->Add(fAntiLambdaResonances);

  fXiEta= new TH1F("fXiEta","fXiEta",1000,-10,10);
  fXiEta->Sumw2();
  foutput->Add(fXiEta);
  fXiPt= new TH1F("fXiPt","fXiPt",1000,0,100);
  fXiPt->Sumw2();
  foutput->Add(fXiPt);
  fXiResonances= new TH1F("fXiResonances","fXiResonances",10000,0.5,10000);
  fXiResonances->Sumw2();
  foutput->Add(fXiResonances);

  fAntiXiEta = new TH1F("fAntiXiEta","fAntiXiEta",1000,-10,10);
  fAntiXiEta->Sumw2();
  foutput->Add(fAntiXiEta);
  fAntiXiPt = new TH1F ("fAntiXiPt","fAntiXiPt",1000,0,100);
  fAntiXiPt->Sumw2();
  foutput->Add(fAntiXiPt);
  fAntiXiResonances= new TH1F("fAntiXiResonances","fAntiXiResonances",10000,0.5,10000);
  fAntiXiResonances->Sumw2();
  foutput->Add(fAntiXiResonances);

  PostData(1,foutput);
}

void AliAnalysisTaskMCFemtoDream::UserExec(Option_t *) {
  AliAODEvent *Event=static_cast<AliAODEvent*>(fInputEvent);
  AliAODHeader *header = dynamic_cast<AliAODHeader*>(Event->GetHeader());
  int mult = header->GetRefMultiplicityComb08();
  AliAODInputHandler *eventHandler =
      dynamic_cast<AliAODInputHandler*>(
          AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliMCEvent* fMC = eventHandler->MCEvent();
  // Loop over the MC tracks
  std::vector<AliFemtoDreamBasePart> AllProtons;
  std::vector<AliFemtoDreamBasePart> AllAntiProtons;
  std::vector<AliFemtoDreamBasePart> AllLambdas;
  std::vector<AliFemtoDreamBasePart> AllAntiLambdas;
  std::vector<AliFemtoDreamBasePart> AllXis;
  std::vector<AliFemtoDreamBasePart> AllAntiXis;
  std::vector<AliFemtoDreamBasePart> AcceptedProtons;
  std::vector<AliFemtoDreamBasePart> AcceptedAntiProtons;
  std::vector<AliFemtoDreamBasePart> AcceptedLambdas;
  std::vector<AliFemtoDreamBasePart> AcceptedAntiLambdas;
  std::vector<AliFemtoDreamBasePart> AcceptedXis;
  std::vector<AliFemtoDreamBasePart> AcceptedAntiXis;

  for(int iPart = 1; iPart < (fMC->GetNumberOfTracks()); iPart++) {
    AliAODMCParticle *mcPart  = (AliAODMCParticle*)fMC->GetTrack(iPart);
    if (mcPart->IsPhysicalPrimary()) {
      fPart->ResetMCInfo();
      fPart->SetMCParticle(mcPart,fMC);
      if (mcPart->PdgCode()==2212 && TMath::Abs(mcPart->Eta())<1e29){
        if (mcPart->GetMother() > 0) {
          fProtonResonances->Fill(TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode()));
//          if (TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode())>10000) {
//            std::cout << "Proton PDG Mother: " << fMC->GetTrack((mcPart->GetMother()))->PdgCode() << std::endl;
//          }
        }
        fProtonEta->Fill(mcPart->Eta());
        fProtonPt->Fill(mcPart->Pt());
        AllProtons.push_back(*fPart);
        fCounter->Fill(1);
        if (!(TMath::Abs(mcPart->Eta())<0.8)) {
          if (!(TMath::Abs(mcPart->Pt())>0.3)) {
            AcceptedProtons.push_back(*fPart);
            fCounter->Fill(2);
          }
        }
      } else if (mcPart->PdgCode()==-2212){
        if (mcPart->GetMother() > 0) {
          fAntiProtonResonances->Fill(TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode()));
//          if (TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode())>10000) {
//            std::cout << "AntiProton PDG Mother: " << fMC->GetTrack((mcPart->GetMother()))->PdgCode() << std::endl;
//          }
        }
        fAntiProtonEta->Fill(mcPart->Eta());
        fAntiProtonPt->Fill(mcPart->Pt());
        AllAntiProtons.push_back(*fPart);
        fCounter->Fill(3);
        if (!(TMath::Abs(mcPart->Eta())<0.8)) {
          if (!(TMath::Abs(mcPart->Pt())>0.3)) {
            AcceptedAntiProtons.push_back(*fPart);
            fCounter->Fill(4);
          }
        }
      } else if (mcPart->PdgCode()==3122){
        if (mcPart->GetMother() > 0) {
          fLambdaResonances->Fill(TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode()));
//          if (TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode())>10000) {
//            std::cout << "Lambda PDG Mother: " << fMC->GetTrack((mcPart->GetMother()))->PdgCode() << std::endl;
//          }
        }
        fLambdaEta->Fill(mcPart->Eta());
        fLambdaPt->Fill(mcPart->Pt());
        AllLambdas.push_back(*fPart);
        fCounter->Fill(5);
        bool pass=true;
        for (int iDaug=0;iDaug<mcPart->GetNDaughters();++iDaug) {
          if (mcPart->GetDaughter(iDaug)<0) {
            pass = false;
            break;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Eta())<0.8)) {
            pass = false;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Pt())>0.3)) {
            pass = false;
          }
        }
        if (pass ) {
          AcceptedLambdas.push_back(*fPart);
          fCounter->Fill(6);
        }
      } else if (mcPart->PdgCode()==-3122){
        if (mcPart->GetMother() > 0) {
          fAntiLambdaResonances->Fill(TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode()));
        }
        fAntiLambdaEta->Fill(mcPart->Eta());
        fAntiLambdaPt->Fill(mcPart->Pt());
        AllAntiLambdas.push_back(*fPart);
        fCounter->Fill(7);
        bool pass=true;
        for (int iDaug=0;iDaug<mcPart->GetNDaughters();++iDaug) {
          if (mcPart->GetDaughter(iDaug)<0) {
            pass = false;
            break;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Eta())<0.8)) {
            pass = false;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Pt())>0.3)) {
            pass = false;
          }
        }
        if (pass ) {
          AcceptedAntiLambdas.push_back(*fPart);
          fCounter->Fill(8);
        }
      } else if (mcPart->PdgCode()==3312){
        if (mcPart->GetMother() > 0) {
          fXiResonances->Fill(TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode()));
//          if (TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode())>10000) {
//            std::cout << "Xi PDG Mother: " << fMC->GetTrack((mcPart->GetMother()))->PdgCode() << std::endl;
//          }
        }
        fXiEta->Fill(mcPart->Eta());
        fXiPt->Fill(mcPart->Pt());
        AllXis.push_back(*fPart);
        fCounter->Fill(9);
        bool pass=true;
        for (int iDaug=0;iDaug<mcPart->GetNDaughters();++iDaug) {
          if (mcPart->GetDaughter(iDaug)<0) {
            pass = false;
            break;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Eta())<0.8)) {
            pass = false;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Pt())>0.3)) {
            pass = false;
          }
        }
        if (pass ) {
          AcceptedXis.push_back(*fPart);
          fCounter->Fill(10);
        }
      } else if (mcPart->PdgCode()==-3312) {
        if (mcPart->GetMother() > 0) {
          fAntiXiResonances->Fill(TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode()));
//          if (TMath::Abs(fMC->GetTrack((mcPart->GetMother()))->PdgCode())>10000) {
//            std::cout << "AntiXi PDG Mother: " << fMC->GetTrack((mcPart->GetMother()))->PdgCode() << std::endl;
//          }
        }
        fAntiXiEta->Fill(mcPart->Eta());
        fAntiXiPt->Fill(mcPart->Pt());
        AllAntiXis.push_back(*fPart);
        fCounter->Fill(11);
        bool pass=true;
        for (int iDaug=0;iDaug<mcPart->GetNDaughters();++iDaug) {
          if (mcPart->GetDaughter(iDaug)<0) {
            pass = false;
            break;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Eta())<0.8)) {
            pass = false;
          }
          if (!(TMath::Abs((fMC->GetTrack(mcPart->GetDaughter(iDaug)))->Pt())>0.3)) {
            pass = false;
          }
        }
        if (pass ) {
          AcceptedAntiXis.push_back(*fPart);
          fCounter->Fill(12);
        }
      }
    }
  }
  std::vector<std::vector<AliFemtoDreamBasePart>> ParticlesAll;
  ParticlesAll.push_back(AllProtons);
  ParticlesAll.push_back(AllAntiProtons);
  ParticlesAll.push_back(AllLambdas);
  ParticlesAll.push_back(AllAntiLambdas);
  ParticlesAll.push_back(AllXis);
  ParticlesAll.push_back(AllAntiXis);
  fPartCollAll->SetEvent(ParticlesAll,1.3,mult,0);
  std::vector<std::vector<AliFemtoDreamBasePart>> ParticlesAccepted;
  ParticlesAccepted.push_back(AcceptedProtons);
  ParticlesAccepted.push_back(AcceptedAntiProtons);
  ParticlesAccepted.push_back(AcceptedLambdas);
  ParticlesAccepted.push_back(AcceptedAntiLambdas);
  ParticlesAccepted.push_back(AcceptedXis);
  ParticlesAccepted.push_back(AcceptedAntiXis);
  fPartCollAccept->SetEvent(ParticlesAccepted,1.3,mult,0);

  PostData(1,foutput);
}
