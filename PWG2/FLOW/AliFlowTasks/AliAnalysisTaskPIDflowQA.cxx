/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial pures is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any pure. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// AliAnalysisTaskPIDflowQA:
// QA for pid
//
//origin: Marek Chojnacki, Marek.Chojnacki@cern.ch
//modified: Mikolaj Krzewicki, Mikolaj.Krzewicki@cern.ch

#include "AliAnalysisTaskPIDflowQA.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliFlowEventCuts.h"
#include "AliFlowTrackCuts.h"

ClassImp( AliAnalysisTaskPIDflowQA)

//________________________________________________________________________
AliAnalysisTaskPIDflowQA:: AliAnalysisTaskPIDflowQA():
  AliAnalysisTaskSE("AliAnalysisTaskPIDflowQA"),
  fESD(NULL),
  fCuts(NULL),
  fEventCuts(NULL),
  fESDpid(NULL),
  fMC(kFALSE),
  fITSsignal(NULL),
  fTPCsignal(NULL),
  fTOFsignal(NULL),
  fITSsignalpip(NULL),
  fTPCsignalpip(NULL),
  fTOFsignalpip(NULL),
  fITSsignalKp(NULL),
  fTPCsignalKp(NULL),
  fTOFsignalKp(NULL),
  fITSsignalpp(NULL),
  fTPCsignalpp(NULL),
  fTOFsignalpp(NULL),
  fITSsignalpiMCp(NULL),
  fTPCsignalpiMCp(NULL),
  fTOFsignalpiMCp(NULL),
  fITSsignalKMCp(NULL),
  fTPCsignalKMCp(NULL),
  fTOFsignalKMCp(NULL),
  fITSsignalpMCp(NULL),
  fTPCsignalpMCp(NULL),
  fTOFsignalpMCp(NULL),
  fTOFsignalBeta(NULL),
  fTOFsignalPiBeta(NULL),
  fTOFsignalKBeta(NULL),
  fTOFsignalPBeta(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fTPCvsGlobalMult(NULL),
  fStandardGlobalCuts(NULL),
  fStandardTPCCuts(NULL),
  fOutputList(NULL)
{
  //def ctor
}

//________________________________________________________________________
AliAnalysisTaskPIDflowQA:: AliAnalysisTaskPIDflowQA(const char *name):
  AliAnalysisTaskSE(name),
  fESD(NULL),
  fCuts(NULL),
  fEventCuts(NULL),
  fESDpid(NULL),
  fMC(kFALSE),
  fITSsignal(NULL),
  fTPCsignal(NULL),
  fTOFsignal(NULL),
  fITSsignalpip(NULL),
  fTPCsignalpip(NULL),
  fTOFsignalpip(NULL),
  fITSsignalKp(NULL),
  fTPCsignalKp(NULL),
  fTOFsignalKp(NULL),
  fITSsignalpp(NULL),
  fTPCsignalpp(NULL),
  fTOFsignalpp(NULL),
  fITSsignalpiMCp(NULL),
  fTPCsignalpiMCp(NULL),
  fTOFsignalpiMCp(NULL),
  fITSsignalKMCp(NULL),
  fTPCsignalKMCp(NULL),
  fTOFsignalKMCp(NULL),
  fITSsignalpMCp(NULL),
  fTPCsignalpMCp(NULL),
  fTOFsignalpMCp(NULL),
  fTOFsignalBeta(NULL),
  fTOFsignalPiBeta(NULL),
  fTOFsignalKBeta(NULL),
  fTOFsignalPBeta(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fTPCvsGlobalMult(NULL),
  fStandardGlobalCuts(NULL),
  fStandardTPCCuts(NULL),
  fOutputList(NULL)
{
  //Constructor
  fESDpid=new AliESDpid();

  //old
  fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,
      2.63394e+01,
      5.04114e-11,
      2.12543e+00,
      4.88663e+00 );
  //new
  //fESDpid->GetTPCResponse().SetBetheBlochParameters(1.28949/50.,
  //                                                  2.74095e+01,
  //                                                  TMath::Exp(-3.21763e+01),
  //                                                  2.44026,
  //                                                  6.58800);

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::UserCreateOutputObjects()
{
  //UserCreateOutputObject
  if (fOutputList) fOutputList->Delete();
  delete fOutputList;
  fOutputList=new TList();

  const  Int_t ndec=2;
  Int_t startvalue=-1;
  const  Int_t npredec=50;
  Double_t tabx[ndec*npredec+1];
  for (Int_t i=0; i<ndec; i++)
  {
    for (Int_t j=0; j<npredec; j++)
    {
      tabx[npredec*i+j]=TMath::Power(10,((Double_t)i)+((Double_t)startvalue)+((Double_t)j)/((Double_t)npredec));
    }
  }
  tabx[ndec*npredec]=TMath::Power(10,ndec+startvalue);

  fITSsignal=new TH2F("fITSsignal",";p [GeV/c];dEdx",ndec*npredec,tabx,900,0,900);
  fOutputList->Add(fITSsignal);
  fTPCsignal=new TH2F("fTPCsignal",";p [GeV/c];dEdx",ndec*npredec,tabx,500,0,500);
  fOutputList->Add(fTPCsignal);
  fTOFsignal=new TH2F("fTOFsignal",";p [GeV/c];t-t_{#pi}",ndec*npredec,tabx,1000,-2000,10000);
  fOutputList->Add(fTOFsignal);

  Int_t kPtBins=60;
  Double_t binsPtDummy[kPtBins+1];
  binsPtDummy[0]=0.0;
  for(int i=1; i<=kPtBins+1; i++)
  {
    if(binsPtDummy[i-1]+0.05<1.01)
      binsPtDummy[i]=binsPtDummy[i-1]+0.05;
    else
      binsPtDummy[i]=binsPtDummy[i-1]+0.1;
  }

  Int_t kPBins=60;
  Double_t binsPDummy[kPBins+1];
  binsPDummy[0]=0.0;
  for(int i=1; i<=kPBins+1; i++)
  {
    if(binsPDummy[i-1]+0.05<1.01)
      binsPDummy[i]=binsPDummy[i-1]+0.05;
    else
      binsPDummy[i]=binsPDummy[i-1]+0.1;
  }

  fITSsignalpip=new TH2F("fITSsignalPions",";p [GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of p for pi+
  fTPCsignalpip=new TH2F("fTPCsignalPions",";p [GeV/c];signal",kPBins,binsPDummy,300,-1,1);//TPC PID signal as function of p for pi+
  fTOFsignalpip=new TH2F("fTOFsignalPions",";p [GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of p for pi+
  fOutputList->Add(fITSsignalpip);
  fOutputList->Add(fTPCsignalpip);
  fOutputList->Add(fTOFsignalpip);

  fITSsignalKp=new TH2F("fITSsignalKaons",";p [GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of p for K+
  fTPCsignalKp=new TH2F("fTPCsignalKaons",";p [GeV/c];signal",kPBins,binsPDummy,300,-1,1);//TPC PID signal as function of p for K+
  fTOFsignalKp=new TH2F("fTOFsignalKaons",";p [GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of p for K+
  fOutputList->Add(fITSsignalKp);
  fOutputList->Add(fTPCsignalKp);
  fOutputList->Add(fTOFsignalKp);

  fITSsignalpp=new TH2F("fITSsignalProtons",";p [GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of p for p
  fTPCsignalpp=new TH2F("fTPCsignalProtons",";p [GeV/c];signal",kPBins,binsPDummy,300,-1,1);//TPC PID signal as function of p for p
  fTOFsignalpp=new TH2F("fTOFsignalProtons",";p [GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of p for p
  fOutputList->Add(fITSsignalpp);
  fOutputList->Add(fTPCsignalpp);
  fOutputList->Add(fTOFsignalpp);

  fTOFsignalBeta=new TH2F("fTOFbeta",";p[GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFsignalPiBeta=new TH2F("fTOFbetaPions",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFsignalKBeta=new TH2F("fTOFbetaKaons",";p [GeV/c];#beta-#beta_{K}",kPBins,binsPDummy,500, -0.25, 0.2);//
  fTOFsignalPBeta=new TH2F("fTOFbetaProtons",";p [GeV/c];#beta-#beta_{p}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fOutputList->Add(fTOFsignalBeta);
  fOutputList->Add(fTOFsignalPiBeta);
  fOutputList->Add(fTOFsignalKBeta);
  fOutputList->Add(fTOFsignalPBeta);

  if(fMC)
  {
    fITSsignalpiMCp=new TH2F("fITSsignalPionsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of pt for pi+
    fTPCsignalpiMCp=new TH2F("fTPCsignalPionsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-1,1);//TPC PID signal as function of pt for pi+
    fTOFsignalpiMCp=new TH2F("fTOFsignalPionsMC",";p [GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt for pi+
    fOutputList->Add(fITSsignalpiMCp);
    fOutputList->Add(fTPCsignalpiMCp);
    fOutputList->Add(fTOFsignalpiMCp);

    fITSsignalKMCp=new TH2F("fITSsignalKaonsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of pt for K+
    fTPCsignalKMCp=new TH2F("fTPCsignalKaonsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-1,1);//TPC PID signal as function of pt for K+
    fTOFsignalKMCp=new TH2F("fTOFsignalKaonsMC",";p [GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt for K+
    fOutputList->Add(fITSsignalKMCp);
    fOutputList->Add(fTPCsignalKMCp);
    fOutputList->Add(fTOFsignalKMCp);

    fITSsignalpMCp=new TH2F("fITSsignalProtonsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of pt for p
    fTPCsignalpMCp=new TH2F("fTPCsignalProtonsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-1,1);//TPC PID signal as function of pt for p
    fTOFsignalpMCp=new TH2F("fTOFsignalProtonsMC",";p [GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt for p
    fOutputList->Add(fITSsignalpMCp);
    fOutputList->Add(fTPCsignalpMCp);
    fOutputList->Add(fTOFsignalpMCp);
  }

  fPvsPt=new TH2F("fPvsPt","p vs p_{t};p [GeV/c];p_{t} [GeV/c]",kPBins,binsPDummy,kPtBins,binsPtDummy);
  fOutputList->Add(fPvsPt);

  fMeanPvsP = new TProfile("fMeanPvsP","Mean P vs P;p [Gev/c];<p> [GeV/c]",kPBins,binsPDummy);
  fOutputList->Add(fMeanPvsP);

  fTPCvsGlobalMult = new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",500,0,2500,500,0,3500);
  fOutputList->Add(fTPCvsGlobalMult);

  fStandardGlobalCuts = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
  fStandardTPCCuts = AliFlowTrackCuts::GetStandardTPCOnlyTrackCuts2010();

  //fOutputList->Add(fESDpid);

  PostData(1,  fOutputList);
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::UserExec(Option_t *)
{
  fESD = dynamic_cast<AliESDEvent*> (InputEvent());
  if (!fESD) return;

  //do the calibration bit
  fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0); // to use T0-TOF 
  fESDpid->MakePID(fESD,kFALSE);

  if(!fCuts || !fEventCuts)
  {
    Printf("No CUTS Defined.........\n");
    PostData(1,  fOutputList);
    return;
  }

  if (!(fEventCuts->IsSelected(fESD)))
  {
    return;
  }

  AliStack* stack=0x0;
  if(fMC)
  {
    AliMCEvent* mcEvent  = (AliMCEvent*) MCEvent();
    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
    stack = mcEvent->Stack();
  }

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  Int_t nTracks=fESD->GetNumberOfTracks();

  AliESDtrack *trackESD=0;

  for(int tr1=0; tr1<nTracks; tr1++)
  {
    trackESD=fESD->GetTrack(tr1);
    if (!trackESD) continue;

    if(!(fCuts->AcceptTrack(trackESD))) continue;

    Int_t label=-1;
    if(fMC) label=trackESD->GetLabel();

    Int_t pdgcode=0;
    if(stack&&fMC)
    {
      TParticle* particle2 = stack->Particle(TMath::Abs(label));
      pdgcode=particle2->GetPdgCode();
    }

    Double_t p=trackESD->GetP();
    Double_t pt=trackESD->Pt();
    fPvsPt->Fill(p,pt);
    fMeanPvsP->Fill(p,p);

    pidITS(trackESD,pdgcode);
    pidTPC(trackESD,pdgcode);
    pidTOF(trackESD,pdgcode);
  }

  //check the correlation between the global and TPConly number of tracks
  fStandardGlobalCuts->SetEvent(fESD);
  fStandardTPCCuts->SetEvent(fESD);
  fTPCvsGlobalMult->Fill(fStandardGlobalCuts->Count(),fStandardTPCCuts->Count());

  // Post output data.
  PostData(1,  fOutputList);
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::Terminate(Option_t *)
{
  //Terminate
  if(fCuts)
    fCuts->Dump();
  if(fMC)
    Printf("MC On\n");

  Printf("AliAnalysisTaskPIDflowQA: end of Terminate");
}

//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidITS(AliESDtrack* t, Int_t pdgcode)
{
  Int_t ngoodSDDSSD=0;
  Double_t sample[4]= {0.0,0.0,0.0,0.0};
  t->GetITSdEdxSamples(sample);
  for(int i=0; i<4; i++)
  {
    if(sample[i]>50.0)
      ngoodSDDSSD++;
  }
  if(ngoodSDDSSD<3) return;

  Float_t dedx=(Float_t)t->GetITSsignal();
  if(dedx<1.0) return;

  Bool_t ifSA=!(t->IsOn(AliESDtrack::kTPCin));
  Float_t p=t->GetP();

  Float_t signalpi=fESDpid->GetITSResponse().Bethe(p,AliPID::ParticleMass(2),ifSA);
  Float_t signalK=fESDpid->GetITSResponse().Bethe(p,AliPID::ParticleMass(3),ifSA);
  Float_t signalp=fESDpid->GetITSResponse().Bethe(p,AliPID::ParticleMass(4),ifSA);
  if(signalpi<1.0||signalK<1.0||signalp<1.0) return;

  fITSsignal->Fill(p,dedx);

  fITSsignalpip->Fill(p,TMath::Log(dedx)-TMath::Log(signalpi));
  fITSsignalKp->Fill(p,TMath::Log(dedx)-TMath::Log(signalK));
  fITSsignalpp->Fill(p,TMath::Log(dedx)-TMath::Log(signalp));

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fITSsignalpiMCp->Fill(p,TMath::Log(dedx)-TMath::Log(signalpi));
    else if(TMath::Abs(pdgcode)==321)
      fITSsignalKMCp->Fill(p,TMath::Log(dedx)-TMath::Log(signalK));
    else if (TMath::Abs(pdgcode)==2212)
      fITSsignalpMCp->Fill(p,TMath::Log(dedx)-TMath::Log(signalp));
  }
}

//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTPC(AliESDtrack* t, Int_t pdgcode)
{
  Double_t pinTPCglobal=t->GetInnerParam()->GetP();
  Float_t sigPion     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kPion);
  Float_t sigKaon     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kKaon);
  Float_t sigProton   = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kProton);
  Float_t p=t->GetP();
  Double_t tpcSignal =t ->GetTPCsignal();
  if(!(sigPion>0.0&&sigKaon>0.0&&sigProton>0.0))
    return;

  fTPCsignal->Fill(pinTPCglobal,tpcSignal);

  fTPCsignalpip->Fill(p,(tpcSignal-sigPion)/sigPion);
  fTPCsignalKp->Fill(p,(tpcSignal-sigKaon)/sigKaon);
  fTPCsignalpp->Fill(p,(tpcSignal-sigProton)/sigProton);

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTPCsignalpiMCp->Fill(p,(tpcSignal-sigPion)/sigPion);
    else if(TMath::Abs(pdgcode)==321)
      fTPCsignalKMCp->Fill(p,(tpcSignal-sigKaon)/sigKaon);
    else if (TMath::Abs(pdgcode)==2212)
      fTPCsignalpMCp->Fill(p,(tpcSignal-sigProton)/sigProton);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTOF(AliESDtrack* t, Int_t pdgcode)
{
  Bool_t goodtrack = (t) &&
                     (t->GetStatus() & AliESDtrack::kTOFpid) &&
                     (t->GetTOFsignal() > 12000) &&
                     (t->GetTOFsignal() < 100000) &&
                     (t->GetIntegratedLength() > 365) &&
                     !(t->GetStatus() & AliESDtrack::kTOFmismatch);

  if (!goodtrack) return;

  const Float_t c = 2.99792457999999984e-02;
  Float_t p=t->GetP();
  Float_t L = t->GetIntegratedLength();
  Float_t fT0track=fESDpid->GetTOFResponse().GetStartTime(p);
  Float_t timeTOF=t->GetTOFsignal()- fT0track;

  //calculate beta for the track
  Float_t beta = L/timeTOF/c;
  
  //2=pion 3=kaon 4=protons
  Double_t inttimes[5]= {-1.0,-1.0,-1.0,-1.0,-1.0};
  t->GetIntegratedTimes(inttimes);
  Float_t betaHypothesis[5] = {0.0,0.0,0.0,0.0,0.0};
  for (Int_t i=0;i<5;i++)
  {
    betaHypothesis[i] = L/inttimes[i]/c;
  }

  fTOFsignal->Fill(p,timeTOF-inttimes[2]);

  //beta part
  fTOFsignalBeta->Fill(p,beta);
  fTOFsignalPiBeta->Fill(p,beta-betaHypothesis[2]);
  fTOFsignalKBeta->Fill(p,beta-betaHypothesis[3]);
  fTOFsignalPBeta->Fill(p,beta-betaHypothesis[4]);

  //P part
  fTOFsignalpip->Fill(p,timeTOF-inttimes[2]);
  fTOFsignalKp->Fill(p,timeTOF-inttimes[3]);
  fTOFsignalpp->Fill(p,timeTOF-inttimes[4]);

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTOFsignalpiMCp->Fill(p,timeTOF-inttimes[2]);
    else if(TMath::Abs(pdgcode)==321)
      fTOFsignalKMCp->Fill(p,timeTOF-inttimes[3]);
    else if (TMath::Abs(pdgcode)==2212)
      fTOFsignalpMCp->Fill(p,timeTOF-inttimes[4]);
  }
}

Float_t AliAnalysisTaskPIDflowQA::Beta(Float_t m, Float_t p) 
{
  //get theoretical beta
  return TMath::Sqrt(1. / (1. + m * m / (p * p)));
}
 
