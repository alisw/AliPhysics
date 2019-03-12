#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include <TROOT.h>

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliCentrality.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskCTauPbPbaod.h"

//extern TROOT *gROOT;

ClassImp(AliAnalysisTaskCTauPbPbaod)


//
//  This is a little task for checking the c*tau of the strange particles 
//

AliAnalysisTaskCTauPbPbaod::AliAnalysisTaskCTauPbPbaod(const char *name) :
AliAnalysisTaskSE(name),
fIsMC(kFALSE),
fCMin(0.),
fCMax(90.),
fCPA(0.9975),
fDCA(1.0),
fTPCcr(70.),
fTPCcrfd(0.8),
fDCApv(0.1),
fRmin(0.9),
fRmax(100.),

fOutput(0),
fMult(0),
fCosPA(0),
fDtrDCA(0),
fTPCrows(0),
fTPCratio(0),
fPrimDCA(0),
fR(0),
fdEdx(0),
fdEdxPid(0),

fK0sM(0),
fK0sSi(0),
fK0sMC(0),
fK0sAs(0),

fLambdaM(0),
fLambdaSi(0),
fLambdaMC(0),
fLambdaAs(0),

fLambdaEff(0),
fLambdaPt(0),

fLambdaBarM(0),
fLambdaBarSi(0),
fLambdaBarMC(0),
fLambdaBarAs(0),

fLambdaBarEff(0),
fLambdaBarPt(0),

fLambdaFromXi(0),
fXiM(0),
fXiSiP(0),

fLambdaBarFromXiBar(0),
fXiBarM(0),
fXiBarSiP(0)
{
  // Constructor. Initialization of pointers
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskCTauPbPbaod::UserCreateOutputObjects()
{
  Int_t    nbins=100;  // number of bins
  Double_t ltMax=100.;
  Double_t ptMax=10.;

  fOutput = new TList(); 
  fOutput->SetOwner();


  fMult=new TH1F("fMult","Multiplicity",1100,0.,3300);
  fMult->GetXaxis()->SetTitle("N tracks"); 
  fOutput->Add(fMult);

  fCosPA=new TH1F("fCosPA","Cos(PA) distribution",50,0.9975,1.0005);
  fCosPA->GetXaxis()->SetTitle("Cos(PA)"); 
  fOutput->Add(fCosPA);

  fDtrDCA=new TH1F("fDtrDCA","DCA between V0 daughters",50,0.0,1.5);
  fDtrDCA->GetXaxis()->SetTitle("DCA (rel. u.)"); 
  fOutput->Add(fDtrDCA);

  fTPCrows=new TH1F("fTPCrows","TPC crossed pad rows",180,0.,180.);
  fTPCrows->GetXaxis()->SetTitle("TPC crossed pad rows"); 
  fOutput->Add(fTPCrows);

  fTPCratio=new TH1F("fTPCratio","TPC crossed/findable pad rows",50,0.0,1.5);
  fTPCratio->GetXaxis()->SetTitle("TPC crossed/findable pad rows"); 
  fOutput->Add(fTPCratio);

  fPrimDCA=new TH1F("fPrimDCA","DCA wrt the primary vertex",50,0.0,1.5);
  fPrimDCA->GetXaxis()->SetTitle("DCA wrt the PV (cm)"); 
  fOutput->Add(fPrimDCA);

  fR=new TH1F("fR","Radius of the V0 vertices",101,0.0,102);
  fR->GetXaxis()->SetTitle("R (cm)"); 
  fOutput->Add(fR);

  fdEdx=new TH2F("fdEdx","dE/dx",50,0.2,3,50,0.,6.);
  fOutput->Add(fdEdx);

  fdEdxPid=new TH2F("fdEdxPid","dE/dx with PID",50,0.2,3,50,0.,6.);
  fOutput->Add(fdEdxPid);

  fK0sM = 
  new TH2F("fK0sM", "Mass for K^{0}_{s}", nbins/2,0.448,0.548,nbins,0.,ptMax);
  fK0sM->GetXaxis()->SetTitle("Mass (GeV/c)"); 
  fOutput->Add(fK0sM);

  fK0sSi = 
  new TH2F("fK0sSi","L_{T} vs p_{T} for K^{0}_{s}, side-band subtracted",
  nbins,0.,ptMax,nbins,0.,ltMax);
  fK0sSi->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sSi->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fK0sSi);

  fK0sMC = 
  new TH2F("fK0sMC","L_{T} vs p_{T} for K^{0}_{s}, from MC stack", 
  nbins,0.,ptMax,nbins,0.,ltMax);
  fK0sMC->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sMC->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fK0sMC);

  fK0sAs = 
  new TH2F("fK0sAs", "L_{T} vs p_{T} for K^{0}_{s}, associated", 
  nbins,0.,ptMax,nbins,0.,ltMax);
  fK0sAs->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fK0sAs->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fK0sAs);

  //----------------------

  fLambdaM = 
  new TH2F("fLambdaM","Mass for \\Lambda", nbins, 1.065, 1.165,nbins,0.,ptMax);
  fLambdaM->GetXaxis()->SetTitle("Mass (GeV/c)"); 
  fOutput->Add(fLambdaM);

  fLambdaSi = 
  new TH2F("fLambdaSi","L_{T} vs p_{T} for \\Lambda, side-band subtructed",
  nbins,0.,ptMax,nbins,0.,ltMax);
  fLambdaSi->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaSi->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaSi);

  fLambdaMC = 
  new TH2F("fLambdaMC","L_{T} vs p_{T} for \\Lambda, from MC stack", 
  nbins,0.,ptMax,nbins,0.,ltMax);
  fLambdaMC->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaMC->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaMC);

  fLambdaAs = 
  new TH2F("fLambdaAs","L_{T} vs p_{T} for \\Lambda, associated",
  nbins,0.,ptMax,nbins,0.,ltMax);
  fLambdaAs->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaAs->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaAs);

  //----------------------

  fLambdaEff=fLambdaAs->ProjectionX();
  fLambdaEff->SetName("fLambdaEff");
  fLambdaEff->SetTitle("Efficiency for #Lambda");
  fOutput->Add(fLambdaEff);

  fLambdaPt=fLambdaAs->ProjectionX();
  fLambdaPt->SetName("fLambdaPt");
  fLambdaPt->SetTitle("Raw #Lambda pT spectrum");
  fOutput->Add(fLambdaPt);

  //----------------------

  fLambdaBarM = 
  new TH2F("fLambdaBarM","Mass for anti-\\Lambda", nbins, 1.065, 1.165,nbins,0.,ptMax);
  fLambdaBarM->GetXaxis()->SetTitle("Mass (GeV/c)"); 
  fOutput->Add(fLambdaBarM);

  fLambdaBarSi = 
  new TH2F("fLambdaBarSi","L_{T} vs p_{T} for anti-\\Lambda, side-band subtructed",
  nbins,0.,ptMax,nbins,0.,ltMax);
  fLambdaBarSi->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaBarSi->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaBarSi);

  fLambdaBarMC = 
  new TH2F("fLambdaBarMC","L_{T} vs p_{T} for anti-\\Lambda, from MC stack", 
  nbins,0.,ptMax,nbins,0.,ltMax);
  fLambdaBarMC->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaBarMC->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaBarMC);

  fLambdaBarAs = 
  new TH2F("fLambdaBarAs","L_{T} vs p_{T} for anti-\\Lambda, associated",
  nbins,0.,ptMax,nbins,0.,ltMax);
  fLambdaBarAs->GetXaxis()->SetTitle("p_{T} (GeV/c)"); 
  fLambdaBarAs->GetYaxis()->SetTitle("L_{T} (cm)"); 
  fOutput->Add(fLambdaBarAs);


  //----------------------

  fLambdaBarEff=fLambdaBarAs->ProjectionX();
  fLambdaBarEff->SetName("fLambdaBarEff");
  fLambdaBarEff->SetTitle("Efficiency for anti-#Lambda");
  fOutput->Add(fLambdaBarEff);

  fLambdaBarPt=fLambdaBarAs->ProjectionX();
  fLambdaBarPt->SetName("fLambdaBarPt");
  fLambdaBarPt->SetTitle("Raw anti-#Lambda pT spectrum");
  fOutput->Add(fLambdaBarPt);

  //----------------------

  fLambdaFromXi=new TH3F("fLambdaFromXi","L_{T} vs p_{T} vs p_{T} of \\Xi for \\Lambda from \\Xi",
  nbins,0.,ptMax,nbins,0.,ltMax,33,0.,ptMax+2);
  fOutput->Add(fLambdaFromXi);

  fXiM  = 
  new TH2F("fXiM", "\\Xi mass distribution", 50, 1.271, 1.371,33,0.,ptMax+2);
  fOutput->Add(fXiM);

  fXiSiP  = new TH1F("fXiSiP", "Pt for \\Xi, side-band subracted",
  33,0.,ptMax+2);
  fOutput->Add(fXiSiP);


  fLambdaBarFromXiBar=new TH3F("fLambdaBarFromXiBar","L_{T} vs p_{T} vs p_{T} of anti-\\Xi for anti-\\Lambda from anti-\\Xi",
  nbins,0.,ptMax,nbins,0.,ltMax,33,0.,ptMax+2);
  fOutput->Add(fLambdaBarFromXiBar);

  fXiBarM  = 
  new TH2F("fXiBarM", "anti-\\Xi mass distribution", 50, 1.271, 1.371,33,0.,ptMax+2);
  fOutput->Add(fXiBarM);

  fXiBarSiP  = new TH1F("fXiBarSiP", "Pt for anti-\\Xi, side-band subracted",
  33,0.,ptMax+2);
  fOutput->Add(fXiBarSiP);


  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskCTauPbPbaod::AcceptTrack(const AliAODTrack *t) {
  if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  //if (t->GetKinkIndex(0)>0) return kFALSE;

  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  if (nCrossedRowsTPC < fTPCcr) return kFALSE;
  Int_t findable=t->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC/findable < fTPCcrfd) return kFALSE;

  if (TMath::Abs(t->Eta()) > 0.8) return kFALSE;

  fTPCrows->Fill(nCrossedRowsTPC);
  fTPCratio->Fill(nCrossedRowsTPC/findable);

  return kTRUE;   
}

Bool_t 
AliAnalysisTaskCTauPbPbaod::AcceptV0(const AliAODv0 *v0,const AliAODEvent *aod)
{
  if (v0->GetOnFlyStatus()) return kFALSE;

  Double_t cpa=v0->CosPointingAngle(aod->GetPrimaryVertex());
  if (cpa < fCPA) return kFALSE;

  Double_t dca=v0->DcaV0Daughters();
  if (dca > fDCA) return kFALSE;

  const AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
  if (!AcceptTrack(ntrack)) return kFALSE;

  const AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);
  if (!AcceptTrack(ptrack)) return kFALSE;

  Float_t xyn=v0->DcaNegToPrimVertex();
  if (TMath::Abs(xyn)<fDCApv) return kFALSE;
  Float_t xyp=v0->DcaPosToPrimVertex();
  if (TMath::Abs(xyp)<fDCApv) return kFALSE;

  Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
  Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
  if (r2<fRmin*fRmin) return kFALSE;
  if (r2>fRmax*fRmax) return kFALSE;

  fCosPA->Fill(cpa);
  fDtrDCA->Fill(dca);
  fPrimDCA->Fill(xyn); fPrimDCA->Fill(xyp);
  fR->Fill(TMath::Sqrt(r2));

  return kTRUE;
}

Bool_t AliAnalysisTaskCTauPbPbaod::AcceptCascade(const AliAODcascade *cs, const AliAODEvent *aod) {

  AliAODVertex *xi = cs->GetDecayVertexXi(); 
  const AliAODTrack *bach=(AliAODTrack *)xi->GetDaughter(0);
  if (!AcceptTrack(bach)) return kFALSE;
   
  Double_t xy=cs->DcaBachToPrimVertex();
  if (TMath::Abs(xy) < 0.03) return kFALSE;

  const AliAODVertex *vtx=aod->GetPrimaryVertex();
  Double_t xv=vtx->GetX(), yv=vtx->GetY(), zv=vtx->GetZ();
  Double_t cpa=cs->CosPointingAngleXi(xv,yv,zv);
  if (cpa<0.999) return kFALSE;

  if (cs->DcaXiDaughters() > 0.3) return kFALSE;

  return kTRUE;
}

static Bool_t AcceptPID(const AliPIDResponse *pidResponse, 
			const AliAODTrack *ptrack, const TClonesArray *stack) {

  const AliAODPid *pid=ptrack->GetDetPid();
  if (!pid) return kTRUE;
  if (pid->GetTPCmomentum() > 1.) return kTRUE;

  if (stack) {
    // MC PID
    Int_t plab=TMath::Abs(ptrack->GetLabel());
    AliAODMCParticle *pp=(AliAODMCParticle*)(*stack)[plab];
    if (!pp) return kTRUE;
    if (pp->Charge() > 0) {
       if (pp->GetPdgCode() == kProton)    return kTRUE;
    } else {
       if (pp->GetPdgCode() == kProtonBar) return kTRUE;
    }
  } else {
    // Real PID
    Double_t nsig=pidResponse->NumberOfSigmasTPC(ptrack,AliPID::kProton);
    if (TMath::Abs(nsig) < 3.) return kTRUE;
  }
  
  return kFALSE; 
}

static AliAODMCParticle*
AssociateV0(const AliAODTrack *ptrack, const AliAODTrack *ntrack,
const TClonesArray *stack, AliAODMCParticle *&mcp) {
  //
  // Try to associate the V0 with the daughters ptrack and ntrack
  // with the Monte Carlo
  //
  if (!stack) return 0;

  Int_t nlab=TMath::Abs(ntrack->GetLabel());
  AliAODMCParticle *n=(AliAODMCParticle*)(*stack)[nlab];
  if (!n) return 0;

  Int_t plab=TMath::Abs(ptrack->GetLabel());
  AliAODMCParticle *p=(AliAODMCParticle*)(*stack)[plab];
  if (!p) return 0;

  Int_t imp=p->GetMother(); //V0 particle, the mother of the pos. track 
  AliAODMCParticle *p0=(AliAODMCParticle*)(*stack)[imp];
  if (!p0) return 0;

  Int_t imn=n->GetMother();
  if (imp != imn) {  // Check decays of the daughters
     return 0; // Fixme
  }

  Int_t code=p0->GetPdgCode();
  if (code != kK0Short)
     if (code != kLambda0)
	if (code != kLambda0Bar) return 0;

  mcp=p;

  return p0;
}


void AliAnalysisTaskCTauPbPbaod::UserExec(Option_t *)
{
  const Double_t yMax=0.5;
  const Double_t pMin=0.0;
  const Double_t lMax=0.001;

  AliAODEvent *aod=(AliAODEvent *)InputEvent();

  if (!aod) {
    Printf("ERROR: aod not available");
    return;
  }

  // Vertex selection
  const AliAODVertex *vtx=aod->GetPrimaryVertex();
  if (vtx->GetNContributors()<3) return;
  Double_t xv=vtx->GetX(), yv=vtx->GetY(), zv=vtx->GetZ();

  if (TMath::Abs(zv) > 10.) return ;   
 

  // Physics selection
  AliAnalysisManager *mgr= AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *hdr=(AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t maskIsSelected = hdr->IsEventSelected();
  Bool_t isSelected = (maskIsSelected & AliVEvent::kMB);
  if (!isSelected) return;

  fMult->Fill(-100); //event counter  

  // Centrality selection
  AliCentrality *cent=aod->GetCentrality();
  if (!cent->IsEventInCentralityClass(fCMin,fCMax,"V0M")) return;

  AliPIDResponse *pidResponse = hdr->GetPIDResponse(); 
 


  //+++++++ MC
  TClonesArray *stack = 0x0;
  Double_t mcXv=0., mcYv=0., mcZv=0.;
   
  if (fIsMC) {
     TList *lst = aod->GetList();
     stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
     if (!stack) {
        Printf("ERROR: stack not available");
        return;
     }
     AliAODMCHeader *
     mcHdr=(AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());

     mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();

     Int_t ntrk=stack->GetEntriesFast();
     while (ntrk--) {
       AliAODMCParticle *p0=(AliAODMCParticle*)stack->UncheckedAt(ntrk);
       Int_t code=p0->GetPdgCode();
       if (code != kK0Short)
	 if (code != kLambda0)
	    if (code != kLambda0Bar) continue;

       Int_t plab=p0->GetDaughterLabel(0), nlab=p0->GetDaughterLabel(1);
       if (nlab==plab) continue;
       AliAODMCParticle *part=(AliAODMCParticle*)(*stack)[plab];
       if (!part) continue;
       Double_t charge=part->Charge();
       if (charge == 0.) continue;

       Double_t pt=p0->Pt();
       if (pt<pMin) continue;
       if (TMath::Abs(p0->Y())>yMax) continue;
    
       Double_t x=p0->Xv(), y=p0->Yv(), z=p0->Zv();
       Double_t dxmc=mcXv-x, dymc=mcYv-y, dzmc=mcZv-z;
       Double_t l=TMath::Sqrt(dxmc*dxmc + dymc*dymc + dzmc*dzmc);

       if (l > lMax) continue; // secondary V0

       x=part->Xv(); y=part->Yv();
       dxmc=mcXv-x; dymc=mcYv-y;
       Double_t lt=TMath::Sqrt(dxmc*dxmc + dymc*dymc);

       switch (code) {
       case kK0Short:
          fK0sMC->Fill(pt,lt);
          break;
       case kLambda0:
          fLambdaMC->Fill(pt,lt);
          break;
       case kLambda0Bar:
          fLambdaBarMC->Fill(pt,lt);
          break;
       default: break;
       }
     }
  }
  //+++++++


  Int_t ntrk1=aod->GetNumberOfTracks();
  Int_t mult=0;
  for (Int_t i=0; i<ntrk1; i++) {
    AliAODTrack *t=dynamic_cast<AliAODTrack*>(aod->GetTrack(i));
    if(!t) { AliFatal("Not a standard AOD"); return; }
    if (t->IsMuonTrack()) continue;
    if (!t->IsOn(AliAODTrack::kTPCrefit)) continue;

    Double_t xyz[3];
    if (t->GetPosition(xyz)) continue;
    if (TMath::Abs(xyz[0])>3.) continue;
    if (TMath::Abs(xyz[1])>3.) continue;

    Double_t pt=t->Pt(),pz=t->Pz();
    if (TMath::Abs(pz/pt)>0.8) continue;

    mult++;

    const AliAODPid *pid=t->GetDetPid();
    if (!pid) continue;

    Double_t p=pid->GetTPCmomentum();
    Double_t dedx=pid->GetTPCsignal()/47.;
    fdEdx->Fill(p,dedx,1);

    Double_t nsig=pidResponse->NumberOfSigmasTPC(t,AliPID::kProton);
    if (TMath::Abs(nsig) < 3.) fdEdxPid->Fill(p,dedx,1);

  }
  fMult->Fill(mult);


  Int_t nv0 = aod->GetNumberOfV0s();
  while (nv0--) {
      AliAODv0 *v0=aod->GetV0(nv0);

      Double_t pt=TMath::Sqrt(v0->Pt2V0());
      if (pt < pMin) continue;

      if (!AcceptV0(v0,aod)) continue;
      
      const AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
      const AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);

      Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
      Double_t dxv=xyz[0]-xv, dyv=xyz[1]-yv;
      Double_t lt=TMath::Sqrt(dxv*dxv + dyv*dyv);

      if (lt/pt > 3*7.89/1.1157) continue;  

      //--- V0 switches
      Bool_t isK0s=kTRUE;
      Bool_t isLambda=kTRUE;
      Bool_t isLambdaBar=kTRUE;

      if (0.4977*lt/pt > 3*2.68) isK0s=kFALSE;
      if (1.1157*lt/pt > 3*7.89) isLambdaBar=isLambda=kFALSE;

      if (!AcceptPID(pidResponse, ptrack, stack)) isLambda=kFALSE;
      if (!AcceptPID(pidResponse, ntrack, stack)) isLambdaBar=kFALSE;

      Double_t yK0s=TMath::Abs(v0->RapK0Short());
      Double_t yLam=TMath::Abs(v0->RapLambda());
      if (yK0s > yMax) isK0s=kFALSE;
      if (yLam > yMax) isLambda=isLambdaBar=kFALSE;
      //---

      Double_t mass=0., m=0., s=0.;
      if (isK0s) {
         mass=v0->MassK0Short();
         fK0sM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
         s=0.0044 + (0.008-0.0044)/(10-1)*(pt - 1.);
         if (TMath::Abs(m-mass) < 3*s) {
            fK0sSi->Fill(pt,lt);
         } else {
	    isK0s=kFALSE;
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fK0sSi->Fill(pt,lt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fK0sSi->Fill(pt,lt,-1);
         }
      }
      
      if (isLambda) {
         mass=v0->MassLambda();
         fLambdaM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
         //s=0.0027 + (0.004-0.0027)/(10-1)*(pt-1);
         //s=0.0015 + (0.002-0.0015)/(2.6-1)*(pt-1);
         s=0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
         if (TMath::Abs(m-mass) < 3*s) {
            fLambdaSi->Fill(pt,lt);
         } else {
	    isLambda=kFALSE;
         } 
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fLambdaSi->Fill(pt,lt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fLambdaSi->Fill(pt,lt,-1);
         }
      }

      if (isLambdaBar) {
         mass=v0->MassAntiLambda();
         fLambdaBarM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kLambda0Bar)->Mass();
         //s=0.0027 + (0.004-0.0027)/(10-1)*(pt-1);
         //s=0.0015 + (0.002-0.0015)/(2.6-1)*(pt-1);
         s=0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
         if (TMath::Abs(m-mass) < 3*s) {
            fLambdaBarSi->Fill(pt,lt);
         } else {
	    isLambdaBar=kFALSE;
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fLambdaBarSi->Fill(pt,lt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fLambdaBarSi->Fill(pt,lt,-1);
         }
      }

      if (!fIsMC) continue;

      //++++++ MC 
      if (!isK0s)
         if (!isLambda)
             if (!isLambdaBar) continue;//check MC only for the accepted V0s 

      AliAODMCParticle *mcp=0;
      AliAODMCParticle *mc0=AssociateV0(ptrack,ntrack,stack,mcp);
      if (!mc0) continue;

      Double_t ptAs=mc0->Pt();
      if (ptAs < pMin) continue;
      Double_t yAs=mc0->Y();
      if (TMath::Abs(yAs) > yMax) continue;

      Int_t code=mc0->GetPdgCode();

      Double_t dx = mcXv - mcp->Xv(), dy = mcYv - mcp->Yv();
      Double_t ltAs=TMath::Sqrt(dx*dx + dy*dy);
 
      Double_t dz=mcZv - mc0->Zv(); dy=mcYv - mc0->Yv(); dx=mcXv - mc0->Xv();
      Double_t l = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
      if (l<lMax) { // Primary V0s
	 switch (code) {
         case kK0Short:
            if (isK0s)       fK0sAs->Fill(ptAs,ltAs);
            break;
	 case kLambda0:
            if (isLambda)    fLambdaAs->Fill(ptAs,ltAs);
            break;
	 case kLambda0Bar:
            if (isLambdaBar) fLambdaBarAs->Fill(ptAs,ltAs);
            break;
         default: break;
	 }
      } else {
	 if (code==kK0Short) continue;

         Int_t nx=mc0->GetMother();
         AliAODMCParticle *xi=(AliAODMCParticle*)(*stack)[nx];
         if (!xi) continue;
         Int_t xcode=xi->GetPdgCode();
         
	 switch (code) {
	 case kLambda0:
            if (isLambda)
	    if ((xcode==kXiMinus) || (xcode==3322)) 
               fLambdaFromXi->Fill(ptAs,ltAs,xi->Pt());
            break;
	 case kLambda0Bar:
            if (isLambdaBar)
	    if ((xcode==kXiPlusBar)||(xcode==-3322)) 
               fLambdaBarFromXiBar->Fill(ptAs,ltAs,xi->Pt());
            break;
         default: break;
	 }
      }
      //++++++
  
  }

  Int_t ncs=aod->GetNumberOfCascades();
  for (Int_t i=0; i<ncs; i++) {
      AliAODcascade *cs=aod->GetCascade(i);
 
      Double_t pt=TMath::Sqrt(cs->Pt2Xi());
      if (pt < pMin) continue;
      if (TMath::Abs(cs->RapXi()) > yMax) continue;
      if (!AcceptCascade(cs,aod)) continue;

      const AliAODv0 *v0=(AliAODv0*)cs;
      if (v0->RapLambda() > yMax) continue;
      if (!AcceptV0(v0,aod)) continue;

      //--- Cascade switches
      Bool_t isXiMinus=kTRUE;
      Bool_t isXiPlusBar=kTRUE;

      const AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);
      if (!AcceptPID(pidResponse, ptrack, stack)) isXiMinus=kFALSE;

      const AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
      if (!AcceptPID(pidResponse, ntrack, stack)) isXiPlusBar=kFALSE;

      Int_t charge=cs->ChargeXi();      
      if (charge > 0) isXiMinus=kFALSE;
      if (charge < 0) isXiPlusBar=kFALSE;
      //---

      if (isXiMinus) {
         Double_t mass=cs->MassXi();
         fXiM->Fill(mass,pt);
         Double_t m=TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass();
         //Double_t s=0.0037;
         Double_t s=0.002 + (0.0032-0.002)/(6-1.5)*(pt-1.5);
         if (TMath::Abs(m-mass) < 3*s) {
            fXiSiP->Fill(pt);
         } else {
            isXiMinus=kFALSE;
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fXiSiP->Fill(pt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fXiSiP->Fill(pt,-1);
         }
      }

      if (isXiPlusBar) {
         Double_t mass=cs->MassXi();
         fXiBarM->Fill(mass,pt);
         Double_t m=TDatabasePDG::Instance()->GetParticle(kXiPlusBar)->Mass();
         //Double_t s=0.0037;
         Double_t s=0.002 + (0.0032-0.002)/(6-1.5)*(pt-1.5);
         if (TMath::Abs(m-mass) < 3*s) {
            fXiBarSiP->Fill(pt);
         } else {
            isXiPlusBar=kFALSE; 
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fXiBarSiP->Fill(pt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fXiBarSiP->Fill(pt,-1);
         }
      }

      if (!fIsMC) continue;

      //++++++ MC 
      if (!isXiMinus)
         if (!isXiPlusBar) continue;//check MC only for the accepted cascades 
      // Here is the future association with MC

  }

}

void AliAnalysisTaskCTauPbPbaod::Terminate(Option_t *)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  fOutput=(TList*)GetOutputData(1);
  if (!fOutput) {
     Printf("ERROR: fOutput not available");
     return;
  }
 
  fMult = dynamic_cast<TH1F*>(fOutput->FindObject("fMult")) ; 
  if (!fMult) {
     Printf("ERROR: fMult not available");
     return;
  }

  fdEdx = dynamic_cast<TH2F*>(fOutput->FindObject("fdEdx")) ; 
  if (!fdEdx) {
     Printf("ERROR: fdEdx not available");
     return;
  }

  fdEdxPid = dynamic_cast<TH2F*>(fOutput->FindObject("fdEdxPid")) ; 
  if (!fdEdxPid) {
     Printf("ERROR: fdEdxPid not available");
     return;
  }


  fK0sMC = dynamic_cast<TH2F*>(fOutput->FindObject("fK0sMC")) ; 
  if (!fK0sMC) {
     Printf("ERROR: fK0sMC not available");
     return;
  }
  TH1D *k0sMcPx=fK0sMC->ProjectionX(); k0sMcPx->Sumw2();
  fK0sAs = dynamic_cast<TH2F*>(fOutput->FindObject("fK0sAs")) ; 
  if (!fK0sAs) {
     Printf("ERROR: fK0sAs not available");
     return;
  }
  TH1D *k0sAsPx=fK0sAs->ProjectionX(); 
  k0sAsPx->Sumw2(); //k0sAsPx->Scale(0.69);



  // Lambda histograms 
  fLambdaFromXi = dynamic_cast<TH3F*>(fOutput->FindObject("fLambdaFromXi")) ; 
  if (!fLambdaFromXi) {
     Printf("ERROR: fLambdaFromXi not available");
     return;
  }
  TH1D *lambdaFromXiPx=fLambdaFromXi->ProjectionX(); lambdaFromXiPx->Sumw2();

  fLambdaMC = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaMC")) ; 
  if (!fLambdaMC) {
     Printf("ERROR: fLambdaMC not available");
     return;
  }
  TH1D *lambdaMcPx=fLambdaMC->ProjectionX(); lambdaMcPx->Sumw2();

  fLambdaAs = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaAs")) ; 
  if (!fLambdaAs) {
     Printf("ERROR: fLambdaAs not available");
     return;
  }
  TH1D *lambdaAsPx=fLambdaAs->ProjectionX(); 
  lambdaAsPx->Sumw2(); //lambdaAsPx->Scale(0.64);

  fLambdaSi = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaSi")) ; 
  if (!fLambdaSi) {
     Printf("ERROR: fLambdaSi not available");
     return;
  }
  TH1D *lambdaSiPx=fLambdaSi->ProjectionX(); 
  lambdaSiPx->SetName("fLambdaPt");
  lambdaSiPx->Sumw2();

  fLambdaEff = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaEff")) ; 
  if (!fLambdaEff) {
     Printf("ERROR: fLambdaEff not available");
     return;
  }
  fLambdaPt = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaPt")) ; 
  if (!fLambdaPt) {
     Printf("ERROR: fLambdaPt not available");
     return;
  }


  // anti-Lambda histograms 
  fLambdaBarFromXiBar = 
  dynamic_cast<TH3F*>(fOutput->FindObject("fLambdaBarFromXiBar")) ; 
  if (!fLambdaBarFromXiBar) {
     Printf("ERROR: fLambdaBarFromXiBar not available");
     return;
  }
  TH1D *lambdaBarFromXiBarPx=
  fLambdaBarFromXiBar->ProjectionX("Bar"); lambdaBarFromXiBarPx->Sumw2();

  fLambdaBarMC = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaBarMC")) ; 
  if (!fLambdaBarMC) {
     Printf("ERROR: fLambdaBarMC not available");
     return;
  }
  TH1D *lambdaBarMcPx=fLambdaBarMC->ProjectionX(); lambdaBarMcPx->Sumw2();

  fLambdaBarAs = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaBarAs")) ; 
  if (!fLambdaBarAs) {
     Printf("ERROR: fLambdaBarAs not available");
     return;
  }
  TH1D *lambdaBarAsPx=fLambdaBarAs->ProjectionX(); 
  lambdaBarAsPx->Sumw2(); //lambdaBarAsPx->Scale(0.64);

  fLambdaBarSi = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaBarSi")) ; 
  if (!fLambdaBarSi) {
     Printf("ERROR: fLambdaBarSi not available");
     return;
  }
  TH1D *lambdaBarSiPx=fLambdaBarSi->ProjectionX(); 
  lambdaBarSiPx->SetName("fLambdaBarPt");
  lambdaBarSiPx->Sumw2();

  fLambdaBarEff = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaBarEff")) ; 
  if (!fLambdaBarEff) {
     Printf("ERROR: fLambdaBarEff not available");
     return;
  }
  fLambdaBarPt = dynamic_cast<TH1D*>(fOutput->FindObject("fLambdaBarPt")) ; 
  if (!fLambdaBarPt) {
     Printf("ERROR: fLambdaBarPt not available");
     return;
  }


  if (!gROOT->IsBatch()) {

    TCanvas *c1 = new TCanvas("c1","Mulitplicity");
    c1->SetLogy();
    fMult->DrawCopy() ;

    new TCanvas("c2","dE/dx");
    fdEdx->DrawCopy() ;

    new TCanvas("c3","dE/dx with PID");
    fdEdxPid->DrawCopy() ;

    if (fIsMC) {
       /*
       TH1D effK(*k0sAsPx); effK.SetTitle("Efficiency for K0s");
       effK.Divide(k0sAsPx,k0sMcPx,1,1,"b");
       new TCanvas("c4","Efficiency for K0s");
       effK.DrawCopy("E") ;
       */

       //+++ Lambdas
       fLambdaEff->Divide(lambdaAsPx,lambdaMcPx,1,1,"b");
       new TCanvas("c5","Efficiency for #Lambda");
       fLambdaEff->DrawCopy("E") ;

       lambdaSiPx->Add(lambdaFromXiPx,-1);
       lambdaSiPx->Divide(fLambdaEff);

       new TCanvas("c6","Corrected #Lambda pt");
       lambdaSiPx->SetTitle("Corrected #Lambda pt");
      *fLambdaPt = *lambdaSiPx; 
       fLambdaPt->SetLineColor(2);
       fLambdaPt->DrawCopy("E");
    
       lambdaMcPx->DrawCopy("same");
 

       //+++ anti-Lambdas
       fLambdaBarEff->Divide(lambdaBarAsPx,lambdaBarMcPx,1,1,"b");
       new TCanvas("c7","Efficiency for anti-#Lambda");
       fLambdaBarEff->DrawCopy("E") ;

       lambdaBarSiPx->Add(lambdaBarFromXiBarPx,-1);
       lambdaBarSiPx->Divide(fLambdaBarEff);

       new TCanvas("c8","Corrected anti-#Lambda pt");
       lambdaBarSiPx->SetTitle("Corrected anti-#Lambda pt");
      *fLambdaBarPt = *lambdaBarSiPx; 
       fLambdaBarPt->SetLineColor(2);
       fLambdaBarPt->DrawCopy("E");
    
       lambdaBarMcPx->DrawCopy("same");
    } else {
       new TCanvas("c6","Raw #Lambda pt");
       lambdaSiPx->SetTitle("Raw #Lambda pt");
      *fLambdaPt = *lambdaSiPx; 
       fLambdaPt->SetLineColor(2);
       fLambdaPt->DrawCopy("E");


       new TCanvas("c7","Raw anti-#Lambda pt");
       lambdaBarSiPx->SetTitle("Raw anti-#Lambda pt");
      *fLambdaBarPt = *lambdaBarSiPx; 
       fLambdaBarPt->SetLineColor(2);
       fLambdaBarPt->DrawCopy("E");
    }
  }
}
