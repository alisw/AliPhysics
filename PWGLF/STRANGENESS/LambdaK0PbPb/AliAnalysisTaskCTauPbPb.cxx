#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticle.h>
#include <TROOT.h>

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"

#include "AliCentrality.h"

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskCTauPbPb.h"

//extern TROOT *gROOT;

ClassImp(AliAnalysisTaskCTauPbPb)


//
//  This is a little task for checking the c*tau of the strange particles 
//

AliAnalysisTaskCTauPbPb::AliAnalysisTaskCTauPbPb(const char *name) :
AliAnalysisTaskSE(name),
fIsMC(kFALSE),
fCMin(0.),
fCMax(90.),
fOutput(0),
fMult(0),
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

void AliAnalysisTaskCTauPbPb::UserCreateOutputObjects()
{
  Int_t    nbins=100;  // number of bins
  Double_t ltMax=100.;
  Double_t ptMax=10.;

  fOutput = new TList(); 
  fOutput->SetOwner();


  fMult=new TH1F("fMult","Multiplicity",1100,0.,3300);
  fMult->GetXaxis()->SetTitle("N tracks"); 
  fOutput->Add(fMult);

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

static Bool_t AcceptTrack(const AliESDtrack *t) {
  if (!t->IsOn(AliESDtrack::kTPCrefit)) return kFALSE;
  if (t->GetKinkIndex(0)>0) return kFALSE;

  //Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  //if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=t->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  //if (nCrossedRowsTPC/findable < 0.8) return kFALSE;

  if (TMath::Abs(t->Eta()) > 0.8) return kFALSE;

  return kTRUE;   
}

static Bool_t AcceptV0(const AliESDv0 *v0, const AliESDEvent *esd) {

  if (v0->GetOnFlyStatus()) return kFALSE;

  Int_t nidx=TMath::Abs(v0->GetNindex());
  AliESDtrack *ntrack=esd->GetTrack(nidx);
  if (!AcceptTrack(ntrack)) return kFALSE;

  Int_t pidx=TMath::Abs(v0->GetPindex());
  AliESDtrack *ptrack=esd->GetTrack(pidx);
  if (!AcceptTrack(ptrack)) return kFALSE;

  Float_t xy,z0;
  ntrack->GetImpactParameters(xy,z0);
  if (TMath::Abs(xy)<0.1) return kFALSE;
  ptrack->GetImpactParameters(xy,z0);
  if (TMath::Abs(xy)<0.1) return kFALSE;

  Double_t dca=v0->GetDcaV0Daughters();
  if (dca>1.0) return kFALSE;

  Double_t cpa=v0->GetV0CosineOfPointingAngle();
  if (cpa<0.998) return kFALSE;

  Double_t xx,yy,zz; v0->GetXYZ(xx,yy,zz);
  Double_t r2=xx*xx + yy*yy;
  if (r2<0.9*0.9) return kFALSE;
  if (r2>100*100) return kFALSE;

  return kTRUE;
}

static Bool_t AcceptCascade(const AliESDcascade *cs, const AliESDEvent *esd) {

  Int_t bidx=TMath::Abs(cs->GetBindex());
  AliESDtrack *btrack=esd->GetTrack(bidx);
  if (!AcceptTrack(btrack)) return kFALSE;

  Float_t xy,z0; 
  btrack->GetImpactParameters(xy,z0);
  if (TMath::Abs(xy)<0.03) return kFALSE;

  const AliESDVertex *vtx=esd->GetPrimaryVertexSPD();
  if (!vtx->GetStatus()) {
     vtx=esd->GetPrimaryVertexTracks();
     if (!vtx->GetStatus()) return kFALSE;
  }
  Double_t xv=vtx->GetXv(), yv=vtx->GetYv(), zv=vtx->GetZv();
  if (cs->GetCascadeCosineOfPointingAngle(xv,yv,zv) < 0.999) return kFALSE;

  if (cs->GetDcaXiDaughters() > 0.3) return kFALSE;

  return kTRUE;
}

static Bool_t AcceptPID(const AliPIDResponse *pidResponse, 
			const AliESDtrack *ptrack, AliStack *stack) {

  const AliExternalTrackParam *par=ptrack->GetInnerParam();
  if (!par) return kTRUE;
  if (par->GetP() > 1.) return kTRUE; 


  if (stack) {
    // MC PID
    Int_t plab=TMath::Abs(ptrack->GetLabel());
    TParticle *pp=stack->Particle(plab);
    if (!pp) return kTRUE;
    if (pp->GetPDG()->Charge() > 0) {
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

static TParticle*
AssociateV0(const AliESDtrack *ptrack,const AliESDtrack *ntrack,AliStack *stack,
TParticle *&mcp) {
  //
  // Try to associate the V0 with the daughters ptrack and ntrack
  // with the Monte Carlo
  //
  if (!stack) return 0;

  Int_t nlab=TMath::Abs(ntrack->GetLabel());
  TParticle *n=stack->Particle(nlab);
  if (!n) return 0;

  Int_t plab=TMath::Abs(ptrack->GetLabel());
  TParticle *p=stack->Particle(plab);
  if (!p) return 0;

  Int_t imp=p->GetFirstMother();
  if (imp<0) return 0;
  if (imp>=stack->GetNtrack()) return 0;
  TParticle *p0=stack->Particle(imp); // V0 particle, mother of the pos. track
  if (!p0) return 0;

  Int_t imn=n->GetFirstMother();
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


void AliAnalysisTaskCTauPbPb::UserExec(Option_t *)
{
  const Double_t yMax=0.5;
  const Double_t pMin=0.0;
  const Double_t lMax=0.001;

  AliESDEvent *esd=(AliESDEvent *)InputEvent();

  if (!esd) {
    Printf("ERROR: esd not available");
    return;
  }

  // Vertex selection
  const AliESDVertex *vtx=esd->GetPrimaryVertex();
  if (!vtx->GetStatus()) return;
  Double_t xv=vtx->GetXv(), yv=vtx->GetYv(), zv=vtx->GetZv();

  if (TMath::Abs(zv) > 10.) return ;   
 

  // Physics selection
  AliAnalysisManager *mgr= AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *hdr=(AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t maskIsSelected = hdr->IsEventSelected();
  Bool_t isSelected = (maskIsSelected & AliVEvent::kMB);
  if (!isSelected) return;


  fMult->Fill(-100); //event counter  


  // Centrality selection
  AliCentrality *cent=esd->GetCentrality();
  if (!cent->IsEventInCentralityClass(fCMin,fCMax,"V0M")) return;

  AliPIDResponse *pidResponse = hdr->GetPIDResponse(); 
 

  //+++++++ MC
  AliStack *stack = 0x0;
  Double_t mcXv=0., mcYv=0., mcZv=0.;

  if (fIsMC) {
     AliMCEvent *mcEvent = MCEvent();
     stack = mcEvent->Stack();
     if (!stack) {
        Printf("ERROR: stack not available");
        return;
     }

     const AliVVertex *mcVtx=mcEvent->GetPrimaryVertex();

     mcXv=mcVtx->GetX(); mcYv=mcVtx->GetY(); mcZv=mcVtx->GetZ();

     Int_t ntrk=stack->GetNtrack();
     while (ntrk--) {
       TParticle *p0=stack->Particle(ntrk);
       Int_t code=p0->GetPdgCode();
       if (code != kK0Short)
	 if (code != kLambda0)
	    if (code != kLambda0Bar) continue;

       Int_t plab=p0->GetFirstDaughter(), nlab=p0->GetLastDaughter();
       if (nlab==plab) continue;
       TParticle *part = stack->Particle(plab);
       if (!part) continue;
       TParticlePDG *partPDG = part->GetPDG();
       if (!partPDG) continue;
       Double_t charge=partPDG->Charge();
       if (charge == 0.) continue;
  
       Double_t pt=p0->Pt();
       if (pt<pMin) continue;
       if (TMath::Abs(p0->Y())>yMax) continue;
    
       Double_t x=p0->Vx(), y=p0->Vy(), z=p0->Vz();
       Double_t dx=mcXv-x, dy=mcYv-y, dz=mcZv-z;
       Double_t l=TMath::Sqrt(dx*dx + dy*dy + dz*dz);

       if (l > lMax) continue; // secondary V0

       x=part->Vx(); y=part->Vy();
       dx=mcXv-x; dy=mcYv-y;
       Double_t lt=TMath::Sqrt(dx*dx + dy*dy);

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


  Int_t ntrk1=esd->GetNumberOfTracks();
  Int_t mult=0;
  for (Int_t i=0; i<ntrk1; i++) {
    AliESDtrack *t=esd->GetTrack(i);
    if (!t->IsOn(AliESDtrack::kTPCrefit)) continue;
    Float_t xy,z0;
    t->GetImpactParameters(xy,z0);
    if (TMath::Abs(xy)>3.) continue;
    if (TMath::Abs(z0)>3.) continue;
    Double_t pt=t->Pt(),pz=t->Pz();
    if (TMath::Abs(pz/pt)>0.8) continue;
    mult++;

    Double_t p=t->GetInnerParam()->GetP();
    Double_t dedx=t->GetTPCsignal()/47.;
    fdEdx->Fill(p,dedx,1);

    Double_t nsig=pidResponse->NumberOfSigmasTPC(t,AliPID::kProton);
    if (TMath::Abs(nsig) < 3.) fdEdxPid->Fill(p,dedx,1);

  }
  fMult->Fill(mult);


  Int_t nv0 = esd->GetNumberOfV0s();
  while (nv0--) {
      AliESDv0 *v0=esd->GetV0(nv0);

      Double_t pt=v0->Pt();
      if (pt < pMin) continue;

      if (!AcceptV0(v0,esd)) continue;

      Int_t nidx=TMath::Abs(v0->GetNindex());
      AliESDtrack *ntrack=esd->GetTrack(nidx);
      Int_t pidx=TMath::Abs(v0->GetPindex());
      AliESDtrack *ptrack=esd->GetTrack(pidx);

      Double_t x,y,z; v0->GetXYZ(x,y,z);
      Double_t dx1=x-xv, dy1=y-yv;
      Double_t lt=TMath::Sqrt(dx1*dx1 + dy1*dy1);

      if (lt/pt > 3*7.89/1.1157) continue;  

      //--- V0 switches
      Bool_t isK0s=kTRUE;
      Bool_t isLambda=kTRUE;
      Bool_t isLambdaBar=kTRUE;

      if (0.4977*lt/pt > 3*2.68) isK0s=kFALSE;
      if (1.1157*lt/pt > 3*7.89) isLambdaBar=isLambda=kFALSE;

      if (v0->PtArmV0() < 0.2*TMath::Abs(v0->AlphaV0())) isK0s=kFALSE;

      if (!AcceptPID(pidResponse, ptrack, stack)) isLambda=kFALSE;
      if (!AcceptPID(pidResponse, ntrack, stack)) isLambdaBar=kFALSE;

      Double_t yK0s=TMath::Abs(v0->RapK0Short());
      Double_t yLam=TMath::Abs(v0->RapLambda());
      if (yK0s > yMax) isK0s=kFALSE;
      if (yLam > yMax) isLambda=isLambdaBar=kFALSE;
      //---

      Double_t mass=0., m=0., s=0.;
      if (isK0s) {
         v0->ChangeMassHypothesis(kK0Short);

         mass=v0->GetEffMass();
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
         v0->ChangeMassHypothesis(kLambda0);

         mass=v0->GetEffMass();
         fLambdaM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
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
         v0->ChangeMassHypothesis(kLambda0Bar);

         mass=v0->GetEffMass();
         fLambdaBarM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kLambda0Bar)->Mass();
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

      TParticle *mcp=0;
      TParticle *mc0=AssociateV0(ptrack,ntrack,stack,mcp);
      if (!mc0) continue;

      Double_t ptAs=mc0->Pt();
      if (ptAs < pMin) continue;
      Double_t yAs=mc0->Y();
      if (TMath::Abs(yAs) > yMax) continue;

      Int_t code=mc0->GetPdgCode();

      Double_t dx = mcXv - mcp->Vx(), dy = mcYv - mcp->Vy();
      Double_t ltAs=TMath::Sqrt(dx*dx + dy*dy);
 
      Double_t dz=mcZv - mc0->Vz(); dy=mcYv - mc0->Vy(); dx=mcXv - mc0->Vx();
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

         Int_t nx=mc0->GetFirstMother();
         if (nx<0) continue;
         if (nx>=stack->GetNtrack()) continue;
         TParticle *xi=stack->Particle(nx);
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

  Double_t kine0;
  Int_t ncs=esd->GetNumberOfCascades();
  for (Int_t i=0; i<ncs; i++) {
      AliESDcascade *cs=esd->GetCascade(i);

      Double_t pt=cs->Pt();
      if (pt < pMin) continue;
      if (TMath::Abs(cs->RapXi()) > yMax) continue;
      if (!AcceptCascade(cs,esd)) continue;

      AliESDv0 *v0 = (AliESDv0*)cs;
      if (v0->Pt() < pMin) continue;
      if (TMath::Abs(v0->RapLambda()) > yMax) continue;
      if (!AcceptV0(v0,esd)) continue;

      //--- Cascade switches
      Bool_t isXiMinus=kTRUE;
      Bool_t isXiPlusBar=kTRUE;

      Int_t pidx=TMath::Abs(v0->GetPindex());
      AliESDtrack *ptrack=esd->GetTrack(pidx);
      if (!AcceptPID(pidResponse, ptrack, stack)) isXiMinus=kFALSE;

      Int_t nidx=TMath::Abs(v0->GetNindex());
      AliESDtrack *ntrack=esd->GetTrack(nidx);
      if (!AcceptPID(pidResponse, ntrack, stack)) isXiPlusBar=kFALSE;

      Int_t charge=cs->Charge();
      if (charge > 0) isXiMinus=kFALSE;
      if (charge < 0) isXiPlusBar=kFALSE;
      //---
      
      if (isXiMinus) {
         cs->ChangeMassHypothesis(kine0,kXiMinus);
         Double_t mass=cs->GetEffMassXi();
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
         cs->ChangeMassHypothesis(kine0,kXiPlusBar);
         Double_t mass=cs->GetEffMassXi();
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

void AliAnalysisTaskCTauPbPb::Terminate(Option_t *)
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
