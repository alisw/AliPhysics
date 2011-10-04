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

extern TROOT *gROOT;

ClassImp(AliAnalysisTaskCTauPbPb)

static Int_t    nbins=102;  // number of bins
static Double_t lMin=0.0, lMax=100.;
static Double_t pMin=0.0, pMax=10.;
static Double_t yMax=0.5;


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

fCPA(0),
fDCA(0),

fLambdaEff(0),
fLambdaPt(0),

fLambdaFromXi(0),
fXiM(0),
fXiSiP(0)
{
  // Constructor. Initialization of pointers
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskCTauPbPb::UserCreateOutputObjects()
{
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
  new TH2F("fK0sM", "Mass for K^{0}_{s}", nbins/2, 0.448, 0.548, 10,pMin,pMax);
  fK0sM->GetXaxis()->SetTitle("Mass [GeV/c]"); 
  fOutput->Add(fK0sM);

  fK0sSi = 
  new TH2F("fK0sSi","L_{T} vs p_{T} for K^{0}_{s}, side-band subtracted",
  nbins,pMin,pMax,nbins,lMin,lMax);
  fK0sSi->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  fK0sSi->GetYaxis()->SetTitle("L_{T} [cm]"); 
  fOutput->Add(fK0sSi);

  fK0sMC = 
  new TH2F("fK0sMC","L_{T} vs p_{T} for K^{0}_{s}, from MC stack", 
  nbins,pMin,pMax,nbins,lMin,lMax);
  fK0sMC->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  fK0sMC->GetYaxis()->SetTitle("L_{T} [cm]"); 
  fOutput->Add(fK0sMC);

  fK0sAs = 
  new TH2F("fK0sAs", "L_{T} vs p_{T} for K^{0}_{s}, associated", 
  nbins,pMin,pMax,nbins,lMin,lMax);
  fK0sAs->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  fK0sAs->GetYaxis()->SetTitle("L_{T} [cm]"); 
  fOutput->Add(fK0sAs);

  //----------------------

  fLambdaM = 
  new TH2F("fLambdaM","Mass for \\Lambda", nbins, 1.065, 1.165,nbins,pMin,pMax);
  //new TH2F("fLambdaM","Mass for \\Lambda", nbins, 1.065, 1.165,10,0.1,1.1);
  fLambdaM->GetXaxis()->SetTitle("Mass [GeV/c]"); 
  fOutput->Add(fLambdaM);

  fLambdaSi = 
  new TH2F("fLambdaSi","L_{T} vs p_{T} for \\Lambda, side-band subtructed",
  nbins,pMin,pMax,nbins,lMin,lMax);
  fLambdaSi->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  fLambdaSi->GetYaxis()->SetTitle("L_{T} [cm]"); 
  fOutput->Add(fLambdaSi);

  fLambdaMC = 
  new TH2F("fLambdaMC","c\\tau for \\Lambda, from MC stack", 
  nbins,pMin,pMax,nbins,lMin,lMax);
  fLambdaMC->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  fLambdaMC->GetYaxis()->SetTitle("L_{T} [cm]"); 
  fOutput->Add(fLambdaMC);

  fLambdaAs = 
  new TH2F("fLambdaAs","c\\tau for \\Lambda, associated",
  nbins,pMin,pMax,nbins,lMin,lMax);
  fLambdaAs->GetXaxis()->SetTitle("p_{T} [GeV/c]"); 
  fLambdaAs->GetYaxis()->SetTitle("L_{T} [cm]"); 
  fOutput->Add(fLambdaAs);

  fCPA=new TH1F("fCPA","Cosine of the pointing angle",30,0.9978,1.);
  fOutput->Add(fCPA);
  fDCA=new TH1F("fDCA","DCA between the daughters",30,0.,1.1);
  fOutput->Add(fDCA);

  fLambdaEff=fLambdaAs->ProjectionX();
  fLambdaEff->SetName("fLambdaEff");
  fLambdaEff->SetTitle("Efficiency for #Lambda");
  fOutput->Add(fLambdaEff);

  fLambdaPt=fLambdaAs->ProjectionX();
  fLambdaPt->SetName("fLambdaPt");
  fLambdaPt->SetTitle("Raw #Lambda pT spectrum");
  fOutput->Add(fLambdaPt);

  //----------------------

  fLambdaFromXi=new TH3F("fLambdaFromXi","L_{T} vs p_{T} vs p_{T} of \\Xi for \\Lambda from Xi",
  nbins,pMin,pMax,nbins,lMin,lMax,33,pMin,pMax+2);
  fOutput->Add(fLambdaFromXi);

  fXiM  = 
  new TH2F("fXiM", "\\Xi mass distribution", 50, 1.271, 1.371,12,pMin,pMax+2);
  fOutput->Add(fXiM);

  fXiSiP  = new TH1F("fXiSiP", "Pt for \\Xi, side-band subracted",
  33,pMin,pMax+2);
  fOutput->Add(fXiSiP);


  PostData(1, fOutput);
}

static Bool_t AcceptTrack(const AliESDtrack *t) {
  if (!t->IsOn(AliESDtrack::kTPCrefit)) return kFALSE;
  if (t->GetKinkIndex(0)>0) return kFALSE;

  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=t->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;

  return kTRUE;   
}

static Bool_t AcceptV0(const AliESDv0 *v0, const AliESDEvent *esd) {

  if (v0->GetOnFlyStatus()) return kFALSE;

  if (v0->Pt() < pMin) return kFALSE;

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
  //if (dca>0.7) return kFALSE;
  //if (dca>0.4) return kFALSE;

  Double_t cpa=v0->GetV0CosineOfPointingAngle();
  if (cpa<0.998) return kFALSE;
  //if (cpa<0.99875) return kFALSE;
  //if (cpa<0.9995) return kFALSE;

  Double_t xx,yy,zz; v0->GetXYZ(xx,yy,zz);
  Double_t r2=xx*xx + yy*yy;
  if (r2<0.9*0.9) return kFALSE;
  if (r2>100*100) return kFALSE;

  return kTRUE;
}

static Bool_t AcceptCascade(const AliESDcascade *cs, const AliESDEvent *esd) {

  if (cs->Pt() < pMin) return kFALSE;

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

void AliAnalysisTaskCTauPbPb::UserExec(Option_t *)
{

  AliESDEvent *esd=(AliESDEvent *)InputEvent();

  if (!esd) {
    Printf("ERROR: esd not available");
    return;
  }

  fMult->Fill(-100); //event counter  

  // Physics selection
  AliAnalysisManager *mgr= AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *hdr=(AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t maskIsSelected = hdr->IsEventSelected();
  Bool_t isSelected = (maskIsSelected & AliVEvent::kMB);
  if (!isSelected) return;

  // Centrality selection
  AliCentrality *cent=esd->GetCentrality();
  if (!cent->IsEventInCentralityClass(fCMin,fCMax,"V0M")) return;

  const AliESDVertex *vtx=esd->GetPrimaryVertexSPD();
  if (!vtx->GetStatus()) {
     vtx=esd->GetPrimaryVertexTracks();
     if (!vtx->GetStatus()) return;
  }
  Double_t xv=vtx->GetXv(), yv=vtx->GetYv(), zv=vtx->GetZv();

  if (TMath::Abs(zv) > 10.) return ;   
 
  AliPIDResponse *pidResponse = hdr->GetPIDResponse(); 
 
  //fMult->Fill(-100); //event counter  

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

     Int_t ntrk=stack->GetNtrack(), ntrk0=ntrk;
     while (ntrk--) {
       TParticle *p0=stack->Particle(ntrk);
       Int_t code=p0->GetPdgCode();
       if (code != kK0Short)
	 if (code != kLambda0) continue;

       Int_t plab=p0->GetFirstDaughter(), nlab=p0->GetLastDaughter();
       if (nlab==plab) continue;
       if (nlab<0) continue;
       if (plab<0) continue;
       if (nlab>=ntrk0) continue;
       if (plab>=ntrk0) continue;
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

       if (l > 0.01) continue; // secondary V0

       x=part->Vx(); y=part->Vy();
       dx=mcXv-x; dy=mcYv-y;
       Double_t lt=TMath::Sqrt(dx*dx + dy*dy);

       if (code == kK0Short) {
          fK0sMC->Fill(pt,lt);
       }
       if (code == kLambda0) {
          fLambdaMC->Fill(pt,lt);
       }
     }
  }


  Int_t ntrk=esd->GetNumberOfTracks();
  Int_t mult=0;
  Double_t nsig;
  for (Int_t i=0; i<ntrk; i++) {
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

    nsig=pidResponse->NumberOfSigmasTPC(t,AliPID::kProton);
    if (TMath::Abs(nsig) < 3.) fdEdxPid->Fill(p,dedx,1);

  }
  fMult->Fill(mult);


  Int_t nv0 = esd->GetNumberOfV0s();
  while (nv0--) {
      AliESDv0 *v0=esd->GetV0(nv0);

      if (!AcceptV0(v0,esd)) continue;

      Int_t nidx=TMath::Abs(v0->GetNindex());
      AliESDtrack *ntrack=esd->GetTrack(nidx);
      Int_t pidx=TMath::Abs(v0->GetPindex());
      AliESDtrack *ptrack=esd->GetTrack(pidx);

      Double_t x,y,z; v0->GetXYZ(x,y,z);
      Double_t dx=x-xv, dy=y-yv;
      Double_t lt=TMath::Sqrt(dx*dx + dy*dy);

      Double_t pt=v0->Pt();

      Bool_t ctK=kTRUE; if (0.4977*lt/pt > 3*2.68) ctK=kFALSE;
      Bool_t ctL=kTRUE; if (1.1157*lt/pt > 3*7.89) ctL=kFALSE;

      //+++++++ MC
      if (stack) {
         Int_t ntrk=stack->GetNtrack();

         Int_t nlab=TMath::Abs(ntrack->GetLabel());
         Int_t plab=TMath::Abs(ptrack->GetLabel());

         if (nlab<0) goto noas;      
         if (nlab>=ntrk) goto noas;      
         if (plab<0) goto noas;      
         if (plab>=ntrk) goto noas;      

         TParticle *np=stack->Particle(nlab);
         TParticle *pp=stack->Particle(plab);
         Int_t i0=pp->GetFirstMother();
         //if (np->GetFirstMother() != i0) goto noas;

         Int_t in0=np->GetFirstMother();
         if (in0<0) goto noas;
         if (in0>=ntrk) goto noas;
         if (in0 != i0) { // did the negative daughter decay ?
            TParticle *nnp=stack->Particle(in0);
            if (nnp->GetFirstMother() != i0) goto noas;
	 }

         if (i0<0) goto noas;
         if (i0>=ntrk) goto noas;
         TParticle *p0=stack->Particle(i0);

         Int_t code=p0->GetPdgCode();
         if (code != kK0Short)
	    if (code != kLambda0) goto noas;

	 if (p0->Pt()<pMin) goto noas;
	 if (TMath::Abs(p0->Y())>yMax ) goto noas;


         Double_t dz=mcZv - p0->Vz(), dy=mcYv - p0->Vy(), dx=mcXv - p0->Vx();
         Double_t l = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

         dx = mcXv - pp->Vx(); dy = mcYv - pp->Vy();
         Double_t ltAs=TMath::Sqrt(dx*dx + dy*dy);
         Double_t ptAs=p0->Pt();

	 if (l > 0.01) { // Secondary V0
	   if (code != kLambda0) goto noas;
           Int_t nx=p0->GetFirstMother();
           if (nx<0) goto noas;
           if (nx>=ntrk) goto noas;
           TParticle *xi=stack->Particle(nx);
           Int_t xcode=xi->GetPdgCode();
           if ( xcode != kXiMinus )
	     if( xcode != 3322 ) goto noas; 
	   fLambdaFromXi->Fill(ptAs,ltAs,xi->Pt());
	 } else {
	   if (code == kLambda0) {
              if (ctL) fLambdaAs->Fill(ptAs,ltAs);
           } else {
              if (ctK)  fK0sAs->Fill(ptAs,ltAs);
	   } 
	 }
 
      }
      //++++++++

  noas:

      Double_t dca=v0->GetDcaV0Daughters();
      Double_t cpa=v0->GetV0CosineOfPointingAngle();

      Double_t mass=0., m=0., s=0.;
      if (ctK)
      if (TMath::Abs(v0->RapK0Short())<yMax) {
         v0->ChangeMassHypothesis(kK0Short);

         mass=v0->GetEffMass();
         fK0sM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
         s=0.0044 + (0.008-0.0044)/(10-1)*(pt - 1.);
         if (TMath::Abs(m-mass) < 3*s) {
            fK0sSi->Fill(pt,lt);
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fK0sSi->Fill(pt,lt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fK0sSi->Fill(pt,lt,-1);
         }
      }
      
      if (ctL)
      if (TMath::Abs(v0->RapLambda())<yMax) {
         Double_t p=ptrack->GetInnerParam()->GetP();
         if (p<1.) {
            nsig=pidResponse->NumberOfSigmasTPC(ptrack,AliPID::kProton);
            if (TMath::Abs(nsig) > 3.) continue;
	 }
         v0->ChangeMassHypothesis(kLambda0);

         mass=v0->GetEffMass();
         fLambdaM->Fill(mass,pt);

         m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
         //s=0.0027 + (0.004-0.0027)/(10-1)*(pt-1);
         //s=0.0015 + (0.002-0.0015)/(2.6-1)*(pt-1);
         s=0.0023 + (0.004-0.0023)/(6-1)*(pt-1);
         if (TMath::Abs(m-mass) < 3*s) {
            fLambdaSi->Fill(pt,lt);
            fCPA->Fill(cpa,1);
            fDCA->Fill(dca,1);
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fLambdaSi->Fill(pt,lt,-1);
            fCPA->Fill(cpa,-1);
            fDCA->Fill(dca,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fLambdaSi->Fill(pt,lt,-1);
            fCPA->Fill(cpa,-1);
            fDCA->Fill(dca,-1);
         }
      }
  }

  Double_t kine0;
  Int_t ncs=esd->GetNumberOfCascades();
  for (Int_t i=0; i<ncs; i++) {
      AliESDcascade *cs=esd->GetCascade(i);

      if (TMath::Abs(cs->RapXi()) > yMax) continue;
      if (!AcceptCascade(cs,esd)) continue;

      AliESDv0 *v0 = (AliESDv0*)cs;
      if (TMath::Abs(v0->RapLambda()) > yMax) continue;
      if (!AcceptV0(v0,esd)) continue;

      Double_t pt=cs->Pt();

      Int_t charge=cs->Charge();      
      if (charge < 0) {         
         Int_t pidx=TMath::Abs(v0->GetPindex());
         AliESDtrack *ptrack=esd->GetTrack(pidx);
         Double_t p=ptrack->GetInnerParam()->GetP();
         if (p<1.) {
            nsig=pidResponse->NumberOfSigmasTPC(ptrack,AliPID::kProton);
            if (TMath::Abs(nsig) > 3.) continue;
	 }
         cs->ChangeMassHypothesis(kine0,kXiMinus);
         Double_t mass=cs->GetEffMassXi();
	 pt=cs->Pt();       
         fXiM->Fill(mass,pt);
         Double_t m=TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass();
         //Double_t s=0.0037;
         Double_t s=0.002 + (0.0032-0.002)/(6-1.5)*(pt-1.5);
         if (TMath::Abs(m-mass) < 3*s) {
            fXiSiP->Fill(pt);
         }
         if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
            fXiSiP->Fill(pt,-1);
         }
         if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
            fXiSiP->Fill(pt,-1);
         }
      }
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
 
    } else {
       new TCanvas("c6","Raw #Lambda pt");
       lambdaSiPx->SetTitle("Raw #Lambda pt");
      *fLambdaPt = *lambdaSiPx; 
       fLambdaPt->SetLineColor(2);
       fLambdaPt->DrawCopy("E");
    }
  }
}
