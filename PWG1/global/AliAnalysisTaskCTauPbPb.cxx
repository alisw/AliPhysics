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

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskCTauPbPb.h"

extern TROOT *gROOT;

ClassImp(AliAnalysisTaskCTauPbPb)

static Int_t    nbins=102;  // number of bins
static Double_t lMin=0.0, lMax=100.;
static Double_t pMin=0.5, pMax=10.;
static Double_t yMax=0.75;


//
//  This is a little task for checking the c*tau of the strange particles 
//

AliAnalysisTaskCTauPbPb::AliAnalysisTaskCTauPbPb(const char *name) :
AliAnalysisTaskSE(name),
fCMin(0.),
fCMax(90.),
fOutput(0),
fMult(0),

fK0sM(0),
fK0sSi(0),
fK0sMC(0),
fK0sAs(0),

fLambdaM(0),
fLambdaSi(0),
fLambdaMC(0),
fLambdaAs(0),

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
  new TH2F("fLambdaM","Mass for \\Lambda", nbins, 1.065, 1.165,10,pMin,pMax);
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
 
  //fMult->Fill(-100); //event counter  

  //+++++++ MC
  AliStack *stack = 0x0;
  const AliVVertex *mcVtx=0x0;
  Double_t mcXv=0., mcYv=0., mcZv=0.;

  AliMCEvent *mcEvent = MCEvent();
  if (mcEvent) {
     stack = mcEvent->Stack();
     if (!stack) {
        Printf("ERROR: stack not available");
        return;
     }

     mcVtx=mcEvent->GetPrimaryVertex();
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
       if (TMath::Abs(p0->Y())>0.75) continue;
    
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
  }
  fMult->Fill(mult);


  Int_t nv0 = esd->GetNumberOfV0s();
  while (nv0--) {
      AliESDv0 *v0=esd->GetV0(nv0);

      if (v0->GetOnFlyStatus()) continue;

      Int_t nidx=TMath::Abs(v0->GetNindex());
      AliESDtrack *ntrack=esd->GetTrack(nidx);
      if (!ntrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (ntrack->GetTPCclusters(0)<70) continue;
      if (ntrack->GetKinkIndex(0)>0) continue;

      Int_t pidx=TMath::Abs(v0->GetPindex());
      AliESDtrack *ptrack=esd->GetTrack(pidx);
      if (!ptrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (ptrack->GetTPCclusters(0)<70) continue;      
      if (ptrack->GetKinkIndex(0)>0) continue;


      Double_t pt=v0->Pt();
      if (pt<pMin) continue;

      // topological selections
      Float_t xy,z0;
      ntrack->GetImpactParameters(xy,z0);
      if (TMath::Abs(xy)<0.1) continue;
      ptrack->GetImpactParameters(xy,z0);
      if (TMath::Abs(xy)<0.1) continue;

      Double_t dca=v0->GetDcaV0Daughters();
      if (dca>1.0) continue;

      Double_t cpa=v0->GetV0CosineOfPointingAngle();
      if (cpa<0.998) continue;

      Double_t xx,yy,zz; v0->GetXYZ(xx,yy,zz);
      Double_t r2=xx*xx + yy*yy;
      if (r2<0.9*0.9) continue;
      if (r2>100*100) continue;
      //

      Double_t x,y,z; v0->GetXYZ(x,y,z);
      Double_t dx=x-xv, dy=y-yv;
      Double_t lt=TMath::Sqrt(dx*dx + dy*dy);


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
         if (np->GetFirstMother() != i0) goto noas;

         if (i0<0) goto noas;
         if (i0>=ntrk) goto noas;
         TParticle *p0=stack->Particle(i0);

         Int_t code=p0->GetPdgCode();
         if (code != kK0Short)
	    if (code != kLambda0) goto noas;

	 if (p0->Pt()<pMin) goto noas;
	 if (TMath::Abs(p0->Y())>0.75 ) goto noas;


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
	   if (code == kLambda0) fLambdaAs->Fill(ptAs,ltAs);
           else  fK0sAs->Fill(ptAs,ltAs);
         }
 
      }
      //++++++++

  noas:

      Double_t mass=0., m=0., s=0.;
      if (TMath::Abs(v0->RapK0Short())<0.75) {
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
      
      if (TMath::Abs(v0->RapLambda())<0.75) {
      v0->ChangeMassHypothesis(kLambda0);

      mass=v0->GetEffMass();
      fLambdaM->Fill(mass,pt);

      m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
      //s=0.0027;
      s=0.0027 + (0.004-0.0027)/(10-1)*(pt-1);
      if (TMath::Abs(m-mass) < 3*s) {
         fLambdaSi->Fill(pt,lt);
      }
      if (TMath::Abs(m-mass + 4.5*s) < 1.5*s) {
         fLambdaSi->Fill(pt,lt,-1);
      }
      if (TMath::Abs(m-mass - 4.5*s) < 1.5*s) {
         fLambdaSi->Fill(pt,lt,-1);
      }
      }
  }

  Double_t kine0;
  Int_t ncs=esd->GetNumberOfCascades();
  for (Int_t i=0; i<ncs; i++) {
      AliESDcascade *cs=esd->GetCascade(i);

      Int_t nidx=TMath::Abs(cs->GetNindex());
      AliESDtrack *ntrack=esd->GetTrack(nidx);
      if (!ntrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (ntrack->GetTPCclusters(0)<70) continue;
      if (ntrack->GetKinkIndex(0)>0) continue;

      Int_t pidx=TMath::Abs(cs->GetPindex());
      AliESDtrack *ptrack=esd->GetTrack(pidx);
      if (!ptrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (ptrack->GetTPCclusters(0)<70) continue;
      if (ptrack->GetKinkIndex(0)>0) continue;

      Int_t bidx=TMath::Abs(cs->GetBindex());
      AliESDtrack *btrack=esd->GetTrack(bidx);
      if (!btrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (btrack->GetTPCclusters(0)<70) continue;
      if (btrack->GetKinkIndex(0)>0) continue;

      Double_t pt=cs->Pt();

      // Xi topological cuts
      if (pt < pMin) continue; 
      if (TMath::Abs(cs->RapXi()) > yMax) continue;

      Float_t xy,z0;
      btrack->GetImpactParameters(xy,z0);
      if (TMath::Abs(xy)<0.03) continue;
      if (cs->GetCascadeCosineOfPointingAngle(xv,yv,zv) < 0.999) continue;
      if (cs->GetDcaXiDaughters() > 0.3) continue;


      // Lambda topological selections
      AliESDv0 *v0 = (AliESDv0*)cs;
      pt=v0->Pt();
      if (pt<pMin) continue;
      if (TMath::Abs(v0->RapLambda()) > yMax) continue;

      ntrack->GetImpactParameters(xy,z0);
      if (TMath::Abs(xy)<0.1) continue;
      ptrack->GetImpactParameters(xy,z0);
      if (TMath::Abs(xy)<0.1) continue;

      Double_t dca=v0->GetDcaV0Daughters();
      if (dca>1.0) continue;

      Double_t cpa=v0->GetV0CosineOfPointingAngle();
      if (cpa<0.998) continue;

      Double_t xx,yy,zz; v0->GetXYZ(xx,yy,zz);
      Double_t r2=xx*xx + yy*yy;
      if (r2<0.9*0.9) continue;
      if (r2>100*100) continue;


      Int_t charge=cs->Charge();      
      if (charge < 0) {
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
 
  /*
  fK0s = dynamic_cast<TH2F*>(fOutput->FindObject("fK0s")) ; 
  if (!fK0s) {
     Printf("ERROR: fK0s not available");
     return;
  }
  */
  fK0sSi = dynamic_cast<TH2F*>(fOutput->FindObject("fK0sSi")) ; 
  if (!fK0sSi) {
     Printf("ERROR: fK0sSi not available");
     return;
  }
  //fK0s->Sumw2(); fK0sSi->Sumw2();
  //fK0s->Add(fK0sSi,-1);  


  /*
  fLambda = dynamic_cast<TH2F*>(fOutput->FindObject("fLambda")) ; 
  if (!fLambda) {
     Printf("ERROR: fLambda not available");
     return;
  }
  */
  fLambdaSi = dynamic_cast<TH2F*>(fOutput->FindObject("fLambdaSi")) ; 
  if (!fLambdaSi) {
     Printf("ERROR: fLambdaSi not available");
     return;
  }
  //fLambda->Sumw2(); fLambdaSi->Sumw2();
  //fLambda->Add(fLambdaSi,-1);  

  if (!gROOT->IsBatch()) {
    /*
    TCanvas *c1 = new TCanvas("c1","Lambda");
    c1->SetLogy();
    c1->SetFillColor(10);
    c1->SetHighLightColor(10);
    fLambda->DrawCopy("E") ;
    */
    TCanvas *c2 = new TCanvas("c2","Lambda side-band subtracted");
    c2->SetFillColor(10);
    c2->SetHighLightColor(10);
    fLambdaSi->DrawCopy("E") ;

    /*
    TCanvas *c3 = new TCanvas("c3","K0s");
    c3->SetLogy();
    c3->SetFillColor(10);
    c3->SetHighLightColor(10);
    fK0s->DrawCopy("E") ;
    */
    TCanvas *c4 = new TCanvas("c4","K0s side-band subtracted");
    c4->SetFillColor(10);
    c4->SetHighLightColor(10);
    fK0sSi->DrawCopy("E") ;
  }
}
