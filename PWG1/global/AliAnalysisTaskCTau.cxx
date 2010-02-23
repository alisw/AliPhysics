#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticle.h>
#include <TROOT.h>

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisTaskCTau.h"

extern TROOT *gROOT;

//
//  This is a little task for checking the c*tau of the strange particles 
//

AliAnalysisTaskCTau::AliAnalysisTaskCTau(const char *name) :
AliAnalysisTaskSE(name),
fOutput(0),
fESD(0), 
fK0s(0),
fK0sMC(0),
fLambdas(0),
fLambdasMC(0),
fLambdaBars(0),
fLambdaBarsMC(0),
fXis(0),
fXisMC(0),
fMass(0),
fMassMC(0)
{
  // Constructor. Initialization of pointers
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskCTau::UserCreateOutputObjects()
{
  fOutput = new TList(); 
  fOutput->SetOwner();

  fK0s     = new TH1F("K0s", "c\\tau for K^{0}_{s}", 100, 0., 25.);//10*2.684);
  fK0s->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
  fOutput->Add(fK0s);
  fK0sMC   = new TH1F("K0sMC", "c\\tau for K^{0}_{s}, MC", 100, 0., 25.);
  fK0sMC->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
  fOutput->Add(fK0sMC);

  fLambdas = new TH1F("Lambdas","c\\tau for \\Lambda", 100, 0.,80.);//10*7.89);
  fLambdas->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
   //fLambdas->Sumw2();
  fOutput->Add(fLambdas);
  fLambdasMC = new TH1F("LambdasMC","c\\tau for \\Lambda, MC", 100, 0.,80.);
  fLambdasMC->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
   //fLambdasMC->Sumw2();
  fOutput->Add(fLambdasMC);

  fLambdaBars = new TH1F("LambdaBars","c\\tau for \\bar{\\Lambda}",100,0.,80.);
  fLambdaBars->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
   //fLambdaBars->Sumw2();
  fOutput->Add(fLambdaBars);
  fLambdaBarsMC=new TH1F("LambdaBarsMC","c\\tau for \\bar{\\Lambda}, MC",100,0.,80.);
  fLambdaBarsMC->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
   //fLambdaBarsMC->Sumw2();
  fOutput->Add(fLambdaBarsMC);

  fXis = new TH1F("Xis","c\\tau for \\Xi^{-} + \\Xi^{+}",50,0.,30.);//10*4.91)
  fXis->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
  fOutput->Add(fXis);
  fXisMC = new TH1F("XisMC","c\\tau for \\Xi^{-} + \\Xi^{+}, MC",50,0.,30.);
  fXisMC->GetXaxis()->SetTitle("L\\frac{m}{p} [cm]"); 
  fOutput->Add(fXisMC);

   //fMass = new TH1F("Mass", "Mass", 50, 0.4477, 0.5477);
  fMass = new TH1F("Mass", "Mass", 50, 1.07, 1.16);
   //fMass = new TH1F("Mass", "Mass", 50, 1.271, 1.371);
  fOutput->Add(fMass);
  fMassMC = new TH1F("MassMC", "Mass, MC", 50, 1.07, 1.16);
  fOutput->Add(fMassMC);

}

void AliAnalysisTaskCTau::UserExec(Option_t *)
{

  fESD=(AliESDEvent *)InputEvent();

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  //Char_t *tname="CINT1B-ABCE-NOPF-ALL";
  //Int_t ok=fESD->IsTriggerClassFired(tname);
  //if (!ok) return;
  
  const AliESDVertex *vtx=fESD->GetPrimaryVertexSPD();
  //if (!vtx->GetStatus()) {
  //   vtx=fESD->GetPrimaryVertexTracks();
  //   if (!vtx->GetStatus()) return;
  //}
  Double_t xv=vtx->GetXv(), yv=vtx->GetYv(), zv=vtx->GetZv();


  //+++++++ MC
  AliStack *stack = 0x0;
  AliMCEvent *mcEvent = MCEvent();
  if (mcEvent) {
     stack = mcEvent->Stack();
     if (!stack) {
        Printf("ERROR: stack not available");
        return;
     }

     const AliVVertex *mcVtx=mcEvent->GetPrimaryVertex();
     Double_t mcXv=mcVtx->GetX(), mcYv=mcVtx->GetY(), mcZv=mcVtx->GetZ();
     Int_t ntrk=stack->GetNtrack();
     while (ntrk--) {
       TParticle *p0=stack->Particle(ntrk);
       Int_t code=p0->GetPdgCode();
       if (code != kK0Short)
       if (code != kLambda0)
       if (code != kLambda0Bar) continue;

       Int_t plab=p0->GetFirstDaughter(), nlab=p0->GetLastDaughter();
       if (nlab==plab) continue;
       if (nlab<0) continue;
       if (plab<0) continue;
       TParticle *part = stack->Particle(plab);
       if (!part) continue;
       TParticlePDG *partPDG = part->GetPDG();
       if (!partPDG) continue;
       Double_t charge=partPDG->Charge();
       if (charge == 0.) continue;
       if (charge < 0.) {
          Int_t i=plab; plab=nlab; nlab=i;
       }
  
       Double_t pz=p0->Pz(), pt=p0->Pt(), p=p0->P();
       if (pt<0.1) continue;
       if (TMath::Abs(pz/pt)>0.7) continue;
       if (p<0.9) continue;
    
       Double_t x=p0->Vx(), y=p0->Vy(), z=p0->Vz();
       Double_t dx=mcXv-x, dy=mcYv-y, dz=mcZv-z;
       Double_t l=TMath::Sqrt(dx*dx + dy*dy + dz*dz);

       if (l > 0.1) continue; // secondary V0

       x=part->Vx(); y=part->Vy(); z=part->Vz();
       dx=mcXv-x; dy=mcYv-y; dz=mcZv-z;
       l=TMath::Sqrt(dx*dx + dy*dy + dz*dz);

       if (code == kK0Short) {
          Double_t m=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
          Double_t ct=l*m/p;    
          fK0sMC->Fill(ct);
       }
       if (code == kLambda0) {
          Double_t m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
          Double_t ct=l*m/p;    
          fLambdasMC->Fill(ct);
       }
       if (code == kLambda0Bar) {
         Double_t m=TDatabasePDG::Instance()->GetParticle(kLambda0Bar)->Mass();
         Double_t ct=l*m/p;    
         fLambdaBarsMC->Fill(ct);
       }
     }
  }

  Int_t nv0 = fESD->GetNumberOfV0s();
  while (nv0--) {
      AliESDv0 *v0=fESD->GetV0(nv0);

      if (v0->GetOnFlyStatus()) continue;

      Int_t nidx=TMath::Abs(v0->GetNindex());
      AliESDtrack *ntrack=fESD->GetTrack(nidx);
      if (!ntrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      //if (!ntrack->HasPointOnITSLayer(5)) continue;

      Int_t pidx=TMath::Abs(v0->GetPindex());
      AliESDtrack *ptrack=fESD->GetTrack(pidx);
      if (!ptrack->IsOn(AliESDtrack::kTPCrefit)) continue;
      //if (!ptrack->HasPointOnITSLayer(5)) continue;


      Double_t pz=v0->Pz(),pt=v0->Pt(),p=v0->P();
      if (pt<0.1) continue;
      if (TMath::Abs(pz/pt)>0.7) continue;
      if (p<0.9) continue;

      Double_t cpa=v0->GetV0CosineOfPointingAngle();
      if (cpa<0.995) continue;
      

      Double_t x,y,z; v0->GetXYZ(x,y,z);
      Double_t dx=x-xv, dy=y-yv, dz=z-zv;
      //Double_t r=TMath::Sqrt(x*x + y*y);
      //if (r > 3.) continue;
      Double_t l=TMath::Sqrt(dx*dx + dy*dy + dz*dz);

      v0->ChangeMassHypothesis(kK0Short);
      Double_t mass=v0->GetEffMass();
      Double_t m=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
      if (TMath::Abs(m-mass)<0.016) {
	 Double_t ct=l*m/p;
         fK0s->Fill(ct);
      }
      
      //if (ntrack->GetP()<0.7) continue;

      v0->ChangeMassHypothesis(kLambda0);
      mass=v0->GetEffMass();
      m=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
      if (TMath::Abs(m-mass)<0.006) {
	 Double_t ct=l*m/p;
         fLambdas->Fill(ct);
      }

      v0->ChangeMassHypothesis(kLambda0Bar);
      mass=v0->GetEffMass();
      m=TDatabasePDG::Instance()->GetParticle(kLambda0Bar)->Mass();
      if (TMath::Abs(m-mass)<0.006) {
	 Double_t ct=l*m/p;
         fLambdaBars->Fill(ct);
      }
      fMass->Fill(mass);


      //+++++++ MC
      Int_t nlab=TMath::Abs(ntrack->GetLabel());
      Int_t plab=TMath::Abs(ptrack->GetLabel());
      
      if ((plab-nlab) != 1)
	if ((nlab-plab) != 1) continue;

      if (!stack) continue;
      TParticle *np=stack->Particle(nlab);
      Int_t i0=np->GetFirstMother();
      TParticle *p0=stack->Particle(i0);
      if (!p0) continue;
      if (p0->GetPdgCode() == kLambda0Bar) fMassMC->Fill(mass);

  }

  Double_t kine0;
  Int_t ncs=fESD->GetNumberOfCascades();
  for (Int_t i=0; i<ncs; i++) {
      AliESDcascade *cs=fESD->GetCascade(i);

      Int_t nidx=TMath::Abs(cs->GetNindex());
      AliESDtrack *ntrack=fESD->GetTrack(nidx);
      if (!ntrack->IsOn(AliESDtrack::kTPCrefit)) continue;

      Int_t pidx=TMath::Abs(cs->GetPindex());
      AliESDtrack *ptrack=fESD->GetTrack(pidx);
      if (!ptrack->IsOn(AliESDtrack::kTPCrefit)) continue;

      Int_t bidx=TMath::Abs(cs->GetBindex());
      AliESDtrack *btrack=fESD->GetTrack(bidx);
      if (!btrack->IsOn(AliESDtrack::kTPCrefit)) continue;

      Int_t charge=cs->Charge();      
      if (charge<0)      
         cs->ChangeMassHypothesis(kine0,kXiMinus);
      else
         cs->ChangeMassHypothesis(kine0,kXiPlusBar);

      Double_t mass=cs->GetEffMassXi();

      if (mass>1.305)
      if (mass<1.335) {
        //TTree *tree=fChain->GetTree();
        //TDirectory *dir=tree->GetDirectory();
        //Int_t n=fESD->GetEventNumberInFile();
        //Printf("%d %d Xi found !!! %s %d\n", ncs, charge, dir->GetName(), n);
        Double_t p=cs->P();
        Double_t m=TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass();
        Double_t x,y,z; cs->GetXYZcascade(x,y,z);
        Double_t dx=x-xv, dy=y-yv, dz=z-zv;
        Double_t l=TMath::Sqrt(dx*dx + dy*dy + dz*dz);

        if (p<0.1) continue;

        Double_t ct=l*m/p;
        fXis->Fill(ct);
      }

  }
  
  PostData(1, fOutput);
}

void AliAnalysisTaskCTau::Terminate(Option_t *)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  fOutput=(TList*)GetOutputData(1);
  if (!fOutput) {
     Printf("ERROR: fOutput not available");
     return;
  }
 
  fMass = dynamic_cast<TH1F*>(fOutput->FindObject("Mass")) ; 
  if (!fMass) {
     Printf("ERROR: fMass not available");
     return;
  }
  fLambdas = dynamic_cast<TH1F*>(fOutput->FindObject("Lambdas")) ; 
  if (!fLambdas) {
     Printf("ERROR: fLambdas not available");
     return;
  }
  fLambdaBars = dynamic_cast<TH1F*>(fOutput->FindObject("LambdaBars")) ; 
  if (!fLambdaBars) {
     Printf("ERROR: fLambdaBars not available");
     return;
  }
  fXis = dynamic_cast<TH1F*>(fOutput->FindObject("Xis")) ; 
  if (!fXis) {
     Printf("ERROR: fXis not available");
     return;
  }
  
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Lambdas");
    c1->SetFillColor(10);
    c1->SetHighLightColor(10);
    fLambdas->DrawCopy("E") ;
    TCanvas *c2 = new TCanvas("c2","LambdaBar");
    c2->SetFillColor(10);
    c2->SetHighLightColor(10);
    fLambdaBars->DrawCopy("E") ;
  }
}
