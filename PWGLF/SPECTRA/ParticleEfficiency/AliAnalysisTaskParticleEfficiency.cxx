#include "AliAnalysisTaskParticleEfficiency.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "TParticle.h"
#include "TObjArray.h"
#include "TList.h"
#include "TH2F.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

ClassImp(AliAnalysisTaskParticleEfficiency)

//_______________________________________________________

AliAnalysisTaskParticleEfficiency::AliAnalysisTaskParticleEfficiency(const Char_t *partName) :
  AliAnalysisTaskSE(partName),
  fParticlePdgCode(0),
  fTrackCuts(NULL)
{

  /* 
   * default constructor
   */

  /* set particle PDG code */
  TParticlePDG *ppdg = TDatabasePDG::Instance()->GetParticle(partName);
  if (ppdg) fParticlePdgCode = ppdg->PdgCode();

  /* init track cuts */
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  fTrackCuts->SetEtaRange(-0.8, 0.8);
  fTrackCuts->SetPtRange(0.15, 1.e10);

  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);
  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskParticleEfficiency::~AliAnalysisTaskParticleEfficiency()
{

  /* 
   * default destructor
   */
  
}

//_______________________________________________________

void
AliAnalysisTaskParticleEfficiency::UserCreateOutputObjects()
{

  /* 
   * user create output objects
   */

  fHistoEvents = new TH2F("hHistoEvents", ";centrality percentile;track multiplicity", 20, 0., 100., 200, 0., 2000.);
  fHistoList->Add(fHistoEvents);

  fHistoGenerated = new TH2F("hHistoGenerated", ";centrality percentile;p_{T} (GeV/c)", 20, 0., 100., 200, 0., 10.);
  fHistoList->Add(fHistoGenerated);
  
  fHistoReconstructed = new TH2F("hHistoReconstructed", ";centrality percentile;p_{T} (GeV/c)", 20, 0., 100., 200, 0., 10.);
  fHistoList->Add(fHistoReconstructed);
  
  for (Int_t idau = 0; idau < 2; idau++) {
    
    fHistoGeneratedDaughter[idau] = new TH2F(Form("hHistoGeneratedDaughter_%d", idau), ";centrality percentile;p_{T} (GeV/c)", 20, 0., 100., 200, 0., 10.);
    fHistoList->Add(fHistoGeneratedDaughter[idau]);
    
    fHistoReconstructedDaughter[idau] = new TH2F(Form("hHistoReconstructedDaughter_%d", idau), ";centrality percentile;p_{T} (GeV/c)", 20, 0., 100., 200, 0., 10.);
    fHistoList->Add(fHistoReconstructedDaughter[idau]);
    
  }  

  PostData(1, fHistoList);
}

//_______________________________________________________

void
AliAnalysisTaskParticleEfficiency::UserExec(Option_t *)
{

  /* 
   * user exec
   */

  /*** INIT ***/

  /* get ESD event */
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) return;
  /* get MC event */
  AliMCEvent *mcEvent = dynamic_cast<AliMCEvent *>(MCEvent());
  if (!mcEvent) return;
  /* get stack */
  AliStack *mcStack = mcEvent->Stack();
  if (!mcStack) return;

  /*** EVENT SELECTION ***/

  /* collision candidate */
  if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;
  /* vertex selection */
  const AliESDVertex *vertex = esdEvent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    vertex = esdEvent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) return;
  }
  /* z-vertex cut */
  if (TMath::Abs(vertex->GetZ()) > 10.) return;
  /* centrality selection */
  AliCentrality *centrality = esdEvent->GetCentrality();
  if (centrality->GetQuality() != 0) return;
  Double_t cent = centrality->GetCentralityPercentileUnchecked("V0M");

  /* fill event histos */
  fHistoEvents->Fill(cent, fTrackCuts->CountAcceptedTracks(esdEvent));

  /*** RECONSTRUCTED TRACKS ***/

  TObjArray recoParticleArray;

  /* loop over ESD tracks */
  for (Int_t itrk = 0; itrk < esdEvent->GetNumberOfTracks(); itrk++) {

    /* get track */
    AliESDtrack *track = esdEvent->GetTrack(itrk);
    if (!track) continue;
    /* check accept track */
    if (!fTrackCuts->AcceptTrack(track)) continue;

    /* get coresponding MC particle */
    TParticle *particle = mcStack->Particle(TMath::Abs(track->GetLabel()));
    if (!particle) continue;

    /* add to reconstructed particle array */
    recoParticleArray.Add(particle);
  }

  /*** MONTECARLO PARTICLES ***/

  /* loop over MC stack */
  for (Int_t ipart = 0; ipart < mcStack->GetNtrack(); ipart++) {
    /* get particle */
    TParticle *particle = mcStack->Particle(ipart);
    if (!particle) continue;
    /* check PDG */
    if (particle->GetPdgCode() != fParticlePdgCode) continue;
    /* check rapidity */
    if (TMath::Abs(particle->Y()) > 0.5) continue;

    /* switch PDG */
    switch (fParticlePdgCode) {
      
    case 211: case -211: /* pions */
    case 321: case -321: /* kaons */
    case 2212: case -2212: /* protons */

      /* check physical primary */
      if (!mcStack->IsPhysicalPrimary(ipart)) continue;

      /* fill particle histos */
      fHistoGenerated->Fill(cent, particle->Pt());
      if (recoParticleArray.Contains(particle))
	fHistoReconstructed->Fill(cent, particle->Pt());
      break;

    case 333: /* phi */
    case 313: /* K*0 */
    
      /* check number of daughters */
      if (particle->GetNDaughters() != 2) continue;
      /* check daughter charge */
      TParticle *daughter0 =  mcStack->Particle(particle->GetDaughter(0));
      TParticle *daughter1 =  mcStack->Particle(particle->GetDaughter(1));
      if (TDatabasePDG::Instance()->GetParticle(daughter0->GetPdgCode())->Charge() == 0.) continue;
      if (TDatabasePDG::Instance()->GetParticle(daughter1->GetPdgCode())->Charge() == 0.) continue;

      /* fill daughter histos */
      fHistoGeneratedDaughter[0]->Fill(cent, daughter0->Pt());
      if (recoParticleArray.Contains(daughter0))
	fHistoReconstructedDaughter[0]->Fill(cent, daughter0->Pt());
      fHistoGeneratedDaughter[1]->Fill(cent, daughter1->Pt());
      if (recoParticleArray.Contains(daughter1))
	fHistoReconstructedDaughter[1]->Fill(cent, daughter1->Pt());

      /* fill particle histos */
      fHistoGenerated->Fill(cent, particle->Pt());
      if (recoParticleArray.Contains(daughter0) &&
	  recoParticleArray.Contains(daughter1)) 
	fHistoReconstructed->Fill(cent, particle->Pt());
      break;
    }
    
  }

  PostData(1, fHistoList);
}

