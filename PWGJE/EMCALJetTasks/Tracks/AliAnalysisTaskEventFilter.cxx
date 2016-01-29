/*
 * AliAnalysisTaskEventFilter.cxx
 *
 *  Created on: Jan 29, 2016
 *      Author: markus
 */
#include <THashList.h>
#include <THistManager.h>
#include <TString.h>

#include <AliAnalysisUtils.h>
#include <AliInputEventHandler.h>
#include <AliVEvent.h>
#include <AliVVertex.h>

#include "AliAnalysisTaskEventFilter.h"

using namespace EMCalTriggerPtAnalysis;

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEventFilter)

AliAnalysisTaskEventFilter::AliAnalysisTaskEventFilter() :
AliAnalysisTaskSE(),
fAnalysisUtils(NULL),
fHistos(NULL)
{
}

AliAnalysisTaskEventFilter::AliAnalysisTaskEventFilter(const char *name) :
AliAnalysisTaskSE(name),
fAnalysisUtils(NULL),
fHistos(NULL)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEventFilter::~AliAnalysisTaskEventFilter() {
  delete fAnalysisUtils;
  delete fHistos;
}

void AliAnalysisTaskEventFilter::UserCreateOutputObjects(){
  fAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager("eventHistos");
  fHistos->CreateTH1("hINT7all", "All INT7 events", 1000, -100, 100);                     // all INT7 events
  fHistos->CreateTH1("hINT7psel", "Physics selected INT7 events", 1000, -100, 100);       // after Physics selection
  fHistos->CreateTH1("hINT7fe", "INT7 no first event in chunk", 1000, -100, 100);         // Not first event in chunk
  fHistos->CreateTH1("hINT7np", "INT7 no pileup event", 1000, -100, 100);                 // Not a pileup event
  fHistos->CreateTH1("hINT7fv", "INT7 after fake vertex cut", 1000, -100, 100);           // survived fake vertex cut
  fHistos->CreateTH1("hINT7tv", "INT7 after true vertex cut", 1000, -100, 100);           // survived true vertex cut

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEventFilter::UserExec(Option_t *){
  PostData(1, fHistos->GetListOfHistograms());
  // Filter event, also run buggy event selection explicitly to check impact
  double vz = InputEvent()->GetPrimaryVertex()->GetZ();
  TString triggerstring(InputEvent()->GetFiredTriggerClasses());
  if(!triggerstring.Contains("INT7-B")) return;
  fHistos->FillTH1("hINT7all", vz);
  if(!(fInputHandler->IsEventSelected() && AliVEvent::kINT7)) return;
  fHistos->FillTH1("hINT7psel", vz);
  if(!fAnalysisUtils->IsFirstEventInChunk(InputEvent())) return;
  fHistos->FillTH1("hINT7fe", vz);
  if(!fAnalysisUtils->IsPileUpEvent(InputEvent())) return;
  fHistos->FillTH1("hINT7np", vz);
  if(FakeVertexSelection2013pA(InputEvent())) return;
  fHistos->FillTH1("hINT7fv", vz);
  if(fAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return;
  fHistos->FillTH1("hINT7tv", vz);
}

Bool_t AliAnalysisTaskEventFilter::FakeVertexSelection2013pA(const AliVEvent *const inputEvent) const {

  Bool_t accept = kFALSE;

  const AliVVertex *trkVtx = dynamic_cast<const AliVVertex*>(inputEvent->GetPrimaryVertex()) ;
  if(!trkVtx || trkVtx->GetNContributors() <= 0){
    accept = kFALSE;
    return accept;
  }

  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return accept;

  Float_t zvtx = trkVtx->GetZ();
  const AliVVertex* spdVtx = dynamic_cast<const AliVVertex*>(inputEvent->GetPrimaryVertexSPD()) ;
  if (spdVtx->GetNContributors()<=0) return accept;

  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (TString(spdVtx->GetTitle()).Contains("vertexer:Z") && (zRes>0.25)) return accept;        // Explicitly done the old way in order to estimate effect of not applied cut
  if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return accept;

  if (TMath::Abs(zvtx) > 10) return accept;

  return kTRUE;
}
