/*
 * AliAnalysisTaskEventFilter.cxx
 *
 *  Created on: Jan 29, 2016
 *      Author: markus
 */
#include <TArray.h>
#include <THashList.h>
#include <THistManager.h>
#include <TMath.h>
#include <TString.h>

#include <AliAnalysisUtils.h>
#include <AliESDtrackCuts.h>
#include <AliESDtrack.h>
#include <AliInputEventHandler.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVVertex.h>

#include <map>
#include <string>

#include "AliAnalysisTaskEventFilter.h"

using namespace PWGJE::EMCALJetTasks;

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEventFilter)

AliAnalysisTaskEventFilter::AliAnalysisTaskEventFilter() :
AliAnalysisTaskSE(),
fAnalysisUtils(NULL),
fTrackCuts(NULL),
fHistos(NULL)
{
}

AliAnalysisTaskEventFilter::AliAnalysisTaskEventFilter(const char *name) :
AliAnalysisTaskSE(name),
fAnalysisUtils(NULL),
fTrackCuts(NULL),
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

  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  fTrackCuts->SetName("Standard Track cuts");
  fTrackCuts->SetMinNCrossedRowsTPC(120);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

  std::map<std::string, std::string> filtersteps;
  filtersteps.insert(std::pair<std::string,std::string>("all", "All INT7 events"));                 // all INT7 events
  filtersteps.insert(std::pair<std::string,std::string>("psel", "Physics selected INT7 events"));   // after Physics selection
  filtersteps.insert(std::pair<std::string,std::string>("fe", "INT7 no first event in chunk"));     // Not first event in chunk
  filtersteps.insert(std::pair<std::string,std::string>("np", "INT7 no pileup event"));             // Not a pileup event
  filtersteps.insert(std::pair<std::string,std::string>("fv", "INT7 after fake vertex cut"));       // survived fake vertex cut
  filtersteps.insert(std::pair<std::string,std::string>("tv", "INT7 after true vertex cut"));       // survived true vertex cut

  TArrayD ptbinning;
  CreatePtBinning(ptbinning);

  fHistos = new THistManager("eventHistos");
  for(std::map<std::string, std::string>::iterator it = filtersteps.begin(); it != filtersteps.end(); ++it){
    fHistos->CreateTH1(Form("hINT7%s", it->first.c_str()), it->second.c_str(), 1000, -100, 100);
  }
  for(std::map<std::string, std::string>::iterator it = filtersteps.begin(); it != filtersteps.end(); ++it){
    fHistos->CreateTH1(Form("hPt%s", it->first.c_str()), it->second.c_str(), ptbinning);
  }
  for(std::map<std::string, std::string>::iterator it = filtersteps.begin(); it != filtersteps.end(); ++it){
    fHistos->CreateTH1(Form("hEta%s", it->first.c_str()), it->second.c_str(), 100, -0.8, 0.8);
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEventFilter::UserExec(Option_t *){
  std::map<int, std::string> eventfilters;
  eventfilters.insert(std::pair<int, std::string>(0, "psel"));
  eventfilters.insert(std::pair<int, std::string>(1, "fe"));
  eventfilters.insert(std::pair<int, std::string>(2, "np"));
  eventfilters.insert(std::pair<int, std::string>(3, "fv"));
  eventfilters.insert(std::pair<int, std::string>(4, "tv"));

  PostData(1, fHistos->GetListOfHistograms());
  // Filter event, also run buggy event selection explicitly to check impact
  double vz = InputEvent()->GetPrimaryVertex()->GetZ();
  TString triggerstring(InputEvent()->GetFiredTriggerClasses());
  if(!triggerstring.Contains("INT7-B")) return;

  UChar_t eventfilter = FilterEvent();
  std::vector<const AliVTrack *> tracks = FilterTracks();
  FillEvent("all", vz);
  FillTracks("all", tracks);
  for(int ifilter = 0; ifilter < 5; ifilter++){
    if(eventfilter & 1 << ifilter){
      std::string filtername = eventfilters.find(ifilter)->second;
      FillEvent(filtername.c_str(), vz);
      FillTracks(filtername.c_str(), tracks);
    }
  }
}

UChar_t AliAnalysisTaskEventFilter::FilterEvent() const {
  UChar_t filterbits = 0;
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return filterbits;
  filterbits |= 1 << 0;
  if(fAnalysisUtils->IsFirstEventInChunk(InputEvent())) return filterbits;
  filterbits |= 1 << 1;
  if(fAnalysisUtils->IsPileUpEvent(InputEvent())) return filterbits;
  filterbits |= 1 << 2;
  if(!FakeVertexSelection2013pA(InputEvent())) return filterbits;
  filterbits |= 1 << 3;
  if(!fAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return filterbits;
  filterbits |= 1 << 4;
  return filterbits;
}

std::vector<const AliVTrack *> AliAnalysisTaskEventFilter::FilterTracks() const {
  std::vector<const AliVTrack *> result;
  AliESDtrack *track(NULL);
  for(Int_t itrk = 0; itrk < InputEvent()->GetNumberOfTracks(); itrk++){
    track = dynamic_cast<AliESDtrack *>(InputEvent()->GetTrack(itrk));
    if(!track) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    if(!fTrackCuts->AcceptTrack(track)) continue;
    result.push_back(track);
  }
  return result;
}

void AliAnalysisTaskEventFilter::FillEvent(const char *filterstep, double vz){
  fHistos->FillTH1(Form("hINT7%s", filterstep), vz);
}

void AliAnalysisTaskEventFilter::FillTracks(const char *filterstep, const std::vector<const AliVTrack *> &tracks){
  for(std::vector<const AliVTrack *>::const_iterator trackiter = tracks.begin(); trackiter != tracks.end(); ++trackiter){
    fHistos->FillTH1(Form("hPt%s", filterstep), TMath::Abs((*trackiter)->Pt()));
    fHistos->FillTH1(Form("hEta%s", filterstep), (*trackiter)->Eta());
  }
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

void AliAnalysisTaskEventFilter::CreatePtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

