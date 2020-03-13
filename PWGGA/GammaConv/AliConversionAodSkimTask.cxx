#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile.h>
#include <TSystem.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include <AliAODMCHeader.h>
#include <AliAODMCParticle.h>
#include <AliAODConversionPhoton.h>
#include <AliAODVertex.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>
#include "AliConversionAodSkimTask.h"
#include "TObjectTable.h"
using namespace std;
ClassImp(AliConversionAodSkimTask)

AliConversionAodSkimTask::AliConversionAodSkimTask(const char* name) :
  AliAodSkimTask(name), fConvMinPt(-1), fConvMinEta(-999.), fConvMaxEta(999.), fConvMinPhi(999.), fConvMaxPhi(999.),fDoBothConvPtAndAcc(0),
  fDoQA(0),fHconvPtBeforeCuts(0), fHconvPtAfterCuts(0), fHconvAccBeforeCuts(0), fHconvAccAfterCuts(0)
{
  if (name) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
  }
}  

AliConversionAodSkimTask::~AliConversionAodSkimTask()
{
  if (fOutputList) {
    delete fOutputList;
  }
  delete fHevs;
  delete fHclus;
  delete fHtrack;
  delete fHconvPtBeforeCuts;
  delete fHconvPtAfterCuts;
  delete fHconvAccBeforeCuts;
  delete fHconvAccAfterCuts;
}

void AliConversionAodSkimTask::UserCreateOutputObjects()
{
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler *oh = (AliAODHandler*)man->GetOutputEventHandler();
  if (oh) {
    TFile *fout = oh->GetTree()->GetCurrentFile();
    fout->SetCompressionLevel(2);
  }

  // Before running skim, check if the V0Reader is running
  TIter next(man->GetTasks());
  TObject *obj;
  while ((obj = next()))
  {
    if (strcmp(obj->ClassName(), "AliV0ReaderV1") == 0)
    {
      AliFatal("Skim task is running but detected V0Reader running in addition! This will cause problems due to relabelling already done now by the V0Reader!");
    }
  }

  // check if convcuts were requested but no gammabranch given
  if ((fConvMinPt > 0) || (fConvMinEta != -999) || (fConvMaxEta != 999) || (fConvMinPhi != -999) || (fConvMaxPhi != 999)){
     if(!fGammaBr || !fGammaBr.CompareTo("")){
       AliFatal(Form("%s: Conversion cuts were requested but no fGammaBr was given", GetName()));
       return;
     }
  }

  fOutputList = new TList;
  fOutputList->SetOwner();
  fHevs = new TH1F("hEvs","",2,-0.5,1.5);
  fOutputList->Add(fHevs);
  fHclus = new TH1F("hClus",";E (GeV)",200,0,100);
  fOutputList->Add(fHclus);
  fHtrack = new TH1F("hTrack",";p_{T} (GeV/c)",200,0,100);
  fOutputList->Add(fHtrack);
  if(fDoQA){
    fHconvPtBeforeCuts = new TH1F("fHconvPtBeforeCuts",";p_{T} (GeV/c)",200,0,100);
    fOutputList->Add(fHconvPtBeforeCuts);
    fHconvPtAfterCuts = new TH1F("fHconvPtAfterCuts",";p_{T} (GeV/c)",200,0,100);
    fOutputList->Add(fHconvPtAfterCuts);
    fHconvAccBeforeCuts = new TH2F("fHconvAccBeforeCuts",";#eta; #phi",100,-0.9,0.9,100,0.,2*TMath::Pi());
    fOutputList->Add(fHconvAccBeforeCuts);
    fHconvAccAfterCuts = new TH2F("fHconvAccAfterCuts",";#eta; #phi",100,-0.9,0.9,100,0.,2*TMath::Pi());
    fOutputList->Add(fHconvAccAfterCuts);
  }
  PostData(1, fOutputList);
}

Bool_t AliConversionAodSkimTask::SelectEvent()
{
  // Accept event only if an EMCal-cluster with a minimum energy of fClusMinE has been found in the event
  Bool_t storeE = kFALSE;
  if (fClusMinE>0) {
    TClonesArray *cls  = fAOD->GetCaloClusters();
    for (Int_t i=0; i<cls->GetEntriesFast(); ++i) {
      AliAODCaloCluster *clus = static_cast<AliAODCaloCluster*>(cls->At(i));
      if (!clus->IsEMCAL())
        continue;
      Double_t e = clus->E();
      fHclus->Fill(e);
      if (e>fClusMinE) {
        storeE = kTRUE;
      }
    }
  } else {
    storeE = kTRUE;
  }

  // Accept event only if an track with a minimum pT of fTrackMinPt has been found in the event
  Bool_t storePt = kFALSE;
  if (fTrackMinPt > 0 ){
    TClonesArray *tracks = fAOD->GetTracks();
    for (Int_t i=0;i<tracks->GetEntries();++i) {
      AliAODTrack *t = static_cast<AliAODTrack*>(tracks->At(i));
      Double_t pt = t->Pt();
      fHtrack->Fill(pt);
      if (pt>fTrackMinPt) {
        storePt = kTRUE;
      }
      if(fTrackMaxPt > 0  &&  fTrackMaxPt < pt){
        storePt = kFALSE;
      } 
    }
  } else {
    storePt = kTRUE;
  }

  // Accept event only if a conversion with a minimum pT of fConvMinPt has been found in the event 
  Bool_t storeConvPt = kFALSE;
  Bool_t storeConvAcc = kFALSE; 
  if ((fConvMinPt > 0) || (fConvMinEta != -999) || (fConvMaxEta != 999) || (fConvMinPhi != -999) || (fConvMaxPhi != 999)){
    // This is only done if any ConvPt or ConvAcceptance cut was given by user    
    TClonesArray *convgammas  = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fGammaBr));
    if(convgammas){
      for(Int_t i=0;i<convgammas->GetEntriesFast();i++){
        AliAODConversionPhoton* g =dynamic_cast<AliAODConversionPhoton*>(convgammas->At(i));
        if(g){
           Double_t pt = g->Pt();
           Double_t eta = g->Eta();
           Double_t phi = g->Phi();
           if(fDoQA){
            fHconvPtBeforeCuts->Fill(pt);
            fHconvAccBeforeCuts->Fill(eta,phi);
           }
           if(fConvMinPt > 0){
              if (pt>fConvMinPt){ // pt cut
                  storeConvPt = kTRUE;
              }
           } else{
              storeConvPt = kTRUE;
           }
           if( (eta>fConvMinEta) && (eta<fConvMaxEta) && (phi>fConvMinPhi) && (phi>fConvMaxPhi)){ // will always be true if no cut was set
             storeConvAcc = kTRUE;
           }
        }
      }
    }
  } else { // neither pt cuts nor acceptance cuts were set
    storeConvPt = kTRUE;
    storeConvAcc = kTRUE;
  }

  Bool_t store = kFALSE;
  if (fDoBothMinTrackAndClus && fClusMinE>0 && fTrackMinPt > 0){
    // request that both conditions are full-filled for propagating the event
    store     = (storeE && storePt);
  } else if (!fDoBothMinTrackAndClus && fClusMinE>0 && fTrackMinPt > 0){
    // request that at least one of the conditions is fullfilled
    store     = (storeE || storePt);
  } else if (fDoBothConvPtAndAcc && (fConvMinPt >0) && (fConvMinEta!=-999 || fConvMaxEta!=999 || fConvMinPhi!=-999 || fConvMaxPhi!=999)){
    store     = (storeConvPt && storeConvAcc);
  } else if (!fDoBothConvPtAndAcc && (fConvMinPt >0) && (fConvMinEta!=-999 || fConvMaxEta!=999 || fConvMinPhi!=-999 || fConvMaxPhi!=999)){
    store     = (storeConvPt || storeConvAcc);
  } else if ( fClusMinE>0 ){
    store     = storeE;
  } else if ( fTrackMinPt>0 ){
    store     = storePt;
  } else if ( fConvMinPt>0 ){
    store     = storeConvPt;
  } else if (fConvMinEta!=-999 || fConvMaxEta!=999 || fConvMinPhi!=-999 || fConvMaxPhi!=999){
    store     = storeConvAcc;
  } else {
    store     = kTRUE;
  }

  if(fDoQA){
    TClonesArray *convgammas  = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fGammaBr));
    if(convgammas){
      for(Int_t i=0;i<convgammas->GetEntriesFast();i++){
        AliAODConversionPhoton* g =dynamic_cast<AliAODConversionPhoton*>(convgammas->At(i));
        if(g){
            Double_t pt = g->Pt();
            Double_t eta = g->Eta();
            Double_t phi = g->Phi();
            fHconvPtAfterCuts->Fill(pt);
            fHconvAccAfterCuts->Fill(eta,phi);
        }
      }
    }
  }
  return store;
}

Bool_t AliConversionAodSkimTask::UserNotify()
{
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();
  Int_t nevents = tree->GetEntriesFast();

  Int_t slevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  Bool_t res = PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);
  gErrorIgnoreLevel=slevel;

  if (res) {
    cout << "AliConversionAodSkimTask " << GetName() << " found xsec info: " << xsection << " " << trials << " " << pthardbin << " " << nevents << endl;
    fPyxsec      = xsection;
    fPytrials    = trials;
    fPypthardbin = pthardbin;
  }

  return res;
}

void AliConversionAodSkimTask::Terminate(Option_t *)
{
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   if (man->GetAnalysisType()!=0)
     return;
   if (fHevs==0)
     return;
   Int_t norm = fHevs->GetEntries();
   if (norm<1)
     norm=1;
   cout << "AliConversionAodSkimTask " << GetName() << " terminated with accepted fraction of events: " << fHevs->GetBinContent(2)/norm
	<< " (" << fHevs->GetBinContent(2) << "/" << fHevs->GetEntries() << ")" << endl;
}

const char *AliConversionAodSkimTask::Str() const
{
  return Form("mine%.2f_%dycut%.2f_%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
              fClusMinE,
              fCutMC,
              fYCutMC,
              fDoCopyHeader,
              fDoCopyVZERO,
              fDoCopyTZERO,
              fDoCopyVertices,
              fDoCopyTOF,
              fDoCopyTracklets,
              fDoCopyTracks,
              fDoCopyTrigger,
              fDoCopyPTrigger,
              fDoCopyCells,
              fDoCopyPCells,
              fDoCopyClusters,
              fDoCopyDiMuons,
              fDoCopyTrdTracks,
              fDoCopyV0s,
              fDoCopyCascades,
              fDoCopyZDC,
              fDoCopyConv,
              fDoCopyMC,
              fDoCopyMCHeader);
}

