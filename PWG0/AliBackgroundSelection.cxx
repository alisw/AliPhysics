#include "AliBackgroundSelection.h"
#include "TH2F.h"
#include "TList.h"
#include "AliLog.h"
#include "TString.h"
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisManager.h"
#include "TTree.h"
#include "../ITS/AliITSRecPoint.h"
#include "AliMultiplicity.h"

ClassImp(AliBackgroundSelection)

AliBackgroundSelection::AliBackgroundSelection():
  AliAnalysisCuts(), fOutputHist(0), fACut(0), fBCut(0)
{
  
  fOutputHist = new TList();
  fOutputHist->SetOwner();
  fACut = 45;
  fBCut = 6.5;
  

}

AliBackgroundSelection::AliBackgroundSelection(const char* name, const char* title):
  AliAnalysisCuts(name,title), fOutputHist(0), fACut(0), fBCut(0)
{

  fOutputHist = new TList();
  fOutputHist->SetOwner();
  fACut = 45;
  fBCut = 6.5;
}

AliBackgroundSelection::AliBackgroundSelection(const AliBackgroundSelection& obj) : AliAnalysisCuts(obj),
fOutputHist(0), fACut(0), fBCut(0)
{

  fOutputHist = obj.fOutputHist;
  fACut       = obj.fACut;
  fBCut       = obj.fBCut;

}

AliBackgroundSelection::~AliBackgroundSelection() {
  if(fOutputHist) delete fOutputHist;

}

Bool_t AliBackgroundSelection::IsSelected(TObject* obj){

  // reset fSelected
  SetSelected(kFALSE);
  // Get rec points
  AliESDInputHandlerRP* handlerRP = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handlerRP)
    AliFatal("Cannot get the AliESDInputHandlerRP");

  TTree* itsClusterTree = handlerRP->GetTreeR("ITS");
  if (!itsClusterTree){
    AliError("Cannot get the ITS Cluster tree");
    return kFALSE;
  }
  //    AliFatal("Cannot get the ITS Cluster tree");

  TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
  TBranch* itsClusterBranch=itsClusterTree->GetBranch("ITSRecPoints");

  itsClusterBranch->SetAddress(&itsClusters);

  Int_t nItsSubs = (Int_t)itsClusterTree->GetEntries();


  AliESDEvent * esdEv = (AliESDEvent*) obj;

  // Get # spd clusters and of tracklets
  Int_t spdClusters=0;


  // loop over the its subdetectors
  for (Int_t iIts=0; iIts < nItsSubs; iIts++) {

    if (!itsClusterTree->GetEvent(iIts))
      continue;

    Int_t nClusters = itsClusters->GetEntriesFast();

    // loop over clusters
    while (nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*) itsClusters->UncheckedAt(nClusters);

      Int_t layer = cluster->GetLayer();

      if (layer < 3) { // SPD
	spdClusters++;
      }
    }
  }

  const AliMultiplicity* mult = esdEv->GetMultiplicity();
  if (!mult){
    AliFatal("No multiplicity object"); // TODO: Should this be fatal?
  }
  Int_t ntracklet = mult->GetNumberOfTracklets();

  Float_t limit = fACut + ntracklet * fBCut;
  
  if (spdClusters > limit) SetSelected(kFALSE);
  else                     SetSelected(kTRUE );


  // Fill control histos for all trigger classes
  TString trgClasses = esdEv->GetFiredTriggerClasses();
  TObjArray * tokens = trgClasses.Tokenize(" ");
  TIter iter(tokens);
  while(TObjString * tok = (TObjString*) iter.Next()){
    TString trg = tok->GetString();
    trg.Strip(TString::kTrailing, ' ');
    trg.Strip(TString::kLeading, ' ');
    //    Printf("TRG: [%s]\n",trg.Data());
    TH2F * hCvsT = GetClusterVsTrackletsHisto(trg.Data());
    if(!hCvsT) {
      // if histo does not exist, book it on the fly (also books accepted histo)
      BookClusterVsTrackletsHisto(trg.Data());
      hCvsT = GetClusterVsTrackletsHisto(trg.Data());
    }
    TH2F * hCvsTa = GetClusterVsTrackletsHistoAccepted(trg.Data());
    hCvsT->Fill(ntracklet,spdClusters);
    if(Selected()) hCvsTa->Fill(ntracklet,spdClusters);

  }
  if(tokens) delete tokens;
  // return decision

  return Selected();
}


void   AliBackgroundSelection::Init(){

  fACut = 45;
  fBCut = 6.5;

}


void AliBackgroundSelection::BookClusterVsTrackletsHisto(const char * trigger_name){

  TH2F * h1 = new TH2F(GetClusterVsTrackletsHistoName(trigger_name),trigger_name, 50, -0.5, 49.5, 1000, -0.5, 999.5);
  h1->SetXTitle("Tracklets");
  h1->SetYTitle("SPD Clusters");
  AliInfo(Form("Creating histos: %s, all and accepted", GetClusterVsTrackletsHistoName(trigger_name)));

  TH2F * h2 = new TH2F(GetClusterVsTrackletsHistoNameAccepted(trigger_name),TString(trigger_name)+ "(accepted)", 
		       50, -0.5, 49.5, 1000, -0.5, 999.5);
  h2->SetXTitle("Tracklets");
  h2->SetYTitle("SPD Clusters");

  fOutputHist->Add(h1);
  fOutputHist->Add(h2);
}

TH2F * AliBackgroundSelection::GetClusterVsTrackletsHisto(const char * trigger_name){
  if(!fOutputHist) {AliError("List of histos not initialized");return 0;}
  return (TH2F*) fOutputHist->FindObject(GetClusterVsTrackletsHistoName(trigger_name));
  
}

TH2F * AliBackgroundSelection::GetClusterVsTrackletsHistoAccepted(const char * trigger_name){

  if(!fOutputHist) {AliError("List of histos not initialized");return 0;}
  return (TH2F*) fOutputHist->FindObject(GetClusterVsTrackletsHistoNameAccepted(trigger_name));
  
}

const char * AliBackgroundSelection::GetClusterVsTrackletsHistoName(const char * trigger_name){
    static TString str;
    str = ("h");
    str = str+GetName()+"_"+trigger_name;
    return str.Data();
}

const char * AliBackgroundSelection::GetClusterVsTrackletsHistoNameAccepted(const char * trigger_name){
    static TString str;
    str = ("h");
    str = str+GetName()+"_"+trigger_name + "_accepted";
    return str.Data();
}

Long64_t AliBackgroundSelection::Merge(TCollection* list)
{
  // Merge a list of AliBackgroundSelection objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  const Int_t nHists = 1;
  TList collections[nHists];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliBackgroundSelection* entry = dynamic_cast<AliBackgroundSelection*> (obj);
    if (entry == 0) 
      continue;

    Int_t n = 0;
    collections[n++].Add(entry->fOutputHist);

    count++;
  }

  Int_t n = 0;
  fOutputHist->Merge(&collections[n++]);
  
  delete iter;

  return count+1;
}


