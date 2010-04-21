// ----------------------------------------------------------------
// AliBackgroundSelection
//
// This class implements to cuts to reject background events from the
// samples to be used in the physics analysis:
// 1. A linear cut on the correlation cluster vs tracklets
// 2. A cut on the delta phi window used by the vertexer Z
// The parameters used in both cuts can be set
// 
// The class also produces control histograms for all and accepted
// events, for each trigger class present in the data independently.
// Histograms are booked on the fly in the UserExec, whenever a new
// trigger class is found.
//
// After the first implementation, it was realized that the deltaphi
// cut is more a quality selection cut than an event selection cut, so
// it is effectively disabled by default.
//
// Author: Michele Floris, CERN
// ----------------------------------------------------------------


#include "AliBackgroundSelection.h"
#include "TH2F.h"
#include "TList.h"
#include "TString.h"
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisManager.h"
#include "TTree.h"
#include "AliMultiplicity.h"
#ifdef PASS1RECO
#include "../ITS/AliITSRecPoint.h"
#endif



ClassImp(AliBackgroundSelection)

AliBackgroundSelection::AliBackgroundSelection():
  AliAnalysisCuts(), fOutputHist(0), fACut(0), fBCut(0), fDeltaPhiCut(10)
{
  // ctor
  fOutputHist = new TList();
  fOutputHist->SetOwner();
  fACut = 65;
  fBCut = 4;
  fDeltaPhiCut = 10; // effectively disabling delta phi cut by default
}

AliBackgroundSelection::AliBackgroundSelection(const char* name, const char* title):
  AliAnalysisCuts(name,title), fOutputHist(0), fACut(0), fBCut(0), fDeltaPhiCut(10)
{
  // ctor
  fOutputHist = new TList();
  fOutputHist->SetOwner();
  fACut = 65;
  fBCut = 4;
  fDeltaPhiCut = 10; //  effectively disabling delta phi cut by default

}

AliBackgroundSelection::AliBackgroundSelection(const AliBackgroundSelection& obj) : AliAnalysisCuts(obj),
fOutputHist(0), fACut(0), fBCut(0), fDeltaPhiCut(0)
{
  // copy ctor
  fOutputHist  = obj.fOutputHist;
  fACut        = obj.fACut;
  fBCut        = obj.fBCut;
  fDeltaPhiCut = obj.fDeltaPhiCut;
}

AliBackgroundSelection::~AliBackgroundSelection() {
  // dtor
  if(fOutputHist) {
    delete fOutputHist;
    fOutputHist = 0;
  }

}

Bool_t AliBackgroundSelection::IsSelected(TObject* const obj) 
{
  // returns false if the event is identifiead as beam background,
  // true otherwise.

  // reset fSelected
  SetSelected(kFALSE);
#ifdef PASS1RECO
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
#endif

  AliESDEvent * esdEv = (AliESDEvent*) obj;

#ifdef PASS1RECO
  Float_t deltaPhi = 0.0; // deltaPhi is not available in pass1

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
#endif

  const AliMultiplicity* mult = esdEv->GetMultiplicity();
  if (!mult){
    AliFatal("No multiplicity object"); // TODO: Should this be fatal?
  }
  Int_t ntracklet = mult->GetNumberOfTracklets();

#ifndef PASS1RECO
  // get deltaphi if vertexer z
  Float_t deltaPhi = 0.0;
  // Get Vertex
  const AliESDVertex * vtxESD = esdEv->GetPrimaryVertexSPD();
  if(vtxESD) {
    if (vtxESD->IsFromVertexerZ()) deltaPhi = vtxESD->GetDispersion(); // dispersion contains deltaphi in case of vertexer Z
  }
  else {
    AliWarning("No Vertex");
  }

  

  // compute number of spd clusters
  Float_t spdClusters = 0;
  for(Int_t ilayer = 0; ilayer < 2; ilayer++){
    spdClusters += mult->GetNumberOfITSClusters(ilayer);
  }
#endif

  // Check cuts
  Bool_t isCvsTOk     = kFALSE;
  Bool_t isDeltaPhiOk = kFALSE;

  Float_t limit = fACut + ntracklet * fBCut;  
  if (spdClusters > limit)        isCvsTOk = kFALSE;
  else                            isCvsTOk = kTRUE ;

  if(deltaPhi > fDeltaPhiCut)     isDeltaPhiOk = kFALSE;
  else                            isDeltaPhiOk = kTRUE ;

  if (!isCvsTOk || !isDeltaPhiOk) SetSelected(kFALSE);
  else                            SetSelected(kTRUE );

  // Fill control histos for all trigger classes
  TString trgClasses = esdEv->GetFiredTriggerClasses();
  TObjArray * tokens = trgClasses.Tokenize(" ");
  TIter iter(tokens);
  while(TObjString * tok = (TObjString*) iter.Next()){
    // clean up trigger name
    TString trg = tok->GetString();
    trg.Strip(TString::kTrailing, ' ');
    trg.Strip(TString::kLeading, ' ');
    
    // cluster vs tracklets
    TH2F * hCvsT = GetClusterVsTrackletsHisto(trg.Data());
    TH2F * hCvsTa = GetClusterVsTrackletsHistoAccepted(trg.Data());
    hCvsT->Fill(ntracklet,spdClusters);
    if(isCvsTOk) hCvsTa->Fill(ntracklet,spdClusters);

    // Delta phi
    TH1F * hDeltaPhi = GetDeltaPhiHisto(trg.Data());
    TH1F * hDeltaPhia = GetDeltaPhiHistoAccepted(trg.Data());
    hDeltaPhi->Fill(deltaPhi);
    if(isDeltaPhiOk) hDeltaPhia->Fill(deltaPhi);
  }
  if(tokens) delete tokens;
  // return decision

#ifdef PASS1RECO
  if(itsClusters) {
    itsClusters->Delete();
    delete itsClusters;
  }
#endif 
  return Selected();
}


void   AliBackgroundSelection::Init(){

  // Set default cut values
  fACut = 65;
  fBCut = 4;

}


void AliBackgroundSelection::BookClusterVsTrackletsHisto(const char * trigger_name){

  // Book control histogram for the cut on the correlation cluster vs tracklets

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

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

  TH1::AddDirectory(oldStatus);

}

void AliBackgroundSelection::BookDeltaPhiHisto(const char * trigger_name){

  // Book control histogram for the cut on the DeltaPhi window used by vertexer Z

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1F * h1 = new TH1F(GetDeltaPhiHistoName(trigger_name),trigger_name, 100,0,0.5);
  h1->SetXTitle("#Delta #phi");
  AliInfo(Form("Creating histos: %s, all and accepted", GetDeltaPhiHistoName(trigger_name)));

  TH1F * h2 = new TH1F(GetDeltaPhiHistoNameAccepted(trigger_name),TString(trigger_name)+ "(accepted)", 100,0,0.5);
  h2->SetXTitle("#Delta #phi");


  fOutputHist->Add(h1);
  fOutputHist->Add(h2);

  TH1::AddDirectory(oldStatus);

}

TH2F * AliBackgroundSelection::GetClusterVsTrackletsHisto(const char * trigger_name){

  // Returns the control histogram corresponding to a given trigger
  // class. If it does not exist, it creates it and adds it to the
  // output list
  // All Events

  if(!fOutputHist) {AliError("List of histos not initialized");return 0;}
  TH2F * h = (TH2F*) fOutputHist->FindObject(GetClusterVsTrackletsHistoName(trigger_name));  
  if(!h) {
    BookClusterVsTrackletsHisto(trigger_name);
    h = (TH2F*) fOutputHist->FindObject(GetClusterVsTrackletsHistoName(trigger_name));  
  }
  return h;
}
TH1F * AliBackgroundSelection::GetDeltaPhiHisto(const char * trigger_name){

  // Returns the control histogram corresponding to a given trigger
  // class. If it does not exist, it creates it and adds it to the
  // output list
  // All Events

  if(!fOutputHist) {AliError("List of histos not initialized");return 0;}
  TH1F * h = (TH1F*) fOutputHist->FindObject(GetDeltaPhiHistoName(trigger_name));  
  if(!h) {
    BookDeltaPhiHisto(trigger_name);
    h  = (TH1F*) fOutputHist->FindObject(GetDeltaPhiHistoName(trigger_name));  
  }
  return h;
}

TH2F * AliBackgroundSelection::GetClusterVsTrackletsHistoAccepted(const char * trigger_name){

  // Returns the control histogram corresponding to a given trigger
  // class. If it does not exist, it creates it and adds it to the
  // output list
  // Events passing the cut only

  if(!fOutputHist) {AliError("List of histos not initialized");return 0;}
  TH2F * h = (TH2F*) fOutputHist->FindObject(GetClusterVsTrackletsHistoNameAccepted(trigger_name));
  if(!h) {
    BookClusterVsTrackletsHisto(trigger_name);
    h = (TH2F*) fOutputHist->FindObject(GetClusterVsTrackletsHistoNameAccepted(trigger_name));  
  }
  return h;
  
}

TH1F * AliBackgroundSelection::GetDeltaPhiHistoAccepted(const char * trigger_name){

  // Returns the control histogram corresponding to a given trigger
  // class. If it does not exist, it creates it and adds it to the
  // output list
  // Events passing the cut only

  if(!fOutputHist) {AliError("List of histos not initialized");return 0;}
  TH1F * h = (TH1F*) fOutputHist->FindObject(GetDeltaPhiHistoNameAccepted(trigger_name));  
  if(!h) {
    BookDeltaPhiHisto(trigger_name);
    h  = (TH1F*) fOutputHist->FindObject(GetDeltaPhiHistoNameAccepted(trigger_name));  
  }
  return h;
  
}

const char * AliBackgroundSelection::GetClusterVsTrackletsHistoName(const char * trigger_name){

  // build up the name of the cluster vs tracklets histo using the trigger class

    static TString str;
    str = ("hCvsT");
    str = str+GetName()+"_"+trigger_name;
    return str.Data();
}

const char * AliBackgroundSelection::GetClusterVsTrackletsHistoNameAccepted(const char * trigger_name){

  // build up the name of the cluster vs tracklets histo using the trigger class (accepted events)
    static TString str;
    str = ("hCvsT");
    str = str+GetName()+"_"+trigger_name + "_accepted";
    return str.Data();
}

const char * AliBackgroundSelection::GetDeltaPhiHistoName(const char * trigger_name){

  // build up the name of the delta phi histo using the trigger class


    static TString str;
    str = ("hDeltaPhi");
    str = str+GetName()+"_"+trigger_name;
    return str.Data();
}

const char * AliBackgroundSelection::GetDeltaPhiHistoNameAccepted(const char * trigger_name){

  // build up the name of the delta phi histo using the trigger class (accepted events)

    static TString str;
    str = ("hDeltaPhi");
    str = str+GetName()+"_"+trigger_name + "_accepted";
    return str.Data();
}

Long64_t AliBackgroundSelection::Merge(TCollection* const list)
{
  // Merge a list of AliBackgroundSelection objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  // We have to make sure that all the list contain the same histos in
  // the same order. We thus also have to sort the list (sorting is
  // done by name in TList).

  AliInfo("Merging");

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
  // 1. Sort this list
  fOutputHist->Sort();
  
  while ((obj = iter->Next())) {
    Bool_t foundDiffinThisIterStep = kFALSE;
    //    Printf("%d - %s",count, obj->GetName());
    AliBackgroundSelection* entry = dynamic_cast<AliBackgroundSelection*> (obj);
    if (entry == 0) 
      continue;

    TList * hlist = entry->fOutputHist;

    // Check if all histos in this fOutputHist are also in the one from entry and viceversa
    // Use getters to automatically book non defined histos    

    Bool_t areListsDifferent=kTRUE;
    Int_t iloop = 0;
    Int_t maxLoops = hlist->GetSize() + fOutputHist->GetSize(); // In the worst case all of the histos will be different...    
    while(areListsDifferent) {
      if(iloop>maxLoops) AliFatal("Infinite Loop?");
      iloop++;
      // sort
      hlist->Sort();
      fOutputHist->Sort();
      // loop over the largest 

      // loop over the largest 
      TObject * hist =0;
      TIterator * iterlist = 0;
      TList * thislist  = 0; // the list over which I'm iterating (i.e. the largest)
      TList * otherlist = 0; // the other list

      if (hlist->GetSize() >= fOutputHist->GetSize()) { 
	thislist  = hlist;
	otherlist = fOutputHist;
      }
      else{
	thislist  = fOutputHist;
	otherlist = hlist;	
      }
      iterlist = thislist->MakeIterator();

      while ((hist= iterlist->Next())){ 
	if(!otherlist->FindObject(hist->GetName())){
	  AliInfo(Form("Adding object %s",hist->GetName()));
	  foundDiffinThisIterStep = kTRUE;
	  TH1 * hclone =  (TH1*) hist->Clone();
	  hclone->Reset();
	  otherlist->Add(hclone);
	}
      }

      // re-sort before checking
      hlist->Sort();
      fOutputHist->Sort();

      // check if everything is fine    
      areListsDifferent=kFALSE;
      if (hlist->GetSize() == fOutputHist->GetSize()) {	
	Int_t nhist =  fOutputHist->GetSize();
	for(Int_t ihist = 0; ihist < nhist; ihist++){
	  if(strcmp(fOutputHist->At(ihist)->GetName(),hlist->At(ihist)->GetName())) areListsDifferent = kTRUE;
	}
      } else {
	areListsDifferent=kTRUE;
      }
    }

    // last check: if something is not ok die loudly 
    if (hlist->GetSize() != fOutputHist->GetSize()) {
      AliFatal("Mismatching size!");
    }
    Int_t nhist =  fOutputHist->GetSize();
    for(Int_t ihist = 0; ihist < nhist; ihist++){
      if(strcmp(fOutputHist->At(ihist)->GetName(),hlist->At(ihist)->GetName())){
	AliFatal(Form("Mismatching histos: %s -> %s", fOutputHist->At(ihist)->GetName(),hlist->At(ihist)->GetName()));
      }
    }

    if (foundDiffinThisIterStep){
      iter->Reset(); // We found a difference: previous lists could
		     // also be affected... We start from scratch
      Int_t n = 0;
      collections[n++].Clear();
      count = 0;
    }
    else {
//       AliInfo("hlist");
//       hlist->Print();
//       AliInfo("fOutputHist");
//       fOutputHist->Print();
      
      Int_t n = 0;
      collections[n++].Add(hlist);
      
      count++;
    }
  }

  Int_t n = 0;
  fOutputHist->Merge(&collections[n++]);
  
  delete iter;

  return count+1;
}


