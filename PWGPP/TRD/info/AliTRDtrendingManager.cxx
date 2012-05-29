////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Trend Value Manager                                                   //
//                                                                        //
//  Mediates interaction with DB (OCDB ?!)                                //                                                                      //                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"

#include "AliLog.h"
#include "AliTRDtrendingManager.h"

ClassImp(AliTRDtrendingManager)

AliTRDtrendingManager* AliTRDtrendingManager::fgInstance=NULL;
Bool_t AliTRDtrendingManager::fgTerminated = kFALSE;

//____________________________________________
AliTRDtrendingManager* AliTRDtrendingManager::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgTerminated != kFALSE) return NULL;

  if (!fgInstance) fgInstance = new AliTRDtrendingManager();

  return fgInstance;
}

//____________________________________________
void AliTRDtrendingManager::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag,
  // instances cannot be requested anymore
  // This function can be called several times.
  //
  
  fgTerminated = kTRUE;
  
  if (fgInstance != NULL) {
    if(TFile::Open("TRD.Trend.root", "RECREATE")){
      if(fEntries) fEntries->Write();
      gFile->Close();
    }
    delete fgInstance; fgInstance = NULL;
  }
}

//____________________________________________
AliTRDtrendingManager::AliTRDtrendingManager() 
  : TObject()
  ,fEntries(NULL)
  ,fValue(NULL)
{
// Constructor
//  fRunRange[0] = 0; fRunRange[1] = AliCDBRunRange::Infinity();
}

//____________________________________________
AliTRDtrendingManager::~AliTRDtrendingManager()
{
// Destructor
  if(fValue) delete fValue;
  if(fEntries) delete fEntries;
}

//____________________________________________
void AliTRDtrendingManager::AddValue(
   const Char_t *name
  ,Double_t mean,Double_t sigm
  ,const Char_t *title
  ,const Char_t *responsible
  ,const Char_t *notifiables
  ,Char_t **messages
  )
{
// Expert Function !!!
// Add a trend value to the map already loaded
// If no map loaded create a new one from scratch
//
// class_name : name of the performance task 
// name       : name of the value to be trended
// title      : description of the value to be trended
// messages   : array of alarm messages for each alarm level
// responsible: name and email of the responsible person. Format "name/email"
// notifiables: name and email of the notifiable persons. Format "name1/email1, name2/email2, etc"
//

  if(!fEntries){ // if no trending map defined create one
    AliDebug(1, "No trending map loaded. Create one from scratch.");
    fEntries = new TObjArray(50);
    fEntries->SetOwner();
    fEntries->SetName("values");
  }

  if(!(fValue = GetValue(name))){
    // create new trending value`
    fValue = new AliTRDtrendValue(name, title?title:"");
    fValue->Set(mean, sigm);
    if(messages) for(Int_t ilevel(AliTRDtrendValue::kNlevels); ilevel--;) if(messages[ilevel]) fValue->SetAlarm(ilevel, messages[ilevel]);
    TObjArray *r(NULL);
    if(responsible){
      TString s(responsible);
      r=s.Tokenize("/");
      if(r->GetEntriesFast()!=2){
        AliWarning("Responsible name/email incorrectly formated.");
      } else {
        fValue->SetResponsible(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
      }
    }
    if(notifiables){
      TString s(notifiables);
      TObjArray *n=s.Tokenize(",");
      for(Int_t in(0); in<TMath::Min(AliTRDtrendValue::kNnotifiable, n->GetEntriesFast()); in++){
        TString ss(((TObjString*)n->At(in))->String());
        r=ss.Tokenize("/");
        if(r->GetEntriesFast()!=2){
          AliWarning(Form("Notifiable person name/email incorrectly formated for [%s].", ss.Data()));
        } else {
          fValue->SetNotifiable(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
        }
      }
    }
  }

  fEntries->AddLast(new AliTRDtrendValue(*fValue));
}

//____________________________________________
AliTRDtrendValue* AliTRDtrendingManager::GetValue(const Char_t *name)
{
// Search trend value list by value "name" formatted according to "class_name"
  if(!fEntries){
    AliError("No trending map defined.");
    return NULL;
  }
  fValue = (AliTRDtrendValue*)fEntries->FindObject(name);
  return fValue;
}

//____________________________________________
Bool_t AliTRDtrendingManager::MakeTrends(const char *fileList)
{
// Make trends with reference to DB for all trend files in "fileList".
// The DB should be loaded
  if(!fEntries){
    AliWarning("Trending map undefined");
    return kFALSE;
  }
  Int_t ntv(fEntries->GetEntries());

  FILE *fp(NULL);
  if(!(fp= fopen(fileList, "rt"))){
    AliWarning(Form("Can not open file list \"%s\"", fileList));
    return kFALSE;
  }

  TGraph **g = new TGraph*[ntv]; memset(g, 0, ntv*sizeof(TGraph*));
  AliTRDtrendValue *TV(NULL), *tv(NULL);
  TString sfp; Int_t run[10000], nr(0);
  while(sfp.Gets(fp)){
    // guess run no from path. Do we need a more robust approach ?
    TObjArray *afp=sfp.Tokenize("/");
    Int_t idx = afp->GetEntries()-2;
    Int_t rno = ((TObjString*)(*afp)[idx])->GetString().Atoi();

    if(!TFile::Open(sfp.Data())) continue;
    run[nr] = rno;
    for(Int_t it(0); it<ntv; it++){
      if(!(TV = (AliTRDtrendValue*)fEntries->At(it))) continue;
      if(!(tv = (AliTRDtrendValue*)gFile->Get(TV->GetName()))) {
        AliWarning(Form("Missing %09d.%s", rno, TV->GetName()));
        continue;
      }
      (*tv)/=(*TV);
      if(!g[it]){
        g[it] = new TGraph();
        g[it]->SetNameTitle(TV->GetName(), TV->GetTitle());
        g[it]->SetMarkerStyle(4);g[it]->SetMarkerSize(0.8);
      }
      g[it]->SetPoint(g[it]->GetN(), nr, tv->GetVal());
    }
    nr++;
  }

  // Draw
  TH1 *hT = new TH1F("hT", ";#bf{RUN};", nr, -0.5, nr-0.5);
  TAxis *ax = hT->GetXaxis(); ax->SetTitleOffset(2.6);ax->CenterTitle(); ax->SetBit(TAxis::kLabelsVert);
  TAxis *ay = hT->GetYaxis(); ay->SetTitleOffset(0.4);ay->CenterTitle(); ay->SetAxisColor(kRed); ay->SetRangeUser(-5, 5);
  for(Int_t ir(0); ir<nr; ir++) ax->SetBinLabel(ir+1, Form("%09d", run[ir]));

  TCanvas *c = new TCanvas("c", "TRD Trend", 1, 1, 1200, 500);
  c->SetLeftMargin(0.03666361);
  c->SetRightMargin(0.005499542);
  c->SetTopMargin(0.02542373);
  c->SetBottomMargin(0.1758475);
  for(Int_t it(0); it<ntv; it++){
    c->Clear();
    ay->SetTitle(Form("#bf{%s [#sigmau]}", g[it]->GetTitle())); hT->Draw("p");
    g[it]->Draw("p");
    c->Modified(); c->Update(); c->SaveAs(Form("Trend_%s.gif", g[it]->GetName()));
    delete g[it];
  }
  delete [] g;
  return kTRUE;
}

//____________________________________________
Bool_t AliTRDtrendingManager::ModifyValue(
  const Char_t *name
  ,const Char_t *title
  ,Double_t mean, Double_t sgm
  ,Char_t **messages
  ,const Char_t *responsible
  ,const Char_t *notifiables
  )
{
// Expert Function !!!
// Modify a trend value in the map already loaded
// see function AddValue() for explanation of input format. 

  if(!fEntries){
    AliError("No trending map loaded.");
    return kFALSE;
  }
  AliWarning("*** EXPERT FUNCTION *** This function is modifying one trending value to the current DB. Continue if you know what yout do!");

  if(!GetValue(name)){
    AliError(Form("Missing trending value %s", name));
    return kFALSE;
  }  
  
  fValue->SetTitle(title);
  fValue->Set(mean, sgm);
  if(messages){ 
    for(Int_t ilevel(AliTRDtrendValue::kNlevels); ilevel--;) fValue->SetAlarm(ilevel, messages[ilevel]);
  }
  TString s;
  TObjArray *r(NULL);
  if(responsible){ 
    s=responsible;
    r=s.Tokenize("/");
    if(r->GetEntriesFast()!=2){ 
      AliWarning("Responsible name/email incorrectly formated.");
    } else { 
      fValue->SetResponsible(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
    }
  }
  if(notifiables){
    s=notifiables;
    TObjArray *n=s.Tokenize(",");
    for(Int_t in(0); in<TMath::Min(AliTRDtrendValue::kNnotifiable, n->GetEntriesFast()); in++){
      TString ss(((TObjString*)n->At(in))->String());
      r=ss.Tokenize("/");
      if(r->GetEntriesFast()!=2){ 
        AliWarning(Form("Notifiable person name/email incorrectly formated for [%s].", ss.Data()));
      } else { 
        fValue->SetNotifiable(((TObjString*)r->At(0))->String().Data(), ((TObjString*)r->At(1))->String().Data());
      }
    }
  }
  return kTRUE;
}

//____________________________________________
void AliTRDtrendingManager::Print(Option_t *o) const
{
// Dump trend value list
  if(!fEntries){
    AliError("No trending map available.");
    return;
  }

  for(Int_t iv(0); iv<fEntries->GetEntriesFast(); iv++){
    ((AliTRDtrendValue*)fEntries->At(iv))->Print(o);
  }
}

//____________________________________________
void AliTRDtrendingManager::Load(const char *fn)
{
// Load TRD trending DB from $ALICE_ROOT/PWGPP/TRD/data.

  AliDebug(1, "Loading TRD trending ...");
  if(!TFile::Open(fn)) return;

  TList *tvList = gFile->GetListOfKeys(); TIterator *iter(tvList->MakeIterator()); AliTRDtrendValue *tv(NULL);
  fEntries = new TObjArray(tvList->GetEntries());
  fEntries->SetOwner();
  TKey *ktv(NULL);
  while((ktv = (TKey*)iter->Next())){
    tv = (AliTRDtrendValue*)gFile->Get(ktv->GetName());
    fEntries->AddLast(new AliTRDtrendValue(*tv));
  }
}
