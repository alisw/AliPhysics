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
#include "TGraphErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"

#include "AliLog.h"
#include "AliTRDtrendingManager.h"

ClassImp(AliTRDtrendingManager)

AliTRDtrendingManager* AliTRDtrendingManager::fgInstance=NULL;
//Bool_t AliTRDtrendingManager::fgTerminated = kFALSE;

//____________________________________________
AliTRDtrendingManager* AliTRDtrendingManager::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
//   if (fgTerminated) return NULL;
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
  
/*  fgTerminated = kTRUE;*/
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
{
// Constructor
}

//____________________________________________
AliTRDtrendingManager::~AliTRDtrendingManager()
{
// Destructor
  if(fEntries){
    fEntries->Delete();
    delete fEntries;
  }
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
    MakeList(1000);
  }

  if(GetValue(name)){
    AliInfo(Form("Trend value \"%s\" already in list. Use Modify function", name));
    return;
  }
  // create new trending value`
  AliTRDtrendValue *fValue = new AliTRDtrendValue(name, title?title:"");
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
    r->Delete(); delete r;
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
      r->Delete(); delete r;
    }
    n->Delete(); delete n;
  }
  fEntries->AddLast(fValue);
}

//____________________________________________
AliTRDtrendValue* AliTRDtrendingManager::GetValue(const Char_t *name)
{
// Search trend value list by value "name" formatted according to "class_name"
  if(!fEntries){
    AliError("No trending map defined.");
    return NULL;
  }
  return (AliTRDtrendValue*)fEntries->FindObject(name);
}

//____________________________________________
TH1* AliTRDtrendingManager::MakeTrends(const char *fileList, TObjArray *dump)
{
// Make trends with reference to DB for all trend files in "fileList".
// The DB should be loaded
  if(!fEntries){
    AliWarning("Trending map undefined");
    return NULL;
  }
  Int_t ntv(fEntries->GetEntries());

  FILE *fp(NULL);
  if(!(fp= fopen(fileList, "rt"))){
    AliWarning(Form("Can not open file list \"%s\"", fileList));
    return NULL;
  }
  Int_t *na = new Int_t[ntv]; memset(na, 0, ntv*sizeof(Int_t));
  Float_t *la = new Float_t[ntv]; memset(la, 0, ntv*sizeof(Float_t));
  Float_t *lm = new Float_t[ntv]; for(Int_t im(0); im<ntv; im++) lm[im] = 1.e5;
  Float_t *lM = new Float_t[ntv]; for(Int_t im(0); im<ntv; im++) lM[im] = -1.e5;
  TGraphErrors **g = new TGraphErrors*[ntv]; memset(g, 0, ntv*sizeof(TGraphErrors*));
  AliTRDtrendValue *TV(NULL), *tv(NULL);
  TString sfp; Int_t run[10000], nr(0);
  while(sfp.Gets(fp)){
    // guess run no from path. Do we need a more robust approach ?
    TObjArray *afp=sfp.Tokenize("/");
    Int_t idx = afp->GetEntries()-2;
    Int_t rno = ((TObjString*)(*afp)[idx])->GetString().Atoi();
    afp->Delete(); delete afp;
    if(!TFile::Open(sfp.Data())) continue;

    run[nr] = rno; Int_t nmiss(0);
    for(Int_t it(0); it<ntv; it++){
      if(!(TV = (AliTRDtrendValue*)fEntries->At(it))) continue;
      if(!(tv = (AliTRDtrendValue*)gFile->Get(TV->GetName()))) {
        AliDebug(1, Form("Missing %09d.%s", rno, TV->GetName()));
        nmiss++;
        continue;
      }
      if(tv->GetVal()<=-998. ||
         (strstr(TV->GetName(), "TRDcheckDET")&&TMath::Abs(tv->GetVal())<1.e-5) ||
         (strstr(TV->GetName(), "TRDefficiency")&&tv->GetVal()<1.e-5) ||
         (!(strcmp(TV->GetName(), "TRDcheckDET_ChargeTracklet"))&&TMath::Abs(tv->GetVal())<1.e1)) continue;
      if(IsRelativeMeanSigma()){
        (*tv)/=(*TV);
        la[it]+=tv->GetVal(); na[it]++;
      } else {
        if(tv->GetVal()<lm[it]) lm[it]=tv->GetVal();
        if(tv->GetVal()>lM[it])lM[it]=tv->GetVal();
      }
      if(!g[it]){
        g[it] = new TGraphErrors();
        g[it]->SetNameTitle(TV->GetName(), TV->GetTitle());
        g[it]->SetMarkerStyle(4);g[it]->SetMarkerSize(1.2);
        g[it]->SetLineStyle(2);g[it]->SetLineWidth(1);
      }
      Int_t ip(g[it]->GetN());
      g[it]->SetPoint(ip, nr, tv->GetVal());
      g[it]->SetPointError(ip, 0., tv->GetErr());
    }
    if(Float_t(nmiss)/ntv>.1) AliWarning(Form("Run[%09d] Missing %6.2f%% values", rno, 1.e2*nmiss/ntv));
    nr++;
  }

  // Draw
  TH1 *hT = new TH1F("hT", ";#bf{RUN};", nr, -0.5, nr-0.5);
  TAxis *ax = hT->GetXaxis(); ax->SetTitleOffset(2.6);ax->CenterTitle(); ax->SetBit(TAxis::kLabelsVert);
  TAxis *ay = hT->GetYaxis(); ay->SetTitleOffset(IsRelativeMeanSigma()?0.4:0.75);ay->CenterTitle(); ay->SetAxisColor(kRed); ay->SetDecimals();
  for(Int_t ir(0); ir<nr; ir++) ax->SetBinLabel(ir+1, Form("%09d", run[ir]));

  TLine *line(NULL);
  TCanvas *c = new TCanvas("c", "TRD Trend", 1, 1, 2400, 1000);
  c->SetLeftMargin(IsRelativeMeanSigma()?0.03666361:0.05685619);
  c->SetRightMargin(0.005499542);
  c->SetTopMargin(0.02542373);
  c->SetBottomMargin(0.1758475);
  for(Int_t it(0); it<ntv; it++){
    if(!g[it]) continue;
    c->Clear();
    if(IsRelativeMeanSigma()){
      ay->SetRangeUser(-5, 5);
      ay->SetTitle(Form("#bf{%s [#sigmau]}", g[it]->GetTitle()));
      line = new TLine(-0.5, na[it]?(la[it]/na[it]):0., nr-0.5, na[it]?(la[it]/na[it]):0.);
      line->SetLineColor(kBlue);
    } else {
      ay->SetRangeUser(lm[it]-0.1*(lM[it]-lm[it]), lM[it]+0.1*(lM[it]-lm[it]));
      ay->SetTitle(Form("#bf{%s}", g[it]->GetTitle()));
    }
    hT->Draw("p");
    g[it]->Draw("ple5");
    if(line) line->Draw();
    c->Modified(); c->Update(); c->SaveAs(Form("Trend_%s.gif", g[it]->GetName()));
    if(dump) dump->Add(g[it]);
    else delete g[it];
    if(line) delete line;
  }
  delete [] g;
  delete [] lm;
  delete [] lM;
  delete [] la;
  delete [] na;
  return hT;
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

  AliTRDtrendValue *fValue(NULL);
  if(!(fValue = GetValue(name))) { 
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
    r->Delete(); delete r;
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
      r->Delete(); delete r;
    }
    n->Delete(); delete n;
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
// Load TRD trending DB from $ALICE_PHYSICS/PWGPP/TRD/data.

  AliDebug(1, "Loading TRD trending ...");
  if(!TFile::Open(fn)) return;

  TList *tvList = gFile->GetListOfKeys(); TIterator *iter(tvList->MakeIterator()); AliTRDtrendValue *tv(NULL);
  MakeList(tvList->GetEntries());
  TKey *ktv(NULL);
  while((ktv = (TKey*)iter->Next())){
    tv = (AliTRDtrendValue*)gFile->Get(ktv->GetName());
    fEntries->AddLast(new AliTRDtrendValue(*tv));
    delete tv;
  }
  gFile->Close();
}

//____________________________________________
void AliTRDtrendingManager::MakeList(Int_t entries)
{
// allocate trending values array
  fEntries = new TObjArray(entries);
  fEntries->SetOwner();
  fEntries->SetName("TrendValues");
}
