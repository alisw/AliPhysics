#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include "KMCDetector.h"
#endif

// Process the output file of testDetKMC.C macro to prepare summary 
// for the "summary class" icl (definition of acceptable track).
// See ProcessSummary comment for details.
// The names of the graphs are prefixed by `pref`.
// If outF is provided, the graphs will be stored there


TObjArray* PrepSummaryKMC(const char* sumf, int cls, const char* pref=0, const char* outF=0);
TObjArray* ProcessSummary(TObjArray* sums, int icl, const char* pref=0);

enum {kSigAD,kSigAZ,kSigAPt,kSigD,kSigZ,kSigPt,kEff,kUpd};
TF1* gs = 0;

TObjArray* PrepSummaryKMC(const char* sumf, int cls, const char* pref, const char* outF)
{
  if (!gROOT->GetClass("KMCDetector")) gROOT->LoadMacro("KMCDetector.cxx+");
  TFile* fl = TFile::Open(sumf);
  if (!fl) {printf("No file %s\n",sumf); return 0;}
  TObjArray* arrs = (TObjArray*)gDirectory->Get("trSum");
  if (!arrs) {printf("No summary in file %s\n",sumf); return 0;}
  //
  TObjArray* sums =  ProcessSummary(arrs,cls,pref);
  if (!sums) return 0;
  //
  if (outF) {
    TFile* flOut = TFile::Open(outF,"update");
    if (!flOut) {printf("Failed to open output file %s\n",outF);}
    else {
      flOut->WriteObject(sums,sums->GetName(),"kSingleKey");
      flOut->Close();
      delete flOut;
      printf("Stored array %s in %s\n",sums->GetName(), outF);
    }
  }
  //
  return sums;
}


TObjArray* ProcessSummary(TObjArray* arrs, int icl, const char* pref)
{
  // Process TObjArray (e.g. for set of pt bins) of TObjArray of KMCTrackSummary objects:
  // pick the KMCTrackSummary for "summary class" icl (definition of acceptable track) and create
  // graphs vs bin.
  // These graphs are returned in new TObjArray
  //
  TString prefs = pref;
  if (!gs) gs = new TF1("gs","gaus",-1,1);
  //
  int nb = arrs->GetEntriesFast();
  TObjArray* sums = (TObjArray*) arrs->At(0);
  int nclass = sums->GetEntriesFast();
  if (icl>=nclass) {printf("summary set has %d classes only, %d requested\n",nclass,icl);return 0;}
  //
  KMCTrackSummary* sm = (KMCTrackSummary*)sums->At(icl);
  //
  TH1* h;
  //
  h = sm->GetHMCSigDCARPhi();  // MC resolution in transverse DCA 
  TGraphErrors * grSigD = 0;
  if (h) {
    grSigD = new TGraphErrors(nb);
    grSigD->SetName(Form("%s%s",prefs.Data(),h->GetName()));
    grSigD->SetTitle(Form("%s%s",prefs.Data(),h->GetTitle()));
  }
  //
  TGraphErrors * grSigZ = 0;
  h = sm->GetHMCSigDCAZ();    // MC resolution in Z DCA
  if (h) {
    grSigZ = new TGraphErrors(nb);
    grSigZ->SetName(Form("%s%s",prefs.Data(),h->GetName()));
    grSigZ->SetTitle(Form("%s%s",prefs.Data(),h->GetTitle()));
  }
  //
  TGraphErrors * grSigAD = 0; // anaitical estimate for resolution in transverse DCA 
  {
    grSigAD = new TGraphErrors(nb);
    grSigAD->SetName(Form("%s%s",prefs.Data(),"sigmaDan"));
    grSigAD->SetTitle(Form("%s%s",prefs.Data(),"#sigmaD an"));
  }
  //
  TGraphErrors * grSigAZ = 0; // anaitical estimate for resolution in Z DCA 
  {
    grSigAZ = new TGraphErrors(nb);
    grSigAZ->SetName(Form("%s%s",prefs.Data(),"sigmaZan"));
    grSigAZ->SetTitle(Form("%s%s",prefs.Data(),"#sigmaZ an"));
  }
  //
  //
  TGraphErrors * grSigPt = 0; // MC res. in pt
  {
    grSigPt = new TGraphErrors(nb);
    grSigPt->SetName(Form("%s%s",prefs.Data(),"sigmaPt"));
    grSigPt->SetTitle(Form("%s%s",prefs.Data(),"#sigmaPt"));
  }
  //
  TGraphErrors * grSigAPt = 0; // analitycal res. in pt
  {
    grSigAPt = new TGraphErrors(nb);
    grSigAPt->SetName(Form("%s%s",prefs.Data(),"sigmaPtan"));
    grSigAPt->SetTitle(Form("%s%s",prefs.Data(),"#sigmaPt an"));
  }

  //
  TGraphErrors * grEff = 0; // MC efficiency
  {
    grEff = new TGraphErrors(nb);
    grEff->SetName(Form("%s_rate",prefs.Data()));
    grEff->SetTitle(Form("%s Rate",prefs.Data()));
  }
  //
  TGraphErrors * grUpd = 0; // number of Kalman track updates
  {
    grUpd = new TGraphErrors(nb);
    grUpd->SetName(Form("%s_updCalls",prefs.Data()));
    grUpd->SetTitle(Form("%s Updates",prefs.Data()));
  }
  //
  for (int ib=0;ib<nb;ib++) {
    sums = (TObjArray*) arrs->At(ib);
    sm = (KMCTrackSummary*)sums->At(icl);
    KMCProbe& prbRef = sm->GetRefProbe();
    KMCProbe& prbAn  = sm->GetAnProbe();
  
    double pt = prbRef.Pt();
    //
    if (grSigAD) {
      grSigAD->SetPoint(ib, pt,prbAn.GetSigmaY2()>0 ? TMath::Sqrt(prbAn.GetSigmaY2()) : 0.);
    }
    //
    if (grSigAZ) {
      grSigAZ->SetPoint(ib, pt,prbAn.GetSigmaZ2()>0 ? TMath::Sqrt(prbAn.GetSigmaZ2()) : 0.);
    }
    //
    if (grSigAPt) {
      double pts = TMath::Sqrt(prbAn.GetSigma1Pt2());
      grSigAPt->SetPoint(ib, pt,pts>0 ? pts*pt : 0.);
    }
    //
    if (grSigPt) {
      h = sm->GetHMCSigPt();
      h->Fit(gs,"0q");
      grSigPt->SetPoint(ib, pt, gs->GetParameter(2));
      grSigPt->SetPointError(ib, 0, gs->GetParError(2));
    }
    //
     if (grSigD) {
      h = sm->GetHMCSigDCARPhi();
      h->Fit(gs,"0q");
      grSigD->SetPoint(ib, pt,gs->GetParameter(2));
      grSigD->SetPointError(ib, 0,gs->GetParError(2));      
    }
    //
    if (grSigZ) {
      h = sm->GetHMCSigDCAZ();
      h->Fit(gs,"0q");
      grSigZ->SetPoint(ib, pt,gs->GetParameter(2));
      grSigZ->SetPointError(ib, 0,gs->GetParError(2));      
    }
    //
    if (grEff) {
      grEff->SetPoint(ib, pt,sm->GetEff());
      grEff->SetPointError(ib, 0,sm->GetEffErr());
    }
    //
    if (grUpd) {
      grUpd->SetPoint(ib, pt,sm->GetUpdCalls());
      grUpd->SetPointError(ib, 0, 0);
    }
  }
  //
  TObjArray* dest = new TObjArray();
  dest->AddAtAndExpand(grSigAD,kSigAD);
  dest->AddAtAndExpand(grSigAZ,kSigAZ);  
  dest->AddAtAndExpand(grSigAPt,kSigAPt);  
  dest->AddAtAndExpand(grSigD,kSigD);
  dest->AddAtAndExpand(grSigZ,kSigZ);  
  dest->AddAtAndExpand(grSigPt,kSigPt);  
  dest->AddAtAndExpand(grEff,kEff);
  dest->AddAtAndExpand(grUpd,kUpd);
  //
  if (!prefs.IsNull()) dest->SetName(pref);
  return dest;
}
