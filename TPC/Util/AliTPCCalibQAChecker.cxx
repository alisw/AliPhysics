/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class to process a tree and create alarms based on thresholds            //
//  origin: jens wiechula: jens.wiechula@cern.ch                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <TObjArray.h>
#include <TString.h>
#include <TObjString.h>
#include <TStyle.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TLine.h>
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TIterator.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStopwatch.h>

#include <AliLog.h>

#include "AliTPCCalibQAChecker.h"

using namespace std;

AliTPCCalibQAChecker::AliTPCCalibQAChecker() :
  TNamed("AliTPCCalibQAChecker","AliTPCCalibQAChecker"),
  fTreePtr(0x0),
  fHistPtr(0x0),
  fGraphPtr(0x0),
  fNumberPtr(0x0),
  fHist(0x0),
  fIterSubCheckers(0x0),
  fArrSubCheckers(0x0),
  fArrAlarmDescriptions(0x0),
  fStrDrawRep(""),
  fStrDrawRepOpt(""),
  fStrDraw(""),
  fStrDrawOpt(""),
  fStrCuts(""),
  fAlarmType(kMean),
  fQualityLevel(kINFO),
  fHistRep(0x0)
{
  //
  // Default ctor
  //
  ResetAlarmThresholds();
}
//_________________________________________________________________________
AliTPCCalibQAChecker::AliTPCCalibQAChecker(const char* name, const char *title) :
  TNamed(name,title),
  fTreePtr(0x0),
  fHistPtr(0x0),
  fGraphPtr(0x0),
  fNumberPtr(0x0),
  fHist(0x0),
  fIterSubCheckers(0x0),
  fArrSubCheckers(0x0),
  fArrAlarmDescriptions(0x0),
  fStrDrawRep(""),
  fStrDrawRepOpt(""),
  fStrDraw(""),
  fStrDrawOpt(""),
  fStrCuts(""),
  fAlarmType(kMean),
  fQualityLevel(kINFO),
  fHistRep(0x0)
{
  //
  // TNamed ctor
  //
  ResetAlarmThresholds();
}
//_________________________________________________________________________
AliTPCCalibQAChecker::~AliTPCCalibQAChecker()
{
  //
  // Default ctor
  //
  if (fHistRep) delete fHistRep;
  if (fIterSubCheckers) delete fIterSubCheckers;
  if (fArrAlarmDescriptions) delete fArrAlarmDescriptions;
  if (fArrSubCheckers) delete fArrSubCheckers;
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::AddSubChecker(AliTPCCalibQAChecker *alarm)
{
  //
  // add a sub checker to this checker
  //
  if (!alarm) return;
  if (!fArrSubCheckers) {
    fArrSubCheckers=new TObjArray;
    fArrSubCheckers->SetOwner();
  }
  fArrSubCheckers->Add(alarm);
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::Process()
{
  //
  // Process the alarm thresholds, decide the alarm level, create the representation histogram
  //

  //reset quality level
  fQualityLevel=kINFO;

  TStopwatch s;
  s.Start();
  //decide which type of checker to use
  if (fArrSubCheckers && fArrSubCheckers->GetEntries()>0) ProcessSub();
  else if (fTreePtr && *fTreePtr) ProcessTree();
  else if (fHistPtr && *fHistPtr) ProcessHist();
  else if (fGraphPtr && *fGraphPtr) ProcessGraph();
  else if (fNumberPtr) ProcessNumber();
  s.Stop();
  AliInfo(Form("Processing Time (%s): %fs",GetName(),s.RealTime()));
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessSub()
{
  //
  // sub checker type checker
  //
  QualityFlag_t quality=kINFO;
  if (fArrSubCheckers && fArrSubCheckers->GetEntries()>0){
    TIter next(fArrSubCheckers);
    TObject *o=0x0;
    while ( (o=next()) ) {
      AliTPCCalibQAChecker *al=(AliTPCCalibQAChecker*)o;
      al->Process();
      QualityFlag_t subQuality=al->GetQuality();
      if (subQuality>quality) quality=subQuality;
    }
  }
  fQualityLevel=quality;
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessTree()
{
  //
  // process tree type checker
  //

  //Create Representation Histogram
  CreateRepresentationHist();
  //
//   if (!fTree) return;
  //chek for the quality

  switch (fAlarmType){
  case kNentries:
    ProcessEntries();
    break;
  case kMean:
  case kBinAny:
  case kBinAll:
    CreateAlarmHist();
    ProcessHist();
    ResetAlarmHist();
    break;
  }
  
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessHist()
{
  //
  // process histogram type checker
  //
  
  if (!(fHistPtr && *fHistPtr)) return;
    
  switch (fAlarmType){
  case kNentries:
  case kMean:
    ProcessMean();
    break;
  case kBinAny:
  case kBinAll:
    ProcessBin();
    break;
  }

  //create representation histogram if we are not in tree mode
  if (!fTreePtr){
    if (!fHistRep) {
      fHistRep=(*fHistPtr)->Clone();
      TH1* h=(TH1*)fHistRep;
      h->SetDirectory(0);
      h->SetNameTitle(Form("h%s",GetName()),GetTitle());
    }
    TH1* h=(TH1*)fHistRep;
    h->Reset();
    h->Add(*fHistPtr);
    if (!fHistRep->InheritsFrom(TH2::Class())) AddQualityLines(h);
  }
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessGraph()
{
  //
  // process graph type checker
  //
  if (!(fGraphPtr && *fGraphPtr)) return;
  Int_t npoints=(*fGraphPtr)->GetN();
  fQualityLevel=GetQuality(npoints,(*fGraphPtr)->GetY());
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessNumber()
{
  //
  // process number type checker
  //
  if (!fNumberPtr) return;
  fQualityLevel=GetQuality(*fNumberPtr);
  //create representation histogram
  if (!fHistRep){
    TH1* h=new TH1F(Form("h%s",GetName()),GetTitle(),1,0,1);
    h->GetXaxis()->SetBinLabel(1,"Value");
    AddQualityLines(h);
    h->SetDirectory(0);
    fHistRep=h;
  }
  TH1 *h=(TH1*)fHistRep;
  TMarker *marker=(TMarker*)h->GetListOfFunctions()->FindObject("TMarker");
  if (!marker) {
    marker=new TMarker(.5,0,20);
    h->GetListOfFunctions()->Add(marker);
  }
  marker->SetMarkerColor(GetQualityColor());
  marker->SetY(*fNumberPtr);
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessEntries()
{
  //
  // Processing function which analyses the number of affected rows of a tree draw
  //
  TString draw=fStrDraw;
  if (draw.IsNull()) return;
  
  TString cuts=fStrCuts;
  
  TString opt=fStrDrawOpt;
  opt+="goff";
  
  //draw and get the histogram
  Int_t res=(*fTreePtr)->Draw(draw.Data(),cuts.Data(),opt.Data());
  fQualityLevel=GetQuality(res);
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessMean()
{
  //
  // Processing function which analyses the mean of the resulting histogram
  //

  TH1* h=(*fHistPtr);
  Double_t value=h->GetMean();
  if (fAlarmType==kNentries) value=h->GetEntries();
  fQualityLevel=GetQuality(value);
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ProcessBin()
{
  //
  // Process a histogram bin by bin and check for thresholds
  //

  //bin quality counters
  Int_t nquality[kNQualityFlags];
  for (Int_t iquality=(Int_t)kINFO; iquality<kNQualityFlags; ++iquality) nquality[iquality]=0;
    
  TH1 *h=(*fHistPtr);
  
  Int_t nbinsX=h->GetNbinsX();
  Int_t nbinsY=h->GetNbinsY();
  Int_t nbinsZ=h->GetNbinsZ();
  Int_t nbinsTotal=nbinsX*nbinsY*nbinsZ;

  //loop over all bins
  for (Int_t ibinZ=1;ibinZ<nbinsZ+1;++ibinZ){
    for (Int_t ibinY=1;ibinY<nbinsY+1;++ibinY){
      for (Int_t ibinX=1;ibinX<nbinsX+1;++ibinX){
        Double_t value = (*fHistPtr)->GetBinContent(ibinX, ibinY, ibinZ);
        QualityFlag_t quality=GetQuality(value);
        nquality[quality]++;
      }
    }
  }

  //loop over Quality levels and set quality
  for (Int_t iquality=(Int_t)kINFO; iquality<kNQualityFlags; ++iquality){
    if (fAlarmType==kBinAny){
      if (nquality[iquality]) fQualityLevel=(QualityFlag_t)iquality;
    } else if (fAlarmType==kBinAll){
      if (nquality[iquality]==nbinsTotal) fQualityLevel=(QualityFlag_t)iquality;
    }
  }
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::CreateRepresentationHist()
{
  //
  // Create the representation histogram which will be shown in the draw function
  //
  ResetRepresentationHist();

  TString draw=fStrDrawRep;
  if (draw.IsNull()) {
    draw=fStrDraw;
    fStrDrawRepOpt=fStrDrawOpt;
  } else {
    draw.ReplaceAll("%alarm%",fStrDraw.Data());
  }
  if (draw.IsNull()) return;
  
  TString cuts=fStrCuts;
  
  TString opt=fStrDrawRepOpt;
  opt+="goff";
  
  Int_t res=(*fTreePtr)->Draw(draw.Data(),cuts.Data(),opt.Data());
  TH1 *hist=(*fTreePtr)->GetHistogram();
  if (res<0 || !hist){
    AliError(Form("Could not create representation histogram of alarm '%s'",GetName()));
    return;
  }
  fHistRep=hist->Clone();
  TH1 *h=(TH1*)fHistRep;
  h->SetDirectory(0);
  h->SetNameTitle(Form("h%s",GetName()),GetTitle());
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::CreateAlarmHist()
{
  //
  // create alarm histogram from the tree
  //
    
  TString draw=fStrDraw;
  if (draw.IsNull()) return;
  
  TString cuts=fStrCuts;
  
  TString opt=fStrDrawOpt;
  opt+="goff";
  
  //draw and get the histogram
  Int_t res=(*fTreePtr)->Draw(draw.Data(),cuts.Data(),opt.Data());
  fHist=(*fTreePtr)->GetHistogram();
  if (res<0 || !fHist){
    AliError(Form("Could not create alarm histogram of alarm '%s'",GetName()));
    return;
  }
  fHist->SetDirectory(0);
  fHistPtr=&fHist;
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ResetAlarmHist()
{
  //
  // delete the alarm histogram and reset the pointer
  //
  if (fHistPtr){
    if (*fHistPtr) delete *fHistPtr;
    fHistPtr=0x0;
  }
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::Draw(Option_t *option)
{
  //
  // object draw function
  // by default the pad backgound color is set to the quality level color
  // use 'nobc' to change this
  //

  //case of a node with subcheckers and no specific representation histogram
  if (HasSubCheckers()&&!fHistRep) {DrawSubNodes(option); return;}
  if (fHistRep) {DrawRepresentationHist(option); return;}
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::Print(Option_t * const option) const
{
  //
  // print the quality status. If we have sub checkers print recursively
  //
  TString sOpt(option);
  cout << sOpt << GetName() << ": " << GetQualityName() << endl;
  if (fArrSubCheckers && fArrSubCheckers->GetEntries()>0){
    sOpt.ReplaceAll("+-","  ");
    sOpt+="+-";
    TIter next(fArrSubCheckers);
    TObject *o=0x0;
    while ( (o=next()) ) o->Print(sOpt.Data());
  }
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::SetAlarmThreshold(const Double_t min, const Double_t max, const QualityFlag_t quality)
{
  //
  //set the alarm thresholds for a specific quality level
  //
  if ((Int_t)quality<(Int_t)kINFO||(Int_t)quality>=kNQualityFlags) return;
  fThresMin[quality]=min;
  fThresMax[quality]=max;
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ResetAlarmThreshold(const QualityFlag_t quality)
{
  //
  //set the alarm thresholds for a specific quality level
  //
  if ((Int_t)quality<(Int_t)kINFO||(Int_t)quality>=kNQualityFlags) return;
  fThresMin[quality]=0;
  fThresMax[quality]=0;
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::ResetAlarmThresholds()
{
  //
  //reset all the alarm thresholds
  //
  for (Int_t i=0;i<kNQualityFlags;++i){
    fThresMin[i]=0;
    fThresMax[i]=0;
  }
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::SetQualityDescription(const char* text, const QualityFlag_t quality)
{
  //
  // set an description for the quality level
  // %min and %max will be replaced by the min and max values of the alarm, when the quality
  // description is queried (see QualityDescription)
  //

  if (quality<kINFO||quality>kFATAL) return;
  if (! fArrAlarmDescriptions ) fArrAlarmDescriptions=new TObjArray(kNQualityFlags);
  TObjString *s=(TObjString*)fArrAlarmDescriptions->At(quality);
  if (!s) fArrAlarmDescriptions->AddAt(s=new TObjString,quality);
  s->SetString(text);
}

//_________________________________________________________________________
const AliTPCCalibQAChecker* AliTPCCalibQAChecker::GetSubChecker(const char* name, Bool_t recursive) const
{
  //
  // Return subnode with 'name'
  //
  TString sname(name);
  if (sname==GetName()) return this;
  if (!fArrSubCheckers || !fArrSubCheckers->GetEntries()) return 0x0;
  const AliTPCCalibQAChecker *al=0x0;
  if (recursive){
    TIter next(fArrSubCheckers);
    TObject *o=0x0;
    while ( (o=next()) ){
      AliTPCCalibQAChecker *sal=(AliTPCCalibQAChecker*)o;
      al=sal->GetSubChecker(name);
      if (al) break;
    }
  }else{
    al=dynamic_cast<AliTPCCalibQAChecker*>(fArrSubCheckers->FindObject(name));
  }
  return al;
}
//_________________________________________________________________________
Int_t AliTPCCalibQAChecker::GetNumberOfSubCheckers(Bool_t recursive) const
{
  //
  // get the number of sub checkers
  // if recursive get total number of non subchecker type sub checkers
  //
  Int_t nsub=0;
  if (recursive){
    if (!fArrSubCheckers) return 1;
    if (!fArrSubCheckers->GetEntries()) return 0;
    TIter next(fArrSubCheckers);
    TObject *o=0x0;
    while ( (o=next()) ){
      AliTPCCalibQAChecker *al=(AliTPCCalibQAChecker*)o;
      nsub+=al->GetNumberOfSubCheckers();
    }
  } else {
    if (fArrSubCheckers) nsub=fArrSubCheckers->GetEntries();
  }
  return nsub;
}
//_________________________________________________________________________
AliTPCCalibQAChecker* AliTPCCalibQAChecker::NextSubChecker()
{
  //
  // loop over sub checkers
  // if recursive, recursively return the pointers of non subchecker type sub checkers
  //
  if (!fArrSubCheckers || !fArrSubCheckers->GetEntries()) return 0;
  if (!fIterSubCheckers) fIterSubCheckers=fArrSubCheckers->MakeIterator();
  AliTPCCalibQAChecker *al=(AliTPCCalibQAChecker*)fIterSubCheckers->Next();
  if (!al){
    delete fIterSubCheckers;
    fIterSubCheckers=0x0;
  }
//   if (recursive && al->GetNumberOfSubCheckers(kFALSE)) al=al->NextSubChecker();
  return al;
}
//_________________________________________________________________________
const char* AliTPCCalibQAChecker::QualityName(const AliTPCCalibQAChecker::QualityFlag_t quality)
{
  //
  // get quality name for quality
  //
  switch (quality){
  case kINFO:
    return "Info";
    break;
  case kWARNING:
    return "Warning";
    break;
  case kERROR:
    return "Error";
    break;
  case kFATAL:
    return "Fatal";
    break;
  default:
    return "";
  }
}
//_________________________________________________________________________
Color_t AliTPCCalibQAChecker::QualityColor(const AliTPCCalibQAChecker::QualityFlag_t quality)
{
  //
  // get quality color for quality
  //
  Color_t info = kSpring-8;
  Color_t warning = kOrange;
  Color_t error = kRed;
  Color_t fatal = kRed+2;
  Color_t none = kWhite;
  
  switch(quality) {
  case kINFO :
    return info;
    break;
  case kWARNING :
    return warning;
    break;
  case kERROR :
    return error;
    break;
  case kFATAL :
    return fatal;
    break;
  default:
    return none;
  }
  return none;
  
}
//_________________________________________________________________________
const char* AliTPCCalibQAChecker::QualityDescription(const QualityFlag_t quality) const
{
  //
  // return description for quality
  //
  if (!fArrAlarmDescriptions || !fArrAlarmDescriptions->At(quality)) return "";
  TString s(fArrAlarmDescriptions->At(quality)->GetName());
  TString min, max;
  min+=fThresMin[quality];
  max+=fThresMax[quality];
  s.ReplaceAll("%min",min);
  s.ReplaceAll("%max",max);
  return s.Data();
}
//_________________________________________________________________________
Int_t AliTPCCalibQAChecker::DrawInPad(TVirtualPad *pad, Int_t sub)
{
  //
  // Draw recursively in 'pad'
  //
  
  if (fArrSubCheckers){
    if (fArrSubCheckers->GetEntries()>0){
      TIter next(fArrSubCheckers);
      TObject *o=0x0;
      while ( (o=next()) ) {
        AliTPCCalibQAChecker *al=(AliTPCCalibQAChecker*)o;
        sub=al->DrawInPad(pad,sub);
      }
    }
  } else {
    pad->cd(sub);
    ++sub;
    Draw();
  }
  return sub;
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::DrawSubNodes(Option_t */*option*/)
{
  //
  // divide the current pad in sub pads and draw all sub nodes
  //
  if (!gPad) new TCanvas;
  TPad *mother=(TPad*)gPad;
  mother->Clear();
  //find number of sub pads needed
  Int_t nPads = GetNumberOfSubCheckers();
  Int_t nCols = (Int_t)TMath::Ceil( TMath::Sqrt(nPads) );
  Int_t nRows = (Int_t)TMath::Ceil( (Double_t)nPads/(Double_t)nCols );
  gPad->Divide(nCols,nRows);
  
  DrawInPad(gPad);
  mother->Update();
  TPad *pad=0;
  Int_t ipad=1;
  while ( (pad=(TPad*)mother->GetPad(ipad)) ){
    TFrame* frame=(TFrame*)pad->GetPrimitive("TFrame");
    if (frame) frame->SetFillColor(kWhite);
    else {printf("no frame\n");}
    pad->Modified();
    ++ipad;
  }
  mother->Update();
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::DrawRepresentationHist(const Option_t *option)
{
  //
  // Draw the representation histogram
  //
  if (!fHistRep) return;
  if (!gPad) new TCanvas;
  Bool_t withBackColor=kTRUE;
  
  TString opt=option;
  opt.ToLower();
  
  if (opt.Contains("nobc")) withBackColor=kFALSE;
  opt.ReplaceAll("nobc","");
  
  if (opt.IsNull()) opt=fStrDrawRepOpt;
  opt.ToLower();
  
  opt.ReplaceAll("prof","");
  
  fHistRep->Draw(opt.Data());
  if (withBackColor) {
    gPad->SetFillColor(GetQualityColor());
  }
  
  gPad->Modified();
}
//_________________________________________________________________________
void AliTPCCalibQAChecker::AddQualityLines(TH1 *hist)
{
  //
  // add lines indicating the quality level to hist
  //

  Double_t xmin=hist->GetXaxis()->GetXmin();
  Double_t xmax=hist->GetXaxis()->GetXmax();
  Double_t min=0;
  Double_t max=0;
  
  for (Int_t i=(Int_t)kINFO; i<kNQualityFlags; ++i){
    if (fThresMin[i]>=fThresMax[i]) continue;
    min=fThresMin[i];
    max=fThresMax[i];
    TLine *lmin=new TLine(xmin,min,xmax,min);
    lmin->SetLineWidth(2);
    lmin->SetLineColor(QualityColor((QualityFlag_t)i));
    TLine *lmax=new TLine(xmin,max,xmax,max);
    lmax->SetLineWidth(2);
    lmax->SetLineColor(QualityColor((QualityFlag_t)i));
    hist->GetListOfFunctions()->Add(lmin);
    hist->GetListOfFunctions()->Add(lmax);
  }
  hist->SetAxisRange(min,max,"Y");
}
