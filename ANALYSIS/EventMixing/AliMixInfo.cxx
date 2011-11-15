//
// Class AliMixInfo
//
// AliMixInfo object contains information about one cut on for event mixing
// available for users containing mixing information
//
// authors:
//          Martin Vala (martin.vala@cern.ch)
//

#include <TList.h>
#include <TPad.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPavesText.h>
#include <TCanvas.h>

#include "AliLog.h"

#include "AliMixInfo.h"
#include "AliMixEventPool.h"
#include "AliMixEventCutObj.h"

ClassImp(AliMixInfo)

//_________________________________________________________________________________________________
AliMixInfo::AliMixInfo(const char *name, const char *title) :
   TNamed(name, title),
   fHistogramList(0)
{
   //
   // Default constructor.
   //
   AliDebug(AliLog::kDebug + 5, "<-");
   AliDebug(AliLog::kDebug + 5, "->");
}
//_________________________________________________________________________________________________
AliMixInfo::AliMixInfo(const AliMixInfo &obj) :
   TNamed(obj),
   fHistogramList(obj.fHistogramList)
{
   //
   // Copy constructor.
   //
   AliDebug(AliLog::kDebug + 5, "<-");
   AliDebug(AliLog::kDebug + 5, "->");
}
//_________________________________________________________________________________________________
AliMixInfo::~AliMixInfo()
{
   //
   // Destructor
   //
}

//_________________________________________________________________________________________________
void AliMixInfo::CreateHistogram(AliMixInfo::EInfoHistorgramType type, Int_t nbins, Int_t min, Int_t max)
{
   //
   // Create mix info histograms
   //
   if (!fHistogramList) {
      fHistogramList = new TList;
      fHistogramList->SetOwner(kTRUE);
   }
   TH1I *hist = (TH1I *) fHistogramList->FindObject(GetNameHistogramByType(type));
   if (hist) return;
   hist = new TH1I(GetNameHistogramByType(type), GetTitleHistogramByType(type), nbins, min, max);
   fHistogramList->Add(hist);
}

//_________________________________________________________________________________________________
void AliMixInfo::FillHistogram(AliMixInfo::EInfoHistorgramType type, Int_t value)
{
   //
   // Create mix info histograms
   //
   if (type == kMixedEvents && value < 0) return;
   if (!fHistogramList) {
      AliError("fHistogramList is null");
      return;
   }
   TH1I *hist = (TH1I *) fHistogramList->FindObject(GetNameHistogramByType(type));
   if (hist) {
      hist->Fill(value);
      AliDebug(AliLog::kDebug, Form("%s was filled with %d sum is %.0f", GetNameHistogramByType(type), value, hist->GetBinContent(value)));
   } else {
      AliError(Form("Problem filling histogram %s", GetNameHistogramByType(type)));
   }
}

//_________________________________________________________________________________________________
const char *AliMixInfo::GetNameHistogramByType(Int_t index) const
{
   //
   // Retruns name of cut
   //
   switch (index) {
      case kMainEvents:
         return "hMainEvents";
      case kMixedEvents:
         return "hMixedEvents";
   }
   return "";
}

//_________________________________________________________________________________________________
const char *AliMixInfo::GetTitleHistogramByType(Int_t index) const
{
   //
   // Retruns name of cut
   //
   switch (index) {
      case kMainEvents:
         return "Main Events";
      case kMixedEvents:
         return "Mixed Events";
   }
   return "";
}

//_________________________________________________________________________________________________
void AliMixInfo::Print(Option_t *option) const
{
   //
   // Print Mix info
   //
   if (!fHistogramList) return;
   if (option)
      AliInfo(Form("Name %s with option is %s", GetName(), option));
   TIter next(fHistogramList);
   TH1I *h = 0;
   for (Int_t i = 0; i < fHistogramList->GetEntries(); i++) {
      h = dynamic_cast<TH1I *>(fHistogramList->At(i));
      if (h) {
         h->Print();
         continue;
      }
   }
}
//_________________________________________________________________________________________________
void AliMixInfo::Draw(Option_t* option)
{
   //
   // Drwas mixi info canvas
   //
   if (!fHistogramList) return;

   // creating main canvas
   TCanvas *cMain = new TCanvas("cMain", "Mixing Info", 500, 500);
   if (!cMain) return;
   cMain->Divide(1, 2, 0.001, 0.001);
   cMain->cd(1);
//     TVirtualPad *upperPad = cMain->cd(1);
//     upperPad->Divide(2,1);
//     TVirtualPad *upperPad1 = gPad->cd(2);
//     upperPad1->Divide(1,2);
//     upperPad1->cd(1);
   TPavesText*text = new TPavesText(0.05, 0.05, 0.95, 0.95, 1);
   text->SetName("mixInfoText");
   text->AddText("Help:");
   text->AddText("Move over histogram to see mix info for different bins");
   text->Draw();

   // gets corresponding histograms
   TH1I *hMain = GetHistogramByType(kMainEvents);
   if (!hMain) {
      AliError("hMain is null");
      return;
   }
   TH1I *hMix =  GetHistogramByType(kMixedEvents);
   if (!hMix) {
      AliError("hMix is null");
      return;
   }


   TH2I *hMixInfo2D = 0;
//     TH1I *hOK=0,*hBad=0;
   AliMixEventPool *evPool = (AliMixEventPool *) GetEventPool("mixEventPool");
   if (evPool) {
      Int_t mixNum = evPool->GetMixNumber();
//         hOK = (TH1I *) hMain->Clone();
//         hBad = (TH1I *) hMain->Clone();
      if (!hMixInfo2D) hMixInfo2D = new TH2I("hMixInfo2D", "hMixInfo2D", hMain->GetXaxis()->GetNbins() + 1, hMain->GetXaxis()->GetXmin() - 1, hMain->GetXaxis()->GetXmax(), 1, 0, 1);
      for (Int_t iBin = 0; iBin < hMain->GetNbinsX() + 1; iBin++) {
         if (!iBin) {
            hMixInfo2D->SetBinContent(iBin + 1, 1, 1);
         } else if (!hMain->GetBinContent(iBin) && !hMix->GetBinContent(iBin)) {
            hMixInfo2D->SetBinContent(iBin + 1, 1, 2);
         } else if (hMix->GetBinContent(iBin) == mixNum * hMain->GetBinContent(iBin)) {
            hMixInfo2D->SetBinContent(iBin + 1, 1, 4);
         } else {
            hMixInfo2D->SetBinContent(iBin + 1, 1, 3);
         }
      }
   }

   TStyle *style = gStyle;
   Int_t cols[4] = { kYellow, kViolet, kRed, kGreen  };
   style->SetPalette(4, cols);
   cMain->cd(2);
//    cMain->SetGrid();
   if (hMixInfo2D) {
      hMixInfo2D->SetMaximum(4);
      hMixInfo2D->SetStats(0);
      hMixInfo2D->SetTitle("");
      hMixInfo2D->GetXaxis()->SetNdivisions(510);
      hMixInfo2D->GetYaxis()->SetNdivisions(0);
   }

   if (hMixInfo2D) hMixInfo2D->Draw(Form("COL %s", option));
//
//     TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
//     legend->AddEntry(hOK,"OK","f");
//     legend->AddEntry(hBad,"NOT OK","f");
//     legend->Draw();

   cMain->cd(2)->AddExec("dynamic", Form("AliMixInfo::DynamicExec((AliMixInfo*)0x%lx)", (ULong_t)this));
}

//_________________________________________________________________________________________________
void AliMixInfo::DynamicExec(AliMixInfo *const mixInfo)
{
   //
   // Function which is run when user move mouse over mix info
   //

   if (!mixInfo) return;

   TObject *select = gPad->GetSelected();
   if (!select) return;
   if (!select->InheritsFrom(TH2I::Class())) {
      gPad->SetUniqueID(0);
      return;
   }

   TH2I *hSelected = (TH2I *) select;
   gPad->GetCanvas()->FeedbackMode(kTRUE);

   //erase old position and draw a line at current position
   Int_t uid = gPad->GetUniqueID();
//     int pxold = gPad->GetUniqueID();
   Int_t px = gPad->GetEventX();
//     Int_t py = gPad->GetEventY();
//     float uxmin = gPad->GetUxmin();
//     float uxmax = gPad->GetUxmax();
//     float uymin = gPad->GetUymin();
//     float uymax = gPad->GetUymax();
   //     Int_t pxmin = gPad->XtoAbsPixel ( uxmin );
   //     Int_t pxmax = gPad->XtoAbsPixel ( uxmax );
   //     Int_t pymin = gPad->YtoAbsPixel ( uymin );
   //     Int_t pymax = gPad->YtoAbsPixel ( uymax );
// //     if(pxold) gVirtualX->DrawLine(pxold,pymin,pxold,pymax);
// //     else gVirtualX->DrawLine(px,pymin,px,pymax);
//     gPad->SetUniqueID ( px );

   Float_t upx = gPad->AbsPixeltoX(px);
//     Float_t upy = gPad->AbsPixeltoY(py);

   Float_t x = gPad->PadtoX(upx);
//     Float_t y = gPad->PadtoY ( upy );

   Int_t binX = hSelected->GetXaxis()->FindBin(x) - 1;
//     Int_t binY = hSelected->GetYaxis()->FindBin(y)-1;



   // return in case of same bin
   if (uid == binX) return;
//     Printf("%d %d",uid,binX);

   //create or set the new canvas cInfo
   TPaveText *text = 0;
   TVirtualPad *padsav = gPad;
   TCanvas *cInfo = (TCanvas *) gROOT->GetListOfCanvases()->FindObject("cMain");
   if (cInfo) {
      text = (TPaveText *)cInfo->GetPrimitive("mixInfoText");
      if (!text) {
         text = new TPavesText(0.05, 0.05, 0.95, 0.95, 1);
      } else {
         text->DeleteText();
      }

   } else   cInfo = new TCanvas("cInfo", "MixInfo Canvas", 510, 0, 350, 150);

   TVirtualPad *upperPad = cInfo->cd(1);
//     TVirtualPad *upperPadL = upperPad->cd(1);
//     TVirtualPad *upperPadR = upperPad->cd(2);
//     TVirtualPad *upperPadR1 = upperPadR->cd(1);
//     TVirtualPad *upperPadR2 = upperPadR->cd(2);
//     TH1I *hMain = 0;
//     TH1I *hMix = 0;


//     mixInfo->Print();
//     return;

   // gets corresponding histograms
   TH1I *hMain = mixInfo->GetHistogramByType(kMainEvents);
   if (!hMain) {
      Printf("hMain is null");
      return;
   }
   TH1I *hMix =  mixInfo->GetHistogramByType(kMixedEvents);
   if (!hMix) {
      Printf("hMix is null");
      return;
   }

   Double_t numMain = hMain->GetBinContent(binX);
   Double_t numMix = hMix->GetBinContent(binX);
   Int_t hist2DValue = (Int_t) hSelected->GetBinContent(binX + 1, 1);

//    Int_t mixNum = 1;
   if (text) {
      if (mixInfo) {
         AliMixEventPool *evPool = (AliMixEventPool *) mixInfo->GetEventPool("mixEventPool");
         if (evPool) {
//             mixNum = evPool->GetMixNumber();
            if (binX - 1 >= 0) {
               if (!evPool->SetCutValuesFromBinIndex(binX - 1)) return;
            }
            text->SetName("mixInfoText");
            text->SetTextAlign(12);
            text->SetToolTipText("Mixing Info about current binX");
            text->SetBorderSize(2);
            text->AddText(Form("binX=%d", binX));
            text->AddText(Form("numMain=%.0f", numMain));
            text->AddText(Form("numMix=%.0f", numMix));
            text->AddText(Form("BINCONTENT=%d", hist2DValue));
            TObjArray *eventCuts = evPool->GetListOfEventCuts();
            if (eventCuts) {
               TObjArrayIter next(eventCuts);
               AliMixEventCutObj *cut;
               while ((cut = (AliMixEventCutObj *) next())) {
                  if (hist2DValue > 1) text->AddText(Form("%s <%.2f,%.2f)", cut->GetCutName(), cut->GetMin(), cut->GetMax()));
                  else text->AddText(Form("%s <Out of Range>", cut->GetCutName()));
               }
            }
         }
      }
      switch (hist2DValue) {
         case 1 :
            text->SetFillColor(kYellow);
            break;
         case 2 :
            text->SetFillColor(kViolet);
            break;
         case 3 :
            text->SetFillColor(kRed);
            break;
         case 4 :
            text->SetFillColor(kGreen);
            break;
         default:
            text->SetFillColor(kWhite);
            break;
      }
      upperPad->cd();
      text->Draw();
//         upperPadR1->cd();
//         TH1D *proj1 = hSelected->ProjectionY("_xxx",binX);
//         proj1->Draw();
//         upperPadR2->cd();
//         TH1D *proj2 = hSelected->ProjectionY("_xxx",binX);
//         proj1->Draw();

   }
   cInfo->Update();
   padsav->cd();

   gPad->SetUniqueID(binX);
}

//_________________________________________________________________________________________________
Long64_t AliMixInfo::Merge(TCollection *list)
{
   //
   // Merge function
   //
   if (!list) return 0;
   TIter nxfc(list);
   AliMixInfo *mi = 0;
   Long64_t counter = 0;
   while ((mi = (AliMixInfo *) nxfc())) {
      // Do not merge with ourself
      if (mi == this) continue;
      // Make sure that it is a AliMixInfo
      if (!mi->InheritsFrom(AliMixInfo::Class())) {
         Error("Merge", "attempt to add object of class: %s to a %s", mi->ClassName(), ClassName());
         return -1;
      }
      // Merge now
      Add(mi);
      counter++;
   }
   // Done
   return counter;
}

TH1I *AliMixInfo::GetHistogramByType(Int_t index) const
{
   //
   // GetHistogramByType
   //
   return (TH1I *) fHistogramList->FindObject(GetNameHistogramByType(index));
}

//_________________________________________________________________________________________________
void AliMixInfo::Add(AliMixInfo *mi)
{
   //
   // adds AliMixInfo
   //

//    AliInfo(Form("Adding %p", mi));
   if (!mi) return;
   if (!fHistogramList) return;
   TH1I *hMain = GetHistogramByType(kMainEvents);
   if (!hMain) {
      AliError("hMain is null");
      return;
   }
   TH1I *hMix =  GetHistogramByType(kMixedEvents);
   if (!hMix) {
      AliError("hMain is null");
      return;
   }
   hMain->Add(mi->GetHistogramByType(kMainEvents));
   hMix->Add(mi->GetHistogramByType(kMixedEvents));
}

//_________________________________________________________________________________________________
void AliMixInfo::SetEventPool(AliMixEventPool *evPool)
{
   //
   // Sets event pool
   //
   if (!evPool) return;

   if (!fHistogramList) return;

   fHistogramList->Add(evPool);
}

//_________________________________________________________________________________________________
AliMixEventPool *AliMixInfo::GetEventPool(const char *name)
{
   //
   // Gets event pool
   //
   if (!fHistogramList) return 0;

   return (AliMixEventPool *) fHistogramList->FindObject(name);
}

