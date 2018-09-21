/// \ingroup STAT/test
/// \brief  test methods of AliPainter
/// Example usage
/*
\code
.L $AliRoot_SRC/STAT/test/AliPainterTest.C+
AliPainterTest();
root.exe -b -q  $AliRoot_SRC/STAT/test/AliPainterTest.C+ 2>&1 | tee AliPainterTest.log
\endcode
*/
// NOTE: if methods in AliPainter.h will be private we should create new class AliPainterTest inherits from Alipainter. @bdrum
// NOTE: we can create class AliTest which could use for testing.
#include "AliPainter.h"
#include "TError.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "Rtypes.h"
#include "AliDrawStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include <cstring>
#include "TH1.h"
#include "THn.h"

void AliPainterTest_GetNextPad();
//void AliPainterTest_DivideTPad();
//void AliPainterTest_SetMultiGraphTimeAxisTest();
//void AliPainterTest_DrawHistogram();
//void AliPainterTest_SetLimits();
void AliPainterTest_GenerateDoxyImages();


void AliPainterTest() {
  AliPainterTest_GetNextPad();
  //AliPainterTest_DivideTPad();
//  AliPainterTest_SetMultiGraphTimeAxisTest();
  //AliPainterTest_DrawHistogram();
  //AliPainterTest_SetLimits();
  AliPainterTest_GenerateDoxyImages();
}

//void AliPainterTest_DivideTPad() {
//  TCanvas *canvasQA = new TCanvas("canvasQATest", "canvasQATest", 1200, 800);
//  AliPainter::DivideTPad(canvasQA, "<horizontal>[1b,1m,1r300px,1lst0.3]");
//  canvasQA->Print("canvasQADivideTPadTest.xml");
//  canvasQA->Print("canvasQADivideTPadTestFixed.xml");
//
//  auto nDiff = gSystem->GetFromPipe("diff canvasQADivideTPadTest.xml $AliRoot_SRC/STAT/test/canvasQADivideTPadTestFixed.xml  | wc -l").Atoi();
//  if (nDiff - 6 <= 0) {
//    ::Info("AliPainterTest","AliPainter::DivideTPad(\"canvasQATest\",\"<horizontal>[1b,1m,1m,1lst0.3]\",\"test\")- IsOK");
//  }else{
//    ::Error("AliPainterTest","AliDrawStyle::DivideTPad(\"canvasQATest\",\"<horizontal>[1b,1m,1m,1lst0.3]\",\"test\")- FAILED");
//  }
//}
//
//void AliPainterTest_DrawHistogram() {
//  TFile::SetCacheFileDir(".");
//  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
//  TTree *tree = (TTree *) finput->Get("hisPtAll");
//  TObjArray *hisArray = new TObjArray();
//  TList *keys = finput->GetListOfKeys();
//  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
//    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
//    hisArray->AddLast(o);
//  }
//  TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
//  AliPainter::DivideTPad(canvasQA, "<horizontal>[1b,1t,1,1]", "Canvas41");
//  canvasQA->cd(1);
//  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
//  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E0,class=PtIts)");
//  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E0,class=Tgl)");
//  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E)");
//  canvasQA->Print("canvasQADrawHistogramTest.xml");
//
//  auto nDiff = gSystem->GetFromPipe("diff canvasQADrawHistogramTest.xml $AliRoot_SRC/STAT/test/canvasQADrawHistogramTestFixed.xml  | wc -l").Atoi();
//  if (nDiff - 6 <= 0) {
//    ::Info("AliPainterTest",
//           "AliPainter::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)\",hisArray)- IsOK");
//  } else {
//    ::Error("AliPainterTest",
//            "AliPainter::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)\",hisArray)- FAILED");
//  }
//}

void AliPainterTest_GetNextPad() {
  TCanvas *canvasQA1 = new TCanvas("canvasQA1", "canvasQA1", 1200, 800);
  canvasQA1->Divide(4,1);
  TPad *pad1 = (TPad *) canvasQA1->cd(1);
  TPad *tPad = AliPainter::GetNextPad(pad1);
  if (std::strncmp(tPad->GetName(), "canvasQA1_2",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", pad1->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", pad1->GetName());
  }
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_3",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_4",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }

  TPad *pad3 = (TPad *) canvasQA1->cd(2);
  pad3->Divide(2,1);
  tPad = (TPad *) pad3->cd(1);
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_2_2",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_3",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }
}

//void AliPainterTest_SetLimits() {
//  TFile::SetCacheFileDir(".");
//  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
//  TTree *tree = (TTree *) finput->Get("hisPtAll");
//  TObjArray *hisArray = new TObjArray();
//  TList *keys = finput->GetListOfKeys();
//  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
//    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
//    hisArray->AddLast(o);
//  }
//  THn *hisN = (THn *) hisArray->FindObject("hisK0DMassQPtTgl");
//  TObjArray *keep = NULL;
//  TObjArray *hisArr = AliPainter::PrepareHistogram(hisN, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0,ylim=[<mean>,<max>])", keep);
//  Long64_t comSize = 0;
//  Double_t *combineArray = AliPainter::GetDataArray(hisArr,comSize,4);
//  Double_t mean = TMath::Mean(comSize, combineArray);
//  Double_t max = TMath::MaxElement(comSize, combineArray);
//  AliPainter::SetLimits(hisArr,4);
//
//  if ((Int_t) mean == 2028 && (Int_t) max == 10730) {
//    ::Info("AliPainterTest",
//           "AliPainter::SetLimits()- IsOK");
//  } else {
//    ::Error("AliPainterTest",
//            "AliPainter::SetLimits()- FAILED");
//  }
//
//}

void AliPainterTest_GenerateDoxyImages() {
  TFile::SetCacheFileDir(".");
  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
  TTree *tree = (TTree *) finput->Get("hisPtAll");
  TObjArray *hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));

  TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1,1,1]", "Canvas41");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E0,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E0,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");

  if (gSystem->GetFromPipe("find $AliRoot_SRC/STAT/ -name imgdoc  | wc -l").Atoi() == 0)
    gSystem->GetFromPipe("mkdir $AliRoot_SRC/STAT/imgdoc");

  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example1.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<vertical>[1,1,1,1]", "Canvas41");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E0,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E0,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example2.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<vertical>[1r,1l,1r,1l]", "Raw,Error", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E0,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E0,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example3.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<vertical>[1,3m,2m,1]", "Raw,Error", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E0,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E0,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=0E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(div=1,class=Mass,drawOpt=E0)");
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example4.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1rpx200,1,1,1]", "Raw,Error", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E0,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E0,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E0,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example5.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl()(0)()()");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example6.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80)(0)()()");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example7.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0,ylim=[0,])");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example8.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0,ylim=[0,], div=1)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example9.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E0,ylim=[<min>+0.5*<mean>,<max>-0.5*<mean>], div=1)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example10.png");
  canvasQA->Clear();
}