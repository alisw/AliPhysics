#include <TGButton.h>
#include <TGComboBox.h>
#include <TEveBrowser.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>
#include <TApplication.h>
#include <TGComboBox.h>
#include <TLatex.h>
#include <TGSlider.h>
#include "TGDoubleSlider.h"
#include <TEvePointSet.h>

/*
class CutsDoubleSlider : public TGMainFrame {

protected:
   TGLabel*		gLabel;
   TGCheckButton*	gApplyCut;
   TGNumberEntry*	gMaxValue;
   TGNumberEntry*	gMinValue;
   TGDoubleSlider*	gDoubleSlider;

public:
   CutsDoubleSlider((const TGWindow* p, Int_t sliderWidth, const char* labelText, Int_t valMin, Int_t valMax, Pixel_t backgroundColor, Pixel_t textColor);

   void GetLabel() { return gLabel; };
   void GetCheckButton() { return gApplyCut; };
   void GetEntryMax() { return gMaxValue; };
   void GetEntryMin() { return gMinValue; };
   void GetSlider() { return gDoubleSlider; };

   ClassDef(CutsDoubleSlider, 0)
};

CutsDoubleSlider::CutsDoubleSlider(const TGWindow* p, Int_t sliderWidth, const char* labelText, Int_t valMin, Int_t valMax, Pixel_t backgroundColor, Pixel_t textColor)
 : TGMainFrame(p, 10, 10, kHorizontalFrame)
{

   TGVerticalFrame *vframe = new TGVerticalFrame(this);

   TGHorizontalFrame* hframe = new TGHorizontalFrame(vframe, 200, 20, kFixedWidth);
   TGLabel *label = 0;

   gLabel = new TGLabel(hframe,labelText);
   gLabel->SetBackgroundColor(backgroundColor);
   gLabel->SetTextColor(textColor);

   gMinValue = new TGNumberEntry(hframe, 0, 6);
   gMinValue->SetNumber(valMin);
   gMinValue->SetBackgroundColor(backgroundColor);

   gMaxValue = new TGNumberEntry(hframe, 0, 6);
   gMaxValue->SetNumber(valMax);
   gMaxValue->SetBackgroundColor(backgroundColor);

   gApplyCut = new TGCheckButton(hframe, "", 10);
   gApplyCut->SetBackgroundColor(backgroundColor);

   hframe->AddFrame(gLabel, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   hframe->AddFrame(gMinValue, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   hframe->AddFrame(gMaxValue, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   hframe->AddFrame(gApplyCut, new TGLayoutHints(kLHintsRight, 4, 0, 0, 0));
   hframe->SetBackgroundColor(backgroundColor);

   vframe->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

   hframe = new TGHorizontalFrame(vframe, 200, 20, kFixedWidth);

   gDoubleSlider = new TGDoubleHSlider(hframe, sliderWidth, 1, -1, kHorizontalFrame, backgroundColor, kFALSE, kFALSE);
   gDoubleSlider->SetRange(valMin, valMax);
   gDoubleSlider->SetPosition(valMin, valMax);
   gDoubleSlider->SetBackgroundColor(backgroundColor);

   hframe->AddFrame(gDoubleSlider, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   hframe->SetBackgroundColor(backgroundColor);

   vframe->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 10));

   AddFrame(vframe, new TGLayoutHints(kLHintsCenterX));

   SetBackgroundColor(backgroundColor);

   Resize();

   MapSubwindows();

   SetWMSizeHints(GetDefaultWidth(), GetDefaultHeight(), 1000, 1000, 0 ,0);
   SetWindowName("Double Slider");
   MapRaised();

}

*/

class SetCutsWindow : public TGMainFrame {

protected:
   TGCheckButton* gDrawV0s;
   TGCheckButton* gDrawCascades;
   TGCheckButton* gDrawKinks;
   TGCheckButton* gDrawVertex;
   TGCheckButton* gDrawTracklets;
   TGCheckButton* gDrawTracks;
   TGCheckButton* gDrawClusters;
   TGCheckButton* gDrawTracksType1;
   TGCheckButton* gDrawTracksType2;
   TGCheckButton* gDrawTracksType3;
   TGCheckButton* gDrawTracksType4;
   TGCheckButton* gDrawTracksType5;
   TGCheckButton* gDrawTracksType6;
   TGCheckButton* gDrawTracksType7;
   TGCheckButton* gCutOnP;
   TGCheckButton* gCutOnPt;
   TGCheckButton* gCutOnEta;
   TGCheckButton* gCutOnMult;
   TGCheckButton* gCutOnCls;
   TEveGDoubleValuator* gPRange;
   TEveGDoubleValuator* gPtRange;
   TEveGDoubleValuator* gEtaRange;
   TGHSlider* gMultRange;
   TGNumberEntry* gMultRangeNE;
   TGHSlider* gClsRange;
   TGNumberEntry* gClsRangeNE;
   TGHSlider* gPMVRange;
   TGLabel* gPMVRangeLabel;
   TGComboBox* gVectorMode;
   TGComboBox* gPosColorList;
   TGTextButton* gPosColorButton;
   TGComboBox* gNegColorList;
   TGTextButton* gNegColorButton;
   TGComboBox* gBkgColorList;
   TGTextButton* gBkgColorButton;
   TGLOverlayButton *gOverlayButton3D;
   TGLOverlayButton *gOverlayButtonRPhi;
   TGLOverlayButton *gOverlayButtonRhoZ;
   Bool_t gDrawHistograms[12];
   
//   CutsDoubleSlider *slajder1;
//   CutsDoubleSlider *slajder2;
//   CutsDoubleSlider *slajder3;

public:
   SetCutsWindow();
   void MultNECallBack();
   void MultSliderCallBack();
   void ClsNECallBack();
   void ClsSliderCallBack();
   void PosTracksCallBack();
   void NegTracksCallBack();
   void BackgroundCallBack();
   void PMVSliderCallBack();
   void AddDescriptions();
   void MuonGeometry();
   void DefaultGeometry();
   void BrightGeometry();
   void TransparentGeometry();
   void YellowGeometry();
   void SaveAllViews();
   void Save3DView();
   void SaveRPhiView();
   void DrawPtHisto();
   void DrawEtaHisto();
   void DrawPhiHisto();
   void DrawPhiPtHisto();
   void DrawPtYHisto();
   void DrawEtaPhiHisto();
   void DrawPtHistoAll();
   void DrawEtaHistoAll();
   void DrawPhiHistoAll();
   void DrawPhiPtHistoAll();
   void DrawPtYHistoAll();
   void DrawEtaPhiHistoAll();
   void DrawHistos();
   void DrawHistosAll();
   void DrawResiduals();
   Int_t GetTrackColorByMomentum(Double_t, Int_t);
   void SetStandardCuts();
   void AddMomentumVectors();
   void SetCuts();
   void CloseTab();
   
   ClassDef(SetCutsWindow, 0)
};

//________________________________________________

namespace
{

   const char *gPictureSaveAsTypes[] = {"PNG Image", "*.png", 0, 0}; //for saving pictures

}

//________________________________________________

SetCutsWindow::SetCutsWindow() : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame)
{

   // Main test window.

   gEve->GetWindowManager()->HideAllEveDecorations();

   SetCleanup(kDeepCleanup);

   TGTextButton *b = 0;
   TGLabel *label = 0;
   TGHorizontalFrame *hframe = 0;
   TGHorizontalFrame *hframeMerge = 0;
   TGHorizontalFrame *hframe1 = 0;
   TGHorizontalFrame *hframe2 = 0;
   TGVerticalFrame *vframe = 0;

   gOverlayButton3D = 0;
   gOverlayButtonRPhi = 0;
   gOverlayButtonRhoZ = 0;

   for(Int_t i = 0; i < 12; i++)
      gDrawHistograms[i] = kFALSE;

   TEveBrowser *browser = gEve->GetBrowser();

   TGShutter *mainShutter = new TGShutter(this, kSunkenFrame);

   TGShutterItem *item1 = new TGShutterItem(mainShutter, new TGHotString("Draw Objects"), 1);

   TGShutterItem *item2 = new TGShutterItem(mainShutter, new TGHotString("Track Cuts"), 2);

   TGShutterItem *item3 = new TGShutterItem(mainShutter, new TGHotString("Colors"), 3);

   TGShutterItem *item4 = new TGShutterItem(mainShutter, new TGHotString("Geometry"), 4);

   TGShutterItem *item5 = new TGShutterItem(mainShutter, new TGHotString("Analysis"), 5);

   TGShutterItem *item6 = new TGShutterItem(mainShutter, new TGHotString("Momentum Vectors"), 6);

   mainShutter->AddItem(item1);

   mainShutter->AddItem(item2);

   mainShutter->AddItem(item3);

   mainShutter->AddItem(item4);

   mainShutter->AddItem(item5);

   mainShutter->AddItem(item6);

   AddFrame(mainShutter, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   TGCompositeFrame *container1 = (TGCompositeFrame *) item1->GetContainer();

   TGCompositeFrame *container2 = (TGCompositeFrame *) item2->GetContainer();

   TGCompositeFrame *container3 = (TGCompositeFrame *) item3->GetContainer();

   TGCompositeFrame *container4 = (TGCompositeFrame *) item4->GetContainer();

   TGCompositeFrame *container5 = (TGCompositeFrame *) item5->GetContainer();

   TGCompositeFrame *container6 = (TGCompositeFrame *) item6->GetContainer();

   // Draw Elements

   TGVerticalFrame *drawElements = new TGVerticalFrame(container1);
   container1->AddFrame(drawElements, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));//, 5, 5, 5, 5));// | kLHintsExpandY

   // Track Cuts
   TGVerticalFrame *trackCuts = new TGVerticalFrame(container2);
   container2->AddFrame(trackCuts, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));//, 5, 5, 5, 5));// | kLHintsExpandY

   // Buttons
   TGVerticalFrame *buttons = new TGVerticalFrame(container3);
   container3->AddFrame(buttons, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));//, 5, 5, 5, 5));// | kLHintsExpandY

   // Geometry
   TGVerticalFrame *geometry = new TGVerticalFrame(container4);
   container4->AddFrame(geometry, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));//, 5, 5, 5, 5));// | kLHintsExpandY

   // Analysis
   TGVerticalFrame *analysis = new TGVerticalFrame(container5);
   container5->AddFrame(analysis, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));//, 5, 5, 5, 5));// | kLHintsExpandY

   // Momentum Vectors
   TGVerticalFrame *momentumVectors = new TGVerticalFrame(container6);
   container6->AddFrame(momentumVectors, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));//, 5, 5, 5, 5));// | kLHintsExpandY

   // DRAW ELEMENTS

   separator = new TGHorizontal3DLine(drawElements);
   drawElements->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(drawElements, "ESD objects");
   drawElements->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(drawElements);
   drawElements->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // V0s

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "V0s");
   gDrawV0s = new TGCheckButton(hframeMerge, "", 10);
   gDrawV0s->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawV0s);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Cascades

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Cascades");
   gDrawCascades = new TGCheckButton(hframeMerge, "", 10);
   gDrawCascades->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawCascades);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Kinks

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Kinks");
   gDrawKinks = new TGCheckButton(hframeMerge, "", 10);
   gDrawKinks->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawKinks);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Primary Vertex

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Primary Vertex");
   gDrawVertex = new TGCheckButton(hframeMerge, "", 10);
   gDrawVertex->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawVertex);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Tracklets

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Tracklets");
   gDrawTracklets = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracklets->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracklets);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Tracks

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Tracks");
   gDrawTracks = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracks->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracks);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Clusters

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Clusters");
   gDrawClusters = new TGCheckButton(hframeMerge, "", 10);
   gDrawClusters->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawClusters);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // TRACK TYPES

   separator = new TGHorizontal3DLine(drawElements);
   drawElements->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(drawElements, "Track types");
   drawElements->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(drawElements);
   drawElements->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // sigma < 3

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "Sigma < 3");
   gDrawTracksType1 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType1->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType1);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // 3 < sigma < 5

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "3 < Sigma < 5");
   gDrawTracksType2 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType2->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType2);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // 5 < Sigma

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "5 < Sigma");
   gDrawTracksType3 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType3->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType3);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // no ITS refit; 

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "no ITS refit");
   gDrawTracksType4 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType4->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType4);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // no TPC refit

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "no TPC refit");
   gDrawTracksType5 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType5->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType5);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // ITS ncl>=3

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "ITS ncl>=3");
   gDrawTracksType6 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType6->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType6);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // ITS others

   hframeMerge = new TGHorizontalFrame(drawElements, 200, 20, kFixedWidth);

   label = new TGLabel(hframeMerge, "ITS others");
   gDrawTracksType7 = new TGCheckButton(hframeMerge, "", 10);
   gDrawTracksType7->SetEnabled(kTRUE);
   hframeMerge->AddFrame(label, new TGLayoutHints(kLHintsExpandX));
   hframeMerge->AddFrame(gDrawTracksType7);

   drawElements->AddFrame(hframeMerge, new TGLayoutHints(kLHintsExpandX));

   // Main menu

   separator = new TGHorizontal3DLine(drawElements);
   drawElements->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(drawElements, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Apply Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   drawElements->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(drawElements, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Close Tab");
   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   drawElements->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // TRACK CUTS

   separator = new TGHorizontal3DLine(trackCuts);
   trackCuts->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   gPRange = new TEveGDoubleValuator(hframe,"P range:", 40, 0);
   gPRange->SetNELength(6);
   gPRange->SetLabelWidth(50);
   gPRange->Build();
   gPRange->GetSlider()->SetWidth(180);
   gPRange->SetLimits(0, 5, TGNumberFormat::kNESRealTwo);
   gPRange->SetValues(0, 5, TGNumberFormat::kNESRealTwo);

   gCutOnP = new TGCheckButton(hframe, "", 10);
   gCutOnP->SetEnabled(kTRUE);

   hframe->AddFrame(gPRange, new TGLayoutHints(kLHintsExpandX, 10, 10, 10, 10));
   hframe->AddFrame(gCutOnP, new TGLayoutHints(kLHintsNormal, 10, 10, 10, 10));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // Pt

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   gPtRange = new TEveGDoubleValuator(hframe,"Pt range:", 40, 0);
   gPtRange->SetNELength(6);
   gPtRange->SetLabelWidth(50);
   gPtRange->Build();
   gPtRange->GetSlider()->SetWidth(180);
   gPtRange->SetLimits(0, 5, TGNumberFormat::kNESRealTwo);
   gPtRange->SetValues(0, 5, TGNumberFormat::kNESRealTwo);

   gCutOnPt = new TGCheckButton(hframe, "", 10);
   gCutOnPt->SetEnabled(kTRUE);

   hframe->AddFrame(gPtRange, new TGLayoutHints(kLHintsExpandX, 10, 10, 10, 10));
   hframe->AddFrame(gCutOnPt, new TGLayoutHints(kLHintsNormal, 10, 10, 10, 10));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // Eta

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   gEtaRange = new TEveGDoubleValuator(hframe,"Eta range:", 40, 0);
   gEtaRange->SetNELength(6);
   gEtaRange->SetLabelWidth(50);
   gEtaRange->Build();
   gEtaRange->GetSlider()->SetWidth(180);
   gEtaRange->SetLimits(-5, 5, TGNumberFormat::kNESRealTwo);
   gEtaRange->SetValues(-5, 5, TGNumberFormat::kNESRealTwo);

   gCutOnEta = new TGCheckButton(hframe, "", 10);
   gCutOnEta->SetEnabled(kTRUE);

   hframe->AddFrame(gEtaRange, new TGLayoutHints(kLHintsExpandX, 10, 10, 10, 10));
   hframe->AddFrame(gCutOnEta, new TGLayoutHints(kLHintsNormal, 10, 10, 10, 10));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // Multiplicity

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "% Tracks displayed:");

   gMultRangeNE = new TGNumberEntry(hframe, 0, 6);
   gMultRangeNE->SetNumber(100);
   gCutOnMult = new TGCheckButton(hframe, "", 10);
   gCutOnMult->SetEnabled(kTRUE);

   hframe->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   hframe->AddFrame(gMultRangeNE, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   hframe->AddFrame(gCutOnMult, new TGLayoutHints(kLHintsRight, 10, 10, 0, 0));//kLHintsNormal

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   gMultRange = new TGHSlider(hframe,180);
   gMultRange->SetRange(0, 100);
   gMultRange->SetPosition(100);
   gMultRange->Connect("PositionChanged(Int_t)", "SetCutsWindow", this, "MultSliderCallBack()");

   gMultRangeNE->Connect("ValueSet(Long_t)", "SetCutsWindow", this, "MultNECallBack()");

   hframe->AddFrame(gMultRange, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 10));

   // TPC Clusters

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "TPC clusters:");

   gClsRangeNE = new TGNumberEntry(hframe, 0, 6);
   gClsRangeNE->SetNumber(0);
   gClsRangeNE->Connect("ValueSet(Long_t)", "SetCutsWindow", this, "ClsNECallBack()");

   gCutOnCls = new TGCheckButton(hframe, "", 10);
   gCutOnCls->SetEnabled(kTRUE);

   hframe->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   hframe->AddFrame(gClsRangeNE, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));
   hframe->AddFrame(gCutOnCls, new TGLayoutHints(kLHintsRight, 10, 10, 0, 0));//kLHintsNormal

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

   hframe = new TGHorizontalFrame(trackCuts, 200, 20, kFixedWidth);

   gClsRange = new TGHSlider(hframe,180);
   gClsRange->SetRange(0, 200);
   gClsRange->SetPosition(0);
   gClsRange->Connect("PositionChanged(Int_t)", "SetCutsWindow", this, "ClsSliderCallBack()");

   hframe->AddFrame(gClsRange, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 10));

   separator = new TGHorizontal3DLine(trackCuts);
   trackCuts->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Standard cuts

   hframe = new TGHorizontalFrame(trackCuts, 100, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Standard Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetStandardCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 10, 10, 10, 10));

   // Main menu

   separator = new TGHorizontal3DLine(trackCuts);
   trackCuts->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(trackCuts, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Apply Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(trackCuts, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Close Tab");
   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   trackCuts->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

//    trackCutsITS->SetMaxDCAToVertexXY(2.4);
//    trackCutsITS->SetMaxDCAToVertexZ(3.2);


   // BUTTONS

   separator = new TGHorizontal3DLine(buttons);
   buttons->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Positive tracks colorset

   hframe = new TGHorizontalFrame(buttons, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "Positive tracks:");
   hframe->AddFrame(label, new TGLayoutHints(kLHintsExpandX));

   gPosColorList = new TGComboBox(hframe);
   gPosColorList->AddEntry("Blue", 1);
   gPosColorList->AddEntry("Red", 2);
   gPosColorList->AddEntry("Mixed", 3);
   gPosColorList->Resize(100,20);
   gPosColorList->Connect("Selected(Int_t)", "SetCutsWindow", this, "PosTracksCallBack()");
   hframe->AddFrame(gPosColorList, new TGLayoutHints(kLHintsExpandX));

   gPosColorButton = new TGTextButton(hframe, "    ");

   hframe->AddFrame(gPosColorButton, new TGLayoutHints(kLHintsNormal, 5, 5, 1, 1));

   buttons->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Negative tracks colorset

   hframe = new TGHorizontalFrame(buttons, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "Negative tracks:");
   hframe->AddFrame(label, new TGLayoutHints(kLHintsExpandX));

   gNegColorList = new TGComboBox(hframe);
   gNegColorList->AddEntry("Blue", 1);
   gNegColorList->AddEntry("Red", 2);
   gNegColorList->AddEntry("Mixed", 3);
   gNegColorList->Resize(100,20);
   gNegColorList->Connect("Selected(Int_t)", "SetCutsWindow", this, "NegTracksCallBack()");
   hframe->AddFrame(gNegColorList, new TGLayoutHints(kLHintsExpandX));

   gNegColorButton = new TGTextButton(hframe, "    ");
   hframe->AddFrame(gNegColorButton, new TGLayoutHints(kLHintsNormal, 5, 5, 1, 1));

   buttons->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Background color

   hframe = new TGHorizontalFrame(buttons, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "Background:");
   hframe->AddFrame(label, new TGLayoutHints(kLHintsExpandX));

   gBkgColorList = new TGComboBox(hframe);
   gBkgColorList->AddEntry("White", 1);
   gBkgColorList->AddEntry("Black", 2);
   gBkgColorList->Resize(100,20);
   gBkgColorList->Connect("Selected(Int_t)", "SetCutsWindow", this, "BackgroundCallBack()");
   hframe->AddFrame(gBkgColorList, new TGLayoutHints(kLHintsExpandX));

   gBkgColorButton = new TGTextButton(hframe, "    ");
   hframe->AddFrame(gBkgColorButton, new TGLayoutHints(kLHintsNormal, 5, 5, 1, 1));

   buttons->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Main menu

   separator = new TGHorizontal3DLine(buttons);
   buttons->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(buttons, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Apply Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   buttons->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(buttons, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Close Tab");
   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   buttons->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // GEOMETRY

   separator = new TGHorizontal3DLine(geometry);
   geometry->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(geometry, "Geometries");
   geometry->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Geometry

   separator = new TGHorizontal3DLine(geometry);
   geometry->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "MUON arm");
   b->Connect("Clicked()", "SetCutsWindow", this, "MuonGeometry()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Default");
   b->Connect("Clicked()", "SetCutsWindow", this, "DefaultGeometry()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Bright");
   b->Connect("Clicked()", "SetCutsWindow", this, "BrightGeometry()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Transparent");
   b->Connect("Clicked()", "SetCutsWindow", this, "TransparentGeometry()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Yellow");
   b->Connect("Clicked()", "SetCutsWindow", this, "YellowGeometry()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

//   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

//   b = new TGTextButton(hframe, "Close Tab");
//   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");
//   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

//   b = new TGTextButton(hframe, "Close Tab");
//   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");
//   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

//   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Snapshots

   separator = new TGHorizontal3DLine(geometry);
   geometry->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(geometry, "Snapshots");
   geometry->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(geometry);
   geometry->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

//   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

//   b = new TGTextButton(hframe, "Descriptions");
//   b->Connect("Clicked()", "SetCutsWindow", this, "AddDescriptions()");
//   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

//   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Save all Views");
   b->Connect("Clicked()", "SetCutsWindow", this, "SaveAllViews()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "3D View");
   b->Connect("Clicked()", "SetCutsWindow", this, "Save3DView()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "RPhi");
   b->Connect("Clicked()", "SetCutsWindow", this, "SaveRPhiView()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "RhoZ");
   b->Connect("Clicked()", "SetCutsWindow", this, "SaveRhoZView()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Main menu

   separator = new TGHorizontal3DLine(geometry);
   geometry->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Apply Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(geometry, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Close Tab");
   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   geometry->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // ANALYSIS

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(analysis, "Single Event Histograms");
   analysis->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Single Event

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Pt");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPtHisto()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Eta");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawEtaHisto()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Phi");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPhiHisto()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Phi-Pt");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPhiPtHisto()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Pt-Y");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPtYHisto()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Eta-Phi");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawEtaPhiHisto()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Draw Histograms");
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawHistos()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // All Events

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(analysis, "All Events Histograms");
   analysis->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Pt");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPtHistoAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Eta");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawEtaHistoAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Phi");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPhiHistoAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Phi-Pt");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPhiPtHistoAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Pt-Y");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawPtYHistoAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   b = new TGTextButton(hframe, "Eta-Phi");
   b->AllowStayDown(kTRUE);
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawEtaPhiHistoAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));//, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Draw Histograms");
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawHistosAll()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(analysis, "Residuals");
   analysis->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(analysis, 200, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Draw Residuals");
   b->Connect("Clicked()", "SetCutsWindow", this, "DrawResiduals()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Main menu

   separator = new TGHorizontal3DLine(analysis);
   analysis->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(analysis, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Apply Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(analysis, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Close Tab");
   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   analysis->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // MOMENTUM VECTORS

   // Draw momentum vectors

   separator = new TGHorizontal3DLine(momentumVectors);
   momentumVectors->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   label = new TGLabel(momentumVectors, "Draw Momentum Vectors");
   momentumVectors->AddFrame(label, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   separator = new TGHorizontal3DLine(momentumVectors);
   momentumVectors->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(momentumVectors, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "Scale:");
   hframe->AddFrame(label, new TGLayoutHints(kLHintsNormal, 5, 5, 0, 0));

   gVectorMode = new TGComboBox(hframe);
   gVectorMode->AddEntry("Log scale", 1);
   gVectorMode->AddEntry("maxP -> 600cm", 2);
   gVectorMode->AddEntry("1GeV -> 100cm", 3);
   gVectorMode->Resize(150,20);
   gVectorMode->Select(2);

   hframe->AddFrame(gVectorMode, new TGLayoutHints(kLHintsNormal, 0, 5, 0, 0));

   momentumVectors->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(momentumVectors, 200, 20, kFixedWidth);

   label = new TGLabel(hframe, "Minimal P: ");
   hframe->AddFrame(label, new TGLayoutHints(kLHintsNormal, 5, 5, 0, 0));

   gPMVRangeLabel = new TGLabel(hframe, "0.0");
   hframe->AddFrame(gPMVRangeLabel, new TGLayoutHints(kLHintsNormal, 0, 0, 0, 0));

   label = new TGLabel(hframe, "GeV/c");
   hframe->AddFrame(label, new TGLayoutHints(kLHintsNormal, 1, 0, 0, 0));

   momentumVectors->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(momentumVectors, 200, 20, kFixedWidth);

   gPMVRange = new TGHSlider(hframe,180);
   gPMVRange->SetRange(0, 50);
   gPMVRange->SetPosition(0);
   gPMVRange->Connect("PositionChanged(Int_t)", "SetCutsWindow", this, "PMVSliderCallBack()");

   hframe->AddFrame(gPMVRange, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   momentumVectors->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(momentumVectors, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Draw");
   b->Connect("Clicked()", "SetCutsWindow", this, "AddMomentumVectors()");
   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 0, 0));

   momentumVectors->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   // Main menu

   separator = new TGHorizontal3DLine(momentumVectors);
   momentumVectors->AddFrame(separator, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(momentumVectors, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Apply Cuts");
   b->Connect("Clicked()", "SetCutsWindow", this, "SetCuts()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX));

   momentumVectors->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   hframe = new TGHorizontalFrame(momentumVectors, 150, 20, kFixedWidth);

   b = new TGTextButton(hframe, "Close Tab");
   b->Connect("Clicked()", "SetCutsWindow", this, "CloseTab()");

   hframe->AddFrame(b, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   momentumVectors->AddFrame(hframe, new TGLayoutHints(kLHintsExpandX));

   // FINAL STUFF

   Resize();

   MapSubwindows();

   SetWMSizeHints(GetDefaultWidth(), GetDefaultHeight(), 1000, 1000, 0 ,0);
   SetWindowName("Pb-Pb set");
   MapRaised();

}

//______________________________________________________________________________

void SetCutsWindow::MultNECallBack()
{

   Double_t entry;

   entry = gMultRangeNE->GetNumber();
   gMultRange->SetPosition(entry);

   return;

}

//______________________________________________________________________________

void SetCutsWindow::MultSliderCallBack()
{

   Double_t entry;

   entry = gMultRange->GetPosition();
   gMultRangeNE->SetNumber(entry);

   return;

}

//______________________________________________________________________________

void SetCutsWindow::ClsNECallBack()
{
M
   Double_t entry;

   entry = gClsRangeNE->GetNumber();
   gClsRange->SetPosition(entry);

   return;

}

//______________________________________________________________________________

void SetCutsWindow::ClsSliderCallBack()
{

   Double_t entry;

   entry = gClsRange->GetPosition();
   gClsRangeNE->SetNumber(entry);

   return;

}

//______________________________________________________________________________

void SetCutsWindow::PosTracksCallBack()
{

   switch(gPosColorList->GetSelected())
   {

      case 1:
         gPosColorButton->SetBackgroundColor(gROOT->GetColor(4)->GetPixel());
         gPosColorButton->SetText("    ");
         break;

      case 2:
         gPosColorButton->SetBackgroundColor(gROOT->GetColor(2)->GetPixel());
         gPosColorButton->SetText("    ");
         break;

      case 3:
         gPosColorButton->SetBackgroundColor(gROOT->GetColor(18)->GetPixel());
         gPosColorButton->SetText(" M ");
         break;

      default:
         gPosColorButton->SetBackgroundColor(gROOT->GetColor(18)->GetPixel());
         gPosColorButton->SetText("    ");
         break;

   }

   return;

}

void SetCutsWindow::NegTracksCallBack()
{

   switch(gNegColorList->GetSelected())
   {

      case 1:
         gNegColorButton->SetBackgroundColor(gROOT->GetColor(4)->GetPixel());
         gNegColorButton->SetText("    ");
         break;

      case 2:
         gNegColorButton->SetBackgroundColor(gROOT->GetColor(2)->GetPixel());
         gNegColorButton->SetText("    ");
         break;

      case 3:
         gNegColorButton->SetBackgroundColor(gROOT->GetColor(18)->GetPixel());
         gNegColorButton->SetText(" M ");
         break;

      default:
         gNegColorButton->SetBackgroundColor(gROOT->GetColor(18)->GetPixel());
         gNegColorButton->SetText("    ");
         break;

   }

   return;

}

void SetCutsWindow::BackgroundCallBack()
{

   switch(gBkgColorList->GetSelected())
   {

      case 1:
         gBkgColorButton->SetBackgroundColor(gROOT->GetColor(0)->GetPixel());
         gBkgColorButton->SetText("    ");

         if(!gEve->GetViewers()->UseLightColorSet())
           gEve->GetViewers()->SwitchColorSet(); //white background

         break;

      case 2:
         gBkgColorButton->SetBackgroundColor(gROOT->GetColor(1)->GetPixel());
         gBkgColorButton->SetText("    ");

         if(gEve->GetViewers()->UseLightColorSet())
            gEve->GetViewers()->SwitchColorSet(); //black background

         break;

      default:
         gBkgColorButton->SetBackgroundColor(gROOT->GetColor(18)->GetPixel());
         gBkgColorButton->SetText("    ");

         break;

   }

   return;

}

//______________________________________________________________________________

void SetCutsWindow::PMVSliderCallBack()
{

//   TGHSlider* gPMVRange;
//   TGLabel* gPMVRangeLabel;


   Double_t entry;

   entry = gPMVRange->GetPosition();
   gPMVRangeLabel->SetText(Form("%2.1f",0.1*entry));

   return;

}

//______________________________________________________________________________

void SetCutsWindow::AddDescriptions()
{
/*
   if(gOverlayButton3D && gOverlayButtonRPhi && gOverlayButtonRhoZ)
   {

      (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("3D View")))->GetGLViewer()->RemoveOverlayElement((TGLOverlayElement*)gOverlayButton3D);
      (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("RPhi View")))->GetGLViewer()->RemoveOverlayElement((TGLOverlayElement*)gOverlayButtonRPhi);
      (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("RhoZ View")))->GetGLViewer()->RemoveOverlayElement((TGLOverlayElement*)gOverlayButtonRhoZ);

   }
   else
   {

      AliEveMultiView *multiView = AliEveMultiView::Instance();

      TGLViewer *glv = (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("3D View")))->GetGLViewer();
      gOverlayButton3D = new TGLOverlayButton(glv,  "ALICE run 123456 event 123", 0, 0, multiView->Get3DView()->GetEveFrame()->GetWidth(), 30);
      gOverlayButton3D->SetAlphaValues(0.5, 0.5);

      glv = (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("RPhi View")))->GetGLViewer();
      gOverlayButtonRPhi = new TGLOverlayButton(glv,  "ALICE run 123456 event 123", 0, 0, multiView->GetRPhiView()->GetEveFrame()->GetWidth(), 30);
      gOverlayButtonRPhi->SetAlphaValues(0.5, 0.5);

      glv = (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("RhoZ View")))->GetGLViewer();
      gOverlayButtonRhoZ = new TGLOverlayButton(glv,  "ALICE run 123456 event 123", 0, 0, multiView->GetRhoZView()->GetEveFrame()->GetWidth(), 30);
      gOverlayButtonRhoZ->SetAlphaValues(0.5, 0.5);

   }

   gEve->Redraw3D(kFALSE, kTRUE);
*/
   return;

}

//______________________________________________________________________________

void SetCutsWindow::MuonGeometry()
{

   if(gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON"))
   {
      if(gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->GetRnrSelf() || gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->GetRnrChildren())
      {
         gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrSelf(kFALSE);
         gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrChildren(kFALSE);
      }
      else
      {
         gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrSelf(kTRUE);
         gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrChildren(kTRUE);
      }  

   }

   gEve->FullRedraw3D();   

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DefaultGeometry()
{

   AliEveMultiView *mv = AliEveMultiView::Instance();

   mv->DestroyAllGeometries(); //destroy RPhi and Rhoz geometries before putting new

   gEve->LoadVizDB("geom_gentle_default.C", kTRUE, kTRUE); //loading geometry

   if(!gEve->GetViewers()->UseLightColorSet())
      gEve->GetViewers()->SwitchColorSet(); //white background

   gEve->FullRedraw3D();   

   return;

}

//______________________________________________________________________________

void SetCutsWindow::BrightGeometry()
{

   AliEveMultiView *mv = AliEveMultiView::Instance();

   mv->DestroyAllGeometries();

   gEve->LoadVizDB("geom_gentle_bright.C", kTRUE, kTRUE);

   if(gEve->GetViewers()->UseLightColorSet())
      gEve->GetViewers()->SwitchColorSet();



   gEve->FullRedraw3D();

   return;

}

//______________________________________________________________________________

void SetCutsWindow::TransparentGeometry()
{

   AliEveMultiView *mv = AliEveMultiView::Instance();

   mv->DestroyAllGeometries();

   gEve->LoadVizDB("geom_gentle_transparent.C", kTRUE, kTRUE);

   if(!gEve->GetViewers()->UseLightColorSet())
      gEve->GetViewers()->SwitchColorSet();

   gEve->FullRedraw3D();

   return;

}

//______________________________________________________________________________

void SetCutsWindow::YellowGeometry()
{

   AliEveMultiView *mv = AliEveMultiView::Instance();

   mv->DestroyAllGeometries();

   gEve->LoadVizDB("geom_gentle_yellow.C", kTRUE, kTRUE);

   if(!gEve->GetViewers()->UseLightColorSet())
      gEve->GetViewers()->SwitchColorSet();

   gEve->FullRedraw3D();    

   return;

}

//______________________________________________________________________________

void SetCutsWindow::SaveAllViews()
{

     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     TString file2(filere[1]);
     TString file3(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += "_3D.png";
       
     if (!file2.EndsWith(".png"))
       file2 += "_RPhi.png";
       
     if (!file3.EndsWith(".png"))
       file3 += "_RhoZ.png";

     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     TEveViewer* view3d = ((TEveViewer*)*i);
     view3d->GetGLViewer()->SavePictureScale(file1,4.0); // getting high resolution
     i++;
     TEveViewer* viewrphi = ((TEveViewer*)*i);
     viewrphi->GetGLViewer()->SavePictureScale(file2,4.0);
     i++;
     TEveViewer* viewrhoz = ((TEveViewer*)*i);
     viewrhoz->GetGLViewer()->SavePictureScale(file3,4.0);
     
     printf("Done.\n"); 

   return;

}

//______________________________________________________________________________

void SetCutsWindow::Save3DView()
{

     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += ".png";

     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     TEveViewer* view3d = ((TEveViewer*)*i);
     view3d->GetGLViewer()->SavePictureScale(file1,4.0);
     
     printf("Done.\n"); 

   return;

}

//______________________________________________________________________________

void SetCutsWindow::SaveRPhiView()
{

     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
     
     if (!file1.EndsWith(".png"))
       file1 += ".png";
     
     gSystem->ChangeDirectory(fi.fIniDir);
      
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     i++;
     TEveViewer* viewrphi = ((TEveViewer*)*i);
     viewrphi->GetGLViewer()->SavePictureScale(file1,4.0);
     
     printf("Done.\n"); 


   return;

}

//______________________________________________________________________________

void SetCutsWindow::SaveRhoZView()
{

     TGFileInfo fi;
     fi.fFileTypes   = gPictureSaveAsTypes;
     fi.fIniDir      = StrDup(""); // current directory
     fi.fFileTypeIdx = 0;
     fi.fOverwrite   = kTRUE;
     new TGFileDialog(gClient->GetDefaultRoot(),
     gEve->GetMainWindow(), kFDSave, &fi);
     if (!fi.fFilename) return;

     TPMERegexp filere(".*/([^/]+$)");
     if (filere.Match(fi.fFilename) != 2)
     {
       Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
       return;
     }

     TString file1(filere[1]);
       
     if (!file1.EndsWith(".png"))
       file1 += ".png";

     gSystem->ChangeDirectory(fi.fIniDir);
     
     printf("Saving...\n");

     TEveViewerList *viewers = gEve->GetViewers();
     TEveElement::List_i i = viewers->BeginChildren();
     i++;
     i++;
     i++;
     TEveViewer* viewrhoz = ((TEveViewer*)*i);
     viewrhoz->GetGLViewer()->SavePictureScale(file1,4.0);
     
     printf("Done.\n"); 


   return;

}
//______________________________________________________________________________

void SetCutsWindow::DrawPtHisto()
{
   if(gDrawHistograms[0])
      gDrawHistograms[0] = kFALSE;
   else
      gDrawHistograms[0] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawEtaHisto()
{
   if(gDrawHistograms[1])
      gDrawHistograms[1] = kFALSE;
   else
      gDrawHistograms[1] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPhiHisto()
{
   if(gDrawHistograms[2])
      gDrawHistograms[2] = kFALSE;
   else
      gDrawHistograms[2] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPhiPtHisto()
{
   if(gDrawHistograms[3])
      gDrawHistograms[3] = kFALSE;
   else
      gDrawHistograms[3] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPtYHisto()
{
   if(gDrawHistograms[4])
      gDrawHistograms[4] = kFALSE;
   else
      gDrawHistograms[4] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawEtaPhiHisto()
{
   if(gDrawHistograms[5])
      gDrawHistograms[5] = kFALSE;
   else
      gDrawHistograms[5] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPtHistoAll()
{
   if(gDrawHistograms[6])
      gDrawHistograms[6] = kFALSE;
   else
      gDrawHistograms[6] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawEtaHistoAll()
{
   if(gDrawHistograms[7])
      gDrawHistograms[7] = kFALSE;
   else
      gDrawHistograms[7] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPhiHistoAll()
{
   if(gDrawHistograms[8])
      gDrawHistograms[8] = kFALSE;
   else
      gDrawHistograms[8] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPhiPtHistoAll()
{
   if(gDrawHistograms[9])
      gDrawHistograms[9] = kFALSE;
   else
      gDrawHistograms[9] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawPtYHistoAll()
{
   if(gDrawHistograms[10])
      gDrawHistograms[10] = kFALSE;
   else
      gDrawHistograms[10] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawEtaPhiHistoAll()
{
   if(gDrawHistograms[11])
      gDrawHistograms[11] = kFALSE;
   else
      gDrawHistograms[11] = kTRUE;

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawHistos()
{

   Int_t nHistos1 = 0;
   Int_t nHistos2 = 0;

   TCanvas* pad1 = 0;
   TCanvas* pad2 = 0;

   if(gDrawHistograms[0])
   {
         nHistos1++;
         TH1D* histPt = new TH1D("Pt\nSingle Event", "AliEve Pt histogram", 100, 0.0, 5.0);
   }

   if(gDrawHistograms[1])
   {
         nHistos1++;
         TH1D* histEta = new TH1D("#eta\nSingle Event", "AliEve #eta histogram", 100, -1.5, 1.5);
   }

   if(gDrawHistograms[2])
   {
         nHistos1++;
         TH1D* histPhi = new TH1D("#phi\nSingle Event", "AliEve #phi histogram", 100, 0.0, 2*TMath::Pi());
   }

   if(gDrawHistograms[3])
   {
         nHistos2++;
         TH2D* histPhiPt = new TH2D("#phi-Pt\nSingle Event", "AliEve #phi-Pt histogram", 100, 0.0, 2*TMath::Pi(), 100, 0.0, 5.0);
   }

   if(gDrawHistograms[4])
   {
         nHistos2++;
         TH2D* histPtY = new TH2D("Pt-Y\nSingle Event", "AliEve Pt-y histogram", 100, 0.0, 5.0, 100, -1.5, 1.5);
   }

   if(gDrawHistograms[5])
   {
         nHistos2++;
         TH2D* histEtaPhi = new TH2D("#eta-#phi\nSingle Event", "AliEve #eta-#phi histogram", 100, -1.5, 1.5, 100, 0.0, 2*TMath::Pi());
   }

   AliESDEvent* esd = AliEveEventManager::AssertESD();

   if(esd->GetNumberOfTracks())
   {

      for(Int_t j = 0; j < esd->GetNumberOfTracks(); j++)
      {

         AliESDtrack* track = esd->GetTrack(j);

         if(gDrawHistograms[0])
            histPt->Fill(track->Pt());
         if(gDrawHistograms[1])
            histEta->Fill(track->Eta());
         if(gDrawHistograms[2])
            histPhi->Fill(track->Phi());
         if(gDrawHistograms[3])
            histPhiPt->Fill(track->Phi(),track->Pt());
         if(gDrawHistograms[4])
            histPtY->Fill(track->Pt(),track->Y());
         if(gDrawHistograms[5])
            histEtaPhi->Fill(track->Eta(),track->Phi());

      }

   }

   switch(nHistos1)
   {
      case 1:
         pad1 = new TCanvas("AliEve 1D histograms - Single Event","AliEve 1D histograms - Single Event",800,600);
         if(gDrawHistograms[0])
            histPt->Draw();
         if(gDrawHistograms[1])
            histEta->Draw();
         if(gDrawHistograms[2])
            histPhi->Draw();
         break;

      case 2:
         pad1 = new TCanvas("AliEve 1D histograms - Single Event","AliEve 1D histograms - Single Event",1200,500);
         pad1->Divide(2,1);
         if(!gDrawHistograms[0])
         {
            pad1->cd(1);
            histEta->Draw();
            pad1->cd(2);
            histPhi->Draw();
         }
         if(!gDrawHistograms[1])
         {
            pad1->cd(1);
            histPt->Draw();
            pad1->cd(2);
            histPhi->Draw();
         }
         if(!gDrawHistograms[2])
         {
            pad1->cd(1);
            histPt->Draw();
            pad1->cd(2);
            histEta->Draw();
         }
         break;

      case 3:
         pad1 = new TCanvas("AliEve 1D histograms - Single Event","AliEve 1D histograms - Single Event",1200,500);
         pad1->Divide(3,1);
         pad1->cd(1);
         histPt->Draw();
         pad1->cd(2);
         histEta->Draw();
         pad1->cd(3);
         histPhi->Draw();
         break;

      default:
         break;
   }

   switch(nHistos2)
   {
      case 1:
         pad2 = new TCanvas("AliEve 2D histograms - Single Event","AliEve 2D histograms - Single Event",800,600);
         if(gDrawHistograms[3])
            histPhiPt->Draw();
         if(gDrawHistograms[4])
            histPtY->Draw();
         if(gDrawHistograms[5])
            histEtaPhi->Draw();
         break;

      case 2:
         pad2 = new TCanvas("AliEve 2D histograms - Single Event","AliEve 2D histograms - Single Event",1200,500);
         pad2->Divide(2,1);
         if(!gDrawHistograms[3])
         {
            pad2->cd(1);
            histPtY->Draw();
            pad2->cd(2);
            histEtaPhi->Draw();
         }
         if(!gDrawHistograms[4])
         {
            pad2->cd(1);
            histPhiPt->Draw();
            pad2->cd(2);
            histEtaPhi->Draw();
         }
         if(!gDrawHistograms[5])
         {
            pad2->cd(1);
            histPhiPt->Draw();
            pad2->cd(2);
            histPtY->Draw();
         }
         break;

      case 3:
         pad2 = new TCanvas("AliEve 2D histograms - Single Event","AliEve 2D histograms - Single Event",1200,500);
         pad2->Divide(3,1);
         pad2->cd(1);
         histPhiPt->Draw();
         pad2->cd(2);
         histPtY->Draw();
         pad2->cd(3);
         histEtaPhi->Draw();
         break;

      default:
         break;
   }

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawHistosAll()
{

   Int_t nHistos1 = 0;
   Int_t nHistos2 = 0;

   TCanvas* pad1 = 0;
   TCanvas* pad2 = 0;

   if(gDrawHistograms[6])
   {
         nHistos1++;
         TH1D* histPt = new TH1D("Pt\nAll Events", "AliEve Pt histogram", 100, 0.0, 5.0);
   }

   if(gDrawHistograms[7])
   {
         nHistos1++;
         TH1D* histEta = new TH1D("#eta\nAll Events", "AliEve #eta histogram", 100, -1.5, 1.5);
   }

   if(gDrawHistograms[8])
   {
         nHistos1++;
         TH1D* histPhi = new TH1D("#phi\nAll Events", "AliEve #phi histogram", 100, 0.0, 2*TMath::Pi());
   }

   if(gDrawHistograms[9])
   {
         nHistos2++;
         TH2D* histPhiPt = new TH2D("#phi-Pt\nAll Events", "AliEve #phi-Pt histogram", 100, 0.0, 2*TMath::Pi(), 100, 0.0, 5.0);
   }

   if(gDrawHistograms[10])
   {
         nHistos2++;
         TH2D* histPtY = new TH2D("Pt-Y\nAll Events", "AliEve Pt-y histogram", 100, 0.0, 5.0, 100, -1.5, 1.5);
   }

   if(gDrawHistograms[11])
   {
         nHistos2++;
         TH2D* histEtaPhi = new TH2D("#eta-#phi\nAll Events", "AliEve #eta-#phi histogram", 100, 1.5, 1.5, 100, 0.0, 2*TMath::Pi());
   }

   Int_t nEvents = AliEveEventManager::GetMaster()->GetMaxEventId();

   AliEveEventManager::GetMaster()->GotoEvent(0);

   for(Int_t i = 0; i <= nEvents; i++)
   {

   AliESDEvent* esd = AliEveEventManager::AssertESD();

   if(esd->GetNumberOfTracks())
   {

      for(Int_t j = 0; j < esd->GetNumberOfTracks(); j++)
      {

         AliESDtrack* track = esd->GetTrack(j);

         if(gDrawHistograms[6])
            histPt->Fill(track->Pt());
         if(gDrawHistograms[7])
            histEta->Fill(track->Eta());
         if(gDrawHistograms[8])
            histPhi->Fill(track->Phi());
         if(gDrawHistograms[9])
            histPhiPt->Fill(track->Phi(),track->Pt());
         if(gDrawHistograms[10])
            histPtY->Fill(track->Pt(),track->Y());
         if(gDrawHistograms[11])
            histEtaPhi->Fill(track->Eta(),track->Phi());

      }

   }

      AliEveEventManager::GetMaster()->NextEvent();

   }

   switch(nHistos1)
   {
      case 1:
         pad1 = new TCanvas("AliEve 1D histograms - All Events","AliEve 1D histograms - All Events",800,600);
         if(gDrawHistograms[6])
            histPt->Draw();
         if(gDrawHistograms[7])
            histEta->Draw();
         if(gDrawHistograms[8])
            histPhi->Draw();
         break;

      case 2:
         pad1 = new TCanvas("AliEve 1D histograms - All Events","AliEve 1D histograms - All Events",1200,500);
         pad1->Divide(2,1);
         if(!gDrawHistograms[6])
         {
            pad1->cd(1);
            histEta->Draw();
            pad1->cd(2);
            histPhi->Draw();
         }
         if(!gDrawHistograms[7])
         {
            pad1->cd(1);
            histPt->Draw();
            pad1->cd(2);
            histPhi->Draw();
         }
         if(!gDrawHistograms[8])
         {
            pad1->cd(1);
            histPt->Draw();
            pad1->cd(2);
            histEta->Draw();
         }
         break;

      case 3:
         pad1 = new TCanvas("AliEve 1D histograms - All Events","AliEve 1D histograms - All Events",1200,500);
         pad1->Divide(3,1);
         pad1->cd(1);
         histPt->Draw();
         pad1->cd(2);
         histEta->Draw();
         pad1->cd(3);
         histPhi->Draw();
         break;

      default:
         break;
   }

   switch(nHistos2)
   {
      case 1:
         pad2 = new TCanvas("AliEve 2D histograms - All Events","AliEve 2D histograms - All Events",800,600);
         if(gDrawHistograms[9])
            histPt->Draw();
         if(gDrawHistograms[10])
            histEta->Draw();
         if(gDrawHistograms[11])
            histPhi->Draw();
         break;

      case 2:
         pad2 = new TCanvas("AliEve 2D histograms - All Events","AliEve 2D histograms - All Events",1200,500);
         pad2->Divide(2,1);
         if(!gDrawHistograms[9])
         {
            pad2->cd(1);
            histPtY->Draw();
            pad2->cd(2);
            histEtaPhi->Draw();
         }
         if(!gDrawHistograms[10])
         {
            pad2->cd(1);
            histPhiPt->Draw();
            pad2->cd(2);
            histEtaPhi->Draw();
         }
         if(!gDrawHistograms[11])
         {
            pad2->cd(1);
            histPhiPt->Draw();
            pad2->cd(2);
            histPtY->Draw();
         }
         break;

      case 3:
         pad2 = new TCanvas("AliEve 2D histograms - All Events","AliEve 2D histograms - All Events",1200,500);
         pad2->Divide(3,1);
         pad2->cd(1);
         histPhiPt->Draw();
         pad2->cd(2);
         histPtY->Draw();
         pad2->cd(3);
         histEtaPhi->Draw();
         break;

      default:
         break;
   }

   return;

}

//______________________________________________________________________________

void SetCutsWindow::DrawResiduals()
{

   TEveUtil::Macro("make_residuals.C");

   return;

}

//______________________________________________________________________________

Int_t SetCutsWindow::GetTrackColorByMomentum(Double_t momentum, Int_t size)
{


   Double_t step = 2.0/size;

   for(Int_t i = 0; i < size; i++)
   {

      if(momentum > i*step && momentum <= (i+1)*step)
      {
//         cout << "tutaj " << i <<endl;
         return i;
      }

   }

//   cout << i <<endl;
   return i;

}

//______________________________________________________________________________

void SetCutsWindow::SetStandardCuts()
{

   gDrawV0s->SetOn(kFALSE,kFALSE);
   gDrawCascades->SetOn(kFALSE,kFALSE);
   gDrawKinks->SetOn(kFALSE,kFALSE);
   gDrawVertex->SetOn(kFALSE,kFALSE);
   gDrawTracklets->SetOn(kFALSE,kFALSE);
   gDrawTracks->SetOn(kTRUE,kFALSE);
   gDrawClusters->SetOn(kFALSE,kFALSE);
   gDrawTracksType1->SetOn(kTRUE,kFALSE);
   gDrawTracksType2->SetOn(kTRUE,kFALSE);
   gDrawTracksType3->SetOn(kTRUE,kFALSE);
   gDrawTracksType4->SetOn(kFALSE,kFALSE);
   gDrawTracksType5->SetOn(kFALSE,kFALSE);
   gDrawTracksType6->SetOn(kFALSE,kFALSE);
   gDrawTracksType7->SetOn(kFALSE,kFALSE);
   gCutOnP->SetOn(kFALSE,kFALSE);
   gCutOnPt->SetOn(kTRUE,kFALSE);
   gCutOnEta->SetOn(kTRUE,kFALSE);
   gCutOnMult->SetOn(kFALSE,kFALSE);
   gCutOnCls->SetOn(kTRUE,kFALSE);
   gPtRange->SetValues(0.15,gPtRange->GetLimitMax());
   gEtaRange->SetValues(-0.9,0.9);
   gClsRangeNE->SetNumber(70);
   gClsRange->SetPosition(70);

}

//______________________________________________________________________________

void SetCutsWindow::AddMomentumVectors()
{

//   Bool_t drawWithTracks = kTRUE;

   if(gEve->GetEventScene()->FirstChild()->FindChild("Momentum Vectors"))
      gEve->GetEventScene()->FirstChild()->FindChild("Momentum Vectors")->Destroy();

   TString str1 = 0;
   TString str2 = 0;

   TEveElement::List_i i = gEve->GetEventScene()->FirstChild()->BeginChildren();
   TEveElement::List_i j = gEve->GetEventScene()->FirstChild()->EndChildren();
   TEveElement::List_i k;

   TEveElementList* momentumVectorList = new TEveElementList("Momentum Vectors");

   Double_t maxMomentum = 0;
   Double_t vectorLength = 600.0;

   Double_t x1 = 0;
   Double_t y1 = 0;
   Double_t z1 = 0;

   Double_t x2 = 0;
   Double_t y2 = 0;
   Double_t z2 = 0;

   Bool_t draw = kFALSE;

   Int_t mode = gVectorMode->GetSelected();
   Double_t cut = 0.1*gPMVRange->GetPosition();

   //==============================================
   // find highest momentum (to normalize)
   //==============================================

   if(mode != 3)
   {
   for(k = i; k != j; k++)
   {
      TEveElement* element = (TEveElement*) *k;

      str1 = element->GetElementName();

      if(str1.Contains("Tracks") || str1.Contains("tracks"))
      {

         TEveElement::List_i m = element->BeginChildren();
         TEveElement::List_i n = element->EndChildren();
         TEveElement::List_i l;

         for(l = m; l != n; l++)
         {
            TEveElement* trackType = (TEveElement*) *l;
            str2 = trackType->GetElementName();

            if(str2.Contains("Sigma < 3"))
            {
               if(trackType->HasChildren())
               {

                  TEveElement::List_i x = trackType->BeginChildren();
                  TEveElement::List_i y = trackType->EndChildren();
                  TEveElement::List_i z;

                  for(z = x; z != y; z++)
                  {

                     AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                     if(trackSingle1->GetESDTrack()->P() > maxMomentum)
                       maxMomentum = trackSingle1->GetESDTrack()->P();

                  }
               }
            }


            if(str2.Contains("3 < Sigma < 5"))
            {

               if(trackType->HasChildren())
               {

                  TEveElement::List_i x = trackType->BeginChildren();
                  TEveElement::List_i y = trackType->EndChildren();
                  TEveElement::List_i z;

                  for(z = x; z != y; z++)
                  {

                     AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                     if(trackSingle1->GetESDTrack()->P() > maxMomentum)
                       maxMomentum = trackSingle1->GetESDTrack()->P();

                  }
               }
            }

            if(str2.Contains("5 < Sigma"))
            {

               if(trackType->HasChildren())
               {

                  TEveElement::List_i x = trackType->BeginChildren();
                  TEveElement::List_i y = trackType->EndChildren();
                  TEveElement::List_i z;

                  for(z = x; z != y; z++)
                  {

                     AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                     if(trackSingle1->GetESDTrack()->P() > maxMomentum)
                       maxMomentum = trackSingle1->GetESDTrack()->P();


                  }
               }
            }
         }
      }
   }
   }
   //==============================================
   // clean the display
   //==============================================
/*
   if(!drawWithTracks)

      for(k = i; k != j; k++)
      {
         TEveElement* element = (TEveElement*) *k;

         str1 = element->GetElementName();

         element->SetRnrSelf(kFALSE);

         if(element->HasChildren())
            element->SetRnrChildren(kFALSE);

      }

   }
*/
   //==============================================
   // draw momentum vectors
   //==============================================

   if(mode == 1 && maxMomentum > 1)
      vectorLength = 100/TMath::Log(100*maxMomentum);
   if(mode == 2 && maxMomentum)
      vectorLength = vectorLength/maxMomentum;
   if(mode == 3)
      vectorLength = 100;

   for(k = i; k != j; k++)
   {
      TEveElement* element = (TEveElement*) *k;

      str1 = element->GetElementName();

      if(str1.Contains("Tracks") || str1.Contains("tracks"))
      {

         TEveElement::List_i m = element->BeginChildren();
         TEveElement::List_i n = element->EndChildren();
         TEveElement::List_i l;

         for(l = m; l != n; l++)
         {
            TEveElement* trackType = (TEveElement*) *l;
            str2 = trackType->GetElementName();

//            trackType->SetRnrSelf(kFALSE);

//            if(trackType->HasChildren())
//               trackType->SetRnrChildren(kFALSE);

            if(str2.Contains("Sigma < 3"))
            {

               if(trackType->HasChildren())
               {

                  TEveElementList* momentumVectorList1 = new TEveElementList("sigma < 3");

                  TEveElement::List_i x = trackType->BeginChildren();
                  TEveElement::List_i y = trackType->EndChildren();
                  TEveElement::List_i z;

                  for(z = x; z != y; z++)
                  {

                     AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                     TEveLine* momentumVector = new TEveLine(TString::Format("Momentum Vector"));

                     if(trackSingle1->GetESDTrack()->P() > cut)
                     {

                        x1 = trackSingle1->GetESDTrack()->Xv();
                        y1 = trackSingle1->GetESDTrack()->Yv();
                        z1 = trackSingle1->GetESDTrack()->Zv();

                        momentumVector->SetPoint(0, x1, y1, z1);

                        if(mode == 1)
                        {
                           x2 = x1+vectorLength*TMath::Log(100*trackSingle1->GetESDTrack()->P())*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*TMath::Log(100*trackSingle1->GetESDTrack()->P())*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*TMath::Log(100*trackSingle1->GetESDTrack()->P())*trackSingle1->GetESDTrack()->Pz();
                        }

                        if(mode == 2)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                        }

                        if(mode == 3)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                        }

                        if((mode != 1 && mode!= 2 && mode != 3) ||
                          ( mode == 3 && trackSingle1->GetESDTrack()->P() <= 0.01))
                           continue;

                        momentumVector->SetPoint(1, x2, y2, z2);

/*
                        if(trackSingle1->GetESDTrack()->Charge() == -1)
                           momentumVector->SetLineColor(kGreen);
                        else
                           momentumVector->SetLineColor(kRed);
*/
                        momentumVector->SetLineColor(kRed);

                        momentumVector->SetLineWidth(1);
                        momentumVector->SetLineStyle(0);
                        momentumVector->SetTitle(Form("%f GeV/c", trackSingle1->GetESDTrack()->P()));

                        momentumVectorList1->AddElement(momentumVector);

                     }

                  }

                  momentumVectorList->AddElement(momentumVectorList1);

                  draw = kTRUE;

               }
            }


            if(str2.Contains("3 < Sigma < 5"))
            {

               if(trackType->HasChildren())
               {

                  TEveElement::List_i x = trackType->BeginChildren();
                  TEveElement::List_i y = trackType->EndChildren();
                  TEveElement::List_i z;

                  TEveElementList* momentumVectorList2 = new TEveElementList("3 < sigma < 5");

                  for(z = x; z != y; z++)
                  {

                     AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                     TEveLine* momentumVector = new TEveLine(TString::Format("Momentum Vector"));

                     if(trackSingle1->GetESDTrack()->P() > cut)
                     {

                        x1 = trackSingle1->GetESDTrack()->Xv();
                        y1 = trackSingle1->GetESDTrack()->Yv();
                        z1 = trackSingle1->GetESDTrack()->Zv();

                        momentumVector->SetPoint(0, x1, y1, z1);

                        if(mode == 1)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px()*TMath::Log(trackSingle1->GetESDTrack()->P());
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py()*TMath::Log(trackSingle1->GetESDTrack()->P());
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz()*TMath::Log(trackSingle1->GetESDTrack()->P());
                        }

                        if(mode == 2)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                        }

                        if(mode == 3)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                        }

                        if(mode != 1 && mode!= 2 && mode != 3)
                           continue;

                        momentumVector->SetPoint(1, x2, y2, z2);

/*
                        if(trackSingle1->GetESDTrack()->Charge() == -1)
                           momentumVector->SetLineColor(kGreen+2);
                        else
                           momentumVector->SetLineColor(kRed+2);
*/
                        momentumVector->SetLineColor(kRed+2);

                        momentumVector->SetLineWidth(1);
                        momentumVector->SetLineStyle(0);
                        momentumVector->SetTitle(Form("%f GeV/c", trackSingle1->GetESDTrack()->P()));

                        momentumVectorList2->AddElement(momentumVector);

                     }

                  }

//                  gEve->AddElement(momentumVectorList2);
                  momentumVectorList->AddElement(momentumVectorList2);

                  draw = kTRUE;

               }
            }

            if(str2.Contains("5 < Sigma"))
            {

               if(trackType->HasChildren())
               {

                  TEveElementList* momentumVectorList3 = new TEveElementList("5 < sigma");

                  TEveElement::List_i x = trackType->BeginChildren();
                  TEveElement::List_i y = trackType->EndChildren();
                  TEveElement::List_i z;

                  for(z = x; z != y; z++)
                  {

                     AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                     TEveLine* momentumVector = new TEveLine(TString::Format("Momentum Vector"));

                     if(trackSingle1->GetESDTrack()->P() > cut)
                     {

                        x1 = trackSingle1->GetESDTrack()->Xv();
                        y1 = trackSingle1->GetESDTrack()->Yv();
                        z1 = trackSingle1->GetESDTrack()->Zv();

                        momentumVector->SetPoint(0, x1, y1, z1);

                        if(mode == 1)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px()*TMath::Log(trackSingle1->GetESDTrack()->P());
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py()*TMath::Log(trackSingle1->GetESDTrack()->P());
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz()*TMath::Log(trackSingle1->GetESDTrack()->P());
                        }

                        if(mode == 2)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                        }

                        if(mode == 3)
                        {
                           x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                           y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                           z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                        }

                        if(mode != 1 && mode!= 2 && mode != 3)
                           continue;

                        momentumVector->SetPoint(1, x2, y2, z2);
/*
                        if(trackSingle1->GetESDTrack()->Charge() == -1)
                           momentumVector->SetLineColor(kGreen+3);
                        else
                           momentumVector->SetLineColor(kRed+3);
*/
                        momentumVector->SetLineColor(kRed+3);

                        momentumVector->SetLineWidth(1);
                        momentumVector->SetLineStyle(0);
                        momentumVector->SetTitle(Form("%f GeV/c", trackSingle1->GetESDTrack()->P()));

                        momentumVectorList3->AddElement(momentumVector);

                     }

                  }

                  //gEve->AddElement(momentumVectorList3);
                  momentumVectorList->AddElement(momentumVectorList3);

                  draw = kTRUE;

               }
            }
         }
      }
   }
 
  gEve->AddElement(momentumVectorList);

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *mv = AliEveMultiView::Instance();

  mv->DestroyEventRPhi();
  mv->DestroyEventRhoZ();

  mv->ImportEventRPhi(top);
  mv->ImportEventRhoZ(top);

  gEve->FullRedraw3D(kFALSE, kTRUE);

}

//______________________________________________________________________________

void SetCutsWindow::SetCuts()
{

   TString str1 = 0;
   TString str2 = 0;

   Int_t posTrackColor= gPosColorList->GetSelected();
   Int_t negTrackColor= gNegColorList->GetSelected();

   TEveElement::List_i i = gEve->GetEventScene()->FirstChild()->BeginChildren();
   TEveElement::List_i j = gEve->GetEventScene()->FirstChild()->EndChildren();
   TEveElement::List_i k;

   Double_t x1, y1, z1, x2, y2, z2;

   for(k = i; k != j; k++)
   {
      TEveElement* element = (TEveElement*) *k;

      str1 = element->GetElementName();

      if(gDrawV0s->IsOn())
      {
         if(str1.Contains("ESD v0") || str1.Contains("ESD V0"))
         {
            element->SetRnrSelf(kTRUE);

            if(element->HasChildren())
               element->SetRnrChildren(kTRUE);

            continue;
         }
      }

      if(gDrawCascades->IsOn())
      {
         if(str1.Contains("ESD cascade") || str1.Contains("ESD Cascade"))
         {
            element->SetRnrSelf(kTRUE);

            if(element->HasChildren())
               element->SetRnrChildren(kTRUE);

            continue;
         }
      }

      if(gDrawKinks->IsOn())
      {
         if(str1.Contains("ESD kink") || str1.Contains("ESD Kink"))
         {
            element->SetRnrSelf(kTRUE);

            if(element->HasChildren())
               element->SetRnrChildren(kTRUE);

            continue;
         }
      }

      if(gDrawVertex->IsOn())
      {
         if(str1.Contains("Primary Vertex") || str1.Contains("primary Vertex") || str1.Contains("Primary vertex") || str1.Contains("primary vertex"))
         {
            element->SetRnrSelf(kTRUE);

            if(element->HasChildren())
               element->SetRnrChildren(kTRUE);

            continue;
         }
      }

      if(gDrawTracklets->IsOn())
      {
         if(str1.Contains("Tracklets") || str1.Contains("tracklets"))
         {
            element->SetRnrSelf(kTRUE);

            if(element->HasChildren())
               element->SetRnrChildren(kTRUE);

            continue;
         }
      }

      if(gDrawTracks->IsOn())
      {

//         Int_t colorPos[10] = {kRed, kBlue, kOrange, kCyan, kGreen, kGray, kViolet, kMagenta, kSpring, kYellow};

Int_t colorNeg[10] = 
{
kCyan-4,
kCyan,
kAzure+10,
kAzure+8,
kAzure+5,
kAzure,
kBlue,
kBlue+1,
kBlue+2,
kBlue+3
};

Int_t colorPos[10] = 
{
kYellow-4,
kYellow,
kOrange+10,
kOrange,
kOrange+7,
kOrange+10,
kRed,
kRed+1,
kRed+2,
kRed+3
};

Int_t colorAll[22] = 
{
kBlue+4,
kBlue+2,
kBlue,
kAzure,
kAzure-3,
kAzure+7,
kCyan,
kCyan-7,
kGreen-7,
kGreen-4,
kGreen,
kSpring,
kSpring+7,
kSpring+8,
kYellow,
kOrange,
kOrange-3,
kOrange+7,
kOrange+4,
kRed,
kRed+2,
kMagenta
};

         if(str1.Contains("Tracks") || str1.Contains("tracks"))
         {
            element->SetRnrSelf(kTRUE);

            if(element->HasChildren())
               element->SetRnrChildren(kTRUE);

            TEveElement::List_i m = element->BeginChildren();
            TEveElement::List_i n = element->EndChildren();
            TEveElement::List_i l;

            for(l = m; l != n; l++)
            {
               TEveElement* trackType = (TEveElement*) *l;
               str2 = trackType->GetElementName();

               trackType->SetRnrSelf(kFALSE);

               (dynamic_cast<TEveTrackList*>trackType)->GetPropagator()->SetMaxR(250);

               if(trackType->HasChildren())
                  trackType->SetRnrChildren(kFALSE);

               if(gDrawTracksType1->IsOn() && str2.Contains("Sigma < 3"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {

                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {

                        AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                        if(trackSingle1->GetESDTrack()->GetSign() > 0)
                        {
                           if(posTrackColor == 1)
                              trackSingle1->SetLineColor(colorNeg[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(posTrackColor == 2)
                              trackSingle1->SetLineColor(colorPos[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(posTrackColor == 3)
                              trackSingle1->SetLineColor(colorAll[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),22)]);
                        }
                        else
                        {
                           if(negTrackColor == 1)
                              trackSingle1->SetLineColor(colorNeg[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(negTrackColor == 2)
                              trackSingle1->SetLineColor(colorPos[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(negTrackColor == 3)
                              trackSingle1->SetLineColor(colorAll[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),22)]);
                        }

                        trackSingle1->SetLineStyle(1);
                        trackSingle1->SetRnrSelf(kTRUE);


                     }

                  }

               }


               if(gDrawTracksType2->IsOn() && str2.Contains("3 < Sigma < 5"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {
                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {

                        AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                        if(trackSingle1->GetESDTrack()->GetSign() > 0)
                        {
                           if(posTrackColor == 1)
                              trackSingle1->SetLineColor(colorNeg[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(posTrackColor == 2)
                              trackSingle1->SetLineColor(colorPos[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(posTrackColor == 3)
                              trackSingle1->SetLineColor(colorAll[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),22)]);
                        }
                        else
                        {
                           if(negTrackColor == 1)
                              trackSingle1->SetLineColor(colorNeg[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(negTrackColor == 2)
                              trackSingle1->SetLineColor(colorPos[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(negTrackColor == 3)
                              trackSingle1->SetLineColor(colorAll[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),22)]);
                        }

                        trackSingle1->SetLineStyle(1);
                        trackSingle1->SetRnrSelf(kTRUE);

                     }
                  }
               }

               if(gDrawTracksType3->IsOn() && str2.Contains("5 < Sigma"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {

                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {

                        AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);

                        if(trackSingle1->GetESDTrack()->GetSign() > 0)
                        {
                           if(posTrackColor == 1)
                              trackSingle1->SetLineColor(colorNeg[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(posTrackColor == 2)
                              trackSingle1->SetLineColor(colorPos[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(posTrackColor == 3)
                              trackSingle1->SetLineColor(colorAll[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),22)]);
                        }
                        else
                        {
                           if(negTrackColor == 1)
                              trackSingle1->SetLineColor(colorNeg[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(negTrackColor == 2)
                              trackSingle1->SetLineColor(colorPos[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),10)]);
                           if(negTrackColor == 3)
                              trackSingle1->SetLineColor(colorAll[GetTrackColorByMomentum(trackSingle1->GetESDTrack()->Pt(),22)]);
                        }

                        trackSingle1->SetLineStyle(1);
                        trackSingle1->SetRnrSelf(kTRUE);

                     }
                  }
               }

               if(gDrawTracksType4->IsOn() && str2.Contains("no ITS refit"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {

                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {
                        TEveElement* trackSingle = (TEveElement*) *z;
                        trackSingle->SetRnrSelf(kTRUE);
                     }
                  }
               }

               if(gDrawTracksType5->IsOn() && str2.Contains("no TPC refit"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {

                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {
                        TEveElement* trackSingle = (TEveElement*) *z;
                        trackSingle->SetRnrSelf(kTRUE);
                     }
                  }
               }

               if(gDrawTracksType6->IsOn() && str2.Contains("ITS ncl>=3"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {

                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {
                        TEveElement* trackSingle = (TEveElement*) *z;
                        trackSingle->SetRnrSelf(kTRUE);
                     }
                  }
               }

               if(gDrawTracksType7->IsOn() && str2.Contains("ITS others"))
               {
                  trackType->SetRnrSelf(kTRUE);

                  if(trackType->HasChildren())
                  {

                     trackType->SetRnrChildren(kTRUE);
                     TEveElement::List_i x = trackType->BeginChildren();
                     TEveElement::List_i y = trackType->EndChildren();
                     TEveElement::List_i z;

                     for(z = x; z != y; z++)
                     {
                        TEveElement* trackSingle = (TEveElement*) *z;
                        trackSingle->SetRnrSelf(kTRUE);
                     }
                  }
               }

            }
            
            continue;
         }

      }

      if(gDrawClusters->IsOn())
      {

         TEvePointSetArray * cc = 0;
         TEvePointSet* clusters = 0;

         if((str1.Contains("Clusters") && !str1.Contains("TPC")) || str1.Contains("Colorized"))
         {

            if(!gCutOnEta->IsOn())
            {            
               element->SetRnrSelf(kTRUE);
               element->SetRnrChildren(kTRUE);

                  if(gEve->GetEventScene()->FirstChild()->FindChild("ITS ClCut"))
                     gEve->GetEventScene()->FirstChild()->FindChild("ITS ClCut")->Destroy();

               continue;
            }
            else
            {            
               element->SetRnrSelf(kFALSE);
               element->SetRnrChildren(kFALSE);

               if(str1.Contains("ITS"))
               {

                  clusters = dynamic_cast<TEvePointSet*>(element);

                  if(gEve->GetEventScene()->FirstChild()->FindChild("ITS ClCut"))
                     gEve->GetEventScene()->FirstChild()->FindChild("ITS ClCut")->Destroy();

                  pointset = new TEvePointSet(clusters->GetLastPoint());
                  pointset->SetMarkerStyle(4);
                  pointset->SetMarkerColor(kBlue);
                  pointset->SetMarkerSize(0.4);
                  pointset->SetName("ITS ClCut");

                  for(Int_t iCluster = 0; iCluster < clusters->GetLastPoint(); iCluster++)
                  {

                     clusters->GetPoint(iCluster, x1, y1, z1);

                     if(TMath::Sqrt(x1*x1 + y1*y1 + z1*z1) != 0 && z1 != 0)
                     {
                        Double_t eta = -TMath::Log(TMath::Tan(0.5*TMath::ACos(z1/TMath::Sqrt(x1*x1 + y1*y1 + z1*z1))));

                        if(eta > gEtaRange->GetMin() && eta < gEtaRange->GetMax())
                        {
                           pointset->SetNextPoint(x1, y1, z1);
                        }
                     }
                  }

                  pointset->SetRnrSelf(kTRUE);
                  pointset->SetRnrChildren(kTRUE);

                  gEve->AddElement(pointset);

               }

               if(str1.Contains("TPC"))
               {

                  cc = new TEvePointSetArray("TPC ClCut");
                  cc->SetMainColor(kRed);
                  cc->SetMarkerStyle(4);
                  cc->SetMarkerSize(0.4);
                  cc->InitBins("Cluster Charge",
                               (dynamic_cast<TEvePointSetArray*>(element))->GetNBins()-2,
                               (dynamic_cast<TEvePointSetArray*>(element))->GetMin(),
                               (dynamic_cast<TEvePointSetArray*>(element))->GetMax());
  
                  cc->GetBin(0)->SetMainColor(kGray);
                  cc->GetBin(0)->SetMarkerSize(0.4);
                  cc->GetBin(1)->SetMainColor(kBlue);
                  cc->GetBin(1)->SetMarkerSize(0.42);
                  cc->GetBin(2)->SetMainColor(kCyan);
                  cc->GetBin(2)->SetMarkerSize(0.44);
                  cc->GetBin(3)->SetMainColor(kGreen);
                  cc->GetBin(3)->SetMarkerSize(0.46);
                  cc->GetBin(4)->SetMainColor(kYellow);
                  cc->GetBin(4)->SetMarkerSize(0.48);
                  cc->GetBin(5)->SetMainColor(kRed);
                  cc->GetBin(5)->SetMarkerSize(0.50);
                  cc->GetBin(6)->SetMainColor(kMagenta);
                  cc->GetBin(6)->SetMarkerSize(0.52);

                  Double_t range = (cc->GetMax()) - (cc->GetMin());

                  if(gEve->GetEventScene()->FirstChild()->FindChild("TPC ClCut"))
                     gEve->GetEventScene()->FirstChild()->FindChild("TPC ClCut")->Destroy();

                  for(Int_t iBin = 0; iBin < cc->GetNBins(); iBin++)
                  {

                     clusters =(dynamic_cast<TEvePointSetArray*>(element))->GetBin(iBin);

                     for(Int_t iCluster = 0; iCluster < clusters->GetLastPoint(); iCluster++)
                     {

                        clusters->GetPoint(iCluster, x1, y1, z1);

                        if(TMath::Sqrt(x1*x1 + y1*y1 + z1*z1) != 0 && z1 != 0)
                        {
                           Double_t eta = -TMath::Log(TMath::Tan(0.5*TMath::ACos(z1/TMath::Sqrt(x1*x1 + y1*y1 + z1*z1))));

                           if(eta > gEtaRange->GetMin() && eta < gEtaRange->GetMax())
                           {
                              cc->Fill(x1, y1, z1,(range/(cc->GetNBins())*iBin)-1);
                           }
                        }
                     }

                  }

                  cc->SetRnrSelf(kTRUE);
                  cc->SetRnrChildren(kTRUE);

                  gEve->AddElement(cc);

               }

               if(str1.Contains("TRD"))
               {

                  clusters = dynamic_cast<TEvePointSet*>(element);

                  if(gEve->GetEventScene()->FirstChild()->FindChild("TRD ClCut"))
                     gEve->GetEventScene()->FirstChild()->FindChild("TRD ClCut")->Destroy();

                  pointset = new TEvePointSet(clusters->GetLastPoint());
                  pointset->SetMarkerStyle(4);
                  pointset->SetMarkerColor(kCyan);
                  pointset->SetMarkerSize(0.4);
                  pointset->SetName("TRD ClCut");

                  for(Int_t iCluster = 0; iCluster < clusters->GetLastPoint(); iCluster++)
                  {

                     clusters->GetPoint(iCluster, x1, y1, z1);

                     if(TMath::Sqrt(x1*x1 + y1*y1 + z1*z1) != 0 && z1 != 0)
                     {
                        Double_t eta = -TMath::Log(TMath::Tan(0.5*TMath::ACos(z1/TMath::Sqrt(x1*x1 + y1*y1 + z1*z1))));

                        if(eta > gEtaRange->GetMin() && eta < gEtaRange->GetMax())
                        {
                           pointset->SetNextPoint(x1, y1, z1);
                        }
                     }
                  }

                  pointset->SetRnrSelf(kTRUE);
                  pointset->SetRnrChildren(kTRUE);

                  gEve->AddElement(pointset);

                  }

                  if(str1.Contains("TOF"))
                  {

                     clusters = dynamic_cast<TEvePointSet*>(element);

                     if(gEve->GetEventScene()->FirstChild()->FindChild("TOF ClCut"))
                        gEve->GetEventScene()->FirstChild()->FindChild("TOF ClCut")->Destroy();

                     pointset = new TEvePointSet(clusters->GetLastPoint());
                     pointset->SetMarkerStyle(4);
                     pointset->SetMarkerColor(kOrange+9);
                     pointset->SetMarkerSize(0.4);
                     pointset->SetName("TOF ClCut");

                     for(Int_t iCluster = 0; iCluster < clusters->GetLastPoint(); iCluster++)
                     {

                        clusters->GetPoint(iCluster, x1, y1, z1);

                        if(TMath::Sqrt(x1*x1 + y1*y1 + z1*z1) != 0 && z1 != 0)
                        {
                           Double_t eta = -TMath::Log(TMath::Tan(0.5*TMath::ACos(z1/TMath::Sqrt(x1*x1 + y1*y1 + z1*z1))));

                           if(eta > gEtaRange->GetMin() && eta < gEtaRange->GetMax())
                           {
                              pointset->SetNextPoint(x1, y1, z1);
                           }
                        }
                     }

                     pointset->SetRnrSelf(kTRUE);
                     pointset->SetRnrChildren(kTRUE);

                     gEve->AddElement(pointset);

                  }

               continue;

            }

         }
/*
         if(str1.Contains("Colorized"))
         {

           cout << "\n\n\n" << (dynamic_cast<TEvePointSetArray*>(element))->GetNBins() << "\n\n\n" << endl;

            if(!gCutOnEta->IsOn())
            {            
               element->SetRnrSelf(kTRUE);
               element->SetRnrChildren(kTRUE);
            }
            else
            {            
               element->SetRnrSelf(kFALSE);
               element->SetRnrChildren(kFALSE);
            }

            continue;

         }
*/
      }
/*
      if(str1.Contains("TPC") && str1.Contains("Clusters") && !str1.Contains("Colorized"))
      {

         element->SetRnrChildren(kFALSE);
         element->SetRnrSelf(kFALSE);

         if(gEve->GetEventScene()->FirstChild()->FindChild("TPC ClCut"))
            gEve->GetEventScene()->FirstChild()->FindChild("TPC ClCut")->Destroy();

         if(gCutOnEta->IsOn())
         {
            clusters = dynamic_cast<TEvePointSet*>(element);

            pointset = new TEvePointSet(clusters->GetLastPoint());
            pointset->SetMarkerStyle(4);
            pointset->SetMarkerColor(kBlue);
            pointset->SetMarkerSize(0.4);
            pointset->SetName("TPC ClCut");

            for(Int_t iCluster = 0; iCluster < clusters->GetLastPoint(); iCluster++)
            {

               clusters->GetPoint(iCluster, x1, y1, z1);

               if(TMath::Sqrt(x1*x1 + y1*y1 + z1*z1) != 0 && z1 != 0)
               {
                  Double_t eta = -TMath::Log(TMath::Tan(0.5*TMath::ACos(z1/TMath::Sqrt(x1*x1 + y1*y1 + z1*z1))));

                  if(eta > gEtaRange->GetMin() && eta < gEtaRange->GetMax())
                  {
                     pointset->SetNextPoint(x1, y1, z1);
                  }
               }
            }

            pointset->SetRnrSelf(kTRUE);
            pointset->SetRnrChildren(kTRUE);

            gEve->AddElement(pointset);

         }

         continue;

      }
*/
      if(!str1.Contains("ClCut"))
      {
         element->SetRnrChildren(kFALSE);
         element->SetRnrSelf(kFALSE);
      }

   }

   if(gDrawTracks->IsOn() || gDrawV0s->IsOn() || gDrawCascades->IsOn() || gDrawKinks->IsOn())
   {

      i = gEve->GetEventScene()->FirstChild()->FindChild("ESD Tracks by category")->BeginChildren();
      j = gEve->GetEventScene()->FirstChild()->FindChild("ESD Tracks by category")->EndChildren();

      for(k = i; k != j; k++)
      {

         TEveElement* trackList = (TEveElement*) *k;

         TEveElement::List_i m = trackList->BeginChildren();
         TEveElement::List_i n = trackList->EndChildren();
         TEveElement::List_i l;

         for(l = m; l != n; l++)
         {

            AliEveTrack *track = dynamic_cast<AliEveTrack*>((TEveElement*) *l);

            if(gCutOnMult->IsOn())
            {

               Double_t draw = gRandom->Rndm();

               if(draw > (gMultRangeNE->GetNumber())/100)
               {

                  track->SetRnrSelf(kFALSE);
                  continue;

               }
            }

            if(gCutOnCls->IsOn())
            {

               if(track->GetESDTrack()->GetTPCNcls() < gClsRangeNE->GetNumber())
               {
                     track->SetRnrSelf(kFALSE);
                     continue;
               }

            }

            if(gCutOnP->IsOn())
            {

               if(gPRange->GetMax() == gPRange->GetLimitMax())
               {

                  if(track->GetESDTrack()->P() < gPRange->GetMin())
                  {
                     track->SetRnrSelf(kFALSE);
                     continue;
                  }

               }

               else
               {
                  if(track->GetESDTrack()->P() < gPRange->GetMin() || track->GetESDTrack()->P() > gPRange->GetMax())
                  {
                     track->SetRnrSelf(kFALSE);
                     continue;
                  }

               }

            }

            if(gCutOnPt->IsOn())
            {

               if(gPtRange->GetMax() == gPtRange->GetLimitMax())
               {

                  if(track->GetESDTrack()->Pt() < gPtRange->GetMin())
                  {
                     track->SetRnrSelf(kFALSE);
                     continue;
                  }

               }

               else
               {
                  if(track->GetESDTrack()->Pt() < gPtRange->GetMin() || track->GetESDTrack()->Pt() > gPtRange->GetMax())
                  {
                     track->SetRnrSelf(kFALSE);
                     continue;
                  }

               }

            }

            if(gCutOnEta->IsOn())
            {

               if(gEtaRange->GetMax() == gEtaRange->GetLimitMax() && gEtaRange->GetMin() == gEtaRange->GetLimitMin())
                  continue;

               if(gEtaRange->GetMax() == gEtaRange->GetLimitMax())
               {

                  if(track->GetESDTrack()->Eta() < gEtaRange->GetMin())
                  {
                     track->SetRnrSelf(kFALSE);
                     continue;
                  }

               }

               if(gEtaRange->GetMin() == gEtaRange->GetLimitMin())
               {

                  if(track->GetESDTrack()->Eta() > gEtaRange->GetMax())
                  {
                     track->SetRnrSelf(kFALSE);
                     continue;
                  }

               }

               if(track->GetESDTrack()->Eta() < gEtaRange->GetMin() || track->GetESDTrack()->Eta() > gEtaRange->GetMax())
               {
                  track->SetRnrSelf(kFALSE);
                  continue;
               }

            }

         }   

      }

/*
      i = gEve->GetEventScene()->FirstChild()->FindChild("ESD v0")->BeginChildren();
      j = gEve->GetEventScene()->FirstChild()->FindChild("ESD v0")->EndChildren();

      for(k = i; k != j; k++)
      {

         AliEveV0 *v0 = dynamic_cast<AliEveV0*>((TEveElement*) *k);
            
         if(gCutOnP->IsOn())
         {

            if(gPRange->GetMax() == gPRange->GetLimitMax())
            {

               if(v0->GetP() < gPRange->GetMin())
               {
                  v0->SetRnrSelf(kFALSE);
                  v0->SetRnrChildren(kFALSE);
                  v0->GetPosTrack()->SetRnrSelf(kFALSE);
                  v0->GetNegTrack()->SetRnrSelf(kFALSE);
                  v0->GetPointingLine()->SetRnrSelf(kFALSE);
                  v0->GetPosTrack()->SetRnrChildren(kFALSE);
                  v0->GetNegTrack()->SetRnrChildren(kFALSE);
                  v0->GetPointingLine()->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            else
            {
               if(v0->GetP() < gPRange->GetMin() || v0->GetP() > gPRange->GetMax())
               {
                  v0->SetRnrSelf(kFALSE);
                  v0->SetRnrChildren(kFALSE);
                  v0->GetPosTrack()->SetRnrSelf(kFALSE);
                  v0->GetNegTrack()->SetRnrSelf(kFALSE);
                  v0->GetPointingLine()->SetRnrSelf(kFALSE);
                  v0->GetPosTrack()->SetRnrChildren(kFALSE);
                  v0->GetNegTrack()->SetRnrChildren(kFALSE);
                  v0->GetPointingLine()->SetRnrChildren(kFALSE);
                  continue;
               }

            }

         }

         if(gCutOnPt->IsOn())
         {

            if(gPtRange->GetMax() == gPtRange->GetLimitMax())
            {

               if(v0->GetPt() < gPtRange->GetMin())
               {
                  v0->SetRnrSelf(kFALSE);
                  v0->SetRnrChildren(kFALSE);
                  v0->GetPosTrack()->SetRnrSelf(kFALSE);
                  v0->GetNegTrack()->SetRnrSelf(kFALSE);
                  v0->GetPointingLine()->SetRnrSelf(kFALSE);
                  v0->GetPosTrack()->SetRnrChildren(kFALSE);
                  v0->GetNegTrack()->SetRnrChildren(kFALSE);
                  v0->GetPointingLine()->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            else
            {
               if(v0->GetPt() < gPtRange->GetMin() || v0->GetPt() > gPtRange->GetMax())
               {
                  v0->SetRnrSelf(kFALSE);
                  v0->SetRnrChildren(kFALSE);
                  v0->GetPosTrack()->SetRnrSelf(kFALSE);
                  v0->GetNegTrack()->SetRnrSelf(kFALSE);
                  v0->GetPointingLine()->SetRnrSelf(kFALSE);
                  v0->GetPosTrack()->SetRnrChildren(kFALSE);
                  v0->GetNegTrack()->SetRnrChildren(kFALSE);
                  v0->GetPointingLine()->SetRnrChildren(kFALSE);
                  continue;
               }

            }

         }

         if(gCutOnEta->IsOn())
         {

            if(gEtaRange->GetMax() == gEtaRange->GetLimitMax() && gEtaRange->GetMax() == gEtaRange->GetLimitMax())
               continue;

            if(gEtaRange->GetMax() == gEtaRange->GetLimitMax())
            {
               if(v0->GetEta() < gEtaRange->GetMin())
               {
                  v0->SetRnrSelf(kFALSE);
                  v0->SetRnrChildren(kFALSE);
                  v0->GetPosTrack()->SetRnrSelf(kFALSE);
                  v0->GetNegTrack()->SetRnrSelf(kFALSE);
                  v0->GetPointingLine()->SetRnrSelf(kFALSE);
                  v0->GetPosTrack()->SetRnrChildren(kFALSE);
                  v0->GetNegTrack()->SetRnrChildren(kFALSE);
                  v0->GetPointingLine()->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            if(gEtaRange->GetMin() == gEtaRange->GetLimitMin())
            {

               if(v0->GetEta() > gEtaRange->GetMax())
               {
                  v0->SetRnrSelf(kFALSE);
                  v0->SetRnrChildren(kFALSE);
                  v0->GetPosTrack()->SetRnrSelf(kFALSE);
                  v0->GetNegTrack()->SetRnrSelf(kFALSE);
                  v0->GetPointingLine()->SetRnrSelf(kFALSE);
                  v0->GetPosTrack()->SetRnrChildren(kFALSE);
                  v0->GetNegTrack()->SetRnrChildren(kFALSE);
                  v0->GetPointingLine()->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            if(v0->GetEta() < gEtaRange->GetMin() || v0->GetEta() > gEtaRange->GetMax())
            {
               v0->SetRnrSelf(kFALSE);
               v0->SetRnrChildren(kFALSE);
               v0->GetPosTrack()->SetRnrSelf(kFALSE);
               v0->GetNegTrack()->SetRnrSelf(kFALSE);
               v0->GetPointingLine()->SetRnrSelf(kFALSE);
               v0->GetPosTrack()->SetRnrChildren(kFALSE);
               v0->GetNegTrack()->SetRnrChildren(kFALSE);
               v0->GetPointingLine()->SetRnrChildren(kFALSE);

               continue;
            }

         }

         v0->SetRnrSelf(kTRUE);
         v0->SetRnrChildren(kTRUE);
               v0->GetPosTrack()->SetRnrSelf(kTRUE);
               v0->GetNegTrack()->SetRnrSelf(kTRUE);
               v0->GetPointingLine()->SetRnrSelf(kTRUE);
               v0->GetPosTrack()->SetRnrChildren(kTRUE);
               v0->GetNegTrack()->SetRnrChildren(kTRUE);
               v0->GetPointingLine()->SetRnrChildren(kTRUE);

      }

cout << "doszlo v0s" << endl;

      i = gEve->GetEventScene()->FirstChild()->FindChild("ESD cascade")->BeginChildren();
      j = gEve->GetEventScene()->FirstChild()->FindChild("ESD cascade")->EndChildren();

      for(k = i; k != j; k++)
      {

         AliEveCascade *cascade = dynamic_cast<AliEveCascade*>((TEveElement*) *k);
            
         if(gCutOnP->IsOn())
         {

            if(gPRange->GetMax() == gPRange->GetLimitMax())
            {

               if(cascade->GetPtot() < gPRange->GetMin())
               {
                  cascade->SetRnrSelf(kFALSE);
                  cascade->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            else
            {
               if(cascade->GetPtot() < gPRange->GetMin() || cascade->GetPtot() > gPRange->GetMax())
               {
                  cascade->SetRnrSelf(kFALSE);
                  cascade->SetRnrChildren(kFALSE);
                  continue;
               }

            }

         }

         if(gCutOnPt->IsOn())
         {

            if(gPtRange->GetMax() == gPtRange->GetLimitMax())
            {

               if(cascade->GetPt() < gPtRange->GetMin())
               {
                  cascade->SetRnrSelf(kFALSE);
                  cascade->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            else
            {
               if(cascade->GetPt() < gPtRange->GetMin() || cascade->GetPt() > gPtRange->GetMax())
               {
                  cascade->SetRnrSelf(kFALSE);
                  cascade->SetRnrChildren(kFALSE);
                  continue;
               }

            }

         }

         if(gCutOnEta->IsOn())
         {

            if(gEtaRange->GetMax() == gEtaRange->GetLimitMax() && gEtaRange->GetMax() == gEtaRange->GetLimitMax())
               continue;

            if(gEtaRange->GetMax() == gEtaRange->GetLimitMax())
            {
               if(cascade->GetEta() < gEtaRange->GetMin())
               {
                  cascade->SetRnrSelf(kFALSE);
                  cascade->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            if(gEtaRange->GetMin() == gEtaRange->GetLimitMin())
            {

               if(cascade->GetEta() > gEtaRange->GetMax())
               {
                  cascade->SetRnrSelf(kFALSE);
                  cascade->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            if(cascade->GetEta() < gEtaRange->GetMin() || cascade->GetEta() > gEtaRange->GetMax())
            {
               cascade->SetRnrSelf(kFALSE);
               cascade->SetRnrChildren(kFALSE);
               continue;
            }

         }

         cascade->SetRnrSelf(kTRUE);
         cascade->SetRnrChildren(kTRUE);

      }


      i = gEve->GetEventScene()->FirstChild()->FindChild("ESD kink")->BeginChildren();
      j = gEve->GetEventScene()->FirstChild()->FindChild("ESD kink")->EndChildren();

      for(k = i; k != j; k++)
      {

         AliEveKink *kink = dynamic_cast<AliEveKink*>((TEveElement*) *k);
            
         if(gCutOnP->IsOn())
         {

            if(gPRange->GetMax() == gPRange->GetLimitMax())
            {

               if(kink->GetESDTrack()->P() < gPRange->GetMin())
               {
                  kink->SetRnrSelf(kFALSE);
                  kink->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            else
            {
               if(kink->GetESDTrack()->P() < gPRange->GetMin() || kink->GetESDTrack()->P() > gPRange->GetMax())
               {
                  kink->SetRnrSelf(kFALSE);
                  kink->SetRnrChildren(kFALSE);
                  continue;
               }

            }

         }

         if(gCutOnPt->IsOn())
         {

            if(gPtRange->GetMax() == gPtRange->GetLimitMax())
            {

               if(kink->GetESDTrack()->Pt() < gPtRange->GetMin())
               {
                  kink->SetRnrSelf(kFALSE);
                  kink->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            else
            {
               if(kink->GetESDTrack()->Pt() < gPtRange->GetMin() || kink->GetESDTrack()->Pt() > gPtRange->GetMax())
               {
                  kink->SetRnrSelf(kFALSE);
                  kink->SetRnrChildren(kFALSE);
                  continue;
               }

            }

         }

         if(gCutOnEta->IsOn())
         {

            if(gEtaRange->GetMax() == gEtaRange->GetLimitMax() && gEtaRange->GetMax() == gEtaRange->GetLimitMax())
               continue;

            if(gEtaRange->GetMax() == gEtaRange->GetLimitMax())
            {
               if(kink->GetESDTrack()->Eta() < gEtaRange->GetMin())
               {
                  kink->SetRnrSelf(kFALSE);
                  kink->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            if(gEtaRange->GetMin() == gEtaRange->GetLimitMin())
            {

               if(kink->GetESDTrack()->Eta() > gEtaRange->GetMax())
               {
                  kink->SetRnrSelf(kFALSE);
                  kink->SetRnrChildren(kFALSE);
                  continue;
               }

            }

            if(kink->GetESDTrack()->Eta() < gEtaRange->GetMin() || kink->GetESDTrack()->Eta() > gEtaRange->GetMax())
            {
               kink->SetRnrSelf(kFALSE);
               kink->SetRnrChildren(kFALSE);
               continue;
            }

         }

         kink->SetRnrSelf(kTRUE);
         kink->SetRnrChildren(kTRUE);

      }
*/

  }

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *mv = AliEveMultiView::Instance();

  mv->DestroyEventRPhi();
  mv->DestroyEventRhoZ();

  mv->ImportEventRPhi(top);
  mv->ImportEventRhoZ(top);

  gEve->FullRedraw3D(kFALSE, kTRUE);

}

//______________________________________________________________________________

void SetCutsWindow::CloseTab()
{

      TEveBrowser *browser = gEve->GetBrowser();
      Int_t current = browser->GetTabLeft()->GetCurrent();

      browser->GetTabLeft()->RemoveTab(current);

}

//______________________________________________________________________________

void alieve_set_cuts()
{

      TEveBrowser *browser = gEve->GetBrowser();

      browser->StartEmbedding(TRootBrowser::kLeft);

      new SetCutsWindow();

      browser->StopEmbedding("Cut Selection");

}
