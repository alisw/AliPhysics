//
//  AliEveCutsWindow.h
//
//  Created by Jeremi Niedziela on 1/12/15
//  main author: Pawel Debski 2010
//

#ifndef __AliEveCutsWindow__
#define __AliEveCutsWindow__

#include <TGButton.h>
#include <TGSlider.h>
#include <TGNumberEntry.h>
#include <TGLOverlayButton.h>
#include <TG3DLine.h>
#include <TGLabel.h>
#include <TGComboBox.h>

#include <TEveGValuators.h>

class AliEveCutsWindow : public TGMainFrame
{
public:
    AliEveCutsWindow();
    ~AliEveCutsWindow(){}
    
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
    void SaveRhoZView();
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
    void SetValues();
    void SaveMacro();
    void LoadMacro();
    void CloseTab();
    void Macro1();
    void Macro2();
    void Macro3();
    void Macro4();
    void Macro5();
    
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
    TGComboBox* gTrackColorScale;
    TGComboBox* gBkgColorList;
    TGTextButton* gBkgColorButton;
    TGLOverlayButton *gOverlayButton3D;
    TGLOverlayButton *gOverlayButtonRPhi;
    TGLOverlayButton *gOverlayButtonRhoZ;
    Bool_t gDrawHistograms[12];
    TGHorizontal3DLine *separator;
    
    const char *gPictureSaveAsTypes[4];
    const char *gMacroSaveAsTypes[4];
    
private:
    AliEveCutsWindow(const AliEveCutsWindow&);
    AliEveCutsWindow& operator=(const AliEveCutsWindow&);
    
    ClassDef(AliEveCutsWindow, 0);
};

#endif
