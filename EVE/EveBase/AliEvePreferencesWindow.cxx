//
//  AliEvePreferencesWindow.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 26/06/15.
//
//

#include "AliEvePreferencesWindow.h"
#include "AliEveInit.h"
#include "AliEveEventManager.h"
#include "AliEveDataSourceOffline.h"
#include "AliEveDataSource.h"

#include <TGLabel.h>
#include <TSystem.h>
#include <TColor.h>

#include <iostream>

using namespace std;

AliEvePreferencesWindow *AliEvePreferencesWindow::fInstance=0;

ClassImp(AliEvePreferencesWindow);


AliEvePreferencesWindow::AliEvePreferencesWindow():
  TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame),
  fTrackWidth(0),
  fDashNoRefit(0),
  fDrawNoRefit(0),
  fTracksByPID(0),
  fTracksByCategory(0),

  fShowV0s(0),
  fShowCascades(0),
  fShowRawData(0),
  fShowPrimaryVertex(0),
  fShowHits(0),
  fShowDigits(0),
  fShowClusters(0),
  fShowKinks(0),
      
  fLogbookHost(0),
  fLogbookPort(0),
  fLogbookDatabase(0),
  fLogbookUser(0),
  fLogbookPassword(0),
      
  fShowMuon(0),
  fShowHLTESDTree(0),
  fOCDBpath(0),
  fAutoload(0),
  fAliceLive(0),

  fSaveAndExitButton(0),
  fCancel(0)
{
    SetName("AliEve Preferences");
    SetLayoutBroken(kTRUE);
    InitWindow();
    ReadFromConfigFile();
}

AliEvePreferencesWindow::~AliEvePreferencesWindow()
{
    delete fTrackWidth;
    delete fDashNoRefit;
    delete fDrawNoRefit;
    delete fTracksByPID;
    delete fTracksByCategory;
    
    delete fShowV0s;
    delete fShowCascades;
    delete fShowRawData;
    delete fShowPrimaryVertex;
    delete fShowHits;
    delete fShowDigits;
    delete fShowClusters;
    delete fShowKinks;
    
    delete fLogbookHost;
    delete fLogbookPort;
    delete fLogbookDatabase;
    delete fLogbookUser;
    delete fLogbookPassword;
    
    delete fShowMuon;
    delete fShowHLTESDTree;
    delete fOCDBpath;
    delete fAutoload;
    delete fAliceLive;
    
    delete fSaveAndExitButton;
    delete fCancel;
}

AliEvePreferencesWindow* AliEvePreferencesWindow::Instance()
{
    if(!fInstance){fInstance = new AliEvePreferencesWindow();}
    return fInstance;
}

void AliEvePreferencesWindow::onExit(bool save)
{
    cout<<"Quitting alieve preferences window";
    if(save){
        cout<<" with saving"<<endl;
        SaveToConfigFile();
        ApplyChanges();
    }
    else{
        cout<<" without saveing"<<endl;
    }
    if(fInstance){delete fInstance;fInstance=0;}
    cout<<" -- OK"<<endl;
}

void AliEvePreferencesWindow::ReadFromConfigFile()
{
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    fTrackWidth->SetNumber(settings.GetValue("tracks.width",2));
    fDashNoRefit->SetOn(settings.GetValue("tracks.noRefit.dash",true));
    fDrawNoRefit->SetOn(settings.GetValue("tracks.noRefit.show",true));
    fTracksByPID->SetOn(settings.GetValue("tracks.byType.show",true));
    fTracksByCategory->SetOn(settings.GetValue("tracks.byCategory.show",false));
    fShowV0s->SetOn(settings.GetValue("V0s.show",false));
    fShowCascades->SetOn(settings.GetValue("cascades.show",false));
    fShowRawData->SetOn(settings.GetValue("rawData.show",false));
    fShowPrimaryVertex->SetOn(settings.GetValue("primary.vertex.show",false));
    fShowHits->SetOn(settings.GetValue("hits.show",false));
    fShowDigits->SetOn(settings.GetValue("digits.show",false));
    fShowClusters->SetOn(settings.GetValue("clusters.show",false));
    fShowKinks->SetOn(settings.GetValue("kinks.show",false));
    fLogbookHost->SetText(settings.GetValue("logbook.host", "hostname"));
    fLogbookPort->SetNumber(settings.GetValue("logbook.port", 0));
    fLogbookDatabase->SetText(settings.GetValue("logbook.db", "database"));
    fLogbookUser->SetText(settings.GetValue("logbook.user", "username"));
    fLogbookPassword->SetText(settings.GetValue("logbook.pass", "password"));
    fShowMuon->SetOn(settings.GetValue("MUON.show", true));
    fShowHLTESDTree->SetOn(settings.GetValue("HLTESDtree.show", false));
    fOCDBpath->SetText(settings.GetValue("OCDB.default.path","local://$ALICE_ROOT/../src/OCDB"));
    fAutoload->SetOn(settings.GetValue("events.autoload.set",false));
    fAliceLive->SetOn(settings.GetValue("ALICE_LIVE.send",false));
    
    fColorSelect2658->SetColor(TColor::Number2Pixel(settings.GetValue("TPC.color",3)));//tpc
    fColorSelect2668->SetColor(TColor::Number2Pixel(settings.GetValue("TOF.color",930)));//tof
    fColorSelect2669->SetColor(TColor::Number2Pixel(settings.GetValue("TRD.color",920)));//trd
    fColorSelect2670->SetColor(TColor::Number2Pixel(settings.GetValue("MUON.color",920)));//muon
    fColorSelect2671->SetColor(TColor::Number2Pixel(settings.GetValue("ITS.SPD.color",924)));//spd
    fColorSelect2672->SetColor(TColor::Number2Pixel(settings.GetValue("ITS.SDD.color",925)));//sdd
    fColorSelect2673->SetColor(TColor::Number2Pixel(settings.GetValue("ITS.SSD.color",926)));//ssd
    fColorSelect2674->SetColor(TColor::Number2Pixel(settings.GetValue("PHOS.color",953)));//phos
    fColorSelect2675->SetColor(TColor::Number2Pixel(settings.GetValue("EMCAL.color",953)));//emcal
    fColorSelect2676->SetColor(TColor::Number2Pixel(settings.GetValue("HMPID.color",5)));//hmpid
    
    fColorSelect2693->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.electron",600)));//electron
    fColorSelect2694->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.muon",416)));//muon
    fColorSelect2695->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.pion",632)));//pion
    fColorSelect2696->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.kaon",400)));//kaon
    fColorSelect2697->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.proton",797)));//proton
    fColorSelect2698->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.deuteron",797)));//deuteron
    fColorSelect2699->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.triton",797)));//triton
    fColorSelect2700->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.he3",797)));//he3
    fColorSelect2701->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.alpha",403)));//alpha
    fColorSelect2702->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.photon",0)));//photon
    fColorSelect2703->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.pi0",616)));//pi0
    fColorSelect2704->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.neutron",900)));//neutron
    fColorSelect2705->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.kaon0",801)));//kaon0
    fColorSelect2706->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.elecon",920)));//elecon
    fColorSelect2707->SetColor(TColor::Number2Pixel(settings.GetValue("tracks.byType.unknown",920)));//unknown

    
    Layout();
}

void AliEvePreferencesWindow::SaveToConfigFile()
{
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    
    settings.SetValue("tracks.width",fTrackWidth->GetNumber());
    
    settings.SetValue("tracks.noRefit.dash",fDashNoRefit->IsOn());
    settings.SetValue("tracks.noRefit.show",fDrawNoRefit->IsOn());
    settings.SetValue("tracks.byType.show",fTracksByPID->IsOn());
    settings.SetValue("tracks.byCategory.show",fTracksByCategory->IsOn());
    
    settings.SetValue("V0s.show",fShowV0s->IsOn());
    settings.SetValue("cascades.show",fShowCascades->IsOn());
    settings.SetValue("rawData.show",fShowRawData->IsOn());
    settings.SetValue("primary.vertex.show",fShowPrimaryVertex->IsOn());
    settings.SetValue("hits.show",fShowHits->IsOn());
    settings.SetValue("digits.show",fShowDigits->IsOn());
    settings.SetValue("clusters.show",fShowClusters->IsOn());
    settings.SetValue("kinks.show",fShowKinks->IsOn());
    
    settings.SetValue("logbook.host",fLogbookHost->GetText());
    settings.SetValue("logbook.port", fLogbookPort->GetNumber());
    settings.SetValue("logbook.db", fLogbookDatabase->GetText());
    settings.SetValue("logbook.user", fLogbookUser->GetText());
    settings.SetValue("logbook.pass", fLogbookPassword->GetText());
    
    settings.SetValue("MUON.show",fShowMuon->IsOn());
    settings.SetValue("HLTESDtree.show", fShowHLTESDTree->IsOn());
    settings.SetValue("OCDB.default.path",fOCDBpath->GetText());
    settings.SetValue("events.autoload.set",fAutoload->IsOn());
    settings.SetValue("ALICE_LIVE.send",fAliceLive->IsOn());
    
    
    settings.SetValue("TPC.color",TColor::GetColor(fColorSelect2658->GetColor()));//tpc
    settings.SetValue("TOF.color",TColor::GetColor(fColorSelect2668->GetColor()));//tof
    settings.SetValue("TRD.color",TColor::GetColor(fColorSelect2669->GetColor()));//trd
    settings.SetValue("MUON.color",TColor::GetColor(fColorSelect2670->GetColor()));//muon
    settings.SetValue("ITS.SPD.color",TColor::GetColor(fColorSelect2671->GetColor()));//spd
    settings.SetValue("ITS.SDD.color",TColor::GetColor(fColorSelect2672->GetColor()));//sdd
    settings.SetValue("ITS.SSD.color",TColor::GetColor(fColorSelect2673->GetColor()));//ssd
    settings.SetValue("PHOS.color",TColor::GetColor(fColorSelect2674->GetColor()));//phos
    settings.SetValue("EMCAL.color",TColor::GetColor(fColorSelect2675->GetColor()));//emcal
    settings.SetValue("HMPID.color",TColor::GetColor(fColorSelect2676->GetColor()));//hmpid
    
    settings.SetValue("tracks.byType.electron",TColor::GetColor(fColorSelect2693->GetColor()));//electron
    settings.SetValue("tracks.byType.muon",TColor::GetColor(fColorSelect2694->GetColor()));//muon
    settings.SetValue("tracks.byType.pion",TColor::GetColor(fColorSelect2695->GetColor()));//pion
    settings.SetValue("tracks.byType.kaon",TColor::GetColor(fColorSelect2696->GetColor()));//kaon
    settings.SetValue("tracks.byType.proton",TColor::GetColor(fColorSelect2697->GetColor()));//proton
    settings.SetValue("tracks.byType.deuteron",TColor::GetColor(fColorSelect2698->GetColor()));//deuteron
    settings.SetValue("tracks.byType.triton",TColor::GetColor(fColorSelect2699->GetColor()));//triton
    settings.SetValue("tracks.byType.he3",TColor::GetColor(fColorSelect2700->GetColor()));//he3
    settings.SetValue("tracks.byType.alpha",TColor::GetColor(fColorSelect2701->GetColor()));//alpha
    settings.SetValue("tracks.byType.photon",TColor::GetColor(fColorSelect2702->GetColor()));//photon
    settings.SetValue("tracks.byType.pi0",TColor::GetColor(fColorSelect2703->GetColor()));//pi0
    settings.SetValue("tracks.byType.neutron",TColor::GetColor(fColorSelect2704->GetColor()));//neutron
    settings.SetValue("tracks.byType.kaon0",TColor::GetColor(fColorSelect2705->GetColor()));//kaon0
    settings.SetValue("tracks.byType.elecon",TColor::GetColor(fColorSelect2706->GetColor()));//elecon
    settings.SetValue("tracks.byType.unknown",TColor::GetColor(fColorSelect2707->GetColor()));//unknown
    
    settings.WriteFile(Form("%s/.eve_config",gSystem->Getenv("HOME")), kEnvAll);
}

void AliEvePreferencesWindow::ApplyChanges()
{
    TEnv settings;
    AliEveInit::GetConfig(&settings);
 
    AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)AliEveEventManager::Instance()->GetDataSourceOffline();
    
//    const Text_t* esdfile = 0;
//    dataSource->SetESDFileName(esdfile);
    
 
    AliEveInit::AddMacros();
 
    AliEveEventManager *man =  AliEveEventManager::Instance();
 
    man->SetAutoLoad(settings.GetValue("events.autoload.set",false));
    
    man->InitOCDB(man->GetCurrentRun());
}

Bool_t AliEvePreferencesWindow::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
    switch (GET_MSG(msg))
    {
        case kC_COMMAND:
            switch (GET_SUBMSG(msg))
        {
            case kCM_BUTTON:
                switch(parm1)
            {
                case 0:onExit(false);break;
                case 1:onExit(true);break;
                default:break;
            }
                break;
            default:break;
        }
            break;
        default:break;
    }
    
    return false;
}

void AliEvePreferencesWindow::InitWindow()
{
    // main frame
//    TGMainFrame *this = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
//    this->SetName("this");
//    this->SetLayoutBroken(kTRUE);
    
    // composite frame
    TGCompositeFrame *fMainFrame2391 = new TGCompositeFrame(this,570,560,kVerticalFrame);
    fMainFrame2391->SetName("fMainFrame2391");
    fMainFrame2391->SetLayoutBroken(kTRUE);
    
    // composite frame
    TGCompositeFrame *fMainFrame2342 = new TGCompositeFrame(fMainFrame2391,580,530,kVerticalFrame);
    fMainFrame2342->SetName("fMainFrame2342");
    fMainFrame2342->SetLayoutBroken(kTRUE);
    
    // tab widget
    TGTab *fTab2632 = new TGTab(fMainFrame2342,570,526);
    
    // container of "Tab1"
    TGCompositeFrame *fCompositeFrame2635;
    fCompositeFrame2635 = fTab2632->AddTab("Appearance");
    fCompositeFrame2635->SetLayoutManager(new TGVerticalLayout(fCompositeFrame2635));
    fCompositeFrame2635->SetLayoutBroken(kTRUE);
    
    // composite frame
    TGCompositeFrame *fCompositeFrame2636 = new TGCompositeFrame(fCompositeFrame2635,570,500,kVerticalFrame);
    fCompositeFrame2636->SetLayoutBroken(kTRUE);
    
    // "Elements to be shown" group frame
    TGGroupFrame *fGroupFrame2637 = new TGGroupFrame(fCompositeFrame2636,"Elements to be shown");
    fGroupFrame2637->SetLayoutBroken(kTRUE);
    fShowV0s = new TGCheckButton(fGroupFrame2637,"V0s");
    fShowV0s->SetTextJustify(36);
    fShowV0s->SetMargins(0,0,0,0);
    fShowV0s->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowV0s, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowV0s->MoveResize(18,18,47,17);
    fShowCascades = new TGCheckButton(fGroupFrame2637,"cascades");
    fShowCascades->SetTextJustify(36);
    fShowCascades->SetMargins(0,0,0,0);
    fShowCascades->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowCascades, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowCascades->MoveResize(18,39,79,17);
    fShowRawData = new TGCheckButton(fGroupFrame2637,"RAW data");
    fShowRawData->SetTextJustify(36);
    fShowRawData->SetMargins(0,0,0,0);
    fShowRawData->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowRawData, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowRawData->MoveResize(18,60,81,17);
    fShowPrimaryVertex = new TGCheckButton(fGroupFrame2637,"primary vertex");
    fShowPrimaryVertex->SetTextJustify(36);
    fShowPrimaryVertex->SetMargins(0,0,0,0);
    fShowPrimaryVertex->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowPrimaryVertex, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowPrimaryVertex->MoveResize(18,81,104,17);
    fShowHits = new TGCheckButton(fGroupFrame2637,"hits");
    fShowHits->SetTextJustify(36);
    fShowHits->SetMargins(0,0,0,0);
    fShowHits->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowHits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowHits->MoveResize(18,102,44,17);
    fShowDigits = new TGCheckButton(fGroupFrame2637,"digits");
    fShowDigits->SetTextJustify(36);
    fShowDigits->SetMargins(0,0,0,0);
    fShowDigits->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowDigits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowDigits->MoveResize(18,123,54,17);
    fShowClusters = new TGCheckButton(fGroupFrame2637,"clusters");
    fShowClusters->SetTextJustify(36);
    fShowClusters->SetMargins(0,0,0,0);
    fShowClusters->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowClusters, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowClusters->MoveResize(18,144,68,17);
    fShowKinks = new TGCheckButton(fGroupFrame2637,"kinks");
    fShowKinks->SetTextJustify(36);
    fShowKinks->SetMargins(0,0,0,0);
    fShowKinks->SetWrapLength(-1);
    fGroupFrame2637->AddFrame(fShowKinks, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowKinks->MoveResize(18,165,53,17);
    
    fGroupFrame2637->SetLayoutManager(new TGVerticalLayout(fGroupFrame2637));
    fGroupFrame2637->Resize(180,200);
    fCompositeFrame2636->AddFrame(fGroupFrame2637, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2637->MoveResize(5,5,180,200);
    
    // "Tracks settings" group frame
    TGGroupFrame *fGroupFrame2646 = new TGGroupFrame(fCompositeFrame2636,"Tracks settings");
    fGroupFrame2646->SetLayoutBroken(kTRUE);
    
    TGFont *ufont;         // will reflect user font changes
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    TGGC   *uGC;           // will reflect user GC changes
    // graphics context changes
    GCValues_t vall2647;
    vall2647.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall2647.fForeground);
    gClient->GetColorByName("#e0e0e0",vall2647.fBackground);
    vall2647.fFillStyle = kFillSolid;
    vall2647.fFont = ufont->GetFontHandle();
    vall2647.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall2647, kTRUE);
    TGLabel *fLabel2647 = new TGLabel(fGroupFrame2646,"width",uGC->GetGC());
    fLabel2647->SetTextJustify(36);
    fLabel2647->SetMargins(0,0,0,0);
    fLabel2647->SetWrapLength(-1);
    fGroupFrame2646->AddFrame(fLabel2647, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2647->MoveResize(18,18,29,16);
    fTrackWidth = new TGNumberEntry(fGroupFrame2646, (Double_t) 0,9,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 0,1,5);
    fGroupFrame2646->AddFrame(fTrackWidth, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTrackWidth->MoveResize(18,38,80,22);
    fDashNoRefit = new TGCheckButton(fGroupFrame2646,"dash no-refit tracks");
    fDashNoRefit->SetTextJustify(36);
    fDashNoRefit->SetMargins(0,0,0,0);
    fDashNoRefit->SetWrapLength(-1);
    fGroupFrame2646->AddFrame(fDashNoRefit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fDashNoRefit->MoveResize(18,64,132,17);
    fDrawNoRefit = new TGCheckButton(fGroupFrame2646,"draw no-refit tracks");
    fDrawNoRefit->SetTextJustify(36);
    fDrawNoRefit->SetMargins(0,0,0,0);
    fDrawNoRefit->SetWrapLength(-1);
    fGroupFrame2646->AddFrame(fDrawNoRefit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fDrawNoRefit->MoveResize(18,85,132,17);
    fTracksByPID = new TGCheckButton(fGroupFrame2646,"tracks by PID");
    fTracksByPID->SetTextJustify(36);
    fTracksByPID->SetMargins(0,0,0,0);
    fTracksByPID->SetWrapLength(-1);
    fGroupFrame2646->AddFrame(fTracksByPID, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTracksByPID->MoveResize(18,106,100,17);
    fTracksByCategory = new TGCheckButton(fGroupFrame2646,"tracks by category");
    fTracksByCategory->SetTextJustify(36);
    fTracksByCategory->SetMargins(0,0,0,0);
    fTracksByCategory->SetWrapLength(-1);
    fGroupFrame2646->AddFrame(fTracksByCategory, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTracksByCategory->MoveResize(18,127,129,17);
    
    fGroupFrame2646->SetLayoutManager(new TGVerticalLayout(fGroupFrame2646));
    fGroupFrame2646->Resize(180,168);
    fCompositeFrame2636->AddFrame(fGroupFrame2646, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2646->MoveResize(5,200,180,168);
    
    // "Detectors colors" group frame
    TGGroupFrame *fGroupFrame2656 = new TGGroupFrame(fCompositeFrame2636,"Detectors colors");
    fGroupFrame2656->SetLayoutBroken(kTRUE);
    TGLabel *fLabel2657 = new TGLabel(fGroupFrame2656,"TPC");
    fLabel2657->SetTextJustify(kTextLeft);
    fLabel2657->SetMargins(0,0,0,0);
    fLabel2657->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2657, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2657->MoveResize(5,20,80,15);
    
    // color select widget
    ULong_t ColPar36;
    gClient->GetColorByName("#ff00ff", ColPar36);
    fColorSelect2658 = new TGColorSelect(fGroupFrame2656, ColPar36, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2658, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2658->MoveResize(90,20,45,23);
    TGLabel *fLabel2659 = new TGLabel(fGroupFrame2656,"TOF");
    fLabel2659->SetTextJustify(kTextLeft);
    fLabel2659->SetMargins(0,0,0,0);
    fLabel2659->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2659, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2659->MoveResize(5,50,80,15);
    TGLabel *fLabel2660 = new TGLabel(fGroupFrame2656,"TRD");
    fLabel2660->SetTextJustify(kTextLeft);
    fLabel2660->SetMargins(0,0,0,0);
    fLabel2660->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2660, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2660->MoveResize(5,80,80,15);
    TGLabel *fLabel2661 = new TGLabel(fGroupFrame2656,"MUON");
    fLabel2661->SetTextJustify(kTextLeft);
    fLabel2661->SetMargins(0,0,0,0);
    fLabel2661->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2661, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2661->MoveResize(5,110,80,15);
    TGLabel *fLabel2662 = new TGLabel(fGroupFrame2656,"SPD");
    fLabel2662->SetTextJustify(kTextLeft);
    fLabel2662->SetMargins(0,0,0,0);
    fLabel2662->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2662, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2662->MoveResize(5,140,80,15);
    TGLabel *fLabel2663 = new TGLabel(fGroupFrame2656,"SDD");
    fLabel2663->SetTextJustify(kTextLeft);
    fLabel2663->SetMargins(0,0,0,0);
    fLabel2663->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2663, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2663->MoveResize(5,170,80,15);
    TGLabel *fLabel2664 = new TGLabel(fGroupFrame2656,"SSD");
    fLabel2664->SetTextJustify(kTextLeft);
    fLabel2664->SetMargins(0,0,0,0);
    fLabel2664->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2664, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2664->MoveResize(5,200,80,15);
    TGLabel *fLabel2665 = new TGLabel(fGroupFrame2656,"PHOS");
    fLabel2665->SetTextJustify(kTextLeft);
    fLabel2665->SetMargins(0,0,0,0);
    fLabel2665->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2665, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2665->MoveResize(5,230,80,15);
    TGLabel *fLabel2666 = new TGLabel(fGroupFrame2656,"EMCAL");
    fLabel2666->SetTextJustify(kTextLeft);
    fLabel2666->SetMargins(0,0,0,0);
    fLabel2666->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2666, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2666->MoveResize(5,260,80,15);
    TGLabel *fLabel2667 = new TGLabel(fGroupFrame2656,"HMPID");
    fLabel2667->SetTextJustify(kTextLeft);
    fLabel2667->SetMargins(0,0,0,0);
    fLabel2667->SetWrapLength(-1);
    fGroupFrame2656->AddFrame(fLabel2667, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2667->MoveResize(5,290,80,15);
    
    // color select widget
    ULong_t ColPar37;
    gClient->GetColorByName("#ff00ff", ColPar37);
    fColorSelect2668 = new TGColorSelect(fGroupFrame2656, ColPar37, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2668, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2668->MoveResize(90,50,45,23);
    
    // color select widget
    ULong_t ColPar38;
    gClient->GetColorByName("#ff00ff", ColPar38);
    fColorSelect2669 = new TGColorSelect(fGroupFrame2656, ColPar38, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2669, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2669->MoveResize(90,80,45,23);
    
    // color select widget
    ULong_t ColPar39;
    gClient->GetColorByName("#ff00ff", ColPar39);
    fColorSelect2670 = new TGColorSelect(fGroupFrame2656, ColPar39, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2670, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2670->MoveResize(90,110,45,23);
    
    // color select widget
    ULong_t ColPar40;
    gClient->GetColorByName("#ff00ff", ColPar40);
    fColorSelect2671 = new TGColorSelect(fGroupFrame2656, ColPar40, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2671, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2671->MoveResize(90,140,45,23);
    
    // color select widget
    ULong_t ColPar41;
    gClient->GetColorByName("#ff00ff", ColPar41);
    fColorSelect2672 = new TGColorSelect(fGroupFrame2656, ColPar41, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2672, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2672->MoveResize(90,170,45,23);
    
    // color select widget
    ULong_t ColPar42;
    gClient->GetColorByName("#ff00ff", ColPar42);
    fColorSelect2673 = new TGColorSelect(fGroupFrame2656, ColPar42, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2673, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2673->MoveResize(90,200,45,23);
    
    // color select widget
    ULong_t ColPar43;
    gClient->GetColorByName("#ff00ff", ColPar43);
    fColorSelect2674 = new TGColorSelect(fGroupFrame2656, ColPar43, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2674, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2674->MoveResize(90,230,45,23);
    
    // color select widget
    ULong_t ColPar44;
    gClient->GetColorByName("#ff00ff", ColPar44);
    fColorSelect2675 = new TGColorSelect(fGroupFrame2656, ColPar44, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2675, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2675->MoveResize(90,260,45,23);
    
    // color select widget
    ULong_t ColPar45;
    gClient->GetColorByName("#ff00ff", ColPar45);
    fColorSelect2676 = new TGColorSelect(fGroupFrame2656, ColPar45, -1);
    
    fGroupFrame2656->AddFrame(fColorSelect2676, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2676->MoveResize(90,290,45,23);
    
    fGroupFrame2656->SetLayoutManager(new TGVerticalLayout(fGroupFrame2656));
    fGroupFrame2656->Resize(180,360);
    fCompositeFrame2636->AddFrame(fGroupFrame2656, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2656->MoveResize(190,5,180,360);
    
    // "Tracks colors" group frame
    TGGroupFrame *fGroupFrame2677 = new TGGroupFrame(fCompositeFrame2636,"Tracks colors");
    fGroupFrame2677->SetLayoutBroken(kTRUE);
    TGLabel *fLabel2678 = new TGLabel(fGroupFrame2677,"electron");
    fLabel2678->SetTextJustify(kTextLeft);
    fLabel2678->SetMargins(0,0,0,0);
    fLabel2678->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2678, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2678->MoveResize(5,20,80,15);
    TGLabel *fLabel2679 = new TGLabel(fGroupFrame2677,"muon");
    fLabel2679->SetTextJustify(kTextLeft);
    fLabel2679->SetMargins(0,0,0,0);
    fLabel2679->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2679, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2679->MoveResize(5,50,80,15);
    TGLabel *fLabel2680 = new TGLabel(fGroupFrame2677,"pion");
    fLabel2680->SetTextJustify(kTextLeft);
    fLabel2680->SetMargins(0,0,0,0);
    fLabel2680->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2680, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2680->MoveResize(5,80,80,15);
    TGLabel *fLabel2681 = new TGLabel(fGroupFrame2677,"kaon");
    fLabel2681->SetTextJustify(kTextLeft);
    fLabel2681->SetMargins(0,0,0,0);
    fLabel2681->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2681, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2681->MoveResize(5,110,80,15);
    TGLabel *fLabel2682 = new TGLabel(fGroupFrame2677,"proton");
    fLabel2682->SetTextJustify(kTextLeft);
    fLabel2682->SetMargins(0,0,0,0);
    fLabel2682->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2682, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2682->MoveResize(5,140,80,15);
    TGLabel *fLabel2683 = new TGLabel(fGroupFrame2677,"deuteron");
    fLabel2683->SetTextJustify(kTextLeft);
    fLabel2683->SetMargins(0,0,0,0);
    fLabel2683->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2683, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2683->MoveResize(5,170,80,15);
    TGLabel *fLabel2684 = new TGLabel(fGroupFrame2677,"triton");
    fLabel2684->SetTextJustify(kTextLeft);
    fLabel2684->SetMargins(0,0,0,0);
    fLabel2684->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2684, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2684->MoveResize(5,200,80,15);
    TGLabel *fLabel2685 = new TGLabel(fGroupFrame2677,"he3");
    fLabel2685->SetTextJustify(kTextLeft);
    fLabel2685->SetMargins(0,0,0,0);
    fLabel2685->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2685, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2685->MoveResize(5,230,80,15);
    TGLabel *fLabel2686 = new TGLabel(fGroupFrame2677,"alpha");
    fLabel2686->SetTextJustify(kTextLeft);
    fLabel2686->SetMargins(0,0,0,0);
    fLabel2686->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2686, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2686->MoveResize(5,260,80,15);
    TGLabel *fLabel2687 = new TGLabel(fGroupFrame2677,"photon");
    fLabel2687->SetTextJustify(kTextLeft);
    fLabel2687->SetMargins(0,0,0,0);
    fLabel2687->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2687, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2687->MoveResize(5,290,80,15);
    TGLabel *fLabel2688 = new TGLabel(fGroupFrame2677,"pi0");
    fLabel2688->SetTextJustify(kTextLeft);
    fLabel2688->SetMargins(0,0,0,0);
    fLabel2688->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2688, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2688->MoveResize(5,320,64,18);
    TGLabel *fLabel2689 = new TGLabel(fGroupFrame2677,"neutron");
    fLabel2689->SetTextJustify(kTextLeft);
    fLabel2689->SetMargins(0,0,0,0);
    fLabel2689->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2689, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2689->MoveResize(5,350,64,18);
    TGLabel *fLabel2690 = new TGLabel(fGroupFrame2677,"kaon0");
    fLabel2690->SetTextJustify(kTextLeft);
    fLabel2690->SetMargins(0,0,0,0);
    fLabel2690->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2690, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2690->MoveResize(5,380,64,18);
    TGLabel *fLabel2691 = new TGLabel(fGroupFrame2677,"elecon");
    fLabel2691->SetTextJustify(kTextLeft);
    fLabel2691->SetMargins(0,0,0,0);
    fLabel2691->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2691, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2691->MoveResize(5,410,64,18);
    TGLabel *fLabel2692 = new TGLabel(fGroupFrame2677,"unknown");
    fLabel2692->SetTextJustify(kTextLeft);
    fLabel2692->SetMargins(0,0,0,0);
    fLabel2692->SetWrapLength(-1);
    fGroupFrame2677->AddFrame(fLabel2692, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2692->MoveResize(5,440,64,18);
    
    // color select widget
    ULong_t ColPar46;
    gClient->GetColorByName("#ff00ff", ColPar46);
    fColorSelect2693 = new TGColorSelect(fGroupFrame2677, ColPar46, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2693, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2693->MoveResize(90,20,45,23);
    
    // color select widget
    ULong_t ColPar47;
    gClient->GetColorByName("#ff00ff", ColPar47);
    fColorSelect2694 = new TGColorSelect(fGroupFrame2677, ColPar47, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2694, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2694->MoveResize(90,50,45,23);
    
    // color select widget
    ULong_t ColPar48;
    gClient->GetColorByName("#ff00ff", ColPar48);
    fColorSelect2695 = new TGColorSelect(fGroupFrame2677, ColPar48, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2695, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2695->MoveResize(90,80,45,23);
    
    // color select widget
    ULong_t ColPar49;
    gClient->GetColorByName("#ff00ff", ColPar49);
    fColorSelect2696 = new TGColorSelect(fGroupFrame2677, ColPar49, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2696, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2696->MoveResize(90,110,45,23);
    
    // color select widget
    ULong_t ColPar50;
    gClient->GetColorByName("#ff00ff", ColPar50);
    fColorSelect2697 = new TGColorSelect(fGroupFrame2677, ColPar50, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2697, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2697->MoveResize(90,140,45,23);
    
    // color select widget
    ULong_t ColPar51;
    gClient->GetColorByName("#ff00ff", ColPar51);
    fColorSelect2698 = new TGColorSelect(fGroupFrame2677, ColPar51, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2698, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2698->MoveResize(90,170,45,23);
    
    // color select widget
    ULong_t ColPar52;
    gClient->GetColorByName("#ff00ff", ColPar52);
    fColorSelect2699 = new TGColorSelect(fGroupFrame2677, ColPar52, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2699, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2699->MoveResize(90,200,45,23);
    
    // color select widget
    ULong_t ColPar53;
    gClient->GetColorByName("#ff00ff", ColPar53);
    fColorSelect2700 = new TGColorSelect(fGroupFrame2677, ColPar53, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2700, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2700->MoveResize(90,230,45,23);
    
    // color select widget
    ULong_t ColPar54;
    gClient->GetColorByName("#ff00ff", ColPar54);
    fColorSelect2701 = new TGColorSelect(fGroupFrame2677, ColPar54, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2701, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2701->MoveResize(90,260,45,23);
    
    // color select widget
    ULong_t ColPar55;
    gClient->GetColorByName("#ff00ff", ColPar55);
    fColorSelect2702 = new TGColorSelect(fGroupFrame2677, ColPar55, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2702, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2702->MoveResize(90,290,45,23);
    
    // color select widget
    ULong_t ColPar56;
    gClient->GetColorByName("#ff00ff", ColPar56);
    fColorSelect2703 = new TGColorSelect(fGroupFrame2677, ColPar56, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2703, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2703->MoveResize(90,320,45,23);
    
    // color select widget
    ULong_t ColPar57;
    gClient->GetColorByName("#ff00ff", ColPar57);
    fColorSelect2704 = new TGColorSelect(fGroupFrame2677, ColPar57, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2704, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2704->MoveResize(90,350,45,23);
    
    // color select widget
    ULong_t ColPar58;
    gClient->GetColorByName("#ff00ff", ColPar58);
    fColorSelect2705 = new TGColorSelect(fGroupFrame2677, ColPar58, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2705, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2705->MoveResize(90,380,45,23);
    
    // color select widget
    ULong_t ColPar59;
    gClient->GetColorByName("#ff00ff", ColPar59);
    fColorSelect2706 = new TGColorSelect(fGroupFrame2677, ColPar59, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2706, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2706->MoveResize(90,410,45,23);
    
    // color select widget
    ULong_t ColPar60;
    gClient->GetColorByName("#ff00ff", ColPar60);
    fColorSelect2707 = new TGColorSelect(fGroupFrame2677, ColPar60, -1);
    
    fGroupFrame2677->AddFrame(fColorSelect2707, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fColorSelect2707->MoveResize(90,440,45,23);
    
    fGroupFrame2677->SetLayoutManager(new TGVerticalLayout(fGroupFrame2677));
    fGroupFrame2677->Resize(180,480);
    fCompositeFrame2636->AddFrame(fGroupFrame2677, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2677->MoveResize(375,5,180,480);
    
    fCompositeFrame2635->AddFrame(fCompositeFrame2636, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fCompositeFrame2636->MoveResize(2,2,570,500);
    
    
    // container of "Tab2"
    TGCompositeFrame *fCompositeFrame2709;
    fCompositeFrame2709 = fTab2632->AddTab("Advanced");
    fCompositeFrame2709->SetLayoutManager(new TGVerticalLayout(fCompositeFrame2709));
    fCompositeFrame2709->SetLayoutBroken(kTRUE);
    
    // "Miscellaneous" group frame
    TGGroupFrame *fGroupFrame2710 = new TGGroupFrame(fCompositeFrame2709,"Miscellaneous");
    fGroupFrame2710->SetLayoutBroken(kTRUE);
    fShowMuon = new TGCheckButton(fGroupFrame2710,"show MUON geometry and tracks");
    fShowMuon->SetTextJustify(36);
    fShowMuon->SetMargins(0,0,0,0);
    fShowMuon->SetWrapLength(-1);
    fGroupFrame2710->AddFrame(fShowMuon, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowMuon->MoveResize(10,20,220,20);
    fShowHLTESDTree = new TGCheckButton(fGroupFrame2710,"show HLT ESD tree");
    fShowHLTESDTree->SetTextJustify(36);
    fShowHLTESDTree->SetMargins(0,0,0,0);
    fShowHLTESDTree->SetWrapLength(-1);
    fGroupFrame2710->AddFrame(fShowHLTESDTree, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowHLTESDTree->MoveResize(10,40,220,20);
    TGLabel *fLabel2713 = new TGLabel(fGroupFrame2710,"OCDB path:");
    fLabel2713->SetTextJustify(kTextLeft);
    fLabel2713->SetMargins(0,0,0,0);
    fLabel2713->SetWrapLength(-1);
    fGroupFrame2710->AddFrame(fLabel2713, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2713->MoveResize(10,100,230,20);
    fAutoload = new TGCheckButton(fGroupFrame2710,"autoload new events");
    fAutoload->SetTextJustify(36);
    fAutoload->SetMargins(0,0,0,0);
    fAutoload->SetWrapLength(-1);
    fGroupFrame2710->AddFrame(fAutoload, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fAutoload->MoveResize(10,60,230,20);
    fAliceLive = new TGCheckButton(fGroupFrame2710,"send pictures to ALICE LIVE");
    fAliceLive->SetTextJustify(36);
    fAliceLive->SetMargins(0,0,0,0);
    fAliceLive->SetWrapLength(-1);
    fGroupFrame2710->AddFrame(fAliceLive, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fAliceLive->MoveResize(10,80,230,20);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valEntry2716;
    valEntry2716.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",valEntry2716.fForeground);
    gClient->GetColorByName("#e0e0e0",valEntry2716.fBackground);
    valEntry2716.fFillStyle = kFillSolid;
    valEntry2716.fFont = ufont->GetFontHandle();
    valEntry2716.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valEntry2716, kTRUE);
    fOCDBpath = new TGTextEntry(fGroupFrame2710, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
    fOCDBpath->SetMaxLength(4096);
    fOCDBpath->SetAlignment(kTextLeft);
    fOCDBpath->SetText("fTextEntry1461");
    fOCDBpath->Resize(300,fOCDBpath->GetDefaultHeight());
    fGroupFrame2710->AddFrame(fOCDBpath, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fOCDBpath->MoveResize(10,120,300,20);
    
    fGroupFrame2710->SetLayoutManager(new TGVerticalLayout(fGroupFrame2710));
    fGroupFrame2710->Resize(320,176);
    fCompositeFrame2709->AddFrame(fGroupFrame2710, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2710->MoveResize(5,5,320,176);
    
    // "Logbook credentials" group frame
    TGGroupFrame *fGroupFrame2717 = new TGGroupFrame(fCompositeFrame2709,"Logbook credentials");
    fGroupFrame2717->SetLayoutBroken(kTRUE);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valEntry2718;
    valEntry2718.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",valEntry2718.fForeground);
    gClient->GetColorByName("#e0e0e0",valEntry2718.fBackground);
    valEntry2718.fFillStyle = kFillSolid;
    valEntry2718.fFont = ufont->GetFontHandle();
    valEntry2718.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valEntry2718, kTRUE);
    fLogbookHost = new TGTextEntry(fGroupFrame2717, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
    fLogbookHost->SetMaxLength(4096);
    fLogbookHost->SetAlignment(kTextLeft);
    fLogbookHost->SetText("fTextEntry1264");
    fLogbookHost->Resize(120,fLogbookHost->GetDefaultHeight());
    fGroupFrame2717->AddFrame(fLogbookHost, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookHost->MoveResize(150,20,120,20);
    fLogbookPort = new TGNumberEntry(fGroupFrame2717, (Double_t) 0,12,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 0,0,99999);
    fGroupFrame2717->AddFrame(fLogbookPort, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookPort->MoveResize(150,45,100,22);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valEntry2723;
    valEntry2723.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",valEntry2723.fForeground);
    gClient->GetColorByName("#e0e0e0",valEntry2723.fBackground);
    valEntry2723.fFillStyle = kFillSolid;
    valEntry2723.fFont = ufont->GetFontHandle();
    valEntry2723.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valEntry2723, kTRUE);
    fLogbookDatabase = new TGTextEntry(fGroupFrame2717, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
    fLogbookDatabase->SetMaxLength(4096);
    fLogbookDatabase->SetAlignment(kTextLeft);
    fLogbookDatabase->SetText("");
    fLogbookDatabase->Resize(120,fLogbookDatabase->GetDefaultHeight());
    fGroupFrame2717->AddFrame(fLogbookDatabase, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookDatabase->MoveResize(150,70,120,20);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valEntry2724;
    valEntry2724.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",valEntry2724.fForeground);
    gClient->GetColorByName("#e0e0e0",valEntry2724.fBackground);
    valEntry2724.fFillStyle = kFillSolid;
    valEntry2724.fFont = ufont->GetFontHandle();
    valEntry2724.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valEntry2724, kTRUE);
    fLogbookUser = new TGTextEntry(fGroupFrame2717, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
    fLogbookUser->SetMaxLength(4096);
    fLogbookUser->SetAlignment(kTextLeft);
    fLogbookUser->SetText("fTextEntry1332");
    fLogbookUser->Resize(120,fLogbookUser->GetDefaultHeight());
    fGroupFrame2717->AddFrame(fLogbookUser, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookUser->MoveResize(150,95,120,20);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valEntry2725;
    valEntry2725.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",valEntry2725.fForeground);
    gClient->GetColorByName("#e0e0e0",valEntry2725.fBackground);
    valEntry2725.fFillStyle = kFillSolid;
    valEntry2725.fFont = ufont->GetFontHandle();
    valEntry2725.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valEntry2725, kTRUE);
    fLogbookPassword = new TGTextEntry(fGroupFrame2717, new TGTextBuffer(14),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kOwnBackground);
    fLogbookPassword->SetMaxLength(4096);
    fLogbookPassword->SetAlignment(kTextLeft);
    fLogbookPassword->SetText("fTextEntry1343");
    fLogbookPassword->Resize(120,fLogbookPassword->GetDefaultHeight());
    fGroupFrame2717->AddFrame(fLogbookPassword, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookPassword->MoveResize(150,120,120,20);
    TGLabel *fLabel2726 = new TGLabel(fGroupFrame2717,"host");
    fLabel2726->SetTextJustify(kTextLeft);
    fLabel2726->SetMargins(0,0,0,0);
    fLabel2726->SetWrapLength(-1);
    fGroupFrame2717->AddFrame(fLabel2726, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2726->MoveResize(20,20,53,16);
    TGLabel *fLabel2727 = new TGLabel(fGroupFrame2717,"port");
    fLabel2727->SetTextJustify(kTextLeft);
    fLabel2727->SetMargins(0,0,0,0);
    fLabel2727->SetWrapLength(-1);
    fGroupFrame2717->AddFrame(fLabel2727, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2727->MoveResize(20,45,53,16);
    TGLabel *fLabel2728 = new TGLabel(fGroupFrame2717,"database");
    fLabel2728->SetTextJustify(kTextLeft);
    fLabel2728->SetMargins(0,0,0,0);
    fLabel2728->SetWrapLength(-1);
    fGroupFrame2717->AddFrame(fLabel2728, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2728->MoveResize(20,70,53,16);
    TGLabel *fLabel2729 = new TGLabel(fGroupFrame2717,"user");
    fLabel2729->SetTextJustify(kTextLeft);
    fLabel2729->SetMargins(0,0,0,0);
    fLabel2729->SetWrapLength(-1);
    fGroupFrame2717->AddFrame(fLabel2729, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2729->MoveResize(20,95,53,16);
    TGLabel *fLabel2730 = new TGLabel(fGroupFrame2717,"password");
    fLabel2730->SetTextJustify(kTextLeft);
    fLabel2730->SetMargins(0,0,0,0);
    fLabel2730->SetWrapLength(-1);
    fGroupFrame2717->AddFrame(fLabel2730, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel2730->MoveResize(20,120,53,16);
    
    fGroupFrame2717->SetLayoutManager(new TGVerticalLayout(fGroupFrame2717));
    fGroupFrame2717->Resize(320,175);
    fCompositeFrame2709->AddFrame(fGroupFrame2717, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2717->MoveResize(5,180,320,175);
    
    
    fTab2632->SetTab(0);
    
    fTab2632->Resize(fTab2632->GetDefaultSize());
    fMainFrame2342->AddFrame(fTab2632, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTab2632->MoveResize(0,0,570,526);
    
    fMainFrame2391->AddFrame(fMainFrame2342, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fMainFrame2342->MoveResize(0,0,580,530);
    fSaveAndExitButton = new TGTextButton(this,"Save and exit",1);
    fSaveAndExitButton->SetTextJustify(36);
    fSaveAndExitButton->SetMargins(0,0,0,0);
    fSaveAndExitButton->SetWrapLength(-1);
    fSaveAndExitButton->Resize(85,22);
    this->AddFrame(fSaveAndExitButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fSaveAndExitButton->MoveResize(480,536,85,22);
    fCancel = new TGTextButton(this,"Cancel",0,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
    fCancel->SetTextJustify(36);
    fCancel->SetMargins(0,0,0,0);
    fCancel->SetWrapLength(-1);
    fCancel->Resize(48,22);
    this->AddFrame(fCancel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fCancel->MoveResize(424,536,48,22);
    
    this->AddFrame(fMainFrame2391, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fMainFrame2391->MoveResize(0,0,570,560);
    
    this->SetMWMHints(kMWMDecorAll,kMWMFuncAll,kMWMInputModeless);
    this->MapSubwindows();
    
    this->Resize(this->GetDefaultSize());
    this->MapWindow();
    this->MoveResize(10,10,574,564);
}