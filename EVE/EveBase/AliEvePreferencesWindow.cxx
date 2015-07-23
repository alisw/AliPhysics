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
    
    settings.WriteFile(Form("%s/eve_config",gSystem->Getenv("HOME")), kEnvAll);
}

void AliEvePreferencesWindow::ApplyChanges()
{
    TEnv settings;
    AliEveInit::GetConfig(&settings);
 
    AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)AliEveEventManager::GetMaster()->GetDataSourceOffline();
    const Text_t* esdfile = 0;
    
    if(settings.GetValue("HLTESDtree.show", false)){
        dataSource->SetESDFileName(esdfile, AliEveDataSourceOffline::kHLTTree);
    }
    else{
        dataSource->SetESDFileName(esdfile, AliEveDataSourceOffline::kOfflineTree);
    }
 
    AliEveInit::AddMacros();
 
    AliEveEventManager *man =  AliEveEventManager::GetMaster();
 
    man->SetESDwidth(settings.GetValue("tracks.width",2));
    man->SetESDdashNoRefit(settings.GetValue("tracks.noRefit.dash",true));
    man->SetESDdrawNoRefit(settings.GetValue("tracks.noRefit.show",true));
    man->SetESDtracksByCategory(settings.GetValue("tracks.byCategory.show",false));
    man->SetESDtracksByType(settings.GetValue("tracks.byType.show",true));
    man->SetSaveViews(settings.GetValue("ALICE_LIVE.send",false));
    man->SetAutoLoad(settings.GetValue("events.autoload.set",false));

//    AliEveEventManager::SetCdbUri(settings.GetValue("OCDB.default.path","local://$ALICE_ROOT/../src/OCDB"));
    
    man->InitOCDB(man->GetCurrentRun());
    
//    man->GetCurrentDataSource()->GotoEvent(man->GetEventId()); // reload event
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
    // "Tracks settings" group frame
    TGGroupFrame *fTracksSettingsGroup = new TGGroupFrame(this,"Tracks settings");
    fTracksSettingsGroup->SetLayoutBroken(kTRUE);
    
    TGLabel *fTrackWidthLabel = new TGLabel(fTracksSettingsGroup,"width");
    fTrackWidthLabel->SetTextJustify(36);
    fTrackWidthLabel->SetMargins(0,0,0,0);
    fTrackWidthLabel->SetWrapLength(-1);
    fTracksSettingsGroup->AddFrame(fTrackWidthLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTrackWidthLabel->MoveResize(8,21,72,20);
    
    fTrackWidth = new TGNumberEntry(fTracksSettingsGroup, (Double_t) 0,8,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 0,1,5);
    fTracksSettingsGroup->AddFrame(fTrackWidth, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTrackWidth->MoveResize(90,20,80,22);
    
    fDashNoRefit = new TGCheckButton(fTracksSettingsGroup,"dash no-refit tracks");
    fDashNoRefit->SetTextJustify(36);
    fDashNoRefit->SetMargins(0,0,0,0);
    fDashNoRefit->SetWrapLength(-1);
    fTracksSettingsGroup->AddFrame(fDashNoRefit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fDashNoRefit->MoveResize(8,45,150,20);
    
    fDrawNoRefit = new TGCheckButton(fTracksSettingsGroup,"draw no-refit tracks");
    fDrawNoRefit->SetTextJustify(36);
    fDrawNoRefit->SetMargins(0,0,0,0);
    fDrawNoRefit->SetWrapLength(-1);
    fTracksSettingsGroup->AddFrame(fDrawNoRefit, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fDrawNoRefit->MoveResize(8,65,150,20);
    
    fTracksByPID = new TGCheckButton(fTracksSettingsGroup,"tracks by PID");
    fTracksByPID->SetTextJustify(36);
    fTracksByPID->SetMargins(0,0,0,0);
    fTracksByPID->SetWrapLength(-1);
    fTracksSettingsGroup->AddFrame(fTracksByPID, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTracksByPID->MoveResize(8,85,150,20);
    
    fTracksByCategory = new TGCheckButton(fTracksSettingsGroup,"tracks by category");
    fTracksByCategory->SetTextJustify(36);
    fTracksByCategory->SetMargins(0,0,0,0);
    fTracksByCategory->SetWrapLength(-1);
    fTracksSettingsGroup->AddFrame(fTracksByCategory, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTracksByCategory->MoveResize(8,105,150,20);
    
    fTracksSettingsGroup->SetLayoutManager(new TGVerticalLayout(fTracksSettingsGroup));
    fTracksSettingsGroup->Resize(180,140);
    AddFrame(fTracksSettingsGroup, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTracksSettingsGroup->MoveResize(8,8,180,140);
    
    // "Elements to be shown" group frame
    TGGroupFrame *fElementsToShowGroup = new TGGroupFrame(this,"Elements to be shown");
    fElementsToShowGroup->SetLayoutBroken(kTRUE);
    
    fShowV0s = new TGCheckButton(fElementsToShowGroup,"V0s");
    fShowV0s->SetTextJustify(36);
    fShowV0s->SetMargins(0,0,0,0);
    fShowV0s->SetWrapLength(-1);
    
    fElementsToShowGroup->AddFrame(fShowV0s, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowV0s->MoveResize(8,60,100,20);
    
    fShowCascades = new TGCheckButton(fElementsToShowGroup,"cascades");
    fShowCascades->SetTextJustify(36);
    fShowCascades->SetMargins(0,0,0,0);
    fShowCascades->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowCascades, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowCascades->MoveResize(8,80,100,20);
    
    fShowRawData = new TGCheckButton(fElementsToShowGroup,"RAW data");
    fShowRawData->SetTextJustify(36);
    fShowRawData->SetMargins(0,0,0,0);
    fShowRawData->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowRawData, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowRawData->MoveResize(8,100,100,20);
    
    fShowPrimaryVertex = new TGCheckButton(fElementsToShowGroup,"primary vertex");
    fShowPrimaryVertex->SetTextJustify(36);
    fShowPrimaryVertex->SetMargins(0,0,0,0);
    fShowPrimaryVertex->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowPrimaryVertex, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowPrimaryVertex->MoveResize(8,120,100,20);
    
    fShowHits = new TGCheckButton(fElementsToShowGroup,"hits");
    fShowHits->SetTextJustify(36);
    fShowHits->SetMargins(0,0,0,0);
    fShowHits->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowHits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowHits->MoveResize(8,140,100,20);
    
    fShowDigits = new TGCheckButton(fElementsToShowGroup,"digits");
    fShowDigits->SetTextJustify(36);
    fShowDigits->SetMargins(0,0,0,0);
    fShowDigits->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowDigits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowDigits->MoveResize(8,160,100,20);
    
    fShowClusters = new TGCheckButton(fElementsToShowGroup,"clusters");
    fShowClusters->SetTextJustify(36);
    fShowClusters->SetMargins(0,0,0,0);
    fShowClusters->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowClusters, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowClusters->MoveResize(8,20,100,20);
    
    fShowKinks = new TGCheckButton(fElementsToShowGroup,"kinks");
    fShowKinks->SetTextJustify(36);
    fShowKinks->SetMargins(0,0,0,0);
    fShowKinks->SetWrapLength(-1);
    fElementsToShowGroup->AddFrame(fShowKinks, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowKinks->MoveResize(8,40,100,20);
    
    fElementsToShowGroup->SetLayoutManager(new TGVerticalLayout(fElementsToShowGroup));
    fElementsToShowGroup->Resize(180,200);
    AddFrame(fElementsToShowGroup, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fElementsToShowGroup->MoveResize(8,145,180,200);
    
    // "Logbook credentials" group frame
    TGGroupFrame *fLogbookCredentialsGroup = new TGGroupFrame(this,"Logbook credentials");
    fLogbookCredentialsGroup->SetLayoutBroken(kTRUE);
    
    fLogbookHost = new TGTextEntry(fLogbookCredentialsGroup, new TGTextBuffer(14),-1);
    fLogbookHost->SetMaxLength(4096);
    fLogbookHost->SetAlignment(kTextLeft);
    fLogbookHost->SetText("fTextEntry1264");
    fLogbookHost->Resize(120,fLogbookHost->GetDefaultHeight());
    fLogbookCredentialsGroup->AddFrame(fLogbookHost, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookHost->MoveResize(80,20,120,20);
    
    fLogbookPort = new TGNumberEntry(fLogbookCredentialsGroup, (Double_t) 0,8,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 0,0,99999);
    fLogbookCredentialsGroup->AddFrame(fLogbookPort, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookPort->MoveResize(80,40,100,22);
    
    fLogbookDatabase = new TGTextEntry(fLogbookCredentialsGroup, new TGTextBuffer(14),-1);
    fLogbookDatabase->SetMaxLength(4096);
    fLogbookDatabase->SetAlignment(kTextLeft);
    fLogbookDatabase->SetText("");
    fLogbookDatabase->Resize(120,fLogbookDatabase->GetDefaultHeight());
    fLogbookCredentialsGroup->AddFrame(fLogbookDatabase, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookDatabase->MoveResize(80,63,120,20);
    
    fLogbookUser = new TGTextEntry(fLogbookCredentialsGroup, new TGTextBuffer(14),-1);
    fLogbookUser->SetMaxLength(4096);
    fLogbookUser->SetAlignment(kTextLeft);
    fLogbookUser->SetText("fTextEntry1332");
    fLogbookUser->Resize(120,fLogbookUser->GetDefaultHeight());
    fLogbookCredentialsGroup->AddFrame(fLogbookUser, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookUser->MoveResize(80,83,120,20);
    
    fLogbookPassword = new TGTextEntry(fLogbookCredentialsGroup, new TGTextBuffer(14),-1);
    fLogbookPassword->SetMaxLength(4096);
    fLogbookPassword->SetAlignment(kTextLeft);
    fLogbookPassword->SetText("fTextEntry1343");
    fLogbookPassword->Resize(120,fLogbookPassword->GetDefaultHeight());
    fLogbookCredentialsGroup->AddFrame(fLogbookPassword, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookPassword->MoveResize(80,103,120,20);
    
    TGLabel *fLogbookHostLabel = new TGLabel(fLogbookCredentialsGroup,"host");
    fLogbookHostLabel->SetTextJustify(36);
    fLogbookHostLabel->SetMargins(0,0,0,0);
    fLogbookHostLabel->SetWrapLength(-1);
    fLogbookCredentialsGroup->AddFrame(fLogbookHostLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookHostLabel->MoveResize(8,20,60,20);
    
    TGLabel *fLogbookPortLabel = new TGLabel(fLogbookCredentialsGroup,"port");
    fLogbookPortLabel->SetTextJustify(36);
    fLogbookPortLabel->SetMargins(0,0,0,0);
    fLogbookPortLabel->SetWrapLength(-1);
    fLogbookCredentialsGroup->AddFrame(fLogbookPortLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookPortLabel->MoveResize(8,42,60,20);
    
    TGLabel *fLogbookDatabaseLabel = new TGLabel(fLogbookCredentialsGroup,"database");
    fLogbookDatabaseLabel->SetTextJustify(36);
    fLogbookDatabaseLabel->SetMargins(0,0,0,0);
    fLogbookDatabaseLabel->SetWrapLength(-1);
    fLogbookCredentialsGroup->AddFrame(fLogbookDatabaseLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookDatabaseLabel->MoveResize(8,62,60,20);
    
    TGLabel *fLogbookUserLabel = new TGLabel(fLogbookCredentialsGroup,"user");
    fLogbookUserLabel->SetTextJustify(36);
    fLogbookUserLabel->SetMargins(0,0,0,0);
    fLogbookUserLabel->SetWrapLength(-1);
    fLogbookCredentialsGroup->AddFrame(fLogbookUserLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookUserLabel->MoveResize(8,82,60,20);
    
    TGLabel *fLogbookPasswordLabel = new TGLabel(fLogbookCredentialsGroup,"password");
    fLogbookPasswordLabel->SetTextJustify(36);
    fLogbookPasswordLabel->SetMargins(0,0,0,0);
    fLogbookPasswordLabel->SetWrapLength(-1);
    fLogbookCredentialsGroup->AddFrame(fLogbookPasswordLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookPasswordLabel->MoveResize(8,102,60,20);
    
    fLogbookCredentialsGroup->SetLayoutManager(new TGVerticalLayout(fLogbookCredentialsGroup));
    fLogbookCredentialsGroup->Resize(260,140);
    AddFrame(fLogbookCredentialsGroup, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLogbookCredentialsGroup->MoveResize(200,205,260,140);
    
    // "Miscellaneous" group frame
    TGGroupFrame *fGeometrySettingsGroup = new TGGroupFrame(this,"Miscellaneous");
    fGeometrySettingsGroup->SetLayoutBroken(kTRUE);
    
    fShowMuon = new TGCheckButton(fGeometrySettingsGroup,"show MUON geometry and tracks");
    fShowMuon->SetTextJustify(36);
    fShowMuon->SetMargins(0,0,0,0);
    fShowMuon->SetWrapLength(-1);
    fGeometrySettingsGroup->AddFrame(fShowMuon, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowMuon->MoveResize(8,21,220,20);
    
    fShowHLTESDTree = new TGCheckButton(fGeometrySettingsGroup,"show HLT ESD tree");
    fShowHLTESDTree->SetTextJustify(36);
    fShowHLTESDTree->SetMargins(0,0,0,0);
    fShowHLTESDTree->SetWrapLength(-1);
    fGeometrySettingsGroup->AddFrame(fShowHLTESDTree, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fShowHLTESDTree->MoveResize(8,40,190,20);
    
    fOCDBpath = new TGTextEntry(fGeometrySettingsGroup, new TGTextBuffer(14),-1);
    fOCDBpath->SetMaxLength(4096);
    fOCDBpath->SetAlignment(kTextLeft);
    fOCDBpath->SetText("fTextEntry1461");
    fOCDBpath->Resize(240,fOCDBpath->GetDefaultHeight());
    fGeometrySettingsGroup->AddFrame(fOCDBpath, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fOCDBpath->MoveResize(8,120,240,22);
    
    TGLabel *fOCDBpathLabel = new TGLabel(fGeometrySettingsGroup,"OCDB path:");
    fOCDBpathLabel->SetTextJustify(36);
    fOCDBpathLabel->SetMargins(0,0,0,0);
    fOCDBpathLabel->SetWrapLength(-1);
    fGeometrySettingsGroup->AddFrame(fOCDBpathLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fOCDBpathLabel->MoveResize(8,100,80,20);
    
    fAutoload = new TGCheckButton(fGeometrySettingsGroup,"autoload new events");
    fAutoload->SetTextJustify(36);
    fAutoload->SetMargins(0,0,0,0);
    fAutoload->SetWrapLength(-1);
    fGeometrySettingsGroup->AddFrame(fAutoload, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fAutoload->MoveResize(8,60,190,20);
    
    fAliceLive = new TGCheckButton(fGeometrySettingsGroup,"send pictures to ALICE LIVE");
    fAliceLive->SetTextJustify(36);
    fAliceLive->SetMargins(0,0,0,0);
    fAliceLive->SetWrapLength(-1);
    fGeometrySettingsGroup->AddFrame(fAliceLive, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fAliceLive->MoveResize(8,80,190,20);
    
    fGeometrySettingsGroup->SetLayoutManager(new TGVerticalLayout(fGeometrySettingsGroup));
    fGeometrySettingsGroup->Resize(260,200);
    AddFrame(fGeometrySettingsGroup, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGeometrySettingsGroup->MoveResize(200,8,260,200);
    
    fSaveAndExitButton = new TGTextButton(this,"Save and exit",1);
    fSaveAndExitButton->SetTextJustify(36);
    fSaveAndExitButton->SetMargins(0,0,0,0);
    fSaveAndExitButton->SetWrapLength(-1);
    fSaveAndExitButton->Resize(100,30);
    AddFrame(fSaveAndExitButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fSaveAndExitButton->MoveResize(360,350,100,30);
    
    fCancel = new TGTextButton(this,"Cancel",0);
    fCancel->SetTextJustify(36);
    fCancel->SetMargins(0,0,0,0);
    fCancel->SetWrapLength(-1);
    fCancel->Resize(100,30);
    AddFrame(fCancel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fCancel->MoveResize(259,350,100,30);
    
    SetMWMHints(kMWMDecorAll,kMWMFuncAll,kMWMInputModeless);
    MapSubwindows();
    
    Resize(GetDefaultSize());
    MapWindow();
    MoveResize(100,100,468,386);
}
