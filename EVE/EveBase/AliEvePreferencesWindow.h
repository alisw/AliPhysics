//
//  AliEvePreferencesWindow.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 26/06/15.
//
//

#ifndef __AliEvePreferencesWindow__
#define __AliEvePreferencesWindow__

#include <TGNumberEntry.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TQObject.h>

class AliEvePreferencesWindow : public TGMainFrame
{
public:
    static AliEvePreferencesWindow* Instance();
    
    void onExit(bool save=true);
    
private:
    AliEvePreferencesWindow();
    ~AliEvePreferencesWindow();
    
    static AliEvePreferencesWindow *fInstance;
    
    void InitWindow();
    void ReadFromConfigFile();
    void SaveToConfigFile();
    void ApplyChanges();
    
    virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
    
    TGNumberEntry *fTrackWidth;
    TGCheckButton *fDashNoRefit;
    TGCheckButton *fDrawNoRefit;
    TGCheckButton *fTracksByPID;
    TGCheckButton *fTracksByCategory;
    
    TGCheckButton *fShowV0s;
    TGCheckButton *fShowCascades;
    TGCheckButton *fShowRawData;
    TGCheckButton *fShowPrimaryVertex;
    TGCheckButton *fShowHits;
    TGCheckButton *fShowDigits;
    TGCheckButton *fShowClusters;
    TGCheckButton *fShowKinks;
    
    TGTextEntry *fLogbookHost;
    TGNumberEntry *fLogbookPort;
    TGTextEntry *fLogbookDatabase;
    TGTextEntry *fLogbookUser;
    TGTextEntry *fLogbookPassword;
    
    TGCheckButton *fShowMuon;
    TGCheckButton *fShowHLTESDTree;
    TGTextEntry *fOCDBpath;
    TGCheckButton *fAutoload;
    TGCheckButton *fAliceLive;
    
    TGTextButton *fSaveAndExitButton;
    TGTextButton *fCancel;
    
    
    AliEvePreferencesWindow(const AliEvePreferencesWindow&);            // Not implemented
    AliEvePreferencesWindow& operator=(const AliEvePreferencesWindow&); // Not implemented
    
    ClassDef(AliEvePreferencesWindow, 0); // Short description.
};


#endif
