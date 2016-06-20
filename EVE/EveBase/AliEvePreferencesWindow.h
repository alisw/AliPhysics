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
#include <TGColorSelect.h>
#include <TGTab.h>

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
    
    TGColorSelect *fColorSelect2658;//tpc
    TGColorSelect *fColorSelect2668;//tof
    TGColorSelect *fColorSelect2669;//trd
    TGColorSelect *fColorSelect2670;//muon
    TGColorSelect *fColorSelect2671;//spd
    TGColorSelect *fColorSelect2672;//sdd
    TGColorSelect *fColorSelect2673;//ssd
    TGColorSelect *fColorSelect2674;//phos
    TGColorSelect *fColorSelect2675;//emcal
    TGColorSelect *fColorSelect2676;//hmpid
    
    TGColorSelect *fColorSelect2693;// electron
    TGColorSelect *fColorSelect2694;// muon
    TGColorSelect *fColorSelect2695;// pion
    TGColorSelect *fColorSelect2696;// kaon
    TGColorSelect *fColorSelect2697;// proton
    TGColorSelect *fColorSelect2698;// deuteron
    TGColorSelect *fColorSelect2699;// triton
    TGColorSelect *fColorSelect2700;// he3
    TGColorSelect *fColorSelect2701;// alpha
    TGColorSelect *fColorSelect2702;// photon
    TGColorSelect *fColorSelect2703;// pi0
    TGColorSelect *fColorSelect2704;// neutron
    TGColorSelect *fColorSelect2705;// kaon0
    TGColorSelect *fColorSelect2706;// elecon
    TGColorSelect *fColorSelect2707;// unknown
    
    
    AliEvePreferencesWindow(const AliEvePreferencesWindow&);            // Not implemented
    AliEvePreferencesWindow& operator=(const AliEvePreferencesWindow&); // Not implemented
    
    ClassDef(AliEvePreferencesWindow, 0); // Short description.
};


#endif
