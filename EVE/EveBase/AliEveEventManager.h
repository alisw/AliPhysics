// $Id: AliEveEventManager.h 64557 2013-10-16 20:03:08Z hristov $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveEventManager_H
#define AliEveEventManager_H

#include <AliEventInfo.h>
#include <AliESDEvent.h>
#include <AliEveSaveViews.h>
#include <AliEveESDTracks.h>
#include "AliEveDataSource.h"

#include <TEveEventManager.h>
#include <TQObject.h>

class AliEveMacroExecutor;
class AliEveEventSelector;
class AliMagF;

//==============================================================================
//
// AliEveEventManager
//
// Interface to ALICE event-data (RunLoader, ESD), magnetic field and
// geometry.
//

class AliEveEventManager : public TEveEventManager, public TQObject
{
public:
    enum EDataSource { kSourceHLT, kSourceOnline, kSourceOffline };
    enum EVisibleESDTrees{ kOfflineTree, kHLTTree };
    
    AliEveEventManager(EDataSource defaultDataSource=kSourceOffline);
    static AliEveEventManager* GetMaster();

    // getters for data from current data source:
//    AliRunLoader*        GetRunLoader()     const { return fCurrentData->fRunLoader; }
    TFile*               GetESDFile()       const { return fCurrentData->fESDFile; }
    TTree*               GetESDTree()       const { return fCurrentData->fESDTree; }
    TTree*               GetHLTESDTree()    const { return fCurrentData->fHLTESDTree; }
    AliESDEvent*         GetESD()           const { return fCurrentData->fESD;     }
//    AliESDfriend*        GetESDfriend()     const { return fCurrentData->fESDfriend; }
//    TFile*               GetAODFile()       const { return fCurrentData->fAODFile; }
//    TTree*               GetAODTree()       const { return fCurrentData->fAODTree; }
//    AliAODEvent*         GetAOD()           const { return fCurrentData->fAOD;     }

    //static getters for drawing macros
    static Int_t  CurrentEventId();
    static Bool_t HasRunLoader();
    static Bool_t HasESD();
    static Bool_t HasESDfriend();
    static Bool_t HasAOD();
    static Bool_t HasRawReader();
    
    static AliRunLoader* AssertRunLoader();
    static AliESDEvent*  AssertESD();
    static AliESDfriend* AssertESDfriend();
    static AliAODEvent*  AssertAOD();
    static AliRawReader* AssertRawReader();
    static AliMagF*      AssertMagField();
    static TGeoManager*  AssertGeometry();
    
    // autoload timer getters and setters
    Double_t      GetAutoLoadTime()        const { return fAutoLoadTime; }
    Bool_t        GetAutoLoad()            const { return fAutoLoad;     }
    bool          GetAutoLoadRunning()     const { return fAutoLoadTimerRunning;}
    
    void          SetAutoLoadTime(Float_t time){fAutoLoadTime = time;}
    void          SetAutoLoad(Bool_t autoLoad);
    void          StartAutoLoadTimer();
    void          StopAutoLoadTimer();
    
    // options' setters
    void          SetSaveViews(bool save){fSaveViews=save;}
    void          SetESDtracksByCategory(bool set){fDrawESDtracksByCategory=set;}
    void          SetESDtracksByType(bool set){fDrawESDtracksByType=set;}
    void          SetESDcolorsByCategory(Color_t colors[9]){fESDdrawer->SetColorsByCategory(colors);}
    void          SetESDwidth(Width_t width){fESDdrawer->SetWidth(width);}
    void          SetESDdashNoRefit(bool dashNoRefit){fESDdrawer->SetDashNoRefit(dashNoRefit);}
    void          SetESDdrawNoRefit(bool drawNoRefit){fESDdrawer->SetDrawNoRefit(drawNoRefit);}
    
    // global and transient elements:
    Bool_t        InsertGlobal(const TString& tag, TEveElement* model);
    Bool_t        InsertGlobal(const TString& tag, TEveElement* model,Bool_t replace, Bool_t update);
    TEveElement*  FindGlobal(const TString& tag);
    static void   RegisterTransient(TEveElement* element);
//    static void   RegisterTransientList(TEveElement* element);
    void          DestroyTransients();
    
    // data sources:
    void ChangeDataSource(EDataSource newSource);
    AliEveDataSource* GetCurrentDataSource(){return fCurrentDataSource;}
    AliEveDataSource* GetDataSourceOnline(){return fDataSourceOnline;}
    AliEveDataSource* GetDataSourceOffline(){return fDataSourceOffline;}
    AliEveDataSource* GetDataSourceHLTZMQ(){return fDataSourceHLTZMQ;}
    

    static void SetCdbUri(TString path);

    // getters and setters for info about events:
    Int_t          GetEventId() const {return fEventId;}
    virtual Int_t  GetMaxEventId(Bool_t refreshESD=kFALSE) const;
    int            GetCurrentRun() {return fCurrentRun;}
//    Bool_t         IsEventAvailable() const {return fHasEvent;}

    void           SetEventId(int eventId)    { fEventId=eventId;}
    void           SetCurrentRun(int run){fCurrentRun = run;}
    void           SetTrigSel(Int_t trig);
    void           SetHasEvent(bool hasEvent){fHasEvent=hasEvent;}
    
    // other public methods:
    void                  ResetMagneticField(){fgMagField=0;}
    virtual void          AfterNewEventLoaded();
    AliEveMacroExecutor*  GetExecutor() const { return fExecutor; }
    AliEveEventSelector*  GetEventSelector() const { return fPEventSelector; }

    Bool_t InitOCDB(int runNo);
    
    // signals:
    void Timeout();             // *SIGNAL*
    void NewEventDataLoaded();  // *SIGNAL*
    void NewEventLoaded();      // *SIGNAL*
    void NoEventLoaded();       // *SIGNAL*
    
    void AutoLoadNextEvent();
private:
    virtual ~AliEveEventManager();
    static AliEveEventManager* fgMaster; // singleton instance of AliEveEventManager
    
    void   InitInternals();
    Bool_t InitGRP();
    
    Int_t         fEventId;		// Id of current event.
    AliEventInfo  fEventInfo;   // Current Event Info
    Bool_t        fHasEvent;    // Is an event available.
    int           fCurrentRun;  // Current run number
    
    AliEveData        fEmptyData;          //just a place holder in case we have no sources
    const AliEveData* fCurrentData;        //current data struct from one of the data sources
    AliEveDataSource* fCurrentDataSource;  // data source in use at the moment
    AliEveDataSource* fDataSourceOnline;   // pointer to online data source
    AliEveDataSource* fDataSourceOffline;  // pointer to offline data source
    AliEveDataSource* fDataSourceHLTZMQ;   // pointer to HLT ZMQ data source
    

    Bool_t   fAutoLoad;              // Automatic loading of events (online)
    Float_t  fAutoLoadTime;          // Auto-load time in seconds
    TTimer*  fAutoLoadTimer;         // Timer for automatic event loading
    Bool_t   fAutoLoadTimerRunning;  // State of auto-load timer.

    
    TMap*             fGlobal;          // Map of global elements
    Bool_t            fGlobalReplace;   // Are global replace
    Bool_t            fGlobalUpdate;    // Are global updates
    TEveElementList*  fTransients;      // Container for additional transient (per event) elements.
    TEveElementList*  fTransientLists;  // Container for lists of transient (per event) elements.
    
    AliEveMacroExecutor*  fExecutor;        // Executor for std macros
    AliEveSaveViews*      fViewsSaver;      // views saver
    AliEveESDTracks*      fESDdrawer;       // drawer of ESD tracks
    AliEveEventSelector*  fPEventSelector;  // Event filter
    
    Bool_t    fgGRPLoaded;     // Global run parameters loaded?
    AliMagF*  fgMagField;      // Global pointer to magnetic field.

    
    bool fSaveViews;
    bool fDrawESDtracksByCategory;
    bool fDrawESDtracksByType;
    bool fFirstEvent;
    
    ClassDef(AliEveEventManager, 0); // Interface for getting all event components in a uniform way.
};

#endif
