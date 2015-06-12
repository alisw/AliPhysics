// $Id: AliEveEventManager.h 64557 2013-10-16 20:03:08Z hristov $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveEventManager_H
#define AliEveEventManager_H

#include <TEveEventManager.h>
#include <TQObject.h>
#include <TObjArray.h>
#include <TThread.h>

#include <AliEventInfo.h>
#include <AliESDEvent.h>
#include <AliEveSaveViews.h>
#include <AliEveESDTracks.h>
#include "AliEveDataSource.h"
#include "AliEveDataSourceOffline.h"

class AliEveMacroExecutor;
class AliEveEventSelector;

class AliRunLoader;
class AliESDEvent;
class AliESDfriend;
class AliAODEvent;
class AliRawReader;

class AliGRPObject;
class AliMagF;

class TEveElement;
class TFile;
class TTree;
class TGeoManager;
class TString;
class TMap;

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

    AliEveEventManager(EDataSource defaultDataSource=kSourceOffline);
    static AliEveEventManager* GetMaster();
    
    enum EVisibleESDTrees{ kOfflineTree, kHLTTree };

    virtual Int_t GetMaxEventId(Bool_t refreshESD=kFALSE) const{
        return fDataSourceOffline?fDataSourceOffline->GetMaxEventId(refreshESD):-1;
    }

    void Timeout(); // * SIGNAL*

    Int_t         GetEventId()         const { return fEventId; }
    void          SetEventId(int eventId)    { fEventId=eventId;}
    AliRunLoader* GetRunLoader()       const { return fCurrentData->fRunLoader; }
    TFile*        GetESDFile()         const { return fCurrentData->fESDFile; }
    TTree*        GetESDTree()         const { return fCurrentData->fESDTree; }
    TTree*        GetHLTESDTree()      const { return fCurrentData->fHLTESDTree; }
    AliESDEvent*  GetESD()             const { return fCurrentData->fESD;     }
    AliESDfriend* GetESDfriend()       const { return fCurrentData->fESDfriend; }
    TFile*        GetAODFile()         const { return fCurrentData->fAODFile; }
    TTree*        GetAODTree()         const { return fCurrentData->fAODTree; }
    AliAODEvent*  GetAOD()             const { return fCurrentData->fAOD;     }
    AliEveEventSelector* GetEventSelector() const { return fPEventSelector; }

    static Int_t  CurrentEventId();
    static Bool_t HasRunLoader();
    static Bool_t HasESD();
    static Bool_t HasESDfriend();
    static Bool_t HasAOD();
    static Bool_t HasRawReader();

    //static getters for drawing macros
    static AliRunLoader* AssertRunLoader();
    static AliESDEvent*  AssertESD();
    static AliESDfriend* AssertESDfriend();
    static AliAODEvent*  AssertAOD();
    static AliRawReader* AssertRawReader();
    static AliMagF*      AssertMagField();
    static TGeoManager*  AssertGeometry();


    static void RegisterTransient(TEveElement* element);
    static void RegisterTransientList(TEveElement* element);

    Double_t      GetAutoLoadTime()        const { return fAutoLoadTime; }
    Bool_t        GetAutoLoad()            const { return fAutoLoad;     }
    bool          GetAutoLoadRunning()     const { return fAutoLoadTimerRunning;}
    
    void          SetAutoLoadTime(Float_t time){fAutoLoadTime = time;}
    void          SetAutoLoad(Bool_t autoLoad);
    void          StartAutoLoadTimer();
    void          StopAutoLoadTimer();
    void          SetTrigSel(Int_t trig);
    void          AutoLoadNextEvent();
    void          SetSaveViews(bool save){fSaveViews=save;}
    
    void          SetESDtracksByCategory(bool set){fDrawESDtracksByCategory=set;}
    void          SetESDtracksByType(bool set){fDrawESDtracksByType=set;}
    void          SetESDcolorsByCategory(Color_t colors[9]){fESDdrawer->SetColorsByCategory(colors);}
    void          SetESDcolorsByType(Color_t colors[15]){fESDdrawer->SetColorsByType(colors);}
    void          SetESDwidth(Width_t width){fESDdrawer->SetWidth(width);}
    void          SetESDdashNoRefit(bool dashNoRefit){fESDdrawer->SetDashNoRefit(dashNoRefit);}
    void          SetESDdrawNoRefit(bool drawNoRefit){fESDdrawer->SetDrawNoRefit(drawNoRefit);}
    
    Bool_t        IsEventAvailable()       const { return fHasEvent;     }
    Bool_t        InsertGlobal(const TString& tag, TEveElement* model);
    Bool_t        InsertGlobal(const TString& tag, TEveElement* model,
                               Bool_t replace, Bool_t update);
    TEveElement*  FindGlobal(const TString& tag);

    virtual void  AfterNewEventLoaded();
    void          NewEventDataLoaded();  // *SIGNAL*
    void          NewEventLoaded();      // *SIGNAL*
    void          NoEventLoaded();      // *SIGNAL*
    void          StorageManagerOk();    // *SIGNAL*
    void          StorageManagerDown();  // *SIGNAL*
    void          EventServerOk();    // *SIGNAL*
    void          EventServerDown();  // *SIGNAL*
    
    AliEveMacroExecutor* GetExecutor() const { return fExecutor; }

    void ChangeDataSource(EDataSource newSource);
    AliEveDataSource* GetDataSourceOnline(){return fDataSourceOnline;}
    AliEveDataSource* GetDataSourceOffline(){return fDataSourceOffline;}
    AliEveDataSource* GetDataSourceHLTZMQ(){return fDataSourceHLTZMQ;}
    
    void DestroyTransients();
    void ResetMagneticField(){fgMagField=0;}
    void SetHasEvent(bool hasEvent){fHasEvent=hasEvent;}
    
    void SetCurrentRun(int run){fCurrentRun = run;}
    int  GetCurrentRun(){return fCurrentRun;}
    
protected:
    virtual ~AliEveEventManager();
    void SetMaster(AliEveEventManager *master);
    
    
    Int_t         fEventId;		// Id of current event.
    int fCurrentRun;

    const AliEveData* fCurrentData; //current data struct from one of the data sources
    AliEveDataSource* fCurrentDataSource; //data source in use at the moment
    AliEveDataSource *fDataSourceOnline;
    AliEveDataSource *fDataSourceOffline;
    AliEveDataSource *fDataSourceHLTZMQ;
    
    //std::map<EDataSource, AliDataSource*> fDataSources; //list of registered data sources (HLT,File,Online)

    AliEventInfo	fEventInfo;		// Current Event Info

    Bool_t        fAutoLoad;              // Automatic loading of events (online)
    Float_t       fAutoLoadTime;          // Auto-load time in seconds
    TTimer       *fAutoLoadTimer;         // Timer for automatic event loading

    Bool_t        fHasEvent;              // Is an event available.

    TMap*         fGlobal;
    Bool_t        fGlobalReplace;         // Are global replace
    Bool_t        fGlobalUpdate;          // Are global updates

    AliEveMacroExecutor *fExecutor;       // Executor for std macros

    TEveElementList     *fTransients;     // Container for additional transient (per event) elements.
    TEveElementList     *fTransientLists; // Container for lists of transient (per event) elements.

    AliEveEventSelector* fPEventSelector; // Event filter

    static Bool_t        fgGRPLoaded;     // Global run parameters loaded?
    static AliMagF      *fgMagField;      // Global pointer to magnetic field.
    static Bool_t        fgUniformField;  // Track with uniform field.
    Bool_t fAutoLoadTimerRunning; // State of auto-load timer.

private:
    static AliEveEventManager* fgMaster;
    
    void InitInternals();
    static Bool_t InitGRP();

    AliEveSaveViews *fViewsSaver;
    AliEveESDTracks *fESDdrawer;

    bool fSaveViews;
    bool fDrawESDtracksByCategory;
    bool fDrawESDtracksByType;
    
    bool fFirstEvent;
    bool fCenterProjectionsAtPrimaryVertex;
    
    AliEveEventManager(const AliEveEventManager&);            // Not implemented
    AliEveEventManager& operator=(const AliEveEventManager&); // Not implemented
    
    ClassDef(AliEveEventManager, 0); // Interface for getting all event components in a uniform way.
};

#endif
