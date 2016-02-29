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
#include <AliEveESDKinks.h>
#include <AliEveESDCascades.h>
#include <AliEveESDV0s.h>
#include <AliEveESDMuonTracks.h>
#include <AliEveESDSPDTracklets.h>
#include <AliEveAODTracks.h>
#include <AliEveDataSource.h>
#include <AliEveMomentumHistograms.h>
#include <AliEvePrimaryVertex.h>
#include <AliEveKineTracks.h>

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
    enum EDataType { kRaw,kHits,kDigits,kClusters,kESD,kAOD };
    
    AliEveEventManager(EDataSource defaultDataSource=kSourceOffline);
    static AliEveEventManager* Instance();
    
    // getters for data from current data source:
    AliRunLoader*  GetRunLoader(){ return fCurrentData->fRunLoader;}
    AliRawReader*  GetRawReader(){ return fCurrentData->fRawReader;}
    TFile*         GetESDFile()  { return fCurrentData->fESDFile;  }
    TTree*         GetESDTree()  { return fCurrentData->fESDTree;  }
    AliESDEvent*   GetESD()      { return fCurrentData->fESD;      }
    AliESDfriend*  GetESDfriend(){ return fCurrentData->fESDfriend;}
    TFile*         GetAODFile()  { return fCurrentData->fAODFile;  }
    TTree*         GetAODTree()  { return fCurrentData->fAODTree;  }
    AliAODEvent*   GetAOD()      { return fCurrentData->fAOD;      }
    
    //static getters for drawing macros
    static Int_t  CurrentEventId();
    static Bool_t HasESD();
    
    static AliRunLoader* AssertRunLoader();
    static AliESDEvent*  AssertESD();
    static AliESDfriend* AssertESDfriend();
    static AliAODEvent*  AssertAOD();
    static AliRawReader* AssertRawReader();
    static AliMagF*      AssertMagField();
    static TGeoManager*  AssertGeometry();
    
    void AddElement(TEveElement *element, TEveElement *parent=0);
    void Redraw3D(){gEve->Redraw3D();}
    void EnableRedraw(){gEve->EnableRedraw();}
    void DisableRedraw(){gEve->DisableRedraw();}
    
    // autoload timer getters and setters
    Double_t      GetAutoLoadTime()        const { return fAutoLoadTime; }
    Bool_t        GetAutoLoad()            const { return fAutoLoad;     }
    bool          GetAutoLoadRunning()     const { return fAutoLoadTimerRunning;}
    
    void          SetAutoLoadTime(Float_t time){fAutoLoadTime = time;}
    void          SetAutoLoad(Bool_t autoLoad);
    void          StartAutoLoadTimer();
    void          StopAutoLoadTimer();
    
    // global and transient elements:
    void          RegisterTransient(TEveElement* element);
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
    
    void           SetEventId(int eventId)    { fEventId=eventId;}
    void           SetCurrentRun(int run){fCurrentRun = run;}
    void           SetTrigSel(Int_t trig);
    void           SetHasEvent(bool hasEvent){fHasEvent=hasEvent;}
    
    // other public methods:
    void                        ResetMagneticField(){fgMagField=0;}
    virtual void                AfterNewEventLoaded();
    AliEveMacroExecutor*        GetExecutor() const { return fExecutor; }
    AliEveEventSelector*        GetEventSelector() const { return fPEventSelector; }
    AliEveMomentumHistograms*   GetMomentumHistogramsDrawer(){return fMomentumHistogramsDrawer;}
    
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
    
    TEveElementList*  fTransients;      // Container for additional transient (per event) elements.
    
    AliEveMacroExecutor*        fExecutor;                  // Executor for std macros
    AliEveSaveViews*            fViewsSaver;                // views saver
    AliEveESDTracks*            fESDTracksDrawer;           // drawer of ESD tracks
    AliEveAODTracks*            fAODTracksDrawer;           // drawer of AOD tracks
    AliEveMomentumHistograms*   fMomentumHistogramsDrawer;  // drawer of momentum histograms
    AliEvePrimaryVertex*        fPrimaryVertexDrawer;       // drawer of primary vertex
    AliEveESDKinks*             fKinksDrawer;               // drawer of ESD kinks
    AliEveESDCascades*          fCascadesDrawer;            // drawer of ESD cascades
    AliEveESDV0s*               fV0sDrawer;                 // drawer of ESD v0s
    AliEveESDMuonTracks*        fMuonTracksDrawer;          // drawer of ESD muon tracks
    AliEveESDSPDTracklets*      fSPDTracklersDrawer;        // drawer of ESD SPD tracklets
    AliEveKineTracks*           fKineTracksDrawer;          // drawer of tracks from Kinematics.root
    AliEveEventSelector*        fPEventSelector;            // Event filter
    
    Bool_t    fgGRPLoaded;     // Global run parameters loaded?
    AliMagF*  fgMagField;      // Global pointer to magnetic field.
    
    ClassDef(AliEveEventManager, 0); // Interface for getting all event components in a uniform way.
};

#endif
