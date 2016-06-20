// $Id: AliEveEventManager.cxx 64557 2013-10-16 20:03:08Z hristov $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveEventManager.h"
#include "AliEveEventSelector.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMultiView.h"
#include "AliEveDataSourceOffline.h"
#include "AliEveInit.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"

#ifdef ZMQ
#include "AliEveDataSourceOnline.h"
#include "AliEveDataSourceHLTZMQ.h"
#endif

#include "AliGRPPreprocessor.h"
#include <TEnv.h>

#include <AliGRPManager.h>
#include <AliLog.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include "AliCDBStorage.h"
#include <AliMagF.h>
#include <AliGeomManager.h>

#include <TEveElement.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TTimeStamp.h>
#include <TROOT.h>
#include <TEveText.h>
#include <TEveTrans.h>
#include <iostream>


using namespace std;

//==============================================================================
//==============================================================================
// AliEveEventManager
//==============================================================================

//______________________________________________________________________________
//
// Provides interface for loading and navigating standard AliRoot data
// (AliRunLoader), ESD, AOD and RAW.
//
// ESDfriend is attached automatically, if the file is found.
//
// AODfriends are not attached automatically as there are several
// possible files involved. To have a specific AODfriend attached, call
// static method
//   AliEveEventManager::AddAODfriend("AliAOD.VertexingHF.root");
// before initializing the event-manager.
//
// Also provides interface to magnetic-field and geometry. Mostly
// intended as wrappers over standard AliRoot functionality for
// convenient use from visualizateion macros.
//
// There can be a single main event-manger, it is stored in private
// data member fgMaster and can be accessed via static member function
// GetMaster().
//
// For event overlaying and embedding one can instantiate additional
// event-managers via static method AddDependentManager(const TString& path).
// This interface is under development.

ClassImp(AliEveEventManager)

AliEveEventManager* AliEveEventManager::fgMaster  = NULL;

AliEveEventManager::AliEveEventManager(EDataSource defaultDataSource) :
TEveEventManager("Event", ""),
fEventId(-1),fEventInfo(),fHasEvent(kFALSE),fCurrentRun(-1),
fCurrentData(&fEmptyData),fCurrentDataSource(NULL),fDataSourceOnline(NULL),fDataSourceOffline(NULL),fDataSourceHLTZMQ(NULL),
fAutoLoad(kFALSE), fAutoLoadTime(5),fAutoLoadTimer(0),fAutoLoadTimerRunning(kFALSE),
fTransients(0),
fExecutor(0),fViewsSaver(0),fESDTracksDrawer(0),fAODTracksDrawer(0),fPrimaryVertexDrawer(0),fKinksDrawer(0),fCascadesDrawer(0),fV0sDrawer(0),fMuonTracksDrawer(0),fSPDTracklersDrawer(0),fKineTracksDrawer(0),  fMomentumHistogramsDrawer(0),fPEventSelector(0),
fgGRPLoaded(false),
fgMagField(0)
{
    InitInternals();
    ChangeDataSource(defaultDataSource);
}

AliEveEventManager* AliEveEventManager::Instance()
{
    // Get master event-manager.
    if(fgMaster){return fgMaster;}
    else{cout<<"FATAL -- Event Manager was not created.\n"<<endl;exit(0);}
}

AliEveEventManager::~AliEveEventManager()
{
    // Destructor.
//    fAutoLoadTimer->Stop();
//    fAutoLoadTimer->Disconnect("Timeout");
//    fAutoLoadTimer->Disconnect("AutoLoadNextEvent");
}

void AliEveEventManager::InitInternals()
{
    // Initialize internal members.
    fgMaster = this;
    
    fAutoLoadTimer = new TTimer;
    fAutoLoadTimer->Connect("Timeout()", "AliEveEventManager", this, "AutoLoadNextEvent()");
    fAutoLoadTimer->Connect("Timeout()", "AliEveEventManager", this, "Timeout()");
    
    fExecutor = new AliEveMacroExecutor;
    
    fTransients = new TEveElementList("Transients", "Transient per-event elements.");
    fTransients->IncDenyDestroy();
    gEve->AddToListTree(fTransients, kFALSE);
    
    fPEventSelector = new AliEveEventSelector(this);
    
    fViewsSaver = new AliEveSaveViews();
    fESDTracksDrawer = new AliEveESDTracks();
    fAODTracksDrawer = new AliEveAODTracks();
    fMomentumHistogramsDrawer = new AliEveMomentumHistograms();
    fPrimaryVertexDrawer = new AliEvePrimaryVertex();
    fKinksDrawer = new AliEveESDKinks();
    fCascadesDrawer = new AliEveESDCascades();
    fV0sDrawer = new AliEveESDV0s();
    fMuonTracksDrawer = new AliEveESDMuonTracks();
    fSPDTracklersDrawer = new AliEveESDSPDTracklets();
    fKineTracksDrawer = new AliEveKineTracks();
    
#ifdef ZMQ
    fDataSourceOnline = new AliEveDataSourceOnline();
    fDataSourceHLTZMQ = new AliEveDataSourceHLTZMQ();    
#endif
    fDataSourceOffline = new AliEveDataSourceOffline();
}

void AliEveEventManager::ChangeDataSource(EDataSource newSource)
{
    //before switching stop autoload timer and process events
    if (fAutoLoadTimerRunning){
        StopAutoLoadTimer();
        gSystem->ProcessEvents();
    }

    fCurrentData = &fEmptyData;
    if(newSource == kSourceOnline)
    {
        fCurrentDataSource = fDataSourceOnline;
    }
    else if(newSource == kSourceOffline)
    {
        fCurrentDataSource = fDataSourceOffline;
    }
    else if(newSource == kSourceHLT)
    {
        fCurrentDataSource = fDataSourceHLTZMQ;
    }
    if (fCurrentDataSource) fCurrentData = fCurrentDataSource->GetData();

    //restore timer
    if (fAutoLoad)
    {
        StartAutoLoadTimer();
    }
}

void AliEveEventManager::DestroyTransients()
{
    TEveManager::TRedrawDisabler rd(gEve);
    gEve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.
    
    gEve->GetViewers()->DeleteAnnotations();
    fTransients->DestroyElements();
    DestroyElements();
    ElementChanged();
}

Int_t AliEveEventManager::GetMaxEventId(Bool_t refreshESD) const
{
    return fDataSourceOffline?fDataSourceOffline->GetMaxEventId(refreshESD):-1;
}

//------------------------------------------------------------------------------
// Static convenience functions, mainly used from macros.
//------------------------------------------------------------------------------

Int_t AliEveEventManager::CurrentEventId()
{
    // Return current event-id.
    
    static const TEveException kEH("AliEveEventManager::CurrentEventId ");
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE)
        throw (kEH + "ALICE event not ready.");
    return fgMaster->GetEventId();
}

Bool_t AliEveEventManager::HasESD()
{
    // Check if AliESDEvent is initialized.
    return fgMaster && fgMaster->fHasEvent && fgMaster->fCurrentData->fESD;
}

AliRunLoader* AliEveEventManager::AssertRunLoader()
{
    // Make sure AliRunLoader is initialized and return it.
    // Static utility for macros.
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE){
        cout<<"AliEveEventManager::AssertRunLoader : ALICE event not ready."<<endl;
        return 0;
    }
    else if (fgMaster->fCurrentData->fRunLoader == 0){
        cout<<"AliEveEventManager::AssertRunLoader : AliRunLoader not initialised."<<endl;
        return 0;
    }
    return fgMaster->fCurrentData->fRunLoader;
}

AliESDEvent* AliEveEventManager::AssertESD()
{
    // Make sure AliESDEvent is initialized and return it.
    // Throws exception in case ESD is not available.
    // Static utility for macros.
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE){
        cout<<"AliEveEventManager::AssertESD() - ALICE event not ready."<<endl;
        return 0;
    }
    if (fgMaster->fCurrentData->fESD == 0)
    {
        cout<<"AliEveEventManager::AssertESD() - AliESD not initialised."<<endl;
        return 0;
    }
    return fgMaster->fCurrentData->fESD;
}

AliESDfriend* AliEveEventManager::AssertESDfriend()
{
    // Make sure AliESDfriend is initialized and return it.
    // Throws exception in case ESDfriend-loader is not available.
    // Static utility for macros.
    
    static const TEveException kEH("AliEveEventManager::AssertESDfriend ");
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE)
        throw (kEH + "ALICE event not ready.");
    if (fgMaster->fCurrentData->fESDfriend == 0)
        throw (kEH + "AliESDfriend not initialised.");
    return fgMaster->fCurrentData->fESDfriend;
}

AliAODEvent* AliEveEventManager::AssertAOD()
{
    // Make sure AliAODEvent is initialized and return it.
    // Throws exception in case AOD is not available.
    // Static utility for macros.
    
    static const TEveException kEH("AliEveEventManager::AssertAOD ");
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE)
        throw (kEH + "ALICE event not ready.");
    if (fgMaster->fCurrentData->fAOD == 0)
        throw (kEH + "AliAOD not initialised.");
    return fgMaster->fCurrentData->fAOD;
}

AliRawReader* AliEveEventManager::AssertRawReader()
{
    return fgMaster->fCurrentData->fRawReader;
}

//==============================================================================

AliMagF* AliEveEventManager::AssertMagField()
{
    // Make sure AliMagF is initialized and returns it.
    // Throws exception in case magnetic field is not available.
    // Static utility for macros.

    static const TEveException kEH("AliEveEventManager::AssertMagField ");
    
    //if we already have a field we're done
    if (fgMaster->fgMagField) return fgMaster->fgMagField;
    
    //if not: init the mag field from the ESD information
    if (fgMaster->AssertESD()) fgMaster->AssertESD()->InitMagneticField();

    //check if we have field from ESD header
    if (TGeoGlobalMagField::Instance()->GetField())
    {
        fgMaster->fgMagField = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
        if (fgMaster->fgMagField) return fgMaster->fgMagField;
    }
    
    //if no field from ESD, try to init from GRP
    if (!fgMaster->fgGRPLoaded)
    {
        if (fgMaster->InitGRP()){
            fgMaster->fgGRPLoaded = kTRUE;
        }
    }
    
    //one last check
    if (!TGeoGlobalMagField::Instance()->GetField())
    {
        AliEveDataSource* dataSource = fgMaster->GetCurrentDataSource();
        if (dataSource) {
          fgMaster->SetCurrentRun(AliCDBManager::Instance()->GetRun());
          dataSource->ReceivePromptRecoParameters(AliCDBManager::Instance()->GetRun());
        }
        if (fgMaster->InitGRP()){
            fgMaster->fgGRPLoaded = kTRUE;
        }
    }
    
    //check if now we have some field from the GRP:
    if (TGeoGlobalMagField::Instance()->GetField())
    {
        fgMaster->fgMagField = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
        if (!fgMaster->fgMagField)
            throw kEH + "Global field set, but it is not AliMagF.";
    }
    else
    {
        throw kEH + "Could not initialize magnetic field.";
    }
    
    return fgMaster->fgMagField;
}

TGeoManager* AliEveEventManager::AssertGeometry()
{
    // Make sure AliGeomManager is initialized and returns the
    // corresponding TGeoManger.
    // gGeoManager is set to the return value.
    // Throws exception if geometry can not be loaded or if it is not
    // available and the TGeoManager is locked.
    // Static utility for macros.
    
    static const TEveException kEH("AliEveEventManager::AssertGeometry ");
    
    if (AliGeomManager::GetGeometry() == 0)
    {
        if (TGeoManager::IsLocked())
            throw (kEH + "geometry is not loaded but TGeoManager is locked.");
        
        gGeoManager = 0;
        AliGeomManager::LoadGeometry();
        if ( ! AliGeomManager::GetGeometry())
        {
            throw (kEH + "can not load geometry.");
        }
        if ( ! AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE"))
        {
            ::Warning(kEH, "mismatch of alignable volumes. Proceeding.");
            // throw (kEH + "could not apply align objs.");
        }
        AliGeomManager::GetGeometry()->DefaultColors();
    }
    
    gGeoManager = AliGeomManager::GetGeometry();
    return gGeoManager;
}

void AliEveEventManager::AddElement(TEveElement *element, TEveElement *parent)
{
    gEve->AddElement(element,parent);
}

void AliEveEventManager::RegisterTransient(TEveElement* element)
{
    fTransients->AddElement(element);
}

//------------------------------------------------------------------------------
// Autoloading of events
//------------------------------------------------------------------------------

void AliEveEventManager::SetAutoLoad(Bool_t autoLoad)
{
    // Set the automatic event loading mode
    static const TEveException kEH("AliEveEventManager::SetAutoLoad ");
    
    cout<<"\n\n setting autoload to:"<<autoLoad<<endl;
    
    if (fAutoLoad == autoLoad)
    {
        Warning(kEH, "Setting autoload to the same value as before - %s. Ignoring.", fAutoLoad ? "true" : "false");
        return;
    }
    
    fAutoLoad = autoLoad;
    if (fAutoLoad)
    {
        //        StorageManagerDown();
        StartAutoLoadTimer();
    }
    else
    {
        //        StorageManagerOk();
        StopAutoLoadTimer();
    }
}

void AliEveEventManager::SetTrigSel(Int_t trig)
{
    static const TEveException kEH("AliEveEventManager::SetTrigSel ");
    
    if (!fCurrentData->fRawReader)
    {
        Warning(kEH, "No Raw-reader exists. Ignoring the call.");
        return;
    }
    else
    {
        ULong64_t trigMask = 0;
        if (trig >= 0) trigMask = (1ull << trig);
        Info(kEH,"Trigger selection: 0x%llx",trigMask);
        fCurrentData->fRawReader->SelectEvents(-1,trigMask,NULL);
    }
}

void AliEveEventManager::StartAutoLoadTimer()
{
    // Start the auto-load timer.
    fAutoLoadTimer->SetTime((Long_t)(1000*fAutoLoadTime));
    fAutoLoadTimer->Reset();
    fAutoLoadTimer->TurnOn();
    fAutoLoadTimerRunning = kTRUE;
}

void AliEveEventManager::StopAutoLoadTimer()
{
    // Stop the auto-load timer.
    fAutoLoadTimerRunning = kFALSE;
    fAutoLoadTimer->TurnOff();
}

void AliEveEventManager::AutoLoadNextEvent()
{
    // Called from auto-load timer, so it has to be public.
    // Do NOT call it directly.
        
    static const TEveException kEH("AliEveEventManager::AutoLoadNextEvent ");
    
    Info(kEH, "called!");
    
    if ( ! fAutoLoadTimerRunning || ! fAutoLoadTimer->HasTimedOut())
    {
        Warning(kEH, "Called unexpectedly - ignoring the call. Should ONLY be called from an internal timer.");
        return;
        }
    
    StopAutoLoadTimer();
    cout<<"Calling NextEvent method on current data source"<<endl;

    fCurrentDataSource->NextEvent();
    
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    
    if(settings.GetValue("save.on.autoload",false))
    {
        if(fCurrentData->fESD->GetNumberOfTracks() != 0){
            fViewsSaver->Save(false,Form("alice_%d.png",fEventId));
        }
    }
    
    if (fAutoLoad){
        StartAutoLoadTimer();
    }
}

//------------------------------------------------------------------------------
// Post event-loading functions
//------------------------------------------------------------------------------

void AliEveEventManager::AfterNewEventLoaded()
{
    // Execute registered macros and commands.
    // At the end emit NewEventLoaded signal.
    //
    // Virtual from TEveEventManager.
  
    Bool_t cdbOK = kTRUE;
    if (fCurrentData->fESD) {
      cdbOK = InitOCDB(fCurrentData->fESD->GetRunNumber());
    }
    if (!cdbOK) {
      printf("CDB not OK! not executing AfterNewEventLoaded\n");
      return;
    }
    else if(fCurrentData->fRawReader){
        InitOCDB(fCurrentData->fRawReader->GetRunNumber());
    }
    AliEveMultiView *mv = AliEveMultiView::Instance();
    mv->DestroyAllEvents();
    

    cout<<"\n\n--------- AliEveEventManager::AfterNewEventLoaded ---------\n\n"<<endl;
    
    ElementChanged();
    NewEventDataLoaded();
    
//    cout<<"\n\n-------- Executing registered macros ---------\n\n"<<endl;
    if (fExecutor){
        fExecutor->ExecMacros();
    }
    else{
        cout<<"no macro executor..."<<endl;
    }
//    cout<<"\n\n-------- Macros have been executed ---------\n\n"<<endl;
    
    TEveEventManager::AfterNewEventLoaded();
    NewEventLoaded();

    TEnv settings;
    AliEveInit::GetConfig(&settings);
    
    if(fCurrentData->fESD && fHasEvent)
    {
        if(settings.GetValue("momentum.histograms.show",false)){fMomentumHistogramsDrawer->Draw();}
        if(settings.GetValue("tracks.primary.vertex.show",false)){fESDTracksDrawer->PrimaryVertexTracks();}
        
        if(settings.GetValue("primary.vertex.global.box",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kGlobal,AliEvePrimaryVertex::kBox);
        }
        if(settings.GetValue("primary.vertex.global.ellipse",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kGlobal,AliEvePrimaryVertex::kEllipse);
        }
        if(settings.GetValue("primary.vertex.global.cross",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kGlobal,AliEvePrimaryVertex::kCross);
        }
        if(settings.GetValue("primary.vertex.tpc.box",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kTPC,AliEvePrimaryVertex::kBox);
        }
        if(settings.GetValue("primary.vertex.tpc.ellipse",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kTPC,AliEvePrimaryVertex::kEllipse);
        }
        if(settings.GetValue("primary.vertex.tpc.cross",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kTPC,AliEvePrimaryVertex::kCross);
        }
        if(settings.GetValue("primary.vertex.spd.box",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kSPD,AliEvePrimaryVertex::kBox);
        }
        if(settings.GetValue("primary.vertex.spd.ellipse",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kSPD,AliEvePrimaryVertex::kEllipse);
        }
        if(settings.GetValue("primary.vertex.spd.cross",false)){
            fPrimaryVertexDrawer->PrimaryVertex(AliEvePrimaryVertex::kSPD,AliEvePrimaryVertex::kCross);
        }
        
        if(settings.GetValue("tracks.byCategory.show",false))fESDTracksDrawer->ByCategory();
        if(settings.GetValue("tracks.byType.show",true))fESDTracksDrawer->ByType();
        if(settings.GetValue("tracks.byPt.show",false))fESDTracksDrawer->ByPt();
        if(settings.GetValue("kinks.show",false)){fKinksDrawer->Draw();}
        if(settings.GetValue("kinks.points.show",false)){fKinksDrawer->DrawPoints();}
        if(settings.GetValue("cascades.show",false)){fCascadesDrawer->Draw();}
        if(settings.GetValue("cascades.points.show",false)){fCascadesDrawer->DrawPoints();}
        if(settings.GetValue("V0s.show",false)){fV0sDrawer->Draw();}
        if(settings.GetValue("V0s.points.offline.show",false)){fV0sDrawer->DrawPointsOffline();}
        if(settings.GetValue("V0s.points.onfly.show",false)){fV0sDrawer->DrawPointsOnfly();}
        if(settings.GetValue("MUON.show",true)){fMuonTracksDrawer->Draw();}
        if(settings.GetValue("SPD.tracklets.show",false)){fSPDTracklersDrawer->Draw();}
        if(settings.GetValue("kine.tracks.show",false)){fKineTracksDrawer->Draw();}
        if(settings.GetValue("HLT.tracks.show",false)){fESDTracksDrawer->HLTTracks();}
        
        Double_t x[3] = { 0, 0, 0 };
        
        AliESDEvent *esd = fCurrentData->fESD;
        esd->GetPrimaryVertex()->GetXYZ(x);
        
        TTimeStamp ts(esd->GetTimeStamp());
        TString win_title("Eve Main Window -- Timestamp: ");
        win_title += ts.AsString("s");
        win_title += "; Event: ";
        win_title += esd->GetEventNumberInFile();
        win_title += "; Run: ";
        win_title += esd->GetRunNumber();
        gEve->GetBrowser()->SetWindowName(win_title);
    }
    if(fCurrentData->fAOD)
    {
        if(settings.GetValue("tracks.aod.byPID.show",false))fAODTracksDrawer->ByPID();
        
        Double_t x[3] = { 0, 0, 0 };
        
        AliAODEvent *aod = fCurrentData->fAOD;
        
        aod->GetPrimaryVertex()->GetXYZ(x);
        
        TTimeStamp ts(aod->GetTimeStamp());
        TString win_title("Eve Main Window -- Timestamp: ");
        win_title += ts.AsString("s");
        win_title += "; Event: ";
        win_title += aod->GetEventNumberInFile();
        win_title += "; Run: ";
        win_title += aod->GetRunNumber();
        gEve->GetBrowser()->SetWindowName(win_title);
    }
    
    TEveElement* top = gEve->GetCurrentEvent();
    mv->ImportEvent(top);
    
    gEve->GetBrowser()->RaiseWindow();
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    
    bool saveViews = settings.GetValue("ALICE_LIVE.send",false);
    
    if(saveViews  && fCurrentData->fESD->GetNumberOfTracks()>0)
    {
        fViewsSaver->SaveForAmore();
        fViewsSaver->SendToAmore();
    }
    
    // animate tracks
    if(settings.GetValue("tracks.byType.animate",false)){
        TEveElementList *tracks = (TEveElementList*)gEve->GetCurrentEvent()->FindChild("ESD Tracks by PID");
        
        vector<TEveTrackPropagator*> propagators;
        
        for(TEveElement::List_i i = tracks->BeginChildren();i != tracks->EndChildren();i++)
        {
            TEveTrackList* trkList = ((TEveTrackList*)*i);
            propagators.push_back(trkList->GetPropagator());
        }
        int animationSpeed = settings.GetValue("tracks.animation.speed",10);
        for(int R=0; R<520; R+=animationSpeed)
        {
            animationSpeed-=10;
            if(animationSpeed<10)animationSpeed=10;
            for(int propIter=0;propIter<propagators.size();propIter++){
                propagators[propIter]->SetMaxR(R);
            }
            gSystem->ProcessEvents();
            gEve->Redraw3D();
        }
    }
}

void AliEveEventManager::Timeout()
{
    Emit("Timeout()");
}
void AliEveEventManager::NewEventDataLoaded()
{
    Emit("NewEventDataLoaded()");
}
void AliEveEventManager::NewEventLoaded()
{
    Emit("NewEventLoaded()");
}
void AliEveEventManager::NoEventLoaded()
{
    Emit("NoEventLoaded()");
}

Bool_t AliEveEventManager::InitGRP()
{
    // Initialization of the GRP entry
    
    static const TEveException kEH("AliEveEventManager::InitGRP ");
    
    AliGRPManager grpMgr;
    if (!grpMgr.ReadGRPEntry()) {
        return kFALSE;
    }
    fgGRPLoaded = kTRUE;
    if (!grpMgr.SetMagField()) {
        throw kEH + "Setting of field failed!";
    }
    
    return kTRUE;
}

Bool_t AliEveEventManager::InitOCDB(int runNo)
{
  Bool_t ok = kTRUE;
  //first check/set the default OCDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if (!cdb->IsDefaultStorageSet())
  {
      TEnv settings;
      AliEveInit::GetConfig(&settings);
      
    TString ocdbStorage = settings.GetValue("OCDB.default.path",Form("local://%s/../src/OCDB",gSystem->Getenv("ALICE_ROOT")));
    if (gSystem->Getenv("ocdbStorage"))
    {
      printf("taking OCDB storage path from env ($ocdbStorage)\n");
      ocdbStorage = gSystem->Getenv("ocdbStorage");
    }
    else
    {
        //now if we don't have a GRP we need to get one from somewhere
        cout<<"AliEveEventManager::InitOCDB : GRP Data from prompt reco params"<<endl;
        fCurrentDataSource->ReceivePromptRecoParameters(runNo);
    }
    
    //on run change destroy the mag field, it will be reinitialized via AssertMagField/InitGRP
    if (runNo != cdb->GetRun())
    {
        delete TGeoGlobalMagField::Instance();
        new TGeoGlobalMagField();
        fgMaster->fgMagField=NULL;
    }
  }
  cdb->SetRun(runNo);
  //check is there is a GRP object for this run
  try {
    AliCDBManager::Instance()->Get("GRP/GRP/Data");
  } catch(...) {
    //now if we don't have a GRP we need to get one from somewhere
    fgMaster->SetCurrentRun(AliCDBManager::Instance()->GetRun());
    ok = fCurrentDataSource->ReceivePromptRecoParameters(runNo);
  }

  //on run change destroy the mag field, it will be reinitialized via AssertMagField/InitGRP
  if (runNo != cdb->GetRun())
  {
    delete TGeoGlobalMagField::Instance();
    new TGeoGlobalMagField();
    fgMaster->fgMagField=NULL;
  }

  return ok;
}

void AliEveEventManager::SetCdbUri(TString path) 
{
  AliCDBManager::Instance()->SetDefaultStorage(path);
}
