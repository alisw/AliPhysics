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
#include "AliEveConfigManager.h"
#include "AliEveVSDCreator.h"
#include "AliEveMultiView.h"

#include <THashList.h>
#include <TEveElement.h>
#include <TEveManager.h>
#include <TEveViewer.h>

#include <AliLog.h>
#include <AliRunLoader.h>
#include <AliRun.h>
#include <AliESDRun.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliAODEvent.h>

#include <AliRecoParam.h>
#include <AliCentralTrigger.h>
#include <AliCDBEntry.h>
#include <AliTriggerClass.h>
#include <AliTriggerConfiguration.h>
#include <AliTriggerCluster.h>
#include <AliDetectorRecoParam.h>

#include <AliDAQ.h>
#include <AliRawEventHeaderBase.h>
#include <AliRawReaderRoot.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliMagF.h>
#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliGRPObject.h>
#include <AliHeader.h>
#include <AliGeomManager.h>
#include <AliGRPManager.h>
#include <AliSysInfo.h>

#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TPRegexp.h>
#include <TError.h>
#include <TEnv.h>
#include <TString.h>
#include <TMap.h>
#include <TROOT.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TEveText.h>
#include <TEveTrans.h>

#ifdef ZMQ
#include "AliEveDataSourceOnline.h"

//#include "AliZMQManager.h"
//#include "AliOnlineReconstructionUtil.h"
//#include "AliGRPPreprocessor.h"
#endif

using std::cout;
using std::endl;
using std::vector;
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

//TString  AliEveEventManager::fgCdbUri;
//TString  AliEveEventManager::fgSpecificCdbUriValue;
//TString  AliEveEventManager::fgSpecificCdbUriPath;

Bool_t   AliEveEventManager::fgGRPLoaded    = kFALSE;
AliMagF* AliEveEventManager::fgMagField     = 0;
AliRecoParam* AliEveEventManager::fgRecoParam = 0;
Bool_t   AliEveEventManager::fgUniformField = kFALSE;

AliEveEventManager* AliEveEventManager::fgMaster  = NULL;

AliEveEventManager::AliEveEventManager(EDataSource defaultDataSource) :
TEveEventManager("Event", ""),
fEventId(-1),
fEventInfo(),
fAutoLoad  (kFALSE), fAutoLoadTime (5),fAutoLoadTimer(0),
fHasEvent(kFALSE),
fGlobal    (0), fGlobalReplace (kTRUE), fGlobalUpdate (kTRUE),
fExecutor    (0), fTransients(0), fTransientLists(0),
fPEventSelector(0),
fAutoLoadTimerRunning(kFALSE),
fViewsSaver(0),
fESDdrawer(0),
fSaveViews(false),
fDrawESDtracksByCategory(false),
fDrawESDtracksByType(false),
fFirstEvent(true),
fCenterProjectionsAtPrimaryVertex(false)
{
    InitInternals();
    
    fDataSourceOnline = new AliEveDataSourceOnline();
    fDataSourceOffline = new AliEveDataSourceOffline();
    
    ChangeDataSource(defaultDataSource);
}

AliEveEventManager::~AliEveEventManager()
{
    // Destructor.
    fAutoLoadTimer->Stop();
    fAutoLoadTimer->Disconnect("Timeout");
    fAutoLoadTimer->Disconnect("AutoLoadNextEvent");
    
    //    fTransients->DecDenyDestroy();
    //    fTransients->Destroy();
    
    //    fTransientLists->DecDenyDestroy();
    //    fTransientLists->Destroy();
    
    //delete fExecutor;
}

void AliEveEventManager::SetMaster(AliEveEventManager *master)
{
    if(fgMaster!=NULL)
    {
        AliFatal("Instance of AliEveEventManager already exists. Cannot create another one!!");
    }
    fgMaster = master;
}

void AliEveEventManager::InitInternals()
{
    // Initialize internal members.
    static const TEveException kEH("AliEveEventManager::InitInternals ");
    
    fgMaster = this;
    
    fAutoLoadTimer = new TTimer;
    fAutoLoadTimer->Connect("Timeout()", "AliEveEventManager", this, "AutoLoadNextEvent()");
    fAutoLoadTimer->Connect("Timeout()", "AliEveEventManager", this, "Timeout()");
    
    fExecutor = new AliEveMacroExecutor;
    
    fTransients = new TEveElementList("Transients", "Transient per-event elements.");
    fTransients->IncDenyDestroy();
    gEve->AddToListTree(fTransients, kFALSE);
    
    fTransientLists = new TEveElementList("Transient Lists", "Containers of transient elements.");
    fTransientLists->IncDenyDestroy();
    gEve->AddToListTree(fTransientLists, kFALSE);
    
    fPEventSelector = new AliEveEventSelector(this);
    
    fGlobal = new TMap; fGlobal->SetOwnerKeyValue();
    
    fViewsSaver = new AliEveSaveViews();
    fESDdrawer = new AliEveESDTracks();
}

void AliEveEventManager::ChangeDataSource(EDataSource newSource)
{
    if(newSource == kSourceOnline)
    {
        fCurrentDataSource = fDataSourceOnline;
    }
    else if(newSource == kSourceOffline)
    {
        fCurrentDataSource = fDataSourceOffline;
    }
    fCurrentData = fCurrentDataSource->GetData();
}

void AliEveEventManager::DestroyTransients()
{
    TEveManager::TRedrawDisabler rd(gEve);
    gEve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.
    
    gEve->GetViewers()->DeleteAnnotations();
    fTransients->DestroyElements();
    for (TEveElement::List_i i = fTransientLists->BeginChildren();
         i != fTransientLists->EndChildren(); ++i)
    {
        (*i)->DestroyElements();
    }
    DestroyElements();
    ElementChanged();
}

void AliEveEventManager::Timeout()
{
    Emit("Timeout()");
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

Bool_t AliEveEventManager::HasRunLoader()
{
    // Check if AliRunLoader is initialized.
    return fgMaster && fgMaster->fHasEvent && fgMaster->fCurrentData->fRunLoader;
}

Bool_t AliEveEventManager::HasESD()
{
    // Check if AliESDEvent is initialized.
    return fgMaster && fgMaster->fHasEvent && fgMaster->fCurrentData->fESD;
}

Bool_t AliEveEventManager::HasESDfriend()
{
    // Check if AliESDfriend is initialized.
    return fgMaster && fgMaster->fHasEvent && fgMaster->fCurrentData->fESDfriend;
}

Bool_t AliEveEventManager::HasAOD()
{
    // Check if AliESDEvent is initialized.
    return fgMaster && fgMaster->fHasEvent && fgMaster->fCurrentData->fAOD;
}

Bool_t AliEveEventManager::HasRawReader()
{
    // Check if raw-reader is initialized.
    return fgMaster && fgMaster->fHasEvent && fgMaster->fCurrentData->fRawReader;
}

AliRunLoader* AliEveEventManager::AssertRunLoader()
{
    // Make sure AliRunLoader is initialized and return it.
    // Throws exception in case run-loader is not available.
    // Static utility for macros.
    
    static const TEveException kEH("AliEveEventManager::AssertRunLoader ");
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE)
        throw (kEH + "ALICE event not ready.");
    if (fgMaster->fCurrentData->fRunLoader == 0)
        throw (kEH + "AliRunLoader not initialised.");
    return fgMaster->fCurrentData->fRunLoader;
}

AliESDEvent* AliEveEventManager::AssertESD()
{
    // Make sure AliESDEvent is initialized and return it.
    // Throws exception in case ESD is not available.
    // Static utility for macros.
    
    static const TEveException kEH("AliEveEventManager::AssertESD ");
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE)
        throw (kEH + "ALICE event not ready.");
    if (fgMaster->fCurrentData->fESD == 0)
        throw (kEH + "AliESD not initialised.");
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
    // Make sure raw-reader is initialized and return it.
    
    static const TEveException kEH("AliEveEventManager::AssertRawReader ");
    
    if (fgMaster == 0 || fgMaster->fHasEvent == kFALSE)
        throw (kEH + "ALICE event not ready.");
    if (fgMaster->fCurrentData->fRawReader == 0)
        throw (kEH + "RawReader not ready.");
    
    return fgMaster->fCurrentData->fRawReader;
}

//==============================================================================

AliMagF* AliEveEventManager::AssertMagField()
{
    // Make sure AliMagF is initialized and returns it.
    // Throws exception in case magnetic field is not available.
    // Static utility for macros.
    
    static const TEveException kEH("AliEveEventManager::AssertMagField ");
    
    if (fgMagField)
    {
        return fgMagField;
    }

    AliEveEventManager::GetMaster()->AssertESD()->InitMagneticField();
    if (TGeoGlobalMagField::Instance()->GetField())
    {
        fgMagField = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
        if (fgMagField == 0)
            throw kEH + "Global field set, but it is not AliMagF.";
        return fgMagField;
    }
    
    if (!fgGRPLoaded)
    {
        InitGRP();
    }
    
    if (TGeoGlobalMagField::Instance()->GetField())
    {
        fgMagField = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
        if (fgMagField == 0)
            throw kEH + "Global field set, but it is not AliMagF.";
    }
    else
    {
        throw kEH + "Could not initialize magnetic field.";
    }
    
    return fgMagField;
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

AliRecoParam* AliEveEventManager::AssertRecoParams()
{
    if(!fgRecoParam)
        InitRecoParam();
    
    return fgRecoParam;
}

Bool_t AliEveEventManager::InitRecoParam()
{
    // This is mostly a reap-off from reconstruction
    // The method accesses OCDB and retrieves all
    // the available reco-param objects from there.
    
    fgRecoParam = new AliRecoParam;
    const Int_t  kNDetectors = 14;
    
    static const TEveException kEH("AliEveEventManager::InitRecoParam");
    
    Bool_t isOK = kTRUE;
    
    if (fgRecoParam->GetDetRecoParamArray(kNDetectors)) {
        ::Info(kEH, "Using custom GRP reconstruction parameters");
    }
    else {
        ::Info(kEH, "Loading GRP reconstruction parameter objects");
        
        AliCDBPath path("GRP","Calib","RecoParam");
        AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
        if(!entry){
            ::Warning(kEH, "Couldn't find GRP RecoParam entry in OCDB");
            isOK = kFALSE;
        }
        else {
            TObject *recoParamObj = entry->GetObject();
            if (dynamic_cast<TObjArray*>(recoParamObj)) {
                // GRP has a normal TobjArray of AliDetectorRecoParam objects
                // Registering them in AliRecoParam
                fgRecoParam->AddDetRecoParamArray(kNDetectors,dynamic_cast<TObjArray*>(recoParamObj));
            }
            else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
                // GRP has only onse set of reco parameters
                // Registering it in AliRecoParam
                ::Info(kEH, "Single set of GRP reconstruction parameters found");
                dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
                fgRecoParam->AddDetRecoParam(kNDetectors,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
            }
            else {
                ::Error(kEH, "No valid GRP RecoParam object found in the OCDB");
                isOK = kFALSE;
            }
            entry->SetOwner(0);
        }
    }
    
    const char* fgkDetectorName[kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE" };
    
    
    for (Int_t iDet = 0; iDet < kNDetectors; iDet++) {
        
        if (fgRecoParam->GetDetRecoParamArray(iDet)) {
            ::Info(kEH, "Using custom reconstruction parameters for detector %s",fgkDetectorName[iDet]);
            continue;
        }
        
        ::Info(kEH, "Loading reconstruction parameter objects for detector %s",fgkDetectorName[iDet]);
        
        AliCDBPath path(fgkDetectorName[iDet],"Calib","RecoParam");
        AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
        if(!entry){
            ::Warning(kEH, "Couldn't find RecoParam entry in OCDB for detector %s",fgkDetectorName[iDet]);
            isOK = kFALSE;
        }
        else {
            TObject *recoParamObj = entry->GetObject();
            if (dynamic_cast<TObjArray*>(recoParamObj)) {
                // The detector has a normal TobjArray of AliDetectorRecoParam objects
                // Registering them in AliRecoParam
                fgRecoParam->AddDetRecoParamArray(iDet,dynamic_cast<TObjArray*>(recoParamObj));
            }
            else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
                // The detector has only onse set of reco parameters
                // Registering it in AliRecoParam
                ::Info(kEH, "Single set of reconstruction parameters found for detector %s",fgkDetectorName[iDet]);
                dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
                fgRecoParam->AddDetRecoParam(iDet,dynamic_cast<AliDetectorRecoParam*>(recoParamObj));
            }
            else {
                ::Error(kEH, "No valid RecoParam object found in the OCDB for detector %s",fgkDetectorName[iDet]);
                isOK = kFALSE;
            }
            entry->SetOwner(0);
            
        }
    }
    
    if(!isOK) {
        delete fgRecoParam;
        fgRecoParam = 0;
    }
    
    return isOK;
}

AliEveEventManager* AliEveEventManager::GetMaster()
{
    // Get master event-manager.
    if(fgMaster)
    {
        return fgMaster;
    }
    else
    {
        cout<<"FATAL -- Event Manager was not created. You must create it first with new AliEveEventManager()\n"<<endl;
        exit(0);
    }
}

void AliEveEventManager::RegisterTransient(TEveElement* element)
{
    GetMaster()->fTransients->AddElement(element);
}

void AliEveEventManager::RegisterTransientList(TEveElement* element)
{
    GetMaster()->fTransientLists->AddElement(element);
}

//------------------------------------------------------------------------------
// Autoloading of events
//------------------------------------------------------------------------------

void AliEveEventManager::SetAutoLoadTime(Float_t time)
{
    // Set the auto-load time in seconds
    fAutoLoadTime = time;
}

void AliEveEventManager::SetAutoLoad(Bool_t autoLoad)
{
    // Set the automatic event loading mode
    static const TEveException kEH("AliEveEventManager::SetAutoLoad ");
    
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
    fCurrentDataSource->NextEvent();
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
    
    static const TEveException kEH("AliEveEventManager::AfterNewEventLoaded ");
    
    Info(kEH, "------------------!!!------------");

    ElementChanged();
    
    NewEventDataLoaded();
    if (fExecutor) fExecutor->ExecMacros();
    
    TEveEventManager::AfterNewEventLoaded();
    NewEventLoaded();
    
    // tests of embedded text
    /*
    TTimeStamp ts(fCurrentData->fESD->GetTimeStamp());
    
    TEveText *txt = new TEveText(Form("Run:%d",fCurrentData->fESD->GetRunNumber()));
    TEveText *txt2 = new TEveText(Form("Timestamp:%s",ts.AsString("s")));
    txt->SetMainColor(kWhite);
    txt->SetFontSize(20);
    txt->SetMainTransparency('1');
    
    txt2->SetMainColor(kWhite);
    txt2->SetFontSize(20);
    txt2->SetMainTransparency('1');
    
    TEveTrans trans = txt->RefMainTrans();
    double *arr = trans.Array();
    arr[12]=0;
    arr[13]=-1200;
    arr[14]=0;
    txt->SetTransMatrix(arr);

    arr[13]+=100;
    txt2->SetTransMatrix(arr);
    
    gEve->AddElement(txt);
    gEve->AddElement(txt2);
     */
    //
    
    if(HasESD())
    {
        if(fDrawESDtracksByCategory)fESDdrawer->ByCategory();
        if(fDrawESDtracksByType)fESDdrawer->ByType();
        
        Double_t x[3] = { 0, 0, 0 };
        
        fCurrentData->fESD->GetPrimaryVertex()->GetXYZ(x);
        
        TTimeStamp ts(fCurrentData->fESD->GetTimeStamp());
        TString win_title("Eve Main Window -- Timestamp: ");
        win_title += ts.AsString("s");
        win_title += "; Event # in ESD file: ";
        win_title += fCurrentData->fESD->GetEventNumberInFile();
        gEve->GetBrowser()->SetWindowName(win_title);
        
        TEveElement* top = gEve->GetCurrentEvent();
        
        AliEveMultiView *mv = AliEveMultiView::Instance();
        
//        mv->DestroyEventRPhi();
        if (fCenterProjectionsAtPrimaryVertex){
            mv->SetCenterRPhi(x[0], x[1], x[2]);
        }
        mv->ImportEventRPhi(top);
        
//        mv->DestroyEventRhoZ();
        if (fCenterProjectionsAtPrimaryVertex){
            mv->SetCenterRhoZ(x[0], x[1], x[2]);
        }
        mv->ImportEventRhoZ(top);
        
        if (fCenterProjectionsAtPrimaryVertex)
            mv->SetCenterMuon(x[0], x[1], x[2]);
        mv->ImportEventMuon(top);
        
 
        gEve->GetBrowser()->RaiseWindow();
        gEve->FullRedraw3D();
        gSystem->ProcessEvents();
        
        if(fFirstEvent)
        {
            gROOT->ProcessLine(".x geom_emcal.C");
            fFirstEvent=false;
        }
    }
    
//    if(fSaveViews  && fCurrentData->fESD->GetNumberOfTracks()>0)
//    {
//        fViewsSaver->Save();
//        fViewsSaver->SendToAmore();
//    }
    
//    if (this == fgMaster && fSubManagers != 0)
//    {
//        TIter next(fSubManagers);
//        while ((fgCurrent = dynamic_cast<AliEveEventManager*>(next())) != 0)
//        {
//            gEve->SetCurrentEvent(fgCurrent);
//            try
//            {
//                fgCurrent->GotoEvent(fEventId);
//            }
//            catch (TEveException& exc)
//            {
//                // !!! Should somehow tag / disable / remove it?
//                Error(kEH, "Getting event %d for sub-event-manager '%s' failed: '%s'.",
//                      fEventId, fgCurrent->GetName(), exc.Data());
//            }
//            Info(kEH, "------------------!!! while() gEve->SetCurrentEvent() ------------");
//        }
//        fgCurrent = fgMaster;
//        Info(kEH, "------------------!!! while() gEve->SetCurrentEvent(MASTER) ------------");
//        gEve->SetCurrentEvent(fgMaster);
//    }
}

void AliEveEventManager::NewEventDataLoaded()
{
    // Emit NewEventDataLoaded signal.
    Emit("NewEventDataLoaded()");
}
void AliEveEventManager::NewEventLoaded()
{
    // Emit NewEventLoaded signal.
    Emit("NewEventLoaded()");
}
void AliEveEventManager::NoEventLoaded()
{
    // Emit NoEventLoaded signal.
    Emit("NoEventLoaded()");
}

//------------------------------------------------------------------------------
// Event info dumpers
//------------------------------------------------------------------------------

const AliEventInfo* AliEveEventManager::GetEventInfo()
{
    // Fill the event info object
    
    AliCentralTrigger *aCTP = NULL;
    if (fCurrentData->fRawReader) {
        fEventInfo.SetEventType(fCurrentData->fRawReader->GetType());
        
        ULong64_t mask = fCurrentData->fRawReader->GetClassMask();
        fEventInfo.SetTriggerMask(mask);
        UInt_t clmask = fCurrentData->fRawReader->GetDetectorPattern()[0];
        fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(clmask));
        
        aCTP = new AliCentralTrigger();
        TString configstr("");
        if (!aCTP->LoadConfiguration(configstr)) { // Load CTP config from OCDB
            AliError("No trigger configuration found in OCDB! The trigger configuration information will not be used!");
            delete aCTP;
            return 0;
        }
        aCTP->SetClassMask(mask);
        aCTP->SetClusterMask(clmask);
        
        if (fCurrentData->fRunLoader) {
            AliCentralTrigger* rlCTP = fCurrentData->fRunLoader->GetTrigger();
            if (rlCTP) {
                rlCTP->SetClassMask(mask);
                rlCTP->SetClusterMask(clmask);
            }
        }
    }
    else {
        fEventInfo.SetEventType(AliRawEventHeaderBase::kPhysicsEvent);
        
        if (fCurrentData->fRunLoader && (!fCurrentData->fRunLoader->LoadTrigger())) {
            aCTP = fCurrentData->fRunLoader->GetTrigger();
            fEventInfo.SetTriggerMask(aCTP->GetClassMask());
            // get inputs from actp - just get
            AliESDHeader* esdheader = fCurrentData->fESD->GetHeader();
            esdheader->SetL0TriggerInputs(aCTP->GetL0TriggerInputs());
            esdheader->SetL1TriggerInputs(aCTP->GetL1TriggerInputs());
            esdheader->SetL2TriggerInputs(aCTP->GetL2TriggerInputs());
            fEventInfo.SetTriggerCluster(AliDAQ::ListOfTriggeredDetectors(aCTP->GetClusterMask()));
        }
        else {
            AliWarning("No trigger can be loaded! The trigger information will not be used!");
            return 0;
        }
    }
    
    AliTriggerConfiguration *config = aCTP->GetConfiguration();
    if (!config) {
        AliError("No trigger configuration has been found! The trigger configuration information will not be used!");
        if (fCurrentData->fRawReader) delete aCTP;
        return 0;
    }
    
    TString declTriggerClasses;
    
    // Load trigger aliases and declare the trigger classes included in aliases
    AliCDBEntry * entry = AliCDBManager::Instance()->Get("GRP/CTP/Aliases");
    if (entry) {
        THashList * lst = dynamic_cast<THashList*>(entry->GetObject());
        if (lst) {
            lst->Sort(kSortDescending); // to avoid problems with substrings
            if (fCurrentData->fRawReader) fCurrentData->fRawReader->LoadTriggerAlias(lst);
            // Now declare all the triggers present in the aliases
            TIter iter(lst);
            TNamed *nmd = 0;
            while((nmd = dynamic_cast<TNamed*>(iter.Next()))){
                declTriggerClasses += " ";
                declTriggerClasses += nmd->GetName();
            }
        }
        else {
            AliError("Cannot cast the object with trigger aliases to THashList!");
        }
    }
    else {
        AliError("No OCDB entry for the trigger aliases!");
    }
    
    // Load trigger classes for this run
    UChar_t clustmask = 0;
    TString trclasses;
    ULong64_t trmask = fEventInfo.GetTriggerMask();
    const TObjArray& classesArray = config->GetClasses();
    Int_t nclasses = classesArray.GetEntriesFast();
    for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
        AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
        if (trclass && trclass->GetMask()>0) {
            Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
            if (fCurrentData->fESD) fCurrentData->fESD->SetTriggerClass(trclass->GetName(),trindex);
            if (fCurrentData->fRawReader) fCurrentData->fRawReader->LoadTriggerClass(trclass->GetName(),trindex);
            if (trmask & (1ull << trindex)) {
                trclasses += " ";
                trclasses += trclass->GetName();
                trclasses += " ";
                clustmask |= trclass->GetCluster()->GetClusterMask();
            }
        }
    }
    fEventInfo.SetTriggerClasses(trclasses);
    
    if (!aCTP->CheckTriggeredDetectors()) {
        if (fCurrentData->fRawReader) delete aCTP;
        return 0;
    }
    
    if (fCurrentData->fRawReader) delete aCTP;
    
    // everything went ok, return pointer
    return (&fEventInfo);
}


TString AliEveEventManager::GetEventInfoHorizontal() const
{
    // Dumps the event-header contents in vertical formatting.
    
    TString rawInfo, esdInfo;
    
    if (!fCurrentData->fRawReader)
    {
        rawInfo = "No raw-data event info is available!\n";
    }
    else
    {
        const UInt_t* attr = fCurrentData->fRawReader->GetAttributes();
        TTimeStamp ts(fCurrentData->fRawReader->GetTimestamp());
        rawInfo.Form("RAW event info: Run#: %d  Event type: %d (%s)  Period: %x  Orbit: %x  BC: %x\n"
                     "Trigger: %llx\nDetectors: %x (%s)\nAttributes:%x-%x-%x  Timestamp: %s\n",
                     fCurrentData->fRawReader->GetRunNumber(),fCurrentData->fRawReader->GetType(),AliRawEventHeaderBase::GetTypeName(fCurrentData->fRawReader->GetType()),
                     fCurrentData->fRawReader->GetPeriod(),fCurrentData->fRawReader->GetOrbitID(),fCurrentData->fRawReader->GetBCID(),
                     fCurrentData->fRawReader->GetClassMask(),
                     *(fCurrentData->fRawReader)->GetDetectorPattern(),AliDAQ::ListOfTriggeredDetectors(*(fCurrentData->fRawReader)->GetDetectorPattern()),
                     attr[0],attr[1],attr[2], ts.AsString("s"));
    }
    
    if (!fCurrentData->fESD)
    {
        esdInfo = "No ESD event info is available!";
    }
    else
    {
        TString acttrclasses   = fCurrentData->fESD->GetESDRun()->GetActiveTriggerClasses();
        TString firedtrclasses = fCurrentData->fESD->GetFiredTriggerClasses();
        TTimeStamp ts(fCurrentData->fESD->GetTimeStamp());
        esdInfo.Form("ESD event info: Run#: %d  Event type: %d (%s)  Period: %x  Orbit: %x  BC: %x\n"
                     "Active trigger classes: %s\nTrigger: %llx (%s)\nEvent# in file: %d  Timestamp: %s, MagField: %.2e",
                     fCurrentData->fESD->GetRunNumber(),
                     fCurrentData->fESD->GetEventType(),AliRawEventHeaderBase::GetTypeName(fCurrentData->fESD->GetEventType()),
                     fCurrentData->fESD->GetPeriodNumber(),fCurrentData->fESD->GetOrbitNumber(),fCurrentData->fESD->GetBunchCrossNumber(),
                     acttrclasses.Data(),
                     fCurrentData->fESD->GetTriggerMask(),firedtrclasses.Data(),
                     fCurrentData->fESD->GetEventNumberInFile(), ts.AsString("s"), fCurrentData->fESD->GetMagneticField());
    }
    
    return rawInfo + esdInfo;
}

TString AliEveEventManager::GetEventInfoVertical() const
{
    // Dumps the event-header contents in vertical formatting.
    
    TString rawInfo, esdInfo;
    
    if (!fCurrentData->fRawReader)
    {
        rawInfo = "No raw-data event info is available!\n";
    }
    else
    {
        const UInt_t* attr = fCurrentData->fRawReader->GetAttributes();
        rawInfo.Form("Raw-data event info:\nRun#: %d\nEvent type: %d (%s)\nPeriod: %x\nOrbit: %x   BC: %x\nTrigger: %llx\nDetectors: %x (%s)\nAttributes:%x-%x-%x\nTimestamp: %x\n",
                     fCurrentData->fRawReader->GetRunNumber(),fCurrentData->fRawReader->GetType(),AliRawEventHeaderBase::GetTypeName(fCurrentData->fRawReader->GetType()),
                     fCurrentData->fRawReader->GetPeriod(),fCurrentData->fRawReader->GetOrbitID(),fCurrentData->fRawReader->GetBCID(),
                     fCurrentData->fRawReader->GetClassMask(),
                     *(fCurrentData->fRawReader)->GetDetectorPattern(),AliDAQ::ListOfTriggeredDetectors(*(fCurrentData->fRawReader)->GetDetectorPattern()),
                     attr[0],attr[1],attr[2],
                     fCurrentData->fRawReader->GetTimestamp());
    }
    
    if (!fCurrentData->fESD)
    {
        esdInfo = "No ESD event info is available!\n";
    }
    else
    {
        TString acttrclasses   = fCurrentData->fESD->GetESDRun()->GetActiveTriggerClasses();
        TString firedtrclasses = fCurrentData->fESD->GetFiredTriggerClasses();
        esdInfo.Form("ESD event info:\nRun#: %d\nActive trigger classes: %s\nEvent type: %d (%s)\nPeriod: %x\nOrbit: %x   BC: %x\nTrigger: %llx (%s)\nEvent# in file:%d\nTimestamp: %x\n",
                     fCurrentData->fESD->GetRunNumber(),
                     acttrclasses.Data(),
                     fCurrentData->fESD->GetEventType(),AliRawEventHeaderBase::GetTypeName(fCurrentData->fESD->GetEventType()),
                     fCurrentData->fESD->GetPeriodNumber(),fCurrentData->fESD->GetOrbitNumber(),fCurrentData->fESD->GetBunchCrossNumber(),
                     fCurrentData->fESD->GetTriggerMask(),firedtrclasses.Data(),
                     fCurrentData->fESD->GetEventNumberInFile(),
                     fCurrentData->fESD->GetTimeStamp());
    }
    
    return rawInfo + "\n" + esdInfo;
}


//==============================================================================
// Reading of GRP and MagneticField.
// This is a reap-off from reconstruction ... should really be a common
// code to do this somewhere in STEER.
//==============================================================================


Bool_t AliEveEventManager::InitGRP()
{
    //------------------------------------
    // Initialization of the GRP entry
    //------------------------------------
    
    static const TEveException kEH("AliEveEventManager::InitGRP ");
    
    AliGRPManager grpMgr;
    if (!grpMgr.ReadGRPEntry()) {
        return kFALSE;
    }
    fgGRPLoaded = kTRUE;
    if (!grpMgr.SetMagField()) {
        throw kEH + "Setting of field failed!";
    }
    
    //*** Get the diamond profiles from OCDB
    // Eventually useful.
    
 
//     entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexSPD");
//     if (entry) {
//     fDiamondProfileSPD = dynamic_cast<AliESDVertex*> (entry->GetObject());
//     } else {
//     ::Error(kEH, "No SPD diamond profile found in OCDB!");
//     }
//     
//     entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertex");
//     if (entry) {
//     fDiamondProfile = dynamic_cast<AliESDVertex*> (entry->GetObject());
//     } else {
//     ::Error(kEH, "No diamond profile found in OCDB!");
//     }
//     
//     entry = AliCDBManager::Instance()->Get("GRP/Calib/MeanVertexTPC");
//     if (entry) {
//     fDiamondProfileTPC = dynamic_cast<AliESDVertex*> (entry->GetObject());
//     } else {
//     ::Error(kEH, "No TPC diamond profile found in OCDB!");
//     }
 
    
    return kTRUE;
}

//------------------------------------
// Global variables management
//------------------------------------

Bool_t AliEveEventManager::InsertGlobal(const TString& tag, TEveElement* model)
{
    // Insert a new visualization-parameter database entry with the default
    return InsertGlobal(tag, model, fGlobalReplace, fGlobalUpdate);
}

Bool_t AliEveEventManager::InsertGlobal(const TString& tag, TEveElement* model,
                                        Bool_t replace, Bool_t update)
{
    TPair* pair = (TPair*) fGlobal->FindObject(tag);
    if (pair)
    {
        if (replace)
        {
            model->IncDenyDestroy();
            model->SetRnrChildren(kFALSE);
            
            TEveElement* old_model = dynamic_cast<TEveElement*>(pair->Value());
            if(!old_model) AliFatal("old_model == 0, dynamic cast failed\n");
            while (old_model->HasChildren())
            {
                TEveElement *el = old_model->FirstChild();
                el->SetVizModel(model);
                if (update)
                {
                    el->CopyVizParams(model);
                    el->PropagateVizParamsToProjecteds();
                }
            }
            old_model->DecDenyDestroy();
            
            pair->SetValue(dynamic_cast<TObject*>(model));
            return kTRUE;
        }
        else
        {
            return kFALSE;
        }
    }
    else
    {
        model->IncDenyDestroy();
        model->SetRnrChildren(kFALSE);
        fGlobal->Add(new TObjString(tag), dynamic_cast<TObject*>(model));
        return kTRUE;
    }
}

TEveElement* AliEveEventManager::FindGlobal(const TString& tag)
{
    return dynamic_cast<TEveElement*>(fGlobal->GetValue(tag));
}



