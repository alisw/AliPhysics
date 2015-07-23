//
//  AliEveDataSourceOffline.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#include "AliEveDataSourceOffline.h"

#include "AliSysInfo.h"
#include "AliEveEventSelector.h"
#include "AliHeader.h"

#include "TSystem.h"
#include "TEnv.h"
#include <TPRegexp.h>

#include <iostream>

using namespace std;

ClassImp(AliEveDataSourceOffline);

AliEveDataSourceOffline::AliEveDataSourceOffline(bool storageManager) :
AliEveDataSource(),
fgESDFileName("AliESDs.root"),
fgESDfriendsFileName("AliESDfriends.root"),
fgAODFileName("AliAOD.root"),
fgGAliceFileName("galice.root"),
fgRawFileName("raw.root"),
fEventManager(0),
fgESDvisibleTrees(kOfflineTree),
fgAssertRunLoader(false),
fgAssertESD(false),
fgAssertAOD(false),
fgAssertRaw(false),
fIsOpen(false),
fESDfriendExists(kFALSE),
fEventInfo(),
fgAODfriends(0),
fgRawFromStandardLoc(false)
{
    cout<<"Constructor of AliEveDataSourceOffline"<<endl;
    cout<<"AliEveData initialized"<<endl;
    
    fEventManager = AliEveEventManager::GetMaster();
    Open();
    cout<<"Files opened"<<endl;
}

AliEveDataSourceOffline::~AliEveDataSourceOffline()
{
    
}


void AliEveDataSourceOffline::SetEvent(AliRunLoader *runLoader, AliRawReader *rawReader, AliESDEvent *esd, AliESDfriend *esdf)
{
    // Set an event from an external source.
    // The method is used in the online visualisation.
    // AOD is not supported.
    
    static const TEveException kEH("AliEveEventManager::SetEvent ");
    
    if (fIsOpen)
    {
        Warning(kEH, "Event-files were open. Closing and switching to external control.");
        Close();
    }
    
    Info(kEH,"setting it!!! ============================");
    
    fCurrentData.fRunLoader = runLoader;
    fCurrentData.fRawReader = rawReader;
    fCurrentData.fESD       = esd;
    fCurrentData.fESDfriend = esdf;
    fCurrentData.fAOD       = 0;
    
    fEventManager->SetEventId(fEventManager->GetEventId()+1);
    fEventManager->SetHasEvent(true);
    
    SetTitle("Online event in memory");
    SetName ("Online Event");
    fEventManager->ElementChanged();
    fEventManager->AfterNewEventLoaded();
    
    if (fEventManager->GetAutoLoad()){
        fEventManager->StartAutoLoadTimer();
    }
    
}

void AliEveDataSourceOffline::GotoEvent(Int_t event)
{
    cout<<"Go to event:"<<event<<endl;
    // Load data for specified event.
    // If event is out of range an exception is thrown and old state
    // is preserved.
    // After successful loading of event, the virtual function
    // AfterNewEventLoaded() is called. This executes commands that
    // were registered via TEveEventManager::AddNewEventCommand().
    //
    // If event is negative, it is subtracted from the number of
    // available events, thus passing -1 will load the last event.
    // This is not supported when raw-data is the only data-source
    // as the number of events is not known.
    
    static const TEveException kEH("AliEveEventManager::GotoEvent ");
    
    if(!fCurrentData.fESD)
    {
        cout<<"No ESD event avaliable. Probably files were not opened."<<endl;
        return;
    }
    
    if(fCurrentData.fESD->GetRunNumber() != fEventManager->GetCurrentRun())
    {
        fEventManager->SetCurrentRun(fCurrentData.fESD->GetRunNumber());
        fEventManager->ResetMagneticField();
    }
    
    if (!fIsOpen){throw (kEH + "Event-files not opened but ED is in offline mode.");}
    
    fEventInfo.Reset();
    fEventManager->SetHasEvent(false);
    
    Int_t maxEvent = 0;
    if ((fCurrentData.fESDTree!=0) || (fCurrentData.fHLTESDTree!=0))
    {
        if(fCurrentData.fESDTree){
            if (event >= fCurrentData.fESDTree->GetEntries())
                fCurrentData.fESDTree->Refresh();
            maxEvent = fCurrentData.fESDTree->GetEntries() - 1;
            if (event < 0)
                event = fCurrentData.fESDTree->GetEntries() + event;
        }
        
        if(fCurrentData.fHLTESDTree){
            if (event >= fCurrentData.fHLTESDTree->GetEntries())
                fCurrentData.fHLTESDTree->Refresh();
            maxEvent = fCurrentData.fHLTESDTree->GetEntries() - 1;
            if (event < 0)
                event = fCurrentData.fHLTESDTree->GetEntries() + event;
            
        }
    }
    else if (fCurrentData.fAODTree)
    {
        maxEvent = fCurrentData.fAODTree->GetEntries() - 1;
        if (event < 0)
            event = fCurrentData.fAODTree->GetEntries() + event;
    }
    else if (fCurrentData.fRunLoader)
    {
        maxEvent = fCurrentData.fRunLoader->GetNumberOfEvents() - 1;
        if (event < 0)
            event = fCurrentData.fRunLoader->GetNumberOfEvents() + event;
    }
    else if (fCurrentData.fRawReader)
    {
        maxEvent = fCurrentData.fRawReader->GetNumberOfEvents() - 1;
        if (maxEvent < 0)
        {
            maxEvent = 10000000;
            if (event < 0) {
                Error(kEH, "current raw-data source does not support direct event access.");
                return;
            }
            Info(kEH, "number of events unknown for current raw-data source, setting max-event id to 10M.");
        }
        else
        {
            if (event < 0)
                event = fCurrentData.fRawReader->GetNumberOfEvents() + event;
        }
    }
    else
    {
        throw (kEH + "neither RunLoader, ESD nor Raw loaded.");
    }
    if (event < 0)
    {
        throw (kEH + Form("event %d not present, available range [%d, %d].",
                          event, 0, maxEvent));
    }
    if (event > maxEvent)
    {
        event=0;
        cout<<"Event number out of range. Going to event 0"<<endl;
    }

    TString sysInfoHeader;
    sysInfoHeader.Form("AliEveEventManager::GotoEvent(%d) - ", event);
    AliSysInfo::AddStamp(sysInfoHeader + "Start");

    fEventManager->DestroyTransients();
    
//    TEveManager::TRedrawDisabler rd(gEve);
//    gEve->Redraw3D(kFALSE, kTRUE); // Enforce drop of all logicals.
//    
//    // !!! MT this is somewhat brutal; at least optionally, one could be
//    // a bit gentler, checking for objs owning their external refs and having
//    // additinal parents.
//    gEve->GetViewers()->DeleteAnnotations();
//    fTransients->DestroyElements();
//    for (TEveElement::List_i i = fTransientLists->BeginChildren();
//         i != fTransientLists->EndChildren(); ++i)
//    {
//        (*i)->DestroyElements();
//    }
//    DestroyElements();
    
    AliSysInfo::AddStamp(sysInfoHeader + "PostDestroy");
    
    if (fCurrentData.fESDTree) {
        if (fCurrentData.fESDTree->GetEntry(event) <= 0)
            throw (kEH + "failed getting required event from ESD.");
        
        if (fESDfriendExists)
            fCurrentData.fESD->SetESDfriend(fCurrentData.fESDfriend);
    }
    if (fCurrentData.fHLTESDTree) {
        if (fCurrentData.fHLTESDTree->GetEntry(event) <= 0)
            throw (kEH + "failed getting required event from HLT ESD.");
        
        if (fESDfriendExists)
            fCurrentData.fESD->SetESDfriend(fCurrentData.fESDfriend);
    }
    
    if (fCurrentData.fAODTree) {
        if (fCurrentData.fAODTree->GetEntry(event) <= 0)
            throw (kEH + "failed getting required event from AOD.");
    }
    
    if (fCurrentData.fRunLoader) {
        if (fCurrentData.fRunLoader->GetEvent(event) != 0)
            throw (kEH + "failed getting required event.");
    }
    if (fCurrentData.fRawReader)
    {
        // AliRawReader::GotoEvent(Int_t) works for AliRawReaderRoot/Chain.
        if (fCurrentData.fRawReader->GotoEvent(event) == kFALSE)
        {
            // Use fallback method - iteration with NextEvent().
            Int_t rawEv = fEventManager->GetEventId();
            if (event < rawEv)
            {
                fCurrentData.fRawReader->RewindEvents();
                rawEv = -1;
            }
            
            while (rawEv < event)
            {
                if ( ! fCurrentData.fRawReader->NextEvent())
                {
                    fCurrentData.fRawReader->RewindEvents();
                    fEventManager->SetEventId(-1);
//                    fEventId = -1;
                    throw (kEH + Form("Error going to next raw-event from event %d.", rawEv));
                }
                ++rawEv;
            }
            Warning(kEH, "Loaded raw-event %d with fallback method.\n", rawEv);
        }
    }
    
    fEventManager->SetHasEvent(true);
    fEventManager->SetEventId(event);
    
    SetName(Form("Event %d", fEventManager->GetEventId()));
    fEventManager->ElementChanged();
    
    AliSysInfo::AddStamp(sysInfoHeader + "PostLoadEvent");
    fEventManager->AfterNewEventLoaded();
    AliSysInfo::AddStamp(sysInfoHeader + "PostUserActions");
}

void AliEveDataSourceOffline::NextEvent()
{
    // Loads next event.
    // Does magick needed for online display when under external event control.
    
    static const TEveException kEH("AliEveEventManager::NextEvent ");
    
//    if (fAutoLoadTimerRunning){throw (kEH + "Event auto-load timer is running.");}
    
    /*
     if ((fCurrentData.fESDTree!=0) || (fCurrentData.fHLTESDTree!=0))
     {
     cout<<"There is ESD or HLTESD tree"<<endl;
     Int_t nextevent=0;
     if (fPEventSelector->FindNext(nextevent))
     {
     cout<<"GotoEvent:"<<nextevent<<endl;
     GotoEvent(nextevent);
     }
     }*/
    /*else*/ if (fEventManager->GetEventId() < GetMaxEventId(kTRUE))
    {
        cout<<"GotoEvent:"<<fEventManager->GetEventId()+1<<endl;
        GotoEvent(fEventManager->GetEventId() + 1);
    }
    else
    {
        cout<<"Going back to event 0"<<endl;
        GotoEvent(0);
    }
    
    gSystem->ProcessEvents();
}

void AliEveDataSourceOffline::PrevEvent()
{
    // Loads previous event.
    static const TEveException kEH("AliEveEventManager::PrevEvent ");
    
    if ((fCurrentData.fESDTree!=0) || (fCurrentData.fHLTESDTree!=0))
    {
        Int_t nextevent=0;
        if (fEventManager->GetEventSelector()->FindPrev(nextevent))
        {
            GotoEvent(nextevent);
        }
    }
    else if (fEventManager->GetEventId() > 0)
    {
        GotoEvent(fEventManager->GetEventId() - 1);
    }
}

void AliEveDataSourceOffline::Open()
{
    // Open event-data from URL specified in path.
    // Attempts to create AliRunLoader() and to open ESD with ESDfriends.
    // Warning is reported if run-loader or ESD is not found.
    // Global data-members fgAssertRunLoader and fgAssertESD can be set
    // to throw exceptions instead.
    
    static const TEveException kEH("AliEveEventManager::Open ");
    if (fIsOpen){throw (kEH + "Event-files already opened.");}
    
    Int_t runNo = -1;

    // Open ESD and ESDfriends
    if ((fCurrentData.fESDFile = TFile::Open(fgESDFileName)))
    {
        fCurrentData.fESD = new AliESDEvent();
        
        switch(fgESDvisibleTrees){
            case kOfflineTree :
                fCurrentData.fESDTree = readESDTree("esdTree", runNo);
                break;
            case kHLTTree :
                fCurrentData.fHLTESDTree = readESDTree("HLTesdTree", runNo);
                break;
            default:
                fCurrentData.fESDTree    = readESDTree("esdTree", runNo);
                fCurrentData.fHLTESDTree = readESDTree("HLTesdTree", runNo);
        }
        
        if(!fCurrentData.fESDTree && !fCurrentData.fHLTESDTree){
            // both ESD trees are == 0
            delete fCurrentData.fESDFile; fCurrentData.fESDFile = 0;
            delete fCurrentData.fESD; fCurrentData.fESD = 0;
        }
        
        
    }
    else{Warning(kEH, "can not read ESD file '%s'.", fgESDFileName.Data());}
    if (fCurrentData.fESDTree == 0 && fCurrentData.fHLTESDTree==0)
    {
        if (fgAssertESD){throw (kEH + "ESD not initialized. Its precence was requested.");}
        else {Warning(kEH, "ESD not initialized.");}
    }
    
    // Open AOD and registered friends
    if ( (fCurrentData.fAODFile = TFile::Open(fgAODFileName)) )
    {
        fCurrentData.fAOD = new AliAODEvent();
        fCurrentData.fAODTree = (TTree*) fCurrentData.fAODFile->Get("aodTree");
        if (fCurrentData.fAODTree != 0)
        {
            // Check if AODfriends exist and attach them.
            TIter       friends(fgAODfriends);
            TObjString *name;
            while ((name = (TObjString*) friends()) != 0)
            {
                TString p(Form("%s/%s", fgAODFileName.Data(), name->GetName()));
                if (fgAODFileName.EndsWith(".zip")) p.Form("%s#%s",fgAODFileName.Data(),name->GetName());
                if (gSystem->AccessPathName(p, kReadPermission) == kFALSE)
                {
                    fCurrentData.fAODTree->AddFriend("aodTree", name->GetName());
                }
            }
            
            fCurrentData.fAOD->ReadFromTree(fCurrentData.fAODTree);
            
            if (fCurrentData.fAODTree->GetEntry(0) <= 0)
            {
                delete fCurrentData.fAODFile; fCurrentData.fAODFile = 0;
                delete fCurrentData.fAOD;     fCurrentData.fAOD     = 0;
                Warning(kEH, "failed getting the first entry from addTree.");
            }
            else if (runNo < 0){runNo = fCurrentData.fAOD->GetRunNumber();}
        }
        else // aodtree == 0
        {
            delete fCurrentData.fAODFile; fCurrentData.fAODFile = 0;
            delete fCurrentData.fAOD;     fCurrentData.fAOD     = 0;
            Warning(kEH, "failed getting the aodTree.");
        }
    }
    else // aod not readable
    {
        Warning(kEH, "can not read AOD file '%s'.", fgAODFileName.Data());
    }
    if (fCurrentData.fAODTree == 0)
    {
        if (fgAssertAOD){throw (kEH + "AOD not initialized. Its precence was requested.");}
        else {Warning(kEH, "AOD not initialized.");}
    }
    // Open RunLoader from galice.root
    //    fgGAliceFileName = "/Users/Jerus/galice.root"; // temp
    
    TFile *gafile = TFile::Open(fgGAliceFileName);
    cout<<"Opening galice"<<endl;
    if (gafile)
    {
        gafile->Close();
        delete gafile;
        cout<<"SETTING RUN LOADER in Open()"<<endl;
        fCurrentData.fRunLoader = AliRunLoader::Open(fgGAliceFileName, fEventManager->GetName());
        if (fCurrentData.fRunLoader)
        {
            TString alicePath(gSystem->DirName(fgGAliceFileName));
            alicePath.Append("/");
            fCurrentData.fRunLoader->SetDirName(alicePath);
            
            if (fCurrentData.fRunLoader->LoadgAlice() != 0){Warning(kEH, "failed loading gAlice via run-loader.");}
            
            if (fCurrentData.fRunLoader->LoadHeader() == 0){
                if(runNo < 0){
                    runNo = fCurrentData.fRunLoader->GetHeader()->GetRun();
                }
            }
            else{
                Warning(kEH, "failed loading run-loader's header.");
                delete fCurrentData.fRunLoader;
                fCurrentData.fRunLoader = 0;
            }
        }
        else{Warning(kEH, "failed opening ALICE run-loader from '%s'.", fgGAliceFileName.Data());}
    }
    else{Warning(kEH, "can not read '%s'.", fgGAliceFileName.Data());}
    
    if (fCurrentData.fRunLoader == 0)
    {
        if (fgAssertRunLoader){throw (kEH + "Bootstraping of run-loader failed. Its precence was requested.");}
        else{Warning(kEH, "Bootstraping of run-loader failed.");}
    }
    // Open raw-data file
    TString rawPath;
    if (fgRawFromStandardLoc)
    {
        if (!fgRawFileName.BeginsWith("alien:")){
            throw kEH + "Standard raw search requested, but the directory is not in AliEn.";
        }
        if (!fgRawFileName.Contains("/ESDs/")){
            throw kEH + "Standard raw search requested, but does not contain 'ESDs' directory.";
        }
        
        TPMERegexp chunk("/([\\d\\.])+/?$");
        Int_t nm = chunk.Match(fgRawFileName);
        if (nm != 2){
            throw kEH + "Standard raw search requested, but the path does not end with chunk-id directory.";
        }
        
        TPMERegexp esdstrip("/ESDs/.*");
        rawPath = fgRawFileName;
        esdstrip.Substitute(rawPath, "/raw/");
        rawPath += chunk[0];
        rawPath += ".root";
        
        Info(kEH, "Standard raw search requested, using the following path:\n  %s\n", rawPath.Data());
    }
    else
    {
        rawPath = fgRawFileName;
    }
    // If i use open directly, raw-reader reports an error but i have
    // no way to detect it.
    // Is this (AccessPathName check) ok for xrootd / alien? Yes, not for http.
    AliLog::EType_t oldLogLevel = (AliLog::EType_t) AliLog::GetGlobalLogLevel();
    if (fgAssertRaw == kFALSE){AliLog::SetGlobalLogLevel(AliLog::kFatal);}
    
    if (gSystem->AccessPathName(rawPath, kReadPermission) == kFALSE){
        fCurrentData.fRawReader = AliRawReader::Create(rawPath);
    }
    else{
        fCurrentData.fRawReader = AliRawReader::Create(fgRawFileName);
    }
    
    if (fgAssertRaw == kFALSE){AliLog::SetGlobalLogLevel(oldLogLevel);}
    
    if (fCurrentData.fRawReader == 0)
    {
        if (fgAssertRaw){throw (kEH + "raw-data not initialized. Its precence was requested.");}
        else{Warning(kEH, "raw-data not initialized.");}
    }
    if (runNo < 0)
    {
        if (fCurrentData.fRawReader)
        {
            if (!fCurrentData.fRawReader->NextEvent()){throw (kEH + "can not go to first event in raw-reader to determine run-id.");}
            runNo = fCurrentData.fRawReader->GetRunNumber();
            Info(kEH, "Determining run-no from raw ... run=%d.", runNo);
            fCurrentData.fRawReader->RewindEvents();
        }
        else
        {
            fEventManager->SetEventId(0);
            return;
        }
    }
    fIsOpen = kTRUE;
}


void AliEveDataSourceOffline::Close()
{
    // Close the event data-files and delete ESD, ESDfriend, run-loader
    // and raw-reader.
    
    cout<<"\n\n\nClose() called!!\n\n\n"<<endl;
    
    static const TEveException kEH("AliEveEventManager::Close ");
    
    if (!fIsOpen)
    {
        throw (kEH + "Event-files not opened.");
    }
    
    if (fEventManager->GetAutoLoadRunning()){
        fEventManager->StopAutoLoadTimer();
    }
    if ((fCurrentData.fESDTree!=0) || (fCurrentData.fHLTESDTree!=0)) {
        delete fCurrentData.fESD;       fCurrentData.fESD       = 0;
        // delete fCurrentData.fESDfriend; // friend tree is deleted with the tree
        fCurrentData.fESDfriend = 0;
        fESDfriendExists = kFALSE;
        
        if(fCurrentData.fESDTree) { delete fCurrentData.fESDTree;   fCurrentData.fESDTree = 0; }
        if(fCurrentData.fHLTESDTree) { delete fCurrentData.fHLTESDTree;   fCurrentData.fHLTESDTree = 0; }
        delete fCurrentData.fESDFile;   fCurrentData.fESDFile = 0;
    }
    
    if (fCurrentData.fAODTree) {
        delete fCurrentData.fAOD;       fCurrentData.fAOD       = 0;
        
        delete fCurrentData.fAODTree;   fCurrentData.fAODTree = 0;
        delete fCurrentData.fAODFile;   fCurrentData.fAODFile = 0;
    }
    
    if (fCurrentData.fRunLoader) {
        delete fCurrentData.fRunLoader; fCurrentData.fRunLoader = 0;
    }
    
    if (fCurrentData.fRawReader) {
        delete fCurrentData.fRawReader; fCurrentData.fRawReader = 0;
    }
    
    fEventManager->SetEventId(-1);
    fIsOpen   = kFALSE;
    fEventManager->SetHasEvent(false);
}

Int_t AliEveDataSourceOffline::GetMaxEventId(Bool_t refreshESD) const
{
    // Returns maximum available event id.
    // If under external control or event is not opened -1 is returned.
    // If raw-data is the only data-source this can not be known
    // and 10,000,000 is returned.
    // If neither data-source is initialised an exception is thrown.
    // If refresh_esd is true and ESD is the primary event-data source
    // its header is re-read from disk.
    
    static const TEveException kEH("AliEveEventManager::GetMaxEventId ");
    
    if (fIsOpen == kFALSE)
    {
        return -1;
    }
    
    if ((fCurrentData.fESDTree!=0) || (fCurrentData.fHLTESDTree!=0))
    {
        if (refreshESD)
        {
            if(fCurrentData.fESDTree!=0) fCurrentData.fESDTree->Refresh();
            if(fCurrentData.fHLTESDTree!=0) fCurrentData.fHLTESDTree->Refresh();
            fEventManager->GetEventSelector()->Update();
        }
        
        Int_t maxEventId=0;
        switch(fgESDvisibleTrees){
            default:
            case AliEveEventManager::kOfflineTree :
                maxEventId = fCurrentData.fESDTree->GetEntries() - 1;
                break;
            case AliEveEventManager::kHLTTree :
                maxEventId = fCurrentData.fHLTESDTree->GetEntries() - 1;
                break;
        }
        
        return maxEventId;
    }
    else if (fCurrentData.fAODTree)
    {
        return fCurrentData.fAODTree->GetEntries() - 1;
    }
    else if (fCurrentData.fRunLoader)
    {
        return fCurrentData.fRunLoader->GetNumberOfEvents() - 1;
    }
    else if (fCurrentData.fRawReader)
    {
        Int_t n = fCurrentData.fRawReader->GetNumberOfEvents() - 1;
        return n > -1 ? n : 10000000;
    }
    else
    {
        throw (kEH + "neither ESD, AOD, RunLoader nor Raw loaded.");
    }
}

/******************************************************************************/

void AliEveDataSourceOffline::SetAssertElements(Bool_t assertRunloader, Bool_t assertEsd,
                                           Bool_t assertAod, Bool_t assertRaw)
{
    // Set global flags that detrmine which parts of the event-data must
    // be present when the event is opened.
    
    fgAssertRunLoader = assertRunloader;
    fgAssertESD = assertEsd;
    fgAssertAOD = assertAod;
    fgAssertRaw = assertRaw;
}

void AliEveDataSourceOffline::SearchRawForCentralReconstruction()
{
    // Enable searching of raw data in standard location. The path passed to
    // Open() is expected to point to a centrally reconstructed run, e.g.:
    // "alien:///alice/data/2009/LHC09c/000101134/ESDs/pass1/09000101134018.10".
    
    fgRawFromStandardLoc = kTRUE;
}

void AliEveDataSourceOffline::SetESDFileName(const TString& esd, EVisibleESDTrees shown)
{
    fgESDvisibleTrees = shown;
    // Set file-name for opening ESD, default "AliESDs.root".
    if (esd.IsNull()) return;
    
    fgESDFileName = esd;
    if (esd.EndsWith(".zip")) fgESDFileName.Form("%s#AliESDs.root",esd.Data());
}

void AliEveDataSourceOffline::SetESDfriendFileName(const TString& esdf)
{
    // Set file-name for opening ESD friend, default "AliESDfriends.root".
    if (esdf.IsNull()) return;
    fgESDfriendsFileName = esdf;
    
    if (esdf.EndsWith(".zip")) fgESDfriendsFileName.Form("%s#AliESDfriends.root",esdf.Data());
}

void AliEveDataSourceOffline::SetAODFileName(const TString& aod)
{
    // Set file-name for opening AOD, default "AliAOD.root".
    if (aod.IsNull()) return;
    fgAODFileName = aod;
    
    if (aod.EndsWith(".zip")) fgAODFileName.Form("%s#AliAOD.root",aod.Data());
}

void AliEveDataSourceOffline::AddAODfriend(const TString& friendFileName)
{
    // Add new AOD friend file-name to be attached when opening AOD.
    // This should include '.root', as in 'AliAOD.VertexingHF.root'.
    
    if (fgAODfriends == 0)
    {
        fgAODfriends = new TList;
        fgAODfriends->SetOwner(kTRUE);
    }
    if (fgAODfriends->FindObject(friendFileName) == 0)
    {
        fgAODfriends->Add(new TObjString(friendFileName));
    }
}

void AliEveDataSourceOffline::SetRawFileName(const TString& raw)
{
    // Set file-name for opening of raw-data, default "raw.root"
    if (raw.IsNull()) return;
    
    fgRawFileName = raw;
}

void AliEveDataSourceOffline::SetGAliceFileName(const TString& galice)
{
    // Set file-name for opening gAlice, default "galice.root".
    
    if ( galice.IsNull()) return;
    fgGAliceFileName = galice;
    
    if (galice.EndsWith(".zip")) fgGAliceFileName.Form("%s#galice.root",galice.Data());
}

void AliEveDataSourceOffline::SetFilesPath(const TString& urlPath)
{
    cout<<"\n\n setting path:"<<urlPath.Data()<<endl;
    
    TString path = urlPath;
    gSystem->ExpandPathName(path);
    if (path.IsNull() || path == ".")
    {
        path = gSystem->WorkingDirectory();
    }
    
    TString sep;
    if(path.EndsWith(".zip")) // if given a path to root_archive.zip
        sep= "#";
    else if(!path.EndsWith("/"))
        sep = "/";
    
    SetESDFileName( TString(Form("%s%sAliESDs.root", path.Data(), sep.Data())) );
    SetESDfriendFileName(  TString(Form("%s%sAliESDfriends.root", path.Data(), sep.Data())) );
    SetAODFileName(  TString(Form("%s%sAliAOD.root", path.Data(), sep.Data())) );
    AddAODfriend(  TString(Form("%s%sAliAOD.VertexingHF.root", path.Data(), sep.Data())) );
    SetGAliceFileName( TString(Form("%s%sgalice.root", path.Data(), sep.Data())) );
    SetRawFileName(TString(Form("%s%sraw.root", path.Data(), sep.Data())));
    
    if(fIsOpen)
    {
        Close(); // close old files
    }
    Open();  // open files with new path
}


TTree* AliEveDataSourceOffline::readESDTree(const char *treeName, int &runNo)
{
    if(!fCurrentData.fESDFile && !fCurrentData.fESD) return 0;
    
    static const TEveException kEH("AliEveEventManager::readESDTree ");
    
    TTree* tempTree = 0;
    
    tempTree =(TTree*) fCurrentData.fESDFile->Get(treeName);
    if (tempTree != 0)
    {
        TFile *esdFriendFile = TFile::Open(fgESDfriendsFileName);
        if (esdFriendFile)
        {
            if (!esdFriendFile->IsZombie())
            {
                esdFriendFile->Close();
                fESDfriendExists = kTRUE;
                tempTree->SetBranchStatus ("ESDfriend*", 1);
            }
            delete esdFriendFile;
        }
        
        fCurrentData.fESD->ReadFromTree(tempTree);
        if (fESDfriendExists)
        {
            fCurrentData.fESDfriend = (AliESDfriend*) fCurrentData.fESD->FindListObject("AliESDfriend");
            Info(kEH, "found and attached ESD friend.");
        }
        else
        {
            Warning(kEH, "ESDfriend not found.");
        }
        
        if (tempTree->GetEntry(0) <= 0)
        {
            Warning(kEH, "failed getting the first entry from tree: %s", treeName);
        }
        else
        {
            if (runNo < 0)
                runNo = fCurrentData.fESD->GetESDRun()->GetRunNumber();
        }
    }
    else // tree == 0
    {
        Warning(kEH, "failed getting the tree:%s", treeName);
    }
    
    return tempTree;
}



