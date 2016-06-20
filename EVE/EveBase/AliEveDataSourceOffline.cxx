//
//  AliEveDataSourceOffline.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#include "AliEveDataSourceOffline.h"
#include "AliEveInit.h"

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
fIsOpen(false),
fESDfriendExists(kFALSE),
fgAODfriends(0),
fgRawFromStandardLoc(false)
{
    fEventManager = AliEveEventManager::Instance();
    Open();
}

AliEveDataSourceOffline::~AliEveDataSourceOffline()
{
    
}

void AliEveDataSourceOffline::NextEvent()
{
    GotoEvent(AliEveEventManager::CurrentEventId()+1);
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
    
    if(!fCurrentData.fESD && !fCurrentData.fAOD && !fCurrentData.fRawReader)
    {
        AliWarning("No ESD event avaliable. Probably files were not opened.");
        return;
    }
    
    if(fCurrentData.fESD)
    {
        if(fCurrentData.fESD->GetRunNumber() != fEventManager->GetCurrentRun())
        {
            fEventManager->ResetMagneticField();
            fEventManager->SetCurrentRun(fCurrentData.fESD->GetRunNumber());
        }
    }
    else if(fCurrentData.fAOD)
    {
        if(fCurrentData.fAOD->GetRunNumber() != fEventManager->GetCurrentRun())
        {
            fEventManager->ResetMagneticField();
            fEventManager->SetCurrentRun(fCurrentData.fAOD->GetRunNumber());
        }
    }
    
    if (!fIsOpen){
        AliFatal("Event-files not opened but ED is in offline mode.");
    }
    
    fEventManager->SetHasEvent(false);
    
    Int_t maxEvent = 0;
    if ((fCurrentData.fESDTree!=0))
    {
        if(fCurrentData.fESDTree){
            if (event >= fCurrentData.fESDTree->GetEntries())
                fCurrentData.fESDTree->Refresh();
            maxEvent = fCurrentData.fESDTree->GetEntries() - 1;
            if (event < 0)
                event = fCurrentData.fESDTree->GetEntries() + event;
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
                AliError("current raw-data source does not support direct event access.");
                return;
            }
            AliInfo("number of events unknown for current raw-data source, setting max-event id to 10M.");
        }
        else
        {
            if (event < 0)
                event = fCurrentData.fRawReader->GetNumberOfEvents() + event;
        }
    }
    else
    {
        AliFatal("neither RunLoader, ESD nor Raw loaded.");
    }
    if (event < 0)
    {
        AliFatal(Form("event %d not present, available range [0, %d].",event, maxEvent));
    }
    if (event > maxEvent)
    {
        AliInfo("Event number out of range. Going to event 0");
        GotoEvent(0);
    }

    TString sysInfoHeader;
    sysInfoHeader.Form("AliEveEventManager::GotoEvent(%d) - ", event);
    AliSysInfo::AddStamp(sysInfoHeader + "Start");

    fEventManager->DestroyTransients();
    
    AliSysInfo::AddStamp(sysInfoHeader + "PostDestroy");
    
    if (fCurrentData.fESDTree) {
        if (fCurrentData.fESDTree->GetEntry(event) <= 0)
            AliFatal("failed getting required event from ESD.");
        
        if (fESDfriendExists)
            fCurrentData.fESD->SetESDfriend(fCurrentData.fESDfriend);
    }
    if (fCurrentData.fAODTree) {
        if (fCurrentData.fAODTree->GetEntry(event) <= 0)
            AliFatal("failed getting required event from AOD.");
    }
    
    if (fCurrentData.fRunLoader) {
        if (fCurrentData.fRunLoader->GetEvent(event) != 0)
            AliFatal("failed getting required event.");
    }
    if (fCurrentData.fRawReader)
    {
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
                    AliFatal(Form("Error going to next raw-event from event %d.", rawEv));
                }
                ++rawEv;
            }
            AliWarning(Form("Loaded raw-event %d with fallback method.\n", rawEv));
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

void AliEveDataSourceOffline::Open()
{
    // Open event-data from URL specified in path.
    // Attempts to create AliRunLoader() and to open ESD with ESDfriends.
    // Warning is reported if run-loader or ESD is not found.
    
    if (fIsOpen){
        AliError("Event-files already opened.");
        return;
    }
    
    Int_t runNo = -1;

    // Open ESD and ESDfriends
    if ((fCurrentData.fESDFile = TFile::Open(fgESDFileName)))
    {
        fCurrentData.fESD = new AliESDEvent();
        fCurrentData.fESDTree = readESDTree("esdTree", runNo);
        
        
        if(!fCurrentData.fESDTree){
            // both ESD trees are == 0
            delete fCurrentData.fESDFile; fCurrentData.fESDFile = 0;
            delete fCurrentData.fESD; fCurrentData.fESD = 0;
        }
        
        
    }
    else
    {
        AliWarning(Form("can not read ESD file '%s'.", fgESDFileName.Data()));
    }
    if (fCurrentData.fESDTree == 0)
    {
        AliWarning("ESD not initialized.");
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
                AliWarning("failed getting the first entry from addTree.");
            }
            else if (runNo < 0){runNo = fCurrentData.fAOD->GetRunNumber();}
        }
        else // aodtree == 0
        {
            delete fCurrentData.fAODFile; fCurrentData.fAODFile = 0;
            delete fCurrentData.fAOD;     fCurrentData.fAOD     = 0;
            AliWarning("failed getting the aodTree.");
        }
    }
    else // aod not readable
    {
        AliWarning(Form("can not read AOD file '%s'.", fgAODFileName.Data()));
    }
    if (fCurrentData.fAODTree == 0)
    {
        AliWarning("AOD not initialized.");
    }
    // Open RunLoader from galice.root
    //    fgGAliceFileName = "/Users/Jerus/galice.root"; // temp
    
    TFile *gafile = TFile::Open(fgGAliceFileName);
    
    if (gafile)
    {
        gafile->Close();
        delete gafile;

        fCurrentData.fRunLoader = AliRunLoader::Open(fgGAliceFileName, fEventManager->GetName());
        if (fCurrentData.fRunLoader)
        {
            TString alicePath(gSystem->DirName(fgGAliceFileName));
            alicePath.Append("/");
            fCurrentData.fRunLoader->SetDirName(alicePath);
            
            if (fCurrentData.fRunLoader->LoadgAlice() != 0)
            {
                AliWarning("failed loading gAlice via run-loader.");
            }
            if (fCurrentData.fRunLoader->LoadHeader() == 0){
                if(runNo < 0){
                    runNo = fCurrentData.fRunLoader->GetHeader()->GetRun();
                }
            }
            else{
                AliWarning("failed loading run-loader's header.");
                delete fCurrentData.fRunLoader;
                fCurrentData.fRunLoader = 0;
            }
        }
        else
        {
            AliWarning(Form("failed opening ALICE run-loader from '%s'.", fgGAliceFileName.Data()));
        }
    }
    else
    {
        AliWarning(Form("can not read '%s'.", fgGAliceFileName.Data()));
    }
    
    if (fCurrentData.fRunLoader == 0)
    {
        AliWarning("Bootstraping of run-loader failed.");
    }
    // Open raw-data file
    TString rawPath;
    if (fgRawFromStandardLoc)
    {
        if (!fgRawFileName.BeginsWith("alien:")){
            AliFatal("Standard raw search requested, but the directory is not in AliEn.");
        }
        if (!fgRawFileName.Contains("/ESDs/")){
            AliFatal("Standard raw search requested, but does not contain 'ESDs' directory.");
        }
        
        TPMERegexp chunk("/([\\d\\.])+/?$");
        Int_t nm = chunk.Match(fgRawFileName);
        if (nm != 2){
            AliFatal("Standard raw search requested, but the path does not end with chunk-id directory.");
        }
        
        TPMERegexp esdstrip("/ESDs/.*");
        rawPath = fgRawFileName;
        esdstrip.Substitute(rawPath, "/raw/");
        rawPath += chunk[0];
        rawPath += ".root";
        
        AliInfo(Form("Standard raw search requested, using the following path:\n  %s\n", rawPath.Data()));
    }
    else
    {
        rawPath = fgRawFileName;
    }
    // If i use open directly, raw-reader reports an error but i have
    // no way to detect it.
    // Is this (AccessPathName check) ok for xrootd / alien? Yes, not for http.
    AliLog::EType_t oldLogLevel = (AliLog::EType_t) AliLog::GetGlobalLogLevel();
    AliLog::SetGlobalLogLevel(AliLog::kFatal);
    
    if (gSystem->AccessPathName(rawPath, kReadPermission) == kFALSE){
        fCurrentData.fRawReader = AliRawReader::Create(rawPath);
    }
    else{
        fCurrentData.fRawReader = AliRawReader::Create(fgRawFileName);
    }
    
    AliLog::SetGlobalLogLevel(oldLogLevel);
    
    if (fCurrentData.fRawReader == 0)
    {
        AliWarning("raw-data not initialized.");
    }
    if (runNo < 0)
    {
        if (fCurrentData.fRawReader)
        {
            if (!fCurrentData.fRawReader->NextEvent())
            {
                AliFatal("can not go to first event in raw-reader to determine run-id.");
            }
            runNo = fCurrentData.fRawReader->GetRunNumber();
            AliInfo(Form("Determining run-no from raw ... run=%d.", runNo));
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
    
    if (!fIsOpen)
    {
        AliError("files were are not opened");
        return;
    }
    
    if (fEventManager->GetAutoLoadRunning()){
        fEventManager->StopAutoLoadTimer();
    }
    if (fCurrentData.fESDTree!=0) {
        delete fCurrentData.fESD;       fCurrentData.fESD       = 0;
        // delete fCurrentData.fESDfriend; // friend tree is deleted with the tree
        fCurrentData.fESDfriend = 0;
        fESDfriendExists = kFALSE;
        
        if(fCurrentData.fESDTree) { delete fCurrentData.fESDTree;   fCurrentData.fESDTree = 0; }
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
    
    if (!fIsOpen){return -1;}
    
    if (fCurrentData.fESDTree!=0)
    {
        if (refreshESD)
        {
            fCurrentData.fESDTree->Refresh();
            fEventManager->GetEventSelector()->Update();
        }
        return fCurrentData.fESDTree->GetEntries() - 1;
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
        AliFatal("neither ESD, AOD, RunLoader nor Raw loaded.");
    }
}

/******************************************************************************/

void AliEveDataSourceOffline::SearchRawForCentralReconstruction()
{
    // Enable searching of raw data in standard location. The path passed to
    // Open() is expected to point to a centrally reconstructed run, e.g.:
    // "alien:///alice/data/2009/LHC09c/000101134/ESDs/pass1/09000101134018.10".
    
    fgRawFromStandardLoc = kTRUE;
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

void AliEveDataSourceOffline::SetFilesPath(const TString& urlPath)
{
    AliInfo(Form("setting path:%s",urlPath.Data()));
    
    TString path = urlPath;
    gSystem->ExpandPathName(path);
    if (path.IsNull() || path == ".")
    {
        path = gSystem->WorkingDirectory();
    }
    
    fgESDFileName        = TString::Format("%s/AliESDs.root", path.Data());
    fgESDfriendsFileName = TString::Format("%s/AliESDfriends.root", path.Data());
    fgAODFileName        = TString::Format("%s/AliAOD.root", path.Data());
    fgRawFileName        = TString::Format("%s/raw.root", path.Data());
    fgGAliceFileName     = TString::Format("%s/galice.root", path.Data());

    AddAODfriend(TString(Form("%s/AliAOD.VertexingHF.root", path.Data())));
    
    if(fIsOpen){Close();} // close old files
    Open();  // open files with new path
}


TTree* AliEveDataSourceOffline::readESDTree(const char *treeName, int &runNo)
{
    if(!fCurrentData.fESDFile && !fCurrentData.fESD) return 0;
    
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
            AliInfo("found and attached ESD friend.");
        }
        else
        {
            AliWarning("ESDfriend not found.");
        }
        
        if (tempTree->GetEntry(0) <= 0)
        {
            AliWarning(Form("failed getting the first entry from tree: %s", treeName));
        }
        else
        {
            if (runNo < 0)
                runNo = fCurrentData.fESD->GetESDRun()->GetRunNumber();
        }
    }
    else // tree == 0
    {
        AliWarning(Form("failed getting the tree:%s", treeName));
    }
    
    return tempTree;
}



