//
//  AliEveDataSourceOffline.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#ifndef __AliEveDataSourceOffline__
#define __AliEveDataSourceOffline__

#include "AliEveEventManager.h"
#include "AliEveDataSource.h"
#include "AliEventInfo.h"

class AliEveDataSourceOffline : public AliEveDataSource
{
public:
    AliEveDataSourceOffline(bool storageManager=false);
    ~AliEveDataSourceOffline();
    
    enum EVisibleESDTrees{ kOfflineTree, kHLTTree };
    
    void SetAssertElements(Bool_t assertRunloader, Bool_t assertEsd,
                                  Bool_t assertAod, Bool_t assertRaw);
    
    void SearchRawForCentralReconstruction();
    
    void SetESDFileName(const TString& esd, EVisibleESDTrees shown=kOfflineTree);
    void SetESDfriendFileName(const TString& esdf);
    void SetAODFileName(const TString& aod);
    void AddAODfriend  (const TString& friendFileName);
    void SetRawFileName(const TString& raw);
    void SetGAliceFileName(const TString& galice);
    void SetFilesPath(const TString& path);

    void GotoEvent(Int_t event);
        
    virtual Int_t GetMaxEventId(Bool_t refreshESD=kFALSE) const;
private:
    void Init(){};// to implement
    void SetEvent(AliRunLoader *runLoader, AliRawReader *rawReader, AliESDEvent *esd, AliESDfriend *esdf);
    void NextEvent();
    void PrevEvent();
    void Open();
    void Close();
    TTree* readESDTree(const char* treeName, int &run);
    
    bool            fIsOpen;         // Are event-files opened.
    AliEventInfo	fEventInfo;		// Current Event Info
    Bool_t        fESDfriendExists;	// Flag specifying if ESDfriend was found during opening of the event-data.
    
    // names of files:
    TString  fgGAliceFileName;        // galice.root file
    TString  fgESDFileName;           // Name by which to open ESD.
    TString  fgESDfriendsFileName;
    TString  fgAODFileName;           // Name by which to open AOD.
    TString  fgRawFileName;           // Name by which to open raw-data file.
    
    AliEveEventManager *fEventManager;
    
    EVisibleESDTrees  fgESDvisibleTrees; // trees to open from ESD
    
    Bool_t   fgAssertRunLoader;	// Global flag specifying if AliRunLoader must be asserted during opening of the event-data.
    
    Bool_t   fgAssertESD;		// Global flag specifying if ESDEvent must be asserted during opening of the event-data.
    Bool_t   fgAssertAOD;		// Global flag specifying if AODEvent must be asserted during opening of the event-data.
    Bool_t   fgAssertRaw;		// Global flag specifying if raw-data presence must be asserted during opening of the event-data.
    
    TList   *fgAODfriends;         // Global list of AOD friend names to be attached during opening of the event-data (empty by default).
    
    Bool_t   fgRawFromStandardLoc; // Global flag to enable looking for raw data in ../../../raw/, as it is stored for central reco.

    AliEveDataSourceOffline(const AliEveDataSourceOffline&);
    AliEveDataSourceOffline& operator=(const AliEveDataSourceOffline&);
    
   ClassDef(AliEveDataSourceOffline, 0);
};


#endif
