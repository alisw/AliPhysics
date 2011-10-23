#ifndef ALIRUNLoader_H
#define ALIRUNLoader_H

//___________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliRunLoader                                             //
//                                                                 //
//  This class aims to be the only one interface for manging data  //
//  It stores Loaders for all modules which knows the filenames    //
//  of the data files to be stored.                                //
//  It aims to substitude AliRun in automatic managing of data     //
//  positioning thus there won't be necessity of loading gAlice    //
//  from file in order to get fast access to the data              //
//                                                                 //
//  Logical place to put the specific Loader to the given          //
//  detector is detector  itself (i.e ITSLoader in ITS).           //
//  But, to load detector object one need to load gAlice, and      //
//  by the way all other detectors with their geometrieces and     //
//  so on. So, if one need to open TPC clusters there is no        //
//  principal nedd to read everything.                             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class TFile;
class TString;
class TFolder;
class TObjArray;
class TTree;
class TParticle;

class AliRun;
class AliLoader;
class AliDetector;
class AliHeader;
class AliStack;
class AliCDBEntry;
class AliCentralTrigger;

#include <TNamed.h>

#include "AliConfig.h"
#include "AliDataLoader.h"

class AliRunLoader: public TNamed
{
  public:
    
    AliRunLoader();
    AliRunLoader(const char* topfoldername);
    AliRunLoader(TFolder* topfolder);
    
    virtual ~AliRunLoader();
    
    static AliRunLoader* Open(const char* filename = "galice.root",
                              const char* eventfoldername = AliConfig::GetDefaultEventFolderName(),
	          Option_t* option = "READ");

    Int_t       GetEventNumber() const {return fCurrentEvent;}

    Int_t       GetEvent(Int_t evno);//sets the event number and reloads data in folders properly
    Int_t       GetNextEvent(){return GetEvent(fCurrentEvent+1);}//gets next event 
    Int_t       SetEventNumber(Int_t evno); //cleans folders and sets the root dirs in files (do not reload data)
    Int_t       SetNextEvent(){return SetEventNumber(fCurrentEvent+1);}
    
    Int_t       GetNumberOfEvents();
    
    AliCDBEntry* GetCDBEntry(const char* name) const;

    void        MakeTree(Option_t *option);
    void        MakeHeader();
    void        MakeTrigger();
    void        MakeStack();
    
    Int_t       LoadgAlice();
    Int_t       LoadHeader();
    Int_t       LoadKinematics(Option_t* option = "READ");
    Int_t       LoadTrigger(Option_t* option = "READ");
    Int_t       LoadTrackRefs(Option_t* option = "READ");
    
    void        UnloadHeader();
    void        UnloadTrigger();
    void        UnloadKinematics();
    void        UnloadgAlice();
    void        UnloadTrackRefs();
    
    void        SetKineFileName(const TString& fname){fKineDataLoader->SetFileName(fname);}
    void        SetTrackRefsFileName(const TString& fname){fTrackRefsDataLoader->SetFileName(fname);}
    
    TTree*      TreeE() const; //returns the tree from folder; shortcut method
    TTree*      TreeCT() const; //returns the tree from folder; shortcut method
    AliHeader*  GetHeader() const;
    AliCentralTrigger*  GetTrigger() const;
    
    AliStack*   Stack() const {return fStack;}
    
    TTree*      TreeK() const; //returns the tree from folder; shortcut method
    TTree*      TreeTR() const; //returns the tree from folder; shortcut method    
    
    AliRun*     GetAliRun()const;
    Int_t       GetRunNumber() const {return fRun;}
    void        SetRunNumber(Int_t run) {fRun=run;}
        
    Int_t       WriteHeader(Option_t* opt="");
    Int_t       WriteTrigger(Option_t* opt="");
    Int_t       WriteAliRun(Option_t* opt="");
    Int_t       WriteKinematics(Option_t* opt="");
    Int_t       WriteTrackRefs(Option_t* opt="");
    Int_t       WriteRunLoader(Option_t* opt="");
    
    Int_t       WriteHits(Option_t* opt=""); 
    Int_t       WriteSDigits(Option_t* opt="");
    Int_t       WriteDigits(Option_t* opt="");
    Int_t       WriteRecPoints(Option_t* opt="");
    Int_t       WriteTracks(Option_t* opt="");
    
    Int_t       LoadHits(Option_t* detectors = "all",Option_t* opt = "READ");
    Int_t       LoadSDigits(Option_t* detectors = "all",Option_t* opt = "READ");
    Int_t       LoadDigits(Option_t* detectors = "all",Option_t* opt = "READ");
    Int_t       LoadRecPoints(Option_t* detectors = "all",Option_t* opt = "READ");
    Int_t       LoadTracks(Option_t* detectors = "all",Option_t* opt = "READ");
    Int_t       LoadRecParticles(Option_t* detectors = "all",Option_t* opt = "READ");
    
    void        UnloadHits(Option_t* detectors = "all");
    void        UnloadSDigits(Option_t* detectors = "all");
    void        UnloadDigits(Option_t* detectors = "all");
    void        UnloadRecPoints(Option_t* detectors = "all");
    void        UnloadTracks(Option_t* detectors = "all");
    void        UnloadRecParticles(Option_t* detectors = "all");
    void        UnloadAll(Option_t* detectors = "all");    
    
    void        AddLoader(AliLoader* loader);
    void        AddLoader(AliDetector* det);
    AliLoader*  GetLoader(const char* detname) const;
    AliLoader*  GetLoader(AliDetector* det) const;
    Int_t       SetEventFolderName(const TString& name = AliConfig::GetDefaultEventFolderName());//sets top folder name for this run; of alread
    void        CleanFolders();//removes all abjects from folder structure
    void        CleanDetectors();
    void        CleanKinematics(){Clean(fgkKineContainerName);}
    void        CleanTrackRefs(){Clean(fgkTrackRefsContainerName);}
    
    void        RemoveEventFolder(); //remove folder structure from top folder 
    void        SetCompressionLevel(Int_t cl);
    void        SetKineComprLevel(Int_t cl);
    void        SetTrackRefsComprLevel(Int_t cl);

    TFolder*    GetEventFolder() const {return fEventFolder;}
    void        CdGAFile();

    void        MakeTrackRefsContainer();
    void        SetDirName(TString& dirname);
    Int_t       GetFileOffset() const;
    void        SetNumberOfEventsPerFile(Int_t nevpf){fNEventsPerFile = nevpf;}
    void        SetNumberOfEventsPerRun(Int_t nevpr) {fNEventsPerRun = nevpr;}
    Int_t       GetNumberOfEventsPerRun()            {return fNEventsPerRun;}
    void        SetDigitsFileNameSuffix(const TString& suffix);//adds the suffix before ".root", 
                                                               //e.g. TPC.Digits.root -> TPC.DigitsMerged.root
                                                               //made on Jiri Chudoba demand
    TString     GetFileName() const;//returns name of galice file
    const TObjArray* GetArrayOfLoaders() const {return fLoaders;}
    void Cd(){fgRunLoader = this;}
    void Synchronize();
    
    AliLoader*    GetDetectorLoader(const char* detname);
    TTree*        GetTreeH(const char* detname, Bool_t maketree);
    TTree*        GetTreeS(const char* detname, Bool_t maketree);
    TTree*        GetTreeD(const char* detname, Bool_t maketree);
    TTree*        GetTreeR(const char* detname, Bool_t maketree);
    TTree*        GetTreeT(const char* detname, Bool_t maketree);
    TTree*        GetTreeP(const char* detname, Bool_t maketree);

  /******************************************/
  /*****   Public S T A T I C Stuff   *******/
  /******************************************/
    static AliRunLoader* GetRunLoader(const char* eventfoldername);
    static AliRunLoader* Instance(){return fgRunLoader;}
    static AliLoader*    GetDetectorLoader(const char* detname, const char* eventfoldername);
    static TTree*        GetTreeH(const char* detname, Bool_t maketree, const char* eventfoldername);
    static TTree*        GetTreeS(const char* detname, Bool_t maketree, const char* eventfoldername);
    static TTree*        GetTreeD(const char* detname, Bool_t maketree, const char* eventfoldername);
    static TTree*        GetTreeR(const char* detname, Bool_t maketree, const char* eventfoldername);
    static TTree*        GetTreeT(const char* detname, Bool_t maketree, const char* eventfoldername);
    static TTree*        GetTreeP(const char* detname, Bool_t maketree, const char* eventfoldername);

    static TString GetRunLoaderName () {return fgkRunLoaderName;}
    static TString GetHeaderContainerName () {return fgkHeaderContainerName;}
    static TString GetTriggerContainerName () {return fgkTriggerContainerName;}
    static TString GetKineContainerName () {return fgkKineContainerName;}
    static TString GetTrackRefsContainerName () {return fgkTrackRefsContainerName;}
    static TString GetHeaderBranchName () {return fgkHeaderBranchName;}
    static TString GetTriggerBranchName () {return fgkTriggerBranchName;}
    static TString GetKineBranchName () {return fgkKineBranchName;}
    static TString GetTriggerFileName() { return fgkDefaultTriggerFileName; }
    static TString GetGAliceName () {return fgkGAliceName;}
     
protected:
    void               SetGAliceFile(TFile* gafile);//sets the pointer to gAlice file
    Int_t              OpenKineFile(Option_t* opt);
    Int_t              OpenTrackRefsFile(Option_t* opt);

    Int_t              OpenDataFile(const TString& filename,TFile*& file,TDirectory*& dir,Option_t* opt,Int_t cl);
    void               SetUnixDir(const TString& udirname);
    const TString      SetFileOffset(const TString& fname);//adds the proper number before .root
    void               SetDetectorAddresses();

    TObjArray         *fLoaders;          //  List of Detectors
    TFolder           *fEventFolder;      //!top folder for this run
    
    Int_t              fRun;               //! Current run number
    Int_t              fCurrentEvent;//!Number of current event
    
    TFile             *fGAFile;//!  pointer to main file with AliRun and Run Loader -> galice.root 
    AliHeader         *fHeader;//!  pointer to header
    AliStack          *fStack; //!  pointer to stack
    AliCentralTrigger *fCTrigger; //! pointer to CEntral Trigger Processor
    
    AliDataLoader     *fKineDataLoader;// kinematics data loader
    AliDataLoader     *fTrackRefsDataLoader;//track reference data loader
    
    Int_t              fNEventsPerFile;  //defines number of events stored per one file
    Int_t              fNEventsPerRun;   //defines number of event per run
    TString            fUnixDirName;     //! name of unix path to directory that contains event
    static const TString   fgkDefaultKineFileName;//default file name with kinamatics
    static const TString   fgkDefaultTrackRefsFileName;//default file name with kinamatics
    static const TString   fgkDefaultTriggerFileName;//default file name with trigger

        
  private:
    AliRunLoader(const AliRunLoader &r);      //Not implemented
    AliRunLoader & operator = (const AliRunLoader &); //Not implemented 
    void  GetListOfDetectors(const char * namelist,TObjArray& pointerarray) const;
    void  CleanHeader(){Clean(fgkHeaderContainerName);}
    void  CleanTrigger(){Clean(fgkTriggerContainerName);}
    void  Clean(const TString& name);
    Int_t SetEvent();
   
    static AliRunLoader*   fgRunLoader; //pointer to the AliRunLoader instance

    static const TString   fgkRunLoaderName;          //default name of the run loader
    static const TString   fgkHeaderContainerName;    //default name of the kinematics container (TREE) name - TreeE
    static const TString   fgkTriggerContainerName;   //default name of the trigger container (TREE) name - TreeCT
    static const TString   fgkKineContainerName;      //default name of the kinematics container (TREE) name - TreeK
    static const TString   fgkTrackRefsContainerName; //default name of the track references container (TREE) name - TreeTR
    static const TString   fgkHeaderBranchName;       //default name of the branch containing the header
    static const TString   fgkTriggerBranchName;      //default name of the branch containing the trigger
    static const TString   fgkKineBranchName;         //default name of the branch with kinematics
    static const TString   fgkGAliceName;             //default name for gAlice file    
    
    ClassDef(AliRunLoader,3)
};

#endif
