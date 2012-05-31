#ifndef ALILOADER_H
#define ALILOADER_H

class TDirectory;
class TFile;
class TString;
class TTree;

#include <TFolder.h>
#include <TObjArray.h>

class AliRunLoader;

#include "AliDataLoader.h"
#include "AliBaseLoader.h"


//___________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliLoader                                                //
//                                                                 //
//  Base class for Loaders.                                        //
//  Loader provides base I/O fascilities for "standard" data.      //
//  Each detector has a laoder data member.                        //
//  Loader is always accessible via folder structure as well       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliLoader: public TNamed
 {
   public:
    AliLoader();
    AliLoader(const Char_t *detname,const Char_t *eventfoldername); //contructor with name of the top folder of the tree
    AliLoader(const Char_t *detname,TFolder* eventfolder);

    virtual ~AliLoader();//----------------- 
 
    virtual Int_t  GetEvent();//changes to root directory
                             //relies on the event number defined in gAlice
    virtual Int_t  SetEvent();//basically the same that GetEvent but do not post data to folders
    virtual void   MakeTree(Option_t* opt);

    //these methods are responsible for reading data
    // 1. Opens the file (accordingly to gAlice->GetEvNumber()/Nevents per File
    // 2. Changes to proper ROOT directory
    // 3. Retrives data and posts to folders
    //options: 
    //   ""       - readonly file
    //   "update" - update file
    //   "recreate" 
    AliDataLoader* GetHitsDataLoader()     const {return (AliDataLoader*)fDataLoaders->At(kHits);}
    AliDataLoader* GetSDigitsDataLoader()  const {return (AliDataLoader*)fDataLoaders->At(kSDigits);}
    AliDataLoader* GetDigitsDataLoader()   const {return (AliDataLoader*)fDataLoaders->At(kDigits);}
    AliDataLoader* GetRecPointsDataLoader()const {return (AliDataLoader*)fDataLoaders->At(kRecPoints);}
    AliDataLoader* GetTracksDataLoader()   const {return (AliDataLoader*)fDataLoaders->At(kTracks);}
    AliDataLoader* GetRecParticlesDataLoader()const {return (AliDataLoader*)fDataLoaders->At(kRecParticles);}
 
    AliDataLoader* GetDataLoader(const char* name);
    void           AddDataLoader(AliDataLoader* dl);
        
    
    Int_t          SetEventFolder(TFolder* eventfolder);//sets the event folder
    Int_t          Register(TFolder* eventFolder);//triggers creation of subfolders for a given detector
      
    TFolder*       GetTopFolder();  //returns top aliroot folder
    TFolder*       GetEventFolder();//returns the folder that event is sitting
    TFolder*       GetDataFolder(); //returns the folder that hits, sdigits, etc are sitting (this contains folders named like detectors)
    TFolder*       GetDetectorDataFolder();//returns the folder that hits, sdigits, etc 
                                     //are sitting for a given detector (subfolder of Data Folder)
    TFolder*       GetModulesFolder();
    
    virtual void   MakeHitsContainer() const {GetHitsDataLoader()->MakeTree();}
    virtual void   MakeSDigitsContainer() const {GetSDigitsDataLoader()->MakeTree();}
    virtual void   MakeDigitsContainer() const {GetDigitsDataLoader()->MakeTree();}
    virtual void   MakeRecPointsContainer() const {GetRecPointsDataLoader()->MakeTree();}
    virtual void   MakeTracksContainer() const {GetTracksDataLoader()->MakeTree();}
    virtual void   MakeRecParticlesContainer() const {GetRecParticlesDataLoader()->MakeTree();}
        
    virtual void   CleanFolders();
    virtual void   CloseFiles();

    TTree*         TreeH() const
      {return GetHitsDataLoader()->Tree();}      //returns the tree from folder; shortcut method
    TTree*         TreeS() const
      {return GetSDigitsDataLoader()->Tree();}   //returns the tree from folder; shortcut method
    TTree*         TreeD() const
      {return GetDigitsDataLoader()->Tree();}    //returns the tree from folder; shortcut method
    TTree*         TreeR() const
      {return GetRecPointsDataLoader()->Tree();} //returns the tree from folder; shortcut method
    TTree*         TreeT() const
      {return GetTracksDataLoader()->Tree();}    //returns the tree from folder; shortcut method
    TTree*         TreeP() const
      {return GetRecParticlesDataLoader()->Tree();} //returns the tree from folder; shortcut method

    Int_t          LoadHits(Option_t* opt=""){
      Int_t status = GetHitsDataLoader()->Load(opt);
      SetTAddrInDet();
      return status;
    }
    Int_t          LoadSDigits(Option_t* opt=""){
      Int_t status = GetSDigitsDataLoader()->Load(opt);
      SetTAddrInDet();
      return status;
    }
    Int_t          LoadDigits(Option_t* opt=""){
      Int_t status = GetDigitsDataLoader()->Load(opt);
      SetTAddrInDet();
      return status;
    }
    

    Int_t          LoadRecPoints(Option_t* opt="") {
      Int_t status = GetRecPointsDataLoader()->Load(opt);
      SetTAddrInDet();
      return status;
    }
    Int_t          LoadTracks(Option_t* opt="") {
      Int_t status = GetTracksDataLoader()->Load(opt);
      SetTAddrInDet();
      return status;
    }
    Int_t          LoadRecParticles(Option_t* opt="") {
      Int_t status = GetRecParticlesDataLoader()->Load(opt);
      SetTAddrInDet();
      return status;
    }

    void           UnloadHits() const {GetHitsDataLoader()->Unload();}
    void           UnloadSDigits() const {GetSDigitsDataLoader()->Unload();}
    void           UnloadDigits() const {GetDigitsDataLoader()->Unload();}
    void           UnloadRecPoints() const{GetRecPointsDataLoader()->Unload();}
    void           UnloadTracks() const {GetTracksDataLoader()->Unload();}
    void           UnloadRecParticles() const {GetRecParticlesDataLoader()->Unload();}
    void           UnloadAll();

    virtual Int_t  ReloadHits() const {return GetHitsDataLoader()->Reload();}  //unload and load again Hits
    virtual Int_t  ReloadSDigits() const {return GetSDigitsDataLoader()->Reload();} //unload and load again 
    virtual Int_t  ReloadDigits() const {return GetDigitsDataLoader()->Reload();} //unload and load again 
    virtual Int_t  ReloadRecPoints() const {return GetRecPointsDataLoader()->Reload();} //unload and load again 
    virtual Int_t  ReloadTracks() const {return GetTracksDataLoader()->Reload();} //unload and load again 
    virtual Int_t  ReloadRecParticles() const {return GetRecParticlesDataLoader()->Reload();} //unload and load again 
    virtual Int_t  ReloadAll(); //unload and load again everything that was loaded 
    
  //these methods writes object from folder to proper file
    virtual Int_t  WriteHits(Option_t* opt="") const;
    virtual Int_t  WriteSDigits(Option_t* opt="") const;
    virtual Int_t  WriteDigits(Option_t* opt="") const {return GetDigitsDataLoader()->WriteData(opt);}
    virtual Int_t  WriteRecPoints(Option_t* opt="") const {return GetRecPointsDataLoader()->WriteData(opt);}
    virtual Int_t  WriteTracks(Option_t* opt="") const {return GetTracksDataLoader()->WriteData(opt);}
    virtual Int_t  WriteRecParticles(Option_t* opt="") const {return GetRecParticlesDataLoader()->WriteData(opt);}

    
  //void SetTopFolder(TString& str){}; //Sets root top folder for the run
    
    void           SetHitsFileName(const TString& fname) const 
      {GetHitsDataLoader()->SetFileName(fname);}
    void           SetSDigitsFileName(const TString& fname) const 
      {GetSDigitsDataLoader()->SetFileName(fname);}
    void           SetDigitsFileName(const TString& fname) const
      {GetDigitsDataLoader()->SetFileName(fname);}
    void           SetRecPointsFileName(const TString& fname) const
      {GetRecPointsDataLoader()->SetFileName(fname);}
    void           SetTracksFileName(const TString& fname) const
      {GetTracksDataLoader()->SetFileName(fname);}
    void           SetRecParticlesFileName(const TString& fname) const
      {GetRecParticlesDataLoader()->SetFileName(fname);}

    const TString& GetHitsFileName() const 
      {return GetHitsDataLoader()->GetFileName();}
    const TString& GetSDigitsFileName() const 
      {return GetSDigitsDataLoader()->GetFileName();}
    const TString& GetDigitsFileName() const
      {return GetDigitsDataLoader()->GetFileName();}
    const TString& GetRecPointsFileName() const
      {return GetRecPointsDataLoader()->GetFileName();}
    const TString& GetTracksFileName() const
      {return GetTracksDataLoader()->GetFileName();}
    const TString& GetRecParticlesFileName() const
      {return GetRecParticlesDataLoader()->GetFileName();}
   
    virtual void  CleanHits() const
      {GetHitsDataLoader()->Clean();}       //cleans hits from folder
    virtual void  CleanSDigits() const
      {GetSDigitsDataLoader()->Clean();}    //cleans digits from folder
    virtual void  CleanDigits() const
      {GetDigitsDataLoader()->Clean();}     //cleans sdigits from folder
    virtual void  CleanRecPoints() const
      {GetRecPointsDataLoader()->Clean();}  //cleans rec. points from folder
    virtual void  CleanTracks() const
      {GetTracksDataLoader()->Clean();}     //cleans tracks from folder

    virtual void  SetHitsFileOption(Option_t* newopt) const
      {GetHitsDataLoader()->SetFileOption(newopt);}          //Sets Hits File Option in open
    virtual void  SetSDigitsFileOption(Option_t* newopt) const
      {GetSDigitsDataLoader()->SetFileOption(newopt);}    //Sets S. Digits File Option in open
    virtual void  SetDigitsFileOption(Option_t* newopt) const
      {GetDigitsDataLoader()->SetFileOption(newopt);}      //Sets Digits File Option in open
    virtual void  SetRecPointsFileOption(Option_t* newopt) const
      {GetRecPointsDataLoader()->SetFileOption(newopt);}//Sets Rec Ponoints File Option in open
    virtual void  SetTracksFileOption(Option_t* newopt) const
      {GetTracksDataLoader()->SetFileOption(newopt);}      //Sets Tracks File Option in open
    virtual void  SetRecParticlesFileOption(Option_t* newopt) const
      {GetRecParticlesDataLoader()->SetFileOption(newopt);}//Sets Rec Ponoints File Option in open
    
    virtual void  SetHitsComprLevel(Int_t cl) const
      {GetHitsDataLoader()->SetCompressionLevel(cl);}
    virtual void  SetDigitsComprLevel(Int_t cl) const
      {GetSDigitsDataLoader()->SetCompressionLevel(cl);}
    virtual void  SetSDigitsComprLevel(Int_t cl) const
      {GetDigitsDataLoader()->SetCompressionLevel(cl);}
    virtual void  SetRecPointsComprLevel(Int_t cl) const
      {GetRecPointsDataLoader()->SetCompressionLevel(cl);}
    virtual void  SetTracksComprLevel(Int_t cl) const
      {GetTracksDataLoader()->SetCompressionLevel(cl);}
    virtual void  SetRecParticlesComprLevel(Int_t cl) const
      {GetRecParticlesDataLoader()->SetCompressionLevel(cl);}
    
    virtual void  SetCompressionLevel(Int_t cl);//Sets compression level in all the files
    void          SetDirName(TString& name);//sets the directory name for all the I/O environment
    
    const TString& GetDetectorName() const{return fDetectorName;}//returns the name of the detector
    AliRunLoader*  GetRunLoader();//gets the run-getter from event folder
    
    void          SetDigitsFileNameSuffix(const TString& suffix) const;//adds the suffix before ".root", 
                                                          //e.g. TPC.Digits.root -> TPC.DigitsMerged.root
                                                              //made on Jiri Chudoba demand
    void Synchronize();

    /**********************************************/
    /***********     P U B L I C     **************/
    /*********       S T A T I C       ************/
    /*********         METHODS         ************/
    /*********     They are used by    ************/
    /*********** AliRunLoader as well**************/
    /**********************************************/
    static TDirectory* ChangeDir(TFile* file, Int_t eventno); //changes the root directory in "file" to directory corresponing to eventno
    static Bool_t      TestFileOption(Option_t* opt);//checks is file is created from scratch
    static Bool_t      IsOptionWritable(const TString& opt);
    
    static void        SetDebug(Int_t deb = 1);//Sets debugging information
    
   protected:

    /*********************************************/
    /************    PROTECTED      **************/
    /*********     M E T H O D S       ***********/
    /*********************************************/
    enum EDataTypes{kHits = 0,kSDigits,kDigits,kRecPoints,kTracks,kRecParticles,kNDataTypes};

    //Opens hits file and jumps to directory cooresponding to current event.
    //If dir does not exists try to create it
    //opt is passed to TFile::Open
    //reads data from the file and posts them into folder


    virtual Int_t PostHits() const
      {return GetHitsDataLoader()->GetBaseLoader(0)->Load();}
    virtual Int_t PostSDigits() const
      {return GetSDigitsDataLoader()->GetBaseLoader(0)->Load();}
    virtual Int_t PostDigits() const
      {return GetDigitsDataLoader()->GetBaseLoader(0)->Load();}
    virtual Int_t PostRecPoints() const
      {return GetRecPointsDataLoader()->GetBaseLoader(0)->Load();}
    virtual Int_t PostTracks() const
      {return GetTracksDataLoader()->GetBaseLoader(0)->Load();}
    virtual Int_t PostRecParticles() const
      {return GetRecParticlesDataLoader()->GetBaseLoader(0)->Load();}
    
    void          Clean();//calls clean for data loaders
    void          Clean(const TString& name);//removes and deletes object from data folder 
    
    
    TString       GetUnixDir() const;
    TObject*      GetDetectorData(const char* name){return GetDetectorDataFolder()->FindObject(name);}
    TObject**     GetDetectorDataRef(TObject* obj);
    
    virtual void InitDefaults();
    void ResetDataInfo();
    
    void SetTAddrInDet();//Call SetTreeAddress for corresponding detector

    /**********************************************/
    /************    PROTECTED      ***************/
    /*********        D A T A          ************/
    /**********************************************/
  
    TObjArray*    fDataLoaders; // array with data loaders each corresponds to 
                                // one data type (Hits, Digits, ...) 

    TString       fDetectorName;//detector name that this loader belongs to
    
    TFolder*      fEventFolder; //! Folder with data that changes from event to event, even in single run
    TFolder*      fDataFolder;  //!  Folder with data (hits, sdigits, etc, grouped in folders named like detectors
    TFolder*      fDetectorDataFolder;//!Folder that contains the detector data
    TFolder*      fModuleFolder;      //!Folder that contains the modules 

 // file option varible was introduced because if TFile is created with "recreate" 
 // stored option in TFile is "CREATE". We need to remeber "recreate" for
 // Max events per file functionality

    static const TString   fgkDefaultHitsContainerName;//default name of conatiner (TREE) for hits
    static const TString   fgkDefaultDigitsContainerName;//default name of conatiner (TREE) for digits
    static const TString   fgkDefaultSDigitsContainerName;//default name of conatiner (TREE) for Sdigits
    static const TString   fgkDefaultRecPointsContainerName;//default name of conatiner (TREE) for Rec Points
    static const TString   fgkDefaultTracksContainerName;//default name of conatiner (TREE) for Tracks
    static const TString   fgkDefaultRecParticlesContainerName;//default name of conatiner (TREE) for Reconstructed Particles
    static const TString   fgkLoaderBaseName;//base name of Loader: fDetectorName+fgkLoaderBaseName. i.e. PHOSLoader
    
   private:
    //descendant classes should
    //use protected interface methods to access these folders
    AliLoader(const AliLoader&); //Not implemented
    AliLoader& operator=(const AliLoader&); //Not implemented

    ClassDef(AliLoader,3)
 };
/******************************************************************/
/************************ I N L I N E S ***************************/
/******************************************************************/

inline TFolder* AliLoader::GetDetectorDataFolder()
 {
 //helper function which returns the folder of name "name" in data folder
   if(!fDetectorDataFolder)
    {
      fDetectorDataFolder = dynamic_cast<TFolder*>(GetDataFolder()->FindObject(fDetectorName.Data()));
      if(!fDetectorDataFolder)
       {
         Fatal("GetDetectorDataFolder",
               "Can not find folder %s in folder %s. Aborting",
               fDetectorName.Data(),GetDataFolder()->GetName());
         return 0x0;
       }
      
    }
   return fDetectorDataFolder;
 }

#endif
