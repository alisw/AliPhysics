#include <AliDataLoader.h>
//__________________________________________
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  class AliDataLoader                                                                    //
//                                                                                         //
//  Container of all data needed for full                                                  //
//  description of each data type                                                          //
//  (Hits, Kine, ...)                                                                      //
//                                                                                         //
//  Each data loader has a basic standard setup of BaseLoaders                             //
//  which can be identuified by indexes (defined by EStdBasicLoaders)                      //
//  Data managed by these standard base loaders has fixed naming convention                //
//  e.g. - tree with hits is always named TreeH                                            //
//                     (defined in AliLoader::fgkDefaultHitsContainerName)                 //
//       - task DtectorName+Name defined                                                   //
//                                                                                         //
//  EStdBasicLoaders   idx     Object Type        Description                              //
//      kData           0    TTree or TObject     main data itself (hits,digits,...)       //
//      kTask           1        TTask            object producing main data               //
//      kQA             2        TTree                quality assurance tree               //
//      kQATask         3        TTask            task producing QA object                 //
//                                                                                         //
//                                                                                         //
//  User can define and add more basic loaders even Run Time.                              //
//  Caution: in order to save information about added base loader                          //
//  user must rewrite Run Loader to galice.file, overwriting old setup                     //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TBits.h>

#include "AliLog.h"
#include "AliRunLoader.h"

ClassImp(AliDataLoader)

AliDataLoader::AliDataLoader():
 fFileName(0),
 fFile(0x0),
 fDirectory(0x0),
 fFileOption(),
 fCompressionLevel(2),
 fNEventsPerFile(0),
 fBaseLoaders(0x0),
 fHasTask(kFALSE),
 fTaskName(),
 fParentalTask(0x0),
 fEventFolder(0x0),
 fFolder(0x0)
{
  
}
/*****************************************************************************/ 

AliDataLoader::AliDataLoader(const char* filename, const char* contname, const char* name, Option_t* opt):
 TNamed(name,name),
 fFileName(filename),
 fFile(0x0),
 fDirectory(0x0),
 fFileOption(0),
 fCompressionLevel(2),
 fNEventsPerFile(0),
 fBaseLoaders(new TObjArray(4)),
 fHasTask(kFALSE),
 fTaskName(),
 fParentalTask(0x0),
 fEventFolder(0x0),
 fFolder(0x0)
{
//constructor
// creates a 0 loader, depending on option, default "T" is specialized loader for trees
// else standard object loader
// trees needs special care, becouse we need to set dir before writing
  AliDebug(1, Form("File name is %s",fFileName.Data()));
   
  TString option(opt);
  AliBaseLoader* bl;
  if (option.CompareTo("T",TString::kIgnoreCase) == 0)
    bl = new AliTreeLoader(contname,this);
  else 
    bl = new AliObjectLoader(contname,this);
  fBaseLoaders->AddAt(bl,kData);
  
}
/*****************************************************************************/ 
AliDataLoader::AliDataLoader(const AliDataLoader& source) : 
  TNamed(source),
  fFileName(source.fFileName),
  fFile(source.fFile),
  fDirectory(source.fDirectory),
  fFileOption(source.fFileOption),
  fCompressionLevel(source.fCompressionLevel),
  fNEventsPerFile(source.fNEventsPerFile),
  fBaseLoaders(source.fBaseLoaders),
  fHasTask(source.fHasTask),
  fTaskName(source.fTaskName),
  fParentalTask(source.fParentalTask),
  fEventFolder(source.fEventFolder),
  fFolder(source.fFolder)
{
  // copy constructor
  AliFatal("Copy constructor not implemented");
}
/*****************************************************************************/ 
AliDataLoader& AliDataLoader::operator=(const AliDataLoader& /*source*/) {
  // Assignment operator
  AliFatal("Assignment operator not implemented");
  return *this;
}
/*****************************************************************************/ 

AliDataLoader::~AliDataLoader()
{
//dtor
 UnloadAll();
}
/*****************************************************************************/ 

Int_t  AliDataLoader::SetEvent()
{
//basically the same that GetEvent but do not post data to folders
 AliRunLoader* rl = GetRunLoader();
 if (rl == 0x0)
  {
    AliError("Can not get RunGettr");
    return 1;
  }
 
 Int_t evno = rl->GetEventNumber();

 TIter next(fBaseLoaders);
 AliBaseLoader* bl;
 while ((bl = (AliBaseLoader*)next()))
  {
    if (bl->DoNotReload() == kFALSE) bl->Clean();
  }

 if(fFile)
  {
    if (CheckReload())
     {
       delete fFile;
       fFile = 0x0;
       AliDebug(1, Form("Reloading new file. File opt is %s",fFileOption.Data()));
       OpenFile(fFileOption);
     }

    fDirectory = AliLoader::ChangeDir(fFile,evno);
    if (fDirectory == 0x0)
      {
        AliError(Form("Can not chage directory in file %s",fFile->GetName()));
        return 1;
      }
   }
 return 0;
}
/*****************************************************************************/ 

Int_t  AliDataLoader::GetEvent()
{
 // posts all loaded data from files to White Board
 // event number is defined in RunLoader
 // 
 //returns:
 //     0  - in case of no error
 //     1  - event not found
 //     
 //for each base laoder post, if was loaded before GetEvent
 
 //call set event to switch to new directory in file


 //post all data that were loaded before 
 // ->SetEvent does not call Unload, but only cleans White Board
 // such IsLoaded flag stays untached
 
 if ( AliLoader::TestFileOption(fFileOption) == kTRUE ) //if file is read or update mode try to post
  {                                                     //in other case there is no sense to post: file is new
   TIter nextbl(fBaseLoaders);
   AliBaseLoader* bl;
   while ((bl = (AliBaseLoader*)nextbl()))
    {
     if (bl->IsLoaded())
      {
        if (bl->DoNotReload() == kFALSE) bl->Post();
      }
    } 
  }
 return 0;
}
/*****************************************************************************/ 

Int_t AliDataLoader::OpenFile(Option_t* opt)
{
//Opens file named 'filename', and assigns pointer to it to 'file'
//jumps to fDirectoryectory corresponding to current event and stores the pointer to it in 'fDirectory'
//option 'opt' is passed to TFile::Open
  if (fFile)
   {
     if(fFile->IsOpen() == kTRUE)
       {
         AliWarning(Form(" File %s already opened. First close it.",fFile->GetName()));
         return 0;
       }
     else
       {
         AliWarning(Form("Pointer to file %s is not null, but file is not opened",
			 fFile->GetName()));
         delete fFile;
         fFile = 0x0;
       }
   }
  
  TString fname(SetFileOffset(fFileName));
  
  fFile = (TFile *)(gROOT->GetListOfFiles()->FindObject(fname));
  if (fFile)
   {
     if(fFile->IsOpen() == kTRUE)
       {
         AliWarning(Form("File %s already opened by sombody else. First close it.",
			 fFile->GetName()));
         return 0;
       }
   }
  
  fFileOption = opt;
  fFile = TFile::Open(fname,fFileOption);//open the file
  if (fFile == 0x0)
   {//file is null
     AliError(Form("Can not open file %s",fname.Data()));
     return 1;
   }
  if (fFile->IsOpen() == kFALSE)
   {//file is null
     AliError(Form("Can not open file %s",fname.Data()));
     return 1;
   }

  fFile->SetCompressionLevel(fCompressionLevel);
  
  AliRunLoader* rg = GetRunLoader();
  if (rg == 0x0)
   {
     AliError("Can not find Run-Loader in folder.");
     return 2;
   }
  Int_t evno = rg->GetEventNumber();
  
  fDirectory = AliLoader::ChangeDir(fFile,evno);
  if (fDirectory == 0x0)
   {
     AliError(Form("Can not chage fDirectory in file %s.",fFile->GetName()));
     return 3; 
   }
  return 0;
}
/*****************************************************************************/ 

void AliDataLoader::Unload()
{
 //unloads main data -  shortcut method 
  GetBaseLoader(0)->Unload();
}
/*****************************************************************************/ 

void AliDataLoader::UnloadAll()
{
//Unloads all data and tasks
 if ( fFile == 0x0 ) return; //nothing loaded
 
 TIter next(fBaseLoaders);
 AliBaseLoader* bl;
 while ((bl = (AliBaseLoader*)next()))
  {
    bl->Unload();
  }
}
/*****************************************************************************/ 

Int_t AliDataLoader::Reload()
{
 //Unloads and loads data again
 if ( fFile == 0x0 ) return 0;
   
 TBits loaded(fBaseLoaders->GetEntries());  
 TIter next(fBaseLoaders);
 AliBaseLoader* bl;

 Int_t i = 0;
 while ((bl = (AliBaseLoader*)next()))
  {
    if (bl->IsLoaded())
     {
       loaded.SetBitNumber(i++,kTRUE);
       bl->Unload();
     }
  }
 
 Int_t retval;
 i = 0;  
 next.Reset();
 while ((bl = (AliBaseLoader*)next()))
  {
    if (loaded.TestBitNumber(i++))
     {
       retval = bl->Load(fFileOption);
       if (retval) 
        {
         AliError(Form("Error occur while loading %s",bl->GetName()));
         return retval;
        }
     }
  }
 

 return 0;
 }
/*****************************************************************************/ 
Int_t AliDataLoader::WriteData(Option_t* opt)
{
//Writes primary data ==  first BaseLoader
  AliDebug(1, Form("Writing %s container for %s data. Option is %s.",
		   GetBaseLoader(0)->GetName(),GetName(),opt));
  return GetBaseLoader(0)->WriteData(opt);
}
/*****************************************************************************/ 

Int_t AliDataLoader::Load(Option_t* opt)
{
//Writes primary data ==  first BaseLoader
  return GetBaseLoader(0)->Load(opt);
}
/*****************************************************************************/ 

Int_t  AliDataLoader::SetEventFolder(TFolder* eventfolder)
{
 //sets the event folder
 if (eventfolder == 0x0)
  {
    AliError("Stupid joke. Argument is NULL");
    return 1;
  }
 AliDebug(1, Form("name = %s Setting Event Folder named %s.",
		  GetName(),eventfolder->GetName()));

 fEventFolder = eventfolder;
 return 0;
}
/*****************************************************************************/ 

Int_t  AliDataLoader::SetFolder(TFolder* folder)
{
  // Sets the folder and the data loaders
 if (folder == 0x0)
  {
    AliError("Stupid joke. Argument is NULL");
    return 1;
  }
 
 AliDebug(1, Form("name = %s Setting folder named %s.",GetName(),folder->GetName()));
 
 fFolder = folder;
 TIter next(fBaseLoaders);
 AliBaseLoader* bl;

 while ((bl = (AliBaseLoader*)next()))
  {
   bl->SetDataLoader(this);
  }  

 return 0;
}
/******************************************************************/

TFolder* AliDataLoader::GetEventFolder()
{
//get EVENT folder (data that are changing from event to event, even in single run)
  AliDebug(1, "EF = %#x");
  return fEventFolder;
}
/*****************************************************************************/ 

AliRunLoader* AliDataLoader::GetRunLoader()
{
//gets the run-loader from event folder
  AliRunLoader* rg = 0x0;
  TFolder* ef = GetEventFolder();
  if (ef == 0x0)
   {
     AliError("Can not get event folder.");
     return 0;
   }
  rg = dynamic_cast<AliRunLoader*>(ef->FindObject(AliRunLoader::GetRunLoaderName()));
  return rg;
}

/*****************************************************************************/ 
void AliDataLoader::CloseFile()
{
  //closes file
  TIter next(fBaseLoaders);
  AliBaseLoader* bl;
  while ((bl = (AliBaseLoader*)next()))
   {
     if (bl->IsLoaded()) return;
   }
  
  AliDebug(1, "Closing and deleting (object) file.");
    
  delete fFile;
  fFile = 0x0;
  fDirectory = 0x0;
}
/*****************************************************************************/ 

void AliDataLoader::Clean()
{
  //Cleans main data
  GetBaseLoader(0)->Clean();
}  
/*****************************************************************************/ 

void AliDataLoader::CleanAll()
{
  //Cleans all folders and tasks
  TIter next(fBaseLoaders);
  AliBaseLoader* bl;
  while ((bl = (AliBaseLoader*)next()))
   {
      bl->Clean();
   }
}
/*****************************************************************************/ 

void AliDataLoader::SetFileNameSuffix(const TString& suffix)
{
  //adds the suffix before ".root", 
  //e.g. TPC.Digits.root -> TPC.DigitsMerged.root
  //made on Jiri Chudoba demand
  AliDebug(1, Form("suffix=%s",suffix.Data()));
  AliDebug(1, Form("   Digits File Name before: %s",fFileName.Data()));
   
  static TString dotroot(".root");
  const TString& suffixdotroot = suffix + dotroot;
  fFileName = fFileName.ReplaceAll(dotroot,suffixdotroot);

  AliDebug(1, Form("                    after : %s",fFileName.Data()));
}
/*****************************************************************************/ 

Bool_t AliDataLoader::CheckReload()
{
//checks if we have to reload given file
 if (fFile == 0x0) return kFALSE;
 TString tmp = SetFileOffset(fFileName);
 if (tmp.CompareTo(fFile->GetName())) return kTRUE;  //file must be reloaded
 return  kFALSE;
}
/*****************************************************************************/ 

const TString AliDataLoader::SetFileOffset(const TString& fname)
{

//return fname;
  Long_t offset = (Long_t)GetRunLoader()->GetFileOffset();
  if (fNEventsPerFile > 0) {
    offset = GetRunLoader()->GetEventNumber()/fNEventsPerFile;
  }
  if (offset < 1) return fname;

  TString soffset;
  soffset += offset;//automatic conversion to string
  TString dotroot(".root");
  const TString& offfsetdotroot = offset + dotroot;
  TString out = fname;
  out = out.ReplaceAll(dotroot,offfsetdotroot);
  AliDebug(1, Form("in=%s  out=%s.",fname.Data(),out.Data()));
  return out;

}
/*****************************************************************************/ 

void AliDataLoader::SetFileOption(Option_t* newopt)
{
  //sets file option
  if (fFileOption.CompareTo(newopt) == 0) return;
  fFileOption = newopt;
  Reload();
}
/*****************************************************************************/ 

void AliDataLoader::SetCompressionLevel(Int_t cl)
{
//sets comression level for data defined by di
  fCompressionLevel = cl;
  if (fFile) fFile->SetCompressionLevel(cl);
}
/*****************************************************************************/ 

void AliDataLoader::MakeTree()
{
  // Makes tree for the current data loader
  AliTreeLoader* tl = dynamic_cast<AliTreeLoader*>(fBaseLoaders->At(0));
  if (tl == 0x0)
   {
     AliError("Can not make a tree because main base loader is not a tree loader");
     return;
   }
  tl->MakeTree();
}
/*****************************************************************************/ 

Bool_t AliDataLoader::IsFileWritable() const
{
//returns true if file is writable
 return (fFile)?fFile->IsWritable():kFALSE;
}
/*****************************************************************************/ 

Bool_t AliDataLoader::IsFileOpen() const
{
//returns true if file is writable
 return (fFile)?fFile->IsOpen():kFALSE;
}
/*****************************************************************************/ 

Bool_t AliDataLoader::IsOptionContrary(const TString& option) const
{
//Checks if passed option is contrary with file open option 
//which is passed option "writable" and existing option not wriable
//in reverse case it is no harm so it is NOT contrary
  if (fFile == 0x0) return kFALSE; //file is not opened - no problem
  
  if ( ( AliLoader::IsOptionWritable(option)      == kTRUE  ) &&     // passed option is writable and 
       ( AliLoader::IsOptionWritable(fFileOption) == kFALSE )    )   // existing one is not writable
    {
      return kTRUE;
    }

  return kFALSE;
}
/*****************************************************************************/ 
void AliDataLoader::AddBaseLoader(AliBaseLoader* bl)
{
//Adds a base loader to lits of base loaders managed by this data loader
//Managed data/task will be stored in proper root directory,
//and posted to 
// - in case of tree/object - data folder connected with detector associated with this data loader
// - in case of task - parental task which defined in this AliTaskLoader 

 if (bl == 0x0)
  {
    AliWarning("Pointer is null.");
    return;
  }
 
 TObject* obj = fBaseLoaders->FindObject(bl->GetName());
 if (obj)
  {
    AliError("Can not add this base loader.");
    AliError(Form("There exists already base loader which manages data named %s for this detector.",obj->GetName()));
    return;
  }
 
 
 fBaseLoaders->Add(bl);
}

/*****************************************************************************/ 

AliBaseLoader* AliDataLoader::GetBaseLoader(const TString& name) const
{
  return dynamic_cast<AliBaseLoader*>(fBaseLoaders->FindObject(name));
}
/*****************************************************************************/ 

AliBaseLoader* AliDataLoader::GetBaseLoader(Int_t n) const
{
  // Gets the n-th base loader (what is n?)
 return dynamic_cast<AliBaseLoader*>(fBaseLoaders->At(n));
}
/*****************************************************************************/ 

TTree* AliDataLoader::Tree() const
{
//returns tree from the main base loader
//it is just shortcut method for comfort of user
//main storage object does not have to be Tree  - 
//that is why first we need to check if it is a TreeLoader 
 AliTreeLoader* tl = dynamic_cast<AliTreeLoader*>(GetBaseLoader(0));
 if (tl == 0x0) return 0x0;
 return tl->Tree();
}
/*****************************************************************************/ 

void  AliDataLoader::SetDirName(TString& dirname)
{
  // Sets the directory name where the files will be stored
  AliDebug(10, Form("FileName before %s",fFileName.Data()));

  Int_t n = fFileName.Last('/');

  AliDebug(10, Form("Slash found on pos %d",n));

  if (n > 0) fFileName = fFileName.Remove(0,n+1);

  AliDebug(10, Form("Core FileName %s",fFileName.Data()));

  fFileName = dirname + "/" + fFileName;

  AliDebug(10, Form("FileName after %s",fFileName.Data()));
}
/*****************************************************************************/ 
AliObjectLoader* AliDataLoader::GetBaseDataLoader()
{
  // Gets the base data loader
 return dynamic_cast<AliObjectLoader*>(GetBaseLoader(kData));
}
/*****************************************************************************/ 
AliTaskLoader* AliDataLoader::GetBaseTaskLoader()
{
  // Gets the base task loader
 return dynamic_cast<AliTaskLoader*>(GetBaseLoader(kTask));
}
/*****************************************************************************/ 
AliBaseLoader* AliDataLoader::GetBaseQALoader()
{
  // Gets the base QA loader
  return GetBaseLoader(kQA);
}
/*****************************************************************************/ 
AliTaskLoader* AliDataLoader::GetBaseQATaskLoader()
{
//returns pointer to QA base loader
 return dynamic_cast<AliTaskLoader*>(GetBaseLoader(kQATask));
}
/*****************************************************************************/ 
void AliDataLoader::SetBaseDataLoader(AliBaseLoader* bl)
{
//sets data base loader
  if (bl == 0x0)
   {
     AliError("Parameter is null");
     return;
   }
  if (GetBaseDataLoader()) delete GetBaseDataLoader();
  fBaseLoaders->AddAt(bl,kData);
}
/*****************************************************************************/ 
void AliDataLoader::SetBaseTaskLoader(AliTaskLoader* bl)
{
//sets Task base loader
  if (bl == 0x0)
   {
     AliError("Parameter is null");
     return;
   }
  if (GetBaseTaskLoader()) delete GetBaseTaskLoader();
  fBaseLoaders->AddAt(bl,kTask);
}
/*****************************************************************************/ 
void AliDataLoader::SetBaseQALoader(AliBaseLoader* bl)
{
//sets QA base loader
  if (bl == 0x0)
   {
     AliError("Parameter is null");
     return;
   }
  if (GetBaseQALoader()) delete GetBaseQALoader();
  fBaseLoaders->AddAt(bl,kQA);
}
/*****************************************************************************/ 
void AliDataLoader::SetBaseQATaskLoader(AliTaskLoader* bl)
{
//sets QA Task base loader
  if (bl == 0x0)
   {
     AliError("Parameter is null");
     return;
   }
  if (GetBaseQATaskLoader()) delete GetBaseQATaskLoader();
  fBaseLoaders->AddAt(bl,kQATask);
}
void AliDataLoader::Synchronize()
{
  //synchrinizes all writtable files 
  if ( fFile ) fFile->Flush();
}

/*****************************************************************************/ 
/*****************************************************************************/ 
/*****************************************************************************/ 
//__________________________________________
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class AliBaseLoader                                                      //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
ClassImp(AliBaseLoader)

AliBaseLoader::AliBaseLoader():
 fIsLoaded(kFALSE),
 fStoreInTopOfFile(kFALSE),
 fDoNotReload(kFALSE),
 fDataLoader(0x0)
{
  //default constructor
}
/*****************************************************************************/ 

AliBaseLoader::AliBaseLoader(const TString& name,  AliDataLoader* dl, Bool_t storeontop):
 TNamed(name,name+" Base Loader"),
 fIsLoaded(kFALSE),
 fStoreInTopOfFile(storeontop),
 fDoNotReload(storeontop),//if stored on top of file - this object is loaded ones pe
 fDataLoader(dl)
{
  //constructor
}

/*****************************************************************************/ 
AliBaseLoader::AliBaseLoader(const AliBaseLoader& source) : 
  TNamed(source),
  fIsLoaded(source.fIsLoaded),
  fStoreInTopOfFile(source.fStoreInTopOfFile),
  fDoNotReload(source.fDoNotReload),
  fDataLoader(source.fDataLoader)
{
  // copy constructor
  AliFatal("Copy constructor not implemented");
}
/*****************************************************************************/ 
AliBaseLoader& AliBaseLoader::operator=(const AliBaseLoader& /*source*/) {
  // Assignment operator
  AliFatal("Assignment operator not implemented");
  return *this;
}
/*****************************************************************************/ 

Int_t AliBaseLoader::Load(Option_t* opt)
{
  // Loads and posts the data
  AliDebug(1, Form("data type = %s, option = %s",GetName(),opt));

  if (Get())
   {
     AliWarning(Form("Data <<%s>> are already loaded. Use ReloadData to force reload. Nothing done",GetName()));
     return 0;
   }
  
  Int_t retval;
  
  if (GetDataLoader()->IsFileOpen() == kTRUE)
   {
     if (GetDataLoader()->IsOptionContrary(opt) == kTRUE)
       {
         AliError(Form("Data Type %s, Container Name %s", GetDataLoader()->GetName(),GetName()));
         AliError("File was already opened in READ-ONLY mode, while now WRITEABLE access is requested.");
         AliError("Use AliDataLoader::SetOption to enforce change of access mode OR");
         AliError("Load previosly loaded data with coherent option.");
         return 10;
       }
   }
  else
   {
     retval = GetDataLoader()->OpenFile(opt);
     if (retval) 
      {
        AliError(Form("Error occured while opening <<%s>> file",GetName()));
        return retval;
      }
   }
  //if file is recreated there is no sense to search for data to post and get Error message
  if (AliLoader::TestFileOption(opt) == kFALSE)
   {
    AliTreeLoader* tl = dynamic_cast<AliTreeLoader*>(this);
    if (tl) tl->MakeTree();
    fIsLoaded = kTRUE;
    return 0;
   }

  retval = Post();
  if (retval)
   {
    AliError(Form("Error occured while posting %s from file to folder.",GetName()));
    return retval;
   }
  
  fIsLoaded = kTRUE;
  return 0;
}
/*****************************************************************************/ 

Int_t AliBaseLoader::Post()
{
//Posts data container to proper folders

  if ( GetDirectory() == 0x0)
   {
     AliError(Form("%s directory is NULL. Load before.",GetDataLoader()->GetName()));
     return 2; 
   }
  
  TObject* data = GetFromDirectory(fName);
  if(data)
   {
     //if such an obejct already exists - remove it first
     return Post(data);
   }
  else
   {
    //check if file is in update mode
    Int_t fileupdate = GetDataLoader()->GetFileOption().CompareTo("update",TString::kIgnoreCase);
    if ( fileupdate == 0)
     { //if it is, it is normal that there is no data yet
       AliDebug(1, Form("Can not find %s in file %s (file is opened in UPDATE mode).",
			GetName(),GetDataLoader()->GetFile()->GetName()));
     }
    else
     {
        AliError(Form("Can not find %s in file %s", GetName(),GetDataLoader()->GetFile()->GetName()));
        return 5;
     }
   }
  return 0;
}
/*****************************************************************************/ 

Int_t AliBaseLoader::Post(TObject* data)
{
//Posts data container to proper folders
 if (data == 0x0)
  {
    AliError("Pointer to object is NULL");
    return 1;
  }

 if ( fName.CompareTo(data->GetName()) != 0)
   {
     AliFatal(Form("Object name is <<%s>> while <<%s>> expected",data->GetName(),GetName()));
     return -1;//pro forma
   }
   
 TObject* obj = Get();
 if (data == obj)
  {
    AliWarning("This object was already posted.");
    return 0;
  }
 if (obj)
  {
    AliWarning(Form("Object named %s already exitsts in data folder. Removing it",GetName()));
    Clean();
  }
 return AddToBoard(data);
}
/*****************************************************************************/ 

Int_t AliBaseLoader::WriteData(Option_t* opt)
{
//Writes data defined by di object
//opt might be "OVERWRITE" in case of forcing overwriting
  AliDebug(1, Form("Writing %s container for %s data. Option is %s.",
		   GetName(),GetDataLoader()->GetName(),opt));
  
  TObject *data = Get();
  if(data == 0x0)
   {//did not get, nothing to write
     AliWarning(Form("Tree named %s not found in folder. Nothing to write.",GetName()));
     return 0;
   }
  
  //check if file is opened
  if (GetDirectory() == 0x0)
   { 
     //if not try to open
     GetDataLoader()->SetFileOption("UPDATE");
     if (GetDataLoader()->OpenFile("UPDATE"))
      {  
        //oops, can not open the file, give an error message and return error code
        AliError(Form("Can not open hits file. %s ARE NOT WRITTEN",GetName()));
        return 1;
      }
   }

  if (GetDataLoader()->IsFileWritable() == kFALSE)
   {
     AliError(Form("File %s is not writable",GetDataLoader()->GetFileName().Data()));
     return 2;
   }
  
  GetDirectory()->cd(); //set the proper directory active

  //see if hits container already exists in this (root) directory
  TObject* obj = GetFromDirectory(GetName());
  if (obj)
   { //if they exist, see if option OVERWRITE is used
     const char *oOverWrite = strstr(opt,"OVERWRITE");
     if(!oOverWrite)
      {//if it is not used -  give an error message and return an error code
        AliError("Tree already exisists. Use option \"OVERWRITE\" to overwrite previous data");
        return 3;
      }
   }
  
  AliDebug(1, Form("DataName = %s, opt = %s, data object name = %s",
                   GetName(),opt,data->GetName()));
  AliDebug(1, Form("File Name = %s, Directory Name = %s Directory's File Name = %s",
                   GetDataLoader()->GetFile()->GetName(),GetDirectory()->GetName(),
                   GetDirectory()->GetFile()->GetName()));
  
  AliDebug(1, "Writing data");
  data->Write(0,TObject::kOverwrite);

  fIsLoaded = kTRUE;  // Just to ensure flag is on. Object can be posted manually to folder structure, not using loader.

  return 0;
 
}
/*****************************************************************************/ 

Int_t AliBaseLoader::Reload()
{
//Unloads and loads datat again - if loaded before
 if (IsLoaded())
  {
    Unload();
    return Load(GetDataLoader()->GetFileOption());
  }
  return 0;
}
/*****************************************************************************/ 

void AliBaseLoader::Clean()
{
//removes objects from folder/task
  AliDebug(1, Form("%s %s",GetName(),GetDataLoader()->GetName()));
  TObject* obj = Get();
  if(obj)
   { 
     AliDebug(1, Form("cleaning %s.",GetName()));
     RemoveFromBoard(obj);
     delete obj;
   }
}
/*****************************************************************************/ 

void AliBaseLoader::Unload()
{
  // Unloads data and closes the files
  Clean();
  fIsLoaded = kFALSE;
  GetDataLoader()->CloseFile();
}
/*****************************************************************************/ 
AliDataLoader* AliBaseLoader::GetDataLoader() const
{
  // Returns pointer to the data loader
 if (fDataLoader == 0x0) 
  {
    AliFatal("Pointer to Data Loader is NULL");
  }
 return fDataLoader;
}
/*****************************************************************************/ 

TDirectory* AliBaseLoader::GetDirectory() const
{
 // returnd TDirectory where data are to be saved
 //if fStoreInTopOfFile flag is true - returns pointer to file
  return (fStoreInTopOfFile)?GetDataLoader()->GetFile():GetDataLoader()->GetDirectory();
}
/*****************************************************************************/ 
/*****************************************************************************/ 
/*****************************************************************************/ 
//__________________________________________
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class AliObjectLoader                                                      //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ClassImp(AliObjectLoader)

AliObjectLoader::AliObjectLoader(const TString& name, AliDataLoader* dl, Bool_t storeontop):
 AliBaseLoader(name,dl,storeontop)
{
//constructor
}
/*****************************************************************************/ 

TFolder* AliObjectLoader::GetFolder() const
{
  // Returns pointer to the object folder
  TFolder* df = GetDataLoader()->GetFolder();
  if (df == 0x0)
   {
     AliFatal("Data Folder is NULL");
   }
  return df;
}
/*****************************************************************************/ 
AliObjectLoader::AliObjectLoader(const AliObjectLoader& source):
  AliBaseLoader(source) {
  // copy constructor
  AliFatal("Copy constructor not implemented");
}
/*****************************************************************************/ 
AliObjectLoader& AliObjectLoader::operator=(const AliObjectLoader& /*source*/) {
  // Assignment operator
  AliFatal("Assignment operator not implemented");
  return *this;
}
/*****************************************************************************/ 

void AliObjectLoader::RemoveFromBoard(TObject* obj)
{
  // Removes "obj" from the board
  GetFolder()->Remove(obj);
}
/*****************************************************************************/ 
Int_t AliObjectLoader::AddToBoard(TObject* obj)
{
  // Adds "obj" to the board
  GetFolder()->Add(obj);
  return 0;
}
/*****************************************************************************/ 

TObject* AliObjectLoader::Get() const
{
  // Returns pointer to the object loader
  return (GetFolder()) ? GetFolder()->FindObject(GetName()) : 0x0;
}

/*****************************************************************************/ 
/*****************************************************************************/ 
/*****************************************************************************/ 
//__________________________________________
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class AliTreeLoader                                                      //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ClassImp(AliTreeLoader)

AliTreeLoader::AliTreeLoader(const TString& name, AliDataLoader* dl,Bool_t storeontop):
 AliObjectLoader(name,dl,storeontop)
{
//constructor
}
/*****************************************************************************/ 
AliTreeLoader::AliTreeLoader(const AliTreeLoader& source):
  AliObjectLoader(source) {
  // copy constructor
  AliFatal("Copy constructor not implemented");
}
/*****************************************************************************/ 
AliTreeLoader& AliTreeLoader::operator=(const AliTreeLoader& /*source*/) {
  // Assignment operator
  AliFatal("Assignment operator not implemented");
  return *this;
}

/*****************************************************************************/ 

Int_t AliTreeLoader::WriteData(Option_t* opt)
{
//Writes data defined by di object
//opt might be "OVERWRITE" in case of forcing overwriting

  AliDebug(1, Form("Writing %s container for %s data. Option is %s.",
		   GetName(),GetDataLoader()->GetName(),opt));

  TObject *data = Get();
  if(data == 0x0)
   {//did not get, nothing to write
     AliWarning(Form("Tree named %s not found in folder. Nothing to write.",GetName()));
     return 0;
   }
  
  //check if file is opened
  if (GetDirectory() == 0x0)
   { 
     //if not try to open
     GetDataLoader()->SetFileOption("UPDATE");
     if (GetDataLoader()->OpenFile("UPDATE"))
      {  
        //oops, can not open the file, give an error message and return error code
        AliError(Form("Can not open hits file. %s ARE NOT WRITTEN",GetName()));
        return 1;
      }
   }

  if (GetDataLoader()->IsFileWritable() == kFALSE)
   {
     AliError(Form("File %s is not writable",GetDataLoader()->GetFileName().Data()));
     return 2;
   }
  
  GetDirectory()->cd(); //set the proper directory active

  //see if hits container already exists in this (root) directory
  TObject* obj = GetFromDirectory(GetName());
  if (obj)
   { //if they exist, see if option OVERWRITE is used
     const char *oOverWrite = strstr(opt,"OVERWRITE");
     if(!oOverWrite)
      {//if it is not used -  give an error message and return an error code
        AliError("Tree already exisists. Use option \"OVERWRITE\" to overwrite previous data");
        return 3;
      }
   }
  
  AliDebug(1, Form("DataName = %s, opt = %s, data object name = %s",
                   GetName(),opt,data->GetName()));
  AliDebug(1, Form("File Name = %s, Directory Name = %s Directory's File Name = %s",
                   GetDataLoader()->GetFile()->GetName(),GetDirectory()->GetName(),
                   GetDirectory()->GetFile()->GetName()));
  
  //if a data object is a tree set the directory
  TTree* tree = dynamic_cast<TTree*>(data);
  if (tree) tree->SetDirectory(GetDirectory()); //forces setting the directory to this directory (we changed dir few lines above)
  
  AliDebug(1, "Writing tree");
  data->Write(0,TObject::kOverwrite);

  fIsLoaded = kTRUE;  // Just to ensure flag is on. Object can be posted manually to folder structure, not using loader.
  
  return 0;
 
}
/*****************************************************************************/ 

void AliTreeLoader::MakeTree()
{
//this virtual method creates the tree in the file
  if (Tree()) 
   {
    AliDebug(1, Form("name = %s, Data Name = %s Tree already exists.",
		     GetName(),GetDataLoader()->GetName()));
    return;//tree already made 
   }
  AliDebug(1, Form("Making Tree named %s.",GetName()));
   
  TString dtypename(GetDataLoader()->GetName());
  TTree* tree = new TTree(GetName(), dtypename + " Container"); //make a tree
  if (tree == 0x0)
   {
     AliError(Form("Can not create %s tree.",GetName()));
     return;
   }
  tree->SetAutoSave(1000000000); //no autosave
  GetFolder()->Add(tree);
  WriteData("OVERWRITE");//write tree to the file
}


/*****************************************************************************/ 
/*****************************************************************************/ 
/*****************************************************************************/ 
//__________________________________________
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class AliTaskLoader                                                      //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ClassImp(AliTaskLoader)

AliTaskLoader::AliTaskLoader(const TString& name, AliDataLoader* dl, TTask* parentaltask, Bool_t storeontop):
 AliBaseLoader(name,dl,storeontop),
 fParentalTask(parentaltask)
{
//constructor
}

/*****************************************************************************/ 
AliTaskLoader::AliTaskLoader(const AliTaskLoader& source):
  AliBaseLoader(source),
  fParentalTask(source.fParentalTask)
{
  // copy constructor
  AliFatal("Copy constructor not implemented");
}
/*****************************************************************************/ 
AliTaskLoader& AliTaskLoader::operator=(const AliTaskLoader& /*source*/) {
  // Assignment operator
  AliFatal("Assignment operator not implemented");
  return *this;
}
/*****************************************************************************/ 
void AliTaskLoader::Clean()
{
//removes tasl from parental task
// DO NOT DELETE OBJECT contrary to BaseLoader
//
  AliDebug(1, Form("Clean","%s %s",GetName(),GetDataLoader()->GetName()));
  TObject* obj = Get();
  if(obj)
   { 
     AliDebug(1, Form("cleaning %s.",GetName()));
     RemoveFromBoard(obj);
   }
}
/*****************************************************************************/ 

void AliTaskLoader::RemoveFromBoard(TObject* obj)
{
  // Removes the task "obj" from the board
  GetParentalTask()->GetListOfTasks()->Remove(obj);
}
/*****************************************************************************/ 

Int_t AliTaskLoader::AddToBoard(TObject* obj)
{
  // Adds task "obj" to the board
  TTask* task = dynamic_cast<TTask*>(obj);
  if (task == 0x0)
   {
     AliError("To TTask board can be added only tasks.");
     return 1;
   }
  GetParentalTask()->Add(task);
  return 0;
}
/*****************************************************************************/ 

TObject* AliTaskLoader::Get() const
{
  // Returns pointer to the current task
  return (GetParentalTask()) ? GetParentalTask()->GetListOfTasks()->FindObject(GetName()) : 0x0;
}
/*****************************************************************************/ 

TTask* AliTaskLoader::GetParentalTask() const
{
//returns parental tasks for this task
  return fParentalTask;
}

/*****************************************************************************/ 

/*****************************************************************************/ 
/*****************************************************************************/ 
/*****************************************************************************/ 


