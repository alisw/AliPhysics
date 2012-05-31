//
//Class (base) responsible for management of data:
//    - opening and closing proper files
//    - posting data to folders
//    - writing data from folders to proper files
//
//Author: Alice Offline Group http://alisoft.cern.ch
//Responsible: Piotr.Skowronski@cern.ch
//
#include "AliLoader.h"

//Root includes
#include <TROOT.h>
#include <TFile.h>
#include <TFolder.h>
#include <TString.h>

//AliRoot includes
#include "AliConfig.h"
#include "AliDetector.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"

const TString AliLoader::fgkDefaultHitsContainerName("TreeH");
const TString AliLoader::fgkDefaultDigitsContainerName = "TreeD";
const TString AliLoader::fgkDefaultSDigitsContainerName = "TreeS";
const TString AliLoader::fgkDefaultRecPointsContainerName = "TreeR";
const TString AliLoader::fgkDefaultTracksContainerName = "TreeT";
const TString AliLoader::fgkDefaultRecParticlesContainerName = "TreeP";
const TString AliLoader::fgkLoaderBaseName("Loader");

ClassImp(AliLoader)
//___________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliLoader                                                //
//                                                                 //
//  Base class for Loaders.                                        //
//  Loader provides base I/O fascilities for "standard" data.      //
//  Each detector has a laoder data member                         //
//  loader is accessible via folder structure as well              //
//                                                                 //
/////////////////////////////////////////////////////////////////////
 
/*****************************************************************************/ 

AliLoader::AliLoader():
 fDataLoaders(0x0),
 fDetectorName(""),
 fEventFolder(0x0),
 fDataFolder(0x0),
 fDetectorDataFolder(0x0),
 fModuleFolder(0x0)
{
  //default constructor
  
}
/******************************************************************/

AliLoader::AliLoader(const Char_t* detname,const Char_t* eventfoldername):
 fDataLoaders(new TObjArray(kNDataTypes)),
 fDetectorName(""),
 fEventFolder(0x0),
 fDataFolder(0x0),
 fDetectorDataFolder(0x0),
 fModuleFolder(0x0)
{
  //ctor
   AliDebug(1, Form("detname = %s eventfoldername = %s",detname,eventfoldername));

   //try to find folder eventfoldername in top alice folder
   //safe because GetTopFolder will abort in case of failure

   fDetectorName = detname;
   fName = fDetectorName+"Loader";
   InitDefaults();

   TObject* fobj = GetTopFolder()->FindObject(eventfoldername);
   fEventFolder = (fobj)?dynamic_cast<TFolder*>(fobj):0x0;//in case FindObject returns NULL dynamic cast cause seg. fault
   SetEventFolder(fEventFolder);
   
 }
/*****************************************************************************/ 

AliLoader::AliLoader(const Char_t * detname,TFolder* eventfolder):
 fDataLoaders(new TObjArray(kNDataTypes)),
 fDetectorName(detname),
 fEventFolder(0x0),
 fDataFolder(0x0),
 fDetectorDataFolder(0x0),
 fModuleFolder(0x0)
{
//constructor
   fDetectorName = detname;
   fName = fDetectorName+"Loader";
   InitDefaults();
   SetEventFolder(eventfolder);
   //fileoption's don't need to initialized because default TString ctor does it correctly
}

/*****************************************************************************/ 
AliLoader::~AliLoader()
{
//detructor
  if (fDataLoaders) fDataLoaders->SetOwner();
  delete fDataLoaders;
}
/*****************************************************************************/ 

void AliLoader::InitDefaults()
{
  // H I T S 
  AliDataLoader* dl;
  dl = new AliDataLoader(fDetectorName + ".Hits.root",fgkDefaultHitsContainerName, "Hits" );
  fDataLoaders->AddAt(dl,kHits);
  
  
  // S U M M A B L E   D I G I T S
  dl = new AliDataLoader(fDetectorName + ".SDigits.root",fgkDefaultSDigitsContainerName, "Summable Digits");
  fDataLoaders->AddAt(dl,kSDigits);

  // D I G I T S  
  dl = new AliDataLoader(fDetectorName + ".Digits.root",fgkDefaultDigitsContainerName, "Digits");
  fDataLoaders->AddAt(dl,kDigits);
  
  // R E C O N S T R U C T E D   P O I N T S aka C L U S T E R S 
  dl = new AliDataLoader(fDetectorName + ".RecPoints.root",fgkDefaultRecPointsContainerName, "Reconstructed Points");
  fDataLoaders->AddAt(dl,kRecPoints);
  
  // T R A C K S
  dl = new AliDataLoader(fDetectorName + ".Tracks.root",fgkDefaultTracksContainerName, "Tracks");
  fDataLoaders->AddAt(dl,kTracks);
  
  // R E C O N S T R U C T E D   P O I N T S aka C L U S T E R S 
  dl = new AliDataLoader(fDetectorName + ".RecParticles.root",fgkDefaultRecParticlesContainerName, "Reconstructed Particles");
  fDataLoaders->AddAt(dl,kRecParticles);

 }
/*****************************************************************************/ 

AliDataLoader* AliLoader::GetDataLoader(const char* name)
{
//returns Data Loader with specified name
  return dynamic_cast<AliDataLoader*>(fDataLoaders->FindObject(name));
}
/*****************************************************************************/ 
void AliLoader::AddDataLoader(AliDataLoader* dl)
{
  //
  // Adds a data loader
  //
  if (dl == 0x0)
   {
     AliError("Pointer is NULL");
     return;
   }
  if (fDataLoaders->FindObject(dl->GetName()))
   {
     AliError("Such a loader exists");
     return;
   }
  fDataLoaders->AddLast(dl);
  dl->SetEventFolder(fEventFolder);
  dl->SetFolder(GetDetectorDataFolder()); //Must exists - ensure register is called before
}
/*****************************************************************************/ 

Int_t AliLoader::SetEvent()
{
 //basically the same that GetEvent but do not post data to folders
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 while ((dl = (AliDataLoader*)next()))
  {
    dl->SetEvent();
  }
 return 0;
}
/******************************************************************/

void AliLoader::UnloadAll()
{
 //calls UnloadAll for all base laoders
 //Unloads everything
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 while ((dl = (AliDataLoader*)next()))
  {
    dl->UnloadAll();
  }
}
/******************************************************************/

Int_t AliLoader::GetEvent()
{
 //changes to proper root  directory and tries to load data from files to folders
 // event number is defined in RunLoader
 // 
 //returns:
 //     0  - in case of no error
 //     1  - event not found
 //     
 
 Int_t retval;   
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 while ((dl = (AliDataLoader*)next()))
  {
    retval = dl->GetEvent();
    if (retval)
     {
       AliError(Form("Error occured while GetEvent for %s",dl->GetName()));
       return retval;
     }
  }

 return 0;
}

/******************************************************************/

TFolder* AliLoader::GetTopFolder()
{
//returns TOP aliroot folder, just a simple interface to AliConfig (gives shorter notation)
 return AliConfig::Instance()->GetTopFolder();
}

/******************************************************************/

TFolder* AliLoader::GetEventFolder()
{
//get EVENT folder (data that are changing from event to event, even in single run)
  return fEventFolder;
}
/******************************************************************/
TFolder* AliLoader::GetDataFolder()
{
//returns the folder speciofic to given detector e.g. /Folders/Event/Data/
 if (!fDataFolder)
  {
   fDataFolder =  dynamic_cast<TFolder*>(GetEventFolder()->FindObject(AliConfig::Instance()->GetDataFolderName()));
   
   if (!fDataFolder)
    {
     AliFatal("Can not find AliRoot data folder. Aborting");
     return 0x0;
    }
  }
  return fDataFolder;
}

/*****************************************************************************/ 

TFolder* AliLoader::GetModulesFolder()
{
  //returns pointer to the folder containing modules
 if (!fModuleFolder)
  {
   fModuleFolder =  dynamic_cast<TFolder*>(GetEventFolder()->FindObjectAny(AliConfig::GetModulesFolderName()));
   
   if (!fModuleFolder)
    {
     AliFatal("Can not find modules folder. Aborting");
     return 0x0;
    }
  }
 return fModuleFolder;
   
}

/*****************************************************************************/ 

TDirectory* AliLoader::ChangeDir(TFile* file, Int_t eventno)
{
//changes the root directory in "file" to "dirname" which corresponds to event 'eventno'
//in case of success returns the pointer to directory
//else NULL
 
 if (file == 0x0)
  {
    AliErrorClass("File is null");
    return 0x0;
  }
 if (file->IsOpen() == kFALSE)
  {
    AliErrorClass("File is not opened");
    return 0x0;
  }

 TString dirname("Event");
 dirname+=eventno;
 AliDebugClass(1, Form("Changing Dir to %s in file %s.",dirname.Data(),file->GetName()));

 Bool_t result;
 
 TDirectory* dir = dynamic_cast<TDirectory*>(file->Get(dirname));

 if (dir == 0x0)
  {
    AliDebugClass(1, Form("Can not find directory %s in file %s, creating...",
		     dirname.Data(),file->GetName()));
    
    if (file->IsWritable() == kFALSE)
     {
       AliErrorClass(Form("Can not create directory. File %s in not writable.",
		     file->GetName()));
       return 0x0;
     }
            
    TDirectory* newdir = file->mkdir(dirname);
    if (newdir == 0x0)
     {
       AliErrorClass(Form("Failed to create new directory in file %s.",
		     file->GetName()));
       return 0x0;
     }
    result = file->cd(dirname);
    if (result == kFALSE)
     {
       return 0x0;
     }
  }
 else
  {
   file->cd();//make a file active 
   file->cd(dirname);//cd event dir
  }

 return gDirectory;
}
/*****************************************************************************/ 

TString AliLoader::GetUnixDir() const
 {
 //This Method will manage jumping through unix directories in case of 
 //run with more events per run than defined in gAlice
 
   TString dir;
   
   return dir;
 }
/*****************************************************************************/ 
/************************************************************/

void AliLoader::MakeTree(Option_t *option)
 {
//Makes a tree depending on option 
//   H: - Hits
//   D: - Digits
//   S: - Summable Digits
//   R: - Reconstructed Points (clusters)
//   T: - Tracks (tracklets)

  const char *oH = strstr(option,"H");
  const char *oD = strstr(option,"D");
  const char *oS = strstr(option,"S");
  const char *oR = strstr(option,"R");
  const char *oT = strstr(option,"T");
  const char *oP = strstr(option,"P");
  
  if (oH) MakeHitsContainer();
  if (oD) MakeDigitsContainer();
  if (oS) MakeSDigitsContainer();
  if (oR) MakeRecPointsContainer();
  if (oT) MakeTracksContainer();
  if (oP) MakeRecParticlesContainer();
 }

/*****************************************************************************/ 
Int_t  AliLoader::WriteHits(Option_t* opt) const
 {
   // Writes hits
   AliDataLoader* dl = GetHitsDataLoader();
   Int_t ret = dl->WriteData(opt);
   return ret;
 }
/*****************************************************************************/ 

Int_t AliLoader::WriteSDigits(Option_t* opt) const
 {
   // Writes summable digits
   AliDataLoader* dl = GetSDigitsDataLoader();
   Int_t ret = dl->WriteData(opt);
   return ret;
 }
 
/*****************************************************************************/ 

TObject** AliLoader::GetDetectorDataRef(TObject *obj)
{
  // Returns pointer to an entry in the list of folders pointing to "obj"
 if (obj == 0x0)
  {
    return 0x0;
  }
 return GetDetectorDataFolder()->GetListOfFolders()->GetObjectRef(obj) ;
}

/*****************************************************************************/ 
void AliLoader::CleanFolders()
 {
 //Cleans all posted objects == removes from folders and deletes them
   TIter next(fDataLoaders);
   AliDataLoader* dl;
   while ((dl = (AliDataLoader*)next()))
    { 
      AliDebug(1, Form("name = %s cleaning",dl->GetName()));
      dl->Clean();
    }
 }

/*****************************************************************************/ 

Int_t AliLoader::ReloadAll()
{
  // Calling Reload function for all the data loaders
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 
 while((dl = (AliDataLoader*)next()))
  {
   Int_t err = dl->Reload();
   if (err)
    {
      AliError(Form("Reload returned error for %s",dl->GetName()));
      return err;
    }
  }
 return 0;
}
/*****************************************************************************/ 

void AliLoader::CloseFiles()
{
//close files for data loaders
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 while((dl = (AliDataLoader*)next()))
  {
   dl->CloseFile();
  }
}
/*****************************************************************************/ 

Int_t  AliLoader::SetEventFolder(TFolder* eventfolder)
{
  //sets the event folder
 if (eventfolder == 0x0)
  {
    AliError("Stupid joke. Argument is NULL");
    return 1;
  }

 fEventFolder = eventfolder;
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 
 while((dl = (AliDataLoader*)next()))
  {
    dl->SetEventFolder(fEventFolder);
    dl->SetFolder(GetDetectorDataFolder()); //Must exists - ensure register is called before
  }

 return 0;
}//sets the event folder
/*****************************************************************************/ 

Int_t AliLoader::Register(TFolder* eventFolder)
{
//triggers creation of subfolders for a given detector
//this method is called when session is read from disk
//
//warning: AliDetector in constructor (not default) calls
//creation of folder structure as well (some detectors needs folders 
//alrady in constructors)

 AliDebug(1, Form("Name is %s.",GetName()));
 if (eventFolder == 0x0)
  {
    AliError("Event folder is not set.");
    return 1;
  }
 Int_t retval = AliConfig::Instance()->AddDetector(eventFolder,fDetectorName,fDetectorName);
 if(retval)
  {
    AliError(Form("Can not create folders for %s. Event folder name is %s",
		  fDetectorName.Data(),eventFolder->GetName()));
    return retval;
  }
 SetEventFolder(eventFolder);
 return 0;
}
/*****************************************************************************/ 
AliRunLoader* AliLoader::GetRunLoader()
{
//gets the run-loader from event folder
  AliRunLoader* rg = 0x0;
  TObject * obj = GetEventFolder()->FindObject(AliRunLoader::GetRunLoaderName());
  if (obj) rg = dynamic_cast<AliRunLoader*>(obj);
  return rg;
}
/*****************************************************************************/ 
Bool_t  AliLoader::TestFileOption(Option_t* opt)
{
//tests the TFile::Option
//if file is truncated at opening moment ("recreate", "new" or "create") returns kFALSE;
//else kTRUE (means opened with "read" or "update" mode)
  TString option(opt);
  if (option.CompareTo("recreate",TString::kIgnoreCase) == 0) return kFALSE;
  if (option.CompareTo("new",TString::kIgnoreCase) == 0) return kFALSE;
  if (option.CompareTo("create",TString::kIgnoreCase) == 0) return kFALSE;
  return kTRUE;
}
/*****************************************************************************/ 
void  AliLoader::SetDirName(TString& dirname)
{
//adds "dirname" to each file. Dirname should end with "#","/", or ":" 
  TIter next(fDataLoaders);
  AliDataLoader* dl;
  while((dl = (AliDataLoader*)next()))
   {
    dl->SetDirName(dirname);
   }
}

/*****************************************************************************/ 

void AliLoader::SetDigitsFileNameSuffix(const TString& suffix) const
{
  //adds the suffix before ".root", 
  //e.g. TPC.Digits.root -> TPC.DigitsMerged.root
  //made on Jiri Chudoba demand
  GetDigitsDataLoader()->SetFileNameSuffix(suffix);
}
/*****************************************************************************/ 

void AliLoader::SetCompressionLevel(Int_t cl)
{
//sets comression level for data defined by di
  TIter next(fDataLoaders);
  AliDataLoader* dl;
  while((dl = (AliDataLoader*)next()))
   {
     dl->SetCompressionLevel(cl);
   }
}
/*****************************************************************************/ 

void AliLoader::Clean()
{
//Cleans all data loaders
  TIter next(fDataLoaders);
  AliDataLoader* dl;
  while((dl = (AliDataLoader*)next()))
   {
     dl->Clean();
   }
}
/*****************************************************************************/

void AliLoader::Clean(const TString& name)
{
  // Removes object with "name" from the detector's data folder
  // and from the memory
  TObject* obj = GetDetectorDataFolder()->FindObject(name);
  if(obj)
   {
     AliDebug(1, Form("name=%s, cleaning %s.",GetName(),name.Data()));
     GetDetectorDataFolder()->Remove(obj);
     delete obj;
   }
}

/*****************************************************************************/ 

Bool_t AliLoader::IsOptionWritable(const TString& opt)
{
  // Returns "true" if the option means also "writable"
  if (opt.CompareTo("recreate",TString::kIgnoreCase)) return kTRUE;
  if (opt.CompareTo("new",TString::kIgnoreCase)) return kTRUE;
  if (opt.CompareTo("create",TString::kIgnoreCase)) return kTRUE;
  if (opt.CompareTo("update",TString::kIgnoreCase)) return kTRUE;
  return kFALSE;
}

/*****************************************************************************/ 

void AliLoader::SetDebug(Int_t deb)
{
  // Sets debug level
  AliLog::SetClassDebugLevel("AliRunLoader", deb);
  AliLog::SetClassDebugLevel("AliLoader", deb);
  AliLog::SetClassDebugLevel("AliDataLoader", deb);
  AliLog::SetClassDebugLevel("AliBaseLoader", deb);
  AliLog::SetClassDebugLevel("AliObjectLoader", deb);
  AliLog::SetClassDebugLevel("AliTreeLoader", deb);
  AliLog::SetClassDebugLevel("AliConfig", deb);
}
/*****************************************************************************/ 

void AliLoader::SetTAddrInDet()
{
  //calls SetTreeAddress for corresponding detector
  AliRunLoader* rl = GetRunLoader();   
  if (rl == 0x0) return;
  AliRun* ar = rl->GetAliRun();
  if (ar == 0x0) return;
  AliDetector* det = ar->GetDetector(fDetectorName);  
  if (det == 0x0) return;
  det->SetTreeAddress();
}
/*****************************************************************************/ 

void AliLoader::Synchronize()
{
  //synchrinizes all writtable files 
 TIter next(fDataLoaders);
 AliDataLoader* dl;
 while ((dl = (AliDataLoader*)next()))
  {
    dl->Synchronize();
  }
  
}
/*****************************************************************************/ 
/*****************************************************************************/ 
/*****************************************************************************/ 

