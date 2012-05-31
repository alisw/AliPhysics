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
//                                                                                         //
//  EStdBasicLoaders   idx     Object Type        Description                              //
//      kData           0    TTree or TObject     main data itself (hits,digits,...)       //
//                                                                                         //
//                                                                                         //
//  User can define and add more basic loaders even Run Time.                              //
//  Caution: in order to save information about added base loader                          //
//  user must rewrite Run Loader to galice.file, overwriting old setup                     //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include "AliDataLoader.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliObjectLoader.h"
#include "AliTreeLoader.h"
#include "AliLog.h"

#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TBits.h>

ClassImp(AliDataLoader)

//______________________________________________________________________________
AliDataLoader::AliDataLoader():
 fFileName(0),
 fFile(0x0),
 fDirectory(0x0),
 fFileOption(),
 fCompressionLevel(2),
 fNEventsPerFile(0),
 fBaseLoaders(0x0),
 fEventFolder(0x0),
 fFolder(0x0)
{
  
}

//______________________________________________________________________________
AliDataLoader::AliDataLoader(const char* filename, const char* contname, 
			     const char* name, Option_t* opt):
 TNamed(name,name),
 fFileName(filename),
 fFile(0x0),
 fDirectory(0x0),
 fFileOption(0),
 fCompressionLevel(2),
 fNEventsPerFile(0),
 fBaseLoaders(new TObjArray(4)),
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

//______________________________________________________________________________
AliDataLoader::~AliDataLoader()
{
  //
  //dtor
  //
  UnloadAll();
}

//______________________________________________________________________________
Int_t  AliDataLoader::SetEvent()
{
  //
  // The same that GetEvent but do not post data to folders
  //
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

//______________________________________________________________________________
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

//______________________________________________________________________________
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
	  TString option1 = fFile->GetOption();
	  if (option1.CompareTo("read",TString::kIgnoreCase) == 0)
	    {
	      AliInfo(Form("File %s already opened in read mode.",fFile->GetName()));
	    }
	  else
	    {
	      TString option2 = opt;
	      if (option2.CompareTo("read",TString::kIgnoreCase) == 0)
		{
		  AliInfo(Form("Open already opened file %s in read mode.",fFile->GetName()));
		}
	      else {
		AliWarning(Form("File %s already opened by sombody else. First close it.",
				fFile->GetName()));
		return 0;
	      }
	    }
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
  
  fFile->SetBit(TFile::kDevNull);
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

//______________________________________________________________________________
void AliDataLoader::Unload()
{
  //
  //unloads main data -  shortcut method 
  //
  GetBaseLoader(0)->Unload();
}

//______________________________________________________________________________
void AliDataLoader::UnloadAll()
{
  //
  // Unloads all data 
  //
  if ( fFile == 0x0 ) return; //nothing loaded
  
  TIter next(fBaseLoaders);
  AliBaseLoader* bl;
  while ((bl = (AliBaseLoader*)next()))
    {
      bl->Unload();
    }
}


//______________________________________________________________________________
Int_t AliDataLoader::Reload()
{
  //
  // Unloads and loads data again
  //
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


//______________________________________________________________________________
Int_t AliDataLoader::WriteData(Option_t* opt)
{
  //
  // Writes primary data ==  first BaseLoader
  //
  AliDebug(1, Form("Writing %s container for %s data. Option is %s.",
		   GetBaseLoader(0)->GetName(),GetName(),opt));
  return GetBaseLoader(0)->WriteData(opt);
}

//______________________________________________________________________________
Int_t AliDataLoader::Load(Option_t* opt)
{
  //
  // Writes primary data ==  first BaseLoader
  //
  return GetBaseLoader(0)->Load(opt);
}

//______________________________________________________________________________
Int_t  AliDataLoader::SetEventFolder(TFolder* eventfolder)
{
  //
  // Sets the event folder
  //
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


//______________________________________________________________________________
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

//______________________________________________________________________________
TFolder* AliDataLoader::GetEventFolder()
{
  //
  // Get EVENT folder
  // Data that are changing from event to event, even in single run
  //
  AliDebug(1, "EF = %#x");
  return fEventFolder;
}

//______________________________________________________________________________
AliRunLoader* AliDataLoader::GetRunLoader()
{
  //
  // Gets the run-loader from event folder
  //
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

//______________________________________________________________________________
void AliDataLoader::CloseFile()
{
  //
  // Closes file
  //
  TIter next(fBaseLoaders);
  AliBaseLoader* bl;
  while ((bl = (AliBaseLoader*)next()))
    {
      if (bl->IsLoaded()) return;
    }
  
  AliDebug(1, "Closing (object) file.");
  
  if (fFile) {
    fFile->Close("R");
    delete fFile;
    fFile = 0x0;
  }
  fDirectory = 0x0;
}


//______________________________________________________________________________
void AliDataLoader::Clean()
{
  //
  // Cleans main data
  //
  GetBaseLoader(0)->Clean();
}

//______________________________________________________________________________
void AliDataLoader::CleanAll()
{
  //
  // Cleans all folders 
  //
  TIter next(fBaseLoaders);
  AliBaseLoader* bl;
  while ((bl = (AliBaseLoader*)next()))
    {
      bl->Clean();
    }
}

//______________________________________________________________________________
void AliDataLoader::SetFileNameSuffix(const TString& suffix)
{
  //
  // adds the suffix before ".root", 
  // e.g. TPC.Digits.root -> TPC.DigitsMerged.root
  // made on Jiri Chudoba demand
  //
  AliDebug(1, Form("suffix=%s",suffix.Data()));
  AliDebug(1, Form("   Digits File Name before: %s",fFileName.Data()));
  
  static const TString dotroot(".root");
  const TString& suffixdotroot = suffix + dotroot;
  fFileName = fFileName.ReplaceAll(dotroot,suffixdotroot);
  
  AliDebug(1, Form("                    after : %s",fFileName.Data()));
}

//______________________________________________________________________________
Bool_t AliDataLoader::CheckReload()
{
  //
  // Checks if we have to reload given file
  //
  if (fFile == 0x0) return kFALSE;
  TString tmp = SetFileOffset(fFileName);
  if (tmp.CompareTo(fFile->GetName())) return kTRUE;  //file must be reloaded
  return  kFALSE;
}

//______________________________________________________________________________
const TString AliDataLoader::SetFileOffset(const TString& fname)
{
  //
  // Return fname
  //
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

//______________________________________________________________________________
void AliDataLoader::SetFileOption(Option_t* newopt)
{
  //
  // Sets file option
  //
  if (fFileOption.CompareTo(newopt) == 0) return;
  fFileOption = newopt;
  Reload();
}

//______________________________________________________________________________
void AliDataLoader::SetCompressionLevel(Int_t cl)
{
  //
  // Sets comression level for data defined by di
  //
  fCompressionLevel = cl;
  if (fFile) fFile->SetCompressionLevel(cl);
}

//______________________________________________________________________________
void AliDataLoader::MakeTree()
{
  //
  // Makes tree for the current data loader
  //
  AliTreeLoader* tl = dynamic_cast<AliTreeLoader*>(fBaseLoaders->At(0));
  if (tl == 0x0)
   {
     AliError("Can not make a tree because main base loader is not a tree loader");
     return;
   }
  tl->MakeTree();
}

//______________________________________________________________________________
Bool_t AliDataLoader::IsFileWritable() const
{
  //
  // Returns true if file is writable
  //
  return (fFile)?fFile->IsWritable():kFALSE;
}

//______________________________________________________________________________
Bool_t AliDataLoader::IsFileOpen() const
{
  //
  // Returns true if file is writable
  //
  return (fFile)?fFile->IsOpen():kFALSE;
}

//______________________________________________________________________________
Bool_t AliDataLoader::IsOptionContrary(const TString& option) const
{
  // Checks if passed option is contrary with file open option 
  // which is passed option "writable" and existing option not wriable
  // in reverse case it is no harm so it is NOT contrary
  if (fFile == 0x0) return kFALSE; //file is not opened - no problem
  
  if ( ( AliLoader::IsOptionWritable(option)      == kTRUE  ) &&     // passed option is writable and 
       ( AliLoader::IsOptionWritable(fFileOption) == kFALSE )    )   // existing one is not writable
    {
      return kTRUE;
    }
  
  return kFALSE;
}


//______________________________________________________________________________
void AliDataLoader::AddBaseLoader(AliBaseLoader* bl)
{
  //Adds a base loader to lits of base loaders managed by this data loader
  //Managed data will be stored in proper root directory,
  //and posted to 
  // - in case of tree/object - data folder connected with detector associated with this data loader
  
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

//______________________________________________________________________________
AliBaseLoader* AliDataLoader::GetBaseLoader(const TString& name) const
{
  //
  // Return pointer to base loader
  //
  return dynamic_cast<AliBaseLoader*>(fBaseLoaders->FindObject(name));
}

//______________________________________________________________________________
AliBaseLoader* AliDataLoader::GetBaseLoader(Int_t n) const
{
  //
  // Gets the n-th base loader (what is n?)
  //
  return dynamic_cast<AliBaseLoader*>(fBaseLoaders->At(n));
}

//______________________________________________________________________________
TTree* AliDataLoader::Tree() const
{
  // Returns tree from the main base loader
  // it is just shortcut method for comfort of user
  // main storage object does not have to be Tree  - 
  // that is why first we need to check if it is a TreeLoader 
  AliTreeLoader* tl = dynamic_cast<AliTreeLoader*>(GetBaseLoader(0));
  if (tl == 0x0) return 0x0;
  return tl->Tree();
}

//______________________________________________________________________________
void  AliDataLoader::SetDirName(TString& dirname)
{
  //
  // Sets the directory name where the files will be stored
  //
  AliDebug(10, Form("FileName before %s",fFileName.Data()));
  Int_t n = fFileName.Last('/');
  AliDebug(10, Form("Slash found on pos %d",n));
  if (n > 0) fFileName = fFileName.Remove(0,n+1);
  AliDebug(10, Form("Core FileName %s",fFileName.Data()));
  fFileName = dirname + fFileName;
  AliDebug(10, Form("FileName after %s",fFileName.Data()));
}

//______________________________________________________________________________
AliObjectLoader* AliDataLoader::GetBaseDataLoader()
{
  //
  // Gets the base data loader
  //
  return dynamic_cast<AliObjectLoader*>(GetBaseLoader(kData));
}

//______________________________________________________________________________
void AliDataLoader::SetBaseDataLoader(AliBaseLoader* bl)
{
  //
  // Sets data base loader
  //
  if (bl == 0x0)
    {
      AliError("Parameter is null");
      return;
    }
  if (GetBaseDataLoader()) delete GetBaseDataLoader();
  fBaseLoaders->AddAt(bl,kData);
}

//______________________________________________________________________________
void AliDataLoader::Synchronize()
{
  //
  // Synchronizes all writable files 
  //
  if ( fFile ) fFile->Flush();
}



