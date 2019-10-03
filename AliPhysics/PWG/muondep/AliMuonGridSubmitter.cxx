/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliMuonGridSubmitter.h"

#include "AliLog.h"
#include "TFile.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include <vector>
#include "AliAnalysisRunList.h"
#include "Riostream.h"


/// Constructor
/// \param jobType type of job to handle (Acc x Eff or QA merging)
/// \param localOnly whether or not the jobs will be completely local
//______________________________________________________________________________
AliMuonGridSubmitter::AliMuonGridSubmitter(AliMuonGridSubmitter::EJobType jobType, Bool_t localOnly)
: TObject(),
fInternalMap(0x0),
fVars(0x0),
fIsValid(kFALSE),
fShouldOverwriteFiles(kFALSE),
fTemplateFileList(0x0),
fLocalFileList(0x0),
fRunList()
{
  if (!gGrid && !localOnly)
  {
    TGrid::Connect("alien://");
    if ( !gGrid )
    {
      AliError("cannot connect to grid");
    }
  }
  
  SetAliPhysicsVersion("VO_ALICE@AliPhysics::v5-08-19-01");
  
  TString basedir = gSystem->ExpandPathName("$ALICE_PHYSICS/PWG/muondep");
  
  TString dir;
  dir.Form("%s/%sTemplates",basedir.Data(),JobTypeName(jobType).Data());
  
  if (!SetTemplateDir(dir.Data()))
  {
    AliError(Form("Could not find %s directory. Please check.",dir.Data()));
  }
  
  SetLocalDir(gSystem->pwd());
}

//______________________________________________________________________________
AliMuonGridSubmitter::~AliMuonGridSubmitter()
{
  // dtor
  delete fTemplateFileList;
  delete fLocalFileList;
  delete fInternalMap;
  delete fVars;
}

/// Add include path for Root ACliC, insuring there's no duplicate
///______________________________________________________________________________
void AliMuonGridSubmitter::AddIncludePath(const char* pathList) const
{
  TObjArray* paths = TString(pathList).Tokenize(" ");
  TIter next(paths);
  TObjString* p;
  TString includePath = gSystem->GetIncludePath();
  
  while ( ( p = static_cast<TObjString*>(next()) ) )
  {
    if ( !includePath.Contains(p->String()) )
    {
      gSystem->AddIncludePath(p->String().Data());
    }
  }
  
  delete paths;
}


/// Add a file to the list of templates
/// and update the local file list if needed
//______________________________________________________________________________
void AliMuonGridSubmitter::AddToTemplateFileList(const char* filename, Bool_t alsoForMerging)
{
  
  TObjArray* a = TemplateFileList();
  
  if ( !a->FindObject(filename) )
  {
    TObjString *o = new TObjString(filename);
    o->SetBit(BIT(23), alsoForMerging);
    a->Add(o);
    UpdateLocalFileList();
  }
}

/// Add a file to the list of local files
//______________________________________________________________________________
void AliMuonGridSubmitter::AddToLocalFileList(const char* filename, Bool_t alsoForMerging)
{
  TObjArray* a = LocalFileList();
  
  if ( !a->FindObject(filename) )
  {
    TObjString *o = new TObjString(filename);
    o->SetBit(BIT(23), alsoForMerging);
    a->Add(o);
  }
}

/// Check whether file can be compiled or not
/// FIXME: use gSystem->TempFileName for tmpfile !
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CheckCompilation(const char* file) const
{
  Bool_t rv(kTRUE);
  
  TString sfile(gSystem->BaseName(file));
  TString tmpfile(Form("tmpfile_%s",sfile.Data()));
  
  gSystem->Exec(Form("cp %s %s",file,tmpfile.Data()));
  
  ReplaceVars(tmpfile.Data());

  if (!gSystem->CompileMacro(Form("%s",tmpfile.Data()), "fc"))
  {
    AliError(Form("macro %s can not be compiled. Please check.",file));
    rv = kFALSE;
  }
  
  gSystem->Exec(Form("rm %s",tmpfile.Data()));
  
  return rv;
}


/// Check whether all required local files are there
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CheckLocal() const
{
  TIter next(LocalFileList());
  TObjString* file;
  
  while ( ( file = static_cast<TObjString*>(next())) )
  {
      if ( gSystem->AccessPathName(file->String().Data()) )
      {
        return kFALSE;
      }
  }
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CheckRemote() const
{
  /// Check whether all required remote files are there
  AliWarning("implement me");
  return kFALSE;
}

/// Clean (remove) local generated files
/// \param excludeList contains a list of comma separated pattern of files
/// to be avoided (i.e. that will _not_ be removed)
//______________________________________________________________________________
void AliMuonGridSubmitter::CleanLocal(const char* excludeList) const
{
  TIter next(LocalFileList());
  TObjString* file;
  TObjArray* excludeArray = TString(excludeList).Tokenize(",");
  
  while ( ( file = static_cast<TObjString*>(next())) )
  {
    Bool_t shouldExclude(kFALSE);
    
    TIter nextExclude(excludeArray);
    TObjString* s;
    
    while ( ( s = static_cast<TObjString*>(nextExclude()))  && !shouldExclude )
    {
      if ( file->String().Contains(s->String()) ) shouldExclude=kTRUE;
    }
    
    if (!shouldExclude)
    {
      gSystem->Unlink(file->String().Data());
    }
  }
  
  delete excludeArray;
}

/// Copy (upload) a local file to remote destination
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CopyFile(const char* localFile)
{
  TString local;
  
  if ( gSystem->IsAbsoluteFileName(localFile) )
  {
    local = localFile;
  }
  else
  {
    local = Form("%s/%s",LocalDir().Data(),gSystem->ExpandPathName(localFile));
  }
  
  if (gSystem->AccessPathName(local.Data()))
  {
    AliErrorClass(Form("Local file %s does not exist",local.Data()));
    return kFALSE;
  }
  
  TString remote;
  
  remote += RemoteDir().Data();
  remote += "/";
  
  if ( gSystem->IsAbsoluteFileName(localFile) )
  {
    TString tmp(localFile);
    tmp.ReplaceAll(GetMapValue("Snapshot"),"");
    tmp.ReplaceAll(GetMapValue("Local"),"");
    remote += tmp;
  }
  else
  {
    remote += localFile;
  }
  
  TString dirName = gSystem->DirName(remote.Data());
  
  Bool_t ok(kTRUE);
  
  if (!RemoteDirectoryExists(dirName.Data()))
  {
    ok = gGrid->Mkdir(dirName.Data(),"-p");
  }
  
  if ( ok )
  {
    AliDebugClass(1,Form("cp %s alien://%s",local.Data(),remote.Data()));
    return TFile::Cp(local.Data(),Form("alien://%s",remote.Data()));
  }
  else
  {
    return kFALSE;
  }
}

/// Check we have a grid connection and that the remote dir exists
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CheckRemoteDir() const
{
  if (RemoteDir().IsNull())
  {
    AliError("you must provide the grid location where to copy the files");
    return kFALSE;
  }
  
  if (!RemoteDirectoryExists(RemoteDir()))
  {
    AliError(Form("directory %s does not exist", RemoteDir().Data()));
    return kFALSE;
  }

  return kTRUE;
}

/// Copy all files necessary to run the simulation into remote directory
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CopyLocalFilesToRemote()
{
  TIter next(LocalFileList());
  TObjString* ftc;
    
  Bool_t allok(kTRUE);
    
  while ( ( ftc = static_cast<TObjString*>(next())) )
  {
    allok = allok && CopyFile(ftc->String());
  }
  return allok;
}

/// Copy (or generate) local files from the template ones
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::CopyTemplateFilesToLocal()
{
  if (!IsValid()) return kFALSE;

  AliDebug(1,"");

  TIter next(TemplateFileList());
  TObjString* file;
  
  Int_t err(0);
  Bool_t potentialProblem(kFALSE);
  
  while ( ( file = static_cast<TObjString*>(next())) )
  {
    if ( file->String().Contains("OCDB") )
    {
      /// OCDB snapshots are not in template
      continue;
    }

    if ( !ShouldOverwriteFiles() && !gSystem->AccessPathName(file->String().Data()) )
    {
      AliError(Form("Local file %s already exists. Remove it first if you want to update overwrite it",file->String().Data()));
      potentialProblem = kTRUE;
    }
    else
    {
      TString stemplate(Form("%s/%s",TemplateDir().Data(),file->String().Data()));
      TString slocal(Form("%s/%s",LocalDir().Data(),file->String().Data()));
      
      AliDebug(1,Form("Copying %s to %s",stemplate.Data(),slocal.Data()));
      
      Int_t c =  gSystem->CopyFile(stemplate.Data(),slocal.Data(),ShouldOverwriteFiles());
      if ( c )
      {
        Bool_t ok(kFALSE);
        if ( stemplate.EndsWith("jdl",TString::kIgnoreCase) )
        {
          ok = Generate(file->String().Data());
        }
        if (!ok)
        {
          AliError(Form("Error %d copying file %s",c,stemplate.Data()));
        }
        else
        {
          c=0;
        }
      }
      else
      {
        if ( HasVars(slocal.Data()) )
        {
          if (!ReplaceVars(slocal.Data()))
          {
            AliError("pb in ReplaceVars");
            c=1;
          }
        }
      }
      err += c;
    }
  }
  
  if ( potentialProblem )
  {
    AliWarning("At least one local file could not be overwritten. Cross-check that the local files are OK before we try to upload them to the Grid !");
    return kFALSE;
  }
  return (err==0);
}

/// Create a (local and empty) JDL file
//______________________________________________________________________________
std::ostream* AliMuonGridSubmitter::CreateJDLFile(const char* name) const
{
  AliDebugClass(1,"");
  
  TString jdl(Form("%s/%s",LocalDir().Data(),name));
  
  if ( !ShouldOverwriteFiles() && !gSystem->AccessPathName(jdl.Data()) )
  {
    AliErrorClass(Form("File %s already exists. Remove it if you want to overwrite it",jdl.Data()));
    return 0x0;
  }
  
  std::ofstream* os = new std::ofstream(gSystem->ExpandPathName(jdl.Data()));
  
  if (os->bad())
  {
    AliErrorClass(Form("Cannot create file %s",jdl.Data()));
    delete os;
    os=0x0;
  }
  
  return os;
}

/// Get the last staging phase already performed
//______________________________________________________________________________
Int_t AliMuonGridSubmitter::GetLastStage(const char* remoteDir)
{
  TGridResult* r = gGrid->Ls(remoteDir);
  Int_t i(0);
  Int_t lastStage(0);
  Int_t offset = TString("Stage_").Length();
  
  while ( r->GetFileName(i) )
  {
    TString file(r->GetFileName(i));
    if  (file.BeginsWith("Stage_") && !file.Contains("xml") )
    {
      Int_t n = TString(file(offset,file.Length()-offset)).Atoi();
      lastStage = TMath::Max(lastStage,n);
    }
    ++i;
  }
  delete r;
  return lastStage;
}

/// Convenience method to access internal map of TObjStrings
//______________________________________________________________________________
TString AliMuonGridSubmitter::GetMapValue(const char* key) const
{
  if (!fInternalMap) return "";
  
  TObject* o = fInternalMap->GetValue(key);
  
  if (o)
  {
    return static_cast<TObjString*>(o)->String();
  }
  
  return "";
}

/// Find the variables in the file
//______________________________________________________________________________
TObjArray* AliMuonGridSubmitter::GetVariables(const char* file) const
{
  std::ifstream in(file);
  char line[1024];
  TObjArray* variables(0x0);
  
  while ( in.getline(line,1023,'\n') )
  {
    TString sline(line);
    while (sline.Contains("VAR_") && !sline.BeginsWith("//") )
    {
      Int_t i1 = sline.Index("VAR_");
      Int_t i2(i1);
      
      while ( ( i2 < sline.Length() ) && ( isalnum(sline[i2]) || sline[i2]=='_' ) ) ++i2;
      
      if (!variables)
      {
        variables = new TObjArray;
        variables->SetOwner(kTRUE);
      }
      
      TString var = sline(i1,i2-i1);
      if ( !variables->FindObject(var) )
      {
        variables->Add(new TObjString(var));
      }
      sline.ReplaceAll(var,"");
    }
  }
  
  in.close();
  
  return variables;
}

/// Whether or not the file contains variable var
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::HasVar(const char* file, const char* var) const
{
  std::ifstream in(file);
  char line[1024];
  while ( in.getline(line,1023,'\n') )
  {
    TString sline(line);
    if (sline.Contains("VAR_") && !sline.BeginsWith("//") && sline.Contains(var))
    {
      return kTRUE;
    }
  }
  return kFALSE;
}

/// Whether or not the file contains variables that have to be substituted
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::HasVars(const char* file) const
{
  std::ifstream in(file);
  char line[1024];
  while ( in.getline(line,1023,'\n') )
  {
    TString sline(line);
    if (sline.Contains("VAR_") && !sline.BeginsWith("//") )
    {
      return kTRUE;
    }
  }
  return kFALSE;
}

//______________________________________________________________________________
TMap* AliMuonGridSubmitter::InternalMap() const
{
  if (!fInternalMap)
  {
    fInternalMap = new TMap;
    fInternalMap->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  return fInternalMap;
}

//______________________________________________________________________________
TString AliMuonGridSubmitter::JobTypeName(AliMuonGridSubmitter::EJobType jobType) const
{
  if ( jobType == kAccEff )
  {
    return "AccEff";
  }
  else if ( jobType == kQAMerge)
  {
    return "QAMerge";
  }
  return "unknown";
}

/// Return (after createing and filling it if needed)
/// the internal file list with paths from the local directory
//______________________________________________________________________________
TObjArray* AliMuonGridSubmitter::LocalFileList() const
{
  if (!fLocalFileList)
  {
    fLocalFileList = static_cast<TObjArray*>(TemplateFileList()->Clone());
  }
  
  return fLocalFileList;
}

/// Number of runs we are dealing with
//______________________________________________________________________________
UInt_t AliMuonGridSubmitter::NofRuns() const
{
  return fRunList.size();
}

/// Return an array where the map's keys are sorted alphabetically
/// the returned array should be deleted by the client
//______________________________________________________________________________
TObjArray* AliMuonGridSubmitter::OrderKeys(const TMap& map) const
{
  TObjArray* keyArray = new TObjArray;
  keyArray->SetOwner(kTRUE);
  TObjString* key;
  
  TIter next(&map);
  while ( ( key = static_cast<TObjString*>(next())) )
  {
    keyArray->Add(new TObjString(key->String()));
  }
  
  keyArray->Sort();
  return keyArray;
}


/// output to ostream of key,{values} pair
//______________________________________________________________________________
void AliMuonGridSubmitter::OutputToJDL(std::ostream& out, const char* key,
                                    const TObjArray& values) const
{
  out << key << " = ";
  
  Int_t n = values.GetEntries();
  
  if ( n > 1 )
  {
    out << "{" << std::endl;
    TIter next(&values);
    TObjString* v;
    
    while ( ( v = static_cast<TObjString*>(next())) )
    {
      if ( v->String().Length() == 0 ) continue;
      
      --n;
      out << "\t\"" << v->String().Data() << "\"";
      if  ( n ) out << ",";
      out << std::endl;
    }
    out << "}";
  }
  else
  {
    TString& v1 = static_cast<TObjString*>(values.At(0))->String();
    
    if ( v1.IsDigit() && !(TString(key).Contains("Price")) && !(TString(key).Contains("SplitMax")) &&
        !(TString(key).Contains("TTL")) && !(TString(key).Contains("Arguments")) )
    {
      out << v1.Atoi();
    }
    else
    {
      out << "\"" << v1.Data() << "\"";
    }
  }
  out << ";" << std::endl;
}

/// output to ostream
//______________________________________________________________________________
void AliMuonGridSubmitter::OutputToJDL(std::ostream& out, const char* key, const char* v1,
                                    const char* v2, const char* v3, const char* v4,
                                    const char* v5, const char* v6, const char* v7,
                                    const char* v8, const char* v9) const
{
  TObjArray values;
  values.SetOwner(kTRUE);
  
  if ( strlen(v1) > 0 ) values.Add(new TObjString(v1));
  if ( strlen(v2) > 0 ) values.Add(new TObjString(v2));
  if ( strlen(v3) > 0 ) values.Add(new TObjString(v3));
  if ( strlen(v4) > 0 ) values.Add(new TObjString(v4));
  if ( strlen(v5) > 0 ) values.Add(new TObjString(v5));
  if ( strlen(v6) > 0 ) values.Add(new TObjString(v6));
  if ( strlen(v7) > 0 ) values.Add(new TObjString(v7));
  if ( strlen(v8) > 0 ) values.Add(new TObjString(v8));
  if ( strlen(v9) > 0 ) values.Add(new TObjString(v9));
  
  OutputToJDL(out,key,values);
}

/// Printout
//______________________________________________________________________________
void AliMuonGridSubmitter::Print(Option_t* opt) const
{
  if (!IsValid())
  {
    std::cout << std::string(80,'*') << std::endl;
    std::cout << "INVALID OBJECT. CHECK BELOW THE CONFIGURATION." << std::endl;
    std::cout << std::string(80,'*') << std::endl;
  }

  TString sopt(opt);
  sopt.ToUpper();

  std::cout << "-- Internals : " << std::endl;
  
  TObjArray* im = OrderKeys(*fInternalMap);
  
  TIter next(im);
  TObjString* key;
  
  while ( ( key = static_cast<TObjString*>(next()) ) )
  {
    TString value = static_cast<TObjString*>(fInternalMap->GetValue(key->String()))->String();
    
    std::cout << key->String() << " : " << value.Data() << std::endl;
  }
  
  delete im;
  
  if ( NofRuns() )
  {
    std::cout << NofRuns() << " run";
    if ( NofRuns() > 1 ) std::cout << "s";
    std::cout << " = ";
    for ( std::vector<int>::size_type i = 0; i < fRunList.size(); ++i )
    {
      std::cout << fRunList[i] << " ";
    }
    std::cout << std::endl;
  }

  PrintVariables(sopt.Contains("ALLVARS"));

  std::cout << std::endl << "-- Files to be uploaded:" << std::endl;
  TIter nextFile(LocalFileList());
  TObjString* sfile;
  while ( ( sfile = static_cast<TObjString*>(nextFile())) )
  {
    std::cout << sfile->String().Data() << std::endl;
  }
}

/// Show the variables 
/// \param all set it to true to show all variables, not just those which exist in the list of template files 
//______________________________________________________________________________
void AliMuonGridSubmitter::PrintVariables(Bool_t all) const
{
  std::cout << std::endl << "-- Variables : " << std::endl;

  TObjArray* iv = OrderKeys(*fVars);
  
  TIter nextVar(iv);
  TObjString* key;

  while ( ( key = static_cast<TObjString*>(nextVar())) )
  {
      Bool_t showIt = kTRUE;

    if ( !all )
    {
        showIt = kFALSE;
        // check first the variable is used somewhere in the template file list
        TIter next(TemplateFileList());
        TObjString* file;
        while ( ( file = static_cast<TObjString*>(next())))
        {
            TString filename = TemplateDir();
            filename += "/";
            filename += file->String();
            if ( HasVar(filename.Data(),key->String().Data())) 
            {
                showIt = kTRUE;
                break;
            }
        }
    }

    if ( showIt )
    {
        TObjString* value = static_cast<TObjString*>(fVars->GetValue(key->String()));
        std::cout << "Variable " << key->String() << " will be replaced by " << value->String() << std::endl;
    }
  }
  
  delete iv;
}

/// Returns true if dirname exists. Can be also a path.
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::RemoteDirectoryExists(const char *dirname) const
{
  if (!gGrid) return kFALSE;
  // Check if dirname is a path
  TString dirstripped = dirname;
  dirstripped = dirstripped.Strip();
  dirstripped = dirstripped.Strip(TString::kTrailing, '/');
  TString dir = gSystem->BaseName(dirstripped);
  dir += "/";
  TString path = gSystem->DirName(dirstripped);
  TGridResult *res = gGrid->Ls(path, "-F");
  if (!res) return kFALSE;
  TIter next(res);
  TMap *map;
  TObject *obj;
  while ((map=dynamic_cast<TMap*>(next()))) {
    obj = map->GetValue("name");
    if (!obj) break;
    if (dir == obj->GetName()) {
      delete res;
      return kTRUE;
    }
  }
  delete res;
  return kFALSE;
}

/// Returns true if lfn exists.
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::RemoteFileExists(const char *lfn)
{
  if (!gGrid) return kFALSE;
  TGridResult *res = gGrid->Ls(lfn);
  if (!res) return kFALSE;
  TMap *map = dynamic_cast<TMap*>(res->At(0));
  if (!map) {
    delete res;
    return kFALSE;
  }
  TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("name"));
  if (!objs || !objs->GetString().Length()) {
    delete res;
    return kFALSE;
  }
  delete res;
  return kTRUE;
}

/// Replace the variables (i.e. things starting by VAR_) found in file
/// by their value
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::ReplaceVars(const char* file) const
{
  AliDebug(1,file);
  
  std::ifstream in(file);
  char line[1024];
  TObjArray lines;
  lines.SetOwner(kTRUE);
  Int_t nvars(0);
  Int_t nreplaced(0);

  TIter next(fVars);

  while ( in.getline(line,1023,'\n') )
  {
    TString sline(line);
    while (sline.Contains("VAR_") && !sline.BeginsWith("//") )
    {
      ++nvars;
      TObjString* key;
      next.Reset();
      
      int n(0);
      
      while ( ( key = static_cast<TObjString*>(next())) )
      {
        AliDebug(1,Form("Does %s contains %s ?",sline.Data(),key->String().Data()));
        
        if ( sline.Contains(key->String()) )
        {
          ++nreplaced;
          ++n;
          TObjString* value = static_cast<TObjString*>(fVars->GetValue(key->String()));
          AliDebug(1,Form("Replacing %s by %s in %s",key->String().Data(),value->String().Data(),sline.Data()));
          sline.ReplaceAll(key->String(),value->String());
          AliDebug(1,"Done");
          break;
        }
      }
      
      if (n==0)
      {
        AliError(Form("Could not find a replacement for variable %s in file %s",sline.Data(),file));
        break;
      }
    }

    lines.Add(new TObjString(sline));
  }
  
  in.close();
  
  if ( nvars > 0 )
  {
    if ( nreplaced != nvars )
    {
      AliError(Form("nvars=%d nreplaced=%d",nvars,nreplaced));
      return kFALSE;
    }
    std::ofstream out(file);
    TIter nextLine(&lines);
    TObjString* s;
    while ( ( s = static_cast<TObjString*>(nextLine()) ) )
    {
      out << s->String().Data() << std::endl;
    }
    out.close();
  }
  
  AliDebug(1,Form("replaced %d vars",nvars));
  
  return kTRUE;
}

/// Return a reference to our runlist
//______________________________________________________________________________
const std::vector<int>& AliMuonGridSubmitter::RunList() const
{
  return fRunList;
}

/// Set the AliPhysics package to be used by the jobs
/// the corresponding aliroot, root, geant3 versions
/// should be set automatically by alien.
//______________________________________________________________________________
void AliMuonGridSubmitter::SetAliPhysicsVersion(const char* aliphysics)
{
  SetMapKeyValue("AliPhysics",aliphysics);
}

/// Set the AliRoot package to be used by the jobs
/// the corresponding root, geant3 versions
/// should be set automatically by alien.
//______________________________________________________________________________
void AliMuonGridSubmitter::SetAliRootVersion(const char* aliroot)
{
  SetMapKeyValue("AliRoot",aliroot);
}

/// Set the generator package to be used by the jobs
//______________________________________________________________________________
void AliMuonGridSubmitter::SetGeneratorPackage(const char* generator)
{
  SetMapKeyValue("Generator",generator);
}

//______________________________________________________________________________
void AliMuonGridSubmitter::SetPackages(const char* aliroot,
                                         const char* root,
                                         const char* geant3,
                                         const char* api)
{
  /// @deprecated
  /// Set the packages to be used by the jobs
  /// If root and geant3 are given (default is to let alien get the correct
  /// values for the requested aliroot version), then they must
  /// correspond to a valid combination, see http://alimonitor.cern.ch/packages/
  
  AliWarning("This method is deprecated and will disappear at some point. Please use SetAliRootVersion or SetAliPhysicsVersion instead");
  
  SetMapKeyValue("AliRoot",aliroot);
  SetMapKeyValue("Root",root);
  SetMapKeyValue("Geant3",geant3);
  SetMapKeyValue("API",api);
}

/// Set the target remote directory (on the grid)
//______________________________________________________________________________
TString AliMuonGridSubmitter::GetRemoteDir(const char* dir, Bool_t create)
{
  if (!RemoteDirectoryExists(dir))
  {
    if (!create)
    {
      AliErrorClass(Form("Remote directory %s does not exist", dir));
      return "";
    }
    else
    {
      AliInfoClass(Form("Remote directory %s does not exist. Trying to create it...",dir));
      if (!gGrid)
      {
          AliErrorClass("cannot connect to grid");
          return "";
      }
      if ( !gGrid->Mkdir(dir,"-p") )
      {
        AliErrorClass(Form("Could not create remote dir. Sorry."));
        return "";
      }
    }
  }
  return dir;
}

//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::SetLocalDirectory(const char* type, const char* path)
{
  if (gSystem->AccessPathName(path)==kFALSE)
  {
    SetMapKeyValue(type,path);
    return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
void AliMuonGridSubmitter::SetMapKeyValue(const char* key, const char* value)
{
  TObjString skey(key);
  InternalMap()->Remove(&skey);
  InternalMap()->Add(new TObjString(key),new TObjString(value));
}

//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::SetRemoteDirectory(const char* type, const char* path)
{
  // Set the merged directory to be used
  TString v = GetRemoteDir(path,kTRUE);
  SetMapKeyValue(type,v);
  return (v.Length()>0);
}

//______________________________________________________________________________
void AliMuonGridSubmitter::SetRunList(const char* runlist)
{
    // set the runlist from a text file
  AliAnalysisRunList rl(runlist);
  fRunList = rl.AsVector();
}

//______________________________________________________________________________
void AliMuonGridSubmitter::SetRunList(int runNumber)
{
  // set the runlist from a text file
  fRunList.clear();
  fRunList.push_back(runNumber);
}

//______________________________________________________________________________
TString AliMuonGridSubmitter::GetVar(const char* key) const
{
  TObjString* o = static_cast<TObjString*>(Vars()->GetValue(key));
  if (o)
  {
    return o->String();
  }
  return "";
}

/// Set the value of a variable
///
/// Pay attention to how string variables should be given here : you have
/// to espace the quotation marks :
///
/// SetVar("VAR_PYTHIA8_SETUP_STRINGS","\"SoftQCD:doubleDiffractive=off\"");
//______________________________________________________________________________
Bool_t AliMuonGridSubmitter::SetVar(const char* varname, const char* value)
{
  TString s(varname);
  s.ToUpper();
  if (!s.BeginsWith("VAR_"))
  {
    AliError("Variable name should start with VAR_");
    return kFALSE;
  }
  if (!fVars)
  {
    fVars = new TMap;
    fVars->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  
  TObject* o = new TObjString(s);
  fVars->Remove(o);
  
  fVars->Add(o,new TObjString(value));
  
  return kTRUE;
}

/// Return (after creating if needed) the list of template files
//______________________________________________________________________________
TObjArray* AliMuonGridSubmitter::TemplateFileList() const
{
  if (!fTemplateFileList)
  {
    fTemplateFileList = new TObjArray;
    fTemplateFileList->SetOwner(kTRUE);
  }
  return fTemplateFileList;
}

/// Insure that local file list contains at least all of the template files
//______________________________________________________________________________
void AliMuonGridSubmitter::UpdateLocalFileList()
{
  TIter next(TemplateFileList());
  TObjString* s;
  
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    AddToLocalFileList(s->String().Data(), s->TestBit(BIT(23)));
  }
}

//______________________________________________________________________________
TMap* AliMuonGridSubmitter::Vars() const
{
  if (!fVars)
  {
    fVars = new TMap;
    fVars->SetOwnerKeyValue(kTRUE,kTRUE);
  }
  return fVars;
}

