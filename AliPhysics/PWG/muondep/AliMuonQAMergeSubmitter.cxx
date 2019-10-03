#include "AliMuonQAMergeSubmitter.h"

#include "AliMuonAccEffSubmitter.h"
#include "AliAnalysisTriggerScalers.h"
#include "AliLog.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TGrid.h"
#include "TGridResult.h"

namespace {
  const Int_t kFinal(999);
}


/// \cond CLASSIMP
ClassImp(AliMuonQAMergeSubmitter)
/// \endcond

//_____________________________________________________________________________
AliMuonQAMergeSubmitter::AliMuonQAMergeSubmitter(const char* period, const char* pass) :
AliMuonGridSubmitter(AliMuonGridSubmitter::kQAMerge),
fPeriod(period),
fPass(pass),
fWhatToMerge("Merged.QA.Data.root"),
fSplitMaxInputFileNumber(50)
{
  if (!fPeriod.BeginsWith("LHC"))
  {
    AliError("Period not starting with LHC !");
  }
  
  AddToTemplateFileList("QAMerge.C");
  AddToTemplateFileList("QAMerge.sh");
  AddToTemplateFileList("validation.sh");
  AddToTemplateFileList(MergeJDLName(kFALSE).Data());
  AddToTemplateFileList(MergeJDLName(kTRUE).Data());
  
  SetVar("VAR_MERGED_OUTPUT_NAME",Form("%s",fWhatToMerge.Data()));

  ShouldOverwriteFiles(kTRUE);

  TString speriod(period);
  
  Int_t year = 2000 + TString(speriod(3,3)).Atoi();

  SetMapKeyValue("DataDirFormat",Form("/alice/data/%d/%s/%%09d/ESDs/%s",year,period,pass));
}

//_____________________________________________________________________________
AliMuonQAMergeSubmitter::~AliMuonQAMergeSubmitter()
{
}

//______________________________________________________________________________
Bool_t AliMuonQAMergeSubmitter::Generate(const char* jdlname) const
{
  /// Create the JDL for merging jobs
  
  AliDebug(1,"");
  
  std::ostream* os = CreateJDLFile(jdlname);
  
  if (!os)
  {
    return kFALSE;
  }
  
  Bool_t final = TString(jdlname).Contains("final",TString::kIgnoreCase);
  
  (*os) << "# Generated merging jdl (production mode)" << std::endl
  << "# $1 = run number" << std::endl
  << "# $2 = merging stage" << std::endl
  << "# Stage_<n>.xml made via: find <OutputDir> *Stage<n-1>/*" << fWhatToMerge.Data() << std::endl;
  
  OutputToJDL(*os,"Packages",
         GetMapValue("AliRoot").Data(),
         GetMapValue("Geant3").Data(),
         GetMapValue("Root").Data(),
         GetMapValue("API").Data());
         
  OutputToJDL(*os,"Executable","QAMerge.sh");
  
  OutputToJDL(*os,"Price","1");
  
  if ( final )
  {
    OutputToJDL(*os,"Jobtag","comment: AliMuonQAMergeSubmitter final merging RUN $1");
  }
  else
  {
    OutputToJDL(*os,"Jobtag","comment: AliMuonQAMergeSubmitter merging RUN $1 stage $2");
  }
  
  OutputToJDL(*os,"Workdirectorysize","5000MB");
  
  OutputToJDL(*os,"Validationcommand",Form("%s/validation.sh",RemoteDir().Data()));
  
  OutputToJDL(*os,"TTL","14400");
  
  OutputToJDL(*os,"OutputArchive",
         "log_archive.zip:stderr,stdout@disk=1",
         Form("root_archive.zip:%s@disk=3",fWhatToMerge.Data())
         );
  
//  OutputToJDL(*os,"Arguments",(final ? "2":"1")); // for QAmerge.sh, 1 means intermediate merging stage, 2 means final merging
  
  if ( !final )
  {
    OutputToJDL(*os,"Arguments","wn.xml");
    OutputToJDL(*os,"InputFile",Form("LF:%s/QAMerge.C",RemoteDir().Data()));
    OutputToJDL(*os,"OutputDir",Form("%s/$1/Stage_$2/#alien_counter_03i#",RemoteDir().Data()));
    OutputToJDL(*os,"Split","se");
    OutputToJDL(*os,"InputDataCollection",Form("LF:%s/$1/Stage_$2.xml,nodownload",RemoteDir().Data()));    
    OutputToJDL(*os,"SplitMaxInputFileNumber",Form("%d",GetSplitMaxInputFileNumber()));
    OutputToJDL(*os,"InputDataListFormat","xml-single");
    OutputToJDL(*os,"InputDataList","wn.xml");
  }
  else
  {
    OutputToJDL(*os,"InputFile",Form("LF:%s/QAMerge.C",RemoteDir().Data()),
           Form("LF:%s/$1/Stage_$2.xml",RemoteDir().Data()));
    OutputToJDL(*os,"Arguments","Stage_$2.xml $1");
    OutputToJDL(*os,"OutputDir",Form("%s/$1",RemoteDir().Data()));
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t AliMuonQAMergeSubmitter::MakeXMLCollectionForRun(Int_t runNumber, Int_t stage)
{
  /// Create a collection named Stage_[stage].xml
  /// with the files from stage-1 (or all the files if stage=0)
  /// Return the total number of files to be merged
  
  if ( stage < 0 )
  {
    AliError(Form("Stage (%d) should be >=0",stage));
    return 0;
  }
  
  TString filename;
  TString sourcedir;

  if ( stage > 1 )
  {
    filename.Form("%s/%d/Stage_%d.xml",LocalDir().Data(),runNumber,stage);
    sourcedir.Form("%s/%d/Stage_%d",RemoteDir().Data(),runNumber,stage-1);
  }
  else if ( stage == 1 )
  {
    sourcedir = GetMapValue("DataDirFormat");
    
    if ( sourcedir.Length() == 0 )
    {
      AliError("Cannot make collections from an empty data dir !");
      return 0;
    }
    filename.Form("%s/%d/Stage_%d.xml",LocalDir().Data(),runNumber,stage);
  }
  else
  {
    AliError("oups");
    return 0;
  }
  
  UInt_t count(0);

  TGridResult* res = gGrid->Query(Form(sourcedir.Data(),runNumber),fWhatToMerge.Data());
  
  Int_t nFiles = res->GetEntries();
  
  if (!nFiles)
  {
    AliError(Form("Got no file for run %d",runNumber));
    return 0;
  }
  
  gSystem->mkdir(Form("%s/",gSystem->DirName(filename.Data())),kTRUE);
  
  AliDebug(1,Form("Creating %s",filename.Data()));
  
  std::ofstream out(filename.Data());
  
  out << Form("<?xml version=\"1.0\"?>\n<alien>\n  <collection name=\"%d-stage-%d\">",runNumber,stage) << std::endl;
  
  Long64_t size(0);
  const Double_t byte2GB(1024*1024*1024);
  
  for (Int_t i = 0; i < nFiles; ++i)
  {
    ++count;
    
    size += TString(res->GetKey(i,"size")).Atoll();
    
    out << Form("    <event name=\"%d\">",count) << std::endl;
    out << Form("      <file name=\"%s\" aclId=\"%s\" broken=\"%s\" ctime=\"%s\" "
                "dir=\"%s\" entryId=\"%s\" expiretime=\"%s\" gowner=\"%s\" "
                "guid=\"%s\" guidtime=\"%s\" lfn=\"%s\" md5=\"%s\" owner=\"%s\" "
                " perm=\"%s\" replicated=\"%s\" size=\"%s\" turl=\"%s\" type=\"%s\" />",
                gSystem->BaseName(res->GetKey(i,"lfn")),
                res->GetKey(i,"aclId"),
                res->GetKey(i,"broken"),
                res->GetKey(i,"ctime"),
                res->GetKey(i,"dir"),
                res->GetKey(i,"entryId"),
                res->GetKey(i,"expiretime"),
                res->GetKey(i,"gowner"),
                res->GetKey(i,"guid"),
                res->GetKey(i,"guidtime"),
                res->GetKey(i,"lfn"),
                res->GetKey(i,"md5"),
                res->GetKey(i,"owner"),
                res->GetKey(i,"perm"),
                res->GetKey(i,"replicated"),
                res->GetKey(i,"size"),
                res->GetKey(i,"turl"),
                res->GetKey(i,"type")) << std::endl;
    out <<      "    </event>" << std::endl;
  }
  
  TString summary(Form("numberoffiles=\"%d\" size=\"%7.2f GB\" ",count,size/byte2GB));
  
  out << Form("  <summary %s />",summary.Data()) << std::endl;
  out << "  </collection>" << std::endl;
  out << "</alien>" << std::endl;
  
  out.close();
  
  delete res;

  Bool_t ok = CopyFile(filename.Data());

  if (!ok) return 0;

  return count;
}

//_____________________________________________________________________________
void AliMuonQAMergeSubmitter::Print(Option_t*) const
{
  std::cout << "SplitMaxInputFileNumber : " << GetSplitMaxInputFileNumber() << std::endl;
  AliMuonGridSubmitter::Print();
}

//_____________________________________________________________________________
Bool_t AliMuonQAMergeSubmitter::Run(const char* mode)
{
  if (!IsValid()) return kFALSE;
  
  TString smode(mode);
  smode.ToUpper();
  
  if ( smode == "FULL")
  {
    return  ( Run("LOCAL") && Run("UPLOAD") && Run("SUBMIT") );
  }
  
  if ( smode == "LOCAL")
  {
    return CopyTemplateFilesToLocal();
  }
  
  if ( smode == "UPLOAD" )
  {
    return CopyLocalFilesToRemote();
  }
  
  if ( smode == "TEST" )
  {
    Bool_t ok = Run("LOCAL") && Run("UPLOAD");
    if ( ok )
    {
      ok = (Submit(kTRUE)>0);
    }
    return ok;
  }
  
  if ( smode == "SUBMIT" )
  {
    return (Submit(kFALSE)>0);
  }
  
  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliMuonQAMergeSubmitter::SetRemoteDir(const char* dir)
{
  if ( AliMuonGridSubmitter::SetRemoteDir(dir) )
  {
    Validate();
    return kTRUE;
  }
  else
  {
    Invalidate();
    return kFALSE;
  }
}

//______________________________________________________________________________
void AliMuonQAMergeSubmitter::ShowStage(Int_t runNumber)
{
  /// Show stage for a given run number
  
  Int_t stage(0);
  
  if ( RemoteFileExists(Form("%s/%d/%s",RemoteDir().Data(),runNumber,fWhatToMerge.Data())) )
  {
    stage = kFinal;
  }
  else
  {
    stage = GetLastStage(Form("%s/%d",RemoteDir().Data(),runNumber));
  }

  std::cout << "RUN " << runNumber << " ";

  if ( stage == kFinal )
  {
    std::cout << "FINAL";
  }
  else
  {
    std::cout << "Stage " << stage;
  }
  std::cout << std::endl;
}

//______________________________________________________________________________
void AliMuonQAMergeSubmitter::ShowStages()
{
  /// Show in remote dir the list of stages we're in for each run
  
  TGridResult* r = gGrid->Ls(RemoteDir());
  Int_t i(0);
  std::map<int,std::vector<int> > stages;
  
  while ( r->GetFileName(i) )
  {
    TString s(r->GetFileName(i));
    if  (s.IsDec())
    {
      Int_t runNumber = s.Atoi();
      Int_t stage(0);
      
      if ( RemoteFileExists(Form("%s/%d/%s",RemoteDir().Data(),runNumber,fWhatToMerge.Data())) )
      {
        stage = kFinal;
      }
      else
      {
        stage = GetLastStage(Form("%s/%d",RemoteDir().Data(),runNumber));
      }

      stages[stage].push_back(runNumber);
    }
    ++i;
  }
  delete r;
  
  std::map<int,std::vector<int> >::const_iterator it;
  
  for ( it = stages.begin(); it != stages.end(); ++it )
  {
    const std::vector<int>& runs = it->second;
    Int_t stage = it->first;
    if ( stage == kFinal )
    {
      std::cout << "FINAL";
    }
    else
    {
      std::cout << "Stage " << stage;
    }
    std::cout << std::endl;
    for ( std::vector<int>::size_type irun = 0; irun < runs.size(); ++irun )
    {
      std::cout << runs[irun] << " ";
    }
    std::cout << std::endl;
  }
}


//_____________________________________________________________________________
Bool_t AliMuonQAMergeSubmitter::Submit(Int_t runNumber, Bool_t dryRun)
{
  /// Submit merging job for one run
  
  TString runDir = GetRemoteDir(Form("%s/%d", RemoteDir().Data(), runNumber),kTRUE);
  
  if (RemoteFileExists(Form("%s/%s", runDir.Data(),fWhatToMerge.Data())))
  {
    AliWarning(" ! final merging already done");
    return kTRUE;
  }
  
  Int_t lastStage = GetLastStage(runDir.Data());
  
  AliDebug(1,Form("lastStage=%d",lastStage));
  
  ++lastStage;
  
  UInt_t n = MakeXMLCollectionForRun(runNumber,lastStage);
  
  Bool_t final  = ( n < GetSplitMaxInputFileNumber() );
  
  TString query;
  TString jdl(MergeJDLName(final));
  
  gGrid->Cd(RemoteDir().Data());
  
  query.Form("submit %s %d %d", jdl.Data(), runNumber, lastStage);

  AliInfo(Form("query=%s",query.Data()));
  
  if (dryRun)
  {
    return kTRUE;
  }

  Bool_t done = kFALSE;
  TGridResult *res = gGrid->Command(query);
  if (res)
  {
    TString cjobId1 = res->GetKey(0,"jobId");
    if (!cjobId1.IsDec())
    {
      AliError(" FAILED");
      gGrid->Stdout();
      gGrid->Stderr();
    }
    else
    {
      AliInfo(Form(" DONE\n   --> the job Id is: %s \n", cjobId1.Data()));
      done = kTRUE;
    }
    delete res;
  }
  else
  {
    AliError(" FAILED");
  }
  
  return done;
}

//_____________________________________________________________________________
Int_t AliMuonQAMergeSubmitter::Submit(Bool_t dryRun)
{
  /// Submit merging jobs

  if (!NofRuns())
  {
    AliError("No run to work with");
    return 0;
  }
  
  const std::vector<int>& runs = RunList();
  
  TString failed;
  
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i )
  {
    Int_t run = runs[i];
    Bool_t ok = Submit(run,dryRun);
  
    if (!ok)
    {
      failed += TString::Format("%d",run);
      failed += " ";
    }
  }
  
  if (failed.Length()>0)
  {
    AliInfo(Form("List of failed runs : %s",failed.Data()));
    return 0;
  }
  
  return 1;
}

