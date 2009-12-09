
//
// AliEn Multi MasterJob Coordination
// Author: K. Read
//

void mylauncher()
{
  TString worksubdir = "work12";  //already existing work directory
//TString worksubdir = "LHC09b4"; //a better name for MultiMode
  // Name of JDL file to upload.  Leave blank for no upload and no submit.
  TString jdlfilename = "anaJete.jdl";//ignored for MultiMode
  // Name of executable to upload.  Leave blank for no upload.
  TString collectfilename = "mycollect.xml";//ignored for MultiMode
  // List any other files to upload in string filenames separated by blanks.
  TString execfilename = "";
  // List any other files to upload in string filenames separated by blanks.
  TString filelist = "run.C validate.sh anaJete.C ConfigJetAnalysisFastJet.C ConfigAnalysisElectron.C mergeout.jdl mergeoutscaled.jdl ANALYSIS.par ANALYSISalice.par AOD.par EMCALUtils.par ESD.par PHOSUtils.par STEERBase.par JETAN.par FASTJETAN.par PWG4PartCorrBase.par PWG4PartCorrDep.par";
  TString filename;

  Bool_t kMultiMode = kFALSE;//uncomment for regular mode
//Bool_t kMultiMode = kTRUE; //uncomment for coordinated multiple masterjobs
//MultiMode Variables.   All ignored for regular mode.
  TString multijdlfilename = "anaJetemulti.jdl";
//TString datatype = "AliESDs.root";//uncomment for ESDs
  TString datatype =  "AliAOD.root";//uncomment for AODs
  TString findlist = "
/alice/sim/LHC09b4/AOD/000
/alice/sim/LHC09b4/AOD/001
/alice/sim/LHC09b4/AOD/002"; //line-feed separated list of desired collections


  gSystem->Load("libNetx.so") ;
  gSystem->Load("libRAliEn.so");

  TGrid::Connect("alien://") ;
  if (gGrid && gGrid->IsConnected()) {
    TString homedir = gGrid->GetHomeDirectory(); // has a trailing slash
    TString workdir = homedir + worksubdir;
    if (gGrid->Cd(workdir)) {

      // Upload files listed
      if (filelist.Length()) {
        arr = filelist.Tokenize(" ");
        TObjString *os;
        TIter next(arr);
        while ((os=(TObjString*)next())) {
          Info("Grid Upload", "Copying %s to your AliEn work directory", os->GetString().Data());
	  if (FileExists(os->GetString())) gGrid->Rm(os->GetString());
	  TFile::Cp(Form("file:%s",os->GetString().Data()), Form("alien://%s/%s", workdir.Data(), os->GetString().Data()));
	}   
	delete arr;   
      }

      // Prepare collection(s)
      if(kMultiMode){
        // Find collections if desired
        if (findlist.Length()) {
          arr = findlist.Tokenize("\n");
          TObjString *os;
          TIter next(arr);
          Int_t count=0;
          while ((os=(TObjString*)next())) {
	    TGridResult *res;
            Info("Find", "Making collection for %s", os->GetString().Data());
            res = gGrid->Command(Form("find -x collection %s %s > mycollect%02d.xml", os->GetString().Data(),datatype.Data(),count));
            if(!res) Info("Find", "Failed making collection for %s \n",os->GetString().Data());
            delete res;
            count++;
          }
        }
      }
      else{
	if (collectfilename.Length()) {
          Info("Grid Upload", "Copying %s to your AliEn work directory", collectfilename.Data());
	  if (FileExists(collectfilename)) gGrid->Rm(collectfilename);
	  TFile::Cp(Form("file:%s",collectfilename.Data()), Form("alien://%s/%s", workdir.Data(), collectfilename.Data()));
        }
      }

      // Upload executable if listed
      if (execfilename.Length()) {
        filename = Form("%sbin/%s", homedir.Data(), execfilename.Data());
        if (FileExists(filename)) gGrid->Rm(filename);
        Info("Grid Upload", "Copying executable file %s to your AliEn bin directory", execfilename.Data());
        TFile::Cp(Form("file:%s",execfilename.Data()), Form("alien://%s", filename.Data()));
      }

      // Upload and submit JDL if listed
      if (kMultiMode) jdlfilename = multijdlfilename;
      if (jdlfilename.Length()) {
        filename = Form("%s/%s", workdir.Data(), jdlfilename.Data());
        if (FileExists(filename)) gGrid->Rm(filename);
        Info("Grid Upload", "Copying JDL file %s to your AliEn work directory", jdlfilename.Data());
        TFile::Cp(Form("file:%s",jdlfilename.Data()), Form("alien://%s", filename.Data()));

        if (!kMultiMode) Int_t count=1;
        for (Int_t index = 0; index < count; index++) {
	  TGridResult *res;
          TString jobID = "";
          if(kMultiMode){
            res = gGrid->Command(Form("submit %s %s %02d", jdlfilename.Data(),worksubdir.Data(),index));
            Info("Launcher:",     "Submitting %s %s %02d", jdlfilename.Data(),worksubdir.Data(),index);
          }
          else{
            res = gGrid->Command(Form("submit %s", jdlfilename.Data()));
            Info("Launcher:",     "Submitting %s", jdlfilename.Data());
          }
          if (res) {
            const char *cjobId = res->GetKey(0,"jobId");
            if (!cjobId) {
              Error("Launcher:", "Your JDL %s could not be submitted", jdlfilename.Data());
              return;
            } 
            else {
              Info("Launcher:", "Your JDL %s was successfully submitted.\n\n\t\t\t THE JOB ID IS: %s\n",
		   jdlfilename.Data(), cjobId);
            }          
            delete res;
          }
        }
      }

      // Launch alien shell
      gSystem->Exec("aliensh");

    }
  }
}

Bool_t FileExists(const char *lfn) const
{
// Returns true if file exists.
   if (!gGrid) {
      Error("FileExists", "No connection to grid");
      return kFALSE;
   }
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
