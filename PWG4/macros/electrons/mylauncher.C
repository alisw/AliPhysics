
//
// Copy files to AliEn space and submit AliEn job
// Author: K. Read
//

void mylauncher()
{
  TString worksubdir = "work12";
  // Name of JDL file to upload.  Leave blank for no upload and no submit.
  TString jdlfilename = "anaJete.jdl";
  // Name of executable to upload.  Leave blank for no upload.
  TString execfilename = "";
  // List any other files to upload in string filenames separated by blanks.
  TString filelist = "run.C validate.sh anaJete.C ConfigJetAnalysisFastJet.C ConfigAnalysisElectron.C mergeout.jdl mergeoutscaled.jdl mycollect.xml ANALYSIS.par ANALYSISalice.par AOD.par ESD.par STEERBase.par JETAN.par FASTJETAN.par PWG4PartCorrBase.par PWG4PartCorrDep.par";
  TString filename;

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

      // Upload executable if listed
      if (execfilename.Length()) {
        filename = Form("%sbin/%s", homedir.Data(), execfilename.Data());
        if (FileExists(filename)) gGrid->Rm(filename);
        Info("Grid Upload", "Copying executable file %s to your AliEn bin directory", execfilename.Data());
        TFile::Cp(Form("file:%s",execfilename.Data()), Form("alien://%s", filename.Data()));
      }

      // Upload and submit JDL if listed
      if (jdlfilename.Length()) {
        filename = Form("%s/%s", workdir.Data(), jdlfilename.Data());
        if (FileExists(filename)) gGrid->Rm(filename);
        Info("Grid Upload", "Copying JDL file %s to your AliEn work directory", jdlfilename.Data());
        TFile::Cp(Form("file:%s",jdlfilename.Data()), Form("alien://%s", filename.Data()));

	TGridResult *res;
        TString jobID = "";
        res = gGrid->Command(Form("submit %s", jdlfilename.Data()));
        Info("Launcher:", "Submitting %s ", jdlfilename.Data());
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
