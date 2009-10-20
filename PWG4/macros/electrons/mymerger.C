
//
// AliEn Initial PYTHIA pt hard bin-by-bin merging to prepare for rescaling
// Author: K. Read
//

void mymerger()
{
  TString filename;
  TString jdlfilename = "mergeoutscaledi.jdl";
  //TString worksubdir = "LHC09b4b";
  //Int_t maxbin = 15; 
  TString worksubdir = "LHC09b2ESDb";
  Int_t maxbin = 16;

  //Later do: aliensh; cd worksubdir/output/merged; cp histoss* file:
  //And then process files locally with MergeFileInBins.C.

  gSystem->Load("libNetx.so") ;
  gSystem->Load("libRAliEn.so");

  TGrid::Connect("alien://") ;
  if (gGrid && gGrid->IsConnected()) {
    TString homedir = gGrid->GetHomeDirectory(); // has a trailing slash
    TString workdir = homedir + worksubdir;
    if (gGrid->Cd(workdir)) {

      // Upload and submit JDL if listed
      if (jdlfilename.Length()) {
        filename = Form("%s/%s", workdir.Data(), jdlfilename.Data());
        if (FileExists(filename)) gGrid->Rm(filename);
        Info("Grid Upload", "Copying JDL file %s to your AliEn work directory", jdlfilename.Data());
        TFile::Cp(Form("file:%s",jdlfilename.Data()), Form("alien://%s", filename.Data()));

        for (Int_t index = 0; index <= maxbin; index++) {
	  TGridResult *res;
          TString jobID = "";
          res = gGrid->Command(Form("submit %s %s %02d", jdlfilename.Data(),worksubdir.Data(),index));
          Info("Launcher:",     "Submitting %s %s %02d", jdlfilename.Data(),worksubdir.Data(),index);
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
