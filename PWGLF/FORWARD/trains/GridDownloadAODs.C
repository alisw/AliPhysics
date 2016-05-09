/**
 * @file   GridDownloadAODs.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue May  3 18:27:13 2016
 * 
 * @brief  Script to download AODs 
 * 
 * 
 */

void
GridDownloadAODs(const char* url,
		 Bool_t      verbose=false,
		 Bool_t      force=false)
{
  if (!TGrid::Connect("alien://")) {
    Error("GridDownloadAODs","Failed to connect to AliEn");
    return;
  }
  Printf("Directory: %s", url);
  Printf("Verbosity: %d", verbose);
  Printf("Force:     %d", force);
	 

  TString      dir(url);
  TString      pat("*/AliAOD.root");
  TGridResult* r = gGrid->Query(dir,pat);
  if (!r) {
    Error("GridDownloadAODs","No result from query");
    return;
  }

  TStopwatch timer;
  timer.Reset();
  Int_t c    = 0;
  Int_t n    = r->GetEntries();
  Int_t save = gErrorIgnoreLevel;
  if (verbose) Printf("=== Got a total of %d AOD files",n);
  for (Int_t i = 0; i < n; i++) {
     TString path(r->GetKey(i, "turl"));
     TString dir(gSystem->DirName(path));
     TString sub(gSystem->BaseName(dir));
     TString subsub(gSystem->BaseName(gSystem->DirName(dir)));
     TString out = TString::Format("AliAOD_%s_%s.root",
                                   subsub.Data(),sub.Data());
     if (!gSystem->AccessPathName(out.Data()) && !force) {
       if (verbose) Printf("=== Already have %s",out.Data());
       continue;
     }
     if (verbose) Printf("=== Getting %s %s (%3d/%3d)",
			 subsub.Data(),sub.Data(),i,n);
     else {
       printf("%5d/%5d %22s ", i+1, n, out.Data());
       fflush(stdout);
     }
     
     gErrorIgnoreLevel = kError;
     if (!TFile::Cp(path, out, verbose)) {
       gErrorIgnoreLevel = save;
       Warning("DownloadAODs","Failed to copy %s -> %s",
               path.Data(), out.Data());
       continue;
     }
     c++;
     gErrorIgnoreLevel = save;
     if (verbose) continue;
     
     Double_t passed = timer.RealTime();
     timer.Continue();
     Double_t perFile = passed / c;
     Double_t remain  = (n-i-1) * perFile;
     printf("%5.1f%%  [%3d copied] ETA: %4d:%02d:%02d (%5.1fs/file)\r",
	    100*float(i+1)/n, c, 
	    Int_t(remain/60/60),
	    Int_t(remain/60)%60,
	    Int_t(remain)%60,
	    perFile);
     fflush(stdout);
	    
   }

}
//
// EOF
// 
