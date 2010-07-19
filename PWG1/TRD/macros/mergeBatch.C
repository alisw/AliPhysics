//______________________________________________________
Char_t* mergeBatch(const Char_t *mark, const Char_t *files, const Int_t nfiles=20, const Int_t first=0, Bool_t kSVN=kTRUE, Bool_t kCLEAR=kFALSE)
{
// Merge files specified in the file list "files" by the token "mark".
// The script will merge "nfiles" files starting from the "first" file. 
// If the file "svnInfo.log" is found together with the files to be merged it is copied locally 
// if option "kSVN". The input files are removed from disk if option "kCLEAR".
//
// On return the name of the merged file is return or NULL in case of failure.
//
  TObjArray arr(nfiles); arr.SetOwner(kTRUE);

  TFileMerger fFM(kTRUE);
  fFM.OutputFile(Form("%s/%d_%s",  gSystem->ExpandPathName("$PWD"), first, mark));
  Int_t iline(0), nbatch(0);
  std::string filename;
  std::ifstream file(files);
  while(getline(file, filename)){
    if(Int_t(filename.find(mark)) < 0) continue;
    iline++;
    if(iline<first) continue;
    if(kSVN){ // download SVN info for trending
      string base=filename.substr(0, filename.find_last_of('/'));
      if(gSystem->Exec(Form("if [ ! -f svnInfo.log ]; then cp -v %s/svnInfo.log %s; fi", base.c_str(), gSystem->ExpandPathName("$PWD"))) == 0) kSVN=kFALSE;
    }
    fFM.AddFile(filename.c_str()); 
    arr.Add(new TObjString(filename.c_str()));
    nbatch++;
    if(nbatch==nfiles) break;
  }
  if(!nbatch){
    Info("mergeBatch.C", "NOTHING TO MERGE"); return NULL;
  }else if(nbatch<nfiles){
    Info("mergeBatch.C", "MERGING INCOMPLETE BATCH %d [%d] ... ", nbatch, nfiles);
  } else {
    Info("mergeBatch.C", "MERGING BATCH STARTING @ %d ... ", first);
  }
  fFM.Merge();

  if(kCLEAR){
    for(Int_t ifile(0); ifile<arr.GetEntries(); ifile++){
      gSystem->Exec(Form("rm -fv %s", ((TObjString*)arr.At(ifile))->GetString().Data()));
    }
  }
  return fFM.GetOutputFileName(); 
}
