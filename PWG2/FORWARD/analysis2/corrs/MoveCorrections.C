/** 
 * 
 * 
 * @param fname 
 * @param n 
 * 
 * @return 
 * @ingroup pwg2_forward_scripts_corr
 */
Bool_t
Backup(const TString& fname, Int_t n=10)
{
  TString fn(fname); fn.Append(Form(".%d", n));
  if (!gSystem->AccessPathName(fn.Data())) {
    Error("Backup", "Last possible backup slot (%d) filled for %s", 
	  n, fname.Data());
    return false;
  }

  for (Int_t i=n-1; i >= 0; i--) { 
    TString fi(fname); 
    if (i > 0) fi.Append(Form(".%d", i));
    if (gSystem->AccessPathName(fi.Data())) continue;

    TString fb(fname); fb.Append(Form(".%d", i+1));

    if (gSystem->Rename(fi.Data(), fb.Data())) { 
      Error("Backup", "Failed to backup %s to %s", fi.Data(), fb.Data());
      return false;
    }
  }
  return true;
}


/** 
 * 
 * 
 * @param what 
 * 
 * @return 
 * @ingroup pwg2_forward_scripts_corr
 */
Bool_t
MoveWhat(UInt_t what) 
{
  TString nWhat;
  switch (what) {
  case AliForwardCorrectionManager::kSecondaryMap:    
    nWhat = "secondary map";    break;
  case AliForwardCorrectionManager::kDoubleHit:	      
    nWhat = "double hit";       break;
  case AliForwardCorrectionManager::kVertexBias:      
    nWhat = "vertex bias";      break;
  case AliForwardCorrectionManager::kMergingEfficiency:
    nWhat = "merging efficiency";break;
  case AliForwardCorrectionManager::kELossFits:       
    nWhat = "energy loss fits";  break;
  }
  Info("MakeWhat", " Copying %s corrections", nWhat.Data());

  AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
  TString dir = gSystem->ExpandPathName(mgr.GetFileDir(what));

  // Make the directory if it doesn't exist
  if (gSystem->AccessPathName(dir.Data())) {
    Info("MakeWhat", " Making directory %s ... ", dir.Data());
    if (gSystem->mkdir(dir.Data(), true)) {
      Error("MakeWhat", "couldn't make directory %s", dir.Data());
      return false;
    }
  }

  TString pattern = mgr.GetFileName(what, 1, 900, 5, false);
  pattern.ReplaceAll("0900GeV", "[0-9]+GeV");
  pattern.ReplaceAll("pp", "[a-zA-Z]+");
  pattern.ReplaceAll("p5kG", "[pm][0-9]+kG");
  pattern.ReplaceAll("real", "[a-zA-Z]+");
  
  TPRegexp regex(pattern);
 
  TSystemDirectory sysdir(".", ".");
  TIter            next(sysdir.GetListOfFiles());
  TSystemFile*     file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) {
    if (file->IsDirectory()) continue;
    TString fname(file->GetName());
 
    if (!regex.Match(fname)) continue;

    // Info("MakeWhat", "  match: %s", fname.Data());
    TString to(gSystem->ConcatFileName(dir.Data(), fname.Data()));
    
    if (!Backup(to)) return false;

    Info("MakeWhat", "  copying %s\n               to %s", 
         fname.Data(), to.Data());
    if (gSystem->CopyFile(fname.Data(), to.Data(), false)) {
      Error("MakeWhat", "Failed to copy %s to %s", fname.Data(), to.Data());
      return false;
    }
  }
  return true;
}
/** 
 * 
 * 
 * @param sec 
 * @param dbl 
 * @param vtx 
 * @param merge 
 * @param eloss 
 * @ingroup pwg2_forward_scripts_corr
 */
void
MoveCorrections(bool sec=true, 
		bool dbl=true, 
		bool vtx=true,
		bool merge=true,
		bool eloss=true)
{
  Info("MoveCorrections", "Loadingl libraries ...");
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  Info("MoveCorrections", "Moving selected corrections ...");
  if (sec)   MoveWhat(AliForwardCorrectionManager::kSecondaryMap);
  if (dbl)   MoveWhat(AliForwardCorrectionManager::kDoubleHit);
  if (vtx)   MoveWhat(AliForwardCorrectionManager::kVertexBias);
  if (merge) MoveWhat(AliForwardCorrectionManager::kMergingEfficiency);
  if (eloss) MoveWhat(AliForwardCorrectionManager::kELossFits);       

}
  
