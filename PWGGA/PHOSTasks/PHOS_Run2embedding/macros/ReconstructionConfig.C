/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Reconstruction steering script
 *
 */

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

enum EReconstruction_t {
  kReconstructionDefault,
  kReconstructionMuon,
  kReconstructionCustom,
  kNReconstructions
};

const Char_t *ReconstructionName[kNReconstructions] = {
  "Default",
  "Muon",
  "Custom"
};

/*****************************************************************/

ReconstructionConfig(AliReconstruction &rec, EReconstruction_t tag, Int_t run)
{

  switch(tag) {

    // Default
  case kReconstructionDefault:
    ReconstructionDefault(rec, run);
    rec.SetRunReconstruction("PHOS");
   return;
    
    // Muon
  case kReconstructionMuon:
    ReconstructionDefault(rec, run);
    rec.SetRunReconstruction("MUON ITS VZERO");
    return;
    
    // Custom
  case kReconstructionCustom:
    if ((gROOT->LoadMacro("ReconstructionCustom.C")) != 0) {
      printf("ERROR: cannot find ReconstructionCustom.C\n");
      abort();
      return;
    }
    ReconstructionCustom(rec, run);
    return;

  }  

}

ReconstructionDefault(AliReconstruction &rec, Int_t run)
{
    //
    //    // set OCDB snapshot mode
//        AliCDBManager *man = AliCDBManager::Instance();
//        man->SetDefaultStorage("alien://Folder=/alice/data/2015/OCDB");
//        man->SetRun(run);
//        man->SetSnapshotMode("OCDBrec.root");
     rec.SetCDBSnapshotMode("Sim/OCDBrec.root");
    //
//     rec.SetCleanESD(kFALSE);
//     rec.SetStopOnError(kFALSE);
//     rec.SetWriteESDfriend();
//     rec.SetWriteAlignmentData();
//     rec.SetFractionFriends(.1);
//     rec.SetRunPlaneEff(kTRUE);
//     rec.SetUseTrackingErrorsForAlignment("ITS");
    //
    rec.SetRunQA(":");
    //
}
