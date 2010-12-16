/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Run first pass of the analysis - that is read in ESD and produce AOD
 * 
 * @param esddir    ESD input directory. Any file matching the pattern 
 *                  *AliESDs*.root are added to the chain 
 * @param nEvents   Number of events to process.  If 0 or less, then 
 *                  all events are analysed
 * @param flags     Job flags. A bit mask of 
 *  - 0x01 (MC)        Monte-Carlo truth handler installed 
 *  - 0x02 (PROOF)     Proof mode 
 *  - 0x04 (FULL)      Run full analysis - including terminate 
 *  - 0x08 (ANALYSE)   Run only analysis - not terminate 
 *  - 0x10 (TERMINATE) Run no analysis - just terminate.  
 * 
 * of these, PROOF, FULL, ANALYSE, and TERMINATE are mutually exclusive. 
 *
 * If PROOF mode is selected, then Terminate will be run on the master node 
 * in any case. 
 * 
 * If FULL is selected, then the full analysis is done.  Note that
 * this can be combined with PROOF but with no effect.
 *
 * ANALYSE cannot be combined with FULL, PROOF, or TERMINATE.  In a
 * local job, the output AnalysisResults.root will still be made, but
 * terminate is not called.
 *
 * In TERMINATE, the file AnalysisResults.root is opened and all
 * containers found are connected to the tasks.  The terminate member
 * function is then called
 * 
 *
 * @ingroup pwg2_forward_scripts
 */
void
Pass1(const char* esddir=".", 
      Int_t       nEvents=1000,
      UShort_t    flags=0x4)
{
  Printf("Flags: 0x%04x", flags);
  Printf("  MC:            %s", flags & 0x01 ? "yes" : "no");
  Printf("  Proof mode     %s", flags & 0x02 ? "yes" : "no");
  Printf("  Full analysis  %s", flags & 0x04 ? "yes" : "no");
  Printf("  Analyse only   %s", flags & 0x08 ? "yes" : "no");
  Printf("  Terminate only %s", flags & 0x10 ? "yes" : "no");

  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/RunManager.C"); 

  RunManager(esddir, nEvents, flags);
}
//
// EOF
//
