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
 * @param proof     Proof mode
 * @param mc        Monte-Carlo truth handler installed 
 *
 * @ingroup pwg2_forward_scripts
 */
void
Pass1(const char* esddir=".", 
      Int_t       nEvents=-1,
      Int_t       proof=0, 
      Bool_t      mc=false)
{
  Printf("  MC:            %s", mc        ? "yes" : "no");
  Printf("  Proof mode     %s", proof > 0 ? "yes" : "no");

  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/MakeAOD.C"); 

  MakeAOD(esddir, nEvents, proof, mc);
}
//
// EOF
//
