/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Read in AOD and generate @f$ dN/d\eta@f$ for the selected 
 * trigger classes and vertex ranges 
 * 
 * @param file     Input file (AOD)
 * @param triggers Triggers to investigate 
 * @param energy   Energy (only used for comparisons)
 * @param vzMin    Minimum interaction point z coordinate
 * @param vzMax    Maximum interaction point z coordinate
 * @param rebin    How many bins to group
 * @param title    Title to put on the plot 
 * @param hhd      Whether to do HHD comparison
 * @param comp     Whether to do comparisons 
 *
 * @ingroup pwg2_forward_scripts
 */
void
Pass2(const char* aoddir=".", 
      Int_t       nEvents=-1,
      const char* triggers="INEL", 
      Double_t    vzMin=-10, 
      Double_t    vzMax=10,
      Int_t       proof=0)
{
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/MakedNdeta.C"); 

  MakedNdeta(aoddir, nEvents, triggers, vzMin, vzMax, proof);
}
//
// EOF
//

      
