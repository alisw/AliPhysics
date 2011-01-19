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
Pass2(const char* file=".", 
      const char* triggers="INEL", 
      Int_t       energy=900, 
      Double_t    vzMin=-10, 
      Double_t    vzMax=10, 
      Int_t       rebin=5, 
      const char* title="",
      bool        hhd=true,
      bool        comp=true)
{
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/Compile.C"); 
  Compile("$ALICE_ROOT/PWG2/FORWARD/analysis2/AnalyseAOD.C","g"); 
  
  Int_t trgMask; 
  TString     trgs(triggers);
  trgs.ToUpper();
  TObjString* trg;
  TIter       next(trgs.Tokenize(" ,|"));
  while ((trg = static_cast<TObjString*>(next()))) { 
    TString s(trg->GetString());
    if      (s.IsNull()) continue;
    if      (s.CompareTo("INEL")  == 0) trgMask = AliAODForwardMult::kInel;
    else if (s.CompareTo("INEL>0")== 0) trgMask = AliAODForwardMult::kInelGt0;
    else if (s.CompareTo("NSD")   == 0) trgMask = AliAODForwardMult::kNSD;
    else 
      Warning("Pass2", "Unknown trigger %s", s.Data());
  }
  if (trgMask == 0) {
    trgMask = 1;
    trgs.Append("INEL");
  }
  TString tit(title);
  tit.ReplaceAll("@", " ");

  printf("--------------------------------------\n"
         "Settings for this:\n"
         "  Input AOD:    %s\n" 
         "  Vertex range: %+4.1f -> %+4.1f cm\n" 
	 "  Rebinning:    %d\n"
	 "  Trigger mask: 0x%02x (%s)\n"
	 "  Energy:       %dGeV\n"
	 "  Title:        %s\n"
	 "  HHD comp.:    %s\n"
	 "  Other comp.:  %s\n"
         "--------------------------------------\n",
	 file, vzMin, vzMax, rebin, trgMask, trgs.Data(), energy, tit.Data(),
	 hhd ? "yes" : "no", comp ? "yes" : "no");
  
  AnalyseAOD dr;
  TStopwatch t;
  t.Start();
  dr.Run(file, vzMin, vzMax, rebin, trgMask, energy, tit.Data(), hhd, comp);
  t.Stop();
  t.Print();  
}
//
// EOF
//

      
