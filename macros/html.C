void html() {
// to run this macro, you must have the correct .rootrc file
// in your galice directory.
// The gAlice classes summary documentation go to directory html
// The gAlice classes source  documentation go to directory html/src
// The example macros documentation go to directory html/examples
   
   gROOT->LoadMacro("loadlibs.C");
   loadlibs();
   THtml html;
   html.MakeAll(1);

   // process some example macros
   html.Convert("aliroot.cxx","The gAlice main program");
   gROOT->LoadMacro("display.C");
   html.Convert("display.C","Macro to start the event display");
   gROOT->LoadMacro("ecut.C");
   html.Convert("ecut.C","Macro to animate the rapidity cut slider");
   gROOT->LoadMacro("HMPID.C");
   html.Convert("HMPID.C","Macro to read HMPID events");
   gROOT->LoadMacro("TPCHits2Digits.C");
   html.Convert("TPCHits2Digits.C","Macro to convert TPC hits to Digits");
   gROOT->LoadMacro("TPCHits2Clusters.C");
   html.Convert("TPCHits2Clusters.C","Macro to convert TPC hits to Clusters");
   gROOT->LoadMacro("TPCDigitsDisplay.C");
   html.Convert("TPCDigitsDisplay.C","Macro to display TPC digits");
   gROOT->LoadMacro("TPCDigits2Clusters.C");
   html.Convert("TPCDigits2Clusters.C","Macro to find TPC clusters");
   gROOT->LoadMacro("anal.C");
   html.Convert("anal.C","Macro to analyse events");
   gROOT->LoadMacro("newanal.C");
   html.Convert("newanal.C","Macro to analyse events");
   gROOT->LoadMacro("PHOS.C");
   html.Convert("PHOS.C","Macro to reconstruct PHOS");
   gROOT->LoadMacro("Config.C");
   html.Convert("Config.C","AliRoot configuration file");
   gROOT->LoadMacro("menu.C");
   html.Convert("menu.C","menu configuration file");
   gROOT->LoadMacro("lego.C");
   html.Convert("lego.C","visualisation of the lego results");
//   gROOT->LoadMacro("html.C");
   html.Convert("html.C","generation of the html doc");
   gROOT->LoadMacro("grun.C");
   html.Convert("grun.C","simple batch macro");
}
