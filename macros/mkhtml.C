void mkhtml (char *macro=0, Int_t force=0) {
// to run this macro, you must have the correct .rootrc file
// in your galice directory.
// The gAlice classes summary documentation go to directory html
// The gAlice classes source  documentation go to directory html/src
// The example macros documentation go to directory html/examples
   
  // gROOT->LoadMacro("loadlibs.C");
  // loadlibs();
  THtml html;
  if(macro) {
    gROOT->LoadMacro(macro);
    html.Convert(macro,"Example Macro");
  } else {
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libRALICE");
    html.MakeAll(force);
  }
}
