DCAmc_std(const Char_t *filename)
{

  gROOT->LoadMacro("MakeLibs.C");
  gROOT->LoadMacro("DCA.C");
  MakeLibs();

  /* start DCAmc */
  DCAmc(filename);

}
