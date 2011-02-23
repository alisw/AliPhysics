aodv1(Bool_t run=0)
{ 
  gROOT->LoadMacro("aodtools.C+g");

  AODProdParameters params;

  LoadLibraries(params);
  gROOT->LoadMacro("runCaloAOD.C+g");
  if (run)
    runCaloAOD(params);
}
