void PathLengthTest() {
  //
  // It's a test macro for path length calculation in AliFastGlauber
  //
  // Andrea Dainese
  //

  AliFastGlauber g;

  g.Init(2);

  // centrality (fraction of geometrical cross section
  g.SetCentralityClass(0.00,0.10); 

  // plot b,ell, and ell1 vs ell2 (back-to-back partons) in that centrality 
  // class
  Int_t nEntries = 100;
  g.PlotBDistr(nEntries);
  g.PlotLengthDistr(nEntries);
  g.PlotLengthB2BDistr(nEntries);


  // examples on getting stuff out of it
  Double_t b,ell,ell1,ell2;

  g.GetRandomBHard(b);
  printf(" Random b in cetrality class (according to hard cross section):\n %f fm\n",b);

  printf(" Length:\n");
  g.GetLength(ell,b);
  printf("   for this b: %f fm\n",ell);
  g.GetLength(ell);
  printf("   for random b : %f fm\n",ell);

  printf(" Lengths for two partons back-to-back:\n");
  g.GetLengthsBackToBack(b,ell1,ell2);
  printf("   for this b: %f fm and %f fm\n",ell1,ell2);
  g.GetLengthsBackToBack(ell1,ell2);
  printf("   for random b: %f fm and %f fm\n",ell1,ell2);


  Double_t phis[3]={2.,4.,6.};
  Double_t ells[3];
  printf(" Lengths for N partons from PYTHIA:\n");
  g.GetLengthsForPythia(3,phis,ells,b);
  printf("   for this b: %f fm, %f fm and %f fm\n",ells[0],ells[1],ells[2]);
  g.GetLengthsForPythia(3,phis,ells);
  printf("   for random b: %f fm, %f fm and %f fm\n",ells[0],ells[1],ells[2]);

  return;
}
