void MakeMatBudPlots(int layer=0)
{
  // Simple interface to GetMaterialBudget.C macro to plot the
  // material budget (for each material separately) for a given layer
  // M.S. 01 Jul 2014

  const int kNLr = 7;
  const int kNLrInner = 3;

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gROOT->LoadMacro("GetMaterialBudget.C");

  enum {kRmn,kRmd,kRmx,kNModPerStave,kPhi0,kNStave,kNPar};
  // Radii are from last TDR (ALICE-TDR-017.pdf Tab. 1.1, rMid is mean value)
  // TO BE KEPT IN SYNC WITH CreateITSUv2.C MACRO!!!
  const double tdr5dat[kNLr][kNPar] = { 
    {2.24, 2.34, 2.67,  9., 16.42, 12}, // for each inner layer: rMin,rMid,rMax,NChip/Stave, phi0, nStaves
    {3.01, 3.15, 3.46,  9., 12.18, 16},
    {3.78, 3.93, 4.21,  9.,  9.55, 20},
    {-1,  19.6 ,   -1,  4.,  0.  , 24},  // for others: -, rMid, -, NMod/HStave, phi0, nStaves
    {-1,  24.55, -1,    4.,  0.  , 30},
    {-1,  34.39, -1,    7.,  0.  , 42},
    {-1,  39.34, -1,    7.,  0.  , 48} 
  };

  if (layer < 0 || layer >= kNLr) {
    printf("Wrong layer number %d - giving up\n",layer);
    return;
  }

  double rmin, rmax;

  if (layer < kNLrInner) { // Inner layers
    rmin = tdr5dat[layer][kRmn] - 0.1;
    rmax = tdr5dat[layer][kRmx] + 0.1;
  } else { // Outer layers
    rmin = tdr5dat[layer][kRmd] - 0.175;
    rmax = tdr5dat[layer][kRmd] + 0.35;
  }

  printf("Drawing material budget for layer %d (Rmin = %f Rmax = %f\n",layer,
	 rmin,rmax);
  DrawMaterialBudget_Splitted(layer, rmin, rmax) ;

}
