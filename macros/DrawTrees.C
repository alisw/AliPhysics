void DrawTrees()
{
   TControlBar *menu = new TControlBar("vertical","DrawTrees menu");
   menu->AddButton("TPC tree",     "((TGeant3*)gMC)->Gdtree(\"tpc\")","Shows the Geant tree for the TPC");
   menu->AddButton("ITS tree",     "((TGeant3*)gMC)->Gdtree(\"ITSV\")","Shows the Geant tree for the ITS");
   menu->AddButton("CASTOR tree",  "((TGeant3*)gMC)->Gdtree(\"OCTA\")","Shows the Geant tree for the CASTOR");
   menu->AddButton("ABSO tree",    "((TGeant3*)gMC)->Gdtree(\"ABSM\")","Shows the Geant tree for the ABSO");
   menu->AddButton("DIPO tree",    "((TGeant3*)gMC)->Gdtree(\"DDIP\")","Shows the Geant tree for the DIPO");
   menu->AddButton("FMD tree",     "((TGeant3*)gMC)->Gdtree(\"IWR3\")","Shows the Geant tree for the FMD");
   menu->AddButton("FRAME tree",   "((TGeant3*)gMC)->Gdtree(\"BFMO\")","Shows the Geant tree for the FRAME");
   menu->AddButton("HALL tree",    "((TGeant3*)gMC)->Gdtree(\"HUWA\")","Shows the Geant tree for the HALL");
   menu->AddButton("MAG tree",     "((TGeant3*)gMC)->Gdtree(\"L3MO\")","Shows the Geant tree for the MAG");
   menu->AddButton("MUON tree",    "((TGeant3*)gMC)->Gdtree(\"CH1A\")","Shows the Geant tree for the MUON");
   menu->AddButton("PHOS tree",    "((TGeant3*)gMC)->Gdtree(\"phos\")","Shows the Geant tree for the PHOS");
   menu->AddButton("PIPE tree",    "((TGeant3*)gMC)->Gdtree(\"QQMO\")","Shows the Geant tree for the PIPE");
   menu->AddButton("PMD tree",     "((TGeant3*)gMC)->Gdtree(\"DPMD\")","Shows the Geant tree for the PMD");
   menu->AddButton("HMPID tree",    "((TGeant3*)gMC)->Gdtree(\"rich\")","Shows the Geant tree for the HMPID");
   menu->AddButton("SHIL tree",    "((TGeant3*)gMC)->Gdtree(\"YMOT\")","Shows the Geant tree for the SHIL");
   menu->AddButton("START tree",    "((TGeant3*)gMC)->Gdtree(\"T0ST\")","Shows the Geant tree for the START");
   menu->AddButton("TOF tree",     "((TGeant3*)gMC)->Gdtree(\"FBAR\")","Shows the Geant tree for the TOF");
   menu->AddButton("TRD tree",     "((TGeant3*)gMC)->Gdtree(\"trd\")","Shows the Geant tree for the TRD");
   menu->AddButton("ZDC tree",     "((TGeant3*)gMC)->Gdtree(\"zdc\")","Shows the Geant tree for the ZDC");
   menu->AddButton("V0R tree", "((TGeant3*)gMC)->Gdtree(\"V0RI\")","Shows the Geant tree for the V0R");
   menu->AddButton("V0L tree", "((TGeant3*)gMC)->Gdtree(\"V0LE\")","Shows the Geant tree for the V0L");
   menu->Show();
}



