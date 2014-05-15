{
  gSystem->Load("libGui");
#if 0
  const char* imgs[] = { 
    "arrow_down.xpm",
    "unchecked_t.xpm",
    0 
  };
  const char** pImg = imgs;
  while (*pImg) { 
    Info("", "Client %p loading %s", gClient, *pImg);
    if (!gClient->GetPicture(*pImg)) { 
      Warning("", "Client %p failed to load %s", gClient, *pImg);
    }
    pImg++;
  }
#endif

  gROOT->LoadMacro("$ALICE_ROOT/../master-src/FMD/scripts/Compile.C");
  Compile("$ALICE_ROOT/../master-src/FMD/scripts/DrawCalib.C","g");
  // gROOT->LoadMacro("DrawCalib.C");
  // Info("", "Client @ %p", gClient);

  DrawCalib();
}
