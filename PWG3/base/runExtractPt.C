void runExtractPt(TH1F* histPtDaug, Bool_t readKineFromNtupla=kTRUE, const char *pathFileName="galice.root") 
 {
  /* Run-macro to extract pt-spectra (and ptMin-spectra) for mothers particles 
     input: 1) pt histogram of daughter particles 
            2) boolean flag: kFALSE -> read Kinematics.root to evaluate correction factors, 
                                       create a TNtuple with kinematic informations of mothers and 
                                       daughters and store it in the file "DecayKine.root"
                             kTRUE  -> read the TNtupla from the file "DecayKine.root" (after it is 
                                       created) to evaluate correction factors  
            3) path of the file "galice.root" to read the Kinematics.root (not needed after the 
               TNtupla is created)
               
     output: file Mothers.root which contains pt-spectra and ptMin-spectra of mothers particles    */
  
  gSystem->Load("libPWG3base.so");

  AliPtMothFromPtDaugh *ptExtr = new AliPtMothFromPtDaugh();
  ptExtr->SetDefaultAnalysis(AliPtMothFromPtDaugh::kBtoJPSI);
  ptExtr->SetBinsPtMoth(0.,10,20,1);
  ptExtr->SetBinsPtMinMoth(0.,10,20,1);
  ptExtr->SetEtaMothers(-1.5,1.5);
  ptExtr->SetEtaDaughter(-1.,1.);
  if(!ptExtr->ReadHistoPtDaught(histPtDaug)) 
     { printf("Daughter pt-Histogram is not defined \n"); return; }
  if(!ptExtr->CreateWeights()) return; 
  ptExtr->SetReadKineFromNtupla(readKineFromNtupla);
  ptExtr->ReadKinematics(pathFileName);
  ptExtr->EvaluatePtMoth();
  ptExtr->WritePtMothHistoToFile();
  return;
 }