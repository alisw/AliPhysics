void readJets( Char_t* fileName = "./analysis/EOR_analyze_100000_kPythia6Jets104_125.root" ) {

  gSystem->Load("libCGAL.so");
  gSystem->Load("${FASTJET}/lib/libfastjet.so");
  
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libJETAN.so");

  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTJET.so");
  
  TFile* f = new TFile(fileName);

  AliHLTJETAnalysisJets* jets = NULL;
  
  jets = static_cast<AliHLTJETAnalysisJets*>(f->Get("AliHLTJETAnalysisJets"));
  
  AliHLTJETAnalysisMerge* merge = new AliHLTJETAnalysisMerge();
  merge->Initialize();

  merge->AddJets( jets ); 
    
  merge->CreateCanvas();
}
