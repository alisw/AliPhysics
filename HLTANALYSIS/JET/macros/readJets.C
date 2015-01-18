void readJets( Char_t* fileName = "./analysis/EOR_analyze_100000_kPythia6Jets104_125.root" ) {

  if ( getenv("FASTJET") ) {
    gSystem->Load("libCGAL");
    gSystem->Load("${FASTJET}/lib/libfastjet");
    gSystem->Load("${FASTJET}/lib/libSISConePlugin");
  }

  gSystem->Load("libTree");
  gSystem->Load("libPhysics");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN");

  gSystem->Load("libAliHLTUtil");
  gSystem->Load("libAliHLTJET");
  
  TFile* f = new TFile(fileName);

  AliHLTJETAnalysisJets* jets = NULL;
  
  jets = static_cast<AliHLTJETAnalysisJets*>(f->Get("AliHLTJETAnalysisJets"));
  
  AliHLTJETAnalysisMerge* merge = new AliHLTJETAnalysisMerge();
  merge->Initialize();
    
  merge->AddJets( jets ); 
    
  merge->CreateCanvas();
}
