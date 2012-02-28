AliFMDAnaCalibBackgroundCorrection*
ViewBGCor(const char* filename)
{
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libPWG0base");
  // gSystem->Load("libPWGLFforward");


  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("ViewBgCor", "Cannot open file %s", filename);
    return 0;
  }

  AliFMDAnaCalibBackgroundCorrection* bg = 
    static_cast<AliFMDAnaCalibBackgroundCorrection*>(file->Get("background"));
  
  return bg;
}

