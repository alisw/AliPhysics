/// \file startGUItime.C
/// \brief Global variables to access the gui components

AliTPCCalibViewerGUItime *guiTime=0x0;
TObjArray *arrGUI=0x0;
AliTPCCalibViewerGUI *gui=0x0;


void startGUItime(const char* file=""){
  TString filename(file);
  if (filename.IsNull()){
    filename=gSystem->ExpandPathName("$GUI_OUTDIR_TIME");
    if (filename=="$GUI_OUTDIR_TIME") return;
    if (!filename.EndsWith("/")) filename+="/";
    filename+="calibTreeTime_*.root";
  }
  arrGUI=AliTPCCalibViewerGUItime::ShowGUI(filename.Data());
  guiTime=(AliTPCCalibViewerGUItime*)arrGUI->At(0);
  gui=(AliTPCCalibViewerGUI*)arrGUI->At(1);
  TString cacheDir=gSystem->ExpandPathName("$GUI_OUTDIR_RUNS");
  if (cacheDir=="$GUI_OUTDIR_RUNS") cacheDir="/tmp";
  guiTime->SetCacheDir(cacheDir.Data());
}
