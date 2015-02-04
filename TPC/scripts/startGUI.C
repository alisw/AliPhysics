/// \file startGUI.C

TObjArray *arrGUI=0x0;
TTree *calPads=0x0;
TTree *calPads2=0x0;
AliTPCCalibViewerGUI *gui=0x0;
AliTPCCalibViewerGUI *gui2=0x0;
AliTPCCalibViewer *viewer=0x0;
AliTPCCalibViewer *viewer2=0x0;

startGUI(char *file=0x0){
  if ( file ) {
    arrGUI=AliTPCCalibViewerGUI::ShowGUI(file);
    gui=((AliTPCCalibViewerGUI*)arrGUI->At(0));
    gui2=((AliTPCCalibViewerGUI*)arrGUI->At(1));
    viewer=gui->GetViewer();
    viewer2=gui2->GetViewer();
    calPads=viewer->GetTree();
    calPads2=viewer2->GetTree();
  }
}
