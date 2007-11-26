void tpc_calib_viewer(const char* file="CalibTree.root")
{
   Reve::RGBrowser* b = gReve->GetBrowser();
   b->StartEmbedding(1);

   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("AliTPCCalibViewer GUI");
   frmMain->SetCleanup(kDeepCleanup);

   AliTPCCalibViewerGUI* calibViewer1 = new AliTPCCalibViewerGUI(frmMain, 1000, 600, (char*)file);
   frmMain->AddFrame(calibViewer1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();

   b->StopEmbedding();
   b->SetTabTitle("TPC CalibViewer", 1);
}
