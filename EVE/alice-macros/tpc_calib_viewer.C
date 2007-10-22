void tpc_calib_viewer(const char* file="CalibTree.root")
{
  Reve::RGBrowser* b = gReve->GetBrowser();
  b->StartEmbedding(1);
  AliTPCCalibViewerGUI::ShowGUI(file);
  b->StopEmbedding();
  b->SetTabTitle("TPC Calib", 1);
}
