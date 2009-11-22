void visscan_raw_raw(const TString& path = ".", Bool_t show_extra_geo = kFALSE)
{
  AliEveEventManager::SearchRawForCentralReconstruction();
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"raw://\", \"%s\", %d)",
				     path.Data(), show_extra_geo));
}
