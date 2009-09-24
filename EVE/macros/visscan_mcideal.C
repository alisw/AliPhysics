void visscan_mcideal(const TString& path = ".", Bool_t show_extra_geo = kFALSE)
{
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"mcideal://\", \"%s\", %d)",
				     path.Data(), show_extra_geo));
}
