void visscan_mcresidual(const TString& path = ".", Bool_t show_extra_geo = kFALSE)
{
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"mcresidual://\", \"%s\", %d)",
				     path.Data(), show_extra_geo));
}
