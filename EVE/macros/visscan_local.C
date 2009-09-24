void visscan_local(const TString& path = ".", Bool_t show_extra_geo = kFALSE)
{
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"local://$ALICE_ROOT/OCDB\", \"%s\", %d)",
				     path.Data(), show_extra_geo));
}
