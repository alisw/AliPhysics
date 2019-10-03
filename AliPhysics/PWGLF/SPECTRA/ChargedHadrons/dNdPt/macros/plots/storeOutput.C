//------------------------------------------------------------------------------
// storeOutput.C
//
// store figures as eps+gif
//------------------------------------------------------------------------------


void storeOutput()
{
if (!SAVE_FIGURES) { return; }

TString outname = outputDir + "/ALICE_";
if (SHOW_PRELIM) { outname += "Preliminary_"; }
TString tmp;


tmp = outname+"Yield_INEL_2PiPt_dNdPt_withPowerLaw_noBinShift_sepRatios.eps";
can23->SaveAs(tmp.Data());

tmp = outname+"ATLAS_CMS_Yield_NSD_2PiPt_dNdPt_noBinShift.eps";
can3->SaveAs(tmp.Data());

tmp = outname+"UA1_InvYield_NSD_Ep_2PiPt_dNdPt_noBinShift.eps";
can4->SaveAs(tmp.Data());

tmp = outname+"PHOJET_PYTHIA_Yield_INEL_2PiPt_dNdPt_noBinShift.eps";
can5->SaveAs(tmp.Data());


tmp = outname+"Yield_INEL_2PiPt_dNdPt_withPowerLaw_noBinShift_sepRatios.gif";
can23->SaveAs(tmp.Data());

tmp = outname+"ATLAS_CMS_Yield_NSD_2PiPt_dNdPt_noBinShift.gif";
can3->SaveAs(tmp.Data());

tmp = outname+"UA1_InvYield_NSD_Ep_2PiPt_dNdPt_noBinShift.gif";
can4->SaveAs(tmp.Data());

tmp = outname+"PHOJET_PYTHIA_Yield_INEL_2PiPt_dNdPt_noBinShift.gif";
can5->SaveAs(tmp.Data());

} 