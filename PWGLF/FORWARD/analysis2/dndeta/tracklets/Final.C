void
Final(Bool_t reweighed=false)
{
  Int_t       run   = 245064;
  const char* mcPer = "LHC15k1_plus21";
  const char* dtPer = "LHC15o";
  const char* mcPst = (reweighed ? "_reweighted" : "");
  const char* prefi = (reweighed ? "reweighed" : "normal");
  TString     dtIn;
  TString     mcIn;
  TString     dOut;
  dtIn.Form("%s_%d_tracklets/root_archive_%09d/trdt.root",dtPer,run,run);
  mcIn.Form("%s_%d_tracklets%s/root_archive_%06d/trmc.root",
	    mcPer,run,mcPst,run);
  dOut.Form("%sResult/TRACKLETS_05023_PbPb.input", prefi);
  TString  rout;
  rout.Form("PbPb5023midRapidity%s.root", reweighed ? "Reweighed" : "Normal");

  gROOT->LoadMacro("SaveCanvas.C+g");
  gROOT->LoadMacro("CorrectSpectraMultiMCBG.C+g");
  gSystem->mkdir("corrFig");
  gSystem->mkdir("corrRes");
    
  CorrectSpectraMultiMCBG(dtIn, mcIn, "PbPb", 9, false);
  
  gROOT->LoadMacro("Extract.C");
  Extract("corrRes/PbPb_9bins_CutEta-2.0_2.0_Zv-15.0_15.0_bg_Shape_wdst_mcLB2_cutSig1.5_cutBg5.0.root","corrRes/TRACKLETS_05023_PbPb.input",reweighed);

  gSystem->Rename("corrFig", Form("%sFigures", prefi));
  gSystem->Rename("corrRes", Form("%sResult", prefi));

  gROOT->LoadMacro("$ANA_SRC/dndeta/QuickDraw.C");
  const char* what[] = {
    // dOut.Data(),
    rout.Data(),
    "../../data_repo/other/fwd/PbPb_05023_CENT.root",
    0 };
  QuickDraw(what);

  gPad->SaveAs(Form("dNdeta_05023_%s.pdf",
		    reweighed ? "Reweighed" : "Normal"));
}
