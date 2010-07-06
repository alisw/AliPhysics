//------------------------------------------------------------------------------
// logoPrelim.C
//
// insert alice logo and/or Preliminiary tag
//------------------------------------------------------------------------------



logoPrelim(TCanvas* can) {

if (SHOW_LOGO) {
    TPad *pad_logo = new TPad("pad_logo","pad_logo",0.754,0.829,0.930,0.957);
    pad_logo->SetMargin(0,0,0,0);
//     pad_logo->SetBorderMode(1);
//     pad_logo->SetBorderSize(1);
    can->cd();
    pad_logo->Draw();
    pad_logo->cd();
    logo->Draw();
}

if (SHOW_PRELIM) {
//     TPad *pad_prel = new TPad("pad_prel","pad_prel",0.587,0.914,0.947,0.964);
    TPad *pad_prel = new TPad("pad_prel","pad_prel",0.200,0.650,0.560,0.700);
    pad_prel->SetMargin(0,0,0,0);
//     pad_prel->SetBorderMode(1);
    can->cd();
    pad_prel->Draw();
    pad_prel->cd();
    TLatex lprel;
    lprel.SetTextColor(kRed);
    lprel.SetTextSizePixels(24);
    lprel.DrawLatex(0.04,0.3,"ALICE Preliminary"); 
}
}