void drawText(Int_t option=1)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  if(option==1)
      latex.DrawLatex(0.5, 0.95, "pp @ 7 TeV (LHC 10 b)");
  else {
    latex.SetTextAngle(90.0);
    latex.DrawLatex(0.95, 0.2, "pp @ 7 TeV (LHC 10 b)");
  }
}
