{
{

  gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(10);

  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(-1);
  gStyle->SetFrameLineWidth(1.2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineWidth(1.2);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadBorderMode(-1);

  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetStatColor(10);
  gStyle->SetTitleColor(10);

  gStyle->SetStatX(0.8);
  gStyle->SetStatW(0.1);
  gStyle->SetStatH(0.06);

//  gStyle->SetTitleFont(42);//org
//  gStyle->SetLabelFont(42,"X");//org
  gStyle->SetLabelFont(72,"Y");//org
  gStyle->SetLabelSize(0.04,"xyz");//no org
  gStyle->SetTitleSize(0.04,"xyz");//no org
  gStyle->SetLabelFont(72,"XY");//no org
  gStyle->SetTitleFont(72);//no org
  gStyle->SetStatFont(72);
  gStyle->SetLabelOffset(0.015,"X");//no org
  gStyle->SetLabelOffset(0.015,"Y");//no org

  gStyle->SetTitleOffset(1.4,"X");
  gStyle->SetTitleOffset(1.7,"Y");

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);


}
}
