//
// Script to draw detail of the FMD
//
void DrawFMD2()
{
  // gAlice->Init("FMD/scripts/ConfigInner.C");
  gAlice->Init("$(ALICE)/FMD/Config.C");
  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  gROOT->LoadMacro("$(ALICE)/FMD/ViewFMD.C");
  gInterpreter->ProcessLine("ViewFMD()");
  gMC->Gsatt("FMD3", "seen", -1);
  gMC->Gsatt("FMD1", "seen", -1);
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  // gMC->Gdraw("alic", 90, 0, 0, 28, 10, .25, .25);
  // gMC->Gdraw("alic", 90, 0, 0, 28, 10, .25, .25);
  gMC->Gdraw("alic", 179, 0, 0, 10, 10, .25, .25);

#if 0
  TArrow* a1 = new TArrow(13.5, 16, 15, 18., .03, "<|");
  a1->SetAngle(30);
  a1->SetFillColor(1);
  a1->Draw();
  
  TLatex* l1 = new TLatex(15, 18, "Honeycomb");
  l1->SetTextAlign(12);
  l1->SetTextFont(132);
  l1->SetTextSize(.04);
  l1->Draw();
  
  a1->DrawArrow(13.4, 14., 15, 15, .03, "<|");
  l1->DrawLatex(15, 15, "Support Leg");
  
  a1->DrawArrow(10.7, 14.2, 15, 13, .03, "<|");
  l1->DrawLatex(15, 13, "Print board");

  a1->DrawArrow(9.7, 12.5, 15, 11, .03, "<|");
  l1->DrawLatex(15, 11, "Silicon sensor");

  a1->DrawArrow(8.3, 12.7, 7, 15, .03, "<|");
  TLatex* l2 = new TLatex(7, 15, "Support Cone");
  l2->SetTextSize(.04);
  l2->SetTextFont(132);
  l2->SetTextAlign(32);
  l2->Draw();

  TLatex* l3 = new TLatex(3, 3, "FMD2");
  l3->SetTextSize(.06);
  l3->SetTextFont(132);
  l3->Draw();  
#endif

  gPad->Modified();
  gPad->cd();
  gPad->Print("FMD2.png");
}
