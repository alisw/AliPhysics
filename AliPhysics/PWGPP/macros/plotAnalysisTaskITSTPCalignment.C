void plotAnalysisTaskITSTPCalignment(const char* option = "b")
{
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetHistFillColor(17);

  TString optionstr(option);
  gROOT->LoadMacro("AliRelAlignerKalmanArray.cxx++");
  TFile f("outputITSTPCalignment.root");

  //////////////////////////////////////////////////////////////////////////////////////////////////
  AliRelAlignerKalmanArray* fArray = dynamic_cast<AliRelAlignerKalmanArray*>(f.Get("outputArrayITSsa"));
  if (!fArray)
  {
    printf("fArray cannot be read!\n");
    return;
  }
  else
  {
    TCanvas* c1 = new TCanvas("c1","psi in time");
    fArray->MakeGraph(0)->Draw("A*");
    c1->SaveAs("graphpsi.eps");

    TCanvas* c2 = new TCanvas("c2","theta in time");
    fArray->MakeGraph(1)->Draw("A*");
    c2->SaveAs("graphtheta.eps");

    TCanvas* c3 = new TCanvas("c3","phi in time");
    fArray->MakeGraph(2)->Draw("A*");
    c3->SaveAs("graphphi.eps");

    TCanvas* c4 = new TCanvas("c4","x in time");
    fArray->MakeGraph(3)->Draw("A*");
    c4->SaveAs("graphx.eps");

    TCanvas* c5 = new TCanvas("c5","y in time");
    fArray->MakeGraph(4)->Draw("A*");
    c5->SaveAs("graphy.eps");

    TCanvas* c6 = new TCanvas("c6","z in time");
    fArray->MakeGraph(5)->Draw("A*");
    c6->SaveAs("graphz.eps");

    TCanvas* c7 = new TCanvas("c7","TPC vd correction in time");
    fArray->MakeGraph(6)->Draw("A*");
    c7->SaveAs("graphvd.eps");

    TCanvas* c8 = new TCanvas("c8","TPC t0 correction in time");
    fArray->MakeGraph(7)->Draw("A*");
    c8->SaveAs("grapht0.eps");

    TCanvas* c9 = new TCanvas("c9","TPC dv/dy in time");
    fArray->MakeGraph(8)->Draw("A*");
    c9->SaveAs("graphdvdy.eps");
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  TList* fList = dynamic_cast<TList*>(f.Get("outputList"));
  TList* pList = NULL;

  TH1F* pMatchingEfficiency = dynamic_cast<TH1F*>(fList->At(3));
  TCanvas* canvasMEff = new TCanvas();
  pMatchingEfficiency->DrawCopy();
  
  
  //------------------------------------------------------------------------------------
  pList = dynamic_cast<TList*>(fList->At(0));
  TH2F* pZYAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(0));
  TH2F* pZYCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(1));
  TH2F* pLPAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(2));
  TH2F* pLPCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(3));
  TH2F* pPhiYAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(4));
  TH2F* pPhiYCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(5));
  TH2F* pPhiZAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(6));
  TH2F* pPhiZCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(7));
  TH2F* pPtYAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(8));
  TH2F* pPtYCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(9));
  TH2F* pPtZAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(10));
  TH2F* pPtZCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(11));
  TH2F* pLowPtYAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(12));
  TH2F* pLowPtYCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(13));
  TH2F* pLowPtZAResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(14));
  TH2F* pLowPtZCResidualsHistBpos = dynamic_cast<TH2F*>(pList->At(15));
  
  TH1D* resycBpos = pZYCResidualsHistBpos->ProjectionY();
  resycBpos->SetTitle("r\\phi residual distribution side C (z<0), Bpos");
  TH1D* reszcBpos = pZYCResidualsHistBpos->ProjectionX();
  reszcBpos->SetTitle("z residual distribution side C (z<0), Bpos");
  TH1D* respcBpos = pLPCResidualsHistBpos->ProjectionY();
  respcBpos->SetTitle("sin(\\phi) residual distribution side C (z<0), Bpos");
  TH1D* reslcBpos = pLPCResidualsHistBpos->ProjectionX();
  reslcBpos->SetTitle("tan(\\lambda) residual distribution side C (z<0), Bpos");
  TProfile* profphiycBpos = pPhiYCResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profphiycBpos->SetTitle("\\phi profile of r\\phi residuals side C (z<0), Bpos");
  profphiycBpos->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profphizcBpos = pPhiZCResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profphizcBpos->SetTitle("\\phi profile of z residuals side C (z<0), Bpos");
  profphizcBpos->SetYTitle("\\deltaz [cm]");
  TProfile* profptycBpos = pPtYCResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profptycBpos->SetTitle("pt profile of r\\phi residuals side C (z<0), Bpos");
  profptycBpos->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profptzcBpos = pPtZCResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profptzcBpos->SetTitle("pt profile of z residuals side C (z<0), Bpos");
  profptzcBpos->SetYTitle("\\deltaz [cm]");
  TProfile* proflowptycBpos = pLowPtYCResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  proflowptycBpos->SetTitle("pt profile of r\\phi residuals side C (z<0), Bpos");
  proflowptycBpos->SetYTitle("\\deltar\\phi [cm]");
  TProfile* proflowptzcBpos = pLowPtZCResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  proflowptzcBpos->SetTitle("pt profile of z residuals side C (z<0), Bpos");
  proflowptzcBpos->SetYTitle("\\deltaz [cm]");

  TH1D* resyaBpos = pZYAResidualsHistBpos->ProjectionY();
  resyaBpos->SetTitle("r\\phi residual distribution side A (z>0), Bpos");
  TH1D* reszaBpos = pZYAResidualsHistBpos->ProjectionX();
  reszaBpos->SetTitle("z residual distribution side A (z>0), Bpos");
  TH1D* respaBpos = pLPAResidualsHistBpos->ProjectionY();
  respaBpos->SetTitle("sin(\\phi) residual distribution side A (z>0), Bpos");
  TH1D* reslaBpos = pLPAResidualsHistBpos->ProjectionX();
  reslaBpos->SetTitle("tan(\\lambda) residual distribution side A (z>0), Bpos");
  TProfile* profphiyaBpos = pPhiYAResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profphiyaBpos->SetTitle("\\phi profile of r\\phi residuals side A (z>0), Bpos");
  profphiyaBpos->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profphizaBpos = pPhiZAResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profphizaBpos->SetTitle("\\phi profile of z residuals side A (z>0), Bpos");
  profphizaBpos->SetYTitle("\\deltaz [cm]");
  TProfile* profptyaBpos = pPtYAResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profptyaBpos->SetTitle("pt profile of r\\phi residuals side A (z>0), Bpos");
  profptyaBpos->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profptzaBpos = pPtZAResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  profptzaBpos->SetTitle("pt profile of z residuals side A (z>0), Bpos");
  profptzaBpos->SetYTitle("\\deltaz [cm]");
  TProfile* proflowptyaBpos = pLowPtYAResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  proflowptyaBpos->SetTitle("pt profile of r\\phi residuals side A (z>0), Bpos");
  proflowptyaBpos->SetYTitle("\\deltar\\phi [cm]");
  TProfile* proflowptzaBpos = pLowPtZAResidualsHistBpos->ProfileX("_pfx",1,-1,"s");
  proflowptzaBpos->SetTitle("pt profile of z residuals side A (z>0), Bpos");
  proflowptzaBpos->SetYTitle("\\deltaz [cm]");

  if (pZYAResidualsHistBpos->GetEntries()>0)
  {
  TCanvas* c1Bpos = new TCanvas("c1Bpos","Residuals ZY 2D, Bpos",800,300);
  c1Bpos->Divide(2,1);
  c1Bpos->cd(1); pZYAResidualsHistBpos->DrawCopy("col");
  c1Bpos->cd(2); pZYCResidualsHistBpos->DrawCopy("col");
  c1Bpos->cd(0);
  c1Bpos->SaveAs("residualsZY2D-Bpos.eps");
  
  TCanvas* c2Bpos = new TCanvas("c2Bpos","Residuals LP 2D, Bpos",800,300);
  c2Bpos->Divide(2,1);
  c2Bpos->cd(1); pLPAResidualsHistBpos->DrawCopy("col");
  c2Bpos->cd(2); pLPCResidualsHistBpos->DrawCopy("col");
  c2Bpos->cd(0);
  c2Bpos->SaveAs("residualsLP2D-Bpos.eps");

  TCanvas* c3Bpos = new TCanvas("c3Bpos","Residuals z 1D, Bpos",800,300);
  c3Bpos->Divide(2,1);
  c3Bpos->cd(1); reszaBpos->DrawCopy(); 
  c3Bpos->cd(2); reszcBpos->DrawCopy(); 
  c3Bpos->cd(0);
  c3Bpos->SaveAs("residualsZ-Bpos.eps");

  TCanvas* c4Bpos = new TCanvas("c4Bpos","Residuals y 1D, Bpos",800,300);
  c4Bpos->Divide(2,1);
  c4Bpos->cd(1); resyaBpos->DrawCopy();
  c4Bpos->cd(2); resycBpos->DrawCopy();
  c4Bpos->cd(0);
  c4Bpos->SaveAs("residualsY-Bpos.eps");

  TCanvas* c5Bpos = new TCanvas("c5Bpos","Residuals phi 1D, Bpos",800,300);
  c5Bpos->Divide(2,1);
  c5Bpos->cd(1); respaBpos->DrawCopy(); 
  c5Bpos->cd(2); respcBpos->DrawCopy(); 
  c5Bpos->cd(0);
  c5Bpos->SaveAs("residualsPhi-Bpos.eps");

  TCanvas* c6Bpos = new TCanvas("c6Bpos","Residuals lambda 1D, Bpos",800,300);
  c6Bpos->Divide(2,1);
  c6Bpos->cd(1); reslaBpos->DrawCopy();
  c6Bpos->cd(2); reslcBpos->DrawCopy();
  c6Bpos->cd(0);
  c6Bpos->SaveAs("residualsLambda-Bpos.eps");

  TCanvas* c7Bpos = new TCanvas("c7Bpos","Profiles z in phi, Bpos",800,300);
  c7Bpos->Divide(2,1);
  c7Bpos->cd(1); profphizaBpos->DrawCopy();
  c7Bpos->cd(2); profphizcBpos->DrawCopy();
  c7Bpos->cd(0);
  c7Bpos->SaveAs("residualsProfilePhiZ-Bpos.eps");

  TCanvas* c8Bpos = new TCanvas("c8Bpos","Profiles y in phi, Bpos",800,300);
  c8Bpos->Divide(2,1);
  c8Bpos->cd(1); profphiyaBpos->DrawCopy();
  c8Bpos->cd(2); profphiycBpos->DrawCopy();
  c8Bpos->cd(0);
  c8Bpos->SaveAs("residualsProfilePhiY-Bpos.eps");

  TCanvas* c9Bpos = new TCanvas("c9Bpos", "Residuals Phi-Z 2D, Bpos",800,300);
  c9Bpos->Divide(2,1);
  c9Bpos->cd(1); pPhiZAResidualsHistBpos->DrawCopy("col");
  c9Bpos->cd(2); pPhiZCResidualsHistBpos->DrawCopy("col");
  c9Bpos->cd(0); 
  c9Bpos->SaveAs("residualsPhiZ2D-Bpos.eps");

  TCanvas* c1Bpos0 = new TCanvas("c1Bpos0", "Residuals Phi-Y 2D, Bpos",800,300);
  c1Bpos0->Divide(2,1);
  c1Bpos0->cd(1); pPhiYAResidualsHistBpos->DrawCopy("col");
  c1Bpos0->cd(2); pPhiYCResidualsHistBpos->DrawCopy("col");
  c1Bpos0->cd(0); 
  c1Bpos0->SaveAs("residualsPhiY2D-Bpos.eps");

  TCanvas* c11Bpos = new TCanvas("c11Bpos", "Residuals Pt-Z 2D, Bpos",800,300);
  c11Bpos->Divide(2,1);
  c11Bpos->cd(1); pPtZAResidualsHistBpos->DrawCopy("col");
  c11Bpos->cd(2); pPtZCResidualsHistBpos->DrawCopy("col");
  c11Bpos->cd(0);
  c11Bpos->SaveAs("residualsPtZ2D-Bpos.eps");

  TCanvas* c12Bpos = new TCanvas("c12Bpos", "Residuals Pt-Y 2D, Bpos",800,300);
  c12Bpos->Divide(2,1);
  c12Bpos->cd(1); pPtYAResidualsHistBpos->DrawCopy("col");
  c12Bpos->cd(2); pPtYCResidualsHistBpos->DrawCopy("col");
  c12Bpos->cd(0);
  c12Bpos->SaveAs("residualsPtY2D-Bpos.eps");

  TCanvas* c13Bpos = new TCanvas("c13Bpos","Profiles Pt-Z, Bpos",800,300);
  c13Bpos->Divide(2,1);
  c13Bpos->cd(1); profptzaBpos->DrawCopy();
  c13Bpos->cd(2); profptzcBpos->DrawCopy();
  c13Bpos->cd(0);
  c13Bpos->SaveAs("residualsProfilePtZ-Bpos.eps");

  TCanvas* c14Bpos = new TCanvas("c14Bpos","Profiles Pt-Y, Bpos",800,300);
  c14Bpos->Divide(2,1);
  c14Bpos->cd(1); profptyaBpos->DrawCopy();
  c14Bpos->cd(2); profptycBpos->DrawCopy();
  c14Bpos->cd(0);
  c14Bpos->SaveAs("residualsProfilePtY-Bpos.eps");

  TCanvas* c15Bpos = new TCanvas("c15Bpos", "Residuals low Pt-Z 2D, Bpos",800,300);
  c15Bpos->Divide(2,1);
  c15Bpos->cd(1); pLowPtZAResidualsHistBpos->DrawCopy("col");
  c15Bpos->cd(2); pLowPtZCResidualsHistBpos->DrawCopy("col");
  c15Bpos->cd(0);
  c15Bpos->SaveAs("residualsLowPtZ2D-Bpos.eps");

  TCanvas* c16Bpos = new TCanvas("c16Bpos", "Residuals in low PtY 2D, Bpos",800,300);
  c16Bpos->Divide(2,1);
  c16Bpos->cd(1); pLowPtYAResidualsHistBpos->DrawCopy("col");
  c16Bpos->cd(2); pLowPtYCResidualsHistBpos->DrawCopy("col");
  c16Bpos->cd(0);
  c16Bpos->SaveAs("residualsLowPtY2D-Bpos.eps");

  TCanvas* c17Bpos = new TCanvas("c17Bpos","Profiles low Pt-Z, Bpos",800,300);
  c17Bpos->Divide(2,1);
  c17Bpos->cd(1); proflowptzaBpos->DrawCopy();
  c17Bpos->cd(2); proflowptzcBpos->DrawCopy();
  c17Bpos->cd(0);
  c17Bpos->SaveAs("residualsProfileLowPtZ-Bpos.eps");

  TCanvas* c18Bpos = new TCanvas("c18Bpos","Profiles low Pt-Y, Bpos",800,300);
  c18Bpos->Divide(2,1);
  c18Bpos->cd(1); proflowptyaBpos->DrawCopy();
  c18Bpos->cd(2); proflowptycBpos->DrawCopy();
  c18Bpos->cd(0);
  c18Bpos->SaveAs("residualsProfileLowPtY-Bpos.eps");
  }

  //------------------------------------------------------------------------------------
  pList = dynamic_cast<TList*>(fList->At(1));
  TH2F* pZYAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(0));
  TH2F* pZYCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(1));
  TH2F* pLPAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(2));
  TH2F* pLPCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(3));
  TH2F* pPhiYAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(4));
  TH2F* pPhiYCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(5));
  TH2F* pPhiZAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(6));
  TH2F* pPhiZCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(7));
  TH2F* pPtYAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(8));
  TH2F* pPtYCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(9));
  TH2F* pPtZAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(10));
  TH2F* pPtZCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(11));
  TH2F* pLowPtYAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(12));
  TH2F* pLowPtYCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(13));
  TH2F* pLowPtZAResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(14));
  TH2F* pLowPtZCResidualsHistBneg = dynamic_cast<TH2F*>(pList->At(15));
  
  TH1D* resycBneg = pZYCResidualsHistBneg->ProjectionY();
  resycBneg->SetTitle("r\\phi residual distribution side C (z<0), Bneg");
  TH1D* reszcBneg = pZYCResidualsHistBneg->ProjectionX();
  reszcBneg->SetTitle("z residual distribution side C (z<0), Bneg");
  TH1D* respcBneg = pLPCResidualsHistBneg->ProjectionY();
  respcBneg->SetTitle("sin(\\phi) residual distribution side C (z<0), Bneg");
  TH1D* reslcBneg = pLPCResidualsHistBneg->ProjectionX();
  reslcBneg->SetTitle("tan(\\lambda) residual distribution side C (z<0), Bneg");
  TProfile* profphiycBneg = pPhiYCResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profphiycBneg->SetTitle("\\phi profile of r\\phi residuals side C (z<0), Bneg");
  profphiycBneg->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profphizcBneg = pPhiZCResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profphizcBneg->SetTitle("\\phi profile of z residuals side C (z<0), Bneg");
  profphizcBneg->SetYTitle("\\deltaz [cm]");
  TProfile* profptycBneg = pPtYCResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profptycBneg->SetTitle("pt profile of r\\phi residuals side C (z<0), Bneg");
  profptycBneg->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profptzcBneg = pPtZCResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profptzcBneg->SetTitle("pt profile of z residuals side C (z<0), Bneg");
  profptzcBneg->SetYTitle("\\deltaz [cm]");
  TProfile* proflowptycBneg = pLowPtYCResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  proflowptycBneg->SetTitle("pt profile of r\\phi residuals side C (z<0), Bneg");
  proflowptycBneg->SetYTitle("\\deltar\\phi [cm]");
  TProfile* proflowptzcBneg = pLowPtZCResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  proflowptzcBneg->SetTitle("pt profile of z residuals side C (z<0), Bneg");
  proflowptzcBneg->SetYTitle("\\deltaz [cm]");

  TH1D* resyaBneg = pZYAResidualsHistBneg->ProjectionY();
  resyaBneg->SetTitle("r\\phi residual distribution side A (z>0), Bneg");
  TH1D* reszaBneg = pZYAResidualsHistBneg->ProjectionX();
  reszaBneg->SetTitle("z residual distribution side A (z>0), Bneg");
  TH1D* respaBneg = pLPAResidualsHistBneg->ProjectionY();
  respaBneg->SetTitle("sin(\\phi) residual distribution side A (z>0), Bneg");
  TH1D* reslaBneg = pLPAResidualsHistBneg->ProjectionX();
  reslaBneg->SetTitle("tan(\\lambda) residual distribution side A (z>0), Bneg");
  TProfile* profphiyaBneg = pPhiYAResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profphiyaBneg->SetTitle("\\phi profile of r\\phi residuals side A (z>0), Bneg");
  profphiyaBneg->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profphizaBneg = pPhiZAResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profphizaBneg->SetTitle("\\phi profile of z residuals side A (z>0), Bneg");
  profphizaBneg->SetYTitle("\\deltaz [cm]");
  TProfile* profptyaBneg = pPtYAResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profptyaBneg->SetTitle("pt profile of r\\phi residuals side A (z>0), Bneg");
  profptyaBneg->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profptzaBneg = pPtZAResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  profptzaBneg->SetTitle("pt profile of z residuals side A (z>0), Bneg");
  profptzaBneg->SetYTitle("\\deltaz [cm]");
  TProfile* proflowptyaBneg = pLowPtYAResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  proflowptyaBneg->SetTitle("pt profile of r\\phi residuals side A (z>0), Bneg");
  proflowptyaBneg->SetYTitle("\\deltar\\phi [cm]");
  TProfile* proflowptzaBneg = pLowPtZAResidualsHistBneg->ProfileX("_pfx",1,-1,"s");
  proflowptzaBneg->SetTitle("pt profile of z residuals side A (z>0), Bneg");
  proflowptzaBneg->SetYTitle("\\deltaz [cm]");

  if (pZYAResidualsHistBneg->GetEntries()>0)
  {
  TCanvas* c1Bneg = new TCanvas("c1Bneg","Residuals ZY 2D, Bneg",800,300);
  c1Bneg->Divide(2,1);
  c1Bneg->cd(1); pZYAResidualsHistBneg->DrawCopy("col");
  c1Bneg->cd(2); pZYCResidualsHistBneg->DrawCopy("col");
  c1Bneg->cd(0);
  c1Bneg->SaveAs("residualsZY2D-Bneg.eps");
  
  TCanvas* c2Bneg = new TCanvas("c2Bneg","Residuals LP 2D, Bneg",800,300);
  c2Bneg->Divide(2,1);
  c2Bneg->cd(1); pLPAResidualsHistBneg->DrawCopy("col");
  c2Bneg->cd(2); pLPCResidualsHistBneg->DrawCopy("col");
  c2Bneg->cd(0);
  c2Bneg->SaveAs("residualsLP2D-Bneg.eps");

  TCanvas* c3Bneg = new TCanvas("c3Bneg","Residuals z 1D, Bneg",800,300);
  c3Bneg->Divide(2,1);
  c3Bneg->cd(1); reszaBneg->DrawCopy(); 
  c3Bneg->cd(2); reszcBneg->DrawCopy(); 
  c3Bneg->cd(0);
  c3Bneg->SaveAs("residualsZ-Bneg.eps");

  TCanvas* c4Bneg = new TCanvas("c4Bneg","Residuals y 1D, Bneg",800,300);
  c4Bneg->Divide(2,1);
  c4Bneg->cd(1); resyaBneg->DrawCopy();
  c4Bneg->cd(2); resycBneg->DrawCopy();
  c4Bneg->cd(0);
  c4Bneg->SaveAs("residualsY-Bneg.eps");

  TCanvas* c5Bneg = new TCanvas("c5Bneg","Residuals phi 1D, Bneg",800,300);
  c5Bneg->Divide(2,1);
  c5Bneg->cd(1); respaBneg->DrawCopy(); 
  c5Bneg->cd(2); respcBneg->DrawCopy(); 
  c5Bneg->cd(0);
  c5Bneg->SaveAs("residualsPhi-Bneg.eps");

  TCanvas* c6Bneg = new TCanvas("c6Bneg","Residuals lambda 1D, Bneg",800,300);
  c6Bneg->Divide(2,1);
  c6Bneg->cd(1); reslaBneg->DrawCopy();
  c6Bneg->cd(2); reslcBneg->DrawCopy();
  c6Bneg->cd(0);
  c6Bneg->SaveAs("residualsLambda-Bneg.eps");

  TCanvas* c7Bneg = new TCanvas("c7Bneg","Profiles z in phi, Bneg",800,300);
  c7Bneg->Divide(2,1);
  c7Bneg->cd(1); profphizaBneg->DrawCopy();
  c7Bneg->cd(2); profphizcBneg->DrawCopy();
  c7Bneg->cd(0);
  c7Bneg->SaveAs("residualsProfilePhiZ-Bneg.eps");

  TCanvas* c8Bneg = new TCanvas("c8Bneg","Profiles y in phi, Bneg",800,300);
  c8Bneg->Divide(2,1);
  c8Bneg->cd(1); profphiyaBneg->DrawCopy();
  c8Bneg->cd(2); profphiycBneg->DrawCopy();
  c8Bneg->cd(0);
  c8Bneg->SaveAs("residualsProfilePhiY-Bneg.eps");

  TCanvas* c9Bneg = new TCanvas("c9Bneg", "Residuals Phi-Z 2D, Bneg",800,300);
  c9Bneg->Divide(2,1);
  c9Bneg->cd(1); pPhiZAResidualsHistBneg->DrawCopy("col");
  c9Bneg->cd(2); pPhiZCResidualsHistBneg->DrawCopy("col");
  c9Bneg->cd(0); 
  c9Bneg->SaveAs("residualsPhiZ2D-Bneg.eps");

  TCanvas* c1Bneg0 = new TCanvas("c1Bneg0", "Residuals Phi-Y 2D, Bneg",800,300);
  c1Bneg0->Divide(2,1);
  c1Bneg0->cd(1); pPhiYAResidualsHistBneg->DrawCopy("col");
  c1Bneg0->cd(2); pPhiYCResidualsHistBneg->DrawCopy("col");
  c1Bneg0->cd(0); 
  c1Bneg0->SaveAs("residualsPhiY2D-Bneg.eps");

  TCanvas* c11Bneg = new TCanvas("c11Bneg", "Residuals Pt-Z 2D, Bneg",800,300);
  c11Bneg->Divide(2,1);
  c11Bneg->cd(1); pPtZAResidualsHistBneg->DrawCopy("col");
  c11Bneg->cd(2); pPtZCResidualsHistBneg->DrawCopy("col");
  c11Bneg->cd(0);
  c11Bneg->SaveAs("residualsPtZ2D-Bneg.eps");

  TCanvas* c12Bneg = new TCanvas("c12Bneg", "Residuals Pt-Y 2D, Bneg",800,300);
  c12Bneg->Divide(2,1);
  c12Bneg->cd(1); pPtYAResidualsHistBneg->DrawCopy("col");
  c12Bneg->cd(2); pPtYCResidualsHistBneg->DrawCopy("col");
  c12Bneg->cd(0);
  c12Bneg->SaveAs("residualsPtY2D-Bneg.eps");

  TCanvas* c13Bneg = new TCanvas("c13Bneg","Profiles Pt-Z, Bneg",800,300);
  c13Bneg->Divide(2,1);
  c13Bneg->cd(1); profptzaBneg->DrawCopy();
  c13Bneg->cd(2); profptzcBneg->DrawCopy();
  c13Bneg->cd(0);
  c13Bneg->SaveAs("residualsProfilePtZ-Bneg.eps");

  TCanvas* c14Bneg = new TCanvas("c14Bneg","Profiles Pt-Y, Bneg",800,300);
  c14Bneg->Divide(2,1);
  c14Bneg->cd(1); profptyaBneg->DrawCopy();
  c14Bneg->cd(2); profptycBneg->DrawCopy();
  c14Bneg->cd(0);
  c14Bneg->SaveAs("residualsProfilePtY-Bneg.eps");

  TCanvas* c15Bneg = new TCanvas("c15Bneg", "Residuals low Pt-Z 2D, Bneg",800,300);
  c15Bneg->Divide(2,1);
  c15Bneg->cd(1); pLowPtZAResidualsHistBneg->DrawCopy("col");
  c15Bneg->cd(2); pLowPtZCResidualsHistBneg->DrawCopy("col");
  c15Bneg->cd(0);
  c15Bneg->SaveAs("residualsLowPtZ2D-Bneg.eps");

  TCanvas* c16Bneg = new TCanvas("c16Bneg", "Residuals in low PtY 2D, Bneg",800,300);
  c16Bneg->Divide(2,1);
  c16Bneg->cd(1); pLowPtYAResidualsHistBneg->DrawCopy("col");
  c16Bneg->cd(2); pLowPtYCResidualsHistBneg->DrawCopy("col");
  c16Bneg->cd(0);
  c16Bneg->SaveAs("residualsLowPtY2D-Bneg.eps");

  TCanvas* c17Bneg = new TCanvas("c17Bneg","Profiles low Pt-Z, Bneg",800,300);
  c17Bneg->Divide(2,1);
  c17Bneg->cd(1); proflowptzaBneg->DrawCopy();
  c17Bneg->cd(2); proflowptzcBneg->DrawCopy();
  c17Bneg->cd(0);
  c17Bneg->SaveAs("residualsProfileLowPtZ-Bneg.eps");

  TCanvas* c18Bneg = new TCanvas("c18Bneg","Profiles low Pt-Y, Bneg",800,300);
  c18Bneg->Divide(2,1);
  c18Bneg->cd(1); proflowptyaBneg->DrawCopy();
  c18Bneg->cd(2); proflowptycBneg->DrawCopy();
  c18Bneg->cd(0);
  c18Bneg->SaveAs("residualsProfileLowPtY-Bneg.eps");
  }

  //------------------------------------------------------------------------------------
  pList = dynamic_cast<TList*>(fList->At(2));
  TH2F* pZYAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(0));
  TH2F* pZYCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(1));
  TH2F* pLPAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(2));
  TH2F* pLPCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(3));
  TH2F* pPhiYAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(4));
  TH2F* pPhiYCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(5));
  TH2F* pPhiZAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(6));
  TH2F* pPhiZCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(7));
  TH2F* pPtYAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(8));
  TH2F* pPtYCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(9));
  TH2F* pPtZAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(10));
  TH2F* pPtZCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(11));
  TH2F* pLowPtYAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(12));
  TH2F* pLowPtYCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(13));
  TH2F* pLowPtZAResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(14));
  TH2F* pLowPtZCResidualsHistBnil = dynamic_cast<TH2F*>(pList->At(15));
  
  TH1D* resycBnil = pZYCResidualsHistBnil->ProjectionY();
  resycBnil->SetTitle("r\\phi residual distribution side C (z<0), Bnil");
  TH1D* reszcBnil = pZYCResidualsHistBnil->ProjectionX();
  reszcBnil->SetTitle("z residual distribution side C (z<0), Bnil");
  TH1D* respcBnil = pLPCResidualsHistBnil->ProjectionY();
  respcBnil->SetTitle("sin(\\phi) residual distribution side C (z<0), Bnil");
  TH1D* reslcBnil = pLPCResidualsHistBnil->ProjectionX();
  reslcBnil->SetTitle("tan(\\lambda) residual distribution side C (z<0), Bnil");
  TProfile* profphiycBnil = pPhiYCResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profphiycBnil->SetTitle("\\phi profile of r\\phi residuals side C (z<0), Bnil");
  profphiycBnil->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profphizcBnil = pPhiZCResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profphizcBnil->SetTitle("\\phi profile of z residuals side C (z<0), Bnil");
  profphizcBnil->SetYTitle("\\deltaz [cm]");
  TProfile* profptycBnil = pPtYCResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profptycBnil->SetTitle("pt profile of r\\phi residuals side C (z<0), Bnil");
  profptycBnil->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profptzcBnil = pPtZCResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profptzcBnil->SetTitle("pt profile of z residuals side C (z<0), Bnil");
  profptzcBnil->SetYTitle("\\deltaz [cm]");
  TProfile* proflowptycBnil = pLowPtYCResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  proflowptycBnil->SetTitle("pt profile of r\\phi residuals side C (z<0), Bnil");
  proflowptycBnil->SetYTitle("\\deltar\\phi [cm]");
  TProfile* proflowptzcBnil = pLowPtZCResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  proflowptzcBnil->SetTitle("pt profile of z residuals side C (z<0), Bnil");
  proflowptzcBnil->SetYTitle("\\deltaz [cm]");

  TH1D* resyaBnil = pZYAResidualsHistBnil->ProjectionY();
  resyaBnil->SetTitle("r\\phi residual distribution side A (z>0), Bnil");
  TH1D* reszaBnil = pZYAResidualsHistBnil->ProjectionX();
  reszaBnil->SetTitle("z residual distribution side A (z>0), Bnil");
  TH1D* respaBnil = pLPAResidualsHistBnil->ProjectionY();
  respaBnil->SetTitle("sin(\\phi) residual distribution side A (z>0), Bnil");
  TH1D* reslaBnil = pLPAResidualsHistBnil->ProjectionX();
  reslaBnil->SetTitle("tan(\\lambda) residual distribution side A (z>0), Bnil");
  TProfile* profphiyaBnil = pPhiYAResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profphiyaBnil->SetTitle("\\phi profile of r\\phi residuals side A (z>0), Bnil");
  profphiyaBnil->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profphizaBnil = pPhiZAResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profphizaBnil->SetTitle("\\phi profile of z residuals side A (z>0), Bnil");
  profphizaBnil->SetYTitle("\\deltaz [cm]");
  TProfile* profptyaBnil = pPtYAResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profptyaBnil->SetTitle("pt profile of r\\phi residuals side A (z>0), Bnil");
  profptyaBnil->SetYTitle("\\deltar\\phi [cm]");
  TProfile* profptzaBnil = pPtZAResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  profptzaBnil->SetTitle("pt profile of z residuals side A (z>0), Bnil");
  profptzaBnil->SetYTitle("\\deltaz [cm]");
  TProfile* proflowptyaBnil = pLowPtYAResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  proflowptyaBnil->SetTitle("pt profile of r\\phi residuals side A (z>0), Bnil");
  proflowptyaBnil->SetYTitle("\\deltar\\phi [cm]");
  TProfile* proflowptzaBnil = pLowPtZAResidualsHistBnil->ProfileX("_pfx",1,-1,"s");
  proflowptzaBnil->SetTitle("pt profile of z residuals side A (z>0), Bnil");
  proflowptzaBnil->SetYTitle("\\deltaz [cm]");

  if (pZYAResidualsHistBnil->GetEntries()>0)
  {
  TCanvas* c1Bnil = new TCanvas("c1Bnil","Residuals ZY 2D, Bnil",800,300);
  c1Bnil->Divide(2,1);
  c1Bnil->cd(1); pZYAResidualsHistBnil->DrawCopy("col");
  c1Bnil->cd(2); pZYCResidualsHistBnil->DrawCopy("col");
  c1Bnil->cd(0);
  c1Bnil->SaveAs("residualsZY2D-Bnil.eps");
  
  TCanvas* c2Bnil = new TCanvas("c2Bnil","Residuals LP 2D, Bnil",800,300);
  c2Bnil->Divide(2,1);
  c2Bnil->cd(1); pLPAResidualsHistBnil->DrawCopy("col");
  c2Bnil->cd(2); pLPCResidualsHistBnil->DrawCopy("col");
  c2Bnil->cd(0);
  c2Bnil->SaveAs("residualsLP2D-Bnil.eps");

  TCanvas* c3Bnil = new TCanvas("c3Bnil","Residuals z 1D, Bnil",800,300);
  c3Bnil->Divide(2,1);
  c3Bnil->cd(1); reszaBnil->DrawCopy(); 
  c3Bnil->cd(2); reszcBnil->DrawCopy(); 
  c3Bnil->cd(0);
  c3Bnil->SaveAs("residualsZ-Bnil.eps");

  TCanvas* c4Bnil = new TCanvas("c4Bnil","Residuals y 1D, Bnil",800,300);
  c4Bnil->Divide(2,1);
  c4Bnil->cd(1); resyaBnil->DrawCopy();
  c4Bnil->cd(2); resycBnil->DrawCopy();
  c4Bnil->cd(0);
  c4Bnil->SaveAs("residualsY-Bnil.eps");

  TCanvas* c5Bnil = new TCanvas("c5Bnil","Residuals phi 1D, Bnil",800,300);
  c5Bnil->Divide(2,1);
  c5Bnil->cd(1); respaBnil->DrawCopy(); 
  c5Bnil->cd(2); respcBnil->DrawCopy(); 
  c5Bnil->cd(0);
  c5Bnil->SaveAs("residualsPhi-Bnil.eps");

  TCanvas* c6Bnil = new TCanvas("c6Bnil","Residuals lambda 1D, Bnil",800,300);
  c6Bnil->Divide(2,1);
  c6Bnil->cd(1); reslaBnil->DrawCopy();
  c6Bnil->cd(2); reslcBnil->DrawCopy();
  c6Bnil->cd(0);
  c6Bnil->SaveAs("residualsLambda-Bnil.eps");

  TCanvas* c7Bnil = new TCanvas("c7Bnil","Profiles z in phi, Bnil",800,300);
  c7Bnil->Divide(2,1);
  c7Bnil->cd(1); profphizaBnil->DrawCopy();
  c7Bnil->cd(2); profphizcBnil->DrawCopy();
  c7Bnil->cd(0);
  c7Bnil->SaveAs("residualsProfilePhiZ-Bnil.eps");

  TCanvas* c8Bnil = new TCanvas("c8Bnil","Profiles y in phi, Bnil",800,300);
  c8Bnil->Divide(2,1);
  c8Bnil->cd(1); profphiyaBnil->DrawCopy();
  c8Bnil->cd(2); profphiycBnil->DrawCopy();
  c8Bnil->cd(0);
  c8Bnil->SaveAs("residualsProfilePhiY-Bnil.eps");

  TCanvas* c9Bnil = new TCanvas("c9Bnil", "Residuals Phi-Z 2D, Bnil",800,300);
  c9Bnil->Divide(2,1);
  c9Bnil->cd(1); pPhiZAResidualsHistBnil->DrawCopy("col");
  c9Bnil->cd(2); pPhiZCResidualsHistBnil->DrawCopy("col");
  c9Bnil->cd(0); 
  c9Bnil->SaveAs("residualsPhiZ2D-Bnil.eps");

  TCanvas* c1Bnil0 = new TCanvas("c1Bnil0", "Residuals Phi-Y 2D, Bnil",800,300);
  c1Bnil0->Divide(2,1);
  c1Bnil0->cd(1); pPhiYAResidualsHistBnil->DrawCopy("col");
  c1Bnil0->cd(2); pPhiYCResidualsHistBnil->DrawCopy("col");
  c1Bnil0->cd(0); 
  c1Bnil0->SaveAs("residualsPhiY2D-Bnil.eps");

  TCanvas* c11Bnil = new TCanvas("c11Bnil", "Residuals Pt-Z 2D, Bnil",800,300);
  c11Bnil->Divide(2,1);
  c11Bnil->cd(1); pPtZAResidualsHistBnil->DrawCopy("col");
  c11Bnil->cd(2); pPtZCResidualsHistBnil->DrawCopy("col");
  c11Bnil->cd(0);
  c11Bnil->SaveAs("residualsPtZ2D-Bnil.eps");

  TCanvas* c12Bnil = new TCanvas("c12Bnil", "Residuals Pt-Y 2D, Bnil",800,300);
  c12Bnil->Divide(2,1);
  c12Bnil->cd(1); pPtYAResidualsHistBnil->DrawCopy("col");
  c12Bnil->cd(2); pPtYCResidualsHistBnil->DrawCopy("col");
  c12Bnil->cd(0);
  c12Bnil->SaveAs("residualsPtY2D-Bnil.eps");

  TCanvas* c13Bnil = new TCanvas("c13Bnil","Profiles Pt-Z, Bnil",800,300);
  c13Bnil->Divide(2,1);
  c13Bnil->cd(1); profptzaBnil->DrawCopy();
  c13Bnil->cd(2); profptzcBnil->DrawCopy();
  c13Bnil->cd(0);
  c13Bnil->SaveAs("residualsProfilePtZ-Bnil.eps");

  TCanvas* c14Bnil = new TCanvas("c14Bnil","Profiles Pt-Y, Bnil",800,300);
  c14Bnil->Divide(2,1);
  c14Bnil->cd(1); profptyaBnil->DrawCopy();
  c14Bnil->cd(2); profptycBnil->DrawCopy();
  c14Bnil->cd(0);
  c14Bnil->SaveAs("residualsProfilePtY-Bnil.eps");

  TCanvas* c15Bnil = new TCanvas("c15Bnil", "Residuals low Pt-Z 2D, Bnil",800,300);
  c15Bnil->Divide(2,1);
  c15Bnil->cd(1); pLowPtZAResidualsHistBnil->DrawCopy("col");
  c15Bnil->cd(2); pLowPtZCResidualsHistBnil->DrawCopy("col");
  c15Bnil->cd(0);
  c15Bnil->SaveAs("residualsLowPtZ2D-Bnil.eps");

  TCanvas* c16Bnil = new TCanvas("c16Bnil", "Residuals in low PtY 2D, Bnil",800,300);
  c16Bnil->Divide(2,1);
  c16Bnil->cd(1); pLowPtYAResidualsHistBnil->DrawCopy("col");
  c16Bnil->cd(2); pLowPtYCResidualsHistBnil->DrawCopy("col");
  c16Bnil->cd(0);
  c16Bnil->SaveAs("residualsLowPtY2D-Bnil.eps");

  TCanvas* c17Bnil = new TCanvas("c17Bnil","Profiles low Pt-Z, Bnil",800,300);
  c17Bnil->Divide(2,1);
  c17Bnil->cd(1); proflowptzaBnil->DrawCopy();
  c17Bnil->cd(2); proflowptzcBnil->DrawCopy();
  c17Bnil->cd(0);
  c17Bnil->SaveAs("residualsProfileLowPtZ-Bnil.eps");

  TCanvas* c18Bnil = new TCanvas("c18Bnil","Profiles low Pt-Y, Bnil",800,300);
  c18Bnil->Divide(2,1);
  c18Bnil->cd(1); proflowptyaBnil->DrawCopy();
  c18Bnil->cd(2); proflowptycBnil->DrawCopy();
  c18Bnil->cd(0);
  c18Bnil->SaveAs("residualsProfileLowPtY-Bnil.eps");
  }
}

