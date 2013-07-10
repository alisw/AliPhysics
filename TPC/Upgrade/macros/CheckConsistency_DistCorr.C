void makeComparisonTree(TString filename, TString addToName)
{
  TFile fn(filename.Data());
  gROOT->cd();
  fTPCCorrection=(AliTPCCorrection*)fn.Get("map");
  fn.Close();

  TString outFile=addToName;
  outFile.Append(".root");
  TTreeSRedirector *sred=new TTreeSRedirector(outFile.Data());
  
  Float_t dx[3]={0,0,0};
  
  for (Float_t iz=-240; iz<=240; iz+=20) {
    Short_t roc=(iz>=0)?0:18;
    for (Float_t ir=86; ir<250; ir+=20) {
      for (Float_t iphi=0; iphi<TMath::TwoPi(); iphi+=10*TMath::DegToRad()){
        Float_t x=ir*(Float_t)TMath::Cos(iphi);
        Float_t y=ir*(Float_t)TMath::Sin(iphi);
        Float_t x3[3]={x,y,iz};
        (*sred) << "t" <<
        "r="   << ir <<
        "phi=" << iphi <<
        "x=" << x <<
        "y=" << y <<
        "z=" << iz;
        
        //distortions
        fTPCCorrection->GetDistortion(x3,roc,dx);
        Float_t xd   = x+dx[0];
        Float_t yd   = y+dx[1];
        Float_t zd   = iz+dx[2];
        Float_t rd   = TMath::Sqrt(xd*xd+yd*yd);
        Float_t phid = TMath::ATan2(yd,xd);
        if (phid<0) phid+=TMath::TwoPi();
        (*sred) << "t" <<
        "xd="   << xd <<
        "yd="   << yd <<
        "zd="   << zd <<
        "rd="   << rd <<
        "phid=" << phid;
        
        // correct back distorted point
        Float_t xd3[3]={xd,yd,zd};
        
        fTPCCorrection->GetCorrection(xd3,roc,dx);
        Float_t xdc   = xd+dx[0];
        Float_t ydc   = yd+dx[1];
        Float_t zdc   = zd+dx[2];
        Float_t rdc   = TMath::Sqrt(xdc*xdc+ydc*ydc);
        Float_t phidc = TMath::ATan2(ydc,xdc);
        if (phidc<0) phidc+=TMath::TwoPi();
        
        (*sred)  << "t" <<
        "xdc="   << xdc <<
        "ydc="   << ydc <<
        "zdc="   << zdc <<
        "rdc="   << rdc <<
        "phidc=" << phidc;
        
        // write current point
        (*sred) << "t" <<
        "\n";
        
      }
    }
  }
  
  delete sred;
}

void makeAllComparisonTrees()
{
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz_precal.lookup.root","LUT_05");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps10_50kHz_precal.lookup.root","LUT_10");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps20_50kHz_precal.lookup.root","LUT_20");
}

TCanvas *GetCanvas(TString addToName);

void makeHistos(TString addToName) {
  TString fileName; //("test_");
  fileName.Append(addToName.Data());
  fileName.Append(".root");
  TFile f(fileName.Data());
  gROOT->cd();
  TTree *t=f.Get("t");
  gStyle->SetTitleX(0.18);
  gStyle->SetTitleW(1-.18-.1);

  t->SetMarkerStyle(20);
  t->SetMarkerSize(.8);

  TCanvas *c=0x0;
  c=GetCanvas(addToName+"_zRes");
  t->Draw("zdc-z:z:r","","colz");
  c->SaveAs(Form("%s_zRes.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_rRes");
  t->Draw("rdc-r:z:r","","colz");
  c->SaveAs(Form("%s_rRes.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_phiRes");
  t->Draw("phidc-phi:z:r","abs(phidc-phi)<1","colz");
  c->SaveAs(Form("%s_phiRes.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_rphiRes");
  t->Draw("(phidc*rdc)-(phi*r):z+(r-84)/(254-84)*18:r","abs(phidc-phi)<1","colz");
  c->SaveAs(Form("%s_rphiRes.png",addToName.Data()));

  f.Close();
}



void makeAllHistos() {
  makeHistos("LUT_05");
  makeHistos("LUT_10");
  makeHistos("LUT_20");

}

TCanvas *GetCanvas(TString addToName)
{
  TString cName(addToName);
  cName.Prepend("c_");
  TCanvas *c=gROOT->GetListOfCanvases()->FindObject(cName.Data());
  if (!c) c=new TCanvas(cName.Data(),addToName.Data());
  c->Clear();
  c->cd();
  return c;
}



