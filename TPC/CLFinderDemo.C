//Begin_Html
/*
<img src="../../../../pic/CLFinderDemo.gif">
*/
//End_Html
void CLFinderDemo()
{

  AliTPCDigitsArray *gDigitsArray= GetDigitsArray();
  AliTPCClustersArray *gExactClusters = GetExactClustersArray();
  
  
  AliH2F * his = gDigitsArray->LoadRow(0,0)->GenerHisto()->GetSubrange2d(440,470,60,90);
  
  AliClusterFinder cf;
  cf.SetThreshold(3);
  cf.SetNoise(1);
  Int_t index[3]={1,0,0};
  cf.SetDetectorIndex(index);
  cf.SetDetectorParam(gTPCParam);
  cf.GetHisto(his);
    
  cf.SetBFit(kFALSE);
  cf.GetMinuit()->SetPrintLevel(-1);
  cf.GetMinuit()->SetMaxIterations(100);
  cf.FindPeaks1()->GetEntriesFast();
    
    
  TCanvas *c0 = new TCanvas("signals", "Cluster finder drawing", 600,1100); 
  c0->Divide(1,3);
  c0->cd(1);
  cf.Draw("lego2");

  c0->cd(2);
  cf.Draw("cont1");
  cf.DrawCluster(4,5,5);  
  c0->cd(3);
  cf.DrawBorders("lego2");
    
}

