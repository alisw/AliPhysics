checkMCmap(Int_t isector=0){
  TCanvas *c = new TCanvas();
  c->Divide(1,2);
  TH2I *h = new TH2I("map",Form("CTTM MAP sector(%i);istrip;half strip",isector),91,0,91,2,0,2);
  TH2I *h2 = new TH2I("map2",Form("Channel MAP sector(%i);istrip;half strip",isector),91,0,91,2,0,2);

  AliTOFTrigger tt;

  Int_t detind[5] = {isector,0,0,0,0};
  Int_t indexLTM[2];

  for(Int_t i=0;i<91;i++){
    for(Int_t j=0;j<2;j++){
      detind[4] = j*24;
      if(i<19){
	detind[1] = 0;
	detind[2] = i;
      }
      else if(i<38){
	detind[1] = 1;
	detind[2] = i-19;
      }
      else if(i<53){
	detind[1] = 2;
	detind[2] = i-38;
      }
      else if(i<72){
	detind[1] = 3;
	detind[2] = i-53;
      }
      else{
	detind[1] = 4;
	detind[2] = i-72;
      }

      tt.GetLTMIndex(detind,indexLTM);
      
      h->SetBinContent(i+1,j+1,indexLTM[0]);
      h2->SetBinContent(i+1,j+1,indexLTM[1]/2);


    }
  }

  TLine *l1 = new TLine(46,1,46,2);
  l1->SetLineWidth(5);
  TLine *l2 = new TLine(45,0,45,1);
  l2->SetLineWidth(5);

  c->cd(1);
  h->Draw("text");
  h->SetStats(0);
  l1->Draw("SAME");
  l2->Draw("SAME");
  c->cd(2);
  h2->Draw("text");
  h2->SetStats(0);
  l1->Draw("SAME");
  l2->Draw("SAME");

}
