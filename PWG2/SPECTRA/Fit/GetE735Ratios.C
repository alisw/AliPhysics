TH1F* GetE735Ratios(Int_t ratio, Int_t energy) {

  Int_t nbins = 10;

  const Float_t binsKpi [] = {0.25,0.30, 0.35, 0.40,0.50,0.60,0.70,0.95,1.15,1.35,1.55};
  const Float_t binsPpi [] = {0.30, 0.35, 0.40,0.50,0.60,0.70,0.95,1.15,1.35,1.55,1.75};
  enum {kEn300,kEn540,kEn1000,kEn1800,kNenergy};

  
  TH1F * hKpi[kNenergy] ;
  TH1F * hPpi[kNenergy] ;
  
  const Int_t markers[] = {24,25,27,28};

  for(Int_t iener = 0; iener < kNenergy; iener++){
    hKpi[iener] = new TH1F(TString("hKpi")+long(iener),TString("hKpi")+long(iener),nbins,binsKpi);
    hPpi[iener] = new TH1F(TString("hPpi")+long(iener),TString("hPpi")+long(iener),nbins,binsPpi);
    hKpi[iener]->SetMarkerStyle(markers[iener]);
    hPpi[iener]->SetMarkerStyle(markers[iener]);    
  }
  


  // 300
  hKpi[kEn300]  ->SetBinContent(0 , 8.1  /100);
  hKpi[kEn300]  ->SetBinContent(1 , 7.5  /100);
  hKpi[kEn300]  ->SetBinContent(2 , 9.5  /100);
  hKpi[kEn300]  ->SetBinContent(3 , 16.1 /100);
  hKpi[kEn300]  ->SetBinContent(4 , 16.4 /100);
  hKpi[kEn300]  ->SetBinContent(5 , 17.9 /100);
  hKpi[kEn300]  ->SetBinContent(6 , 19.4 /100);
  hKpi[kEn300]  ->SetBinContent(7 , 23.8 /100);
  hKpi[kEn300]  ->SetBinContent(8 , 25.7 /100);
  hKpi[kEn300]  ->SetBinContent(9 , 0    /100);
  hKpi[kEn300]  ->SetBinError  (0 , 1.7/100);
  hKpi[kEn300]  ->SetBinError  (1 , 1.4/100);
  hKpi[kEn300]  ->SetBinError  (2 , 1.6/100);
  hKpi[kEn300]  ->SetBinError  (3 , 1.5/100);
  hKpi[kEn300]  ->SetBinError  (4 , 2.5/100);
  hKpi[kEn300]  ->SetBinError  (5 , 3.5/100);
  hKpi[kEn300]  ->SetBinError  (6 , 2.4/100);
  hKpi[kEn300]  ->SetBinError  (7 , 4.9/100);
  hKpi[kEn300]  ->SetBinError  (8 , 7.5/100);
  hKpi[kEn300]  ->SetBinError  (9 , 0  /100);

  hKpi[kEn540]  ->SetBinContent(0 , 3.3 /100);
  hKpi[kEn540]  ->SetBinContent(1 , 12.7/100);
  hKpi[kEn540]  ->SetBinContent(2 , 10.6/100);
  hKpi[kEn540]  ->SetBinContent(3 , 13.1/100);
  hKpi[kEn540]  ->SetBinContent(4 , 16.9/100);
  hKpi[kEn540]  ->SetBinContent(5 , 17.3/100);
  hKpi[kEn540]  ->SetBinContent(6 , 20.2/100);
  hKpi[kEn540]  ->SetBinContent(7 , 21.2/100);
  hKpi[kEn540]  ->SetBinContent(8 , 25.0/100);
  hKpi[kEn540]  ->SetBinContent(9 , 31.2/100);
  hKpi[kEn540]  ->SetBinError  (0 , 0.9/100);
  hKpi[kEn540]  ->SetBinError  (1 , 1.1/100);
  hKpi[kEn540]  ->SetBinError  (2 , 0.9/100);
  hKpi[kEn540]  ->SetBinError  (3 , 0.7/100);
  hKpi[kEn540]  ->SetBinError  (4 , 1.6/100);
  hKpi[kEn540]  ->SetBinError  (5 , 1.7/100);
  hKpi[kEn540]  ->SetBinError  (6 , 1.4/100);
  hKpi[kEn540]  ->SetBinError  (7 , 2.1/100);
  hKpi[kEn540]  ->SetBinError  (8 , 3.2/100);
  hKpi[kEn540]  ->SetBinError  (9 , 4.9/100);

  hKpi[kEn1000] ->SetBinContent(0 , 9.2 /100);
  hKpi[kEn1000] ->SetBinContent(1 , 10.4/100);
  hKpi[kEn1000] ->SetBinContent(2 , 9.3 /100);
  hKpi[kEn1000] ->SetBinContent(3 , 11.8/100);
  hKpi[kEn1000] ->SetBinContent(4 , 15.3/100);
  hKpi[kEn1000] ->SetBinContent(5 , 19.8/100);
  hKpi[kEn1000] ->SetBinContent(6 , 19.8/100);
  hKpi[kEn1000] ->SetBinContent(7 , 25.9/100);
  hKpi[kEn1000] ->SetBinContent(8 , 27.8/100);
  hKpi[kEn1000] ->SetBinContent(9 , 32.8/100);
  hKpi[kEn1000] ->SetBinError  (0 , 0.7/100);
  hKpi[kEn1000] ->SetBinError  (1 , 0.7/100);
  hKpi[kEn1000] ->SetBinError  (2 , 0.7/100);
  hKpi[kEn1000] ->SetBinError  (3 , 0.7/100);
  hKpi[kEn1000] ->SetBinError  (4 , 1.2/100);
  hKpi[kEn1000] ->SetBinError  (5 , 1.5/100);
  hKpi[kEn1000] ->SetBinError  (6 , 1.2/100);
  hKpi[kEn1000] ->SetBinError  (7 , 1.9/100);
  hKpi[kEn1000] ->SetBinError  (8 , 2.6/100);
  hKpi[kEn1000] ->SetBinError  (9 , 3.9/100);

  hKpi[kEn1800] ->SetBinContent(0 , 7.8 /100);
  hKpi[kEn1800] ->SetBinContent(1 , 9.1 /100);
  hKpi[kEn1800] ->SetBinContent(2 , 9.7 /100);
  hKpi[kEn1800] ->SetBinContent(3 , 13.9/100);
  hKpi[kEn1800] ->SetBinContent(4 , 16.2/100);
  hKpi[kEn1800] ->SetBinContent(5 , 18.1/100);
  hKpi[kEn1800] ->SetBinContent(6 , 21.6/100);
  hKpi[kEn1800] ->SetBinContent(7 , 24.6/100);
  hKpi[kEn1800] ->SetBinContent(8 , 25.8/100);
  hKpi[kEn1800] ->SetBinContent(9 , 25.2/100);
  hKpi[kEn1800] ->SetBinError  (0 , 0.3/100);
  hKpi[kEn1800] ->SetBinError  (1 , 0.3/100);
  hKpi[kEn1800] ->SetBinError  (2 , 0.3/100);
  hKpi[kEn1800] ->SetBinError  (3 , 0.3/100);
  hKpi[kEn1800] ->SetBinError  (4 , 0.6/100);
  hKpi[kEn1800] ->SetBinError  (5 , 0.7/100);
  hKpi[kEn1800] ->SetBinError  (6 , 0.7/100);
  hKpi[kEn1800] ->SetBinError  (7 , 0.8/100);
  hKpi[kEn1800] ->SetBinError  (8 , 0.9/100);
  hKpi[kEn1800] ->SetBinError  (9 , 1.0/100);

  if (ratio == 0) return hKpi[energy];
  else            return hPpi[energy];

}
