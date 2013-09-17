void testAxis()
{
  const Int_t   nCenter = 56;
  const Int_t   nSat    = 10;
  const Int_t   nBins   = 2*nSat + nCenter;
  Float_t bins[nBins+1];
  
  for (Int_t i = 0; i < nSat;i++) {
    bins[i] = (i-nSat-.5) * 37.5;
    printf("bins[%2d]=%+6.2f\n", i, bins[i]);
  }
  for (Int_t i = nSat; i < nSat+nCenter+1; i++) {
    bins[i] = -nCenter + (i-nSat) * 2;
    printf("bins[%2d]=%+6.2f\n", i, bins[i]);
  }
  for (Int_t i = nSat+nCenter+1; i < 2*nSat+nCenter; i++) {
    bins[i] = (i-nSat-nCenter +.5) * 37.5;
    printf("bins[%2d]=%+6.2f\n", i, bins[i]);
  }
  bins[nBins] = (nSat + .5) * 37.5;
  
  // printf("\n");

  TH1* h = new TH1F("h", "h", nBins, bins);
  h->SetFillColor(kRed+1);
  h->SetFillStyle(3001);

  TAxis* a = h->GetXaxis();
  for (Int_t i = 1; i <= nBins; i++) 
    printf("%2d/%2d: %+6.2f - %+6.2f: %+6.2f\n", 
	   i,nBins,
	   a->GetBinLowEdge(i), 
	   a->GetBinUpEdge(i), 
	   a->GetBinCenter(i));
  
  for (Int_t i = 1; i <= nSat; i++) 
    h->SetBinContent(i, i);
  for (Int_t i = nSat+1; i <= nSat+nCenter; i++) 
    h->SetBinContent(i, nSat+nCenter/2 - TMath::Abs(i-nCenter/2-nSat));
  for (Int_t i = nSat+nCenter+1; i <= 2*nSat+nCenter; i++) 
    h->SetBinContent(i, (2*nSat+nCenter+1)-i);

  
  h->Draw();
}

  
    
