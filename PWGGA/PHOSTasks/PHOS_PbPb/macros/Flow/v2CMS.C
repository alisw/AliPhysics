TGraphErrors* v2CMS2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 6;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

  Double_t _x[] = {1.79, 2.24, 2.74, 3.46, 4.89, 6.75};
  Double_t _y[] = {0.161, 0.154, 0.169, 0.166, 0.131, 0.107};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.033, 0.031, 0.022, 0.033, 0.029, 0.034};


  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* v2CMS3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 6;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

  Double_t _x[] = {1.79, 2.24, 2.74, 3.46, 4.89, 6.75};
  Double_t _y[] = {0.178, 0.175, 0.192, 0.180, 0.137, 0.118};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.016, 0.017, 0.023, 0.024, 0.024, 0.018};


  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* v2CMS4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 6;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

  Double_t _x[] = {1.79, 2.24, 2.74, 3.46, 4.89, 6.75};
  Double_t _y[] = {0.192, 0.189, 0.194, 0.182, 0.153, 0.108};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.012, 0.013, 0.020, 0.017, 0.016, 0.016};


  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

