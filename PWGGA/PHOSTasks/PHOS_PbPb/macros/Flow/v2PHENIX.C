//source:   Phys. Rev. Lett. 105, 142301 (2010)
//http://www.phenix.bnl.gov/phenix/WWW/info/data/ppg110/fig1.txt

TGraphErrors* v2PHENIX010(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

   Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17};
   Double_t _y[] = {0.0523,0.0617,0.0697,0.0777,0.0746,0.0695,0.0612,0.0615,0.061,0.0391,0.0526,0.0467,0.0419,0.0055,0.0615,0.1435,0.1409};
Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t _yerr[] = {0.0031,0.0019,0.0018,0.0022,0.0029,0.0034,0.0042,0.0054,0.0057,0.0088,0.0126,0.0184,0.0261,0.0283,0.0512,0.1019,0.1372};

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* v2PHENIX1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

   Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17};
   Double_t _y[] = {0.0996,0.1208,0.132,0.1386,0.1373,0.1221,0.1174,0.1071,0.0949,0.0855,0.0887,0.0783,0.0991,0.0648,0.1224,0.1507,-0.1153};
Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t _yerr[] = {0.002,0.0012,0.0012,0.0014,0.0019,0.0022,0.003,0.004,0.0041,0.0066,0.01,0.0148,0.0211,0.0236,0.044,0.0753,0.1492};

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* v2PHENIX2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

   Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17};
   Double_t _y[] = {0.1346,0.1572,0.1704,0.175,0.1726,0.1646,0.1493,0.1333,0.1242,0.1081,0.1043,0.0796,0.0741,0.0577,0.0294,0.0161,0.1917};
Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t _yerr[] = {0.0018,0.0011,0.0011,0.0014,0.0018,0.0023,0.0031,0.0041,0.0039,0.0066,0.0103,0.0156,0.0222,0.0249,0.0455,0.0801,0.1424};

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2PHENIX3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

   Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17};
   Double_t _y[] = {0.1534,0.1788,0.1898,0.192,0.1823,0.1791,0.164,0.1584,0.1406,0.1252,0.1224,0.109,0.1163,0.1112,0.1214,0.1556,0.1111};
Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t _yerr[] = {0.0019,0.001,0.0011,0.0014,0.0019,0.0024,0.0032,0.0044,0.0044,0.0076,0.0122,0.0183,0.0278,0.0299,0.056,0.1043,0.194};

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2PHENIX4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

   Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17};
   Double_t _y[] = {0.1676,0.1862,0.1935,0.1989,0.1989,0.175,0.1736,0.1671,0.1473,0.1253,0.12,0.1098,0.1383,0.1269,0.054,0.0532,-0.0181};
Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t _yerr[] = {0.0022,0.0012,0.0013,0.0017,0.0022,0.0028,0.0039,0.0054,0.0058,0.0102,0.0166,0.0256,0.0389,0.0426,0.0879,0.1456,0.2846};

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* v2PHENIX5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

   Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17};
   Double_t _y[] = {0.165,0.1841,0.1894,0.1888,0.1825,0.1763,0.1838,0.1747,0.1734,0.1771,0.1153,0.1617,0.1986,0.1398,0.2379,0.08,0.287};
Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Double_t _yerr[] = {0.003,0.0016,0.0017,0.0022,0.0031,0.0042,0.0061,0.0085,0.0092,0.0161,0.0271,0.0419,0.0601,0.07,0.1227,0.2248,0.4645};

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

