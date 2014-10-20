Bool_t kSyst = kTRUE;
Bool_t kStat = kTRUE;

Double_t v0Cres[] = {0.448977,0.615462,0.712913,0.735905,0.697669,0.610942,0.480438,0.326472,0.185398};
Double_t v0Cres3[] = {0.296482,0.303404,0.294064,0.261364,0.220765,0.176069,0.112545,0.056252,0.009874};


void PrintPoints(TGraphErrors *g1,TGraphErrors *g2,const char *title="QM11 preliminary plots 10-20"){

  Int_t n1 = g1->GetN();
  Int_t n2 = g2->GetN();


  printf("\n%s\n",title);

  if(n1==n2){
    for(Int_t i=0;i<n1;i++){
      Float_t x1=g1->GetX()[i];
      Float_t x2=g2->GetX()[i];

      if(x1==x2){
	printf("%5.2f %e %e %e\n",x1,g1->GetY()[i],g1->GetEY()[i],g2->GetEY()[i]);
      }
    }
  }
  else{
    printf("The TGraphErrors have different number of points\nNothing done\n");
  }
}

void PrintPoints(TGraphAsymmErrors *g1,TGraphAsymmErrors *g2,const char *title="high pt flow paper 10-20"){

  Int_t n1 = g1->GetN();
  Int_t n2 = g2->GetN();


  printf("\n%s\n",title);

  if(n1==n2){
    for(Int_t i=0;i<n1;i++){
      Float_t x1=g1->GetX()[i];
      Float_t x2=g2->GetX()[i];

      if(x1==x2){
	printf("%5.2f %e %e %e %e\n",x1,g1->GetY()[i],g1->GetEYlow()[i],g2->GetEYlow()[i],g2->GetEYhigh()[i]);
      }
    }
  }
  else{
    printf("The TGraphErrors have different number of points\nNothing done\n");
  }
}



Double_t Eval(TGraphErrors *g,Double_t x){
    Int_t n = g->GetN();
    Int_t i = 0;
    while(x > g->GetX()[i]){
	i++;
    }
    if(i == 0) i = 1;

    Float_t distAll = g->GetX()[i] - g->GetX()[i-1];
    Float_t dist = x - g->GetX()[i-1];

    printf("%i\n",i);

    return g->GetY()[i-1] + (g->GetY()[i] - g->GetY()[i-1])/distAll*dist;
}

Double_t Eval(TGraphAsymmErrors *g,Double_t x){
    Int_t n = g->GetN();
    Int_t i = 0;
    while(x > g->GetX()[i]){
	i++;
    }
    if(i == 0) i = 1;

    Float_t distAll = g->GetX()[i] - g->GetX()[i-1];
    Float_t dist = x - g->GetX()[i-1];

    printf("%i\n",i);

    return g->GetY()[i-1] + (g->GetY()[i] - g->GetY()[i-1])/distAll*dist;
}

Double_t EvalError(TGraphErrors *g,Float_t x){
    Int_t n = g->GetN();
    Int_t i = 0;
    while(x < g->GetX()[i]){
	i++;
    }
    if(i == 0) i = 1;

    Float_t distAll = g->GetX()[i] - g->GetX()[i-1];
    Float_t dist = x - g->GetX()[i-1];

    return g->GetEY()[i-1] + (g->GetEY()[i] - g->GetEY()[i-1])/distAll*dist;
}

Double_t EvalError(TGraphAsymmErrors *g,Float_t x){
    Int_t n = g->GetN();
    Int_t i = 0;
    while(x < g->GetX()[i]){
	i++;
    }
    if(i == 0) i = 1;

    Float_t distAll = g->GetX()[i] - g->GetX()[i-1];
    Float_t dist = x - g->GetX()[i-1];

    return g->GetEY()[i-1] + (g->GetEYhigh()[i] - g->GetEYhigh()[i-1])/distAll*dist;
}

TGraph *GetV2(Int_t isp=0,Int_t cmin,Int_t cmax){
  // isp = 0(pi) 1(K) 2(pbar) 3(k0s) 4(lambda) 5(phi) 6(Xi) 7 (omega)
  if(cmin==0 && cmax==5){
    switch(isp){
    case 0:
      return v2Pion0005(1,20,4,36);
      break;
    case 1:
      return v2Kaon0005(1,20,5,36);
      break;
    case 2:
      return v2Antiproton0005(1,20,6,36);
      break;
    case 3:
      return Kv2_05_SPVZE();
      break;
    case 4:
      return Lv2_05_SPVZE();
      break;
    case 5:
      break;
    case 6:
      break;
    case 7:
      break;
    }
  }
  else if(cmin==5 && cmax==10){
    switch(isp){
    case 0:
      return v2Pion0510(1,20,4,36);
      break;
    case 1:
      return v2Kaon0510(1,20,5,36);
      break;
    case 2:
      return v2Antiproton0510(1,20,6,36);
      break;
    case 3:
      return Kv2_510_SPVZE();
      break;
    case 4:
      return Lv2_510_SPVZE();
      break;
    case 5:
      break;
    case 6:
      return v2Xi0510();
      break;
    case 7:
      return v2Omega0510();
      break;
    }
  }
  else if(cmin==10 && cmax==20){
    switch(isp){
    case 0:
      return v2Pion1020(1,20,4,36);
      break;
    case 1:
      return v2Kaon1020(1,20,5,36);
      break;
    case 2:
      return v2Antiproton1020(1,20,6,36);
      break;
    case 3:
      return Kv2_1020_SPVZE();
      break;
    case 4:
      return Lv2_1020_SPVZE();
      break;
    case 5:
      return v2EPPhi1020();
      break;
    case 6:
      return v2Xi1020();
      break;
    case 7:
      return v2Omega1020();
      break;
    }
  }
  else if(cmin==20 && cmax==30){
    switch(isp){
    case 0:
      return v2Pion2030(1,20,4,36);
      break;
    case 1:
      return v2Kaon2030(1,20,5,36);
      break;
    case 2:
      return v2Antiproton2030(1,20,6,36);
      break;
    case 3:
      return Kv2_2030_SPVZE();
      break;
    case 4:
      return Lv2_2030_SPVZE();
      break;
    case 5:
      return v2EPPhi2030();
      break;
    case 6:
      return v2Xi2030();
      break;
    case 7:
      return v2Omega2030();
      break;
    }
  }
  else if(cmin==30 && cmax==40){
    switch(isp){
    case 0:
      return v2Pion3040(1,20,4,36);
      break;
    case 1:
      return v2Kaon3040(1,20,5,36);
      break;
    case 2:
      return v2Antiproton3040(1,20,6,36);
      break;
    case 3:
      return Kv2_3040_SPVZE();
      break;
    case 4:
      return Lv2_3040_SPVZE();
      break;
    case 5:
      return v2EPPhi3040();
      break;
    case 6:
      return v2Xi3040();
      break;
    case 7:
      return v2Omega3040();
      break;
    }
  }
  else if(cmin==40 && cmax==50){
    switch(isp){
    case 0:
      return v2Pion4050(1,20,4,36);
      break;
    case 1:
      return v2Kaon4050(1,20,5,36);
      break;
    case 2:
      return v2Antiproton4050(1,20,6,36);
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      return v2Xi4050();
      break;
    case 7:
      return v2Omega4050();
      break;
    }
  }
  else if(cmin==50 && cmax==60){
    switch(isp){
    case 0:
      return v2Pion5060(1,20,4,36);
      break;
    case 1:
      return v2Kaon5060(1,20,5,36);
      break;
    case 2:
      return v2Antiproton5060(1,20,6,36);
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      return v2Xi5060();
      break;
    case 7:
      return v2Omega5060();
      break;
    }
  }
  else if(cmin==60 && cmax==70){
    switch(isp){
    case 0:
      return v2Pion6070(1,20,4,36);
      break;
    case 1:
      return v2Kaon6070(1,20,5,36);
      break;
    case 2:
      return v2Antiproton6070(1,20,6,36);
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      break;
    case 7:
      break;
    }
  }
  else if(cmin==70 && cmax==80){
    switch(isp){
    case 0:
      return v2Pion7080(1,20,4,36);
      break;
    case 1:
      return v2Kaon7080(1,20,5,36);
      break;
    case 2:
      return v2Antiproton7080(1,20,6,36);
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      break;
    case 7:
      break;
    }
  }

  return NULL;
}

// pion merged bins
TGraphAsymmErrors *v2Pion(Int_t cmin,Int_t cmax){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphAsymmErrors *gStat[10];
  TGraphAsymmErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gStat[ngr] = v2Pion0005(1,20,4,36);
      break;
    case 1:
      gStat[ngr] = v2Pion0510(1,20,4,36);
      break;
    case 2:
      gStat[ngr] = v2Pion1020(1,20,4,36);
      break;
    case 3:
      gStat[ngr] = v2Pion2030(1,20,4,36);
      break;
    case 4:
      gStat[ngr] = v2Pion3040(1,20,4,36);
      break;
    case 5:
      gStat[ngr] = v2Pion4050(1,20,4,36);
      break;
    case 6:
      gStat[ngr] = v2Pion5060(1,20,4,36);
      break;
    case 7:
      gStat[ngr] = v2Pion6070(1,20,4,36);
      break;
    case 8:
      gStat[ngr] = v2Pion7080(1,20,4,36);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gSyst[ngr] = v2Pion0005(1,20,4,36);
      break;
    case 1:
      gSyst[ngr] = v2Pion0510(1,20,4,36);
      break;
    case 2:
      gSyst[ngr] = v2Pion1020(1,20,4,36);
      break;
    case 3:
      gSyst[ngr] = v2Pion2030(1,20,4,36);
      break;
    case 4:
      gSyst[ngr] = v2Pion3040(1,20,4,36);
      break;
    case 5:
      gSyst[ngr] = v2Pion4050(1,20,4,36);
      break;
    case 6:
      gSyst[ngr] = v2Pion5060(1,20,4,36);
      break;
    case 7:
      gSyst[ngr] = v2Pion6070(1,20,4,36);
      break;
    case 8:
      gSyst[ngr] = v2Pion7080(1,20,4,36);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEYlow()[j] * gStat[i]->GetEYlow()[j] * v0Cres[imin+i]*v0Cres[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEYlow()[j]*gStat[i]->GetEYlow()[j] * weight * weight;
	eval2 += gSyst[i]->GetEYlow()[j] * weight;
	eval3 += gSyst[i]->GetEYhigh()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphAsymmErrors *gRes = new TGraphAsymmErrors(np,x,y,ex,ex,ey,ey2);
  gRes->SetMarkerStyle(20);
  gRes->SetMarkerColor(4);
  gRes->SetLineColor(4);

  return gRes;
}

// pion merged bins
TGraphAsymmErrors *v2PionAlex(Int_t cmin,Int_t cmax){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphAsymmErrors *gStat[10];
  TGraphAsymmErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gStat[ngr] = v2PionAlex0005(1,20);
      break;
    case 1:
      gStat[ngr] = v2PionAlex0510(1,20);
      break;
    case 2:
      gStat[ngr] = v2PionAlex1020(1,20);
      break;
    case 3:
      gStat[ngr] = v2PionAlex2030(1,20);
      break;
    case 4:
      gStat[ngr] = v2PionAlex3040(1,20);
      break;
    case 5:
      gStat[ngr] = v2PionAlex4050(1,20);
      break;
    case 6:
      gStat[ngr] = v2PionAlex5060(1,20);
      break;
    case 7:
      gStat[ngr] = v2PionAlex6070(1,20);
      break;
    case 8:
      gStat[ngr] = v2PionAlex7080(1,20);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gSyst[ngr] = v2PionAlex0005(1,20);
      break;
    case 1:
      gSyst[ngr] = v2PionAlex0510(1,20);
      break;
    case 2:
      gSyst[ngr] = v2PionAlex1020(1,20);
      break;
    case 3:
      gSyst[ngr] = v2PionAlex2030(1,20);
      break;
    case 4:
      gSyst[ngr] = v2PionAlex3040(1,20);
      break;
    case 5:
      gSyst[ngr] = v2PionAlex4050(1,20);
      break;
    case 6:
      gSyst[ngr] = v2PionAlex5060(1,20);
      break;
    case 7:
      gSyst[ngr] = v2PionAlex6070(1,20);
      break;
    case 8:
      gSyst[ngr] = v2PionAlex7080(1,20);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEYlow()[j] * gStat[i]->GetEYlow()[j] * v0Cres[imin+i]*v0Cres[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEYlow()[j]*gStat[i]->GetEYlow()[j] * weight * weight;
	eval2 += gSyst[i]->GetEYlow()[j] * weight;
	eval3 += gSyst[i]->GetEYhigh()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphAsymmErrors *gRes = new TGraphAsymmErrors(np,x,y,ex,ex,ey,ey2);
  gRes->SetMarkerStyle(20);
  gRes->SetMarkerColor(4);
  gRes->SetLineColor(4);

  return gRes;
}

TGraphErrors *v2PionQM11(Int_t cmin,Int_t cmax){
  if(cmin < 10) cmin = 10;
  if(cmax >60) cmax = 60;
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphErrors *gStat[10];
  TGraphErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      break;
    case 1:
      break;
    case 2:
      gStat[ngr] = v22_etagap10_1020_pion(1,20,3,31);
      break;
    case 3:
      gStat[ngr] = v22_etagap10_2030_pion(1,20,3,31);
      break;
    case 4:
      gStat[ngr] = v22_etagap10_3040_pion(1,20,3,31);
      break;
    case 5:
      gStat[ngr] = v22_etagap10_4050_pion(1,20,3,31);
      break;
    case 6:
      gStat[ngr] = v22_etagap10_5060_pion(1,20,3,31);
      break;
    case 7:
      break;
    case 8:
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
     switch(i){
    case 0:
      break;
    case 1:
      break;
    case 2:
      gSyst[ngr] = v22_etagap10_1020_pion(1,20,3,31);
      break;
    case 3:
      gSyst[ngr] = v22_etagap10_2030_pion(1,20,3,31);
      break;
    case 4:
      gSyst[ngr] = v22_etagap10_3040_pion(1,20,3,31);
      break;
    case 5:
      gSyst[ngr] = v22_etagap10_4050_pion(1,20,3,31);
      break;
    case 6:
      gSyst[ngr] = v22_etagap10_5060_pion(1,20,3,31);
      break;
    case 7:
      break;
    case 8:
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEY()[j] * gStat[i]->GetEY()[j];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEY()[j]*gStat[i]->GetEY()[j] * weight * weight;
	eval2 += gSyst[i]->GetEY()[j] * weight;
	eval3 += gSyst[i]->GetEY()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphErrors *gRes = new TGraphErrors(np,x,y,ex,ey);
  gRes->SetMarkerStyle(20);
  gRes->SetMarkerColor(4);
  gRes->SetLineColor(4);

  return gRes;
}

TGraphErrors *v3Pion(Int_t cmin,Int_t cmax,Int_t rebin=1){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphErrors *gStat[10];
  TGraphErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
	gStat[ngr] = v3Pion0005(1,20,4,36,rebin);
      break;
    case 1:
	gStat[ngr] = v3Pion0510(1,20,4,36,rebin);
      break;
    case 2:
	gStat[ngr] = v3Pion1020(1,20,4,36,rebin);
      break;
    case 3:
	gStat[ngr] = v3Pion2030(1,20,4,36,rebin);
      break;
    case 4:
	gStat[ngr] = v3Pion3040(1,20,4,36,rebin);
      break;
    case 5:
	gStat[ngr] = v3Pion4050(1,20,4,36,rebin);
      break;
    case 6:
	gStat[ngr] = v3Pion5060(1,20,4,36,rebin);
      break;
    case 7:
	gStat[ngr] = v3Pion6070(1,20,4,36,rebin);
      break;
    case 8:
	gStat[ngr] = v3Pion7080(1,20,4,36,rebin);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
	gSyst[ngr] = v3Pion0005(1,20,4,36,rebin);
      break;
    case 1:
	gSyst[ngr] = v3Pion0510(1,20,4,36,rebin);
      break;
    case 2:
	gSyst[ngr] = v3Pion1020(1,20,4,36,rebin);
      break;
    case 3:
	gSyst[ngr] = v3Pion2030(1,20,4,36,rebin);
      break;
    case 4:
	gSyst[ngr] = v3Pion3040(1,20,4,36,rebin);
      break;
    case 5:
	gSyst[ngr] = v3Pion4050(1,20,4,36,rebin);
      break;
    case 6:
	gSyst[ngr] = v3Pion5060(1,20,4,36,rebin);
      break;
    case 7:
	gSyst[ngr] = v3Pion6070(1,20,4,36,rebin);
      break;
    case 8:
	gSyst[ngr] = v3Pion7080(1,20,4,36,rebin);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t evalNR = 0;
    Double_t eval2 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEY()[j] * gStat[i]->GetEY()[j] * v0Cres3[imin+i]*v0Cres3[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEY()[j]*gStat[i]->GetEY()[j] * weight * weight;
	eval2 += gSyst[i]->GetEY()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
    }
    if(oldflag2){
      ey[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);

  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphErrors *gRes = new TGraphErrors(np,x,y,ex,ey);
  gRes->SetMarkerStyle(20);
  gRes->SetMarkerColor(4);
  gRes->SetLineColor(4);

  return gRes;
}

// antiproton merged bins
TGraphAsymmErrors *v2Antiproton(Int_t cmin,Int_t cmax){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphAsymmErrors *gStat[10];
  TGraphAsymmErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gStat[ngr] = v2Antiproton0005(1,20,6,40);
      break;
    case 1:
      gStat[ngr] = v2Antiproton0510(1,20,6,40);
      break;
    case 2:
      gStat[ngr] = v2Antiproton1020(1,20,6,40);
      break;
    case 3:
      gStat[ngr] = v2Antiproton2030(1,20,6,40);
      break;
    case 4:
      gStat[ngr] = v2Antiproton3040(1,20,6,40);
      break;
    case 5:
      gStat[ngr] = v2Antiproton4050(1,20,6,40);
      break;
    case 6:
      gStat[ngr] = v2Antiproton5060(1,20,6,40);
      break;
    case 7:
      gStat[ngr] = v2Antiproton6070(1,20,6,40);
      break;
    case 8:
      gStat[ngr] = v2Antiproton7080(1,20,6,40);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gSyst[ngr] = v2Antiproton0005(1,20,6,40);
      break;
    case 1:
      gSyst[ngr] = v2Antiproton0510(1,20,6,40);
      break;
    case 2:
      gSyst[ngr] = v2Antiproton1020(1,20,6,40);
      break;
    case 3:
      gSyst[ngr] = v2Antiproton2030(1,20,6,40);
      break;
    case 4:
      gSyst[ngr] = v2Antiproton3040(1,20,6,40);
      break;
    case 5:
      gSyst[ngr] = v2Antiproton4050(1,20,6,40);
      break;
    case 6:
      gSyst[ngr] = v2Antiproton5060(1,20,6,40);
      break;
    case 7:
      gSyst[ngr] = v2Antiproton6070(1,20,6,40);
      break;
    case 8:
      gSyst[ngr] = v2Antiproton7080(1,20,6,40);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEYlow()[j] * gStat[i]->GetEYlow()[j] * v0Cres[imin+i]*v0Cres[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEYlow()[j]*gStat[i]->GetEYlow()[j] * weight * weight;
	eval2 += gSyst[i]->GetEYlow()[j] * weight;
	eval3 += gSyst[i]->GetEYhigh()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphAsymmErrors *gRes = new TGraphAsymmErrors(np,x,y,ex,ex,ey,ey2);
  gRes->SetMarkerStyle(22);
  gRes->SetMarkerColor(2);
  gRes->SetLineColor(2);

  return gRes;
}

TGraphAsymmErrors *v2AntiprotonAlex(Int_t cmin,Int_t cmax){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphAsymmErrors *gStat[10];
  TGraphAsymmErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gStat[ngr] = v2AntiprotonAlex0005(1,20);
      break;
    case 1:
      gStat[ngr] = v2AntiprotonAlex0510(1,20);
      break;
    case 2:
      gStat[ngr] = v2AntiprotonAlex1020(1,20);
      break;
    case 3:
      gStat[ngr] = v2AntiprotonAlex2030(1,20);
      break;
    case 4:
      gStat[ngr] = v2AntiprotonAlex3040(1,20);
      break;
    case 5:
      gStat[ngr] = v2AntiprotonAlex4050(1,20);
      break;
    case 6:
      gStat[ngr] = v2AntiprotonAlex5060(1,20);
      break;
    case 7:
      gStat[ngr] = v2AntiprotonAlex6070(1,20);
      break;
    case 8:
      gStat[ngr] = v2AntiprotonAlex7080(1,20);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gSyst[ngr] = v2AntiprotonAlex0005(1,20);
      break;
    case 1:
      gSyst[ngr] = v2AntiprotonAlex0510(1,20);
      break;
    case 2:
      gSyst[ngr] = v2AntiprotonAlex1020(1,20);
      break;
    case 3:
      gSyst[ngr] = v2AntiprotonAlex2030(1,20);
      break;
    case 4:
      gSyst[ngr] = v2AntiprotonAlex3040(1,20);
      break;
    case 5:
      gSyst[ngr] = v2AntiprotonAlex4050(1,20);
      break;
    case 6:
      gSyst[ngr] = v2AntiprotonAlex5060(1,20);
      break;
    case 7:
      gSyst[ngr] = v2AntiprotonAlex6070(1,20);
      break;
    case 8:
      gSyst[ngr] = v2AntiprotonAlex7080(1,20);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEYlow()[j] * gStat[i]->GetEYlow()[j] * v0Cres[imin+i]*v0Cres[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEYlow()[j]*gStat[i]->GetEYlow()[j] * weight * weight;
	eval2 += gSyst[i]->GetEYlow()[j] * weight;
	eval3 += gSyst[i]->GetEYhigh()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphAsymmErrors *gRes = new TGraphAsymmErrors(np,x,y,ex,ex,ey,ey2);
  gRes->SetMarkerStyle(22);
  gRes->SetMarkerColor(2);
  gRes->SetLineColor(2);

  return gRes;
}

TGraphErrors *v2AntiprotonQM11(Int_t cmin,Int_t cmax){
  if(cmin < 10) cmin = 10;
  if(cmax >60) cmax = 60;
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphErrors *gStat[10];
  TGraphErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      break;
    case 1:
      break;
    case 2:
      gStat[ngr] = v22_etagap10_1020_antiproton(1,20,5,35);
      break;
    case 3:
      gStat[ngr] = v22_etagap10_2030_antiproton(1,20,5,35);
      break;
    case 4:
      gStat[ngr] = v22_etagap10_3040_antiproton(1,20,5,35);
      break;
    case 5:
      gStat[ngr] = v22_etagap10_4050_antiproton(1,20,5,35);
      break;
    case 6:
      gStat[ngr] = v22_etagap10_5060_antiproton(1,20,5,35);
      break;
    case 7:
      break;
    case 8:
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
     switch(i){
    case 0:
      break;
    case 1:
      break;
    case 2:
      gSyst[ngr] = v22_etagap10_1020_antiproton(1,20,5,35);
      break;
    case 3:
      gSyst[ngr] = v22_etagap10_2030_antiproton(1,20,5,35);
      break;
    case 4:
      gSyst[ngr] = v22_etagap10_3040_antiproton(1,20,5,35);
      break;
    case 5:
      gSyst[ngr] = v22_etagap10_4050_antiproton(1,20,5,35);
      break;
    case 6:
      gSyst[ngr] = v22_etagap10_5060_antiproton(1,20,5,35);
      break;
    case 7:
      break;
    case 8:
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEY()[j] * gStat[i]->GetEY()[j];
      if(weight > 0.0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEY()[j]*gStat[i]->GetEY()[j] * weight * weight;
	eval2 += gSyst[i]->GetEY()[j] * weight;
	eval3 += gSyst[i]->GetEY()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphErrors *gRes = new TGraphErrors(np,x,y,ex,ey);
  gRes->SetMarkerStyle(20);
  gRes->SetMarkerColor(4);
  gRes->SetLineColor(4);

  return gRes;
}

TGraphErrors *v3Antiproton(Int_t cmin,Int_t cmax,Int_t rebin=1){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphErrors *gStat[10];
  TGraphErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
	gStat[ngr] = v3Antiproton0005(2,22,6,36,rebin);
      break;
    case 1:
	gStat[ngr] = v3Antiproton0510(2,22,6,36,rebin);
      break;
    case 2:
	gStat[ngr] = v3Antiproton1020(2,22,6,36,rebin);
      break;
    case 3:
	gStat[ngr] = v3Antiproton2030(2,22,6,36,rebin);
      break;
    case 4:
	gStat[ngr] = v3Antiproton3040(2,22,6,36,rebin);
      break;
    case 5:
	gStat[ngr] = v3Antiproton4050(2,22,6,36,rebin);
      break;
    case 6:
	gStat[ngr] = v3Antiproton5060(2,22,6,36,rebin);
      break;
    case 7:
	gStat[ngr] = v3Antiproton6070(2,22,6,36,rebin);
      break;
    case 8:
	gStat[ngr] = v3Antiproton7080(2,22,6,36,rebin);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
	gSyst[ngr] = v3Antiproton0005(1,20,6,36,rebin);
      break;
    case 1:
	gSyst[ngr] = v3Antiproton0510(1,20,6,36,rebin);
      break;
    case 2:
	gSyst[ngr] = v3Antiproton1020(1,20,6,36,rebin);
      break;
    case 3:
	gSyst[ngr] = v3Antiproton2030(1,20,6,36,rebin);
      break;
    case 4:
	gSyst[ngr] = v3Antiproton3040(1,20,6,36,rebin);
      break;
    case 5:
	gSyst[ngr] = v3Antiproton4050(1,20,6,36,rebin);
      break;
    case 6:
	gSyst[ngr] = v3Antiproton5060(1,20,6,36,rebin);
      break;
    case 7:
	gSyst[ngr] = v3Antiproton6070(1,20,6,36,rebin);
      break;
    case 8:
	gSyst[ngr] = v3Antiproton7080(1,20,6,36,rebin);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEY()[j] * gStat[i]->GetEY()[j] * v0Cres3[imin+i]*v0Cres3[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEY()[j]*gStat[i]->GetEY()[j] * weight * weight;
	eval2 += gSyst[i]->GetEY()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
    }
    if(oldflag2){
      ey[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphErrors *gRes = new TGraphErrors(np,x,y,ex,ey);
  gRes->SetMarkerStyle(22);
  gRes->SetMarkerColor(2);
  gRes->SetLineColor(2);

  return gRes;
}

// kaon merged bins
TGraphAsymmErrors *v2Kaon(Int_t cmin,Int_t cmax){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphAsymmErrors *gStat[10];
  TGraphAsymmErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gStat[ngr] = v2Kaon0005(1,20,5,36);
      break;
    case 1:
      gStat[ngr] = v2Kaon0510(1,20,5,36);
      break;
    case 2:
      gStat[ngr] = v2Kaon1020(1,20,5,36);
      break;
    case 3:
      gStat[ngr] = v2Kaon2030(1,20,5,36);
      break;
    case 4:
      gStat[ngr] = v2Kaon3040(1,20,5,36);
      break;
    case 5:
      gStat[ngr] = v2Kaon4050(1,20,5,36);
      break;
    case 6:
      gStat[ngr] = v2Kaon5060(1,20,5,36);
      break;
    case 7:
      gStat[ngr] = v2Kaon6070(1,20,5,36);
      break;
    case 8:
      gStat[ngr] = v2Kaon7080(1,20,5,36);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      gSyst[ngr] = v2Kaon0005(1,20,5,36);
      break;
    case 1:
      gSyst[ngr] = v2Kaon0510(1,20,5,36);
      break;
    case 2:
      gSyst[ngr] = v2Kaon1020(1,20,5,36);
      break;
    case 3:
      gSyst[ngr] = v2Kaon2030(1,20,5,36);
      break;
    case 4:
      gSyst[ngr] = v2Kaon3040(1,20,5,36);
      break;
    case 5:
      gSyst[ngr] = v2Kaon4050(1,20,5,36);
      break;
    case 6:
      gSyst[ngr] = v2Kaon5060(1,20,5,36);
      break;
    case 7:
      gSyst[ngr] = v2Kaon6070(1,20,5,36);
      break;
    case 8:
      gSyst[ngr] = v2Kaon7080(1,20,5,36);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEYlow()[j] * gStat[i]->GetEYlow()[j] * v0Cres[imin+i]*v0Cres[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEYlow()[j]*gStat[i]->GetEYlow()[j] * weight * weight;
	eval2 += gSyst[i]->GetEYlow()[j] * weight;
	eval3 += gSyst[i]->GetEYhigh()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphAsymmErrors *gRes = new TGraphAsymmErrors(np,x,y,ex,ex,ey,ey2);
  gRes->SetMarkerStyle(21);
  gRes->SetMarkerColor(1);
  gRes->SetLineColor(1);

  return gRes;
}

TGraphErrors *v2KaonQM11(Int_t cmin,Int_t cmax){
  if(cmin < 10) cmin = 10;
  if(cmax >60) cmax = 60;
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphErrors *gStat[10];
  TGraphErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
      break;
    case 1:
      break;
    case 2:
      gStat[ngr] = v22_etagap10_1020_kaon(1,20,4,23);
      break;
    case 3:
      gStat[ngr] = v22_etagap10_2030_kaon(1,20,4,23);
      break;
    case 4:
      gStat[ngr] = v22_etagap10_3040_kaon(1,20,4,23);
      break;
    case 5:
      gStat[ngr] = v22_etagap10_4050_kaon(1,20,4,23);
      break;
    case 6:
      gStat[ngr] = v22_etagap10_5060_kaon(1,20,4,23);
      break;
    case 7:
      break;
    case 8:
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
     switch(i){
    case 0:
      break;
    case 1:
      break;
    case 2:
      gSyst[ngr] = v22_etagap10_1020_kaon(1,20,4,23);
      break;
    case 3:
      gSyst[ngr] = v22_etagap10_2030_kaon(1,20,4,23);
      break;
    case 4:
      gSyst[ngr] = v22_etagap10_3040_kaon(1,20,4,23);
      break;
    case 5:
      gSyst[ngr] = v22_etagap10_4050_kaon(1,20,4,23);
      break;
    case 6:
      gSyst[ngr] = v22_etagap10_5060_kaon(1,20,4,23);
      break;
    case 7:
      break;
    case 8:
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t eval3 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEY()[j] * gStat[i]->GetEY()[j];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEY()[j]*gStat[i]->GetEY()[j] * weight * weight;
	eval2 += gSyst[i]->GetEY()[j] * weight;
	eval3 += gSyst[i]->GetEY()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval3 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    ey2[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
      ey2[j] += eval3*eval3;
    }
    if(oldflag2){
      ey[j] += eval*eval;
      ey2[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
    ey2[j] = TMath::Sqrt(ey2[j]);
  }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphErrors *gRes = new TGraphErrors(np,x,y,ex,ey);
  gRes->SetMarkerStyle(20);
  gRes->SetMarkerColor(4);
  gRes->SetLineColor(4);

  return gRes;
}

TGraphErrors *v3Kaon(Int_t cmin,Int_t cmax,Int_t rebin=1){
  Int_t icentr[] = {0,5,10,20,30,40,50,60,70,80};
  Int_t imin=100;
  Int_t imax = -1;
  for(Int_t j=8;j>=0;j--)
    if(cmin < icentr[j+1]) imin = j;
  for(Int_t j=0;j<8;j++)
    if(cmax > icentr[j]) imax = j;
  
  Bool_t oldflag1 = kSyst;
  Bool_t oldflag2 = kStat;
  
  Int_t ngr = 0;
  TGraphErrors *gStat[10];
  TGraphErrors *gSyst[10];
  
  kSyst=0;
  kStat=1;
  
  char name[100];

  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
	gStat[ngr] = v3Kaon0005(1,20,5,36,rebin);
      break;
    case 1:
	gStat[ngr] = v3Kaon0510(1,20,5,36,rebin);
      break;
    case 2:
	gStat[ngr] = v3Kaon1020(1,20,5,36,rebin);
      break;
    case 3:
	gStat[ngr] = v3Kaon2030(1,20,5,36,rebin);
      break;
    case 4:
	gStat[ngr] = v3Kaon3040(1,20,5,36,rebin);
      break;
    case 5:
	gStat[ngr] = v3Kaon4050(1,20,5,36,rebin);
      break;
    case 6:
	gStat[ngr] = v3Kaon5060(1,20,5,36,rebin);
      break;
    case 7:
	gStat[ngr] = v3Kaon6070(1,20,5,36,rebin);
      break;
    case 8:
	gStat[ngr] = v3Kaon7080(1,20,5,36,rebin);
      break;
    }
    sprintf(name,"Stat%i",ngr);
    gStat[ngr]->SetName(name);
    ngr++;
  }
  
  kSyst=1;
  kStat=0;
  ngr = 0;
  
  for(Int_t i=imin;i <= imax;i++){
    switch(i){
    case 0:
	gSyst[ngr] = v3Kaon0005(1,20,5,36,rebin);
      break;
    case 1:
	gSyst[ngr] = v3Kaon0510(1,20,5,36,rebin);
      break;
    case 2:
	gSyst[ngr] = v3Kaon1020(1,20,5,36,rebin);
      break;
    case 3:
	gSyst[ngr] = v3Kaon2030(1,20,5,36,rebin);
      break;
    case 4:
	gSyst[ngr] = v3Kaon3040(1,20,5,36,rebin);
      break;
    case 5:
	gSyst[ngr] = v3Kaon4050(1,20,5,36,rebin);
      break;
    case 6:
	gSyst[ngr] = v3Kaon5060(1,20,5,36,rebin);
      break;
    case 7:
	gSyst[ngr] = v3Kaon6070(1,20,5,36,rebin);
      break;
    case 8:
	gSyst[ngr] = v3Kaon7080(1,20,5,36,rebin);
      break;
    }
    sprintf(name,"Syst%i",ngr);
    gSyst[ngr]->SetName(name);
    ngr++;
  }
  
  Double_t x[100],ex[100],y[100],ey[100],ey2[100];
  
  Int_t np = gSyst[0]->GetN();
  for(Int_t j=0;j<np;j++){
    Double_t val = 0;
    Double_t eval = 0;
    Double_t eval2 = 0;
    Double_t sumw = 0;
    for(Int_t i=0;i < ngr;i++){
      Double_t weight = gStat[i]->GetEY()[j] * gStat[i]->GetEY()[j] * v0Cres3[imin+i]*v0Cres3[imin+i];
      if(weight > 0){
	weight = 1./weight;
	val += gStat[i]->GetY()[j] * weight;
	sumw += weight;
	eval += gStat[i]->GetEY()[j]*gStat[i]->GetEY()[j] * weight * weight;
	eval2 += gSyst[i]->GetEY()[j] * weight;
      }
    }
    val /= sumw;
    eval2 /= sumw;
    eval = TMath::Sqrt(eval)/sumw;
    
    x[j] = gStat[0]->GetX()[j];
    ex[j] = 0;
    y[j] = val;
    ey[j] = 0;
    if(oldflag1){
      ey[j] += eval2*eval2;
    }
    if(oldflag2){
      ey[j] += eval*eval;
    }
    else{
	ex[j] = 0.05;
    }
    ey[j] = TMath::Sqrt(ey[j]);
   }
  
  kSyst = oldflag1;
  kStat = oldflag2;

  TGraphErrors *gRes = new TGraphErrors(np,x,y,ex,ey);
  gRes->SetMarkerStyle(21);
  gRes->SetMarkerColor(1);
  gRes->SetLineColor(1);

  return gRes;
}


// ***********************  K0 and Lambda flow ************************************************************************

void k0lambda() {
  Kv2_3040_QC2()->Draw("APS");
  K2030_QC2()->Draw("PSAME");
  K1020_QC2()->Draw("PSAME");
  Kv2_510_QC2()->Draw("PSAME");
  Kv2_05_QC2()->Draw("PSAME");

  new TCanvas;
  Lv2_3040_QC2()->Draw("APS");
  L2030_QC2()->Draw("PSAME");
  L1020_QC2()->Draw("PSAME");
  Lv2_510_QC2()->Draw("PSAME");
  Lv2_05_QC2()->Draw("PSAME");
}

TGraphErrors* Kv2_05_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.004252, 0.005188, 0.008694, 0.021937, 0.027455, 0.034104, 0.039625, 0.044746, 0.048867, 0.050852, 0.056576, 0.057511, 0.061194, 0.064164, 0.043417, 0.034294, 0.045038, 0.055057, 0.054040};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.007064, 0.002132, 0.001111, 0.000885, 0.000831, 0.000879, 0.001009, 0.001196, 0.001423, 0.001768, 0.002138, 0.002071, 0.003061, 0.004600, 0.006617, 0.007092, 0.013762, 0.017200, 0.034379};
  Double_t _ysys[] = {0.000407, 0.000469, 0.000679, 0.001100, 0.001374, 0.001748, 0.001981, 0.002328, 0.002508, 0.002882, 0.002882, 0.003027, 0.003337, 0.003455, 0.003552, 0.002006, 0.002474, 0.009195, 0.010798};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_510_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.007000, 0.011162, 0.021922, 0.035766, 0.049923, 0.061198, 0.072225, 0.079441, 0.087496, 0.091636, 0.097925, 0.103443, 0.098518, 0.112176, 0.095534, 0.085755, 0.081923, 0.056758, 0.037123};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004525, 0.001531, 0.000857, 0.000685, 0.000641, 0.000694, 0.000797, 0.000935, 0.001034, 0.001370, 0.001699, 0.001582, 0.002280, 0.003377, 0.004748, 0.004655, 0.009180, 0.012872, 0.023541};
  Double_t _ysys[] = {0.001626, 0.000605, 0.000954, 0.001506, 0.002063, 0.002470, 0.002897, 0.003226, 0.003512, 0.003752, 0.004013, 0.004140, 0.004256, 0.004645, 0.004638, 0.003733, 0.008835, 0.005530, 0.006966};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_1020_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.006448, 0.017702, 0.036487, 0.056327, 0.073865, 0.093749, 0.106839, 0.118791, 0.128331, 0.138452, 0.143128, 0.150689, 0.152251, 0.144184, 0.136654, 0.134234, 0.122976, 0.090678, 0.062086};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004170, 0.001288, 0.000779, 0.000640, 0.000620, 0.000674, 0.000751, 0.000886, 0.001046, 0.001270, 0.001556, 0.001467, 0.002149, 0.003077, 0.004349, 0.004403, 0.007759, 0.011108, 0.026570};
  Double_t _ysys[] = {0.000588, 0.000584, 0.001122, 0.001690, 0.002222, 0.002921, 0.003219, 0.003614, 0.003882, 0.004232, 0.004320, 0.004618, 0.004829, 0.004410, 0.004243, 0.004993, 0.003940, 0.004961, 0.019516};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_2030_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.010006, 0.024794, 0.052177, 0.078209, 0.102262, 0.123894, 0.140711, 0.156015, 0.169287, 0.179601, 0.183815, 0.192834, 0.198249, 0.191741, 0.176891, 0.166178, 0.138895, 0.110028, 0.153199};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004230, 0.001405, 0.000881, 0.000736, 0.000726, 0.000765, 0.000863, 0.000990, 0.001184, 0.001466, 0.001765, 0.001639, 0.002407, 0.003441, 0.004740, 0.004662, 0.008490, 0.011135, 0.021290};
  Double_t _ysys[] = {0.000948, 0.001279, 0.001622, 0.002386, 0.003122, 0.003718, 0.004253, 0.004744, 0.005222, 0.005395, 0.005579, 0.005818, 0.006045, 0.005821, 0.005667, 0.005197, 0.005507, 0.010924, 0.021326};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_3040_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.016004, 0.032272, 0.063662, 0.094458, 0.122230, 0.147693, 0.165353, 0.182370, 0.192551, 0.203017, 0.209774, 0.215077, 0.205180, 0.207687, 0.193207, 0.177072, 0.138215, 0.095648, 0.094620};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003860, 0.001294, 0.000829, 0.000700, 0.000694, 0.000746, 0.000854, 0.000996, 0.001182, 0.001428, 0.001702, 0.001573, 0.002237, 0.003101, 0.004298, 0.004143, 0.007344, 0.010269, 0.019498};
  Double_t _ysys[] = {0.001196, 0.001172, 0.001913, 0.002838, 0.003787, 0.004457, 0.004986, 0.005593, 0.005961, 0.006311, 0.006323, 0.006607, 0.006218, 0.006364, 0.005809, 0.005347, 0.004869, 0.006404, 0.011155};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_05_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.007244, 0.000766, 0.002918, 0.010512, 0.020161, 0.028609, 0.035194, 0.042263, 0.048493, 0.054680, 0.063401, 0.062718, 0.074361, 0.069395, 0.054816, 0.035160, 0.065080};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.007859, 0.003415, 0.001886, 0.001502, 0.001318, 0.001249, 0.001266, 0.001380, 0.001507, 0.001271, 0.001736, 0.002462, 0.003518, 0.004229, 0.011940, 0.025584, 0.064570};
  Double_t _ysys[] = {0.010257, 0.003898, 0.000858, 0.000705, 0.001076, 0.001437, 0.001761, 0.002292, 0.002449, 0.002735, 0.003459, 0.003155, 0.003993, 0.003532, 0.004840, 0.013337, 0.037902};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_510_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.000448, -0.002030, 0.007388, 0.022971, 0.036594, 0.048758, 0.064942, 0.076033, 0.084732, 0.100738, 0.118755, 0.127093, 0.137373, 0.131762, 0.119606, 0.081458, 0.061535};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.005760, 0.002502, 0.001393, 0.001129, 0.001021, 0.000977, 0.000998, 0.001047, 0.001164, 0.000988, 0.001377, 0.001973, 0.002857, 0.003306, 0.008045, 0.016394, 0.044974};
  Double_t _ysys[] = {0.001294, 0.002166, 0.001195, 0.001337, 0.001570, 0.001979, 0.002601, 0.003128, 0.003399, 0.004032, 0.004824, 0.005133, 0.005533, 0.005522, 0.004800, 0.016161, 0.022762};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_1020_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.003309, 0.008249, 0.024390, 0.038326, 0.057884, 0.077701, 0.098417, 0.115582, 0.130611, 0.152156, 0.176539, 0.188140, 0.202555, 0.200892, 0.178190, 0.158896, 0.154960};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.005072, 0.001860, 0.001218, 0.000996, 0.000922, 0.000912, 0.000944, 0.001027, 0.001121, 0.000962, 0.001303, 0.001840, 0.002665, 0.002998, 0.007228, 0.013583, 0.037967};
  Double_t _ysys[] = {0.000509, 0.000711, 0.001270, 0.001621, 0.001764, 0.002576, 0.002954, 0.003583, 0.003920, 0.004565, 0.005324, 0.005754, 0.006182, 0.006043, 0.006355, 0.006942, 0.025241};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_2030_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.005534, 0.016922, 0.040005, 0.063529, 0.089775, 0.113536, 0.140542, 0.157994, 0.179955, 0.205746, 0.228441, 0.245757, 0.250396, 0.253334, 0.212181, 0.144409, 0.127257};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004080, 0.001874, 0.001310, 0.001097, 0.001033, 0.001037, 0.001084, 0.001178, 0.001300, 0.001124, 0.001515, 0.002140, 0.003040, 0.003451, 0.007882, 0.014386, 0.040230};
  Double_t _ysys[] = {0.002796, 0.002731, 0.002055, 0.002177, 0.002961, 0.003688, 0.004239, 0.004790, 0.005451, 0.006265, 0.006935, 0.007398, 0.007591, 0.007668, 0.007851, 0.004401, 0.035996};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_3040_SPVZE(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.008101, 0.028716, 0.056285, 0.085861, 0.113487, 0.142217, 0.168666, 0.190902, 0.214881, 0.236344, 0.259644, 0.277207, 0.284432, 0.269306, 0.230249, 0.170618, 0.111317};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003465, 0.001650, 0.001186, 0.001001, 0.001016, 0.000989, 0.001040, 0.001151, 0.001294, 0.001105, 0.001520, 0.002137, 0.002988, 0.003283, 0.007364, 0.013282, 0.033430};
  Double_t _ysys[] = {0.002027, 0.000896, 0.001954, 0.002742, 0.003461, 0.004299, 0.005087, 0.005752, 0.006511, 0.007091, 0.007937, 0.008339, 0.008937, 0.008115, 0.007644, 0.007032, 0.022812};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* Kv2_05_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.001041, 0.001124, 0.008947, 0.020038, 0.029261, 0.038507, 0.045613, 0.052333, 0.056686, 0.062340, 0.065821, 0.065011, 0.070707, 0.072572, 0.063348, 0.069030, 0.051120, 0.056467, 0.023178};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004180, 0.001479, 0.000893, 0.000735, 0.000683, 0.000732, 0.000832, 0.000975, 0.001149, 0.001386, 0.001683, 0.001578, 0.002352, 0.003466, 0.004909, 0.004808, 0.008910, 0.011868, 0.022452};
  Double_t _ysys[] = {0.000454, 0.000997, 0.000792, 0.001007, 0.001464, 0.001947, 0.002286, 0.002674, 0.002867, 0.003308, 0.003299, 0.003312, 0.003542, 0.003664, 0.003686, 0.003527, 0.003099, 0.003576, 0.009997};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_510_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.003497, 0.010205, 0.020795, 0.036825, 0.051873, 0.063821, 0.075186, 0.085009, 0.093230, 0.099266, 0.106385, 0.110677, 0.107768, 0.111739, 0.108875, 0.099813, 0.097247, 0.090101, 0.047768};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003794, 0.001375, 0.000857, 0.000698, 0.000671, 0.000712, 0.000790, 0.000930, 0.001042, 0.001321, 0.001621, 0.001466, 0.002187, 0.003218, 0.004481, 0.004463, 0.008598, 0.010870, 0.018896};
  Double_t _ysys[] = {0.001352, 0.000521, 0.000942, 0.001586, 0.002090, 0.002557, 0.003016, 0.003424, 0.003743, 0.003981, 0.004308, 0.004506, 0.004379, 0.004541, 0.004788, 0.004033, 0.005403, 0.003820, 0.003682};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_1020_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.007051, 0.016100, 0.037009, 0.057360, 0.076636, 0.096443, 0.110868, 0.124309, 0.134059, 0.143032, 0.148407, 0.156759, 0.161398, 0.156367, 0.148881, 0.142463, 0.135274, 0.108699, 0.109107};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003485, 0.001355, 0.000849, 0.000704, 0.000690, 0.000741, 0.000824, 0.000949, 0.001078, 0.001346, 0.001635, 0.001515, 0.002217, 0.003181, 0.004494, 0.004505, 0.008140, 0.011124, 0.020627};
  Double_t _ysys[] = {0.001480, 0.000588, 0.001159, 0.001734, 0.002305, 0.002927, 0.003347, 0.003831, 0.004036, 0.004330, 0.004453, 0.004743, 0.005010, 0.005319, 0.004618, 0.004824, 0.004091, 0.003268, 0.013869};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_2030_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.002309, 0.024107, 0.051778, 0.079557, 0.105023, 0.127597, 0.144243, 0.161524, 0.174976, 0.186361, 0.189533, 0.197314, 0.204868, 0.196847, 0.192251, 0.177620, 0.169643, 0.122656, 0.139651};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003981, 0.001516, 0.000978, 0.000819, 0.000813, 0.000852, 0.000954, 0.001078, 0.001302, 0.001583, 0.001903, 0.001755, 0.002527, 0.003632, 0.005061, 0.004993, 0.008796, 0.013641, 0.022930};
  Double_t _ysys[] = {0.000531, 0.000796, 0.001567, 0.002411, 0.003173, 0.003837, 0.004342, 0.004852, 0.005271, 0.005592, 0.005688, 0.005956, 0.006215, 0.006195, 0.005855, 0.005481, 0.005122, 0.005476, 0.011500};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_3040_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.014508, 0.034328, 0.065797, 0.096833, 0.125000, 0.150947, 0.169894, 0.186345, 0.199821, 0.210376, 0.220423, 0.222461, 0.217520, 0.221469, 0.207910, 0.196929, 0.152205, 0.135415, 0.145111};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003475, 0.001344, 0.000838, 0.000757, 0.000755, 0.000803, 0.000902, 0.001042, 0.001223, 0.001458, 0.001762, 0.001633, 0.002270, 0.003212, 0.004436, 0.004225, 0.007481, 0.009788, 0.019137};
  Double_t _ysys[] = {0.000908, 0.001215, 0.001981, 0.002925, 0.003915, 0.004574, 0.005100, 0.005639, 0.006018, 0.006414, 0.006650, 0.006688, 0.006582, 0.006824, 0.006319, 0.006235, 0.004788, 0.006151, 0.008116};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_05_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.012915, -0.017239, -0.004504, 0.006377, 0.016216, 0.027976, 0.037605, 0.043651, 0.052304, 0.061816, 0.072772, 0.080232, 0.087232, 0.090775, 0.075730, 0.093122, 0.142425};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.005322, 0.002200, 0.001487, 0.001218, 0.001104, 0.001044, 0.001061, 0.001133, 0.001215, 0.001046, 0.001348, 0.001905, 0.002767, 0.003223, 0.008212, 0.015866, 0.041348};
  Double_t _ysys[] = {0.001355, 0.001125, 0.000820, 0.000597, 0.001260, 0.001425, 0.002446, 0.002519, 0.003251, 0.003270, 0.004089, 0.004892, 0.004666, 0.004660, 0.003883, 0.015552, 0.020457};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_510_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.014544, -0.007924, 0.004837, 0.018192, 0.033375, 0.049417, 0.065286, 0.076631, 0.089073, 0.104481, 0.125041, 0.133032, 0.146169, 0.147679, 0.136833, 0.119569, 0.125120};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004719, 0.001979, 0.001368, 0.001133, 0.001031, 0.000997, 0.001002, 0.001055, 0.001146, 0.000975, 0.001325, 0.001863, 0.002759, 0.003156, 0.007594, 0.013760, 0.037775};
  Double_t _ysys[] = {0.002082, 0.002344, 0.000918, 0.000859, 0.001522, 0.002101, 0.002649, 0.003080, 0.003580, 0.004204, 0.005062, 0.005391, 0.005913, 0.005966, 0.007519, 0.005709, 0.009710};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_1020_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.010600, 0.001161, 0.020148, 0.037486, 0.055986, 0.077790, 0.099118, 0.117230, 0.133895, 0.156136, 0.181408, 0.195720, 0.205453, 0.207396, 0.189198, 0.182420, 0.141741};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004043, 0.001840, 0.001298, 0.001073, 0.000994, 0.000989, 0.001018, 0.001092, 0.001192, 0.001025, 0.001357, 0.001917, 0.002753, 0.003091, 0.007417, 0.014128, 0.036993};
  Double_t _ysys[] = {0.003305, 0.000447, 0.000611, 0.001148, 0.001861, 0.002359, 0.003017, 0.003581, 0.004034, 0.004722, 0.005444, 0.005884, 0.006169, 0.006594, 0.005812, 0.005919, 0.007949};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_2030_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.001192, 0.013454, 0.038257, 0.061408, 0.088608, 0.115157, 0.141186, 0.160925, 0.184930, 0.211126, 0.234058, 0.251488, 0.259255, 0.264109, 0.231330, 0.190554, 0.203114};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004142, 0.001977, 0.001425, 0.001193, 0.001140, 0.001137, 0.001178, 0.001275, 0.001414, 0.001197, 0.001612, 0.002259, 0.003224, 0.003610, 0.008382, 0.014979, 0.039175};
  Double_t _ysys[] = {0.003300, 0.000466, 0.001315, 0.001868, 0.002704, 0.003535, 0.004245, 0.004835, 0.005564, 0.006438, 0.007037, 0.007627, 0.007836, 0.007952, 0.007945, 0.006502, 0.021830};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_3040_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.001546, 0.027875, 0.057020, 0.084991, 0.114717, 0.144647, 0.171366, 0.195803, 0.217127, 0.242040, 0.268601, 0.285450, 0.292867, 0.283628, 0.252365, 0.208233, 0.184286};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003424, 0.001730, 0.001253, 0.001064, 0.001076, 0.001074, 0.001095, 0.001202, 0.001344, 0.001144, 0.001557, 0.002167, 0.003051, 0.003422, 0.007545, 0.013201, 0.030984};
  Double_t _ysys[] = {0.001534, 0.000858, 0.001718, 0.002558, 0.003445, 0.004341, 0.005153, 0.005881, 0.006542, 0.007262, 0.008076, 0.008579, 0.009650, 0.008542, 0.007900, 0.007809, 0.014849};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_510_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.024448, 0.006490, 0.015154, 0.028703, 0.041548, 0.051786, 0.059059, 0.063980, 0.069818, 0.079014, 0.072393, 0.081744, 0.072246, 0.083706, 0.064711, 0.086149, 0.101979, -0.020145, -0.030268};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.014226, 0.005388, 0.003665, 0.003072, 0.003015, 0.003076, 0.003523, 0.004001, 0.004486, 0.005411, 0.006457, 0.005767, 0.007998, 0.012337, 0.016395, 0.016698, 0.033892, 0.039218, 0.068058};
  Double_t _ysys[] = {0.009051, 0.001662, 0.000632, 0.001870, 0.001916, 0.002614, 0.002448, 0.002588, 0.003495, 0.005612, 0.002927, 0.004893, 0.006262, 0.009444, 0.012260, 0.005168, 0.021448, 0.025203, 0.039657};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_1020_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.010304, 0.016204, 0.033282, 0.050568, 0.068130, 0.082275, 0.098983, 0.107410, 0.118314, 0.122915, 0.127885, 0.135673, 0.133027, 0.125030, 0.119996, 0.114883, 0.103795, 0.056606, 0.073532};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.008246, 0.003348, 0.002248, 0.001928, 0.001887, 0.001998, 0.002171, 0.002489, 0.002762, 0.003426, 0.004034, 0.003694, 0.005368, 0.007697, 0.010785, 0.010791, 0.019094, 0.026001, 0.046208};
  Double_t _ysys[] = {0.003572, 0.000550, 0.001082, 0.001520, 0.002068, 0.003303, 0.003411, 0.003818, 0.003650, 0.003819, 0.003868, 0.004158, 0.005506, 0.003772, 0.003747, 0.003794, 0.003375, 0.008648, 0.025117};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_2030_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.004850, 0.027542, 0.046313, 0.075131, 0.094502, 0.115978, 0.130658, 0.141709, 0.153858, 0.163352, 0.166724, 0.172142, 0.171170, 0.154985, 0.148088, 0.141224, 0.112097, 0.083910, 0.066892};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.007987, 0.003216, 0.002194, 0.001899, 0.001856, 0.001936, 0.002187, 0.002423, 0.002887, 0.003383, 0.004106, 0.003696, 0.005245, 0.007658, 0.010429, 0.010345, 0.018054, 0.024491, 0.045947};
  Double_t _ysys[] = {0.000285, 0.000849, 0.001456, 0.002505, 0.003076, 0.003492, 0.003941, 0.004288, 0.004863, 0.005408, 0.005098, 0.005306, 0.005200, 0.004955, 0.005220, 0.004885, 0.004737, 0.012374, 0.012547};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_3040_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.006512, 0.034943, 0.062114, 0.088630, 0.109348, 0.133985, 0.149426, 0.161983, 0.172958, 0.178089, 0.179392, 0.185962, 0.166788, 0.179056, 0.169265, 0.143349, 0.121933, 0.108497, 0.112296};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.007580, 0.003105, 0.002090, 0.001853, 0.001904, 0.001991, 0.002207, 0.002507, 0.002931, 0.003438, 0.004063, 0.003738, 0.005244, 0.007411, 0.009942, 0.009662, 0.016750, 0.021741, 0.039454};
  Double_t _ysys[] = {0.002840, 0.001382, 0.001869, 0.002727, 0.003337, 0.004088, 0.004513, 0.004862, 0.005472, 0.005347, 0.005438, 0.005909, 0.005198, 0.006915, 0.008132, 0.004807, 0.003987, 0.010813, 0.004271};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_510_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.012776, -0.006532, 0.002407, 0.013908, 0.034298, 0.043741, 0.053164, 0.061653, 0.066822, 0.083646, 0.099757, 0.113607, 0.130599, 0.088641, 0.108932, 0.091642, 0.003873};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.018501, 0.008241, 0.005990, 0.004994, 0.004546, 0.004265, 0.004395, 0.004456, 0.004907, 0.003900, 0.005136, 0.007247, 0.010837, 0.011661, 0.029174, 0.053332, 0.123189};
  Double_t _ysys[] = {0.002945, 0.007251, 0.002181, 0.001572, 0.001461, 0.001965, 0.002908, 0.002978, 0.004778, 0.004521, 0.004007, 0.007012, 0.006891, 0.013885, 0.014025, 0.004857, 0.022835};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_1020_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.018354, 0.002501, 0.019160, 0.035782, 0.053545, 0.072513, 0.088244, 0.107401, 0.116970, 0.139628, 0.158563, 0.168384, 0.175704, 0.168032, 0.176285, 0.129312, 0.232092};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.009524, 0.004461, 0.003259, 0.002770, 0.002591, 0.002571, 0.002613, 0.002780, 0.002965, 0.002566, 0.003380, 0.004681, 0.006765, 0.007527, 0.017929, 0.032354, 0.089079};
  Double_t _ysys[] = {0.005604, 0.000546, 0.002036, 0.001541, 0.001710, 0.002222, 0.002798, 0.004193, 0.003757, 0.004209, 0.005259, 0.005308, 0.006049, 0.006522, 0.010668, 0.008480, 0.015839};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_2030_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.006141, 0.016013, 0.033080, 0.061927, 0.081133, 0.103266, 0.128730, 0.147669, 0.165525, 0.186983, 0.207332, 0.221861, 0.227361, 0.236137, 0.179255, 0.127896, 0.158990};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.008335, 0.004076, 0.003074, 0.002625, 0.002498, 0.002492, 0.002622, 0.002774, 0.003122, 0.002550, 0.003460, 0.004807, 0.006813, 0.007646, 0.017565, 0.032766, 0.074537};
  Double_t _ysys[] = {0.001704, 0.000755, 0.001575, 0.002334, 0.002515, 0.003214, 0.004683, 0.004838, 0.005278, 0.005797, 0.006703, 0.008488, 0.008067, 0.009027, 0.005602, 0.011053, 0.009362};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv2_3040_QC4(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.012902, 0.025954, 0.049970, 0.076305, 0.105496, 0.130038, 0.153960, 0.176067, 0.194724, 0.213066, 0.226722, 0.241776, 0.237693, 0.227219, 0.182900, 0.143968, -0.042035};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.007691, 0.003874, 0.003023, 0.002608, 0.002564, 0.002662, 0.002693, 0.002911, 0.003235, 0.002711, 0.003660, 0.005053, 0.007198, 0.007945, 0.017207, 0.031596, 0.071902};
  Double_t _ysys[] = {0.003312, 0.003793, 0.002313, 0.002394, 0.004059, 0.004057, 0.005027, 0.005388, 0.006191, 0.006417, 0.006924, 0.007293, 0.007340, 0.007175, 0.005544, 0.017290, 0.003166};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv3_05_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.003647, 0.000666, 0.002345, 0.010674, 0.020633, 0.029478, 0.037807, 0.044923, 0.052115, 0.059175, 0.065205, 0.074218, 0.077106, 0.082393, 0.089078, 0.081828, 0.073522, 0.075938, 0.067269};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004758, 0.001683, 0.001029, 0.000836, 0.000794, 0.000838, 0.000957, 0.001115, 0.001294, 0.001584, 0.001933, 0.001800, 0.002665, 0.003967, 0.005624, 0.005425, 0.010326, 0.013241, 0.023159};
  Double_t _ysys[] = {0.001013, 0.000245, 0.000494, 0.000555, 0.001042, 0.001487, 0.001891, 0.002255, 0.002627, 0.003002, 0.003343, 0.003748, 0.003877, 0.004135, 0.004890, 0.004358, 0.003707, 0.004124, 0.004060};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv3_510_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.003787, -0.000149, 0.005713, 0.014075, 0.023652, 0.033837, 0.044640, 0.053538, 0.061587, 0.069731, 0.074860, 0.085709, 0.088923, 0.096413, 0.091740, 0.085116, 0.106610, 0.068623, 0.022712};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004620, 0.001699, 0.001043, 0.000858, 0.000833, 0.000879, 0.000972, 0.001140, 0.001362, 0.001633, 0.001969, 0.001819, 0.002577, 0.003814, 0.005396, 0.005238, 0.010101, 0.013129, 0.021670};
  Double_t _ysys[] = {0.000336, 0.000398, 0.000445, 0.000581, 0.000977, 0.001378, 0.001838, 0.002218, 0.002526, 0.003519, 0.002998, 0.003431, 0.003621, 0.005050, 0.003793, 0.003749, 0.005098, 0.004616, 0.006951};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv3_3040_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.001683, 0.006079, 0.012886, 0.025565, 0.040869, 0.051072, 0.064341, 0.074004, 0.084896, 0.091381, 0.098828, 0.101845, 0.108973, 0.108196, 0.110746, 0.105394, 0.094440, 0.129151, 0.079003};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004967, 0.001948, 0.001272, 0.001085, 0.001066, 0.001154, 0.001290, 0.001509, 0.001743, 0.002074, 0.002519, 0.002180, 0.003337, 0.004463, 0.006051, 0.005854, 0.010565, 0.013474, 0.026326};
  Double_t _ysys[] = {0.001352, 0.000662, 0.000479, 0.000802, 0.001418, 0.001544, 0.001965, 0.002444, 0.002565, 0.002998, 0.003136, 0.003062, 0.003630, 0.003348, 0.004205, 0.003653, 0.004495, 0.006762, 0.009228};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv3_05_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.027886, -0.017589, -0.008723, -0.003043, 0.005581, 0.015099, 0.027081, 0.040961, 0.049351, 0.063924, 0.081395, 0.097965, 0.108813, 0.114558, 0.103125, 0.093272, 0.087033};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.006162, 0.002556, 0.001698, 0.001399, 0.001249, 0.001188, 0.001214, 0.001288, 0.001370, 0.001168, 0.001509, 0.002123, 0.003149, 0.003627, 0.009132, 0.018091, 0.040165};
  Double_t _ysys[] = {0.003723, 0.002771, 0.000452, 0.000446, 0.001228, 0.001403, 0.002048, 0.002379, 0.002743, 0.003654, 0.004349, 0.005028, 0.005496, 0.005970, 0.005193, 0.006945, 0.037907};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv3_510_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.026865, -0.015319, -0.006049, 0.000713, 0.012147, 0.022344, 0.034518, 0.048523, 0.061018, 0.076231, 0.096640, 0.114272, 0.127109, 0.134074, 0.122842, 0.092986, 0.076467};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.005845, 0.002422, 0.001705, 0.001400, 0.001268, 0.001206, 0.001249, 0.001277, 0.001410, 0.001197, 0.001530, 0.002242, 0.003363, 0.003781, 0.009326, 0.017323, 0.044707};
  Double_t _ysys[] = {0.003573, 0.003004, 0.001444, 0.002206, 0.000804, 0.001032, 0.001700, 0.002355, 0.002558, 0.003161, 0.003905, 0.004778, 0.005111, 0.005989, 0.009185, 0.005342, 0.039017};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Lv3_3040_QC2(Int_t color=1, Int_t marker=20) {
  Int_t _nPoints = 17;
  Double_t _x[] = {0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {-0.007458, -0.002458, 0.007903, 0.022297, 0.035848, 0.047811, 0.065085, 0.079516, 0.091341, 0.107211, 0.128847, 0.139651, 0.152615, 0.162387, 0.162078, 0.146706, 0.152210};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.004968, 0.002495, 0.001817, 0.001582, 0.001429, 0.001516, 0.001591, 0.001717, 0.001936, 0.001642, 0.002285, 0.003113, 0.004521, 0.004913, 0.010693, 0.018956, 0.045095};
  Double_t _ysys[] = {0.004997, 0.002936, 0.001453, 0.000976, 0.001325, 0.002050, 0.003174, 0.002808, 0.002827, 0.003269, 0.003971, 0.004223, 0.004751, 0.004884, 0.005438, 0.010324, 0.008929};
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// ***********************  Xi flow ************************************************************************

// v2 for Xi in 5-10%
TGraphAsymmErrors* v2Xi0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 9;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
	Double_t _y[] = {0.00478264,0.0415018,0.0640852,0.0894661,0.0863706,0.120973,0.110222,0.176181,0.172742};
	Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t _yerrH[] =  {0.0083636,0.00470284,0.00441949,0.00540378,0.00771696,0.0112068,0.0160397,0.0262699,0.0328755};
	Double_t _yerrL[] = {0.0083636,0.00470284,0.00441949,0.00540378,0.00771696, 0.0112068,0.0160397,0.0262699,0.0328755};
	Double_t sysErrH[] = {0.00241906,0.00874917,0.00977195,0.0133772,0.0121558,0.0153083,0.0154306,0.0269707,0.0268289};
	Double_t sysErrL[] = {0.00245092,0.00886439,0.00977195,0.0133772,0.0121558,0.0153083,0.0154306,0.0269707,0.0268289};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerrH[i] = 0;
			_yerrL[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
		  _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
		  _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
		}
	}
	
	TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
							 &_x[first], 
							 &_y[first], 
							 &_xerr[first], 
							 &_xerr[first], 
							 &_yerrL[first], 
							 &_yerrH[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	return graph;
}

// v2 for Xi in 10-20%
TGraphAsymmErrors* v2Xi1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme                                                                                                            
  Int_t _nPoints = 9;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0298374,0.0681879,0.102952,0.134715,0.161414,0.178835,0.210738,0.177071,0.160964};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.00522021,0.00310098,0.00292487,0.00356788,0.00486366,0.00762432,0.0114765,0.0186645,0.0224286};
  Double_t _yerrL[] = {0.00522021,0.00310098,0.00292487,0.00356788,0.00486366,0.00762432,0.0114765,0.0186645,0.0224286};
  Double_t sysErrH[] = {0.00447339,0.00880626,0.00549531,0.00735239,0.0089939,0.00965687,0.0118862,0.00935573,0.00947006};
  Double_t sysErrL[] = {0.00466358,0.00918066,0.00549531,0.00735239,0.0089939,0.00965687,0.0118862,0.00935573,0.00947006};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
						   &_x[first],
						   &_y[first],
						   &_xerr[first],
						   &_xerr[first],
						   &_yerrL[first],
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Xi in 20-30%                                                                                                        
TGraphAsymmErrors* v2Xi2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme             
  Int_t _nPoints = 9;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0459866,0.098373,0.144052,0.186874,0.215865,0.228814,0.235308,0.247873,0.204186};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.00544938,0.00285116,0.00328415,0.00396176,0.00555323,0.0084036,0.0125012,0.0192451,0.023887};
  Double_t _yerrL[] = {0.00544938,0.00285116,0.00328415,0.00396176,0.00555323,0.0084036,0.0125012,0.0192451,0.023887};
  Double_t sysErrH[] = {0.00377003,0.0082164,0.00604054,0.00768675,0.00902968,0.00927681,0.00914397,0.0102107,0.00819201};
  Double_t sysErrL[] = {0.00415033,0.00904521,0.0075187,0.00956775,0.0112393,0.0115469,0.0113816,0.0127093,0.0101966};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Xi in 30-40%                                                                                                        
TGraphAsymmErrors* v2Xi3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme                                                                                                                 
  Int_t _nPoints = 9;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0696847,0.127436,0.184481,0.224152,0.236251,0.256007,0.263712,0.262631,0.216026};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.00645888,0.00417658,0.00551857,0.00668778,0.00707181,0.0104058,0.0155499,0.0230136,0.0283509};
  Double_t _yerrL[] = {0.00645888,0.00417658,0.00551857,0.00668778,0.00707181,0.0104058,0.0155499,0.0230136,0.0283509};
  Double_t sysErrH[] = {0.00505421,0.00985009,0.00727407,0.0089936,0.00952512,0.0104066,0.00995713,0.00941505,0.00842769};
  Double_t sysErrL[] = {0.00768463,0.0149765,0.017209,0.0212771,0.0225345,0.02462,0.0235566,0.0222741,0.0199382};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Xi in 40-50%                                                                                                        
TGraphAsymmErrors* v2Xi4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme                                                                                                                 
  Int_t _nPoints = 9;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.108667,0.143941,0.210207,0.246373,0.261297,0.28388,0.256455,0.233381,0.232871};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.00880908,0.00592602,0.00575493,0.00739944,0.0102572,0.0151811,0.0221249,0.0321652,0.0365146};
  Double_t _yerrL[] = {0.00880908,0.00592602,0.00575493,0.00739944,0.0102572,0.0151811,0.0221249,0.0321652,0.0365146};
  Double_t sysErrH[] = {0.00613324,0.0075205,0.00791397,0.00895048,0.00998192,0.0109442,0.00965148,0.00990504,0.00967035};
  Double_t sysErrL[] = {0.0180145,0.0220892,0.0320963,0.0363,0.0404832,0.0443857,0.039143,0.0401713,0.0392195};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Xi in 50-60%                                                                                                        
TGraphAsymmErrors* v2Xi5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme                                                                                                                 
  Int_t _nPoints = 9;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.119091,0.164749,0.214061,0.24545,0.250945,0.212644,0.195534,0.218019,0.190839};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0140518,0.00994419,0.0102977,0.0130384,0.0180408,0.0258517,0.037823,0.0520584,0.0613603};
  Double_t _yerrL[] = {0.0140518,0.00994419,0.0102977,0.0130384,0.0180408,0.0258517,0.037823,0.0520584,0.0613603};
  Double_t sysErrH[] = {0.00481909,0.0069128,0.00797333,0.00918592,0.00972359,0.00809423,0.00807699,0.00934366,0.0119026};
  Double_t sysErrL[] = {0.0179392,0.025733,0.0332072,0.0382574,0.0404967,0.0337108,0.0336389,0.0389143,0.0495717};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Omega in 5-10%
TGraphAsymmErrors* v2Omega0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme 
  Int_t _nPoints = 8;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0235771,0.0673959,0.0958733,0.111819,0.10579,0.0908653,0.159274,0.123938};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0320571,0.0227004,0.0261064,0.0291098,0.0344799,0.049148,0.0704155,0.0831405};
  Double_t _yerrL[] = {0.0320571,0.0227004,0.0261064,0.0291098,0.0344799,0.049148,0.0704155,0.0831405};
  Double_t sysErrH[] = {0.00865259,0.0128786,0.01763,0.0160203,0.0201751,0.0203363,0.0355451,0.0353583};
  Double_t sysErrL[] = {0.00876908,0.0128786,0.01763,0.0160203,0.0201751,0.0203363,0.0355451,0.0353583};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Omega in 10-20%
TGraphAsymmErrors* v2Omega1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme                                                                                                                 
  Int_t _nPoints = 8;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0100348,0.107418,0.121454,0.137174,0.171493,0.221296,0.179937,0.217591};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0188218,0.0151658,0.014954,0.0180363,0.0217451,0.0300963,0.0446354,0.0511276};
  Double_t _yerrL[] = {0.0188218,0.0151658,0.014954,0.0180363,0.0217451,0.0300963,0.0446354,0.0511276};
  Double_t sysErrH[] = {0.010783,0.010599,0.0141808,0.0173468,0.0186255,0.0229252,0.0180447,0.0182652};
  Double_t sysErrL[] = {0.0110909,0.0109134,0.0146015,0.0178614,0.0191781,0.0236054,0.01858,0.0188071};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Omega in 20-30%
TGraphAsymmErrors* v2Omega2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme             
  Int_t _nPoints = 8;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0318731,0.124177,0.182995,0.199765,0.203559,0.257703,0.219963,0.300351};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0196181,0.0147677,0.0153109,0.0191081,0.0234568,0.0314174,0.04477,0.0510465};
  Double_t _yerrL[] = {0.0196181,0.0147677,0.0153109,0.0191081,0.0234568,0.0314174,0.04477,0.0510465};
  Double_t sysErrH[] = {0.00994006,0.0116037,0.014766,0.0173458,0.0178205,0.0175653,0.0196145,0.0157366};
  Double_t sysErrL[] = {0.0106354,0.012622,0.0160618,0.0188679,0.0193843,0.0191067,0.0213357,0.0171176};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}


// v2 for Omega in 30-40%
TGraphAsymmErrors* v2Omega3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme             
  Int_t _nPoints = 8;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0711479,0.126749,0.223198,0.244707,0.21216,0.265079,0.25481,0.248865};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0229248,0.0179945,0.018397,0.022049,0.0283539,0.0350356,0.0575597,0.0646626};
  Double_t _yerrL[] = {0.0229248,0.0179945,0.018397,0.022049,0.0283539,0.0350356,0.0575597,0.0646626};
  Double_t sysErrH[] = {0.0115736,0.0140511,0.0173726,0.0183994,0.0201022,0.0192339,0.0181868,0.0162795};
  Double_t sysErrL[] = {0.0161623,0.0211234,0.0261168,0.0276603,0.0302202,0.0289149,0.0273407,0.0244735};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Omega in 40-50% 
TGraphAsymmErrors* v2Omega4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme             
  Int_t _nPoints = 8;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.146527,0.199316,0.212732,0.337707,0.328127,0.238654,0.327706,0.173636};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0305757,0.0249314,0.0278559,0.0320731,0.0397579,0.0543356,0.0688925,0.0809091};
  Double_t _yerrL[] = {0.0305757,0.0249314,0.0278559,0.0320731,0.0397579,0.0543356,0.0688925,0.0809091};
  Double_t sysErrH[] = {0.0115225,0.0162654,0.0183957,0.0205156,0.0224933,0.0198364,0.0203575,0.0198752};
  Double_t sysErrL[] = {0.0237517,0.0351013,0.0396986,0.0442734,0.0485414,0.0428078,0.0439324,0.0428915};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// v2 for Omega in 50-60%
TGraphAsymmErrors* v2Omega5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme             
  Int_t _nPoints = 7;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.5};
  Double_t _y[] = {0.0778981,0.236571,0.250775,0.240466,0.254392,0.262685,0.159428};
  Double_t _xerr[] = {0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t _yerrH[] = {0.0489679,0.044233,0.0456819,0.0582985,0.0735127,0.0952028,0.134597};
  Double_t _yerrL[] = {0.0489679,0.044233,0.0456819,0.0582985,0.0735127,0.0952028,0.134597};
  Double_t sysErrH[] = {0.013064,0.016725,0.0192686,0.0203964,0.0169786,0.0169424,0.0195994};
  Double_t sysErrL[] = {0.0280191,0.0363163,0.0418393,0.0442882,0.0368669,0.0367884,0.0425577};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = 0;
      _yerrL[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + sysErrH[i]*sysErrH[i]);
      _yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + sysErrL[i]*sysErrL[i]);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1,
                                                   &_x[first],
                                                   &_y[first],
                                                   &_xerr[first],
                                                   &_xerr[first],
                                                   &_yerrL[first],
                                                   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}


// ***********************  phi meson flow ************************************************************************

TGraphErrors* v2QC2Phi1060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 6;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5};
	Double_t _y[] = {0.0435127, 0.0967571, 0.153809, 0.185704, 0.195696, 0.154209 };
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] =  {0.00399773, 0.00367556, 0.00406245, 0.00486993, 0.00531555, 0.00979028, 0.0168998 };
	Double_t systPhi[] = {0.0107045, 0.0111588, 0.0171032, 0.0182152, 0.0206451, 0.0265506};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	return graph;
}

TGraphAsymmErrors* v2SPPhi1060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _y[] = {0.0497706,0.125828,0.191562,0.233368,0.248731,0.197478};
	Double_t _yerr[] = {0.00325214,0.00300283,0.00355632,0.00429377,0.00470054,0.00866158};
	Double_t _yerr1[] = {0.00325214,0.00300283,0.00355632,0.00429377,0.00470054,0.00866158};
	Double_t systPhi[] = {0.0108287,0.0104627,0.0161681,0.0171251,0.0196162,0.0282542};
	Double_t syst1Phi[] ={0.0114964,0.0114986,0.0195793,0.0275249,0.0329337,0.0465338};

	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

TGraphAsymmErrors* v2QC2Phi1040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] = {0.0314334,0.0834544,0.136862,0.17338,0.192154,0.152361,0.146152};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] =  {0.00399773, 0.00367556, 0.00406245, 0.00486993, 0.00531555, 0.00979028, 0.0168998 };
	Double_t _yerr1[] =  {0.00399773, 0.00367556, 0.00406245, 0.00486993, 0.00531555, 0.00979028, 0.0168998 };
	Double_t systPhi[] = {0.00779561, 0.00866222, 0.0160667, 0.0175538, 0.0183183, 0.0261145, 0.0520927 };
	Double_t syst1Phi[] = {0.00841917,0.00946453, 0.0184, 0.0249112, 0.0286591, 0.0399144, 0.0649789 };
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

TGraphAsymmErrors* v2SPPhi1040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _y[] = {0.0269576,0.0794705,0.13547,0.173661,0.193061,0.15142,0.146814};
	Double_t _yerr[] = {0.00386249,0.00363086,0.00413961,0.00500739,0.00549182,0.0101523,0.0175617};
	Double_t _yerr1[] = {0.00386249,0.00363086,0.00413961,0.00500739,0.00549182,0.0101523,0.0175617};
	Double_t systPhi[] = {0.00819262,0.00763137,0.015061,0.0161453,0.0166553,0.0276003,0.0579412};
	Double_t syst1Phi[] ={0.00878341,0.00855589,0.0175993,0.0239641,0.0276759,0.040985,0.0697352};

	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}


// v2{2} for phi meson in 10-20%
TGraphAsymmErrors* v2QC2Phi1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] = {0.0150641,0.0592553,0.103549,0.136659,0.164012,0.131764,0.0942402};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] =  {0.00660891,0.00627844,0.00642448,0.00797092,0.00880354,0.0164602,0.0284161};
	Double_t _yerr1[] =  {0.00660891,0.00627844,0.00642448,0.00797092,0.00880354,0.0164602,0.0284161};
	Double_t systPhi[] = {0.00397099, 0.00621952, 0.0107314, 0.0137233, 0.0183169, 0.0198932, 0.0471365};
	Double_t syst1Phi[] = {0.00468544, 0.00691432, 0.0128854, 0.0197752, 0.0253878, 0.0314116, 0.0562815};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{SP} for phi meson in 10-20%
TGraphAsymmErrors* v2SPPhi1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};	
	Double_t _y[] = {0.0117566,0.0562341,0.102551,0.136863,0.165251,0.131885,0.0927932 };
	Double_t _yerr[] = {0.00623673,0.00608246,0.00651475,0.00818594,0.00910491,0.0171075,0.0296194 };
	Double_t _yerr1[] = {0.00623673,0.00608246,0.00651475,0.00818594,0.00910491,0.0171075,0.0296194 };
	Double_t systPhi[] = {0.00420059, 0.00595573, 0.0106467, 0.013754, 0.0185174, 0.0204594, 0.0511227};
	Double_t syst1Phi[] = {0.00488155, 0.00667803, 0.012815, 0.0197966, 0.0255328, 0.0317732, 0.0596598};

	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{EP} for phi meson in 10-20%
TGraphErrors* v2EPPhi1020(Int_t color=3, Int_t marker=22, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] =  {0.0150641,0.0592553,0.103549,0.136659,0.164012,0.131764,0.0942402};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] =  {0.00660891,0.00627844,0.00642448,0.00797092,0.00880354,0.0164602,0.0284161};
	Double_t systPhi[] = {0.00495865, 0.00618128, 0.0109778, 0.012896, 0.0177644, 0.0173466, 0.0675537};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	return graph;
}

// v2{2} for phi meson in 20-30%
TGraphAsymmErrors* v2QC2Phi2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] =  {0.0403995,0.0863463,0.15077,0.188012,0.199488,0.182664,0.154738 };
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] = {0.0071221,0.0067618,0.00767483,0.00899532,0.00975384,0.0169804,0.0328156 };
	Double_t _yerr1[] = {0.0071221,0.0067618,0.00767483,0.00899532,0.00975384,0.0169804,0.0328156 };
	Double_t systPhi[] = {0.00829679, 0.00985329, 0.0224768, 0.0177579, 0.0135778, 0.0290853, 0.049472};
	Double_t syst1Phi[] = {0.00883775, 0.0105244, 0.0241132, 0.0248834, 0.0254462, 0.0416129, 0.0621651 };
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{SP} for phi meson in 20-30%
TGraphAsymmErrors* v2SPPhi2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _y[] = {0.0366716,0.0834632,0.149996,0.189005,0.200436,0.181567,0.157942 };
	Double_t _yerr[] = {0.00689717,0.0066824,0.00780865,0.00922471,0.0100461,0.0175647,0.0337462 };
	Double_t _yerr1[] = {0.00689717,0.0066824,0.00780865,0.00922471,0.0100461,0.0175647,0.0337462 };
	Double_t systPhi[] = {0.00880338, 0.0103119, 0.0232535, 0.017869, 0.0137857, 0.0331623, 0.0568229};
	Double_t syst1Phi[] = {0.00931496, 0.010955, 0.0248389, 0.0249628, 0.0255578, 0.044558, 0.0681606};

	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{EP} for phi meson in 20-30%
TGraphErrors* v2EPPhi2030(Int_t color=3, Int_t marker=22, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] = {0.0355455,0.0788878,0.14843,0.18069,0.190489,0.164176,0.131192 };
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] =  {0.00772047,0.00747238,0.00858801,0.0102612,0.0110171,0.0189772,0.0373178 };
	Double_t systPhi[] = {0.00795919, 0.0114993, 0.0211859, 0.0176374, 0.017292, 0.0191694, 0.0493063};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	return graph;
}


// ====================================================
// v2{2} for phi meson in 30-40%
TGraphAsymmErrors* v2QC2Phi3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] = {0.0413533,0.104016,0.166315,0.20161,0.216033,0.14348,0.188358 };
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] =  {0.00707788,0.00611069,0.00718185,0.00843123,0.0091371,0.0174769,0.0273762};
	Double_t _yerr1[] =  {0.00707788,0.00611069,0.00718185,0.00843123,0.0091371,0.0174769,0.0273762};
	Double_t systPhi[] = {0.0116873, 0.0100034, 0.0171211, 0.0216601, 0.0224799, 0.029981, 0.0585166};
	Double_t syst1Phi[] = {0.0122882, 0.0110147, 0.0202885, 0.0306819, 0.0350024, 0.0477006, 0.0750098};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{SP} for phi meson in 30-40%
TGraphAsymmErrors* v2SPPhi3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _y[] = {0.0361503,0.0998508,0.164688,0.201586,0.21652,0.141248,0.189121 };
	Double_t _yerr[] = {0.00701879,0.00615103,0.00737338,0.0087015,0.00945538,0.0181249,0.0285779 };
	Double_t _yerr1[] = {0.00701879,0.00615103,0.00737338,0.0087015,0.00945538,0.0181249,0.0285779 };
	Double_t systPhi[] = {0.0126161, 0.00707381, 0.013411, 0.0173136, 0.0171892, 0.0296933, 0.0650907};
	Double_t syst1Phi[] ={0.0131747, 0.0084436, 0.0172727, 0.0277845, 0.0318636, 0.0475204, 0.0802438 };
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}	
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{EP} for phi meson in 30-40%
TGraphErrors* v2EPPhi3040(Int_t color=3, Int_t marker=22, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 7;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5, 5.5};
	Double_t _y[] =  {0.0341134,0.0980865,0.156852,0.194519,0.20455,0.14186,0.166363};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] = {0.00790964,0.00696526,0.00832211,0.00974027,0.0105543,0.0203945,0.0317642};
	Double_t systPhi[] = {0.00629562, 0.00599028, 0.0120576, 0.0139109, 0.0142922, 0.0230755, 0.0427984};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	return graph;
}


// v2{2} for phi meson in 40-50%
TGraphAsymmErrors* v2QC2Phi4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 6;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5};
	Double_t _y[] = {0.069958,0.119243,0.202656,0.212522,0.209619,0.163328 };
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] = {0.00749632,0.006472,0.00826359,0.0098238,0.0104237,0.0197831};	
	Double_t _yerr1[] = {0.00749632,0.006472,0.00826359,0.0098238,0.0104237,0.0197831};	
	Double_t systPhi[] = {0.015809, 0.0146439, 0.0220097, 0.0231634, 0.0289709, 0.0232595};
	Double_t syst1Phi[] = {0.0165427, 0.0157946, 0.0260705, 0.0362577, 0.0450041, 0.0530007};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}
// v2{SP} for phi meson in 40-50%
TGraphAsymmErrors* v2SPPhi4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 6;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
	Double_t _y[] = {0.0619342,0.113185,0.20116,0.211573,0.208769,0.160066};
	Double_t _yerr[] ={0.00762882,0.00666748,0.00862832,0.0102957,0.0109352,0.020797};
	Double_t _yerr1[] ={0.00762882,0.00666748,0.00862832,0.0102957,0.0109352,0.020797};
	Double_t systPhi[] = {0.0161593, 0.0144677, 0.0219467, 0.0231281, 0.0294595, 0.0237431};
	Double_t syst1Phi[] ={0.0168778, 0.0156313, 0.0260173, 0.0362352, 0.0453202, 0.0532147};

	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}

// v2{2} for phi meson in 50-60%
TGraphAsymmErrors* v2QC2Phi5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 6;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5};
	Double_t _y[] = {0.0681429,0.127928,0.190599,0.227551,0.195058,0.151621 };
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t _yerr[] = {0.00929287,0.00834279,0.0108673,0.0131563,0.0150857,0.0253655};
	Double_t _yerr1[] = {0.00929287,0.00834279,0.0108673,0.0131563,0.0150857,0.0253655};
	Double_t systPhi[] = {0.018578, 0.0182298, 0.0160349, 0.0141674, 0.0219474, 0.0348885 }; 
	Double_t syst1Phi[] = {0.0196864, 0.0198723, 0.0246167, 0.0398876, 0.0509999, 0.072594};
	
	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}
// v2{SP} for phi meson in 50-60%
TGraphAsymmErrors* v2SPPhi5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
	//commentme
	Int_t _nPoints = 6;
	if (last>_nPoints-1) last=_nPoints-1;
	if (last<0 && first<0) last=_nPoints-1;
	if (last<0) last=_nPoints-1+last;
	if (first<0) first=0;
	Double_t _x[] = {0.9, 1.5, 2.1, 2.7, 3.5, 4.5};
	Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};	
	Double_t _y[] = {0.0527229,0.117624,0.185003,0.225434,0.190129,0.142005};
	Double_t _yerr[] = {0.00983516,0.00892543,0.0117201,0.0142501,0.016352,0.0275762};
	Double_t _yerr1[] = {0.00983516,0.00892543,0.0117201,0.0142501,0.016352,0.0275762};
	Double_t systPhi[] = {0.0190608, 0.020395, 0.0143803, 0.0135605, 0.0238556, 0.0410104};
	Double_t syst1Phi[] ={0.0201427, 0.0218755, 0.0235725, 0.0396761, 0.0518497, 0.0757267};

	if(!kStat){
		for(Int_t i=0;i<_nPoints;i++){
			_yerr[i] = 0;
			_yerr1[i] = 0;
			_xerr[i] = 0.05;
		}
	}
	
	if(kSyst){
		for(Int_t i=0;i<_nPoints;i++){
			Float_t systerr = systPhi[i];
			Float_t systerr1 = syst1Phi[i];
			_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
			_yerr1[i] = TMath::Sqrt(_yerr1[i]*_yerr1[i] + systerr1*systerr1);
		}
	}
	
	TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(marker);
	
	TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr1[first],&_yerr[first]);
	graph1->SetLineColor(color);
	graph1->SetMarkerColor(color);
	graph1->SetMarkerStyle(marker);
	return graph1;
}
//************************  phi meson flow  END *******************************************************************************************


TGraphAsymmErrors* v2Pion0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00824246, 0.0104024, 0.0129201, 0.0151121, 0.0175799, 0.0201892, 0.0226423, 0.0267943, 0.0305399, 0.0336959, 0.0375974, 0.0402311, 0.0433197, 0.0458848, 0.0478647, 0.0497621, 0.0523315, 0.0544167, 0.0563532, 0.0548481, 0.0567523, 0.0601806, 0.0608681, 0.0608455, 0.061474, 0.0584937, 0.0615577, 0.0570016, 0.057145, 0.0422802, 0.0320468, 0.035704, 0.0288706, 0.0226647, 0.0400384, 0.0106084, 0.0443823, 0.0543296, -0.0050952, 0.0338315, 0.028349};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000133001, 0.000130299, 0.000131827, 0.000135471, 0.000140745, 0.000147232, 0.000166315, 0.000177255, 0.000194932, 0.000216528, 0.000242008, 0.000271456, 0.000304427, 0.000342121, 0.00038542, 0.000435245, 0.000493359, 0.000558285, 0.000632774, 0.000717061, 0.00081164, 0.000633477, 0.000853291, 0.00113628, 0.001482, 0.00191709, 0.00243193, 0.00305123, 0.0037621, 0.00352847, 0.0049445, 0.00658436, 0.00840823, 0.010414, 0.012692, 0.0113534, 0.0153486, 0.0197452, 0.0192072, 0.0251356, 0.0357553};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.000133001, 0.000130299, 0.000131827, 0.000135471, 0.000140745, 0.000147232, 0.000166315, 0.000177255, 0.000194932, 0.000216528, 0.000242008, 0.000271456, 0.000304427, 0.000342121, 0.00038542, 0.000435245, 0.000493359, 0.000558285, 0.000632774, 0.000717061, 0.00081164, 0.000633477, 0.000853291, 0.00113628, 0.001482, 0.00191709, 0.00243193, 0.00305123, 0.0037621, 0.00352847, 0.0049445, 0.00658436, 0.00840823, 0.010414, 0.012692, 0.0113534, 0.0153486, 0.0197452, 0.0192072, 0.0251356, 0.0357553};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.23867300000000003e-03;
	Float_t pol1 = 6.71160999999999980e-04;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(0, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphAsymmErrors* v2Kaon0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00134701, 0.00317163, 0.0031598, 0.00390976, 0.00703514, 0.0114981, 0.013297, 0.0192762, 0.0226932, 0.0285625, 0.0326775, 0.0361104, 0.0381141, 0.0414397, 0.0461846, 0.0489079, 0.0510579, 0.053841, 0.0571077, 0.056352, 0.0609899, 0.0638829, 0.0650988, 0.0588761, 0.0665654, 0.0715115, 0.0578553, 0.0541069, 0.0748944, 0.0532155, 0.0590017, 0.0630081, 0.0316938, 0.0783863, 0.263692, 0.128943, 0.286008, 0.13549, -0.199536, 0.342328};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.000987431, 0.000805962, 0.000767625, 0.000932952, 0.00134873, 0.000760668, 0.00064691, 0.000599764, 0.00058296, 0.000590131, 0.000608036, 0.000638975, 0.000678305, 0.000727843, 0.000792013, 0.000861015, 0.000949487, 0.00105108, 0.00116845, 0.00130555, 0.00100112, 0.00134176, 0.00180046, 0.00244084, 0.00325403, 0.00431594, 0.00572013, 0.00755246, 0.00794723, 0.0139484, 0.0246449, 0.0407051, 0.0661392, 0.0996223, 0.128813, 0.245527, 0.355049, 0.312221, 0.442529, 0.973744};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.000987431, 0.000805962, 0.000767625, 0.000932952, 0.00134873, 0.000760668, 0.00064691, 0.000599764, 0.00058296, 0.000590131, 0.000608036, 0.000638975, 0.000678305, 0.000727843, 0.000792013, 0.000861015, 0.000949487, 0.00105108, 0.00116845, 0.00130555, 0.00100112, 0.00134176, 0.00180046, 0.00244084, 0.00325403, 0.00431594, 0.00572013, 0.00755246, 0.00794723, 0.0139484, 0.0246449, 0.0407051, 0.0661392, 0.0996223, 0.128813, 0.245527, 0.355049, 0.312221, 0.442529, 0.973744};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
  	Float_t pol0 = 1.23867300000000003e-03;
	Float_t pol1 = 6.71160999999999980e-04;
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(0, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00268584, 0.000878408, 0.00217129, 0.00143351, 0.00322161, 0.00453038, 0.00681432, 0.00762821, 0.0112393, 0.013895, 0.0172801, 0.021795, 0.026308, 0.0281495, 0.0325376, 0.0381697, 0.0407043, 0.0433402, 0.0510857, 0.052863, 0.0633415, 0.0675941, 0.0756583, 0.0766736, 0.0811324, 0.0861403, 0.0740058, 0.0820468, 0.0858504, 0.0493029, 0.0626479, 0.0850436, 0.138466, 0.0540207, -0.151327, 0.0230493, -0.0698608, -0.418291, -0.24384};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00320162, 0.00255264, 0.00209768, 0.00182089, 0.00111836, 0.000987526, 0.00135817, 0.00123459, 0.00116575, 0.0011279, 0.00110913, 0.00110802, 0.00112118, 0.00115314, 0.00119301, 0.00125462, 0.00133187, 0.00142342, 0.00153028, 0.00110658, 0.00137104, 0.00172252, 0.00218898, 0.00193816, 0.00251692, 0.00327637, 0.0042846, 0.00444436, 0.00754896, 0.0122858, 0.0193753, 0.0297557, 0.0469546, 0.0523895, 0.0995941, 0.179351, 0.229041, 0.222489, 0.288131};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00320162, 0.00255264, 0.00209768, 0.00182089, 0.00111836, 0.000987526, 0.00135817, 0.00123459, 0.00116575, 0.0011279, 0.00110913, 0.00110802, 0.00112118, 0.00115314, 0.00119301, 0.00125462, 0.00133187, 0.00142342, 0.00153028, 0.00110658, 0.00137104, 0.00172252, 0.00218898, 0.00193816, 0.00251692, 0.00327637, 0.0042846, 0.00444436, 0.00754896, 0.0122858, 0.0193753, 0.0297557, 0.0469546, 0.0523895, 0.0995941, 0.179351, 0.229041, 0.222489, 0.288131};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.69829837126647719e-03;
	Float_t pol1 = 1.21632530760273721e-03;
 	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(0, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0137263, 0.0172372, 0.0214492, 0.0254373, 0.029454, 0.0333164, 0.0379873, 0.0447695, 0.0512336, 0.0575635, 0.0629713, 0.0682518, 0.0737371, 0.0779552, 0.0820896, 0.0861183, 0.0897984, 0.0936959, 0.0968309, 0.0985123, 0.102473, 0.105277, 0.107769, 0.110146, 0.109389, 0.10991, 0.105635, 0.106773, 0.101545, 0.0922347, 0.0793453, 0.077956, 0.0672527, 0.0487211, 0.0664732, 0.0510957, 0.0581614, 0.0307352, 0.0329276, 0.0530418, 0.0541353};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000103947, 0.000102006, 0.000103226, 0.000106107, 0.000110274, 0.000115385, 0.000129675, 0.000138283, 0.000152102, 0.00016898, 0.000188771, 0.000211627, 0.000237187, 0.000266274, 0.000299556, 0.000337933, 0.0003821, 0.000431893, 0.000488473, 0.000552876, 0.000625098, 0.000486471, 0.000652335, 0.000863765, 0.00112109, 0.001441, 0.00182259, 0.00226753, 0.00279327, 0.0025993, 0.00362121, 0.00481019, 0.00607083, 0.00754175, 0.00920136, 0.00840368, 0.0113048, 0.0143551, 0.0142097, 0.0190686, 0.0262873};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.000103947, 0.000102006, 0.000103226, 0.000106107, 0.000110274, 0.000115385, 0.000129675, 0.000138283, 0.000152102, 0.00016898, 0.000188771, 0.000211627, 0.000237187, 0.000266274, 0.000299556, 0.000337933, 0.0003821, 0.000431893, 0.000488473, 0.000552876, 0.000625098, 0.000486471, 0.000652335, 0.000863765, 0.00112109, 0.001441, 0.00182259, 0.00226753, 0.00279327, 0.0025993, 0.00362121, 0.00481019, 0.00607083, 0.00754175, 0.00920136, 0.00840368, 0.0113048, 0.0143551, 0.0142097, 0.0190686, 0.0262873};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.34651999999999992e-03;
	Float_t pol1 = 7.29597999999999922e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(1, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00444172, 0.00461262, 0.00605512, 0.00803394, 0.0115242, 0.0168087, 0.0236315, 0.0320238, 0.0387647, 0.0459735, 0.0531368, 0.0608119, 0.0661568, 0.0723506, 0.0772225, 0.0821846, 0.0838383, 0.0900321, 0.093199, 0.0960424, 0.103874, 0.11007, 0.116071, 0.12213, 0.117554, 0.119993, 0.105736, 0.113478, 0.0964476, 0.106342, 0.110875, 0.139353, 0.158666, 0.215898, 0.104821, 0.0470767, -0.170325, -0.442802, -0.282548, -1.62346};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.000768421, 0.000626958, 0.000596685, 0.000720582, 0.00104589, 0.000591218, 0.000503924, 0.000467582, 0.000455346, 0.000460641, 0.000474365, 0.000497623, 0.000527587, 0.000565326, 0.000614276, 0.000666982, 0.000734439, 0.000811081, 0.00089994, 0.00100514, 0.000767714, 0.00102462, 0.00136677, 0.00184566, 0.00244838, 0.00323692, 0.00423034, 0.0056044, 0.00581909, 0.0101702, 0.0175893, 0.029834, 0.0477719, 0.0717025, 0.0865885, 0.161362, 0.234659, 0.209886, 0.275306, 0.161189};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.000768421, 0.000626958, 0.000596685, 0.000720582, 0.00104589, 0.000591218, 0.000503924, 0.000467582, 0.000455346, 0.000460641, 0.000474365, 0.000497623, 0.000527587, 0.000565326, 0.000614276, 0.000666982, 0.000734439, 0.000811081, 0.00089994, 0.00100514, 0.000767714, 0.00102462, 0.00136677, 0.00184566, 0.00244838, 0.00323692, 0.00423034, 0.0056044, 0.00581909, 0.0101702, 0.0175893, 0.029834, 0.0477719, 0.0717025, 0.0865885, 0.161362, 0.234659, 0.209886, 0.275306, 0.161189};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.34651999999999992e-03;
	Float_t pol1 = 7.29597999999999922e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(1, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.00246409, 0.00296217, 0.00441478, 0.00368168, 0.00503813, 0.0066585, 0.0104877, 0.0131502, 0.0179621, 0.0259345, 0.0298481, 0.0373793, 0.0422295, 0.0498634, 0.0559304, 0.0645674, 0.0691012, 0.077218, 0.0842499, 0.0947646, 0.10716, 0.122472, 0.132405, 0.133811, 0.140163, 0.148689, 0.147146, 0.144192, 0.135299, 0.129794, 0.137321, 0.121362, 0.0528114, 0.0732102, 0.0807746, -0.00812459, -0.0918688, -0.121837, 0.123071};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00245065, 0.00195951, 0.00161286, 0.00140387, 0.000861586, 0.000761485, 0.00104237, 0.00095234, 0.000900714, 0.000872152, 0.000860265, 0.000859144, 0.000872256, 0.000897719, 0.00093013, 0.000978436, 0.00103962, 0.00111159, 0.00119467, 0.000863984, 0.00106895, 0.00134145, 0.00170067, 0.00149682, 0.00193438, 0.00250201, 0.0032605, 0.00335482, 0.00558796, 0.00920953, 0.0145869, 0.0223702, 0.0344418, 0.0392677, 0.0736332, 0.147416, 0.17009, 0.225622, 0.217007};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00245065, 0.00195951, 0.00161286, 0.00140387, 0.000861586, 0.000761485, 0.00104237, 0.00095234, 0.000900714, 0.000872152, 0.000860265, 0.000859144, 0.000872256, 0.000897719, 0.00093013, 0.000978436, 0.00103962, 0.00111159, 0.00119467, 0.000863984, 0.00106895, 0.00134145, 0.00170067, 0.00149682, 0.00193438, 0.00250201, 0.0032605, 0.00335482, 0.00558796, 0.00920953, 0.0145869, 0.0223702, 0.0344418, 0.0392677, 0.0736332, 0.147416, 0.17009, 0.225622, 0.217007};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.84616637865581655e-03;
	Float_t pol1 = 1.32222872399611561e-03;
 	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(1, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0201224, 0.0252062, 0.0311323, 0.0369631, 0.0428155, 0.0483589, 0.0550696, 0.0651319, 0.0746561, 0.0834061, 0.0919292, 0.0993882, 0.106362, 0.11341, 0.119441, 0.125327, 0.130489, 0.135768, 0.140111, 0.1436, 0.147624, 0.15197, 0.157217, 0.15934, 0.159981, 0.157703, 0.15448, 0.151198, 0.142915, 0.132599, 0.113408, 0.10271, 0.0943079, 0.0879974, 0.0790374, 0.0795734, 0.0614621, 0.067106, 0.0609654, 0.0515049, 0.0223675};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 7.16967e-05, 7.04768e-05, 7.14158e-05, 7.34838e-05, 7.64088e-05, 8.00149e-05, 8.9464e-05, 9.55187e-05, 0.000105152, 0.000116857, 0.000130562, 0.000146276, 0.000163809, 0.000183749, 0.00020652, 0.000232603, 0.000262503, 0.000296109, 0.000334244, 0.000377637, 0.00042637, 0.00033071, 0.000441212, 0.000580034, 0.000749139, 0.000958424, 0.00120529, 0.00149142, 0.00182397, 0.0016818, 0.00233162, 0.00308243, 0.00393626, 0.00487286, 0.0058437, 0.00539979, 0.00725315, 0.00952612, 0.0093786, 0.0123107, 0.0176851};
  Double_t _yerr2[] = {0, 0, 0, 0, 7.16967e-05, 7.04768e-05, 7.14158e-05, 7.34838e-05, 7.64088e-05, 8.00149e-05, 8.9464e-05, 9.55187e-05, 0.000105152, 0.000116857, 0.000130562, 0.000146276, 0.000163809, 0.000183749, 0.00020652, 0.000232603, 0.000262503, 0.000296109, 0.000334244, 0.000377637, 0.00042637, 0.00033071, 0.000441212, 0.000580034, 0.000749139, 0.000958424, 0.00120529, 0.00149142, 0.00182397, 0.0016818, 0.00233162, 0.00308243, 0.00393626, 0.00487286, 0.0058437, 0.00539979, 0.00725315, 0.00952612, 0.0093786, 0.0123107, 0.0176851};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.36284000000000027e-04;
	Float_t pol1 = 8.43209000000000037e-04;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00640075, 0.00890777, 0.0114231, 0.0121152, 0.0193763, 0.0262718, 0.0373136, 0.04773, 0.0585364, 0.0685442, 0.0781197, 0.0884599, 0.0959829, 0.104197, 0.110977, 0.118551, 0.124407, 0.130268, 0.137113, 0.141952, 0.149588, 0.158874, 0.165006, 0.168625, 0.168325, 0.16697, 0.173543, 0.161645, 0.15523, 0.14532, 0.128752, 0.142204, 0.0843617, 0.269387, 0.0231167, 0.431226, 0.0167267, 0.315711, 0.235569, -0.000340595};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.000525284, 0.000429881, 0.000408989, 0.000491132, 0.000713549, 0.000405664, 0.000346767, 0.000322585, 0.000314851, 0.000318631, 0.00032842, 0.000344056, 0.000364314, 0.00039006, 0.000423204, 0.000459161, 0.000504151, 0.000555309, 0.000614542, 0.000683822, 0.000519796, 0.000688975, 0.00091857, 0.00122243, 0.00160751, 0.00211499, 0.00275193, 0.0035775, 0.00370145, 0.00631253, 0.0106874, 0.0182738, 0.0313279, 0.0476483, 0.0611653, 0.106198, 0.207285, 0.182173, 0.191574, 0.224581};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.000525284, 0.000429881, 0.000408989, 0.000491132, 0.000713549, 0.000405664, 0.000346767, 0.000322585, 0.000314851, 0.000318631, 0.00032842, 0.000344056, 0.000364314, 0.00039006, 0.000423204, 0.000459161, 0.000504151, 0.000555309, 0.000614542, 0.000683822, 0.000519796, 0.000688975, 0.00091857, 0.00122243, 0.00160751, 0.00211499, 0.00275193, 0.0035775, 0.00370145, 0.00631253, 0.0106874, 0.0182738, 0.0313279, 0.0476483, 0.0611653, 0.106198, 0.207285, 0.182173, 0.191574, 0.224581};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.36284000000000027e-04;
	Float_t pol1 = 8.43209000000000037e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00514464, 0.00344225, 0.00468868, 0.00546, 0.0082079, 0.012222, 0.0165934, 0.0233094, 0.0288221, 0.0374016, 0.0454486, 0.0563712, 0.064759, 0.075207, 0.0867709, 0.0963642, 0.106439, 0.116716, 0.12656, 0.138558, 0.157754, 0.173931, 0.185026, 0.199048, 0.206145, 0.208087, 0.213577, 0.213699, 0.21574, 0.192142, 0.152677, 0.172787, 0.122158, 0.11594, 0.170479, 0.143031, 0.0113144, 0.188084, 0.909081};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00163153, 0.00131053, 0.00108878, 0.000950215, 0.000582294, 0.000516432, 0.000703097, 0.000643592, 0.000611608, 0.000594501, 0.000588832, 0.000590229, 0.000600458, 0.000619565, 0.00064412, 0.00067877, 0.000720346, 0.000770869, 0.000828332, 0.000600635, 0.000741676, 0.000928701, 0.00116982, 0.00102807, 0.00131593, 0.00169927, 0.00218807, 0.00222766, 0.00363962, 0.00585315, 0.00929629, 0.013941, 0.020436, 0.0231194, 0.046922, 0.0895645, 0.111115, 0.118019, 0.100696};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00163153, 0.00131053, 0.00108878, 0.000950215, 0.000582294, 0.000516432, 0.000703097, 0.000643592, 0.000611608, 0.000594501, 0.000588832, 0.000590229, 0.000600458, 0.000619565, 0.00064412, 0.00067877, 0.000720346, 0.000770869, 0.000828332, 0.000600635, 0.000741676, 0.000928701, 0.00116982, 0.00102807, 0.00131593, 0.00169927, 0.00218807, 0.00222766, 0.00363962, 0.00585315, 0.00929629, 0.013941, 0.020436, 0.0231194, 0.046922, 0.0895645, 0.111115, 0.118019, 0.100696};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.13364693747489361e-03;
	Float_t pol1 = 1.52812297971200361e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0264906, 0.0328913, 0.0406349, 0.0482191, 0.0559519, 0.063266, 0.072363, 0.0854648, 0.0978324, 0.109596, 0.120324, 0.130441, 0.139475, 0.147956, 0.1559, 0.163043, 0.169642, 0.176537, 0.180972, 0.185598, 0.189913, 0.1954, 0.201295, 0.202612, 0.200368, 0.197602, 0.193443, 0.185868, 0.174288, 0.163585, 0.148547, 0.135952, 0.126668, 0.106048, 0.112211, 0.0935685, 0.100671, 0.080048, 0.0687841, 0.0681896, 0.0283596};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 8.33757e-05, 8.2413e-05, 8.37108e-05, 8.6308e-05, 8.99009e-05, 9.42972e-05, 0.000104942, 0.000112303, 0.000123869, 0.000137888, 0.000154129, 0.000172772, 0.00019344, 0.000216903, 0.000243535, 0.000274009, 0.000308524, 0.000347492, 0.000391614, 0.000441228, 0.000496805, 0.000384006, 0.000509506, 0.000665912, 0.000855374, 0.00107254, 0.0013387, 0.00164823, 0.00200117, 0.00183046, 0.00250297, 0.00330227, 0.00420942, 0.00522087, 0.00628707, 0.00569617, 0.00777439, 0.0102038, 0.00997685, 0.013336, 0.0187935};
  Double_t _yerr2[] = {0, 0, 0, 0, 8.33757e-05, 8.2413e-05, 8.37108e-05, 8.6308e-05, 8.99009e-05, 9.42972e-05, 0.000104942, 0.000112303, 0.000123869, 0.000137888, 0.000154129, 0.000172772, 0.00019344, 0.000216903, 0.000243535, 0.000274009, 0.000308524, 0.000347492, 0.000391614, 0.000441228, 0.000496805, 0.000384006, 0.000509506, 0.000665912, 0.000855374, 0.00107254, 0.0013387, 0.00164823, 0.00200117, 0.00183046, 0.00250297, 0.00330227, 0.00420942, 0.00522087, 0.00628707, 0.00569617, 0.00777439, 0.0102038, 0.00997685, 0.013336, 0.0187935};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.90101399999999992e-03;
	Float_t pol1 = 1.03004500000000005e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(3, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0100568, 0.0124986, 0.0168451, 0.0175805, 0.0283886, 0.0371263, 0.0509478, 0.0644455, 0.0780852, 0.0932833, 0.104516, 0.115791, 0.127519, 0.137186, 0.146611, 0.155235, 0.161983, 0.169526, 0.174815, 0.182918, 0.190393, 0.199666, 0.204089, 0.206915, 0.210964, 0.20384, 0.199675, 0.195183, 0.176372, 0.172127, 0.137577, 0.149046, 0.183246, 0.155882, 0.0653847, -0.00956166, 0.314558, 0.285004, -0.0182406, 0.30176};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00060424, 0.000495522, 0.000474726, 0.000575723, 0.000823512, 0.000470094, 0.000404114, 0.0003779, 0.000371012, 0.000376148, 0.000388761, 0.000407363, 0.000431318, 0.000462215, 0.000501103, 0.000543433, 0.00059585, 0.0006553, 0.000723921, 0.000804308, 0.000610234, 0.000801757, 0.00105659, 0.00139631, 0.00186122, 0.00243238, 0.00315255, 0.00411179, 0.00426523, 0.00727213, 0.0126056, 0.0218825, 0.03676, 0.0584615, 0.0738541, 0.148133, 0.190062, 0.255908, 0.221697, 0.310706};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.00060424, 0.000495522, 0.000474726, 0.000575723, 0.000823512, 0.000470094, 0.000404114, 0.0003779, 0.000371012, 0.000376148, 0.000388761, 0.000407363, 0.000431318, 0.000462215, 0.000501103, 0.000543433, 0.00059585, 0.0006553, 0.000723921, 0.000804308, 0.000610234, 0.000801757, 0.00105659, 0.00139631, 0.00186122, 0.00243238, 0.00315255, 0.00411179, 0.00426523, 0.00727213, 0.0126056, 0.0218825, 0.03676, 0.0584615, 0.0738541, 0.148133, 0.190062, 0.255908, 0.221697, 0.310706};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.90101399999999992e-03;
	Float_t pol1 = 1.03004500000000005e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(3, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.0055522, 0.00689601, 0.00484755, 0.00818012, 0.0119217, 0.0176936, 0.0247658, 0.0330758, 0.0432781, 0.0555955, 0.0660769, 0.0792523, 0.0923913, 0.105875, 0.118953, 0.130771, 0.144687, 0.157169, 0.167856, 0.184878, 0.207554, 0.228001, 0.243, 0.25184, 0.262881, 0.265347, 0.264864, 0.266241, 0.261999, 0.254697, 0.223834, 0.22018, 0.135383, 0.162749, 0.148328, 0.165155, 0.256436, -0.0440721, 0.430871};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00181005, 0.0014695, 0.00122215, 0.00107081, 0.000656891, 0.000584984, 0.00079404, 0.000732572, 0.000699401, 0.000684536, 0.000681312, 0.00068764, 0.00070337, 0.000728803, 0.000762285, 0.000804511, 0.0008565, 0.000919188, 0.000989172, 0.000718546, 0.000889519, 0.00111315, 0.00139867, 0.00121871, 0.00155384, 0.0019896, 0.00255073, 0.00256814, 0.00413703, 0.00660142, 0.0100799, 0.0156344, 0.0222041, 0.0265145, 0.0477444, 0.102528, 0.104764, 0.134924, 0.17457};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00181005, 0.0014695, 0.00122215, 0.00107081, 0.000656891, 0.000584984, 0.00079404, 0.000732572, 0.000699401, 0.000684536, 0.000681312, 0.00068764, 0.00070337, 0.000728803, 0.000762285, 0.000804511, 0.0008565, 0.000919188, 0.000989172, 0.000718546, 0.000889519, 0.00111315, 0.00139867, 0.00121871, 0.00155384, 0.0019896, 0.00255073, 0.00256814, 0.00413703, 0.00660142, 0.0100799, 0.0156344, 0.0222041, 0.0265145, 0.0477444, 0.102528, 0.104764, 0.134924, 0.17457};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.60641395006084133e-03;
	Float_t pol1 = 1.86671983155921230e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(3, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0297602, 0.0376167, 0.0463315, 0.0556559, 0.0642205, 0.0730444, 0.0836154, 0.0989402, 0.113114, 0.126486, 0.139145, 0.150305, 0.16061, 0.170721, 0.178498, 0.186109, 0.194151, 0.198377, 0.205201, 0.208348, 0.213488, 0.218346, 0.22189, 0.223529, 0.221045, 0.21347, 0.209264, 0.200217, 0.187252, 0.17979, 0.166723, 0.151286, 0.135703, 0.131906, 0.124924, 0.120519, 0.098863, 0.109076, 0.0705744, 0.0527536, 0.0317052};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000109358, 0.000107746, 0.000109743, 0.000113428, 0.000118448, 0.000124546, 0.000138355, 0.000148586, 0.000164345, 0.000183373, 0.000205473, 0.000230642, 0.000258532, 0.000289823, 0.000325585, 0.000365999, 0.000411749, 0.000463224, 0.000521013, 0.000585997, 0.000658582, 0.000506735, 0.000667446, 0.000864468, 0.0011026, 0.00137946, 0.00170617, 0.00207612, 0.00250699, 0.00227451, 0.00308872, 0.00405663, 0.0051683, 0.00635124, 0.00778926, 0.00704853, 0.00964995, 0.0125646, 0.0127765, 0.0168288, 0.0240879};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.000109358, 0.000107746, 0.000109743, 0.000113428, 0.000118448, 0.000124546, 0.000138355, 0.000148586, 0.000164345, 0.000183373, 0.000205473, 0.000230642, 0.000258532, 0.000289823, 0.000325585, 0.000365999, 0.000411749, 0.000463224, 0.000521013, 0.000585997, 0.000658582, 0.000506735, 0.000667446, 0.000864468, 0.0011026, 0.00137946, 0.00170617, 0.00207612, 0.00250699, 0.00227451, 0.00308872, 0.00405663, 0.0051683, 0.00635124, 0.00778926, 0.00704853, 0.00964995, 0.0125646, 0.0127765, 0.0168288, 0.0240879};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.36398699999999969e-03;
	Float_t pol1 = 1.28090200000000000e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(4, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0118116, 0.0163667, 0.0202925, 0.0260998, 0.0355875, 0.0478286, 0.0614344, 0.0776, 0.0945835, 0.110036, 0.123176, 0.137558, 0.147839, 0.158324, 0.168516, 0.177569, 0.186743, 0.193229, 0.197969, 0.205318, 0.211439, 0.22, 0.224539, 0.224465, 0.223164, 0.217814, 0.215515, 0.207766, 0.195001, 0.178179, 0.19419, 0.165809, 0.177231, 0.135574, 0.0795671, 0.0127847, 0.0560518, -0.21127, -0.236194, -0.389479};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.000774866, 0.000639173, 0.000617921, 0.000754178, 0.00106575, 0.00061229, 0.000529924, 0.000498954, 0.000492243, 0.000501011, 0.000520013, 0.000546239, 0.000579176, 0.000621522, 0.000674318, 0.00073272, 0.000802985, 0.000882957, 0.000973762, 0.00107984, 0.000815824, 0.00106284, 0.00138328, 0.00179768, 0.00234626, 0.00300534, 0.00382029, 0.00482682, 0.00474745, 0.00752479, 0.0117531, 0.0179539, 0.0272536, 0.0425466, 0.0499865, 0.0932052, 0.144085, 0.200123, 0.562735, 0.572957};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.000774866, 0.000639173, 0.000617921, 0.000754178, 0.00106575, 0.00061229, 0.000529924, 0.000498954, 0.000492243, 0.000501011, 0.000520013, 0.000546239, 0.000579176, 0.000621522, 0.000674318, 0.00073272, 0.000802985, 0.000882957, 0.000973762, 0.00107984, 0.000815824, 0.00106284, 0.00138328, 0.00179768, 0.00234626, 0.00300534, 0.00382029, 0.00482682, 0.00474745, 0.00752479, 0.0117531, 0.0179539, 0.0272536, 0.0425466, 0.0499865, 0.0932052, 0.144085, 0.200123, 0.562735, 0.572957};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.36398699999999969e-03;
	Float_t pol1 = 1.28090200000000000e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(4, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.0134608, 0.00872409, 0.0100999, 0.0141023, 0.0197051, 0.0258617, 0.0336928, 0.0463434, 0.0576032, 0.0717588, 0.0843552, 0.100229, 0.116334, 0.129007, 0.147308, 0.163198, 0.172034, 0.189326, 0.202834, 0.216141, 0.244193, 0.259947, 0.279023, 0.284821, 0.291713, 0.296845, 0.287988, 0.283412, 0.260332, 0.250927, 0.227238, 0.208844, 0.19487, 0.168319, 0.125243, 0.249898, 0.220953, 0.359189, 0.46984};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00223486, 0.00181593, 0.00152581, 0.00133558, 0.000823213, 0.000739234, 0.00100455, 0.000930685, 0.000897515, 0.000884604, 0.000889995, 0.00090428, 0.00093094, 0.000969855, 0.00101916, 0.00108241, 0.00115601, 0.00124551, 0.0013468, 0.000981188, 0.00121854, 0.00152492, 0.00192614, 0.00166776, 0.00212082, 0.00269735, 0.00344015, 0.00344034, 0.00548854, 0.00841199, 0.0127197, 0.0186227, 0.0259636, 0.0277095, 0.0473207, 0.0821222, 0.104711, 0.170825, 0.259822};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00223486, 0.00181593, 0.00152581, 0.00133558, 0.000823213, 0.000739234, 0.00100455, 0.000930685, 0.000897515, 0.000884604, 0.000889995, 0.00090428, 0.00093094, 0.000969855, 0.00101916, 0.00108241, 0.00115601, 0.00124551, 0.0013468, 0.000981188, 0.00121854, 0.00152492, 0.00192614, 0.00166776, 0.00212082, 0.00269735, 0.00344015, 0.00344034, 0.00548854, 0.00841199, 0.0127197, 0.0186227, 0.0259636, 0.0277095, 0.0473207, 0.0821222, 0.104711, 0.170825, 0.259822};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 3.24118098903184706e-03;
	Float_t pol1 = 2.32134148520698243e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(4, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0312325, 0.0393881, 0.0492528, 0.0594062, 0.0688574, 0.0779215, 0.0897294, 0.105849, 0.121535, 0.13572, 0.148727, 0.160756, 0.17149, 0.180516, 0.189661, 0.196275, 0.203114, 0.209617, 0.21372, 0.217062, 0.222655, 0.223982, 0.227325, 0.224833, 0.220611, 0.215738, 0.206053, 0.199035, 0.194861, 0.178505, 0.15725, 0.147485, 0.131276, 0.140082, 0.129746, 0.119138, 0.0962541, 0.110656, 0.103842, 0.107556, 0.0562334};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000160203, 0.000158458, 0.000161984, 0.000167948, 0.000175992, 0.000185633, 0.000206362, 0.000222676, 0.000247387, 0.000276984, 0.000311291, 0.000350389, 0.000393386, 0.000441718, 0.000496267, 0.000558361, 0.000627807, 0.000705629, 0.00079231, 0.000890586, 0.000997622, 0.0007646, 0.000999208, 0.00128014, 0.00161834, 0.00199908, 0.00245428, 0.00297957, 0.00354308, 0.00320161, 0.00431837, 0.00563286, 0.0070497, 0.00882622, 0.010657, 0.00970365, 0.0132692, 0.0175344, 0.0174432, 0.0233582, 0.0345269};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.000160203, 0.000158458, 0.000161984, 0.000167948, 0.000175992, 0.000185633, 0.000206362, 0.000222676, 0.000247387, 0.000276984, 0.000311291, 0.000350389, 0.000393386, 0.000441718, 0.000496267, 0.000558361, 0.000627807, 0.000705629, 0.00079231, 0.000890586, 0.000997622, 0.0007646, 0.000999208, 0.00128014, 0.00161834, 0.00199908, 0.00245428, 0.00297957, 0.00354308, 0.00320161, 0.00431837, 0.00563286, 0.0070497, 0.00882622, 0.010657, 0.00970365, 0.0132692, 0.0175344, 0.0174432, 0.0233582, 0.0345269};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 3.02312299999999996e-03;
	Float_t pol1 = 1.63804399999999996e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(5, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0164446, 0.0194303, 0.025841, 0.0330678, 0.0434366, 0.0534827, 0.0699271, 0.0880427, 0.105501, 0.120385, 0.132267, 0.146892, 0.157761, 0.169855, 0.180112, 0.188618, 0.197751, 0.199262, 0.20814, 0.212515, 0.220558, 0.221804, 0.22802, 0.229711, 0.218398, 0.217427, 0.212665, 0.192197, 0.187403, 0.150519, 0.199072, 0.150158, 0.0439847, 0.228831, 0.226831, 0.0184782, -0.116235, -0.368718, 0.9956, 0.0201774};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00111871, 0.000922607, 0.000885857, 0.00106002, 0.00155242, 0.000899587, 0.000786825, 0.000747841, 0.000742576, 0.000760796, 0.000793007, 0.000836375, 0.000888282, 0.000956738, 0.00104026, 0.00113218, 0.00124419, 0.00137324, 0.00151438, 0.0016807, 0.00126703, 0.00164921, 0.00215552, 0.00277288, 0.00361221, 0.00456292, 0.00576781, 0.00718472, 0.00701046, 0.0109007, 0.0168378, 0.0258657, 0.0386923, 0.0541827, 0.0655031, 0.145958, 0.230034, 0.441662, 0.564042, 0};
 Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.00111871, 0.000922607, 0.000885857, 0.00106002, 0.00155242, 0.000899587, 0.000786825, 0.000747841, 0.000742576, 0.000760796, 0.000793007, 0.000836375, 0.000888282, 0.000956738, 0.00104026, 0.00113218, 0.00124419, 0.00137324, 0.00151438, 0.0016807, 0.00126703, 0.00164921, 0.00215552, 0.00277288, 0.00361221, 0.00456292, 0.00576781, 0.00718472, 0.00701046, 0.0109007, 0.0168378, 0.0258657, 0.0386923, 0.0541827, 0.0655031, 0.145958, 0.230034, 0.441662, 0.564042, 0};

  if(!kStat){
     for(Int_t i=0;i<_nPoints;i++){
	 _yerr[i] = 0;
	 _yerr2[i] = 0;
	 _xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 3.02312299999999996e-03;
	Float_t pol1 = 1.63804399999999996e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(5, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00886425, 0.0141768, 0.0101568, 0.0161863, 0.0219319, 0.0334622, 0.0443905, 0.0592395, 0.0751561, 0.0864854, 0.103567, 0.120373, 0.132621, 0.149678, 0.166081, 0.179562, 0.198764, 0.205123, 0.220172, 0.237277, 0.25903, 0.27644, 0.284633, 0.292908, 0.297426, 0.28959, 0.292919, 0.290167, 0.269848, 0.261622, 0.229602, 0.212424, 0.295423, 0.22624, 0.194105, 0.0526927, -0.076895, 0.217367, 0.983948};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00307721, 0.00251481, 0.00211521, 0.0018645, 0.00115595, 0.00104369, 0.00142809, 0.0013364, 0.00130398, 0.00129796, 0.00131775, 0.00135219, 0.00140202, 0.00147393, 0.00156009, 0.00166608, 0.00179217, 0.00193757, 0.00209925, 0.00153881, 0.00192408, 0.00241141, 0.00305588, 0.00264725, 0.0033566, 0.00426608, 0.00547565, 0.00538587, 0.00848915, 0.0131123, 0.0189988, 0.0279199, 0.0386729, 0.0415593, 0.0663031, 0.108596, 0.116259, 0.193609, 0.221554};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00307721, 0.00251481, 0.00211521, 0.0018645, 0.00115595, 0.00104369, 0.00142809, 0.0013364, 0.00130398, 0.00129796, 0.00131775, 0.00135219, 0.00140202, 0.00147393, 0.00156009, 0.00166608, 0.00179217, 0.00193757, 0.00209925, 0.00153881, 0.00192408, 0.00241141, 0.00305588, 0.00264725, 0.0033566, 0.00426608, 0.00547565, 0.00538587, 0.00848915, 0.0131123, 0.0189988, 0.0279199, 0.0386729, 0.0415593, 0.0663031, 0.108596, 0.116259, 0.193609, 0.221554};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 4.14488622049151260e-03;
	Float_t pol1 = 2.96857730797800597e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(5, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0303956, 0.0398194, 0.0497382, 0.0592532, 0.0692265, 0.0789018, 0.0906535, 0.107528, 0.122981, 0.137459, 0.150972, 0.16168, 0.170871, 0.180953, 0.186814, 0.195909, 0.200081, 0.203635, 0.209983, 0.21007, 0.213946, 0.213234, 0.216813, 0.215134, 0.2081, 0.207871, 0.185115, 0.194539, 0.179635, 0.169884, 0.165801, 0.145236, 0.12241, 0.119177, 0.123804, 0.111177, 0.101511, 0.139219, 0.0612838, 0.152114, 0.14425};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000275012, 0.000271644, 0.000278885, 0.000290482, 0.000305581, 0.000323479, 0.000360388, 0.000391261, 0.000437084, 0.000491624, 0.000554574, 0.000626577, 0.000704307, 0.000792713, 0.000891942, 0.00100363, 0.00113002, 0.00126999, 0.00142248, 0.00159975, 0.00178952, 0.00136889, 0.00177929, 0.00227244, 0.00284049, 0.00351991, 0.00429153, 0.0051599, 0.00614605, 0.0054959, 0.00736442, 0.00956336, 0.0120611, 0.0151024, 0.0177286, 0.0166113, 0.0232109, 0.0302416, 0.0304321, 0.0412011, 0.0613213};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.000275012, 0.000271644, 0.000278885, 0.000290482, 0.000305581, 0.000323479, 0.000360388, 0.000391261, 0.000437084, 0.000491624, 0.000554574, 0.000626577, 0.000704307, 0.000792713, 0.000891942, 0.00100363, 0.00113002, 0.00126999, 0.00142248, 0.00159975, 0.00178952, 0.00136889, 0.00177929, 0.00227244, 0.00284049, 0.00351991, 0.00429153, 0.0051599, 0.00614605, 0.0054959, 0.00736442, 0.00956336, 0.0120611, 0.0151024, 0.0177286, 0.0166113, 0.0232109, 0.0302416, 0.0304321, 0.0412011, 0.0613213};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 4.01257199999999984e-03;
	Float_t pol1 = 2.17416299999999979e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(6, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0190181, 0.0200425, 0.0286688, 0.03481, 0.0451018, 0.0588837, 0.0782952, 0.0947995, 0.108637, 0.1221, 0.135844, 0.148488, 0.157787, 0.169447, 0.173264, 0.188723, 0.18833, 0.195569, 0.199314, 0.205223, 0.211733, 0.208414, 0.213651, 0.208302, 0.19844, 0.196712, 0.200806, 0.15693, 0.18652, 0.115381, 0.133751, 0.0982841, 0.107246, 0.0417931, 0.147759, -0.0735589, 0.571869, 0.488767, -0.315304, 1.71937};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00188333, 0.00156073, 0.00151031, 0.00181199, 0.00266245, 0.00155702, 0.00137462, 0.00132048, 0.00131871, 0.00136175, 0.00142741, 0.00150997, 0.00161189, 0.00173747, 0.00189366, 0.00206161, 0.00226274, 0.00248796, 0.00273639, 0.00303486, 0.00227613, 0.00293909, 0.00376378, 0.00482482, 0.00608811, 0.00752018, 0.00942054, 0.0115408, 0.0108768, 0.0160409, 0.0235301, 0.0337122, 0.0466674, 0.0674531, 0.0743068, 0.11877, 0.253677, 0.235641, 0.512443, 0.266431};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.00188333, 0.00156073, 0.00151031, 0.00181199, 0.00266245, 0.00155702, 0.00137462, 0.00132048, 0.00131871, 0.00136175, 0.00142741, 0.00150997, 0.00161189, 0.00173747, 0.00189366, 0.00206161, 0.00226274, 0.00248796, 0.00273639, 0.00303486, 0.00227613, 0.00293909, 0.00376378, 0.00482482, 0.00608811, 0.00752018, 0.00942054, 0.0115408, 0.0108768, 0.0160409, 0.0235301, 0.0337122, 0.0466674, 0.0674531, 0.0743068, 0.11877, 0.253677, 0.235641, 0.512443, 0.266431};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 4.01257199999999984e-03;
	Float_t pol1 = 2.17416299999999979e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(6, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00919056, 0.0178105, 0.0134286, 0.0295189, 0.0277064, 0.0434262, 0.0535665, 0.0715688, 0.0842733, 0.106768, 0.117998, 0.137741, 0.150551, 0.161897, 0.178599, 0.192242, 0.20269, 0.219659, 0.227496, 0.242855, 0.257706, 0.265325, 0.281052, 0.2813, 0.271633, 0.267676, 0.272005, 0.277673, 0.273667, 0.196607, 0.24668, 0.181115, 0.31068, 0.114003, 0.333805, 0.276964, -0.0473593, 0.207914, 0.12068};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00490617, 0.00402616, 0.00341963, 0.00302967, 0.00188597, 0.00172265, 0.00237703, 0.00226014, 0.00222069, 0.00223993, 0.00230435, 0.00239039, 0.00250534, 0.00265745, 0.00283418, 0.00304524, 0.00330038, 0.00358151, 0.00389282, 0.00287646, 0.00361525, 0.00453627, 0.00570268, 0.00498584, 0.00629115, 0.00797278, 0.010008, 0.00986531, 0.0151841, 0.0225399, 0.0331056, 0.0484569, 0.0652464, 0.0676624, 0.111366, 0.15898, 0.211355, 0.256669, 0.378091};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00490617, 0.00402616, 0.00341963, 0.00302967, 0.00188597, 0.00172265, 0.00237703, 0.00226014, 0.00222069, 0.00223993, 0.00230435, 0.00239039, 0.00250534, 0.00265745, 0.00283418, 0.00304524, 0.00330038, 0.00358151, 0.00389282, 0.00287646, 0.00361525, 0.00453627, 0.00570268, 0.00498584, 0.00629115, 0.00797278, 0.010008, 0.00986531, 0.0151841, 0.0225399, 0.0331056, 0.0484569, 0.0652464, 0.0676624, 0.111366, 0.15898, 0.211355, 0.256669, 0.378091};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 5.50148450680349593e-03;
	Float_t pol1 = 3.94017620709329313e-03;

 	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(6, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0275772, 0.0353857, 0.0451589, 0.0569504, 0.0645817, 0.0744759, 0.0846387, 0.103224, 0.116637, 0.12955, 0.141348, 0.150827, 0.163529, 0.166875, 0.175975, 0.177624, 0.186353, 0.18951, 0.191246, 0.197078, 0.195822, 0.194678, 0.19574, 0.182454, 0.186385, 0.192537, 0.171284, 0.176511, 0.15328, 0.170476, 0.165746, 0.142293, 0.153372, 0.152722, 0.132789, 0.103282, 0.0631696, 0.142603, 0.126161, 0.119172, -0.0503365};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000573254, 0.0005678, 0.000585761, 0.000612834, 0.00064773, 0.000688254, 0.000769859, 0.000841415, 0.000944491, 0.00106887, 0.00120962, 0.00137073, 0.00154395, 0.00174052, 0.00195924, 0.00220447, 0.00247752, 0.00278203, 0.00312082, 0.00349313, 0.00390327, 0.00297633, 0.00383881, 0.00485695, 0.00607907, 0.00740945, 0.00894688, 0.0107216, 0.0126575, 0.0112992, 0.0149453, 0.0195003, 0.0243844, 0.0308124, 0.036297, 0.0342476, 0.0454655, 0.0641991, 0.063416, 0.0867013, 0.127476};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.000573254, 0.0005678, 0.000585761, 0.000612834, 0.00064773, 0.000688254, 0.000769859, 0.000841415, 0.000944491, 0.00106887, 0.00120962, 0.00137073, 0.00154395, 0.00174052, 0.00195924, 0.00220447, 0.00247752, 0.00278203, 0.00312082, 0.00349313, 0.00390327, 0.00297633, 0.00383881, 0.00485695, 0.00607907, 0.00740945, 0.00894688, 0.0107216, 0.0126575, 0.0112992, 0.0149453, 0.0195003, 0.0243844, 0.0308124, 0.036297, 0.0342476, 0.0454655, 0.0641991, 0.063416, 0.0867013, 0.127476};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 5.49768000000000045e-03;
	Float_t pol1 = 2.97884999999999990e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(7, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0166047, 0.0248478, 0.030513, 0.0411593, 0.0380696, 0.0639456, 0.0758458, 0.0897174, 0.102605, 0.116338, 0.133685, 0.140557, 0.141469, 0.15527, 0.167424, 0.171481, 0.179479, 0.187214, 0.177038, 0.178377, 0.187908, 0.181605, 0.20429, 0.173224, 0.18915, 0.193766, 0.195795, 0.189594, 0.17632, 0.158214, 0.154361, 0.0684557, 0.263961, 0.0779351, 0.0523301, -0.0182867, 0.264014, 0.135541, 0.027261, -3.29017};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00388478, 0.00326808, 0.00326331, 0.00427722, 0.00561727, 0.0033143, 0.00296018, 0.00287265, 0.00289959, 0.00299629, 0.00316135, 0.00335928, 0.00358978, 0.0038917, 0.00425027, 0.00463553, 0.0050853, 0.0056094, 0.00617222, 0.00685316, 0.00509904, 0.00653626, 0.00835353, 0.010359, 0.0134302, 0.0166309, 0.0204693, 0.0247671, 0.0232361, 0.0352802, 0.0483286, 0.0662901, 0.0989231, 0.12961, 0.143373, 0.313321, 0.49827, 0.675383, 0, 0.348563};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.00388478, 0.00326808, 0.00326331, 0.00427722, 0.00561727, 0.0033143, 0.00296018, 0.00287265, 0.00289959, 0.00299629, 0.00316135, 0.00335928, 0.00358978, 0.0038917, 0.00425027, 0.00463553, 0.0050853, 0.0056094, 0.00617222, 0.00685316, 0.00509904, 0.00653626, 0.00835353, 0.010359, 0.0134302, 0.0166309, 0.0204693, 0.0247671, 0.0232361, 0.0352802, 0.0483286, 0.0662901, 0.0989231, 0.12961, 0.143373, 0.313321, 0.49827, 0.675383, 0, 0.348563};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 5.49768000000000045e-03;
	Float_t pol1 = 2.97884999999999990e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(7, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.0130981, 0.0179731, 0.00939326, 0.0248886, 0.029577, 0.0550703, 0.0636149, 0.0798523, 0.104272, 0.112925, 0.115704, 0.146685, 0.149302, 0.179475, 0.181385, 0.184513, 0.208168, 0.213384, 0.225983, 0.232891, 0.231734, 0.241272, 0.270741, 0.246586, 0.266346, 0.227266, 0.245789, 0.208304, 0.234844, 0.194027, -0.0336696, 0.267591, 0.2404, 0.241992, 0.125957, 0.288086, -0.0589748, 0.311092, 0.477676};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00943873, 0.00785252, 0.00670716, 0.00596963, 0.00376045, 0.00346777, 0.00486688, 0.00467992, 0.00468118, 0.00479058, 0.00497572, 0.00521532, 0.00552991, 0.00590979, 0.00633209, 0.00685181, 0.00747093, 0.00816394, 0.00889636, 0.00662104, 0.00834123, 0.0105143, 0.0133603, 0.0113902, 0.014318, 0.0179343, 0.0226294, 0.0223227, 0.0342436, 0.0511045, 0.0730614, 0.097839, 0.132614, 0.150997, 0.218463, 0.348699, 0.362474, 0.57912, 0.796176};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00943873, 0.00785252, 0.00670716, 0.00596963, 0.00376045, 0.00346777, 0.00486688, 0.00467992, 0.00468118, 0.00479058, 0.00497572, 0.00521532, 0.00552991, 0.00590979, 0.00633209, 0.00685181, 0.00747093, 0.00816394, 0.00889636, 0.00662104, 0.00834123, 0.0105143, 0.0133603, 0.0113902, 0.014318, 0.0179343, 0.0226294, 0.0223227, 0.0342436, 0.0511045, 0.0730614, 0.097839, 0.132614, 0.150997, 0.218463, 0.348699, 0.362474, 0.57912, 0.796176};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.53765601848897532e-03;
	Float_t pol1 = 5.39848705646177145e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(7, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0195733, 0.0309692, 0.0356384, 0.0462504, 0.0543824, 0.0647011, 0.0749217, 0.085954, 0.0999237, 0.117407, 0.12294, 0.137339, 0.140553, 0.145824, 0.156141, 0.15101, 0.159541, 0.163759, 0.171796, 0.179879, 0.167389, 0.166001, 0.176266, 0.179223, 0.18775, 0.185847, 0.206443, 0.148862, 0.175494, 0.189032, 0.233424, 0.0839163, 0.263219, 0.0629226, 0.320809, 0.134834, 0.190157, 0.0960499, 0.136835, -0.0485677, 0.389939};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.00157008, 0.00155061, 0.00160592, 0.00168889, 0.0017938, 0.00191324, 0.00215022, 0.00236305, 0.00267145, 0.00303387, 0.00344413, 0.00390738, 0.00441634, 0.00498694, 0.00562308, 0.00632594, 0.00710876, 0.00796577, 0.00892046, 0.00999595, 0.0111319, 0.00845832, 0.0108732, 0.0136163, 0.0167438, 0.0206408, 0.0247138, 0.0295918, 0.0351168, 0.0307851, 0.0417363, 0.0528758, 0.0669722, 0.0817049, 0.0994534, 0.0947033, 0.126801, 0.177027, 0.176152, 0.245739, 0.391682};
  Double_t _yerr2[] = {0, 0, 0, 0, 0.00157008, 0.00155061, 0.00160592, 0.00168889, 0.0017938, 0.00191324, 0.00215022, 0.00236305, 0.00267145, 0.00303387, 0.00344413, 0.00390738, 0.00441634, 0.00498694, 0.00562308, 0.00632594, 0.00710876, 0.00796577, 0.00892046, 0.00999595, 0.0111319, 0.00845832, 0.0108732, 0.0136163, 0.0167438, 0.0206408, 0.0247138, 0.0295918, 0.0351168, 0.0307851, 0.0417363, 0.0528758, 0.0669722, 0.0817049, 0.0994534, 0.0947033, 0.126801, 0.177027, 0.176152, 0.245739, 0.391682};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 8.29529000000000030e-03;
	Float_t pol1 = 5.60073999999999989e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(8, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00726463, 0.0144253, 0.0334631, 0.0384761, 0.0494532, 0.0693552, 0.0762569, 0.0719228, 0.0937563, 0.0955609, 0.105749, 0.117711, 0.129626, 0.133837, 0.128721, 0.132298, 0.161963, 0.138618, 0.139166, 0.18394, 0.169776, 0.181286, 0.173951, 0.14043, 0.115846, 0.214474, 0.218691, 0.103945, 0.37536, -0.0865605, 0.228924, -0.115581, 0.110384, 0.184704, 0.509455, -0.485498, 0.0837004, -0.813161, -0.0297345, 0.0323336};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.0106177, 0.00884964, 0.00836979, 0.00948119, 0.0156272, 0.00929314, 0.00842046, 0.00821347, 0.0083515, 0.00868388, 0.00921062, 0.00978082, 0.0105426, 0.0114431, 0.0124756, 0.0136669, 0.0149571, 0.0164684, 0.0181788, 0.0200314, 0.014819, 0.0190758, 0.0243153, 0.0308707, 0.0377526, 0.0457948, 0.0569713, 0.0684405, 0.0634441, 0.0935611, 0.129774, 0.182993, 0.241246, 0.29603, 0.342445, 0.523882, 0.918355, 1.40502, 0, 0};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.0106177, 0.00884964, 0.00836979, 0.00948119, 0.0156272, 0.00929314, 0.00842046, 0.00821347, 0.0083515, 0.00868388, 0.00921062, 0.00978082, 0.0105426, 0.0114431, 0.0124756, 0.0136669, 0.0149571, 0.0164684, 0.0181788, 0.0200314, 0.014819, 0.0190758, 0.0243153, 0.0308707, 0.0377526, 0.0457948, 0.0569713, 0.0684405, 0.0634441, 0.0935611, 0.129774, 0.182993, 0.241246, 0.29603, 0.342445, 0.523882, 0.918355, 1.40502, 0, 0};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 8.29529000000000030e-03;
	Float_t pol1 = 5.60073999999999989e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(8, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Antiproton7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.0243855, 4.59379e-05, 0.00505587, 0.0407476, 0.0360989, 0.0659522, 0.0654352, 0.0672759, 0.0947919, 0.0946055, 0.108647, 0.108389, 0.120068, 0.154908, 0.201656, 0.179421, 0.179193, 0.209699, 0.209577, 0.190403, 0.267222, 0.257467, 0.17816, 0.252145, 0.270207, 0.311401, 0.147204, 0.197523, 0.149844, 0.175334, 0.0566661, 0.322004, -0.123699, 0.1207, 0.101863, -0.448652, -1.05087, 2.292, 0.286693};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.0243279, 0.020212, 0.0174078, 0.0154737, 0.00988768, 0.00926115, 0.0131909, 0.0128732, 0.0130524, 0.0135259, 0.0142208, 0.0150952, 0.0160249, 0.0173254, 0.0187753, 0.0203282, 0.0222585, 0.0243728, 0.0268049, 0.0199116, 0.0251891, 0.0316377, 0.0397848, 0.0345062, 0.0433967, 0.0542165, 0.0687569, 0.0658095, 0.0953791, 0.140737, 0.210719, 0.265241, 0.339048, 0.375181, 0.545855, 1.04449, 0.904123, 0.837606, 2.36757};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.0243279, 0.020212, 0.0174078, 0.0154737, 0.00988768, 0.00926115, 0.0131909, 0.0128732, 0.0130524, 0.0135259, 0.0142208, 0.0150952, 0.0160249, 0.0173254, 0.0187753, 0.0203282, 0.0222585, 0.0243728, 0.0268049, 0.0199116, 0.0251891, 0.0316377, 0.0397848, 0.0345062, 0.0433967, 0.0542165, 0.0687569, 0.0658095, 0.0953791, 0.140737, 0.210719, 0.265241, 0.339048, 0.375181, 0.545855, 1.04449, 0.904123, 0.837606, 2.36757};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = 0;
	_yerr2[i] = 0;
	_xerr[i] = 0.05;
   }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.13733682403531055e-02;
	Float_t pol1 = 8.14563321585857515e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(8, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* ptscaling(TGraphErrors *g,Int_t nq){
  Int_t n = g->GetN();
  Double_t x[100],y[100],ex[100],ey[100];
  Double_t *xo,*yo,*eyo,*exo;
  xo = g->GetX();
  yo = g->GetY();
  eyo = g->GetEY();
  exo = g->GetEX();
  for(Int_t i=0;i < n;i++){
    x[i] = xo[i]/nq;
    y[i] = yo[i]/nq;
    ey[i] = eyo[i]/nq;
    ex[i] = exo[i]/2;
  }

  TGraphErrors *g2 = new TGraphErrors(n,x,y,ex,ey);

  g2->SetMarkerStyle(g->GetMarkerStyle());
  g2->SetMarkerColor(g->GetMarkerColor());
  g2->SetLineColor(g->GetLineColor());
  
  return g2;
}

TGraphAsymmErrors* ptscaling(TGraphAsymmErrors *g,Int_t nq){
  Int_t n = g->GetN();
  Double_t x[100],y[100],ex[100],ey[100],ex2[100],ey2[100];
  Double_t *xo,*yo,*eyo,*eyo2,*exo;
  xo = g->GetX();
  yo = g->GetY();
  exo = g->GetEXlow();
  eyo = g->GetEYlow();
  eyo2 = g->GetEYhigh();
  for(Int_t i=0;i < n;i++){
    x[i] = xo[i]/nq;
    y[i] = yo[i]/nq;
    ey[i] = eyo[i]/nq;
    ey2[i] = eyo2[i]/nq;
    ex[i] = exo[i]/2;
    ex2[i] = exo[i]/2;
  }

  TGraphAsymmErrors *g2 = new TGraphAsymmErrors(n,x,y,ex,ex2,ey,ey2);

  g2->SetMarkerStyle(g->GetMarkerStyle());
  g2->SetMarkerColor(g->GetMarkerColor());
  g2->SetLineColor(g->GetLineColor());
  
  return g2;
}

TGraphErrors* mtscaling(TGraphErrors *g,Float_t mass,Int_t nq){
  Int_t n = g->GetN();
  Double_t x[100],y[100],ex[100],ey[100];
  Double_t *xo,*yo,*eyo,*exo;
  xo = g->GetX();
  yo = g->GetY();
  eyo = g->GetEY();
  exo = g->GetEX();
  for(Int_t i=0;i < n;i++){
    x[i] = (TMath::Sqrt(xo[i]*xo[i] + mass*mass) - mass)/nq;
    y[i] = yo[i]/nq;
    ey[i] = eyo[i]/nq;
    ex[i] = exo[i]/2;
  }

  TGraphErrors *g2 = new TGraphErrors(n,x,y,ex,ey);
  
  g2->SetMarkerStyle(g->GetMarkerStyle());
  g2->SetMarkerColor(g->GetMarkerColor());
  g2->SetLineColor(g->GetLineColor());

  return g2;
}

TGraphAsymmErrors* mtscaling(TGraphAsymmErrors *g,Float_t mass,Int_t nq){
  Int_t n = g->GetN();
  Double_t x[100],y[100],ex[100],ey[100],ex2[100],ey2[100];
  Double_t *xo,*yo,*eyo,*eyo2,*exo;
  xo = g->GetX();
  yo = g->GetY();
  eyo = g->GetEYlow();
  eyo2 = g->GetEYhigh();
  exo = g->GetEXhigh();
  for(Int_t i=0;i < n;i++){
    x[i] = (TMath::Sqrt(xo[i]*xo[i] + mass*mass) - mass)/nq;
    y[i] = yo[i]/nq;
    ey[i] = eyo[i]/nq;
    ey2[i] = eyo2[i]/nq;
    ex[i] = exo[i]/2;
    ex2[i] = exo[i]/2;
  }

  TGraphAsymmErrors *g2 = new TGraphAsymmErrors(n,x,y,ex,ex2,ey,ey2);
  
  g2->SetMarkerStyle(g->GetMarkerStyle());
  g2->SetMarkerColor(g->GetMarkerColor());
  g2->SetLineColor(g->GetLineColor());

  return g2;
}

Float_t nonflow(Int_t ic,Float_t pt){ // non flow systematics [%]
  Float_t mult[] = {55,47.5,35,22,14,8,5,3,2};

  return (pt/300 + 0.01)*TMath::Sqrt(mult[ic]/35);
}

Float_t NUO(Int_t ic,Float_t pt){ // dE/dx modulation [%]
  if(pt < 2) return 0.0;
  if(pt > 3) return 0.06;

  return (pt - 2)*0.06;
}

Float_t PIDpi(Int_t ic,Float_t pt){ // syst pi PID [%]
  if(pt < 3) return 0.0;

  return 0.02/6*pt;
}

Float_t PIDka(Int_t ic,Float_t pt){ // syst K PID [%]

  return 0.004*(pt-2.3)*(pt-2.3);
}

Float_t PIDpr(Int_t ic,Float_t pt){ // syst p PID [%]

  return 0.004*(pt-2.3)*(pt-2.3);
}

Float_t NUOtofall(Int_t ic,Float_t pt){ // syst pi PID NUO TOF [%]
  if(pt < 0.4) return 0;
  return 0.01;
}

Float_t NUOtofka(Int_t ic,Float_t pt){ // syst K PID NUO TOF [%]
  if(pt < 0.4) return 0;
  return 0.01/pt/pt;
}

Float_t NUOtofpr(Int_t ic,Float_t pt){ // syst p PID NUO TOF [%]
  if(pt < 0.8) return 0;
  return 0.03/pt/pt;
}

Float_t FeedDown(Int_t ic,Float_t pt){ // syst p PID NUO TOFfeed down [%]
  return 0.05/pt/pt;
}

Float_t systPi(Int_t ic,Float_t pt){
  Float_t e0 = 0.005; // event plane res
  Float_t e1 = 0;//nonflow(ic,pt);
  Float_t e2 = NUO(ic,pt);
  Float_t e3 = PIDpi(ic,pt);
  Float_t e4 = NUOtofall(ic,pt);

  return TMath::Sqrt(e0*e0+e1*e1+e2*e2+e3*e3+e4*e4 + 0.02*0.02);
}

Float_t systKa(Int_t ic,Float_t pt){
  Float_t e0 = 0.005; // event plane res
  Float_t e1 = 0;//nonflow(ic,pt);
  Float_t e2 = NUO(ic,pt);
  Float_t e3 = PIDka(ic,pt);
  Float_t e4 = NUOtofall(ic,pt);
  Float_t e5 = NUOtofka(ic,pt);

  return TMath::Sqrt(e0*e0+e1*e1+e2*e2+e3*e3+e4*e4+e5*e5 + 0.02*0.02);
}

Float_t systPr(Int_t ic,Float_t pt){
  Float_t e0 = 0.005; // event plane res
  Float_t e1 = 0;//nonflow(ic,pt);
  Float_t e2 = NUO(ic,pt);
  Float_t e3 = PIDpr(ic,pt);
  Float_t e4 = NUOtofall(ic,pt);
  Float_t e5 = NUOtofpr(ic,pt);
  Float_t e6 = FeedDown(ic,pt);

  return TMath::Sqrt(e0*e0+e1*e1+e2*e2+e3*e3+e4*e4+e5*e5+e6*e6 + 0.02*0.02);
}

TGraphErrors* v3Pion0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00336251, 0.0042964, 0.00591363, 0.00745722, 0.00917811, 0.0114605, 0.0140557, 0.0179659, 0.0219244, 0.0262026, 0.0296011, 0.0327218, 0.036855, 0.0404912, 0.0438127, 0.045798, 0.0512975, 0.0549897, 0.059595, 0.058824, 0.0625847, 0.0662299, 0.0730445, 0.073744, 0.0756572, 0.0801984, 0.0794842, 0.0772571, 0.06879, 0.0540573, 0.0427378, 0.0531661, 0.0420874, 0.0124048, 0.0326177, 0.0529338, 0.080694, 0.00958407, -0.0393533, -0.0185034, 0.0631424};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000210522, 0.000206387, 0.000214394, 0.000242337, 0.000275482, 0.000302358, 0.00018653, 0.000199449, 0.000219671, 0.000242186, 0.000265206, 0.000305424, 0.000353422, 0.000419006, 0.000494644, 0.000562017, 0.00074144, 0.000805938, 0.000880169, 0.000956853, 0.00103696, 0.000785844, 0.00102401, 0.00134318, 0.00173053, 0.00221324, 0.00278667, 0.00347647, 0.00426881, 0.00398661, 0.00557066, 0.0074162, 0.00947343, 0.0116492, 0.0142187, 0.0127081, 0.0171297, 0.0219476, 0.0216085, 0.028262, 0.0392477};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.23867300000000003e-03;
	Float_t pol1 = 6.71160999999999980e-04;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(0, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0015475, 0.00157659, 0.00184087, 0.00126819, 0.00306379, 0.00554225, 0.00878697, 0.0108394, 0.0140058, 0.0198785, 0.0221055, 0.029062, 0.0333972, 0.0361533, 0.0402771, 0.0465897, 0.0483209, 0.0515488, 0.0555625, 0.0597794, 0.0641036, 0.0744021, 0.0742929, 0.0796722, 0.0763208, 0.074773, 0.0773214, 0.107016, 0.0676339, 0.061203, 0.0238502, 0.0728963, -0.00581352, 0.136265, 0.23944, 0.351119, -0.296351, 0.224582, -0.408449, -0.145808};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00156347, 0.00127827, 0.00124381, 0.00157615, 0.0020528, 0.000885764, 0.000846416, 0.00088859, 0.000918216, 0.000933995, 0.000962436, 0.00101051, 0.00106814, 0.0011348, 0.00121794, 0.00128274, 0.00141546, 0.00158458, 0.00179081, 0.00203398, 0.00158186, 0.00212304, 0.00284786, 0.00386322, 0.0051468, 0.00683846, 0.0090759, 0.0119298, 0.0125558, 0.0220895, 0.0390313, 0.064347, 0.106476, 0.156529, 0.201189, 0.401203, 0.546208, 0.484561, 0.692492, 1.58935};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 1.23867300000000003e-03;
	  Float_t pol1 = 6.71160999999999980e-04;
	  
	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = _y[i]*systKa(0, _x[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.00967204, 0.00535552, -0.000611655, 0.00156903, -0.000214743, 0.000930684, 0.00628733, 0.00346132, 0.00752828, 0.00776737, 0.00689048, 0.0117864, 0.0160731, 0.0216243, 0.0242861, 0.0317626, 0.036016, 0.0461436, 0.0447851, 0.0532914, 0.0656139, 0.0774954, 0.093497, 0.0932595, 0.097548, 0.104816, 0.103702, 0.104301, 0.113606, 0.0802858, 0.0906418, 0.0615356, 0.0523877, -0.0314623, 0.295411, 0.364149, 0.426971, 0.252544, 1.09528};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00507777, 0.00409405, 0.00348366, 0.00320362, 0.00211831, 0.00200809, 0.00152647, 0.00138163, 0.00133943, 0.00133443, 0.00143251, 0.00154796, 0.00168496, 0.00181804, 0.00188818, 0.00198545, 0.00210711, 0.00225372, 0.00242099, 0.00175226, 0.00216888, 0.00270259, 0.00337084, 0.00294624, 0.00377967, 0.00485491, 0.00629107, 0.00645914, 0.0108552, 0.017606, 0.0278862, 0.0433882, 0.0699234, 0.0828253, 0.157932, 0.276799, 0.359305, 0.35556, 0.390246};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.69829837126647719e-03;
	Float_t pol1 = 1.21632530760273721e-03;
 	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(0, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00420951, 0.00533606, 0.00727045, 0.00880758, 0.0112568, 0.0140611, 0.0168202, 0.0213167, 0.0259051, 0.0307, 0.0352693, 0.0392127, 0.0438769, 0.0474563, 0.0511212, 0.0543492, 0.061692, 0.0649522, 0.0673741, 0.0720042, 0.0725524, 0.0783617, 0.0818539, 0.0921856, 0.0887707, 0.0907123, 0.0907136, 0.0818291, 0.0862718, 0.0811945, 0.0746985, 0.0519574, 0.0399439, 0.0504592, 0.0295955, 0.0107689, 0.0380752, 0.0217176, 0.0331197, 0.0475655, 0.0421606};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000228508, 0.000224354, 0.000233169, 0.000264022, 0.000300728, 0.000330391, 0.000201905, 0.000215871, 0.000237699, 0.000262189, 0.000287106, 0.000330327, 0.000382108, 0.000452355, 0.000533182, 0.0006056, 0.000799207, 0.000868309, 0.00094664, 0.00102796, 0.00111185, 0.00083961, 0.0010888, 0.00141987, 0.00181911, 0.00231573, 0.0029012, 0.0035933, 0.00441343, 0.00407112, 0.00565933, 0.00751697, 0.00952897, 0.0118462, 0.0142275, 0.013119, 0.0175258, 0.0226027, 0.0220753, 0.02934, 0.0409402};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.34651999999999992e-03;
	Float_t pol1 = 7.29597999999999922e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
     Float_t systerr = _y[i]*systPi(1, _x[i]);
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, -0.00139525, 0.000403738, 0.00194544, 0.000362192, 0.00255818, 0.00623812, 0.00852105, 0.0157652, 0.0178431, 0.0237543, 0.0294554, 0.0330147, 0.0386089, 0.0435374, 0.0484218, 0.052563, 0.0576409, 0.063316, 0.0665299, 0.0712927, 0.0704316, 0.0830872, 0.0888531, 0.0855421, 0.0955255, 0.100167, 0.108107, 0.0858922, 0.0867649, 0.102963, 0.102207, 0.0689617, 0.137493, 0.232975, 0.014938, 0.109416, -0.11437, -0.0219904, -1.05656, 0.480083};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00168914, 0.00138095, 0.00134266, 0.00169063, 0.00220897, 0.000961223, 0.000916024, 0.00096254, 0.000996154, 0.00101286, 0.0010432, 0.00109366, 0.00115359, 0.00122476, 0.00131136, 0.00138, 0.00151985, 0.00169741, 0.00191603, 0.00217312, 0.00168575, 0.0022548, 0.00301049, 0.00406655, 0.00539095, 0.00712795, 0.00931687, 0.0123094, 0.0128226, 0.0224061, 0.0389769, 0.0662474, 0.103818, 0.157602, 0.189825, 0.433172, 0.5729, 0.481099, 0.629613, 2.73969};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

   if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.34651999999999992e-03;
	Float_t pol1 = 7.29597999999999922e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(1, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00179347, -0.00174397, 0.00562031, -0.000854717, 0.00102079, 0.00268305, 0.00484297, 0.0095663, 0.0100514, 0.0122446, 0.010814, 0.0159408, 0.0223776, 0.0276956, 0.0328181, 0.0362557, 0.0455214, 0.049074, 0.0556873, 0.0661178, 0.0797023, 0.0959318, 0.10964, 0.112015, 0.118346, 0.119038, 0.140106, 0.14157, 0.10431, 0.116439, 0.10492, 0.0871708, 0.0568561, -0.0455202, 0.198776, 0.426306, 0.107729, 0.370331, 0.0994008};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00538956, 0.00437022, 0.00372889, 0.00343717, 0.00227295, 0.00215975, 0.0016284, 0.00148015, 0.00143673, 0.00143295, 0.00154079, 0.00166453, 0.00181946, 0.00196447, 0.00204408, 0.00215335, 0.00228741, 0.00244586, 0.00262781, 0.00190225, 0.00235167, 0.00292554, 0.0036415, 0.00316742, 0.0040373, 0.00516761, 0.00665627, 0.00678022, 0.0111978, 0.0183614, 0.0288174, 0.0448455, 0.0723913, 0.0836158, 0.15176, 0.3484, 0.344686, 0.400118, 0.447911};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

 if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.84616637865581655e-03;
	Float_t pol1 = 1.32222872399611561e-03;
 	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(1, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00461794, 0.00642023, 0.00852456, 0.0103066, 0.0130283, 0.0157136, 0.0194803, 0.0248819, 0.0296363, 0.0344324, 0.0398318, 0.0442839, 0.0488051, 0.0535568, 0.0587408, 0.0608659, 0.0691928, 0.0715776, 0.0753062, 0.0791578, 0.0809953, 0.0870097, 0.0914991, 0.0965167, 0.100661, 0.0974432, 0.102311, 0.0922769, 0.0863138, 0.0857585, 0.0742174, 0.058852, 0.0443585, 0.0521248, 0.0365875, 0.0124116, 0.0456829, 0.00609852, 0.0251501, 0.0762395, 0.0729517};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000192925, 0.00018977, 0.000197506, 0.000224275, 0.000256092, 0.000281821, 0.00017049, 0.000182355, 0.000200813, 0.000221589, 0.000243065, 0.00027971, 0.000324058, 0.000382571, 0.000450527, 0.00051265, 0.000674555, 0.000731225, 0.000795309, 0.00086146, 0.000932844, 0.000702476, 0.000905643, 0.0011698, 0.00149335, 0.00189126, 0.00235292, 0.0029003, 0.00353017, 0.00324194, 0.00448472, 0.00592329, 0.00753801, 0.00932365, 0.0111649, 0.0103037, 0.0139087, 0.0180768, 0.0179247, 0.0232962, 0.0328189};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0; 
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.36284000000000027e-04;
	Float_t pol1 = 8.43209000000000037e-04;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00103604, -0.000214098, 0.00273713, 0.00117868, 0.00364223, 0.00963693, 0.0110632, 0.013632, 0.0215956, 0.0263831, 0.0334408, 0.0397796, 0.0429066, 0.0491487, 0.0559746, 0.059679, 0.0665367, 0.0664515, 0.0734363, 0.0706062, 0.0827271, 0.0854827, 0.0950996, 0.0986797, 0.0978439, 0.0942852, 0.0926222, 0.0934745, 0.107198, 0.0728447, 0.0950053, 0.044093, 0.013894, -0.247519, 0.362324, -0.215948, -0.210557, -0.496652, 0.187874, 1.75622};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00141303, 0.00115861, 0.00112637, 0.00141086, 0.00183742, 0.000797225, 0.000770476, 0.000812609, 0.000843751, 0.000858169, 0.000884989, 0.000926432, 0.000975764, 0.0010329, 0.00110305, 0.00115771, 0.00127264, 0.00141801, 0.00159964, 0.00181237, 0.0014013, 0.00186508, 0.00248639, 0.00331328, 0.00435613, 0.00572852, 0.00745468, 0.00969517, 0.00999537, 0.0169686, 0.0289983, 0.0484386, 0.0833753, 0.129986, 0.161529, 0.305458, 0.448385, 0.512319, 0.50547, 0.546504};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.36284000000000027e-04;
	Float_t pol1 = 8.43209000000000037e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.0162324, -0.00252875, -5.90144e-06, -0.00476797, -0.00229146, 0.00398366, 0.00388682, 0.0067282, 0.00975148, 0.0139623, 0.0196024, 0.0240234, 0.0274354, 0.0338669, 0.0383077, 0.0439452, 0.0512254, 0.0570444, 0.0664354, 0.0752314, 0.0935012, 0.100361, 0.121053, 0.122332, 0.132495, 0.14685, 0.13494, 0.144166, 0.147784, 0.143769, 0.119115, 0.0753315, 0.101949, 0.109027, -0.0198193, 0.195004, 0.532964, 0.159071, 0.304631};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00439881, 0.00358002, 0.00308529, 0.00284691, 0.00188604, 0.00179879, 0.00134515, 0.00121996, 0.00118715, 0.00119346, 0.00128806, 0.00139995, 0.00153443, 0.00166203, 0.00173607, 0.00183116, 0.00194285, 0.00208101, 0.00223829, 0.00162341, 0.0020054, 0.00249063, 0.00308347, 0.00267283, 0.00338056, 0.00430813, 0.00549979, 0.00551871, 0.00893073, 0.014277, 0.0224559, 0.0335484, 0.0527298, 0.0616519, 0.128009, 0.242888, 0.324704, 0.327696, 0.410988};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.13364693747489361e-03;
	Float_t pol1 = 1.52812297971200361e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00519639, 0.00686685, 0.00875694, 0.0115813, 0.0149347, 0.0174483, 0.0218948, 0.0279923, 0.0332986, 0.0395437, 0.0448705, 0.0495976, 0.0549734, 0.0596968, 0.0646702, 0.0694171, 0.0761951, 0.0800467, 0.0820417, 0.0882722, 0.0909226, 0.0949481, 0.0998517, 0.0981284, 0.103912, 0.105908, 0.0958954, 0.0959066, 0.091531, 0.0858565, 0.0824873, 0.0789383, 0.0546646, 0.0323943, 0.0398512, 0.058909, 0.0315019, 0.0458475, 0.0152419, 0.0401339, 0.00171052};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000264306, 0.000261423, 0.000272779, 0.000310985, 0.000356496, 0.000393309, 0.000235694, 0.000252646, 0.000278736, 0.000308116, 0.000338422, 0.0003898, 0.000452198, 0.000538239, 0.000644435, 0.000738943, 0.000937402, 0.00101476, 0.0010981, 0.00118455, 0.00128004, 0.000961587, 0.00123273, 0.00158511, 0.00201169, 0.00249688, 0.00308622, 0.00378193, 0.00456781, 0.00416223, 0.00567462, 0.00744567, 0.00947753, 0.0117554, 0.014131, 0.0128392, 0.0175674, 0.0227266, 0.022365, 0.0298124, 0.0421548};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.90101399999999992e-03;
	Float_t pol1 = 1.03004500000000005e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(3, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00380575, 0.00351231, 0.00042128, 0.00337534, 0.00625768, 0.0112856, 0.0125639, 0.0191966, 0.0242974, 0.0322568, 0.0358068, 0.0447696, 0.0519689, 0.054212, 0.0615544, 0.0622953, 0.0705619, 0.0726914, 0.0826837, 0.0818793, 0.09019, 0.090596, 0.0953699, 0.100226, 0.108469, 0.102736, 0.0902474, 0.0857363, 0.0690082, 0.0664356, 0.0393149, -0.00799427, 0.0412672, 0.303438, 0.00582671, 0.820174, -0.580176, 1.08132, 1.61094, 0.836442};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00191587, 0.00157322, 0.00153987, 0.00194804, 0.00249873, 0.00108422, 0.00105534, 0.00112328, 0.00117285, 0.00119472, 0.00123609, 0.00129458, 0.00136514, 0.00144791, 0.00154846, 0.00162865, 0.00178834, 0.00198891, 0.0022366, 0.00252929, 0.00194906, 0.00257001, 0.0033886, 0.00447902, 0.00597657, 0.0077999, 0.010111, 0.0131682, 0.0136283, 0.0232701, 0.040096, 0.0695593, 0.11547, 0.188365, 0.224126, 0.453195, 0.735385, 0.612546, 0.702625, 0.84789};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.90101399999999992e-03;
	Float_t pol1 = 1.03004500000000005e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(3, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.00162998, 0.00135553, 0.00268008, -2.88326e-05, 0.00717294, 0.00338302, 0.0098108, 0.0128466, 0.0121947, 0.0198757, 0.0265091, 0.0285274, 0.0366512, 0.0396301, 0.0494639, 0.0573915, 0.060973, 0.0692938, 0.0779223, 0.0859542, 0.105306, 0.123225, 0.124215, 0.133615, 0.145517, 0.14513, 0.152834, 0.154076, 0.149022, 0.158324, 0.141536, 0.0908961, 0.15833, 0.0454533, 0.195242, -0.368495, -0.0699533, -0.320681, 1.18474};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00574103, 0.00472569, 0.00408722, 0.00379418, 0.00251462, 0.00241207, 0.00178918, 0.00162871, 0.00159259, 0.00161675, 0.00174712, 0.00191264, 0.00211679, 0.0023059, 0.00242345, 0.00256207, 0.00273003, 0.00293429, 0.00315888, 0.00229873, 0.00285069, 0.00355799, 0.00442805, 0.00378705, 0.00477178, 0.00603747, 0.00764458, 0.00760955, 0.0121272, 0.0190318, 0.0291333, 0.0447506, 0.0639781, 0.0781551, 0.149547, 0.31162, 0.336667, 0.414484, 0.628442};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.60641395006084133e-03;
	Float_t pol1 = 1.86671983155921230e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(3, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00514108, 0.00672007, 0.00967282, 0.012633, 0.0156524, 0.0188224, 0.0237045, 0.0299925, 0.0361635, 0.0422306, 0.0474108, 0.0534242, 0.0594919, 0.0656865, 0.0688581, 0.0724254, 0.0791741, 0.0818226, 0.0878019, 0.088663, 0.0937712, 0.0953295, 0.0966741, 0.10527, 0.0982133, 0.100045, 0.0955392, 0.103438, 0.0861584, 0.0826865, 0.0631282, 0.0609323, 0.0296266, 0.0360597, 0.0527848, 0.0153368, -0.00640239, 0.0210628, -0.113878, -0.0420435, 0.0684772};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.000387513, 0.000382062, 0.000399717, 0.000457563, 0.000526693, 0.000582927, 0.000347414, 0.000373975, 0.000414497, 0.00045968, 0.000506558, 0.000583763, 0.000678084, 0.000797031, 0.000941195, 0.00107603, 0.00138966, 0.00149598, 0.00160715, 0.00172442, 0.00186579, 0.00139376, 0.00178089, 0.00227312, 0.00286831, 0.00355367, 0.00438117, 0.0053019, 0.00637583, 0.00578408, 0.00781772, 0.0102514, 0.0130533, 0.0160274, 0.0196597, 0.0177473, 0.0243562, 0.0317803, 0.0318796, 0.0419259, 0.0606142};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.36398699999999969e-03;
	Float_t pol1 = 1.28090200000000000e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(4, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, -0.000448657, -0.000397745, 0.00223929, 0.0020443, 0.00766057, 0.0147414, 0.0166677, 0.0204658, 0.0273593, 0.0331384, 0.041485, 0.0472078, 0.0543817, 0.063844, 0.0672682, 0.0634426, 0.0772583, 0.0792137, 0.0800184, 0.0878916, 0.0868313, 0.101706, 0.0976556, 0.108331, 0.0940934, 0.0924703, 0.1065, 0.0925592, 0.0792342, 0.0528274, 0.0691005, -0.00104652, 0.0222572, 0.273817, -0.309606, 0.0212317, 0.47615, -0.542112, 2.11997, -1.47365};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00274573, 0.00226854, 0.00224066, 0.00285215, 0.00353177, 0.0015494, 0.0015083, 0.00161636, 0.0017295, 0.00178111, 0.00184845, 0.00193832, 0.00203435, 0.0021395, 0.00227017, 0.00237759, 0.00261493, 0.00292321, 0.00330223, 0.00375077, 0.00290497, 0.0038074, 0.0049577, 0.00644356, 0.00841452, 0.0107389, 0.0136556, 0.0171925, 0.0169879, 0.0269182, 0.0418678, 0.0637165, 0.098592, 0.147218, 0.172146, 0.330455, 0.509198, 0.723025, 1.22814, 1.76288};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.36398699999999969e-03;
	Float_t pol1 = 1.28090200000000000e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(4, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.00134837, 0.0129444, -0.00455304, -0.00508216, 0.00360254, 0.00671058, 0.0103081, 0.0153634, 0.0182269, 0.021164, 0.0300812, 0.0303042, 0.0409975, 0.0510105, 0.0536346, 0.0642107, 0.0722623, 0.0766332, 0.0832146, 0.0961881, 0.114689, 0.117237, 0.12622, 0.137184, 0.150702, 0.152045, 0.15773, 0.17032, 0.119951, 0.14929, 0.0586223, 0.1223, 0.0986238, 0.0153053, 0.133439, 0.254889, -0.0343679, 0.519122, 0.741424};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00793538, 0.00654322, 0.00570034, 0.00530356, 0.00353453, 0.00341486, 0.00254129, 0.00232278, 0.00226879, 0.00231184, 0.00250502, 0.00275635, 0.00307742, 0.00341543, 0.00362731, 0.00385972, 0.00412541, 0.00445048, 0.00482212, 0.00351727, 0.00437345, 0.00543228, 0.00678885, 0.0057682, 0.00720011, 0.00902802, 0.0113323, 0.0111222, 0.0173232, 0.0261736, 0.0385091, 0.0554513, 0.078188, 0.0815409, 0.138491, 0.263223, 0.350662, 0.635402, 0.904071};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 3.24118098903184706e-03;
	Float_t pol1 = 2.32134148520698243e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(4, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00470199, 0.00668197, 0.00976864, 0.0125668, 0.0165945, 0.0168669, 0.0223222, 0.0298682, 0.0368004, 0.0414957, 0.0466974, 0.0505491, 0.0584346, 0.0627035, 0.069766, 0.0730282, 0.0751266, 0.0784967, 0.0821991, 0.0868877, 0.0912784, 0.0915794, 0.0868381, 0.0954055, 0.0915013, 0.0850078, 0.078301, 0.0715817, 0.0778199, 0.0596224, 0.0704722, 0.0530317, 0.0880876, 0.018892, -0.0141094, 0.0397421, 0.0203508, 0.0266409, -0.0898402, 0.13245, 0.0615592};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.00061778, 0.000611419, 0.000642008, 0.00073809, 0.000853287, 0.000947983, 0.000563933, 0.000610585, 0.000681148, 0.000759122, 0.000837411, 0.000960857, 0.0011128, 0.00131825, 0.0015757, 0.0018105, 0.00228684, 0.00245037, 0.00261375, 0.00278801, 0.00300688, 0.00224118, 0.00285223, 0.00359793, 0.00451726, 0.00555495, 0.00680936, 0.00822186, 0.00977507, 0.00881085, 0.0118381, 0.0154212, 0.0192629, 0.0240823, 0.0291823, 0.0264709, 0.0360962, 0.0478299, 0.0477297, 0.0641824, 0.0930879};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 3.02312299999999996e-03;
	Float_t pol1 = 1.63804399999999996e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(5, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, -0.000835722, -0.00111267, 0.00157752, 0.00530462, 0.0048155, 0.0162994, 0.0174401, 0.0175551, 0.0268205, 0.030134, 0.0433559, 0.0506271, 0.0540468, 0.05611, 0.0721392, 0.0651572, 0.0658194, 0.0746101, 0.0681309, 0.0737592, 0.0816147, 0.0780384, 0.0877877, 0.12141, 0.098673, 0.0775791, 0.103087, 0.0484423, 0.131149, 0.121459, -0.00237364, -0.0762335, 0.16375, -0.0902876, 0.255197, 0.320249, 0.501913, 0.241749, 0.982189, 0.100093};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00431504, 0.00356388, 0.00349553, 0.00436078, 0.00555725, 0.00247842, 0.00240435, 0.00259563, 0.00282875, 0.00294149, 0.00306727, 0.00322682, 0.00339458, 0.00358301, 0.0038133, 0.00400849, 0.00441698, 0.00494144, 0.00558109, 0.0063388, 0.00489129, 0.00641331, 0.00839007, 0.0107902, 0.014018, 0.0177179, 0.0222906, 0.0281227, 0.0272413, 0.0424294, 0.0655903, 0.100751, 0.148188, 0.20562, 0.253328, 0.548881, 0.746892, 1.45581, 2.2303, 0};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 3.02312299999999996e-03;
	Float_t pol1 = 1.63804399999999996e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(5, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.00160495, -0.00375938, 0.00225361, -0.00138298, 0.00672974, 0.0146752, 0.0157013, 0.0181489, 0.0126964, 0.0253593, 0.0215832, 0.029294, 0.034828, 0.0494889, 0.0590901, 0.059864, 0.0638435, 0.0902087, 0.0877695, 0.0938779, 0.113694, 0.112781, 0.1215, 0.114081, 0.128978, 0.141657, 0.0967336, 0.113659, 0.0814051, 0.11728, 0.0545098, 0.302757, 0.121751, -0.0440605, -0.285575, -0.236346, -0.156475, -0.66292, 0.0716732};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.0118804, 0.00987968, 0.00863847, 0.00805332, 0.00541803, 0.00525946, 0.00392364, 0.00362233, 0.00357483, 0.00366657, 0.00398395, 0.00440916, 0.00496336, 0.00561059, 0.00604309, 0.00646088, 0.00695276, 0.00752836, 0.0081738, 0.00599113, 0.00750449, 0.00934324, 0.0117725, 0.0100423, 0.0125607, 0.0156953, 0.0197364, 0.0191932, 0.0295074, 0.0447571, 0.0636396, 0.0897913, 0.122813, 0.128055, 0.203871, 0.344707, 0.409889, 0.559407, 0.842927};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

 if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 4.14488622049151260e-03;
	Float_t pol1 = 2.96857730797800597e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(5, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00348915, 0.00756235, 0.0107608, 0.00988755, 0.0121196, 0.0204889, 0.0220417, 0.030525, 0.0316641, 0.0401874, 0.0432457, 0.0534317, 0.0556032, 0.0610604, 0.0638937, 0.0631415, 0.0672388, 0.0790177, 0.0733123, 0.0786108, 0.0841621, 0.0848912, 0.0758076, 0.0703749, 0.0833439, 0.0778629, 0.0876323, 0.0388327, 0.0640147, 0.0516244, 0.0726542, 0.0812565, 0.0449663, -0.0553643, 0.0332307, 0.203045, 0.133275, 0.169177, 0.056018, -0.26421, -0.440226};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.00124936, 0.00123461, 0.0013016, 0.00150388, 0.00174677, 0.00194891, 0.00115884, 0.00126079, 0.00141183, 0.00157959, 0.00175432, 0.00202357, 0.00235303, 0.00278943, 0.00333422, 0.0038491, 0.00485018, 0.00519756, 0.00553165, 0.00591171, 0.00637311, 0.00473058, 0.00597128, 0.00750991, 0.00931667, 0.0114766, 0.0139393, 0.016774, 0.0199227, 0.0177383, 0.0237252, 0.0309059, 0.0389086, 0.0485901, 0.0578131, 0.0533689, 0.0743548, 0.0972631, 0.0985038, 0.13383, 0.197948};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 4.01257199999999984e-03;
	Float_t pol1 = 2.17416299999999979e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(6, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00792801, 0.00266082, 0.00470314, -0.00978977, -0.00565417, 0.0175624, 0.0186389, 0.0245958, 0.0396072, 0.0361056, 0.0526261, 0.0490507, 0.0426577, 0.0705987, 0.0541062, 0.0705978, 0.0609208, 0.07822, 0.0705974, 0.0903856, 0.0910678, 0.0668338, 0.0912261, 0.0751016, 0.0704388, 0.113234, 0.101156, 0.100512, 0.00977196, 0.0365165, 0.0487332, 0.223454, -0.042882, -0.278059, -0.434993, -0.428813, 3.15745, 0.0246815, 5.204, 1.15718};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.00854556, 0.00710258, 0.00701836, 0.00878387, 0.0111126, 0.00502841, 0.00493816, 0.00538809, 0.00590501, 0.00619556, 0.00648759, 0.00683347, 0.00718756, 0.00754521, 0.0079971, 0.00837843, 0.00924004, 0.0103362, 0.0116509, 0.0132779, 0.0102584, 0.0134035, 0.017172, 0.0220316, 0.0277159, 0.0343268, 0.0430676, 0.0526704, 0.0495564, 0.0733204, 0.106265, 0.150402, 0.212986, 0.312029, 0.337181, 0.54166, 1.02797, 1.06409, 2.42013, 6.94904};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 4.01257199999999984e-03;
	Float_t pol1 = 2.17416299999999979e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(6, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.0192241, 0.0260949, 0.0182832, -0.00275147, 0.0187805, 0.00631485, 0.0211281, 0.0259732, 0.0239513, 0.0400317, 0.0491792, 0.0397261, 0.0519205, 0.0519361, 0.0752147, 0.0842076, 0.07053, 0.0909382, 0.0816053, 0.0891785, 0.110612, 0.077322, 0.144557, 0.109386, 0.119593, 0.0787762, 0.0846064, 0.104219, 0.0807165, 0.0144733, 0.269359, -0.0582495, 0.247274, -0.709632, 0.214561, -0.705862, 1.36533, -0.323964, -0.366863};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.0223046, 0.0186193, 0.0164312, 0.0154014, 0.0104109, 0.0102522, 0.00769725, 0.00723415, 0.00717281, 0.00744624, 0.00819774, 0.00914984, 0.0103872, 0.0118351, 0.0128888, 0.0138632, 0.0150433, 0.0163627, 0.0178012, 0.013148, 0.0165023, 0.0206064, 0.0255859, 0.0222043, 0.0275954, 0.0343904, 0.0425941, 0.0410592, 0.0613099, 0.0888349, 0.127822, 0.177487, 0.232908, 0.244446, 0.385938, 0.593035, 0.727183, 0.763612, 1.33937};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 5.50148450680349593e-03;
	Float_t pol1 = 3.94017620709329313e-03;

 	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(6, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.00532019, 0.00818378, 0.00724118, 0.00516368, 0.0128364, 0.0226381, 0.0233119, 0.0242261, 0.0404504, 0.0300591, 0.0367223, 0.0392187, 0.0631186, 0.0534841, 0.0547624, 0.0693861, 0.0707018, 0.0330524, 0.0641408, 0.0326427, 0.0693289, 0.0732374, 0.056526, 0.0186195, 0.13461, -0.00250564, -0.0807247, -0.0230955, 0.0770975, -0.074968, 0.114935, -0.0557285, 0.0529142, -0.115263, -0.265437, 0.025273, 0.163827, 0.357759, -0.314276, 0.815976, 0.15717};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.00354324, 0.00351164, 0.00371981, 0.0043181, 0.00504426, 0.00564883, 0.00336921, 0.00369073, 0.00416028, 0.00467844, 0.00520676, 0.0060233, 0.00700758, 0.00830575, 0.00992362, 0.0114483, 0.0143172, 0.0152843, 0.01618, 0.0171645, 0.018547, 0.0137816, 0.0173176, 0.0216648, 0.0269479, 0.0326113, 0.0392359, 0.0470618, 0.0553329, 0.0496574, 0.0657208, 0.0853741, 0.107096, 0.134301, 0.158415, 0.150224, 0.200109, 0.282743, 0.282242, 0.388098, 0.600491};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 5.49768000000000045e-03;
	Float_t pol1 = 2.97884999999999990e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(7, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.0197483, 0.0180058, -0.0171612, 0.0491509, 0.0332618, 0.0129175, 0.00558811, -0.00699739, 0.0317913, 0.0267255, 0.0401873, 0.0714799, 0.0342715, 0.0109885, 0.0215154, 0.0293005, 0.0669465, 0.0296294, 0.0647696, 0.129599, 0.0300014, 0.113043, 0.0821695, -0.0340911, 0.0442565, 0.165024, -0.0487434, -0.0213451, 0.0183962, 0.17259, 0.6626, -0.785362, 0.671771, 0.448768, -0.0679629, 0.861771, 4.27246, -2.42423, 0.39535, 4.5634};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.0239986, 0.0202105, 0.0206189, 0.0280924, 0.034233, 0.0146649, 0.0143844, 0.0157827, 0.0175924, 0.0185466, 0.0195311, 0.0206602, 0.0217104, 0.0228648, 0.0242888, 0.025484, 0.0282024, 0.0315494, 0.035468, 0.0404523, 0.0311365, 0.0404745, 0.0517413, 0.0641674, 0.082981, 0.103051, 0.126052, 0.154289, 0.144323, 0.217221, 0.303621, 0.412496, 0.605089, 0.792056, 0.900304, 1.79174, 2.68906, 5.38102, 0, 12.1419};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
     //printf("%f %f %f\n",xrb[i],yrb[i],erb[i]);
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 5.49768000000000045e-03;
	Float_t pol1 = 2.97884999999999990e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(7, _x[i]);
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, -0.0158129, 0.0656096, 0.0124522, 0.00267904, -0.024072, -0.0102232, 0.0401282, 0.0782716, 0.0132832, 0.0439476, 0.0226909, 0.0658784, 0.0467835, 0.0536947, 0.0234852, 0.0946235, 0.0983099, 0.0460588, 0.098392, 0.127964, 0.0978314, 0.268099, 0.138069, 0.18096, 0.0490962, -0.0279796, 0.33413, 0.312404, 0.0128821, 0.41463, 0.321723, 0.166142, -0.92738, -1.15805, -1.08636, -1.1947, 1.47742, -1.20375, 1.79054};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.0586293, 0.0494651, 0.0436259, 0.0414706, 0.0282552, 0.0281392, 0.0213093, 0.0203187, 0.0205207, 0.0215662, 0.0237977, 0.0268047, 0.0308095, 0.0356603, 0.0391952, 0.042407, 0.0463146, 0.0506837, 0.0550799, 0.0410941, 0.0513746, 0.0644241, 0.0820897, 0.0689592, 0.0852705, 0.105254, 0.130132, 0.12491, 0.185621, 0.265899, 0.381963, 0.506455, 0.65307, 0.689405, 0.991451, 1.50079, 1.63196, 3.34345, 4.67234};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.53765601848897532e-03;
	Float_t pol1 = 5.39848705646177145e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(7, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Pion7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, -0.0277715, 0.0200781, -0.147342, -0.0910199, 0.138053, -0.135879, 0.13, 0.00721135, 0.169047, 0.118246, 0.150001, 0.154404, 0.0494712, 0.289205, 0.0683699, 0.221813, 0.323357, 0.723906, 0.195442, 0.0600025, -0.551415, -0.0983314, -0.0333938, 0.391585, -0.0025392, 0.751968, -0.266348, -0.381468, 1.40929, -1.07184, -0.402591, 1.12346, -0.341909, 0.393444, 0.544109, 0.695897, -6.6012, -5.45095, -2.99495, -5.37993, -14.2584};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.0585941, 0.0578992, 0.0615916, 0.0718761, 0.0843703, 0.0949434, 0.0568047, 0.0626236, 0.0711161, 0.0805044, 0.0899988, 0.103695, 0.120911, 0.143206, 0.170331, 0.195689, 0.24354, 0.261, 0.274216, 0.291895, 0.315833, 0.234549, 0.294646, 0.363354, 0.446777, 0.548363, 0.657545, 0.786572, 0.929265, 0.817846, 1.09839, 1.40458, 1.78799, 2.15641, 2.71042, 2.50982, 3.42489, 4.76495, 4.69896, 6.3885, 11.134};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 8.29529000000000030e-03;
	Float_t pol1 = 5.60073999999999989e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(8, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Kaon7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, -0.0288316, -0.0592418, 0.369199, 0.575608, 0.0766474, -0.1637, 0.278165, 0.177826, -0.0377868, 0.63419, 0.0407312, -0.0816281, 0.255651, -0.0857342, -0.758954, 0.43096, 0.361551, 0.228773, 0.0348557, 1.36953, -0.264354, 0.774332, -0.0470566, 0.259239, 1.14279, -0.28058, -1.18877, -1.3862, 0.258608, -4.08418, -3.95216, 2.38053, 1.97998, -7.55988, 3.81951, -12.191, -18.1662, 39.2848, -0.47591, 2.15941};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.396338, 0.330962, 0.320101, 0.378465, 0.517482, 0.248631, 0.244079, 0.266511, 0.301912, 0.324328, 0.34321, 0.361752, 0.380482, 0.395796, 0.417082, 0.440378, 0.482559, 0.541749, 0.61127, 0.699858, 0.539332, 0.711999, 0.905082, 1.15617, 1.39784, 1.72269, 2.12391, 2.55149, 2.41808, 3.44732, 4.75497, 7.01377, 9.18817, 11.0011, 12.3568, 20.7339, 35.7767, 47.8611, 0, 0};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 8.29529000000000030e-03;
	Float_t pol1 = 5.60073999999999989e-03;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(8, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3Antiproton7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1,Int_t rebin=1)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.781035, 0.265145, -0.629633, 0.067364, -0.224935, -0.248775, -0.0243976, -0.0547503, -0.206082, -0.320101, 0.0315955, 0.472627, -1.13402, 0.461416, 0.0404126, 0.763037, -0.155347, 0.107662, -0.543443, 0.176083, 0.833912, -0.640997, -0.489657, 1.29716, -0.302417, -0.145165, 3.6357, -2.42992, 1.71601, 9.63302, -5.31106, -2.49473, 15.5048, -2.05208, 28.5906, -63.943, -57.8212, 13.7189, -51.4378};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.907785, 0.767146, 0.684745, 0.650857, 0.450223, 0.454094, 0.348794, 0.337984, 0.345745, 0.366159, 0.405809, 0.461621, 0.529708, 0.62245, 0.700744, 0.759008, 0.834023, 0.909172, 1.00179, 0.743342, 0.938799, 1.16354, 1.42633, 1.24775, 1.5373, 1.85767, 2.31073, 2.13039, 3.05418, 4.28533, 6.58975, 8.00336, 9.78783, 10.9298, 15.4178, 23.7608, 22.9556, 34.1062, 44.7479};

  // rebinning
  Int_t npRB = 0;
  Double_t xrb[100];
  Double_t yrb[100];
  Double_t erb[100];

  if(rebin > 1){
    for(Int_t i=first;i<=last;){
      Double_t sumw = 0;
      xrb[npRB] = 0;
      yrb[npRB] = 0;
      erb[npRB] = 0;
      for(Int_t j=0;j<rebin&&i <=last;j++){
	Double_t weight = _yerr[i] * _yerr[i];
	if(weight > 0){
	  weight = 1./weight;
	  sumw += weight;
	  xrb[npRB] += weight * _x[i];
	  yrb[npRB] += weight * _y[i];
	  erb[npRB] += weight * weight * _yerr[i] * _yerr[i];
	}
	i++;
      }
      if(sumw > 0){
	xrb[npRB] /= sumw;
	yrb[npRB] /= sumw;
	erb[npRB] = TMath::Sqrt(erb[npRB]) / sumw;
      }
      npRB++;
    }

    _nPoints = npRB;
    first = 0;
    last = npRB-1;
    for(Int_t i=0; i < npRB;i++){
      _x[i] = xrb[i];
      _y[i] = yrb[i];
      _yerr[i] = erb[i];
    }
  }

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
       _xerr[i] = 0.05;
   }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 1.13733682403531055e-02;
	Float_t pol1 = 8.14563321585857515e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(8, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr + nonflow*nonflow);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

Float_t systPiV3(Int_t ic,Float_t pt){
  Float_t e0 = 0.04; // event plane res
  Float_t e1 = 0;//nonflow(ic,pt);
  Float_t e2 = 0;//NUO(ic,pt);
  Float_t e3 = PIDpi(ic,pt);
  Float_t e4 = NUOtofall(ic,pt);

  return TMath::Sqrt(e0*e0+e1*e1+e2*e2+e3*e3+e4*e4);
}

Float_t systKaV3(Int_t ic,Float_t pt){
  Float_t e0 = 0.04; // event plane res
  Float_t e1 = 0;//nonflow(ic,pt);
  Float_t e2 = 0;//NUO(ic,pt);
  Float_t e3 = PIDka(ic,pt);
  Float_t e4 = NUOtofall(ic,pt);
  Float_t e5 = NUOtofka(ic,pt);

  return TMath::Sqrt(e0*e0+e1*e1+e2*e2+e3*e3+e4*e4+e5*e5);
}

Float_t systPrV3(Int_t ic,Float_t pt){
  Float_t e0 = 0.04; // event plane res
  Float_t e1 = 0;//nonflow(ic,pt);
  Float_t e2 = 0;//NUO(ic,pt);
  Float_t e3 = PIDpr(ic,pt);
  Float_t e4 = NUOtofall(ic,pt);
  Float_t e5 = NUOtofpr(ic,pt);
  Float_t e6 = FeedDown(ic,pt);

  return TMath::Sqrt(e0*e0+e1*e1+e2*e2+e3*e3+e4*e4+e5*e5+e6*e6);
}


// phi v3 in 10-20% with LHC11h data
TGraphErrors* v3QC2Phi1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1){
  Double_t _x[ ] = {0.6, 1.2, 1.8, 2.4, 3, 4, 5};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const Double_t v3QC2[] = {0.0102678, 0.0124233, 0.048642, 0.0545929, 0.0748703, 0.0682869 };
  const Double_t v3QC2err[] = {0.0059634, 0.00477554, 0.00424124, 0.00477529, 0.00597224, 0.0122696 };
  Int_t _nPoints = 6;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &v3QC2[first], &_xerr[first], &v3QC2err[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;

}

TGraphErrors* v3SPPhi1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1){
  const Double_t v3SP[] = {0.0031926, 0.00884663, 0.046138, 0.0569548, 0.0814731, 0.0734251};
  const Double_t v3SPerr[] = {0.0125239, 0.0101962, 0.00935323, 0.0106264, 0.0133925, 0.0277158 };
  Double_t _x[ ] = {0.6, 1.2, 1.8, 2.4, 3, 4, 5};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //commentme
  Int_t _nPoints = 6;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &v3SP[first], &_xerr[first], &v3SPerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2PionAlex0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0572121, 0.0503675, 0.0375937, 0.0268805, 0.0266244, 0.00924442, 0.0159123, 0.0441546, -0.0248131, 0.0214811, 0.0730123};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00201358, 0.00304022, 0.00437662, 0.0060147, 0.00619144, 0.00936874, 0.0130334, 0.0137247, 0.0214058, 0.0253415, 0.0478934};
  Double_t _yerrH[] = {0.00201358, 0.00304022, 0.00437662, 0.0060147, 0.00619144, 0.00936874, 0.0130334, 0.0137247, 0.0214058, 0.0253415, 0.0478934};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 1.23867300000000003e-03;
  Float_t pol1 = 6.71160999999999980e-04;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.100286, 0.0951594, 0.0870573, 0.0786874, 0.0653645, 0.0541408, 0.0457221, 0.0482269, 0.0519502, 0.0084316, 0.01122};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.0015204, 0.00226955, 0.00323721, 0.00443409, 0.00452888, 0.00689104, 0.00966223, 0.0102573, 0.0160961, 0.019199, 0.0347209};
  Double_t _yerrH[] = {0.0015204, 0.00226955, 0.00323721, 0.00443409, 0.00452888, 0.00689104, 0.00966223, 0.0102573, 0.0160961, 0.019199, 0.0347209};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 1.34651999999999992e-03;
  Float_t pol1 = 7.29597999999999922e-04;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.148558, 0.13737, 0.127354, 0.107954, 0.0923903, 0.0770206, 0.0654796, 0.0618298, 0.064932, 0.0533334, 0.0544499};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00101876, 0.00150302, 0.00211816, 0.00287069, 0.00290956, 0.00443253, 0.00620384, 0.00660053, 0.0104934, 0.012534, 0.02299};
  Double_t _yerrH[] = {0.00101876, 0.00150302, 0.00211816, 0.00287069, 0.00290956, 0.00443253, 0.00620384, 0.00660053, 0.0104934, 0.012534, 0.02299};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 7.36284000000000027e-04;
  Float_t pol1 = 8.43209000000000037e-04;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.184662, 0.172692, 0.158474, 0.133359, 0.121578, 0.0985604, 0.0726206, 0.0936145, 0.0556499, 0.0463445, 0.0218997};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00115285, 0.00167411, 0.00233323, 0.00311742, 0.00316527, 0.00477094, 0.00668583, 0.00706099, 0.0114757, 0.0134119, 0.0255548};
  Double_t _yerrH[] = {0.00115285, 0.00167411, 0.00233323, 0.00311742, 0.00316527, 0.00477094, 0.00668583, 0.00706099, 0.0114757, 0.0134119, 0.0255548};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 1.90101399999999992e-03;
  Float_t pol1 = 1.03004500000000005e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.201648, 0.188309, 0.170016, 0.158705, 0.132905, 0.107812, 0.104791, 0.0941046, 0.0588267, 0.033593, 0.0528952};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.0014898, 0.0021307, 0.00292622, 0.00390727, 0.0038914, 0.00588152, 0.00829931, 0.00880313, 0.0140397, 0.0171821, 0.031647};
  Double_t _yerrH[] = {0.0014898, 0.0021307, 0.00292622, 0.00390727, 0.0038914, 0.00588152, 0.00829931, 0.00880313, 0.0140397, 0.0171821, 0.031647};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 2.36398699999999969e-03;
  Float_t pol1 = 1.28090200000000000e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.20446, 0.185559, 0.177906, 0.163602, 0.141989, 0.127367, 0.105003, 0.100395, 0.110096, 0.11953, 0.0644066};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00220519, 0.00310504, 0.00420546, 0.00550637, 0.00551497, 0.00824085, 0.0116214, 0.0123683, 0.0201189, 0.0243315, 0.044922};
  Double_t _yerrH[] = {0.00220519, 0.00310504, 0.00420546, 0.00550637, 0.00551497, 0.00824085, 0.0116214, 0.0123683, 0.0201189, 0.0243315, 0.044922};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 3.02312299999999996e-03;
  Float_t pol1 = 1.63804399999999996e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.193713, 0.190999, 0.162508, 0.161684, 0.129319, 0.139517, 0.108649, 0.131493, 0.0577825, 0.153995, 0.0963159};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00385242, 0.00531271, 0.00718802, 0.00946356, 0.00929931, 0.0136617, 0.0196584, 0.0210118, 0.0344036, 0.0421224, 0.0798216};
  Double_t _yerrH[] = {0.00385242, 0.00531271, 0.00718802, 0.00946356, 0.00929931, 0.0136617, 0.0196584, 0.0210118, 0.0344036, 0.0421224, 0.0798216};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 4.01257199999999984e-03;
  Float_t pol1 = 2.17416299999999979e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.182112, 0.159495, 0.126087, 0.169947, 0.119372, 0.0960817, 0.121692, 0.080404, 0.0940447, 0.0126981, 0.0345412};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00863382, 0.0117304, 0.0157993, 0.0203184, 0.0202819, 0.0300107, 0.0437674, 0.0456399, 0.0741986, 0.0879486, 0.175398};
  Double_t _yerrH[] = {0.00863382, 0.0117304, 0.0157993, 0.0203184, 0.0202819, 0.0300107, 0.0437674, 0.0456399, 0.0741986, 0.0879486, 0.175398};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 5.49768000000000045e-03;
  Float_t pol1 = 2.97884999999999990e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2PionAlex7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.17601, 0.18709, 0.0865992, 0.180606, 0.157716, 0.123429, 0.294631, 0.275269, 0.238633, 0.299389, 0.21161};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.0256103, 0.0351521, 0.0468882, 0.0601322, 0.0589281, 0.0875188, 0.129327, 0.138888, 0.218522, 0.274498, 0.540804};
  Double_t _yerrH[] = {0.0256103, 0.0351521, 0.0468882, 0.0601322, 0.0589281, 0.0875188, 0.129327, 0.138888, 0.218522, 0.274498, 0.540804};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 8.29529000000000030e-03;
  Float_t pol1 = 5.60073999999999989e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0803428, 0.0795488, 0.0823786, 0.0722152, 0.0553335, 0.079474, 0.0727801, -9.92707e-05, 0.0413631, 0.0473471, -0.0285734};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00239983, 0.00330892, 0.0047863, 0.00697627, 0.00807388, 0.0146574, 0.0229816, 0.0263452, 0.0433039, 0.0502205, 0.0823823};
  Double_t _yerrH[] = {0.00239983, 0.00330892, 0.0047863, 0.00697627, 0.00807388, 0.0146574, 0.0229816, 0.0263452, 0.0433039, 0.0502205, 0.0823823};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 1.69829837126647719e-03;
  Float_t pol1 = 1.21632530760273721e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.152701, 0.154497, 0.144623, 0.149111, 0.138367, 0.10334, 0.103953, 0.123783, 0.152709, 0.0364577, 0.11831};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00193571, 0.0026302, 0.00376345, 0.00542586, 0.00623612, 0.0112432, 0.0177243, 0.0202246, 0.0345635, 0.0402301, 0.0615073};
  Double_t _yerrH[] = {0.00193571, 0.0026302, 0.00376345, 0.00542586, 0.00623612, 0.0112432, 0.0177243, 0.0202246, 0.0345635, 0.0402301, 0.0615073};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 1.84616637865581655e-03;
  Float_t pol1 = 1.32222872399611561e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.215715, 0.220329, 0.213535, 0.20526, 0.179853, 0.164986, 0.144028, 0.121122, 0.08483, 0.0873156, 0.142706};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00134312, 0.0017999, 0.0025414, 0.00364127, 0.00414403, 0.00743048, 0.01166, 0.0136733, 0.0240196, 0.0273348, 0.0425547};
  Double_t _yerrH[] = {0.00134312, 0.0017999, 0.0025414, 0.00364127, 0.00414403, 0.00743048, 0.01166, 0.0136733, 0.0240196, 0.0273348, 0.0425547};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 2.13364693747489361e-03;
  Float_t pol1 = 1.52812297971200361e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.267297, 0.271003, 0.264474, 0.253482, 0.228673, 0.175124, 0.134638, 0.135447, 0.14589, 0.13513, 0.267056};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00158214, 0.00209045, 0.0029266, 0.00409535, 0.00461357, 0.00827324, 0.0130239, 0.0151731, 0.0255994, 0.0331381, 0.0517058};
  Double_t _yerrH[] = {0.00158214, 0.00209045, 0.0029266, 0.00409535, 0.00461357, 0.00827324, 0.0130239, 0.0151731, 0.0255994, 0.0331381, 0.0517058};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 2.60641395006084133e-03;
  Float_t pol1 = 1.86671983155921230e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.29359, 0.299738, 0.285779, 0.261864, 0.227056, 0.187394, 0.167713, 0.167319, 0.0600002, 0.186913, 0.185933};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00212147, 0.00278448, 0.00381831, 0.00531457, 0.0059205, 0.0103991, 0.016121, 0.018671, 0.0338814, 0.0395379, 0.0737746};
  Double_t _yerrH[] = {0.00212147, 0.00278448, 0.00381831, 0.00531457, 0.0059205, 0.0103991, 0.016121, 0.018671, 0.0338814, 0.0395379, 0.0737746};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 3.24118098903184706e-03;
  Float_t pol1 = 2.32134148520698243e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.29195, 0.277526, 0.285899, 0.258206, 0.232818, 0.208127, 0.149166, 0.165656, 0.132402, 0.0834634, 0.139781};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.0032595, 0.00422106, 0.00576192, 0.00787184, 0.00847629, 0.0143821, 0.0225136, 0.0266895, 0.0474944, 0.0639686, 0.097342};
  Double_t _yerrH[] = {0.0032595, 0.00422106, 0.00576192, 0.00787184, 0.00847629, 0.0143821, 0.0225136, 0.0266895, 0.0474944, 0.0639686, 0.097342};
 
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 4.14488622049151260e-03;
  Float_t pol1 = 2.96857730797800597e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.260472, 0.267568, 0.260439, 0.272361, 0.211514, 0.19635, 0.130636, 0.153759, 0.179847, 0.155681, -0.199146};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.00584754, 0.00745765, 0.00994934, 0.0133111, 0.0142698, 0.0240255, 0.0377183, 0.0435218, 0.0758897, 0.0932454, 0.157862};
  Double_t _yerrH[] = {0.00584754, 0.00745765, 0.00994934, 0.0133111, 0.0142698, 0.0240255, 0.0377183, 0.0435218, 0.0758897, 0.0932454, 0.157862};
 
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 5.50148450680349593e-03;
  Float_t pol1 = 3.94017620709329313e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.265336, 0.217832, 0.223959, 0.179655, 0.143289, 0.157325, 0.115741, 0.197337, 0.139223, 0.220918, -1.13372};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.0133475, 0.0168104, 0.0220827, 0.0299551, 0.0316788, 0.0525262, 0.0819623, 0.0973168, 0.159925, 0.207026, 0.346889};
  Double_t _yerrH[] = {0.0133475, 0.0168104, 0.0220827, 0.0299551, 0.0316788, 0.0525262, 0.0819623, 0.0973168, 0.159925, 0.207026, 0.346889};
 
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 7.53765601848897532e-03;
  Float_t pol1 = 5.39848705646177145e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* v2AntiprotonAlex7080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.288452, 0.272661, 0.24954, 0.282237, 0.254303, 0.428618, 0.627356, 0.254672, -0.245046, 1.56461, 0.350998};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrL[] = {0.0411726, 0.0515451, 0.0685134, 0.090632, 0.0937251, 0.155169, 0.242241, 0.298679, 0.660934, 0.787338, 1.0507};
  Double_t _yerrH[] = {0.0411726, 0.0515451, 0.0685134, 0.090632, 0.0937251, 0.155169, 0.242241, 0.298679, 0.660934, 0.787338, 1.0507};
 
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerrL[i] = 0;
      _yerrH[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  Float_t pol0 = 1.13733682403531055e-02;
  Float_t pol1 = 8.14563321585857515e-03;

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t systerr = TMath::Sqrt(0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.003*0.003 + 0.01*0.01 + 0.05*0.05) * _y[i];
	_yerrH[i] = TMath::Sqrt(_yerrH[i]*_yerrH[i] + systerr*systerr);

	systerr = TMath::Sqrt(systerr*systerr + (pol0 + pol1*_x[i])*(pol0 + pol1*_x[i]));
	_yerrL[i] = TMath::Sqrt(_yerrL[i]*_yerrL[i] + systerr*systerr);
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, 
						   &_x[first], 
						   &_y[first], 
						   &_xerr[first], 
						   &_xerr[first], 
						   &_yerrL[first], 
						   &_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);


  return graph;
}

TGraphAsymmErrors* mergeGraph(TGraphErrors *g1,TGraphAsymmErrors *g2,Float_t pt=3)
{
  Int_t np=0;
  Double_t _x[100];
  Double_t _xerr[100];
  Double_t _y[100];
  Double_t _yerrL[100];
  Double_t _yerrH[100];
  if(g1){
    for(Int_t i=0;i < g1->GetN();i++){
      if(g1->GetX()[i] < pt){
	_x[np] = g1->GetX()[i];
	_xerr[np] = g1->GetEX()[i];
	_y[np] = g1->GetY()[i];
	_yerrL[np] = g1->GetEY()[i];
	_yerrH[np] = g1->GetEY()[i];
	np++;
      }
    }
  }

  if(g2){
    for(Int_t i=0;i < g2->GetN();i++){
      if(g2->GetX()[i] > pt){
	    _x[np] = g2->GetX()[i];
	    _xerr[np] = g2->GetEXlow()[i];
	    _y[np] = g2->GetY()[i];
	    _yerrL[np] = g2->GetEYlow()[i];
	    _yerrH[np] = g2->GetEYhigh()[i];
	    np++;
	}
    }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(np, 
						   _x, 
						   _y, 
						   _xerr, 
						   _xerr, 
						   _yerrL, 
						   _yerrH);
  return graph;
}

TGraphAsymmErrors* mergeGraph(TGraphAsymmErrors *g1,TGraphAsymmErrors *g2,Float_t pt=3)
{
    Int_t np=0;
    Double_t _x[100];
    Double_t _xerr[100];
    Double_t _y[100];
    Double_t _yerrL[100];
    Double_t _yerrH[100];

    if(g1){
      for(Int_t i=0;i < g1->GetN();i++){
	if(g1->GetX()[i] < pt){
	  _x[np] = g1->GetX()[i];
	  _xerr[np] = g1->GetEXlow()[i];
	  _y[np] = g1->GetY()[i];
	  _yerrL[np] = g1->GetEYlow()[i];
	  _yerrH[np] = g1->GetEYhigh()[i];
	  np++;
	}
      }
    }      

    if(g2){
      for(Int_t i=0;i < g2->GetN();i++){
	if(g2->GetX()[i] > pt){
	  _x[np] = g2->GetX()[i];
	  _xerr[np] = g2->GetEXlow()[i];
	  _y[np] = g2->GetY()[i];
	  _yerrL[np] = g2->GetEYlow()[i];
	  _yerrH[np] = g2->GetEYhigh()[i];
	  np++;
	}
      }
    }
    
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(np, 
						   _x, 
						   _y, 
						   _xerr, 
						   _xerr, 
						   _yerrL, 
						   _yerrH);
    return graph;
}

TGraphErrors* v3PionAlex0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0739456, 0.0727898, 0.0471398, 0.0381677, 0.0479825, 0.0139417, 0.0179028, 0.0395375, -0.0112012, -0.00243147, 0.0181051};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00311161, 0.00470213, 0.0067746, 0.00927911, 0.00956509, 0.0145144, 0.02021, 0.0211796, 0.0333082, 0.039139, 0.0719018};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 1.23867300000000003e-03;
	  Float_t pol1 = 6.71160999999999980e-04;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.094798, 0.0852534, 0.0775351, 0.0629343, 0.048353, 0.0236994, 0.0133263, -0.00409467, 0.0381892, 0.00921291, 0.125329};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00323655, 0.00483726, 0.00689185, 0.00941221, 0.00965958, 0.0146739, 0.020628, 0.0217555, 0.0338957, 0.0406934, 0.0737968};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 1.34651999999999992e-03;
	  Float_t pol1 = 7.29597999999999922e-04;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0973472, 0.0790525, 0.0917722, 0.0810776, 0.0541923, 0.0569199, 0.0228887, 0.0146595, 0.0613978, 0.0654468, 0.0554691};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0026594, 0.00392309, 0.00551918, 0.007477, 0.0075928, 0.0115107, 0.0161408, 0.0171886, 0.0274932, 0.0326007, 0.0591033};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 7.36284000000000027e-04;
	  Float_t pol1 = 8.43209000000000037e-04;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.101055, 0.0875675, 0.0764201, 0.0896778, 0.0745688, 0.0508416, 0.0388729, 0.0392451, 0.0131711, -0.00475471, 0.163426};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00350531, 0.00508545, 0.00707063, 0.00944129, 0.0095415, 0.0143669, 0.0201868, 0.0214139, 0.0346041, 0.0406055, 0.0773843};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 1.90101399999999992e-03;
	  Float_t pol1 = 1.03004500000000005e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0931397, 0.0940591, 0.0925068, 0.0748423, 0.0456136, 0.015763, 0.00281953, 0.019022, -0.0626833, 0.00124118, -0.0710237};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00505401, 0.00722202, 0.00992587, 0.0132183, 0.0131525, 0.019794, 0.0278607, 0.0295971, 0.0471327, 0.0572401, 0.103944};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 2.36398699999999969e-03;
	  Float_t pol1 = 1.28090200000000000e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0737769, 0.0729844, 0.0669981, 0.0561869, 0.054891, -0.000311302, 0.00534945, 0.0473785, -0.0374926, -0.0609766, -0.0159348};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00821781, 0.0115396, 0.0156427, 0.0204538, 0.0204816, 0.030527, 0.0428393, 0.0459644, 0.0759178, 0.0894202, 0.162667};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 3.02312299999999996e-03;
	  Float_t pol1 = 1.63804399999999996e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.092788, 0.0365816, 0.0197854, 0.107325, 0.0648891, 0.0662658, 0.0233037, -0.0236122, -0.0682863, -0.207984, -0.336516};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0171159, 0.0236563, 0.0318208, 0.0416273, 0.0413343, 0.0610424, 0.0869106, 0.093073, 0.153723, 0.187896, 0.346795};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 4.01257199999999984e-03;
	  Float_t pol1 = 2.17416299999999979e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3PionAlex6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0244232, -0.109517, -0.0827906, 0.161598, 0.0307236, -0.0999002, -0.392744, 0.145035, 0.152917, 0.0490188, 0.723248};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0426231, 0.0581957, 0.0778244, 0.100537, 0.0994573, 0.147997, 0.21388, 0.223613, 0.376758, 0.442241, 0.906303};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 5.49768000000000045e-03;
	  Float_t pol1 = 2.97884999999999990e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.0977845, 0.101319, 0.103091, 0.11351, 0.0722551, 0.10129, -0.00459366, 0.0529236, 0.173699, 0.0607847, 0.28877};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0037089, 0.00511687, 0.00740435, 0.0107892, 0.0124665, 0.0225462, 0.0361307, 0.0404895, 0.0672709, 0.0779826, 0.129084};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 1.69829837126647719e-03;
	  Float_t pol1 = 1.21632530760273721e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.110978, 0.119958, 0.134291, 0.111447, 0.114645, 0.0959741, 0.0778946, 0.0935617, 0.139158, 0.0898417, 0.035706};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00413053, 0.00561351, 0.0080384, 0.0115829, 0.0133399, 0.0239995, 0.0377082, 0.0437155, 0.0735175, 0.0854179, 0.133933};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 1.84616637865581655e-03;
	  Float_t pol1 = 1.32222872399611561e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.125172, 0.144867, 0.148934, 0.139593, 0.130398, 0.0804742, 0.0631514, 0.120612, -0.00968759, 0.0707184, -0.110791};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00353228, 0.00473819, 0.00666899, 0.00956468, 0.0108567, 0.0195176, 0.0305967, 0.0360753, 0.0642099, 0.0707663, 0.115075};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 2.13364693747489361e-03;
	  Float_t pol1 = 1.52812297971200361e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.137062, 0.133346, 0.138529, 0.143784, 0.12754, 0.101273, 0.0387196, 0.0968617, 0.110227, -0.103164, 0.252989};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00486502, 0.00644471, 0.00899533, 0.0125946, 0.0141578, 0.0251028, 0.0389365, 0.0461861, 0.0771875, 0.102432, 0.157139};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 2.60641395006084133e-03;
	  Float_t pol1 = 1.86671983155921230e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.131965, 0.139946, 0.158363, 0.0903046, 0.146713, 0.0872603, 0.113731, -0.124344, -0.000589159, 0.190572, -0.124482};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00729178, 0.00957906, 0.0131514, 0.01819, 0.0202096, 0.035221, 0.0537048, 0.0634878, 0.114044, 0.137072, 0.240824};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 3.24118098903184706e-03;
	  Float_t pol1 = 2.32134148520698243e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.138793, 0.145757, 0.0673552, 0.0754698, 0.165909, 0.14957, 0.0758562, -0.290112, -0.152356, -0.232926, 0.106716};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0122641, 0.015819, 0.0218204, 0.0294662, 0.0317047, 0.0533635, 0.0841067, 0.0986234, 0.179252, 0.226761, 0.366702};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 4.14488622049151260e-03;
	  Float_t pol1 = 2.96857730797800597e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.113443, 0.147393, 0.0654217, 0.155842, 0.0597213, -0.0416118, -0.14306, 0.273004, -0.328444, 0.192514, -0.127087};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0260424, 0.0332128, 0.0445912, 0.0593885, 0.0636087, 0.10658, 0.165977, 0.195681, 0.34504, 0.445066, 0.718467};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 5.50148450680349593e-03;
	  Float_t pol1 = 3.94017620709329313e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v3AntiprotonAlex6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 11;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 9, 11, 14, 18};
  Double_t _y[] = {0.083703, 0.083255, 0.0590565, 0.0624293, 0.0867053, 0.247706, 0.576466, 0.0501761, 0.429465, 0.326308, 1.31989};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0663444, 0.0834688, 0.109026, 0.146809, 0.155555, 0.258989, 0.405504, 0.462926, 0.780885, 0.994817, 2.2086};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t pol0 = 7.53765601848897532e-03;
	  Float_t pol1 = 5.39848705646177145e-03;

	  Float_t nonflow = (pol0 + pol1*_x[i]);
	  Float_t systerr = TMath::Sqrt( nonflow*nonflow + (0.01*0.01 + 0.015*0.015 + 0.01*0.01 + 0.04*0.04 + 0.04*0.04 + 0.03*0.03 + 0.05*0.05)*_y[i]*_y[i]);
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// QM (SP eta gap > 1) by Mikolaj
TGraphErrors* v22_etagap10_5060_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, -0.0873219, 0.0597375, 0.0743915, 0.0931236, 0.109237, 0.122847, 0.135165, 0.150852, 0.160262, 0.166853, 0.178234, 0.192645, 0.192824, 0.196815, 0.200117, 0.202778, 0.1995, 0.210411, 0.214184, 0.223324, 0.217727, 0.213249, 0.208803, 0.209717, 0.178609, 0.226572, 0.230321, 0.243208, 0.177178, 0.197599, 0.263397, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0.204799, 0.00488736, 0.0036102, 0.00391356, 0.00438879, 0.00490408, 0.00543342, 0.00607888, 0.00663274, 0.00717516, 0.0078937, 0.00876489, 0.00942504, 0.0102981, 0.0113314, 0.0125538, 0.0139245, 0.0158814, 0.0182328, 0.0208889, 0.0240269, 0.0277764, 0.0317578, 0.0359011, 0.0412365, 0.0445151, 0.0492254, 0.0546197, 0.0355352, 0.0709352, 0.465996, 0.110652, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0.0999384, 0.00245063, 0.00153724, 0.00147518, 0.00156441, 0.00172669, 0.00192613, 0.00215318, 0.00241725, 0.00272092, 0.00306104, 0.00347236, 0.00392824, 0.0044454, 0.00506733, 0.00579985, 0.00667986, 0.0077041, 0.00901047, 0.0104883, 0.012436, 0.0144056, 0.0167576, 0.0190382, 0.0220889, 0.0237475, 0.0258547, 0.0288486, 0.0189844, 0.0390867, 0.250826, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

 TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_4050_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, -0.000724975, 0.0552354, 0.0732448, 0.0906073, 0.107367, 0.120443, 0.134765, 0.148343, 0.156962, 0.171834, 0.182009, 0.186255, 0.197146, 0.201325, 0.207395, 0.212444, 0.223133, 0.227715, 0.215306, 0.227535, 0.229997, 0.230652, 0.246553, 0.213974, 0.225505, 0.228866, 0.20663, 0.213873, 0.211377, 0.159562, 0.269798, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0.12403, 0.00318662, 0.0027737, 0.0031637, 0.00364762, 0.00407577, 0.00455525, 0.00503037, 0.00539388, 0.00594982, 0.00639986, 0.00674716, 0.00729324, 0.00774568, 0.00832146, 0.00898665, 0.00986877, 0.0108603, 0.0118255, 0.013404, 0.0151841, 0.017233, 0.0197056, 0.0218697, 0.0247173, 0.0273353, 0.0300092, 0.0323873, 0.0218155, 0.0442529, 0.253422, 0.0522819, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0.0661206, 0.00149826, 0.000936119, 0.000890212, 0.000936912, 0.00102774, 0.00114045, 0.00126956, 0.0014225, 0.00160352, 0.00179964, 0.002034, 0.00229325, 0.002606, 0.00297516, 0.0034027, 0.00390643, 0.00455883, 0.00535362, 0.00623856, 0.00738427, 0.00865204, 0.0100195, 0.0113015, 0.0129887, 0.014458, 0.0159583, 0.0168897, 0.0114366, 0.0240444, 0.139903, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_3040_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0.0134225, 0.0519218, 0.0665359, 0.0842276, 0.100815, 0.114049, 0.127699, 0.140248, 0.150695, 0.16182, 0.172197, 0.18134, 0.187407, 0.192411, 0.201391, 0.211336, 0.209113, 0.217154, 0.220327, 0.219292, 0.213611, 0.223944, 0.231156, 0.222101, 0.210007, 0.230393, 0.229598, 0.231779, 0.221272, 0.216799, 0.303872, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0.0872183, 0.00247355, 0.00232267, 0.00276805, 0.00324864, 0.00365973, 0.00409325, 0.00450541, 0.00486738, 0.00526119, 0.00564843, 0.00602156, 0.00633633, 0.00665612, 0.00712018, 0.0076721, 0.00800948, 0.00870467, 0.00945645, 0.0103166, 0.0113037, 0.012766, 0.0144866, 0.0160066, 0.0177591, 0.019504, 0.0214444, 0.0231496, 0.0161895, 0.0316968, 0.202708, 0.0396501, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0.0464175, 0.00105187, 0.000653968, 0.000618086, 0.000645372, 0.000704275, 0.000779312, 0.00086849, 0.000970154, 0.00108892, 0.00122422, 0.00138357, 0.00156242, 0.00177261, 0.00201701, 0.00231225, 0.00266135, 0.0030997, 0.00363381, 0.00428479, 0.00504451, 0.00586784, 0.00689493, 0.00788004, 0.00900252, 0.00993424, 0.010906, 0.0118421, 0.00803899, 0.0172747, 0.114742, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_2030_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0.0675354, 0.0453406, 0.0597064, 0.0729257, 0.0856194, 0.099743, 0.11087, 0.12233, 0.133139, 0.142766, 0.150798, 0.15726, 0.164193, 0.17088, 0.180412, 0.185874, 0.187267, 0.190105, 0.19678, 0.193095, 0.203633, 0.189152, 0.197994, 0.206382, 0.212811, 0.19176, 0.195015, 0.197307, 0.202853, 0.1925, 0.252937, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0.065632, 0.00201808, 0.00201334, 0.00235429, 0.00272524, 0.00315311, 0.0035035, 0.0038686, 0.00422263, 0.00455118, 0.00484569, 0.00510894, 0.00540221, 0.00571347, 0.00612565, 0.00647902, 0.00677918, 0.00722552, 0.00786669, 0.00838991, 0.00940688, 0.0101343, 0.0114481, 0.0127734, 0.014505, 0.0155823, 0.0168149, 0.0184895, 0.0131523, 0.0251873, 0.161027, 0.0348376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0.0344332, 0.000816139, 0.0005061, 0.000476092, 0.00049599, 0.000538675, 0.000594465, 0.000658739, 0.000735579, 0.000824272, 0.000930675, 0.00104848, 0.00118504, 0.00134617, 0.00153402, 0.00176201, 0.00202544, 0.00237128, 0.00279278, 0.00327205, 0.0038541, 0.00454915, 0.0053078, 0.00607066, 0.00707512, 0.00785901, 0.00858144, 0.0094215, 0.00634443, 0.0135698, 0.0903221, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_1020_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, -0.0244474, 0.0361825, 0.0458614, 0.0568154, 0.0653387, 0.0761904, 0.0855912, 0.0931017, 0.100605, 0.108806, 0.115534, 0.122009, 0.12765, 0.133345, 0.139141, 0.140422, 0.142585, 0.149623, 0.153545, 0.157209, 0.158014, 0.159055, 0.161039, 0.163363, 0.160728, 0.173582, 0.151455, 0.159009, 0.162018, 0.148968, 0.552532, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0.0549447, 0.00166051, 0.00157917, 0.00185472, 0.00210364, 0.00243219, 0.00272718, 0.0029743, 0.00322744, 0.00350813, 0.00375677, 0.00400998, 0.00425557, 0.00452476, 0.00482201, 0.00504417, 0.0053393, 0.00583931, 0.00636429, 0.00699075, 0.00768206, 0.00858417, 0.00960373, 0.0107919, 0.0120125, 0.0133585, 0.0142786, 0.0157435, 0.0111965, 0.0217748, 0.107326, 0.0304973, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0.0297487, 0.000693128, 0.000430282, 0.000403147, 0.000418835, 0.000453801, 0.00049934, 0.000553943, 0.000617715, 0.000693185, 0.00078017, 0.000880407, 0.000998026, 0.0011363, 0.00129749, 0.00149084, 0.00171888, 0.00200843, 0.00237001, 0.00278848, 0.0032784, 0.00386639, 0.00452296, 0.005239, 0.00598587, 0.00667402, 0.00732162, 0.00816298, 0.00547164, 0.0118032, 0.0561392, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_510_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v22_etagap10_510_pion_nPoints = 32;
  if (last>v22_etagap10_510_pion_nPoints-1) last=v22_etagap10_510_pion_nPoints-1;
  if (last<0 && first<0) last=v22_etagap10_510_pion_nPoints-1;
  if (last<0) last=v22_etagap10_510_pion_nPoints-1+last;
  if (first<0) first=0;
  Double_t v22_etagap10_510_pion_x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6};
  Double_t v22_etagap10_510_pion_y[] = {0.,0.,0.,0.026881, 0.0316749, 0.0381939, 0.0462715, 0.0527852, 0.0589098, 0.0674904, 0.0711453, 0.0749032, 0.078903, 0.0833461, 0.0894323, 0.0928401, 0.097527, 0.0963054, 0.103171, 0.108102, 0.104966, 0.104179, 0.111621, 0.110965, 0.118524, 0.119141, 0.123375, 0.118692, 0.10753, 0.112698, 0.103193, 0.0966573};
  Double_t v22_etagap10_510_pion_xerr[] = {0.,0.,0.,0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v22_etagap10_510_pion_yerr[] = {0.,0.,0.,0.00214347, 0.00155063, 0.00162521, 0.00183525, 0.00205243, 0.00228006, 0.00258007, 0.00278512, 0.00301629, 0.00327933, 0.00358512, 0.00395689, 0.0043373, 0.00479985, 0.00524212, 0.00591262, 0.00673476, 0.00761607, 0.00876874, 0.0102147, 0.0118854, 0.013785, 0.0157642, 0.0180509, 0.0199667, 0.0219864, 0.0245517, 0.016479, 0.0353244};
  TGraphErrors* graph = new TGraphErrors(last-first+1, &v22_etagap10_510_pion_x[first], &v22_etagap10_510_pion_y[first], &v22_etagap10_510_pion_xerr[first], &v22_etagap10_510_pion_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_5060_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.0592854, 0.0581809, 0.0528559, 0.0758825, 0.0948136, 0.103104, 0.118047, 0.139711, 0.14803, 0.14964, 0.172821, 0.173198, 0.191228, 0.191121, 0.206676, 0.188007, 0.180187, 0.200772, 0.21883, 0.179475, 0.187122, 0.218256, 0.250577, 0.223592, 0.292999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.0327399, 0.00952958, 0.00630952, 0.00536493, 0.0051041, 0.00505363, 0.00518118, 0.00541767, 0.00568562, 0.0060391, 0.00652267, 0.00706427, 0.00770117, 0.00846333, 0.00957273, 0.0109742, 0.0128912, 0.0151542, 0.0179697, 0.0223891, 0.0282206, 0.0392386, 0.0554708, 0.10338, 0.262149, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0623433, 0.0177094, 0.0118733, 0.0103362, 0.0100447, 0.0100381, 0.0104779, 0.0111214, 0.0116894, 0.0123778, 0.0135523, 0.0144775, 0.0158768, 0.017185, 0.0192743, 0.0216116, 0.0248209, 0.0293869, 0.0347478, 0.0427324, 0.0536825, 0.0742418, 0.104119, 0.180642, 0.555504, 0.0456015, 0.0717251, 0.202558, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_4050_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.0619661, 0.0400546, 0.0619503, 0.0785653, 0.0916647, 0.105258, 0.118573, 0.140119, 0.148773, 0.157981, 0.167909, 0.176368, 0.190243, 0.197994, 0.192225, 0.201264, 0.208214, 0.217569, 0.215782, 0.209895, 0.205265, 0.214319, 0.231061, 0.257986, 0.383446, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.0201743, 0.00584789, 0.00382593, 0.00323562, 0.00303369, 0.00301693, 0.00305804, 0.00316891, 0.00329939, 0.00352794, 0.00377958, 0.00406431, 0.00447333, 0.00491623, 0.00554596, 0.00633044, 0.00740304, 0.0086718, 0.0104701, 0.0128577, 0.0165705, 0.0233379, 0.0316965, 0.0604997, 0.142815, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0368872, 0.010646, 0.00721477, 0.00638506, 0.00625038, 0.00641893, 0.00670512, 0.00724411, 0.00762642, 0.00811159, 0.00868884, 0.00928689, 0.0101344, 0.0109697, 0.0118869, 0.0132718, 0.015169, 0.0173855, 0.0205779, 0.025083, 0.031145, 0.0430324, 0.0585126, 0.102606, 0.255642, 0.0254592, 0.0349152, 0.0887712, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_3040_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.000553789, 0.0404834, 0.0500686, 0.0659864, 0.0809581, 0.0979551, 0.111338, 0.125787, 0.138353, 0.152205, 0.162979, 0.171183, 0.178498, 0.189554, 0.194077, 0.195814, 0.204214, 0.210309, 0.215834, 0.206115, 0.227138, 0.20424, 0.216186, 0.203484, 0.1687, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.0142928, 0.00413782, 0.00269443, 0.00225036, 0.00211002, 0.00207131, 0.00208791, 0.00214672, 0.00224472, 0.00238231, 0.00255078, 0.00274784, 0.00298966, 0.00331916, 0.00372056, 0.00426177, 0.00498107, 0.00590833, 0.00707412, 0.0086493, 0.0114182, 0.0157553, 0.0219151, 0.0404066, 0.143742, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0259527, 0.00764997, 0.00519162, 0.00462467, 0.00462972, 0.0048717, 0.00515435, 0.00553526, 0.00593861, 0.00641789, 0.0068783, 0.00731213, 0.00780885, 0.00848646, 0.00915576, 0.00997916, 0.0112493, 0.0128153, 0.0148045, 0.017474, 0.022579, 0.0298356, 0.0427476, 0.0757909, 0.238668, 0.0178573, 0.0263078, 0.0551329, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_2030_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.0212288, 0.0259616, 0.0411856, 0.0524764, 0.0668708, 0.0808397, 0.0942864, 0.106513, 0.119792, 0.130224, 0.139591, 0.15116, 0.15632, 0.16214, 0.1692, 0.170107, 0.183808, 0.194388, 0.191664, 0.192513, 0.192104, 0.197077, 0.200502, 0.182628, 0.147351, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.0107508, 0.00320911, 0.00207722, 0.00172773, 0.00161598, 0.00158404, 0.00158521, 0.00162636, 0.00169633, 0.00179334, 0.00191295, 0.00206472, 0.00224106, 0.00247143, 0.00279206, 0.00319885, 0.00375698, 0.00443323, 0.00534009, 0.00661798, 0.00862949, 0.0118678, 0.0174641, 0.030209, 0.100816, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0197657, 0.00587757, 0.00400434, 0.00356801, 0.0036119, 0.00382296, 0.00410052, 0.00441878, 0.00480124, 0.00516269, 0.00552864, 0.00597777, 0.00631951, 0.00674724, 0.00731513, 0.00790242, 0.00897728, 0.0101926, 0.0115865, 0.0137535, 0.0170439, 0.0226759, 0.0325735, 0.0548696, 0.18938, 0.0135812, 0.0203971, 0.0475816, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_1020_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.026086, 0.0200684, 0.0287286, 0.0413946, 0.0518688, 0.0605028, 0.0735013, 0.0829178, 0.0920512, 0.0970675, 0.109182, 0.112781, 0.121785, 0.124268, 0.131138, 0.140937, 0.137646, 0.148586, 0.147178, 0.159076, 0.153502, 0.152984, 0.164356, 0.185492, 0.199013, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.00909186, 0.00275495, 0.00176528, 0.00146439, 0.00136497, 0.00132616, 0.00132809, 0.00135637, 0.00141312, 0.00149465, 0.00159557, 0.00171572, 0.00187653, 0.00206897, 0.00234094, 0.00268047, 0.00314097, 0.00371948, 0.00447599, 0.00550291, 0.00716883, 0.0100772, 0.0149141, 0.0263166, 0.0846981, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0163953, 0.00499284, 0.00332771, 0.00295984, 0.00296082, 0.00305623, 0.00330957, 0.00354906, 0.00382349, 0.00403788, 0.00443445, 0.0046715, 0.00507456, 0.00539324, 0.00589863, 0.00656852, 0.00718652, 0.00826722, 0.00947331, 0.0113341, 0.0140336, 0.0190078, 0.0276842, 0.0481405, 0.142308, 0.0116544, 0.0178366, 0.0362519, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_510_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v22_etagap10_510_kaon_nPoints = 26;
  if (last>v22_etagap10_510_kaon_nPoints-1) last=v22_etagap10_510_kaon_nPoints-1;
  if (last<0 && first<0) last=v22_etagap10_510_kaon_nPoints-1;
  if (last<0) last=v22_etagap10_510_kaon_nPoints-1+last;
  if (first<0) first=0;
  Double_t v22_etagap10_510_kaon_x[] = {0.05, 0.15, 0.25,0.35,0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55};
  Double_t v22_etagap10_510_kaon_y[] = {0.,0.,0.,0.,0.00966259, 0.0207807, 0.0256179, 0.0355533, 0.042711, 0.0502604, 0.0580054, 0.0591413, 0.0659689, 0.0769689, 0.076443, 0.0875564, 0.0921437, 0.0983687, 0.0969186, 0.0992941, 0.101035, 0.102316, 0.102327, 0.0984255, 0.119298, 0.138143};
  Double_t v22_etagap10_510_kaon_xerr[] = {0.,0.,0.,0.,0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v22_etagap10_510_kaon_yerr[] = {0.,0.,0.,0.,0.00784647, 0.00508557, 0.00427572, 0.00407675, 0.00403995, 0.00413015, 0.00430378, 0.0044739, 0.00477018, 0.00519501, 0.00552012, 0.00606902, 0.00663695, 0.00745315, 0.00835876, 0.009673, 0.0113138, 0.0136154, 0.0165971, 0.0213286, 0.0291213, 0.0427121};
  TGraphErrors* graph = new TGraphErrors(last-first+1, &v22_etagap10_510_kaon_x[first], &v22_etagap10_510_kaon_y[first], &v22_etagap10_510_kaon_xerr[first], &v22_etagap10_510_kaon_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_5060_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 1.04828, 0.0216201, 0.0288792, 0.0475747, 0.0452145, 0.0510448, 0.0916252, 0.0800232, 0.0999849, 0.133603, 0.151354, 0.166441, 0.164859, 0.189033, 0.191063, 0.224191, 0.2235, 0.208137, 0.262829, 0.234087, 0.250857, 0.259322, 0.221027, 0.287725, 0.297794, 0.245738, 0.296381, 0.228615, 0.265392, 0.1527, 0.702683, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.00267392, 0.0480364, 0.0204747, 0.0141307, 0.0116684, 0.0105368, 0.0103799, 0.0102472, 0.0103013, 0.0106565, 0.0109217, 0.0113955, 0.0122391, 0.0130204, 0.0138631, 0.0148995, 0.0163369, 0.0172412, 0.0191165, 0.0208819, 0.0221608, 0.0248779, 0.0268306, 0.0292929, 0.0327247, 0.0364108, 0.0390226, 0.0242139, 0.0370034, 0.0681912, 0.36359, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0318153, 0.0863738, 0.0349737, 0.024328, 0.0204729, 0.0186663, 0.0184856, 0.0183288, 0.0185378, 0.0193167, 0.0201325, 0.0210981, 0.0223835, 0.0240191, 0.0257134, 0.0278661, 0.0303526, 0.0317838, 0.0355249, 0.0383376, 0.0412807, 0.0459884, 0.0500093, 0.0540488, 0.0595738, 0.0658101, 0.0698881, 0.0445687, 0.0687029, 0.124069, 0.596873, 0.0761726, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_4050_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0, -0.00253575, 0.0196403, 0.0296492, 0.0558592, 0.0585807, 0.0789371, 0.0829314, 0.106347, 0.109234, 0.130626, 0.135786, 0.154787, 0.183858, 0.196135, 0.210303, 0.204343, 0.231183, 0.231211, 0.249107, 0.221292, 0.273605, 0.266969, 0.257376, 0.27806, 0.280984, 0.292262, 0.26905, 0.274387, 0.211313, 0.324339, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0.0313968, 0.013031, 0.00897124, 0.00746666, 0.00669911, 0.00643486, 0.00625644, 0.00625102, 0.0063233, 0.00648134, 0.00679227, 0.00711397, 0.00756203, 0.00803456, 0.00861781, 0.00930825, 0.00993977, 0.0107897, 0.0119757, 0.012844, 0.014168, 0.0151713, 0.0166566, 0.0184211, 0.0202448, 0.0218659, 0.0140385, 0.0214924, 0.0430906, 0.178145, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.137912, 0.0550066, 0.0223057, 0.0155486, 0.0131162, 0.0119509, 0.0115947, 0.0113617, 0.0116088, 0.011777, 0.0123083, 0.0128657, 0.0135363, 0.0147092, 0.0156493, 0.0168069, 0.0178454, 0.0193603, 0.0208125, 0.0228027, 0.0242045, 0.0268854, 0.029038, 0.031391, 0.0345494, 0.0377226, 0.040746, 0.027034, 0.040086, 0.0771034, 0.360776, 0.0543393, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_3040_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.338705, 0.0518824, 0.0390561, 0.0276461, 0.0288275, 0.0450128, 0.0541352, 0.0696097, 0.078628, 0.0945815, 0.114324, 0.121956, 0.143111, 0.16049, 0.163625, 0.187258, 0.196448, 0.198742, 0.20469, 0.213722, 0.243047, 0.24687, 0.247373, 0.248214, 0.278438, 0.27313, 0.276065, 0.28618, 0.275928, 0.255872, 0.115524, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.000645641, 0.0219796, 0.00970872, 0.00665225, 0.00546807, 0.0048868, 0.00460971, 0.00445534, 0.00443915, 0.00440264, 0.00449212, 0.00469103, 0.00489613, 0.00509739, 0.00537234, 0.00577706, 0.00619733, 0.00656813, 0.00723524, 0.00786358, 0.00854128, 0.00919734, 0.0100571, 0.0107599, 0.0121505, 0.0133044, 0.0145959, 0.00927871, 0.0144753, 0.0271428, 0.101974, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.0103031, 0.0399007, 0.0167902, 0.0116564, 0.00967558, 0.0087891, 0.00842169, 0.00826439, 0.00832211, 0.00845603, 0.00884102, 0.00924786, 0.00986379, 0.010438, 0.0109558, 0.0119274, 0.0127354, 0.0135081, 0.0146049, 0.0157881, 0.0172627, 0.0184926, 0.020007, 0.0214871, 0.0239718, 0.0257944, 0.0281455, 0.0192021, 0.0280399, 0.0496494, 0.182599, 0.0387317, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_2030_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.662474, 0.0185994, 0.0117562, 0.0259402, 0.0214211, 0.0304214, 0.0394749, 0.0546118, 0.0597476, 0.0794185, 0.0866745, 0.0970393, 0.10449, 0.129477, 0.136203, 0.149914, 0.158318, 0.181555, 0.188754, 0.194446, 0.197162, 0.223073, 0.227114, 0.23305, 0.245297, 0.240345, 0.245049, 0.259105, 0.26242, 0.287932, 0.243008, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.364261, 0.0177086, 0.00772697, 0.00537977, 0.00442251, 0.00393434, 0.0036613, 0.00350753, 0.00343486, 0.00345618, 0.00351172, 0.00356278, 0.00370219, 0.0038888, 0.00407012, 0.00430605, 0.00459197, 0.00497422, 0.0053417, 0.00573583, 0.00627158, 0.00680288, 0.00749132, 0.00815469, 0.00883801, 0.00985937, 0.0107836, 0.00687874, 0.0105095, 0.020312, 0.0866113, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 0.786284, 0.0311215, 0.0133156, 0.00936461, 0.00774258, 0.00698169, 0.00658302, 0.00644905, 0.00637499, 0.00661681, 0.00679006, 0.00704375, 0.00736515, 0.00801826, 0.00841806, 0.00901493, 0.00960231, 0.0105292, 0.0112683, 0.0120311, 0.0128863, 0.0141925, 0.0152828, 0.0165364, 0.0178314, 0.0193874, 0.02123, 0.0149339, 0.0209609, 0.0385035, 0.152991, 0.0307381, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_1020_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 48;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 9.9};
  Double_t _y[] = {0, 0, 0, 0.0255935, 0.0326516, 0.00631974, 0.0156774, 0.0138975, 0.0183874, 0.0223836, 0.0376923, 0.0473871, 0.0504273, 0.0598239, 0.0756639, 0.0858968, 0.094869, 0.105419, 0.109821, 0.124009, 0.132104, 0.135292, 0.14661, 0.154133, 0.151872, 0.16556, 0.172129, 0.1835, 0.193187, 0.206297, 0.192187, 0.209371, 0.236101, 0.215908, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0.823453, 0.0150733, 0.00662216, 0.0047245, 0.00385231, 0.00342087, 0.00319976, 0.00302894, 0.00295194, 0.00293546, 0.00294981, 0.00303874, 0.00312892, 0.00323678, 0.00340637, 0.00359611, 0.00385576, 0.00412362, 0.00442889, 0.00484538, 0.00520886, 0.00565058, 0.00614164, 0.00677307, 0.00739301, 0.00811282, 0.00894112, 0.0057491, 0.00900794, 0.0178728, 0.0777444, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerrSyst[] = {0, 0, 0, 1.28647, 0.0263428, 0.0113809, 0.00814393, 0.00674068, 0.00603507, 0.00567695, 0.00548342, 0.00543364, 0.0054437, 0.00558493, 0.00589599, 0.00618781, 0.00649389, 0.00690933, 0.00729146, 0.00789278, 0.00847596, 0.0090342, 0.00983544, 0.0105641, 0.0112837, 0.0122616, 0.0133408, 0.0146642, 0.0159207, 0.017622, 0.0120255, 0.017693, 0.0329639, 0.141832, 0.026807, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v22_etagap10_510_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v22_etagap10_510_antiproton_nPoints = 33;
  if (last>v22_etagap10_510_antiproton_nPoints-1) last=v22_etagap10_510_antiproton_nPoints-1;
  if (last<0 && first<0) last=v22_etagap10_510_antiproton_nPoints-1;
  if (last<0) last=v22_etagap10_510_antiproton_nPoints-1+last;
  if (first<0) first=0;
  Double_t v22_etagap10_510_antiproton_x[] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.2, 3.6, 4};
  Double_t v22_etagap10_510_antiproton_y[] = {0.,0.,0.,0.,0.,0.0107737, 0.0128734, 0.00917658, 0.0105757, 0.0115926, 0.0329602, 0.0261356, 0.0313384, 0.0452816, 0.0526354, 0.0618944, 0.0574567, 0.0701679, 0.0701348, 0.0875208, 0.0878755, 0.104836, 0.10294, 0.110327, 0.11429, 0.106938, 0.115266, 0.107332, 0.133372, 0.130432, 0.121658, 0.154643, 0.105422};
  Double_t v22_etagap10_510_antiproton_xerr[] = {0.,0.,0.,0.,0.,0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v22_etagap10_510_antiproton_yerr[] = {0.,0.,0.,0.,0.,0.0187085, 0.0131334, 0.0109682, 0.00970316, 0.00905275, 0.00863548, 0.00845393, 0.00840426, 0.00854155, 0.00874073, 0.00902577, 0.00933217, 0.00988566, 0.0104503, 0.011238, 0.0120288, 0.0130014, 0.0139594, 0.0152531, 0.0166315, 0.0177388, 0.0196984, 0.0216449, 0.0237533, 0.025945, 0.016823, 0.0267045, 0.0516681};
  TGraphErrors* graph = new TGraphErrors(last-first+1, &v22_etagap10_510_antiproton_x[first], &v22_etagap10_510_antiproton_y[first], &v22_etagap10_510_antiproton_xerr[first], &v22_etagap10_510_antiproton_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// QM11 Ks and lambda
TGraphErrors* v2KsQM114050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 22;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.2, 3.6, 4, 4.4, 4.8};
  Double_t _y[] = {0.0149204, 0.0586289, 0.0649511, 0.0869282, 0.095823, 0.118562, 0.129187, 0.143078, 0.151888, 0.144034, 0.165162, 0.201767, 0.179882, 0.20807, 0.193481, 0.201477, 0.190785, 0.218286, 0.201261, 0.195918, 0.200895, 0.162648, 0.202368, 0.19096, 0.143372};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0300707, 0.021654, 0.0182861, 0.0190635, 0.0186501, 0.0204063, 0.0204742, 0.02106, 0.0211224, 0.0196578, 0.021203, 0.0238578, 0.0215306, 0.0235206, 0.0196785, 0.0196917, 0.0200961, 0.0221347, 0.0230998, 0.0261922, 0.0234749, 0.0306236, 0.0433844, 0.0608945, 0.0726564};
  Double_t _yerrS[] = {0.0300707, 0.021654, 0.0182861, 0.0190635, 0.0186501, 0.0204063, 0.0204742, 0.02106, 0.0211224, 0.0196578, 0.021203, 0.0238578, 0.0215306, 0.0235206, 0.0196785, 0.0196917, 0.0200961, 0.0221347, 0.0230998, 0.0261922, 0.0234749, 0.0306236, 0.0433844, 0.0608945, 0.0726564};

  _yerr[1-1]=0.0298;_yerrS[1-1]=0.0043;
  _yerr[2-1]=0.0170;_yerrS[2-1]=0.0134;
  _yerr[3-1]=0.0129;_yerrS[3-1]=0.0129;
  _yerr[4-1]=0.0111;_yerrS[4-1]=0.0155;
  _yerr[5-1]=0.0103;_yerrS[5-1]=0.0156;
  _yerr[6-1]=0.0101;_yerrS[6-1]=0.0177;
  _yerr[7-1]=0.0099;_yerrS[7-1]=0.0179;
  _yerr[8-1]=0.0100;_yerrS[8-1]=0.0185;
  _yerr[9-1]=0.0104;_yerrS[9-1]=0.0184;
  _yerr[10-1]=0.0109;_yerrS[10-1]=0.0164;
  _yerr[11-1]=0.0117;_yerrS[11-1]=0.0177;
  _yerr[12-1]=0.0125;_yerrS[12-1]=0.0203;
  _yerr[13-1]=0.0131;_yerrS[13-1]=0.0171;
  _yerr[14-1]=0.0143;_yerrS[14-1]=0.0187;
  _yerr[15-1]=0.0115;_yerrS[15-1]=0.0159;
  _yerr[16-1]=0.0129;_yerrS[16-1]=0.0148;
  _yerr[17-1]=0.0157;_yerrS[17-1]=0.0126;
  _yerr[18-1]=0.0181;_yerrS[18-1]=0.0128;
  _yerr[19-1]=0.0206;_yerrS[19-1]=0.0105;
  _yerr[20-1]=0.0246;_yerrS[20-1]=0.0090;
  _yerr[21-1]=0.0222;_yerrS[21-1]=0.0076;
  _yerr[22-1]=0.0303;_yerrS[22-1]=0.0046;
  _yerr[23-1]=0.0432;_yerrS[23-1]=0.0039;
  _yerr[24-1]=0.0609;_yerrS[24-1]=0.0022;
  _yerr[25-1]=0.0727;_yerrS[25-1]=0.0006;

 if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2LambdaQM114050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 10;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.1, 1.3, 1.5, 1.7, 1.9, 2.15, 2.45, 2.75, 2.95, 3.5, 4.5};
  Double_t _y[] = {0.070888, 0.119342, 0.149007, 0.174366, 0.179891, 0.219007, 0.246098, 0.280332, 0.233308, 0.254311, 0.30027};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0193462, 0.026561, 0.0322, 0.036856, 0.0379916, 0.0464809, 0.0503634, 0.0963763, 0.0946555, 0.0545124, 0.0964832};
  Double_t _yerrS[] = {0.0193462, 0.026561, 0.0322, 0.036856, 0.0379916, 0.0464809, 0.0503634, 0.0963763, 0.0946555, 0.0545124, 0.0964832};

  _yerr[1-1]=0.0134;_yerrS[1-1]=0.0139;
  _yerr[2-1]=0.0119;_yerrS[2-1]=0.0237;
  _yerr[3-1]=0.0116;_yerrS[3-1]=0.0300;
  _yerr[4-1]=0.0121;_yerrS[4-1]=0.0348;
  _yerr[5-1]=0.0133;_yerrS[5-1]=0.0356;
  _yerr[6-1]=0.0150;_yerrS[6-1]=0.0440;
  _yerr[7-1]=0.0139;_yerrS[7-1]=0.0484;
  _yerr[8-1]=0.0882;_yerrS[8-1]=0.0389;
  _yerr[9-1]=0.0740;_yerrS[9-1]=0.0590;
  _yerr[10-1]=0.0213;_yerrS[10-1]=0.0502;
  _yerr[11-1]=0.0826;_yerrS[11-1]=0.0498;


 if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2KsQM111020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 22;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.2, 3.6, 4, 4.4, 4.8};
  Double_t _y[] = {0.037872, 0.0427513, 0.040969, 0.042276, 0.0585498, 0.0594143, 0.0734317, 0.0825041, 0.0995797, 0.104865, 0.103593, 0.114106, 0.121967, 0.134521, 0.133603, 0.143376, 0.153956, 0.168784, 0.155898, 0.155901, 0.165062, 0.165088, 0.108553, 0.132339, 0.0308495};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0229751, 0.0146732, 0.0115235, 0.01008, 0.0113093, 0.0106629, 0.0116316, 0.0120922, 0.0133899, 0.0135076, 0.0127666, 0.0133961, 0.0137166, 0.0143667, 0.0126623, 0.0129384, 0.0137158, 0.0146097, 0.0157186, 0.0171808, 0.0151266, 0.0197935, 0.0268463, 0.0391686, 0.0471636};
  Double_t _yerrS[] = {0.0229751, 0.0146732, 0.0115235, 0.01008, 0.0113093, 0.0106629, 0.0116316, 0.0120922, 0.0133899, 0.0135076, 0.0127666, 0.0133961, 0.0137166, 0.0143667, 0.0126623, 0.0129384, 0.0137158, 0.0146097, 0.0157186, 0.0171808, 0.0151266, 0.0197935, 0.0268463, 0.0391686, 0.0471636};

  _yerr[1-1]=0.0202;_yerrS[1-1]=0.0109;
  _yerr[2-1]=0.0109;_yerrS[2-1]=0.0098;
  _yerr[3-1]=0.0082;_yerrS[3-1]=0.0081;
  _yerr[4-1]=0.0067;_yerrS[4-1]=0.0075;
  _yerr[5-1]=0.0061;_yerrS[5-1]=0.0095;
  _yerr[6-1]=0.0059;_yerrS[6-1]=0.0089;
  _yerr[7-1]=0.0056;_yerrS[7-1]=0.0102;
  _yerr[8-1]=0.0057;_yerrS[8-1]=0.0107;
  _yerr[9-1]=0.0058;_yerrS[9-1]=0.0121;
  _yerr[10-1]=0.0063;_yerrS[10-1]=0.0119;
  _yerr[11-1]=0.0063;_yerrS[11-1]=0.0111;
  _yerr[12-1]=0.0069;_yerrS[12-1]=0.0115;
  _yerr[13-1]=0.0073;_yerrS[13-1]=0.0116;
  _yerr[14-1]=0.0078;_yerrS[14-1]=0.0121;
  _yerr[15-1]=0.0063;_yerrS[15-1]=0.0110;
  _yerr[16-1]=0.0075;_yerrS[16-1]=0.0106;
  _yerr[17-1]=0.0092;_yerrS[17-1]=0.0101;
  _yerr[18-1]=0.0107;_yerrS[18-1]=0.0099;
  _yerr[19-1]=0.0135;_yerrS[19-1]=0.0081;
  _yerr[20-1]=0.0156;_yerrS[20-1]=0.0072;
  _yerr[21-1]=0.0138;_yerrS[21-1]=0.0062;
  _yerr[22-1]=0.0192;_yerrS[22-1]=0.0046;
  _yerr[23-1]=0.0268;_yerrS[23-1]=0.0021;
  _yerr[24-1]=0.0391;_yerrS[24-1]=0.0015;
  _yerr[25-1]=0.0472;_yerrS[25-1]=0.0001;


 if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2LambdaQM111020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 10;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.1, 1.3, 1.5, 1.7, 1.9, 2.15, 2.45, 2.75, 2.95, 3.5, 4.5};
  Double_t _y[] = {0.0329731, 0.0598229, 0.0582877, 0.0862622, 0.105079, 0.11425, 0.131433, 0.165193, 0.167494, 0.206327, 0.156109};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0128034, 0.0148904, 0.0134535, 0.0186384, 0.0223601, 0.0240979, 0.0275226, 0.0479114, 0.041216, 0.0458347, 0.0765811};
  Double_t _yerrS[] = {0.0128034, 0.0148904, 0.0134535, 0.0186384, 0.0223601, 0.0240979, 0.0275226, 0.0479114, 0.041216, 0.0458347, 0.0765811};

  _yerr[1-1]=0.0111;_yerrS[1-1]=0.0065;
  _yerr[2-1]=0.0082;_yerrS[2-1]=0.0124;
  _yerr[3-1]=0.0073;_yerrS[3-1]=0.0113;
  _yerr[4-1]=0.0070;_yerrS[4-1]=0.0173;
  _yerr[5-1]=0.0074;_yerrS[5-1]=0.0211;
  _yerr[6-1]=0.0080;_yerrS[6-1]=0.0227;
  _yerr[7-1]=0.0072;_yerrS[7-1]=0.0266;
  _yerr[8-1]=0.0408;_yerrS[8-1]=0.0252;
  _yerr[9-1]=0.0282;_yerrS[9-1]=0.0300;
  _yerr[10-1]=0.0283;_yerrS[10-1]=0.0361;
  _yerr[11-1]=0.0637;_yerrS[11-1]=0.0425;


 if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}


// Kink with VZERO SP method
TGraphErrors* v2Kink0005(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.00278912, 0.00150957, 0.00969046, 0.0104301, 0.0188253, 0.0227717, 0.0309487, 0.0341199, 0.0505573, 0.0482094, 0.0587934, 0.0527458, 0.058417, 0.0243267, -0.00413263, 0.049706, 0.124202};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00265733, 0.0021355, 0.00200978, 0.00210868, 0.0023306, 0.00266326, 0.00233616, 0.00321753, 0.00440557, 0.00612543, 0.00843734, 0.0115667, 0.0125789, 0.0219909, 0.0384996, 0.0573552, 0.0832343};
  Double_t _yerrSyst[] = {0.000115953, 0.000744338, 0.000801153, 0.001446, 0.00174913, 0.00237721, 0.0026208, 0.00388338, 0.00370303, 0.00451601, 0.00405148, 0.0044871, 0.00232061, -0.000394227, 0.00761978, 0.0190399, 0.0190399};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* v2Kink0510(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.0046225, 0.00633764, 0.0161419, 0.0231465, 0.031357, 0.0387739, 0.0525418, 0.0582742, 0.0787954, 0.0953374, 0.0912154, 0.0912694, 0.100176, 0.115121, 0.080661, 0.124693, 0.0479271};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00226415, 0.00180494, 0.0016965, 0.00176181, 0.0019341, 0.00219713, 0.0019242, 0.00263926, 0.00362214, 0.00502443, 0.00690064, 0.00938605, 0.0103232, 0.0181839, 0.0305616, 0.0480648, 0.0734767};
  Double_t _yerrSyst[] = {0.000486803, 0.00123988, 0.00177792, 0.00240858, 0.00297828, 0.00403581, 0.00447613, 0.00605239, 0.00732301, 0.00700639, 0.00701053, 0.00769464, 0.0109818, 0.00769457, 0.019115, 0.00734709, 0.00734709};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2Kink1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.00749225, 0.0170844, 0.0243086, 0.0381418, 0.052451, 0.0672117, 0.0814359, 0.101393, 0.117803, 0.132852, 0.150789, 0.147415, 0.148044, 0.154442, 0.173818, 0.22373, 0.105331};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00162693, 0.00129069, 0.00120872, 0.00124744, 0.0013639, 0.00154052, 0.00134308, 0.00184132, 0.00253054, 0.0034831, 0.00483312, 0.00652163, 0.00717336, 0.0125484, 0.0208724, 0.0309835, 0.0502633};
  Double_t _yerrSyst[] = {0.00131228, 0.00186718, 0.00292973, 0.00402884, 0.00516263, 0.00625521, 0.00778814, 0.00904864, 0.0102046, 0.0115823, 0.0113232, 0.0113715, 0.0147328, 0.0165812, 0.0342972, 0.016147, 0.016147};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2Kink2030(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.0072888, 0.0241204, 0.038972, 0.0557968, 0.0706379, 0.0905339, 0.110115, 0.137367, 0.154966, 0.173271, 0.187214, 0.197233, 0.196449, 0.204328, 0.190329, 0.198035, 0.215135};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00185852, 0.00147397, 0.00138144, 0.00142629, 0.00155702, 0.00175274, 0.00152882, 0.00209354, 0.00289508, 0.00399616, 0.00545056, 0.00742193, 0.00797779, 0.0139342, 0.0227111, 0.036162, 0.053161};
  Double_t _yerrSyst[] = {0.00185272, 0.00299349, 0.00428583, 0.0054258, 0.00695404, 0.0084581, 0.0105514, 0.0119032, 0.0133092, 0.0143802, 0.0151497, 0.0150896, 0.0194916, 0.0181562, 0.0303583, 0.0329796, 0.0329796};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2Kink3040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.0198535, 0.0269752, 0.0489322, 0.0685815, 0.0871186, 0.106312, 0.130012, 0.153869, 0.181945, 0.191183, 0.200373, 0.217543, 0.22113, 0.227982, 0.152184, 0.204723, 0.10483};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00252247, 0.00200528, 0.00189274, 0.00196066, 0.00214449, 0.00242964, 0.0021233, 0.0029313, 0.00405387, 0.00556614, 0.00764924, 0.0105114, 0.0113762, 0.0191464, 0.0310815, 0.0459809, 0.063268};
  Double_t _yerrSyst[] = {0.00207201, 0.00375855, 0.00526785, 0.0066917, 0.008166, 0.00998641, 0.0118189, 0.0139754, 0.0146851, 0.0153909, 0.0167098, 0.0169853, 0.0217481, 0.0145174, 0.0313835, 0.0160701, 0.0160701};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2Kink4050(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.0184006, 0.039312, 0.0560155, 0.0744549, 0.105477, 0.117306, 0.14288, 0.168705, 0.187757, 0.211535, 0.21727, 0.230958, 0.232356, 0.211124, 0.215568, 0.261458, 0.29527};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00340412, 0.00271914, 0.00258583, 0.00270762, 0.00298639, 0.00338847, 0.00298335, 0.00410873, 0.00569788, 0.00789697, 0.0109218, 0.0148735, 0.0154399, 0.0265265, 0.045388, 0.0613099, 0.0871169};
  Double_t _yerrSyst[] = {0.00301961, 0.00430263, 0.00571899, 0.00810186, 0.00901046, 0.0109748, 0.0129584, 0.0144219, 0.0162483, 0.0166888, 0.0177403, 0.0178476, 0.0201399, 0.0205639, 0.0400808, 0.045264, 0.045264};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2Kink5060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.0184053, 0.0421066, 0.0636918, 0.093209, 0.109357, 0.122793, 0.135598, 0.167734, 0.179235, 0.202347, 0.223503, 0.204778, 0.198038, 0.199902, 0.208692, 0.258469, 0.126576};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00520242, 0.00422132, 0.0040556, 0.00427233, 0.00476878, 0.00545619, 0.00489075, 0.00675269, 0.00942665, 0.0131378, 0.0176433, 0.0228161, 0.025738, 0.0438358, 0.0596872, 0.0951616, 0.126892};
  Double_t _yerrSyst[] = {0.00323427, 0.00489226, 0.00715952, 0.00839988, 0.00943188, 0.0104155, 0.0128839, 0.0137673, 0.0155425, 0.0171676, 0.0157293, 0.0152116, 0.0190695, 0.0199079, 0.0396226, 0.0194037, 0.0194037};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v2Kink6070(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 17;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375, 2.75, 3.25, 3.75, 4.25, 4.75};
  Double_t _y[] = {0.0160665, 0.0369738, 0.0569239, 0.0756837, 0.10208, 0.116458, 0.12661, 0.128646, 0.146764, 0.171994, 0.180963, 0.183109, 0.139119, 0.19161, 0.108737, 0.177263, 0.189475};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00988879, 0.00800023, 0.00786797, 0.00846542, 0.00954037, 0.0110158, 0.00974575, 0.0135531, 0.0191612, 0.0260883, 0.0352833, 0.049656, 0.0488693, 0.0751465, 0.130463, 0.176409, 0.232558};
  Double_t _yerrSyst[] = {0.00284001, 0.00437241, 0.00581338, 0.00784093, 0.00894528, 0.00972507, 0.00988149, 0.0112731, 0.0132111, 0.0139, 0.0140648, 0.010686, 0.0182784, 0.0103728, 0.0271739, 0.0290459, 0.0290459};

  if(!kStat){
      for(Int_t i=0;i<_nPoints;i++){
	  _yerr[i] = 0.0;
	  _xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<_nPoints;i++){
	  Float_t systerr = _yerrSyst[i];
	  _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
      }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// Xi and Omega approved for QM
TGraphAsymmErrors* v2XiSQM110020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2XiSQM110020_nPoints = 8;
  if (last>v2XiSQM110020_nPoints-1) last=v2XiSQM110020_nPoints-1;
  if (last<0 && first<0) last=v2XiSQM110020_nPoints-1;
  if (last<0) last=v2XiSQM110020_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2XiSQM110020_x[] = {1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 5.25};
  Double_t v2XiSQM110020_y[] = {0.0122607, 0.0467335, 0.0787376, 0.103407, 0.127762, 0.14248, 0.138947, 0.158982};
  Double_t v2XiSQM110020_xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
  Double_t v2XiSQM110020_yerrL[] = {0.000674581, 0.00295378, 0.00649944, 0.0107775, 0.016228, 0.0214342, 0.0242101, 0.0239034};
  Double_t v2XiSQM110020_yerrH[] = {0.00196832, 0.00551615, 0.00760563, 0.0086593, 0.00962982, 0.00992244, 0.00911361, 0.00672358};
  Double_t v2XiSQM110020_yerr[] = {0.00695454, 0.00390572, 0.00370856, 0.00461745, 0.00655036, 0.0100333, 0.015608, 0.0198221};

  if(!kStat){
      for(Int_t i=0;i<v2XiSQM110020_nPoints;i++){
	  v2XiSQM110020_yerr[i] = 0.0;
	  v2XiSQM110020_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2XiSQM110020_nPoints;i++){
	  v2XiSQM110020_yerrL[i] = TMath::Sqrt(v2XiSQM110020_yerr[i]*v2XiSQM110020_yerr[i] + v2XiSQM110020_yerrL[i]*v2XiSQM110020_yerrL[i]);
	  v2XiSQM110020_yerrH[i] = TMath::Sqrt(v2XiSQM110020_yerr[i]*v2XiSQM110020_yerr[i] + v2XiSQM110020_yerrH[i]*v2XiSQM110020_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2XiSQM110020_nPoints;i++){
	v2XiSQM110020_yerrL[i] =  v2XiSQM110020_yerr[i];
	v2XiSQM110020_yerrH[i] =  v2XiSQM110020_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2XiSQM110020_x[first], &v2XiSQM110020_y[first], &v2XiSQM110020_xerr[first], &v2XiSQM110020_xerr[first], &v2XiSQM110020_yerrL[first], &v2XiSQM110020_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2XiSQM112040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2XiSQM112040_nPoints = 8;
  if (last>v2XiSQM112040_nPoints-1) last=v2XiSQM112040_nPoints-1;
  if (last<0 && first<0) last=v2XiSQM112040_nPoints-1;
  if (last<0) last=v2XiSQM112040_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2XiSQM112040_x[] = {1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 5.25};
  Double_t v2XiSQM112040_y[] = {0.0712601, 0.112841, 0.172442, 0.213023, 0.244154, 0.264343, 0.301539, 0.250904};
  Double_t v2XiSQM112040_xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
  Double_t v2XiSQM112040_yerrL[] = {0.00177606, 0.00378591, 0.00775254, 0.0122552, 0.0172471, 0.0222139, 0.0294304, 0.0213675};
  Double_t v2XiSQM112040_yerrH[] = {0.00490055, 0.00593952, 0.00764103, 0.00840135, 0.00887289, 0.00905297, 0.00988339, 0.00559612};
  Double_t v2XiSQM112040_yerr[] = {0.00681536, 0.0042213, 0.00417748, 0.00523159, 0.00731666, 0.0110357, 0.0169457, 0.0196909};

  if(!kStat){
      for(Int_t i=0;i<v2XiSQM112040_nPoints;i++){
	  v2XiSQM112040_yerr[i] = 0.0;
	  v2XiSQM112040_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2XiSQM112040_nPoints;i++){
	  v2XiSQM112040_yerrL[i] = TMath::Sqrt(v2XiSQM112040_yerr[i]*v2XiSQM112040_yerr[i] + v2XiSQM112040_yerrL[i]*v2XiSQM112040_yerrL[i]);
	  v2XiSQM112040_yerrH[i] = TMath::Sqrt(v2XiSQM112040_yerr[i]*v2XiSQM112040_yerr[i] + v2XiSQM112040_yerrH[i]*v2XiSQM112040_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2XiSQM112040_nPoints;i++){
	v2XiSQM112040_yerrL[i] =  v2XiSQM112040_yerr[i];
	v2XiSQM112040_yerrH[i] =  v2XiSQM112040_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2XiSQM112040_x[first], &v2XiSQM112040_y[first], &v2XiSQM112040_xerr[first], &v2XiSQM112040_xerr[first], &v2XiSQM112040_yerrL[first], &v2XiSQM112040_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2XiSQM114080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2XiSQM114080_nPoints = 8;
  if (last>v2XiSQM114080_nPoints-1) last=v2XiSQM114080_nPoints-1;
  if (last<0 && first<0) last=v2XiSQM114080_nPoints-1;
  if (last<0) last=v2XiSQM114080_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2XiSQM114080_x[] = {1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 5.25};
  Double_t v2XiSQM114080_y[] = {0.107555, 0.165226, 0.226318, 0.28282, 0.278062, 0.316736, 0.297305, 0.256274};
  Double_t v2XiSQM114080_xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
  Double_t v2XiSQM114080_yerrL[] = {0.0063979, 0.0126659, 0.0216612, 0.0327278, 0.0378896, 0.049769, 0.0529833, 0.0387182};
  Double_t v2XiSQM114080_yerrH[] = {0.00641104, 0.0096893, 0.0131776, 0.0164053, 0.0160933, 0.0183042, 0.0171653, 0.0104354};
  Double_t v2XiSQM114080_yerr[] = {0.011815, 0.00792612, 0.00819157, 0.01045, 0.014573, 0.0213005, 0.0330642, 0.0357233};

  if(!kStat){
    for(Int_t i=0;i<v2XiSQM114080_nPoints;i++){
	  v2XiSQM114080_yerr[i] = 0.0;
	  v2XiSQM114080_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2XiSQM114080_nPoints;i++){
	  v2XiSQM114080_yerrL[i] = TMath::Sqrt(v2XiSQM114080_yerr[i]*v2XiSQM114080_yerr[i] + v2XiSQM114080_yerrL[i]*v2XiSQM114080_yerrL[i]);
	  v2XiSQM114080_yerrH[i] = TMath::Sqrt(v2XiSQM114080_yerr[i]*v2XiSQM114080_yerr[i] + v2XiSQM114080_yerrH[i]*v2XiSQM114080_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2XiSQM114080_nPoints;i++){
	v2XiSQM114080_yerrL[i] =  v2XiSQM114080_yerr[i];
	v2XiSQM114080_yerrH[i] =  v2XiSQM114080_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2XiSQM114080_x[first], &v2XiSQM114080_y[first], &v2XiSQM114080_xerr[first], &v2XiSQM114080_xerr[first], &v2XiSQM114080_yerrL[first], &v2XiSQM114080_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2XiSQM110080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2XiSQM110080_nPoints = 8;
  if (last>v2XiSQM110080_nPoints-1) last=v2XiSQM110080_nPoints-1;
  if (last<0 && first<0) last=v2XiSQM110080_nPoints-1;
  if (last<0) last=v2XiSQM110080_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2XiSQM110080_x[] = {1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 5.25};
  Double_t v2XiSQM110080_y[] = {0.0464933, 0.081694, 0.122039, 0.153568, 0.177207, 0.196881, 0.204888, 0.201544};
  Double_t v2XiSQM110080_xerr[] = {0, 0, 0, 0, 0, 0, 0, 0};
  Double_t v2XiSQM110080_yerrL[] = {0.00215042, 0.00311255, 0.00400468, 0.00489034, 0.00630052, 0.00841482, 0.0106497, 0.00987109};
  Double_t v2XiSQM110080_yerrH[] = {0.00595234, 0.00876142, 0.0112838, 0.0124532, 0.0126323, 0.012266, 0.0110847, 0.00597443};
  Double_t v2XiSQM110080_yerr[] = {0.00463512, 0.00276248, 0.00269241, 0.00337086, 0.00475049, 0.00721294, 0.0111561, 0.0134438};

  if(!kStat){
      for(Int_t i=0;i<v2XiSQM110080_nPoints;i++){
	  v2XiSQM110080_yerr[i] = 0.0;
	  v2XiSQM110080_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2XiSQM110080_nPoints;i++){
	  v2XiSQM110080_yerrL[i] = TMath::Sqrt(v2XiSQM110080_yerr[i]*v2XiSQM110080_yerr[i] + v2XiSQM110080_yerrL[i]*v2XiSQM110080_yerrL[i]);
	  v2XiSQM110080_yerrH[i] = TMath::Sqrt(v2XiSQM110080_yerr[i]*v2XiSQM110080_yerr[i] + v2XiSQM110080_yerrH[i]*v2XiSQM110080_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2XiSQM110080_nPoints;i++){
	v2XiSQM110080_yerrL[i] =  v2XiSQM110080_yerr[i];
	v2XiSQM110080_yerrH[i] =  v2XiSQM110080_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2XiSQM110080_x[first], &v2XiSQM110080_y[first], &v2XiSQM110080_xerr[first], &v2XiSQM110080_xerr[first], &v2XiSQM110080_yerrL[first], &v2XiSQM110080_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2OmegaSQM110020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2OmegaSQM110020_nPoints = 5;
  if (last>v2OmegaSQM110020_nPoints-1) last=v2OmegaSQM110020_nPoints-1;
  if (last<0 && first<0) last=v2OmegaSQM110020_nPoints-1;
  if (last<0) last=v2OmegaSQM110020_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2OmegaSQM110020_x[] = {1.75, 2.25, 2.75, 3.25, 4.25};
  Double_t v2OmegaSQM110020_y[] = {0.0519195, 0.0834295, 0.107549, 0.123927, 0.118231};
  Double_t v2OmegaSQM110020_xerr[] = {0, 0, 0, 0, 0};
  Double_t v2OmegaSQM110020_yerrL[] = {0.00640247, 0.0102868, 0.0132677, 0.0152817, 0.00918066};
  Double_t v2OmegaSQM110020_yerrH[] = {0.00875293, 0.0119141, 0.0138805, 0.0150067, 0.00861946};
  Double_t v2OmegaSQM110020_yerr[] = {0.025909, 0.0176272, 0.0178206, 0.0234692, 0.0248951};

  if(!kStat){
      for(Int_t i=0;i<v2OmegaSQM110020_nPoints;i++){
	  v2OmegaSQM110020_yerr[i] = 0.0;
	  v2OmegaSQM110020_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2OmegaSQM110020_nPoints;i++){
	  v2OmegaSQM110020_yerrL[i] = TMath::Sqrt(v2OmegaSQM110020_yerr[i]*v2OmegaSQM110020_yerr[i] + v2OmegaSQM110020_yerrL[i]*v2OmegaSQM110020_yerrL[i]);
	  v2OmegaSQM110020_yerrH[i] = TMath::Sqrt(v2OmegaSQM110020_yerr[i]*v2OmegaSQM110020_yerr[i] + v2OmegaSQM110020_yerrH[i]*v2OmegaSQM110020_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2OmegaSQM110020_nPoints;i++){
	v2OmegaSQM110020_yerrL[i] =  v2OmegaSQM110020_yerr[i];
	v2OmegaSQM110020_yerrH[i] =  v2OmegaSQM110020_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2OmegaSQM110020_x[first], &v2OmegaSQM110020_y[first], &v2OmegaSQM110020_xerr[first], &v2OmegaSQM110020_xerr[first], &v2OmegaSQM110020_yerrL[first], &v2OmegaSQM110020_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2OmegaSQM112040(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2OmegaSQM112040_nPoints = 5;
  if (last>v2OmegaSQM112040_nPoints-1) last=v2OmegaSQM112040_nPoints-1;
  if (last<0 && first<0) last=v2OmegaSQM112040_nPoints-1;
  if (last<0) last=v2OmegaSQM112040_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2OmegaSQM112040_x[] = {1.75, 2.25, 2.75, 3.25, 4.25};
  Double_t v2OmegaSQM112040_y[] = {0.0814009, 0.12127, 0.185807, 0.195268, 0.283225};
  Double_t v2OmegaSQM112040_xerr[] = {0, 0, 0, 0, 0};
  Double_t v2OmegaSQM112040_yerrL[] = {0.00965916, 0.0143777, 0.0219579, 0.0231921, 0.0196319};
  Double_t v2OmegaSQM112040_yerrH[] = {0.00830306, 0.012386, 0.0189283, 0.0200527, 0.0169729};
  Double_t v2OmegaSQM112040_yerr[] = {0.0234382, 0.0177883, 0.0182779, 0.0242621, 0.023935};

  if(!kStat){
      for(Int_t i=0;i<v2OmegaSQM112040_nPoints;i++){
	  v2OmegaSQM112040_yerr[i] = 0.0;
	  v2OmegaSQM112040_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2OmegaSQM112040_nPoints;i++){
	  v2OmegaSQM112040_yerrL[i] = TMath::Sqrt(v2OmegaSQM112040_yerr[i]*v2OmegaSQM112040_yerr[i] + v2OmegaSQM112040_yerrL[i]*v2OmegaSQM112040_yerrL[i]);
	  v2OmegaSQM112040_yerrH[i] = TMath::Sqrt(v2OmegaSQM112040_yerr[i]*v2OmegaSQM112040_yerr[i] + v2OmegaSQM112040_yerrH[i]*v2OmegaSQM112040_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2OmegaSQM112040_nPoints;i++){
	v2OmegaSQM112040_yerrL[i] =  v2OmegaSQM112040_yerr[i];
	v2OmegaSQM112040_yerrH[i] =  v2OmegaSQM112040_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2OmegaSQM112040_x[first], &v2OmegaSQM112040_y[first], &v2OmegaSQM112040_xerr[first], &v2OmegaSQM112040_xerr[first], &v2OmegaSQM112040_yerrL[first], &v2OmegaSQM112040_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2OmegaSQM114080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2OmegaSQM114080_nPoints = 5;
  if (last>v2OmegaSQM114080_nPoints-1) last=v2OmegaSQM114080_nPoints-1;
  if (last<0 && first<0) last=v2OmegaSQM114080_nPoints-1;
  if (last<0) last=v2OmegaSQM114080_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2OmegaSQM114080_x[] = {1.75, 2.25, 2.75, 3.25, 4.25};
  Double_t v2OmegaSQM114080_y[] = {0.142785, 0.198855, 0.245025, 0.324906, 0.251377};
  Double_t v2OmegaSQM114080_xerr[] = {0, 0, 0, 0, 0};
  Double_t v2OmegaSQM114080_yerrL[] = {0.0144081, 0.0200663, 0.0247293, 0.0327983, 0.0142125};
  Double_t v2OmegaSQM114080_yerrH[] = {0.0149509, 0.0207167, 0.0254687, 0.0337345, 0.0145248};
  Double_t v2OmegaSQM114080_yerr[] = {0.0382144, 0.033091, 0.033689, 0.0448856, 0.0454132};

  if(!kStat){
      for(Int_t i=0;i<v2OmegaSQM114080_nPoints;i++){
	  v2OmegaSQM114080_yerr[i] = 0.0;
	  v2OmegaSQM114080_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2OmegaSQM114080_nPoints;i++){
	  v2OmegaSQM114080_yerrL[i] = TMath::Sqrt(v2OmegaSQM114080_yerr[i]*v2OmegaSQM114080_yerr[i] + v2OmegaSQM114080_yerrL[i]*v2OmegaSQM114080_yerrL[i]);
	  v2OmegaSQM114080_yerrH[i] = TMath::Sqrt(v2OmegaSQM114080_yerr[i]*v2OmegaSQM114080_yerr[i] + v2OmegaSQM114080_yerrH[i]*v2OmegaSQM114080_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2OmegaSQM114080_nPoints;i++){
	v2OmegaSQM114080_yerrL[i] =  v2OmegaSQM114080_yerr[i];
	v2OmegaSQM114080_yerrH[i] =  v2OmegaSQM114080_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2OmegaSQM114080_x[first], &v2OmegaSQM114080_y[first], &v2OmegaSQM114080_xerr[first], &v2OmegaSQM114080_xerr[first], &v2OmegaSQM114080_yerrL[first], &v2OmegaSQM114080_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2OmegaSQM110080(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v2OmegaSQM110080_nPoints = 5;
  if (last>v2OmegaSQM110080_nPoints-1) last=v2OmegaSQM110080_nPoints-1;
  if (last<0 && first<0) last=v2OmegaSQM110080_nPoints-1;
  if (last<0) last=v2OmegaSQM110080_nPoints-1+last;
  if (first<0) first=0;
  Double_t v2OmegaSQM110080_x[] = {1.75, 2.25, 2.75, 3.25, 4.25};
  Double_t v2OmegaSQM110080_y[] = {0.0728013, 0.105837, 0.141493, 0.162115, 0.17853};
  Double_t v2OmegaSQM110080_xerr[] = {0, 0, 0, 0, 0};
  Double_t v2OmegaSQM110080_yerrL[] = {0.00476832, 0.00692521, 0.00925901, 0.0106143, 0.00707741};
  Double_t v2OmegaSQM110080_yerrH[] = {0.00806897, 0.0106351, 0.0134897, 0.0150059, 0.00969192};
  Double_t v2OmegaSQM110080_yerr[] = {0.0170676, 0.0124432, 0.0127703, 0.0168486, 0.0173784};

  if(!kStat){
      for(Int_t i=0;i<v2OmegaSQM110080_nPoints;i++){
	  v2OmegaSQM110080_yerr[i] = 0.0;
	  v2OmegaSQM110080_xerr[i] = 0.05;
      }
  }

  if(kSyst){
      for(Int_t i=0;i<v2OmegaSQM110080_nPoints;i++){
	  v2OmegaSQM110080_yerrL[i] = TMath::Sqrt(v2OmegaSQM110080_yerr[i]*v2OmegaSQM110080_yerr[i] + v2OmegaSQM110080_yerrL[i]*v2OmegaSQM110080_yerrL[i]);
	  v2OmegaSQM110080_yerrH[i] = TMath::Sqrt(v2OmegaSQM110080_yerr[i]*v2OmegaSQM110080_yerr[i] + v2OmegaSQM110080_yerrH[i]*v2OmegaSQM110080_yerrH[i]);
      }
  }
  else{
      for(Int_t i=0;i<v2OmegaSQM110080_nPoints;i++){
	v2OmegaSQM110080_yerrL[i] =  v2OmegaSQM110080_yerr[i];
	v2OmegaSQM110080_yerrH[i] =  v2OmegaSQM110080_yerr[i];
      }
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &v2OmegaSQM110080_x[first], &v2OmegaSQM110080_y[first], &v2OmegaSQM110080_xerr[first], &v2OmegaSQM110080_xerr[first], &v2OmegaSQM110080_yerrL[first], &v2OmegaSQM110080_yerrH[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// new phenix results
//run7_pion_v2_0_10_new_sys.txt
TGraphErrors* v2PionPHENIX0010(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 16;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.08603, 1.28304, 1.48005, 1.67805, 1.87506, 2.07207, 2.27106, 2.46906, 2.66607, 2.86506, 3.16899, 3.66399, 4.16493, 4.6629, 5.16186, 5.66379};
  Double_t _y[] = {0.0447245, 0.0512504, 0.0574795, 0.0617296, 0.0659824, 0.067367, 0.0678192, 0.0700349, 0.0671657, 0.0659291, 0.0683715, 0.0631313, 0.0549763, 0.0721619, 0.020939, 0.0417612};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000222701, 0.000298047, 0.000398824, 0.000531302, 0.000704507, 0.000921466, 0.00122938, 0.0016748, 0.00236021, 0.00362083, 0.00317037, 0.00518788, 0.00802641, 0.0112411, 0.0160144, 0.0189578};
  Double_t _yerrS[] = {0.00232395, 0.00266305, 0.00298672, 0.00320756, 0.00342855, 0.00350049, 0.00352399, 0.00363912, 0.00349003, 0.00342578, 0.00355269, 0.0032804, 0.00315815, 0.00510262, 0.00148061, 0.0033927};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

//run7_pion_v2_10_20_new_sys.txt
TGraphErrors* v2PionPHENIX1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 16;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.08055, 1.27656, 1.47258, 1.66958, 1.86559, 2.0616, 2.25959, 2.45659, 2.65261, 2.85059, 3.15298, 3.64548, 4.14389, 4.63935, 5.13579, 5.63519};
  Double_t _y[] = {0.0838016, 0.0968691, 0.106677, 0.114742, 0.121115, 0.125837, 0.128677, 0.127835, 0.131393, 0.131373, 0.123985, 0.123396, 0.115131, 0.0932285, 0.0945603, 0.0750595};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000190244, 0.000253536, 0.000337772, 0.000447675, 0.000590113, 0.000766905, 0.00101727, 0.00137521, 0.00193639, 0.00295942, 0.00247985, 0.0040054, 0.00611423, 0.00866822, 0.0126883, 0.0155469};
  Double_t _yerrS[] = {0.00435446, 0.00503347, 0.0055431, 0.00596217, 0.00629332, 0.00653868, 0.00668625, 0.0066425, 0.00682738, 0.00682634, 0.00644245, 0.00641184, 0.00661377, 0.00659225, 0.00668642, 0.00609786};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

//run7_pion_v2_40_60_new_sys.txt
TGraphErrors* v2PionPHENIX4060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 16;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.07396, 1.26878, 1.46361, 1.65941, 1.85423, 2.04905, 2.24583, 2.44163, 2.63645, 2.83323, 3.13378, 3.62328, 4.11865, 4.61109, 5.10451, 5.60086};
  Double_t _y[] = {0.146249, 0.164058, 0.177614, 0.187203, 0.19349, 0.197953, 0.199169, 0.200115, 0.192651, 0.195101, 0.193315, 0.175433, 0.181669, 0.170734, 0.138554, 0.166739};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000319726, 0.000427822, 0.000567271, 0.000746729, 0.000975056, 0.00125177, 0.00162971, 0.00214835, 0.0029677, 0.00413081, 0.00363466, 0.00567607, 0.00854216, 0.011881, 0.0183544, 0.0223662};
  Double_t _yerrS[] = {0.00485053, 0.00544119, 0.00589079, 0.00620882, 0.00641734, 0.00656536, 0.00660569, 0.00663706, 0.00638951, 0.00647077, 0.00641153, 0.00581845, 0.0061446, 0.00580747, 0.0081134, 0.0104448};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// run7_kaon_v2_0_10_new_sys.txt
TGraphErrors* v2KaonPHENIX0010(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 12;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.08822, 1.28563, 1.48403, 1.68144, 1.87885, 2.07626, 2.27565, 2.47405, 2.67146, 2.87085, 3.17539, 3.67139};
  Double_t _y[] = {0.0357656, 0.0445811, 0.0510559, 0.056433, 0.0611297, 0.0651863, 0.0664281, 0.0659084, 0.0686639, 0.0618081, 0.0691171, 0.070925};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000486801, 0.000587744, 0.000725444, 0.000904875, 0.00111919, 0.00136158, 0.00169396, 0.00216426, 0.00273625, 0.00346302, 0.00727094, 0.0150273};
  Double_t _yerrS[] = {0.00185844, 0.0023165, 0.00265294, 0.00293234, 0.00317639, 0.00338718, 0.00345171, 0.0034247, 0.00356788, 0.00321164, 0.0037857, 0.00459647};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// run7_kaon_v2_10_20_new_sys.txt
TGraphErrors* v2KaonPHENIX1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 12;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.08164, 1.27786, 1.47506, 1.67127, 1.86748, 2.0637, 2.26188, 2.45908, 2.6553, 2.85348, 3.15619, 3.64919};
  Double_t _y[] = {0.0692813, 0.0828747, 0.0962352, 0.105263, 0.112471, 0.120031, 0.125193, 0.124436, 0.125841, 0.131877, 0.120636, 0.105541};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.0004194, 0.000504879, 0.000621801, 0.000772546, 0.000952328, 0.00115052, 0.00142343, 0.00179927, 0.00225918, 0.002827, 0.00553934, 0.0113803};
  Double_t _yerrS[] = {0.00359996, 0.0043063, 0.00500053, 0.00546963, 0.00584416, 0.00623699, 0.00650522, 0.00646588, 0.00653889, 0.00685253, 0.00660751, 0.00683984};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
// run7_kaon_v2_40_60_new_sys.txt
TGraphErrors* v2KaonPHENIX4060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 12;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.07506, 1.27008, 1.4651, 1.6611, 1.85612, 2.05114, 2.24812, 2.44412, 2.63914, 2.83612, 3.13698, 3.62698};
  Double_t _y[] = {0.124107, 0.144895, 0.161632, 0.171836, 0.182351, 0.184241, 0.188441, 0.189411, 0.188685, 0.182045, 0.178534, 0.211813};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000736268, 0.00089208, 0.0011, 0.00136603, 0.00167614, 0.00200908, 0.00243513, 0.00303598, 0.0037641, 0.00468641, 0.00828064, 0.0164834};
  Double_t _yerrS[] = {0.00411616, 0.00480562, 0.00536073, 0.00569916, 0.0060479, 0.00611058, 0.00624988, 0.00628205, 0.00625797, 0.00603775, 0.00668013, 0.00792532};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// run7_proton_v2_0_10_new_sys.txt
TGraphErrors* v2AntiprotonPHENIX0010(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 16;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.08932, 1.28693, 1.48454, 1.68314, 1.88074, 2.07835, 2.27794, 2.47654, 2.67415, 2.87374, 3.17859, 3.67509, 4.17755, 4.67703, 5.1775, 5.68095};
  Double_t _y[] = {0.0210294, 0.0324285, 0.0433725, 0.054565, 0.0644553, 0.0715762, 0.0779503, 0.084824, 0.0906549, 0.0951113, 0.0942727, 0.115555, 0.106061, 0.126451, 0.0902255, 0.120818};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000434053, 0.000480734, 0.0005626, 0.000684719, 0.000846039, 0.00106283, 0.00136236, 0.00173081, 0.00220325, 0.00280525, 0.00415681, 0.00777055, 0.013864, 0.0230551, 0.036464, 0.0527073};
  Double_t _yerrS[] = {0.00109272, 0.00168503, 0.0022537, 0.00283528, 0.0033492, 0.00371921, 0.00405042, 0.00440758, 0.00471057, 0.00494213, 0.00516353, 0.00632921, 0.00609274, 0.00779496, 0.00637991, 0.0098153};


  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// run7_proton_v2_10_20_new_sys.txt
TGraphErrors* v2AntiprotonPHENIX1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 16;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.08274, 1.27915, 1.47557, 1.67297, 1.86938, 2.06579, 2.26418, 2.46158, 2.65799, 2.85638, 3.15939, 3.65289, 4.15231, 4.64877, 5.14622, 5.64663};
  Double_t _y[] = {0.0463211, 0.0659061, 0.0868213, 0.103327, 0.120261, 0.133348, 0.145277, 0.15602, 0.164359, 0.169368, 0.170511, 0.179104, 0.158896, 0.159868, 0.166116, 0.130806};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.00037414, 0.000415018, 0.000486394, 0.000591698, 0.000730528, 0.00091517, 0.00116734, 0.00147657, 0.00186563, 0.00235995, 0.00329999, 0.00605461, 0.0105693, 0.0175421, 0.0284934, 0.0418104};
  Double_t _yerrS[] = {0.00240691, 0.00342458, 0.00451137, 0.00536903, 0.00624894, 0.00692897, 0.00754881, 0.00810704, 0.00854034, 0.00880062, 0.00933927, 0.00980993, 0.00912788, 0.00985493, 0.0117462, 0.0106267};


  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }
  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

// run7_proton_v2_40_60_new_sys.txt
TGraphErrors* v2AntiprotonPHENIX4060(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t _nPoints = 16;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {1.07616, 1.27138, 1.4666, 1.6628, 1.85801, 2.05323, 2.25041, 2.44661, 2.64183, 2.83901, 3.14018, 3.63068, 4.12707, 4.62051, 5.11493, 5.6123};
  Double_t _y[] = {0.111728, 0.142026, 0.168396, 0.192156, 0.212232, 0.230724, 0.24573, 0.248557, 0.254284, 0.262583, 0.26775, 0.267792, 0.26648, 0.226501, 0.217841, 0.20696};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0.000627176, 0.00072307, 0.00087183, 0.0010793, 0.00135044, 0.00169897, 0.00217077, 0.00272808, 0.00343839, 0.00428991, 0.00552279, 0.00973688, 0.0163398, 0.0270507, 0.0416056, 0.0595448};
  Double_t _yerrS[] = {0.0037056, 0.00471047, 0.00558506, 0.00637309, 0.00703894, 0.00765225, 0.00814994, 0.0082437, 0.00843365, 0.00870889, 0.00937125, 0.00987931, 0.010307, 0.00902884, 0.012994, 0.0142622};

  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _xerr[i] = 0.05;
    }
  }

  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _yerrS[i]*_yerrS[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &_x[first], &_y[first], &_xerr[first], &_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}


// v3 from QM11 QC{2}
TGraphErrors* v32_1020_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //  commentme
    Int_t v32_1020_pion_nPoints = 30;
  if (last>v32_1020_pion_nPoints-1) last=v32_1020_pion_nPoints-1;
  if (last<0 && first<0) last=v32_1020_pion_nPoints-1;
  if (last<0) last=v32_1020_pion_nPoints-1+last;
  if (first<0) first=0;
  Double_t v32_1020_pion_x[] = {0.33, 0.43, 0.53, 0.63, 0.73, 0.83, 0.93, 1.03, 1.13, 1.23, 1.33, 1.43, 1.53, 1.63, 1.73, 1.83, 1.93, 2.18, 2.58, 2.98, 3.38, 3.78, 4.18, 4.58, 4.98, 5.38, 5.78, 6.18, 6.58, 6.98};
  Double_t v32_1020_pion_y[] = {0.0134229, 0.0167432, 0.02107, 0.0262901, 0.0312443, 0.0362567, 0.0422287, 0.0474885, 0.0524595, 0.0567917, 0.0625859, 0.0673443, 0.0725898, 0.0765344, 0.0815771, 0.0849598, 0.08939, 0.0941734, 0.107589, 0.11532, 0.112264, 0.114684, 0.116536, 0, 0, 0, 0, 0, 0, 0};
  Double_t v32_1020_pion_xerr[] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v32_1020_pion_yerr[] = {0.000816086, 0.000658693, 0.000742595, 0.000887094, 0.0010371, 0.00119473, 0.00138196, 0.00155296, 0.00171972, 0.00187557, 0.00207941, 0.00226254, 0.00247083, 0.0026593, 0.00289646, 0.00311194, 0.00338275, 0.00316069, 0.00400211, 0.00508831, 0.00656859, 0.00987805, 0.0189298, 0, 0, 0, 0, 0, 0, 0};

  Double_t v32_1020_pion_statError_yerr[] = {0.000390094, 0.00023236, 0.000217013, 0.000228925, 0.00025102, 0.000280217, 0.000313372, 0.000350632, 0.000393107, 0.000444735, 0.000506531, 0.00057781, 0.000661905, 0.00076183, 0.000878565, 0.00101342, 0.00117191, 0.000806488, 0.00134727, 0.00212538, 0.0032084, 0.00527795, 0.0105326, 0, 0, 0, 0, 0, 0, 0};

  if(kStat && !kSyst){
    for(Int_t i=0;i < v32_1020_pion_nPoints;i++){
      v32_1020_pion_yerr[i] =   v32_1020_pion_statError_yerr[i];
    }
  }
  else if(!kStat && kSyst){
    for(Int_t i=0;i < v32_1020_pion_nPoints;i++){
      v32_1020_pion_yerr[i] =  v32_1020_pion_yerr[i]*v32_1020_pion_yerr[i] -  v32_1020_pion_statError_yerr[i]* v32_1020_pion_statError_yerr[i];
      v32_1020_pion_yerr[i] = TMath::Sqrt(v32_1020_pion_yerr[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &v32_1020_pion_x[first], &v32_1020_pion_y[first], &v32_1020_pion_xerr[first], &v32_1020_pion_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v32_1020_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v32_1020_kaon_nPoints = 25;
  if (last>v32_1020_kaon_nPoints-1) last=v32_1020_kaon_nPoints-1;
  if (last<0 && first<0) last=v32_1020_kaon_nPoints-1;
  if (last<0) last=v32_1020_kaon_nPoints-1+last;
  if (first<0) first=0;
  Double_t v32_1020_kaon_x[] = {0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.2, 2.6, 3, 3.4, 3.8, 4.2, 4.6, 5, 5.4};
  Double_t v32_1020_kaon_y[] = {0.00302743, 0.00670068, 0.00960944, 0.0140959, 0.0210006, 0.0261218, 0.0331506, 0.0390413, 0.0432323, 0.0498768, 0.0535018, 0.0609539, 0.0663136, 0.0701533, 0.0746404, 0.0797591, 0.0864142, 0.0969327, 0.117939, 0, 0, 0, 0, 0, 0};
  Double_t v32_1020_kaon_xerr[] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v32_1020_kaon_yerr[] = {0.00291973, 0.00189558, 0.00156194, 0.00147753, 0.00151666, 0.00159995, 0.00174912, 0.00190683, 0.0020656, 0.00228462, 0.00246573, 0.00273975, 0.00300067, 0.00325875, 0.00355464, 0.00392612, 0.00336279, 0.00513077, 0.0134291, 0, 0, 0, 0, 0, 0};

  Double_t v32_1020_kaon_statError_yerr[] = {0.00162595, 0.00105604, 0.000863028, 0.000799207, 0.000779598, 0.000790079, 0.000816257, 0.000853837, 0.000912805, 0.000980786, 0.00106406, 0.00115936, 0.00127764, 0.00141682, 0.00157074, 0.00177344, 0.00121994, 0.00240937, 0.00737108, 0, 0, 0, 0, 0, 0};

  if(kStat && !kSyst){
    for(Int_t i=0;i < v32_1020_kaon_nPoints;i++){
      v32_1020_kaon_yerr[i] =   v32_1020_kaon_statError_yerr[i];
    }
  }
  else if(!kStat && kSyst){
    for(Int_t i=0;i < v32_1020_kaon_nPoints;i++){
      v32_1020_kaon_yerr[i] =  v32_1020_kaon_yerr[i]*v32_1020_kaon_yerr[i] -  v32_1020_kaon_statError_yerr[i]* v32_1020_kaon_statError_yerr[i];
      v32_1020_kaon_yerr[i] = TMath::Sqrt(v32_1020_kaon_yerr[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &v32_1020_kaon_x[first], &v32_1020_kaon_y[first], &v32_1020_kaon_xerr[first], &v32_1020_kaon_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v32_1020_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v32_1020_antiproton_nPoints = 28;
  if (last>v32_1020_antiproton_nPoints-1) last=v32_1020_antiproton_nPoints-1;
  if (last<0 && first<0) last=v32_1020_antiproton_nPoints-1;
  if (last<0) last=v32_1020_antiproton_nPoints-1+last;
  if (first<0) first=0;
  Double_t v32_1020_antiproton_x[] = {0.57, 0.67, 0.77, 0.87, 0.97, 1.07, 1.17, 1.27, 1.37, 1.47, 1.57, 1.67, 1.77, 1.87, 1.97, 2.22, 2.62, 3.02, 3.42, 3.82, 4.22, 4.62, 5.02, 5.42, 5.82, 6.22, 6.62, 7.02};
  Double_t v32_1020_antiproton_y[] = {0.00186436, 7.05025e-05, 0.00226371, 0.0025971, 0.00426813, 0.00596588, 0.0109206, 0.0170278, 0.0240831, 0.0275234, 0.0358475, 0.0415758, 0.0496074, 0.0573969, 0.0630463, 0.0785925, 0.102838, 0.120118, 0.137493, 0.142522, 0.151444, 0.144125, 0.192062, 0, 0, 0, 0, 0};
  Double_t v32_1020_antiproton_xerr[] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v32_1020_antiproton_yerr[] = {0.00661632, 0.00481615, 0.00405812, 0.00364377, 0.00344754, 0.0033137, 0.0032578, 0.00326806, 0.0033384, 0.00343628, 0.00359836, 0.00377852, 0.00402085, 0.00431005, 0.00458806, 0.00340917, 0.00454227, 0.00595661, 0.00811159, 0.0117587, 0.0190361, 0.0338669, 0.134234, 0, 0, 0, 0, 0};

  Double_t v32_1020_antiproton_statError_yerr[] = {0.00377499, 0.00276915, 0.00234124, 0.00210598, 0.00199562, 0.00191919, 0.00188437, 0.00187885, 0.00189756, 0.00194448, 0.00200543, 0.0020857, 0.00218184, 0.00231353, 0.00244858, 0.00144405, 0.00195883, 0.00279008, 0.00411264, 0.00646947, 0.0109628, 0.0199308, 0.0781864, 0, 0, 0, 0, 0};

  if(kStat && !kSyst){
    for(Int_t i=0;i < v32_1020_antiproton_nPoints;i++){
      v32_1020_antiproton_yerr[i] =   v32_1020_antiproton_statError_yerr[i];
    }
  }
  else if(!kStat && kSyst){
    for(Int_t i=0;i < v32_1020_antiproton_nPoints;i++){
      v32_1020_antiproton_yerr[i] =  v32_1020_antiproton_yerr[i]*v32_1020_antiproton_yerr[i] -  v32_1020_antiproton_statError_yerr[i]* v32_1020_antiproton_statError_yerr[i];
      v32_1020_antiproton_yerr[i] = TMath::Sqrt(v32_1020_antiproton_yerr[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &v32_1020_antiproton_x[first], &v32_1020_antiproton_y[first], &v32_1020_antiproton_xerr[first], &v32_1020_antiproton_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v32_4050_pion(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v32_4050_pion_nPoints = 30;
  if (last>v32_4050_pion_nPoints-1) last=v32_4050_pion_nPoints-1;
  if (last<0 && first<0) last=v32_4050_pion_nPoints-1;
  if (last<0) last=v32_4050_pion_nPoints-1+last;
  if (first<0) first=0;
  Double_t v32_4050_pion_x[] = {0.33, 0.43, 0.53, 0.63, 0.73, 0.83, 0.93, 1.03, 1.13, 1.23, 1.33, 1.43, 1.53, 1.63, 1.73, 1.83, 1.93, 2.18, 2.58, 2.98, 3.38, 3.78, 4.18, 4.58, 4.98, 5.38, 5.78, 6.18, 6.58, 6.98};
  Double_t v32_4050_pion_y[] = {0.0169421, 0.0226956, 0.027289, 0.0344464, 0.0401389, 0.0483862, 0.0538993, 0.0615618, 0.0668474, 0.0724796, 0.0775622, 0.0842534, 0.0894518, 0.0928943, 0.0977908, 0.106364, 0.107759, 0.114188, 0.121071, 0.121902, 0.129527, 0.134133, 0.180864, 0, 0, 0, 0, 0, 0, 0};
  Double_t v32_4050_pion_xerr[] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v32_4050_pion_yerr[] = {0.00182104, 0.00124715, 0.00128606, 0.00147443, 0.0016766, 0.00195826, 0.00219382, 0.00249919, 0.00277752, 0.00307774, 0.0033925, 0.00375772, 0.00411319, 0.00447598, 0.00490496, 0.00544854, 0.00595246, 0.00477098, 0.00662221, 0.00943227, 0.0134523, 0.0209874, 0.038741, 0, 0, 0, 0, 0, 0, 0};

  Double_t v32_4050_pion_statError_yerr[] = {0.000970601, 0.000579259, 0.00055439, 0.000591078, 0.000657755, 0.000742217, 0.000836199, 0.000950415, 0.00108627, 0.00123124, 0.00139699, 0.00157447, 0.00176577, 0.00198297, 0.00222574, 0.00250362, 0.00283409, 0.00188237, 0.00314542, 0.00493338, 0.00728232, 0.0117328, 0.0214989, 0, 0, 0, 0, 0, 0, 0};


  if(kStat && !kSyst){
    for(Int_t i=0;i < v32_4050_pion_nPoints;i++){
      v32_4050_pion_yerr[i] =   v32_4050_pion_statError_yerr[i];
    }
  }
  else if(!kStat && kSyst){
    for(Int_t i=0;i < v32_4050_pion_nPoints;i++){
      v32_4050_pion_yerr[i] =  v32_4050_pion_yerr[i]*v32_4050_pion_yerr[i] -  v32_4050_pion_statError_yerr[i]* v32_4050_pion_statError_yerr[i];
      v32_4050_pion_yerr[i] = TMath::Sqrt(v32_4050_pion_yerr[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &v32_4050_pion_x[first], &v32_4050_pion_y[first], &v32_4050_pion_xerr[first], &v32_4050_pion_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v32_4050_kaon(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v32_4050_kaon_nPoints = 25;
  if (last>v32_4050_kaon_nPoints-1) last=v32_4050_kaon_nPoints-1;
  if (last<0 && first<0) last=v32_4050_kaon_nPoints-1;
  if (last<0) last=v32_4050_kaon_nPoints-1+last;
  if (first<0) first=0;
  Double_t v32_4050_kaon_x[] = {0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.2, 2.6, 3, 3.4, 3.8, 4.2, 4.6, 5, 5.4};
  Double_t v32_4050_kaon_y[] = {0.00800406, 0.00824057, 0.0184803, 0.0236651, 0.0334773, 0.0394258, 0.042902, 0.0526186, 0.056656, 0.0609829, 0.0704108, 0.0727429, 0.0780242, 0.0906113, 0.0879405, 0.0954966, 0.101673, 0.110859, 0.12648, 0, 0, 0, 0, 0, 0};
  Double_t v32_4050_kaon_xerr[] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v32_4050_kaon_yerr[] = {0.00647827, 0.00440681, 0.00382398, 0.00365843, 0.00370453, 0.0038315, 0.00398085, 0.00424222, 0.00450969, 0.00481112, 0.00520842, 0.00558284, 0.00606004, 0.00670006, 0.00727355, 0.00812917, 0.00608355, 0.0109859, 0.0319517, 0, 0, 0, 0, 0, 0};

  Double_t v32_4050_kaon_statError_yerr[] = {0.00360466, 0.0024692, 0.00212939, 0.0020259, 0.00201447, 0.00206226, 0.00213368, 0.0022294, 0.00236744, 0.00252395, 0.0026997, 0.00291806, 0.00317155, 0.00347944, 0.00385033, 0.00432545, 0.00299177, 0.00595115, 0.0180812, 0, 0, 0, 0, 0, 0};

  if(kStat && !kSyst){
    for(Int_t i=0;i < v32_4050_kaon_nPoints;i++){
      v32_4050_kaon_yerr[i] =   v32_4050_kaon_statError_yerr[i];
    }
  }
  else if(!kStat && kSyst){
    for(Int_t i=0;i < v32_4050_kaon_nPoints;i++){
      v32_4050_kaon_yerr[i] =  v32_4050_kaon_yerr[i]*v32_4050_kaon_yerr[i] -  v32_4050_kaon_statError_yerr[i]* v32_4050_kaon_statError_yerr[i];
      v32_4050_kaon_yerr[i] = TMath::Sqrt(v32_4050_kaon_yerr[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &v32_4050_kaon_x[first], &v32_4050_kaon_y[first], &v32_4050_kaon_xerr[first], &v32_4050_kaon_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphErrors* v32_4050_antiproton(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1)
{
  //commentme
  Int_t v32_4050_antiproton_nPoints = 28;
  if (last>v32_4050_antiproton_nPoints-1) last=v32_4050_antiproton_nPoints-1;
  if (last<0 && first<0) last=v32_4050_antiproton_nPoints-1;
  if (last<0) last=v32_4050_antiproton_nPoints-1+last;
  if (first<0) first=0;
  Double_t v32_4050_antiproton_x[] = {0.57, 0.67, 0.77, 0.87, 0.97, 1.07, 1.17, 1.27, 1.37, 1.47, 1.57, 1.67, 1.77, 1.87, 1.97, 2.22, 2.62, 3.02, 3.42, 3.82, 4.22, 4.62, 5.02, 5.42, 5.82, 6.22, 6.62, 7.02};
  Double_t v32_4050_antiproton_y[] = {0.0088111, -0.00168695, 0.010351, 0.00911283, 0.0276992, 0.0240616, 0.0326823, 0.0380326, 0.0510652, 0.0538593, 0.0628152, 0.074492, 0.0831343, 0.0866683, 0.094179, 0.102877, 0.128451, 0.145513, 0.151148, 0.171606, 0.194494, 0.147614, -0.143488, 0, 0, 0, 0, 0};
  Double_t v32_4050_antiproton_xerr[] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04};
  Double_t v32_4050_antiproton_yerr[] = {0.013943, 0.0098319, 0.0083454, 0.00762029, 0.00734332, 0.00721975, 0.00728574, 0.00737796, 0.0076756, 0.00797667, 0.00841186, 0.00894512, 0.00951378, 0.0101552, 0.0109422, 0.0070593, 0.00971959, 0.0135803, 0.0190688, 0.0278606, 0.045911, 0.0816557, 0.29937, 0, 0, 0, 0, 0};

  Double_t v32_4050_antiproton_statError_yerr[] = {0.00790802, 0.00562316, 0.00479201, 0.00440497, 0.00422239, 0.00416259, 0.00418901, 0.00423755, 0.00437205, 0.00454936, 0.00478116, 0.00504625, 0.00534515, 0.00573639, 0.00617799, 0.00370405, 0.00521801, 0.00753314, 0.0108532, 0.0161127, 0.0265911, 0.0483794, 0.176792, 0, 0, 0, 0, 0};

  if(kStat && !kSyst){
    for(Int_t i=0;i < v32_4050_antiproton_nPoints;i++){
      v32_4050_antiproton_yerr[i] =   v32_4050_antiproton_statError_yerr[i];
    }
  }
  else if(!kStat && kSyst){
    for(Int_t i=0;i < v32_4050_antiproton_nPoints;i++){
      v32_4050_antiproton_yerr[i] =  v32_4050_antiproton_yerr[i]*v32_4050_antiproton_yerr[i] -  v32_4050_antiproton_statError_yerr[i]* v32_4050_antiproton_statError_yerr[i];
      v32_4050_antiproton_yerr[i] = TMath::Sqrt(v32_4050_antiproton_yerr[i]);
    }
  }

  TGraphErrors* graph = new TGraphErrors(last-first+1, &v32_4050_antiproton_x[first], &v32_4050_antiproton_y[first], &v32_4050_antiproton_xerr[first], &v32_4050_antiproton_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}


