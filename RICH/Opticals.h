#ifdef __CINT__
void Opticals()
{
  gROOT->Reset();
#endif   
  int i;
  const Int_t kNbins=30;
  
  Float_t aPckov[kNbins];  for(i=0;i<kNbins;i++) aPckov[i]=(0.1*i+5.5)*1e-9;//Photons energy bins 5.5 eV - 8.5 eV step 0.1 eV
    
  
  Float_t aIdxC6F14[kNbins];  
  Float_t aIdxCH4[kNbins];
  Float_t aIdxGrid[kNbins];
  Float_t aIdxSiO2[kNbins];
  Float_t aIdxOpSiO2[kNbins];
        
  for (i=0;i<kNbins;i++){
    aIdxC6F14[i] = (Float_t)AliRICHParam::IndOfRefC6F14(aPckov[i]*1e9);
    aIdxSiO2[i]  = (Float_t)AliRICHParam::IndOfRefSiO2(aPckov[i]*1e9);
    aIdxCH4[i]   = (Float_t)AliRICHParam::IndOfRefCH4();
    aIdxGrid[i]    =1;
    aIdxOpSiO2[i]  =1;
  } 
    
  Float_t aAbsC6F14[kNbins]={//New values from A.DiMauro 28.10.03 total 31
    32701.4219, 17996.1141, 10039.7281, 1799.1230, 1799.1231, 1799.1231, 1241.4091, 179.0987, 179.0986, 179.0987,
      179.0987,   118.9800,    39.5058,   23.7244,   11.1283,    7.1573,    3.6249,   2.1236,   0.7362,   0.5348,
        0.3387,     0.3074,     0.3050,    0.0001,    0.0001,    0.0001,    0.0001,   0.0001,   0.0001,   0.0001};    
  Float_t aAbsSiO2[kNbins]={//New values from A.DiMauro 28.10.03 total 31
       34.4338, 30.5424, 30.2584, 31.4928, 31.7868, 17.8397, 9.3410, 6.4492, 6.1128, 5.8128,
        5.5589,  5.2877,  5.0162,  4.7999,  4.5734,  4.2135, 3.7471, 2.6033, 1.5223, 0.9658,
        0.4242,  0.2500,  0.1426,  0.0863,  0.0793,  0.0724, 0.0655, 0.0587, 0.0001, 0.0001};
    
  Float_t aAbsOpSiO2[kNbins];
  Float_t aAbsCH4[kNbins];
  Float_t aAbsGrid[kNbins];
  Float_t aAbsCsI[kNbins];
    for (i=0;i<kNbins;i++){
      aAbsCH4[i]         =AbsoCH4(aPckov[i]*1e9); 
      aAbsOpSiO2[i]      =1e-5; 
      aAbsCsI[i]         =1e-4; 
      aAbsGrid[i]        =1e-4; 
    }
    
  Float_t aQeCsI[kNbins] = {//New values from A.DiMauro 28.10.03 total 31
                            0.0002, 0.0006, 0.0007, 0.0010, 0.0049, 0.0073, 0.0104, 0.0519, 0.0936, 0.1299,
                            0.1560, 0.1768, 0.1872, 0.1976, 0.2142, 0.2288, 0.2434, 0.2599, 0.2673, 0.2808,
                            0.2859, 0.2954, 0.3016, 0.3120, 0.3172, 0.3224, 0.3266, 0.3328, 0.3359, 0.3390};
//                             0.3431};

  
  Float_t aQeAll[kNbins];
  for(i=0;i<kNbins;i++){
    aQeCsI[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeAll[i]=1; //QE for all other materials except for PC must be 1.
  }
       
                    
#ifdef __CINT__

  Float_t aAbsC6F14old[kNbins]={//previous values 26 in total added 0.0001 to the 30
      179.0987, 179.0987, 179.0987, 179.0987, 179.0987, 179.0987, 179.0987, 179.0987, 179.0987, 142.9206, 
       56.6496,  25.5862,  13.9529,  12.0391,  10.4295,   8.8042,   7.0690,   4.4611,   2.0284,   1.2930, 
        0.5772,   0.4074,   0.3350,   0.0001,   0.0001,   0.0001,   0.0001,   0.0001,   0.0001,   0.0001};
  Float_t aAbsSiO2old[kNbins]={//previous values 26 in total added 0.0001 to the 30
      105.8000,  65.5200,  48.5800,  42.8500,  35.7900,  31.2620,  28.5980,  27.5270,  25.0070,  22.8150, 
       21.0040,  19.2660,  17.5250,  15.8780,  14.1770,  11.7190,   9.2820,   6.6200,   4.0930,   2.6010, 
        1.1490,   0.6670,   0.3630,   0.1920,   0.1500,   0.1090,   0.0001,   0.0001,   0.0001,   0.0001};
  Float_t aQeCsIold[kNbins]={
    0.0002, 0.0006, 0.0007, 0.0050, 0.0075, 0.0101, 0.0243, 0.0405, 0.0689, 0.1053, 
    0.1215, 0.1417, 0.1579, 0.1620, 0.1661, 0.1677, 0.1743, 0.1768, 0.1793, 0.1826,
    0.1859, 0.1876, 0.1892, 0.1909, 0.2075, 0.2158, 0.0001, 0.0001, 0.0001, 0.0001 };      

//Now plot all the thigs  
//Freon, Quartz, Opaque ,Methane,CsI,Grid   
  const Int_t kC6F14M= 24; const Int_t kC6F14C=kRed;  
  const Int_t kCH4M= 25;   const Int_t kCH4Color=kGreen;  
  const Int_t kSiO2M= 26;  const Int_t kSiO2C=kBlue;  
  const Int_t kCsIMarker = 2;   const Int_t kCsIColor =kMagenta;  
  
  TCanvas *pC=new TCanvas("c1","RICH optics to check",1100,900);

  pC->Divide(3,2);           
           
  pC->cd(1);                 
  TGraph *pIdxC6F14G=new TGraph(kNbins,aPckov,aIdxC6F14);pIdxC6F14G->SetMarkerStyle(kC6F14M); pIdxC6F14G->SetMarkerColor(kC6F14C);
  TGraph *pIdxSiO2G =new TGraph(kNbins,aPckov,aIdxSiO2); pIdxSiO2G->SetMarkerStyle(kSiO2M);   pIdxSiO2G->SetMarkerColor(kSiO2C);  
  TGraph *pIdxCH4G  =new TGraph(kNbins,aPckov,aIdxCH4);  pIdxCH4G->SetMarkerStyle(kCH4M);     pIdxCH4G->SetMarkerColor(kCH4Color);
  TMultiGraph *pIdxMG=new TMultiGraph();  TLegend *pIdxLe=new TLegend(0.6,0.2,0.85,0.4);
  pIdxMG->Add(pIdxC6F14G);     pIdxLe->AddEntry(pIdxC6F14G,   "C6F14","p");            
  pIdxMG->Add(pIdxSiO2G);      pIdxLe->AddEntry(pIdxSiO2G,    "SiO2","p");          
  pIdxMG->Add(pIdxCH4G);       pIdxLe->AddEntry(pIdxCH4G,     "CH4","p");          
  pIdxMG->Draw("APL");pIdxMG->GetXaxis()->SetTitle("ref index versus energy, GeV");pIdxMG->Draw("APL");
  pIdxLe->Draw();      


  pC->cd(2);
  TGraph *pAbsC6F14G=new TGraph(kNbins,aPckov,aAbsC6F14);pAbsC6F14G->SetMarkerStyle(kC6F14M); pAbsC6F14G->SetMarkerColor(kC6F14C);
  TGraph *pAbsSiO2G =new TGraph(kNbins,aPckov,aAbsSiO2); pAbsSiO2G->SetMarkerStyle(kSiO2M);   pAbsSiO2G->SetMarkerColor(kSiO2C);
  TMultiGraph *pAbsMG=new TMultiGraph();  TLegend *pAbsLegend=new TLegend(0.6,0.3,0.85,0.5);
  pAbsMG->Add(pAbsC6F14G);      pAbsLegend->AddEntry(pAbsC6F14G,  "C6F14","p"); 
  pAbsMG->Add(pAbsSiO2G);       pAbsLegend->AddEntry(pAbsSiO2G,  "SiO2", "p"); 
  pAbsMG->Draw("APL"); pAbsMG->GetXaxis()->SetTitle("absorbtion length, cm versus energy, GeV");  pAbsMG->Draw("APL");
  pAbsLegend->Draw();   

  pC->cd(3);      
  TGraph *pAbsCH4Gr=new TGraph(kNbins,aPckov,aAbsCH4);
  pAbsCH4Gr->SetMarkerStyle(kCH4M); pAbsCH4Gr->SetMarkerColor(kCH4Color);
  pAbsCH4Gr->Draw("APL");
  pAbsCH4Gr->GetXaxis()->SetTitle("energy, GeV");
  pAbsCH4Gr->SetTitle("CH4 absorption length, cm");
  pAbsCH4Gr->Draw("APL");
  
  pC->cd(4);
  TGraph *pQeCsIG=new TGraph(kNbins,aPckov,aQeCsI);
  pQeCsIG->SetMarkerStyle(kCsIMarker); pQeCsIG->SetMarkerColor(kCsIColor);
  pQeCsIG->Draw("APL");
  pQeCsIG->GetXaxis()->SetTitle("energy, GeV");
  pQeCsIG->SetTitle("QE");
  pQeCsIG->Draw("APL");
  
//transmission  
  pC->cd(5);
  Float_t aTrC6F14[kNbins],aTrSiO2[kNbins],aTrCH4[kNbins];
  Float_t aTotTr[kNbins];
  for(Int_t i=0;i<kNbins;i++){
    aTrC6F14[i]=TMath::Exp(-15*mm/(aAbsC6F14[i]+0.0001));
    aTrSiO2[i] =TMath::Exp(-5*mm/(aAbsSiO2[i] +0.0001));
    aTrCH4[i]  =TMath::Exp(-8*cm/(aAbsCH4[i]  +0.0001));    
    aTotTr[i]    =aTrC6F14[i]*aTrSiO2[i]*aTrCH4[i]*aQeCsI[i];
  }
  TGraph *pTrC6F14G=new TGraph(kNbins,aPckov,aTrC6F14);pTrC6F14G->SetMarkerStyle(kC6F14M);pTrC6F14G->SetMarkerColor(kC6F14C);  
  TGraph *pTrSiO2G=new TGraph(kNbins,aPckov,aTrSiO2);pTrSiO2G->SetMarkerStyle(kSiO2M); pTrSiO2G->SetMarkerColor(kSiO2C);  
  TGraph *pTrCH4G=new TGraph(kNbins,aPckov,aTrCH4);pTrCH4G->SetMarkerStyle(kCH4M); pTrCH4G->SetMarkerColor(kCH4Color);   
  TGraph *pTotTrG=new TGraph(kNbins,aPckov,aTotTr);pTotTrG->SetMarkerStyle(30);pTotTrG->SetMarkerColor(kYellow);
  TMultiGraph *pTrMG=new TMultiGraph();
  TLegend *pTrLegend=new TLegend(0.2,0.4,0.35,0.6);
  pTrMG->Add(pQeCsIG);       pTrLegend->AddEntry(pQeCsIG,   "CsI QE", "p");  
  pTrMG->Add(pTrC6F14G);     pTrLegend->AddEntry(pTrC6F14G, "C6F14", "p");  
  pTrMG->Add(pTrSiO2G);      pTrLegend->AddEntry(pTrSiO2G,  "SiO2",  "p");          
  pTrMG->Add(pTrCH4G);       pTrLegend->AddEntry(pTrCH4G,   "CH4",   "p");          
  pTrMG->Add(pTotTrG);       pTrLegend->AddEntry(pTotTrG,   "total", "p");          
  pTrMG->Draw("APL");pTrMG->GetXaxis()->SetTitle("transmission versus energy, GeV");pTrMG->Draw("APL");
  pTrLegend->Draw();      

  pC->cd(6);
  
  TMultiGraph *pCompMG=new TMultiGraph;
  pCompMG->Add(pQeCsIG);
  pCompMG->Add(new TGraph(kNbins,aPckov,aQeCsIold));
  pCompMG->Draw("APL");pCompMG->GetXaxis()->SetTitle("comparison of QE versus energy, GeV");pCompMG->Draw("APL");
    
  return;
  TCanvas *pQeC=new TCanvas("pQeC","CsI QE currently all the same",800,900); 
  pQeC->Divide(2,4);
  for(int i=1;i<=7;i++){
    pQeC->cd(i);  
    switch(i){
      case 1: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 1");break;
      case 2: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 2");break;
      case 3: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 3");break;
      case 4: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 4");break;
      case 5: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 5");break;
      case 6: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 6");break;
      case 7: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI);pQeCsIGr->SetTitle("Module 7");break;
    }
    pQeCsIGr->SetMarkerStyle(kCsIMarker); pQeCsIGr->SetMarkerColor(kCsIColor);  
    pQeCsIGr->Draw("APL");
    pQeCsIGr->GetXaxis()->SetTitle("energy, GeV");
    pQeCsIGr->Draw("APL");
  }
}//main

//__________________________________________________________________________________________________
Float_t AbsoCH4(Float_t x)
{//Evaluate the absorbtion lenght of CH4
  Float_t sch4[9] = {.12,.16,.23,.38,.86,2.8,7.9,28.,80.};              //MB X 10^22
  Float_t em[9] = {8.1,8.158,8.212,8.267,8.322,8.378,8.435,8.493,8.55};
  const Float_t kLoschmidt=2.686763e19;                                      // LOSCHMIDT NUMBER IN CM-3
  const Float_t kPressure=750.,kTemperature=283.;                                      
  const Float_t pn=kPressure/760.;
  const Float_t tn=kTemperature/273.16;
  const Float_t c0=-1.655279e-1;
  const Float_t c1=6.307392e-2;
  const Float_t c2=-8.011441e-3;
  const Float_t c3=3.392126e-4;
    		
  Float_t crossSection=0;                        
  if (x<7.75) 
    crossSection=.06e-22;
  else if(x>=7.75 && x<=8.1){                 //------ METHANE CROSS SECTION cm-2 ASTROPH. J. 214, L47 (1978)                                               
	crossSection=(c0+c1*x+c2*x*x+c3*x*x*x)*1.e-18;
  }else if (x> 8.1){
    Int_t j=0;
    while (x<=em[j] || x>=em[j+1]){
      j++;
      Float_t a=(sch4[j+1]-sch4[j])/(em[j+1]-em[j]);
      crossSection=(sch4[j]+a*(x-em[j]))*1e-22;
    }
  }//if
    
    Float_t density=kLoschmidt*pn/tn; //CH4 molecular density 1/cm-3
    return 1./(density*crossSection);
}//AbsoCH4()
//__________________________________________________________________________________________________
Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
{//ENE(EV), PDOTI=COS(INC.ANG.), PDOTR=COS(POL.PLANE ROT.ANG.)
    
    Float_t en[36] = {5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,
		      6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,
		      7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5};
     

    Float_t csin[36] = {2.14,2.21,2.33,2.48,2.76,2.97,2.99,2.59,2.81,3.05,
			2.86,2.53,2.55,2.66,2.79,2.96,3.18,3.05,2.84,2.81,2.38,2.11,
			2.01,2.13,2.39,2.73,3.08,3.15,2.95,2.73,2.56,2.41,2.12,1.95,
			1.72,1.53};
      
    Float_t csik[36] = {0.,0.,0.,0.,0.,0.196,0.408,0.208,0.118,0.49,0.784,0.543,
	 		0.424,0.404,0.371,0.514,0.922,1.102,1.139,1.376,1.461,1.253,0.878,
			0.69,0.612,0.649,0.824,1.347,1.571,1.678,1.763,1.857,1.824,1.824,
			1.714,1.498};
    Float_t xe=ene;
    Int_t  j=Int_t(xe*10)-49;
    Float_t cn=csin[j]+((csin[j+1]-csin[j])/0.1)*(xe-en[j]);
    Float_t ck=csik[j]+((csik[j+1]-csik[j])/0.1)*(xe-en[j]);

    //FORMULAE FROM HANDBOOK OF OPTICS, 33.23 OR
    //W.R. HUNTER, J.O.S.A. 54 (1964),15 , J.O.S.A. 55(1965),1197

    Float_t sinin=TMath::Sqrt(1-pdoti*pdoti);
    Float_t tanin=sinin/pdoti;

    Float_t c1=cn*cn-ck*ck-sinin*sinin;
    Float_t c2=4*cn*cn*ck*ck;
    Float_t aO=TMath::Sqrt(0.5*(TMath::Sqrt(c1*c1+c2)+c1));
    Float_t b2=0.5*(TMath::Sqrt(c1*c1+c2)-c1);
    
    Float_t rs=((aO-pdoti)*(aO-pdoti)+b2)/((aO+pdoti)*(aO+pdoti)+b2);
    Float_t rp=rs*((aO-sinin*tanin)*(aO-sinin*tanin)+b2)/((aO+sinin*tanin)*(aO+sinin*tanin)+b2);
    

    //CORRECTION FACTOR FOR SURFACE ROUGHNESS
    //B.J. STAGG  APPLIED OPTICS, 30(1991),4113

    Float_t sigraf=18.;
    Float_t lamb=1240/ene;
    Float_t fresn;
 
    Float_t  rO=TMath::Exp(-(4*TMath::Pi()*pdoti*sigraf/lamb)*(4*TMath::Pi()*pdoti*sigraf/lamb));

    if(pola)
    {
	Float_t pdotr=0.8;                                 //DEGREE OF POLARIZATION : 1->P , -1->S
	fresn=0.5*(rp*(1+pdotr)+rs*(1-pdotr));
    }
    else
	fresn=0.5*(rp+rs);
      
    fresn = fresn*rO;
    return(fresn);
}//Fresnel()
#endif
