#ifdef __CINT__
void Opticals()
{
  gROOT->Reset();
#endif   
  int i;
  const Int_t kNbins=30;
  Float_t aPckov[kNbins];
  for(i=0;i<kNbins;i++){    //Photons energy bins 5.5 eV - 8.5 eV step 0.1 eV
    aPckov[i]=(0.1*i+5.5)*1e-9;
  }
  
  Float_t aIndexFreon[kNbins];  
  Float_t aIndexSiO2[kNbins];
  Float_t aIndexOpaqueQuartz[kNbins];
  Float_t aIndexCH4[kNbins];
  Float_t aIndexGrid[kNbins];
        
  Float_t  e1= 10.666;Float_t  e2= 18.125;  Float_t  f1= 46.411; Float_t  f2= 228.71;//RICH TDR page 35 
  for (i=0;i<kNbins;i++){
    aIndexFreon[i]        = aPckov[i] * .0172 * 1e9 + 1.177;
    Float_t ene=aPckov[i]*1e9;
    aIndexSiO2[i]        = TMath::Sqrt(1. + f1/(e1*e1 - ene*ene) + f2/(e2*e2 - ene*ene) );
    aIndexOpaqueQuartz[i]  =1;
    aIndexCH4[i]      =1.000444;
    aIndexGrid[i]         =1;
  } 
    
  Float_t aAbsFreon[kNbins]={179.0987,  179.0987,    179.0987,  179.0987,   179.0987,   //Previous values from code
                             179.0987,  179.0987,    179.0987,  179.0987,   142.9206, 
                              56.64957,  25.58622,    13.95293,  12.03905,   10.42953, 
                               8.804196,  7.069031,    4.461292,  2.028366,   1.293013, 
                               0.577267,  0.40746,     0.334964,  0.00001,        0.00001,       
                               0.00001};
    
//   Float_t aAbsFreon[kNbins]={32701.4219, 17996.1141, 10039.7281, 1799.1231, 1799.1231,//New values from A.DiMauro 28.10.03
//                               1799.1231,  1241.4091,   179.0987,  179.0987,  179.0987,  179.0987,   179.0987,  
//                              179.0987,  179.0987,    179.0987,  179.0987,   142.9206, 
//                               56.64957,  25.58622,    13.95293,  12.03905,   10.42953, 
//                                8.804196,  7.069031,    4.461292,  2.028366,   1.293013, 
//                                0.577267,  0.40746,     0.334964,  0.0,        0.0,       
//                                0.0000, 0.0000, 0.0000, 0.0000, 0.0000};

  Float_t aAbsSiO2[kNbins]={105.8,   65.52,  48.58,  42.85,   35.79,
                               31.262, 28.598, 27.527, 25.007,  22.815, 
                               21.004, 19.266, 17.525, 15.878,  14.177, 
                               11.719, 9.282,   6.62,   4.0925,  2.601, 
                               1.149,  0.667,   0.3627, 0.192,   0.1497, 
                               0.10857};
    
  Float_t aAbsOpaqueQuartz[kNbins];
  Float_t aAbsCH4[kNbins];
  Float_t aAbsGrid[kNbins];
  Float_t aAbsCsI[kNbins];
    for (i=0;i<kNbins;i++){
      aAbsCH4[i]    =AbsoCH4(aPckov[i]*1e9); 
      aAbsOpaqueQuartz[i]=1e-5; 
      aAbsCsI[i]        =1e-4; 
      aAbsGrid[i]       =1e-4; 
    }
    
  Float_t aQeCsI1[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  Float_t aQeCsI2[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  Float_t aQeCsI3[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  Float_t aQeCsI4[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  Float_t aQeCsI5[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  Float_t aQeCsI6[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  Float_t aQeCsI7[kNbins] = {0.0002, 0.0003, 0.0005, 0.0010, 0.0020, 0.0020, 0.0099, 0.0500, 0.1800, 0.1850, 
                             0.1899, 0.1920, 0.1979, 0.2099, 0.2179, 0.2249, 0.2349, 0.2599, 0.2800, 0.3000, 
                             0.3100, 0.3249, 0.3400, 0.3549, 0.3600, 0.3700, 0.3799, 0.3899, 0.4000, 0.4099};
  
  Float_t aQeAll[kNbins];
  for(i=0;i<kNbins;i++){
    aQeCsI1[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeCsI2[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeCsI3[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeCsI4[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeCsI5[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeCsI6[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeCsI7[i]/= (1.0-Fresnel(aPckov[i]*1e9,1.0,0)); //FRESNEL LOSS CORRECTION
    aQeAll[i]=1; //QE for all other materials except for PC must be 1.
  }
       
                    
#ifdef __CINT__

//Now plot all the thigs  
//Freon, Quartz, Opaque ,Methane,CsI,Grid   
  const Int_t kFreonMarker= 24; const Int_t kFreonColor=kRed;  
  const Int_t kCH4Marker= 25;   const Int_t kCH4Color=kGreen;  
  const Int_t kSiO2Marker= 26;  const Int_t kSiO2Color=kBlue;  
  const Int_t kCsIMarker = 2;   const Int_t kCsIColor =kMagenta;  
  
  TCanvas *pC=new TCanvas("c1","RICH optics to check",1100,900);

  pC->Divide(3,2);           
           
  pC->cd(1);                 
  TGraph *pAbsFreonGr=new TGraph(kNbins,aPckov,aAbsFreon);
  pAbsFreonGr->SetMarkerStyle(kFreonMarker); pAbsFreonGr->SetMarkerColor(kFreonColor);
  TGraph *pAbsSiO2Gr=new TGraph(kNbins,aPckov,aAbsSiO2);
  pAbsSiO2Gr->SetMarkerStyle(kSiO2Marker); pAbsSiO2Gr->SetMarkerColor(kSiO2Color);  
  TMultiGraph *pAbsMG=new TMultiGraph();
  TLegend *pAbsLegend=new TLegend(0.6,0.3,0.85,0.5);
  pAbsMG->Add(pAbsFreonGr);     pAbsLegend->AddEntry(pAbsFreonGr,   "freon","p");            //1
  pAbsMG->Add(pAbsSiO2Gr);      pAbsLegend->AddEntry(pAbsSiO2Gr,  "SiO2","p");          //2
  pAbsMG->Draw("APL");                       		
  pAbsMG->GetXaxis()->SetTitle("energy, GeV");
  pAbsMG->GetYaxis()->SetTitle("absorption length, cm");
  pAbsMG->Draw("APL");
  pAbsLegend->Draw();   


  pC->cd(2);
  TGraph *pIndexFreonGr=new TGraph(kNbins,aPckov,aIndexFreon);
  pIndexFreonGr->SetMarkerStyle(kFreonMarker); pIndexFreonGr->SetMarkerColor(kFreonColor);  
  TGraph *pIndexSiO2Gr=new TGraph(kNbins,aPckov,aIndexSiO2);
  pIndexSiO2Gr->SetMarkerStyle(kSiO2Marker); pIndexSiO2Gr->SetMarkerColor(kSiO2Color);  
  TGraph *pIndexCH4Gr=new TGraph(kNbins,aPckov,aIndexCH4);
  pIndexCH4Gr->SetMarkerStyle(kCH4Marker); pIndexCH4Gr->SetMarkerColor(kCH4Color);   
  TMultiGraph *pIndexMG=new TMultiGraph();
  TLegend *pIndexLegend=new TLegend(0.6,0.2,0.85,0.4);
  pIndexMG->Add(pIndexFreonGr);     pIndexLegend->AddEntry(pIndexFreonGr,   "C6F14","p");            //1
  pIndexMG->Add(pIndexSiO2Gr);      pIndexLegend->AddEntry(pIndexSiO2Gr,  "SiO2","p");          
  pIndexMG->Add(pIndexCH4Gr);       pIndexLegend->AddEntry(pIndexCH4Gr,   "CH4","p");          
  pIndexMG->Draw("APL");                       		
  pIndexMG->GetXaxis()->SetTitle("energy, GeV");
  pIndexMG->GetYaxis()->SetTitle("refraction index");
  pIndexMG->Draw("APL");
  pIndexLegend->Draw();      

  pC->cd(3);      
  TGraph *pAbsCH4Gr=new TGraph(kNbins,aPckov,aAbsCH4);
  pAbsCH4Gr->SetMarkerStyle(kCH4Marker); pAbsCH4Gr->SetMarkerColor(kCH4Color);
  pAbsCH4Gr->Draw("APL");
  pAbsCH4Gr->GetXaxis()->SetTitle("energy, GeV");
  pAbsCH4Gr->SetTitle("CH4 absorption length, cm");
  pAbsCH4Gr->Draw("APL");
  
  pC->cd(4);
  TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI1);
  pQeCsIGr->SetMarkerStyle(kCsIMarker); pQeCsIGr->SetMarkerColor(kCsIColor);
  pQeCsIGr->Draw("APL");
  pQeCsIGr->GetXaxis()->SetTitle("energy, GeV");
  pQeCsIGr->SetTitle("QE");
  pQeCsIGr->Draw("APL");
  
  
  pC->cd(5);
  Float_t aTrFreon[kNbins],aTrSiO2[kNbins],aTrCH4[kNbins];
  Float_t aTotTr[kNbins];
  for(Int_t i=0;i<kNbins;i++){
    aTrFreon[i]=TMath::Exp(-AliRICHParam::FreonThickness() /(aAbsFreon[i]+0.0001));
    aTrSiO2[i] =TMath::Exp(-AliRICHParam::QuartzThickness()/(aAbsSiO2[i] +0.0001));
    aTrCH4[i]  =TMath::Exp(-AliRICHParam::GapThickness()   /(aAbsCH4[i]  +0.0001));    
    aTotTr[i]    =aTrFreon[i]*aTrSiO2[i]*aTrCH4[i]*aQeCsI1[i];
  }
  TGraph *pTrFreonG=new TGraph(kNbins,aPckov,aTrFreon);pTrFreonG->SetMarkerStyle(kFreonMarker);pTrFreonG->SetMarkerColor(kFreonColor);  
  TGraph *pTrSiO2G=new TGraph(kNbins,aPckov,aTrSiO2);pTrSiO2G->SetMarkerStyle(kSiO2Marker); pTrSiO2G->SetMarkerColor(kSiO2Color);  
  TGraph *pTrCH4G=new TGraph(kNbins,aPckov,aTrCH4);pTrCH4G->SetMarkerStyle(kCH4Marker); pTrCH4G->SetMarkerColor(kCH4Color);   
  TGraph *pTotTrG=new TGraph(kNbins,aPckov,aTotTr);pTotTrG->SetMarkerStyle(30);
  TMultiGraph *pTrMG=new TMultiGraph();
  TLegend *pTrLegend=new TLegend(0.6,0.2,0.85,0.4);
  pTrMG->Add(pTrFreonG);     pTrLegend->AddEntry(pTrFreonG, "freon","p");            //1
  pTrMG->Add(pTrSiO2G);      pTrLegend->AddEntry(pTrSiO2G,  "SiO2","p");          
  pTrMG->Add(pTrCH4G);       pTrLegend->AddEntry(pTrCH4G,   "CH4","p");          
  pTrMG->Add(pTotTrG);       pTrLegend->AddEntry(pTotTrG,   "total","p");          
  pTrMG->Draw("APL");                       		
  pTrMG->GetXaxis()->SetTitle("energy, GeV");
  pTrMG->GetYaxis()->SetTitle("transmission");
  pTrMG->Draw("APL");
  pTrLegend->Draw();      

  
  return;
  TCanvas *pQeC=new TCanvas("pQeC","CsI QE currently all the same",800,900); 
  pQeC->Divide(2,4);
  for(int i=1;i<=7;i++){
    pQeC->cd(i);  
    switch(i){
      case 1: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI1);pQeCsIGr->SetTitle("Module 1");break;
      case 2: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI2);pQeCsIGr->SetTitle("Module 2");break;
      case 3: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI3);pQeCsIGr->SetTitle("Module 3");break;
      case 4: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI4);pQeCsIGr->SetTitle("Module 4");break;
      case 5: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI5);pQeCsIGr->SetTitle("Module 5");break;
      case 6: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI6);pQeCsIGr->SetTitle("Module 6");break;
      case 7: TGraph *pQeCsIGr=new TGraph(kNbins,aPckov,aQeCsI7);pQeCsIGr->SetTitle("Module 7");break;
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
}//Fresnel(...)
#endif
