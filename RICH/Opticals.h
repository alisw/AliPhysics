#ifdef __CINT__
void Opticals()
{
  gROOT->Reset();
#endif   
  int i;
  const Int_t kNbins=30;//number of photon energy points
  
  Float_t aPckov[kNbins];  for(i=0;i<kNbins;i++) aPckov[i]=1e-9*AliRICHParam::PhotonEnergy(i);
    
  
  Float_t aIdxC6F14[kNbins];  //Freon ref index
  Float_t aIdxCH4[kNbins];    //Methane ref index
  Float_t aIdxSiO2[kNbins];   //Quartz ref index
  Float_t aIdxMetal[kNbins];  //Metal ref index always 0
  Float_t aIdxGel[kNbins];    //Aerogel ref index
        
  for (i=0;i<kNbins;i++){
    aIdxC6F14[i] = AliRICHParam::RefIdxC6F14(aPckov[i]*1e9);
    aIdxSiO2[i]  = AliRICHParam::RefIdxSiO2(aPckov[i]*1e9);
    aIdxCH4[i]   = AliRICHParam::RefIdxCH4(5.5);
    aIdxGel[i]   = AliRICHParam::RefIdxGel(5.5);
    aIdxMetal[i]  =0;//metal
  } 
    
  Float_t aAbsC6F14[kNbins]={//New values from A.DiMauro 28.10.03 total 30
    32701.4219, 17996.1141, 10039.7281, 1799.1230, 1799.1231, 1799.1231, 1241.4091, 179.0987, 179.0986, 179.0987,
      179.0987,   118.9800,    39.5058,   23.7244,   11.1283,    7.1573,    3.6249,   2.1236,   0.7362,   0.5348,
        0.3387,     0.3074,     0.3050,    0.0001,    0.0001,    0.0001,    0.0001,   0.0001,   0.0001,   0.0001};    
  Float_t aAbsSiO2[kNbins]={//New values from A.DiMauro 28.10.03 total 30
       34.4338, 30.5424, 30.2584, 31.4928, 31.7868, 17.8397, 9.3410, 6.4492, 6.1128, 5.8128,
        5.5589,  5.2877,  5.0162,  4.7999,  4.5734,  4.2135, 3.7471, 2.6033, 1.5223, 0.9658,
        0.4242,  0.2500,  0.1426,  0.0863,  0.0793,  0.0724, 0.0655, 0.0587, 0.0001, 0.0001};
    
  Float_t aAbsOpSiO2[kNbins];
  Float_t aAbsCH4[kNbins];
  Float_t aAbsGrid[kNbins];
  Float_t aAbsCsI[kNbins];
  Float_t aAbsGel[kNbins];
  Float_t aAbsRef[kNbins];//abs prob on reflector
  
    for (i=0;i<kNbins;i++){
      aAbsCH4[i]         =AliRICHParam::AbsCH4(AliRICHParam::PhotonEnergy(i)); 
      aAbsOpSiO2[i]      =1e-5; 
      aAbsCsI[i]         =1e-4; 
      aAbsGrid[i]        =1e-4;       
      aAbsGel[i]         =AliRICHParam::AbsGel(AliRICHParam::PhotonEnergy(i)); 
      aAbsRef[i]         =0.01;//1% reflector abs prob 
    }
    
  Float_t aQeCsI[kNbins] = {//New values from A.DiMauro 28.10.03 total 31
                            0.0002, 0.0006, 0.0007, 0.0010, 0.0049, 0.0073, 0.0104, 0.0519, 0.0936, 0.1299,
                            0.1560, 0.1768, 0.1872, 0.1976, 0.2142, 0.2288, 0.2434, 0.2599, 0.2673, 0.2808,
                            0.2859, 0.2954, 0.3016, 0.3120, 0.3172, 0.3224, 0.3266, 0.3328, 0.3359, 0.3390};
//                             0.3431};

  
  Float_t aQeAll[kNbins];
  for(i=0;i<kNbins;i++){
    aQeCsI[i]/= (1.0-Fresnel(AliRICHParam::PhotonEnergy(i),1.0,0)); //FRESNEL LOSS CORRECTION
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
//Markers and colors for different curves:
    
  const Double_t kWidth=0.25,kHeight=0.2;  
  const Int_t kC6F14Marker=24 , kC6F14Color=kRed;  
  const Int_t kCH4Marker  =25 , kCH4Color  =kGreen;  
  const Int_t kSiO2M      =26 , kSiO2Color =kBlue;  
  const Int_t kCsIMarker  = 2 , kCsIColor  =kMagenta;  
  const Int_t kGelMarker  =27 , kGelColor  =46;  
  const Int_t kRefMarker  =28 , kRefColor  =47;  
  
  TCanvas *pC=new TCanvas("c1","RICH optics to check",1100,900);

  pC->Divide(3,2);           
//Ref index           
  pC->cd(1);                 
  TGraph *pIdxC6F14=new TGraph(kNbins,aPckov,aIdxC6F14);pIdxC6F14->SetMarkerStyle(kC6F14Marker);pIdxC6F14->SetMarkerColor(kC6F14Color);
  TGraph *pIdxSiO2 =new TGraph(kNbins,aPckov,aIdxSiO2); pIdxSiO2 ->SetMarkerStyle(kSiO2M)      ;pIdxSiO2 ->SetMarkerColor(kSiO2Color);  
  TGraph *pIdxCH4  =new TGraph(kNbins,aPckov,aIdxCH4);  pIdxCH4  ->SetMarkerStyle(kCH4Marker)       ;pIdxCH4  ->SetMarkerColor(kCH4Color);
  TGraph *pIdxGel  =new TGraph(kNbins,aPckov,aIdxGel);  pIdxGel  ->SetMarkerStyle(kGelMarker)  ;pIdxGel  ->SetMarkerColor(kGelColor);
  
  TMultiGraph *pIdxMG=new TMultiGraph();  TLegend *pIdxLe=new TLegend(0.5,0.21,0.5+kWidth,0.21+kHeight);
  pIdxMG->Add(pIdxC6F14); pIdxLe->AddEntry(pIdxC6F14,"C6F14"  ,"p");            
  pIdxMG->Add(pIdxSiO2) ; pIdxLe->AddEntry(pIdxSiO2 ,"SiO2"   ,"p");          
  pIdxMG->Add(pIdxCH4)  ; pIdxLe->AddEntry(pIdxCH4  ,"CH4"    ,"p");          
  pIdxMG->Add(pIdxGel)  ; pIdxLe->AddEntry(pIdxGel  ,"Aerogel","p");          
  pIdxMG->Draw("APL");
  pIdxMG->SetTitle("Refractive index");  pIdxMG->GetXaxis()->SetTitle("energy, GeV");
  pIdxMG->Draw("APL");
  pIdxLe->Draw();      
//Absorbtion
  pC->cd(2);
  gPad->SetLogy();
  TGraph *pAbsC6F14=new TGraph(kNbins,aPckov,aAbsC6F14);pAbsC6F14->SetMarkerStyle(kC6F14Marker); pAbsC6F14->SetMarkerColor(kC6F14Color);
  TGraph *pAbsSiO2 =new TGraph(kNbins,aPckov,aAbsSiO2) ;pAbsSiO2 ->SetMarkerStyle(kSiO2M)      ; pAbsSiO2 ->SetMarkerColor(kSiO2Color);
  TGraph *pAbsCH4  =new TGraph(kNbins,aPckov,aAbsCH4)  ;pAbsCH4  ->SetMarkerStyle(kCH4Marker)  ; pAbsCH4  ->SetMarkerColor(kCH4Color);
  TGraph *pAbsGel  =new TGraph(kNbins,aPckov,aAbsGel)  ;pAbsGel  ->SetMarkerStyle(kGelMarker)  ; pAbsGel  ->SetMarkerColor(kGelColor);
  TGraph *pAbsRef  =new TGraph(kNbins,aPckov,aAbsRef)  ;pAbsRef  ->SetMarkerStyle(kRefMarker)  ; pAbsRef  ->SetMarkerColor(kRefColor);
  
  TMultiGraph *pAbsMG=new TMultiGraph();  TLegend *pAbsLe=new TLegend(0.2,0.15,0.2+kWidth,0.15+kHeight);
  pAbsMG->Add(pAbsC6F14);      pAbsLe->AddEntry(pAbsC6F14,  "C6F14"    ,"p"); 
  pAbsMG->Add(pAbsSiO2);       pAbsLe->AddEntry(pAbsSiO2 ,  "SiO2"     ,"p"); 
  pAbsMG->Add(pAbsCH4);        pAbsLe->AddEntry(pAbsCH4  ,  "CH4"      ,"p"); 
  pAbsMG->Add(pAbsGel);        pAbsLe->AddEntry(pAbsGel  ,  "Aerogel"  ,"p"); 
  pAbsMG->Add(pAbsRef);        pAbsLe->AddEntry(pAbsRef  ,  "Reflector","p"); 
  pAbsMG->Draw("APL"); 
  pAbsMG->SetTitle("Absorbtion length,cm"); pAbsMG->GetXaxis()->SetTitle("energy, GeV");  
  pAbsMG->Draw("APL");
  pAbsLe->Draw();   
//QE  
  pC->cd(4);
  TGraph *pQeCsI=new TGraph(kNbins,aPckov,aQeCsI); pQeCsI->SetMarkerStyle(kCsIMarker); pQeCsI->SetMarkerColor(kCsIColor);
  pQeCsI->Draw("APL");
  pQeCsI->SetTitle("QE");  pQeCsI->GetXaxis()->SetTitle("energy, GeV");
  pQeCsI->Draw("APL");  
//transmission  
  pC->cd(5);
  Float_t mm =0.1;Float_t cm=1.; 
  Float_t aTrC6F14[kNbins],aTrSiO2[kNbins],aTrCH4[kNbins];
  Float_t aTotTr[kNbins];
  for(Int_t i=0;i<kNbins;i++){//calculate probability for photon to survive during transversing a volume of material with absorption length  
    aTrC6F14[i]=TMath::Exp(-15*mm/(aAbsC6F14[i]+0.0001));
    aTrSiO2[i] =TMath::Exp(-5*mm/(aAbsSiO2[i] +0.0001));
    aTrCH4[i]  =TMath::Exp(-8*cm/(aAbsCH4[i]  +0.0001));    
    aTotTr[i]    =aTrC6F14[i]*aTrSiO2[i]*aTrCH4[i]*aQeCsI[i];
  }
  TGraph *pTrC6F14=new TGraph(kNbins,aPckov,aTrC6F14);pTrC6F14->SetMarkerStyle(kC6F14Marker);pTrC6F14->SetMarkerColor(kC6F14Color);  
  TGraph *pTrSiO2 =new TGraph(kNbins,aPckov,aTrSiO2) ;pTrSiO2 ->SetMarkerStyle(kSiO2M)      ;pTrSiO2 ->SetMarkerColor(kSiO2Color);  
  TGraph *pTrCH4  =new TGraph(kNbins,aPckov,aTrCH4)  ;pTrCH4  ->SetMarkerStyle(kCH4Marker)  ;pTrCH4  ->SetMarkerColor(kCH4Color);   
  TGraph *pTotTr  =new TGraph(kNbins,aPckov,aTotTr)  ;pTotTr  ->SetMarkerStyle(30)          ;pTotTr  ->SetMarkerColor(kYellow);
  
  TMultiGraph *pTrMG=new TMultiGraph();  TLegend *pTrLe=new TLegend(0.2,0.4,0.2+kWidth,0.4+kHeight);
  pTrMG->Add(pQeCsI);       pTrLe->AddEntry(pQeCsI,   "CsI QE", "p");  
  pTrMG->Add(pTrC6F14);     pTrLe->AddEntry(pTrC6F14, "C6F14", "p");  
  pTrMG->Add(pTrSiO2);      pTrLe->AddEntry(pTrSiO2,  "SiO2",  "p");          
  pTrMG->Add(pTrCH4);       pTrLe->AddEntry(pTrCH4,   "CH4",   "p");          
  pTrMG->Add(pTotTr);       pTrLe->AddEntry(pTotTr,   "total", "p");          
  pTrMG->Draw("APL");
  pTrMG->SetTitle("Transmission");  pTrMG->GetXaxis()->SetTitle("energy, GeV");
  pTrMG->Draw("APL");
  pTrLe->Draw();      
//comparison of new and old staff
  pC->cd(6);
  TGraph *pQeCsIold=new TGraph(kNbins,aPckov,aQeCsIold); pQeCsIold->SetMarkerStyle(kC6F14Marker);pQeCsIold->SetMarkerColor(kC6F14Color);  
  
  TMultiGraph *pCompMG=new TMultiGraph;  TLegend *pCompLe=new TLegend(0.2,0.6,0.2+kWidth,0.6+kHeight);
  pCompMG->Add(pQeCsI);       pCompLe->AddEntry(pQeCsI,    "QE new 30.10.03", "p");  
  pCompMG->Add(pQeCsIold);    pCompLe->AddEntry(pQeCsIold, "QE old 01.01.02", "p");
  pCompMG->Draw("APL");
  pCompMG->SetTitle("Comparison of new and old staff");
  pCompMG->GetXaxis()->SetTitle("energy, GeV");
  pCompMG->Draw("APL");
  pCompLe->Draw();      
    
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
