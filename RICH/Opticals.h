#ifdef __CINT__
void Opticals()
{
  gROOT->Reset();
#endif   
  int i;
  const Int_t kNbins=26;
  Float_t aPckov[kNbins];
  for(i=0;i<kNbins;i++){    //Photons energy intervals
    aPckov[i]=(0.1*i+5.5)*1e-9;
  }
  
  Float_t aIndexFreon[kNbins];  
  Float_t aIndexQuartz[kNbins];
  Float_t aIndexOpaqueQuartz[kNbins];
  Float_t aIndexCH4[kNbins];
  Float_t aIndexGrid[kNbins];
        
  Float_t  e1= 10.666;Float_t  e2= 18.125;  Float_t  f1= 46.411; Float_t  f2= 228.71;//RICH TDR page 35 
  for (i=0;i<kNbins;i++){
    aIndexFreon[i]        = aPckov[i] * .0172 * 1e9 + 1.177;
    Float_t ene=aPckov[i]*1e9;
    aIndexQuartz[i]        = TMath::Sqrt(1. + f1/(e1*e1 - ene*ene) + f2/(e2*e2 - ene*ene) );
    aIndexOpaqueQuartz[i]  =1;
    aIndexCH4[i]      =1.000444;
    aIndexGrid[i]         =1;
  } 
    
  Float_t aAbsFreon[kNbins]={179.0987,  179.0987,    179.0987,  179.0987,   179.0987,  
                             179.0987,  179.0987,    179.0987,  179.0987,   142.9206, 
                              56.64957,  25.58622,    13.95293,  12.03905,   10.42953, 
                               8.804196,  7.069031,    4.461292,  2.028366,   1.293013, 
                               0.577267,  0.40746,     0.334964,  0.0,        0.0,       
                               0};
    

  Float_t aAbsQuartz[kNbins]={105.8,   65.52,  48.58,  42.85,   35.79,
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
    
  Float_t aQeCsI1[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
  Float_t aQeCsI2[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
  Float_t aQeCsI3[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
  Float_t aQeCsI4[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
  Float_t aQeCsI5[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
  Float_t aQeCsI6[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
  Float_t aQeCsI7[kNbins] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983,
                            0.010125,	    0.0242999997,   0.0405000001,   0.0688500032,  0.105299994, 
                            0.121500008,    0.141749993,    0.157949999,    0.162,         0.166050002, 
                            0.167669997,    0.174299985,    0.176789999,    0.179279998,   0.182599992, 
                            0.18592,        0.187579989,    0.189239994,    0.190899998,   0.207499996, 
                            0.215799987};
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

  
//Freon, Quartz, Opaque Quartz,Methane,CsI,Grid   
  const Int_t kFreonMarker= 24; const Int_t kFreonColor=kRed;  
  const Int_t kCH4Marker= 25;   const Int_t kCH4Color=kGreen;  
  const Int_t kSiO2Marker= 26;  const Int_t kSiO2Color=kBlue;  
  const Int_t kCsIMarker = 2;   const Int_t kCsIColor =kMagenta;  
  
  TCanvas *pC=new TCanvas("c1","RICH optics to check",800,900);

  pC->Divide(2,2);           
           
  pC->cd(1);                 
  TGraph *pAbsFreonGr=new TGraph(kNbins,aPckov,aAbsFreon);
  pAbsFreonGr->SetMarkerStyle(kFreonMarker); pAbsFreonGr->SetMarkerColor(kFreonColor);
  TGraph *pAbsSiO2Gr=new TGraph(kNbins,aPckov,aAbsQuartz);
  pAbsSiO2Gr->SetMarkerStyle(kSiO2Marker); pAbsSiO2Gr->SetMarkerColor(kSiO2Color);  
  TMultiGraph *pAbsMG=new TMultiGraph();
  TLegend *pAbsLegend=new TLegend(0.6,0.3,0.85,0.5);
  pAbsMG->Add(pAbsFreonGr);     pAbsLegend->AddEntry(pAbsFreonGr,   "freon","p");            //1
  pAbsMG->Add(pAbsSiO2Gr);      pAbsLegend->AddEntry(pAbsSiO2Gr,  "quartz","p");          //2
  pAbsMG->Draw("APL");                       		
  pAbsMG->GetXaxis()->SetTitle("energy, GeV");
  pAbsMG->GetYaxis()->SetTitle("absorption length, cm");
  pAbsMG->Draw("APL");
  pAbsLegend->Draw();   


  pC->cd(2);
  TGraph *pIndexFreonGr=new TGraph(kNbins,aPckov,aIndexFreon);
  pIndexFreonGr->SetMarkerStyle(kFreonMarker); pIndexFreonGr->SetMarkerColor(kFreonColor);  
  TGraph *pIndexSiO2Gr=new TGraph(kNbins,aPckov,aIndexQuartz);
  pIndexSiO2Gr->SetMarkerStyle(kSiO2Marker); pIndexSiO2Gr->SetMarkerColor(kSiO2Color);  
  TGraph *pIndexCH4Gr=new TGraph(kNbins,aPckov,aIndexCH4);
  pIndexCH4Gr->SetMarkerStyle(kCH4Marker); pIndexCH4Gr->SetMarkerColor(kCH4Color);   
  TMultiGraph *pIndexMG=new TMultiGraph();
  TLegend *pIndexLegend=new TLegend(0.6,0.2,0.85,0.4);
  pIndexMG->Add(pIndexFreonGr);     pIndexLegend->AddEntry(pIndexFreonGr,   "freon","p");            //1
  pIndexMG->Add(pIndexSiO2Gr);      pIndexLegend->AddEntry(pIndexSiO2Gr,  "quartz","p");          
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

Float_t AbsoCH4(Float_t x)
{

    //KLOSCH,SCH4(9),WL(9),EM(9),ALENGTH(31)
    Float_t sch4[9] = {.12,.16,.23,.38,.86,2.8,7.9,28.,80.};              //MB X 10^22
    //Float_t wl[9] = {153.,152.,151.,150.,149.,148.,147.,146.,145};
    Float_t em[9] = {8.1,8.158,8.212,8.267,8.322,8.378,8.435,8.493,8.55};
    const Float_t kLosch=2.686763E19;                                      // LOSCHMIDT NUMBER IN CM-3
    const Float_t kIgas1=100, kIgas2=0, kOxy=10., kWater=5., kPressure=750.,kTemperature=283.;  
    Float_t pn=kPressure/760.;
    Float_t tn=kTemperature/273.16;
    
	
// ------- METHANE CROSS SECTION -----------------
// ASTROPH. J. 214, L47 (1978)
	
    Float_t sm=0;
    if (x<7.75) 
	sm=.06e-22;
    
    if(x>=7.75 && x<=8.1)
    {
	Float_t c0=-1.655279e-1;
	Float_t c1=6.307392e-2;
	Float_t c2=-8.011441e-3;
	Float_t c3=3.392126e-4;
	sm=(c0+c1*x+c2*x*x+c3*x*x*x)*1.e-18;
    }
    
    if (x> 8.1)
    {
	Int_t j=0;
	while (x<=em[j] && x>=em[j+1])
	{
	    j++;
	    Float_t a=(sch4[j+1]-sch4[j])/(em[j+1]-em[j]);
	    sm=(sch4[j]+a*(x-em[j]))*1e-22;
	}
    }
    
    Float_t dm=(kIgas1/100.)*(1.-((kOxy+kWater)/1.e6))*kLosch*pn/tn;
    Float_t abslm=1./sm/dm;
    
//    ------- ISOBUTHANE CROSS SECTION --------------
//     i-C4H10 (ai) abs. length from curves in
//     Lu-McDonald paper for BARI RICH workshop .
//     -----------------------------------------------------------
    
    Float_t ai;
    Float_t absli;
    if (kIgas2 != 0) 
    {
	if (x<7.25)
	    ai=100000000.;
	
	if(x>=7.25 && x<7.375)
	    ai=24.3;
	
	if(x>=7.375)
	    ai=.0000000001;
	
	Float_t si = 1./(ai*kLosch*273.16/293.);                    // ISOB. CRO.SEC.IN CM2
	Float_t di=(kIgas2/100.)*(1.-((kOxy+kWater)/1.e6))*kLosch*pn/tn;
	absli =1./si/di;
    }
    else
	absli=1.e18;
//    ---------------------------------------------------------
//
//       transmission of O2
//
//       y= path in cm, x=energy in eV
//       so= cross section for UV absorption in cm2
//       do= O2 molecular density in cm-3
//    ---------------------------------------------------------
    
    Float_t abslo;
    Float_t so=0;
    if(x>=6.0)
    {
	if(x>=6.0 && x<6.5)
	{
	    so=3.392709e-13 * TMath::Exp(2.864104 *x);
	    so=so*1e-18;
	}
	
	if(x>=6.5 && x<7.0) 
	{
	    so=2.910039e-34 * TMath::Exp(10.3337*x);
	    so=so*1e-18;
	}
	    

	if (x>=7.0) 
	{
	    Float_t a0=-73770.76;
	    Float_t a1=46190.69;
	    Float_t a2=-11475.44;
	    Float_t a3=1412.611;
	    Float_t a4=-86.07027;
	    Float_t a5=2.074234;
	    so= a0+(a1*x)+(a2*x*x)+(a3*x*x*x)+(a4*x*x*x*x)+(a5*x*x*x*x*x);
	    so=so*1e-18;
	}
	
	Float_t dox=(kOxy/1e6)*kLosch*pn/tn;
	abslo=1./so/dox;
    }
    else
	abslo=1.e18;
//     ---------------------------------------------------------
//
//       transmission of H2O
//
//       y= path in cm, x=energy in eV
//       sw= cross section for UV absorption in cm2
//       dw= H2O molecular density in cm-3
//     ---------------------------------------------------------
    
    Float_t abslw;
    
    Float_t b0=29231.65;
    Float_t b1=-15807.74;
    Float_t b2=3192.926;
    Float_t b3=-285.4809;
    Float_t b4=9.533944;
    
    if(x>6.75)
    {    
	Float_t sw= b0+(b1*x)+(b2*x*x)+(b3*x*x*x)+(b4*x*x*x*x);
	sw=sw*1e-18;
	Float_t dw=(kWater/1e6)*kLosch*pn/tn;
	abslw=1./sw/dw;
    }
    else
    	abslw=1.e18;
	    
//    ---------------------------------------------------------
    
    Float_t alength=1./(1./abslm+1./absli+1./abslo+1./abslw);
    return (alength);
}//AbsoCH4


Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
{

    //ENE(EV), PDOTI=COS(INC.ANG.), PDOTR=COS(POL.PLANE ROT.ANG.)
    
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
