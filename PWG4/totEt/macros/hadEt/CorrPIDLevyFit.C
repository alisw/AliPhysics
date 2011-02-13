//Christine Nattrass, University of Tennessee at Knoxville
//This macro is to calculate the contributions to Et from various particles based on Levy fits to ALICE data at 900 GeV
//It uses d^2N/dpTdy from the papers, weighs it by Et, and integrates over all pT. 
//A=0 for mesons
//A=1 for antinucleons since they would anihilate in the calorimeter if they were observed by a calorimeter
//A=-1 for nucleons since they would not anihilate and their mass would not be measured
//At the end this is used to calculate the pid calculation correction
class AliAnalysisLevyPtModified{
public:
  virtual ~AliAnalysisLevyPtModified(){;};
  

  Double_t Evaluate(Double_t *pt, Double_t *par)
  {
    Double_t ldNdy  = par[0];
    Double_t l2pi   = 2*TMath::Pi();
    Double_t lTemp = par[1];
    Double_t lPower = par[2];
    Double_t lMass  = par[3];
    Double_t A = par[4];
    Double_t lMassEt  = par[5];//only used for calculating Et
    //this is the Et we would calculate if we had identified the particle as having a mass lMassEt
    Double_t Et = TMath::Sqrt(pt[0]*pt[0]+lMassEt*lMassEt)+A*lMassEt;
    
    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);
    //This is the density of particles times the Et weight
    return ldNdy *Et* pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
  }
  ClassDef(AliAnalysisLevyPtModified, 1);
};
class AliAnalysisLevyPtModifiedStrangeness{
public:
  virtual ~AliAnalysisLevyPtModifiedStrangeness(){;};
  

  Double_t Evaluate(Double_t *pt, Double_t *par)
  {
    Double_t ldNdy  = par[0];
    Double_t l2pi   = 2*TMath::Pi();
    Double_t lTemp = par[1];
    Double_t lPower = par[2];
    Double_t lMass  = par[3];
    Double_t A = par[4];
    Double_t lMassEt  = par[5];//only used for calculating Et
    //this is the Et we would calculate if we had identified the particle as having a mass lMassEt
    Double_t Et = TMath::Sqrt(pt[0]*pt[0]+lMassEt*lMassEt)+A*lMassEt;
    
    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);
    Double_t strangeness = (0.4333333*pt[0]-0.0666666)/(0.2333333333*pt[0]-0.01666666);
    if(strangeness<1.0 || pt[0]<0.1) strangeness = 1.0;
    //This is the density of particles times the Et weight
    return ldNdy *strangeness* Et* pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
  }
  ClassDef(AliAnalysisLevyPtModifiedStrangeness, 1);
};
class AliAnalysisLevyPtModifiedBaryonEnhanced{
public:
  virtual ~AliAnalysisLevyPtModifiedBaryonEnhanced(){;};
  

  Double_t Evaluate(Double_t *pt, Double_t *par)
  {
    Double_t ldNdy  = par[0];
    Double_t l2pi   = 2*TMath::Pi();
    Double_t lTemp = par[1];
    Double_t lPower = par[2];
    Double_t lMass  = par[3];
    Double_t A = par[4];
    Double_t lMassEt  = par[5];//only used for calculating Et
    Double_t lBary0 = par[6];
    Double_t lBary1 = par[7];
    Double_t lBary2 = par[8];
    Double_t lBary3 = par[9];
    Double_t lBary4 = par[10];
    Double_t lBary5 = par[11];
    //this is the Et we would calculate if we had identified the particle as having a mass lMassEt
    Double_t Et = TMath::Sqrt(pt[0]*pt[0]+lMassEt*lMassEt)+A*lMassEt;
    
    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);
    //this is the baryon enhancement factor
    Double_t baryonEnhancement = lBary0*pow(pt[0],lBary1)*exp(-pow(pt[0]/lBary2,lBary3))/(lBary4+lBary5*pt[0]);
    if(baryonEnhancement<1.0) baryonEnhancement = 1.0;
    Double_t strangeness = (0.4333333*pt[0]-0.0666666)/(0.2333333333*pt[0]-0.01666666);
    if(strangeness<1.0 || pt[0]<0.1) strangeness = 1.0;
    //This is the density of particles times the Et weight
    return ldNdy *Et*strangeness* baryonEnhancement*pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);  }
  ClassDef(AliAnalysisLevyPtModifiedBaryonEnhanced, 1);
};
class AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced{
public:
  virtual ~AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced(){;};
  

  Double_t Evaluate(Double_t *pt, Double_t *par)
  {
    Double_t ldNdy  = par[0];
    Double_t l2pi   = 2*TMath::Pi();
    Double_t lTemp = par[1];
    Double_t lPower = par[2];
    Double_t lMass  = par[3];
    Double_t A = par[4];
    Double_t lMassEt  = par[5];//only used for calculating Et
    Double_t lBary0 = par[6];
    Double_t lBary1 = par[7];
    Double_t lBary2 = par[8];
    Double_t lBary3 = par[9];
    Double_t lBary4 = par[10];
    Double_t lBary5 = par[11];
    //this is the Et we would calculate if we had identified the particle as having a mass lMassEt
    Double_t Et = TMath::Sqrt(pt[0]*pt[0]+lMassEt*lMassEt)+A*lMassEt;
    
    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);
    //this is the baryon enhancement factor
    Double_t baryonEnhancement = lBary0*pow(pt[0],lBary1)*exp(-pow(pt[0]/lBary2,lBary3))/(lBary4+lBary5*pt[0]);
    if(baryonEnhancement<1.0) baryonEnhancement = 1.0;
    //This is the density of particles times the Et weight
    return ldNdy *Et* baryonEnhancement*pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);  }
  ClassDef(AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced, 1);
};

void CorrPIDLevyFit(bool hadronic = false){
  //pt cuts where we can identify things
  float protoncut = 0.9+0.25;
  float protoncuterr = 0.25;
  float kaoncut = 0.45+0.25;
  float kaoncuterr = 0.25;
  //pt cut
  float ptcut = 0.1;
// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $\pi^{+}+\pi^{-}$ & 2.977  $\pm$ 0.15  & 0.126 $\pm$ 0.001 & 7.82 $\pm$ 0.1  &                  &   &     \\
  AliAnalysisLevyPtModified *function = new AliAnalysisLevyPtModified();
  fPion = new TF1("fPion",function, &AliAnalysisLevyPtModified::Evaluate,0,50,6,"AliAnalysisLevyPtModified","Evaluate");
  fPion->SetParameter(0,2.977);//dN/dy
  fPion->SetParameter(1,0.126);//T
  fPion->SetParameter(2,7.82);//n
  fPion->SetParameter(3,0.13957);//mass
  fPion->SetParameter(4,0);//A=0 for pions
  fPion->SetParameter(5,0.13957);//mass et
  float integralPion = fPion->Integral(ptcut,50);
  float integralErrPion = integralPion*0.007/2.977;
  float myerrorPionT = 0.0;
  float myerrorPionn = 0.0;
  float tmpPion;
  float T = 0.126;
  float Terr = 0.001;
  float n = 7.82;
  float nerr = 0.1;
  fPion->SetParameter(1,T+Terr);//T
  tmpPion = fPion->Integral(ptcut,50);
  myerrorPionT = TMath::Abs(integralPion-tmpPion);
  fPion->SetParameter(1,T-Terr);//T
  tmpPion = fPion->Integral(ptcut,50);
  if(TMath::Abs(integralPion-tmpPion)>myerrorPionT) myerrorPionT = TMath::Abs(integralPion-tmpPion);
  fPion->SetParameter(1,T);//T
  fPion->SetParameter(2,n+nerr);//n
  tmpPion = fPion->Integral(ptcut,50);
  myerrorPionn = TMath::Abs(integralPion-tmpPion);
  fPion->SetParameter(2,n-nerr);//n
  tmpPion = fPion->Integral(ptcut,50);
  if(TMath::Abs(integralPion-tmpPion)>myerrorPionn) myerrorPionn = TMath::Abs(integralPion-tmpPion);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  //To get the correct error one would have to fit the spectra data to get the covariance matrix...
  cout<<"Errors: dN/dy "<<integralErrPion<<" T "<<myerrorPionT<<" n "<<myerrorPionn<<endl;
  integralErrPion = TMath::Sqrt(TMath::Power(integralErrPion,2)+TMath::Power(myerrorPionT,2)+TMath::Power(myerrorPionn,2));
  cout<<"Pion Et = "<<integralPion<<"$\\pm$"<<integralErrPion<<endl;




// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $K^{+}+K^{-}$     & 0.366 $\pm$ 0.03   & 0.160 $\pm$ 0.006 & 6.087 $\pm$ 0.4 &                  &   &     \\
  fKaon = new TF1("fKaon",function, &AliAnalysisLevyPtModified::Evaluate,0,50,6,"AliAnalysisLevyPtModified","Evaluate");
  fKaon->SetParameter(0,0.366);//dN/dy
  fKaon->SetParameter(1,0.160);//T
  fKaon->SetParameter(2,6.087);//n
  fKaon->SetParameter(3,0.493677);//mass
  fKaon->SetParameter(4,0);//A=0 for kaons
  fKaon->SetParameter(5,0.493677);//mass et
  fKaonStrange = new TF1("fKaonStrange",function, &AliAnalysisLevyPtModifiedStrangeness::Evaluate,0,50,6,"AliAnalysisLevyPtModifiedStrangeness","Evaluate");
  for(int i=0;i<6;i++){fKaonStrange->SetParameter(i,fKaon->GetParameter(i));}
  float integralKaon = fKaon->Integral(ptcut,50);
  float integralKaonStrange = fKaonStrange->Integral(ptcut,50);
  float integralErrKaon = integralKaon*0.006/0.366;
  float myerrorKaonT = 0.0;
  float myerrorKaonn = 0.0;
  float tmpKaon;
  float T = 0.160;
  float Terr = 0.006;
  float n = 6.087;
  float nerr = 0.4;
  fKaon->SetParameter(1,T+Terr);//T
  tmpKaon = fKaon->Integral(ptcut,50);
  myerrorKaonT = TMath::Abs(integralKaon-tmpKaon);
  fKaon->SetParameter(1,T-Terr);//T
  tmpKaon = fKaon->Integral(ptcut,50);
  if(TMath::Abs(integralKaon-tmpKaon)>myerrorKaonT) myerrorKaonT = TMath::Abs(integralKaon-tmpKaon);
  fKaon->SetParameter(1,T);//T
  fKaon->SetParameter(2,n+nerr);//n
  tmpKaon = fKaon->Integral(ptcut,50);
  myerrorKaonn = TMath::Abs(integralKaon-tmpKaon);
  fKaon->SetParameter(2,n-nerr);//n
  tmpKaon = fKaon->Integral(ptcut,50);
  if(TMath::Abs(integralKaon-tmpKaon)>myerrorKaonn) myerrorKaonn = TMath::Abs(integralKaon-tmpKaon);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  cout<<"Errors: dN/dy "<<integralErrKaon<<" T "<<myerrorKaonT<<" n "<<myerrorKaonn<<endl;
  integralErrKaon = TMath::Sqrt(TMath::Power(integralErrKaon,2)+TMath::Power(myerrorKaonT,2)+TMath::Power(myerrorKaonn,2));
  cout<<"Kaon Et = "<<integralKaon<<"$\\pm$"<<integralErrKaon<<endl;
  cout<<"Kaon Et Strangeness Enhanced = "<<integralKaonStrange<<endl;
  //now calculate the integral for kaons we can identify:
  float integralKaonLowPtMean = fKaon->Integral(ptcut,kaoncut);
  float integralKaonLowPtPlus = fKaon->Integral(ptcut,kaoncut+kaoncuterr);
  float integralKaonLowPtMinus = fKaon->Integral(ptcut,kaoncut-kaoncuterr);
  float integralKaonStrangeLowPtMinus = fKaonStrange->Integral(ptcut,kaoncut-kaoncuterr);
  //now we switch the kaon mass to the pion mass
  fKaon->SetParameter(5,0.13957);//mass et
  fKaonStrange->SetParameter(5,0.13957);//mass et
  //and do the integrals
  float integralKaonHighPtMean = fKaon->Integral(kaoncut,50);
  float integralKaonHighPtPlus = fKaon->Integral(kaoncut+kaoncuterr,50);
  float integralKaonHighPtMinus = fKaon->Integral(kaoncut-kaoncuterr,50);
  float integralKaonStrangeHighPtMinus = fKaonStrange->Integral(kaoncut-kaoncuterr,50);
  float integralKaonNoID = fKaon->Integral(ptcut,50);
  cout<<"Kaon Et as pion = "<<integralKaonHighPtMean+integralKaonLowPtMean<<endl;
  cout<<"Kaon Et Strange as pion = "<<integralKaonStrangeHighPtMinus+integralKaonStrangeLowPtMinus<<endl;
  TF1 *fKaonPion = fKaon->Clone("fKaonPion");
  //and change it back just to avoid confusion
  fKaon->SetParameter(5,0.493677);//mass et
  fKaonStrange->SetParameter(5,0.493677);//mass et

// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $p + \bar{p}$     &  0.157 $\pm$ 0.012 & 0.196 $\pm$ 0.009 & 8.6 $\pm$ 1.1   &                  &   &     \\
  fProton = new TF1("fProton",function, &AliAnalysisLevyPtModified::Evaluate,0,50,6,"AliAnalysisLevyPtModified","Evaluate");
  fProton->SetParameter(0,0.08);//dN/dy
  fProton->SetParameter(1,0.196);//T
  fProton->SetParameter(2,8.6);//n
  fProton->SetParameter(3,0.938272);//mass
  fProton->SetParameter(4,-1);//A=0 for protons
  fProton->SetParameter(5,0.938272);//mass et
  fProtonEnhanced = new TF1("fProtonEnhanced",function, &AliAnalysisLevyPtModifiedBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedBaryonEnhanced","Evaluate");
  for(int i=0;i<6;i++){fProtonEnhanced->SetParameter(i,fProton->GetParameter(i));}//set all of the spectra parameters to their normal values
  //and now set the baryon enhancement parameters
  fProtonEnhanced->SetParameter(6,0.900878);
  fProtonEnhanced->SetParameter(7,1.38882);
  fProtonEnhanced->SetParameter(8,2.6361);
  fProtonEnhanced->SetParameter(9,1.37751);
  fProtonEnhanced->SetParameter(10,0.5);
  fProtonEnhanced->SetParameter(11,-0.03);
  float integralProton = fProton->Integral(ptcut,50);
  float integralErrProton = integralProton*0.002/0.08;
  float myerrorProtonT = 0.0;
  float myerrorProtonn = 0.0;
  float tmpProton;
  float T = 0.196;
  float Terr = 0.009;
  float n = 8.6;
  float nerr = 1.1;
  fProton->SetParameter(1,T+Terr);//T
  tmpProton = fProton->Integral(ptcut,50);
  myerrorProtonT = TMath::Abs(integralProton-tmpProton);
  fProton->SetParameter(1,T-Terr);//T
  tmpProton = fProton->Integral(ptcut,50);
  if(TMath::Abs(integralProton-tmpProton)>myerrorProtonT) myerrorProtonT = TMath::Abs(integralProton-tmpProton);
  fProton->SetParameter(1,T);//T
  fProton->SetParameter(2,n+nerr);//n
  tmpProton = fProton->Integral(ptcut,50);
  myerrorProtonn = TMath::Abs(integralProton-tmpProton);
  fProton->SetParameter(2,n-nerr);//n
  tmpProton = fProton->Integral(ptcut,50);
  if(TMath::Abs(integralProton-tmpProton)>myerrorProtonn) myerrorProtonn = TMath::Abs(integralProton-tmpProton);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  cout<<"Errors: dN/dy "<<integralErrProton<<" T "<<myerrorProtonT<<" n "<<myerrorProtonn<<endl;
  integralErrProton = TMath::Sqrt(TMath::Power(integralErrProton,2)+TMath::Power(myerrorProtonT,2)+TMath::Power(myerrorProtonn,2));
  cout<<"Proton Et = "<<integralProton<<"$\\pm$"<<integralErrProton<<endl;
  //now calculate the integral for kaons we can identify:
  float integralProtonLowPtMean = fProton->Integral(ptcut,protoncut);
  float integralProtonLowPtPlus = fProton->Integral(ptcut,protoncut+protoncuterr);
  float integralProtonLowPtMinus = fProton->Integral(ptcut,protoncut-protoncuterr);
  float integralProtonEnhancedTrue =integralProtonLowPtMinus+ fProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fProtonEnhanced->SetParameter(6,1.2*fProtonEnhanced->GetParameter(6));
  float integralProtonReallyEnhancedTrue = integralProtonLowPtMinus+ fProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fProtonEnhanced->SetParameter(6,1.0/1.2*fProtonEnhanced->GetParameter(6));
  //now we switch the kaon mass to the pion mass
  fProton->SetParameter(5,0.13957);//mass et
  fProton->SetParameter(4,0);//A=0 for pions
  //and do the integrals
  float integralProtonHighPtMean = fProton->Integral(protoncut,50);
  float integralProtonHighPtPlus = fProton->Integral(protoncut+protoncuterr,50);
  float integralProtonHighPtMinus = fProton->Integral(protoncut-protoncuterr,50);

  //now we switch the kaon mass to the pion mass
  fProtonEnhanced->SetParameter(5,0.13957);//mass et
  fProtonEnhanced->SetParameter(4,0);//A=0 for pions
  float integralProtonHighPtMinusEnhanced = fProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fProtonEnhanced->SetParameter(6,1.2*fProtonEnhanced->GetParameter(6));
  float integralProtonHighPtMinusReallyEnhanced = fProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fProtonEnhanced->SetParameter(6,1.0/1.2*fProtonEnhanced->GetParameter(6));

  float integralProtonNoID = fProton->Integral(ptcut,50);
  float integralProtonEnhancedNoID = integralProtonHighPtMinusEnhanced+integralProtonLowPtMinus;
  float integralProtonReallyEnhancedNoID = integralProtonHighPtMinusReallyEnhanced+integralProtonLowPtMinus;
  cout<<"Proton Et as pion = "<<integralProtonHighPtMean+integralProtonLowPtMean<<endl;
  cout<<"Proton Et really enhanced = "<<integralProtonReallyEnhancedTrue<<"("<<integralProtonHighPtMinusReallyEnhanced+integralProtonLowPtMinus<<")"<<endl;
  TF1 *fProtonPion = fProton->Clone("fProtonPion");
  //and change it back just to avoid confusion
  fProton->SetParameter(5,0.938272);//mass et
  fProton->SetParameter(4,-1);//A=0 for protons
  fProtonEnhanced->SetParameter(5,0.938272);//mass et
  fProtonEnhanced->SetParameter(4,-1);//A=0 for protons



  //Antiprotons...
  fProton->SetParameter(0,0.077);//dN/dy
  fProton->SetParameter(2,n);//n
  fProton->SetParameter(4,1);//A=0 for protons
  fAntiProtonEnhanced = new TF1("fAntiProtonEnhanced",function, &AliAnalysisLevyPtModifiedBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedBaryonEnhanced","Evaluate");
  for(int i=0;i<6;i++){fAntiProtonEnhanced->SetParameter(i,fProton->GetParameter(i));}//set all of the spectra parameters to their normal values
  //and now set the baryon enhancement parameters
  fAntiProtonEnhanced->SetParameter(6,0.900878);
  fAntiProtonEnhanced->SetParameter(7,1.38882);
  fAntiProtonEnhanced->SetParameter(8,2.6361);
  fAntiProtonEnhanced->SetParameter(9,1.37751);
  fAntiProtonEnhanced->SetParameter(10,0.5);
  fAntiProtonEnhanced->SetParameter(11,-0.03);
  float integralAntiProton = fProton->Integral(ptcut,50);
  float integralErrAntiProton = integralAntiProton*0.002/0.077;
  float myerrorAntiProtonT = 0.0;
  float myerrorAntiProtonn = 0.0;
  float tmpAntiProton;
  float T = 0.196;
  float Terr = 0.009;
  float n = 8.6;
  float nerr = 1.1;
  fProton->SetParameter(1,T+Terr);//T
  tmpAntiProton = fProton->Integral(ptcut,50);
  myerrorAntiProtonT = TMath::Abs(integralAntiProton-tmpAntiProton);
  fProton->SetParameter(1,T-Terr);//T
  tmpAntiProton = fProton->Integral(ptcut,50);
  if(TMath::Abs(integralAntiProton-tmpAntiProton)>myerrorAntiProtonT) myerrorAntiProtonT = TMath::Abs(integralAntiProton-tmpAntiProton);
  fProton->SetParameter(1,T);//T
  fProton->SetParameter(2,n+nerr);//n
  tmpAntiProton = fProton->Integral(ptcut,50);
  myerrorAntiProtonn = TMath::Abs(integralAntiProton-tmpAntiProton);
  fProton->SetParameter(2,n-nerr);//n
  tmpAntiProton = fProton->Integral(ptcut,50);
  if(TMath::Abs(integralAntiProton-tmpAntiProton)>myerrorAntiProtonn) myerrorAntiProtonn = TMath::Abs(integralAntiProton-tmpAntiProton);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  cout<<"Errors: dN/dy "<<integralErrAntiProton<<" T "<<myerrorAntiProtonT<<" n "<<myerrorAntiProtonn<<endl;
  integralErrAntiProton = TMath::Sqrt(TMath::Power(integralErrAntiProton,2)+TMath::Power(myerrorAntiProtonT,2)+TMath::Power(myerrorAntiProtonn,2));
  cout<<"AntiProton Et = "<<integralAntiProton<<"$\\pm$"<<integralErrAntiProton<<endl;
  TF1 *fAntiProton = fProton->Clone("fAntiProton");
  //now calculate the integral for kaons we can identify:
  float integralAntiProtonLowPtMean = fProton->Integral(ptcut,protoncut);
  float integralAntiProtonLowPtPlus = fProton->Integral(ptcut,protoncut+protoncuterr);
  float integralAntiProtonLowPtMinus = fProton->Integral(ptcut,protoncut-protoncuterr);
  float integralAntiProtonEnhancedTrue =integralAntiProtonLowPtMinus+ fAntiProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fAntiProtonEnhanced->SetParameter(6,1.2*fAntiProtonEnhanced->GetParameter(6));
  float integralAntiProtonReallyEnhancedTrue = integralAntiProtonLowPtMinus+ fAntiProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fAntiProtonEnhanced->SetParameter(6,1.0/1.2*fAntiProtonEnhanced->GetParameter(6));
  //now we switch the kaon mass to the pion mass
  fProton->SetParameter(5,0.13957);//mass et
  fProton->SetParameter(4,0);//A=0 for pions
  //and do the integrals
  float integralAntiProtonHighPtMean = fProton->Integral(protoncut,50);
  float integralAntiProtonHighPtPlus = fProton->Integral(protoncut+protoncuterr,50);
  float integralAntiProtonHighPtMinus = fProton->Integral(protoncut-protoncuterr,50);

  //now we switch the kaon mass to the pion mass
  fAntiProtonEnhanced->SetParameter(5,0.13957);//mass et
  fAntiProtonEnhanced->SetParameter(4,0);//A=0 for pions
  float integralAntiProtonHighPtMinusEnhanced = fAntiProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fAntiProtonEnhanced->SetParameter(6,1.2*fAntiProtonEnhanced->GetParameter(6));
  float integralAntiProtonHighPtMinusReallyEnhanced = fAntiProtonEnhanced->Integral(protoncut-protoncuterr,50);
  fAntiProtonEnhanced->SetParameter(6,1.0/1.2*fAntiProtonEnhanced->GetParameter(6));

  float integralAntiProtonNoID = fProton->Integral(ptcut,50);
  cout<<"AntiProton Et as pion = "<<integralAntiProtonHighPtMean+integralAntiProtonLowPtMean<<endl;
  cout<<"AntiProton Et really enhanced = "<<integralAntiProtonReallyEnhancedTrue<<"("<<integralAntiProtonHighPtMinusReallyEnhanced+integralAntiProtonLowPtMinus<<")"<<endl;
  TF1 *fAntiProtonPion = fProton->Clone("fAntiProtonPion");
  //and change it back just to avoid confusion
  fProton->SetParameter(5,0.938272);//mass
  //cout<<"CHECK "<<fProton->GetParameter(5)<<endl;
  fProton->SetParameter(4,-1);//A=0 for protons
  fAntiProtonEnhanced->SetParameter(5,0.938272);//mass et
  fAntiProtonEnhanced->SetParameter(4,-1);//A=0 for protons

  float totTrue = integralPion + integralKaon + integralProton + integralAntiProton;
  float totTrueErr = TMath::Sqrt(TMath::Power(integralErrAntiProton,2)+TMath::Power(integralErrProton,2)+TMath::Power(integralErrKaon,2)+TMath::Power(integralErrPion,2));
  float totTrueEnhanced = integralPion + integralKaon + integralProtonEnhancedTrue + integralAntiProtonEnhancedTrue;
  float totTrueReallyEnhanced = integralPion + integralKaon + integralProtonReallyEnhancedTrue + integralAntiProtonReallyEnhancedTrue;
  float totTrueStrange = integralPion + integralKaonStrange + integralProton + integralAntiProton;
  float totTrueStrangeEnhanced = integralPion + integralKaonStrange + integralProtonReallyEnhancedTrue + integralAntiProtonReallyEnhancedTrue;
  cout<<"totEt "<<totTrue<<"+/-"<<totTrueErr<<endl;

  float measEt = integralPion+integralKaonLowPtMean+integralKaonHighPtMean+integralProtonLowPtMean+integralProtonHighPtMean+integralAntiProtonLowPtMean+integralAntiProtonHighPtMean;
  float measEtPlus = integralPion+integralKaonLowPtPlus+integralKaonHighPtPlus+integralProtonLowPtPlus+integralProtonHighPtPlus+integralAntiProtonLowPtPlus+integralAntiProtonHighPtPlus;
  float measEtMinus = integralPion+integralKaonLowPtMinus+integralKaonHighPtMinus+integralProtonLowPtMinus+integralProtonHighPtMinus+integralAntiProtonLowPtMinus+integralAntiProtonHighPtMinus;
  float measEtMinusStrange = integralPion+integralKaonStrangeLowPtMinus+integralKaonStrangeHighPtMinus+integralProtonLowPtMinus+integralProtonHighPtMinus+integralAntiProtonLowPtMinus+integralAntiProtonHighPtMinus;
  float measEtEnhanced = integralPion+integralKaonLowPtMinus+integralKaonHighPtMinus+integralProtonLowPtMinus+integralProtonHighPtMinusEnhanced+integralAntiProtonLowPtMinus+integralAntiProtonHighPtMinusEnhanced;
  float measEtReallyEnhanced = integralPion+integralKaonLowPtMinus+integralKaonHighPtMinus+integralProtonLowPtMinus+integralProtonHighPtMinusReallyEnhanced+integralAntiProtonLowPtMinus+integralAntiProtonHighPtMinusReallyEnhanced;
  cout<<"measEt "<<measEt<<" "<<measEtPlus<<" "<<measEtMinus<<endl;
  float measEtStrangeEnhanced = integralPion+integralKaonStrangeLowPtMinus+integralKaonStrangeHighPtMinus+integralProtonLowPtMinus+integralProtonHighPtMinusReallyEnhanced+integralAntiProtonLowPtMinus+integralAntiProtonHighPtMinusReallyEnhanced;
  cout<<"measEt "<<measEt<<" "<<measEtPlus<<" "<<measEtMinus<<endl;
  float measNoID = integralPion+integralKaonNoID+integralProtonNoID+integralAntiProtonNoID;

  float fpid = measEt/totTrue;
  float fpiderr = fpid*totTrueErr/totTrue;
  float fpidsyserr = TMath::Abs(measEtMinus-measEt)/totTrue;
  float fnopid = measNoID/totTrue;
  float fnopiderr = fnopid*totTrueErr/totTrue;
  if(TMath::Abs(measEtPlus-measEt)> fpidsyserr) fpidsyserr = TMath::Abs(measEtPlus-measEt);
  cout<<"fpid "<<fpid<<" +/- "<<fpiderr<<" +/- "<<fpidsyserr<<endl;
  cout<<"fpid (no pt uncertainty) "<<measEtMinus/totTrue<<" +/- "<<(measEtMinus/totTrue)*totTrueErr/totTrue<<endl;
  cout<<"For no id fpid "<<fnopid<<" +/- "<<fnopiderr<<endl;
  cout<<"fpid enhanced "<<measEtEnhanced/totTrueEnhanced<<endl;
  cout<<"fpid strange "<<measEtMinusStrange/totTrueStrange<<endl;
  cout<<"fpid strange enhanced "<<measEtStrangeEnhanced/totTrueStrangeEnhanced<<endl;
  cout<<"fpid really enhanced "<<measEtReallyEnhanced/totTrueReallyEnhanced<<endl;



  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c1 = new TCanvas("c1","c1",500,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLogy();
  fPion->SetLineColor(1);
  fKaon->SetLineColor(2);
  fKaonPion->SetLineColor(2);
  fKaonPion->SetLineStyle(2);
  fProton->SetLineColor(4);
  fProtonPion->SetLineColor(4);
  fProtonPion->SetLineStyle(2);
  fAntiProton->SetLineColor(TColor::kCyan);
  fAntiProtonPion->SetLineColor(TColor::kCyan);
  fAntiProtonPion->SetLineStyle(2);
//   cout<<"Kaon pion mass "<<Form("%2.4f",fKaonPion->GetParameter(5))<<" A "<<fKaonPion->GetParameter(4)<<endl;
//   cout<<"Proton pion mass "<<fProtonPion->GetParameter(5)<<" A "<<fProtonPion->GetParameter(4)<<endl;
//   cout<<"AntiProton pion mass "<<fAntiProtonPion->GetParameter(5)<<" A "<<fAntiProtonPion->GetParameter(4)<<endl;
//   cout<<"Kaon mass "<<Form("%2.4f",fKaon->GetParameter(5))<<" A "<<fKaon->GetParameter(4)<<endl;
//   cout<<"Proton mass "<<fProton->GetParameter(5)<<" A "<<fProton->GetParameter(4)<<endl;
//   cout<<"AntiProton mass "<<fAntiProton->GetParameter(5)<<" A "<<fAntiProton->GetParameter(4)<<endl;
  TH1F *frame = new TH1F("frame","frame",1,0,2);
  frame->GetYaxis()->SetTitle("1/N_{eve}d^{2}NE_{T}/dydp_{T}");
  frame->GetXaxis()->SetTitle("p_{T}");
  //fPion->SetRange(0,2);
  frame->SetMinimum(1e-4);
  frame->SetMaximum(0.5);
  frame->Draw();
  fPion->Draw("same");
  fKaon->Draw("same");
  fProton->Draw("same");
  fAntiProton->Draw("same");
  fKaonPion->Draw("same");
  fProtonPion->Draw("same");
  fAntiProtonPion->Draw("same");
  TLegend *leg = new  TLegend(0.782258,0.774194,0.925403,0.962366);
  leg->AddEntry(fPion,"#pi");
  leg->AddEntry(fProton,"p");
  leg->AddEntry(fKaon,"K");
  leg->AddEntry(fAntiProton,"#bar{p}");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
 leg->Draw();
 TLine *lineProton = new TLine(protoncut-protoncuterr,frame->GetMinimum(),protoncut-protoncuterr,frame->GetMaximum());
 lineProton->SetLineColor(fProton->GetLineColor());
 lineProton->SetLineStyle(3);
 lineProton->SetLineWidth(2);
 lineProton->Draw();
 TLine *lineKaon = new TLine(kaoncut-kaoncuterr,frame->GetMinimum(),kaoncut-kaoncuterr,frame->GetMaximum());
 lineKaon->SetLineColor(fKaon->GetLineColor());
 lineKaon->SetLineStyle(3);
 lineKaon->SetLineWidth(2);
 lineKaon->Draw();
 

  TCanvas *c2 = new TCanvas("c2","c2",500,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  c2->SetLogy();
  fPionSpectrum = (TF1*) fPion->Clone("fPionSpectrum");
  fKaonSpectrum = (TF1*) fKaon->Clone("fKaonSpectrum");
  fProtonSpectrum = (TF1*) fProton->Clone("fProtonSpectrum");
  fAntiProtonSpectrum = (TF1*) fAntiProton->Clone("fAntiProtonSpectrum");
  fPionSpectrum->SetLineColor(1);
  fKaonSpectrum->SetLineColor(2);
  fProtonSpectrum->SetLineColor(4);
  fAntiProtonSpectrum->SetLineColor(TColor::kCyan);
  fKaonSpectrum->SetParameter(5,0.0);//mass et
  fProtonSpectrum->SetParameter(5,0.0);//mass et
  fPionSpectrum->SetParameter(5,0.0);//mass et
  fAntiProtonSpectrum->SetParameter(5,0.0);//mass et
  TH1F *frame2 = new TH1F("frame2","frame2",1,0,2);
  frame2->GetYaxis()->SetTitle("1/N_{eve}d^{2}N/dydp_{T}");
  frame2->GetXaxis()->SetTitle("p_{T}");
  //fPionSpectrum->SetRange(0,2);
  frame2->SetMinimum(1e-5);
  frame2->SetMaximum(1);
  frame2->Draw();
  fPionSpectrum->Draw("same");
  fKaonSpectrum->Draw("same");
  fProtonSpectrum->Draw("same");
  //fAntiProtonSpectrum->Draw("same");
  TLegend *leg2 = new  TLegend(0.782258,0.774194,0.925403,0.962366);
  leg2->AddEntry(fPionSpectrum,"#pi");
  leg2->AddEntry(fProtonSpectrum,"p,#bar{p}");
  leg2->AddEntry(fKaonSpectrum,"K");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.06);
 leg2->Draw();

  TCanvas *c3 = new TCanvas("c3","c3",500,400);
  c3->SetTopMargin(0.02);
  c3->SetRightMargin(0.02);
  c3->SetBorderSize(0);
  c3->SetFillColor(0);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillColor(0);
  c3->SetFrameBorderMode(0);
  //c3->SetLogy();
  TH1F *frame3 = new TH1F("frame3","frame3",1,0,2);
  frame3->GetYaxis()->SetTitle("E_{T}");
  frame3->GetXaxis()->SetTitle("p_{T}");
  //fPionSpectrum->SetRange(0,2);
  frame3->SetMinimum(0.0);
  frame3->SetMaximum(4);
  frame3->Draw();

  TF1 *fPionEt = new TF1("fPionEt","sqrt(x*x+[0]*[0])",0,2);
  fPionEt->SetParameter(0,0.13957);
  fPionEt->Draw("same");
  TF1 *fKaonEt = new TF1("fKaonEt","sqrt(x*x+[0]*[0])",0,2);
  fKaonEt->SetParameter(0,0.493677);
  fKaonEt->Draw("same");
  TF1 *fProtonEt = new TF1("fProtonEt","sqrt(x*x+[0]*[0])-[0]",0,2);
  fProtonEt->SetParameter(0,0.938272);
  fProtonEt->Draw("same");
  TF1 *fAntiProtonEt = new TF1("fAntiProtonEt","sqrt(x*x+[0]*[0])+[0]",0,2);
  fAntiProtonEt->SetParameter(0,0.938272);
  fAntiProtonEt->Draw("same");
  fPionEt->SetLineColor(1);
  fKaonEt->SetLineColor(2);
  fProtonEt->SetLineColor(4);
  fAntiProtonEt->SetLineColor(TColor::kCyan);
  TLegend *leg3 = new  TLegend(0.782258,0.774194,0.925403,0.962366);
  leg3->AddEntry(fPionEt,"#pi");
  leg3->AddEntry(fProtonEt,"p");
  leg3->AddEntry(fAntiProtonEt,"#bar{p}");
  leg3->AddEntry(fKaonEt,"K");
  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.06);
  leg3->Draw();

  c1->SaveAs("pics/PIDEtweight.eps");
  c2->SaveAs("pics/PIDSpectra.eps");
  c3->SaveAs("pics/PIDEt.eps");


}

