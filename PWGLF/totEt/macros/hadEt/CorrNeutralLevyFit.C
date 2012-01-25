//Christine Nattrass, University of Tennessee at Knoxville
//This macro is to calculate the contributions to Et from various particles based on Levy fits to ALICE data at 900 GeV
//It uses d^2N/dpTdy from the papers, weighs it by Et, and integrates over all pT. 
//A=0 for mesons
//A=1 for antinucleons since they would anihilate in the calorimeter if they were observed by a calorimeter
//A=-1 for nucleons since they would not anihilate and their mass would not be measured
//At the end this is used to calculate the neutral energy correction
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
    Double_t Et = TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)+A*lMass;
    
    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);
    
    return ldNdy *Et* pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
  }
  ClassDef(AliAnalysisLevyPtModified, 1);
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
    Double_t lMassEt  = par[3];//only used for calculating Et
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
  ClassDef(AliAnalysisLevyPtModifiedBaryonEnhanced, 1);
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
    Double_t lMassEt  = par[3];//only used for calculating Et
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
    Double_t lMassEt  = par[3];//only used for calculating Et
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
  ClassDef(AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced, 1);
};


void CorrNeutralLevyFit(bool hadronic = false){

  float factor = 0.0;
  float factorerr = 0.0;
  if(hadronic){
    factor = 0.548;
    factorerr = 0.003;
  }

// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $\pi^{+}+\pi^{-}$ & 2.977  $\pm$ 0.15  & 0.126 $\pm$ 0.001 & 7.82 $\pm$ 0.1  &                  &   &     \\
  AliAnalysisLevyPtModified *function = new AliAnalysisLevyPtModified();
  fPion = new TF1("fPion",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fPion->SetParameter(0,2.977);//dN/dy
  fPion->SetParameter(1,0.126);//T
  fPion->SetParameter(2,7.82);//n
  fPion->SetParameter(3,0.13957);//mass
  fPion->SetParameter(4,0);//A=0 for pions
  float integralPion = fPion->Integral(0,50);
  float integralErrPion = integralPion*0.15/2.977;
  float myerrorPionT = 0.0;
  float myerrorPionn = 0.0;
  float tmpPion;
  float T = 0.126;
  float Terr = 0.001;
  float n = 7.82;
  float nerr = 0.1;
  fPion->SetParameter(1,T+Terr);//T
  tmpPion = fPion->Integral(0,50);
  myerrorPionT = TMath::Abs(integralPion-tmpPion);
  fPion->SetParameter(1,T-Terr);//T
  tmpPion = fPion->Integral(0,50);
  if(TMath::Abs(integralPion-tmpPion)>myerrorPionT) myerrorPionT = TMath::Abs(integralPion-tmpPion);
  fPion->SetParameter(1,T);//T
  fPion->SetParameter(2,n+nerr);//n
  tmpPion = fPion->Integral(0,50);
  myerrorPionn = TMath::Abs(integralPion-tmpPion);
  fPion->SetParameter(2,n-nerr);//n
  tmpPion = fPion->Integral(0,50);
  if(TMath::Abs(integralPion-tmpPion)>myerrorPionn) myerrorPionn = TMath::Abs(integralPion-tmpPion);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  //To get the correct error one would have to fit the spectra data to get the covariance matrix...
  integralErrPion = TMath::Sqrt(TMath::Power(integralErrPion,2)+TMath::Power(myerrorPionT,2)+TMath::Power(myerrorPionn,2));
  cout<<"Pion Et = "<<integralPion<<"$\\pm$"<<integralErrPion<<endl;

// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $K^{+}+K^{-}$     & 0.366 $\pm$ 0.03   & 0.160 $\pm$ 0.006 & 6.087 $\pm$ 0.4 &                  &   &     \\
  fKaon = new TF1("fKaon",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fKaon->SetParameter(0,0.366);//dN/dy
  fKaon->SetParameter(1,0.160);//T
  fKaon->SetParameter(2,6.087);//n
  fKaon->SetParameter(3,0.493677);//mass
  fKaon->SetParameter(4,0);//A=0 for kaons
  float integralKaon = fKaon->Integral(0,50);
  float integralErrKaon = integralKaon*0.03/0.366;
  fKaonStrange = new TF1("fKaonStrange",function, &AliAnalysisLevyPtModifiedStrangeness::Evaluate,0,50,5,"AliAnalysisLevyPtModifiedStrangeness","Evaluate");
  for(int i=0; i<5;i++){fKaonStrange->SetParameter(i,fKaon->GetParameter(i));}
  float integralKaonStrange = fKaonStrange->Integral(0,50);
  float myerrorKaonT = 0.0;
  float myerrorKaonn = 0.0;
  float tmpKaon;
  float T = 0.160;
  float Terr = 0.006;
  float n = 6.087;
  float nerr = 0.4;
  fKaon->SetParameter(1,T+Terr);//T
  tmpKaon = fKaon->Integral(0,50);
  myerrorKaonT = TMath::Abs(integralKaon-tmpKaon);
  fKaon->SetParameter(1,T-Terr);//T
  tmpKaon = fKaon->Integral(0,50);
  if(TMath::Abs(integralKaon-tmpKaon)>myerrorKaonT) myerrorKaonT = TMath::Abs(integralKaon-tmpKaon);
  fKaon->SetParameter(1,T);//T
  fKaon->SetParameter(2,n+nerr);//n
  tmpKaon = fKaon->Integral(0,50);
  myerrorKaonn = TMath::Abs(integralKaon-tmpKaon);
  fKaon->SetParameter(2,n-nerr);//n
  tmpKaon = fKaon->Integral(0,50);
  if(TMath::Abs(integralKaon-tmpKaon)>myerrorKaonn) myerrorKaonn = TMath::Abs(integralKaon-tmpKaon);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  integralErrKaon = TMath::Sqrt(TMath::Power(integralErrKaon,2)+TMath::Power(myerrorKaonT,2)+TMath::Power(myerrorKaonn,2));
  cout<<"Kaon Et = "<<integralKaon<<"$\\pm$"<<integralErrKaon<<endl;
  cout<<"Kaon Et Strange = "<<integralKaonStrange<<endl;

// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $p + \bar{p}$     &  0.157 $\pm$ 0.012 & 0.196 $\pm$ 0.009 & 8.6 $\pm$ 1.1   &                  &   &     \\
  fProton = new TF1("fProton",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fProton->SetParameter(0,0.157/2.0);//dN/dy
  fProton->SetParameter(1,0.196);//T
  fProton->SetParameter(2,8.6);//n
  fProton->SetParameter(3,0.938272);//mass
  fProton->SetParameter(4,-1);//A=0 for protons
  fProtonEnhanced = new TF1("fProtonEnhanced",function, &AliAnalysisLevyPtModifiedBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedBaryonEnhanced","Evaluate");
  for(int i=0;i<6;i++){fProtonEnhanced->SetParameter(i,fProton->GetParameter(i));}//set all of the spectra parameters to their normal values
  //and now set the baryon enhancement parameters
  fProtonEnhanced->SetParameter(6,0.900878*1.2);
  fProtonEnhanced->SetParameter(7,1.38882);
  fProtonEnhanced->SetParameter(8,2.6361);
  fProtonEnhanced->SetParameter(9,1.37751);
  fProtonEnhanced->SetParameter(10,0.5);
  fProtonEnhanced->SetParameter(11,-0.03);
  float integralProtonEnhanced = fProtonEnhanced->Integral(0,50);
  float integralProton = fProton->Integral(0,50);
  float integralErrProton = integralProton*0.012/0.157;
  float myerrorProtonT = 0.0;
  float myerrorProtonn = 0.0;
  float tmpProton;
  float T = 0.196;
  float Terr = 0.009;
  float n = 8.6;
  float nerr = 1.1;
  fProton->SetParameter(1,T+Terr);//T
  tmpProton = fProton->Integral(0,50);
  myerrorProtonT = TMath::Abs(integralProton-tmpProton);
  fProton->SetParameter(1,T-Terr);//T
  tmpProton = fProton->Integral(0,50);
  if(TMath::Abs(integralProton-tmpProton)>myerrorProtonT) myerrorProtonT = TMath::Abs(integralProton-tmpProton);
  fProton->SetParameter(1,T);//T
  fProton->SetParameter(2,n+nerr);//n
  tmpProton = fProton->Integral(0,50);
  myerrorProtonn = TMath::Abs(integralProton-tmpProton);
  fProton->SetParameter(2,n-nerr);//n
  tmpProton = fProton->Integral(0,50);
  if(TMath::Abs(integralProton-tmpProton)>myerrorProtonn) myerrorProtonn = TMath::Abs(integralProton-tmpProton);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  integralErrProton = TMath::Sqrt(TMath::Power(integralErrProton,2)+TMath::Power(myerrorProtonT,2)+TMath::Power(myerrorProtonn,2));
  cout<<"Proton Et = "<<integralProton<<"$\\pm$"<<integralErrProton<<endl;





  //Antiprotons...
  fProton->SetParameter(0,0.157/2.0);//dN/dy
  fProton->SetParameter(1,0.196);//T
  fProton->SetParameter(2,8.6);//n
  fProton->SetParameter(3,0.938272);//mass
  fProton->SetParameter(2,n);//n
  fProton->SetParameter(4,1);//A=0 for protons
  float integralAntiProton = fProton->Integral(0,50);
  float integralErrAntiProton = integralAntiProton*0.012/0.157;
  fAntiProtonEnhanced = new TF1("fAntiProtonEnhanced",function, &AliAnalysisLevyPtModifiedBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedBaryonEnhanced","Evaluate");
  for(int i=0;i<6;i++){fAntiProtonEnhanced->SetParameter(i,fProton->GetParameter(i));}//set all of the spectra parameters to their normal values
  fAntiProtonEnhanced->SetParameter(2,n);//n
  fAntiProtonEnhanced->SetParameter(4,1);//A=0 for protons
  //and now set the baryon enhancement parameters
  fAntiProtonEnhanced->SetParameter(6,0.900878*1.2);
  fAntiProtonEnhanced->SetParameter(7,1.38882);
  fAntiProtonEnhanced->SetParameter(8,2.6361);
  fAntiProtonEnhanced->SetParameter(9,1.37751);
  fAntiProtonEnhanced->SetParameter(10,0.5);
  fAntiProtonEnhanced->SetParameter(11,-0.03);
  float integralAntiProtonEnhanced = fAntiProtonEnhanced->Integral(0,50);
  float myerrorAntiProtonT = 0.0;
  float myerrorAntiProtonn = 0.0;
  float tmpAntiProton;
  float T = 0.196;
  float Terr = 0.009;
  float n = 8.6;
  float nerr = 1.1;
  fProton->SetParameter(1,T+Terr);//T
  tmpAntiProton = fProton->Integral(0,50);
  myerrorAntiProtonT = TMath::Abs(integralAntiProton-tmpAntiProton);
  fProton->SetParameter(1,T-Terr);//T
  tmpAntiProton = fProton->Integral(0,50);
  if(TMath::Abs(integralAntiProton-tmpAntiProton)>myerrorAntiProtonT) myerrorAntiProtonT = TMath::Abs(integralAntiProton-tmpAntiProton);
  fProton->SetParameter(1,T);//T
  fProton->SetParameter(2,n+nerr);//n
  tmpAntiProton = fProton->Integral(0,50);
  myerrorAntiProtonn = TMath::Abs(integralAntiProton-tmpAntiProton);
  fProton->SetParameter(2,n-nerr);//n
  tmpAntiProton = fProton->Integral(0,50);
  if(TMath::Abs(integralAntiProton-tmpAntiProton)>myerrorAntiProtonn) myerrorAntiProtonn = TMath::Abs(integralAntiProton-tmpAntiProton);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  integralErrAntiProton = TMath::Sqrt(TMath::Power(integralErrAntiProton,2)+TMath::Power(myerrorAntiProtonT,2)+TMath::Power(myerrorAntiProtonn,2));
  cout<<"AntiProton Et = "<<integralAntiProton<<"$\\pm$"<<integralErrAntiProton<<endl;
  cout<<"AntiProton enhanced "<<integralAntiProtonEnhanced<<", "<<integralAntiProton/integralAntiProtonEnhanced*100.0<<endl;
  cout<<"Proton enhanced "<<integralProtonEnhanced<<", "<<integralProton/integralProtonEnhanced*100.0<<endl;
  cout<<"Proton and antiproton enhanced "<< (integralProtonEnhanced+integralAntiProtonEnhanced) <<", "
      << (integralProtonEnhanced+integralAntiProtonEnhanced) / (integralProton+integralAntiProton)*100.0 
      <<endl;


// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// \kzeroshort       &0.184 $\pm$ 0.006   & 0.168 $\pm$ 0.005 & 6.6 $\pm$ 0.3   &                  &   &     \\
  fK0S = new TF1("fK0S",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fK0S->SetParameter(0,0.184);//dN/dy
  fK0S->SetParameter(1,0.168);//T
  fK0S->SetParameter(2,6.6);//n
  fK0S->SetParameter(3,.497614);//mass
  fK0S->SetParameter(4,0);//A=0 for kaons
  float integralK0S = fK0S->Integral(0,50);
  float integralErrK0S = integralK0S*0.006/0.184;
  fK0Strange = new TF1("fK0Strange",function, &AliAnalysisLevyPtModifiedStrangeness::Evaluate,0,50,5,"AliAnalysisLevyPtModifiedStrangeness","Evaluate");
  for(int i=0; i<5;i++){fK0Strange->SetParameter(i,fK0S->GetParameter(i));}
  float integralK0SStrange = fK0Strange->Integral(0,50);
  float myerrorK0ST = 0.0;
  float myerrorK0Sn = 0.0;
  float tmpK0S;
  float T = 0.184;
  float Terr = 0.006;
  float n = 6.6;
  float nerr = 0.3;
  fK0S->SetParameter(1,T+Terr);//T
  tmpK0S = fK0S->Integral(0,50);
  myerrorK0ST = TMath::Abs(integralK0S-tmpK0S);
  fK0S->SetParameter(1,T-Terr);//T
  tmpK0S = fK0S->Integral(0,50);
  if(TMath::Abs(integralK0S-tmpK0S)>myerrorK0ST) myerrorK0ST = TMath::Abs(integralK0S-tmpK0S);
  fK0S->SetParameter(1,T);//T
  fK0S->SetParameter(2,n+nerr);//n
  tmpK0S = fK0S->Integral(0,50);
  myerrorK0Sn = TMath::Abs(integralK0S-tmpK0S);
  fK0S->SetParameter(2,n-nerr);//n
  tmpK0S = fK0S->Integral(0,50);
  if(TMath::Abs(integralK0S-tmpK0S)>myerrorK0Sn) myerrorK0Sn = TMath::Abs(integralK0S-tmpK0S);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  integralErrK0S = TMath::Sqrt(TMath::Power(integralErrK0S,2)+TMath::Power(myerrorK0ST,2)+TMath::Power(myerrorK0Sn,2));
  cout<<"K0S Et = "<<integralK0S<<"$\\pm$"<<integralErrK0S<<endl;
  cout<<"K0S Et Strange = "<<integralK0SStrange<<endl;

// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// \lam              &0.048 $\pm$ 0.004   & 0.229 $\pm$ 0.015 & 10.8 $\pm$ 2.0  &                  &   &     \\
  fLambda = new TF1("fLambda",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fLambda->SetParameter(0,0.048);//dN/dy
  fLambda->SetParameter(1,0.229);//T
  fLambda->SetParameter(2,10.8);//n
  fLambda->SetParameter(3,1.115683);//mass
  fLambda->SetParameter(4,0);//A=0 for kaons
  float integralLambda = fLambda->Integral(0,50);
  float integralErrLambda = integralLambda*0.004/0.048;
  fLambdaEnhanced = new TF1("fLambdaEnhanced",function, &AliAnalysisLevyPtModifiedBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedBaryonEnhanced","Evaluate");
  for(int i=0;i<6;i++){fLambdaEnhanced->SetParameter(i,fLambda->GetParameter(i));}//set all of the spectra parameters to their normal values
  //and now set the baryon enhancement parameters
  fLambdaEnhanced->SetParameter(6,0.900878);
  fLambdaEnhanced->SetParameter(7,1.38882);
  fLambdaEnhanced->SetParameter(8,2.6361);
  fLambdaEnhanced->SetParameter(9,1.37751);
  fLambdaEnhanced->SetParameter(10,0.5);
  fLambdaEnhanced->SetParameter(11,-0.03);
  float integralLambdaEnhanced = fLambdaEnhanced->Integral(0,50);
//   fLambdaEnhanced->Draw();
//   return;
  fLambdaStrange = new TF1("fLambdaStrange",function, &AliAnalysisLevyPtModifiedStrangeness::Evaluate,0,50,5,"AliAnalysisLevyPtModifiedStrangeness","Evaluate");
  for(int i=0; i<5;i++){fLambdaStrange->SetParameter(i,fLambda->GetParameter(i));}
  //for(int i=0;i<12;i++){cout<<"Lambda a"<<i<<" "<<fLambdaStrange->GetParameter(i)<<" b"<<i<<" "<<fLambda->GetParameter(i)<<endl;}
  float integralLambdaStrange = fLambdaStrange->Integral(0,50);
  //cout<<"test "<<integralLambdaStrange<<endl;
  fLambdaStrangeEnhanced = new TF1("fLambdaStrangeEnhanced",function, &AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced","Evaluate");
  for(int i=0; i<12;i++){fLambdaStrangeEnhanced->SetParameter(i,fLambdaEnhanced->GetParameter(i));}
//   //and now set the baryon enhancement parameters
//   fLambdaStrangeEnhanced->SetParameter(6,0.900878);
//   fLambdaStrangeEnhanced->SetParameter(7,1.38882);
//   fLambdaStrangeEnhanced->SetParameter(8,2.6361);
//   fLambdaStrangeEnhanced->SetParameter(9,1.37751);
//   fLambdaStrangeEnhanced->SetParameter(10,0.5);
//   fLambdaStrangeEnhanced->SetParameter(11,-0.03);
  //for(int i=0;i<12;i++){cout<<"a"<<i<<" "<<fLambdaStrangeEnhanced->GetParameter(i)<<" b"<<i<<" "<<fLambdaEnhanced->GetParameter(i)<<" c"<<i<<" "<<fLambda->GetParameter(i)<<endl;}
  float integralLambdaStrangeEnhanced = fLambdaStrangeEnhanced->Integral(0,50);
  cout<<"Lambda enhanced "<<integralLambdaEnhanced<<", ";
  if(integralLambdaEnhanced>0.0) cout<<integralLambda/integralLambdaEnhanced*100.0;
  cout<<endl;
  float myerrorLambdaT = 0.0;
  float myerrorLambdan = 0.0;
  float tmpLambda;
  float T = 0.229;
  float Terr = 0.015;
  float n = 10.8;
  float nerr = 2.0;
  fLambda->SetParameter(1,T+Terr);//T
  tmpLambda = fLambda->Integral(0,50);
  myerrorLambdaT = TMath::Abs(integralLambda-tmpLambda);
  fLambda->SetParameter(1,T-Terr);//T
  tmpLambda = fLambda->Integral(0,50);
  if(TMath::Abs(integralLambda-tmpLambda)>myerrorLambdaT) myerrorLambdaT = TMath::Abs(integralLambda-tmpLambda);
  fLambda->SetParameter(1,T);//T
  fLambda->SetParameter(2,n+nerr);//n
  tmpLambda = fLambda->Integral(0,50);
  myerrorLambdan = TMath::Abs(integralLambda-tmpLambda);
  fLambda->SetParameter(2,n-nerr);//n
  tmpLambda = fLambda->Integral(0,50);
  if(TMath::Abs(integralLambda-tmpLambda)>myerrorLambdan) myerrorLambdan = TMath::Abs(integralLambda-tmpLambda);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  integralErrLambda = TMath::Sqrt(TMath::Power(integralErrLambda,2)+TMath::Power(myerrorLambdaT,2)+TMath::Power(myerrorLambdan,2));
  cout<<"Lambda Et = "<<integralLambda<<"$\\pm$"<<integralErrLambda<<endl;
  cout<<"Lambda Et Strange = "<<integralLambdaStrange<<endl;
  cout<<"Lambda Et Strange Enhanced = "<<integralLambdaStrangeEnhanced<<endl;


// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// \alam             & 0.047 $\pm$ 0.005  & 0.210 $\pm$ 0.015 & 9.2 $\pm$ 1.4   &                  &   &     \\                  &   &     \\
  fAntiLambda = new TF1("fAntiLambda",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fAntiLambda->SetParameter(0,0.047);//dN/dy
  fAntiLambda->SetParameter(1,0.210);//T
  fAntiLambda->SetParameter(2,9.2);//n
  fAntiLambda->SetParameter(3,1.115683);//mass
  fAntiLambda->SetParameter(4,0);//A=0 for kaons
  float integralAntiLambda = fAntiLambda->Integral(0,50);
  float integralErrAntiLambda = integralAntiLambda*0.005/0.047;
  fAntiLambdaEnhanced = new TF1("fAntiLambdaEnhanced",function, &AliAnalysisLevyPtModifiedBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedBaryonEnhanced","Evaluate");
  for(int i=0;i<6;i++){fAntiLambdaEnhanced->SetParameter(i,fAntiLambda->GetParameter(i));}//set all of the spectra parameters to their normal values
  //and now set the baryon enhancement parameters
  fAntiLambdaEnhanced->SetParameter(6,0.900878*1.2);
  fAntiLambdaEnhanced->SetParameter(7,1.38882);
  fAntiLambdaEnhanced->SetParameter(8,2.6361);
  fAntiLambdaEnhanced->SetParameter(9,1.37751);
  fAntiLambdaEnhanced->SetParameter(10,0.5);
  fAntiLambdaEnhanced->SetParameter(11,-0.03);
  float integralAntiLambdaEnhanced = fAntiLambdaEnhanced->Integral(0,50);
  cout<<"AntiLambda enhanced "<<integralAntiLambdaEnhanced<<", "<<integralAntiLambda/integralAntiLambdaEnhanced*100.0<<endl;
  fAntiLambdaStrange = new TF1("fAntiLambdaStrange",function, &AliAnalysisLevyPtModifiedStrangeness::Evaluate,0,50,5,"AliAnalysisLevyPtModifiedStrangeness","Evaluate");
  for(int i=0; i<5;i++){fAntiLambdaStrange->SetParameter(i,fAntiLambda->GetParameter(i));}
  float integralAntiLambdaStrange = fAntiLambdaStrange->Integral(0,50);

  fAntiLambdaStrangeEnhanced = new TF1("fAntiLambdaStrangeEnhanced",function, &AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced::Evaluate,0,50,12,"AliAnalysisLevyPtModifiedStrangenessBaryonEnhanced","Evaluate");
  for(int i=0; i<12;i++){fAntiLambdaStrangeEnhanced->SetParameter(i,fAntiLambdaEnhanced->GetParameter(i));}
  float integralAntiLambdaStrangeEnhanced = fAntiLambdaStrangeEnhanced->Integral(0,50);

  float myerrorAntiLambdaT = 0.0;
  float myerrorAntiLambdan = 0.0;
  float tmpAntiLambda;
  float T = 0.210;
  float Terr = 0.015;
  float n = 9.2;
  float nerr = 1.4;
  fAntiLambda->SetParameter(1,T+Terr);//T
  tmpAntiLambda = fAntiLambda->Integral(0,50);
  myerrorAntiLambdaT = TMath::Abs(integralAntiLambda-tmpAntiLambda);
  fAntiLambda->SetParameter(1,T-Terr);//T
  tmpAntiLambda = fAntiLambda->Integral(0,50);
  if(TMath::Abs(integralAntiLambda-tmpAntiLambda)>myerrorAntiLambdaT) myerrorAntiLambdaT = TMath::Abs(integralAntiLambda-tmpAntiLambda);
  fAntiLambda->SetParameter(1,T);//T
  fAntiLambda->SetParameter(2,n+nerr);//n
  tmpAntiLambda = fAntiLambda->Integral(0,50);
  myerrorAntiLambdan = TMath::Abs(integralAntiLambda-tmpAntiLambda);
  fAntiLambda->SetParameter(2,n-nerr);//n
  tmpAntiLambda = fAntiLambda->Integral(0,50);
  if(TMath::Abs(integralAntiLambda-tmpAntiLambda)>myerrorAntiLambdan) myerrorAntiLambdan = TMath::Abs(integralAntiLambda-tmpAntiLambda);
  //This isn't strictly correct because the errors on the parameters should be correlated but it's close
  integralErrAntiLambda = TMath::Sqrt(TMath::Power(integralErrAntiLambda,2)+TMath::Power(myerrorAntiLambdaT,2)+TMath::Power(myerrorAntiLambdan,2));
  cout<<"AntiLambda Et = "<<integralAntiLambda<<"$\\pm$"<<integralErrAntiLambda<<endl;
  cout<<"AntiLambda Et Strange = "<<integralAntiLambdaStrange<<endl;
  cout<<"AntiLambda Et Strange Enhanced = "<<integralAntiLambdaStrangeEnhanced<<endl;

  //now we apply various assumptions
  float integralK0L = integralK0S;
  float integralErrK0L = integralErrK0S;
  float integralNeutron = integralProton;
  float integralErrNeutron = integralErrProton;

  float totalEt = (1.0+factor)*integralPion+integralKaon+2.0*integralProton+2.0*integralAntiProton+2.0*integralK0S+integralLambda+integralAntiLambda;
  float measuredEt = integralPion+integralKaon+1.0*integralProton+1.0*integralAntiProton;
  float fneutral = measuredEt/totalEt;
  cout<<"fneutral = "<<fneutral<<endl;

  float totalEtEnhanced = (1.0+factor)*integralPion+integralKaon+2.0*integralProtonEnhanced+2.0*integralAntiProtonEnhanced+2.0*integralK0S+integralLambdaEnhanced+integralAntiLambdaEnhanced;
  float measuredEtEnhanced = integralPion+integralKaon+1.0*integralProtonEnhanced+1.0*integralAntiProtonEnhanced;
  float fneutralEnhanced = measuredEtEnhanced/totalEtEnhanced;
  cout<<"fneutralEnhanced = "<<fneutralEnhanced<<endl;


  float totalEtStrange = (1.0+factor)*integralPion+integralKaonStrange+2.0*integralProton+2.0*integralAntiProton+2.0*integralK0SStrange+integralLambdaStrange+integralAntiLambdaStrange;
  float measuredEtStrange = integralPion+integralKaonStrange+1.0*integralProton+1.0*integralAntiProton;
  float fneutralStrange = measuredEtStrange/totalEtStrange;
  cout<<"fneutralStrange = "<<fneutralStrange<<endl;



  float totalEtStrangeEnhanced = (1.0+factor)*integralPion+integralKaonStrange+2.0*integralProtonEnhanced+2.0*integralAntiProtonEnhanced+2.0*integralK0SStrange+integralLambdaStrangeEnhanced+integralAntiLambdaStrangeEnhanced;
  float measuredEtStrangeEnhanced = integralPion+integralKaonStrange+1.0*integralProtonEnhanced+1.0*integralAntiProtonEnhanced;
  float fneutralStrangeEnhanced = measuredEtStrangeEnhanced/totalEtStrangeEnhanced;
  cout<<"fneutralStrangeEnhanced = "<<fneutralStrangeEnhanced<<endl;

  //this is from ugly derivative taking for error propagation
  //form 1:  pions, kaons, protons
  //for f=(xa+A)/(ya+B)
  //df/da=(xf-yf^2)/(xa+A)
  //x=y=1; xa+A=measuredEt
  float errPion = (fneutral-fneutral*fneutral)/measuredEt*integralErrPion;
  if(hadronic){
    //then we have x=1, y = 1+factor
    //the error on the extra bit of et from pi0s
    float extraerror = integralPion*factorerr;
    integralErrPion = TMath::Sqrt(TMath::Power(integralErrPion,2)+TMath::Power(extraerror,2));
    errPion = (fneutral-(1.0+factor)*fneutral*fneutral)/measuredEt*integralErrPion;
  }
  float errKaon = (fneutral-fneutral*fneutral)/measuredEt*integralErrKaon;
  //x=1,y=2
  float errProton = (fneutral-2.0*fneutral*fneutral)/measuredEt*integralErrProton;
  float errAntiProton = (fneutral-2.0*fneutral*fneutral)/measuredEt*integralErrAntiProton;
  //form 2: K0S, lambda, antilambda
  //for f=A/(B+xa)
  //df/da=-xf^2/A
  //x=1,A=measuredEt
  float errLambda = fneutral*fneutral/measuredEt*integralErrLambda;
  float errAntiLambda = fneutral*fneutral/measuredEt*integralErrAntiLambda;
  //x=2,A=measuredEt
  float errK0S = 2.0*fneutral*fneutral/measuredEt*integralErrK0S;
  cout<<"Errors Pion "<<errPion<<" Kaon "<<errKaon<<" Proton "<<errProton<<" Lambda "<<errLambda<<" AntiLambda "<<errAntiLambda<<" K0S "<<errK0S<<endl;
  float totalErr = TMath::Sqrt( TMath::Power(errPion,2) + TMath::Power(errKaon,2) + TMath::Power(errProton,2)+ TMath::Power(errAntiProton,2) + TMath::Power(errLambda,2) + TMath::Power(errAntiLambda,2) + TMath::Power(errK0S,2) );
  
  cout<<"fneutral = "<<fneutral<<"$\\pm$"<<totalErr<<endl;
  cout<<"1/fneutral = "<<1.0/fneutral<<"$\\pm$"<<1.0/fneutral*totalErr/fneutral<<endl;
  cout<<"Percentage error "<<totalErr/fneutral*100.0<<endl;

  float fneutralPb = (fneutralStrangeEnhanced+fneutral)/2.0;
  float fneutralPbErr = TMath::Sqrt( TMath::Power((fneutralStrangeEnhanced-fneutral)/2.0,2.0) + TMath::Power(totalErr,2.0));
  cout<<"fneutral strange enhanced "<<fneutralPb<<"$\\pm$"<<fneutralPbErr<<endl;
    cout<<"1/fneutral strange enhanced = "<<1.0/fneutralPb<<"$\\pm$"<<1.0/fneutralPb*fneutralPbErr/fneutralPb<<endl;
  cout<<"Percentage error "<<fneutralPbErr/fneutralPb*100.0<<endl;

}

