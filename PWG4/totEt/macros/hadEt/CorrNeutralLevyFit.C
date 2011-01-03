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

// particle          & $\frac{dN}{dy}$    & T (GeV)           & n               & $\frac{dE_T}{dy}$& a &\ET \\ \hline
// $p + \bar{p}$     &  0.157 $\pm$ 0.012 & 0.196 $\pm$ 0.009 & 8.6 $\pm$ 1.1   &                  &   &     \\
  fProton = new TF1("fProton",function, &AliAnalysisLevyPtModified::Evaluate,0,50,5,"AliAnalysisLevyPtModified","Evaluate");
  fProton->SetParameter(0,0.157/2.0);//dN/dy
  fProton->SetParameter(1,0.196);//T
  fProton->SetParameter(2,8.6);//n
  fProton->SetParameter(3,0.938272);//mass
  fProton->SetParameter(4,-1);//A=0 for protons
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
  fProton->SetParameter(2,n);//n
  fProton->SetParameter(4,1);//A=0 for protons
  float integralAntiProton = fProton->Integral(0,50);
  float integralErrAntiProton = integralAntiProton*0.012/0.157;
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

  //now we apply various assumptions
  float integralK0L = integralK0S;
  float integralErrK0L = integralErrK0S;
  float integralNeutron = integralProton;
  float integralErrNeutron = integralErrProton;

  float totalEt = (1.0+factor)*integralPion+integralKaon+2.0*integralProton+2.0*integralAntiProton+2.0*integralK0S+integralLambda+integralAntiLambda;
  float measuredEt = integralPion+integralKaon+1.0*integralProton+1.0*integralAntiProton;
  float fneutral = measuredEt/totalEt;
  cout<<"fneutral = "<<fneutral<<endl;
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
}

