/*
  Formulas:

  z  = s* (z0 - vd*(t-t0))

  s  - side -1 and +1 
  t0 - time 0
  vd - nominal drift velocity
  zs - miscalibrated position

  zs = s*(z0 - vd*(1+vr)*(t-(t0+dt))
  vr  - relative change of the drift velocity
  dzt - vd*dt
  dr  = zz0-s*z
  ..
  ==>
  zs ~ z - s*vr*(z0-s*z)+s*dzt
  --------------------------------
  1. Correction function vr constant:


  dz = zs-z = -s*vr *(z0-s*z)+s*dzt         
  dzs/dl = dz/dl +s*s*vr*dz/dl 
  d(dz/dl) = vr*dz/dl 

*/




gSystem->Load("libSTAT.so");
TStatToolkit toolkit;
Double_t chi2;
TVectorD fitParam;
TMatrixD covMatrix;
Int_t npoints;



TFile f("driftitsTPC.root");
TTree * tree = (TTree*)f.Get("Test");

tree->SetAlias("side","(-1+(pTPC.fP[1]>0)*2)");     //side
tree->SetAlias("z","(pTPC.fP[1]+0.0)");     //z position
tree->SetAlias("dr","(1-abs(pTPC.fP[1])/250.)");    //norm drift length
tree->SetAlias("tl","(pTPC.fP[3]+pITS.fP[3])*0.5"); //tan lampbda
tree->SetAlias("sa","sin(pTPC.fAlpha+0.)");         //sin alpha
tree->SetAlias("ca","cos(pTPC.fAlpha+0.)");         //cos alpha

tree->SetAlias("dz","(pTPC.fP[1]-pITS.fP[1])");      //z delta
tree->SetAlias("dtl","(pTPC.fP[3]-pITS.fP[3])");     //delta tan lampbda
tree->SetAlias("etl","sqrt(pITS.fC[9]+0.)");         //error tan lampbda
tree->SetAlias("ez","sqrt(pITS.fC[2]+0.)");          //error z

TCut cut0("abs(pITS.fP[0]-pTPC.fP[0])<2&&abs(pTPC.fP[1])>10");
//TCut cutOut="(side<0)"
TCut cutA=cut0;



TCut cut38009("run==38009");  // 524
TCut cut38010("run==38010");  // 60

TCut cut38286("run==38286");
TCut cut38532("run==38532");
TCut cut38576("run==38576");
TCut cut38586("run==38586");
TCut cut38591("run==38591");
TCut cut45639("run==45639");
TCut cut46993("run==46993");
TCut cut47164("run==47164");
TCut cut47175("run==47175");
TCut cut47180("run==47180");
TCut cut47274("run==47274");

TCut cutA = cut38591+cut0;
TCut cut1 = "1";
TCut cut2 = "1";

void Fit1(){
  TString fstringTl1="";
  fstringTl1+="(side)++";
  fstringTl1+="(tl)++";
  fstringTl1+="side*(tl)++";
  
  //
  //
  TString fstringZ1="";
  fstringZ1+="(side)++";
  fstringZ1+="(-side*(250.0-side*z))++";
  fstringZ1+="((250.-side*z))++";
  // systematic  shift
  fstringZ1+="sa++";   
  fstringZ1+="ca++";
  fstringZ1+="side*sa++";
  fstringZ1+="side*ca++";
  //

  TString *strTl1 = toolkit.FitPlane(tree,"dtl:etl",fstringTl1->Data(), cutA+cut1, chi2,npoints,fitParam,covMatrix);
  strTl1->Tokenize("+")->Print();
  printf("Tl1: Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitTl1",strTl1->Data());
  tree->Draw("dtl-fitTl1",cutA+cut1);
  //
  TString *strZ1 = toolkit.FitPlane(tree,"dz:ez",fstringZ1->Data(), cutA+cut1, chi2,npoints,fitParam,covMatrix);
  strZ1->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitZ1",strZ1->Data());
  tree->Draw("dz-fitZ1",cutA+cut1);
  cut1 =  "abs(dz-fitZ1)<2  && abs(dtl-fitTl1)<0.02"; 
}



void Fit2(){
  TString fstringTl2="";
  fstringTl2+="(side)++";
  fstringTl2+="(tl)++";
  fstringTl2+="(side*tl)++";
  //
  fstringTl2+="(sa)++";
  fstringTl2+="(side*sa)++";
  fstringTl2+="(ca)++";
  fstringTl2+="(side*ca)++";
  //
  fstringTl2+="(sa*tl)++";
  fstringTl2+="(side*sa*tl)++";
  fstringTl2+="(ca*tl)++";
  fstringTl2+="(side*ca*tl)++";
  //
  //
  //
  //
  TString fstringZ2="";
  fstringZ2+="(side)++";
  fstringZ2+="(-side*(250.0-side*z))++";
  fstringZ2+="((250.-side*z))++";
  // systematic  shift
  fstringZ2+="sa++";   
  fstringZ2+="ca++";
  fstringZ2+="side*sa++";
  fstringZ2+="side*ca++";
  //

  TString *strTl2 = toolkit.FitPlane(tree,"dtl:etl",fstringTl2->Data(), cutA+cut1+cut2, chi2,npoints,fitParam,covMatrix);
  strTl2->Tokenize("+")->Print();
  printf("Tl2: Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitTl2",strTl2->Data());
  tree->Draw("dtl-fitTl2",cutA+cut1);
  //
  TString *strZ2 = toolkit.FitPlane(tree,"dz:ez",fstringZ2->Data(), cutA+cut1+cut2, chi2,npoints,fitParam,covMatrix);
  strZ2->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitZ2",strZ2->Data());
  tree->Draw("dz-fitZ2",cutA+cut1);
  cut2 =  "abs(dz-fitZ2)<2  && abs(dtl-fitTl2)<0.02"; 
}
















TString fstringP3="";
fstringP3+="tl++";
fstringP3+="sa++";
fstringP3+="ca++";
fstringP3+="side*tl++";
fstringP3+="side*sa++";
fstringP3+="side*ca++";

TString *strP3 = toolkit.FitPlane(tree,"dtl:etl",fstringP3->Data(), cutA+cut38009, chi2,npoints,fitParam,covMatrix,0.9);
strP3->Tokenize("+")->Print();
printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
tree->SetAlias("fitP3",strP3->Data());
tree->Draw("dtl-fitP3",cut0+cut38009);



TString fstringP3="";
//
fstringP3+="side++";
fstringP3+="tl*(run<40000)++";
fstringP3+="tl*(abs(run-45600)<200)++";
fstringP3+="tl*(run>46000)++";
//
fstringP3+="sa++";
fstringP3+="ca++";
fstringP3+="sa*tl++";
fstringP3+="ca*tl++";
//
fstringP3+="side*tl++";
fstringP3+="side*sa++";
fstringP3+="side*ca++";
fstringP3+="side*sa*tl++";
fstringP3+="side*ca*tl++";



TString *strP3 = toolkit.FitPlane(tree,"dtl:etl",fstringP3->Data(), cutA, chi2,npoints,fitParam,covMatrix);
strP3->Tokenize("+")->Print();
printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
tree->SetAlias("fitP3",strP3->Data());
tree->Draw("dtl-fitP3",cut0);



TString fstringP1="";
//
fstringP1+="side++";
fstringP1+="dr*(run<40000)++";
fstringP1+="dr*(abs(run-45600)<200)++";
fstringP1+="dr*(run>46000)++";
fstringP1+="side*dr*(run<40000)++";
fstringP1+="side*dr*(abs(run-45600)<200)++";
fstringP1+="side*dr*(run>46000)++";
//
fstringP1+="tl*(run<40000)++";
fstringP1+="tl*(abs(run-45600)<200)++";
fstringP1+="tl*(run>46000)++";
fstringP1+="side*tl*(run<40000)++";
fstringP1+="side*tl*(abs(run-45600)<200)++";
fstringP1+="side*tl*(run>46000)++";
//
fstringP1+="sa++";
fstringP1+="ca++";
fstringP1+="sa*tl++";
fstringP1+="ca*tl++";
//
fstringP1+="side*sa++";
fstringP1+="side*ca++";
fstringP1+="side*sa*tl++";
fstringP1+="side*ca*tl++";

TCut cuttl="abs(dtl-fitP3)<0.02"

TString *strP1 = toolkit.FitPlane(tree,"dz",fstringP1->Data(),cuttl+cutA, chi2,npoints,fitParam,covMatrix);
strP1->Tokenize("+")->Print();
printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
tree->SetAlias("fitP1",strP1->Data());
tree->Draw("dz-fitP1",cuttl+cut0);

