/// \file driftITSTPC.C

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


void Init(){

  gSystem->Load("libSTAT");
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam, fitParam1,fitParam2;
  TVectorD errParam, errParam1,errParam2;
  TMatrixD covMatrix,covMatrix1,covMatrix2;
  
  Int_t npoints;
  //
  TFile f("driftitsTPC.root");
  TTree * tree = (TTree*)f.Get("Test");
  tree->SetAlias("side","(-1+(pTPC.fP[1]>0)*2)");     //side
  tree->SetAlias("z","(pTPC.fP[1]+0.0)");     //z position
  tree->SetAlias("dr","(1-abs(pTPC.fP[1])/250.)");    //norm drift length
  tree->SetAlias("tl","(pTPC.fP[3]+pITS.fP[3])*0.5"); //tan lampbda
  tree->SetAlias("sa","sin(pTPC.fAlpha+0.)");         //sin alpha
  tree->SetAlias("ca","cos(pTPC.fAlpha+0.)");         //cos alpha
  
  tree->SetAlias("dz","(pTPC.fP[1]-pITS.fP[1])");      //z delta
  tree->SetAlias("dy","(pTPC.fP[0]-pITS.fP[0])");      //z delta
  tree->SetAlias("dtl","(pTPC.fP[3]-pITS.fP[3])");     //delta tan lampbda
  tree->SetAlias("etl","sqrt(pITS.fC[9]+0.)");         //error tan lampbda
  tree->SetAlias("ez","sqrt(pITS.fC[2]+0.)");          //error z
}



void InitCuts(){
  TCut cut0("abs(dy)<1&&abs(pTPC.fP[1])>10");
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
}

//
// variable array
//
npoints=0;
TVectorD vside;
TVectorD vz;
TVectorD vdr;
TVectorD vtl;
TVectorD vsa;
TVectorD vca;
TVectorD vdz;
TVectorD vdtl;
TVectorD vez;
TVectorD vetl;

TVectorD errors;


void FillVar(){
  ///

  npoints =tree->Draw("side",cutA+cut1+cut2);
  vside.ResizeTo(npoints);  vside.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("z",cutA+cut1+cut2);
  vz.ResizeTo(npoints);  vz.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("dr",cutA+cut1+cut2);
  vdr.ResizeTo(npoints);  vdr.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("tl",cutA+cut1+cut2);
  vtl.ResizeTo(npoints);  vtl.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("sa",cutA+cut1+cut2);
  vsa.ResizeTo(npoints);  vsa.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("ca",cutA+cut1+cut2);
  vca.ResizeTo(npoints);  vca.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("dz",cutA+cut1+cut2);
  vdz.ResizeTo(npoints);  vdz.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("dtl",cutA+cut1+cut2);
  vdtl.ResizeTo(npoints);  vdtl.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("ez",cutA+cut1+cut2);
  vez.ResizeTo(npoints);  vez.SetElements(tree->GetV1());
  //
  npoints =tree->Draw("etl*0.5",cutA+cut1+cut2);
  vetl.ResizeTo(npoints);  vetl.SetElements(tree->GetV1()); 
  //
  npoints = tree->Draw("run","run>0");
  Int_t *np=new Int_t[npoints]; for (Int_t i=0;i<npoints;i++) np[i]=tree->GetV1()[i];
  Int_t *nsort=new Int_t[npoints];
  Int_t nruns = toolkit.Freq(npoints,np,nsort,kTRUE); 
}


void FitI1(){
  /// Independent fit 1

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

void Fit2I(){
  /// Independent fit

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


void Fit1(){
  ///  dz = zs-z = -s*vr *(z0-s*z)+s*dzt
  ///  dzs/dl = dz/dl +s*s*vr*dz/dl
  ///  d(dz/dl) = vr*dz/dl

  TLinearFitter fitter1(5, "hyp4");
  // parameters
  // 0 - P3   offset
  // 1 - P1   offset
  // 2 - vdr  factor
  //
  // 3 - side offset P3
  // 4 - side offset P1

  //  dz = zs-z = -s*vr *(z0-s*z)+s*dzt = [1] + [2]*(-s*(z0-s*z))+ [4]*s         
  //  d(dz/dl) = vr*dz/dl               = [0] + [2]*dz/dl        + [3]*s

  Double_t xxx[100];
  for (Int_t i=0;i<npoints;i++){
    for (Int_t jc =0;jc<10;jc++) xxx[jc]=0;
    // P1
    xxx[0] = 1;
    xxx[1] = -vside[i]*(250.0-vside[i]*vz[i]);
    xxx[3] = vside[i];
    fitter1.AddPoint(xxx,vdz[i],vez[i]);
    for (Int_t jc =0;jc<10;jc++) xxx[jc]=0;
    xxx[1] = vtl[i];
    xxx[2] = vside[i];
    fitter1.AddPoint(xxx,vdtl[i],vetl[i]);
  }
  fitter1.Eval();
  fitter1.GetParameters(fitParam1);  
  fitter1.GetErrors(errParam1);  
  //  dz = zs-z = -s*vr *(z0-s*z)+s*dzt = [1] + [2]*(-s*(z0-s*z))+ [4]*s         
  //  d(dz/dl) = vr*dz/dl               = [0] + [2]*dz/dl        + [3]*s

  TString fstrP1f1 ="";
  TString fstrP3f1 ="";
  fstrP1f1+=fitParam1[1];fstrP1f1+="-side*(250.0-side*z)*(";fstrP1f1+=fitParam1[2];
  fstrP1f1+=")+side*("; fstrP1f1+=fitParam1[4]; fstrP1f1+=")";
  //
  fstrP3f1+=fitParam1[0];fstrP3f1+="+tl*(";fstrP3f1+=fitParam1[2];
  fstrP3f1+=")+side*("; fstrP3f1+=fitParam1[3]; fstrP3f1+=")";
  //
  tree->SetAlias("fP1f1",fstrP1f1->Data());
  tree->SetAlias("fP3f1",fstrP3f1->Data());
  //
}



void Fit2(){
  ///  dz = zs-z = -s*vr *(z0-s*z)+s*dzt
  ///  dzs/dl = dz/dl +s*s*vr*dz/dl
  ///  d(dz/dl) = vr*dz/dl

  TLinearFitter fitter2(7, "hyp6");
  // parameters
  // 0 - P3   offset
  // 1 - P1   offset
  // 2 - vdr  factor
  //
  // 3 - side offset P3
  // 4 - side offset P1
  //
  // 5 - vdr difference -A side - C side
  // 6 - vdr gradient   -sin(alpha)*tl dependence 


  //  dz = zs-z = -s*vr *(z0-s*z)+s*dzt = [1] + [2]*(-s*(z0-s*z))+ [4]*s +[5]* (-s*s*(z0-s*z))      
  //  dz+= [6]*(-s*(z0-s*z))*sa
  //    
  //  d(dz/dl) = vr*dz/dl               = [0] + [2]*dz/dl        + [3]*s +[5]* s*dz/dl
  //   d(dz/dl)+= [6]*dz/dl*sa
  Double_t xxx[100];
  for (Int_t i=0;i<npoints;i++){
    for (Int_t jc =0;jc<10;jc++) xxx[jc]=0;
    // P1
    xxx[0] = 1;
    xxx[1] = -vside[i]*(250.0-vside[i]*vz[i]);
    xxx[3] = vside[i];
    xxx[4] = -vside[i]*vside[i]*(250.0-vside[i]*vz[i]);
    xxx[5] = -vside[i]*(250.0-vside[i]*vz[i])*vsa[i];
    fitter2.AddPoint(xxx,vdz[i],vez[i]);
    //P3
    for (Int_t jc =0;jc<10;jc++) xxx[jc]=0;
    xxx[1] = vtl[i];
    xxx[2] = vside[i];
    xxx[4] = vside[i]*vtl[i];
    xxx[5] = vtl[i]*vsa[i];
    fitter2.AddPoint(xxx,vdtl[i],vetl[i]);
  }
  fitter2.Eval();
  fitter2.GetParameters(fitParam2);  
  fitter2.GetErrors(errParam2);
  //  dz = zs-z = -s*vr *(z0-s*z)+s*dzt = [1] + [2]*(-s*(z0-s*z))+ [4]*s +[5]* (-s*s*(z0-s*z))      
  //  dz+= [6]*(-s*(z0-s*z))*sa
  //    
  //  d(dz/dl) = vr*dz/dl               = [0] + [2]*dz/dl        + [3]*s +[5]* s*dz/dl
  //   d(dz/dl)+= [6]*dz/dl*sa


  TString fstrP1f2 ="";
  TString fstrP3f2 ="";
  fstrP1f2+=fitParam2[1];fstrP1f2+="-side*(250.0-side*z)*(";fstrP1f2+=fitParam2[2];
  fstrP1f2+=")+side*("; fstrP1f2+=fitParam2[4];
  fstrP1f2+=")-side*side*(250.0-side*z)*("; fstrP1f2+=fitParam2[5];
  fstrP1f2+=")-side*(250.0-side*z)*sa*tl*(";fstrP1f2+=fitParam2[6];fstrP1f2+=")";
  //
  fstrP3f2+=fitParam2[0];fstrP3f2+="+tl*(";fstrP3f2+=fitParam2[2];
  fstrP3f2+=")+side*("; fstrP3f2+=fitParam2[3];
  fstrP3f2+=")+side*tl*("; fstrP3f2+=fitParam2[5];
  fstrP3f2+=")*tl*sa*("; fstrP3f2+=fitParam2[6];fstrP3f2+=")";
 
  //
  tree->SetAlias("fP1f2",fstrP1f2->Data());
  tree->SetAlias("fP3f2",fstrP3f2->Data());
  //
}






















TString fstringP3="";
fstringP3+="tl++";
fstringP3+="sa++";
fstringP3+="ca++";
fstringP3+="side*tl++";
fstringP3+="side*sa++";
fstringP3+="side*ca++";

TString *strP3 = toolkit.FitPlane(tree,"dtl:etl",fstringP3->Data(), cutA+cut38009, chi2,npoints,fitParam,covMatrix);
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

*/
