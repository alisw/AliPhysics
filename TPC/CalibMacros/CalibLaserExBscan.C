/// \file CalibLaserExBscan.C
///
/// ~~~
/// .x ~/rootlogon.C
/// .L $ALICE_ROOT/TPC/CalibMacros/CalibLaserExBscan.C
/// ~~~
///
/// 0. Make a calibration
/// 1. Make a laser scan list
///     e.g in TPC workscape
///
/// 2. Define a reference data  e.g:
/// ~~~
/// for a in `cat laserScan.txt`; do echo `pwd`/../mergerunMag0.list/laserMean.root; done >laserScanRef.txt
/// ~~~
///
/// ~~~{.cpp}
/// Init();         // init chain
/// MakeAliases();  // make alaises for variables
/// ~~~

TCut  cutE="eY.fElements<0.02&&Rr.eY.fElements<0.02";         // error cut 200 microns 
TCut  cutN="nCl.fElements>10&&Rr.nCl.fElements>10";           // number of clusters cut
TCut  cutEd="(abs(Rr.Y.fElements)-Rr.X.fElements*0.155<-1)";   // edge cut
TCut  cutB="LTr.fVecLX.fElements>0";                          // local x cut
TCut  cutR="(LTr.fVecLX.fElements%3)==1&&((abs(LTr.fVecLZ.fElements)+abs(LTr.fVecLY.fElements))%3)==1";   // cut - skip points
TCut  cutA=cutE+cutN+cutEd+cutB;
//
TCut  cutAM5  =cutA+"LTr.fP[1]>0&&abs(bz+5)<0.1";
TCut  cutAP2  =cutA+"LTr.fP[1]>0&&abs(bz-2)<0.1";
TCut  cutAP5  =cutA+"LTr.fP[1]>0&&abs(bz-5)<0.1";
//
TCut  cutCM5  =cutA+"LTr.fP[1]<0&&abs(bz+5)<0.1";
TCut  cutCP2  =cutA+"LTr.fP[1]<0&&abs(bz-2)<0.1";
TCut  cutCP5  =cutA+"LTr.fP[1]<0&&abs(bz-5)<0.1";


TChain * chain =0;

TVectorD fitParamB[6];      // magnetic fit param
TVectorD fitParamBT[6];     // + electric field tilt
TVectorD fitParamBT0[6];    // + electric field close to ROC              
TVectorD fitParamA[6];      // + electric field rotation
//
TMatrixD covarB[6];         // covariance matrices
TMatrixD covarBT[6];        //
TMatrixD covarBT0[6];       //
TMatrixD covarA[6];         //
//
TMatrixD *errB[6];         // covariance matrices
TMatrixD *errBT[6];        //
TMatrixD *errBT0[6];       //
TMatrixD *errA[6];         //
//
Double_t chi2B[6];          // chi2
Double_t chi2BT[6];         // chi2
Double_t chi2BT0[6];        //
Double_t chi2A[6];          //


void Init(){
  ///

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib"); 
  gSystem->Load("libSTAT");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool;
  chain = tool.MakeChainRandom("laserScan.txt","Mean",0,10200);
  chain->Lookup();
  chainRef = tool.MakeChain("laserScanRef.txt","Mean",0,10200);
  chainRef->Lookup();
  chain->AddFriend(chainRef,"Rr");
}



void MakeAliases(){
  /// shortcuts for variables to fit
  ///
  /// bserved distrotions

  chain->SetAlias("dy","(dY.fElements-Rr.dY.fElements)");              // y ~ (r-phi) distortion
  chain->SetAlias("ctany","(dY.fElements-Rr.dY.fElements)/(250*dr)");  // mean distortion (angle) over drift length
  //
  chain->SetAlias("dr","(250.-abs(LTr.fP[1]))");                       // drift length 0-250
  chain->SetAlias("dr1","(1.-abs(LTr.fP[1]/250.))");                       // drift length 0-250
  chain->SetAlias("phiL","(LTr.fVecLY.fElements/LTr.fVecLX.fElements)");       // tann of local phi
  chain->SetAlias("Esign","(1.-LTr.fSide*2.)");          // E field sign:    +1 for  A side : -1 for c side 
  //
  // Br and Brfi  
  //
  chain->SetAlias("ablx","(iblx.fElements/(ibz.fElements))");
  chain->SetAlias("ably","(ibly.fElements/(ibz.fElements))");
  chain->SetAlias("abr","(ibr.fElements/(ibz.fElements))");
  chain->SetAlias("abrf","(ibrphi.fElements/(ibz.fElements))");


  chain->SetAlias("r","(R.fElements-165.5)/80.3");          // local X - 0 at middle   -1 first +1 last padrow
  chain->SetAlias("ky","(kY.fElements)");                   // local inclination angle
  chain->SetAlias("sec","LTr.fVecSec.fElements");           // sector number
  chain->SetAlias("ca","cos(LTr.fVecPhi.fElements+0.)");    // cos of global phi position
  chain->SetAlias("sa","sin(LTr.fVecPhi.fElements+0.)");    // sin of global phi position  
}


TMatrixD * MakeErrVector(TMatrixD & mat){
  /// get error vector

  Int_t    nrows=mat.GetNrows();
  TMatrixD *err = new TMatrixD(nrows,1);
  for (Int_t i=0; i<nrows;i++) (*err)(i,0)=TMath::Sqrt(mat(i,i));
  return err;
}


void MakeFit(Int_t i, TCut cutI, TString  aName){
  ///

  Int_t  ntracks=3000000;
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  //
  TString  fstringB="";   // magnetic part 
  TString  fstringT="";    // E file dtilting part
  TString  fstring0="";   // ROC      part 
  TString  fstringL="";   // linear part
  //
  //  Magnetic field map part
  // 
  //                           // 0
  fstringB+="ablx*dr++";       // 1
  fstringB+="ablx*ky*dr++";    // 2
  fstringB+="ably*dr++";       // 3
  fstringB+="ably*ky*dr++";    // 4
  //
  // Electric field tilting part
  // 
  //
  fstringT+="(dr1)++";             // 5     - ex
  fstringT+="(dr1)*ky++";          // 6     - ex
  fstringT+="(dr1)*phiL++";        // 7     - ey
  fstringT+="(dr1)*phiL*ky++";     // 8     - ey
  //
  // E field close to the ROC - radius independent part  
  //                          
  // 
  fstring0+="ca++";           // 9
  fstring0+="sa++";           // 10
  fstring0+="ky++";           // 11 
  fstring0+="ky*ca++";        // 12
  fstring0+="ky*sa++";        // 13
  //
  // E field close to the ROC - radius dependent part
  //
  fstring0+="r++";              // 14
  fstring0+="ca*r++";           // 15
  fstring0+="sa*r++";           // 16
  fstring0+="ky*r++";           // 17
  fstring0+="ky*ca*r++";        // 18
  fstring0+="ky*sa*r++";        // 19
  //
  // E field rotation - drift length linear part
  //  
  fstringL+="ca*dr1++";         //20
  fstringL+="sa*dr1++";         //21
  fstringL+="ky*ca*dr1++";      //22 
  fstringL+="ky*sa*dr1++";      //23
  //
  TString fstringBT   = fstringB +fstringT;
  TString fstringBT0  = fstringBT+fstring0; 
  TString fstringA    = fstringBT0+fstringL; 
  //
  //
  TCut cutF=cutI;
  TString * strFit = 0;
  //
  //
  cutF=cutR+cutI;
  //
  strFit = TStatToolkit::FitPlane(chain,"dy:0.1", fstringB.Data(),cutF, chi2,npoints,fitParamB[i],covarB[i],1,0, ntracks);
  chain->SetAlias(aName+"_FB",strFit->Data());
  cutF=cutR+cutI + Form("abs(dy-%s)<0.2",(aName+"_FB").Data());
  //
  //
  //
  strFit = TStatToolkit::FitPlane(chain,"dy:0.1", fstringB.Data(),cutF, chi2,npoints,fitParamB[i],covarB[i],1,0, ntracks);
  chain->SetAlias(aName+"_FB",strFit->Data()); 
  chi2B[i]=TMath::Sqrt(chi2/npoints);
  printf("B fit sqrt(Chi2/npoints)=%f\n",chi2B[i]);
  //
  strFit = TStatToolkit::FitPlane(chain,"dy:0.1", fstringBT.Data(),cutF, chi2,npoints,fitParamBT[i],covarBT[i],1,0, ntracks);
  chain->SetAlias(aName+"_FBT",strFit->Data()); 
  chi2BT[i]=TMath::Sqrt(chi2/npoints);
  printf("BT fit sqrt(Chi2/npoints)=%f\n",chi2BT[i]);
  //
  strFit = TStatToolkit::FitPlane(chain,"dy:0.1", fstringBT0.Data(),cutF, chi2,npoints,fitParamBT0[i],covarBT0[i],1,0, ntracks);
  chain->SetAlias(aName+"_FBT0",strFit->Data()); 
  chi2BT0[i]=TMath::Sqrt(chi2/npoints);
  printf("BT0 Fit: sqrt(Chi2/npoints)=%f\n",chi2BT0[i]);
  //
  //
  strFit = TStatToolkit::FitPlane(chain,"dy:0.1", fstringA.Data(),cutF, chi2,npoints,fitParamA[i],covarA[i],1, 0, ntracks);
  chain->SetAlias(aName+"_FA",strFit->Data()); 
  chi2A[i]=TMath::Sqrt(chi2/npoints);
  printf("A fit: sqrt(Chi2/npoints)=%f\n",chi2A[i]);
  errB[i]   =(MakeErrVector(covarB[i]));
  errBT[i]  =(MakeErrVector(covarBT[i]));
  errBT0[i] =(MakeErrVector(covarBT0[i]));
  errA[i]   =(MakeErrVector(covarA[i]));
}





void MakeFitDy(){
  ///


  MakeFit(0,cutAM5,"dyAM5");
  MakeFit(1,cutAP2,"dyAP2");
  MakeFit(2,cutAP5,"dyAP5");
  //
  MakeFit(3,cutCM5,"dyCM5");
  MakeFit(4,cutCP2,"dyCP2");
  //  MakeFit(5,cutCP5,"dyAP5");
  DumpFit();
}


void DrawPhi(){
  ///

  chain->Draw("(dy-dyAM5_FA):LTr.fVecPhi.fElements>>hisAM5(60,-3.14,3.14,100,-0.2,0.2)",cutAM5,"");
  chain->Draw("(dy-dyAP5_FA):LTr.fVecPhi.fElements>>hisAP5(60,-3.14,3.14,100,-0.2,0.2)",cutAP5,"");
  chain->Draw("(dy-dyAP2_FA):LTr.fVecPhi.fElements>>hisAP2(60,-3.14,3.14,100,-0.2,0.2)",cutAP2,"");
  hisAM5->FitSlicesY(0,0,-1,20);
  hisAP5->FitSlicesY(0,0,-1,20);
  hisAP2->FitSlicesY(0,0,-1,20);
  hisAM5_1->SetMinimum(-0.2);
  hisAM5_1->SetMaximum(0.2);
  hisAM5_1->SetMarkerStyle(20);
  hisAP5_1->SetMarkerStyle(21);
  hisAP2_1->SetMarkerStyle(22);
  hisAM5_1->SetMarkerColor(1);
  hisAP5_1->SetMarkerColor(2);
  hisAP2_1->SetMarkerColor(4);
  hisAM5_1->Draw();
  hisAP5_1->Draw("same");
  hisAP2_1->Draw("same");
}

void MakeGraphs(){
  Double_t bz[3] ={-5,2,5};
  Double_t p0A[3]={fitParam[0][0],fitParam[1][0],fitParam[2][0]};
  Double_t p1A[3]={fitParam[0][1],fitParam[1][1],fitParam[2][1]};
  Double_t p2A[3]={fitParam[0][2],fitParam[1][2],fitParam[2][2]};
  Double_t p3A[3]={fitParam[0][3],fitParam[1][3],fitParam[2][3]};
  Double_t p4A[3]={fitParam[0][4],fitParam[1][4],fitParam[2][4]};
  Double_t p5A[3]={fitParam[0][5],fitParam[1][5],fitParam[2][5]};
  Double_t p6A[3]={fitParam[0][6],fitParam[1][6],fitParam[2][6]};
  Double_t p7A[3]={fitParam[0][7],fitParam[1][7],fitParam[2][7]};
  Double_t p8A[3]={fitParam[0][8],fitParam[1][8],fitParam[2][8]};
  Double_t p9A[3]={fitParam[0][9],fitParam[1][9],fitParam[2][9]};
  TGraph grA0(3,bz,p0A);
  TGraph grA1(3,bz,p1A);
  TGraph grA2(3,bz,p2A);
  TGraph grA3(3,bz,p3A);
  TGraph grA4(3,bz,p4A);
  TGraph grA5(3,bz,p5A);
  TGraph grA6(3,bz,p6A);
  TGraph grA7(3,bz,p7A);
  TGraph grA8(3,bz,p8A);
  TGraph grA9(3,bz,p9A);
  TF1 f1("f1","[0]*x/(1+([0]*x)^2)");
  TF1 f2("f2","([0]*x)^2/(1+([0]*x)^2)");
  TMatrixD matB(5,5);
  TMatrixD matBT(6,5);
  TMatrixD matBT0(6,5);
  for (Int_t i=0; i<5;i++) for (Int_t j=0; j<5;j++) matB[i][j]=fitParamB[j][i];
  for (Int_t i=0; i<5;i++) for (Int_t j=0; j<5;j++) matBT[i][j]=fitParamBT[j][i];
  for (Int_t i=0; i<6;i++) for (Int_t j=0; j<5;j++) matB0[i][j]=fitParamBT0[j][i];
}


void DumpFit(){
  ///

  TTreeSRedirector *pcstream = new TTreeSRedirector("exbFits.root");

  for (Int_t i=0; i<5;i++){
    Double_t bz=0;
    Double_t side=0;
    if (i==0) { bz=-5; side=1;}
    if (i==1) { bz=2;  side=1;}
    if (i==2) { bz=5;  side=1;}
    if (i==3) { bz=-5; side=-1;}
    if (i==4) { bz=2;  side=-1;}
    (*pcstream)<<"fit"<<
      "bz="<<bz<<
      "side="<<side<<
      "pb.="<<&fitParamB[i]<<
      "pbT.="<<&fitParamBT[i]<<
      "pbT0.="<<&fitParamBT0[i]<<
      "pbA.="<<&fitParamA[i]<<
      "eb.="<<errB[i]<<
      "ebT.="<<errBT[i]<<
      "ebT0.="<<errBT0[i]<<
      "ebA.="<<errA[i]<<
      "\n";
  }
  delete pcstream;
}

void MakeFitPic(){
  TFile f("exbFits.root");

}
