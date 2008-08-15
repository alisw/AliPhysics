void CalibExB(){
  // 
  // Macro to test ExB correction versus laser B filed scan data
  //
  // Before running it be sure the file laserScan.root exist
  // 



  //AliMagF* field = new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 2);
  AliTracker::SetFieldMap(field,1);
  AliTPCExBFirst *exbfirst = new  AliTPCExBFirst(field,2.64000000000000000e+04,50,50,50);
  AliTPCExB::SetInstance(exbfirst);
  //
  // Create this file from scan - See AliTPCCalibLaser.C
  //
  TFile fscan("laserScan.root");
  TTree * treeT = (TTree*)fscan.Get("Mean");
  gSystem->Load("libSTAT.so");
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam;
  TMatrixD covMatrix;
  Int_t npoints;
  
  TCut cutF0("abs(gphi1-(pphi0+pphi1*bz))<0.05");
  TCut cutF1("abs(gphiP1-(pphiP0+pphiP1*bz))<0.0005");
  TCut cutN("entries>2");
  
  TCut cutA = cutF0+cutF1+cutN;
  
  treeT->SetAlias("side","(-1+(LTr.fP[2]>0)*2)");       // side
  treeT->SetAlias("dr","(abs(LTr.fP[1]/250.))");
  treeT->SetAlias("sa","sin(atan2(lx1+0.0,lx0+0.0))");
  treeT->SetAlias("ca","cos(atan2(lx1+0.0,lx0+0.0))");
  treeT->SetAlias("ta","tan(asin(LTr.fP[2]+0.0))");

  TString fstring="";
  //
  fstring+="((dr)^3-1)*bz++";           //1
  fstring+="((dr)^3-1)*bz*ta++";        //2
  fstring+="((dr)^3-1)*bz*sa++";        //3
  fstring+="((dr)^3-1)*bz*sa*ta++";     //4
  fstring+="((dr)^3-1)*bz*ca++";        //5
  fstring+="((dr)^3-1)*bz*ca*ta++";     //6

  fstring+="side*((dr)^3-1)*bz++";           //1
  fstring+="side*((dr)^3-1)*bz*ta++";        //2
  fstring+="side*((dr)^3-1)*bz*sa++";        //3
  fstring+="side*((dr)^3-1)*bz*sa*ta++";     //4
  fstring+="side*((dr)^3-1)*bz*ca++";        //5
  fstring+="side*((dr)^3-1)*bz*ca*ta++";     //6

  fstring+="side*((dr)^1-1)*bz++";           //7
  fstring+="side*((dr)^1-1)*bz*ta++";        //8 
  fstring+="side*((dr)^1-1)*bz*sa++";        //9
  fstring+="side*((dr)^1-1)*bz*sa*ta++";     //10
  fstring+="side*((dr)^1-1)*bz*ca++";        //11
  fstring+="side*((dr)^1-1)*bz*ca*ta++";     //12




  TString *strq0 = toolkit.FitPlane(treeT,"gphi1-pphi0",fstring->Data(), cutA, chi2,npoints,fitParam,covMatrix);
  strq0->Tokenize("+")->Print();
  treeT->SetAlias("fit",strq0->Data());
  
  
  TString fstringeb="";
  //
  fstringeb+="AliTPCExB::GetDrphi(260,atan2(lx1,lx0),LTr.fP[1])++";       //1
  fstringeb+="AliTPCExB::GetDr(260,atan2(lx1,lx0),LTr.fP[1])*ta++";       //1
  // fstringeb+="bz*bz*AliTPCExB::GetDrphi(260,atan2(lx1,lx0),LTr.fP[1])++";              //1
  //fstringeb+="bz*bz*AliTPCExB::GetDr(260,atan2(lx1,lx0),LTr.fP[1])*ta++";       //1
  
  //  fstringeb+="side*((dr)^1-1)*bz++";           //3
  //  fstringeb+="side*((dr)^1-1)*bz*ta++";        //4 
  //  fstringeb+="side*((dr)^1-1)*bz*sa++";        //7
  //  fstringeb+="side*((dr)^1-1)*bz*sa*ta++";     //8   
  //  fstringeb+="side*((dr)^1-1)*bz*ca++";        //11
  //  fstringeb+="side*((dr)^1-1)*bz*ca*ta++";     //12
  
  
  
  TString *strExB = toolkit.FitPlane(treeT,"gphi1-pphi0",fstringeb->Data(), "abs(bz+0.4)<0.05"+cutA, chi2,npoints,fitParam,covMatrix);
  strExB->Tokenize("+")->Print();
  treeT->SetAlias("fitEB",strExB->Data());
}
