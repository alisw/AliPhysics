void Init(){
  //
  // Initialize
  //
  AliMagF* field = new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG);
  AliTPCExB::RegisterField(0,field);
  AliMagF* fieldC0 = new AliMagF("Maps","Maps", 2, 1, 1, 10., AliMagF::k5kG);
  AliTPCExB::RegisterField(1,fieldC0);
  AliMagF* fieldC1 = new AliMagF("Maps","Maps", 2, 1, 1, 10., AliMagF::k5kG,kTRUE,"$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  AliTPCExB::RegisterField(2,fieldC1);

  gSystem->Load("libSTAT.so");
  AliTPCExBFirst *exbfirst1  = new  AliTPCExBFirst(fieldC1,0.88*2.6400e+04,50,50,50);
  AliTPCExB::SetInstance(exbfirst1);

}

void FitField(){
  //
  //
  //
  AliTPCExB * exb =  AliTPCExB::Instance();
  exb->TestExB("field.root");
  TFile f("field.root");
  TTree * tree = (TTree*)f.Get("positions");
  //
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam;
  TMatrixD covMatrix;
  Int_t npoints;
  //
  // SetAliases
  //
  tree->SetAlias("sa","sin(phi+0.0)");
  tree->SetAlias("ca","cos(phi+0.0)");
  tree->SetAlias("sa2","sin(phi*2+0.0)");
  tree->SetAlias("ca2","cos(phi*2+0.0)");
  tree->SetAlias("zn","(x2/250.)");
  tree->SetAlias("rn","(r/250.)");
    
  TString fstringSym="";
  //  
  fstringSym+="zn++";
  fstringSym+="rn++";
  fstringSym+="zn*rn++";
  fstringSym+="zn*zn++";
  fstringSym+="zn*zn*rn++";
  fstringSym+="zn*rn*rn++";
  //
  fstringSym+="sa++";
  fstringSym+="ca++";  
  fstringSym+="ca2++";
  fstringSym+="sa2++";
  fstringSym+="ca*zn++";
  fstringSym+="sa*zn++";
  fstringSym+="ca2*zn++";
  fstringSym+="sa2*zn++";
  fstringSym+="ca*zn*zn++";
  fstringSym+="sa*zn*zn++";
  fstringSym+="ca*zn*rn++";
  fstringSym+="sa*zn*rn++";

  // BrBz
  TString *strBrBz = toolkit.FitPlane(tree,"br/bz",fstringSym->Data(), "abs(r)<250&&abs(r)>90", chi2,npoints,fitParam,covMatrix);
  strBrBz->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitBrBz",strBrBz->Data());
  exb->fMatBrBz = new TVectorD(fitParam);
  // BrfiBz
  TString *strBrfiBz = toolkit.FitPlane(tree,"brfi/bz",fstringSym->Data(), "abs(r)<250&&abs(r)>90", chi2,npoints,fitParam,covMatrix);
  strBrfiBz->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitBrfiBz",strBrfiBz->Data());
  exb->fMatBrfiBz = new TVectorD(fitParam);

  
  //
  // BR integral parameterization
  //
  TString *strBRI0 = toolkit.FitPlane(tree,"bri",fstringSym->Data(), "abs(r)<250&&abs(r)>90&&x2>0&&x2<250", chi2,npoints,fitParam,covMatrix);
  strBRI0->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitBRI0",strBRI0->Data());
  exb->fMatBrBzI0 = new TVectorD(fitParam);

  TString *strBRI1 = toolkit.FitPlane(tree,"bri",fstringSym->Data(), "abs(r)<250&&abs(r)>90&&x2<-10&&x2>-250", chi2,npoints,fitParam,covMatrix);
  strBRI1->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitBRI1",strBRI1->Data());
  exb->fMatBrBzI1 = new TVectorD(fitParam);
  //

  TString *strBRFII0 = toolkit.FitPlane(tree,"brfii",fstringSym->Data(), "abs(r)<250&&abs(r)>90&&x2>0&&x2<250", chi2,npoints,fitParam,covMatrix);
  strBRFII0->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitBRFII0",strBRFII0->Data());
  exb->fMatBrfiBzI0 = new TVectorD(fitParam);
 
  TString *strBRFII1 = toolkit.FitPlane(tree,"brfii",fstringSym->Data(), "abs(r)<250&&abs(r)>90&&x2<-10&&x2>-240", chi2,npoints,fitParam,covMatrix);
  strBRFII1->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  tree->SetAlias("fitBRFII1",strBRFII1->Data());
  exb->fMatBrfiBzI1 = new TVectorD(fitParam);

  TFile * fout = new TFile("bfit.root","recreate");
  exb->Write("bfit");
  fout->Close();

}






void CalibExB(){
  // 
  // Macro to test ExB correction versus laser B filed scan data
  //
  // Before running it be sure the file laserScan.root exist
  //
  // Create this file from scan - See AliTPCCalibLaser.C
  //
  TFile fscan("laserScan.root");
  TTree * treeT = (TTree*)fscan.Get("Mean");
  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParam;
  TMatrixD covMatrix;
  Int_t npoints;
  
  TCut cutF0("abs(gphi1-(pphi0+pphi1*bz))<0.05");
  TCut cutF1("abs(gphiP1-(pphiP0+pphiP1*bz))<0.0005");
  TCut cutN("entries>2");
  
  TCut cutA = cutF0+cutF1+cutN;
  
  treeT->SetAlias("side","(-1+(LTr.fP[1]>0)*2)");       // side
  treeT->SetAlias("dr","(abs(LTr.fP[1]/250.))");
  treeT->SetAlias("sa","sin(atan2(lx1+0.0,lx0+0.0))");
  treeT->SetAlias("ca","cos(atan2(lx1+0.0,lx0+0.0))");
  treeT->SetAlias("ta","tan(asin(LTr.fP[2]+0.0))");

  TString fstring="";
  //
  fstring+="((dr)^3-1)*bz++";             //1
  fstring+="((dr)^3-1)*bz*sa++";        //3
  fstring+="((dr)^3-1)*bz*ca++";        //5
  //
  fstring+="((dr)^1-1)*bz++";           //1
  fstring+="((dr)^1-1)*bz*sa++";        //3
  fstring+="((dr)^1-1)*bz*ca++";        //5
  //
  fstring+="side*((dr)^3-1)*bz++";             //1
  fstring+="side*((dr)^3-1)*bz*sa++";        //3
  fstring+="side*((dr)^3-1)*bz*ca++";        //5
  //
  fstring+="side*((dr)^1-1)*bz++";           //7
  fstring+="side*((dr)^1-1)*bz*sa++";        //9
  fstring+="side*((dr)^1-1)*bz*ca++";        //11
  //
  fstring+="((dr)^3-1)*bz^2*ta++";        //2
  fstring+="((dr)^3-1)*bz^2*sa*ta++";     //4
  fstring+="((dr)^3-1)*bz^2*ca*ta++";     //6

  fstring+="((dr)^1-1)*bz^2*ta++";        //2
  fstring+="((dr)^1-1)*bz^2*sa*ta++";     //4
  fstring+="((dr)^1-1)*bz^2*ca*ta++";     //6
  //  
  fstring+="side*((dr)^1-1)*bz^2*ta++";        //8 
  fstring+="side*((dr)^1-1)*bz^2*sa*ta++";     //10
  fstring+="side*((dr)^1-1)*bz^2*ca*ta++";     //12




  TString *strq0 = toolkit.FitPlane(treeT,"gphi1-pphi0",fstring->Data(), cutA, chi2,npoints,fitParam,covMatrix);
  strq0->Tokenize("+")->Print();
  treeT->SetAlias("fit",strq0->Data());
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));


  
  treeT->SetAlias("bir1",  "-AliTPCExB::GetBrI(254,atan2(lx1,lx0+0),LTr.fP[1],2)");
  treeT->SetAlias("birfi1","-AliTPCExB::GetBrfiI(254,atan2(lx1,lx0+0),LTr.fP[1],2)");

  treeT->SetAlias("bir0",  "AliTPCExB::GetBrI(254,atan2(lx1,lx0+0),LTr.fP[1],2)");
  treeT->SetAlias("birfi0","AliTPCExB::GetBrfiI(254,atan2(lx1,lx0+0),LTr.fP[1],2)");

  treeT->SetAlias("fbz00", "(bir0+birfi0*ta)");
  treeT->SetAlias("fbz02", "(birfi0+bir0*ta)");
  treeT->SetAlias("fbz10", "(bir1+birfi1*ta)");
  treeT->SetAlias("fbz12", "(birfi1+bir1*ta)");
  //
  treeT->SetAlias("fbz0", "((fSide==0)*fbz00+(fSide==1)*fbz10)");
  treeT->SetAlias("fbz2", "((fSide==0)*fbz02+(fSide==1)*fbz12)");

  //
  TString fstringeb="";
  //
  fstringeb+="bz*fbz0++";                //1
  fstringeb+="sign(bz)*bz^2*fbz2++";     //2
  fstringeb+="((dr)-1)*bz*sa++";        //9
  fstringeb+="side*((dr)-1)*bz*sa++";        //9


  fstringeb+="((dr)-1)*bz*ca++";        //9
  fstringeb+="side*((dr)^3-1)*bz*sa++";        //9
  fstringeb+="side*((dr)^3-1)*bz++";           //7

  fstringeb+="((dr)^3-1)*bz++";           //7
  fstringeb+="side*((dr)^3-1)*bz*ca++";        //11
  
  //  
  TString *strExB = toolkit.FitPlane(treeT,"gphi1-pphi0",fstringeb->Data(), "abs(gphi1-pphi0-fit)<0.06&&abs(bz)>0.1"+cutA, chi2,npoints,fitParam,covMatrix);
  strExB->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  treeT->SetAlias("fitEB",strExB->Data());



  TString fstringeb="";
  //
  fstringeb+="AliTPCExB::GetDrphi(254,atan2(lx1,lx0),LTr.fP[1],bz*10)++";        //1
  fstringeb+="AliTPCExB::GetDr(254,atan2(lx1,lx0),LTr.fP[1],-bz*10)*ta++";        //2
  //
  fstringeb+="side*((dr)^1-1)*bz++";        //9
  fstringeb+="side*((dr)^1-1)*bz*sa++";        //9
  //
  fstringeb+="side*((dr)^1-1)*bz*ta++";        //9
  fstringeb+="side*((dr)^1-1)*bz*sa*ta++";        //9
  fstringeb+="side*((dr)^1-1)*bz*ca*ta++";        //9

  TString *strExB = toolkit.FitPlane(treeT,"fit",fstringeb->Data(), "abs(gphi1-pphi0-fit)<0.06&&abs(bz)>0.1"+cutA, chi2,npoints,fitParam,covMatrix);
  strExB->Tokenize("+")->Print();
  printf("Chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  treeT->SetAlias("fitEB",strExB->Data());



 



  {
    Float_t otcor0    =0.7;
    Float_t omegatau0 =0.35;
    Float_t rms0      =10;
    for (Float_t otcor=0.6; otcor<1.1;otcor+=0.1)
      for (Float_t ometau=0.4; ometau<0.45;ometau+=0.01){
	char hname[100];
	sprintf(hname,"hist_ometau_%f(50,-0.10,0.10)",ometau);
	char expr[100];
	sprintf(expr,"gphi1-pphi0-bz*fbz0*%f-sign(bz)*(%f*%f*bz)^2*fbz2>>%s",ometau,ometau,otcor,hname);
	treeT->Draw(expr,"abs(gphi1-pphi0-fit)<0.06&&abs(bz)>0.1"+cutA);    
	printf("Ometau=%f\tCor=%f\tRMS=%f\n",ometau,otcor,treeT->GetHistogram()->GetRMS());
	if (rms0>treeT->GetHistogram()->GetRMS()){
	  otcor0=otcor;
	  omegatau0=ometau;
	  rms0=treeT->GetHistogram()->GetRMS();
	}
      }
  }
 
}




void MakePic(){
  //
  //
  //

}

