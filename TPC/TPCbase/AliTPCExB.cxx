#include "AliTPCExB.h"
#include "TMath.h"
//#include "TTreeStream.h"
#include "AliMagF.h"
#include "TLinearFitter.h"
#include "AliTPCcalibDB.h"


//
// Abstract class for ExB effect parameterization
// 
//
// 
// The ExB correction map is stored in the calib DB
// The lookup can be dumped to the tree:
/*

   //
  char *storage = "local://OCDBres"
  Int_t RunNumber=0;
  AliCDBManager::Instance()->SetDefaultStorage(storage);
  AliCDBManager::Instance()->SetRun(RunNumber) 
  AliTPCExBFirst * exb = AliTPCcalibDB::Instance()->GetExB();
  //
  // See example macro $ALICE_ROOT/TPC/macros/AliTPCExBdraw.C 
  //
  .L $ALICE_ROOT/TPC/macros/AliTPCExBdraw.C 
  Draw(0)




*/

AliTPCExB* AliTPCExB::fgInstance = 0;

TObjArray   AliTPCExB::fgArray;



ClassImp(AliTPCExB)

AliTPCExB::AliTPCExB():
  TObject(),
  fMatBrBz(0),       //param matrix Br/Bz
  fMatBrfiBz(0),     //param matrix Br/Bz
  fMatBrBzI0(0),     //param matrix Br/Bz integral  z>0 
  fMatBrBzI1(0),     //param matrix Br/Bz integral  z<0 
  fMatBrfiBzI0(0),   //param matrix Br/Bz integral  z>0 
  fMatBrfiBzI1(0)    //param matrix Br/Bz integral  z<0
{
  //
  // default constructor
  //
}

AliTPCExB::AliTPCExB(const AliTPCExB& exb):
  TObject(exb),
  fMatBrBz(new TVectorD(*(exb.fMatBrBz))),       //param matrix Br/Bz
  fMatBrfiBz(new TVectorD(*(exb.fMatBrfiBz))),     //param matrix Br/Bz
  fMatBrBzI0(new TVectorD(*(exb.fMatBrBzI0))),     //param matrix Br/Bz integral  z>0 
  fMatBrBzI1(new TVectorD(*(exb.fMatBrBzI1))),     //param matrix Br/Bz integral  z<0 
  fMatBrfiBzI0(new TVectorD(*(exb.fMatBrfiBzI0))),   //param matrix Br/Bz integral  z>0 
  fMatBrfiBzI1(new TVectorD(*(exb.fMatBrfiBzI1)))    //param matrix Br/Bz integral  z<0
{
  //
  // copy constructor
  //
}

AliTPCExB& AliTPCExB::operator=(const AliTPCExB &/*exb*/)
{
  //
  // Dummy  assignment
  //
  return *this;
}




void AliTPCExB::TestExB(const char* fileName) {
  //
  // Test ExB  - 
  // Dump the filed and corrections to the tree in file fileName  
  //
  // 
  TTreeSRedirector ts(fileName);
  Double_t x[3];
  for (x[0]=-250.;x[0]<=250.;x[0]+=10.)
    for (x[1]=-250.;x[1]<=250.;x[1]+=10.)
      for (x[2]=-250.;x[2]<=250.;x[2]+=20.) {
	Double_t r=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
	if (r<20) continue;
	if (r>260) continue;
	Double_t z = x[2];
	Double_t d[3];
	Correct(x,d);
	Double_t rd=TMath::Sqrt(d[0]*d[0]+d[1]*d[1]);
	Double_t dr=r-rd;
	Double_t phi=TMath::ATan2(x[1],x[0]);
	Double_t phid=TMath::ATan2(d[1],d[0]);
	Double_t dphi=phi-phid;
	if (dphi<0.) dphi+=TMath::TwoPi();
	if (dphi>TMath::Pi()) dphi=TMath::TwoPi()-dphi;
	Double_t drphi=r*dphi;
	Double_t dx=x[0]-d[0];
	Double_t dy=x[1]-d[1];
	Double_t dz=x[2]-d[2];
	//
	Double_t bx    = GetBx(r,phi,z,0);
	Double_t by    = GetBy(r,phi,z,0);
	Double_t bz    = GetBz(r,phi,z,0);
	Double_t br    = GetBr(r,phi,z,0);
	Double_t brfi  = GetBrfi(r,phi,z,0);
	//
	Double_t bxi    = GetBxI(r,phi,z,0);
	Double_t byi    = GetByI(r,phi,z,0);
	Double_t bzi    = GetBzI(r,phi,z,0);
	Double_t bri    = GetBrI(r,phi,z,0);
	Double_t brfii  = GetBrfiI(r,phi,z,0);

	ts<<"positions"<<
	  "x0="<<x[0]<<
	  "x1="<<x[1]<<
	  "x2="<<x[2]<<
	  "dx="<<dx<<
	  "dy="<<dy<<
	  "dz="<<dz<<
	  "r="<<r<<
	  "phi="<<phi<<
	  "dr="<<dr<<
	  "drphi="<<drphi<<
	  //
	  // B-Field
	  //
	  "bx="<<bx<<
	  "by="<<by<<
	  "bz="<<bz<<
	  "br="<< br<<
	  "brfi="<<brfi<<
	  // B-field integ
	  "bxi="<<bxi<<
	  "byi="<<byi<<
	  "bzi="<<bzi<<
	  "bri="<< bri<<
	  "brfii="<<brfii<<
	  "\n";
      }
}



Double_t AliTPCExB::GetDr(Double_t r, Double_t phi, Double_t z, Double_t bz){
  //
  // Static function
  // Posibble to us it for visualization 
  // 
  //
  AliTPCExB *exb = Instance();
  if (!exb) exb = AliTPCcalibDB::GetExB(bz,kFALSE);
  if (!exb) return 0;
  Double_t pos0[3] = {r*TMath::Cos(phi), r*TMath::Sin(phi),z};
  Double_t pos1[3];
  exb->Correct(pos0,pos1);
  Double_t dx=pos1[0]-pos0[0];
  Double_t dy=pos1[1]-pos0[1];
  //  Double_t dz=pos1[2]-pos0[2];
  // return TMath::Sqrt(dx*dx+dy*dy);  
  Float_t dr = (dx*pos0[0]+dy*pos0[1])/r;
  return dr;
}


Double_t AliTPCExB::GetDrphi(Double_t r, Double_t phi, Double_t z, Double_t bz){
  //
  //
  //
  AliTPCExB *exb = Instance();
  if (!exb) exb = AliTPCcalibDB::GetExB(bz,kFALSE);
  if (!exb) return 0;
  Double_t pos0[3] = {r*TMath::Cos(phi), r*TMath::Sin(phi),z};
  Double_t pos1[3];
  exb->Correct(pos0,pos1);
  Double_t dphi=TMath::ATan2(pos1[1],pos1[0])-TMath::ATan2(pos0[1],pos0[0]);
  if (dphi>TMath::Pi()) dphi-=TMath::TwoPi();
  if (dphi<-TMath::Pi()) dphi+=TMath::TwoPi();
  return r*dphi;  

}


Double_t AliTPCExB::GetDphi(Double_t r, Double_t phi, Double_t z, Double_t bz){
  //
  //
  // 
  AliTPCExB *exb = Instance();
  if (!exb) exb = AliTPCcalibDB::GetExB(bz,kFALSE);
  if (!exb) return 0;
  Double_t pos0[3] = {r*TMath::Cos(phi), r*TMath::Sin(phi),z};
  Double_t pos1[3];
  exb->Correct(pos0,pos1);
  Double_t dphi=TMath::ATan2(pos1[1],pos1[0])-TMath::ATan2(pos0[1],pos0[0]);
  return dphi;  

}

Double_t AliTPCExB::GetDz(Double_t r, Double_t phi, Double_t z, Double_t bz){
  //
  //
  //
  AliTPCExB *exb = Instance();
  if (!exb) exb = AliTPCcalibDB::GetExB(bz,kFALSE);
  if (!exb) return 0;
  Double_t pos0[3] = {r*TMath::Cos(phi), r*TMath::Sin(phi),z};
  Double_t pos1[3];
  exb->Correct(pos0,pos1);
  Double_t dz=pos1[2]-pos0[2];
  return dz;  
}

//
// Magnetic field
//




void AliTPCExB::RegisterField(Int_t index, AliMagF * magf){
  //
  // add the filed to the list
  //
  fgArray.AddAt(magf,index);
}



Double_t AliTPCExB::GetBx(Double_t r, Double_t phi, Double_t z,Int_t index){
  //
  // 
  //
  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Double_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //  xyz[1]+=30;
  Double_t bxyz[3];
  mag->Field(xyz,bxyz);
  return bxyz[0];
}  

Double_t AliTPCExB::GetBy(Double_t r, Double_t phi, Double_t z,Int_t index){
  //
  // 
  //
  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Double_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //  xyz[1]+=30;
  Double_t bxyz[3];
  mag->Field(xyz,bxyz);
  return bxyz[1];
}  

Double_t AliTPCExB::GetBz(Double_t r, Double_t phi, Double_t z,Int_t index){
  //
  // 
  //
  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Double_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //  xyz[1]+=30;
  Double_t bxyz[3];
  mag->Field(xyz,bxyz);
  return bxyz[2];
}  



Double_t AliTPCExB::GetBr(Double_t r, Double_t phi, Double_t z,Int_t index){
  //
  // 
  //
  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Double_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //xyz[1]+=30;
  Double_t bxyz[3];
  mag->Field(xyz,bxyz);
  if (r==0) return 0;
  Double_t br = (bxyz[0]*xyz[0]+bxyz[1]*xyz[1])/r;
  return br;
}  

Double_t AliTPCExB::GetBrfi(Double_t r, Double_t phi, Double_t z,Int_t index){
  //
  // 
  //
  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Double_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //xyz[1]+=30;
  Double_t bxyz[3];
  mag->Field(xyz,bxyz);
  if (r==0) return 0;
  Double_t br = (-bxyz[0]*xyz[1]+bxyz[1]*xyz[0])/r;
  return br;
}  




Double_t AliTPCExB::GetBxI(Double_t r, Double_t phi, Double_t z,Int_t index)
{
  Double_t sumf  =0;
  if (z>0 &&z<250){
    for (Float_t zi=z;zi<250;zi+=5){
      sumf+=GetBx(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  if (z<0 &&z>-250){
    for (Float_t zi=z;zi>-250;zi-=5){
      sumf+=GetBx(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  return sumf*5;
}

Double_t AliTPCExB::GetByI(Double_t r, Double_t phi, Double_t z,Int_t index)
{
  Double_t sumf =0;
  if (z>0 &&z<250){
    for (Float_t zi=z;zi<250;zi+=5){
      sumf+=GetBy(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  if (z<0 &&z>-250){
    for (Float_t zi=z;zi>-250;zi-=5){
      sumf+=GetBy(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  return sumf*5;
}

Double_t AliTPCExB::GetBzI(Double_t r, Double_t phi, Double_t z,Int_t index)
{
  Double_t sumf =0;
  if (z>0 &&z<250){
    for (Float_t zi=z;zi<250;zi+=5){
      sumf+=GetBz(r,phi,zi,index);
    }
  }
  if (z<0 &&z>-250){
    for (Float_t zi=z;zi>-250;zi-=5){
      sumf+=GetBz(r,phi,zi,index);
    }
  }
  return sumf*5;
}


Double_t AliTPCExB::GetBrI(Double_t r, Double_t phi, Double_t z,Int_t index)
{
  Double_t sumf =0;
  if (z>0 &&z<250){
    for (Float_t zi=z;zi<250;zi+=5){
      sumf+=GetBr(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  if (z<0 &&z>-250){
    for (Float_t zi=z;zi>-250;zi-=5){
      sumf+=GetBr(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  return sumf*5;
}

Double_t AliTPCExB::GetBrfiI(Double_t r, Double_t phi, Double_t z,Int_t index)
{
  Double_t sumf =0;
  if (z>0 &&z<250){
    for (Float_t zi=z;zi<250;zi+=5.){
      sumf+=GetBrfi(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  if (z<0 &&z>-250){
    for (Float_t zi=z;zi>-250;zi-=5){
      sumf+=GetBrfi(r,phi,zi,index)/GetBz(r,phi,zi,index);
    }
  }
  return sumf*5;
}


Double_t AliTPCExB::Eval(Int_t type, Double_t r, Double_t phi, Double_t z){
  //
  // Evaluate parameterization
  //
  // br integral param 
  if (type==0) {
    if (z>0 && fMatBrBzI0) return EvalMat(*fMatBrBzI0,r,phi,z);
    if (z<0 && fMatBrBzI1) return EvalMat(*fMatBrBzI1,r,phi,z);
  }
  // brfi integral param   
  if (type==1) {
    if (z>0 && fMatBrfiBzI0) return EvalMat(*fMatBrfiBzI0,r,phi,z);
    if (z<0 && fMatBrfiBzI1) return EvalMat(*fMatBrfiBzI1,r,phi,z);
  }
  // brbz param
  if (type==2 && fMatBrBz) return EvalMat(*fMatBrBz,r,phi,z);
  // brfibz param
  if (type==3 && fMatBrfiBz) return EvalMat(*fMatBrfiBz,r,phi,z);
  return 0;
}


Double_t AliTPCExB::EvalMat(const TVectorD &vec, Double_t r, Double_t phi, Double_t z){
  //
  // Evaluate taylor expansion in r,phi,z
  //
  // Variables  
  //tree->SetAlias("sa","sin(phi+0.0)");
  //tree->SetAlias("ca","cos(phi+0.0)");
  //tree->SetAlias("sa2","sin(phi*2+0.0)");
  //tree->SetAlias("ca2","cos(phi*2+0.0)");
  //tree->SetAlias("zn","(x2/250.)");
  //tree->SetAlias("rn","(r/250.)")
  //   TString fstringSym="";
  //   //  
  //   fstringSym+="zn++";
  //   fstringSym+="rn++";
  //   fstringSym+="zn*rn++";
  //   fstringSym+="zn*zn++";
  //   fstringSym+="zn*zn*rn++";
  //   fstringSym+="zn*rn*rn++";
  //   //
  //   fstringSym+="sa++";
  //   fstringSym+="ca++";  
  //   fstringSym+="ca2++";
  //   fstringSym+="sa2++";
  //   fstringSym+="ca*zn++";
  //   fstringSym+="sa*zn++";
  //   fstringSym+="ca2*zn++";
  //   fstringSym+="sa2*zn++";
  //   fstringSym+="ca*zn*zn++";
  //   fstringSym+="sa*zn*zn++";
  //   fstringSym+="ca*zn*rn++";
  //   fstringSym+="sa*zn*rn++";


  Double_t sa  = TMath::Sin(phi);
  Double_t ca  = TMath::Cos(phi);
  Double_t sa2 = TMath::Sin(phi*2);
  Double_t ca2 = TMath::Cos(phi*2);
  Double_t zn  = z/250.;
  Double_t rn  = r/250.;
  Int_t  ipoint=0;
  Double_t res = vec[ipoint++];
  res+=vec[ipoint++]*zn;
  res+=vec[ipoint++]*rn;
  res+=vec[ipoint++]*zn*rn;
  res+=vec[ipoint++]*zn*zn;
  res+=vec[ipoint++]*zn*zn*rn;
  res+=vec[ipoint++]*zn*rn*rn;
  //
  res+=vec[ipoint++]*sa;
  res+=vec[ipoint++]*ca;  
  res+=vec[ipoint++]*ca2;
  res+=vec[ipoint++]*sa2;
  res+=vec[ipoint++]*ca*zn;
  res+=vec[ipoint++]*sa*zn;
  res+=vec[ipoint++]*ca2*zn;
  res+=vec[ipoint++]*sa2*zn;
  res+=vec[ipoint++]*ca*zn*zn;
  res+=vec[ipoint++]*sa*zn*zn;
  res+=vec[ipoint++]*ca*zn*rn;
  res+=vec[ipoint++]*sa*zn*rn;
  return res;
}














/*
  
 AliTPCExB draw;
 draw.RegisterField(0,new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG));
 draw.RegisterField(1,new AliMagFMaps("Maps","Maps", 2, 1., 10., 2));

 TF2 fbz_rz_0pi("fbz_rz_0pi","AliTPCExB::GetBz(x,0*pi,y)",0,250,-250,250);
 fbz_rz_0pi->Draw("surf2");
 
 TF1 fbz_z_90_00pi("fbz_z_90_00pi","AliTPCExB::GetBz(90,0*pi,x)",-250,250);
  TF1 fbz_z_90_05pi("fbz_z_90_05pi","AliTPCExB::GetBz(90,0.5*pi,x)",-250,250);
  TF1 fbz_z_90_10pi("fbz_z_90_10pi","AliTPCExB::GetBz(90,1.0*pi,x)",-250,250);
  TF1 fbz_z_90_15pi("fbz_z_90_15pi","AliTPCExB::GetBz(90,1.5*pi,x)",-250,250);
  fbz_z_90_00pi->SetLineColor(2);
  fbz_z_90_05pi->SetLineColor(3);
  fbz_z_90_10pi->SetLineColor(4);
  fbz_z_90_15pi->SetLineColor(5);
  fbz_z_90_00pi->Draw()
  fbz_z_90_05pi->Draw("same")
  fbz_z_90_15pi->Draw("same")
  fbz_z_90_10pi->Draw("same")
  

  TF1 fbr_z_90_00pi("fbz_z_90_00pi","AliTPCExB::GetBr(90,0*pi,x)",-250,250);
  TF1 fbr_z_90_05pi("fbz_z_90_05pi","AliTPCExB::GetBr(90,0.5*pi,x)",-250,250);
  TF1 fbr_z_90_10pi("fbz_z_90_10pi","AliTPCExB::GetBr(90,1.0*pi,x)",-250,250);
  TF1 fbr_z_90_15pi("fbz_z_90_15pi","AliTPCExB::GetBr(90,1.5*pi,x)",-250,250);
  fbr_z_90_00pi->SetLineColor(2);
  fbr_z_90_05pi->SetLineColor(3);
  fbr_z_90_10pi->SetLineColor(4);
  fbr_z_90_15pi->SetLineColor(5);
  fbr_z_90_00pi->Draw()
  fbr_z_90_05pi->Draw("same")
  fbr_z_90_15pi->Draw("same")
  fbr_z_90_10pi->Draw("same")

  //
  TF2 fbz_xy_0z("fbz_xy_0z","AliTPCExB::GetBz(sqrt(x^2+y^2),atan2(y,x),0)",-250,250,-250,250);
  fbz_xy_0z.SetNpy(100);
  fbz_xy_0z.SetNpx(100);
  fbz_xy_0z->Draw("colz");
  //
  TF2 fbz_xy_250z("fbz_xy_250z","AliTPCExB::GetBz(sqrt(x^2+y^2),atan2(y,x),250)",-250,250,-250,250);
  fbz_xy_250z.SetNpy(100);
  fbz_xy_250z.SetNpx(100)
  fbz_xy_250z->Draw("colz");
  //
   TF2 fbz_xy_m250z("fbz_xy_m250z","AliTPCExB::GetBz(sqrt(x^2+y^2),atan2(y,x),-250)",-250,250,-250,250);
  fbz_xy_m250z.SetNpy(100);
  fbz_xy_m250z.SetNpx(100)
  fbz_xy_m250z->Draw("colz");
  //





*/


