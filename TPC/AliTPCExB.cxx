#include "AliTPCExB.h"
#include "TMath.h"
#include "TTreeStream.h"


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
  Draw0(0)




*/

AliTPCExB* AliTPCExB::fgInstance = 0;

ClassImp(AliTPCExB)


void AliTPCExB::TestExB(const char* fileName) {
  //
  // well, as the name sais...
  //
  TTreeSRedirector ts(fileName);
  Double_t x[3];
  for (x[0]=-250.;x[0]<=250.;x[0]+=10.)
    for (x[1]=-250.;x[1]<=250.;x[1]+=10.)
      for (x[2]=-250.;x[2]<=250.;x[2]+=10.) {
	Double_t d[3];
	Correct(x,d);
	Double_t r=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
	Double_t rd=TMath::Sqrt(d[0]*d[0]+d[1]*d[1]);
	Double_t dr=r-rd;
	Double_t phi=TMath::ATan2(x[0],x[1]);
	Double_t phid=TMath::ATan2(d[0],d[1]);
	Double_t dphi=phi-phid;
	if (dphi<0.) dphi+=TMath::TwoPi();
	if (dphi>TMath::Pi()) dphi=TMath::TwoPi()-dphi;
	Double_t drphi=r*dphi;
	Double_t dx=x[0]-d[0];
	Double_t dy=x[1]-d[1];
	Double_t dz=x[2]-d[2];
	ts<<"positions"
	  <<"x0="<<x[0]
	  <<"x1="<<x[1]
	  <<"x2="<<x[2]
	  <<"dx="<<dx
	  <<"dy="<<dy
	  <<"dz="<<dz
	  <<"r="<<r
	  <<"phi="<<phi
	  <<"dr="<<dr
	  <<"drphi="<<drphi
	  <<"\n";
      }
}



Double_t AliTPCExB::GetDr(Double_t r, Double_t phi, Double_t z){
  //
  // Static function
  // Posibble to us it for visualization 
  // 
  //
  if (!fgInstance) return 0;
  Double_t pos0[3] = {r*TMath::Sin(phi), r*TMath::Cos(phi),z};
  Double_t pos1[3];
  fgInstance->Correct(pos0,pos1);
  Double_t dx=pos1[0]-pos0[0];
  Double_t dy=pos1[1]-pos0[1];
  //  Double_t dz=pos1[2]-pos0[2];
  return TMath::Sqrt(dx*dx+dy*dy);  
}


Double_t AliTPCExB::GetDrphi(Double_t r, Double_t phi, Double_t z){
  //
  //
  //
  if (!fgInstance) return 0;
  Double_t pos0[3] = {r*TMath::Sin(phi), r*TMath::Cos(phi),z};
  Double_t pos1[3];
  fgInstance->Correct(pos0,pos1);
  Double_t dphi=TMath::ATan2(pos1[1],pos1[0])-TMath::ATan2(pos0[1],pos0[0]);
  return r*dphi;  

}


Double_t AliTPCExB::GetDphi(Double_t r, Double_t phi, Double_t z){
  //
  //
  //
  if (!fgInstance) return 0;
  Double_t pos0[3] = {r*TMath::Sin(phi), r*TMath::Cos(phi),z};
  Double_t pos1[3];
  fgInstance->Correct(pos0,pos1);
  Double_t dphi=TMath::ATan2(pos1[1],pos1[0])-TMath::ATan2(pos0[1],pos0[0]);
  return dphi;  

}

Double_t AliTPCExB::GetDz(Double_t r, Double_t phi, Double_t z){
  //
  //
  //
  if (!fgInstance) return 0;
  Double_t pos0[3] = {r*TMath::Sin(phi), r*TMath::Cos(phi),z};
  Double_t pos1[3];
  fgInstance->Correct(pos0,pos1);
  Double_t dz=pos1[2]-pos0[2];
  return dz;  
}
