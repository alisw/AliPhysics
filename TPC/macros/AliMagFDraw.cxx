/// \class AliMagFDraw
///
/// ~~~{.cpp}
/// .L $ALICE_ROOT/TPC/macros/AliMagFDraw.cxx+
/// AliMagFDraw draw;
/// draw.RegisterField(0,new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG));
/// draw.RegisterField(1,new AliMagFMaps("Maps","Maps", 2, 1., 10., 2));
///
/// TF2 fbz_rz_0pi("fbz_rz_0pi","AliMagFDraw::GetBz(x,0*pi,y)",0,250,-250,250);
/// fbz_rz_0pi->Draw("surf2");
/// 
/// TF1 fbz_z_90_00pi("fbz_z_90_00pi","AliMagFDraw::GetBz(90,0*pi,x)",-250,250);
/// TF1 fbz_z_90_05pi("fbz_z_90_05pi","AliMagFDraw::GetBz(90,0.5*pi,x)",-250,250);
/// TF1 fbz_z_90_10pi("fbz_z_90_10pi","AliMagFDraw::GetBz(90,1.0*pi,x)",-250,250);
/// TF1 fbz_z_90_15pi("fbz_z_90_15pi","AliMagFDraw::GetBz(90,1.5*pi,x)",-250,250);
/// fbz_z_90_00pi->SetLineColor(2);
/// fbz_z_90_05pi->SetLineColor(3);
/// fbz_z_90_10pi->SetLineColor(4);
/// fbz_z_90_15pi->SetLineColor(5);
/// fbz_z_90_00pi->Draw()
/// fbz_z_90_05pi->Draw("same")
/// fbz_z_90_15pi->Draw("same")
/// fbz_z_90_10pi->Draw("same")
/// 
/// TF1 fbr_z_90_00pi("fbz_z_90_00pi","AliMagFDraw::GetBr(90,0*pi,x)",-250,250);
/// TF1 fbr_z_90_05pi("fbz_z_90_05pi","AliMagFDraw::GetBr(90,0.5*pi,x)",-250,250);
/// TF1 fbr_z_90_10pi("fbz_z_90_10pi","AliMagFDraw::GetBr(90,1.0*pi,x)",-250,250);
/// TF1 fbr_z_90_15pi("fbz_z_90_15pi","AliMagFDraw::GetBr(90,1.5*pi,x)",-250,250);
/// fbr_z_90_00pi->SetLineColor(2);
/// fbr_z_90_05pi->SetLineColor(3);
/// fbr_z_90_10pi->SetLineColor(4);
/// fbr_z_90_15pi->SetLineColor(5);
/// fbr_z_90_00pi->Draw()
/// fbr_z_90_05pi->Draw("same")
/// fbr_z_90_15pi->Draw("same")
/// fbr_z_90_10pi->Draw("same")
///
/// TF2 fbz_xy_0z("fbz_xy_0z","AliMagFDraw::GetBz(sqrt(x^2+y^2),atan2(y,x),0)",-250,250,-250,250);
/// fbz_xy_0z.SetNpy(100);
/// fbz_xy_0z.SetNpx(100);
/// fbz_xy_0z->Draw("colz");
///
/// TF2 fbz_xy_250z("fbz_xy_250z","AliMagFDraw::GetBz(sqrt(x^2+y^2),atan2(y,x),250)",-250,250,-250,250);
/// fbz_xy_250z.SetNpy(100);
/// fbz_xy_250z.SetNpx(100)
/// fbz_xy_250z->Draw("colz");
///
///  TF2 fbz_xy_m250z("fbz_xy_m250z","AliMagFDraw::GetBz(sqrt(x^2+y^2),atan2(y,x),-250)",-250,250,-250,250);
/// fbz_xy_m250z.SetNpy(100);
/// fbz_xy_m250z.SetNpx(100)
/// fbz_xy_m250z->Draw("colz");
/// ~~~

#include "TObjArray.h"
#include "TMath.h"
#include "AliMagF.h"
#include "TLinearFitter.h"
#include "TString.h"

class AliMagFDraw  : public AliMagF
{
public:
  AliMagFDraw():AliMagF(){}
  static  void RegisterField(Int_t index, AliMagF * magf);
  static  Double_t GetBx(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBy(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBz(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBr(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  Double_t GetBrfi(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  //static  Double_t GetBr2(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  //static  Double_t GetBrfi2(Double_t r, Double_t phi, Double_t z,Int_t index=0);
  static  TObjArray *Fit(const char *formula, Int_t index=0);
public:
  static TObjArray   fgArray;
  ClassDef(AliMagFDraw,2) 
};


/// \cond CLASSIMP
ClassImp(AliMagFDraw)
/// \endcond


TObjArray   AliMagFDraw::fgArray;

void AliMagFDraw::RegisterField(Int_t index, AliMagF * magf){
  /// add the filed to the list

  fgArray.AddAt(magf,index);
}

Double_t AliMagFDraw::GetBz(Double_t r, Double_t phi, Double_t z,Int_t index){
  ///

  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Float_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //  xyz[1]+=30;
  Float_t bxyz[3];
  mag->Field(xyz,bxyz);
  return bxyz[2];
}  

Double_t AliMagFDraw::GetBy(Double_t r, Double_t phi, Double_t z,Int_t index){
  ///

  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Float_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //  xyz[1]+=30;
  Float_t bxyz[3];
  mag->Field(xyz,bxyz);
  return bxyz[1];
}  


Double_t AliMagFDraw::GetBx(Double_t r, Double_t phi, Double_t z,Int_t index){
  ///

  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Float_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //  xyz[1]+=30;
  Float_t bxyz[3];
  mag->Field(xyz,bxyz);
  return bxyz[0];
}  




Double_t AliMagFDraw::GetBr(Double_t r, Double_t phi, Double_t z,Int_t index){
  ///

  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Float_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //xyz[1]+=30;
  Float_t bxyz[3];
  mag->Field(xyz,bxyz);
  if (r==0) return 0;
  Float_t br = bxyz[0]*xyz[0]/r+bxyz[1]*xyz[1]/r;
  return br;
}  

Double_t AliMagFDraw::GetBrfi(Double_t r, Double_t phi, Double_t z,Int_t index){
  ///

  AliMagF *mag = (AliMagF*)fgArray.At(index);
  if (!mag) return 0;
  Float_t xyz[3]={r*TMath::Cos(phi),r*TMath::Sin(phi),z};
  //xyz[1]+=30;
  Float_t bxyz[3];
  mag->Field(xyz,bxyz);
  if (r==0) return 0;
  Float_t br = -bxyz[0]*xyz[1]/r+bxyz[1]*xyz[0]/r;
  return br;
}  


TObjArray * AliMagFDraw::Fit(const char *formula, Int_t index){
  /// formula=`1++x+x^2++cos(y)++cos(y)^2++z++z^2`
  /// index=0
  TObjArray *fstrings = TString(formula).Tokenize("++");
  Int_t ndim = fstrings->GetEntries();
  TObjArray *formulas = new TObjArray(ndim);
  for (Int_t i=0;i<ndim;i++){
    formulas->AddAt(new TFormula(Form("fff_%d",i),fstrings->At(i)->GetName()),i);		    
  }  
  TLinearFitter * fitR   = new TLinearFitter(ndim+1,Form("hyp%d",ndim));
  TLinearFitter * fitRFI = new TLinearFitter(ndim+1,Form("hyp%d",ndim));
  TLinearFitter * fitZ   = new TLinearFitter(ndim+1,Form("hyp%d",ndim));
  Double_t x[ndim];    
  for (Float_t r=20; r<250;r+=20){
    for (Float_t fi=0; fi<TMath::Pi()*2;fi+=0.2){
      for (Float_t z=-250; z<250;z+=20){
	for (Int_t ifor=0;ifor<ndim;ifor++){
	  x[ifor]= ((TFormula*)formulas->At(ifor))->Eval(r/250.,fi,z/250.);
	}
	fitR->AddPoint(x,AliMagFDraw::GetBr(r,fi,z,index));
	fitRFI->AddPoint(x,AliMagFDraw::GetBrfi(r,fi,z,index));
	fitZ->AddPoint(x,AliMagFDraw::GetBz(r,fi,z,index));     
      }
    }
  }
  fitR->Eval();  
  fitRFI->Eval();  
  fitZ->Eval();  
  TObjArray *res = new TObjArray; 
  res->AddAt(fitR,0);
  res->AddAt(fitRFI,1);
  res->AddAt(fitZ,2);
  printf("\tchi2\tn\tRMS\n");
  printf("\t%f\t%d\t%f\n",fitR->GetChisquare(),fitR->GetNpoints(),TMath::Sqrt(fitR->GetChisquare()/fitR->GetNpoints()));
  printf("\t%f\t%d\t%f\n",fitRFI->GetChisquare(),fitRFI->GetNpoints(),TMath::Sqrt(fitRFI->GetChisquare()/fitRFI->GetNpoints()));
  printf("\t%f\t%d\t%f\n",fitZ->GetChisquare(),fitZ->GetNpoints(),TMath::Sqrt(fitZ->GetChisquare()/fitZ->GetNpoints()));

  TFormula * funBZ = new TFormula("funBZ",formula);
  TFormula * funBR = new TFormula("funBR",formula);
  TFormula * funBRFI = new TFormula("funBRFI",formula);
  TVectorD vec;
  fitR->GetParameters(vec);
  funBR->SetParameters(vec.GetMatrixArray());
  fitRFI->GetParameters(vec);
  funBRFI->SetParameters(vec.GetMatrixArray());
  fitZ->GetParameters(vec);
  funBZ->SetParameters(vec.GetMatrixArray());

  return res;
}


TString MakeString(){
  TString str="";
  {
    Int_t counter=0;
    for (Int_t ix=0;ix<3;ix++)
      for (Int_t iz=0;iz<3;iz++){
	if (ix+iz>0) {
	  str+=Form("x^%d*z^%d++",ix,iz);
	  printf("x^%d*z^%d++\n",ix,iz);
	}
	for (Int_t iy=1;iy<4;iy++){
	  if (ix+iz+iy==0) continue;
	  str+=Form("x^%d*z^%d*sin(y*%d)++",ix,iz,iy);
	  str+=Form("x^%d*z^%d*cos(y*%d)++",ix,iz,iy);
	  printf(   "x^%d*z^%d*sin(y*%d)++\n",ix,iz,iy);
	  printf(   "x^%d*z^%d*cos(y*%d)++\n",ix,iz,iy);
	  counter++;
	}
      }
    str[str.Length()-2]=0;
  }
  return str;
}



/*

TObjArray * array = AliMagFDraw::Fit("x",0)
//
chi2    n       RMS
15.534116       9600    0.040226
4.422160        9600    0.021463
8.378622        9600    0.029543

TObjArray * array = AliMagFDraw::Fit("x++z++z^2++x*sin(y)++x*cos(y)++z*x*sin(y)++z*x*cos(y)++z^2*x*sin(y)++z^2*x*cos(y)++z^2*x*sin(y)^2++z^2*x*cos(y)^2++z^3*x*sin(y)^2++z^3*x*cos(y)^2",0);
chi2    n       RMS
1.842443        9600    0.013854
0.957056        9600    0.009985
0.097799        9600    0.003192

TObjArray * array = AliMagFDraw::Fit("x++z++z^2++x*sin(y)++x*cos(y)++z*x*sin(y)++z*x*cos(y)++z^2*x*sin(y)++z^2*x*cos(y)++z^2*x*sin(y)^2++z^2*x*cos(y)^2++z^3*x*sin(y)^2++z^3*x*cos(y)++z^3*sin(y)++z^3*cos(y)++z^3*x*cos(y)++z^2*sin(y)++z^2*cos(y)++z*sin(y)++z*cos(y)++sin(y)++cos(y)",0);

chi2    n       RMS
1.564687        9600    0.012767
0.002291        9600    0.000489
0.097063        9600    0.003180


TObjArray * array = AliMagFDraw::Fit(MakeString(),0);

chi2    n       RMS
0.000303        9600    0.000178
0.000066        9600    0.000083
0.003110        9600    0.000569


}
*/

