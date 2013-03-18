/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliDielectronBtoJPSItoEleCDFfitFCNfitter.h"
#include "TArrayD.h"
#include "TFormula.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

AliDielectronBtoJPSItoEleCDFfitFCNfitter::AliDielectronBtoJPSItoEleCDFfitFCNfitter():
 fInvMass(0x0),
 fX(0x0),
 fX2D(0x0),
 fFCN(0x0),
 fFitOpt("")
{
 //
 // ctor
 //
 fFCN = new AliDielectronBtoJPSItoEleCDFfitFCN();
 fFCN->SetCrystalBallFunction(kTRUE);
 for(Int_t i=0; i<kInvMassTotal; i++) fParameterInvMassToFix[i]=kFALSE;
 for(Int_t i=0; i<kPseudo; i++)  fParameterXToFix[i]=kFALSE;
 for(Int_t i=0; i<kPseudoBkg; i++)  fParameterXbkgToFix[i]=kFALSE;
}
//__________________________________________________________________________________________
AliDielectronBtoJPSItoEleCDFfitFCNfitter::~AliDielectronBtoJPSItoEleCDFfitFCNfitter()
{
 //
 // dtor
 //
 if(fInvMass) delete fInvMass;
 if(fX) delete fX;
 if(fX2D) delete fX2D;
 if(fFCN) delete fFCN;
}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::SetInvMassSignalParameters(const Double_t massPar[5])
{
 //
 // Setter of the inv mass signal distribution parameters
 //
fFCN->SetCrystalBallMmean(massPar[0]); //9
fFCN->SetCrystalBallNexp(massPar[1]); //10
fFCN->SetCrystalBallSigma(massPar[2]); //11
fFCN->SetCrystalBallAlpha(massPar[3]); //12
fFCN->SetCrystalBallNorm(massPar[4]); //13
 //for(Int_t i=0; i<kInvMassSignal; i++) printf(" par%i %f \n",i,massPar[i]);
}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::SetInvMassParameters(const Double_t massPar[9])
{
 //
 // Setter of the inv mass total distribution parameters
 //
 Double_t massParSig[5];
 Double_t massParBkg[4];

 for(Int_t iPar=0; iPar<9; iPar++){
  if(iPar<5) massParSig[iPar]=massPar[iPar];
  else massParBkg[iPar-5]=massPar[iPar];
 }
 SetInvMassSignalParameters(massParSig);
 SetInvMassBkgParameters(massParBkg);
}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::SetInvMassBkgParameters(const Double_t massPar[4])
{
 //
 // Setter of the Inv Mass background distribution parameters
 //
fFCN->SetBkgInvMassNorm(massPar[0]);//14
fFCN->SetBkgInvMassMean(massPar[1]);//15
fFCN->SetBkgInvMassSlope(massPar[2]);//16
fFCN->SetBkgInvMassConst(massPar[3]);//17
}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::GetInvMassParameters(TArrayD &massPar)
{
 //
 // Getter of the inv mass total distribution parameters
 //
 massPar.Reset();
 massPar.Set(kInvMassTotal);
 massPar.AddAt(fFCN->GetCrystalBallMmean(),0); //9
 massPar.AddAt(fFCN->GetCrystalBallNexp(),1); //10
 massPar.AddAt(fFCN->GetCrystalBallSigma(),2); //11
 massPar.AddAt(fFCN->GetCrystalBallAlpha(),3); //12
 massPar.AddAt(fFCN->GetCrystalBallNorm(),4); //13
 massPar.AddAt(fFCN->GetBkgInvMassNorm(),5);//14
 massPar.AddAt(fFCN->GetBkgInvMassMean(),6);//15
 massPar.AddAt(fFCN->GetBkgInvMassSlope(),7); //16
 massPar.AddAt(fFCN->GetBkgInvMassConst(),8); //17
}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::GetInvMassSignalParameters(TArrayD &massPar)
{
 //
 // Getter of the Inv Mass signal distribution parameters
 //
 massPar.Reset();
 massPar.Set(kInvMassSignal);
 massPar.AddAt(fFCN->GetCrystalBallMmean(),0); //9
 massPar.AddAt(fFCN->GetCrystalBallNexp(),1); //10
 massPar.AddAt(fFCN->GetCrystalBallSigma(),2); //11
 massPar.AddAt(fFCN->GetCrystalBallAlpha(),3); //12
 massPar.AddAt(fFCN->GetCrystalBallNorm(),4); //13

}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::GetInvMassBkgParameters(TArrayD &massPar)
{
 //
 // Getter of the Inv Mass background distribution parameters
 //
 massPar.Reset();
 massPar.Set(kInvMassBkg);
 massPar.AddAt(fFCN->GetBkgInvMassNorm(),0); //14
 massPar.AddAt(fFCN->GetBkgInvMassMean(),1); //15
 massPar.AddAt(fFCN->GetCrystalBallSigma(),2); //11
 massPar.AddAt(fFCN->GetCrystalBallAlpha(),3); //12
}

//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::SetPseudoProperBkgParameters(Double_t xBkgPar[7])
{
 //
 // Setter of the psudoproper background distribution parameters
 //
 fFCN->SetResWeight(xBkgPar[0]);//0
 fFCN->SetFPlus(xBkgPar[1]);//1
 fFCN->SetFMinus(xBkgPar[2]);//2
 fFCN->SetFSym(xBkgPar[3]);//3
 fFCN->SetLamPlus(xBkgPar[4]); //4
 fFCN->SetLamMinus(xBkgPar[5]); // 5	
 fFCN->SetLamSym(xBkgPar[6]); // 6

}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::GetPseudoProperBkgParameters(TArrayD &xBkgPar)
{
 //
 // Getter of the psudoproper background distribution parameters
 //
 xBkgPar.Reset();
 xBkgPar.Set(kPseudoBkg);
 xBkgPar.AddAt(fFCN->GetResWeight(),0);
 xBkgPar.AddAt(fFCN->GetFPlus(),1);
 xBkgPar.AddAt(fFCN->GetFMinus(),2);
 xBkgPar.AddAt(fFCN->GetFSym(),3);
 xBkgPar.AddAt(fFCN->GetLamPlus(),4);
 xBkgPar.AddAt(fFCN->GetLamMinus(),5);
 xBkgPar.AddAt(fFCN->GetLamSym(),6);
}
//__________________________________________________________________________________________
TH1F* AliDielectronBtoJPSItoEleCDFfitFCNfitter::FitInvMassSignal(Double_t norm,Double_t mMin, Double_t mMax)
{
 //
 // Fit on the invariant mass signal distribution ( Crystal ball params ) 
 //
 if(!fInvMass){
  printf("no histogram...exiting\n");
  return 0x0;
 }

 TArrayD pars;
 GetInvMassSignalParameters(pars);

 TF1 * invMass = new TF1("invMassSignal",this,&AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFInvMassSignal,mMin,mMax,kInvMassSignal+1);
 for(Int_t i=0; i<kInvMassSignal ; i++) {
  if(fParameterInvMassToFix[i]) invMass->FixParameter(i,pars.At(i));
  else invMass->SetParameter(i,pars.At(i));
 }
 invMass->SetParameter(kInvMassSignal,norm);
 fInvMass->Fit("invMassSignal",Form("R%s",fFitOpt.Data()),"",mMin,mMax);
 return fInvMass;
}
//__________________________________________________________________________________________
TH1F * AliDielectronBtoJPSItoEleCDFfitFCNfitter::FitInvMass(Double_t norm[2],Double_t mMin, Double_t mMax)
{
 //
 // Fit on the invariant mass total distribution
 //
 if(!fInvMass){
  printf("no histogram...exiting\n");
  return 0x0;
 }

 TArrayD pars;
 GetInvMassParameters(pars);

 TF1 * invMassTot = new TF1("invMassTot",this,&AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFInvMassTotal,mMin,mMax,11);

 for(Int_t i=0; i<kInvMassTotal ; i++) {
  if(fParameterInvMassToFix[i]) invMassTot->FixParameter(i,pars.At(i));
  else invMassTot->SetParameter(i,pars.At(i));
 }

 if(!fParameterInvMassToFix[5]) invMassTot->SetParLimits(5,0,999999999);


 invMassTot->SetParameter(9,norm[0]);
 invMassTot->SetParLimits(9,0,999999999);
 invMassTot->SetParameter(10,norm[1]);
 invMassTot->SetParLimits(10,0,999999999);

 fInvMass->Fit("invMassTot",Form("R%s",fFitOpt.Data()),"",mMin,mMax);

 //printf("Inv Mass Params : ");
 //for(Int_t i=0; i< 11; i++) printf("p%i %f ",i,invMassTot->GetParameter(i));
 //printf("\n");

 return fInvMass;
}
//__________________________________________________________________________________________
void AliDielectronBtoJPSItoEleCDFfitFCNfitter::GetPseudoProperParameters(TArrayD &xPar, Int_t type)
{
 //
 // Getter of the psudoproper background distribution parameters per dielectron type (FF, FS, SS)
 //
 xPar.Reset();
 xPar.Set(kPseudo);
 xPar.AddAt(fFCN->GetNormGaus1ResFunc(type),0); 
 xPar.AddAt(fFCN->GetNormGaus2ResFunc(type),3); 
 xPar.AddAt(fFCN->GetResMean1(type),1); 
 xPar.AddAt(fFCN->GetResSigma1(type),2); 
 xPar.AddAt(fFCN->GetResMean2(type),4); 
 xPar.AddAt(fFCN->GetResSigma2(type),5); 
 xPar.AddAt(fFCN->GetResAlfa(type),6); 
 xPar.AddAt(fFCN->GetResLambda(type),7); 
 xPar.AddAt(fFCN->GetResNormExp(type),8); 
}
//__________________________________________________________________________________________
TH1F * AliDielectronBtoJPSItoEleCDFfitFCNfitter::FitResolutionFunction(Double_t xMin, Double_t xMax, Int_t type, Double_t norm)
{
 //
 // Fit on the resolution function distribution per dielectron type
 //
 if(!fX){
  printf("no Pseudoproper decayhistogram...exiting\n");
  return 0x0;
 }

 TArrayD pars;
 GetPseudoProperParameters(pars,type);

 TF1 * xPrompt = new TF1("xPrompt",this,&AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFResolutionFunction,xMin,xMax,11);
 for(Int_t i=0; i<kPseudo ; i++) {
  if(fParameterXToFix[i]) xPrompt->FixParameter(i,pars.At(i));
  else xPrompt->SetParameter(i,pars.At(i));
  printf("%f ",pars.At(i));
 }
 printf("\n");
 if(!fParameterXToFix[3]) xPrompt->SetParLimits(3,0,999999999);
 xPrompt->FixParameter(9,type);
 xPrompt->SetParameter(10,norm);
 fX->Fit("xPrompt",Form("R%s",fFitOpt.Data()),"",xMin,xMax);

// printf("Res Func type %i : ",type);
// for(Int_t i=0; i< 11; i++) printf("p%i %6.3f ",i,xPrompt->GetParameter(i));
// printf("\n");
 return fX;
}

//__________________________________________________________________________________________
TH2F * AliDielectronBtoJPSItoEleCDFfitFCNfitter::FitBkgPsudoProper(Double_t xMin, Double_t xMax, Double_t norm)
{
 //
 // Fit on the resolution function distribution of the background ( in 2 D, whose one D is the dielectron
 // type)
 //

 if(!fX2D){
  printf("no BKG Pseudoproper decay 2D histogram ...exiting\n");
  return 0x0;
 }

 TArrayD bkgParam;
 GetPseudoProperBkgParameters(bkgParam);

 TF2 *psProperBkgnd = new TF2("bkgPseudoProper",this,
   &AliDielectronBtoJPSItoEleCDFfitFCNfitter::PsProperBackFitFunc,xMin,xMax,-0.5,2.5,8,"AliDielectronBtoJPSItoEleCDFfitFCNfitter","PsProperBackFitFunc");

 for(Int_t i=0; i< kPseudoBkg ; i++) {
  if(fParameterXbkgToFix[i]) psProperBkgnd->FixParameter(i,bkgParam.At(i));
  else psProperBkgnd->SetParameter(i,bkgParam.At(i));
 }
 psProperBkgnd->SetParameter(kPseudoBkg,norm);

 fX2D->Fit("bkgPseudoProper",Form("R%s",fFitOpt.Data()),"LEGO2",xMin,xMax);

 return fX2D;
}

//__________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFInvMassSignal(const Double_t *x,const Double_t *par)
{
 //
 // function
 //
 SetInvMassSignalParameters(par);
 return par[5]*fFCN->EvaluateCDFInvMassSigDistr(x[0]);
}
//__________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFInvMassBkg(const Double_t *x,const Double_t *par)
{
 //
 // function
 //
 SetInvMassBkgParameters(par);
 return fFCN->EvaluateCDFInvMassBkgDistr(x[0]);
}
//__________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFInvMassTotal(const Double_t *x,const Double_t *par)
{
 //
 // function (it doesn't use the function integral->pay attention to it)
 //
 SetInvMassParameters(par);
 //return ((EvaluateCDFInvMassSigDistr(x[0])))*par[9] + (EvaluateCDFInvMassBkgDistr(x[0]))*par[10];
 return (((fFCN->EvaluateCDFInvMassSigDistr(x[0])))*par[9] + (1-par[9])*(fFCN->EvaluateCDFInvMassBkgDistr(x[0])))*par[10];
}
//__________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCNfitter::CDFResolutionFunction(const Double_t *x,const Double_t *par)
{
 //
 // function 
 //
//  Double_t resParam[kPseudo];
//  for(Int_t i=0; i<kPseudo; i++) {
//   resParam[i] = par[i];
//  }
 fFCN->SetResolutionConstants(par,(Int_t)par[9]);
 return par[10]*fFCN->ResolutionFunc(x[0],(Int_t)(par[9]));
}
//__________________________________________________________________________________________
Double_t AliDielectronBtoJPSItoEleCDFfitFCNfitter::PsProperBackFitFunc(const Double_t* x, const Double_t* par)
{	
 //
 // function 
 //
 Int_t type = TMath::Nint(x[1]);
 Double_t params[7];
 for(Int_t i=0; i<7; i++) params[i]=par[i]; 
 SetPseudoProperBkgParameters(params);		
 return fFCN->EvaluateCDFDecayTimeBkgDistr(x[0],type)*par[7]*fFCN->GetResWeight(type);
}
