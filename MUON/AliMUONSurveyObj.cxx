/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------
/// \class AliMUONSurveyObj
/// Base class for the survey processing of the ALICE DiMuon spectrometer 
///
/// This base object provides methods to process the survey+photogrammetry
/// data of the Chambers (frames) and Detection Elements of the DiMuon
/// Spectrometer and calculate their misalignments. 
///
/// \author Javier Castillo
//-----------------------------------------------------------------------------

#include <fstream>

#include "TMath.h"
#include "TVector3.h"
#include "TGeoMatrix.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TString.h"
#include "TH2.h"
#include "TF2.h"
#include "TGraph2DErrors.h"
#include "TArrayD.h"

#include "AliLog.h"
#include "AliSurveyPoint.h"

#include "AliMUONSurveyObj.h"
#include "AliMUONSurveyUtil.h"

void SurveyFcn(int &npar, double *g, double &f, double *par, int iflag);

/// \cond CLASSIMP
ClassImp(AliMUONSurveyObj)
/// \endcond

AliMUONSurveyObj::AliMUONSurveyObj() 
  : TObject() 
  , fSTargets(0x0)
  , fGBTargets(0x0)
  , fLBTargets(0x0)
  , fLocalTrf(0x0)
  , fAlignTrf(0x0)
  , fBaseTrf(0x0)
  , fOwnerLocalTrf(kFALSE)
  , fOwnerAlignTrf(kTRUE)
  , fOwnerBaseTrf(kFALSE)
  , fUseCM(kTRUE)
  , fPlane(0x0)
  , fFitter(0x0)
  , fXMin(-400.)
  , fXMax(400.)
  , fYMin(-400.)
  , fYMax(400.)
  , fZMin(-2000.)
  , fZMax(2000.)
{
/// Default constructor

  fSTargets = new TObjArray();  
  fSTargets->SetOwner(kFALSE);
  fGBTargets = new TObjArray();  
  fGBTargets->SetOwner(kFALSE);
  fLBTargets = new TObjArray();  
  fLBTargets->SetOwner(kFALSE);  

  fAlignTrf = new TGeoCombiTrans();

  fFitter = new TFitter(100);
}

AliMUONSurveyObj::~AliMUONSurveyObj() {
  /// Destructor
  if(fSTargets) {
    fSTargets->Delete();
    fSTargets = 0x0;
  }
  if(fGBTargets) {
    fGBTargets->Delete();
    fGBTargets = 0x0;
  }
  if(fLBTargets) {
    fLBTargets->Delete();
    fLBTargets = 0x0;
  }
  if (fPlane) {
    fPlane->Delete();
    fPlane = 0x0;
  }
  if(fOwnerLocalTrf && fLocalTrf) {
    fLocalTrf->Delete();
    fLocalTrf = 0x0;
  }
  if(fOwnerAlignTrf && fAlignTrf) {
    fAlignTrf->Delete();
    fAlignTrf = 0x0;
  }
  if(fOwnerBaseTrf && fBaseTrf) {
    fBaseTrf->Delete();
    fBaseTrf = 0x0;
  }
  if (fFitter){
    fFitter->Delete();
    fFitter = 0x0;
  }
}

void AliMUONSurveyObj::AddStickerTarget(AliSurveyPoint *stPoint){
  /// Add sticker target
  if (fUseCM) {
    fSTargets->Add(ConvertPointUnits(stPoint,0.1));
  } else {
    fSTargets->Add(stPoint);
  }
}

void AliMUONSurveyObj::AddGButtonTarget(AliSurveyPoint *btPoint){
  /// Add global button target
  if (fUseCM) {
    fGBTargets->Add(ConvertPointUnits(btPoint,0.1));
  } else {
    fGBTargets->Add(btPoint);
  }  
}

void AliMUONSurveyObj::AddLButtonTarget(AliSurveyPoint *btPoint){
  /// Add local button target target; AliSurveyPoint
  if (fUseCM) {
    fLBTargets->Add(ConvertPointUnits(btPoint,0.1));
  } else {
    fLBTargets->Add(btPoint);
  }  
}

void AliMUONSurveyObj::AddLButtonTarget(TVector3 *btVector){
  /// Add local button target target; TVector3
  fLBTargets->Add(btVector);
}

Int_t AliMUONSurveyObj::AddStickerTargets(TObjArray *pArray, TString stBaseName, Int_t lTargetMax){
  /// Add a maximum of lTargetMax sticker targets with stBaseName from pArray 
  if (!pArray) {
    AliError(Form("Survey points array is empty %p!",pArray));
    return 0;
  }
  if (stBaseName.IsNull()) {
    AliError(Form("Need base name for sticker targets %s!",stBaseName.Data()));
    return 0;
  }
  
  Int_t stIndex = 0;
  AliSurveyPoint *pointSST = 0x0;

  TString stNumber;
  
  for (int iPoint=0; iPoint<lTargetMax; iPoint++) {
    TString stFullName(stBaseName);
    stNumber = Form("%d",iPoint+1);
    if(lTargetMax>9&&iPoint+1<10) {
      stFullName+="0";
    }
    stFullName+=stNumber;

    pointSST = (AliSurveyPoint *)pArray->FindObject(stFullName.Data());

    if(pointSST) {
      AddStickerTarget(pointSST);
      AliInfo(Form("Added survey sticker target %s at index %d",pointSST->GetName(),stIndex));
      stIndex++;
    }
  }

  AliInfo(Form("Found %d sticker targets with base name %s",fSTargets->GetEntries(),stBaseName.Data()));
  return stIndex;
}

Int_t AliMUONSurveyObj::AddGButtonTargets(TObjArray *pArray, TString btBaseName, Int_t lTargetMax){
  /// Add a maximum of lTargetMax global button targets with stBaseName from pArray 
  printf("%s \n",btBaseName.Data());
  if (!pArray) {
    AliError(Form("Survey points array is empty %p!",pArray));
    return 0;
  }
  if (btBaseName.IsNull()) {
    AliError(Form("Need base name for button targets %s!",btBaseName.Data()));
    return 0;
  }
  
  Int_t btIndex = 0;
  AliSurveyPoint *pointSBT = 0x0;

  TString btNumber;
  
  for (int iPoint=0; iPoint<lTargetMax; iPoint++) {
    TString btFullName(btBaseName);
    btNumber = Form("%d",iPoint+1);
    if(lTargetMax>9&&iPoint+1<10) {
      btFullName+="0";
    }
    btFullName+=btNumber;
    printf("%s \n",btFullName.Data());
    pointSBT = (AliSurveyPoint *)pArray->FindObject(btFullName.Data());

    if(pointSBT) {
      AddGButtonTarget(pointSBT);
      AliInfo(Form("Added survey button target %s at index %d",pointSBT->GetName(),btIndex));
      btIndex++;
    }
  }

  AliInfo(Form("Found %d button targets with base name %s",fGBTargets->GetEntries(),btBaseName.Data()));
  return btIndex;
}

Int_t AliMUONSurveyObj::AddLButtonTargets(TObjArray *pArray, TString btBaseName, Int_t lTargetMax){
  /// Add a maximum of lTargetMax local button targets with stBaseName from pArray 
  printf("%s \n",btBaseName.Data());
  if (!pArray) {
    AliError(Form("Local points array is empty %p!",pArray));
    return 0;
  }
  if (btBaseName.IsNull()) {
    AliError(Form("Need base name for button targets %s!",btBaseName.Data()));
    return 0;
  }
  
  Int_t btIndex = 0;
  AliSurveyPoint *pointSBT = 0x0;

  TString btNumber;
  
  for (int iPoint=0; iPoint<lTargetMax; iPoint++) {
    TString btFullName(btBaseName);
    btNumber = Form("%d",iPoint+1);
    if(lTargetMax>9&&iPoint+1<10) {
      btFullName+="0";
    }
    btFullName+=btNumber;
    printf("%s \n",btFullName.Data());
    pointSBT = (AliSurveyPoint *)pArray->FindObject(btFullName.Data());

    if(pointSBT) {
      AddLButtonTarget(pointSBT);
      AliInfo(Form("Added local button target %s at index %d",pointSBT->GetName(),btIndex));
      btIndex++;
    }
  }

  AliInfo(Form("Found %d local button targets with base name %s",fLBTargets->GetEntries(),btBaseName.Data()));
  return btIndex;
}

Int_t AliMUONSurveyObj::GetNStickerTargets() {
  /// return number of sticker targets
  return fSTargets->GetEntriesFast();
}

AliSurveyPoint* AliMUONSurveyObj::GetStickerTarget(Int_t stIndex){
  /// return sticker target at stIndex
  if (stIndex<0||stIndex>=fSTargets->GetEntriesFast()) {
    AliError(Form("No sticker target at index %d",stIndex));
    return 0x0;
  }
  else {
    return (AliSurveyPoint*)fSTargets->At(stIndex);
  }
}

Int_t AliMUONSurveyObj::GetNGButtonTargets() {
  /// return number of global button targets
  return fGBTargets->GetEntriesFast();
}

AliSurveyPoint* AliMUONSurveyObj::GetGButtonTarget(Int_t btIndex){
  /// return global button target at btIndex
  if (btIndex<0||btIndex>=fGBTargets->GetEntriesFast()) {
    AliError(Form("No surveyed button target at index %d",btIndex));
    return 0x0;
  }
  else {
    return (AliSurveyPoint*)fGBTargets->At(btIndex);
  }
}

Int_t AliMUONSurveyObj::GetNLButtonTargets() {
  /// return number of local button targets
  return fGBTargets->GetEntriesFast();
}

AliSurveyPoint* AliMUONSurveyObj::GetLButtonTarget(Int_t btIndex){
  /// return local button target at btIndex
  if (btIndex<0||btIndex>=fLBTargets->GetEntriesFast()) {
    AliError(Form("No surveyed button target at index %d",btIndex));
    return 0x0;
  }
  else {
    if(fLBTargets->At(btIndex)->IsA()==TVector3::Class()){
      TVector3 *lBT = (TVector3*)fLBTargets->At(btIndex);
      TString str("B");
      return (new AliSurveyPoint(TString("local"),(float)lBT->X(),(float)lBT->Y(),(float)lBT->Z(),(float)0.,(float)0.,(float)0.,'B',kTRUE));
    } else if(fLBTargets->At(btIndex)->IsA()==AliSurveyPoint::Class()) {
      AliSurveyPoint *lBT = (AliSurveyPoint*)fLBTargets->At(btIndex);
      return lBT;
    } else {
      AliError(Form("Unexpected class %s ! Valid classes are TVector3 or AliSurveyPoint",fLBTargets->At(btIndex)->ClassName()));
      return 0;
    }
  }
}

void AliMUONSurveyObj::SetPlane(TString pName, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax){
  /// Set the plane function for the plane fitting
  if(fPlane) {
    fPlane->Delete();
    fPlane = 0x0;
  }
  fPlane = new TF2(pName,this,&AliMUONSurveyObj::EqPlane,xMin,xMax,yMin,yMax,3,"AliMUONSurveyObj","EqPlane");
}

void AliMUONSurveyObj::SetPlaneParameters(Double_t p0, Double_t p1, Double_t p2) {
  /// Set the parameters of plane function for the plane fitting
  if (!fPlane) {
    AliError("Must use SetPlane before SetPlaneParameters!!!");
  }
  else {
    fPlane->SetParameter(0,p0);
    fPlane->SetParameter(1,p1);
    fPlane->SetParameter(2,p2);
  }
}

void AliMUONSurveyObj::DrawSTargets() {
  /// Draw a graph of the sticker targets
  TGraph2DErrors *gST = new TGraph2DErrors(3);
  AliSurveyPoint *pST = 0x0;
  for (Int_t iPoint=0; iPoint<GetNStickerTargets(); iPoint++) {
    pST = GetStickerTarget(iPoint);
    //    pST->PrintPoint();
    gST->SetPoint(iPoint,pST->GetX(),pST->GetY(),pST->GetZ());
    gST->SetPointError(iPoint,pST->GetPrecisionX(),pST->GetPrecisionY(),pST->GetPrecisionZ());
  }
  gST->DrawClone("P0");

  if (gST) gST->Delete();
}

Double_t AliMUONSurveyObj::FitPlane() {
  /// Fit plane to sticker targets
  if (!fPlane) {
    AliError("Must use SetPlane before FitPlane!!!");
    return 0.;
  }
  if (fSTargets->GetEntriesFast()<3) {
    AliError("Not enough sticker targets (%d) for plane fitting!!!");
    return 0.;
  }

  Double_t pl[3] = {0};
  Double_t pg[3] = {0};
    
  TGraph2DErrors *gST = new TGraph2DErrors(3);
  AliSurveyPoint *pST = 0x0;
  for (Int_t iPoint=0; iPoint<GetNStickerTargets(); iPoint++) {
    pST = GetStickerTarget(iPoint);
    //    pST->PrintPoint();
    pg[0] =  pST->GetX(); 
    pg[1] =  pST->GetY(); 
    pg[2] =  pST->GetZ(); 
    fBaseTrf->MasterToLocal(pg,pl);
    gST->SetPoint(iPoint,pl[0],pl[1],pl[2]);
    printf("%d %f %f %f\n",iPoint,pl[0],pl[1],pl[2]);
    gST->SetPointError(iPoint,pST->GetPrecisionX(),pST->GetPrecisionY(),pST->GetPrecisionZ());
  }
  gST->Fit(fPlane);

  if (gST) gST->Delete();

  return fPlane->GetChisquare();
}

Double_t AliMUONSurveyObj::SurveyChi2(Double_t *par){
  /// Returns the chisquare between local2global transform of local button targets and their surveyed position
  TGeoTranslation transTemp;
  TGeoRotation rotTemp;
  TGeoCombiTrans trfTemp;

  Double_t lChi2=0.;

  trfTemp.SetTranslation(transTemp);
  trfTemp.SetRotation(rotTemp);
  trfTemp.Clear();
  trfTemp.RotateZ(TMath::RadToDeg()*par[5]);
  trfTemp.RotateY(TMath::RadToDeg()*par[4]);
  trfTemp.RotateX(TMath::RadToDeg()*par[3]);
  trfTemp.SetTranslation(par[0],par[1],par[2]);

  TGeoHMatrix matGlo = (*fBaseTrf)*trfTemp;
  TGeoCombiTrans trfGlo(matGlo);
		
  Double_t pl[3] = {0};
  Double_t pg[3] = {0};
  
  for(Int_t iPoint=0; iPoint<fGBTargets->GetEntries(); iPoint++){
    AliSurveyPoint *gBT = (AliSurveyPoint*)fGBTargets->At(iPoint);
    if(fLBTargets->At(iPoint)->IsA()==TVector3::Class()){
      TVector3 *lBT = (TVector3*)fLBTargets->At(iPoint);
      pl[0]=lBT->X();
      pl[1]=lBT->Y();
      pl[2]=lBT->Z();
    } else if(fLBTargets->At(iPoint)->IsA()==AliSurveyPoint::Class()) {
      AliSurveyPoint *lBT = (AliSurveyPoint*)fLBTargets->At(iPoint);
      pl[0]=lBT->GetX();
      pl[1]=lBT->GetY();
      pl[2]=lBT->GetZ();
    } else {
      AliError(Form("Unexpected class %s ! Valid classes are TVector3 or AliSurveyPoint",fLBTargets->At(iPoint)->ClassName()));
      return 0;
    }

    trfGlo.LocalToMaster(pl, pg);
    //    printf("%d %f %f %f\n",iPoint,pg[0],pg[1],pg[2]);
    if(fLBTargets->At(iPoint)->IsA()==TVector3::Class()){
      lChi2 += (pg[0]-gBT->GetX())*(pg[0]-gBT->GetX())/(gBT->GetPrecisionX()*gBT->GetPrecisionX());
      lChi2 += (pg[1]-gBT->GetY())*(pg[1]-gBT->GetY())/(gBT->GetPrecisionY()*gBT->GetPrecisionY());
      lChi2 += (pg[2]-gBT->GetZ())*(pg[2]-gBT->GetZ())/(gBT->GetPrecisionZ()*gBT->GetPrecisionZ());
    } else if(fLBTargets->At(iPoint)->IsA()==AliSurveyPoint::Class()) {
      AliSurveyPoint *lBT = (AliSurveyPoint*)fLBTargets->At(iPoint);
      lChi2 += (pg[0]-gBT->GetX())*(pg[0]-gBT->GetX())/(gBT->GetPrecisionX()*gBT->GetPrecisionX()+lBT->GetPrecisionX()*lBT->GetPrecisionX());
      lChi2 += (pg[1]-gBT->GetY())*(pg[1]-gBT->GetY())/(gBT->GetPrecisionY()*gBT->GetPrecisionY()+lBT->GetPrecisionY()*lBT->GetPrecisionY());
      lChi2 += (pg[2]-gBT->GetZ())*(pg[2]-gBT->GetZ())/(gBT->GetPrecisionZ()*gBT->GetPrecisionZ()+lBT->GetPrecisionZ()*lBT->GetPrecisionZ());
    } else {
      AliError(Form("Unexpected class %s ! Valid classes are TVector3 or AliSurveyPoint",fLBTargets->At(iPoint)->ClassName()));
      return 0;
    }
  }

  return lChi2;
}

//_____________________________________________________________________________
void SurveyFcn(int &npar, double *g, double &f, double *par, int iflag) {
  /// 
  /// Standard function as needed by Minuit-like minimization procedures. 
  /// For the set of parameters par calculates and returns chi-squared.
  ///

  // smuggle a C++ object into a C function
  AliMUONSurveyObj *aSurveyObj = (AliMUONSurveyObj*) gMinuit->GetObjectFit(); 

  f = aSurveyObj->SurveyChi2(par);
  if (iflag==3) {}
  if (npar) {} 
  if (g) {} // no warnings about unused stuff...

}

//_____________________________________________________________________________
Int_t AliMUONSurveyObj::SurveyToAlign(TGeoCombiTrans &quadTransf, Double_t *parErr, Double_t psi, Double_t tht, Double_t epsi, Double_t etht) {
  /// Main function to obtain the misalignments from the surveyed position of the button targets; 
  if (fGBTargets->GetEntries()!=fLBTargets->GetEntries()){
    AliError(Form("Different number of button targets: %d survey points and %d local coord!",
		  fGBTargets->GetEntries(),fLBTargets->GetEntries()));
    return 0;
  }

  TFitter fitter(100);
  gMinuit->SetObjectFit(this);
  fitter.SetFCN(SurveyFcn);
  fitter.SetParameter(0,"dx",0,0.1,-20,20);
  fitter.SetParameter(1,"dy",0,0.1,-20,20);
  fitter.SetParameter(2,"dz",0,0.1,0,0);
  fitter.SetParameter(3,"rx",psi,0.0001,-0.05,0.05);
  fitter.SetParameter(4,"ry",tht,0.0001,-0.05,0.05);
//   fitter.SetParameter(3,"rx",psi,0.0001,psi-5*epsi,psi+5*epsi);
//   fitter.SetParameter(4,"ry",tht,0.0001,tht-5*etht,tht+5*etht);
  fitter.SetParameter(5,"rz",0,0.0001,0,0);

  if(psi) fitter.FixParameter(3);
  if(tht) fitter.FixParameter(4);

  double arglist[100];
  arglist[0] = 2;
  fitter.ExecuteCommand("SET PRINT", arglist, 1);
  fitter.ExecuteCommand("SET ERR", arglist, 1);
  arglist[0]=0;
  //fitter.ExecuteCommand("SIMPLEX", arglist, 1);
  //  fitter.ExecuteCommand("MINIMIZE", arglist, 1);
  fitter.ExecuteCommand("MIGRAD", arglist, 1);
  fitter.ExecuteCommand("IMPROVE", arglist, 1);
  //  fitter.ExecuteCommand("MINOS", arglist, 1);
  //  fitter.ExecuteCommand("CALL 3", arglist,0);

  for (int j=0; j<3; j++) printf("%10.3f ",fitter.GetParameter(j));   
  for (int j=3; j<6; j++) printf("%10.3f ",1000*fitter.GetParameter(j));   
  printf("\n");
  for (int j=0; j<3; j++) printf("%10.3f ",fitter.GetParError(j));
  for (int j=3; j<6; j++) printf("%10.3f ",1000*fitter.GetParError(j));
  printf("\n");

  quadTransf.Clear();
  quadTransf.RotateZ(TMath::RadToDeg()*fitter.GetParameter(5));
  quadTransf.RotateY(TMath::RadToDeg()*fitter.GetParameter(4));
  quadTransf.RotateX(TMath::RadToDeg()*fitter.GetParameter(3));
  quadTransf.SetTranslation(fitter.GetParameter(0),fitter.GetParameter(1),fitter.GetParameter(2));

  for(Int_t iPar=0; iPar<6; iPar++){
    parErr[iPar] = fitter.GetParError(iPar);
  }
  if(epsi) parErr[3] = epsi;
  if(etht) parErr[4] = etht;

  return 1;

}

//_____________________________________________________________________________
Int_t AliMUONSurveyObj::SurveyToAlign(Double_t psi, Double_t tht, Double_t epsi, Double_t etht) {
  /// Main function to obtain the misalignments from the surveyed position of the button targets; 
  if (fGBTargets->GetEntries()!=fLBTargets->GetEntries()){
    AliError(Form("Different number of button targets: %d survey points and %d local coord!",
		  fGBTargets->GetEntries(),fLBTargets->GetEntries()));
    return 0;
  }

  //  TFitter fitter(100);
  gMinuit->SetObjectFit(this);
  fFitter->SetFCN(SurveyFcn);
  fFitter->SetParameter(0,"dx",0,0.1,-20,20);
  fFitter->SetParameter(1,"dy",0,0.1,-20,20);
  fFitter->SetParameter(2,"dz",0,0.1,0,0);
  if(psi)
      fFitter->SetParameter(3,"rx",psi,epsi,psi-5*epsi,psi+5*epsi);
  else
    fFitter->SetParameter(3,"rx",psi,0.0001,-0.05,0.05);
  if(tht)
    fFitter->SetParameter(4,"ry",tht,etht,tht-5*etht,tht+5*etht);
  else
    fFitter->SetParameter(4,"ry",tht,0.0001,-0.05,0.05);
//   fFitter->SetParameter(3,"rx",psi,0.0001,psi-5*epsi,psi+5*epsi);
//   fFitter->SetParameter(4,"ry",tht,0.0001,tht-5*etht,tht+5*etht);
  fFitter->SetParameter(5,"rz",0,0.0001,0,0);

  if(psi) fFitter->FixParameter(3);
  if(tht) fFitter->FixParameter(4);

  double arglist[100];
  arglist[0] = 2;
  fFitter->ExecuteCommand("SET PRINT", arglist, 1);
  fFitter->ExecuteCommand("SET ERR", arglist, 1);
  arglist[0]=0;
  //fFitter->ExecuteCommand("SIMPLEX", arglist, 1);
  //  fFitter->ExecuteCommand("MINIMIZE", arglist, 1);
  fFitter->ExecuteCommand("MIGRAD", arglist, 1);
  fFitter->ExecuteCommand("IMPROVE", arglist, 1);
//   fFitter->ExecuteCommand("MINOS", arglist, 1);
//   fFitter->ExecuteCommand("CALL 3", arglist,0);

  for (int j=0; j<3; j++) printf("%10.3f ",fFitter->GetParameter(j));   
  for (int j=3; j<6; j++) printf("%10.3f ",1000*fFitter->GetParameter(j));   
  printf("\n");
  for (int j=0; j<3; j++) printf("%10.3f ",fFitter->GetParError(j));
  for (int j=3; j<6; j++) printf("%10.3f ",1000*fFitter->GetParError(j));
  printf("\n");

  fAlignTrf->Clear();
  fAlignTrf->RotateZ(TMath::RadToDeg()*fFitter->GetParameter(5));
  fAlignTrf->RotateY(TMath::RadToDeg()*fFitter->GetParameter(4));
  fAlignTrf->RotateX(TMath::RadToDeg()*fFitter->GetParameter(3));
  fAlignTrf->SetTranslation(fFitter->GetParameter(0),fFitter->GetParameter(1),fFitter->GetParameter(2));

  if(epsi) fFitter->ReleaseParameter(3); // To get error
  if(etht) fFitter->ReleaseParameter(4); // To get error

  TGeoCombiTrans lGlobalTrf = TGeoCombiTrans((*fBaseTrf)*(*fAlignTrf));
  AliSurveyPoint *pointGBT;
  AliSurveyPoint *pointLBT;
  Double_t pl[3] = {0};
  Double_t pg[3] = {0};
  for (int iPoint=0; iPoint<GetNGButtonTargets(); iPoint++){
    pointGBT=GetGButtonTarget(iPoint);
    pointLBT=GetLButtonTarget(iPoint);
    pl[0] = pointLBT->GetX();
    pl[1] = pointLBT->GetY();
    pl[2] = pointLBT->GetZ();
    lGlobalTrf.LocalToMaster(pl,pg);
    printf("Point %d  local: %.3f %.3f %.3f\n",iPoint,pl[0],pl[1],pl[2]);
    printf("Point %d global: %.3f %.3f %.3f\n",iPoint,pg[0],pg[1],pg[2]);
    printf("Point %d survey: %.3f %.3f %.3f\n",iPoint,pointGBT->GetX(),pointGBT->GetY(),pointGBT->GetZ());
  }

  return 1;

}

Double_t AliMUONSurveyObj::EvalFunction(const TF2 *lFunction, Int_t iP1, Int_t iP2, const Char_t *lCoord) {
  /// Evaluate the given function at the given points for the given coordinate
  if (!lFunction) {
    AliError("No function given!!!");
    return 0;
  }
  AliSurveyPoint *gP1 = GetGButtonTarget(iP1);
  AliSurveyPoint *gP2 = GetGButtonTarget(iP2);
  
  if(!gP1||!gP2){
    AliError("Missing global button target!!!");
    return 0;
  }

  //  AliInfo(Form("Function %s parameters %f %f %f %f %f %f",lFunction->GetName(),lFunction->GetParameter(0),lFunction->GetParameter(1),lFunction->GetParameter(2),lFunction->GetParameter(3),lFunction->GetParameter(4),lFunction->GetParameter(5)));
  Double_t pl1[3] = {0};
  Double_t pl2[3] = {0};
  Double_t pg1[3] = {0};
  Double_t pg2[3] = {0};
    
  pg1[0] =  gP1->GetX(); 
  pg1[1] =  gP1->GetY(); 
  pg1[2] =  gP1->GetZ(); 
  pg2[0] =  gP2->GetX(); 
  pg2[1] =  gP2->GetY(); 
  pg2[2] =  gP2->GetZ(); 

  fBaseTrf->MasterToLocal(pg1,pl1);
  fBaseTrf->MasterToLocal(pg2,pl2);

  Double_t lVal = 0.;
  switch (lCoord[0]) {
  case 'X':
    {
      lVal = lFunction->Eval(pl1[0],pl2[0]);
      //      lVal = lFunction->Eval(gP1->GetX(),gP2->GetX());
      //      AliInfo(Form("case X, lVal = %f",lVal));
      return lVal;
    }
  case 'Y':
    {
      lVal = lFunction->Eval(pl1[1],pl2[1]);
      //      lVal = lFunction->Eval(gP1->GetY(),gP2->GetY());
      //      AliInfo(Form("case Y, lVal = %f",lVal));
      return lVal;
    }
  case 'Z':
    {
      lVal = lFunction->Eval(pl1[2],pl2[2]);
      //      lVal = lFunction->Eval(gP1->GetZ(),gP2->GetZ());
      //      AliInfo(Form("case Z, lVal = %f",lVal));
      return lVal;
    }
  default:
    {
      AliError(Form("Coordinate %c is not valid, options are X Y Z",lCoord));
      return 0;
    }
  }
}

void AliMUONSurveyObj::CalculateTranslation(TF2 *xFunc, TF2 *yFunc, TF2 *zFunc, Int_t iP1, Int_t iP2, Double_t *lCenTemp) {
  /// Calculate the center translation using analytic functions
  lCenTemp[0] = EvalFunction(xFunc,iP1,iP2,"X");
  lCenTemp[1] = EvalFunction(yFunc,iP1,iP2,"Y");
  lCenTemp[2] = EvalFunction(zFunc,iP1,iP2,"Z");

}

Double_t AliMUONSurveyObj::CalculateGlobalDiff(TGeoCombiTrans &lTransf, Int_t nPoints, TArrayD &lDiff){
  /// By hand computation of distance between local2global transform of target position and its surveyed position
  if (nPoints > GetNGButtonTargets()) {
    nPoints = GetNGButtonTargets();
  }

  for(Int_t iVal=0; iVal<nPoints*(3+1)+1; iVal++){
    lDiff[iVal] = 0.;
  }

  Double_t pl[3] = {0};
  Double_t pg[3] = {0};
  Double_t pml[3] = {0};
  Double_t pmg[3] = {0};
  AliSurveyPoint *gBT = 0x0;
  for(Int_t iPoint=0; iPoint<nPoints; iPoint++){
    gBT = GetGButtonTarget(iPoint);
    if(!gBT||!fLBTargets->At(iPoint)){
      AliError(Form("The local or global target %d is missing!",iPoint));
      lDiff[nPoints*(3+1)] = 1.e7;
      return lDiff[nPoints*(3+1)];
    }
    if(fLBTargets->At(iPoint)->IsA()==TVector3::Class()){
      TVector3 *lBT = (TVector3*)fLBTargets->At(iPoint);
      pl[0]=lBT->X();
      pl[1]=lBT->Y();
      pl[2]=lBT->Z();
    } else if(fLBTargets->At(iPoint)->IsA()==AliSurveyPoint::Class()) {
      AliSurveyPoint *lBT = (AliSurveyPoint*)fLBTargets->At(iPoint);
      pl[0]=lBT->GetX();
      pl[1]=lBT->GetY();
      pl[2]=lBT->GetZ();
    } else {
      AliError(Form("Unexpected class %s ! Valid classes are TVector3 or AliSurveyPoint",fLBTargets->At(iPoint)->ClassName()));
      return 0;
    }

    lTransf.LocalToMaster(pl,pg);
    pmg[0] = gBT->GetX();
    pmg[1] = gBT->GetY();
    pmg[2] = gBT->GetZ();
    fBaseTrf->MasterToLocal(pmg,pml);
//     printf("l %d %f %f %f\n",iPoint,pl[0],pl[1],pl[2]);
//     printf("g %d %f %f %f\n",iPoint,pg[0],pg[1],pg[2]);
//     printf("ml %d %f %f %f\n",iPoint,pml[0],pml[1],pml[2]);
//     printf("mg %d %f %f %f\n",iPoint,gBT->GetX(),gBT->GetY(),gBT->GetZ());
    lDiff[iPoint*(3+1)+0] = (pml[0]-pg[0]);
    lDiff[iPoint*(3+1)+1] = (pml[1]-pg[1]);
    lDiff[iPoint*(3+1)+2] = (pml[2]-pg[2]);
    
    lDiff[iPoint*(3+1)+3] = TMath::Sqrt(lDiff[iPoint*(3+1)+0]*lDiff[iPoint*(3+1)+0]+
					lDiff[iPoint*(3+1)+1]*lDiff[iPoint*(3+1)+1]+
					lDiff[iPoint*(3+1)+2]*lDiff[iPoint*(3+1)+2]);

    lDiff[nPoints*(3+1)] += lDiff[iPoint*(3+1)+3]*lDiff[iPoint*(3+1)+3];
  }
		
  lDiff[nPoints*(3+1)] = TMath::Sqrt(lDiff[nPoints*(3+1)]);
  return lDiff[nPoints*(3+1)];
}

Int_t AliMUONSurveyObj::CalculateBestTransf(Int_t iP1, Int_t iP2, Double_t *lXYZ, Double_t *lPTP) {
  /// By hand calculation of the best local to global transform using 2 button targets
  Double_t lPsi = lPTP[0];
  Double_t lTht = lPTP[1];

  Double_t pl1[3] = {0};
  Double_t pl2[3] = {0};

  if(!fLBTargets->At(iP1)||!fLBTargets->At(iP2)){
    AliError(Form("Local target %d or %d is missing!",iP1,iP2));
    return 0;
  }

  if(fLBTargets->At(iP1)->IsA()==TVector3::Class()){
    TVector3 *lBT1 = (TVector3*)fLBTargets->At(iP1);
    pl1[0]=lBT1->X();
    pl1[1]=lBT1->Y();
    pl1[2]=lBT1->Z();
  } else if(fLBTargets->At(iP1)->IsA()==AliSurveyPoint::Class()) {
    AliSurveyPoint *lBT1 = (AliSurveyPoint*)fLBTargets->At(iP1);
    pl1[0]=lBT1->GetX();
    pl1[1]=lBT1->GetY();
    pl1[2]=lBT1->GetZ();
  } else {
    AliError(Form("Unexpected class %s ! Valid classes are TVector3 or AliSurveyPoint",fLBTargets->At(iP1)->ClassName()));
    return 0;
  }
  if(fLBTargets->At(iP2)->IsA()==TVector3::Class()){
    TVector3 *lBT2 = (TVector3*)fLBTargets->At(iP2);
    pl2[0]=lBT2->X();
    pl2[1]=lBT2->Y();
    pl2[2]=lBT2->Z();
  } else if(fLBTargets->At(iP2)->IsA()==AliSurveyPoint::Class()) {
    AliSurveyPoint *lBT2 = (AliSurveyPoint*)fLBTargets->At(iP2);
    pl2[0]=lBT2->GetX();
    pl2[1]=lBT2->GetY();
    pl2[2]=lBT2->GetZ();
  } else {
    AliError(Form("Unexpected class %s ! Valid classes are TVector3 or AliSurveyPoint",fLBTargets->At(iP2)->ClassName()));
    return 0;
  }

  
  AliMUONSurveyUtil *surveyUtil = AliMUONSurveyUtil::Instance();

  // Xcenter functions
  const char *fxcName = "fXcn00"; 
  TF2 **fXc = new TF2*[2];
  fxcName = "fXcn";
  fXc[0] = new TF2(fxcName,surveyUtil,&AliMUONSurveyUtil::XnCenter,fXMin,fXMax,fYMin,fYMax,7,"AliMUONSurveyUtil","XnCenter");
  fxcName = "fXcp";
  fXc[1] = new TF2(fxcName,surveyUtil,&AliMUONSurveyUtil::XpCenter,fXMin,fXMax,fYMin,fYMax,7,"AliMUONSurveyUtil","XpCenter");

  // Ycenter functions
  const char *fycName = "fYcn00"; 
  TF2 **fYc = new TF2*[2];
  fycName = "fYcn";
  fYc[0] = new TF2(fycName,surveyUtil,&AliMUONSurveyUtil::YnCenter,fYMin,fYMax,fYMin,fYMax,8,"AliMUONSurveyUtil","YnCenter");
  fycName = "fYcp";
  fYc[1] = new TF2(fycName,surveyUtil,&AliMUONSurveyUtil::YpCenter,fYMin,fYMax,fYMin,fYMax,8,"AliMUONSurveyUtil","YpCenter");   

  // Zcenter functions
  const char *fzcName = "fZcn00"; 
  TF2 **fZc = new TF2*[2];
  fzcName = "fZcn";
  fZc[0] = new TF2(fzcName,surveyUtil,&AliMUONSurveyUtil::ZnCenter,fZMin,fZMax,fZMin,fZMax,8,"AliMUONSurveyUtil","ZnCenter");
  fzcName = "fZcp";
  fZc[1] = new TF2(fzcName,surveyUtil,&AliMUONSurveyUtil::ZpCenter,fZMin,fZMax,fZMin,fZMax,8,"AliMUONSurveyUtil","ZpCenter");   

  // Phi rotation using xglobal coords functions
  const char *fphixName = "fPhiXnn00"; 
  TF2 ***fPhiX = new TF2**[2];
  for (Int_t iX =0; iX<2; iX++) {
    fPhiX[iX] = new TF2*[2];
  }
  fphixName = "fPhiXnn";
  fPhiX[0][0] = new TF2(fphixName,surveyUtil,&AliMUONSurveyUtil::PhiXnn,fXMin,fXMax,fXMin,fXMax,7,"AliMUONSurveyUtil","PhiXnn");
  fphixName = "fPhiXnp";
  fPhiX[0][1] = new TF2(fphixName,surveyUtil,&AliMUONSurveyUtil::PhiXnp,fXMin,fXMax,fXMin,fXMax,7,"AliMUONSurveyUtil","PhiXnp");   
  fphixName = "fPhiXpn";
  fPhiX[1][0] = new TF2(fphixName,surveyUtil,&AliMUONSurveyUtil::PhiXpn,fXMin,fXMax,fXMin,fXMax,7,"AliMUONSurveyUtil","PhiXpn");
  fphixName = "fPhiXpp";
  fPhiX[1][1] = new TF2(fphixName,surveyUtil,&AliMUONSurveyUtil::PhiXpp,fXMin,fXMax,fXMin,fXMax,7,"AliMUONSurveyUtil","PhiXpp");   

  // Phi rotation using yglobal coords functions
  const char *fphiyName = "fPhiYnn00"; 
  TF2 ***fPhiY = new TF2**[2];
  for (Int_t iY =0; iY<2; iY++) {
    fPhiY[iY] = new TF2*[2];
  }
  fphiyName = "fPhiYnn";
  fPhiY[0][0] = new TF2(fphiyName,surveyUtil,&AliMUONSurveyUtil::PhiYnn,fYMin,fYMax,fYMin,fYMax,8,"AliMUONSurveyUtil","PhiYnn");
  fphiyName = "fPhiYnp";
  fPhiY[0][1] = new TF2(fphiyName,surveyUtil,&AliMUONSurveyUtil::PhiYnp,fYMin,fYMax,fYMin,fYMax,8,"AliMUONSurveyUtil","PhiYnp");   
  fphiyName = "fPhiYpn";
  fPhiY[1][0] = new TF2(fphiyName,surveyUtil,&AliMUONSurveyUtil::PhiYpn,fYMin,fYMax,fYMin,fYMax,8,"AliMUONSurveyUtil","PhiYpn");
  fphiyName = "fPhiYpp";
  fPhiY[1][1] = new TF2(fphiyName,surveyUtil,&AliMUONSurveyUtil::PhiYpp,fYMin,fYMax,fYMin,fYMax,8,"AliMUONSurveyUtil","PhiYpp");   
  

  // Set Parameters of functions
  for(Int_t iS=0; iS<2; iS++){
    fXc[iS]->SetParameters(pl1[0],pl1[1],pl1[2],pl2[0],pl2[1],pl2[2],lTht);
    fYc[iS]->SetParameters(pl1[0],pl1[1],pl1[2],pl2[0],pl2[1],pl2[2],lPsi,lTht);
    fZc[iS]->SetParameters(pl1[0],pl1[1],pl1[2],pl2[0],pl2[1],pl2[2],lPsi,lTht);
//     fXc[iS]->SetParameters(lBT1->X(),lBT1->Y(),lBT1->Z(),lBT2->X(),lBT2->Y(),lBT2->Z(),lTht);
//     fYc[iS]->SetParameters(lBT1->X(),lBT1->Y(),lBT1->Z(),lBT2->X(),lBT2->Y(),lBT2->Z(),lPsi,lTht);
//     fZc[iS]->SetParameters(lBT1->X(),lBT1->Y(),lBT1->Z(),lBT2->X(),lBT2->Y(),lBT2->Z(),lPsi,lTht);
    for(Int_t jS=0; jS<2; jS++){
      fPhiX[iS][jS]->SetParameters(pl1[0],pl1[1],pl1[2],pl2[0],pl2[1],pl2[2],lTht);
      fPhiY[iS][jS]->SetParameters(pl1[0],pl1[1],pl1[2],pl2[0],pl2[1],pl2[2],lPsi,lTht);
//     fPhiX[iS][jS]->SetParameters(lBT1->X(),lBT1->Y(),lBT1->Z(),lBT2->X(),lBT2->Y(),lBT2->Z(),lTht);
//     fPhiY[iS][jS]->SetParameters(lBT1->X(),lBT1->Y(),lBT1->Z(),lBT2->X(),lBT2->Y(),lBT2->Z(),lPsi,lTht);
    }
  }

  Double_t lCenTemp[3];
  Double_t lRotTemp[3];

  TGeoCombiTrans trfTemp;

  Int_t nPoints = GetNGButtonTargets();

  TArrayD lDiffTemp(nPoints*(3+1)+1);
  TArrayD lDiffMin(nPoints*(3+1)+1);

  for(Int_t i=0; i<nPoints*(3+1)+1; i++){
    lDiffMin[i]=1000000.;
    lDiffTemp[i]=0.;
  }

  //
  // Calculate Detection Element Center from button targets
  //	
  
  // Trying 2x*2y*2z*(2phi+2phi) possibilities
  for(Int_t iX=0; iX<2; iX++){
    for(Int_t iY=0; iY<2; iY++){
      for(Int_t iZ=0; iZ<2; iZ++){
	CalculateTranslation(fXc[iX],fYc[iY],fZc[iZ],iP1,iP2,lCenTemp);
	
	lRotTemp[0] = lPsi;
	lRotTemp[1] = lTht;
	for(Int_t iP=0; iP<2; iP++){
	  lRotTemp[2] = EvalFunction(fPhiX[iX][iP],iP1,iP2,"X");
	  
	  trfTemp.Clear();
	  trfTemp.RotateZ(TMath::RadToDeg()*lRotTemp[2]);
	  trfTemp.RotateY(TMath::RadToDeg()*lRotTemp[1]);
	  trfTemp.RotateX(TMath::RadToDeg()*lRotTemp[0]);
	  trfTemp.SetTranslation(lCenTemp);
	  
	  if(CalculateGlobalDiff(trfTemp,nPoints,lDiffTemp)<lDiffMin[nPoints*(3+1)]){
	    printf("Diffs");
	    for(Int_t i=0; i<nPoints*(3+1)+1; i++){
	      printf(" %f",lDiffTemp[i]); 
	    }
	    printf("\n");
	    printf(" : mycenX%dY%dZ%d(%f,%f,%f); rotx%d(%f,%f,%f)\n",iX,iY,iZ,lCenTemp[0],lCenTemp[1],lCenTemp[2],iP,lRotTemp[0],lRotTemp[1],lRotTemp[2]);
	    printf("Transformation improved ...\n");
	    for (int i=0; i<3; i++) {
	      lXYZ[i] = lCenTemp[i];
	    } 
	    lPTP[2] = lRotTemp[2];
	    for(Int_t i=0; i<nPoints*(3+1)+1; i++){
	      lDiffMin[i]=lDiffTemp[i];
	    }
	  }
	}
	for(Int_t iP=0; iP<2; iP++){
	  lRotTemp[2] = EvalFunction(fPhiY[iY][iP],iP1,iP2,"Y");
	  
	  trfTemp.Clear();
	  trfTemp.RotateZ(TMath::RadToDeg()*lRotTemp[2]);
	  trfTemp.RotateY(TMath::RadToDeg()*lRotTemp[1]);
	  trfTemp.RotateX(TMath::RadToDeg()*lRotTemp[0]);
	  trfTemp.SetTranslation(lCenTemp);
	  
	  if(CalculateGlobalDiff(trfTemp,nPoints,lDiffTemp)<lDiffMin[nPoints*(3+1)]){
	    printf("Diffs");
	    for(Int_t i=0; i<nPoints*(3+1)+1; i++){
	      printf(" %f",lDiffTemp[i]); 
	    }
	    printf("\n");
	    printf(" : mycenX%dY%dZ%d(%f,%f,%f); roty%d(%f,%f,%f)\n",iX,iY,iZ,lCenTemp[0],lCenTemp[1],lCenTemp[2],iP,lRotTemp[0],lRotTemp[1],lRotTemp[2]);
	    printf("Transformation improved ...\n");
	    for (int i=0; i<3; i++) {
	      lXYZ[i] = lCenTemp[i];
	    } 
	    lPTP[2] = lRotTemp[2];
	    for(Int_t i=0; i<nPoints*(3+1)+1; i++){
	      lDiffMin[i]=lDiffTemp[i];
	    }
	  }
	}
      }
    }
  }

  for (Int_t i=0; i<2; i++) {
    delete fXc[i];
    delete fYc[i];
    delete fZc[i];
    for (Int_t j=0; j<2; j++) {
      delete fPhiX[i][j];
      delete fPhiY[i][j];
    }
    delete[] fPhiX[i];
    delete[] fPhiY[i];
  }
  delete[] fXc;
  delete[] fYc;
  delete[] fZc;
  delete[] fPhiX;
  delete[] fPhiY;

  if (lDiffMin[nPoints*(3+1)]>20) return 0;

  return 1;
}

void AliMUONSurveyObj::CalculateMeanTransf(Double_t *lXYZ, Double_t *lPTP) {
  /// By hand calculation of the mean (for nPairs of targets) of the best local to global transform using 2 button targets
    Double_t xce=0.;
    Double_t yce=0.;
    Double_t zce=0.;
    Double_t phi=0.;
    
    Int_t nPairs = 0;
    Int_t nPoints = GetNGButtonTargets();
    // Loop over all possible pairs of button tragets
    for(Int_t iP1=0; iP1<nPoints; iP1++){
      for(Int_t iP2=iP1+1; iP2<nPoints; iP2++){
 	printf("%d and %d\n",iP1,iP2);

	if(CalculateBestTransf(iP1,iP2,lXYZ,lPTP)) {
	  nPairs++;
	
	  xce+=lXYZ[0];
	  yce+=lXYZ[1];
	  zce+=lXYZ[2];
	  phi+=lPTP[2];
	}
      }
    }

    if (!nPairs) return;
    
    lXYZ[0]=xce/nPairs;
    lXYZ[1]=yce/nPairs;
    lXYZ[2]=zce/nPairs;
    lPTP[2]=phi/nPairs;
}

void AliMUONSurveyObj::PrintLocalTrf() {
  /// Print the local transformation
  Double_t lRotTemp[3];
  AliMUONSurveyUtil::MatrixToAngles(fLocalTrf->GetRotationMatrix(),lRotTemp);
  printf("(%.3f %.3f %.3f), (%.6f %.6f %.6f)\n",fLocalTrf->GetTranslation()[0],fLocalTrf->GetTranslation()[1],fLocalTrf->GetTranslation()[2],lRotTemp[0],lRotTemp[1],lRotTemp[2]);
}

void AliMUONSurveyObj::PrintAlignTrf() {
  /// Print the alignment transformation
  Double_t lRotTemp[3];
  AliMUONSurveyUtil::MatrixToAngles(fAlignTrf->GetRotationMatrix(),lRotTemp);
  printf("(%.3f %.3f %.3f), (%.6f %.6f %.6f)\n",fAlignTrf->GetTranslation()[0],fAlignTrf->GetTranslation()[1],fAlignTrf->GetTranslation()[2],lRotTemp[0],lRotTemp[1],lRotTemp[2]);
}

void AliMUONSurveyObj::FillSTHistograms(TString baseNameC, TH2 *hSTc, TString baseNameA, TH2 *hSTa) {
  /// Fill sticker target histograms for monitoring
  if(baseNameC.IsNull()||!hSTc){
    AliError("Need base name for points on side C and/or a histogram for them!");
    return;
  }
  AliSurveyPoint *pointST = 0x0;
  for (Int_t iPoint=0; iPoint<GetNStickerTargets(); iPoint++) {
    pointST = GetStickerTarget(iPoint);
    if (!pointST) continue;
    if (pointST->GetPointName().Contains(baseNameC)){
      hSTc->Fill(pointST->GetX(),pointST->GetY(),-pointST->GetZ());
    } else if ((!baseNameA.IsNull()) && 
	      (pointST->GetPointName().Contains(baseNameA))) {
      if (!hSTa){
	AliError("Base name for points on side A provided but no histogram for them!");
	continue;
      }
      hSTa->Fill(pointST->GetX(),pointST->GetY(),-pointST->GetZ());
    }
  }
}

Double_t AliMUONSurveyObj::GetAlignResX() {
  /// Returns the uncertainty of the x translation parameter 
  if(!fFitter) {
    AliError("There is no fitter for this object! X resolution will be 0.");
    return 0.;
  }
  return fFitter->GetParError(0);
}

Double_t AliMUONSurveyObj::GetAlignResY() {
  /// Returns the uncertainty of the y translation parameter 
  if(!fFitter) {
    AliError("There is no fitter for this object! Y resolution will be 0.");
    return 0.;
  }
  return fFitter->GetParError(1);
}

AliSurveyPoint* AliMUONSurveyObj::ConvertPointUnits(AliSurveyPoint *stPoint, Float_t lFactor) {
  /// Return the AliSurveyPoint with new units. Default is from mm -> cm 
  return new AliSurveyPoint(stPoint->GetPointName(),
			    lFactor*stPoint->GetX(),lFactor*stPoint->GetY(),lFactor*stPoint->GetZ(),
			    lFactor*stPoint->GetPrecisionX(),lFactor*stPoint->GetPrecisionY(),lFactor*stPoint->GetPrecisionZ(),
			    stPoint->GetType(), stPoint->GetTarget());
}
