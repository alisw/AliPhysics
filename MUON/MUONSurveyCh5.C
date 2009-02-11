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

/* $Id:$ */

/// \ingroup macros
/// \file MUONSurveyCh8L.C
/// \brief Macro to process survey and photogrammetry data of chamber 5
/// 
/// Macro loads the survey data from .txt file using AliSurveyObj.
/// Macro uses AliMUONSurvey... classes.
/// The transformations of the detection elements or chambers can be obtained 
/// in three ways:
/// - A: Fit of local to global transformation using the fixed button targets.
/// - B: Fit a plane to the sticker targets -> psi, theta
///      Use above psi and theta and fit remaining 4 parameters using the fixed 
///      button targets
/// - C: Fit a plane to the sticker targets -> psi, theta
///      Use above psi and theta to calculate xc, yc, zc and phi by solving 
///      the equations from a local to global transformation of the fixed button 
///      targets
///
/// Methods A and B are prefered to C, and B is better if sticker targets are 
/// available and lie on a plane!
/// For slats only methods B and C can be used.
/// Various histograms are filled and printed for monitoring.
/// MisAlignment object is then created.
///
/// \author Javier Castillo

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryMisAligner.h"
#include "AliMUONSurveyChamber.h"
#include "AliMUONSurveyDetElem.h"
#include "AliMUONSurveyUtil.h"

#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliAlignObjMatrix.h"
#include "AliAlignObj.h"

#include <TROOT.h>
#include <TGeoManager.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TArrayD.h>
#include <TMath.h>
#include <TString.h>
#include <TFitter.h>
#include <TH2.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TPostScript.h>
#include <TPaveLabel.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMatrixDSym.h>

#include <fstream>

#endif

void MUONSurveyCh5() {
  
  TString str;
  TString title;
  
  Bool_t bMonitor = kTRUE;
  Bool_t bOCDB = kTRUE;
  Bool_t saveps = kFALSE;
  const int cWidth = (int)(700*(29./21.));
  const int cHeight = 700;
  const int filetype  = 112; // landscape  

  Int_t chId = 4;
  Int_t nChs = 2;
  Int_t nDetElems = 18;
  Int_t nDetElemsI = 9;
  //  Int_t nDetElemsO = 9;
  Int_t iDetElemToDetElemId[18] = {514,506,516,508,500,510,502,512,504,513,503,511,501,509,517,507,515,505};
  Int_t iDetElemPseudoIdToDetElem[18] = {4,12,6,10,8,17,1,15,3,13,5,11,7,9,0,16,2,14};
  Int_t iDetElemsOfChamber[2][9] = {{0,16,2,14,4,12,6,10,8},
				    {9,7,11,5,13,3,15,1,17}};
  
  TObjArray *myChamberArray = new TObjArray(); 
  myChamberArray->Add(new AliMUONSurveyChamber(chId));  
  myChamberArray->Add(new AliMUONSurveyChamber(chId));  

  AliMUONSurveyChamber *myChamberI = 0x0;
  AliMUONSurveyChamber *myChamberO = 0x0;
  AliMUONSurveyChamber *myChamber = 0x0;
  AliMUONSurveyDetElem *myDetElem = 0x0;

  myChamberI = (AliMUONSurveyChamber*)myChamberArray->At(0);
  myChamberI->GetSurveyObj()->FillFromLocalFile("../Reports/AliceSt3_TC5_1381b.txt");
  myChamberO = (AliMUONSurveyChamber*)myChamberArray->At(1);
  myChamberO->GetSurveyObj()->FillFromLocalFile("../Reports/AliceSt3_TC5_1381b.txt");

  myChamber = myChamberI;
  myChamber->PrintSurveyReport();

  // Chamber & DetElems button targets local coordinates
  AliSurveyObj *lSO = new AliSurveyObj();    
  lSO->FillFromLocalFile("$ALICE_ROOT/MUON/data/MUONTrackingTargetsLocal.txt");

  // Set survey targets for chambers 
  myChamberI->AddGButtonTargets(Form("%dC_1000",chId+1),4);
  myChamberI->AddStickerTargets(Form("%dA_702",chId+1),9);
  myChamberI->AddStickerTargets(Form("%dA_703",chId+1),9);
  myChamberI->AddStickerTargets(Form("%dC_700",chId+1),9);
  myChamberI->AddStickerTargets(Form("%dC_701",chId+1),9);
  myChamberI->AddLButtonTargets(lSO->GetData(),Form("%d_1000",chId+1),4);
  myChamberO->AddGButtonTargets(Form("%dC_1010",chId+1),4);
  myChamberO->AddStickerTargets(Form("%dA_700",chId+1),9);
  myChamberO->AddStickerTargets(Form("%dA_701",chId+1),9);
  myChamberO->AddStickerTargets(Form("%dC_702",chId+1),9);
  myChamberO->AddStickerTargets(Form("%dC_703",chId+1),9);
  myChamberO->AddLButtonTargets(lSO->GetData(),Form("%d_1010",chId+1),4);

  // Set survey targets for detection elements 
  for (int iCh =0; iCh<=1; iCh++) {
    myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);
    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      Int_t iDetElem = iDetElemsOfChamber[iCh][iDetElemI];
      myChamber->AddSurveyDetElem(iDetElemToDetElemId[iDetElem]);
      TString baseName;
      myDetElem =  myChamber->GetDetElem(iDetElemI);
      if (myDetElem) {
	if(iDetElem<9) baseName = Form("%dA_50%d",chId+1,iDetElem+1);	
	else baseName = Form("%dC_50%d",chId+1,iDetElem+1-9);
	myDetElem->AddStickerTargets(baseName,9);
	if(iDetElem<9) baseName = Form("%dA_60%d",chId+1,iDetElem+1);
	else baseName = Form("%dC_60%d",chId+1,iDetElem+1-9);
	myDetElem->AddGButtonTargets(baseName,2);
	if(iDetElemToDetElemId[iDetElem]-500<10) baseName = Form("%d_50%d",chId+1,iDetElemToDetElemId[iDetElem]-500);
	else baseName = Form("%d_5%d",chId+1,iDetElemToDetElemId[iDetElem]-500);
	myDetElem->AddLButtonTargets(lSO->GetData(),baseName,2);
      }
    }
  }
  printf("All targets added! \n");

  Double_t **lCenDetElem = new Double_t*[nDetElems+nChs];
  Double_t **lRotDetElem = new Double_t*[nDetElems+nChs];
  Double_t **lDiffCenDetElem0 = new Double_t*[nDetElems+nChs];
  Double_t **lDiffRotDetElem0 = new Double_t*[nDetElems+nChs];
  Double_t **lDiffThCenDetElem0 = new Double_t*[nDetElems+nChs];
  Double_t **lDiffThRotDetElem0 = new Double_t*[nDetElems+nChs];
  Double_t **lDeltaDiffCenDetElem0 = new Double_t*[nDetElems+nChs];
  Double_t **lDeltaDiffRotDetElem0 = new Double_t*[nDetElems+nChs];

  for (int iDetElem=0; iDetElem<nDetElems+nChs; iDetElem++) {
    lCenDetElem[iDetElem] = new Double_t[3];
    lRotDetElem[iDetElem] = new Double_t[3];
    lDiffCenDetElem0[iDetElem] = new Double_t[3];
    lDiffRotDetElem0[iDetElem] = new Double_t[3];
    lDiffThCenDetElem0[iDetElem] = new Double_t[3];
    lDiffThRotDetElem0[iDetElem] = new Double_t[3];
    lDeltaDiffCenDetElem0[iDetElem] = new Double_t[3];
    lDeltaDiffRotDetElem0[iDetElem] = new Double_t[3];
  }

  TGeoCombiTrans dtrfDetElem[nDetElems+nChs];
  TGeoCombiTrans localTrfDetElem[nDetElems+nChs];
  TGeoCombiTrans localTrfThDetElem[nDetElems+nChs];

  // Import TGeo geometry 
  char* geoFilename = "geometry.root";
  if ( ! AliGeomManager::GetGeometry() ) {
    AliGeomManager::LoadGeometry(geoFilename);
    if (! AliGeomManager::GetGeometry() ) {
      printf("MUONSurveyCh%d: getting geometry from file %s failed\n", chId+1,geoFilename);
      return;
    }
    cout << "geometry imported" << endl;
  }

  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer();
  transform->LoadGeometryData();

  TGeoCombiTrans trfThChamber;
  TGeoCombiTrans trfThDetElem;

  for (int iCh =0; iCh<=1; iCh++) {
    myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);

    trfThChamber = TGeoCombiTrans(*transform->GetModuleTransformerByDEId(iDetElemToDetElemId[iCh*nDetElemsI])->GetTransformation());
    trfThChamber.Print();
    myChamber->SetLocalTransformation(new TGeoCombiTrans(trfThChamber),kTRUE);

    // Set Chamber plane function
    cout << "Setting plane for Chamber" << iCh+1 << " ..." << endl;
    myChamber->SetPlane(Form("fChamber%d",iCh+1));
    myChamber->SetPlaneParameters(0.,0.,0.);

    // Fit a plane to sticker targets
    Double_t lCChi2 = myChamber->FitPlane();
    cout << "... chi2 = " << lCChi2 << " ..." << endl; 
    
    // Get best transformation from fitting method 
    // (local to global transformation)
    cout << "Trying fitting method for chamber " << iCh << endl;
    myChamber->SurveyToAlign();    
    //    myChamber->SurveyToAlign(myChamber->GetPlane()->GetParameter(0),myChamber->GetPlane()->GetParameter(1),myChamber->GetPlane()->GetParError(0),myChamber->GetPlane()->GetParError(1));    
    
    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      myDetElem =  myChamber->GetDetElem(iDetElemI);
      Int_t iDetElem = iDetElemsOfChamber[iCh][iDetElemI];

      trfThDetElem.Clear();
      trfThDetElem = TGeoCombiTrans(*transform->GetDetElement(iDetElemToDetElemId[iDetElem])->GetLocalTransformation());
      trfThDetElem.Print();
      myDetElem->SetLocalTransformation(new TGeoCombiTrans(trfThDetElem),kTRUE);

      for (int iCor=0; iCor<3; iCor++){
	lCenDetElem[iDetElem][iCor] = 0;
	lRotDetElem[iDetElem][iCor] = 0;
      }

      if (bMonitor){
	// MONITOR: Draw graph with sticker targets for plane fitting
	myDetElem->DrawSTargets();
	gPad->Update();
      }

      // Get DetElem transformation. 
      // Psi and Theta are obtained by fitting a plane to the sticker targets.
      // Then Xc, Yc, Zc and Phi are obtained by solving the equations to the ref.
      // syst. transformation of the button targets

      // Set DetElem plane function
      cout << "Setting plane for DetElem" << iDetElem+1 << " ..." << endl;
      myDetElem->SetPlane(Form("fDetElem%d",iDetElem+1));
      myDetElem->SetPlaneParameters(0.,0.,3.)
;
      // Fit a plane to sticker targets
      Double_t lChi2 = myDetElem->FitPlane();
      cout << "... chi2 = " << lChi2 << " ..." << endl; 
      
      lRotDetElem[iDetElem][0]=TMath::ATan(myDetElem->GetPlane()->GetParameter(0));
      lRotDetElem[iDetElem][1]=TMath::ATan(myDetElem->GetPlane()->GetParameter(1));
      
      // Calculate Mean transformation using previous plane fit 
      // and pairs of button targets
      myDetElem->CalculateMeanTransf(lCenDetElem[iDetElem],lRotDetElem[iDetElem]);
    
      cout << "DetElem" << iDetElem+1 << " : mycen(" << lCenDetElem[iDetElem][0] << "," << lCenDetElem[iDetElem][1] << "," << lCenDetElem[iDetElem][2] << "); rot(" << lRotDetElem[iDetElem][0] << "," << lRotDetElem[iDetElem][1] << "," << lRotDetElem[iDetElem][2] << ")"  << endl;  	
      
     
      // Get best transformation from fitting method 
      // (local to global transformation)
      cout << "Trying fitting method for DetElemId " << iDetElemToDetElemId[iDetElem] << endl;
      //      myDetElem->SurveyToAlign();     
      myDetElem->SurveyToAlign(lRotDetElem[iDetElem][0],lRotDetElem[iDetElem][1],myDetElem->GetPlane()->GetParError(0),myDetElem->GetPlane()->GetParError(1));    
    }
  }

  // Print found transformation of Detection Element (plane fit + loop)        
  for (int iCh =0; iCh<=1; iCh++) {
    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      Int_t iDetElem = iDetElemsOfChamber[iCh][iDetElemI];
      cout << "DetElem" << iDetElem+1 << " : mycen(" << lCenDetElem[iDetElem][0] << "," << lCenDetElem[iDetElem][1] << "," << lCenDetElem[iDetElem][2] << "); rot(" << lRotDetElem[iDetElem][0] << "," << lRotDetElem[iDetElem][1] << "," << lRotDetElem[iDetElem][2] << ")"  << endl;  	
    }
  }

  // Print Theoretical transformation of Detection Element
  for (int iCh =0; iCh<=1; iCh++) {
    myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);
    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      myChamber->GetDetElem(iDetElemI)->PrintLocalTrf();
    }
  }

  // Print found delta transformation of Detection Element
  for (int iCh =0; iCh<=1; iCh++) {
    myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);      
    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      myChamber->GetDetElem(iDetElemI)->PrintAlignTrf();
    }
  }

  //
  // Compare transformations to expected ones 
  //
  Int_t iDetElemToPos[18] = {0, 16, 2, 14, 4, 12, 6, 10, 8, 9, 7, 11, 5, 13, 3, 15, 1, 17};

  TGraphErrors *gDeltaDiffCenXDetElem0 = new TGraphErrors(nDetElems);
  TGraphErrors *gDeltaDiffCenYDetElem0 = new TGraphErrors(nDetElems);
  TGraphErrors *gDeltaDiffCenZDetElem0 = new TGraphErrors(nDetElems);
  TGraphErrors *gDeltaDiffPsiDetElem0 = new TGraphErrors(nDetElems);
  TGraphErrors *gDeltaDiffThtDetElem0 = new TGraphErrors(nDetElems);
  TGraphErrors *gDeltaDiffPhiDetElem0 = new TGraphErrors(nDetElems);

  for (int iCh =0; iCh<=1; iCh++) {
    myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);      
    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      myChamber->GetDetElem(iDetElemI)->GetAlignTrf()->Print();
    }
  }

  for (int iCh =0; iCh<=1; iCh++) {
    myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);
    // Store delta transformations to use for CDB entry creation
    dtrfDetElem[nDetElems+iCh].Clear();
    dtrfDetElem[nDetElems+iCh] = *(myChamber->GetAlignTrf());

    // Those are for checks and visualizations
    localTrfDetElem[nDetElems+iCh].Clear();
    localTrfDetElem[nDetElems+iCh] = (*(myChamber->GetLocalTrf())*(*(myChamber->GetAlignTrf())));
    localTrfDetElem[nDetElems+iCh].Print();
    lDiffCenDetElem0[nDetElems+iCh] = (Double_t*)localTrfDetElem[nDetElems+iCh].GetTranslation();
    AliMUONSurveyUtil::MatrixToAngles(localTrfDetElem[nDetElems+iCh].GetRotationMatrix(),lDiffRotDetElem0[nDetElems+iCh]);
    
    localTrfThDetElem[nDetElems+iCh].Clear();
    localTrfThDetElem[nDetElems+iCh] = (*(myChamber->GetLocalTrf()));
    localTrfThDetElem[nDetElems+iCh].Print();
    lDiffThCenDetElem0[nDetElems+iCh] = (Double_t*)localTrfThDetElem[nDetElems+iCh].GetTranslation();
    AliMUONSurveyUtil::MatrixToAngles(localTrfThDetElem[nDetElems+iCh].GetRotationMatrix(),lDiffThRotDetElem0[nDetElems+iCh]);

    for (int iCor=0; iCor<3; iCor++){
      lDeltaDiffCenDetElem0[nDetElems+iCh][iCor] = lDiffCenDetElem0[nDetElems+iCh][iCor]-lDiffThCenDetElem0[nDetElems+iCh][iCor];
      lDeltaDiffRotDetElem0[nDetElems+iCh][iCor] = lDiffRotDetElem0[nDetElems+iCh][iCor]-lDiffThRotDetElem0[nDetElems+iCh][iCor];
      if (lDeltaDiffRotDetElem0[nDetElems+iCh][iCor]>TMath::Pi()) lDeltaDiffRotDetElem0[nDetElems+iCh][iCor]-=TMath::TwoPi();
      if (lDeltaDiffRotDetElem0[nDetElems+iCh][iCor]<-TMath::Pi()) lDeltaDiffRotDetElem0[nDetElems+iCh][iCor]+=TMath::TwoPi();
    }      

    for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
      Int_t iDetElem = iDetElemsOfChamber[iCh][iDetElemI];
      myDetElem =  myChamber->GetDetElem(iDetElemI);
      // Store delta transformations to use for CDB entry creation
      dtrfDetElem[iDetElem].Clear();
      dtrfDetElem[iDetElem] = *(myDetElem->GetAlignTrf());

      // Those are for checks and visualizations
      localTrfDetElem[iDetElem].Clear();
      localTrfDetElem[iDetElem] = (*(myDetElem->GetLocalTrf())*(*(myDetElem->GetAlignTrf())));
      //      localTrfDetElem[iDetElem] = (*(myDetElem->GetBaseTrf())*(*(myDetElem->GetAlignTrf())));
      localTrfDetElem[iDetElem].Print();
      lDiffCenDetElem0[iDetElem] = (Double_t*)localTrfDetElem[iDetElem].GetTranslation();
      AliMUONSurveyUtil::MatrixToAngles(localTrfDetElem[iDetElem].GetRotationMatrix(),lDiffRotDetElem0[iDetElem]);
//       lDiffCenDetElem0[iDetElem] = lCenDetElem[iDetElem];
//       lDiffRotDetElem0[iDetElem] = lRotDetElem[iDetElem];

      localTrfThDetElem[iDetElem].Clear();
      localTrfThDetElem[iDetElem] = (*(myDetElem->GetLocalTrf()));
      //      localTrfThDetElem[iDetElem] = (*(myDetElem->GetBaseTrf()));
      localTrfThDetElem[iDetElem].Print();
      lDiffThCenDetElem0[iDetElem] = (Double_t*)localTrfThDetElem[iDetElem].GetTranslation();
      AliMUONSurveyUtil::MatrixToAngles(localTrfThDetElem[iDetElem].GetRotationMatrix(),lDiffThRotDetElem0[iDetElem]);
      
      for (int iCor=0; iCor<3; iCor++){
	lDeltaDiffCenDetElem0[iDetElem][iCor] = lDiffCenDetElem0[iDetElem][iCor]-lDiffThCenDetElem0[iDetElem][iCor];
	lDeltaDiffRotDetElem0[iDetElem][iCor] = lDiffRotDetElem0[iDetElem][iCor]-lDiffThRotDetElem0[iDetElem][iCor];
	if (lDeltaDiffRotDetElem0[iDetElem][iCor]>TMath::Pi()) lDeltaDiffRotDetElem0[iDetElem][iCor]-=TMath::TwoPi();
	if (lDeltaDiffRotDetElem0[iDetElem][iCor]<-TMath::Pi()) lDeltaDiffRotDetElem0[iDetElem][iCor]+=TMath::TwoPi();
      }

      gDeltaDiffCenXDetElem0->SetPoint(iDetElem,1e1*lDeltaDiffCenDetElem0[iDetElem][0],iDetElemToPos[iDetElem]+1);
      gDeltaDiffCenYDetElem0->SetPoint(iDetElem,1e1*lDeltaDiffCenDetElem0[iDetElem][1],iDetElemToPos[iDetElem]+1);
      gDeltaDiffCenZDetElem0->SetPoint(iDetElem,1e1*lDeltaDiffCenDetElem0[iDetElem][2],iDetElemToPos[iDetElem]+1);
      gDeltaDiffPsiDetElem0->SetPoint(iDetElem,1e3*lDeltaDiffRotDetElem0[iDetElem][0],iDetElemToPos[iDetElem]+1);
      gDeltaDiffThtDetElem0->SetPoint(iDetElem,1e3*lDeltaDiffRotDetElem0[iDetElem][1],iDetElemToPos[iDetElem]+1);
      gDeltaDiffPhiDetElem0->SetPoint(iDetElem,1e3*lDeltaDiffRotDetElem0[iDetElem][2],iDetElemToPos[iDetElem]+1);
      gDeltaDiffCenXDetElem0->SetPointError(iDetElem,1e1*myDetElem->GetFitter()->GetParError(0),0.);
      gDeltaDiffCenYDetElem0->SetPointError(iDetElem,1e1*myDetElem->GetFitter()->GetParError(1),0.);
      gDeltaDiffCenZDetElem0->SetPointError(iDetElem,1e1*myDetElem->GetFitter()->GetParError(2),0.);
      gDeltaDiffPsiDetElem0->SetPointError(iDetElem,1e3*myDetElem->GetFitter()->GetParError(3),0.);
      gDeltaDiffThtDetElem0->SetPointError(iDetElem,1e3*myDetElem->GetFitter()->GetParError(4),0.);
      gDeltaDiffPhiDetElem0->SetPointError(iDetElem,1e3*myDetElem->GetFitter()->GetParError(5),0.);  
    }
  }

  // Apply the found misalignments to the geometry
  AliMUONGeometryTransformer *newTransform = AliMUONSurveyUtil::ReAlign(transform,chId,nDetElems,iDetElemPseudoIdToDetElem,dtrfDetElem,true); 
  newTransform->WriteTransformations(Form("transform2ReAlignSurveyCh%d.dat",chId+1));

  if(bOCDB){
    // Generate realigned data in local cdb
    const TClonesArray* array = newTransform->GetMisAlignmentData();
    
    // Set the alignment resolution in the align objects for this chamber
    Double_t chResX = myChamberI->GetAlignResX();
    Double_t chResY = myChamberI->GetAlignResY();
    Double_t deResX = (myChamberI->GetMeanDetElemAlignResX()+myChamberO->GetMeanDetElemAlignResX())/2.;
    Double_t deResY = (myChamberI->GetMeanDetElemAlignResY()+myChamberO->GetMeanDetElemAlignResY())/2.;
    printf("Chamber alignment resolutions: resX=%f , resY=%f\n",chResX,chResY); 
    printf("Detection Elements alignment resolutions: resX=%f , resY=%f\n",deResX,deResY); 
    chResX = TMath::Sqrt(0.1*0.1+chResX*chResX);
    chResY = TMath::Sqrt(0.1*0.1+chResY*chResY);
    AliMUONSurveyUtil::SetAlignmentResolution(array,chId,chResX,chResY,deResX,deResY);
    
    // CDB manager
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    cdbManager->SetDefaultStorage(Form("local://ReAlignSurveyCh%dCDB",chId+1));
    
    AliCDBMetaData* cdbData = new AliCDBMetaData();
    cdbData->SetResponsible("Dimuon Offline project");
    cdbData->SetComment("MUON alignment objects with survey realignment");
    AliCDBId id("MUON/Align/Data", 0, AliCDBRunRange::Infinity()); 
    cdbManager->Put(const_cast<TClonesArray*>(array), id, cdbData);
  }

  for(Int_t iCor=0; iCor<3; iCor++){
    for(Int_t iDetElem=0; iDetElem<nDetElems; iDetElem++){
      cout << lDeltaDiffCenDetElem0[iDetElem][iCor] << " " << lDiffCenDetElem0[iDetElem][iCor] << " " << lDiffThCenDetElem0[iDetElem][iCor] << endl;
    }
    cout << endl;
  }
  cout << endl;
  for(Int_t iCor=0; iCor<3; iCor++){
    for(Int_t iDetElem=0; iDetElem<nDetElems; iDetElem++){
      cout << lDeltaDiffRotDetElem0[iDetElem][iCor] << " " << lDiffRotDetElem0[iDetElem][iCor] << " " << lDiffThRotDetElem0[iDetElem][iCor] << endl;
    }
    cout << endl;
  }

  TH1F *myDetElemDeltaDiffCenX = new TH1F("myDetElemDeltaDiffCenX","myDetElemDeltaDiffCenX",100,-10,10);
  myDetElemDeltaDiffCenX->SetMaximum(nDetElems+1);
  myDetElemDeltaDiffCenX->SetMinimum(0);
  TH1F *myDetElemDeltaDiffCenY = new TH1F("myDetElemDeltaDiffCenY","myDetElemDeltaDiffCenY",100,-10,10);
  myDetElemDeltaDiffCenY->SetMaximum(nDetElems+1);
  myDetElemDeltaDiffCenY->SetMinimum(0);
  TH1F *myDetElemDeltaDiffCenZ = new TH1F("myDetElemDeltaDiffCenZ","myDetElemDeltaDiffCenZ",100,-30,30);
  myDetElemDeltaDiffCenZ->SetMaximum(nDetElems+1);
  myDetElemDeltaDiffCenZ->SetMinimum(0);

  TH1F *myDetElemDeltaDiffRotX = new TH1F("myDetElemDeltaDiffRotX","myDetElemDeltaDiffRotX",100,-15,15);
  myDetElemDeltaDiffRotX->SetMaximum(nDetElems+1);
  myDetElemDeltaDiffRotX->SetMinimum(0);
  TH1F *myDetElemDeltaDiffRotY = new TH1F("myDetElemDeltaDiffRotY","myDetElemDeltaDiffRotY",100,-15,15);
  myDetElemDeltaDiffRotY->SetMaximum(nDetElems+1);
  myDetElemDeltaDiffRotY->SetMinimum(0);
  TH1F *myDetElemDeltaDiffRotZ = new TH1F("myDetElemDeltaDiffRotZ","myDetElemDeltaDiffRotZ",100,-5,5);
  myDetElemDeltaDiffRotZ->SetMaximum(nDetElems+1);
  myDetElemDeltaDiffRotZ->SetMinimum(0);

  //
  // ******** Starting plots 
  //
  TCanvas *canvas;
  TPad *pad;
  TPaveLabel *theTitle;
  gStyle->SetPalette(1);

  TPostScript *ps = 0;
  
  if( saveps ){
    ps = new TPostScript(Form("surveyChamber%d",chId+1),filetype); 
    ps->NewPage();
  }

  // Observed misalignments
  str = Form("Chamber %d",chId+1);
  TCanvas *cvn2 = new TCanvas("cvn2",str,cWidth,cHeight);
  canvas = cvn2;
  canvas->Range(0,0,21,29);
  
  title = Form(" MisAlignments Chamber %d - PL2G - In Frame ",chId+1);
  TPaveLabel *theTitle2 = new TPaveLabel(3,27.0,18,28.5,title,"br");
  theTitle = theTitle2;
  theTitle->SetFillColor(18);
  theTitle->SetTextFont(32);
  theTitle->SetTextSize(0.4);
  theTitle->SetTextColor(1);
  theTitle->Draw();
 
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.98,0.91,0);
  pad = pad2;
  pad->Draw();
  TLine *ch0Line = new TLine(0,1,0,2);
  TLine *ch1Line = new TLine(0,1,0,2);
  ch1Line->SetLineStyle(2);
  pad->Divide(3,2);
  
  pad->cd(1);
  myDetElemDeltaDiffCenX->Draw();
  myDetElemDeltaDiffCenX->SetXTitle("#Delta[xc_{i}^{m}-xc_{i}^{th}] (mm)");
  myDetElemDeltaDiffCenX->SetYTitle("DetElem arbitrary ordering");
  gDeltaDiffCenXDetElem0->SetMarkerStyle(20);
  gDeltaDiffCenXDetElem0->Draw("P");
  ch0Line->DrawLine(1e1*lDeltaDiffCenDetElem0[nDetElems+0][0],0.5,1e1*lDeltaDiffCenDetElem0[nDetElems+0][0],9.5);
  ch1Line->DrawLine(1e1*lDeltaDiffCenDetElem0[nDetElems+1][0],9.5,1e1*lDeltaDiffCenDetElem0[nDetElems+1][0],18.5);

  pad->cd(2);
  myDetElemDeltaDiffCenY->Draw();
  myDetElemDeltaDiffCenY->SetXTitle("#Delta[yc_{i}^{m}-yc_{i}^{th}] (mm)");
  myDetElemDeltaDiffCenY->SetYTitle("DetElem arbitrary ordering");
  gDeltaDiffCenYDetElem0->SetMarkerStyle(20);
  gDeltaDiffCenYDetElem0->Draw("P");
  ch0Line->DrawLine(1e1*lDeltaDiffCenDetElem0[nDetElems+0][1],0.5,1e1*lDeltaDiffCenDetElem0[nDetElems+0][1],9.5);
  ch1Line->DrawLine(1e1*lDeltaDiffCenDetElem0[nDetElems+1][1],9.5,1e1*lDeltaDiffCenDetElem0[nDetElems+1][1],18.5);

  pad->cd(3);
  myDetElemDeltaDiffCenZ->Draw();
  myDetElemDeltaDiffCenZ->SetXTitle("#Delta[zc_{i}^{m}-zc_{i}^{th}] (mm)");
  myDetElemDeltaDiffCenZ->SetYTitle("DetElem arbitrary ordering");
  gDeltaDiffCenZDetElem0->SetMarkerStyle(20);
  gDeltaDiffCenZDetElem0->Draw("P");
  ch0Line->DrawLine(1e1*lDeltaDiffCenDetElem0[nDetElems+0][2],0.5,1e1*lDeltaDiffCenDetElem0[nDetElems+0][2],9.5);
  ch1Line->DrawLine(1e1*lDeltaDiffCenDetElem0[nDetElems+1][2],9.5,1e1*lDeltaDiffCenDetElem0[nDetElems+1][2],18.5);

  pad->cd(4);
  myDetElemDeltaDiffRotX->Draw();
  myDetElemDeltaDiffRotX->SetXTitle("#Delta[#psi_{i}^{m}-#psi_{i}^{th}] (mrad)");
  myDetElemDeltaDiffRotX->SetYTitle("DetElem arbitrary ordering");
  gDeltaDiffPsiDetElem0->SetMarkerStyle(20);
  gDeltaDiffPsiDetElem0->Draw("P");
  ch0Line->DrawLine(1e3*lDeltaDiffRotDetElem0[nDetElems+0][0],0.5,1e3*lDeltaDiffRotDetElem0[nDetElems+0][0],9.5);
  ch1Line->DrawLine(1e3*lDeltaDiffRotDetElem0[nDetElems+1][0],9.5,1e3*lDeltaDiffRotDetElem0[nDetElems+1][0],18.5);

  pad->cd(5);
  myDetElemDeltaDiffRotY->Draw();
  myDetElemDeltaDiffRotY->SetXTitle("#Delta[#theta_{i}^{m}-#theta_{i}^{th}] (mrad)");
  myDetElemDeltaDiffRotY->SetYTitle("DetElem arbitrary ordering");
  gDeltaDiffThtDetElem0->SetMarkerStyle(20);
  gDeltaDiffThtDetElem0->Draw("P");
  ch0Line->DrawLine(1e3*lDeltaDiffRotDetElem0[nDetElems+0][1],0.5,1e3*lDeltaDiffRotDetElem0[nDetElems+0][1],9.5);
  ch1Line->DrawLine(1e3*lDeltaDiffRotDetElem0[nDetElems+1][1],9.5,1e3*lDeltaDiffRotDetElem0[nDetElems+1][1],18.5);

  pad->cd(6);
  myDetElemDeltaDiffRotZ->Draw();
  myDetElemDeltaDiffRotZ->SetXTitle("#Delta[#phi_{i}^{m}-#phi_{i}^{th}] (mrad)");
  myDetElemDeltaDiffRotZ->SetYTitle("DetElem arbitrary ordering");
  gDeltaDiffPhiDetElem0->SetMarkerStyle(20);
  gDeltaDiffPhiDetElem0->Draw("P");
  ch0Line->DrawLine(1e3*lDeltaDiffRotDetElem0[nDetElems+0][2],0.5,1e3*lDeltaDiffRotDetElem0[nDetElems+0][2],9.5);
  ch1Line->DrawLine(1e3*lDeltaDiffRotDetElem0[nDetElems+1][2],9.5,1e3*lDeltaDiffRotDetElem0[nDetElems+1][2],18.5);

  pad->Update();
 
  if(bMonitor){    
    // MONITOR: Histograms for monitoring
    TH2F *hCPSTa = new TH2F("hCPSTa","hCPSTa",100,-200,200,100,-200,200);
    TH2F *hCPSTc = new TH2F("hCPSTc","hCPSTc",100,-200,200,100,-200,200);
    TH2F *hSSTa = new TH2F("hSSTa","hSSTa",200,-200,200,200,-200,200);
    TH2F *hSSTc = new TH2F("hSSTc","hSSTc",200,-200,200,200,-200,200);
    TH2F *hSSTap = new TH2F("hSSTap","hSSTap",800,-200,200,800,-200,200);
    TH2F *hSSTcp = new TH2F("hSSTcp","hSSTcp",800,-200,200,800,-200,200);
    
    // MONITOR: Fill histograms with chambers and slats sticker target positions
    for (int iCh =0; iCh<=1; iCh++) {
      myChamber = (AliMUONSurveyChamber*)myChamberArray->At(iCh);            
      cout << "Filling histograms of sticker target for chamber" << iCh+1 << " ..." << endl;
      myChamber->FillCPSTHistograms(TString("C"),hCPSTc,TString("A"),hCPSTa);
      myChamber->FillDESTHistograms(TString("C"),hSSTc,TString("A"),hSSTa);

      for (int iDetElemI=0; iDetElemI<nDetElemsI; iDetElemI++){
	//	Int_t iDetElem = iDetElemsOfChamber[iCh][iDetElemI];
	myDetElem =  myChamber->GetDetElem(iDetElemI);
	// MONITOR: Fill slat plane for fit monitor.
	Double_t pl[3] = {0};
	Double_t pg[3] = {0};
	AliSurveyPoint *pointSBT0 = myDetElem->GetLButtonTarget(0);
	AliSurveyPoint *pointSBT1 = myDetElem->GetLButtonTarget(1);
	if(pointSBT0&&pointSBT1) {
	  if (pointSBT0->GetX()>pointSBT1->GetX()){
	    pointSBT0=myDetElem->GetLButtonTarget(1);
	    pointSBT1=myDetElem->GetLButtonTarget(0);
	  }
	  Double_t lX = pointSBT0->GetX();
	  while(lX<pointSBT1->GetX()){
	    Double_t lY = pointSBT0->GetY()-20;
	    while(lY<pointSBT0->GetY()+20){
	      pl[0] = lX;  pl[1] = lY;  pl[2] = 0.;
	      (TGeoCombiTrans((*(myDetElem->GetBaseTrf()))*(*(myDetElem->GetAlignTrf())))).LocalToMaster(pl,pg);
	      if(myDetElem->GetGButtonTarget(0)->GetPointName().Contains("A")){
		if (hSSTap->GetBinContent(hSSTap->FindBin(pg[0],pg[1]))==0)
		  hSSTap->Fill(pg[0],pg[1],-pg[2]);
	      }
	      else {
		if (hSSTcp->GetBinContent(hSSTcp->FindBin(pg[0],pg[1]))==0)
		  hSSTcp->Fill(pg[0],pg[1],-pg[2]);
	      }
	      lY+=hSSTap->GetYaxis()->GetBinWidth(1);
	    }
	    lX+=hSSTap->GetXaxis()->GetBinWidth(1);
	  }
	}
      }
    }

    if( saveps ){
      ps->NewPage();
    }

    // View from side A
    str = Form("Chamber %d - side A",chId+1);
    TCanvas *cvn0 = new TCanvas("cvn0",str,cWidth,cHeight);
    canvas = cvn0;
    canvas->Range(0,0,21,29);
  
    title = Form(" Deformations of chamber %d - View from side A ",chId+1);
    TPaveLabel *theTitle0 = new TPaveLabel(3,27.0,18,28.5,title,"br");
    theTitle = theTitle0;
    theTitle->SetFillColor(18);
    theTitle->SetTextFont(32);
    theTitle->SetTextSize(0.4);
    theTitle->SetTextColor(1);
    theTitle->Draw();
 
    TPad *pad0 = new TPad("pad0","pad0",0.01,0.01,0.98,0.91,0);
    pad = pad0;
    pad->Draw();
    pad->Divide(2,2);

    Double_t lMin, lMax;

    pad->cd(1);
    lMin = hCPSTa->GetMinimum(0);
    hCPSTa->SetMinimum(TMath::Floor(lMin));
    lMax = hCPSTa->GetMaximum();
    hCPSTa->SetMaximum(TMath::Ceil(lMax));
    hCPSTa->Draw("lego2z");
    
    pad->cd(2);
    lMin = hSSTa->GetMinimum(0);
    hSSTa->SetMinimum(TMath::Floor(lMin));
    lMax = hSSTa->GetMaximum();
    hSSTa->SetMaximum(TMath::Ceil(lMax));
    hSSTa->Draw("lego2z");
    
    pad->cd(3);
    
    pad->cd(4);
    lMin = hSSTap->GetMinimum(0);
    hSSTap->SetMinimum(TMath::Floor(lMin));
    lMax = hSSTap->GetMaximum();
    hSSTap->SetMaximum(TMath::Ceil(lMax));
    hSSTap->Draw("surf2z");
    
    pad->Update();
    if(saveps){
      ps->NewPage();
    }

    // Inv Mass, Multiplicity
    str = Form("Chamber %d - side C",chId+1);
    TCanvas *cvn1 = new TCanvas("cvn1",str,cWidth,cHeight);
    canvas = cvn1;
    canvas->Range(0,0,21,29);

    title = Form(" Deformations of chamber %d - View from side C ",chId+1);
    TPaveLabel *theTitle1 = new TPaveLabel(3,27.0,18,28.5,title,"br");
    theTitle = theTitle1;
    theTitle->SetFillColor(18);
    theTitle->SetTextFont(32);
    theTitle->SetTextSize(0.4);
    theTitle->SetTextColor(1);
    theTitle->Draw();
 
    TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.98,0.91,0);
    pad = pad1;
    pad->Draw();
    pad->Divide(2,2);

    pad->cd(1);
    lMin = hCPSTc->GetMinimum(0);
    hCPSTc->SetMinimum(TMath::Floor(lMin));
    lMax = hCPSTc->GetMaximum();
    hCPSTc->SetMaximum(TMath::Ceil(lMax));
    hCPSTc->Draw("lego2z");

    pad->cd(2);
    lMin = hSSTc->GetMinimum(0);
    hSSTc->SetMinimum(TMath::Floor(lMin));
    lMax = hSSTc->GetMaximum();
    hSSTc->SetMaximum(TMath::Ceil(lMax));
    hSSTc->Draw("lego2z");

    pad->cd(3);

    pad->cd(4);
    lMin = hSSTcp->GetMinimum(0);
    hSSTcp->SetMinimum(TMath::Floor(lMin));
    lMax = hSSTcp->GetMaximum();
    hSSTcp->SetMaximum(TMath::Ceil(lMax));
    hSSTcp->Draw("surf2z");
  }

  pad->Update();
  if( saveps ){
    ps->Close();
  }

  TFile *hFile = new TFile(Form("spCH%d_PL2G_IF.root",chId+1),"RECREATE");
  gDeltaDiffCenXDetElem0->Write("gDeltaDiffCenXDetElem0");
  gDeltaDiffCenYDetElem0->Write("gDeltaDiffCenYDetElem0");
  gDeltaDiffCenZDetElem0->Write("gDeltaDiffCenZDetElem0");
  gDeltaDiffPsiDetElem0->Write("gDeltaDiffPsiDetElem0");
  gDeltaDiffThtDetElem0->Write("gDeltaDiffThtDetElem0");
  gDeltaDiffPhiDetElem0->Write("gDeltaDiffPhiDetElem0");
  myDetElemDeltaDiffCenX->Write();
  myDetElemDeltaDiffCenY->Write();
  myDetElemDeltaDiffCenZ->Write();
  myDetElemDeltaDiffRotX->Write();
  myDetElemDeltaDiffRotY->Write();
  myDetElemDeltaDiffRotZ->Write();
  hFile->Close();
  hFile->Delete();
}
