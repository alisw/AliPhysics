#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliITSRecPoint.h"
#include "AliESDEvent.h"
#include "AliTrackPointArray.h"
#include "AliITSgeomTGeo.h"
#include "AliITSTPArrayFit.h"
#include "AliESDfriend.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliGeomManager.h"
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TChain.h>
#include <TGeoGlobalMagField.h>
#include "AliESDInputHandlerRP.h"

/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysiTaskITSAlignQA
// AliAnalysisTaskSE to extract from ESD + ESDfriends 
// the track-to-point residuals and dE/dx vs, time for SDD modules
//
// Author: F. Prino, prino@to.infn.it
//*************************************************************************


#include "AliAnalysisTaskITSAlignQA.h"

ClassImp(AliAnalysisTaskITSAlignQA)
//______________________________________________________________________________
AliAnalysisTaskITSAlignQA::AliAnalysisTaskITSAlignQA() : AliAnalysisTaskSE("SDD Calibration"), 
  fOutput(0),
  fHistNEvents(0),
  fHistPtAccept(0),
  fDoSPDResiduals(kTRUE),
  fDoSDDResiduals(kTRUE),
  fDoSSDResiduals(kTRUE),
  fDoSDDdEdxCalib(kTRUE),
  fUseITSsaTracks(kFALSE),
  fLoadGeometry(kFALSE),
  fMinITSpts(3),
  fMinTPCpts(70),
  fMinPt(0.5),
  fNPtBins(8),
  fFitter(0),
  fRunNb(0),
  fOCDBLocation("local://$ALICE_ROOT/OCDB")
{
  //
  fFitter = new AliITSTPArrayFit(5);
  Double_t xbins[9]={0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};
  SetPtBinLimits(8,xbins);
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskITSAlignQA::~AliAnalysisTaskITSAlignQA(){
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fHistNEvents){
    delete fHistNEvents;
    fHistNEvents=0;    
  }
  for(Int_t i=0;i<kNSDDmods;i++){
    if(fHistSDDResidXvsX[i]){ delete fHistSDDResidXvsX[i]; fHistSDDResidXvsX[i]=0;}
    if(fHistSDDResidXvsZ[i]){ delete fHistSDDResidXvsX[i]; fHistSDDResidXvsX[i]=0;}
    if(fHistSDDResidZvsX[i]){ delete fHistSDDResidXvsX[i]; fHistSDDResidXvsX[i]=0;}
    if(fHistSDDResidZvsZ[i]){ delete fHistSDDResidXvsX[i]; fHistSDDResidXvsX[i]=0;}
    if(fHistSDDdEdxvsDrTime[i]){ delete fHistSDDdEdxvsDrTime[i]; fHistSDDdEdxvsDrTime[i]=0;}
    if(fHistSDDDrTimeAll[i]){ delete fHistSDDDrTimeAll[i]; fHistSDDDrTimeAll[i]=0;}
    if(fHistSDDDrTimeExtra[i]){ delete fHistSDDDrTimeExtra[i]; fHistSDDDrTimeExtra[i]=0;}
    if(fHistSDDDrTimeAttac[i]){ delete fHistSDDDrTimeAttac[i]; fHistSDDDrTimeAttac[i]=0;}
    for (int ix=2;ix--;) {
      if (fHProfSDDResidXvsXD[i][ix]) delete fHProfSDDResidXvsXD[i][ix]; fHProfSDDResidXvsXD[i][ix] = 0;
      if (fHProfSDDResidXvsZ[i][ix])  delete fHProfSDDResidXvsZ[i][ix];  fHProfSDDResidXvsZ[i][ix]  = 0;
      if (fHProfSDDDrTimevsXD[i][ix]) delete fHProfSDDDrTimevsXD[i][ix]; fHProfSDDDrTimevsXD[i][ix] = 0;
      if (fHProfSDDDrTimevsZ[i][ix])  delete fHProfSDDDrTimevsZ[i][ix];  fHProfSDDDrTimevsZ[i][ix]  = 0;
    }
  }
  if(fFitter){
    delete fFitter;
    fFitter=0;
  }
}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::UserCreateOutputObjects() {
  //

  if(fLoadGeometry) LoadGeometryFromOCDB();

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistPtAccept = new TH1F("hPtAccept","Pt distrib of accepted tracks",50,0.,5.);
  fHistPtAccept->Sumw2();
  fHistPtAccept->SetMinimum(0);
  fOutput->Add(fHistPtAccept);

  if(fDoSPDResiduals) CreateSPDHistos();
  if(fDoSDDResiduals || fDoSDDdEdxCalib) CreateSDDHistos();
  if(fDoSSDResiduals) CreateSSDHistos();

  PostData(1,fOutput);
}

//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CreateSPDHistos(){
  // Histos for SPD
  

  for(Int_t iMod=0; iMod<kNSPDmods; iMod++){
    fHistSPDResidX[iMod] = new TH2F(Form("hSPDResidX%d",iMod),
				    Form("hSPDResidX%d",iMod),
				    fNPtBins,fPtBinLimits,
				    250,-0.05,0.05);
    fHistSPDResidX[iMod]->Sumw2();
    fOutput->Add(fHistSPDResidX[iMod]);

    fHistSPDResidZ[iMod] = new TH2F(Form("hSPDResidZ%d",iMod),
				    Form("hSPDResidZ%d",iMod),
				    fNPtBins,fPtBinLimits,
				    250,-0.1,0.1);
    fHistSPDResidZ[iMod]->Sumw2();
    fOutput->Add(fHistSPDResidZ[iMod]);
  }
  return;
}

//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CreateSDDHistos(){
  // Histos for SDD

  for(Int_t iMod=0; iMod<kNSDDmods; iMod++){
    if(fDoSDDResiduals){
      fHistSDDResidX[iMod] = new TH2F(Form("hSDDResidX%d",iMod+kNSPDmods),
				      Form("hSDDResidX%d",iMod+kNSPDmods),
				      fNPtBins,fPtBinLimits,
				      300,-0.15,0.15);
      fHistSDDResidX[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidX[iMod]);
      
      fHistSDDResidZ[iMod] = new TH2F(Form("hSDDResidZ%d",iMod+kNSPDmods),
				      Form("hSDDResidZ%d",iMod+kNSPDmods),
				      fNPtBins,fPtBinLimits,
				      200,-0.1,0.1);
      fHistSDDResidZ[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidZ[iMod]);
      
      fHistSDDResidXvsX[iMod] = new TH2F(Form("hSDDResidXvsX%d",iMod+kNSPDmods),
					 Form("hSDDResidXvsX%d",iMod+kNSPDmods),
					 40,-3.5,3.5,300,-0.15,0.15);   
      fHistSDDResidXvsX[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidXvsX[iMod]);
      
      fHistSDDResidXvsZ[iMod] = new TH2F(Form("hSDDResidXvsZ%d",iMod+kNSPDmods),
					 Form("hSDDResidXvsZ%d",iMod+kNSPDmods),
					 10,-3.8,3.8,300,-0.15,0.15);   
      fHistSDDResidXvsZ[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidXvsZ[iMod]);
      
      fHistSDDResidZvsX[iMod] = new TH2F(Form("hSDDResidZvsX%d",iMod+kNSPDmods),
				      Form("hSDDResidZvsX%d",iMod+kNSPDmods),
				      40,-3.5,3.5,200,-0.1,0.1);   
      fHistSDDResidZvsX[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidZvsX[iMod]);
      
      fHistSDDResidZvsZ[iMod] = new TH2F(Form("hSDDResidZvsZ%d",iMod+kNSPDmods),
					 Form("hSDDResidZvsZ%d",iMod+kNSPDmods),
					 10,-3.8,3.8,200,-0.1,0.1);   
      fHistSDDResidZvsZ[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidZvsZ[iMod]);
      //
      for (int ix=0;ix<2;ix++) { // profile histos per side
	//
	char* hnm = Form("hpSDDResXvsXD%d_%d",iMod+kNSPDmods,ix);
	int nbX = 50, nbZ = 20, nbXOffs = 2, nbZOffs = 2;
	double xRange = 3.5085e4, zRange = 7.5264e4, xOffs = nbXOffs*xRange/nbX, zOffs = nbZOffs*zRange/nbZ;
	fHProfSDDResidXvsXD[iMod][ix] = new TProfile(hnm, hnm, nbX+2*nbXOffs, -xOffs, xRange + xOffs);
	fHProfSDDResidXvsXD[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDResidXvsXD[iMod][ix]);
	//
	hnm = Form("hpSDDDrTimevsXD%d_%d",iMod+kNSPDmods,ix);
	fHProfSDDDrTimevsXD[iMod][ix] = new TProfile(hnm, hnm, nbX+2*nbXOffs, -xOffs, xRange + xOffs);
	fHProfSDDDrTimevsXD[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDDrTimevsXD[iMod][ix]);
	//
	hnm = Form("hpSDDResXvsZ%d_%d",iMod+kNSPDmods,ix);
	fHProfSDDResidXvsZ[iMod][ix] = new TProfile(hnm, hnm, nbZ+2*nbZOffs, -(0.5*zRange+zOffs),(0.5*zRange+zOffs));
	fHProfSDDResidXvsZ[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDResidXvsZ[iMod][ix]);
	//
	hnm = Form("hpSDDDrTimevsZ%d_%d",iMod+kNSPDmods,ix);
	fHProfSDDDrTimevsZ[iMod][ix] = new TProfile(hnm, hnm, nbZ+2*nbZOffs, -(0.5*zRange+zOffs),(0.5*zRange+zOffs));
	fHProfSDDDrTimevsZ[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDDrTimevsZ[iMod][ix]);
	//
      }
    }

    if(fDoSDDdEdxCalib){
      fHistSDDdEdxvsDrTime[iMod] = new TH2F(Form("hSDDdEdxvsDrTime%d",iMod+kNSPDmods),
					    Form("hSDDdEdxvsDrTime%d",iMod+kNSPDmods),
					    16,0.,6400.,100,0.,300.);
      fHistSDDdEdxvsDrTime[iMod]->Sumw2();
      fOutput->Add(fHistSDDdEdxvsDrTime[iMod]);
    }

    fHistSDDDrTimeAll[iMod]=new TH1F(Form("hSDDDrTimeAll%d",iMod+kNSPDmods),
				     Form("hSDDDrTimeAll%d",iMod+kNSPDmods),
				     3200,0.,6400.);
    fHistSDDDrTimeAll[iMod]->Sumw2();
    fHistSDDDrTimeAll[iMod]->SetMinimum(0.);
    fOutput->Add(fHistSDDDrTimeAll[iMod]);

    fHistSDDDrTimeExtra[iMod]=new TH1F(Form("hSDDDrTimeExtra%d",iMod+kNSPDmods),
				       Form("hSDDDrTimeExtra%d",iMod+kNSPDmods),
				       3200,0.,6400.);
    fHistSDDDrTimeExtra[iMod]->Sumw2();
    fHistSDDDrTimeExtra[iMod]->SetMinimum(0.);
    fOutput->Add(fHistSDDDrTimeExtra[iMod]);

    fHistSDDDrTimeAttac[iMod]=new TH1F(Form("hSDDDrTimeAttac%d",iMod+kNSPDmods),
				       Form("hSDDDrTimeAttac%d",iMod+kNSPDmods),
				       3200,0.,6400.);
    fHistSDDDrTimeAttac[iMod]->Sumw2();
    fHistSDDDrTimeAttac[iMod]->SetMinimum(0.);
    fOutput->Add(fHistSDDDrTimeAttac[iMod]);

  }
  return;

}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CreateSSDHistos(){
  // Histos for SSD
  for(Int_t iMod=0; iMod<kNSSDmods; iMod++){
    fHistSSDResidX[iMod] = new TH2F(Form("hSSDResidX%d",iMod+kNSPDmods+kNSDDmods),
				    Form("hSSDResidX%d",iMod+kNSPDmods+kNSDDmods),
				    fNPtBins,fPtBinLimits,
				    250,-0.1,0.1);
    fHistSSDResidX[iMod]->Sumw2();
    fOutput->Add(fHistSSDResidX[iMod]);

    fHistSSDResidZ[iMod] = new TH2F(Form("hSSDResidZ%d",iMod+kNSPDmods+kNSDDmods),
				    Form("hSSDResidZ%d",iMod+kNSPDmods+kNSDDmods),
				    fNPtBins,fPtBinLimits,
				    250,-1.,1.);
    fHistSSDResidZ[iMod]->Sumw2();
    fOutput->Add(fHistSSDResidZ[iMod]);
  }
  return;
}
//______________________________________________________________________________
void AliAnalysisTaskITSAlignQA::UserExec(Option_t *)
{
  //
  AliESDEvent *esd = (AliESDEvent*) (InputEvent());

  if(!esd) {
    printf("AliAnalysisTaskITSAlignQA::Exec(): bad ESD\n");
    return;
  } 


  if(!ESDfriend()) {
    printf("AliAnalysisTaskITSAlignQA::Exec(): bad ESDfriend\n");
    return;
  }

  fHistNEvents->Fill(0);
  fFitter->SetBz(esd->GetMagneticField());

  const AliTrackPointArray *array = 0;
  Int_t ntracks = esd->GetNumberOfTracks();

  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    AliESDtrack * track = esd->GetTrack(itrack);
    if(!track) continue;
    if(!AcceptTrack(track)) continue;
    array = track->GetTrackPointArray();
    if(!array) continue;

    Int_t npts=array->GetNPoints();
    if(fDoSPDResiduals){ 
      FitAndFillSPD(1,array,npts,track);
      FitAndFillSPD(2,array,npts,track);
    }
    if(fDoSDDResiduals || fDoSDDdEdxCalib){
      FitAndFillSDDrphi(array,npts,track);
      FitAndFillSDDz(3,array,npts,track);
      FitAndFillSDDz(4,array,npts,track);
    }
    if(fDoSSDResiduals){ 
      FitAndFillSSD(5,array,npts,track);
      FitAndFillSSD(6,array,npts,track);
    }
  }

  PostData(1,fOutput);
  
}
//___________________________________________________________________________
Bool_t AliAnalysisTaskITSAlignQA::AcceptTrack(AliESDtrack * track){
  // track selection cuts
  Bool_t accept=kTRUE;
  if(fUseITSsaTracks){ 
      if(track->GetNcls(1)>0) accept=kFALSE;
  }else{
    if(track->GetNcls(1)<fMinTPCpts) accept=kFALSE;
  }
  if(track->GetNcls(0) < fMinITSpts) accept=kFALSE;
  Int_t trstatus=track->GetStatus();
  if(!(trstatus&AliESDtrack::kITSrefit)) accept=kFALSE;
  Float_t pt=track->Pt();
  if(pt<fMinPt) accept=kFALSE;
  if(accept) fHistPtAccept->Fill(pt);
  return accept;
}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::FitAndFillSPD(Int_t iLayer, const AliTrackPointArray *array, Int_t npts,AliESDtrack * track){
  // fit track and fills histos for SPD
  fFitter->AttachPoints(array,0, npts-1); 
  Int_t iPtSPD[4],modIdSPD[4];
  Int_t nPtSPD=0;
  Double_t resGlo[3],resLoc[3];
  for(Int_t ipt=0; ipt<npts; ipt++) {
    AliTrackPoint point;
    Int_t modId;
    array->GetPoint(point,ipt);
    Int_t volId = point.GetVolumeID();
    Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==iLayer){
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      iPtSPD[nPtSPD] = ipt;
      modIdSPD[nPtSPD] = modId;
      ++nPtSPD;
      fFitter->SetCovIScale(ipt,1e-4); 
    }
  }
  if(nPtSPD>0){
    fFitter->Fit(track->Charge(),track->Pt(),0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSPD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSPD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSPD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSPD[ip];
      fHistSPDResidX[index]->Fill(track->Pt(),resLoc[0]);
      fHistSPDResidZ[index]->Fill(track->Pt(),resLoc[2]);
    }
  }    
}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::FitAndFillSDDrphi(const AliTrackPointArray *array, Int_t npts, AliESDtrack * track){
  // fit track and fills histos for SDD along rphi (drift coord.)
  Double_t dedx[4];
  track->GetITSdEdxSamples(dedx);

  fFitter->AttachPoints(array,0, npts-1); 
  Int_t iPtSDD[4],modIdSDD[4],modSide[4];
  Double_t xLocSDD[4],zLocSDD[4],drTime[4];
  Int_t nPtSDD=0;
  Int_t nPtSSDSPD=0;
  Double_t resGlo[3],resLoc[3];
  Float_t  posGloF[3];
  Double_t posGlo[3],posLoc[3];

  for(Int_t ipt=0; ipt<npts; ipt++) {
    AliTrackPoint point;
    Int_t modId;
    array->GetPoint(point,ipt);
    Int_t volId = point.GetVolumeID();
    Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==3 || layerId==4){
      drTime[nPtSDD] = point.GetDriftTime();
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      Int_t index=modId-kNSPDmods;
      fHistSDDDrTimeAll[index]->Fill(drTime[nPtSDD]);
      if(point.IsExtra()) fHistSDDDrTimeExtra[index]->Fill(drTime[nPtSDD]);
      else fHistSDDDrTimeAttac[index]->Fill(drTime[nPtSDD]);
      Float_t dedxLay=dedx[layerId-3];
      if(dedxLay>1.) fHistSDDdEdxvsDrTime[index]->Fill(drTime[nPtSDD],dedxLay);
      iPtSDD[nPtSDD] = ipt;
      modIdSDD[nPtSDD] = modId;
      modSide[nPtSDD] = point.GetClusterType()&BIT(16) ? 0:1; 
      point.GetXYZ(posGloF);
      for(Int_t icoor=0;icoor<3;icoor++) posGlo[icoor]=posGloF[icoor];
      AliITSgeomTGeo::GlobalToLocal(modId,posGlo,posLoc);
      xLocSDD[nPtSDD]=posLoc[0];
      zLocSDD[nPtSDD]=posLoc[2];
      ++nPtSDD;
      fFitter->SetCovIScale(ipt,1e-4); // scaling for inverted errors of SDD
    }else{
      ++nPtSSDSPD;
    }
  }
  if(nPtSDD>0 && nPtSSDSPD>=2){
    fFitter->Fit(track->Charge(),track->Pt(),0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSDD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSDD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSDD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSDD[ip]-kNSPDmods;
      fHistSDDResidX[index]->Fill(track->Pt(),resLoc[0]);
      fHistSDDResidXvsX[index]->Fill(xLocSDD[ip],resLoc[0]);
      fHistSDDResidXvsZ[index]->Fill(zLocSDD[ip],resLoc[0]);
      //
      double cf = modSide[ip] ? 1.e4:-1.e4;
      double xMeas = cf*xLocSDD[0];            // measured coordinate in microns
      double xRes  = cf*resLoc[0];             // X residual in microns
      double xDriftTrue  = 3.5085e4 - (xMeas + xRes);   // "true" drift distance
      //
      fHProfSDDResidXvsXD[index][modSide[ip]]->Fill(xDriftTrue, xRes);
      fHProfSDDResidXvsZ[index][modSide[ip]]->Fill(zLocSDD[ip], xRes);
      fHProfSDDDrTimevsXD[index][modSide[ip]]->Fill(xDriftTrue, drTime[ip]);
      fHProfSDDDrTimevsZ[index][modSide[ip]]->Fill(zLocSDD[ip], drTime[ip]);      
    }
  }
}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::FitAndFillSDDz(Int_t iLayer, const AliTrackPointArray *array, Int_t npts, AliESDtrack * track){
  // fit track and fills histos for SDD along z

  fFitter->AttachPoints(array,0, npts-1); 
  Int_t iPtSDD[4],modIdSDD[4];
  Double_t xLocSDD[4],zLocSDD[4];
  Int_t nPtSDD=0;
  Double_t resGlo[3],resLoc[3];
  Float_t  posGloF[3];
  Double_t posGlo[3],posLoc[3];
  for(Int_t ipt=0; ipt<npts; ipt++) {
    AliTrackPoint point;
    Int_t modId;
    array->GetPoint(point,ipt);
    Int_t volId = point.GetVolumeID();
    Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==iLayer){
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      iPtSDD[nPtSDD] = ipt;
      modIdSDD[nPtSDD] = modId;
      point.GetXYZ(posGloF);
      for(Int_t icoor=0;icoor<3;icoor++) posGlo[icoor]=posGloF[icoor];
      AliITSgeomTGeo::GlobalToLocal(modId,posGlo,posLoc);
      xLocSDD[nPtSDD]=posLoc[0];
      zLocSDD[nPtSDD]=posLoc[2];
      ++nPtSDD;
      fFitter->SetCovIScale(ipt,1e-4); // scaling for inverted errors of SDD
    }
  }
  if(nPtSDD>0){
    fFitter->Fit(track->Charge(),track->Pt(),0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSDD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSDD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSDD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSDD[ip]-kNSPDmods;
      fHistSDDResidZ[index]->Fill(track->Pt(),resLoc[2]);
      fHistSDDResidZvsX[index]->Fill(xLocSDD[ip],resLoc[2]);
      fHistSDDResidZvsZ[index]->Fill(zLocSDD[ip],resLoc[2]);
    }
  }
}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::FitAndFillSSD(Int_t iLayer, const AliTrackPointArray *array, Int_t npts,AliESDtrack * track){
  // fit track and fills histos for SSD
  fFitter->AttachPoints(array,0, npts-1); 
  Int_t iPtSSD[4],modIdSSD[4];
  Int_t nPtSSD=0;
  Double_t resGlo[3],resLoc[3];
  for(Int_t ipt=0; ipt<npts; ipt++) {
    AliTrackPoint point;
    Int_t modId;
    array->GetPoint(point,ipt);
    Int_t volId = point.GetVolumeID();
    Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==iLayer){
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      iPtSSD[nPtSSD] = ipt;
      modIdSSD[nPtSSD] = modId;
      ++nPtSSD;
      fFitter->SetCovIScale(ipt,1e-4); 
    }  
  }
  if(nPtSSD>0){
    fFitter->Fit(track->Charge(),track->Pt(),0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSSD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSSD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSSD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSSD[ip]-kNSPDmods-kNSDDmods;
      fHistSSDResidX[index]->Fill(track->Pt(),resLoc[0]);
      fHistSSDResidZ[index]->Fill(track->Pt(),resLoc[2]);
    }
  }
}
//______________________________________________________________________________
void AliAnalysisTaskITSAlignQA::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(2));
  return;
}


//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::LoadGeometryFromOCDB(){
  //method to get the gGeomanager
  // it is called at the CreatedOutputObject stage
  // to comply with the CAF environment
  AliInfo("Loading geometry");

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(fOCDBLocation.Data());
  man->SetRun(fRunNb);
  AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
  AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
  AliGeomManager::GetNalignable("ITS");
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");
}






