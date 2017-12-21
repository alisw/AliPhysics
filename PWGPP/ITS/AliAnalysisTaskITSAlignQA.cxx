#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliITSRecPoint.h"
#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliDAQ.h"
#include "AliTrackPointArray.h"
#include "AliITSgeomTGeo.h"
#include "AliITSTPArrayFit.h"
#include "AliESDfriend.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliGeomManager.h"
#include "AliMultiplicity.h"
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TChain.h>
#include <TGeoGlobalMagField.h>
#include "AliESDInputHandlerRP.h"
#include "AliITSSumTP.h"
#include "AliMagF.h"

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
  fDoSDDVDriftCalib(kTRUE),
  fDoSDDDriftTime(kTRUE),
  fDoFillTPTree(kFALSE),
  fUseITSsaTracks(kFALSE),
  fLoadGeometry(kFALSE),
  fUseVertex(kFALSE),
  fUseVertexForZOnly(kFALSE),
  fUseTPCMomentum(kFALSE),
  fUseTOFTiming(10.),
  fMinVtxContributors(5),
  fRemovePileupWithSPD(kTRUE),
  fMinITSpts(3),
  fMinTPCpts(70),
  fMinPt(0.5),
  fNPtBins(8),
  fMinMult(0),
  fMaxMult(1e9),
  fCutDCAXY(100.),
  fCutDCAZ(100.),
  fFitter(0),
  fITSSumTP(),
  fTPTree(),
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
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    delete fITSSumTP;
  }
  delete fFitter;
  delete fTPTree;
  //
}
//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::UserCreateOutputObjects() {
  //

  if(fLoadGeometry) LoadGeometryFromOCDB();

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",kNEvStatBins,-0.5,kNEvStatBins-0.5);
  //fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(kEvAll+1,"All Events");
  fHistNEvents->GetXaxis()->SetBinLabel(kEvCnt+1,"After Centrality cut");
  fHistNEvents->GetXaxis()->SetBinLabel(kEvVtx+1,"After Vertex cut");
  fHistNEvents->GetXaxis()->SetBinLabel(kEvPlp+1,"After Pileup cut");
  fHistNEvents->GetXaxis()->SetBinLabel(kNTracks+1,"Tracks Accepted");
  fOutput->Add(fHistNEvents);

  fHistPtAccept = new TH1F("hPtAccept","Pt distrib of accepted tracks",50,0.,5.);
  //fHistPtAccept->Sumw2();
  fHistPtAccept->SetMinimum(0);
  fOutput->Add(fHistPtAccept);

  if(fDoSPDResiduals) CreateSPDHistos();
  if(fDoSDDResiduals || fDoSDDdEdxCalib || fDoSDDVDriftCalib || fDoSDDDriftTime) CreateSDDHistos();
  if(fDoSSDResiduals) CreateSSDHistos();
  //
  if (fDoFillTPTree) {
    TFile* troutf = OpenFile(2);
    if (!troutf) {
      AliFatal("Failed to open output file for AliITSSumTP tree");
      exit(1);
    }
    fITSSumTP = new AliITSSumTP();
    fTPTree = new TTree("ITSSumTP","ITS TP Summary");
    fTPTree->Branch("AliITSSumTP","AliITSSumTP",&fITSSumTP);
    PostData(2,fTPTree);
  }
  //
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
    //fHistSPDResidX[iMod]->Sumw2();
    fOutput->Add(fHistSPDResidX[iMod]);

    fHistSPDResidZ[iMod] = new TH2F(Form("hSPDResidZ%d",iMod),
				    Form("hSPDResidZ%d",iMod),
				    fNPtBins,fPtBinLimits,
				    250,-0.1,0.1);
    //fHistSPDResidZ[iMod]->Sumw2();
    fOutput->Add(fHistSPDResidZ[iMod]);
  }
  return;
}

//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CreateSDDHistos(){
  // Histos for SDD

  for(Int_t iMod=0; iMod<kNSDDmods; iMod++){
    if (fDoSDDResiduals) {
      fHistSDDResidX[iMod] = new TH2F(Form("hSDDResidX%d",iMod+kNSPDmods),
				      Form("hSDDResidX%d",iMod+kNSPDmods),
				      fNPtBins,fPtBinLimits,
				      300,-0.15,0.15);
      //fHistSDDResidX[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidX[iMod]);
      
      fHistSDDResidZ[iMod] = new TH2F(Form("hSDDResidZ%d",iMod+kNSPDmods),
				      Form("hSDDResidZ%d",iMod+kNSPDmods),
				      fNPtBins,fPtBinLimits,
				      200,-0.1,0.1);
      //fHistSDDResidZ[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidZ[iMod]);
      
      fHistSDDResidXvsX[iMod] = new TH2F(Form("hSDDResidXvsX%d",iMod+kNSPDmods),
					 Form("hSDDResidXvsX%d",iMod+kNSPDmods),
					 40,-3.5,3.5,300,-0.15,0.15);   
      //fHistSDDResidXvsX[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidXvsX[iMod]);
      
      fHistSDDResidXvsZ[iMod] = new TH2F(Form("hSDDResidXvsZ%d",iMod+kNSPDmods),
					 Form("hSDDResidXvsZ%d",iMod+kNSPDmods),
					 10,-3.8,3.8,300,-0.15,0.15);   
      //fHistSDDResidXvsZ[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidXvsZ[iMod]);
      
      fHistSDDResidZvsX[iMod] = new TH2F(Form("hSDDResidZvsX%d",iMod+kNSPDmods),
				      Form("hSDDResidZvsX%d",iMod+kNSPDmods),
				      40,-3.5,3.5,200,-0.1,0.1);   
      //fHistSDDResidZvsX[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidZvsX[iMod]);
      
      fHistSDDResidZvsZ[iMod] = new TH2F(Form("hSDDResidZvsZ%d",iMod+kNSPDmods),
					 Form("hSDDResidZvsZ%d",iMod+kNSPDmods),
					 10,-3.8,3.8,200,-0.1,0.1);   
      //fHistSDDResidZvsZ[iMod]->Sumw2();
      fOutput->Add(fHistSDDResidZvsZ[iMod]);
      //
    }
    //
    if (fDoSDDVDriftCalib) {
      for (int ix=0;ix<2;ix++) { // profile histos per side
	//
	char* hnm = Form("hpSDDResXvsXD%d_%d",iMod+kNSPDmods,ix);
	int nbX = 50, nbZ = 20, nbXOffs = 2, nbZOffs = 2;
	double xRange = 3.5085e4, zRange = 7.5264e4, xOffs = nbXOffs*xRange/nbX, zOffs = nbZOffs*zRange/nbZ;
	fHProfSDDResidXvsXD[iMod][ix] = new TProfile(hnm, hnm, nbX+2*nbXOffs, -xOffs, xRange + xOffs);
	//fHProfSDDResidXvsXD[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDResidXvsXD[iMod][ix]);
	//
	hnm = Form("hpSDDDrTimevsXD%d_%d",iMod+kNSPDmods,ix);
	fHProfSDDDrTimevsXD[iMod][ix] = new TProfile(hnm, hnm, nbX+2*nbXOffs, -xOffs, xRange + xOffs);
	//fHProfSDDDrTimevsXD[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDDrTimevsXD[iMod][ix]);
	//
	hnm = Form("hpSDDResXvsZ%d_%d",iMod+kNSPDmods,ix);
	fHProfSDDResidXvsZ[iMod][ix] = new TProfile(hnm, hnm, nbZ+2*nbZOffs, -(0.5*zRange+zOffs),(0.5*zRange+zOffs));
	//fHProfSDDResidXvsZ[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDResidXvsZ[iMod][ix]);
	//
	hnm = Form("hpSDDDrTimevsZ%d_%d",iMod+kNSPDmods,ix);
	fHProfSDDDrTimevsZ[iMod][ix] = new TProfile(hnm, hnm, nbZ+2*nbZOffs, -(0.5*zRange+zOffs),(0.5*zRange+zOffs));
	//fHProfSDDDrTimevsZ[iMod][ix]->Sumw2();
	fOutput->Add(fHProfSDDDrTimevsZ[iMod][ix]);
	//
      }
    }
    
    if(fDoSDDdEdxCalib){
      fHistSDDdEdxvsDrTime[iMod] = new TH2F(Form("hSDDdEdxvsDrTime%d",iMod+kNSPDmods),
					    Form("hSDDdEdxvsDrTime%d",iMod+kNSPDmods),
					    16,0.,6400.,100,0.,300.);
      //fHistSDDdEdxvsDrTime[iMod]->Sumw2();
      fOutput->Add(fHistSDDdEdxvsDrTime[iMod]);
    }
    //
    if (fDoSDDDriftTime) {
      fHistSDDDrTimeAll[iMod]=new TH1F(Form("hSDDDrTimeAll%d",iMod+kNSPDmods),
				       Form("hSDDDrTimeAll%d",iMod+kNSPDmods),
				       3200,0.,6400.);
      //fHistSDDDrTimeAll[iMod]->Sumw2();
      fHistSDDDrTimeAll[iMod]->SetMinimum(0.);
      fOutput->Add(fHistSDDDrTimeAll[iMod]);
      
      fHistSDDDrTimeExtra[iMod]=new TH1F(Form("hSDDDrTimeExtra%d",iMod+kNSPDmods),
					 Form("hSDDDrTimeExtra%d",iMod+kNSPDmods),
					 3200,0.,6400.);
      //fHistSDDDrTimeExtra[iMod]->Sumw2();
      fHistSDDDrTimeExtra[iMod]->SetMinimum(0.);
      fOutput->Add(fHistSDDDrTimeExtra[iMod]);
      
      fHistSDDDrTimeAttac[iMod]=new TH1F(Form("hSDDDrTimeAttac%d",iMod+kNSPDmods),
					 Form("hSDDDrTimeAttac%d",iMod+kNSPDmods),
					 3200,0.,6400.);
      //fHistSDDDrTimeAttac[iMod]->Sumw2();
      fHistSDDDrTimeAttac[iMod]->SetMinimum(0.);
      fOutput->Add(fHistSDDDrTimeAttac[iMod]);
    }
  }
  return;
  //
}

//___________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CreateSSDHistos(){
  // Histos for SSD
  for(Int_t iMod=0; iMod<kNSSDmods; iMod++){
    fHistSSDResidX[iMod] = new TH2F(Form("hSSDResidX%d",iMod+kNSPDmods+kNSDDmods),
				    Form("hSSDResidX%d",iMod+kNSPDmods+kNSDDmods),
				    fNPtBins,fPtBinLimits,
				    250,-0.1,0.1);
    //fHistSSDResidX[iMod]->Sumw2();
    fOutput->Add(fHistSSDResidX[iMod]);

    fHistSSDResidZ[iMod] = new TH2F(Form("hSSDResidZ%d",iMod+kNSPDmods+kNSDDmods),
				    Form("hSSDResidZ%d",iMod+kNSPDmods+kNSDDmods),
				    fNPtBins,fPtBinLimits,
				    250,-1.,1.);
    //fHistSSDResidZ[iMod]->Sumw2();
    fOutput->Add(fHistSSDResidZ[iMod]);
  }
  return;
}
//______________________________________________________________________________
void AliAnalysisTaskITSAlignQA::UserExec(Option_t *)
{
  //
  static AliTrackPointArray* arrayITS = 0;
  AliTrackPointArray* arrayITSNoVtx = 0;
  //
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (fITSSumTP) fITSSumTP->Reset();

  if(!esd) {
    AliInfo("No ESD");
    return;
  } 
  AliESDfriend* fr = dynamic_cast<AliESDfriend*>(ESDfriend());
  if(!fr || !fr->GetNumberOfTracks()) {
    AliInfo("No or empty ESDfriend");
    return;
  }
  //
  static Bool_t firstCheck = kTRUE;
  if (firstCheck) {
    //    
    if (TMath::Abs(esd->GetCurrentL3())<300) { // no field
      SetMinPt(0.005);
      AliInfo("No magnetic field: eliminating pt cut");
    }
    const AliESDRun *esdrn = esd->GetESDRun();
    if (!esdrn) return;
    Int_t activeDetectors = esdrn->GetDetectorsInReco();
    if ( !(activeDetectors & AliDAQ::kTPC) ) {
      AliInfo("No TPC, suppress TPC points request");
      SetUseITSstandaloneTracks(kTRUE);
      SetUseTPCMomentum(kFALSE);
    }
    firstCheck = kFALSE;
  }
  //
  fHistNEvents->Fill(kEvAll);
  //
  if (!AcceptCentrality(esd)) return;
  fHistNEvents->Fill(kEvCnt);

  const AliESDVertex* vtx=0,*vtxSPD=0;
  vtx    = esd->GetPrimaryVertex();
  vtxSPD = esd->GetPrimaryVertexSPD();
  //
  if (fUseVertex) {  // check the vertex if it is requested as an extra point
    if (!AcceptVertex(vtx,vtxSPD)) return;
  }

  fHistNEvents->Fill(kEvVtx);
  if (fRemovePileupWithSPD){
    // skip events tagged by SPD as pileup
    if(esd->IsPileupFromSPD()) return;
  }
  fHistNEvents->Fill(kEvPlp);

  //
  fFitter->SetBz(esd->GetMagneticField());

  const AliTrackPointArray *array = 0;
  Int_t ntracks = esd->GetNumberOfTracks();

  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    //
    if (arrayITS) {delete arrayITS; arrayITS = 0;}  // reset points from previous tracks 
    arrayITSNoVtx = 0;
    //
    AliESDtrack * track = esd->GetTrack(itrack);
    if(!track) continue;
    array = track->GetTrackPointArray();
    if(!array) continue;
    if(!AcceptTrack(track, vtx)) continue;
    arrayITS = PrepareTrack(array, vtx);
    if (fITSSumTP) {
      arrayITSNoVtx = PrepareTrack(array, 0);
      arrayITSNoVtx->SetUniqueID(itrack);
      fITSSumTP->AddTrack(arrayITSNoVtx);
    }
    //
    fHistNEvents->Fill(kNTracks);
    //
    Int_t npts  = arrayITS->GetNPoints();
    Int_t npts1 = fUseVertexForZOnly ? npts-1 : npts;
    //
    if(fDoSPDResiduals){ 
      FitAndFillSPD(1,arrayITS,npts1,track);
      FitAndFillSPD(2,arrayITS,npts1,track);
    }
    if(fDoSDDResiduals || fDoSDDdEdxCalib || fDoSDDVDriftCalib || fDoSDDDriftTime) {
      FitAndFillSDDrphi(arrayITS,npts,track);
      if (fDoSDDResiduals) {
	FitAndFillSDDz(3,arrayITS,npts1,track);
	FitAndFillSDDz(4,arrayITS,npts1,track);
      }
    }
    if(fDoSSDResiduals){ 
      FitAndFillSSD(5,arrayITS,npts1,track);
      FitAndFillSSD(6,arrayITS,npts1,track);
    }
  }
  //
  if (fITSSumTP) { // store vertex and mometum info
    double xyz[3]={0,0,0};
    fITSSumTP->SetVertex(vtx);
    TObjArray& tps = fITSSumTP->GetTracks();
    int ntp = tps.GetEntriesFast();
    fITSSumTP->BookNTracks(ntp);
    for (int it=ntp;it--;) {
      AliTrackPointArray* tp = (AliTrackPointArray*)tps[it];
      if (!tp) continue;
      AliESDtrack* esdTr = esd->GetTrack(tp->GetUniqueID());
      double crv =  esdTr->GetC(esd->GetMagneticField());
      double crve = TMath::Sqrt(esdTr->GetSigma1Pt2()) * esd->GetMagneticField()*kB2C;
      fITSSumTP->SetCrvGlo(it,crv);
      fITSSumTP->SetCrvGloErr(it,crve);
      const AliExternalTrackParam* inTPC =  esdTr->GetTPCInnerParam(); // TPC track at vtx
      if (inTPC) {
	 crv =  inTPC->GetC(esd->GetMagneticField());
	 crve = TMath::Sqrt(inTPC->GetSigma1Pt2()) * TMath::Abs(esd->GetMagneticField()*kB2C);
	 fITSSumTP->SetCrvTPC(it,crv);
	 fITSSumTP->SetCrvTPCErr(it,crve);
      }
      inTPC = esdTr->GetInnerParam();  // TPC track at the inner wall
      if (inTPC) {
	inTPC->GetXYZ(xyz);
	fITSSumTP->SetTPCInnerXYZ(it,xyz);
      }
    }
    fITSSumTP->SetUniqueID(fCurrentRunNumber);
    if (ntp) fTPTree->Fill();
    CopyUserInfo();
  }

  //
  PostData(1,fOutput);
  
}

//___________________________________________________________________________
Bool_t AliAnalysisTaskITSAlignQA::AcceptTrack(const AliESDtrack * track, const AliESDVertex* vtx)
{
  // track selection cuts
  Bool_t accept=kTRUE;
  if(fUseITSsaTracks){ 
    if(track->GetNcls(1)>0) accept=kFALSE;
  }else{
    if (!track->IsOn(AliESDtrack::kTPCrefit) || track->GetNcls(1)<fMinTPCpts) accept=kFALSE;
    if (fUseTOFTiming>0) {
      if (TMath::Abs(track->GetTOFExpTDiff())>fUseTOFTiming) accept=kFALSE;
    }
  }
  if(track->GetNcls(0) < fMinITSpts) accept=kFALSE;
  Int_t trstatus=track->GetStatus();
  if(!(trstatus&AliESDtrack::kITSrefit)) accept=kFALSE;
  Float_t pt = 0;
  if (fUseTPCMomentum && track->IsOn(AliESDtrack::kTPCrefit)) pt = track->GetTPCInnerParam()->Pt();
  else pt = track->Pt();
  //
  if(pt<fMinPt) accept=kFALSE;
  //
  // if vertex constraint is used, apply soft DCA cut
  if (vtx) {
    Double_t dz[2],cov[3];
    AliExternalTrackParam trc = *track;
    if (!trc.PropagateToDCA(vtx, fFitter->GetBz(), 3.0, dz, cov)) accept=kFALSE;
    else {
      if (dz[0]*dz[0]/(1e-4+cov[0])>fCutDCAXY) accept=kFALSE;
      if (dz[1]*dz[1]/(4e-4+cov[2])>fCutDCAZ)  accept=kFALSE;
    }
  }
  //
  if(accept) fHistPtAccept->Fill(pt);
  return accept;
}

//___________________________________________________________________________
Bool_t AliAnalysisTaskITSAlignQA::AcceptVertex(const AliESDVertex * vtx, const AliESDVertex * vtxSPD) {
  // vertex selection cuts
  if (!vtx || vtx->GetStatus()<1) return kFALSE;
  if (!vtxSPD || vtxSPD->GetStatus()<1) return kFALSE;
  if (vtx->GetNContributors()<fMinVtxContributors) return kFALSE;
  if (TMath::Abs(vtx->GetZ()-vtxSPD->GetZ())>0.3) return kFALSE;
  return kTRUE;
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
    if (volId == kVtxSensVID) continue; // this is a vertex constraint
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
    double pt = (fUseTPCMomentum && track->IsOn(AliESDtrack::kTPCrefit)) ? track->GetTPCInnerParam()->Pt() : track->Pt();
    fFitter->Fit(track->Charge(),pt,0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSPD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSPD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSPD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSPD[ip];
      fHistSPDResidX[index]->Fill(pt,resLoc[0]);
      fHistSPDResidZ[index]->Fill(pt,resLoc[2]);
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
    if (volId == kVtxSensVID) continue; // this is a vertex constraint
    Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if(layerId==3 || layerId==4){
      drTime[nPtSDD] = point.GetDriftTime();
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      Int_t index=modId-kNSPDmods;
      if (fDoSDDDriftTime) {
	fHistSDDDrTimeAll[index]->Fill(drTime[nPtSDD]);
	if(point.IsExtra()) fHistSDDDrTimeExtra[index]->Fill(drTime[nPtSDD]);
	else fHistSDDDrTimeAttac[index]->Fill(drTime[nPtSDD]);
      }
      if (fDoSDDdEdxCalib) {
	Float_t dedxLay=dedx[layerId-3];
	if(dedxLay>1.) fHistSDDdEdxvsDrTime[index]->Fill(drTime[nPtSDD],dedxLay);
      }
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
    double pt = (fUseTPCMomentum && track->IsOn(AliESDtrack::kTPCrefit)) ? track->GetTPCInnerParam()->Pt() : track->Pt();
    fFitter->Fit(track->Charge(),pt,0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSDD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSDD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSDD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSDD[ip]-kNSPDmods;
      if (fDoSDDResiduals) {
	fHistSDDResidX[index]->Fill(pt,resLoc[0]);
	fHistSDDResidXvsX[index]->Fill(xLocSDD[ip],resLoc[0]);
	fHistSDDResidXvsZ[index]->Fill(zLocSDD[ip],resLoc[0]);
      }
      //
      if (fDoSDDVDriftCalib) {
	double cf = modSide[ip] ? 1.e4:-1.e4;
	double xMeas = cf*xLocSDD[ip];            // measured coordinate in microns
	double xRes  = cf*resLoc[0];             // X residual in microns
	double xDriftTrue  = 3.5085e4 - (xMeas + xRes);   // "true" drift distance
	//
	fHProfSDDResidXvsXD[index][modSide[ip]]->Fill(xDriftTrue, xRes);
	fHProfSDDResidXvsZ[index][modSide[ip]]->Fill(zLocSDD[ip]*1e4, xRes);
	fHProfSDDDrTimevsXD[index][modSide[ip]]->Fill(xDriftTrue, drTime[ip]);
	fHProfSDDDrTimevsZ[index][modSide[ip]]->Fill(zLocSDD[ip]*1e4, drTime[ip]);      
      }
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
    if (volId == kVtxSensVID) continue; // this is a vertex constraint
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
    double pt = (fUseTPCMomentum && track->IsOn(AliESDtrack::kTPCrefit)) ? track->GetTPCInnerParam()->Pt() : track->Pt();
    fFitter->Fit(track->Charge(),pt,0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSDD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSDD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSDD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSDD[ip]-kNSPDmods;
      fHistSDDResidZ[index]->Fill(pt,resLoc[2]);
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
    if (volId == kVtxSensVID) continue; // this is a vertex constraint
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
    double pt = (fUseTPCMomentum && track->IsOn(AliESDtrack::kTPCrefit)) ? track->GetTPCInnerParam()->Pt() : track->Pt();
    fFitter->Fit(track->Charge(),pt,0.);
    Double_t chi2=fFitter->GetChi2NDF();
    if ( chi2<0 || chi2>1e4 ) return; // fit failed, abandon this track
    for (Int_t ip=0; ip<nPtSSD;ip++) {
      fFitter->GetResiduals(resGlo,iPtSSD[ip]);
      TGeoHMatrix *mcurr = AliITSgeomTGeo::GetMatrix(modIdSSD[ip]);
      mcurr->MasterToLocalVect(resGlo,resLoc);
      Int_t index=modIdSSD[ip]-kNSPDmods-kNSDDmods;
      fHistSSDResidX[index]->Fill(pt,resLoc[0]);
      fHistSSDResidZ[index]->Fill(pt,resLoc[2]);
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
  if(fHistNEvents){
    AliInfo(Form("Number of analyzed events = %d, %d tracks accepted",
		 (Int_t)fHistNEvents->GetBinContent(kEvAcc+1),(Int_t)fHistNEvents->GetBinContent(kNTracks+1)));
  }else{
    printf("Warning: pointer to fHistNEvents is NULL\n");
  }
  //
  if (fDoFillTPTree) CreateUserInfo();
  //
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
  if(obj){
    AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
    AliGeomManager::GetNalignable("ITS");
    AliGeomManager::ApplyAlignObjsFromCDB("ITS");
  }
  else AliFatal("Geometry object not found in OCDB");
}


//______________________________________________________________________________________
AliTrackPointArray* AliAnalysisTaskITSAlignQA::PrepareTrack(const AliTrackPointArray* inp, const AliESDVertex* vtx)
{
  // Extract from the global TrackPointArray the ITS part and optionally add vertex as the last measured point
  //
  int npts = inp->GetNPoints();
  int modID=0,nptITS = 0;
  int itsRefs[24];
  const UShort_t *vids = inp->GetVolumeID();
  for(int ipt=0; ipt<npts; ipt++) { // count ITS points
    if (vids[ipt]<=0) continue;
    int layerId = AliGeomManager::VolUIDToLayer(vids[ipt],modID);
    if(layerId<1 || layerId>6) continue;
    itsRefs[nptITS++] = ipt;
  }
  //
  AliTrackPointArray *trackCopy = new AliTrackPointArray(nptITS + (vtx ? 1:0)); // reserve extra space if vertex provided
  AliTrackPoint point;
  for(int ipt=0; ipt<nptITS; ipt++) {
    inp->GetPoint(point,itsRefs[ipt]);
    trackCopy->AddPoint(ipt,&point);
  }
  //
  if (vtx) {
    PrepareVertexConstraint(vtx,point);
    trackCopy->AddPoint(nptITS,&point); // add vertex constraint as a last point
  }
  return trackCopy;
}

//_______________________________________________________________________________________
void AliAnalysisTaskITSAlignQA::PrepareVertexConstraint(const AliESDVertex* vtx, AliTrackPoint &point)
{
  // convert vertex to measured point with dummy VID
  if (!vtx) return;
  //
  double cmat[6];
  float cmatF[6];
  point.SetVolumeID(kVtxSensVID);
  //
  vtx->GetCovMatrix(cmat);
  cmatF[0] = cmat[0]; // xx
  cmatF[1] = cmat[1]; // xy
  cmatF[2] = cmat[3]; // xz
  cmatF[3] = cmat[2]; // yy
  cmatF[4] = cmat[4]; // yz
  cmatF[5] = cmat[5]; // zz
  point.SetXYZ(vtx->GetX(),vtx->GetY(),vtx->GetZ(), cmatF);
}


//_______________________________________________________________________________________
Bool_t AliAnalysisTaskITSAlignQA::AcceptCentrality(const AliESDEvent *esd) const
{
  // check if events is in the required multiplicity range
  //
  const AliMultiplicity *alimult = esd->GetMultiplicity();
  Int_t nclsSPDouter=0;
  if(alimult) nclsSPDouter = alimult->GetNumberOfITSClusters(1);
  if(nclsSPDouter<fMinMult || nclsSPDouter>fMaxMult) return kFALSE;
  //
  return kTRUE;
}

//_______________________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CreateUserInfo()
{
  // if needed, set user info of the output tree
  if (!fTPTree) {
    AliError("TrackPoints summary tree does not exist"); 
    return;
  }
  TList* uInfo = fTPTree->GetUserInfo();
  TMap  *cdbMapCopy  = (TMap*)uInfo->FindObject("cdbMap");
  TList *cdbListCopy = (TList*)uInfo->FindObject("cdbList"); 
  TList *bzList      = (TList*)uInfo->FindObject("BzkGauss"); 
  if (cdbMapCopy && cdbListCopy && bzList) return; //already done
  //
  const TMap *cdbMap = AliCDBManager::Instance()->GetStorageMap();	 
  const TList *cdbList = AliCDBManager::Instance()->GetRetrievedIds();	 
  //
  cdbMapCopy = new TMap(cdbMap->GetEntries());	 
  cdbMapCopy->SetOwner(1);	 
  cdbMapCopy->SetName("cdbMap");	 
  TIter iter(cdbMap->GetTable());	 
  //
  TPair* pair = 0;	 
  while((pair = dynamic_cast<TPair*> (iter.Next()))){	 
    TObjString* keyStr = dynamic_cast<TObjString*> (pair->Key());	 
    TObjString* valStr = dynamic_cast<TObjString*> (pair->Value());
    if (keyStr && valStr) {
      cdbMapCopy->Add(new TObjString(keyStr->GetName()), new TObjString(valStr->GetName()));	 
      AliInfo(Form("Add %s : %s to cdbMap of ITSTPUserInfo",keyStr->GetName(),valStr->GetName()));
    }
  }	 
  //
  cdbListCopy = new TList();	 
  cdbListCopy->SetOwner(1);	 
  cdbListCopy->SetName("cdbList");	 
  // 
  TIter iter2(cdbList);	  	 
  AliCDBId* id=0;
  while((id = dynamic_cast<AliCDBId*> (iter2.Next()))){	 
    cdbListCopy->Add(new TObjString(id->ToString().Data()));	 
    AliInfo(Form("Add %s to cdbList of ITSTPUserInfo",id->ToString().Data()));
  }	 
  // 
  uInfo->Add(cdbMapCopy);	 
  uInfo->Add(cdbListCopy);  
  //
  AliMagF *fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = fld ? fld->SolenoidField() : 0;
  TString bzString; bzString+=bz;
  TObjString *bzObjString = new TObjString(bzString);
  bzList = new TList();	 
  bzList->SetOwner(1);	 
  bzList->SetName("BzkGauss");	 
  bzList->Add(bzObjString);
  uInfo->Add(bzList);
  //
}

//_______________________________________________________________________________________
void AliAnalysisTaskITSAlignQA::CopyUserInfo()
{
  // if available, copy the UserInfo from the ESDtree to the output tree
  static Bool_t done = kFALSE;
  if (done) return;
  if (!fTPTree) {
    AliError("TrackPoints summary tree does not exist"); 
    return;
  }
  AliESDInputHandler *handler = (AliESDInputHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  TTree* esdTree = 0;
  if (!handler || !(esdTree=handler->GetTree())) return;
  if (esdTree->InheritsFrom(TChain::Class())) esdTree = esdTree->GetTree();
  TList* uInfoSrc = esdTree->GetUserInfo();
  const TMap  *cdbMapSrc  = (TMap*)uInfoSrc->FindObject("cdbMap");
  const TList *cdbListSrc = (TList*)uInfoSrc->FindObject("cdbList"); 
  if (!cdbMapSrc || !cdbListSrc) return;
  //
  AliInfo("Create ITSTPUserInfo from esdTree");
  TList* uInfoDst = fTPTree->GetUserInfo();
  //
  TMap *cdbMapCopy = new TMap(cdbMapSrc->GetEntries());	 
  cdbMapCopy->SetOwner(1);	 
  cdbMapCopy->SetName("cdbMap");	 
  TIter iter(cdbMapSrc->GetTable());	 
  TPair* pair = 0;	 
  while((pair = dynamic_cast<TPair*> (iter.Next()))){	 
    TObjString* keyStr = dynamic_cast<TObjString*> (pair->Key());	 
    TObjString* valStr = dynamic_cast<TObjString*> (pair->Value());
    if (keyStr && valStr) {
      cdbMapCopy->Add(new TObjString(keyStr->GetName()), new TObjString(valStr->GetName()));	 
      AliInfo(Form("Add %s : %s to cdbMap of ITSTPUserInfo",keyStr->GetName(),valStr->GetName()));
    }
  }	 
  //  
  TList *cdbListCopy = new TList();	 
  cdbListCopy->SetOwner(1);	 
  cdbListCopy->SetName("cdbList");	 
  // 
  TIter iter2(cdbListSrc);	  	 
  TObjString* id=0;
  while((id = dynamic_cast<TObjString*> (iter2.Next()))){	 
    cdbListCopy->Add(new TObjString(*id));	 
    AliInfo(Form("Add %s to cdbList of ITSTPUserInfo",id->GetName()));
  }	 
  // 
  uInfoDst->Add(cdbMapCopy);	 
  uInfoDst->Add(cdbListCopy);  
  //
  AliMagF *fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = fld ? fld->SolenoidField() : 0;
  TString bzString; bzString+=bz;
  TObjString *bzObjString = new TObjString(bzString);
  TList *bzList = new TList();	 
  bzList->SetOwner(1);	 
  bzList->SetName("BzkGauss");	 
  bzList->Add(bzObjString);
  uInfoDst->Add(bzList);
  //
  done = kTRUE;
}
