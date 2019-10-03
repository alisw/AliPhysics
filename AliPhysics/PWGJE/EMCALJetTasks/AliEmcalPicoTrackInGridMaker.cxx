// $Id$
//
// Class to put collection of tracks into grid of PicoTracks
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TH3F.h>

#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliJetContainer.h"

#include "AliEmcalPicoTrackInGridMaker.h"

ClassImp(AliEmcalPicoTrackInGridMaker)

//________________________________________________________________________
AliEmcalPicoTrackInGridMaker::AliEmcalPicoTrackInGridMaker() : 
AliAnalysisTaskEmcalJet("AliEmcalPicoTrackInGridMaker",kTRUE),
  fTracksOutName("PicoTracksInGrid"),
  fTracksOut(0),
  fL1Slide(0),
  fCellSize(0.0145),
  fMinCellE(0.15),
  fExclLeadingPatch(0),
  fPatchSub(3),
  fRhoMean(184.),
  fNCells(-1),
  fNCellsEMCal(-1),
  fNCellsDCal(-1),
  fMultVsRho(0)
{
  // Constructor.

  for(Int_t i = 0; i<2; i++) {
    fCellGrid[i] = 0;
    fMiniPatchGrid[i] = 0;
    fActiveAreaMP[i] = 0;
  }

  fPhiMin[0] = 1.405;
  fPhiMax[0] = 1.405+TMath::DegToRad()*110.;
  fPhiMin[1] = 4.547;
  fPhiMax[1] = 5.71;
  fEtaMin[0] = -0.7;
  fEtaMax[0] = 0.7;
  fEtaMin[1] = -0.7;
  fEtaMax[1] = 0.7;

  for(Int_t i = 0; i<5; i++) {
    fNPatchesEMCal[i] = 0;
    fh1RhoEmcal[i] = 0;
    fh1RhoDcal[i] = 0;
    fPatchEnVsActivityEmcal[i] = 0;
    fPatchEnVsActivityDcal[i]  = 0;

    for(Int_t j = 0; j<2; j++) {
      fPatchGrid[j][i] = 0;
      fActiveAreaMPP[j][i] = 0;
      fActiveAreaCP[j][i]  = 0;
      fPatchECorr[j][i] = 0;
      fPatchECorrPar[j][i] = 0;
      fPatchERaw[j][i] = 0;
      fPatchECorrRho[j][i] = 0;
      fPatchECorrECorrRho[j][i] = 0;
      fh2JetPtPatchECorr[j][i] = 0;
    }
    fh2PatchEtaPhiEmcal[i] = 0;
    fh2PatchEtaPhiDcal[i]  = 0;
  }

  for(Int_t i = 0; i<3; i++) {
    fh2MedianTypeEmcal[i] = 0;
    fh2MedianTypeDcal[i] = 0;
    fpMedianTypeEmcal[i] = 0;
    fpMedianTypeDcal[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliEmcalPicoTrackInGridMaker::AliEmcalPicoTrackInGridMaker(const char *name) : 
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fTracksOutName("PicoTracksInGrid"),
  fTracksOut(0),
  fL1Slide(0),
  fCellSize(0.0145),
  fMinCellE(0.15),
  fExclLeadingPatch(0),
  fPatchSub(3),
  fRhoMean(184.),
  fNCells(-1),
  fNCellsEMCal(-1),
  fNCellsDCal(-1),
  fMultVsRho(0)
{
  // Constructor.

  for(Int_t i = 0; i<2; i++) {
    fCellGrid[i] = 0;
    fMiniPatchGrid[i] = 0;
    fActiveAreaMP[i] = 0;
  }

  fPhiMin[0] = 1.405;
  fPhiMax[0] = 1.405+TMath::DegToRad()*110.;
  fPhiMin[1] = 4.547;
  fPhiMax[1] = 5.71;
  fEtaMin[0] = -0.7;
  fEtaMax[0] = 0.7;
  fEtaMin[1] = -0.7;
  fEtaMax[1] = 0.7;

  for(Int_t i = 0; i<5; i++) {
    fNPatchesEMCal[i] = 0;
    fh1RhoEmcal[i] = 0;
    fh1RhoDcal[i] = 0;
    fPatchEnVsActivityEmcal[i] = 0;
    fPatchEnVsActivityDcal[i]  = 0;

    for(Int_t j = 0; j<2; j++) {
      fPatchGrid[j][i] = 0;
      fActiveAreaMPP[j][i] = 0;
      fActiveAreaCP[j][i]  = 0;
      fPatchECorr[j][i] = 0;
      fPatchECorrPar[j][i] = 0;
      fPatchERaw[j][i] = 0;
      fPatchECorrRho[j][i] = 0;
      fPatchECorrECorrRho[j][i] = 0;
      fh2JetPtPatchECorr[j][i] = 0;
    }
    fh2PatchEtaPhiEmcal[i] = 0;
    fh2PatchEtaPhiDcal[i]  = 0;
  }

  for(Int_t i = 0; i<3; i++) {
    fh2MedianTypeEmcal[i] = 0;
    fh2MedianTypeDcal[i] = 0;
    fpMedianTypeEmcal[i] = 0;
    fpMedianTypeDcal[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliEmcalPicoTrackInGridMaker::~AliEmcalPicoTrackInGridMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::UserCreateOutputObjects()
{
  // Create my user objects.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Int_t nBinsMed = 500;
  Double_t minMed = 0.;
  Double_t maxMed = 500.;

  Int_t nBinsPhiEmcal = 132+64;
  Double_t phiMinEmcal = 1.436931 - 32.*fCellSize;
  Double_t phiMaxEmcal = 3.292931 + 32.*fCellSize;

  Int_t nBinsPhiDcal = 80+64;
  Double_t phiMinDcal = 4.664500 - 32.*fCellSize;
  Double_t phiMaxDcal = 5.592500 + 32.*fCellSize;

  Int_t nBinsEta = 96+64;
  Double_t etaMin = -0.696 - 32.*fCellSize;
  Double_t etaMax =  0.696 + 32.*fCellSize;

  for(Int_t i = 0; i<3; i++) {
    fh2MedianTypeEmcal[i] = new TH2F(Form("fh2MedianTypeEmcalAreaType%d",i),Form("fh2MedianTypeEmcalAreaType%d",i),5,0.5,5.5,nBinsMed,minMed,maxMed);
    fOutput->Add(fh2MedianTypeEmcal[i]);

    fh2MedianTypeDcal[i] = new TH2F(Form("fh2MedianTypeDcalAreaType%d",i),Form("fh2MedianTypeDcalAreaType%d",i),5,0.5,5.5,nBinsMed,minMed,maxMed);
    fOutput->Add(fh2MedianTypeDcal[i]);

    fpMedianTypeEmcal[i] = new TProfile(Form("fpMedianTypeEmcalAreaType%d",i),Form("fpMedianTypeEmcalAreaType%d",i),5,0.5,5.5,"s");
    fOutput->Add(fpMedianTypeEmcal[i]);

    fpMedianTypeDcal[i] = new TProfile(Form("fpMedianTypeDcalAreaType%d",i),Form("fpMedianTypeDcalAreaType%d",i),5,0.5,5.5,"s");
    fOutput->Add(fpMedianTypeDcal[i]);
  }

  TString det[2] = {"Emcal","Dcal"};
  for(Int_t i = 0; i<5; i++) { //loop over patch types
    fh1RhoEmcal[i] = new TH1F(Form("fh1RhoEmcal_%d",i),Form("fh1RhoEmcal_%d",i),500,0.,1000.);
    fOutput->Add(fh1RhoEmcal[i]);
    fh1RhoDcal[i] = new TH1F(Form("fh1RhoDcal_%d",i),Form("fh1RhoDcal_%d",i),500,0.,1000.);
    fOutput->Add(fh1RhoDcal[i]);

    fPatchEnVsActivityEmcal[i] = new TH2F(Form("fh2PatchEnVsActivityEmcal_%d",i),Form("fh2PatchEnVsActivityEmcal_%d",i),300,0.,300.,150,-0.5,149.5);
    fOutput->Add(fPatchEnVsActivityEmcal[i]);

    fPatchEnVsActivityDcal[i] = new TH2F(Form("fh2PatchEnVsActivityDcal_%d",i),Form("fh2PatchEnVsActivityDcal_%d",i),300,0.,300.,150,-0.5,149.5);
    fOutput->Add(fPatchEnVsActivityDcal[i]);

    for(Int_t j = 0; j<2; j++) {
      fPatchECorr[j][i] = new TH1F(Form("fPatchECorr%s_%d",det[j].Data(),i),Form("fPatchECorr%s_%d;#it{E}_{patch}^{corr}",det[j].Data(),i),250,-50.,200.);
      fOutput->Add(fPatchECorr[j][i]);

      fPatchECorrPar[j][i] = new TH1F(Form("fPatchECorrPar%s_%d",det[j].Data(),i),Form("fPatchECorrPar%s_%d;#it{E}_{patch}^{corr}",det[j].Data(),i),250,-50.,200.);
      fOutput->Add(fPatchECorrPar[j][i]);  

      fPatchERaw[j][i] = new TH1F(Form("fPatchERaw%s_%d",det[j].Data(),i),Form("fPatchERaw%s_%d;#it{E}_{patch}^{corr}",det[j].Data(),i),250,-50.,200.);
      fOutput->Add(fPatchERaw[j][i]);

      fPatchECorrRho[j][i] = new TH2F(Form("fPatchECorrRho%s_%d",det[j].Data(),i),Form("fPatchECorrRho%s_%d;#it{E}_{patch}^{corr};#rho",det[j].Data(),i),250,-50.,200.,500,0.,500.);
      fOutput->Add(fPatchECorrRho[j][i]);  

      fPatchECorrECorrRho[j][i] = new TH3F(Form("fPatchECorrECorrRho%s_%d",det[j].Data(),i),Form("fPatchECorrECorrRho%s_%d;#it{E}_{patch,det1}^{corr};#it{E}_{patch,det2}^{corr};#rho",det[j].Data(),i),210,-30.,180.,210,-30.,180.,250,0.,250.);
      fOutput->Add(fPatchECorrECorrRho[j][i]);

      fh2JetPtPatchECorr[j][i] = new TH2F(Form("fh2JetPtPatchECorr%s_%d",det[j].Data(),i),Form("fh2JetPtPatchECorr%s_%d",det[j].Data(),i),250,-50.,200.,250,-50.,200.);
      fOutput->Add(fh2JetPtPatchECorr[j][i]);
    }

    fh2PatchEtaPhiEmcal[i] = new TH2F(Form("fh2PatchEtaPhiEmcal_%d",i),Form("fh2PatchEtaPhiEmcal_%d;#eta;#phi",i),nBinsEta,etaMin,etaMax,nBinsPhiEmcal,phiMinEmcal,phiMaxEmcal);
    fOutput->Add(fh2PatchEtaPhiEmcal[i]);

    fh2PatchEtaPhiDcal[i] = new TH2F(Form("fh2PatchEtaPhiDcal_%d",i),Form("fh2PatchEtaPhiDcal_%d;#eta;#phi",i),nBinsEta,etaMin,etaMax,nBinsPhiDcal,phiMinDcal,phiMaxDcal);
    fOutput->Add(fh2PatchEtaPhiDcal[i]);
  }

  fMultVsRho = new TH2F("fMultVsRho","fMultVsRho",3000,0,3000,400,0,400);
  fOutput->Add(fMultVsRho);

  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::Run()
{
  // Main loop, called for each event.

  Bool_t b = CreateGridCells();
  if(!b) return kFALSE;
  b = CreateGridMiniPatches();
  if(!b) return kFALSE;

  //L0 single shower trigger
  CreateGridPatches(4,0);
  //  return kTRUE;

  //L1 triggers: sliding window
  CreateGridPatches(4,1);
  CreateGridPatches(8,1);
  CreateGridPatches(16,1);
  CreateGridPatches(32,1);


  Double_t medL0 = CalculateMedian(0,0);
  fh2MedianTypeEmcal[0]->Fill(0.5,medL0);
  fpMedianTypeEmcal[0]->Fill(0.5,medL0);
  medL0 = CalculateMedian(0,1);
  fh2MedianTypeDcal[0]->Fill(0.5,medL0);
  fpMedianTypeDcal[0]->Fill(0.5,medL0);

  Double_t medL1[4][2];
  for(Int_t i = 0; i<4; i++) { //patches
    for(Int_t type = 0; type<2; type++) { //EMCal or DCal
      for(Int_t areaT = 0; areaT<1; areaT++) { //areay type (passive vs active)
	medL1[i][type] = CalculateMedian(i+1,type,areaT);
	if(type==0) {
	  fh2MedianTypeEmcal[areaT]->Fill((Double_t)(i+1)+0.5,medL1[i][type]);
	  fpMedianTypeEmcal[areaT]->Fill((Double_t)(i+1)+0.5,medL1[i][type]);
	}
	if(type==1) {
	  fh2MedianTypeDcal[areaT]->Fill((Double_t)(i+1)+0.5,medL1[i][type]);
	  fpMedianTypeDcal[areaT]->Fill((Double_t)(i+1)+0.5,medL1[i][type]);
	}
      }
    }
  }

  // subtract energy density and store energy distributions of corrected patches in histo
  Int_t EleadID[5][2];
  Double_t Elead[5][2];
  Double_t EleadRaw[5][2];
  for(Int_t i = 1; i<5; i++) { //patch types
    Int_t stepSize = GetTriggerPatchIdStepSizeNoOverlap(GetPatchDim(i));
    for(Int_t type = 0; type<2; type++) {
      EleadID[i][type] = -1;
      Elead[i][type] = -1e6;
      EleadRaw[i][type] = -1e6;
      Int_t subType = 1;
      if(type==1) subType = 0;
      //      for(Int_t j = 0; j<(fPatchGrid[type][i].GetSize()-stepSize+1); j+=stepSize) { //patches
      for(Int_t k = 0; k<GetNColTriggerPatches(type,GetPatchDim(i),i); k+=stepSize) {
	for(Int_t l = 0; l<GetNRowTriggerPatches(type,GetPatchDim(i),i); l+=stepSize) {
	  Int_t id = GetTriggerPatchID(k,l,type,GetPatchDim(i),i);
	  //	  if(type==1 && i==4) Printf("id: %d/%d k: %d/%d l: %d/%d",id,fPatchGrid[type][i].GetSize(),k,GetNColTriggerPatches(type,GetPatchDim(i),i),l,GetNRowTriggerPatches(type,GetPatchDim(i),i));
	  if(fPatchGrid[type][i].At(id)>0.) { //don't do anything with empty patches
	    Double_t sub = medL1[fPatchSub-1][subType]*GetPatchArea(i);
	    fPatchECorr[type][i]->Fill(fPatchGrid[type][i].At(id) - sub);
	    fPatchECorrPar[type][i]->Fill(fPatchGrid[type][i].At(id) - fRhoMean*GetPatchArea(i));
	
	    //Bookkeep leading patches
	    if((fPatchGrid[type][i].At(id)-sub)>Elead[i][type]) {
	      EleadID[i][type] = id;
	      Elead[i][type] = fPatchGrid[type][i].At(id)-sub;
	    }
	    if(fPatchGrid[type][i].At(id)>EleadRaw[i][type])
	      EleadRaw[i][type] = fPatchGrid[type][i].At(id);
	  }
	}//cols
      }//rows
    }//type

    AliJetContainer *cont = GetJetContainer(0);
    for(Int_t k = 0; k<2; k++) {
      Int_t subType = 1;
      if(k==1) subType=0;
      fPatchECorrRho[k][i]->Fill(Elead[i][k],medL1[fPatchSub-1][subType]);
      fPatchECorrECorrRho[k][i]->Fill(Elead[i][k],Elead[i][subType],medL1[fPatchSub-1][subType]);
      fPatchERaw[k][i]->Fill(EleadRaw[i][k]);
      //Save jet spectra within EMCal and DCal fiducial acceptance
      //jet pT vs highest energy patch for preferred trigger patch types
      if(cont) {
	Double_t r = cont->GetJetRadius();
	cont->SetJetEtaLimits(fEtaMin[k]+r,fEtaMax[k]-r);
	cont->SetJetPhiLimits(fPhiMin[k]+r,fPhiMax[k]-r);
	Double_t rho = cont->GetRhoVal();
	AliEmcalJet *jet = NULL;
	cont->ResetCurrentID();
	while((jet = cont->GetNextAcceptJet())) {
	  Double_t jetpt = jet->Pt() - rho*jet->Area();
	  fh2JetPtPatchECorr[k][i]->Fill(jetpt,Elead[i][k]);
	}
      }
    }
    Double_t eta = 0.; Double_t phi = 0.;
    GetEtaPhiFromTriggerPatchID(EleadID[i][0],0,GetPatchDim(i),1,eta,phi);
    fh2PatchEtaPhiEmcal[i]->Fill(eta,phi);

    GetEtaPhiFromTriggerPatchID(EleadID[i][1],1,GetPatchDim(i),1,eta,phi);
    fh2PatchEtaPhiDcal[i]->Fill(eta,phi);
  } //patch types

  fMultVsRho->Fill(GetParticleContainer(0)->GetNParticles(),medL1[3][0]);
  return kTRUE;
}

//________________________________________________________________________
AliEmcalJet* AliEmcalPicoTrackInGridMaker::GetClosestJet(const Double_t eta, const Double_t phi, const Int_t icont) const {

  AliJetContainer *cont = GetJetContainer(icont);
  if(!cont) return NULL;
  
  Double_t closest_dr = 1e6;
  Int_t closest_id = -1;
  AliEmcalJet *jet = NULL;
  cont->ResetCurrentID();
  while((jet = cont->GetNextAcceptJet())) {
    Double_t dphi = jet->Phi() - phi;
    Double_t deta = jet->Eta() - eta;
    dphi = TVector2::Phi_mpi_pi(dphi);
    Double_t dr = TMath::Sqrt ( dphi * dphi + deta * deta );
    if(dr<closest_dr) {
      closest_dr = dr;
      closest_id = cont->GetCurrentID();
    }
  }
  jet = cont->GetJet(closest_id);
  return jet;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::CalculateSum(const Int_t patchType) const {
  //calc total energy of all patches
  Int_t n = fPatchGrid[0][patchType].GetSize();
  if(n<1) return -1.;

  Double_t sum = 0.;
  Int_t count = 0;
  Int_t stepSize = GetTriggerPatchIdStepSizeNoOverlap(GetPatchDim(patchType));
  for(Int_t type = 0; type<2; type++) {
    for(Int_t i = 0; i<fPatchGrid[type][patchType].GetSize(); i+=stepSize) {
      if(fPatchGrid[type][patchType].At(i)>0.) count++;
      sum+=fPatchGrid[type][patchType].At(i);
    }
  }
  return sum;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::CalculateMedian(const Int_t patchType, const Int_t type, const Int_t areaType) {
  //areaType:
  //0: passive area
  //1: active area mini patches
  //2: active arear cells

  Int_t n = fPatchGrid[type][patchType].GetSize();
  if(n<1) return -1.;

  Int_t level = 0;
  if(patchType>0) level = 1;
  Int_t dim = GetPatchDim(patchType);
  Double_t area = GetPatchArea(patchType);
  Int_t stepSize = GetTriggerPatchIdStepSizeNoOverlap(GetPatchDim(patchType),level);
  //Printf("patchType: %d dim: %d stepSizeNoOverlap: %d ",patchType,GetPatchDim(patchType),stepSize);

  static Double_t arr[999];
  Int_t c = 0;

  //find patch with highest energy
  Int_t imax = -1;
  Double_t max = 0.;
  for(Int_t i = 0; i<(GetNColTriggerPatches(type,dim,patchType)); i+=stepSize) {
    for(Int_t j = 0; j<(GetNRowTriggerPatches(type,dim,patchType)); j+=stepSize) {
      Int_t id = GetTriggerPatchID(i,j,type,dim,patchType);
      //      Printf("id: %d/%d i: %d/%d j:%d/%d",id,fPatchGrid[type][patchType].GetSize(),i,GetNColTriggerPatches(type,dim,patchType),j,GetNRowTriggerPatches(type,dim,patchType));
      if(fPatchGrid[type][patchType].At(id)>max) {
	imax = id;
	max = fPatchGrid[type][patchType].At(id);
      }
    }//cols
  }//rows

  for(Int_t i = 0; i<GetNColTriggerPatches(type,dim,patchType); i+=stepSize) {
    for(Int_t j = 0; j<GetNRowTriggerPatches(type,dim,patchType); j+=stepSize) {
      Int_t id = GetTriggerPatchID(i,j,type,dim,patchType);
      if(fExclLeadingPatch>0 && id==imax) continue;
      if(fPatchGrid[type][patchType].At(id)>0.) {
	Int_t active = 99;
	if(areaType==1) active = fActiveAreaMPP[type][patchType].At(id);
	else if(areaType==2) active = fActiveAreaCP[type][patchType].At(id);
	if(areaType>0) area = GetPatchAreaActive(id,type,patchType,areaType-1);

	if(area>0. && active>1) {
	  arr[c] = fPatchGrid[type][patchType].At(id)/area;
	  c++;
	}
	if(type==0 && areaType==0) {
	  fh1RhoEmcal[patchType]->Fill(arr[c-1]);
	  fPatchEnVsActivityEmcal[patchType]->Fill(fPatchGrid[type][patchType].At(id),fActiveAreaCP[type][patchType].At(id));
	}
	if(type==1 && areaType==0) {
	  fh1RhoDcal[patchType]->Fill(arr[c-1]);
	  fPatchEnVsActivityDcal[patchType]->Fill(fPatchGrid[type][patchType].At(id),fActiveAreaCP[type][patchType].At(id));
	}
      }
    }
  }
  Double_t med = TMath::Median(c,arr);
  return med;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::CreateGridCells() {
  //create cells from track input
  if(!InitCells()) return kFALSE;

  AliVParticle *track = NULL;
  AliParticleContainer *trackCont = GetParticleContainer(0);
  if(!trackCont) return kFALSE;
  trackCont->ResetCurrentID();
  while((track = trackCont->GetNextAcceptParticle())) {
    if(track->Pt()<fMinCellE) continue;
    Int_t id = GetGridID(track);
    Int_t type = GetCellType(track);
    if(id>-1)
      fCellGrid[type].AddAt(fCellGrid[type].At(id)+track->Pt(),id);
    }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::CreateGridMiniPatches() {
  //create mini patches (2x2 cells)
  if(!InitMiniPatches()) return kFALSE;

  //loop over edges of mini patches
  Int_t nm = 0; //mini patch number
  for(Int_t type = 0; type<2; type++) {
    nm = 0;
    for(Int_t i = 0; i<(GetNCellsCol(type)-1); i+=2) {
      for(Int_t j = 0; j<(GetNCellsRow(type)-1); j+=2) {
	//loop over cells in mini patch
	for(Int_t k = 0; k<2; k++) {
	  for(Int_t l = 0; l<2; l++) {
	    Int_t id = GetGridID(i+k,j+l,type);
	    fMiniPatchGrid[type].AddAt(fMiniPatchGrid[type].At(nm)+fCellGrid[type].At(id),nm);
	    if(fCellGrid[type].At(id)>0.)
	      fActiveAreaMP[type].AddAt(fActiveAreaMP[type].At(nm)+1,nm);
	  }
	}
	nm++;
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNRowMiniPatches(const Int_t type) const {
  //returns number of rows of mini patches in detector of type X (0: EMCal 1: DCal)
  Int_t nRows = TMath::FloorNint(0.5*GetNCellsCol(type));
  return nRows;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNColMiniPatches(const Int_t type) const {
  //returns number of rows of mini patches in detector of type X (0: EMCal 1: DCal)
  Int_t nCols = TMath::FloorNint(0.5*GetNCellsRow(type));
  return nCols;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::CreateGridPatches(const Int_t dim, const Int_t level) {
  //create trigger patches
  if(!InitPatches(dim,level)) return kFALSE;

  Int_t pt = GetPatchType(dim,level);
  if(pt<0) return kFALSE;
  Int_t nm = (Int_t)(dim/2.);    //size of trigger patch in number of mini patches
  Int_t stepm = (Int_t)(dim/2.); //step size through grid in mini patches
  if(level==1 && fL1Slide)  stepm = GetSlidingStepSizeMiniPatches(dim,level);
  //loop over edges of mini patches
  for(Int_t type = 0; type<2; type++) {
    Int_t np = 0; //patch number
      for(Int_t j = 0; j<=(GetNColMiniPatches(type)-nm); j+=stepm) {
    for(Int_t i = 0; i<=(GetNRowMiniPatches(type)-nm); i+=stepm) {
      //      for(Int_t j = 0; j<=(GetNColMiniPatches(type)-nm); j+=stepm) {
	//loop over mini patches in patch
	for(Int_t k = 0; k<nm; k++) {
	  for(Int_t l = 0; l<nm; l++) {
	    Int_t row = i+k;
	    Int_t col = j+l;
	    Int_t id = GetMiniPatchID(row,col,type);
	    fPatchGrid[type][pt].AddAt(fPatchGrid[type][pt].At(np)+fMiniPatchGrid[type].At(id),np);
	    if(fMiniPatchGrid[type].At(id)>0.) {
	      fActiveAreaMPP[type][pt].AddAt(fActiveAreaMPP[type][pt].At(np)+1,np);
	      fActiveAreaCP[type][pt].AddAt(fActiveAreaCP[type][pt].At(np)+fActiveAreaMP[type].At(id),np);
	    }
	  }
	}
	np++;
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetTriggerPatchID(const Int_t row, const Int_t col, const Int_t type, const Int_t dim, const Int_t patchType) const {
  Int_t id = row*GetNRowTriggerPatches(type,dim,patchType) + col;
  return id;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetMiniPatchID(const Int_t row, const Int_t col, const Int_t type) const {
  Int_t id = row*GetNColMiniPatches(type) + col;
  return id;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetCellType(const Double_t eta, const Double_t phi) const {
  //cell in EMCal (0) or DCal (1)
  for(Int_t i = 0; i<2; i++) {
    if(eta>fEtaMin[i] && eta<fEtaMax[i] && phi>fPhiMin[i] && phi<fPhiMax[i]) return i;
  }
  return -1;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetGridID(const Int_t row, const Int_t col, const Int_t type) const {
  Int_t id = row*GetNCellsRow(type) + col;
  return id;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetGridID(const Double_t eta, const Double_t phi) const {
  
  Int_t type = GetCellType(eta,phi);
  if(type<0 || type>1) return -1; //position is not in EMCal or DCal

  // grid ID convention:
  // upper left corner (min phi, min eta) is first ID
  // then walk through grid from upper left to lower right accross the rows in phi

  Int_t id = -1;
  Double_t etaRel = eta-fEtaMin[type];
  Double_t phiRel = phi-fPhiMin[type];
  Int_t row = TMath::FloorNint(etaRel/fCellSize);
  Int_t col = TMath::FloorNint(phiRel/fCellSize);
  id = GetGridID(row,col,type);

  if(id>=fNCells) {
    Printf("Got too large id %d %d type: %d",id,fNCells,type);
    Printf("eta: %f phi: %f",eta,phi);
    Printf("etaRel: %f phiRel: %f",etaRel,phiRel);
    Printf("row: %d col: %d  ->  %d + %d = %d",row,col,row*GetNCellsRow(type) + col,fNCellsEMCal,id);
    Printf("n cells row: %d",GetNCellsRow(type));
    Printf("\n");
  }
  return id;
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::GetEtaPhiFromGridID(const Int_t id, const Int_t type, Double_t &eta, Double_t &phi) const {
  //returns eta phi of cell at lower right edge (lowest eta, lowest phi)
  Int_t row = TMath::FloorNint(id/GetNCellsRow(type));
  Int_t col = id - row * GetNCellsRow(type);
  eta = fEtaMin[type] + row*fCellSize;
  phi = fPhiMin[type] + col*fCellSize;
  AliDebug(2,Form("id: %d type: %d row: %d col: %d eta: %f phi: %f",id,type,row,col,eta,phi));
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::GetEtaPhiFromMiniPatchID(const Int_t id, const Int_t type, Double_t &eta, Double_t &phi) const {
  //returns eta phi of mini patch at lower right edge (lowest eta, lowest phi)
  Int_t row = TMath::FloorNint(id/GetNColMiniPatches(type));
  Int_t col = id - row * GetNColMiniPatches(type);
  eta = fEtaMin[type] + row*2.*fCellSize;
  phi = fPhiMin[type] + col*2.*fCellSize;
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::GetEtaPhiFromTriggerPatchID(const Int_t id, const Int_t type, const Int_t dim, const Int_t level, Double_t &eta, Double_t &phi) const {
  //returns eta phi of mini patch at lower right edge (lowest eta, lowest phi)
  Int_t step = GetSlidingStepSizeCells(dim,level); //id: 8/96 k(row): 0/8 l(col): 8/12
  //  Int_t offset = dim/step;
  Int_t row = TMath::FloorNint(id/GetNColTriggerPatches(type,dim,GetPatchType(dim,level)));
  Int_t col = id - row * GetNColTriggerPatches(type,dim,GetPatchType(dim,level));
  eta = fEtaMin[type] + row*step*fCellSize;
  phi = fPhiMin[type] + col*step*fCellSize;
  // if(dim==32) {
  //   Printf("dim: %d step: %d offset: %d",dim,step,offset);
  //   Printf("dim: %d id: %d row: %d col: %d eta: %f phi: %f",dim,id,row,col,eta,phi);
  // }
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNColTriggerPatches(const Int_t type, const Int_t dim, const Int_t patchType) const {
  //returns number of trigger patch columns
  Int_t level = 0;
  if(patchType>0) level = 1;
  Double_t stepmp = (Double_t)GetSlidingStepSizeMiniPatches(dim,level);
  Int_t nmp = GetNColMiniPatches(type);
  Int_t ntc = TMath::FloorNint(nmp/stepmp);
  //  Printf("dim: %d stepmp: %f nmp: %d ntc: %d",dim,stepmp,nmp,ntc);
  return ntc;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNRowTriggerPatches(const Int_t type, const Int_t dim, const Int_t patchType) const {
  //returns number of trigger patch rows
  Int_t level = 0;
  if(patchType>0) level = 1;
  Double_t stepmp = (Double_t)GetSlidingStepSizeMiniPatches(dim,level);
  Int_t nmp = GetNRowMiniPatches(type);
  Int_t ntr = TMath::FloorNint(nmp/stepmp);
  return ntr;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNCellsCol(const Int_t type) const {
  //returns number of cells in column
  Double_t deta = fEtaMax[type] - fEtaMin[type];
  Int_t nCellsCol = TMath::FloorNint(deta/fCellSize);
  return nCellsCol;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNCellsRow(const Int_t type) const {
  //returns number of cells in row
  Double_t dPhi = fPhiMax[type] - fPhiMin[type];
  Int_t nCellsRow = TMath::FloorNint(dPhi/fCellSize);
  return nCellsRow;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::InitCells() {
  //initialize cells array
  CheckEdges();
  if(!CheckEdges()) return kFALSE;

  //number of cells in EMCal acceptance
  Int_t nCellsPhiE = GetNCellsCol(0);
  Int_t nCellsEtaE = GetNCellsRow(0);
  fNCellsEMCal = nCellsPhiE*nCellsEtaE;
  
  //number of cells in DCal acceptance
  Int_t nCellsPhiD = GetNCellsCol(1);
  Int_t nCellsEtaD = GetNCellsRow(1);
  fNCellsDCal = nCellsPhiD*nCellsEtaD;

  //total number of cells
  fNCells = fNCellsEMCal + fNCellsDCal;
  
  AliDebug(2,Form("EMCal: %d x %d",nCellsEtaE,nCellsPhiE));
  AliDebug(2,Form("DCal: row: %d x col: %d",nCellsEtaD,nCellsPhiD));
  AliDebug(2,Form("fNCells: %d fNCellsE: %d fnCellsD: %d",fNCells,fNCellsEMCal,fNCellsDCal));
  
  fCellGrid[0].Set(fNCellsEMCal);
  fCellGrid[1].Set(fNCellsDCal);
  fCellGrid[0].Reset(0);
  fCellGrid[1].Reset(0);
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::InitMiniPatches() {
  //initialize mini patch array
  if(fCellGrid[0].GetSize()<0) return kFALSE;
  if(fCellGrid[1].GetSize()<0) return kFALSE;
  Double_t conv = 0.25; //dimension of mini patch is 2x2 cells
  Int_t nMiniPatches[2];
  nMiniPatches[0] = (Int_t)(conv*fNCellsEMCal);
  nMiniPatches[1] = (Int_t)(conv*fNCellsDCal);
  for(Int_t i = 0; i<2; i++) {
    fMiniPatchGrid[i].Set(nMiniPatches[i]);
    fMiniPatchGrid[i].Reset(0.);

    fActiveAreaMP[i].Set(nMiniPatches[i]);
    fActiveAreaMP[i].Reset(0);
  }
  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNTriggerPatches(const Int_t type, const Int_t dim, const Int_t level) const {
  //get number of trigger patches in EMCal or DCal
  Double_t dimd = (Double_t)dim;
  Double_t conv = 1./(dimd*dimd);
  if(level==1) {
    Double_t step = (Double_t)GetSlidingStepSizeCells(dim);
    conv = 1./(step*step);
  }
  Int_t nPatches = 0;
  if(type==0) nPatches = (Int_t)(conv*fNCellsEMCal);
  else if(type==1) nPatches = (Int_t)(conv*fNCellsDCal);
  return nPatches;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::InitPatches(const Int_t dim, const Int_t level) {
  // dimensions in cell units
  // if level==1: sliding window will be applied
  // L1 4x4: slide by 2 cells: 1 mini patch
  // L1 8x8: slide by 4 cells: 2 mini patches
  // L1 16x16: slide by 4 cells: 2 mini patches
  // L1 32x32: slide by 8 cells: 4 mini patches
  
  if(fCellGrid[0].GetSize()<0) return kFALSE;
  if(fCellGrid[1].GetSize()<0) return kFALSE;

  Int_t type = GetPatchType(dim,level);
  if(type<0 || type>4) return kFALSE;

  Int_t nPatches[2]; //number of trigger patches in EMCal and DCal
  for(Int_t i = 0; i<2; i++)
    nPatches[i] = GetNTriggerPatches(i,dim,level);
  //total number of trigger patches
  Int_t nPatchesT = nPatches[0] + nPatches[1];

  fNPatchesEMCal[type] = nPatches[0];
  AliDebug(2,Form("Create trigger patch of type %d with dim %d and with %d patches EMCAL: %d DCAL: %d",type,dim,nPatchesT,nPatches[0],nPatches[1]));
  //Printf("Create trigger patch of type %d with dim %d and with %d patches EMCAL: %d DCAL: %d",type,dim,nPatchesT,nPatches[0],nPatches[1]);
  for(Int_t i = 0; i<2; i++) {
    fPatchGrid[i][type].Set(nPatches[i]);
    fPatchGrid[i][type].Reset(0);
    fActiveAreaMPP[i][type].Set(nPatches[i]);
    fActiveAreaMPP[i][type].Reset(0);
    fActiveAreaCP[i][type].Set(nPatches[i]);
    fActiveAreaCP[i][type].Reset(0);
  }

  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetPatchDim(const Int_t ipatch) const {
  //returns total area of patch
  Int_t ncell = 4;
  if(ipatch==0) ncell = 4;
  if(ipatch==1) ncell = 4;
  if(ipatch==2) ncell = 8;
  if(ipatch==3) ncell = 16;
  if(ipatch==4) ncell = 32;

  return ncell;
}
 
//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::GetPatchArea(const Int_t ipatch) const {
  //returns total area of patch
  Double_t ncell = (Double_t)GetPatchDim(ipatch);
  Double_t area = ncell*ncell*fCellSize*fCellSize;
  return area;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::GetPatchAreaActive(const Int_t id, const Int_t type, const Int_t ipatch, const Int_t atype) const {
  //atype = 0 : active area from mini patches
  //atype = 1 : active area from cells

  Int_t active = 0;
  if(atype==0) active = fActiveAreaMPP[type][ipatch].At(id);
  else if(atype==1) active = fActiveAreaCP[type][ipatch].At(id);
  else return -1;

  Double_t fac = 1.;
  if(atype==0) fac = 4.;
  Double_t area = active*fac*fCellSize*fCellSize;
  return area;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetSlidingStepSizeCells(const Int_t dim, const Int_t level) const {
  //get step size for mock-up L1 trigger
  if(!fL1Slide || level==0) return dim;

  if(dim==4) return 2;
  else if(dim==8) return 4;
  else if(dim==16) return 4;
  else if(dim==32) return 8;
  else return -1;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetSlidingStepSizeMiniPatches(const Int_t dim, const Int_t level) const {
  //returns step size in mini patches
  Double_t step = (Double_t)GetSlidingStepSizeCells(dim,level);
  if(step<0) return step;
  Int_t stepMiniPatch = (Int_t)(step/2.);
  return stepMiniPatch;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetTriggerPatchIdStepSizeNoOverlap(const Int_t dim, const Int_t level) const {
  //return step for trigger patch id's without overlapping patches
  if(!fL1Slide) return 1;
  
  Int_t cellStep = GetSlidingStepSizeCells(dim,level);
  Int_t step = TMath::FloorNint((Double_t)(dim/cellStep));
  return step;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetPatchType(const Int_t dim, const Int_t level) const {
  //type of trigger patch (size)
  Int_t type = -1;
  if(level==0 && dim==4) type = 0;
  if(level==1) {
    if(dim==4)  type = 1;
    if(dim==8)  type = 2;
    if(dim==16) type = 3;
    if(dim==32) type = 4;
  }
  return type;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::CheckEdges() {
  //Check if defined edges of EMCal and DCal make sense
  if(fPhiMin[0]<0. || fPhiMax[0]<fPhiMin[0]) {
    AliDebug(11,Form("EMCal phi edges not defined %f-%f",fPhiMin[0],fPhiMax[0]));
    return kFALSE;
  }

  if(fPhiMin[1]<0. || fPhiMax[1]<fPhiMin[1]) {
    AliDebug(11,Form("DCal phi edges not defined %f-%f",fPhiMin[1],fPhiMax[1]));
    return kFALSE;
  }

  if(fEtaMin[0]<-10. || fEtaMax[0]<fEtaMin[1]) {
    AliDebug(11,Form("EMCal eta edges not well defined %f-%f",fEtaMin[0],fEtaMax[0]));
    return kFALSE;
  }

  if(fEtaMin[1]<-10. || fEtaMax[1]<fEtaMin[1]) {
    AliDebug(11,Form("DCal eta edges not well defined %f-%f",fEtaMin[1],fEtaMax[1]));
    return kFALSE;
  }

  for(Int_t type = 0; type<2; type++) {
    Double_t dphi = fPhiMax[type] - fPhiMin[type];
    Double_t nPatchPhi32 = TMath::Floor(dphi/(32.*fCellSize));
    Double_t nCellsColExact = nPatchPhi32 * 32. ;
    
    Double_t deta = fEtaMax[type] - fEtaMin[type];
    Double_t nPatchEta32 = TMath::Floor(deta/(32.*fCellSize));
    Double_t nCellsRowExact = nPatchEta32 * 32. ;

    Double_t col_extra = dphi/fCellSize - TMath::Floor(nCellsColExact);
    Double_t phi_extra = col_extra*fCellSize;
    fPhiMin[type] += phi_extra*0.5;
    fPhiMax[type] -= phi_extra*0.5;

    Double_t row_extra = deta/fCellSize - TMath::Floor(nCellsRowExact);
    Double_t eta_extra = row_extra*fCellSize;
    fEtaMin[type] += eta_extra*0.5;
    fEtaMax[type] -= eta_extra*0.5;
    AliDebug(2,Form("type: %d  exact: col: %f row: %f",type,nCellsColExact,nCellsRowExact));
    //    Printf("type: %d  exact: col: %f row: %f",type,nCellsColExact,nCellsRowExact);
    //  PrintAcceptance();
  }
  return kTRUE;
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::PrintAcceptance() const {
  Printf("EMCal");
  Printf("phi: %f-%f",fPhiMin[0],fPhiMax[0]);
  Printf("eta: %f-%f",fEtaMin[0],fEtaMax[0]);

  Printf("DCal");
  Printf("phi: %f-%f",fPhiMin[1],fPhiMax[1]);
  Printf("eta: %f-%f",fEtaMin[1],fEtaMax[1]);
}
