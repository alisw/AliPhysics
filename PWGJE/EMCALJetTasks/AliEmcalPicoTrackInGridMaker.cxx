// $Id$
//
// Class to put collection of tracks into grid of PicoTracks
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TH3F.h>
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
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
  fCellGrid(),
  fMiniPatchGrid(),
  fActiveAreaMP(),
  fMultVsRho(0)
{
  // Constructor.

  fPhiMin[0] = 1.405;
  fPhiMax[0] = 1.405+TMath::DegToRad()*110.;
  fPhiMin[1] = 4.547;
  fPhiMax[1] = 5.71;
  fEtaMin[0] = -0.7;
  fEtaMax[0] = 0.7;
  fEtaMin[1] = -0.7;
  fEtaMax[1] = 0.7;

  for(Int_t i = 0; i<5; i++) {
    fPatchGrid[i] = 0;
    fNPatchesEMCal[i] = 0;
    fActiveAreaMPP[i] = 0;
    fActiveAreaCP[i]  = 0;

    fh1RhoEmcal[i] = 0;
    fh1RhoDcal[i] = 0;
    fPatchEnVsActivityEmcal[i] = 0;
    fPatchEnVsActivityDcal[i]  = 0;

    for(Int_t j = 0; j<2; j++) {
      fPatchECorr[j][i] = 0;
      fPatchECorrPar[j][i] = 0;
      fPatchECorrRho[j][i] = 0;
      fPatchECorrRhoDijet[j][i] = 0;
      fPatchECorrECorrRho[j][i] = 0;
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
  fCellGrid(),
  fMiniPatchGrid(),
  fActiveAreaMP(),
  fMultVsRho(0)
{
  // Constructor.

  fPhiMin[0] = 1.405;
  fPhiMax[0] = 1.405+TMath::DegToRad()*110.;
  fPhiMin[1] = 4.547;
  fPhiMax[1] = 5.71;
  fEtaMin[0] = -0.7;
  fEtaMax[0] = 0.7;
  fEtaMin[1] = -0.7;
  fEtaMax[1] = 0.7;

  for(Int_t i = 0; i<5; i++) {
    fPatchGrid[i] = 0;
    fNPatchesEMCal[i] = 0;
    fActiveAreaMPP[i] = 0;
    fActiveAreaCP[i]  = 0;

    fh1RhoEmcal[i] = 0;
    fh1RhoDcal[i] = 0;

    fPatchEnVsActivityEmcal[i] = 0;
    fPatchEnVsActivityDcal[i]  = 0;

    for(Int_t j = 0; j<2; j++) {
      fPatchECorr[j][i] = 0;
      fPatchECorrPar[j][i] = 0;
      fPatchECorrRho[j][i] = 0;
      fPatchECorrRhoDijet[j][i] = 0;
      fPatchECorrECorrRho[j][i] = 0;
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

      fPatchECorrRho[j][i] = new TH2F(Form("fPatchECorrRho%s_%d",det[j].Data(),i),Form("fPatchECorrRho%s_%d;#it{E}_{patch}^{corr};#rho",det[j].Data(),i),250,-50.,200.,500,0.,500.);
      fOutput->Add(fPatchECorrRho[j][i]);  

      fPatchECorrRhoDijet[j][i] = new TH2F(Form("fPatchECorrRhoDijet%s_%d",det[j].Data(),i),Form("fPatchECorrRhoDijet%s_%d;#it{E}_{patch}^{corr};#rho",det[j].Data(),i),250,-50.,200.,500,0.,500.);
      fOutput->Add(fPatchECorrRhoDijet[j][i]);

      fPatchECorrECorrRho[j][i] = new TH3F(Form("fPatchECorrECorrRho%s_%d",det[j].Data(),i),Form("fPatchECorrECorrRho%s_%d;#it{E}_{patch,det1}^{corr};#it{E}_{patch,det2}^{corr};#rho",det[j].Data(),i),210,-30.,180.,210,-30.,180.,250,0.,250.);
      fOutput->Add(fPatchECorrECorrRho[j][i]);
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
  for(Int_t i = 1; i<5; i++) { //patch types
    Int_t EleadID[2] = {-1,-1};
    Double_t Elead[2] = {0.,0.};
    for(Int_t j = 0; j<fPatchGrid[i].GetSize(); j++) { //patches
      if(fPatchGrid[i].At(j)>0.) { //don't do anything with empty patches
	Int_t type = 0; //EMCal
	Int_t subType = 1;
	if(j>=fNPatchesEMCal[i]) {
	  type = 1; //DCal
	  subType = 0;
	}
	Double_t sub = medL1[fPatchSub-1][subType]*GetPatchArea(i);
	fPatchECorr[type][i]->Fill(fPatchGrid[i].At(j) - sub);
	fPatchECorrPar[type][i]->Fill(fPatchGrid[i].At(j) - fRhoMean*GetPatchArea(i));
	
	//Bookkeep leading patches
	if((fPatchGrid[i].At(j)-sub)>Elead[type]) {
	  EleadID[type] = j;
	  Elead[type] = fPatchGrid[i].At(j)-sub;
	}
      }
    }

    for(Int_t k = 0; k<2; k++) {
      Int_t subType = 1;
      if(k==1) subType=0;
      fPatchECorrRho[k][i]->Fill(Elead[k],medL1[fPatchSub-1][subType]);
      fPatchECorrECorrRho[k][i]->Fill(Elead[k],Elead[subType],medL1[fPatchSub-1][subType]);
      if(Elead[subType]>30.) {
    	fPatchECorrRhoDijet[k][i]->Fill(Elead[k],medL1[fPatchSub-1][subType]);
      }
    }
    Double_t eta = 0.; Double_t phi = 0.;
    GetEtaPhiFromTriggerPatchID(EleadID[0],GetPatchDim(i),1,eta,phi);
    fh2PatchEtaPhiEmcal[i]->Fill(eta,phi);

    GetEtaPhiFromTriggerPatchID(EleadID[1],GetPatchDim(i),1,eta,phi);
    fh2PatchEtaPhiDcal[i]->Fill(eta,phi);
  }

  fMultVsRho->Fill(GetParticleContainer(0)->GetNParticles(),medL1[3][0]);
  return kTRUE;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::CalculateSum(const Int_t patchType) const {
  //calc total energy of all patches
  Int_t n = fPatchGrid[patchType].GetSize();
  if(n<1) return -1.;

  Double_t sum = 0.;
  Int_t count = 0;
  for(Int_t i = 0; i<fPatchGrid[patchType].GetSize(); i++) {
    if(fPatchGrid[patchType].At(i)>0.) count++;
    sum+=fPatchGrid[patchType].At(i);
  }
  return sum;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::CalculateMedian(const Int_t patchType, const Int_t type, const Int_t areaType) {
  //areaType:
  //0: passive area
  //1: active area mini patches
  //2: active arear cells

  Int_t n = fPatchGrid[patchType].GetSize();
  if(n<1) return -1.;
  //  Double_t *arr = fPatchGrid[patchType].GetArray();
  
  static Double_t arr[999];
  Int_t c = 0;
  Int_t start = 0;
  Int_t last = n;
  if(type==0) {
    start = 0;
    last = fNPatchesEMCal[patchType];
  } else if(type==1) {
    start = fNPatchesEMCal[patchType];
    last = n;
  }

  Double_t area = GetPatchArea(patchType);

  //find patch with highest energy
  Int_t imax = -1;
  Double_t max = 0.;
  for(Int_t i = start; i<last; i++) {
    if(fPatchGrid[patchType].At(i)>max) {
      imax = i;
      max = fPatchGrid[patchType].At(i);
    }
  }
  for(Int_t i = start; i<last; i++) {
    if(fExclLeadingPatch>0 && i==imax) continue;
    if(fPatchGrid[patchType].At(i)>0.) {
      Int_t active = 99;
      if(areaType==1) active = fActiveAreaMPP[patchType].At(i);
      else if(areaType==2) active = fActiveAreaCP[patchType].At(i);
      if(areaType>0) area = GetPatchAreaActive(i,patchType,areaType-1);

      if(area>0. && active>1) {
	arr[c] = fPatchGrid[patchType].At(i)/area;
	c++;
      }
      if(type==0 && areaType==0) {
	fh1RhoEmcal[patchType]->Fill(arr[c-1]);
	fPatchEnVsActivityEmcal[patchType]->Fill(fPatchGrid[patchType].At(i),fActiveAreaCP[patchType].At(i));
	if(fActiveAreaCP[patchType].At(i)==0) Printf("WARNING activity: %d while En=%f active from MP: %d",fActiveAreaCP[patchType].At(i),fPatchGrid[patchType].At(i),fActiveAreaMPP[patchType].At(i));
      }
      if(type==1 && areaType==0) {
	fh1RhoDcal[patchType]->Fill(arr[c-1]);
	fPatchEnVsActivityDcal[patchType]->Fill(fPatchGrid[patchType].At(i),fActiveAreaCP[patchType].At(i));
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
    if(id>-1)
      fCellGrid.AddAt(fCellGrid.At(id)+track->Pt(),id);
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
    for(Int_t i = 0; i<(GetNCellsCol(type)-1); i+=2) {
      for(Int_t j = 0; j<(GetNCellsRow(type)-1); j+=2) {
	//loop over cells in mini patch
	for(Int_t k = 0; k<2; k++) {
	  for(Int_t l = 0; l<2; l++) {
	    Int_t id = GetGridID(i+k,j+l,type);
	    fMiniPatchGrid.AddAt(fMiniPatchGrid.At(nm)+fCellGrid.At(id),nm);
	    if(fCellGrid.At(id)>0.)
	      fActiveAreaMP.AddAt(fActiveAreaMP.At(nm)+1,nm);
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
  Int_t nm = (Int_t)(dim/2.);
  if(level==1 && fL1Slide)  nm = GetSlidingStepSizeMiniPatches(dim);
  //loop over edges of mini patches
  Int_t np = 0; //patch number
  for(Int_t type = 0; type<2; type++) {
    for(Int_t i = 0; i<=(GetNColMiniPatches(type)-nm); i+=nm) {
      for(Int_t j = 0; j<=(GetNRowMiniPatches(type)-nm); j+=nm) {
	//loop over mini patches in patch
	for(Int_t k = 0; k<nm; k++) {
	  for(Int_t l = 0; l<nm; l++) {
	    Int_t row = i+k;
	    Int_t col = j+l;
	    Int_t id = GetMiniPatchID(row,col,type);
	    fPatchGrid[pt].AddAt(fPatchGrid[pt].At(np)+fMiniPatchGrid.At(id),np);
	    if(fMiniPatchGrid.At(id)>0.) {
	      fActiveAreaMPP[pt].AddAt(fActiveAreaMPP[pt].At(np)+1,np);
	      fActiveAreaCP[pt].AddAt(fActiveAreaCP[pt].At(np)+fActiveAreaMP.At(id),np);
	    }
	  }
	}
	Double_t eta = 0.; Double_t phi = 0.;
	GetEtaPhiFromTriggerPatchID(np,dim,level,eta,phi);
	//Sanity check
	if(fActiveAreaMPP[pt].At(np)>1 && fActiveAreaCP[pt].At(np)==0) {
	  Printf("Activity MP: %d Cells: %d",fActiveAreaMPP[pt].At(np),fActiveAreaCP[pt].At(np));
	  Printf("pt: %d np: %d en: %f",pt,np,fPatchGrid[pt].At(np));
	}
	np++;
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetMiniPatchID(const Int_t row, const Int_t col, const Int_t type) const {
  Int_t id = row*GetNColMiniPatches(type) + col;
  //if in DCal move ID to after all EMCal IDs
  if(type==1) {
    Int_t ne = GetNColMiniPatches(0) * GetNRowMiniPatches(0);
    id+=ne;
  }
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
  //if in DCal move ID to after all EMCal IDs
  if(type==1) id+=fNCellsEMCal;
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
void AliEmcalPicoTrackInGridMaker::GetEtaPhiFromGridID(const Int_t id, Double_t &eta, Double_t &phi) const {
  //returns eta phi of cell at lower right edge (lowest eta, lowest phi)
  Int_t type = 0;
  Int_t useID = id;
  if(id>=fNCellsEMCal) {
    type = 1;
    useID = id-fNCellsEMCal;
  }

  Int_t row = TMath::FloorNint(useID/GetNCellsRow(type));
  Int_t col = useID - row * GetNCellsRow(type);
  eta = fEtaMin[type] + row*fCellSize;
  phi = fPhiMin[type] + col*fCellSize;

  AliDebug(2,Form("id: %d type: %d useID: %d row: %d col: %d eta: %f phi: %f",id,type,useID,row,col,eta,phi));
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::GetEtaPhiFromMiniPatchID(const Int_t id, Double_t &eta, Double_t &phi) const {
  //returns eta phi of mini patch at lower right edge (lowest eta, lowest phi)
  Int_t nMiniPatchesE = (Int_t)(0.25*fNCellsEMCal);
  Int_t type = 0;
  Int_t useID = id;

  if(id>=nMiniPatchesE) {
    type = 1;
    useID = id-nMiniPatchesE;
  }

  Int_t row = TMath::FloorNint(useID/GetNColMiniPatches(type));
  Int_t col = useID - row * GetNColMiniPatches(type);
  eta = fEtaMin[type] + row*2.*fCellSize;
  phi = fPhiMin[type] + col*2.*fCellSize;
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::GetEtaPhiFromTriggerPatchID(const Int_t id, const Int_t dim, const Int_t level, Double_t &eta, Double_t &phi) const {
  //returns eta phi of mini patch at lower right edge (lowest eta, lowest phi)
  Int_t nPatchesE = GetNTriggerPatches(0,dim,level);
  Int_t type = 0;
  Int_t useID = id;
  if(id>=nPatchesE) {
    type = 1;
    useID = id-nPatchesE;
  }
  Int_t row = TMath::FloorNint(useID/GetNColTriggerPatches(type,dim));
  Int_t col = useID - row * GetNColTriggerPatches(type,dim);
  eta = fEtaMin[type] + row*dim*fCellSize;
  phi = fPhiMin[type] + col*dim*fCellSize;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNColTriggerPatches(const Int_t type, const Int_t dim) const {
  //returns number of trigger patch columns
  Double_t stepmp = 0.5*(Double_t)GetSlidingStepSizeCells(dim);
  Int_t nmp = GetNColMiniPatches(type);
  Int_t ntc = TMath::FloorNint(nmp/stepmp);
  return ntc;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNRowTriggerPatches(const Int_t type, const Int_t dim) const {
  //returns number of trigger patch rows
  Double_t stepmp = 0.5*(Double_t)GetSlidingStepSizeCells(dim);
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
  
  fCellGrid.Set(fNCells);
  fCellGrid.Reset(0);
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::InitMiniPatches() {
  //initialize mini patch array
  if(fCellGrid.GetSize()<0) return kFALSE;
  Double_t conv = 0.25; //dimension of mini patch is 2x2 cells
  Int_t nMiniPatchesE = (Int_t)(conv*fNCellsEMCal);
  Int_t nMiniPatchesD = (Int_t)(conv*fNCellsDCal);

  //total number of mini patches
  Int_t nMiniPatches = nMiniPatchesE + nMiniPatchesD;

  fMiniPatchGrid.Set(nMiniPatches);
  fMiniPatchGrid.Reset(0.);

  fActiveAreaMP.Set(nMiniPatches);
  fActiveAreaMP.Reset(0);
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
  
  if(fCellGrid.GetSize()<0) return kFALSE;

  Int_t nPatchesE = GetNTriggerPatches(0,dim,level);
  Int_t nPatchesD = GetNTriggerPatches(1,dim,level);

  //total number of mini patches
  Int_t nPatches = nPatchesE + nPatchesD;

  Int_t type = GetPatchType(dim,level);
  if(type<0 || type>4) return kFALSE;
  fNPatchesEMCal[type] = nPatchesE;
  AliDebug(2,Form("Create trigger patch of type %d with dim %d and with %d patches EMCAL: %d DCAL: %d",type,dim,nPatches,nPatchesE,nPatchesD));
  fPatchGrid[type].Set(nPatches);
  fPatchGrid[type].Reset(0);
  fActiveAreaMPP[type].Set(nPatches);
  fActiveAreaMPP[type].Reset(0);
  fActiveAreaCP[type].Set(nPatches);
  fActiveAreaCP[type].Reset(0);

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
Double_t AliEmcalPicoTrackInGridMaker::GetPatchAreaActive(const Int_t id, const Int_t ipatch, const Int_t type) const {
  //type = 0 : active area from mini patches
  //type = 1 : active area from cells

  Int_t active = 0;
  if(type==0) active = fActiveAreaMPP[ipatch].At(id);
  else if(type==1) active = fActiveAreaCP[ipatch].At(id);
  else return -1;

  Double_t fac = 1.;
  if(type==0) fac = 4.;
  Double_t area = active*fac*fCellSize*fCellSize;
   return area;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetSlidingStepSizeCells(const Int_t dim) const {
  //get step size for mock-up L1 trigger
  if(!fL1Slide) return dim;

  if(dim==4) return 2;
  else if(dim==8) return 4;
  else if(dim==16) return 4;
  else if(dim==32) return 8;
  else return -1;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetSlidingStepSizeMiniPatches(const Int_t dim) const {
  //returns step size in mini patches
  Double_t step = (Double_t)GetSlidingStepSizeCells(dim);
  if(step<0) return step;
  Int_t stepMiniPatch = (Int_t)(step/2.);
  return stepMiniPatch;
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
