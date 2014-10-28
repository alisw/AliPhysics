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
#include "AliEmcalPicoTrackInGridMaker.h"

ClassImp(AliEmcalPicoTrackInGridMaker)

//________________________________________________________________________
AliEmcalPicoTrackInGridMaker::AliEmcalPicoTrackInGridMaker() : 
AliAnalysisTaskSE("AliEmcalPicoTrackInGridMaker"),
  fTracksOutName("PicoTracksInGrid"),
  fTracksInName("tracks"),
  fTracksIn(0),
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
  fMultVsRho(0),
  fHistList(0)
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
  }

  for(Int_t i = 0; i<3; i++) {
    fh2MedianTypeEmcal[i] = 0;
    fh2MedianTypeDcal[i] = 0;
    fpMedianTypeEmcal[i] = 0;
    fpMedianTypeDcal[i] = 0;
  }

}

//________________________________________________________________________
AliEmcalPicoTrackInGridMaker::AliEmcalPicoTrackInGridMaker(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksOutName("PicoTracksInGrid"),
  fTracksInName("tracks"),
  fTracksIn(0),
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
  fMultVsRho(0),
  fHistList(0)
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
  }

  for(Int_t i = 0; i<3; i++) {
    fh2MedianTypeEmcal[i] = 0;
    fh2MedianTypeDcal[i] = 0;
    fpMedianTypeEmcal[i] = 0;
    fpMedianTypeDcal[i] = 0;
  }

  // Output slot #1 write into a TList
  DefineOutput(1, TList::Class());
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

  OpenFile(1);
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);

  Int_t nBinsMed = 500;
  Double_t minMed = 0.;
  Double_t maxMed = 500.;

  // Int_t nBinsMult = 150;
  // Double_t minMult = 0.;
  // Double_t maxMult = 4000.;

  for(Int_t i = 0; i<3; i++) {
    fh2MedianTypeEmcal[i] = new TH2F(Form("fh2MedianTypeEmcalAreaType%d",i),Form("fh2MedianTypeEmcalAreaType%d",i),5,0.5,5.5,nBinsMed,minMed,maxMed);
    fHistList->Add(fh2MedianTypeEmcal[i]);

    fh2MedianTypeDcal[i] = new TH2F(Form("fh2MedianTypeDcalAreaType%d",i),Form("fh2MedianTypeDcalAreaType%d",i),5,0.5,5.5,nBinsMed,minMed,maxMed);
    fHistList->Add(fh2MedianTypeDcal[i]);

    fpMedianTypeEmcal[i] = new TProfile(Form("fpMedianTypeEmcalAreaType%d",i),Form("fpMedianTypeEmcalAreaType%d",i),5,0.5,5.5,"s");
    fHistList->Add(fpMedianTypeEmcal[i]);

    fpMedianTypeDcal[i] = new TProfile(Form("fpMedianTypeDcalAreaType%d",i),Form("fpMedianTypeDcalAreaType%d",i),5,0.5,5.5,"s");
    fHistList->Add(fpMedianTypeDcal[i]);
  }

  TString det[2] = {"Emcal","Dcal"};
  for(Int_t i = 0; i<5; i++) { //loop over patch types
    fh1RhoEmcal[i] = new TH1F(Form("fh1RhoEmcal_%d",i),Form("fh1RhoEmcal_%d",i),500,0.,1000.);
    fHistList->Add(fh1RhoEmcal[i]);
    fh1RhoDcal[i] = new TH1F(Form("fh1RhoDcal_%d",i),Form("fh1RhoDcal_%d",i),500,0.,1000.);
    fHistList->Add(fh1RhoDcal[i]);

    fPatchEnVsActivityEmcal[i] = new TH2F(Form("fh2PatchEnVsActivityEmcal_%d",i),Form("fh2PatchEnVsActivityEmcal_%d",i),300,0.,300.,150,-0.5,149.5);
    fHistList->Add(fPatchEnVsActivityEmcal[i]);

    fPatchEnVsActivityDcal[i] = new TH2F(Form("fh2PatchEnVsActivityDcal_%d",i),Form("fh2PatchEnVsActivityDcal_%d",i),300,0.,300.,150,-0.5,149.5);
    fHistList->Add(fPatchEnVsActivityDcal[i]);

    for(Int_t j = 0; j<2; j++) {
      fPatchECorr[j][i] = new TH1F(Form("fPatchECorr%s_%d",det[j].Data(),i),Form("fPatchECorr%s_%d;#it{E}_{patch}^{corr}",det[j].Data(),i),250,-50.,200.);
      fHistList->Add(fPatchECorr[j][i]);

      fPatchECorrPar[j][i] = new TH1F(Form("fPatchECorrPar%s_%d",det[j].Data(),i),Form("fPatchECorrPar%s_%d;#it{E}_{patch}^{corr}",det[j].Data(),i),250,-50.,200.);
      fHistList->Add(fPatchECorrPar[j][i]);  

      fPatchECorrRho[j][i] = new TH2F(Form("fPatchECorrRho%s_%d",det[j].Data(),i),Form("fPatchECorrRho%s_%d;#it{E}_{patch}^{corr};#rho",det[j].Data(),i),250,-50.,200.,500,0.,500.);
      fHistList->Add(fPatchECorrRho[j][i]);  

      fPatchECorrRhoDijet[j][i] = new TH2F(Form("fPatchECorrRhoDijet%s_%d",det[j].Data(),i),Form("fPatchECorrRhoDijet%s_%d;#it{E}_{patch}^{corr};#rho",det[j].Data(),i),250,-50.,200.,500,0.,500.);
      fHistList->Add(fPatchECorrRhoDijet[j][i]);

      fPatchECorrECorrRho[j][i] = new TH3F(Form("fPatchECorrECorrRho%s_%d",det[j].Data(),i),Form("fPatchECorrECorrRho%s_%d;#it{E}_{patch,det1}^{corr};#it{E}_{patch,det2}^{corr};#rho",det[j].Data(),i),210,-30.,180.,210,-30.,180.,250,0.,250.);
      fHistList->Add(fPatchECorrECorrRho[j][i]);
    }
  }

  fMultVsRho = new TH2F("fMultVsRho","fMultVsRho",3000,0,3000,400,0,400);
  fHistList->Add(fMultVsRho);

  PostData(1, fHistList);
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // retrieve tracks from input.
  if (!fTracksIn) { 
    fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));
    if (!fTracksIn) {
      AliError(Form("Could not retrieve tracks %s!", fTracksInName.Data())); 
      return;
    }
    if (!fTracksIn->GetClass()->GetBaseClass("AliVParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTracksInName.Data())); 
      return;
    }
  }

  Bool_t b = CreateGridCells();
  if(!b) return;
  b = CreateGridMiniPatches();
  if(!b) return;

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
    //    Int_t EleadID[2] = {-1,-1};
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
	  //	  EleadID[type] = j;
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
  }

  fMultVsRho->Fill(fTracksIn->GetEntriesFast(),medL1[3][0]);
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::CalculateSum(Int_t patchType) {

  Int_t n = fPatchGrid[patchType].GetSize();
  if(n<1) return -1.;

  Double_t sum = 0.;
  Int_t count = 0;
  for(Int_t i = 0; i<fPatchGrid[patchType].GetSize(); i++) {
    if(fPatchGrid[patchType].At(i)>0.) count++;
    sum+=fPatchGrid[patchType].At(i);
  }
  //  Double_t avg = -1.;
  //  if(count>0) avg = sum/((Double_t)count);
  return sum;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::CalculateMedian(Int_t patchType, Int_t type, Int_t areaType) {
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

  if(!InitCells()) return kFALSE;

  // loop over tracks
  const Int_t ntr = fTracksIn->GetEntriesFast();
  for (Int_t itr = 0; itr < ntr; ++itr) {

    AliVParticle *track = static_cast<AliVParticle*>(fTracksIn->At(itr));
    if (!track)
      continue;

    if(track->Pt()<fMinCellE) continue;

    Int_t id = GetGridID(track);
    if(id>-1) {
      //      Double_t old = fCellGrid.At(id);
      fCellGrid.AddAt(fCellGrid.At(id)+track->Pt(),id);
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::CreateGridMiniPatches() {

  if(!InitMiniPatches()) return kFALSE;

  //Int_t nMiniPatches = fMiniPatchGrid.GetSize();
  //Int_t nMiniPatchesE = (Int_t)(0.25*fNCellsEMCal);
  
  Int_t type = 0;
  // Int_t nMiniPatchesRow = TMath::FloorNint(0.5*GetNCellsRow(type));
  //loop over edges of mini patches
  Int_t nm = 0; //mini patch number
  for(Int_t i = 0; i<(GetNCellsCol(type)-1); i+=2) {
    for(Int_t j = 0; j<(GetNCellsRow(type)-1); j+=2) {
      //loop over cells in mini patch
      for(Int_t k = 0; k<2; k++) {
	for(Int_t l = 0; l<2; l++) {
	  Int_t row = i+k;
	  Int_t col = j+l;
	  Int_t id = GetGridID(row,col,type);
	  fMiniPatchGrid.AddAt(fMiniPatchGrid.At(nm)+fCellGrid.At(id),nm);
	  if(fCellGrid.At(id)>0.)
	    fActiveAreaMP.AddAt(fActiveAreaMP.At(nm)+1,nm);
	}
      }
      nm++;
    }
  }

  //DCal mini patches
  type = 1;
  for(Int_t i = 0; i<(GetNCellsCol(type)-1); i+=2) {
    for(Int_t j = 0; j<(GetNCellsRow(type)-1); j+=2) {
      //loop over cells in mini patch
      for(Int_t k = 0; k<2; k++) {
	for(Int_t l = 0; l<2; l++) {
	  Int_t row = i+k;
	  Int_t col = j+l;
	  Int_t id = GetGridID(row,col,type);
	  fMiniPatchGrid.AddAt(fMiniPatchGrid.At(nm)+fCellGrid.At(id),nm);
	  if(fCellGrid.At(id)>0.) 
	    fActiveAreaMP.AddAt(fActiveAreaMP.At(nm)+1,nm);
	}
      }
      nm++;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNRowMiniPatches(Int_t type) {
  //returns number of rows of mini patches in detector of type X (0: EMCal 1: DCal)
  Int_t nRows = TMath::FloorNint(0.5*GetNCellsCol(type));
  return nRows;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNColMiniPatches(Int_t type) {
  //returns number of rows of mini patches in detector of type X (0: EMCal 1: DCal)
  Int_t nCols = TMath::FloorNint(0.5*GetNCellsRow(type));
  return nCols;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::CreateGridPatches(Int_t dim, Int_t level) {

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
Int_t AliEmcalPicoTrackInGridMaker::GetMiniPatchID(Int_t row, Int_t col, Int_t type) {
  Int_t id = row*GetNColMiniPatches(type) + col;
  //if in DCal move ID to after all EMCal IDs
  if(type==1) {
    Int_t ne = GetNColMiniPatches(0) * GetNRowMiniPatches(0);
    id+=ne;
  }
  return id;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetCellType(Double_t eta, Double_t phi) {

  for(Int_t i = 0; i<2; i++) {
    if(eta>fEtaMin[i] && eta<fEtaMax[i] && phi>fPhiMin[i] && phi<fPhiMax[i]) return i;
  }
  return -1;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetGridID(Int_t row, Int_t col, Int_t type) {
  Int_t id = row*GetNCellsRow(type) + col;
  //if in DCal move ID to after all EMCal IDs
  if(type==1) id+=fNCellsEMCal;

  return id;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetGridID(Double_t eta, Double_t phi) {
  
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
Int_t AliEmcalPicoTrackInGridMaker::GetNCellsCol(Int_t type) {

  Double_t deta = fEtaMax[type] - fEtaMin[type];
  Int_t nCellsCol = TMath::FloorNint(deta/fCellSize);
  return nCellsCol;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetNCellsRow(Int_t type) {

  Double_t dPhi = fPhiMax[type] - fPhiMin[type];
  Int_t nCellsRow = TMath::FloorNint(dPhi/fCellSize);
  return nCellsRow;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::InitCells() {

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

  /*
    Printf("EMCal: %d x %d",nCellsEtaE,nCellsPhiE);
    Printf("DCal: row: %d x col: %d",nCellsEtaD,nCellsPhiD);
    Printf("fNCells: %d fNCellsE: %d fnCellsD: %d",fNCells,fNCellsEMCal,fNCellsDCal);
  */
  fCellGrid.Set(fNCells);
  fCellGrid.Reset(0);
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalPicoTrackInGridMaker::InitMiniPatches() {
  
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
Bool_t AliEmcalPicoTrackInGridMaker::InitPatches(Int_t dim, Int_t level) {
  // dimensions in cell units
  // if level==1: sliding window will be applied
  // L1 4x4: slide by 2 cells: 1 mini patch
  // L1 8x8: slide by 4 cells: 2 mini patches
  // L1 16x16: slide by 4 cells: 2 mini patches
  // L1 32x32: slide by 8 cells: 4 mini patches
  
  if(fCellGrid.GetSize()<0) return kFALSE;

  Double_t dimd = (Double_t)dim;
  Double_t conv = 1./(dimd*dimd);
  if(level==1) {
    Double_t step = (Double_t)GetSlidingStepSizeCells(dim);
    conv = 1./(step*step);
  }
  Int_t nPatchesE = (Int_t)(conv*fNCellsEMCal);
  Int_t nPatchesD = (Int_t)(conv*fNCellsDCal);

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
Double_t AliEmcalPicoTrackInGridMaker::GetPatchArea(Int_t ipatch) {

  Double_t ncell = 4.;
  if(ipatch==0) ncell = 4.;
  if(ipatch==1) ncell = 4.;
  if(ipatch==2) ncell = 8.;
  if(ipatch==3) ncell = 16.;
  if(ipatch==4) ncell = 32.;
  Double_t area = ncell*ncell*fCellSize*fCellSize;
  return area;
}

//________________________________________________________________________
Double_t AliEmcalPicoTrackInGridMaker::GetPatchAreaActive(Int_t id, Int_t ipatch, Int_t type) {
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
Int_t AliEmcalPicoTrackInGridMaker::GetSlidingStepSizeCells(Int_t dim) {

  if(!fL1Slide) return dim;

  if(dim==4) return 2;
  else if(dim==8) return 4;
  else if(dim==16) return 4;
  else if(dim==32) return 8;
  else return -1;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetSlidingStepSizeMiniPatches(Int_t dim) {
  //returns step size in mini patches
  Double_t step = (Double_t)GetSlidingStepSizeCells(dim);
  if(step<0) return step;
  Int_t stepMiniPatch = (Int_t)(step/2.);
  return stepMiniPatch;
}

//________________________________________________________________________
Int_t AliEmcalPicoTrackInGridMaker::GetPatchType(Int_t dim, Int_t level) {
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
  //
  //Check if defined edges of EMCal and DCal make sense
  //
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
    Double_t nCellsColExact = dphi/fCellSize;
    
    Double_t deta = fEtaMax[type] - fEtaMin[type];
    Double_t nCellsRowExact = deta/fCellSize;

    Double_t col_extra = nCellsColExact - TMath::Floor(nCellsColExact);
    Double_t phi_extra = col_extra*fCellSize;
    fPhiMin[type] += phi_extra*0.5;
    fPhiMax[type] -= phi_extra*0.5;

    Double_t row_extra = nCellsRowExact - TMath::Floor(nCellsRowExact);
    Double_t eta_extra = row_extra*fCellSize;
    fEtaMin[type] += eta_extra*0.5;
    fEtaMax[type] -= eta_extra*0.5;
    AliDebug(2,Form("type: %d  exact: col: %f row: %f",type,nCellsColExact,nCellsRowExact));
    //PrintAcceptance();
  }
  return kTRUE;
}

//________________________________________________________________________
void AliEmcalPicoTrackInGridMaker::PrintAcceptance() {

  Printf("EMCal");
  Printf("phi: %f-%f",fPhiMin[0],fPhiMax[0]);
  Printf("eta: %f-%f",fEtaMin[0],fEtaMax[0]);

  Printf("DCal");
  Printf("phi: %f-%f",fPhiMin[1],fPhiMax[1]);
  Printf("eta: %f-%f",fEtaMin[1],fEtaMax[1]);

}
