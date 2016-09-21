// $Id$
//
// Class to put cells into trigger patches
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TH3F.h>

#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALTriggerConstants.h"

#include "AliEmcalPatchFromCellMaker.h"

ClassImp(AliEmcalPatchFromCellMaker)

//________________________________________________________________________
AliEmcalPatchFromCellMaker::AliEmcalPatchFromCellMaker() : 
  AliAnalysisTaskEmcal("AliEmcalPatchFromCellMaker",kTRUE),
  fCaloTriggersOutName("EmcalPatches32x32"),
  fCaloTriggersOut(0),
  fPatchDim(32),
  fMinCellE(0.05),
  fCellTimeMin(485e-9),
  fCellTimeMax(685e-9),
  fL1Slide(0),
  fTriggerBitConfig(0x0),
  fh3EEtaPhiCell(0),
  fh2CellEnergyVsTime(0),
  fh1CellEnergySum(0)
{
  // Constructor.
 for (Int_t i = 0; i < kPatchCols; i++) {
    for (Int_t j = 0; j < kPatchRows; j++) {
      fPatchADCSimple[i][j] = 0.;
      fPatchESimple[i][j] = 0.;
    }
 }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliEmcalPatchFromCellMaker::AliEmcalPatchFromCellMaker(const char *name) : 
  AliAnalysisTaskEmcal(name,kTRUE),
  fCaloTriggersOutName("EmcalPatches32x32"),
  fCaloTriggersOut(0),
  fPatchDim(32),
  fMinCellE(0.05),
  fCellTimeMin(485e-9),
  fCellTimeMax(685e-9),
  fL1Slide(0),
  fTriggerBitConfig(0x0),
  fh3EEtaPhiCell(0),
  fh2CellEnergyVsTime(0),
  fh1CellEnergySum(0)
{
  // Constructor.
 for (Int_t i = 0; i < kPatchCols; i++) {
    for (Int_t j = 0; j < kPatchRows; j++) {
      fPatchADCSimple[i][j] = 0.;
      fPatchESimple[i][j] = 0.;
    }
 }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliEmcalPatchFromCellMaker::~AliEmcalPatchFromCellMaker()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalPatchFromCellMaker::ExecOnce() 
{
  // Init the analysis.

  AliAnalysisTaskEmcal::ExecOnce();

  if (!fLocalInitialized)
    return;

  if (!fCaloTriggersOutName.IsNull()) {
    fCaloTriggersOut = new TClonesArray("AliEMCALTriggerPatchInfo");
    fCaloTriggersOut->SetName(fCaloTriggersOutName);

    if (!(InputEvent()->FindListObject(fCaloTriggersOutName))) {
      InputEvent()->AddObject(fCaloTriggersOut);
    }
    else {
      fLocalInitialized = kFALSE;
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fCaloTriggersOutName.Data()));
      return;
    }
  }

  if(!fTriggerBitConfig)
    fTriggerBitConfig = new AliEMCALTriggerBitConfigNew();

}

//________________________________________________________________________
void AliEmcalPatchFromCellMaker::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  Int_t fgkNPhiBins = 18*8;
  Float_t kMinPhi   = 0.;
  Float_t kMaxPhi   = 2.*TMath::Pi();
  Double_t *binsPhi = new Double_t[fgkNPhiBins+1];
  for(Int_t i=0; i<=fgkNPhiBins; i++) binsPhi[i]=(Double_t)kMinPhi + (kMaxPhi-kMinPhi)/fgkNPhiBins*(Double_t)i ;

  Int_t fgkNEtaBins = 100;
  Float_t fgkEtaMin = -1.;
  Float_t fgkEtaMax =  1.;
  Double_t *binsEta=new Double_t[fgkNEtaBins+1];
  for(Int_t i=0; i<=fgkNEtaBins; i++) binsEta[i]=(Double_t)fgkEtaMin + (fgkEtaMax-fgkEtaMin)/fgkNEtaBins*(Double_t)i ;

  Int_t fgkNTimeBins = 600;
  Float_t kMinTime   = -200.;
  Float_t kMaxTime   = 1000;
  Double_t *binsTime = new Double_t[fgkNTimeBins+1];
  for(Int_t i=0; i<=fgkNTimeBins; i++) binsTime[i]=(Double_t)kMinTime + (kMaxTime-kMinTime)/fgkNTimeBins*(Double_t)i ;

  Double_t enBinEdges[3][2];
  enBinEdges[0][0] = 1.; //10 bins
  enBinEdges[0][1] = 0.1;
  enBinEdges[1][0] = 5.; //8 bins
  enBinEdges[1][1] = 0.5;
  enBinEdges[2][0] = 100.;//95 bins
  enBinEdges[2][1] = 1.;

  const Float_t enmin1 =  0;
  const Float_t enmax1 =  enBinEdges[0][0];
  const Float_t enmin2 =  enmax1 ;
  const Float_t enmax2 =  enBinEdges[1][0];
  const Float_t enmin3 =  enmax2 ;
  const Float_t enmax3 =  enBinEdges[2][0];//fgkEnMax;
  const Int_t nbin11 = (int)((enmax1-enmin1)/enBinEdges[0][1]);
  const Int_t nbin12 = (int)((enmax2-enmin2)/enBinEdges[1][1])+nbin11;
  const Int_t nbin13 = (int)((enmax3-enmin3)/enBinEdges[2][1])+nbin12;

  Int_t fgkNEnBins=nbin13;
  Double_t *binsEn=new Double_t[fgkNEnBins+1];
  for(Int_t i=0; i<=fgkNEnBins; i++) {
    if(i<=nbin11) binsEn[i]=(Double_t)enmin1 + (enmax1-enmin1)/nbin11*(Double_t)i ;
    if(i<=nbin12 && i>nbin11) binsEn[i]=(Double_t)enmin2 + (enmax2-enmin2)/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;
    if(i<=nbin13 && i>nbin12) binsEn[i]=(Double_t)enmin3 + (enmax3-enmin3)/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;
  }

  fh3EEtaPhiCell = new TH3F("fh3EEtaPhiCell","fh3EEtaPhiCell;E_{cell};#eta;#phi",fgkNEnBins,binsEn,fgkNEtaBins,binsEta,fgkNPhiBins,binsPhi);
  fOutput->Add(fh3EEtaPhiCell);

  fh2CellEnergyVsTime = new TH2F("fh2CellEnergyVsTime","fh2CellEnergyVsTime;E_{cell};time",fgkNEnBins,binsEn,fgkNTimeBins,binsTime);
  fOutput->Add(fh2CellEnergyVsTime);

  fh1CellEnergySum = new TH1F("fh1CellEnergySum","fh1CellEnergySum;E_{cell};time",fgkNEnBins,binsEn);
  fOutput->Add(fh1CellEnergySum);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.

  if(binsEn)                delete [] binsEn;
  if(binsPhi)               delete [] binsPhi;
  if(binsEta)               delete [] binsEta;
  if(binsTime)              delete [] binsTime;
}

//________________________________________________________________________
Bool_t AliEmcalPatchFromCellMaker::Run() 
{
  // Main loop, called for each event.
  
  fCaloTriggersOut->Delete();

  if (!fCaloCells) {
    AliError(Form("Calo cells container %s not available.", fCaloCellsName.Data()));
    return kFALSE;
  }

 for (Int_t i = 0; i < kPatchCols; i++) {
    for (Int_t j = 0; j < kPatchRows; j++) {
      fPatchADCSimple[i][j] = 0.;
      fPatchESimple[i][j] = 0.;
    }
 }

  if(!FillPatchADCSimple()) {
    AliError(Form("%s Could not create simple ADC patches",GetName()));
    return kFALSE;
  }

  RunSimpleOfflineTrigger();

  Double_t sum = 0.;
  for (Int_t i = 0; i < kPatchCols; i++) {
    for (Int_t j = 0; j < kPatchRows; j++) {
      sum+=fPatchESimple[i][j];
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalPatchFromCellMaker::FillPatchADCSimple()
{

  // fill the array for offline trigger processing

  // fill the patch ADCs from cells
  Double_t sum = 0.;
  Int_t nCell = fCaloCells->GetNumberOfCells();
  for(Int_t iCell = 0; iCell < nCell; ++iCell) {
    // get the cell info, based in index in array
    Short_t cellId = fCaloCells->GetCellNumber(iCell);

    Double_t cellT = fCaloCells->GetCellTime(cellId); 
    Double_t amp = fCaloCells->GetAmplitude(iCell);
    fh2CellEnergyVsTime->Fill(amp,cellT*1e9);

    //timing cuts
    if(cellT<fCellTimeMin || cellT>fCellTimeMax) continue;
    //energy cut
    if(amp<fMinCellE) continue;
    sum+=amp;

    // get position
    Int_t absId=-1;
    fGeom->GetFastORIndexFromCellIndex(cellId, absId);
    Int_t globCol=-1, globRow=-1;
    fGeom->GetPositionInEMCALFromAbsFastORIndex(absId, globCol, globRow);
    // add
    fPatchADCSimple[globCol][globRow] += amp/EMCALTrigger::kEMCL1ADCtoGeV;
    fPatchESimple[globCol][globRow] += amp;

    TVector3 pos;
    fGeom->GetGlobal(cellId, pos);
    TLorentzVector lv(pos,amp);
    Double_t cellEta = lv.Eta();
    Double_t cellPhi = lv.Phi();
    if(cellPhi<0.) cellPhi+=TMath::TwoPi();
    if(cellPhi>TMath::TwoPi()) cellPhi-=TMath::TwoPi();
    fh3EEtaPhiCell->Fill(amp,cellEta,cellPhi); 
  }
  fh1CellEnergySum->Fill(sum);

  return kTRUE;
}

//________________________________________________________________________
void AliEmcalPatchFromCellMaker::RunSimpleOfflineTrigger()
{
  // Runs a simple offline trigger algorithm.
  // It creates separate patches with dimension fPatchDim

  // run the trigger algo, stepping by stepsize (in trigger tower units)
  Int_t itrig = 0;
  Int_t patchSize = GetDimFastor();
  Int_t stepSize = GetSlidingStepSizeFastor();
  Int_t maxCol = kPatchCols - patchSize;
  //  Int_t maxRow = kPatchRows - patchSize;
  Int_t maxRow = fGeom->GetNTotalTRU()*2 - patchSize; //apparently this is the total TRU in phi ??
  //  Printf("fGeom->GetNTotalTRU(): %d %d = %d x %d ;  %d x %d",fGeom->GetNTRU(),fGeom->GetNTotalTRU(),fGeom->GetNTRUEta(),fGeom->GetNTRUPhi(),fGeom->GetNModulesInTRUEta(),fGeom->GetNModulesInTRUPhi());

  for (Int_t i = 0; i <= maxCol; i += stepSize) {
    for (Int_t j = 0; j <= maxRow; j += stepSize) {
      // get the trigger towers composing the patch
      Int_t   adcAmp = 0;
      Double_t enAmp = 0.;
      // window
      for (Int_t k = 0; k < patchSize; ++k) {
        for (Int_t l = 0; l < patchSize; ++l) {
	  // add amplitudes
          adcAmp += (ULong64_t)fPatchADCSimple[i+k][j+l];
	  enAmp += fPatchESimple[i+k][j+l];
	}
      }

      if (adcAmp == 0) {
	AliDebug(2,"EMCal trigger patch with 0 ADC counts.");
	continue;
      }

      Int_t absId=-1;
      Int_t cellAbsId[4]={-1,-1,-1,-1};

      // get low left edge (eta max, phi min)
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(i, j, absId);
      // convert to the 4 absId of the cells composing the trigger channel
      fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
      TVector3 edge1;
      fGeom->GetGlobal(cellAbsId[0], edge1);

      // get up right edge (eta min, phi max)
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(i+patchSize-1, j+patchSize-1, absId);
      fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
      TVector3 edge2;
      fGeom->GetGlobal(cellAbsId[3], edge2);

      // get the center of the patch
      Int_t offsetCenter = TMath::FloorNint(0.5*patchSize);
      fGeom->GetAbsFastORIndexFromPositionInEMCAL(i+offsetCenter-1, j+offsetCenter-1, absId);
      fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
      TVector3 center1;
      fGeom->GetGlobal(cellAbsId[3], center1);

      fGeom->GetAbsFastORIndexFromPositionInEMCAL(i+offsetCenter, j+offsetCenter, absId);
      fGeom->GetCellIndexFromFastORIndex(absId, cellAbsId);
      TVector3 center2;
      fGeom->GetGlobal(cellAbsId[0], center2);

      TVector3 centerGeo(center1);
      centerGeo += center2;
      centerGeo *= 0.5;

      // save the trigger object
      AliEMCALTriggerPatchInfo *trigger =
          new ((*fCaloTriggersOut)[itrig]) AliEMCALTriggerPatchInfo();
      itrig++;
      trigger->SetTriggerBitConfig(fTriggerBitConfig);
      trigger->SetCenterGeo(centerGeo, enAmp);
      trigger->SetEdge1(edge1, enAmp);
      trigger->SetEdge2(edge2, enAmp);
      trigger->SetADCAmp(adcAmp);
      trigger->SetEdgeCell(i*2, j*2); // from triggers to cells
    }
  } // trigger algo
  AliDebug(2,Form("Created %d trigger patches (%d) in this event",itrig,patchSize));

}

//________________________________________________________________________
Int_t AliEmcalPatchFromCellMaker::GetDimFastor() const {

  Int_t dim = TMath::FloorNint((Double_t)(fPatchDim/2.));
  return dim;
}

//________________________________________________________________________
Int_t AliEmcalPatchFromCellMaker::GetSlidingStepSizeFastor() const {

  Int_t dim = GetDimFastor();
  if(!fL1Slide) return dim;

  if(dim==2) return 2;
  else if(dim==4) return 4;
  else if(dim==8) return 4;
  else if(dim==16) return 8;
  else return -1;
}




