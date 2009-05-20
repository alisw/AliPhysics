/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Boris Polishchuk                                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//---------------------------------------------------------------------------// 
// Fill histograms with two-cluster invariant mass                           //
// using calibration coefficients of the previous iteration.                 //
//---------------------------------------------------------------------------//

#include "AliAnalysisTaskPi0CalibSelection.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSPIDv1.h"
#include "TRefArray.h"
#include "TList.h"
#include "TH1.h"

ClassImp(AliAnalysisTaskPi0CalibSelection)

AliAnalysisTaskPi0CalibSelection::AliAnalysisTaskPi0CalibSelection() :
AliAnalysisTaskSE(),fOutputContainer(0x0),fRecoParam(0x0),fPhosGeo(0x0),fHmgg(0x0),
  fEmin(0.)
{
  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }
  
}

AliAnalysisTaskPi0CalibSelection::AliAnalysisTaskPi0CalibSelection(const char* name) :
  AliAnalysisTaskSE(name),fOutputContainer(0x0),fRecoParam(0x0),fPhosGeo(0x0),fHmgg(0x0),
  fEmin(0.)
{
  DefineOutput(1,TList::Class());

  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }

}

AliAnalysisTaskPi0CalibSelection::~AliAnalysisTaskPi0CalibSelection()
{
  
  if(fOutputContainer){
    fOutputContainer->Clear() ; 
    delete fOutputContainer ;
  }
}

void AliAnalysisTaskPi0CalibSelection::UserCreateOutputObjects()
{
  fOutputContainer = new TList();
  
  char hname[128], htitl[128];
  
  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	sprintf(hname,"%d_%d_%d",iMod,iX,iZ);
	sprintf(htitl,"Two-gamma inv. mass for mod %d, cell (%d,%d)",iMod,iX,iZ);
	fHmpi0[iMod][iX][iZ] = new TH1F(hname,htitl,100,0.,300.);
	fOutputContainer->Add(fHmpi0[iMod][iX][iZ]);
      }
    }
  }

  fHmgg = new TH1F("hmgg","2-cluster invariant mass",100,0.,300.);
  fOutputContainer->Add(fHmgg);
}

void AliAnalysisTaskPi0CalibSelection::UserExec(Option_t* /* option */)
{
  //Analysis per event.
  
  AliVEvent* event = InputEvent();
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event) ;

  Double_t v[3] ; //vertex ;
  esd->GetVertex()->GetXYZ(v) ;
  TVector3 vtx(v);
  printf("Vertex: (%.3f,%.3f,%.3f)\n",vtx.X(),vtx.Y(),vtx.Z());
 
  Int_t runNum = esd->GetRunNumber();
  printf("Run number: %d\n",runNum);
  
  if(!fPhosGeo) {
    
    AliCDBEntry* entryGeo = AliCDBManager::Instance()->Get("GRP/Geometry/Data",runNum);
    if(!entryGeo) AliFatal("No Geometry entry found in OCDB!!");
    
    TGeoManager* geoManager = dynamic_cast<TGeoManager*>(entryGeo->GetObject());
    if(!geoManager) AliFatal("No valid TGeoManager object found in the OCDB.");
    
    gGeoManager = geoManager;
    fPhosGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
  }
  
  //Calibrations from previous iteration
  AliPHOSCalibData calibData(runNum);
  
  //Get RecoParameters. See AliQADataMakerRec::InitRecoParams().
  if(!fRecoParam) {
    
    AliCDBEntry* entryRecoPar = AliCDBManager::Instance()->Get("PHOS/Calib/RecoParam",runNum);
    if(!entryRecoPar) AliFatal("No RecoParam entry in OCDB!!");
    
    TObject * recoParamObj = entryRecoPar->GetObject() ;
    AliDetectorRecoParam* rprm = 0;
    
    if (dynamic_cast<TObjArray*>(recoParamObj)) {
      TObjArray *recoParamArray = dynamic_cast<TObjArray*>(recoParamObj) ;
      for (Int_t iRP=0; iRP<recoParamArray->GetEntriesFast(); iRP++) {
	rprm = dynamic_cast<AliDetectorRecoParam*>(recoParamArray->At(iRP)) ;
	if (rprm->IsDefault()) break;
      }
    } 
    
    if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
      dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
      rprm = dynamic_cast<AliDetectorRecoParam*>(recoParamObj) ;
    }
    
    if(!rprm) AliFatal("No valid RecoParam object found in the OCDB.");
    fRecoParam = dynamic_cast<AliPHOSRecoParam*>(rprm);
    if(!fRecoParam) AliFatal("recoparams are _NOT_ of type AliPHOSRecoParam!!");
  }
  
  Float_t logWeight = fRecoParam->GetEMCLogWeight();
  printf("Will use logWeight %.3f .\n",logWeight);

  AliPHOSPIDv1 pid;

  Int_t mod1 = -1;
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p12;
  Int_t relid[4] ;
  Int_t maxId;

  TRefArray * caloClustersArr  = new TRefArray();
  esd->GetPHOSClusters(caloClustersArr);
  
  const Int_t kNumberOfPhosClusters   = caloClustersArr->GetEntries() ;
  printf("CaloClusters: %d\n", kNumberOfPhosClusters);

  // PHOS cells
  AliESDCaloCells *phsCells = esd->GetPHOSCells();
  
  // loop over PHOS clusters
  for(Int_t iClu=0; iClu<kNumberOfPhosClusters; iClu++) {
    
    AliESDCaloCluster *c1 = (AliESDCaloCluster *) caloClustersArr->At(iClu);
    if(!c1->IsPHOS()) continue; // EMCAL cluster!

    Float_t E1_i = c1->E();   // cluster energy before correction
    if(E1_i<fEmin) continue;

    AliPHOSEsdCluster clu1(*c1);
    clu1.Recalibrate(&calibData, phsCells);
    clu1.EvalAll(logWeight,vtx);
    clu1.EnergyCorrection(&pid) ;

    clu1.GetMomentum(p1,v);

    MaxEnergyCellPos(phsCells,&clu1,maxId);
    fPhosGeo->AbsToRelNumbering(maxId, relid);

    mod1 = relid[0]-1;        // module
    Int_t iX = relid[2]-1;    // cluster X-coord
    Int_t iZ = relid[3]-1;    // cluster Z-coord
    Float_t E1_ii = clu1.E(); // cluster energy after correction
    
    for (Int_t jClu=iClu; jClu<kNumberOfPhosClusters; jClu++) {
      AliESDCaloCluster *c2 = (AliESDCaloCluster *) caloClustersArr->At(jClu);
      if(!c2->IsPHOS())   continue; // EMCAL cluster!
      if(c2->IsEqual(c1)) continue;

      Float_t E2_i = c2->E();
      if(E2_i<fEmin) continue;
      
      AliPHOSEsdCluster clu2(*c2);
      clu2.Recalibrate(&calibData, phsCells);
      clu2.EvalAll(logWeight,vtx);
      clu2.EnergyCorrection(&pid) ;
        
      clu2.GetMomentum(p2,v);
      Float_t E2_ii = clu2.E();

      p12 = p1+p2;

      fHmgg->Fill(p12.M()*1000); 
      fHmpi0[mod1][iX][iZ]->Fill(p12.M()*1000);
      
      printf("Mass in (mod%d,%d,%d): %.3f GeV  E1_i=%f E1_ii=%f  E2_i=%f E2_ii=%f\n",
	     mod1,iX,iZ,p12.M(),E1_i,E1_ii,E2_i,E2_ii);
    }
    
  } // end of loop over PHOS clusters
  
  delete caloClustersArr;
  PostData(1,fOutputContainer);
}

void AliAnalysisTaskPi0CalibSelection::MaxEnergyCellPos(AliESDCaloCells *cells, AliESDCaloCluster* clu, Int_t& maxId)
{
  
  Double_t eMax = -111;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    Int_t cellAbsId = clu->GetCellAbsId(iDig);
    Double_t eCell = cells->GetCellAmplitude(cellAbsId)*clu->GetCellAmplitudeFraction(iDig);
    if(eCell>eMax)  { 
      eMax = eCell; 
      maxId = cellAbsId;
    }
  }

}
