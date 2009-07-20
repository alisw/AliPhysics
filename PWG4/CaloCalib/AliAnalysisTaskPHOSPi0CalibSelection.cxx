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
//                                                                           //
// Fill histograms (one per cell) with two-cluster invariant mass            //
// using calibration coefficients of the previous iteration.                 //
// Histogram for a given cell is filled if the most energy of one cluster    //
// is deposited in this cell and the other cluster could be anywhere in PHOS.//
//                                                                           //
//---------------------------------------------------------------------------//

#include <cstdlib>

// Root 
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRefArray.h"
#include "TList.h"

// AliRoot
#include "AliAnalysisTaskPHOSPi0CalibSelection.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSPIDv1.h"

ClassImp(AliAnalysisTaskPHOSPi0CalibSelection)

AliAnalysisTaskPHOSPi0CalibSelection::AliAnalysisTaskPHOSPi0CalibSelection() :
AliAnalysisTaskSE(),fOutputContainer(0x0),fPhosGeo(0x0),fCalibData(0x0),fHmgg(0x0),
  fEmin(0.), fLogWeight(4.5)
{
  //Default constructor.

  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }
		
}

AliAnalysisTaskPHOSPi0CalibSelection::AliAnalysisTaskPHOSPi0CalibSelection(const char* name) :
  AliAnalysisTaskSE(name),fOutputContainer(0x0),fPhosGeo(0x0),fCalibData(0x0),fHmgg(0x0),
  fEmin(0.), fLogWeight(4.5)
{
  //Named constructor which should be used.
  
  DefineOutput(1,TList::Class());
  
  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }

}

AliAnalysisTaskPHOSPi0CalibSelection::~AliAnalysisTaskPHOSPi0CalibSelection()
{
  //Destructor.
  
  if(fOutputContainer){
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
	
  if(fCalibData) delete fCalibData;
  if(fPhosGeo)   delete fPhosGeo;
	
}

void AliAnalysisTaskPHOSPi0CalibSelection::UserCreateOutputObjects()
{
  //Create output container
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
	
  fCalibData = new AliPHOSCalibData();
  fPhosGeo =  AliPHOSGeometry::GetInstance("IHEP") ;	

}

void AliAnalysisTaskPHOSPi0CalibSelection::UserExec(Option_t* /* option */)
{
  //Analysis per event.
  if(DebugLevel() > 1) printf("AliAnalysisTaskPHOSPi0CalibSelection <<< Event %d >>>\n",(Int_t)Entry());
  AliAODEvent* aod = 0x0;
  if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) aod = dynamic_cast<AliAODEvent*>(InputEvent());
  else  if(!strcmp(InputEvent()->GetName(),"AliESDEvent")) aod = AODEvent();
  else {
	  printf("AliAnalysisTaskPHOSPi0CalibSelection: Unknown event type, STOP!\n");
	  abort();
  }	

  Double_t v[] = {aod->GetVertex(0)->GetX(),aod->GetVertex(0)->GetY(),aod->GetVertex(0)->GetZ()}; //to check!!
  //aod->GetVertex()->GetXYZ(v) ;
  TVector3 vtx(v); //Check
	
  if(DebugLevel() > 1) printf("AliAnalysisTaskPHOSPi0CalibSelection Vertex: (%.3f,%.3f,%.3f)\n",vtx.X(),vtx.Y(),vtx.Z());
 
  Int_t runNum = aod->GetRunNumber();
  if(DebugLevel() > 1) printf("Run number: %d\n",runNum);
	
  //Get the matrix with geometry information
  //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
  if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) {
	  if(DebugLevel() > 1) 
		  printf("AliAnalysisTaskPHOSPi0CalibSelection Use ideal geometry, values geometry matrix not kept in AODs.\n");
  }
  else{	
	  if(DebugLevel() > 1) printf("AliAnalysisTaskPHOSPi0CalibSelection Load Misaligned matrices. \n");
	  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
	  for(Int_t mod=0; mod<5; mod++){ 
		if(esd->GetPHOSMatrix(mod))
			fPhosGeo->SetMisalMatrix(esd->GetPHOSMatrix(mod),mod) ;
		}
  }
	
  if(DebugLevel() > 1) printf("AliAnalysisTaskPHOSPi0CalibSelection Will use fLogWeight %.3f .\n",fLogWeight);

  AliPHOSPIDv1 pid;

  Int_t mod1 = -1;
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p12;
  Int_t relid[4] ;
  Int_t maxId;

  TRefArray * caloClustersArr  = new TRefArray();
  aod->GetPHOSClusters(caloClustersArr);
  
  const Int_t kNumberOfPhosClusters   = caloClustersArr->GetEntries() ;
  if(DebugLevel() > 1) printf("AliAnalysisTaskPHOSPi0CalibSelection CaloClusters: %d\n", kNumberOfPhosClusters);

  // PHOS cells
  AliAODCaloCells *phsCells = aod->GetPHOSCells();
  
  // loop over PHOS clusters
  for(Int_t iClu=0; iClu<kNumberOfPhosClusters; iClu++) {
    
    AliAODCaloCluster *c1 = (AliAODCaloCluster *) caloClustersArr->At(iClu);
    if(!c1->IsPHOSCluster()) continue; // EMCAL cluster!

    Float_t e1i = c1->E();   // cluster energy before correction
    if(e1i<fEmin) continue;

    AliPHOSAodCluster clu1(*c1);
    clu1.Recalibrate(fCalibData, phsCells);
    clu1.EvalAll(fLogWeight,vtx);
    clu1.EnergyCorrection(&pid) ;

    clu1.GetMomentum(p1,v);

    MaxEnergyCellPos(phsCells,&clu1,maxId);
    fPhosGeo->AbsToRelNumbering(maxId, relid);

    mod1 = relid[0]-1;        // module
    Int_t iX = relid[2]-1;    // cluster X-coord
    Int_t iZ = relid[3]-1;    // cluster Z-coord
	Float_t e1ii = clu1.E(); // cluster energy after correction
    
    for (Int_t jClu=iClu; jClu<kNumberOfPhosClusters; jClu++) {
      AliAODCaloCluster *c2 = (AliAODCaloCluster *) caloClustersArr->At(jClu);
      if(!c2->IsPHOSCluster())   continue; // EMCAL cluster!
      if(c2->IsEqual(c1)) continue;

      Float_t e2i = c2->E();
      if(e2i<fEmin) continue;
      
      AliPHOSAodCluster clu2(*c2);
      clu2.Recalibrate(fCalibData, phsCells);
      clu2.EvalAll(fLogWeight,vtx);
      clu2.EnergyCorrection(&pid) ;
        
      clu2.GetMomentum(p2,v);
	  Float_t e2ii = clu2.E();

      p12 = p1+p2;

      fHmgg->Fill(p12.M()*1000); 
      fHmpi0[mod1][iX][iZ]->Fill(p12.M()*1000);
      
	  if(DebugLevel() > 1) printf("AliAnalysisTaskPHOSPi0CalibSelection Mass in (mod%d,%d,%d): %.3f GeV  E1_i=%f E1_ii=%f  E2_i=%f E2_ii=%f\n",
	     mod1,iX,iZ,p12.M(),e1i,e1ii,e2i,e2ii);
    }
    
  } // end of loop over PHOS clusters
  
  delete caloClustersArr;
  PostData(1,fOutputContainer);
}

void AliAnalysisTaskPHOSPi0CalibSelection::MaxEnergyCellPos(AliAODCaloCells *cells, AliAODCaloCluster* clu, Int_t& maxId)
{
  //For a given CaloCluster calculates the absId of the cell 
  //with maximum energy deposit.
  
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

void AliAnalysisTaskPHOSPi0CalibSelection::SetCalibCorrections(AliPHOSCalibData* cdata)
{
  //Set new correction factors (~1) to calibration coefficients, delete previous.

   if(fCalibData) delete fCalibData;
   fCalibData = cdata;
	
}

