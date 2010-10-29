/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*****************************************************
//   Class AliCentralitySelectionTask
//   Class to analyze determine centrality            
//   author: Alberica Toia
//*****************************************************

#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <iostream>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliESDCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliCentralitySelectionTask.h"

ClassImp(AliCentralitySelectionTask)


//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask():
AliAnalysisTaskSE(),
  fDebug(0),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fFile(0),
  fFile2(0),
  fCentfilename(0),
  fCentfilename2(0),
  fCentV0M(0),
  fCentFMD(0),
  fCentTRK(0),
  fCentTKL(0),
  fCentCL0(0),
  fCentV0MvsFMD(0),
  fCentTKLvsV0M(0),
  fCentZEMvsZDC(0),
  fHtempV0M(0),
  fHtempFMD(0),
  fHtempTRK(0),
  fHtempTKL(0),
  fHtempCL0(0),
  fHtempV0MvsFMD(0),
  fHtempTKLvsV0M(0),
  fHtempZEMvsZDC(0)
{   
  // Default constructor
}   

//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fDebug(0),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fFile(0),
  fFile2(0),
  fCentfilename(0),
  fCentfilename2(0),
  fCentV0M(0),
  fCentFMD(0),
  fCentTRK(0),
  fCentTKL(0),
  fCentCL0(0),
  fCentV0MvsFMD(0),
  fCentTKLvsV0M(0),
  fCentZEMvsZDC(0),
  fHtempV0M(0),
  fHtempFMD(0),
  fHtempTRK(0),
  fHtempTKL(0),
  fHtempCL0(0),
  fHtempV0MvsFMD(0),
  fHtempTKLvsV0M(0),
  fHtempZEMvsZDC(0)
{
  // Default constructor
}

//________________________________________________________________________
AliCentralitySelectionTask& AliCentralitySelectionTask::operator=(const AliCentralitySelectionTask& c)
{
  // Assignment operator
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
  }
  return *this;
}

//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask(const AliCentralitySelectionTask& ana):
  AliAnalysisTaskSE(ana),
  fDebug(ana.fDebug),	  
  fAnalysisInput(ana.fDebug),
  fIsMCInput(ana.fIsMCInput),
  fFile(ana.fFile),
  fFile2(ana.fFile2),
  fCentfilename(ana.fCentfilename),
  fCentfilename2(ana.fCentfilename2),
  fCentV0M(ana.fCentV0M),
  fCentFMD(ana.fCentFMD),
  fCentTRK(ana.fCentTRK),
  fCentTKL(ana.fCentTKL),
  fCentCL0(ana.fCentCL0),
  fCentV0MvsFMD(ana.fCentV0MvsFMD),
  fCentTKLvsV0M(ana.fCentTKLvsV0M),
  fCentZEMvsZDC(ana.fCentZEMvsZDC),
  fHtempV0M(ana.fHtempV0M),
  fHtempFMD(ana.fHtempFMD),
  fHtempTRK(ana.fHtempTRK),
  fHtempTKL(ana.fHtempTKL),
  fHtempCL0(ana.fHtempCL0),
  fHtempV0MvsFMD(ana.fHtempV0MvsFMD),
  fHtempTKLvsV0M(ana.fHtempTKLvsV0M),
  fHtempZEMvsZDC(ana.fHtempZEMvsZDC)
{
  // Copy Constructor	
}
 
//________________________________________________________________________
 AliCentralitySelectionTask::~AliCentralitySelectionTask()
 {
   // Destructor
 }  

//________________________________________________________________________
void AliCentralitySelectionTask::UserCreateOutputObjects()
{  
  // Create the output containers
  if(fDebug>1) printf("AnalysisCentralitySelectionTask::UserCreateOutputObjects() \n");
}

//________________________________________________________________________
void AliCentralitySelectionTask::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliCentralitySelectionTask::UserExec() \n");
  
  Float_t  zncEnergy;               //  ZNC Energy
  Float_t  zpcEnergy;               //  ZPC Energy
  Float_t  znaEnergy;               //  ZNA Energy
  Float_t  zpaEnergy;               //  ZPA Energy
  Float_t  zem1Energy = 0.;         //  ZEM1 Energy
  Float_t  zem2Energy = 0.;         //  ZEM2 Energy
  
  Int_t    nTracks    = 0;          //  no. tracks
  Int_t    nTracklets = 0;          //  no. tracklets
  Int_t    nClusters[6];            //  no. clusters on 6 ITS layers
  Int_t    nChips[2];               //  no. chips on 2 SPD layers

  Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
  Float_t  multFMDA = 0;            //  multiplicity from FMD on detector A
  Float_t  multFMDC = 0;            //  multiplicity from FMD on detector C

  AliESDCentrality *esdCent = 0;

  if(fAnalysisInput.CompareTo("ESD")==0){

    AliVEvent* event = InputEvent();
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);

    esdCent = esd->GetCentrality();

    // ***** V0 info    
    AliESDVZERO* esdV0 = esd->GetVZEROData();
    multV0A=esdV0->GetMTotV0A();
    multV0C=esdV0->GetMTotV0C();
    
    // ***** CB info (tracklets, clusters, chips)
    nTracks    = event->GetNumberOfTracks();     

    const AliMultiplicity *mult = esd->GetMultiplicity();

    nTracklets = mult->GetNumberOfTracklets();

    for(Int_t ilay=0; ilay<6; ilay++){
      nClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
    }

    for(Int_t ilay=0; ilay<2; ilay++){
      nChips[ilay] = mult->GetNumberOfFiredChips(ilay);
    }
    

    // ***** FMD info
    AliESDFMD *fmd = esd->GetFMDData();
    Float_t totalMultA = 0;
    Float_t totalMultC = 0;
    const Float_t fFMDLowCut = 0.4;
    
    for(UShort_t det=1;det<=3;det++) {
      Int_t nRings = (det==1 ? 1 : 2);
      for (UShort_t ir = 0; ir < nRings; ir++) {	  
	Char_t   ring = (ir == 0 ? 'I' : 'O');
	UShort_t nsec = (ir == 0 ? 20  : 40);
	UShort_t nstr = (ir == 0 ? 512 : 256);
	for(UShort_t sec =0; sec < nsec;  sec++)  {
	  for(UShort_t strip = 0; strip < nstr; strip++) {
	    
	    Float_t FMDmult = fmd->Multiplicity(det,ring,sec,strip);
	    if(FMDmult == 0 || FMDmult == AliESDFMD::kInvalidMult) continue;
	    
	    Float_t nParticles=0;
	    
	    if(FMDmult > fFMDLowCut) {
	      nParticles = 1.;
	    }
	    
	    if (det<3) totalMultA = totalMultA + nParticles;
	    else totalMultC = totalMultC + nParticles;
	    
	  }
	}
      }
    }
    multFMDA = totalMultA;
    multFMDC = totalMultC;
    
    // ***** ZDC info
    AliESDZDC *esdZDC = esd->GetESDZDC();
    zncEnergy = (Float_t) (esdZDC->GetZDCN1Energy());
    zpcEnergy = (Float_t) (esdZDC->GetZDCP1Energy());
    znaEnergy = (Float_t) (esdZDC->GetZDCN2Energy());
    zpaEnergy = (Float_t) (esdZDC->GetZDCP2Energy());
    zem1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0));
    zem2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1));
    
  }   
  else if(fAnalysisInput.CompareTo("AOD")==0){
    //AliAODEvent *aod =  dynamic_cast<AliAODEvent*> (InputEvent());
    // to be implemented
    printf("  AOD analysis not yet implemented!!!\n\n");
    return;
  }

  // ***** Centrality Selection
  fCentV0M = fHtempV0M->GetBinContent(fHtempV0M->FindBin((multV0A+multV0C)));
  fCentFMD = fHtempFMD->GetBinContent(fHtempFMD->FindBin((multFMDA+multFMDC)));
  fCentTRK = fHtempTRK->GetBinContent(fHtempTRK->FindBin(nTracks));
  fCentTKL = fHtempTKL->GetBinContent(fHtempTKL->FindBin(nTracklets));
  fCentCL0 = fHtempCL0->GetBinContent(fHtempCL0->FindBin(nClusters[0]));
  
  fCentV0MvsFMD = fHtempV0MvsFMD->GetBinContent(fHtempV0MvsFMD->FindBin((multV0A+multV0C)));
  fCentTKLvsV0M = fHtempTKLvsV0M->GetBinContent(fHtempTKLvsV0M->FindBin(nTracklets));
  fCentZEMvsZDC = fHtempZEMvsZDC->GetBinContent(fHtempZEMvsZDC->FindBin((zem1Energy+zem2Energy)/1000.));
  
  esdCent->SetCentralityV0M(fCentV0M);
  esdCent->SetCentralityFMD(fCentFMD);
  esdCent->SetCentralityTRK(fCentTRK);
  esdCent->SetCentralityTKL(fCentTKL);
  esdCent->SetCentralityCL0(fCentCL0);
  esdCent->SetCentralityV0MvsFMD(fCentV0MvsFMD);
  esdCent->SetCentralityTKLvsV0M(fCentTKLvsV0M);
  esdCent->SetCentralityZEMvsZDC(fCentZEMvsZDC);
}

//________________________________________________________________________
void AliCentralitySelectionTask::ReadCentralityHistos() 
{
  fFile  = new TFile(fCentfilename);
  fHtempV0M  = (TH1D*) (fFile->Get("hmultV0_percentile")); 
  fHtempFMD  = (TH1D*) (fFile->Get("hmultFMD_percentile")); 
  fHtempTRK  = (TH1D*) (fFile->Get("hNtracks_percentile")); 
  fHtempTKL  = (TH1D*) (fFile->Get("hNtracklets_percentile")); 
  fHtempCL0  = (TH1D*) (fFile->Get("hNclusters0_percentile")); 
}  
  
void AliCentralitySelectionTask::ReadCentralityHistos2() 
{
  fFile2  = new TFile(fCentfilename2);
  fHtempV0MvsFMD  = (TH1D*) (fFile2->Get("hmultV0vsmultFMD_all_percentile")); 
  fHtempTKLvsV0M  = (TH1D*) (fFile2->Get("hNtrackletsvsmultV0_all_percentile")); 
  fHtempZEMvsZDC  = (TH1D*) (fFile2->Get("hEzemvsEzdc_all_percentile"));  
}

void AliCentralitySelectionTask::SetPercentileFile(TString filename) 
{
  fCentfilename = filename;
  ReadCentralityHistos();
}

void AliCentralitySelectionTask::SetPercentileFile2(TString filename) 
{
  fCentfilename2 = filename;
  ReadCentralityHistos2();
}

//________________________________________________________________________
void AliCentralitySelectionTask::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fFile->Close();  
  fFile2->Close();  
}
//________________________________________________________________________

