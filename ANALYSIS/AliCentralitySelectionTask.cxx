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

#include "AliCentralitySelectionTask.h"

#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TFile.h>
#include <TObjString.h>
#include <TString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <iostream>

#include "AliAnalysisManager.h"
#include "AliHeader.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliESDFMD.h"
#include "AliESDVZERO.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODTracklets.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliESDUtils.h"

ClassImp(AliCentralitySelectionTask)


//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask():
AliAnalysisTaskSE(),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fPass(0),
  fFile(0),
  fFile2(0),
  fCurrentRun(-1),
  fRunNo(-1),
  fLowRunN(0),
  fHighRunN(0),
  fUseScaling(0),
  fUseCleaning(0),
  fTrackCuts(0),
  fZVCut(10),
  fOutliersCut(5),
  fQuality(999),
  fCentV0M(0),
  fCentFMD(0),
  fCentTRK(0),
  fCentTKL(0),
  fCentCL0(0),
  fCentCL1(0),
  fCentV0MvsFMD(0),
  fCentTKLvsV0M(0),
  fCentZEMvsZDC(0),
  fHtempV0M(0),
  fHtempFMD(0),
  fHtempTRK(0),
  fHtempTKL(0),
  fHtempCL0(0),
  fHtempCL1(0),
  fHtempV0MvsFMD(0),
  fHtempTKLvsV0M(0),
  fHtempZEMvsZDC(0),
  fOutputList(0),
  fHOutCentV0M     (0),
  fHOutCentFMD     (0),
  fHOutCentTRK     (0),
  fHOutCentTKL     (0),
  fHOutCentCL0     (0),
  fHOutCentCL1     (0),
  fHOutCentV0MvsFMD(0),
  fHOutCentTKLvsV0M(0),
  fHOutCentZEMvsZDC(0),
  fHOutCentV0MvsCentCL1(0),
  fHOutCentV0MvsCentTRK(0),
  fHOutCentTRKvsCentCL1(0),
  fHOutCentV0MvsCentZDC(0),
  fHOutMultV0M(0),
  fHOutMultV0R(0),
  fHOutMultFMD(0),
  fHOutMultTRK(0),
  fHOutMultTKL(0),
  fHOutMultCL0(0),
  fHOutMultCL1(0),
  fHOutMultV0MvsZDN(0),
  fHOutMultZEMvsZDN(0),
  fHOutMultV0MvsZDC(0),
  fHOutMultZEMvsZDC(0),
  fHOutMultZEMvsZDCw(0),
  fHOutMultV0MvsCL1(0),
  fHOutMultV0MvsTRK(0),
  fHOutMultTRKvsCL1(0),
  fHOutCentV0Mqual1(0),
  fHOutCentTRKqual1(0),
  fHOutCentCL1qual1(0),
  fHOutMultV0MvsCL1qual1(0),
  fHOutMultV0MvsTRKqual1(0),
  fHOutMultTRKvsCL1qual1(0),
  fHOutCentV0Mqual2(0),
  fHOutCentTRKqual2(0),
  fHOutCentCL1qual2(0),
  fHOutMultV0MvsCL1qual2(0),
  fHOutMultV0MvsTRKqual2(0),
  fHOutMultTRKvsCL1qual2(0),
  fHOutQuality(0),
  fHOutVertex(0)
{   
  // Default constructor
  AliInfo("Centrality Selection enabled.");
  fLowRunN =136851;
  fHighRunN=139517;

  for (Int_t i=0; i < 2667; i++) {
    fV0MScaleFactor[i]=0.0;
    fSPDScaleFactor[i]=0.0;
    fTPCScaleFactor[i]=0.0;
    fV0MScaleFactorMC[i]=0.0;
  }
  fUseScaling=kTRUE;
  fUseCleaning=kTRUE;
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,AliESDZDC.,AliESDFMD.,AliESDVZERO."
               ",SPDVertex.,TPCVertex.,PrimaryVertex.,AliMultiplicity.,Tracks ";
}   

//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fPass(0),
  fFile(0),
  fFile2(0),
  fCurrentRun(-1),
  fRunNo(-1),
  fLowRunN(0),
  fHighRunN(0),
  fUseScaling(0),
  fUseCleaning(0),
  fTrackCuts(0),
  fZVCut(10),
  fOutliersCut(5),
  fQuality(999),
  fCentV0M(0),
  fCentFMD(0),
  fCentTRK(0),
  fCentTKL(0),
  fCentCL0(0),
  fCentCL1(0),
  fCentV0MvsFMD(0),
  fCentTKLvsV0M(0),
  fCentZEMvsZDC(0),
  fHtempV0M(0),
  fHtempFMD(0),
  fHtempTRK(0),
  fHtempTKL(0),
  fHtempCL0(0),
  fHtempCL1(0),
  fHtempV0MvsFMD(0),
  fHtempTKLvsV0M(0),
  fHtempZEMvsZDC(0),
  fOutputList(0),
  fHOutCentV0M     (0),
  fHOutCentFMD     (0),
  fHOutCentTRK     (0),
  fHOutCentTKL     (0),
  fHOutCentCL0     (0),
  fHOutCentCL1     (0),
  fHOutCentV0MvsFMD(0),
  fHOutCentTKLvsV0M(0),
  fHOutCentZEMvsZDC(0),
  fHOutCentV0MvsCentCL1(0),
  fHOutCentV0MvsCentTRK(0),
  fHOutCentTRKvsCentCL1(0),
  fHOutCentV0MvsCentZDC(0),
  fHOutMultV0M(0),
  fHOutMultV0R(0),
  fHOutMultFMD(0),
  fHOutMultTRK(0),
  fHOutMultTKL(0),
  fHOutMultCL0(0),
  fHOutMultCL1(0),
  fHOutMultV0MvsZDN(0),
  fHOutMultZEMvsZDN(0),
  fHOutMultV0MvsZDC(0),
  fHOutMultZEMvsZDC(0),
  fHOutMultZEMvsZDCw(0),
  fHOutMultV0MvsCL1(0),
  fHOutMultV0MvsTRK(0),
  fHOutMultTRKvsCL1(0),
  fHOutCentV0Mqual1(0),
  fHOutCentTRKqual1(0),
  fHOutCentCL1qual1(0),
  fHOutMultV0MvsCL1qual1(0),
  fHOutMultV0MvsTRKqual1(0),
  fHOutMultTRKvsCL1qual1(0),
  fHOutCentV0Mqual2(0),
  fHOutCentTRKqual2(0),
  fHOutCentCL1qual2(0),
  fHOutMultV0MvsCL1qual2(0),
  fHOutMultV0MvsTRKqual2(0),
  fHOutMultTRKvsCL1qual2(0),
  fHOutQuality(0),
  fHOutVertex(0)
{
    // Default constructor
    fLowRunN =136851;
    fHighRunN=139517;

  AliInfo("Centrality Selection enabled.");
  DefineOutput(1, TList::Class());
  for (Int_t i=0; i<2667; i++) {
    fV0MScaleFactor[i]=0.0;
    fSPDScaleFactor[i]=0.0;
    fTPCScaleFactor[i]=0.0;
    fV0MScaleFactorMC[i]=0.0;
  }
  fUseScaling=kTRUE;
  fUseCleaning=kTRUE;
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,AliESDZDC.,AliESDFMD.,AliESDVZERO."
               ",SPDVertex.,TPCVertex.,PrimaryVertex.,AliMultiplicity.,Tracks ";
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
  fAnalysisInput(ana.fDebug),
  fIsMCInput(ana.fIsMCInput),
  fPass(ana.fPass),
  fFile(ana.fFile),
  fFile2(ana.fFile2),
  fCurrentRun(ana.fCurrentRun),
  fRunNo(ana.fRunNo),
  fLowRunN(ana.fLowRunN),
  fHighRunN(ana.fHighRunN),
  fUseScaling(ana.fUseScaling),
  fUseCleaning(ana.fUseCleaning),
  fTrackCuts(ana.fTrackCuts),
  fZVCut(ana.fZVCut),
  fOutliersCut(ana.fOutliersCut),
  fQuality(ana.fQuality),
  fCentV0M(ana.fCentV0M),
  fCentFMD(ana.fCentFMD),
  fCentTRK(ana.fCentTRK),
  fCentTKL(ana.fCentTKL),
  fCentCL0(ana.fCentCL0),
  fCentCL1(ana.fCentCL1),
  fCentV0MvsFMD(ana.fCentV0MvsFMD),
  fCentTKLvsV0M(ana.fCentTKLvsV0M),
  fCentZEMvsZDC(ana.fCentZEMvsZDC),
  fHtempV0M(ana.fHtempV0M),
  fHtempFMD(ana.fHtempFMD),
  fHtempTRK(ana.fHtempTRK),
  fHtempTKL(ana.fHtempTKL),
  fHtempCL0(ana.fHtempCL0),
  fHtempCL1(ana.fHtempCL1),
  fHtempV0MvsFMD(ana.fHtempV0MvsFMD),
  fHtempTKLvsV0M(ana.fHtempTKLvsV0M),
  fHtempZEMvsZDC(ana.fHtempZEMvsZDC),
  fOutputList(ana.fOutputList),
  fHOutCentV0M     (ana.fHOutCentV0M     ),
  fHOutCentFMD     (ana.fHOutCentFMD     ),
  fHOutCentTRK     (ana.fHOutCentTRK     ),
  fHOutCentTKL     (ana.fHOutCentTKL     ),
  fHOutCentCL0     (ana.fHOutCentCL0     ),
  fHOutCentCL1     (ana.fHOutCentCL1     ),
  fHOutCentV0MvsFMD(ana.fHOutCentV0MvsFMD),
  fHOutCentTKLvsV0M(ana.fHOutCentTKLvsV0M),
  fHOutCentZEMvsZDC(ana.fHOutCentZEMvsZDC),
  fHOutCentV0MvsCentCL1(ana.fHOutCentV0MvsCentCL1),
  fHOutCentV0MvsCentTRK(ana.fHOutCentV0MvsCentTRK),
  fHOutCentTRKvsCentCL1(ana.fHOutCentTRKvsCentCL1),
  fHOutCentV0MvsCentZDC(ana.fHOutCentV0MvsCentZDC),
  fHOutMultV0M(ana.fHOutMultV0M),
  fHOutMultV0R(ana.fHOutMultV0R),
  fHOutMultFMD(ana.fHOutMultFMD),
  fHOutMultTRK(ana.fHOutMultTRK),
  fHOutMultTKL(ana.fHOutMultTKL),
  fHOutMultCL0(ana.fHOutMultCL0),
  fHOutMultCL1(ana.fHOutMultCL1),
  fHOutMultV0MvsZDN(ana.fHOutMultV0MvsZDN),
  fHOutMultZEMvsZDN(ana.fHOutMultZEMvsZDN),
  fHOutMultV0MvsZDC(ana.fHOutMultV0MvsZDC),
  fHOutMultZEMvsZDC(ana.fHOutMultZEMvsZDC),
  fHOutMultZEMvsZDCw(ana.fHOutMultZEMvsZDCw),
  fHOutMultV0MvsCL1(ana.fHOutMultV0MvsCL1),
  fHOutMultV0MvsTRK(ana.fHOutMultV0MvsTRK),
  fHOutMultTRKvsCL1(ana.fHOutMultTRKvsCL1),
  fHOutCentV0Mqual1(ana.fHOutCentV0Mqual1),
  fHOutCentTRKqual1(ana.fHOutCentTRKqual1),
  fHOutCentCL1qual1(ana.fHOutCentCL1qual1),
  fHOutMultV0MvsCL1qual1(ana.fHOutMultV0MvsCL1qual1),
  fHOutMultV0MvsTRKqual1(ana.fHOutMultV0MvsTRKqual1),
  fHOutMultTRKvsCL1qual1(ana.fHOutMultTRKvsCL1qual1),
  fHOutCentV0Mqual2(ana.fHOutCentV0Mqual2),
  fHOutCentTRKqual2(ana.fHOutCentTRKqual2),
  fHOutCentCL1qual2(ana.fHOutCentCL1qual2),
  fHOutMultV0MvsCL1qual2(ana.fHOutMultV0MvsCL1qual2),
  fHOutMultV0MvsTRKqual2(ana.fHOutMultV0MvsTRKqual2),
  fHOutMultTRKvsCL1qual2(ana.fHOutMultTRKvsCL1qual2),
  fHOutQuality(ana.fHOutQuality),
  fHOutVertex(ana.fHOutVertex)
{
  // Copy Constructor	
    for (Int_t i=0; i<2667; i++) {
	fV0MScaleFactor[i]=0.0;
	fSPDScaleFactor[i]=0.0;
	fTPCScaleFactor[i]=0.0;
	fV0MScaleFactorMC[i]=0.0;
    }

}
 
//________________________________________________________________________
AliCentralitySelectionTask::~AliCentralitySelectionTask()
{
  // Destructor  
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;
  if (fTrackCuts) delete fTrackCuts;
}  

//________________________________________________________________________
void AliCentralitySelectionTask::UserCreateOutputObjects()
{  
  // Create the output containers
  if(fDebug>1) printf("AnalysisCentralitySelectionTask::UserCreateOutputObjects() \n");
  AliLog::SetClassDebugLevel("AliCentralitySelectionTask", AliLog::kInfo);

  fOutputList = new TList();
  fOutputList->SetOwner();
  fHOutCentV0M     = new TH1F("fHOutCentV0M","fHOutCentV0M; Centrality V0",505,0,101);
  fHOutCentFMD     = new TH1F("fHOutCentFMD","fHOutCentFMD; Centrality FMD",505,0,101);
  fHOutCentTRK     = new TH1F("fHOutCentTRK","fHOutCentTRK; Centrality TPC",505,0,101);
  fHOutCentTKL     = new TH1F("fHOutCentTKL","fHOutCentTKL; Centrality tracklets",505,0,101);
  fHOutCentCL0     = new TH1F("fHOutCentCL0","fHOutCentCL0; Centrality SPD inner",505,0,101);
  fHOutCentCL1     = new TH1F("fHOutCentCL1","fHOutCentCL1; Centrality SPD outer",505,0,101);
  fHOutCentV0MvsFMD= new TH1F("fHOutCentV0MvsFMD","fHOutCentV0MvsFMD; Centrality V0 vs FMD",505,0,101);
  fHOutCentTKLvsV0M= new TH1F("fHOutCentTKLvsV0M","fHOutCentTKLvsV0M; Centrality tracklets vs V0",505,0,101);
  fHOutCentZEMvsZDC= new TH1F("fHOutCentZEMvsZDC","fHOutCentZEMvsZDC; Centrality ZEM vs ZDC",505,0,101);
  fHOutCentV0MvsCentCL1= new TH2F("fHOutCentV0MvsCentCL1","fHOutCentV0MvsCentCL1; Cent V0; Cent SPD",505,0,101,505,0,101);
  fHOutCentV0MvsCentTRK= new TH2F("fHOutCentV0MvsCentTRK","fHOutCentV0MvsCentTRK; Cent V0; Cent TPC",505,0,101,505,0,101);
  fHOutCentTRKvsCentCL1= new TH2F("fHOutCentTRKvsCentCL1","fHOutCentTRKvsCentCL1; Cent TPC; Cent SPD",505,0,101,505,0,101);
  fHOutCentV0MvsCentZDC= new TH2F("fHOutCentV0MvsCentZDC","fHOutCentV0MvsCentZDC; Cent V0; Cent ZDC",505,0,101,505,0,101);

  fHOutMultV0M = new TH1F("fHOutMultV0M","fHOutMultV0M; Multiplicity V0",25000,0,25000);
  fHOutMultV0R = new TH1F("fHOutMultV0R","fHOutMultV0R; Multiplicity V0",30000,0,30000);
  fHOutMultFMD = new TH1F("fHOutMultFMD","fHOutMultFMD; Multiplicity FMD",24000,0,24000);
  fHOutMultTRK = new TH1F("fHOutMultTRK","fHOutMultTRK; Multiplicity TPC",4000,0,4000);
  fHOutMultTKL = new TH1F("fHOutMultTKL","fHOutMultTKL; Multiplicity tracklets",5000,0,5000);
  fHOutMultCL0 = new TH1F("fHOutMultCL0","fHOutMultCL0; Multiplicity SPD inner",7000,0,7000);
  fHOutMultCL1 = new TH1F("fHOutMultCL1","fHOutMultCL1; Multiplicity SPD outer",7000,0,7000);
  fHOutMultV0MvsZDN = new TH2F("fHOutMultV0MvsZDN","fHOutMultV0MvsZDN; Multiplicity V0; Energy ZDC-N",500,0,25000,500,0,180000);
  fHOutMultZEMvsZDN = new TH2F("fHOutMultZEMvsZDN","fHOutMultZEMvsZDN; Energy ZEM; Energy ZDC-N",500,0,2500,500,0,180000);
  fHOutMultV0MvsZDC = new TH2F("fHOutMultV0MvsZDC","fHOutMultV0MvsZDC; Multiplicity V0; Energy ZDC",500,0,25000,500,0,200000);
  fHOutMultZEMvsZDC = new TH2F("fHOutMultZEMvsZDC","fHOutMultZEMvsZDC; Energy ZEM; Energy ZDC",500,0,2500,500,0,200000);
  fHOutMultZEMvsZDCw = new TH2F("fHOutMultZEMvsZDCw","fHOutMultZEMvsZDCw; Energy ZEM; Energy ZDC (weigthed with V0 percentile)",500,0,2500,500,0,200000);
  fHOutMultV0MvsCL1 = new TH2F("fHOutMultV0MvsCL1","fHOutMultV0MvsCL1; Multiplicity V0; Multiplicity SPD outer",2500,0,25000,700,0,7000);
  fHOutMultV0MvsTRK = new TH2F("fHOutMultV0MvsTRK","fHOutMultV0MvsTRK; Multiplicity V0; Multiplicity TPC",2500,0,25000,400,0,4000);
  fHOutMultTRKvsCL1 = new TH2F("fHOutMultTRKvsCL1","fHOutMultTRKvsCL1; Multiplicity TPC; Multiplicity SPD outer",400,0,4000,700,0,7000);

  fHOutCentV0Mqual1 = new TH1F("fHOutCentV0M_qual1","fHOutCentV0M_qual1; Centrality V0",505,0,101);
  fHOutCentTRKqual1 = new TH1F("fHOutCentTRK_qual1","fHOutCentTRK_qual1; Centrality TPC",505,0,101);
  fHOutCentCL1qual1 = new TH1F("fHOutCentCL1_qual1","fHOutCentCL1_qual1; Centrality SPD outer",505,0,101);
  fHOutMultV0MvsCL1qual1 = new TH2F("fHOutMultV0MvsCL1_qual1","fHOutMultV0MvsCL1_qual1; Multiplicity V0; Multiplicity SPD outer",2500,0,25000,700,0,7000);
  fHOutMultV0MvsTRKqual1 = new TH2F("fHOutMultV0MvsTRK_qual1","fHOutMultV0MvsTRK_qual1; Multiplicity V0; Multiplicity TPC",2500,0,25000,400,0,4000);
  fHOutMultTRKvsCL1qual1 = new TH2F("fHOutMultTRKvsCL1_qual1","fHOutMultTRKvsCL1_qual1; Multiplicity TPC; Multiplicity SPD outer",400,0,4000,700,0,7000);

  fHOutCentV0Mqual2 = new TH1F("fHOutCentV0M_qual2","fHOutCentV0M_qual2; Centrality V0",505,0,101);
  fHOutCentTRKqual2 = new TH1F("fHOutCentTRK_qual2","fHOutCentTRK_qual2; Centrality TPC",505,0,101);
  fHOutCentCL1qual2 = new TH1F("fHOutCentCL1_qual2","fHOutCentCL1_qual2; Centrality SPD outer",505,0,101);
  fHOutMultV0MvsCL1qual2 = new TH2F("fHOutMultV0MvsCL1_qual2","fHOutMultV0MvsCL1_qual2; Multiplicity V0; Multiplicity SPD outer",2500,0,25000,700,0,7000);
  fHOutMultV0MvsTRKqual2 = new TH2F("fHOutMultV0MvsTRK_qual2","fHOutMultV0MvsTRK_qual2; Multiplicity V0; Multiplicity TPC",2500,0,25000,400,0,4000);
  fHOutMultTRKvsCL1qual2 = new TH2F("fHOutMultTRKvsCL1_qual2","fHOutMultTRKvsCL1_qual2; Multiplicity TPC; Multiplicity SPD outer",400,0,4000,700,0,7000);

  fHOutQuality = new TH1F("fHOutQuality", "fHOutQuality", 100,-0.5,99.5);
  fHOutVertex  = new TH1F("fHOutVertex", "fHOutVertex", 100,-20,20);

  fOutputList->Add(  fHOutCentV0M     );
  fOutputList->Add(  fHOutCentFMD     );
  fOutputList->Add(  fHOutCentTRK     );
  fOutputList->Add(  fHOutCentTKL     );
  fOutputList->Add(  fHOutCentCL0     );
  fOutputList->Add(  fHOutCentCL1     );
  fOutputList->Add(  fHOutCentV0MvsFMD);
  fOutputList->Add(  fHOutCentTKLvsV0M);
  fOutputList->Add(  fHOutCentZEMvsZDC);
  fOutputList->Add(  fHOutCentV0MvsCentCL1);
  fOutputList->Add(  fHOutCentV0MvsCentTRK);
  fOutputList->Add(  fHOutCentTRKvsCentCL1);
  fOutputList->Add(  fHOutCentV0MvsCentZDC);
  fOutputList->Add(  fHOutMultV0M); 
  fOutputList->Add(  fHOutMultV0R); 
  fOutputList->Add(  fHOutMultFMD); 
  fOutputList->Add(  fHOutMultTRK); 
  fOutputList->Add(  fHOutMultTKL); 
  fOutputList->Add(  fHOutMultCL0); 
  fOutputList->Add(  fHOutMultCL1); 
  fOutputList->Add(  fHOutMultV0MvsZDN);
  fOutputList->Add(  fHOutMultZEMvsZDN);
  fOutputList->Add(  fHOutMultV0MvsZDC);
  fOutputList->Add(  fHOutMultZEMvsZDC);
  fOutputList->Add(  fHOutMultZEMvsZDCw);
  fOutputList->Add(  fHOutMultV0MvsCL1);
  fOutputList->Add(  fHOutMultV0MvsTRK);
  fOutputList->Add(  fHOutMultTRKvsCL1);
  fOutputList->Add(  fHOutCentV0Mqual1 );
  fOutputList->Add(  fHOutCentTRKqual1 );
  fOutputList->Add(  fHOutCentCL1qual1 );                   
  fOutputList->Add(  fHOutMultV0MvsCL1qual1);
  fOutputList->Add(  fHOutMultV0MvsTRKqual1);
  fOutputList->Add(  fHOutMultTRKvsCL1qual1);
  fOutputList->Add(  fHOutCentV0Mqual2 );
  fOutputList->Add(  fHOutCentTRKqual2 );
  fOutputList->Add(  fHOutCentCL1qual2 );
  fOutputList->Add(  fHOutMultV0MvsCL1qual2);
  fOutputList->Add(  fHOutMultV0MvsTRKqual2);
  fOutputList->Add(  fHOutMultTRKvsCL1qual2);
  fOutputList->Add(  fHOutQuality );
  fOutputList->Add(  fHOutVertex );


  fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

  PostData(1, fOutputList); 

  MyInitScaleFactor();
  if (fIsMCInput) MyInitScaleFactorMC();

  if (fPass==0) AliFatal("Which pass are you analyzing? You should set it via taskCentrality->SetPass(N)");
}

//________________________________________________________________________
void AliCentralitySelectionTask::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliCentralitySelectionTask::UserExec() \n");
  
  Float_t  zncEnergy = 0.;          //  ZNC Energy
  Float_t  zpcEnergy = 0.;          //  ZPC Energy
  Float_t  znaEnergy = 0.;          //  ZNA Energy
  Float_t  zpaEnergy = 0.;          //  ZPA Energy
  Float_t  zem1Energy = 0.;         //  ZEM1 Energy
  Float_t  zem2Energy = 0.;         //  ZEM2 Energy
  Bool_t   zdcEnergyCal = kFALSE;   // if zdc is calibrated (in pass2)

  Int_t    nTracks = 0;             //  no. tracks
  Int_t    nTracklets = 0;          //  no. tracklets
  Int_t    nClusters[6] = {0};      //  no. clusters on 6 ITS layers
  Int_t    nChips[2];               //  no. chips on 2 SPD layers
  Float_t  spdCorr =0;              //  corrected spd2 multiplicity

  Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
  Float_t  multFMDA = 0;            //  multiplicity from FMD on detector A
  Float_t  multFMDC = 0;            //  multiplicity from FMD on detector C

  Short_t v0Corr = 0;               // corrected V0 multiplicity
  Short_t v0CorrResc = 0;           // corrected and rescaled V0 multiplicity

  Float_t zvtx =0;                  // z-vertex SPD
  Int_t zvtxNcont =0;               // contributors to z-vertex SPD

 
  AliCentrality *esdCent = 0;

  if(fAnalysisInput.CompareTo("ESD")==0){

    AliVEvent* event = InputEvent();
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (!esd) {
	AliError("No ESD Event");
	return;
    }

    LoadBranches();
    
    if (fRunNo<=0) {
	if (SetupRun(esd)<0) {
	    AliError("Centrality File not available for this run");
	    return;
	}
    }

    esdCent = esd->GetCentrality();

    // ***** V0 info    
    AliESDVZERO* esdV0 = esd->GetVZEROData();
    multV0A=esdV0->GetMTotV0A();
    multV0C=esdV0->GetMTotV0C();

    Float_t v0CorrR;
    v0Corr = (Short_t)AliESDUtils::GetCorrV0(esd,v0CorrR);
    v0CorrResc = (Short_t)v0CorrR;

    // ***** Vertex Info
    const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
    zvtx        = vtxESD->GetZ(); 
    zvtxNcont   = vtxESD->GetNContributors();

    // ***** CB info (tracklets, clusters, chips)
    //nTracks    = event->GetNumberOfTracks();     
    nTracks    = fTrackCuts ? (Short_t)fTrackCuts->GetReferenceMultiplicity(esd,kTRUE):-1;

    const AliMultiplicity *mult = esd->GetMultiplicity();

    nTracklets = mult->GetNumberOfTracklets();

    for(Int_t ilay=0; ilay<6; ilay++){
      nClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
    }

    for(Int_t ilay=0; ilay<2; ilay++){
      nChips[ilay] = mult->GetNumberOfFiredChips(ilay);
    }

    spdCorr = AliESDUtils::GetCorrSPD2(nClusters[1],zvtx);    

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
	    
	    Float_t fmdMult = fmd->Multiplicity(det,ring,sec,strip);
	    if(fmdMult == 0 || fmdMult == AliESDFMD::kInvalidMult) continue;
	    
	    Float_t nParticles=0;
	    
	    if(fmdMult > fFMDLowCut) {
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
    zdcEnergyCal = esdZDC->AliESDZDC::TestBit(AliESDZDC::kEnergyCalibratedSignal);
    if (zdcEnergyCal) {      
      zncEnergy = (Float_t) (esdZDC->GetZDCN1Energy());
      zpcEnergy = (Float_t) (esdZDC->GetZDCP1Energy());
      znaEnergy = (Float_t) (esdZDC->GetZDCN2Energy());
      zpaEnergy = (Float_t) (esdZDC->GetZDCP2Energy());
    } else {
      zncEnergy = (Float_t) (esdZDC->GetZDCN1Energy())/8.;
      zpcEnergy = (Float_t) (esdZDC->GetZDCP1Energy())/8.;
      znaEnergy = (Float_t) (esdZDC->GetZDCN2Energy())/8.;
      zpaEnergy = (Float_t) (esdZDC->GetZDCP2Energy())/8.;
    }
    zem1Energy = (Float_t) (esdZDC->GetZDCEMEnergy(0))/8.;
    zem2Energy = (Float_t) (esdZDC->GetZDCEMEnergy(1))/8.;
  
  }   

  else if(fAnalysisInput.CompareTo("AOD")==0){
    //AliAODEvent *aod =  dynamic_cast<AliAODEvent*> (InputEvent());
    // to be implemented
    printf("  AOD analysis not yet implemented!!!\n\n");
    return;
  } 

  // ***** Scaling
  // ***** Scaling for pass2 
  if (fPass==2) {
    fUseScaling=kFALSE;
  }
  // ***** Scaling for MC
  if (fIsMCInput) {
    fUseScaling=kFALSE;
    Float_t tempScalefactorV0M = MyGetScaleFactorMC(fCurrentRun);
    v0Corr  = Short_t((multV0A+multV0C)  * tempScalefactorV0M);
  }
  // ***** Scaling for Data
  if (fUseScaling) {
    Float_t tempScalefactorV0M = MyGetScaleFactor(fCurrentRun,0);
    Float_t tempScalefactorSPD = MyGetScaleFactor(fCurrentRun,1);
    Float_t tempScalefactorTPC = MyGetScaleFactor(fCurrentRun,2);
    v0Corr  = Short_t(v0Corr / tempScalefactorV0M);
    spdCorr = spdCorr / tempScalefactorSPD;
    nTracks = Int_t(nTracks / tempScalefactorTPC);
  }

  // ***** Centrality Selection
  if(fHtempV0M) fCentV0M = fHtempV0M->GetBinContent(fHtempV0M->FindBin((v0Corr)));
  if(fHtempFMD) fCentFMD = fHtempFMD->GetBinContent(fHtempFMD->FindBin((multFMDA+multFMDC)));
  if(fHtempTRK) fCentTRK = fHtempTRK->GetBinContent(fHtempTRK->FindBin(nTracks));
  if(fHtempTKL) fCentTKL = fHtempTKL->GetBinContent(fHtempTKL->FindBin(nTracklets));
  if(fHtempCL0) fCentCL0 = fHtempCL0->GetBinContent(fHtempCL0->FindBin(nClusters[0]));
  if(fHtempCL1) fCentCL1 = fHtempCL1->GetBinContent(fHtempCL1->FindBin(spdCorr));
  
  if(fHtempV0MvsFMD) fCentV0MvsFMD = fHtempV0MvsFMD->GetBinContent(fHtempV0MvsFMD->FindBin((multV0A+multV0C)));
  if(fHtempTKLvsV0M) fCentTKLvsV0M = fHtempTKLvsV0M->GetBinContent(fHtempTKLvsV0M->FindBin(nTracklets));
  if(fHtempZEMvsZDC) fCentZEMvsZDC = fHtempZEMvsZDC->GetBinContent(fHtempZEMvsZDC->FindBin(zem1Energy+zem2Energy,zncEnergy+znaEnergy+zpcEnergy+zpaEnergy));

  // ***** Cleaning
  if (fUseCleaning) {
      fQuality=0;
      fZVCut=10;
      fOutliersCut=6;
      
      // ***** vertex
      if (TMath::Abs(zvtx)>fZVCut || zvtxNcont<1) fQuality += 1;   

      // ***** outliers, skip in case of MC input
      if (!fIsMCInput) {
	// **** V0 vs SPD
	if (IsOutlierV0MSPD(spdCorr, v0Corr, int(fCentV0M))) fQuality  += 2;
	// ***** V0 vs TPC
	if (IsOutlierV0MTPC(nTracks, v0Corr, int(fCentV0M))) fQuality  += 4;
	// ***** V0 vs ZDC
	if (IsOutlierV0MZDC((zncEnergy+znaEnergy+zpcEnergy+zpaEnergy), v0Corr) &&
	    (zdcEnergyCal==kFALSE) && !(fIsMCInput)) fQuality  += 8;
	if (IsOutlierV0MZDCECal((zncEnergy+znaEnergy+zpcEnergy+zpaEnergy), v0Corr) &&
	    ((zdcEnergyCal==kTRUE) || (fIsMCInput))) fQuality  += 8;
      }
  } else {
      fQuality = 0;
  }

        
  if (esdCent) {
      esdCent->SetQuality(fQuality);
      esdCent->SetCentralityV0M(fCentV0M);
      esdCent->SetCentralityFMD(fCentFMD);
      esdCent->SetCentralityTRK(fCentTRK);
      esdCent->SetCentralityTKL(fCentTKL);
      esdCent->SetCentralityCL0(fCentCL0);
      esdCent->SetCentralityCL1(fCentCL1);
      esdCent->SetCentralityV0MvsFMD(fCentV0MvsFMD);
      esdCent->SetCentralityTKLvsV0M(fCentTKLvsV0M);
      esdCent->SetCentralityZEMvsZDC(fCentZEMvsZDC);
  }

  fHOutQuality->Fill(fQuality);
  fHOutVertex->Fill(zvtx);
  
  if (fQuality==0) {  
    fHOutCentV0M->Fill(fCentV0M);
    fHOutCentFMD->Fill(fCentFMD);
    fHOutCentTRK->Fill(fCentTRK);
    fHOutCentTKL->Fill(fCentTKL);
    fHOutCentCL0->Fill(fCentCL0);
    fHOutCentCL1->Fill(fCentCL1);
    fHOutCentV0MvsFMD->Fill(fCentV0MvsFMD);
    fHOutCentTKLvsV0M->Fill(fCentTKLvsV0M);
    fHOutCentZEMvsZDC->Fill(fCentZEMvsZDC);
    fHOutCentV0MvsCentCL1->Fill(fCentV0M,fCentCL1);
    fHOutCentV0MvsCentTRK->Fill(fCentV0M,fCentTRK);
    fHOutCentTRKvsCentCL1->Fill(fCentTRK,fCentCL1);
    fHOutCentV0MvsCentZDC->Fill(fCentV0M,fCentZEMvsZDC);
    fHOutMultV0M->Fill(v0Corr);
    fHOutMultV0R->Fill(multV0A+multV0C);
    fHOutMultFMD->Fill((multFMDA+multFMDC));
    fHOutMultTRK->Fill(nTracks);
    fHOutMultTKL->Fill(nTracklets);
    fHOutMultCL0->Fill(nClusters[0]);
    fHOutMultCL1->Fill(spdCorr);
    fHOutMultV0MvsZDN->Fill(v0Corr,(zncEnergy+znaEnergy));
    fHOutMultZEMvsZDN->Fill((zem1Energy+zem2Energy),(zncEnergy+znaEnergy));
    fHOutMultV0MvsZDC->Fill(v0Corr,(zncEnergy+znaEnergy+zpcEnergy+zpaEnergy));
    fHOutMultZEMvsZDC->Fill((zem1Energy+zem2Energy),(zncEnergy+znaEnergy+zpcEnergy+zpaEnergy));
    fHOutMultZEMvsZDCw->Fill((zem1Energy+zem2Energy),(zncEnergy+znaEnergy+zpcEnergy+zpaEnergy),fCentV0M);
    fHOutMultV0MvsCL1->Fill(v0Corr,spdCorr);
    fHOutMultV0MvsTRK->Fill(v0Corr,nTracks);
    fHOutMultTRKvsCL1->Fill(nTracks,spdCorr);
  } else if (fQuality%2 == 0) {
    fHOutCentV0Mqual1->Fill(fCentV0M);
    fHOutCentTRKqual1->Fill(fCentTRK);
    fHOutCentCL1qual1->Fill(fCentCL1);
    fHOutMultV0MvsCL1qual1->Fill(v0Corr,spdCorr);
    fHOutMultV0MvsTRKqual1->Fill(v0Corr,nTracks);
    fHOutMultTRKvsCL1qual1->Fill(nTracks,spdCorr);
  } else {
    fHOutCentV0Mqual2->Fill(fCentV0M);
    fHOutCentTRKqual2->Fill(fCentTRK);
    fHOutCentCL1qual2->Fill(fCentCL1);
    fHOutMultV0MvsCL1qual2->Fill(v0Corr,spdCorr);
    fHOutMultV0MvsTRKqual2->Fill(v0Corr,nTracks);
    fHOutMultTRKvsCL1qual2->Fill(nTracks,spdCorr);
  }

  PostData(1, fOutputList); 
}

//________________________________________________________________________
void AliCentralitySelectionTask::ReadCentralityHistos(TString fCentfilename) 
{
  //  Read centrality histograms
  TDirectory *owd = gDirectory;
  // Check if the file is present
  TString path = gSystem->ExpandPathName(fCentfilename.Data());
  if (gSystem->AccessPathName(path)) {
     AliError(Form("File %s does not exist", path.Data()));
     return;
  }
  fFile  = TFile::Open(fCentfilename);
  owd->cd();
  fHtempV0M  = (TH1F*) (fFile->Get("hmultV0_percentile"));
  fHtempFMD  = (TH1F*) (fFile->Get("hmultFMD_percentile"));
  fHtempTRK  = (TH1F*) (fFile->Get("hNtracks_percentile"));
  fHtempTKL  = (TH1F*) (fFile->Get("hNtracklets_percentile"));
  fHtempCL0  = (TH1F*) (fFile->Get("hNclusters0_percentile"));
  fHtempCL1  = (TH1F*) (fFile->Get("hNclusters1_percentile"));

  if (!fHtempV0M) AliWarning(Form("Calibration for V0M does not exist in %s", path.Data()));
  if (!fHtempFMD) AliWarning(Form("Calibration for FMD does not exist in %s", path.Data()));
  if (!fHtempTRK) AliWarning(Form("Calibration for TRK does not exist in %s", path.Data()));
  if (!fHtempTKL) AliWarning(Form("Calibration for TKL does not exist in %s", path.Data()));
  if (!fHtempCL0) AliWarning(Form("Calibration for CL0 does not exist in %s", path.Data()));
  if (!fHtempCL1) AliWarning(Form("Calibration for CL1 does not exist in %s", path.Data()));
  
  owd->cd();
}  

//________________________________________________________________________
void AliCentralitySelectionTask::ReadCentralityHistos2(TString fCentfilename2) 
{
  //  Read centrality histograms
  TDirectory *owd = gDirectory;
  TString path = gSystem->ExpandPathName(fCentfilename2.Data());
  if (gSystem->AccessPathName(path)) {
     AliError(Form("File %s does not exist", path.Data()));
     return;
  }   
  fFile2  = TFile::Open(fCentfilename2);
  owd->cd();
  fHtempV0MvsFMD =  (TH1F*) (fFile2->Get("hmultV0vsmultFMD_all_percentile"));
  fHtempTKLvsV0M  = (TH1F*) (fFile2->Get("hNtrackletsvsmultV0_all_percentile"));
  fHtempZEMvsZDC  = (TH2F*) (fFile2->Get("hEzemvsEzdc_all_percentile"));

  if (!fHtempV0MvsFMD) AliWarning(Form("Calibration for V0MvsFMD does not exist in %s", path.Data()));
  if (!fHtempTKLvsV0M) AliWarning(Form("Calibration for TKLvsV0M does not exist in %s", path.Data()));
  if (!fHtempZEMvsZDC) AliWarning(Form("Calibration for ZEMvsZDC does not exist in %s", path.Data()));
  
  owd->cd();
}

//________________________________________________________________________
void AliCentralitySelectionTask::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  if (fFile && fFile->IsOpen())
    fFile->Close();  
  if (fFile2 && fFile2->IsOpen())
    fFile2->Close();  
}
//________________________________________________________________________
Int_t AliCentralitySelectionTask::SetupRun(AliESDEvent* const esd)
{
  // Setup files for run

  if (!esd)
    return -1;

  // check if something to be done
  if (fCurrentRun == esd->GetRunNumber())
    return 0;
  else
    fCurrentRun = esd->GetRunNumber();
  
  AliInfo(Form("Setup Centrality Selection for run %d\n",fCurrentRun));

  // CHANGE HERE FOR RUN RANGES
  switch (fPass) {
    case 1:
      if ( fCurrentRun <= 137165 ) fRunNo = 137161;
      else fRunNo = 137366;
      break;
    case 2:
      if      ( fCurrentRun >= 136851  && fCurrentRun <= 137165 ) fRunNo = 137161;
      else if ( fCurrentRun >= 137230  && fCurrentRun <= 137531 ) fRunNo = 137366;
      else if ( fCurrentRun >= 137539  && fCurrentRun <= 137848 ) fRunNo = 137722;
      else if ( fCurrentRun >= 138125  && fCurrentRun <= 138154 ) fRunNo = 138150;
      else if ( fCurrentRun >= 138190  && fCurrentRun <= 138275 ) fRunNo = 138200;
      else fRunNo = 139172;
      break;       
    default:
      AliError(Form("Run %d not known to centrality selection!", fCurrentRun));
      return -1;
  }    
  TString fileName =(Form("%s/COMMON/CENTRALITY/data/pass%d/AliCentralityBy1D_%d.root", AliAnalysisManager::GetOADBPath(), fPass, fRunNo));
  TString fileName2=(Form("%s/COMMON/CENTRALITY/data/pass%d/AliCentralityByFunction_%d.root", AliAnalysisManager::GetOADBPath(), fPass, fRunNo));
  
  AliInfo(Form("Centrality Selection for run %d and pass %d is initialized with %s", fCurrentRun, fPass, fileName.Data()));
  ReadCentralityHistos(fileName.Data());
  ReadCentralityHistos2(fileName2.Data());
  if (!fFile && !fFile2) {
     AliError(Form("Run %d not known to centrality selection!", fCurrentRun));       
     return -1;
  }   
  return 0;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MSPD(Float_t spd, Float_t v0, Int_t cent) const
{
// Clean outliers
  Float_t val= -0.579579 +  0.273949 * v0;
  Float_t spdSigma[101]={216.473, 194.377, 190.432, 185.019, 179.787, 175.99, 170.695, 167.276, 162.229, 159.016, 155.592, 151.909, 148.967, 146.367, 142.752, 139.718, 136.838, 134.133, 131.374, 128.616, 126.195, 123.49, 120.76, 118.368, 115.671, 113.354, 110.986, 108.761, 105.816, 103.373, 101.525, 99.3941, 96.8891, 95.1021, 92.728, 90.6434, 88.6393, 86.4624, 84.19, 82.3937, 80.4803, 78.5633, 76.5256, 74.7059, 73.0146, 70.9321, 69.1671, 67.2923, 65.5859, 63.9027, 61.7826, 60.1855, 58.6887, 56.8523, 55.1695, 53.422, 51.9776, 50.3506, 48.5304, 47.1517, 45.4786, 43.9914, 42.4579, 41.0889, 39.6201, 38.3046, 36.8624, 35.3697, 33.9223, 32.4637, 30.859, 29.9707, 28.764, 27.5677, 26.3437, 24.8656, 23.5253, 22.9878, 21.5172, 20.4645, 19.7614, 19.0363, 18.0875, 17.3675, 16.7313, 16.4101, 15.8235, 15.4795, 14.897, 14.578, 14.3506, 14.3506, 14.3506, 14.3506, 14.3506, 14.3506, 14.3506, 14.3506, 14.3506, 22.5965, 22.5965};

  if ( TMath::Abs(spd-val) > fOutliersCut*spdSigma[cent] ) 
    return kTRUE;
  else 
    return kFALSE;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MTPC(Int_t tracks, Float_t v0, Int_t cent) const
{
// Clean outliers
  Float_t val = -1.03873 +  0.125691 * v0;
  Float_t tpcSigma[101]={105.1, 95.7594, 94.5777, 91.7369, 89.655, 87.9469, 85.2747, 83.8137, 81.5663, 79.9871, 78.3491, 76.664, 75.1352, 73.9838, 72.0122, 70.7627, 69.1832, 67.8707, 66.3974, 65.0222, 64.1065, 62.7248, 61.498, 60.11, 58.7405, 57.681, 56.2979, 55.3153, 53.983, 52.6813, 51.5321, 50.5531, 49.3309, 48.2269, 47.1665, 46.1924, 45.1143, 44.0394, 42.8571, 41.8471, 40.9594, 39.8298, 38.8128, 37.8845, 36.8788, 35.9896, 34.9477, 34.1293, 33.1267, 32.161, 31.1132, 30.3521, 29.5237, 28.5992, 27.6283, 26.8885, 26.0719, 25.1334, 24.2125, 23.5376, 22.5383, 21.8396, 21.0384, 20.2942, 19.4746, 18.5545, 18.26, 17.2997, 16.7109, 15.9324, 15.2755, 14.764, 14.2358, 13.6728, 13.2606, 12.846, 12.2636, 11.722, 11.1491, 10.5952, 10.235, 9.87852, 9.48122, 9.17571, 8.90995, 8.81675, 8.62595, 8.59371, 8.52936, 8.58058, 8.5543, 8.5543, 8.5543, 8.5543, 8.5543, 8.5543, 8.5543, 8.5543, 8.5543, 14.2584, 14.2584};

  if ( TMath::Abs(tracks-val) > fOutliersCut*tpcSigma[cent] ) 
    return kTRUE;
  else 
    return kFALSE;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MZDC(Float_t zdc, Float_t v0) const
{
// Clean outliers
  Float_t val = 235000. - 9.5 * v0;
  if (zdc >  val) 
    return kTRUE;
  else 
    return kFALSE;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MZDCECal(Float_t /*zdc*/, Float_t /*v0*/) const
{
// Clean outliers
    return kFALSE;
}

//________________________________________________________________________  
Float_t AliCentralitySelectionTask::MyGetScaleFactor(Int_t runnumber, Int_t flag) const
{
// Get scaling factor
  if (! (runnumber >= fLowRunN && runnumber <=fHighRunN)) {
    cout << "MyGetScaleFactor error in run number range " << runnumber << endl;
    return 0.0;
  }

  Float_t scalefactor=0.0;
  if (flag==0)
    scalefactor = fV0MScaleFactor[runnumber - fLowRunN]; // subtracting reference offset index
  else if (flag==1)
    scalefactor = fSPDScaleFactor[runnumber - fLowRunN]; // subtracting reference offset index
  else if (flag==2)
    scalefactor = fTPCScaleFactor[runnumber - fLowRunN]; // subtracting reference offset index

  return scalefactor;

}

//________________________________________________________________________  
Float_t AliCentralitySelectionTask::MyGetScaleFactorMC(Int_t runnumber) const
{
// Get MC scaling factor
  if (! (runnumber >= fLowRunN && runnumber <=fHighRunN)) {
    cout << "MyGetScaleFactor error in run number range " << runnumber << endl;
    return 0.0;
  }

  Float_t scalefactor= fV0MScaleFactorMC[runnumber - fLowRunN]; // subtracting reference offset index
  return scalefactor;

}

//________________________________________________________________________  
void AliCentralitySelectionTask::MyInitScaleFactor () 
{
// Initialize the scaling factors
  for (int i=0; i<=(fHighRunN-fLowRunN); i++) fV0MScaleFactor[i] = 0.0;
  for (int i=0; i<=(fHighRunN-fLowRunN); i++) fSPDScaleFactor[i] = 0.0;
  for (int i=0; i<=(fHighRunN-fLowRunN); i++) fTPCScaleFactor[i] = 0.0;
  
  // scale factors determined from <V0 charge> on a run-by-run basis
   fV0MScaleFactor[310] = 1.;
   fV0MScaleFactor[311] = 1.;
   fV0MScaleFactor[514] = 1.0046;
   fV0MScaleFactor[515] = 0.983535;
   fV0MScaleFactor[579] = 0.988185;
   fV0MScaleFactor[580] = 0.983351;
   fV0MScaleFactor[581] = 0.989013;
   fV0MScaleFactor[583] = 0.990056;
   fV0MScaleFactor[588] = 0.974438;
   fV0MScaleFactor[589] = 0.981572;
   fV0MScaleFactor[590] = 0.989316;
   fV0MScaleFactor[592] = 0.98345;
   fV0MScaleFactor[688] = 0.993647;
   fV0MScaleFactor[690] = 0.994758;
   fV0MScaleFactor[693] = 0.989569;
   fV0MScaleFactor[698] = 0.993119;
   fV0MScaleFactor[744] = 0.989583;
   fV0MScaleFactor[757] = 0.990377;
   fV0MScaleFactor[787] = 0.990176;
   fV0MScaleFactor[788] = 0.98723;
   fV0MScaleFactor[834] = 1.00403;
   fV0MScaleFactor[835] = 0.994376;
   fV0MScaleFactor[840] = 0.99759;
   fV0MScaleFactor[842] = 1.01031;
   fV0MScaleFactor[900] = 0.996216;
   fV0MScaleFactor[901] = 0.994205;
   fV0MScaleFactor[992] = 0.998479;
   fV0MScaleFactor[997] = 1.00078;
   fV0MScaleFactor[1299] = 1.00515;
   fV0MScaleFactor[1303] = 1.00094;
   fV0MScaleFactor[1339] = 0.986596;
   fV0MScaleFactor[1346] = 0.972226;
   fV0MScaleFactor[1349] = 0.960358;
   fV0MScaleFactor[1350] = 0.970023;
   fV0MScaleFactor[1374] = 1.00575;
   fV0MScaleFactor[1545] = 1.00471;
   fV0MScaleFactor[1587] = 1.00611;
   fV0MScaleFactor[1588] = 1.00976;
   fV0MScaleFactor[1618] = 1.00771;
   fV0MScaleFactor[1682] = 1.01622;
   fV0MScaleFactor[1727] = 1.01305;
   fV0MScaleFactor[1728] = 1.00388;
   fV0MScaleFactor[1731] = 1.00673;
   fV0MScaleFactor[1732] = 1.00916;
   fV0MScaleFactor[1770] = 1.0092;
   fV0MScaleFactor[1773] = 1.00728;
   fV0MScaleFactor[1786] = 1.01655;
   fV0MScaleFactor[1787] = 1.00672;
   fV0MScaleFactor[1801] = 0.983339;
   fV0MScaleFactor[1802] = 1.00754;
   fV0MScaleFactor[1811] = 1.00608;
   fV0MScaleFactor[1815] = 1.01227;
   fV0MScaleFactor[1879] = 0.99907;
   fV0MScaleFactor[1881] = 0.995696;
   fV0MScaleFactor[2186] = 1.00559;
   fV0MScaleFactor[2187] = 1.00631;
   fV0MScaleFactor[2254] = 1.01512;
   fV0MScaleFactor[2256] = 0.998727;
   fV0MScaleFactor[2321] = 1.00701;
   fSPDScaleFactor[310] = 1.00211;
   fSPDScaleFactor[311] = 1.00067;
   fSPDScaleFactor[514] = 1.02268;
   fSPDScaleFactor[515] = 0.994902;
   fSPDScaleFactor[579] = 1.00215;
   fSPDScaleFactor[580] = 0.993421;
   fSPDScaleFactor[581] = 1.00129;
   fSPDScaleFactor[583] = 1.00242;
   fSPDScaleFactor[588] = 0.984762;
   fSPDScaleFactor[589] = 0.994355;
   fSPDScaleFactor[590] = 1.00073;
   fSPDScaleFactor[592] = 0.995889;
   fSPDScaleFactor[688] = 0.994532;
   fSPDScaleFactor[690] = 0.998307;
   fSPDScaleFactor[693] = 0.994052;
   fSPDScaleFactor[698] = 0.993224;
   fSPDScaleFactor[744] = 0.993279;
   fSPDScaleFactor[757] = 0.992494;
   fSPDScaleFactor[787] = 0.992678;
   fSPDScaleFactor[788] = 0.996563;
   fSPDScaleFactor[834] = 1.01116;
   fSPDScaleFactor[835] = 0.993108;
   fSPDScaleFactor[840] = 0.997574;
   fSPDScaleFactor[842] = 1.01829;
   fSPDScaleFactor[900] = 0.999438;
   fSPDScaleFactor[901] = 0.995849;
   fSPDScaleFactor[992] = 0.999227;
   fSPDScaleFactor[997] = 1.00575;
   fSPDScaleFactor[1299] = 0.99877;
   fSPDScaleFactor[1303] = 0.999682;
   fSPDScaleFactor[1339] = 0.978198;
   fSPDScaleFactor[1346] = 0.964178;
   fSPDScaleFactor[1349] = 0.959439;
   fSPDScaleFactor[1350] = 0.956945;
   fSPDScaleFactor[1374] = 0.994434;
   fSPDScaleFactor[1545] = 1.0016;
   fSPDScaleFactor[1587] = 1.00153;
   fSPDScaleFactor[1588] = 1.00698;
   fSPDScaleFactor[1618] = 1.00554;
   fSPDScaleFactor[1682] = 1.0123;
   fSPDScaleFactor[1727] = 1.011;
   fSPDScaleFactor[1728] = 1.00143;
   fSPDScaleFactor[1731] = 1.00486;
   fSPDScaleFactor[1732] = 1.00646;
   fSPDScaleFactor[1770] = 1.00515;
   fSPDScaleFactor[1773] = 1.00485;
   fSPDScaleFactor[1786] = 1.01215;
   fSPDScaleFactor[1787] = 1.00419;
   fSPDScaleFactor[1801] = 0.983327;
   fSPDScaleFactor[1802] = 1.00529;
   fSPDScaleFactor[1811] = 1.00367;
   fSPDScaleFactor[1815] = 1.01045;
   fSPDScaleFactor[1879] = 0.996374;
   fSPDScaleFactor[1881] = 0.988827;
   fSPDScaleFactor[2186] = 1.00354;
   fSPDScaleFactor[2187] = 1.00397;
   fSPDScaleFactor[2254] = 1.01138;
   fSPDScaleFactor[2256] = 0.996641;
   fSPDScaleFactor[2321] = 1.00357;
   fTPCScaleFactor[310] = 1.00434;
   fTPCScaleFactor[311] = 1.0056;
   fTPCScaleFactor[514] = 1.02185;
   fTPCScaleFactor[515] = 1.0036;
   fTPCScaleFactor[579] = 1.00607;
   fTPCScaleFactor[580] = 1.00183;
   fTPCScaleFactor[581] = 1.00693;
   fTPCScaleFactor[583] = 1.00746;
   fTPCScaleFactor[588] = 0.990524;
   fTPCScaleFactor[589] = 0.998582;
   fTPCScaleFactor[590] = 1.00618;
   fTPCScaleFactor[592] = 1.00088;
   fTPCScaleFactor[688] = 1.00598;
   fTPCScaleFactor[690] = 1.00658;
   fTPCScaleFactor[693] = 1.00144;
   fTPCScaleFactor[698] = 1.00279;
   fTPCScaleFactor[744] = 1.00122;
   fTPCScaleFactor[757] = 1.002;
   fTPCScaleFactor[787] = 0.997818;
   fTPCScaleFactor[788] = 0.994583;
   fTPCScaleFactor[834] = 1.01508;
   fTPCScaleFactor[835] = 1.00218;
   fTPCScaleFactor[840] = 1.00569;
   fTPCScaleFactor[842] = 1.01789;
   fTPCScaleFactor[900] = 1.00739;
   fTPCScaleFactor[901] = 1.00462;
   fTPCScaleFactor[992] = 1.00502;
   fTPCScaleFactor[997] = 1.00943;
   fTPCScaleFactor[1299] = 1.00438;
   fTPCScaleFactor[1303] = 0.996701;
   fTPCScaleFactor[1339] = 0.978641;
   fTPCScaleFactor[1346] = 0.968906;
   fTPCScaleFactor[1349] = 0.954311;
   fTPCScaleFactor[1350] = 0.958764;
   fTPCScaleFactor[1374] = 0.997899;
   fTPCScaleFactor[1545] = 0.992;
   fTPCScaleFactor[1587] = 0.992635;
   fTPCScaleFactor[1588] = 1.00207;
   fTPCScaleFactor[1618] = 1.00126;
   fTPCScaleFactor[1682] = 1.00324;
   fTPCScaleFactor[1727] = 1.00042;
   fTPCScaleFactor[1728] = 0.978881;
   fTPCScaleFactor[1731] = 0.999818;
   fTPCScaleFactor[1732] = 1.00109;
   fTPCScaleFactor[1770] = 0.99935;
   fTPCScaleFactor[1773] = 0.998531;
   fTPCScaleFactor[1786] = 0.999125;
   fTPCScaleFactor[1787] = 0.998479;
   fTPCScaleFactor[1801] = 0.9775;
   fTPCScaleFactor[1802] = 0.999095;
   fTPCScaleFactor[1811] = 0.998197;
   fTPCScaleFactor[1815] = 1.00413;
   fTPCScaleFactor[1879] = 0.990916;
   fTPCScaleFactor[1881] = 0.987241;
   fTPCScaleFactor[2186] = 1.00048;
   fTPCScaleFactor[2187] = 1.00057;
   fTPCScaleFactor[2254] = 1.00588;
   fTPCScaleFactor[2256] = 0.991503;
   fTPCScaleFactor[2321] = 1.00057;


  // set all missing values to the value of the run before it ....
  for (int i=0; i<=(fHighRunN-fLowRunN); i++) {    
    if (fV0MScaleFactor[i] == 0.0) {     
      if (i==0) {
	fV0MScaleFactor[i] = 1.00;
      } else {
	// search for last run number with non-zero value
	for (int j=i-1;j>=0;j--) {
	  if (fV0MScaleFactor[j] != 0.0) {
	    fV0MScaleFactor[i] = fV0MScaleFactor[j];
	    break;
	  }
	}
      }
    }
  } // end loop over checking all run-numbers 

  for (int i=0; i<=(fHighRunN-fLowRunN); i++) {    
    if (fSPDScaleFactor[i] == 0.0) {     
      if (i==0) {
	fSPDScaleFactor[i] = 1.00;
      } else {
	for (int j=i-1;j>=0;j--) {
	  if (fSPDScaleFactor[j] != 0.0) {
	    fSPDScaleFactor[i] = fSPDScaleFactor[j];
	    break;
	  }
	}
      }
    }
  } 

  for (int i=0; i<=(fHighRunN-fLowRunN); i++) {    
    if (fTPCScaleFactor[i] == 0.0) {     
      if (i==0) {
	fTPCScaleFactor[i] = 1.00;
      } else {
	for (int j=i-1;j>=0;j--) {
	  if (fTPCScaleFactor[j] != 0.0) {
	    fTPCScaleFactor[i] = fTPCScaleFactor[j];
	    break;
	  }
	}
      }
    }
  } 
  

  //    for (int i=0; i<=(fHighRunN-fLowRunN); i++) cout << "Scale Factor = " << fV0MScaleFactor[i] << " for Run " << i+fLowRunN << endl;
  
  return;

}


//________________________________________________________________________  
void AliCentralitySelectionTask::MyInitScaleFactorMC() 
{
// Initialize the MC scaling factors
  for (int i=0; i<(fHighRunN-fLowRunN); i++) fV0MScaleFactorMC[i] = 0.0;
  // scale factors determined from <V0 charge> on a run-by-run basis
  switch (fPass) {
  case 1:
    fV0MScaleFactorMC[0] = 0.75108;
    break;
  case 2:
    fV0MScaleFactorMC[0] = 0.8;
    break;
  default:
    AliError(Form("Unknown reconstruction pass (%d)! Setting MC scale in V0 to 1", fPass));
    fV0MScaleFactorMC[0] = 1.0;
    break;
  }

  // set all missing values to the value of the run before it ....
  for (int i=0; i<=(fHighRunN-fLowRunN); i++) {    
    if (fV0MScaleFactorMC[i] == 0.0) {     
      if (i==0) {
	fV0MScaleFactorMC[i] = 1.00;
      } else {
	// search for last run number with non-zero value
	for (int j=i-1;j>=0;j--) {
	  if (fV0MScaleFactorMC[j] != 0.0) {
	    fV0MScaleFactorMC[i] = fV0MScaleFactorMC[j];
	    break;
	  }
	}
      }
    }
  } // end loop over checking all run-numbers 


  //    for (int i=0; i<=(fHighRunN-fLowRunN); i++) cout << "Scale Factor = " << fV0MScaleFactorMC[i] << " for Run " << i+fLowRunN << endl;
  
  return;

}

