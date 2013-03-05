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
#include "AliESDTZERO.h"
#include "AliESDVZERO.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliOADBCentrality.h"
#include "AliOADBContainer.h"
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
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
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
  fCurrentRun(-1),
  fUseScaling(0),
  fUseCleaning(0),
  fFillHistos(0),
  fV0MScaleFactor(0),
  fSPDScaleFactor(0),
  fTPCScaleFactor(0),
  fV0MScaleFactorMC(0),
  fV0MSPDOutlierPar0(0),  
  fV0MSPDOutlierPar1(0),  
  fV0MTPCOutlierPar0(0),  
  fV0MTPCOutlierPar1(0),  
  fV0MSPDSigmaOutlierPar0(0),  
  fV0MSPDSigmaOutlierPar1(0),  
  fV0MSPDSigmaOutlierPar2(0),  
  fV0MTPCSigmaOutlierPar0(0),  
  fV0MTPCSigmaOutlierPar1(0),  
  fV0MTPCSigmaOutlierPar2(0),   
  fV0MZDCOutlierPar0(0),	    
  fV0MZDCOutlierPar1(0),	    
  fV0MZDCEcalOutlierPar0(0),   
  fV0MZDCEcalOutlierPar1(0),   
  fTrackCuts(0),
  fEsdTrackCuts(0),
  fEsdTrackCutsExtra1(0),
  fEsdTrackCutsExtra2(0),
  fZVCut(10),
  fOutliersCut(5),
  fQuality(999),
  fIsSelected(0),
  fMSL(0),
  fMSH(0),
  fMUL(0),
  fMLL(0),
  fEJE(0),
  fEGA(0),
  fPHS(0),
  fMB(0),
  fCVHN(0),
  fCVLN(0),
  fCVHNbit(0),
  fCVLNbit(0),
  fCCENT(0),
  fCSEMI(0),
  fCCENTbit(0),
  fCSEMIbit(0),
  fCentV0M(0),
  fCentV0A(0),
  fCentV0C(0),
  fCentV0MEq(0),
  fCentV0AEq(0),
  fCentV0CEq(0),
  fCentFMD(0),
  fCentTRK(0),
  fCentTKL(0),
  fCentCL0(0),
  fCentCL1(0),
  fCentCND(0),
  fCentZNA(0),
  fCentZNC(0),
  fCentNPA(0),
  fCentV0MvsFMD(0),
  fCentTKLvsV0M(0),
  fCentZEMvsZDC(0),
  fCentV0Mtrue(0),
  fCentV0Atrue(0),
  fCentV0Ctrue(0),
  fCentV0MEqtrue(0),
  fCentV0AEqtrue(0),
  fCentV0CEqtrue(0),
  fCentFMDtrue(0),
  fCentTRKtrue(0),
  fCentTKLtrue(0),
  fCentCL0true(0),
  fCentCL1true(0),
  fCentCNDtrue(0),
  fCentZNAtrue(0),
  fCentZNCtrue(0),
  fHtempV0M(0),
  fHtempV0A(0),
  fHtempV0C(0),
  fHtempV0MEq(0),
  fHtempV0AEq(0),
  fHtempV0CEq(0),
  fHtempFMD(0),
  fHtempTRK(0),
  fHtempTKL(0),
  fHtempCL0(0),
  fHtempCL1(0),
  fHtempCND(0),
  fHtempZNA(0),
  fHtempZNC(0),
  fHtempV0MvsFMD(0),
  fHtempTKLvsV0M(0),
  fHtempZEMvsZDC(0),
  fHtempNPA(0),
  fHtempV0Mtrue(0),
  fHtempV0Atrue(0),
  fHtempV0Ctrue(0),
  fHtempV0MEqtrue(0),
  fHtempV0AEqtrue(0),
  fHtempV0CEqtrue(0),
  fHtempFMDtrue(0),
  fHtempTRKtrue(0),
  fHtempTKLtrue(0),
  fHtempCL0true(0),
  fHtempCL1true(0),
  fHtempCNDtrue(0),
  fHtempZNAtrue(0),
  fHtempZNCtrue(0),
  fOutputList(0),
  fHOutCentV0M(0),
  fHOutCentV0A(0),
  fHOutCentV0C(0),
  fHOutCentV0MEq(0),
  fHOutCentV0AEq(0),
  fHOutCentV0CEq(0),
  fHOutCentV0MCVHN(0),
  fHOutCentV0MCVLN(0),
  fHOutCentV0MCVHNinMB(0),
  fHOutCentV0MCVLNinMB(0),
  fHOutCentV0MCCENT(0),
  fHOutCentV0MCSEMI(0),
  fHOutCentV0MCCENTinMB(0),
  fHOutCentV0MCSEMIinMB(0),
  fHOutCentV0MMSL(0),
  fHOutCentV0MMSH(0),
  fHOutCentV0MMUL(0),
  fHOutCentV0MMLL(0),
  fHOutCentV0MEJE(0),
  fHOutCentV0MEGA(0),
  fHOutCentV0MPHS(0),
  fHOutCentV0MMSLinMB(0),
  fHOutCentV0MMSHinMB(0),
  fHOutCentV0MMULinMB(0),
  fHOutCentV0MMLLinMB(0),
  fHOutCentV0MEJEinMB(0),
  fHOutCentV0MEGAinMB(0),
  fHOutCentV0MPHSinMB(0),
  fHOutCentFMD(0),
  fHOutCentTRK(0),
  fHOutCentTKL(0),
  fHOutCentCL0(0),
  fHOutCentCL1(0),
  fHOutCentCND(0),
  fHOutCentNPA(0),
  fHOutCentZNA(0),
  fHOutCentZNC(0),
  fHOutCentV0MvsFMD(0),
  fHOutCentTKLvsV0M(0),
  fHOutCentZEMvsZDC(0),
  fHOutCentV0MvsCentCL1(0),
  fHOutCentV0MvsCentTRK(0),
  fHOutCentTRKvsCentCL1(0),
  fHOutCentV0MvsCentZDC(0),
  fHOutCentV0AvsCentV0C(0),
  fHOutCentV0AvsCentTRK(0),
  fHOutCentV0AvsCentCND(0),
  fHOutCentV0AvsCentCL1(0),
  fHOutCentV0CvsCentTRK(0),
  fHOutCentV0CvsCentCND(0),
  fHOutCentV0CvsCentCL1(0),
  fHOutCentNPAvsCentV0A(0),
  fHOutCentNPAvsCentV0C(0),
  fHOutCentNPAvsCentTRK(0),
  fHOutCentNPAvsCentCND(0),
  fHOutCentNPAvsCentCL1(0),
  fHOutCentZNAvsCentV0A(0),
  fHOutCentZNAvsCentV0C(0),
  fHOutCentZNAvsCentTRK(0),
  fHOutCentZNAvsCentCND(0),
  fHOutCentZNAvsCentCL1(0),
  fHOutMultV0AC(0),
  fHOutMultV0M(0),
  fHOutMultV0A(0),
  fHOutMultV0C(0),
  fHOutMultV0MEq(0),
  fHOutMultV0AEq(0),
  fHOutMultV0CEq(0),
  fHOutMultV0Mnc(0),
  fHOutMultV0Anc(0),
  fHOutMultV0Cnc(0),
  fHOutMultV0O(0),
  fHOutMultV0Cells(0),
  fHOutMultFMD(0),
  fHOutMultTRK(0),
  fHOutMultTKL(0),
  fHOutMultCL0(0),
  fHOutMultCL1(0),
  fHOutMultCND(0),
  fHOutMultNPA(0),
  fHOutMultZNA(0),
  fHOutMultZNC(0),
  fHOutMultV0MvsZDN(0),
  fHOutMultZEMvsZDN(0),
  fHOutMultV0MvsZDC(0),
  fHOutMultZEMvsZDC(0),
  fHOutMultZEMvsZDCw(0),
  fHOutMultV0MvsCL1(0),
  fHOutMultV0MvsTRK(0),
  fHOutMultTRKvsCL1(0),
  fHOutMultV0MvsV0O(0),
  fHOutMultV0OvsCL1(0),
  fHOutMultV0OvsTRK(0),
  fHOutMultCL1vsTKL(0),
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
  fHOutVertex(0),
  fHOutVertexT0(0)
{   
  // Default constructor
  AliInfo("Centrality Selection enabled.");

  fUseScaling=kTRUE;
  fUseCleaning=kTRUE;
  fFillHistos=kFALSE;
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,AliESDZDC.,AliESDFMD.,AliESDVZERO.,AliESDTZERO."
    ",SPDVertex.,TPCVertex.,PrimaryVertex.,AliMultiplicity.,Tracks ";
}   

//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fCurrentRun(-1),
  fUseScaling(0),
  fUseCleaning(0),
  fFillHistos(0),
  fV0MScaleFactor(0),
  fSPDScaleFactor(0),
  fTPCScaleFactor(0),
  fV0MScaleFactorMC(0),
  fV0MSPDOutlierPar0(0),  
  fV0MSPDOutlierPar1(0),  
  fV0MTPCOutlierPar0(0),  
  fV0MTPCOutlierPar1(0),  
  fV0MSPDSigmaOutlierPar0(0),  
  fV0MSPDSigmaOutlierPar1(0),  
  fV0MSPDSigmaOutlierPar2(0),  
  fV0MTPCSigmaOutlierPar0(0),  
  fV0MTPCSigmaOutlierPar1(0),  
  fV0MTPCSigmaOutlierPar2(0),   
  fV0MZDCOutlierPar0(0),	    
  fV0MZDCOutlierPar1(0),	    
  fV0MZDCEcalOutlierPar0(0),   
  fV0MZDCEcalOutlierPar1(0),   
  fTrackCuts(0),
  fEsdTrackCuts(0),
  fEsdTrackCutsExtra1(0),
  fEsdTrackCutsExtra2(0),
  fZVCut(10),
  fOutliersCut(5),
  fQuality(999),
  fIsSelected(0),
  fMSL(0),
  fMSH(0),
  fMUL(0),
  fMLL(0),
  fEJE(0),
  fEGA(0),
  fPHS(0),
  fMB(0),
  fCVHN(0),
  fCVLN(0),
  fCVHNbit(0),
  fCVLNbit(0),
  fCCENT(0),
  fCSEMI(0),
  fCCENTbit(0),
  fCSEMIbit(0),
  fCentV0M(0),
  fCentV0A(0),
  fCentV0C(0),
  fCentV0MEq(0),
  fCentV0AEq(0),
  fCentV0CEq(0),
  fCentFMD(0),
  fCentTRK(0),
  fCentTKL(0),
  fCentCL0(0),
  fCentCL1(0),
  fCentCND(0),
  fCentZNA(0),
  fCentZNC(0),
  fCentNPA(0),
  fCentV0MvsFMD(0),
  fCentTKLvsV0M(0),
  fCentZEMvsZDC(0),
  fCentV0Mtrue(0),
  fCentV0Atrue(0),
  fCentV0Ctrue(0),
  fCentV0MEqtrue(0),
  fCentV0AEqtrue(0),
  fCentV0CEqtrue(0),
  fCentFMDtrue(0),
  fCentTRKtrue(0),
  fCentTKLtrue(0),
  fCentCL0true(0),
  fCentCL1true(0),
  fCentCNDtrue(0),
  fCentZNAtrue(0),
  fCentZNCtrue(0),
  fHtempV0M(0),
  fHtempV0A(0),
  fHtempV0C(0),
  fHtempV0MEq(0),
  fHtempV0AEq(0),
  fHtempV0CEq(0),
  fHtempFMD(0),
  fHtempTRK(0),
  fHtempTKL(0),
  fHtempCL0(0),
  fHtempCL1(0),
  fHtempCND(0),
  fHtempZNA(0),
  fHtempZNC(0),
  fHtempV0MvsFMD(0),
  fHtempTKLvsV0M(0),
  fHtempZEMvsZDC(0),
  fHtempNPA(0),
  fHtempV0Mtrue(0),
  fHtempV0Atrue(0),
  fHtempV0Ctrue(0),
  fHtempV0MEqtrue(0),
  fHtempV0AEqtrue(0),
  fHtempV0CEqtrue(0),
  fHtempFMDtrue(0),
  fHtempTRKtrue(0),
  fHtempTKLtrue(0),
  fHtempCL0true(0),
  fHtempCL1true(0),
  fHtempCNDtrue(0),
  fHtempZNAtrue(0),
  fHtempZNCtrue(0),
  fOutputList(0),
  fHOutCentV0M(0),
  fHOutCentV0A(0),
  fHOutCentV0C(0),
  fHOutCentV0MEq(0),
  fHOutCentV0AEq(0),
  fHOutCentV0CEq(0),
  fHOutCentV0MCVHN(0),
  fHOutCentV0MCVLN(0),
  fHOutCentV0MCVHNinMB(0),
  fHOutCentV0MCVLNinMB(0),
  fHOutCentV0MCCENT(0),
  fHOutCentV0MCSEMI(0),
  fHOutCentV0MCCENTinMB(0),
  fHOutCentV0MCSEMIinMB(0),
  fHOutCentV0MMSL(0),
  fHOutCentV0MMSH(0),
  fHOutCentV0MMUL(0),
  fHOutCentV0MMLL(0),
  fHOutCentV0MEJE(0),
  fHOutCentV0MEGA(0),
  fHOutCentV0MPHS(0),
  fHOutCentV0MMSLinMB(0),
  fHOutCentV0MMSHinMB(0),
  fHOutCentV0MMULinMB(0),
  fHOutCentV0MMLLinMB(0),
  fHOutCentV0MEJEinMB(0),
  fHOutCentV0MEGAinMB(0),
  fHOutCentV0MPHSinMB(0),
  fHOutCentFMD(0),
  fHOutCentTRK(0),
  fHOutCentTKL(0),
  fHOutCentCL0(0),
  fHOutCentCL1(0),
  fHOutCentCND(0),
  fHOutCentNPA(0),
  fHOutCentZNA(0),
  fHOutCentZNC(0),
  fHOutCentV0MvsFMD(0),
  fHOutCentTKLvsV0M(0),
  fHOutCentZEMvsZDC(0),
  fHOutCentV0MvsCentCL1(0),
  fHOutCentV0MvsCentTRK(0),
  fHOutCentTRKvsCentCL1(0),
  fHOutCentV0MvsCentZDC(0),
  fHOutCentV0AvsCentV0C(0),
  fHOutCentV0AvsCentTRK(0),
  fHOutCentV0AvsCentCND(0),
  fHOutCentV0AvsCentCL1(0),
  fHOutCentV0CvsCentTRK(0),
  fHOutCentV0CvsCentCND(0),
  fHOutCentV0CvsCentCL1(0),
  fHOutCentNPAvsCentV0A(0),
  fHOutCentNPAvsCentV0C(0),
  fHOutCentNPAvsCentTRK(0),
  fHOutCentNPAvsCentCND(0),
  fHOutCentNPAvsCentCL1(0),
  fHOutCentZNAvsCentV0A(0),
  fHOutCentZNAvsCentV0C(0),
  fHOutCentZNAvsCentTRK(0),
  fHOutCentZNAvsCentCND(0),
  fHOutCentZNAvsCentCL1(0),
  fHOutMultV0AC(0),
  fHOutMultV0M(0),
  fHOutMultV0A(0),
  fHOutMultV0C(0),
  fHOutMultV0MEq(0),
  fHOutMultV0AEq(0),
  fHOutMultV0CEq(0),
  fHOutMultV0Mnc(0),
  fHOutMultV0Anc(0),
  fHOutMultV0Cnc(0),
  fHOutMultV0O(0),
  fHOutMultV0Cells(0),
  fHOutMultFMD(0),
  fHOutMultTRK(0),
  fHOutMultTKL(0),
  fHOutMultCL0(0),
  fHOutMultCL1(0),
  fHOutMultCND(0),
  fHOutMultNPA(0),
  fHOutMultZNA(0),
  fHOutMultZNC(0),
  fHOutMultV0MvsZDN(0),
  fHOutMultZEMvsZDN(0),
  fHOutMultV0MvsZDC(0),
  fHOutMultZEMvsZDC(0),
  fHOutMultZEMvsZDCw(0),
  fHOutMultV0MvsCL1(0),
  fHOutMultV0MvsTRK(0),
  fHOutMultTRKvsCL1(0),
  fHOutMultV0MvsV0O(0),
  fHOutMultV0OvsCL1(0),
  fHOutMultV0OvsTRK(0),
  fHOutMultCL1vsTKL(0),
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
  fHOutVertex(0),
  fHOutVertexT0(0)
{
  // Default constructor
  AliInfo("Centrality Selection enabled.");
  //DefineOutput(1, TList::Class());
  fUseScaling=kTRUE;
  fUseCleaning=kTRUE;
  fFillHistos=kFALSE;
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,AliESDZDC.,AliESDFMD.,AliESDVZERO.,AliESDTZERO."
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
  fAnalysisInput(ana.fAnalysisInput),
  fIsMCInput(ana.fIsMCInput),
  fCurrentRun(ana.fCurrentRun),
  fUseScaling(ana.fUseScaling),
  fUseCleaning(ana.fUseCleaning),
  fFillHistos(ana.fFillHistos),
  fV0MScaleFactor(ana.fV0MScaleFactor),
  fSPDScaleFactor(ana.fSPDScaleFactor),
  fTPCScaleFactor(ana.fTPCScaleFactor),
  fV0MScaleFactorMC(ana.fV0MScaleFactorMC),
  fV0MSPDOutlierPar0(ana.fV0MSPDOutlierPar0),  
  fV0MSPDOutlierPar1(ana.fV0MSPDOutlierPar1),  
  fV0MTPCOutlierPar0(ana.fV0MTPCOutlierPar0),  
  fV0MTPCOutlierPar1(ana.fV0MTPCOutlierPar1),  
  fV0MSPDSigmaOutlierPar0(ana.fV0MSPDSigmaOutlierPar0),  
  fV0MSPDSigmaOutlierPar1(ana.fV0MSPDSigmaOutlierPar1),  
  fV0MSPDSigmaOutlierPar2(ana.fV0MSPDSigmaOutlierPar2),  
  fV0MTPCSigmaOutlierPar0(ana.fV0MTPCSigmaOutlierPar0),  
  fV0MTPCSigmaOutlierPar1(ana.fV0MTPCSigmaOutlierPar1),  
  fV0MTPCSigmaOutlierPar2(ana.fV0MTPCSigmaOutlierPar2), 
  fV0MZDCOutlierPar0(ana.fV0MZDCOutlierPar0),	    
  fV0MZDCOutlierPar1(ana.fV0MZDCOutlierPar1),	    
  fV0MZDCEcalOutlierPar0(ana.fV0MZDCEcalOutlierPar0),   
  fV0MZDCEcalOutlierPar1(ana.fV0MZDCEcalOutlierPar1),   
  fTrackCuts(ana.fTrackCuts),
  fEsdTrackCuts(ana.fEsdTrackCuts),
  fEsdTrackCutsExtra1(ana.fEsdTrackCutsExtra1),
  fEsdTrackCutsExtra2(ana.fEsdTrackCutsExtra2),
  fZVCut(ana.fZVCut),
  fOutliersCut(ana.fOutliersCut),
  fQuality(ana.fQuality),
  fIsSelected(ana.fIsSelected),
  fMSL(ana.fMSL),
  fMSH(ana.fMSH),
  fMUL(ana.fMUL),
  fMLL(ana.fMLL),
  fEJE(ana.fEJE),
  fEGA(ana.fEGA),
  fPHS(ana.fPHS),
  fMB(ana.fMB),
  fCVHN(ana.fCVHN),
  fCVLN(ana.fCVLN),
  fCVHNbit(ana.fCVHNbit),
  fCVLNbit(ana.fCVLNbit),
  fCCENT(ana.fCCENT),
  fCSEMI(ana.fCSEMI),
  fCCENTbit(ana.fCCENTbit),
  fCSEMIbit(ana.fCSEMIbit),
  fCentV0M(ana.fCentV0M),
  fCentV0A(ana.fCentV0A),
  fCentV0C(ana.fCentV0C),
  fCentV0MEq(ana.fCentV0MEq),
  fCentV0AEq(ana.fCentV0AEq),
  fCentV0CEq(ana.fCentV0CEq),
  fCentFMD(ana.fCentFMD),
  fCentTRK(ana.fCentTRK),
  fCentTKL(ana.fCentTKL),
  fCentCL0(ana.fCentCL0),
  fCentCL1(ana.fCentCL1),
  fCentCND(ana.fCentCND),
  fCentZNA(ana.fCentZNA),
  fCentZNC(ana.fCentZNC),
  fCentNPA(ana.fCentNPA),
  fCentV0MvsFMD(ana.fCentV0MvsFMD),
  fCentTKLvsV0M(ana.fCentTKLvsV0M),
  fCentZEMvsZDC(ana.fCentZEMvsZDC),
  fCentV0Mtrue(ana.fCentV0Mtrue),
  fCentV0Atrue(ana.fCentV0Atrue),
  fCentV0Ctrue(ana.fCentV0Ctrue),
  fCentV0MEqtrue(ana.fCentV0MEqtrue),
  fCentV0AEqtrue(ana.fCentV0AEqtrue),
  fCentV0CEqtrue(ana.fCentV0CEqtrue),
  fCentFMDtrue(ana.fCentFMDtrue),
  fCentTRKtrue(ana.fCentTRKtrue),
  fCentTKLtrue(ana.fCentTKLtrue),
  fCentCL0true(ana.fCentCL0true),
  fCentCL1true(ana.fCentCL1true),
  fCentCNDtrue(ana.fCentCNDtrue),
  fCentZNAtrue(ana.fCentZNAtrue),
  fCentZNCtrue(ana.fCentZNCtrue),
  fHtempV0M(ana.fHtempV0M),
  fHtempV0A(ana.fHtempV0A),
  fHtempV0C(ana.fHtempV0C),
  fHtempV0MEq(ana.fHtempV0MEq),
  fHtempV0AEq(ana.fHtempV0AEq),
  fHtempV0CEq(ana.fHtempV0CEq),
  fHtempFMD(ana.fHtempFMD),
  fHtempTRK(ana.fHtempTRK),
  fHtempTKL(ana.fHtempTKL),
  fHtempCL0(ana.fHtempCL0),
  fHtempCL1(ana.fHtempCL1),
  fHtempCND(ana.fHtempCND),
  fHtempZNA(ana.fHtempZNA),
  fHtempZNC(ana.fHtempZNC),
  fHtempV0MvsFMD(ana.fHtempV0MvsFMD),
  fHtempTKLvsV0M(ana.fHtempTKLvsV0M),
  fHtempZEMvsZDC(ana.fHtempZEMvsZDC),
  fHtempNPA(ana.fHtempNPA),
  fHtempV0Mtrue(ana.fHtempV0Mtrue),
  fHtempV0Atrue(ana.fHtempV0Atrue),
  fHtempV0Ctrue(ana.fHtempV0Ctrue),
  fHtempV0MEqtrue(ana.fHtempV0MEqtrue),
  fHtempV0AEqtrue(ana.fHtempV0AEqtrue),
  fHtempV0CEqtrue(ana.fHtempV0CEqtrue),
  fHtempFMDtrue(ana.fHtempFMDtrue),
  fHtempTRKtrue(ana.fHtempTRKtrue),
  fHtempTKLtrue(ana.fHtempTKLtrue),
  fHtempCL0true(ana.fHtempCL0true),
  fHtempCL1true(ana.fHtempCL1true),
  fHtempCNDtrue(ana.fHtempCNDtrue),
  fHtempZNAtrue(ana.fHtempZNAtrue),
  fHtempZNCtrue(ana.fHtempZNCtrue),
  fOutputList(ana.fOutputList),
  fHOutCentV0M(ana.fHOutCentV0M),
  fHOutCentV0A(ana.fHOutCentV0A),
  fHOutCentV0C(ana.fHOutCentV0C),
  fHOutCentV0MEq(ana.fHOutCentV0MEq),
  fHOutCentV0AEq(ana.fHOutCentV0AEq),
  fHOutCentV0CEq(ana.fHOutCentV0CEq),
  fHOutCentV0MCVHN(ana.fHOutCentV0MCVHN),
  fHOutCentV0MCVLN(ana.fHOutCentV0MCVLN),
  fHOutCentV0MCVHNinMB(ana.fHOutCentV0MCVHNinMB),
  fHOutCentV0MCVLNinMB(ana.fHOutCentV0MCVLNinMB),
  fHOutCentV0MCCENT(ana.fHOutCentV0MCCENT),
  fHOutCentV0MCSEMI(ana.fHOutCentV0MCSEMI),
  fHOutCentV0MCCENTinMB(ana.fHOutCentV0MCCENTinMB),
  fHOutCentV0MCSEMIinMB(ana.fHOutCentV0MCSEMIinMB),
  fHOutCentV0MMSL(ana.fHOutCentV0MMSL),
  fHOutCentV0MMSH(ana.fHOutCentV0MMSH),
  fHOutCentV0MMUL(ana.fHOutCentV0MMUL),
  fHOutCentV0MMLL(ana.fHOutCentV0MMLL),
  fHOutCentV0MEJE(ana.fHOutCentV0MEJE),
  fHOutCentV0MEGA(ana.fHOutCentV0MEGA),
  fHOutCentV0MPHS(ana.fHOutCentV0MPHS),
  fHOutCentV0MMSLinMB(ana.fHOutCentV0MMSLinMB),
  fHOutCentV0MMSHinMB(ana.fHOutCentV0MMSHinMB),
  fHOutCentV0MMULinMB(ana.fHOutCentV0MMULinMB),
  fHOutCentV0MMLLinMB(ana.fHOutCentV0MMLLinMB),
  fHOutCentV0MEJEinMB(ana.fHOutCentV0MEJEinMB),
  fHOutCentV0MEGAinMB(ana.fHOutCentV0MEGAinMB),
  fHOutCentV0MPHSinMB(ana.fHOutCentV0MPHSinMB),
  fHOutCentFMD(ana.fHOutCentFMD),
  fHOutCentTRK(ana.fHOutCentTRK),
  fHOutCentTKL(ana.fHOutCentTKL),
  fHOutCentCL0(ana.fHOutCentCL0),
  fHOutCentCL1(ana.fHOutCentCL1),
  fHOutCentCND(ana.fHOutCentCND),
  fHOutCentNPA(ana.fHOutCentNPA),
  fHOutCentZNA(ana.fHOutCentZNA),
  fHOutCentZNC(ana.fHOutCentZNC),
  fHOutCentV0MvsFMD(ana.fHOutCentV0MvsFMD),
  fHOutCentTKLvsV0M(ana.fHOutCentTKLvsV0M),
  fHOutCentZEMvsZDC(ana.fHOutCentZEMvsZDC),
  fHOutCentV0MvsCentCL1(ana.fHOutCentV0MvsCentCL1),
  fHOutCentV0MvsCentTRK(ana.fHOutCentV0MvsCentTRK),
  fHOutCentTRKvsCentCL1(ana.fHOutCentTRKvsCentCL1),
  fHOutCentV0MvsCentZDC(ana.fHOutCentV0MvsCentZDC),
  fHOutCentV0AvsCentV0C(ana.fHOutCentV0AvsCentV0C),
  fHOutCentV0AvsCentTRK(ana.fHOutCentV0AvsCentTRK),
  fHOutCentV0AvsCentCND(ana.fHOutCentV0AvsCentCND),
  fHOutCentV0AvsCentCL1(ana.fHOutCentV0AvsCentCL1),
  fHOutCentV0CvsCentTRK(ana.fHOutCentV0CvsCentTRK),
  fHOutCentV0CvsCentCND(ana.fHOutCentV0CvsCentCND),
  fHOutCentV0CvsCentCL1(ana.fHOutCentV0CvsCentCL1),
  fHOutCentNPAvsCentV0A(ana.fHOutCentNPAvsCentV0A),
  fHOutCentNPAvsCentV0C(ana.fHOutCentNPAvsCentV0C),
  fHOutCentNPAvsCentTRK(ana.fHOutCentNPAvsCentTRK),
  fHOutCentNPAvsCentCND(ana.fHOutCentNPAvsCentCND),
  fHOutCentNPAvsCentCL1(ana.fHOutCentNPAvsCentCL1),
  fHOutCentZNAvsCentV0A(ana.fHOutCentZNAvsCentV0A),
  fHOutCentZNAvsCentV0C(ana.fHOutCentZNAvsCentV0C),
  fHOutCentZNAvsCentTRK(ana.fHOutCentZNAvsCentTRK),
  fHOutCentZNAvsCentCND(ana.fHOutCentZNAvsCentCND),
  fHOutCentZNAvsCentCL1(ana.fHOutCentZNAvsCentCL1),
  fHOutMultV0AC(ana.fHOutMultV0AC),
  fHOutMultV0M(ana.fHOutMultV0M),
  fHOutMultV0A(ana.fHOutMultV0A),
  fHOutMultV0C(ana.fHOutMultV0C),
  fHOutMultV0MEq(ana.fHOutMultV0MEq),
  fHOutMultV0AEq(ana.fHOutMultV0AEq),
  fHOutMultV0CEq(ana.fHOutMultV0CEq),
  fHOutMultV0Mnc(ana.fHOutMultV0Mnc),
  fHOutMultV0Anc(ana.fHOutMultV0Anc),
  fHOutMultV0Cnc(ana.fHOutMultV0Cnc),
  fHOutMultV0O(ana.fHOutMultV0O),
  fHOutMultV0Cells(ana.fHOutMultV0Cells),
  fHOutMultFMD(ana.fHOutMultFMD),
  fHOutMultTRK(ana.fHOutMultTRK),
  fHOutMultTKL(ana.fHOutMultTKL),
  fHOutMultCL0(ana.fHOutMultCL0),
  fHOutMultCL1(ana.fHOutMultCL1),
  fHOutMultCND(ana.fHOutMultCND),
  fHOutMultNPA(ana.fHOutMultNPA),
  fHOutMultZNA(ana.fHOutMultZNA),
  fHOutMultZNC(ana.fHOutMultZNC),
  fHOutMultV0MvsZDN(ana.fHOutMultV0MvsZDN),
  fHOutMultZEMvsZDN(ana.fHOutMultZEMvsZDN),
  fHOutMultV0MvsZDC(ana.fHOutMultV0MvsZDC),
  fHOutMultZEMvsZDC(ana.fHOutMultZEMvsZDC),
  fHOutMultZEMvsZDCw(ana.fHOutMultZEMvsZDCw),
  fHOutMultV0MvsCL1(ana.fHOutMultV0MvsCL1),
  fHOutMultV0MvsTRK(ana.fHOutMultV0MvsTRK),
  fHOutMultTRKvsCL1(ana.fHOutMultTRKvsCL1),
  fHOutMultV0MvsV0O(ana.fHOutMultV0MvsV0O),
  fHOutMultV0OvsCL1(ana.fHOutMultV0OvsCL1),
  fHOutMultV0OvsTRK(ana.fHOutMultV0OvsTRK),
  fHOutMultCL1vsTKL(ana.fHOutMultCL1vsTKL),
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
  fHOutVertex(ana.fHOutVertex),
  fHOutVertexT0(ana.fHOutVertexT0)
{
  // Copy Constructor	

}

//________________________________________________________________________
AliCentralitySelectionTask::~AliCentralitySelectionTask()
{
  // Destructor  
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fOutputList;
  if (fTrackCuts) delete fTrackCuts;
  if (fEsdTrackCuts) delete fEsdTrackCuts;
  if (fEsdTrackCutsExtra1) delete fEsdTrackCutsExtra1;
  if (fEsdTrackCutsExtra2) delete fEsdTrackCutsExtra2;
}  

//________________________________________________________________________
void AliCentralitySelectionTask::UserCreateOutputObjects()
{  
  // Create the output containers
  if(fDebug>1) printf("AnalysisCentralitySelectionTask::UserCreateOutputObjects() \n");
  AliLog::SetClassDebugLevel("AliCentralitySelectionTask", AliLog::kInfo);

  if (fFillHistos) {    
    fOutputList = new TList();
    fOutputList->SetOwner();
    fHOutCentV0M     = new TH1F("fHOutCentV0M","fHOutCentV0M; Centrality V0",505,0,101);
    fHOutCentV0A    = new TH1F("fHOutCentV0A","fHOutCentV0A; Centrality V0A",505,0,101);
    fHOutCentV0C    = new TH1F("fHOutCentV0C","fHOutCentV0C; Centrality V0C",505,0,101);
    fHOutCentV0MEq     = new TH1F("fHOutCentV0MEq","fHOutCentV0MEq; Centrality V0 equalized",505,0,101);
    fHOutCentV0AEq    = new TH1F("fHOutCentV0AEq","fHOutCentV0AEq; Centrality V0A equalized",505,0,101);
    fHOutCentV0CEq    = new TH1F("fHOutCentV0CEq","fHOutCentV0CEq; Centrality V0C equalized",505,0,101);
    fHOutCentV0MCVHN= new TH1F("fHOutCentV0M_CVHN","fHOutCentV0M_CVHN; Centrality V0",505,0,101);
    fHOutCentV0MCVLN= new TH1F("fHOutCentV0M_CVLN","fHOutCentV0M_CVLN; Centrality V0",505,0,101);
    fHOutCentV0MCVHNinMB= new TH1F("fHOutCentV0M_CVHNinMB","fHOutCentV0M_CVHN; Centrality V0",505,0,101);
    fHOutCentV0MCVLNinMB= new TH1F("fHOutCentV0M_CVLNinMB","fHOutCentV0M_CVLN; Centrality V0",505,0,101);
    fHOutCentV0MCCENT= new TH1F("fHOutCentV0M_CCENT","fHOutCentV0M_CCENT; Centrality V0",505,0,101);
    fHOutCentV0MCSEMI= new TH1F("fHOutCentV0M_CSEMI","fHOutCentV0M_CSEMI; Centrality V0",505,0,101);
    fHOutCentV0MCCENTinMB= new TH1F("fHOutCentV0M_CCENTinMB","fHOutCentV0M_CCENT; Centrality V0",505,0,101);
    fHOutCentV0MCSEMIinMB= new TH1F("fHOutCentV0M_CSEMIinMB","fHOutCentV0M_CSEMI; Centrality V0",505,0,101);
    fHOutCentV0MMSL= new TH1F("fHOutCentV0M_MSL","fHOutCentV0M_MSL; Centrality V0",505,0,101);
    fHOutCentV0MMSH= new TH1F("fHOutCentV0M_MSH","fHOutCentV0M_MSH; Centrality V0",505,0,101);
    fHOutCentV0MMUL= new TH1F("fHOutCentV0M_MUL","fHOutCentV0M_MUL; Centrality V0",505,0,101);
    fHOutCentV0MMLL= new TH1F("fHOutCentV0M_MLL","fHOutCentV0M_MLL; Centrality V0",505,0,101);
    fHOutCentV0MEJE= new TH1F("fHOutCentV0M_EJE","fHOutCentV0M_EJE; Centrality V0",505,0,101);
    fHOutCentV0MEGA= new TH1F("fHOutCentV0M_EGA","fHOutCentV0M_EGA; Centrality V0",505,0,101);
    fHOutCentV0MPHS= new TH1F("fHOutCentV0M_PHS","fHOutCentV0M_PHS; Centrality V0",505,0,101);
    fHOutCentV0MMSLinMB= new TH1F("fHOutCentV0M_MSLinMB","fHOutCentV0M_MSLinMB; Centrality V0",505,0,101);
    fHOutCentV0MMSHinMB= new TH1F("fHOutCentV0M_MSHinMB","fHOutCentV0M_MSHinMB; Centrality V0",505,0,101);
    fHOutCentV0MMULinMB= new TH1F("fHOutCentV0M_MULinMB","fHOutCentV0M_MULinMB; Centrality V0",505,0,101);
    fHOutCentV0MMLLinMB= new TH1F("fHOutCentV0M_MLLinMB","fHOutCentV0M_MLLinMB; Centrality V0",505,0,101);
    fHOutCentV0MEJEinMB= new TH1F("fHOutCentV0M_EJEinMB","fHOutCentV0M_EJEinMB; Centrality V0",505,0,101);
    fHOutCentV0MEGAinMB= new TH1F("fHOutCentV0M_EGAinMB","fHOutCentV0M_EGAinMB; Centrality V0",505,0,101);
    fHOutCentV0MPHSinMB= new TH1F("fHOutCentV0M_PHSinMB","fHOutCentV0M_PHSinMB; Centrality V0",505,0,101);
    fHOutCentFMD     = new TH1F("fHOutCentFMD","fHOutCentFMD; Centrality FMD",505,0,101);
    fHOutCentTRK     = new TH1F("fHOutCentTRK","fHOutCentTRK; Centrality TPC",505,0,101);
    fHOutCentTKL     = new TH1F("fHOutCentTKL","fHOutCentTKL; Centrality tracklets",505,0,101);
    fHOutCentCL0     = new TH1F("fHOutCentCL0","fHOutCentCL0; Centrality SPD inner",505,0,101);
    fHOutCentCL1     = new TH1F("fHOutCentCL1","fHOutCentCL1; Centrality SPD outer",505,0,101);
    fHOutCentCND     = new TH1F("fHOutCentCND","fHOutCentCND; Centrality candle",505,0,101);
    fHOutCentNPA     = new TH1F("fHOutCentNPA","fHOutCentNPA; Centrality Npart",505,0,101);
    fHOutCentZNA     = new TH1F("fHOutCentZNA","fHOutCentZNA; Centrality ZNA",505,0,101);
    fHOutCentZNC     = new TH1F("fHOutCentZNC","fHOutCentZNC; Centrality ZNC",505,0,101);
    fHOutCentV0MvsFMD= new TH1F("fHOutCentV0MvsFMD","fHOutCentV0MvsFMD; Centrality V0 vs FMD",505,0,101);
    fHOutCentTKLvsV0M= new TH1F("fHOutCentTKLvsV0M","fHOutCentTKLvsV0M; Centrality tracklets vs V0",505,0,101);
    fHOutCentZEMvsZDC= new TH1F("fHOutCentZEMvsZDC","fHOutCentZEMvsZDC; Centrality ZEM vs ZDC",505,0,101);
    fHOutCentV0MvsCentCL1= new TH2F("fHOutCentV0MvsCentCL1","fHOutCentV0MvsCentCL1; Cent V0; Cent SPD",505,0,101,505,0,101);
    fHOutCentV0MvsCentTRK= new TH2F("fHOutCentV0MvsCentTRK","fHOutCentV0MvsCentTRK; Cent V0; Cent TPC",505,0,101,505,0,101);
    fHOutCentTRKvsCentCL1= new TH2F("fHOutCentTRKvsCentCL1","fHOutCentTRKvsCentCL1; Cent TPC; Cent SPD",505,0,101,505,0,101);
    fHOutCentV0MvsCentZDC= new TH2F("fHOutCentV0MvsCentZDC","fHOutCentV0MvsCentZDC; Cent V0; Cent ZDC",505,0,101,505,0,101);
    fHOutCentV0AvsCentV0C= new TH2F("fHOutCentV0AvsCentV0C","fHOutCentV0AvsCentV0C; Cent V0A; Cent V0C;", 505,0,101,505,0,101);
    fHOutCentV0AvsCentTRK= new TH2F("fHOutCentV0AvsCentTRK","fHOutCentV0AvsCentTRK; Cent V0A; Cent TRK;", 505,0,101,505,0,101);
    fHOutCentV0AvsCentCND= new TH2F("fHOutCentV0AvsCentCND","fHOutCentV0AvsCentCND; Cent V0A; Cent CND;", 505,0,101,505,0,101);
    fHOutCentV0AvsCentCL1= new TH2F("fHOutCentV0AvsCentCL1","fHOutCentV0AvsCentCL1; Cent V0A; Cent CL1;", 505,0,101,505,0,101);
    fHOutCentV0CvsCentTRK= new TH2F("fHOutCentV0CvsCentTRK","fHOutCentV0CvsCentTRK; Cent V0C; Cent TRK;", 505,0,101,505,0,101);
    fHOutCentV0CvsCentCND= new TH2F("fHOutCentV0CvsCentCND","fHOutCentV0CvsCentCND; Cent V0C; Cent CND;", 505,0,101,505,0,101);
    fHOutCentV0CvsCentCL1= new TH2F("fHOutCentV0CvsCentCL1","fHOutCentV0CvsCentCL1; Cent V0C; Cent CL1;", 505,0,101,505,0,101);
    fHOutCentNPAvsCentV0A= new TH2F("fHOutCentNPAvsCentV0A","fHOutCentNPAvsCentV0A; Cent NPA; Cent V0A;", 505,0,101,505,0,101);
    fHOutCentNPAvsCentV0C= new TH2F("fHOutCentNPAvsCentV0C","fHOutCentNPAvsCentV0C; Cent NPA; Cent V0C;", 505,0,101,505,0,101);
    fHOutCentNPAvsCentTRK= new TH2F("fHOutCentNPAvsCentTRK","fHOutCentNPAvsCentTRK; Cent NPA; Cent TRK;", 505,0,101,505,0,101);
    fHOutCentNPAvsCentCND= new TH2F("fHOutCentNPAvsCentCND","fHOutCentNPAvsCentCND; Cent NPA; Cent CND;", 505,0,101,505,0,101);
    fHOutCentNPAvsCentCL1= new TH2F("fHOutCentNPAvsCentCL1","fHOutCentNPAvsCentCL1; Cent NPA; Cent CL1;", 505,0,101,505,0,101);
    fHOutCentZNAvsCentV0A= new TH2F("fHOutCentZNAvsCentV0A","fHOutCentZNAvsCentV0A; Cent ZNA; Cent V0A;", 505,0,101,505,0,101);
    fHOutCentZNAvsCentV0C= new TH2F("fHOutCentZNAvsCentV0C","fHOutCentZNAvsCentV0C; Cent ZNA; Cent V0C;", 505,0,101,505,0,101);
    fHOutCentZNAvsCentTRK= new TH2F("fHOutCentZNAvsCentTRK","fHOutCentZNAvsCentTRK; Cent ZNA; Cent TRK;", 505,0,101,505,0,101);
    fHOutCentZNAvsCentCND= new TH2F("fHOutCentZNAvsCentCND","fHOutCentZNAvsCentCND; Cent ZNA; Cent CND;", 505,0,101,505,0,101);
    fHOutCentZNAvsCentCL1= new TH2F("fHOutCentZNAvsCentCL1","fHOutCentZNAvsCentCL1; Cent ZNA; Cent CL1;", 505,0,101,505,0,101);

    fHOutMultV0AC = new TH2F("fHOutMultV0AC","fHOutMultV0AC; Multiplicity V0A; Multiplicity V0C",1000,0,1000,1000,0,1000);
    fHOutMultV0M  = new TH1F("fHOutMultV0M","fHOutMultV0M; Multiplicity V0",25000,0,25000);
    fHOutMultV0A  = new TH1F("fHOutMultV0A","fHOutMultV0A; Multiplicity V0",25000,0,25000);
    fHOutMultV0C  = new TH1F("fHOutMultV0C","fHOutMultV0C; Multiplicity V0",25000,0,25000);
    fHOutMultV0Mnc= new TH1F("fHOutMultV0Mnc","fHOutMultV0Mnc; Multiplicity V0",25000,0,25000);
    fHOutMultV0Anc= new TH1F("fHOutMultV0Anc","fHOutMultV0Anc; Multiplicity V0",25000,0,25000);
    fHOutMultV0Cnc= new TH1F("fHOutMultV0Cnc","fHOutMultV0Cnc; Multiplicity V0",25000,0,25000);
    fHOutMultV0O  = new TH1F("fHOutMultV0O","fHOutMultV0O; Multiplicity V0",40000,0,40000);
    fHOutMultV0Cells = new TH2F("fHOutMultV0Cells","fHOutMultV0Cells",33,-0.5,32.5,33,-0.5,32.5); 
    fHOutMultFMD = new TH1F("fHOutMultFMD","fHOutMultFMD; Multiplicity FMD",24000,0,24000);
    fHOutMultTRK = new TH1F("fHOutMultTRK","fHOutMultTRK; Multiplicity TPC",4000,0,4000);
    fHOutMultTKL = new TH1F("fHOutMultTKL","fHOutMultTKL; Multiplicity tracklets",5000,0,5000);
    fHOutMultCL0 = new TH1F("fHOutMultCL0","fHOutMultCL0; Multiplicity SPD inner",7000,0,7000);
    fHOutMultCL1 = new TH1F("fHOutMultCL1","fHOutMultCL1; Multiplicity SPD outer",7000,0,7000);
    fHOutMultCND = new TH1F("fHOutMultCND","fHOutMultCND; Multiplicity candle",4000,0,4000);
    fHOutMultNPA = new TH1F("fHOutMultNPA","fHOutMultNPA; Nparticipants",450,0,450);
    fHOutMultZNA = new TH1F("fHOutMultZNA","fHOutMultZNA; ZNA Energy",500,0,2000);

    fHOutMultV0MvsZDN = new TH2F("fHOutMultV0MvsZDN","fHOutMultV0MvsZDN; Multiplicity V0; Energy ZDC-N",500,0,30000,500,0,180000);
    fHOutMultZEMvsZDN = new TH2F("fHOutMultZEMvsZDN","fHOutMultZEMvsZDN; Energy ZEM; Energy ZDC-N",500,0,2500,500,0,180000);
    fHOutMultV0MvsZDC = new TH2F("fHOutMultV0MvsZDC","fHOutMultV0MvsZDC; Multiplicity V0; Energy ZDC",500,0,30000,500,0,200000);
    fHOutMultZEMvsZDC = new TH2F("fHOutMultZEMvsZDC","fHOutMultZEMvsZDC; Energy ZEM; Energy ZDC",500,0,2500,500,0,200000);
    fHOutMultZEMvsZDCw = new TH2F("fHOutMultZEMvsZDCw","fHOutMultZEMvsZDCw; Energy ZEM; Energy ZDC (weigthed with V0 percentile)",500,0,2500,500,0,200000);
    fHOutMultV0MvsCL1 = new TH2F("fHOutMultV0MvsCL1","fHOutMultV0MvsCL1; Multiplicity V0; Multiplicity SPD outer",2500,0,30000,700,0,7000);
    fHOutMultV0MvsTRK = new TH2F("fHOutMultV0MvsTRK","fHOutMultV0MvsTRK; Multiplicity V0; Multiplicity TPC",2500,0,30000,400,0,4000);
    fHOutMultTRKvsCL1 = new TH2F("fHOutMultTRKvsCL1","fHOutMultTRKvsCL1; Multiplicity TPC; Multiplicity SPD outer",400,0,4000,700,0,7000);
    fHOutMultV0MvsV0O = new TH2F("fHOutMultV0MvsV0O","fHOutMultV0MvsV0O; Multiplicity V0; Multiplicity V0 Online",500,0,30000,500,0,40000);
    fHOutMultV0OvsCL1 = new TH2F("fHOutMultV0OvsCL1","fHOutMultV0OvsCL1; Multiplicity V0; Multiplicity SPD outer",500,0,40000,700,0,7000);
    fHOutMultV0OvsTRK = new TH2F("fHOutMultV0OvsTRK","fHOutMultV0OvsTRK; Multiplicity V0; Multiplicity TPC",500,0,40000,400,0,4000);
    fHOutMultV0MvsV0O = new TH2F("fHOutMultV0MvsV0O","fHOutMultV0MvsV0O; Multiplicity V0; Multiplicity V0 Online",500,0,30000,500,0,30000);
    fHOutMultV0OvsCL1 = new TH2F("fHOutMultV0OvsCL1","fHOutMultV0OvsCL1; Multiplicity V0; Multiplicity SPD outer",2500,0,30000,700,0,7000);
    fHOutMultV0OvsTRK = new TH2F("fHOutMultV0OvsTRK","fHOutMultV0OvsTRK; Multiplicity V0; Multiplicity TPC",2500,0,30000,400,0,4000);
    fHOutMultCL1vsTKL = new TH2F ("fHOutMultCL1vsTKL","fHOutMultCL1vsTKL; Multiplicity SPD outer; Multiplicity tracklets",700,0,7000,700,0,7000);
    
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
    fHOutVertexT0  = new TH1F("fHOutVertexT0", "fHOutVertexT0", 100,-20,20);
  
    fOutputList->Add(fHOutCentV0M);
    fOutputList->Add(fHOutCentV0A);
    fOutputList->Add(fHOutCentV0C);
    fOutputList->Add(fHOutCentV0MEq);
    fOutputList->Add(fHOutCentV0AEq);
    fOutputList->Add(fHOutCentV0CEq);
    fOutputList->Add(fHOutCentV0MCVHN);
    fOutputList->Add(fHOutCentV0MCVLN);
    fOutputList->Add(fHOutCentV0MCVHNinMB);
    fOutputList->Add(fHOutCentV0MCVLNinMB);
    fOutputList->Add(fHOutCentV0MCCENT);
    fOutputList->Add(fHOutCentV0MCSEMI);
    fOutputList->Add(fHOutCentV0MCCENTinMB);
    fOutputList->Add(fHOutCentV0MCSEMIinMB);
    fOutputList->Add(fHOutCentV0MMSL);
    fOutputList->Add(fHOutCentV0MMSH);
    fOutputList->Add(fHOutCentV0MMUL);
    fOutputList->Add(fHOutCentV0MMLL);
    fOutputList->Add(fHOutCentV0MEJE);
    fOutputList->Add(fHOutCentV0MEGA);
    fOutputList->Add(fHOutCentV0MPHS);
    fOutputList->Add(fHOutCentV0MMSLinMB);
    fOutputList->Add(fHOutCentV0MMSHinMB);
    fOutputList->Add(fHOutCentV0MMULinMB);
    fOutputList->Add(fHOutCentV0MMLLinMB);
    fOutputList->Add(fHOutCentV0MEJEinMB);
    fOutputList->Add(fHOutCentV0MEGAinMB);
    fOutputList->Add(fHOutCentV0MPHSinMB);
    fOutputList->Add(fHOutCentFMD);
    fOutputList->Add(fHOutCentTRK);
    fOutputList->Add(fHOutCentTKL);
    fOutputList->Add(fHOutCentCL0);
    fOutputList->Add(fHOutCentCL1);
    fOutputList->Add(fHOutCentCND);
    fOutputList->Add(fHOutCentNPA);
    fOutputList->Add(fHOutCentZNA);
    fOutputList->Add(fHOutCentZNC);
    fOutputList->Add(fHOutCentV0MvsFMD);
    fOutputList->Add(fHOutCentTKLvsV0M);
    fOutputList->Add(fHOutCentZEMvsZDC);
    fOutputList->Add(fHOutCentV0MvsCentCL1);
    fOutputList->Add(fHOutCentV0MvsCentTRK);
    fOutputList->Add(fHOutCentTRKvsCentCL1);
    fOutputList->Add(fHOutCentV0MvsCentZDC);
    fOutputList->Add(fHOutCentV0AvsCentV0C);
    fOutputList->Add(fHOutCentV0AvsCentTRK);
    fOutputList->Add(fHOutCentV0AvsCentCND);
    fOutputList->Add(fHOutCentV0AvsCentCL1);
    fOutputList->Add(fHOutCentV0CvsCentTRK);
    fOutputList->Add(fHOutCentV0CvsCentCND);
    fOutputList->Add(fHOutCentV0CvsCentCL1);
    fOutputList->Add(fHOutCentNPAvsCentV0A);
    fOutputList->Add(fHOutCentNPAvsCentV0C);
    fOutputList->Add(fHOutCentNPAvsCentTRK);
    fOutputList->Add(fHOutCentNPAvsCentCND);
    fOutputList->Add(fHOutCentNPAvsCentCL1);
    fOutputList->Add(fHOutCentZNAvsCentV0A);
    fOutputList->Add(fHOutCentZNAvsCentV0C);
    fOutputList->Add(fHOutCentZNAvsCentTRK);
    fOutputList->Add(fHOutCentZNAvsCentCND);
    fOutputList->Add(fHOutCentZNAvsCentCL1);

    fOutputList->Add(fHOutMultV0AC); 
    fOutputList->Add(fHOutMultV0M); 
    fOutputList->Add(fHOutMultV0A); 
    fOutputList->Add(fHOutMultV0C); 
    fOutputList->Add(fHOutMultV0Mnc); 
    fOutputList->Add(fHOutMultV0Anc); 
    fOutputList->Add(fHOutMultV0Cnc); 
    fOutputList->Add(fHOutMultV0O);
    fOutputList->Add(fHOutMultV0Cells) ;   
    fOutputList->Add(fHOutMultFMD); 
    fOutputList->Add(fHOutMultTRK); 
    fOutputList->Add(fHOutMultTKL); 
    fOutputList->Add(fHOutMultCL0); 
    fOutputList->Add(fHOutMultCL1); 
    fOutputList->Add(fHOutMultCND); 
    fOutputList->Add(fHOutMultNPA); 
    fOutputList->Add(fHOutMultZNA); 
    fOutputList->Add(fHOutMultV0MvsZDN);
    fOutputList->Add(fHOutMultZEMvsZDN);
    fOutputList->Add(fHOutMultV0MvsZDC);
    fOutputList->Add(fHOutMultZEMvsZDC);
    fOutputList->Add(fHOutMultZEMvsZDCw);
    fOutputList->Add(fHOutMultV0MvsCL1);
    fOutputList->Add(fHOutMultV0MvsTRK);
    fOutputList->Add(fHOutMultTRKvsCL1);
    fOutputList->Add(fHOutMultV0MvsV0O);
    fOutputList->Add(fHOutMultV0OvsCL1);
    fOutputList->Add(fHOutMultV0OvsTRK);
    fOutputList->Add(fHOutMultCL1vsTKL);
    fOutputList->Add(fHOutCentV0Mqual1);
    fOutputList->Add(fHOutCentTRKqual1);
    fOutputList->Add(fHOutCentCL1qual1);                   
    fOutputList->Add(fHOutMultV0MvsCL1qual1);
    fOutputList->Add(fHOutMultV0MvsTRKqual1);
    fOutputList->Add(fHOutMultTRKvsCL1qual1);
    fOutputList->Add(fHOutCentV0Mqual2);
    fOutputList->Add(fHOutCentTRKqual2);
    fOutputList->Add(fHOutCentCL1qual2);
    fOutputList->Add(fHOutMultV0MvsCL1qual2);
    fOutputList->Add(fHOutMultV0MvsTRKqual2);
    fOutputList->Add(fHOutMultTRKvsCL1qual2);
    fOutputList->Add(fHOutQuality);
    fOutputList->Add(fHOutVertex);
    fOutputList->Add(fHOutVertexT0);
  
    PostData(1, fOutputList); 
  }
  
  fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  // Add SPD requirement
  fEsdTrackCutsExtra1 = new AliESDtrackCuts("SPD", "Require 1 cluster in SPD");
  fEsdTrackCutsExtra1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  // Add SDD requirement
  fEsdTrackCutsExtra2 = new AliESDtrackCuts("SDD", "Require 1 cluster in first layer SDD");
  fEsdTrackCutsExtra2->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kFirst);  
}

//________________________________________________________________________
void AliCentralitySelectionTask::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliCentralitySelectionTask::UserExec() \n");

  Int_t    runType = 0;             // 0:PbPb, 1:pPb or Pbp 

  Float_t  zncEnergy = 0.;          //  ZNC Energy
  Float_t  zpcEnergy = 0.;          //  ZPC Energy
  Float_t  znaEnergy = 0.;          //  ZNA Energy
  Float_t  zpaEnergy = 0.;          //  ZPA Energy
  Float_t  zem1Energy = 0.;         //  ZEM1 Energy
  Float_t  zem2Energy = 0.;         //  ZEM2 Energy
  Bool_t   zdcEnergyCal = kFALSE;   // if zdc is calibrated (in pass2)
  Double_t znaTower = 0.;           // common PMT of ZNA 
  Double_t zncTower = 0.;
  Bool_t   znaFired = kFALSE;
  Bool_t   zncFired = kFALSE;

  Int_t    nTracks = 0;             //  no. tracks
  Int_t    nTracklets = 0;          //  no. tracklets
  Int_t    nClusters[6] = {0};      //  no. clusters on 6 ITS layers
  Int_t    nChips[2] = {0,0};       //  no. chips on 2 SPD layers
  Float_t  spdCorr =0;              //  corrected spd2 multiplicity
  Int_t    multCND = 0;             //  no. tracks (candle condition)

  Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
  Float_t  multV0AEq  = 0;          //  multiplicity from V0 reco side A
  Float_t  multV0CEq  = 0;          //  multiplicity from V0 reco side C
  Float_t  multV0ACorr  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0CCorr  = 0;            //  multiplicity from V0 reco side C
  Short_t  multV0AOnline  = 0;      //  multiplicity from V0 reco side A
  Short_t  multV0COnline  = 0;      //  multiplicity from V0 reco side C
  Float_t  v0Corr = 0;              // corrected V0 multiplicity (used for MC)
  Int_t nV0A = 0;
  Int_t nV0C = 0;

  Float_t  multFMDA = 0;            //  multiplicity from FMD on detector A
  Float_t  multFMDC = 0;            //  multiplicity from FMD on detector C

  Float_t zvtx =0;                  // z-vertex SPD
  Int_t zvtxNcont =0;               // contributors to z-vertex SPD

  Float_t zvtxT0 =0;                // z-vertex T0

  Int_t Npart =0;                   // N. of participants (true MC)

  AliCentrality *esdCent = 0;

  AliVEvent *event = InputEvent();
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aod =  dynamic_cast<AliAODEvent*>(event);
  if(fAnalysisInput.CompareTo("ESD")==0){
    if (!esd) {
      AliError("No ESD Event");
      return;
    }
  } else if(fAnalysisInput.CompareTo("AOD")==0){
    if (!aod) {
      AliError("No AOD Event");
      return;
    }
  }
  LoadBranches();

  if (SetupRun(event)<0) {
    AliError("Centrality File not available for this run");
    return;
  }

  if (esd) {
    if (strcmp(esd->GetESDRun()->GetBeamType(), "A-A") == 0) runType=0;
    else runType=1;  
  } else {
    Int_t runNumber = event->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593))    // LHC11h
      runType=0;
  }

  esdCent = event->GetCentrality();

  // ***** Vertex Info
  if (esd) {
    const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
    zvtx        = vtxESD->GetZ(); 
    zvtxNcont   = vtxESD->GetNContributors();
  } else {
    const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
    zvtx        = spdVtx->GetZ(); 
    zvtxNcont   = spdVtx->GetNContributors();
  }

    // ***** V0 info    
  AliVVZERO* esdV0 = event->GetVZEROData();
  if (!esdV0) {
    AliError("AliVVZERO not available");
    return;
  }

  multV0A=esdV0->GetMTotV0A();
  multV0C=esdV0->GetMTotV0C();

  multV0ACorr = AliESDUtils::GetCorrV0A(multV0A,zvtx);    
  multV0CCorr = AliESDUtils::GetCorrV0C(multV0C,zvtx);    

  v0Corr = multV0A+multV0C; // Todo: C.L. not clear why here we do not use the sum of the corrected values?

  multV0AOnline=esdV0->GetTriggerChargeA(); 
  multV0COnline=esdV0->GetTriggerChargeC(); 

  // Count V0 flags
  for(Int_t i = 0; i < 32; ++i) {
    if (esdV0->GetBBFlag(i)) nV0C++;
    if (esdV0->GetBBFlag(i+32)) nV0A++;
  }
    
  // Equalized signals
  multV0AEq=0.;
  multV0CEq=0.;
  for(Int_t iCh = 4; iCh < 7; ++iCh) {
    Double_t mult = esd->GetVZEROEqMultiplicity(iCh);
    multV0AEq += mult;
  }
  for(Int_t iCh = 0; iCh < 3; ++iCh) {
    Double_t mult = esd->GetVZEROEqMultiplicity(iCh);
    multV0CEq += mult;
  }
  
  Bool_t kT0BB = kFALSE;    
  if (esd) {
    // ***** T0 info    
    const AliESDTZERO* esdT0 = esd->GetESDTZERO();
    if (!esdT0)
      {
	AliError("AliESDTZERO not available");
	return;
      }
    Int_t trig=esdT0->GetT0Trig();
    if(trig&1) kT0BB=kTRUE;
    zvtxT0=esdT0->GetT0zVertex();
  } else {
    const AliAODTZERO* esdT0 = aod->GetTZEROData();
    if (!esdT0)
      {
	AliError("AliAODTZERO not available");
	return;
      }
    Int_t trig=1;//esdT0->GetT0Trig(); //* Todo: C.L. This info is not in AOD? */
    if(trig&1) kT0BB=kTRUE;
    zvtxT0=esdT0->GetT0zVertex();
  }

  // ***** Trigger info    
  fIsSelected = ((esdV0->GetV0ADecision()==1) && (esdV0->GetV0CDecision()==1));
  TString trigStr;
  if (esd)
    trigStr = esd->GetFiredTriggerClasses();
  else
    trigStr = aod->GetFiredTriggerClasses();
    
  fMB=kFALSE;
  fCVHN=kFALSE; fCVLN=kFALSE; fCCENT=kFALSE; fCSEMI=kFALSE; 
  fMSL=kFALSE;  fMSH=kFALSE;  fMUL=kFALSE;   fMLL=kFALSE;
  fEJE=kFALSE;  fEGA=kFALSE;  fPHS=kFALSE;

  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI")) && (fIsSelected)) 
    fMB=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CVHN")) && (fIsSelected)) 
    fCVHN=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CVLN")) && (fIsSelected))
    fCVLN=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CCENT")) && (fIsSelected)) 
    fCCENT=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CSEMI")) && (fIsSelected))
    fCSEMI=kTRUE;
    
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1MSL")) && (fIsSelected))
    fMSL=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1MSH")) && (fIsSelected))
    fMSH=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1MUL")) && (fIsSelected))
    fMUL=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1MLL")) && (fIsSelected))
    fMLL=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1EJE")) && (fIsSelected))
    fEJE=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1EGA")) && (fIsSelected))
    fEGA=kTRUE;
  if ( (trigStr.Contains("-B-")) &&  (trigStr.Contains("CPBI1PHS")) && (fIsSelected))
    fPHS=kTRUE;
  
  fCVHNbit=kFALSE;    fCVLNbit=kFALSE;       fCCENTbit=kFALSE;    fCSEMIbit=kFALSE; 
  if (esdV0->GetTriggerBits() & (1<<8)) 
    fCVHNbit=kTRUE;
  if (esdV0->GetTriggerBits() & (1<<6)) 
    fCVLNbit=kTRUE;
    
  if (kT0BB && fCVHNbit)
    fCCENTbit=kTRUE;
  if (kT0BB && fCVLNbit)
    fCSEMIbit=kTRUE;

  if (esd) {  
    // ***** CB info (tracklets, clusters, chips)
    //nTracks    = event->GetNumberOfTracks();     
    nTracks    = fTrackCuts ? (Short_t)fTrackCuts->GetReferenceMultiplicity(esd,kTRUE):-1;
  } else {
    AliAODHeader *h = aod->GetHeader();
    nTracks    = h!=0 ? (Short_t)h->GetTPConlyRefMultiplicity():-1;
  }

  if (esd) {
    Short_t nTrTPCcandle = 0;
    for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {

      AliESDtrack* track = esd->GetTrack(iTracks);
      if (!track) continue;
      
      if (! fEsdTrackCuts->IsSelected(track) )continue;
      
      if (fEsdTrackCutsExtra1 && fEsdTrackCutsExtra2 &&
       	  !fEsdTrackCutsExtra1->IsSelected(track) &&
       	  !fEsdTrackCutsExtra2->IsSelected(track)) continue;
      
      if (track->Pt() > 0.4 && TMath::Abs(track->Eta()) < 0.9)	nTrTPCcandle++;
    } 
    multCND = nTrTPCcandle;
  } else {
    Short_t nTrTPCcandle = 0;
    for (Int_t iTracks = 0; iTracks < aod->GetNumberOfTracks(); iTracks++) {

      AliAODTrack* track = aod->GetTrack(iTracks);

      if (!track) continue;
      if (!track->TestFilterBit(1<<5) && 
	  !track->TestFilterBit(1<<6)) continue;

      if (track->Pt() > 0.4 && TMath::Abs(track->Eta()) < 0.9)	nTrTPCcandle++;
    }
    multCND = nTrTPCcandle;
  }

  if (esd) {
    const AliMultiplicity *mult = esd->GetMultiplicity();
    nTracklets = mult->GetNumberOfTracklets();

    for(Int_t ilay=0; ilay<6; ilay++){
      nClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
    }

    for(Int_t ilay=0; ilay<2; ilay++){
      nChips[ilay] = mult->GetNumberOfFiredChips(ilay);
    }
  } else {
    AliAODTracklets *mult = aod->GetTracklets();
    nTracklets = mult->GetNumberOfTracklets();
    AliAODHeader *h = aod->GetHeader();
    for(Int_t ilay=0; ilay<6; ilay++){
      nClusters[ilay] = h->GetNumberOfITSClusters(ilay);
    }
  }
  spdCorr = AliESDUtils::GetCorrSPD2(nClusters[1],zvtx);
  
  if (esd) {
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
  }

  if (esd) {
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

    const Double_t *ZNAtower = esdZDC->GetZN2TowerEnergy(); 
    znaTower = ZNAtower[0];
    const Double_t *ZNCtower = esdZDC->GetZN1TowerEnergy();
    zncTower = ZNCtower[0];
    
    for (Int_t j = 0; j < 4; ++j) 
      if (esdZDC->GetZDCTDCData(12,j) != 0) 
	znaFired = kTRUE;
    
    for (Int_t j = 0; j < 4; ++j) 
      if (esdZDC->GetZDCTDCData(10,j) != 0) 
	zncFired = kTRUE;   

  } else {
    AliAODHeader *h = aod->GetHeader();
    zncEnergy = (Float_t) (h->GetZDCN1Energy());
    zpcEnergy = (Float_t) (h->GetZDCP1Energy());
    znaEnergy = (Float_t) (h->GetZDCN2Energy());
    zpaEnergy = (Float_t) (h->GetZDCP2Energy());
    zem1Energy = (Float_t) (h->GetZDCEMEnergy(0))/8.; //Todo: C.L. Should we devide here by 8? It is done in the ESD case!
    zem2Energy = (Float_t) (h->GetZDCEMEnergy(1))/8.;

    AliAODZDC *aodZDC = aod->GetZDCData();
    const Double_t *ZNAtower = aodZDC->GetZNATowerEnergy(); 
    znaTower = ZNAtower[0];

    znaFired = kFALSE; // trick because info is not stored in AOD
    if (esdCent->GetCentralityPercentile("ZNA") != 101)
      znaFired = kTRUE;
  }

  if (esd) {
    // ***** MC info
    AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
    AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
    AliStack*    stack=0;
    AliMCEvent*  mcEvent=0;
    if (fIsMCInput && eventHandler && (mcEvent=eventHandler->MCEvent()) && (stack=mcEvent->Stack())) {
      AliGenHijingEventHeader* hHijing=0;
      AliGenDPMjetEventHeader* dpmHeader=0;
      
      AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
      if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())) 
	hHijing = (AliGenHijingEventHeader*)mcGenH;      
      else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
	TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
	hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing"));
	if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing pPb_0"));
      }
      else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
	dpmHeader = (AliGenDPMjetEventHeader*)mcGenH;
      }
      if(hHijing)   Npart = hHijing->ProjectileParticipants()+hHijing->TargetParticipants();
      if(dpmHeader) Npart = dpmHeader->ProjectileParticipants()+ dpmHeader->TargetParticipants();
    }
  } 

  // ***** Scaling for MC
  if (fIsMCInput) {
    fUseScaling=kFALSE;
    v0Corr  = Short_t((multV0A+multV0C)  * fV0MScaleFactorMC);
    multV0A  = multV0A * fV0MScaleFactorMC;
    multV0C  = multV0C * fV0MScaleFactorMC;
  }
  // ***** Scaling for Data 
  if (fUseScaling) {
    v0Corr  = Short_t(v0Corr /   fV0MScaleFactor);
    spdCorr = spdCorr / fSPDScaleFactor;
    nTracks = Int_t(nTracks /   fTPCScaleFactor);
  }

  // ***** Centrality Selection
  if(fHtempV0M) fCentV0M = fHtempV0M->GetBinContent(fHtempV0M->FindBin((v0Corr)));
  if(fHtempV0A) fCentV0A = fHtempV0A->GetBinContent(fHtempV0A->FindBin((multV0ACorr)));
  if(fHtempV0C) fCentV0C = fHtempV0C->GetBinContent(fHtempV0C->FindBin((multV0CCorr)));
  if(fHtempV0MEq) fCentV0MEq = fHtempV0MEq->GetBinContent(fHtempV0MEq->FindBin((multV0AEq+multV0CEq)));
  if(fHtempV0AEq) fCentV0AEq = fHtempV0AEq->GetBinContent(fHtempV0AEq->FindBin((multV0AEq)));
  if(fHtempV0CEq) fCentV0CEq = fHtempV0CEq->GetBinContent(fHtempV0CEq->FindBin((multV0CEq)));
  if(fHtempFMD) fCentFMD = fHtempFMD->GetBinContent(fHtempFMD->FindBin((multFMDA+multFMDC)));
  if(fHtempTRK) fCentTRK = fHtempTRK->GetBinContent(fHtempTRK->FindBin(nTracks));
  if(fHtempTKL) fCentTKL = fHtempTKL->GetBinContent(fHtempTKL->FindBin(nTracklets));
  if(fHtempCL0) fCentCL0 = fHtempCL0->GetBinContent(fHtempCL0->FindBin(nClusters[0]));
  if(fHtempCL1) fCentCL1 = fHtempCL1->GetBinContent(fHtempCL1->FindBin(spdCorr));
  if(fHtempCND) fCentCND = fHtempCND->GetBinContent(fHtempCND->FindBin(multCND));
  if(fHtempZNA) {
    if(znaFired) fCentZNA = fHtempZNA->GetBinContent(fHtempZNA->FindBin(znaTower));
    else fCentZNA = 101;
  }
  if(fHtempZNC) {
    if(zncFired) fCentZNC = fHtempZNC->GetBinContent(fHtempZNC->FindBin(zncTower));
    else fCentZNC = 101;
  }
  if(fHtempV0MvsFMD) fCentV0MvsFMD = fHtempV0MvsFMD->GetBinContent(fHtempV0MvsFMD->FindBin((multV0A+multV0C)));
  if(fHtempTKLvsV0M) fCentTKLvsV0M = fHtempTKLvsV0M->GetBinContent(fHtempTKLvsV0M->FindBin(nTracklets));
  if(fHtempZEMvsZDC) fCentZEMvsZDC = fHtempZEMvsZDC->GetBinContent(fHtempZEMvsZDC->FindBin(zem1Energy+zem2Energy,zncEnergy+znaEnergy+zpcEnergy+zpaEnergy));

  if(fHtempNPA) fCentNPA = fHtempNPA->GetBinContent(fHtempNPA->FindBin(Npart));
  if(fHtempV0Mtrue) fCentV0Mtrue = fHtempV0Mtrue->GetBinContent(fHtempV0Mtrue->FindBin((multV0ACorr+multV0CCorr)));
  if(fHtempV0Atrue) fCentV0Atrue = fHtempV0Atrue->GetBinContent(fHtempV0Atrue->FindBin((multV0ACorr)));
  if(fHtempV0Ctrue) fCentV0Ctrue = fHtempV0Ctrue->GetBinContent(fHtempV0Ctrue->FindBin((multV0CCorr)));
  if(fHtempV0MEqtrue) fCentV0MEqtrue = fHtempV0MEqtrue->GetBinContent(fHtempV0MEqtrue->FindBin((multV0AEq+multV0CEq)));
  if(fHtempV0AEqtrue) fCentV0AEqtrue = fHtempV0AEqtrue->GetBinContent(fHtempV0AEqtrue->FindBin((multV0AEq)));
  if(fHtempV0CEqtrue) fCentV0CEqtrue = fHtempV0CEqtrue->GetBinContent(fHtempV0CEqtrue->FindBin((multV0CEq)));
  if(fHtempFMDtrue) fCentFMDtrue = fHtempFMDtrue->GetBinContent(fHtempFMDtrue->FindBin((multFMDA+multFMDC)));
  if(fHtempTRKtrue) fCentTRKtrue = fHtempTRKtrue->GetBinContent(fHtempTRKtrue->FindBin(nTracks));
  if(fHtempTKLtrue) fCentTKLtrue = fHtempTKLtrue->GetBinContent(fHtempTKLtrue->FindBin(nTracklets));
  if(fHtempCL0true) fCentCL0true = fHtempCL0true->GetBinContent(fHtempCL0true->FindBin(nClusters[0]));
  if(fHtempCL1true) fCentCL1true = fHtempCL1true->GetBinContent(fHtempCL1true->FindBin(spdCorr));
  if(fHtempCNDtrue) fCentCNDtrue = fHtempCNDtrue->GetBinContent(fHtempCNDtrue->FindBin(multCND));
  if(fHtempZNAtrue) fCentZNAtrue = fHtempZNAtrue->GetBinContent(fHtempZNAtrue->FindBin(znaTower));
  if(fHtempZNCtrue) fCentZNCtrue = fHtempZNCtrue->GetBinContent(fHtempZNCtrue->FindBin(zncTower));
   

  // ***** Cleaning
  if (fUseCleaning) {
    fQuality=0;
    
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
	  (zdcEnergyCal==kFALSE) ) fQuality  += 8;
      if (IsOutlierV0MZDCECal((zncEnergy+znaEnergy+zpcEnergy+zpaEnergy), v0Corr) &&
	  (zdcEnergyCal==kTRUE) ) fQuality  += 8;
    }
  } else {
    fQuality = 0;
  }

      
  if (esdCent) {
    if (aod&&(fDebug>1)) {
      Double_t v0m = esdCent->GetCentralityPercentile("V0M");
      Double_t cl1 = esdCent->GetCentralityPercentile("CL1");
      Double_t trk = esdCent->GetCentralityPercentile("TRK");
      Double_t cnd = esdCent->GetCentralityPercentile("CND");
      Double_t zna = esdCent->GetCentralityPercentile("ZNA");
      printf("AOD: v0m %.2f %.2f (%.2f) cl1 %.2f %.2f (%.2f) trk %.2f %.2f (%.2f) cnd %.2f %.2f (%.2f) zna %.2f %.2f (%.2f)\n", 
	     v0m, fCentV0M, fCentV0M!=0?v0m/fCentV0M:1, cl1, fCentCL1, fCentCL1!=0?cl1/fCentCL1:1, trk, fCentTRK, 
	     fCentTRK!=0?trk/fCentTRK:1, cnd, fCentCND, fCentCND!=0?cnd/fCentCND:1, zna, fCentZNA, fCentZNA!=0?zna/fCentZNA:1);
    }
    esdCent->SetQuality(fQuality);
    esdCent->SetCentralityV0M(fCentV0M);
    esdCent->SetCentralityV0A(fCentV0A);
    esdCent->SetCentralityV0C(fCentV0C);
    esdCent->SetCentralityV0MEq(fCentV0MEq);
    esdCent->SetCentralityV0AEq(fCentV0AEq);
    esdCent->SetCentralityV0CEq(fCentV0CEq);
    esdCent->SetCentralityFMD(fCentFMD);
    esdCent->SetCentralityTRK(fCentTRK);
    esdCent->SetCentralityTKL(fCentTKL);
    esdCent->SetCentralityCL0(fCentCL0);
    esdCent->SetCentralityCL1(fCentCL1);
    esdCent->SetCentralityCND(fCentCND);
    esdCent->SetCentralityNPA(fCentNPA);
    esdCent->SetCentralityZNA(fCentZNA);
    esdCent->SetCentralityZNC(fCentZNC);
    esdCent->SetCentralityV0MvsFMD(fCentV0MvsFMD);
    esdCent->SetCentralityTKLvsV0M(fCentTKLvsV0M);
    esdCent->SetCentralityZEMvsZDC(fCentZEMvsZDC);
  }

  // filling QA histograms
  if (fFillHistos) {    

    if (fIsMCInput) { // fill histo with true centrality for simulations
      fCentV0M = fCentV0Mtrue;        
      fCentV0A = fCentV0Atrue;        
      fCentV0C = fCentV0Ctrue;        
      fCentV0MEq = fCentV0MEqtrue;        
      fCentV0AEq = fCentV0AEqtrue;        
      fCentV0CEq = fCentV0CEqtrue;        
      fCentFMD = fCentFMDtrue;        
      fCentTRK = fCentTRKtrue;        
      fCentTKL = fCentTKLtrue;        
      fCentCL0 = fCentCL0true;        
      fCentCL1 = fCentCL1true;        
      fCentCND = fCentCNDtrue;        
      fCentZNA = fCentZNAtrue;        
      fCentZNC = fCentZNCtrue;        
    }


    if ((fMB) && (abs(zvtx)<10))	fHOutMultCL1vsTKL->Fill(spdCorr,nTracklets);

    if (fCVHN)   fHOutCentV0MCVHN->Fill(fCentV0M);
    if (fCVLN)   fHOutCentV0MCVLN->Fill(fCentV0M);
    if (fCCENT)  fHOutCentV0MCCENT->Fill(fCentV0M);
    if (fCSEMI)  fHOutCentV0MCSEMI->Fill(fCentV0M);
    if (fMSL) fHOutCentV0MMSL->Fill(fCentV0M);
    if (fMSH) fHOutCentV0MMSH->Fill(fCentV0M);
    if (fMUL) fHOutCentV0MMUL->Fill(fCentV0M);
    if (fMLL) fHOutCentV0MMLL->Fill(fCentV0M);
    if (fEJE) fHOutCentV0MEJE->Fill(fCentV0M);
    if (fEGA) fHOutCentV0MEGA->Fill(fCentV0M);
    if (fPHS) fHOutCentV0MPHS->Fill(fCentV0M);

    if (((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB) && (runType==0)) ||
	((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7) && (runType==1))) { // fill the QA histograms only for MB events!
      
      fHOutQuality->Fill(fQuality);
      fHOutVertex->Fill(zvtx);
      fHOutVertexT0->Fill(zvtxT0);
      
      if (fQuality==0) {  
	fHOutCentV0M->Fill(fCentV0M);
	fHOutCentV0A->Fill(fCentV0A);
	fHOutCentV0C->Fill(fCentV0C);
	fHOutCentV0MEq->Fill(fCentV0MEq);
	fHOutCentV0AEq->Fill(fCentV0AEq);
	fHOutCentV0CEq->Fill(fCentV0CEq);
	
	if (fCVHNbit)  fHOutCentV0MCVHNinMB->Fill(fCentV0M);
	if (fCVLNbit)  fHOutCentV0MCVLNinMB->Fill(fCentV0M);
	if (fCCENTbit) fHOutCentV0MCCENTinMB->Fill(fCentV0M);
	if (fCSEMIbit) fHOutCentV0MCSEMIinMB->Fill(fCentV0M);
	if (fMSL) fHOutCentV0MMSLinMB->Fill(fCentV0M);
	if (fMSH) fHOutCentV0MMSHinMB->Fill(fCentV0M);
	if (fMUL) fHOutCentV0MMULinMB->Fill(fCentV0M);
	if (fMLL) fHOutCentV0MMLLinMB->Fill(fCentV0M);
	if (fEJE) fHOutCentV0MEJEinMB->Fill(fCentV0M);
	if (fEGA) fHOutCentV0MEGAinMB->Fill(fCentV0M);
	if (fPHS) fHOutCentV0MPHSinMB->Fill(fCentV0M);
       	
	fHOutCentFMD->Fill(fCentFMD);
	fHOutCentTRK->Fill(fCentTRK);
	fHOutCentTKL->Fill(fCentTKL);
	fHOutCentCL0->Fill(fCentCL0);
	fHOutCentCL1->Fill(fCentCL1);
	fHOutCentCND->Fill(fCentCND);
	fHOutCentNPA->Fill(fCentNPA);
	fHOutCentZNA->Fill(fCentZNA);
	fHOutCentZNC->Fill(fCentZNC);
	fHOutCentV0MvsFMD->Fill(fCentV0MvsFMD);
	fHOutCentTKLvsV0M->Fill(fCentTKLvsV0M);
	fHOutCentZEMvsZDC->Fill(fCentZEMvsZDC);
	fHOutCentV0MvsCentCL1->Fill(fCentV0M,fCentCL1);
	fHOutCentV0MvsCentTRK->Fill(fCentV0M,fCentTRK);
	fHOutCentTRKvsCentCL1->Fill(fCentTRK,fCentCL1);
	fHOutCentV0MvsCentZDC->Fill(fCentV0M,fCentZEMvsZDC);
	fHOutCentV0AvsCentV0C->Fill(fCentV0A,fCentV0C);
	fHOutCentV0AvsCentTRK->Fill(fCentV0A,fCentTRK);
	fHOutCentV0AvsCentCND->Fill(fCentV0A,fCentCND);
	fHOutCentV0AvsCentCL1->Fill(fCentV0A,fCentCL1);
	fHOutCentV0CvsCentTRK->Fill(fCentV0C,fCentTRK);
	fHOutCentV0CvsCentCND->Fill(fCentV0C,fCentCND);
	fHOutCentV0CvsCentCL1->Fill(fCentV0C,fCentCL1);
	fHOutCentNPAvsCentV0A->Fill(fCentNPA,fCentV0A);
	fHOutCentNPAvsCentV0C->Fill(fCentNPA,fCentV0C);
	fHOutCentNPAvsCentTRK->Fill(fCentNPA,fCentTRK);
	fHOutCentNPAvsCentCND->Fill(fCentNPA,fCentCND);
	fHOutCentNPAvsCentCL1->Fill(fCentNPA,fCentCL1);
	fHOutCentZNAvsCentV0A->Fill(fCentZNA,fCentV0A);
	fHOutCentZNAvsCentV0C->Fill(fCentZNA,fCentV0C);
	fHOutCentZNAvsCentTRK->Fill(fCentZNA,fCentTRK);
	fHOutCentZNAvsCentCND->Fill(fCentZNA,fCentCND);
	fHOutCentZNAvsCentCL1->Fill(fCentZNA,fCentCL1);

	fHOutMultV0AC->Fill(multV0A,multV0C);
	fHOutMultV0M->Fill(multV0ACorr+multV0CCorr);
	fHOutMultV0A->Fill(multV0ACorr);
	fHOutMultV0C->Fill(multV0CCorr);
	fHOutMultV0MEq->Fill(multV0AEq+multV0CEq);
	fHOutMultV0AEq->Fill(multV0AEq);
	fHOutMultV0CEq->Fill(multV0CEq);
	fHOutMultV0Mnc->Fill(multV0A+multV0C);
	fHOutMultV0Anc->Fill(multV0A);
	fHOutMultV0Cnc->Fill(multV0C);
	fHOutMultV0O->Fill(multV0AOnline+multV0COnline);
	fHOutMultV0Cells->Fill(nV0A,nV0C); 
	fHOutMultFMD->Fill(multFMDA+multFMDC);
	fHOutMultTRK->Fill(nTracks);
	fHOutMultTKL->Fill(nTracklets);
	fHOutMultCL0->Fill(nClusters[0]);
	fHOutMultCL1->Fill(spdCorr);
	fHOutMultCND->Fill(multCND);
	fHOutMultNPA->Fill(Npart);
	if(znaFired)fHOutMultZNA->Fill(znaTower);
	if(zncFired)fHOutMultZNC->Fill(zncTower);

	fHOutMultV0MvsZDN->Fill(v0Corr,(zncEnergy+znaEnergy));
	fHOutMultZEMvsZDN->Fill((zem1Energy+zem2Energy),(zncEnergy+znaEnergy));
	fHOutMultV0MvsZDC->Fill(v0Corr,(zncEnergy+znaEnergy+zpcEnergy+zpaEnergy));
	fHOutMultZEMvsZDC->Fill((zem1Energy+zem2Energy),(zncEnergy+znaEnergy+zpcEnergy+zpaEnergy));
	fHOutMultZEMvsZDCw->Fill((zem1Energy+zem2Energy),(zncEnergy+znaEnergy+zpcEnergy+zpaEnergy),fCentV0M);
	fHOutMultV0MvsCL1->Fill(v0Corr,spdCorr);
	fHOutMultV0MvsTRK->Fill(v0Corr,nTracks);
	fHOutMultTRKvsCL1->Fill(nTracks,spdCorr);
	fHOutMultV0MvsV0O->Fill(v0Corr,(multV0AOnline+multV0COnline));
	fHOutMultV0OvsCL1->Fill((multV0AOnline+multV0COnline),spdCorr);
	fHOutMultV0OvsTRK->Fill((multV0AOnline+multV0COnline),nTracks);
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
    }
    PostData(1, fOutputList); 
  }
}
//________________________________________________________________________
void AliCentralitySelectionTask::Terminate(Option_t */*option*/)
{
  // Terminate analysis
}
//________________________________________________________________________
Int_t AliCentralitySelectionTask::SetupRun(const AliVEvent* const esd)
{
  // Setup files for run

  if (!esd)
    return -1;

  // check if something to be done
  if (fCurrentRun == esd->GetRunNumber())
    return 0;
  else
    fCurrentRun = esd->GetRunNumber();

  TString fileName =(Form("%s/COMMON/CENTRALITY/data/centrality.root", AliAnalysisManager::GetOADBPath()));
  AliInfo(Form("Setup Centrality Selection for run %d with file %s\n",fCurrentRun,fileName.Data()));

  AliOADBContainer *con = new AliOADBContainer("OADB");
  con->InitFromFile(fileName,"Centrality");

  AliOADBCentrality*  centOADB = 0;
  centOADB = (AliOADBCentrality*)(con->GetObject(fCurrentRun));
  if (!centOADB) {
    AliWarning(Form("Centrality OADB does not exist for run %d, using Default \n",fCurrentRun ));
    centOADB  = (AliOADBCentrality*)(con->GetDefaultObject("oadbDefault"));
  }

  Bool_t isHijing=kFALSE;
  Bool_t isDpmjet=kFALSE;
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
  AliMCEvent*  mcEvent=0;
  if (fIsMCInput && eventHandler && (mcEvent=eventHandler->MCEvent()) ) {     
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())) isHijing=kTRUE;
    else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) isHijing=kTRUE;
    else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) isDpmjet=kTRUE;
  }


  // modes
  fUseScaling   = centOADB->UseScaling();
  fUseCleaning  = centOADB->UseCleaning();

  // cuts
  fZVCut        = centOADB->ZVCut();
  fOutliersCut  = centOADB->OutliersCut(); 

  // centrality histos
  fHtempV0M       = centOADB->V0hist(); 
  fHtempV0A       = centOADB->V0Ahist(); 
  fHtempV0C       = centOADB->V0Chist(); 
  fHtempV0MEq     = centOADB->V0Eqhist(); 
  fHtempV0AEq     = centOADB->V0AEqhist(); 
  fHtempV0CEq     = centOADB->V0CEqhist(); 
  fHtempTRK       = centOADB->TPChist();
  fHtempCL1       = centOADB->SPDhist();
  fHtempCND       = centOADB->CNDhist();
  fHtempFMD       = centOADB->FMDhist();
  fHtempZNA       = centOADB->ZNAhist();
  fHtempZNC       = centOADB->ZNChist();
  fHtempZEMvsZDC  = centOADB->ZEMvsZDChist();

   if (isHijing) {       
     fHtempNPA           = centOADB->NPAhist();
     fHtempV0Mtrue       = centOADB->V0histtrue(); 
     fHtempV0Atrue       = centOADB->V0Ahisttrue(); 
     fHtempV0Ctrue       = centOADB->V0Chisttrue(); 
     fHtempV0MEqtrue     = centOADB->V0Eqhisttrue(); 
     fHtempV0AEqtrue     = centOADB->V0AEqhisttrue(); 
     fHtempV0CEqtrue     = centOADB->V0CEqhisttrue(); 
     fHtempTRKtrue       = centOADB->TPChisttrue();
     fHtempCL1true       = centOADB->SPDhisttrue();
     fHtempCNDtrue       = centOADB->CNDhisttrue();
     fHtempFMDtrue       = centOADB->FMDhisttrue();
     fHtempZNAtrue       = centOADB->ZNAhisttrue();
     fHtempZNCtrue       = centOADB->ZNChisttrue();
   }   else if (isDpmjet)   {
     fHtempNPA           = centOADB->NPAhistDPM();
     fHtempV0Mtrue       = centOADB->V0histtrueDPM(); 
     fHtempV0Atrue       = centOADB->V0AhisttrueDPM(); 
     fHtempV0Ctrue       = centOADB->V0ChisttrueDPM(); 
     fHtempV0MEqtrue     = centOADB->V0EqhisttrueDPM(); 
     fHtempV0AEqtrue     = centOADB->V0AEqhisttrueDPM(); 
     fHtempV0CEqtrue     = centOADB->V0CEqhisttrueDPM(); 
     fHtempTRKtrue       = centOADB->TPChisttrueDPM();
     fHtempCL1true       = centOADB->SPDhisttrueDPM();
     fHtempCNDtrue       = centOADB->CNDhisttrueDPM();
     fHtempFMDtrue       = centOADB->FMDhisttrueDPM();
     fHtempZNAtrue       = centOADB->ZNAhisttrueDPM();
     fHtempZNCtrue       = centOADB->ZNChisttrueDPM();
   }


    TString path = gSystem->ExpandPathName(fileName.Data());
  if (!fHtempV0M) AliWarning(Form("Calibration for V0M does not exist in %s", path.Data()));
  if (!fHtempV0A) AliWarning(Form("Calibration for V0A does not exist in %s", path.Data()));
  if (!fHtempV0C) AliWarning(Form("Calibration for V0C does not exist in %s", path.Data()));
  if (!fHtempV0MEq) AliWarning(Form("Calibration for V0MEq does not exist in %s", path.Data()));
  if (!fHtempV0AEq) AliWarning(Form("Calibration for V0AEq does not exist in %s", path.Data()));
  if (!fHtempV0CEq) AliWarning(Form("Calibration for V0CEq does not exist in %s", path.Data()));
  if (!fHtempTRK) AliWarning(Form("Calibration for TRK does not exist in %s", path.Data()));
  if (!fHtempCL1) AliWarning(Form("Calibration for CL1 does not exist in %s", path.Data()));
  if (!fHtempCND) AliWarning(Form("Calibration for CND does not exist in %s", path.Data()));
  if (!fHtempZNA) AliWarning(Form("Calibration for ZNA does not exist in %s", path.Data()));
  if (!fHtempZNC) AliWarning(Form("Calibration for ZNC does not exist in %s", path.Data()));
  if (!fHtempFMD) AliWarning(Form("Calibration for FMD does not exist in %s", path.Data()));
  if (!fHtempZEMvsZDC) AliWarning(Form("Calibration for ZEMvsZDC does not exist in %s", path.Data()));
  if (!fHtempNPA) AliWarning(Form("Calibration for NPA does not exist in %s", path.Data()));

  if (!fHtempV0Mtrue) AliWarning(Form("Calibration for V0Mtrue does not exist in %s", path.Data()));
  if (!fHtempV0Atrue) AliWarning(Form("Calibration for V0Atrue does not exist in %s", path.Data()));
  if (!fHtempV0Ctrue) AliWarning(Form("Calibration for V0Ctrue does not exist in %s", path.Data()));
  if (!fHtempV0MEqtrue) AliWarning(Form("Calibration for V0MEqtrue does not exist in %s", path.Data()));
  if (!fHtempV0AEqtrue) AliWarning(Form("Calibration for V0AEqtrue does not exist in %s", path.Data()));
  if (!fHtempV0CEqtrue) AliWarning(Form("Calibration for V0CEqtrue does not exist in %s", path.Data()));
  if (!fHtempTRKtrue) AliWarning(Form("Calibration for TRKtrue does not exist in %s", path.Data()));
  if (!fHtempCL1true) AliWarning(Form("Calibration for CL1true does not exist in %s", path.Data()));
  if (!fHtempCNDtrue) AliWarning(Form("Calibration for CNDtrue does not exist in %s", path.Data()));
  if (!fHtempZNAtrue) AliWarning(Form("Calibration for ZNAtrue does not exist in %s", path.Data()));
  if (!fHtempZNCtrue) AliWarning(Form("Calibration for ZNCtrue does not exist in %s", path.Data()));
  if (!fHtempFMDtrue) AliWarning(Form("Calibration for FMDtrue does not exist in %s", path.Data()));


  // scale factors
  fV0MScaleFactor    = centOADB->V0MScaleFactor();
  fSPDScaleFactor    = centOADB->SPDScaleFactor();
  fTPCScaleFactor    = centOADB->TPCScaleFactor();
  fV0MScaleFactorMC  = centOADB->V0MScaleFactorMC();

  // outliers parameters
  fV0MSPDOutlierPar0 = centOADB->V0MSPDOutlierPar0();      
  fV0MSPDOutlierPar1 = centOADB->V0MSPDOutlierPar1();     
  fV0MTPCOutlierPar0 = centOADB->V0MTPCOutlierPar0();      
  fV0MTPCOutlierPar1 = centOADB->V0MTPCOutlierPar1();     
  			   			   
  fV0MSPDSigmaOutlierPar0 = centOADB->V0MSPDSigmaOutlierPar0(); 
  fV0MSPDSigmaOutlierPar1 = centOADB->V0MSPDSigmaOutlierPar1(); 
  fV0MSPDSigmaOutlierPar2 = centOADB->V0MSPDSigmaOutlierPar2();
  fV0MTPCSigmaOutlierPar0 = centOADB->V0MTPCSigmaOutlierPar0(); 
  fV0MTPCSigmaOutlierPar1 = centOADB->V0MTPCSigmaOutlierPar1(); 
  fV0MTPCSigmaOutlierPar2 = centOADB->V0MTPCSigmaOutlierPar2(); 
  			    
  fV0MZDCOutlierPar0 =      centOADB->V0MZDCOutlierPar0();      
  fV0MZDCOutlierPar1 =      centOADB->V0MZDCOutlierPar1();      
  fV0MZDCEcalOutlierPar0 =  centOADB->V0MZDCEcalOutlierPar0();  
  fV0MZDCEcalOutlierPar1 =  centOADB->V0MZDCEcalOutlierPar1();  

  return 0;
}



//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MSPD(Float_t spd, Float_t v0, Int_t cent) const
{
  // Clean outliers
  Float_t val = fV0MSPDOutlierPar0 +  fV0MSPDOutlierPar1 * v0;
  Float_t spdSigma = fV0MSPDSigmaOutlierPar0 + fV0MSPDSigmaOutlierPar1*cent + fV0MSPDSigmaOutlierPar2*cent*cent;
  if ( TMath::Abs(spd-val) > fOutliersCut*spdSigma ) 
    return kTRUE;
  else 
    return kFALSE;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MTPC(Int_t tracks, Float_t v0, Int_t cent) const
{
  // Clean outliers
  Float_t val = fV0MTPCOutlierPar0 +  fV0MTPCOutlierPar1 * v0;
  Float_t tpcSigma = fV0MTPCSigmaOutlierPar0 + fV0MTPCSigmaOutlierPar1*cent + fV0MTPCSigmaOutlierPar2*cent*cent;
  if ( TMath::Abs(tracks-val) > fOutliersCut*tpcSigma ) 
    return kTRUE;
  else 
    return kFALSE;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MZDC(Float_t zdc, Float_t v0) const
{
  // Clean outliers
  Float_t val = fV0MZDCOutlierPar0 + fV0MZDCOutlierPar1 * v0;
  if (zdc >  val) 
    return kTRUE;
  else 
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCentralitySelectionTask::IsOutlierV0MZDCECal(Float_t zdc, Float_t v0) const
{
  // Clean outliers
  Float_t val = fV0MZDCEcalOutlierPar0 + fV0MZDCEcalOutlierPar1 * v0;
  if (zdc >  val) 
    return kTRUE;
  else 
    return kFALSE;
}
