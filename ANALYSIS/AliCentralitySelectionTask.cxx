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
#include <TObjString.h>
#include <TString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TDirectory.h>
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
  fCentfilename(""),
  fCentfilename2(""),
  fFileList(new TList),
  fFileList2(new TList),
  fCurrentRun(-1),
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
  fHtempZEMvsZDC(0)
{   
  // Default constructor
  fFileList->SetOwner();
  fFileList2->SetOwner();

  AliInfo("Centrality Selection enabled.");
}   

//________________________________________________________________________
AliCentralitySelectionTask::AliCentralitySelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fDebug(0),
  fAnalysisInput("ESD"),
  fIsMCInput(kFALSE),
  fFile(0),
  fFile2(0),
  fCentfilename(""),
  fCentfilename2(""),
  fFileList(new TList),
  fFileList2(new TList),
  fCurrentRun(-1),
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
  fHtempZEMvsZDC(0)
{
  // Default constructor
  fFileList->SetOwner();
  fFileList2->SetOwner();

  AliInfo("Centrality Selection enabled.");
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
  fFileList(ana.fFileList),
  fFileList2(ana.fFileList2),
  fCurrentRun(ana.fCurrentRun),
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
  fHtempZEMvsZDC(ana.fHtempZEMvsZDC)
{
  // Copy Constructor	
}
 
//________________________________________________________________________
AliCentralitySelectionTask::~AliCentralitySelectionTask()
{
  // Destructor
  
  if (fFileList) {
    fFileList->Clear();
    delete fFileList;
  }
  fFileList = NULL;

  if (fFileList2) {
    fFileList2->Clear();
    delete fFileList2;
  }
  fFileList2 = NULL;
}  

//________________________________________________________________________
void AliCentralitySelectionTask::UserCreateOutputObjects()
{  
  // Create the output containers
  if(fDebug>1) printf("AnalysisCentralitySelectionTask::UserCreateOutputObjects() \n");

  if (fFileList->GetEntries() < 1) {
    AliError("No Inputfiles Added");
  }

  AliLog::SetClassDebugLevel("AliCentralitySelectionTask", AliLog::kInfo);
}

//________________________________________________________________________
void AliCentralitySelectionTask::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliCentralitySelectionTask::UserExec() \n");
  
  Float_t  zncEnergy;          //  ZNC Energy
  Float_t  zpcEnergy;          //  ZPC Energy
  Float_t  znaEnergy;          //  ZNA Energy
  Float_t  zpaEnergy;          //  ZPA Energy
  Float_t  zem1Energy = 0.;         //  ZEM1 Energy
  Float_t  zem2Energy = 0.;         //  ZEM2 Energy
  
  Int_t    nTracks = 0;             //  no. tracks
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

    if (SetupRun(esd))   
      return;

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
  if(fHtempV0M)  fCentV0M = fHtempV0M->GetBinContent(fHtempV0M->FindBin((multV0A+multV0C)));
  ///  else     printf("  Centrality by V0 not available!!!\n\n");
  if(fHtempFMD) fCentFMD = fHtempFMD->GetBinContent(fHtempFMD->FindBin((multFMDA+multFMDC)));
  //  else     printf("  Centrality by FMD not available!!!\n\n");
  if(fHtempTRK) fCentTRK = fHtempTRK->GetBinContent(fHtempTRK->FindBin(nTracks));
  //  else     printf("  Centrality by TRK not available!!!\n\n");
  if(fHtempTKL) fCentTKL = fHtempTKL->GetBinContent(fHtempTKL->FindBin(nTracklets));
  //  else     printf("  Centrality by TKL not available!!!\n\n");
  if(fHtempCL0) fCentCL0 = fHtempCL0->GetBinContent(fHtempCL0->FindBin(nClusters[0]));
  //  else     printf("  Centrality by CL0 not available!!!\n\n");
  if(fHtempCL1) fCentCL1 = fHtempCL1->GetBinContent(fHtempCL1->FindBin(nClusters[1]));
  ///  else     printf("  Centrality by CL1 not available!!!\n\n");
  
  if(fHtempV0MvsFMD) fCentV0MvsFMD = fHtempV0MvsFMD->GetBinContent(fHtempV0MvsFMD->FindBin((multV0A+multV0C)));
  //  else     printf("  Centrality by V0 vs FMD not available!!!\n\n");
  if(fHtempTKLvsV0M) fCentTKLvsV0M = fHtempTKLvsV0M->GetBinContent(fHtempTKLvsV0M->FindBin(nTracklets));
  //  else     printf("  Centrality by V0 vs TKL not available!!!\n\n");
  if(fHtempZEMvsZDC) fCentZEMvsZDC = fHtempZEMvsZDC->GetBinContent(fHtempZEMvsZDC->FindBin((zem1Energy+zem2Energy)/1000.));
  //  else     printf("  Centrality by ZEM vs ZDC not available!!!\n\n");

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

//________________________________________________________________________
void AliCentralitySelectionTask::ReadCentralityHistos() 
{
//  Read centrality histograms
    TDirectory *owd = gDirectory;
    fFile  = TFile::Open(fCentfilename);
    owd->cd();
    fHtempV0M  = (TH1F*) (fFile->Get("hmultV0_percentile"));
    fHtempFMD  = (TH1F*) (fFile->Get("hmultFMD_percentile"));
    fHtempTRK  = (TH1F*) (fFile->Get("hNtracks_percentile"));
    fHtempTKL  = (TH1F*) (fFile->Get("hNtracklets_percentile"));
    fHtempCL0  = (TH1F*) (fFile->Get("hNclusters0_percentile"));
    fHtempCL1  = (TH1F*) (fFile->Get("hNclusters1_percentile"));
}  

//________________________________________________________________________
void AliCentralitySelectionTask::ReadCentralityHistos2() 
{
//  Read centrality histograms
    TDirectory *owd = gDirectory;
    fFile2  = TFile::Open(fCentfilename2);
    owd->cd();
    fHtempV0MvsFMD =  (TH1F*) (fFile2->Get("hmultV0vsmultFMD_all_percentile"));
    fHtempTKLvsV0M  = (TH1F*) (fFile2->Get("hNtrackletsvsmultV0_all_percentile"));
    fHtempZEMvsZDC  = (TH1F*) (fFile2->Get("hEzemvsEzdc_all_percentile"));
}

//________________________________________________________________________
void AliCentralitySelectionTask::SetPercentileFile(TString filename) 
{
// Set the percentile file name
  fCentfilename = filename;
  ReadCentralityHistos();
}

//________________________________________________________________________
void AliCentralitySelectionTask::SetPercentileFile2(TString filename) 
{
// Set the percentile file name
  fCentfilename2 = filename;
  ReadCentralityHistos2();
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
Int_t AliCentralitySelectionTask::SetupRun(AliESDEvent* esd)
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

  Int_t runNo = fCurrentRun;

  // CHANGE HERE FOR RUN RANGES
  if ( runNo == 137162 ) runNo = 137161;
  else if ( runNo == 137365 ) runNo = 137366;
  // CHANGE HERE FOR RUN RANGES

  TString runName(Form("%d", runNo));
  TString fileName("");
  Bool_t isRunKnown = kFALSE;

  // Check if run is in fileList
  // if not, take the last name in the list
  for ( Int_t idx=0 ; idx < fFileList->GetEntries(); ++idx ) {

    TString str((dynamic_cast<TObjString*>(fFileList->At(idx)))->GetString());
    if (str.Contains(runName)) {
      fileName += str;
      isRunKnown = kTRUE;
      break;
    }
  }

  if (!isRunKnown) {
    if (fFileList->Last()) {
      fileName += (dynamic_cast<TObjString*>(fFileList->Last()))->GetString();
      AliError(Form("Run %d not known to centrality selection!", fCurrentRun));
    }
  }

  if (fileName.Contains(".root")) {
    AliInfo(Form("Centrality Selection for run %d is initialized with %s", fCurrentRun, fileName.Data()));
    SetPercentileFile(fileName.Data());
    return 0;
  }
  
  return -1;
}
