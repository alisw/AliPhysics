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
#include "AliESDVertex.h"
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
  fRunNo(-1),
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
  fRunNo(-1),
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
  fRunNo(ana.fRunNo),
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
  Float_t  spdCorr =0;              //  corrected spd2 multiplicity

  Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
  Float_t  multFMDA = 0;            //  multiplicity from FMD on detector A
  Float_t  multFMDC = 0;            //  multiplicity from FMD on detector C

  Short_t v0Corr = 0;               // corrected V0 multiplicity
  Short_t v0CorrResc = 0;           // corrected and rescaled V0 multiplicity

  Float_t zvtx =0;                  // z-vertex SPD
 
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

    float v0CorrR;
    v0Corr = (Short_t) (GetCorrV0(esd,v0CorrR));
    v0CorrResc = (Short_t)v0CorrR;

    // ***** Vertex Info
    const AliESDVertex* vtxESD = esd->GetPrimaryVertexSPD();
    zvtx        = vtxESD->GetZ(); 

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

    spdCorr = GetCorrSPD2(nClusters[1],zvtx);    

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
  if(fHtempV0M)  fCentV0M = fHtempV0M->GetBinContent(fHtempV0M->FindBin((v0Corr)));
  ///  else     printf("  Centrality by V0 not available!!!\n\n");
  if(fHtempFMD) fCentFMD = fHtempFMD->GetBinContent(fHtempFMD->FindBin((multFMDA+multFMDC)));
  //  else     printf("  Centrality by FMD not available!!!\n\n");
  if(fHtempTRK) fCentTRK = fHtempTRK->GetBinContent(fHtempTRK->FindBin(nTracks));
  //  else     printf("  Centrality by TRK not available!!!\n\n");
  if(fHtempTKL) fCentTKL = fHtempTKL->GetBinContent(fHtempTKL->FindBin(nTracklets));
  //  else     printf("  Centrality by TKL not available!!!\n\n");
  if(fHtempCL0) fCentCL0 = fHtempCL0->GetBinContent(fHtempCL0->FindBin(nClusters[0]));
  //  else     printf("  Centrality by CL0 not available!!!\n\n");
  if(fHtempCL1) fCentCL1 = fHtempCL1->GetBinContent(fHtempCL1->FindBin(spdCorr));
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

  fRunNo = fCurrentRun;

  // CHANGE HERE FOR RUN RANGES
  if ( fRunNo == 137162 ) fRunNo = 137161;
  else if ( fRunNo == 137365 ) fRunNo = 137366;
  else if ( fRunNo > 137366 ) fRunNo = 137366;
  // CHANGE HERE FOR RUN RANGES

  TString runName(Form("%d", fRunNo));
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
//________________________________________________________________________
Float_t AliCentralitySelectionTask::GetCorrV0(const AliESDEvent* esd, float &v0CorrResc) 
{
  // correct V0 non-linearity, prepare a version rescaled to SPD2 corr
  Double_t *par0;
  Double_t *par1;
  Double_t *par2;
  
  Double_t par0_137161[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
			5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
			5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
			6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
			5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
			6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
			2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
			2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
			4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
			4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
			4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  Double_t par1_137161[64] = { -6.68e-05 , -7.78e-05 , -6.88e-05 , -5.92e-05 , 
			-2.43e-05 , -3.54e-05 , -2.91e-05 , -1.99e-05 , -1.40e-05 , -4.01e-05 , 
			-2.29e-05 , -3.68e-05 , -2.53e-05 , -2.44e-06 , -9.22e-06 , -1.51e-05 , 
			-2.80e-05 , -2.34e-05 , -1.72e-05 , -1.81e-05 , -1.29e-05 , -2.65e-05 , 
			-1.61e-05 , -2.86e-05 , -1.74e-05 , -4.23e-05 , -3.41e-05 , -1.05e-04 , 
			-2.76e-05 , -4.71e-05 , -3.06e-05 , -2.32e-05 , -1.55e-06 , 2.15e-05 , 
			1.40e-05 , 2.16e-05 , 1.21e-05 , 3.05e-06 , 1.67e-05 , -3.84e-06 , 
			3.09e-06 , 1.50e-05 , 3.47e-06 , 4.87e-06 , -3.71e-07 , -1.75e-06 , 
			-1.80e-06 , 9.99e-06 , -6.46e-06 , -4.91e-06 , 1.33e-05 , -2.52e-07 , 
			-3.85e-06 , 4.94e-06 , -2.48e-07 , -1.20e-05 , 2.07e-06 , 6.12e-06 , 
			-1.18e-06 , 4.54e-06 , -1.54e-05 , -1.25e-05 , 1.46e-06 , -6.67e-06 };
  Double_t par2_137161[64] = { 1.29e-08 , 1.51e-08 , 1.43e-08 , 1.11e-08 , 
			5.04e-09 , 6.99e-09 , 5.58e-09 , 4.15e-09 , 4.00e-09 , 8.22e-09 , 
			4.97e-09 , 7.66e-09 , 4.91e-09 , 1.10e-09 , 2.64e-09 , 3.64e-09 , 
			5.76e-09 , 5.46e-09 , 3.38e-09 , 3.47e-09 , 2.43e-09 , 4.13e-09 , 
			2.80e-09 , 5.80e-09 , 3.86e-09 , 7.46e-09 , 5.98e-09 , 2.58e-08 , 
			5.50e-09 , 8.72e-09 , 5.23e-09 , 4.37e-09 , 2.33e-09 , -6.01e-10 , 
			3.99e-11 , -2.02e-10 , 7.67e-10 , 2.03e-09 , 1.17e-10 , 2.56e-09 , 
			1.16e-09 , -4.75e-10 , 1.28e-09 , 1.23e-09 , 1.62e-09 , 1.61e-09 , 
			1.93e-09 , 2.97e-10 , 2.21e-09 , 2.16e-09 , 5.22e-10 , 1.03e-09 , 
			1.56e-09 , 5.00e-10 , 1.01e-09 , 2.93e-09 , 1.05e-09 , 9.96e-11 , 
			1.21e-09 , 7.45e-10 , 3.07e-09 , 2.31e-09 , 6.70e-10 , 1.89e-09 };

  Double_t par0_137366[64] = { 7.12e-02 , 7.34e-02 , 7.39e-02 , 6.54e-02 , 6.11e-02 , 6.31e-02 , 6.15e-02 , 
			       6.00e-02 , 6.10e-02 , 6.49e-02 , 6.17e-02 , 6.33e-02 , 6.00e-02 , 5.48e-02 , 
			       5.44e-02 , 5.81e-02 , 6.49e-02 , 7.07e-02 , 5.91e-02 , 6.18e-02 , 4.82e-02 , 
			       5.67e-02 , 5.36e-02 , 6.60e-02 , 6.37e-02 , 6.78e-02 , 6.31e-02 , 1.04e-01 , 
			       6.91e-02 , 7.32e-02 , 6.61e-02 , 6.16e-02 , 2.64e-02 , 2.81e-02 , 2.64e-02 , 
			       2.85e-02 , 2.87e-02 , 2.18e-02 , 2.19e-02 , 2.43e-02 , 2.81e-02 , 4.37e-02 , 
			       3.90e-02 , 4.66e-02 , 4.24e-02 , 4.09e-02 , 4.21e-02 , 3.88e-02 , 4.83e-02 , 
			       5.23e-02 , 5.44e-02 , 4.85e-02 , 4.42e-02 , 4.58e-02 , 4.74e-02 , 3.14e-02 , 
			       6.31e-02 , 5.30e-02 , 5.01e-02 , 5.33e-02 , 5.70e-02 , 3.95e-02 , 4.98e-02 , 5.31e-02 };
  Double_t par1_137366[64] = { -6.99e-05 , -6.99e-05 , -6.94e-05 , -6.55e-05 , -3.55e-05 , -4.50e-05 , 
			       -3.10e-05 , -2.81e-05 , -2.29e-05 , -3.89e-05 , -2.53e-05 , -4.25e-05 ,
			       -1.87e-05 , -2.01e-05 , -1.53e-05 , -2.14e-05 , -2.86e-05 , -4.70e-05 ,
			       -2.23e-05 , -3.30e-05 ,-9.74e-06 , -2.62e-05 , -1.76e-05 , -2.38e-05 , 
			       -2.40e-05 , -3.43e-05 , -2.75e-05 , -6.86e-05 ,-2.35e-05 , -4.45e-05 , 
			       -2.51e-05 , -2.20e-05 , -1.25e-16 , -2.04e-17 , -2.06e-17 , -3.74e-19 ,
			       -1.18e-18 , -2.02e-15 , -3.78e-06 , -1.26e-06 , -2.71e-06 , -6.23e-17 , 
			       -7.39e-08 , -1.76e-16 , -8.98e-06 , -4.10e-18 , -1.34e-05 , -1.06e-16 , 
			       -3.34e-06 , -1.04e-05 , -5.28e-06 , -7.34e-06 , -1.05e-05 , -7.68e-06 ,
			       -1.78e-05 , -1.19e-05 , -1.78e-05 , -1.34e-06 , -9.23e-06 , -3.34e-06 ,
			       -8.02e-06 , -1.39e-05 , -1.38e-05 , -1.40e-05 };
  Double_t par2_137366[64] = { 1.41e-08 , 1.47e-08 , 1.48e-08 , 1.24e-08 , 6.82e-09 , 8.73e-09 , 6.26e-09 , 
			       5.53e-09 , 5.40e-09 , 7.93e-09 , 5.49e-09 , 8.77e-09 , 4.21e-09 , 3.93e-09 , 
			       3.60e-09 , 4.67e-09 , 5.59e-09 , 8.81e-09 , 3.89e-09 , 6.19e-09 , 1.97e-09 , 
			       4.38e-09 , 3.26e-09 , 5.00e-09 , 4.58e-09 , 6.39e-09 , 5.03e-09 , 1.30e-08 , 
			       4.95e-09 , 8.26e-09 , 4.57e-09 , 4.10e-09 , 2.35e-09 , 2.30e-09 , 2.15e-09 , 
			       2.27e-09 , 2.17e-09 , 2.27e-09 , 2.97e-09 , 2.25e-09 , 1.69e-09 , 1.44e-09 , 
			       1.66e-09 , 1.75e-09 , 2.88e-09 , 1.82e-09 , 3.64e-09 , 1.80e-09 , 1.71e-09 , 
			       2.66e-09 , 3.01e-09 , 1.95e-09 , 2.64e-09 , 2.42e-09 , 3.68e-09 , 2.66e-09 , 
			       3.92e-09 , 1.18e-09 , 2.26e-09 , 1.57e-09 , 2.02e-09 , 2.71e-09 , 2.99e-09 , 3.04e-09 }; 
  
  
  if (esd->GetRunNumber() <= 137165) {
    par0=par0_137161;
    par1=par1_137161;
    par2=par2_137161;
  }  else  {
    par0=par0_137366;
    par1=par1_137366;
    par2=par2_137366;
 }
  //
  Float_t multCorr = 0;
  Float_t multCorr2 = 0;
  Float_t multChCorr[64];
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  for(Int_t i = 0; i < 64; ++i) {
    Double_t b = (esdV0->GetMultiplicity(i)*par1[i]-par0[i]);
    Double_t s = (b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i));
    Double_t n = (s<0) ? -b : (-b + TMath::Sqrt(s));
    multChCorr[i] = 2.*esdV0->GetMultiplicity(i)/n*par0[i];
    multCorr += multChCorr[i];
    multCorr2 += (multChCorr[i]/par0[i]/64.);
  }
  v0CorrResc =  multCorr2;
  return multCorr;
}

//____________________________________________________________________
Float_t AliCentralitySelectionTask::GetCorrSPD2(Float_t spd2raw,Float_t zv)
{
  // renormalize N spd2 clusters at given Zv to acceptance at Zv=0
  const double pars[] = {8.10030e-01,-2.80364e-03,-7.19504e-04};
  zv -= pars[0];
  float corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? spd2raw/corr : -1;
}
