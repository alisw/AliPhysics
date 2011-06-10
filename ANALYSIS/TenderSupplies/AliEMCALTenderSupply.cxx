/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  EMCAL tender, apply corrections to EMCAl clusters                        //
//  and do track matching                                                    //                                                                           
//  Author : Deepa Thomas (Utrecht University)                               //                      
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"
//#include "TAlienFile.h"
#include "TGrid.h"
#include "TTree.h"
#include "TInterpreter.h"
#include "TObjArray.h"

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include "AliOADBContainer.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliESDCaloCluster.h"
#include "AliEMCALTenderSupply.h"

//EMCAL
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

AliEMCALTenderSupply::AliEMCALTenderSupply() :
	AliTenderSupply()
	,fEMCALGeo(0x0)
	,fEMCALGeoName("EMCAL_FIRSTYEARV1")
	,fEMCALRecoUtils(new AliEMCALRecoUtils)
	,fConfigName("")
	,fDebugLevel(0)
	,fNonLinearFunc(AliEMCALRecoUtils::kNoCorrection) 
	,fNonLinearThreshold(30)      	
	,fReCalibCluster(kFALSE)	
	,fReCalibCell(kFALSE)	
	,fRecalClusPos(kFALSE)
	,fFiducial(kFALSE) 
	,fNCellsFromEMCALBorder(1)	
	,fRecalDistToBadChannels(kFALSE)	
	,fInputTree(0)	
	,fInputFile(0)
	,fFilepass(0) 
	,fMass(0.139)
	,fStep(1)
	,fRcut(0.05)	
{
	//
	// default ctor
	//
	for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
}

//_____________________________________________________
AliEMCALTenderSupply::AliEMCALTenderSupply(const char *name, const AliTender *tender) :
	AliTenderSupply(name,tender)
	,fEMCALGeo(0x0)
	,fEMCALGeoName("EMCAL_FIRSTYEARV1")
	,fEMCALRecoUtils(new AliEMCALRecoUtils)
	,fConfigName("") 
	,fDebugLevel(0)
	,fNonLinearFunc(AliEMCALRecoUtils::kNoCorrection)      	
	,fNonLinearThreshold(30)      	
	,fReCalibCluster(kFALSE)	
	,fReCalibCell(kFALSE)	
	,fRecalClusPos(kFALSE)
	,fFiducial(kFALSE) 
	,fNCellsFromEMCALBorder(1)	
	,fRecalDistToBadChannels(kFALSE)	
	,fInputTree(0)	
	,fInputFile(0)
	,fFilepass(0)
	,fMass(0.139)
	,fStep(1)
	,fRcut(0.05)
{
	//
	// named ctor
	//
	for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
}

//_____________________________________________________
AliEMCALTenderSupply::~AliEMCALTenderSupply()
{
	//Destructor
	delete fEMCALGeo;
	delete fEMCALRecoUtils;
	delete fInputTree;
	delete fInputFile;
}

//_____________________________________________________
void AliEMCALTenderSupply::Init()
{
	//
	// Initialise EMCAL tender
	//

	if (fDebugLevel>0) AliInfo("Init EMCAL Tender supply\n");	

	AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();

	fInputTree = mgr->GetTree();

	if(gROOT->LoadMacro(fConfigName) >=0){
		AliEMCALTenderSupply *tender = (AliEMCALTenderSupply*)gInterpreter->ProcessLine("ConfigEMCALTenderSupply()");
		fDebugLevel         = tender->fDebugLevel;
		fEMCALGeoName       = tender->fEMCALGeoName; 
		fEMCALRecoUtils     = tender->fEMCALRecoUtils; 
		fConfigName         = tender->fConfigName;
		fNonLinearFunc      = tender->fNonLinearFunc;
		fNonLinearThreshold = tender->fNonLinearThreshold;
		fReCalibCluster     = tender->fReCalibCluster;
		fReCalibCell        = tender->fReCalibCell;
		fRecalClusPos       = tender->fRecalClusPos;
		fFiducial	    = tender->fFiducial;
		fNCellsFromEMCALBorder = tender->fNCellsFromEMCALBorder;
		fRecalDistToBadChannels = tender->fRecalDistToBadChannels;    
		fMass               = tender->fMass;
		fStep               = tender->fStep;
		fRcut               = tender->fRcut;
	}

	// Init goemetry	
	fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;

	fEMCALRecoUtils = new AliEMCALRecoUtils();

	//Initialising Non linearity parameters
	fEMCALRecoUtils->SetNonLinearityThreshold(fNonLinearThreshold);

	//Setting mass, step size and residual cut 
	fEMCALRecoUtils->SwitchOnCutEtaPhiSum(); 
	fEMCALRecoUtils->SetCutR(fRcut);
	fEMCALRecoUtils->SetMass(fMass);
	fEMCALRecoUtils->SetStep(fStep);

	if(fDebugLevel>1) fEMCALRecoUtils->Print("");

}

//_____________________________________________________
void AliEMCALTenderSupply::ProcessEvent()
{
	//Event loop	
	AliESDEvent *event=fTender->GetEvent();
	if (!event) return;

	if(fTender->RunChanged()){ 
		//Initialising parameters once per run number
		if (!InitBadChannels()) return;
		if (fRecalClusPos){ if (!InitMisalignMatrix()) return;}
		if (fReCalibCluster || fReCalibCell){ if (!InitRecalib()) return;}
	}

	AliESDCaloCells *cells= event->GetEMCALCells();

	//------------Good clusters------------
	TClonesArray *clusArr;

	clusArr = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
	if(!clusArr) clusArr = dynamic_cast<TClonesArray*>(event->FindListObject("CaloClusters"));
	if(!clusArr) return;

	for (Int_t icluster=0; icluster < clusArr->GetEntries(); icluster++ ){
		AliVCluster *clust = static_cast<AliVCluster*>(clusArr->At(icluster));
		if(!clust) continue;
		if (!clust->IsEMCAL()) continue;
		if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo, clust->GetCellsAbsId(), clust->GetNCells() )){
			delete clusArr->RemoveAt(icluster);
			continue;
		}
		if(fFiducial){
			if(!fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo, clust, cells)){
				delete clusArr->RemoveAt(icluster);
				continue;
			}
		}
		fEMCALRecoUtils->CorrectClusterEnergyLinearity(clust);
		if(fRecalDistToBadChannels) fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, cells, clust);  
		if(fReCalibCluster) fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, clust, cells);
		if(fRecalClusPos) fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, cells, clust);
		//fEMCALRecoUtils->SetTimeDependentCorrections(event->GetRunNumber());
	}
	clusArr->Compress();

	//------------ EMCAL cells loop ------------
	
	//Recalibrate cells
	if(fReCalibCell) RecalibrateCells();
	
	//-------- Track Matching ------------------

	//set magnetic field
	Double_t  magF = event->GetMagneticField();
	Double_t magSign = 1.0;
	if(magF<0)magSign = -1.0;

	if (!TGeoGlobalMagField::Instance()->GetField()) 
	{
		AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
		TGeoGlobalMagField::Instance()->SetField(field);
	}


	fEMCALRecoUtils->FindMatches(event,0x0,fEMCALGeo);

	Int_t nTracks = event->GetNumberOfTracks();

	if(nTracks>0){
		SetClusterMatchedToTrack(event);
		SetTracksMatchedToCluster(event);
	}
	
}

//_____________________________________________________
void AliEMCALTenderSupply::SetClusterMatchedToTrack(AliESDEvent *event)
{
	//checks if tracks are matched to EMC clusters and set the matched EMCAL cluster index to ESD track. 

	Int_t matchClusIndex = -1;
	Int_t nTracks = event->GetNumberOfTracks();

	// Track loop 
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) 
	{
		AliESDtrack* track = event->GetTrack(iTrack);
		if (!track) {
			printf("ERROR: Could not receive track %d\n", iTrack);
			continue;
		}

		matchClusIndex = fEMCALRecoUtils->GetMatchedClusterIndex(iTrack);		   
		track->SetEMCALcluster(matchClusIndex); //sets -1 if track not matched within residual
	}

	if (fDebugLevel>2) AliInfo("Track matched to closest cluster\n");	
}

//_____________________________________________________
void AliEMCALTenderSupply::SetTracksMatchedToCluster(AliESDEvent *event)
{
	//checks if EMC clusters are matched to ESD track.
	//Adds track indexes of all the tracks matched to a cluster withing residuals in ESDCalocluster

	Int_t matchTrackIndex = -1;
	Int_t nTracks = event->GetNumberOfTracks();

	// Cluster loop
	for (Int_t iClus=0; iClus < event->GetNumberOfCaloClusters(); iClus++)
	{
		AliESDCaloCluster * cluster = event->GetCaloCluster(iClus);
		if (!cluster->IsEMCAL()) continue;

		Int_t nMatched =0;
		TArrayI arrayTrackMatched(nTracks);

		//get the closest track matched to the cluster
		matchTrackIndex = fEMCALRecoUtils->GetMatchedTrackIndex(iClus); 
		arrayTrackMatched[nMatched] = matchTrackIndex;
		nMatched++;

		//get all other tracks matched to the cluster
		//track loop
		for(Int_t iTrk=0; iTrk<nTracks; iTrk++){
			AliESDtrack* track = event->GetTrack(iTrk);
			if(iTrk == matchTrackIndex) continue;

			if(track->GetEMCALcluster() == iClus){
				arrayTrackMatched[nMatched] = iTrk;
				nMatched++;
			}
		}

		arrayTrackMatched.Set(nMatched);
		cluster->AddTracksMatched(arrayTrackMatched);

		Float_t eta= -999, phi = -999;
		if(matchTrackIndex != -1) fEMCALRecoUtils->GetMatchedResiduals(iClus, eta, phi);
		cluster->SetTrackDistance(phi, eta);
	}
	if (fDebugLevel>2) AliInfo("Cluster matched to tracks \n");	
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitMisalignMatrix()
{
	//Initialising Misalignment matrix

	AliESDEvent *event=fTender->GetEvent();
	if (!event) return kFALSE;
	
	if (fDebugLevel>0) AliInfo("Initialising Misalignment matrix \n");	

	if(fInputTree){ 
		fInputFile = fInputTree->GetCurrentFile();
		if(fInputFile){
			const char *fileName = fInputFile->GetName();
			TString FileName = TString(fileName);
			if     (FileName.Contains("pass1")) fFilepass = TString("pass1");
			else if(FileName.Contains("pass2")) fFilepass = TString("pass2");
			else if(FileName.Contains("pass3")) fFilepass = TString("pass3");
			else AliError("pass number not found");
		}
		else AliError("File not found");
	}
	else AliError("Tree not found");

	Int_t runGM = 0;
	runGM = event->GetRunNumber();

	fEMCALRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

	if(runGM <=140000){ //2010 data

		Double_t rotationMatrix[4][9] = {{-0.014587, -0.999892, -0.002031, 0.999892, -0.014591,  0.001979, -0.002009, -0.002002,  0.999996},
			{-0.014587,  0.999892,  0.002031, 0.999892,  0.014591, -0.001979, -0.002009,  0.002002, -0.999996},
			{-0.345864, -0.938278, -0.003412, 0.938276, -0.345874,  0.003010, -0.004004, -0.002161,  0.999990},
			{-0.345861,  0.938280,  0.003412, 0.938276,  0.345874, -0.003010, -0.004004,  0.002161, -0.999990}};

		Double_t translationMatrix[4][3] = {{0.351659,    447.576446,  176.269742},
			{1.062577,    446.893974, -173.728870},
			{-154.213287, 419.306156,  176.753692},
			{-153.018950, 418.623681, -173.243605}};

		for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
		{
			//if(DebugLevel() > 1)  fEMCALMatrix[mod]->Print();
			fEMCALMatrix[mod]->SetRotation(rotationMatrix[mod]);
			fEMCALMatrix[mod]->SetTranslation(translationMatrix[mod]);		
			fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod); 
		}
	}


	if(runGM>140000 && runGM <148531 && (fFilepass = "pass1"))
	{ // 2011 LHC11a pass1 data
		AliOADBContainer emcalgeoCont(Form("emcal2011"));
		emcalgeoCont.InitFromFile(Form("BetaGood.root"),Form("AliEMCALgeo"));

		TObjArray *mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byS");

		for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
			fEMCALMatrix[mod] = (TGeoHMatrix*) mobj->At(mod);
			fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod); 

			fEMCALMatrix[mod]->Print();
		}
	}

	else AliInfo("MISALLIGNMENT NOT APPLIED");

	return kTRUE;
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitBadChannels()
{
	//Initialising Bad channel maps

	AliESDEvent *event=fTender->GetEvent();
	if (!event) return kFALSE;

	Int_t fRunBC= 0;
	fRunBC = event->GetRunNumber();

	if (fDebugLevel>0) AliInfo("Initialising Bad channel map \n");	

	if(fFiducial){
		fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(fNCellsFromEMCALBorder);
		fEMCALRecoUtils->SwitchOnNoFiducialBorderInEMCALEta0();
	}

	fEMCALRecoUtils->SwitchOnBadChannelsRemoval();
	if(fRecalDistToBadChannels) fEMCALRecoUtils->SwitchOnDistToBadChannelRecalculation();

	TFile *fbad;
	//2010 
	if(fRunBC <=140000){
		fbad = new TFile("BadChannels.root","read");
		if (fbad->IsZombie()){
			if (fDebugLevel>1) AliInfo("Connecting to alien to get BadChannels.root \n");	
			TGrid::Connect("alien://");
			fbad = TFile::Open("alien:///alice/cern.ch/user/g/gconesab/BadChannelsDB/BadChannels.root");
		}

		TH2I * hb0 = ( TH2I *)fbad->Get("EMCALBadChannelMap_Mod0");
		TH2I * hb1 = ( TH2I *)fbad->Get("EMCALBadChannelMap_Mod1");
		TH2I * hb2 = ( TH2I *)fbad->Get("EMCALBadChannelMap_Mod2");
		TH2I * hb3 = ( TH2I *)fbad->Get("EMCALBadChannelMap_Mod3"); 
		fEMCALRecoUtils->SetEMCALChannelStatusMap(0,hb0);
		fEMCALRecoUtils->SetEMCALChannelStatusMap(1,hb1);
		fEMCALRecoUtils->SetEMCALChannelStatusMap(2,hb2);
		fEMCALRecoUtils->SetEMCALChannelStatusMap(3,hb3); 
	}

	//2011
	if(fRunBC>140000){

		const Int_t nTowers=89;
		Int_t hotChannels[nTowers]={74, 103, 152, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 368, 369, 370, 371, 372, 373, 374,375, 376, 377, 378, 379, 380, 381, 382, 383, 917, 1275, 1288, 1519, 1595, 1860, 1967, 2022, 2026, 2047, 2117, 2298, 2540, 2776, 3135, 3764, 6095, 6111, 6481, 6592, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7371, 7425, 7430, 7457, 7491, 7709, 8352, 8353, 8356, 8357, 8808, 8810, 8812, 8814, 9056, 9769, 9815, 9837};

		Int_t nSupMod=-1, nModule=-1, nIphi=-1, nIeta=-1, iphi=-1, ieta=-1;
		for(Int_t i=0; i<nTowers; i++)
		{
			fEMCALGeo->GetCellIndex(hotChannels[i],nSupMod,nModule,nIphi,nIeta);

			fEMCALGeo->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta,iphi,ieta);
			fEMCALRecoUtils->SetEMCALChannelStatus(nSupMod, ieta, iphi);
		}

	}

	return kTRUE;
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitRecalib()
{
	//Initialising Recalibration Factors

	AliESDEvent *event=fTender->GetEvent();
	if (!event) return kFALSE;

	if (fDebugLevel>0) AliInfo("Initialising Recalibration factors \n");
	
	if(fInputTree){ 
		fInputFile = fInputTree->GetCurrentFile();
		if(fInputFile){
			const char *fileName = fInputFile->GetName();
			TString FileName = TString(fileName);
			if     (FileName.Contains("pass1")) fFilepass = TString("pass1");
			else if(FileName.Contains("pass2")) fFilepass = TString("pass2");
			else if(FileName.Contains("pass3")) fFilepass = TString("pass3");
			else AliError("pass number not found");
		}
		else AliError("File not found");
	}
	else AliError("Tree not found");
	//        else {cout << "Tree not found " <<endl; return kFALSE;}
  
	Int_t runRC = event->GetRunNumber();

  //if (event->GetRunNumber()==runRC)
  //     return kFALSE;

	fEMCALRecoUtils->SwitchOnRecalibration();

	TFile* fRecalib;

	if(runRC <=140000){
		fRecalib = new TFile("RecalibrationFactors.root","read");
		if (fRecalib->IsZombie()){
			if (fDebugLevel>1) AliInfo("Connecting to alien to get RecalibrationFactors.root \n");	

			if(     (runRC >= 114737 && runRC <= 117223 && (fFilepass = "pass2")) || //LHC10b pass2 
					(runRC >= 118503 && runRC <= 121040 && ((fFilepass = "pass2")||(fFilepass = "pass3"))) || //LHC10c pass2, LHC10c pass3
					(runRC >= 122195 && runRC <= 126437 && (fFilepass = "pass1")) || //LHC10d pass1
					(runRC >= 127712 && runRC <= 130850 && (fFilepass = "pass1")) || //LHC10e pass1
					(runRC >= 133004 && runRC < 134657  && (fFilepass = "pass1")))// LHC10f pass1 <134657
			{
				TGrid::Connect("alien://");
				fRecalib = TFile::Open("alien:///alice/cern.ch/user/g/gconesab/RecalDB/summer_december_2010/RecalibrationFactors.root");
			}

			else if((runRC >= 122195 && runRC <= 126437 && (fFilepass = "pass2")) || //LHC10d pass2
					(runRC >= 134657 && runRC <= 135029 && (fFilepass = "pass1")) || //LHC10f pass1 >= 134657
					(runRC >= 135654 && runRC <= 136377 && (fFilepass = "pass1")) || //LHC10g pass1
					(runRC >= 136851 && runRC < 137231  && (fFilepass = "pass1"))) //LHC10h pass1 untill christmas
			{
				TGrid::Connect("alien://");
				fRecalib = TFile::Open("alien:///alice/cern.ch/user/g/gconesab/RecalDB/december2010/RecalibrationFactors.root");
			}

			else if(runRC >= 137231 && runRC <= 139517 && (fFilepass = "pass1")) //LHC10h pass1 from christmas
			{
				TGrid::Connect("alien://");
				fRecalib = TFile::Open("alien:///alice/cern.ch/user/g/gconesab/RecalDB/summer2010/RecalibrationFactors.root");
			}
			else {
				AliError("Run number or pass number not found; RECALIBRATION NOT APPLIED");
				return kTRUE;
			}

		}

		TH2F * r0 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM0");
		TH2F * r1 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM1");
		TH2F * r2 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM2");
		TH2F * r3 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM3");
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(0,r0);
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(1,r1);
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(2,r2);
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(3,r3);
	}

	if(runRC > 140000){
		fRecalib = new TFile("RecalibrationFactors.root","read");
		if (fRecalib->IsZombie()){
			if (fDebugLevel>1) AliInfo("Connecting to alien to get RecalibrationFactors.root \n");

			TGrid::Connect("alien://");
			fRecalib = TFile::Open("alien:///alice/cern.ch/user/g/gconesab/RecalDB/2011_v0/RecalibrationFactors.root");
			if(!fRecalib) AliError("Recalibration file not found");
			return kFALSE;
		}

		TH2F * r0 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM0");
		TH2F * r1 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM1");
		TH2F * r2 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM2");
		TH2F * r3 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM3");
		TH2F * r4 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM4");
		TH2F * r5 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM5");
		TH2F * r6 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM6");
		TH2F * r7 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM7");
		TH2F * r8 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM8");
		TH2F * r9 = ( TH2F *)fRecalib->Get("EMCALRecalFactors_SM9");

		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(0,r0);
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(1,r1);
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(2,r2);
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(3,r3); 
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(4,r4); 
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(5,r5); 
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(6,r6); 
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(7,r7); 
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(8,r8); 
		fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(9,r9); 

	}

	return kTRUE;
}

//_____________________________________________________
void AliEMCALTenderSupply::RecalibrateCells()
{
	AliESDCaloCells *cells = fTender->GetEvent()->GetEMCALCells();

	Int_t nEMCcell = cells->GetNumberOfCells();
	Double_t calibFactor = 1.;

	for(Int_t icell=0; icell<nEMCcell; icell++){
		Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
		fEMCALGeo->GetCellIndex(cells->GetCellNumber(icell),imod,iTower,iIphi,iIeta);
		fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
		calibFactor = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
		cells->SetCell(icell,cells->GetCellNumber(icell),cells->GetAmplitude(icell)*calibFactor,cells->GetTime(icell));	
	}	

}

