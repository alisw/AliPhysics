/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  M.Siciliano Aug 2008 QA RecPoints 
//  INFN Torino

// --- ROOT system ---

#include <TProfile2D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TBranch.h>
#include <TTree.h>
#include <TMath.h>
//#include <TObjArray.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQASDDDataMakerRec.h"
#include "AliQAv1.h"
#include "AliRawReader.h"
#include "AliITSRawStream.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSdigit.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeomTGeo.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"

class TGaxis;
class TF1;
class TSystem;
class AliLog;
class AliQAChecker;
class AliITSRawStreamSDDCompressed;
class AliCDBStorage;
class Riostream;
class AliITSdigitSDD;
class AliITS;
class AliRunLoader;
class AliITSLoader;
class AliITSDetTypeRec;



ClassImp(AliITSQASDDDataMakerRec)

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSDDhRawsTask(0),
fSDDhDigitsTask(0),
fSDDhRecPointsTask(0),
fOnlineOffsetRaws(0),
fOnlineOffsetRecPoints(0),
fGenRawsOffset(0),
fGenDigitsOffset(0),
fGenRecPointsOffset(0),
fTimeBinSize(1),
fNEvent(0),
fNEventRP(0),
fDDLModuleMap(0),
fCalibration(0),
fHistoCalibration(0)
{
  //ctor used to discriminate OnLine-Offline analysis
  if(fLDC < 0 || fLDC > 6) {
	AliError("Error: LDC number out of range; return\n");
  }
	fGenRawsOffset = new Int_t[AliRecoParam::kNSpecies];
	fGenRecPointsOffset = new Int_t[AliRecoParam::kNSpecies];
	fGenDigitsOffset = new Int_t[AliRecoParam::kNSpecies];
	for(Int_t i=0; i<AliRecoParam::kNSpecies; i++) {
		fGenRawsOffset[i] = 0;
		fGenRecPointsOffset[i] = 0;
		fGenDigitsOffset[i]=0;
	}

	InitCalibrationArray();
}

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(const AliITSQASDDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSDDhRawsTask(qadm.fSDDhRawsTask),
fSDDhDigitsTask(qadm.fSDDhDigitsTask),
fSDDhRecPointsTask(qadm.fSDDhRecPointsTask),
fOnlineOffsetRaws(qadm.fOnlineOffsetRaws),
fOnlineOffsetRecPoints(qadm.fOnlineOffsetRecPoints),
fGenRawsOffset(qadm.fGenRawsOffset),
fGenDigitsOffset(qadm.fGenDigitsOffset),
fGenRecPointsOffset(qadm.fGenRecPointsOffset),
fTimeBinSize(qadm.fTimeBinSize),
fNEvent(qadm.fNEvent),
fNEventRP(qadm.fNEventRP),
fDDLModuleMap(qadm.fDDLModuleMap),
fCalibration(qadm.fCalibration),
fHistoCalibration(qadm.fHistoCalibration)
{
  //copy ctor 
  fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
  fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
  //fDDLModuleMap=NULL;
}

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::~AliITSQASDDDataMakerRec(){
  // destructor
  //if(fDDLModuleMap) delete fDDLModuleMap;
  if(fHistoCalibration){delete fHistoCalibration; fHistoCalibration=NULL;}
}
//__________________________________________________________________
AliITSQASDDDataMakerRec& AliITSQASDDDataMakerRec::operator = (const AliITSQASDDDataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQASDDDataMakerRec();
  new(this) AliITSQASDDDataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::StartOfDetectorCycle()
{

  //Start of a QA cycle

  AliDebug(AliQAv1::GetQADebugLevel(),Form("Start of SDD Cycle with event specie %s for task %s\n",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie()),AliQAv1::GetTaskName(fAliITSQADataMakerRec->GetTaskIndexSelected()).Data()));
  if(!fCalibration) {CreateTheCalibration();}

  if(fAliITSQADataMakerRec->GetEventSpecie()==0) return;//not the kDefault EventSpecie
    //Detector specific actions at start of cycle
    if(fAliITSQADataMakerRec->GetTaskIndexSelected()==AliQAv1::kRAWS){
      AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SDD Cycle\n");
      if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kFALSE)return;

	AliDebug(AliQAv1::GetQADebugLevel(),Form("Reset of Raw Data normalized histograms with eventspecie %s ",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie())));

	fAliITSQADataMakerRec->GetRawsData(3+ fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
	fAliITSQADataMakerRec->GetRawsData(4+ fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
	fAliITSQADataMakerRec->GetRawsData(5+ fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
      
    }
    if(fAliITSQADataMakerRec->GetTaskIndexSelected()==AliQAv1::kRECPOINTS){
      if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRECPOINTS)==kFALSE)return;

	AliDebug(AliQAv1::GetQADebugLevel(),Form("Reset of RecPoints normalized histograms with eventspecie %s ",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie())));

	fAliITSQADataMakerRec->GetRecPointsData(9+  fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
	fAliITSQADataMakerRec->GetRecPointsData(10+ fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
	fAliITSQADataMakerRec->GetRecPointsData(11+ fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
      }
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray* /*list*/)
{
  //end of a QA cycle
	AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
	if(task==AliQAv1::kRAWS){
	  //	  printf("fNevent %d \n",fNEvent);
	  if(fNEvent!=0){
	    ((TH1D*)fAliITSQADataMakerRec->GetRawsData(3 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH1D*)fAliITSQADataMakerRec->GetRawsData(0 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH1D*) (fHistoCalibration->At(0))),1.,(Double_t)fNEvent);
	  
	    ((TH2D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH2D*)fAliITSQADataMakerRec->GetRawsData(1 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH2D*)(fHistoCalibration->At(1))),1.,(Double_t)fNEvent);
	    
	    ((TH2D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH2D*)fAliITSQADataMakerRec->GetRawsData(2 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH2D*)(fHistoCalibration->At(2))),1.,(Double_t)fNEvent);
	  }	  

		Int_t xbin3 = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		Int_t ybin3 = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		
		for(Int_t i=0; i<xbin3; i++) {
			for(Int_t j=0; j<ybin3; j++) {
				((TH1D*)fAliITSQADataMakerRec->GetRawsData(6 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Fill(((TH1D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetBinContent(i,j));
			}
		}
		
		Int_t xbin4 = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		Int_t ybin4 = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		
		for(Int_t i=0; i<xbin4; i++) {
			for(Int_t j=0; j<ybin4; j++) {
				((TH1D*)fAliITSQADataMakerRec->GetRawsData(7 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Fill(((TH1D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetBinContent(i,j));
			}
		}
			
	}//end raws
	
	if(task==AliQAv1::kRECPOINTS){
	  //	  printf("fNeventRP %d \n",fNEventRP);
	  if(fNEventRP!=0){
	    ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(9 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(6 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH1D*) (fHistoCalibration->At(0))),1.,(Double_t)fNEventRP);
	    
		  ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(10+ fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(7 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH2D*)(fHistoCalibration->At(1))),1.,(Double_t)fNEventRP);
		  
		  ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(11+ fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(8 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH2D*)(fHistoCalibration->At(2))),1.,(Double_t)fNEventRP);

		  //  ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(21+ fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(7 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH2D*)(fHistoCalibration->At(1))),1.,(Double_t)fNEventRP);
		  
		  // ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(22+ fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(8 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]),((TH2D*)(fHistoCalibration->At(2))),1.,(Double_t)fNEventRP);
	  }
		Int_t xbin3 = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(10 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		Int_t ybin3 = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(10 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		
		for(Int_t i=0; i<xbin3; i++) {
			for(Int_t j=0; j<ybin3; j++) {
				((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(19 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Fill(((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(10 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetBinContent(i,j));
			}
		}
		
		Int_t xbin4 = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(11 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		Int_t ybin4 = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(11 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		
		for(Int_t i=0; i<xbin4; i++) {
			for(Int_t j=0; j<ybin4; j++) {
				((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(20 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Fill(((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(11 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetBinContent(i,j));
			}
		}

		// RecPoints 2 Raws Ratio
		if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kTRUE)
		  {
		    Int_t xbin3RP = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(21 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		    Int_t ybin3RP = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(21 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		    Int_t xbin3R  = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		    Int_t ybin3R  = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		    if((xbin3RP == xbin3R) && (ybin3RP == ybin3R)) {
		      ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(21 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide(((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(10 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])),((TH2D*)fAliITSQADataMakerRec->GetRawsData(4 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])));
		      for(Int_t i=0; i<xbin3R; i++) {
			for(Int_t j=0; j<ybin3R; j++) {
			  ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(23 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Fill(((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(21 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetBinContent(i,j));
			}
		      }
		    } else 
		      AliWarning("Number of bins for Raws and RecPoints (Layer 3) do not match\n");
		    
		    Int_t xbin4RP = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(22 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		    Int_t ybin4RP = ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(22 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		    Int_t xbin4R = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsX();
		    Int_t ybin4R = ((TH1D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetNbinsY();
		    if((xbin4RP == xbin4R) && (ybin4RP == ybin4R)) {
 ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(22 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide(((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(11 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])),((TH2D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])));
 //		      ((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(22 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Divide(((TH2D*)fAliITSQADataMakerRec->GetRawsData(5 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])));
		      for(Int_t i=0; i<xbin4R; i++) {
			for(Int_t j=0; j<ybin4R; j++) {
			  ((TH1D*)fAliITSQADataMakerRec->GetRecPointsData(24 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->Fill(((TH2D*)fAliITSQADataMakerRec->GetRecPointsData(22 + fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]))->GetBinContent(i,j));
			}
		      }
		    } else 
		      AliWarning("Number of bins for Raws and RecPoints (Layer 4) do not match\n");
		  }
		else{AliWarning("Ratio between RecPoints and Raws not executed because the raw list has not been created\n");}
	}//end recpoints
	
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitRaws()
{ 

  // Initialization for RAW data - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  Int_t rv = 0 ; 
  Int_t lay, lad, det;
  Int_t indexlast = 0;
  Int_t index1 = 0;

	fSDDhRawsTask = 0;
	if(fkOnline){AliInfo("Book Online Histograms for SDD\n");}
  else {AliInfo("Book Offline Histograms for SDD\n ");}
  TH1D *h0 = new TH1D("SDDModPattern","HW Modules pattern",fgknSDDmodules,239.5,499.5); //0
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  rv = fAliITSQADataMakerRec->Add2RawsList(h0,0+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2D *hphil3 = new TH2D("SDDphizL3","SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);//1
  hphil3->GetXaxis()->SetTitle("z[Module Number L3 ]");
  hphil3->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil3,1+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2D *hphil4 = new TH2D("SDDphizL4","SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5); //2
  hphil4->GetXaxis()->SetTitle("z[Module Number L4]");
  hphil4->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
   rv = fAliITSQADataMakerRec->Add2RawsList(hphil4,2+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr); 
  fSDDhRawsTask++;
  
  //normalized histograms
  TH1D *h0norm = new TH1D("SDDModPatternNORM","NORM HW Modules pattern",fgknSDDmodules,239.5,499.5); //3
  h0norm->GetXaxis()->SetTitle("Module Number");
  h0norm->GetYaxis()->SetTitle("Counts");
  rv = fAliITSQADataMakerRec->Add2RawsList(h0norm,3+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2D *hphil3norm = new TH2D("SDDphizL3NORM","NORM SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);//4
  hphil3norm->GetXaxis()->SetTitle("z[Module Number L3 ]");
  hphil3norm->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil3norm,4+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2D *hphil4norm = new TH2D("SDDphizL4NORM","NORM SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5); //5
  hphil4norm->GetXaxis()->SetTitle("z[Module Number L4]");
  hphil4norm->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
   rv = fAliITSQADataMakerRec->Add2RawsList(hphil4norm,5+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image, !saveCorr); 
  fSDDhRawsTask++;

	
	Float_t hMax = 0.2;
	
	TH1F *oL3 = new TH1F("SDDL3_RelativeOccupancy","Layer 3 Relative Occupancy",200,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RawsList(oL3,6+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr); //6  
	fSDDhRawsTask++;
	
	TH1F *oL4 = new TH1F("SDDL4_RelativeOccupancy","Layer 4 Relative Occupancy",200,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RawsList(oL4,7+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr); //7   
	fSDDhRawsTask++;
	
	fOnlineOffsetRaws = fSDDhRawsTask;
	if(fkOnline){
      //DDL Pattern 
      TH2D *hddl = new TH2D("SDDDDLPattern","SDD DDL Pattern ",24,-0.5,11.5,24,-0.5,23.5); //8
      hddl->GetXaxis()->SetTitle("Channel");
      hddl->GetYaxis()->SetTitle("DDL Number");
      rv = fAliITSQADataMakerRec->Add2RawsList(hddl,fOnlineOffsetRaws+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);
      fSDDhRawsTask++;
      Int_t indexlast1 = 0;
  
      fTimeBinSize = 4;
      indexlast = 0;
      index1 = 0;
      indexlast1 = fSDDhRawsTask;
      for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
		  for(Int_t iside=0;iside<fgknSide;iside++){
			  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
			  TProfile2D *fModuleChargeMapFSE = new TProfile2D(Form("SDDchargeMapFSE_L%d_%d_%d_%d",lay,lad,det,iside),Form("SDDChargeMapForSingleEvent_L%d_%d_%d_%d",lay,lad,det,iside)  ,256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
			  fModuleChargeMapFSE->GetXaxis()->SetTitle("Time Bin");
			  fModuleChargeMapFSE->GetYaxis()->SetTitle("Anode");
			  rv = fAliITSQADataMakerRec->Add2RawsList(fModuleChargeMapFSE,indexlast1 + index1 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr);	  
			  fSDDhRawsTask++;
			  index1++;	 
		  }
      }
      
      for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
		  for(Int_t iside=0;iside<fgknSide;iside++){
			  AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
			  TProfile2D *fModuleChargeMap = new TProfile2D(Form("SDDchargeMap_L%d_%d_%d_%d",lay,lad,det,iside),Form("SDDChargeMap_L%d_%d_%d_%d",lay,lad,det,iside),256/fTimeBinSize,-0.5,255.5,256,-0.5,255.5);
			  fModuleChargeMap->GetXaxis()->SetTitle("Time Bin");
			  fModuleChargeMap->GetYaxis()->SetTitle("Anode Number");
			  rv = fAliITSQADataMakerRec->Add2RawsList(fModuleChargeMap,indexlast1 + index1 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr); 
			  fSDDhRawsTask++;
			  index1++;	 
		  }
      }
      
      //Event Size 
      TH1F *hsize = new TH1F("SDDEventSize","SDD Event Size ",500,-0.5,199.5); 
      hsize->SetBit(TH1::kCanRebin);
      hsize->GetXaxis()->SetTitle("Event Size [kB]");
      hsize->GetYaxis()->SetTitle("Entries");
      rv = fAliITSQADataMakerRec->Add2RawsList(hsize,indexlast1 + index1 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image, !saveCorr); 
      fSDDhRawsTask++;
	  
    }  // kONLINE
	
	cout << fSDDhRawsTask << " SDD Raws histograms booked" << endl;
	
	AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Raws histograms booked\n",fSDDhRawsTask));
	return rv ; 
}


//____________________________________________________________________________
Int_t AliITSQASDDDataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW - SDD -
  Int_t rv = 0;
  // Check id histograms already created for this Event Specie
  fNEvent++;
  if(!fDDLModuleMap){CreateTheMap();}
  if(rawReader->GetType() != 7) return rv;  // skips non physical triggers
  AliDebug(AliQAv1::GetQADebugLevel(),"entering MakeRaws\n");                 
  rawReader->Reset();       
  rawReader->Select("ITSSDD");
  AliITSRawStream *stream=AliITSRawStreamSDD::CreateRawStreamSDD(rawReader);
   stream->SetDDLModuleMap(fDDLModuleMap);
  
  Int_t lay, lad, det; 
  
  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0;iside<fgknSide;iside++) {
		if(fSDDhRawsTask > fOnlineOffsetRaws + 1 + index) fAliITSQADataMakerRec->GetRawsData(fOnlineOffsetRaws + 1 + index +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();   
	// 4  because the 2D histos for single events start after the fourth position
	index++;
      }
    }
  }
  
  Int_t cnt = 0;
  Int_t ildcID = -1;
  Int_t iddl = -1;
  Int_t isddmod = -1;
  Int_t coord1, coord2, signal, moduleSDD, activeModule, index1; 
  //if(fkOnline)
  //{
      Int_t prevDDLID = -1;
      UInt_t size = 0;
      int totalddl=static_cast<int>(fDDLModuleMap->GetNDDLs());
      Bool_t *ddldata=new Bool_t[totalddl];
      for(Int_t jddl=0;jddl<totalddl;jddl++){ddldata[jddl]=kFALSE;}
      //}
  while(stream->Next()) {
    ildcID = rawReader->GetLDCId();
    iddl = rawReader->GetDDLID();// - fgkDDLIDshift;
    if(iddl<0)isddmod=-1;
    //printf("----------------------iddl %i\n",iddl);
    else isddmod = fDDLModuleMap->GetModuleNumber(iddl,stream->GetCarlosId());

    if(isddmod==-1){
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Found module with iddl: %d, stream->GetCarlosId: %d \n",iddl,stream->GetCarlosId()));
      continue;
    }
    if(stream->IsCompletedModule()) {
      AliDebug(AliQAv1::GetQADebugLevel(),Form("IsCompletedModule == KTRUE\n"));
      continue;
    } 
    if(stream->IsCompletedDDL()) {
      if(fkOnline){
	if ((rawReader->GetDDLID() != prevDDLID)&&(ddldata[iddl])==kFALSE){
	  size += rawReader->GetDataSize();//in bytes
	  prevDDLID = rawReader->GetDDLID();
	  ddldata[iddl]=kTRUE;
	}
      }
      AliDebug(AliQAv1::GetQADebugLevel(),Form("IsCompletedDDL == KTRUE\n"));
      continue;
    } 
    
    coord1 = stream->GetCoord1();
    coord2 = stream->GetCoord2();
    signal = stream->GetSignal();
    
    moduleSDD = isddmod - fgkmodoffset;
    
    if(isddmod <fgkmodoffset|| isddmod>fgknSDDmodules+fgkmodoffset-1) {
      AliDebug(AliQAv1::GetQADebugLevel(),Form( "Module SDD = %d, resetting it to 1 \n",isddmod));
      isddmod = 1;
    }
    
    AliITSgeomTGeo::GetModuleId(isddmod, lay, lad, det);
    Short_t iside = stream->GetChannel();

    //printf(" \n%i %i %i %i \n ",lay, lad, det,iside );
    fAliITSQADataMakerRec->GetRawsData( 0 + fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()] )->Fill(isddmod);   
    if(lay==3) fAliITSQADataMakerRec->GetRawsData(1+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(det+0.5*iside-0.25,lad); 
    if(lay==4) fAliITSQADataMakerRec->GetRawsData(2+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(det+0.5*iside-0.25,lad);
 
    if(fkOnline) {

      fAliITSQADataMakerRec->GetRawsData(fOnlineOffsetRaws+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill((stream->GetCarlosId())+0.5*iside -0.5,iddl);
      //  printf("content ddlmap %d, %d = %f \n",(stream->GetCarlosId()+0.5*iside -0.5),iddl,fAliITSQADataMakerRec->GetRawsData(3+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetBinContent(1+(det-1)*2;lad));
      //printf("content ddlmap %d, %d = %f \n",(stream->GetCarlosId())+0.5*iside -0.5,iddl,fAliITSQADataMakerRec->GetRawsData(3+fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetBinContent(1+(stream->GetCarlosId()-1)*2,iddl));
      activeModule = moduleSDD;
      index1 = activeModule * 2 + iside;
      
      if(index1<0){
        AliDebug(AliQAv1::GetQADebugLevel(),Form("Wrong index number %d - patched to 0\n",index1));
		  index1 = 0;
      }      

      if(fSDDhRawsTask > fOnlineOffsetRaws +1 + index1) {                                  
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(fOnlineOffsetRaws +1 + index1 +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])))->Fill(coord2, coord1, signal);     
        ((TProfile2D *)(fAliITSQADataMakerRec->GetRawsData(fOnlineOffsetRaws +1 + index1 + 260*2 +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])))->Fill(coord2, coord1, signal); 
      }
    }//online
    cnt++;
    if(!(cnt%10000)) AliDebug(AliQAv1::GetQADebugLevel(),Form(" %d raw digits read",cnt));
  }//end next()
  if(fkOnline){((TH1F*)(fAliITSQADataMakerRec->GetRawsData(fOnlineOffsetRaws +1 + 260*4 +fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])))->Fill(size/1024.);//KB
  }
	
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Event completed, %d raw digits read",cnt)); 
  delete stream;
  stream = NULL; 

//	if(fkOnline) {
//		AnalyseBNG(); // Analyse Baseline, Noise, Gain
//		AnalyseINJ(); // Analyse Injectors
//	}

  delete []ddldata;
  ddldata=NULL;
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitDigits()
{ 
  //  if(!fHistoCalibration)InitCalibrationArray();
  // Initialization for DIGIT data - SDD -  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
//  fGenDigitsOffset = (fAliITSQADataMakerRec->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  TH1F* h0=new TH1F("SDD DIGITS Module Pattern","SDD DIGITS Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h0,fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  // printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  TH1F* h1=new TH1F("SDD Anode Distribution","DIGITS Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h1,1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  TH1F* h2=new TH1F("SDD Tbin Distribution","DIGITS Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h2,2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  TH1F* h3=new TH1F("SDD ADC Counts Distribution","DIGITS ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h3,3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->GetName() , fSDDhDigitsTask );
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Digits histograms booked\n",fSDDhDigitsTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASDDDataMakerRec::MakeDigits(TTree * digits)
{ 

  // Fill QA for DIGIT - SDD -
  //AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  //fITS->SetTreeAddress();
  //TClonesArray *iITSdigits  = fITS->DigitsAddress(1);


  Int_t rv = 0 ; 

  TBranch *branchD = digits->GetBranch("ITSDigitsSDD");

  if (!branchD) {
    AliError("can't get the branch with the ITS SDD digits !");
    return rv ;
  }
  // Check id histograms already created for this Event Specie
//  if ( ! fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset) )
//    rv = InitDigits() ;
  
  static TClonesArray statDigits("AliITSdigitSDD");
  TClonesArray *iITSdigits = &statDigits;
  branchD->SetAddress(&iITSdigits);

  for(Int_t i=0; i<260; i++){
    Int_t nmod=i+240;
    digits->GetEvent(nmod);
    Int_t ndigits = iITSdigits->GetEntries();
    fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(nmod,ndigits);

    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      fAliITSQADataMakerRec->GetDigitsData(1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(iz);
      fAliITSQADataMakerRec->GetDigitsData(2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(ix);
      fAliITSQADataMakerRec->GetDigitsData(3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(sig);
    }
  }
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitRecPoints()
{
  //if(!fHistoCalibration)InitCalibrationArray();
	//AliInfo("Initialize SDD recpoints histos\n");
  // Initialization for RECPOINTS - SDD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
//  fGenRecPointsOffset = (fAliITSQADataMakerRec->fRecPointsQAList[AliRecoParam::kDefault])->GetEntries();

  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  TH1F *h0 = new TH1F("SDDLay3TotCh","Layer 3 total charge",250,-0.5, 499.5); //position number 0
  //h0->SetBit(TH1::kCanRebin);
  h0->GetXaxis()->SetTitle("KeV");
  h0->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h0, 0 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);//NON expert image
  fSDDhRecPointsTask++;
 
  TH1F *h1 = new TH1F("SDDLay4TotCh","Layer 4 total charge",250,-0.5, 499.5);//position number 1
  //h1->SetBit(TH1::kCanRebin);
  h1->GetXaxis()->SetTitle("Kev");
  h1->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h1, 1 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);//NON expert image
  fSDDhRecPointsTask++;

  TH2F *h2 = new TH2F("SDDGlobalCoordDistribYX","YX Global Coord Distrib",56,-28,28,56,-28,28);//position number 2
  h2->GetYaxis()->SetTitle("Y[cm]");
  h2->GetXaxis()->SetTitle("X[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h2,2+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);// NON expert image
  fSDDhRecPointsTask++;

  TH2F *h3 = new TH2F("SDDGlobalCoordDistribRZ","RZ Global Coord Distrib",128,-32,32,56,12,26);//position number 3
  h3->GetYaxis()->SetTitle("R[cm]");
  h3->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h3,3+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);// NON expert image
  fSDDhRecPointsTask++;
  
  TH2F *h4 = new TH2F("SDDGlobalCoordDistribL3PHIZ","#varphi Z Global Coord Distrib L3",96,-23,23,112,-TMath::Pi(),TMath::Pi());//position number 4
  h4->GetYaxis()->SetTitle("#phi[rad]");
  h4->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h4,4+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);//NON expert image
  fSDDhRecPointsTask++;

  TH2F *h5 = new TH2F("SDDGlobalCoordDistribL4PHIZ","#varphi Z Global Coord Distrib L4",128,-31,31,176,-TMath::Pi(),TMath::Pi());//position number 5
  h5->GetYaxis()->SetTitle("#phi[rad]");
  h5->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h5,5+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);//NON expert image
  fSDDhRecPointsTask++;
  
  TH1D *h6 = new TH1D("SDDModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); //position number 6
  h6->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h6->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h6,6 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  
  TH2D *h7 = new TH2D("SDDModPatternL3RP","Modules pattern L3 RP",12,0.5,6.5,14,0.5,14.5);  //position number 7
  h7->GetXaxis()->SetTitle("z[#Module L3 ]");
  h7->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h7,7 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  TH2D *h8 = new TH2D("SDDModPatternL4RP","Modules pattern L4 RP",16,0.5,8.5,22,0.5,22.5); //position number 8
  h8->GetXaxis()->SetTitle("[#Module L3 ]");
  h8->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h8,8 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  //------------------------norm--------------------------//


  TH1D *h9 = new TH1D("SDDModPatternRPNORM","Modules pattern RP NORM",fgknSDDmodules,239.5,499.5); //position number 9
  h9->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h9->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h9,9 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  
  TH2D *h10 = new TH2D("SDDModPatternL3RPNORM","Modules pattern L3 RP NORM",12,0.5,6.5,14,0.5,14.5);  //position number 10
  h10->GetXaxis()->SetTitle("z[#Module L3 ]");
  h10->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h10,10 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO  image
  fSDDhRecPointsTask++;

  TH2D *h11 = new TH2D("SDDModPatternL4RPNORM","Modules pattern L4 RP NORM",16,0.5,8.5,22,0.5,22.5); //position number 11
  h11->GetXaxis()->SetTitle("[#Module L4 ]");
  h10->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h11,11 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);//  expert NO image
  fSDDhRecPointsTask++;

  //--------------------------------------------------------//

  TH2F *h12 = new TH2F("SDDLocalCoordDistrib","Local Coord Distrib",160,-4,4,160,-4,4);//position number 12
  h12->GetXaxis()->SetTitle("X local coord, drift, cm");
  h12->GetYaxis()->SetTitle("Z local coord, anode, cm");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h12,12 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);//  expert  NO image
  fSDDhRecPointsTask++;
  
  //AliInfo("Create SDD recpoints histos\n");
  
  TH1F *h13 = new TH1F("SDDrdistrib_Layer3" ,"SDD r distribution Layer3" ,100,14.,16.5);//position number 13 (L3)
  h13->SetBit(TH1::kCanRebin);
  h13->GetXaxis()->SetTitle("r[cm]");
  h13->GetXaxis()->CenterTitle();
  h13->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h13,13 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
  fSDDhRecPointsTask++;
  
  TH1F *h14 = new TH1F("SDDrdistrib_Layer4" ,"SDD r distribution Layer4" ,100,23.,25.);// and position number 14 (L4)
  h14->SetBit(TH1::kCanRebin);
  h14->GetXaxis()->SetTitle("r[cm]");
  h14->GetXaxis()->CenterTitle();
  h14->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h14,14 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
  fSDDhRecPointsTask++;
  
  for(Int_t iLay=0; iLay<=1; iLay++){
    TH1F *h15 = new TH1F(Form("SDDphidistrib_Layer%d",iLay+3),Form("SDDphidistrib_Layer%d",iLay+3) ,180,-TMath::Pi(),TMath::Pi());//position number 15 (L3) and position number 16 (L4)
    h15->GetXaxis()->SetTitle("#varphi[rad]");
    h15->GetXaxis()->CenterTitle();
    h15->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(h15,iLay+15+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
    fSDDhRecPointsTask++;
  }
  
  for(Int_t iLay=0; iLay<=1; iLay++){
    TH1F *h17 = new TH1F(Form("SDDdrifttime_Layer%d",iLay+3),Form("SDDdrifttime_Layer%d",iLay+3),45,-0.5,4499.5);//position number 17 (L3) and position number 18 (L4)
    h17->SetBit(TH1::kCanRebin);
    h17->GetXaxis()->SetTitle("drift time[ns]");
    h17->GetXaxis()->CenterTitle();
    h17->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(h17,iLay+17+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);// NON expert  image
    fSDDhRecPointsTask++;
  }
  
	Float_t hMax = 0.2;
	
	TH1F *oL3 = new TH1F("SDDL3_RelativeOccupancy","Layer 3 Relative Occupancy (RecPoints)",200,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(oL3,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); // 19
	fSDDhRecPointsTask++;
	
	TH1F *oL4 = new TH1F("SDDL4_RelativeOccupancy","Layer 4 Relative Occupancy (RecPoints)",200,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(oL4,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); // 20
	fSDDhRecPointsTask++;
	
	
	TH2D *h21 = new TH2D("SDDL3_Rec2Raw_2D","L3 RecPoints to Raws 2D",12,0.5,6.5,14,0.5,14.5);  //position number 21
	h21->GetXaxis()->SetTitle("z[#Module L3 ]");
	h21->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
	rv = fAliITSQADataMakerRec->Add2RecPointsList(h21,21 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO  image
	fSDDhRecPointsTask++;
	
	TH2D *h22 = new TH2D("SDDL4_Rec2Raw_2D","L4 RecPoints to Raws 2D",16,0.5,8.5,22,0.5,22.5); //position number 22
	h22->GetXaxis()->SetTitle("[#Module L4 ]");
	h22->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
	rv = fAliITSQADataMakerRec->Add2RecPointsList(h22,22 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);//  expert NO image
	fSDDhRecPointsTask++;

        hMax = 0.3;	
	TH1F *R2RL3 = new TH1F("SDDL3_Rec2Raw","L3 RecPoints to Raws ratio",150,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(R2RL3,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); // 23
	fSDDhRecPointsTask++;
	
	TH1F *R2RL4 = new TH1F("SDDL4_Rec2Raw","L4 RecPoints to Raws ratio",150,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(R2RL4,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); // 24
	fSDDhRecPointsTask++;

	TH1F *dedxL3 = new TH1F("SDDL3_dedx","L3 dE/dX",100,0.,1.);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(dedxL3,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); // 25
	fSDDhRecPointsTask++;
	
	TH1F *dedxL4 = new TH1F("SDDL4_dedx","L4 dE/dX",100,0.,1.);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(dedxL4,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); // 26
	fSDDhRecPointsTask++;
	
	fOnlineOffsetRecPoints = fSDDhRecPointsTask;
	if(fkOnline){
      TH2F *h19 = new TH2F("SDDGlobalCoordDistribYXFSE","YX Global Coord Distrib FSE",112,-28,28,112,-28,28);//position number 27
      h19->GetYaxis()->SetTitle("Y[cm]");
      h19->GetXaxis()->SetTitle("X[cm]");
      rv = fAliITSQADataMakerRec->Add2RecPointsList(h19,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
      fSDDhRecPointsTask++;
      
      TH2F *h20 = new TH2F("SDDGlobalCoordDistribRZFSE","RZ Global Coord Distrib FSE",128,-32,32,56,12,26);//position number 28
      h20->GetYaxis()->SetTitle("R[cm]");
      h20->GetXaxis()->SetTitle("Z[cm]");
      rv = fAliITSQADataMakerRec->Add2RecPointsList(h20,fSDDhRecPointsTask+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);// expert NO image
      fSDDhRecPointsTask++;      
    }//online
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Recs histograms booked\n",fSDDhRecPointsTask));
  
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::MakeRecPoints(TTree * clustersTree)
{
 // Fill QA for RecPoints - SDD -
  Int_t rv = 0 ;
  fNEventRP++; 
  Int_t lay, lad, det; 
  //AliInfo("get the branch with the ITS clusters !\n");
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  TClonesArray *recpoints=NULL; 
  if(fkOnline){recpoints = rpcont->FetchClusters(0,clustersTree,fAliITSQADataMakerRec->GetEventNumber());}
  else{recpoints = rpcont->FetchClusters(0,clustersTree);}
  AliDebug(10,Form("Fetched RecPoints for %d SDD modules",recpoints->GetEntriesFast()));
  if(!rpcont->GetStatusOK() || !rpcont->IsSDDActive()){
    AliError("can't get SDD clusters !");
    return rv;
  }

  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  if(fkOnline){
      for(Int_t i=27;i<29;i++){
	  fAliITSQADataMakerRec->GetRecPointsData(i+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
	}
    }
  // AliITSgeomTGeo::GetModuleIndex() issues an error in case the arguments
  // are illegal and returns -1
  Int_t firMod=TMath::Max(0,AliITSgeomTGeo::GetModuleIndex(3,1,1));
  Int_t lasMod=AliITSgeomTGeo::GetModuleIndex(5,1,1);
  for(Int_t module=firMod; module<lasMod;module++){
    //AliInfo(Form("Module %d\n",module));
    recpoints = rpcont->UncheckedGetClusters(module);
    npoints += recpoints->GetEntries();
    for(Int_t j=0;j<recpoints->GetEntries();j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j); 
      Int_t index = recp->GetDetectorIndex();
      lay=recp->GetLayer();
		if(lay < 2 || lay > 3) continue;
		Int_t modnumb=index+AliITSgeomTGeo::GetModuleIndex(lay+1,1,1);
		AliITSgeomTGeo::GetModuleId(modnumb, lay, lad, det);  
//		AliInfo(Form("modnumb %d, lay %d, lad %d, det %d \n",module, lay, lad, det));
      fAliITSQADataMakerRec->GetRecPointsData(6 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(modnumb);//modpatternrp
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      Float_t drifttime=recp->GetDriftTime();
		fAliITSQADataMakerRec->GetRecPointsData(12 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetDetLocalX(),recp->GetDetLocalZ());//local distribution
		fAliITSQADataMakerRec->GetRecPointsData(2 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[0],cluglo[1]);//global distribution YX
		fAliITSQADataMakerRec->GetRecPointsData(3 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],rad);//global distribution rz
		if(fkOnline) {
			fAliITSQADataMakerRec->GetRecPointsData(27 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[0],cluglo[1]);//global distribution YX FSE
			fAliITSQADataMakerRec->GetRecPointsData(28 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],rad);//global distribution rz FSE
		}
		Int_t iside=recp->GetDriftSide();
                lay=recp->GetLayer();
		if(lay == 2) {
			fAliITSQADataMakerRec->GetRecPointsData(0  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetQ()) ;//total charge of layer 3
			fAliITSQADataMakerRec->GetRecPointsData(7  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(det+0.5*iside-0.5,lad);//mod pattern layer 3
			fAliITSQADataMakerRec->GetRecPointsData(13 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(rad);//r distribution layer 3
			fAliITSQADataMakerRec->GetRecPointsData(15 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(phi);// phi distribution layer 3
			fAliITSQADataMakerRec->GetRecPointsData(4  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],phi);// zphi distribution layer
			fAliITSQADataMakerRec->GetRecPointsData(17  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(drifttime);// time distribution layer 3
			fAliITSQADataMakerRec->GetRecPointsData(25  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetdEdX());// charge distribution layer 3
		} else if(lay == 3) {
			fAliITSQADataMakerRec->GetRecPointsData(1  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetQ()) ;//total charge layer 4
			fAliITSQADataMakerRec->GetRecPointsData(8  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(det+0.5*iside-0.5,lad);//mod pattern layer 4
			fAliITSQADataMakerRec->GetRecPointsData(14 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(rad);//r distribution
			fAliITSQADataMakerRec->GetRecPointsData(16 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(phi);//phi distribution
			fAliITSQADataMakerRec->GetRecPointsData(5  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(cluglo[2],phi);// zphi distribution layer 4
			fAliITSQADataMakerRec->GetRecPointsData(18  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(drifttime);// time distribution layer 4
			fAliITSQADataMakerRec->GetRecPointsData(26  +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Fill(recp->GetdEdX());// charge distribution layer 4
      }
    }
  }
  return rv ; 
}

//_______________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task, Int_t specie)const
{
  // Returns offset number according to the specified task
  Int_t offset=0;
  if( task == AliQAv1::kRAWS ){offset=fGenRawsOffset[specie];}
  else if(task == AliQAv1::kDIGITSR ){offset=fGenDigitsOffset[specie];}
  else if( task == AliQAv1::kRECPOINTS ){offset=fGenRecPointsOffset[specie];}
  return offset;
}

//_______________________________________________________________

void AliITSQASDDDataMakerRec::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset, Int_t specie) {
  // Set offset number according to the specified task
  if( task == AliQAv1::kRAWS ) {fGenRawsOffset[specie]=offset;}
  else if( task == AliQAv1::kDIGITSR ) {fGenDigitsOffset[specie]=offset;}
  else if( task == AliQAv1::kRECPOINTS ) {fGenRecPointsOffset[specie]=offset;}
}

//_______________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetTaskHisto(AliQAv1::TASKINDEX_t task)
{
  //return the number of histo booked for a given Task
  Int_t histotot=0;
  if( task == AliQAv1::kRAWS ){ histotot=fSDDhRawsTask ;}
  else if(task == AliQAv1::kDIGITSR) { histotot=fSDDhDigitsTask;}
  else if( task == AliQAv1::kRECPOINTS ){ histotot=fSDDhRecPointsTask;}
  else {AliInfo("No task has been selected. TaskHisto set to zero.\n");}
  return histotot;
}


//_______________________________________________________________


void AliITSQASDDDataMakerRec::CreateTheMap()
{
  //Create the SDD DDL Module Map
  AliCDBEntry *ddlMapSDD = AliCDBManager::Instance()->Get("ITS/Calib/DDLMapSDD");
  Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
  if(!ddlMapSDD){
      AliError("Calibration object retrieval failed! SDD will not be processed");
      fDDLModuleMap = NULL;
      //return rv;
    }
  else{
    fDDLModuleMap = (AliITSDDLModuleMapSDD*)ddlMapSDD->GetObject();
    if(!cacheStatus)ddlMapSDD->SetObject(NULL);
    ddlMapSDD->SetOwner(kTRUE);
    if(!cacheStatus){ delete ddlMapSDD;}
    AliInfo("DDL Map Created\n ");
  }
}

//_______________________________________________________________


void AliITSQASDDDataMakerRec::CreateTheCalibration()
{
  //Take from the OCDB the calibration information for the SDD 
    AliCDBEntry *calibSDD = AliCDBManager::Instance()->Get("ITS/Calib/CalibSDD");
    Bool_t cacheStatus = AliCDBManager::Instance()->GetCacheFlag();
    if(!calibSDD)
      {
	AliError("Calibration object retrieval failed! SDD will not be processed");
	fCalibration = NULL;
      }
    else{
      fCalibration = (TObjArray *)calibSDD->GetObject();
      
      if(!cacheStatus)calibSDD->SetObject(NULL);
      calibSDD->SetOwner(kTRUE);
      if(!cacheStatus){delete calibSDD;}
      
      AliITSCalibrationSDD * cal=NULL;
      for(Int_t imod=0;imod<fgknSDDmodules;imod++)
	{
	  //cal=NULL;
	  Int_t fillmodhisto1=fgkTotalNumberSDDAnodes;
	  Int_t fillmodhisto2side0=fgkNumberOfSDDAnodesperSide;
	  Int_t fillmodhisto2side1=fgkNumberOfSDDAnodesperSide;
	  Int_t fillmodhisto3side0=fgkNumberOfSDDAnodesperSide;
	  Int_t fillmodhisto3side1=fgkNumberOfSDDAnodesperSide;
	  
	  Int_t badmodhisto1=0;
	  Int_t badmodhisto2side0=0;
	  Int_t badmodhisto2side1=0;
	  Int_t badmodhisto3side0=0;
	  Int_t badmodhisto3side1=0;
	  //printf("imod %i\t ==== \t",imod);
	  Int_t module=imod + 240;
	  //printf("module %i\t ==== \t",module);
	  cal=(AliITSCalibrationSDD*)fCalibration->At(imod);
	  Int_t lay,lad,det;
	  AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
	  Int_t index=1+(det-1)*2;
	  if(cal==0){continue;}
	  if (cal->IsBad()){continue;}//bad module check
	  else{
	    for(Int_t i=0;i<8;i++) //check on bad chips in good modules
	      {
		if(lay==3){
		  if(cal->IsChipBad(i)){
		    if(i<4){badmodhisto2side0+=64;}
		    if(i>=4){badmodhisto2side1+=64;}
		  }//end if chip
		}//end if  layer3
		else if(lay==4){
		  if(cal->IsChipBad(i)){
		    if(i<4){badmodhisto3side0+=64;}
		    if(i>=4){badmodhisto3side1+=64;}		 
		  }//end if  chip
		}//ens if layer4
	      }//end for  chip
	    for(Int_t iAn=0; iAn<512; iAn++){//anodes loop 
	      Int_t ic=cal->GetChip(iAn);//chip with this anode number
	      if(!cal->IsChipBad(ic) && !cal->IsBad() && cal->IsBadChannel(iAn)){// good chip   good module   bad channel 
		if(lay==3){
		  if(ic<4) badmodhisto2side0++;
		  else if(ic>=4)badmodhisto2side1++;
		}//end if layer 3
		else if(lay==4){
		  if(ic<4) badmodhisto3side0++;
		  else if(ic>=4)badmodhisto3side1++;
		}//end if layer 4
	      }//end if chip module channel
	    }//end for anodes
	    if(lay==3){
	      badmodhisto1=badmodhisto2side0+badmodhisto2side1;
	      fillmodhisto1-=badmodhisto1;
	      fillmodhisto2side0-=badmodhisto2side0;
	      fillmodhisto2side1-=badmodhisto2side1;
	      ((TH1D*)(fHistoCalibration->At(0)))->SetBinContent(imod+1,fillmodhisto1);
	      ((TH2D*)(fHistoCalibration->At(1)))->SetBinContent(index,lad,fillmodhisto2side0);
	      ((TH2D*)(fHistoCalibration->At(1)))->SetBinContent(index+1,lad,fillmodhisto2side1);
	    }//end layer 3
	    else if(lay==4){
	      badmodhisto1=badmodhisto3side0+badmodhisto3side1;
	      fillmodhisto1-=badmodhisto1;
	      fillmodhisto3side0-=badmodhisto3side0;
	      fillmodhisto3side1-=badmodhisto3side1;
	      ((TH1D*)(fHistoCalibration->At(0)))->SetBinContent(imod+1,fillmodhisto1);
	      ((TH2D*)(fHistoCalibration->At(2)))->SetBinContent(index,lad,fillmodhisto3side0);
	      ((TH2D*)(fHistoCalibration->At(2)))->SetBinContent(index+1,lad,fillmodhisto2side1);
	    }//end layer 4
	  }//end else bad module
	}//end module for
    }
}

//____________________________________________________________________

void AliITSQASDDDataMakerRec::InitCalibrationArray()
{
  //create the histograms with the calibration informations. The histograms are stored in a TObjArray
    TH1D *pattern1  = new TH1D("CALSDDModPattern","Calibration HW Modules pattern",fgknSDDmodules,239.5,499.5);
    pattern1->SetDirectory(0) ;
    TH2D *patternl3 = new TH2D("CALSDDphizL3","Calibration SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);
    patternl3->SetDirectory(0) ;
    TH2D *patternl4 = new TH2D("CALSDDphizL4"," Calibration SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5);
    patternl4->SetDirectory(0) ;

    if(!fHistoCalibration)fHistoCalibration = new TObjArray(3);
    fHistoCalibration->AddAtAndExpand(pattern1,0);
    fHistoCalibration->AddAtAndExpand(patternl3,1);
    fHistoCalibration->AddAtAndExpand(patternl4,2);
    fHistoCalibration->SetOwner(kTRUE); 
    //    printf("Calibration Histograms created!\n");
}

//____________________________________________________________________

void AliITSQASDDDataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{
  //reset the SDD calibration histograms
  AliInfo(Form("Reset detector in SDD called for task index %i", task));
  if(task== AliQAv1::kRAWS ){
    fDDLModuleMap=NULL;
  }

  fCalibration=NULL;

  ((TH1D*)(fHistoCalibration->At(0)))->Reset();
  ((TH2D*)(fHistoCalibration->At(1)))->Reset();
  ((TH2D*)(fHistoCalibration->At(2)))->Reset();
  //delete fHistoCalibration;
  //fHistoCalibration=NULL;
  
}

