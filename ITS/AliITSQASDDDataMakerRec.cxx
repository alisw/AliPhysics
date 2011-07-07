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
#include <TH2F.h>
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
  //
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Start of SDD Cycle with event specie %s for task %s\n",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie()),AliQAv1::GetTaskName(fAliITSQADataMakerRec->GetTaskIndexSelected()).Data()));
  if(!fCalibration) {
    //fAliITSQADataMakerRec->GetRawsData(9+ fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()])->Reset();
    CreateTheCalibration();
  }

  if(fAliITSQADataMakerRec->GetEventSpecie()==0) return;//not the kDefault EventSpecie
  //Detector specific actions at start of cycle
  if(fAliITSQADataMakerRec->GetTaskIndexSelected()==AliQAv1::kRAWS){
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SDD Cycle\n");
    if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kFALSE)return;
    //
    int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
    AliDebug(AliQAv1::GetQADebugLevel(),Form("Reset of Raw Data normalized histograms with eventspecie %s ",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie())));
    fAliITSQADataMakerRec->ResetRawsData(10+offsRW );
    fAliITSQADataMakerRec->ResetRawsData(3+ offsRW);
    fAliITSQADataMakerRec->ResetRawsData(4+ offsRW);
    fAliITSQADataMakerRec->ResetRawsData(5+ offsRW);
    //
  }
  if(fAliITSQADataMakerRec->GetTaskIndexSelected()==AliQAv1::kRECPOINTS){
    if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRECPOINTS)==kFALSE)return;
    
    AliDebug(AliQAv1::GetQADebugLevel(),Form("Reset of RecPoints normalized histograms with eventspecie %s ",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie())));
    int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
    fAliITSQADataMakerRec->ResetRecPointsData(9+  offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(10+ offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(11+ offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(27+ offsRP);
    
  }
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** /*list*/)
{
  //end of a QA cycle

  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  
  Double_t entriescalibration1= 0.;
  Double_t entriescalibrationL3=0.;
  Double_t entriescalibrationL4=0.;
  int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  //
  if(fHistoCalibration){
    entriescalibration1= ((TH1F*)(fHistoCalibration->At(0)))->GetEntries();
    entriescalibrationL3=((TH2F*)(fHistoCalibration->At(1)))->GetEntries();
    entriescalibrationL4=((TH2F*)(fHistoCalibration->At(2)))->GetEntries();
  }
  else{ AliWarning("Calibration TObjArray is NULL! No Normalization and calibtaion plot will be filled\n");}
  //
  for (int trCl=-1;trCl<fAliITSQADataMakerRec->GetNTrigClasses();trCl++) { // RS Loop over all trigger classes (+ global counter, -1)
    //
    if(task==AliQAv1::kRAWS) {
      TObjArray &harr = *fAliITSQADataMakerRec->GetRawsDataOfTrigClass(trCl);
      int nEvent = GetNumberOfEvents(AliQAv1::kRAWS,trCl);
      //	  printf("fNevent %d \n",nEvent);
      TH1* h10 = (TH1*)harr[10+offsRW];
      if (h10) h10->SetBinContent(1,nEvent);
      //
      if(nEvent!=0) {
	if (h10 && harr[0+offsRW] && harr[1+offsRW] && harr[2]) { // were they all cloned?
	  //
	  Double_t entriesmodpattern = ((TH1*)harr[0+offsRW])->GetEntries();
	  Double_t entriesL3         = ((TH1*)harr[1+offsRW])->GetEntries();
	  Double_t entriesL4         = ((TH1*)harr[2+offsRW])->GetEntries();
	  //	  
	  Double_t normentriesmodpattern=0.; 
	  Double_t normentriesL3= 0.;       
	  Double_t normentriesL4= 0.;     
	  
	  if(entriescalibration1!=0.)normentriesmodpattern=entriesmodpattern/(entriescalibration1*nEvent);
	  if(entriesL3!=0.)   normentriesL3= entriesL3/(entriescalibrationL3*nEvent);
	  if(entriesL4!=0.)   normentriesL4= entriesL4/(entriescalibrationL4*nEvent);
	  //
	  h10->SetBinContent(2,normentriesmodpattern);
	  h10->SetBinContent(3,normentriesL3);
	  h10->SetBinContent(4,normentriesL4);
	  //	
	  if(fHistoCalibration){
	    if (harr[3+offsRW] && harr[0+offsRW]) ((TH1*)harr[3+offsRW])->Divide( (TH1*)harr[0+offsRW],((TH1F*)(fHistoCalibration->At(0))),1.,nEvent);
	    if (harr[4+offsRW] && harr[1+offsRW]) ((TH1*)harr[4+offsRW])->Divide( (TH1*)harr[1+offsRW],((TH1F*)(fHistoCalibration->At(1))),1.,nEvent);
	    if (harr[5+offsRW] && harr[2+offsRW]) ((TH1*)harr[5+offsRW])->Divide( (TH1*)harr[2+offsRW],((TH1F*)(fHistoCalibration->At(2))),1.,nEvent);
	  }
	}	  
      }
      //
      TH1* htmp4 = (TH1*)harr[4+offsRW];
      Int_t xbin3 = htmp4 ? ((TH1*)harr[4+offsRW])->GetNbinsX() : 0;
      Int_t ybin3 = htmp4 ? ((TH1*)harr[4+offsRW])->GetNbinsY() : 0;
      //
      TH1* htmp6 = ((TH1*)harr[6+offsRW]);
      if (htmp6 && htmp4) for(Int_t i=0; i<xbin3; i++) for(Int_t j=0; j<ybin3; j++) htmp6->Fill(htmp4->GetBinContent(i,j));
      //
      TH1* htmp5 = (TH1*)harr[5+offsRW];
      Int_t xbin4 = htmp5 ? htmp5->GetNbinsX() : 0;
      Int_t ybin4 = htmp5 ? htmp5->GetNbinsY() : 0;
      //
      TH1* htmp7 = (TH1*)harr[7+offsRW];
      if (htmp7 && htmp5) for(Int_t i=0; i<xbin4; i++) for(Int_t j=0; j<ybin4; j++) htmp7->Fill(htmp5->GetBinContent(i,j));
      //
      if(fHistoCalibration){
	for(Int_t inumb=8;inumb<10;inumb++) {
	  TH1* htmp = (TH1*)harr[inumb + offsRW];
	  if (!htmp) continue;
	  TH1* hcal = (TH1*)fHistoCalibration->At(inumb-7);
	  for(Int_t ii=1; ii<htmp->GetNbinsX()+1;ii++) {
	    for(Int_t jj=1; jj<htmp->GetNbinsY()+1;jj++) {	      
	      if( hcal->GetBinContent(ii,jj) != 0. ) htmp->SetBinContent(ii,jj,1);
	      else if( hcal->GetBinContent(ii,jj) == 0. ) htmp->SetBinContent(ii,jj,0);
	    }//end x axis
	  }//end y axis
	}//end for inumb
      }
      //      
    }//end raws
    //
    //
    if(task==AliQAv1::kRECPOINTS) {
      //
      TObjArray &harr = *fAliITSQADataMakerRec->GetRecPointsDataOfTrigClass(trCl);
      //
      TH1* htmp0  = (TH1*)harr[0+offsRP];
      TH1* htmp1  = (TH1*)harr[1+offsRP];
      TH1* htmp13 = (TH1*)harr[13+offsRP];
      TH1* htmp14 = (TH1*)harr[14+offsRP];
      //
      int nEventRP = GetNumberOfEvents(AliQAv1::kRECPOINTS,trCl);
      if (htmp0 && htmp1 && htmp13 && htmp14) {
	Double_t chargeL3    = htmp0->GetMean();
	Double_t errchargeL3 = htmp0->GetMeanError();
	//
	Double_t chargeL4    = htmp1->GetMean();
	Double_t errchargeL4 = htmp1->GetMeanError();
	//
	Double_t radiusL3    = htmp13->GetMean();
	Double_t errradiusL3 = htmp13->GetMeanError();
	
	Double_t radiusL4= htmp14->GetMean();
	Double_t errradiusL4= htmp14->GetMeanError();
	//
	// printf("fNeventRP %d \n",nEventRP);
	TH1* htmp27 = (TH1*)harr[27+offsRP];
	if (htmp27) {
	  htmp27->SetBinContent(1,nEventRP);
	  
	  htmp27->SetBinContent(31,chargeL3);
	  htmp27->SetBinError(31,errchargeL3);
	  
	  htmp27->SetBinContent(47,chargeL4);
	  htmp27->SetBinError(47,errchargeL4);
	  
	  htmp27->SetBinContent(32,radiusL3);
	  htmp27->SetBinError(32,errradiusL3);
	  
	  htmp27->SetBinContent(48,radiusL4);
	  htmp27->SetBinError(48,errradiusL4);
	}
	//
	TH1 *htmp6=(TH1*)harr[6+offsRP],*htmp7=(TH1*)harr[7+offsRP],*htmp8=(TH1*)harr[8+offsRP];
	Double_t entriesmodpattern = htmp6 ? htmp6->GetEntries() : 0;
	Double_t entriesL3         = htmp7 ? htmp7->GetEntries() : 0;
	Double_t entriesL4         = htmp8 ? htmp8->GetEntries() : 0;
    
	//printf("entries modpatternrp = %02f \t L3= %02f \t L4 =%02f \n",entriesmodpattern,entriesL3,entriesL4 );
    
	if(nEventRP!=0) {
	  Double_t normentriesmodpattern=0.; 
	  Double_t normentriesL3= 0.;       
	  Double_t normentriesL4= 0.;     
	  
	  if(entriescalibration1!=0.) normentriesmodpattern=entriesmodpattern/(entriescalibration1*nEventRP);
	  if(entriesL3!=0.)           normentriesL3= entriesL3/(entriescalibrationL3*nEventRP);
	  if(entriesL4!=0.)           normentriesL4= entriesL4/(entriescalibrationL4*nEventRP);
      
	  if (htmp27) {
	    htmp27->SetBinContent(2,normentriesmodpattern);
	    htmp27->SetBinContent(3,normentriesL3);
	    htmp27->SetBinContent(4,normentriesL4);	      
	  }
	  // printf("NORM entries modpatternrp = %02f \t L3= %02f \t L4 =%02f \n",normentriesmodpattern,normentriesL3,normentriesL4 );
      
	  if(fHistoCalibration) {
	    if (harr[9 +offsRP] && htmp6) ((TH1*)harr[9 +offsRP])->Divide(htmp6, (TH1*)fHistoCalibration->At(0),nEventRP);
	    if (harr[10+offsRP] && htmp7) ((TH1*)harr[10+offsRP])->Divide(htmp7, (TH1*)fHistoCalibration->At(1),nEventRP);
	    if (harr[11+offsRP] && htmp8) ((TH1*)harr[11+offsRP])->Divide(htmp8, (TH1*)fHistoCalibration->At(2),nEventRP);
	    //if (harr[21+offsRP] && htmp7) ((TH1*)harr[21+offsRP])->Divide(htmp7, (TH1*)fHistoCalibration->At(1),nEventRP);
	    //if (harr[22+offsRP] && htmp8) ((TH1*)harr[22+offsRP])->Divide(htmp8, (TH1*)fHistoCalibration->At(2),nEventRP);
	  }
	}
	TH1* htmp10 = (TH1*)harr[10+offsRP];
	Int_t xbin3 = htmp10 ? htmp10->GetNbinsX() : 0;
	Int_t ybin3 = htmp10 ? htmp10->GetNbinsY() : 0;
	//
	TH1* htmp19 = (TH1*)harr[19+offsRP];	
	if (htmp19&&htmp10) for(Int_t i=0; i<xbin3; i++) for(Int_t j=0; j<ybin3; j++) htmp19->Fill(htmp10->GetBinContent(i,j));
	//
	TH1* htmp11 = (TH1*)harr[11+offsRP];
	Int_t xbin4 = htmp11 ? htmp11->GetNbinsX() : 0;
	Int_t ybin4 = htmp11 ? htmp11->GetNbinsY() : 0;
	//
	TH1* htmp20 = (TH1*)harr[20+offsRP];
	if (htmp11&&htmp20) for(Int_t i=0; i<xbin4; i++) for(Int_t j=0; j<ybin4; j++) htmp20->Fill(htmp11->GetBinContent(i,j));
      }
      //
      
      // RecPoints 2 Raws Ratio
      if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kTRUE) {
	//
	TH1* htmp21  = (TH1*)harr[21+offsRP];
	TH1* htmp23  = (TH1*)harr[23+offsRP];
	TH1* htmp10  = (TH1*)harr[10+offsRP];
	TH1* htmp4rw =  fAliITSQADataMakerRec->GetRawsData(4 + offsRW,trCl);
	if (htmp21 && htmp23 && htmp10 && htmp4rw) {
	  Int_t xbin3RP = htmp21->GetNbinsX();
	  Int_t ybin3RP = htmp21->GetNbinsY();
	  Int_t xbin3R  = htmp4rw->GetNbinsX();
	  Int_t ybin3R  = htmp4rw->GetNbinsY();
	  if((xbin3RP == xbin3R) && (ybin3RP == ybin3R)) {
	    htmp21->Divide(htmp10,htmp4rw);
	    for(Int_t i=0; i<xbin3R; i++) for(Int_t j=0; j<ybin3R; j++) htmp23->Fill(htmp21->GetBinContent(i,j));
	  }
	  else AliWarning("Number of bins for Raws and RecPoints (Layer 3) do not match\n");
	}
	//
	TH1* htmp22  = (TH1*)harr[22+offsRP];
	TH1* htmp24  = (TH1*)harr[24+offsRP];
	TH1* htmp11  = (TH1*)harr[11+offsRP];
	TH1* htmp5rw =  fAliITSQADataMakerRec->GetRawsData(5 + offsRW,trCl);
	if (htmp22 && htmp24 && htmp11 && htmp5rw) {
	  Int_t xbin4RP = htmp22->GetNbinsX();
	  Int_t ybin4RP = htmp22->GetNbinsY();
	  Int_t xbin4R  = htmp5rw->GetNbinsX();
	  Int_t ybin4R  = htmp5rw->GetNbinsY();
	  if((xbin4RP == xbin4R) && (ybin4RP == ybin4R)) {
	    htmp22->Divide(htmp11,htmp5rw);
	    for(Int_t i=0; i<xbin4R; i++) for(Int_t j=0; j<ybin4R; j++) htmp24->Fill(htmp22->GetBinContent(i,j));
	  } 
	  else AliWarning("Number of bins for Raws and RecPoints (Layer 4) do not match\n");
	}
      }
      else{AliWarning("Ratio between RecPoints and Raws not executed because the raw list has not been created\n");}
    }//end recpoints
    //
  } // RS Loop over all trigger classes (+ global counter, -1)
  //
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
  int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  fSDDhRawsTask = 0;
  if(fkOnline){AliInfo("Book Online Histograms for SDD\n");}
  else {AliInfo("Book Offline Histograms for SDD\n ");}
  TH1F *h0 = new TH1F("SDDModPattern","HW Modules pattern",fgknSDDmodules,239.5,499.5); //0
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  h0->SetOption("bar1");
  h0->SetBarOffset(0.01);
  h0->SetBarWidth(0.95);
  h0->SetFillColor(45);
  rv = fAliITSQADataMakerRec->Add2RawsList(h0,0+offsRW, expert, !image, !saveCorr);
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2F *hphil3 = new TH2F("SDDphizL3","SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);//1
  hphil3->GetXaxis()->SetTitle("z[Module Number L3 ]");
  hphil3->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  hphil3->SetStats(0);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil3,1+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2F *hphil4 = new TH2F("SDDphizL4","SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5); //2
  hphil4->GetXaxis()->SetTitle("z[Module Number L4]");
  hphil4->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
  hphil4->SetStats(0);
   rv = fAliITSQADataMakerRec->Add2RawsList(hphil4,2+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;
  
  //normalized histograms
  TH1F *h0norm = new TH1F("SDDModPatternNORM","NORM HW Modules pattern",fgknSDDmodules,239.5,499.5); //3
  h0norm->GetXaxis()->SetTitle("Module Number");
  h0norm->GetYaxis()->SetTitle("Counts");
  h0norm->SetOption("bar1");
  h0norm->SetBarOffset(0.01);
  h0norm->SetBarWidth(0.95);
  h0norm->SetFillColor(46);
  //h0norm->SetStats(0);
  rv = fAliITSQADataMakerRec->Add2RawsList(h0norm,3+offsRW, expert, !image, !saveCorr);
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2F *hphil3norm = new TH2F("SDDphizL3NORM","NORM SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);//4
  hphil3norm->GetXaxis()->SetTitle("z[Module Number L3 ]");
  hphil3norm->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  hphil3norm->SetStats(0);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil3norm,4+offsRW, expert, !image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2F *hphil4norm = new TH2F("SDDphizL4NORM","NORM SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5); //5
  hphil4norm->GetXaxis()->SetTitle("z[Module Number L4]");
  hphil4norm->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
  hphil4norm->SetStats(0);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil4norm,5+offsRW, expert, !image, !saveCorr); 
  fSDDhRawsTask++;

	
  Float_t hMax = 0.2;
  
  TH1F *oL3 = new TH1F("SDDL3_RelativeOccupancy","Layer 3 Relative Occupancy",200,0.,hMax);
  rv = fAliITSQADataMakerRec->Add2RawsList(oL3,6+offsRW, expert, !image, !saveCorr); //6  
  fSDDhRawsTask++;
  
  TH1F *oL4 = new TH1F("SDDL4_RelativeOccupancy","Layer 4 Relative Occupancy",200,0.,hMax);
  rv = fAliITSQADataMakerRec->Add2RawsList(oL4,7+offsRW, expert, !image, !saveCorr); //7   
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2F *hcalibl3 = new TH2F("SDDphizCalibL3","SDDCalibL3 ",12,0.5,6.5,14,0.5,14.5);//8
  hcalibl3->GetXaxis()->SetTitle("z[Module Number L3]");
  hcalibl3->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  hcalibl3->SetStats(0);
  hcalibl3->SetMaximum(2);
  rv = fAliITSQADataMakerRec->Add2RawsList(hcalibl3,8+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2F *hcalibl4 = new TH2F("SDDphizCalibL4","SDDCalibL4 ",16,0.5,8.5,22,0.5,22.5); //9
  hcalibl4->GetXaxis()->SetTitle("z[Module Number L4]");
  hcalibl4->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
  hcalibl4->SetStats(0);
  hcalibl4->SetMaximum(2);
  rv = fAliITSQADataMakerRec->Add2RawsList(hcalibl4,9+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;

  TH1F *hsummarydata = new TH1F("SDDRawDataCheck","SDDRawDataCheck",46,-0.5,45.5);//10 summary of raw data checks
  hsummarydata->GetXaxis()->SetLabelSize(0.02);
  hsummarydata->GetXaxis()->SetTickLength(0.01);
  hsummarydata->GetXaxis()->SetNdivisions(110);
  hsummarydata->GetXaxis()->SetTicks("-");
  hsummarydata->GetYaxis()->SetLabelSize(0.02);
  hsummarydata->GetYaxis()->SetTitleSize(0.02);
  hsummarydata->GetYaxis()->SetTitleOffset(1.5);
  hsummarydata->GetYaxis()->SetTitle("#events(norm) or #modules (m) or #drift regions (dr)");
  hsummarydata->SetStats(0);
  hsummarydata->SetMaximum(272);
  hsummarydata->SetMinimum(0.);

  hsummarydata->SetOption("textbar1");
  hsummarydata->SetBarOffset(0.01);
  hsummarydata->SetBarWidth(0.95);
  hsummarydata->SetFillColor(8);

  //information about the events
  hsummarydata->GetXaxis()->SetBinLabel(1,"Events");
  hsummarydata->GetXaxis()->SetBinLabel(2,"Ent_NORM");
  hsummarydata->GetXaxis()->SetBinLabel(3,"Ent_NORML3");
  hsummarydata->GetXaxis()->SetBinLabel(4,"Ent_NORML4");
  //all

  hsummarydata->GetXaxis()->SetBinLabel(5,"m_act");
  hsummarydata->GetXaxis()->SetBinLabel(6,"m_fil");
  hsummarydata->GetXaxis()->SetBinLabel(7,"dr_act");
  hsummarydata->GetXaxis()->SetBinLabel(8,"dr_fil");
  hsummarydata->GetXaxis()->SetBinLabel(9,"m_exc");
  hsummarydata->GetXaxis()->SetBinLabel(10,"m_emp");
  hsummarydata->GetXaxis()->SetBinLabel(11,"dr_exc");
  hsummarydata->GetXaxis()->SetBinLabel(12,"dr_emp");
  hsummarydata->GetXaxis()->SetBinLabel(13,"m_exc-act");
  hsummarydata->GetXaxis()->SetBinLabel(14,"m_act-emp");
  hsummarydata->GetXaxis()->SetBinLabel(15,"dr_exc-act");
  hsummarydata->GetXaxis()->SetBinLabel(16,"dr_act-emp");
  hsummarydata->GetXaxis()->SetBinLabel(17,"m_overth");
  hsummarydata->GetXaxis()->SetBinLabel(18,"dr_overth");

  //l3

  hsummarydata->GetXaxis()->SetBinLabel(19,"m_actl3");
  hsummarydata->GetXaxis()->SetBinLabel(20,"m_fill3");
  hsummarydata->GetXaxis()->SetBinLabel(21,"dr_actl3");
  hsummarydata->GetXaxis()->SetBinLabel(22,"dr_fill3");
  hsummarydata->GetXaxis()->SetBinLabel(23,"m_excl3");
  hsummarydata->GetXaxis()->SetBinLabel(24,"m_empl3");
  hsummarydata->GetXaxis()->SetBinLabel(25,"dr_excl3");
  hsummarydata->GetXaxis()->SetBinLabel(26,"dr_empl3");
  hsummarydata->GetXaxis()->SetBinLabel(27,"m_exc-actl3");
  hsummarydata->GetXaxis()->SetBinLabel(28,"m_act-empl3");
  hsummarydata->GetXaxis()->SetBinLabel(29,"dr_exc-actl3");
  hsummarydata->GetXaxis()->SetBinLabel(30,"dr_act-empl3");
  hsummarydata->GetXaxis()->SetBinLabel(31,"m_overthl3");
  hsummarydata->GetXaxis()->SetBinLabel(32,"dr_overthl3");

  //l4

  hsummarydata->GetXaxis()->SetBinLabel(33,"m_actl4");
  hsummarydata->GetXaxis()->SetBinLabel(34,"m_fill4");
  hsummarydata->GetXaxis()->SetBinLabel(35,"dr_actl4");
  hsummarydata->GetXaxis()->SetBinLabel(36,"dr_fill4");
  hsummarydata->GetXaxis()->SetBinLabel(37,"m_excl4");
  hsummarydata->GetXaxis()->SetBinLabel(38,"m_empl4");
  hsummarydata->GetXaxis()->SetBinLabel(39,"dr_excl4");
  hsummarydata->GetXaxis()->SetBinLabel(40,"dr_empl4");
  hsummarydata->GetXaxis()->SetBinLabel(41,"m_exc-actl4");
  hsummarydata->GetXaxis()->SetBinLabel(42,"m_act-empl4");
  hsummarydata->GetXaxis()->SetBinLabel(43,"dr_exc-actl4");
  hsummarydata->GetXaxis()->SetBinLabel(44,"dr_act-empl4");
  hsummarydata->GetXaxis()->SetBinLabel(45,"m_overthl4");
  hsummarydata->GetXaxis()->SetBinLabel(46,"dr_overthl4");
  hsummarydata->GetXaxis()->LabelsOption("v");

  rv = fAliITSQADataMakerRec->Add2RawsList(hsummarydata,10+offsRW, expert, !image, !saveCorr); 
  fSDDhRawsTask++;
  fOnlineOffsetRaws = fSDDhRawsTask;

  //online part
  if(fkOnline){
      //DDL Pattern 
    TH2F *hddl = new TH2F("SDDDDLPattern","SDD DDL Pattern ",24,-0.5,11.5,24,-0.5,23.5); //11
    hddl->GetXaxis()->SetTitle("Channel");
    hddl->GetYaxis()->SetTitle("DDL Number");
    rv = fAliITSQADataMakerRec->Add2RawsList(hddl,fOnlineOffsetRaws+offsRW, expert, !image, !saveCorr);
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
	rv = fAliITSQADataMakerRec->Add2RawsList(fModuleChargeMapFSE,indexlast1 + index1 + offsRW, expert, !image, !saveCorr);	  
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
	rv = fAliITSQADataMakerRec->Add2RawsList(fModuleChargeMap,indexlast1 + index1 + offsRW, expert, !image, !saveCorr); 
	fSDDhRawsTask++;
	index1++;	 
      }
    }
    
    //Event Size 
    TH1F *hsize = new TH1F("SDDEventSize","SDD Event Size ",500,-0.5,199.5); 
    hsize->SetBit(TH1::kCanRebin);
    hsize->GetXaxis()->SetTitle("Event Size [kB]");
    hsize->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RawsList(hsize,indexlast1 + index1 + offsRW, expert, !image, !saveCorr); 
    fSDDhRawsTask++;
    
  }  // kONLINE
  
  cout << fSDDhRawsTask << " SDD Raws histograms booked" << endl;
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Raws histograms booked\n",fSDDhRawsTask));
  //
  return rv ; 
}


//____________________________________________________________________________
Int_t AliITSQASDDDataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW - SDD -
  Int_t rv = 0;
  int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  // Check id histograms already created for this Event Specie
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
	if(fSDDhRawsTask > fOnlineOffsetRaws + 1 + index) fAliITSQADataMakerRec->ResetRawsData(fOnlineOffsetRaws + 1 + index +offsRW);   
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
    fAliITSQADataMakerRec->FillRawsData( 0 + offsRW,isddmod);   
    if(lay==3) fAliITSQADataMakerRec->FillRawsData(1+offsRW,det+0.5*iside-0.25,lad); //phiz l3 not norm
    if(lay==4) fAliITSQADataMakerRec->FillRawsData(2+offsRW,det+0.5*iside-0.25,lad); //phiz l4 not norm
 
    if(fkOnline) {

      fAliITSQADataMakerRec->FillRawsData(fOnlineOffsetRaws+offsRW, (stream->GetCarlosId())+0.5*iside -0.5,iddl);
      //  printf("content ddlmap %d, %d = %f \n",(stream->GetCarlosId()+0.5*iside -0.5),iddl,fAliITSQADataMakerRec->GetRawsData(3+offsRW)->GetBinContent(1+(det-1)*2;lad));
      //printf("content ddlmap %d, %d = %f \n",(stream->GetCarlosId())+0.5*iside -0.5,iddl,fAliITSQADataMakerRec->GetRawsData(3+offsRW)->GetBinContent(1+(stream->GetCarlosId()-1)*2,iddl));
      activeModule = moduleSDD;
      index1 = activeModule * 2 + iside;
      
      if(index1<0){
        AliDebug(AliQAv1::GetQADebugLevel(),Form("Wrong index number %d - patched to 0\n",index1));
		  index1 = 0;
      }      

      if(fSDDhRawsTask > fOnlineOffsetRaws +1 + index1) {                                  
        fAliITSQADataMakerRec->FillRawsData(fOnlineOffsetRaws +1 + index1 +offsRW        , coord2, coord1, signal);     
        fAliITSQADataMakerRec->FillRawsData(fOnlineOffsetRaws +1 + index1 + 260*2 +offsRW, coord2, coord1, signal); 
      }
    }//online
    cnt++;
    if(!(cnt%10000)) AliDebug(AliQAv1::GetQADebugLevel(),Form(" %d raw digits read",cnt));
  }//end next()
  if(fkOnline){
    fAliITSQADataMakerRec->FillRawsData(fOnlineOffsetRaws +1 + 260*4 +offsRW,size/1024.);//KB
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
  //
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
  int offsD = fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
//  fGenDigitsOffset = (fAliITSQADataMakerRec->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSDDhTask must be incremented by one unit every time a histogram is ADDED to the QA List
  TH1F* h0=new TH1F("SDD DIGITS Module Pattern","SDD DIGITS Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h0,offsD, !expert, image);
  fSDDhDigitsTask ++;
  // printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(offsD)->GetName() , fSDDhDigitsTask );
  TH1F* h1=new TH1F("SDD Anode Distribution","DIGITS Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h1,1+offsD, !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(1+offsD)->GetName() , fSDDhDigitsTask );
  TH1F* h2=new TH1F("SDD Tbin Distribution","DIGITS Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h2,2+offsD, !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(2+offsD)->GetName() , fSDDhDigitsTask );
  TH1F* h3=new TH1F("SDD ADC Counts Distribution","DIGITS ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h3,3+offsD, !expert, image);
  fSDDhDigitsTask ++;
  //printf("Add %s \t the task offset is %i\n",fAliITSQADataMakerRec->GetDigitsData(3+offsD)->GetName() , fSDDhDigitsTask );
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Digits histograms booked\n",fSDDhDigitsTask));
  //
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
  int offsD = fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
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
    fAliITSQADataMakerRec->FillDigitsData(offsD,nmod,ndigits);

    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      fAliITSQADataMakerRec->FillDigitsData(1+offsD,iz);
      fAliITSQADataMakerRec->FillDigitsData(2+offsD,ix);
      fAliITSQADataMakerRec->FillDigitsData(3+offsD,sig);
    }
  }
  //
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
  int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
//  fGenRecPointsOffset = (fAliITSQADataMakerRec->fRecPointsQAList[AliRecoParam::kDefault])->GetEntries();

  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("offsRP %d\n",offsRP));
  TH1F *h0 = new TH1F("SDDLay3TotCh","Layer 3 total charge",250,-0.5, 499.5); //position number 0
  //h0->SetBit(TH1::kCanRebin);
  h0->GetXaxis()->SetTitle("KeV");
  h0->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h0, 0 +offsRP, !expert, image);//NON expert image
  fSDDhRecPointsTask++;
 
  TH1F *h1 = new TH1F("SDDLay4TotCh","Layer 4 total charge",250,-0.5, 499.5);//position number 1
  //h1->SetBit(TH1::kCanRebin);
  h1->GetXaxis()->SetTitle("Kev");
  h1->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h1, 1 +offsRP, !expert, image);//NON expert image
  fSDDhRecPointsTask++;

  TH2F *h2 = new TH2F("SDDGlobalCoordDistribYX","YX Global Coord Distrib",56,-28,28,56,-28,28);//position number 2
  h2->GetYaxis()->SetTitle("Y[cm]");
  h2->GetXaxis()->SetTitle("X[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h2,2+offsRP, !expert, image);// NON expert image
  fSDDhRecPointsTask++;

  TH2F *h3 = new TH2F("SDDGlobalCoordDistribRZ","RZ Global Coord Distrib",128,-32,32,56,12,26);//position number 3
  h3->GetYaxis()->SetTitle("R[cm]");
  h3->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h3,3+offsRP, !expert, image);// NON expert image
  fSDDhRecPointsTask++;
  
  TH2F *h4 = new TH2F("SDDGlobalCoordDistribL3PHIZ","#varphi Z Global Coord Distrib L3",96,-23,23,112,-TMath::Pi(),TMath::Pi());//position number 4
  h4->GetYaxis()->SetTitle("#phi[rad]");
  h4->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h4,4+offsRP, !expert, image);//NON expert image
  fSDDhRecPointsTask++;

  TH2F *h5 = new TH2F("SDDGlobalCoordDistribL4PHIZ","#varphi Z Global Coord Distrib L4",128,-31,31,176,-TMath::Pi(),TMath::Pi());//position number 5
  h5->GetYaxis()->SetTitle("#phi[rad]");
  h5->GetXaxis()->SetTitle("Z[cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h5,5+offsRP, !expert, image);//NON expert image
  fSDDhRecPointsTask++;
  
  TH1F *h6 = new TH1F("SDDModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5); //position number 6
  h6->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h6->GetYaxis()->SetTitle("Entries");
  h6->SetOption("bar1");
  h6->SetBarOffset(0.01);
  h6->SetBarWidth(0.95);
  h6->SetFillColor(39);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h6,6 +offsRP, expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  
  TH2F *h7 = new TH2F("SDDModPatternL3RP","Modules pattern L3 RP",12,0.5,6.5,14,0.5,14.5);  //position number 7
  h7->GetXaxis()->SetTitle("z[#Module L3 ]");
  h7->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h7,7 +offsRP, expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  TH2F *h8 = new TH2F("SDDModPatternL4RP","Modules pattern L4 RP",16,0.5,8.5,22,0.5,22.5); //position number 8
  h8->GetXaxis()->SetTitle("[#Module L3 ]");
  h8->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h8,8 +offsRP, expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  //------------------------norm--------------------------//


  TH1F *h9 = new TH1F("SDDModPatternRPNORM","Modules pattern RP NORM",fgknSDDmodules,239.5,499.5); //position number 9
  h9->GetXaxis()->SetTitle("Module number"); //spd offset = 240
  h9->GetYaxis()->SetTitle("Entries");
  h9->SetOption("bar1");
  h9->SetBarOffset(0.01);
  h9->SetBarWidth(0.95);
  h9->SetFillColor(49);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h9,9 +offsRP, expert, !image);// expert NO image
  fSDDhRecPointsTask++;

  
  TH2F *h10 = new TH2F("SDDModPatternL3RPNORM","Modules pattern L3 RP NORM",12,0.5,6.5,14,0.5,14.5);  //position number 10
  h10->GetXaxis()->SetTitle("z[#Module L3 ]");
  h10->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h10,10 +offsRP, expert, !image);// expert NO  image
  fSDDhRecPointsTask++;

  TH2F *h11 = new TH2F("SDDModPatternL4RPNORM","Modules pattern L4 RP NORM",16,0.5,8.5,22,0.5,22.5); //position number 11
  h11->GetXaxis()->SetTitle("[#Module L4 ]");
  h10->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h11,11 +offsRP, expert, !image);//  expert NO image
  fSDDhRecPointsTask++;

  //--------------------------------------------------------//

  TH2F *h12 = new TH2F("SDDLocalCoordDistrib","Local Coord Distrib",160,-4,4,160,-4,4);//position number 12
  h12->GetXaxis()->SetTitle("X local coord, drift, cm");
  h12->GetYaxis()->SetTitle("Z local coord, anode, cm");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h12,12 +offsRP, expert, !image);//  expert  NO image
  fSDDhRecPointsTask++;
  
  //AliInfo("Create SDD recpoints histos\n");
  
  TH1F *h13 = new TH1F("SDDrdistrib_Layer3" ,"SDD r distribution Layer3" ,100,14.,16.5);//position number 13 (L3)
  h13->SetBit(TH1::kCanRebin);
  h13->GetXaxis()->SetTitle("r[cm]");
  h13->GetXaxis()->CenterTitle();
  h13->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h13,13 +offsRP, expert, !image);// expert NO image
  fSDDhRecPointsTask++;
  
  TH1F *h14 = new TH1F("SDDrdistrib_Layer4" ,"SDD r distribution Layer4" ,100,23.,25.);// and position number 14 (L4)
  h14->SetBit(TH1::kCanRebin);
  h14->GetXaxis()->SetTitle("r[cm]");
  h14->GetXaxis()->CenterTitle();
  h14->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h14,14 +offsRP, expert, !image);// expert NO image
  fSDDhRecPointsTask++;
  
  for(Int_t iLay=0; iLay<=1; iLay++){
    TH1F *h15 = new TH1F(Form("SDDphidistrib_Layer%d",iLay+3),Form("SDDphidistrib_Layer%d",iLay+3) ,180,-TMath::Pi(),TMath::Pi());//position number 15 (L3) and position number 16 (L4)
    h15->GetXaxis()->SetTitle("#varphi[rad]");
    h15->GetXaxis()->CenterTitle();
    h15->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(h15,iLay+15+offsRP, expert, !image);// expert NO image
    fSDDhRecPointsTask++;
  }
  
  for(Int_t iLay=0; iLay<=1; iLay++){
    TH1F *h17 = new TH1F(Form("SDDdrifttime_Layer%d",iLay+3),Form("SDDdrifttime_Layer%d",iLay+3),45,-0.5,4499.5);//position number 17 (L3) and position number 18 (L4)
    h17->SetBit(TH1::kCanRebin);
    h17->GetXaxis()->SetTitle("drift time[ns]");
    h17->GetXaxis()->CenterTitle();
    h17->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(h17,iLay+17+offsRP, !expert, image);// NON expert  image
    fSDDhRecPointsTask++;
  }
  
	Float_t hMax = 0.2;
	
	TH1F *oL3 = new TH1F("SDDL3_RelativeOccupancy","Layer 3 Relative Occupancy (RecPoints)",200,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(oL3,fSDDhRecPointsTask+offsRP, expert, !image); // 19
	fSDDhRecPointsTask++;
	
	TH1F *oL4 = new TH1F("SDDL4_RelativeOccupancy","Layer 4 Relative Occupancy (RecPoints)",200,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(oL4,fSDDhRecPointsTask+offsRP, expert, !image); // 20
	fSDDhRecPointsTask++;
	
	
	TH2F *h21 = new TH2F("SDDL3_Rec2Raw_2D","L3 RecPoints to Raws 2D",12,0.5,6.5,14,0.5,14.5);  //position number 21
	h21->GetXaxis()->SetTitle("z[#Module L3 ]");
	h21->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
	rv = fAliITSQADataMakerRec->Add2RecPointsList(h21,21 +offsRP, expert, !image);// expert NO  image
	fSDDhRecPointsTask++;
	
	TH2F *h22 = new TH2F("SDDL4_Rec2Raw_2D","L4 RecPoints to Raws 2D",16,0.5,8.5,22,0.5,22.5); //position number 22
	h22->GetXaxis()->SetTitle("[#Module L4 ]");
	h22->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
	rv = fAliITSQADataMakerRec->Add2RecPointsList(h22,22 +offsRP, expert, !image);//  expert NO image
	fSDDhRecPointsTask++;

        hMax = 0.3;	
	TH1F *R2RL3 = new TH1F("SDDL3_Rec2Raw","L3 RecPoints to Raws ratio",150,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(R2RL3,fSDDhRecPointsTask+offsRP, expert, !image); // 23
	fSDDhRecPointsTask++;
	
	TH1F *R2RL4 = new TH1F("SDDL4_Rec2Raw","L4 RecPoints to Raws ratio",150,0.,hMax);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(R2RL4,fSDDhRecPointsTask+offsRP, expert, !image); // 24
	fSDDhRecPointsTask++;

	TH1F *dedxL3 = new TH1F("SDDL3_dedx","L3 dE/dX",100,0.,1.);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(dedxL3,fSDDhRecPointsTask+offsRP, expert, !image); // 25
	fSDDhRecPointsTask++;
	
	TH1F *dedxL4 = new TH1F("SDDL4_dedx","L4 dE/dX",100,0.,1.);
	rv = fAliITSQADataMakerRec->Add2RecPointsList(dedxL4,fSDDhRecPointsTask+offsRP, expert, !image); // 26
	fSDDhRecPointsTask++;
	

  TH1F *hsummarydatarp = new TH1F("SDDRecPointCheck","SDDRecPointCheck",50,-0.5,49.5);//27 summary of recpoint checks
  hsummarydatarp->GetXaxis()->SetLabelSize(0.02);
  hsummarydatarp->GetXaxis()->SetTickLength(0.01);
  hsummarydatarp->GetXaxis()->SetNdivisions(110);
  hsummarydatarp->GetXaxis()->SetTicks("-");
  hsummarydatarp->GetYaxis()->SetLabelSize(0.02);
  hsummarydatarp->GetYaxis()->SetTitleSize(0.02);
  hsummarydatarp->GetYaxis()->SetTitleOffset(1.5);
  hsummarydatarp->GetYaxis()->SetTitle("#events(norm) or #modules (m) or #drift regions (dr)");
  hsummarydatarp->SetStats(0);
  hsummarydatarp->SetMaximum(272);

  hsummarydatarp->SetOption("text bar1 hist");
  hsummarydatarp->SetBarOffset(0.05);
  hsummarydatarp->SetBarWidth(0.95);
  hsummarydatarp->SetFillColor(32);
  hsummarydatarp->SetMinimum(0.);

  //information about the events
  hsummarydatarp->GetXaxis()->SetBinLabel(1,"Events");
  hsummarydatarp->GetXaxis()->SetBinLabel(2,"Ent_NORM");
  hsummarydatarp->GetXaxis()->SetBinLabel(3,"Ent_NORML3");
  hsummarydatarp->GetXaxis()->SetBinLabel(4,"Ent_NORML4");
  //all

  hsummarydatarp->GetXaxis()->SetBinLabel(5, "m_act");
  hsummarydatarp->GetXaxis()->SetBinLabel(6, "m_fil");
  hsummarydatarp->GetXaxis()->SetBinLabel(7, "dr_act");
  hsummarydatarp->GetXaxis()->SetBinLabel(8, "dr_fil");
  hsummarydatarp->GetXaxis()->SetBinLabel(9, "m_exc");
  hsummarydatarp->GetXaxis()->SetBinLabel(10,"m_emp");
  hsummarydatarp->GetXaxis()->SetBinLabel(11,"dr_exc");
  hsummarydatarp->GetXaxis()->SetBinLabel(12,"dr_emp");
  hsummarydatarp->GetXaxis()->SetBinLabel(13,"m_exc-act");
  hsummarydatarp->GetXaxis()->SetBinLabel(14,"m_act-emp");
  hsummarydatarp->GetXaxis()->SetBinLabel(15,"dr_exc-act");
  hsummarydatarp->GetXaxis()->SetBinLabel(16,"dr_act-emp");
  hsummarydatarp->GetXaxis()->SetBinLabel(17,"m_overth");
  hsummarydatarp->GetXaxis()->SetBinLabel(18,"dr_overth");

  //l3

  hsummarydatarp->GetXaxis()->SetBinLabel(19,"m_actl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(20,"m_fill3");
  hsummarydatarp->GetXaxis()->SetBinLabel(21,"dr_actl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(22,"dr_fill3");
  hsummarydatarp->GetXaxis()->SetBinLabel(23,"m_excl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(24,"m_empl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(25,"dr_excl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(26,"dr_empl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(27,"m_exc-actl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(28,"m_act-empl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(29,"dr_exc-actl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(30,"dr_act-empl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(31,"meanch_l3");
  hsummarydatarp->GetXaxis()->SetBinLabel(32,"meanrad_l3");
  hsummarydatarp->GetXaxis()->SetBinLabel(33,"m_overthl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(34,"dr_overthl3");

  //l4

  hsummarydatarp->GetXaxis()->SetBinLabel(35,"m_actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(36,"m_fill4");
  hsummarydatarp->GetXaxis()->SetBinLabel(37,"dr_actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(38,"dr_fill4");
  hsummarydatarp->GetXaxis()->SetBinLabel(39,"m_excl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(40,"m_empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(41,"dr_excl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(42,"dr_empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(43,"m_exc-actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(44,"m_act-empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(45,"dr_exc-actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(46,"dr_act-empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(47,"meanch_l4");
  hsummarydatarp->GetXaxis()->SetBinLabel(48,"meanrad_l4");
  hsummarydatarp->GetXaxis()->SetBinLabel(49,"m_overthl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(50,"dr_overthl4");

  hsummarydatarp->GetXaxis()->LabelsOption("v");

  rv = fAliITSQADataMakerRec->Add2RecPointsList(hsummarydatarp,fSDDhRecPointsTask+offsRP, expert, !image);// expert NO image



  fSDDhRecPointsTask++;

  fOnlineOffsetRecPoints = fSDDhRecPointsTask;
  if(fkOnline){
    TH2F *h19 = new TH2F("SDDGlobalCoordDistribYXFSE","YX Global Coord Distrib FSE",112,-28,28,112,-28,28);//position number 27
    h19->GetYaxis()->SetTitle("Y[cm]");
    h19->GetXaxis()->SetTitle("X[cm]");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(h19,fSDDhRecPointsTask+offsRP, expert, !image);// expert NO image
    fSDDhRecPointsTask++;
    
    TH2F *h20 = new TH2F("SDDGlobalCoordDistribRZFSE","RZ Global Coord Distrib FSE",128,-32,32,56,12,26);//position number 28
    h20->GetYaxis()->SetTitle("R[cm]");
    h20->GetXaxis()->SetTitle("Z[cm]");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(h20,fSDDhRecPointsTask+offsRP, expert, !image);// expert NO image
    fSDDhRecPointsTask++;      
  }//online
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Recs histograms booked\n",fSDDhRecPointsTask));
  //
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::MakeRecPoints(TTree * clustersTree)
{
 // Fill QA for RecPoints - SDD -
  Int_t rv = 0 ;
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
  int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  if(fkOnline){
    for(Int_t i=28;i<30;i++) fAliITSQADataMakerRec->ResetRecPointsData(i+offsRP);
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
      fAliITSQADataMakerRec->FillRecPointsData(6 +offsRP,modnumb);//modpatternrp
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      Float_t drifttime=recp->GetDriftTime();
      fAliITSQADataMakerRec->FillRecPointsData(12 +offsRP,recp->GetDetLocalX(),recp->GetDetLocalZ());//local distribution
      fAliITSQADataMakerRec->FillRecPointsData(2 +offsRP,cluglo[0],cluglo[1]);//global distribution YX
      fAliITSQADataMakerRec->FillRecPointsData(3 +offsRP,cluglo[2],rad);//global distribution rz
      if(fkOnline) {
	fAliITSQADataMakerRec->FillRecPointsData(28 +offsRP,cluglo[0],cluglo[1]);//global distribution YX FSE
	fAliITSQADataMakerRec->FillRecPointsData(29 +offsRP,cluglo[2],rad);//global distribution rz FSE
      }
      Int_t iside=recp->GetDriftSide();
      lay=recp->GetLayer();
      if(lay == 2) {
	fAliITSQADataMakerRec->FillRecPointsData(0  +offsRP, recp->GetQ()) ;//total charge of layer 3
	fAliITSQADataMakerRec->FillRecPointsData(7  +offsRP, det+0.5*iside-0.5,lad);//mod pattern layer 3
	fAliITSQADataMakerRec->FillRecPointsData(13 +offsRP, rad);//r distribution layer 3
	fAliITSQADataMakerRec->FillRecPointsData(15 +offsRP, phi);// phi distribution layer 3
	fAliITSQADataMakerRec->FillRecPointsData(4  +offsRP, cluglo[2],phi);// zphi distribution layer
	fAliITSQADataMakerRec->FillRecPointsData(17  +offsRP, drifttime);// time distribution layer 3
	fAliITSQADataMakerRec->FillRecPointsData(25  +offsRP, recp->GetdEdX());// charge distribution layer 3
      } else if(lay == 3) {
	fAliITSQADataMakerRec->FillRecPointsData(1  +offsRP, recp->GetQ()) ;//total charge layer 4
	fAliITSQADataMakerRec->FillRecPointsData(8  +offsRP, det+0.5*iside-0.5,lad);//mod pattern layer 4
	fAliITSQADataMakerRec->FillRecPointsData(14 +offsRP, rad);//r distribution
	fAliITSQADataMakerRec->FillRecPointsData(16 +offsRP, phi);//phi distribution
	fAliITSQADataMakerRec->FillRecPointsData(5  +offsRP, cluglo[2],phi);// zphi distribution layer 4
	fAliITSQADataMakerRec->FillRecPointsData(18  +offsRP, drifttime);// time distribution layer 4
	fAliITSQADataMakerRec->FillRecPointsData(26  +offsRP, recp->GetdEdX());// charge distribution layer 4
      }
    }
  }
  //
  return rv ; 
}

//_______________________________________________________________
Int_t AliITSQASDDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task, Int_t specie) const
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
	      ((TH1F*)(fHistoCalibration->At(0)))->SetBinContent(imod+1,fillmodhisto1);
	      ((TH2F*)(fHistoCalibration->At(1)))->SetBinContent(index,lad,fillmodhisto2side0);
	      ((TH2F*)(fHistoCalibration->At(1)))->SetBinContent(index+1,lad,fillmodhisto2side1);
	    }//end layer 3
	    else if(lay==4){
	      badmodhisto1=badmodhisto3side0+badmodhisto3side1;
	      fillmodhisto1-=badmodhisto1;
	      fillmodhisto3side0-=badmodhisto3side0;
	      fillmodhisto3side1-=badmodhisto3side1;
	      ((TH1F*)(fHistoCalibration->At(0)))->SetBinContent(imod+1,fillmodhisto1);
	      ((TH2F*)(fHistoCalibration->At(2)))->SetBinContent(index,lad,fillmodhisto3side0);
	      ((TH2F*)(fHistoCalibration->At(2)))->SetBinContent(index+1,lad,fillmodhisto3side1);
	    }//end layer 4
	  }//end else bad module
	}//end module for
    }

}

//____________________________________________________________________
void AliITSQASDDDataMakerRec::InitCalibrationArray()
{
  //create the histograms with the calibration informations. The histograms are stored in a TObjArray
    TH1F *pattern1  = new TH1F("CALSDDModPattern","Calibration HW Modules pattern",fgknSDDmodules,239.5,499.5);
    pattern1->SetDirectory(0) ;
    TH2F *patternl3 = new TH2F("CALSDDphizL3","Calibration SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);
    patternl3->SetDirectory(0) ;
    TH2F *patternl4 = new TH2F("CALSDDphizL4"," Calibration SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5);
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

  ((TH1F*)(fHistoCalibration->At(0)))->Reset();
  ((TH2F*)(fHistoCalibration->At(1)))->Reset();
  ((TH2F*)(fHistoCalibration->At(2)))->Reset();
  //delete fHistoCalibration;
  //fHistoCalibration=NULL;
  
}

//____________________________________________________________________

Int_t AliITSQASDDDataMakerRec::GetNumberOfEvents(AliQAv1::TASKINDEX_t task, Int_t trigCl)
{
  //return the number of the processed events for a given task and trigger class (-1 for all)
  return fAliITSQADataMakerRec->GetEvCountTotal(task, trigCl);
}
