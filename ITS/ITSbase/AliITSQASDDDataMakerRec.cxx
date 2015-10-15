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
//  last review: F. Prino Apr 2015
//  INFN Torino

// --- ROOT system ---

#include <RVersion.h>
#include <TProfile2D.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TBranch.h>
#include <TTree.h>
#include <TMath.h>
#include <TPaveText.h>
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
#include "AliITSQADataMakerRec.h"
#include "AliLog.h"
#include <iostream>

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

using std::endl;
using std::cout;
ClassImp(AliITSQASDDDataMakerRec)

//____________________________________________________________________________ 
AliITSQASDDDataMakerRec::AliITSQASDDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc) : TObject(),
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
  fTimeBinSize(4),
  fDDLModuleMap(0),
  fCalibration(0),
  fHistoCalibration(0),
  fPulserRun(99999999)
{
  //ctor used to discriminate OnLine-Offline analysis
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
  fHistoCalibration(qadm.fHistoCalibration),
  fPulserRun(qadm.fPulserRun)
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
    CreateTheCalibration();
  }

  if(fAliITSQADataMakerRec->GetEventSpecie()==0) return;

  //Detector specific actions at start of cycle
  if(fAliITSQADataMakerRec->GetTaskIndexSelected()==AliQAv1::kRAWS){
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SDD Cycle\n");
    if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kFALSE)return;
    //
    int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
    AliDebug(AliQAv1::GetQADebugLevel(),Form("Reset of Raw Data normalized histograms with eventspecie %s ",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie())));
    fAliITSQADataMakerRec->ResetRawsData(kSDDRawDataCheck+offsRW );
    fAliITSQADataMakerRec->ResetRawsData(kSDDRawModPatternNorm+offsRW);
    fAliITSQADataMakerRec->ResetRawsData(kSDDRawLadModLay3Norm+offsRW);
    fAliITSQADataMakerRec->ResetRawsData(kSDDRawLadModLay4Norm+offsRW);
  }

  if(fAliITSQADataMakerRec->GetTaskIndexSelected()==AliQAv1::kRECPOINTS){
    if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRECPOINTS)==kFALSE)return;
    
    AliDebug(AliQAv1::GetQADebugLevel(),Form("Reset of RecPoints normalized histograms with eventspecie %s ",AliRecoParam::GetEventSpecieName(fAliITSQADataMakerRec->GetEventSpecie())));
    int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
    fAliITSQADataMakerRec->ResetRecPointsData(kSDDRecpModPatternNorm+offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(kSDDRecpLadModLay3Norm+offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(kSDDRecpLadModLay4Norm+offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(kSDDRecpDataCheck+offsRP);
  }
}

//____________________________________________________________________________ 
void AliITSQASDDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** /*list*/)
{
  //end of a QA cycle

  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  
  Double_t nGoodAnodesTot= 0.;
  Double_t nGoodAnodesL3=0.;
  Double_t nGoodAnodesL4=0.;
  int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  //
  if(fHistoCalibration){
    nGoodAnodesTot= ((TH1F*)(fHistoCalibration->At(0)))->Integral();
    nGoodAnodesL3=((TH2F*)(fHistoCalibration->At(1)))->Integral();
    nGoodAnodesL4=((TH2F*)(fHistoCalibration->At(2)))->Integral();
    AliInfo(Form("Number of good anodes: Lay3=%.0f   Lay4=%.0f   Tot=%.0f",nGoodAnodesL3,nGoodAnodesL4,nGoodAnodesTot));
  }
  else{ AliWarning("Calibration TObjArray is NULL! No Normalization and calibtaion plot will be filled\n");}
  //
  for (int trCl=-1;trCl<fAliITSQADataMakerRec->GetNTrigClasses();trCl++) { // RS Loop over all trigger classes (+ global counter, -1)
    //
    if(task==AliQAv1::kRAWS) {
      TObjArray &harr = *fAliITSQADataMakerRec->GetRawsDataOfTrigClass(trCl);
      Int_t nEvent = GetNumberOfEvents(AliQAv1::kRAWS,trCl);
      TH1* h10 = (TH1*)harr[kSDDRawDataCheck+offsRW];
      if(nEvent==0) if(h10) nEvent=h10->GetBinContent(1);
      Double_t countsMod=0.;
      Double_t countsLay3=0.;
      Double_t countsLay4=0.;
      if(harr[kSDDRawModPattern+offsRW]) countsMod=((TH1*)harr[kSDDRawModPattern+offsRW])->GetEntries();
      if(harr[kSDDRawLadModLay3+offsRW]) countsLay3=((TH1*)harr[kSDDRawLadModLay3+offsRW])->GetEntries();
      if(harr[kSDDRawLadModLay4+offsRW]) countsLay4=((TH1*)harr[kSDDRawLadModLay4+offsRW])->GetEntries();
      Double_t normCountsMod=0.;
      Double_t normCountsLay3=0.;
      Double_t normCountsLay4=0.;
      if(nEvent!=0) {
	if(nGoodAnodesTot!=0.) normCountsMod = countsMod/(nGoodAnodesTot*nEvent);
	if(nGoodAnodesL3!=0.) normCountsLay3 = countsLay3/(nGoodAnodesL3*nEvent);
	if(nGoodAnodesL4!=0.) normCountsLay4 = countsLay4/(nGoodAnodesL4*nEvent);
	if(fHistoCalibration){
	  if (harr[kSDDRawModPattern+offsRW] && harr[kSDDRawModPatternNorm+offsRW]){
	    ((TH1*)harr[kSDDRawModPatternNorm+offsRW])->Divide( (TH1*)harr[kSDDRawModPattern+offsRW],((TH1F*)(fHistoCalibration->At(0))),1.,nEvent);
	  }
	  if (harr[kSDDRawLadModLay3+offsRW] && harr[kSDDRawLadModLay3Norm+offsRW]){
	    ((TH2*)harr[kSDDRawLadModLay3Norm+offsRW])->Divide( (TH2*)harr[kSDDRawLadModLay3+offsRW],((TH2F*)(fHistoCalibration->At(1))),1.,nEvent);
	  }
	  if (harr[kSDDRawLadModLay4+offsRW] && harr[kSDDRawLadModLay4Norm+offsRW]){
	    ((TH2*)harr[kSDDRawLadModLay4Norm+offsRW])->Divide( (TH2*)harr[kSDDRawLadModLay4+offsRW],((TH2F*)(fHistoCalibration->At(2))),1.,nEvent);
	  }
	}
      }
      if (h10){
	h10->SetBinContent(1,nEvent);
	h10->SetBinContent(2,normCountsMod);
	h10->SetBinContent(3,normCountsLay3);
	h10->SetBinContent(4,normCountsLay4);
      }

      //
      if(fHistoCalibration){
	TH1* hact3 = (TH1*)harr[kActiveModLay3+offsRW];
	if(hact3){
	  TH1* hcal3 = (TH1*)fHistoCalibration->At(1);
	  for(Int_t ii=1; ii<hact3->GetNbinsX()+1;ii++) {
	    for(Int_t jj=1; jj<hact3->GetNbinsY()+1;jj++) {	      
	      if( hcal3->GetBinContent(ii,jj) > 0. ) hact3->SetBinContent(ii,jj,1);
	      else hact3->SetBinContent(ii,jj,0);
	    }
	  }
	  hact3->SetTitle(Form("SDDCalibL3 - Pulser Run %d",fPulserRun));
	}
	TH1* hact4 = (TH1*)harr[kActiveModLay4+offsRW];
	if(hact4){
	  TH1* hcal4 = (TH1*)fHistoCalibration->At(2);
	  for(Int_t ii=1; ii<hact4->GetNbinsX()+1;ii++) {
	    for(Int_t jj=1; jj<hact4->GetNbinsY()+1;jj++) {	      
	      if( hcal4->GetBinContent(ii,jj) > 0. ) hact4->SetBinContent(ii,jj,1);
	      else hact4->SetBinContent(ii,jj,0);
	    }
	  }
	  hact4->SetTitle(Form("SDDCalibL4 - Pulser Run %d",fPulserRun));
	}
      }
    }//end raws


    if(task==AliQAv1::kRECPOINTS) {
      TObjArray &harr = *fAliITSQADataMakerRec->GetRecPointsDataOfTrigClass(trCl);
      Int_t nEventRP = GetNumberOfEvents(AliQAv1::kRECPOINTS,trCl);
      TH1* htmp27 = (TH1*)harr[kSDDRecpDataCheck+offsRP];
      if(nEventRP==0) if(htmp27) nEventRP=htmp27->GetBinContent(1);
      Double_t nrecpMod=0.;
      Double_t nrecpLay3=0.;
      Double_t nrecpLay4=0.;
      if(harr[kSDDRecpModPattern+offsRP]) nrecpMod=((TH1*)harr[kSDDRecpModPattern+offsRP])->GetEntries();
      if(harr[kSDDRecpLadModLay3+offsRP]) nrecpLay3=((TH1*)harr[kSDDRecpLadModLay3+offsRP])->GetEntries();
      if(harr[kSDDRecpLadModLay4+offsRP]) nrecpLay4=((TH1*)harr[kSDDRecpLadModLay4+offsRP])->GetEntries();
      Double_t normNrecpMod=0.;
      Double_t normNrecpLay3=0.;
      Double_t normNrecpLay4=0.;
      if(nEventRP!=0) {
	if(nGoodAnodesTot!=0.) normNrecpMod = nrecpMod/(nGoodAnodesTot*nEventRP);
	if(nGoodAnodesL3!=0.) normNrecpLay3 = nrecpLay3/(nGoodAnodesL3*nEventRP);
	if(nGoodAnodesL4!=0.) normNrecpLay4 = nrecpLay4/(nGoodAnodesL4*nEventRP);
	if(fHistoCalibration){
	  if (harr[kSDDRecpModPattern+offsRP] && harr[kSDDRecpModPatternNorm+offsRP]){
	    ((TH1*)harr[kSDDRecpModPatternNorm+offsRP])->Divide( (TH1*)harr[kSDDRecpModPattern+offsRP],((TH1F*)(fHistoCalibration->At(0))),1.,nEventRP);
	  }
	  if (harr[kSDDRecpLadModLay3+offsRP] && harr[kSDDRecpLadModLay3Norm+offsRP]){
	    ((TH2*)harr[kSDDRecpLadModLay3Norm+offsRP])->Divide( (TH2*)harr[kSDDRecpLadModLay3+offsRP],((TH2F*)(fHistoCalibration->At(1))),1.,nEventRP);
	  }
	  if (harr[kSDDRecpLadModLay4+offsRP] && harr[kSDDRecpLadModLay4Norm+offsRP]){
	    ((TH2*)harr[kSDDRecpLadModLay4Norm+offsRP])->Divide( (TH2*)harr[kSDDRecpLadModLay4+offsRP],((TH2F*)(fHistoCalibration->At(2))),1.,nEventRP);
	  }
	}
      }
      if (htmp27) {
	htmp27->SetBinContent(1,nEventRP);
	htmp27->SetBinContent(2,normNrecpMod);
	htmp27->SetBinContent(3,normNrecpLay3);
	htmp27->SetBinContent(4,normNrecpLay4);
      }

      FillRelativeOccupancyHistos((TH2*)harr[kSDDRecpLadModLay3Norm+offsRP],(TH1*)harr[kSDDRecpRelOccLay3+offsRP]);
      FillRelativeOccupancyHistos((TH2*)harr[kSDDRecpLadModLay4Norm+offsRP],(TH1*)harr[kSDDRecpRelOccLay4+offsRP]);
      
      // RecPoints 2 Raws Ratio
      if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kTRUE) {
	TH2* hratio3  = (TH2*)harr[kSDDRecpToRawLay3+offsRP];
	TH2* hrp3  = (TH2*)harr[kSDDRecpLadModLay3+offsRP];
	TH2* hrw3 =  (TH2*)fAliITSQADataMakerRec->GetRawsData(kSDDRawLadModLay3 + offsRW,trCl);
	FillRecToRaw(hrp3,hrw3,hratio3);
	TH2* hratio4  = (TH2*)harr[kSDDRecpToRawLay4+offsRP];
	TH2* hrp4  = (TH2*)harr[kSDDRecpLadModLay4+offsRP];
	TH2* hrw4 =  (TH2*)fAliITSQADataMakerRec->GetRawsData(kSDDRawLadModLay4 + offsRW,trCl);
	FillRecToRaw(hrp4,hrw4,hratio4);
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

  int offsRW = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  fSDDhRawsTask = 0;
  if(fkOnline){AliInfo("Book Online Histograms for SDD\n");}
  else {AliInfo("Book Offline Histograms for SDD\n ");}

  TPaveText *paveText0=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
  paveText0->AddText("");
  paveText0->SetFillColor(kMagenta+2);
  paveText0->SetTextColor(kWhite);
  paveText0->SetBorderSize(1);
  paveText0->SetLineWidth(1);	
  
  TH1D *h0 = new TH1D("SDDModPattern","HW Modules pattern",fgknSDDmodules,239.5,499.5); //0
  h0->GetXaxis()->SetTitle("Module Number");
  h0->GetYaxis()->SetTitle("Counts");
  h0->SetOption("bar1");
  h0->SetBarOffset(0.01);
  h0->SetBarWidth(0.95);
  h0->SetFillColor(45);
  h0->GetListOfFunctions()->Add(paveText0);
	
  rv = fAliITSQADataMakerRec->Add2RawsList(h0,kSDDRawModPattern+offsRW, expert, !image, !saveCorr);
  fSDDhRawsTask++;
        
  TPaveText *paveText1=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
  paveText1->AddText("");
  paveText1->SetFillColor(kMagenta+2);
  paveText1->SetTextColor(kWhite);
  paveText1->SetBorderSize(1);
  paveText1->SetLineWidth(1);

  //zPhi distribution using ladder and modules numbers
  TH2D *hphil3 = new TH2D("SDDphizL3","SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);//1
  hphil3->GetXaxis()->SetTitle("z[Module Number L3 ]");
  hphil3->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  hphil3->SetStats(0);
  hphil3->GetListOfFunctions()->Add(paveText1);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil3,kSDDRawLadModLay3+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;
 
  TPaveText *paveText2=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
  paveText2->AddText("");
  paveText2->SetFillColor(kMagenta+2);
  paveText2->SetTextColor(kWhite);
  paveText2->SetBorderSize(1);
  paveText2->SetLineWidth(1);

  TH2D *hphil4 = new TH2D("SDDphizL4","SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5); //2
  hphil4->GetXaxis()->SetTitle("z[Module Number L4]");
  hphil4->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
  hphil4->SetStats(0);
  hphil4->GetListOfFunctions()->Add(paveText2);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil4,kSDDRawLadModLay4+offsRW, !expert, image, !saveCorr); 
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
  rv = fAliITSQADataMakerRec->Add2RawsList(h0norm,kSDDRawModPatternNorm+offsRW, expert, !image, !saveCorr);
  fSDDhRawsTask++;
  
  //zPhi distribution using ladder and modules numbers
  TH2F *hphil3norm = new TH2F("SDDphizL3NORM","NORM SDD #varphiz Layer3 ",12,0.5,6.5,14,0.5,14.5);//4
  hphil3norm->GetXaxis()->SetTitle("z[Module Number L3 ]");
  hphil3norm->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  hphil3norm->SetStats(0);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil3norm,kSDDRawLadModLay3Norm+offsRW, expert, !image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2F *hphil4norm = new TH2F("SDDphizL4NORM","NORM SDD #varphiz Layer4 ",16,0.5,8.5,22,0.5,22.5); //5
  hphil4norm->GetXaxis()->SetTitle("z[Module Number L4]");
  hphil4norm->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
  hphil4norm->SetStats(0);
  rv = fAliITSQADataMakerRec->Add2RawsList(hphil4norm,kSDDRawLadModLay4Norm+offsRW, expert, !image, !saveCorr); 
  fSDDhRawsTask++;
	
  // Histos with number of digits distributions (filled once per event)
  Float_t digbins[201];
  digbins[0]=-0.5;
  for(Int_t ib=1; ib<=200; ib++) digbins[ib]=digbins[ib-1]+(1<<(ib/15));
  Float_t modBins[fgknSDDmodules+1];
  for(Int_t ib=0; ib<=fgknSDDmodules; ib++) modBins[ib]=239.5+ib;
  
  TH1F *hTotDigits = new TH1F("SDDnOfDigits","SDD Total number of digits ; # of digits",200,digbins);
  rv = fAliITSQADataMakerRec->Add2RawsList(hTotDigits,kSDDNofDigits+offsRW, expert, !image, !saveCorr); //6  
  fSDDhRawsTask++;
  
  TH2F *hDigitsPerModule = new TH2F("SDDnOfDigitsPerModule","SDD Digits vs module ; module ID ; # of digits", fgknSDDmodules,modBins,200,digbins);
  rv = fAliITSQADataMakerRec->Add2RawsList(hDigitsPerModule,kSDDNofDigitsVsMod+offsRW, expert, !image, !saveCorr); //7   
  fSDDhRawsTask++;
  
  // active modules and drift regions
  TH2F *hcalibl3 = new TH2F("SDDphizCalibL3","SDDCalibL3 ",12,0.5,6.5,14,0.5,14.5);//8
  hcalibl3->GetXaxis()->SetTitle("z[Module Number L3]");
  hcalibl3->GetYaxis()->SetTitle("#varphi[ Ladder Number L3]");
  hcalibl3->SetStats(0);
  hcalibl3->SetMaximum(2);
  rv = fAliITSQADataMakerRec->Add2RawsList(hcalibl3,kActiveModLay3+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;
  
  TH2F *hcalibl4 = new TH2F("SDDphizCalibL4","SDDCalibL4 ",16,0.5,8.5,22,0.5,22.5); //9
  hcalibl4->GetXaxis()->SetTitle("z[Module Number L4]");
  hcalibl4->GetYaxis()->SetTitle("#varphi[Ladder Number L4]");
  hcalibl4->SetStats(0);
  hcalibl4->SetMaximum(2);
  rv = fAliITSQADataMakerRec->Add2RawsList(hcalibl4,kActiveModLay4+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++;

  TH1F *hsummarydata = new TH1F("SDDRawDataCheck","SDDRawDataCheck",46,-0.5,45.5);//10 summary of raw data checks Non expert and image
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

  rv = fAliITSQADataMakerRec->Add2RawsList(hsummarydata,kSDDRawDataCheck+offsRW, !expert, image, !saveCorr); 
  fSDDhRawsTask++; 
  fOnlineOffsetRaws = fSDDhRawsTask;

  //online part
  if(fkOnline){
    //DDL Pattern 
    TH2F *hddl = new TH2F("SDDDDLPattern","SDD DDL Pattern ",24,-0.5,11.5,24,-0.5,23.5); //11
    hddl->GetXaxis()->SetTitle("Channel");
    hddl->GetYaxis()->SetTitle("DDL Number");
    rv = fAliITSQADataMakerRec->Add2RawsList(hddl,kSDDOnlineDDLPattern+offsRW, expert, !image, !saveCorr);
    fSDDhRawsTask++;

    //Event Size 
    TH1F *hsize = new TH1F("SDDEventSize","SDD Event Size ",500,-0.5,199.5);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
    hsize->SetBit(TH1::kCanRebin);
#endif
    hsize->GetXaxis()->SetTitle("Event Size [kB]");
    hsize->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RawsList(hsize,kSDDDataSize+offsRW, expert, !image, !saveCorr); 
    fSDDhRawsTask++;

    Int_t nTimeBins= 256/fTimeBinSize;
    Int_t index1 = 0;

    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0;iside<fgknSide;iside++){
	AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	TProfile2D *fModuleChargeMapFSE = new TProfile2D(Form("SDDchargeMapFSE_L%d_%d_%d_%d",lay,lad,det,iside),Form("SDDChargeMapForSingleEvent_L%d_%d_%d_%d",lay,lad,det,iside),nTimeBins,-0.5,255.5,256,-0.5,255.5);
	fModuleChargeMapFSE->GetXaxis()->SetTitle("Time Bin");
	fModuleChargeMapFSE->GetYaxis()->SetTitle("Anode");
	rv = fAliITSQADataMakerRec->Add2RawsList(fModuleChargeMapFSE,kChargeMapFirstMod + index1 + offsRW, expert, !image, !saveCorr);	  
	index1++;
      }
    }
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0;iside<fgknSide;iside++){
	AliITSgeomTGeo::GetModuleId(moduleSDD+fgkmodoffset, lay, lad, det);
	TProfile2D *fModuleChargeMap = new TProfile2D(Form("SDDchargeMap_L%d_%d_%d_%d",lay,lad,det,iside),Form("SDDChargeMap_L%d_%d_%d_%d",lay,lad,det,iside),nTimeBins,-0.5,255.5,256,-0.5,255.5);
	fModuleChargeMap->GetXaxis()->SetTitle("Time Bin");
	fModuleChargeMap->GetYaxis()->SetTitle("Anode Number");
	rv = fAliITSQADataMakerRec->Add2RawsList(fModuleChargeMap,kChargeMapFirstMod + index1 + offsRW, expert, !image, !saveCorr); 
	index1++;
      }
    }
    fSDDhRawsTask+=index1;
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

  if(!fDDLModuleMap) CreateTheMap();
  if(rawReader->GetType() != 7) return rv;  // skips non physical triggers
  AliDebug(AliQAv1::GetQADebugLevel(),"entering MakeRaws\n");                 
  rawReader->Reset();       
  rawReader->Select("ITSSDD");
  AliITSRawStream *stream=AliITSRawStreamSDD::CreateRawStreamSDD(rawReader);
  stream->SetDDLModuleMap(fDDLModuleMap);
  
  Int_t index=0;
  if(fkOnline) {
    for(Int_t moduleSDD =0; moduleSDD<fgknSDDmodules; moduleSDD++){
      for(Int_t iside=0;iside<fgknSide;iside++) {
	fAliITSQADataMakerRec->ResetRawsData(kChargeMapFirstMod+index+offsRW);   
	index++;
      }
    }
  }
  
  Int_t lay, lad, det;   
  Int_t cnt = 0;
  Int_t cntMod[fgknSDDmodules];
  for(Int_t jmod=0; jmod<fgknSDDmodules; jmod++) cntMod[jmod]=0;
  
  Int_t iddl = -1;
  Int_t isddmod = -1;
  Int_t coord1, coord2, signal, moduleSDD; 
  Int_t prevDDLID = -1;
  UInt_t size = 0;
  Int_t totalddl=static_cast<int>(fDDLModuleMap->GetNDDLs());
  Bool_t *ddldata=new Bool_t[totalddl];
  for(Int_t jddl=0;jddl<totalddl;jddl++) ddldata[jddl]=kFALSE;

  while(stream->Next()) {
    iddl = rawReader->GetDDLID();// - fgkDDLIDshift;
    if(iddl<0)isddmod=-1;
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
    fAliITSQADataMakerRec->FillRawsData( kSDDRawModPattern + offsRW,isddmod);   
    if(lay==3) fAliITSQADataMakerRec->FillRawsData(kSDDRawLadModLay3+offsRW,det+0.5*iside-0.25,lad); //phiz l3 not norm
    else if(lay==4) fAliITSQADataMakerRec->FillRawsData(kSDDRawLadModLay4+offsRW,det+0.5*iside-0.25,lad); //phiz l4 not norm
 
    if(fkOnline) {

      fAliITSQADataMakerRec->FillRawsData(kSDDOnlineDDLPattern+offsRW, (stream->GetCarlosId())+0.5*iside -0.5,iddl);
      Int_t index1 = moduleSDD * 2 + iside;
      
      if(index1<0){
        AliDebug(AliQAv1::GetQADebugLevel(),Form("Wrong index number %d - patched to 0\n",index1));
	index1 = 0;
      }      

      fAliITSQADataMakerRec->FillRawsData(kChargeMapFirstMod + index1 +offsRW, coord2, coord1, signal);     
      fAliITSQADataMakerRec->FillRawsData(kChargeMapFirstMod + index1 + fgknSDDmodules*fgknSide +offsRW, coord2, coord1, signal); 

    }//online
    cnt++;
    if(moduleSDD>=0 && moduleSDD<fgknSDDmodules) cntMod[moduleSDD]++;
    if(!(cnt%10000)) AliDebug(AliQAv1::GetQADebugLevel(),Form(" %d raw digits read",cnt));
  }//end next()
  fAliITSQADataMakerRec->FillRawsData(kSDDNofDigits+offsRW,cnt);
  for(Int_t jmod=0; jmod<fgknSDDmodules; jmod++){
    fAliITSQADataMakerRec->FillRawsData(kSDDNofDigitsVsMod+offsRW,jmod+fgkmodoffset,cntMod[jmod]);
  }
  if(fkOnline){
    fAliITSQADataMakerRec->FillRawsData(kSDDDataSize+offsRW,size/1024.);//KB
  }
	
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Event completed, %d raw digits read",cnt)); 
  delete stream;
  delete [] ddldata;
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
  fSDDhDigitsTask=0;
  int offsD = fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()];

  TH1F* h0=new TH1F("SDD DIGITS Module Pattern","SDD DIGITS Module Pattern",260,239.5,499.5);       //hmod
  h0->GetXaxis()->SetTitle("SDD Module Number");
  h0->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h0,kSDDDigModPattern+offsD, !expert, image);
  fSDDhDigitsTask++;

  TH1F* h1=new TH1F("SDD Anode Distribution","DIGITS Anode Distribution",512,-0.5,511.5);      //hanocc
  h1->GetXaxis()->SetTitle("Anode Number");
  h1->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h1,kSDDDigAnode+offsD, !expert, image);
  fSDDhDigitsTask++;

  TH1F* h2=new TH1F("SDD Tbin Distribution","DIGITS Tbin Distribution",256,-0.5,255.5);      //htbocc
  h2->GetXaxis()->SetTitle("Tbin Number");
  h2->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h2,kSDDDigTimeBin+offsD, !expert, image);
  fSDDhDigitsTask++;

  TH1F* h3=new TH1F("SDD ADC Counts Distribution","DIGITS ADC Counts Distribution",200,0.,1024.);          //hsig
  h3->GetXaxis()->SetTitle("ADC Value");
  h3->GetYaxis()->SetTitle("# DIGITS");
  rv = fAliITSQADataMakerRec->Add2DigitsList(h3,kSDDDigADC+offsD, !expert, image);
  fSDDhDigitsTask++;

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SDD Digits histograms booked\n",fSDDhDigitsTask));
  //
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASDDDataMakerRec::MakeDigits(TTree * digits)
{ 

  // Fill QA for DIGIT - SDD -

  Int_t rv = 0 ; 
  int offsD = fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  TBranch *branchD = digits->GetBranch("ITSDigitsSDD");

  if (!branchD) {
    AliError("can't get the branch with the ITS SDD digits !");
    return rv ;
  }
  
  static TClonesArray statDigits("AliITSdigitSDD");
  TClonesArray *iITSdigits = &statDigits;
  branchD->SetAddress(&iITSdigits);

  for(Int_t i=0; i<260; i++){
    Int_t nmod=i+240;
    digits->GetEvent(nmod);
    Int_t ndigits = iITSdigits->GetEntries();
    fAliITSQADataMakerRec->FillDigitsData(kSDDDigModPattern+offsD,nmod,ndigits);

    for (Int_t idig=0; idig<ndigits; idig++) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t iz=dig->GetCoord1();  // cell number z
      Int_t ix=dig->GetCoord2();  // cell number x
      Int_t sig=dig->GetSignal();
      fAliITSQADataMakerRec->FillDigitsData(kSDDDigAnode+offsD,iz);
      fAliITSQADataMakerRec->FillDigitsData(kSDDDigTimeBin+offsD,ix);
      fAliITSQADataMakerRec->FillDigitsData(kSDDDigADC+offsD,sig);
    }
  }

  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASDDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SDD -

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];

  TH1F *h0 = new TH1F("SDDLay3TotCh","Layer 3 total charge",250,-0.5, 499.5); //position number 0
  //h0->SetBit(TH1::kCanRebin);
  h0->GetXaxis()->SetTitle("Charge (keV)");
  h0->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h0, kSDDRecpChargeLay3+offsRP, !expert, image);//NON expert image
  fSDDhRecPointsTask++;
 
  TH1F *h1 = new TH1F("SDDLay4TotCh","Layer 4 total charge",250,-0.5, 499.5);
  //h1->SetBit(TH1::kCanRebin);
  h1->GetXaxis()->SetTitle("Charge (keV)");
  h1->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h1, kSDDRecpChargeLay4+offsRP, !expert, image);
  fSDDhRecPointsTask++;

  TH2F *h2 = new TH2F("SDDGlobalCoordDistribYX","YX Global Coord Distrib",56,-28,28,56,-28,28);
  h2->GetYaxis()->SetTitle("Y (cm)");
  h2->GetXaxis()->SetTitle("X (cm)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h2,kSDDRecpGloXY+offsRP, expert, image);
  fSDDhRecPointsTask++;

  TH2F *h3 = new TH2F("SDDGlobalCoordDistribRZ","RZ Global Coord Distrib",128,-32,32,56,12,26);
  h3->GetYaxis()->SetTitle("R (cm)");
  h3->GetXaxis()->SetTitle("Z (cm)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h3,kSDDRecpGloRZ+offsRP, expert, image);
  fSDDhRecPointsTask++;
  
  TH2F *h4 = new TH2F("SDDGlobalCoordDistribL3PHIZ","#varphi Z Global Coord Distrib L3",96,-23,23,112,-TMath::Pi(),TMath::Pi());
  h4->GetYaxis()->SetTitle("#varphi (rad)");
  h4->GetXaxis()->SetTitle("Z (cm)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h4,kSDDRecpPhiZLay3+offsRP, expert, image);
  fSDDhRecPointsTask++;

  TH2F *h5 = new TH2F("SDDGlobalCoordDistribL4PHIZ","#varphi Z Global Coord Distrib L4",128,-31,31,176,-TMath::Pi(),TMath::Pi());
  h5->GetYaxis()->SetTitle("#varphi (rad)");
  h5->GetXaxis()->SetTitle("Z (cm)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h5,kSDDRecpPhiZLay4+offsRP, expert, image);
  fSDDhRecPointsTask++;
  
  TH1F *h6 = new TH1F("SDDModPatternRP","Modules pattern RP",fgknSDDmodules,239.5,499.5);
  h6->GetXaxis()->SetTitle("Module number");
  h6->GetYaxis()->SetTitle("Entries");
  h6->SetOption("bar1");
  h6->SetBarOffset(0.01);
  h6->SetBarWidth(0.95);
  h6->SetFillColor(39);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h6,kSDDRecpModPattern+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  TPaveText *paveText7=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
  paveText7->AddText("");
  paveText7->SetFillColor(kMagenta+2);
  paveText7->SetTextColor(kWhite);
  paveText7->SetBorderSize(1);
  paveText7->SetLineWidth(1);	
  	
  TH2F *h7 = new TH2F("SDDModPatternL3RP","Modules pattern L3 RP",12,0.5,6.5,14,0.5,14.5);
  h7->GetXaxis()->SetTitle("z[#Module L3 ]");
  h7->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  h7->GetListOfFunctions()->Add(paveText7);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h7,kSDDRecpLadModLay3+offsRP, !expert, image);
  fSDDhRecPointsTask++;

  TPaveText *paveText8=new TPaveText(0.3,0.88,0.9,0.99,"NDC");
  paveText8->AddText("");
  paveText8->SetFillColor(kMagenta+2);
  paveText8->SetTextColor(kWhite);
  paveText8->SetBorderSize(1);
  paveText8->SetLineWidth(1);	
  
  TH2F *h8 = new TH2F("SDDModPatternL4RP","Modules pattern L4 RP",16,0.5,8.5,22,0.5,22.5); 
  h8->GetXaxis()->SetTitle("[#Module L3 ]");
  h8->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  h8->GetListOfFunctions()->Add(paveText8);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h8,kSDDRecpLadModLay4+offsRP, !expert, image);
  fSDDhRecPointsTask++;

  //------------------------norm--------------------------//


  TH1F *h9 = new TH1F("SDDModPatternRPNORM","Modules pattern RP NORM",fgknSDDmodules,239.5,499.5); 
  h9->GetXaxis()->SetTitle("Module number"); 
  h9->GetYaxis()->SetTitle("Entries");
  h9->SetOption("bar1");
  h9->SetBarOffset(0.01);
  h9->SetBarWidth(0.95);
  h9->SetFillColor(49);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h9,kSDDRecpModPatternNorm+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  
  TH2F *h10 = new TH2F("SDDModPatternL3RPNORM","Modules pattern L3 RP NORM",12,0.5,6.5,14,0.5,14.5);
  h10->GetXaxis()->SetTitle("z[#Module L3 ]");
  h10->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h10,kSDDRecpLadModLay3Norm+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  TH2F *h11 = new TH2F("SDDModPatternL4RPNORM","Modules pattern L4 RP NORM",16,0.5,8.5,22,0.5,22.5);
  h11->GetXaxis()->SetTitle("[#Module L4 ]");
  h10->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h11,kSDDRecpLadModLay4Norm+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  //--------------------------------------------------------//

  TH2F *h12 = new TH2F("SDDLocalCoordDistrib","Local Coord Distrib",160,-4,4,160,-4,4);
  h12->GetXaxis()->SetTitle("X local coord, drift (cm)");
  h12->GetYaxis()->SetTitle("Z local coord, anode (cm)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h12,kSDDRecpLocalCoord+offsRP, expert, !image);
  fSDDhRecPointsTask++;
  
  TH1F *h13 = new TH1F("SDDrdistrib_Layer3" ,"SDD r distribution Layer3" ,100,14.,16.5);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  h13->SetBit(TH1::kCanRebin);
#endif
  h13->GetXaxis()->SetTitle("r (cm)");
  h13->GetXaxis()->CenterTitle();
  h13->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h13,kSDDRecpRLay3+offsRP, expert, !image);
  fSDDhRecPointsTask++;
  
  TH1F *h14 = new TH1F("SDDrdistrib_Layer4" ,"SDD r distribution Layer4" ,100,23.,25.);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  h14->SetBit(TH1::kCanRebin);
#endif
  h14->GetXaxis()->SetTitle("r (cm)");
  h14->GetXaxis()->CenterTitle();
  h14->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h14,kSDDRecpRLay4+offsRP, expert, !image);
  fSDDhRecPointsTask++;
  
  TH1F *h15 = new TH1F("SDDphidistrib_Layer3","SDDphidistrib_Layer3" ,180,-TMath::Pi(),TMath::Pi());
  h15->GetXaxis()->SetTitle("#varphi (rad)");
  h15->GetXaxis()->CenterTitle();
  h15->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h15,kSDDRecpPhiLay3+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  TH1F *h16 = new TH1F("SDDphidistrib_Layer4","SDDphidistrib_Layer4" ,180,-TMath::Pi(),TMath::Pi());
  h16->GetXaxis()->SetTitle("#varphi (rad)");
  h16->GetXaxis()->CenterTitle();
  h16->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h16,kSDDRecpPhiLay4+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  
  TH1F *h17 = new TH1F("SDDdrifttime_Layer3","SDDdrifttime_Layer3",45,-0.5,4499.5);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  h17->SetBit(TH1::kCanRebin);
#endif
  h17->GetXaxis()->SetTitle("Drift time (ns)");
  h17->GetXaxis()->CenterTitle();
  h17->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h17,kSDDRecpDriftTimeLay3+offsRP, !expert, image);
  fSDDhRecPointsTask++;
  
  TH1F *h18 = new TH1F("SDDdrifttime_Layer4","SDDdrifttime_Layer4",45,-0.5,4499.5);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,0)
  h18->SetBit(TH1::kCanRebin);
#endif
  h18->GetXaxis()->SetTitle("Drift time (ns)");
  h18->GetXaxis()->CenterTitle();
  h18->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h18,kSDDRecpDriftTimeLay4+offsRP, !expert, image);
  fSDDhRecPointsTask++;
  
  TH1F *h19 = new TH1F("SDDL3_RelativeOccupancy","Layer 3 Relative Occupancy (RecPoints)",200,0.,0.2);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h19,kSDDRecpRelOccLay3+offsRP, expert, !image); 
  fSDDhRecPointsTask++;
  
  TH1F *h20 = new TH1F("SDDL4_RelativeOccupancy","Layer 4 Relative Occupancy (RecPoints)",200,0.,0.2);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h20,kSDDRecpRelOccLay4+offsRP, expert, !image);
  fSDDhRecPointsTask++;
		
  TH2F *h21 = new TH2F("SDDL3_Rec2Raw_2D","L3 RecPoints to Raws 2D",12,0.5,6.5,14,0.5,14.5);  
  h21->GetXaxis()->SetTitle("z[#Module L3 ]");
  h21->GetYaxis()->SetTitle("#varphi[#Ladder L3]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h21,kSDDRecpToRawLay3+offsRP, expert, !image);
  fSDDhRecPointsTask++;
  
  TH2F *h22 = new TH2F("SDDL4_Rec2Raw_2D","L4 RecPoints to Raws 2D",16,0.5,8.5,22,0.5,22.5);
  h22->GetXaxis()->SetTitle("[#Module L4 ]");
  h22->GetYaxis()->SetTitle("#varphi[#Ladder L4]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h22,kSDDRecpToRawLay4+offsRP, expert, !image);
  fSDDhRecPointsTask++;

  TH2F *h23 = new TH2F("SDDL3_clusizAnode","L3 CluSiz Anode",50,0.,6400.,21,-0.5,20.5);
  h23->GetXaxis()->SetTitle("Drift time (ns)");
  h23->GetYaxis()->SetTitle("Cluster Size (anodes)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h23,kSDDRecpCluSizAnLay3+offsRP, expert, !image); // 25
  fSDDhRecPointsTask++;
  
  TH2F *h24 = new TH2F("SDDL4_clusizAnode","L4 CluSiz Anode",50,0.,6400.,21,-0.5,20.5);
  h24->GetXaxis()->SetTitle("Drift time (ns)");
  h24->GetYaxis()->SetTitle("Cluster Size (anodes)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h24,kSDDRecpCluSizAnLay4+offsRP, expert, !image); // 26
  fSDDhRecPointsTask++;

  TH2F *h25 = new TH2F("SDDL3_clusizTB","L3 CluSiz TimeBin",50,0.,6400.,21,-0.5,20.5);
  h25->GetXaxis()->SetTitle("Drift time (ns)");
  h25->GetYaxis()->SetTitle("Cluster Size (time bins)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h25,kSDDRecpCluSizTbLay3+offsRP, expert, !image); // 25
  fSDDhRecPointsTask++;
  
  TH2F *h26 = new TH2F("SDDL4_clusizTB","L4 CluSiz TimeBin",50,0.,6400.,21,-0.5,20.5);
  h26->GetXaxis()->SetTitle("Drift time (ns)");
  h26->GetYaxis()->SetTitle("Cluster Size (time bins)");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(h26,kSDDRecpCluSizTbLay4+offsRP, expert, !image); // 26
  fSDDhRecPointsTask++;
	

  TH1F *hsummarydatarp = new TH1F("SDDRecPointCheck","SDDRecPointCheck",46,-0.5,45.5);//27 summary of recpoint checks
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

  hsummarydatarp->SetOption("text bar1");
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
  hsummarydatarp->GetXaxis()->SetBinLabel(31,"m_overthl3");
  hsummarydatarp->GetXaxis()->SetBinLabel(32,"dr_overthl3");

  //l4

  hsummarydatarp->GetXaxis()->SetBinLabel(33,"m_actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(34,"m_fill4");
  hsummarydatarp->GetXaxis()->SetBinLabel(35,"dr_actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(36,"dr_fill4");
  hsummarydatarp->GetXaxis()->SetBinLabel(37,"m_excl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(38,"m_empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(39,"dr_excl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(40,"dr_empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(41,"m_exc-actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(42,"m_act-empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(43,"dr_exc-actl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(44,"dr_act-empl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(45,"m_overthl4");
  hsummarydatarp->GetXaxis()->SetBinLabel(46,"dr_overthl4");

  hsummarydatarp->GetXaxis()->LabelsOption("v");

  rv = fAliITSQADataMakerRec->Add2RecPointsList(hsummarydatarp,kSDDRecpDataCheck+offsRP, !expert, image);
  fSDDhRecPointsTask++;

  fOnlineOffsetRecPoints = fSDDhRecPointsTask;
  if(fkOnline){
    TH2F *hxy1ev = new TH2F("SDDGlobalCoordDistribYXFSE","YX Global Coord Distrib FSE",112,-28,28,112,-28,28);
    hxy1ev->GetYaxis()->SetTitle("Y (cm)");
    hxy1ev->GetXaxis()->SetTitle("X (cm)");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hxy1ev,kSDDRecpOnlineGloXYSingEv+offsRP, expert, !image);
    fSDDhRecPointsTask++;
    
    TH2F *hrz1ev = new TH2F("SDDGlobalCoordDistribRZFSE","RZ Global Coord Distrib FSE",128,-32,32,56,12,26);
    hrz1ev->GetYaxis()->SetTitle("R (cm)");
    hrz1ev->GetXaxis()->SetTitle("Z (cm)");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hrz1ev,kSDDRecpOnlineGloRZSingEv+offsRP, expert, !image);
    fSDDhRecPointsTask++;      
  }
  
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

  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  TClonesArray *recpoints=NULL; 
  if(fkOnline){
    recpoints = rpcont->FetchClusters(0,clustersTree,fAliITSQADataMakerRec->GetEventNumber());
  }else{
    recpoints = rpcont->FetchClusters(0,clustersTree);
  }
  AliDebug(10,Form("Fetched RecPoints for %d SDD modules",recpoints->GetEntriesFast()));
  if(!rpcont->GetStatusOK() || !rpcont->IsSDDActive()){
    AliError("can't get SDD clusters !");
    return rv;
  }
  int offsRP = fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  if(fkOnline){
    fAliITSQADataMakerRec->ResetRecPointsData(kSDDRecpOnlineGloXYSingEv+offsRP);
    fAliITSQADataMakerRec->ResetRecPointsData(kSDDRecpOnlineGloRZSingEv+offsRP);
  }

  Int_t firMod=TMath::Max(0,AliITSgeomTGeo::GetModuleIndex(3,1,1));
  Int_t lasMod=AliITSgeomTGeo::GetModuleIndex(5,1,1);
  for(Int_t module=firMod; module<lasMod;module++){
    recpoints = rpcont->UncheckedGetClusters(module);
    npoints += recpoints->GetEntries();
    for(Int_t j=0;j<recpoints->GetEntries();j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j); 
      Int_t index = recp->GetDetectorIndex();
      lay=recp->GetLayer();
      if(lay < 2 || lay > 3) continue;
      Int_t modnumb=index+AliITSgeomTGeo::GetModuleIndex(lay+1,1,1);
      AliITSgeomTGeo::GetModuleId(modnumb, lay, lad, det);  
      fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpModPattern+offsRP,modnumb);
      recp->GetGlobalXYZ(cluglo);
      Float_t rad=TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi=TMath::ATan2(cluglo[1],cluglo[0]);
      Float_t drifttime=recp->GetDriftTime();
      fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpLocalCoord+offsRP,recp->GetDetLocalX(),recp->GetDetLocalZ());
      fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpGloXY+offsRP,cluglo[0],cluglo[1]);
      fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpGloRZ+offsRP,cluglo[2],rad);
      if(fkOnline) {
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpOnlineGloXYSingEv+offsRP,cluglo[0],cluglo[1]);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpOnlineGloRZSingEv+offsRP,cluglo[2],rad);
      }
      Int_t iside=recp->GetDriftSide();
      lay=recp->GetLayer();
      if(lay == 2) {
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpChargeLay3+offsRP, recp->GetQ());
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpLadModLay3+offsRP, det+0.5*iside-0.5,lad);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpRLay3+offsRP, rad);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpPhiLay3+offsRP, phi);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpPhiZLay3+offsRP, cluglo[2],phi);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpDriftTimeLay3+offsRP, drifttime);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpCluSizAnLay3+offsRP, drifttime,recp->GetNz());
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpCluSizTbLay3+offsRP, drifttime,recp->GetNy());
      } else if(lay == 3) {
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpChargeLay4+offsRP, recp->GetQ());
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpLadModLay4+offsRP, det+0.5*iside-0.5,lad);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpRLay4+offsRP, rad);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpPhiLay4+offsRP, phi);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpPhiZLay4+offsRP, cluglo[2],phi);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpDriftTimeLay4+offsRP, drifttime);
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpCluSizAnLay4+offsRP, drifttime, recp->GetNz());
	fAliITSQADataMakerRec->FillRecPointsData(kSDDRecpCluSizTbLay4+offsRP, drifttime, recp->GetNy());
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
  if(!calibSDD){
    AliError("Calibration object retrieval failed! SDD will not be processed");
    fCalibration = NULL;
  }else{
    fCalibration = (TObjArray *)calibSDD->GetObject();      
    if(!cacheStatus) calibSDD->SetObject(NULL);
    calibSDD->SetOwner(kTRUE);
    if(!cacheStatus) delete calibSDD;

    AliCDBId cdbid=calibSDD->GetId();
    fPulserRun=cdbid.GetFirstRun();

    AliITSCalibrationSDD * cal=NULL;
    for(Int_t imod=0;imod<fgknSDDmodules;imod++){
      Int_t module=imod + 240;
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(module,lay,lad,det);
      Int_t index=1+(det-1)*2;
      Int_t nGoodAnPerSide[2]={0,0};
      cal=(AliITSCalibrationSDD*)fCalibration->At(imod);
      if(cal && !cal->IsBad()){
	for(Int_t iAn=0; iAn<512; iAn++){
	  Int_t iwing=cal->GetWing(iAn);
	  Int_t ichip=cal->GetChip(iAn);
	  if(!cal->IsChipBad(ichip) && !cal->IsBadChannel(iAn)) nGoodAnPerSide[iwing]++;	    
	}
      }
      Int_t totGoodAn=nGoodAnPerSide[0]+nGoodAnPerSide[1];
      ((TH1F*)(fHistoCalibration->At(0)))->SetBinContent(imod+1,totGoodAn);
      ((TH2F*)(fHistoCalibration->At(lay-2)))->SetBinContent(index,lad,nGoodAnPerSide[0]);
      ((TH2F*)(fHistoCalibration->At(lay-2)))->SetBinContent(index+1,lad,nGoodAnPerSide[1]);
    }
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
//____________________________________________________________________
void AliITSQASDDDataMakerRec::FillRelativeOccupancyHistos(TH2* hnormc, TH1* hrelocc) const{
  // fill histograms with relative occupancy
  if (hnormc && hrelocc){
    for(Int_t i=0; i<hnormc->GetNbinsX(); i++){
      for(Int_t j=0; j<hnormc->GetNbinsY(); j++){
	hrelocc->Fill(hnormc->GetBinContent(i,j));
      }
    }
  }
}
//____________________________________________________________________
void AliITSQASDDDataMakerRec::FillRecToRaw(TH2* hrpLay,TH2* hrwLay,TH2* hratioLay) const{
  // fill histogram with rec to raw ratios
  if (hrpLay && hrwLay && hratioLay){
    Int_t nbx1=hrpLay->GetNbinsX();
    Int_t nbx2=hrwLay->GetNbinsX();
    Int_t nbx3=hratioLay->GetNbinsX();
    Int_t nby1=hrpLay->GetNbinsY();
    Int_t nby2=hrwLay->GetNbinsY();
    Int_t nby3=hratioLay->GetNbinsY();
    if(nbx1==nbx2 && nbx1==nbx3 && nby1==nby2 && nby1==nby3){
      hratioLay->Divide(hrpLay,hrwLay);
    }
    else AliWarning("Number of bins for Raws and RecPoints do not match\n");
  }
}
