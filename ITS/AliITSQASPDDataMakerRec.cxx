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
/* $Id$  */
//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino
//  M. Nicassio D. Elia INFN Bari March 2008
//  maria.nicassio@ba.infn.it
    

// --- ROOT system ---
#include <TTree.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASPDDataMakerRec.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSPDErrorLog.h"
#include "AliITSdigitSPD.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"

ClassImp(AliITSQASPDDataMakerRec)

//____________________________________________________________________________
AliITSQASPDDataMakerRec::AliITSQASPDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc, AliITSRawStreamSPDErrorLog *aliITSRawStreamSPDErrorLog) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSPDhRawsTask(0),
fSPDhDigitsTask(0),
fSPDhRecPointsTask(0),
fGenRawsOffset(0),
fGenDigitsOffset(0),
fGenRecPointsOffset(0),
fAdvLogger(aliITSRawStreamSPDErrorLog) 
{
  //ctor used to discriminate OnLine-Offline analysis  
  //AliInfo(Form("AliRecoParam::kNSpecies %d\n",AliRecoParam::kNSpecies));
	fGenRawsOffset = new Int_t[AliRecoParam::kNSpecies];
	fGenRecPointsOffset = new Int_t[AliRecoParam::kNSpecies];
	fGenDigitsOffset = new Int_t[AliRecoParam::kNSpecies];
	for(Int_t i=0; i<AliRecoParam::kNSpecies;i++) {
		fGenRawsOffset[i] = 0;
		fGenRecPointsOffset[i] = 0;
		fGenDigitsOffset[i]=0;
	}
}

//____________________________________________________________________________ 
AliITSQASPDDataMakerRec::AliITSQASPDDataMakerRec(const AliITSQASPDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSPDhRawsTask(qadm.fSPDhRawsTask),
fSPDhDigitsTask(qadm.fSPDhDigitsTask),
fSPDhRecPointsTask(qadm.fSPDhRecPointsTask),
fGenRawsOffset(qadm.fGenRawsOffset),
fGenDigitsOffset(qadm.fGenDigitsOffset),
fGenRecPointsOffset(qadm.fGenRecPointsOffset),
fAdvLogger(qadm.fAdvLogger)
{
  //copy ctor 
  fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
  fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
  }

//__________________________________________________________________
AliITSQASPDDataMakerRec::~AliITSQASPDDataMakerRec(){
  // destructor
//  delete fAdvLogger;
}
//__________________________________________________________________

AliITSQASPDDataMakerRec& AliITSQASPDDataMakerRec::operator = (const AliITSQASPDDataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQASPDDataMakerRec();
  new(this) AliITSQASPDDataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SPD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list)
{
  // launch the QA checking
  //
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  if(!list){
    AliError(" Histogram list is NULL");
    return;
  }
  //
  Int_t shift = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  for (int trCl=-1;trCl<fAliITSQADataMakerRec->GetNTrigClasses();trCl++) { // RS Loop over all trigger classes (+ global counter, -1)
    if(task == AliQAv1::kRAWS) {
      TObjArray &harr = *fAliITSQADataMakerRec->GetRawsDataOfTrigClass(trCl);  // extract the clones for given trigger class
      //
      if (harr[5+shift]) ((TH1*)harr[5+shift])->Reset();  // clean up MEB histo (needed at the first cycle for small statistics)
      if (harr[7+shift]&&harr[1+shift]&&harr[0+shift]) ((TH1*)harr[7+shift])->Divide((TH1*)harr[1+shift],(TH1*)harr[0+shift]);
      if (harr[8+shift]&&harr[2+shift]&&harr[0+shift]) ((TH1*)harr[8+shift])->Divide((TH1*)harr[2+shift],(TH1*)harr[0+shift]); // missing FO ratio (per event) 
      if (harr[9+shift]&&harr[3+shift]&&harr[4+shift]) ((TH1*)harr[9+shift])->Divide((TH1*)harr[3+shift],(TH1*)harr[4+shift]); // noisy FO ratio   (per event)
      //
      TH1 *h0 = (TH1*)harr[0+shift];
      TH1 *h4 = (TH1*)harr[4+shift];
      TH2 *h5 = (TH2*)harr[5+shift];
      TH1 *h8 = (TH1*)harr[8+shift];
      TH1 *h9 = (TH1*)harr[9+shift];
      //
      if (h0 && h4 && h5 && h8 && h9) {
	for(Int_t i=0; i<1200; i++) {
	  // threshold for protection in case of technical runs (->few entries per chip)
	  Short_t thre=20; // 20 is ok in run 104792 (where the problem occured).
	  if (h0->GetBinContent(i+1)<thre) continue; // expected FO yield
	  if (h4->GetBinContent(i+1)<thre) continue; // the total FO yield.
	  //
	  if(h8->GetBinContent(i+1)>0.5 && h9->GetBinContent(i+1)>0.5) {
	    Int_t eq=i/60;
	    Int_t hs=(i%60)/10;
	    Int_t chip = (i%60)%10;
	    Int_t sect = 999; 	Int_t chipbin = 999; 
	    if(eq<10) {sect = eq; chipbin = 19-chip;} 
	    else {sect = eq-10; chipbin = chip;}
	    h5->Fill(sect*6+hs,chipbin);
	  }
	}
      }
    } // QA Task Index == RAWS
  } // RS: loop over trigger classes
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerRec::InitRaws()
{ 
 
  // Initialization for RAW data - SPD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  
  if(!fAdvLogger) fAdvLogger = new AliITSRawStreamSPDErrorLog();  
  AliDebug(AliQAv1::GetQADebugLevel(), "Book Offline Histograms for SPD\n ");

  Char_t name[50];
  Char_t title[50];
  // offset for online histogram numbering
  Int_t shift = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  
 
  Float_t range[2] = {-0.5,1199.5};

// **********  online histo booking (shift is added)   *********************

  //0
  TH1F *hFiredChips = new TH1F("SPDFiredChips_OnlineSPD","FiredChips - SPD",fgkSPDchips,range[0],range[1]);
  hFiredChips->GetXaxis()->SetTitle("chip index (eq*60 + hs*10 + chip)");
  hFiredChips->GetYaxis()->SetTitle("Fired Chip yield");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFiredChips, 0+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
  // 1
  TH1F *hFastOrFiredChips = new TH1F("SPDFastOrFiredChips_OnlineSPD","FastOr-Fired Chips (if pixel hit present) - SPD",fgkSPDchips,range[0],range[1]);
  hFastOrFiredChips->GetXaxis()->SetTitle("chip index (eq*60 + hs*10 + chip)");
  hFastOrFiredChips->GetYaxis()->SetTitle("FastOr-Fired Chip yield (per event)");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrFiredChips, 1+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 2
  TH1F *hFastOrMissing = new TH1F("SPDFastOrMissing_OnlineSPD","Missing FastOr signal - SPD",fgkSPDchips,range[0],range[1]);
  hFastOrMissing->GetXaxis()->SetTitle("chip index (eq*60 + hs*10 + chip)");
  hFastOrMissing->GetYaxis()->SetTitle("Missing Fast Or yield");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrMissing, 2+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 3
  TH1F *hFastOrNoisy = new TH1F("SPDFastOrNoisy_OnlineSPD","Noisy (no pixel hit present) FastOr signal - SPD",fgkSPDchips,range[0],range[1]);
  hFastOrNoisy->GetXaxis()->SetTitle("chipkey");
  hFastOrNoisy->GetYaxis()->SetTitle("Noisy Fast Or");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrNoisy, 3+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 4
  TH1F *hFastOrCumulative = new TH1F("SPDFastOrCumulative_OnlineSPD","Cumulative FastOr signal - SPD",fgkSPDchips,range[0],range[1]);
  hFastOrCumulative->GetXaxis()->SetTitle("chipkey");
  hFastOrCumulative->GetYaxis()->SetTitle("Cumulative Fast Or yield");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrCumulative, 4+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 5
  
   TH2I *hSPDChipsMEB = new TH2I("SPDChipsMEB_OnlineSPD","Chips with MEB problem - SPD",60,-0.5,59.5,20,-0.2,19.5);
  // hSPDChipsMEB->GetXaxis()->SetTitle("Stave");
   hSPDChipsMEB->GetXaxis()->SetNdivisions(60,kFALSE);
   hSPDChipsMEB->GetYaxis()->SetTitle("SIDE C   ->   SIDE A           Chip");
   hSPDChipsMEB->GetYaxis()->SetNdivisions(20,kFALSE);
   //hSPDChipsMEB->SetOption("COLZ");
   hSPDChipsMEB->SetOption("text");
   hSPDChipsMEB->UseCurrentStyle();
   for(Int_t ibinx =0; ibinx< hSPDChipsMEB->GetNbinsX(); ibinx++){
   if(ibinx%6==0) hSPDChipsMEB->GetXaxis()->SetBinLabel(ibinx+1,Form("Sector %i__%i",ibinx/6,ibinx%6));
   else hSPDChipsMEB->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",ibinx%6));
   }
   for(Int_t ibiny =0; ibiny< hSPDChipsMEB->GetNbinsY(); ibiny++){
   if(ibiny < 10) hSPDChipsMEB->GetYaxis()->SetBinLabel(ibiny+1,Form("%i",ibiny));
   else hSPDChipsMEB->GetYaxis()->SetBinLabel(ibiny+1,Form("%i",19-ibiny));
   }
   rv = fAliITSQADataMakerRec->Add2RawsList(hSPDChipsMEB, 5+shift, !expert, image, !saveCorr);
   fSPDhRawsTask++;  
// 6
  TH2F *hFastOrCorrelation = new TH2F("SPDFastOrCorrelation_OnlineSPD","Fast Or multiplicity correlation - SPD",100,0.,100.,100,0,100);
  hFastOrCorrelation->GetXaxis()->SetTitle("Layer 1");
  hFastOrCorrelation->GetYaxis()->SetTitle("Layer 2");
  hFastOrCorrelation->SetOption("COLZ");
  hFastOrCorrelation->UseCurrentStyle();
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrCorrelation, 6+shift, !expert, image, !saveCorr);
  fSPDhRawsTask++;
// 7
  TH1F *hFastOrEfficiency, *hFastOrMissingRatio, *hFastOrNoisyRatio;
  
  hFastOrEfficiency = new TH1F("SPDFastOrEfficiency_OnlineSPD","FastOr Efficiency : FastOr / fired chips (per event) - SPD",fgkSPDchips,range[0],range[1]);
  hFastOrEfficiency->SetFillColor(kBlue);
  hFastOrEfficiency->SetMaximum(1.05);
  hFastOrEfficiency->GetXaxis()->SetTitle("chip index [eq*60+hs*10+chip]");
  hFastOrEfficiency->GetYaxis()->SetTitle("FastOr Efficiency (per event)");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrEfficiency, 7+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 8 
  hFastOrMissingRatio = new TH1F("SPDFastOrMissingRatio_OnlineSPD"," Missing FastOr / fired chips (per event) - SPD)",fgkSPDchips,range[0],range[1]);
  hFastOrMissingRatio->SetFillColor(kBlue);
  hFastOrMissingRatio->SetMaximum(1.05);
  hFastOrMissingRatio->GetXaxis()->SetTitle("chip index [eq*60+hs*10+chip]");
  hFastOrMissingRatio->GetYaxis()->SetTitle("ratio of Missing FastOr (per event)");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrMissingRatio, 8+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 9
  hFastOrNoisyRatio = new TH1F("SPDFastOrNoisyRatio_OnlineSPD","Noisy Fast Or / fired chips (per event) - SPD",fgkSPDchips,range[0],range[1]);
  hFastOrNoisyRatio->SetFillColor(kBlue);
  hFastOrNoisyRatio->SetMaximum(1.05);
  hFastOrNoisyRatio->GetXaxis()->SetTitle("chip index [eq*60+hs*10+chip]");
  hFastOrNoisyRatio->GetYaxis()->SetTitle("ratio of Noisy FastOr (per event)");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrNoisyRatio, 9+shift, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 10
  TH2F *herrorsAll = new TH2F("SPDErrorsAll_OnlineSPD","Error codes - SPD",20,-0.5,19.5,22,-0.5,21.5);
  herrorsAll->GetXaxis()->SetTitle("DDL");
  herrorsAll->GetYaxis()->SetTitle("Error Type");
//  herrorsAll->SetOption("COLZ");
  herrorsAll->SetOption("text");
  herrorsAll->UseCurrentStyle();
  rv = fAliITSQADataMakerRec->Add2RawsList(herrorsAll, kAmoreFoOffset+shift, !expert, image, !saveCorr);
  fSPDhRawsTask++;
//11-30
  TH1F * herrors[20];
  for (Int_t iEq=0; iEq<20; iEq++) {
    snprintf(name,49,"SPDErrors_Eq%d_OnlineSPD",iEq+1);
    snprintf(title,49,"Error codes - SPD Eq %d",iEq+1);
    herrors[iEq] = new TH1F (name,title,fAdvLogger->GetNrErrorCodes(),0,fAdvLogger->GetNrErrorCodes());
    herrors[iEq]->SetXTitle("Error Code");
    herrors[iEq]->SetYTitle("Nr of errors");
    rv = fAliITSQADataMakerRec->Add2RawsList(herrors[iEq], 11+shift+iEq, expert, !image, !saveCorr);
    fSPDhRawsTask++;
  }
//   *********   offline histo booking  (offset is added) ****************************  

 // offset for offline histogram numbering
  Int_t offset = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + kAmoreFoOffset + kAmoreErrorsOffset;
  
 // printf("now booking offline raw data : genrawoffset %i, kAmoreOffset %i , kAmoreErrorsOffset %i -> total %i , list numbering %i\n",fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()],(Int_t)kAmoreFoOffset,(Int_t)kAmoreErrorsOffset, offset,fSPDhRawsTask);
// 0
  TH1F *hlayer = new TH1F("SPDLayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RawsList(hlayer, 0+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;

  TH1F * hmod[2];
  TH2F * hhitMap[20];

  
// 1-2
  for (Int_t iLay=0; iLay<2; iLay++) {
    snprintf(name,50,"SPDModPattern_SPD%d",iLay+1);
    snprintf(title,50,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,fgknSPDmodules,0,fgknSPDmodules);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RawsList(hmod[iLay], 1+iLay+offset, expert, !image, !saveCorr);
    fSPDhRawsTask++;
  }
// 3
  TH2F *hHitMapHalfStaveChipSideA 
     = new TH2F("SPDHitMapHalfStaveChipSideA_SPD","Hit map per HalfStave per Chip Side A - SPD",60,0.,60.,10,0.,10.);
  hHitMapHalfStaveChipSideA->GetXaxis()->SetTitle("HalfStave");
  hHitMapHalfStaveChipSideA->GetYaxis()->SetTitle("Chip");
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapHalfStaveChipSideA, 3+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 4
  TH2F *hHitMapHalfStaveChipSideC 
     = new TH2F("SPDHitMapHalfStaveChipSideC_SPD","Hit map per HalfStave per Chip Side C - SPD",60,0.,60.,10,0.,10.);
  hHitMapHalfStaveChipSideC->GetXaxis()->SetTitle("HalfStave");
  hHitMapHalfStaveChipSideC->GetYaxis()->SetTitle("Chip");
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapHalfStaveChipSideC, 4+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
 //5-24
  for (Int_t iDDL=0; iDDL<20; iDDL++) {
    snprintf(name,49,"SPDHitMap_SPD_DDL%d",iDDL+1);
    snprintf(title,49,"Hit map - SPD DDL %d",iDDL+1);
    hhitMap[iDDL]=new TH2F(name,title,320,0,10*32,1536,0,6*256);
    hhitMap[iDDL]->GetXaxis()->SetTitle("Column");
    hhitMap[iDDL]->GetYaxis()->SetTitle("Row");
    rv = fAliITSQADataMakerRec->Add2RawsList(hhitMap[iDDL], 5+iDDL+offset, expert, !image, !saveCorr);
    fSPDhRawsTask++;
    }
// 25-26
  TH1F * hMultSPDhits[2];
  for (Int_t iLay=0; iLay<2; iLay++) {
    snprintf(name,49,"SPDHitsMultiplicity_SPD%d",iLay+1);
    snprintf(title,49,"Hit multiplicity - SPD Layer %d",iLay+1);
    hMultSPDhits[iLay]=new TH1F(name,title,200,0.,200.);
    hMultSPDhits[iLay]->GetXaxis()->SetTitle("Hit multiplicity");
    hMultSPDhits[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RawsList(hMultSPDhits[iLay], 25+iLay+offset, expert, !image, !saveCorr);
    fSPDhRawsTask++;
  }
// 27
  TH2F *hMultSPDhits2MultSPDhits1 
         = new TH2F("SPDHitMultCorrelation_SPD","Hit multiplicity correlation - SPD",200,0.,200.,200,0.,200.);
  hMultSPDhits2MultSPDhits1->GetXaxis()->SetTitle("Hit multiplicity (Layer 1)");
  hMultSPDhits2MultSPDhits1->GetYaxis()->SetTitle("Hit multiplicity (Layer 2)");
  rv = fAliITSQADataMakerRec->Add2RawsList(hMultSPDhits2MultSPDhits1, 27+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 28
  TH2F *hFastOrMapStaveChip 
         = new TH2F("SPDFastOrMapStaveChip_SPD","FastOr map per Stave per Chip - SPD",60,-0.5,59.5,20,-0.5,19.5);
//  hFastOrMapStaveChip->GetXaxis()->SetTitle("Stave");
  hFastOrMapStaveChip->GetYaxis()->SetTitle("SIDE C   ->   SIDE A           Chip");
  for(Int_t ibinx =0; ibinx< hFastOrMapStaveChip->GetNbinsX(); ibinx++){
  if(ibinx%6==0) hFastOrMapStaveChip->GetXaxis()->SetBinLabel(ibinx+1,Form("Sector %i__%i",ibinx/6,ibinx%6));
  else hFastOrMapStaveChip->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",ibinx%6));
  }
  for(Int_t ibiny =0; ibiny< hFastOrMapStaveChip->GetNbinsY(); ibiny++){
  if(ibiny < 10) hFastOrMapStaveChip->GetYaxis()->SetBinLabel(ibiny+1,Form("%i",ibiny));
  else hFastOrMapStaveChip->GetYaxis()->SetBinLabel(ibiny+1,Form("%i",19-ibiny));
  }
  hFastOrMapStaveChip->SetOption("COLZ");
  hFastOrMapStaveChip->UseCurrentStyle();
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrMapStaveChip, 28+offset, !expert, image, !saveCorr);
  fSPDhRawsTask++;
// 29
  TH1F *hFastOrFiredMap = new TH1F("SPDFastOrPattern_SPD","FastOrFiredChip map - SPD",1200,0.,1200.);
  hFastOrFiredMap->GetXaxis()->SetTitle("Chip number");
  hFastOrFiredMap->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RawsList(hFastOrFiredMap, 29+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 30
  TH2F *hHitMapHalfStaveChipInner 
     = new TH2F("SPDHitMapStaveChipInner_SPD","Hit map per Stave per Chip Inner Layer- SPD",20,0.,20.,20,0.,20.);
  hHitMapHalfStaveChipInner->GetXaxis()->SetTitle("SIDE C                                 SIDE A                   Chip");
  //hHitMapHalfStaveChipInner->GetYaxis()->SetTitle("Stave in Sector S");
  for(Int_t ibinx =0; ibinx< hHitMapHalfStaveChipInner->GetNbinsX(); ibinx++){
  if(ibinx < 10) hHitMapHalfStaveChipInner->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",ibinx));
  else hHitMapHalfStaveChipInner->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",20-(ibinx+1)));
  }
  for(Int_t ibiny =0; ibiny< hHitMapHalfStaveChipInner->GetNbinsY(); ibiny++){
  if(ibiny%2==1) hHitMapHalfStaveChipInner->GetYaxis()->SetBinLabel(ibiny+1,Form(" S %i - %i",ibiny/2,ibiny%2));
  else hHitMapHalfStaveChipInner->GetYaxis()->SetBinLabel(ibiny+1,Form("%i",ibiny%2));
  hHitMapHalfStaveChipInner->GetYaxis()->SetTitleOffset(1.4);
  }
  hHitMapHalfStaveChipInner->SetOption("COLZ");
  hHitMapHalfStaveChipInner->UseCurrentStyle();
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapHalfStaveChipInner, 30+offset, !expert, image, !saveCorr);
  fSPDhRawsTask++;
// 31
  TH2F *hHitMapHalfStaveChipOuter 
     = new TH2F("SPDHitMapStaveChipOuter_SPD","Hit map per Stave per Chip Outer Layer - SPD",20,0.,20.,40,0.,40.);
  hHitMapHalfStaveChipOuter->GetXaxis()->SetTitle("SIDE C                                 SIDE A                   Chip");
  //hHitMapHalfStaveChipOuter->GetYaxis()->SetTitle("Stave in Sector S");
  for(Int_t ibinx =0; ibinx< hHitMapHalfStaveChipOuter->GetNbinsX(); ibinx++){
  if(ibinx < 10) hHitMapHalfStaveChipOuter->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",ibinx));
  else hHitMapHalfStaveChipOuter->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",20-(ibinx+1)));
  }
  for(Int_t ibiny =0; ibiny< hHitMapHalfStaveChipOuter->GetNbinsY(); ibiny++){
  if(ibiny%4==3) hHitMapHalfStaveChipOuter->GetYaxis()->SetBinLabel(ibiny+1,Form(" S %i - %i",ibiny/4,ibiny%4+2));
  else hHitMapHalfStaveChipOuter->GetYaxis()->SetBinLabel(ibiny+1,Form("%i",ibiny%4+2));
  hHitMapHalfStaveChipOuter->GetYaxis()->SetTitleOffset(1.4);
  }
  hHitMapHalfStaveChipOuter->SetOption("COLZ");
  hHitMapHalfStaveChipOuter->UseCurrentStyle();
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapHalfStaveChipOuter, 31+offset, !expert, image, !saveCorr);
  fSPDhRawsTask++;
// 32
  TH1F *hHitMapChipInnerZ = new TH1F("SPDHitMapChipInnerZ_SPD","Hit map per ChipZ Inner - SPD",20,0.,20.);
  hHitMapChipInnerZ->GetXaxis()->SetTitle("Chip");
  hHitMapChipInnerZ->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapChipInnerZ, 32+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 33
  TH1F *hHitMapChipOuterZ = new TH1F("SPDHitMapChipOuterZ_SPD","Hit map per ChipZ Outer - SPD",20,0.,20.);
  hHitMapChipOuterZ->GetXaxis()->SetTitle("Chip");
  hHitMapChipOuterZ->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapChipOuterZ, 33+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 34
  TH1F *hHitMapStaveInnerPhi = new TH1F("SPDHitMapChipInnerPhi_SPD","Hit map per Stave in Phi Inner - SPD",20,0.,20.);
  //hHitMapStaveInnerPhi->GetXaxis()->SetTitle("Stave");
  hHitMapStaveInnerPhi->GetYaxis()->SetTitle("Entries");
  for(Int_t ibinx =0; ibinx< hHitMapStaveInnerPhi->GetNbinsX(); ibinx++){
  if(ibinx%2==0) hHitMapStaveInnerPhi->GetXaxis()->SetBinLabel(ibinx+1,Form("%i___Sector %i",ibinx%2,ibinx/2));
  else hHitMapStaveInnerPhi->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",ibinx%2));
  }
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapStaveInnerPhi, 34+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
// 35
  TH1F *hHitMapStaveOuterPhi = new TH1F("SPDHitMapChipOuterPhi_SPD","Hit map per Stave in Phi Outer - SPD",40,0.,40.);
  //hHitMapStaveOuterPhi->GetXaxis()->SetTitle("Stave");
  hHitMapStaveOuterPhi->GetYaxis()->SetTitle("Entries");
  for(Int_t ibinx =0; ibinx< hHitMapStaveOuterPhi->GetNbinsX(); ibinx++){
  if(ibinx%4==0) hHitMapStaveOuterPhi->GetXaxis()->SetBinLabel(ibinx+1,Form("%i___Sector %i ",ibinx%4+2,ibinx/4));
  else hHitMapStaveOuterPhi->GetXaxis()->SetBinLabel(ibinx+1,Form("%i",ibinx%4+2));
  }
  rv = fAliITSQADataMakerRec->Add2RawsList(hHitMapStaveOuterPhi, 35+offset, expert, !image, !saveCorr);
  fSPDhRawsTask++;
   
  //AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SPD Raws histograms booked\n",fSPDhRawsTask));
  //printf("------------------         %d SPD Raws histograms booked                \n",fSPDhRawsTask);
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASPDDataMakerRec::MakeRaws(AliRawReader* rawReader)
{ 
  // Fill QA for RAW - SPD -
  Int_t rv = 0 ; 
  if(!rawReader) {
   AliError("rawReader is NULL"); 
   return -1;
  } 
  rawReader->Reset();
  AliITSRawStreamSPD rawStreamSPD(rawReader);
  rawStreamSPD.ActivateAdvancedErrorLog(kTRUE,fAdvLogger);
  // shift for online histos
  Int_t shift = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];
  
   // shift for offline histos
  Int_t offset = fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + kAmoreFoOffset+kAmoreErrorsOffset;

  Int_t nDigitsL1 = 0;
  Int_t nDigitsL2 = 0;
  Int_t iEq;
  Int_t iLayer;
  Int_t iHalfStave, iChip;
  Int_t chipKey;
  Int_t col, row; 
  UInt_t module, colM, rowM;
  Bool_t isFiredChip[1200];
  for(Int_t ichK=0; ichK<1200; ichK++) isFiredChip[ichK] = kFALSE;
  Bool_t isOnlineFiredChip[1200];
  for(Int_t iOnlineKey=0; iOnlineKey<1200; iOnlineKey++) isOnlineFiredChip[iOnlineKey] = kFALSE;
  UInt_t nFastOr[2]={0,0};

  while(rawStreamSPD.Next()) {

    iEq = rawReader->GetDDLID();
    if (iEq>=0 && iEq<20) {
      iHalfStave = rawStreamSPD.GetHalfStaveNr();
      iChip = rawStreamSPD.GetChipAddr();
      col   = rawStreamSPD.GetChipCol();
      row   = rawStreamSPD.GetChipRow();
      isOnlineFiredChip[iEq*60+iHalfStave*10+iChip] = kTRUE;
      chipKey = rawStreamSPD.GetOfflineChipKeyFromOnline(iEq,iHalfStave,iChip);
      isFiredChip[chipKey] = kTRUE;

      rawStreamSPD.OnlineToOffline(iEq, iHalfStave, iChip, col, row, module, colM, rowM);

      if (iHalfStave>=0 && iHalfStave<2) iLayer=0;
      else iLayer=1;
      
      fAliITSQADataMakerRec->FillRawsData(0+offset,iLayer);
      if (iLayer==0) {
        fAliITSQADataMakerRec->FillRawsData(1+offset,module);
        nDigitsL1++;
      } else {
        fAliITSQADataMakerRec->FillRawsData(2+offset,module);
        nDigitsL2++;
      }
      
      if(iEq<10) {
         fAliITSQADataMakerRec->FillRawsData(3+offset,iHalfStave+iEq*6,iChip);
      } 
      else       {
         fAliITSQADataMakerRec->FillRawsData(4+offset,iHalfStave+(iEq-10)*6,iChip);
      }

      if(iLayer==0) {
         if(iEq<10)  { 
            fAliITSQADataMakerRec->FillRawsData(30+offset,19-iChip,iHalfStave+iEq*2);
            fAliITSQADataMakerRec->FillRawsData(32+offset,19-iChip);
            fAliITSQADataMakerRec->FillRawsData(34+offset,iHalfStave+iEq*2);
         }
         else {
            fAliITSQADataMakerRec->FillRawsData(30+offset,iChip,iHalfStave+(iEq-10)*2);
            fAliITSQADataMakerRec->FillRawsData(32+offset,iChip);
            fAliITSQADataMakerRec->FillRawsData(34+offset,iHalfStave+(iEq-10)*2);
         }
      }
      else         {   
         if(iEq<10)  { 
            fAliITSQADataMakerRec->FillRawsData(31+offset,19-iChip,iHalfStave-2+iEq*4);
            fAliITSQADataMakerRec->FillRawsData(33+offset,19-iChip);
            fAliITSQADataMakerRec->FillRawsData(35+offset,iHalfStave-2+iEq*4);
         }
         else {
            fAliITSQADataMakerRec->FillRawsData(31+offset,iChip,iHalfStave-2+(iEq-10)*4);
            fAliITSQADataMakerRec->FillRawsData(33+offset,iChip);
            fAliITSQADataMakerRec->FillRawsData(35+offset,iHalfStave-2+(iEq-10)*4);
         }
      }
      fAliITSQADataMakerRec->FillRawsData(5+iEq+offset,colM+(module%2)*160,rowM+iHalfStave*256); 
    }
  }

  UInt_t nErrorsDDL[20];
  for (Int_t ieq=0; ieq<20; ieq++) {
    nErrorsDDL[ieq] = 0;
    if (!fAdvLogger) continue;
    for (UInt_t ierr=0; ierr<fAdvLogger->GetNrErrorCodes(); ierr++) {
      fAliITSQADataMakerRec->FillRawsData(ieq+(kAmoreFoOffset+1)+shift,ierr,fAdvLogger->GetNrErrors(ierr,ieq));
      if(ierr>0) {
       nErrorsDDL[ieq] = nErrorsDDL[ieq] + fAdvLogger->GetNrErrors(ierr,ieq); 
     }
      fAliITSQADataMakerRec->FillRawsData(kAmoreFoOffset+shift,ieq,ierr,fAdvLogger->GetNrErrors(ierr,ieq));
    } 
    for (Int_t ihs=0; ihs<6; ihs++) {
      for (Int_t ichip=0; ichip<10; ichip++) {
      if(isOnlineFiredChip[ieq*60+ihs*10+ichip]) fAliITSQADataMakerRec->FillRawsData(0+shift,ieq*60+ihs*10+ichip); // online
       if(rawStreamSPD.GetFastOrSignal(ieq,ihs,ichip)) fAliITSQADataMakerRec->FillRawsData(4+shift,ieq*60+ihs*10+ichip); // online
       // now filling the 3 possibile combinations
      if(rawStreamSPD.GetFastOrSignal(ieq,ihs,ichip) && isOnlineFiredChip[ieq*60+ihs*10+ichip]) fAliITSQADataMakerRec->FillRawsData(1+shift,ieq*60+ihs*10+ichip); // online
      if(!rawStreamSPD.GetFastOrSignal(ieq,ihs,ichip) && isOnlineFiredChip[ieq*60+ihs*10+ichip]) fAliITSQADataMakerRec->FillRawsData(2+shift,ieq*60+ihs*10+ichip); // online
      if(rawStreamSPD.GetFastOrSignal(ieq,ihs,ichip) && !isOnlineFiredChip[ieq*60+ihs*10+ichip]) fAliITSQADataMakerRec->FillRawsData(3+shift,ieq*60+ihs*10+ichip); // online       
      
        chipKey = rawStreamSPD.GetOfflineChipKeyFromOnline(ieq,ihs,ichip);
        Int_t sect = 999;
	Int_t chipbin = 999;
	 
        if(rawStreamSPD.GetFastOrSignal(ieq,ihs,ichip)) {
          if(ihs <2) nFastOr[0]++; // online
	  else nFastOr[1]++;       // online
	  if(ieq<10) {sect = ieq; chipbin = 19-ichip;} 
	  else {sect = ieq-10; chipbin = ichip;}
	  
	  fAliITSQADataMakerRec->FillRawsData(28+offset,sect*6+ihs,chipbin);
          fAliITSQADataMakerRec->FillRawsData(29+offset,chipKey);
        }
      }
    } 


  }
  fAliITSQADataMakerRec->FillRawsData(6+shift,nFastOr[0],nFastOr[1]); // online
 
  fAdvLogger->Reset();
  fAliITSQADataMakerRec->FillRawsData(25+offset,nDigitsL1);
  fAliITSQADataMakerRec->FillRawsData(26+offset,nDigitsL2);
  fAliITSQADataMakerRec->FillRawsData(27+offset,nDigitsL1,nDigitsL2);
  //
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Event completed, %d raw digits read",nDigitsL1+nDigitsL2));
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerRec::InitDigits()
{ 
  // Initialization for DIGIT data - SPD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ;
  Int_t rv = 0 ; 
//  fGenDigitsOffset = (fAliITSQADataMakerRec->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  //fSPDhDigitsTask must be incremented by one unit every time a histogram is ADDED to the QA List
  
  Char_t name[50];
  Char_t title[50];
  
  TH1F *hlayer = new TH1F("SPDLayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2DigitsList(hlayer,fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  fSPDhDigitsTask++;
  
  TH1F * hmod[2];
  for (Int_t iLay=0; iLay<2; iLay++) {
    snprintf(name,49,"SPDModPattern_SPD%d",iLay+1);
    snprintf(title,49,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,240,0,240);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2DigitsList(hmod[iLay],1+iLay+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
    fSPDhDigitsTask++;
  }
  
  TH1F *hcolumns = new TH1F("SPDColumns_SPD","Columns - SPD",160,0.,160.);
  hcolumns->GetXaxis()->SetTitle("Column number");
  hcolumns->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2DigitsList(hcolumns,3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  fSPDhDigitsTask++;
  
  TH1F *hrows = new TH1F("SPDRows_SPD","Rows - SPD",256,0.,256.);
  hrows->GetXaxis()->SetTitle("Row number");
  hrows->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2DigitsList(hrows,4+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  fSPDhDigitsTask++;
  
  TH1F * hMultSPDdigits[2];
  for (Int_t iLay=0; iLay<2; ++iLay) {
    snprintf(name,49,"SPDDigitMultiplicity_SPD%d",iLay+1);
    snprintf(title,49,"Digit multiplicity - SPD Layer %d",iLay+1);
    hMultSPDdigits[iLay]=new TH1F(name,title,200,0.,200.);
    hMultSPDdigits[iLay]->GetXaxis()->SetTitle("Digit multiplicity");
    hMultSPDdigits[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2DigitsList(hMultSPDdigits[iLay], 5+iLay+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
    fSPDhDigitsTask++;
  }
  
  TH2F *hMultSPDdig2MultSPDdig1 
    = new TH2F("SPDDigitMultCorrelation_SPD","Digit multiplicity correlation - SPD",200,0.,200.,200,0.,200.);
  hMultSPDdig2MultSPDdig1->GetXaxis()->SetTitle("Digit multiplicity (Layer 1)");
  hMultSPDdig2MultSPDdig1->GetYaxis()->SetTitle("Digit multiplicity (Layer 2)");
  rv = fAliITSQADataMakerRec->Add2DigitsList(hMultSPDdig2MultSPDdig1,7+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSPDhDigitsTask++;
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SPD Digits histograms booked\n",fSPDhDigitsTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASPDDataMakerRec::MakeDigits(TTree *digits)
{ 
  // Fill QA for DIGIT - SPD -
  
  Int_t rv = 0 ; 
  //  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  //  fITS->SetTreeAddress();
  //  TClonesArray *iITSdigits  = fITS->DigitsAddress(0);  // 0->SPD
  TBranch *branchD = digits->GetBranch("ITSDigitsSPD");
  if (!branchD) { 
    AliError("can't get the branch with the SPD ITS digits !");
    return rv;
  }
  static TClonesArray statDigits("AliITSdigitSPD");
  TClonesArray *iITSdigits = &statDigits;
  branchD->SetAddress(&iITSdigits);  
  Int_t nDigitsL1=0;
  Int_t nDigitsL2=0;
  
  for (Int_t imod=0; imod<240; ++imod){
    digits->GetEvent(imod);
    Int_t ndigits = iITSdigits->GetEntries();
    if (imod<80) {
      fAliITSQADataMakerRec->FillDigitsData(0+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],0.5,ndigits);
      fAliITSQADataMakerRec->FillDigitsData(1+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],imod,ndigits);
      nDigitsL1+=ndigits;
    }
    else {
      fAliITSQADataMakerRec->FillDigitsData(0+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],1,ndigits);
      fAliITSQADataMakerRec->FillDigitsData(2+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],imod,ndigits);
      nDigitsL2+=ndigits;
    }
    for (Int_t idig=0; idig<ndigits; ++idig) {
      AliITSdigit *dig=(AliITSdigit*)iITSdigits->UncheckedAt(idig);
      Int_t col=dig->GetCoord1();  // cell number z
      Int_t row=dig->GetCoord2();  // cell number x
      fAliITSQADataMakerRec->FillDigitsData(3+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],col);
      fAliITSQADataMakerRec->FillDigitsData(4+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],row);
    }
  }
  fAliITSQADataMakerRec->FillDigitsData(5+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],nDigitsL1);
  fAliITSQADataMakerRec->FillDigitsData(6+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],nDigitsL2);
  fAliITSQADataMakerRec->FillDigitsData(7+fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()],nDigitsL1,nDigitsL2);
  //
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SPD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
//  fGenRecPointsOffset = (fAliITSQADataMakerRec->fRecPointsQAList[AliRecoParam::kDefault])->GetEntries();
  TH1F* hlayer= new TH1F("SPDLayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(hlayer, 0+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); 
  fSPDhRecPointsTask++;

  TH1F * hmod[2];
  TH1F * hxl[2];
  TH1F * hzl[2];
  TH1F * hxg[2];
  TH1F * hyg[2];
  TH1F * hzg[2];
  TH1F * hr[2];
  TH1F * hphi[2];
  TH1F * hMultSPDcl[2];
  TH2F * hNyNz[2];  // y and z cluster length
  TH1F * hNpixels[2];  // cluster size in number of pixels
  TH1F * hType[2];  // cluster type according to conventional table
  TH2F * hPhiZ[2];

  Float_t xlim[2]={4.5,8.};
  Float_t zlim[2]={15.,15.};

  Char_t name[50];
  Char_t title[50];
  for (Int_t iLay=0;iLay<2;iLay++) {
    snprintf(name,49,"SPDModPattern_SPD%d",iLay+1);
    snprintf(title,49,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,fgknSPDmodules,0,fgknSPDmodules);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hmod[iLay], 1+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); 
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDxLoc_SPD%d",iLay+1);
    snprintf(title,49,"Local x coordinate - SPD Layer %d",iLay+1);
    hxl[iLay]=new TH1F(name,title,100,-4.,4.);
    hxl[iLay]->GetXaxis()->SetTitle("Local x [cm]");
    hxl[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hxl[iLay], 2+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDzLoc_SPD%d",iLay+1);
    snprintf(title,49,"Local z coordinate - SPD Layer %d",iLay+1);
    hzl[iLay]=new TH1F(name,title,100,-4.,4.);
    hzl[iLay]->GetXaxis()->SetTitle("Local z [cm]");
    hzl[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hzl[iLay], 3+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); 
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDxGlob_SPD%d",iLay+1);
    snprintf(title,49,"Global x coordinate - SPD Layer %d",iLay+1);
    hxg[iLay]=new TH1F(name,title,100,-xlim[iLay],xlim[iLay]);
    hxg[iLay]->GetXaxis()->SetTitle("Global x [cm]");
    hxg[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hxg[iLay],4+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);  
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDyGlob_SPD%d",iLay+1);
    snprintf(title,49,"Global y coordinate - SPD Layer %d",iLay+1);
    hyg[iLay]=new TH1F(name,title,100,-xlim[iLay],xlim[iLay]);
    hyg[iLay]->GetXaxis()->SetTitle("Global y [cm]");
    hyg[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hyg[iLay], 5+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); 
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDzGlob_SPD%d",iLay+1);
    snprintf(title,49,"Global z coordinate - SPD Layer %d",iLay+1);
    hzg[iLay]=new TH1F(name,title,150,-zlim[iLay],zlim[iLay]);
    hzg[iLay]->GetXaxis()->SetTitle("Global z [cm]");
    hzg[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hzg[iLay], 6+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); 
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDr_SPD%d",iLay+1);
    snprintf(title,49,"Radius - SPD Layer %d",iLay+1);
    hr[iLay]=new TH1F(name,title,100,0.,10.);
    hr[iLay]->GetXaxis()->SetTitle("r [cm]");
    hr[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hr[iLay], 7+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);  
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDphi_SPD%d",iLay+1);
    snprintf(title,49,"#varphi - SPD Layer %d",iLay+1);
    hphi[iLay]=new TH1F(name,title,1000,0.,2*TMath::Pi());
    hphi[iLay]->GetXaxis()->SetTitle("#varphi [rad]");
    hphi[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hphi[iLay], 8+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
    fSPDhRecPointsTask++;
    
    snprintf(name,49,"SPDSizeYvsZ_SPD%d",iLay+1);
    snprintf(title,49,"Cluster dimension - SPD Layer %d",iLay+1);
    hNyNz[iLay]=new TH2F(name,title,100,0.,100.,100,0.,100.);
    hNyNz[iLay]->GetXaxis()->SetTitle("z length");
    hNyNz[iLay]->GetYaxis()->SetTitle("y length");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hNyNz[iLay], 9+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image); 
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDSizeTot_SPD%d",iLay+1);
    snprintf(title,49,"Cluster size - SPD Layer %d",iLay+1);
    hNpixels[iLay]=new TH1F(name,title,100,0.,100.);
    hNpixels[iLay]->GetXaxis()->SetTitle("Cluster size");
    hNpixels[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hNpixels[iLay], 10+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);  
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDType_SPD%d",iLay+1);
    snprintf(title,49,"Cluster type - SPD Layer %d",iLay+1);
    hType[iLay]=new TH1F(name,title,20,0.,20.);
    hType[iLay]->GetXaxis()->SetTitle("Cluster type");
    hType[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hType[iLay], 11+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);  
    fSPDhRecPointsTask++;

    snprintf(name,49,"SPDphi_z_SPD%d",iLay+1);
    snprintf(title,49,"#varphi vs z - SPD Layer %d",iLay+1);
    hPhiZ[iLay]=new TH2F(name,title,150,-zlim[iLay],zlim[iLay],200,0.,2*TMath::Pi());
    hPhiZ[iLay]->GetXaxis()->SetTitle("Global z [cm]");
    hPhiZ[iLay]->GetYaxis()->SetTitle("#varphi [rad]");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hPhiZ[iLay], 12+(12*iLay)+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
    fSPDhRecPointsTask++;

  }

  TH2F *hrPhi=new TH2F("SPDr_phi_SPD","#varphi vs r - SPD",100,0.,10.,100,0.,2*TMath::Pi());
  hrPhi->GetXaxis()->SetTitle("r [cm]");
  hrPhi->GetYaxis()->SetTitle("#varphi [rad]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(hrPhi, 25+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], expert, !image);
  fSPDhRecPointsTask++;

  TH2F *hxy=new TH2F("SPDx_y_SPD","Global y vs x - SPD",200,-10.,10.,200,-10.,10.);
  hxy->GetXaxis()->SetTitle("Global x [cm]");
  hxy->GetYaxis()->SetTitle("Global y [cm]");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(hxy, 26+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSPDhRecPointsTask++;

  for (Int_t iLay=0;iLay<2;iLay++) {
    snprintf(name,49,"SPDMultiplicity_SPD%d",iLay+1);
    snprintf(title,49,"Cluster multiplicity - SPD Layer %d",iLay+1);
    hMultSPDcl[iLay]=new TH1F(name,title,200,0.,200.);
    hMultSPDcl[iLay]->GetXaxis()->SetTitle("Cluster multiplicity");
    hMultSPDcl[iLay]->GetYaxis()->SetTitle("Entries");
    rv = fAliITSQADataMakerRec->Add2RecPointsList(hMultSPDcl[iLay], 27+iLay+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
    fSPDhRecPointsTask++;
  } 

  TH2F *hMultSPDcl2MultSPDcl1 =
        new TH2F("SPDMultCorrelation_SPD","Cluster multiplicity correlation - SPD",200,0.,200.,200,0.,200.);
  hMultSPDcl2MultSPDcl1->GetXaxis()->SetTitle("Clusters multiplicity (Layer 1)");
  hMultSPDcl2MultSPDcl1->GetYaxis()->SetTitle("Clusters multiplicity (Layer 2)"); 
  rv = fAliITSQADataMakerRec->Add2RecPointsList(hMultSPDcl2MultSPDcl1, 29+fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()], !expert, image);
  fSPDhRecPointsTask++;

  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SPD Recs histograms booked\n",fSPDhRecPointsTask));

  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASPDDataMakerRec::MakeRecPoints(TTree * clusterTree)
{
  // Fill QA for RecPoints - SPD -

  Int_t rv = 0 ;
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  TClonesArray *recpoints = NULL;
  rpcont->FetchClusters(0,clusterTree);
  if(!rpcont->GetStatusOK() || !rpcont->IsSPDActive()){
    AliError("can't get SPD clusters !");
    return rv;
  }

  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  Int_t nSPDmod = AliITSgeomTGeo::GetModuleIndex(3,1,1);

  Float_t cluGlo[3] = {0.,0.,0.};
  Int_t nClusters[2] = {0,0};

  for (Int_t iIts=0; iIts < nSPDmod; iIts++) {
    recpoints = rpcont->UncheckedGetClusters(iIts);
    Int_t nCluster = recpoints->GetEntriesFast();
    if(nCluster == 0)continue;
    // loop over clusters
    while(nCluster--) {
      AliITSRecPoint* cluster =
                      (AliITSRecPoint*)recpoints->UncheckedAt(nCluster);
      if (cluster->GetLayer()>1)continue;
      Int_t lay=cluster->GetLayer();
      fAliITSQADataMakerRec->FillRecPointsData(0 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],lay);
      cluster->GetGlobalXYZ(cluGlo);
      Float_t rad=TMath::Sqrt(cluGlo[0]*cluGlo[0]+cluGlo[1]*cluGlo[1]);
        Float_t phi= TMath::Pi() + TMath::ATan2(-cluGlo[1],-cluGlo[0]);
        if (lay==0) {
          fAliITSQADataMakerRec->FillRecPointsData(1 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],iIts);
          fAliITSQADataMakerRec->FillRecPointsData(2 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetDetLocalX());
          fAliITSQADataMakerRec->FillRecPointsData(3 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetDetLocalZ());
          fAliITSQADataMakerRec->FillRecPointsData(4 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[0]);
          fAliITSQADataMakerRec->FillRecPointsData(5 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[1]);
          fAliITSQADataMakerRec->FillRecPointsData(6 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[2]);
          fAliITSQADataMakerRec->FillRecPointsData(7 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],rad);
          fAliITSQADataMakerRec->FillRecPointsData(8 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],phi);
          fAliITSQADataMakerRec->FillRecPointsData(9 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetNz(),cluster->GetNy());
          fAliITSQADataMakerRec->FillRecPointsData(10 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetNpixels());
          fAliITSQADataMakerRec->FillRecPointsData(11 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetSPDclusterType());
          fAliITSQADataMakerRec->FillRecPointsData(12 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[2],phi);
        } else  {
          fAliITSQADataMakerRec->FillRecPointsData(13 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],iIts);
          fAliITSQADataMakerRec->FillRecPointsData(14 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetDetLocalX());
          fAliITSQADataMakerRec->FillRecPointsData(15 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetDetLocalZ());
          fAliITSQADataMakerRec->FillRecPointsData(16 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[0]);
          fAliITSQADataMakerRec->FillRecPointsData(17 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[1]);
          fAliITSQADataMakerRec->FillRecPointsData(18 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[2]);
          fAliITSQADataMakerRec->FillRecPointsData(19 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],rad);
          fAliITSQADataMakerRec->FillRecPointsData(20 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],phi);
          fAliITSQADataMakerRec->FillRecPointsData(21 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetNz(),cluster->GetNy());
          fAliITSQADataMakerRec->FillRecPointsData(22 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetNpixels());
          fAliITSQADataMakerRec->FillRecPointsData(23 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluster->GetSPDclusterType());
          fAliITSQADataMakerRec->FillRecPointsData(24 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[2],phi);
        }
        fAliITSQADataMakerRec->FillRecPointsData(25 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],rad,phi);
        fAliITSQADataMakerRec->FillRecPointsData(26 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],cluGlo[0],cluGlo[1]);

        nClusters[lay]++;
    } // end of cluster loop
  } // end of its "subdetector" loop

  for (Int_t iLay=0; iLay<2; iLay++)
    fAliITSQADataMakerRec->FillRecPointsData(27 +iLay +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],nClusters[iLay]);

  fAliITSQADataMakerRec->FillRecPointsData(29 +fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()],nClusters[0],nClusters[1]);
  //
  return rv ;
}



//_______________________________________________________________

Int_t AliITSQASPDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task,Int_t specie) const {
  // Returns offset number according to the specified task
  Int_t offset=0;
  if( task == AliQAv1::kRAWS ) {
    offset=fGenRawsOffset[specie];
  }
  else if( task == AliQAv1::kDIGITSR ) {
    offset=fGenDigitsOffset[specie];
  }
  else if( task == AliQAv1::kRECPOINTS ) {
    offset=fGenRecPointsOffset[specie];
  }

  return offset;
}

//_______________________________________________________________

void AliITSQASPDDataMakerRec::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset, Int_t specie) {
  // Returns offset number according to the specified task
  if( task == AliQAv1::kRAWS ) {
    fGenRawsOffset[specie]=offset;
  }
  else if( task == AliQAv1::kDIGITSR ) {
    fGenDigitsOffset[specie]=offset;
  }
  else if( task == AliQAv1::kRECPOINTS ) {
    fGenRecPointsOffset[specie]=offset;
  }
}

//_______________________________________________________________

Int_t AliITSQASPDDataMakerRec::GetTaskHisto(AliQAv1::TASKINDEX_t task) const {
  // Returns the number of histograms associated to the specified task

  Int_t histotot=0;

  if( task == AliQAv1::kRAWS ) {
    histotot=fSPDhRawsTask;
  }
  else if( task == AliQAv1::kDIGITSR ) {
    histotot=fSPDhDigitsTask;
  }
  else if( task == AliQAv1::kRECPOINTS ){
    histotot=fSPDhRecPointsTask;
  }

  return histotot;
}

