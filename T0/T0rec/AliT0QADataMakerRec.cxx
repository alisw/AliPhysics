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


/* $Id$ */

//---
//  Produces the data needed to calculate the quality assurance. 
//  Alla.Maevskaya@cern.ch
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h> 
#include <TDirectory.h>
#include <TMath.h>
// --- Standard library ---

// --- AliRoot header files ---

#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliT0digit.h" 
#include "AliT0hit.h"
#include "AliT0RecPoint.h"
#include "AliT0QADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliT0RawReader.h"
#include "AliT0RecoParam.h"
#include "AliQAThresholds.h"
#include "AliDAQ.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "THnSparse.h"

#include "TFitResultPtr.h"

#include "Riostream.h"
ClassImp(AliT0QADataMakerRec)
           
//____________________________________________________________________________ 
  AliT0QADataMakerRec::AliT0QADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kT0), 
		  "T0 Quality Assurance Data Maker"),
    fMeanRawVertexParam(0),
    fMeanORAParam(0),
    fMeanORCParam(0),
    fCFDEffSubRangeLowParam(0),
    fCFDEffSubRangeHighParam(3000),
    fLEDEffSubRangeLowParam(0),
    fLEDEffSubRangeHighParam(3000)
//  fnEventCal(0),
//  fnEventPhys(0)
{
  // ctor
  // RS: There is some inconsistency here: the separation of physics and calib. events/histos is done by
  // fEventSpecie. Why do we book separate histos on different slots for calib and physics ? 
  // I am changing this in such way that we don't need local counters like fNumTriggers (the corresponding
  // histos now incremented in the MakeRaws, and for the normalization I will use the framework's counters
  // AliQADataMaker::GetEvCountCycle(...), AliQADataMaker::GetEvCountTotal(...)
  // All these fTrEff.. feff.. will by directly filled in corresponding histos

  for(Int_t i=0; i<24; i++){
    fMeanCFDFromGoodRunParam[i]=0; 
  }
}


//____________________________________________________________________________ 
AliT0QADataMakerRec::AliT0QADataMakerRec(const AliT0QADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fMeanRawVertexParam(qadm.fMeanRawVertexParam),
  fMeanORAParam(qadm.fMeanORAParam),
  fMeanORCParam(qadm.fMeanORCParam),
  fCFDEffSubRangeLowParam(qadm.fCFDEffSubRangeLowParam),
  fCFDEffSubRangeHighParam(qadm.fCFDEffSubRangeHighParam),
  fLEDEffSubRangeLowParam(qadm.fLEDEffSubRangeLowParam),
  fLEDEffSubRangeHighParam(qadm.fLEDEffSubRangeHighParam)
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle());
  for(Int_t i=0; i<24; i++){
    fMeanCFDFromGoodRunParam[i]=qadm.fMeanCFDFromGoodRunParam[i]; 
  }
}

//__________________________________________________________________
AliT0QADataMakerRec& AliT0QADataMakerRec::operator = (const AliT0QADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliT0QADataMakerRec();
  new(this) AliT0QADataMakerRec(qadm);
  return *this;
}
//__________________________________________________________________
AliT0QADataMakerRec::~AliT0QADataMakerRec()
{
  //destructor
}
//____________________________________________________________________________
void AliT0QADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliInfo(Form("Task: %d",task));
  ResetEventTrigClasses();
  
  TH1* hcounter = 0;
  TH1* heff = 0;
  TH1* htmp = 0;
  TH1F* hEventCounter=NULL;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    //
    // RS: There is some inconsistency here: the separation of physics and calib. events/histos is done by
    // fEventSpecie. Why do we book separate histos on different slots for calib and physics ? 
    // I am changing this in such way that we don't need local counters like fNumTriggers (the corresponding
    // histos now incremented in the MakeRaws, and for the normalization I will use the framework's counters
    // AliQADataMaker::GetEvCountCycle(...), AliQADataMaker::GetEvCountTotal(...)
    //
    // I think the histos xx+250 should be suppressed (the xx calib histos of specie==calibration will be 
    // used automatically)
    //
    if (! IsValidEventSpecie(specie, list)) continue;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    //
    for (int itc=-1;itc<GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class
      // 
      if ( task == AliQAv1::kRAWS ) {
	//
	float nEvent = GetEvCountCycleRaws(itc);   // counted events for given trigger class
	if(nEvent>0) { 
          Float_t numberOfEventsAllCycles = 0.0;
          if((hEventCounter=(TH1F*) GetRawsData(240,itc))){
            numberOfEventsAllCycles = hEventCounter->Integral() + nEvent;// count all events upto now
            hEventCounter->SetBinContent(1,numberOfEventsAllCycles); // increase counter 
          } 

          SetEfficiency(169, 241, itc, numberOfEventsAllCycles);
          SetEfficiency(207, 242, itc, numberOfEventsAllCycles);
          SetEfficiency(208, 243, itc, numberOfEventsAllCycles);
          SetEfficiency(237, 244, itc, numberOfEventsAllCycles);
          SetEfficiency(238, 245, itc, numberOfEventsAllCycles);


          //fk// orA and orC for given specie and trigger class
          Float_t  numberOfORAEndOfCycle = 0.0;
          Float_t  numberOfORCEndOfCycle = 0.0;
	  if((htmp=GetRawsData(172,itc))) numberOfORAEndOfCycle = htmp->Integral(); //ORA     
	  if((htmp=GetRawsData(173,itc))) numberOfORCEndOfCycle = htmp->Integral(); //ORC     
          
          if((heff=GetRawsData(209,itc))){ //QTC efficiency
            if((hcounter=GetRawsData(246,itc))){ //QTC counter
              if(numberOfORCEndOfCycle>0){
                for(int ipmt=0; ipmt<12; ipmt++){
                  Float_t val = hcounter->GetBinContent(ipmt+1); //first bin has consequtive number 1 
                  heff->SetBinContent(ipmt+1,val/numberOfORCEndOfCycle);
                                 
               }
              }else{
                for(int ipmt=0;ipmt<12; ipmt++)
                   heff->SetBinContent(ipmt+1,0);
              }
              if(numberOfORAEndOfCycle>0){
                for(int ipmt=12;ipmt<24; ipmt++){
                  Float_t val = hcounter->GetBinContent(ipmt+1);
                  heff->SetBinContent(ipmt+1,val/numberOfORAEndOfCycle);
                }
              }else{
                for(int ipmt=0;ipmt<12; ipmt++)
                  heff->SetBinContent(ipmt+1,0);
              }
            }
          }
	}//Evt >0
      } // kRAWS
    } // RS: loop over eventual clones per trigger class
  } // loop over species
  //
  AliQAChecker::Instance()->Run(AliQAv1::kT0, task, list); //FK
}


//____________________________________________________________________________
void AliT0QADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliCDBManager* man = AliCDBManager::Instance();
  //man->SetDefaultStorage(gSystem->Getenv("AMORE_CDB_URI"));
  if(!man) return; 
  AliCDBEntry* entry = man->Get("GRP/Calib/QAThresholds");
  if(!entry) return;
  TObjArray* t0branch = (TObjArray*) entry->GetObject();
  AliQAThresholds*  thresholds = (AliQAThresholds*) t0branch->FindObject("T00");
  // here you should test that you got a non-null pointer


if(!thresholds) return;
  if(AliDAQ::DetectorID("T0")!= thresholds->GetDetectorId()){
    AliInfo(Form("DETECTOR ID %d DOES NOT MATCH TO TZERO",thresholds->GetDetectorId()));
    return;
  }
  
  int iparam = 0; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fMeanRawVertexParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  }

  iparam = 76; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fMeanORAParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  } 
 
  iparam = 77; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fMeanORCParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  } 
 
  iparam = 78; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fCFDEffSubRangeLowParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  } 
 
  iparam = 79; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fCFDEffSubRangeHighParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  } 
    
  iparam = 80; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fLEDEffSubRangeLowParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  } 
    
  iparam = 81; 
  if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean raw vertex 
    fLEDEffSubRangeHighParam = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
  } 

  for(int ipmt=0; ipmt<24;ipmt++){ 
    iparam = ipmt + 1; //current consecutive number of parameter
    if((TParameter<float>*) thresholds->GetThreshold(iparam)){ // mean CFD from a good run 
      fMeanCFDFromGoodRunParam[ipmt] = ((TParameter<float>*) thresholds->GetThreshold(iparam))->GetVal();
    }
  }
}
 
//____________________________________________________________________________ 
void AliT0QADataMakerRec::InitRaws()
{

  // create Raw histograms in Raw subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Float_t low[500];
  Float_t high[500];
  //triggers
  const Char_t *triggers[6] = {"T0 OR", "vertex","ORA","ORC","central","semi-central"};
   
  
  for (Int_t i=0; i<500; i++){
    low[i] = 0;
    high[i] = 30000;

  }

  TString timename, ampname, qtcname, ledname;
  TString timeCalname, ampCalname, ledCalname, qtcCalname;
  TString qt1name, qt0name, qt1Calname, qt0Calname;
  TString nhits;

  TH1F* hRefPoint = new TH1F("hRefPoint","Ref Point", 10000, 0 ,50000);
  hRefPoint->SetLabelSize(0.02);
  Add2RawsList( hRefPoint,0, expert, !image, !saveCorr);

  TH1F *hRefPointcal = new TH1F("hRefPointcal","Ref Point laser", 5000, 0 ,20000);
  Add2RawsList( hRefPointcal,249, expert, !image, !saveCorr);

  TH1F *hRawCFD[24]; 
  TH1F *hRawLEDamp[24];
  TH1F *hRawQTC[24]; TH1F *hRawLED[24];
  TH1F *hRawQT1[24]; TH1F *hRawQT0[24];
  TH1F* hRawNhits[24];
  for(Int_t i=0; i<24; i++){
      timename ="CFD/hRawCFD";
      ledname = "LED/hRawLED";
      qtcname = "QTC/hRawQTC";
      qt0name = "QTC/start/hRawQT0_";
      qt1name = "QTC/stop/hRawQT1_";
      ampname = "LEDminCFD/hRawLEDminCFD";
      nhits = "Hits/hRawNhits";
      timename += i+1;
      ampname += i+1;
      qtcname += i+1;
      qt0name += i+1;
      qt1name += i+1;
      ledname += i+1;
      nhits   += i+1;
      
      hRawCFD[i] = new TH1F(timename.Data(), Form("%s;CFD [#channels]; Counts", timename.Data()),Int_t((high[i+1]-low[i+1])/4),low[i+1],high[i+1]);
      //      ForbidCloning(hRawCFD[i]);       //RS I don't know how histos 1-24 should be processed in MakeRaws, for the moment forbidding the cloning
      Add2RawsList( hRawCFD[i],i+1, expert, !image, !saveCorr);
      hRawLED[i] = new TH1F(ledname.Data(),  Form("%s;LED [#channels]; Counts", ledname.Data()),Int_t((high[i+25]-low[i+25])/4),low[i+25],high[i+25]);
      Add2RawsList( hRawLED[i],i+25, expert, !image, !saveCorr);
      hRawLEDamp[i] = new TH1F(ampname.Data(),  Form("%s;LED-CFD [#channels]; Counts", ampname.Data()),1000,0,1000);
      Add2RawsList( hRawLEDamp[i],i+49, expert, !image, !saveCorr);
      hRawQTC[i] = new TH1F(qtcname.Data(),  Form("%s;QTC[#channels]; Counts", qtcname.Data()), 2500,0,10000); //fk
      //QT0
      Add2RawsList( hRawQTC[i],i+73, expert, !image, !saveCorr);
      hRawQT0[i] = new TH1F(qt0name.Data(),  Form("%s; QT0 [#channels]; Counts", qt0name.Data()),Int_t((high[97+i]-low[97+i])/4),low[97+i],high[97+i]);
      Add2RawsList( hRawQT0[i],97+i, expert, !image, !saveCorr);
      //QT1
      hRawQT1[i] = new TH1F(qt1name.Data(),  Form("%s; QT1 [#channels]; Counts", qt1name.Data()),Int_t((high[121+i]-low[121+i])/4),low[121+i],high[121+i]);
      Add2RawsList( hRawQT1[i],121+i, expert, !image, !saveCorr);
      
      hRawNhits[i] = new TH1F(nhits.Data(),  Form("%s;#Hits;Events", nhits.Data()),20, 0, 20);
      Add2RawsList( hRawNhits[i],176+i, expert, !image, !saveCorr);
    }
      //new QTC
  Int_t ihist=0;
  TH1F* hallhist[220];
  TString namech[4]=   {"00", "01", "10", "11"};
  TString namehist;
  for (Int_t i=0; i<12; i++)
    {      
      for (Int_t ih=0; ih<4; ih++) {
	 namehist = Form("newQT/StartStop/hnewRawQT%s_C%i",namech[ih].Data(),i+1);
	 hallhist[ihist]=new TH1F(namehist.Data(),  Form("%s; #channels;Events",namehist.Data()),1000, 0, 30000);
	 Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
	 ihist++;
       }
    }    
  for (Int_t i=12; i<24; i++)
    {      
      for (Int_t ih=0; ih<4; ih++) {
	namehist = Form("newQT/StartStop/hnewRawQT%s_A%i",namech[ih].Data(),i+1-12);
	 hallhist[ihist]=new TH1F(namehist.Data(),  Form("%s; #channels;Events",namehist.Data()),1000, 0, 30000);
	 Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
	 ihist++;
      }
    }    

  for (Int_t i=0; i<24; i++)
    {      
      hallhist[ihist] = new TH1F(Form("newQT/hnewRawQTC0_%i_diff",i+1),  Form("hRawQTC new %s - %s ch %i ;#channels;Events",namech[0].Data(), namech[1].Data(),i+1),1200, -100, 1100);
 	Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
	ihist++;
    }
  for (Int_t i=0; i<24; i++)
    {      
	hallhist[ihist] = new TH1F(Form("newQT/hnewRawQTC1_%i_diff",i+1),  Form("hRawQTC new %s - %s ch  %i ;#channels;Events",namech[2].Data(), namech[3].Data(),i+1),1200, -100, 1100);
	Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
	ihist++;
    }
  
  // new mult QTC
  TString namediff[4] = {"C_00min01","C_10min11", "A_00min01", "A_10min11"};  
  for (Int_t i=0; i<4; i++) {
    hallhist[ihist]  = new TH1F(Form("newMPD/hnewRawMultC_%s",namech[i].Data()),  Form("new C sum mult %s; #channels;Events",namech[i].Data()), 1000, 0, 30000);
    Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
    ihist++;
  }
  for (Int_t i=0; i<4; i++) {
    hallhist[ihist] = new TH1F(Form("newMPD/hnewRawMultA_%s",namech[i].Data()), Form("new A sum mult %s; #channels;Events",namech[i].Data()), 1000, 0, 30000);
    Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
    ihist++;
  }
  
  for (Int_t i=0; i<4; i++) {
    TString namempd=Form("newMPD/hnewRawMPD_%s_diff",namediff[i].Data());
    hallhist[ihist] =  new TH1F ( namempd.Data(),  namempd.Data(),
				  1200, -100, 1100) ; 
    Add2RawsList(hallhist[ihist],250+ihist, expert, !image, !saveCorr);
    ihist++;
  }
  // end new QTC
 
  TH1F* hRawTrigger = new TH1F("hRawTrigger"," triggers;Trigger ;Counts",6,0,6);
  for (Int_t itr=0; itr<6; itr++) hRawTrigger->Fill(triggers[itr], 0); // RS Modified to allow cloning (no fNumTriggers member anymore)
  Add2RawsList(hRawTrigger ,169, !expert, image, !saveCorr);
  TH1F* hRawMean = new TH1F("Triggers/hRawMean","online timer mean signal, physics event;",Int_t((high[170]-low[170])/4),low[170],high[170]);
  Add2RawsList( hRawMean,170, expert, !image, !saveCorr);

  TH1F* hRawVertex = new TH1F("Triggers/hRawVertex","online 0TVX vertex signal; counts",Int_t((high[171]-low[171])/4),low[171],high[171]);
  Add2RawsList( hRawVertex,171, expert, !image, !saveCorr);//FK

  TH1F* hRawORA = new TH1F("Triggers/hRawORA","online OR A; counts",Int_t((high[172]-low[172])/4),low[172],high[172]);
  Add2RawsList( hRawORA,172, expert, !image, !saveCorr);
  TH1F* hRawORC = new TH1F("Triggers/hRawORC","online OR C;counts",Int_t(( high[173]-low[173])/4),low[173],high[173]);
  Add2RawsList( hRawORC,173, expert, !image, !saveCorr);
  TH1F* hMultCentr = new TH1F("Triggers/hMultCentr","online trigger Central;counts ",Int_t(( high[174]-low[174])/4),low[174],high[174]);
  Add2RawsList( hMultCentr,174, expert, !image, !saveCorr);
  TH1F* hMultSeCentr = new TH1F("Triggers/hMultSemiCentr","online trigger SemiCentral;counts ",Int_t(( high[175]-low[175])/4),low[175],high[175]);
  Add2RawsList( hMultSeCentr,175, expert, !image, !saveCorr);

  TH1F* hMultA = new TH1F("Triggers/hMultA","full mulltiplicity A side;Multiplicity;Entries", Int_t((high[201]-low[201])/4) ,low[201],high[201]);
  Add2RawsList( hMultA,201, expert, !image, !saveCorr );//FK
  
  TH1F* hMultAS = new TH1F("Triggers/hMultASemi","full multiplicity with semi-central trigger A side ;Multiplicity;Entries",
			    Int_t((high[202]-low[202])/4),low[202],high[202] );
  Add2RawsList( hMultAS, 202, expert, !image, !saveCorr);
  TH1F* hMultAC = new TH1F("Triggers/hMultACentr","full multiplicity with central trigger;Multiplicity;Entries", 
			    Int_t((high[203]-low[203])/4),low[203],high[203]);
  Add2RawsList( hMultAC, 203, expert, !image, !saveCorr);
  
  
  //side C
   TH1F* hMultC = new TH1F("Triggers/hMultC","full mulltiplicity C side;Multiplicity;Entries", Int_t(high[204]-low[204]/4) ,low[204],high[204]);
  Add2RawsList( hMultC,204, expert, !image, !saveCorr );//FK
  TH1F* hMultCS = new TH1F("Triggers/hMultCSemi","full multiplicity with semi-central trigger C side;Multiplicity;Entries",
			    Int_t((high[205]-low[205])/4),low[205],high[205] );
  Add2RawsList( hMultCS,205, expert, !image, !saveCorr);
  TH1F* hMultCC = new TH1F("Triggers/hMultCCentr","full multiplicity with central trigger C side;Multiplicity;Entries", 
			    Int_t((high[206]-low[206])/4),low[206],high[206]);
  Add2RawsList( hMultCC,206, expert, !image, !saveCorr);
  
  
  //efficiency
  TH1F* hCFDeff= new TH1F("hCFDeff"," CFD efficiency; #PMT; #CFD counts/nEvents",24, 0 ,24);  
  hCFDeff->SetMinimum(0);
  hCFDeff->SetMaximum(2);
  hCFDeff->SetMarkerStyle(20);//fk
  hCFDeff->SetMarkerColor(2);//fk
  hCFDeff->SetOption("p");//fk
  Add2RawsList( hCFDeff, 207, expert, image, !saveCorr);//FK   
  TH1F* hEffLED = new TH1F("hEffLED","LED efficiency; #PMT; #LED counts/nEvent",24, 0 ,24);
  hEffLED ->SetMinimum(0);
  hEffLED->SetMaximum(2);
  hEffLED->SetMarkerStyle(28);//fk
  hEffLED->SetMarkerColor(1);//fk
  hEffLED->SetOption("p,same");//fk
  Add2RawsList( hEffLED, 208, expert, !image, !saveCorr);//FK is published attahced to the CFD efficiency 
  
  TH1F* hEffQTC = new TH1F("hEffQTC","QTC efficiency; #PMT; QTC efficiency%s;",24, 0 ,24);
  hEffQTC->SetMinimum(0);
  hEffQTC->SetMaximum(2);
  Add2RawsList( hEffQTC,209, !expert, image, !saveCorr);
   
  TH2F* hCFD = new TH2F("hCFD","CFD ; #PMT; CFD {#channnels}", 24, 0 , 24,Int_t((high[210]-low[210])/4),low[210],high[210]);
  hCFD->SetOption("COLZ");
  Add2RawsList( hCFD,210, expert, !image, !saveCorr);//fk
    
  TH2F* hLED = new TH2F("hLED","LED ; #PMT; LED [#channnels]", 24, 0 , 24,Int_t((high[211]-low[211])/4),low[211],high[211]);
  hLED->SetOption("COLZ");
  Add2RawsList( hLED,211, expert, !image, !saveCorr);//fk

  TH2F* hQTC = new TH2F("hQTC","QTC ; #PMT; QTC [#channnels]", 24, 0, 24,Int_t( (high[212]-low[212])/4),low[212],high[212]);
  hQTC->SetOption("COLZ");
  Add2RawsList( hQTC,212, expert, !image, !saveCorr);//fk
  
  TH1F* hNumPMTA= new TH1F("hNumPMTA","number of PMT hitted per event A side",13, 0 ,13);
  Add2RawsList(hNumPMTA ,213, expert, image, !saveCorr);
  
  TH1F* hNumPMTC= new TH1F("hNumPMTC","number of PMT hitted per event C side",13, 0 ,13);
  Add2RawsList(hNumPMTC ,214, expert, image, !saveCorr);
  
  TH1F* hHitsOrA= new TH1F("hHitsOrA","T0_OR A hit multiplicity",20, 0 ,20);
  Add2RawsList( hHitsOrA,215, expert, !image, !saveCorr);
  
  TH1F* hHitsOrC= new TH1F("hHitsOrC","T0_OR C hit multiplicity",20, 0 ,20);
  Add2RawsList(hHitsOrC ,216, expert, !image, !saveCorr);
  
  
  TH1F* hOrCminOrA= new TH1F("Beam/hOrCminOrA","T0_OR C - T0_OR A [cm]",10000,-5000,5000);
  Add2RawsList( hOrCminOrA,219, expert, !image, !saveCorr); //FK

  TH1F* hOrCminOrATvdcOn= new TH1F("Beam/hOrCminOrATvdcOn","T0_OR C - T0_OR A TVDC on [cm]",10000,-5000,5000);
  Add2RawsList( hOrCminOrATvdcOn,217, expert, !image, !saveCorr);//FK
  

  TH1F* hOrCminOrATvdcOff= new TH1F("Beam/hOrCminOrATvdcOff","T0_OR C - T0_OR A TVDC off [cm]",10000,-5000,5000);
  Add2RawsList( hOrCminOrATvdcOff,218, expert, !image, !saveCorr);//FK

   //satellite  & beam background
  TH2F* hBeam = new TH2F("Beam/hBeam", "Mean vs Vertex from 1st hit", 120, -30, 30, 120, -30, 30);
  hBeam->SetOption("COLZ");
  hBeam->GetXaxis()->SetTitle("(T0C-T0A)/2, ns from 1st"); //vtx
  hBeam->GetYaxis()->SetTitle("(T0C+T0A)/2, ns"); //time
  Add2RawsList( hBeam,220, !expert, image, !saveCorr);

  TH2F* hBeamTVDCon = new TH2F("Beam/hBeamTVDCon", "Mean vs Vertex TVDC on from 1st hit",50, -5, 5, 50, -5, 5);//FK
  hBeamTVDCon->SetOption("COLZ");
  hBeamTVDCon->GetXaxis()->SetTitle("(T0C-T0A)/2, ns from 1st hit");
  hBeamTVDCon->GetYaxis()->SetTitle("(T0C+T0A)/2, ns");
  Add2RawsList( hBeamTVDCon,221, expert, image, !saveCorr);

  TH2F* hBeamTVDCoff = new TH2F("Beam/hBeamTVDCoff", "Mean vs Vertex TVDC off from 1st hit", 120, -30, 30, 120, -30, 30);
  hBeamTVDCoff->GetXaxis()->SetTitle("(T0C-T0A)/2, ns from 1st hit");
  hBeamTVDCoff->GetYaxis()->SetTitle("(T0C+T0A)/2, ns");
  hBeamTVDCoff->SetOption("COLZ");
  Add2RawsList( hBeamTVDCoff,222, expert, image, !saveCorr);

  //vertex 1st
  // TH1F* hVertex1stTVDCon = new TH1F("Beam/hVertex1stTVDCon", "(T0A-T0C)/2, ps, from 1st hit TVDC on", 200, -2000, 2000); //FK
  TH1F* hVertex1stTVDCon = new TH1F("Beam/hVertex1stTVDCon", "(T0A-T0C)/2, cm, from 1st hit TVDC on", 200, -100, 100); //alla
   Add2RawsList(hVertex1stTVDCon ,223, !expert, image, !saveCorr);
   //  TH1F* hVertex1stTVDCoff = new TH1F("Beam/hVertex1stTVDCoff", "(T0A-T0C)/2, ps, from 1st hit TVDC off", 500, -2000, 2000);//FK
  TH1F* hVertex1stTVDCoff = new TH1F("Beam/hVertex1stTVDCoff", "(T0A-T0C)/2, cm, from 1st hit TVDC off", 500, -100, 100);//alla
  Add2RawsList( hVertex1stTVDCoff,225, !expert, image, !saveCorr);
  TH1F* hMean1stTVDCon  = new TH1F("Beam/hMean1stTVDCon", "(T0A+T0C)/2, ps, from 1st hit TVDC on", 200, -2000, 2000);//FK
  Add2RawsList( hMean1stTVDCon,  226, !expert, image, !saveCorr);
  TH1F* hMean1stTVDCoff = new TH1F("Beam/hMean1stTVDCoff", "(T0A+T0C)/2, ps, from 1st hit TVDC off", 200, -2000, 2000);//FK
  Add2RawsList( hMean1stTVDCoff, 227, !expert, image, !saveCorr);

   
  //FK histograms start from 230
  TH1F* hRawVertexMinMean = new TH1F("hRawVertexMinMean","online 0TVX vertex signal minus mean; channels",200,-200,200);
  Add2RawsList(hRawVertexMinMean,230, expert, image, !saveCorr);//FK

  TH1F* hCFDSubtrMean = new TH1F("hCFDSubtrMean","CFD minus mean; #PMT; CFD - mean {#channnels}", 24, 0, 24);
  hCFDSubtrMean->SetMarkerStyle(20);
  hCFDSubtrMean->SetOption("p");
  Add2RawsList( hCFDSubtrMean,231, !expert, image, !saveCorr);//fk filled in Checker
    
  TH1F* hLEDSubtrMean = new TH1F("hLEDSubtrMean","LED minus mean; #PMT; LED - mean [#channnels]", 24, 0, 24);
  hLEDSubtrMean->SetMarkerStyle(20);
  hLEDSubtrMean->SetOption("p");
  Add2RawsList( hLEDSubtrMean,232, expert, image, !saveCorr);//fk filled in Checker

  TH1F* hQTCSubtrMean = new TH1F("hQTCSubtrMean","QTC minus mean; #PMT; QTC - mean [#channnels]", 24, 0, 24);
  hQTCSubtrMean->SetMarkerStyle(20);
  hQTCSubtrMean->SetOption("p");
  Add2RawsList( hQTCSubtrMean,233, expert, image, !saveCorr);//fk filled in Checker
 
  TH2F* hDiffOrCVersusDiffOrATvdcOn= new TH2F("Beam/hDiffOrCVersusDiffOrATvdcOn","ORC-meanORC versus ORA-meanORA (TVDC on)",50,-200,200,50,-200,200);
  hDiffOrCVersusDiffOrATvdcOn->SetOption("COLZ");
  hDiffOrCVersusDiffOrATvdcOn->GetXaxis()->SetTitle("ORA - mean ORA [channel]");
  hDiffOrCVersusDiffOrATvdcOn->GetYaxis()->SetTitle("ORC - mean ORC [channel]");
  Add2RawsList(hDiffOrCVersusDiffOrATvdcOn, 234, expert, image, !saveCorr);//FK
  
  TH2F* hDiffOrCVersusDiffOrATvdcOff= new TH2F("Beam/hDiffOrCVersusDiffOrATvdcOff","ORC-meanORC vetsus ORA-meanORA (TVDC off)",50,-200,200,50,-200,200);
  hDiffOrCVersusDiffOrATvdcOff->SetOption("COLZ");
  hDiffOrCVersusDiffOrATvdcOff->GetXaxis()->SetTitle("ORA - mean ORA [channel]");
  hDiffOrCVersusDiffOrATvdcOff->GetYaxis()->SetTitle("ORC - mean ORC [channel]");
  Add2RawsList(hDiffOrCVersusDiffOrATvdcOff, 235, expert, image, !saveCorr);//FK
 
  TH2F* hBCID = new TH2F("hBCID", "header BCID vs TRM BC ID ", 500, 0, 5000, 500, 0, 5000);
  hBCID->SetOption("COLZ");
  hBCID->GetXaxis()->SetTitle("TRM BC ID");
  hBCID->GetYaxis()->SetTitle("event header BC ID");
  Add2RawsList(hBCID ,236, !expert, image, !saveCorr);

  //CFD and LED efficiency in range ~2000- ~3000 
  TH1F* hCFDeffSubRange = new TH1F("hCFDeffSubRange"," CFD eff in subrange; #PMT; #CFD counts/nEvents",24, 0 ,24);  
  Add2RawsList( hCFDeffSubRange, 237, expert, !image, !saveCorr);//FK  

 
  TH1F* hEffLEDSubRange = new TH1F("hEffLEDSubRange","LED eff in subrange; #PMT; #LED counts/nEvent",24, 0 ,24);
  Add2RawsList( hEffLEDSubRange,238, expert, !image, !saveCorr);//FK
  // ratio CDF eff /LEF eff in subragne 
  TH1F* hRatioCFDLEDeff = new TH1F("hRatioCFDLEDeff","Ratio CFD/LED eff in subrange; #PMT; ratio CDF/LED eff",24, 0 ,24);  
  hRatioCFDLEDeff->SetMinimum(0);
  hRatioCFDLEDeff->SetMaximum(2);
  Add2RawsList( hRatioCFDLEDeff, 239, expert, image, !saveCorr);//FK   
 
  TH1F* hEventCounter = new TH1F("hEventCounter","Event counter for eff histos; X; number of events",1, 0 ,1);  
  Add2RawsList( hEventCounter, 240, expert, !image, !saveCorr);//FK   

  //counters 
  TH1F* hRawTriggerCounter = new TH1F("hRawTriggerCounter"," triggers;Trigger ;Counts",6,0,6);
  for (Int_t itr=0; itr<6; itr++) hRawTriggerCounter->Fill(triggers[itr], 0);
  Add2RawsList(hRawTriggerCounter ,241, expert, !image, !saveCorr);
 
  TH1F* hCFDCounter= new TH1F("hCFDCounter"," CFD counter #PMT; #CFD counts",24, 0 ,24);  
  Add2RawsList( hCFDCounter, 242, expert, !image, !saveCorr);//FK   
  TH1F* hLEDCounter = new TH1F("hLEDCounter","LED counter; #PMT; #LED counts",24, 0 ,24);
  Add2RawsList( hLEDCounter, 243, expert, !image, !saveCorr);//FK
 
  TH1F* hCFDeffSubRangeCounter = new TH1F("hCFDeffSubRangeCounter"," CFD eff in subrange counter; #PMT; #CFD counts",24, 0 ,24);  
  Add2RawsList( hCFDeffSubRangeCounter, 244, expert, !image, !saveCorr);//FK  
  TH1F* hEffLEDSubRangeCounter = new TH1F("hEffLEDSubRangeCounter","LED eff in subrange counter; #PMT; #LED counts",24, 0 ,24);
  Add2RawsList( hEffLEDSubRangeCounter,245, expert, !image, !saveCorr);//FK
 
  TH1F* hQTCCounter = new TH1F("hQTCCounter","QTC counter; #PMT; QTC counts;",24, 0 ,24);
  Add2RawsList( hQTCCounter,246, expert, !image, !saveCorr);
 
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}
  
//____________________________________________________________________________ 
void AliT0QADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH2F * fhDigCFD = new TH2F("fhDigCFD", " CFD digits; #PMT; CFD digits[#channels]",25,-0.5,24.5,100,0,1000);
  fhDigCFD->SetOption("COLZ");
  Add2DigitsList( fhDigCFD,0, !expert, image);
  TH2F *fhDigLEDamp = new TH2F("fhDigLEDamp", " LED-CFD digits; #PMT; LED-CFD amplitude ",25,-0.5,24.5,100,100,1000);
  fhDigLEDamp->SetOption("COLZ");
  Add2DigitsList( fhDigLEDamp,1, !expert, !image);
  TH2F * fhDigQTC = new TH2F("fhDigQTC", " QTC digits; #PMT; QTC amplitude",25,-0.5,24.5,100,100,10000);
  fhDigQTC->SetOption("COLZ");
  Add2DigitsList( fhDigQTC,2, !expert, !image);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 

void AliT0QADataMakerRec::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH2F* fhRecCFD = new TH2F("hRecCFD"," CFD time;#PMT; CFD Time [ns];",24, 0 ,24, 
			      100,-50,50);
  fhRecCFD->SetOption("COLZ");
  Add2RecPointsList ( fhRecCFD,0, !expert, image);

  TH2F* fhRecAmpDiff = new TH2F("hRecAmpDiff"," LED-CFD  min QTC amplitude;#PMT; difference [MIPs];",
				24, 0 ,24, 200,-10,10);
  fhRecAmpDiff->SetOption("COLZ");
  Add2RecPointsList (fhRecAmpDiff, 1, !expert, image);
  
  TH1F *fhMean = new TH1F("hMean","online - rec mean;online - rec mean[#channels];",2000, -1000, 1000);
  Add2RecPointsList ( fhMean,2, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//____________________________________________________________________________
void AliT0QADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *fhESDMean = new TH1F("hESDmean"," ESD mean; mean time[%channels]",1000, -5, 5);
  Add2ESDsList(fhESDMean, 0, expert, !image) ;
  TH1F * fhESDVertex = new TH1F("hESDvertex","ESDvertex; vertex[cm];",82,-30,30);
  Add2ESDsList(fhESDVertex, 1, expert, !image) ;
  
  TH1F * fhESDResolution = new TH1F("hESDResolution","(T0A-T0C)/2 corrected by SPD vertex; ns",800,-2,2);
  Add2ESDsList(fhESDResolution, 2, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRaws( AliRawReader* rawReader)
{
  //indices in lookup table lookUpTable_tanay.txt 
  enum { kTZeroRefPoint=0, kTZeroFirstCfdC=1, kTZeroFirstLedC=13, kTZeroFirstQT0C=25,kTZeroFirstQT1C=26,
         kTZeroVertex=50, kTZeroOrA=51, kTZeroOrC=52, kT0multAQ0=53, kT0multAQ1=54, kTZeroMultCent=55, kTZeroMultSemi=56,
         kTZeroFirstCfdA=57, kTZeroFirstLedA=69, kTZeroFirstQT0A=81,kTZeroFirstQT1A=82,
         kT0multCQ0=105, kT0multCQ1=106, kT0meaner=49 
  }; 

  Int_t  time[24] ;
  for(Int_t i=0; i<24; i++) time[i] = 0;	  
  rawReader->Reset() ; 
  //fills QA histos for RAW
  //Int_t shift=0;
  // Int_t refPointParam = GetRecoParam()->GetRefPoint();
  Int_t refpoint = 0;
  Int_t refPointParam = 0;
  
  AliT0RawReader *start = new AliT0RawReader(rawReader);
  
  if (! start->Next()) {
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" no raw data found!!"));
    delete start;
    return;
  }
  UInt_t type =rawReader->GetType();
  if (GetEventSpecie()==AliRecoParam::kCalib && type!=8) {
     delete start;
    return;
  }
  //
  // RS: Don't use custom counters, they create problems with trigger cloning
  //     Use instead framework counters, incremented in the end of this routine
  // RS: There is some inconsistency here: the separation of physics and calib. events/histos is done by
  // fEventSpecie. Why do we book separate histos on different slots for calib and physics ? 
  // I am changing this in such way that we don't need local counters like fNumTriggers (the corresponding
  // histos now incremented in the MakeRaws, and for the normalization I will use the framework's counters
  // AliQADataMaker::GetEvCountCycle(...), AliQADataMaker::GetEvCountTotal(...)
  //
  // I think the histos xx+250 should be suppressed (the xx calib histos of specie==calibration will be 
  // used automatically)
      
 //
  //BC ID
  //  if (type == 7){
  UInt_t bcid = rawReader->GetBCID();
  UInt_t	trmbcid = start->GetTRMBunchID();
  
  FillRawsData(236,trmbcid, bcid);
  //FillRawsData(236,100, bcid); fake error for testing
 
 //  }    
  //    if (type == 7){ shift=1;   fnEventPhys++;}
  Int_t allData[220][5];
  for(Int_t i0=0; i0<220; i0++){
    for(Int_t j0=0; j0<5; j0++){
      allData[i0][j0]=0;
    } 
  }

  for(Int_t i=0; i<=211; i++){ 
    for(Int_t iHit=0; iHit<5; iHit++){
      allData[i][iHit]= start->GetData(i,iHit);
    }
  }

  if( allData[kTZeroRefPoint][0] > 0  /*&& (type == 7)*/){
    FillRawsData(0, allData[kTZeroRefPoint][0]); //Reference point
  }
  refpoint = allData[refPointParam][0];
  if(refPointParam <  0 ) refpoint=0; 
  if(refPointParam == 0 ) refpoint = allData[kTZeroRefPoint][0] - 5000;
  
  Int_t offsetCDF, offsetLED, offsetQT0, offsetQT1;
  Int_t numPmtC=0;    
  Int_t numPmtA=0;   
 
  for(Int_t ik = 0; ik<24; ik++){
    Int_t ipmt = ik; //C side
    if(ik<12) {
      offsetCDF = kTZeroFirstCfdC;
      offsetLED = kTZeroFirstLedC;
      offsetQT0 = kTZeroFirstQT0C;
      offsetQT1 = kTZeroFirstQT1C;
      if(allData[ipmt+offsetCDF][0]>0 /*&& type == 7 */ )  numPmtC++;
    }else{
      ipmt = ik - 12; //A side 
      offsetCDF = kTZeroFirstCfdA;
      offsetLED = kTZeroFirstLedA;
      offsetQT0 = kTZeroFirstQT0A;
      offsetQT1 = kTZeroFirstQT1A;
      if(allData[ipmt + offsetCDF][0]>0 /*&& type == 7 */) numPmtA++;
    }
    Int_t nhitsPMT=0; //count hits for this pmt
    Bool_t  tvdcon=kFALSE;
    Bool_t orcon=kFALSE;
    Bool_t oraon=kFALSE;

    for (Int_t iHt=0; iHt<5; iHt++) {
      //cfd
      if(allData[ipmt+offsetCDF][iHt]>0){
	FillRawsData(ik+1, allData[ipmt+offsetCDF][iHt]);  //CFD for each PMT
	FillRawsData(210, ik, allData[ipmt+offsetCDF][iHt]); //CFD vs PMT
	FillRawsData(242,ik,1.); // CFD counter for efficiency  
	if( fCFDEffSubRangeLowParam<allData[ipmt+offsetCDF][iHt] && allData[ipmt+offsetCDF][iHt]<fCFDEffSubRangeHighParam){
          FillRawsData(244,ik,1.); //count CDF entries in given subrange  for   CDF/LED eff ratio
        }
	AliDebug(50,Form("%i CFD %i  data %s",ik, ipmt+offsetCDF,  GetRawsData(ik+1)->GetName()));
	nhitsPMT++;
      }
      //led
      if(allData[ipmt+offsetLED][iHt] > 0){ 
	FillRawsData(ik+25,allData[ipmt+offsetLED][iHt]);
	FillRawsData(211,ik, allData[ipmt+offsetLED][iHt]);
	FillRawsData(243,ik,1.); //LED counter for LED efficiency 
        if(fLEDEffSubRangeLowParam < allData[ipmt+offsetLED][iHt] && allData[ipmt+offsetLED][iHt]<fLEDEffSubRangeHighParam){
	  FillRawsData(245,ik,1.); //count LED entries in given subrange for   CDF/LED eff ratio
        } 
	AliDebug(50,Form("%i LED %i  data %s",ik, ipmt+offsetLED,  GetRawsData(ik+25)->GetName()));
      }

      //led -cfd
      if(allData[ipmt+offsetLED][iHt] > 0 && allData[ipmt+offsetCDF][iHt] > 0 )
	FillRawsData(ik+49, allData[ipmt+offsetLED][iHt]-allData[ipmt+offsetCDF][iHt]);
      
      //qtc
      if(allData[2*ipmt+offsetQT0][iHt] > 0 &&
	 allData[2*ipmt+offsetQT1][iHt] > 0 ) {
        
	FillRawsData(ik+73, allData[2*ipmt+offsetQT0][iHt]-allData[2*ipmt+offsetQT1][iHt]); //QTC = QT0 - QT1 for each channel
	FillRawsData(212,ik, allData[2*ipmt+offsetQT0][iHt]-allData[2*ipmt+offsetQT1][iHt]); //QTC vs pmt

	FillRawsData(246,ik,1.); // QTC counter for QTC efficiency

	AliDebug(50,Form("%i QTC %i  data %s",ik, 2*ipmt+offsetQT0, GetRawsData(ik+73)->GetName()));
	
      }
      if(allData[2*ipmt+offsetQT0][iHt] > 0) { //QT0
	AliDebug(50,Form("%i QT0 %i  data %s",ik, 2*ipmt+offsetQT0, GetRawsData(ik+97)->GetName()));
	FillRawsData(ik+97,allData[2*ipmt+offsetQT0][iHt]);
      }
      if(allData[2*ipmt+offsetQT1][iHt] > 0) {//QT1
	AliDebug(50,Form("%i QT1 %i  data %s",ik, 2*ipmt+offsetQT1, GetRawsData(ik+121)->GetName()));
	FillRawsData(ik+121,allData[2*ipmt+offsetQT1][iHt]);
      }
    }
      
    FillRawsData(ik+176, nhitsPMT);
  }
  FillRawsData(213, numPmtA);
  FillRawsData(214, numPmtC);
     
 
  Int_t trChannel[6] = {kT0meaner, kTZeroVertex, kTZeroOrA, kTZeroOrC, kTZeroMultCent, kTZeroMultSemi};
  Float_t ch2cm = 24.4*0.029979;     
  Int_t nhitsOrA=0;
  Int_t nhitsOrC=0;
  
  for (Int_t iHt = 0; iHt < 5; iHt++) {
     for (Int_t itr = 0; itr < 6; itr++) { //T0_MEAN,TO_VERTX,ORA,ORC,T0_mult,T0_mult
       if (allData[trChannel[itr]][iHt] > 0) {
	 // FillWeighedRawsData(169 + shift, itr, 1.); //hRawTrigger RS: increment counters
	 FillRawsData(241, itr, 1.); // fill trigger counter
	 // printf(" triggers %i  data %i\n", itr, iHt);
	 FillRawsData(170 + itr, allData[trChannel[itr]][iHt]);
	 
	 if (trChannel[itr] == kTZeroVertex) //T0_VERTEX minus mean from config files
	   FillRawsData(230, allData[kTZeroVertex][iHt] - fMeanRawVertexParam);
       }
     }
  }
  // ORC-mean   ORA -mean  //Alla 
  Bool_t  tvdcon=kFALSE;
  Bool_t orcon=kFALSE;
  Bool_t oraon=kFALSE;
  Int_t orAch, orCch; 
  Float_t  diffORA=-999999,  diffORC=-999999;
  for (Int_t iHt = 0; iHt < 5; iHt++) {
    if (allData[kTZeroOrA][iHt] > fMeanORAParam-400 && 
	allData[kTZeroOrA][iHt] < fMeanORAParam+400) {
      diffORA =  allData[kTZeroOrA][iHt] - fMeanORAParam;
      orAch = allData[kTZeroOrA][iHt];
      oraon=kTRUE;
      nhitsOrA++;
    }
    if ( allData[kTZeroOrC][iHt] > fMeanORCParam-400 && 
	 allData[kTZeroOrC][iHt] < fMeanORCParam+400 ) { 
	diffORC =  allData[kTZeroOrC][iHt] - fMeanORCParam;
	orCch = allData[kTZeroOrC][iHt];
	orcon=kTRUE; 
	nhitsOrC++;
    }
    if (allData[kTZeroVertex][iHt]>fMeanRawVertexParam-400 &&
	allData[kTZeroVertex][iHt]<fMeanRawVertexParam+400 ) tvdcon=kTRUE;
  }
  if (oraon&&orcon) {
    FillRawsData(219, (orCch - orAch) * ch2cm);
    if(tvdcon)   {//TVDC on
      FillRawsData(234, diffORA, diffORC);     
      FillRawsData(217, (orCch - orAch) * ch2cm);
    }
    else  {//TVDC off 
      FillRawsData(235, diffORA, diffORC);
      FillRawsData(218, (orCch - orAch) * ch2cm);
    }
  }
  //mult trigger signals phys
  //A side
  for (Int_t iHt=0; iHt<5; iHt++) {
    if(allData[kT0multAQ0][iHt]>0 && allData[kT0multAQ1][iHt]>0) {
      AliDebug(50,Form(" mpdA %i  data %s", 201,  GetRawsData(201)->GetName()));
      
      FillRawsData(201,allData[kT0multAQ0][iHt]-allData[kT0multAQ1][iHt]);
      if(allData[kTZeroMultSemi][iHt]>0) FillRawsData(202,allData[kT0multAQ0][iHt]-allData[kT0multAQ1][iHt]);
      if(allData[kTZeroMultCent][iHt]>0) FillRawsData(203,allData[kT0multAQ0][iHt]-allData[kT0multAQ1][iHt]);
    }
    
    //C side 
    if(allData[kT0multCQ0][iHt]>0 && allData[kT0multCQ1][iHt]>0) {
      AliDebug(50,Form(" mpdC %i  data %s", 204,  GetRawsData(204)->GetName()));
      
      FillRawsData(204,allData[kT0multCQ0][iHt]-allData[kT0multCQ1][iHt]);
      if(allData[kTZeroMultSemi][iHt]>0) FillRawsData(205,allData[kT0multCQ0][iHt]-allData[kT0multCQ1][iHt]);
      if(allData[kTZeroMultCent][iHt]>0) FillRawsData(206,allData[kT0multCQ0][iHt]-allData[kT0multCQ1][iHt]);
    }
  }
  
    FillRawsData(215,nhitsOrA);
    FillRawsData(216,nhitsOrC);

  // new QTC 
  Float_t diff[4];
  Int_t pmt;
  for(Int_t iHt = 0; iHt<5; iHt++) {
    for(int id=0; id<4; id++) diff[id] = 0;
    //new QTC C side
    for (Int_t  ik=0; ik<56; ik++)
      {
	if(ik<48) {
	  pmt=ik/4;    
	  if  (allData[107+ik][iHt]!=0) 
	    FillRawsData(ik+250, allData[107+ik][iHt]);	  
	}
	else
	  if  (allData[107+ik][iHt]!=0) FillRawsData(250+ik+144-48, allData[107+ik][iHt]);
      }
    
    for (Int_t ik=0; ik<48; ik+=4)
      {
	pmt=ik/4;    
	diff[0]=allData[107+pmt*4][iHt] - allData[107+pmt*4+1][iHt];
	diff[1]=allData[107+pmt*4+2][iHt] - allData[107+pmt*4+3][iHt];
	if(diff[0] != 0)   FillRawsData(250+pmt+96, diff[0]);
	if(diff[1] != 0)   FillRawsData(250+pmt+120,  diff[1]); //!!!
	
      }
    //new MPD ch 48+
    diff[0] = allData[107+48][iHt] - allData[107+48+1][iHt];
    diff[1] = allData[107+48+2][iHt] - allData[107+48+3][iHt];
    diff[2] = allData[107+48+4][iHt] - allData[107+48+5][iHt];
    diff[3] = allData[107+48+6][iHt] - allData[107+48+7][iHt];
    for (Int_t  i=0; i<4; i++) 
      if (diff[i] !=0) FillRawsData(250+152+i, diff[i]);
    
    //new QTC A
    for (Int_t  ik=56; ik<106; ik++)
      {
	pmt=(ik-8)/4;
	if  (allData[107+ik][iHt]!=0) {
	  FillRawsData(250+ik-8, allData[107+ik][iHt]);	  
	}
      }
    for (Int_t ik=56; ik<106; ik+=4)
      {
	pmt=(ik-8)/4;
	diff[0]=allData[107+ik][iHt] - allData[107+ik+1][iHt];
	diff[1]=allData[107+ik+2][iHt] - allData[107+ik+3][iHt];
	if(diff[0] != 0 ) {
	  FillRawsData(250+pmt+96, diff[0]);
	}
	if(diff[1] != 0 ) {
	  FillRawsData(250+pmt+120, diff[1]); //!!!
	}
      } 
  } //iHit
  //end new QTC
  //draw satellite
  for (int itr=-1;itr<GetNEventTrigClasses();itr++) { //RS loop over all active trigger classes, including the global one
    int itrID = itr==-1 ? -1 : int( GetEventTrigClass(itr)->GetUniqueID());

    Float_t c = 0.0299792458; // cm/ps
    TH2 *hBeam             = (TH2*)GetRawsData(220,itrID);
    TH2 *hBeamTVDCon       = (TH2*)GetRawsData(221,itrID);
    TH2 *hBeamTVDCoff      = (TH2*)GetRawsData(222,itrID);
    TH1 *hVertex1stTVDCon  = (TH1*)GetRawsData(223,itrID);
    TH1 *hVertex1stTVDCoff = (TH1*)GetRawsData(225,itrID);
    TH1 *hMean1stTVDCon    = (TH1*)GetRawsData(226,itrID);
    TH1 *hMean1stTVDCoff   = (TH1*)GetRawsData(227,itrID);
    if(hBeam || hBeamTVDCon || hBeamTVDCoff || hVertex1stTVDCon || hVertex1stTVDCoff || hMean1stTVDCon || hMean1stTVDCoff){
     
      Int_t time1stA=9999999;
      Int_t time1stC=9999999;
      for(Int_t ipmt=0; ipmt<12; ipmt++){
	if(allData[ipmt+kTZeroFirstCfdC][0] > 1 ) {
	  time[ipmt] = allData[ipmt+kTZeroFirstCfdC][0] - (Int_t) fMeanCFDFromGoodRunParam[ipmt]; //fk
	  if(time[ipmt] < time1stC)  time1stC=time[ipmt]; //timeC
	}
      }
      for( Int_t ipmt=12; ipmt<24; ipmt++){
	if(allData[ipmt-12+kTZeroFirstCfdA][0] > 0) {
	  time[ipmt] = allData[ipmt-12+kTZeroFirstCfdA][0] - (Int_t) fMeanCFDFromGoodRunParam[ipmt];//fk
	  if(time[ipmt] < time1stA)  time1stA=time[ipmt]; //timeC
	}
      }
      if(time1stA<99999 && time1stC< 99999) { //From First
	Float_t t01st  = 24.4 * (Float_t) (( time1stA + time1stC)/2.0);
	Float_t ver1st = 24.4 * (Float_t) (( time1stC - time1stA)/2.0);
	if(hBeam) hBeam->Fill(0.001*ver1st, 0.001*(t01st)); //Mean versus vertex
	if(allData[kTZeroVertex][0] > 0){//TVDC on
          if(hBeamTVDCon)      hBeamTVDCon->Fill(0.001*ver1st, 0.001*(t01st));//Mean versus  TVDC on from first
          //if(hVertex1stTVDCon) hVertex1stTVDCon->Fill(ver1st);
	  if(hVertex1stTVDCon) hVertex1stTVDCon->Fill(c*ver1st);   //alla ps2cm
          if(hMean1stTVDCon)   hMean1stTVDCon->Fill(t01st);
       }else{//TVDC off
          if(hBeamTVDCoff)      hBeamTVDCoff->Fill(0.001*ver1st, 0.001*(t01st));//FK// TVDC off from first
          if(hVertex1stTVDCoff) hVertex1stTVDCoff->Fill(c*ver1st); //alla ps2cm
          if(hMean1stTVDCoff)   hMean1stTVDCoff->Fill(t01st);
        }
      }	  
    } //
  } // RS loop over all active trigger classes, including the global one
    //
    IncEvCountCycleRaws();
    IncEvCountTotalRaws();
    //
    delete start;
}
  
//____________________________________________________________________________
void AliT0QADataMakerRec::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits

  TArrayI *digCFD = new TArrayI(24);
  TArrayI *digLED = new TArrayI(24);
  TArrayI *digQT0 = new TArrayI(24);
  TArrayI *digQT1 = new TArrayI(24);
  Int_t refpoint=0;
  
  TBranch *brDigits=digitsTree->GetBranch("T0");
  AliT0digit *fDigits = new AliT0digit() ;
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError(Form("EXEC Branch T0 digits not found"));
    delete digCFD;
    delete digLED;
    delete digQT0;
    delete digQT1;
    return;
  }
  digitsTree->GetEvent(0);
  digitsTree->GetEntry(0);
  brDigits->GetEntry(0);
  fDigits->GetTimeCFD(*digCFD);
  fDigits->GetTimeLED(*digLED);
  fDigits->GetQT0(*digQT0);
  fDigits->GetQT1(*digQT1);
  refpoint = fDigits->RefPoint();
  for (Int_t i=0; i<24; i++) {
    if (digCFD->At(i)>0) {
      Int_t cfd=digCFD->At(i)- refpoint;
      FillDigitsData(0,i,cfd);
      FillDigitsData(1,i, (digLED->At(i) - digCFD->At(i)));
      FillDigitsData(2,i, (digQT1->At(i) - digQT0->At(i)));
    }
  }  
  
  delete digCFD;
  delete digLED;
  delete digQT0;
  delete digQT1;
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits(); 
  delete fDigits;
  //  
}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  //fills QA histos for clusters

  AliT0RecPoint* frecpoints= new AliT0RecPoint ();
  if (!frecpoints) {
    AliError(":MakeRecPoints >> no recpoints found");
    return;
  }
  TBranch *brRec =clustersTree ->GetBranch("T0");
  if (brRec) {
    brRec->SetAddress(&frecpoints);
  }else{
    AliError(Form("EXEC Branch T0 rec not found "));
    return;
  } 

  brRec->GetEntry(0);
  
  for ( Int_t i=0; i<24; i++) {
    if(i<12)
      FillRecPointsData(0, i, frecpoints -> GetTime(i) - frecpoints -> GetTime(0)); 
    if(i>11)
      FillRecPointsData(0, i,  frecpoints -> GetTime(i) - frecpoints -> GetTime(12)); 
    FillRecPointsData(1, i, frecpoints -> GetAmp(i) - frecpoints->AmpLED(i));
  }
  Double_t mmm=frecpoints->GetOnlineMean()- frecpoints->GetMeanTime();
  FillRecPointsData(2,mmm);
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //  
  delete frecpoints;
}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD
  
  const Double32_t  *mean;
  mean = esd->GetT0TOF();
  Double32_t t0time= 0.001*mean[0];
  Double32_t orA= 0.001*mean[1];
  Double32_t orC=0.001* mean[2];

  if (t0time<99)   FillESDsData(0,t0time);
  if( esd->GetT0zVertex() <99) FillESDsData(1, esd->GetT0zVertex());
  if( orA<99 && orC<99) FillESDsData(2,(orA-orC)/2.);
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
}
//____________________________________________________________________________
//____________________________________________________________________________
void AliT0QADataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{

 //reset the detector histograms for a given task
  AliQADataMakerRec::ResetDetector(task);

  for(int ih=0; ih<=250; ih++){
    for(int itr=-1; itr < GetNEventTrigClasses(); itr++){ 
      int itrID = itr==-1 ? -1 : int( GetEventTrigClass(itr)->GetUniqueID());

      TH1 *htmp = (TH1*) GetRawsData(ih,itrID);
      if(htmp) htmp->Reset();
    } 
  }
}


/*
void AliT0QADataMakerRec::GetMeanAndSigma(TH1F* hist, Float_t &mean, Float_t &sigma) 
{

  const double window = 3.;  //fit window 
 
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  sigmaEstimate = hist->GetRMS();
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","RQ","Q");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

  delete fit;
}
*/

void AliT0QADataMakerRec::SetEfficiency(Int_t idxEffHisto, Int_t idxCounterHisto, Int_t trigger, Float_t totNumOfEvts){
  //calculate efficiency =  counts/number of events
  TH1* heff     = GetRawsData(idxEffHisto,trigger);
  TH1* hcounter = GetRawsData(idxCounterHisto,trigger);
  if(heff && hcounter && (totNumOfEvts>0.0)){
    int nb = heff->GetNbinsX();
    for(int ib=1;ib<=nb;ib++){
      heff->SetBinContent(ib,((Float_t) hcounter->GetBinContent(ib))/((Float_t) totNumOfEvts));
    }
  }
  return;
}
