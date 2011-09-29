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
#include "THnSparse.h"

#include "Riostream.h"
ClassImp(AliT0QADataMakerRec)
           
//____________________________________________________________________________ 
  AliT0QADataMakerRec::AliT0QADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kT0), 
		  "T0 Quality Assurance Data Maker")
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

  for(Int_t i=0; i<24; i++) fMeans[i]=2500;
}


//____________________________________________________________________________ 
AliT0QADataMakerRec::AliT0QADataMakerRec(const AliT0QADataMakerRec& qadm) :
  AliQADataMakerRec()

{
  //copy ctor 
 SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
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
  AliQAChecker::Instance()->Run(AliQAv1::kT0, task, list) ;
  //  const Char_t *triggers[6] = {"mean", "vertex","ORA","ORC","central","semi-central"};
  //
  TH1* htmp = 0;
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
	if (nEvent>0) {
	  if ( (htmp=GetRawsData(169,itc))     ) htmp->Scale(1./nEvent); // RS 
	  if ( (htmp=GetRawsData(207,itc))     ) htmp->Scale(1./nEvent); // RS  // former feffPhysC count. 
	  if ( (htmp=GetRawsData(208,itc)) ) htmp->Scale(1./nEvent); // RS  // former feffA count. Remove?
	  if ( (htmp=GetRawsData(209,itc)) ) htmp->Scale(1./nEvent); // RS  // former feffqtc count. 
	}
      } // kRAWS
    } // RS: loop over eventual clones per trigger class
  } // loop over species
  //
}


//____________________________________________________________________________
void AliT0QADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

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
  const Char_t *triggers[6] = {"mean", "vertex","ORA","ORC","central","semi-central"};
   
  
  for (Int_t i=0; i<500; i++){
    low[i] = 0;
    high[i] = 30000;

  }

  TString timename, ampname, qtcname, ledname;
  TString timeCalname, ampCalname, ledCalname, qtcCalname;
  TString qt1name, qt0name, qt1Calname, qt0Calname;
  TString nhits;

  TH1F* fhRefPoint = new TH1F("hRefPoint","Ref Point", 10000, 0 ,50000);
  fhRefPoint->SetLabelSize(0.02);
  Add2RawsList( fhRefPoint,0, expert, !image, !saveCorr);

  TH1F* fhRefPointcal = new TH1F("hRefPointcal","Ref Point laser", 5000, 0 ,20000);
  Add2RawsList( fhRefPointcal,250, expert, !image, !saveCorr);

  TH1F *fhRawCFD[24]; 
  TH1F * fhRawLEDamp[24];
  TH1F *fhRawQTC[24]; TH1F * fhRawLED[24];
  TH1F *fhRawQT1[24]; TH1F *fhRawQT0[24];
  TH1F* fhRawNhits[24];
 
  for (Int_t i=0; i<24; i++)
    {
      timename ="hRawCFD";
      ledname = "hRawLED";
      qtcname = "hRawQTC";
      qt0name = "hRawQT0_";
      qt1name = "hRawQT1_";
      ampname = "hRawLEDminCFD";
      nhits = "hRawNhits";
      timename += i+1;
      ampname += i+1;
      qtcname += i+1;
      qt0name += i+1;
      qt1name += i+1;
      ledname += i+1;
      nhits   += i+1;
      fhRawCFD[i] = new TH1F(timename.Data(), Form("%s;CFD [#channels];Counts", timename.Data()),Int_t((high[i+1]-low[i+1])/4),low[i+1],high[i+1]);
      //      ForbidCloning(fhRawCFD[i]);       //RS I don't know how histos 1-24 should be processed in MakeRaws, for the moment forbidding the cloning
      Add2RawsList( fhRawCFD[i],i+1, expert, !image, !saveCorr);
      fhRawLED[i] = new TH1F(ledname.Data(),  Form("%s;LED[#channels];Counts", ledname.Data()),Int_t((high[i+25]-low[i+25])/4),low[i+25],high[i+25]);
      Add2RawsList( fhRawLED[i],i+25, expert, !image, !saveCorr);
      fhRawLEDamp[i] = new TH1F(ampname.Data(),  Form("%s;LED-CFD [#channels];Counts", ampname.Data()),1000,0,1000);
      Add2RawsList( fhRawLEDamp[i],i+49, expert, !image, !saveCorr);
      fhRawQTC[i] = new TH1F(qtcname.Data(),  Form("%s;QTC[#channels];Counts", qtcname.Data()),10000,0,10000);
      Add2RawsList( fhRawQTC[i],i+73, expert, !image, !saveCorr);
      fhRawQT1[i] = new TH1F(qt1name.Data(),  Form("%s;QT1[#channels];Counts", qt1name.Data()),Int_t((high[97+i]-low[97+i])/4),low[97+i],high[97+i]);
      Add2RawsList( fhRawQT1[i],97+i, expert, !image, !saveCorr);
      fhRawQT0[i] = new TH1F(qt0name.Data(),  Form("%s;QT0[#channels];Counts", qt0name.Data()),Int_t((high[121+i]-low[121+i])/4),low[121+i],high[121+i]);
      Add2RawsList( fhRawQT0[i],121+i, expert, !image, !saveCorr);
      
      fhRawNhits[i] = new TH1F(nhits.Data(),  Form("%s;#Hits;Events", nhits.Data()),20, 0, 20);
      Add2RawsList( fhRawNhits[i],176+i, expert, !image, !saveCorr);
      
    }
  
      
  TH1F* fhRawTrigger = new TH1F("hRawTrigger"," triggers;Trigger ;Counts",6,0,6);
  for (Int_t itr=0; itr<6; itr++) fhRawTrigger->Fill(triggers[itr], 0); // RS Modified to allow cloning (no fNumTriggers member anymore)
  fhRawTrigger->SetMinimum(0);
  fhRawTrigger->SetMaximum(2);
   Add2RawsList(fhRawTrigger ,169, !expert, image, !saveCorr);
  
  TH1F* fhRawMean = new TH1F("hRawMean","online mean signal, physics event;",Int_t((high[170]-low[170])/4),low[170],high[170]);
  Add2RawsList( fhRawMean,170, expert, !image, !saveCorr);

  TH1F* fhRawVertex = new TH1F("hRawVertex","online vertex signal; counts",Int_t((high[171]-low[171])/4),low[171],high[171]);
  Add2RawsList( fhRawVertex,171, expert, image, !saveCorr);

  TH1F* fhRawORA = new TH1F("hRawORA","online OR A; counts",Int_t((high[172]-low[172])/4),low[172],high[172]);
  Add2RawsList( fhRawORA,172, expert, !image, !saveCorr);
  TH1F* fhRawORC = new TH1F("hRawORC","online OR C;counts",Int_t(( high[173]-low[173])/4),low[173],high[173]);
  Add2RawsList( fhRawORC,173, expert, !image, !saveCorr);
  TH1F* fhMultCentr = new TH1F("hMultCentr","online trigger Central;counts ",Int_t(( high[174]-low[174])/4),low[174],high[174]);
  Add2RawsList( fhMultCentr,174, expert, !image, !saveCorr);
  TH1F* fhMultSeCentr = new TH1F("hMultSemiCentr","online trigger SemiCentral;counts ",Int_t(( high[175]-low[175])/4),low[175],high[175]);
  Add2RawsList( fhMultSeCentr,175, expert, !image, !saveCorr);

  TH1F* fhMultA = new TH1F("hMultA","full mulltiplicity A side;Multiplicity;Entries", Int_t((high[201]-low[201])/4) ,low[201],high[201]);
  Add2RawsList( fhMultA,201, expert, image, !saveCorr );
  
  TH1F* fhMultAS = new TH1F("hMultASemi","full multiplicity with semi-central trigger A side ;Multiplicity;Entries",
			    Int_t((high[202]-low[202])/4),low[202],high[202] );
  Add2RawsList( fhMultAS, 202, expert, !image, !saveCorr);
  TH1F* fhMultAC = new TH1F("hMultACentr","full multiplicity with central trigger;Multiplicity;Entries", 
			    Int_t((high[203]-low[203])/4),low[203],high[203]);
  Add2RawsList( fhMultAC, 203, expert, !image, !saveCorr);
  
  
  //side C
   TH1F* fhMultC = new TH1F("hMultC","full mulltiplicity C side;Multiplicity;Entries", Int_t(high[204]-low[204]/4) ,low[204],high[204]);
  Add2RawsList( fhMultC,204, expert, image, !saveCorr );
  TH1F* fhMultCS = new TH1F("hMultCSemi","full multiplicity with semi-central trigger C side;Multiplicity;Entries",
			    Int_t((high[205]-low[205])/4),low[205],high[205] );
  Add2RawsList( fhMultCS,205, expert, !image, !saveCorr);
  TH1F* fhMultCC = new TH1F("hMultCCentr","full multiplicity with central trigger C side;Multiplicity;Entries", 
			    Int_t((high[206]-low[206])/4),low[206],high[206]);
  Add2RawsList( fhMultCC,206, expert, !image, !saveCorr);
  
  
  //efficiency
  TH1F* fhCFDeff= new TH1F("fhCFDeff"," CFD efficiency; #PMT; #CFD counts/nEvents",24, 0 ,24);  
  fhCFDeff->SetMinimum(0);
  fhCFDeff->SetMaximum(2);
  Add2RawsList( fhCFDeff, 207, !expert, image, !saveCorr);
  TH1F* fhEffLED = new TH1F("hEffLED","LEDefficiecy; #PMT; #LED counts/nEvent",24, 0 ,24);
  fhEffLED ->SetMinimum(0);
  fhEffLED->SetMaximum(2);
  Add2RawsList( fhEffLED,208, !expert, image, !saveCorr);
  
  TH1F* fhEffQTC = new TH1F("hEffQTC","QTC efficiecy; #PMT; QTC efficiency%s;",24, 0 ,24);
  fhEffQTC->SetMinimum(0);
  fhEffQTC->SetMaximum(2);
  Add2RawsList( fhEffQTC,209, !expert, image, !saveCorr);
   
  TH2F* fhCFD = new TH2F("hCFD","CFD ; #PMT; CFD {#channnels}",25, 0 ,25,Int_t((high[210]-low[210])/4),low[210],high[210]);
  fhCFD->SetOption("COLZ");
  Add2RawsList( fhCFD,210, expert, image, !saveCorr);
    
  TH2F* fhLED = new TH2F("hLED","LED ; #PMT; LED [#channnels]",25, 0 ,25,Int_t((high[211]-low[211])/4),low[211],high[211]);
  fhLED->SetOption("COLZ");
  Add2RawsList( fhLED,211, expert, image, !saveCorr);

  TH2F* fhQTC = new TH2F("hQTC","QTC ; #PMT; QTC [#channnels]",25, 0 ,25,Int_t( high[212]-low[212]),low[212],high[212]);
  fhQTC->SetOption("COLZ");
  Add2RawsList( fhQTC,212, expert, image, !saveCorr);
  
  TH1F* fhNumPMTA= new TH1F("hNumPMTA","number of PMT hitted per event A side",13, 0 ,13);
  Add2RawsList(fhNumPMTA ,213, expert, image, !saveCorr);
  
  TH1F* fhNumPMTC= new TH1F("hNumPMTC","number of PMT hitted per event C side",13, 0 ,13);
  Add2RawsList(fhNumPMTC ,214, expert, image, !saveCorr);
  
  TH1F* fhHitsOrA= new TH1F("fhHitsOrA","T0_OR A hit multiplicitie",20, 0 ,20);
  Add2RawsList( fhHitsOrA,215, expert, !image, !saveCorr);
  
  TH1F* fhHitsOrC= new TH1F("fhHitsOrC","T0_OR C hit multiplicitie",20, 0 ,20);
  Add2RawsList(fhHitsOrC ,216, expert, !image, !saveCorr);
  
  
  TH1F* fhOrCminOrA= new TH1F("fhOrCminOrA","T0_OR C - T0_OR A",10000,-5000,5000);
  Add2RawsList( fhOrCminOrA,219, expert, !image, !saveCorr);

  TH1F* fhOrCminOrATvdcOn= new TH1F("fhOrCminOrATvdcOn","T0_OR C - T0_OR A TVDC on",10000,-5000,5000);
  Add2RawsList( fhOrCminOrATvdcOn,217, expert, !image, !saveCorr);
  

  TH1F* fhOrCminOrATvdcOff= new TH1F("fhOrCminOrATvdcOff","T0_OR C - T0_OR A TVDC off",10000,-5000,5000);
  Add2RawsList( fhOrCminOrATvdcOff,218, expert, !image, !saveCorr);

   //satellite  & beam background
  TH2F* fhBeam = new TH2F("fhBeam", " Mean vs Vertex ", 120, -30, 30, 120, -30, 30);
  fhBeam->SetOption("COLZ");
  fhBeam->GetXaxis()->SetTitle("(T0C-T0A/2, ns");
  fhBeam->GetYaxis()->SetTitle("(T0C+T0A/2, ns");
  Add2RawsList( fhBeam,220, !expert, image, !saveCorr);
   TH2F* fhBeamTVDCon = new TH2F("fhBeamTVDCon", " Mean vs Vertex TVDC on ",50, -5, 5, 50, -5, 5) ;
   fhBeamTVDCon->SetOption("COLZ");
   fhBeamTVDCon->GetXaxis()->SetTitle("(T0C-T0A/2, ns");
   fhBeamTVDCon->GetYaxis()->SetTitle("(T0C+T0A/2, ns");
   Add2RawsList( fhBeamTVDCon,221, expert, image, !saveCorr);
   TH2F* fhBeamTVDCoff = new TH2F("fhBeamTVDCoff", " Mean vs Vertex TVDC off", 120, -30, 30, 120, -30, 30);
   fhBeamTVDCoff->GetXaxis()->SetTitle("(T0C-T0A/2, ns");
   fhBeamTVDCoff->GetYaxis()->SetTitle("(T0C+T0A/2, ns");
   fhBeamTVDCoff->SetOption("COLZ");
   Add2RawsList( fhBeamTVDCoff,222, expert, image, !saveCorr);
   TH1F* fhMean = new TH1F("fhMean", " (T0A+T0C)/2, ps ", 200, -2000, 2000);
   Add2RawsList( fhMean,223, !expert, image, !saveCorr);
   //vertex 1st 
   TH1F* fhVertex1st = new TH1F("fhVertex1st", " (T0A-T0C)/2, ps 1st particle", 200, -2000, 2000);
   Add2RawsList( fhVertex1st ,225, !expert, image, !saveCorr);
   TH1F* fhVertexBest = new TH1F("fhVertexBest", " (T0A-T0C)/2, ps , best particle", 200, -2000, 2000);
   Add2RawsList( fhVertexBest ,226, !expert, image, !saveCorr);
   TH1F* fhMean1st = new TH1F("fhMean1st", " (T0A + T0C)/2, ps , 1st particle", 200, -2000, 2000);
   Add2RawsList( fhMean1st ,227, expert, image, !saveCorr);

   // 
   TH2F* fhBCID = new TH2F("fhBCID", "header BCID vs TRM BC ID ", 500, 0, 5000, 500, 0, 5000);
   fhBCID->SetOption("COLZ");
   fhBCID->GetXaxis()->SetTitle("TRM BC ID");
   fhBCID->GetYaxis()->SetTitle("event header BC ID");
   Add2RawsList(fhBCID ,224, !expert, image, !saveCorr);

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
  Int_t  time[24] ;
  for(Int_t i=0; i<24; i++) time[i] = 0;	  
  
  rawReader->Reset() ; 
  //fills QA histos for RAW
  Int_t shift=0;
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
  FillRawsData(224,trmbcid, bcid);
  //  }    
  //    if (type == 7){ shift=1;   fnEventPhys++;}
  Int_t allData[110][5];
  for (Int_t i0=0; i0<110; i0++)
    {
      for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0;
    }
  for (Int_t i=0; i<107; i++) 
    for (Int_t iHit=0; iHit<5; iHit++)
      allData[i][iHit]= start->GetData(i,iHit);
  
  //if ( allData[0][0] > 0  && (type == 7))
    FillRawsData(0, allData[0][0]);
  refpoint = allData[refPointParam][0];
  if (refPointParam <  0 ) refpoint=0; 
  if (refPointParam == 0 ) refpoint = allData[0][0] - 5000;
  
  Int_t sideshift, sideshiftqtc;
  Int_t numPmtC=0;    
  Int_t numPmtA=0;    
  for (Int_t ik = 0; ik<24; ik++) {
    if(ik<12) {
      sideshift=1;
      sideshiftqtc=1;
      if(allData[ik+sideshift][0]>0 /*&& type == 7 */ )  numPmtC++;
    }
    else {
      if(allData[ik+45][0]>0 /*&& type == 7 */) numPmtA++;
      sideshift=45;
      sideshiftqtc=33;	
    }
    Int_t nhitsPMT=0;
    
    for (Int_t iHt=0; iHt<5; iHt++) {
      //cfd
      if(allData[ik+sideshift][iHt]>0) {
	FillRawsData(shift+ik+1, allData[ik+sideshift][iHt]);
	FillRawsData(210+shift, ik+1, allData[ik+sideshift][iHt]);
	FillRawsData(207+shift,ik,1.); // instead of incrementing former feff's, increment histo directly
	AliDebug(50,Form("%i CFD %i  data %s",ik, ik+sideshift,  GetRawsData(shift+ik+1)->GetName()));
	nhitsPMT++;
	  
      }
      //led
      if(allData[ik+12+sideshift][iHt] > 0) { 
	FillRawsData(shift+ik+24+1,allData[ik+12+sideshift][iHt]);
	FillRawsData(211+shift,ik+1, allData[ik+12+sideshift][iHt]);
	AliDebug(50,Form("%i LED %i  data %s",ik, ik+12+sideshift,  GetRawsData(shift+ik+1+24)->GetName()));
	FillRawsData(208+shift,ik,1.); // instead of incrementing former feff's, increment histo directly

      }
      //led -cfd
      
      if(allData[ik+12+sideshift][iHt] > 0 && allData[ik+sideshift][iHt] >0 )
	FillRawsData(shift+ik+48+1, allData[ik+12+sideshift][iHt]-allData[ik+sideshift][iHt]);
      
      //qtc
      if(allData[2*ik+sideshiftqtc+24][iHt] > 0 &&
	 allData[2*ik+sideshiftqtc+25][iHt] > 0) {
	FillRawsData(shift+ik+72+1, allData[2*ik+sideshiftqtc+24][iHt]-allData[2*ik+sideshiftqtc+25][iHt]);
	FillRawsData(212+shift,ik+1, allData[2*ik+sideshiftqtc+24][iHt]-allData[2*ik+sideshiftqtc+25][iHt]);

	FillRawsData(209+shift,ik,1.); // instead of incrementing former feff's, increment histo directly

	AliDebug(50,Form("%i QTC %i  data %s",ik, 2*ik+sideshiftqtc+24,  GetRawsData(shift+ik+1+72)->GetName()));
	
      }
      if(allData[2*ik+sideshiftqtc+24][iHt] > 0) {
	AliDebug(50,Form("%i QT0 %i  data %s",ik, 2*ik+sideshiftqtc+24,  GetRawsData(shift+ik+1+96)->GetName()));
	FillRawsData(shift+ik+96+1,allData[2*ik+sideshiftqtc+24][iHt]);
      }
      if(allData[2*ik+sideshiftqtc+25][iHt] > 0) {
	AliDebug(50,Form("%i QT0 %i  data %s",ik, 2*ik+sideshiftqtc+25,  GetRawsData(shift+ik+1+120)->GetName()));
	FillRawsData(shift+ik+120+1,allData[2*ik+sideshiftqtc+25][iHt]);
      }
    }
      
    FillRawsData(ik+176,nhitsPMT);
    FillRawsData(213,numPmtA);
    FillRawsData(214,numPmtC);
  }   
  
  Int_t trChannel[6] = {49,50,51,52,55,56};  
  Float_t ch2cm = 24.4*0.029979;     
  Int_t nhitsOrA=0;
  Int_t nhitsOrC=0;
  for (Int_t iHt=0; iHt<5; iHt++) {
    
    //orA-orC phys tvdc 1 
    if((allData[51][iHt]>0 && allData[52][iHt]>0) && allData[50][iHt]>0) {
      AliDebug(10,Form("orA-orC phys tvdc 1  %i  data %s", 217+shift,  GetRawsData(shift+217)->GetName()));
      
      FillRawsData(217+shift,(allData[52][iHt] - allData[51][iHt])*ch2cm);
    }
    //orA-orC phys tvdc 0 
    if((allData[51][iHt]>0 && allData[52][iHt]>0) && allData[50][iHt]<=0) {
      AliDebug(10,Form("orA-orC phys tvdc 0  %i  data %s", 218+shift,  GetRawsData(shift+218)->GetName()));
	
      FillRawsData(218+shift,(allData[52][iHt] - allData[51][iHt])*ch2cm);
    }
    if(allData[51][iHt]>0 && allData[52][iHt]>0) {
      AliDebug(50,Form("orA-orC phys tvdc all  %i  data %s", 219+shift,  GetRawsData(shift+219)->GetName()));
      FillRawsData(219+shift,(allData[52][iHt] - allData[51][iHt])*ch2cm);
    }
    for (Int_t itr=0; itr<6; itr++) {
      if (allData[trChannel[itr]][iHt] >0) {
	//
	// RS instead of incremented custom counters, fill directly the specie-specific histos
	//  	FillRawsData(169+shift, 0.5+itr, 1.);  // RS: increment counters
	FillRawsData(169+shift, itr, 1.);  // RS: increment counters
	AliDebug(50,Form(" triggers %i  data %s", 170+itr+shift,  GetRawsData(170+itr+shift)->GetName()));
	
	FillRawsData(170+itr+shift,allData[trChannel[itr]][iHt]);
      }
    }
    
    /*   if(type == 7) */if(allData[51][iHt] >0) nhitsOrA++;
    /* if(type == 7) */if(allData[52][iHt] >0) nhitsOrC++;
    
    //mult trigger signals phys
    //C side
    if(allData[53][iHt]>0 && allData[54][iHt]>0) {
      AliDebug(50,Form(" mpdA %i  data %s", 201+shift,  GetRawsData(201+shift)->GetName()));
      
      FillRawsData(201+shift,allData[53][iHt]-allData[54][iHt]);
      if(allData[56][iHt]>0) FillRawsData(202+shift,allData[53][iHt]-allData[54][iHt]);
      if(allData[55][iHt]>0) FillRawsData(203+shift,allData[53][iHt]-allData[54][iHt]);
    }
    
    //A side 
    if(allData[105][iHt]>0 && allData[106][iHt]>0) {
      AliDebug(50,Form(" mpdC %i  data %s", 204+shift,  GetRawsData(204+shift)->GetName()));
	
      FillRawsData(204+shift,allData[105][iHt]-allData[106][iHt]);
      if(allData[56][iHt]>0) FillRawsData(205+shift,allData[105][iHt]-allData[106][iHt]);
      if(allData[55][iHt]>0) FillRawsData(206+shift,allData[105][iHt]-allData[106][iHt]);
    }
  }
    
  FillRawsData(215,nhitsOrA);
  FillRawsData(216,nhitsOrC);
    
  //draw satellite
  for (int itr=-1;itr<GetNEventTrigClasses();itr++) { //RS loop over all active trigger classes, including the global one
    int itrID = itr==-1 ? -1 : GetEventTrigClass(itr)->GetUniqueID();
    int nEvent = GetEvCountCycleRaws(itrID);
    // 
    //   if(type == 7) {	 // RS Do we need here event type check, specie should do the job
    if(nEvent == 1000) {
      for (Int_t ik=0; ik<24; ik++) {
	//
	TH1* hik = (TH1*) GetRawsData(ik+1,itrID); if (!hik) continue; 
	hik->GetXaxis()->SetRangeUser(2000, 3000);
	int  maxBin  =  hik->GetMaximumBin(); 
	double   meanEstimate  =  hik->GetBinCenter( maxBin);
	TF1* fit= new TF1("fit","gaus", meanEstimate - 40, meanEstimate + 40);
	fit->SetParameters (hik->GetBinContent(maxBin), meanEstimate, 80);
	hik->Fit("fit","RQ","Q",  meanEstimate-40,  meanEstimate+40);
	fMeans[ik]= (Int_t) fit->GetParameter(1);
	hik->GetXaxis()->SetRangeUser(0, 30000);
	fit->Delete();
      }
    }
    //
    TH2 *h220 = (TH2*)GetRawsData(220,itrID);
    TH2 *h221 = (TH2*)GetRawsData(221,itrID);
    TH2 *h222 = (TH2*)GetRawsData(222,itrID);
    TH1 *h223 = (TH1*)GetRawsData(223,itrID);
    TH1 *h225 = (TH1*)GetRawsData(225,itrID);
    TH1 *h226 = (TH1*)GetRawsData(226,itrID);
    TH1 *h227 = (TH1*)GetRawsData(227,itrID);
    if( nEvent>1000 && (h220 || h221 || h222 || h223 || h225 || h226|| h227) ) 
      {
      Int_t besttimeA=9999999;
      Int_t besttimeC=9999999;
      Int_t time1stA=9999999;
      Int_t time1stC=9999999;
      for (Int_t ipmt=0; ipmt<12; ipmt++){
	if(allData[ipmt+1][0] > 1 ) {
	  //	      time[ipmt] = allData[ipmt+1][0] - Int_t(GetRawsData(1)->GetMean());
	  time[ipmt] = allData[ipmt+1][0] - fMeans[ipmt];
	  if(TMath::Abs(time[ipmt]) < TMath::Abs(besttimeC))
	    besttimeC=time[ipmt]; //timeC
	  if(time[ipmt] < time1stC)  time1stC=time[ipmt]; //timeC
	}
      }
      for ( Int_t ipmt=12; ipmt<24; ipmt++){
	if(allData[ipmt+45][0] > 0) {
	  time[ipmt] = allData[ipmt+45][0] - fMeans[ipmt] ;
	  if(TMath::Abs(time[ipmt]) < TMath::Abs(besttimeA) )
	    besttimeA=time[ipmt]; //timeA
	  if(time[ipmt] < time1stA)  time1stA=time[ipmt]; //timeC
	}
      }
      if(besttimeA<99999 &&besttimeC< 99999) {
	Float_t t0 =  24.4 * (Float_t( besttimeA+besttimeC)/2. );
	Float_t ver = 24.4 * Float_t( besttimeC-besttimeA)/2.;
	if (h220) h220->Fill(0.001*ver, 0.001*(t0));
	if(allData[50][0] > 0)  if (h221) h221->Fill(0.001*ver, 0.001*(t0));
	if(allData[50][0] <= 0) if (h222) h222->Fill(0.001*ver, 0.001*(t0));
	if (h223) h223->Fill(t0);
	if (h226) h226->Fill(ver);
      }	  
      if(time1stA<99999 &&time1stC< 99999) {
	Float_t t01st =  24.4 * (Float_t( time1stA + time1stC)/2. );
	Float_t ver1st = 24.4 * Float_t( time1stC - time1stA)/2.;
	if (h225) h225->Fill(ver1st);
	if (h227) h227->Fill(t01st);
      }	  
    } //event >100
      //   } //type 7
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
