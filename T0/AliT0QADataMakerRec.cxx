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
		  "T0 Quality Assurance Data Maker"),
  fnEventCal(0),
  fnEventPhys(0)

{
  // ctor
  for (Int_t i=0; i<6; i++) {
    fNumTriggers[i]=0;
    fNumTriggersCal[i]=0;
    fTrEffCal[i] = 0;
    fTrEffPhys[i] = 0;

  }
  for (Int_t i=0; i<24; i++)
    {
      feffC[i]=0;
      feffA[i]=0;
      feffqtc[i]=0;
   }
}

//____________________________________________________________________________ 
AliT0QADataMakerRec::AliT0QADataMakerRec(const AliT0QADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fnEventCal(0),
  fnEventPhys(0)

  
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
//____________________________________________________________________________
void AliT0QADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kT0, task, list) ;
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! IsValidEventSpecie(specie, list)) 
      continue ;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie)) ; 
    if ( task == AliQAv1::kRAWS ) {
      GetRawsData(0)->SetLabelSize(0.02);
      const Char_t *triggers[6] = {"mean", "vertex","ORA","ORC","central","semi-central"};
      for (Int_t itr=0; itr<6; itr++) {
	if ( fnEventCal>0) 
	  fTrEffCal[itr] = Float_t (fNumTriggersCal[itr])/Float_t (fnEventCal);
	if ( fnEventPhys>0) 
	  fTrEffPhys[itr] = Float_t (fNumTriggers[itr])/Float_t (fnEventPhys);

        GetRawsData(420)->Fill(triggers[itr], fTrEffCal[itr]);
        GetRawsData(420)->SetBinContent(itr+1, fTrEffCal[itr]);
        GetRawsData(97)->Fill(triggers[itr], fTrEffPhys[itr]);
        GetRawsData(97)->SetBinContent(itr+1, fTrEffPhys[itr]);
      } 
      Float_t effic=0;
      for(Int_t ik=0; ik<24; ik++)
	{  
	  effic=0;
	  if ( fnEventCal>0) effic = Float_t(feffC[ik])/Float_t(fnEventCal);
	  GetRawsData(205)->SetBinContent(ik+1,effic) ;
	  effic=0;
	  if ( fnEventCal>0) effic = Float_t(feffA[ik])/Float_t(fnEventCal);
	  GetRawsData(206)->SetBinContent(ik+1,effic );
	  effic=0;
	  if ( fnEventCal>0) effic = Float_t(feffqtc[ik])/Float_t(fnEventCal);
	  GetRawsData(207)->SetBinContent(ik+1, effic);
	}
      
      
    }
    
  }

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

  
  for (Int_t i=0; i<500; i++){
    //    low[i]  = GetRecoParam()->GetLow(i);
    ///   high[i]  = GetRecoParam()->GetHigh(i);
    low[i] = 0;
    high[i] = 10000;

  }

  TString timename, ampname, qtcname, ledname;
  TString timeCalname, ampCalname, ledCalname, qtcCalname;
  TString qt1name, qt0name, qt1Calname, qt0Calname;
  TString nhits;

  TH1F* fhRefPoint = new TH1F("hRefPoint","Ref Point", 20000, 25000 ,45000);
  Add2RawsList( fhRefPoint,0, expert, !image, !saveCorr);

  TH1F* fhRefPointcal = new TH1F("hRefPointcal","Ref Point laser", 10000, 18000 ,28000);
  Add2RawsList( fhRefPointcal,358, expert, !image, !saveCorr);

  TH1F *fhRawCFD[24]; 
  TH1F * fhRawLEDamp[24];
  TH1F *fhRawQTC[24]; TH1F * fhRawLED[24];
  TH1F *fhRawCFDcal[24]; TH1F * fhRawLEDampcal[24]; 
  TH1F *fhRawCFDcalpmt[24];  
  TH1F *fhRawCFDpmt[24];
  TH1F *fhRawQTCcal[24];  TH1F * fhRawLEDcal[24];
  TH1F *fhRawQT0cal[24]; TH1F *fhRawQT1cal[24];
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
      // fhRawCFD[i] = new THnSparseF(timename.Data(), Form("%s;CFD [#channels];Counts", timename.Data()),1,spabins,spaxmin,spaxmax);
      fhRawCFD[i] = new TH1F(timename.Data(), Form("%s;CFD [#channels];Counts", timename.Data()),Int_t(high[i+1]-low[i+1]),low[i+1],high[i+1]);
        Add2RawsList( fhRawCFD[i],i+1, expert, !image, !saveCorr);
	fhRawLED[i] = new TH1F(ledname.Data(),  Form("%s;LED[#channels];Counts", ledname.Data()),Int_t(high[i+24+1]-low[i+24+1]),low[i+24+1],high[i]);
      Add2RawsList( fhRawLED[i],i+24+1, expert, !image, !saveCorr);
      fhRawLEDamp[i] = new TH1F(ampname.Data(),  Form("%s;LED-CFD [#channels];Counts", ampname.Data()),1000,0,1000);
     Add2RawsList( fhRawLEDamp[i],i+48+1, expert, !image, !saveCorr);
      fhRawQTC[i] = new TH1F(qtcname.Data(),  Form("%s;QTC[#channels];Counts", qtcname.Data()),10000,0,10000);
      Add2RawsList( fhRawQTC[i],i+72+1, expert, !image, !saveCorr);
      fhRawQT1[i] = new TH1F(qt1name.Data(),  Form("%s;QT1[#channels];Counts", qtcname.Data()),Int_t(high[270+i]-low[270+i]),low[270+i],high[270+i]);
       Add2RawsList( fhRawQT1[i],270+i, expert, !image, !saveCorr);
       fhRawQT0[i] = new TH1F(qt0name.Data(),  Form("%s;QT0[#channels];Counts", qtcname.Data()),Int_t(high[270+24+i]-low[270+24+i]),low[270+24+i],high[270+24+i]);
     Add2RawsList( fhRawQT0[i],270+24+i, expert, !image, !saveCorr);

      fhRawNhits[i] = new TH1F(nhits.Data(),  Form("%s;#Hits;Events", nhits.Data()),10, 0, 10);
     Add2RawsList( fhRawNhits[i],244+i, expert, !image, !saveCorr);

    }
  TH1F* fhRawTrigger = new TH1F("hRawTrigger"," phys triggers;Trigger #;Counts",5,0,5);
  Add2RawsList(fhRawTrigger ,97, expert, !image, !saveCorr);

  TH1F* fhRawMean = new TH1F("hRawMean","online mean signal, physics event;",Int_t(high[98]-low[98]),low[98],high[98]);
  Add2RawsList( fhRawMean,98, expert, !image, !saveCorr);

  TH1F* fhRawVertex = new TH1F("hRawVertex","online vertex signal; counts",Int_t(high[100]-low[99]),low[99],high[99]);
  Add2RawsList( fhRawVertex,99, expert, !image, !saveCorr);

  TH1F* fhRawORA = new TH1F("hRawORA","online OR A; counts",Int_t(high[100]-low[100]),low[100],high[100]);
  Add2RawsList( fhRawORA,100, expert, !image, !saveCorr);
  TH1F* fhRawORC = new TH1F("hRawORC","online OR C;counts",Int_t( high[101]-low[101]),low[101],high[101]);
  Add2RawsList( fhRawORC,101, expert, !image, !saveCorr);


  for (Int_t i=0; i<24; i++)
    {
      // for events with trigger CALIBRATION_EVENT
      timeCalname ="hRawCFDcal";
      ledCalname = "hRawLEDcal";
      ampCalname = "hRawLEDminCFDcal";
      qtcCalname = "hRawQTCcal";
      qt0Calname = "hRawQT0cal";
      qt1Calname = "hRawQT1cal";
      timeCalname += i+1;
      ledCalname += i+1;
      ampCalname += i+1;
      qtcCalname += i+1;
      qt0Calname += i+1;
      qt1Calname += i+1;
 
      fhRawCFDcal[i] = new TH1F(timeCalname.Data(),  Form("%s;Time ;Counts", timeCalname.Data()),Int_t(high[101+i+1]-low[101+i+1]),low[101+i+1],high[101+i+1]);
      Add2RawsList( fhRawCFDcal[i],101+i+1, expert, !image, !saveCorr);
 
      fhRawLEDcal[i] = new TH1F(ledCalname.Data(),  Form("%s;Time ;Counts", ledCalname.Data()),Int_t(high[101+i+24+1]-low[101+i+24+1]),low[101+i+24+1],high[101+i+24+1]);
      Add2RawsList( fhRawLEDcal[i],101+i+24+1, expert, !image, !saveCorr);

      fhRawLEDampcal[i] = new TH1F(ampCalname.Data(), Form("%s;Amplitude [ADC counts];Counts", ampCalname.Data()),1000,0,1000);
      Add2RawsList( fhRawLEDampcal[i],101+i+48+1, expert, !image, !saveCorr);

      fhRawQTCcal[i] = new TH1F(qtcCalname.Data(), Form("%s;Charge ;Counts",qtcCalname.Data()),10000,0,10000);
      Add2RawsList( fhRawQTCcal[i],101+i+72+1, expert, !image, !saveCorr);

      fhRawQT0cal[i] = new TH1F(qt0Calname.Data(), Form("%s;Charge ;Counts",qt0Calname.Data()),Int_t(high[371+i]-low[371+i]),low[371+i],high[371+i]);
      Add2RawsList( fhRawQT0cal[i],371+i, expert, !image, !saveCorr);

      fhRawQT1cal[i] = new TH1F(qt1Calname.Data(), Form("%s;Charge ;Counts",qt1Calname.Data()),Int_t(high[i+371+24]-low[i+371+24]),low[i+371+24],high[i+371+24]);
      Add2RawsList( fhRawQT1cal[i],i+371+24, expert, !image, !saveCorr);

      }

  //from PMT1 (equalizing)
  for (Int_t i=0; i<24; i++)
    {
      // for events with trigger CALIBRATION_EVENT
      timeCalname ="hRawCFDcalpmt";
      timename ="hRawCFDpmt";
      timeCalname += i+1;
      ledCalname += i+1;
      timename += i+1;
      ledname += i+1;
      fhRawCFDcalpmt[i] = new TH1F(timeCalname.Data(),  Form("%s;Time;Counts", timeCalname.Data()),2000,-1000,1000);
       Add2RawsList( fhRawCFDcalpmt[i],321+i , expert, !image, !saveCorr);

      fhRawCFDpmt[i] = new TH1F(timename.Data(),  Form("%s;Time;Counts", timename.Data()),2000,-1000,1000);
      Add2RawsList( fhRawCFDpmt[i],220+i, expert, !image, !saveCorr);
     }

  TH1F* fhRawTriggerCal = new TH1F("hRawTriggerCal"," laser triggers",6,0,6);
  Add2RawsList(fhRawTriggerCal ,420 , !expert, image, !saveCorr);
 
  TH1F* fhRawMeanCal = new TH1F("hRawMeanCal","online mean signal, calibration event",Int_t(high[198]-low[198]),low[198],high[198]);
  Add2RawsList( fhRawMeanCal,198, expert, !image, !saveCorr);
 
  TH1F* fhRawVertexCal = new TH1F("hRawVertexCal","online vertex signal, calibration event ",Int_t( high[199]-low[199]),low[199],high[199] );
   Add2RawsList( fhRawVertexCal,199, expert, !image, !saveCorr);
 
 
   TH1F* fhRawORAcal = new TH1F("hRawORAcal","laser OR A; counts",Int_t( high[200]-low[200]),low[200],high[200]);
   Add2RawsList( fhRawORAcal,200, expert, !image, !saveCorr );
 
 
   TH1F* fhRawORCcal = new TH1F("hRawORCcal","laserOR C;counts ",Int_t( high[201]-low[201]),low[201],high[201]);
   Add2RawsList( fhRawORCcal,201, expert, !image, !saveCorr);
 

  //multiplicity trigger
   TH1F* fhMultcal = new TH1F("hMultcal","full mulltiplicity;Multiplicity;Entries",Int_t( high[202]-low[202]),low[202],high[202]);
  Add2RawsList( fhMultcal,202, expert, !image, !saveCorr );
   TH1F* fhMultScal = new TH1F("hMultSemical","full multiplicity with semi-central trigger;Multiplicity;Entries",
			       Int_t( high[203]-low[203]),low[203],high[203] );
  Add2RawsList( fhMultScal,203, expert, !image, !saveCorr);
  TH1F* fhMultCcal = new TH1F("hMultCentrcal","full multiplicity with central trigger;Multiplicity;Entries", 
			      Int_t( high[204]-low[204]),low[204],high[204]);
  Add2RawsList( fhMultCcal,204, expert, !image, !saveCorr);
 
  TH1F* fhMult = new TH1F("hMult","full mulltiplicity;Multiplicity;Entries", high[216]-low[216],low[216],high[216]);
  Add2RawsList( fhMult,216, expert, !image, !saveCorr );
   TH1F* fhMultS = new TH1F("hMultSemi","full multiplicity with semi-central trigger;Multiplicity;Entries",
			    Int_t( high[217]-low[217]),low[217],high[217] );
  Add2RawsList( fhMultS,217, expert, !image, !saveCorr);
  TH1F* fhMultC = new TH1F("hMultCentr","full multiplicity with central trigger;Multiplicity;Entries", 
			     high[218]-low[218],low[218],high[218]);
  Add2RawsList( fhMultC,218, expert, !image, !saveCorr);

  TH1F* fhEffCFD = new TH1F("hEffCFD","#PMT; #CFD counts/nEvents",24, 0 ,24); 
  Add2RawsList( fhEffCFD,205, !expert, image, !saveCorr);
 
  TH1F* fhEffLED = new TH1F("hEffLED","#PMT; #LED counts/nEvent",24, 0 ,24);
  Add2RawsList( fhEffLED,206, !expert, image, !saveCorr);

  TH1F* fhEffQTC = new TH1F("hEffQTC","#PMT; QTC efficiency%s;",24, 0 ,24);
  Add2RawsList( fhEffQTC,207, !expert, image, !saveCorr);

  
  TH2F* fhCFDcal = new TH2F("hCFDcal","CFD laser; #PMT; CFD {#channnels}",25, 0 ,25,Int_t(high[208]-low[208]),low[208],high[208]);
  fhCFDcal->SetOption("COLZ");
  Add2RawsList( fhCFDcal,208, !expert, image, !saveCorr);


  TH2F* fhLEDcal = new TH2F("hLEDcal","LED laser; #PMT; LED [#channnels]",25, 0 ,25,Int_t( high[209]-low[209]),low[209],high[209]);
  fhLEDcal->SetOption("COLZ");
  Add2RawsList( fhLEDcal,209, !expert, image, !saveCorr);

  TH2F* fhQTCcal = new TH2F("hQTCcal","QTC laser; #PMT; QTC [#channnels]",25, 0 ,25,Int_t( high[210]-low[210]),low[210],high[210]);
  fhQTCcal->SetOption("COLZ");
  Add2RawsList( fhQTCcal,210, !expert, image, !saveCorr);


  TH1F* fhNumPMTA= new TH1F("hNumPMTA","number of PMT hitted per event",13, 0 ,13);
  Add2RawsList(fhNumPMTA ,211, expert, !image, !saveCorr);

  TH1F* fhNumPMTC= new TH1F("hNumPMTC","number of PMT hitted per event",13, 0 ,13);
  Add2RawsList(fhNumPMTC ,212, expert, !image, !saveCorr);

  TH1F* fhHitsOrA= new TH1F("fhHitsOrA","T0_OR A hit multiplicitie",10, 0 ,10);
  Add2RawsList( fhHitsOrA,213, expert, !image, !saveCorr);

  TH1F* fhHitsOrC= new TH1F("fhHitsOrC","T0_OR C hit multiplicitie",10, 0 ,10);
  Add2RawsList(fhHitsOrC ,214, expert, !image, !saveCorr);


  TH1F* fhOrCminOrA= new TH1F("fhOrCminOrA","T0_OR C - T0_OR A",10000,-5000,5000);
  Add2RawsList( fhOrCminOrA,215, expert, !image, !saveCorr);

  TH1F* fhOrCminOrAcal= new TH1F("fhOrCminOrAcal","T0_OR C - T0_OR A",10000,-5000,5000);
  Add2RawsList( fhOrCminOrAcal,219, expert, !image, !saveCorr);

  TH1F* fhOrCminOrATvdcOn= new TH1F("fhOrCminOrATvdcOn","T0_OR C - T0_OR A TVDC on",10000,-5000,5000);
  Add2RawsList( fhOrCminOrATvdcOn,350, expert, !image, !saveCorr);

  TH1F* fhOrCminOrATvdcOncal= new TH1F("fhOrCminOrATvdcOncal","T0_OR C - T0_OR A TVDC on laser",10000,-5000,5000);
  Add2RawsList( fhOrCminOrATvdcOncal,351, expert, !image, !saveCorr);

  TH1F* fhOrCminOrATvdcOff= new TH1F("fhOrCminOrATvdcOff","T0_OR C - T0_OR A TVDC off",10000,-5000,5000);
  Add2RawsList( fhOrCminOrATvdcOff,352, expert, !image, !saveCorr);

  TH1F* fhOrCminOrATvdcOffcal= new TH1F("fhOrCminOrATvdcOffcal","T0_OR C - T0_OR ATVDC off laser",10000,-5000,5000);
  Add2RawsList( fhOrCminOrATvdcOffcal,353, expert, !image, !saveCorr);

  TH1F* fhCFD1abs= new TH1F("fhCFD1abs"," CFD 1 pure ( L1) ", 20000, 0, 50000 );
  Add2RawsList(fhCFD1abs ,354, expert, !image, !saveCorr);

  TH1F* fhCFD2abs= new TH1F("fhCFD2abs"," CFD 2 pure ( L1)  ", 20000, 0, 50000 );
  Add2RawsList( fhCFD2abs,355, expert, !image, !saveCorr);

  TH1F* fhCFD1abscal= new TH1F("fhCFD1abscal"," CFD 1 pure ( L1)laser ", 20000,0 ,50000 );
  Add2RawsList(fhCFD1abscal,356, expert, !image, !saveCorr);

  TH1F* fhCFD2abscal= new TH1F("fhCFD2abscal"," CFD 2 pure ( L1)laser ", 20000, 0 ,50000 );
  Add2RawsList( fhCFD2abscal,357, expert, !image, !saveCorr);

 const Char_t *triggers[6] = {"mean", "vertex","ORA","ORC","central","semi-central"};
  for (Int_t itr=0; itr<6; itr++) {
    GetRawsData(420)->Fill(triggers[itr], fNumTriggersCal[itr]);
    GetRawsData(420)->SetBinContent(itr+1, fNumTriggersCal[itr]);
    GetRawsData(97)->Fill(triggers[itr], fNumTriggers[itr]);
    GetRawsData(97)->SetBinContent(itr+1, fNumTriggers[itr]);
  }  

  TH1F* fhT0meancal= new TH1F("fhT0meancal"," (OrA+OrC)/2 in ns Laser",  1000,-100,100) ;
  Add2RawsList(fhT0meancal,421, expert, !image, !saveCorr);
  TH1F* fhOrAnscal= new TH1F("fhOrAnscal"," OrA in ns Laser",  1000,-100,100) ;
  Add2RawsList(fhOrAnscal,422, expert, !image, !saveCorr);
  TH1F* fhOrCnscal= new TH1F("fhOrCnscal"," OrC in ns Laser",  1000,-100,100) ;
  Add2RawsList(fhOrCnscal,423, expert, !image, !saveCorr);

  TH1F* fhT0meanphys= new TH1F("fhT0meanphys"," (OrA+OrC)/2 in ns Physics", 1000,-100,100) ;
  Add2RawsList(fhT0meanphys,424, expert, !image, !saveCorr);
  TH1F* fhOrAnsphys= new TH1F("fhOrAnsphys"," OrA in ns Physics", 1000,-100,100) ;
  Add2RawsList(fhOrAnsphys,425, expert, !image, !saveCorr);
  TH1F* fhOrCnsphys= new TH1F("fhOrCnsphys"," OrC in ns Physics",  1000,-100,100) ;
  Add2RawsList(fhOrCnsphys,426, expert, !image, !saveCorr);

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

}

//____________________________________________________________________________
void AliT0QADataMakerRec::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *fhESDMean = new TH1F("hESDmean"," ESD mean; mean time[%channels]",1000,0,1000);
  Add2ESDsList(fhESDMean, 0, expert, !image) ;
  TH1F * fhESDVertex = new TH1F("hESDvertex","ESDvertex; vertex[cm];",82,-30,30);
  Add2ESDsList(fhESDVertex, 1, expert, !image) ;
  
  TH1F * fhESDResolution = new TH1F("hESDResolution","(T0A-T0C)/2 corrected by SPD vertex; ns",400,-5,5);
  Add2ESDsList(fhESDResolution, 2, !expert, image) ;

}

//____________________________________________________________________________
void AliT0QADataMakerRec::MakeRaws( AliRawReader* rawReader)
{

  Float_t latencyHPTDC = 9000;
  Float_t latencyL1    = 8.19754744000000028e+03;
  Float_t latencyL1C   = 8.19707652000000053e+03;
  Float_t latencyL1A   = 8.19801836000000003e+03;
  /*
// 2009 latency
  Float_t latencyHPTDC = 22000;
  Float_t latencyL1    = 7.76597940000000017e+03;
  Float_t latencyL1C    = 7.76528399999999965e+03;
  Float_t latencyL1A    = 7.76667480000000069e+03;
  */


  rawReader->Reset() ; 
  //fills QA histos for RAW
  Int_t shift=0;
  // Int_t refPointParam = GetRecoParam()->GetRefPoint();
  Int_t refpoint = 0;
  Int_t refPointParam = 0;

  AliT0RawReader *start = new AliT0RawReader(rawReader);
  
  if (! start->Next())
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" no raw data found!!"));
  else
    {  
      
      UInt_t type =rawReader->GetType();
      if (type == 8){ shift=101; fnEventCal++;} 
      if (type == 7){ shift=0;   fnEventPhys++;}
      Int_t allData[110][5];
      for (Int_t i0=0; i0<105; i0++)
	{
	  for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0;
	}
      for (Int_t i=0; i<105; i++) 
	for (Int_t iHit=0; iHit<5; iHit++)
	  allData[i][iHit]= start->GetData(i,iHit);
      
      if ( type == 7  && allData[0][0] > 0 ) GetRawsData(0) -> Fill( allData[0][0]);
      if ( type == 8  && allData[0][0] > 0 ) GetRawsData(358) -> Fill( allData[0][0]);
 
      refpoint = allData[refPointParam][0];
      if (refPointParam <  0 ) refpoint=0; 
      if (refPointParam == 0 ) refpoint = allData[0][0] - 5000;

      Int_t numPmtC=0;    
      for (Int_t ik = 0; ik<12; ik++)
	{
	//	for (Int_t iHt=0; iHt<1; iHt++){

	Int_t nhitsPMT=0;
	if(allData[ik+1][0]>0 && type == 7  ) numPmtC++;
	for (Int_t iHt=0; iHt<5; iHt++){
	  //cfd
	  if(allData[ik+1][iHt]>0) {
	    GetRawsData(shift+ik+1) -> Fill(allData[ik+1][iHt]-refpoint);
	    if(allData[1][0]>0) 
	    GetRawsData(shift+ik+220) -> Fill(allData[ik+1][iHt]-allData[1][0]);
	    //	    cout<<"C  cfd "<<ik<<" "<<iHt<<" "<<allData[ik+1][iHt]<<" -RF "<<refpoint<<"  "<<allData[ik+1][iHt]-refpoint<<endl;
	    if(type == 8 && iHt==0 ) {
	      feffC[ik]++;
	      GetRawsData(208)->Fill(ik+1, allData[ik+1][iHt]-refpoint);
	    } 
	    if(type == 7  )  nhitsPMT++;
	    
	  }
	  //led
	  if(allData[ik+13][iHt] > 0) { 
	    GetRawsData(shift+ik+24+1)->  Fill(allData[ik+13][iHt]-refpoint);
	    //	    cout<<"C  led "<<ik<<" "<<iHt<<" "<<allData[ik+1][iHt]<<" -RF "<<refpoint<<"  "<<allData[ik+13][iHt]-refpoint<<endl;
	    if(type == 8   && iHt==0) {
	      feffA[ik]++;
	      GetRawsData(209)->Fill(ik+1, allData[ik+13][iHt]-refpoint);
	    }
	  }
	  //led -cfd

	  if(allData[ik+13][iHt] > 0 && allData[ik+1][iHt] >0 )
	    GetRawsData(shift+ik+48+1)->
	      Fill(allData[ik+13][iHt]-allData[ik+1][iHt]);

	  //qtc
	  if(allData[2*ik+25][iHt] > 0 || allData[2*ik+26][iHt] > 0) {
	    GetRawsData(shift+ik+72+1)->
	      Fill(allData[2*ik+25][iHt]-allData[2*ik+26][iHt]);
	    GetRawsData(shift+ik+270)->Fill(allData[2*ik+26][iHt] - refpoint);
	    GetRawsData(shift+ik+24+270)->Fill(allData[2*ik+25][iHt] - refpoint);

	    if(type == 8   ) {
	      feffqtc[ik]++;
	      GetRawsData(210)->Fill(ik+1, allData[2*ik+25][iHt]-allData[2*ik+26][iHt]);
	    }
	  }
	}
	if(type == 7  ) GetRawsData(ik+244)->Fill(nhitsPMT);
      }
      if(type == 7  ) GetRawsData(212)->Fill(numPmtC);
      Int_t numPmtA=0;    
       for (Int_t ik = 12; ik<24; ik++) {
	//	for (Int_t iHt=0; iHt<1; iHt++) {
	if(allData[ik+45][0]>0 && type == 7 ) numPmtA++;
	Int_t nhitsPMT=0;
	for (Int_t iHt=0; iHt<5; iHt++) {
	  if(allData[ik+45][iHt]>0) {
	    //cfd
	    GetRawsData(shift+ik+1)->
	      Fill(allData[ik+45][iHt]-refpoint);
	    if(allData[57][0]>0) 
	    GetRawsData(shift+ik+220) -> Fill(allData[ik+45][iHt]-allData[57][0]);

	    //	    cout<<"A  cfd "<<ik<<" "<<iHt<<" "<<allData[ik+1][iHt]<<" -RF "<<refpoint<<"  "<<allData[ik+1][iHt]-refpoint<<endl;
	    if(type == 8  ) {
	      feffC[ik]++;
	      GetRawsData(208)->Fill(ik+1, allData[ik+45][iHt]-refpoint);
	    } 
	    if(type == 7  )   nhitsPMT++;
	    

	  }
	  //led
	  if(allData[ik+57][iHt] > 0 ) {
	    GetRawsData(shift+ik+24+1)->Fill(allData[ik+57][iHt]-refpoint);
	    //   cout<<"C  led "<<ik<<" "<<iHt<<" "<<allData[ik+1][iHt]<<" -RF "<<refpoint<<"  "<<allData[ik+13][iHt]-refpoint<<endl;
	    if(type == 8  ) {
	      feffA[ik]++;
	      GetRawsData(209)->Fill(ik+1, allData[ik+57][iHt]-refpoint);
	    }
	  }
	    //qtc	  
	  if(allData[2*ik+57][iHt]>0 || allData[2*ik+58][iHt]>0)
	    {
	      GetRawsData(shift+ik+72+1)-> Fill(allData[2*ik+57][iHt]-allData[2*ik+58][iHt]);
	      GetRawsData(shift+ik+270)->Fill(allData[2*ik+58][iHt]-refpoint);
	      GetRawsData(shift+ik+24+270)->Fill(allData[2*ik+57][iHt]-refpoint);
	      if(type == 8  ){
		feffqtc[ik]++;
		GetRawsData(210)->Fill(ik+1, allData[2*ik+57][iHt]-allData[2*ik+58][iHt]);
	      }
	    } 
	  //led-cfd
	  if(allData[ik+57][iHt] > 0 &&allData[ik+45][iHt]>0)
	    GetRawsData(shift+ik+48+1)->
	      Fill(allData[ik+57][iHt]-allData[ik+45][iHt]);
	}
	if(type == 7  ) GetRawsData(ik+244)->Fill(nhitsPMT);

      }
       if(type == 7  ) GetRawsData(211)->Fill(numPmtA);

      Int_t trChannel[6] = {49,50,51,52,55,56};  
 
      Float_t ch2cm = 24.4*0.029979;     
      if(type == 7)
	{
	  Int_t nhitsOrA=0;
	  Int_t nhitsOrC=0;
	  for (Int_t iHt=0; iHt<5; iHt++) {
	    // pure CFD1 & CFD2 for trigger needs
	    if(allData[1][iHt]>0) GetRawsData(354)->Fill(allData[1][iHt]);
	    if(allData[2][iHt]>0) GetRawsData(355)->Fill(allData[2][iHt]);

	    //orA-orC phys tvdc 1 
	    if((allData[51][iHt]>0 && allData[52][iHt]>0) && allData[50][iHt]>0)
	      GetRawsData(350)->Fill((allData[52][iHt]-allData[51][iHt])*ch2cm);
	    //orA-orC phys tvdc 0 
	    if((allData[51][iHt]>0 && allData[52][iHt]>0) && allData[50][iHt]<=0)
	      GetRawsData(352)->Fill((allData[52][iHt]-allData[51][iHt])*ch2cm);
	    for (Int_t itr=0; itr<6; itr++) {
	      if (allData[trChannel[itr]][iHt] >0) {
		fNumTriggers[itr]++;
		GetRawsData(98+itr)->Fill(allData[trChannel[itr]][iHt]-refpoint);
	      }
	    }
	    if(allData[51][iHt]>0 && allData[52][iHt]>0) {
	      GetRawsData(424) ->Fill(24.4*0.001*(allData[51][iHt]+allData[52][iHt])/2.-latencyHPTDC+latencyL1);
	    }
	    if(allData[51][iHt]>0 ) {
	       GetRawsData(425) ->Fill(24.4*0.001*allData[51][iHt]-latencyHPTDC+latencyL1A);
	    }
	    if(allData[52][iHt]>0 ) { 
		  GetRawsData(426) ->Fill(24.4*0.001*allData[52][iHt]-latencyHPTDC + latencyL1C);
	    }
		  

	    if(allData[51][iHt] >0) nhitsOrA++;
	    if(allData[52][iHt] >0) nhitsOrC++;
	    
	    if(allData[53][iHt]>0 && allData[54][iHt]>0) {
	      GetRawsData(216)->Fill(allData[53][iHt]-allData[54][iHt]);
	      if(allData[55][iHt]>0) GetRawsData(216)->Fill(allData[53][iHt]-allData[54][iHt]);
	      if(allData[56][iHt]>0) GetRawsData(218)->Fill(allData[53][iHt]-allData[54][iHt]);
	    }
	    if(allData[51][iHt]>0 && allData[52][iHt]>0)
	      GetRawsData(215)->Fill((allData[52][iHt]-allData[51][iHt])*ch2cm);
	    
	  }
	  GetRawsData(213)->Fill(nhitsOrA);
	  GetRawsData(214)->Fill(nhitsOrC);
	  
	}
      
      if(type == 8)
	  {
	    for (Int_t iHt=0; iHt<5; iHt++) {
	      //pure CFD 1 & 2
	    if(allData[1][iHt]>0) GetRawsData(356)->Fill(allData[1][iHt]);
	    if(allData[2][iHt]>0) GetRawsData(357)->Fill(allData[2][iHt]);

	      //orA-orC laser
	    if(allData[51][iHt]>0 && allData[52][iHt]>0)
	      GetRawsData(219)->Fill((allData[52][iHt]-allData[51][iHt])*ch2cm);
	    
	    //orA-orC laser tvdc 1 
	    if((allData[51][iHt]>0 && allData[52][iHt]>0) && allData[50][iHt]>0)
	      GetRawsData(351)->Fill((allData[52][iHt]-allData[51][iHt])*ch2cm);
	    //orA-orC laser tvdc 0 
	    if((allData[51][iHt]>0 && allData[52][iHt]>0) && allData[50][iHt]<=0)
	      GetRawsData(353)->Fill((allData[52][iHt]-allData[51][iHt])*ch2cm);
	    //trigger laser
	      for (Int_t itr=0; itr<6; itr++) {
		if(allData[trChannel[itr]][iHt]>0)
		  {
		    GetRawsData(198+itr)->
		      Fill(allData[trChannel[itr]][iHt]-refpoint);
		    fNumTriggersCal[itr]++;
		  }
	      }
	      //T0Tof in ns
	    if(allData[51][iHt]>0 && allData[52][iHt]>0)
	      GetRawsData(421) ->Fill(24.4*0.001*(allData[51][iHt]+allData[52][iHt])/2.-latencyHPTDC+latencyL1);
	      if(allData[51][iHt]>0 ) 
		GetRawsData(422) ->Fill(24.4*0.001*allData[51][iHt]-latencyHPTDC+latencyL1A);
	      if(allData[52][iHt]>0 ) 
		GetRawsData(423) ->Fill(24.4*0.001*allData[52][iHt]-latencyHPTDC+latencyL1A);
	      //mult trigger signals laser
	      if(allData[53][iHt]>0 && allData[54][iHt]>0) {
		GetRawsData(202)->Fill(allData[53][iHt]-allData[54][iHt]);
		if(allData[55][iHt]>0) GetRawsData(202)->Fill(allData[53][iHt]-allData[54][iHt]);
		if(allData[56][iHt]>0) GetRawsData(204)->Fill(allData[53][iHt]-allData[54][iHt]);
	      }
	    }
	  } 
      
           delete start;
    }
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
  for (Int_t i=0; i<24; i++)
    {
    if (digCFD->At(i)>0) {
      Int_t cfd=digCFD->At(i)- refpoint;
      GetDigitsData(0) ->Fill(i,cfd);
      GetDigitsData(1) -> Fill(i, (digLED->At(i) - digCFD->At(i)));
      GetDigitsData(2) -> Fill(i, (digQT1->At(i) - digQT0->At(i)));
    }
    }  
  
  delete digCFD;
  delete digLED;
  delete digQT0;
  delete digQT1;
  
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
      GetRecPointsData(0) -> Fill(i, frecpoints -> GetTime(i) - frecpoints -> GetTime(0)); 
    if(i>11)
      GetRecPointsData(0) -> Fill(i,  frecpoints -> GetTime(i) - frecpoints -> GetTime(12)); 
    GetRecPointsData(1) -> Fill( i, frecpoints -> GetAmp(i) - frecpoints->AmpLED(i));
  }
  Double_t mmm=frecpoints->GetOnlineMean()- frecpoints->GetMeanTime();
  GetRecPointsData(2) ->Fill(mmm);
  
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

  GetESDsData(0) -> Fill(t0time);
  GetESDsData(1)-> Fill(esd->GetT0zVertex());
  GetESDsData(2)-> Fill((orA-orC)/2.);
  
}
