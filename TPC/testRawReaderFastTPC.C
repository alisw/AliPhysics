/// \file testRawReaderFastTPC.C

#include <stdio.h>
#include <TString.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>

#include "AliLog.h"

#include "AliRawReader.h"
#include "AliRawReaderRoot.h"

#include "AliTPCRawStream.h"
#include "AliTPCRawStreamFast.h"

#include "AliAltroRawStream.h"
#include "AliAltroRawStreamFast.h"

#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"

/*
    .L testRawReaderFastTPC.C+
    testRawReaderFastTPC()
*/


AliTPCCalPad *testRawReaderFastTPC(const Char_t *file="/data.local/data/06000002142000.1A.root")
{

    AliLog::SetGlobalDebugLevel(0) ;
    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    //TString filename("/d/alice05/testtpc/raw/pulser/06000002142000.1A.root");  //nfs
    //TString filename("root://lxfs35.gsi.de:1094//alice/testtpc/raw2006/06000001537001.001.root");

    TString filename(file);  //local
    // on castor: /castor/cern.ch/alice/data/2006/09/18/15/06000002142000.1A.root


    printf("File: %s\n", filename.Data());

    AliRawReader *rawReader = new AliRawReaderRoot(filename);
    if ( !rawReader ) return 0x0;
    rawReader->RewindEvents();

    AliTPCRawStreamFast *sf = new AliTPCRawStreamFast(rawReader);
    AliTPCRawStream      *s = new AliTPCRawStream(rawReader);

    s->SetNoAltroMapping(kFALSE);

    Int_t ievent = 0;
    Int_t ievent2 = 0;
    Int_t count=0;

    AliTPCCalPad *padOld = new AliTPCCalPad("old","old");
    AliTPCCalPad *padNew = new AliTPCCalPad("new","new");
    AliTPCCalPad *padSum = new AliTPCCalPad("sum","sum");
    AliTPCCalPad *padhadd = new AliTPCCalPad("sum","sum");

    for (Int_t sec=0; sec<72; sec++){
	AliTPCCalROC *roc = padhadd->GetCalROC(sec);
	for (UInt_t ch=0; ch<roc->GetNchannels(); ch++)
	    roc->SetValue(ch,-1);
    }

    TStopwatch timer1;
    TStopwatch timer2;
    TStopwatch timer3;

    TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
    if ( !c1 ) c1 = new TCanvas("c1","c1");
    c1->Clear();
    c1->Divide(2,2);


    while (rawReader->NextEvent()){
	//old algorithm
	printf("\nevent: %d (%d)\n",ievent, ievent2);

	Bool_t input=kFALSE;

//	if ( ievent != 19 ) input=kTRUE;
//	else{
	    timer1.Start();timer2.Start(kFALSE);

	    Int_t sum1 = 0;
            Int_t sum2 = 0;
	while ( s->Next() ){
            AliTPCCalROC *roc = padhadd->GetCalROC(s->GetSector());

            //check if hwaddress gets overwritten
	    if ( roc->GetValue(s->GetRow(), s->GetPad()) == -1 )
	        roc->SetValue(s->GetRow(), s->GetPad(), s->GetHWAddress());
	    else
		if ( roc->GetValue(s->GetRow(), s->GetPad()) != s->GetHWAddress()){
		    printf("#%.2d: %.2d.%.3d.%.2d,%.4d - old [%.4d] now[%.4d]",
			   ievent, s->GetSector(), s->GetRow(),s->GetPad(), s->GetTime(),
			   (Int_t)roc->GetValue(s->GetRow(), s->GetPad()), s->GetHWAddress());
		}


	    Float_t val=padOld->GetCalROC(s->GetSector())->GetValue(s->GetRow(), s->GetPad());
	    padOld->GetCalROC(s->GetSector())->SetValue(s->GetRow(), s->GetPad(), s->GetSignal()+val);

/*	    if ( ievent == 19 && s->GetSector() == 25 && s->GetRow() == 00 && s->GetPad()==9 ){
		printf("old         | (%.2d.%.3d.%.2d.%.4d | %.3d ): %.3d\n",
		       s->GetSector(), s->GetRow(), s->GetPad(), s->GetTime(), s->GetHWAddress(), s->GetSignal());
                sum1+=s->GetSignal();
	    }
*/

	    input=kTRUE;
	    count++;
	}
	timer1.Stop();timer2.Stop();
	printf("old  --  Time: %.4f (%.4f)\n", timer1.RealTime(), timer1.CpuTime());

	rawReader->Reset();

	//new algorithm
        timer1.Start();timer3.Start(kFALSE);
	while ( sf->NextDDL() ){
	    while ( sf->NextChannel() ){
		UInt_t signal=0;
		while ( sf->NextBunch() ){
		    for (UInt_t timebin=sf->GetStartTimeBin(); timebin<sf->GetEndTimeBin(); timebin++){

			AliTPCCalROC *roc = padhadd->GetCalROC(sf->GetSector());
//			if ( roc->GetValue(sf->GetRow(), sf->GetPad()) == -1 )
//			    roc->SetValue(sf->GetRow(), sf->GetPad(), sf->GetHWAddress());
//			else
			    if ( roc->GetValue(sf->GetRow(), sf->GetPad()) != sf->GetHWAddress()){
				printf("#%.2d: %.2d.%.3d.%.2d,%.4d - old [%.4d] now[%.4d]",
				       ievent, sf->GetSector(), sf->GetRow(),sf->GetPad(), timebin+1,
				       (Int_t)roc->GetValue(sf->GetRow(), sf->GetPad()), sf->GetHWAddress());
			    }
                            /*
			    Int_t sig = sf->GetSignals()[timebin-sf->GetStartTimeBin()];
			    if ( ievent == 19 && sf->GetSector() == 25 && sf->GetRow() == 00 && sf->GetPad()==9  ){
				printf("new         | (%.2d.%.3d.%.2d.%.4d | %.3d ): %.3d\n",
				       sf->GetSector(), sf->GetRow(), sf->GetPad(), timebin+1, sf->GetHWAddress(), sig );
				sum2+=sig;
			    }
                            */
			signal+=sf->GetSignals()[timebin-sf->GetStartTimeBin()];
		    }
		    padNew->GetCalROC(sf->GetSector())->SetValue(sf->GetRow(),sf->GetPad(),signal);
		}
	    }
	}
	timer1.Stop();timer3.Stop();
	printf("new  --  Time: %.4f (%.4f)\n", timer1.RealTime(), timer1.CpuTime());

	printf("sum1: %d, sum2: %d, diff: %d\n",sum1,sum2,sum1-sum2);

	AliTPCCalPad com(*padOld);
	com.Add(padNew, -1);

	c1->cd(1);
	padOld->MakeHisto2D(1)->Draw("colz");
	c1->cd(2);
	padNew->MakeHisto2D(1)->Draw("colz");
	c1->cd(3);
	com.MakeHisto2D(1)->Draw("colz");

        //loop over all sectors, rows, pads
	for ( UInt_t iSec=0; iSec<72; iSec++ )
	    for ( UInt_t iRow=0; iRow<com.GetCalROC(iSec)->GetNrows(); iRow++ )
		for ( UInt_t iPad=0; iPad<com.GetCalROC(iSec)->GetNPads(iRow); iPad++){
		    Float_t val = com.GetCalROC(iSec)->GetValue(iRow, iPad);
		    Float_t valo = padNew->GetCalROC(iSec)->GetValue(iRow, iPad);
                    Float_t valn = padOld->GetCalROC(iSec)->GetValue(iRow, iPad);
                    //check if values for old and new algorithm differ
//		    if ( val != 0 ){
		    if ( (Int_t)valo != (Int_t)valn ){
                        Float_t hadd = padhadd->GetCalROC(iSec)->GetValue(iRow,iPad);
			printf("Event: %.2d | (%.2d.%.3d.%.2d | %f ): %.3f\n", ievent, iSec, iRow, iPad, hadd, val);
			padSum->GetCalROC(iSec)->SetValue(iRow, iPad, padSum->GetCalROC(iSec)->GetValue(iRow, iPad)+1);
		    }
		    com.GetCalROC(iSec)->SetValue(iRow, iPad, 0);
		    padOld->GetCalROC(iSec)->SetValue(iRow, iPad, 0);
		    padNew->GetCalROC(iSec)->SetValue(iRow, iPad, 0);
		}


	c1->cd(4);
	padSum->MakeHisto2D(1)->Draw("colz");

	c1->Modified();
	c1->Update();
  //      }// end event sel
	if (input) ievent++;
	ievent2++;
    }
    TFile f("output.root","recreate");
    padSum->Write();
    f.Save();
    f.Close();

    printf("total old  --  Time: %.4f (%.4f)\n", timer2.RealTime(), timer2.CpuTime());
    printf("total new  --  Time: %.4f (%.4f)\n", timer3.RealTime(), timer3.CpuTime());

    delete rawReader;
    return padSum;
}

