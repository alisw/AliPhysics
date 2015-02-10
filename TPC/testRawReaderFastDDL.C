/// \file testRawReaderFastDDL.C
///
/// compare old and new AltroRawStream algorithms for each pad and timebin
///
/// check if bins are filled twice (which should not be the case!!!)

#include <stdio.h>
#include <TString.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TH2I.h>
#include <TFile.h>

#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliLog.h"

#include "AliAltroRawStreamFast.h"
#include "AliAltroRawStream.h"


void testRawReaderFastDDL(const Char_t *file="/data.local/data/06000002142000.1A.root")
{
    // set minimal screen output
    AliLog::SetGlobalDebugLevel(0) ;
    AliLog::SetGlobalLogLevel(AliLog::kFatal);

//    TString filename("/d/alice05/testtpc/raw/pulser/06000002142000.1A.root");
    TString filename(file);

    printf("File: %s\n", filename.Data());

    AliRawReader *rawReader = new AliRawReaderRoot(filename);
    if ( !rawReader ) return;
    rawReader->RewindEvents();

    AliAltroRawStreamFast *sf = new AliAltroRawStreamFast(rawReader);
    AliAltroRawStream      *s = new AliAltroRawStream(rawReader);

    s->SetNoAltroMapping(kFALSE);
    s->SelectRawData("TPC");

    Int_t ievent = 0;
    Int_t count=0;

    TH2I *h2ddlT1[216];
    TH2I *h2ddlT2[216];

    for ( Int_t i=0; i<216; i++ ){
	h2ddlT1[i] = 0x0;
	h2ddlT2[i] = 0x0;
    }

    TStopwatch timer1;
    TStopwatch timer2;
    TStopwatch timer3;

    while (rawReader->NextEvent()){
	printf("\nevent: %d\n",ievent);
	Bool_t input=kFALSE;


	//old algorithm
	timer1.Start();timer2.Start(kFALSE);
	while ( s->Next() ){
            if ( !h2ddlT1[s->GetDDLNumber()] ) h2ddlT1[s->GetDDLNumber()] = new TH2I(Form("hddl1_%d",s->GetDDLNumber()),"h2c1",3584,0,3584,1024,0,1024);
	    TH2I *hist = h2ddlT1[s->GetDDLNumber()];

	    //fast filling, TH1::Fill takes awfully long
	    Int_t bin=(s->GetTime()+1)*(3584+2)+s->GetHWAddress()+1;
	    // check if this bin was allready filled
	    if ( hist->GetArray()[bin] > 0 )
		printf(" not 0:   |  %.3d : %.3d (%.3d)\n",
		        s->GetHWAddress(), s->GetSignal(), hist->GetArray()[bin]);
	    else
		hist->GetArray()[bin]=s->GetSignal();


	    input=kTRUE;
	    count++;
	}
	timer1.Stop();timer2.Stop();
	printf("old  --  Time: %.4f (%.4f)\n", timer1.RealTime(), timer1.CpuTime());
	// END OLD

	rawReader->Reset();

	//new algorithm
        timer1.Start();timer3.Start(kFALSE);
	while ( sf->NextDDL() ){
	    if ( !h2ddlT2[sf->GetDDLNumber()] ) h2ddlT2[sf->GetDDLNumber()] = new TH2I(Form("hddl2_%d",s->GetDDLNumber()),"h2c1",3584,0,3584,1024,0,1024);
	    TH2I *hist = h2ddlT2[sf->GetDDLNumber()];
	    while ( sf->NextChannel() ){
		UInt_t signal=0;
		while ( sf->NextBunch() ){
		    for (UInt_t timebin=sf->GetStartTimeBin(); timebin<sf->GetEndTimeBin(); timebin++){
			signal=sf->GetSignals()[timebin-sf->GetStartTimeBin()];

			//fast filling, TH1::Fill takes awfully long
			Int_t bin=(timebin+1+1)*(3584+2)+sf->GetHWAddress()+1; // timebins of old and new algorithm differ by 1!!!
                        // check if this bin was allready filled
			if ( hist->GetArray()[bin] > 0 )
			    printf(" not 0:   |  %.3d : %.3d (%.3d)\n",
				    sf->GetHWAddress(), signal, hist->GetArray()[bin]);
			else
			    hist->GetArray()[bin]=signal;
		    }
		}
	    }
	}
	timer1.Stop();timer3.Stop();
	printf("new  --  Time: %.4f (%.4f)\n", timer1.RealTime(), timer1.CpuTime());
	// END NEW

        //check if all data are the same for both algorithms
	for ( Int_t ddl=0; ddl<216; ddl++ ){
	    TH2I *hist = h2ddlT1[ddl];
	    if ( !hist ) continue;
	    TH2I *hist2 = h2ddlT2[ddl];
	    for ( Int_t hadd=0; hadd<3584; hadd++ )
		for ( Int_t time=0; time<1024; time++ ){
		    Int_t bin=(time+1)*(3584+2)+hadd+1;
		    Int_t val1 = hist->GetArray()[bin];
		    Int_t val2 = hist2->GetArray()[bin];
		    if ( val1 != val2 )
			printf("%.2d. %.3d %.4d %.4d: %d - %d = %d\n", ievent, ddl, hadd, time, val1, val2, val1-val2);

                    //reset for the next event
		    hist->GetArray()[bin]=0;
		    hist2->GetArray()[bin]=0;
		}
	}

	if (input) ievent++;
    }

    printf("total old  --  Time: %.4f (%.4f)\n", timer2.RealTime(), timer2.CpuTime());
    printf("total new  --  Time: %.4f (%.4f)\n", timer3.RealTime(), timer3.CpuTime());

    delete rawReader;
    delete [] h2ddlT1;
    delete [] h2ddlT2;
}

