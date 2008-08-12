/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#define USE_ROOT
#include "AliHLTHOMERReader.h"
#include "AliHLTMessage.h"
#include "Riostream.h"
#include "Rtypes.h"
#include "TH1D.h"
#include "TObject.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TCanvas.h"
#include <errno.h>
#else
#error You must compile this macro. Try the following command "> aliroot MonitorRawData.C++"
#endif


TH1D* gErrorHist = NULL;
TH1D* gManuHist = NULL;
TH1D* gSignalHist = NULL;


const char* DataTypeToString(homer_uint64 type)
{
	union
	{
		homer_uint64 val;
		char bytes[8];
	};
	val = type;
	
	static char str[9];
	for (int i = 0; i < 8; i++)
	{
		str[i] = bytes[7-i];
	}
	str[8] = '\0'; // Null terminate the string.
	return str;
}


const char* OriginToString(homer_uint32 origin)
{
	union
	{
		homer_uint32 val;
		char bytes[4];
	};
	val = origin;
	
	static char str[5];
	for (int i = 0; i < 4; i++)
	{
		str[i] = bytes[3-i];
	}
	str[4] = '\0'; // Null terminate the string.
	return str;
}


bool UpdateHists(AliHLTHOMERReader& homerReader, const char* hostname, UShort_t port, bool addHists)
{
	/// This routine just reads all the data blocks from the HOMER reader interface,
	/// and checks if they are histogram objects we can deal with.
	/// Any histogram objects with the following strings in their names:
	/// "rawDataErrors", "manuDistrib" or "signalDistrib"
	/// will be added to the global cumulative histograms.

	int result = homerReader.ReadNextEvent(3000000);  // 3 second timeout.
	TTimeStamp now;
	if (result != 0)
	{
		if (result != ETIMEDOUT)
		{
			cerr << "ERROR: Could not read another event from HLT on "
				<< hostname << ":" << port << "\t" << now.AsString() << endl;
			return false;
		}
		else
		{
			cerr << "Timed out when trying to read from HLT on "
				<< hostname << ":" << port << "\t" << now.AsString() << endl;
			return true;
		}
	}	
	
	cout << "Reading stats up to event: " << homerReader.GetEventID()
		<< "\t" << now.AsString() << endl;
	
	for (unsigned long n = 0; n < homerReader.GetBlockCnt(); n++)
	{
		char* buffer = new char[homerReader.GetBlockDataLength(n)];
		memcpy(buffer, homerReader.GetBlockData(n), homerReader.GetBlockDataLength(n));
		AliHLTMessage msg(buffer, homerReader.GetBlockDataLength(n));
		TClass* objclass = msg.GetClass();
		if (objclass == NULL)
		{
			cerr << "WARNING: Do not know how to handle block " << n
				<< ", block type = " << DataTypeToString(homerReader.GetBlockDataType(n))
				<< ", origin = " << OriginToString(homerReader.GetBlockDataOrigin(n))
				<< ", but unknown class type. Skipping block." << endl;
			//delete buffer;
			continue;
		}
		TObject* obj = msg.ReadObject(objclass);
		if (obj == NULL)
		{
			cerr << "WARNING: Could not read object from block " << n
				<< ", class type = " << objclass->GetName()
				<< ". Skipping block." << endl;
			delete buffer;
			continue;
		}
		if (obj->IsA() != TH1D::Class())
		{
			cerr << "WARNING: Do not know how to handle object of type "
				<< obj->ClassName() << " recevied in block " << n
				<< ". Skipping block." << endl;
			delete buffer;
			continue;
		}
		
		TH1D* hist = static_cast<TH1D*>(obj);
		TString name = hist->GetName();
		if (name.Contains("rawDataErrors"))
		{
			if (! addHists) gErrorHist->Reset("M");
			gErrorHist->Add(hist);
		}
		else if (name.Contains("manuDistrib"))
		{
			if (! addHists) gManuHist->Reset("M");
			gManuHist->Add(hist);
		}
		else if (name.Contains("signalDistrib"))
		{
			if (! addHists) gSignalHist->Reset("M");
			gSignalHist->Add(hist);
		}
		else
		{
			cerr << "WARNING: Do not know how to handle histogram " << name
				<< " found in data block " << n << endl;
		}
		delete buffer;
	}
	
	return true;
}


void MonitorRawData(const char* hostname = "alihlt-dcs0", UShort_t port = 58784, bool addHists = false)
{
	if (gErrorHist == NULL && gDirectory->Get("errorHist") == NULL)
	{
		gErrorHist = new TH1D("errorHist", "Error codes found in raw data", 40, 0.5, 40.5);
		gErrorHist->SetXTitle("Error code");
		gErrorHist->SetYTitle("Number of errors");
	}
	if (gManuHist == NULL && gDirectory->Get("manuHist") == NULL)
	{
		gManuHist = new TH1D("manuHist", "Distribution of signals found per MANU", 2048, -0.5, 2047.5);
		gManuHist->SetXTitle("MANU number (as seen in raw data)");
		gManuHist->SetYTitle("Number of signals received.");
	}
	if (gSignalHist == NULL && gDirectory->Get("signalHist") == NULL)
	{
		gSignalHist = new TH1D("signalHist", "Distribution of ADC signal values", 4096, -0.5, 4095.5);
		gSignalHist->SetXTitle("Channels");
		gSignalHist->SetYTitle("dN/dChannel");
	}
	
	TCanvas* c1 = new TCanvas("errorCanvas", "Error histogram", 0, 0, 600, 450);
	gErrorHist->Draw();
	TCanvas* c2 = new TCanvas("manuCanvas", "MANU histogram", 610, 00, 600, 450);
	if (gManuHist->GetEntries() > 0) c2->SetLogy();
	gManuHist->Draw();
	TCanvas* c3 = new TCanvas("signalCanvas", "Signal histogram", 0, 480, 600, 450);
	if (gSignalHist->GetEntries() > 0) c3->SetLogy();
	gSignalHist->Draw();
	
	for (;;)
	{
		bool updateOk = true;

		AliHLTHOMERReader homerReader(hostname, port);
		int status = homerReader.GetConnectionStatus();
		if (status == 0)
		{
			try
			{
				updateOk = UpdateHists(homerReader, hostname, port, addHists);
			}
			catch (...)
			{
				cerr << "ERROR: exception occured. Trying to recover and continue..." << endl;
			}
		}
		else
		{
			cerr << "ERROR: Could not connect to HLT on " << hostname << ":" << port << endl;
		}

		// Clear histograms on errors. This is usualy due to end of run.
		if (! updateOk || status != 0)
		{
			gErrorHist->Reset("M");
			gManuHist->Reset("M");
			gSignalHist->Reset("M");
			c1->cd();
			gErrorHist->Draw();
			c2->cd();
			c2->SetLogy(kFALSE);
			gManuHist->Draw();
			c3->cd();
			c3->SetLogy(kFALSE);
			gSignalHist->Draw();
		}

		if (gManuHist->GetEntries() > 0) c2->SetLogy();
		if (gSignalHist->GetEntries() > 0) c3->SetLogy();

		c1->Update();
		c2->Update();
		c3->Update();

		// Effectively wait for 2 seconds but dispatch ROOT events while waiting.
		for (int i = 0; i < 200; i++)
		{
			gSystem->DispatchOneEvent(kTRUE);
			gSystem->Sleep(10);
		}
	}
}

