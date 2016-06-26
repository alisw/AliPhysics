/// \ingroup macros
/// \file MUONBusPachEvolution.C
/// \brief Macro to help display the result of the MCHBPEVOda
///
/// The result (an AliMergeableCollection, aka histogram collection) is either
/// the root file output directly from the DA or the corresponding OCDB object
///
/// Typical usage :
///
/// .L MUONBusPatchEvolution.C+
/// AliMergeableCollection* bp = BPEVO(runNumber,"output.root"); /// using OCDB object
///
/// Alternatively, in order to test the DA (expert mode), one can use :
///
/// AliMergeableCollection* bp = BPEVO("daoutput.root","output.root");
///
/// And then :
///
/// PlotStationOccupancies(*bp);
/// PlotDDLOccupancies(*bp);
/// PlotChamberOccupancies(*bp,stationNumber);
/// PlotDEOccupancies(*bp,chamberNumber);
/// PlotBusPatchOccupancies(*bp,ddlNumber);
///
/// \author Laurent Aphecetche, Subatech

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDAQ.h"
#include "AliMergeableCollection.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpDDL.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMUONBusPatchEvolution.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TSystem.h"
#include <AliMpDDLStore.h>
#include <cassert>
#include <limits>
#include <map>
#include <set>

//_________________________________________________________________________________________________
void AssertMapping(Int_t runNumber=0) {

	AliCDBManager* cdbm = AliCDBManager::Instance();

	if (!cdbm->IsDefaultStorageSet()) {
		cdbm->SetDefaultStorage(
				"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
	}

	cdbm->SetRun(runNumber);

	AliMpCDB::LoadAll();
}

//_________________________________________________________________________________________________
AliMergeableCollection* BPEVO(Int_t runNumber, const char* output = "mchbepevo.complete.root")
{
	AliCDBManager* cdb = AliCDBManager::Instance();


	if (!cdb->IsDefaultStorageSet()) {
		cdb->SetDefaultStorage(
				"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
	}

	AssertMapping(runNumber);

	AliCDBEntry* e = cdb->Get("MUON/Calib/BPEVO");

	AliMergeableCollection* hc = static_cast<AliMergeableCollection*>(e->GetObject()->Clone());

	if (!hc) return 0x0;

	AliMUONBusPatchEvolution bpe(*hc);

	bpe.Augment();

	if (output)
	{
		TFile* fout = TFile::Open(gSystem->ExpandPathName(output), "RECREATE");
		hc->Write("bpevo");
		delete fout;
	}

	return hc;
}

//_________________________________________________________________________________________________
AliMergeableCollection* BPEVO(const char* daoutput, const char* output = "mchbepevo.complete.root")
{
	/// Get an augmented version of the mergeable collection found in the DA output file

	TFile* f = TFile::Open(gSystem->ExpandPathName(daoutput));

	AliMergeableCollection* hc = static_cast<AliMergeableCollection*>(f->Get(
			"bpevo"))->Clone();

	delete f;

	AliMUONBusPatchEvolution bpe(*hc);

	bpe.Augment();

	if (output)
	{
		TFile* fout = TFile::Open(gSystem->ExpandPathName(output), "RECREATE");
		hc->Write("bpevo");
		delete fout;
	}

	return hc;
}

//_________________________________________________________________________________________________
void MakeConfigFileForBPEVOda(const char* filename = "mchbpevo.conf") {
	/// Generate a suitable configuration file for the MCHBPEVO DA
	/// to be stored in the DAQ detector database

	std::ofstream out(gSystem->ExpandPathName(filename));

	AssertMapping();

	out << "# Histogram the hit counts in 10 and 60 seconds bins" << std::endl
			<< std::endl;
	out << "+timeResolutions: 10s" << std::endl;
	out << "+timeResolutions: 60s" << std::endl << std::endl;

	out << "# maxDuration set to 5 hours" << std::endl << std::endl;
	out << "maxDuration: 18000" << std::endl << std::endl;

	out << "# number of events needed to take a decision" << std::endl
			<< std::endl;
	out << "nofEventsRequiredForDecision: 1000" << std::endl << std::endl;

	out << "# occupancy threshold for alarm" << std::endl << std::endl;
	out << "occupancyThreshold: 0.5" << std::endl << std::endl;

	TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
	AliMpBusPatch* bp;

	std::map<int, int> buspatches;

	while ((bp = static_cast<AliMpBusPatch*>(next()))) {
		int npads(0);

		Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(bp->GetId());

		AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(
				detElemId);

		for (Int_t i = 0; i < bp->GetNofManus(); ++i) {
			Int_t manuId = bp->GetManuId(i);

			npads += de->NofChannelsInManu(manuId);
		}

		buspatches[bp->GetId()] = npads;
	}

	std::map<int, int>::const_iterator it;
	int total(0);

	out << "# +bp: buspatchId:nofPadsInThatBusPatch" << std::endl << std::endl;

	for (it = buspatches.begin(); it != buspatches.end(); ++it) {
		total += it->second;
		out << "+bp: " << it->first << ":" << it->second << std::endl;
	}

	assert(total == 1064008);
}

//_________________________________________________________________________________________________
void MakeConfigCodeForBPEVOda(const char* buspatchmapname = "buspatches") {
	/// Generate a piece of code to make the default buspatch map
	/// used in the MCHBPEVO DA

	AssertMapping();

	TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
	AliMpBusPatch* bp;

	int total(0);

	while ((bp = static_cast<AliMpBusPatch*>(next()))) {
		int npads(0);

		Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(bp->GetId());

		AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(
				detElemId);

		for (Int_t i = 0; i < bp->GetNofManus(); ++i) {
			Int_t manuId = bp->GetManuId(i);

			npads += de->NofChannelsInManu(manuId);
		}

		std::cout << Form("%s[%d]=%d;", buspatchmapname, bp->GetId(), npads)
						<< std::endl;

		total += npads;
	}
	assert(total == 1064008);
}


//_________________________________________________________________________________________________
void plot(const std::vector<TH1*>& v, Double_t min, Double_t max, Double_t factor)
{
	if ( !v.size() ) return;

	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	if (!v[0]) return;

	TString title = v[0]->GetName();
	title.ReplaceAll("RIGHT","");
	title.ReplaceAll("LEFT","");

	TLegend* legend = new TLegend(0.8,0.1,0.99,0.90,title.Data());
	legend->SetTextSize(0.08);
	gPad->SetRightMargin(0.22);

	for ( std::vector<TH1*>::size_type i = 0; i < v.size(); ++i )
	{
		if (!v[i]) continue;

		TString tf = v[i]->GetXaxis()->GetTimeFormat();

		tf.ReplaceAll(":%s","");
		v[i]->GetXaxis()->SetTimeFormat(tf.Data());

		v[i]->Scale(factor);
		v[i]->SetMaximum(max);
		v[i]->SetMinimum(min);
		v[i]->GetYaxis()->SetLabelSize(0.10);
		v[i]->GetXaxis()->SetLabelSize(0.10);
		v[i]->GetYaxis()->SetRange(min,max);
		Double_t integralError;
		Double_t integral = v[i]->IntegralAndError(1,v[i]->GetNbinsX(),integralError);
		integralError /= integral;
		integral /= v[i]->GetNbinsX();
		integralError *= integral;
		TLegendEntry* entry = legend->AddEntry(v[i],Form("Mean %5.2g #pm %5.2g %%",integral,integralError));
		entry->SetTextSize(0.04);
		if (i)
		{
			v[i]->Draw("same");
		}
		else
		{
			v[i]->Draw();
		}
	}
	legend->Draw();

}

//_________________________________________________________________________________________________
void PlotStationOccupancies(AliMergeableCollection& hc,
		int timeResolution = 60) {

	TCanvas* c = new TCanvas("station", "station");
	c->Divide(1, 5);

	for (int i = 1; i <= 5; ++i) {

		c->cd(i);

		std::vector<TH1*> v;

		v.push_back(hc.Histo(Form("/STATION/OCC/%ds/STATION%dRIGHT", timeResolution, i)));
		v.push_back(hc.Histo(Form("/STATION/OCC/%ds/STATION%dLEFT", timeResolution, i)));
		v.push_back(hc.Histo(Form("/STATION/OCC/%ds/STATION%d", timeResolution, i)));

		gPad->SetLogy();


		plot(v,1E-2,10,100.0);

	}
}

//_________________________________________________________________________________________________
void PlotChamberOccupancies(AliMergeableCollection& hc,
														int stationId = 3,
														int timeResolution = 60) {

	TCanvas* c = new TCanvas("chamber", "chamber");

	c->Divide(1, 2);

	for (int i = 1; i <= 2; ++i) {

		c->cd(i);

		std::vector<TH1*> v;

		int chamberId = (stationId-1)*2 + i;

		v.push_back(hc.Histo(Form("/CHAMBER/OCC/%ds/CHAMBER%dRIGHT", timeResolution, chamberId)));
		v.push_back(hc.Histo(Form("/CHAMBER/OCC/%ds/CHAMBER%dLEFT", timeResolution, chamberId)));
		v.push_back(hc.Histo(Form("/CHAMBER/OCC/%ds/CHAMBER%d", timeResolution, chamberId)));

		gPad->SetLogy();

		plot(v,1E-2,10,100.0);
	}
}

//_________________________________________________________________________________________________
void PlotDDLOccupancies(AliMergeableCollection& hc,
														int timeResolution = 60) {

	TCanvas* c = new TCanvas("ddl", "ddl");

	c->Divide(4, 5);

	for ( int i = 1; i < 20; ++i ) {

		c->cd(i);

		std::vector<TH1*> v;

		v.push_back(hc.Histo(Form("/DDL/OCC/%ds/DDL%04d", timeResolution, 2560+i)));

		gPad->SetLogy();

		plot(v,1E-2,10,100.0);
	}
}

//_________________________________________________________________________________________________
void PlotDEOccupancies(AliMergeableCollection& hc,
														int chamberId = 5,
														int timeResolution = 60) {

	AssertMapping();

	TCanvas* c = new TCanvas("de", "de");

	c->Divide(5, 5);

	AliMpDEIterator it;

	it.First(chamberId-1);

	int i(1);

	while (!it.IsDone()) {

		c->cd(i);

		int detElemId = it.CurrentDEId();

		std::vector<TH1*> v;

		v.push_back(hc.Histo(Form("/DE/OCC/%ds/DE%04d", timeResolution, detElemId)));

		it.Next();

		gPad->SetLogy();

		plot(v,1E-2,10,100.0);

		++i;
	}
}

//_________________________________________________________________________________________________
void PlotBusPatchOccupancies(AliMergeableCollection& hc,
														 int ddlId = 2575,
														  int timeResolution = 60) {

	AssertMapping();

	TCanvas* c = new TCanvas("bp", "bp");

	AliMpDDL* ddl = AliMpDDLStore::Instance()->GetDDL(ddlId - AliDAQ::DdlIDOffset("MUONTRK"));

	int nbuspatches = ddl->GetNofBusPatches();

	c->DivideSquare(nbuspatches);

	for ( int i = 1; i < nbuspatches; ++i ) {

		c->cd(i);

		std::vector<TH1*> v;

		v.push_back(hc.Histo(Form("/BUSPATCH/OCC/%ds/BP%04d", timeResolution, ddl->GetBusPatchId(i))));

		gPad->SetLogy();

		plot(v,1E-1,100,100.0);
	}
}

