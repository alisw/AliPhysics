//===========================================================
// AliJEbePercentile.h
//  
//   J
//===========================================================

#ifndef ALIJEBEPERCENTILE_H
#define ALIJEBEPERCENTILE_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>
#include  "AliJConst.h"
#include  "AliJCard.h"

class AliJCard;

const int NBVn = 100;

using namespace std;

class AliJEbePercentile {

	public:
		AliJEbePercentile();
		virtual ~AliJEbePercentile(){delete fFile;}    //destructor
		AliJEbePercentile(const AliJEbePercentile& obj);
		AliJEbePercentile(AliJCard *fCard, TString file);
		AliJEbePercentile& operator=(const AliJEbePercentile& obj);

		double GetEbeFlowPercentile(int cBin, int ih, double vn);


	protected:
		TFile *fFile;
		AliJCard *fCard;
		double vnLimit[kMaxNoCentrBin][kNHarmonics][NBVn];
		double vnPercentile[NBVn+1];
		TH1D *vnobs[kMaxNoCentrBin][kNHarmonics];
};

#endif






















