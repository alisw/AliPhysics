#include <TGrid.h>
#include  "AliJEbePercentile.h"

//______________________________________________________________________________
AliJEbePercentile::AliJEbePercentile():
	fFile(NULL),
	fCard(NULL)
{
}
//______________________________________________________________________________
AliJEbePercentile::AliJEbePercentile(const AliJEbePercentile& obj):
	fFile(obj.fFile),
	fCard(obj.fCard)
{
	// copy constructor
}

AliJEbePercentile::AliJEbePercentile(AliJCard *fcard, TString input):
	fFile(NULL),
	fCard(NULL)
{

	fCard = fcard;
	cout << input << endl;
	
	if (TString(input).BeginsWith("alien:"))  TGrid::Connect("alien:");
	fFile = TFile::Open(input);

	fFile->Print();

	for(int ibv = 0; ibv < NBVn ; ibv++){
		vnPercentile[ibv] = ibv*1.0/NBVn;
	}

	TString name;
	int NH = 4;
	for(int ic = 0; ic<fCard->GetNoOfBins(kCentrType);ic++){
		for(int ih=2;ih<NH;ih++) {
			name = "hVnObsVector";
			vnobs[ic][ih] = (TH1D*)fFile->Get(name.Append(Form("%02d%02d",ic,ih)));
		}
	}

	double entr;
	double sum;
	int found[NBVn];

	for(int ic = 0 ; ic < fCard->GetNoOfBins(kCentrType); ic++){
		for(int ih = 2; ih < NH ; ih++){
			vnLimit[ic][ih][0]=1;
			vnLimit[ic][ih][NBVn-1]=0;
			for(int ibv = 1 ; ibv < NBVn ; ibv++){
				found[ibv] = 0;
			}
			sum = 0;
			entr = vnobs[ic][ih]->GetEntries();
			//cout << "Centrality " << ic << " of " << NC ;
			//cout << " Harmonic " << ih << " of " << NH << endl;
			for(int ib = 1; ib < vnobs[ic][ih]->GetNbinsX() +1 ; ib++){
				sum = sum + vnobs[ic][ih]->GetBinContent(ib);
				for(int ibv = 1; ibv < NBVn ; ibv++){
					if( sum/entr > 1- vnPercentile[ibv] && found[ibv] == 0){
						//cout << "top " << vnpercentile[ibv]*100 << "% limit: " << vnobs[ih][ic]->GetBinCenter(ib) << endl;
						found[ibv] = 1;
						vnLimit[ic][ih][ibv] = vnobs[ic][ih]->GetBinCenter(ib);
					}
				}

			}
		}
	}
}




//______________________________________________________________________________
AliJEbePercentile& AliJEbePercentile::operator=(const AliJEbePercentile& obj){
	// copy constructor
	JUNUSED(obj);
	return *this;
}


double AliJEbePercentile::GetEbeFlowPercentile(int cBin, int ih, double vn){
	double pecentile = -1;
	for(int ibv = 0 ; ibv < NBVn ; ibv++){
		if(vn <= vnLimit[cBin][ih][ibv]){
			pecentile = vnPercentile[ibv];
		}
	}
	return pecentile;
}

