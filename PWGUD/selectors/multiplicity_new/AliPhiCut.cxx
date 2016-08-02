#include "AliPhiCut.h"
#include <TMath.h>
#include <cstdio>
#include <iostream>

ClassImp ( AliPhiCut )

const Double_t phiconv = TMath::Pi() / 180. ;

AliPhiCut::AliPhiCut( Int_t energy ) : TObject() {
	if ( (energy > -1 ) && ( energy < 4 ) ) CreateDefaultCut( energy ); else std::cout<<"Empty cut Initialized"<<std::endl;
}

AliPhiCut::~AliPhiCut() {

}

void AliPhiCut::Print ( Option_t* option ) const {
	TObject::Print(option);
	std::cout<<"AliPhiCut:"<<std::endl;
	for (Int_t i = 0 ; i < edges.GetSize(); ++i) {
		std::cout<<edges.At(i) / phiconv <<" - "<< ( edges.At(i) + widths.At(i) ) / phiconv<<std::endl;
	}
}

AliPhiCut::AliPhiCut ( Double_t* e, Double_t* w, Int_t n ) {
	if ( e && w && ( n > 0 ) ) {
		edges.Set(n, e);
		widths.Set(n, w);
	} else {
		std::cerr<<"Input data incomplete, cut not initialized"<<std::endl;
	}
}

Bool_t AliPhiCut::CheckCut ( const Double_t phi ) {
	for (Int_t i = 0 ; i < edges.GetSize(); ++i) {
		if ( ( phi >= edges.At(i) ) && ( phi < edges.At(i) + widths.At(i)) ) {
			return kTRUE;
		}
	}
	return kFALSE;
}

void AliPhiCut::CreateDefaultCut ( Int_t energy ) {
	switch (energy) {
		case 0:	{			// 0.9 TeV
				const Int_t ncut = 3;
				Double_t e[ncut] = {0, 135 * phiconv, 205 * phiconv};
				Double_t w[ncut] = {90 * phiconv, 25 * phiconv, 15 * phiconv};
				edges.Set(ncut, e);
				widths.Set(ncut, w);
				std::cout<<"Initialized phi cut for 0.9 TeV data"<<std::endl;
			}
			break;
		case 1:	{			// 2.76 TeV
				const Int_t ncut = 4;
				Double_t e[ncut] = {0, 135 * phiconv, 205 * phiconv, 295 * phiconv};
				Double_t w[ncut] = {70 * phiconv, 13 * phiconv, 10 * phiconv, 25 * phiconv};
				edges.Set(ncut, e);
				widths.Set(ncut, w);
				std::cout<<"Initialized phi cut for 2.76 TeV data"<<std::endl;
			}
			break;
		case 2:	{			// 7 TeV
				const Int_t ncut = 4;
				Double_t e[ncut] = {5, 135 * phiconv, 205 * phiconv, 295 * phiconv};
				Double_t w[ncut] = {65 * phiconv, 25 * phiconv, 10 * phiconv, 25 * phiconv};
				edges.Set(ncut, e);
				widths.Set(ncut, w);
				std::cout<<"Initialized phi cut for 7 TeV data"<<std::endl;
			}
			break;
		case 3:	{			// 8 TeV
				const Int_t ncut = 5;
				Double_t e[ncut] = {0, 135 * phiconv, 205 * phiconv, 255 * phiconv, 325 * phiconv};
				Double_t w[ncut] = {90 * phiconv, 50 * phiconv, 30 * phiconv, 45 * phiconv, 15 * phiconv};
				edges.Set(ncut, e);
				widths.Set(ncut, w);
				std::cout<<"Initialized phi cut for 8 TeV data"<<std::endl;
			}
			break;
		default:
			std::cerr<<"Energy not set, skipping cut"<<std::endl;
			break;
	}
}
