/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Martin Wilde, Friederike Bock, Daniel Lohner, Svein Lindal     *
 * based on previous version by Kenneth Aamodt and Ana Marin			  *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class for handling of background calculation
//---------------------------------------------
////////////////////////////////////////////////

#include "AliGammaConversionAODBGHandler.h"
#include "AliKFParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"

using namespace std;

ClassImp(AliGammaConversionAODBGHandler)

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler() :
	TObject(),
	fNEvents(10),
	fBGEventCounter(NULL),
	fBGEventENegCounter(NULL),
	fBGEventMesonCounter(NULL),
	fBGProbability(NULL),
	fBGEventVertex(NULL),
	fNBinsZ(0),
	fNBinsMultiplicity(0),
	fBinLimitsArrayZ(NULL),
	fBinLimitsArrayMultiplicity(NULL),
	fBGEvents(),
	fBGEventsENeg(),
	fBGEventsMeson()
{
	// constructor
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler(Int_t binsZ,Int_t binsMultiplicity,Int_t nEvents) :
	TObject(),
	fNEvents(nEvents),
	fBGEventCounter(NULL),
	fBGEventENegCounter(NULL),
	fBGEventMesonCounter(NULL),
	fBGProbability(NULL),
	fBGEventVertex(NULL),
	fNBinsZ(binsZ),
	fNBinsMultiplicity(binsMultiplicity),
	fBinLimitsArrayZ(NULL),
	fBinLimitsArrayMultiplicity(NULL),
	fBGEvents(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents))),
	fBGEventsENeg(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents))),
	fBGEventsMeson(binsZ,AliGammaConversionMotherMultipicityVector(binsMultiplicity,AliGammaConversionMotherBGEventVector(nEvents)))
{
	// constructor
}


//_____________________________________________________________________________________________________________________________
AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler(Int_t collisionSystem, Int_t centMin, Int_t centMax,
                                                               Int_t nEvents, Bool_t useTrackMult, Int_t mode, Int_t binsZ, Int_t binsMultiplicity) :
	TObject(),
	fNEvents(nEvents),
	fBGEventCounter(NULL),
	fBGEventENegCounter(NULL),
	fBGEventMesonCounter(NULL),
	fBGProbability(NULL),
	fBGEventVertex(NULL),
	fNBinsZ(binsZ),
	fNBinsMultiplicity(binsMultiplicity),
	fBinLimitsArrayZ(NULL),
	fBinLimitsArrayMultiplicity(NULL),
	fBGEvents(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents))),
	fBGEventsENeg(binsZ,AliGammaConversionMultipicityVector(binsMultiplicity,AliGammaConversionBGEventVector(nEvents))),
	fBGEventsMeson(binsZ,AliGammaConversionMotherMultipicityVector(binsMultiplicity,AliGammaConversionMotherBGEventVector(nEvents)))
{
	// constructor
    if(fNBinsZ>8) fNBinsZ = 8;
    if(fNBinsMultiplicity>5) fNBinsMultiplicity = 5;

	// Initializing z vertex bins
	fBinLimitsArrayZ= new Double_t[fNBinsZ] ;
	if(collisionSystem > 0 && collisionSystem < 8){ // PbPb
		Double_t fBinLimitsArrayZPbPb[8] = 	{-50, 	-5.5, 	-2.9, 	-0.65,
											 1.45, 	3.65, 	6.15, 	50};
		for (Int_t i = 0; i < fNBinsZ; i++){
			fBinLimitsArrayZ[i] =  fBinLimitsArrayZPbPb[i];
		}	
	} else if(collisionSystem == 0){				// pp
		Double_t fBinLimitsArrayZpp[8] = 	{-50, 	-3.375, -1.605, -0.225, 
											 1.065, 2.445, 	4.245, 	50};
		for (Int_t i = 0; i < fNBinsZ; i++){
			fBinLimitsArrayZ[i] =  fBinLimitsArrayZpp[i];
		}	
	} else { 										// pPb
		Double_t fBinLimitsArrayZpPb[8] = 	{-50, 	-5.85, 	-3.35, 	-1.15, 
											 0.85, 	2.95, 	5.55, 	50};
		for (Int_t i = 0; i < fNBinsZ; i++){
			fBinLimitsArrayZ[i] =  fBinLimitsArrayZpPb[i];
		}	
	}

	// Initializing multiplicity bins 
	fBinLimitsArrayMultiplicity= new Double_t[fNBinsMultiplicity];
	if(useTrackMult){ // multiplicity binning based on number of good global tracks
		// default pp values
		Double_t fBinLimitsArrayMultiplicitypp[5] = 	{0., 	8.5, 	16.5, 	27.5, 	200.};
		for (Int_t i = 0; i < fNBinsMultiplicity; i++){
			fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypp[i];
		}	
		if(collisionSystem > 0 && collisionSystem < 8){ // PbPb values
			if(centMin == 0 && centMax == 5){ // 0-5% central
				Double_t fBinLimitsArrayMultiplicityPbPb0005[5] = 	{0., 1540., 1665., 1780., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0005[i];
				}	
			} else if(centMin == 0 && centMax == 10){ // 0-10% central
				Double_t fBinLimitsArrayMultiplicityPbPb0010[5] = 	{0., 1360., 1520., 1685., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0010[i];
				}	
			} else if(centMin == 0 && centMax == 20){ // 0-20% central
				Double_t fBinLimitsArrayMultiplicityPbPb0020[5] = 	{0., 1110., 1360., 1600., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0020[i];
				}	
			} else if(centMin == 0 && centMax == 80){ // 0-80% central
				Double_t fBinLimitsArrayMultiplicityPbPb0080[5] = 	{0., 890., 1240., 1540., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0080[i];
				}	
			} else if(centMin == 5 && centMax == 10){ // 5-10% central
				Double_t fBinLimitsArrayMultiplicityPbPb0510[5] = 	{0., 1250., 1345., 1445., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0510[i];
				}	
			} else if(centMin == 10 && centMax == 20){ // 10-20% central
				Double_t fBinLimitsArrayMultiplicityPbPb1020[5] = 	{0., 915., 1020., 1130., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb1020[i];
				}	
			} else if(centMin == 20 && centMax == 40){ // 20-40% central
				Double_t fBinLimitsArrayMultiplicityPbPb2040[5] = 	{0., 510., 625., 730., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb2040[i];
				}	
			} else if(centMin == 40 && centMax == 80){ // 40-80% central
				Double_t fBinLimitsArrayMultiplicityPbPb4080[5] = 	{0., 185., 250., 300., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb4080[i];
				}	
			} else if(centMin == 60 && centMax == 80){ // 60-80% central
				Double_t fBinLimitsArrayMultiplicityPbPb6080[5] = 	{0., 55., 80., 100., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb6080[i];
				}	
			} else { // all other centrality classes 
				Double_t fBinLimitsArrayMultiplicityPbPb[5] = 	{0., 510., 625., 730., 5000};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb[i];
				}	
			}
		} else if(collisionSystem == 8 || collisionSystem == 9){ // pPb
			Double_t fBinLimitsArrayMultiplicitypPb[5] = 	{0., 7.5, 16.5, 29.5, 500};
			for (Int_t i = 0; i < fNBinsMultiplicity; i++){
				fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypPb[i];
			}	
			if(centMin == 0 && centMax == 20){ // pPb 0-20 %
				Double_t fBinLimitsArrayMultiplicitypPb0020[5] = 	{0., 31.5, 40.5, 50.5, 500};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypPb0020[i];
				}	
			} else if(centMin == 20 && centMax == 40){
				Double_t fBinLimitsArrayMultiplicitypPb2040[5] = 	{0., 19.5, 25.5, 32.5, 500};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypPb2040[i];
				}	
			} else if(centMin == 40 && centMax == 60){
				Double_t fBinLimitsArrayMultiplicitypPb4060[5] = 	{0., 12.5, 16.5, 22.5, 500};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypPb4060[i];
				}	
			} else if(centMin == 60 && centMax == 80){ 
				Double_t fBinLimitsArrayMultiplicitypPb6080[5] = 	{0., 5.5, 9.5, 13.5, 500};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypPb6080[i];
				}	
			} else if(centMin == 60 && centMax == 100){
				Double_t fBinLimitsArrayMultiplicitypPb60100[5] = 	{0., 2.5, 6.5, 11.5, 500};
				for (Int_t i = 0; i < fNBinsMultiplicity; i++){
					fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypPb60100[i];
				}	
			}  
		}
	} else{	// Initializing Multiplicity binning with photon Mult 
		if (mode == 0 || mode == 1) { // settings for Conv-Conv && Conv-Dalitz
			// pp & pPb defaults
			Double_t fBinLimitsArrayMultiplicitypp[5] = 	{2., 3., 4., 5., 9999.};
			for (Int_t i = 0; i < fNBinsMultiplicity; i++){
				fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypp[i];
			}	
			if(collisionSystem > 0 && collisionSystem < 8){ // settings PbPb
				if(centMin == 0 && centMax == 5){ 			// 0-5% 
					Double_t fBinLimitsArrayMultiplicityPbPb0005[5] = 	{0., 27., 31., 36., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0005[i];
					}	
				} else if(centMin == 0 && centMax == 10){	// 0-10%
					Double_t fBinLimitsArrayMultiplicityPbPb0010[5] = 	{0., 25., 30., 36., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0010[i];
					}	
				} else if(centMin == 0 && centMax == 20){	// 0-20%
					Double_t fBinLimitsArrayMultiplicityPbPb0020[5] = 	{0., 22., 27., 33., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0020[i];
					}	
				} else if(centMin == 0 && centMax == 80){	// 0-80%
					Double_t fBinLimitsArrayMultiplicityPbPb0080[5] = 	{0., 18., 25., 32., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0080[i];
					}	
				} else if(centMin == 5 && centMax == 10){ 	// 5-10%
					Double_t fBinLimitsArrayMultiplicityPbPb0510[5] = 	{0., 23., 27., 32., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0510[i];
					}	
				} else if(centMin == 10 && centMax == 20){	//10-20%
					Double_t fBinLimitsArrayMultiplicityPbPb1020[5] = 	{0., 18., 22., 27., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb1020[i];
					}	
				} else if(centMin == 20 && centMax == 40){	// 20-40%
					Double_t fBinLimitsArrayMultiplicityPbPb2040[5] = 	{0., 11., 14., 18., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb2040[i];
					}	
				} else if(centMin == 40 && centMax == 80){ // 40-80%
					Double_t fBinLimitsArrayMultiplicityPbPb4080[5] = 	{0., 5., 7., 11., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb4080[i];
					}	
				} else if(centMin == 60 && centMax == 80){ // 60-80%
					Double_t fBinLimitsArrayMultiplicityPbPb6080[5] = 	{0., 2., 3., 5., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb6080[i];
					}	
				} else{ // default PbPb
					Double_t fBinLimitsArrayMultiplicityPbPb[5] = 	{0., 11., 14., 18., 100.};
					for (Int_t i = 0; i < fNBinsMultiplicity; i++){
						fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb[i];
					}	
				}
			}
		} else if (mode == 2 || mode == 3 || mode == 4 || mode == 5){ // setting for EMCAL-Conv, PHOS-Conv, EMCAL-EMCAL, PHOS-PHOS
            if(collisionSystem > 0 && collisionSystem < 8){ // settings PbPb
                if(centMin == 0 && centMax == 5){ 			// 0-5%
                    Double_t fBinLimitsArrayMultiplicityPbPb0005[5] = 	{0., 27., 31., 36., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0005[i];
                    }
                } else if(centMin == 0 && centMax == 10){	// 0-10%
                    Double_t fBinLimitsArrayMultiplicityPbPb0010[5] = 	{0., 25., 30., 36., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0010[i];
                    }
                } else if(centMin == 0 && centMax == 20){	// 0-20%
                    Double_t fBinLimitsArrayMultiplicityPbPb0020[5] = 	{0., 22., 27., 33., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0020[i];
                    }
                } else if(centMin == 0 && centMax == 80){	// 0-80%
                    Double_t fBinLimitsArrayMultiplicityPbPb0080[5] = 	{0., 18., 25., 32., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0080[i];
                    }
                } else if(centMin == 5 && centMax == 10){ 	// 5-10%
                    Double_t fBinLimitsArrayMultiplicityPbPb0510[5] = 	{0., 23., 27., 32., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb0510[i];
                    }
                } else if(centMin == 10 && centMax == 20){	//10-20%
                    Double_t fBinLimitsArrayMultiplicityPbPb1020[5] = 	{0., 18., 22., 27., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb1020[i];
                    }
                } else if(centMin == 20 && centMax == 40){	// 20-40%
                    Double_t fBinLimitsArrayMultiplicityPbPb2040[5] = 	{0., 11., 14., 18., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb2040[i];
                    }
                } else if(centMin == 40 && centMax == 80){ // 40-80%
                    Double_t fBinLimitsArrayMultiplicityPbPb4080[5] = 	{0., 5., 7., 11., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb4080[i];
                    }
                } else if(centMin == 60 && centMax == 80){ // 60-80%
                    Double_t fBinLimitsArrayMultiplicityPbPb6080[5] = 	{0., 2., 3., 5., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb6080[i];
                    }
                } else{ // default PbPb
                    Double_t fBinLimitsArrayMultiplicityPbPb[5] = 	{0., 11., 14., 18., 100.};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicityPbPb[i];
                    }
                }
            }
            else
            {
                //seperate settings for ConvCalo and CaloCalo for pp/pPb
                if (mode == 2 || mode == 3){ //ConvCalo
                    Double_t fBinLimitsArrayMultiplicitypp_pPbConvCalo[5] = {1., 2., 3., 4., 9999};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypp_pPbConvCalo[i];
                    }
                }
                else{ //CaloCalo
                    Double_t fBinLimitsArrayMultiplicitypp_pPbCaloCalo[5] = {2., 3., 4., 5., 9999};
                    for (Int_t i = 0; i < fNBinsMultiplicity; i++){
                        fBinLimitsArrayMultiplicity[i] =  fBinLimitsArrayMultiplicitypp_pPbCaloCalo[i];
                    }
                }
            }
		}	
	} 
	
	Initialize(fBinLimitsArrayZ,fBinLimitsArrayMultiplicity);
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODBGHandler::AliGammaConversionAODBGHandler(const AliGammaConversionAODBGHandler & original) :
	TObject(original),
	fNEvents(original.fNEvents),
	fBGEventCounter(original.fBGEventCounter),
	fBGEventENegCounter(original.fBGEventENegCounter),
	fBGEventMesonCounter(original.fBGEventMesonCounter),
	fBGProbability(original.fBGProbability),
	fBGEventVertex(original.fBGEventVertex),
	fNBinsZ(original.fNBinsZ),
	fNBinsMultiplicity(original.fNBinsMultiplicity),
	fBinLimitsArrayZ(original.fBinLimitsArrayZ),
	fBinLimitsArrayMultiplicity(original.fBinLimitsArrayMultiplicity),
	fBGEvents(original.fBGEvents),
	fBGEventsENeg(original.fBGEventsENeg),
	fBGEventsMeson(original.fBGEventsMeson)
{
	//copy constructor	
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODBGHandler & AliGammaConversionAODBGHandler::operator = (const AliGammaConversionAODBGHandler & /*source*/)
{
	// assignment operator
	return *this;
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODBGHandler::~AliGammaConversionAODBGHandler(){

	if(fBGEventCounter){
		for(Int_t z=0;z<fNBinsZ;z++){
			delete[] fBGEventCounter[z];
		}
		delete[] fBGEventCounter;
		fBGEventCounter = NULL;
	}

	if(fBGEventVertex){
		for(Int_t z=0;z<fNBinsZ;z++){
			for(Int_t m=0;m<fNBinsMultiplicity;m++){
				delete [] fBGEventVertex[z][m];
			}
			delete [] fBGEventVertex[z];
		}
		delete [] fBGEventVertex;
	}

	if(fBGEventENegCounter){
		for(Int_t z=0;z<fNBinsZ;z++){
			delete[] fBGEventENegCounter[z];
		}
		delete[] fBGEventENegCounter;
		fBGEventENegCounter = NULL;
	}

	if(fBGEventMesonCounter){
		for(Int_t z=0;z<fNBinsZ;z++){
			delete[] fBGEventMesonCounter[z];
		}
		delete[] fBGEventMesonCounter;
		fBGEventMesonCounter = NULL;
	}

	if(fBinLimitsArrayZ){
		delete[] fBinLimitsArrayZ;
	}

	if(fBinLimitsArrayMultiplicity){
		delete[] fBinLimitsArrayMultiplicity;
	}
}

//_____________________________________________________________________________________________________________________________
void AliGammaConversionAODBGHandler::Initialize(Double_t * const zBinLimitsArray, Double_t * const multiplicityBinLimitsArray){
  // see header file for documantation  

	if(zBinLimitsArray){
		fBinLimitsArrayZ = zBinLimitsArray;
	}
	else{
		//Print warning
	}
	
	if(multiplicityBinLimitsArray){
		fBinLimitsArrayMultiplicity = multiplicityBinLimitsArray ;
	}
	else{
		//Print warning
	}
	if(fBGEventCounter == NULL){
		fBGEventCounter= new Int_t*[fNBinsZ];
	}
	for(Int_t z=0;z<fNBinsZ;z++){
		fBGEventCounter[z]=new Int_t[fNBinsMultiplicity];
	}
	
	for(Int_t z=0;z<fNBinsZ;z++){
		for(Int_t m=0;m<fNBinsMultiplicity;m++){
			fBGEventCounter[z][m]=0;
		}
	}

	if(fBGEventMesonCounter == NULL){
		fBGEventMesonCounter= new Int_t*[fNBinsZ];
	}
	for(Int_t z=0;z<fNBinsZ;z++){
		fBGEventMesonCounter[z]=new Int_t[fNBinsMultiplicity];
	}
	
	for(Int_t z=0;z<fNBinsZ;z++){
		for(Int_t m=0;m<fNBinsMultiplicity;m++){
			fBGEventMesonCounter[z][m]=0;
		}
	}

	
	if(fBGEventVertex == NULL){
		fBGEventVertex = new GammaConversionVertex**[fNBinsZ];
	}
	for(Int_t z=0; z < fNBinsZ; z++){
		fBGEventVertex[z]= new GammaConversionVertex*[fNBinsMultiplicity];
	}
	for(Int_t z=0;z<fNBinsZ;z++){
		for(Int_t m=0;m<fNBinsMultiplicity; m++){
			fBGEventVertex[z][m]= new GammaConversionVertex[fNEvents];
		}
	}
	if( fBGEventENegCounter == NULL){
		fBGEventENegCounter = new Int_t*[fNBinsZ];
	}

	for(Int_t z=0; z < fNBinsZ; z++){
		fBGEventENegCounter[z] = new Int_t[fNBinsMultiplicity];
	}

	for(Int_t z=0;z<fNBinsZ;z++){
		for(Int_t m=0;m<fNBinsMultiplicity; m++){
			fBGEventENegCounter[z][m] = 0;
		}
	}

	if(fBGProbability == NULL){
		fBGProbability = new Double_t*[fNBinsZ];
	}
	for(Int_t z=0; z < fNBinsZ; z++){
		fBGProbability[z] = new Double_t[fNBinsMultiplicity];
	}
    Double_t BGProbabilityLookup[7][4] =
           {
             {0.243594,0.279477,0.305104,0.315927},
             {0.241964,0.272995,0.307165,0.292248},
             {0.241059,0.27509,0.283657,0.310512},
             {0.23888,0.283418,0.297232,0.348188},
             {0.245555,0.281218,0.317236,0.323495},
             {0.244572,0.259498,0.278383,0.284696},
             {0.24703, 0.275265,0.284004,0.343584}
           };
	for(Int_t z=0;z<fNBinsZ;z++){
		for(Int_t m=0;m<fNBinsMultiplicity; m++){
            if((z<7)&&(m<4)){
                fBGProbability[z][m] = BGProbabilityLookup[z][m];
            }else{
                fBGProbability[z][m] = 1;
            }
		}
	}
	
}

//_____________________________________________________________________________________________________________________________
Int_t AliGammaConversionAODBGHandler::GetZBinIndex(Double_t zvalue) const{
	// see header file for documantation
	if(fNBinsZ<2 || zvalue<=fBinLimitsArrayZ[0]){
		return 0;
	}

	for(Int_t i=0; i<fNBinsZ-1 ;i++){
		if(zvalue >= fBinLimitsArrayZ[i] && zvalue <= fBinLimitsArrayZ[i+1]){
			return i;
		}
	}
	return fNBinsZ-1;
}

//_____________________________________________________________________________________________________________________________
Int_t AliGammaConversionAODBGHandler::GetMultiplicityBinIndex(Int_t multiplicity) const{
	// see header file for documantation  
	if(fNBinsMultiplicity<2){
		return 0;
	}

	for(Int_t i=0; i<fNBinsMultiplicity-1 ;i++){
		if(multiplicity >= fBinLimitsArrayMultiplicity[i] && multiplicity < fBinLimitsArrayMultiplicity[i+1]){
			return i;
		}
	}
	return fNBinsMultiplicity-1;
}

//_____________________________________________________________________________________________________________________________
void AliGammaConversionAODBGHandler::AddEvent(TList* const eventGammas,Double_t xvalue, Double_t yvalue, Double_t zvalue, Int_t multiplicity, Double_t epvalue){

	// see header file for documantation  

	//  cout<<"Entering the AddEvent function"<<endl;

	Int_t z = GetZBinIndex(zvalue);
	Int_t m = GetMultiplicityBinIndex(multiplicity);

	if(fBGEventCounter[z][m] >= fNEvents){
		fBGEventCounter[z][m]=0;
	}
	Int_t eventCounter=fBGEventCounter[z][m];
	
	/*
	if(fBGEventVertex[z][m][eventCounter]){
		delete fBGEventVertex[z][m][eventCounter];
	}
	*/
	fBGEventVertex[z][m][eventCounter].fX = xvalue;
	fBGEventVertex[z][m][eventCounter].fY = yvalue;
	fBGEventVertex[z][m][eventCounter].fZ = zvalue;
	fBGEventVertex[z][m][eventCounter].fEP = epvalue;

	//first clear the vector
	// cout<<"Size of vector: "<<fBGEvents[z][m][eventCounter].size()<<endl;
	//  cout<<"Checking the entries: Z="<<z<<", M="<<m<<", eventCounter="<<eventCounter<<endl;

	//  cout<<"The size of this vector is: "<<fBGEvents[z][m][eventCounter].size()<<endl;
    for(Int_t d=0;d<fBGEvents[z][m][eventCounter].size();d++){
		delete (AliAODConversionPhoton*)(fBGEvents[z][m][eventCounter][d]);
	}
	fBGEvents[z][m][eventCounter].clear();
	
	// add the gammas to the vector
	for(Int_t i=0; i< eventGammas->GetEntries();i++){
		//    AliKFParticle *t = new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i)));
		fBGEvents[z][m][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
	}
	fBGEventCounter[z][m]++;
}

//_____________________________________________________________________________________________________________________________
void AliGammaConversionAODBGHandler::AddMesonEvent(TList* const eventMothers, Double_t xvalue, Double_t yvalue, Double_t zvalue, Int_t multiplicity, Double_t epvalue){

	// see header file for documantation  
	//  cout<<"Entering the AddEvent function"<<endl;
	Int_t z = GetZBinIndex(zvalue);
	Int_t m = GetMultiplicityBinIndex(multiplicity);

	if(fBGEventMesonCounter[z][m] >= fNEvents){
		fBGEventMesonCounter[z][m]=0;
	}
	Int_t eventCounter=fBGEventMesonCounter[z][m];
	
	fBGEventVertex[z][m][eventCounter].fX = xvalue;
	fBGEventVertex[z][m][eventCounter].fY = yvalue;
	fBGEventVertex[z][m][eventCounter].fZ = zvalue;
	fBGEventVertex[z][m][eventCounter].fEP = epvalue;

	//first clear the vector
  for(Int_t d=0;d<fBGEventsMeson[z][m][eventCounter].size();d++){
		delete (AliAODConversionMother*)(fBGEventsMeson[z][m][eventCounter][d]);
	}
	fBGEventsMeson[z][m][eventCounter].clear();
	
	// add the gammas to the vector
	for(Int_t i=0; i< eventMothers->GetEntries();i++){
		fBGEventsMeson[z][m][eventCounter].push_back(new AliAODConversionMother(*(AliAODConversionMother*)(eventMothers->At(i))));
	}
	fBGEventMesonCounter[z][m]++;
}

void AliGammaConversionAODBGHandler::AddMesonEvent(const std::vector<AliAODConversionMother> &eventMother, Double_t xvalue, Double_t yvalue, Double_t zvalue, Int_t multiplicity, Double_t epvalue){
  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetMultiplicityBinIndex(multiplicity);

  if(fBGEventMesonCounter[z][m] >= fNEvents){
    fBGEventMesonCounter[z][m]=0;
  }
  Int_t eventCounter=fBGEventMesonCounter[z][m];

  fBGEventVertex[z][m][eventCounter].fX = xvalue;
  fBGEventVertex[z][m][eventCounter].fY = yvalue;
  fBGEventVertex[z][m][eventCounter].fZ = zvalue;
  fBGEventVertex[z][m][eventCounter].fEP = epvalue;

  //first clear the vector
  for(Int_t d=0;d<fBGEvents[z][m][eventCounter].size();d++){
    delete (AliAODConversionMother*)(fBGEventsMeson[z][m][eventCounter][d]);
  }
  fBGEventsMeson[z][m][eventCounter].clear();

  // add the gammas to the vector
  for(const auto &mother : eventMother){
    fBGEventsMeson[z][m][eventCounter].push_back(new AliAODConversionMother(mother));
  }
  fBGEventMesonCounter[z][m]++;
}

//_____________________________________________________________________________________________________________________________
void AliGammaConversionAODBGHandler::AddElectronEvent(TClonesArray* const eventENeg, Double_t zvalue, Int_t multiplicity){

	Int_t z = GetZBinIndex(zvalue);
	Int_t m = GetMultiplicityBinIndex(multiplicity);

	if(fBGEventENegCounter[z][m] >= fNEvents){
		fBGEventENegCounter[z][m]=0;
	}
	Int_t eventENegCounter=fBGEventENegCounter[z][m];
	
	//first clear the vector
	// cout<<"Size of vector: "<<fBGEvents[z][m][eventCounter].size()<<endl;
	//  cout<<"Checking the entries: Z="<<z<<", M="<<m<<", eventCounter="<<eventCounter<<endl;

	//  cout<<"The size of this vector is: "<<fBGEvents[z][m][eventCounter].size()<<endl;
    for(Int_t d=0;d<fBGEventsENeg[z][m][eventENegCounter].size();d++){
		delete (AliAODConversionPhoton*)(fBGEventsENeg[z][m][eventENegCounter][d]);
	}

	fBGEventsENeg[z][m][eventENegCounter].clear();

	// add the electron to the vector
	for(Int_t i=0; i< eventENeg->GetEntriesFast();i++){
		//    AliKFParticle *t = new AliKFParticle(*(AliKFParticle*)(eventGammas->At(i)));
		fBGEventsENeg[z][m][eventENegCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventENeg->At(i))));
	}
	fBGEventENegCounter[z][m]++;
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODVector* AliGammaConversionAODBGHandler::GetBGGoodV0s(Int_t zbin, Int_t mbin, Int_t event){
	//see headerfile for documentation
	return &(fBGEvents[zbin][mbin][event]);
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionMotherAODVector* AliGammaConversionAODBGHandler::GetBGGoodMesons(Int_t zbin, Int_t mbin, Int_t event){
	//see headerfile for documentation
	return &(fBGEventsMeson[zbin][mbin][event]);
}

//_____________________________________________________________________________________________________________________________
Int_t AliGammaConversionAODBGHandler::GetNBackgroundEventsInBuffer(Int_t binz, int binMult) const {
  return fBGEventsMeson[binz][binMult].size();
}

//_____________________________________________________________________________________________________________________________
AliGammaConversionAODVector* AliGammaConversionAODBGHandler::GetBGGoodENeg(Int_t event, Double_t zvalue, Int_t multiplicity){
	//see headerfile for documentation
	Int_t z = GetZBinIndex(zvalue);
	Int_t m = GetMultiplicityBinIndex(multiplicity);
	return &(fBGEventsENeg[z][m][event]);
}

//_____________________________________________________________________________________________________________________________
void AliGammaConversionAODBGHandler::PrintBGArray(){
	//see headerfile for documentation
	for(Int_t z=0;z<fNBinsZ;z++){
		if(z==2){
			cout<<"Getting the data for z bin: "<<z<<endl;
			for(Int_t multiplicity=0;multiplicity<fNBinsMultiplicity;multiplicity++){
				if(multiplicity==2){
					cout<<"Getting the data for multiplicity bin: "<<multiplicity<<endl;	
					for(Int_t event=0;event<fNEvents;event++){
						if(fBGEvents[z][multiplicity][event].size()>0){
						cout<<"Event: "<<event<<" has: "<<fBGEvents[z][multiplicity][event].size()<<endl;
						}
					}
				}
			}
		}
	}
}
