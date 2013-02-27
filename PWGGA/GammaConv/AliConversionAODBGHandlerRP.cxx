#if !defined( __CINT__) || defined(__MAKECINT__)

#include <exception>
#include <iostream>
#include "AliLog.h"
#include "AliConversionAODBGHandlerRP.h"
using namespace std;
#endif


// Author Daniel Lohner (Daniel.Lohner@cern.ch)


ClassImp(AliConversionAODBGHandlerRP);

//________________________________________________________________________
AliConversionAODBGHandlerRP::AliConversionAODBGHandlerRP(Bool_t IsHeavyIon,Bool_t UseChargedTrackMult,Int_t NEvents) : TObject(),
    fIsHeavyIon(IsHeavyIon),
    fUseChargedTrackMult(UseChargedTrackMult),
    fNEvents(NEvents),
    fBGEventCounter(NULL),
    fNBGEvents(NULL),
    fNBinsZ(7),
    fNBinsMultiplicity(5+Int_t(fUseChargedTrackMult)),
    fBinLimitsArrayZ(NULL),
    fBinLimitsArrayMultiplicity(NULL),
    fBGPool(fNBinsZ,AliGammaConversionMultiplicityVector(fNBinsMultiplicity,AliGammaConversionBGEventVector(fNEvents)))
{
    // Vertex Z Binning

    fBinLimitsArrayZ = new Double_t[fNBinsZ+1];
    fBinLimitsArrayZ[0] = -50.00;
    fBinLimitsArrayZ[1] = -3.375;
    fBinLimitsArrayZ[2] = -1.605;
    fBinLimitsArrayZ[3] = -0.225;
    fBinLimitsArrayZ[4] = 1.065;
    fBinLimitsArrayZ[5] = 2.445;
    fBinLimitsArrayZ[6] = 4.245;
    fBinLimitsArrayZ[7] = 50.00;

    // MultiplicityBins

    fBinLimitsArrayMultiplicity= new Double_t[fNBinsMultiplicity+1];

    if(fUseChargedTrackMult){
        // Use Charged Particle Multiplicity

	fBinLimitsArrayMultiplicity[0] = 0;
	fBinLimitsArrayMultiplicity[1] = 8.5;
	fBinLimitsArrayMultiplicity[2] = 16.5;
	fBinLimitsArrayMultiplicity[3] = 27.5;
	fBinLimitsArrayMultiplicity[4] = 41.5;
	fBinLimitsArrayMultiplicity[5] = 200.;

	if(fIsHeavyIon){
	    fBinLimitsArrayMultiplicity[0] = 0;
	    fBinLimitsArrayMultiplicity[1] = 200.;
	    fBinLimitsArrayMultiplicity[2] = 500.;
	    fBinLimitsArrayMultiplicity[3] = 1000.;
	    fBinLimitsArrayMultiplicity[4] = 1500.;
	    fBinLimitsArrayMultiplicity[5] = 5000.;
	}
    }
    else{
	// Use V0 Multiplicity

	fBinLimitsArrayMultiplicity[0] = 2;
	fBinLimitsArrayMultiplicity[1] = 3;
	fBinLimitsArrayMultiplicity[2] = 4;
	fBinLimitsArrayMultiplicity[3] = 5;
	fBinLimitsArrayMultiplicity[4] = 9999;

	if(fIsHeavyIon){
	    fBinLimitsArrayMultiplicity[0] = 2;
	    fBinLimitsArrayMultiplicity[1] = 10;
	    fBinLimitsArrayMultiplicity[2] = 30;
	    fBinLimitsArrayMultiplicity[3] = 50;
	    fBinLimitsArrayMultiplicity[4] = 9999;
	}
    }
    Initialize();
}

//________________________________________________________________________
AliConversionAODBGHandlerRP::~AliConversionAODBGHandlerRP()
{
    if(fBinLimitsArrayZ){
	delete[] fBinLimitsArrayZ;
	fBinLimitsArrayZ=0x0;}
    if(fBinLimitsArrayMultiplicity){
	delete[] fBinLimitsArrayMultiplicity;
	fBinLimitsArrayMultiplicity=0x0;}

    if(fBGEventCounter){

	for(Int_t z=0;z<fNBinsZ;z++){
	    delete[] fBGEventCounter[z];
	}
	delete[] fBGEventCounter;
	fBGEventCounter = NULL;
    }

    // Delete pool

    for(Int_t z=0;z<fNBinsZ;z++){
	for(Int_t m=0;m<fNBinsMultiplicity;m++){
	    for(Int_t eventCounter=0;eventCounter<fNBGEvents[z][m]&&eventCounter<fNEvents;eventCounter++){

		for(UInt_t d=0;d<fBGPool[z][m][eventCounter].size();d++){
		    delete (AliAODConversionPhoton*)(fBGPool[z][m][eventCounter][d]);
		}
	    }
	}
    }

    if(fNBGEvents){
	for(Int_t z=0;z<fNBinsZ;z++){
	    delete[] fNBGEvents[z];
	}
	delete[] fNBGEvents;
	fNBGEvents = NULL;
    }

}

//________________________________________________________________________
void AliConversionAODBGHandlerRP::Initialize(){

    // Counter

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

    if(fNBGEvents == NULL){
	fNBGEvents= new Int_t*[fNBinsZ];
    }
    for(Int_t z=0;z<fNBinsZ;z++){
	fNBGEvents[z]=new Int_t[fNBinsMultiplicity];
    }
    for(Int_t z=0;z<fNBinsZ;z++){
	for(Int_t m=0;m<fNBinsMultiplicity;m++){
	    fNBGEvents[z][m]=0;
	}
    }
}

//-------------------------------------------------------------

Int_t AliConversionAODBGHandlerRP::GetZBinIndex(Double_t zvalue) const{

    if(zvalue<=fBinLimitsArrayZ[0]){
	return -1;
    }

    if(fNBinsZ<2){return 0;}

    for(Int_t i=0; i<fNBinsZ ;i++){
	if(zvalue >= fBinLimitsArrayZ[i] && zvalue <= fBinLimitsArrayZ[i+1]){
	    return i;
	}
    }
    return -1;
}
//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetMultiplicityBinIndex(Int_t multiplicity) const{
  if(fNBinsMultiplicity<2){
    return 0;
  }

  for(Int_t i=0; i<fNBinsMultiplicity ;i++){
    if(multiplicity >= fBinLimitsArrayMultiplicity[i] && multiplicity < fBinLimitsArrayMultiplicity[i+1]){
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------
Bool_t AliConversionAODBGHandlerRP::FindBins(TObjArray * const eventGammas,AliVEvent *fInputEvent,Int_t &zbin,Int_t &mbin){
    Double_t vertexz=fInputEvent->GetPrimaryVertex()->GetZ();
    zbin = GetZBinIndex(vertexz);

    Int_t multiplicity=0;
    if(fUseChargedTrackMult){
	multiplicity=fInputEvent->GetNumberOfTracks();
    }
    else{
	multiplicity=eventGammas->GetEntries();
    }
    mbin = GetMultiplicityBinIndex(multiplicity);

    if(zbin<fNBinsZ&&mbin<fNBinsMultiplicity){
	if(zbin>=0&&mbin>=0){
	    return kTRUE;
	}
    }
    //cout<<Form("Requested BG pool does not exist:  z %i m %i",zbin,mbin)<<endl;
    return kFALSE;
}
//-------------------------------------------------------------
Bool_t AliConversionAODBGHandlerRP::FindBins(TList * const eventGammas,AliVEvent *fInputEvent,Int_t &zbin,Int_t &mbin){
    Double_t vertexz=fInputEvent->GetPrimaryVertex()->GetZ();
    zbin = GetZBinIndex(vertexz);

    Int_t multiplicity=0;
    if(fUseChargedTrackMult){
	multiplicity=fInputEvent->GetNumberOfTracks();
    }
    else{
	multiplicity=eventGammas->GetEntries();
    }
    mbin = GetMultiplicityBinIndex(multiplicity);

    if(zbin<fNBinsZ&&mbin<fNBinsMultiplicity){
	if(zbin>=0&&mbin>=0){
	    return kTRUE;
	}
    }
    //cout<<Form("Requested BG pool does not exist:  z %i m %i",zbin,mbin)<<endl;
    return kFALSE;
}


//-------------------------------------------------------------
void AliConversionAODBGHandlerRP::AddEvent(TObjArray * const eventGammas,AliVEvent *fInputEvent){

    if(eventGammas->GetEntriesFast()==0)return;

    Int_t z;
    Int_t m;
    if(FindBins(eventGammas,fInputEvent,z,m)){
  
   
	// If Event Stack is full, replace the first entry (First in first out)
	if(fBGEventCounter[z][m] >= fNEvents){
	    fBGEventCounter[z][m]=0;
	}

	// Update number of Events stored
	if(fNBGEvents[z][m] < fNEvents){
	    fNBGEvents[z][m]++;
	}

	Int_t eventCounter=fBGEventCounter[z][m];

	//clear the vector for old gammas
	for(UInt_t d=0;d<fBGPool[z][m][eventCounter].size();d++){
	    delete (AliAODConversionPhoton*)(fBGPool[z][m][eventCounter][d]);
	}

	fBGPool[z][m][eventCounter].clear();

	// add the gammas to the vector

	for(Int_t i=0; i< eventGammas->GetEntriesFast();i++){
	    fBGPool[z][m][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
	}

	fBGEventCounter[z][m]++;
    }
}
//-------------------------------------------------------------
void AliConversionAODBGHandlerRP::AddEvent(TList * const eventGammas,AliVEvent *fInputEvent){
    if(eventGammas->GetEntries()==0)return;

    Int_t z;
    Int_t m;
    if(FindBins(eventGammas,fInputEvent,z,m)){
  
   
	// If Event Stack is full, replace the first entry (First in first out)
	if(fBGEventCounter[z][m] >= fNEvents){
	    fBGEventCounter[z][m]=0;
	}

	// Update number of Events stored
	if(fNBGEvents[z][m] < fNEvents){
	    fNBGEvents[z][m]++;
	}

	Int_t eventCounter=fBGEventCounter[z][m];

	//clear the vector for old gammas
	for(UInt_t d=0;d<fBGPool[z][m][eventCounter].size();d++){
	    delete (AliAODConversionPhoton*)(fBGPool[z][m][eventCounter][d]);
	}

	fBGPool[z][m][eventCounter].clear();

	// add the gammas to the vector

	for(Int_t i=0; i< eventGammas->GetEntries();i++){
	    fBGPool[z][m][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
	}

	fBGEventCounter[z][m]++;
    }
}

//-------------------------------------------------------------
AliGammaConversionPhotonVector* AliConversionAODBGHandlerRP::GetBGGoodGammas(TObjArray * const eventGammas,AliVEvent *fInputEvent,Int_t event){
    Int_t zbin;
    Int_t mbin;
    if(FindBins(eventGammas,fInputEvent,zbin,mbin)){
	return &(fBGPool[zbin][mbin][event]);
    }
    return NULL;
}
//-------------------------------------------------------------
AliGammaConversionPhotonVector* AliConversionAODBGHandlerRP::GetBGGoodGammas(TList * const eventGammas,AliVEvent *fInputEvent,Int_t event){
    Int_t zbin;
    Int_t mbin;
    if(FindBins(eventGammas,fInputEvent,zbin,mbin)){
	return &(fBGPool[zbin][mbin][event]);
    }
    return NULL;
}

//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetNBGEvents(TObjArray * const eventGammas,AliVEvent *fInputEvent){
    Int_t zbin;
    Int_t mbin;
    if(FindBins(eventGammas,fInputEvent,zbin,mbin)){
	return fNBGEvents[zbin][mbin];
    }
    return 0;
}
//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetNBGEvents(TList * const eventGammas,AliVEvent *fInputEvent){
   Int_t zbin;
   Int_t mbin;
   if(FindBins(eventGammas,fInputEvent,zbin,mbin)){
      return fNBGEvents[zbin][mbin];
   }
   return 0;
}

