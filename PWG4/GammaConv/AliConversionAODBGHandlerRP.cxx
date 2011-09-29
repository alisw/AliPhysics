#if !defined( __CINT__) || defined(__MAKECINT__)

#include <exception>
#include <iostream>
#include "AliLog.h"
#include "AliConversionAODBGHandlerRP.h"
using namespace std;
#endif


// Author Daniel Lohner (Daniel.Lohner@cern.ch)


ClassImp(AliConversionAODBGHandlerRP)

AliConversionAODBGHandlerRP::AliConversionAODBGHandlerRP() :
    fNEvents(20),
    fBGEventCounter(NULL),
    fNBGEvents(NULL),
    fNBinsZ(7),
    fNBinsCentrality(6),
    fNBinsRP(6),
    fBinLimitsArrayZ(NULL),
    fBinLimitsArrayCentrality(NULL),
    fBinLimitsArrayRP(NULL),
    fBGEvents(fNBinsRP,AliGammaConversionVertexPositionVector(fNBinsZ,AliGammaConversionCentralityVector(fNBinsCentrality,AliGammaConversionBGEventVector(fNEvents))))
{
    // Default Constructor

    // CentralityBins
    fBinLimitsArrayCentrality= new Double_t[fNBinsCentrality+1];

    fBinLimitsArrayCentrality[0] = 0;
    fBinLimitsArrayCentrality[1] = 5;
    fBinLimitsArrayCentrality[2] = 10;
    fBinLimitsArrayCentrality[3] = 20;
    fBinLimitsArrayCentrality[4] = 40;
    fBinLimitsArrayCentrality[5] = 60;
    fBinLimitsArrayCentrality[6] = 90;

    // Vertex Z Binning

    fBinLimitsArrayZ = new Double_t[fNBinsZ+1];

    fBinLimitsArrayZ[0] = -10;
    fBinLimitsArrayZ[1] = -5.4;
    fBinLimitsArrayZ[2] = -2.8;
    fBinLimitsArrayZ[3] = -0.6;
    fBinLimitsArrayZ[4] = 1.4;
    fBinLimitsArrayZ[5] = 3.5;
    fBinLimitsArrayZ[6] = 6.1;
    fBinLimitsArrayZ[7] = 10;

    // Initialize

    Initialize();

}

AliConversionAODBGHandlerRP::AliConversionAODBGHandlerRP(Int_t NEvents,Int_t NBinsZ,Double_t *BinLimitsZ,Int_t NBinsCentrality,Double_t *BinLimitsCentrality,Int_t NBinsRP) :
fNEvents(NEvents),
fBGEventCounter(NULL),
fNBGEvents(NULL),
fNBinsZ(NBinsZ),
fNBinsCentrality(NBinsCentrality),
fNBinsRP(NBinsRP),
fBinLimitsArrayZ(BinLimitsZ),
fBinLimitsArrayCentrality(BinLimitsCentrality),
fBinLimitsArrayRP(NULL),
fBGEvents(fNBinsRP,AliGammaConversionVertexPositionVector(fNBinsZ,AliGammaConversionCentralityVector(fNBinsCentrality,AliGammaConversionBGEventVector(fNEvents))))
{
    Initialize();
}


AliConversionAODBGHandlerRP::~AliConversionAODBGHandlerRP()
{
    if(fBinLimitsArrayZ){
    delete[] fBinLimitsArrayZ;
    fBinLimitsArrayZ=0x0;}
     if(fBinLimitsArrayCentrality){
    delete[] fBinLimitsArrayCentrality;
    fBinLimitsArrayCentrality=0x0;}
     if(fBinLimitsArrayRP){
    delete[] fBinLimitsArrayRP;
    fBinLimitsArrayRP=0x0;}

    if(fBGEventCounter){

	for(Int_t psi=0;psi<fNBinsRP;psi++){
	    for(Int_t z=0;z<fNBinsZ;z++){
		delete[] fBGEventCounter[psi][z];
	    }
	    delete[] fBGEventCounter[psi];
	}
	delete[] fBGEventCounter;
	fBGEventCounter = NULL;
    }

    // Delete pool

    for(Int_t psi=0;psi<fNBinsRP;psi++){
	for(Int_t z=0;z<fNBinsZ;z++){
	    for(Int_t m=0;m<fNBinsCentrality;m++){
		for(Int_t eventCounter=0;eventCounter<fNBGEvents[psi][z][m]&&eventCounter<fNEvents;eventCounter++){

		    for(UInt_t d=0;d<fBGEvents[psi][z][m][eventCounter].size();d++){
		    delete (AliAODConversionPhoton*)(fBGEvents[psi][z][m][eventCounter][d]);
		    }
		}
	    }
	}
    }

    if(fNBGEvents){
	for(Int_t psi=0;psi<fNBinsRP;psi++){
	    for(Int_t z=0;z<fNBinsZ;z++){
		delete[] fNBGEvents[psi][z];
	    }
        delete[] fNBGEvents[psi];
	}
	delete[] fNBGEvents;
	fNBGEvents = NULL;
    }


}


void AliConversionAODBGHandlerRP::Initialize(){

    // RP Binning

    fBinLimitsArrayRP=new Double_t[fNBinsRP+1];

    for(Int_t i=0;i<fNBinsRP+1;i++){
	fBinLimitsArrayRP[i] = i*TMath::Pi()/fNBinsRP;
    }

    // Counter

    if(fBGEventCounter == NULL){
	fBGEventCounter= new Int_t**[fNBinsRP];
    }
    for(Int_t psi=0;psi<fNBinsRP;psi++){
	fBGEventCounter[psi]= new Int_t*[fNBinsZ];

	for(Int_t z=0;z<fNBinsZ;z++){
	    fBGEventCounter[psi][z]=new Int_t[fNBinsCentrality];
	}
    }
    for(Int_t psi=0;psi<fNBinsRP;psi++){
	for(Int_t z=0;z<fNBinsZ;z++){
	    for(Int_t m=0;m<fNBinsCentrality;m++){
		fBGEventCounter[psi][z][m]=0;
	    }
	}
    }


    if(fNBGEvents == NULL){
	fNBGEvents= new Int_t**[fNBinsRP];
    }
    for(Int_t psi=0;psi<fNBinsRP;psi++){
	fNBGEvents[psi]= new Int_t*[fNBinsZ];

	for(Int_t z=0;z<fNBinsZ;z++){
	    fNBGEvents[psi][z]=new Int_t[fNBinsCentrality];
	}
    }
    for(Int_t psi=0;psi<fNBinsRP;psi++){
	for(Int_t z=0;z<fNBinsZ;z++){
	    for(Int_t m=0;m<fNBinsCentrality;m++){
		fNBGEvents[psi][z][m]=0;
	    }
	}
    }

}

//-------------------------------------------------------------

Int_t AliConversionAODBGHandlerRP::GetZBinIndex(Double_t zvalue) const{
  // see header file for documantation
  if(fNBinsZ<2 || zvalue<=fBinLimitsArrayZ[0]){
    return 0;
  }

  for(Int_t i=0; i<fNBinsZ ;i++){
    if(zvalue >= fBinLimitsArrayZ[i] && zvalue <= fBinLimitsArrayZ[i+1]){
      return i;
    }
  }
  return -1;
}
//-------------------------------------------------------------
Int_t AliConversionAODBGHandlerRP::GetCentralityBinIndex(Int_t multiplicity) const{
  if(fNBinsCentrality<2){
    return 0;
  }

  for(Int_t i=0; i<fNBinsCentrality ;i++){
    if(multiplicity >= fBinLimitsArrayCentrality[i] && multiplicity < fBinLimitsArrayCentrality[i+1]){
      return i;
    }
  }
  return -1;
}

//-------------------------------------------------------------

Int_t AliConversionAODBGHandlerRP::GetRPBinIndex(Double_t psivalue) const{
 
  if(fNBinsRP<2 || psivalue<=fBinLimitsArrayRP[0]){
    return 0;
  }

  for(Int_t i=0; i<fNBinsRP ;i++){
    if(psivalue >= fBinLimitsArrayRP[i] && psivalue <= fBinLimitsArrayRP[i+1]){
      return i;
    }
  }
  return -1;
}


//-------------------------------------------------------------
void AliConversionAODBGHandlerRP::AddEvent(TClonesArray * const eventGammas,Double_t psivalue,Double_t zvalue, Int_t multiplicity){


  Int_t z = GetZBinIndex(zvalue);
  Int_t m = GetCentralityBinIndex(multiplicity);
  Int_t psi=GetRPBinIndex(psivalue);
 
  // If Event Stack is full, replace the first entry (First in first out)
  if(fBGEventCounter[psi][z][m] >= fNEvents){
     fBGEventCounter[psi][z][m]=0;
  }

  // Update number of Events stored
  if(fNBGEvents[psi][z][m] < fNEvents){
  fNBGEvents[psi][z][m]++;
  }


  Int_t eventCounter=fBGEventCounter[psi][z][m];

  //clear the vector for old gammas
  for(UInt_t d=0;d<fBGEvents[psi][z][m][eventCounter].size();d++){
   delete (AliAODConversionPhoton*)(fBGEvents[psi][z][m][eventCounter][d]);
  }

  fBGEvents[psi][z][m][eventCounter].clear();


  // add the gammas to the vector
  
    for(Int_t i=0; i< eventGammas->GetEntriesFast();i++){
       fBGEvents[psi][z][m][eventCounter].push_back(new AliAODConversionPhoton(*(AliAODConversionPhoton*)(eventGammas->At(i))));
    }

  fBGEventCounter[psi][z][m]++;


}

//-------------------------------------------------------------
AliGammaConversionPhotonVector* AliConversionAODBGHandlerRP::GetBGGoodGammas(Int_t psi,Int_t zbin, Int_t mbin, Int_t event){
  //see headerfile for documentation      
  return &(fBGEvents[psi][zbin][mbin][event]);
}

