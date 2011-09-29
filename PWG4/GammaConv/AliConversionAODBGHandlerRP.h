#ifndef __ALICONVERSIONAODBGHANDLERRP_H__
#define __ALICONVERSIONAODBGHANDLERRP_H__

#include "AliLog.h"
#include "TObject.h"
#include "AliAODConversionPhoton.h"
#include "TClonesArray.h"
using namespace std;

typedef vector<AliAODConversionPhoton*> AliGammaConversionPhotonVector;

class AliConversionAODBGHandlerRP: public TObject{

public:

    AliConversionAODBGHandlerRP();
    AliConversionAODBGHandlerRP(Int_t NEvents,Int_t NBinsZ,Double_t *BinLimitsZ,Int_t NBinsCentrality,Double_t *BinLimitsCentrality,Int_t NBinsRP);

    
virtual ~AliConversionAODBGHandlerRP();


Int_t GetZBinIndex(Double_t z) const;
Int_t GetCentralityBinIndex(Int_t mult) const;
Int_t GetRPBinIndex(Double_t psi) const;
void Initialize();

typedef vector<AliGammaConversionPhotonVector> AliGammaConversionBGEventVector;       // Event contains vector of gammas (AliConversionPhotons)
typedef vector<AliGammaConversionBGEventVector> AliGammaConversionCentralityVector;  // Centrality classes containing event vectors
typedef vector<AliGammaConversionCentralityVector> AliGammaConversionVertexPositionVector;       // z vertex position ...
typedef vector<AliGammaConversionVertexPositionVector> AliGammaConversionBGVector;       // RP angle


AliGammaConversionPhotonVector* GetBGGoodGammas(Int_t psi,Int_t zbin, Int_t mbin, Int_t event);
void AddEvent(TClonesArray * const eventGammas, Double_t psi,Double_t zvalue, Int_t multiplicity);
Int_t GetNBGEvents()const {return fNEvents;} // Size of the Pool (20)
Int_t GetNBGEvents(Int_t psi,Int_t z,Int_t m)const{if(z<fNBinsZ&&m<fNBinsCentrality){return fNBGEvents[psi][z][m];}else{AliError(Form("Requested BG pool does not exist:  z %i m %i psi %i",z,m,psi));return 0;}};
Int_t GetNZBins()const{return fNBinsZ;};
Int_t GetNRPBins()const{return fNBinsRP;};
Int_t GetNCentralityBins()const{return fNBinsCentrality;};

private:

 Int_t fNEvents; 
 Int_t *** fBGEventCounter; //! bg counter
 Int_t *** fNBGEvents;
 Int_t fNBinsZ; //n z bins
 Int_t fNBinsCentrality; //n bins multiplicity
 Int_t fNBinsRP;
 Double_t *fBinLimitsArrayZ;//! bin limits z array
 Double_t *fBinLimitsArrayCentrality;//! bin limit multiplicity array
 Double_t *fBinLimitsArrayRP;
 AliGammaConversionBGVector fBGEvents; //background events


 ClassDef(AliConversionAODBGHandlerRP,0);

};
#endif
