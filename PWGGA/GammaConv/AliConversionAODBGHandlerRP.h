#ifndef __ALICONVERSIONAODBGHANDLERRP_H__
#define __ALICONVERSIONAODBGHANDLERRP_H__

#include "AliLog.h"
#include "TObject.h"
#include "AliAODConversionPhoton.h"
#include "TObjArray.h"
#include "TList.h"
#include <vector>

using namespace std;

// typedef vector<AliAODConversionPhoton*> AliGammaConversionPhotonVector;   				// Vector containing photons
// typedef vector<AliGammaConversionPhotonVector> AliGammaConversionBGEventVector;			// Event contains vector of gammas (AliConversionPhotons)
// typedef vector<AliGammaConversionBGEventVector> AliGammaConversionMultiplicityVector;	// Multiplicity classes containing event vectors
// typedef vector<AliGammaConversionMultiplicityVector> AliGammaConversionBGVector;		// z vertex position ...
typedef vector<AliAODConversionPhoton*> AliGammaConversionPhotonVector;                 // Vector containing photons
typedef vector<AliGammaConversionPhotonVector> AliGammaConversionBGEventVector;         // Event contains vector of gammas (AliConversionPhotons)
typedef vector<AliGammaConversionBGEventVector> AliGammaConversionVertexPositionVector;       // z vertex position ...
typedef vector<AliGammaConversionVertexPositionVector> AliGammaConversionBGVector;       // RP angle



class AliConversionAODBGHandlerRP: public TObject{

  public:

    AliConversionAODBGHandlerRP                     ( Bool_t IsHeavyIon=kFALSE,
                                                      Bool_t UseChargedTrackMult=kTRUE,
                                                      Int_t NEvents=10 );
    
    virtual ~AliConversionAODBGHandlerRP();

    Int_t GetRPBinIndex                             (Double_t psi) const;
    Int_t GetZBinIndex                              ( Double_t z ) const;
    Int_t GetMultiplicityBinIndex                   ( Int_t mult ) const;
    void Initialize                                 ();
    Bool_t FindBins                                 ( TObjArray * const eventGammas, 
                                                      AliVEvent *fInputEvent,
                                                      Int_t &psibin,
                                                      Int_t &zbin);
    Bool_t FindBins                                 ( TList * const eventGammas, 
                                                      AliVEvent *fInputEvent,
                                                      Int_t &psibin,
                                                      Int_t &zbin);

  
    AliGammaConversionPhotonVector* GetBGGoodGammas ( TObjArray * const eventGammas,
                                                      AliVEvent *fInputEvent,
                                                      Int_t event );
    AliGammaConversionPhotonVector* GetBGGoodGammas ( TList * const eventGammas, 
                                                      AliVEvent *fInputEvent,
                                                      Int_t event );
    void AddEvent                                   ( TObjArray * const eventGammas,
                                                      AliVEvent *fInputEvent );
    void AddEvent                                   ( TList * const eventGammas, 
                                                      AliVEvent *fInputEvent );
    Int_t GetNBGEvents                              ()const                                         { return fNEvents                             ;} // Size of the Pool (20)
    Int_t GetNBGEvents                              ( TObjArray * const eventGammas,
                                                      AliVEvent *fInputEvent );
    Int_t GetNBGEvents                              ( TList * const eventGammas,
                                                      AliVEvent *fInputEvent );
    Int_t GetNRPBins                                ()const                                         { return fNBinsRP                             ;}
    Int_t GetNZBins                                 ()const                                         { return fNBinsZ                              ;}
    Int_t GetNMultiplicityBins                      ()const                                         { return fNBinsMultiplicity                   ;}

  private:
    Bool_t                      fIsHeavyIon;                      // flag for heavy ion
    Bool_t                      fUseChargedTrackMult;             // flag for multiplicity switch
    Int_t                       fNEvents;                         // number of events
    Int_t**                     fBGEventCounter;                  //! bg counter
    Int_t**                     fNBGEvents;                       //! n BG events
    Int_t                       fNBinsRP;                         //n RP bins
    Int_t                       fNBinsZ;                          //n z bins
    Int_t                       fNBinsMultiplicity;               //n bins multiplicity
    Double_t*                   fBinLimitsArrayRP;                //! bin limits RP array    
    Double_t*                   fBinLimitsArrayZ;                 //! bin limits z array
    Double_t*                   fBinLimitsArrayMultiplicity;      //! bin limit multiplicity array
    AliGammaConversionBGVector  fBGEvents;                        //background events
//     AliGammaConversionBGVector  fBGPool;                          //background events

    AliConversionAODBGHandlerRP(AliConversionAODBGHandlerRP &original);
    AliConversionAODBGHandlerRP &operator=(const AliConversionAODBGHandlerRP &ref);

  ClassDef(AliConversionAODBGHandlerRP,1);

};
#endif
