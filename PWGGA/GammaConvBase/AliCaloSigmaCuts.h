#ifndef ALICALOSIGMACUTS_H
#define ALICALOSIGMACUTS_H

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "AliAODMCParticle.h"
#include "AliCaloPhotonCuts.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;


/**
 * @class AliCaloSigmaCuts
 * @brief Class handling all kinds of selection cuts for Sigma analysis
 * @author Steven Merkel
 * @author Adrian Mechler
 * @ingroup GammaConv
 *
 * The cut configuration is set as a string with an 14 digit number.
 * Each digit in the string corresponds to a certain cut type, while
 * its values represent the cut values. The cut configuration is listed here:
 *
 * | Position in the cut string                | Cut type                 |
 * |-------------------------------------------|--------------------------|
 * |                  0                        | PID-Variation            |
 * |                  1                        | FilterBit                |
 * |                  2                        | N Cluster TPC            |
 * |                  3                        | chi2 TPC                 |
 * |                  4                        | N Cluster ITS            |
 * |                  5                        | chi2 ITS                 |
 * |                  6                        | Invert chi2 ITS          |
 * |                  7                        | Min DCA XY               |
 * |                  8                        | Min DCA Z                |
 * |                  9                        | Low p n Sigma TPC Signal |
 * |                  10                       | n Sigma TPC Signal       |
 * |                  11                       | n Sigma TOF Signal       |
 * |                  12                       | Pion Mass lower Cut      |
 * |                  13                       | Pion Mass upper Cut      |
 * |                  14                       | Podolanski Cut           |
 * |                  15                       | Opening Angle Cut        |
 * |                  16                       | Background estimation    |
*/
class AliCaloSigmaCuts : public AliAnalysisCuts {

  public:


    enum cutIds {
      kPIDVariation,
      kFilterBit,
      kNTPCCluster,
      kChi2TPC,
      kNITSCluster,
      kChi2ITS,
      kInvertChi2ITS,
      kMinDCAXY,
      kMinDCAZ,
      kLowPNSigmaTPC,
      kNSigmaTPC,
      kNSigmaTOF,
      kPionMassLower,
      kPionMassUpper,
      kAmenterosCut,
      kOpeningAngleCut,
      kBackgroundEstimation,
      kNCuts
    };

    Bool_t  SetCutIds(TString cutString);
    Int_t   fCuts[kNCuts];
    Bool_t  SetCut(cutIds cutID, Int_t cut);
    Bool_t  UpdateCutString();

    static const char * fgkCutNames[kNCuts];

    Bool_t  InitializeCutsFromCutString(const TString analysisCutSelection);

    AliCaloSigmaCuts(const char *name="SigmaCuts", const char * title="Sigma Cuts");
    AliCaloSigmaCuts(const AliCaloSigmaCuts&);
    AliCaloSigmaCuts& operator=(const AliCaloSigmaCuts&);

    virtual ~AliCaloSigmaCuts();                            //virtual destructor

    virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

    TString GetCutNumber();

     // Cut Selection
    Bool_t PionIsSelectedByMassCut (Double_t    pionMass);
    Bool_t ArmenterosLikeQtCut(Double_t    alpha, Double_t    qT);
    Bool_t SigmaDaughtersOpeningangleCut(Double_t    openingangle);
    Double_t UseRotationmethod(){return fBackgroundestimation;};
    Bool_t TrackIsSelected(AliAODTrack* track, AliPIDResponse* fPIDResponse);
    Bool_t TrackIsSelectedByDCACut(AliAODTrack* track);


    // Set Individual Cuts
    Bool_t SetPIDVariationCut(Int_t PIDVariationCut);
    Bool_t SetFilterBitCut(Int_t FilterBitCut);
    Bool_t SetNClusterTPCCut(Int_t NClusterTPCCut);
    Bool_t SetChi2TPCCut(Int_t Chi2TPCCut);
    Bool_t SetNClusterITSCut(Int_t NClusterITSCut);
    Bool_t SetChi2ITSCut(Int_t Chi2ITSCut);
    Bool_t SetInvertChi2ITSCut(Int_t InvertChi2ITSCut);
    Bool_t SetDCAXYCut(Int_t DCAXYCut);
    Bool_t SetDCAZCut(Int_t DCAZCut);
    Bool_t SetLowPNSigmaTPCCut(Int_t LowPNSigmaTPCCut);
    Bool_t SetNSigmaTPCCut(Int_t NSigmaTPCCut);
    Bool_t SetNSigmaTOFCut(Int_t NSigmaTOFCut);
    Bool_t SetMinPionMassCut(Int_t PionMinMassCut);
    Bool_t SetMaxPionMassCut(Int_t PionMaxMassCut);
    Bool_t SetAmenterosCut(Int_t AmenterosCut);
    Bool_t SetOpeningAngleCut(Int_t OpeningAngleCut);
    Bool_t SetBackgroundEstimation(Int_t BackgroundEstimation);

    // Cut Histogramms
    void    InitCutHistograms(  TString name="");
    TList*  GetCutHistograms() { return fHistograms; }
    void    SetFillCutHistograms( TString name="")  { if(!fHistograms){ InitCutHistograms(name);} return; }
 
    

  protected:
    TList*      fHistograms;                    
    TObjString* fCutString;                     // cut number used for analysis
    TString     fCutStringRead;                 
   

    TF1*        fAmenterosCut;                 


    Int_t       fPIDVariation;                  // min N Cluster ITS
    UInt_t      fFilterBit;                     // FilterBit
    UInt_t      fNClusterTPC;                   // min N Cluster TPC
    Double_t    fChi2TPC;                       // max Chi2 TPC
    Int_t       fNClusterITS;                   // min N Cluster ITS
    Double_t    fChi2ITS;                       // max Chi2 ITS
    Int_t       fInvertChi2ITS;                 // max Chi2 ITS
    Double_t    fDCAXY;                         // min DCA in xy-Richtung
    Double_t    fDCAZ;                          // min DCA in z-Richtung
    Double_t    fLowPNSigmaTPC;                 // max n sigma TPC
    Double_t    fNSigmaTPC;                     // max n sigma TPC
    Double_t    fNSigmaTOF;                     // max N sigma ITS
    Double_t    fMaxPionMass;                   // max pion mass
    Double_t    fMinPionMass;                   // min pion mass
    Double_t    fMaxAlpha;                      // max alpha
    Double_t    fMinAlpha;                      // min alpha
    Double_t    fQt;                            // qT cut
    Double_t    fMaxOpeningAngle;               // max opneningangle
    Double_t    fMinOpeningAngle;               // min openingangle
    Double_t    fBackgroundestimation;          // min openingangle

    //Histogramme
    TH2F*       fHistDEDx;                     
    TH2F*       fHistTOFBeta;                   
    TH2F*       fHistTPCSignal;                 
    
  private:

    /// \cond CLASSIMP
    ClassDef(AliCaloSigmaCuts,5)
    /// \endcond
};


#endif
