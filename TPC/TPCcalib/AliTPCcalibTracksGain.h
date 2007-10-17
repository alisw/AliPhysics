#ifndef AliTPCCALIBTRACKSGAIN_H
#define AliTPCCALIBTRACKSGAIN_H

#include <TChain.h>
#include <TNamed.h>


#include <TObjArray.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include <iostream>
using namespace std;

class TTreeSRedirector;
class TH3F; 
class TLinearFitter;

class AliTPCClusterParam; 
class AliTPCParamSR; 
class AliTPCCalROC; 
class AliTPCseed; 
class AliTPCclusterMI; 
class AliTrackPointArray;
class TTreeStream;

class AliTPCcalibTracksGain : public TNamed {
public :
   enum {
      kShortPads = 0,
      kMediumPads = 1,
      kLongPads = 2
   };
   enum {
      kSimpleFitter = 0,
      kSqrtFitter = 1,
      kLogFitter = 2
   };
   
   AliTPCcalibTracksGain(const char* name = 0, const char* title = 0);
   virtual ~AliTPCcalibTracksGain();
   static Bool_t   AcceptTrack(AliTPCseed * track);
   void            DumpTrack(AliTPCseed * track);
  Bool_t GetDedx(AliTPCseed * track, Int_t padType, Int_t *rows);
   void            AddTrack(AliTPCseed* seed);
   void            AddCluster(AliTPCclusterMI* cluster);
   Int_t           Evaluate(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD* fitParam = 0, TVectorD* fitError = 0, Double_t* redChi2 = 0, Bool_t robust = kFALSE);
   void            GetParameters(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitParam);
   void            GetErrors(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitError);
   Double_t        GetRedChi2(UInt_t segment, UInt_t padType, UInt_t fitType);
   void            GetCovarianceMatrix(UInt_t segment, UInt_t padType, UInt_t fitType, TMatrixD& covMatrix);
   AliTPCCalROC*   CreateFitCalROC(UInt_t sector, UInt_t padType, TVectorD &fitParam, Int_t undoTransformation = -1, Bool_t normalizeToPadSize = kFALSE);
   AliTPCCalROC*   CreateCombinedCalROC(const AliTPCCalROC* roc1, const AliTPCCalROC* roc2);
   TLinearFitter*  GetFitter(UInt_t segment, UInt_t padType, UInt_t fitType);
   Double_t        GetPadLength(Double_t lx);
   Int_t           GetPadType(Double_t lx);
  //
  //
  static Float_t   TPCBetheBloch(Float_t p, Float_t mass=0.1057);



private:
  TTreeSRedirector   *fDebugStream;       //! debug stream for

   TObjArray*      fShortFitter;          // simple fitter for short pads
   TObjArray*      fMediumFitter;         // simple fitter for medium pads
   TObjArray*      fLongFitter;           // simple fitter for long pads
   
   TObjArray*      fSqrtShortFitter;      // sqrt fitter for short pads
   TObjArray*      fSqrtMediumFitter;     // sqrt fitter for medium pads
   TObjArray*      fSqrtLongFitter;       // sqrt fitter for long pads
   
   TObjArray*      fLogShortFitter;       // log fitter for short pads
   TObjArray*      fLogMediumFitter;      // log fitter for medium pads
   TObjArray*      fLogLongFitter;        // log fitter for long pads
   
   UInt_t          fNShortClusters[36];   // number of clusters registered on short pads
   UInt_t          fNMediumClusters[36];  // number of clusters registered on medium pads
   UInt_t          fNLongClusters[36];    // number of clusters registered on medium pads
   AliTPCParamSR*  fTPCparam;             //! helper object for geometry related operations

   static const Double_t fgkM;            // value used in the transformation of the charge values for the logarithmic fitter
       
   ClassDef(AliTPCcalibTracksGain, 1);
};

#endif
