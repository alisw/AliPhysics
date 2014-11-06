#ifndef AliTPCCALIBTRACKSGAIN_H
#define AliTPCCALIBTRACKSGAIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include <TVectorDfwd.h>
#include <TMatrixDfwd.h>

#include <AliTPCcalibBase.h>

using namespace std;

class TTreeSRedirector;
class TObjString;
class TLinearFitter;
class TProfile; 
class TProfile2D; 
class TH1F; 

class AliTPCClusterParam; 
class AliTPCParamSR; 
class AliTPCCalROC;
class AliTPCCalPad;
class AliTPCseed; 
class AliTPCclusterMI; 
class AliTrackPointArray;
class TTreeStream;
class AliTPCcalibTracksCuts;
class AliTPCFitPad;
class TGraph;
class AliTPCCClusterParam;

class AliTPCcalibTracksGain : public AliTPCcalibBase {
public:
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
   
  AliTPCcalibTracksGain();
  AliTPCcalibTracksGain(const AliTPCcalibTracksGain& obj);
  AliTPCcalibTracksGain(const char* name, const char* title, AliTPCcalibTracksCuts* cuts);
  virtual ~AliTPCcalibTracksGain();
  AliTPCcalibTracksGain& operator=(const AliTPCcalibTracksGain& rhs);
  virtual Long64_t        Merge(TCollection *list);
  virtual void            Process(AliTPCseed* seed);
  virtual void            Terminate();
  virtual void            Analyze();
  //
  // Tracks and cluster manipulation
  //
  Float_t         GetGain(AliTPCclusterMI *cl);
  Float_t         GetMaxNorm(AliTPCclusterMI * cl,  Float_t ky, Float_t kz);
  Float_t         GetQNorm(AliTPCclusterMI * cl,  Float_t ky, Float_t kz);
  void            AddTrack(AliTPCseed* seed);
  void            DumpTrack(AliTPCseed* track);
  Bool_t          GetDedx(AliTPCseed* track, Int_t padType,  Int_t*rows,
			  Int_t &sector, Int_t& npoints, 
			  TVectorD &dedxM, TVectorD &dedxQ, 
			  TVectorD &parY, TVectorD &parZ, TVectorD & meanPos);
  void            AddCluster(AliTPCclusterMI* cluster, Float_t momenta, Float_t mdedx, Int_t padType, Float_t xcenter, TVectorD &dedxQ, TVectorD & dedxM, Float_t fraction, Float_t fraction2, Float_t dedge, TVectorD& parY, TVectorD& parZ, TVectorD& meanPos);
  void            AddTracklet(UInt_t sector, UInt_t padType,TVectorD &dedxQ, TVectorD &dedxM,TVectorD& parY, TVectorD& parZ, TVectorD& meanPos);

  void            Add(AliTPCcalibTracksGain* cal);

  //
  // Get Derived results - gain maps
  //
  AliTPCCalPad*   CreateFitCalPad(UInt_t fitType, Bool_t undoTransformation = kFALSE, Bool_t normalizeToPadSize = kFALSE);
  AliTPCCalROC*   CreateFitCalROC(UInt_t sector, UInt_t fitType, Bool_t undoTransformation = kFALSE, Bool_t normalizeToPadSize = kFALSE);
  AliTPCCalROC*   CreateFitCalROC(UInt_t sector, UInt_t padType, TVectorD &fitParam, UInt_t fitType, Bool_t undoTransformation = kFALSE, Bool_t normalizeToPadSize = kFALSE);
  AliTPCCalROC*   CreateCombinedCalROC(const AliTPCCalROC* roc1, const AliTPCCalROC* roc2);
  //
  void            Evaluate(Bool_t robust = kFALSE, Double_t frac = -1.);
  Bool_t            GetParameters(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitParam);
  void            GetErrors(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitError);
  void            GetCovarianceMatrix(UInt_t segment, UInt_t padType, UInt_t fitType, TMatrixD& covMatrix);
  //
  //
  void            UpdateClusterParam(AliTPCClusterParam *param);
  TGraph *        CreateAmpGraph(Int_t ipad, Bool_t qmax);


  TLinearFitter*  GetFitter(UInt_t segment, UInt_t padType, UInt_t fitType);  
  //void     Process(AliESDEvent *event) {AliTPCcalibBase::Process(event);}
  void     Process(AliVEvent *event) {AliTPCcalibBase::Process(event);}
  //void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliVTrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  
public:
  //
  // Helper function
  //
  static Double_t GetPadLength(Double_t lx);
  static Int_t    GetPadType(Double_t lx);  
  //
  //
  AliTPCcalibTracksCuts* fCuts;            // cuts that are used for sieving the tracks used for calibration
  //
  // kalman fit - alignment of dEdx
  //
  //
  // Fitters
  //
  AliTPCFitPad*     fSimpleFitter;         // simple fitter for short pads
  AliTPCFitPad*     fSqrtFitter;           // sqrt fitter for medium pads
  AliTPCFitPad*     fLogFitter;            // log fitter for long pads
  //
  TLinearFitter*    fFitter0M;          // fitting of the atenuation, angular correction, and mean chamber gain
  TLinearFitter*    fFitter1M;          // fitting of the atenuation, angular correction, and mean chamber gain
  TLinearFitter*    fFitter2M;          // fitting of the atenuation, angular correction, and mean chamber gain
  TLinearFitter*    fFitter0T;          // fitting of the atenuation, angular correction, and mean chamber gain
  TLinearFitter*    fFitter1T;          // fitting of the atenuation, angular correction, and mean chamber gain
  TLinearFitter*    fFitter2T;          // fitting of the atenuation, angular correction, and mean chamber gain
  //
  // angular adn diffusion effect fitter
  // 
  TLinearFitter*    fDFitter0M;          // fitting of the atenuation, angular correction
  TLinearFitter*    fDFitter1M;          // fitting of the atenuation, angular correction
  TLinearFitter*    fDFitter2M;          // fitting of the atenuation, angular correction
  TLinearFitter*    fDFitter0T;          // fitting of the atenuation, angular correction 
  TLinearFitter*    fDFitter1T;          // fitting of the atenuation, angular correction
  TLinearFitter*    fDFitter2T;          // fitting of the atenuation, angular correction
  //
  AliTPCFitPad*     fSingleSectorFitter;   // just for debugging  
  //
  // Conters
  //
  UInt_t          fTotalTracks;         // just for debugging
  UInt_t          fAcceptedTracks;      // just for debugging
  //
  //
  //
  // Setup
  //
   static       AliTPCParamSR* fgTPCparam;              //! helper object for geometry related operations
   static const Double_t       fgkM;                    // value used in the transformation of the charge values for the logarithmic fitter
   static const char*          fgkDebugStreamFileName;  // filename of the debug stream file
   static const Bool_t         fgkUseTotalCharge;       // whether to use the cluster's total or maximum charge

   ClassDef(AliTPCcalibTracksGain, 1);
};

#endif
