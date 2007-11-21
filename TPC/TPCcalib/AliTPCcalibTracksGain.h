#ifndef AliTPCCALIBTRACKSGAIN_H
#define AliTPCCALIBTRACKSGAIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TChain.h>
#include <TNamed.h>


#include <TObjArray.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include <iostream>
#include <TH1F.h>

using namespace std;

class TTreeSRedirector;
class TH3F; 
class TLinearFitter;

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

class AliTPCcalibTracksGain : public TNamed {
private:
   class PreProcess {
   public:
      PreProcess(AliTPCseed* seed);
      virtual ~PreProcess();
      
      AliTPCseed* GetSeed() { return fSeed; }
      Double_t    GetAngleTrackPadrow(UInt_t segment, UInt_t padType) { return fAngleTrackPadrow[segment + 36*padType]; }
      Double_t    GetAngleTrackBeam(UInt_t segment, UInt_t padType) { return fAngleTrackBeam[segment + 36*padType]; }
      Double_t    GetMeanCharge(UInt_t segment, UInt_t padType) { return fMeanCharge[segment + 36*padType]; }
      Double_t    GetMeanZ(UInt_t segment, UInt_t padType) { return fMeanZ[segment + 36*padType]; }
      Bool_t      IsClusterAccepted(UInt_t clusterIndex) { return fAcceptedClusters[clusterIndex]; }
      
   private:
      AliTPCseed*  fSeed;                   // the seed used for preprocessing
      Double_t*    fAngleTrackPadrow;       // angle between the track and the padrows of one TPC pad region, calculated using a line fit through the appropriate clusters
      Double_t*    fAngleTrackBeam;         // angle between the track and the beam axis, calculated using a line fit through the appropriate clusters
      Double_t*    fMeanCharge;             // the mean charges for this seed (for each pad region)
      Double_t*    fMeanZ;                  // the mean z position for this seed (for each pad region)
      Bool_t       fAcceptedClusters[159];  // contains kFALSE for not accepted clusters, kTRUE for accepted ones (index is padrow)
   };
   
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
   AliTPCcalibTracksGain(const char* name, const char* title, AliTPCcalibTracksCuts* cuts, TNamed* debugStreamPrefix = 0, AliTPCcalibTracksGain* prevIter = 0);
   virtual ~AliTPCcalibTracksGain();
   AliTPCcalibTracksGain& operator=(const AliTPCcalibTracksGain& rhs);
   
   static void testSeed(char* file = "/d/alice07/sma/v4-06-Rev-03/pp/0201/AliESDs.root", Int_t entry = 87, Int_t track = 0); // LOESCH MICH ODER DU WIRST UNTERGEHEN
   static void     AddInfo(TChain* chain, char* debugStreamPrefix = 0, char* prevIterFileName = 0);
   Int_t          AcceptTrack(AliTPCseed* track);
   void            Terminate();
   //Bool_t          AcceptCluster(AliTPCclusterMI* cluster);
   void            Add(AliTPCcalibTracksGain* cal);
   void            AddTrack(AliTPCseed* seed);
   void            AddCluster(AliTPCclusterMI* cluster, AliTPCcalibTracksGain::PreProcess& preProc);
   void            AddCluster2(AliTPCclusterMI* cluster, Float_t momenta, Float_t mdedx, Int_t padType, Float_t xcenter, TVectorD dedxQ, TVectorD dedxM, Float_t fraction, Float_t fraction2, Float_t dedge, TVectorD parY, TVectorD parZ, TVectorD meanPos);
   void            Process(AliTPCseed* seed);
   Long64_t        Merge(TCollection *list);
   void            Evaluate(Bool_t robust = kFALSE, Double_t frac = -1.);
   Int_t           Evaluate(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD* fitParam = 0, TVectorD* fitError = 0, Double_t* redChi2 = 0, Bool_t robust = kFALSE, Double_t frac = -1.);
   void            GetParameters(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitParam);
   void            GetErrors(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitError);
   Double_t        GetRedChi2(UInt_t segment, UInt_t padType, UInt_t fitType);
   void            GetCovarianceMatrix(UInt_t segment, UInt_t padType, UInt_t fitType, TMatrixD& covMatrix);
   AliTPCCalPad*   CreateFitCalPad(UInt_t fitType, Bool_t undoTransformation = kFALSE, Bool_t normalizeToPadSize = kFALSE);
   AliTPCCalROC*   CreateFitCalROC(UInt_t sector, UInt_t fitType, Bool_t undoTransformation = kFALSE, Bool_t normalizeToPadSize = kFALSE);
   AliTPCCalROC*   CreateFitCalROC(UInt_t sector, UInt_t padType, TVectorD &fitParam, UInt_t fitType, Bool_t undoTransformation = kFALSE, Bool_t normalizeToPadSize = kFALSE);
   AliTPCCalROC*   CreateCombinedCalROC(const AliTPCCalROC* roc1, const AliTPCCalROC* roc2);
   TLinearFitter*  GetFitter(UInt_t segment, UInt_t padType, UInt_t fitType);
   static Double_t GetPadLength(Double_t lx);
   static Int_t    GetPadType(Double_t lx);
   void            DumpTrack(AliTPCseed* track);
   Bool_t          GetDedx(AliTPCseed* track, Int_t padType, Int_t* rows);
   
   static Bool_t   GetRowPad(Double_t lx, Double_t ly, Int_t& row, Int_t& pad); // just for debugging, remove me
   UInt_t          fTotalTracks;         // just for debugging, remove me
   UInt_t          fAcceptedTracks;      // just for debugging, remove me
   AliTPCCalPad*   fDebugCalPadRaw;      // just for debugging, remove me
   AliTPCCalPad*   fDebugCalPadCorr;     // just for debugging, remove me

public:
   TTreeSRedirector* fDebugStream;          //! debug stream for debugging
   
   AliTPCFitPad*     fSimpleFitter;         // simple fitter for short pads
   AliTPCFitPad*     fSqrtFitter;           // sqrt fitter for medium pads
   AliTPCFitPad*     fLogFitter;            // log fitter for long pads
   AliTPCFitPad*     fSingleSectorFitter;   // just for debugging

   AliTPCcalibTracksGain* fPrevIter;        // the calibration object in its previous iteration (will not be owned by the new object, don't forget to delete it!)
   
   UInt_t            fNShortClusters[36];   // number of clusters registered on short pads
   UInt_t            fNMediumClusters[36];  // number of clusters registered on medium pads
   UInt_t            fNLongClusters[36];    // number of clusters registered on long pads
   TObjString*       fDebugStreamPrefix;    // pathname of the final location of the debug stream file (may also be an xrootd directory)
   AliTPCcalibTracksCuts* fCuts;            // cuts that are used for sieving the tracks used for calibration

   static       AliTPCParamSR* fgTPCparam;              //! helper object for geometry related operations
   static const Double_t       fgkM;                    // value used in the transformation of the charge values for the logarithmic fitter
   static const char*          fgkDebugStreamFileName;  // filename of the debug stream file
   static const Bool_t         fgkUseTotalCharge;       // whether to use the cluster's total or maximum charge

   ClassDef(AliTPCcalibTracksGain, 1);
};

#endif
