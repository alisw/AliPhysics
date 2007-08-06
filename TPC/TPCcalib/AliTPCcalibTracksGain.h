#ifndef AliTPCCALIBTRACKSGAIN_H
#define AliTPCCALIBTRACKSGAIN_H

#include <TChain.h>
#include <TNamed.h>


#include <TObjArray.h>
#include <TH2D.h>
#include <TVectorD.h>

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


class AliTPCcalibTracksGain : public TNamed {
public :
   // List of branches
   AliTPCcalibTracksGain(const char* name = 0, const char* title = 0);
   virtual ~AliTPCcalibTracksGain();
   
   static Bool_t   AcceptTrack(AliTPCseed * track);
   
   void            AddTrack(AliTPCseed* seed);
   void            AddCluster(AliTPCclusterMI* cluster);
   Int_t           Evaluate(UInt_t padType, TObjArray* fitParam = 0, TObjArray* fitError = 0, Double_t* redChi2 = 0, Bool_t robust = kFALSE);
   void            GetParameters(TVectorD &fitParam, UInt_t padType, UInt_t fitType);
   void            GetErrors(TVectorD &fitError, UInt_t padType, UInt_t fitType);
   Double_t        GetRedChi2(UInt_t padType, UInt_t fitType);
   AliTPCCalROC*   CreateFitCalROC(UInt_t sector, UInt_t padType, TVectorD &fitParam, Int_t undoTransformation = -1, Bool_t normalizeToPadSize = kFALSE);
   AliTPCCalROC*   CreateCombinedCalROC(const AliTPCCalROC* roc1, const AliTPCCalROC* roc2);
   Double_t        GetPadLength(Double_t lx);
   Int_t           GetPadType(Double_t lx);

private:
   //TTreeSRedirector   *fDebugStream;  //! debug stream for
   //TList          *fOutput;            //output list
   
   TLinearFitter*  fShortFitter;       // simple fitter for short pads
   TLinearFitter*  fMediumFitter;      // simple fitter for medium pads
   TLinearFitter*  fLongFitter;        // simple fitter for long pads
   
   TLinearFitter*  fSqrtShortFitter;   // sqrt fitter for short pads
   TLinearFitter*  fSqrtMediumFitter;  // sqrt fitter for medium pads
   TLinearFitter*  fSqrtLongFitter;    // sqrt fitter for long pads
   
   TLinearFitter*  fLogShortFitter;    // log fitter for short pads
   TLinearFitter*  fLogMediumFitter;   // log fitter for medium pads
   TLinearFitter*  fLogLongFitter;     // log fitter for long pads
   
   UInt_t          fNShortClusters;    // number of clusters registered on short pads
   UInt_t          fNMediumClusters;   // number of clusters registered on medium pads
   UInt_t          fNLongClusters;     // number of clusters registered on medium pads
   AliTPCParamSR*  fTPCparam;          //! helper object for geometry related operations

   static const Double_t fgkM;         // value used in the transformation of the charge values for the logarithmic fitter
       
   ClassDef(AliTPCcalibTracksGain, 1);
};

#endif
