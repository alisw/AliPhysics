#ifndef ALITPCCALIBGAINMULT_H
#define ALITPCCALIBGAINMULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "AliTPCCalPad.h"
#include "TH3F.h"
#include "TF1.h"
#include "THnSparse.h"
#include "THn.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
class TH1F;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliTPCseed;

#include "TTreeStream.h"


class AliTPCcalibGainMult:public AliTPCcalibBase {
public:
  AliTPCcalibGainMult(); 
  AliTPCcalibGainMult(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibGainMult();
  void SetBBParam(TVectorD * param) {fBBParam=param;}
  //  virtual void Terminate();  
  //
  virtual void           Process(AliESDEvent *event);
  virtual void           ProcessV0s(AliESDEvent *event);
  virtual void           ProcessCosmic(const AliESDEvent *event);
  virtual void           ProcessKinks(const AliESDEvent *event);
  virtual void           ProcessTOF(const AliESDEvent *event);  
  virtual void           DumpHPT(const AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze();
  void                   DumpTrack(AliESDtrack * track, AliESDfriendTrack *ftrack, AliTPCseed * seed, Int_t index);
  static Double_t GetTruncatedMeanPosition(Double_t q0, Double_t q1, Double_t q2, Int_t ntracks, Int_t tuneIndex=0, TTreeSRedirector *pcstream=0);
  //
  TH1F   *          GetHistNTracks() const {return fHistNTracks;};
  TH1F   *          GetHistClusterShape() const {return fHistClusterShape;};
  TH3F   *          GetHistQA() const {return fHistQA;};
  //
  THnSparseF *      GetHistGainSector() const {return fHistGainSector;};
  THnSparseF *      GetHistPadEqual() const {return fHistPadEqual;};
  THnSparseF *      GetHistGainMult() const {return fHistGainMult;};  
  THnF       *      GetHistTopology() const {return fHistTopology;};
  //
  THnSparseF * GetHistdEdxMap() const { return fHistdEdxMap;}      // 4D dedx histogram
  THnSparseF * GetHistdEdxMax() const { return fHistdEdxMax;}      // 4D dedx histogram
  THnSparseF * GetHistdEdxTot() const { return fHistdEdxTot;}      // 4D dedx histogram
  TTree *      GetdEdxTree() const {return fdEdxTree;}         // tree for the later minimization

  TGraphErrors* GetGainPerChamber(Int_t padRegion=1, Bool_t plotQA=kFALSE);
  TGraphErrors* GetGainPerChamberRobust(Int_t padRegion=1, Bool_t plotQA=kFALSE, TObjArray *arrQA=0x0, Bool_t normQA=kTRUE);
  //
  const TString&   GetTimeGainID()      const { return fTimeGainID;      }
  const TString&   GetTimeGainStorage() const { return fTimeGainStorage; }
  const TObjArray* GetTimeGainObjects() const { return fTimeGainObjects; }
  //
  void SetMIPvalue(Float_t mip){fMIP = mip;};
  void SetLowerTrunc(Float_t lowerTrunc){fLowerTrunc = lowerTrunc;};
  void SetUpperTrunc(Float_t upperTrunc){fUpperTrunc = upperTrunc;};
  void SetUseMax(Bool_t useMax){fUseMax = useMax;};
  //
  void SetCutMinCrossRows(Int_t crossRows){fCutCrossRows = crossRows;};
  void SetCutMaxEta(Float_t maxEta){fCutEtaWindow = maxEta;};
  void SetCutRequireITSrefit(Bool_t requireItsRefit = kFALSE){fCutRequireITSrefit = requireItsRefit;};
  void SetCutMaxDcaXY(Float_t maxXY){fCutMaxDcaXY = maxXY;};
  void SetCutMaxDcaZ(Float_t maxZ){fCutMaxDcaZ = maxZ;};
  //
  void    SetMinTPCsignalN(Float_t minSignalN) { fMinTPCsignalN=minSignalN; }
  Float_t GetMinTPCsignalN() const             { return fMinTPCsignalN;     }
  //
  void SetMinMomentumMIP(Float_t minMom = 0.4){fMinMomentumMIP = minMom;};
  void SetMaxMomentumMIP(Float_t maxMom = 0.6){fMaxMomentumMIP = maxMom;};
  void SetAlephParameters(Float_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j];};
  //
  //
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
  //
  void     MakeLookup(THnSparse * hist, Char_t * outputFile);
  //
  void     UpdateGainMap();
  void     UpdateClusterParam();

  Double_t GetEntries() const {return fHistGainSector->GetEntries();}

  static void SetMergeEntriesCut(Double_t c) {fgMergeEntriesCut=c;}

private:
  static Double_t fgMergeEntriesCut;  //maximal number of entries for merging  -can be modified via setter

  //
  // parameter specifications
  //
  Float_t fMIP;                  // MIP position to be in fMIP
  Float_t fLowerTrunc;           // lower truncation for dEdx
  Float_t fUpperTrunc;           // upper truncation for dEdx
  //
  Bool_t fUseMax;                 // flag if Qmax or Qtot should be used
  //
  // track cuts
  //
  Int_t   fCutCrossRows;                // minimum number of crossed rows 
  Float_t fCutEtaWindow;                // maximum eta of tracks
  Bool_t  fCutRequireITSrefit;          // if ITSrefit should be required (dangerous in cpass0)
  Float_t fCutMaxDcaXY;                 // max dca_xy (only TPConly resolution is guaranteed!)
  Float_t fCutMaxDcaZ;                  // max dca_z  (dangerous if vDrift is not calibrated)
  Float_t fMinTPCsignalN;               // minimum number of PID clusters
  //
  // definition of MIP window
  //
  Float_t fMinMomentumMIP;              // minimum momentum of MIP region, e.g. 400 MeV
  Float_t fMaxMomentumMIP;              // maximum momentum of MIP region, e.g. 600 MeV
  Float_t fAlephParameters[5];          // parameters for equalization in MIP window, parameter set should be =1 at MIP
  //
  // Store timeGain calibration that was used to create this object
  //
  TString fTimeGainID;            // ID of timeGain object that was used to create this calibration
  TString fTimeGainStorage;       // Storage of time gain object
  TObjArray* fTimeGainObjects;     // Time gain calibration objects used
  //
  // histograms
  //
  TH1F  *fHistNTracks;            //  histogram showing number of ESD tracks per event
  TH1F  *fHistClusterShape;       //  histogram to check the cluster shape
  TH3F  *fHistQA;                 //  dE/dx histogram showing the final spectrum
  //
  //
  THnSparseF * fHistGainSector;   //  histogram which shows MIP peak for each of the 3x36 sectors (pad region)
  THnSparseF * fHistPadEqual;     //  histogram for the equalization of the gain in the different pad regions -> pass0
  THnSparseF * fHistGainMult;     //  histogram which shows decrease of MIP signal as a function
  THnF       * fHistTopology;     //  histogram for topological corrections of signal - dip angle theta and curvature (1/pT)
  TMatrixD *fPIDMatrix;           //! custom PID matrix
  //
  THnSparseF * fHistdEdxMap;      // 4D dedx histogram - per sector/phi
  THnSparseF * fHistdEdxMax;      // 5D dedx histogram - per 1/dedx, tan(theta), tan(phi), mult, pad-type
  THnSparseF * fHistdEdxTot;      // 5D dedx histogram - per 1/dedx, tan(theta), tan(phi), mult, pad-type
  TTree *      fdEdxTree;         // tree for the later minimization
  TVectorD    *fBBParam;          // BetheBloch parameterization used for the dedx expected calculation
  //
  AliTPCcalibGainMult(const AliTPCcalibGainMult&); 
  AliTPCcalibGainMult& operator=(const AliTPCcalibGainMult&); 

  ClassDef(AliTPCcalibGainMult, 6);
};

#endif


