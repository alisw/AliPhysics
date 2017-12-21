#ifndef ALIDIELECTRONPID_H
#define ALIDIELECTRONPID_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronPID                     #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <TBits.h>

#include <AliPID.h>
#include <AliAnalysisCuts.h>
#include <AliTRDPIDResponse.h>

class TF1;
class TList;
class AliVTrack;
class TGraph;
class THnBase;
class AliPIDResponse;
class AliDielectronVarManager;
class AliDielectronVarCuts;

class AliDielectronPID : public AliAnalysisCuts {
public:
  enum DetType {kITS, kTPC, kTRD, kTRD2D, kTRD3D, kTRD7D, kTRDeleEff, kTRDeleEff2D, kTRDeleEff3D, kTRDeleEff7D, kTOF, kEMCAL};
  enum PIDbitType {kIgnore=0, kRequire, kIfAvailable};
  
  AliDielectronPID();
  AliDielectronPID(const char*name, const char* title);

  virtual ~AliDielectronPID();

  void AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp=-99999.,
              Double_t min=0, Double_t max=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire, 
	      Int_t var=-1);

  void AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, TF1 * const funUp,
              Double_t min=0, Double_t max=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire, 
	      Int_t var=-1);

  void AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, Double_t nSigmaUp,
              Double_t min=0, Double_t max=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire, 
	      Int_t var=-1);

  void AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, TF1 * const funUp,
              Double_t min=0, Double_t max=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire, 
	      Int_t var=-1);
  void AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp, Double_t min, Double_t max, Bool_t exclude, UInt_t pidBitType,              TF1 * const funSigma);
  void AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp,
  	      AliDielectronVarCuts *varcuts, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire );

  void AddCut(DetType det, AliPID::EParticleType type, THnBase * const histLow, Double_t nSigmaUp,
	       Double_t min=0, Double_t max=0, Bool_t exclude=kFALSE,
	       UInt_t pidBitType=AliDielectronPID::kRequire, Int_t var=-1);

  void SetDefaults(Int_t def);

  Int_t GetNCuts() { return fNcuts;}
  //
  //Analysis cuts interface
  //const
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}

  static void SetCorrGraph(TGraph * const gr) { fgFitCorr=gr; }
  static TGraph *GetCorrGraph()  { return fgFitCorr; }
  
  static void SetCorrVal(Double_t run);
  static Double_t GetCorrVal()   { return fgCorr; }
  static Double_t GetCorrValdEdx()   { return fgCorrdEdx; }
  
  static void SetCorrGraphdEdx(TGraph * const gr) { fgdEdxRunCorr=gr; }
  static TGraph *GetCorrGraphdEdx()  { return fgdEdxRunCorr; }

  static void SetEtaCorrFunction(TF1 *fun) {fgFunEtaCorr=fun;}
  static TF1* GetEtaCorrFunction() { return fgFunEtaCorr; }
  static void SetCentroidCorrFunction(TH1 *fun) { fgFunCntrdCorr=fun; }
  static void SetWidthCorrFunction(TH1 *fun) { fgFunWdthCorr=fun; }
  static void SetCentroidCorrFunctionITS(TH1 *fun) { fgFunCntrdCorrITS=fun; }
  static void SetWidthCorrFunctionITS(TH1 *fun) { fgFunWdthCorrITS=fun; }
  static void SetCentroidCorrFunctionTOF(TH1 *fun) { fgFunCntrdCorrTOF=fun; }
  static void SetWidthCorrFunctionTOF(TH1 *fun) { fgFunWdthCorrTOF=fun; }

  static Double_t GetEtaCorr(const AliVTrack *track);
  static Double_t GetCntrdCorr(const AliVTrack *track) { return (fgFunCntrdCorr ? GetPIDCorr(track,fgFunCntrdCorr) : 0.0); }
  static Double_t GetWdthCorr(const AliVTrack *track)  { return (fgFunWdthCorr  ? GetPIDCorr(track,fgFunWdthCorr)  : 1.0); }
  static Double_t GetCntrdCorrITS(const AliVTrack *track) { return (fgFunCntrdCorrITS ? GetPIDCorr(track,fgFunCntrdCorrITS) : 0.0); }
  static Double_t GetWdthCorrITS(const AliVTrack *track)  { return (fgFunWdthCorrITS  ? GetPIDCorr(track,fgFunWdthCorrITS)  : 1.0); }
  static Double_t GetCntrdCorrTOF(const AliVTrack *track) { return (fgFunCntrdCorrTOF ? GetPIDCorr(track,fgFunCntrdCorrTOF) : 0.0); }
  static Double_t GetWdthCorrTOF(const AliVTrack *track)  { return (fgFunWdthCorrTOF  ? GetPIDCorr(track,fgFunWdthCorrTOF)  : 1.0); }

private:
  enum {kNmaxPID=30};

  TBits     *fUsedVars;            // list of used variables
  DetType  fDetType[kNmaxPID];    //detector type of nsigma cut
  AliPID::EParticleType fPartType[kNmaxPID]; //particle type
  Float_t  fNsigmaLow[kNmaxPID];  //lower nsigma bound
  Float_t  fNsigmaUp[kNmaxPID];   //upper nsigma bound
  Double_t fmin[kNmaxPID];        //lower cut limit for range
  Double_t fmax[kNmaxPID];        //upper cut limit for range
  Bool_t   fExclude[kNmaxPID];    //use as exclusion band
  TF1     *fFunUpperCut[kNmaxPID];//use function as upper cut
  TF1     *fFunLowerCut[kNmaxPID];//use function as lower cut
  UChar_t  fNcuts;                //number of cuts
  UChar_t  fRequirePIDbit[kNmaxPID]; //How to make use of the pid bit (see)
  UShort_t fActiveCuts[kNmaxPID]; // list of activated cuts
  Double_t fSigmaFunLow[kNmaxPID]; // lower bound for fFunSigma
  Double_t fSigmaFunUp[kNmaxPID];  // upper bound for fFunSigma
  TF1      *fFunSigma[kNmaxPID];   // use function as cut range
  AliDielectronVarCuts *fVarCuts[kNmaxPID]; // varcuts

  AliPIDResponse *fPIDResponse;   //! pid response object
  
  static TGraph *fgFitCorr;       //spline fit object to correct the nsigma deviation in the TPC electron band
  static Double_t fgCorr;         //!correction value for current run. Set if fgFitCorr is set and SetCorrVal(run)
                                  // was called
  static Double_t fgCorrdEdx;     //!dEdx correction value for current run. Set if fgFitCorr is set and SetCorrVal(run)
                                  // was called
  static TF1    *fgFunEtaCorr;    //function for eta correction of electron sigma
  static TH1    *fgFunCntrdCorr;  //function for correction of electron sigma (centroid) in TPC
  static TH1    *fgFunWdthCorr;   //function for correction of electron sigma (width) in TPC
  static TH1    *fgFunCntrdCorrITS;  //function for correction of electron sigma (centroid) in ITS
  static TH1    *fgFunWdthCorrITS;   //function for correction of electron sigma (width) in ITS
  static TH1    *fgFunCntrdCorrTOF;  //function for correction of electron sigma (centroid) in TOF
  static TH1    *fgFunWdthCorrTOF;   //function for correction of electron sigma (width) in TOF
  static TGraph *fgdEdxRunCorr;   //run by run correction for dEdx

  static Double_t GetPIDCorr(const AliVTrack *track, TH1 *hist);
  
  THnBase* fMapElectronCutLow[kNmaxPID];  //map for the electron lower cut in units of n-sigma widths 1 centered to zero
  Bool_t IsSelectedITS(AliVTrack * const part, Int_t icut);
  Bool_t IsSelectedTPC(AliVTrack * const part, Int_t icut, Double_t *values);
	Bool_t IsSelectedTRD(AliVTrack * const part, Int_t icut, AliTRDPIDResponse::ETRDPIDMethod PIDmethod);
  Bool_t IsSelectedTRDeleEff(AliVTrack * const part, Int_t icut, AliTRDPIDResponse::ETRDPIDMethod PIDmethod=AliTRDPIDResponse::kLQ1D);
  Bool_t IsSelectedTOF(AliVTrack * const part, Int_t icut);
  Bool_t IsSelectedEMCAL(AliVTrack * const part, Int_t icut);

  AliDielectronPID(const AliDielectronPID &c);
  AliDielectronPID &operator=(const AliDielectronPID &c);

  ClassDef(AliDielectronPID,8)         // Dielectron PID
};

#endif
