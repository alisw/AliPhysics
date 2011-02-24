#ifndef ALIDIELECTRONPID_H
#define ALIDIELECTRONPID_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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

#include <AliPID.h>
#include <AliAODTrack.h>
#include <AliAODPid.h>

#include <AliAnalysisCuts.h>

class TF1;
class TList;
class AliVTrack;
class TGraph;
class AliESDpid;
class AliAODpidUtil;

class AliDielectronPID : public AliAnalysisCuts {
public:
  enum DetType {kITS, kTPC, kTRD, kTOF};
  enum PIDbitTupe {kIgnore=0, kRequire, kIfAvailable};
  
  AliDielectronPID();
  AliDielectronPID(const char*name, const char* title);

  virtual ~AliDielectronPID();

  void AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp=-99999.,
              Double_t pMin=0, Double_t pMax=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire);

  void AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, TF1 * const funUp,
              Double_t pMin=0, Double_t pMax=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire);

  void AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, Double_t nSigmaUp,
              Double_t pMin=0, Double_t pMax=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire);

  void AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, TF1 * const funUp,
              Double_t pMin=0, Double_t pMax=0, Bool_t exclude=kFALSE, UInt_t pidBitType=AliDielectronPID::kRequire);
  
  void SetDefaults(Int_t def);

  //
  //Analysis cuts interface
  //const
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}

  static void SetCorrGraph(TGraph * const gr) { fgFitCorr=gr; }
  static void SetCorrVal(Double_t run);
  static Double_t GetCorrVal()   { return fgCorr; }
  static TGraph *GetCorrGraph()  { return fgFitCorr; }
  
private:
  enum {kNmaxPID=10};
  
  DetType  fDetType[kNmaxPID];    //detector type of nsigma cut
  AliPID::EParticleType fPartType[kNmaxPID]; //particle type
  Float_t  fNsigmaLow[kNmaxPID];  //lower nsigma bound
  Float_t  fNsigmaUp[kNmaxPID];   //upper nsigma bound
  Double_t fPmin[kNmaxPID];       //lower momentum
  Double_t fPmax[kNmaxPID];       //upper momentum
  Bool_t   fExclude[kNmaxPID];    //use as exclusion band
  TF1     *fFunUpperCut[kNmaxPID];//use function as upper cut
  TF1     *fFunLowerCut[kNmaxPID];//use function as lower cut
  UChar_t  fNcuts;                //number of cuts
  UChar_t  fRequirePIDbit[kNmaxPID]; //How to make use of the pid bit (see)

  AliESDpid *fESDpid;             //! esd pid object
  AliAODpidUtil *fAODpidUtil;     //! AOD pid object
  
  static TGraph *fgFitCorr;       //spline fit object to correct the nsigma deviation in the TPC electron band
  static Double_t fgCorr;         //!correction value for current run. Set if fgFitCorr is set and SetCorrVal(run)
                                  // was called
  
  Bool_t IsSelectedITS(AliVTrack * const part, Int_t icut);
  Bool_t IsSelectedTPC(AliVTrack * const part, Int_t icut);
  Bool_t IsSelectedTRD(AliVTrack * const part, Int_t icut);
  Bool_t IsSelectedTOF(AliVTrack * const part, Int_t icut);

  AliDielectronPID(const AliDielectronPID &c);
  AliDielectronPID &operator=(const AliDielectronPID &c);

  ClassDef(AliDielectronPID,3)         // Dielectron PID
};

#endif
