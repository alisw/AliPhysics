#ifndef ALITPCTENDERSUPPLY_H
#define ALITPCTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  TPC tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <TString.h>

#include <AliTenderSupply.h>

class TObjArray;
class AliESDpid;
class AliSplineFit;
class AliGRPObject;
class TGraphErrors;
class AliAnalysisManager;
class TF1;

class AliTPCTenderSupply: public AliTenderSupply {
  
public:
  AliTPCTenderSupply();
  AliTPCTenderSupply(const char *name, const AliTender *tender=NULL);
  
  virtual ~AliTPCTenderSupply(){;}

  void SetGainCorrection(Bool_t gainCorr) {fGainCorrection=gainCorr;}
  void SetAttachmentCorrection(Bool_t attCorr) {fAttachmentCorrection=attCorr;}
  void SetDebugLevel(Int_t level)         {fDebugLevel=level;}
  void SetMip(Double_t mip)               {fMip=mip;}
  void SetResponseFunctions(TObjArray *arr) {fArrPidResponseMaster=arr;}
  Double_t GetMultiplicityCorrectionMean(Double_t tpcMulti);
  Double_t GetMultiplicityCorrectionSigma(Double_t tpcMulti);

  void AddSpecificStorage(const char* cdbPath, const char* storage);

  virtual void              Init();
  virtual void              ProcessEvent();
  
private:
  AliESDpid          *fESDpid;         //! ESD pid object
  AliSplineFit       *fGainNew;        //! New gain correction
  AliSplineFit       *fGainOld;        //! Old gain correction
  TGraphErrors       *fGainAttachment; //! new attachment correction

  Bool_t fIsMC;                      //! if or not we have MC analysis
  Bool_t fGainCorrection;            //  Perform gain correction
  Bool_t fAttachmentCorrection;      //  Perform attachment correction
  Bool_t fPcorrection;               //!Perform pressure correction
  Bool_t fMultiCorrection;           //!Perform multiplicity correction
  TObjArray *fArrPidResponseMaster;  //array with gain curves
  TF1 *fMultiCorrMean;               //!multiplicity correction for mean
  TF1 *fMultiCorrSigma;              //!multiplicity correction for resolution
  TObjArray *fSpecificStorages;      //array with specific storages
  
  Int_t fDebugLevel;                 //debug level
  Double_t fMip;                     //mip position
  
  AliGRPObject *fGRP;                //!GRP for pressure temperature correction

  TString fBeamType;                 //! beam type (PP) or (PBPB)
  TString fLHCperiod;                //! LHC period
  TString fMCperiod;                 //! corresponding MC period to use for the splines
  Int_t   fRecoPass;                 //! reconstruction pass

  void SetSplines();
  Double_t GetGainCorrection();

  Double_t GetTPCMultiplicityBin();

  void SetBeamType();
  
  void SetParametrisation();
  void SetRecoInfo();
  
  AliTPCTenderSupply(const AliTPCTenderSupply&c);
  AliTPCTenderSupply& operator= (const AliTPCTenderSupply&c);
  
  ClassDef(AliTPCTenderSupply, 2);  // TPC tender task
};


#endif

