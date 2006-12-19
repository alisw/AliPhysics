#ifndef AliEMCALPID_H
#define AliEMCALPID_H

/* $Id$ */
/* History of cvs commits:
 *
 * $Log$
 *
 */

///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALPID
///////////////////////////////////////////////////////////////////////////////

#include "TTask.h"
#include "TArrayD.h"
#include "AliESD.h"
#include "AliPID.h" 

class AliEMCALPID : public TTask {

public:
  
  AliEMCALPID();
  virtual ~AliEMCALPID() { }
  
  void     RunPID(AliESD *esd);
  void     ComputePID(Double_t energy, Double_t lambda0); // give the PID of a cluster
  TArrayD  DistLambda0(Double_t energy, Int_t nature); // compute lambda0 distributions
  
  Double_t GetPID(Int_t idx) const {if (idx>=0&&idx<3) return fPID[idx]; else return 0.;}
  Double_t GetPIDFinal(Int_t idx) const {if (idx>=0&&idx<AliPID::kSPECIESN) return fPIDFinal[idx]; else return 0.;}
  Double_t GetPIDWeight(Int_t idx) const {if (idx>=0&&idx<3) return fPIDWeight[idx]; else return 0.;}
  
  void     SetPID(Double_t val, Int_t idx) {if (idx>=0&&idx<3) fPID[idx] = val;}
  void     SetPIDFinal(Double_t val, Int_t idx) {if (idx>=0&&idx<AliPID::kSPECIESN) fPIDFinal[idx] = val;}
  void     SetPIDWeight(Double_t val, Int_t idx) {if (idx>=0&&idx<3) fPIDWeight[idx] = val;}
  void     SetPrintInfo(Bool_t yesno) {fPrintInfo = yesno;}
   void     SetReconstructor(Bool_t yesno) {fReconstructor = yesno;}
 private:
  
  Double_t Polynomial(Double_t x, Double_t *params);
  
  Bool_t   fPrintInfo;          // flag to decide if details about PID must be printed
  
  Double_t fGamma[6][6];        // Parameter to Compute PID
  Double_t fHadron[6][6];		  // Parameter to Compute PID
  Double_t fPiZero5to10[6][6];  // Parameter to Compute PID
  Double_t fPiZero10to60[6][6]; // Parameter to Compute PID
  
  Float_t fPID[3];
  
  Float_t fPIDFinal[AliPID::kSPECIESN];  // final PID format
  Float_t fPIDWeight[3];                 // order: gamma, pi0, hadrons,
  Double_t fProbGamma;	                // probility to be a Gamma
  Double_t fProbPiZero;	                // probility to be a PiO
  Double_t fProbHadron;	                // probility to be a Hadron
  Bool_t    fReconstructor;               //Fill esdcalocluster when called from EMCALReconstructor
  
  ClassDef(AliEMCALPID, 0)
};

#endif // ALIEMCALPID_H
