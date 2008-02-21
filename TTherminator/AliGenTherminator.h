#ifndef ALIGENTHERMINATOR_H
#define ALIGENTHERMINATOR_H

// ALICE event generator based on the THERMINATOR model
// It invokes the event generation and then puts the resulting
// event on the stack. All the parameters of the model are 
// accessible via this interface
// Author: Adam.Kisiel@cern.ch

class TF1;
class TClonesArray;

#include "AliGenerator.h"
#include "AliDecayer.h"
#include "AliGenMC.h"
#include "TTherminator.h"

class AliGenTherminator : public AliGenMC
{
 public:

  AliGenTherminator();
  AliGenTherminator(Int_t npart);
  virtual ~AliGenTherminator();
  virtual void Generate();
  virtual void Init();
  // Set model parameters

  void SetFileName(const char *infilename);
  void SetEventNumberInFile(int evnum);

  void SetTemperature(Double_t temp);       
  void SetMiuI(Double_t miu);              
  void SetMiuS(Double_t miu);              
  void SetMiuB(Double_t miu);              
  void SetAlfaRange(Double_t range);         
  void SetRapRange(Double_t range);          

  void SetRhoMax(Double_t rho);  
  void SetTau(Double_t tau);     
  void SetBWA(Double_t bwa);     
  void SetBWVt(Double_t bwv);    
  void SetBWDelay(Double_t bwd); 

  void SetModel(const char *model);
  void SetLhyquidSet(const char *set);
  void SetLhyquidInputDir(const char *inputdir);
  void ReadShareParticleTable();

 protected:
  void     CreateTherminatorInputFile();

  Int_t    fNt;                // CurrentTrack;
  Int_t    fEventNumber;       // Number of the event to read
  TString  fFileName;          // FileName of the file with events

  TString  fFreezeOutModel;    // FreezeOut model to use
  TString  fFOHSlocation;      // Directory where Lhyquid input files are

  // Parameters common for all models
  Double_t fTemperature;       // Freeze-out temperature [GeV/c]
  Double_t fMiuI;              // Isospin chemical potential
  Double_t fMiuS;              // Strance chemical potential
  Double_t fMiuB;              // Baryonic chemical potential
  Double_t fAlfaRange;         // Spatial rapidity range for primordial particles
  Double_t fRapRange;          // Momentum rapidity range for primordial particles
  
  // Model dependent paramters
  Double_t fRhoMax;            // Transverse bound of the system
  Double_t fTau;               // Proper freeze-out time parameter
  Double_t fBWA;               // Blast-wave freeze-out surface slope in t-rho plane
  Double_t fBWVt;              // Blast-wave velocity profile parameter
  Double_t fBWDelay;           // Blast-wave time delay parameter

 private:

  ClassDef(AliGenTherminator,1) // Hijing parametrisation generator
};
#endif

inline void AliGenTherminator::SetTemperature(Double_t temp) { fTemperature = temp; }
inline void AliGenTherminator::SetMiuI(Double_t miu) { fMiuI = miu; }
inline void AliGenTherminator::SetMiuS(Double_t miu) { fMiuS = miu; }
inline void AliGenTherminator::SetMiuB(Double_t miu) { fMiuB = miu; }
inline void AliGenTherminator::SetAlfaRange(Double_t range) { fAlfaRange = range; }
inline void AliGenTherminator::SetRapRange(Double_t range) { fRapRange = range; }

inline void AliGenTherminator::SetRhoMax(Double_t rho) { fRhoMax = rho; }
inline void AliGenTherminator::SetTau(Double_t tau) { fTau = tau; }
inline void AliGenTherminator::SetBWA(Double_t bwa) { fBWA = bwa; }
inline void AliGenTherminator::SetBWVt(Double_t bwv) { fBWVt = bwv; }
inline void AliGenTherminator::SetBWDelay(Double_t bwd) { fBWDelay = bwd; }








