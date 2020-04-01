#ifndef ALIANALYSISEMCALJETHELPEREA_H
#define ALIANALYSISEMCALJETHELPEREA_H

#include <TObject.h>
#include <TString.h>
#include <TArrayD.h>
#include <TArrayI.h>



// ANALYSIS OF EVENT ACTIVITY WITH HIGH PT HADRON BIAS
// Author Filip Krizek   (8 AUG 2019)
namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisEmcalJetHelperEA : public TObject {
  
   public:

   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisEmcalJetHelperEA();
   virtual  ~AliAnalysisEmcalJetHelperEA(){;}

   
  // ######### SETTERS/GETTERS
  Double_t  GetV0A(Int_t runnumber) const; 
  Double_t  GetV0C(Int_t runnumber) const; 
  Double_t  GetV0M(Int_t runnumber) const; 
  Int_t     GetNRuns() const { return  fnRun;}
  Int_t     GetRun(Int_t i) const { return  fRuns[i];}

  Double_t  GetV0APartLevel() const { return fMeanV0A_PartLevel;} 
  Double_t  GetV0CPartLevel() const { return fMeanV0C_PartLevel;} 
  Double_t  GetV0MPartLevel() const { return fMeanV0M_PartLevel;} 

  Double_t  GetV0ADetLevel() const { return fMeanV0A_DetLevel;} 
  Double_t  GetV0CDetLevel() const { return fMeanV0C_DetLevel;} 
  Double_t  GetV0MDetLevel() const { return fMeanV0M_DetLevel;} 

  void      SetV0MeanForMCWithDeltaElectronBug();

 private:


   TArrayD  fMeanV0A;                              // mean V0A signal in incl. MB  run by run 
   TArrayD  fMeanV0C;                              // mean V0C signal in incl. MB  run by run
   TArrayD  fMeanV0M;                              // mean V0C signal in incl. MB  run by run
   TArrayI  fRuns;                                 // run numbers
   Int_t    fnRun;                                 //  the number of runs 
                     
   Double_t fMeanV0A_PartLevel;                    // mean V0A signal in incl. MB particle level 
   Double_t fMeanV0C_PartLevel;                    // mean V0C signal in incl. MB particle level 
   Double_t fMeanV0M_PartLevel;                    // mean V0M signal in incl. MB particle level 
 
   Double_t fMeanV0A_DetLevel;                     // mean V0A signal in incl. MB detector level 
   Double_t fMeanV0C_DetLevel;                     // mean V0C signal in incl. MB detector level 
   Double_t fMeanV0M_DetLevel;                     // mean V0M signal in incl. MB detector level 
 

   AliAnalysisEmcalJetHelperEA(const AliAnalysisEmcalJetHelperEA&);
   AliAnalysisEmcalJetHelperEA& operator=(const AliAnalysisEmcalJetHelperEA&);

   /// \cond CLASSIMP
   ClassDef(AliAnalysisEmcalJetHelperEA, 3); 
   /// \endcond

};
}
}
#endif
