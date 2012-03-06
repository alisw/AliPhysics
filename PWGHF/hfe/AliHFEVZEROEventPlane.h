//
// VZERO event plane task for 2010
//
#ifndef AliHFEVZEROEVENTPLANE_H
#define AliHFEVZEROEVENTPLANE_H

#include "TProfile.h"
#include <AliESDEvent.h>
#include "TString.h"
#include "TH1F.h"
#include "TList.h"


class AliHFEVZEROEventPlane : public TNamed {
 public:
  AliHFEVZEROEventPlane();
  AliHFEVZEROEventPlane(const char *name, const Char_t *title);
  AliHFEVZEROEventPlane(const AliHFEVZEROEventPlane &ref);
  AliHFEVZEROEventPlane& operator=(const AliHFEVZEROEventPlane &ref);
  virtual void Copy(TObject &o) const;
  ~AliHFEVZEROEventPlane();

  void ProcessEvent(AliESDEvent *event);
  
  void  SetNameFile(TString namefile) {fnamefile = namefile;};
  Bool_t OpenInfoCalbration(Int_t run);

  Double_t GetEventPlaneV0A() const {return fEventPlaneV0A;};
  Double_t GetEventPlaneV0C() const {return fEventPlaneV0C;};
  Double_t GetEventPlaneV0()  const {return fEventPlaneV0;};

  TList *GetOutputList() const {return fOutputList;};

 private:
  virtual void Analyze(AliESDEvent* esdEvent); 

  Double_t fEventPlaneV0A;          // Corrected event plane V0A
  Double_t fEventPlaneV0C;          // Corrected event plane V0C
  Double_t fEventPlaneV0;           // Corrected event plane V0

  AliESDEvent* fESD;                //! ESD object
  Int_t        fRun;                // Run number
  TProfile *fMultV0;                //! fMultiplicityV0
  Float_t fV0Cpol,fV0Apol;          // fV0Cpol, fV0Apol
  static const Int_t fgknCentrBin = 9; // Centrality bins
  Float_t fMeanQ[fgknCentrBin][2][2];  // mean for centering
  Float_t fWidthQ[fgknCentrBin][2][2]; // rms for centering
  TString fnamefile;                // name of the file with the coefficient
  TList       *fOutputList;         //! Output list
  TProfile    *fMultV0Before;       //! fMultiplicityV0 Before
  TProfile    *fMultV0After;        //! fMultiplicityV0 After
  TH1F *fQBefore[fgknCentrBin][2][2]; //! Q centering before
  TH1F *fQAfter[fgknCentrBin][2][2];  //! Q centering after
 

  ClassDef(AliHFEVZEROEventPlane, 1);    //Analysis task for high pt analysis 
};

#endif
