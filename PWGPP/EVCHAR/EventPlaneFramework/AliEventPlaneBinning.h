/*
***********************************************************
  Binning class used for THnF creation
  Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  *********************************************************
*/

#ifndef ALIEVENTPLANEBINNING_H
#define ALIEVENTPLANEBINNING_H

//#include <TClonesArray.h>
//#include <TBits.h>
//#include <TMath.h>
//#include <TList.h>
//#include <iostream>
//#include <Rtypes.h>
//#include <TArrayS.h>
//#include <THn.h>
//#include <TVector.h>
#include <TObject.h>
#include <TArrayI.h>
#include <TArrayD.h>
#include <TAxis.h>
//#include <TProfile.h>

//const Int_t fgkEPMaxHarmonics = 6;
//const Int_t fgkEPMaxDetectors = 20;


//_________________________________________________________________________
class AliEventPlaneBinning : public TObject {

 public:
  //AliEventPlaneBinning();
  AliEventPlaneBinning(Int_t dim=1);
  AliEventPlaneBinning(AliEventPlaneBinning* EPbinning);
  ~AliEventPlaneBinning();

  // setters
  void SetDim(Int_t dim)  {fDim=dim;}
  void SetVar(Int_t dim, Int_t* var) {fVar=TArrayI(dim, var);}
  void SetNbins(Int_t dim, Int_t* nbins) {fNbins=TArrayI(dim, nbins);}
  void AddBinAxis(Int_t axis, Int_t var, Int_t nwidths, Int_t * nbins, Double_t * edges);
  void SetNchannels(Int_t nchan);

  // getters
  Int_t Dim()  const {return fDim;}
  Int_t* Var()  {return fVar.GetArray();}
  Int_t Var(Int_t var)  {return fVar.GetArray()[var];}
  Int_t* Nbins()  {return fNbins.GetArray();}
  //TArrayD* Array()  {return fArray;}
  TAxis* Axes() {return fAxes;}
  TAxis Axis(Int_t ax) {return fAxes[ax];}

 private:
  Int_t  fDim;
  TArrayI fVar;
  TArrayI fNbins;
  //TArrayD* fArray;
  TAxis fAxes[10];


  AliEventPlaneBinning(const AliEventPlaneBinning &c);
  AliEventPlaneBinning& operator= (const AliEventPlaneBinning &c);

  ClassDef(AliEventPlaneBinning, 1);
};

#endif
