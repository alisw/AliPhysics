#ifndef ALIQNCORRECTIONS_AXES_H
#define ALIQNCORRECTIONS_AXES_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2015                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/
 
 

//#include <TClonesArray.h>
//#include <TBits.h>
//#include <TMath.h>
//#include <TList.h>
//#include <iostream>
//#include <Rtypes.h>
//#include <THn.h>
//#include <TVector.h>
#include "AliQnCorrectionsConstants.h"
#include <TObject.h>
#include <TAxis.h>
//#include <TProfile.h>

//const Int_t fgkEPMaxHarmonics = 6;
//const Int_t AliQnCorrectionsConstants::nQnConfigurations = 20;

//const Int_t AliQnCorrectionsConstants::nHistogramDimensions=10;

//_________________________________________________________________________
class AliQnCorrectionsAxes : public TObject {

 public:
  //AliQnCorrectionsAxes();
  AliQnCorrectionsAxes(Int_t dim=1);
  AliQnCorrectionsAxes(const AliQnCorrectionsAxes &c);
  ~AliQnCorrectionsAxes();

  // setters
  void SetDim(Int_t dim)  {fDim=dim;}
  void SetVar(Int_t dim, Int_t* var) {for(Int_t i=0; i<dim; i++) fVar[i]=var[i];}
  void SetAxis(Int_t axis, Int_t var, Double_t binArray[][2], TString label);
  void SetAxis(Int_t axis, Int_t var, Int_t nwidths, Int_t * nbins, Double_t * edges, TString label="");
  void SetAxis(Int_t axis, Int_t var, TAxis ax, TString label="");
  void SetAxisLabel(Int_t axis, TString label) {fAxesLabels[axis]=label;}
  void SetNchannels(Int_t nchan);

  // getters
  Int_t Dim()  const {return fDim;}
  const Int_t* Var() const  {return fVar;}
  Int_t Var(Int_t var) const  {return fVar[var];}
  Int_t  Nbins(Int_t ax) const  {return fAxes[ax].GetNbins();}
  TAxis* Axes() {return fAxes;}
  TAxis Axis(Int_t ax) const {return fAxes[ax];}
  const Double_t* Bins(Int_t ax)const  {return fAxes[ax].GetXbins()->GetArray();}
  TString AxisLabel(Int_t ax) const {return fAxesLabels[ax];}

  Double_t GetLowEdge(Int_t ax) const {return Bins(ax)[0];}
  Double_t GetUpEdge(Int_t ax) const {return Bins(ax)[Nbins(ax)];}
  static TAxis MakeAxis(Double_t binArray[][2]);

 private:


  Int_t  fDim;
  Int_t fVar[AliQnCorrectionsConstants::nHistogramDimensions];
  TAxis fAxes[AliQnCorrectionsConstants::nHistogramDimensions];
  TString fAxesLabels[AliQnCorrectionsConstants::nHistogramDimensions];

  ClassDef(AliQnCorrectionsAxes, 1);
};

#endif
