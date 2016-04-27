/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2015                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
 
 
  

 
#include "AliQnCorrectionsConstants.h"
#include "AliQnCorrectionsAxes.h"

#include <TObject.h>
#include <TMath.h>
//#include <TList.h>
//#include <TClonesArray.h>
//#include <TRandom3.h>
#include <TAxis.h>
#include <iostream>

ClassImp(AliQnCorrectionsAxes)


//____________________________________________________________________________
AliQnCorrectionsAxes::AliQnCorrectionsAxes(Int_t dim) :
  // Constructor
  TObject(),
  fDim(dim)
{
  for(Int_t i=0; i<AliQnCorrectionsConstants::nHistogramDimensions; i++) {fAxes[i]=TAxis();fVar[i]=-1;fAxesLabels[i]="";}

}

//____________________________________________________________________________
AliQnCorrectionsAxes::AliQnCorrectionsAxes(const AliQnCorrectionsAxes &c) :
  // Copy constructor
  TObject()
{
  fDim = c.Dim();
  for(Int_t i=0; i<fDim; i++) {
    fVar[i]=c.Var(i);
    fAxes[i]=TAxis(c.Axis(i));
    fAxesLabels[i]=c.AxisLabel(i);
  }
  
}




//_______________________________________________________________________________
AliQnCorrectionsAxes::~AliQnCorrectionsAxes()
{
  //
  // De-Constructor
  //
}





//____________________________________________________________________________________
void AliQnCorrectionsAxes::SetAxis(Int_t axis, Int_t var, TAxis ax, TString label) {

  
  fAxesLabels[axis]=label;
  fVar[axis]=var;
  //fNbins.SetAt(ax.GetNbins(),axis);
  fAxes[axis] = TAxis(ax);
  fAxes[axis].SetTitle(label);

  return;
}

//____________________________________________________________________________________
void AliQnCorrectionsAxes::SetAxis(Int_t axis, Int_t var, Int_t nwidths, Int_t * nbins, Double_t * edges, TString label) {


  Int_t Nnewbins=0;
  for(Int_t i=0; i<nwidths; i++) Nnewbins+=nbins[i];
  
  Double_t newbins[Nnewbins+1];
  
  newbins[0] = edges[0];
  
  Int_t ibin = 1;
  
  for(Int_t iw=0; iw<nwidths; iw++) {
  	
  	Double_t xwidth = (edges[iw+1]-edges[iw])/nbins[iw];
  
  	for(Int_t ib=(iw==0 ? 1 : 0); ib<(iw==0 ? nbins[iw]+1 : nbins[iw]); ib++){ 
  		newbins[ibin] = newbins[ibin-1] + xwidth;
  		ibin++;
  	
  }}
  
  fAxesLabels[axis]=label;
  fVar[axis]=var;
  TAxis ax = TAxis(ibin-1, newbins);
  fAxes[axis] = TAxis(ax);
  fAxes[axis].SetTitle(label);

  return;
}




//____________________________________________________________________________________
void AliQnCorrectionsAxes::SetAxis(Int_t axis, Int_t var, Double_t binArray[][2], TString label) {

  fAxesLabels[axis]=label;
  fVar[axis]=var;
  TAxis ax = MakeAxis(binArray);
  fAxes[axis] = TAxis(ax);
  fAxes[axis].SetTitle(label);

  return;
}



//____________________________________________________________________________________
TAxis AliQnCorrectionsAxes::MakeAxis(Double_t binArray[][2]){


  Int_t Nnewbins=0;
  for(Int_t i=1; i<=(binArray[0][1]-1); i++) Nnewbins+=binArray[i][1];
  
  Double_t newbins[Nnewbins+1];
  
  newbins[0] = binArray[0][0];
  
  Int_t ibin = 1;
  
  for(Int_t iw=0; iw<(binArray[0][1]-1); iw++) {
  	
  	Double_t xwidth = (binArray[iw+1][0]-binArray[iw][0])/binArray[iw+1][1];
  
  	for(Int_t ib=(iw==0 ? 1 : 0); ib<(iw==0 ? binArray[iw+1][1]+1 : binArray[iw+1][1]); ib++){ 
  		newbins[ibin] = newbins[ibin-1] + xwidth;
  		ibin++;
  	
  }}
  
  return TAxis(ibin-1, newbins);

}





//____________________________________________________________________________________
void AliQnCorrectionsAxes::SetNchannels(Int_t nchannels){
    
  Double_t channelArray[][2] = {{-0.5, 2}, {-0.5+nchannels, (Double_t)nchannels}};

  SetAxis(fDim-1, 0, channelArray, "Channel number");
  
}





