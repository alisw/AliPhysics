/*
***********************************************************
  Binning class used for THnF creation
  Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  *********************************************************
*/


#ifndef ALIEVENTPLANEBINNING_H
#include "AliEventPlaneBinning.h"
#endif

#include <TObject.h>
#include <TMath.h>
#include <TArrayI.h>
#include <TArrayD.h>
//#include <TList.h>
//#include <TClonesArray.h>
//#include <TRandom3.h>
//#include <TArrayS.h>
#include <TAxis.h>
#include <iostream>

ClassImp(AliEventPlaneBinning)


//____________________________________________________________________________
AliEventPlaneBinning::AliEventPlaneBinning(Int_t dim) :
  // Constructor
  TObject(),
  fDim(dim)
  //fArray(0x0)
{
  fVar = TArrayI(dim);
  fNbins = TArrayI(dim);
  //fArray = new TArrayD [dim];
  //fAxes = new TAxis [dim];

  //for(Int_t i=0; i<dim; i++) fAxes[i]=0x0;

}

//____________________________________________________________________________
AliEventPlaneBinning::AliEventPlaneBinning(AliEventPlaneBinning* EPbinning) :
  // Copy constructor
  TObject()
{
  fDim = EPbinning->Dim();
  fVar = TArrayI(fDim);
  fNbins = TArrayI(fDim);
  //fArray = new TArrayD [dim];

  for(Int_t i=0; i<fDim; i++) fAxes[i]=TAxis(EPbinning->Axis(i));

}


//_______________________________________________________________________________
AliEventPlaneBinning::~AliEventPlaneBinning()
{
  //
  // De-Constructor
  //
}





//____________________________________________________________________________________
void AliEventPlaneBinning::AddBinAxis(Int_t axis, Int_t var, Int_t nwidths, Int_t * nbins, Double_t * edges) {

  //if(!fArray)  fArray = new TArrayD [fDim];
  //std::cout<<axis<<"  "<<fArray<<std::endl;

  Int_t Nnewbins=0;
  for(Int_t i=0; i<nwidths; i++) Nnewbins+=nbins[i];
  
  Double_t * newbins;
  newbins = new Double_t [Nnewbins+1];
  
  newbins[0] = edges[0];
  
  Int_t ibin = 1;
  
  for(Int_t iw=0; iw<nwidths; iw++) {
  	
  	Double_t xwidth = (edges[iw+1]-edges[iw])/nbins[iw];
  
  	for(Int_t ib=(iw==0 ? 1 : 0); ib<(iw==0 ? nbins[iw]+1 : nbins[iw]); ib++){ 
  		newbins[ibin] = newbins[ibin-1] + xwidth;
  			ibin++;
  	
  }}
  
  fVar.SetAt(var,axis);
  fNbins.SetAt(ibin-1,axis);
  //fArray[axis] = TArrayD(ibin, newbins);
  fAxes[axis] = TAxis(ibin-1, newbins);
  //std::cout<<axis<<"  "<<fAxes<<std::endl;
  //delete [] newbins;
  return;
}


//____________________________________________________________________________________
void AliEventPlaneBinning::SetNchannels(Int_t nchannels){
    
  const Int_t nChannelwidths = 1;
  Double_t Channeledges[nChannelwidths+1] = {-0.5,nchannels-0.5};
  Int_t Channelbins[nChannelwidths] = {nchannels};

  AddBinAxis(fDim-1, 0, nChannelwidths, Channelbins, Channeledges);
  
}


