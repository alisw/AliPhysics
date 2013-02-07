/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//_________________________________________________________________________
//  Class to perform flattening of the event plane distribution
//       It stores necessary parameterizations and apply when requested
//
//*-- Author: Dmitri Peressounko (RRC KI)

// --- ROOT system ---

#include "TMath.h"
#include "TH2.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSEPFlattener.h"



ClassImp(AliPHOSEPFlattener)

//____________________________________________________________________________
AliPHOSEPFlattener::AliPHOSEPFlattener() :
  TNamed(),
  fNCentrBins(0), 
  fNHarmonics(0),
  fNparam(0),
  fParam(0x0)
{
  // default ctor
//  Nothing to initialize
}

//____________________________________________________________________________
AliPHOSEPFlattener::AliPHOSEPFlattener(const char * name) :
  TNamed(name,"PHOSEPFlattener"),
  fNCentrBins(0), 
  fNHarmonics(0),
  fNparam(0),
  fParam(0x0)
{
//  Nothing to initialize
}
//____________________________________________________________________________
AliPHOSEPFlattener::AliPHOSEPFlattener(const AliPHOSEPFlattener & fl):
  TNamed(fl.GetName(),"PHOSEPFlattener"),
  fNCentrBins(0), 
  fNHarmonics(0),
  fNparam(0),
  fParam(0x0)

{
  
  fNCentrBins = fl.fNCentrBins ;
  fNHarmonics = fl.fNHarmonics ;
  fNparam     = fl.fNparam ;
  fParam = new Double32_t[fNparam] ;
  for(Int_t i=0; i<fNparam; i++)
    fParam[i]=fl.fParam[i] ;  
}  

//____________________________________________________________________________
AliPHOSEPFlattener & AliPHOSEPFlattener::operator = (const AliPHOSEPFlattener & fl)
{
  if(this== &fl)
    return *this ;
  SetName(fl.GetName()) ;
  fNCentrBins = fl.fNCentrBins ;
  fNHarmonics = fl.fNHarmonics ;
  fNparam     = fl.fNparam ;
  if(fParam) delete [] fParam ;
  fParam = new Double32_t[fNparam] ;
  for(Int_t i=0; i<fNparam; i++)
    fParam[i]=fl.fParam[i] ;
  
  return *this;
}
//____________________________________________________________________________
AliPHOSEPFlattener::~AliPHOSEPFlattener() 
{
// Nothing to delete
  fNCentrBins = 0 ; 
  fNHarmonics = 0 ;
  fNparam     = 0 ;
  if(fParam)
    delete [] fParam ;
  fParam = 0x0 ;
  
}
//____________________________________________________________________________
Double_t AliPHOSEPFlattener::MakeFlat(Double_t oldPhi,Double_t centrality)const
{
  //Apply flattening using existing parameterizations
  if(fNCentrBins==0) return oldPhi; //No correction encoded
  
  Double_t result = oldPhi ;  
  //Centrality bin
  Int_t icen=(Int_t) (centrality*fNCentrBins/100.) ;
  if(icen>=fNCentrBins)icen = fNCentrBins-1 ;
  
  //Offset in the array
  icen=icen*fNHarmonics ;

  for(Int_t i = 1; i<=fNHarmonics/2; i++){
    Double_t c = 2./i*fParam[icen+2*i-2] ;  //fParam==Mean cos(n*phi) for a given centrality
    Double_t s = 2./i*fParam[icen+2*i-1];   //fParam==Mean sin(n*phi) for a given centrality
    result += c*TMath::Sin(i*oldPhi)-s*TMath::Cos(i*oldPhi) ;      
  }
  return result ;
}
//____________________________________________________________________________
void AliPHOSEPFlattener::SetParameterization(TH2 * h){
 //Fill parameterizations
 //We expect histogram with <cos(i*phi)>, <sin(i*phi)> with centrality bins in x axis
 //and harmonics: <cos(x)>, <sin(x)>, <cos(2x)>, <sin(2x)>, .... in y axis
 
  fNCentrBins = h->GetXaxis()->GetNbins() ;
  fNHarmonics = h->GetYaxis()->GetNbins() ;
  fNparam     = fNCentrBins*fNHarmonics ;
  if(fParam) delete [] fParam ;
  fParam = new Double32_t[fNparam] ;
  for(Int_t i=0; i<fNCentrBins; i++)
    for(Int_t j=0; j<fNHarmonics; j++)
      fParam[i*fNHarmonics+j]=h->GetBinContent(i+1,j+1) ;
}

