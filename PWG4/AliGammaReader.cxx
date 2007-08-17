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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for reading data in order to do prompt gamma correlations
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

//---- ANALYSIS system ----
#include "Riostream.h"
#include "AliLog.h"
#include "AliGammaReader.h"

ClassImp(AliGammaReader)


//____________________________________________________________________________
AliGammaReader::AliGammaReader() : 
  TObject(), fDataType(0),
  fCTSEtaCut(0.), fEMCALEtaCut(0.), fPHOSEtaCut(0.),
  fNeutralPtCut(0.),
  fChargedPtCut(0.)
{
  //Ctor

  fPhiEMCALCut[0]=0.;
  fPhiEMCALCut[1]=0.;
  fPhiPHOSCut[0]=0.;
  fPhiPHOSCut[1]=0.;
 
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliGammaReader::AliGammaReader(const AliGammaReader & g) :   
  TObject(g), fDataType(g.fDataType),
  fCTSEtaCut(g.fCTSEtaCut),  fEMCALEtaCut(g.fEMCALEtaCut),  fPHOSEtaCut(g.fPHOSEtaCut),
  fNeutralPtCut(g.fNeutralPtCut),
  fChargedPtCut(g.fChargedPtCut)
{
  // cpy ctor

  fPhiEMCALCut[0]=g.fPhiEMCALCut[0];
  fPhiEMCALCut[1]=g.fPhiEMCALCut[1];
  fPhiPHOSCut[0]=g.fPhiPHOSCut[0];
  fPhiPHOSCut[1]=g.fPhiPHOSCut[1];
}

//_________________________________________________________________________
AliGammaReader & AliGammaReader::operator = (const AliGammaReader & source)
{
  // assignment operator

  if(&source == this) return *this;

  fDataType = source.fDataType ;
  fCTSEtaCut = source.fCTSEtaCut;  
  fEMCALEtaCut = source.fEMCALEtaCut;  
  fPHOSEtaCut = source.fPHOSEtaCut;
  fNeutralPtCut = source.fNeutralPtCut;
  fChargedPtCut = source.fChargedPtCut; 

  fPhiEMCALCut[0]=source.fPhiEMCALCut[0];
  fPhiEMCALCut[1]=source.fPhiEMCALCut[1];
  fPhiPHOSCut[0]=source.fPhiPHOSCut[0];
  fPhiPHOSCut[1]=source.fPhiPHOSCut[1];

  return *this;

}

//_______________________________________________________________
void AliGammaReader::InitParameters()
{
 
  //Initialize the parameters of the analysis.
  fDataType = kData ;
  fCTSEtaCut         = 0.7 ;  
  fEMCALEtaCut         = 0.7 ;  
  fPHOSEtaCut         = 0.12 ;
  fPhiEMCALCut[0] = 80 *TMath::DegToRad();
  fPhiEMCALCut[1] = 190*TMath::DegToRad();
  fPhiPHOSCut[0] = 220. *TMath::DegToRad();
  fPhiPHOSCut[1] = 320.*TMath::DegToRad();
  fNeutralPtCut   = 0.5 ;
  fChargedPtCut   = 0.5 ;

}


//________________________________________________________________
void AliGammaReader::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  printf("Data type           : %d\n", fDataType) ;
  printf("CTS Eta cut           : %f\n", fCTSEtaCut) ;
  printf("EMCAL Eta cut           : %f\n", fEMCALEtaCut) ;
  printf("PHOS Eta cut           : %f\n", fPHOSEtaCut) ;
  printf("Phi EMCAL cut           : [%f, %f]\n", fPhiEMCALCut[0],fPhiEMCALCut[1]) ;
  printf("Phi PHOS cut           : [%f, %f]\n", fPhiPHOSCut[0],fPhiPHOSCut[1]) ;
  printf("pT neutral cut           : %f GeV/c\n", fNeutralPtCut) ;
  printf("pT charged cut           : %f GeV/c\n", fChargedPtCut) ;

} 


