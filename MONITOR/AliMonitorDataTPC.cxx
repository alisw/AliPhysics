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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Data structure for the TPC monitor tree branch                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliLog.h"
#include "AliMonitorDataTPC.h"


ClassImp(AliMonitorDataTPC) 


//_____________________________________________________________________________
AliMonitorDataTPC::AliMonitorDataTPC()
{
// default constructor

  fPt = fEta = fPhi = NULL;
  fSize = 0;
}

//_____________________________________________________________________________
AliMonitorDataTPC::AliMonitorDataTPC(const AliMonitorDataTPC& data) :
  TObject(data)
{
  AliFatal("copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorDataTPC& AliMonitorDataTPC::operator = (const AliMonitorDataTPC& 
						  /*data*/)
{
  AliFatal("assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliMonitorDataTPC::AliMonitorDataTPC(Int_t size)
{
// constructor with given size

  fPt = new Float_t[size];
  fEta = new Float_t[size];
  fPhi = new Float_t[size];
  fSize = size;
}

//_____________________________________________________________________________
AliMonitorDataTPC::~AliMonitorDataTPC()
{
// destructor: free allocated memory

  delete[] fPt;
  delete[] fEta;
  delete[] fPhi;
}

//_____________________________________________________________________________
void AliMonitorDataTPC::SetSize(Int_t size)
{
// set a new array size and allocate memory if necessary

  if (size > fSize) {
    delete[] fPt;
    delete[] fEta;
    delete[] fPhi;
    fPt = new Float_t[size];
    fEta = new Float_t[size];
    fPhi = new Float_t[size];
    fSize = size;
  }
}

//_____________________________________________________________________________
void AliMonitorDataTPC::SetNTracks(Int_t nTracks)
{
// set the number of tracks

  SetSize(nTracks);
  fNTracks = nTracks;
}

//_____________________________________________________________________________
void AliMonitorDataTPC::SetData(Int_t i, Float_t pt, Float_t eta, Float_t phi)
{
// set the data of the i-th track

  if ((i < 0) || (i >= fSize)) {
    AliError(Form("index %d out of range (0-%d)", i, fSize));
    return;
  }
  fPt[i] = pt;
  fEta[i] = eta;
  fPhi[i] = phi;
}
