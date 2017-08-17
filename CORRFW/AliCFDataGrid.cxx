/* $Id$ */
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
//--------------------------------------------------------------------//
//                                                                    //
// AliCFDataGrid Class                                                //
// Class to handle observed data and correct them                     // 
//                                                                    //
// -- Author : S.Arcelli                                              //
//                                                                    //
// substantially modified by r. vernet                                //
//                                                                    //
//--------------------------------------------------------------------//
//
//
#include "TMath.h"
#include "AliLog.h"
#include "AliCFDataGrid.h"

//____________________________________________________________________
ClassImp(AliCFDataGrid)

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid() : 
  AliCFGridSparse(),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // default constructor
  //
}

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn) :
  AliCFGridSparse(name,title,nVarIn,nBinIn),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // main constructor
  //
  SumW2();// errors saved
}

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const Char_t* name, const Char_t* title, const AliCFContainer &c, Int_t step) :  
  AliCFGridSparse(name,title),
  fSelData(step),
  fContainer(&c)
{
  //
  // main constructor
  // assign directly the selection step
  //

  //simply clones the container's data at specified step
  fData = (THnSparse*) fContainer->GetGrid(fSelData)->GetGrid()->Clone();
  SumW2();
  AliInfo(Form("retrieving measured data from Container %s at selection step %i.",fContainer->GetName(),fSelData));
}

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const AliCFDataGrid& data) : 
  AliCFGridSparse(data),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // copy constructor
  //
  ((AliCFDataGrid &)data).Copy(*this);
}

//____________________________________________________________________
AliCFDataGrid::~AliCFDataGrid()
{
  //
  // destructor
  //
}

//____________________________________________________________________
AliCFDataGrid &AliCFDataGrid::operator=(const AliCFDataGrid &c)
{
  //
  // assigment operator
  //
  if (this != &c) c.Copy(*this);
  return *this;
} 

//____________________________________________________________________
void AliCFDataGrid::ApplyEffCorrection(const AliCFEffGrid &c)
{

  //
  // Apply the efficiency correction
  //
  if (c.GetNVar()!=GetNVar()) {
    AliInfo("Different number of variables, cannot apply correction");
    return;
  }
  Divide(&c);
  AliInfo(Form("correction applied on data grid %s with efficiency %s.",GetName(),c.GetName()));
}
//____________________________________________________________________
void AliCFDataGrid::ApplyBGCorrection(const AliCFDataGrid &c)
{

  //
  // Apply correction for background
  //
  if (c.GetNVar()!=GetNVar()) {
    AliInfo("Different number of variables, cannot apply correction");
    return;
  }
  Add(&c,-1);
  AliInfo(Form("background %s subtracted from data %s.",c.GetName(),GetName()));
}

//____________________________________________________________________
void AliCFDataGrid::Copy(TObject& c) const
{
  // copy function
  AliCFGridSparse::Copy(c);
  AliCFDataGrid& target = (AliCFDataGrid &) c;
  target.fContainer=fContainer;
  target.fSelData=fSelData;
}
