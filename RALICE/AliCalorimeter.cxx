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

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliCalorimeter
// Description of a modular calorimeter system.
// A generic 2D geometry is used in which a module is identified by (row,col).
// Obviously this geometry can be a matrix, but also any other regular
// structure is supported, provided the user has adopted a proper convention
// to uniquely address a module via the (row,col) indices.  
// Note : First module is identified as (1,1).
//
// This is the way to define and enter signals into a calorimeter :
//
//   AliCalorimeter cal;
//
//   cal.AddSignal(5,7,85.4);
//   cal.AddSignal(5,7,25.9);
//   cal.AddSignal(3,5,1000);
//   cal.SetSignal(5,7,10.3);
//   cal.Reset(3,5);             // Reset module (3,5) as being 'not fired'
//                               // All module data are re-initialised.
//   cal.SetEdgeOn(1,1);         // Declare module (1,1) as an 'edge module'
//   cal.SetDead(8,3);
//   cal.SetGain(2,8,3.2);
//
//   Float_t vec[3]={6,1,20};
//   cal.SetPosition(2,8,vec,"car");
//
//   AliSignal s;
//   Float_t loc[3]={-1,12,3};
//   s.SetPosition(loc,"car");
//   s.SetSignal(328);
//   cal.AddVetoSignal(s); // Associate (extrapolated) signal as a veto
//
//   cal.Group(2);      // Group 'fired' modules into clusters
//                      // Perform grouping over 2 rings around the center
//   cal.Reset();       // Reset the complete calorimeter
//                      // Normally to prepare for the next event data
//                      // Note : Module gain, offset, edge and dead flags remain
//
//--- Author: Nick van Eijndhoven 13-jun-1997 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliCalorimeter.h"
#include "Riostream.h"

ClassImp(AliCalorimeter) // Class implementation to enable ROOT I/O
 
AliCalorimeter::AliCalorimeter() : AliDevice()
{
// Default constructor, all parameters set to 0.
// Create a calorimeter module matrix with fixed row and column size.
// Note : Due to the dynamic size extension when signals are set,
//        the  "edge modules" can NOT be marked automatically.
//        This has to be done manually by the user via the SetEdgeOn() 
//        memberfunction.
 fNrows=0;
 fNcolumns=0;
 fSwap=0;
 fMatrix=0;
 fClusters=0;
 fHmodules=0;
 fHclusters=0;
 fVetos=0;
 fAttributes=0;
 fPositions=0;
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter::~AliCalorimeter()
{
// Destructor to delete memory allocated to the various arrays and matrices
 if (fClusters)
 {
  delete fClusters;
  fClusters=0;
 }
 if (fVetos)
 {
  delete fVetos;
  fVetos=0;
 }
 if (fHmodules)
 {
  delete fHmodules;
  fHmodules=0;
 }
 if (fHclusters)
 {
  delete fHclusters;
  fHclusters=0;
 }
 if (fMatrix)
 {
  delete fMatrix;
  fMatrix=0;
 }
 if (fPositions)
 {
  delete fPositions;
  fPositions=0;
 }
 if (fAttributes)
 {
  delete fAttributes;
  fAttributes=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter::AliCalorimeter(Int_t nrow,Int_t ncol) : AliDevice()
{
// Create a calorimeter module matrix with fixed row and column size.
// The modules at the edges are automatically marked as "edge modules".
 fNrows=nrow;
 fNcolumns=ncol;
 fClusters=0;

 fSwap=0;
 fMatrix=0;
 fPositions=0;

 fAttributes=new TObjArray(nrow);
 fAttributes->SetOwner();
 
 // Mark the edge modules
 for (Int_t row=1; row<=nrow; row++)
 {
  AliAttribObj* a=new AliAttribObj();
  if (row==1 || row==nrow)
  {
   for (Int_t col=1; col<=ncol; col++)
   {
    a->SetEdgeOn(col);
   }
  }
  else
  {
   a->SetEdgeOn(1);
   a->SetEdgeOn(ncol);
  }
  fAttributes->Add(a);
 }
 
 fHmodules=0;
 fHclusters=0;

 fVetos=0;
}
///////////////////////////////////////////////////////////////////////////
AliCalorimeter::AliCalorimeter(const AliCalorimeter& c) : AliDevice(c)
{
// Copy constructor
 fClusters=0;
 fVetos=0;

 fAttributes=0;

 fHmodules=0;
 fHclusters=0;

 fMatrix=0;
 fPositions=0;

 fNrows=c.fNrows;
 fNcolumns=c.fNcolumns;

 fSwap=c.fSwap;

 if (c.fPositions)
 {
  Int_t nrows=(c.fPositions)->GetMaxRow();
  Int_t ncols=(c.fPositions)->GetMaxColumn();
  for (Int_t irow=1; irow<=nrows; irow++)
  {
   for (Int_t icol=1; icol<=ncols; icol++)
   {
    AliPosition* p=(AliPosition*)(c.fPositions->GetObject(irow,icol));
    if (p) SetPosition(irow,icol,*p);
   } 
  }  
 }

 Int_t size=0;
 if (c.fAttributes) size=c.fAttributes->GetSize();
 if (size)
 {
  fAttributes=new TObjArray(size);
  fAttributes->SetOwner();
  for (Int_t ia=0; ia<size; ia++)
  {
   AliAttribObj* a=(AliAttribObj*)(c.fAttributes->At(ia));
   if (a) fAttributes->AddAt(new AliAttribObj(*a),ia);
  }
 }

 Int_t n=0;
 n=c.GetNclusters();
 if (n)
 {
  fClusters=new TObjArray();
  fClusters->SetOwner();
  for (Int_t icl=1; icl<=n; icl++)
  {
   AliCalcluster* cl=c.GetCluster(icl);
   if (cl) fClusters->Add(new AliCalcluster(*cl));
  }
 }

 n=c.GetNvetos();
 for (Int_t iv=1; iv<=n; iv++)
 {
  AliSignal* s=c.GetVetoSignal(iv);
  if (s) AddVetoSignal(s);
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNrows()
{
// Provide the number of rows for the calorimeter module matrix
 Int_t nrows=fNrows;
 if (!fMatrix) LoadMatrix();
 if (fMatrix && !nrows) nrows=fMatrix->GetMaxRow();
 return nrows;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNcolumns()
{
// Provide the number of columns for the calorimeter module matrix
 Int_t ncols=fNcolumns;
 if (!fMatrix) LoadMatrix();
 if (fMatrix && !ncols) ncols=fMatrix->GetMaxColumn();
 return ncols;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetSignal(Int_t row,Int_t col,Float_t sig)
{
// Set the signal for a certain calorimeter module.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 if (!fMatrix) LoadMatrix();

 if (!fMatrix)
 {
  fMatrix=new AliObjMatrix();
  fMatrix->SetSwapMode(fSwap);
 }

 AliCalmodule* mx=GetModule(row,col);
 if (mx) // Existing module
 {
  mx->SetSignal(sig);
 }
 else // Initialise for a new module
 {
  AliCalmodule m;
  m.SetRow(row);
  m.SetColumn(col);
  m.SetSignal(sig);
  AliPosition* r=0;
  if (fPositions) r=(AliPositionObj*)fPositions->GetObject(row,col);
  if (r) m.SetPosition(*r);
  if (fAttributes)
  {
   AliAttribObj* a=0;
   if (row <= fAttributes->GetSize()) a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (a->GetGainFlag(col)) m.SetGain(a->GetGain(col));
    if (a->GetOffsetFlag(col)) m.SetOffset(a->GetOffset(col));
    if (a->GetDeadValue(col)) m.SetDead();
    if (a->GetEdgeValue(col)) m.SetEdgeValue(a->GetEdgeValue(col));
   }
  }
  AddHit(m);
  fMatrix->EnterObject(row,col,fHits->Last());
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddSignal(Int_t row, Int_t col, Float_t sig)
{
// Add the signal to a certain calorimeter module.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::AddSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }
 
  AliCalmodule* m=GetModule(row,col);
  if (!m) // initialise for new modules
  {
   SetSignal(row,col,sig);
  }
  else
  {
   m->AddSignal(sig);
  }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddSignal(AliCalmodule* mod)
{
// Add the signal of module mod to the current calorimeter data.
// This enables mixing of calorimeter data of various events.
//
// Note : The position and attributes according to the user provided data
//        for the corresponding (row,col) location will be used.
//        In case there is no user provided data present, the position and
//        attributes of the first module added to the corresponding (row,col)
//        location will be taken, except for the "edge" and "dead" indicators.
//        The latter will then both be set to 0.

 if (!mod) return;
 
 Int_t row=mod->GetRow();
 Int_t col=mod->GetColumn();
 Float_t sig=mod->GetSignal();

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::AddSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 if (!fMatrix) LoadMatrix();
 
 if (!fMatrix)
 {
  fMatrix=new AliObjMatrix();
  fMatrix->SetSwapMode(fSwap);
 }

 AliCalmodule* mx=GetModule(row,col);
 if (!mx) // No module existed yet at this position
 {
  AliCalmodule m(*mod);
  AliPosition* r=0;
  if (fPositions) r=(AliPositionObj*)fPositions->GetObject(row,col);
  if (r) m.SetPosition(*r);
  // Don't take the dead and edge attributes from this module,
  // but from the calorimeter dbase, if present.
  m.SetEdgeOff();
  m.SetAlive();
  if (fAttributes)
  {
   AliAttribObj* a=0;
   if (row <= fAttributes->GetSize()) a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (a->GetGainFlag(col)) m.SetGain(a->GetGain(col));
    if (a->GetOffsetFlag(col)) m.SetOffset(a->GetOffset(col));
    if (a->GetDeadValue(col)) m.SetDead();
    if (a->GetEdgeValue(col)) m.SetEdgeValue(a->GetEdgeValue(col));
   }
  }
  AddHit(m);
  fMatrix->EnterObject(row,col,fHits->Last());
 }
 else
 {
  mx->AddSignal(sig);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Reset(Int_t row,Int_t col)
{
// Reset the signal for a certain calorimeter module.
// Note : Module position and attributes remain unchanged.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::Reset* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }
 
 AliCalmodule* m=GetModule(row,col);
 if (m)
 {
  RemoveHit(m);
  fMatrix->RemoveObject(row,col);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Reset(Int_t mode)
{
// Reset the signals for the complete calorimeter.
// Normally this is done to prepare for the data of the next event.
//
// mode = 0 : Swap mode, module positions and attributes remain unchanged.
//        1 : Swap mode, module positions and attributes are cleared.
//
// The default is mode=0.
//
// Note : In the case of reading AliCalorimeter objects from a data file,
//        one has to reset the AliCalorimeter object with mode=1
//        (or explicitly delete it) before reading-in the next object
//        in order to prevent memory leaks.

 if (mode<0 || mode>1)
 {
  cout << " *AliCalorimeter::Reset* Wrong argument. mode = " << mode << endl;
  return;
 }

 AliDevice::Reset(mode);

 if (fClusters)
 {
  delete fClusters;
  fClusters=0;
 }

 if (fVetos)
 {
  delete fVetos;
  fVetos=0;
 }

 if (mode==1)
 {
  if (fMatrix)
  {
   delete fMatrix;
   fMatrix=0;
  }
  if (fPositions)
  {
   delete fPositions;
   fPositions=0;
  }
 }
 else
 {
  if (fMatrix) fMatrix->Reset();
 }

 // Free memory allocated for the various arrays.
 if (mode==1)
 {
  if (fAttributes)
  {
   delete fAttributes;
   fAttributes=0;
  }
  if (fHmodules)
  {
   delete fHmodules;
   fHmodules=0;
  }
  if (fHclusters)
  {
   delete fHclusters;
   fHclusters=0;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetSignal(Int_t row,Int_t col,Int_t mode)
{
// Provide the signal of a certain calorimeter module.
// In case the module was marked dead, 0 is returned.
//
// mode = 0 : Just the module signal is returned
//        1 : The module signal is corrected for the gain and offset.
//            In case the gain value was not set, gain=1 will be assumed.
//            In case the gain value was 0, a signal value of 0 is returned.
//            In case the offset value was not set, offset=0 will be assumed.
//
// The corrected signal (sigc) is determined as follows :
//
//              sigc=(signal/gain)-offset 
//
// The default is mode=0.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
 
 Float_t signal=0;
 Float_t gain=1;
 Float_t offset=0;
 AliCalmodule* m=GetModule(row,col);
 if (m)
 {
  Int_t dead=m->GetDeadValue();
  if (!dead) signal=m->GetSignal();

  if (mode==0 || dead) return signal;

  // Correct the signal for the gain and offset
  if (GetGainFlag(row,col))
  {
   gain=GetGain(row,col);
  }
  else
  {
   if (m->GetGainFlag()) gain=m->GetGain();
  }

  if (GetOffsetFlag(row,col))
  {
   offset=GetOffset(row,col);
  }
  else
  {
   if (m->GetOffsetFlag()) offset=m->GetOffset();
  }

  if (fabs(gain)>0.)
  {
   signal=(signal/gain)-offset;
  }
  else
  {
   signal=0;
  }
 }
 return signal;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetEdgeOn(Int_t row,Int_t col)
{
// Indicate a certain calorimeter module as 'edge module'.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetEdgeOn* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 if (!fAttributes)
 {
  fAttributes=new TObjArray(row);
  fAttributes->SetOwner();
 }
 else
 {
  if (row > fAttributes->GetSize()) fAttributes->Expand(row);
 }

 AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
 if (a)
 {
  a->SetEdgeOn(col);
 }
 else
 {
  a=new AliAttribObj();
  a->SetEdgeOn(col);
  fAttributes->AddAt(a,row-1);
 } 

 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetEdgeOn();
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetEdgeOff(Int_t row,Int_t col)
{
// Indicate a certain calorimeter module as 'non-edge module'.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetEdgeOff* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 // Only action on fAttributes in case an attribute is present at (row,col),
 // since by default a module has edge=0 unless explicitly set otherwise.
 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a) a->SetEdgeOff(col);
  }
 } 

 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetEdgeOff();
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetDead(Int_t row,Int_t col)
{
// Indicate a certain calorimeter module as 'dead module'

 // Check for (row,col) boundaries in case of a fixed size calorimeter
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetDead* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 // Make Attributes storage 1 row (and also 1 column) larger than needed
 // because the 'edge value' of the (future) surrounding modules has
 // to be updated as well.
 if (!fAttributes)
 {
  fAttributes=new TObjArray(row+1);
  fAttributes->SetOwner();
 }
 else
 {
  if (row >= fAttributes->GetSize()) fAttributes->Expand(row+1);
 }

 AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
 if (a)
 {
  a->SetDead(col);
 }
 else
 {
  a=new AliAttribObj();
  a->SetDead(col);
  fAttributes->AddAt(a,row-1);
 } 

 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetDead();
 
 // Increase the 'edge value' of surrounding modules
 Int_t rlow=row-1;
 Int_t rup=row+1;
 Int_t clow=col-1;
 Int_t cup=col+1;
 
 if (rlow < 1) rlow=row;
 if (clow < 1) clow=col;
 
 for (Int_t i=rlow; i<=rup; i++)
 {
  for (Int_t j=clow; j<=cup; j++)
  {
   if (i!=row || j!=col) // No increase of edge value for the 'dead' module itself
   {
    AliAttribObj* a=(AliAttribObj*)fAttributes->At(i-1);
    if (a)
    {
     a->IncreaseEdgeValue(j);
    }
    else
    {
    a=new AliAttribObj();
    a->SetEdgeOn(j);
    fAttributes->AddAt(a,i-1);
    } 

    AliCalmodule* m=GetModule(i,j);
    if (m) m->IncreaseEdgeValue();
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetAlive(Int_t row,Int_t col)
{
// Indicate a certain calorimeter module as 'active module'.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetAlive* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 // Only action on fAttributes in case an attribute is present at (row,col),
 // since by default a module has dead=0 unless explicitly set otherwise.
 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a) a->SetAlive(col);
  }
 } 

 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetAlive();
 
 // Decrease the 'edge value' of surrounding modules
 Int_t rlow=row-1;
 Int_t rup=row+1;
 Int_t clow=col-1;
 Int_t cup=col+1;
 
 if (rlow < 1) rlow=row;
 if (clow < 1) clow=col;

 for (Int_t i=rlow; i<=rup; i++)
 {
  for (Int_t j=clow; j<=cup; j++)
  {
   if (i!=row || j!=col) // No decrease of edge value for the 'alive' module itself
   {
    if (i <= fAttributes->GetSize())
    {
     AliAttribObj* a=(AliAttribObj*)fAttributes->At(i-1);
     if (a) a->DecreaseEdgeValue(j);
    }
    AliCalmodule* m=GetModule(i,j);
    if (m) m->DecreaseEdgeValue();
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetGain(Int_t row,Int_t col,Float_t gain)
{
// Set the gain value for a certain calorimeter module.
// See the memberfunction GetSignal() for a definition of the gain value.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetGain* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 if (!fAttributes)
 {
  fAttributes=new TObjArray(row);
  fAttributes->SetOwner();
 }
 else
 {
  if (row > fAttributes->GetSize()) fAttributes->Expand(row);
 }

 AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
 if (a)
 {
  a->SetGain(gain,col);
 }
 else
 {
  a=new AliAttribObj();
  a->SetGain(gain,col);
  fAttributes->AddAt(a,row-1);
 } 

 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetGain(gain);
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetOffset(Int_t row,Int_t col,Float_t offset)
{
// Set the offset value for a certain calorimeter module.
// See the memberfunction GetSignal() for a definition of the offset value.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetOffset* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 if (!fAttributes)
 {
  fAttributes=new TObjArray(row);
  fAttributes->SetOwner();
 }
 else
 {
  if (row > fAttributes->GetSize()) fAttributes->Expand(row);
 }

 AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
 if (a)
 {
  a->SetOffset(offset,col);
 }
 else
 {
  a=new AliAttribObj();
  a->SetOffset(offset,col);
  fAttributes->AddAt(a,row-1);
 } 

 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetOffset(offset);
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetPosition(Int_t row,Int_t col,Float_t* vec,TString f)
{
// Set the position in user coordinates for a certain calorimeter module
 Ali3Vector r;
 r.SetVector(vec,f);
 SetPosition(row,col,r);
} 
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetPosition(Int_t row,Int_t col,Ali3Vector& r)
{
// Set the position for a certain calorimeter module

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::SetPosition* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return;
 }

 if (!fPositions)
 {
  fPositions=new AliObjMatrix();
  fPositions->SetOwner();
  fPositions->SetSwapMode(fSwap);
 }

 AliPositionObj* p=(AliPositionObj*)fPositions->GetObject(row,col);

 if (p)
 {
  p->Load(r);
 }
 else
 {
  p=new AliPositionObj();
  p->Load(r);
  fPositions->EnterObject(row,col,p);
 }

 // Update the position of the calorimeter module itself as well if it exists
 AliCalmodule* m=GetModule(row,col);
 if (m) m->SetPosition(r);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetEdgeValue(Int_t row,Int_t col)
{
// Provide the value of the edge flag of a certain module.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetEdgeValue* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 Int_t edge=0;

 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (col <= a->GetNcalflags())
    {
     edge=a->GetEdgeValue(col);
     return edge;
    }
   }
  }
 }

 AliCalmodule* m=GetModule(row,col);
 if (m) edge=m->GetEdgeValue();
 return edge;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetDeadValue(Int_t row,Int_t col)
{
// Provide the value of the dead flag of a certain module

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetDeadValue* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 Int_t dead=0;

 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (col <= a->GetNcalflags())
    {
     dead=a->GetDeadValue(col);
     return dead;
    }
   }
  }
 }

 AliCalmodule* m=GetModule(row,col);
 if (m) dead=m->GetDeadValue();
 return dead;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetGainFlag(Int_t row,Int_t col)
{
// Provide the value of the gain flag of a certain module.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetGainFlag* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 Int_t gf=0;

 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (col <= a->GetNcalflags())
    {
     gf=a->GetGainFlag(col);
     return gf;
    }
   }
  }
 }

 AliCalmodule* m=GetModule(row,col);
 if (m) gf=m->GetGainFlag();
 return gf;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetOffsetFlag(Int_t row,Int_t col)
{
// Provide the value of the offset flag of a certain module.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetOffsetFlag* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 Int_t of=0;

 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (col <= a->GetNcalflags())
    {
     of=a->GetOffsetFlag(col);
     return of;
    }
   }
  }
 }

 AliCalmodule* m=GetModule(row,col);
 if (m) of=m->GetOffsetFlag();
 return of;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetGain(Int_t row,Int_t col)
{
// Provide the gain value of a certain module.
// See the memberfunction GetSignal() for a definition of the gain value.
//
// In case the gain value is unknown, the value 0 will be returned. 

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetGain* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 Float_t gain=0;

 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (col <= a->GetNcalflags())
    {
     if (a->GetGainFlag(col))
     {
      gain=a->GetGain(col);
      return gain;
     }
    }
   }
  }
 }

 AliCalmodule* m=GetModule(row,col);
 if (m)
 {
  if (m->GetGainFlag())
  {
   gain=m->GetGain();
  }
 }
 return gain;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetOffset(Int_t row,Int_t col)
{
// Provide the offset value of a certain module.
// See the memberfunction GetSignal() for a definition of the offset value.
//
// In case the offset value is unknown, the value 0 will be returned. 

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetOffset* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 Float_t offset=0;

 if (fAttributes)
 {
  if (row <= fAttributes->GetSize())
  {
   AliAttribObj* a=(AliAttribObj*)fAttributes->At(row-1);
   if (a)
   {
    if (col <= a->GetNcalflags())
    {
     if (a->GetOffsetFlag(col))
     {
      offset=a->GetOffset(col);
      return offset;
     }
    }
   }
  }
 }

 AliCalmodule* m=GetModule(row,col);
 if (m)
 {
  if (m->GetOffsetFlag())
  {
   offset=m->GetOffset();
  }
 }
 return offset;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::GetPosition(Int_t row,Int_t col,Float_t* vec,TString f)
{
// Return the position in user coordinates for a certain calorimeter module
 vec[0]=0;
 vec[1]=0;
 vec[2]=0;

 AliPosition* p=GetPosition(row,col);
 if (p) p->GetVector(vec,f);
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliCalorimeter::GetPosition(Int_t row,Int_t col)
{
// Access to the position of a certain calorimeter module.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetPosition* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }

 AliPositionObj* po=0;
 if (fPositions) po=(AliPositionObj*)fPositions->GetObject(row,col);
 if (po) return po;

 AliCalmodule* m=GetModule(row,col);
 return m;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalorimeter::GetClusteredSignal(Int_t row,Int_t col)
{
// Provide the module signal after clustering.

 // Check for (row,col) boundaries.
 if (row<1 || col<1 || (fNrows && fNcolumns && (row>fNrows || col>fNcolumns)))
 {
  cout << " *AliCalorimeter::GetClusteredSignal* row,col : " << row << "," << col
       << " out of range." << endl;
  if (fNrows && fNcolumns) cout << " Nrows,Ncols = " << fNrows << "," << fNcolumns << endl;
  return 0;
 }
 
 Float_t sig=0;

 AliCalmodule* m=GetModule(row,col);
 if (m) sig=m->GetClusteredSignal();

 return sig;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNsignals() const
{
// Provide the number of modules that contain a signal
// Note : The number of modules marked 'dead' but which had a signal
//        are included.
 return GetNhits();
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Group(Int_t n,Int_t mode)
{
// Group the individual modules into clusters.
// Module signals of n rings around the central module will be grouped.
// The grouping process will start with the module containing the highest signal
// in an iterative way.
// For this all fired modules are ordered w.r.t. decreasing signal.
// The search mode for the module signal hierarchy can be specified by the user.
//
// mode = 1 : Search performed via the (row,col) structure of the matrix (SortM)
//        2 : Search performed via the linear array of fired modules (SortA)
//
// See the docs of the memberfunctions SortM and SortA for additional details.
//
// Default values : n=1 mode=1.

 if (mode<1 || mode>2)
 {
  cout << " *AliCalorimeter::Group* Invalid mode : " << mode << endl;
  cout << " Default value mode=1 will be used." << endl;
  mode=1;
 }

 if (fClusters)
 {
  delete fClusters;
  fClusters=0;
 }

 if (!fMatrix) LoadMatrix();

 if (!fMatrix) return;
 
 Int_t nsignals=GetNsignals();
 if (nsignals > 0) // Only continue if there are fired modules
 {
  if (GetNclusters() > 0) Ungroup(); // Restore unclustered situation if needed
 
  // Order the modules with decreasing signal
  if (mode==1) SortM(); 
  if (mode==2) SortA();

  Int_t nord=0;
  if (fOrdered) nord=fOrdered->GetEntries();
 
  // Clustering of modules. Start with the highest signal.
  fClusters=new TObjArray();
  fClusters->SetOwner();
  Int_t row=0;
  Int_t col=0;
  AliCalcluster* c=0;
  for (Int_t i=0; i<nord; i++)
  {
   AliCalmodule* m=(AliCalmodule*)fOrdered->At(i);
   if (!m) continue;

   row=m->GetRow();    // row number of cluster center
   col=m->GetColumn(); // column number of cluster center

   // only use modules not yet used in a cluster
   if (m->GetClusteredSignal() > 0.)
   {
    Int_t edge=GetEdgeValue(row,col);
    c=new AliCalcluster();
    if (!edge) c->Start(*m);   // module to start the cluster if not on edge
    if (c->GetNmodules() > 0)  // cluster started successfully (no edge)
    {
     fClusters->Add(c);
     AddRing(row,col,n); // add signals of n rings around the center
    }
    else
    {
     if (c) delete c;
     c=0;
    }
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SortM()
{
// Order the modules with decreasing signal by looping over the (row,col) grid
// of the matrix.
// Modules which were declared as "Dead" will be rejected.
// The gain etc... corrected module signals will be used in the ordering process.
//
// Note : This method may become slow for large, very finely granulated calorimeters.
//
// Very specific case :
// ====================
// In case of various overlapping showers of which the central modules have
// EXACTLY the same signal this ordering procedure may have the following
// advantages and disadvantages.
//
// Advantages :
// ------------
// * In case of multi-overlapping showers, the central shower will NOT
//   be "eaten-up" from both sides, resulting in a slightly more accurate
//   cluster signal.
// * This method produces re-producable results, irrespective of the filling
//   order of the matrix modules.
//
// Disadvantages :
// ---------------
// * In case of a very high occupancy, there might be a slight effect on the
//   cluster signals depending on the geometrical location in the detector matrix.

 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }

 Int_t nrows=fMatrix->GetMaxRow();
 Int_t ncols=fMatrix->GetMaxColumn();

 Float_t signal=0.;
 Int_t nord=0;
 for (Int_t irow=1; irow<=nrows; irow++) // loop over all modules of the matrix
 {
  for (Int_t icol=1; icol<=ncols; icol++)
  {
   AliCalmodule* m=(AliCalmodule*)fMatrix->GetObject(irow,icol);
   if (!m) continue;

   signal=m->GetSignal(1,1);   // get the gain etc... corrected signal
   if (signal <= 0.) continue; // only take alive modules with a signal
 
   if (nord == 0) // store the first module with a signal at the first ordered position
   {
    if (!fOrdered)
    {
     Int_t nhits=GetNhits();
     fOrdered=new TObjArray(nhits);
    }
    nord++;
    fOrdered->AddAt(m,nord-1);
    continue;
   }
 
   for (Int_t j=0; j<=nord; j++) // put module in the right ordered position
   {
    if (j == nord) // module has smallest signal seen so far
    {
     nord++;
     fOrdered->AddAt(m,j); // add module at the end
     break; // go for next matrix module
    }
 
    if (signal < ((AliCalmodule*)fOrdered->At(j))->GetSignal(1,1)) continue;
 
    nord++;
    for (Int_t k=nord-1; k>j; k--) // create empty position
    {
     fOrdered->AddAt(fOrdered->At(k-1),k);
    }
    fOrdered->AddAt(m,j); // put module at empty position
    break; // go for next matrix module
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SortA()
{
// Order the modules with decreasing signal by looping over the linear array
// of fired modules.
// Modules which were declared as "Dead" will be rejected.
// The gain etc... corrected module signals will be used in the ordering process.
//
// Note : This method is rather fast even for large, very finely granulated calorimeters.
//
// Very specific case :
// ====================
// In case of various overlapping showers of which the central modules have
// EXACTLY the same signal this ordering procedure may have the following
// advantages and disadvantages.
//
// Advantages :
// ------------
// * Even in case of a very high occupancy, the resulting cluster signals
//   will in general NOT depend on the geometrical location in the detector matrix.
//
// Disadvantages :
// ---------------
// * In case of multi-overlapping showers, the central shower might be
//   "eaten-up" from both sides, resulting in a slightly too low value
//   of the resulting cluster signal.
// * This method might produce results depending on the filling
//   order of the matrix modules.

 SortHits();
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddRing(Int_t row, Int_t col, Int_t n)
{
// Add module signals of 1 ring around (row,col) to current cluster.
// The gain etc... corrected module signals will be used in this process.
// The parameter n denotes the maximum number of rings around cluster center.
// Note : This function is used recursively.

 if (!fMatrix) return;

 Int_t nrows=fMatrix->GetMaxRow();
 Int_t ncols=fMatrix->GetMaxColumn();
 
 if (n >= 1) // Check if any rings left for recursive calls
 {
  Float_t signal=GetSignal(row,col,1); // Gain etc... corrected signal of (row,col) module
 
  Int_t lrow=row-1; if (lrow < 1) lrow=1;         // row lowerbound for ring
  Int_t urow=row+1; if (urow > nrows) urow=nrows; // row upperbound for ring
  Int_t lcol=col-1; if (lcol < 1) lcol=1;         // col lowerbound for ring
  Int_t ucol=col+1; if (ucol > ncols) ucol=ncols; // row upperbound for ring
 
  for (Int_t i=lrow; i<=urow; i++)
  {
   for (Int_t j=lcol; j<=ucol; j++)
   {
    // add module(i,j) to cluster if the signal <= signal(row,col)
    if (GetSignal(i,j,1) <= signal)
    {
     AliCalmodule* m=(AliCalmodule*)fMatrix->GetObject(i,j);
     if (m) ((AliCalcluster*)fClusters->At(GetNclusters()-1))->Add(*m);
    }
    AddRing(i,j,n-1); // Go for ring of modules around this (i,j) one
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNclusters() const
{
// Provide the number of clusters
 Int_t nclu=0;
 if (fClusters) nclu=fClusters->GetEntries();
 return nclu;
}
///////////////////////////////////////////////////////////////////////////
AliCalcluster* AliCalorimeter::GetCluster(Int_t j) const
{
// Provide cluster number j
// Note : j=1 denotes the first cluster

 if (!fClusters) return 0;

 if ((j >= 1) && (j <= GetNclusters()))
 {
  return (AliCalcluster*)fClusters->At(j-1);
 }
 else
 {
  cout << " *AliCalorimeter::GetCluster* cluster number : " << j
       << " out of range ==> 0 returned." << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule* AliCalorimeter::GetModule(Int_t j) const
{
// Provide 'fired' module number j
// Note : j=1 denotes the first 'fired' module

 AliCalmodule* m=(AliCalmodule*)GetHit(j);
 return m;
}
///////////////////////////////////////////////////////////////////////////
AliCalmodule* AliCalorimeter::GetModule(Int_t row,Int_t col)
{
// Provide access to module (row,col).
// Note : first module is at (1,1).

 AliCalmodule* m=0;
 if (!fMatrix) LoadMatrix();
 if (fMatrix) m=(AliCalmodule*)fMatrix->GetObject(row,col);
 return m;
}
///////////////////////////////////////////////////////////////////////////
TH2F* AliCalorimeter::DrawModules(Float_t thresh,Int_t mode)
{
// Provide a lego plot of the module signals.
// The input parameter mode (default mode=0) has the same meaning as
// specified in the memberfunction GetSignal(row,col,mode).
// Only modules with a (corrected) signal value above the threshold
// (default thresh=0) will be displayed.

 Int_t nrows=fNrows;
 Int_t ncols=fNcolumns;

 if (!fMatrix) LoadMatrix();

 if (fMatrix && !nrows && !ncols)
 {
  nrows=fMatrix->GetMaxRow();
  ncols=fMatrix->GetMaxColumn();
 }
 
 if (fHmodules)
 {
  fHmodules->Reset();
 }
 else
 {
  fHmodules=new TH2F("fHmodules","Module signals",
            ncols,0.5,float(ncols)+0.5,nrows,0.5,float(nrows)+0.5);
 
  fHmodules->SetDirectory(0); // Suppress global character of histo pointer
 }
 
 Int_t nmods=GetNsignals();

 Int_t row,col;
 Float_t signal;
 Int_t dead;
 for (Int_t i=1; i<=nmods; i++)
 {
  AliCalmodule* m=(AliCalmodule*)fMatrix->GetObject(i);
  if (m)
  {
   row=m->GetRow();
   col=m->GetColumn();
   dead=m->GetDeadValue();
   signal=0;
   if (!dead) signal=GetSignal(row,col,mode);
   if (signal>thresh) fHmodules->Fill(float(col),float(row),signal);
  }
 }
 
 fHmodules->Draw("lego");
 return fHmodules;
}
///////////////////////////////////////////////////////////////////////////
TH2F* AliCalorimeter::DrawClusters(Float_t thresh)
{
// Provide a lego plot of the cluster signals.
// Only clusters with a signal value above the threshold (default thresh=0)
// will be displayed.

 Int_t nrows=fNrows;
 Int_t ncols=fNcolumns;

 if (!fMatrix) LoadMatrix();

 if (fMatrix && !nrows && !ncols)
 {
  nrows=fMatrix->GetMaxRow();
  ncols=fMatrix->GetMaxColumn();
 }
 
 if (fHclusters)
 {
  fHclusters->Reset();
 }
 else
 {
  fHclusters=new TH2F("fHclusters","Cluster signals",
            ncols,0.5,float(ncols)+0.5,nrows,0.5,float(nrows)+0.5);
 
  fHclusters->SetDirectory(0); // Suppress global character of histo pointer
 }
 
 AliCalcluster* c;
 Int_t row,col;
 Float_t signal;
 for (Int_t i=0; i<GetNclusters(); i++)
 {
  c=(AliCalcluster*)fClusters->At(i);
  if (c)
  {
   row=c->GetRow();
   col=c->GetColumn();
   signal=c->GetSignal();
   if (signal>thresh) fHclusters->Fill(float(col),float(row),signal);
  }
 }
 
 fHclusters->Draw("lego");
 return fHclusters;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::Ungroup()
{
// Set the module signals back to the non-clustered situation

 if (!fMatrix) LoadMatrix();

 if (!fMatrix) return;
 
 Int_t nsig=GetNsignals();

 Float_t signal=0;
 for (Int_t j=1; j<=nsig; j++)
 {
  AliCalmodule* m=(AliCalmodule*)fMatrix->GetObject(j);
  if (m)
  {
   signal=m->GetSignal();
   m->SetClusteredSignal(signal);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::AddVetoSignal(AliSignal& s)
{
// Associate an (extrapolated) AliSignal as veto to the calorimeter.
 if (!fVetos)
 {
  fVetos=new TObjArray();
  fVetos->SetOwner();
 } 

 AliSignal* sx=new AliSignal(s);

 fVetos->Add(sx);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetNvetos() const
{
// Provide the number of veto signals associated to the calorimeter.
 Int_t nvetos=0;
 if (fVetos) nvetos=fVetos->GetEntries();
 return nvetos;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliCalorimeter::GetVetoSignal(Int_t i) const
{
// Provide access to the i-th veto signal of this calorimeter
// Note : The first hit corresponds to i=1

 if (i>0 && i<=GetNvetos())
 {
  return (AliSignal*)fVetos->At(i-1);
 }
 else
 {
  cout << " *AliCalorimeter::GetVetoSignal* Signal number " << i
       << " out of range ==> 0 returned." << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::SetMatrixSwapMode(Int_t swap)
{
// Set the swap mode for the module and position matrices.
// At invokation of this memberfunction the default argument is swap=1.
// For further details see the documentation of AliObjMatrix.
 if (swap==0 || swap==1)
 {
  fSwap=swap;
 }
 else
 {
  cout << " *AliCalorimeter::SetMatrixSwapMode* Invalid argument : swap = " << swap << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalorimeter::GetMatrixSwapMode() const
{
// Provide the swap mode for the module and position matrices.
// For further details see the documentation of AliObjMatrix.
 return fSwap;
}
///////////////////////////////////////////////////////////////////////////
void AliCalorimeter::LoadMatrix()
{
// Load the matrix lookup table of module pointers from the linear hit array.
 Int_t nhits=GetNhits();
 
 if (!nhits) return;

 fMatrix=new AliObjMatrix();
 fMatrix->SetSwapMode(fSwap);

 Int_t row=0;
 Int_t col=0;
 for (Int_t i=1; i<=nhits; i++)
 {
  AliCalmodule* m=(AliCalmodule*)GetHit(i);
  if (m)
  {
   row=m->GetRow();
   col=m->GetColumn();
   fMatrix->EnterObject(row,col,m);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliCalorimeter::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers like AliEvent when adding objects in case the
// container owns the objects.

 AliCalorimeter* cal=new AliCalorimeter(*this);
 if (name)
 {
  if (strlen(name)) cal->SetName(name);
 } 
 return cal;
}
///////////////////////////////////////////////////////////////////////////
