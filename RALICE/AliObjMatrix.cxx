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
// Class AliObjMatrix
// Handling of a matrix structure of objects.
// All objects which are derived from TObject may be entered into the matrix 
// structure. This means that also TObjArray objects can be entered,
// which implies an increase of the dimension of the resulting structure.
//
// Example :
// =========
//
// AliObjMatrix* matrix=new AliObjMatrix();
// matrix->SetOwner();
// matrix->SetSwapMode();
//
// Float_t pos[3];
//
// AliSignal* s=0;
//
// s=new AliSignal();
// s->SetSignal(135);
// pos[0]=-120.4
// pos[1]=78.25
// pos[3]=12.93
// s->SetPosition(pos,"car");
// matrix->EnterObject(6,21,s);
//
// s=new AliSignal();
// s->SetSignal(25.84);
// pos[0]=68.7
// pos[1]=-53.88
// pos[3]=22.69
// s->SetPosition(pos,"car");
// matrix->EnterObject(8,13,s);
//
// s=new AliSignal();
// s->SetSignal(87.25);
// pos[0]=154.8
// pos[1]=932.576
// pos[3]=-1382.754
// s->SetPosition(pos,"car");
// matrix->EnterObject(64,3,s);
//
// Int_t nrows=matrix->GetMaxRow();
// Int_t ncols=matrix->GetMaxColumn();
//
//  cout << " Maxrow : " << nrows << " Maxcol : " << ncols
//       << " Nobjects : " << matrix->GetNobjects() << endl;
//
//  for (Int_t i=1; i<=nrows; i++)
//  {
//   for (Int_t j=1; j<=ncols; j++)
//   {
//    s=(AliSignal*)matrix->GetObject(i,j);
//    if (s) cout << " At (" << i << "," << j << ") Signal : " << s->GetSignal() << endl;
//   }
//  }
//
//--- Author: Nick van Eijndhoven 23-jan-2003 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliObjMatrix.h"
#include "Riostream.h"
 
ClassImp(AliObjMatrix) // Class implementation to enable ROOT I/O
 
AliObjMatrix::AliObjMatrix() : TObject()
{
// Default constructor.
// Note : The owner and swap mode flags will be initialised to 0.
//        See the memberfunctions SetOwner() and SetSwapMode() for further
//        details. 
 fRows=0;
 fOwn=0;
 fSwap=0;
 fMaxrow=0;
 fMaxcol=0;
 fObjects=0;
}
///////////////////////////////////////////////////////////////////////////
AliObjMatrix::~AliObjMatrix()
{
// Default destructor.
 if (fRows)
 {
  delete fRows;
  fRows=0;
 }
 if (fObjects)
 {
  delete fObjects;
  fObjects=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliObjMatrix::Reset()
{
// Reset the whole matrix structure.
// Note : The values of the owner and swap mode flags will not be modified.
//        To modify the ownership, use the memberfunction SetOwner(). 
//        To modify the swap mode, use the memberfunction SetSwapMode(). 
 if (fRows)
 {
  delete fRows;
  fRows=0;
 }
 if (fObjects)
 {
  delete fObjects;
  fObjects=0;
 }

 fMaxrow=0;
 fMaxcol=0;
}
///////////////////////////////////////////////////////////////////////////
void AliObjMatrix::SetOwner(Int_t own)
{
// Set the owner flag (0/1) for the stored objects.
// When the owner flag is set to 1, all entered objects are owned by the
// matrix structure.
// At invokation of this memberfunction the default argument is own=1.
//
 fOwn=own;

 if (!fRows) return;

 for (Int_t irow=0; irow<fRows->GetSize(); irow++)
 {
  TObjArray* mrow=(TObjArray*)fRows->At(irow);
  if (mrow)
  {
   if (own)
   {
    mrow->SetOwner(kTRUE);
   }
   else
   {
    mrow->SetOwner(kFALSE);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetOwner()
{
// Provide the owner flag for the stored objects.
 return fOwn;
}
///////////////////////////////////////////////////////////////////////////
void AliObjMatrix::SetSwapMode(Int_t swap)
{
// Set the swap mode flag (0/1) for the internal matrix storage.
// In case the number of rows differs considerably from the number of columns,
// it might be more efficient (w.r.t. memory usage and/or output file size)
// to internally store the matrix with the rows and colums swapped.
// This swapping is only related with the internal storage and as such
// is completely hidden for the user.
// At invokation of this memberfunction the default argument is swap=1.
//
// Note : The swap mode can only be set as long as no objects have
//        been stored in the matrix structure (i.e. a new instance
//        of AliObjMatrix or after invokation of the Reset() function). 
//
 if (!fRows)
 {
  fSwap=swap;
 }
 else
 {
  cout << " *AliObjMatrix::SetSwapMode* Matrix not empty ==> No action." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetSwapMode()
{
// Provide the swap mode flag for this matrix.
 return fSwap;
}
///////////////////////////////////////////////////////////////////////////
void AliObjMatrix::EnterObject(Int_t row,Int_t col,TObject* obj)
{
// Enter an object to the matrix structure at location (row,col).
// In case the location already contained an object, the existing object
// will first be removed before the new object is stored.
// According to the status of the owner flag (see the SetOwner() function)
// the existing object will also be deleted.
// Note : The first location in the matrix is indicated as (1,1).
 if (row<1 || col<1)
 {
  cout << " *AliObjMatrix::AddObject* Invalid argument(s) (row,col) : ("
       << row << "," << col << ")" << endl;
  return;
 }

 if (row>fMaxrow) fMaxrow=row;
 if (col>fMaxcol) fMaxcol=col;

 Int_t rowx=row;
 if (fSwap) rowx=col;
 Int_t colx=col;
 if (fSwap) colx=row;

 if (!fRows)
 {
  fRows=new TObjArray(rowx);
  fRows->SetOwner();
 }
 else
 {
  if (rowx > fRows->GetSize()) fRows->Expand(rowx);
 }

 TObjArray* mrow=(TObjArray*)fRows->At(rowx-1);

 if (!mrow)
 {
  TObjArray* columns=new TObjArray(colx);
  if (fOwn) columns->SetOwner();
  fRows->AddAt(columns,rowx-1);
  mrow=columns;
 }
 else
 {
  if (colx > mrow->GetSize()) mrow->Expand(colx);
 }

 TObject* old=(TObject*)mrow->At(colx-1);
 if (old)
 {
  fObjects->Remove(old);
  fObjects->Compress();
  if (fOwn) delete old;
 }

 mrow->AddAt(obj,colx-1);

 if (!fObjects) fObjects=new TObjArray();
 fObjects->Add(obj);
}
///////////////////////////////////////////////////////////////////////////
void AliObjMatrix::RemoveObject(Int_t row,Int_t col)
{
// Remove the object stored at the matrix location (row,col).
// In case the object was owned by the matrix, it will be deleted.
//
// Note : The first location in the matrix is indicated as (1,1).

 TObject* obj=0;

 if (!fRows || row<1 || col<1) return;


 Int_t rowx=row;
 if (fSwap) rowx=col;
 Int_t colx=col;
 if (fSwap) colx=row;

 TObjArray* mrow=0;
 if (rowx <= fRows->GetSize()) mrow=(TObjArray*)fRows->At(rowx-1);

 if (!mrow) return;

 if (colx <= mrow->GetSize()) obj=(TObject*)mrow->At(colx-1);

 if (obj)
 {
  fObjects->Remove(obj);
  fObjects->Compress();
  mrow->Remove(obj);
  if (fOwn) delete obj;
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliObjMatrix::GetObject(Int_t row,Int_t col)
{
// Provide a pointer to the object stored at the matrix location (row,col).
// In case no object was stored at the indicated location or the location
// would reside outside the matrix boundaries, a value 0 will be returned.
// Note : The first location in the matrix is indicated as (1,1).

 TObject* obj=0;

 if (!fRows || row<1 || col<1) return obj;


 Int_t rowx=row;
 if (fSwap) rowx=col;
 Int_t colx=col;
 if (fSwap) colx=row;

 TObjArray* mrow=0;
 if (rowx <= fRows->GetSize()) mrow=(TObjArray*)fRows->At(rowx-1);

 if (!mrow) return obj;

 if (colx <= mrow->GetSize()) obj=(TObject*)mrow->At(colx-1);

 return obj;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliObjMatrix::GetObject(Int_t j)
{
// Provide a pointer to the j-th stored object.
// In case the index j is invalid, a value 0 will be returned.
// The first stored object is indicated as j=1.
//
// Note : Do NOT delete the object.
//        To remove an object, the memberfunction RemoveObject()
//        should be used.

 TObject* obj=0;
 Int_t nobj=0;
 if (fObjects) nobj=fObjects->GetSize();

 if (j>0 && j<=nobj) obj=(TObject*)fObjects->At(j-1);

 return obj;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliObjMatrix::GetObjects()
{
// Provide references to all the stored objects.
// In case no objects are present, a value 0 will be returned.
//
// Note : Do NOT make any changes to the reference array apart from
//        changing the order of the pointers of the various objects.
//        For addition or removal of objects, the memberfunctions
//        EnterObject() and RemoveObject() should be used.

 return fObjects;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetMaxRow()
{
// Provide the maximum row number index.
 return fMaxrow;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetMaxColumn()
{
// Provide the maximum column number index.
 return fMaxcol;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetNobjects()
{
// Provide the number of stored objects.
 Int_t nobj=0;
 if (fObjects) nobj=fObjects->GetEntries();

 return nobj;
}
///////////////////////////////////////////////////////////////////////////
