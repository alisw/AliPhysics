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
 
AliObjMatrix::AliObjMatrix() : TNamed()
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
AliObjMatrix::AliObjMatrix(const AliObjMatrix& m) : TNamed(m)
{
// Copy constructor

 fRows=0;
 fMaxrow=0;
 fMaxcol=0;
 fObjects=0;

 fOwn=m.fOwn;
 fSwap=m.fSwap;

 Int_t maxrow=m.GetMaxRow();
 Int_t maxcol=m.GetMaxColumn();
 for (Int_t irow=1; irow<=maxrow; irow++)
 {
  for (Int_t icol=1; icol<=maxcol; icol++)
  {
   TObject* obj=m.GetObject(irow,icol);
   if (obj)
   {
    if (!fOwn)
    {
     EnterObject(irow,icol,obj);
    }
    else
    {
     EnterObject(irow,icol,obj->Clone());
    }
   }
  }
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
Int_t AliObjMatrix::GetOwner() const
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
Int_t AliObjMatrix::GetSwapMode() const
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
void AliObjMatrix::RemoveObjects(TObject* obj,Int_t row,Int_t col)
{
// Remove object(s) from the matrix according to user specified selections.
// In case the object was owned by the matrix, it will be deleted.
//
// An object is only removed from the matrix if the stored reference matches
// the argument "obj".
// In case obj=0 no check on the matching of the stored reference is performed
// and the stored object is always removed in accordance with the other
// selection criteria.
//
// In case the argument "row" is specified, only the object references from
// that matrix row will be deleted.
// In case row=0 (default) no checking on the row index is performed.
//
// In case the argument "col" is specified, only the object references from
// that matrix column will be deleted.
// In case col=0 (default) no checking on the column index is performed.
//
// So, invokation of RemoveObjects(obj) will remove all references to the
// object "obj" from the total matrix, whereas RemoveObjects(obj,0,col)
// will remove all references to the object "obj" only from column "col".
//
// Notes :
// -------
// The first location in the matrix is indicated as (1,1).
//
// Invokation of RemoveObjects(0,row,col) is equivalent to invoking the
// memberfunction RemoveObject(row,col).
// Invoking the latter directly is slightly faster.
//
// Invokation of RemoveObjects(0) is equivalent to invoking Reset().
// Invoking the latter directly is slightly faster.
//
 TArrayI rows;
 TArrayI cols;
 Int_t nrefs=0;

 if (row && col)
 {
  if (!obj)
  {
   RemoveObject(row,col);
  }
  else
  {
   TObject* objx=GetObject(row,col);
   if (objx==obj) RemoveObject(row,col);
  }
  return;
 }

 if (!row && !col)
 {
  if (!obj)
  {
   Reset();
   return;
  }
  else
  {
   nrefs=GetIndices(obj,rows,cols);
  }
 }

 if (row && !col) nrefs=GetIndices(obj,row,cols);
 if (!row && col) nrefs=GetIndices(obj,rows,col);

 // Remove the selected objects based on the obtained row and column indices
 Int_t irow,icol;
 for (Int_t i=0; i<nrefs; i++)
 {
  irow=row;
  if (!irow) irow=rows.At(i);
  icol=col;
  if (!icol) icol=cols.At(i);
  RemoveObject(irow,icol);
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliObjMatrix::GetObject(Int_t row,Int_t col) const
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
TObject* AliObjMatrix::GetObject(Int_t j) const
{
// Provide a pointer to the j-th stored object.
// In case the index j is invalid, a value 0 will be returned.
// The first stored object is indicated as j=1.
//
// Note : Do NOT delete the object.
//        To remove an object, the memberfunction RemoveObject() or
//        RemoveObjects() should be used.

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
//        EnterObject(), RemoveObject() or RemoveObjects() should be used.

 return fObjects;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetMaxRow() const
{
// Provide the maximum row number index.
 return fMaxrow;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetMaxColumn() const
{
// Provide the maximum column number index.
 return fMaxcol;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetNobjects() const
{
// Provide the number of stored objects.
 Int_t nobj=0;
 if (fObjects) nobj=fObjects->GetEntries();

 return nobj;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetNrefs(TObject* obj) const
{
// Provide the number of stored references to the specified object.
// If obj=0 the total number of stored references for all objects is returned.
 Int_t nobjs=GetNobjects();

 if (!obj) return nobjs;

 Int_t nrefs=0;
 for (Int_t i=1; i<=nobjs; i++)
 {
  TObject* objx=GetObject(i);
  if (objx==obj) nrefs++;
 }
 return nrefs;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetIndices(TObject* obj,TArrayI& rows,TArrayI& cols) const
{
// Provide the (row,col) indices of all the storage locations of the
// specified object.
// The row and column indices are returned in the two separate TArrayI arrays
// from which the (row,col) pairs can be obtained from the corresponding
// array indices like (row,col)=(rows.At(j),cols.At(j)).
// The integer return argument represents the number of (row,col) pairs which
// were encountered for the specified object.
//
// If obj=0 no object selection is performed and all (row,col) indices
// of the stored references for all objects are returned.
//
// Notes :
// -------
// As usual the convention is that row and column numbering starts at 1.
// 
// This memberfunction always resets the two TArrayI arrays at the start.
//
// This memberfunction can only be used to obtain the (row,col) indices
// of the object as stored via the EnterObject() memberfunction.
// This means that in case the user has entered a TObjArray as object
// (to increase the dimension of the resulting structure), the (row,col)
// indices of that TObjArray are obtained and NOT the indices of the
// actual objects contained in that TObjArray structure.
//
 Int_t nrefs=GetNrefs(obj);
 rows.Reset();
 cols.Reset();
 rows.Set(nrefs);
 cols.Set(nrefs);
 if (!nrefs) return 0;

 Int_t irow,icol;
 Int_t jref=0;
 for (Int_t i=0; i<fRows->GetSize(); i++)
 {
  TObjArray* columns=(TObjArray*)fRows->At(i);
  if (!columns) continue;

  for (Int_t j=0; j<columns->GetSize(); j++)
  {
   TObject* objx=(TObject*)columns->At(j);
   if (objx && (objx==obj || !obj))
   {
    irow=i+1;
    if (fSwap) irow=j+1;
    icol=j+1;
    if (fSwap) icol=i+1;
    rows.AddAt(irow,jref);
    cols.AddAt(icol,jref);
    jref++;
   }
   // All references found ==> Done
   if (jref==nrefs) break;
  }
  // All references found ==> Done
  if (jref==nrefs) break;
 }
 return nrefs;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetIndices(TObject* obj,Int_t row,TArrayI& cols) const
{
// Provide the column indices of all the storage locations of the
// specified object in the specified row of the matrix.
// The column indices are returned in the TArrayI array.
// The integer return argument represents the number of storage locations which
// were encountered for the specified object in the specified matrix row.
//
// If obj=0 no object selection is performed and all column indices
// of the stored references for all objects in this specified matrix row
// are returned.
//
// If row=0 all rows will be scanned and all column indices matching the
// object selection are returned.
// Note that in this case multiple appearances of the same column index
// will only be recorded once in the returned TArrayI array.
//
// Notes :
// -------
// As usual the convention is that row and column numbering starts at 1.
// 
// This memberfunction always resets the TArrayI array at the start.
//
// This memberfunction can only be used to obtain the column indices
// of the object as stored via the EnterObject() memberfunction.
// This means that in case the user has entered a TObjArray as object
// (to increase the dimension of the resulting structure), the column
// indices of that TObjArray are obtained and NOT the indices of the
// actual objects contained in that TObjArray structure.
//
 cols.Reset();

 if (row<0 || row>GetMaxRow()) return 0;

 Int_t nrefs=GetNrefs(obj);
 cols.Set(nrefs);
 if (!nrefs) return 0;

 Int_t irow,icol;
 Int_t jref=0;

 // No specific row selection
 if (!row)
 {
  TArrayI ar;
  TArrayI ac;
  Int_t n=GetIndices(obj,ar,ac);
  Int_t found=0;
  for (Int_t idx=0; idx<n; idx++)
  {
   icol=ac.At(idx);
   found=0;
   for (Int_t k=0; k<jref; k++)
   {
    if (icol==cols.At(k)) found=1;
   }
   if (!found)
   {
    cols.AddAt(icol,jref);
    jref++;
   }
  }
  // Set the array size to the actual number of different column indices
  cols.Set(jref);

  return jref;
 }

 // Specific row selection
 for (Int_t i=0; i<fRows->GetSize(); i++)
 {
  TObjArray* columns=(TObjArray*)fRows->At(i);
  if (!columns) continue;

  for (Int_t j=0; j<columns->GetSize(); j++)
  {
   TObject* objx=(TObject*)columns->At(j);
   if (objx && (objx==obj || !obj))
   {
    irow=i+1;
    if (fSwap) irow=j+1;
    icol=j+1;
    if (fSwap) icol=i+1;
    if (irow==row)
    {
     cols.AddAt(icol,jref);
     jref++;
    }
   }
   // All references found ==> Done
   if (jref==nrefs) break;
  }
  // All references found ==> Done
  if (jref==nrefs) break;
 }
 // Set the array size to the actual number of found occurrences
 cols.Set(jref);

 return jref;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliObjMatrix::GetIndices(TObject* obj,TArrayI& rows,Int_t col) const
{
// Provide the row indices of all the storage locations of the
// specified object in the specified column of the matrix.
// The row indices are returned in the TArrayI array.
// The integer return argument represents the number of storage locations which
// were encountered for the specified object in the specified matrix column.
//
// If obj=0 no object selection is performed and all row indices
// of the stored references for all objects in this specified matrix column
// are returned.
//
// If col=0 all columns will be scanned and all row indices matching the
// object selection are returned.
// Note that in this case multiple appearances of the same row index
// will only be recorded once in the returned TArrayI array.
//
// Notes :
// -------
// As usual the convention is that row and column numbering starts at 1.
// 
// This memberfunction always resets the TArrayI array at the start.
//
// This memberfunction can only be used to obtain the row indices
// of the object as stored via the EnterObject() memberfunction.
// This means that in case the user has entered a TObjArray as object
// (to increase the dimension of the resulting structure), the row
// indices of that TObjArray are obtained and NOT the indices of the
// actual objects contained in that TObjArray structure.
//
 rows.Reset();

 if (col<0 || col>GetMaxColumn()) return 0;

 Int_t nrefs=GetNrefs(obj);
 rows.Set(nrefs);
 if (!nrefs) return 0;

 Int_t irow,icol;
 Int_t jref=0;

 // No specific column selection
 if (!col)
 {
  TArrayI ar;
  TArrayI ac;
  Int_t n=GetIndices(obj,ar,ac);
  Int_t found=0;
  for (Int_t idx=0; idx<n; idx++)
  {
   irow=ar.At(idx);
   found=0;
   for (Int_t k=0; k<jref; k++)
   {
    if (irow==rows.At(k)) found=1;
   }
   if (!found)
   {
    rows.AddAt(irow,jref);
    jref++;
   }
  }
  // Set the array size to the actual number of different row indices
  rows.Set(jref);

  return jref;
 }

 // Specific column selection
 for (Int_t i=0; i<fRows->GetSize(); i++)
 {
  TObjArray* columns=(TObjArray*)fRows->At(i);
  if (!columns) continue;

  for (Int_t j=0; j<columns->GetSize(); j++)
  {
   TObject* objx=(TObject*)columns->At(j);
   if (objx && (objx==obj || !obj))
   {
    irow=i+1;
    if (fSwap) irow=j+1;
    icol=j+1;
    if (fSwap) icol=i+1;
    if (icol==col)
    {
     rows.AddAt(irow,jref);
     jref++;
    }
   }
   // All references found ==> Done
   if (jref==nrefs) break;
  }
  // All references found ==> Done
  if (jref==nrefs) break;
 }
 // Set the array size to the actual number of found occurrences
 rows.Set(jref);

 return jref;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliObjMatrix::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.

 AliObjMatrix* m=new AliObjMatrix(*this);
 if (name)
 {
  if (strlen(name)) m->SetName(name);
 }
 return m;
}
///////////////////////////////////////////////////////////////////////////
