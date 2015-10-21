/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEmcalTriggerDataGrid.h"

/// \cond CLASSIMP
templateClassImp(AliEmcalTriggerDataGrid)
/// \endcond

template<typename T>
AliEmcalTriggerDataGrid<T>::AliEmcalTriggerDataGrid():
    fNCols(0),
    fNRows(0),
    fValues(NULL)
{
}

template<typename T>
AliEmcalTriggerDataGrid<T>::AliEmcalTriggerDataGrid(Int_t ncols, Int_t nrows):
    fNCols(ncols),
    fNRows(nrows),
    fValues(NULL)
{
  fValues = new T[fNCols * fNRows];
  memset(fValues, 0, sizeof(T) * fNCols * fNRows);
}

template<typename T>
AliEmcalTriggerDataGrid<T>::AliEmcalTriggerDataGrid(const AliEmcalTriggerDataGrid &ref):
    fNCols(ref.fNCols),
    fNRows(ref.fNRows),
    fValues(NULL)
{
  fValues = new T[fNCols * fNRows];
  memcpy(fValues, ref.fValues, sizeof(T) * fNCols * fNRows);
}

template<typename T>
AliEmcalTriggerDataGrid<T> &AliEmcalTriggerDataGrid<T>::operator=(const AliEmcalTriggerDataGrid<T> &ref){
  TObject::operator =(ref);
  if(this != &ref){
    if(fValues) delete[] fValues;
    fNCols = ref.fNCols;
    fNRows = ref.fNRows;
    fValues = new T[fNCols * fNRows];
    memcpy(fValues, ref.fValues, sizeof(T) * fNCols * fNRows);
  }
  return *this;
}

template<typename T>
const T &AliEmcalTriggerDataGrid<T>::operator()(Int_t col, Int_t row) const{
  if(!fValues) throw UninitException();
  if(row >= fNRows) throw OutOfBoundsException(OutOfBoundsException::kRowDir, row, fNRows);
  if(col >= fNCols) throw OutOfBoundsException(OutOfBoundsException::kColDir, col, fNCols);
  return fValues[GetIndex(col, row)];
}

template<typename T>
T &AliEmcalTriggerDataGrid<T>::operator()(Int_t col, Int_t row) {
  if(!fValues) throw UninitException();
  if(row >= fNRows) throw OutOfBoundsException(OutOfBoundsException::kRowDir, row, fNRows);
  if(col >= fNCols) throw OutOfBoundsException(OutOfBoundsException::kColDir, col, fNCols);
  return fValues[GetIndex(col, row)];
}

template<typename T>
AliEmcalTriggerDataGrid<T>::~AliEmcalTriggerDataGrid() {
  if(fValues) delete[] fValues;
}

template<typename T>
void AliEmcalTriggerDataGrid<T>::Allocate(Int_t ncols, Int_t nrows){
  if(fValues){
    T *tempstorage = new T[ncols * nrows];
    memset(tempstorage, 0, sizeof(T) * ncols * nrows);
    // Copy old content in new array
    for(int irow = 0; irow < fNRows; irow++){
      for(int icol = 0; icol < fNCols; icol++){
        int indexnew = irow * ncols + icol,
            indexold = irow * fNCols + icol;
        tempstorage[indexnew] = fValues[indexold];
      }
    }
    delete[] fValues;
    fValues = tempstorage;
  } else {
    fValues = new T[ncols * nrows];
    memset(fValues, 0, sizeof(T) * ncols * nrows);
  }
  fNCols = ncols;
  fNRows = nrows;
}

template<typename T>
void AliEmcalTriggerDataGrid<T>::SetADC(Int_t col, Int_t row, const T &adc) {
  if(!fValues)
    throw UninitException();
  if(row >= fNRows)
    throw OutOfBoundsException(OutOfBoundsException::kRowDir, row, fNRows);
  if(col >= fNCols)
    throw OutOfBoundsException(OutOfBoundsException::kColDir, col, fNCols);
  fValues[GetIndex(col, row)] = adc;
}

template<typename T>
void AliEmcalTriggerDataGrid<T>::Reset() {
  memset(fValues, 0, sizeof(T) * fNCols * fNRows);
}

template<typename T>
const T &AliEmcalTriggerDataGrid<T>::GetADC(Int_t col, Int_t row) const {
  if(!fValues)
    throw UninitException();
  if(row >= fNRows)
    throw OutOfBoundsException(OutOfBoundsException::kRowDir, row, fNRows);
  if(col >= fNCols)
    throw OutOfBoundsException(OutOfBoundsException::kColDir, col, fNCols);
  return fValues[GetIndex(col, row)];
}

template<typename T>
Int_t AliEmcalTriggerDataGrid<T>::GetIndex(Int_t col, Int_t row) const {
  return row * fNCols + col;
}

template class AliEmcalTriggerDataGrid<int>;
template class AliEmcalTriggerDataGrid<double>;
template class AliEmcalTriggerDataGrid<float>;
template class AliEmcalTriggerDataGrid<char>;
