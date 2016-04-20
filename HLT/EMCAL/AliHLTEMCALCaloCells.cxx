/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <AliHLTEMCALCaloCells.h>

/// \cond CLASSIMP
ClassImp(AliHLTEMCALCaloCells);
/// \endcond

/// Default constructor.
AliHLTEMCALCaloCells::AliHLTEMCALCaloCells() :
  AliVCaloCells(), fNCells(0), fCellNumber(0),
  fAmplitude(0), fType(kUndef), fIsSorted(kTRUE), fCapacity(0),
  fSwapCellNumber(0), fSwapAmplitude(0), fSwapIndexArray(0), fSwapCapacity(0)
{
}

/// Named constructor.
/// \param name Name of the object
/// \param title Title of the object
/// \param ttype Cell type (EMCal, PHOS)
AliHLTEMCALCaloCells::AliHLTEMCALCaloCells(const char* name, const char* title, VCells_t ttype) :
   AliVCaloCells(name, title), fNCells(0), fCellNumber(0),
   fAmplitude(0), fType(ttype), fIsSorted(kTRUE), fCapacity(0),
   fSwapCellNumber(0), fSwapAmplitude(0), fSwapIndexArray(0), fSwapCapacity(0)
{
}

/// Copy constructor.
/// \param c Const reference to copy from
AliHLTEMCALCaloCells::AliHLTEMCALCaloCells(const AliHLTEMCALCaloCells& c) :
  AliVCaloCells(c), fNCells(c.fNCells),  fCellNumber(0),
  fAmplitude(0), fType(c.fType), fIsSorted(c.fIsSorted), fCapacity(0),
  fSwapCellNumber(0), fSwapAmplitude(0), fSwapIndexArray(0), fSwapCapacity(0)
{
  fCellNumber = new Short_t[fNCells];
  fAmplitude  = new Double32_t[fNCells];

  for(Int_t i = 0; i < fNCells; i++) {
    fCellNumber[i]    = c.fCellNumber[i];
    fAmplitude[i]     = c.fAmplitude[i];
  }
}


/// Assignment operator.
/// \param source Const reference to copy from
/// \return Reference to this
AliHLTEMCALCaloCells & AliHLTEMCALCaloCells::operator =(const AliHLTEMCALCaloCells& source)
{
  if(this != &source)
  {
    AliVCaloCells::operator=(source);
    
    if (fCapacity != source.fCapacity) {
      delete [] fCellNumber;
      delete [] fAmplitude;
  
      fCellNumber = new Short_t[fCapacity];
      fAmplitude  = new Double32_t[fCapacity];
    }
    
    fNCells = source.fNCells;

    memcpy(fCellNumber,source.fCellNumber,fNCells*sizeof(Short_t));
    memcpy(fAmplitude, source.fAmplitude, fNCells*sizeof(Double32_t));

    fIsSorted = source.fIsSorted;
    fType = source.fType;
  }

  return *this;
}

/// This overwrites the virtual TObject::Copy()
/// to allow run time copying without casting
void AliHLTEMCALCaloCells::Copy(TObject &obj) const
{
  if(this==&obj)return;
  AliHLTEMCALCaloCells *robj = dynamic_cast<AliHLTEMCALCaloCells*>(&obj);
  if(!robj)return; // not an AliHLTCaloCells
  *robj = *this;
}

/// Copy the calo cells into a new object. If option all=FALSE, just the object type, 
/// for mixing.
AliVCaloCells* AliHLTEMCALCaloCells::CopyCaloCells(Bool_t all = kTRUE) const
{  
  AliVCaloCells *obj = new AliHLTEMCALCaloCells();
  
  if(all){
    obj->SetName (GetName()) ; 
    obj->SetTitle(GetTitle()) ; 
    obj->SetType (GetType()) ; 
    
    obj->SetNumberOfCells(fNCells);
    for (Short_t i = 0; i < fNCells; i++) {
      obj->SetCell(i,fCellNumber[i],fAmplitude[i],-1,-1,-1);
    }
  }

  return obj;
}

/// Destructor.
AliHLTEMCALCaloCells::~AliHLTEMCALCaloCells()
{
  DeleteContainer();
}

/// Clear number of cells.
void AliHLTEMCALCaloCells::Clear(Option_t*)
{
  for (Int_t i = 0; i < fCapacity; i++) {
    fCellNumber[i] = 0;
    fAmplitude[i] = 0;
  }
  for (Int_t i = 0; i < fSwapCapacity; i++) {
    fSwapCellNumber[i] = 0;
    fSwapAmplitude[i] = 0;
    fSwapIndexArray[i] = 0;
  }

  fNCells = 0;
}

/// Function that creates container to store calorimeter cell data.
/// The current content of the container is reset.
/// The container is only enlarged: if the container already exists and
/// it is bigger or the same size of nCells than the content container will only
/// be reset.
/// \param nCells capacity of the container
void AliHLTEMCALCaloCells::CreateContainer(Short_t nCells)
{
  fNCells = 0;
  
  if (nCells <= 0) {
    fCapacity = 0;
    return;
  }

  if (nCells > fCapacity) {
    DeleteContainer();
    fCellNumber = new Short_t[nCells];
    fAmplitude  = new Double32_t[nCells];
    fCapacity = nCells;
  }

  for (Int_t i = 0; i < fCapacity; i++) {
    fCellNumber[i] = 0;
    fAmplitude[i] = 0;
  }
}

/// Deletes allocated memory
void AliHLTEMCALCaloCells::DeleteContainer()
{
  if (fCellNumber) {
    delete[] fCellNumber;
    fCellNumber = 0;
  }

  if (fAmplitude) {
    delete[] fAmplitude;
    fAmplitude = NULL;
  }

  if (fSwapCellNumber) {
    delete[] fSwapCellNumber;
    fSwapCellNumber = 0;
  }

  if (fSwapAmplitude) {
    delete[] fSwapAmplitude;
    fSwapAmplitude = NULL;
  }

  if (fSwapIndexArray) {
    delete[] fSwapIndexArray;
    fSwapIndexArray = NULL;
  }

  fSwapCapacity = 0;
  fCapacity = 0;
  fNCells = 0;
  fIsSorted = kFALSE;
}

/// Sort the cell array by cell number.
void AliHLTEMCALCaloCells::Sort()
{
  if (fSwapCapacity < fCapacity) {
    delete[] fSwapCellNumber;
    delete[] fSwapAmplitude;
    fSwapCellNumber = new Short_t[fCapacity];
    fSwapAmplitude = new Double32_t[fCapacity];
    fSwapIndexArray = new Int_t[fCapacity];
    fSwapCapacity = fCapacity;

    for (Int_t i = 0; i < fSwapCapacity; i++) {
      fSwapCellNumber[i] = 0;
      fSwapAmplitude[i] = 0;
      fSwapIndexArray[i] = -1;
    }
  }

  TMath::Sort(fNCells, fCellNumber, fSwapIndexArray, kFALSE);
  
  for (Int_t i = 0; i < fNCells; i++) {
    fSwapCellNumber[i] = fCellNumber[fSwapIndexArray[i]];
    fSwapAmplitude[i] = fAmplitude[fSwapIndexArray[i]];
  }
  
  Short_t *tempCellNumber = fSwapCellNumber;
  Double32_t *tempAmplitude = fSwapAmplitude;
  Int_t tempCapacity = fSwapCapacity;

  fSwapCellNumber = fCellNumber;
  fSwapAmplitude = fAmplitude;
  fSwapCapacity = fCapacity;

  fCellNumber = tempCellNumber;
  fAmplitude = tempAmplitude;
  fCapacity = tempCapacity;
  
  fIsSorted = kTRUE;
} 

/// Sets a cell at the given position.
/// \param pos: cell position in array.
/// \param cellNumber: Cell absolute Id. number.
/// \param amplitude: Cell signal (GeV).
/// \return true on success
Bool_t AliHLTEMCALCaloCells::SetCell(Short_t pos, Short_t cellNumber, Double32_t amplitude,
    Double32_t /*time*/, Int_t /*mclabel*/, Double32_t /*efrac*/, Bool_t /*isHG*/)
{
  if (pos < 0) return kFALSE;
  if (pos >= fNCells) SetNumberOfCells(pos+1);

  fCellNumber[pos] = cellNumber;
  fAmplitude[pos]  = amplitude;

  fIsSorted = kFALSE;
  return kTRUE;
}

/// Set the number of cells. If needed, enlarges the container (container is never shrunk)
/// \param n Number of cells
void AliHLTEMCALCaloCells::SetNumberOfCells(Int_t n)
{
  Int_t min = 0;
  Int_t max = 0;

  if (fNCells == n) {
    return;
  }
  else if (fNCells > n) {
    min = n;
    max = fNCells;
  }
  else {
    min = fNCells;
    max = n;
  }

  if (fCapacity < n) {
    // We need to enlarge the container
    Short_t *tempCellNumber = fCellNumber;
    Double32_t *tempAmplitude = fAmplitude;
    // Prevent deallocation in CreateContainer
    fCellNumber = 0;
    fAmplitude = 0;
    CreateContainer(n);
    for (Int_t i = 0; i < min; i++) {
      fCellNumber[i] = tempCellNumber[i];
      fAmplitude[i] = tempAmplitude[i];
    }
    delete[] tempCellNumber;
    delete[] tempAmplitude;
  }
  else {
    // The container is big enough

    // Set to zero extra cells
    for (Int_t i = min; i < max; i++){
      fCellNumber[i] = 0;
      fAmplitude[i] = 0;
    }
  }

  fNCells = n;
}
