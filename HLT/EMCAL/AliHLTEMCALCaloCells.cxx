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
  AliVCaloCells(), fCapacity(0), fNCells(0),
  fAmplitude(0), fType(kUndef)
{
}

/// Named constructor.
/// \param name Name of the object
/// \param title Title of the object
/// \param ttype Cell type (EMCal, PHOS)
AliHLTEMCALCaloCells::AliHLTEMCALCaloCells(const char* name, const char* title, VCells_t ttype) :
   AliVCaloCells(name, title), fCapacity(0), fNCells(0),
   fAmplitude(0), fType(kUndef)
{
}

/// Copy constructor.
/// \param c Const reference to copy from
AliHLTEMCALCaloCells::AliHLTEMCALCaloCells(const AliHLTEMCALCaloCells& c) :
  AliVCaloCells(c), fCapacity(c.fCapacity), fNCells(c.fNCells),
  fAmplitude(0), fType(kUndef)
{
  fAmplitude = new Double_t[fCapacity];

  memcpy(fAmplitude, c.fAmplitude, sizeof(Double_t)*fCapacity);
}


/// Assignment operator.
/// \param source Const reference to copy from
/// \return Reference to this
AliHLTEMCALCaloCells & AliHLTEMCALCaloCells::operator =(const AliHLTEMCALCaloCells& source)
{
  if (this != &source) {
    AliVCaloCells::operator=(source);
    
    if (fCapacity != source.fCapacity) {
      delete[] fAmplitude;
      fCapacity = source.fCapacity;

      fAmplitude = new Double_t[fCapacity];
    }

    fNCells = source.fNCells;
    memcpy(fAmplitude, source.fAmplitude, sizeof(Double_t)*fCapacity);
  }

  return *this;
}

/// This overwrites the virtual TObject::Copy()
/// to allow run time copying without casting
void AliHLTEMCALCaloCells::Copy(TObject &obj) const
{
  if(this==&obj)return;
  AliHLTEMCALCaloCells *robj = dynamic_cast<AliHLTEMCALCaloCells*>(&obj);
  if (!robj) return; // not an AliHLTCaloCells
  *robj = *this;
}

/// Copy the calo cells into a new object. If option all=FALSE, just the object type, 
/// for mixing.
AliVCaloCells* AliHLTEMCALCaloCells::CopyCaloCells(Bool_t all = kTRUE) const
{  
  AliVCaloCells *obj = new AliHLTEMCALCaloCells(*this);

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
  TObject::Clear();

  memset(fAmplitude, 0, sizeof(Double_t)*fCapacity);

  fNCells = 0;
}

/// Function that creates container to store calorimeter cell data.
/// The current content of the container is reset.
/// \param nCells capacity of the container
void AliHLTEMCALCaloCells::CreateContainer(Short_t nCells)
{
  if (fCapacity != nCells) {
    delete[] fAmplitude;
    fCapacity = nCells;
    fAmplitude = new Double_t[fCapacity];
  }

  Clear();
}

/// Deletes allocated memory
void AliHLTEMCALCaloCells::DeleteContainer()
{
  delete[] fAmplitude;
  fAmplitude = 0;
  fCapacity = 0;
  fNCells = 0;
}

/// Sets a cell at the given position.
/// \param pos: cell position in array.
/// \param amplitude: Cell signal (GeV).
/// \return true on success
Bool_t AliHLTEMCALCaloCells::SetCell(Short_t pos, Short_t /*cellNumber*/, Double_t amplitude,
    Double_t /*time*/, Int_t /*mclabel*/, Double_t /*efrac*/, Bool_t /*isHG*/)
{
  if (pos < 0 || pos >= fCapacity) return kFALSE;

  if (amplitude == 0 && fAmplitude[pos] > 0 && fNCells > 0) {
    fNCells--;
  }
  else if (amplitude > 0 && fAmplitude[pos] == 0) {
    fNCells++;
  }

  fAmplitude[pos] = amplitude;

  return kTRUE;
}

/// Add energy to a cell at the given position.
/// \param pos: cell position in array.
/// \param amplitude: Cell signal (GeV).
/// \return true on success
Bool_t AliHLTEMCALCaloCells::AddCell(Short_t pos, Double_t amplitude)
{
  if (pos < 0 || pos >= fCapacity) {
    Printf(Form("ERROR-AliHLTEMCALCaloCells::AddCell - Could not add energy of cell %d (capacity is %d)", pos, fCapacity));
    return kFALSE;
  }

  if (amplitude > 0 && fAmplitude[pos] == 0) {
    fNCells++;
  }

  fAmplitude[pos] += amplitude;

  return kTRUE;
}

/// Function that creates container to store calorimeter cell data.
/// The current content of the container is reset.
/// \param nCells capacity of the container
void AliHLTEMCALCaloCells::SetNumberOfCells(Int_t n)
{
  CreateContainer(n);
}
