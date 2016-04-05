#ifndef ALIHLTCALOCELLS_H
#define ALIHLTCALOCELLS_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
/// \class AliHLTEMCALCaloCells
/// \brief Class for calorimeter cell HLT data handling
///
/// HLT class to store calorimeter cell data.
/// This is a lightweigth version of AliESDCaloCells.
///
/// Data is stored in different arrays, each entry of the array corresponding to a cell.
/// The data stored is the cell energy, absolute id number.
///
///  \author Salvatore Aiola, <salvatore.aiola@cern.ch>, Yale University
///
//-------------------------------------------------------------------------

#include <AliVCaloCells.h>
#include <TMath.h>

class AliHLTEMCALCaloCells : public AliVCaloCells
{
 public:

  AliHLTEMCALCaloCells();
  AliHLTEMCALCaloCells(const char* name, const char* title, VCells_t ttype=kUndef);
  AliHLTEMCALCaloCells(const AliHLTEMCALCaloCells & cells);
  AliHLTEMCALCaloCells & operator=(const AliHLTEMCALCaloCells& source);
  virtual ~AliHLTEMCALCaloCells();
  
  virtual AliVCaloCells * CopyCaloCells(Bool_t all) const;
  virtual void    Copy(TObject &obj) const;
  void            Clear(Option_t* option = "");
  void            CreateContainer(Short_t nCells);
  void            DeleteContainer();
  void            Sort();
  
  Bool_t          IsEMCAL()  const { return (fType == kEMCALCell); }
  Bool_t          IsPHOS()   const { return (fType == kPHOSCell) ; }
  Char_t          GetType()  const { return  fType               ; }
  void            SetType(Char_t t){ fType = t                   ; }
  
  inline Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, Double_t &time, Int_t &mclabel,      Double_t &efrac) const;
  Bool_t          SetCell(Short_t pos, Short_t  cellNumber, Double_t  amplitude, Double_t  time, Int_t  mclabel = -1, Double_t  efrac = 0., Bool_t isHG=kFALSE);

  Short_t         GetNumberOfCells() const  { return fNCells ; }
  void            SetNumberOfCells(Int_t n);
  
  inline Double_t GetCellAmplitude(Short_t cellNumber);
  inline Bool_t   GetCellHighGain(Short_t cellNumber) { return kFALSE; }  //is this cell High Gain
  inline Short_t  GetCellPosition(Short_t cellNumber);
  inline Double_t GetCellTime(Short_t cellNumber) { return -1; }
  
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Bool_t   GetHighGain(Short_t pos) const { return kFALSE; }
  inline Double_t GetTime(Short_t pos) const { return -1; }
  inline Short_t  GetCellNumber(Short_t pos) const;

  // MC & embedding
  inline Int_t    GetCellMCLabel(Short_t cellNumber) { return -1; }
  inline Int_t    GetMCLabel(Short_t pos) const { return -1; }
  
  inline Double_t GetCellEFraction(Short_t cellNumber) { return -1; }
  inline Double_t GetEFraction(Short_t pos) const { return -1; }
  
  inline void     SetEFraction    (Short_t /*pos*/, Double32_t /*efrac*/) { ; }
  inline void     SetCellEFraction(Short_t /*cellNumber*/, Double32_t /*efrac*/) {;}
  
 protected:
  
  Int_t       fNCells;       ///< Number of cells
  
  /// Array of cell absolute Id. numbers.
  Short_t    *fCellNumber;   //[fNCells]
  
  /// Array with cell amplitudes (= energy!).
  Double32_t *fAmplitude;    //[fNCells][0.,0.,16]
  
  Char_t      fType;         ///< Cell type.

  Int_t       fCapacity;     //!<! Capacity of the containers
  Bool_t      fIsSorted;     //!<! True if cell arrays are sorted by index.

  Short_t    *fSwapCellNumber; //!<! Array used as a swap for sorting operations
  Double32_t *fSwapAmplitude;  //!<! Array used as a swap for sorting operations
  Int_t      *fSwapIndexArray; //!<! Array used as a swap for sorting operations
  Int_t       fSwapCapacity;   //!<! Size of the above arrays used as a swap for sorting operations

  /// \cond CLASSIMP
  ClassDef(AliHLTEMCALCaloCells, 1) ;
  /// \endcond
};

///
/// Given the position index in the array, return the cell parameters.
///
/// \param pos: Index of cell in array
/// \param cellNumber: Absolute cell Id. number
/// \param amplitude: Cell energy
/// \param time: Cell time
/// \param mclabel: MC particle index in kine tree
/// \param efrac: Fraction of energy (embedding)
///
/// \return True if pos is correct not negative or larger than expected.
///
Bool_t AliHLTEMCALCaloCells::GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude,
    Double_t & time, Int_t & mclabel, Double_t & efrac) const
{ 
  if (pos>=0 && pos<fNCells) {
    cellNumber = fCellNumber[pos];
    amplitude  = fAmplitude[pos];
    
    time       = -1;
    mclabel    = -1;
    efrac      = -1;
    return kTRUE;

  }
  else {
    return kFALSE;
  }
}

///
/// \return Cell amplitude (GeV).
/// \param cellNumber: Cell absolute Id.
///
Double_t AliHLTEMCALCaloCells::GetCellAmplitude(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted = kTRUE;
  }

  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber ) {
    return fAmplitude[pos];
  }
  else {
    return 0.;
  }
}

///
/// \return Cell amplitude (GeV).
/// \param pos: Cell position in array.
///
Double_t AliHLTEMCALCaloCells::GetAmplitude(Short_t pos) const
{ 
  if (pos>=0 && pos<fNCells) {
    return fAmplitude[pos];
  }
  else {
    return 0.;
  }
}

///
/// \return Cell absolute Id. number.
/// \param pos: Cell position in array.
///
Short_t AliHLTEMCALCaloCells::GetCellNumber(Short_t pos) const
{ 
  if (pos>=0 && pos<fNCells) {
    return fCellNumber[pos];
  }
  else {
    return fNCells;
  }
}

///
/// \param cellNumber: Cell absolute Id. number.
/// \return Cell position in array.
///
Short_t AliHLTEMCALCaloCells::GetCellPosition(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Int_t nabove, nbelow, middle;
  Short_t pos = -1;
  
  nabove = fNCells + 1;
  nbelow = 0;
  while (nabove - nbelow > 1) {
    middle = (nabove + nbelow) / 2;
    if (cellNumber == fCellNumber[middle-1]) {
      pos =   middle - 1;
      break;
    }
    if (cellNumber  < fCellNumber[middle-1]) nabove = middle;
    else                                     nbelow = middle;
  }
  
  return pos;
}

#endif
