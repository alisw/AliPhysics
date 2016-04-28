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
  void            Sort() {;}
  
  Bool_t          IsEMCAL()  const { return (fType == kEMCALCell); }
  Bool_t          IsPHOS()   const { return (fType == kPHOSCell) ; }
  Char_t          GetType()  const { return  fType               ; }
  void            SetType(Char_t t){ fType = t                   ; }
  
  inline Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, Double_t &time, Int_t &mclabel,      Double_t &efrac) const;
  Bool_t          SetCell(Short_t pos, Short_t  cellNumber, Double_t  amplitude, Double_t  time, Int_t  mclabel = -1, Double_t  efrac = 0., Bool_t isHG=kFALSE);
  Bool_t          AddCell(Short_t pos, Double_t  amplitude);

  Short_t         GetNumberOfCells() const  { return fNCells ; }
  void            SetNumberOfCells(Int_t n);
  
  inline Double_t GetCellAmplitude(Short_t cellNumber) { return GetAmplitude(cellNumber); }
  inline Bool_t   GetCellHighGain(Short_t cellNumber) { return kFALSE; }  //is this cell High Gain
  inline Short_t  GetCellPosition(Short_t cellNumber) { return cellNumber; }
  inline Double_t GetCellTime(Short_t cellNumber) { return -1; }
  
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Bool_t   GetHighGain(Short_t pos) const { return kFALSE; }
  inline Double_t GetTime(Short_t pos) const { return -1; }
  inline Short_t  GetCellNumber(Short_t pos) const { return pos; }

  // MC & embedding
  inline Int_t    GetCellMCLabel(Short_t cellNumber) { return -1; }
  inline Int_t    GetMCLabel(Short_t pos) const { return -1; }
  
  inline Double_t GetCellEFraction(Short_t cellNumber) { return -1; }
  inline Double_t GetEFraction(Short_t pos) const { return -1; }
  
  inline void     SetEFraction    (Short_t /*pos*/, Double_t /*efrac*/) { ; }
  inline void     SetCellEFraction(Short_t /*cellNumber*/, Double_t /*efrac*/) {;}
  
 protected:
  
  Int_t       fCapacity;     ///< Total number of cells
  Int_t       fNCells;       ///< Number of cells
  
  /// Array with cell amplitudes (= energy!).
  Double_t   *fAmplitude;    //[fCapacity][0.,0.,16]
  
  Char_t      fType;         ///< Cell type.

  /// \cond CLASSIMP
  ClassDef(AliHLTEMCALCaloCells, 2) ;
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
  if (pos >= 0 && pos < fCapacity) {
    cellNumber = pos;
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
/// \param pos: Cell position in array.
///
Double_t AliHLTEMCALCaloCells::GetAmplitude(Short_t pos) const
{ 
  if (pos >= 0 && pos < fNCells) {
    return fAmplitude[pos];
  }
  else {
    return 0.;
  }
}

#endif
