/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
/// \class AliAODCaloCells
/// \brief Class for calorimeter cell AOD data handling
///
/// AOD class to store calorimeter cell data
///
/// Data is stored in different arrays, each entry of the array corresponding to a cell.
/// The data stored is the cell energy, time, high gain bool, absolute id number, 
/// MC label that deposited most energy, and a container for the MC embedded energy.
///
/// \author Markus Oldenburg, CERN.
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-Grenoble
///
//-------------------------------------------------------------------------

#ifndef ALIAODCELLS_H
#define ALIAODCELLS_H

#include <AliVCaloCells.h>
#include <TMath.h>

class AliAODCaloCells : public AliVCaloCells
{
 public:
  AliAODCaloCells();
  AliAODCaloCells(const char* name, const char* title, VCells_t ttype=kUndef);
  AliAODCaloCells(const AliAODCaloCells& cells); 
  AliAODCaloCells& operator=(const AliAODCaloCells& cells);
  virtual ~AliAODCaloCells();

  virtual AliVCaloCells* CopyCaloCells(Bool_t all) const;
  virtual void    Copy(TObject &obj) const;
  void            Clear(const Option_t*);
  void            CreateContainer(Short_t nCells);
  void            DeleteContainer();
  void            Sort();
  
  inline Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude,  Double_t &time, Int_t &mclabel,      Double_t &efrac) const ;
  Bool_t          SetCell(Short_t pos, Short_t  cellNumber, Double_t  amplitude, Double_t  time, Int_t  mclabel = -1, Double_t  efrac = 0., Bool_t isHG=kFALSE);
  
  Short_t         GetNumberOfCells() const  { return fNCells ; }
  void            SetNumberOfCells(Int_t n) { fNCells = n    ; }
  
  inline Double_t GetCellAmplitude(Short_t cellNumber);
  inline Bool_t   GetCellHighGain(Short_t cellNumber);  //is this cell High Gain
  inline Short_t  GetCellPosition(Short_t cellNumber);
  inline Double_t GetCellTime(Short_t cellNumber);
  
  inline Double_t GetAmplitude(Short_t pos) const;
  inline Bool_t   GetHighGain(Short_t pos) const;
  inline Short_t  GetCellNumber(Short_t pos) const;
  inline Double_t GetTime(Short_t pos) const;
  
  Bool_t          IsEMCAL() const { return (fType == kEMCALCell); }
  Bool_t          IsPHOS()  const { return (fType == kPHOSCell) ; }
  
  Char_t          GetType() const { return fType;}
  void            SetType(Char_t ttype) { fType=ttype; }
  
  // MC & embedding
  inline Int_t    GetCellMCLabel(Short_t cellNumber) ;
  inline Int_t    GetMCLabel(Short_t pos) const ;
  
  inline Double_t GetCellEFraction(Short_t cellNumber) ;
  inline Double_t GetEFraction(Short_t pos) const ;  
  
  inline void     SetEFraction    (Short_t pos,         Double32_t efrac) ;
  inline void     SetCellEFraction(Short_t cellNumber,  Double32_t efrac) ;
  
 protected:
  
  Int_t       fNCells;       ///< Number of cells.
  
  /// If Cell is High Gain or Low Gain
  Bool_t     *fHGLG;         //[fNCells] 
  
  /// Array of cell absolute Id. numbers.
  Short_t    *fCellNumber;   //[fNCells] 
  
  /// Array with cell amplitudes (= energy!).
  Double32_t *fAmplitude;    //[fNCells][0.,0.,16] 
  
  /// Array with cell times.
  Double32_t *fTime;         //[fNCells][0.,0.,16] 
  
  /// Array with fraction of MC energy and data - for embedding.
  Double32_t *fEFraction;    //[fNCells][0.,0.,16] 
  
  ///  Array of MC labels, each label is the highest contributor to the cell signal.
  Int_t      *fMCLabel;      //[fNCells]
  
  Bool_t      fIsSorted;     //!<! True if cell arrays are sorted by index.
  
  Char_t      fType;         ///< Cell type
  
  /// \cond CLASSIMP
  ClassDef(AliAODCaloCells, 5) ;
  /// \endcond

};

///
/// Given the position index in the array, return the cell parameters
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
Bool_t AliAODCaloCells::GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, 
                                Double_t &time, Int_t & mclabel, Double_t & efrac) const 
{ 
  if (pos>=0 && pos<fNCells) 
  {
    cellNumber = fCellNumber[pos];
    amplitude  = fAmplitude[pos];
    
    if(fTime)      time    = fTime[pos];
    else           time    =-1.;
    if(fMCLabel)   mclabel = fMCLabel[pos];
    else           mclabel =-1 ; 
    if(fEFraction) efrac   = fEFraction[pos];
    else           efrac   = 0 ;
    
    return kTRUE;
    
  } else {
    return kFALSE;
  }
}

///
/// \return Cell amplitude (GeV)
/// \param cellNumber: Cell absolute Id.
///
Double_t AliAODCaloCells::GetCellAmplitude(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }

  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && fCellNumber[pos] == cellNumber) {
    return fAmplitude[pos];
  } else {
    return 0.;
  }
}

///
/// \return True for high gain, False for low gain.
/// \param cellNumber: Cell absolute Id.
///
Bool_t AliAODCaloCells::GetCellHighGain(Short_t cellNumber)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }

  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber ) {
    if(fHGLG)
      return fHGLG[pos];
    else{
      if(fMCLabel) //old version of AOD, 
         return !(fMCLabel[pos]==-2) ;
      else
	 return kTRUE ;
    }
  } else {
    return kFALSE;
  }
}

///
/// \return Cell time (s).
/// \param cellNumber: Cell absolute Id.
///
Double_t AliAODCaloCells::GetCellTime(Short_t cellNumber)
{ 
  if(!fTime) return -1;
  
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber) {
    return fTime[pos];
  } else {
    return -1.;
  }
}

///
/// \return Cell amplitude (GeV).
/// \param pos: Cell position in array.
///
Double_t AliAODCaloCells::GetAmplitude(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fAmplitude[pos];
  } else {
    return 0.;
  }
}

///
/// \return True for high gain, False for low gain.
/// \param pos: Cell position in array.
///
Bool_t AliAODCaloCells::GetHighGain(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    if(fHGLG)
      return fHGLG[pos];
    else{
      if(fMCLabel)    //Old version of AOD store this flag in MCLabel
        return !(fMCLabel[pos]==-2) ;
      else
	return kTRUE ;
    }
  } else {
    return kFALSE;
  }
}

///
/// \return Cell time (s).
/// \param pos: Cell position in array.
///
Double_t AliAODCaloCells::GetTime(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells && fTime) {
    return fTime[pos];
  } else {
    return -1.;
  }
}

///
/// \return Cell absolute Id. number.
/// \param pos: Cell position in array.
///
Short_t AliAODCaloCells::GetCellNumber(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells) {
    return fCellNumber[pos];
  } else {
    return fNCells;
  }
}

///
/// \param cellNumber: Cell absolute Id. number.
/// \return Cell position in array.
///
Short_t AliAODCaloCells::GetCellPosition(Short_t cellNumber)
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

///
/// \return MC label of highest contributor particle depositing energy in cell.
/// \param pos: Cell position in array.
///
Int_t AliAODCaloCells::GetMCLabel(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells && fMCLabel) {
    return fMCLabel[pos];
  } else {
    return 0;
  }
}

///
/// \return Fraction of energy in cell from embedding.
/// \param pos: Cell position in array.
///
Double_t AliAODCaloCells::GetEFraction(Short_t pos) const 
{ 
  if (pos>=0 && pos<fNCells && fEFraction) {
    return fEFraction[pos];
  } else {
    return 0.;
  }
}

///
/// \return MC label of highest contributor particle depositing energy in cell.
/// \param cellNumber: Cell position in array.
///
Int_t AliAODCaloCells::GetCellMCLabel(Short_t cellNumber)
{ 
  if(!fMCLabel) return -1;
  
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && fCellNumber[pos] == cellNumber) {
    return fMCLabel[pos];
  } else {
    return 0;
  }
}

///
/// \return Fraction of energy in cell from embedding.
/// \param cellNumber: Absolute Id number of cell.
///
Double_t AliAODCaloCells::GetCellEFraction(Short_t cellNumber)
{ 
  if(!fEFraction) return 0;

  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber) {
    return fEFraction[pos];
  } else {
    return -1.;
  }
}

///
/// Set Fraction of energy in cell from embedding.
/// \param pos: Cell position in array.
/// \param efrac: fraction of energy.
///
void AliAODCaloCells::SetEFraction(Short_t pos,  Double32_t efrac)
{  
  if (pos>=0 && pos < fNCells) 
  {
    if(!fEFraction) fEFraction = new Double32_t[fNCells];
    fEFraction[pos]  = efrac;
  } 
}

///
/// Set fraction of energy in cell from embedding.
/// \param cellNumber: Absolute cell Id. number of cell.
/// \param efrac: fraction of energy.
///
void AliAODCaloCells::SetCellEFraction(Short_t cellNumber, Double32_t efrac)
{ 
  if (!fIsSorted) {
    Sort();
    fIsSorted=kTRUE;
  }
  
  Short_t pos = TMath::BinarySearch(fNCells, fCellNumber, cellNumber);
  if (pos>=0 && pos < fNCells && fCellNumber[pos] == cellNumber) 
  {
    if(!fEFraction) fEFraction = new Double32_t[fNCells];
    fEFraction[pos] = efrac;
  } 
}


#endif
