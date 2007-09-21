// -*- mode: C++ -*-
/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
/** @file 
    @brief Declaration of an Axis in a Flow "histogram" */
#ifndef ALIFMDFLOWAXIS_H
#define ALIFMDFLOWAXIS_H
#ifndef ROOT_TObject
# include <TObject.h>
#endif
//  Axis object for the 
//  AliFMDFlowBinned1D and 
//  AliFMDFlowBinned2D 
//  "histograms" of  objects 
//  of class AliFMDFlowBin. 
//

//______________________________________________________
/** @class AliFMDFlowAxis flow/AliFMDFlowAxis.h <flow/AliFMDFlowAxis.h>
    @brief Axis object for the AliFMDFlowBinned1D and
    AliFMDFlowBinned2D "histograms" of  objects of class AliFMDFlowBin. 
    @ingroup c_binned 
*/
class AliFMDFlowAxis : public TObject
{
public:
  /** Constructor 
      @param n    Number of bins 
      @param bins Bin limits (@a n+1 entries) */
  AliFMDFlowAxis(UShort_t n, Double_t* bins);
  /** Constructor 
      @param n    Number of bins 
      @param min  Minimum 
      @param max  Maximum */
  AliFMDFlowAxis(UShort_t n, Double_t min, Double_t max);
  /** Copy constructor 
      @param other Object to copy from */
  AliFMDFlowAxis(const AliFMDFlowAxis& other);
  /** Assignement operator 
      @param other Object to assign from */ 
  AliFMDFlowAxis& operator=(const AliFMDFlowAxis& other);
  /** Destructor */ 
  virtual ~AliFMDFlowAxis();
  /** Find a bin 
      @param v Value to look for 
      @return bin number of @a v */
  Int_t FindBin(Double_t v) const;
  /** Get the width of the @a i th bin 
      @param i  Bin to get width of */
  Double_t BinWidth(UShort_t i) const;
  /** Get the center of a bin 
      @param i Bin to get center of 
      @return Center of the @a i th bin */
  Double_t BinCenter(UShort_t i) const;
  /** Get the lower limit of a bin 
      @param i Bin to get the lower limit of 
      @return lower limit of bin @a i */
  Double_t BinLower(UShort_t i) const;
  /** Get the upper limit of a bin 
      @param i Bin to get the upper limit of 
      @return upper limit of bin @a i */
  Double_t BinUpper(UShort_t i) const;
  /** Get a pointer to the bins 
      @return pointer to the bins */ 
  Double_t* Bins() const { return fBins; }
  /** Get the number of bins */ 
  UShort_t N() const { return fN; }
  /** Print the axis */ 
  void Print(Option_t* option="") const; //*MENU*
protected:
  /** Number of bins */ 
  UShort_t fN; // Number of bins
  /** Borders of the bins */  
  Double_t* fBins; //[fN+1] Bin boundaries
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowAxis,1);
};  


#endif
//
// EOF
//
