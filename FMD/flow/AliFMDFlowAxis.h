// -*- mode: C++ -*-
/** @file 
    @brief Declaration of an Axis in a Flow "histogram" */
#ifndef FLOW_AXIS_H
#define FLOW_AXIS_H
#ifndef ROOT_TObject
# include <TObject.h>
#endif

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
protected:
  /** Number of bins */ 
  UShort_t fN;
  /** Borders of the bins */ 
  Double_t* fBins; //[fN]
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowAxis,1);
};  


#endif
//
// EOF
//
