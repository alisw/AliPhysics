// -*- mode: C++ -*-
/** @file 
    @brief Declaration of a 2-dimensional Flow "histogram" */
#ifndef FLOW_BINNED2D_H
#define FLOW_BINNED2D_H
#include <flow/AliFMDFlowAxis.h>
#include <TObject.h>

// Forward declaration 
class AliFMDFlowBin;
class TBrowser;

//______________________________________________________
/** @class AliFMDFlowBinned2D flow/AliFMDFlowBinned2D.h <flow/AliFMDFlowBinned2D.h>
    @brief A 2 dimensional "histogram" of objects of class AliFMDFlowBin. 
    @ingroup c_binned */
class AliFMDFlowBinned2D : public TObject
{
public:
  /** Constructor 
      @param order    Order of the harmonic 
      @param nxbins   Number of X bins.
      @param xbins    Borders of X bins (@a nxbins+1 entries) 
      @param nybins   Number of Y bins.
      @param ybins    Borders of Y bins (@a nybins+1 entries) */
  AliFMDFlowBinned2D(UShort_t order, 
		     UShort_t nxbins, Double_t* xbins,
		     UShort_t nybins, Double_t* ybins);
  /** Constructor 
      @param order Order of the harmonic 
      @param xaxis Axis object 
      @param yaxis Axis object */
  AliFMDFlowBinned2D(UShort_t order, const AliFMDFlowAxis& xaxis, 
		     const AliFMDFlowAxis& yaxis);
      
  /** Copy constructor */
  AliFMDFlowBinned2D(const AliFMDFlowBinned2D& other);
  /** Copy constructor */
  AliFMDFlowBinned2D& operator=(const AliFMDFlowBinned2D& other);
  /** Destructor */
  virtual ~AliFMDFlowBinned2D();
  
  /** Called at the beginning of an event */
  virtual void Begin();
  /** Called at the end of an event */ 
  virtual void End();
  /** Called to add a contribution to the event plane 
      @param x   Bin to fill into 
      @param y   Bin to fill into 
      @param phi The angle @f$ \varphi@f$ in radians 
      @param w   Weight
      @param a   If true, add to sub-event A, otherwise sub-event B */ 
  virtual Bool_t AddToEventPlane(Double_t x, Double_t y, Double_t phi, 
				 Double_t w, Bool_t a);
  /** Called to add a contribution to the harmonic
      @param x   Bin to fill into 
      @param y   Bin to fill into 
      @param phi The angle @f$ \varphi@f$ in radians */
  virtual Bool_t AddToHarmonic(Double_t x, Double_t y, Double_t phi);
  /** Process a full event. 
      @param phis List of @a n @f$ \varphi=[0,2\pi]@f$ angles 
      @param xs   List of @a n @f$ x@f$ values. 
      @param ys   List of @a n @f$ y@f$ values. 
      @param ws   Weights
      @param n    Size of @a phis and @a xs */
  virtual void Event(Double_t* phis, Double_t* xs, Double_t* ys, 
		     Double_t* ws, ULong_t n);
  /** Called at the end of a run */
  virtual void Finish();
  /** Get the bin at @a x 
      @param x The bin value to find a flow object for. 
      @param y The bin value to find a flow object for. 
      @return The flow object at @a x,y or 0 if out of range */ 
  virtual AliFMDFlowBin* GetBin(Double_t x, Double_t y) const;
  /** Get the bin at @a x 
      @param i The bin index to find a flow object for. 
      @param j The bin index to find a flow object for. 
      @return The flow object at @a i,j or 0 if out of range */ 
  virtual AliFMDFlowBin* GetBin(UShort_t i, UShort_t j) const;
  /** Whether this is to be considered a folder */
  Bool_t IsFolder() const { return kTRUE; }
  /** Browse this object */ 
  void Browse(TBrowser* b);
protected:
  /** X axis */ 
  AliFMDFlowAxis fXAxis;
  /** Y axis */ 
  AliFMDFlowAxis fYAxis;
  /** Array of the flow objects */ 
  AliFMDFlowBin** fBins;
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowBinned2D,1);
};


#endif
//
// EOF
//

