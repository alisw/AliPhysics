// -*- mode: C++ -*-
/** @file 
    @brief Declaration of a 1-dimensional Flow "histogram" */
#ifndef FLOW_BINNED1D_H
#define FLOW_BINNED1D_H
#include <flow/AliFMDFlowAxis.h>

// Forward declaration 
class AliFMDFlowBin;
class TBrowser;

//______________________________________________________
/** @class AliFMDFlowBinned1D flow/AliFMDFlowBinned1D.h <flow/AliFMDFlowBinned1D.h>
    @brief A 1 dimensional "histogram" of objects of class AliFMDFlowBin. 
    @ingroup c_binned 
    @example test_flow.cxx 
    @example ana_flow.cxx 
*/
class AliFMDFlowBinned1D : public TObject
{
public:
  /** Flags for print */
  enum {
    /** Show details */ 
    details = 0x1, 
    /** Show summary */ 
    summary = 0x2
  };
  /** Constructor 
      @param order    Order of the harmonic
      @param nxbins   Number of X bins.
      @param xbins    Borders of X bins (@a nxbins+1 entries) */
  AliFMDFlowBinned1D(UShort_t order, UShort_t nxbins, Double_t* xbins);
  /** Constructor 
      @param order    Order of the harmonic
      @param xaxis Axis object  */
  AliFMDFlowBinned1D(UShort_t order, const AliFMDFlowAxis& xaxis);
  /** Copy constructor */
  AliFMDFlowBinned1D(const AliFMDFlowBinned1D& other);
  /** Copy constructor */
  AliFMDFlowBinned1D& operator=(const AliFMDFlowBinned1D& other);
  /** Destructor */
  virtual ~AliFMDFlowBinned1D();
  
  /** Called at the beginning of an event */
  virtual void Begin();
  /** Called at the end of an event */ 
  virtual void End();
  /** Called to add a contribution to the event plane 
      @param x   Bin to fill into 
      @param w   Weight
      @param phi The angle @f$ \varphi@f$ in radians 
      @param a   If true, add to sub-event A, otherwise sub-event B */ 
  virtual Bool_t AddToEventPlane(Double_t x, Double_t phi, 
				 Double_t w, Bool_t a);
  /** Called to add a contribution to the harmonic
      @param x   Bin to fill into 
      @param phi The angle @f$ \varphi@f$ in radians */
  virtual Bool_t AddToHarmonic(Double_t x, Double_t phi);
  /** Process a full event. 
      @param phis List of @a n @f$ \varphi=[0,2\pi]@f$ angles 
      @param xs   List of @a n @f$ x@f$ values. 
      @param ws   Weights
      @param n    Size of @a phis and @a xs */
  virtual void Event(Double_t* phis, Double_t* xs, Double_t* ws, 
		     ULong_t n);
  /** Called at the end of a run */
  virtual void Finish();
  /** Get the bin at @a x 
      @param x The bin value to find a flow object for. 
      @return The flow object at @a x or 0 if out of range */ 
  virtual AliFMDFlowBin* GetBin(Double_t x) const;
  /** Get the bin @a i 
      @param i The bin number to get
      @return The flow object in bin @a i or 0 if out of range */ 
  virtual AliFMDFlowBin* GetBin(UShort_t i) const;
  /** Print to standard out */ 
  virtual void Print(Option_t* option="s") const;
  /** Whether this is to be considered a folder */
  Bool_t IsFolder() const { return kTRUE; }
  /** Browse this object */ 
  void Browse(TBrowser* b);
protected:
  /** X axis */ 
  AliFMDFlowAxis fXAxis;
  /** Array of the flow objects */ 
  AliFMDFlowBin** fBins;
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowBinned1D,1);
};

#endif
//
// EOF
//
