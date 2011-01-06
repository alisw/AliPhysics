// This class calculates the exclusive charged particle density
// in each for the 5 FMD rings. 
//
#ifndef ALIFMDCORRECTIONS_H
#define ALIFMDCORRECTIONS_H
#include <TNamed.h>
#include <TList.h>
#include "AliForwardUtil.h"
class TH2D;

/** 
 * @defgroup pwg2_forward_algo Algorithms 
 *
 * @ingroup pwg2_forward 
 */
/** 
 * This class calculates the exclusive charged particle density
 * in each for the 5 FMD rings. 
 *
 * @par Input:
 *   - 5 RingHistos objects - each with a number of vertex dependent 
 *     2D histograms of the inclusive charge particle density 
 *
 * @par Output:
 *   - 5 RingHistos objects - each with a number of vertex dependent 
 *     2D histograms of the exclusive charge particle density 
 * 
 * @par Corrections used: 
 *   - AliFMDCorrSecondaryMap;
 *   - AliFMDCorrVertexBias
 *   - AliFMDCorrMergingEfficiency
 *
 * @ingroup pwg2_forward_algo 
 */
class AliFMDCorrections : public TNamed
{
public:
  /** 
   * Constructor 
   */
  AliFMDCorrections();
  /** 
   * Constructor 
   * 
   * @param name Name of object
   */
  AliFMDCorrections(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDCorrections(const AliFMDCorrections& o);
  /** 
   * Destructor 
   */
  virtual ~AliFMDCorrections();
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  AliFMDCorrections& operator=(const AliFMDCorrections&);
  /** 
   * Do the calculations 
   * 
   * @param hists    Cache of histograms 
   * @param vtxBin   Vertex bin 
   * 
   * @return true on successs 
   */
  virtual Bool_t Correct(AliForwardUtil::Histos& hists, UShort_t vtxBin);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param dir     Where the output is stored
   * @param nEvents Number of events 
   */
  virtual void ScaleHistograms(TList* dir, Int_t nEvents);
  /** 
   * Output diagnostic histograms to directory 
   * 
   * @param dir List to write in
   */  
  virtual void DefineOutput(TList* dir);
  /** 
   * Set the debug level.  The higher the value the more output 
   * 
   * @param dbg Debug level 
   */
  void SetDebug(Int_t dbg=1) { fDebug = dbg; }
  /** 
   * Print information
   * 
   * @param option Not used 
   */
  void Print(Option_t* option="") const;
protected:
  /** 
   * Internal data structure to keep track of the histograms
   */
  struct RingHistos : public AliForwardUtil::RingHistos 
  { 
    /** 
     * Default CTOR
     */
    RingHistos();
    /** 
     * Constructor
     * 
     * @param d detector
     * @param r ring 
     */
    RingHistos(UShort_t d, Char_t r);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    RingHistos(const RingHistos& o);
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this 
     */
    RingHistos& operator=(const RingHistos& o);
    /** 
     * Destructor 
     */
    ~RingHistos();
    /** 
     * Make output 
     * 
     * @param dir Where to put it 
     */
    void Output(TList* dir);
    /** 
     * Scale the histograms to the total number of events 
     * 
     * @param dir     where the output is stored
     * @param nEvents Number of events 
     */
    void ScaleHistograms(TList* dir, Int_t nEvents);
    TH2D*     fDensity;      // Distribution primary Nch
    ClassDef(RingHistos,1);
  };
  /** 
   * Get the ring histogram container 
   * 
   * @param d Detector
   * @param r Ring 
   * 
   * @return Ring histogram container 
   */
  RingHistos* GetRingHistos(UShort_t d, Char_t r) const;

  TList    fRingHistos;    // List of histogram containers
  Int_t    fDebug;         //  Debug level 

  ClassDef(AliFMDCorrections,1); // Calculate Nch density 
};

#endif
// Local Variables:
//   mode: C++
// End:

