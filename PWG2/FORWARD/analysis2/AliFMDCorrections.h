#ifndef ALIROOT_PWG2_FORWARD_ANALYSIS2_ALIFMDCORRECTIONS_H
#define ALIROOT_PWG2_FORWARD_ANALYSIS2_ALIFMDCORRECTIONS_H
#include <TNamed.h>
#include <TList.h>
#include "AliForwardUtil.h"
class AliESDFMD;
class TH2D;

/** 
 * This class calculates the inclusive charged particle density
 * in each for the 5 FMD rings. 
 *
 * @par Input:
 *   - AliESDFMD object possibly corrected for sharing
 *
 * @par Output:
 *   - 5 RingHistos objects - each with a number of vertex dependent 
 *     2D histograms of the inclusive charge particle density 
 * 
 * @par Corrections used: 
 *   - AliFMDAnaCalibBackgroundCorrection
 *   - AliFMDAnaCalibEventSelectionEfficiency
 *   - AliFMDAnaCalibSharingEfficiency
 *
 * @ingroup pwg2_forward_analysis 
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
  virtual Bool_t Correct(AliForwardUtil::Histos& hists, Int_t vtxBin);
  /** 
   * Scale the histograms to the total number of events 
   * 
   * @param nEvents Number of events 
   */
  void ScaleHistograms(TList* dir, Int_t nEvents);
  /** 
   * Output diagnostic histograms to directory 
   * 
   * @param dir List to write in
   */  
  void DefineOutput(TList* dir);
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
    void Output(TList* dir);
    /** 
     * Scale the histograms to the total number of events 
     * 
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
  Double_t fMultCut;       // Low cut on scaled energy loss

  ClassDef(AliFMDCorrections,1); // Calculate Nch density 
};

#endif
// Local Variables:
//   mode: C++
// End:

