/**
 * @file   AliDisplacedVertexSelection.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 00:53:03 2013
 * 
 * @brief  Selection of events from satellite interactions 
 * 
 * @ingroup pwglf_forward_aod
 * 
 */
#ifndef ALIDISPLACEDVERTEXSELECTION_H
#define ALIDISPLACEDVERTEXSELECTION_H
#include <TObject.h>
class AliESDEvent;

/** 
 * Selection of events from satellite interactions 
 *
 * @ingroup pwglf_forward_aod
 */
class AliDisplacedVertexSelection : public TObject 
{
public:
  /** 
   * Constructor 
   */
  AliDisplacedVertexSelection();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliDisplacedVertexSelection(const AliDisplacedVertexSelection& o);
  /** 
   * Assingment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliDisplacedVertexSelection& operator=(const AliDisplacedVertexSelection& o);
  /** 
   * Define the output 
   * 
   * @param l     List to add output to
   * @param name  Name of the list 
   */
  void CreateOutputObjects(TList* l, const char* name=0) const;
  /** 
   * Print information 
   * 
   * @param option  Not used 
   */
  void Print(Option_t* option="") const;
  /** 
   * Process an ESD event to get the information 
   * 
   * @param esd ESD event 
   * 
   * @return true on success
   */
  Bool_t Process(const AliESDEvent* esd);
  /** 
   * Get the interaction point Z-coordinate from ZDC timing. 
   * 
   * 
   * @return Interaction point Z-coordinate
   */
  Double_t GetVertexZ() const { return fVertexZ; }
  /** 
   * Return the centrality percentile 
   * 
   * 
   * @return Centrality percentile (ZDC vs ZEM)
   */
  Double_t GetCentralityPercentile() const { return fCent; }
  /** 
   * Check for displaced vertices (M.Guilbaud) 
   * 
   * @param esd  Event 
   * 
   * @return displaced vertex
   */
  Double_t CheckDisplacedVertex(const AliESDEvent* esd) const;
   /** 
   * Calculate Centrality for displaced vertices (M.Guilbaud) 
   * 
   * @param esd  Event 
   * 
   * @return displaced vertex centrality
   */
  Double_t CalculateDisplacedVertexCent(const AliESDEvent* esd) const;
  
protected:
  Double_t fVertexZ; // Interaction point Z-coordinate
  Double_t fCent;    // Centrality percentile
  
  ClassDef(AliDisplacedVertexSelection,2); // Satelitte collisions 
};

#endif
// Local Variables: 
//  mode: C++
// End:
