#ifndef ALIDISPLACEDVERTEXSELECTION_H
#define ALIDISPLACEDVERTEXSELECTION_H
#include <TObject.h>
class AliESDEvent;

/** 
 * Selection of events from satellite interactions 
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

  Double_t GetVertexZ() const { return fVertexZ; }
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
  Double_t fVertexZ;
  Double_t fCent;
  
  ClassDef(AliDisplacedVertexSelection,2); // Cuts on ESD Mult 
};

#endif
// Local Variables: 
//  mode: C++
// End:
