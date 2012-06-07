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
  void Output(TList* l, const char* name=0) const;
  /** 
   * Print information 
   * 
   * @param option  Not used 
   */
  void Print(Option_t* option="") const;
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
  
  
  ClassDef(AliDisplacedVertexSelection,1); // Cuts on ESD Mult 
};

#endif
// Local Variables: 
//  mode: C++
// End:
