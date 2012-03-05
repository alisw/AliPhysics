#ifndef ALIDISPLACEDVERTEXSELECTION_H
#define ALIDISPLACEDVERTEXSELECTION_H
#include <TObject.h>
class AliESDEvent;
class AliDisplacedVertexSelection : public TObject 
{
public:
  AliDisplacedVertexSelection();
  AliDisplacedVertexSelection(const AliDisplacedVertexSelection& o);
  AliDisplacedVertexSelection& operator=(const AliDisplacedVertexSelection& o);
  
  void Output(TList* l, const char* name=0) const;
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
