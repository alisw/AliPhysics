// Author: Svein Lindal <slindal@fys.uio.no>

#ifndef ALIHLTEVEANY_H
#define ALIHLTEVEANY_H

#include "AliHLTEveBase.h"
class AliHLTHOMERBlockDesc;

class AliHLTEveAny : public AliHLTEveBase {

public:
  
  /** Constructor  **/
  AliHLTEveAny();

  /** Destructor **/
 ~AliHLTEveAny();

  /** Inherited form AliHLTEveBase */
  void ProcessBlock(AliHLTHOMERBlockDesc * block);

  /** inherited from AliHLTEveBase */
  void UpdateElements();
  
  /** inherited from AliHLTEveBase */
  void ResetElements();

private:
  
  /** copy constructor prohibited */
  AliHLTEveAny(const AliHLTEveAny&);
  /** assignment operator prohibited */
  AliHLTEveAny& operator = (const AliHLTEveAny &);

  void ProcessHistogram(AliHLTHOMERBlockDesc * block );

  ClassDef(AliHLTEveAny, 0);
};

#endif
