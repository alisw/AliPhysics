#ifndef ALIHBTWEIGHTS_H
#define ALIHBTWEIGHTS_H

#include <TObject.h>

class AliHBTPair;

class AliHBTWeights: public TObject
 {
   public:
      virtual ~AliHBTWeights();
      static Double_t Weight(AliHBTPair* partpair){return (fgWeights)?fgWeights->GetWeight(partpair):0.0;}
      virtual Double_t GetWeight(AliHBTPair* partpair) = 0;
      virtual void Set() = 0;
      static AliHBTWeights* Instance() {return fgWeights;}
      
   protected:
      static AliHBTWeights* fgWeights;
      
      ClassDef(AliHBTWeights,2)
 };

#endif
