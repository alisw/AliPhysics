#ifndef RICHSegResCkv_H
#define RICHSegResCkv_H

#include "AliRICHSegResV0.h"

class AliRICHresponseCkv : public AliRICHresponseV0 {
    
 public:
    AliRICHresponseCkv(){}
    virtual ~AliRICHresponseCkv(){}
    
    virtual Float_t IntPH();
    
    ClassDef(AliRICHresponseCkv,1)
	
	
	};
	
#endif
