#ifndef MUONSegResV04_H
#define MUONSegResV04_H
/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "AliMUONSegResV01.h"
#include "TArrayF.h"
#include "TArrayI.h"
class AliMUONsegmentationV04 :
public AliMUONsegmentationV01 {
 public:
    AliMUONsegmentationV04(){}
    virtual ~AliMUONsegmentationV04(){}
    // Initialisation
    virtual void Init(AliMUONchamber*);
    // Test points for auto calibration
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    ClassDef(AliMUONsegmentationV04,1)
};
	

#endif






