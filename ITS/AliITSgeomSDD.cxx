#include "AliITSgeomSDD.h"

ClassImp(AliITSgeomSDD)
AliITSgeomSDD::AliITSgeomSDD(){
    //
    // default constructor
    //
    fShapeSDD = new TBRIK("ActiveSDD","Active volume of SDD","SDD SI CHIP",
			    3.0E-2/2.,7.25/2.,7.53/2.);
}
