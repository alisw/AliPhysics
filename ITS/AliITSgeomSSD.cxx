#include "AliITSgeomSSD.h"

ClassImp(AliITSgeomSSD)
AliITSgeomSSD::AliITSgeomSSD(){
    //
    // default constructor
    //
    fShapeSSD = new TBRIK("ActiveSSD","Active volume of SSD","SI",
		  	    3.0E-2/2.,1.73/2.,4.0/2.);
}
