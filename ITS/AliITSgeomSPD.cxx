#include "AliITSgeomSPD.h"

ClassImp(AliITSgeomSPD)
AliITSgeomSPD::AliITSgeomSPD(){
    //
    // default constructor
    //
    fShapeSPD = new TBRIK("ActiveSPD","Active volume of SPD","SPD SI CHIP",
			  2.5E-2/2.0,1.38/2.0,8.2/2.0);
}
