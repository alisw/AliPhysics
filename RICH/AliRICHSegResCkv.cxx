#include "AliRICHSegResCkv.h"
#include "TMath.h"
#include "TRandom.h"


//___________________________________________
ClassImp(AliRICHresponseCkv)
//___________________________________________	
Float_t AliRICHresponseCkv::IntPH(Float_t)
{
    
    Float_t charge = -fChslope*TMath::Log(gRandom->Rndm());
    return charge;
}

