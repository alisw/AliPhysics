#include <TGClient.h>
#include <TGPicture.h>
#include <TSystem.h>
#include <TString.h>

#include "AliEveUtil.h"

TGPicturePool* gAliEvePicturePool=0;
AliEveUtil* AliEveUtil::fgAliEveUtil=0;

ClassImp(AliEveUtil)
AliEveUtil::AliEveUtil()
{
    gAliEvePicturePool = GetPicturePool();
}

void AliEveUtil::Init()
{
    if(fgAliEveUtil) return;

    fgAliEveUtil = new AliEveUtil;
}


TGPicturePool* AliEveUtil::GetPicturePool()
{
    if(gAliEvePicturePool) return gAliEvePicturePool;

    TString iconSearchPath(gSystem->Getenv("ALICE_ROOT") );
    iconSearchPath.Append("/EVE/icons/");

    gAliEvePicturePool = new TGPicturePool(gClient, iconSearchPath.Data() );

    return gAliEvePicturePool;
}

