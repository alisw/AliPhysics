#include <stdio.h>
#ifdef __APPLE__

#include "AliTPCTrackHitsInterfaces.h"


AliClassAliTrackHitsInfo & g_AliClassAliTrackHitsInfo()
{
  static AliClassAliTrackHitsInfo p;
  //  printf("AliClassAliTrackHitsInfo\n");
  return p;
}

AliClassAliTrackHitsParam & g_AliTrackHitsParam()
{
  static AliClassAliTrackHitsParam p;
  //  printf("AliClassAliTrackHitsParam\n");
  return p;
}

AliClassAliHitInfo & g_AliHitInfo()
{
  static AliClassAliHitInfo p;
  //  printf("AliClassAliHitInfo\n");
  return p;
}
#endif
