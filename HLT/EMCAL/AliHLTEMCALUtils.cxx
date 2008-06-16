#include "AliHLTEMCALUtils.h"

#include "TClass.h"

#include "AliEMCALRecParam.h"
#include "AliEMCALReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "AliEMCALGeometry.h"
#include "AliEMCALClusterizerv1.h" 
#include "AliEMCALRawUtils.h"      

ClassImp(AliHLTEMCALUtils);

AliEMCALGeometry*      AliHLTEMCALUtils::fgGeom        = NULL;
AliEMCALClusterizerv1* AliHLTEMCALUtils::fgClusterizer = NULL;
AliEMCALRecParam*      AliHLTEMCALUtils::fgRecParam    = NULL;
AliEMCALRawUtils*      AliHLTEMCALUtils::fgRawUtils    = NULL;

AliHLTEMCALUtils::AliHLTEMCALUtils()
  : TObject()
{
  //
  // Default constructor
  //
  ;
}

AliHLTEMCALUtils::~AliHLTEMCALUtils()
{
  //
  // Destructor
  //
  ;
}

AliHLTEMCALUtils::AliHLTEMCALUtils(const AliHLTEMCALUtils & /*t*/)  
  : TObject()
{
  //
  // copy ctor not to be used
  //
  AliFatal("May not use.");  
}
  
AliHLTEMCALUtils& AliHLTEMCALUtils::operator = (const AliHLTEMCALUtils & /*t*/)  
{
  //
  // assignement operator not to be used
  //
  AliFatal("May not use.") ;
  return *this ; 
}

void AliHLTEMCALUtils::InitRecParam()
{
  //
  // Please check the AliEMCALReconstructor for comparison
  // Check if the instance of AliEMCALRecParam exists, 
  // if not, get it from OCDB if available, otherwise create a default one
  //

  fgRecParam = (AliEMCALRecParam*) AliEMCALReconstructor::GetRecParam();
  if (fgRecParam)
    return;
  
 if (!fgRecParam  && (AliCDBManager::Instance()->IsDefaultStorageSet())) 
   {
     AliCDBEntry *entry = (AliCDBEntry*) 
       AliCDBManager::Instance()->Get("EMCAL/Config/RecParam");
     if (entry) fgRecParam = (AliEMCALRecParam*) entry->GetObject();
   }
 
  if(!fgRecParam)
    {
      AliWarningClass("The Reconstruction parameters initialized to default.");
      fgRecParam = new AliEMCALRecParam;
    }

  if (!fgRecParam)
    {
      AliErrorClass("Unable to init the reco params. Something is really wrong. Memory?");
    }
  
  AliEMCALReconstructor::SetRecParam(fgRecParam);
}

const AliEMCALRecParam* AliHLTEMCALUtils::GetRecParam()
{
  //
  // Init the parameters and reuse the Reconstructor
  //
  AliHLTEMCALUtils::InitRecParam();
  return fgRecParam;
}

const AliEMCALRawUtils*      AliHLTEMCALUtils::GetRawUtils()
{
  //
  // Init EMCAL raw utils
  //
  if (fgRawUtils == NULL)
    {
      if (AliCDBManager::Instance()->IsDefaultStorageSet())
	{
	  fgRawUtils = new AliEMCALRawUtils;
	  AliInfoClass("Raw Utils initialized.");
	}
      else
	{
	  AliErrorClass("OCDB not initialized. Unable to init raw utils.");
	}
    }

  return fgRawUtils;
}

const AliEMCALClusterizerv1* AliHLTEMCALUtils::GetClusterizer()
{
  //
  // Init EMCAL clusterizer
  //
  if (fgClusterizer == NULL)
    {
      fgClusterizer = new AliEMCALClusterizerv1;
      AliInfoClass("ClusterizerV1 initialized.");
   }

  return fgClusterizer;
}

const AliEMCALGeometry*      AliHLTEMCALUtils::GetGeometry()
{
  //
  // Init EMCAL geometry
  //
  if (fgGeom == NULL)
    {
      AliInfoClass(Form("Using default geometry"));
      fgGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    }
  return fgGeom;
}
