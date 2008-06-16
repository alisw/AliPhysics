#ifndef ALIHLTEMCALUTILS_H
#define ALIHLTEMCALUTILS_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTEMCALUtils.h
    @author m. ploskon
    @date   
    @brief  Utils for EMCAL-HLT
*/

#include <TObject.h>

/**
 * @class AliHLTEMCALUtils
 */

class AliEMCALGeometry; 
class AliEMCALClusterizerv1;
class AliEMCALRecParam;    
class AliEMCALRawUtils;

class AliHLTEMCALUtils : public TObject
{
 public:
  AliHLTEMCALUtils();
  AliHLTEMCALUtils(const AliHLTEMCALUtils & /*t*/);
  AliHLTEMCALUtils & operator = (const AliHLTEMCALUtils & /*t*/);
  virtual ~AliHLTEMCALUtils();

  static void                         InitRecParam();
  static const AliEMCALRecParam*      GetRecParam();

  static const AliEMCALRawUtils*      GetRawUtils();
  static const AliEMCALClusterizerv1* GetClusterizer();
  static const AliEMCALGeometry*      GetGeometry();

 protected:
 private:

  static AliEMCALGeometry*      fgGeom;        // EMCAL geometry
  static AliEMCALClusterizerv1* fgClusterizer; // ECMAL clusterizer
  static AliEMCALRecParam*      fgRecParam;    // EMCAL reconstruction parameters
  static AliEMCALRawUtils*      fgRawUtils;    // EMCAL raw utilities 
  ClassDef(AliHLTEMCALUtils, 0)
};

#endif
