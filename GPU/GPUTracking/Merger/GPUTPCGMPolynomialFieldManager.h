//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUTPCGMPolynomialFieldManager.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCGMPOLYNOMIALFIELDMANAGER_H
#define GPUTPCGMPOLYNOMIALFIELDMANAGER_H

#include "GPUCommonDef.h"
class AliMagF;

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCGMPolynomialField;
}
} // namespace GPUCA_NAMESPACE

/**
 * @class GPUTPCGMPolynomialFieldManager
 *
 */

class GPUTPCGMPolynomialFieldManager
{
 public:
  enum StoredField_t { kUnknown,
                       kUniform,
                       k2kG,
                       k5kG }; // known fitted polynomial fields, stored in constants

  GPUTPCGMPolynomialFieldManager() CON_DEFAULT;

  /* Get appropriate pre-calculated polynomial field for the given field value nominalFieldkG
 */
  static int GetPolynomialField(float nominalFieldkG, GPUCA_NAMESPACE::gpu::GPUTPCGMPolynomialField& field);

#if defined(GPUCA_ALIROOT_LIB) & !defined(GPUCA_GPUCODE)

  /* Get pre-calculated polynomial field for the current ALICE field (if exists)
 */
  static int GetPolynomialField(GPUCA_NAMESPACE::gpu::GPUTPCGMPolynomialField& field);

  /* Fit given field for TPC
 */
  static int FitFieldTpc(AliMagF* fld, GPUCA_NAMESPACE::gpu::GPUTPCGMPolynomialField& field, double step = 1.);

  /* Fit given field for TRD
 */
  static int FitFieldTrd(AliMagF* fld, GPUCA_NAMESPACE::gpu::GPUTPCGMPolynomialField& field, double step = 1.);

  /* Fit given field for ITS
 */
  static int FitFieldIts(AliMagF* fld, GPUCA_NAMESPACE::gpu::GPUTPCGMPolynomialField& field, double step = 1.);

#endif

  /* Get pre-calculated polynomial field of type "type", scaled with respect to nominalFieldkG
 */
  static int GetPolynomialField(StoredField_t type, float nominalFieldkG, GPUCA_NAMESPACE::gpu::GPUTPCGMPolynomialField& field);
};

#endif
