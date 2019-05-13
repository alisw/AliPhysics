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

/// \file GPUTPCGMPolynomialField.cxx
/// \author Sergey Gorbunov, David Rohr

#include "GPUTPCGMPolynomialField.h"
using namespace GPUCA_NAMESPACE::gpu;

#if defined(GPUCA_ALIROOT_LIB) & !defined(GPUCA_GPUCODE)

#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

void GPUTPCGMPolynomialField::Print() const
{
  const double kCLight = 0.000299792458;
  typedef std::numeric_limits<float> flt;
  cout << std::scientific;
#if __cplusplus >= 201103L
  cout << std::setprecision(flt::max_digits10 + 2);
#endif
  cout << " nominal field " << mNominalBz << " [kG * (2.99792458E-4 GeV/c/kG/cm)]"
       << " == " << mNominalBz / kCLight << " [kG]" << endl;

  cout << " TpcBx[NTPCM] = { ";
  for (int i = 0; i < NTPCM; i++) {
    cout << mTpcBx[i];
    if (i < NTPCM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TpcBy[NTPCM] = { ";
  for (int i = 0; i < NTPCM; i++) {
    cout << mTpcBy[i];
    if (i < NTPCM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TpcBz[NTPCM] = { ";
  for (int i = 0; i < NTPCM; i++) {
    cout << mTpcBz[i];
    if (i < NTPCM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << "TRD field: \n"
       << endl;

  cout << " TrdBx[NTRDM] = { ";
  for (int i = 0; i < NTRDM; i++) {
    cout << mTrdBx[i];
    if (i < NTRDM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TrdBy[NTRDM] = { ";
  for (int i = 0; i < NTRDM; i++) {
    cout << mTrdBy[i];
    if (i < NTRDM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TrdBz[NTRDM] = { ";
  for (int i = 0; i < NTRDM; i++) {
    cout << mTrdBz[i];
    if (i < NTRDM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << "ITS field: \n"
       << endl;

  cout << " ItsBx[NITSM] = { ";
  for (int i = 0; i < NITSM; i++) {
    cout << mItsBx[i];
    if (i < NITSM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " ItsBy[NITSM] = { ";
  for (int i = 0; i < NITSM; i++) {
    cout << mItsBy[i];
    if (i < NITSM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " ItsBz[NITSM] = { ";
  for (int i = 0; i < NITSM; i++) {
    cout << mItsBz[i];
    if (i < NITSM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }
}

#else

void GPUTPCGMPolynomialField::Print() const
{
  // do nothing
}

#endif
