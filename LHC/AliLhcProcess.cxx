/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include "AliLhcProcess.h"
#include "AliLHC.h"

ClassImp(AliLhcProcess)

AliLhcProcess::AliLhcProcess(AliLHC *lhc, const char* name, const char* title)
    :TNamed(name,title),
     fAccelerator(lhc)
{
// Constructor
}

AliLhcProcess::AliLhcProcess(const AliLhcProcess& process)
    : TNamed(process), AliLhcMonitor(process), fAccelerator(0)
{
// copy constructor
}

AliLhcProcess::~AliLhcProcess()
{
// Destructor

}

void AliLhcProcess::Evolve(Float_t dt)
{
    printf("\n Here process %s %f:", GetName(), dt);
}

AliLhcProcess& AliLhcProcess::operator=(const  AliLhcProcess & /*rhs*/)
{
// Assignment operator
    return *this;
}


