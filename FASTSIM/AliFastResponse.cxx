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

/*
$Log$
Revision 1.1  2002/09/20 13:32:51  morsch
Base classes for fast simulation. First commit.

*/


#include "AliFastResponse.h"
#include "AliFastParticle.h"

ClassImp(AliFastResponse)


Float_t AliFastResponse::Evaluate(AliFastParticle* part)
{
//
// Basic implementation of this method 
//
    Float_t theta = part->Theta();
    Float_t phi   = part->Phi();
    Float_t pt    = part->Pt();
    Float_t eff   = Evaluate(pt, theta, phi);
    return eff;
}

void AliFastResponse::Evaluate(Float_t   p,  Float_t  theta , Float_t   phi,
			       Float_t& pS,  Float_t& thetaS, Float_t&  phiS)
{
//    
// Basic implementation of this method 
//
    pS     = p;
    thetaS = theta;
    phiS   = phi;
}

