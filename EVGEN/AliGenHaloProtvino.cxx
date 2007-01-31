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

// Read background particles from a boundary source
// Very specialized generator to simulate background from beam halo.
// The input file is a text file specially prepared 
// for this purpose.
// Author: andreas.morsch@cern.ch

#include <stdlib.h>

#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPDGCode.h>
#include <TSystem.h>

#include "AliGenHaloProtvino.h"
#include "AliRun.h"

ClassImp(AliGenHaloProtvino)

AliGenHaloProtvino::AliGenHaloProtvino()
    :AliGenerator(-1), 
     fFile(0),
     fFileName(0),
     fSide(1),
     fRunPeriod(kY3D90),
     fTimePerEvent(1.e-4),
     fNskip(0),
     fZ1(0),
     fZ2(0),
     fG1(0),
     fG2(0),
     fGPASize(0)
{
// Constructor
    
    fName  = "HaloProtvino";
    fTitle = "Halo from LHC Tunnel";
//
//  Read all particles
    fNpart = -1;
    SetAnalog(0);
}

AliGenHaloProtvino::AliGenHaloProtvino(Int_t npart)
    :AliGenerator(npart),
     fFile(0),
     fFileName(0),
     fSide(1),
     fRunPeriod(kY3D90),
     fTimePerEvent(1.e-4),
     fNskip(0),
     fZ1(0),
     fZ2(0),
     fG1(0),
     fG2(0),
     fGPASize(0)
{
// Constructor
    fName = "Halo";
    fTitle= "Halo from LHC Tunnel";
//
    fNpart   = npart;
//
    SetAnalog(0);
}


/*
# Title:    README file for the sources of IR8 machine induced background
# Author:   Vadim Talanov <Vadim.Talanov@cern.ch>
# Modified: 12-12-2000 

0. Overview

	There are three files, named ring.one.beta.[01,10,50].m, which
	contain the lists of background particles, induced by proton losses
	upstream of IP8 in the LHC ring one, for the beta* values of 1, 10
	and 50 m, respectively.

1. File contents

	Each line in the files contains the coordinates of particle track
	crossing with the infinite plane, positioned at z=-1m, together with
	the physical properties of corresponding particle, namely:

	S  - S coordinate of the primary interaction vertex, cm;
	N  - type of the gas nuclei at interaction, 1 is H, 2 - C and 3 - O;
	I  - particle ID in PDG particle numbering scheme;
	W  - particle weight;
	E  - particle kinetic energy, GeV;
	X  - x coordinate of the crossing point, cm;
	Y  - y coordinate of the crossing point, cm;
	Dx - x direction cosine;
	Dy - y direction cosine.

2. Normalisation

	Each file is given per unity of linear density of proton inelastic
	interactions with the gas nuclei, [1 inelastic interaction/m].

# ~/vtalanov/public/README.mib: the end.

*/




