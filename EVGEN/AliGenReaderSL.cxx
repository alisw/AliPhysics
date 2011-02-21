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
//
// Realisations of the AliGenReader interface to be used with AliGenExFile.
// NextEvent() loops over events 
// and NextParticle() loops over particles. 
// This implementation reads various StarLight output formats 
// Author: andreas.morsch@cern.ch

#include <TVirtualMC.h> 
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliGenReaderSL.h"
#include "AliRun.h"
#include "AliStack.h"


ClassImp(AliGenReaderSL)


AliGenReaderSL& AliGenReaderSL::operator=(const  AliGenReaderSL& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

void AliGenReaderSL::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}

void AliGenReaderSL::Init()
{
    // Initialisation
    if( !(fFile = fopen(fFileName,"r")) ) {
	printf("Couldn't open input file: %s \n", fFileName);
    } else {
	printf("File %s opened \n", fFileName);
    }
    
    
}

Int_t AliGenReaderSL::NextEvent()
{
// Read the next event

    if (fFormat == 0) {
	fNParticles = 4;
	return (fNParticles);
    }
    
// Example 
// EVENT: 7 23 1
// GAMMAENERGIES: 27.6431
// VERTEX: 0 0 0 0 1 0 0 23
    char linelabel[20];
    int i1 = 0;
    int i2 = 0;
    int i3 = 0;


    double x1 = 0.0;
    double x2 = 0.0;
    double x3 = 0.0;
    double x4 = 0.0;

    int ntrk = 0;
    int nvtx = 0;
    //
    Float_t eGamma1, eGamma2; 
    eGamma2 = 0.;
    Int_t nb;
    // Event line 
    nb = fscanf(fFile,"%6s %d %d %d ",linelabel, &i1, &ntrk, &i2);
    if (nb == 0) return (0);
    fNParticles = ntrk;
    //printf("Event line: %s i1 = %5d ntrk = %5d i2 = %5d \n", linelabel, i1, ntrk, i2);

    // Gamma line
    if (fFormat == 1) {
	nb = fscanf(fFile, "%14s %f  ",linelabel, &eGamma1); 
    } else if (fFormat == 2) {
	nb = fscanf(fFile, "%14s %f %f ",linelabel, &eGamma1, &eGamma2); 
    }
    if (nb == 0) return (0);
//    printf("Gamma line: %s Egamma1 = %13.3f Egamma2 = %13.3f \n", linelabel, eGamma1, eGamma2); 

    // Vertex line 
    nb = fscanf(fFile, "%7s %lf %lf %lf %lf %d %d %d %d",
	   linelabel, &x1, &x2, &x3, &x4, &i1, &i2, &i3, &nvtx);
//    printf("Vertex line: %s (x = %13.3f, y =  %13.3f, z =  %13.3f, t =  %13.3f) i1 = %5d i2 = %5d i3 = %5d nvtx = %5d \n", 
//	   linelabel, x1, x2, x3, x4, i1, i2, i3, nvtx);
    if (nb == 0) return (0);
    if(ntrk != nvtx) printf("ERROR: ntrk = %5d  nvtx = %5d \n", ntrk, nvtx);
    
    return (fNParticles);
    
}

TParticle* AliGenReaderSL::NextParticle()
{
    // Read next particle

    Float_t px, py, pz;
    Int_t pdg;
    static TParticle particle;
    Int_t nb;
    if (fFormat == 0) {
	Int_t ievent;
	Int_t ipart;
	nb = fscanf(fFile, "%d %d %d %f %f %f ", &ievent, &ipart, &pdg, &px, &py, &pz); 
//	printf("%5d %5d %5d %13.3f %13.3f %13.3f \n", ievent, ipart, pdg, px, py, pz);
	
	if (pdg == 8) pdg =  211;
	if (pdg == 9) pdg = -211;

    } else {
	char tracklabel[20];
	int i1 = 0;
	int i2 = 0;
	int i3 = 0;
	int i4 = 0;
	nb = fscanf(fFile,"%6s %d %f %f %f %d %d %d %d",
	       tracklabel, &i1, &px, &py, &pz, &i2, &i3, &i4, &pdg);
//	printf("Particle %5d %13.3f %13.3f %13.3f \n",  pdg, px, py, pz);
    }
    if (nb == 0) return (0x0);

    Double_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
    Float_t e = TMath::Sqrt(px * px + py * py + pz * pz + mass * mass);
    particle.SetMomentum(px, py, pz, e);
    particle.SetPdgCode(pdg);
    particle.SetFirstMother(-1);
    particle.SetLastMother(-1);
    particle.SetBit(kTransportBit);
    return &particle;
}

void AliGenReaderSL::RewindEvent()
{
//    fFile->Rewind();
}
