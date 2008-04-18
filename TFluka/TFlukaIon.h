#ifndef TFLUKAION_H
#define TFLUKAION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// Class that gives access to properties of ions used as primary particles   //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>

class TFlukaIon : public TNamed
{

public:
    TFlukaIon();
    TFlukaIon(const char* name, Int_t z, Int_t a, Int_t q, Double_t exE, Double_t mass = 0.);
    Int_t    GetZ()                const  {return fZ;}
    Int_t    GetA()                const  {return fA;}

    Int_t    GetQ()                const  {return fQ;}
    Double_t GetExcitationEnergy() const  {return fExEnergy;}
    Double_t GetMass()             const  {return fMass;}
    Int_t    GetPdgCode()          const  {return GetIonPdg(fZ, fA);}
    //
    void     WriteUserInputCard(FILE* file) const;
    //
    static void  AddIon(Int_t a, Int_t z);
    static void  AddIon(const char* name, Int_t z, Int_t a, Int_t q,
			Double_t exE, Double_t mass);
    static Int_t GetIonPdg(Int_t z, Int_t a, Int_t i = 0);
    static Int_t    GetZ(Int_t pdg);
    static Int_t    GetA(Int_t pdg);

 protected:
    Int_t    fZ;         // Z
    Int_t    fA;         // A
    Int_t    fQ;         // Q
    Double_t fExEnergy;  // Excitation energy
    Double_t fMass;      // Mass
 private:
    // Copy constructor and operator= declared but not implemented (-Weff++ flag)
    TFlukaIon(const TFlukaIon&);
    TFlukaIon& operator=(const TFlukaIon&);
    
    ClassDef(TFlukaIon, 1)          // Ion Properties
};
	
#endif
	
