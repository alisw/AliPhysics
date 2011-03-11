#ifndef ALIITSMISALIGNER_H
#define ALIITSMISALIGNER_H

// Class building the alignment objects for ITS (SPD, SDD, SSD)
// It derives from AliMisAligner, thus providing the methods
// MakeAlObjsArray (builds and returns the array of alignment objects)
// and GetCDBMetaData (returns the metadata for the OCDB entry)
//

/* $Id$ */

#include <TString.h>
#include "AliMisAligner.h"
#include <TRandom3.h>

class TClonesArray;
class AliAlignObjParams;

//-------------------------------------------------------------------------
class AliITSMisAligner : public AliMisAligner{
    public:
	AliITSMisAligner();
	~AliITSMisAligner() {};

	TClonesArray* MakeAlObjsArray();
	AliCDBMetaData* GetCDBMetaData() const;

	void  SetSeed(Int_t seed) {fRnd.SetSeed(seed); return;}

	void SetWholeITSPars(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fWholeITS[0] = dx;   fWholeITS[1] = dy;     fWholeITS[2] = dz;
	    fWholeITS[3] = dpsi; fWholeITS[4] = dtheta; fWholeITS[5] = dphi;
	}

	void SetSPDSectorSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDSector[0] = dx;   fSPDSector[1] = dy;     fSPDSector[2] = dz;
	    fSPDSector[3] = dpsi; fSPDSector[4] = dtheta; fSPDSector[5] = dphi;
	}

	void SetSPDHSSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDHS[0] = dx;   fSPDHS[1] = dy;     fSPDHS[2] = dz;
	    fSPDHS[3] = dpsi; fSPDHS[4] = dtheta; fSPDHS[5] = dphi;
	}

	void SetSPDLadderSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDLadder[0] = dx;   fSPDLadder[1] = dy;     fSPDLadder[2] = dz;
	    fSPDLadder[3] = dpsi; fSPDLadder[4] = dtheta; fSPDLadder[5] = dphi;
	}

	void SetSPDHBSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDHB[0] = dx;   fSPDHB[1] = dy;     fSPDHB[2] = dz;
	    fSPDHB[3] = dpsi; fSPDHB[4] = dtheta; fSPDHB[5] = dphi;
	}

	void SetSPDBarrelSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDBarrel[0] = dx;   fSPDBarrel[1] = dy;     fSPDBarrel[2] = dz;
	    fSPDBarrel[3] = dpsi; fSPDBarrel[4] = dtheta; fSPDBarrel[5] = dphi;
	}

	void SetSDDLayerSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSDDLayer[0] = dx;   fSDDLayer[1] = dy;     fSDDLayer[2] = dz;
	    fSDDLayer[3] = dpsi; fSDDLayer[4] = dtheta; fSDDLayer[5] = dphi;
	}

	void SetSDDBarrelSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSDDBarrel[0] = dx;   fSDDBarrel[1] = dy;     fSDDBarrel[2] = dz;
	    fSDDBarrel[3] = dpsi; fSDDBarrel[4] = dtheta; fSDDBarrel[5] = dphi;
	}

	void SetSDDLadderSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSDDLadder[0] = dx;   fSDDLadder[1] = dy;     fSDDLadder[2] = dz;
	    fSDDLadder[3] = dpsi; fSDDLadder[4] = dtheta; fSDDLadder[5] = dphi;
	}

	void SetSDDModuleSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSDDModule[0] = dx;   fSDDModule[1] = dy;     fSDDModule[2] = dz;
	    fSDDModule[3] = dpsi; fSDDModule[4] = dtheta; fSDDModule[5] = dphi;
	}

	void SetSSDBarrelPars(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSSDBarrel[0] = dx;   fSSDBarrel[1] = dy;     fSSDBarrel[2] = dz;
	    fSSDBarrel[3] = dpsi; fSSDBarrel[4] = dtheta; fSSDBarrel[5] = dphi;
	}

	void SetSSDLadderSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSSDLadder[0] = dx;   fSSDLadder[1] = dy;     fSSDLadder[2] = dz;
	    fSSDLadder[3] = dpsi; fSSDLadder[4] = dtheta; fSSDLadder[5] = dphi;
	}

	void SetSSDModuleSigmas(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSSDModule[0] = dx;   fSSDModule[1] = dy;     fSSDModule[2] = dz;
	    fSSDModule[3] = dpsi; fSSDModule[4] = dtheta; fSSDModule[5] = dphi;
	}

	void SetSPDLadderShiftT(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDLadderShiftT[0] = dx;   fSPDLadderShiftT[1] = dy;     fSPDLadderShiftT[2] = dz;
	    fSPDLadderShiftT[3] = dpsi; fSPDLadderShiftT[4] = dtheta; fSPDLadderShiftT[5] = dphi;
	}

	void SetSPDLadderShiftB(Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
	{
	    fSPDLadderShiftB[0] = dx;   fSPDLadderShiftB[1] = dy;     fSPDLadderShiftB[2] = dz;
	    fSPDLadderShiftB[3] = dpsi; fSPDLadderShiftB[4] = dtheta; fSPDLadderShiftB[5] = dphi;
	}

	void SetSPDLadderShiftT(Double_t pars[6])
	{
	    for(Int_t ii=0; ii<6; ii++)
		fSPDLadderShiftT[ii] = pars[ii];
	}

	void SetSPDLadderShiftB(Double_t pars[6])
	{
	    for(Int_t ii=0; ii<6; ii++)
		fSPDLadderShiftB[ii] = pars[ii];
	}

	void  SetWholeITSMisAlignment();
	void  SetSPDMisAlignment();
	void  SetSDDMisAlignment();
	void  SetSSDMisAlignment();

	Double_t GetUnif(Double_t x1,Double_t x2) {return fRnd.Uniform(x1,x2);}

	Bool_t AddAlignObj(char* name,Double_t dx,Double_t dy,Double_t dz,
		Double_t dpsi,Double_t dtheta,Double_t dphi,
		const char* distrib);

	Bool_t AddAlignObj(Int_t lay,Double_t dx,Double_t dy,Double_t dz,
		Double_t dpsi,Double_t dtheta,Double_t dphi,
		Bool_t unif);

	Bool_t AddAlignObj(Int_t lay,Int_t ladd,Double_t dx,Double_t dy,Double_t dz,
		Double_t dpsi,Double_t dtheta,Double_t dphi,
		Double_t xShift,Double_t yShift,Double_t zShift,
		Double_t psiShift,Double_t thetaShift,Double_t phiShift,
		Bool_t unif);

	Bool_t AddSectorAlignObj(Int_t sectMin,Int_t sectMax,
		Double_t dx,Double_t dy,Double_t dz,
		Double_t dpsi,Double_t dtheta,Double_t dphi,
		Double_t xShift,Double_t yShift,Double_t zShift,
		Double_t psiShift,Double_t thetaShift,Double_t phiShift,
		Bool_t unif);

	void  ShiftAlignObj(AliAlignObjParams &alObj,Double_t dx,Double_t dy,Double_t dz,Double_t dpsi,Double_t dtheta,Double_t dphi);
	void  SmearAlignObj(AliAlignObjParams &alObj,Double_t sx,Double_t sy,Double_t sz,Double_t spsi,Double_t stheta,Double_t sphi);

	const char*  GetSymbName(Int_t layer) const;
	const char*  GetSymbName(Int_t layer,Int_t ladd) const;
	const char*  GetSymbName(Int_t layer,Int_t ladd,Int_t mod) const;
	const char*  GetHalfStaveLadderSymbName(Int_t layer,Int_t ladd,Int_t halfStave) const;
	static const char*  GetParentSymName(const char* symname) ;
	// const char*  GetSistersSymName(const char* symname) const;
	static Bool_t GetLayerAndLevel(const char* symname, Int_t &layer, Int_t &level);

	static Int_t GetNLayers()             {return kNLayers;}
	static Int_t GetNLadders(Int_t lay)   {return fgkNLadders[lay];}
	static Int_t GetNDetectors(Int_t lay) {return fgkNDetectors[lay];}
	static Int_t GetNSisters(const char* symname);
	static Int_t GetNDaughters(const char* symname);

    protected:
	TRandom3     fRnd; // TRandom3 object
	//TRandom     fRnd; // TRandom object
	Int_t        fInd; // index of current AliAlignObjParams in fAlignObjArray
	TClonesArray *fAlignObjArray; // array of AliAlignObjParams
	TString      fStrSPD; // name of SPD
	TString      fStrSDD; // name of SDD
	TString      fStrSSD; // name of SSD
	TString      fStrStave; // name of SPD stave
	TString      fStrHalfStave; // name of SPD half-stave
	TString      fStrLadder; // name of SPD ladder
	TString      fStrSector; // name of SPD sector
	TString      fStrSensor; // name of sensitive volume

    private:
	AliITSMisAligner(const AliITSMisAligner &mAligner);
	AliITSMisAligner &operator= (const AliITSMisAligner &mAligner);
	enum {kNLayers = 6}; // The number of layers.
	static const Int_t  fgkNLadders[kNLayers];  // Array of the number of ladders/layer(layer)
	static const Int_t  fgkNDetectors[kNLayers];// Array of the number of detector/ladder(layer)

	// Parameters setting the misalignment at all SPD/SDD/SSD levels
	Double_t fWholeITS[6];   // parameters for the whole TIS
	Double_t fSPDSector[6];  // sectors
	Double_t fSPDHB[6];      // SPD half barrel
	Double_t fSPDBarrel[6];  // SPD barrel
	Double_t fSPDHS[6];      // SPD half stave
	Double_t fSPDLadder[6];  // SPD ladder
	Double_t fSDDLayer[6];   // SPD layer
	Double_t fSDDBarrel[6];  // SDD barrel
	Double_t fSDDLadder[6];  // SDD ladder
	Double_t fSDDModule[6];  // SDD module
	Double_t fSSDBarrel[6];  // SSD barrel
	Double_t fSSDLayer[6];   // SSD layer
	Double_t fSSDLadder[6];  // SSD ladder
	Double_t fSSDModule[6];  // SSD module

	// Parameters setting common shifts (used for "full" misalignment)
	Double_t fSPDLadderShiftT[6]; // for top half-barrel ladders
	Double_t fSPDLadderShiftB[6]; // for bottom half-barrel ladders
	Double_t fSDDLadderShift1[6]; // for ladder first SDD layer
	Double_t fSDDLadderShift2[6]; // for ladder second SDD layer

	// Choice between uniform (kTRUE) or gaussian (kFALSE) distribution in the smearing
	Bool_t fUnifSPDSector, fUnifSPDHS, fUnifSDDLadder, fUnifSSDLadder, fUnifSPDLadder, fUnifSDDModule, fUnifSSDModule;

	ClassDef(AliITSMisAligner,0)   //ITS MisAligner
};


#endif
