#ifndef AliGRPOBJECT_H
#define AliGRPOBJECT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TMap.h>

#include "AliDCSSensor.h"
#include "AliSplineFit.h"
#include "AliLog.h"

class AliGRPObject : public TObject {
 public:

	enum Stats {kMean = 0, kTruncMedian = 1, kMedian = 2, kSDMean = 3, kSDMedian = 4};

	enum DP_HallProbes { 
		 k_HP_L3_BSF17_H1= 0 , k_HP_L3_BSF17_H2, k_HP_L3_BSF17_H3, k_HP_L3_BSF17_Temperature, 
		 k_HP_L3_BSF4_H1, k_HP_L3_BSF4_H2, k_HP_L3_BSF4_H3, k_HP_L3_BSF4_Temperature, 
		 k_HP_L3_BKF17_H1, k_HP_L3_BKF17_H2, k_HP_L3_BKF17_H3, k_HP_L3_BKF17_Temperature, 
		 k_HP_L3_BKF4_H1, k_HP_L3_BKF4_H2, k_HP_L3_BKF4_H3, k_HP_L3_BKF4_Temperature, 
		 k_HP_L3_BSF13_H1, k_HP_L3_BSF13_H2, k_HP_L3_BSF13_H3, k_HP_L3_BSF13_Temperature,
		 k_HP_L3_BSF8_H1, k_HP_L3_BSF8_H2, k_HP_L3_BSF8_H3, k_HP_L3_BSF8_Temperature,
		 k_HP_L3_BKF13_H1, k_HP_L3_BKF13_H2, k_HP_L3_BKF13_H3, k_HP_L3_BKF13_Temperature,
		 k_HP_L3_BKF8_H1, k_HP_L3_BKF8_H2, k_HP_L3_BKF8_H3, k_HP_L3_BKF8_Temperature,
		 k_HP_Dipole_Inside_H1, k_HP_Dipole_Inside_H2, k_HP_Dipole_Inside_H3, k_HP_Dipole_Inside_Temperature,
		 k_HP_Dipole_Outside_H1, k_HP_Dipole_Outside_H2, k_HP_Dipole_Outside_H3, k_HP_Dipole_Outside_Temperature};


	AliGRPObject();
	AliGRPObject(const AliGRPObject & obj);
	AliGRPObject& operator=(const AliGRPObject & obj);
	~AliGRPObject();

	// getters

	time_t    GetTimeStart() const {return fTimeStart;}
	time_t    GetTimeEnd() const {return fTimeEnd;}
	Float_t   GetBeamEnergy() const {return fBeamEnergy;}
	TString   GetBeamType() const {return fBeamType;}
	Char_t    GetNumberOfDetectors() const {return fNumberOfDetectors;}
	UInt_t     GetDetectorMask() const {return fDetectorMask;}
	TString   GetLHCPeriod() const {return fLHCPeriod;}
	TString   GetRunType() const {return fRunType;}
	TString   GetLHCState() const {return fLHCState;}
	Float_t*  GetLHCLuminosity() {return fLHCLuminosity;}
	Float_t   GetLHCLuminosity(Stats stat) const {return fLHCLuminosity[stat];}
	AliSplineFit*  GetLHCLuminositySplineFit() const {return fLHCLuminositySplineFit;}
	Float_t*  GetBeamIntensity() {return fBeamIntensity;}
	Float_t   GetBeamIntensity(Stats stat) const {return fBeamIntensity[stat];}
	AliSplineFit*  GetBeamIntensitySplineFit() const {return fBeamIntensitySplineFit;}
	Char_t    GetL3Polarity() const {return fL3Polarity;}
	Char_t    GetDipolePolarity() const {return fDipolePolarity;}
	Float_t*  GetL3Current() {return fL3Current;}
	Float_t   GetL3Current(Stats stat) const {return fL3Current[stat];}
	Float_t*  GetDipoleCurrent() {return fDipoleCurrent;}
	Float_t   GetDipoleCurrent(Stats stat) const {return fDipoleCurrent[stat];}
	Float_t*  GetCavernTemperature() {return fCavernTemperature;}
	Float_t   GetCavernTemperature(Stats stat) const {return fCavernTemperature[stat];}
	//	Float_t*  GetCavernAtmosPressure() {return fCavernAtmosPressure;}
	//Float_t   GetCavernAtmosPressure(Stats stat) const {return fCavernAtmosPressure[stat];}
	AliDCSSensor*   GetCavernAtmosPressure() const {return fCavernAtmosPressure;}
	AliDCSSensor*   GetSurfaceAtmosPressure() const {return fSurfaceAtmosPressure;}

	Float_t*  GetHallProbes(DP_HallProbes hp);
	Float_t   GetHallProbes(Int_t hp) const {return fHallProbes[hp];}
	Float_t   GetHallProbes(DP_HallProbes hp, Stats stat) const {return fHallProbes[hp*fPoints+stat];}

	Int_t    GetPoints() const {return fPoints;}
	Int_t    GetDimension() const {return fDimension;}

	// setters

	void SetTimeStart(time_t timeStart)  {fTimeStart = timeStart;}
	void SetTimeEnd(time_t timeEnd)  {fTimeEnd = timeEnd;}
	void SetBeamEnergy(Float_t beamEnergy)  {fBeamEnergy = beamEnergy;}
	void SetBeamType(TString beamType)  {fBeamType = beamType;}
	void SetNumberOfDetectors(Char_t numberOfDetectors)  {fNumberOfDetectors = numberOfDetectors;}
	void SetDetectorMask(UInt_t detectorMask)  {fDetectorMask = detectorMask;}
	void SetLHCPeriod(TString lhcPeriod)  {fLHCPeriod = lhcPeriod;}
	void SetRunType(TString runType)  {fRunType = runType;}
	void SetLHCState(TString lhcState)  {fLHCState = lhcState;}
	void SetLHCLuminosity(Float_t* lhcLuminosity)  {
  		for (Int_t i = 0;i<fPoints;i++) fLHCLuminosity[i] = lhcLuminosity[i];
	}
	void SetLHCLuminosity(Float_t lhcLuminosity, Stats stat)  {fLHCLuminosity[stat] = lhcLuminosity;}
	void SetLHCLuminositySplineFit(AliSplineFit* lhcLuminositySplineFit)  {fLHCLuminositySplineFit = lhcLuminositySplineFit;}
	void SetBeamIntensity(Float_t* beamIntensity)  {
  		for (Int_t i = 0;i<fPoints;i++) fBeamIntensity[i] = beamIntensity[i];
	}
	void SetBeamIntensity(Float_t beamIntensity, Stats stat)  {fBeamIntensity[stat] = beamIntensity;}
	void SetBeamIntensitySplineFit(AliSplineFit* beamIntensitySplineFit)  {fBeamIntensitySplineFit = beamIntensitySplineFit;}
	void SetL3Polarity(Char_t l3Polarity)  {fL3Polarity = l3Polarity;}
	void SetDipolePolarity(Char_t dipolePolarity)  {fDipolePolarity = dipolePolarity;}
	void SetL3Current(Float_t* l3Current)  {
		for (Int_t i = 0;i<fPoints;i++) fL3Current[i] = l3Current[i];
	}
	void SetL3Current(Float_t l3Current, Stats stat)  {fL3Current[stat] = l3Current;}
	void SetDipoleCurrent(Float_t* dipoleCurrent) {
  		for (Int_t i = 0;i<fPoints;i++) fDipoleCurrent[i] = dipoleCurrent[i];
	}
	void SetDipoleCurrent(Float_t dipoleCurrent, Stats stat)  {fDipoleCurrent[stat] = dipoleCurrent;}
	void SetCavernTemperature(Float_t* cavernTemperature)  {
  		for (Int_t i = 0;i<fPoints;i++) fCavernTemperature[i] = cavernTemperature[i];
	}
	void SetCavernTemperature(Float_t cavernTemperature, Stats stat)  {fCavernTemperature[stat] = cavernTemperature;}
	//	void SetCavernAtmosPressure(Float_t* cavernAtmosPressure)  {
	//		for (Int_t i = 0;i<fPoints;i++) fCavernAtmosPressure[i] = cavernAtmosPressure[i];
	//}
//	void SetCavernAtmosPressure(Float_t cavernAtmosPressure, Stats stat)  {fCavernAtmosPressure[stat] = cavernAtmosPressure;}
	void SetCavernAtmosPressure(AliDCSSensor* cavernAtmosPressure)  {fCavernAtmosPressure = cavernAtmosPressure;}
	void SetSurfaceAtmosPressure(AliDCSSensor* surfacePressure)  {fSurfaceAtmosPressure = surfacePressure;}

	void SetHallProbes(DP_HallProbes hp, Float_t hall_probe, Stats stat)  {fHallProbes[hp*fPoints+stat] = hall_probe;}
	void SetHallProbes(Float_t *hall_probe){
		for (Int_t i = 0; i< fDimension; i++) fHallProbes[i] =  hall_probe[i];}

	void SetHallProbes(DP_HallProbes hp, Float_t* hall_probe);  
	void SetPoints(Int_t points) {fPoints = points;}
	void SetDimension(Int_t dimension) {fDimension = dimension;}

	// getters for "invalid" flags

	static Float_t GetInvalidFloat() {return fgkInvalidFloat;}
	static TString GetInvalidString() {return fgkInvalidString;}
	static Int_t GetInvalidInt() {return fgkInvalidInt;}
	static Int_t GetInvalidUInt() {return fgkInvalidUInt;}
	static Char_t GetInvalidChar() {return fgkInvalidChar;}
	static Int_t GetNumberOfHP() {return fgknDCSDP_HallProbes;}
	static const char* GetHPDP(Int_t indexHP) {return fgkDCSDataPoints_HallProbes[indexHP];}

	// to read old GRP object in TMap format

	void ReadValuesFromMap(TMap* map);	

 private:

	static const Float_t fgkInvalidFloat;   // value to identify invalid data - float
	static const TString fgkInvalidString;  // value to identify invalid data - string
	static const Char_t fgkInvalidChar;     // value to identify invalid data - char
	static const Int_t fgkInvalidInt;       // value to identify invalid data - int
	static const Int_t fgkInvalidUInt;       // value to identify invalid data - uint
	static const Int_t   fgknDCSDP_HallProbes;               //! number of dcs dps
	static const char*   fgkDCSDataPoints_HallProbes[];      //! names of dcs dps

	Int_t fPoints;                    // number of statistical quantities to be stored
	Int_t fDimension;                 // dimension of Hall Probes array

	time_t   fTimeStart;              // DAQ_time_start entry from DAQ logbook
 	time_t   fTimeEnd;                // DAQ_time_end entry from DAQ logbook
	Float_t  fBeamEnergy;             // beamEnergy entry from DAQ logbook
	TString  fBeamType;               // beamType entry from DAQ logbook
	Char_t   fNumberOfDetectors;      // numberOfDetectors entry from DAQ logbook
	UInt_t   fDetectorMask;           // detectorMask entry from DAQ logbook
	TString  fLHCPeriod;              // LHCperiod entry from DAQ logbook 
	TString  fRunType;                // RunType entry from DAQ logbook 
	TString  fLHCState;               // LHCState entry from DCS DB
	Float_t*  fLHCLuminosity;         // [fPoints]
	                                  // LHCLuminosity entry from DCS DB
	AliSplineFit*  fLHCLuminositySplineFit;       // LHCLuminosity SplineFit from DCS DB
	Float_t*  fBeamIntensity   ;      // [fPoints]
	                                  // BeamIntensity entry from DCS DB
	AliSplineFit*  fBeamIntensitySplineFit;       // BeamIntensity SplineFit from DCS DB
	Char_t    fL3Polarity;             // L3Polarity entry from DCS DB
	Char_t    fDipolePolarity;         // [fPoints]
	                                  // DipolePolarity entry from DCS DB
	Float_t*  fL3Current;             // [fPoints]
	                                  // L3Current entry from DCS DB
	Float_t*  fDipoleCurrent;         // [fPoints]
                                          // DipoleCurrent entry from DCS DB
	Float_t*  fCavernTemperature;     // [fPoints]
                                          // CavernTemperature entry from DCS DB
	//	Float_t*  fCavernAtmosPressure;   // [fPoints]
                                          // CavernAtmosPressure entry from DCS DB
	AliDCSSensor*  fCavernAtmosPressure;    // CavernAtmosPressure entry from DCS DB
	AliDCSSensor*  fSurfaceAtmosPressure;   // SurfaceAtmosPressure entry from DCS DB

	// Hall Probes

	Float_t* fHallProbes;       //[fDimension] 
	                            // array containg the values for the Hall Probes

	ClassDef(AliGRPObject,1)

};

#endif
