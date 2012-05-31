#ifndef ALIGRPOBJECT_H
#define ALIGRPOBJECT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// AliGRPObject
// class to store the information
// coming from the GRP preprocessor
// 

#include <time.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TArrayD.h>

class TMap;

class AliDCSSensor;
class AliSplineFit;
class AliLog;

class AliGRPObject : public TObject {
 public:

	enum Stats {kMean = 0, kTruncMean = 1, kMedian = 2, kSDMean = 3, kSDMedian = 4};
 	
 	enum StatusBits {kPolConvLHC = BIT(14), kUniformBMap = BIT(15),kConvSqrtSHalfGeV = BIT(16)};

	enum DP_HallProbes { 
		 khpL3bsf17H1= 0 , khpL3bsf17H2, khpL3bsf17H3, khpL3bsf17Temperature, 
		 khpL3bsf4H1, khpL3bsf4H2, khpL3bsf4H3, khpL3bsf4Temperature, 
		 khpL3bkf17H1, khpL3bkf17H2, khpL3bkf17H3, khpL3bkf17Temperature, 
		 khpL3bkf4H1, khpL3bkf4H2, khpL3bkf4H3, khpL3bkf4Temperature, 
		 khpL3bsf13H1, khpL3bsf13H2, khpL3bsf13H3, khpL3bsf13Temperature,
		 khpL3bsf8H1, khpL3bsf8H2, khpL3bsf8H3, khpL3bsfy8Temperature,
		 khpL3bkf13H1, khpL3bkf13H2, khpL3bkf13H3, khpL3bkf13Temperature,
		 khpL3bkf8H1, khpL3bkf8H2, khpL3bkf8H3, khpL3bkf8Temperature,
		 khpDipoleInsideH1, khpDipoleInsideH2, khpDipoleInsideH3, khpDipoleInsideTemperature,
		 khpDipoleOutsideH1, khpDipoleOutsideH2, khpDipoleOutsideH3, khpDipoleOutsideTemperature};


	AliGRPObject();
	AliGRPObject(const AliGRPObject & obj);
	AliGRPObject& operator=(const AliGRPObject & obj);
	~AliGRPObject();

	// getters
	Bool_t    IsBeamEnergyIsSqrtSHalfGeV() const {return TestBit(kConvSqrtSHalfGeV);}
 	Bool_t    IsPolarityConventionLHC() const {return TestBit(kPolConvLHC);}
 	Bool_t    IsUniformBMap() const {return TestBit(kUniformBMap);}
	time_t    GetTimeStart() const {return fTimeStart;}
	time_t    GetTimeEnd() const {return fTimeEnd;}
	Float_t   GetBeamEnergy() const;
	TString   GetBeamType() const {return fBeamType;}
	Char_t    GetNumberOfDetectors() const {return fNumberOfDetectors;}
	UInt_t    GetDetectorMask() const {return fDetectorMask;}
	TString   GetLHCPeriod() const {return fLHCPeriod;}
	TString   GetRunType() const {return fRunType;}
	TString   GetLHCState() const {return fLHCState;}
	TString   GetMachineMode() const {return fMachineMode;}
	TObjArray*   GetLHCStateArray() const {return fLHCStateArray;}
	TObjArray*   GetMachineModeArray() const {return fMachineModeArray;}
	Char_t    GetL3Polarity() const {return fL3Polarity;}
	Char_t    GetDipolePolarity() const {return fDipolePolarity;}
	Float_t*  GetL3Current() const {return fL3Current;}
	Float_t   GetL3Current(Stats stat) const {return fL3Current[stat];}
	Float_t*  GetDipoleCurrent() const {return fDipoleCurrent;}
	Float_t   GetDipoleCurrent(Stats stat) const {return fDipoleCurrent[stat];}
	Float_t*  GetCavernTemperature() const {return fCavernTemperature;}
	Float_t   GetCavernTemperature(Stats stat) const {return fCavernTemperature[stat];}
	AliDCSSensor*   GetCavernAtmosPressure() const {return fCavernAtmosPressure;}
	AliDCSSensor*   GetCavernAtmosPressure2() const {return fCavernAtmosPressure2;}
	AliDCSSensor*   GetSurfaceAtmosPressure() const {return fSurfaceAtmosPressure;}
	AliDCSSensor*   GetBestCavernAtmosPressure() const ;
	AliDCSSensor*   GetBestCavernAtmosPressure(const TTimeStamp& time) const;
	static AliDCSSensor* GetBestCavernAtmosPressure(AliDCSSensor* cavern1, 
	                   AliDCSSensor* cavern2, AliDCSSensor* surface, const TTimeStamp& time);


	Float_t*  GetHallProbesArray(DP_HallProbes hp) const;
	Float_t   GetHallProbes(Int_t hp) const {return fHallProbes[hp];}
	Float_t   GetHallProbes(DP_HallProbes hp, Stats stat) const {return fHallProbes[hp*fPoints+stat];}

	Int_t    GetPoints() const {return fPoints;}
	Int_t    GetDimension() const {return fDimension;}

	Double_t GetMaxTimeLHCValidity() const {return fMaxTimeLHCValidity;}
	//
	TObjArray* GetQATrigClasses()    const {return (TObjArray*)fQATrigClasses;}
	TObjArray* GetQACloningRequest() const {return (TObjArray*)fQACloningRequest;}
	//
	// setters
 	void SetBeamEnergyIsSqrtSHalfGeV(Bool_t v=kTRUE) {SetBit(kConvSqrtSHalfGeV,v);}
 	void SetPolarityConventionLHC(Bool_t v=kTRUE) {return SetBit(kPolConvLHC,v);}
 	void SetUniformBMap(Bool_t v=kTRUE) {return SetBit(kUniformBMap,v);}
	void SetTimeStart(time_t timeStart)  {fTimeStart = timeStart;}
	void SetTimeEnd(time_t timeEnd)  {fTimeEnd = timeEnd;}
	void SetBeamEnergy(Float_t beamEnergy)  {fBeamEnergy = beamEnergy;}
	void SetBeamType(TString beamType)  {fBeamType = beamType;}
	void SetNumberOfDetectors(Char_t numberOfDetectors)  {fNumberOfDetectors = numberOfDetectors;}
	void SetDetectorMask(UInt_t detectorMask)  {fDetectorMask = detectorMask;}
	void SetLHCPeriod(TString lhcPeriod)  {fLHCPeriod = lhcPeriod;}
	void SetRunType(TString runType)  {fRunType = runType;}
	void SetLHCState(TString lhcState)  {fLHCState = lhcState;}
	void SetMachineMode(TString machineMode)  {fMachineMode = machineMode;}
	void SetLHCStateArray(TObjArray* lhcStateArray)  {fLHCStateArray = lhcStateArray;}
	void SetMachineModeArray(TObjArray* machineModeArray)  {fMachineModeArray = machineModeArray;}
	void SetL3Polarity(Char_t l3Polarity)  {fL3Polarity = l3Polarity;}
	void SetDipolePolarity(Char_t dipolePolarity)  {fDipolePolarity = dipolePolarity;}
	void SetL3Current(const Float_t* l3Current)  {
		for (Int_t i = 0;i<fPoints;i++) fL3Current[i] = l3Current[i];
	}
	void SetL3Current(Float_t l3Current, Stats stat)  {fL3Current[stat] = l3Current;}
	void SetDipoleCurrent(const Float_t* dipoleCurrent) {
  		for (Int_t i = 0;i<fPoints;i++) fDipoleCurrent[i] = dipoleCurrent[i];
	}
	void SetDipoleCurrent(Float_t dipoleCurrent, Stats stat)  {fDipoleCurrent[stat] = dipoleCurrent;}
	void SetCavernTemperature(const Float_t* cavernTemperature)  {
  		for (Int_t i = 0;i<fPoints;i++) fCavernTemperature[i] = cavernTemperature[i];
	}
	void SetCavernTemperature(Float_t cavernTemperature, Stats stat)  {fCavernTemperature[stat] = cavernTemperature;}
	void SetCavernAtmosPressure(AliDCSSensor* const cavernAtmosPressure)  {fCavernAtmosPressure = cavernAtmosPressure;}
	void SetCavernAtmosPressure2(AliDCSSensor* const cavernAtmosPressure)  {fCavernAtmosPressure2 = cavernAtmosPressure;}
	void SetSurfaceAtmosPressure(AliDCSSensor* const surfacePressure)  {fSurfaceAtmosPressure = surfacePressure;}

	void SetHallProbes(DP_HallProbes hp, Float_t hall_probe, Stats stat)  {fHallProbes[hp*fPoints+stat] = hall_probe;}
	void SetHallProbes(const Float_t* hall_probe){
		for (Int_t i = 0; i< fDimension; i++) fHallProbes[i] =  hall_probe[i];}

	void SetHallProbes(DP_HallProbes hp, const Float_t* hall_probe);  
	void SetPoints(Int_t points) {fPoints = points;}
	void SetDimension(Int_t dimension) {fDimension = dimension;}
	void SetMaxTimeLHCValidity(Double_t maxTimeLHCValidity) {fMaxTimeLHCValidity = maxTimeLHCValidity;}
	//
	void SetQATrigClasses(TObjArray* arr)    {fQATrigClasses = arr;}
	void SetQACloningRequest(TObjArray* arr) {fQACloningRequest = arr;}
	//
	// getters for "invalid" flags

	static Float_t GetInvalidFloat() {return fgkInvalidFloat;}
	static TString GetInvalidString() {return fgkInvalidString;}
	static Int_t GetInvalidInt() {return fgkInvalidInt;}
	static Int_t GetInvalidUInt() {return fgkInvalidUInt;}
	static Char_t GetInvalidChar() {return fgkInvalidChar;}
	static Int_t GetNumberOfHP() {return fgknDCSDPHallProbes;}
	static const char* GetHPDP(Int_t indexHP) {return fgkDCSDataPointsHallProbes[indexHP];}

        Double_t EvalCavernPressure(const TTimeStamp& time, Bool_t& inside) const;
	static Double_t EvalCavernPressure(AliDCSSensor* cavern1, 
	              AliDCSSensor* cavern2, AliDCSSensor* surface, 
		      const TTimeStamp& time, Bool_t& inside);

	// to read old GRP object in TMap format

	void ReadValuesFromMap(const TMap* map);	

	void SetSingleBeamType(Int_t ibeamType, TString beamType)  {fSeparateBeamType[ibeamType] = beamType;}
	TString   GetSingleBeamType(Int_t ibeamType) const {return fSeparateBeamType[ibeamType];}

	void SetNFalseDataQualityFlag(Int_t nFalses) {fNFalseDataQualityFlag = nFalses;}
	Int_t GetNFalseDataQualityFlag() const {return fNFalseDataQualityFlag;}

	void SetFalseDataQualityFlagPeriods(Double_t* falses);
	TArrayD* GetFalseDataQualityFlagPeriods() const {return fFalseDataQualityFlag;}

	Double_t GetStartFalseDataQualityFlag(Int_t iperiod) const;
	Double_t GetEndFalseDataQualityFlag(Int_t iperiod) const;

	void SetBeamTypeFromLHC(TString beamTypeFromLHC)  {fBeamTypeFromLHC = beamTypeFromLHC;}
	TString   GetBeamTypeFromLHC() const {return fBeamTypeFromLHC;}

 private:

	static const Float_t fgkInvalidFloat;   // value to identify invalid data - float
	static const TString fgkInvalidString;  // value to identify invalid data - string
	static const Char_t fgkInvalidChar;     // value to identify invalid data - char
	static const Int_t fgkInvalidInt;       // value to identify invalid data - int
	static const Int_t fgkInvalidUInt;       // value to identify invalid data - uint
	static const Int_t   fgknDCSDPHallProbes;               //! number of dcs dps
	static const char*   fgkDCSDataPointsHallProbes[];      //! names of dcs dps

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
	Char_t    fL3Polarity;            // L3Polarity entry from DCS DB
	Char_t    fDipolePolarity;        // DipolePolarity entry from DCS DB 	                                  
	Float_t*  fL3Current;             // [fPoints]
	                                  // L3Current entry from DCS DB
	Float_t*  fDipoleCurrent;         // [fPoints]
                                          // DipoleCurrent entry from DCS DB
	Float_t*  fCavernTemperature;     // [fPoints]
                                          // CavernTemperature entry from DCS DB
	AliDCSSensor*  fCavernAtmosPressure;    // CavernAtmosPressure entry from DCS DB
	AliDCSSensor*  fCavernAtmosPressure2;    // CavernAtmosPressure2 entry from DCS DB
	AliDCSSensor*  fSurfaceAtmosPressure;   // SurfaceAtmosPressure entry from DCS DB

	// Hall Probes

	Float_t* fHallProbes;       //[fDimension] 
	                            // array containg the values for the Hall Probes

	TString  fMachineMode;      // Machine Mode from LHC
	TObjArray* fLHCStateArray;     // Array of values for the LHC State
	TObjArray* fMachineModeArray;  // Array of values for the LHC State
	TObjArray* fQATrigClasses;     // RS: Array of trigger class to watch in QA
	TObjArray* fQACloningRequest;  // RS: Array of cloning requests for QA histos
	Double_t fMaxTimeLHCValidity;    // time until which the LHC Data Machine Mode and Beam Mode didn't change 
	TString  fSeparateBeamType[2];   // separate beam Types from LHC
	Int_t fNFalseDataQualityFlag;    // number of times the data quality flag turned to FALSE
	TArrayD* fFalseDataQualityFlag;  // array of starts (even positions) and ends (odd poistions) of the periods
	                                 // when the data quality flag was FALSE
	TString fBeamTypeFromLHC;        // string containing the information about the beam types AS SENT BY LHC (in the form "beam1-beam2")
	
	ClassDef(AliGRPObject,10)

};

#endif
