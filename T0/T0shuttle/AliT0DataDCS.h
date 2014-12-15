#ifndef ALIT0DATADCS_H
#define ALIT0DATADCS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TMap.h>
//#include <TClonesArray.h>
//#include <TGraph.h>
#include <TObject.h>
class TString;

// AliT0DataDCS class
// fetching T0 data points from DCS, calculating mean values for the run
// and storing the result to Reference DB

class AliT0DataDCS : public TObject {
public:
	enum {kNAliases=201, kHV=12, kLV=2, kCFD=12, kScalers=32, kTRM=20, kDRM=5, kAtten=1}; // 201 T0 aliases

	AliT0DataDCS();
	AliT0DataDCS(Int_t nRun, UInt_t startTime, UInt_t endTime, UInt_t startTimeDCSQuery, UInt_t endTimeDCSQuery);
	AliT0DataDCS(const AliT0DataDCS & data);
  	AliT0DataDCS& operator=(const AliT0DataDCS & data);
       	~AliT0DataDCS();

	void SetRun(Int_t run) {fRun = run;}
	void SetStartTime(Int_t startTime) {fStartTime = startTime;}
	void SetEndTime(Int_t endTime) {fEndTime = endTime;}
	void SetStartTimeDCSQuery(Int_t startTimeDCSQuery) {fStartTimeDCSQuery = startTimeDCSQuery;}
	void SetEndTimeDCSQuery(Int_t endTimeDCSQuery) {fEndTimeDCSQuery = endTimeDCSQuery;}


        Int_t GetRun() const {return fRun;}
	Int_t GetStartTime() const {return fStartTime;}
	Int_t GetEndTime() const {return fEndTime;}
	Int_t GetStartTimeDCSQuery() const {return fStartTimeDCSQuery;}
  	Int_t GetEndTimeDCSQuery() const {return fEndTimeDCSQuery;}	

        Int_t GetTVDCtop() const {return    fTVDCtop;}     
	Int_t GetTVDCbottom() const {return fTVDCbottom;} 

        Float_t GetMeanPMTCurrentASide(Int_t pmt)   const {return fHViA[pmt];} 
        Float_t GetMeanPMTVoltageASide(Int_t pmt)   const {return fHVvA[pmt];} 
        Float_t GetMeanLVCurrentASide (Int_t index) const {return fLViA[index];} 
        Float_t GetMeanLVVoltageASide (Int_t index) const {return fLVvA[index];}

        Float_t GetMeanPMTCurrentCSide(Int_t pmt)   const {return fHViC[pmt];} 
        Float_t GetMeanPMTVoltageCSide(Int_t pmt)   const {return fHVvC[pmt];} 
        Float_t GetMeanLVCurrentCSide (Int_t index) const {return fLViC[index];} 
        Float_t GetMeanLVVoltageCSide (Int_t index) const {return fLVvC[index];}

        Float_t GetMeanCFDThresholdASide(Int_t pmt) const{return fCFDtA[pmt];}
	Float_t GetMeanCFDWalkASide   (Int_t pmt)   const{return fCFDwA[pmt];}	
        Float_t GetMeanCFDThresholdCSide(Int_t pmt) const{return fCFDtC[pmt];}
	Float_t GetMeanCFDWalkCSide   (Int_t pmt)   const{return fCFDwC[pmt];}	

        UInt_t GetMeanValueOfT0ScalerFromEntireRun(Int_t scaler) const{return fScalerMean[scaler];}         
        UInt_t GetMeanValueOfT0ScalerPerSecond    (Int_t scaler) const{return fScalerSecMean[scaler];} 

        Float_t GetMeanTemperatureTRM (Int_t trm)   const {return fTRM[trm];}			
	Float_t GetMeanTemperatureDRM (Int_t drm)   const {return fDRM[drm];}	
	Float_t GetLaserAmplitude() const { return  fAtten;}	
	Int_t  GetMultDiscriminatorCentralASide() const     {return fMPDcentA;}
	Int_t  GetMultDiscriminatorCentralCSide() const     {return fMPDcentC;}     
        Int_t  GetMultDiscriminatorSemiCentralASide() const {return fMPDsemiCentA;}
        Int_t  GetMultDiscriminatorSemiCentralCSide() const {return fMPDsemiCentC;}          
	Int_t  GetMultDiscriminatorMode() const {return fMPDmode;}


	Bool_t ProcessData(TMap& aliasMap);

	const char* GetAliasName(Int_t pos) const {return pos<kNAliases ? fAliasNames[pos].Data() : 0;}

        void PrintT0Data() const;

private:
	void Init();
	void Introduce(UInt_t numAlias, const TObjArray* aliasArr)const;
        void PrintfArray(const char *label, const Float_t *array, Int_t numElements) const;

	Int_t fRun;				// Run number
	UInt_t fStartTime;			// Start time
	UInt_t fEndTime;			// End time
	Int_t fStartTimeDCSQuery; 		// Begin time DCSQuery
  	Int_t fEndTimeDCSQuery;   		// End time DCSQuery
	Float_t fHViA[kHV];			// Mean value of HV current in uA on A-side
	Float_t fHVvA[kHV]; 			// Mean value of HV voltage in V on A-side
	Float_t fLViA[kLV]; 			// Mean value of LV current in uA on A-side
	Float_t fLVvA[kLV];			// Mean value of LV voltage in V on A-side
	Float_t fHViC[kHV];			// Mean value of HV current in uA on C-side
        Float_t fHVvC[kHV];			// Mean value of HV voltage in V on C-side
        Float_t fLViC[kLV];			// Mean value of LV current in uA on C-side
        Float_t fLVvC[kLV];			// Mean value of LV voltage in V on C-side
 	Float_t fCFDtA[kCFD];			// Mean threshold on CFD in V on A-side	
	Float_t fCFDwA[kCFD];			// Mean walk on CFD in V on A-side
	Float_t fCFDtC[kCFD];			// Mean threshold on CFD in V on C-side
        Float_t fCFDwC[kCFD];			// Mean walk on CFD in V on C-side
        UInt_t fScalerMean[kScalers];  	// Mean value of T0 scaler counts from the entire run
        UInt_t fScalerSecMean[kScalers]; 	// Mean value of T0 scaler counts per second
	Float_t fTRM[kTRM];			// Mean temperature on TRM in degrees of Celsius
	Float_t fDRM[kDRM];			// Mean temperature on DRM in degrees of Celsius
	Float_t fAtten;				// Laser amplitude in MIPs
	Int_t fMPDcentA;			// Multiplicity Discriminator central on A-side
	Int_t fMPDcentC;                	// Multiplicity Discriminator central on C-side
        Int_t fMPDsemiCentA;                	// Multiplicity Discriminator semi-central on A-side
        Int_t fMPDsemiCentC;                	// Multiplicity Discriminator semi-central on C-side
	Int_t fTVDCtop;				// T0 Vertex Unit top
	Int_t fTVDCbottom;			// T0 Vertex Unit bottom
	Int_t fMPDmode;                       // Multiplicity Discriminator on C-side only, A-side only, or both sides
	TString fAliasNames[kNAliases];		// T0 data points aliases  
	Bool_t fIsProcessed;			// status - was processing data successful
	ClassDef(AliT0DataDCS, 2)

protected:
};

#endif
