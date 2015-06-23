#ifndef ALIADTRIGGERSIMULATOR_H
#define ALIADTRIGGERSIMULATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// 
// Class AliADTriggerSimulator
// ------------------------------
//  Simulate the AD Trigger response
// Use FEE parameters stored in Database
// Can work on real data or in simulation
//

#include <TObject.h>
#include "AliADConst.h"

class AliADLogicalSignal;
class AliADCalibData;
class TTree;
class TClonesArray;

class AliADTriggerSimulator : public TObject {
public:
	AliADTriggerSimulator();
	AliADTriggerSimulator(TTree * digitsTree, TClonesArray* digits);
	~AliADTriggerSimulator();
	
	AliADCalibData * GetCalibData() const {return fCalibData;};
	
	Bool_t GetBBAandBBC() const		{return (fTriggerWord & 0x1);};
	Bool_t GetBBAorBBC() const		{return (fTriggerWord>>1 & 0x1);};
	Bool_t GetBGAandBBC() const		{return (fTriggerWord>>2 & 0x1);};
	Bool_t GetBGA() const			{return (fTriggerWord>>3 & 0x1);};
	Bool_t GetBGCandBBA() const		{return (fTriggerWord>>4 & 0x1);};
	Bool_t GetBGC() const			{return (fTriggerWord>>5 & 0x1);};
	Bool_t GetCTA1andCTC1() const	{return (fTriggerWord>>6 & 0x1);};
	Bool_t GetCTA1orCTC1() const	{return (fTriggerWord>>7 & 0x1);};
	Bool_t GetCTA2andCTC2() const	{return (fTriggerWord>>8 & 0x1);};
	Bool_t GetCTA2orCTC2() const	{return (fTriggerWord>>9 & 0x1);};
	Bool_t GetMTAandMTC() const		{return (fTriggerWord>>10 & 0x1);};
	Bool_t GetMTAorMTC() const		{return (fTriggerWord>>11 & 0x1);};
	Bool_t GetBBA() const			{return (fTriggerWord>>12 & 0x1);};
	Bool_t GetBBC() const			{return (fTriggerWord>>13 & 0x1);};
	Bool_t GetBGAorBGC() const		{return (fTriggerWord>>14 & 0x1);};
	Bool_t GetBeamGas() const		{return (fTriggerWord>>15 & 0x1);};
	
	void SetBBAandBBC()		{ (fTriggerWord += 0x1);};
	void SetBBAorBBC()		{ (fTriggerWord += 0x1<<1);};
	void SetBGAandBBC()		{ (fTriggerWord += 0x1<<2);};
	void SetBGA()			{ (fTriggerWord += 0x1<<3);};
	void SetBGCandBBA()		{ (fTriggerWord += 0x1<<4);};
	void SetBGC()			{ (fTriggerWord += 0x1<<5);};
	void SetCTA1andCTC1()	{ (fTriggerWord += 0x1<<6);};
	void SetCTA1orCTC1()	{ (fTriggerWord += 0x1<<7);};
	void SetCTA2andCTC2()	{ (fTriggerWord += 0x1<<8);};
	void SetCTA2orCTC2()	{ (fTriggerWord += 0x1<<9);};
	void SetMTAandMTC()		{ (fTriggerWord += 0x1<<10);};
	void SetMTAorMTC()		{ (fTriggerWord += 0x1<<11);};	
	void SetBBA()			{ (fTriggerWord += 0x1<<12);};
	void SetBBC()			{ (fTriggerWord += 0x1<<13);};
	void SetBGAorBGC()		{ (fTriggerWord += 0x1<<14);};
	void SetBeamGas()		{ (fTriggerWord += 0x1<<15);};
	
	void Run();
	void FillFlags(Bool_t *bbFlag, Bool_t *bgFlag, Float_t time[16]);
	virtual void Print(Option_t* /* opt */) const;
	
private:
	// Private methods
	AliADTriggerSimulator(const AliADTriggerSimulator &/*triggerSim*/);
	AliADTriggerSimulator& operator= (const AliADTriggerSimulator & /*triggerSim*/);
	AliADCalibData * LoadCalibData() const ;
	void                  LoadClockOffset();
	void GenerateBBWindows();
	void GenerateBGWindows();
	void GenerateBCMask();
	Bool_t AreGatesOpen() const;
	
	// Members
	AliADLogicalSignal * fBBGate[kNCIUBoards];  // BB Observation window
	AliADLogicalSignal * fBGGate[kNCIUBoards];  // BG Observation window
	AliADLogicalSignal * fBCMask[kNCIUBoards];  // BC Mask

	AliADCalibData *fCalibData; // Object holding the trigger configuration parameters
	Float_t fWindowOffset[kNCIUBoards]; // TDC clock offset including roll-over, trig count and L0->L1 delay
	
	TTree* fDigitsTree; //Pointer to AD digit tree
	TClonesArray* fDigits; //Pointer to AD digit array
	
	Float_t fTime[16];
	Bool_t fBBFlags[16]; // Individual BB Flags
	Bool_t fBGFlags[16]; // Individual BG Flags
	Float_t  fCharges[16]; // Individual Charge
	
	UShort_t fTriggerWord; // Word holding the 16 triggers return by the FEE
		
	ClassDef( AliADTriggerSimulator, 4 )  

};


#endif // ALIADTRIGGERSIMULATOR_H


