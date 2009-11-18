#ifndef ALIVZEROTRIGGERDATA_H
#define ALIVZEROTRIGGERDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */

// Class AliVZEROTriggerData
// -------------------------
// Retrieves and hold the FEE parameters
// The parameters are recieved from the shuttle 
// AliVZEROTriggerData is then used in the AliVZEROTriggerSimulator
//

#include <TNamed.h>

#include "AliVZERODataFEE.h"


class AliVZEROTriggerData : public TNamed {
public:
	AliVZEROTriggerData();
	AliVZEROTriggerData(Int_t nRun, UInt_t startTime, UInt_t endTime);
	~AliVZEROTriggerData();
	
	void FillData(AliVZERODataFEE * data);

	// ----- Setters -----

	void SetClk1Win1(UShort_t* clks);
	void SetClk1Win1(UShort_t clk, Int_t board);
	void SetClk2Win1(UShort_t* clks);
	void SetClk2Win1(UShort_t clk, Int_t board);
	
	void SetClk1Win2(UShort_t* clks);
	void SetClk1Win2(UShort_t clk, Int_t board);
	void SetClk2Win2(UShort_t* clks);
	void SetClk2Win2(UShort_t clk, Int_t board);
	
	void SetDelayClk1Win1(UShort_t* delays);
	void SetDelayClk1Win1(UShort_t delay, Int_t board);
	void SetDelayClk2Win1(UShort_t* delays);
	void SetDelayClk2Win1(UShort_t delay, Int_t board);
	
	void SetDelayClk1Win2(UShort_t* delays);
	void SetDelayClk1Win2(UShort_t delay, Int_t board);
	void SetDelayClk2Win2(UShort_t* delays);
	void SetDelayClk2Win2(UShort_t delay, Int_t board);
	
	void SetLatchWin1(UShort_t *latchs);
	void SetLatchWin1(UShort_t latch, Int_t board);
	void SetLatchWin2(UShort_t *latchs);
	void SetLatchWin2(UShort_t latch, Int_t board);
	
	void SetResetWin1(UShort_t *resets);
	void SetResetWin1(UShort_t reset, Int_t board);
	void SetResetWin2(UShort_t *resets);
	void SetResetWin2(UShort_t reset, Int_t board);
	
	void SetPedestalSubtraction(Bool_t *peds);
	void SetPedestalSubtraction(Bool_t ped, Int_t board);
	
	void SetBBAThreshold(UShort_t th) {fBBAThreshold = th;};
	void SetBBCThreshold(UShort_t th) {fBBCThreshold = th;};

	void SetBGAThreshold(UShort_t th) {fBGAThreshold = th;};
	void SetBGCThreshold(UShort_t th) {fBGCThreshold = th;};

	void SetBBAForBGThreshold(UShort_t th) {fBBAForBGThreshold = th;};
	void SetBBCForBGThreshold(UShort_t th) {fBBCForBGThreshold = th;};
	
	void SetCentralityV0AThrLow(UShort_t th) {fCentralityVOAThrLow = th;};
	void SetCentralityV0AThrHigh(UShort_t th) {fCentralityVOAThrHigh = th;};
	
	void SetCentralityV0CThrLow(UShort_t th) {fCentralityVOCThrLow = th;};
	void SetCentralityV0CThrHigh(UShort_t th) {fCentralityVOCThrHigh = th;};
	
	void SetMultV0AThrLow(UShort_t th) {fMultV0AThrLow = th;};
	void SetMultV0AThrHigh(UShort_t th) {fMultV0AThrHigh = th;};
	
	void SetMultV0CThrLow(UShort_t th) {fMultV0CThrLow = th;};
	void SetMultV0CThrHigh(UShort_t th) {fMultV0CThrHigh = th;};
	
	void SetTriggerSelected(UShort_t trigger, Int_t output);
	
	void SetEnableCharge(Bool_t val, Int_t board, Int_t channel);
	void SetEnableTiming(Bool_t val, Int_t board, Int_t channel);
	void SetDiscriThr(UShort_t val, Int_t board, Int_t channel);
	void SetDelayHit(UShort_t val, Int_t board, Int_t channel);
	void SetPedestal(UShort_t val, Int_t integrator, Int_t board, Int_t channel);
	void SetPedestalCut(UShort_t val, Int_t integrator, Int_t board, Int_t channel);

	
	// ----- Getters -----
	
	UShort_t * GetClk1Win1() const {return (UShort_t*)fClk1Win1;};
	UShort_t GetClk1Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk1Win1[board]:0);};
	UShort_t * GetClk2Win1() const {return (UShort_t*)fClk2Win1;};
	UShort_t GetClk2Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk2Win1[board]:0);};

	UShort_t * GetClk1Win2() const {return (UShort_t*)fClk1Win2;};
	UShort_t GetClk1Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk1Win2[board]:0);};
	UShort_t * GetClk2Win2() const {return (UShort_t*)fClk2Win2;};
	UShort_t GetClk2Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fClk2Win2[board]:0);};

	UShort_t * GetDelayClk1Win1() const {return (UShort_t*)fDelayClk1Win1;};
	UShort_t GetDelayClk1Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk1Win1[board]:0);};
	UShort_t * GetDelayClk2Win1() const {return (UShort_t*)fDelayClk2Win1;};
	UShort_t GetDelayClk2Win1(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk2Win1[board]:0);};
	
	UShort_t * GetDelayClk1Win2() const {return (UShort_t*)fDelayClk1Win2;};
	UShort_t GetDelayClk1Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk1Win2[board]:0);};
	UShort_t * GetDelayClk2Win2() const {return (UShort_t*)fDelayClk2Win2;};
	UShort_t GetDelayClk2Win2(Int_t board ) const {return ((board>=0 && board<kNCIUBoards)?fDelayClk2Win2[board]:0);};
	
	UShort_t * GetLatchWin1() const {return (UShort_t*)fLatchWin1;};
	UShort_t GetLatchWin1(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fLatchWin1[board]:0);};
	UShort_t * GetLatchWin2() const {return (UShort_t*)fLatchWin2;};
	UShort_t GetLatchWin2(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fLatchWin2[board]:0);};
	
	UShort_t * GetResetWin1() const {return (UShort_t*)fResetWin1;};
	UShort_t GetResetWin1(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fResetWin1[board]:0);};
	UShort_t * GetResetWin2() const {return (UShort_t*)fResetWin2;};
	UShort_t GetResetWin2(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fResetWin2[board]:0);};
	
	Bool_t * GetPedestalSubtraction() const {return (Bool_t*) fPedestalSubtraction;};
	Bool_t GetPedestalSubtraction(Int_t board ) const  {return ((board>=0 && board<kNCIUBoards)?fPedestalSubtraction[board]:0);};

	UShort_t GetBBAThreshold() const {return fBBAThreshold;};
	UShort_t GetBBCThreshold() const {return fBBCThreshold;};

	UShort_t GetBGAThreshold() const {return fBGAThreshold;};
	UShort_t GetBGCThreshold() const {return fBGCThreshold;};

	UShort_t GetBBAForBGThreshold() const {return fBBAForBGThreshold;};
	UShort_t GetBBCForBGThreshold() const {return fBBCForBGThreshold;};
	
	UShort_t GetCentralityV0AThrLow() const {return fCentralityVOAThrLow;};
	UShort_t GetCentralityV0AThrHigh() const {return fCentralityVOAThrHigh;};
	
	UShort_t GetCentralityV0CThrLow() const {return fCentralityVOCThrLow;};
	UShort_t GetCentralityV0CThrHigh() const {return fCentralityVOCThrHigh;};

	UShort_t GetMultV0AThrLow() const {return fMultV0AThrLow;};
	UShort_t GetMultV0AThrHigh() const {return fMultV0AThrHigh;};

	UShort_t GetMultV0CThrLow() const {return fMultV0CThrLow;};
	UShort_t GetMultV0CThrHigh() const {return fMultV0CThrHigh;};

	UShort_t GetTriggerSelected(Int_t output) const {return ((output>=0 && output<kNTriggerOutputs)?fTriggerSelected[output]:0);};
	
	Bool_t GetEnableCharge(Int_t board, Int_t channel);
	Bool_t GetEnableTiming(Int_t board, Int_t channel);
	UShort_t GetDiscriThr(Int_t board, Int_t channel);
	UShort_t GetDelayHit(Int_t board, Int_t channel);
	UShort_t GetPedestal(Int_t integrator, Int_t board, Int_t channel);
	UShort_t GetPedestalCut(Int_t integrator, Int_t board, Int_t channel);

	enum {
		kNCIUBoards = AliVZERODataFEE::kNCIUBoards,
		kNAliases = AliVZERODataFEE::kNAliases,
		kNTriggerOutputs = 5,
		kNChannels = 8
	};
	
private:
	AliVZEROTriggerData(const AliVZEROTriggerData &/*triggerData*/);
	AliVZEROTriggerData& operator= (const AliVZEROTriggerData &/*triggerData*/);
	
	UShort_t fClk1Win1[kNCIUBoards]; //Profil of the Clock 1  of the Window 1 (BB window)
	UShort_t fClk2Win1[kNCIUBoards]; //Profil of the Clock 2  of the Window 1 (BB window)
	UShort_t fClk1Win2[kNCIUBoards]; //Profil of the Clock 1  of the Window 2 (BG window)
	UShort_t fClk2Win2[kNCIUBoards]; //Profil of the Clock 2  of the Window 2 (BG window)
	UShort_t fDelayClk1Win1[kNCIUBoards]; // Delays of the Clock 1  of the Window 1 (BB window)
	UShort_t fDelayClk2Win1[kNCIUBoards]; // Delays of the Clock 2 of the Window 1 (BB window)
	UShort_t fDelayClk1Win2[kNCIUBoards]; // Delays of the Clock 1  of the Window 2 (BG window)
	UShort_t fDelayClk2Win2[kNCIUBoards]; // Delays of the Clock 2 of the Window 2 (BG window)
	UShort_t fLatchWin1[kNCIUBoards]; //Profil of the Clock of the Latch signal of Window 1 (BB window)
	UShort_t fLatchWin2[kNCIUBoards]; //Profil of the Clock of the Latch signal of Window 2 (BG window)
	UShort_t fResetWin1[kNCIUBoards]; //Profil of the Clock of the Reset signal of Window 1 (BB window)
	UShort_t fResetWin2[kNCIUBoards]; //Profil of the Clock of the Reset signal of Window 2 (BG window)
	Bool_t	 fPedestalSubtraction[kNCIUBoards]; // Flag to en(dis)able pedestal subtraction before centrality trigger calculation
	UShort_t fBBAThreshold;  // Minimum bias Threshold in number of channel hit for V0A
	UShort_t fBBCThreshold;  // Minimum bias Threshold in number of channel hit for V0C
	UShort_t fBGAThreshold;  // Beam Gas Threshold in number of channel hit for V0A
	UShort_t fBGCThreshold;  // Beam Gas Threshold in number of channel hit for V0C
	UShort_t fBBAForBGThreshold;  // BBA threshold for Beam Gas triggers (i.e. BBA and BGC)
	UShort_t fBBCForBGThreshold;  // BBC threshold for Beam Gas triggers (i.e. BBC and BGA)
	UShort_t fCentralityVOAThrLow;  // Threshold used for centrality triggers (i.e. CTA1 and CTC1)
	UShort_t fCentralityVOAThrHigh; // Threshold used for centrality triggers (i.e. CTA2 and CTC2)
	UShort_t fCentralityVOCThrLow;  // Threshold used for centrality triggers (i.e. CTA1 and CTC1)
	UShort_t fCentralityVOCThrHigh; // Threshold used for centrality triggers (i.e. CTA2 and CTC2)
	UShort_t fMultV0AThrLow;  // Threshold used for multiplicity triggers (i.e. MTA and MTC)
	UShort_t fMultV0AThrHigh; // Threshold used for multiplicity triggers (i.e. MTA and MTC)
	UShort_t fMultV0CThrLow;  // Threshold used for multiplicity triggers (i.e. MTA and MTC)
	UShort_t fMultV0CThrHigh; // Threshold used for multiplicity triggers (i.e. MTA and MTC)
	UShort_t fTriggerSelected[kNTriggerOutputs]; // Triggers selected on the 5 outputs to CTP
	Bool_t	 fEnableCharge[kNCIUBoards][kNChannels]; // Flag to know is a channel is participating to the Charge triggers
	Bool_t	 fEnableTiming[kNCIUBoards][kNChannels]; // Flag to know is a channel is participating to the Timing triggers
	UShort_t fDiscriThr[kNCIUBoards][kNChannels]; // Threshold of each discriminator
	UShort_t fDelayHit[kNCIUBoards][kNChannels]; // Individual delays of each channel 
	UShort_t fPedestalOdd[kNCIUBoards][kNChannels]; // Pedestals for the Odd integrators
	UShort_t fPedestalEven[kNCIUBoards][kNChannels]; // Pedestals for the Even integrators
	UShort_t fPedestalCutOdd[kNCIUBoards][kNChannels]; // Pedestals Cut for the Odd integrators
	UShort_t fPedestalCutEven[kNCIUBoards][kNChannels]; // Pedestals Cut for the Even integrators

	Int_t fRun;       // Run number
	Int_t fStartTime; // Start time
	Int_t fEndTime;   // End time
	TString fAliasNames[kNAliases];	// aliases for DCS data
	Bool_t fIsProcessed; // bool to know processing status
	
	Bool_t	IsClkValid(UShort_t clock) const;
	void SetParameter(TString name, Int_t val);
	
	
	ClassDef( AliVZEROTriggerData, 2 )  

};

#endif // ALIVZEROTRIGGERDATA_H


