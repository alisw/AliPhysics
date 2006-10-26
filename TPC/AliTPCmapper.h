#ifndef AliTPCmapper_H
#define AliTPCmapper_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class AliTPCmapper : public TObject{


public:

    AliTPCmapper();
    virtual ~AliTPCmapper();

    void Init();

    Int_t ReadMapping();


    Int_t GetRow(Int_t rcu, Int_t branch, Int_t fec, Int_t altro, Int_t channel) const{
	return fAddressToRow[rcu][branch][fec][altro][channel];}

    Int_t GetPad(Int_t rcu, Int_t branch, Int_t fec, Int_t altro, Int_t channel) const {
	return fAddressToPad[rcu][branch][fec][altro][channel];}

    Int_t GetPadSector(Int_t rcu, Int_t branch, Int_t fec, Int_t altro, Int_t channel) const {
	return fRowPadToPadsec[fAddressToPad[rcu][branch][fec][altro][channel]]
	    [fAddressToPad[rcu][branch][fec][altro][channel]];}


    Int_t GetRow(Int_t altroaddr) const;

    Int_t GetPad(Int_t altroaddr) const;

    Int_t GetPadFromPadSector(Int_t padsector) const{
	return fPadsecToPad[padsector];}

    Int_t GetRowFromPadSector(Int_t padsector) const {
	return fPadsecToRow[padsector];}

    Int_t GetPadSector(Int_t row,Int_t pad) const{
	return fRowPadToPadsec[row][pad];}

    Int_t GetPadsInRowS(Int_t row) const;

    Double_t GetPadXlocalS (Int_t row, Int_t pad) const;
    Double_t GetPadXlocalS (Int_t padsector) const;
    Double_t GetPadYlocalS (Int_t row, Int_t pad) const;
    Double_t GetPadYlocalS (Int_t padsector) const;
    Double_t GetPadXglobalS(Int_t row, Int_t pad, Int_t sector) const;
    Double_t GetPadYglobalS(Int_t row, Int_t pad, Int_t sector) const;
    Double_t GetPadWidthS  (Int_t row) const;
    Double_t GetPadLengthS (Int_t row) const;


    Int_t GetRCUs(Int_t row, Int_t pad) const{
	return fRowPadToRCU[row][pad];}

    Int_t GetBranchS(Int_t row, Int_t pad) const {
	return fRowPadToBranch[row][pad];}

    Int_t GetFECs(Int_t row, Int_t pad) const {
	return fRowPadToFEC[row][pad];}

    Int_t GetAltroS(Int_t row, Int_t pad) const{
	return fRowPadToAltro[row][pad];}

    Int_t GetChannelS(Int_t row, Int_t pad) const {
	return fRowPadToChannel[row][pad];}

    void PrintRBFACinfo(Int_t row, Int_t pad);

    void PrintAddressArray(Int_t row, Int_t pad);


    Int_t GetAltroAddrwPatch(const Int_t row, const Int_t pad) const;

    Int_t GetAltroAddrwPatch(const Int_t padsector) const;

    //for aliroot compatibility (sector == roc)
    Int_t GetPadsInRowS(Int_t row, Int_t sector) const {
	return GetPadsInRowS(row+(sector/36)*63);}
    Double_t GetPadXlocal  (Int_t row, Int_t pad, Int_t sector) const {
	return GetPadXlocalS(row+(sector/36)*63,pad);}
    Double_t GetPadYlocal  (Int_t row, Int_t pad, Int_t sector) const {
	return GetPadYlocalS(row+(sector/36)*63,pad);}
    Double_t GetPadXglobal (Int_t row, Int_t pad, Int_t sector) const {
    	return GetPadXlocalS(row+(sector/36)*63,pad);}
    Double_t GetPadYglobal (Int_t row, Int_t pad, Int_t sector) const{
    	return GetPadYlocalS(row+(sector/36)*63,pad);}
    Int_t GetRCU(Int_t row, Int_t pad, Int_t sector) const {
    	return GetRCUs(row+(sector/36)*63,pad);}
    Int_t GetBranch(Int_t row, Int_t pad, Int_t sector) const {
    	return GetBranchS(row+(sector/36)*63,pad);}
    Int_t GetFEC(Int_t row, Int_t pad, Int_t sector) const {
    	return GetFECs(row+(sector/36)*63,pad);}
    Int_t GetAltro(Int_t row, Int_t pad, Int_t sector) const {
    	return GetAltroS(row+(sector/36)*63,pad);}
    Int_t GetChannel(Int_t row, Int_t pad, Int_t sector) const {
    	return GetChannelS(row+(sector/36)*63,pad);}
    Double_t GetPadWidth   (Int_t row, Int_t sector) const {
    	return GetPadWidthS(row+(sector/36)*63);}
    Double_t GetPadLength  (Int_t row, Int_t sector) const{
    	return GetPadLengthS(row+(sector/36)*63);}

private:
    enum {
	kNrcu       = 6,
	kNbranch    = 2,
	kNfecMax    = 13,
	kNaltro     = 8,
	kNchannel   = 16,
	kNpadrow    = 159,
	kNpadMax    = 140,
	kNaddrSize  = 20,
	kNpadSector = 15488
    };

    Int_t fAddressToRow[kNrcu][kNbranch][kNfecMax][kNaltro][kNchannel]; //fAddressToRow
    Int_t fAddressToPad[kNrcu][kNbranch][kNfecMax][kNaltro][kNchannel]; //fAddressToPad

    Int_t fRowPadToRCU[kNpadrow][kNpadMax]; //fRowPadToRCU
    Int_t fRowPadToBranch[kNpadrow][kNpadMax]; //fRowPadToBranch
    Int_t fRowPadToFEC[kNpadrow][kNpadMax];   //fRowPadToFEC
    Int_t fRowPadToAltro[kNpadrow][kNpadMax]; // fRowPadToAltro
    Int_t fRowPadToChannel[kNpadrow][kNpadMax]; //RowPadToChannel

    Int_t fPadsecToRow[kNpadSector];  //PadsecToRow
    Int_t fPadsecToPad[kNpadSector];  //PadsecToPad

    Int_t fRowPadToPadsec[kNpadrow][kNpadMax]; //RowPadToPadsec

    Char_t fMapfileName[255];                  //MapfileName



	ClassDef(AliTPCmapper,0)
};

#endif
