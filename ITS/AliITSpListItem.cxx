/**************************************************************************
 * Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
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
#include <Riostream.h>
#include <TMath.h>
#include "AliLog.h"
#include "AliITSpListItem.h"
// ************************************************************************
// the data member "fa" of the AliITSpList class
// is a TObjectArray of AliITSpListItem objects
// each AliITSpListItem object contains information related to
// the digitization such signal, noise, module number,...
// plus some info related to the simulation (hits/track)
// in order to allow efficiency studies
//************************************************************************
ClassImp(AliITSpListItem)
//______________________________________________________________________
AliITSpListItem::AliITSpListItem():
fmodule(-1),
findex(-1),
fTsignal(0.0),
fNoise(0.0),
fSignalAfterElect(0.0){
    // Default constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A zeroed/empty AliITSpListItem class.

    for(Int_t i=0;i<this->fgksize;i++){
        this->fTrack[i]  = -2;
        this->fHits[i]   = -1;
        this->fSignal[i] = 0.0;
    } // end if i
}
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(Int_t module,Int_t index,Double_t noise):
fmodule(module),
findex(index),
fTsignal(0.0),
fNoise(noise),
fSignalAfterElect(0.0){
    // Standard noise constructor
    // Inputs:
    //    Int_t module   The module where this noise occurred
    //    Int_t index    The cell index where this noise occurred
    //    Double_t noise The value of the noise.
    // Outputs:
    //    none.
    // Return:
    //    A setup and noise filled AliITSpListItem class.

    for(Int_t i=0;i<this->fgksize;i++){
        this->fTrack[i]  = -2;
        this->fSignal[i] = 0.0;
        this->fHits[i]   = -1;
    } // end if i
}
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(Int_t track,Int_t hit,Int_t module,
                               Int_t index,Double_t signal):
fmodule(module),
findex(index),
fTsignal(signal),
fNoise(0.0),
fSignalAfterElect(0.0){
    // Standard signal constructor
    // Inputs:
    //    Int_t track     The track number which produced this signal
    //    Int_t hit       The hit number which produced this signal
    //    Int_t module    The module where this signal occurred
    //    Int_t index     The cell index where this signal occurred
    //    Double_t signal The value of the signal (ionization)
    // Outputs:
    //    none.
    // Return:
    //    A setup and signal filled  AliITSpListItem class.

    this->fTrack[0]  = track;
    this->fHits[0]   = hit;
    this->fSignal[0] = signal;
    for(Int_t i=1;i<this->fgksize;i++){
        this->fTrack[i]  = -2;
        this->fSignal[i] = 0.0;
        this->fHits[i]   = -1;
    } // end if i
}
//______________________________________________________________________
AliITSpListItem::~AliITSpListItem(){
    // Destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A properly destroyed AliITSpListItem class.

}
//______________________________________________________________________
AliITSpListItem& AliITSpListItem::operator=(const AliITSpListItem &source){
    // = operator
    // Inputs:
    //    AliITSpListItem &source   A AliITSpListItem Object
    // Outputs:
    //    none.
    // Return:
    //    A copied AliITSpListItem object
  this->~AliITSpListItem();
  new(this) AliITSpListItem(source);
  return *this;

}
//______________________________________________________________________
AliITSpListItem::AliITSpListItem(const AliITSpListItem &source) : 
TObject(source),
fmodule(source.fmodule),
findex(source.findex),
fTsignal(source.fTsignal),
fNoise(source.fNoise),
fSignalAfterElect(source.fSignalAfterElect){
    // Copy operator
    // Inputs:
    //    AliITSpListItem &source   A AliITSpListItem Object
    // Outputs:
    //    none.
    // Return:
    //    A copied AliITSpListItem object
  
    for(Int_t i=0;i<this->fgksize;i++){
      this->fTrack[i]  = source.fTrack[i];
      this->fSignal[i] = source.fSignal[i];
      this->fHits[i]   = source.fHits[i];
    } // end if i

}
//______________________________________________________________________
void AliITSpListItem::AddSignal(Int_t track,Int_t hit,Int_t module,
                               Int_t index,Double_t signal){
    // Adds this track number and signal to the pList and orders them
    // Inputs:
    //    Int_t track     The track number which produced this signal
    //    Int_t hit       The hit number which produced this signal
    //    Int_t module    The module where this signal occurred
    //    Int_t index     The cell index where this signal occurred
    //    Double_t signal The value of the signal (ionization)
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t    i,j,trk,hts;
    Double_t sig;
    Bool_t   flg=kFALSE;

    if (TMath::Abs(signal)>2147483647.0) {
      //PH 2147483647 is the max. integer
      //PH This apparently is a problem which needs investigation
      AliWarning(Form("Too big or too small signal value %f",signal));
      signal = TMath::Sign((Double_t)2147483647,signal);
    }
    if(findex!=index || fmodule!=module) 
        Warning("AddSignal","index=%d != findex=%d or module=%d != fmodule=%d",
                 index,findex,module,fmodule);
    fTsignal += signal; // Keep track of sum signal.

    //    for(i=0;i<fgksize;i++) if( track==fTrack[i] && hit==fHits[i] ){
    for(i=0;i<fgksize;i++) if( track==fTrack[i]  ){
        fSignal[i] += signal;
        flg = kTRUE;
    } // end for i & if.
    //cout << "track="<<track<<endl;
    if(flg){ // resort arrays.  
        for(i=1;i<fgksize;i++){
            j = i;
            while(j>0 && fSignal[j]>fSignal[j-1]){
                trk = fTrack[j-1];
                hts = fHits[j-1];
                sig = fSignal[j-1];
                fTrack[j-1]  = fTrack[j];
                fHits[j-1]   = fHits[j];
                fSignal[j-1] = fSignal[j];                
                fTrack[j]  = trk;
                fHits[j]   = hts;
                fSignal[j] = sig;
		//cout << "#fTrack["<<j-1<<"]="<<fTrack[j-1]<< " fTrack["<<
		// j<<"]="<<fTrack[j]<<endl;
                j--;
            } // end while
        } // end if i
        return;
    } // end if added to existing and resorted array

    // new entry add it in order.
    // if this signal is <= smallest then don't add it.
    if(signal <= fSignal[fgksize-1]) return;
    for(i=fgksize-2;i>=0;i--){
        if(signal > fSignal[i]){
            fSignal[i+1] = fSignal[i];
            fTrack[i+1]  = fTrack[i];
            fHits[i+1]   = fHits[i];
        }else{
            fSignal[i+1] = signal;
            fTrack[i+1]  = track;
            fHits[i+1]   = hit;
            return; // put it in the right place, now exit.
        } //  end if
	//cout << "$fTrack["<<i+1<<"]="<<fTrack[i+1]<< " fTrack["<<i<<"]="
	//<<fTrack[i]<< " fHits["<<i+1<<"]="<<fHits[i+1]<< " fHits["<<i<<"]="
	//<<fHits[i]<< " fSignal["<<i+1<<"]="<<fSignal[i+1]<< " fSignal["<<i
	//<<"]="<<fSignal[i]<<endl;
    } // end if; end for i
    // Still haven't found the right place. Must be at top of list.
    fSignal[0] = signal;
    fTrack[0]  = track;
    fHits[0]   = hit;
    //cout << "$fTrack["<<0<<"]="<<fTrack[0]<<" fHits["<<0<<"]="<<fHits[0]
    //<<" fSignal["<<0<<"]="<<fSignal[0]<<endl;
    return;
}
//______________________________________________________________________
void AliITSpListItem::AddNoise(Int_t module,Int_t index,Double_t noise){
    // Adds noise to this existing list.
    // Inputs:
    //    Int_t module   The module where this noise occurred
    //    Int_t index    The cell index where this noise occurred
    //    Double_t noise The value of the noise.
    // Outputs:
    //    none.
    // Return:
    //    none.

    if(findex!=index || fmodule!=module) 
        Warning("AddNoise","index=%d != findex=%d or module=%d != fmodule=%d",
            index,findex,module,fmodule);
    fNoise += noise; // Keep track of sum signal.
}
//______________________________________________________________________
void AliITSpListItem::AddSignalAfterElect(Int_t module,Int_t index,Double_t signal){
    // Adds signal after electronics to this existing list.
    // Inputs:
    //    Int_t module   The module where this noise occurred
    //    Int_t index    The cell index where this noise occurred
    //    Double_t signal The value of the signal.
    // Outputs:
    //    none.
    // Return:
    //    none.

    if(findex!=index || fmodule!=module) 
        Warning("AddSignalAfterElect","index=%d != findex=%d or module=%d "
		"!= fmodule=%d",index,findex,module,fmodule);
    fSignalAfterElect += signal; // Keep track of sum signal.
}
//______________________________________________________________________
void AliITSpListItem::Add(AliITSpListItem *pl){
    // Adds the contents of pl to this
    // pl could come from different module and index 
    // Inputs:
    //    AliITSpListItem *pl  an AliITSpListItem to be added to this class.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;
    Double_t sig  = 0.0;
    Double_t sigT = 0.0;

    for(i=0;i<pl->GetNsignals();i++){
        sig = pl->GetSignal(i); 
        if( sig <= 0.0 ) break; // no more signals
        AddSignal(pl->GetTrack(i),pl->GetHit(i),fmodule,findex,sig);
        sigT += sig;
    } // end for i
    fTsignal += (pl->fTsignal - sigT);
    fNoise   += pl->fNoise;
    return;
}
//______________________________________________________________________
void AliITSpListItem::AddTo(Int_t fileIndex,AliITSpListItem *pl){
    // Adds the contents of pl to this with track number off set given by
    // fileIndex.
    // Inputs:
    //    Int_t fileIndex      track number offset value
    //    AliITSpListItem *pl  an AliITSpListItem to be added to this class.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,trk;
    Double_t sig  = 0.0;

    Int_t module = pl->GetModule();
    Int_t index  = pl->GetIndex();
    for(i=0;i<pl->GetNsignals();i++){
        sig = pl->GetSignal(i); 
        if( sig <= 0.0 ) break; // no more signals
        trk = pl->GetTrack(i);
        trk += fileIndex; 
        AddSignal(trk,pl->GetHit(i),module,index,sig);
    } // end for i
    fSignalAfterElect += (pl->fSignalAfterElect + pl->fNoise - fNoise);
    fNoise = pl->fNoise;
    return;
}
//______________________________________________________________________
Int_t AliITSpListItem::ShiftIndex(Int_t in,Int_t trk) const {
    // Shift an index number to occupy the upper four bits. No longer used.
    // Inputs:
    //    Int_t in   The file number
    //    Int_t trk  The track number
    // Outputs:
    //    none.
    // Return:
    //    Int_t The track number with the file number in the upper bits.
    Int_t si = sizeof(Int_t) * 8;
    UInt_t uin,utrk; // use UInt_t to avoid interger overflow-> goes negative.

    uin = in;
    utrk = trk;
    for(Int_t i=0;i<si-4;i++) uin *= 2;
    uin += utrk;
    in = uin;
    return in;
}
//______________________________________________________________________
void AliITSpListItem::Print(ostream *os) const {
    //Standard output format for this class
    // Inputs:
    //    ostream *os  The output stream
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;

    *os << fmodule <<","<<findex<<",";
    *os << fgksize <<",";
    for(i=0;i<fgksize;i++) *os << fTrack[i] <<",";
    for(i=0;i<fgksize;i++) *os << fHits[i] <<",";
    for(i=0;i<fgksize;i++) *os << fSignal[i] <<",";
    *os << fTsignal <<","<< fNoise << "," << fSignalAfterElect;
}
//______________________________________________________________________
void AliITSpListItem::Read(istream *is){
    // Standard output streaming function.
    // Inputs:
    //    istream *is The input stream
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,iss;

    *is >> fmodule >> findex;
    *is >> iss; // read in fgksize
    for(i=0;i<fgksize&&i<iss;i++) *is >> fTrack[i];
    for(i=0;i<fgksize&&i<iss;i++) *is >> fHits[i];
    for(i=0;i<fgksize&&i<iss;i++) *is >> fSignal[i];
    *is >> fTsignal >> fNoise >> fSignalAfterElect;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSpListItem &source){
    // Standard output streaming function.
    // Inputs:
    //    ostream &os             The output stream
    //    AliITSpListItem &source The AliITSpListItem object to be written out.
    // Outputs:
    //    none.
    // Return:
    //    ostream  The output stream

    source.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &os,AliITSpListItem &source){
    // Standard output streaming function.
    // Inputs:
    //    istream os              The input stream
    //    AliITSpListItem &source The AliITSpListItem object to be inputted
    // Outputs:
    //    none.
    // Return:
    //    istream The input stream.

    source.Read(&os);
    return os;
}
