#ifndef ALI_MUON_ST1_RESPONSE_H
#define ALI_MUON_ST1_RESPONSE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1Response
// ----------------------------
// Response class for station 1 including electronics and detector response. 
// Individual pedestals or noise levels can be controlled separately. 
// The current pulse height responses do not contain any physics

#include <TString.h>
#include <TList.h>
#include "AliMUONResponseV0.h"
#include "AliMUONSt1Types.h"
#include "AliMUONSt1ElectronicElement.h"

class AliMpPlane;
class AliMpPlaneSegmentation;
class AliMpZone;
class AliMpSector;
class TArrayF;
class TObjArray;
class AliMUONSt1ResponseParameter;

class AliMUONSt1Response : public AliMUONResponseV0 
{
public:
    AliMUONSt1Response(Int_t chamber=1);
    virtual ~AliMUONSt1Response();
    
    //
    // Configuration methods
    //
    void SetIniFileName(Int_t plane,const TString& fileName);

    virtual Float_t  IntPH(Float_t eloss);

    // Noise, zero-suppression, adc saturation
    virtual Int_t DigitResponse(Int_t digit,AliMUONTransientDigit* where);
    void PrintStatistics() const;


private:
    //private constants
    static const Int_t fgkNofZones=4;
    static const TString fgkTopDir;
    static const TString fgkDataDir;
    static const TString fgkConfigBaseName;
    static const TString fgkStandardIniFileName;    

    // static names
    static const TString fgkBaseName ;
    static const TString fgkIncludeName ;
    static const TString fgkParameterName ;
    static const TString fgkRegionName ;
    static const TString fgkRuleName ;
    static const TString fgkNameName ;
    static const TString fgkPedestalName ;
    static const TString fgkNoiseName ;
    static const TString fgkStateName ;
    static const TString fgkMName ;
    static const TString fgkMGName ;
    static const TString fgkMGCName ;
    static const TString fgkIJName ;
    static const TString fgkXYName ;
    static const TString fgkZoneName ;
    static const TString fgkStickyOnName ;
    static const TString fgkStickyOffName ;
    static const TString fgkFileName ;
    static const TString fgkValueName ;
    static const TString fgkGausName ;
    static const TString fgkNotName ;
    static const TString fgkNofSigmaName ;

    //protected methods
    AliMpZone* FindZone(AliMpSector* sector,Int_t posId); // to be moved in AliMpSector::
    void ReadFiles();
    void ReadIniFile(Int_t plane,const TString& fileName,Bool_t rdParam,Bool_t rdRegion,Bool_t rdRule);
    void ReadIniFile(Int_t plane);
    void ReadCouplesOfIntRanges(const string& value,TList* list,AliMUONSt1ElectronicElement::TDescription descr);
    void ReadCouplesOfFloatRanges(const string& value,TList* list);
    void SetPairToParam(const string& name,const string& value,AliMUONSt1ResponseParameter* param) const;
    void SetPairToListElem(const string& name,const string& value,TList* list);


    //data members
    AliMpPlane* fPlane[2];       // !The mapping planes
    AliMpPlaneSegmentation* fPlaneSegmentation[2]; // !The mapping plane segmentation
    TString fIniFileName[2];// file names for initialisation of each cathode

    AliMUONSt1ResponseParameter* fDefaultParameters[2][fgkNofZones]; // !Response for each zone
    TList fRulesList[2]; //! list of special rules

    Int_t fCountNofCalls;    // number of calls to DigitResponse()
    Int_t fCountUnknownZone; // ntimes the DigitResponse was called in an unknown zone
    Int_t fCountUnknownIndices; // ntimes the DigitResponse was called with unknown indices

    Int_t fChamber;                // The chamber number

    TParamsMap fParams;  //! internal parameter list
    TListMap   fRegions; //! internal list of regions
    TList      fTrashList; //!internal trash list 

  ClassDef(AliMUONSt1Response,1) // Overall detector response
};

#endif //ALI_MUON_ST1_RESPONSE_H
