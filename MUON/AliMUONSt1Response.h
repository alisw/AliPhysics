#ifndef ALI_MUON_ST1_RESPONSE_H
#define ALI_MUON_ST1_RESPONSE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1Response
// ----------------------------
// Response class for station 1 including electronics and detector response. 
// Individual pedestals or noise levels can be controlled separately. 
// The current pulse height responses do not contain any physics

#include <map>
#ifndef __HP_aCC
  using std::map;
#endif

#include <TString.h>
#include <TList.h>

#include "AliMUONResponseV0.h"
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

  protected:
    AliMUONSt1Response(const AliMUONSt1Response& rhs);

    // operators
    AliMUONSt1Response& operator=(const AliMUONSt1Response & rhs);

  private:
    // typedefs
    typedef map<string, AliMUONSt1ResponseParameter*> ParamsMap;
    typedef map<string, TList*>  ListMap;

    // private methods
    AliMpZone* FindZone(AliMpSector* sector,Int_t posId) const; // to be moved in AliMpSector::
    void ReadFiles();
    void ReadIniFile(Int_t plane,const TString& fileName,Bool_t rdParam,Bool_t rdRegion,Bool_t rdRule);
    void ReadIniFile(Int_t plane);
    void ReadCouplesOfIntRanges(const string& value,TList* list,AliMUONSt1ElectronicElement::TDescription descr);
    void ReadCouplesOfFloatRanges(const string& value,TList* list);
    void SetPairToParam(const string& name,const string& value,AliMUONSt1ResponseParameter* param) const;
    void SetPairToListElem(const string& name,const string& value,TList* list);

    // private constants
    static const Int_t fgkNofZones=4;           // number of zones
    static const TString fgkTopDir;             // top directory path
    static const TString fgkDataDir;            // data directory path
    static const TString fgkConfigBaseName;     // config file base name
    static const TString fgkStandardIniFileName;// standard ini file name    

    // static names
    static const TString fgkBaseName ;          // base name
    static const TString fgkIncludeName ;       // include name
    static const TString fgkParameterName ;     // parameter name
    static const TString fgkRegionName ;        // region name
    static const TString fgkRuleName ;          // rule name
    static const TString fgkNameName ;          // name name
    static const TString fgkPedestalName ;      // pedestal name
    static const TString fgkNoiseName ;         // noise name
    static const TString fgkStateName ;         // state name
    static const TString fgkMName ;             // M name
    static const TString fgkMGName ;            // MG name
    static const TString fgkMGCName ;           // MGC name
    static const TString fgkIJName ;            // i,j name
    static const TString fgkXYName ;            // x,y name
    static const TString fgkZoneName ;          // zone name
    static const TString fgkStickyOnName ;      // sticky on name
    static const TString fgkStickyOffName ;     // sticky off
    static const TString fgkFileName ;          // file name
    static const TString fgkValueName ;         // value name
    static const TString fgkGausName ;          // gauss name
    static const TString fgkNotName ;           // not name
    static const TString fgkNofSigmaName ;      // nof sigma name

    // data members
    AliMpPlane* fPlane[2];  // !The mapping planes
    AliMpPlaneSegmentation* fPlaneSegmentation[2]; // !The mapping plane segmentation
    TString fIniFileName[2];// file names for initialisation of each cathode

    AliMUONSt1ResponseParameter* fDefaultParameters[2][fgkNofZones]; // !Response for each zone
    TList fRulesList[2]; //! list of special rules

    Int_t fCountNofCalls;    // number of calls to DigitResponse()
    Int_t fCountUnknownZone; // ntimes the DigitResponse was called in an unknown zone
    Int_t fCountUnknownIndices; // ntimes the DigitResponse was called with unknown indices

    Int_t fChamber;                // The chamber number

    ParamsMap  fParams;    //! internal parameter list
    ListMap    fRegions;   //! internal list of regions
    TList      fTrashList; //! internal trash list 

  ClassDef(AliMUONSt1Response,1) // Overall detector response
};

#endif //ALI_MUON_ST1_RESPONSE_H
