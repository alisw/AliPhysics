// -*- mode: C++ -*- 
#ifndef ALIAODMCHEADER_H
#define ALIAODMCHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Class AliAODMCHeader
//   Some MC specific inforamtion for filtering KINE infomration to the AOD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

class AliGenEventHeader;

#include "AliVHeader.h"
#include "TList.h"

class AliAODMCHeader: public AliVHeader {
public:
  AliAODMCHeader();
  virtual ~AliAODMCHeader();
  AliAODMCHeader(const AliAODMCHeader& header);
  AliAODMCHeader& operator=(const AliAODMCHeader& header);
  virtual void Copy(TObject &obj) const;

  virtual void      SetGeneratorName(const char* c){fGenerator = c;}
  virtual void      AddGeneratorName(const char* c);
  virtual const char* GetGeneratorName() const {return fGenerator.Data();}

  virtual void SetVertex(Double_t *vtx){
    fVertex[0] = vtx[0]; fVertex[1] = vtx[1]; fVertex[2] = vtx[2];
  }
  virtual void SetVertex(Double_t x,Double_t y,Double_t z){
    fVertex[0] = x; fVertex[1] = y; fVertex[2] = z;
  }
  virtual void GetVertex(Double_t *vtx) const {
    vtx[0] = fVertex[0]; vtx[1] = fVertex[1]; vtx[2] = fVertex[2];
  }

  virtual Double_t GetVtxX() const { return fVertex[0]; }
  virtual Double_t GetVtxY() const { return fVertex[1]; }
  virtual Double_t GetVtxZ() const { return fVertex[2]; }

  
  virtual void      SetImpactParameter(Double_t b){fImpactPar = b;}
  virtual Double_t  GetImpactParameter() const {return fImpactPar;}

  virtual void      SetPtHard(Double_t f){fPtHard = f;}
  virtual Double_t  GetPtHard() const {return fPtHard;}

  virtual void      SetCrossSection(Double_t f){fXsection = f;}
  virtual Double_t  GetCrossSection() const {return fXsection;}

  virtual void      AddTrial(Int_t i) {fTrials+=i;}
  virtual void      SetTrials(Int_t f){fTrials = f;}
  virtual Int_t     GetTrials() const {return fTrials;}

  virtual void      SetReactionPlaneAngle(Double_t b){fReactionPlaneAngle = b;}
  virtual Double_t  GetReactionPlaneAngle() const {return fReactionPlaneAngle;}

  virtual void      SetEventType(UInt_t eventType){fEventType = eventType;}
  virtual UInt_t    GetEventType() const {return fEventType;}

  virtual void      Reset();
  virtual void      Print(const Option_t *opt=0) const;

  // needed to make class non virtual
  virtual UShort_t  GetBunchCrossNumber()   const {return 0;}
  virtual UInt_t    GetOrbitNumber()        const {return 0;}
  virtual UInt_t    GetPeriodNumber()       const {return 0;}
  virtual ULong64_t GetTriggerMask()        const {return 0;}
  virtual UChar_t   GetTriggerCluster()     const {return 0;}
  virtual UInt_t    GetTimeStamp()          const {return 0;}
  // 
  
  // Access to header informations

  virtual void AddCocktailHeader(const AliGenEventHeader* header);
  virtual void AddCocktailHeaders(AliGenEventHeader* header);
  virtual AliGenEventHeader* GetCocktailHeader(Int_t i);
  virtual TList* GetCocktailHeaders(){return fHeaders;}
  virtual UInt_t GetNCocktailHeaders(){
    if(fHeaders)return fHeaders->GetEntries();
    return 0;
  }

  static const char* StdBranchName(){return fgkStdBranchName.Data();}

private:

  static TString fgkStdBranchName;      // Standard branch name

  // General event information

  TString      fGenerator;         // Name of the generator, combination of names in case of gen cocktail 
  Double32_t   fVertex[3];         // MC vertex
  Double32_t   fImpactPar;         // Impact parameter in case of Pb+Pb
  Double32_t   fPtHard;            // [0,0,12] Pt hard for jet events
  Double32_t   fXsection;          // Cross section for particlar process
  UInt_t       fTrials;            // Number of trials
  UInt_t       fEventType;         // MC Process Type of Event
  Double32_t   fReactionPlaneAngle;// MC Reaction Plane Angle

  // more details in the headers
  TList  *fHeaders;                // List of all MC Headers 

  ClassDef(AliAODMCHeader,6)

};

#endif
