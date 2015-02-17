#ifndef AliDCSSensor_H
#define AliDCSSensor_H
/* Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



////////////////////////////////////////////////////////////////////////////
//              Container class for DCS sensors            
////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "AliSplineFit.h"

class TObjArray;
class TClonesArray;
class TTree;
class TFile;
class TString;
class TGraph;
class TVector3;


////////////////////////////////////////////////////////////////////////
//              Class AliDCSSensor
////////////////////////////////////////////////////////////////////////

class AliDCSSensor : public TNamed {

public:
  AliDCSSensor();
  AliDCSSensor(const AliDCSSensor& source);
  virtual ~AliDCSSensor();
  AliDCSSensor& operator=(const AliDCSSensor& source);
  virtual void Print(const Option_t* option="") const;
  virtual void Draw(Option_t* option="") ;
  
  Int_t       GetId()     const {return fId;     }
  Int_t       GetIdDCS()  const {return fIdDCS;     }
  const TString& GetStringID() const {return fStringID; }

  Double_t    GetX()   	   const {return fX;      }
  Double_t    GetY()   	   const {return fY;      }
  Double_t    GetZ()   	   const {return fZ;      }
  Double_t    GetR()   	   const {return TMath::Sqrt(fX*fX+fY*fY);}
  Double_t    GetPhi() 	   const {return TMath::ATan2(fY,fX);    }
  
  UInt_t      GetStartTime() const {return fStartTime;} 
  UInt_t      GetEndTime() const { return fEndTime; }
  TGraph*     GetGraph()   const {return fGraph; }
  AliSplineFit* GetFit()   const {return fFit; }

  void SetId     (Int_t id)        {fId     = id;    }
  void SetIdDCS  (Int_t iddcs)     {fIdDCS  = iddcs;    }
  void SetStringID (const TString& stringID)  {fStringID = stringID; }
  
  void SetX       (Double_t x)       {fX     = x;      }
  void SetY       (Double_t y)       {fY     = y;      }
  void SetZ       (Double_t z)       {fZ     = z;      }

  void SetGraph   (TGraph *gr)      {fGraph = gr; }
  void SetFit     (AliSplineFit *f) {fFit = f; }
  void SetStartTime (UInt_t stime)  {fStartTime = stime; }
  void SetStartTime (TTimeStamp time)  {fStartTime = time.GetSec(); }
  void SetEndTime (UInt_t stime)  {fEndTime = stime; }
  void SetEndTime (TTimeStamp time)  {fEndTime = time.GetSec(); }
  Double_t GetValue(UInt_t timeSec);   // timeSec offset for individual sensor map
  Double_t GetValue(TTimeStamp time);  // use this for absolute time stamp
  Double_t Eval(const TTimeStamp& time, Bool_t& inside) const;
  Double_t EvalGraph(const TTimeStamp& time, Bool_t& inside) const;
  Double_t EvalGraph(const Double_t& timeHour) const;
  TGraph *MakeGraph (Int_t nPoints=100, Bool_t debug=kFALSE) const;
  static TClonesArray *  ReadTree(TTree *tree);
  

protected:
  Int_t fId;         // Internal number of sensor id  (array index)
  Int_t fIdDCS;      // ID number in DCS
  TString  fStringID; // Amanda String ID
  UInt_t   fStartTime;  // start time for DCS map/fit
  UInt_t   fEndTime;    // end time for DCS map/fit
  TGraph * fGraph;      // graph with values
  AliSplineFit *fFit;   // fit to graph values
  Double_t fX;      //X-position of the sensor
  Double_t fY;      //Y-position of the sensor
  Double_t fZ;      //Z-position of the sensor

  ClassDef(AliDCSSensor,3);
};
#endif

