#ifndef ALIFMDALIGNFAKER_H
#define ALIFMDALIGNFAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */

//____________________________________________________________________
//
//  Class to make fake alignration parameters 
//
#ifndef ROOT_TTask
# include <TTask.h>
#endif
#ifndef ROOT_TVector3
# include <TVector3.h>
#endif
class TClonesArray;
class TString;

class AliFMDAlignFaker : public TTask
{
public:
  enum EWhat {
    kSensors =  1, 
    kHalves
  };
  enum {
    kAll             = (1<<kHalves|1<<kSensors)
  };
  AliFMDAlignFaker(Int_t mask=kAll, 
		   const char* geo="geometry.root",
		   const char* loc="local://cdb");
  virtual ~AliFMDAlignFaker() {}
  void AddAlign(EWhat w) { SETBIT(fMask, w); }
  void RemoveAlign(EWhat w) { SETBIT(fMask, w); }
  void SetAlign(Int_t mask) { fMask = mask; }
  void SetSensorDisplacement(Double_t x1=0,   Double_t y1=0,   Double_t z1=0,
			     Double_t x2=.01, Double_t y2=.01, Double_t z2=0);
  void SetSensorRotation(Double_t x1=0,  Double_t y1=0,  Double_t z1=0,
			 Double_t x2=.5, Double_t y2=.5, Double_t z2=.5);
  void SetHalfDisplacement(Double_t x1=0,   Double_t y1=0,   Double_t z1=0,
			   Double_t x2=.05, Double_t y2=.05, Double_t z2=.05);
  void SetHalfRotation(Double_t x1=0, Double_t y1=0, Double_t z1=0,
		       Double_t x2=0, Double_t y2=0, Double_t z2=0);
  void SetOutput(const char* file) { SetTitle(file); }
  void SetGeometryFile(const char* file) { SetName(file); }
  void Exec(Option_t* option="");
protected:
  Bool_t MakeAlign(const TString& path, Int_t volID, 
		   Double_t transX, Double_t transY, Double_t transZ,
		   Double_t rotX, Double_t rotY, Double_t rotZ);
  Bool_t MakeAlignSensor(const TString& path, Int_t id);
  Bool_t MakeAlignHalf(const TString& path, Int_t id);
  void   WriteToCDB();
  void   WriteToFile();
  
  Long_t        fMask;            // What to write 
  TVector3      fSensorTransMin;
  TVector3      fSensorTransMax;
  TVector3      fSensorRotMin;
  TVector3      fSensorRotMax;
  TVector3      fHalfTransMin;
  TVector3      fHalfTransMax;
  TVector3      fHalfRotMin;
  TVector3      fHalfRotMax;
  Int_t         fRunMin;
  Int_t         fRunMax;
  TClonesArray* fArray;
  
  ClassDef(AliFMDAlignFaker,0)
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//

