#ifndef ALITRDALIGNMENT_H
#define ALITRDALIGNMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTRDalignment class is an instrument for reading, writing, and       //
// manipulating of the TRD alignment data.                                //
// D.Miskowiec, November 2006                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TRandom.h>
#include <TObjString.h>
#include <AliGeomManager.h>
class TClonesArray;
class AliSurveyObj;

class AliTRDalignment : public TObject {
  
 public:

  AliTRDalignment();
  AliTRDalignment(const AliTRDalignment& source);
  AliTRDalignment& operator=(const AliTRDalignment& source);  
  AliTRDalignment& operator*=(double fac);  
  AliTRDalignment& operator+=(const AliTRDalignment& source);  
  AliTRDalignment& operator-=(const AliTRDalignment& source);  
  Bool_t operator==(const AliTRDalignment& source) const;  
  virtual ~AliTRDalignment() {};

  // setters 

  void SetSmZero();                                  // reset to zero supermodule data
  void SetChZero();                                  // reset to zero chamber data
  void SetSm(int sm, const double a[6])              {for (int i = 0; i < 6; i++) fSm[sm][i] = a[i];}
  void SetCh(int ch, const double a[6])              {for (int i = 0; i < 6; i++) fCh[ch][i] = a[i];}
  void SetSmRandom(double a[6]);                     // generate random gaussians with sigmas a
  void SetChRandom(double a[6]);                     // generate random gaussians with sigmas a
  void SetSmFull();                                  // set supermodule data to initial aka full 
  void SetChFull();                                  // set chamber data to initial aka full 
  void SetSmResidual();                              // set supermodule data to final aka residual
  void SetChResidual();                              // set chamber data to final aka residual
  void SetZero()                                     {SetSmZero(); SetChZero();}
  void SetIdeal()                                    {SetZero();}
  void SetFull()                                     {SetSmFull(); SetChFull();}
  void SetResidual()                                 {SetSmResidual(); SetChResidual();}
  void SetComment(char *s)                           {fComment.SetString(s);} 

  // simple getters (for other getters see below)

  void GetSm(int sm, double * const a) const         {for (int i = 0; i < 6; i++) a[i] = fSm[sm][i];}
  void GetCh(int ch, double * const a) const         {for (int i = 0; i < 6; i++) a[i] = fCh[ch][i];}
  
  // dumping on screen

  void PrintSm(int sm, FILE * const fp = stdout) const;   // print data of a supermodule
  void PrintCh(int ch, FILE * const fp = stdout) const;   // print data of a chamber
  void PrintSm(FILE * const fp = stdout) const       {for (int i = 0; i <  18; i++)  PrintSm(i,fp);}
  void PrintCh(FILE * const fp = stdout) const       {for (int i = 0; i < 540; i++)  PrintCh(i,fp);}
  void Print(FILE * const fp = stdout) const         {PrintSm(fp); PrintCh(fp);                    }
  void Print(Option_t *) const                       {Print();                                     } 

  // reading-in from file

  void ReadAscii(const char * const filename);           // read from ascii file
  void ReadCurrentGeo();                                 // read from currently loaded geometry
  void ReadRoot(const char * const filename);            // read from root file
  void ReadDB(const char * const filename);              // read from DB file
  void ReadDB(const char * const db, const char * const path, int run, int version=-1, int subversion=-1);
  Bool_t DecodeSurveyPointName(TString pna, Int_t &sm, Int_t &iz,Int_t &ir, Int_t &iphi);
  void ReadSurveyReport(const char * const filename);    // read from survey report
  void ReadSurveyReport(const AliSurveyObj * const so);  // read from survey object 
  void ReadAny(const char * const  filename);            // read from any kind of file

  // writing on file

  void WriteAscii(const char * const filename) const;                    // store data on ascii file
  void WriteRoot(const char * const filename);                           // store data on root file
  void WriteDB(const char * const fina, int r0, int r1, int v, int s);   // store data on a local DB-like file
  void WriteDB(char * const db, const char * const pa, int r0, int r1);  // store data on DB file
  void WriteGeo(char *filename);                                         // apply misalignment and store geometry 

  // geometry and symbolic names getters

  // phi-sector number of chamber ch, 0-17
  int GetSec(int ch) const           {return ch/30;}
  // stack number, 0-4
  int GetSta(int ch) const           {return ch%30/6;}
  // plane number, 0-5 
  int GetPla(int ch) const           {return ch%30%6;}
  // module number, 0-89
  int GetMod(int ch) const           {return 5*GetSec(ch)+GetSta(ch);} 
  // layer number, 9-14
  int GetLay(int ch) const           {return AliGeomManager::kTRD1+GetPla(ch);}
  // volume id
  UShort_t GetVoi(int ch) const      {return AliGeomManager::LayerToVolUID(GetLay(ch),GetMod(ch));}
  // symbolic name of a supermodule
  char *GetSmName(int sm) const      {return Form("TRD/sm%02d",sm);}
  // symbolic name of a chamber
  char *GetChName(int ch) const      {return Form("TRD/sm%02d/st%d/pl%d",GetSec(ch),GetSta(ch),GetPla(ch));}
  // index of a supermodule 
  int GetSmIndex(const char *name)   {for (int i=0; i<18; i++) if (strcmp(name,GetSmName(i))==0) return i; return -1;}
  // index of a chamber
  int GetChIndex(const char *name)   {for (int i=0; i<540; i++) if (strcmp(name,GetChName(i))==0) return i; return -1;}

  // data analysis

  double GetSmRMS(int xyz) const;                    // calculate rms fSm[*][xyz]
  double GetChRMS(int xyz) const;                    // calculate rms fCh[*][xyz]
  void   PrintSmRMS() const;                         // print rms of fSm
  void   PrintChRMS() const;                         // print rms of fCh
  void   PrintRMS() const                            {PrintSmRMS(); PrintChRMS();}

  double SurveyChi2(int i, const double * const a);  // compare survey with ideal, return chi2
  double SurveyChi2(const double * const a)          {return SurveyChi2(fIbuffer[0],a);}
  void   SurveyToAlignment(int i, const char * const flag);  // determine alignment of supermodule i based on survey
  void   SurveyToAlignment(const char * const flag)  {for (int i=0; i<18; i++) SurveyToAlignment(i,flag);}

 protected:

  void   ArToNumbers(TClonesArray * const ar);       // read ar and fill fSm and fCh
  void   NumbersToAr(TClonesArray * const ar);       // build ar using fSm and fCh data
  int    IsGeoLoaded();                              // check if geometry is loaded

 protected:

  double     fSm[18][6];                             // supermodule data
  double     fCh[540][6];                            // chamber data 
  TObjString fComment;                               // info concerning origin of the data etc.
  TRandom    fRan;                                   // random generator for fake alignment data

  // Temporary storage for ideal position of the survey points and the survey data.
  // The survey data are in master frame and in cm. Each supermodule has 8 survey marks. 
  // The indices are sm number, z-end, radius, phi. 
  // The ideal positions of survey points are in local frame of supermodule and in cm. 
  // The indices are z-end, radius, phi. 
  // The processed survey results are stored in fSm.
  double fSurveyX[18][2][2][2];                      // supermodule survey point X
  double fSurveyY[18][2][2][2];                      // supermodule survey point Y
  double fSurveyZ[18][2][2][2];                      // supermodule survey point Z
  double fSurveyEX[18][2][2][2];                     // supermodule survey point X error
  double fSurveyEY[18][2][2][2];                     // supermodule survey point Y error
  double fSurveyEZ[18][2][2][2];                     // supermodule survey point Z error
  double fSurveyX0[2][2][2];                         // ideal X position of the survey marks
  double fSurveyY0[2][2][2];                         // ideal Y position of the survey marks
  double fSurveyZ0[2][2][2];                         // ideal Z position of the survey marks
  int    fIbuffer[1000];                             // generic buffer for misc. operations
  double fDbuffer[1000];                             // generic buffer for misc. operations

  ClassDef(AliTRDalignment,1)    

};
    
#endif
