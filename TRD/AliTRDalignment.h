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

#include <AliAlignObj.h>

class AliTRDalignment : public TObject {
  
 public:

  AliTRDalignment();
  AliTRDalignment(const AliTRDalignment& source);
  AliTRDalignment& operator=(const AliTRDalignment& source);  
  AliTRDalignment& operator+=(const AliTRDalignment& source);  
  AliTRDalignment& operator-=(const AliTRDalignment& source);  
  Bool_t operator==(const AliTRDalignment& source) const;  
  virtual ~AliTRDalignment() {};                     // destructor

  // setting 

  void SetSmZero();                                  // reset to zero supermodule data
  void SetChZero();                                  // reset to zero chamber data
  void SetZero() {SetSmZero(); SetChZero();}         // reset to zero both
  void SetSm(Int_t sm, const Double_t a[6])          { for (int i = 0; i < 6; i++) fSm[sm][i] = a[i]; }
  void SetCh(Int_t ch, const Double_t a[6])          { for (int i = 0; i < 6; i++) fCh[ch][i] = a[i]; }
  void SetSmRandom(Double_t a[6]);                   // generate random gaussians with sigmas a
  void SetChRandom(Double_t a[6]);                   // generate random gaussians with sigmas a
  void SetSmFull();                                  // set supermodule data to initial aka full 
  void SetChFull();                                  // set chamber data to initial aka full 
  void SetSmResidual();                              // set supermodule data to final aka residual
  void SetChResidual();                              // set chamber data to final aka residual
  void SetFull() {SetSmFull(); SetChFull();}
  void SetResidual() {SetSmResidual(); SetChResidual();}

  // dumping on screen

  void PrintSm(Int_t sm, FILE *fp = stdout) const;   // print data of a supermodule
  void PrintCh(Int_t sm, FILE *fp = stdout) const;   // print data of a chamber
  void PrintSm(FILE *fp = stdout) const              { for (int i = 0; i <  18; i++)  PrintSm(i,fp);  }
  void PrintCh(FILE *fp = stdout) const              { for (int i = 0; i < 540; i++)  PrintCh(i,fp);  }
  void Print(FILE *fp = stdout) const                { PrintSm(fp); PrintCh(fp);                      }
  void Print(Option_t *) const                       { Print();                                       } 

  // reading-in from file

  void ReadAscii(char *filename);                    // read from ascii file
  void ReadRoot(char *filename);                     // read from root file
  void ReadDB(char *filename);                       // read from DB file
  void ReadDB(char *db, char *path, Int_t run, Int_t version=-1, Int_t subversion=-1);
  void ReadGeo(char *misaligned);                    // read from misaligned_geometry.root
  void ReadSurveyReport(char *filename);             // read from survey report
  void ReadAny(char *filename);                      // read from any kind of file

  // writing on file

  void WriteAscii(char *filename) const;             // store data on ascii file
  void WriteRoot(char *filename);                    // store data on root file
  void WriteDB(char *filename, char *comment, Int_t run0, Int_t run1);
                                                     // store data on a local DB-like file
  void WriteDB(char *db, char *path, char *comment, Int_t run0, Int_t run1); 
                                                     // store data on DB file
  void WriteGeo(char *filename);                     // apply misalignment and store geometry 

  // geometry and symbolic names getters

  // phi-sector number of chamber ch, 0-17
  Int_t    GetSec(Int_t ch) const                    { return ch/30;   }
  // stack number, 0-4
  Int_t    GetSta(Int_t ch) const                    { return ch%30/6; }
  // plane number, 0-5
  Int_t    GetPla(Int_t ch) const                    { return ch%30%6; }
  // module number, 0-89
  Int_t    GetMod(Int_t ch) const                    { return 5*GetSec(ch)+GetSta(ch);                           } 
  // layer number, 9-14
  Int_t    GetLay(Int_t ch) const                    { return AliAlignObj::kTRD1+GetPla(ch);                     }
  // volume id
  UShort_t GetVoi(Int_t ch) const                    { return AliAlignObj::LayerToVolUID(GetLay(ch),GetMod(ch)); }
  char    *GetSmName(Int_t sm) const                 { return Form("TRD/sm%02d",sm);                             }
  char    *GetChName(Int_t ch) const                 { return Form("TRD/sm%02d/st%d/pl%d",GetSec(ch),GetSta(ch),GetPla(ch)); }

  // data analysis

  Double_t GetSmRMS(Int_t xyz) const;                // calculate rms fSm[*][xyz]
  Double_t GetChRMS(Int_t xyz) const;                // calculate rms fCh[*][xyz]
  void     PrintSmRMS() const;                       // print rms of fSm
  void     PrintChRMS() const;                       // print rms of fCh
  void     PrintRMS() const                          { PrintSmRMS(); PrintChRMS();}

 protected:

  void ArToNumbers(TClonesArray *ar);                // read ar and fill fSm and fCh
  void NumbersToAr(TClonesArray *ar);                // build ar using fSm and fCh data
  void LoadIdealGeometry(char *filename);            // load ideal geometry from file
  void LoadIdealGeometry() {LoadIdealGeometry("geometry.root");} 

 protected:

  Double_t fSm[18][6];                               // supermodule data
  Double_t fCh[540][6];                              // chamber data 
  TRandom  fRan;                                     // random generator for fake alignment data

  ClassDef(AliTRDalignment,1)    

};
    
#endif
