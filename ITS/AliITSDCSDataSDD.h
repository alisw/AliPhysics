#ifndef ALIITSDCSDATASDD_H
#define ALIITSDCSDATASDD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
// Class to define object containing SDD DCS data                //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//         V.Pospisil, CTU Prague, gdermog@seznam.cz             //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<TArrayF.h>
#include<TArrayI.h>
#include<TArrayC.h>
#include<TArrayS.h>

class AliITSDCSDataSDD : public TObject
{

 public:
  AliITSDCSDataSDD( void );
                        // Default constructor
  ~AliITSDCSDataSDD( void ){};
                        // Destructor is void

/* There are allowed ranges of temperatures and medium voltages

     MV ....... 0.0 - 65.535 V
     TL, TR ... 0.0 - 655.35 C

   because these variables are stores in UShort_t arrays (it halves
   needed space). If value of any variable exceed allowed range,
   something very stupid would be stored.
*/

  void SetNPointsTempLeft( Int_t npts );
  void SetNPointsTempRight( Int_t npts );
  void SetNPointsHV( Int_t npts );
  void SetNPointsMV( Int_t npts );
  void SetNPointsStatus( Int_t npts );
                        // Sets sizes of the DCS variable arrays

  void SetValueTempLeft(Int_t time, Float_t temperature );
  void SetValueTempRight(Int_t time, Float_t temperature );
  void SetValueHV(Int_t time, Float_t voltage );
  void SetValueMV(Int_t time, Float_t voltage );
  void SetValueStatus(Int_t time, Char_t status );
                        // Inserts value of a DCS variable into the appropriate array.
                        //  Resizes and sorts array if necessary.

  void Compress();      // Minimize array sizes

  Int_t GetTempLeftRecords()  const {return fTempLeftSetPoints;}
  Int_t GetTempRightRecords() const {return fTempRightSetPoints;}
  Int_t GetHVRecords()        const {return fHVSetPoints;}
  Int_t GetMVRecords()        const {return fMVSetPoints;}
  Int_t GetStatusRecords()    const {return fStatusSetPoints;}
                       // Returns number of stored values of specific DCS variable

  Int_t GetTempLeftSize()  const {return fTempLeftMaxPoints;}
  Int_t GetTempRightSize() const {return fTempRightMaxPoints;}
  Int_t GetHVSize()        const {return fHVMaxPoints;}
  Int_t GetMVSize()        const {return fMVMaxPoints;}
  Int_t GetStatusSize()    const {return fStatusMaxPoints;}
                       // Returns size of selected array

  Float_t GetTempLeftIdx( Int_t index )  const { UShort_t val = (UShort_t)fTempLeft.At(index);
                                                 return (Float_t)val / fgkTPrec; };
  Float_t GetTempRightIdx( Int_t index ) const { UShort_t val = (UShort_t)fTempRight.At(index);
                                                 return (Float_t)val / fgkTPrec;};
  Float_t GetHVIdx( Int_t index )        const { return fHV.At(index);};
  Float_t GetMVIdx( Int_t index )        const { UShort_t val = (UShort_t)fMV.At(index);
                                                 return (Float_t)val / fgkMVPrec;};
  Char_t  GetStatusIdx( Int_t index )    const { return fStatus.At(index);};
                       // Returns value of specific DCS variable by index in the array

  Bool_t GetOKStatIdx( Int_t index )        const { return (Bool_t)(fStatus.At(index) & 1 ); };
  Bool_t GetTempLeftStatIdx( Int_t index )  const { return (Bool_t)(fStatus.At(index) & 2 ); };
  Bool_t GetTempRightStatIdx( Int_t index ) const { return (Bool_t)(fStatus.At(index) & 4 ); };
                       // Return status of a readout device by index in the array

  Int_t   GetTempLeftTimeIdx( Int_t index )  const {return fTempLeftTimeStamp.At(index);};
  Int_t   GetTempRightTimeIdx( Int_t index ) const {return fTempRightTimeStamp.At(index);};
  Int_t   GetHVTimeIdx( Int_t index )        const {return fHVTimeStamp.At(index);};
  Int_t   GetMVTimeIdx( Int_t index )        const {return fMVTimeStamp.At(index);};
  Int_t   GetStatusTimeIdx( Int_t index )    const {return fStatusTimeStamp.At(index);};
                       // Returns time stamp of specific DCS variable by index in the array

  Float_t GetTempLeft( Int_t time )  const { Int_t i = FindIndex( time, fTempLeftTimeStamp, fTempLeftSetPoints );
                                             return ( i < 0 ) ? -1.0 : GetTempLeftIdx( i ); }
  Float_t GetTempRight( Int_t time ) const { Int_t i = FindIndex( time, fTempRightTimeStamp, fTempRightSetPoints );
                                             return ( i < 0 ) ? -1.0 : GetTempRightIdx( i ); }
  Float_t GetHV( Int_t time )        const { Int_t i = FindIndex( time, fHVTimeStamp, fHVSetPoints );
                                             return ( i < 0 ) ? -1.0 : GetHVIdx( i ); }
  Float_t GetMV( Int_t time )        const { Int_t i = FindIndex( time, fMVTimeStamp, fMVSetPoints );
                                             return ( i < 0 ) ? -1.0 : GetMVIdx( i ); }
  Char_t  GetStatus( Int_t time )    const { Int_t i = FindIndex( time, fStatusTimeStamp, fStatusSetPoints );
                                             return ( i < 0 ) ? -1 : GetStatusIdx( i ); }
                       // Returns value of specific DCS variable by time stamp

  Bool_t GetOKStat( Int_t time )        const { return (Bool_t)( GetStatus( time ) & 1 ); };
  Bool_t GetTempLeftStat( Int_t time )  const { return (Bool_t)( GetStatus( time ) & 2 ); };
  Bool_t GetTempRightStat( Int_t time ) const { return (Bool_t)( GetStatus( time ) & 4 ); };
                       // Return status of a readout device in given time

  Float_t GetDriftField( Int_t timeStamp ) const;
                       // Returns drift field counted for specific time

  Float_t GetDriftSpeed( Int_t /*timeStamp*/ ) const;  /* --- DUMMY --- */
                       // Returns drift speed counted for specific time. Calculation is based on temerature 
                       //  taken  from DCS. This metod is not dedicated for normal usage, it should be used
                       //  only in cases that the injectors for given module fails. 
                       //
                       // Presently only a prototype, returns -1.0.

  void PrintValues( FILE *output = stdout ) const;
                       // Displays stored DCS varable values or writes it into a text file
 private:

  TArrayS fTempLeft;           // Temperature on left side. If there is stored a negative value
                               //  something wrong happend with the temperature chip.
                               //  Temperatures and medium voltages are stored as UShort_t, which 
                               //  takes half memory as Float_t. It makes ranges 
                               //
                               //  MV ....... 0.0 - 65.535 volts
                               //  TL, TR ... 0.0 - 655.35 C
                               //
                               //  which should be enough.
                               //
  TArrayI fTempLeftTimeStamp;  // Time stamps of the temperatures on left side
  Int_t   fTempLeftMaxPoints;  // Size of the arrays
  Int_t   fTempLeftSetPoints;  // Number of filled array cells (number of set values)


  TArrayS fTempRight;          // Temperature on right side. If there is stored a negative value
                               //  something wrong happend with the temperature chip.
  TArrayI fTempRightTimeStamp; // Time stamps of temperatures on right side
  Int_t   fTempRightMaxPoints; // Size of the arrays
  Int_t   fTempRightSetPoints; // Number of filled array cells (number of set values)


  TArrayF fHV;                 // High voltage on SDD
  TArrayI fHVTimeStamp;        // Time stamps of HV
  Int_t   fHVMaxPoints;        // Size of the arrays
  Int_t   fHVSetPoints;        // Number of filled array cells (number of set values)


  TArrayS fMV;                 // Medium voltage on SDD
  TArrayI fMVTimeStamp;        // Time stamps of MV
  Int_t   fMVMaxPoints;        // Size of the arrays
  Int_t   fMVSetPoints;        // Number of filled array cells (number of set values)

  TArrayC fStatus;             // Status of temperature and voltage readout
                               //
                               // bit 0 - _OK
                               // bit 1 - _TEMP_L_STATE
                               // bit 2 - _TEMP_R_STATE

  TArrayI fStatusTimeStamp;    // Time stamps of MV
  Int_t   fStatusMaxPoints;    // Size of the arrays
  Int_t   fStatusSetPoints;    // Number of filled array cells (number of set values)


  static const Float_t fgkTPrec;  // Number of temperature decimal places stored
  static const Float_t fgkMVPrec; // Number of medium voltage decimal places stored
                               //  There are three possibilities :
                               //  10.0 ..... one decimal place
                               //  100.0 .... two decimal places
                               //  1000.0 ... three decimal places
                               //  Values are set in AliITSDCSDataSDD.cxx, by default 
                               //  it is fgkTPrec = 100.0, fgkMVPrec = 1000.0

  Int_t   FindIndex( Int_t timeStamp, const TArrayI &timeStampArray, Int_t n ) const;
                        // Provides binary search in the time array. Returns index in the array of time 
                        //  stamps by selected value. Returns -1 if the time is less than time stamp in 
                        //  the timeArray[0]

  ClassDef(AliITSDCSDataSDD, 3)

}; /*class AliITSDCSDataSDD*/

#endif
