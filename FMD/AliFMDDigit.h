// -*- mode: c++ -*- 
#ifndef ALIFMDDIGIT_H
#define ALIFMDDIGIT_H

//////////////////////////////////////////////////////////////////////
//
//  Digits classes for the FMD                
//
//////////////////////////////////////////////////////////////////////
#ifndef ROOT_TObject
# include <TObject.h>
#endif

//____________________________________________________________________
class AliFMDBaseDigit : public TObject 
{
protected:
  UShort_t fDetector;  // (Sub) Detector # (1,2, or 3)
  Char_t   fRing;      // Ring ID ('I' or 'O')
  UShort_t fSector;    // Sector # (phi division)
  UShort_t fStrip;     // Strip # (radial division)
public: 
  AliFMDBaseDigit();
  AliFMDBaseDigit(UShort_t detector, 
		  Char_t   ring='\0', 
		  UShort_t sector=0, 
		  UShort_t strip=0);
  virtual ~AliFMDBaseDigit() {}
  UShort_t     Detector()	   const { return fDetector; }
  Char_t       Ring()	           const { return fRing;     }
  UShort_t     Sector()	           const { return fSector;   }
  UShort_t     Strip()	           const { return fStrip;    }
  virtual void Print(Option_t* opt="") const;
  ClassDef(AliFMDBaseDigit, 1) // Base class for FMD digits 
};

//____________________________________________________________________
class AliFMDDigit : public AliFMDBaseDigit
{
protected:
  UShort_t fCount1;     // Digital signal 
  Short_t  fCount2;     // Digital signal (-1 if not used)
  Short_t  fCount3;     // Digital signal (-1 if not used)
public:
  AliFMDDigit();
  AliFMDDigit(UShort_t detector, 
	      Char_t   ring='\0', 
	      UShort_t sector=0, 
	      UShort_t strip=0, 
	      UShort_t count=0, 
	      Short_t  count2=-1, 
	      Short_t  count3=-1);
  virtual ~AliFMDDigit() {}
  UShort_t Count1()                const { return fCount1;   }
  Short_t  Count2()                const { return fCount2;   }
  Short_t  Count3()                const { return fCount3;   }
  UShort_t Counts()                const;
  void     Print(Option_t* opt="") const;
  ClassDef(AliFMDDigit,1)     // Normal FMD digit
};

inline UShort_t 
AliFMDDigit::Counts() const 
{
  return fCount1 
    + (fCount2 >= 0 ? fCount2 : 0)
    + (fCount3 >= 0 ? fCount3 : 0);
}

//____________________________________________________________________
class AliFMDSDigit : public AliFMDBaseDigit
{
protected:
  Float_t  fEdep;       // Energy deposited 
  UShort_t fCount1;     // Digital signal 
  Short_t  fCount2;     // Digital signal (-1 if not used)
  Short_t  fCount3;     // Digital signal (-1 if not used)
public:
  AliFMDSDigit();
  AliFMDSDigit(UShort_t detector, 
	       Char_t   ring='\0', 
	       UShort_t sector=0, 
	       UShort_t strip=0, 
	       Float_t  edep=0,
	       UShort_t count=0, 
	       Short_t  count2=-1, 
	       Short_t  count3=-1);
  virtual ~AliFMDSDigit() {}
  UShort_t Count1()                const { return fCount1;   }
  Short_t  Count2()                const { return fCount2;   }
  Short_t  Count3()                const { return fCount3;   }
  Float_t  Edep()                  const { return fEdep;     }
  UShort_t Counts()                const;
  void     Print(Option_t* opt="") const;
  ClassDef(AliFMDSDigit,1)     // Summable FMD digit
};
  
inline UShort_t 
AliFMDSDigit::Counts() const 
{
  return fCount1 
    + (fCount2 >= 0 ? fCount2 : 0)
    + (fCount3 >= 0 ? fCount3 : 0);
}


#endif
//____________________________________________________________________
//
// EOF
//
