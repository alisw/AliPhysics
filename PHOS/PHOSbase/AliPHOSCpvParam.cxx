/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

// Class provides a set of static methods to convert absolute number of pad to pair (X,Y) 
// and vice versa
// and some other
// Author - Mikhail Stolpovskiy, IHEP Protvino (2013)

#include "AliPHOSCpvParam.h" 
// #include "AliPHOSCpvRawStream.h"

ClassImp(AliPHOSCpvParam);

using namespace std;

//====================================================================================
AliPHOSCpv3GConnection AliPHOSCpvParam::fConnection;
//====================================================================================
Bool_t AliPHOSCpvParam::DecodeRawWord(Int_t ddl,Int_t rWord, Int_t & abs, Int_t & q, Int_t & eType) {

  //std::cout<<"ddl = "<<ddl<<", mod = "<<DDL2Mod(ddl)<<std::endl;
  
  if (((rWord >> 27) & 1)) { //check if it's end_of_event word for 3gassiplex card.
    eType = -1;
    return kFALSE;
  }

  //std::cout<<"AliPHOSCpvParam::DecodeRawWord(): passed ((rWord >> 27) & 1)"<<std::endl;

  UInt_t statusControlRow = 0x32a8;
  UInt_t rowControlWord = rWord & 0xfbff;
  if(rowControlWord == statusControlRow) {
    eType = -1;
    return kFALSE;
  }

  //std::cout<<"AliPHOSCpvParam::DecodeRawWord(): passed rowControlWord == statusControlRow"<<std::endl;


  abs = rWord>>12;
  q = rWord & 0xfff;
  abs |= (DDL2Mod(ddl))<<15;
  if(!IsValidAbs(abs)) {
    eType = -1;
    return kFALSE;
  }
  //std::cout<<"AliPHOSCpvParam::DecodeRawWord(): passed IsValidAbs(abs)"<<std::endl;
  if(A2Pad(abs) >= kNPadAdd) {
    // eType = AliPHOSCpvRawStream::kWrongPadErr; // Eliminate circular dependence AliPHOSCpvParam <-> AliPHOSCpvRawStream!
    eType = 4;
    return kFALSE;
  }

  //Printf("AliPHOSCpvParam::DecodeRawWord: mod, cc, 3g, pad , q = %d, %d, %d, %d, %d",A2Mod(abs),A2CC(abs),A23G(abs),A2Pad(abs), q);

  return kTRUE;
}
//====================================================================================
Int_t AliPHOSCpvParam::Abs(Int_t ddl,Int_t columnCtrl,Int_t gassiplex3,Int_t pad) {
   if(ddl<0 || ddl>=kNDDL ||
     columnCtrl<0 || columnCtrl>=kNRows  ||
     gassiplex3<0 || gassiplex3>=kN3GAdd ||
     pad<0        || pad>=kNPadAdd) return -1;
   Int_t module = DDL2Mod(ddl); // module is number of CPV module (1-5) accordance to the PHOS offline 
   if(module == -1) return -1;
   return module<<15 
     | (columnCtrl+1)<<10
     | (gassiplex3+1)<<6
     | pad;
}
//====================================================================================
Bool_t AliPHOSCpvParam::IsValidAbs(Int_t abs) {
  //  Printf("AliPHOSCpvParam::IsValidAbs: abs = %d",abs);
  Int_t mod = A2Mod(abs),
    cc = A2CC(abs),
    g3 = A23G(abs),
    pad = A2Pad(abs);
  //  Printf("mod = %d, cc = %d, g3 = %d, pad = %d",mod,cc,g3,pad);
  if(mod<1 || mod>kNModules ||
     cc <0 || cc >=kNRows ||
     g3 <0 || g3 >=kN3GAdd ||
     pad<0 || pad>=kNPadAdd) return kFALSE;
  return kTRUE;
}
//====================================================================================
Int_t AliPHOSCpvParam::A2DDL(Int_t abs) { return Mod2DDL(A2Mod(abs));}
//====================================================================================
Int_t AliPHOSCpvParam::A2Mod(Int_t abs) { /*cout << "module is" << (abs>>15) << endl;*/ return abs>>15; }
//====================================================================================
Int_t AliPHOSCpvParam::DDL2Mod(Int_t ddl) {
  switch(ddl) {
  case (0) : return 5; break;
  case (2) : return 4; break;
  case (4) : return 3; break;
  case (6) : return 2; break;
  case (8) : return 1; break;
  default : return -1; break;
  }
}
//====================================================================================
Int_t AliPHOSCpvParam::Mod2DDL(Int_t mod) {
  switch(mod) {
  case (1) : return 8; break;
  case (2) : return 6; break;
  case (3) : return 4; break;
  case (4) : return 2; break;
  case (5) : return 0; break;
  default : return -1; break;
  }
}
//====================================================================================
Int_t AliPHOSCpvParam::A2CC (Int_t abs) { return ((abs >> 10) & 0x1f) - 1; }
//====================================================================================
Int_t AliPHOSCpvParam::A23G (Int_t abs) { return ((abs >> 6 ) & 0xf) - 1; }
//====================================================================================
Int_t AliPHOSCpvParam::A2Pad(Int_t abs) { return  abs & 0x3f; }
//====================================================================================
Int_t AliPHOSCpvParam::A2X  (Int_t abs) {
  if(!IsValidAbs(abs)) {
    //    Printf("AliPHOSCpvParam::A2X: abs is not valid!");
    return -1;
  }
  return (kNRows - 1 - A2CC(abs))*(kPadPcX/kNRows) + ( fConnection.Pad2X(A2Pad(abs)));
}
//====================================================================================
Int_t AliPHOSCpvParam::A2Y  (Int_t abs) {
  if(!IsValidAbs(abs)) {
    //    Printf("AliPHOSCpvParam::A2Y: abs is not valid!");
    return -1;
  }
  //return A23G(abs)*(kPadPcY/kN3GAdd) + connection.Ch2Y(A2Pad(abs));
  //return (kN3GAdd - 1 - A23G(abs))*(kPadPcY/kN3GAdd) + (kPadPcY/kN3GAdd - 1 - fConnection.pad2Y(A2Pad(abs)));
  return (A23G(abs))*(kPadPcY/kN3GAdd) + (5-fConnection.Pad2Y(A2Pad(abs)));
  //return (kN3GAdd - 1 - A23G(abs))*(kPadPcY/kN3GAdd) + connection.Ch2Y(A2Pad(abs));
}
//====================================================================================
Int_t AliPHOSCpvParam::XY2A (Int_t ddl, Int_t x, Int_t y) {
  if(x<kMinPx || x>kMaxPx || y<kMinPy || y>kMaxPy) return -1;
  return Abs(ddl,X2CC(x),Y23G(y),XY2Pad(x,y));
} // XY2A
//====================================================================================
Int_t AliPHOSCpvParam::X2CC (Int_t x) {
  if(x<kMinPx|| x>kMaxPx) return -1;
  return kNRows - 1 - x/(kPadPcX/kNRows);
} // X2CC
//====================================================================================
Int_t AliPHOSCpvParam::Y23G (Int_t y) {
  if(y<kMinPy || y>kMaxPy) return -1;
  return y/(kPadPcY/kN3GAdd);
} // Y23G
//====================================================================================
Int_t AliPHOSCpvParam::XY2Pad(Int_t x, Int_t y) {
  if(x<kMinPx || x>kMaxPx || y<kMinPy || y>kMaxPy) return -1;
  Int_t xPad = x - (kNRows - 1 - X2CC(x))*(kPadPcX/kNRows);
  //Int_t yPad = y - (kN3GAdd- 1 - Y23G(y))*(kPadPcY/kN3GAdd);
  Int_t yPad = y - (Y23G(y))*(kPadPcY/kN3GAdd);
  //return connection.XY2Ch(xPad,yPad);
  //return fConnection.XY2pad(xPad,kPadPcY/kN3GAdd - 1 - yPad);
 return fConnection.XY2Pad(xPad,5-yPad);
} // XY2Pad 
//====================================================================================
Bool_t AliPHOSCpvParam::GetLimOfCConX( Int_t cc, Int_t &xmin, Int_t &xmax) {
  //cout<<"cc="<<cc;
  if(cc < 0 || cc > kNRows) return kFALSE;
  Int_t a1 = Abs(0,cc,0,fConnection.XY2Pad(0,0)),
        a2 = Abs(0,cc,0,fConnection.XY2Pad(7,0));
  if(!(IsValidAbs(a1) && IsValidAbs(a2))) return kFALSE;
  xmin = A2X(a1);
  xmax = A2X(a2);
  //cout<<": xmin = "<<xmin<<" xmax = "<<xmax<<endl;
  return kTRUE;
} // GetLimOfCConX
//====================================================================================
Bool_t AliPHOSCpvParam::GetLimOf3GonY( Int_t g3, Int_t &ymin, Int_t &ymax) {
  if(g3 < 0 || g3 > kN3GAdd) return kFALSE;
  Int_t a1 = Abs(0,0,g3,fConnection.XY2Pad(0,0)),
        a2 = Abs(0,0,g3,fConnection.XY2Pad(0,5)); 
  if(!(IsValidAbs(a1) && IsValidAbs(a2))) return kFALSE;
  ymin = A2Y(a1);
  ymax = A2Y(a2);
  return kTRUE;
} // GetLimOf3GonY
//====================================================================================
Int_t AliPHOSCpvParam::A2fId(Int_t abs) {
  Int_t fId = 17920;
  fId += kPadPcX*kPadPcY*(A2DDL(abs)-1);
  fId += (kPadPcX/kNRows )*A2CC(abs);
  fId += (kPadPcY/kN3GAdd)*A23G(abs);
  fId += A2Pad(abs);
  return fId;
} // A2fId
