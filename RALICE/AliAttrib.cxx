/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliAttrib
// Generic handling of detector signal (calibration) attributes.
// Normally this class is only used as a base class to provide the various
// attributes to the derived class. An example of this is AliSignal.
// However, one can of course also use this class on its own as shown
// in the simple example hereafter.
//
// Example :
// ---------
// AliAttrib a;
// a.SetSlotName("PMT amplitude in Volt");
// a.SetGain(250.7);
// a.SetSlotName("Time of flight in ns",2);
// a.SetOffset(-22.5,2);
// a.SetSlotName("PMT amplitude in ADC",3);
// a.SetGain(1340,3);
// a.SetSlotName("TDC",4);
// a.SetOffset(10.75,"TDC");
// a.SetEdgeOn(3);
// a.SetDead(1);
// a.List();
//
//--- Author: Nick van Eijndhoven 18-sep-2003 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliAttrib.h"
#include "Riostream.h"
 
ClassImp(AliAttrib) // Class implementation to enable ROOT I/O
 
AliAttrib::AliAttrib()
{
// Creation of an AliAttrib object and initialisation of parameters.
// Several values of the same type (e.g. gain) can be stored in different slots.
// If needed, the storage for values will be expanded automatically
// when entering values.
 fGains=0;
 fOffsets=0;
 fCalflags=0;
 fNames=0;
}
///////////////////////////////////////////////////////////////////////////
AliAttrib::~AliAttrib()
{
// Destructor to delete dynamically allocated memory
 if (fGains)
 {
  delete fGains;
  fGains=0;
 }
 if (fOffsets)
 {
  delete fOffsets;
  fOffsets=0;
 }
 if (fCalflags)
 {
  delete fCalflags;
  fCalflags=0;
 }
 if (fNames)
 {
  delete fNames;
  fNames=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliAttrib::AliAttrib(const AliAttrib& a)
{
// Copy constructor
 fGains=0;
 fOffsets=0;
 fCalflags=0;
 fNames=0;

 Int_t n=0;
 Double_t val=0;

 n=a.GetNgains();
 for (Int_t ig=1; ig<=n; ig++)
 {
  val=a.GetGain(ig);
  if (a.GetGainFlag(ig)) SetGain(val,ig);
 }

 n=a.GetNoffsets();
 for (Int_t io=1; io<=n; io++)
 {
  val=a.GetOffset(io);
  if (a.GetOffsetFlag(io)) SetOffset(val,io);
 }

 n=a.GetNcalflags();
 for (Int_t ic=1; ic<=n; ic++)
 {
  SetEdgeValue(a.GetEdgeValue(ic),ic);
  if (a.GetDeadValue(ic)) SetDead(ic);
 }

 n=a.GetNnames();
 TString s;
 for (Int_t in=1; in<=n; in++)
 {
  s=a.GetSlotName(in);
  if (s!="") SetSlotName(s,in);
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetNgains() const
{
// Provide the number of specified gains for this attribute.
 Int_t n=0;
 if (fGains) n=fGains->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetNoffsets() const
{
// Provide the number of specified offsets for this attribute.
 Int_t n=0;
 if (fOffsets) n=fOffsets->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetNcalflags() const
{
// Provide the number of specified calib. flags for this attribute.
 Int_t n=0;
 if (fCalflags) n=fCalflags->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetNnames() const
{
// Provide the maximum number of specified names for this attribute.
 Int_t n=0;
 if (fNames) n=fNames->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetGain(Double_t gain,Int_t j)
{
// Store gain value of the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for gain values, the number of reserved slots for the gain
// values is increased automatically.

 if (j<1) 
 {
  cout << " *AliAttrib::SetGain* Invalid argument j = " << j << endl;
  return;
 }

 if (!fGains)
 {
  fGains=new TArrayF(j);
 }

 Int_t size=fGains->GetSize();

 if (j>size)
 {
  fGains->Set(j);
 }

 fGains->AddAt(float(gain),j-1);

 Int_t oflag=GetOffsetFlag(j);

 SetCalFlags(1,oflag,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetGain(Double_t gain,TString name)
{
// Store gain value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetGain(gain,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetOffset(Double_t off,Int_t j)
{
// Store offset value of the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for offset values, the number of reserved slots for the offset
// values is increased automatically.

 if (j<1) 
 {
  cout << " *AliAttrib::GetOffset* Invalid argument j = " << j << endl;
  return;
 }

 if (!fOffsets)
 {
  fOffsets=new TArrayF(j);
 }

 Int_t size=fOffsets->GetSize();

 if (j>size)
 {
  fOffsets->Set(j);
 }

 fOffsets->AddAt(float(off),j-1);

 Int_t gflag=GetGainFlag(j);

 SetCalFlags(gflag,1,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetOffset(Double_t off,TString name)
{
// Store offset value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetOffset(off,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetCalFlags(Int_t gainflag,Int_t offsetflag,Int_t j)
{
// Store calibration flags of the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the calib. flags, the number of reserved slots for the calib.
// flags is increased automatically.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::GetCalFlags* Invalid argument j = " << j << endl;
  return;
 }

 if (!fCalflags)
 {
  fCalflags=new TArrayI(j);
 }

 Int_t size=fCalflags->GetSize();

 if (j>size)
 {
  fCalflags->Set(j);
 }

 Int_t edge=GetEdgeValue(j);
 Int_t dead=GetDeadValue(j);

 Int_t word=1000*edge+100*dead+10*gainflag+offsetflag;
 
 fCalflags->AddAt(word,j-1);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetGainFlag(Int_t j) const
{
// Provide gain flag of the j-th (default j=1) attribute slot.
//
// flag = 1 : Gain was set
//        0 : Gain was not set
//
// Note : The first attribute slot is at j=1.
// In case j is invalid, 0 is returned.

 if (j<1) 
 {
  cout << " *AliAttrib::GetGainFlag* Invalid argument j = " << j << endl;
  return 0;
 }
 Int_t gflag=0;
 if (fCalflags)
 {
  if (j>0 && j<=(fCalflags->GetSize()))
  {
   Int_t word=fCalflags->At(j-1);
   word=word%100;
   gflag=word/10;
  }
 }
 return gflag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetGainFlag(TString name) const
{
// Provide gain flag of the name-specified attribute slot.
//
// flag = 1 : Gain was set
//        0 : Gain was not set
//
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t flag=0;
 if (j>0) flag=GetGainFlag(j);
 return flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetOffsetFlag(Int_t j) const
{
// Provide offset flag of the j-th (default j=1) attribute slot.
//
// flag = 1 : Offset was set
//        0 : Offset was not set
//
// Note : The first attribute slot is at j=1.
// In case j is invalid, 0 is returned.

 if (j<1) 
 {
  cout << " *AliAttrib::GetOffsetFlag* Invalid argument j = " << j << endl;
  return 0;
 }

 Int_t oflag=0;
 if (fCalflags)
 {
  if (j>0 && j<=(fCalflags->GetSize()))
  {
   Int_t word=fCalflags->At(j-1);
   oflag=word%10;
  }
 }
 return oflag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetOffsetFlag(TString name) const
{
// Provide ofset flag of the name-specified attribute slot.
//
// flag = 1 : Offset was set
//        0 : Offset was not set
//
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t flag=0;
 if (j>0) flag=GetOffsetFlag(j);
 return flag;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliAttrib::GetGain(Int_t j) const
{
// Provide gain value of the j-th (default j=1) attribute slot.
// The first attribute slot is at j=1.
// In case no gain value was set or the argument j is invalid, 0 is returned.
// Note : Use GetGainFlag(j) to check whether this gain was set or not.

 if (j<1) 
 {
  cout << " *AliAttrib::GetGain* Invalid argument j = " << j << endl;
  return 0;
 }

 Float_t gain=0;
 if (fGains)
 {
  if (j>0 && j<=(fGains->GetSize()))
  {
   if (GetGainFlag(j)) gain=fGains->At(j-1);
  }
 }
 return gain;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliAttrib::GetGain(TString name) const
{
// Provide gain value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Float_t gain=0;
 if (j>0) gain=GetGain(j);
 return gain;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliAttrib::GetOffset(Int_t j) const
{
// Provide offset value of the j-th (default j=1) attribute slot.
// The first attribute slot at j=1.
// In case no offset value was set or the argument j is invalid, 0 is returned.
// Note : Use GetOffsetFlag(j) to check whether this offset was set or not.

 if (j<1) 
 {
  cout << " *AliAttrib::GetOffset* Invalid argument j = " << j << endl;
  return 0;
 }

 Float_t offset=0;
 if (fOffsets)
 {
  if (j>0 && j<=(fOffsets->GetSize()))
  {
   if (GetOffsetFlag(j)) offset=fOffsets->At(j-1);
  }
 }
 return offset;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliAttrib::GetOffset(TString name) const
{
// Provide offset value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Float_t offset=0;
 if (j>0) offset=GetOffset(j);
 return offset;
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::ResetGain(Int_t j)
{
// Reset the gain value of the j-th (default j=1) attribute slot.
// Notes : The first attribute slot is at j=1.
//         j=0 ==> All gain values will be reset.
 
 if (!fGains) return;

 Int_t size=fGains->GetSize();

 if ((j>=0) && (j<=size))
 {
  if (j)
  {
   fGains->AddAt(0,j-1);
   Int_t oflag=GetOffsetFlag(j);
   SetCalFlags(0,oflag,j);
  }
  else
  {
   for (Int_t i=0; i<size; i++)
   {
    fGains->AddAt(0,i);
    Int_t oflag=GetOffsetFlag(i);
    SetCalFlags(0,oflag,i);
   }
  }
 }
 else
 {
  cout << " *AliAttrib::ResetGain* Index j = " << j << " invalid." << endl;
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::ResetGain(TString name)
{
// Reset the gain value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) ResetGain(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::ResetOffset(Int_t j)
{
// Reset the offset value of the j-th (default j=1) attribute slot.
// Notes : The first attribute slot is at j=1.
//         j=0 ==> All offset values will be reset.
 
 if (!fOffsets) return;

 Int_t size=fOffsets->GetSize();

 if ((j>=0) && (j<=size))
 {
  if (j)
  {
   fOffsets->AddAt(0,j-1);
   Int_t gflag=GetGainFlag(j);
   SetCalFlags(gflag,0,j);
  }
  else
  {
   for (Int_t i=0; i<size; i++)
   {
    fOffsets->AddAt(0,i);
    Int_t gflag=GetGainFlag(i);
    SetCalFlags(gflag,0,i);
   }
  }
 }
 else
 {
  cout << " *AliAttrib::ResetOffset* Index j = " << j << " invalid." << endl;
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::ResetOffset(TString name)
{
// Reset the offset value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) ResetOffset(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::DeleteCalibrations(Int_t mode)
{
// User selected delete of all gains and/or offsets.
// mode = 0 : All attributes (names, gains, offsets, edge and dead values) are deleted.
//        1 : Only the gains are deleted.
//        2 : Only the offsets are deleted.
//        3 : Both gains and offsets are deleted, but names, edge and dead values are kept.
//
// The default when invoking DeleteCalibrations() corresponds to mode=0.

 if (mode<0 || mode>3)
 {
  cout << " *AliAttrib::DeleteCalibrations* Unknown mode : " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 if (mode==0 || mode==3)
 {
  ResetGain(0);
  if (fGains)
  {
   delete fGains;
   fGains=0;
  }
  ResetOffset(0);
  if (fOffsets)
  {
   delete fOffsets;
   fOffsets=0;
  }
  if (fCalflags && mode==0)
  {
   delete fCalflags;
   fCalflags=0;
  }
  if (fNames && mode==0)
  {
   delete fNames;
   fNames=0;
  }
  return;
 }

 if (mode==1)
 {
  ResetGain(0);
  if (fGains)
  {
   delete fGains;
   fGains=0;
  }
 }
 else
 {
  ResetOffset(0);
  if (fOffsets)
  {
   delete fOffsets;
   fOffsets=0;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetDead(Int_t j)
{
// Set the dead flag to 1 for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, the number of reserved slots for the flags
// is increased automatically.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::SetDead* Invalid argument j = " << j << endl;
  return;
 }

 if (!fCalflags)
 {
  fCalflags=new TArrayI(j);
 }

 Int_t size=fCalflags->GetSize();

 if (j>size)
 {
  fCalflags->Set(j);
 }

 Int_t dead=1;
 Int_t oflag=GetOffsetFlag(j);
 Int_t gflag=GetGainFlag(j);
 Int_t edge=GetEdgeValue(j);

 Int_t word=1000*edge+100*dead+10*gflag+oflag;
 
 fCalflags->AddAt(word,j-1);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetDead(TString name)
{
// Set the dead flag to 1 for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetDead(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetAlive(Int_t j)
{
// Set the dead flag to 0 for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, no action is taken since by default the dead flag is 0.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::SetAlive* Invalid argument j = " << j << endl;
  return;
 }

 if (!fCalflags || j>fCalflags->GetSize()) return;

 Int_t dead=0;
 Int_t oflag=GetOffsetFlag(j);
 Int_t gflag=GetGainFlag(j);
 Int_t edge=GetEdgeValue(j);

 Int_t word=1000*edge+100*dead+10*gflag+oflag;
 
 fCalflags->AddAt(word,j-1);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetAlive(TString name)
{
// Set the dead flag to 0 for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetAlive(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetEdgeOn(Int_t j)
{
// Set the edge value to 1 for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, the number of reserved slots for the flags
// is increased automatically.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::SetEdgeOn* Invalid argument j = " << j << endl;
  return;
 }

 SetEdgeValue(1,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetEdgeOn(TString name)
{
// Set the edge value to 1 for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetEdgeOn(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetEdgeOff(Int_t j)
{
// Set the edge value to 0 for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, no action is taken since by default the edge flag is 0.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::SetEdgeOff* Invalid argument j = " << j << endl;
  return;
 }

 if (!fCalflags || j>fCalflags->GetSize()) return;

 SetEdgeValue(0,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetEdgeOff(TString name)
{
// Set the edge value to 0 for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetEdgeOff(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetEdgeValue(Int_t val,Int_t j)
{
// Set the edge value to "val" for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, the number of reserved slots for the flags
// is increased automatically.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::SetEdgeValue* Invalid argument j = " << j << endl;
  return;
 }

 if (!fCalflags)
 {
  fCalflags=new TArrayI(j);
 }

 Int_t size=fCalflags->GetSize();

 if (j>size)
 {
  fCalflags->Set(j);
 }

 Int_t edge=val;
 Int_t dead=GetDeadValue(j);
 Int_t gflag=GetGainFlag(j);
 Int_t oflag=GetOffsetFlag(j);

 Int_t word=1000*edge+100*dead+10*gflag+oflag;
 
 fCalflags->AddAt(word,j-1);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetEdgeValue(Int_t val,TString name)
{
// Set the edge value to "val" for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) SetEdgeValue(val,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::IncreaseEdgeValue(Int_t j)
{
// Increase the edge value by 1 for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, the number of reserved slots for the flags
// is increased automatically.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::IncreaseEdgeValue* Invalid argument j = " << j << endl;
  return;
 }

 Int_t edge=GetEdgeValue();
 SetEdgeValue(edge+1,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::IncreaseEdgeValue(TString name)
{
// Increase the edge value by 1 for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) IncreaseEdgeValue(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::DecreaseEdgeValue(Int_t j)
{
// Decrease the edge value by 1 for the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, the number of reserved slots for the flags
// is increased automatically.
// The value stored is : 1000*edge + 100*dead + 10*gainflag + offsetflag.

 if (j<1) 
 {
  cout << " *AliAttrib::DecreaseEdgeValue* Invalid argument j = " << j << endl;
  return;
 }

 Int_t edge=GetEdgeValue();
 SetEdgeValue(edge-1,j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::DecreaseEdgeValue(TString name)
{
// Decrease the edge value by 1 for the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) DecreaseEdgeValue(j);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetEdgeValue(Int_t j) const
{
// Provide edge value of the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case j is invalid, 0 is returned.

 if (j<1) 
 {
  cout << " *AliAttrib::GetEdgeValue* Invalid argument j = " << j << endl;
  return 0;
 }

 Int_t edge=0;
 if (fCalflags)
 {
  if (j>0 && j<=(fCalflags->GetSize()))
  {
   Int_t word=fCalflags->At(j-1);
   edge=word/1000;
  }
 }
 return edge;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetEdgeValue(TString name) const
{
// Provide edge value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t val=0;
 if (j>0) val=GetEdgeValue(j);
 return val;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetDeadValue(Int_t j) const
{
// Provide dead value of the j-th (default j=1) attribute slot.
// Note : The first attribute slot is at j=1.
// In case j is invalid, 0 is returned.

 if (j<1) 
 {
  cout << " *AliAttrib::GetDeadValue* Invalid argument j = " << j << endl;
  return 0;
 }

 Int_t dead=0;
 if (fCalflags)
 {
  if (j>0 && j<=(fCalflags->GetSize()))
  {
   Int_t word=fCalflags->At(j-1);
   word=word%1000;
   dead=word/100;
  }
 }
 return dead;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetDeadValue(TString name) const
{
// Provide dead value of the name-specified attribute slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t val=0;
 if (j>0) val=GetDeadValue(j);
 return val;
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::SetSlotName(TString s,Int_t j)
{
// Set a user defined name for the j-th (default j=1) slot. 
// Note : The first attribute slot is at j=1.

 if (j<1) 
 {
  cout << " *AliAttrib::SetSlotName* Invalid argument j = " << j << endl;
  return;
 }

 if (!fNames)
 {
  fNames=new TObjArray(j);
  fNames->SetOwner();
 }

 if (j>fNames->GetSize()) fNames->Expand(j);

 TObjString* so=(TObjString*)fNames->At(j-1);
 if (!so)
 {
  so=new TObjString(s.Data());
  fNames->AddAt(so,j-1);
 }
 else
 {
  so->SetString(s);
 }
}
///////////////////////////////////////////////////////////////////////////
TString AliAttrib::GetSlotName(Int_t j) const
{
// Provide the user defined name for the j-th (default j=1) slot. 
// Note : The first attribute slot is at j=1.

 TString s="";
 if (j<1) 
 {
  cout << " *AliAttrib::GetSlotName* Invalid argument j = " << j << endl;
  return s;
 }

 if (fNames)
 {
  if (j<=fNames->GetSize())
  {
   TObjString* so=(TObjString*)fNames->At(j-1);
   if (so) s=so->GetString();
  }
 }
 return s;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliAttrib::GetSlotIndex(TString name) const
{
// Provide the slot index for the matching name.
// If no matching name is found, 0 is returned.
// Note : The first attribute slot is at j=1.

 Int_t index=0;

 if (fNames)
 {
  TString s;
  Int_t size=fNames->GetSize();
  for (Int_t i=0; i<size; i++)
  {
   TObjString* so=(TObjString*)fNames->At(i);
   if (so) s=so->GetString();
   if (s==name)
   {
    index=i+1;
    break;
   }
  }
 }
 return index;
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::List(Int_t j) const
{
// Provide attribute information for the j-th slot.
// The first slot is at j=1.
// In case j=0 (default) the data of all slots will be listed.

 if (j<0) 
 {
  cout << " *AliAttrib::Data* Invalid argument j = " << j << endl;
  return;
 }

 if (j>0)
 {
  if (GetGainFlag(j)) cout << " gain : " << GetGain(j);
  if (GetOffsetFlag(j)) cout << " offset : " << GetOffset(j);
  if (GetEdgeValue(j)) cout << " edge : " << GetEdgeValue(j);
  if (GetDeadValue(j)) cout << " dead : " << GetDeadValue(j);
  TString s=GetSlotName(j);
  if (s!="") cout << " name : " << s.Data();
 }
 else
 {
  Int_t ng=GetNgains();
  Int_t no=GetNoffsets();
  Int_t nf=0;
  if (fCalflags) nf=fCalflags->GetSize();
  Int_t nn=GetNnames();
  Int_t n=ng;
  if (n<no) n=no;
  if (n<nn) n=nn;
  if (n<nf) n=nf;
  Int_t printf=0;
  TString s;
  for (Int_t i=1; i<=n; i++)
  {
   printf=0;
   if (GetGainFlag(i))   {cout << " gain : " << GetGain(i); printf=1;}
   if (GetOffsetFlag(i)) {cout << " offset : " << GetOffset(i); printf=1;}
   if (GetEdgeValue(i))  {cout << " edge : " << GetEdgeValue(i); printf=1;}
   if (GetDeadValue(i))  {cout << " dead : " << GetDeadValue(i); printf=1;}
   s=GetSlotName(i);
   if (s!="") {cout << " name : " << s.Data(); printf=1;}
   if (printf) cout << endl;
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliAttrib::List(TString name) const
{
// Provide attribute information for the name-specified slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) List(j);
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::Load(AliAttrib& a,Int_t j)
{
// Load attributes of the j-th slot of the input AliAttrib into this AliAttrib object.
// 
// Note : if j=0, then all attributes of all slots are loaded
//
// The default is j=0.

 if (j<0) 
 {
  cout << " *AliAttrib::Load* Invalid argument j = " << j << endl;
  return;
 }

 Int_t n=0;

 if (j==0) // load attributes for all slots
 {
  n=a.GetNgains();
  for (Int_t ig=1; ig<=n; ig++)
  {
   if (a.GetGainFlag(ig))
   {
    SetGain(a.GetGain(ig),ig);
   }
   else
   {
    ResetGain(ig);
   }
  }
  n=a.GetNoffsets();
  for (Int_t io=1; io<=n; io++)
  {
   if (a.GetOffsetFlag(io))
   {
    SetOffset(a.GetOffset(io),io);
   }
   else
   {
    ResetOffset(io);
   }
  }
  n=a.GetNcalflags();
  for (Int_t ic=1; ic<=n; ic++)
  {
   SetEdgeValue(a.GetEdgeValue(ic),ic);
   if (a.GetDeadValue(ic))
   {
    SetDead(ic);
   }
   else
   {
    SetAlive(ic);
   }
  }
  n=a.GetNnames();
  {
   TString s;
   for (Int_t in=1; in<=n; in++)
   {
    s=a.GetSlotName(in);
    if (s!="") SetSlotName(s,in);
   }
  }
 }
 else // load attributes for specified j-th slot only
 {
  n=a.GetNgains();
  if (j<=n)
  {
   if (a.GetGainFlag(j))
   {
    SetGain(a.GetGain(j),j);
   }
   else
   {
    ResetGain(j);
   }
  }
  n=a.GetNoffsets();
  if (j<=n)
  {
   if (a.GetOffsetFlag(j))
   {
    SetOffset(a.GetOffset(j),j);
   }
   else
   {
    ResetOffset(j);
   } 
  }
  n=a.GetNcalflags();
  if (j<=n)
  {
   SetEdgeValue(a.GetEdgeValue(j),j);
   if (a.GetDeadValue(j))
   {
    SetDead(j);
   }
   else
   {
    SetAlive(j);
   }
  }
  n=a.GetNnames();
  {
   TString s;
   if (j<=n)
   {
    s=a.GetSlotName(j);
    if (s!="") SetSlotName(s,j);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliAttrib::Load(AliAttrib& a,TString name)
{
// Load attributes of the name-specified slot of the input AliAttrib into
// this AliAttrib object.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) Load(a,j);
}
///////////////////////////////////////////////////////////////////////////
