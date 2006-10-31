// Author: Cvetan Cheshkov  10/10/2005

/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawEventHeaderBase                                                //
// This a new versioning scheme for raw data root-ification and reading //
// For details look at offline weekly meeting 20/10/2005                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <unistd.h>

#include <TClass.h>
#include <TDataMember.h>
#include <TMethodCall.h>

#include "AliLog.h"
#include "AliRawEventHeaderBase.h"


ClassImp(AliRawEventHeaderBase)

//______________________________________________________________________________
AliRawEventHeaderBase::AliRawEventHeaderBase():
fSize(0),
fMagic(0),
fHeadSize(0),
fVersion(0),
fExtendedDataSize(0),
fExtendedData(NULL),
fIsSwapped(kFALSE)
{
  // Default constructor
}

//______________________________________________________________________________
void *AliRawEventHeaderBase::HeaderBegin()
{
  // Returns the pointer to the first data member
  // beyond the base class data members

  TList *datalist = IsA()->GetListOfDataMembers();
  TIter next(datalist);                           
  TDataMember *member = (TDataMember *)next();

  if(!strcmp(member->GetTypeName(),"TClass"))
    member = (TDataMember *)next();

  return (void *)((char *)this+member->GetOffset());
}

//______________________________________________________________________________
Int_t AliRawEventHeaderBase::HeaderSize() const
{
  // Returns the size of the data members list
  // beyond the base class data members

  Int_t size = 0;

  TList *datalist = IsA()->GetListOfDataMembers();
  TIter next(datalist);                           
  TDataMember *member;
  while ((member=(TDataMember *)next()) != 0x0) {
    if (!strcmp(member->GetTypeName(),"TClass")) continue;
    UInt_t unitsize = member->GetUnitSize();
    UInt_t ndim = member->GetArrayDim();
    if (ndim == 0)
      size += unitsize;
    else
      for(UInt_t i=0;i<ndim;i++) size += member->GetMaxIndex(i)*unitsize;
  }

  return size;
}

//______________________________________________________________________________
void AliRawEventHeaderBase::Swap()
{
   // Swap base header data.
   // Update the fIsSwapped flag which
   // is then use to copy in an appropriate way
   // the rest of the header data from the raw data stream

   if (IsSwapped()) {
      fIsSwapped    = kTRUE;
      fSize         = net2host(fSize);
      fMagic        = net2host(fMagic);
      fHeadSize     = net2host(fHeadSize);
      fVersion      = net2host(fVersion);
   }
}

//______________________________________________________________________________
const char *AliRawEventHeaderBase::GetTypeName()
{
   // Get event type as a string.
   // Will fail in case data header
   // does not contain eventType field
   Int_t eventType = Get("Type");

   switch (eventType) {
      case kStartOfRun:
         return "START_OF_RUN";
         break;
      case kEndOfRun:
         return "END_OF_RUN";
         break;
      case kStartOfRunFiles:
         return "START_OF_RUN_FILES";
         break;
      case kEndOfRunFiles:
         return "END_OF_RUN_FILES";
         break;
      case kStartOfBurst:
         return "START_OF_BURST";
         break;
      case kEndOfBurst:
         return "END_OF_BURST";
         break;
      case kPhysicsEvent:
         return "PHYSICS_EVENT";
         break;
      case kCalibrationEvent:
         return "CALIBRATION_EVENT";
         break;
      case kFormatError:
         return "EVENT_FORMAT_ERROR";
         break;
      case kStartOfData:
	 return "START_OF_DATA";
	 break;
      case kEndOfData:
	 return "END_OF_DATA";
	 break;
      case kSystemSoftwareTriggerEvent:
	 return "SYSTEM_SOFTWARE_TRIGGER_EVENT";
	 break;
      case kDetectorSoftwareTriggerEvent:
	 return "DETECTOR_SOFTWARE_TRIGGER_EVENT";
	 break;
      default:
	 return "UNKNOWN EVENT TYPE NUMBER";
         break;
   }
}

//______________________________________________________________________________
AliRawEventHeaderBase* AliRawEventHeaderBase::Create(char*& data)
{
  // Static method to create AliRawEventHeaderVX object
  // The actual header class version is taken from the
  // raw data

  // First create AlirawEVentHeaderBase class object
  AliRawEventHeaderBase header;

  // Copy the first common part of the raw data header
  memcpy(header.HeaderBaseBegin(), data, header.HeaderBaseSize());
 
    // Swap header data if needed
  if (header.IsSwapped())
    header.Swap();

  // Is header valid...
  if (!header.IsValid()) {
    AliFatalClass("Invalid header format!");
    // try recovery... how?
    return 0x0;
  }

  if (header.GetEventSize() < (UInt_t)header.HeaderBaseSize()) {
    AliFatalClass("Invalid header base size!");
    // try recovery... how?
    return 0x0;
  }

  // Now check the DATE version and create the corresponding header
  // class object
  UInt_t version = header.GetVersion();
  UInt_t majorversion = (version>>16)&0x0000ffff;
  UInt_t minorversion = version&0x0000ffff;
  TString classname;
  classname.Form("AliRawEventHeaderV%d_%d",majorversion,minorversion);
    
  TClass *tcl = TClass::GetClass(classname.Data());
  if (!tcl) {
    AliFatalClass(Form("Unknown header version (%s)!",classname.Data()));
    return 0x0;
  }

  //  header.Dump(); tcl->Dump();

  AliRawEventHeaderBase *hdr = (AliRawEventHeaderBase *)tcl->New();
  if (!hdr) {
    AliFatalClass(Form("Can not create object of class %s",classname.Data()));
    return 0x0;
  }

  // Copy the base header data members and initialize other data members
  memcpy(hdr->HeaderBaseBegin(),header.HeaderBaseBegin(), header.HeaderBaseSize());
  memset(hdr->HeaderBegin(),0, hdr->HeaderSize());
  hdr->fIsSwapped = header.fIsSwapped;

  // Consistency check
  if (hdr->GetEventSize() < ((UInt_t)hdr->HeaderBaseSize() + (UInt_t)hdr->HeaderSize())) {
    AliFatalClass(Form("Invalid header size (%d < %d +%d)!",
		       hdr->GetEventSize(),hdr->HeaderBaseSize(),hdr->HeaderSize()));
    // try recovery... how?
    return 0x0;
  }

  // Check for the presence of header extension and its size
  Int_t extsize = (Int_t)hdr->GetHeadSize() - (hdr->HeaderBaseSize() + hdr->HeaderSize());
  if (extsize < 0) {
    AliFatalClass(Form("Invalid header size (%d < %d +%d)!",
		       hdr->GetHeadSize(),hdr->HeaderBaseSize(),hdr->HeaderSize()));
    // try recovery... how?
    return 0x0;
  }
  else {
    if (extsize > 0) {
      hdr->SetExtendedDataSize(extsize);
      char *extdata = new char[extsize];
      memset(extdata,0,extsize);
      hdr->SetExtendedData(extdata);
    }
  }

  return hdr;
}

//______________________________________________________________________________
Int_t AliRawEventHeaderBase::ReadHeader(char*& data)
{
  // Read header info from DATE data stream.
  // Returns bytes read

  Long_t start = (Long_t)data;
  // Swap header data if needed
  if (DataIsSwapped()) {
    swab(data,HeaderBaseBegin(), HeaderBaseSize());
    data += HeaderBaseSize();
    swab(data, HeaderBegin(), HeaderSize());
    data += HeaderSize();
    if(GetExtendedDataSize()>0) {
      swab(data, GetExtendedData(), GetExtendedDataSize());
      data += GetExtendedDataSize();
    }
  }
  else {
    memcpy(HeaderBaseBegin(), data, HeaderBaseSize());
    data += HeaderBaseSize();
    memcpy(HeaderBegin(), data, HeaderSize());
    data += HeaderSize();
    if(GetExtendedDataSize()>0) {
      memcpy(GetExtendedData(), data, GetExtendedDataSize());
      data += GetExtendedDataSize();
    }
  }

  return (Int_t)((Long_t)data - start);
}

//______________________________________________________________________________
UInt_t AliRawEventHeaderBase::Get(const char *datamember) const
{
  // The method to get a data member from the header object
  // Except for the data members of the base class, all the
  // other header data should be retrieved ONLY by this method
  // The name of the data member should be supplied without "f"
  // in front

  char buf[256] = "f";
  strcat(buf,datamember);

  TDataMember *member = IsA()->GetDataMember(buf);
  if (!member) {
    AliFatal(Form("No data member %s is found! Check the raw data version!",buf));
    return 0;
  }

  if (member->GetArrayDim() != 0) {
    AliFatal(Form("Member %s is an array! Use the GetP() method!",buf));
    return 0;
  }

  if (strcmp(member->GetTypeName(),"UInt_t") != 0) {
    AliFatal(Form("Member %s is not of type UInt_t!",buf));
    return 0;
  }

  const void *pointer = (char *)this+member->GetOffset();

  return *((UInt_t *)pointer);
}

//______________________________________________________________________________
const UInt_t* AliRawEventHeaderBase::GetP(const char *datamember) const
{
  // The method to get a data member from the header object
  // Except for the data members of the base class, all the
  // other header data should be retrieved ONLY by this method
  // The name of the data member should be supplied without "f"
  // in front

  char buf[256] = "f";
  strcat(buf,datamember);

  TDataMember *member = IsA()->GetDataMember(buf);
  if (!member) {
    AliFatal(Form("No data member %s is found! Check the raw data version!",buf));
    return 0;
  }

  //  if (member->GetArrayDim() == 0) {
  //    AliFatal(Form("Member %s is not an array! Use the Get() method!",buf));
  //    return 0;
  //  }

  if (strcmp(member->GetTypeName(),"UInt_t") != 0) {
    AliFatal(Form("Member %s is not of type UInt_t*!",buf));
    return 0;
  }

  const void *pointer = (char *)this+member->GetOffset();

  return (const UInt_t*)pointer;
}
