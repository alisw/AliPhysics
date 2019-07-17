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

#include "AliDirList.h"

#include <iostream>

#include <TClass.h>
#include <TDirectoryFile.h>
#include <TKey.h>

//______________________________________________________________________________
TObject *AliDirList::PHolder::Get(TDirectory *dir)
{
   // Return stored object if any, otherwise read from directory.
   if (fObj) return fObj;
   if (dir) fObj = dir->Get(fName.Data());
   return fObj;
}

const AliDirList::PHolder &AliDirList::PHolder::operator=(const AliDirList::PHolder &other)
{
   if (&other != this) {
      TNamed::operator=(other);
      fObj = other.fObj;
      fClassName = other.fClassName;
   }
   return *this;
}

//______________________________________________________________________________
AliDirList::AliDirList(const char *name, TCollection *list)
           :TNamed(name, "")
{
   // Create a directory list from a collection
   TIter next(list);
   TObject *obj;
   while ((obj = next())) Add(obj);
}

//______________________________________________________________________________
AliDirList::~AliDirList()
{
   Clear();
}

//______________________________________________________________________________
TObject *AliDirList::At(int i)
{
   // Return element at position i.
   if (!fDir || i > fList.size() - 1) return nullptr;
   PHolder &h = fList[i];
   return h.Get(fDir);
}

//______________________________________________________________________________
void AliDirList::Clear(Option_t *)
{
   if (fOwner) {
      for (auto &h : fList)
         h.ClearObject(fOwner);
   }
   fList.clear();
}

//______________________________________________________________________________
TObject *AliDirList::FindObject(const char *name) const {
   for (auto &h : fList) {
      if (!strcmp(h.GetName(), name))
         return ((PHolder&)h).Get(fDir);
   }
   return nullptr;
}

//______________________________________________________________________________
void AliDirList::ls(Option_t *)  const
{
   for (auto &h : fList) {
      TObject *obj = h.GetObject();
      if (obj)
         std::cout << "  MEM:    " << "\t" << h.GetClassName() << "\t" << h.GetName() << "\t" << h.GetTitle() << std::endl;
      else
         std::cout << "  ONFILE: " << "\t" << h.GetClassName() << "\t" << h.GetName() << "\t" << h.GetTitle() << std::endl;
   }
}

//______________________________________________________________________________
void AliDirList::Print(Option_t *option)  const
{
   std::cout << "Directory list: " << GetName() << std::endl;
   ls(option);
}

//______________________________________________________________________________
AliDirList *AliDirList::CreateFrom(const char *dirname, Bool_t mem)
{
   TDirectory *dir = gDirectory->GetDirectory(dirname);
   if (!dir) return nullptr;
   AliDirList *list = new AliDirList(dirname);
   list->ReadFrom("", mem);
   return list;
}

//______________________________________________________________________________
Int_t AliDirList::ReadFrom(const char *dirname, Bool_t mem)
{
   // Read keys from the current directory and create placeholders.
   // Does not actually read the objects. These are read when calling At function.
   if (dirname && dirname[0] != '\0') SetName(dirname);
   fDir = gDirectory->GetDirectory(GetName());
   if (!fDir) return 0;
   TIter next(fDir->GetListOfKeys());
   TKey *key;
   while ((key = (TKey*)next())) {
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (!cl->InheritsFrom(TDirectory::Class())) {
         // Good key, create placeholder for the key name
         TObject *obj = (mem) ? key->ReadObj() : nullptr;
         fList.push_back(PHolder(key->GetName(), key->GetTitle(), key->GetClassName(), obj));
      }
   }
   return 0;
}

//______________________________________________________________________________
Int_t AliDirList::Write(const char *name, Int_t, Int_t)
{
   // Write content to directory. If directory does not exist, create it in the 
   // current directory with the name of the list.
   Int_t sum = 0;
   if (name && name[0] != '\0')
      SetName(name);
   TDirectory *crtdir = gDirectory;
   fDir = crtdir->GetDirectory(GetName());
   if (!fDir)
      fDir = crtdir->mkdir(GetName());
   if (!fDir->InheritsFrom(TDirectoryFile::Class())) {
      Error("Write", "Current directory not connected to a file. List not written");
      crtdir->cd();
      return 0;
   }
   fDir->cd();
   for (auto &h : fList)
      sum += h.Get()->Write();
   crtdir->cd();
   return sum;
}
