#ifndef ALIDIRLIST_H
#define ALIDIRLIST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliDirList
/// \author Mihaela Gheata
/// \brief AliDirList
/// A special list handling objects to be written or retrieved from a ROOT directory.
/// The list itself is not persistent, writing it to a file will create first a directory
/// and then call the method Write for all objects in the list. Getting back AliDirList
/// from a file requires having gDirectory pointing to the parent directory in the file.
/// The objects in the list are internally hold into a smart placeholder object. AliDirList
/// will scan the content of the folder to be loaded and create empty placeholders for each
/// object. When the Get() method is invoked, the placeholder will do the actual reading from
/// file and hand over the read object. This is to read in memory only the objects requested
// by the user. A AliDirList cannot contain other AliDirList, all directories from the read
/// directory will be skipped.
/// \date 21/01/2018

#include <vector>

#include <TNamed.h>

class TDirectory;
class TDirectoryFile;
class TCollection;

class AliDirList : public TNamed
{
   // Internal placeholder for AliDirList items
   class PHolder : public TNamed {
      TObject *fObj = nullptr;    //! pointer to object
      TString  fClassName;        //! class name

   public:
      PHolder(const char *name, const char *title = "", const char *classname = "", TObject *obj = nullptr)
         : TNamed(name, title), fObj(obj), fClassName(classname) {}
      PHolder(TObject *obj) : TNamed(obj->GetName(), obj->ClassName()), fObj(obj) {}
      PHolder(const PHolder &other) : TNamed(other), fObj(other.fObj), fClassName(other.fClassName) {}
      ~PHolder() {}

      const PHolder &operator=(const PHolder &other);

      inline void     ClearObject(bool del) { if (del) delete fObj; fObj = nullptr; }
      TObject        *Get(TDirectory *dir = nullptr);
      inline TObject *GetObject() const { return fObj; }
      inline void     Set(TObject *obj) { fObj = obj; }
      const char     *GetClassName() const { return (fObj) ? fObj->ClassName() : fClassName.Data(); }
      void            SetClassName(const char *name) { fClassName = name; }
   };

private:
   Bool_t fOwner = kFALSE;                //! Container owns the content
   TDirectory *fDir = nullptr;           //! Directory containing the the objects to read
   std::vector<PHolder> fList;           //! List content

public:
   AliDirList(const char *name = "") : TNamed(name, "") {}
   AliDirList(const char *name, TCollection *list);
   virtual ~AliDirList();

   // Compatibility methods (TCollection, TList, TObjArray)
   void           Add(TObject *obj) { fList.push_back(PHolder(obj)); }
   TObject       *At(Int_t i);
   TObject       *operator[](Int_t i) { return At(i); }
   virtual void   Clear(Option_t *option = "");
   TObject       *First() { return At(0); }
   TObject       *Last() { return At(Size() - 1); }
   int            GetSize() const { return (Int_t)fList.size(); }
   int            GetEntries() const { return (Int_t)fList.size(); }
   int            GetEntriesFast() const { return (Int_t)fList.size(); }
   TObject       *FindObject(const char *name) const;
   static AliDirList *CreateFrom(const char *dirname, bool mem = kFALSE);

   virtual void   ls (Option_t *option="") const;
   virtual void   Print(Option_t *option="") const;
   void           SetOwner(Bool_t flag = true) { fOwner = flag; }
   size_t         Size() const { return fList.size(); }
   Int_t          Write (const char *dir="", Int_t option=0, Int_t bufsize=0);
   Int_t          ReadFrom(const char *dir="", bool mem = kFALSE);

   ClassDef(AliDirList, 0)  // Class representing a special list of objects from a ROOT directory
};
#endif // ALIDIRLIST_H
