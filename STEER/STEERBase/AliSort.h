#ifndef ALISORT_H
#define ALISORT_H
#include <TClonesArray.h>

// fast implementations of sorting for TClonesArray by making
// use of element type information directly
// created 3.5.2016; sandro.wenzel@cern.ch

namespace AliSort {

// anonymous namespace to prevent using this helper class from outside
namespace {
  // An *internal* class to gain access to fLast and fKeep of a TClonesArray via inheritance
  class AliClonesArrayWrapper : public TClonesArray {
    public:
      TObjArray *GetKeep() const {return fKeep;}
      void SetLast(Int_t l) {fLast = l;}
      void SetSorted(bool b) {fSorted = b;}
      Int_t GetAbsLastM() const {return this->GetAbsLast();}
      Int_t GetLowerBound() const {return fLowerBound;}
  };
}

  // a helper function for the QSortT algorithm
  template <typename T>
  inline Int_t ObjCompareT(TObject *a, TObject *b)
  {
     if (a == nullptr && b == nullptr) return 0;
     if (a == nullptr) return 1;
     if (b == nullptr) return -1;
     return ((T*)a)->T::Compare(b);
  }

  template <typename T>
  inline void QSortT(TObject **a, Int_t nBs, TObject ***b, Int_t first, Int_t last)
  {
     // in contrast to ROOT originals implementation we do not protect for thread safety here
     // since AliROOT does not use threads

     static TObject *tmp1, **tmp2;
     static int i; // "static" to save stack space
     int j,k;

     static int depth = 0;
     if (depth == 0 && nBs > 0) tmp2 = new TObject*[nBs];
     depth++;

     while (last - first > 1) {
        i = first;
        j = last;
        for (;;) {
           while (++i < last && ObjCompareT<T>(a[i], a[first]) < 0) {}
           while (--j > first && ObjCompareT<T>(a[j], a[first]) > 0) {}
           if (i >= j) break;

           tmp1 = a[i]; for(k=0;k<nBs;k++) tmp2[k] = b[k][i];
           a[i] = a[j]; for(k=0;k<nBs;k++) b[k][i] = b[k][j];
           a[j] = tmp1; for(k=0;k<nBs;k++) b[k][j] = tmp2[k];
        }
        if (j == first) {
           ++first;
           continue;
        }
        tmp1 = a[first]; for(k=0;k<nBs;k++) tmp2[k] = b[k][first];
        a[first] = a[j]; for(k=0;k<nBs;k++) b[k][first] = b[k][j];
        a[j] = tmp1; for(k=0;k<nBs;k++) b[k][j] = tmp2[k];
        if (j - first < last - (j + 1)) {
           QSortT<T>(a, nBs, b, first, j);
           first = j + 1; // QSort(j + 1, last);
        } else {
           QSortT<T>(a, nBs, b, j + 1, last);
           last = j; // QSort(first, j);
        }
     }
     depth--;

     if (depth == 0 && nBs > 0) delete [] tmp2;
  }

  template <typename T>
  inline void QSortT(TObject **a, TObject **b, Int_t first, Int_t last) { QSortT<T>(a, 1, &b, first, last); }

  // a fast template replacement for ROOT's TClonesArray Sort
  // sorts a TClonesArray fast by using information about the type that
  // the TClonesArray is storing
  // this implementation is taken from ROOT
  template <typename T>
  inline void TClonesArraySort(TClonesArray *a, Int_t upto=kMaxInt){
      AliClonesArrayWrapper *wrapper = static_cast<AliClonesArrayWrapper *>(a);
      Int_t nentries = wrapper->GetAbsLastM()+1;
      if (nentries <= 0 || wrapper->TClonesArray::IsSorted()) return;

      TObject **rawcontainer = wrapper->GetObjectRef();

      // this test is a quite expensive runtime test
      // which is essentially useless for TClonesArray since by definition
      // all elements have the same type
#ifndef NDEBUG
      for (Int_t i = 0; i < wrapper->TClonesArray::GetSize(); ++i)
         if (rawcontainer[i]) {
            if (!((T*)rawcontainer[i])->T::IsSortable()) {
              // Error("Sort", "objects in array are not sortable");
               return;
            }
         }
#endif

      Int_t lowerbound = wrapper->GetLowerBound();

      TObject **othercontainer = wrapper->GetKeep()->GetObjectRef();
      QSortT<T>(rawcontainer, othercontainer, 0, std::min(nentries, upto-lowerbound));

      wrapper->SetLast(-2);
      wrapper->SetSorted(kTRUE);
  }

}


#endif // ALISORT_H
