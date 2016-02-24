#include "TFile.h"
#include "AliMUONDataInterface.h"
#include "AliMUONVDigitStore.h"
#include "TTree.h"
#include "AliLog.h"
#include "TMath.h"
#include "Riostream.h"
#include "TKey.h"
#include "AliMUONVClusterStore.h"

/// \ingroup macros
/// \file CompareDigitsOrClusters.C
/// \brief Macro to quickly check if two MUON.Digits.root or two MUON.RecPoints are identical
///
/// Usage :
///
/// .L CompareDigitsOrClusters.C+
/// CompareClusters("dir1/MUON.RecPoints.root","dir2/MUON.RecPoints.root");
/// CompareDigits("dir1/MUON.Digits.root","dir2/MUON.Digits.root");
///
/// If the two files are identical, returns 0, otherwise will create a pair of text file (file1.#event.txt
/// and file2.#event.txt for each event where the information differs,
/// containing a printout of the corresponding stores (either cluster or digits).
///
/// \author: Laurent Aphecetche, Subatech

using namespace std;

template<class T>
T* Store(TFile& file, Long64_t event)
{
	Bool_t ok = file.cd(Form("Event%lld",event));
	if (!ok)
	{
		AliErrorGeneral("Store",Form("Cannot cd to Event%lld",event));
		return 0x0;
	}

	TObject* key = gDirectory->GetListOfKeys()->First();

	TTree* tree = static_cast<TTree*>(gDirectory->Get(key->GetName()));

	if (!tree)
	{
		AliErrorGeneral("Store",Form("Cannot get %s from file for event %lld",key->GetName(),event));
		return 0x0;
	}

	T* store = T::Create(*tree);

	if (!store)
	{
		AliErrorGeneral("Store","Cannot get store from tree");
		return 0x0;
	}

	store->Connect(*tree);

	tree->GetEntry(0);

	return store;
}

template<class T>
Int_t Compare(const T& ds1, const T& ds2)
{
	if ( ds1.GetSize() != ds2.GetSize() )
	{
		return 1;
	}

	TIter next1(ds1.CreateIterator());
	TIter next2(ds2.CreateIterator());

	TObject* d1;

	while ( ( d1 = next1()))
	{
		TObject* d2 = next2();
		if ( d1->Compare(d2) ) return 2;
	}

	return 0;
}

TFile* Open(const char* filename, Long64_t& nevents)
{
	TFile* file = TFile::Open(filename);
	if (!file)
	{
		AliErrorGeneral("Open",Form("Cannot open file %s",filename));
		return 0x0;
	}
	TIter next(file->GetListOfKeys());
	TKey* key;

	nevents = 0;
	while ( ( key = static_cast<TKey*>(next())))
	{
		if ( TString(key->GetName()).BeginsWith("Event") )
		{
			++nevents;
		}
	}

	return file;
}

void DumpSorted(const char* filename, Long64_t event, const AliMUONVStore& store)
{
  /// Dump the given store, in sorted order

  TIter next(store.CreateIterator());
  TObject* object;
  TList list;
  list.SetOwner(kFALSE);

  while ( ( object = next() ) )
  {
    list.Add(object);
  }

  list.Sort();

  std::ofstream out(Form("%s.%lld.txt",filename,event));

  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to filename

  list.Print();

  std::cout.rdbuf(coutbuf); //reset to standard output again
}

template<class T>
Int_t Compare(const char* file1, const char* file2, Long64_t maxEvents=0)
{
	Long64_t n1;
	Long64_t n2;

	TFile* f1 = Open(file1,n1);
	TFile* f2 = Open(file2,n2);

	if ( n1 != n2 )
	{
		AliWarningGeneral("CompareDigits",Form("file1 has %lld events while file2 has %lld ...",
				n1,n2));
	}

	if ( !maxEvents )
	{
		maxEvents = TMath::Min(n1,n2);
	}

	Long64_t n = TMath::Min(maxEvents,n1);

	n = TMath::Min(n,n2);

	for ( Long64_t i = 0; i < n; ++i )
	{
		T* store1 = Store<T>(*f1,i);
		if (!store1)
		{
			AliErrorGeneral("Compare",Form("Cannot get store from file1 for event %lld",i));
			continue;
		}
		T* store2 = Store<T>(*f2,i);
		if (!store2)
		{
			AliErrorGeneral("Compare",Form("Cannot get store from file2 for event %lld",i));
			continue;
		}

		if ( Compare<T>(*store1,*store2) )
		{
			AliErrorGeneral("Compare",Form("differ for event %lld",i));
			DumpSorted("file1",i,*store1);
			DumpSorted("file2",i,*store2);
		}

		delete store1;
		delete store2;
	}

	return 0;
}

Int_t CompareDigits(const char* file1, const char* file2, Long64_t maxEvents=0)
{
	return Compare<AliMUONVDigitStore>(file1,file2,maxEvents);
}

Int_t CompareClusters(const char* file1, const char* file2, Long64_t maxEvents=0)
{
	return Compare<AliMUONVClusterStore>(file1,file2,maxEvents);
}
