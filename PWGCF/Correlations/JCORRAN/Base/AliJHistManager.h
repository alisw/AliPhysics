// @(#)utils:$Id$
// Author: Beonsu Chang   26/12/94
#ifndef UTILS_ALIJARRAY_H
#define UTILS_ALIJARRAY_H

#include <vector>
#include <map>
#include <TString.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TPRegexp.h>
#include <TVector.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>
#include <iostream>
#include <TClass.h>

#define JERROR(x)  {std::cout<<"!!! JERROR : "<<x<<" "<<__LINE__<<" "<<__FILE__<<" "<<std::endl , gSystem->Exit(100); }
#define JDEBUG(x,y)  if(x<100){std::cout<<"JDEBUG : "<<#x<<" : "<<(y)<<" "<<__LINE__<<" "<<__FILE__<<" "<<std::endl;}

//class AliJMap;
class AliJArrayBase;
class AliJBin;
class AliJArrayAlgorithm;
class AliJArrayAlgorithmSimple;
class AliJTH1;
class AliJHistManager;
template<typename t> class AliJTH1Derived;
template<typename t> class AliJTH1DerivedPlayer;

//////////////////////////////////////////////////////
//  Utils
//////////////////////////////////////////////////////
std::vector<TString> Tokenize( TString s, TString d, int quote=1 );
TString Join( std::vector<TString>& ss , TString del=" ");
bool OutOf( int, int, int );

typedef std::vector<int> ArrayInt;
typedef std::vector<double> ArrayDouble;
//typedef AliJArray<void*> ArrayVoid*;

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJNamed                                                              //
//                                                                      //
// 
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class AliJNamed{
    public:
        AliJNamed( TString name, TString title, TString opt, int mode );
        virtual ~AliJNamed();
        TString GetName(){ return fName; }
        TString GetTitle(){ return fTitle; }
        TString GetOption(){ return fOption; }
        TString GetOption(TString key);
        bool    HasOption(TString key){ return GetOption(key)!=UndefinedOption(); }
        bool    HasOption(TString key, TString val){ return GetOption(key)==val; } // TODO sensitive?
        void SetName( const char * s ) { fName=s; }
        void SetTitle( const char * s ) { fTitle=s; }
        void SetNameTitle( TString n, TString t ){ SetName(n);SetTitle(t); }
        void SetFullOption( const char * s ) { fOption=s; }
        void SetOption( TString key, TString value="" );
        void RemoveOption( TString key ); // TODO
        //void SetOptionWithString( TString s );// TODO
        static TString UndefinedOption();
    protected:
        TString  fName;
        TString  fTitle;
        TString  fOption;
        int      fMode;
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJBin                                                              //
//                                                                      //
// 
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class AliJBin : public AliJNamed {
    public:
        enum { kSingle, kRange, kString , kNMode};
        AliJBin();
        AliJBin(TString config, AliJHistManager * hmg);
        AliJBin(const AliJBin& obj);
        AliJBin& operator=(const AliJBin& obj);
        AliJBin & Set( TString name, TString iname, TString Title, int mode=kRange);
        void      AddToManager( AliJHistManager* hmg );
        AliJBin & SetBin( const int n, const float * v );
        AliJBin & SetBin( const int n, const double * v );
        AliJBin & SetBin( TVector * v );
        AliJBin & SetBin( const TString  v );
        AliJBin & SetBin( const int n );

        int GetBin( double x );

        double GetMin(){ return fBinD[0]; }
        double GetMax(){ return fBinD[RawSize()-1]; }

        TString BuildTitle( int i );
        int RawSize(){ return fBinD.size(); }
        int Size(){ return fMode==kRange?fBinD.size()-1:fBinD.size(); }
        double At(int i){ return fBinD[i]; }
        TString GetIndexName(){ return fIndexName; }

        void Print();
        TString GetString();

        operator int(){ return Size(); }

        static TString GetModeString(int i);
        static int GetMode( TString mode );
    private:
        void AddBin( const TString & v );
        void AddBin( float v );
        virtual void FixBin();

        std::vector<double>  fBinD;
        std::vector<TString> fBinStr;
        bool    fIsFixedBin;
        TString fIndexName;
        AliJHistManager * fHMG;
};



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJArrayBase                                                        //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//________________________________________________________________________
class AliJArrayBase : public AliJNamed{
    public:
        enum { kNormal, kSingle };
        virtual ~AliJArrayBase();
        AliJArrayBase& operator=(const AliJArrayBase& obj);

        int AddDim( int i){ fDim.push_back(i);return fDim.size();}
        int Dimension(){ return fDim.size(); }
        int GetEntries(){ return fArraySize; }
        int SizeOf(int i) { return fDim.at(i); }

        ArrayInt& Index(){ return fIndex; }
        int  Index( int d );
        void SetIndex( int i, int d );
        void ClearIndex(){ fIndex.clear();fIndex.resize( Dimension(), 0 ); }

        void * GetItem();
        void * GetSingleItem();

        ///void LockBin(bool is=true){}//TODO
        //bool IsBinLocked(){ return fIsBinLocked; }

        virtual void FixBin();
        bool IsBinFixed(){ return fIsBinFixed; }

        bool OutOfDim( int d ){ return OutOf( d, 0, Dimension()-1 ); }
        bool OutOfSize( int i, int d ){ return OutOfDim(d) || OutOf( i, 0, SizeOf(d)-1); }


        // Virtual 
        virtual void *  BuildItem()=0;
        virtual TString BuildName()=0;
        virtual TString BuildTitle()=0;
        virtual void  Print()=0;
        virtual TString GetString()=0;

        //int Resize( int size, int dim=-1 ); // NextStep 
        void    InitIterator();
        bool    Next(void *& item ); 

    protected:
        AliJArrayBase(); // Prevent direct creation of AliJArrayBase
        AliJArrayBase(const AliJArrayBase& obj);

        ArrayInt        fDim;           // Comment test
        ArrayInt        fIndex;         /// Comment test
        int         fArraySize;         /// Comment test3
        int         fNGenerated;
        bool        fIsBinFixed;
        bool        fIsBinLocked;
        AliJArrayAlgorithm * fAlg;
        friend class AliJArrayAlgorithm;
};


//________________________________________________________________________
class AliJArrayAlgorithm {
    public:
        AliJArrayAlgorithm(AliJArrayBase * cmd); //TODO Move to private
        AliJArrayAlgorithm(const AliJArrayAlgorithm& obj);
        AliJArrayAlgorithm& operator=(const AliJArrayAlgorithm& obj);
        virtual ~AliJArrayAlgorithm();
        int Dimension(){ return fCMD->Dimension(); }
        int SizeOf(int i){ return fCMD->SizeOf(i); }
        int GetEntries(){ return fCMD->GetEntries(); }
        int Index(int i){ return fCMD->Index(i); }
        virtual int BuildArray()=0;
        virtual void * GetItem()=0;
        virtual void SetItem(void * item)=0;
        virtual void InitIterator()=0;
        virtual bool Next(void *& item) = 0;
        virtual void ** GetRawItem()=0;
        virtual void * GetPosition()=0;
        virtual bool IsCurrentPosition(void * pos)=0;
        virtual void SetPosition(void * pos )=0;
        virtual void DeletePosition( void * pos ) =0;
    protected:
        AliJArrayBase * fCMD;
};

//________________________________________________________________________
class AliJArrayAlgorithmSimple : public AliJArrayAlgorithm {
    public:
        AliJArrayAlgorithmSimple( AliJArrayBase * cmd);
        AliJArrayAlgorithmSimple(const AliJArrayAlgorithmSimple& obj);
        AliJArrayAlgorithmSimple& operator=(const AliJArrayAlgorithmSimple& obj);
        virtual ~AliJArrayAlgorithmSimple();
        virtual int BuildArray();
        int  GlobalIndex();
        void ReverseIndex(int iG );
        virtual void * GetItem();
        virtual void SetItem(void * item);
        virtual void InitIterator(){ fPos = 0; }
        virtual void ** GetRawItem(){ return &fArray[GlobalIndex()]; }
        virtual bool Next(void *& item){
            item = fPos<GetEntries()?(void*)fArray[fPos]:NULL;
            if( fPos<GetEntries() ) ReverseIndex(fPos);
            return fPos++<GetEntries();
        }
        virtual void * GetPosition(){ return  static_cast<void*> (new int( fPos )); }
        virtual bool IsCurrentPosition(void * pos){ return *static_cast<int*>(pos)==fPos; }
        virtual void  SetPosition(void *pos){ fPos=*static_cast<int*>(pos);ReverseIndex(fPos); } 
        virtual void  DeletePosition(void *pos){ delete static_cast<int*>(pos); }
    private:
        ArrayInt    fDimFactor;
        void    **fArray;
        int     fPos;
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJTH1                                                              //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class AliJTH1 : public AliJArrayBase{
    public:
        AliJTH1();
        AliJTH1(TString config, AliJHistManager * hmg);
        AliJTH1(const AliJTH1& obj);
        AliJTH1& operator=(const AliJTH1& obj);
        virtual ~AliJTH1();

        int AddDim( int i){ fDim.push_back(i);return fDim.size();}
        int AddDim( AliJBin * bin );
        int AddDim( TString v);
        void AddToManager( AliJHistManager * hmg );
        AliJBin* GetBinPtr(int i){ return fBins.at(i); }

        // Virtual from AliJArrayBase
        virtual void * BuildItem() ;
        virtual TString GetString();
        virtual void  Print();
        virtual void  FixBin();
        // Virtual from this
        virtual Int_t Write();
        //virtual Int_t WriteAll();
        virtual const char * ClassName(){ return "AliJTH1"; }

        // Not Virtual
        virtual TString BuildName();
        virtual TString BuildTitle();
        bool    IsLoadMode(); 
        void    SetTemplate(TH1* h);
        TH1*    GetTemplatePtr(){ return fTemplate; }



    protected:

        TDirectory      *fDirectory;
        TDirectory      *fSubDirectory;
        AliJHistManager *fHMG;
        TH1             *fTemplate;
        std::vector<AliJBin*>  fBins;
};
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJTH1Derived                                                       //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
template< typename T >
class AliJTH1Derived : public AliJTH1 {
    protected:
    public:
        AliJTH1Derived();
        AliJTH1Derived(TString config, AliJHistManager *hmg):
            AliJTH1(config, hmg),fPlayer(this){}
        virtual ~AliJTH1Derived();

        AliJTH1DerivedPlayer<T> & operator[](int i){ fPlayer.Init();fPlayer[i];return fPlayer; }
        T * operator->(){ return static_cast<T*>(GetSingleItem()); }
        operator T*(){ return static_cast<T*>(GetSingleItem()); }
        // Virtual from AliJArrayBase

        // Virtual from AliJTH1
        virtual const char * ClassName(){ return Form("AliJ%s",T::Class()->GetName()); }

        AliJTH1Derived<T>& operator<<(int i){ AddDim(i);return *this; }
        AliJTH1Derived<T>& operator<<(AliJBin& v){ AddDim(&v);return *this; }
        AliJTH1Derived<T>& operator<<(TString v){ AddDim(v);return *this; }
        AliJTH1Derived<T>& operator<<(T v){ SetTemplate(&v);return *this; }
        void SetWith( AliJTH1Derived<T>& v, TString name, TString title="" ){
          SetTemplate( v.GetTemplatePtr() );
          fName = name;
          fTitle = title;
          GetTemplatePtr()->SetName( name );
          GetTemplatePtr()->SetTitle( title );
          for( int i=0;i<v.Dimension();i++ ) AddDim( v.GetBinPtr(i) );
          AddDim("END");
        }
        void SetWith( AliJTH1Derived<T>& v, T tem ){
          SetTemplate( &tem );
          for( int i=0;i<v.Dimension();i++ ) AddDim( v.GetBinPtr(i) );
          AddDim("END");
        }
    protected:
        AliJTH1DerivedPlayer<T> fPlayer;

};


//////////////////////////////////////////////////////////////////////////
// AliJTH1DerivedPlayer                                                 //
//////////////////////////////////////////////////////////////////////////
template< typename T>
class AliJTH1DerivedPlayer {
    public:
        AliJTH1DerivedPlayer( AliJTH1Derived<T> * cmd ):fLevel(0),fCMD(cmd){};
        AliJTH1DerivedPlayer<T>& operator[](int i){
            if( fLevel > fCMD->Dimension() ) { JERROR("Exceed Dimension"); }
            if( OutOf( i, 0,  fCMD->SizeOf(fLevel)-1) ){ JERROR(Form("wrong Index %d of %dth in ",i, fLevel)+fCMD->GetName()); }
            fCMD->SetIndex(i, fLevel++);
            return *this;
        }
        void Init(){ fLevel=0;fCMD->ClearIndex(); }
        T* operator->(){ return static_cast<T*>(fCMD->GetItem()); } 
        operator T*(){ return static_cast<T*>(fCMD->GetItem()); } 
        operator TObject*(){ return static_cast<TObject*>(fCMD->GetItem()); } 
        operator TH1*(){ return static_cast<TH1*>(fCMD->GetItem()); } 
    private:
        int fLevel;
        AliJTH1Derived<T> * fCMD;
};

typedef AliJTH1Derived<TH1D> AliJTH1D;
typedef AliJTH1Derived<TH2D> AliJTH2D;
typedef AliJTH1Derived<TH3D> AliJTH3D;
typedef AliJTH1Derived<TProfile> AliJTProfile;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJHistManager                                                       //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class AliJHistManager: public AliJNamed{
    public:
        AliJHistManager(TString name, TString dirname="" );
        AliJHistManager(const AliJHistManager& obj);
        AliJHistManager& operator=(const AliJHistManager& obj);
        void Add( AliJBin * o );
        void Add( AliJTH1 * o );

        int GetNBin(){ return fBin.size()>fBinNames.size()?fBin.size():fBinNames.size(); } // TODO
        int GetNHist(){ return fHist.size()>fHistNames.size()?fHist.size():fHistNames.size(); } // TODO
        void Print();
        int  LoadConfig();
        TDirectory * GetDirectory(){ return fDirectory;} 
        void SetDirectory(TDirectory* d){ fDirectory = d; }
        static AliJHistManager* GlobalManager();
        static AliJHistManager* CurrentManager( AliJHistManager * hmg=NULL);
        AliJHistManager * cd(){ return AliJHistManager::CurrentManager(this); }
        void SetLoadMode(bool b=true){ fIsLoadMode = b; }
        bool IsLoadMode(){ return fIsLoadMode; }
        TString GetString(){
            TString st;
            for( int i=0;i<GetNBin();i++ ) st+=fBin[i]->GetString()+"\n";
            for( int i=0;i<GetNHist();i++ ){
                st+=fHist[i]->GetString()+"\n";
            }
            return st;
        }
        void Write();
        void WriteConfig();

        AliJBin * GetBin( TString name); 
        AliJBin * GetBuiltBin( TString name); 
        AliJTH1 * GetTH1( TString name );
        AliJTH1 * GetBuiltTH1( TString name );
        AliJTProfile& GetTProfile( TString name){ return dynamic_cast<AliJTProfile&>(*GetTH1(name)); }
        AliJTH1D& GetTH1D( TString name){ return dynamic_cast<AliJTH1D&>(*GetTH1(name)); }
        AliJTH2D& GetTH2D( TString name){ return dynamic_cast<AliJTH2D&>(*GetTH1(name)); }
        AliJTH3D& GetTH3D( TString name){ return dynamic_cast<AliJTH3D&>(*GetTH1(name)); }
        bool fIsLoadMode;

        TString GetHistName(int i){ return fHistNames[i]; }

        AliJTH1 * GetAliJTH1(int i){ return GetTH1(fHistNames[i]); }
        int GetNAliJTH1(){ return fHistNames.size(); }
        bool HistogramExists(TString name);


    private:
        TDirectory *fDirectory;
        TString                     fConfigStr;
        std::vector<AliJBin*>       fBin;
        std::vector<AliJTH1*>  fHist;
        std::vector<AliJHistManager*> fManager;
        std::vector<TString>        fBinNames;
        std::vector<TString>        fBinConfigs;
        std::vector<TString>        fHistNames;
        std::vector<TString>        fHistConfigs;
  
};

#endif
