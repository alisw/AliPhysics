#include "AliJHistManager.h"
#include <TMath.h>
using namespace std;
//////////////////////////////////////////////////////
//  AliJBin
//////////////////////////////////////////////////////

AliJNamed::AliJNamed(TString name, TString title, TString opt, int mode) :
  fName(name),
  fTitle(title),
  fOption(opt),
  fMode(mode)
{
    // constructor
}

AliJNamed::~AliJNamed(){
    // virtual destructor for base class
}

TString AliJNamed::GetOption(TString key){
    TPMERegexp a("&"+key+"=?([^&]*)","i");
    int nMatch = a.Match(fOption);
    if( nMatch < 2 ) return UndefinedOption();
    return a[1];
}
void AliJNamed::SetOption( TString key, TString value){
    TPMERegexp a("&"+key+"=?[^&]*","i");
    int nMatch = a.Match(fOption);
    TString newOpt = "&"+key +( value.Length()?"="+value:"");
    if( value == UndefinedOption() ) newOpt = "";
    if( nMatch < 1 ) fOption += newOpt;
    else fOption.ReplaceAll( a[0], newOpt );
}
void AliJNamed::RemoveOption( TString key ){
    SetOption( key, UndefinedOption() );
}
TString AliJNamed::UndefinedOption(){
    //static TString undefinedOption = "Undefined";
    //return undefinedOption;
    return "Undefined";
}

//////////////////////////////////////////////////////
//  AliJBin
//////////////////////////////////////////////////////

//_____________________________________________________
AliJBin::AliJBin():
    AliJNamed("AliJBin","%.2f-%2.f", "&Mode=Range", kRange),
    fBinD(0),
    fBinStr(0),
    fIsFixedBin(false), 
    fIndexName("H"), 
    fHMG(NULL)
{;}
//_____________________________________________________
AliJBin::AliJBin(TString config, AliJHistManager * hmg):
    AliJNamed("AliJBin","%.2f-%2.f", "&Mode=Range", kRange),
    fBinD(0),
    fBinStr(0),
    fIsFixedBin(false), 
    fIndexName("H"), 
    fHMG(NULL)
{
    cout<< config<<endl;
    std::vector<TString> t = Tokenize(config, " \t,");
    TString type = t[0];
    SetName( t[1] );
    fIndexName =  t[2] ;
    SetTitle( t[3] );
    fTitle.ReplaceAll("\"","" );
    SetFullOption( t[4] );
    fMode = GetMode( GetOption("mode") );
    AddToManager( hmg );
    TString s;
    for( int i=5;i<int(t.size());i++ ) s+=" "+t[i];
    SetBin(s);
}

//_____________________________________________________
AliJBin::AliJBin(const AliJBin& obj) :
    AliJNamed(obj.fName,obj.fTitle,obj.fOption,obj.fMode),
    fBinD(obj.fBinD),
    fBinStr(obj.fBinStr),
    fIsFixedBin(obj.fIsFixedBin), 
    fIndexName(obj.fIndexName), 
    fHMG(obj.fHMG)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
AliJBin& AliJBin::operator=(const AliJBin& obj)
{
  // assignment operator
  if(this != &obj){
    // TODO: proper implementation
  }
  return *this;
}

//_____________________________________________________
void AliJBin::FixBin(){ 
    if(fIsFixedBin ) return;
    fIsFixedBin = true;
    if(!fHMG) AddToManager( AliJHistManager::CurrentManager());
    
}
//_____________________________________________________
void AliJBin::AddToManager( AliJHistManager* hmg ){ 
    hmg->Add(this); 
}
//_____________________________________________________
AliJBin & AliJBin::Set( TString name, TString iname, TString title, int mode){
    SetNameTitle( name, title );
    fIndexName = iname;
    fMode = mode;
    SetOption("mode",GetModeString(mode));
    return *this;
}
//_____________________________________________________
TString AliJBin::GetModeString(int i){
    static TString mode[] = { "Single","Range","String" };
    if( i<0 || i>2 ) return "";
    return mode[i];
}
int AliJBin::GetMode( TString mode ){
    for( int i=0;i<kNMode;i++ ) if( mode == GetModeString(i) ) return i;
    return -1;
}
//_____________________________________________________
AliJBin & AliJBin::SetBin( const int n, const float * v ){
    for( int i=0;i<n;i++ ) {
        AddBin( v[i] );
    }
    FixBin();
    return *this;
}
//_____________________________________________________
AliJBin & AliJBin::SetBin( const int n, const double * v ){
    for( int i=0;i<n;i++ ) {
        AddBin( v[i] );
    }
    FixBin();
    return *this;
}
AliJBin & AliJBin::SetBin( TVector *v ){
    for( int i=0;i<v->GetNrows();i++ ) {
        AddBin( (v->GetMatrixArray())[i] );
    }
    FixBin();
    return *this;
}
//_____________________________________________________
AliJBin& AliJBin::SetBin(const TString  v){
    std::vector<TString> ar = Tokenize( v, "\t ,");
    for( UInt_t i=0; i<ar.size();i++ ) {
        AddBin( ar[i] );
    }
    FixBin();
    return *this;
}

//_____________________________________________________
AliJBin& AliJBin::SetBin(const int  n){
    for( UInt_t i=0; i<UInt_t(n);i++ ) {
        AddBin( i );
    }
    FixBin();
    return *this;
}
//_____________________________________________________
void AliJBin::AddBin( const TString& v ){
    if( fIsFixedBin ) { JERROR( "You can't Add Bini "+GetName()); }
    fBinStr.push_back( (v=="_")?"":v );
    fBinD.push_back( v.Atof() );
}
//_____________________________________________________
void AliJBin::AddBin( float v ){
    if( fIsFixedBin ) { JERROR( "You can't Add Bin "+GetName()); }
    fBinD.push_back( v ); 
    fBinStr.push_back(Form("%f",v));
}
//_____________________________________________________
TString AliJBin::BuildTitle( int i ){
    if( i < 0 || i > Size() ) return "";
    if( fMode == kSingle )
        return TString(Form(fTitle.Data(), fBinD[i] ));
    if( fMode == kRange )
        return TString(Form(fTitle.Data(), fBinD[i], fBinD[i+1]));
    if( fMode == kString )
        return TString( Form(fTitle.Data(), fBinStr[i].Data()) );
    JERROR( TString("Bad Mode of AliJBin type ") + char(fMode) + " in " + fName+ "!!!" );
    return "";
}
//_____________________________________________________
TString AliJBin::GetString(){
    SetOption( "mode",GetModeString(fMode) );
    return "AliJBin\t"+fName+"\t"+fIndexName+"\t\""+fTitle+"\""+"\t"+fOption+"\t"+Join(fBinStr," ");

}
//_____________________________________________________
void AliJBin::Print(){
    std::cout<<"*"+GetString()<<std::endl;
}

int AliJBin::GetBin(double x){
  int i =  TMath::BinarySearch( fBinD.size(), &fBinD[0], x ); 
  if( fMode == kRange && i+1 >= int(fBinD.size()) ) return -1;
  return i;
}

//////////////////////////////////////////////////////
// AliJArrayBase 
//////////////////////////////////////////////////////

//_____________________________________________________
AliJArrayBase::AliJArrayBase():
    AliJNamed("AliJArayBase","","&Dir=default&Lazy",0),
    //AliJNamed("AliJArayBase","","&Dir=default&LessLazy",0),
    fDim(0),
    fIndex(0),
    fArraySize(0),
    fNGenerated(0),
    fIsBinFixed(false),
    fIsBinLocked(false),
    fAlg(NULL)
{
  // constrctor
}
//_____________________________________________________
AliJArrayBase::~AliJArrayBase(){
    //destructor
    if(fAlg) delete fAlg;
}

//_____________________________________________________
AliJArrayBase::AliJArrayBase(const AliJArrayBase& obj) :
    AliJNamed(obj.fName,obj.fTitle,obj.fOption,obj.fMode),
    fDim(obj.fDim),
    fIndex(obj.fIndex),
    fArraySize(obj.fArraySize),
    fNGenerated(obj.fNGenerated),
    fIsBinFixed(obj.fIsBinFixed),
    fIsBinLocked(obj.fIsBinLocked),
    fAlg(obj.fAlg)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
AliJArrayBase& AliJArrayBase::operator=(const AliJArrayBase& obj)
{
  // assignment operator
  if(this != &obj){
    // TODO: proper implementation
  }
  return *this;
}

//_____________________________________________________
void* AliJArrayBase::GetItem(){
    void * item = fAlg->GetItem();
    if( !item ){ 
        BuildItem() ; 
        item = fAlg->GetItem();
    }
    return item;
}
//_____________________________________________________
void* AliJArrayBase::GetSingleItem(){
    if(fMode == kSingle )return GetItem();
    JERROR("This is not single array");
    return NULL;
}
//_____________________________________________________
void AliJArrayBase::FixBin(){
    if( Dimension() == 0 ){
        AddDim(1);SetOption("Single");
        fMode = kSingle;
        if( HasOption("dir","default")) RemoveOption("dir");
    }
    ClearIndex();
    fAlg = new AliJArrayAlgorithmSimple(this);
    fArraySize = fAlg->BuildArray();
}
//_____________________________________________________
int AliJArrayBase::Index(int d){
    if( OutOfDim(d) ) JERROR("Wrong Dim");
    return fIndex[d];
}
void AliJArrayBase::SetIndex(int i, int d ){
    if( OutOfSize( i, d ) ) JERROR( "Wrong Index" );
    fIndex[d] = i;
}

void AliJArrayBase::InitIterator(){ fAlg->InitIterator(); }
bool AliJArrayBase::Next(void *& item){ return fAlg->Next(item); }


//////////////////////////////////////////////////////
//  AliJArrayAlgorithm
//////////////////////////////////////////////////////

//_____________________________________________________
AliJArrayAlgorithm::AliJArrayAlgorithm(AliJArrayBase * cmd):
    fCMD(cmd)
{
  // constructor
}
//_____________________________________________________
AliJArrayAlgorithm::~AliJArrayAlgorithm(){
  // destructor
}

//_____________________________________________________
AliJArrayAlgorithm::AliJArrayAlgorithm(const AliJArrayAlgorithm& obj) :
    fCMD(obj.fCMD)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
AliJArrayAlgorithm& AliJArrayAlgorithm::operator=(const AliJArrayAlgorithm& obj)
{
  // assignment operator
  if(this != &obj){
    *fCMD = *(obj.fCMD);
  }
  return *this;
}

//////////////////////////////////////////////////////
//  AliJArrayAlgorithmSimple
//////////////////////////////////////////////////////

//_____________________________________________________
AliJArrayAlgorithmSimple::AliJArrayAlgorithmSimple(AliJArrayBase * cmd):
    AliJArrayAlgorithm(cmd),
    fDimFactor(0),
    fArray(NULL),
    fPos(0)
{
  // constructor
}
//_____________________________________________________
AliJArrayAlgorithmSimple::~AliJArrayAlgorithmSimple(){
    // Dimension, GetEntries, SizeOf
    if( fArray ) delete [] (void**)fArray;
}

//_____________________________________________________
AliJArrayAlgorithmSimple::AliJArrayAlgorithmSimple(const AliJArrayAlgorithmSimple& obj) :
    AliJArrayAlgorithm(obj.fCMD),
    fDimFactor(obj.fDimFactor),
    fArray(obj.fArray),
    fPos(obj.fPos)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
AliJArrayAlgorithmSimple& AliJArrayAlgorithmSimple::operator=(const AliJArrayAlgorithmSimple& obj)
{
  // assignment operator TODO: proper implementation
  if(this != &obj){
    *fCMD = *(obj.fCMD);
  }
  return *this;
}
//_____________________________________________________
int AliJArrayAlgorithmSimple::BuildArray(){
    fDimFactor.resize( Dimension(), 1 );
    for( int i=Dimension()-2; i>=0; i-- ){
        fDimFactor[i] = fDimFactor[i+1] * SizeOf(i+1);
    } // TODO split to BuildArray and lazyArray in GetItem
    int arraySize = fDimFactor[0] * SizeOf(0);
    fArray = new void*[arraySize];
    for( int i=0;i<arraySize;i++ ) fArray[i] = NULL;
    return arraySize;
}
//_____________________________________________________
int  AliJArrayAlgorithmSimple::GlobalIndex(){
    int iG = 0;
    for( int i=0;i<Dimension();i++ ) // Index is checked by fCMD
        iG+= Index(i)*fDimFactor[i];
    // TODO check iG
    return iG;

}
void AliJArrayAlgorithmSimple::ReverseIndex( int iG ){
    int n = iG;
    for( int i=0;i<Dimension();i++ ){
        int n1 = int(n/fDimFactor[i]);
        fCMD->SetIndex( n1 , i );
        n-=n1*fDimFactor[i];
    }
}
void * AliJArrayAlgorithmSimple::GetItem(){
    return fArray[GlobalIndex()];

}
void AliJArrayAlgorithmSimple::SetItem(void * item){
    fArray[GlobalIndex()] = item;
}


//////////////////////////////////////////////////////
//  AliJTH1
//////////////////////////////////////////////////////
//_____________________________________________________
AliJTH1::AliJTH1():
    fDirectory(NULL),
    fSubDirectory(NULL),
    fHMG(NULL),
    fTemplate(NULL),
    fBins(0)
{
    // default constructor
    fName="AliJTH1";
}

//_____________________________________________________
AliJTH1::AliJTH1(TString config, AliJHistManager * hmg):
    fDirectory(NULL),
    fSubDirectory(NULL),
    fHMG(NULL),
    fTemplate(NULL),
    fBins(0)
{
    // constructor
    std::vector<TString> t = Tokenize(config, " \t,");
    TString type = t[0];
    SetName( t[1] );
    SetTitle( t[2] );
    fTitle.ReplaceAll("\"","");
    SetFullOption( t[3] );
    fMode = HasOption("mode","Single")?kSingle:kNormal;
    AddToManager( hmg );
    TString s;
    for( int i=4;i<int(t.size());i++ ) s+=" "+t[i];
    AddDim( s );
    FixBin();
}
//_____________________________________________________
AliJTH1::~AliJTH1(){
    // destructor
    if( fNGenerated == 0 && fTemplate ) delete fTemplate;
}

//_____________________________________________________
AliJTH1::AliJTH1(const AliJTH1& obj) :
    AliJArrayBase(),
    fDirectory(obj.fDirectory),
    fSubDirectory(obj.fSubDirectory),
    fHMG(obj.fHMG),
    fTemplate(obj.fTemplate),
    fBins(obj.fBins)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
AliJTH1& AliJTH1::operator=(const AliJTH1& obj)
{
  // assignment operator
  if(this != &obj){
    // TODO: proper implementation
  }
  return *this;
}

//_____________________________________________________
int AliJTH1::AddDim( AliJBin * bin){
    int ndim = this->AliJArrayBase::AddDim( bin->Size() );
    fBins.resize( ndim, NULL );
    fBins[ndim-1] = bin;
    return ndim;
}

int AliJTH1::AddDim(TString v) {
    if( v == "END" ) { FixBin(); }
    else{
        std::vector<TString> o= Tokenize(v, "\t ,");
        for( UInt_t i=0;i<o.size();i++ ){
            TString & s = o[i];
            if(s.Length() == 0 ) continue;
            if( s.IsFloat() ) { // TODO IsInt? IsDigit?
                AddDim( s.Atoi() );
                continue;
            }
            AliJBin * b = NULL;
            if( fHMG ) b = fHMG->GetBin(s);
            if( b ) this->AddDim(b);
            else {JERROR("Wrong terminator of Array : \"" + s+"\" in " + fName ); }
        }
    }
    return Dimension();
}
//_____________________________________________________
Int_t AliJTH1::Write(){
    TDirectory *owd = (TDirectory*) gDirectory;
    InitIterator();
    void * item;
    if( fSubDirectory ) fSubDirectory->cd();
    //else fDirectory->cd();
    while( Next(item) ){
        if( !item ) continue;
        TH1 * obj = static_cast<TH1*>(item);
        obj->Write( );
        //obj->Write( 0, TObject::kOverwrite );
    }
    if( owd != gDirectory ) owd->cd();
    return 0;
}
//_____________________________________________________
TString AliJTH1::GetString( ){
    TString s = Form( "%s\t%s\t\"%s\"\t%s\t", 
            ClassName(), fName.Data(), fTitle.Data(), fOption.Data() );
    for( int i=0;i<Dimension();i++ ){
        if( int(fBins.size()) > i && fBins[i] != NULL ){
            s+= " "+fBins[i]->GetName();
        }else{
            s+= TString(" ")+Form("%d", SizeOf(i));
        }
    }
    return s;
}
//_____________________________________________________
void AliJTH1::FixBin(){
    this->AliJArrayBase::FixBin();

    if(!fHMG) {
        AddToManager( AliJHistManager::CurrentManager() );
    }
    if(!fDirectory) fDirectory = fHMG->GetDirectory();
}

//_____________________________________________________
void AliJTH1::AddToManager(AliJHistManager *hmg){
    if(fHMG) return; // TODO handle error
    fHMG = hmg;
    hmg->Add(this);
}
//_____________________________________________________
void AliJTH1::Print(){
    std::cout<<"*"<<GetString()<<std::endl;
    // TODO more details.
}
//_____________________________________________________
void AliJTH1::SetTemplate(TH1 *h){
    if( fTemplate ) return; /// TDOO give error
    fTemplate = (TH1*)h->Clone();
    fTemplate->Sumw2();
    fTemplate->SetDirectory(0);
    fName = h->GetName();
    fTitle = h->GetTitle();
}
//_____________________________________________________
TString AliJTH1::BuildName(){
    TString name = fName;
    if( !HasOption("Single") )
        for( int i=0;i<Dimension();i++ ){
            name+=((int(fBins.size()) > i && fBins[i] != NULL)?fBins[i]->GetIndexName():"H")
                +Form("%02d",Index(i));
        }
    return name;
}
//_____________________________________________________
TString AliJTH1::BuildTitle(){
    TString title = fTitle;
    for( int i=0;i<Dimension();i++ )
        title+=((int(fBins.size()) > i && fBins[i] != NULL)?" "+fBins[i]->BuildTitle(Index(i)):"")
            +Form("%02d",Index(i));
    return title;
}
//_____________________________________________________
void * AliJTH1::BuildItem(){
    TDirectory * owd = (TDirectory*) gDirectory;
    gROOT->cd();
    TString name = BuildName();
    TH1 * item = NULL;
    if( !fSubDirectory ){
        if( !HasOption("dir") ) {
            fSubDirectory = fDirectory;
        }
        else {
            fSubDirectory = fDirectory->GetDirectory(fName);
            if( !fSubDirectory && !IsLoadMode() ){
                fSubDirectory = fDirectory->mkdir(fName);
            }
        }
    }
    if( IsLoadMode() ){
        //if( fSubDirectory ) JDEBUG(2, fSubDirectory->GetName() );
        if( fSubDirectory )
            item = dynamic_cast<TH1*>(fSubDirectory->Get( name ));
        if( !item ){
            void ** rawItem = fAlg->GetRawItem();
            InitIterator();
            void * tmp;
            while( Next(tmp) ){
                item = (TH1*)fSubDirectory->Get(BuildName());
                if( item ) break;
            }
            if( item ) {
                item = dynamic_cast<TH1*>((static_cast<TH1*>(item))->Clone(name));
                if( !item ){ JERROR("Any of "+fName+" doesn't exists. I need at least one"); return NULL;}
                item->Reset();
                item->SetTitle( BuildTitle() );
                item->SetDirectory(0);
                *rawItem = (void*)item;
            }
        }
        if( !item ){ JERROR("Any of "+fName+" doesn't exists. I need at least one"); return NULL;}
    }
    else{ //  Gen Mode
        TH1 * titem = NULL;
        if(fNGenerated == 0 ) {
            titem = fTemplate;
        }
        else titem =(TH1*) fTemplate->Clone();
        titem->SetDirectory( fSubDirectory );
        titem->Reset();
        titem->SetName( BuildName() );
        titem->SetTitle( BuildTitle() );
        fNGenerated++;
        item=titem;
    }
    if( item ) fAlg->SetItem(item);
    owd->cd();
    return (void*)item;
}
//_____________________________________________________
bool AliJTH1::IsLoadMode(){
    return fHMG->IsLoadMode();
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJTH1Derived                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
template< typename T>
AliJTH1Derived<T>::AliJTH1Derived():
    AliJTH1(), fPlayer(this)
{
}
template< typename T>
AliJTH1Derived<T>::~AliJTH1Derived(){
}




//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliJHistManager                                                       //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
AliJHistManager::AliJHistManager(TString name, TString dirname):
    AliJNamed(name,"","",0),
    fIsLoadMode(false),
    fDirectory(gDirectory),
    fConfigStr(),
    fBin(0),
    fHist(0),
    fManager(0),
    fBinNames(0),
    fBinConfigs(0),
    fHistNames(0),
    fHistConfigs(0)
{
    // constructor
    if( dirname.Length() == 0 ) dirname = name;
    if( dirname.Length() > 0 ) {
        fDirectory = (TDirectory*)gDirectory->Get(dirname);
        if( fDirectory ){
            std::cout<<"JWARNING : "<<Form("Hist directory %s exists", dirname.Data() )<<std::endl;
            // gSystem->Exit(1);  // We might actually want the directory to exist, so no exit
        }
        if( !fDirectory ){
            fDirectory = gDirectory->mkdir( dirname );
        }
    }
    if( !fDirectory ){
        std::cout<<"JERROR : "<<Form("Fail to generate Hist directory %s", dirname.Data() )<<std::endl;
        gSystem->Exit(1);
    }
    this->cd();
}

//_____________________________________________________
AliJHistManager::AliJHistManager(const AliJHistManager& obj) :
    AliJNamed(obj.fName,obj.fTitle,obj.fOption,obj.fMode),
    fIsLoadMode(obj.fIsLoadMode),
    fDirectory(obj.fDirectory),
    fConfigStr(obj.fConfigStr),
    fBin(obj.fBin),
    fHist(obj.fHist),
    fManager(obj.fManager),
    fBinNames(obj.fBinNames),
    fBinConfigs(obj.fBinConfigs),
    fHistNames(obj.fHistNames),
    fHistConfigs(obj.fHistConfigs)
{
    // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
AliJHistManager& AliJHistManager::operator=(const AliJHistManager& obj)
{
    // assignment operator
    if(this != &obj){
        // TODO: proper implementation
    }
    return *this;
}

AliJHistManager* AliJHistManager::GlobalManager(){
    static AliJHistManager* singleton = new AliJHistManager("GlobalHistManager");
    return singleton;
}

AliJHistManager* AliJHistManager::CurrentManager( AliJHistManager * hmg){
    static AliJHistManager* currentManager = NULL;//;AliJHistManager::GlobalManager();
    if( hmg ) currentManager = hmg; 
    return currentManager;
}

AliJBin* AliJHistManager::GetBuiltBin(TString s ){
    for( int i=0;i<int(fBin.size());i++ )
        if( fBin[i]->GetName() == s ) return fBin[i];
    return NULL;
}
AliJBin* AliJHistManager::GetBin(TString s ){
    AliJBin* h = GetBuiltBin(s);
    if(h) return h;
    for( int i=0;i<GetNBin();i++ )
        if( fBinNames[i] == s ){
            return new AliJBin( fBinConfigs[i],this );
        }
    return NULL;
}
AliJTH1 * AliJHistManager::GetBuiltTH1(TString s ){
    for( int i=0;i<int(fHist.size());i++ )
        if( fHist[i]->GetName() == s ) return fHist[i];
    return NULL;
}
// Note: Returning NULL crashes the code, something should be done about this.
// The error given by compiler is: non-const lvalue reference to type 'AliJTH1D' (aka 'AliJTH1Derived<TH1D>') cannot bind to a temporary of type 'long'
// The reoson for crash is that NULL cannot be used in dynamic_cast<AliJTH1D&>
AliJTH1 * AliJHistManager::GetTH1(TString s ){
    AliJTH1 * h = GetBuiltTH1(s);
    if( h ) return h;
    for( int i=0;i<GetNHist();i++ )
        if( fHistNames[i] == s ){
            if( fHistConfigs[i].BeginsWith("AliJTH1D")) return new AliJTH1D( fHistConfigs[i], this );
            if( fHistConfigs[i].BeginsWith("AliJTH2D")) return new AliJTH2D( fHistConfigs[i], this );
            if( fHistConfigs[i].BeginsWith("AliJTH3D")) return new AliJTH3D( fHistConfigs[i], this );
            if( fHistConfigs[i].BeginsWith("AliJTProfile")) return new AliJTProfile( fHistConfigs[i], this );
        }
    return NULL;
}
void AliJHistManager::Add(AliJBin *o ){
    if( !o ) return;
    if( GetBuiltBin( o->GetName() ) ) return; // TODO error handle
    fBin.push_back( o ); 
}
void AliJHistManager::Add(AliJTH1 *o ){
    if( !o ) return;
    if( GetBuiltTH1( o->GetName() ) ) return; // TODO error handle
    fHist.push_back( o ); 
}
void AliJHistManager::Print(){
    if( IsLoadMode() ) {
        cout<<fConfigStr<<endl;
        return;
    }
    cout<<"============ AliJHistManager : "<<fName<<" ==================="<<endl;
    cout<<endl;
    cout<<"---- AliJBin ----"<<endl;
    for( int i=0;i<GetNBin();i++ ){
        fBin[i]->Print();
    }
    cout<<endl;
    cout<<"---- AliJTH1 ----"<<endl;
    for( int i=0;i<GetNHist();i++ ){
        fHist[i]->Print();
    }
}
void AliJHistManager::Write(){
    for( int i=0;i<GetNHist();i++ )
        fHist[i]->Write();
}

void AliJHistManager::WriteConfig(){
    TDirectory *owd = fDirectory;
    //cout<<"DEBUG_T1: "<<fDirectory<<endl;
    //cout<<"DEBUG_T2: "<<GetName()<<"\t"<<fDirectory->GetName()<<endl;
    //exit(1);
    // TODO 1.Error Check 2.Duplicaition check
    TDirectory * fHistConfigDir = fDirectory->mkdir("HistManager");
    fHistConfigDir->cd();
    TObjString * config = new TObjString(GetString().Data());
    config->Write("Config");
    owd->cd();
}

int AliJHistManager::LoadConfig(){
    SetLoadMode(true);
    TObjString *strobj = (TObjString*)fDirectory->Get("HistManager/Config");
    if( !strobj ) return 0; //TODO
    TString config = strobj->String();
    fConfigStr = config;
    vector<TString> lines = Tokenize(config, "\n");
    cout<< Form("Read Config.%d objects found\n", (int)lines.size() );
    for( UInt_t i=0;i < lines.size();i++ ){
        TString line = lines.at(i);
        std::vector<TString> t = Tokenize(line, " \t,");
        if(line.BeginsWith("AliJBin")) {
            fBinNames.push_back( t[1] );
            fBinConfigs.push_back( line );
        }else if( line.BeginsWith("AliJ")){
            fHistNames.push_back( t[1] );
            fHistConfigs.push_back( line );
        }
    }
    return 1;
}

bool AliJHistManager::HistogramExists(TString name){
  for(int i = 0; i < fHistNames.size(); i++){
    if(fHistNames[i] == name) return true;
  }
  return false;
}


//////////////////////////////////////////////////////
//  Utils
//////////////////////////////////////////////////////
vector<TString> Tokenize( TString s, TString d, int quote ){
    //int nd = d.Length();
    bool flagBeforeToken = 0;
    bool inQuote = 0;
    TString tok="";
    vector<TString> toks;
    s += d[0];
    for( int i=0;i<s.Length();i++ ){
        if( quote == 1 && s[i] == '\"' ){ inQuote = ! inQuote; }
        if( d.First(s[i]) != kNPOS && !inQuote ){
            if( flagBeforeToken == 0 && tok.Length()>0){
                toks.push_back(tok);
                tok.Clear();
                flagBeforeToken = 1;
            }
        }else{
            tok+=s[i];
            flagBeforeToken = 0;
        }
    }
    return toks;
}


TString Join( vector<TString>& ss , TString del){
    if( ss.size() < 1 ) return "";
    TString s = ss[0];
    for( UInt_t i=1;i<ss.size();i++ ) s+=del+ss[i];
    return s;
}

template class AliJTH1Derived<TH1D>;
template class AliJTH1Derived<TH2D>;
template class AliJTH1Derived<TH3D>;
template class AliJTH1Derived<TProfile>;

bool OutOf( int i, int x, int y ){ return ( i<x || i>y ); }

#include <TFile.h>

void ttestAliJArray(){
    AliJHistManager * fHMG;
    AliJBin fCentBin;
    AliJBin fVtxBin;
    AliJBin fPTtBin;
    AliJBin fPTaBin;
    AliJBin fXEBin;
    AliJBin fKLongBin;
    AliJBin fRGapBin;
    AliJBin fEtaGapBin;
    AliJBin fPhiGapBin;
    AliJBin fMassBin;
    AliJBin fTypBin;
    AliJBin fTypBin3;
    AliJBin fPairPtBin;
    AliJTH1D fhTriggPtBin;
    AliJTH1D fhTriggMult;
    AliJTH1D fhIphiTrigg;
    AliJTH1D fhIetaTrigg;
    AliJTH2D test1;
    AliJTProfile test2;

    TFile * f = new TFile("test.root","RECREATE");
    fHMG = AliJHistManager::GlobalManager();
    fCentBin   .Set("Cent",   "C", "C %2.0f-%2.0f%%" ).SetBin( "0 100");
    fVtxBin    .Set("Vtx",    "V", "Vtx %.0f-%.0f" ).SetBin("-10 10");
    fPTtBin    .Set("PTt",    "T", "p_{Tt} %.1f-%.1f").SetBin("3 5 8 10 15 20");
    fPTaBin    .Set("PTa",    "A", "p_{Tt} %.1f-%.1f").SetBin("3 5 8 10 15 20");

    fhTriggMult
        << TH1D( "hTriggMult", "",  100, -0.5, 99.5) 
        <<  fCentBin << fPTtBin  << "END";
    fhIphiTrigg
        << TH1D( "fhIphiTrigg", "",  3, -0.1, 0.1 ) 
        <<  fCentBin << fPTtBin  << "END";
    fhIetaTrigg
        << TH1D( "hIetaTrigg", "",  80, -5, 5 ) 
        <<  fCentBin << fPTtBin  << "END";// inclusive eta
    fhTriggPtBin
        << TH1D( "hTriggPtBin", "", 10,0,10 )
        <<  fCentBin << fVtxBin << fPTtBin  << "END";

    fhTriggMult[0][0]->Fill(1);
    fHMG->Print();

    f->Write();
    fHMG->Write();
    fHMG->WriteConfig();

}
