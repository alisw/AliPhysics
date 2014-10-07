//====================================
//last modified FK 6.NOV 2009
//====================================
//blah
// blah

#include "AliJCard.h"

//ClassImp(AliJCard);

AliJCard::AliJCard():
    AliJBaseCard(),
    fhCorr(0),
    feventV3kv(0),
    fIndexVector(0),
    fpi0massbin(0)
{   
    //constructor
}

AliJCard::AliJCard(const char *filename):
    AliJBaseCard(filename),
    fhCorr(0),
    feventV3kv(0),
    fIndexVector(0),
    fpi0massbin(0)
{  
    //constructor
    InitCard();
    MakeFastCorrTypeIndex();
}



AliJCard::AliJCard(const AliJCard& obj) : 
    AliJBaseCard(),
    fhCorr(obj.fhCorr),
    feventV3kv(obj.feventV3kv),
    fIndexVector(obj.fIndexVector),
    fpi0massbin(obj.fpi0massbin)
{
    // copy constructor
    JUNUSED(obj);
}


AliJCard& AliJCard::operator=(const AliJCard& obj){
    // copy constructor
    JUNUSED(obj);
    return *this;
}

AliJCard::~AliJCard(){
    // destructor
    if( fpi0massbin ){
        for( int i = 0; i < GetNoOfBins( kCentrType ); i++ ){
            delete [] fpi0massbin[i];
        }    // destructor
        delete [] fpi0massbin;
    }
    if( fhCorr ){
        delete fhCorr;
    }
}

void AliJCard::MakeFastCorrTypeIndex(){
    // make fast findex array
    for( unsigned int i=0;i<fKeyWordVector.size();i++ ){
        int corrIndex = GetCorrType( fKeyWordVector[i] );
        if( corrIndex != kNoType )
            fIndexVector[corrIndex] = i;
    }

}

int AliJCard::IsLessThanUpperPairPtCut(double inPairPt){ 
    // pt cut
    int nB = GetN("UpperPairPtCut");
    if(inPairPt == -999) return nB;
    if(-999 < inPairPt && inPairPt <=0) return int(Get("UpperPairPtCut",-int(inPairPt)));
    int bin=-1, i=0;
    while(i<nB && inPairPt>Get("UpperPairPtCut",i)) i++;
    if(i<nB) bin=i;
    //cout<<" i="<<i<<" bin="<<bin<<" inval="<<inPairPt<<endl;;
    return bin;
}



TString AliJCard::GetKeyWord(corrType ctype){
    // get keyword
    TString kw;

    switch(ctype){
        case kTriggType:  kw = "TriggPtBorders"; break;
        case kAssocType:  kw = "AssocPtBorders"; break;
        case kXeType:     kw = "xEBorders"; break;
        case kLongType:   kw = "KlongBorders"; break;
        case kCentrType:  kw = "CentBinBorders"; break;
        case kZVertType:  kw = "zVertBins"; break;
        case kMassType:   kw = "PairInvariantMassBins"; break;
        case kEtaGapType: kw = "EtaGapThresholds"; break;
        case kDiJetType:  kw = "DiJetMassBorders"; break;
        case kRGapType:   kw = "RGapThresholds"; break;
                          //case kEPType:  kw = "EPBorders"; break;
        default : cout<<"ERROR: kNoType on input to AliJCard::GetKeyWord"<<endl; exit(1); break;
    }
    return kw; 
}

corrType AliJCard::GetCorrType( TString inStr ){
    // get corr type
    TString kw;

    int i;
    TString s;

    // go through all corrType to check if the parameter should be added
    // into fast findex array
    for( i = 0; i < kNcorrType; i++ ){
        if( i == kNoType )
            continue;

        s = GetKeyWord( (corrType)i );
        if( ! strcmp( s.Data(), inStr.Data() ))
            return (corrType)i;
    }

    return kNoType; 
}

int AliJCard::GetN(corrType ctype){
    //returns size of TVector
    return GetN(GetKeyWord(ctype));
}

int AliJCard::GetNFast(corrType ctype){
    //returns size of TVector
    int findex = fIndexVector[ctype];
    if( findex > -1 ){
        return (int) fValuesVector[findex].GetNrows();
    }else{
        cout<<"ERROR: fValuesVector fast findex out of range "<< findex << " " << (GetKeyWord( ctype )).Data() << endl;
        exit(1);
    }
}


float AliJCard::Get(corrType ctype, int VectorComponent){
    //returns VectorComponent Component of  fValuesVecto`uor TVector for given keyword
    return Get(GetKeyWord(ctype), VectorComponent);
}

float AliJCard::GetFast(corrType ctype, int VectorComponent){
    // fast get

    if(0<=VectorComponent && VectorComponent<GetNFast(ctype) && fIndexVector[ctype] > -1 ){
        return fValuesVector[fIndexVector[ctype]][VectorComponent+1];
    }else{
        cout<<"ERROR: fValuesVector fast findex out of range "<< (GetKeyWord(ctype)).Data()<<endl;
        exit(1);
    }
}

void AliJCard::ReCompile(){
    InitCard();
    MakeFastCorrTypeIndex();
}


void AliJCard::InitCard(){
    // Init card
    cout<<"Init of AliJCard"<<endl;
    // set the length of fIndexVector and disable all indices
    fIndexVector.resize( kNcorrType  );
    for( int i = 0; i < kNcorrType; i++ )
        fIndexVector[i] = -1;
}

void AliJCard::FinishCard(){
    // Finish loading of card
    AliJBaseCard::FinishCard();
    // recompute fast idices

    fpi0massbin = new Double_t*[GetNoOfBins( kCentrType )];

    char centstr[10];

    // read the pi0 mass bin parameters
    for( int i = 0; i < GetNoOfBins( kCentrType ); i++ ){
        fpi0massbin[i] = new Double_t[8];

        for( int j = 0; j < 8; j++ ){
            sprintf( centstr, "cent%d", i );
            fpi0massbin[i][j] = Get( centstr, j );
        }
    }

}

void AliJCard::PrintOut(){
    // Print out contents of card
    AliJBaseCard::PrintOut();
    cout << "----- fast array ----------" << endl;
    for( unsigned int ii = 0; ii < fIndexVector.size(); ii++ ){
        if( ii == kNoType )
            continue;
        cout << (GetKeyWord( (corrType)ii )).Data() << " " << fIndexVector[ii];
        if( fIndexVector[ii] > -1 ){
            cout<<" (dim ="<<fValuesVector[fIndexVector[ii]].GetNrows()<<") ";//print size of TVector
            for(int j=1; j<=fValuesVector[fIndexVector[ii]].GetNrows(); j++){
                cout<<fValuesVector[fIndexVector[ii]][j]<<" ";//TVector components
            }
        }
        else
            cout << " no link!";

        cout << endl;  
    }
}

int AliJCard::GetBin(corrType ctype, float val){
    // get bin

    if(ctype == kNoType) return 0;

    TVector * v = &fValuesVector[fIndexVector[ctype]];
    int iBin2 = TMath::BinarySearch( v->GetNrows(), v->GetMatrixArray(), val );
    if( iBin2 >= v->GetNrows()-1 ) iBin2 = -1;

    return iBin2;
}

int AliJCard::GetBinFast(corrType ctype, float val){
    // fast get bin

    if(ctype == kNoType) return 0;

    for(int i=0; i<(GetNFast(ctype)-1); i++)
        if(GetFast(ctype,i)<=val && val<GetFast(ctype,i+1))
            return i;

    return -1;
}



bool AliJCard::IsGoodRun(int runID){
    // run quality
    bool isgood = true;
    for(int i=0; i<GetN("badRuns");i++) if(((int) Get("badRuns",i)) == runID) isgood = false;

    return isgood;
}

bool AliJCard::MbTrigger(int triggin) const {
    // MB trigger
    return (triggin & (1<<kMinBiasTriggerBitJCorran)) > 0; //masking on >0 masking off ==0 desired trigger
}


bool AliJCard::SimilarVertZ(float Z1, float Z2){
    // z vert
    static double v = Get("maxMixDZ");
    if( v < 0 ) return 1;
    return fabs(Z1-Z2) < v ;
}

bool AliJCard::SimilarMultiplicity(float mult1, float mult2){
    // multi
    static double v = Get("maxMixDMult");
    if( v <0) return 1;
    return fabs(mult1-mult2) < v;
}


bool AliJCard::SimilarCentrality(float c1, float c2, int cbin){
    // centra
    static TVector *v = GetVector("maxDCent");
    if(v==NULL) return 1;
    return fabs(c1-c2) < (*v)[cbin+1]; // TODO
}


//--------- P H E N I X    C G L --------------

bool AliJCard::InPhiRange(float Phi){
    // phi
    bool isIn = false; 
    for(int i=1; i<GetN("phiDCRange"); i+=2 ) 
        isIn = isIn || (Get("phiDCRange",i-1)<Phi && Phi<Get("phiDCRange",i));
    return isIn;
} 



bool  AliJCard::CheckTrackParamsInTPC(int NClustersTPC,float Chi2PerClusterTPC){
    // tpc pars
    bool isGoodTrack = true;
    if(NClustersTPC < Get("MinNClustersTPC"))           isGoodTrack = false;
    if(Chi2PerClusterTPC > Get("MaxChi2PerClusterTPC")) isGoodTrack = false;
    return isGoodTrack;
}

bool  AliJCard::CheckMinNumTPCClustPt(int NClustersTPC, float fpt){
    // tpc pars
    float minNTPCCls =  Get("ParMinNClustTPCPt",0)*log(Get("ParMinNClustTPCPt",1)*fpt + Get("ParMinNClustTPCPt",2)); 
    if(NClustersTPC > minNTPCCls) return true; //track had enough clusters
    else return false;  //track did not have sufficient number of clusters
}


bool  AliJCard::CheckTrackImpact(float xyIm, float zIm, float fpt){
    // tpc pars
    bool isGoodTrack = true;
    if(TMath::Abs(xyIm) > Get("ParMaxDCAToVertexXYPtDep",0)+Get("ParMaxDCAToVertexXYPtDep",1)/pow(fpt,Get("ParMaxDCAToVertexXYPtDep",2))){
        isGoodTrack = false;
    }
    if(TMath::Abs(zIm) > Get("MaxDCAToVertexZ")) isGoodTrack = false;

    return isGoodTrack;
}



bool AliJCard::DeltaEtaCheck(const AliJBaseTrack *ftk1, const AliJBaseTrack *ftk2) {
    // deta
    double delta = fabs( log(tan(ftk1->Theta()/2.0)) -  log(tan(ftk2->Theta()/2.0)) );
    if(Get("etaCut")>0) return delta < Get("etaCut");
    if(Get("etaCut")<0) return delta > Get("etaCut");

    return true;
}





















