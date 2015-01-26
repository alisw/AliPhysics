//====================================
//last modified FK 6.NOV 2009
//====================================
//blah
// blah

#include <TPRegexp.h>
#include "AliJConst.h"
#include "AliJBaseCard.h"

//ClassImp(AliJBaseCard);

AliJBaseCard::AliJBaseCard() :
  fnentry(0),
  fKeyWordVector(0),
  fValuesVector(0),
  fValueString(0),
  fKeyTable(0)
{   
  //constructor
}

AliJBaseCard::AliJBaseCard(const char *filename) :
  fnentry(0),
  fKeyWordVector(0),
  fValuesVector(0),
  fValueString(0),
  fKeyTable(0)
{  
  //constructor

  strcpy(fcardname, filename);  //needed in PrintOut()

  if( strlen( filename ) > 0 ){

    //----  r e a d   t h e   c a r d ----
    ReadInputCard();//read config file fill Tvectors
  }
}


AliJBaseCard& AliJBaseCard::operator=(const AliJBaseCard& obj){
  // equal operator
  JUNUSED(obj);
  return *this;
}

AliJBaseCard::~AliJBaseCard(){
  // destructor
}


unsigned int AliJBaseCard::GetTVectorIndex(TString keyword, int tol){

  //   //returns findex of a TVector according to its position in std::hash_map 
  //   std::hash_map< TString, unsigned int >::iterator iter = MapKeyWordToTVector.begin();
  //   iter = MapKeyWordToTVector.find(keyword);
  //   if(iter != MapKeyWordToTVector.end()){
  //     return (unsigned int) iter->second;
  //   }else{
  //     cout << "ERROR: \""<<keyword.Data()<<"\" must be defined  "<< endl;
  //     exit(1);
  //   }

  UInt_t i;
  Int_t ind;
  i = 0;
  ind = -1;

  //cout<<"ALIJBASECARD_SEARCH_MODE_HASHLIST"<<endl;
  TNamed * ko = (TNamed*)fKeyTable.FindObject( keyword.Data() ); 
  //cout<<ko<<endl;
  if(ko) ind = ko->GetUniqueID();
  if( ind == -1 ){
    for( UInt_t ii=0;ii<fKeyWordVector.size();ii++ ){
      if( fKeyWordVector[ii] == keyword ) return i;
    }
    if( tol == 0 ){
      cout << "ERROR: \""<<keyword.Data()<<"\" must be defined  "<< endl;
      exit(1);
    }else if ( tol == 1 ){
      cout << "Warning: \""<<keyword.Data()<<"\" is not exist. return default value  "<< endl;
      return -1;
    }else{
        return -1;
    }
  }

  return ind;
}

int AliJBaseCard::GetN(TString keyword){
  //returns size of TVector
  unsigned int findex = GetTVectorIndex(keyword);
  return (int) fValuesVector[findex].GetNrows();
}

TVector *  AliJBaseCard::GetVector(TString keyword ){
  int findex = GetTVectorIndex(keyword);
  return &fValuesVector[findex];
}

float AliJBaseCard::Get(TString keyword, int VectorComponent){
  //returns VectorComponent Component of  fValuesVector TVector for given keyword
  int findex = GetTVectorIndex(keyword);
  if(0<=VectorComponent && VectorComponent<GetNwithIndex(findex)){
    return fValuesVector[findex](VectorComponent+1);
  }else{
    cout<<"ERROR: fValuesVector findex out of range "<<keyword.Data()<<endl;
    cout << "   Max findex: " << GetN(keyword) -  1<< " Asked: " <<  VectorComponent << endl;
    exit(1);
  }
}

TString AliJBaseCard::GetStr(TString keyword ){
  int findex = GetTVectorIndex(keyword, 1);
  if( findex < 0  ) return TString("");
  return fValueString[findex];
}


void AliJBaseCard::InitCard(){
  // set the length of fIndexVector and disable all indices
}

void AliJBaseCard::FinishCard(){
  // recompute fast idices
}

void AliJBaseCard::ReadInputCard(){
  // read card

  char buffer[kMaxDimBuffer];
  ifstream incard;

  cout << "Reading fcard from file: " << fcardname << endl;
  incard.open(fcardname,ios::in);

  if(!incard){
    cout<<"ERROR: Config file <"<<fcardname<<"> not found!"<<endl;
    exit(1);
  }

  InitCard();

  while(!incard.eof()){ //loop over the input fcard

    incard.getline(buffer,kMaxDimBuffer); //read a line

    if(fnentry > 1000){//is the file reasonably long?
      cout<<"Maximum number of 1000 lines reached in AliJBaseCard.C"<<endl;
      exit(1);
    }

    ReadInputLine( buffer );

    fnentry++;
  }//while eof

  FinishCard();

  return;
}

void AliJBaseCard::ReadLine( const char * buffer ){
    TString tstr(buffer);
    TPMERegexp rsp(";");
    TPMERegexp csp1("=");
    TPMERegexp csp2(",");
    int nrow = rsp.Split( tstr );
    for( int i=0;i<nrow;i++ ){
        TString row = rsp[i];
        int nst = csp1.Split(row);
        if( nst!=2 ) continue; // TODO Error or warning
        TString key = csp1[0];
        TString val = csp1[1];
        int nc = csp2.Split( val );

        vector< float > items;//auxiliary vector

        for(int j=0; j<nc; j++){ //loop over the numbers 
            TString token = csp2[j];//read a string

            if(token.IsFloat()){
                items.push_back(token.Atof());//if string is float number store it to vector
            }else{
                items.push_back(0);
                // cout<<"ERROR: char "<<token.Data()<<" among numbers"<<endl;
                // exit(1);
            }
        }//end of the for loop


        //Fill TVectors and Map 
        int index =  GetTVectorIndex( key, 2 );
        if(  index > -1 ){
            //fKeyWordVector[index] = key;
            //fValuesVector[index] = TVector( 1, items.size(), &items[0]);
            fValuesVector[index].ResizeTo( 1, items.size()) ;
            //fValuesVector[index] = TVector( 1, items.size(), &items[0]);
            //fValuesVector[index].SetElements( &items[0] );
            for( unsigned int ii=0;ii< items.size(); ii ++ ){
                fValuesVector[index][ii+1] = items[ii];
            }


            fValueString[index] = val;
        }else{
            fKeyWordVector.push_back( key.Data() );//put the new keyword at the end of the array


            fValuesVector.push_back( TVector( 1, items.size(), &items[0]) );//store TVector to array
            fValueString.push_back( val );
            //       MapKeyWordToTVector.insert(pair<TString, unsigned int>(entryname.Data(),fKeyWordVector.size()-1)); 
            AddToKeyTable( key, fValuesVector.size()-1 ); 
        }

    }

}

void AliJBaseCard::ReadInputLine( const char *buffer ){
    // parse a line



    TString tstr(buffer); //convert the line in the buffer to TString

    if( tstr.BeginsWith("#") ) return;//skipp comments
    tstr.ReplaceAll("\t"," ");//get rid of tabelators

    //remove comment in line
    Ssiz_t startOFcomment = tstr.First('#');
    if(startOFcomment>0){
        tstr.Remove(startOFcomment,tstr.Length() - startOFcomment);
    }

    //remove white spaces from the begining
    if(tstr.BeginsWith(" ")){
        Ssiz_t startOFkeyword = 0;
        while(1){
            TString s = tstr[startOFkeyword];
            if(s.CompareTo(" ")) break;
            startOFkeyword++;
        }
        tstr.Replace(0,startOFkeyword,"",0);
    }

    //separate inputs 
    TObjArray *lineContents = tstr.Tokenize(" ");

    if(lineContents->GetEntriesFast() < 1) return;//skipp empty lines

    //----- Read a keyword -----
    TString entryname = ((TObjString*)(lineContents->At(0)))->String(); //read a key word

    if(lineContents->GetEntriesFast() == 1){
        cout<<"WARNING: single keyword "<<entryname.Data()<<" on line"<<endl;
    }else{


        //----- Read parameters -----
        vector< float > items;//auxiliary vector

        for(int i=1; i<lineContents->GetEntriesFast(); i++){ //loop over the numbers 
            TString token = ((TObjString*)(lineContents->At(i)))->String();//read a string

            if(token.IsFloat()){
                items.push_back(token.Atof());//if string is float number store it to vector
            }else{
                items.push_back(0);
                // cout<<"ERROR: char "<<token.Data()<<" among numbers"<<endl;
                // exit(1);
            }
        }//end of the for loop


        //Fill TVectors and Map 
        fKeyWordVector.push_back( entryname.Data() );//put the new keyword at the end of the array


        fValuesVector.push_back( TVector( 1, items.size(), &items[0]) );//store TVector to array
        fValueString.push_back( ((TObjString*)(lineContents->At(1)))->String() );

        //       MapKeyWordToTVector.insert(pair<TString, unsigned int>(entryname.Data(),fKeyWordVector.size()-1)); 
        AddToKeyTable( entryname, fValuesVector.size()-1 ); 


    }//else

    lineContents->~TObjArray();//remove array from heap
}

void AliJBaseCard::PrintOut(){
    // echo
    cout<<endl<<"======== "<<fcardname<<" ========="<<endl;
    for(unsigned int i=0; i<fValuesVector.size();i++){
        cout<<Form("%15s",fKeyWordVector[i].Data());//print keyword
        cout<<" (dim ="<<fValuesVector[i].GetNrows()<<") ";//print size of TVector
        for(int j=1; j<=fValuesVector[i].GetNrows(); j++){
            cout<<fValuesVector[i][j]<<" ";//TVector components
        }
        cout<<endl;
    }
}


void AliJBaseCard::WriteCard(TDirectory *file){
    // write
    cout<<endl<<"====== Writing into file ========="<<endl;

    if(!file->GetDirectory("JCard")) {
        file->mkdir("JCard");//directory to store input parameters
    }
    file->cd("JCard");
    for(unsigned int i=0;i<fValuesVector.size();i++){ 
        fValuesVector[i].Write(fKeyWordVector[i]);
    }
}




















