#ifndef AliACORDEDataModule_H
#define AliACORDEDataModule_H

#include "TObject.h"
#include "TNamed.h"


                                                                                                                                       
class AliACORDEDataModule : public TNamed {  //Set and get the name of the module
                                                                                
public:
                                                                                
       AliACORDEDataModule();
       AliACORDEDataModule(Float_t value,Bool_t fStatus,const char* name);
       //AliACORDEDataModule(const AliACORDEDataModule & Data);
       //AliACORDEDataModule& operator=(const AliACORDEDataModule & Data);
      ~AliACORDEDataModule();
       Float_t GetRate(); //Get the rate value
       Bool_t GetStatus();//Get the Module status 
       void SetRate(Float_t value); //Set the rate for modules
       void SetStatus(Bool_t status){fStatus=status;} // give the status 0 or 1 for modules
       
                                                
                                                                                
                                                                                
private:
                                                                                
      Float_t fRate; //Module Rate 
      Bool_t fStatus; //Module Status


ClassDef(AliACORDEDataModule, 2);
};
#endif
