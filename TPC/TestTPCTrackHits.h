void ConvertHits(const char * benchmark="0", Bool_t debug=kFALSE);
void CompareHits(const char * benchmark="1", Bool_t debug=kFALSE);

class AliTPChitD : public AliTPChit {
public:
  AliTPChit * GetDelta() {return &fDelta;}
private:
  AliTPChit fDelta;     //delta of hit information 
  ClassDef(AliTPChitD,1) 
};
