#ifndef ALIPHOSRAWDIGIPRODUCER_H
#define ALIPHOSRAWDIGIPRODUCER_H

class AliPHOSRawDecoder;

class AliPHOSRawDigiProducer {

public:

  AliPHOSRawDigiProducer() {}
  virtual ~AliPHOSRawDigiProducer() {}

  virtual void MakeDigits(TClonesArray *digits, AliPHOSRawDecoder* decoder);

  ClassDef(AliPHOSRawDigiProducer,1)
};

#endif
