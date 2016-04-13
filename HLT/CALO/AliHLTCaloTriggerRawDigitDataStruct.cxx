#include <AliHLTCaloTriggerRawDigitDataStruct.h>
#include <cstring>


void InitializeRawDigit(AliHLTCaloTriggerRawDigitDataStruct &rawdigit){
  rawdigit.fID = 0;
  rawdigit.fTriggerBits = 0;
  rawdigit.fL1TimeSum = 0;
  rawdigit.fNL0Times = 0;
  memset(rawdigit.fL0Times, 0, sizeof(rawdigit.fL0Times));
  rawdigit.fNTimeSamples = 0;
  memset(rawdigit.fTimeSamples, 0, sizeof(rawdigit.fTimeSamples));
}

void SetRawDigitID(AliHLTCaloTriggerRawDigitDataStruct &digit, Int_t id){
  digit.fID = id;
}

void SetTriggerBit(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t bit, Int_t mode){
  dig.fTriggerBits |= 1 << (bit + mode * kTriggerTypeEnd);
}

void SetL0Time(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t i) {
  bool hasfound(false);
  for (Int_t j = 0; j < dig.fNL0Times; j++) {
    if (i == dig.fL0Times[j]) {
      hasfound = true;
      break;
    }
  }
  if(hasfound) return;
  dig.fNL0Times++;
  if (dig.fNL0Times > 9) return;
  dig.fL0Times[dig.fNL0Times - 1] = i;
}

void SetL1TimeSum(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t l1timeSum){
  dig.fL1TimeSum = l1timeSum;
}

void SetL1SubRegion(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t l1subregion){
  dig.fL1SubRegion = l1subregion;
}

void SetTimeSamples(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t nsamples, Int_t *samples){
  memcpy(dig.fTimeSamples, samples, sizeof(Int_t) * 15);
  dig.fNTimeSamples = nsamples;
}

Int_t GetRawDigitID(const AliHLTCaloTriggerRawDigitDataStruct &dig) {
  return dig.fID;
}

Bool_t GetL0Time(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t i, Int_t& time) {
  if (i < 0 || i > dig.fNL0Times) {
    return kFALSE;
  }
  time = dig.fL0Times[i];
  return kTRUE;
}

void GetL0Times(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t times[]) {
  for (Int_t i = 0; i < dig.fNL0Times; i++) times[i] = dig.fL0Times[i];
}

Bool_t GetTimeSample(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t iSample, Int_t& timeBin, Int_t& amp) {
  if (iSample > dig.fNTimeSamples || iSample < 0) return kFALSE;
  amp     = (Short_t)(dig.fTimeSamples[iSample] & 0xFFFF);
  timeBin = (Short_t)(dig.fTimeSamples[iSample] >> 16 );
  return kTRUE;
}


Int_t GetL0TimeSum(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t time) {
  Int_t value = 0;
  for (Int_t i = 0; i < dig.fNTimeSamples; i++) {
    Int_t timeBin, amp;
    GetTimeSample(dig, i, timeBin, amp);

    if (timeBin >= time && timeBin < time + 4) value += amp;
  }
  return value;
}

Int_t GetL1SubRegion(const AliHLTCaloTriggerRawDigitDataStruct &dig){
  return dig.fL1SubRegion;
}

Int_t GetTriggerBit(const AliHLTCaloTriggerRawDigitDataStruct &dig, const TriggerType_t type, const Int_t mode) {
  Int_t shift = kTriggerTypeEnd * mode;
  Int_t mask  = 1 << type;
  return ((dig.fTriggerBits >> shift) & mask);
}

Bool_t GetRawDigitMaximumAmplitude(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t& amplitude, Int_t& time) {
  if (!dig.fNTimeSamples) return kFALSE;

  amplitude = 0;
  for (Int_t i = 0; i < dig.fNTimeSamples; i++) {
    Int_t t, a;
    if (GetTimeSample(dig, i, t, a)) {
      if (a > amplitude)
      {
        amplitude = a;
        time      = t;
      }
    }
  }
  return kTRUE;
}

void PrintRawDigit(const AliHLTCaloTriggerRawDigitDataStruct &dig) {
  printf("===\n| Digit  id: %4d /  %d Time Samples: \n", dig.fID, dig.fNTimeSamples);
  for (Int_t i=0; i < dig.fNTimeSamples; i++) {
    Int_t timeBin, amp;
    GetTimeSample(dig, i, timeBin, amp);
    printf("| (%d,%d) ",timeBin,amp);
  }
  printf("\n");
  printf("| L0: (%d,%d) / %d Time(s): \n",GetTriggerBit(dig, kL0,1),GetTriggerBit(dig, kL0,0), dig.fNL0Times);
  for (Int_t i = 0; i < dig.fNL0Times; i++) {
    Int_t time;
    if (GetL0Time(dig, i, time)) printf("| %d ",time);
  }
  printf("\n");
  printf("| L1: g high (%d,%d) g low (%d,%d) j high (%d,%d) j low (%d,%d) / Time sum: %d\n",
       GetTriggerBit(dig, kL1GammaHigh,1),GetTriggerBit(dig, kL1GammaHigh,0),GetTriggerBit(dig, kL1GammaLow,1),GetTriggerBit(dig, kL1GammaLow,0),
       GetTriggerBit(dig, kL1JetHigh,1),  GetTriggerBit(dig, kL1JetHigh,0),  GetTriggerBit(dig, kL1JetLow,1),  GetTriggerBit(dig, kL1JetLow,0),
       dig.fL1TimeSum);
}

void SPrintRawDigit(const AliHLTCaloTriggerRawDigitDataStruct &dig, char *outputbuffer) {
  char *currentpos = outputbuffer, tempstr[1024];
  sprintf(tempstr, "===\n| Digit  id: %4d /  %d Time Samples: \n", dig.fID, dig.fNTimeSamples);
  strcpy(currentpos, tempstr);
  currentpos += strlen(tempstr)/sizeof(char);
  for (Int_t i=0; i < dig.fNTimeSamples; i++) {
    Int_t timeBin, amp;
    GetTimeSample(dig, i, timeBin, amp);
    sprintf(tempstr, "| (%d,%d) ",timeBin,amp);
    strcpy(currentpos, tempstr);
    currentpos += strlen(tempstr)/sizeof(char);
  }
  sprintf(tempstr, "| L0: (%d,%d) / %d Time(s): \n",GetTriggerBit(dig, kL0,1),GetTriggerBit(dig, kL0,0), dig.fNL0Times);
  strcpy(currentpos, tempstr);
  currentpos += strlen(tempstr)/sizeof(char);
  for (Int_t i = 0; i < dig.fNL0Times; i++) {
    Int_t time;
    if (GetL0Time(dig, i, time)) sprintf(tempstr, "| %d ",time);
    strcpy(currentpos, tempstr);
    currentpos += strlen(tempstr)/sizeof(char);
  }
  sprintf(currentpos, "| L1: g high (%d,%d) g low (%d,%d) j high (%d,%d) j low (%d,%d) / Time sum: %d\n",
       GetTriggerBit(dig, kL1GammaHigh,1),GetTriggerBit(dig, kL1GammaHigh,0),GetTriggerBit(dig, kL1GammaLow,1),GetTriggerBit(dig, kL1GammaLow,0),
       GetTriggerBit(dig, kL1JetHigh,1),  GetTriggerBit(dig, kL1JetHigh,0),  GetTriggerBit(dig, kL1JetLow,1),  GetTriggerBit(dig, kL1JetLow,0),
       dig.fL1TimeSum);
}
