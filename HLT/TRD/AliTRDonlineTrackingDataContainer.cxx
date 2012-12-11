#include <cstdlib>
#include "TString.h"
#include "AliHLTLogging.h"
#include "AliTRDonlineTrackingDataContainer.h"


ClassImp(AliTRDonlineTrackingDataContainer)

const Float_t AliTRDonlineTrackingDataContainer::fgkBinWidthY   = 160e-4; // 160 um

AliTRDonlineTrackingDataContainer::AliTRDonlineTrackingDataContainer() : AliHLTLogging(),
  fGtuPtMultiplier(1.),
  fStoreGtuInfo(kTRUE),
  fLogPrefix(NULL)
{
  fLogPrefix = new TString("");
  Clear();
}

AliTRDonlineTrackingDataContainer::AliTRDonlineTrackingDataContainer(const AliTRDonlineTrackingDataContainer& /*cont*/) : AliHLTLogging(),
  fGtuPtMultiplier(1.),
  fStoreGtuInfo(kTRUE),
  fLogPrefix(NULL)
{
  fLogPrefix = new TString("");
  Clear();
}

AliTRDonlineTrackingDataContainer::~AliTRDonlineTrackingDataContainer(){
  if (fLogPrefix)
    delete fLogPrefix;
  fLogPrefix = NULL;
}

void AliTRDonlineTrackingDataContainer::Clear(const Option_t*) {
  memset(fNumTracklets, 0, sizeof(UInt_t)*fgkNumChambers);
  memset(fNumTracks, 0, sizeof(UInt_t)*fgkNumStacks);

  memset(fSectorTrgWords, 0, sizeof(UInt_t)*fgkNumSectors);
  memset(fStackTrgWords, 0, sizeof(ULong64_t)*fgkNumSectors*fgkNumStacksPerSector);

  //## todo: only for debugging, may be removed without harm
//  memset(fTrackletWords, 0, sizeof(UInt_t)*fgkNumChambers*fgkMaxTrackletsPerChamber);
//  memset(fTrackletHCId, 0, sizeof(UInt_t)*fgkNumChambers*fgkMaxTrackletsPerChamber);
//  memset(fTrackWords, 0, sizeof(ULong64_t)*fgkNumStacks*fgkMaxTracksPerStack);
//  memset(fTrackExtWords, 0, sizeof(ULong64_t)*fgkNumStacks*fgkMaxTracksPerStack);
//  memset(fTrackTrackletWords, 0, sizeof(UInt_t)*fgkNumStacks*fgkMaxTracksPerStack*fgkNumLayers);
}

Int_t AliTRDonlineTrackingDataContainer::GetNumTracklets() {
  Int_t count = 0;
  for (UShort_t det = 0; det < fgkNumChambers; ++det)
    count += fNumTracklets[det];
  return count;
}

Int_t AliTRDonlineTrackingDataContainer::GetNumTracklets(UInt_t det) {
  return fNumTracklets[det];
}

Int_t AliTRDonlineTrackingDataContainer::GetNumTracks() {
  Int_t count = 0;
  for (UShort_t stack = 0; stack < fgkNumStacks; ++stack)
    count += fNumTracks[stack];
  return count;
}

Int_t AliTRDonlineTrackingDataContainer::GetNumTracks(UShort_t stack){
  return fNumTracks[stack];
}

Int_t AliTRDonlineTrackingDataContainer::GetTrackletBinY(UInt_t det, UInt_t trackletIndex) {
  UInt_t trackletWord = fTrackletWords[det][trackletIndex];
  if (trackletWord & 0x1000) {
    return -((~(trackletWord - 1)) & 0x1fff);
  }
  else {
    return (trackletWord & 0x1fff);
  }
}

Int_t AliTRDonlineTrackingDataContainer::GetTrackletBinDy(UInt_t det, UInt_t trackletIndex) {
  UInt_t trackletWord = fTrackletWords[det][trackletIndex];
  if (trackletWord & (1 << 19))
    return -((~((trackletWord >> 13) - 1)) & 0x7f);
  else
    return ((trackletWord >> 13) & 0x7f);
};

Int_t AliTRDonlineTrackingDataContainer::GetTrackletBinZ(UInt_t det, UInt_t trackletIndex) {
  return ((fTrackletWords[det][trackletIndex] >> 20) & 0xf);
}

Int_t AliTRDonlineTrackingDataContainer::GetTrackletPID(UInt_t det, UInt_t trackletIndex) {
  return ((fTrackletWords[det][trackletIndex] >> 24) & 0xff);
};

Float_t AliTRDonlineTrackingDataContainer::GetTrackletLocalY(UInt_t det, UInt_t trackletIndex) {
  return GetTrackletBinY(det, trackletIndex) * fgkBinWidthY;
}

AliESDTrdTracklet* AliTRDonlineTrackingDataContainer::GetTracklet(UInt_t det, UInt_t trackletIndex) {
  AliESDTrdTracklet* trkl = NULL;
  if ((det < fgkNumChambers) && (trackletIndex < fNumTracklets[det])){
    trkl = new AliESDTrdTracklet(fTrackletWords[det][trackletIndex], fTrackletHCId[det][trackletIndex], -1);
  }
  return trkl;
}

AliESDTrdTrack* AliTRDonlineTrackingDataContainer::GetTrack(UInt_t stack, UInt_t trackIndex, Bool_t constructTracklets){
  AliESDTrdTrack* trk = NULL;
  if ((stack < fgkNumStacks) && (trackIndex < fNumTracks[stack])){
    trk = new AliESDTrdTrack();
    ULong64_t tw = fTrackWords[stack][trackIndex];
    ULong64_t etw = fTrackWords[stack][trackIndex];
    trk->SetLayerMask(GetTrackLayerMask(stack, trackIndex));
    trk->SetA(GetTrackA(stack, trackIndex));
    trk->SetB( (((tw >> 20) & 0x3ffff) ^ 0x20000) - 0x20000);
    trk->SetC( (((tw >> 8)  &  0xffff) ^  0x8000) -  0x8000);
    trk->SetPID(GetTrackPID(stack, trackIndex));
    trk->SetSector(stack/5);
    trk->SetStack(stack%5);
    trk->SetLabel(-3);
    trk->SetFlags((etw >> 52) & 0x7ff);
    trk->SetReserved((etw >> 49) & 0x7);
    trk->SetY((etw >> 36) & 0x1fff);
    trk->SetTrackletIndex((etw >>  0) & 0x3f, 0);
    trk->SetTrackletIndex((etw >>  6) & 0x3f, 1);
    trk->SetTrackletIndex((etw >> 12) & 0x3f, 2);
    trk->SetTrackletIndex((etw >> 18) & 0x3f, 3);
    trk->SetTrackletIndex((etw >> 24) & 0x3f, 4);
    trk->SetTrackletIndex((etw >> 30) & 0x3f, 5);

    if (constructTracklets) {
      for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer){
	AliESDTrdTracklet * trkl = new AliESDTrdTracklet(GetTrackTrackletWord(stack, trackIndex, iLayer), 2*(stack*6 + iLayer));
	trk->AddTrackletReference(trkl, iLayer);
      }
    }

  } else {
    HLTError("invalid stack (%d) or track index (%d) in GetTrack", stack, trackIndex);
  }
  return trk;
}

Double_t AliTRDonlineTrackingDataContainer::GetTrackPt(UInt_t stack, UInt_t trackIndex){

  // calculate pt from a as done in hardware
  const Int_t maskIdLut[64] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,
    -1, -1, -1, -1, -1, -1, -1,  1, -1, -1, -1,  2, -1,  3,  4,  5,
    -1, -1, -1, -1, -1, -1, -1,  6, -1, -1, -1,  7, -1,  8,  9, 10,
    -1, -1, -1, 11, -1, 12, 13, 14, -1, 15, 16, 17, 18, 19, 20, 21
  };

  const Int_t c1Lut[32] = {
    -2371, -2474, -2474, -2474, -2563, -2448, -2578, -2578,
    -2578, -2670, -2557, -2578, -2578, -2670, -2557, -2578,
    -2670, -2557, -2763, -2557, -2644, -2523,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1
  };

  Int_t a = GetTrackA(stack, trackIndex);
  UShort_t lm = GetTrackLayerMask(stack, trackIndex);

  if (a != 0) {
    Int_t layerMaskId = maskIdLut[lm];
    Int_t c1 = c1Lut[layerMaskId];
    Int_t c1Ext = c1 << 8;
    Int_t ptRawStage4 = c1Ext / ((a >> 2) != 0 ? (a >> 2) : 1 );
    Int_t ptRawComb4 = ptRawStage4;
    Int_t ptExtComb4 = (ptRawComb4 > 0) ? ptRawComb4 + 33 : ptRawComb4 - 30;

    return ((-ptExtComb4/2) / 128.) * fGtuPtMultiplier;
  }
  else
    return 0.;

}

UShort_t AliTRDonlineTrackingDataContainer::GetTrackLayerNum(UInt_t stack, UInt_t trackIndex) {
  UShort_t num = 0;
  UShort_t layerMask = GetTrackLayerMask(stack, trackIndex);
  for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer)
    if ((layerMask >> iLayer) & 1)
      num++;
  return num;
}

UInt_t AliTRDonlineTrackingDataContainer::GetTrackTrackletWord(UInt_t stack, UInt_t trackIndex, UShort_t layer) {
  return fTrackTrackletWords[stack][trackIndex][layer];
}


Int_t AliTRDonlineTrackingDataContainer::GetTrackTrackletBinY(UInt_t stack, UInt_t trackIndex, UShort_t layer) {
  UInt_t trackletWord = fTrackTrackletWords[stack][trackIndex][layer];
  if (trackletWord & 0x1000) {
    return -((~(trackletWord - 1)) & 0x1fff);
  }
  else {
    return (trackletWord & 0x1fff);
  }
}

Int_t AliTRDonlineTrackingDataContainer::GetTrackAddInfo(UShort_t stack, UInt_t trackIndex) {
  return fTrackAddInfo[stack][trackIndex];
}

Float_t AliTRDonlineTrackingDataContainer::GetTrackTrackletLocalY(UInt_t stack, UInt_t trackIndex, UShort_t layer) {
  return GetTrackTrackletBinY(stack, trackIndex, layer) * fgkBinWidthY;
}

Int_t AliTRDonlineTrackingDataContainer::GetTrackTrackletBinZ(UInt_t stack, UInt_t trackIndex, UShort_t layer) {
  return ((fTrackTrackletWords[stack][trackIndex][layer] >> 20) & 0xf);
}

UShort_t AliTRDonlineTrackingDataContainer::GetTrackTrackletPID(UInt_t stack, UInt_t trackIndex, UShort_t layer) {
  return ((fTrackTrackletWords[stack][trackIndex][layer] >> 24) & 0xff);
}

Int_t AliTRDonlineTrackingDataContainer::AddTracklet(UInt_t HCId, UInt_t trackletWord) {
  UShort_t det = HCId/2;
  Int_t pos = fNumTracklets[det]++;
  fTrackletWords[det][pos] = trackletWord;
  fTrackletHCId[det][pos] = HCId;
  return pos;
}

Int_t AliTRDonlineTrackingDataContainer::AddTracklet(const AliESDTrdTracklet* tracklet) {
  return AddTracklet(tracklet->GetHCId(), tracklet->GetTrackletWord());
}

Int_t AliTRDonlineTrackingDataContainer::AddTrack(UShort_t stack,
						  ULong64_t trackWord, ULong64_t extTrackWord,
						  const UInt_t trackletWords[6], const Int_t addInfo){
  Int_t pos = fNumTracks[stack]++;
  fTrackWords[stack][pos] = trackWord;
  fTrackExtWords[stack][pos] = extTrackWord;
  fTrackAddInfo[stack][pos] = addInfo;
  for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer){
    fTrackTrackletWords[stack][pos][iLayer] = trackletWords[iLayer];
  }
  return pos;
}

Int_t AliTRDonlineTrackingDataContainer::AddTrack(const AliESDTrdTrack* track, Int_t addInfo) {

  UInt_t trackletWords[fgkNumLayers];
  UShort_t lm = track->GetLayerMask();
  Bool_t checkOk = kTRUE;
  for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer){
    if ((lm >> iLayer) & 1){
      if (track->GetTracklet(iLayer))
	trackletWords[iLayer] = track->GetTracklet(iLayer)->GetTrackletWord();
      else
	checkOk = kFALSE;  // inconsistency between layer mask and tracklet pointers: do not use this track
    } else
      trackletWords[iLayer] = fgkInvalidTrackletWord;
  }

  if (checkOk){
    return AddTrack(track->GetSector()*5 + track->GetStack(),
		    track->GetTrackWord(0), track->GetExtendedTrackWord(0), trackletWords, addInfo);
  } else {
    // DbgLog("ERROR", "Ignoring GTU track with inconsistency between layer mask and tracklet pointers");
    printf("<ERROR> Ignoring GTU track with inconsistency between layer mask and tracklet pointers\n");
    return -1;
  }

}

void AliTRDonlineTrackingDataContainer::SetTrackAddInfo(UShort_t stack, UInt_t trackIndex, Int_t addInfo) {
  fTrackAddInfo[stack][trackIndex] = addInfo;
}

void AliTRDonlineTrackingDataContainer::SetGtuPtMultiplierFromMagField(Double_t magField) {
  if (magField > 0)
    fGtuPtMultiplier = -1.;
  else
    fGtuPtMultiplier = 1.;
}

void AliTRDonlineTrackingDataContainer::PrintBuffer(const void* buffer, UInt_t sizeInBytes, const char* identifier) {

  UInt_t* buf = (UInt_t*)buffer;
  UInt_t currPos = 0;
  TString str("");
  TString ident(identifier);

  HLTDebug("BUFFER DUMP for <%s>", ident.Data());

  while (currPos < sizeInBytes/4){

    str = Form("%06d: 0x%08x  ", currPos, buf[currPos]);

    if (currPos == 0){  // leading magic constant

      if (buf[currPos++] == fgkLeadingMagicConst)
	str += "correct leading magic constant";
      else
	str += Form("incorrect leading magic constant, should be 0x%08x",
		    fgkLeadingMagicConst);
    }

    if (currPos == sizeInBytes/4 - 1){

      if (buf[sizeInBytes/4 - 1] == fgkTrailingMagicConst)
	str += "correct trailing magic constant";
      else
	str += Form("incorrect trailing magic constant, should be 0x%08x",
		    fgkTrailingMagicConst);
    }

    currPos++;

    HLTDebug(str.Data());

  }

}

void AliTRDonlineTrackingDataContainer::PrintSummary(const char* identifier){

  TString ident(identifier);

  HLTDebug("TRDGM AliTRDonlineTrackingDataContainer <%s> summary: %5d tracklets, %2d tracks, %6ld bytes comp. mem size",
	   ident.Data(), GetNumTracklets(), GetNumTracks(), DataWordsNeeded()*sizeof(UInt_t));

}

inline UInt_t AliTRDonlineTrackingDataContainer::MakeDataHeader(){
  if (fStoreGtuInfo)
    return (fgkHeaderTypeData << 28) | (1 << 16) | DataWordsNeeded();
  else
    return (fgkHeaderTypeData << 28) | DataWordsNeeded();
}

inline UInt_t AliTRDonlineTrackingDataContainer::MakeStackHeader(UShort_t stack,
								 Bool_t trackletsPresent,
								 Bool_t tracksPresent,
								 UInt_t size){
  return
    (fgkHeaderTypeStack << 28) |
    (((tracksPresent) ? 1 : 0) << 27) |
    (((trackletsPresent) ? 1 : 0) << 26) |
    ((stack & 0xff) << 16) | (size & 0xffff);
}

inline Int_t AliTRDonlineTrackingDataContainer::StackFromStackHeader(UInt_t stackHeader) {
  return (stackHeader >> 16) & 0xff;
}

inline UInt_t AliTRDonlineTrackingDataContainer::SizeFromStackHeader(UInt_t stackHeader) {
  return stackHeader & 0xffff;
}

inline UInt_t AliTRDonlineTrackingDataContainer::MakeTrackletHeader(UShort_t halfChamber){
  return (fgkHeaderTypeTracklet << 28) | ((halfChamber & 0xfff) << 16) | 0x2;
}

inline Int_t AliTRDonlineTrackingDataContainer::HCIdFromTrackletHeader(UInt_t trackletHeader) {
  return (trackletHeader >> 16) & 0xfff;
}

inline UInt_t AliTRDonlineTrackingDataContainer::MakeTrackHeader(UShort_t stack){
  return (fgkHeaderTypeTrack << 28) | ((stack & 0xff) << 16) | 12;
}

inline UInt_t AliTRDonlineTrackingDataContainer::MakeGtuHeader(Bool_t storeInfo){
  if (storeInfo)
    return (fgkHeaderTypeGtu << 28) | (1 << 16) | (1 + 18 + 2*90);
  else
    return (fgkHeaderTypeGtu << 28) | 1;
}

inline Bool_t AliTRDonlineTrackingDataContainer::StoreGtuInfoFromDataHeader(UInt_t dataHeader){
  return ((dataHeader >> 16) & 1);
}

inline Bool_t AliTRDonlineTrackingDataContainer::GtuInfoPresentFromGtuHeader(UInt_t gtuHeader){
  return ((gtuHeader >> 16) & 1);
}

inline UInt_t AliTRDonlineTrackingDataContainer::DataWordsNeeded() {

  UInt_t size = 0;

  size += 1;                              // leading magic word
  size += 1;                              // overall data header
  size += 90;                             // stack headers
  size += GetNumTracklets()*(1 + 1);      // tracklets (mark + trackword)
  size += GetNumTracks()*(1 + 4 + 1 + 6); // GTU tracks (mark + trackword[2] + exttrackword[2] + addInfo + 6*trackletword)
  size += 1;                              // crosscheck word
  size += 1;                              // trailing magic word

  if (fStoreGtuInfo)
    size += 1 + 18 + 2*90;                // trigger header + 18 sector trigger words, 90 stack tracking done words[2]
  else
    size += 1;                            // trigger header only

  return size;
}

inline UInt_t calcCheck(UInt_t prevCrc, UInt_t c) {
  // CCITT 16 bit (X^16 + X^12 + X^5 + 1).
  UInt_t crc = prevCrc;
  crc  = (unsigned char)(crc >> 8) | (crc << 8);
  crc ^= c;
  crc ^= (unsigned char)(crc & 0xff) >> 4;
  crc ^= (crc << 8) << 4;
  crc ^= ((crc & 0xff) << 4) << 1;
  return (crc);
}

Bool_t AliTRDonlineTrackingDataContainer::Compress(void* &buffer, UInt_t &sizeInBytes){

  // TString ptrInfo("COMPRESS CALLED: [");

  UInt_t bufSize = sizeof(UInt_t)*DataWordsNeeded();
  UInt_t* buf = (UInt_t*)malloc(bufSize);
  UInt_t crosscheck = fgkCrosscheckSeed;
  UInt_t lastStackHeaderPos;

  if (!buf){
    HLTError("Could not allocate %d bytes for buffer", bufSize);
    return kFALSE;
  }

  memset(buf, 0, bufSize);
  // ptrInfo += Form(" memset(%p, 0, %d) -> %p - %p, ", buf, bufSize, buf, (char*)buf + bufSize);
  UInt_t currPos = 0;
  UInt_t det;

  // ptrInfo += Form(" lowest write addr: %p, ", &(buf[currPos]));
  buf[currPos++] = fgkLeadingMagicConst;
  buf[currPos++] = MakeDataHeader();
  for (UShort_t iStack = 0; iStack < 90; ++iStack){

    // add stack info
    lastStackHeaderPos = currPos;
    buf[currPos++] = MakeStackHeader(iStack, 1, 1, 0);

    // add tracklet infos
    for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer){
      det = iStack*fgkNumLayers + iLayer;
      for (UInt_t iTrkl = 0; iTrkl < fNumTracklets[det]; ++iTrkl){
	buf[currPos++] = MakeTrackletHeader(fTrackletHCId[det][iTrkl]);
	buf[currPos++] = fTrackletWords[det][iTrkl];
      }
    }

    // add track infos
    for (UInt_t iTrk = 0; iTrk < fNumTracks[iStack]; ++iTrk){
      buf[currPos++] = MakeTrackHeader(iStack);
      buf[currPos++] = fTrackWords[iStack][iTrk] & 0xffffffff;
      buf[currPos++] = (fTrackWords[iStack][iTrk] >> 32) & 0xffffffff;
      buf[currPos++] = fTrackExtWords[iStack][iTrk] & 0xffffffff;
      buf[currPos++] = (fTrackExtWords[iStack][iTrk] >> 32) & 0xffffffff;
      buf[currPos++] = fTrackAddInfo[iStack][iTrk];
      for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer)
	buf[currPos++] = fTrackTrackletWords[iStack][iTrk][iLayer];
    }

    buf[lastStackHeaderPos] = MakeStackHeader(iStack,
					      ((GetNumTracklets() > 0) ? 1 : 0),
					      ((GetNumTracks() > 0) ? 1 : 0),
					      currPos - lastStackHeaderPos); // update size  //## todo: check off-by-one issues

  } // loop over all stacks

  // add Gtu header info
  buf[currPos++] = MakeGtuHeader(fStoreGtuInfo);
  if (fStoreGtuInfo){
    // store trigger info from GTU headers
    for (UShort_t iSector = 0; iSector < fgkNumSectors; ++iSector)
      buf[currPos++] = fSectorTrgWords[iSector];
    for (UShort_t iSector = 0; iSector < fgkNumSectors; ++iSector){
      for (UShort_t iStack = 0; iStack < fgkNumStacksPerSector; ++iStack){
	buf[currPos++] = fStackTrgWords[iSector][iStack] & 0xffffffff;
	buf[currPos++] = (fStackTrgWords[iSector][iStack] >> 32) & 0xffffffff;
      }
    }
  }

  // calc crosscheck
  for (UInt_t ofs = 1; ofs < currPos; ++ofs)
    crosscheck = calcCheck(crosscheck, buf[ofs]);

  buf[currPos++] = crosscheck;
  buf[currPos++] = fgkTrailingMagicConst;
  // ptrInfo += Form(" highest write addr: %p, %d]", &(buf[currPos - 1]), (currPos-1)*4);

  if (sizeof(UInt_t)*currPos != bufSize){
    HLTError("inconsistent memory layout! (%ld %d)", sizeof(UInt_t)*currPos, bufSize);
  }

//  for (UInt_t ofs = 0; ofs < bufSize/4; ++ofs)
//    printf("%04d: 0x%08x\n", ofs, buf[ofs]);

  buffer = buf;
  sizeInBytes = bufSize;

  // DbgLog("", ptrInfo.Data());

  return kFALSE;
}

Bool_t AliTRDonlineTrackingDataContainer::Decompress(const void* buffer, UInt_t sizeInBytes, Bool_t cumulative, Bool_t verbose) {

  // TString ptrInfo(Form("DECOMPRESS CALLED: [buf: %p, size: %d - %p - %p, ", buffer, sizeInBytes, buffer, (char*)buffer + sizeInBytes));

  UInt_t* buf = (UInt_t*)buffer;
  UInt_t currPos = 0;
  UInt_t crosscheck = fgkCrosscheckSeed;
  UInt_t stackHeader;
  UInt_t gtuInfoHeader = 0;
  Int_t stack = 0;
  UInt_t size = 0;

  if (!cumulative)
    Clear();

  if (buf[currPos++] != fgkLeadingMagicConst){
    HLTError("invalid data: leading mark should be 0x%08x, is 0x%08x", fgkLeadingMagicConst, buf[currPos - 1]);
    return kFALSE;
  } else if (verbose)
    HLTError("0x%05d: 0x%08x  correct leading mark", currPos, buf[currPos - 1]);

  UInt_t overallDataHeader = buf[currPos++];
  if ((overallDataHeader >> 28) != fgkHeaderTypeData){
    HLTError("invalid data header: 0x%08x", overallDataHeader);
    return kFALSE;
  } else {
    fStoreGtuInfo = StoreGtuInfoFromDataHeader(overallDataHeader);
  }

  if (buf[sizeInBytes/4 - 1] != fgkTrailingMagicConst){
    HLTError("invalid data: trailing mark should be 0x%08x, is 0x%08x", fgkTrailingMagicConst, buf[sizeInBytes/4 - 1]);
    return kFALSE;
  } else if (verbose){
    HLTDebug("0x%05d: 0x%08x  correct trailing mark", sizeInBytes/4 - 1, buf[sizeInBytes/4 - 1]);
  }

  while (currPos < (sizeInBytes/4) - 2) { // stack-level loop

    stackHeader = buf[currPos++];

    // stack header + encapsulated tracklet and track data
    if (((stackHeader >> 28) & 0xf) == fgkHeaderTypeStack){
      stack = StackFromStackHeader(stackHeader);
      size = SizeFromStackHeader(stackHeader);

      if (verbose){
	HLTDebug("STACK HEADER: 0x%08x - S%02d-%d, size: %d  [checksum: 0x%08x]", stackHeader, stack/5, stack%5, size, buf[sizeInBytes/4 - 2]);
      }

      while (currPos < sizeInBytes/4 - 2){

	if (((buf[currPos] >> 28) & 0xf) == fgkHeaderTypeTracklet){
	  UInt_t trklHdr = buf[currPos++];
	  UInt_t trklWord = buf[currPos++];
	  AddTracklet(HCIdFromTrackletHeader(trklHdr), trklWord);
	  if (verbose){
	    HLTDebug("Tracklet: 0x%08x 0x%08x", trklHdr, trklWord);
	  }
	} else if (((buf[currPos] >> 28) & 0xf) == fgkHeaderTypeTrack){
	  UInt_t trkHdr = buf[currPos++];
	  if (trkHdr == 0)
	    HLTError("ERROR", "invalid track header");
	  ULong64_t trackWord = buf[currPos++];
	  trackWord |= ((ULong64_t)buf[currPos++] << 32);
	  ULong64_t extTrackWord = buf[currPos++];
	  extTrackWord |= ((ULong64_t)buf[currPos++] << 32);
	  UInt_t addInfo = buf[currPos++];
	  UInt_t trackletWords[6];
	  for (UShort_t iLayer = 0; iLayer < fgkNumLayers; ++iLayer){
	    trackletWords[iLayer] = buf[currPos++];
	  }
	  AddTrack(stack, trackWord, extTrackWord, trackletWords, addInfo);
	  if (verbose){
	    HLTDebug("GTU track: 0x%016llx 0x%016llx", trackWord, extTrackWord);
	  }
	} else if (((buf[currPos] >> 28) & 0xf) == fgkHeaderTypeStack){
	//printf("next stack header\n");
	  break;
	} else if (((buf[currPos] >> 28) & 0xf) == fgkHeaderTypeGtu){
	  gtuInfoHeader = buf[currPos];
	  break;
	} else {
	  HLTError("unknown data block: 0x%08x", buf[currPos]);
	  break;
	}
      }

      if (gtuInfoHeader)
	break;

    } else {
      HLTError("invalid header while analyzing tracking data block: 0x%08x - S%02d-%d, size: %d",
			   stackHeader, stack/5, stack%5, size);
      break;
    }

  } // stack-level loop

  // GTU header data loop
  if (((gtuInfoHeader >> 28) & 0xf) == fgkHeaderTypeGtu) {
    UInt_t gtuHeader = buf[currPos++];
    HLTInfo("gtu header: 0x%08x", gtuHeader);
    if (GtuInfoPresentFromGtuHeader(gtuHeader)){
      for (UShort_t iSector = 0; iSector < fgkNumSectors; ++iSector)
	fSectorTrgWords[iSector] = buf[currPos++];
      for (UShort_t iSector = 0; iSector < fgkNumSectors; ++iSector){
	for (UShort_t iStack = 0; iStack < fgkNumStacksPerSector; ++iStack){
	  UInt_t low = buf[currPos++];
	  UInt_t high = buf[currPos++];
	  fStackTrgWords[iSector][iStack] = ((ULong64_t)high << 32) | low;
	}
      }
    }
  } else {
    HLTError("expected GtuInfoHeader at position %d, but is 0x%08x",
			 currPos, buf[currPos]);
  }

  if (currPos != (sizeInBytes/4) - 2){
    HLTError("invalid read position after analyzing tracking data block.");
    return kFALSE;
  }

  for (UInt_t ofs = 1; ofs < sizeInBytes/4 - 2; ++ofs)
    crosscheck = calcCheck(crosscheck, buf[ofs]);

  if (crosscheck != buf[sizeInBytes/4 - 2]){  // compare recalculated checksum with the one in the data
    HLTError("decompress checksum mismatch: should be 0x%08x, is 0x%08x",
	   crosscheck, buf[sizeInBytes/4 - 2]);
    return kFALSE;
  }

  // HLTDebug(ptrInfo.Data());

  return kTRUE;

}
