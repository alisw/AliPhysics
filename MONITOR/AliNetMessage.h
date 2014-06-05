#ifndef AliNetMessage_H
#define AliNetMessage_H

// adapted from AliHLTMessage

#include <TBufferFile.h>

#ifndef ROOT_MessageTypes
#include <MessageTypes.h>
#endif
#ifndef ROOT_TBits
#include <TBits.h>
#endif

class AliNetMessage : public TBufferFile
{
public:
	AliNetMessage(UInt_t what = kMESS_ANY);
	AliNetMessage(void *buf, Int_t bufsize);
	virtual ~AliNetMessage();
	
	 void     ForceWriteInfo(TVirtualStreamerInfo *info, Bool_t force);
   void     Forward();
   void 		TagStreamerInfo(TVirtualStreamerInfo *info);
   TClass  *GetClass() const { return fClass;}
   void     IncrementLevel(TVirtualStreamerInfo* info);
   void     Reset();
   void     Reset(UInt_t what) { SetWhat(what); Reset(); }
   UInt_t   What() const { return fWhat; }
   void			SetLength() const;
   void     SetWhat(UInt_t what);

   void     EnableSchemaEvolution(Bool_t enable = kTRUE) { fEvolution = enable; }
   Bool_t   UsesSchemaEvolution() const { return fEvolution; }
 
   void     WriteObject(const TObject *obj);

   static void   EnableSchemaEvolutionForAll(Bool_t enable = kTRUE);
   static Bool_t UsesSchemaEvolutionForAll();

   const TList* GetStreamerInfos() const {return fInfos;}

private:
	AliNetMessage(const AliNetMessage &); 
	void operator=(const AliNetMessage &);
	
   UInt_t   fWhat;        //!Message type
   TClass  *fClass;       //!If message is kMESS_OBJECT pointer to object's class
   char    *fBufUncompressed; //!Uncompressed buffer
   TList   *fInfos;       //Array of TStreamerInfo used in WriteObject
   Bool_t   fEvolution;   //True if support for schema evolution required

   static Bool_t fgEvolution;  //True if global support for schema evolution required
  
	ClassDef(AliNetMessage, 1);
};
#endif
