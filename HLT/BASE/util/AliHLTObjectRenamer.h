#ifndef AliHLTObjectRenamer_h
#define AliHLTObjectRenamer_h

#include "AliHLTProcessor.h"
#include <string>
#include "TCollection.h"

class AliHLTObjectRenamer : public AliHLTProcessor
{
public:
	AliHLTObjectRenamer() : AliHLTProcessor(), fSuffix("") {}
	virtual ~AliHLTObjectRenamer() {}

	virtual const char* GetComponentID() { return "ObjectRenamer"; }
	virtual void GetInputDataTypes( vector<AliHLTComponentDataType>& list) { list.push_back(kAliHLTAnyDataType); } 
	virtual AliHLTComponentDataType GetOutputDataType() { return kAliHLTAnyDataType; }
	
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
	{
		constBase = 1024*1024;
		inputMultiplier = 1;
	}
	
	virtual AliHLTComponent* Spawn() { return new AliHLTObjectRenamer; }

protected:
	
	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit() { return 0; }

	virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

	using AliHLTProcessor::DoEvent;
	
private:
	// Do not allow copying of this class.
	AliHLTObjectRenamer(const AliHLTObjectRenamer&);
	AliHLTObjectRenamer& operator=(const AliHLTObjectRenamer&);
	
	TString fSuffix;   /// Suffix to appent to object names.

  template <typename ListType>
    TObject* CloneList(const TObject* obj, const TString& suffix)
    {
      if (obj->IsA() == ListType::Class())
      {
        const ListType* oldArray = static_cast<const ListType*>(obj);
        ListType* newArray = new ListType;
        newArray->SetOwner(kTRUE);
        TIter next(oldArray);
        const TObject* tmpObj = NULL;
        while ((tmpObj = next()) != NULL)
        {
          TString name = tmpObj->GetName();
          name += suffix;
          TObject* newObj = tmpObj->Clone(name);
          newArray->Add(newObj);
        }
        return newArray;
      }
      return NULL;
    }

  ClassDef(AliHLTObjectRenamer, 0)
};
#endif
