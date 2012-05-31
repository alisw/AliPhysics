
#ifndef ALIHLTFMDRECONSTRUCTIONCOMPONENT_H
#define ALIHLTFMDRECONSTRUCTIONCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliFMDReconstructor.h"
#include "AliFMDDigit.h"
#include "AliESDEvent.h"

class AliHLTFMDReconstructionComponent : public AliHLTProcessor
{
    public:
      AliHLTFMDReconstructionComponent();
	virtual ~AliHLTFMDReconstructionComponent();
	
	const char* GetComponentID()  { return "FMDReconstruction"; }
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();
	Int_t GetRunNumber() {return fRunNumber;}
    protected:
	
	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	
	using AliHLTProcessor::DoEvent;
		
	class AliHLTFMDReconstructor : public AliFMDReconstructor {
	public:
	  AliHLTFMDReconstructor();
	  ~AliHLTFMDReconstructor() {};
	  AliESDFMD* GetFMDData() { return fESDObj; }
	  void ReconstructDigits(TClonesArray* digitArray) {ProcessDigits(digitArray);}
	private:
	  
	};
		
	
    private:
	
	unsigned long int fRunNumber;
	ClassDef(AliHLTFMDReconstructionComponent, 0)

    };
#endif
