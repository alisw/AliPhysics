#ifndef ALIJESE_H
#define ALIJESE_H

class AliJESE{
public:
	AliJESE();
	~AliJESE();
	typedef unsigned int uint;
	bool Initialize();
	void Destroy();
	double Getqc2Perc(class AliAODEvent *, float, uint);
	enum OBJECT{
		OBJECT_MULTV0,
		OBJECT_QXA2M,
		OBJECT_QXA2S,
		OBJECT_QYA2M,
		OBJECT_QYA2S,
		OBJECT_QXC2M,
		OBJECT_QXC2S,
		OBJECT_QYC2M,
		OBJECT_QYC2S,
		OBJECT_COUNT
	};
	class TFile *poadbf, *psplf;
	class AliOADBContainer *poadb[OBJECT_COUNT];
	class TSpline3 *psplineQ2c[90];
};

#endif

