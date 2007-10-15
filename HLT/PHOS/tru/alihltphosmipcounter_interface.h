#ifndef __ALIHLTPHOSMIPCOUNTER_INTERFACE_H
#define __ALIHLTPHOSMIPCOUNTER_INTERFACE_H

class AliHLTPHOSMIPCounterInterface
{
public:
    AliHLTPHOSMIPCounterInterface() {}
    virtual ~AliHLTPHOSMIPCounterInterface() {}


private:
    AliHLTPHOSMIPCounterInterface( const AliHLTPHOSMIPCounterInterface& source );
    void operator = ( const AliHLTPHOSMIPCounterInterface& source );
};


#endif // __ALIHLTPHOSMIPCOUNTER_INTERFACE_H
