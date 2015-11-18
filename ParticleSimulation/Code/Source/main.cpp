#include "../Header/D3DApp.h"

#include <crtdbg.h>
//Get AVX instrinsics
#include <immintrin.h>
//Get CPUID capability
#include <intrin.h>
//Code for detection of AVX support found at https://insufficientlycomplicated.wordpress.com/2011/11/07/detecting-intel-advanced-vector-extensions-avx-in-visual-studio/
// 2015-04-21, 15:50
bool DetectAXVSupport()
{
	bool avxSupported = false;
 
    // If Visual Studio 2010 SP1 or later
#if ( _MSC_FULL_VER >= 160040219 )
    // Checking for AVX requires 3 things:
    // 1) CPUID indicates that the OS uses XSAVE and XRSTORE
    //     instructions (allowing saving YMM registers on context
    //     switch)
    // 2) CPUID indicates support for AVX
    // 3) XGETBV indicates the AVX registers will be saved and
    //     restored on context switch
    //
    // Note that XGETBV is only available on 686 or later CPUs, so
    // the instruction needs to be conditionally run.
    int cpuInfo[4]; 
    __cpuid( cpuInfo, 1 );
 
    bool osUsesXSAVE_XRSTORE = cpuInfo[2] & (1 << 27) || false;
    bool cpuAVXSuport = cpuInfo[2] & (1 << 28) || false;
 
    if ( osUsesXSAVE_XRSTORE && cpuAVXSuport )
    {
        // Check if the OS will save the YMM registers
        unsigned long long xcrFeatureMask = _xgetbv( _XCR_XFEATURE_ENABLED_MASK );
        //avxSupported = (xcrFeatureMask & 0x6) || false;
		//avxSupported = (xcrFeatureMask & 0x6) == 6;
		
		avxSupported = ( (xcrFeatureMask & 0x6) == 6 ) && ( (xcrFeatureMask & 0x6) || false );
    }
#endif
 
    if ( avxSupported )
		//AVX is supported
        return true;

    else
		//AVX is NOT supported
        return false;
 
    return 0;

 }


int WINAPI wWinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPWSTR lpCmdLine, int nCmdShow )
{

#if defined(DEBUG) | defined(_DEBUG)

	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_CrtSetBreakAlloc(200); // Break at specific memory allocation point

#endif

	//First check for AVX support
	if( !DetectAXVSupport() )
	{
		MessageBox( 0, "System does not support AVX.", "Error!", MB_OK );
		return 0;
	}
	
	int retVal = 0;

	D3DApp* app = new D3DApp();

	if( FAILED( retVal = app->InitializeWindow( hInstance, nCmdShow ) ) )
		return retVal;

	if( FAILED( retVal = app->InitializeDirectX11() ) )
		return retVal;

	if( FAILED( retVal = app->InitializeSimulationComponents() ) )
		return retVal;


	retVal = app->Run();

	app->Release();

	delete app;

	return retVal;
}