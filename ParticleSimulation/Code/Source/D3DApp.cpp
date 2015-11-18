#include "../Header/D3DApp.h"
#include <sstream>
#include <string>

void D3DApp::SetFPSAndDeltaTimeInWindow()
{
	std::wstringstream wss;
	std::wstring ws;

	wss.clear();
			
	wss << "FPS: " << mTimer.GetFPS() << "  DeltaTime: " << mTimer.GetDeltaTime();
	ws = wss.str();
	
	SetWindowTextW( mHWnd, ws.c_str() );
	ws.clear();
	wss.str(ws);
}

void D3DApp::Render()
{
    // Clear Back Buffer
    static float clearColor[4] = { 0.3f, 0.2f, 0.5f, 1.0f };
    mDeviceContext->ClearRenderTargetView( mRenderTargetView, clearColor );
 
    // Clear Depth Buffer
    mDeviceContext->ClearDepthStencilView( mDepthStencilView, D3D11_CLEAR_DEPTH, 1.0f, 0 );
 
	// Set Rasterizer State
	mDeviceContext->RSSetState( mCurrentRasterizerState );
	
	//SetViewPort();
	mDeviceContext->OMSetRenderTargets( 1, &mRenderTargetView, mDepthStencilView );

	//Render scene
	mScene.Render();

    // Swap Front and Back Bufffer
    mSwapChain->Present( 0, 0 );
}

LRESULT CALLBACK D3DApp::WndProc( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam )
{
    PAINTSTRUCT ps;
    HDC hdc;
 
    switch ( msg )
    {
        case WM_PAINT:
            hdc = BeginPaint( hWnd, &ps );
            EndPaint( hWnd, &ps );
            break;
 
        case WM_DESTROY:
            PostQuitMessage( 0 );
            break;
 
        case WM_KEYDOWN:
            switch ( wParam )
            {
                case VK_ESCAPE:
                    PostQuitMessage( 0 );
                    break;
            }
            break;
         
        default:
            return DefWindowProc( hWnd, msg, wParam, lParam );
    }
    return 0;
}

D3DApp::D3DApp()
{
    mHInstance     = NULL;
    mHWnd          = NULL;
 
    mDevice            = nullptr;
    mDeviceContext     = nullptr;
    mSwapChain         = nullptr;
    mRenderTargetView  = nullptr;
    mDepthStencilView  = nullptr;
}

D3DApp::~D3DApp()
{}

HRESULT D3DApp::InitializeWindow( HINSTANCE hInstance,  int nCmdShow )
{
    HRESULT hr = S_OK;
 
    mHInstance = hInstance;
 
    //----------------------
    // Register Window Class
    //----------------------
    WNDCLASSEX wc;
    wc.cbSize           = sizeof( WNDCLASSEX );
    wc.style            = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc      = WndProc;
    wc.cbClsExtra       = 0;
    wc.cbWndExtra       = 0;
    wc.hInstance        = mHInstance;
    wc.hIcon            = 0;
    wc.hCursor          = LoadCursor( NULL, IDC_HAND );
    wc.hbrBackground    = (HBRUSH)( COLOR_WINDOW + 1 );
    wc.lpszMenuName     = NULL;
    wc.lpszClassName    = "ParticleSimulation";
    wc.hIconSm          = 0;
 
    if( FAILED( hr = RegisterClassEx( &wc ) ) )
        return hr;
 
    //-----------------------
    // Adjust & Create Window
    //-----------------------
    RECT rc = { 0, 0, SCREEN_WIDTH , SCREEN_HEIGHT };
    AdjustWindowRect( &rc, WS_OVERLAPPEDWINDOW, FALSE );
 
    if( !( mHWnd = CreateWindow(	"ParticleSimulation",
                                    "ParticleSimulation",
                                    WS_OVERLAPPEDWINDOW,
                                    CW_USEDEFAULT,
                                    CW_USEDEFAULT,
                                    rc.right - rc.left,
                                    rc.bottom - rc.top,
                                    NULL,
                                    NULL,
                                    mHInstance,
                                    NULL ) ) )
    {
        return E_FAIL;
    }
 
    ShowWindow( mHWnd, nCmdShow );
    ShowCursor( TRUE );

 
    return hr;
}

HRESULT D3DApp::InitializeDirectX11()
{
	HRESULT hr = E_FAIL;
 
    RECT rc;
    GetClientRect( mHWnd, &rc );
 
    int width   = rc.right - rc.left;
    int height  = rc.bottom - rc.top;
 
 
    //-------------------------------------------
    // Create Swap Chain, Device & Device Context
    //------------------------------------------- 
    D3D_DRIVER_TYPE driverTypes[] = { D3D_DRIVER_TYPE_HARDWARE,
                                      D3D_DRIVER_TYPE_WARP,
                                      D3D_DRIVER_TYPE_REFERENCE };
 
    DXGI_SWAP_CHAIN_DESC sd;
    memset( &sd, 0, sizeof( sd ) );
    sd.BufferCount              = 1;
    sd.BufferDesc.Width                     = width;
    sd.BufferDesc.Height                    = height;
    sd.BufferDesc.Format                    = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.RefreshRate.Numerator     = 0;
    sd.BufferDesc.RefreshRate.Denominator   = 1;
    sd.BufferUsage                          = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.OutputWindow                         = mHWnd;
    sd.SampleDesc.Count                     = 1;
    sd.Windowed                             = TRUE;
 
    D3D_FEATURE_LEVEL featureLevelsToTry[] = { D3D_FEATURE_LEVEL_11_0 };
    D3D_FEATURE_LEVEL initiatedFeatureLevel;

	UINT creationFlags = D3D11_CREATE_DEVICE_BGRA_SUPPORT;
	#if defined(_DEBUG)
        // If the project is in a debug build, enable the debug layer.
        creationFlags |= D3D11_CREATE_DEVICE_DEBUG;
	#endif
 
    for ( UINT driverTypeIndex = 0; driverTypeIndex < ARRAYSIZE( driverTypes ) && FAILED( hr ); driverTypeIndex++ )
    {

        hr = D3D11CreateDeviceAndSwapChain( nullptr,
                                            driverTypes[driverTypeIndex],
                                            NULL,
											creationFlags,
                                            featureLevelsToTry,
                                            ARRAYSIZE( featureLevelsToTry ),
                                            D3D11_SDK_VERSION,
                                            &sd,
                                            &mSwapChain,
                                            &mDevice,
                                            &initiatedFeatureLevel,
                                            &mDeviceContext );
    }
 
    if( FAILED( hr ) )
        return hr;
 
 
    //--------------------------
    // Create Render Target View
    //--------------------------
    ID3D11Texture2D* pBackBuffer;
    if ( SUCCEEDED( hr = mSwapChain->GetBuffer( 0, __uuidof( ID3D11Texture2D ), (void**)&pBackBuffer ) ) )
    {
        hr = mDevice->CreateRenderTargetView( pBackBuffer, nullptr, &mRenderTargetView );
        SAFE_RELEASE( pBackBuffer );
    }
     
    if( FAILED( hr ) )
        return hr;
 
 
    //--------------------------
    // Create Depth Stencil View
    //-------------------------- 
    D3D11_TEXTURE2D_DESC dsd;
    dsd.Width               = width;
    dsd.Height              = height;
    dsd.MipLevels           = 1;
    dsd.ArraySize           = 1;
    dsd.Format              = DXGI_FORMAT_D32_FLOAT;
    dsd.SampleDesc.Count    = 1;
    dsd.SampleDesc.Quality  = 0;
    dsd.Usage               = D3D11_USAGE_DEFAULT;
    dsd.BindFlags           = D3D11_BIND_DEPTH_STENCIL;
    dsd.CPUAccessFlags      = 0;
    dsd.MiscFlags           = 0;
 
	ID3D11Texture2D* depthStencil = nullptr;

    if( FAILED( hr = mDevice->CreateTexture2D( &dsd, nullptr, &depthStencil ) ) )
        return hr;
 
    D3D11_DEPTH_STENCIL_VIEW_DESC dsvd;
    ZeroMemory( &dsvd, sizeof( dsvd ) );
    dsvd.Format             = dsd.Format;
    dsvd.ViewDimension      = D3D11_DSV_DIMENSION_TEXTURE2D;
    dsvd.Texture2D.MipSlice = 0;
 
    if( FAILED( hr = mDevice->CreateDepthStencilView( depthStencil, &dsvd, &mDepthStencilView ) ) )
        return hr;
 
    //-------------
    // Set Viewport
    //-------------
    D3D11_VIEWPORT vp;
    vp.Width    = (float)width;
    vp.Height   = (float)height;
    vp.MinDepth = 0.0f;
    vp.MaxDepth = 1.0f;
    vp.TopLeftX = 0;
    vp.TopLeftY = 0;
    mDeviceContext->RSSetViewports( 1, &vp );
 
    mDeviceContext->OMSetRenderTargets( 1, &mRenderTargetView, mDepthStencilView );
 
	//Set Rasterizer States
	//=====================
	//------------------------------
	// Create Solid Rasterizer State
	//------------------------------
	D3D11_RASTERIZER_DESC solidRazDesc;
	memset( &solidRazDesc, 0, sizeof( solidRazDesc ) );
	solidRazDesc.FillMode			= D3D11_FILL_SOLID;
	solidRazDesc.CullMode			= D3D11_CULL_NONE;
	solidRazDesc.DepthClipEnable	= true;

	hr = mDevice->CreateRasterizerState( &solidRazDesc, &mRasterizerStateSolid );

	// Set current rasterizerstate to solid
	mCurrentRasterizerState = mRasterizerStateSolid;

	//----------------------------------
	// Create Wireframe Rasterizer State
	//----------------------------------
	D3D11_RASTERIZER_DESC wiredRazDesc;
	memset( &wiredRazDesc, 0, sizeof( wiredRazDesc ) );
	wiredRazDesc.FillMode			= D3D11_FILL_WIREFRAME;
	wiredRazDesc.CullMode			= D3D11_CULL_NONE;
	wiredRazDesc.DepthClipEnable	= true;

	hr = mDevice->CreateRasterizerState( &wiredRazDesc, &mRasterizerStateWired );


	//Initialize Timer
	mTimer.Initialize();

    return hr;
}

HRESULT D3DApp::InitializeSimulationComponents()
{
	if( FAILED( mScene.Initialize( mDevice, mDeviceContext ) ) )
		return E_FAIL;

	return S_OK;
}

int D3DApp::Run()
{
	MSG msg = {0};
 
    while ( WM_QUIT != msg.message )
    {
        if( PeekMessage( &msg, NULL, 0, 0, PM_REMOVE ) )
        {
            TranslateMessage( &msg );
            DispatchMessage ( &msg );
        }
        else
        {
			mTimer.Update();
			mScene.Update( mTimer.GetDeltaTime() );
			Render();

			SetFPSAndDeltaTimeInWindow();
        }
    }
 
    return (int)msg.wParam;
}

void D3DApp::Release()
{
	SAFE_RELEASE( mDevice );
	SAFE_RELEASE( mDeviceContext );
	SAFE_RELEASE( mSwapChain );
	SAFE_RELEASE( mRenderTargetView );
	SAFE_RELEASE( mDepthStencilView );

	SAFE_RELEASE( mRasterizerStateSolid	);
	SAFE_RELEASE( mRasterizerStateWired	);
	SAFE_RELEASE( mCurrentRasterizerState );
	
	mTimer.Release();
	mScene.Release();
}