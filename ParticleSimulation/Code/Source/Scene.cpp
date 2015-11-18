#include "../Header/Scene.h"

HRESULT Scene::CompileShader( char* shaderFile, char* pEntrypoint, char* pTarget, D3D10_SHADER_MACRO* pDefines, ID3DBlob** pCompiledShader )
{
    HRESULT hr = S_OK;
 
    DWORD dwShaderFlags =  D3DCOMPILE_DEBUG;// D3DCOMPILE_ENABLE_STRICTNESS |
                           // D3DCOMPILE_IEEE_STRICTNESS;
 
    std::string shader_code;
    std::ifstream in( shaderFile, std::ios::in | std::ios::binary );
 
    if ( in )
    {
        in.seekg( 0, std::ios::end );
        shader_code.resize( (unsigned int)in.tellg() );
        in.seekg( 0, std::ios::beg );
        in.read( &shader_code[0], shader_code.size() );
        in.close();
    }
 
    ID3DBlob* pErrorBlob = nullptr;
    hr = D3DCompile( shader_code.data(),
                             shader_code.size(),
                             NULL,
                             pDefines,
                             nullptr,
                             pEntrypoint,
                             pTarget,
                             dwShaderFlags,
                             NULL,
                             pCompiledShader,
                             &pErrorBlob );
 
    if( pErrorBlob )
    {
        OutputDebugStringA( (char*)pErrorBlob->GetBufferPointer() );
        SAFE_RELEASE( pErrorBlob );
    }
 
 
    return hr;
}

HRESULT Scene::InitializeShaders( ID3D11Device* device )
{
	//TEST 
	UINT yOffset = NUM_PARTICLES;
	UINT zOffset = NUM_PARTICLES * 2;


    HRESULT hr = S_OK;
 
    //-------------------------------
    // Compile Particle Vertex Shader
    //-------------------------------
    ID3DBlob* pvs = nullptr;

	 
     
	if ( SUCCEEDED( hr = CompileShader( "Code/Shader/ParticleShader.hlsl", "VS_main", "vs_5_0", nullptr, &pvs ) ) )
    {
        if( SUCCEEDED( hr = device->CreateVertexShader( pvs->GetBufferPointer(),
                                                          pvs->GetBufferSize(),
                                                          nullptr,
                                                          &mVertexShader ) ) )
        {
            D3D11_INPUT_ELEMENT_DESC inputDesc[] = {                
                { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
            };
 
            hr = device->CreateInputLayout( inputDesc,
                                            ARRAYSIZE( inputDesc ),
                                            pvs->GetBufferPointer(),
                                            pvs->GetBufferSize(),
                                            &mInputLayout );
        }
 
        SAFE_RELEASE( pvs );
    }
 
 
    //---------------------------------
    // Compile Particle Geometry Shader
    //-------------------------------- 
    ID3DBlob* pgs = nullptr;
 
    if( SUCCEEDED( hr = CompileShader( "Code/Shader/ParticleShader.hlsl", "GS_main", "gs_5_0", nullptr, &pgs ) ) )
    {
        hr = device->CreateGeometryShader( pgs->GetBufferPointer(),
                                           pgs->GetBufferSize(),
                                           nullptr,
                                           &mGeometryShader );
        SAFE_RELEASE( pgs );
    }
 
 
    //------------------------------
    // Compile Particle Pixel Shader
    //------------------------------
    ID3DBlob* pps = nullptr;
 
    if( SUCCEEDED( hr = CompileShader( "Code/Shader/ParticleShader.hlsl", "PS_main", "ps_5_0", nullptr, &pps ) ) )
    {
        hr = device->CreatePixelShader( pps->GetBufferPointer(),
                                        pps->GetBufferSize(),
                                        nullptr,
                                        &mPixelShader );
        SAFE_RELEASE( pps );
    }
 
 
    return hr;
}

HRESULT Scene::UpdateCBPerFrame()
{
	HRESULT hr = S_OK;

	// Retrive view and projection matrices from camera to per frame buffer data
	XMStoreFloat4x4( &mPerFrameData.viewMatrix, XMMatrixTranspose( XMLoadFloat4x4( &mCamera.GetViewMatrix() ) ) );
	XMStoreFloat4x4( &mPerFrameData.projMatrix, XMMatrixTranspose( XMLoadFloat4x4( &mCamera.GetProjectionMatrix() ) ) );

	// Get camera eye position and add to per frame data.
	// Must have to perform Geometry shader billboarding
	mPerFrameData.eyePosition = mCamera.GetEyePosition();

	D3D11_MAPPED_SUBRESOURCE mappedResource;
	hr = mDeviceContext->Map( mPerFrameCBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource );

	if( SUCCEEDED( hr ) )
	{
		memcpy( mappedResource.pData, &mPerFrameData, sizeof(PerFrameData) );	
		mDeviceContext->Unmap( mPerFrameCBuffer, 0 );

		// Set constant buffer to Geometry shader stage
		mDeviceContext->GSSetConstantBuffers( 0, 1, &mPerFrameCBuffer );
	}

	return hr;
}

void Scene::Update( float deltaTime )
{
	UpdateCBPerFrame();

	mParticleSystem.Update( deltaTime );
}

void Scene::Render()
{
	//Set shader stages
	mDeviceContext->VSSetShader( mVertexShader, nullptr, 0 );
	mDeviceContext->HSSetShader( nullptr, nullptr, 0 );
	mDeviceContext->DSSetShader( nullptr, nullptr, 0 );
	mDeviceContext->GSSetShader( mGeometryShader, nullptr, 0 );
	mDeviceContext->PSSetShader( mPixelShader, nullptr, 0 );

	// Set position vertex description
	mDeviceContext->IASetInputLayout( mInputLayout );

	// Set topology
	mDeviceContext->IASetPrimitiveTopology( D3D11_PRIMITIVE_TOPOLOGY_POINTLIST );

	mParticleSystem.Render( mDeviceContext );

}

HRESULT Scene::Initialize( ID3D11Device* device, ID3D11DeviceContext* deviceContext )
{
	HRESULT hr = S_OK;

	mDeviceContext = deviceContext;

	if( FAILED( hr = mParticleSystem.Initialize( device, XMFLOAT3( 0.0f, 2.0f, 0.0f ) ) ) )
		return hr;

	if( FAILED( hr = InitializeShaders( device ) ) )
		return hr;

	D3D11_BUFFER_DESC cbDesc;
	cbDesc.ByteWidth			= sizeof( mPerFrameData );
	cbDesc.Usage				= D3D11_USAGE_DYNAMIC;
	cbDesc.BindFlags			= D3D11_BIND_CONSTANT_BUFFER;
	cbDesc.CPUAccessFlags		= D3D11_CPU_ACCESS_WRITE;
	cbDesc.MiscFlags			= 0;
	cbDesc.StructureByteStride	= 0;

	hr = device->CreateBuffer( &cbDesc, nullptr, &mPerFrameCBuffer );

	return hr;
}

void Scene::Release()
{
	mParticleSystem.Release();

	SAFE_RELEASE( mVertexShader );
	SAFE_RELEASE( mGeometryShader );
	SAFE_RELEASE( mPixelShader );
	SAFE_RELEASE( mInputLayout );
	SAFE_RELEASE( mPerFrameCBuffer );
}

Scene::Scene()
{
	mDeviceContext			= nullptr;

	mVertexShader			= nullptr;
	mGeometryShader			= nullptr;
	mPixelShader			= nullptr;
	mInputLayout			= nullptr;

	mPerFrameCBuffer		= nullptr;
}

Scene::~Scene()
{}