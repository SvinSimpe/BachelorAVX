#include "../Header/ParticleSystem.h"
#include <time.h>
#include <sstream>
#include <string>


HRESULT ParticleSystem::UpdateAndSetBuffer( ID3D11DeviceContext* deviceContext )
{
	D3D11_MAPPED_SUBRESOURCE mappedResource;
	HRESULT hr = deviceContext->Map( mVertexBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource );

	if( SUCCEEDED( hr ) )
	{
		size_t counter	= 0;
		//for (int baseIndex = 0; baseIndex < 125; baseIndex++)
		for (int baseIndex = 0; baseIndex < NUM_PARTICLES * 0.125f; baseIndex++) // NUM_PARTICLES * 0.125f == NUM_PARTICLES / 8
		{
			//Mapping from SOA-pattern to AOS-pattern courtesy of Intel
			//https://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx

			//Load
			__m256 xReg = _mm256_load_ps( &mXPosition[baseIndex * 8] );
			__m256 yReg = _mm256_load_ps( &mYPosition[baseIndex * 8] );
			__m256 zReg = _mm256_load_ps( &mZPosition[baseIndex * 8] );


			//Shuffle
			__m256 xyReg = _mm256_shuffle_ps( xReg, yReg, _MM_SHUFFLE( 2,0,2,0 ) );
			__m256 yzReg = _mm256_shuffle_ps( yReg, zReg, _MM_SHUFFLE( 3,1,3,1 ) );
			__m256 zxReg = _mm256_shuffle_ps( zReg, xReg, _MM_SHUFFLE( 3,1,2,0 ) );

			__m256 reg03 = _mm256_shuffle_ps( xyReg, zxReg, _MM_SHUFFLE( 2, 0, 2, 0 ) );
			__m256 reg14 = _mm256_shuffle_ps( yzReg, xyReg, _MM_SHUFFLE( 3, 1, 2, 0 ) );
			__m256 reg25 = _mm256_shuffle_ps( zxReg, yzReg, _MM_SHUFFLE( 3, 1, 3, 1 ) );
			


			
			//Map, xyz
			__m128* vertexRegAOS = (__m128*)mTempPtr;

			vertexRegAOS[0] = _mm256_castps256_ps128( reg03 );	// x8,y8,z8,x7
			vertexRegAOS[1] = _mm256_castps256_ps128( reg14 );	// y7,z7,x6,y6
			vertexRegAOS[2] = _mm256_castps256_ps128( reg25 );	// z6,x5,y5,z5

			vertexRegAOS[3] = _mm256_extractf128_ps( reg03, 1 );	// x4,y4,z4,x3
			//vertexRegAOS[4] = _mm256_extractf128_ps( reg14, 1 );	// y3,z3,x2,y2
			//vertexRegAOS[5] = _mm256_extractf128_ps( reg25, 1 );	// z2,x1,y1,z1

			//CORDES========================================================
			__m256 reg45 = _mm256_permute2f128_ps (reg14, reg25, 1|(3<<4) );
			_mm256_storeu_ps( (float*)(vertexRegAOS + 4), reg45);
			//CORDES========================================================

			/*memcpy( mappedResource.pData, &vertexRegAOS[0], sizeof( ParticleVertex12 ) * NUM_PARTICLES );	*/

			

			for ( int index = 0, subIndex = 0 ; index < 6; index++ )
			{
				mVertices[counter++] = vertexRegAOS[index].m128_f32[(subIndex++) % 4];
				mVertices[counter++] = vertexRegAOS[index].m128_f32[(subIndex++) % 4];
				mVertices[counter++] = vertexRegAOS[index].m128_f32[(subIndex++) % 4];
				mVertices[counter++] = vertexRegAOS[index].m128_f32[(subIndex++) % 4];
			}

		}
		


		/// ================ EASY APPROACH ===================
		//int counter = 0;
		//for (int i = 0; i < NUM_PARTICLES / 8; i++)
		//{
		//	__m256 xPositionReg = _mm256_load_ps( &mXPosition[i * 8] );
		//	__m256 yPositionReg = _mm256_load_ps( &mYPosition[i * 8] );
		//	__m256 zPositionReg = _mm256_load_ps( &mZPosition[i * 8] );

		//	for (size_t i = 0; i < 8; i++)
		//	{
		//		mVertices[counter++] = xPositionReg.m256_f32[i];
		//		mVertices[counter++] = yPositionReg.m256_f32[i];
		//		mVertices[counter++] = zPositionReg.m256_f32[i];
		//	}	

		//	_mm256_zeroall();
		//}
		/// ================ EASY APPROACH ===================




		memcpy( mappedResource.pData, mVertices, sizeof( ParticleVertex12 ) * NUM_PARTICLES );
		deviceContext->Unmap( mVertexBuffer, 0 );
	}

	return hr;
}

HRESULT ParticleSystem::UpdateAndSetPerObjectData( ID3D11DeviceContext* deviceContext, float deltaTime )
{
	HRESULT hr = S_OK;

	XMStoreFloat4x4( &mPerObjectData.worldMatrix, XMMatrixTranspose( XMMatrixIdentity() ) );

	D3D11_MAPPED_SUBRESOURCE mappedResource;
	hr = deviceContext->Map( mPerObjectCBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource );

	if( SUCCEEDED( hr ) )
	{
		memcpy( mappedResource.pData, &mPerObjectData, sizeof( mPerObjectData ) );	
		deviceContext->Unmap( mPerObjectCBuffer, 0 );

		// Set constant buffer to shader stages
		deviceContext->GSSetConstantBuffers( 1, 1, &mPerObjectCBuffer );	
	}

	return hr;
}

void ParticleSystem::UpdateParticleLogic( float deltaTime )
{
	const __m256 accelerationReg	= _mm256_set1_ps( -9.82f );
	const __m256 deltaReg			= _mm256_set1_ps( deltaTime * 100.0f );

	for ( int i = 0; i < NUM_PARTICLES / 8; i++ )
	{
		//Load positions into registers
		__m256 xPositionReg = _mm256_load_ps( &mXPosition[i*8] );
		__m256 yPositionReg = _mm256_load_ps( &mYPosition[i*8] );
		__m256 zPositionReg = _mm256_load_ps( &mZPosition[i*8] );

		//Load velocity into registers
		__m256 xVelocityReg = _mm256_load_ps( &mXVelocity[i*8] );
		__m256 yVelocityReg = _mm256_load_ps( &mYVelocity[i*8] );
		__m256 zVelocityReg = _mm256_load_ps( &mZVelocity[i*8] );

		// Calculate new Y-velocity based on acceleration and timestep (deltaTime)
		__m256 newYVelocityReg = _mm256_add_ps( _mm256_mul_ps( deltaReg , accelerationReg ),  yVelocityReg );

		// Store new Y-velocity
		yVelocityReg = newYVelocityReg;

		// Calculate new positions   |  oldPosition + ( velocity * timeStep )
		__m256 newXPositionReg = _mm256_add_ps( _mm256_mul_ps( deltaReg , xVelocityReg ), xPositionReg  );
		__m256 newYPositionReg = _mm256_add_ps( _mm256_mul_ps( deltaReg , yVelocityReg ), yPositionReg );
		__m256 newZPositionReg = _mm256_add_ps( _mm256_mul_ps( deltaReg , zVelocityReg ), zPositionReg );

		// Store new positions back to memory
		_mm256_store_ps( &mXPosition[i*8], newXPositionReg );
		_mm256_store_ps( &mYPosition[i*8], newYPositionReg );
		_mm256_store_ps( &mZPosition[i*8], newZPositionReg );

		// Store new velocity back to memory
		_mm256_store_ps( &mYVelocity[i*8], yVelocityReg );

	}

	

	// -9.2559631349317831e+061 - 1024 particles

	// 0.00010142520525020329 - 10 024 particles

	// 0.00082951328579630550 - 100 024 particles

	// 0.0061422621322355262 - 10^6 + 24 particles


	int k = 4;
}

void ParticleSystem::SetRandomVelocity( int index )
{
	XMFLOAT3 randomVelocity[8];
	float randomSpreadAngle = 0.0f;

	for ( int i = 0; i < 8; i++ )
	{
		randomSpreadAngle = (float)( ( rand() % 3000 * 2 ) - 3000 ) * 0.01f; // 3000 == 100 * 30 degree spread angle
		XMVECTOR randomAimingDirection	= XMVector3TransformCoord( XMLoadFloat3( &XMFLOAT3( 0.0f, 1.0f, 0.0f ) ), XMMatrixRotationX( XMConvertToRadians( randomSpreadAngle ) ) );
		
		// Test a magnitude
		int magnitude = ( ( rand() % 15 ) + 38 );// - 5;
		randomAimingDirection *= magnitude;

		XMStoreFloat3( &randomVelocity[i], XMVector3TransformCoord( randomAimingDirection, XMMatrixRotationZ( XMConvertToRadians( randomSpreadAngle ) ) ) );	
		
	}

	//Set random velocity to registers
	__m256 xVelReg = _mm256_set_ps( randomVelocity[0].x, randomVelocity[1].x,
									randomVelocity[2].x, randomVelocity[3].x,
									randomVelocity[4].x, randomVelocity[5].x,
									randomVelocity[6].x, randomVelocity[7].x );

	__m256 yVelReg = _mm256_set_ps( randomVelocity[0].y, randomVelocity[1].y,
			 						randomVelocity[2].y, randomVelocity[3].y,
			  						randomVelocity[4].y, randomVelocity[5].y,
									randomVelocity[6].y, randomVelocity[7].y );
	
	__m256 zVelReg = _mm256_set_ps( randomVelocity[0].z, randomVelocity[1].z,
			 						randomVelocity[2].z, randomVelocity[3].z,
			  						randomVelocity[4].z, randomVelocity[5].z,
									randomVelocity[6].z, randomVelocity[7].z );

	//Store random velocity in Velocity members
	_mm256_store_ps( &mXVelocity[index], xVelReg );
	_mm256_store_ps( &mYVelocity[index], yVelReg );
	_mm256_store_ps( &mZVelocity[index], zVelReg );

	//Store random velocity in Initial Velocity members
	_mm256_store_ps( &mInitialXVelocity[index], xVelReg );
	_mm256_store_ps( &mInitialYVelocity[index], yVelReg );
	_mm256_store_ps( &mInitialZVelocity[index], zVelReg );
}

void ParticleSystem::CheckDeadParticles()
{
	const __m256 zeroReg	= _mm256_set1_ps( -200.0f ); 
	//const __m256 zeroReg	= _mm256_setzero_ps();                 // const vector of zeros
	
	for ( int i = 0; i + 8 <= NUM_PARTICLES; i += 8 )
	{
	    __m256 yPositionReg	= _mm256_loadu_ps( &mYPosition[i] );					// load 8 x floats
	    __m256 cmpReg		= _mm256_cmp_ps( yPositionReg, zeroReg, _CMP_LE_OS );	// compare for <= 0
	    int	mask			= _mm256_movemask_ps( cmpReg );							// get MS bits from comparison result
	    if ( mask != 0 )															// if any bits set
	    {																			// then we have 1 or more elements <= 0
	        for ( int k = 0; k < 8; ++k )											// test each element in vector
	        {																		// using scalar code...
	            if ( ( mask & 1 ) != 0 )
	            {
	                // found element at index i + k
	                // do something with it...
					ResetParticle( i + k );
	            }
	            mask >>= 1;
	        }
	    }
	}
	//// deal with any remaining elements in case where n is not a multiple of 8
	//for ( int j = i; j < n; ++j )
	//{
	//    if ( mYPosition[j] <= 0.0f )
	//    {
	//        // found element at index j
	//        // do something with it...
	//		ResetParticle( j );
	//    }
	//}
}

void ParticleSystem::ResetParticle( int index )
{
	mXPosition[index] = 0.0f;
	mYPosition[index] = 0.0f;
	mZPosition[index] = 0.0f;

	mXVelocity[index] = mInitialXVelocity[index];
	mYVelocity[index] = mInitialYVelocity[index];
	mZVelocity[index] = mInitialZVelocity[index];
}

void ParticleSystem::Update( float deltaTime )
{
	deltaTime *= 0.01f;

	if( mTimer.IsRunning() )
	{
		// Print Timer data
		std::wstringstream stream;
		stream << "Time per frame:  " << mTimer.StopTimer() << " ms";
		OutputDebugStringW( stream.str().c_str() );
		OutputDebugStringA( "\n" );
	}
	mTimer.StartTimer();


	UpdateParticleLogic( deltaTime );
	CheckDeadParticles();
}

void ParticleSystem::Render( ID3D11DeviceContext* deviceContext )
{
	UpdateAndSetBuffer( deviceContext );
	UpdateAndSetPerObjectData( deviceContext, 0.0f );

	deviceContext->VSSetConstantBuffers( 1, 1, &mPerObjectCBuffer);
	deviceContext->GSSetConstantBuffers( 1, 1, &mPerObjectCBuffer);
	deviceContext->PSSetConstantBuffers( 1, 1, &mPerObjectCBuffer);

	UINT32 offset				 	= 0;
	UINT32 stride					= sizeof(ParticleVertex12);
	ID3D11Buffer* buffersToSet[]	= { mVertexBuffer };
	deviceContext->IASetVertexBuffers( 0, 1, buffersToSet, &stride, &offset );

	deviceContext->Draw( NUM_PARTICLES, 0 );
}

HRESULT ParticleSystem::Initialize( ID3D11Device* device, XMFLOAT3 emitterPosition )
{
	//Allocating aligned memory for each of the attribute component. Allocated aligned
	//memory should be deallocated using '_aligned_free( void *memblock )'
	mXPosition = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );
	mYPosition = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );
	mZPosition = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );

	mXVelocity = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );
	mYVelocity = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );
	mZVelocity = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );

	mInitialXVelocity = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );
	mInitialYVelocity = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );
	mInitialZVelocity = (float*) _aligned_malloc( NUM_PARTICLES * sizeof(float), 32 );

	srand(time(NULL));

	//Zeroes all YMM registers
	_mm256_zeroall();

	//Since the AVX instructions used below is fetching 32bytes or a set of 8 single-precision
	//floating points from memory to registers, we increment loop index with 8 per iteration
	for ( size_t i = 0; i < NUM_PARTICLES; i+=8 )
	{
		// Zero position members
		const __m256 zeroReg = _mm256_setzero_ps();

		_mm256_store_ps( &mXPosition[i], zeroReg );
		_mm256_store_ps( &mYPosition[i], zeroReg );
		_mm256_store_ps( &mZPosition[i], zeroReg );


		//Generate eight random velocities each call. The functions utilizes
		//AVX-intrinsics for filling YMM-registers with random velocities and
		//store into members for both Velocity and Initial Velocity
		SetRandomVelocity( i );
	}

	//Build Vertex Buffer
	D3D11_BUFFER_DESC desc;
	ZeroMemory( &desc, sizeof( desc ) );

    desc.BindFlags				= D3D11_BIND_VERTEX_BUFFER;
    desc.Usage					= D3D11_USAGE_DYNAMIC;
    desc.ByteWidth				= sizeof(ParticleVertex12) * NUM_PARTICLES;
	desc.StructureByteStride	= sizeof(ParticleVertex12);
	desc.CPUAccessFlags			= D3D11_CPU_ACCESS_WRITE;
	
	int k = ( NUM_PARTICLES * 3 ) * sizeof(float);

	//Allocating aligned memory for array used for maping vertices to buffer
	mVertices = (float*) _aligned_malloc( ( NUM_PARTICLES * 3 ) * sizeof(float), 32 );
	mTempPtr  = new float[24];

	D3D11_SUBRESOURCE_DATA subData;
	subData.pSysMem = mVertices;

	if( FAILED( device->CreateBuffer( &desc, &subData, &mVertexBuffer ) ) )
		return E_FAIL;


	D3D11_BUFFER_DESC cbDesc;
	cbDesc.ByteWidth			= sizeof( PerObjectData );
	cbDesc.Usage				= D3D11_USAGE_DYNAMIC;
	cbDesc.BindFlags			= D3D11_BIND_CONSTANT_BUFFER;
	cbDesc.CPUAccessFlags		= D3D11_CPU_ACCESS_WRITE;
	cbDesc.MiscFlags			= 0;
	cbDesc.StructureByteStride	= 0;


	if( FAILED( device->CreateBuffer( &cbDesc, nullptr, &mPerObjectCBuffer ) ) )
		return E_FAIL;

	// Timer
	mTimer.Initialize();

	return S_OK;
}

void ParticleSystem::Release()
{
	_aligned_free( mXPosition );
	_aligned_free( mYPosition );
	_aligned_free( mZPosition );

	_aligned_free( mXVelocity );
	_aligned_free( mYVelocity );
	_aligned_free( mZVelocity );

	_aligned_free( mInitialXVelocity );
	_aligned_free( mInitialYVelocity );
	_aligned_free( mInitialZVelocity );

	SAFE_RELEASE( mVertexBuffer );
	_aligned_free( mVertices );
	SAFE_DELETE( mTempPtr );

	//Zeroes all YMM registers
	_mm256_zeroall();

	SAFE_RELEASE( mPerObjectCBuffer );
}

ParticleSystem::ParticleSystem()
{
	mXPosition		= nullptr;
	mYPosition		= nullptr;
	mZPosition		= nullptr;
			  
	mXVelocity		= nullptr;
	mYVelocity		= nullptr;
	mZVelocity		= nullptr;

	mVertexBuffer	= nullptr;
	mVertices		= nullptr;
	mTempPtr		= nullptr;

	mPerObjectCBuffer	= nullptr;
}

ParticleSystem::~ParticleSystem()
{}

