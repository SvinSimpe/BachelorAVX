#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_

#include "../Header/3DLibs.h"
#include <vector>

//Get AVX instrinsics
#include <immintrin.h>


#include "../Header/Timer.h" 

class ParticleSystem
{
	private:
		float*		mXPosition;
		float*		mYPosition;
		float*		mZPosition;

		float*		mXVelocity;
		float*		mYVelocity;
		float*		mZVelocity;

		float*		mInitialXVelocity;
		float*		mInitialYVelocity;
		float*		mInitialZVelocity;

		ID3D11Buffer*	mVertexBuffer;
		float*			mVertices;		// Used for maping, SOA -> AOS-pattern
		float*			mTempPtr;		// Used for shuffle mapping

		ID3D11Buffer*	mPerObjectCBuffer;
		PerObjectData	mPerObjectData;

		// Timer
		Timer mTimer;

	private:
		HRESULT		UpdateAndSetBuffer( ID3D11DeviceContext* deviceContext );
		HRESULT		UpdateAndSetPerObjectData( ID3D11DeviceContext* deviceContext, float deltaTime );
		void		UpdateParticleLogic( float deltaTime );

		void		SetRandomVelocity( int index );
		void		CheckDeadParticles();
		void		ResetParticle( int index );

	public:
		void Update( float deltaTime );
		void Render( ID3D11DeviceContext* deviceContext );

		HRESULT Initialize( ID3D11Device* device, XMFLOAT3 emitterPosition );
		void Release();

		ParticleSystem();
		~ParticleSystem();
};
#endif