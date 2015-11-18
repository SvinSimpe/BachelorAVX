#ifndef _SCENE_H_
#define _SCENE_H_

#include "../Header/ParticleSystem.h"
#include "../Header/Camera.h"
#include <fstream>

class Scene
{
	private:
		ID3D11DeviceContext*		mDeviceContext;
		ParticleSystem				mParticleSystem;
		Camera						mCamera;

		ID3D11VertexShader*			mVertexShader;
		ID3D11GeometryShader*		mGeometryShader;
		ID3D11PixelShader*			mPixelShader;
		ID3D11InputLayout*			mInputLayout;	

		PerFrameData				mPerFrameData;
		ID3D11Buffer*				mPerFrameCBuffer;

	private:
		HRESULT CompileShader( char* shaderFile, char* pEntrypoint, char* pTarget, D3D10_SHADER_MACRO* pDefines, ID3DBlob** pCompiledShader );
		HRESULT InitializeShaders( ID3D11Device* device );
		HRESULT UpdateCBPerFrame();

	public:		
		void Update( float deltaTime );
		void Render();

		HRESULT	Initialize( ID3D11Device* device, ID3D11DeviceContext* deviceContext );
		void	Release();

		Scene();
		~Scene();
};
#endif