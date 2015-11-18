#ifndef _D3DAPP_H_
#define _D3DAPP_H_

#include "../Header/Timer.h"
#include "../Header/Scene.h"

const unsigned int SCREEN_WIDTH		= 1920;
const unsigned int SCREEN_HEIGHT	= 1080;

class D3DApp
{
	private:
		//Window
		HINSTANCE	mHInstance;
		HWND		mHWnd;

		//DirectX
		ID3D11Device*				mDevice;
		ID3D11DeviceContext*		mDeviceContext;

		IDXGISwapChain*				mSwapChain;
		ID3D11RenderTargetView*		mRenderTargetView;
		ID3D11DepthStencilView*		mDepthStencilView;

		ID3D11RasterizerState*		mRasterizerStateSolid;
		ID3D11RasterizerState*		mRasterizerStateWired;
		ID3D11RasterizerState*		mCurrentRasterizerState;

		//Simulation
		Timer mTimer;
		Scene mScene;

	private:
		void SetFPSAndDeltaTimeInWindow();

	public:
		void Render();

		static LRESULT CALLBACK WndProc( HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam );

		D3DApp();
		~D3DApp();

		HRESULT InitializeWindow( HINSTANCE hInstance,  int nCmdShow );
		HRESULT InitializeDirectX11();
		HRESULT InitializeSimulationComponents();

		int Run();
		void Release();
};
#endif