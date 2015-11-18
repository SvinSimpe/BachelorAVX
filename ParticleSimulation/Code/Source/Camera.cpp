#include "../Header/Camera.h"

Camera::Camera()
{
	mEyePosition		= XMFLOAT3( -360.0f, 0.0f, -360.0f );
	mFocusPosition		= XMFLOAT3( 0.0f, 0.0f,   1.0f );
	mUpDirection		= XMFLOAT3( 0.0f, 1.0f,   0.0f );
	mRightDirection		= XMFLOAT3( 1.0f, 0.0f,   0.0f );

	XMMATRIX tempViewMatrix = XMMatrixLookAtLH( XMLoadFloat3( &mEyePosition ), 
												XMLoadFloat3( &mFocusPosition ),
												XMLoadFloat3( &mUpDirection ) );
	XMStoreFloat4x4( &mViewMatrix, tempViewMatrix );


	XMMATRIX tempProjectionMatrix = XMMatrixPerspectiveFovLH( 0.75f, 1920.0f / 1080.0f, 0.5f, 5000.0f );
	XMStoreFloat4x4( &mProjectionMatrix, tempProjectionMatrix );
}
Camera::~Camera()
{

}

XMFLOAT4X4 Camera::GetViewMatrix() const
{
	return mViewMatrix;
}

XMFLOAT4X4 Camera::GetProjectionMatrix() const
{
	return mProjectionMatrix;
}

XMFLOAT3 Camera::GetEyePosition() const
{
	return mEyePosition;
}

void Camera::UpdateViewMatrix()
{
	XMMATRIX tempViewMatrix = XMMatrixLookAtLH( XMLoadFloat3( &mEyePosition ), 
												XMLoadFloat3( &mFocusPosition ),
												XMLoadFloat3( &mUpDirection ) );
	XMStoreFloat4x4( &mViewMatrix, tempViewMatrix );
}