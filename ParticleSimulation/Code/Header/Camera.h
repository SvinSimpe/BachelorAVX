#ifndef _CAMERA_H_
#define _CAMERA_H_
 
#include "../Header/3DLibs.h"
 
class Camera
{
    private:
        XMFLOAT3 mEyePosition;
        XMFLOAT3 mFocusPosition;
        XMFLOAT3 mUpDirection;
        XMFLOAT3 mRightDirection;

        XMFLOAT4X4 mViewMatrix;
        XMFLOAT4X4 mProjectionMatrix;
 
    public:
        Camera();
        ~Camera();
 
        XMFLOAT4X4 GetViewMatrix() const;
        XMFLOAT4X4 GetProjectionMatrix() const;
        XMFLOAT3 GetEyePosition() const;
 
        void UpdateViewMatrix();
};
 
#endif
