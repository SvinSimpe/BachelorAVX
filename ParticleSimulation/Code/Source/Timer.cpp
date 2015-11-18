#include "../Header/Timer.h"

float Timer::GetDeltaTime()
{
	return mDeltaTime;
}

float Timer::GetFPS()
{
	return mFps;
}

void Timer::StartTimer()
{
	mStartTime = 0;
	QueryPerformanceCounter( (LARGE_INTEGER*)&mStartTime );
	mIsRunning		= true;
}

double Timer::StopTimer()
{
	mStopTime = 0;
	QueryPerformanceCounter( (LARGE_INTEGER*)&mStopTime );
	mIsRunning	= false;

	return (double)( ( mStopTime - mStartTime ) * mSecsPerCnt );
}

bool Timer::IsRunning() const
{
	return mIsRunning;
}

void Timer::Update()
{
	__int64 currTimeStamp = 0;
	QueryPerformanceCounter( (LARGE_INTEGER*)&currTimeStamp );

	mDeltaTime = (float)( ( currTimeStamp - mPrevTimeStamp ) * mSecsPerCnt );

	mPrevTimeStamp = currTimeStamp;

	mFrameTime -= mDeltaTime;
	if( mFrameTime <= 0.0f )
	{	
		mFps = (float)mFrameCnt / (float)( currTimeStamp - mFrameStart );
		mFps = (float)mFrameCnt;
		mFrameCnt = 0; // Reset frame count
		QueryPerformanceCounter( (LARGE_INTEGER*)&mFrameStart ); // Store frame start time
		mFrameTime = 1.0f;
	}
	else
		mFrameCnt++;
}

HRESULT Timer::Initialize()
{
	mCntsPerSec = 0;
	QueryPerformanceFrequency( (LARGE_INTEGER*)&mCntsPerSec );

	mSecsPerCnt = 1.0 / (double)mCntsPerSec;

	mPrevTimeStamp = 0;
	QueryPerformanceCounter( (LARGE_INTEGER*)&mPrevTimeStamp );

	mFps		= 0.0f;
	mFrameTime  = 1.0f; // 1 sec
	
	QueryPerformanceCounter( (LARGE_INTEGER*)&mFrameStart );

	

	mFrameCnt = 0;

	mIsRunning	= false;
	
	return S_OK;
}

void Timer::Release()
{
	mCntsPerSec		= 0;
	mSecsPerCnt		= 0.0f;
	mPrevTimeStamp	= 0;
	mFrameCnt		= 0;
	mFps			= 0.0f;
	mFrameTime		= 0.0f;
	mFrameStart		= 0;
	mIsRunning		= false;
}

Timer::Timer()
{
	mCntsPerSec		= 0;
	mSecsPerCnt		= 0.0f;
	mPrevTimeStamp	= 0;
	mFrameCnt		= 0;
	mFps			= 0.0f;
	mFrameTime		= 0.0f;
	mFrameStart		= 0;
	mIsRunning		= false;
}

Timer::~Timer()
{
}