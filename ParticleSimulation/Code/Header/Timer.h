#ifndef _TIMER_H_
#define _TIMER_H_

#include "../Header/3DLibs.h"

class Timer
{
	private:
		__int64		mCntsPerSec;
		double		mSecsPerCnt;
		__int64		mPrevTimeStamp;
		int			mFrameCnt;
		float		mFps;
		float		mDeltaTime;
		float		mFrameTime;
		__int64		mFrameStart;
		bool		mIsRunning;

		//Timer 
		__int64		mStartTime;
		__int64		mStopTime;

	public:
		float		GetDeltaTime();
		float		GetFPS();

		void		StartTimer();
		double		StopTimer();
		bool		IsRunning() const;
		void		Update();

		HRESULT		Initialize();
		void		Release();
					Timer();
		virtual		~Timer();
};
#endif