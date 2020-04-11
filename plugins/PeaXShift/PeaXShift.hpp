// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "../PeaXCommon/PeaXCommon.hpp"

namespace Trax {

class PeaXShift : public SCUnit {
public:
    PeaXShift();
    ~PeaXShift();

private:
    // Calc function
    void next(int nSamples);
    void init(SndBuf* buf, int numbins, float rOversampling);
		void transpose(PeaX_Data* data, SCPolarBuf* a, SCPolarBuf* dstBuf);

		PeaX_Common* common;

    float* m_tempbuf;
    PeaX_Data* data_b;
    double *m_inFreqs;
    float *m_inPrevPhases;
    double *m_outPrevPhases;

    // in0(0) is bufferA
    // in0(1) is bufferB
    // in0(2) is oversampling
		// in0(3) is rootFreq
		int m_rootFreq, m_rootBinB;
		// in0(4) is ampThreshold
		float m_ampThr;
		// in0(5) is maxPeaks
		int m_maxPeaks;
};

} // namespace Trax
