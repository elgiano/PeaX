// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "../PeaXCommon/PeaXCommon.hpp"

namespace Trax {

class PeaX : public SCUnit {
public:
    PeaX();

    ~PeaX();

private:
    // Calc function
    void next(int nSamples);
    void init(SndBuf* buf, int numbins, float rOversampling);
		void combinePeaks(PeaX_Data* a, PeaX_Data* b, SCPolarBuf* dstBuf);

		PeaX_Common* common;

    float* m_tempbuf;
    PeaX_Data* data_a;
    PeaX_Data* data_b;
    double *m_outPrevPhases;

		// in0(0) is bufferA
    // in0(1) is bufferB
    // in0(2) is oversampling
		// in0(3-4) is rootFreqA and B
		int m_rootFreqA, m_rootBinA,
		    m_rootFreqB, m_rootBinB;
		// in0(5-6) is ampThresholdA and B
		float m_ampThrA, m_ampThrB;
		// in0(7-8) is maxPeaksA and B
		int m_maxPeaksA, m_maxPeaksB;
};

} // namespace Trax
