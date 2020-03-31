// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace Trax {

//freq in this case is omega, angular frequency in radians per sample = 2*PI*f/SR
struct PeaXShift2_peak {
	int bin;
	float mag, phase;
	double freq;
};

struct PeaXShift2_Data{

	//use buffer swapping of pointer as needed
	PeaXShift2_peak *peaks;
	int numpeaks;

  World *world;

  float *prevPhases;
  float ampThr;
	float maxMag, maxFreq, rootFreq;

  PeaXShift2_Data(World* world, int numbins):
    world(world), numpeaks(0), ampThr(0)
  {
    peaks = (PeaXShift2_peak*)RTAlloc(world, numbins * sizeof(PeaXShift2_peak));
    prevPhases = (float*)RTAlloc(world, numbins * sizeof(float));
    memset(prevPhases, 0, numbins * sizeof(float));
  }

  ~PeaXShift2_Data(){
    if(peaks) RTFree(world, peaks);
    if(prevPhases) RTFree(world, prevPhases);
  }
};

class PeaXShift2 : public SCUnit {
public:
    PeaXShift2();

    ~PeaXShift2();

private:
    // Calc function
    void next(int nSamples);
    void init(SndBuf* buf, int numbins);
    void findPeaks(PeaXShift2_Data* data, SCPolarBuf* fft);
		void combinePeaks(PeaXShift2_Data* a, PeaXShift2_Data* b, SCPolarBuf* dstBuf);
    double getTrueFreq(int bin_i, float phaseDiff);
    double getPhaseCorrection(int bin_i, double freq);
		void getFreqs(SCPolarBuf* fft, float* prevPhases);
		void transpose(PeaXShift2_Data* data, SCPolarBuf* a, SCPolarBuf* dstBuf);

    float* m_tempbuf;
    PeaXShift2_Data* data_a;
    PeaXShift2_Data* data_b;
		double *m_freqsB;
    float *m_prevPhases;
    float *m_phaseAccumulator;

    // in0(0) is buf
    float m_rOversampling;
    int m_maxpeaks;
    int m_numpeaksrequested;
    float m_tolerance;
    // in0(5) is ampThr

		int m_rootFreq;
		int m_rootBinB;

    int m_numbins;
    int m_windowsize;
    int m_hopsize;
    double m_freqPerBin, m_rFreqPerBin;
    double m_expct;
    double m_ampmult;
    double m_rAmpmult;

};

} // namespace Trax
