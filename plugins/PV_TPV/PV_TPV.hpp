// PluginTrax.hpp
// Gianluca Elia (elgiano@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace Trax {

//cubic interpolation of phase parameters for formula (37) where t is from 0 to 1 as interpolation parameter
struct PartialTrack {
	float theta1, omega1, theta2, omega2, alpha, beta; //cubic interpolation of phase
	float amp1, amp2; //linear interpolation of amplitude
};

//freq in this case is omega, angular frequency in radians per sample = 2*PI*f/SR
struct TPVPeak {
	float mag, phase;
	double freq; //improve frequency estimate by interpolation over amplitude of nearby points, or by time-frequency reassignment
};

struct PV_TPV_Data{
  //or for each partial to be rendered persistent data for rendering need phasek, angfreqk (omegak), alphak, betak as per (37)
	PartialTrack *tracks; //space for double maxpeaks if all birth and die at once!
	int numtracks;

	//use buffer swapping of pointer as needed
	TPVPeak *prevpeaks;
	TPVPeak *newpeaks;
	int numprevpeaks;
	int numnewpeaks;
  int maxpeaks;

  World *world;

  float *prevPhases;
  float ampThr;

  PV_TPV_Data(World* world, int maxpeaks, int numbins):
    world(world), maxpeaks(maxpeaks),
    numprevpeaks(0), numnewpeaks(0), numtracks(0),
    ampThr(0)
  {
    tracks = (PartialTrack*)RTAlloc(world, 2*maxpeaks * sizeof(PartialTrack));
    prevpeaks = (TPVPeak*)RTAlloc(world, maxpeaks * sizeof(TPVPeak));
    newpeaks = (TPVPeak*)RTAlloc(world, maxpeaks * sizeof(TPVPeak));
    prevPhases = (float*)RTAlloc(world, numbins * sizeof(float));
    memset(prevPhases, 0, numbins * sizeof(float));
  }

  ~PV_TPV_Data(){
    if(tracks) RTFree(world, tracks);
    if(prevpeaks) RTFree(world, prevpeaks);
    if(newpeaks) RTFree(world, newpeaks);
    if(prevPhases) RTFree(world, prevPhases);
  }
};

class PV_TPV : public SCUnit {
public:
    PV_TPV();

    ~PV_TPV();

private:
    // Calc function
    void next(int nSamples);
    void init(SndBuf* buf, int numbins);
    void findPeaks(PV_TPV_Data* data, SCPolarBuf* fft);
    void calcTracks(PV_TPV_Data* data);
    void writeTracks(int numtracks, PartialTrack* tracks, SCPolarBuf* fft);

    double getTrueFreq(int bin_i, float phaseDiff);
    double getPhaseCorrection(int bin_i, double freq);

    float* m_tempbuf;
    PV_TPV_Data* data;
    float *m_prevPhases;
    float *m_phaseAccumulator;

    // in0(0) is buf
    float m_rOversampling;
    int m_maxpeaks;
    int m_numpeaksrequested;
    float m_tolerance;
    // in0(5) is ampThr

    int m_numbins;
    int m_windowsize;
    int m_hopsize;
    double m_freqPerBin, m_rFreqPerBin;
    double m_expct;
    double m_ampmult;
    double m_rAmpmult;

};

} // namespace Trax
