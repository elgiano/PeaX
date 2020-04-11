// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "FFT_UGens.h"
#include "PeaX.hpp"

namespace Trax {

PeaX::PeaX(): m_tempbuf(0) {
    out0(0) = in0(0);
    mCalcFunc = make_calc_function<PeaX, &PeaX::next>();
    next(1);
}

void PeaX::init(SndBuf* buf, int numbins, float rOversampling){
	common = new PeaX_Common(mWorld->mFullRate.mSampleRate, numbins, buf->samples,rOversampling);

  m_tempbuf = (float*)RTAlloc(mWorld, buf->samples * sizeof(float));
  m_outPrevPhases = (double*)RTAlloc(mWorld, numbins * sizeof(double));
  data_a = new PeaX_Data(mWorld, numbins);
  data_b = new PeaX_Data(mWorld, numbins);

	memset(m_outPrevPhases, 0, numbins * sizeof(double));
}

PeaX::~PeaX(){
	if(m_tempbuf) RTFree(mWorld, m_tempbuf);
  if(m_outPrevPhases) RTFree(mWorld, m_outPrevPhases);

	if(data_a) delete data_a;
	if(data_b) delete data_b;
	if(common) delete common;
}

void PeaX::next(int nSamples) {
  auto* unit = this;
  PV_GET_BUF2
	// in0(2): rOversampling
  if(!m_tempbuf) init(buf1,numbins, in0(2));

	m_rootFreqA = in0(3);
	m_rootFreqB = in0(4);
	m_ampThrA = in0(5);
	m_ampThrB = in0(6);
	m_maxPeaksA = in0(7);
	if(m_maxPeaksA < 0 || m_maxPeaksA > numbins) m_maxPeaksA = numbins ;
	m_maxPeaksB = in0(8);
	if(m_maxPeaksB < 0 || m_maxPeaksB > numbins) m_maxPeaksB = numbins ;

  SCPolarBuf* fft_a = ToPolarApx(buf1);
  SCPolarBuf* fft_b = ToPolarApx(buf2);
  SCPolarBuf* t = (SCPolarBuf*) m_tempbuf;
	memset(m_tempbuf, 0, buf1->samples * sizeof(float));

  // Print("findPeaks\n");
	// Print("[a] ");
  common->findPeaks(data_a, fft_a, m_ampThrA, m_maxPeaksA);
	// Print("[b] ");
  common->findPeaks(data_b, fft_b, m_ampThrB, m_maxPeaksB);
  // Print("combine\n");
  combinePeaks(data_a, data_b, t);

  memcpy(fft_a, t, buf1->samples * sizeof(float));
}

void PeaX::combinePeaks(PeaX_Data* a, PeaX_Data* b, SCPolarBuf* dstBuf){

	int numpeaks_a = a->numpeaks;
	int numpeaks_b = b->numpeaks;
	// Print("[a] %d peaks, max: %fHz %f\n",numpeaks_a,a->maxFreq,a->maxMag);
	// Print("[b] %d peaks, max: %fHz %f\n",numpeaks_b,b->maxFreq,b->maxMag);

	if(numpeaks_a*numpeaks_b==0) return;
	// Print("[c] %d\n", numpeaks_a*numpeaks_b);
	PeaX_Peak *peaks_a = a->peaks;
	PeaX_Peak *peaks_b = b->peaks;

  // float minPhase, maxPhase, minMag, maxMag;
	float magRatio = a->maxMag/b->maxMag;
	float currentMagRatio;
	// Print("magratio: %f\n",magRatio);

	double rFreqPerBin = common->constants.rFreqPerBin;
	int numbins = common->constants.numbins;

	for (int i=0; i < numpeaks_a; ++i) {
		float freqRatio = a->peaks[i].freq/m_rootFreqB;
		currentMagRatio = peaks_a[i].mag * magRatio;
		// Print("frqratio %f\n", freqRatio);
		for (int j=0; j< numpeaks_b; ++j) {
			double freq = b->peaks[j].freq * freqRatio;
			int index = (int) freq * rFreqPerBin;
			if( index < 0 || index >= numbins) continue;
			float phase = common->getPhaseCorrection(index, freq, m_outPrevPhases);
			float mag = b->peaks[j].mag * currentMagRatio;
			dstBuf->bin[index].mag += mag;
			dstBuf->bin[index].phase = phase;
		}
	}

	double rAmpmult = common->constants.rAmpmult;
	for(int i=0; i< numbins; i++){
		dstBuf->bin[i].mag *= rAmpmult;
	}
}

} // namespace Trax

InterfaceTable* ft;

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::PeaX>(ft, "PeaX");
}
