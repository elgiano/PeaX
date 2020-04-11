// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "FFT_UGens.h"
#include "PeaXShift2.hpp"

namespace Trax {

PeaXShift2::PeaXShift2(): m_tempbuf(0) {
    out0(0) = in0(0);
    mCalcFunc = make_calc_function<PeaXShift2, &PeaXShift2::next>();
    next(1);
}

void PeaXShift2::init(SndBuf* buf, int numbins, float rOversampling){
  common = new PeaX_Common(mWorld->mFullRate.mSampleRate, numbins, buf->samples,rOversampling);

	m_tempbuf = (float*)RTAlloc(mWorld, buf->samples * sizeof(float));
	m_inPrevPhases = (float*)RTAlloc(mWorld, numbins * sizeof(float));
	m_outPrevPhases = (double*)RTAlloc(mWorld, numbins * sizeof(double));
	m_inFreqs = (double*)RTAlloc(mWorld, numbins * sizeof(double));
	data_a = new PeaX_Data(mWorld, numbins);

	memset(m_inPrevPhases, 0, numbins * sizeof(float));
	memset(m_outPrevPhases, 0, numbins * sizeof(double));
}

PeaXShift2::~PeaXShift2(){
  if(m_tempbuf) RTFree(mWorld, m_tempbuf);
  if(m_inPrevPhases) RTFree(mWorld, m_inPrevPhases);
  if(m_outPrevPhases) RTFree(mWorld, m_outPrevPhases);
  if(m_inFreqs) RTFree(mWorld, m_inFreqs);
	if(data_a) delete data_a;
	if(common) delete common;
}

void PeaXShift2::next(int nSamples) {
  auto* unit = this;
  PV_GET_BUF2
	// in0(2): rOversampling
  if(!m_tempbuf) init(buf1,numbins, in0(2));

	m_rootFreq = in0(3);
	if(m_rootFreq<=0) m_rootFreq = 1;
	m_rootBinB = m_rootFreq * common->constants.rFreqPerBin;
  m_ampThr = in0(4);
  m_maxPeaks = in0(5);
	if(m_maxPeaks < 0 || m_maxPeaks > numbins) m_maxPeaks = numbins ;

  SCPolarBuf* fft_a = ToPolarApx(buf1);
  SCPolarBuf* fft_b = ToPolarApx(buf2);
  SCPolarBuf* t = (SCPolarBuf*) m_tempbuf;
	memset(m_tempbuf, 0, buf1->samples * sizeof(float));

  // Print("getFreqs\n");
  common->getFreqs(fft_b, m_inPrevPhases, m_inFreqs);
	// Print("find Peaks");
  common->findPeaks(data_a, fft_a, m_ampThr, m_maxPeaks);
  // Print("combine\n");
  transpose(data_a, fft_b, t);

  memcpy(fft_a, t, buf1->samples * sizeof(float));
}

void PeaXShift2::transpose(PeaX_Data* data, SCPolarBuf* b, SCPolarBuf* dstBuf){
	int numpeaks = data->numpeaks;
	int numbins = common->constants.numbins;
	double ampmult = common->constants.ampmult, rAmpmult = common->constants.rAmpmult;

	if(numpeaks==0) return;
	//Print("[c] %d\n", numpeaks);
	PeaX_Peak *peaks = data->peaks;
	double avgFreq, freqRatio, mag, sourceMag;
	int numFreqs, transposingBin;

	for (int i=0; i < numbins; ++i) {
		avgFreq = 0; numFreqs = 0; mag = 0;
		for (int j=0; j< numpeaks; ++j) {
			freqRatio =  m_rootFreq/data->peaks[j].freq;
			transposingBin = (i  - data->peaks[j].bin) * freqRatio + m_rootBinB;
		  // Print("b %d , r %f\n", transposingBin , rFreqRatio);
			if( transposingBin >= numbins ) break;
			else if( transposingBin < 0 ) continue;
		  // Print("b %f , r %f", m_inFreqs[transposingBin] , rFreqRatio);
			avgFreq += m_inFreqs[transposingBin] / freqRatio;
			sourceMag = (double) b->bin[transposingBin].mag * ampmult;
			mag += (double) data->peaks[j].mag / (double) data->maxMag *  sourceMag ;
			++numFreqs;
		}
		mag *= rAmpmult;
		if(numFreqs==0 || mag == 0){
			dstBuf->bin[i].phase = 	dstBuf->bin[i].mag = m_outPrevPhases[i] = 0;;
			continue;
		}
		avgFreq = (double) avgFreq / numFreqs;
		// Print("freq %d %f %f\n", i, avgFreq, mag );
		dstBuf->bin[i].phase = (float) common->getPhaseCorrection(i, avgFreq, m_outPrevPhases);
		dstBuf->bin[i].mag = mag;
	}
}

} // namespace Trax

InterfaceTable* ft;

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::PeaXShift2>(ft, "PeaXShift2");
}
