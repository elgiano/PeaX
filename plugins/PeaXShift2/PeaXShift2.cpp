// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "FFT_UGens.h"
#include "PeaXShift2.hpp"

namespace Trax {

///////////////////////////////////////////////////

double PeaXShift2::getTrueFreq(int bin_i, float phaseDiff){
	double diff = phaseDiff - (double)bin_i*m_expct;
	// map delta phase into +/- Pi interval
	long qpd = diff/pi;
	if (qpd >= 0) qpd += qpd&1;
	else qpd -= qpd&1;
	diff -= pi*(double)qpd;
	// get deviation from bin frequency from the +/- Pi interval
	diff /= m_rOversampling*twopi;
	// compute the k-th partials' true frequency
	return (double) bin_i*m_freqPerBin + diff*m_freqPerBin;
}

double PeaXShift2::getPhaseCorrection(int bin_i, double freq){
    // subtract bin mid frequency
		double tmp = freq - (double)bin_i*m_freqPerBin;
    // get bin deviation from freq deviation
		tmp *= m_rFreqPerBin;
    // take osamp into account
		tmp *= twopi*m_rOversampling;
    // add the overlap phase advance back in
		tmp += (double)bin_i*m_expct;
    // accumulate delta phase to get bin phase
		m_phaseAccumulator[bin_i] += tmp;
		return m_phaseAccumulator[bin_i];
}

//////////////////////////////////////////////////

PeaXShift2::PeaXShift2(): m_tempbuf(0) {
    out0(0) = in0(0);
    mCalcFunc = make_calc_function<PeaXShift2, &PeaXShift2::next>();
    next(1);
}

void PeaXShift2::init(SndBuf* buf, int numbins){
  m_numbins = numbins;
  m_ampmult = (double) 2.f/(numbins+1); //*(1.0/unit->m_maxpeaks);
  m_rAmpmult = (double) (numbins+1)/2.f; //*(1.0/unit->m_maxpeaks);
  m_windowsize = buf->samples;

  m_rOversampling = in0(2);
  m_hopsize = m_rOversampling * m_windowsize;

  m_freqPerBin = (double) mWorld->mFullRate.mSampleRate / m_windowsize;
  m_rFreqPerBin = (double) m_windowsize / mWorld->mFullRate.mSampleRate;
  m_expct = twopi * (double) m_rOversampling;

  m_tempbuf = (float*)RTAlloc(mWorld, buf->samples * sizeof(float));
  m_prevPhases = (float*)RTAlloc(mWorld, numbins * sizeof(float));
  m_phaseAccumulator = (float*)RTAlloc(mWorld, numbins * sizeof(float));
  m_freqsB = (double*)RTAlloc(mWorld, numbins * sizeof(double));
  data_a = new PeaXShift2_Data(mWorld, numbins);
  data_b = new PeaXShift2_Data(mWorld, numbins);

  memset(m_prevPhases, 0, numbins * sizeof(float));
	memset(m_phaseAccumulator, 0, numbins * sizeof(float));

	if(mWorld->mFullRate.mBufLength!=64)
    Print("PeaXShift2 complains: block size not 64, you have %d\n", mWorld->mFullRate.mBufLength);
}

PeaXShift2::~PeaXShift2(){
  if(data_a) delete data_a;
  if(data_b) delete data_b;
  if(m_tempbuf) RTFree(mWorld, m_tempbuf);
  if(m_prevPhases) RTFree(mWorld, m_prevPhases);
  if(m_phaseAccumulator) RTFree(mWorld, m_phaseAccumulator);
  if(m_freqsB) RTFree(mWorld, m_freqsB);
}

void PeaXShift2::next(int nSamples) {
  auto* unit = this;
  PV_GET_BUF2

  if (!m_tempbuf)
      init(buf1,numbins);
  else if (numbins != m_numbins)
      return;

  m_maxpeaks = in0(3);
	data_a->rootFreq = in0(4);
	m_rootFreq = in0(5);
	m_rootBinB = m_rootFreq * m_rFreqPerBin;
  m_tolerance = in0(6);
  data_a->ampThr = in0(7);
  data_b->ampThr = in0(7);

  SCPolarBuf* fft_a = ToPolarApx(buf1);
  SCPolarBuf* fft_b = ToPolarApx(buf2);
  SCPolarBuf* t = (SCPolarBuf*) m_tempbuf;

	memset(m_tempbuf, 0, m_windowsize * sizeof(float));

  // Print("findPeaks\n");
	// Print("[a] ");
  // findPeaks(data_a, fft_a);
	getFreqs(fft_b, data_b->prevPhases);
	// Print("[b] ");
  findPeaks(data_a, fft_a);
  // Print("combine\n");
  // combinePeaks(data_a, data_b, t);
	transpose(data_a, fft_b, t);
  memcpy(fft_a, t, m_windowsize * sizeof(float));
}

void PeaXShift2::getFreqs(SCPolarBuf* fft, float* prevPhases){
	float currPhase;
	for(int i = 0; i < m_numbins; i++){
		currPhase = fft->bin[i].phase;
		m_freqsB[i] = getTrueFreq(i, currPhase - prevPhases[i]);
		prevPhases[i] = currPhase;
	}
}

void PeaXShift2::findPeaks(PeaXShift2_Data* data, SCPolarBuf* fft){

	//swap new peaks to old; current now safe to overwrite;
	PeaXShift2_peak *peaks = data->peaks;
	int numpeaks = 0;

  // angular frequency is pi*(i/nover2)
  // double angmult = pi/(m_numbins+1);
  float ampcheck = data->ampThr * m_rAmpmult;

	float phase, prevmag, mag, nextmag;

	prevmag = fft->bin[0].mag;
	mag = fft->bin[1].mag;
	float maxMag, maxFreq;
	if(prevmag > mag && prevmag > ampcheck){
		peaks[0].mag = prevmag * m_ampmult;
		peaks[0].freq = getTrueFreq(0, phase - data->prevPhases[0]);
		peaks[0].bin = 0;
		maxMag = prevmag;
		maxFreq = peaks[0].freq;
		++numpeaks;
	}else if(mag > fft->bin[2].mag && mag > ampcheck){
		peaks[0].mag = mag * m_ampmult;
		peaks[0].freq = getTrueFreq(1, phase - data->prevPhases[1]);
		peaks[0].bin = 1;
		maxMag = mag;
		maxFreq = peaks[0].freq;
		++numpeaks;
	}

	//could restrict not to go above nover4 (numbins/2)!
	for (int i=1; i< m_numbins; ++i) {

		phase= fft->bin[i].phase;
		nextmag= fft->bin[i].mag;

		if (
      (prevmag<mag) && (nextmag<mag) && (mag>ampcheck)
    ) {
			// found a peak
			peaks[numpeaks].mag = mag * m_ampmult;
      float phaseDev = phase - data->prevPhases[i];
      double freq = getTrueFreq(i, phaseDev);
			peaks[numpeaks].freq = freq;
			peaks[numpeaks].bin = i;
			if(peaks[numpeaks].mag>maxMag){
				maxMag = peaks[numpeaks].mag;
				maxFreq = freq;
			}
			++numpeaks;
		}

		data->prevPhases[i] = phase;
		prevmag=mag;
		mag=nextmag;
	}

	// Print("%d peaks, max: %fHz %f\n",numpeaks,maxFreq,maxMag);

	data->maxMag = maxMag;
	data->maxFreq = maxFreq;
	data->numpeaks = numpeaks;
}

void PeaXShift2::transpose(PeaXShift2_Data* data, SCPolarBuf* b, SCPolarBuf* dstBuf){
	int numpeaks = data->numpeaks;
	if(numpeaks==0) return;
	//Print("[c] %d\n", numpeaks);
	PeaXShift2_peak *peaks = data->peaks;
	double rRootFreq = 1.f/m_rootFreq;
	float sourceMag;
	double avgFreq, rFreqRatio, mag;
	int numFreqs, transposingBin;

	for (int i=0; i < m_numbins; ++i) {
		avgFreq = 0; numFreqs = 0; mag = 0;
		for (int j=0; j< numpeaks; ++j) {
			rFreqRatio =  data->peaks[j].freq * rRootFreq;
			transposingBin = (i  - data->peaks[j].bin) * rFreqRatio + m_rootBinB;
			// Print("b %d , r %f\n", transposingBin , rFreqRatio);
			if( transposingBin < 0 || transposingBin >= m_numbins) continue;
			// Print("b %f , r %f", m_freqsB[transposingBin] , rFreqRatio);
			avgFreq += m_freqsB[transposingBin] * rFreqRatio;
			sourceMag = b->bin[transposingBin].mag * m_ampmult;
			mag += (double) data->peaks[j].mag / (double) data->maxMag * (double) sourceMag;
			++numFreqs;
		}
		avgFreq = (double) avgFreq / numFreqs;
		// Print("freq %d %f %f\n", i, avgFreq, mag * m_rAmpmult);
		dstBuf->bin[i].phase = (float) getPhaseCorrection(i, avgFreq);
		dstBuf->bin[i].mag = mag * m_rAmpmult;
	}
}

} // namespace Trax

InterfaceTable* ft;

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::PeaXShift2>(ft, "PeaXShift2");
}
