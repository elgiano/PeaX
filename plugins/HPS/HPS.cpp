// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "FFT_UGens.h"
#include "HPS.hpp"

#define PV_GET_BUF_ANAL                                                                                                \
    float fbufnum = in0(0);                                                                                           \
    if (fbufnum < 0.f) {                                                                                               \
        out0(0) = outval; return;                                                                               \
        return;                                                                                                        \
    }                                                                                                                  \
    uint32 ibufnum = (uint32)fbufnum;                                                                                  \
    World* world = mWorld;                                                                                       \
    SndBuf* buf;                                                                                                       \
    if (ibufnum >= world->mNumSndBufs) {                                                                               \
        int localBufNum = ibufnum - world->mNumSndBufs;                                                                \
        Graph* parent = mParent;                                                                                 \
        if (localBufNum <= parent->localBufNum) {                                                                      \
            buf = parent->mLocalSndBufs + localBufNum;                                                                 \
        } else {                                                                                                       \
            buf = world->mSndBufs;                                                                                     \
        }                                                                                                              \
    } else {                                                                                                           \
        buf = world->mSndBufs + ibufnum;                                                                               \
    }                                                                                                                  \
    LOCK_SNDBUF(buf);                                                                                                  \
    int numbins = (buf->samples - 2) >> 1;


namespace Trax {

HPS::HPS(): m_prevPhases(0), m_init(false) {
    mCalcFunc = make_calc_function<HPS, &HPS::next>();
    next(1);
}
HPS::~HPS(){
    if(m_prevPhases != 0) RTFree(mWorld, m_prevPhases);
}

void HPS::next(int nSamples) {
  PV_GET_BUF_ANAL

  if(!m_init){
    //Print("Init\n");
    int windowsize = buf->samples;
    m_oversampling = in0(1);
    m_freqPerBin = world->mFullRate.mSampleRate / windowsize;
    m_rFreqPerBin = windowsize / world->mFullRate.mSampleRate;
    m_expct = twopi * (double) m_oversampling;
    m_prevPhases = (float*) RTAlloc(world, numbins * sizeof(float));
    memset(m_prevPhases, 0, numbins * sizeof(float));
    m_init = true;
  }

  SCPolarBuf* fft = ToPolarApx(buf);

  int harmonics = in0(2);
  int minbin =  in0(3) * m_rFreqPerBin;
  int maxbin =  in0(4) * m_rFreqPerBin;
  int thr =  in0(5);
  bool correction = in0(6);
  if(minbin<=0) minbin = 0;
  if(maxbin<=0 || maxbin >= numbins) maxbin = (numbins-1);
  // Print("args: %d -> %d (%d)\n",minbin,maxbin, harmonics);
  // Print("calc: %f (%f1)\n",m_freqPerBin, m_oversampling);

  int fund_bin = hps(fft->bin, numbins, minbin, maxbin, harmonics, thr, correction);
  float bin_freq = (float) fund_bin*m_freqPerBin;
  float freq =  (float) getTrueFreq(fund_bin, fft->bin[fund_bin].phase - m_prevPhases[fund_bin]);
  for(int i=0; i < numbins; i++){
    m_prevPhases[i] = fft->bin[i].phase;
  }
  //Print("%d -> %f (%f)\n",fund_bin,freq,bin_freq);

  out0(0) = outval = freq ;
}

double HPS::getTrueFreq(int bin_i, float phaseDiff){
	double diff = phaseDiff - (double)bin_i*m_expct;
	// map delta phase into +/- Pi interval
	long qpd = diff/pi;
	if (qpd >= 0) qpd += qpd&1;
	else qpd -= qpd&1;
	diff -= pi*(double)qpd;
	// get deviation from bin frequency from the +/- Pi interval
	diff /= m_oversampling*twopi;
	// compute the k-th partials' true frequency
	return (double) bin_i*m_freqPerBin + diff*m_freqPerBin;
}

int HPS::hps(SCPolar* spectrum, int numbins, int minbin, int maxbin, int harmonics=8, float thr=0, bool correction=1){
  // set the maximum index to search for a pitch
   int i, j, maxHIndex;
   maxHIndex = numbins/harmonics;
   if (maxbin < maxHIndex) maxHIndex = maxbin;

   // generate the Harmonic Product Spectrum values and keep track of the
   // maximum amplitude value to assign to a pitch.

   int maxLocation = minbin;
   for (j=minbin; j<=maxbin; j++) {
		  if(spectrum[j].mag < thr) continue;
      for (i=2; i<=harmonics; i++) {
         if(j*i>maxbin) continue;
         spectrum[j].mag += spectrum[j*i].mag;
      }
      if (spectrum[j].mag > spectrum[maxLocation].mag) {
         maxLocation = j;
      }
   }

   // Correct for octave too high errors.  If the maximum subharmonic
   // of the measured maximum is approximately 1/2 of the maximum
   // measured frequency, AND if the ratio of the sub-harmonic
   // to the total maximum is greater than 0.2, THEN the pitch value
   // is assigned to the subharmonic.

   if(correction){
     int max2 = minbin;
     int maxsearch = maxLocation * 3 / 4;
     for (i=minbin+1; i<maxsearch; i++) {
        if (spectrum[i].mag > spectrum[max2].mag) {
           max2 = i;
        }
     }
     if (abs(max2 * 2 - maxLocation) < 4) {
        if (spectrum[max2].mag / spectrum[maxLocation].mag > 0.2) {
           maxLocation = max2;
        }
     }
   }

   return maxLocation;
}

} // namespace Trax

InterfaceTable* ft;

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::HPS>(ft, "HPS");
}
