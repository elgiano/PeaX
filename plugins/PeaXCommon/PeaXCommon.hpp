#pragma once

namespace Trax{

struct PeaX_Peak {
	int bin;
	float mag, phase;
	double freq;
};

struct PeaX_Data{

	int numpeaks;

	float maxMag, maxFreq, rootFreq;
	double peakMul;

	float *prevPhases;
	PeaX_Peak *peaks;
	World *world;

  PeaX_Data(World* world, int numbins):
    world(world), numpeaks(0)
  {
    peaks = (PeaX_Peak*)RTAlloc(world, numbins * sizeof(PeaX_Peak));
    prevPhases = (float*)RTAlloc(world, numbins * sizeof(float));
    memset(prevPhases, 0, numbins * sizeof(float));
  }

  ~PeaX_Data(){
		if(peaks) RTFree(world, peaks);
		if(prevPhases) RTFree(world, prevPhases);
  }
};

struct PeaX_Constants{

  int numbins;
  float rOversampling;
  double freqPerBin, rFreqPerBin,
         ampmult, rAmpmult,
         expct;


  PeaX_Constants(int sampleRate, int numbins, int windowsize, float rOversampling):
  rOversampling(rOversampling), numbins(numbins){

    ampmult = (double) 2.f/(numbins+1); //*(1.0/unit->m_maxpeaks);
    rAmpmult = (double) (numbins+1)/2.f; //*(1.0/unit->m_maxpeaks);
    freqPerBin = (double) sampleRate / windowsize;
    rFreqPerBin = (double) windowsize / sampleRate;
    expct = twopi * (double) rOversampling;

  }
};

class PeaX_Common{

public:

  PeaX_Constants constants;

  PeaX_Common(int sampleRate, int numbins, int windowsize, float rOversampling):
  constants(sampleRate, numbins, windowsize, rOversampling){

  };

  double getTrueFreq(int bin_i, float phaseDiff){
  	double diff = phaseDiff - (double)bin_i*constants.expct;
  	// map delta phase into +/- Pi interval
  	long qpd = diff/pi;
  	if (qpd >= 0) qpd += qpd&1;
  	else qpd -= qpd&1;
  	diff -= pi*(double)qpd;
  	// get deviation from bin frequency from the +/- Pi interval
  	diff /= constants.rOversampling*twopi;
  	// compute the k-th partials' true frequency
  	return (double) bin_i*constants.freqPerBin + diff*constants.freqPerBin;
  }

  double getPhaseCorrection(int bin_i, double freq, double* outPrevPhases){
      // subtract bin mid frequency
  		double tmp = freq - (double)bin_i*constants.freqPerBin;
      // get bin deviation from freq deviation
  		tmp *= constants.rFreqPerBin;
      // take osamp into account
  		tmp *= twopi*constants.rOversampling;
      // add the overlap phase advance back in
  		tmp += (double)bin_i*constants.expct;
      // accumulate delta phase to get bin phase
  		outPrevPhases[bin_i] += tmp;
  		return outPrevPhases[bin_i];
  }

  void getFreqs(SCPolarBuf* fft, float* prevPhases, double* freqs){
  	float currPhase;
  	for(int i = 0; i < constants.numbins; i++){
  		currPhase = fft->bin[i].phase;
  		freqs[i] = getTrueFreq(i, currPhase - prevPhases[i]);
  		prevPhases[i] = currPhase;
  	}
  }

  void findPeaks(PeaX_Data* data, SCPolarBuf* fft, float ampThr, int maxPeaks){

  	//swap new peaks to old; current now safe to overwrite;
  	PeaX_Peak *peaks = data->peaks;
  	int numpeaks = 0;

  	data->peakMul = 1;

    // angular frequency is pi*(i/nover2)
    // double angmult = pi/(m_numbins+1);
    double ampcheck = (double) ampThr * constants.rAmpmult;
  	float maxMag, maxFreq;

  	float phase = fft->bin[0].phase,
  				prevmag = fft->bin[0].mag,
  				mag = fft->bin[1].mag,
  				nextmag;

  	if(prevmag > mag && prevmag > ampcheck){
  		peaks[0].mag = prevmag * constants.ampmult;
  		peaks[0].freq = getTrueFreq(0, phase - data->prevPhases[0]);
  		peaks[0].bin = 0;
  		maxMag = prevmag;
  		maxFreq = peaks[0].freq;
  		data->peakMul *= peaks[0].mag;
  		++numpeaks;
  		data->prevPhases[0] = phase;
  	}

  	//could restrict not to go above nover4 (numbins/2)!
  	for (int i=2; i < constants.numbins; ++i) {

  		phase= fft->bin[i-1].phase;
  		nextmag= fft->bin[i].mag;

  		if (
        (prevmag<mag) && (nextmag<mag) && (mag>ampcheck) && numpeaks < maxPeaks
      ) {
  			// found a peak
  			peaks[numpeaks].mag = mag * constants.ampmult;
        double freq = getTrueFreq(i-1, phase - data->prevPhases[i-1]);
  			peaks[numpeaks].freq = freq;
  			peaks[numpeaks].bin = i-1;
  			if(peaks[numpeaks].mag > maxMag){
  				maxMag = peaks[numpeaks].mag;
  				maxFreq = freq;
  			}
  			// Print("%f ",peaks[numpeaks].mag);
  			data->peakMul *= peaks[numpeaks].mag;
  			++numpeaks;
  		}

  		data->prevPhases[i-1] = phase;
  		prevmag=mag;
  		mag=nextmag;
  	}
  	// Print("\n-----\n");
    // Print("%d peaks, max: %fHz %f tot: %f\n",numpeaks,maxFreq,maxMag, data->peakMul);

  	data->maxMag = maxMag;
  	data->maxFreq = maxFreq;
  	data->numpeaks = numpeaks;
  }



};

}
