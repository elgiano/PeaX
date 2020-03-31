// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "FFT_UGens.h"
#include "PV_TPV.hpp"

namespace Trax {

///////////////////////////////////////////////////

double PV_TPV::getTrueFreq(int bin_i, float phaseDiff){
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

double PV_TPV::getPhaseCorrection(int bin_i, double freq){
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

PV_TPV::PV_TPV(): m_tempbuf(0) {
    out0(0) = in0(0);
    mCalcFunc = make_calc_function<PV_TPV, &PV_TPV::next>();
    next(1);
}

void PV_TPV::init(SndBuf* buf, int numbins){
  m_numbins = numbins;
  m_ampmult = (double) 2.f/(numbins+1); //*(1.0/unit->m_maxpeaks);
  m_rAmpmult = (double) (numbins+1)/2.f; //*(1.0/unit->m_maxpeaks);
  m_windowsize = buf->samples;

  m_rOversampling = in0(1);
  m_maxpeaks = in0(2);
  m_hopsize = m_rOversampling * m_windowsize;

  m_freqPerBin = (double) mWorld->mFullRate.mSampleRate / m_windowsize;
  m_rFreqPerBin = (double) m_windowsize / mWorld->mFullRate.mSampleRate;
  m_expct = twopi * (double) m_rOversampling;

  m_tempbuf = (float*)RTAlloc(mWorld, buf->samples * sizeof(float));
  m_prevPhases = (float*)RTAlloc(mWorld, numbins * sizeof(float));
  m_phaseAccumulator = (float*)RTAlloc(mWorld, numbins * sizeof(float));
  data = new PV_TPV_Data(mWorld, m_maxpeaks, numbins);

  memset(m_prevPhases, 0, numbins * sizeof(float));
	memset(m_phaseAccumulator, 0, numbins * sizeof(float));

	if(mWorld->mFullRate.mBufLength!=64)
    Print("PV_TPV complains: block size not 64, you have %d\n", mWorld->mFullRate.mBufLength);
}

PV_TPV::~PV_TPV(){
  if(data) delete data;
  if(m_tempbuf) RTFree(mWorld, m_tempbuf);
  if(m_prevPhases) RTFree(mWorld, m_prevPhases);
  if(m_phaseAccumulator) RTFree(mWorld, m_phaseAccumulator);
}

void PV_TPV::next(int nSamples) {
  auto* unit = this;
  PV_GET_BUF

  if (!m_tempbuf)
      init(buf,numbins);
  else if (numbins != m_numbins)
      return;

  m_numpeaksrequested = in0(3);
  m_tolerance = in0(4);
  data->ampThr = in0(5);

  SCPolarBuf* fft = ToPolarApx(buf);
  SCPolarBuf* t = (SCPolarBuf*) m_tempbuf;

  // Print("findPeaks\n");
  findPeaks(data, fft);
  // Print("calcTracks\n");
  calcTracks(data);
  // Print("writeTracks\n");
  writeTracks(data->numtracks, data->tracks, t);


  memcpy(fft, t, m_windowsize * sizeof(float));
}

void PV_TPV::writeTracks(int numtracks, PartialTrack* tracks, SCPolarBuf* fft){
  memset(fft->bin, 0, m_numbins * sizeof(SCPolar));
  for (int i=0; i<numtracks; ++i) {
    PartialTrack *track = &tracks[i];
		int index = sc_floor(track->omega2 * m_rFreqPerBin);
    if(index < 0 || index > m_numbins){
      Print("BIN out of range\n");
      continue;
    }
    double phase = getPhaseCorrection(index, track->omega2);
		fft->bin[index].phase = phase;
		fft->bin[index].mag += track->amp2;
  }

  for(int i= 0; i<m_numbins; ++i)
    fft->bin[i].mag *= m_rAmpmult;

}

void PV_TPV::findPeaks(PV_TPV_Data* data, SCPolarBuf* fft){

	//swap new peaks to old; current now safe to overwrite;
	TPVPeak *prevpeaks = data->prevpeaks;
	TPVPeak *newpeaks = data->newpeaks;
	int numprevpeaks = data->numprevpeaks;
	int numnewpeaks = data->numnewpeaks;
	//Print("prev pointer %p new pointer %p \n",prevpeaks, newpeaks);
	//ditch old
	numprevpeaks = numnewpeaks;
	numnewpeaks = 0;
	//swap pointers ready to write new peaks
	TPVPeak *temp = prevpeaks;
	prevpeaks = newpeaks;
	newpeaks = temp;
	//Print("prev pointer %p new pointer %p temp %p \n",prevpeaks, newpeaks, temp);

  int maxpeaks = sc_min(data->maxpeaks,m_numpeaksrequested);

  // angular frequency is pi*(i/nover2)
  // double angmult = pi/(m_numbins+1);
  float ampcheck = data->ampThr * m_rAmpmult;

	float phase, prevmag, mag, nextmag;

	//bin 1 can't be pick since no interpolation possible! dc should be ignored
	//test each if peak candidate; if so, add to list and add to peaks total
	prevmag = fft->bin[0].mag; //this is at analysis frequency, not dc
	mag = fft->bin[1].mag;

  // float minPhase, maxPhase, minMag, maxMag;

	//could restrict not to go above nover4 (numbins/2)!
	for (int i=1; i< m_numbins; ++i) {

		phase= fft->bin[i].phase;
		nextmag= fft->bin[i].mag;

		if (
      (prevmag<mag) && (nextmag<mag) &&
      (mag>ampcheck) && (numnewpeaks<maxpeaks)
    ) {
			// found a peak
			// could use cubic interpolation/successive parabolic interpolation to refine peak location; or should have zero padded
			newpeaks[numnewpeaks].mag = mag * m_ampmult; //must divide by fftsize before resynthesis!
      float phaseDev = phase - data->prevPhases[i];
      double freq = getTrueFreq(i, phaseDev);
			newpeaks[numnewpeaks].freq = freq;
			newpeaks[numnewpeaks].phase = phase;	//is this in range -pi to pi? more like -1 to 5 or so, but hey, is in radians
			// Print("newpeak %d bin: %d p: %f pp: %f f: %f\n",numnewpeaks, i , phase, data->prevPhases[i], newpeaks[numnewpeaks].freq, fft->bin[i-1].phase);
      /*if(numnewpeaks == 0 ){
        minPhase = phase; maxPhase = phase;
        minMag = mag*m_ampmult; maxMag = mag*m_ampmult;
      }else{
        maxPhase = sc_max(maxPhase,phase);
        minPhase = sc_min(minPhase,phase);
        maxMag = sc_max(maxMag,mag*m_ampmult);
        minMag = sc_min(minMag,mag*m_ampmult);
      }*/
			++numnewpeaks;
		}

		data->prevPhases[i] = phase;
		prevmag=mag;
    Print("mag %f\n",mag);
		mag=nextmag;
	}
  // Print("phis: %f %f\n",minPhase,maxPhase);
  // Print("mags: %f %f\n",minMag,maxMag);

	data->prevpeaks = prevpeaks;
	data->newpeaks = newpeaks;
	data->numprevpeaks = numprevpeaks;
	data->numnewpeaks = numnewpeaks;
}

void PV_TPV::calcTracks(PV_TPV_Data* data){
  //now peak matching algorithm
  //int leftsort=0;
  int rightsort=0;
  bool flag= true;
  //float rightfreq= newpeaks[0].freq;

  PartialTrack * tracks = data->tracks;
  int numtracks = 0; //unit->m_numtracks;
  //increase tolerance
  float tolerance= m_tolerance;//*angmult;

  float testfreq;

  float T = m_hopsize;
  int i,j;
  TPVPeak *prevpeaks = data->prevpeaks;
  TPVPeak *newpeaks = data->newpeaks;
  int numprevpeaks = data->numprevpeaks;
  int numnewpeaks = data->numnewpeaks;

  //	SEEMS OK
  //	Print("numprevpeaks %d numnewpeaks %d \n",numprevpeaks, numnewpeaks);
  //	//print all and look for junk data
  //	for (i=0; i<numnewpeaks; ++i)
  //	Print("new i %d amp %f freq %f phase %f \n",i,newpeaks[i].mag,newpeaks[i].freq,newpeaks[i].phase);
  //
  //	for (i=0; i<numprevpeaks; ++i)
  //	Print("prev i %d amp %f freq %f phase %f \n",i,prevpeaks[i].mag,prevpeaks[i].freq,prevpeaks[i].phase);
  //
  //

	// ASSUMES BOTH PEAKS LISTS ARE IN ORDER OF INCREASING FREQUENCY
	// while right less than left-tolerance then birth on right
	// if right within tolerance, find closest; if less than, match, else must check next on left whether better match. If not, match, else, check previous on right. If within tolerance, match, else death on right.

	//step through prevpeaks
	for (i=0; i<numprevpeaks; ++i) {
		float freqnow= prevpeaks[i].freq;

		flag=true;
		while(flag) {

			if(rightsort>=numnewpeaks) {flag=false;} else {
				testfreq= newpeaks[rightsort].freq;

				if((testfreq*tolerance)<freqnow) {
					//birth on right
					tracks[numtracks].omega1=newpeaks[rightsort].freq;
          tracks[numtracks].omega2=newpeaks[rightsort].freq; //match to itself
          tracks[numtracks].theta1=newpeaks[rightsort].phase - (T*(newpeaks[rightsort].freq)); //should really be current phase + freq*hopsize
					tracks[numtracks].theta2=newpeaks[rightsort].phase;
					tracks[numtracks].amp1=0.0;
					tracks[numtracks].amp2=newpeaks[rightsort].mag;
					++numtracks;
					++rightsort;

				} else flag=false;
			}
		}

		flag=false; //whether match process fails
		if(rightsort>=numnewpeaks) {flag=true;} else {
			//Print("testfreq %f freqnow %f tolerance %f \n ", testfreq, freqnow, tolerance);
			//assumption that testfreq already valid;
			if (testfreq>(freqnow*tolerance)) {flag=true;} else {

				//now have a candidate. search for closest amongst remaining; as soon as no closer, break
				//Print("candidate! \n ");

				float bestsofar = fabs(freqnow- testfreq);
				int bestindex = rightsort;

				for (j=(rightsort+1); j<numnewpeaks; ++j) {
					float newcandidate= newpeaks[j].freq;
					float newproximity= fabs(newcandidate-freqnow);
					//must keep getting closer, else no use
					if(newproximity<bestsofar) {bestindex= j; bestsofar= newproximity;}
					else break; //nothing better
				}
				//now have closest estimate. If less than freqnow have match
				float closest= newpeaks[bestindex].freq;
				bool havematch=false;
				//Print("closest! %f bestindex %d rightsort %d \n ", closest, bestindex, rightsort);

				if(closest<freqnow || (i==(numprevpeaks-1))) havematch=true;
				else { //test next i as available in this case

					float competitor = prevpeaks[i+1].freq;
					if (fabs(competitor-closest)<bestsofar) {
						//if no alternative
						if (bestindex==rightsort) flag= true; //failure to match anything
						else {bestindex= rightsort-1;
							havematch=true;
						}
					} else
					havematch=true;
				}

				if(havematch) {
					//int newrightsort= bestindex;
					//if() newrightsort=

					//TIDY UP ANY CANIDATES MISSED OUT BY THIS PROCESS
					for (j=rightsort; j<=(bestindex-1);++j) {
						//BIRTHS ON RIGHT
						tracks[numtracks].omega1=newpeaks[j].freq;
						tracks[numtracks].theta2=newpeaks[j].phase;
						tracks[numtracks].omega2=newpeaks[j].freq; //match to itself
						tracks[numtracks].theta1=sc_wrap(newpeaks[j].phase - (T*(newpeaks[j].freq)),0.0f,(float)twopi); //backcalculate starting phase
						tracks[numtracks].amp1=0.0;
						tracks[numtracks].amp2=newpeaks[j].mag;
						++numtracks;
						++rightsort;
					}

					//Print("match! \n ");
					//MATCH!
					tracks[numtracks].theta1=prevpeaks[i].phase;
					tracks[numtracks].omega1=prevpeaks[i].freq;
					tracks[numtracks].theta2=newpeaks[rightsort].phase; //match to itself; should really be current phase + freq*hopsize
					tracks[numtracks].omega2=newpeaks[rightsort].freq; //match to itself
					tracks[numtracks].amp1=prevpeaks[i].mag;
					tracks[numtracks].amp2=newpeaks[rightsort].mag;

					//yes, OK
					//Print("amp check i %d amp1 %f amp2 %f source1 %f source2 %f\n",i,tracks[numtracks].amp1, tracks[numtracks].amp2, prevpeaks[i].mag, newpeaks[rightsort].mag);
					++numtracks;
					++rightsort;
					//rightsort=bestindex+1;
				}
				//if was flag==true, then none missed out, still on rightsort
			}
		}

		//match failed, death on left
		if (flag==true) {
			//DEATH ON LEFT
			//death on left
			tracks[numtracks].omega1=prevpeaks[i].freq;
      tracks[numtracks].omega2=prevpeaks[i].freq; //match to itself
      tracks[numtracks].theta1=prevpeaks[i].phase;
			tracks[numtracks].theta2=sc_wrap(prevpeaks[i].phase + (T*prevpeaks[i].freq),0.0f,(float)twopi); //match to itself; should really be current phase + freq*hopsize
			tracks[numtracks].amp1=prevpeaks[i].mag;
			tracks[numtracks].amp2=0.0;
			++numtracks;
			//ADDCODE
			//++leftsort;
		}
	}
	//rightsort should equal numnewpeaks!
	data->numtracks = numtracks;
}

} // namespace Trax

InterfaceTable* ft;

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::PV_TPV>(ft, "PV_TPV");
}
