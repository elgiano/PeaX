HPS : UGen {
	*new { arg buffer, oversampling=0.5, harmonics=8, minfreq=50, maxfreq=10000, ampThr=0.01, correction=1;
		^this.multiNew('control', buffer, oversampling, harmonics, minfreq, maxfreq, ampThr, correction)
	}
	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}
