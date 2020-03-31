PV_TPV : UGen {
	*new { arg buffer, hop, maxPeaks=80, numPeaks=80, tolerance=1, ampThr=0;
		^this.multiNew('control', buffer, hop, maxPeaks, numPeaks, tolerance, ampThr)
	}

	checkInputs { ^this.checkValidInputs }

}
