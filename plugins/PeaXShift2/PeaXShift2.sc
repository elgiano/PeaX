PeaXShift2 : UGen {
	*new { arg bufferA, bufferB, hop, rootFreq=440, ampThr=0, maxPeaks=(-1);
		^this.multiNew('control', bufferA, bufferB, hop, rootFreq, ampThr, maxPeaks)
	}

	checkInputs { ^this.checkValidInputs }

}
