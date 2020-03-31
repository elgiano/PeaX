PeaXShift2 : UGen {
	*new { arg bufferA, bufferB, hop, maxPeaks=80, rootA=440, rootB=440, tolerance=1, ampThr=0;
		^this.multiNew('control', bufferA, bufferB, hop, maxPeaks, rootA, rootB, tolerance, ampThr)
	}

	checkInputs { ^this.checkValidInputs }

}
