PV_ShiftConv : PV_ChainUGen {
	*new { arg bufferA, bufferB, minFreqA=0, maxFreqA=(-1), minFreqB=0, maxFreqB=(-1), thr=0, complex=0, average=0, stretch=0, f0=440, interp=1;
		^this.multiNew('control', bufferA, bufferB,
			stretch, complex,
		  minFreqA, maxFreqA, minFreqB, maxFreqB, thr, average,
			f0, interp
		)
	}

	*stretch { arg bufferA, bufferB, minFreqA=0, maxFreqA=(-1), minFreqB=0, maxFreqB=(-1), thr=0, f0=440, interp=1, complex=0, average=0;
		^this.multiNew('control', bufferA, bufferB,
	  	1, complex,
			minFreqA, maxFreqA, minFreqB, maxFreqB, thr, average,
		  f0, interp
	)}

	checkInputs { ^this.checkValidInputs }

}


PV_MagShiftStretchConv : PV_ChainUGen {
	*new { arg bufferA, bufferB, minFreqA=0, maxFreqA=(-1), minFreqB=0, maxFreqB=(-1), thr=0, f0=440, interp=1, complex=0, average=0;
		^PV_ShiftConv.stretch(
			bufferA, bufferB, minFreqA, maxFreqA, minFreqB, maxFreqB, thr, f0, interp, complex, average
		)
	}
}
//PV_ShiftStretchConv : PV_MagShiftStretchConv {}
