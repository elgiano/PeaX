PitchShiftFollow : UGen {
  *new{ |source, template, clarThr=0.9, ampLag=0.1, pitchTracker=\Pitch|
      ^this.multiNew(source, template, clarThr, ampLag, pitchTracker);
  }

  *new1{ |source, template, clarThr=0.9, ampLag=0.1, pitchTracker=\Pitch|
      var sourceAmp = Amplitude.kr(source).lag(ampLag);
      var pitch = pitchTracker.asClass.kr([source,template], clar:1);
      pitch = pitch.collect{|p|Latch.kr(p[0],p[1]>clarThr)};

      ^Normalizer.ar(
          PitchShift.ar(template,0.1,
              pitch[0]/pitch[1],0.01,0.01
          )
		  )*sourceAmp;
  }

}
