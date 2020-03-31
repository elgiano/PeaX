// PluginTrax.cpp
// Gianluca Elia (elgiano@gmail.com)

#include "SC_PlugIn.hpp"
#include "FFT_UGens.h"
#include "PV_ShiftConv.hpp"

namespace Trax {

PV_ShiftConv::PV_ShiftConv(): m_tempbuf(0) {
    out0(0) = in0(0);
    mCalcFunc = make_calc_function<PV_ShiftConv, &PV_ShiftConv::next>();
    next(1);
}

void PV_ShiftConv::next(int nSamples) {
  auto* unit = this;
  PV_GET_BUF2

  if (!m_tempbuf) {
      m_tempbuf = (float*) RTAlloc(mWorld, buf1->samples * sizeof(float));
      m_binCounter = (float*) RTAlloc(mWorld, numbins * sizeof(float));
      m_numbins = numbins;
      m_ampmult = (double) 2.f/(numbins+1); //*(1.0/unit->m_maxpeaks);
      m_rAmpmult = (double) (numbins+1)/2.f; //*(1.0/unit->m_maxpeaks);
  } else if (numbins != m_numbins)
      return;

  // Controls
  bool stretch = (bool)in0(2);
  bool complex = (bool)in0(3);

  memset(m_binCounter, 1, numbins*sizeof(float));

  if(stretch){
    combineFunc = &PV_ShiftConv::shiftStretch;
  }else{
    if(complex)
      combineFunc = &PV_ShiftConv::shift_complex;
    else
      combineFunc = &PV_ShiftConv::shift;
  };


  fromBinA = (int)(in0(4) * buf1->samples / mWorld->mFullRate.mSampleRate);
  toBinA = (int)(in0(5) * buf1->samples / mWorld->mFullRate.mSampleRate);
  fromBinB = (int)(in0(6) * buf1->samples / mWorld->mFullRate.mSampleRate);
  toBinB = (int)(in0(7) * buf1->samples / mWorld->mFullRate.mSampleRate);
  if(fromBinA<=0 || fromBinA >=numbins ) fromBinA = 0;
  if(fromBinB<=0 || fromBinB >=numbins ) fromBinB = 0;
  if(toBinA<=0 || toBinA >= numbins) toBinA = numbins-1;
  if(toBinB<=0 || toBinB >= numbins) toBinB = numbins-1;
  thr = in0(8);

  int doAverage = (int) in0(9);

  f0Bin = (int) (in0(10) * (float) buf1->samples / mWorld->mFullRate.mSampleRate);
  if(f0Bin<=0) f0Bin = 1;
  if(f0Bin>=m_numbins) f0Bin = m_numbins-1;
  interp = (int) in0(11);

  (this->*combineFunc)(buf1,buf2);

  SCPolarBuf* t = (SCPolarBuf*) m_tempbuf;

  if(doAverage){
    if(complex){
      SCComplexBuf* complex = (SCComplexBuf*) m_tempbuf;
      for(int i=0; i < numbins; ++i){
        complex->bin[i].real /= m_binCounter[i];
        complex->bin[i].imag /= m_binCounter[i];
      }
    }else{
      for(int i=0; i < numbins; ++i){
        t->bin[i].mag /= m_binCounter[i];
      }
    }
  }


  memcpy(ToPolarApx(buf1)->bin, t->bin, m_numbins * sizeof(SCPolar));

}

void PV_ShiftConv::shift(SndBuf* buf1, SndBuf* buf2){

  SCPolar* a = ToPolarApx(buf1)->bin;
  SCPolar* b = ToPolarApx(buf2)->bin;
  SCPolar* t = ((SCPolarBuf*) m_tempbuf)->bin;

  // zero tmp buf
  memset(t,0, m_numbins * sizeof(SCPolar));

  int numCopies, maxB, i, j;
  float aMag, bMag;

  for(i=fromBinA; i < toBinA ; ++i){
    t[i].phase = a[i].phase;
    numCopies = sc_min(i - fromBinA, toBinB - fromBinB) + 1;
    maxB = sc_min(fromBinB+numCopies-1, toBinB);
    for(j=0; j < numCopies; ++j){
      aMag = a[i-j].mag* m_ampmult;
      bMag = b[fromBinB+j].mag* m_ampmult;
      if( aMag < thr || bMag < thr) continue;
      t[i].mag += aMag * bMag;
      m_binCounter[i]++;
    };
    t[i].mag *= m_rAmpmult;
  }

}

void PV_ShiftConv::shift_complex(SndBuf* buf1, SndBuf* buf2){

  SCComplex* a = ToComplexApx(buf1)->bin;
  SCComplex* b = ToComplexApx(buf2)->bin;
  SCComplex* t = ((SCComplexBuf*) m_tempbuf)->bin;

  int numCopies, maxB, i, j;

  // zero tmp buf
  memset(t,0,m_numbins*sizeof(SCComplex));

  for(i=fromBinA; i < toBinA ; ++i){
    if(a[i].ToPolarApx().mag < thr) continue;
    numCopies = i - fromBinA;
    maxB = sc_min(fromBinB+i, toBinB);
    for(j=0; j <= numCopies, maxB-j >= 0; ++j){
      if(b[maxB-j].ToPolarApx().mag < thr) continue;
      t[i] += a[fromBinA + j] * b[maxB-j];
      m_binCounter[i]++;
    };
  }

}

void PV_ShiftConv::shiftStretch(SndBuf* buf1, SndBuf* buf2){

  SCPolar* a = ToPolarApx(buf1)->bin;
  SCPolar* b = ToPolarApx(buf2)->bin;
  SCPolar* t = ((SCPolarBuf*) m_tempbuf)->bin;

  int i, j, sourceBin, shift;
  double stretch, fpos;
  float sourceMag;

  // init tmp buf: 0 magnitudes and copy phases
  /*for(i=0; i<m_numbins; ++i) {
    t[i].mag = 0; t[i].phase = a[i].phase;
  }*/
  // zero tmp buf
  memset(t,0, m_numbins * sizeof(SCPolar));

  // shift stretch sums
  for( sourceBin=fromBinA; sourceBin < toBinA ; ++sourceBin ){
    sourceMag = a[sourceBin].mag * m_ampmult;
    if(sourceMag < thr) continue;
    shift = sourceBin - f0Bin;
    stretch = (double) sourceBin / (double) f0Bin;
    //Print("amag %d %f f0 %d\n",sourceBin, sourceMag, f0Bin);
    //Print("sh %d f0 %d st %f\n",shift,f0Bin,stretch);

    if (interp <= 0) {
      for (i = fromBinB/*, fpos = i + shift*/; i < toBinB; ++i/*, fpos += stretch*/) {
          // if(i < fromBinB) continue;
          float bMag = b[i].mag * m_ampmult;
          //Print("bmag %f\n", bMag);
          if(bMag < thr) continue;
          int binDiff = (i - f0Bin);
          fpos = (binDiff*stretch) + sourceBin;
          int32 pos = (int32)(fpos + 0.5);
          //Print("bin %d f0 %d fpos %f sh %d st %f\n",i,f0Bin,fpos,sourceBin,stretch);
          if (pos >= 0 && pos < m_numbins) {
              t[pos].mag += sourceMag * bMag;
              m_binCounter[pos]++;
          }
      }
    } else {
      for (i = fromBinB/*, fpos = i + shift*/; i < toBinB; ++i/*, fpos += stretch*/) {
          //if(i<fromBinB) continue;
          float bMag = b[i].mag * m_ampmult;
          if(bMag < thr) continue;
          int binDiff = (i - f0Bin);
          fpos = (binDiff*stretch) + sourceBin ;
          //Print("bmag %f\n", bMag);
          //Print("bin %d f0 %d fpos %f sh %d st %f\n",i,f0Bin,fpos,sourceBin,stretch);
          int32 fpos0 = (int32)std::floor(fpos);
          int32 fpos1 = fpos0 + 1;
          float beta = fpos - std::floor(fpos);
          float alpha = 1.0f - beta;
          if (fpos0 >= 0 && fpos0 < m_numbins) {
              t[fpos0].mag += sourceMag * alpha * b[i].mag;
              m_binCounter[fpos0]++;
          }
          if (fpos1 >= 0 && fpos1 < m_numbins) {
              t[fpos1].mag += sourceMag * beta * b[i].mag;
              m_binCounter[fpos1]++;
          }
      }
    }


  }
  for(i=0; i < m_numbins; i++){
    if(m_binCounter[i]>0) t[i].mag *= m_rAmpmult;
    t[i].phase = a[i].phase;
  }

}

// TODO: COPY FIXES FROM shiftStretch
void PV_ShiftConv::shiftStretch_complex(SndBuf *buf1, SndBuf *buf2){

  SCComplex* a = ToComplexApx(buf1)->bin;
  SCComplex* b = ToComplexApx(buf2)->bin;
  SCComplex* t = ((SCComplexBuf*) m_tempbuf)->bin;

  int i, j, sourceBin, destBin, coeffBin;
  float shift, stretch, sourceMag, fpos;

  int f0Bin = (int) (in0(10) * buf1->samples / mWorld->mFullRate.mSampleRate);
  if(f0Bin<=0) f0Bin = 1;
  if(f0Bin>=m_numbins) f0Bin = m_numbins-1;
  int interp = (int) in0(11);

  // zero tmp buf
  memset(t,0,m_numbins*sizeof(SCComplex));

  // shift stretch complex sums
  for( sourceBin=fromBinA; sourceBin < toBinA ; ++sourceBin ){
      shift = sourceBin - f0Bin;
      stretch = sourceBin/f0Bin;

      if (interp > 0) {
            for (i = 0, fpos = shift; i < toBinB; ++i, fpos += stretch) {
                if(i<fromBinB) continue;
                int32 fpos0 = (int32)std::floor(fpos);
                int32 fpos1 = fpos0 + 1;
                float beta = fpos - std::floor(fpos);
                float alpha = 1.0f - beta;
                if (fpos0 >= 0 && fpos0 < m_numbins) {
                    t[fpos0] += a[sourceBin] * alpha * b[i];
                }
                if (fpos1 >= 0 && fpos1 < m_numbins) {
                    t[fpos1] += a[sourceBin] * beta * b[i];
                }
            }
        } else {
            for (i = 0, fpos = shift; i < toBinB; ++i, fpos += stretch) {
                if(i<fromBinB) continue;
                int32 pos = (int32)(fpos + 0.5);
                if (pos >= 0 && pos < m_numbins) {
                    t[pos] += a[sourceBin] *b[i];
                }
            }
        }

  }

}

} // namespace Trax

InterfaceTable* ft;

PluginLoad(TraxUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<Trax::PV_ShiftConv>(ft, "PV_ShiftConv");
}
