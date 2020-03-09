import numpy as np
import matplotlib.pyplot as plt

# Data
Df=75e3 # Frequency deviation relative to the carrier frequency - 75kHz - ENACOM
fp=19e3 # 19kHz - ENACOM
fc=fp*30.0 # Frecuencia de la portadora
# fc=0
Ap=10.0 # 10v
Ac=1.0
fs=fp*100.0 # fs=kHz    Ts=5.263us
Df_Sp=Df*0.1 # 10% of total deviation
h=Df_Sp/fp # Modulation index - h=0.4
K=Df_Sp/Ap # [Hz/V]

print("Frecuencia de portadora: {}Hz".format(fc))
print(h)
# Signal creation
t=np.arange(0,160000)/fs
Sp=Ap*np.cos(2*np.pi*fp*t)
FFTsp=np.fft.fft(Sp)/len(Sp)
FFTspdB=20*np.log10(np.abs(FFTsp))
f=np.fft.fftfreq(len(Sp),1/fs)

# PLOT
# plt.plot(t,Sp)
# plt.grid()
# plt.show()
# plt.plot(f,FFTspdB)
# plt.grid()
# plt.show()

IntSp=Ap/(2*np.pi*fp)*np.cos(2*np.pi*fp*t)
# plt.plot(t,IntSp)
# plt.grid()
# plt.show()
SFM=Ac*np.cos(2*np.pi*fc*t+K*2*np.pi*IntSp)
FFTsfm=np.fft.fft(SFM)/len(SFM)
FFTsfmdB=20.0*np.log10(np.abs(FFTsfm))
# plt.plot(t,SFM)
# plt.grid()
# plt.show()
plt.subplot(2,1,1)
plt.plot(f,FFTsfmdB)
plt.grid()
plt.title("Espectro FM con fórmula")
# plt.show()

# Bessel coefficients

J=[0.9604,0.196,0.0197,0.0197,0.196] # Considering harmonics of relative amplitude 0.01 (>-40dB)
SFMB=np.zeros(len(t))
for i in range(-2,3):
    SFMB+=J[i]*np.cos(2*np.pi*fc*t+i*2*np.pi*fp*t+i*np.pi/2.0)
FFTsfmb=np.fft.fft(SFMB)/len(SFMB)
FFTsfmbdB=20.0*np.log10(np.abs(FFTsfmb))
f=np.fft.fftfreq(len(SFMB),1/fs)
# plt.plot(t,SFMB)
# plt.grid()
# plt.show()
plt.subplot(2,1,2)
plt.plot(f,FFTsfmbdB)
plt.grid()
plt.title("Espectro FM con Bessel")
plt.show()

for i in range(len(f)):
    if (f[i]%19e3)==0 and f[i]<(fc+3*fp) and f[i]>(fc-3*fp):
        print("SFM (f={:.3f})={:.3f}dB\t{:.3}\tfc  {:2}*fp\t{:.3f}dB".format(f[i],FFTsfmdB[i],FFTsfm[i],int((f[i]-fc)/fp),(FFTsfmdB[i]-FFTsfmdB[i-1600])))
        print("SFMB(f={:.3f})={:.3f}dB\t{:.3}\tfc  {:2}*fp\t{:.3f}dB".format(f[i],FFTsfmbdB[i],FFTsfmb[i],int((f[i]-fc)/fp),(FFTsfmbdB[i]-FFTsfmbdB[i-1600])))