# WORLD - a high-quality speech analysis, manipulation and synthesis system

WORLD is free software for high-quality speech analysis, manipulation and synthesis.
It can estimate Fundamental frequency (F0), aperiodicity and spectral envelope and also generate the speech like input speech with only estimated parameters.

This source code is released under the modified-BSD license.
There is no patent in all algorithms in WORLD.

## Important notice (2018/01/21)
FFT functions were modified.
If you have used another FFT function, you must change the sign of imaginary part in the spectrum for the compatibility.
This version outputs the same result compared with major FFT functions (e.g. FFTW and MATLAB).

## References
When you cite the latest version of WORLD in your paper, please use the sentence "WORLD \[1\] (D4C edition [2])" and cite the following papers.  
[1] M. Morise, F. Yokomori, and K. Ozawa: WORLD: a vocoder-based high-quality speech synthesis system for real-time applications, IEICE transactions on information and systems, vol. E99-D, no. 7, pp. 1877-1884, 2016.  
[2] M. Morise: D4C, a band-aperiodicity estimator for high-quality speech synthesis, Speech Communication, vol. 84, pp. 57-65, Nov. 2016. http://www.sciencedirect.com/science/article/pii/S0167639316300413  

In CheapTrick, you can refer the following references.  
[3] M. Morise: CheapTrick, a spectral envelope estimator for high-quality speech synthesis, Speech Communication, vol. 67, pp. 1-7, March 2015. http://www.sciencedirect.com/science/article/pii/S0167639314000697  
[4] M. Morise: Error evaluation of an F0-adaptive spectral envelope estimator in robustness against the additive noise and F0 error, IEICE transactions on information and systems, vol. E98-D, no. 7, pp. 1405-1408, July 2015.  

In DIO, you can refer the following reference.  
[5] M. Morise, H. Kawahara and H. Katayose: Fast and reliable F0 estimation method based on the period extraction of vocal fold vibration of singing voice and speech, AES 35th International Conference, CD-ROM Proceeding, Feb. 2009.

In Harvest, you can refer the following reference.  
[6] M. Morise: Harvest: A high-performance fundamental frequency estimator from speech signals, in Proc. INTERSPEECH 2017, pp. 2321–2325, 2017. http://www.isca-speech.org/archive/Interspeech_2017/abstracts/0068.html

In the codec of spectral envelope, you can refer the following reference.  
[7] M. Morise, G. Miyashita and K. Ozawa: Low-dimensional representation of spectral envelope without deterioration for full-band speech analysis/synthesis system, in Proc. INTERSPEECH 2017, pp. 409-413, 2017. http://www.isca-speech.org/archive/Interspeech_2017/abstracts/0067.html

A paper was published to demonstrate that the current version of WORLD was superior to the similar vocoders in the sound quality of re-synthesized speech. This paper also includes the detailed information in the D4C LoveTrain used in the latest version.  
[8] M. Morise and Y. Watanabe: Sound quality comparison among high-quality vocoders by using re-synthesized speech, Acoust. Sci. & Tech., vol. 39, no. 3, pp. 263-265, May 2018. https://www.jstage.jst.go.jp/article/ast/39/3/39_E1779/_article/-char/en

## Introduction of WORLD family

I introduce useful software in WORLD. If you want to introduce your project in WORLD, please contact me.

PyWorldVocoder (https://github.com/JeremyCCHsu/Python-Wrapper-for-World-Vocoder) is a Python wrapper for World Vocoder.

Python-WORLD (https://github.com/tuanad121/Python-WORLD) is line-by-line implementation of WORLD vocoder (Matlab, C++) in python.

world-class (https://github.com/yukara-ikemiya/world-class) is a C++ library of WORLD.

Note: To avoid making the project complicated, I decided not to merge it to my repository and introduce your project here. The other reason is that I can't support some computer languages.
