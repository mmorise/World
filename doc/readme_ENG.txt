WORLD - a high-quality speech analysis, manipulation and synthesis system

WORLD is free software for high-quality speech analysis, manipulation and synthesis.
It can estimate Fundamental frequency (F0), aperiodicity and spectral envelope
and also generate the speech like input speech with only estimated parameters.

2. Usage
Please see test.cpp.

(1) F0 estimation by Dio()
(1-1) F0 is refined by StoneMask() if you need more accurate result.
(2) Spectral envelope estimation by CheapTrick()
(3) Aperiodicity estimation by D4C().
(4) You can manipulation these parameters in this phase. 
(5) Voice synthesis by Synthesis().

English document is written by a Japanese poor editor.
I willingly accept your kind indication on my English text.

3. License
WORLD is free software, and you can redistribute it and 
modify it under the terms of the modified BSD License.
Please see copying.txt for more information.
You can use this program for business, while I hope that 
you contact me after you developed software with WORLD.
This information is crucial to obtain a grant to develop WORLD.

4. Contacts
WORLD was written by Masanori Morise.
You can contact him via e-mail (mmorise [at] yamanashi.ac.jp)
or Twitter: @m_morise