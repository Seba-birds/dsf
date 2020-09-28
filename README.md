\mainpage Puredata DSF external

\section intro Discrete summation formula

When synthesizing sounds digitally,
classic waveforms with rich spectra
inevitably create unwanted alias-frequencies.
This results from the fact that the overtone
series goes to infinity while a digital system
can only realize a limited bandwith. The
surplus energy is being reflected back into
the created frequency spectrum as alias-frequencies.

One way around this problem is to create the
overtones piece by piece up to the maximal frequency
(which is the Nyquist frequency or half the sample
rate at which the digital system operates).






