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
Especially with low frequencies, this is not
feasible, because the number of sine waves to
keep track of can grow up to 1000 for low
bass notes.

DSF synthesis creates sine and cosine waves
by rotating phasors. Instead of generating
sine waves up to the Nyquist frequency,
the complex numbers which represent the 
position of the phasors are being put into
a geometric series for which a closed
formula exists. This formula can be computed
in near-linear time, so that computational
costs stay low.

One main factor in this is the efficient
calculation of a complex power function
as implemented here: power_complex().
Using a "divide-and-conquer" approach,
the number of calculations are reduced
to the 2-logarithm of the number of calculations
in the naive implementation power_complex_naiv().



