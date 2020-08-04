# k-mer Pattern Partition Prediction

A rust library for reading and querying point mutation probabilities.

For example, when you have an input file like this:
```
A->C MNAANVR 7.589541872637903e-06
A->C RNKAVKN 7.839974716364649e-06
A->G NAAATNN 2.6344844835403154e-05
A->G NRTASWD 3.8916665703091594e-05
A->T NNAATNV 5.161581369396004e-06
A->T NNHACNN 6.636144968007699e-06
C->A ANACAAN 1.9414406870708827e-05
C->A BBACVAN 1.783626464296544e-05
C->A HNRCWKN 1.1882887322469946e-05
C->G NADCYDD 1.9273637543165863e-05
C->G NCHCKNC 2.0450439690457853e-05
C->G NDNCGNN 2.3540691630409924e-05
C->G NNWCAAH 1.3384209291990454e-05
C->T DNSCYTN 5.0466135573024847e-05
C->T NCACGVN 0.0006964474916140189
C->T NHACMRA 3.95017705645649e-05
C->T NVSCGNW 0.0006052619765592556
C->A NNRCGKN 3.8108173013492654e-05
```

The first row tells you that a point mutation from an Adenine to Cytosine with
the [IUPAC](https://www.bioinformatics.org/sms/iupac.html) context `MNA` (left)
and `NVR` (right) has a probability of 7.589541872637903e-06.

This library makes it more conventient to access this information from a
papapred file (the file whose's contents are shown above).
