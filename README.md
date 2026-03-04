# ISMRMRD.jl

Fast Julia reader for the [ISMRMRD](https://ismrmrd.readthedocs.io/en/latest/) MRI raw data format.

Compared to the reader in [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), this package reads the entire acquisition dataset in a **single HDF5 call** (O(1) instead of O(N)), deserialises headers via direct byte-offset access (no `IOBuffer`, no reflection), and constructs profiles in parallel with `Threads.@threads`.

Non-registered so add with url. After some more testing it could be submitted as a PR to MRIReco.jl.