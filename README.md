# Spatially-coupled-LDPC
C simulation for spatially coupled LDPC (SC LDPC) codes on binary erasure channel (BEC) with belief propagation (BP) decoding.

The specific construction of the SC LDPC codes follows from [this paper by Kudekar et al](https://arxiv.org/abs/1001.1826). We start from the (l,r)-regular protograph, perform spatial coupling with parameter L, then randomly lift M times. The BP decoder is particularly easy to implement in the case of BEC.

Below is a plot of the EXIT curve vs the erasure rate, along with the BP and MAP threshold for the (l,r) regular ensemble. Here the SC LDPC ensemble is (l,r,L), with l=5, r=10, and L=100, 200, 500. Also, M = 100. We vary the number of BP iterations from 50 to 500.

![](/spatiallycoupledLDPC.jpg)

The SC LDPC codes move away from the BP threshold of the regular ensemble and approach closely to the MAP threshold. The MAP threshold is very close to 0.5, which is the Shannon limit on BEC for any code of rate 0.5. In the limit where M, L and (l,r) tend to infinity in that order, we can approach the Shannon capacity on BEC. The paper Kudekar et al gave a proof of this fact.
