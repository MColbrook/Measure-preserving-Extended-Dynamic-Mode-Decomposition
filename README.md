# Measure-preserving-Extended-Dynamic-Mode-Decomposition
General purpose code and examples for the mpEDMD algorithm. Computation of spectral properties of Koopman operators associated with measure-preserving dynamical systems. Highlights include: convergence to spectral properties of Koopman operator, measure-preserving discretization, and increased robustness to noise. The algorithm can work with any choice of dictionary or shape/format of data set.

The reference paper can be found here:<br>
[https://arxiv.org/pdf/2209.02244.pdf<br>](https://epubs.siam.org/doi/abs/10.1137/22M1521407?journalCode=sjnaam)
Please send comments, bugs and suggestions to: m.colbrook@damtp.cam.ac.uk

The code includes a main routine **mpEDMD** as well as the three examples from the paper. Other examples will be added in the future. If you have an example and want to add it, please get in touch! 

**Datasets** (needed for turbulent flow example) can be found here: https://www.dropbox.com/sh/xj59e5in7dfsobi/AAAfkxqa1x9WFSTgrvqoqqRqa?dl=0

Some of the code for setting up the examples makes use of Chebfun, which can be found at https://www.chebfun.org/.
