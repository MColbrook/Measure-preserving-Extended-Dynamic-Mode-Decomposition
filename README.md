# Measure-preserving-Extended-Dynamic-Mode-Decomposition
General-purpose code and examples for the mpEDMD algorithm. Computation of spectral properties of Koopman operators associated with arbitrary measure-preserving dynamical systems. Highlights include: convergence to spectral properties of Koopman operator, measure-preserving discretization, enhanced robustness to noise, and improved prediction and qualitative behavior of the Koopman mode decomposition. The algorithm is versatile and can work with any choice of dictionary or any shape/format of the dataset.

The reference paper can be found here:<br>
[The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving Dynamical Systems](https://epubs.siam.org/doi/abs/10.1137/22M1521407?journalCode=sjnaam)<br>
Please send comments, bugs and suggestions to: m.colbrook@damtp.cam.ac.uk

The code includes a main routine **mpEDMD** as well as the three examples from the paper. Other examples will be added in the future. If you have an example and want to add it, please get in touch! 

**Datasets** (needed for turbulent flow example) can be found here: https://www.dropbox.com/sh/xj59e5in7dfsobi/AAAfkxqa1x9WFSTgrvqoqqRqa?dl=0

Some of the code for setting up the examples makes use of Chebfun, which can be found at https://www.chebfun.org/.
