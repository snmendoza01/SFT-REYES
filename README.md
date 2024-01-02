# SFT-REYES
A repository containing the content of the REYES summer project I was involved in.


The program was a wonderful learning experience, and so I have decided to share my work through this repository. 5 different problem sets can be found in this repository (labelled by `week_i_Mendoza.ipynb`, with $i \in \{1, \ldots, 5\}$), as well as a Final Project that we made during the last week, here called `anisotropic_heisenberg_chain`. 

Each problem set was made using a template that the mentors provided that contain the quesitons, with each following cell having my answer to the respective problem. Mistakes were probably made, but I tried my best to reach a satisfactory answer each time :)


The final project contains the code I used, a report I wrote on it, and a poster which I presented to my mentors. I think the end result was really nice, but I wish I had more time to expand on the more theoretical/mathematical side because it could have yielded some pretty cool stuff. Sadly, coding/debugging took most of my time, but I still mention what I initially intended to achieve through this project. (It  also didn't help that I initially tried to work in Jax, but the randomization process and conditionals made the conversion to tracer object extremely complicated (arrays were not of fixed size), so I had to redo everything using standard NumPy techniques).

The content available is mainly:

- **Week 1**: Ising system, entropy, and partition function in a lattice. Pretty standard material, but fun to work through nonetheless.
- **Week 2**: Once again the Ising system, although going a bit more in depth. I did some hefty work to compute the 1-2 point correlation functions, but I must say it was quite rewarding. We had some symmetries/symmetry breaking, but it was very elemental
- **Week 3**: This one became more challenging. The first parts were pretty ok, with phase transitions and correlation functions being at the center of the stage. The last exercise had me do some computations to analyze the (classical) interaction between two vortices (initially not necessarily a vortex-antivortex pair). It was a little overboard, but with the help of some online resources (Tong's Lecture notes mostly) I was able to pull through and obtain what seems to be a satisfactory result. Mathematicians *beware*, I was not rigorous of course, although I did try to use some differential geometry (Stokes Theorem and considered using forms at some point) to work through the problem.
- **Week 4**: This set was my absolute favorite. I was finally able to see QM as a 1-dim QFT (and, more generally, Statistical Field Theories with their Quantum Analogs), and it was absolutely beautiful. I had a completely different perspective for the spin-$\frac{1}{2}$ particle after seeing it as the Wick rotation of an Ising chain, and solving the QHO through the path integral was absolutely magestic. 
- **Week 5**: This final week we went a bit more in depth into the numerical methods, and worked through some discretization and clustering algorithms. It was alright, although we had to replicate the results given in <a href="https://arxiv.org/pdf/hep-lat/9206004.pdf" target="_blank">this paper</a> in the end (the code I wrote was not very organized/sophisticated, but I didn't have enough time to beautify it).
