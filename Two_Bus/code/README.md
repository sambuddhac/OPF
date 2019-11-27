An ADMM-based solver for dynamic energy management
==================================================

This is the source code for our solver. It is written in C++. There is also a 
`cvxgen` dependency which ships with this version.

This version is not meant to be public (mostly because of the `cvxgen`
dependency).

In addition, you'll need `ruby` and `protobuffers` to run the code.

To make the program, you can simply run

`make`

It will put the program in the `bin` folder.

There are also the options to compile the separate modules:

`make proto`
`make cvxgen`
`make d_opf`
`make create_network`

The first will compile the protobuffers we use in the code. The second will 
compile the `cvxgen` models. These two options will only produce `.o`
object files (meaning, their result is not very useful). The third 
will make the full binary for dynamic OPF, and the fourth makes a binary
to create sample networks.

To clean up the directory, there are two options:

`make clean`
`make fullclean`

The first just cleans up the object files used to build the binaries. The 
second also removes the generated source files from CVXGEN and Protobuffers.

TODO:
rename `opsp_solver` class to `d_opf_solver` or something related.
rename `OPF::Bus` to `DOPF::Device` or something.
