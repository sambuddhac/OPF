# CVXGEN devices
DEVICES		= generator battery deferrable_load
# C++ devices
DEVICESRC	= opsp_lossy_dc_line.cpp opsp_battery.cpp opsp_curtailable_load.cpp opsp_deferrable_load.cpp opsp_fixed_load.cpp opsp_gen_load.cpp opsp_generator.cpp opsp_lossless_line.cpp opsp_lossless_quadratic_line.cpp opsp_lossy_line.cpp

# directories
BUILDDIR	= build
BINDIR		= bin
SRCDIR		= src
BUILDCVXDIRS	= $(addprefix $(BUILDDIR)/,$(DEVICES))
PROTODIR	= proto
CVXGENDIR	= rb_cvxgen

# compiler
CC			= g++

INCLUDES	= -I $(SRCDIR)/$(PROTODIR) $(addprefix -I $(SRCDIR)/,$(DEVICES))
LDFLAGS		= -lgomp -lprotobuf
CFLAGS  	= -g -p -Wall -mavx -march=corei7-avx -ffast-math -O3 -fopenmp

# CVXGEN source
CVXGENSRC	= solver.c ldl.c matrix_support.c
CVXGENOBJ	= $(foreach DIR, $(BUILDCVXDIRS), $(addprefix $(DIR)/, $(CVXGENSRC:.c=.o)))

# PROTOBUF source
PROTOSRC	= power_network.pb.cc
PROTOOBJ	= $(BUILDDIR)/$(PROTODIR)/$(PROTOSRC:.cc=.o)

# D_OPF source
OPFSRC		= device.cpp terminal.cpp net.cpp opsp_solver.cpp $(addprefix devices/, $(DEVICESRC))
OPFOBJ		= $(addprefix $(BUILDDIR)/, $(OPFSRC:.cpp=.o))

# device headers (we omit dependencies on generated header files)
HEADERS = $(wildcard $(SRCDIR)/*.hpp) $(wildcard $(SRCDIR)/*.h) $(wildcard $(SRCDIR)/devices/*.hpp)

.PHONY: all
all: proto cvxgen d_opf create_network perturb_network

# ADMM based optimal power flow solver
.PHONY: cvxgen
cvxgen: $(BUILDDIR) $(BUILDCVXDIRS) $(addprefix $(SRCDIR)/,$(DEVICES))  $(CVXGENOBJ)

.PHONY: proto
proto: $(BUILDDIR) $(BUILDDIR)/$(PROTODIR) $(SRCDIR)/$(PROTODIR)/$(PROTOSRC) $(PROTOOBJ)

.PHONY: d_opf
d_opf: $(BUILDDIR) $(BUILDDIR)/devices $(BINDIR) cvxgen proto $(OPFOBJ) $(BINDIR)/d_opf

.PHONY: create_network
create_network: $(BUILDDIR) $(BUILDDIR)/devices $(BINDIR) proto $(BINDIR)/create_network

.PHONY: perturb_network
perturb_network: $(BUILDDIR) $(BUILDDIR)/devices $(BINDIR) proto $(BINDIR)/perturb_network

$(BINDIR)/%: $(BUILDDIR)/%.o $(CVXGENOBJ) $(OPFOBJ) $(PROTOOBJ)
	$(CC) $< $(CVXGENOBJ) $(OPFOBJ) $(PROTOOBJ) -o $@ $(LDFLAGS) -pg

%: $(CVXGENDIR)/%.cvx
	# rm -rf $@
	# ruby $(SRCDIR)/$(CVXGENDIR)/cvxgen_api.rb $<
	# ruby $(SRCDIR)/$(CVXGENDIR)/cvxgen_convert_to_class.rb $(@:$(SRCDIR)/%=%)
	# mv cvxgen $@

%.pb.cc: %.proto
	protoc -I=$(dir $@) --cpp_out=$(dir $@) $<

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDES) $< -pg

$(BUILDDIR)/%.pb.o: $(SRCDIR)/%.pb.cc
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDES) $< -pg

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDES) $< -pg

$(BUILDDIR)/%:
	mkdir $@

$(BUILDDIR):
	mkdir $(BUILDDIR)

$(BINDIR):
	mkdir $(BINDIR)

# Cleanup tasks
.PHONY: clean
clean:
	rm -rf $(BINDIR) $(BUILDDIR)

.PHONY: fullclean
fullclean:
	rm -rf $(BINDIR) $(BUILDDIR) $(addprefix $(SRCDIR)/,$(DEVICES)) $(SRCDIR)/$(PROTODIR)/*.cc $(SRCDIR)/$(PROTODIR)/*.h

