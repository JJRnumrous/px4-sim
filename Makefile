CC=gcc
CFLAGS=-I.
DEPS = utils.h quad_model.h sensors_model.h quad_parameters.h px4_sim_communication.h px4_quad_sim.h
LIBS=-lm

_OBJ = sensors_model.o utils.o quad_model.o px4_sim_communication.o px4_quad_sim.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

ODIR=obj

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

px4_quad_sim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o