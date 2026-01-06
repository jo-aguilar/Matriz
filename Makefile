CXX = gcc
LIBS = matriz.c dmpc.c teste.c
HEADER = matriz.h dmpc.h -lm 
ALVO = teste
LINKS = -lm

all: $(ALVO)

$(ALVO): $(LIBS) $(HEADER)
	$(CXX) $(LIBS) -o $(ALVO) $(LINKS)

clean: 
	rm -f $(ALVO)

run:
	clear && ./$(ALVO)

