CXXFLAGS += -g -Wall -Werror --pedantic
SRC = generate_meshes.cpp
OBJ = generate_meshes

all: $(OBJ)
	@echo Compiled Generate Meshes...

$(OBJ): $(SRC) MeshGenerator.hpp
	$(CC) $(CXXFLAGS) $< -o $@
