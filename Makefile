CC    = g++
SRC1  = src/angle_calculation.cpp
SRC2  = src/distance_calculation.cpp
EXE1  = angle_calculation
EXE2  = distance_calculation

all: $(EXE1) $(EXE2)

$(EXE1): $(SRC1)
	$(CC) $(SRC1) -o $(EXE1)

$(EXE2): $(SRC2)
	$(CC) $(SRC2) -o $(EXE2)
