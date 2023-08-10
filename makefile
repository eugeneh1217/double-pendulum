CPPFLAGS=-std=c++20 -Wall -Wextra -lSDL2 -lSDL2_ttf -o dp
GPP = g++

dp: dp.cpp
	g++ $^ $(CPPFLAGS)

clean:
	rm dp
