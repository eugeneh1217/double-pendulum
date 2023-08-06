CPPFLAGS=-w -lSDL2 -o dp
GPP = g++

dp: dp.cpp
	g++ $^ $(CPPFLAGS)

clean:
	rm dp
