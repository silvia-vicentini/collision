# collision
To execute analysis.cpp:
Compile from ROOT whith the following commands:
```bash
$ .L particletype.cpp
$ .L resonancetype.cpp
$ .L particle.cpp
$ .L simulatin.cpp
$ simulation()
$ .L analysis.cpp
$ analysis()
```

To execute test.cpp:
Compile from SHELL whith the following commands:
```bash
$ g++  -Wall -Wextra particletype.hpp particletype.cpp resonancetype.hpp resonancetype.cpp particle.hpp particle.cpp test.cpp -o name_of_progect
$ ./name_of_progect
```
