# ACTio-ReACTio
Repo for Mathematica and C++ codes used in the EFTofDE+PPF to halo model reaction map

GtoPT.nb : Notebook going from Horndeski functions to Poisson Equation modifications 

Horndeski_term_check.nb : Notebook deriving terms found in Appendix B of 1606.02520


## Docker

I have included a Dockerfile here which I've tested. Once you've installed docker [https://www.docker.com/] you can just place this file into a folder and build the images into a container

```
docker build -t mybuild . 
```

Once this builds, you can jump into the container (which has ReACT installed) using the following command

```
docker run -v /Users/bbose/Desktop/myfolder:/home -i -t mybuild
```

This will also automatically take all the local files in `/Users/bbose/Desktop/myfolder` to the `/home` folder within the container. Anything placed in the home folder from within the container will then automatically show up locally. This lets you transfer ReACT output produced within the container to the local system. 

I've tested this works on the tests and example files within the `reactions/examples` and `reactions/tests` directories. You can use the following command to run them

```
g++ -I/ReACT/reactions/include -L/ReACT/reactions/lib  spt.cpp -lgsl -lcopter
```
```
./a.out
```

To exit the container simply use the exit command 

```
exit 
```
