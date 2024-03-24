# VLS-mapper: Long Read Mapper Based On Variable Length Seed


## Usage
Ensure that the executable files for Minimap2, Gem-mapper, and Samtools are included in the same directory as the code.  

**Minimap2:** [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)  
**Gem-mapper:** [https://github.com/smarco/gem3-mapper](https://github.com/smarco/gem3-mapper)  
**Samtools:** [https://github.com/samtools/samtools](https://github.com/samtools/samtools)  
  

Then use g++ to compile our code:
```
$ g++ -o vlsmapper -std=c++11 main.cpp Basic.h Basic.cpp output.h output.cpp edlib/include/edlib.h edlib/src/edlib.cpp -lhts
```

Use the following command to map the read dataset **read.fq** to the reference genome **reference.fq**:

```
$ ./vls-mapper -i read.fq -I reference.fq 
```




