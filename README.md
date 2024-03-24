# VLS-mapper: Long Read Mapper Based On Variable Length Seed

Release Date: March 24, 2024

Author

	~Changyong Yu (Northeastern University in CHINA)
	~Shenzhao Zheng (Northeastern University in CHINA)  
  
1.Introduction
--  
VLS mapper is a long read mapper based on variable length seeds, which can efficiently and accurately map noisy long reads generated by third-generation sequencing technology to corresponding positions on the reference genome

2.Building Notes
--
Ensure that the executable files for Minimap2, Gem-mapper, and Samtools are included in the same directory as the code.  

**Minimap2:** [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)  
**Gem-mapper:** [https://github.com/smarco/gem3-mapper](https://github.com/smarco/gem3-mapper)  
**Samtools:** [https://github.com/samtools/samtools](https://github.com/samtools/samtools)  
  

Then use g++ to compile our code:
```
$ g++ -o vlsmapper -std=c++11 main.cpp Basic.h Basic.cpp output.h output.cpp edlib/include/edlib.h edlib/src/edlib.cpp -lhts
```

3.Usage Notes
--
Use the following command to map the read dataset **read.fq** to the reference genome **reference.fq**:

```
$ ./vls-mapper -i read.fq -I reference.fq 
```

4.Parameter Settings
--

@parameter (-r,<ref_path>)
The format of a parameter of VLS-mapper in the command line is a pair of strings, here we denote the pair as (-p, [q]) or (-p,<q>). String p is the name of the parameter. String q is the value of the parameter input in the command line. [q] represents that the parameter is a optional parameter. <q> represents that the parameter is a necessary parameter.

@parameter (-i,<read_path>)

  Parameter 'i' gives the path of a text file which saves the filenames of the read dataset files. For example 'read_path'="dataset1", than dataset1 is a text file which saves the filename "read.fa".

@parameter (-I,<ref_path>)

  Parameter 'I' gives the path of a text file which saves the filenames of the reference genome files. For example 'ref_path'="dataset2", than dataset2 is a text file which saves the filename "reference.fa".

@parameter (-o,<result_path>)

  Parameter 'o' gives the path of a text file which saves the filenames of the mapping result files. For example 'result_path'="result", than result is a text file which saves the filename "map_res.sam".

@parameter (-m,[mapping_rate])

  The parameter 'm' is optional, the default value for that parameter is 0.3, which means that the probability of mapping variable length seeds to the reference genome is 30%. The input 'm' must be between 0 and 1, otherwise the default value will be used.

@parameter (-c,[chain_method])

  The parameter 'c' is optional, the default value for that parameter is 'AnchorDegree', which means that VLS-mapper uses a chain generation algorithm based on anchor degree when generating chains. The input 'c' must be 'AnchorDegree' or 'MaximalClique', otherwise the default value will be used.

5.Contacts
--

	Please e-mail your feedback at cyyneu@126.com.




