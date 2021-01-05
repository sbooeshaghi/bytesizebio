<a href="https://colab.research.google.com/github/sbooeshaghi/bytesizebio/blob/main/tutorials/notebooks/terminal_common_fastq_tasks.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Common FASTQ tasks

FASTQ files are text files that contain biological sequences, usually nucleotides, and their associated quality scores. Often we'd like to perform various tasks on the files such as:



1.   finding the reverse complement,
2.   performing advanced sequence searching,
3.   converting FASTQ files to other formats,
4.    

----------------
**tl;dr** In this tutorial we introduce more advanced ways to manipulate and interact with FASTQ files


## Downloading FASTQ files with `curl`

We will use the `curl` command to download FASTQ files to our machine. [`curl`](https://curl.se/) is a command line tool that allows us to download files. It takes a URL as an argument and downloads the file to your local computer in the current directory.

FASTQ files usually come in pairs, known as "read 1" and "read 2". This naming convention originates from the order in which the molecules in a sequencing machine are sequenced.


```bash
%%bash
curl -Ls -o read1.fastq https://caltech.box.com/shared/static/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz
curl -Ls -o read2.fastq https://caltech.box.com/shared/static/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
ls -lht
```

    total 434M
    -rw-r--r-- 1 root root 247M Jan  5 19:34 read2.fastq
    -rw-r--r-- 1 root root 188M Jan  5 19:34 read1.fastq
    drwxr-xr-x 1 root root 4.0K Dec 21 17:29 sample_data


The `curl` command uses three options: 
1.    `-L` tells curl to follow any redirects in case the file was moved,
2.    `-s` tells curl to not print any logging information (i.e. the command is silent),
3.    and the `-o` specifies the new name of the file to be saved to our machine.

## Printing out only sequences

Often times we care only about the sequences in the FASTQ files and nothing else. We are tasked with converting the FASTQ file into just a file of sequences.

Recall that a FASTQ file is comprised of sequences of molecules and associated metadata. Each molecule is represented by four lines:

```
@SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
ACATCNGTCGAGAACGATCGTGTCCG
+SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
AAAAA#EEEEEEEEEEEEEEEEEEEE
```
The first line is the molecule identifier, the second line is the actual sequence of the molecule, the third line starts with a + sign and provides optional information, and the fourth line is the quality score of the sequence. Note that the quality score must be the same length of the molecule. 

So we need to figure out a way to print the sequence (line 2) and every 4 lines after that.

### The `awk` command

The `awk` command is a commandline tool for advanced processing of text files. It can perform many tasks and we will use it for printing out only the sequences from the FASTQ files. `awk` operates on a line-by-line basis and treats each line as a collection of "columns" separated by a "separator" such as a space, a comma, or a tab. In the simplest case we can print out every line input to awk by printing the `0th` column: 

```
awk '{print $0}' file.txt
```

Let's use it on our FASTQ files. If you recall, a fastq file that has not been compressed with `gzip` is simply a text file. Since we don't want the entire contents of the file to be printed to the screen, we can limit the number of lines with the `head` command.


```bash
%%bash
awk '{print $0}' read1.fastq | head -n 4
```

    @SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    ACATCNGTCGAGAACGATCGTGTCCG
    +SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    AAAAA#EEEEEEEEEEEEEEEEEEEE


The `awk` command also has a built-in variable called `NR`. This variable stores the current line number. We can also print it out:


```bash
%%bash
awk '{print NR, $0}' read1.fastq | head -n 12
```

    1 @SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    2 ACATCNGTCGAGAACGATCGTGTCCG
    3 +SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    4 AAAAA#EEEEEEEEEEEEEEEEEEEE
    5 @SRR8611943.2 NS500272:478:HVL5HBGX5:1:11101:20079:1043 length=26
    6 TCAGCNCCAACTGCTAGTCTTTCCCT
    7 +SRR8611943.2 NS500272:478:HVL5HBGX5:1:11101:20079:1043 length=26
    8 AAAAA#EEEEEE6EEEEEEEEEEEEE
    9 @SRR8611943.3 NS500272:478:HVL5HBGX5:1:11101:10651:1043 length=26
    10 TGAAANATCACGCGGTTCATCAGTAG
    11 +SRR8611943.3 NS500272:478:HVL5HBGX5:1:11101:10651:1043 length=26
    12 AAAAA#EEEAEEEEEEEAEEEEEEEE


### The modulo (`%`) operator 
Notice how the sequence of interest occurs on the 2nd line and then again on the 6th line and then again on the 10th line. The sequences are spaced four lines apart and start on line 2. We can be clever and print out only those lines using the modulo `%` operator. The modulo operator is the same as seen in python. Given two numbers `a` and `b`, the operation `a%b` first divides `a` by `b` and then returns the remainder. For example:


```
%%python
a = 17
b = 4
print(f"{a}/{b} = {a//b} with a remainder of {a%b}.")
print(f"In other words: {b}*{a//b} + {a%b} = {a}.")
```

    17/4 = 4 with a remainder of 1.
    In other words: 4*4 + 1 = 17.


To get every 4th line of the file we want to print out the lines when `NR%4==0`. We can declare the condition to print the line by writing the condition before the `{print}` in the following way:

```
awk '(condition){print $0}' file.txt
```

Let's try it on our file:


```bash
%%bash
awk 'NR%4==0{print NR, $0}' read1.fastq | head -n 6
```

    4 AAAAA#EEEEEEEEEEEEEEEEEEEE
    8 AAAAA#EEEEEE6EEEEEEEEEEEEE
    12 AAAAA#EEEAEEEEEEEAEEEEEEEE
    16 AAAAA#EEEEEEEEEEEEEEEEEEEE
    20 AAAAA#EEEEEEEEEEAEEEEEEEEE
    24 AAA6A#EEEEAEEEAEE6EEEEEEEE


Notice that this gives us the rows `4, 8, 12...`, i.e. the quality score. We need to shift the condition by 2 so that we get the sequence instead:



```bash
%%bash
awk 'NR%4==2{print NR, $0}' read1.fastq | head -n 6
```

    2 ACATCNGTCGAGAACGATCGTGTCCG
    6 TCAGCNCCAACTGCTAGTCTTTCCCT
    10 TGAAANATCACGCGGTTCATCAGTAG
    14 GCGCANCAGGCTACGATCGATCATGG
    18 CTGATNCTCCTTTCTCATGGAGAATG
    22 CGAGCNAAGTGAAGAGTGCGAAAGCC


### Saving the sequences to a file

Now if we want to save this as a new file we can direct the stream of text to a new file like so:


```bash
%%bash
awk 'NR%4==2{print $0}' read1.fastq > read1.txt
```

## Finding the reverse complement

DNA (Dexoy-riboNucleic Acid) is a two-stranded biopolymer where each monormer is one of the four nucleotides

1.    Adenosine (`A`)
2.    Thymine   (`T`)
3.    Guanine   (`G`)
4.    Cytosine  (`C`).

The strands have a directionality, commonly called to as the `3'` and `5'` ends which refers to the position of the carbon atoms in the pentose-sugar-ring. The `5'`-carbon usually has a phosphate group attached to it and the `3'`-carbon end usually has a hydroxyl (`OH`) group group attached to it. The nucleotides bind complementarily and run anti-parallel to each other. Here is an example:

```
5' - A A T G G A C C T A C A C T G T C A - 3'
     | | | | | | | | | | | | | | | | | |
3' - T T A C C T G G A T G T G A C A G T - 5'
```

Often bioinformaticians want to know the the reverse complement of a given sequence.

### Relevant commands

To find the reverse complement we will use two commands, `echo` and `tr`.

#### The `echo` command

The `echo` command is a commandline tool that allows us to "print" whatever argument was supplied to it. For example:


```bash
%%bash
echo "Hello world"
```

    Hello world


#### The `tr` command

The `tr` command is a commandline tool that is commonly referred to as the "translate" command. It will take in a set of letters and replace them with the new set. A common use case is to convert lowercase letters to uppercase letters. 


**Note**: `tr` takes the text that you will convert from a pipe. We don't cover pipes in this video but here is an [instructive video](https://www.youtube.com/watch?v=bKzonnwoR2I) on them.




```bash
%%bash
echo "aaaaaBBBCCCC" | tr "a" "A" 
```

    AAAAABBBCCCC


#### The `rev` command

The `rev` command is a commandline tool that reverses the input lines. For example:


```bash
%%bash
echo "aaaaaBBBCCCC" | rev
```

    CCCCBBBaaaaa


### Using `tr`, `rev`, and `echo` together

We can now `echo` the sequence we want to convert into the `tr` command to find the complementary base and then pass it into the `rev` command to reverse the sequence.


```bash
%%bash
echo "AATGGACCTACACTGTCA" | tr "ATGC" "TACG"
```

    TTACCTGGATGTGACAGT


### Processing all sequences

We are now armed with the tools to process all of the sequences in our `read1.txt` file that we generated earlier.


```bash
%%bash
cat read1.txt | tr "ATGC" "TACG" | rev | head
```

    CGGACACGATCGTTCTCGACNGATGT
    AGGGAAAGACTAGCAGTTGGNGCTGA
    CTACTGATGAACCGCGTGATNTTTCA
    CCATGATCGATCGTAGCCTGNTGCGC
    CATTCTCCATGAGAAAGGAGNATCAG
    GGCTTTCGCACTCTTCACTTNGCTCG
    AGCTCTAACTTATGCCCTGCNATTCT
    ACCGAAACGTAGTCGTAACTNGCTGA
    TTTTGGAAACTGGATGACTGNTGTCA
    GCCGATGTCCGCGCTAGACTNAAGGT


And we can save these seqences as a new file:


```bash
%%bash
cat read1.txt | tr "ATGC" "TACG" | rev > read1_revcomp.txt
```

## Summary

In summary we learned how to use the `awk` command along with the modulo (`%`) operator in order to convert a FASTQ file into a text file of just sequences. We then learned how to use the `echo`, `tr`, and `rev` command to find the reverse complement of the sequences in this file.
