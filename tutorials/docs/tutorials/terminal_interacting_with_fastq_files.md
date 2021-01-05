<a href="https://colab.research.google.com/github/sbooeshaghi/bytesizebio/blob/main/tutorials/staging/terminal_interacting_with_fastq_files.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Interacting with FASTQ files

FASTQ files are text files that contain biological sequences, usually nucleotides, and their associated quality scores. The quality score is a measure of how accurate the base-call is for each letter in the sequence. For more information please see the [wikipedia entry](https://en.wikipedia.org/wiki/FASTQ_format).

----------

**tl;dr** In this tutorial we introduce ways to interact with FASTQ files by the command line.

## The `curl` command

We will use the `curl` command to download FASTQ files to our machine. [`curl`](https://curl.se/) is a command line tool that allows us to download files. It takes a URL as an argument and downloads the file to your local computer in the current directory. For most cases when downloading data with `curl` we will use the following syntax:

```
$ curl -Ls -o outputfile.txt http://www.somewebsite.com/myfile.txt
```

Where the `-L` tells the command to follow "redirects" in case the file as been moved, `-s` tells the command to be silent and the `-o` tells the command the name of the file to be saved. 

## Downloading FASTQ files

FASTQ files usually come in pairs. They are coloquially called "read 1" and "read 2". This naming refers to the order in which the molecules were sequenced. We won't cover the many ways in which FASTQ files can be generated in this tutorial but for more details please check out these resources: [Illumina documentation](https://web.archive.org/web/20200414070028/https://www.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-06.pdf), [Next Generation Sequencing Core @UPenn](https://ngsc.med.upenn.edu/faqs-mini/FASTQ-Files.html).


```
!curl -Ls -o read1.fastq https://caltech.box.com/shared/static/fh81mkceb8ydwma3tlrqfgq22z4kc4nt.gz
!curl -Ls -o read2.fastq https://caltech.box.com/shared/static/ycxkluj5my7g3wiwhyq3vhv71mw5gmj5.gz
```


```
!ls -lht
```

    total 434M
    -rw-r--r-- 1 root root 247M Jan  4 23:49 read2.fastq
    -rw-r--r-- 1 root root 188M Jan  4 23:48 read1.fastq
    drwxr-xr-x 1 root root 4.0K Dec 21 17:29 sample_data


## Interacting with FASTQ files

FASTQ files can be quite large, often on the order of tens to hundreds of gigabytes, therefore it is impractical to open them up and view them in a text editor. For this reason FASTQ files are often "compressed" to save storage space. The type of compression is colloquially called "gzip" compression to refer to the command line tool that compresses files.

### Useful commands

#### The `gzip` command and the `gunzip` command

Often times FASTQ files end with `.fastq.gz`. This means that they have been compressed with the [`gzip`](https://www.gnu.org/software/gzip/) command. They can be decompressed them with the [`gunzip`](https://linux.die.net/man/1/gunzip) command. Though most of the time we don't want to decompress them since that would require a lot of space. The syntax for decompressing and compressing a file in most cases will look like the following:

```
$ gunzip file.fastq.gz 
$ gzip file.fastq
```

Our files are currently decompressed. Let's see how much space we save by compressing them.


```
!ls -lht
!gzip read1.fastq
!gzip read2.fastq
!ls -lht
```

    total 434M
    -rw-r--r-- 1 root root 247M Jan  4 23:49 read2.fastq
    -rw-r--r-- 1 root root 188M Jan  4 23:48 read1.fastq
    drwxr-xr-x 1 root root 4.0K Dec 21 17:29 sample_data
    total 74M
    -rw-r--r-- 1 root root  47M Jan  4 23:49 read2.fastq.gz
    -rw-r--r-- 1 root root  28M Jan  4 23:48 read1.fastq.gz
    drwxr-xr-x 1 root root 4.0K Dec 21 17:29 sample_data


The `read1` FASTQ file went from 188 Megabytes in size (decompressed, `read1.fastq`) to 28 Megabytes in size (compressed, `read1.fastq.gz`) and the `read2` FASTQ file went from 247 Megabytes (decompressed, `read2.fastq`) to 47 Megabytes (compressed, `read2.fastq.gz`). The `read1` file was compressed by a factor of 6.7 and the `read2` file was compressed by a factor of 5.2.


#### The `zcat` command

Since we almost always want to keep our FASTQ files compressed but we still want to be able to read them we will use a modified version of the `cat` command but for files that have been compressed with `gzip`. The [`zcat`](https://linux.die.net/man/1/zcat) command does exactly the same thing to a file as the `cat`, namely it prints the contents of the file to the screen. The difference, however is that the `zcat` command expects a compressed file. The syntax looks like this:

```
$ zcat file.fastq.gz
```

### Reading FASTQ files

If we simply use `zcat` as above it will print the whole file to the screen so instead we "pipe" the output of the zcat command to the `head` command. We don't conver pipes in this video but there are a lot of great resources covering them. Here is a [nice video](https://www.youtube.com/watch?v=bKzonnwoR2I) explaining pipes.

The `head` command just stops the process of printing after a certain number of lines. The default is 10 lines.


```
!zcat read1.fastq.gz | head
```

    @SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    ACATCNGTCGAGAACGATCGTGTCCG
    +SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    AAAAA#EEEEEEEEEEEEEEEEEEEE
    @SRR8611943.2 NS500272:478:HVL5HBGX5:1:11101:20079:1043 length=26
    TCAGCNCCAACTGCTAGTCTTTCCCT
    +SRR8611943.2 NS500272:478:HVL5HBGX5:1:11101:20079:1043 length=26
    AAAAA#EEEEEE6EEEEEEEEEEEEE
    @SRR8611943.3 NS500272:478:HVL5HBGX5:1:11101:10651:1043 length=26
    TGAAANATCACGCGGTTCATCAGTAG


We can change the number of lines printed by `head` by specifying a `-n` option:


```
!zcat read1.fastq.gz | head -n 4
```

    @SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    ACATCNGTCGAGAACGATCGTGTCCG
    +SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
    AAAAA#EEEEEEEEEEEEEEEEEEEE


### FASTQ Format

Every four lines in a FASTQ file represent one molecule that was sequenced.

```
@SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
ACATCNGTCGAGAACGATCGTGTCCG
+SRR8611943.1 NS500272:478:HVL5HBGX5:1:11101:9611:1040 length=26
AAAAA#EEEEEEEEEEEEEEEEEEEE
```

The first line is the molecule identifier, the second line is the actual sequence of the molecule, the third line starts with a `+` sign and provides optional information, and the fourth line is the quality score of the sequence. Note that the quality score must be the same length of the molecule. To learn more about quality scores please check out the [Illumina documentation](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm).

### Counting the number of molecules (reads)

A common task for bioinformaticians is to count the number of reads (also known as the number of molecules) in a FASTQ file. This means we need to count the number of lines in the FASTQ file and then divide by four, since a read is represented by four lines.

We will use the [`wc`](https://linuxize.com/post/linux-wc-command/) command to count the number of lines in a file. The `wc` command is a commandline tool to count the number of words, characters, and lines in a file. The syntax for counting lines in a file is the following:

```
$ wc -l myfile.txt
```

We will use the `wc` command to count the number of lines in the FASTQ file and then we will divide by four to get the number of reads.


```
!zcat read1.fastq.gz | wc -l
!zcat read2.fastq.gz | wc -l
```

    4000000
    4000000


Both read1 and read2 have 4,000,000 lines, this means that they have 1 million reads. For "paired-end" sequencing data where read 1 and read 2 are given, both FASTQ files will always have the same number of reads since sequneces from a single molecule are represented in both files.

### Finding an exact sequence

Another common task for bioinformaticians is to find an exact sequence in the FASTQ files. We will use the `grep` command along with the `zcat` command. The [`grep`](https://man7.org/linux/man-pages/man1/grep.1.html) command is a commandline tool that enables searching. The `grep` command can get complicated but most use cases can be covered by the following syntax:

```
$ grep "find this" myfile.txt
```

Let's find the sequence `ATTAGGAGCCG` in `read1.fastq.gz`.


```
!zcat read2.fastq.gz | grep "ATTAGGAGCCG"
```

    GAAGATGTTGTCGTGGATACTGAAATGCGTCGTCAAAAATTAGGAGCCGTTCTTTTG
    GTCAAAAATTAGGAGCCGTTCTTTTGAAGACTCTTGTTTCTCTTGGAAAGTCTCTCG
    GCGTCGTCAAAAATTAGGAGCCGTTGGTTTGAAGAATCTTTTTTCTAATGGAAAGTG
    GCGTCGTCAAAAATTAGGAGCCGTTCTTTTGAAGACTCTTGTTTCTCTTGGAAAGTC
    CAAAAATTAGGAGCCGTTCTTTTGAAGACTCTTGTTTCTCTTGGAAAGTCTCTCGGA


### Enumerating lines of a FASTQ file

In the above example we found five reads that contain our sequence of interest but we don't know what line number they came from. We can figure this out by using the `nl` command. The `nl` command is a commandline tool to add line numbers to the file. The syntax is as follows:

```
$ nl file.txt
```

Let's find out the line numbers where we found our sequence.


```
!zcat read2.fastq.gz | nl | grep "ATTAGGAGCCG"
```

    535046	GAAGATGTTGTCGTGGATACTGAAATGCGTCGTCAAAAATTAGGAGCCGTTCTTTTG
    745522	GTCAAAAATTAGGAGCCGTTCTTTTGAAGACTCTTGTTTCTCTTGGAAAGTCTCTCG
    818990	GCGTCGTCAAAAATTAGGAGCCGTTGGTTTGAAGAATCTTTTTTCTAATGGAAAGTG
    1654078	GCGTCGTCAAAAATTAGGAGCCGTTCTTTTGAAGACTCTTGTTTCTCTTGGAAAGTC
    1956810	CAAAAATTAGGAGCCGTTCTTTTGAAGACTCTTGTTTCTCTTGGAAAGTCTCTCGGA


## Summary

In this tutorial we learned how to interact with FASTQ files using built in shell commands to perform simple tasks. 


```

```
