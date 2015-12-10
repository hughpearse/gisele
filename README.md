# gisele

## Description
Automated Data Model and State Machine Inference of Previously Unobserved Network Protocols.


## Introduction
The aim of this project is to create a tool that can analyse a packet capture file which contains plaintext data in the OSI layer 7 of each packet, and then output an xml document describing the structure of the protocol. This could be used to analyse previously unobserved protocols, which would be a common occurence for a malware analyst.

This tool uses machine learning techniques to produce the xml document. The first step is to apply a multiple sequence alignment algorithm commonly used for aligning multiple DNA sequences. This ensures the data is optimally aligned for analysis. Then from the sequence alignment algorithm the distances (not the similarities) between each packet are used to populate a matrix where each row/column is a packet and each cell contains the distance value. Next a multidimensional scaling algorithm uses the distance matrix to plot each packet on an XY plane. The points can then be clustered using traditional clustering algorithms. Each cluster is then re-aligned using the multiple sequence alignment algorithm and packets are recursively merged to form the protocol structure. The final result is an xml document that follows the Peach pit xml file structure. This means the process of describing and fuzzing a protocol for vulnerabilities can be fully automated from start to finish.

## Instructions
The program is called graph.py

The third party libraries can be found in the "lib" folder.

The environment setup for the eclipse project can be found in the 2 screenshots "program-arguments.png" and "third-party-libs.png" in the "INSTRUCTIONS" folder.


