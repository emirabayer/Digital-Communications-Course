# Digital Communications Course (EEE431)

This repository contains the computational assignments and reports for the Digital Communications course at Bilkent University. The work focuses on the mathematical modeling and simulation of digital communication systems, from source coding to modulation and error analysis in noise.

***course link:*** https://catalog.bilkent.edu.tr/course/c14431.html

<br>

## Computational Assignment 1 - Source Coding, Quantization & Image Processing
This project explores the conversion of continuous sources into digital formats through various quantization techniques. 
- **Scalar Quantization:** Generated a triangularly distributed source using **Inverse Transform Sampling**. Implemented both **Uniform** and **Lloyd-Max** (iterative centroid) quantizers, demonstrating up to a 1 dB improvement in SQNR with the optimal Lloyd-Max algorithm[cite: 2200, 2295, 2300].
- **Image Quantization:** Analyzed image statistics by fitting 4th-order polynomial PDFs to histograms[cite: 2324]. 
- **Frequency Domain Analysis:** Utilized 2D-FFT to visualize high-frequency quantization noise and implemented rectangular **Low-Pass Filters (LPF)** to effectively suppress noise in quantized images[cite: 2409, 2412].

<br>

## Computational Assignment 2 - Modulation, Optimal Receivers & Error Performance
This assignment focuses on the design and performance evaluation of receivers for AWGN channels using MATLAB Monte Carlo simulations ($10^6$ realizations)[cite: 2480].
- **Binary Antipodal Modulation:** Designed an optimal correlation receiver for a parabolic pulse shape. Derived and verified the **Theoretical Symbol Error Probability** against simulated performance for both equal and unequal priors (MAP vs. ML rules)[cite: 2477, 2542, 2586].
- **M-ary Constellations:** Analyzed two distinct 4-symbol constellations in 2D signal space. Derived **Union Bounds** for SEP and demonstrated how constellation geometry affects performance[cite: 2621, 2676].
- **Coding Schemes:** Evaluated Bit Error Rate (BER) performance comparing **Natural Coding** versus **Gray Coding**, verifying that Gray coding significantly minimizes bit errors by mapping adjacent symbols to a Hamming distance of 1[cite: 2729, 2757].
