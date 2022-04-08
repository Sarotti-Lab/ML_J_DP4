# ML_JDP4 Package

## An integrated quantum mechanics-machine learning approach for ultra-fast NMR structural elucidation

![GA](https://user-images.githubusercontent.com/101182775/162478494-d9d52f2a-65f4-4b30-97b1-8922bd877c54.png)

This repository contains all codes and data required to run ML-J-DP4 calculations. 

### Description
ML-J-DP4 is a machine learning-based Python program to compute the ML-J-DP4 probabilities using the GIAO NMR calculations conduced at the RHF/STO-3G//MMFF level as input. The script feeds on the Gaussian output files (including *.log and *.out extensions), and creates an input matrix by computing different local descriptors from the 3D geometries and NMR/NBO data. The input matrix is transformed into refined chemical shifts using a KRR trained ML. In paralell, the script calculates the 3JHH coupling constants using the HLA formalism. The chemical shifts and coupling constants are Boltzmann averaged, and correlated with the experimental data provided to obtain the ML-J-DP4 probabilities for each candidate isomer. 

### Installation Requirements

**ML_JDP4.py** needs python 3.8 or later to work. The module can be installed by console using:
`pip3 install ml-jdp4`

Usage: `ml_jdp4`

### User Guide

To run ML_JDP4 it is required that the information is located in a folder containing the following files: 

      1. All the Gaussian outputs from the NMR/NBO calculations (all conformers for all isomers). 
      
      2. An Excel file containing the experimental data and the labels of each nucleus associated with each experimental value.
      
 ### Technical requirements
 
**1) The output files:** must be named following the following convention: number_*.log or .out, where number identifies the ith isomer, ranging from 1 to N (where N is the number of candidate isomers under study). For example: 
 
       1_NewNatProd_c01.log (Conformer 1 of isomer 1 of a compound named NewNatProd)

       1_NewNatProd_c02.log (Conformer 2 of isomer 1 of a compound named NewNatProd)

       2_NewNatProd_c01.log (Conformer 1 of isomer 2 of a compound named NewNatProd)
       
       2_NewNatProd_c02.log (Conformer 2 of isomer 2 of a compound named NewNatProd)

*The script allows the use of outputs from Gaussian 03, 09 and 16.*

**2) The input Excel file:** The experimental data and the labels of the candidate structures must be provided in an Excel file which must be made according the following rules. The Excel file should be constituted by two sheets; one containing the data for the coupling constants (named ‘J’) and the other with the NMR chemical shifts (named ‘shifts’).

**“J” sheet:** the first column *“exp”* contains the experimental <sup>3</sup>J coupling constants. In case of interchangeable values, they should be arranged “upside-down” (that is, the larger value first). The second column *“exchange”* serves to indicate **0** (not interchangeable value) or **1** (experimental data interchangeable with its following value). The third and fourth columns are intended to place the labels of the coupled protons. For cases when there are isomers with different labels, there should be two columns for each isomer as indicated below.  

![image](https://user-images.githubusercontent.com/101136961/161282945-682190b8-2f04-4e53-bcbd-7e54b5dd9908.png)

**“shifts” sheet:** the first column *“nuclei”* contains the identity of the atom ‘c or C’ for <sup>13</sup>C and ‘h or H’ for hydrogen atoms. The second column *“exp_data”* contains the experimental chemical shifts. In case of interchangeable values, they should be arranged *“upside-down”* (that is, the larger value first). The third column *“exchange”* serves to indicate **0** (not interchangeable value) or **1** (experimental data interchangeable with its following value). The following columns are intended to place the labels of the nuclei associated to the corresponding chemical shift. If 2 or more values are added in that region, the isotropic shielding values will be averaged (as in the case of methyl groups or equivalent methylene groups). For cases when there are isomers with different labels, there should be three columns for each isomer as indicated below.

![image](https://user-images.githubusercontent.com/101136961/161283203-35f3f2df-e6a3-43d4-b8b4-87eb0c7bca18.png)


**3) The Output excel file:** once the ML_JDP4.py is executed, a filed named *‘Results_ML_J_DP4.xlsx’* is created in the same folder containing the Gaussian output files and the Excel input file. The Excel output file contains five sheets: 


**DP4 sheet:**  the DP4 probabilities are shown for each isomer considering the information of H, C and J individually, and altogether. Although the high accuracy in the ML predictions, it must be emphasized that some environments might not be correctly reproduced leading to large unscaled errors that would affect the scaling procedure and the concomitant J-DP4 values. Hence, to avoid potential misassignments, the following sheets contain information regrding the scaled and unscaled chemical shifts, and the corresponding errors (differences with the experimental values). In case all isomers display alarmingly high errors for a given nucleous, it would be advisable to re-compute J-DP4 after removing or revising the conflicting signal. 

![dp4](https://user-images.githubusercontent.com/101182775/162080130-6fe10393-efaf-4871-b2df-2b1516adc063.png)


**Unscaled Chemical Shifts sheet:** this sheet is labeled as *“Shifts_Unsc”*, and displays the experimental chemical shifts (<sup>13</sup>C and <sup>1</sup>H) and the coupling constants (<sup>3</sup>J), in that order, along with the predicted Boltzmann-averaged unscaled values computed for each isomer.

![unsc](https://user-images.githubusercontent.com/101182775/162080256-50f92f90-7d89-4ecf-8c9d-2bd405b456d5.png)


**Scaled Chemical Shifts sheet:** this sheet is labeled as *“Shifts”*, and displays the experimental chemical shifts (<sup>13</sup>C and <sup>1</sup>H) and the coupling constants (<sup>3</sup>J), in that order, along with the corresponding scaled values after for each isomer. The scaling procedure is done according to δs = (δu − b)/m, where m and b are the slope and intercept, respectively, resulting from a linear regression calculation on a plot of δu against δexp.

![sc](https://user-images.githubusercontent.com/101182775/162080321-678e630b-ed2d-44be-b0a7-d9313a3f296e.png)


**Unscaled Errors sheet:** this sheet is labeled as *“Unsc_Errors”*, and displays the differences (absolute value) between unscaled and experimental chemical shifts and coupling constants. It is necessary to check this sheet to verify the absence of large outliers that could affect the assignment. In case all isomers display alarmingly high errors for a given nucleous, it would be advisable to re-compute J-DP4 after removing or revising the conflicting signal. 

![unErr](https://user-images.githubusercontent.com/101182775/162080373-90b639d9-c309-4e47-966c-cdda836836e8.png)

**Errors sheet:** it displays the differences (absolute value) between scaled and experimental chemical shifts and coupling constants. 

![ScErr](https://user-images.githubusercontent.com/101182775/162080396-15ba6550-044b-4666-a0dd-1d7891fcbb59.png)

# Case study: (-)-Menthol

In order to illustrate the ML-dJ-DP4 and the ML-iJ/dJ-DP4 workflows, we present the analysis of (-)-menthol following these two approaches. As indicated in the Figure 1, there are four possible diastereoisomers. 
 
![menthol](https://user-images.githubusercontent.com/101182775/162082954-d6a71d19-8e3d-4520-91d2-edc0104c8729.jpg)
**Figure 1**

Following the recommended computational procedure, a total number of 123 conformers were found at the MMFF force field. Each structure was submitted to NMR calcualtions at the RHF/STO-3G level (with the pop=nbo option). The corresponding output files are provided in the Folder **“menthol_ML_dJ-DP4”**. According to Gaussian numbering scheme, the labels corresponding to each nuclei are given in Figure 2. 

![menthol-numeracion](https://user-images.githubusercontent.com/101182775/162085447-f6a68302-acd4-4545-9b9e-e1e370971017.jpg)
**Figure 2.** Carbon labels, and proton lavels (in parenthesis). 

## ML-dJ-DP4

The input excel file required for running the ML_JDP4.py script must be filled as follows. The Excel file is also provided in the folder *“menthol_ML_dJ-DP4”*.

**J sheet**

![Jdir (2)](https://user-images.githubusercontent.com/101182775/162087115-c72f178b-2cfe-4b12-8694-1e47e4e49ca8.png)

**Shifts sheet**

![shiftdir (2)](https://user-images.githubusercontent.com/101182775/162087211-25112778-554f-42fd-b56a-0f872a915f68.png)

Once installed, the script is run through the console as `ml_jdp4`. A pop-up window is opened, requesting to select the directory in which the Gaussian output files and the Excel input file are located. Once the folder is selected, a new pop-up window is opened requesting to select the Excel input file with the experimental data and labels. Once the script is run, the resulting excel report file **“Results_ML_J_DP4”** will be generated. The output Excel file is provided in the folder **“menthol_ML_dJ-DP4”**.

**DP4 sheet:**

![dp4](https://user-images.githubusercontent.com/101182775/162087380-d2445cec-3e2c-4ec6-a206-e402060628fe.png)

## ML-iJ/dJ-DP4

This approach involves removing unsuitable conformations that are incompatible with selected experimental <sup>3</sup>J<sub>HH</sub> values. In this case, we selected the J<sub>16,19</sub>= 10 Hz to constrain the sampling, leading to 34 conformations.

![mentholRemove](https://user-images.githubusercontent.com/101182775/162088111-620e4fd3-ea72-4adf-9860-0ba770069aae.jpg)
 
**Figure 3.** J<sub>16,19</sub> (10 Hz) was chosen to constrain the conformations

Once removed the unsuitable conformations, the remaining conformations are submitted to GIAO NMR calculations at the RHF/STO-3G level with the pop=nbo option. The resulting output files are given in the folder **“menthol_ML_iJ-dJ-DP4”**. It is important to note that in this method, the selection of the conformations should be done before running the `ML_JDP4.py script`. The Excel file with the experimental data is the same as in the dJ-DP4 analysis, but leaving out the information of the <sup>3</sup>J<sub>HH</sub> used to guide the conformational sampling. 

**J sheet:**

![iJJota (2)](https://user-images.githubusercontent.com/101182775/162088779-68c544f3-8361-4ce7-819f-93c7183321f2.png)

Once the script is run (in the same way discussed above), the resulting excel report file **“Results_ML_J_DP4” will be generated**. The output Excel file is provided in the folder *“menthol_ML_iJ-dJ-DP4”*.

**DP4 sheet:**

![iJdp4 (2)](https://user-images.githubusercontent.com/101182775/162088878-4d97c35c-305b-4d2d-b14c-c5a6fc9d430c.png)





